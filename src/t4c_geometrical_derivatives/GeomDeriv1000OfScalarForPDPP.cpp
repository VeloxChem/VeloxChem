#include "GeomDeriv1000OfScalarForPDPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_pdpp_0(CSimdArray<double>& buffer_1000_pdpp,
                     const CSimdArray<double>& buffer_sdpp,
                     const CSimdArray<double>& buffer_ddpp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_pdpp.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_pdpp

    auto g_x_0_0_0_x_xx_x_x = buffer_1000_pdpp[0];

    auto g_x_0_0_0_x_xx_x_y = buffer_1000_pdpp[1];

    auto g_x_0_0_0_x_xx_x_z = buffer_1000_pdpp[2];

    auto g_x_0_0_0_x_xx_y_x = buffer_1000_pdpp[3];

    auto g_x_0_0_0_x_xx_y_y = buffer_1000_pdpp[4];

    auto g_x_0_0_0_x_xx_y_z = buffer_1000_pdpp[5];

    auto g_x_0_0_0_x_xx_z_x = buffer_1000_pdpp[6];

    auto g_x_0_0_0_x_xx_z_y = buffer_1000_pdpp[7];

    auto g_x_0_0_0_x_xx_z_z = buffer_1000_pdpp[8];

    auto g_x_0_0_0_x_xy_x_x = buffer_1000_pdpp[9];

    auto g_x_0_0_0_x_xy_x_y = buffer_1000_pdpp[10];

    auto g_x_0_0_0_x_xy_x_z = buffer_1000_pdpp[11];

    auto g_x_0_0_0_x_xy_y_x = buffer_1000_pdpp[12];

    auto g_x_0_0_0_x_xy_y_y = buffer_1000_pdpp[13];

    auto g_x_0_0_0_x_xy_y_z = buffer_1000_pdpp[14];

    auto g_x_0_0_0_x_xy_z_x = buffer_1000_pdpp[15];

    auto g_x_0_0_0_x_xy_z_y = buffer_1000_pdpp[16];

    auto g_x_0_0_0_x_xy_z_z = buffer_1000_pdpp[17];

    auto g_x_0_0_0_x_xz_x_x = buffer_1000_pdpp[18];

    auto g_x_0_0_0_x_xz_x_y = buffer_1000_pdpp[19];

    auto g_x_0_0_0_x_xz_x_z = buffer_1000_pdpp[20];

    auto g_x_0_0_0_x_xz_y_x = buffer_1000_pdpp[21];

    auto g_x_0_0_0_x_xz_y_y = buffer_1000_pdpp[22];

    auto g_x_0_0_0_x_xz_y_z = buffer_1000_pdpp[23];

    auto g_x_0_0_0_x_xz_z_x = buffer_1000_pdpp[24];

    auto g_x_0_0_0_x_xz_z_y = buffer_1000_pdpp[25];

    auto g_x_0_0_0_x_xz_z_z = buffer_1000_pdpp[26];

    auto g_x_0_0_0_x_yy_x_x = buffer_1000_pdpp[27];

    auto g_x_0_0_0_x_yy_x_y = buffer_1000_pdpp[28];

    auto g_x_0_0_0_x_yy_x_z = buffer_1000_pdpp[29];

    auto g_x_0_0_0_x_yy_y_x = buffer_1000_pdpp[30];

    auto g_x_0_0_0_x_yy_y_y = buffer_1000_pdpp[31];

    auto g_x_0_0_0_x_yy_y_z = buffer_1000_pdpp[32];

    auto g_x_0_0_0_x_yy_z_x = buffer_1000_pdpp[33];

    auto g_x_0_0_0_x_yy_z_y = buffer_1000_pdpp[34];

    auto g_x_0_0_0_x_yy_z_z = buffer_1000_pdpp[35];

    auto g_x_0_0_0_x_yz_x_x = buffer_1000_pdpp[36];

    auto g_x_0_0_0_x_yz_x_y = buffer_1000_pdpp[37];

    auto g_x_0_0_0_x_yz_x_z = buffer_1000_pdpp[38];

    auto g_x_0_0_0_x_yz_y_x = buffer_1000_pdpp[39];

    auto g_x_0_0_0_x_yz_y_y = buffer_1000_pdpp[40];

    auto g_x_0_0_0_x_yz_y_z = buffer_1000_pdpp[41];

    auto g_x_0_0_0_x_yz_z_x = buffer_1000_pdpp[42];

    auto g_x_0_0_0_x_yz_z_y = buffer_1000_pdpp[43];

    auto g_x_0_0_0_x_yz_z_z = buffer_1000_pdpp[44];

    auto g_x_0_0_0_x_zz_x_x = buffer_1000_pdpp[45];

    auto g_x_0_0_0_x_zz_x_y = buffer_1000_pdpp[46];

    auto g_x_0_0_0_x_zz_x_z = buffer_1000_pdpp[47];

    auto g_x_0_0_0_x_zz_y_x = buffer_1000_pdpp[48];

    auto g_x_0_0_0_x_zz_y_y = buffer_1000_pdpp[49];

    auto g_x_0_0_0_x_zz_y_z = buffer_1000_pdpp[50];

    auto g_x_0_0_0_x_zz_z_x = buffer_1000_pdpp[51];

    auto g_x_0_0_0_x_zz_z_y = buffer_1000_pdpp[52];

    auto g_x_0_0_0_x_zz_z_z = buffer_1000_pdpp[53];

    auto g_x_0_0_0_y_xx_x_x = buffer_1000_pdpp[54];

    auto g_x_0_0_0_y_xx_x_y = buffer_1000_pdpp[55];

    auto g_x_0_0_0_y_xx_x_z = buffer_1000_pdpp[56];

    auto g_x_0_0_0_y_xx_y_x = buffer_1000_pdpp[57];

    auto g_x_0_0_0_y_xx_y_y = buffer_1000_pdpp[58];

    auto g_x_0_0_0_y_xx_y_z = buffer_1000_pdpp[59];

    auto g_x_0_0_0_y_xx_z_x = buffer_1000_pdpp[60];

    auto g_x_0_0_0_y_xx_z_y = buffer_1000_pdpp[61];

    auto g_x_0_0_0_y_xx_z_z = buffer_1000_pdpp[62];

    auto g_x_0_0_0_y_xy_x_x = buffer_1000_pdpp[63];

    auto g_x_0_0_0_y_xy_x_y = buffer_1000_pdpp[64];

    auto g_x_0_0_0_y_xy_x_z = buffer_1000_pdpp[65];

    auto g_x_0_0_0_y_xy_y_x = buffer_1000_pdpp[66];

    auto g_x_0_0_0_y_xy_y_y = buffer_1000_pdpp[67];

    auto g_x_0_0_0_y_xy_y_z = buffer_1000_pdpp[68];

    auto g_x_0_0_0_y_xy_z_x = buffer_1000_pdpp[69];

    auto g_x_0_0_0_y_xy_z_y = buffer_1000_pdpp[70];

    auto g_x_0_0_0_y_xy_z_z = buffer_1000_pdpp[71];

    auto g_x_0_0_0_y_xz_x_x = buffer_1000_pdpp[72];

    auto g_x_0_0_0_y_xz_x_y = buffer_1000_pdpp[73];

    auto g_x_0_0_0_y_xz_x_z = buffer_1000_pdpp[74];

    auto g_x_0_0_0_y_xz_y_x = buffer_1000_pdpp[75];

    auto g_x_0_0_0_y_xz_y_y = buffer_1000_pdpp[76];

    auto g_x_0_0_0_y_xz_y_z = buffer_1000_pdpp[77];

    auto g_x_0_0_0_y_xz_z_x = buffer_1000_pdpp[78];

    auto g_x_0_0_0_y_xz_z_y = buffer_1000_pdpp[79];

    auto g_x_0_0_0_y_xz_z_z = buffer_1000_pdpp[80];

    auto g_x_0_0_0_y_yy_x_x = buffer_1000_pdpp[81];

    auto g_x_0_0_0_y_yy_x_y = buffer_1000_pdpp[82];

    auto g_x_0_0_0_y_yy_x_z = buffer_1000_pdpp[83];

    auto g_x_0_0_0_y_yy_y_x = buffer_1000_pdpp[84];

    auto g_x_0_0_0_y_yy_y_y = buffer_1000_pdpp[85];

    auto g_x_0_0_0_y_yy_y_z = buffer_1000_pdpp[86];

    auto g_x_0_0_0_y_yy_z_x = buffer_1000_pdpp[87];

    auto g_x_0_0_0_y_yy_z_y = buffer_1000_pdpp[88];

    auto g_x_0_0_0_y_yy_z_z = buffer_1000_pdpp[89];

    auto g_x_0_0_0_y_yz_x_x = buffer_1000_pdpp[90];

    auto g_x_0_0_0_y_yz_x_y = buffer_1000_pdpp[91];

    auto g_x_0_0_0_y_yz_x_z = buffer_1000_pdpp[92];

    auto g_x_0_0_0_y_yz_y_x = buffer_1000_pdpp[93];

    auto g_x_0_0_0_y_yz_y_y = buffer_1000_pdpp[94];

    auto g_x_0_0_0_y_yz_y_z = buffer_1000_pdpp[95];

    auto g_x_0_0_0_y_yz_z_x = buffer_1000_pdpp[96];

    auto g_x_0_0_0_y_yz_z_y = buffer_1000_pdpp[97];

    auto g_x_0_0_0_y_yz_z_z = buffer_1000_pdpp[98];

    auto g_x_0_0_0_y_zz_x_x = buffer_1000_pdpp[99];

    auto g_x_0_0_0_y_zz_x_y = buffer_1000_pdpp[100];

    auto g_x_0_0_0_y_zz_x_z = buffer_1000_pdpp[101];

    auto g_x_0_0_0_y_zz_y_x = buffer_1000_pdpp[102];

    auto g_x_0_0_0_y_zz_y_y = buffer_1000_pdpp[103];

    auto g_x_0_0_0_y_zz_y_z = buffer_1000_pdpp[104];

    auto g_x_0_0_0_y_zz_z_x = buffer_1000_pdpp[105];

    auto g_x_0_0_0_y_zz_z_y = buffer_1000_pdpp[106];

    auto g_x_0_0_0_y_zz_z_z = buffer_1000_pdpp[107];

    auto g_x_0_0_0_z_xx_x_x = buffer_1000_pdpp[108];

    auto g_x_0_0_0_z_xx_x_y = buffer_1000_pdpp[109];

    auto g_x_0_0_0_z_xx_x_z = buffer_1000_pdpp[110];

    auto g_x_0_0_0_z_xx_y_x = buffer_1000_pdpp[111];

    auto g_x_0_0_0_z_xx_y_y = buffer_1000_pdpp[112];

    auto g_x_0_0_0_z_xx_y_z = buffer_1000_pdpp[113];

    auto g_x_0_0_0_z_xx_z_x = buffer_1000_pdpp[114];

    auto g_x_0_0_0_z_xx_z_y = buffer_1000_pdpp[115];

    auto g_x_0_0_0_z_xx_z_z = buffer_1000_pdpp[116];

    auto g_x_0_0_0_z_xy_x_x = buffer_1000_pdpp[117];

    auto g_x_0_0_0_z_xy_x_y = buffer_1000_pdpp[118];

    auto g_x_0_0_0_z_xy_x_z = buffer_1000_pdpp[119];

    auto g_x_0_0_0_z_xy_y_x = buffer_1000_pdpp[120];

    auto g_x_0_0_0_z_xy_y_y = buffer_1000_pdpp[121];

    auto g_x_0_0_0_z_xy_y_z = buffer_1000_pdpp[122];

    auto g_x_0_0_0_z_xy_z_x = buffer_1000_pdpp[123];

    auto g_x_0_0_0_z_xy_z_y = buffer_1000_pdpp[124];

    auto g_x_0_0_0_z_xy_z_z = buffer_1000_pdpp[125];

    auto g_x_0_0_0_z_xz_x_x = buffer_1000_pdpp[126];

    auto g_x_0_0_0_z_xz_x_y = buffer_1000_pdpp[127];

    auto g_x_0_0_0_z_xz_x_z = buffer_1000_pdpp[128];

    auto g_x_0_0_0_z_xz_y_x = buffer_1000_pdpp[129];

    auto g_x_0_0_0_z_xz_y_y = buffer_1000_pdpp[130];

    auto g_x_0_0_0_z_xz_y_z = buffer_1000_pdpp[131];

    auto g_x_0_0_0_z_xz_z_x = buffer_1000_pdpp[132];

    auto g_x_0_0_0_z_xz_z_y = buffer_1000_pdpp[133];

    auto g_x_0_0_0_z_xz_z_z = buffer_1000_pdpp[134];

    auto g_x_0_0_0_z_yy_x_x = buffer_1000_pdpp[135];

    auto g_x_0_0_0_z_yy_x_y = buffer_1000_pdpp[136];

    auto g_x_0_0_0_z_yy_x_z = buffer_1000_pdpp[137];

    auto g_x_0_0_0_z_yy_y_x = buffer_1000_pdpp[138];

    auto g_x_0_0_0_z_yy_y_y = buffer_1000_pdpp[139];

    auto g_x_0_0_0_z_yy_y_z = buffer_1000_pdpp[140];

    auto g_x_0_0_0_z_yy_z_x = buffer_1000_pdpp[141];

    auto g_x_0_0_0_z_yy_z_y = buffer_1000_pdpp[142];

    auto g_x_0_0_0_z_yy_z_z = buffer_1000_pdpp[143];

    auto g_x_0_0_0_z_yz_x_x = buffer_1000_pdpp[144];

    auto g_x_0_0_0_z_yz_x_y = buffer_1000_pdpp[145];

    auto g_x_0_0_0_z_yz_x_z = buffer_1000_pdpp[146];

    auto g_x_0_0_0_z_yz_y_x = buffer_1000_pdpp[147];

    auto g_x_0_0_0_z_yz_y_y = buffer_1000_pdpp[148];

    auto g_x_0_0_0_z_yz_y_z = buffer_1000_pdpp[149];

    auto g_x_0_0_0_z_yz_z_x = buffer_1000_pdpp[150];

    auto g_x_0_0_0_z_yz_z_y = buffer_1000_pdpp[151];

    auto g_x_0_0_0_z_yz_z_z = buffer_1000_pdpp[152];

    auto g_x_0_0_0_z_zz_x_x = buffer_1000_pdpp[153];

    auto g_x_0_0_0_z_zz_x_y = buffer_1000_pdpp[154];

    auto g_x_0_0_0_z_zz_x_z = buffer_1000_pdpp[155];

    auto g_x_0_0_0_z_zz_y_x = buffer_1000_pdpp[156];

    auto g_x_0_0_0_z_zz_y_y = buffer_1000_pdpp[157];

    auto g_x_0_0_0_z_zz_y_z = buffer_1000_pdpp[158];

    auto g_x_0_0_0_z_zz_z_x = buffer_1000_pdpp[159];

    auto g_x_0_0_0_z_zz_z_y = buffer_1000_pdpp[160];

    auto g_x_0_0_0_z_zz_z_z = buffer_1000_pdpp[161];

    auto g_y_0_0_0_x_xx_x_x = buffer_1000_pdpp[162];

    auto g_y_0_0_0_x_xx_x_y = buffer_1000_pdpp[163];

    auto g_y_0_0_0_x_xx_x_z = buffer_1000_pdpp[164];

    auto g_y_0_0_0_x_xx_y_x = buffer_1000_pdpp[165];

    auto g_y_0_0_0_x_xx_y_y = buffer_1000_pdpp[166];

    auto g_y_0_0_0_x_xx_y_z = buffer_1000_pdpp[167];

    auto g_y_0_0_0_x_xx_z_x = buffer_1000_pdpp[168];

    auto g_y_0_0_0_x_xx_z_y = buffer_1000_pdpp[169];

    auto g_y_0_0_0_x_xx_z_z = buffer_1000_pdpp[170];

    auto g_y_0_0_0_x_xy_x_x = buffer_1000_pdpp[171];

    auto g_y_0_0_0_x_xy_x_y = buffer_1000_pdpp[172];

    auto g_y_0_0_0_x_xy_x_z = buffer_1000_pdpp[173];

    auto g_y_0_0_0_x_xy_y_x = buffer_1000_pdpp[174];

    auto g_y_0_0_0_x_xy_y_y = buffer_1000_pdpp[175];

    auto g_y_0_0_0_x_xy_y_z = buffer_1000_pdpp[176];

    auto g_y_0_0_0_x_xy_z_x = buffer_1000_pdpp[177];

    auto g_y_0_0_0_x_xy_z_y = buffer_1000_pdpp[178];

    auto g_y_0_0_0_x_xy_z_z = buffer_1000_pdpp[179];

    auto g_y_0_0_0_x_xz_x_x = buffer_1000_pdpp[180];

    auto g_y_0_0_0_x_xz_x_y = buffer_1000_pdpp[181];

    auto g_y_0_0_0_x_xz_x_z = buffer_1000_pdpp[182];

    auto g_y_0_0_0_x_xz_y_x = buffer_1000_pdpp[183];

    auto g_y_0_0_0_x_xz_y_y = buffer_1000_pdpp[184];

    auto g_y_0_0_0_x_xz_y_z = buffer_1000_pdpp[185];

    auto g_y_0_0_0_x_xz_z_x = buffer_1000_pdpp[186];

    auto g_y_0_0_0_x_xz_z_y = buffer_1000_pdpp[187];

    auto g_y_0_0_0_x_xz_z_z = buffer_1000_pdpp[188];

    auto g_y_0_0_0_x_yy_x_x = buffer_1000_pdpp[189];

    auto g_y_0_0_0_x_yy_x_y = buffer_1000_pdpp[190];

    auto g_y_0_0_0_x_yy_x_z = buffer_1000_pdpp[191];

    auto g_y_0_0_0_x_yy_y_x = buffer_1000_pdpp[192];

    auto g_y_0_0_0_x_yy_y_y = buffer_1000_pdpp[193];

    auto g_y_0_0_0_x_yy_y_z = buffer_1000_pdpp[194];

    auto g_y_0_0_0_x_yy_z_x = buffer_1000_pdpp[195];

    auto g_y_0_0_0_x_yy_z_y = buffer_1000_pdpp[196];

    auto g_y_0_0_0_x_yy_z_z = buffer_1000_pdpp[197];

    auto g_y_0_0_0_x_yz_x_x = buffer_1000_pdpp[198];

    auto g_y_0_0_0_x_yz_x_y = buffer_1000_pdpp[199];

    auto g_y_0_0_0_x_yz_x_z = buffer_1000_pdpp[200];

    auto g_y_0_0_0_x_yz_y_x = buffer_1000_pdpp[201];

    auto g_y_0_0_0_x_yz_y_y = buffer_1000_pdpp[202];

    auto g_y_0_0_0_x_yz_y_z = buffer_1000_pdpp[203];

    auto g_y_0_0_0_x_yz_z_x = buffer_1000_pdpp[204];

    auto g_y_0_0_0_x_yz_z_y = buffer_1000_pdpp[205];

    auto g_y_0_0_0_x_yz_z_z = buffer_1000_pdpp[206];

    auto g_y_0_0_0_x_zz_x_x = buffer_1000_pdpp[207];

    auto g_y_0_0_0_x_zz_x_y = buffer_1000_pdpp[208];

    auto g_y_0_0_0_x_zz_x_z = buffer_1000_pdpp[209];

    auto g_y_0_0_0_x_zz_y_x = buffer_1000_pdpp[210];

    auto g_y_0_0_0_x_zz_y_y = buffer_1000_pdpp[211];

    auto g_y_0_0_0_x_zz_y_z = buffer_1000_pdpp[212];

    auto g_y_0_0_0_x_zz_z_x = buffer_1000_pdpp[213];

    auto g_y_0_0_0_x_zz_z_y = buffer_1000_pdpp[214];

    auto g_y_0_0_0_x_zz_z_z = buffer_1000_pdpp[215];

    auto g_y_0_0_0_y_xx_x_x = buffer_1000_pdpp[216];

    auto g_y_0_0_0_y_xx_x_y = buffer_1000_pdpp[217];

    auto g_y_0_0_0_y_xx_x_z = buffer_1000_pdpp[218];

    auto g_y_0_0_0_y_xx_y_x = buffer_1000_pdpp[219];

    auto g_y_0_0_0_y_xx_y_y = buffer_1000_pdpp[220];

    auto g_y_0_0_0_y_xx_y_z = buffer_1000_pdpp[221];

    auto g_y_0_0_0_y_xx_z_x = buffer_1000_pdpp[222];

    auto g_y_0_0_0_y_xx_z_y = buffer_1000_pdpp[223];

    auto g_y_0_0_0_y_xx_z_z = buffer_1000_pdpp[224];

    auto g_y_0_0_0_y_xy_x_x = buffer_1000_pdpp[225];

    auto g_y_0_0_0_y_xy_x_y = buffer_1000_pdpp[226];

    auto g_y_0_0_0_y_xy_x_z = buffer_1000_pdpp[227];

    auto g_y_0_0_0_y_xy_y_x = buffer_1000_pdpp[228];

    auto g_y_0_0_0_y_xy_y_y = buffer_1000_pdpp[229];

    auto g_y_0_0_0_y_xy_y_z = buffer_1000_pdpp[230];

    auto g_y_0_0_0_y_xy_z_x = buffer_1000_pdpp[231];

    auto g_y_0_0_0_y_xy_z_y = buffer_1000_pdpp[232];

    auto g_y_0_0_0_y_xy_z_z = buffer_1000_pdpp[233];

    auto g_y_0_0_0_y_xz_x_x = buffer_1000_pdpp[234];

    auto g_y_0_0_0_y_xz_x_y = buffer_1000_pdpp[235];

    auto g_y_0_0_0_y_xz_x_z = buffer_1000_pdpp[236];

    auto g_y_0_0_0_y_xz_y_x = buffer_1000_pdpp[237];

    auto g_y_0_0_0_y_xz_y_y = buffer_1000_pdpp[238];

    auto g_y_0_0_0_y_xz_y_z = buffer_1000_pdpp[239];

    auto g_y_0_0_0_y_xz_z_x = buffer_1000_pdpp[240];

    auto g_y_0_0_0_y_xz_z_y = buffer_1000_pdpp[241];

    auto g_y_0_0_0_y_xz_z_z = buffer_1000_pdpp[242];

    auto g_y_0_0_0_y_yy_x_x = buffer_1000_pdpp[243];

    auto g_y_0_0_0_y_yy_x_y = buffer_1000_pdpp[244];

    auto g_y_0_0_0_y_yy_x_z = buffer_1000_pdpp[245];

    auto g_y_0_0_0_y_yy_y_x = buffer_1000_pdpp[246];

    auto g_y_0_0_0_y_yy_y_y = buffer_1000_pdpp[247];

    auto g_y_0_0_0_y_yy_y_z = buffer_1000_pdpp[248];

    auto g_y_0_0_0_y_yy_z_x = buffer_1000_pdpp[249];

    auto g_y_0_0_0_y_yy_z_y = buffer_1000_pdpp[250];

    auto g_y_0_0_0_y_yy_z_z = buffer_1000_pdpp[251];

    auto g_y_0_0_0_y_yz_x_x = buffer_1000_pdpp[252];

    auto g_y_0_0_0_y_yz_x_y = buffer_1000_pdpp[253];

    auto g_y_0_0_0_y_yz_x_z = buffer_1000_pdpp[254];

    auto g_y_0_0_0_y_yz_y_x = buffer_1000_pdpp[255];

    auto g_y_0_0_0_y_yz_y_y = buffer_1000_pdpp[256];

    auto g_y_0_0_0_y_yz_y_z = buffer_1000_pdpp[257];

    auto g_y_0_0_0_y_yz_z_x = buffer_1000_pdpp[258];

    auto g_y_0_0_0_y_yz_z_y = buffer_1000_pdpp[259];

    auto g_y_0_0_0_y_yz_z_z = buffer_1000_pdpp[260];

    auto g_y_0_0_0_y_zz_x_x = buffer_1000_pdpp[261];

    auto g_y_0_0_0_y_zz_x_y = buffer_1000_pdpp[262];

    auto g_y_0_0_0_y_zz_x_z = buffer_1000_pdpp[263];

    auto g_y_0_0_0_y_zz_y_x = buffer_1000_pdpp[264];

    auto g_y_0_0_0_y_zz_y_y = buffer_1000_pdpp[265];

    auto g_y_0_0_0_y_zz_y_z = buffer_1000_pdpp[266];

    auto g_y_0_0_0_y_zz_z_x = buffer_1000_pdpp[267];

    auto g_y_0_0_0_y_zz_z_y = buffer_1000_pdpp[268];

    auto g_y_0_0_0_y_zz_z_z = buffer_1000_pdpp[269];

    auto g_y_0_0_0_z_xx_x_x = buffer_1000_pdpp[270];

    auto g_y_0_0_0_z_xx_x_y = buffer_1000_pdpp[271];

    auto g_y_0_0_0_z_xx_x_z = buffer_1000_pdpp[272];

    auto g_y_0_0_0_z_xx_y_x = buffer_1000_pdpp[273];

    auto g_y_0_0_0_z_xx_y_y = buffer_1000_pdpp[274];

    auto g_y_0_0_0_z_xx_y_z = buffer_1000_pdpp[275];

    auto g_y_0_0_0_z_xx_z_x = buffer_1000_pdpp[276];

    auto g_y_0_0_0_z_xx_z_y = buffer_1000_pdpp[277];

    auto g_y_0_0_0_z_xx_z_z = buffer_1000_pdpp[278];

    auto g_y_0_0_0_z_xy_x_x = buffer_1000_pdpp[279];

    auto g_y_0_0_0_z_xy_x_y = buffer_1000_pdpp[280];

    auto g_y_0_0_0_z_xy_x_z = buffer_1000_pdpp[281];

    auto g_y_0_0_0_z_xy_y_x = buffer_1000_pdpp[282];

    auto g_y_0_0_0_z_xy_y_y = buffer_1000_pdpp[283];

    auto g_y_0_0_0_z_xy_y_z = buffer_1000_pdpp[284];

    auto g_y_0_0_0_z_xy_z_x = buffer_1000_pdpp[285];

    auto g_y_0_0_0_z_xy_z_y = buffer_1000_pdpp[286];

    auto g_y_0_0_0_z_xy_z_z = buffer_1000_pdpp[287];

    auto g_y_0_0_0_z_xz_x_x = buffer_1000_pdpp[288];

    auto g_y_0_0_0_z_xz_x_y = buffer_1000_pdpp[289];

    auto g_y_0_0_0_z_xz_x_z = buffer_1000_pdpp[290];

    auto g_y_0_0_0_z_xz_y_x = buffer_1000_pdpp[291];

    auto g_y_0_0_0_z_xz_y_y = buffer_1000_pdpp[292];

    auto g_y_0_0_0_z_xz_y_z = buffer_1000_pdpp[293];

    auto g_y_0_0_0_z_xz_z_x = buffer_1000_pdpp[294];

    auto g_y_0_0_0_z_xz_z_y = buffer_1000_pdpp[295];

    auto g_y_0_0_0_z_xz_z_z = buffer_1000_pdpp[296];

    auto g_y_0_0_0_z_yy_x_x = buffer_1000_pdpp[297];

    auto g_y_0_0_0_z_yy_x_y = buffer_1000_pdpp[298];

    auto g_y_0_0_0_z_yy_x_z = buffer_1000_pdpp[299];

    auto g_y_0_0_0_z_yy_y_x = buffer_1000_pdpp[300];

    auto g_y_0_0_0_z_yy_y_y = buffer_1000_pdpp[301];

    auto g_y_0_0_0_z_yy_y_z = buffer_1000_pdpp[302];

    auto g_y_0_0_0_z_yy_z_x = buffer_1000_pdpp[303];

    auto g_y_0_0_0_z_yy_z_y = buffer_1000_pdpp[304];

    auto g_y_0_0_0_z_yy_z_z = buffer_1000_pdpp[305];

    auto g_y_0_0_0_z_yz_x_x = buffer_1000_pdpp[306];

    auto g_y_0_0_0_z_yz_x_y = buffer_1000_pdpp[307];

    auto g_y_0_0_0_z_yz_x_z = buffer_1000_pdpp[308];

    auto g_y_0_0_0_z_yz_y_x = buffer_1000_pdpp[309];

    auto g_y_0_0_0_z_yz_y_y = buffer_1000_pdpp[310];

    auto g_y_0_0_0_z_yz_y_z = buffer_1000_pdpp[311];

    auto g_y_0_0_0_z_yz_z_x = buffer_1000_pdpp[312];

    auto g_y_0_0_0_z_yz_z_y = buffer_1000_pdpp[313];

    auto g_y_0_0_0_z_yz_z_z = buffer_1000_pdpp[314];

    auto g_y_0_0_0_z_zz_x_x = buffer_1000_pdpp[315];

    auto g_y_0_0_0_z_zz_x_y = buffer_1000_pdpp[316];

    auto g_y_0_0_0_z_zz_x_z = buffer_1000_pdpp[317];

    auto g_y_0_0_0_z_zz_y_x = buffer_1000_pdpp[318];

    auto g_y_0_0_0_z_zz_y_y = buffer_1000_pdpp[319];

    auto g_y_0_0_0_z_zz_y_z = buffer_1000_pdpp[320];

    auto g_y_0_0_0_z_zz_z_x = buffer_1000_pdpp[321];

    auto g_y_0_0_0_z_zz_z_y = buffer_1000_pdpp[322];

    auto g_y_0_0_0_z_zz_z_z = buffer_1000_pdpp[323];

    auto g_z_0_0_0_x_xx_x_x = buffer_1000_pdpp[324];

    auto g_z_0_0_0_x_xx_x_y = buffer_1000_pdpp[325];

    auto g_z_0_0_0_x_xx_x_z = buffer_1000_pdpp[326];

    auto g_z_0_0_0_x_xx_y_x = buffer_1000_pdpp[327];

    auto g_z_0_0_0_x_xx_y_y = buffer_1000_pdpp[328];

    auto g_z_0_0_0_x_xx_y_z = buffer_1000_pdpp[329];

    auto g_z_0_0_0_x_xx_z_x = buffer_1000_pdpp[330];

    auto g_z_0_0_0_x_xx_z_y = buffer_1000_pdpp[331];

    auto g_z_0_0_0_x_xx_z_z = buffer_1000_pdpp[332];

    auto g_z_0_0_0_x_xy_x_x = buffer_1000_pdpp[333];

    auto g_z_0_0_0_x_xy_x_y = buffer_1000_pdpp[334];

    auto g_z_0_0_0_x_xy_x_z = buffer_1000_pdpp[335];

    auto g_z_0_0_0_x_xy_y_x = buffer_1000_pdpp[336];

    auto g_z_0_0_0_x_xy_y_y = buffer_1000_pdpp[337];

    auto g_z_0_0_0_x_xy_y_z = buffer_1000_pdpp[338];

    auto g_z_0_0_0_x_xy_z_x = buffer_1000_pdpp[339];

    auto g_z_0_0_0_x_xy_z_y = buffer_1000_pdpp[340];

    auto g_z_0_0_0_x_xy_z_z = buffer_1000_pdpp[341];

    auto g_z_0_0_0_x_xz_x_x = buffer_1000_pdpp[342];

    auto g_z_0_0_0_x_xz_x_y = buffer_1000_pdpp[343];

    auto g_z_0_0_0_x_xz_x_z = buffer_1000_pdpp[344];

    auto g_z_0_0_0_x_xz_y_x = buffer_1000_pdpp[345];

    auto g_z_0_0_0_x_xz_y_y = buffer_1000_pdpp[346];

    auto g_z_0_0_0_x_xz_y_z = buffer_1000_pdpp[347];

    auto g_z_0_0_0_x_xz_z_x = buffer_1000_pdpp[348];

    auto g_z_0_0_0_x_xz_z_y = buffer_1000_pdpp[349];

    auto g_z_0_0_0_x_xz_z_z = buffer_1000_pdpp[350];

    auto g_z_0_0_0_x_yy_x_x = buffer_1000_pdpp[351];

    auto g_z_0_0_0_x_yy_x_y = buffer_1000_pdpp[352];

    auto g_z_0_0_0_x_yy_x_z = buffer_1000_pdpp[353];

    auto g_z_0_0_0_x_yy_y_x = buffer_1000_pdpp[354];

    auto g_z_0_0_0_x_yy_y_y = buffer_1000_pdpp[355];

    auto g_z_0_0_0_x_yy_y_z = buffer_1000_pdpp[356];

    auto g_z_0_0_0_x_yy_z_x = buffer_1000_pdpp[357];

    auto g_z_0_0_0_x_yy_z_y = buffer_1000_pdpp[358];

    auto g_z_0_0_0_x_yy_z_z = buffer_1000_pdpp[359];

    auto g_z_0_0_0_x_yz_x_x = buffer_1000_pdpp[360];

    auto g_z_0_0_0_x_yz_x_y = buffer_1000_pdpp[361];

    auto g_z_0_0_0_x_yz_x_z = buffer_1000_pdpp[362];

    auto g_z_0_0_0_x_yz_y_x = buffer_1000_pdpp[363];

    auto g_z_0_0_0_x_yz_y_y = buffer_1000_pdpp[364];

    auto g_z_0_0_0_x_yz_y_z = buffer_1000_pdpp[365];

    auto g_z_0_0_0_x_yz_z_x = buffer_1000_pdpp[366];

    auto g_z_0_0_0_x_yz_z_y = buffer_1000_pdpp[367];

    auto g_z_0_0_0_x_yz_z_z = buffer_1000_pdpp[368];

    auto g_z_0_0_0_x_zz_x_x = buffer_1000_pdpp[369];

    auto g_z_0_0_0_x_zz_x_y = buffer_1000_pdpp[370];

    auto g_z_0_0_0_x_zz_x_z = buffer_1000_pdpp[371];

    auto g_z_0_0_0_x_zz_y_x = buffer_1000_pdpp[372];

    auto g_z_0_0_0_x_zz_y_y = buffer_1000_pdpp[373];

    auto g_z_0_0_0_x_zz_y_z = buffer_1000_pdpp[374];

    auto g_z_0_0_0_x_zz_z_x = buffer_1000_pdpp[375];

    auto g_z_0_0_0_x_zz_z_y = buffer_1000_pdpp[376];

    auto g_z_0_0_0_x_zz_z_z = buffer_1000_pdpp[377];

    auto g_z_0_0_0_y_xx_x_x = buffer_1000_pdpp[378];

    auto g_z_0_0_0_y_xx_x_y = buffer_1000_pdpp[379];

    auto g_z_0_0_0_y_xx_x_z = buffer_1000_pdpp[380];

    auto g_z_0_0_0_y_xx_y_x = buffer_1000_pdpp[381];

    auto g_z_0_0_0_y_xx_y_y = buffer_1000_pdpp[382];

    auto g_z_0_0_0_y_xx_y_z = buffer_1000_pdpp[383];

    auto g_z_0_0_0_y_xx_z_x = buffer_1000_pdpp[384];

    auto g_z_0_0_0_y_xx_z_y = buffer_1000_pdpp[385];

    auto g_z_0_0_0_y_xx_z_z = buffer_1000_pdpp[386];

    auto g_z_0_0_0_y_xy_x_x = buffer_1000_pdpp[387];

    auto g_z_0_0_0_y_xy_x_y = buffer_1000_pdpp[388];

    auto g_z_0_0_0_y_xy_x_z = buffer_1000_pdpp[389];

    auto g_z_0_0_0_y_xy_y_x = buffer_1000_pdpp[390];

    auto g_z_0_0_0_y_xy_y_y = buffer_1000_pdpp[391];

    auto g_z_0_0_0_y_xy_y_z = buffer_1000_pdpp[392];

    auto g_z_0_0_0_y_xy_z_x = buffer_1000_pdpp[393];

    auto g_z_0_0_0_y_xy_z_y = buffer_1000_pdpp[394];

    auto g_z_0_0_0_y_xy_z_z = buffer_1000_pdpp[395];

    auto g_z_0_0_0_y_xz_x_x = buffer_1000_pdpp[396];

    auto g_z_0_0_0_y_xz_x_y = buffer_1000_pdpp[397];

    auto g_z_0_0_0_y_xz_x_z = buffer_1000_pdpp[398];

    auto g_z_0_0_0_y_xz_y_x = buffer_1000_pdpp[399];

    auto g_z_0_0_0_y_xz_y_y = buffer_1000_pdpp[400];

    auto g_z_0_0_0_y_xz_y_z = buffer_1000_pdpp[401];

    auto g_z_0_0_0_y_xz_z_x = buffer_1000_pdpp[402];

    auto g_z_0_0_0_y_xz_z_y = buffer_1000_pdpp[403];

    auto g_z_0_0_0_y_xz_z_z = buffer_1000_pdpp[404];

    auto g_z_0_0_0_y_yy_x_x = buffer_1000_pdpp[405];

    auto g_z_0_0_0_y_yy_x_y = buffer_1000_pdpp[406];

    auto g_z_0_0_0_y_yy_x_z = buffer_1000_pdpp[407];

    auto g_z_0_0_0_y_yy_y_x = buffer_1000_pdpp[408];

    auto g_z_0_0_0_y_yy_y_y = buffer_1000_pdpp[409];

    auto g_z_0_0_0_y_yy_y_z = buffer_1000_pdpp[410];

    auto g_z_0_0_0_y_yy_z_x = buffer_1000_pdpp[411];

    auto g_z_0_0_0_y_yy_z_y = buffer_1000_pdpp[412];

    auto g_z_0_0_0_y_yy_z_z = buffer_1000_pdpp[413];

    auto g_z_0_0_0_y_yz_x_x = buffer_1000_pdpp[414];

    auto g_z_0_0_0_y_yz_x_y = buffer_1000_pdpp[415];

    auto g_z_0_0_0_y_yz_x_z = buffer_1000_pdpp[416];

    auto g_z_0_0_0_y_yz_y_x = buffer_1000_pdpp[417];

    auto g_z_0_0_0_y_yz_y_y = buffer_1000_pdpp[418];

    auto g_z_0_0_0_y_yz_y_z = buffer_1000_pdpp[419];

    auto g_z_0_0_0_y_yz_z_x = buffer_1000_pdpp[420];

    auto g_z_0_0_0_y_yz_z_y = buffer_1000_pdpp[421];

    auto g_z_0_0_0_y_yz_z_z = buffer_1000_pdpp[422];

    auto g_z_0_0_0_y_zz_x_x = buffer_1000_pdpp[423];

    auto g_z_0_0_0_y_zz_x_y = buffer_1000_pdpp[424];

    auto g_z_0_0_0_y_zz_x_z = buffer_1000_pdpp[425];

    auto g_z_0_0_0_y_zz_y_x = buffer_1000_pdpp[426];

    auto g_z_0_0_0_y_zz_y_y = buffer_1000_pdpp[427];

    auto g_z_0_0_0_y_zz_y_z = buffer_1000_pdpp[428];

    auto g_z_0_0_0_y_zz_z_x = buffer_1000_pdpp[429];

    auto g_z_0_0_0_y_zz_z_y = buffer_1000_pdpp[430];

    auto g_z_0_0_0_y_zz_z_z = buffer_1000_pdpp[431];

    auto g_z_0_0_0_z_xx_x_x = buffer_1000_pdpp[432];

    auto g_z_0_0_0_z_xx_x_y = buffer_1000_pdpp[433];

    auto g_z_0_0_0_z_xx_x_z = buffer_1000_pdpp[434];

    auto g_z_0_0_0_z_xx_y_x = buffer_1000_pdpp[435];

    auto g_z_0_0_0_z_xx_y_y = buffer_1000_pdpp[436];

    auto g_z_0_0_0_z_xx_y_z = buffer_1000_pdpp[437];

    auto g_z_0_0_0_z_xx_z_x = buffer_1000_pdpp[438];

    auto g_z_0_0_0_z_xx_z_y = buffer_1000_pdpp[439];

    auto g_z_0_0_0_z_xx_z_z = buffer_1000_pdpp[440];

    auto g_z_0_0_0_z_xy_x_x = buffer_1000_pdpp[441];

    auto g_z_0_0_0_z_xy_x_y = buffer_1000_pdpp[442];

    auto g_z_0_0_0_z_xy_x_z = buffer_1000_pdpp[443];

    auto g_z_0_0_0_z_xy_y_x = buffer_1000_pdpp[444];

    auto g_z_0_0_0_z_xy_y_y = buffer_1000_pdpp[445];

    auto g_z_0_0_0_z_xy_y_z = buffer_1000_pdpp[446];

    auto g_z_0_0_0_z_xy_z_x = buffer_1000_pdpp[447];

    auto g_z_0_0_0_z_xy_z_y = buffer_1000_pdpp[448];

    auto g_z_0_0_0_z_xy_z_z = buffer_1000_pdpp[449];

    auto g_z_0_0_0_z_xz_x_x = buffer_1000_pdpp[450];

    auto g_z_0_0_0_z_xz_x_y = buffer_1000_pdpp[451];

    auto g_z_0_0_0_z_xz_x_z = buffer_1000_pdpp[452];

    auto g_z_0_0_0_z_xz_y_x = buffer_1000_pdpp[453];

    auto g_z_0_0_0_z_xz_y_y = buffer_1000_pdpp[454];

    auto g_z_0_0_0_z_xz_y_z = buffer_1000_pdpp[455];

    auto g_z_0_0_0_z_xz_z_x = buffer_1000_pdpp[456];

    auto g_z_0_0_0_z_xz_z_y = buffer_1000_pdpp[457];

    auto g_z_0_0_0_z_xz_z_z = buffer_1000_pdpp[458];

    auto g_z_0_0_0_z_yy_x_x = buffer_1000_pdpp[459];

    auto g_z_0_0_0_z_yy_x_y = buffer_1000_pdpp[460];

    auto g_z_0_0_0_z_yy_x_z = buffer_1000_pdpp[461];

    auto g_z_0_0_0_z_yy_y_x = buffer_1000_pdpp[462];

    auto g_z_0_0_0_z_yy_y_y = buffer_1000_pdpp[463];

    auto g_z_0_0_0_z_yy_y_z = buffer_1000_pdpp[464];

    auto g_z_0_0_0_z_yy_z_x = buffer_1000_pdpp[465];

    auto g_z_0_0_0_z_yy_z_y = buffer_1000_pdpp[466];

    auto g_z_0_0_0_z_yy_z_z = buffer_1000_pdpp[467];

    auto g_z_0_0_0_z_yz_x_x = buffer_1000_pdpp[468];

    auto g_z_0_0_0_z_yz_x_y = buffer_1000_pdpp[469];

    auto g_z_0_0_0_z_yz_x_z = buffer_1000_pdpp[470];

    auto g_z_0_0_0_z_yz_y_x = buffer_1000_pdpp[471];

    auto g_z_0_0_0_z_yz_y_y = buffer_1000_pdpp[472];

    auto g_z_0_0_0_z_yz_y_z = buffer_1000_pdpp[473];

    auto g_z_0_0_0_z_yz_z_x = buffer_1000_pdpp[474];

    auto g_z_0_0_0_z_yz_z_y = buffer_1000_pdpp[475];

    auto g_z_0_0_0_z_yz_z_z = buffer_1000_pdpp[476];

    auto g_z_0_0_0_z_zz_x_x = buffer_1000_pdpp[477];

    auto g_z_0_0_0_z_zz_x_y = buffer_1000_pdpp[478];

    auto g_z_0_0_0_z_zz_x_z = buffer_1000_pdpp[479];

    auto g_z_0_0_0_z_zz_y_x = buffer_1000_pdpp[480];

    auto g_z_0_0_0_z_zz_y_y = buffer_1000_pdpp[481];

    auto g_z_0_0_0_z_zz_y_z = buffer_1000_pdpp[482];

    auto g_z_0_0_0_z_zz_z_x = buffer_1000_pdpp[483];

    auto g_z_0_0_0_z_zz_z_y = buffer_1000_pdpp[484];

    auto g_z_0_0_0_z_zz_z_z = buffer_1000_pdpp[485];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_x_0_0_0_x_xx_x_x, g_x_0_0_0_x_xx_x_y, g_x_0_0_0_x_xx_x_z, g_xx_xx_x_x, g_xx_xx_x_y, g_xx_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xx_x_x[i] = -g_0_xx_x_x[i] + 2.0 * g_xx_xx_x_x[i] * a_exp;

        g_x_0_0_0_x_xx_x_y[i] = -g_0_xx_x_y[i] + 2.0 * g_xx_xx_x_y[i] * a_exp;

        g_x_0_0_0_x_xx_x_z[i] = -g_0_xx_x_z[i] + 2.0 * g_xx_xx_x_z[i] * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_x_0_0_0_x_xx_y_x, g_x_0_0_0_x_xx_y_y, g_x_0_0_0_x_xx_y_z, g_xx_xx_y_x, g_xx_xx_y_y, g_xx_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xx_y_x[i] = -g_0_xx_y_x[i] + 2.0 * g_xx_xx_y_x[i] * a_exp;

        g_x_0_0_0_x_xx_y_y[i] = -g_0_xx_y_y[i] + 2.0 * g_xx_xx_y_y[i] * a_exp;

        g_x_0_0_0_x_xx_y_z[i] = -g_0_xx_y_z[i] + 2.0 * g_xx_xx_y_z[i] * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_x_0_0_0_x_xx_z_x, g_x_0_0_0_x_xx_z_y, g_x_0_0_0_x_xx_z_z, g_xx_xx_z_x, g_xx_xx_z_y, g_xx_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xx_z_x[i] = -g_0_xx_z_x[i] + 2.0 * g_xx_xx_z_x[i] * a_exp;

        g_x_0_0_0_x_xx_z_y[i] = -g_0_xx_z_y[i] + 2.0 * g_xx_xx_z_y[i] * a_exp;

        g_x_0_0_0_x_xx_z_z[i] = -g_0_xx_z_z[i] + 2.0 * g_xx_xx_z_z[i] * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_x_0_0_0_x_xy_x_x, g_x_0_0_0_x_xy_x_y, g_x_0_0_0_x_xy_x_z, g_xx_xy_x_x, g_xx_xy_x_y, g_xx_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xy_x_x[i] = -g_0_xy_x_x[i] + 2.0 * g_xx_xy_x_x[i] * a_exp;

        g_x_0_0_0_x_xy_x_y[i] = -g_0_xy_x_y[i] + 2.0 * g_xx_xy_x_y[i] * a_exp;

        g_x_0_0_0_x_xy_x_z[i] = -g_0_xy_x_z[i] + 2.0 * g_xx_xy_x_z[i] * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_x_0_0_0_x_xy_y_x, g_x_0_0_0_x_xy_y_y, g_x_0_0_0_x_xy_y_z, g_xx_xy_y_x, g_xx_xy_y_y, g_xx_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xy_y_x[i] = -g_0_xy_y_x[i] + 2.0 * g_xx_xy_y_x[i] * a_exp;

        g_x_0_0_0_x_xy_y_y[i] = -g_0_xy_y_y[i] + 2.0 * g_xx_xy_y_y[i] * a_exp;

        g_x_0_0_0_x_xy_y_z[i] = -g_0_xy_y_z[i] + 2.0 * g_xx_xy_y_z[i] * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_x_0_0_0_x_xy_z_x, g_x_0_0_0_x_xy_z_y, g_x_0_0_0_x_xy_z_z, g_xx_xy_z_x, g_xx_xy_z_y, g_xx_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xy_z_x[i] = -g_0_xy_z_x[i] + 2.0 * g_xx_xy_z_x[i] * a_exp;

        g_x_0_0_0_x_xy_z_y[i] = -g_0_xy_z_y[i] + 2.0 * g_xx_xy_z_y[i] * a_exp;

        g_x_0_0_0_x_xy_z_z[i] = -g_0_xy_z_z[i] + 2.0 * g_xx_xy_z_z[i] * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_x_0_0_0_x_xz_x_x, g_x_0_0_0_x_xz_x_y, g_x_0_0_0_x_xz_x_z, g_xx_xz_x_x, g_xx_xz_x_y, g_xx_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xz_x_x[i] = -g_0_xz_x_x[i] + 2.0 * g_xx_xz_x_x[i] * a_exp;

        g_x_0_0_0_x_xz_x_y[i] = -g_0_xz_x_y[i] + 2.0 * g_xx_xz_x_y[i] * a_exp;

        g_x_0_0_0_x_xz_x_z[i] = -g_0_xz_x_z[i] + 2.0 * g_xx_xz_x_z[i] * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_x_0_0_0_x_xz_y_x, g_x_0_0_0_x_xz_y_y, g_x_0_0_0_x_xz_y_z, g_xx_xz_y_x, g_xx_xz_y_y, g_xx_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xz_y_x[i] = -g_0_xz_y_x[i] + 2.0 * g_xx_xz_y_x[i] * a_exp;

        g_x_0_0_0_x_xz_y_y[i] = -g_0_xz_y_y[i] + 2.0 * g_xx_xz_y_y[i] * a_exp;

        g_x_0_0_0_x_xz_y_z[i] = -g_0_xz_y_z[i] + 2.0 * g_xx_xz_y_z[i] * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_x_0_0_0_x_xz_z_x, g_x_0_0_0_x_xz_z_y, g_x_0_0_0_x_xz_z_z, g_xx_xz_z_x, g_xx_xz_z_y, g_xx_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xz_z_x[i] = -g_0_xz_z_x[i] + 2.0 * g_xx_xz_z_x[i] * a_exp;

        g_x_0_0_0_x_xz_z_y[i] = -g_0_xz_z_y[i] + 2.0 * g_xx_xz_z_y[i] * a_exp;

        g_x_0_0_0_x_xz_z_z[i] = -g_0_xz_z_z[i] + 2.0 * g_xx_xz_z_z[i] * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_x_0_0_0_x_yy_x_x, g_x_0_0_0_x_yy_x_y, g_x_0_0_0_x_yy_x_z, g_xx_yy_x_x, g_xx_yy_x_y, g_xx_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yy_x_x[i] = -g_0_yy_x_x[i] + 2.0 * g_xx_yy_x_x[i] * a_exp;

        g_x_0_0_0_x_yy_x_y[i] = -g_0_yy_x_y[i] + 2.0 * g_xx_yy_x_y[i] * a_exp;

        g_x_0_0_0_x_yy_x_z[i] = -g_0_yy_x_z[i] + 2.0 * g_xx_yy_x_z[i] * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_x_0_0_0_x_yy_y_x, g_x_0_0_0_x_yy_y_y, g_x_0_0_0_x_yy_y_z, g_xx_yy_y_x, g_xx_yy_y_y, g_xx_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yy_y_x[i] = -g_0_yy_y_x[i] + 2.0 * g_xx_yy_y_x[i] * a_exp;

        g_x_0_0_0_x_yy_y_y[i] = -g_0_yy_y_y[i] + 2.0 * g_xx_yy_y_y[i] * a_exp;

        g_x_0_0_0_x_yy_y_z[i] = -g_0_yy_y_z[i] + 2.0 * g_xx_yy_y_z[i] * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_x_0_0_0_x_yy_z_x, g_x_0_0_0_x_yy_z_y, g_x_0_0_0_x_yy_z_z, g_xx_yy_z_x, g_xx_yy_z_y, g_xx_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yy_z_x[i] = -g_0_yy_z_x[i] + 2.0 * g_xx_yy_z_x[i] * a_exp;

        g_x_0_0_0_x_yy_z_y[i] = -g_0_yy_z_y[i] + 2.0 * g_xx_yy_z_y[i] * a_exp;

        g_x_0_0_0_x_yy_z_z[i] = -g_0_yy_z_z[i] + 2.0 * g_xx_yy_z_z[i] * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_x_0_0_0_x_yz_x_x, g_x_0_0_0_x_yz_x_y, g_x_0_0_0_x_yz_x_z, g_xx_yz_x_x, g_xx_yz_x_y, g_xx_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yz_x_x[i] = -g_0_yz_x_x[i] + 2.0 * g_xx_yz_x_x[i] * a_exp;

        g_x_0_0_0_x_yz_x_y[i] = -g_0_yz_x_y[i] + 2.0 * g_xx_yz_x_y[i] * a_exp;

        g_x_0_0_0_x_yz_x_z[i] = -g_0_yz_x_z[i] + 2.0 * g_xx_yz_x_z[i] * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_x_0_0_0_x_yz_y_x, g_x_0_0_0_x_yz_y_y, g_x_0_0_0_x_yz_y_z, g_xx_yz_y_x, g_xx_yz_y_y, g_xx_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yz_y_x[i] = -g_0_yz_y_x[i] + 2.0 * g_xx_yz_y_x[i] * a_exp;

        g_x_0_0_0_x_yz_y_y[i] = -g_0_yz_y_y[i] + 2.0 * g_xx_yz_y_y[i] * a_exp;

        g_x_0_0_0_x_yz_y_z[i] = -g_0_yz_y_z[i] + 2.0 * g_xx_yz_y_z[i] * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_x_0_0_0_x_yz_z_x, g_x_0_0_0_x_yz_z_y, g_x_0_0_0_x_yz_z_z, g_xx_yz_z_x, g_xx_yz_z_y, g_xx_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yz_z_x[i] = -g_0_yz_z_x[i] + 2.0 * g_xx_yz_z_x[i] * a_exp;

        g_x_0_0_0_x_yz_z_y[i] = -g_0_yz_z_y[i] + 2.0 * g_xx_yz_z_y[i] * a_exp;

        g_x_0_0_0_x_yz_z_z[i] = -g_0_yz_z_z[i] + 2.0 * g_xx_yz_z_z[i] * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_x_0_0_0_x_zz_x_x, g_x_0_0_0_x_zz_x_y, g_x_0_0_0_x_zz_x_z, g_xx_zz_x_x, g_xx_zz_x_y, g_xx_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_zz_x_x[i] = -g_0_zz_x_x[i] + 2.0 * g_xx_zz_x_x[i] * a_exp;

        g_x_0_0_0_x_zz_x_y[i] = -g_0_zz_x_y[i] + 2.0 * g_xx_zz_x_y[i] * a_exp;

        g_x_0_0_0_x_zz_x_z[i] = -g_0_zz_x_z[i] + 2.0 * g_xx_zz_x_z[i] * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_x_0_0_0_x_zz_y_x, g_x_0_0_0_x_zz_y_y, g_x_0_0_0_x_zz_y_z, g_xx_zz_y_x, g_xx_zz_y_y, g_xx_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_zz_y_x[i] = -g_0_zz_y_x[i] + 2.0 * g_xx_zz_y_x[i] * a_exp;

        g_x_0_0_0_x_zz_y_y[i] = -g_0_zz_y_y[i] + 2.0 * g_xx_zz_y_y[i] * a_exp;

        g_x_0_0_0_x_zz_y_z[i] = -g_0_zz_y_z[i] + 2.0 * g_xx_zz_y_z[i] * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_x_0_0_0_x_zz_z_x, g_x_0_0_0_x_zz_z_y, g_x_0_0_0_x_zz_z_z, g_xx_zz_z_x, g_xx_zz_z_y, g_xx_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_zz_z_x[i] = -g_0_zz_z_x[i] + 2.0 * g_xx_zz_z_x[i] * a_exp;

        g_x_0_0_0_x_zz_z_y[i] = -g_0_zz_z_y[i] + 2.0 * g_xx_zz_z_y[i] * a_exp;

        g_x_0_0_0_x_zz_z_z[i] = -g_0_zz_z_z[i] + 2.0 * g_xx_zz_z_z[i] * a_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_x_0_0_0_y_xx_x_x, g_x_0_0_0_y_xx_x_y, g_x_0_0_0_y_xx_x_z, g_xy_xx_x_x, g_xy_xx_x_y, g_xy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xx_x_x[i] = 2.0 * g_xy_xx_x_x[i] * a_exp;

        g_x_0_0_0_y_xx_x_y[i] = 2.0 * g_xy_xx_x_y[i] * a_exp;

        g_x_0_0_0_y_xx_x_z[i] = 2.0 * g_xy_xx_x_z[i] * a_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_x_0_0_0_y_xx_y_x, g_x_0_0_0_y_xx_y_y, g_x_0_0_0_y_xx_y_z, g_xy_xx_y_x, g_xy_xx_y_y, g_xy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xx_y_x[i] = 2.0 * g_xy_xx_y_x[i] * a_exp;

        g_x_0_0_0_y_xx_y_y[i] = 2.0 * g_xy_xx_y_y[i] * a_exp;

        g_x_0_0_0_y_xx_y_z[i] = 2.0 * g_xy_xx_y_z[i] * a_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_x_0_0_0_y_xx_z_x, g_x_0_0_0_y_xx_z_y, g_x_0_0_0_y_xx_z_z, g_xy_xx_z_x, g_xy_xx_z_y, g_xy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xx_z_x[i] = 2.0 * g_xy_xx_z_x[i] * a_exp;

        g_x_0_0_0_y_xx_z_y[i] = 2.0 * g_xy_xx_z_y[i] * a_exp;

        g_x_0_0_0_y_xx_z_z[i] = 2.0 * g_xy_xx_z_z[i] * a_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_x_0_0_0_y_xy_x_x, g_x_0_0_0_y_xy_x_y, g_x_0_0_0_y_xy_x_z, g_xy_xy_x_x, g_xy_xy_x_y, g_xy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xy_x_x[i] = 2.0 * g_xy_xy_x_x[i] * a_exp;

        g_x_0_0_0_y_xy_x_y[i] = 2.0 * g_xy_xy_x_y[i] * a_exp;

        g_x_0_0_0_y_xy_x_z[i] = 2.0 * g_xy_xy_x_z[i] * a_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_x_0_0_0_y_xy_y_x, g_x_0_0_0_y_xy_y_y, g_x_0_0_0_y_xy_y_z, g_xy_xy_y_x, g_xy_xy_y_y, g_xy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xy_y_x[i] = 2.0 * g_xy_xy_y_x[i] * a_exp;

        g_x_0_0_0_y_xy_y_y[i] = 2.0 * g_xy_xy_y_y[i] * a_exp;

        g_x_0_0_0_y_xy_y_z[i] = 2.0 * g_xy_xy_y_z[i] * a_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_x_0_0_0_y_xy_z_x, g_x_0_0_0_y_xy_z_y, g_x_0_0_0_y_xy_z_z, g_xy_xy_z_x, g_xy_xy_z_y, g_xy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xy_z_x[i] = 2.0 * g_xy_xy_z_x[i] * a_exp;

        g_x_0_0_0_y_xy_z_y[i] = 2.0 * g_xy_xy_z_y[i] * a_exp;

        g_x_0_0_0_y_xy_z_z[i] = 2.0 * g_xy_xy_z_z[i] * a_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_0_0_0_y_xz_x_x, g_x_0_0_0_y_xz_x_y, g_x_0_0_0_y_xz_x_z, g_xy_xz_x_x, g_xy_xz_x_y, g_xy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xz_x_x[i] = 2.0 * g_xy_xz_x_x[i] * a_exp;

        g_x_0_0_0_y_xz_x_y[i] = 2.0 * g_xy_xz_x_y[i] * a_exp;

        g_x_0_0_0_y_xz_x_z[i] = 2.0 * g_xy_xz_x_z[i] * a_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_0_0_0_y_xz_y_x, g_x_0_0_0_y_xz_y_y, g_x_0_0_0_y_xz_y_z, g_xy_xz_y_x, g_xy_xz_y_y, g_xy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xz_y_x[i] = 2.0 * g_xy_xz_y_x[i] * a_exp;

        g_x_0_0_0_y_xz_y_y[i] = 2.0 * g_xy_xz_y_y[i] * a_exp;

        g_x_0_0_0_y_xz_y_z[i] = 2.0 * g_xy_xz_y_z[i] * a_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_0_0_0_y_xz_z_x, g_x_0_0_0_y_xz_z_y, g_x_0_0_0_y_xz_z_z, g_xy_xz_z_x, g_xy_xz_z_y, g_xy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xz_z_x[i] = 2.0 * g_xy_xz_z_x[i] * a_exp;

        g_x_0_0_0_y_xz_z_y[i] = 2.0 * g_xy_xz_z_y[i] * a_exp;

        g_x_0_0_0_y_xz_z_z[i] = 2.0 * g_xy_xz_z_z[i] * a_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_x_0_0_0_y_yy_x_x, g_x_0_0_0_y_yy_x_y, g_x_0_0_0_y_yy_x_z, g_xy_yy_x_x, g_xy_yy_x_y, g_xy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yy_x_x[i] = 2.0 * g_xy_yy_x_x[i] * a_exp;

        g_x_0_0_0_y_yy_x_y[i] = 2.0 * g_xy_yy_x_y[i] * a_exp;

        g_x_0_0_0_y_yy_x_z[i] = 2.0 * g_xy_yy_x_z[i] * a_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_x_0_0_0_y_yy_y_x, g_x_0_0_0_y_yy_y_y, g_x_0_0_0_y_yy_y_z, g_xy_yy_y_x, g_xy_yy_y_y, g_xy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yy_y_x[i] = 2.0 * g_xy_yy_y_x[i] * a_exp;

        g_x_0_0_0_y_yy_y_y[i] = 2.0 * g_xy_yy_y_y[i] * a_exp;

        g_x_0_0_0_y_yy_y_z[i] = 2.0 * g_xy_yy_y_z[i] * a_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_x_0_0_0_y_yy_z_x, g_x_0_0_0_y_yy_z_y, g_x_0_0_0_y_yy_z_z, g_xy_yy_z_x, g_xy_yy_z_y, g_xy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yy_z_x[i] = 2.0 * g_xy_yy_z_x[i] * a_exp;

        g_x_0_0_0_y_yy_z_y[i] = 2.0 * g_xy_yy_z_y[i] * a_exp;

        g_x_0_0_0_y_yy_z_z[i] = 2.0 * g_xy_yy_z_z[i] * a_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_x_0_0_0_y_yz_x_x, g_x_0_0_0_y_yz_x_y, g_x_0_0_0_y_yz_x_z, g_xy_yz_x_x, g_xy_yz_x_y, g_xy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yz_x_x[i] = 2.0 * g_xy_yz_x_x[i] * a_exp;

        g_x_0_0_0_y_yz_x_y[i] = 2.0 * g_xy_yz_x_y[i] * a_exp;

        g_x_0_0_0_y_yz_x_z[i] = 2.0 * g_xy_yz_x_z[i] * a_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_x_0_0_0_y_yz_y_x, g_x_0_0_0_y_yz_y_y, g_x_0_0_0_y_yz_y_z, g_xy_yz_y_x, g_xy_yz_y_y, g_xy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yz_y_x[i] = 2.0 * g_xy_yz_y_x[i] * a_exp;

        g_x_0_0_0_y_yz_y_y[i] = 2.0 * g_xy_yz_y_y[i] * a_exp;

        g_x_0_0_0_y_yz_y_z[i] = 2.0 * g_xy_yz_y_z[i] * a_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_x_0_0_0_y_yz_z_x, g_x_0_0_0_y_yz_z_y, g_x_0_0_0_y_yz_z_z, g_xy_yz_z_x, g_xy_yz_z_y, g_xy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yz_z_x[i] = 2.0 * g_xy_yz_z_x[i] * a_exp;

        g_x_0_0_0_y_yz_z_y[i] = 2.0 * g_xy_yz_z_y[i] * a_exp;

        g_x_0_0_0_y_yz_z_z[i] = 2.0 * g_xy_yz_z_z[i] * a_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_x_0_0_0_y_zz_x_x, g_x_0_0_0_y_zz_x_y, g_x_0_0_0_y_zz_x_z, g_xy_zz_x_x, g_xy_zz_x_y, g_xy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_zz_x_x[i] = 2.0 * g_xy_zz_x_x[i] * a_exp;

        g_x_0_0_0_y_zz_x_y[i] = 2.0 * g_xy_zz_x_y[i] * a_exp;

        g_x_0_0_0_y_zz_x_z[i] = 2.0 * g_xy_zz_x_z[i] * a_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_x_0_0_0_y_zz_y_x, g_x_0_0_0_y_zz_y_y, g_x_0_0_0_y_zz_y_z, g_xy_zz_y_x, g_xy_zz_y_y, g_xy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_zz_y_x[i] = 2.0 * g_xy_zz_y_x[i] * a_exp;

        g_x_0_0_0_y_zz_y_y[i] = 2.0 * g_xy_zz_y_y[i] * a_exp;

        g_x_0_0_0_y_zz_y_z[i] = 2.0 * g_xy_zz_y_z[i] * a_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_x_0_0_0_y_zz_z_x, g_x_0_0_0_y_zz_z_y, g_x_0_0_0_y_zz_z_z, g_xy_zz_z_x, g_xy_zz_z_y, g_xy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_zz_z_x[i] = 2.0 * g_xy_zz_z_x[i] * a_exp;

        g_x_0_0_0_y_zz_z_y[i] = 2.0 * g_xy_zz_z_y[i] * a_exp;

        g_x_0_0_0_y_zz_z_z[i] = 2.0 * g_xy_zz_z_z[i] * a_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_x_0_0_0_z_xx_x_x, g_x_0_0_0_z_xx_x_y, g_x_0_0_0_z_xx_x_z, g_xz_xx_x_x, g_xz_xx_x_y, g_xz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xx_x_x[i] = 2.0 * g_xz_xx_x_x[i] * a_exp;

        g_x_0_0_0_z_xx_x_y[i] = 2.0 * g_xz_xx_x_y[i] * a_exp;

        g_x_0_0_0_z_xx_x_z[i] = 2.0 * g_xz_xx_x_z[i] * a_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_x_0_0_0_z_xx_y_x, g_x_0_0_0_z_xx_y_y, g_x_0_0_0_z_xx_y_z, g_xz_xx_y_x, g_xz_xx_y_y, g_xz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xx_y_x[i] = 2.0 * g_xz_xx_y_x[i] * a_exp;

        g_x_0_0_0_z_xx_y_y[i] = 2.0 * g_xz_xx_y_y[i] * a_exp;

        g_x_0_0_0_z_xx_y_z[i] = 2.0 * g_xz_xx_y_z[i] * a_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_x_0_0_0_z_xx_z_x, g_x_0_0_0_z_xx_z_y, g_x_0_0_0_z_xx_z_z, g_xz_xx_z_x, g_xz_xx_z_y, g_xz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xx_z_x[i] = 2.0 * g_xz_xx_z_x[i] * a_exp;

        g_x_0_0_0_z_xx_z_y[i] = 2.0 * g_xz_xx_z_y[i] * a_exp;

        g_x_0_0_0_z_xx_z_z[i] = 2.0 * g_xz_xx_z_z[i] * a_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_x_0_0_0_z_xy_x_x, g_x_0_0_0_z_xy_x_y, g_x_0_0_0_z_xy_x_z, g_xz_xy_x_x, g_xz_xy_x_y, g_xz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xy_x_x[i] = 2.0 * g_xz_xy_x_x[i] * a_exp;

        g_x_0_0_0_z_xy_x_y[i] = 2.0 * g_xz_xy_x_y[i] * a_exp;

        g_x_0_0_0_z_xy_x_z[i] = 2.0 * g_xz_xy_x_z[i] * a_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_x_0_0_0_z_xy_y_x, g_x_0_0_0_z_xy_y_y, g_x_0_0_0_z_xy_y_z, g_xz_xy_y_x, g_xz_xy_y_y, g_xz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xy_y_x[i] = 2.0 * g_xz_xy_y_x[i] * a_exp;

        g_x_0_0_0_z_xy_y_y[i] = 2.0 * g_xz_xy_y_y[i] * a_exp;

        g_x_0_0_0_z_xy_y_z[i] = 2.0 * g_xz_xy_y_z[i] * a_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_x_0_0_0_z_xy_z_x, g_x_0_0_0_z_xy_z_y, g_x_0_0_0_z_xy_z_z, g_xz_xy_z_x, g_xz_xy_z_y, g_xz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xy_z_x[i] = 2.0 * g_xz_xy_z_x[i] * a_exp;

        g_x_0_0_0_z_xy_z_y[i] = 2.0 * g_xz_xy_z_y[i] * a_exp;

        g_x_0_0_0_z_xy_z_z[i] = 2.0 * g_xz_xy_z_z[i] * a_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_x_0_0_0_z_xz_x_x, g_x_0_0_0_z_xz_x_y, g_x_0_0_0_z_xz_x_z, g_xz_xz_x_x, g_xz_xz_x_y, g_xz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xz_x_x[i] = 2.0 * g_xz_xz_x_x[i] * a_exp;

        g_x_0_0_0_z_xz_x_y[i] = 2.0 * g_xz_xz_x_y[i] * a_exp;

        g_x_0_0_0_z_xz_x_z[i] = 2.0 * g_xz_xz_x_z[i] * a_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_x_0_0_0_z_xz_y_x, g_x_0_0_0_z_xz_y_y, g_x_0_0_0_z_xz_y_z, g_xz_xz_y_x, g_xz_xz_y_y, g_xz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xz_y_x[i] = 2.0 * g_xz_xz_y_x[i] * a_exp;

        g_x_0_0_0_z_xz_y_y[i] = 2.0 * g_xz_xz_y_y[i] * a_exp;

        g_x_0_0_0_z_xz_y_z[i] = 2.0 * g_xz_xz_y_z[i] * a_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_x_0_0_0_z_xz_z_x, g_x_0_0_0_z_xz_z_y, g_x_0_0_0_z_xz_z_z, g_xz_xz_z_x, g_xz_xz_z_y, g_xz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xz_z_x[i] = 2.0 * g_xz_xz_z_x[i] * a_exp;

        g_x_0_0_0_z_xz_z_y[i] = 2.0 * g_xz_xz_z_y[i] * a_exp;

        g_x_0_0_0_z_xz_z_z[i] = 2.0 * g_xz_xz_z_z[i] * a_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_x_0_0_0_z_yy_x_x, g_x_0_0_0_z_yy_x_y, g_x_0_0_0_z_yy_x_z, g_xz_yy_x_x, g_xz_yy_x_y, g_xz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yy_x_x[i] = 2.0 * g_xz_yy_x_x[i] * a_exp;

        g_x_0_0_0_z_yy_x_y[i] = 2.0 * g_xz_yy_x_y[i] * a_exp;

        g_x_0_0_0_z_yy_x_z[i] = 2.0 * g_xz_yy_x_z[i] * a_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_x_0_0_0_z_yy_y_x, g_x_0_0_0_z_yy_y_y, g_x_0_0_0_z_yy_y_z, g_xz_yy_y_x, g_xz_yy_y_y, g_xz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yy_y_x[i] = 2.0 * g_xz_yy_y_x[i] * a_exp;

        g_x_0_0_0_z_yy_y_y[i] = 2.0 * g_xz_yy_y_y[i] * a_exp;

        g_x_0_0_0_z_yy_y_z[i] = 2.0 * g_xz_yy_y_z[i] * a_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_x_0_0_0_z_yy_z_x, g_x_0_0_0_z_yy_z_y, g_x_0_0_0_z_yy_z_z, g_xz_yy_z_x, g_xz_yy_z_y, g_xz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yy_z_x[i] = 2.0 * g_xz_yy_z_x[i] * a_exp;

        g_x_0_0_0_z_yy_z_y[i] = 2.0 * g_xz_yy_z_y[i] * a_exp;

        g_x_0_0_0_z_yy_z_z[i] = 2.0 * g_xz_yy_z_z[i] * a_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_x_0_0_0_z_yz_x_x, g_x_0_0_0_z_yz_x_y, g_x_0_0_0_z_yz_x_z, g_xz_yz_x_x, g_xz_yz_x_y, g_xz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yz_x_x[i] = 2.0 * g_xz_yz_x_x[i] * a_exp;

        g_x_0_0_0_z_yz_x_y[i] = 2.0 * g_xz_yz_x_y[i] * a_exp;

        g_x_0_0_0_z_yz_x_z[i] = 2.0 * g_xz_yz_x_z[i] * a_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_x_0_0_0_z_yz_y_x, g_x_0_0_0_z_yz_y_y, g_x_0_0_0_z_yz_y_z, g_xz_yz_y_x, g_xz_yz_y_y, g_xz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yz_y_x[i] = 2.0 * g_xz_yz_y_x[i] * a_exp;

        g_x_0_0_0_z_yz_y_y[i] = 2.0 * g_xz_yz_y_y[i] * a_exp;

        g_x_0_0_0_z_yz_y_z[i] = 2.0 * g_xz_yz_y_z[i] * a_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_x_0_0_0_z_yz_z_x, g_x_0_0_0_z_yz_z_y, g_x_0_0_0_z_yz_z_z, g_xz_yz_z_x, g_xz_yz_z_y, g_xz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yz_z_x[i] = 2.0 * g_xz_yz_z_x[i] * a_exp;

        g_x_0_0_0_z_yz_z_y[i] = 2.0 * g_xz_yz_z_y[i] * a_exp;

        g_x_0_0_0_z_yz_z_z[i] = 2.0 * g_xz_yz_z_z[i] * a_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_x_0_0_0_z_zz_x_x, g_x_0_0_0_z_zz_x_y, g_x_0_0_0_z_zz_x_z, g_xz_zz_x_x, g_xz_zz_x_y, g_xz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_zz_x_x[i] = 2.0 * g_xz_zz_x_x[i] * a_exp;

        g_x_0_0_0_z_zz_x_y[i] = 2.0 * g_xz_zz_x_y[i] * a_exp;

        g_x_0_0_0_z_zz_x_z[i] = 2.0 * g_xz_zz_x_z[i] * a_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_x_0_0_0_z_zz_y_x, g_x_0_0_0_z_zz_y_y, g_x_0_0_0_z_zz_y_z, g_xz_zz_y_x, g_xz_zz_y_y, g_xz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_zz_y_x[i] = 2.0 * g_xz_zz_y_x[i] * a_exp;

        g_x_0_0_0_z_zz_y_y[i] = 2.0 * g_xz_zz_y_y[i] * a_exp;

        g_x_0_0_0_z_zz_y_z[i] = 2.0 * g_xz_zz_y_z[i] * a_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_x_0_0_0_z_zz_z_x, g_x_0_0_0_z_zz_z_y, g_x_0_0_0_z_zz_z_z, g_xz_zz_z_x, g_xz_zz_z_y, g_xz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_zz_z_x[i] = 2.0 * g_xz_zz_z_x[i] * a_exp;

        g_x_0_0_0_z_zz_z_y[i] = 2.0 * g_xz_zz_z_y[i] * a_exp;

        g_x_0_0_0_z_zz_z_z[i] = 2.0 * g_xz_zz_z_z[i] * a_exp;
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_xy_xx_x_x, g_xy_xx_x_y, g_xy_xx_x_z, g_y_0_0_0_x_xx_x_x, g_y_0_0_0_x_xx_x_y, g_y_0_0_0_x_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xx_x_x[i] = 2.0 * g_xy_xx_x_x[i] * a_exp;

        g_y_0_0_0_x_xx_x_y[i] = 2.0 * g_xy_xx_x_y[i] * a_exp;

        g_y_0_0_0_x_xx_x_z[i] = 2.0 * g_xy_xx_x_z[i] * a_exp;
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_xy_xx_y_x, g_xy_xx_y_y, g_xy_xx_y_z, g_y_0_0_0_x_xx_y_x, g_y_0_0_0_x_xx_y_y, g_y_0_0_0_x_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xx_y_x[i] = 2.0 * g_xy_xx_y_x[i] * a_exp;

        g_y_0_0_0_x_xx_y_y[i] = 2.0 * g_xy_xx_y_y[i] * a_exp;

        g_y_0_0_0_x_xx_y_z[i] = 2.0 * g_xy_xx_y_z[i] * a_exp;
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_xy_xx_z_x, g_xy_xx_z_y, g_xy_xx_z_z, g_y_0_0_0_x_xx_z_x, g_y_0_0_0_x_xx_z_y, g_y_0_0_0_x_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xx_z_x[i] = 2.0 * g_xy_xx_z_x[i] * a_exp;

        g_y_0_0_0_x_xx_z_y[i] = 2.0 * g_xy_xx_z_y[i] * a_exp;

        g_y_0_0_0_x_xx_z_z[i] = 2.0 * g_xy_xx_z_z[i] * a_exp;
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_xy_xy_x_x, g_xy_xy_x_y, g_xy_xy_x_z, g_y_0_0_0_x_xy_x_x, g_y_0_0_0_x_xy_x_y, g_y_0_0_0_x_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xy_x_x[i] = 2.0 * g_xy_xy_x_x[i] * a_exp;

        g_y_0_0_0_x_xy_x_y[i] = 2.0 * g_xy_xy_x_y[i] * a_exp;

        g_y_0_0_0_x_xy_x_z[i] = 2.0 * g_xy_xy_x_z[i] * a_exp;
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_xy_xy_y_x, g_xy_xy_y_y, g_xy_xy_y_z, g_y_0_0_0_x_xy_y_x, g_y_0_0_0_x_xy_y_y, g_y_0_0_0_x_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xy_y_x[i] = 2.0 * g_xy_xy_y_x[i] * a_exp;

        g_y_0_0_0_x_xy_y_y[i] = 2.0 * g_xy_xy_y_y[i] * a_exp;

        g_y_0_0_0_x_xy_y_z[i] = 2.0 * g_xy_xy_y_z[i] * a_exp;
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_xy_xy_z_x, g_xy_xy_z_y, g_xy_xy_z_z, g_y_0_0_0_x_xy_z_x, g_y_0_0_0_x_xy_z_y, g_y_0_0_0_x_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xy_z_x[i] = 2.0 * g_xy_xy_z_x[i] * a_exp;

        g_y_0_0_0_x_xy_z_y[i] = 2.0 * g_xy_xy_z_y[i] * a_exp;

        g_y_0_0_0_x_xy_z_z[i] = 2.0 * g_xy_xy_z_z[i] * a_exp;
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_xy_xz_x_x, g_xy_xz_x_y, g_xy_xz_x_z, g_y_0_0_0_x_xz_x_x, g_y_0_0_0_x_xz_x_y, g_y_0_0_0_x_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xz_x_x[i] = 2.0 * g_xy_xz_x_x[i] * a_exp;

        g_y_0_0_0_x_xz_x_y[i] = 2.0 * g_xy_xz_x_y[i] * a_exp;

        g_y_0_0_0_x_xz_x_z[i] = 2.0 * g_xy_xz_x_z[i] * a_exp;
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_xy_xz_y_x, g_xy_xz_y_y, g_xy_xz_y_z, g_y_0_0_0_x_xz_y_x, g_y_0_0_0_x_xz_y_y, g_y_0_0_0_x_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xz_y_x[i] = 2.0 * g_xy_xz_y_x[i] * a_exp;

        g_y_0_0_0_x_xz_y_y[i] = 2.0 * g_xy_xz_y_y[i] * a_exp;

        g_y_0_0_0_x_xz_y_z[i] = 2.0 * g_xy_xz_y_z[i] * a_exp;
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_xy_xz_z_x, g_xy_xz_z_y, g_xy_xz_z_z, g_y_0_0_0_x_xz_z_x, g_y_0_0_0_x_xz_z_y, g_y_0_0_0_x_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xz_z_x[i] = 2.0 * g_xy_xz_z_x[i] * a_exp;

        g_y_0_0_0_x_xz_z_y[i] = 2.0 * g_xy_xz_z_y[i] * a_exp;

        g_y_0_0_0_x_xz_z_z[i] = 2.0 * g_xy_xz_z_z[i] * a_exp;
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_xy_yy_x_x, g_xy_yy_x_y, g_xy_yy_x_z, g_y_0_0_0_x_yy_x_x, g_y_0_0_0_x_yy_x_y, g_y_0_0_0_x_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yy_x_x[i] = 2.0 * g_xy_yy_x_x[i] * a_exp;

        g_y_0_0_0_x_yy_x_y[i] = 2.0 * g_xy_yy_x_y[i] * a_exp;

        g_y_0_0_0_x_yy_x_z[i] = 2.0 * g_xy_yy_x_z[i] * a_exp;
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_xy_yy_y_x, g_xy_yy_y_y, g_xy_yy_y_z, g_y_0_0_0_x_yy_y_x, g_y_0_0_0_x_yy_y_y, g_y_0_0_0_x_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yy_y_x[i] = 2.0 * g_xy_yy_y_x[i] * a_exp;

        g_y_0_0_0_x_yy_y_y[i] = 2.0 * g_xy_yy_y_y[i] * a_exp;

        g_y_0_0_0_x_yy_y_z[i] = 2.0 * g_xy_yy_y_z[i] * a_exp;
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_xy_yy_z_x, g_xy_yy_z_y, g_xy_yy_z_z, g_y_0_0_0_x_yy_z_x, g_y_0_0_0_x_yy_z_y, g_y_0_0_0_x_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yy_z_x[i] = 2.0 * g_xy_yy_z_x[i] * a_exp;

        g_y_0_0_0_x_yy_z_y[i] = 2.0 * g_xy_yy_z_y[i] * a_exp;

        g_y_0_0_0_x_yy_z_z[i] = 2.0 * g_xy_yy_z_z[i] * a_exp;
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_xy_yz_x_x, g_xy_yz_x_y, g_xy_yz_x_z, g_y_0_0_0_x_yz_x_x, g_y_0_0_0_x_yz_x_y, g_y_0_0_0_x_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yz_x_x[i] = 2.0 * g_xy_yz_x_x[i] * a_exp;

        g_y_0_0_0_x_yz_x_y[i] = 2.0 * g_xy_yz_x_y[i] * a_exp;

        g_y_0_0_0_x_yz_x_z[i] = 2.0 * g_xy_yz_x_z[i] * a_exp;
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_xy_yz_y_x, g_xy_yz_y_y, g_xy_yz_y_z, g_y_0_0_0_x_yz_y_x, g_y_0_0_0_x_yz_y_y, g_y_0_0_0_x_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yz_y_x[i] = 2.0 * g_xy_yz_y_x[i] * a_exp;

        g_y_0_0_0_x_yz_y_y[i] = 2.0 * g_xy_yz_y_y[i] * a_exp;

        g_y_0_0_0_x_yz_y_z[i] = 2.0 * g_xy_yz_y_z[i] * a_exp;
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_xy_yz_z_x, g_xy_yz_z_y, g_xy_yz_z_z, g_y_0_0_0_x_yz_z_x, g_y_0_0_0_x_yz_z_y, g_y_0_0_0_x_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yz_z_x[i] = 2.0 * g_xy_yz_z_x[i] * a_exp;

        g_y_0_0_0_x_yz_z_y[i] = 2.0 * g_xy_yz_z_y[i] * a_exp;

        g_y_0_0_0_x_yz_z_z[i] = 2.0 * g_xy_yz_z_z[i] * a_exp;
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_xy_zz_x_x, g_xy_zz_x_y, g_xy_zz_x_z, g_y_0_0_0_x_zz_x_x, g_y_0_0_0_x_zz_x_y, g_y_0_0_0_x_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_zz_x_x[i] = 2.0 * g_xy_zz_x_x[i] * a_exp;

        g_y_0_0_0_x_zz_x_y[i] = 2.0 * g_xy_zz_x_y[i] * a_exp;

        g_y_0_0_0_x_zz_x_z[i] = 2.0 * g_xy_zz_x_z[i] * a_exp;
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_xy_zz_y_x, g_xy_zz_y_y, g_xy_zz_y_z, g_y_0_0_0_x_zz_y_x, g_y_0_0_0_x_zz_y_y, g_y_0_0_0_x_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_zz_y_x[i] = 2.0 * g_xy_zz_y_x[i] * a_exp;

        g_y_0_0_0_x_zz_y_y[i] = 2.0 * g_xy_zz_y_y[i] * a_exp;

        g_y_0_0_0_x_zz_y_z[i] = 2.0 * g_xy_zz_y_z[i] * a_exp;
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_xy_zz_z_x, g_xy_zz_z_y, g_xy_zz_z_z, g_y_0_0_0_x_zz_z_x, g_y_0_0_0_x_zz_z_y, g_y_0_0_0_x_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_zz_z_x[i] = 2.0 * g_xy_zz_z_x[i] * a_exp;

        g_y_0_0_0_x_zz_z_y[i] = 2.0 * g_xy_zz_z_y[i] * a_exp;

        g_y_0_0_0_x_zz_z_z[i] = 2.0 * g_xy_zz_z_z[i] * a_exp;
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_y_0_0_0_y_xx_x_x, g_y_0_0_0_y_xx_x_y, g_y_0_0_0_y_xx_x_z, g_yy_xx_x_x, g_yy_xx_x_y, g_yy_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xx_x_x[i] = -g_0_xx_x_x[i] + 2.0 * g_yy_xx_x_x[i] * a_exp;

        g_y_0_0_0_y_xx_x_y[i] = -g_0_xx_x_y[i] + 2.0 * g_yy_xx_x_y[i] * a_exp;

        g_y_0_0_0_y_xx_x_z[i] = -g_0_xx_x_z[i] + 2.0 * g_yy_xx_x_z[i] * a_exp;
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_y_0_0_0_y_xx_y_x, g_y_0_0_0_y_xx_y_y, g_y_0_0_0_y_xx_y_z, g_yy_xx_y_x, g_yy_xx_y_y, g_yy_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xx_y_x[i] = -g_0_xx_y_x[i] + 2.0 * g_yy_xx_y_x[i] * a_exp;

        g_y_0_0_0_y_xx_y_y[i] = -g_0_xx_y_y[i] + 2.0 * g_yy_xx_y_y[i] * a_exp;

        g_y_0_0_0_y_xx_y_z[i] = -g_0_xx_y_z[i] + 2.0 * g_yy_xx_y_z[i] * a_exp;
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_y_0_0_0_y_xx_z_x, g_y_0_0_0_y_xx_z_y, g_y_0_0_0_y_xx_z_z, g_yy_xx_z_x, g_yy_xx_z_y, g_yy_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xx_z_x[i] = -g_0_xx_z_x[i] + 2.0 * g_yy_xx_z_x[i] * a_exp;

        g_y_0_0_0_y_xx_z_y[i] = -g_0_xx_z_y[i] + 2.0 * g_yy_xx_z_y[i] * a_exp;

        g_y_0_0_0_y_xx_z_z[i] = -g_0_xx_z_z[i] + 2.0 * g_yy_xx_z_z[i] * a_exp;
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_y_0_0_0_y_xy_x_x, g_y_0_0_0_y_xy_x_y, g_y_0_0_0_y_xy_x_z, g_yy_xy_x_x, g_yy_xy_x_y, g_yy_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xy_x_x[i] = -g_0_xy_x_x[i] + 2.0 * g_yy_xy_x_x[i] * a_exp;

        g_y_0_0_0_y_xy_x_y[i] = -g_0_xy_x_y[i] + 2.0 * g_yy_xy_x_y[i] * a_exp;

        g_y_0_0_0_y_xy_x_z[i] = -g_0_xy_x_z[i] + 2.0 * g_yy_xy_x_z[i] * a_exp;
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_y_0_0_0_y_xy_y_x, g_y_0_0_0_y_xy_y_y, g_y_0_0_0_y_xy_y_z, g_yy_xy_y_x, g_yy_xy_y_y, g_yy_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xy_y_x[i] = -g_0_xy_y_x[i] + 2.0 * g_yy_xy_y_x[i] * a_exp;

        g_y_0_0_0_y_xy_y_y[i] = -g_0_xy_y_y[i] + 2.0 * g_yy_xy_y_y[i] * a_exp;

        g_y_0_0_0_y_xy_y_z[i] = -g_0_xy_y_z[i] + 2.0 * g_yy_xy_y_z[i] * a_exp;
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_y_0_0_0_y_xy_z_x, g_y_0_0_0_y_xy_z_y, g_y_0_0_0_y_xy_z_z, g_yy_xy_z_x, g_yy_xy_z_y, g_yy_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xy_z_x[i] = -g_0_xy_z_x[i] + 2.0 * g_yy_xy_z_x[i] * a_exp;

        g_y_0_0_0_y_xy_z_y[i] = -g_0_xy_z_y[i] + 2.0 * g_yy_xy_z_y[i] * a_exp;

        g_y_0_0_0_y_xy_z_z[i] = -g_0_xy_z_z[i] + 2.0 * g_yy_xy_z_z[i] * a_exp;
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_y_0_0_0_y_xz_x_x, g_y_0_0_0_y_xz_x_y, g_y_0_0_0_y_xz_x_z, g_yy_xz_x_x, g_yy_xz_x_y, g_yy_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xz_x_x[i] = -g_0_xz_x_x[i] + 2.0 * g_yy_xz_x_x[i] * a_exp;

        g_y_0_0_0_y_xz_x_y[i] = -g_0_xz_x_y[i] + 2.0 * g_yy_xz_x_y[i] * a_exp;

        g_y_0_0_0_y_xz_x_z[i] = -g_0_xz_x_z[i] + 2.0 * g_yy_xz_x_z[i] * a_exp;
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_y_0_0_0_y_xz_y_x, g_y_0_0_0_y_xz_y_y, g_y_0_0_0_y_xz_y_z, g_yy_xz_y_x, g_yy_xz_y_y, g_yy_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xz_y_x[i] = -g_0_xz_y_x[i] + 2.0 * g_yy_xz_y_x[i] * a_exp;

        g_y_0_0_0_y_xz_y_y[i] = -g_0_xz_y_y[i] + 2.0 * g_yy_xz_y_y[i] * a_exp;

        g_y_0_0_0_y_xz_y_z[i] = -g_0_xz_y_z[i] + 2.0 * g_yy_xz_y_z[i] * a_exp;
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_y_0_0_0_y_xz_z_x, g_y_0_0_0_y_xz_z_y, g_y_0_0_0_y_xz_z_z, g_yy_xz_z_x, g_yy_xz_z_y, g_yy_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xz_z_x[i] = -g_0_xz_z_x[i] + 2.0 * g_yy_xz_z_x[i] * a_exp;

        g_y_0_0_0_y_xz_z_y[i] = -g_0_xz_z_y[i] + 2.0 * g_yy_xz_z_y[i] * a_exp;

        g_y_0_0_0_y_xz_z_z[i] = -g_0_xz_z_z[i] + 2.0 * g_yy_xz_z_z[i] * a_exp;
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_y_0_0_0_y_yy_x_x, g_y_0_0_0_y_yy_x_y, g_y_0_0_0_y_yy_x_z, g_yy_yy_x_x, g_yy_yy_x_y, g_yy_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yy_x_x[i] = -g_0_yy_x_x[i] + 2.0 * g_yy_yy_x_x[i] * a_exp;

        g_y_0_0_0_y_yy_x_y[i] = -g_0_yy_x_y[i] + 2.0 * g_yy_yy_x_y[i] * a_exp;

        g_y_0_0_0_y_yy_x_z[i] = -g_0_yy_x_z[i] + 2.0 * g_yy_yy_x_z[i] * a_exp;
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_y_0_0_0_y_yy_y_x, g_y_0_0_0_y_yy_y_y, g_y_0_0_0_y_yy_y_z, g_yy_yy_y_x, g_yy_yy_y_y, g_yy_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yy_y_x[i] = -g_0_yy_y_x[i] + 2.0 * g_yy_yy_y_x[i] * a_exp;

        g_y_0_0_0_y_yy_y_y[i] = -g_0_yy_y_y[i] + 2.0 * g_yy_yy_y_y[i] * a_exp;

        g_y_0_0_0_y_yy_y_z[i] = -g_0_yy_y_z[i] + 2.0 * g_yy_yy_y_z[i] * a_exp;
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_y_0_0_0_y_yy_z_x, g_y_0_0_0_y_yy_z_y, g_y_0_0_0_y_yy_z_z, g_yy_yy_z_x, g_yy_yy_z_y, g_yy_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yy_z_x[i] = -g_0_yy_z_x[i] + 2.0 * g_yy_yy_z_x[i] * a_exp;

        g_y_0_0_0_y_yy_z_y[i] = -g_0_yy_z_y[i] + 2.0 * g_yy_yy_z_y[i] * a_exp;

        g_y_0_0_0_y_yy_z_z[i] = -g_0_yy_z_z[i] + 2.0 * g_yy_yy_z_z[i] * a_exp;
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_y_0_0_0_y_yz_x_x, g_y_0_0_0_y_yz_x_y, g_y_0_0_0_y_yz_x_z, g_yy_yz_x_x, g_yy_yz_x_y, g_yy_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yz_x_x[i] = -g_0_yz_x_x[i] + 2.0 * g_yy_yz_x_x[i] * a_exp;

        g_y_0_0_0_y_yz_x_y[i] = -g_0_yz_x_y[i] + 2.0 * g_yy_yz_x_y[i] * a_exp;

        g_y_0_0_0_y_yz_x_z[i] = -g_0_yz_x_z[i] + 2.0 * g_yy_yz_x_z[i] * a_exp;
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_y_0_0_0_y_yz_y_x, g_y_0_0_0_y_yz_y_y, g_y_0_0_0_y_yz_y_z, g_yy_yz_y_x, g_yy_yz_y_y, g_yy_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yz_y_x[i] = -g_0_yz_y_x[i] + 2.0 * g_yy_yz_y_x[i] * a_exp;

        g_y_0_0_0_y_yz_y_y[i] = -g_0_yz_y_y[i] + 2.0 * g_yy_yz_y_y[i] * a_exp;

        g_y_0_0_0_y_yz_y_z[i] = -g_0_yz_y_z[i] + 2.0 * g_yy_yz_y_z[i] * a_exp;
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_y_0_0_0_y_yz_z_x, g_y_0_0_0_y_yz_z_y, g_y_0_0_0_y_yz_z_z, g_yy_yz_z_x, g_yy_yz_z_y, g_yy_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yz_z_x[i] = -g_0_yz_z_x[i] + 2.0 * g_yy_yz_z_x[i] * a_exp;

        g_y_0_0_0_y_yz_z_y[i] = -g_0_yz_z_y[i] + 2.0 * g_yy_yz_z_y[i] * a_exp;

        g_y_0_0_0_y_yz_z_z[i] = -g_0_yz_z_z[i] + 2.0 * g_yy_yz_z_z[i] * a_exp;
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_y_0_0_0_y_zz_x_x, g_y_0_0_0_y_zz_x_y, g_y_0_0_0_y_zz_x_z, g_yy_zz_x_x, g_yy_zz_x_y, g_yy_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_zz_x_x[i] = -g_0_zz_x_x[i] + 2.0 * g_yy_zz_x_x[i] * a_exp;

        g_y_0_0_0_y_zz_x_y[i] = -g_0_zz_x_y[i] + 2.0 * g_yy_zz_x_y[i] * a_exp;

        g_y_0_0_0_y_zz_x_z[i] = -g_0_zz_x_z[i] + 2.0 * g_yy_zz_x_z[i] * a_exp;
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_y_0_0_0_y_zz_y_x, g_y_0_0_0_y_zz_y_y, g_y_0_0_0_y_zz_y_z, g_yy_zz_y_x, g_yy_zz_y_y, g_yy_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_zz_y_x[i] = -g_0_zz_y_x[i] + 2.0 * g_yy_zz_y_x[i] * a_exp;

        g_y_0_0_0_y_zz_y_y[i] = -g_0_zz_y_y[i] + 2.0 * g_yy_zz_y_y[i] * a_exp;

        g_y_0_0_0_y_zz_y_z[i] = -g_0_zz_y_z[i] + 2.0 * g_yy_zz_y_z[i] * a_exp;
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_y_0_0_0_y_zz_z_x, g_y_0_0_0_y_zz_z_y, g_y_0_0_0_y_zz_z_z, g_yy_zz_z_x, g_yy_zz_z_y, g_yy_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_zz_z_x[i] = -g_0_zz_z_x[i] + 2.0 * g_yy_zz_z_x[i] * a_exp;

        g_y_0_0_0_y_zz_z_y[i] = -g_0_zz_z_y[i] + 2.0 * g_yy_zz_z_y[i] * a_exp;

        g_y_0_0_0_y_zz_z_z[i] = -g_0_zz_z_z[i] + 2.0 * g_yy_zz_z_z[i] * a_exp;
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_y_0_0_0_z_xx_x_x, g_y_0_0_0_z_xx_x_y, g_y_0_0_0_z_xx_x_z, g_yz_xx_x_x, g_yz_xx_x_y, g_yz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xx_x_x[i] = 2.0 * g_yz_xx_x_x[i] * a_exp;

        g_y_0_0_0_z_xx_x_y[i] = 2.0 * g_yz_xx_x_y[i] * a_exp;

        g_y_0_0_0_z_xx_x_z[i] = 2.0 * g_yz_xx_x_z[i] * a_exp;
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_y_0_0_0_z_xx_y_x, g_y_0_0_0_z_xx_y_y, g_y_0_0_0_z_xx_y_z, g_yz_xx_y_x, g_yz_xx_y_y, g_yz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xx_y_x[i] = 2.0 * g_yz_xx_y_x[i] * a_exp;

        g_y_0_0_0_z_xx_y_y[i] = 2.0 * g_yz_xx_y_y[i] * a_exp;

        g_y_0_0_0_z_xx_y_z[i] = 2.0 * g_yz_xx_y_z[i] * a_exp;
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_y_0_0_0_z_xx_z_x, g_y_0_0_0_z_xx_z_y, g_y_0_0_0_z_xx_z_z, g_yz_xx_z_x, g_yz_xx_z_y, g_yz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xx_z_x[i] = 2.0 * g_yz_xx_z_x[i] * a_exp;

        g_y_0_0_0_z_xx_z_y[i] = 2.0 * g_yz_xx_z_y[i] * a_exp;

        g_y_0_0_0_z_xx_z_z[i] = 2.0 * g_yz_xx_z_z[i] * a_exp;
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_y_0_0_0_z_xy_x_x, g_y_0_0_0_z_xy_x_y, g_y_0_0_0_z_xy_x_z, g_yz_xy_x_x, g_yz_xy_x_y, g_yz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xy_x_x[i] = 2.0 * g_yz_xy_x_x[i] * a_exp;

        g_y_0_0_0_z_xy_x_y[i] = 2.0 * g_yz_xy_x_y[i] * a_exp;

        g_y_0_0_0_z_xy_x_z[i] = 2.0 * g_yz_xy_x_z[i] * a_exp;
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_y_0_0_0_z_xy_y_x, g_y_0_0_0_z_xy_y_y, g_y_0_0_0_z_xy_y_z, g_yz_xy_y_x, g_yz_xy_y_y, g_yz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xy_y_x[i] = 2.0 * g_yz_xy_y_x[i] * a_exp;

        g_y_0_0_0_z_xy_y_y[i] = 2.0 * g_yz_xy_y_y[i] * a_exp;

        g_y_0_0_0_z_xy_y_z[i] = 2.0 * g_yz_xy_y_z[i] * a_exp;
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_y_0_0_0_z_xy_z_x, g_y_0_0_0_z_xy_z_y, g_y_0_0_0_z_xy_z_z, g_yz_xy_z_x, g_yz_xy_z_y, g_yz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xy_z_x[i] = 2.0 * g_yz_xy_z_x[i] * a_exp;

        g_y_0_0_0_z_xy_z_y[i] = 2.0 * g_yz_xy_z_y[i] * a_exp;

        g_y_0_0_0_z_xy_z_z[i] = 2.0 * g_yz_xy_z_z[i] * a_exp;
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_y_0_0_0_z_xz_x_x, g_y_0_0_0_z_xz_x_y, g_y_0_0_0_z_xz_x_z, g_yz_xz_x_x, g_yz_xz_x_y, g_yz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xz_x_x[i] = 2.0 * g_yz_xz_x_x[i] * a_exp;

        g_y_0_0_0_z_xz_x_y[i] = 2.0 * g_yz_xz_x_y[i] * a_exp;

        g_y_0_0_0_z_xz_x_z[i] = 2.0 * g_yz_xz_x_z[i] * a_exp;
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_y_0_0_0_z_xz_y_x, g_y_0_0_0_z_xz_y_y, g_y_0_0_0_z_xz_y_z, g_yz_xz_y_x, g_yz_xz_y_y, g_yz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xz_y_x[i] = 2.0 * g_yz_xz_y_x[i] * a_exp;

        g_y_0_0_0_z_xz_y_y[i] = 2.0 * g_yz_xz_y_y[i] * a_exp;

        g_y_0_0_0_z_xz_y_z[i] = 2.0 * g_yz_xz_y_z[i] * a_exp;
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_y_0_0_0_z_xz_z_x, g_y_0_0_0_z_xz_z_y, g_y_0_0_0_z_xz_z_z, g_yz_xz_z_x, g_yz_xz_z_y, g_yz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xz_z_x[i] = 2.0 * g_yz_xz_z_x[i] * a_exp;

        g_y_0_0_0_z_xz_z_y[i] = 2.0 * g_yz_xz_z_y[i] * a_exp;

        g_y_0_0_0_z_xz_z_z[i] = 2.0 * g_yz_xz_z_z[i] * a_exp;
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_y_0_0_0_z_yy_x_x, g_y_0_0_0_z_yy_x_y, g_y_0_0_0_z_yy_x_z, g_yz_yy_x_x, g_yz_yy_x_y, g_yz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yy_x_x[i] = 2.0 * g_yz_yy_x_x[i] * a_exp;

        g_y_0_0_0_z_yy_x_y[i] = 2.0 * g_yz_yy_x_y[i] * a_exp;

        g_y_0_0_0_z_yy_x_z[i] = 2.0 * g_yz_yy_x_z[i] * a_exp;
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_y_0_0_0_z_yy_y_x, g_y_0_0_0_z_yy_y_y, g_y_0_0_0_z_yy_y_z, g_yz_yy_y_x, g_yz_yy_y_y, g_yz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yy_y_x[i] = 2.0 * g_yz_yy_y_x[i] * a_exp;

        g_y_0_0_0_z_yy_y_y[i] = 2.0 * g_yz_yy_y_y[i] * a_exp;

        g_y_0_0_0_z_yy_y_z[i] = 2.0 * g_yz_yy_y_z[i] * a_exp;
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_y_0_0_0_z_yy_z_x, g_y_0_0_0_z_yy_z_y, g_y_0_0_0_z_yy_z_z, g_yz_yy_z_x, g_yz_yy_z_y, g_yz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yy_z_x[i] = 2.0 * g_yz_yy_z_x[i] * a_exp;

        g_y_0_0_0_z_yy_z_y[i] = 2.0 * g_yz_yy_z_y[i] * a_exp;

        g_y_0_0_0_z_yy_z_z[i] = 2.0 * g_yz_yy_z_z[i] * a_exp;
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_y_0_0_0_z_yz_x_x, g_y_0_0_0_z_yz_x_y, g_y_0_0_0_z_yz_x_z, g_yz_yz_x_x, g_yz_yz_x_y, g_yz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yz_x_x[i] = 2.0 * g_yz_yz_x_x[i] * a_exp;

        g_y_0_0_0_z_yz_x_y[i] = 2.0 * g_yz_yz_x_y[i] * a_exp;

        g_y_0_0_0_z_yz_x_z[i] = 2.0 * g_yz_yz_x_z[i] * a_exp;
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_y_0_0_0_z_yz_y_x, g_y_0_0_0_z_yz_y_y, g_y_0_0_0_z_yz_y_z, g_yz_yz_y_x, g_yz_yz_y_y, g_yz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yz_y_x[i] = 2.0 * g_yz_yz_y_x[i] * a_exp;

        g_y_0_0_0_z_yz_y_y[i] = 2.0 * g_yz_yz_y_y[i] * a_exp;

        g_y_0_0_0_z_yz_y_z[i] = 2.0 * g_yz_yz_y_z[i] * a_exp;
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_y_0_0_0_z_yz_z_x, g_y_0_0_0_z_yz_z_y, g_y_0_0_0_z_yz_z_z, g_yz_yz_z_x, g_yz_yz_z_y, g_yz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yz_z_x[i] = 2.0 * g_yz_yz_z_x[i] * a_exp;

        g_y_0_0_0_z_yz_z_y[i] = 2.0 * g_yz_yz_z_y[i] * a_exp;

        g_y_0_0_0_z_yz_z_z[i] = 2.0 * g_yz_yz_z_z[i] * a_exp;
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_y_0_0_0_z_zz_x_x, g_y_0_0_0_z_zz_x_y, g_y_0_0_0_z_zz_x_z, g_yz_zz_x_x, g_yz_zz_x_y, g_yz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_zz_x_x[i] = 2.0 * g_yz_zz_x_x[i] * a_exp;

        g_y_0_0_0_z_zz_x_y[i] = 2.0 * g_yz_zz_x_y[i] * a_exp;

        g_y_0_0_0_z_zz_x_z[i] = 2.0 * g_yz_zz_x_z[i] * a_exp;
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_y_0_0_0_z_zz_y_x, g_y_0_0_0_z_zz_y_y, g_y_0_0_0_z_zz_y_z, g_yz_zz_y_x, g_yz_zz_y_y, g_yz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_zz_y_x[i] = 2.0 * g_yz_zz_y_x[i] * a_exp;

        g_y_0_0_0_z_zz_y_y[i] = 2.0 * g_yz_zz_y_y[i] * a_exp;

        g_y_0_0_0_z_zz_y_z[i] = 2.0 * g_yz_zz_y_z[i] * a_exp;
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_y_0_0_0_z_zz_z_x, g_y_0_0_0_z_zz_z_y, g_y_0_0_0_z_zz_z_z, g_yz_zz_z_x, g_yz_zz_z_y, g_yz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_zz_z_x[i] = 2.0 * g_yz_zz_z_x[i] * a_exp;

        g_y_0_0_0_z_zz_z_y[i] = 2.0 * g_yz_zz_z_y[i] * a_exp;

        g_y_0_0_0_z_zz_z_z[i] = 2.0 * g_yz_zz_z_z[i] * a_exp;
    }
    // integrals block (324-327)

    #pragma omp simd aligned(g_xz_xx_x_x, g_xz_xx_x_y, g_xz_xx_x_z, g_z_0_0_0_x_xx_x_x, g_z_0_0_0_x_xx_x_y, g_z_0_0_0_x_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xx_x_x[i] = 2.0 * g_xz_xx_x_x[i] * a_exp;

        g_z_0_0_0_x_xx_x_y[i] = 2.0 * g_xz_xx_x_y[i] * a_exp;

        g_z_0_0_0_x_xx_x_z[i] = 2.0 * g_xz_xx_x_z[i] * a_exp;
    }
    // integrals block (327-330)

    #pragma omp simd aligned(g_xz_xx_y_x, g_xz_xx_y_y, g_xz_xx_y_z, g_z_0_0_0_x_xx_y_x, g_z_0_0_0_x_xx_y_y, g_z_0_0_0_x_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xx_y_x[i] = 2.0 * g_xz_xx_y_x[i] * a_exp;

        g_z_0_0_0_x_xx_y_y[i] = 2.0 * g_xz_xx_y_y[i] * a_exp;

        g_z_0_0_0_x_xx_y_z[i] = 2.0 * g_xz_xx_y_z[i] * a_exp;
    }
    // integrals block (330-333)

    #pragma omp simd aligned(g_xz_xx_z_x, g_xz_xx_z_y, g_xz_xx_z_z, g_z_0_0_0_x_xx_z_x, g_z_0_0_0_x_xx_z_y, g_z_0_0_0_x_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xx_z_x[i] = 2.0 * g_xz_xx_z_x[i] * a_exp;

        g_z_0_0_0_x_xx_z_y[i] = 2.0 * g_xz_xx_z_y[i] * a_exp;

        g_z_0_0_0_x_xx_z_z[i] = 2.0 * g_xz_xx_z_z[i] * a_exp;
    }
    // integrals block (333-336)

    #pragma omp simd aligned(g_xz_xy_x_x, g_xz_xy_x_y, g_xz_xy_x_z, g_z_0_0_0_x_xy_x_x, g_z_0_0_0_x_xy_x_y, g_z_0_0_0_x_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xy_x_x[i] = 2.0 * g_xz_xy_x_x[i] * a_exp;

        g_z_0_0_0_x_xy_x_y[i] = 2.0 * g_xz_xy_x_y[i] * a_exp;

        g_z_0_0_0_x_xy_x_z[i] = 2.0 * g_xz_xy_x_z[i] * a_exp;
    }
    // integrals block (336-339)

    #pragma omp simd aligned(g_xz_xy_y_x, g_xz_xy_y_y, g_xz_xy_y_z, g_z_0_0_0_x_xy_y_x, g_z_0_0_0_x_xy_y_y, g_z_0_0_0_x_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xy_y_x[i] = 2.0 * g_xz_xy_y_x[i] * a_exp;

        g_z_0_0_0_x_xy_y_y[i] = 2.0 * g_xz_xy_y_y[i] * a_exp;

        g_z_0_0_0_x_xy_y_z[i] = 2.0 * g_xz_xy_y_z[i] * a_exp;
    }
    // integrals block (339-342)

    #pragma omp simd aligned(g_xz_xy_z_x, g_xz_xy_z_y, g_xz_xy_z_z, g_z_0_0_0_x_xy_z_x, g_z_0_0_0_x_xy_z_y, g_z_0_0_0_x_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xy_z_x[i] = 2.0 * g_xz_xy_z_x[i] * a_exp;

        g_z_0_0_0_x_xy_z_y[i] = 2.0 * g_xz_xy_z_y[i] * a_exp;

        g_z_0_0_0_x_xy_z_z[i] = 2.0 * g_xz_xy_z_z[i] * a_exp;
    }
    // integrals block (342-345)

    #pragma omp simd aligned(g_xz_xz_x_x, g_xz_xz_x_y, g_xz_xz_x_z, g_z_0_0_0_x_xz_x_x, g_z_0_0_0_x_xz_x_y, g_z_0_0_0_x_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xz_x_x[i] = 2.0 * g_xz_xz_x_x[i] * a_exp;

        g_z_0_0_0_x_xz_x_y[i] = 2.0 * g_xz_xz_x_y[i] * a_exp;

        g_z_0_0_0_x_xz_x_z[i] = 2.0 * g_xz_xz_x_z[i] * a_exp;
    }
    // integrals block (345-348)

    #pragma omp simd aligned(g_xz_xz_y_x, g_xz_xz_y_y, g_xz_xz_y_z, g_z_0_0_0_x_xz_y_x, g_z_0_0_0_x_xz_y_y, g_z_0_0_0_x_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xz_y_x[i] = 2.0 * g_xz_xz_y_x[i] * a_exp;

        g_z_0_0_0_x_xz_y_y[i] = 2.0 * g_xz_xz_y_y[i] * a_exp;

        g_z_0_0_0_x_xz_y_z[i] = 2.0 * g_xz_xz_y_z[i] * a_exp;
    }
    // integrals block (348-351)

    #pragma omp simd aligned(g_xz_xz_z_x, g_xz_xz_z_y, g_xz_xz_z_z, g_z_0_0_0_x_xz_z_x, g_z_0_0_0_x_xz_z_y, g_z_0_0_0_x_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xz_z_x[i] = 2.0 * g_xz_xz_z_x[i] * a_exp;

        g_z_0_0_0_x_xz_z_y[i] = 2.0 * g_xz_xz_z_y[i] * a_exp;

        g_z_0_0_0_x_xz_z_z[i] = 2.0 * g_xz_xz_z_z[i] * a_exp;
    }
    // integrals block (351-354)

    #pragma omp simd aligned(g_xz_yy_x_x, g_xz_yy_x_y, g_xz_yy_x_z, g_z_0_0_0_x_yy_x_x, g_z_0_0_0_x_yy_x_y, g_z_0_0_0_x_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yy_x_x[i] = 2.0 * g_xz_yy_x_x[i] * a_exp;

        g_z_0_0_0_x_yy_x_y[i] = 2.0 * g_xz_yy_x_y[i] * a_exp;

        g_z_0_0_0_x_yy_x_z[i] = 2.0 * g_xz_yy_x_z[i] * a_exp;
    }
    // integrals block (354-357)

    #pragma omp simd aligned(g_xz_yy_y_x, g_xz_yy_y_y, g_xz_yy_y_z, g_z_0_0_0_x_yy_y_x, g_z_0_0_0_x_yy_y_y, g_z_0_0_0_x_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yy_y_x[i] = 2.0 * g_xz_yy_y_x[i] * a_exp;

        g_z_0_0_0_x_yy_y_y[i] = 2.0 * g_xz_yy_y_y[i] * a_exp;

        g_z_0_0_0_x_yy_y_z[i] = 2.0 * g_xz_yy_y_z[i] * a_exp;
    }
    // integrals block (357-360)

    #pragma omp simd aligned(g_xz_yy_z_x, g_xz_yy_z_y, g_xz_yy_z_z, g_z_0_0_0_x_yy_z_x, g_z_0_0_0_x_yy_z_y, g_z_0_0_0_x_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yy_z_x[i] = 2.0 * g_xz_yy_z_x[i] * a_exp;

        g_z_0_0_0_x_yy_z_y[i] = 2.0 * g_xz_yy_z_y[i] * a_exp;

        g_z_0_0_0_x_yy_z_z[i] = 2.0 * g_xz_yy_z_z[i] * a_exp;
    }
    // integrals block (360-363)

    #pragma omp simd aligned(g_xz_yz_x_x, g_xz_yz_x_y, g_xz_yz_x_z, g_z_0_0_0_x_yz_x_x, g_z_0_0_0_x_yz_x_y, g_z_0_0_0_x_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yz_x_x[i] = 2.0 * g_xz_yz_x_x[i] * a_exp;

        g_z_0_0_0_x_yz_x_y[i] = 2.0 * g_xz_yz_x_y[i] * a_exp;

        g_z_0_0_0_x_yz_x_z[i] = 2.0 * g_xz_yz_x_z[i] * a_exp;
    }
    // integrals block (363-366)

    #pragma omp simd aligned(g_xz_yz_y_x, g_xz_yz_y_y, g_xz_yz_y_z, g_z_0_0_0_x_yz_y_x, g_z_0_0_0_x_yz_y_y, g_z_0_0_0_x_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yz_y_x[i] = 2.0 * g_xz_yz_y_x[i] * a_exp;

        g_z_0_0_0_x_yz_y_y[i] = 2.0 * g_xz_yz_y_y[i] * a_exp;

        g_z_0_0_0_x_yz_y_z[i] = 2.0 * g_xz_yz_y_z[i] * a_exp;
    }
    // integrals block (366-369)

    #pragma omp simd aligned(g_xz_yz_z_x, g_xz_yz_z_y, g_xz_yz_z_z, g_z_0_0_0_x_yz_z_x, g_z_0_0_0_x_yz_z_y, g_z_0_0_0_x_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yz_z_x[i] = 2.0 * g_xz_yz_z_x[i] * a_exp;

        g_z_0_0_0_x_yz_z_y[i] = 2.0 * g_xz_yz_z_y[i] * a_exp;

        g_z_0_0_0_x_yz_z_z[i] = 2.0 * g_xz_yz_z_z[i] * a_exp;
    }
    // integrals block (369-372)

    #pragma omp simd aligned(g_xz_zz_x_x, g_xz_zz_x_y, g_xz_zz_x_z, g_z_0_0_0_x_zz_x_x, g_z_0_0_0_x_zz_x_y, g_z_0_0_0_x_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_zz_x_x[i] = 2.0 * g_xz_zz_x_x[i] * a_exp;

        g_z_0_0_0_x_zz_x_y[i] = 2.0 * g_xz_zz_x_y[i] * a_exp;

        g_z_0_0_0_x_zz_x_z[i] = 2.0 * g_xz_zz_x_z[i] * a_exp;
    }
    // integrals block (372-375)

    #pragma omp simd aligned(g_xz_zz_y_x, g_xz_zz_y_y, g_xz_zz_y_z, g_z_0_0_0_x_zz_y_x, g_z_0_0_0_x_zz_y_y, g_z_0_0_0_x_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_zz_y_x[i] = 2.0 * g_xz_zz_y_x[i] * a_exp;

        g_z_0_0_0_x_zz_y_y[i] = 2.0 * g_xz_zz_y_y[i] * a_exp;

        g_z_0_0_0_x_zz_y_z[i] = 2.0 * g_xz_zz_y_z[i] * a_exp;
    }
    // integrals block (375-378)

    #pragma omp simd aligned(g_xz_zz_z_x, g_xz_zz_z_y, g_xz_zz_z_z, g_z_0_0_0_x_zz_z_x, g_z_0_0_0_x_zz_z_y, g_z_0_0_0_x_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_zz_z_x[i] = 2.0 * g_xz_zz_z_x[i] * a_exp;

        g_z_0_0_0_x_zz_z_y[i] = 2.0 * g_xz_zz_z_y[i] * a_exp;

        g_z_0_0_0_x_zz_z_z[i] = 2.0 * g_xz_zz_z_z[i] * a_exp;
    }
    // integrals block (378-381)

    #pragma omp simd aligned(g_yz_xx_x_x, g_yz_xx_x_y, g_yz_xx_x_z, g_z_0_0_0_y_xx_x_x, g_z_0_0_0_y_xx_x_y, g_z_0_0_0_y_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xx_x_x[i] = 2.0 * g_yz_xx_x_x[i] * a_exp;

        g_z_0_0_0_y_xx_x_y[i] = 2.0 * g_yz_xx_x_y[i] * a_exp;

        g_z_0_0_0_y_xx_x_z[i] = 2.0 * g_yz_xx_x_z[i] * a_exp;
    }
    // integrals block (381-384)

    #pragma omp simd aligned(g_yz_xx_y_x, g_yz_xx_y_y, g_yz_xx_y_z, g_z_0_0_0_y_xx_y_x, g_z_0_0_0_y_xx_y_y, g_z_0_0_0_y_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xx_y_x[i] = 2.0 * g_yz_xx_y_x[i] * a_exp;

        g_z_0_0_0_y_xx_y_y[i] = 2.0 * g_yz_xx_y_y[i] * a_exp;

        g_z_0_0_0_y_xx_y_z[i] = 2.0 * g_yz_xx_y_z[i] * a_exp;
    }
    // integrals block (384-387)

    #pragma omp simd aligned(g_yz_xx_z_x, g_yz_xx_z_y, g_yz_xx_z_z, g_z_0_0_0_y_xx_z_x, g_z_0_0_0_y_xx_z_y, g_z_0_0_0_y_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xx_z_x[i] = 2.0 * g_yz_xx_z_x[i] * a_exp;

        g_z_0_0_0_y_xx_z_y[i] = 2.0 * g_yz_xx_z_y[i] * a_exp;

        g_z_0_0_0_y_xx_z_z[i] = 2.0 * g_yz_xx_z_z[i] * a_exp;
    }
    // integrals block (387-390)

    #pragma omp simd aligned(g_yz_xy_x_x, g_yz_xy_x_y, g_yz_xy_x_z, g_z_0_0_0_y_xy_x_x, g_z_0_0_0_y_xy_x_y, g_z_0_0_0_y_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xy_x_x[i] = 2.0 * g_yz_xy_x_x[i] * a_exp;

        g_z_0_0_0_y_xy_x_y[i] = 2.0 * g_yz_xy_x_y[i] * a_exp;

        g_z_0_0_0_y_xy_x_z[i] = 2.0 * g_yz_xy_x_z[i] * a_exp;
    }
    // integrals block (390-393)

    #pragma omp simd aligned(g_yz_xy_y_x, g_yz_xy_y_y, g_yz_xy_y_z, g_z_0_0_0_y_xy_y_x, g_z_0_0_0_y_xy_y_y, g_z_0_0_0_y_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xy_y_x[i] = 2.0 * g_yz_xy_y_x[i] * a_exp;

        g_z_0_0_0_y_xy_y_y[i] = 2.0 * g_yz_xy_y_y[i] * a_exp;

        g_z_0_0_0_y_xy_y_z[i] = 2.0 * g_yz_xy_y_z[i] * a_exp;
    }
    // integrals block (393-396)

    #pragma omp simd aligned(g_yz_xy_z_x, g_yz_xy_z_y, g_yz_xy_z_z, g_z_0_0_0_y_xy_z_x, g_z_0_0_0_y_xy_z_y, g_z_0_0_0_y_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xy_z_x[i] = 2.0 * g_yz_xy_z_x[i] * a_exp;

        g_z_0_0_0_y_xy_z_y[i] = 2.0 * g_yz_xy_z_y[i] * a_exp;

        g_z_0_0_0_y_xy_z_z[i] = 2.0 * g_yz_xy_z_z[i] * a_exp;
    }
    // integrals block (396-399)

    #pragma omp simd aligned(g_yz_xz_x_x, g_yz_xz_x_y, g_yz_xz_x_z, g_z_0_0_0_y_xz_x_x, g_z_0_0_0_y_xz_x_y, g_z_0_0_0_y_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xz_x_x[i] = 2.0 * g_yz_xz_x_x[i] * a_exp;

        g_z_0_0_0_y_xz_x_y[i] = 2.0 * g_yz_xz_x_y[i] * a_exp;

        g_z_0_0_0_y_xz_x_z[i] = 2.0 * g_yz_xz_x_z[i] * a_exp;
    }
    // integrals block (399-402)

    #pragma omp simd aligned(g_yz_xz_y_x, g_yz_xz_y_y, g_yz_xz_y_z, g_z_0_0_0_y_xz_y_x, g_z_0_0_0_y_xz_y_y, g_z_0_0_0_y_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xz_y_x[i] = 2.0 * g_yz_xz_y_x[i] * a_exp;

        g_z_0_0_0_y_xz_y_y[i] = 2.0 * g_yz_xz_y_y[i] * a_exp;

        g_z_0_0_0_y_xz_y_z[i] = 2.0 * g_yz_xz_y_z[i] * a_exp;
    }
    // integrals block (402-405)

    #pragma omp simd aligned(g_yz_xz_z_x, g_yz_xz_z_y, g_yz_xz_z_z, g_z_0_0_0_y_xz_z_x, g_z_0_0_0_y_xz_z_y, g_z_0_0_0_y_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xz_z_x[i] = 2.0 * g_yz_xz_z_x[i] * a_exp;

        g_z_0_0_0_y_xz_z_y[i] = 2.0 * g_yz_xz_z_y[i] * a_exp;

        g_z_0_0_0_y_xz_z_z[i] = 2.0 * g_yz_xz_z_z[i] * a_exp;
    }
    // integrals block (405-408)

    #pragma omp simd aligned(g_yz_yy_x_x, g_yz_yy_x_y, g_yz_yy_x_z, g_z_0_0_0_y_yy_x_x, g_z_0_0_0_y_yy_x_y, g_z_0_0_0_y_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yy_x_x[i] = 2.0 * g_yz_yy_x_x[i] * a_exp;

        g_z_0_0_0_y_yy_x_y[i] = 2.0 * g_yz_yy_x_y[i] * a_exp;

        g_z_0_0_0_y_yy_x_z[i] = 2.0 * g_yz_yy_x_z[i] * a_exp;
    }
    // integrals block (408-411)

    #pragma omp simd aligned(g_yz_yy_y_x, g_yz_yy_y_y, g_yz_yy_y_z, g_z_0_0_0_y_yy_y_x, g_z_0_0_0_y_yy_y_y, g_z_0_0_0_y_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yy_y_x[i] = 2.0 * g_yz_yy_y_x[i] * a_exp;

        g_z_0_0_0_y_yy_y_y[i] = 2.0 * g_yz_yy_y_y[i] * a_exp;

        g_z_0_0_0_y_yy_y_z[i] = 2.0 * g_yz_yy_y_z[i] * a_exp;
    }
    // integrals block (411-414)

    #pragma omp simd aligned(g_yz_yy_z_x, g_yz_yy_z_y, g_yz_yy_z_z, g_z_0_0_0_y_yy_z_x, g_z_0_0_0_y_yy_z_y, g_z_0_0_0_y_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yy_z_x[i] = 2.0 * g_yz_yy_z_x[i] * a_exp;

        g_z_0_0_0_y_yy_z_y[i] = 2.0 * g_yz_yy_z_y[i] * a_exp;

        g_z_0_0_0_y_yy_z_z[i] = 2.0 * g_yz_yy_z_z[i] * a_exp;
    }
    // integrals block (414-417)

    #pragma omp simd aligned(g_yz_yz_x_x, g_yz_yz_x_y, g_yz_yz_x_z, g_z_0_0_0_y_yz_x_x, g_z_0_0_0_y_yz_x_y, g_z_0_0_0_y_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yz_x_x[i] = 2.0 * g_yz_yz_x_x[i] * a_exp;

        g_z_0_0_0_y_yz_x_y[i] = 2.0 * g_yz_yz_x_y[i] * a_exp;

        g_z_0_0_0_y_yz_x_z[i] = 2.0 * g_yz_yz_x_z[i] * a_exp;
    }
    // integrals block (417-420)

    #pragma omp simd aligned(g_yz_yz_y_x, g_yz_yz_y_y, g_yz_yz_y_z, g_z_0_0_0_y_yz_y_x, g_z_0_0_0_y_yz_y_y, g_z_0_0_0_y_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yz_y_x[i] = 2.0 * g_yz_yz_y_x[i] * a_exp;

        g_z_0_0_0_y_yz_y_y[i] = 2.0 * g_yz_yz_y_y[i] * a_exp;

        g_z_0_0_0_y_yz_y_z[i] = 2.0 * g_yz_yz_y_z[i] * a_exp;
    }
    // integrals block (420-423)

    #pragma omp simd aligned(g_yz_yz_z_x, g_yz_yz_z_y, g_yz_yz_z_z, g_z_0_0_0_y_yz_z_x, g_z_0_0_0_y_yz_z_y, g_z_0_0_0_y_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yz_z_x[i] = 2.0 * g_yz_yz_z_x[i] * a_exp;

        g_z_0_0_0_y_yz_z_y[i] = 2.0 * g_yz_yz_z_y[i] * a_exp;

        g_z_0_0_0_y_yz_z_z[i] = 2.0 * g_yz_yz_z_z[i] * a_exp;
    }
    // integrals block (423-426)

    #pragma omp simd aligned(g_yz_zz_x_x, g_yz_zz_x_y, g_yz_zz_x_z, g_z_0_0_0_y_zz_x_x, g_z_0_0_0_y_zz_x_y, g_z_0_0_0_y_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_zz_x_x[i] = 2.0 * g_yz_zz_x_x[i] * a_exp;

        g_z_0_0_0_y_zz_x_y[i] = 2.0 * g_yz_zz_x_y[i] * a_exp;

        g_z_0_0_0_y_zz_x_z[i] = 2.0 * g_yz_zz_x_z[i] * a_exp;
    }
    // integrals block (426-429)

    #pragma omp simd aligned(g_yz_zz_y_x, g_yz_zz_y_y, g_yz_zz_y_z, g_z_0_0_0_y_zz_y_x, g_z_0_0_0_y_zz_y_y, g_z_0_0_0_y_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_zz_y_x[i] = 2.0 * g_yz_zz_y_x[i] * a_exp;

        g_z_0_0_0_y_zz_y_y[i] = 2.0 * g_yz_zz_y_y[i] * a_exp;

        g_z_0_0_0_y_zz_y_z[i] = 2.0 * g_yz_zz_y_z[i] * a_exp;
    }
    // integrals block (429-432)

    #pragma omp simd aligned(g_yz_zz_z_x, g_yz_zz_z_y, g_yz_zz_z_z, g_z_0_0_0_y_zz_z_x, g_z_0_0_0_y_zz_z_y, g_z_0_0_0_y_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_zz_z_x[i] = 2.0 * g_yz_zz_z_x[i] * a_exp;

        g_z_0_0_0_y_zz_z_y[i] = 2.0 * g_yz_zz_z_y[i] * a_exp;

        g_z_0_0_0_y_zz_z_z[i] = 2.0 * g_yz_zz_z_z[i] * a_exp;
    }
    // integrals block (432-435)

    #pragma omp simd aligned(g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_z_0_0_0_z_xx_x_x, g_z_0_0_0_z_xx_x_y, g_z_0_0_0_z_xx_x_z, g_zz_xx_x_x, g_zz_xx_x_y, g_zz_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xx_x_x[i] = -g_0_xx_x_x[i] + 2.0 * g_zz_xx_x_x[i] * a_exp;

        g_z_0_0_0_z_xx_x_y[i] = -g_0_xx_x_y[i] + 2.0 * g_zz_xx_x_y[i] * a_exp;

        g_z_0_0_0_z_xx_x_z[i] = -g_0_xx_x_z[i] + 2.0 * g_zz_xx_x_z[i] * a_exp;
    }
    // integrals block (435-438)

    #pragma omp simd aligned(g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_z_0_0_0_z_xx_y_x, g_z_0_0_0_z_xx_y_y, g_z_0_0_0_z_xx_y_z, g_zz_xx_y_x, g_zz_xx_y_y, g_zz_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xx_y_x[i] = -g_0_xx_y_x[i] + 2.0 * g_zz_xx_y_x[i] * a_exp;

        g_z_0_0_0_z_xx_y_y[i] = -g_0_xx_y_y[i] + 2.0 * g_zz_xx_y_y[i] * a_exp;

        g_z_0_0_0_z_xx_y_z[i] = -g_0_xx_y_z[i] + 2.0 * g_zz_xx_y_z[i] * a_exp;
    }
    // integrals block (438-441)

    #pragma omp simd aligned(g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_z_0_0_0_z_xx_z_x, g_z_0_0_0_z_xx_z_y, g_z_0_0_0_z_xx_z_z, g_zz_xx_z_x, g_zz_xx_z_y, g_zz_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xx_z_x[i] = -g_0_xx_z_x[i] + 2.0 * g_zz_xx_z_x[i] * a_exp;

        g_z_0_0_0_z_xx_z_y[i] = -g_0_xx_z_y[i] + 2.0 * g_zz_xx_z_y[i] * a_exp;

        g_z_0_0_0_z_xx_z_z[i] = -g_0_xx_z_z[i] + 2.0 * g_zz_xx_z_z[i] * a_exp;
    }
    // integrals block (441-444)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_z_0_0_0_z_xy_x_x, g_z_0_0_0_z_xy_x_y, g_z_0_0_0_z_xy_x_z, g_zz_xy_x_x, g_zz_xy_x_y, g_zz_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xy_x_x[i] = -g_0_xy_x_x[i] + 2.0 * g_zz_xy_x_x[i] * a_exp;

        g_z_0_0_0_z_xy_x_y[i] = -g_0_xy_x_y[i] + 2.0 * g_zz_xy_x_y[i] * a_exp;

        g_z_0_0_0_z_xy_x_z[i] = -g_0_xy_x_z[i] + 2.0 * g_zz_xy_x_z[i] * a_exp;
    }
    // integrals block (444-447)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_z_0_0_0_z_xy_y_x, g_z_0_0_0_z_xy_y_y, g_z_0_0_0_z_xy_y_z, g_zz_xy_y_x, g_zz_xy_y_y, g_zz_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xy_y_x[i] = -g_0_xy_y_x[i] + 2.0 * g_zz_xy_y_x[i] * a_exp;

        g_z_0_0_0_z_xy_y_y[i] = -g_0_xy_y_y[i] + 2.0 * g_zz_xy_y_y[i] * a_exp;

        g_z_0_0_0_z_xy_y_z[i] = -g_0_xy_y_z[i] + 2.0 * g_zz_xy_y_z[i] * a_exp;
    }
    // integrals block (447-450)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_z_0_0_0_z_xy_z_x, g_z_0_0_0_z_xy_z_y, g_z_0_0_0_z_xy_z_z, g_zz_xy_z_x, g_zz_xy_z_y, g_zz_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xy_z_x[i] = -g_0_xy_z_x[i] + 2.0 * g_zz_xy_z_x[i] * a_exp;

        g_z_0_0_0_z_xy_z_y[i] = -g_0_xy_z_y[i] + 2.0 * g_zz_xy_z_y[i] * a_exp;

        g_z_0_0_0_z_xy_z_z[i] = -g_0_xy_z_z[i] + 2.0 * g_zz_xy_z_z[i] * a_exp;
    }
    // integrals block (450-453)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_z_0_0_0_z_xz_x_x, g_z_0_0_0_z_xz_x_y, g_z_0_0_0_z_xz_x_z, g_zz_xz_x_x, g_zz_xz_x_y, g_zz_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xz_x_x[i] = -g_0_xz_x_x[i] + 2.0 * g_zz_xz_x_x[i] * a_exp;

        g_z_0_0_0_z_xz_x_y[i] = -g_0_xz_x_y[i] + 2.0 * g_zz_xz_x_y[i] * a_exp;

        g_z_0_0_0_z_xz_x_z[i] = -g_0_xz_x_z[i] + 2.0 * g_zz_xz_x_z[i] * a_exp;
    }
    // integrals block (453-456)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_z_0_0_0_z_xz_y_x, g_z_0_0_0_z_xz_y_y, g_z_0_0_0_z_xz_y_z, g_zz_xz_y_x, g_zz_xz_y_y, g_zz_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xz_y_x[i] = -g_0_xz_y_x[i] + 2.0 * g_zz_xz_y_x[i] * a_exp;

        g_z_0_0_0_z_xz_y_y[i] = -g_0_xz_y_y[i] + 2.0 * g_zz_xz_y_y[i] * a_exp;

        g_z_0_0_0_z_xz_y_z[i] = -g_0_xz_y_z[i] + 2.0 * g_zz_xz_y_z[i] * a_exp;
    }
    // integrals block (456-459)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_z_0_0_0_z_xz_z_x, g_z_0_0_0_z_xz_z_y, g_z_0_0_0_z_xz_z_z, g_zz_xz_z_x, g_zz_xz_z_y, g_zz_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xz_z_x[i] = -g_0_xz_z_x[i] + 2.0 * g_zz_xz_z_x[i] * a_exp;

        g_z_0_0_0_z_xz_z_y[i] = -g_0_xz_z_y[i] + 2.0 * g_zz_xz_z_y[i] * a_exp;

        g_z_0_0_0_z_xz_z_z[i] = -g_0_xz_z_z[i] + 2.0 * g_zz_xz_z_z[i] * a_exp;
    }
    // integrals block (459-462)

    #pragma omp simd aligned(g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_z_0_0_0_z_yy_x_x, g_z_0_0_0_z_yy_x_y, g_z_0_0_0_z_yy_x_z, g_zz_yy_x_x, g_zz_yy_x_y, g_zz_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yy_x_x[i] = -g_0_yy_x_x[i] + 2.0 * g_zz_yy_x_x[i] * a_exp;

        g_z_0_0_0_z_yy_x_y[i] = -g_0_yy_x_y[i] + 2.0 * g_zz_yy_x_y[i] * a_exp;

        g_z_0_0_0_z_yy_x_z[i] = -g_0_yy_x_z[i] + 2.0 * g_zz_yy_x_z[i] * a_exp;
    }
    // integrals block (462-465)

    #pragma omp simd aligned(g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_z_0_0_0_z_yy_y_x, g_z_0_0_0_z_yy_y_y, g_z_0_0_0_z_yy_y_z, g_zz_yy_y_x, g_zz_yy_y_y, g_zz_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yy_y_x[i] = -g_0_yy_y_x[i] + 2.0 * g_zz_yy_y_x[i] * a_exp;

        g_z_0_0_0_z_yy_y_y[i] = -g_0_yy_y_y[i] + 2.0 * g_zz_yy_y_y[i] * a_exp;

        g_z_0_0_0_z_yy_y_z[i] = -g_0_yy_y_z[i] + 2.0 * g_zz_yy_y_z[i] * a_exp;
    }
    // integrals block (465-468)

    #pragma omp simd aligned(g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_z_0_0_0_z_yy_z_x, g_z_0_0_0_z_yy_z_y, g_z_0_0_0_z_yy_z_z, g_zz_yy_z_x, g_zz_yy_z_y, g_zz_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yy_z_x[i] = -g_0_yy_z_x[i] + 2.0 * g_zz_yy_z_x[i] * a_exp;

        g_z_0_0_0_z_yy_z_y[i] = -g_0_yy_z_y[i] + 2.0 * g_zz_yy_z_y[i] * a_exp;

        g_z_0_0_0_z_yy_z_z[i] = -g_0_yy_z_z[i] + 2.0 * g_zz_yy_z_z[i] * a_exp;
    }
    // integrals block (468-471)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_z_0_0_0_z_yz_x_x, g_z_0_0_0_z_yz_x_y, g_z_0_0_0_z_yz_x_z, g_zz_yz_x_x, g_zz_yz_x_y, g_zz_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yz_x_x[i] = -g_0_yz_x_x[i] + 2.0 * g_zz_yz_x_x[i] * a_exp;

        g_z_0_0_0_z_yz_x_y[i] = -g_0_yz_x_y[i] + 2.0 * g_zz_yz_x_y[i] * a_exp;

        g_z_0_0_0_z_yz_x_z[i] = -g_0_yz_x_z[i] + 2.0 * g_zz_yz_x_z[i] * a_exp;
    }
    // integrals block (471-474)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_z_0_0_0_z_yz_y_x, g_z_0_0_0_z_yz_y_y, g_z_0_0_0_z_yz_y_z, g_zz_yz_y_x, g_zz_yz_y_y, g_zz_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yz_y_x[i] = -g_0_yz_y_x[i] + 2.0 * g_zz_yz_y_x[i] * a_exp;

        g_z_0_0_0_z_yz_y_y[i] = -g_0_yz_y_y[i] + 2.0 * g_zz_yz_y_y[i] * a_exp;

        g_z_0_0_0_z_yz_y_z[i] = -g_0_yz_y_z[i] + 2.0 * g_zz_yz_y_z[i] * a_exp;
    }
    // integrals block (474-477)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_z_0_0_0_z_yz_z_x, g_z_0_0_0_z_yz_z_y, g_z_0_0_0_z_yz_z_z, g_zz_yz_z_x, g_zz_yz_z_y, g_zz_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yz_z_x[i] = -g_0_yz_z_x[i] + 2.0 * g_zz_yz_z_x[i] * a_exp;

        g_z_0_0_0_z_yz_z_y[i] = -g_0_yz_z_y[i] + 2.0 * g_zz_yz_z_y[i] * a_exp;

        g_z_0_0_0_z_yz_z_z[i] = -g_0_yz_z_z[i] + 2.0 * g_zz_yz_z_z[i] * a_exp;
    }
    // integrals block (477-480)

    #pragma omp simd aligned(g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_z_0_0_0_z_zz_x_x, g_z_0_0_0_z_zz_x_y, g_z_0_0_0_z_zz_x_z, g_zz_zz_x_x, g_zz_zz_x_y, g_zz_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_zz_x_x[i] = -g_0_zz_x_x[i] + 2.0 * g_zz_zz_x_x[i] * a_exp;

        g_z_0_0_0_z_zz_x_y[i] = -g_0_zz_x_y[i] + 2.0 * g_zz_zz_x_y[i] * a_exp;

        g_z_0_0_0_z_zz_x_z[i] = -g_0_zz_x_z[i] + 2.0 * g_zz_zz_x_z[i] * a_exp;
    }
    // integrals block (480-483)

    #pragma omp simd aligned(g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_z_0_0_0_z_zz_y_x, g_z_0_0_0_z_zz_y_y, g_z_0_0_0_z_zz_y_z, g_zz_zz_y_x, g_zz_zz_y_y, g_zz_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_zz_y_x[i] = -g_0_zz_y_x[i] + 2.0 * g_zz_zz_y_x[i] * a_exp;

        g_z_0_0_0_z_zz_y_y[i] = -g_0_zz_y_y[i] + 2.0 * g_zz_zz_y_y[i] * a_exp;

        g_z_0_0_0_z_zz_y_z[i] = -g_0_zz_y_z[i] + 2.0 * g_zz_zz_y_z[i] * a_exp;
    }
    // integrals block (483-486)

    #pragma omp simd aligned(g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_z_0_0_0_z_zz_z_x, g_z_0_0_0_z_zz_z_y, g_z_0_0_0_z_zz_z_z, g_zz_zz_z_x, g_zz_zz_z_y, g_zz_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_zz_z_x[i] = -g_0_zz_z_x[i] + 2.0 * g_zz_zz_z_x[i] * a_exp;

        g_z_0_0_0_z_zz_z_y[i] = -g_0_zz_z_y[i] + 2.0 * g_zz_zz_z_y[i] * a_exp;

        g_z_0_0_0_z_zz_z_z[i] = -g_0_zz_z_z[i] + 2.0 * g_zz_zz_z_z[i] * a_exp;
    }
}

} // t4c_geom namespace

