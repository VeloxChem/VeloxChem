#include "GeomDeriv1010OfScalarForPDSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_pdsp_0(CSimdArray<double>& buffer_1010_pdsp,
                     const CSimdArray<double>& buffer_sdpp,
                     const CSimdArray<double>& buffer_ddpp,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_pdsp.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1010_pdsp

    auto g_x_0_x_0_x_xx_0_x = buffer_1010_pdsp[0];

    auto g_x_0_x_0_x_xx_0_y = buffer_1010_pdsp[1];

    auto g_x_0_x_0_x_xx_0_z = buffer_1010_pdsp[2];

    auto g_x_0_x_0_x_xy_0_x = buffer_1010_pdsp[3];

    auto g_x_0_x_0_x_xy_0_y = buffer_1010_pdsp[4];

    auto g_x_0_x_0_x_xy_0_z = buffer_1010_pdsp[5];

    auto g_x_0_x_0_x_xz_0_x = buffer_1010_pdsp[6];

    auto g_x_0_x_0_x_xz_0_y = buffer_1010_pdsp[7];

    auto g_x_0_x_0_x_xz_0_z = buffer_1010_pdsp[8];

    auto g_x_0_x_0_x_yy_0_x = buffer_1010_pdsp[9];

    auto g_x_0_x_0_x_yy_0_y = buffer_1010_pdsp[10];

    auto g_x_0_x_0_x_yy_0_z = buffer_1010_pdsp[11];

    auto g_x_0_x_0_x_yz_0_x = buffer_1010_pdsp[12];

    auto g_x_0_x_0_x_yz_0_y = buffer_1010_pdsp[13];

    auto g_x_0_x_0_x_yz_0_z = buffer_1010_pdsp[14];

    auto g_x_0_x_0_x_zz_0_x = buffer_1010_pdsp[15];

    auto g_x_0_x_0_x_zz_0_y = buffer_1010_pdsp[16];

    auto g_x_0_x_0_x_zz_0_z = buffer_1010_pdsp[17];

    auto g_x_0_x_0_y_xx_0_x = buffer_1010_pdsp[18];

    auto g_x_0_x_0_y_xx_0_y = buffer_1010_pdsp[19];

    auto g_x_0_x_0_y_xx_0_z = buffer_1010_pdsp[20];

    auto g_x_0_x_0_y_xy_0_x = buffer_1010_pdsp[21];

    auto g_x_0_x_0_y_xy_0_y = buffer_1010_pdsp[22];

    auto g_x_0_x_0_y_xy_0_z = buffer_1010_pdsp[23];

    auto g_x_0_x_0_y_xz_0_x = buffer_1010_pdsp[24];

    auto g_x_0_x_0_y_xz_0_y = buffer_1010_pdsp[25];

    auto g_x_0_x_0_y_xz_0_z = buffer_1010_pdsp[26];

    auto g_x_0_x_0_y_yy_0_x = buffer_1010_pdsp[27];

    auto g_x_0_x_0_y_yy_0_y = buffer_1010_pdsp[28];

    auto g_x_0_x_0_y_yy_0_z = buffer_1010_pdsp[29];

    auto g_x_0_x_0_y_yz_0_x = buffer_1010_pdsp[30];

    auto g_x_0_x_0_y_yz_0_y = buffer_1010_pdsp[31];

    auto g_x_0_x_0_y_yz_0_z = buffer_1010_pdsp[32];

    auto g_x_0_x_0_y_zz_0_x = buffer_1010_pdsp[33];

    auto g_x_0_x_0_y_zz_0_y = buffer_1010_pdsp[34];

    auto g_x_0_x_0_y_zz_0_z = buffer_1010_pdsp[35];

    auto g_x_0_x_0_z_xx_0_x = buffer_1010_pdsp[36];

    auto g_x_0_x_0_z_xx_0_y = buffer_1010_pdsp[37];

    auto g_x_0_x_0_z_xx_0_z = buffer_1010_pdsp[38];

    auto g_x_0_x_0_z_xy_0_x = buffer_1010_pdsp[39];

    auto g_x_0_x_0_z_xy_0_y = buffer_1010_pdsp[40];

    auto g_x_0_x_0_z_xy_0_z = buffer_1010_pdsp[41];

    auto g_x_0_x_0_z_xz_0_x = buffer_1010_pdsp[42];

    auto g_x_0_x_0_z_xz_0_y = buffer_1010_pdsp[43];

    auto g_x_0_x_0_z_xz_0_z = buffer_1010_pdsp[44];

    auto g_x_0_x_0_z_yy_0_x = buffer_1010_pdsp[45];

    auto g_x_0_x_0_z_yy_0_y = buffer_1010_pdsp[46];

    auto g_x_0_x_0_z_yy_0_z = buffer_1010_pdsp[47];

    auto g_x_0_x_0_z_yz_0_x = buffer_1010_pdsp[48];

    auto g_x_0_x_0_z_yz_0_y = buffer_1010_pdsp[49];

    auto g_x_0_x_0_z_yz_0_z = buffer_1010_pdsp[50];

    auto g_x_0_x_0_z_zz_0_x = buffer_1010_pdsp[51];

    auto g_x_0_x_0_z_zz_0_y = buffer_1010_pdsp[52];

    auto g_x_0_x_0_z_zz_0_z = buffer_1010_pdsp[53];

    auto g_x_0_y_0_x_xx_0_x = buffer_1010_pdsp[54];

    auto g_x_0_y_0_x_xx_0_y = buffer_1010_pdsp[55];

    auto g_x_0_y_0_x_xx_0_z = buffer_1010_pdsp[56];

    auto g_x_0_y_0_x_xy_0_x = buffer_1010_pdsp[57];

    auto g_x_0_y_0_x_xy_0_y = buffer_1010_pdsp[58];

    auto g_x_0_y_0_x_xy_0_z = buffer_1010_pdsp[59];

    auto g_x_0_y_0_x_xz_0_x = buffer_1010_pdsp[60];

    auto g_x_0_y_0_x_xz_0_y = buffer_1010_pdsp[61];

    auto g_x_0_y_0_x_xz_0_z = buffer_1010_pdsp[62];

    auto g_x_0_y_0_x_yy_0_x = buffer_1010_pdsp[63];

    auto g_x_0_y_0_x_yy_0_y = buffer_1010_pdsp[64];

    auto g_x_0_y_0_x_yy_0_z = buffer_1010_pdsp[65];

    auto g_x_0_y_0_x_yz_0_x = buffer_1010_pdsp[66];

    auto g_x_0_y_0_x_yz_0_y = buffer_1010_pdsp[67];

    auto g_x_0_y_0_x_yz_0_z = buffer_1010_pdsp[68];

    auto g_x_0_y_0_x_zz_0_x = buffer_1010_pdsp[69];

    auto g_x_0_y_0_x_zz_0_y = buffer_1010_pdsp[70];

    auto g_x_0_y_0_x_zz_0_z = buffer_1010_pdsp[71];

    auto g_x_0_y_0_y_xx_0_x = buffer_1010_pdsp[72];

    auto g_x_0_y_0_y_xx_0_y = buffer_1010_pdsp[73];

    auto g_x_0_y_0_y_xx_0_z = buffer_1010_pdsp[74];

    auto g_x_0_y_0_y_xy_0_x = buffer_1010_pdsp[75];

    auto g_x_0_y_0_y_xy_0_y = buffer_1010_pdsp[76];

    auto g_x_0_y_0_y_xy_0_z = buffer_1010_pdsp[77];

    auto g_x_0_y_0_y_xz_0_x = buffer_1010_pdsp[78];

    auto g_x_0_y_0_y_xz_0_y = buffer_1010_pdsp[79];

    auto g_x_0_y_0_y_xz_0_z = buffer_1010_pdsp[80];

    auto g_x_0_y_0_y_yy_0_x = buffer_1010_pdsp[81];

    auto g_x_0_y_0_y_yy_0_y = buffer_1010_pdsp[82];

    auto g_x_0_y_0_y_yy_0_z = buffer_1010_pdsp[83];

    auto g_x_0_y_0_y_yz_0_x = buffer_1010_pdsp[84];

    auto g_x_0_y_0_y_yz_0_y = buffer_1010_pdsp[85];

    auto g_x_0_y_0_y_yz_0_z = buffer_1010_pdsp[86];

    auto g_x_0_y_0_y_zz_0_x = buffer_1010_pdsp[87];

    auto g_x_0_y_0_y_zz_0_y = buffer_1010_pdsp[88];

    auto g_x_0_y_0_y_zz_0_z = buffer_1010_pdsp[89];

    auto g_x_0_y_0_z_xx_0_x = buffer_1010_pdsp[90];

    auto g_x_0_y_0_z_xx_0_y = buffer_1010_pdsp[91];

    auto g_x_0_y_0_z_xx_0_z = buffer_1010_pdsp[92];

    auto g_x_0_y_0_z_xy_0_x = buffer_1010_pdsp[93];

    auto g_x_0_y_0_z_xy_0_y = buffer_1010_pdsp[94];

    auto g_x_0_y_0_z_xy_0_z = buffer_1010_pdsp[95];

    auto g_x_0_y_0_z_xz_0_x = buffer_1010_pdsp[96];

    auto g_x_0_y_0_z_xz_0_y = buffer_1010_pdsp[97];

    auto g_x_0_y_0_z_xz_0_z = buffer_1010_pdsp[98];

    auto g_x_0_y_0_z_yy_0_x = buffer_1010_pdsp[99];

    auto g_x_0_y_0_z_yy_0_y = buffer_1010_pdsp[100];

    auto g_x_0_y_0_z_yy_0_z = buffer_1010_pdsp[101];

    auto g_x_0_y_0_z_yz_0_x = buffer_1010_pdsp[102];

    auto g_x_0_y_0_z_yz_0_y = buffer_1010_pdsp[103];

    auto g_x_0_y_0_z_yz_0_z = buffer_1010_pdsp[104];

    auto g_x_0_y_0_z_zz_0_x = buffer_1010_pdsp[105];

    auto g_x_0_y_0_z_zz_0_y = buffer_1010_pdsp[106];

    auto g_x_0_y_0_z_zz_0_z = buffer_1010_pdsp[107];

    auto g_x_0_z_0_x_xx_0_x = buffer_1010_pdsp[108];

    auto g_x_0_z_0_x_xx_0_y = buffer_1010_pdsp[109];

    auto g_x_0_z_0_x_xx_0_z = buffer_1010_pdsp[110];

    auto g_x_0_z_0_x_xy_0_x = buffer_1010_pdsp[111];

    auto g_x_0_z_0_x_xy_0_y = buffer_1010_pdsp[112];

    auto g_x_0_z_0_x_xy_0_z = buffer_1010_pdsp[113];

    auto g_x_0_z_0_x_xz_0_x = buffer_1010_pdsp[114];

    auto g_x_0_z_0_x_xz_0_y = buffer_1010_pdsp[115];

    auto g_x_0_z_0_x_xz_0_z = buffer_1010_pdsp[116];

    auto g_x_0_z_0_x_yy_0_x = buffer_1010_pdsp[117];

    auto g_x_0_z_0_x_yy_0_y = buffer_1010_pdsp[118];

    auto g_x_0_z_0_x_yy_0_z = buffer_1010_pdsp[119];

    auto g_x_0_z_0_x_yz_0_x = buffer_1010_pdsp[120];

    auto g_x_0_z_0_x_yz_0_y = buffer_1010_pdsp[121];

    auto g_x_0_z_0_x_yz_0_z = buffer_1010_pdsp[122];

    auto g_x_0_z_0_x_zz_0_x = buffer_1010_pdsp[123];

    auto g_x_0_z_0_x_zz_0_y = buffer_1010_pdsp[124];

    auto g_x_0_z_0_x_zz_0_z = buffer_1010_pdsp[125];

    auto g_x_0_z_0_y_xx_0_x = buffer_1010_pdsp[126];

    auto g_x_0_z_0_y_xx_0_y = buffer_1010_pdsp[127];

    auto g_x_0_z_0_y_xx_0_z = buffer_1010_pdsp[128];

    auto g_x_0_z_0_y_xy_0_x = buffer_1010_pdsp[129];

    auto g_x_0_z_0_y_xy_0_y = buffer_1010_pdsp[130];

    auto g_x_0_z_0_y_xy_0_z = buffer_1010_pdsp[131];

    auto g_x_0_z_0_y_xz_0_x = buffer_1010_pdsp[132];

    auto g_x_0_z_0_y_xz_0_y = buffer_1010_pdsp[133];

    auto g_x_0_z_0_y_xz_0_z = buffer_1010_pdsp[134];

    auto g_x_0_z_0_y_yy_0_x = buffer_1010_pdsp[135];

    auto g_x_0_z_0_y_yy_0_y = buffer_1010_pdsp[136];

    auto g_x_0_z_0_y_yy_0_z = buffer_1010_pdsp[137];

    auto g_x_0_z_0_y_yz_0_x = buffer_1010_pdsp[138];

    auto g_x_0_z_0_y_yz_0_y = buffer_1010_pdsp[139];

    auto g_x_0_z_0_y_yz_0_z = buffer_1010_pdsp[140];

    auto g_x_0_z_0_y_zz_0_x = buffer_1010_pdsp[141];

    auto g_x_0_z_0_y_zz_0_y = buffer_1010_pdsp[142];

    auto g_x_0_z_0_y_zz_0_z = buffer_1010_pdsp[143];

    auto g_x_0_z_0_z_xx_0_x = buffer_1010_pdsp[144];

    auto g_x_0_z_0_z_xx_0_y = buffer_1010_pdsp[145];

    auto g_x_0_z_0_z_xx_0_z = buffer_1010_pdsp[146];

    auto g_x_0_z_0_z_xy_0_x = buffer_1010_pdsp[147];

    auto g_x_0_z_0_z_xy_0_y = buffer_1010_pdsp[148];

    auto g_x_0_z_0_z_xy_0_z = buffer_1010_pdsp[149];

    auto g_x_0_z_0_z_xz_0_x = buffer_1010_pdsp[150];

    auto g_x_0_z_0_z_xz_0_y = buffer_1010_pdsp[151];

    auto g_x_0_z_0_z_xz_0_z = buffer_1010_pdsp[152];

    auto g_x_0_z_0_z_yy_0_x = buffer_1010_pdsp[153];

    auto g_x_0_z_0_z_yy_0_y = buffer_1010_pdsp[154];

    auto g_x_0_z_0_z_yy_0_z = buffer_1010_pdsp[155];

    auto g_x_0_z_0_z_yz_0_x = buffer_1010_pdsp[156];

    auto g_x_0_z_0_z_yz_0_y = buffer_1010_pdsp[157];

    auto g_x_0_z_0_z_yz_0_z = buffer_1010_pdsp[158];

    auto g_x_0_z_0_z_zz_0_x = buffer_1010_pdsp[159];

    auto g_x_0_z_0_z_zz_0_y = buffer_1010_pdsp[160];

    auto g_x_0_z_0_z_zz_0_z = buffer_1010_pdsp[161];

    auto g_y_0_x_0_x_xx_0_x = buffer_1010_pdsp[162];

    auto g_y_0_x_0_x_xx_0_y = buffer_1010_pdsp[163];

    auto g_y_0_x_0_x_xx_0_z = buffer_1010_pdsp[164];

    auto g_y_0_x_0_x_xy_0_x = buffer_1010_pdsp[165];

    auto g_y_0_x_0_x_xy_0_y = buffer_1010_pdsp[166];

    auto g_y_0_x_0_x_xy_0_z = buffer_1010_pdsp[167];

    auto g_y_0_x_0_x_xz_0_x = buffer_1010_pdsp[168];

    auto g_y_0_x_0_x_xz_0_y = buffer_1010_pdsp[169];

    auto g_y_0_x_0_x_xz_0_z = buffer_1010_pdsp[170];

    auto g_y_0_x_0_x_yy_0_x = buffer_1010_pdsp[171];

    auto g_y_0_x_0_x_yy_0_y = buffer_1010_pdsp[172];

    auto g_y_0_x_0_x_yy_0_z = buffer_1010_pdsp[173];

    auto g_y_0_x_0_x_yz_0_x = buffer_1010_pdsp[174];

    auto g_y_0_x_0_x_yz_0_y = buffer_1010_pdsp[175];

    auto g_y_0_x_0_x_yz_0_z = buffer_1010_pdsp[176];

    auto g_y_0_x_0_x_zz_0_x = buffer_1010_pdsp[177];

    auto g_y_0_x_0_x_zz_0_y = buffer_1010_pdsp[178];

    auto g_y_0_x_0_x_zz_0_z = buffer_1010_pdsp[179];

    auto g_y_0_x_0_y_xx_0_x = buffer_1010_pdsp[180];

    auto g_y_0_x_0_y_xx_0_y = buffer_1010_pdsp[181];

    auto g_y_0_x_0_y_xx_0_z = buffer_1010_pdsp[182];

    auto g_y_0_x_0_y_xy_0_x = buffer_1010_pdsp[183];

    auto g_y_0_x_0_y_xy_0_y = buffer_1010_pdsp[184];

    auto g_y_0_x_0_y_xy_0_z = buffer_1010_pdsp[185];

    auto g_y_0_x_0_y_xz_0_x = buffer_1010_pdsp[186];

    auto g_y_0_x_0_y_xz_0_y = buffer_1010_pdsp[187];

    auto g_y_0_x_0_y_xz_0_z = buffer_1010_pdsp[188];

    auto g_y_0_x_0_y_yy_0_x = buffer_1010_pdsp[189];

    auto g_y_0_x_0_y_yy_0_y = buffer_1010_pdsp[190];

    auto g_y_0_x_0_y_yy_0_z = buffer_1010_pdsp[191];

    auto g_y_0_x_0_y_yz_0_x = buffer_1010_pdsp[192];

    auto g_y_0_x_0_y_yz_0_y = buffer_1010_pdsp[193];

    auto g_y_0_x_0_y_yz_0_z = buffer_1010_pdsp[194];

    auto g_y_0_x_0_y_zz_0_x = buffer_1010_pdsp[195];

    auto g_y_0_x_0_y_zz_0_y = buffer_1010_pdsp[196];

    auto g_y_0_x_0_y_zz_0_z = buffer_1010_pdsp[197];

    auto g_y_0_x_0_z_xx_0_x = buffer_1010_pdsp[198];

    auto g_y_0_x_0_z_xx_0_y = buffer_1010_pdsp[199];

    auto g_y_0_x_0_z_xx_0_z = buffer_1010_pdsp[200];

    auto g_y_0_x_0_z_xy_0_x = buffer_1010_pdsp[201];

    auto g_y_0_x_0_z_xy_0_y = buffer_1010_pdsp[202];

    auto g_y_0_x_0_z_xy_0_z = buffer_1010_pdsp[203];

    auto g_y_0_x_0_z_xz_0_x = buffer_1010_pdsp[204];

    auto g_y_0_x_0_z_xz_0_y = buffer_1010_pdsp[205];

    auto g_y_0_x_0_z_xz_0_z = buffer_1010_pdsp[206];

    auto g_y_0_x_0_z_yy_0_x = buffer_1010_pdsp[207];

    auto g_y_0_x_0_z_yy_0_y = buffer_1010_pdsp[208];

    auto g_y_0_x_0_z_yy_0_z = buffer_1010_pdsp[209];

    auto g_y_0_x_0_z_yz_0_x = buffer_1010_pdsp[210];

    auto g_y_0_x_0_z_yz_0_y = buffer_1010_pdsp[211];

    auto g_y_0_x_0_z_yz_0_z = buffer_1010_pdsp[212];

    auto g_y_0_x_0_z_zz_0_x = buffer_1010_pdsp[213];

    auto g_y_0_x_0_z_zz_0_y = buffer_1010_pdsp[214];

    auto g_y_0_x_0_z_zz_0_z = buffer_1010_pdsp[215];

    auto g_y_0_y_0_x_xx_0_x = buffer_1010_pdsp[216];

    auto g_y_0_y_0_x_xx_0_y = buffer_1010_pdsp[217];

    auto g_y_0_y_0_x_xx_0_z = buffer_1010_pdsp[218];

    auto g_y_0_y_0_x_xy_0_x = buffer_1010_pdsp[219];

    auto g_y_0_y_0_x_xy_0_y = buffer_1010_pdsp[220];

    auto g_y_0_y_0_x_xy_0_z = buffer_1010_pdsp[221];

    auto g_y_0_y_0_x_xz_0_x = buffer_1010_pdsp[222];

    auto g_y_0_y_0_x_xz_0_y = buffer_1010_pdsp[223];

    auto g_y_0_y_0_x_xz_0_z = buffer_1010_pdsp[224];

    auto g_y_0_y_0_x_yy_0_x = buffer_1010_pdsp[225];

    auto g_y_0_y_0_x_yy_0_y = buffer_1010_pdsp[226];

    auto g_y_0_y_0_x_yy_0_z = buffer_1010_pdsp[227];

    auto g_y_0_y_0_x_yz_0_x = buffer_1010_pdsp[228];

    auto g_y_0_y_0_x_yz_0_y = buffer_1010_pdsp[229];

    auto g_y_0_y_0_x_yz_0_z = buffer_1010_pdsp[230];

    auto g_y_0_y_0_x_zz_0_x = buffer_1010_pdsp[231];

    auto g_y_0_y_0_x_zz_0_y = buffer_1010_pdsp[232];

    auto g_y_0_y_0_x_zz_0_z = buffer_1010_pdsp[233];

    auto g_y_0_y_0_y_xx_0_x = buffer_1010_pdsp[234];

    auto g_y_0_y_0_y_xx_0_y = buffer_1010_pdsp[235];

    auto g_y_0_y_0_y_xx_0_z = buffer_1010_pdsp[236];

    auto g_y_0_y_0_y_xy_0_x = buffer_1010_pdsp[237];

    auto g_y_0_y_0_y_xy_0_y = buffer_1010_pdsp[238];

    auto g_y_0_y_0_y_xy_0_z = buffer_1010_pdsp[239];

    auto g_y_0_y_0_y_xz_0_x = buffer_1010_pdsp[240];

    auto g_y_0_y_0_y_xz_0_y = buffer_1010_pdsp[241];

    auto g_y_0_y_0_y_xz_0_z = buffer_1010_pdsp[242];

    auto g_y_0_y_0_y_yy_0_x = buffer_1010_pdsp[243];

    auto g_y_0_y_0_y_yy_0_y = buffer_1010_pdsp[244];

    auto g_y_0_y_0_y_yy_0_z = buffer_1010_pdsp[245];

    auto g_y_0_y_0_y_yz_0_x = buffer_1010_pdsp[246];

    auto g_y_0_y_0_y_yz_0_y = buffer_1010_pdsp[247];

    auto g_y_0_y_0_y_yz_0_z = buffer_1010_pdsp[248];

    auto g_y_0_y_0_y_zz_0_x = buffer_1010_pdsp[249];

    auto g_y_0_y_0_y_zz_0_y = buffer_1010_pdsp[250];

    auto g_y_0_y_0_y_zz_0_z = buffer_1010_pdsp[251];

    auto g_y_0_y_0_z_xx_0_x = buffer_1010_pdsp[252];

    auto g_y_0_y_0_z_xx_0_y = buffer_1010_pdsp[253];

    auto g_y_0_y_0_z_xx_0_z = buffer_1010_pdsp[254];

    auto g_y_0_y_0_z_xy_0_x = buffer_1010_pdsp[255];

    auto g_y_0_y_0_z_xy_0_y = buffer_1010_pdsp[256];

    auto g_y_0_y_0_z_xy_0_z = buffer_1010_pdsp[257];

    auto g_y_0_y_0_z_xz_0_x = buffer_1010_pdsp[258];

    auto g_y_0_y_0_z_xz_0_y = buffer_1010_pdsp[259];

    auto g_y_0_y_0_z_xz_0_z = buffer_1010_pdsp[260];

    auto g_y_0_y_0_z_yy_0_x = buffer_1010_pdsp[261];

    auto g_y_0_y_0_z_yy_0_y = buffer_1010_pdsp[262];

    auto g_y_0_y_0_z_yy_0_z = buffer_1010_pdsp[263];

    auto g_y_0_y_0_z_yz_0_x = buffer_1010_pdsp[264];

    auto g_y_0_y_0_z_yz_0_y = buffer_1010_pdsp[265];

    auto g_y_0_y_0_z_yz_0_z = buffer_1010_pdsp[266];

    auto g_y_0_y_0_z_zz_0_x = buffer_1010_pdsp[267];

    auto g_y_0_y_0_z_zz_0_y = buffer_1010_pdsp[268];

    auto g_y_0_y_0_z_zz_0_z = buffer_1010_pdsp[269];

    auto g_y_0_z_0_x_xx_0_x = buffer_1010_pdsp[270];

    auto g_y_0_z_0_x_xx_0_y = buffer_1010_pdsp[271];

    auto g_y_0_z_0_x_xx_0_z = buffer_1010_pdsp[272];

    auto g_y_0_z_0_x_xy_0_x = buffer_1010_pdsp[273];

    auto g_y_0_z_0_x_xy_0_y = buffer_1010_pdsp[274];

    auto g_y_0_z_0_x_xy_0_z = buffer_1010_pdsp[275];

    auto g_y_0_z_0_x_xz_0_x = buffer_1010_pdsp[276];

    auto g_y_0_z_0_x_xz_0_y = buffer_1010_pdsp[277];

    auto g_y_0_z_0_x_xz_0_z = buffer_1010_pdsp[278];

    auto g_y_0_z_0_x_yy_0_x = buffer_1010_pdsp[279];

    auto g_y_0_z_0_x_yy_0_y = buffer_1010_pdsp[280];

    auto g_y_0_z_0_x_yy_0_z = buffer_1010_pdsp[281];

    auto g_y_0_z_0_x_yz_0_x = buffer_1010_pdsp[282];

    auto g_y_0_z_0_x_yz_0_y = buffer_1010_pdsp[283];

    auto g_y_0_z_0_x_yz_0_z = buffer_1010_pdsp[284];

    auto g_y_0_z_0_x_zz_0_x = buffer_1010_pdsp[285];

    auto g_y_0_z_0_x_zz_0_y = buffer_1010_pdsp[286];

    auto g_y_0_z_0_x_zz_0_z = buffer_1010_pdsp[287];

    auto g_y_0_z_0_y_xx_0_x = buffer_1010_pdsp[288];

    auto g_y_0_z_0_y_xx_0_y = buffer_1010_pdsp[289];

    auto g_y_0_z_0_y_xx_0_z = buffer_1010_pdsp[290];

    auto g_y_0_z_0_y_xy_0_x = buffer_1010_pdsp[291];

    auto g_y_0_z_0_y_xy_0_y = buffer_1010_pdsp[292];

    auto g_y_0_z_0_y_xy_0_z = buffer_1010_pdsp[293];

    auto g_y_0_z_0_y_xz_0_x = buffer_1010_pdsp[294];

    auto g_y_0_z_0_y_xz_0_y = buffer_1010_pdsp[295];

    auto g_y_0_z_0_y_xz_0_z = buffer_1010_pdsp[296];

    auto g_y_0_z_0_y_yy_0_x = buffer_1010_pdsp[297];

    auto g_y_0_z_0_y_yy_0_y = buffer_1010_pdsp[298];

    auto g_y_0_z_0_y_yy_0_z = buffer_1010_pdsp[299];

    auto g_y_0_z_0_y_yz_0_x = buffer_1010_pdsp[300];

    auto g_y_0_z_0_y_yz_0_y = buffer_1010_pdsp[301];

    auto g_y_0_z_0_y_yz_0_z = buffer_1010_pdsp[302];

    auto g_y_0_z_0_y_zz_0_x = buffer_1010_pdsp[303];

    auto g_y_0_z_0_y_zz_0_y = buffer_1010_pdsp[304];

    auto g_y_0_z_0_y_zz_0_z = buffer_1010_pdsp[305];

    auto g_y_0_z_0_z_xx_0_x = buffer_1010_pdsp[306];

    auto g_y_0_z_0_z_xx_0_y = buffer_1010_pdsp[307];

    auto g_y_0_z_0_z_xx_0_z = buffer_1010_pdsp[308];

    auto g_y_0_z_0_z_xy_0_x = buffer_1010_pdsp[309];

    auto g_y_0_z_0_z_xy_0_y = buffer_1010_pdsp[310];

    auto g_y_0_z_0_z_xy_0_z = buffer_1010_pdsp[311];

    auto g_y_0_z_0_z_xz_0_x = buffer_1010_pdsp[312];

    auto g_y_0_z_0_z_xz_0_y = buffer_1010_pdsp[313];

    auto g_y_0_z_0_z_xz_0_z = buffer_1010_pdsp[314];

    auto g_y_0_z_0_z_yy_0_x = buffer_1010_pdsp[315];

    auto g_y_0_z_0_z_yy_0_y = buffer_1010_pdsp[316];

    auto g_y_0_z_0_z_yy_0_z = buffer_1010_pdsp[317];

    auto g_y_0_z_0_z_yz_0_x = buffer_1010_pdsp[318];

    auto g_y_0_z_0_z_yz_0_y = buffer_1010_pdsp[319];

    auto g_y_0_z_0_z_yz_0_z = buffer_1010_pdsp[320];

    auto g_y_0_z_0_z_zz_0_x = buffer_1010_pdsp[321];

    auto g_y_0_z_0_z_zz_0_y = buffer_1010_pdsp[322];

    auto g_y_0_z_0_z_zz_0_z = buffer_1010_pdsp[323];

    auto g_z_0_x_0_x_xx_0_x = buffer_1010_pdsp[324];

    auto g_z_0_x_0_x_xx_0_y = buffer_1010_pdsp[325];

    auto g_z_0_x_0_x_xx_0_z = buffer_1010_pdsp[326];

    auto g_z_0_x_0_x_xy_0_x = buffer_1010_pdsp[327];

    auto g_z_0_x_0_x_xy_0_y = buffer_1010_pdsp[328];

    auto g_z_0_x_0_x_xy_0_z = buffer_1010_pdsp[329];

    auto g_z_0_x_0_x_xz_0_x = buffer_1010_pdsp[330];

    auto g_z_0_x_0_x_xz_0_y = buffer_1010_pdsp[331];

    auto g_z_0_x_0_x_xz_0_z = buffer_1010_pdsp[332];

    auto g_z_0_x_0_x_yy_0_x = buffer_1010_pdsp[333];

    auto g_z_0_x_0_x_yy_0_y = buffer_1010_pdsp[334];

    auto g_z_0_x_0_x_yy_0_z = buffer_1010_pdsp[335];

    auto g_z_0_x_0_x_yz_0_x = buffer_1010_pdsp[336];

    auto g_z_0_x_0_x_yz_0_y = buffer_1010_pdsp[337];

    auto g_z_0_x_0_x_yz_0_z = buffer_1010_pdsp[338];

    auto g_z_0_x_0_x_zz_0_x = buffer_1010_pdsp[339];

    auto g_z_0_x_0_x_zz_0_y = buffer_1010_pdsp[340];

    auto g_z_0_x_0_x_zz_0_z = buffer_1010_pdsp[341];

    auto g_z_0_x_0_y_xx_0_x = buffer_1010_pdsp[342];

    auto g_z_0_x_0_y_xx_0_y = buffer_1010_pdsp[343];

    auto g_z_0_x_0_y_xx_0_z = buffer_1010_pdsp[344];

    auto g_z_0_x_0_y_xy_0_x = buffer_1010_pdsp[345];

    auto g_z_0_x_0_y_xy_0_y = buffer_1010_pdsp[346];

    auto g_z_0_x_0_y_xy_0_z = buffer_1010_pdsp[347];

    auto g_z_0_x_0_y_xz_0_x = buffer_1010_pdsp[348];

    auto g_z_0_x_0_y_xz_0_y = buffer_1010_pdsp[349];

    auto g_z_0_x_0_y_xz_0_z = buffer_1010_pdsp[350];

    auto g_z_0_x_0_y_yy_0_x = buffer_1010_pdsp[351];

    auto g_z_0_x_0_y_yy_0_y = buffer_1010_pdsp[352];

    auto g_z_0_x_0_y_yy_0_z = buffer_1010_pdsp[353];

    auto g_z_0_x_0_y_yz_0_x = buffer_1010_pdsp[354];

    auto g_z_0_x_0_y_yz_0_y = buffer_1010_pdsp[355];

    auto g_z_0_x_0_y_yz_0_z = buffer_1010_pdsp[356];

    auto g_z_0_x_0_y_zz_0_x = buffer_1010_pdsp[357];

    auto g_z_0_x_0_y_zz_0_y = buffer_1010_pdsp[358];

    auto g_z_0_x_0_y_zz_0_z = buffer_1010_pdsp[359];

    auto g_z_0_x_0_z_xx_0_x = buffer_1010_pdsp[360];

    auto g_z_0_x_0_z_xx_0_y = buffer_1010_pdsp[361];

    auto g_z_0_x_0_z_xx_0_z = buffer_1010_pdsp[362];

    auto g_z_0_x_0_z_xy_0_x = buffer_1010_pdsp[363];

    auto g_z_0_x_0_z_xy_0_y = buffer_1010_pdsp[364];

    auto g_z_0_x_0_z_xy_0_z = buffer_1010_pdsp[365];

    auto g_z_0_x_0_z_xz_0_x = buffer_1010_pdsp[366];

    auto g_z_0_x_0_z_xz_0_y = buffer_1010_pdsp[367];

    auto g_z_0_x_0_z_xz_0_z = buffer_1010_pdsp[368];

    auto g_z_0_x_0_z_yy_0_x = buffer_1010_pdsp[369];

    auto g_z_0_x_0_z_yy_0_y = buffer_1010_pdsp[370];

    auto g_z_0_x_0_z_yy_0_z = buffer_1010_pdsp[371];

    auto g_z_0_x_0_z_yz_0_x = buffer_1010_pdsp[372];

    auto g_z_0_x_0_z_yz_0_y = buffer_1010_pdsp[373];

    auto g_z_0_x_0_z_yz_0_z = buffer_1010_pdsp[374];

    auto g_z_0_x_0_z_zz_0_x = buffer_1010_pdsp[375];

    auto g_z_0_x_0_z_zz_0_y = buffer_1010_pdsp[376];

    auto g_z_0_x_0_z_zz_0_z = buffer_1010_pdsp[377];

    auto g_z_0_y_0_x_xx_0_x = buffer_1010_pdsp[378];

    auto g_z_0_y_0_x_xx_0_y = buffer_1010_pdsp[379];

    auto g_z_0_y_0_x_xx_0_z = buffer_1010_pdsp[380];

    auto g_z_0_y_0_x_xy_0_x = buffer_1010_pdsp[381];

    auto g_z_0_y_0_x_xy_0_y = buffer_1010_pdsp[382];

    auto g_z_0_y_0_x_xy_0_z = buffer_1010_pdsp[383];

    auto g_z_0_y_0_x_xz_0_x = buffer_1010_pdsp[384];

    auto g_z_0_y_0_x_xz_0_y = buffer_1010_pdsp[385];

    auto g_z_0_y_0_x_xz_0_z = buffer_1010_pdsp[386];

    auto g_z_0_y_0_x_yy_0_x = buffer_1010_pdsp[387];

    auto g_z_0_y_0_x_yy_0_y = buffer_1010_pdsp[388];

    auto g_z_0_y_0_x_yy_0_z = buffer_1010_pdsp[389];

    auto g_z_0_y_0_x_yz_0_x = buffer_1010_pdsp[390];

    auto g_z_0_y_0_x_yz_0_y = buffer_1010_pdsp[391];

    auto g_z_0_y_0_x_yz_0_z = buffer_1010_pdsp[392];

    auto g_z_0_y_0_x_zz_0_x = buffer_1010_pdsp[393];

    auto g_z_0_y_0_x_zz_0_y = buffer_1010_pdsp[394];

    auto g_z_0_y_0_x_zz_0_z = buffer_1010_pdsp[395];

    auto g_z_0_y_0_y_xx_0_x = buffer_1010_pdsp[396];

    auto g_z_0_y_0_y_xx_0_y = buffer_1010_pdsp[397];

    auto g_z_0_y_0_y_xx_0_z = buffer_1010_pdsp[398];

    auto g_z_0_y_0_y_xy_0_x = buffer_1010_pdsp[399];

    auto g_z_0_y_0_y_xy_0_y = buffer_1010_pdsp[400];

    auto g_z_0_y_0_y_xy_0_z = buffer_1010_pdsp[401];

    auto g_z_0_y_0_y_xz_0_x = buffer_1010_pdsp[402];

    auto g_z_0_y_0_y_xz_0_y = buffer_1010_pdsp[403];

    auto g_z_0_y_0_y_xz_0_z = buffer_1010_pdsp[404];

    auto g_z_0_y_0_y_yy_0_x = buffer_1010_pdsp[405];

    auto g_z_0_y_0_y_yy_0_y = buffer_1010_pdsp[406];

    auto g_z_0_y_0_y_yy_0_z = buffer_1010_pdsp[407];

    auto g_z_0_y_0_y_yz_0_x = buffer_1010_pdsp[408];

    auto g_z_0_y_0_y_yz_0_y = buffer_1010_pdsp[409];

    auto g_z_0_y_0_y_yz_0_z = buffer_1010_pdsp[410];

    auto g_z_0_y_0_y_zz_0_x = buffer_1010_pdsp[411];

    auto g_z_0_y_0_y_zz_0_y = buffer_1010_pdsp[412];

    auto g_z_0_y_0_y_zz_0_z = buffer_1010_pdsp[413];

    auto g_z_0_y_0_z_xx_0_x = buffer_1010_pdsp[414];

    auto g_z_0_y_0_z_xx_0_y = buffer_1010_pdsp[415];

    auto g_z_0_y_0_z_xx_0_z = buffer_1010_pdsp[416];

    auto g_z_0_y_0_z_xy_0_x = buffer_1010_pdsp[417];

    auto g_z_0_y_0_z_xy_0_y = buffer_1010_pdsp[418];

    auto g_z_0_y_0_z_xy_0_z = buffer_1010_pdsp[419];

    auto g_z_0_y_0_z_xz_0_x = buffer_1010_pdsp[420];

    auto g_z_0_y_0_z_xz_0_y = buffer_1010_pdsp[421];

    auto g_z_0_y_0_z_xz_0_z = buffer_1010_pdsp[422];

    auto g_z_0_y_0_z_yy_0_x = buffer_1010_pdsp[423];

    auto g_z_0_y_0_z_yy_0_y = buffer_1010_pdsp[424];

    auto g_z_0_y_0_z_yy_0_z = buffer_1010_pdsp[425];

    auto g_z_0_y_0_z_yz_0_x = buffer_1010_pdsp[426];

    auto g_z_0_y_0_z_yz_0_y = buffer_1010_pdsp[427];

    auto g_z_0_y_0_z_yz_0_z = buffer_1010_pdsp[428];

    auto g_z_0_y_0_z_zz_0_x = buffer_1010_pdsp[429];

    auto g_z_0_y_0_z_zz_0_y = buffer_1010_pdsp[430];

    auto g_z_0_y_0_z_zz_0_z = buffer_1010_pdsp[431];

    auto g_z_0_z_0_x_xx_0_x = buffer_1010_pdsp[432];

    auto g_z_0_z_0_x_xx_0_y = buffer_1010_pdsp[433];

    auto g_z_0_z_0_x_xx_0_z = buffer_1010_pdsp[434];

    auto g_z_0_z_0_x_xy_0_x = buffer_1010_pdsp[435];

    auto g_z_0_z_0_x_xy_0_y = buffer_1010_pdsp[436];

    auto g_z_0_z_0_x_xy_0_z = buffer_1010_pdsp[437];

    auto g_z_0_z_0_x_xz_0_x = buffer_1010_pdsp[438];

    auto g_z_0_z_0_x_xz_0_y = buffer_1010_pdsp[439];

    auto g_z_0_z_0_x_xz_0_z = buffer_1010_pdsp[440];

    auto g_z_0_z_0_x_yy_0_x = buffer_1010_pdsp[441];

    auto g_z_0_z_0_x_yy_0_y = buffer_1010_pdsp[442];

    auto g_z_0_z_0_x_yy_0_z = buffer_1010_pdsp[443];

    auto g_z_0_z_0_x_yz_0_x = buffer_1010_pdsp[444];

    auto g_z_0_z_0_x_yz_0_y = buffer_1010_pdsp[445];

    auto g_z_0_z_0_x_yz_0_z = buffer_1010_pdsp[446];

    auto g_z_0_z_0_x_zz_0_x = buffer_1010_pdsp[447];

    auto g_z_0_z_0_x_zz_0_y = buffer_1010_pdsp[448];

    auto g_z_0_z_0_x_zz_0_z = buffer_1010_pdsp[449];

    auto g_z_0_z_0_y_xx_0_x = buffer_1010_pdsp[450];

    auto g_z_0_z_0_y_xx_0_y = buffer_1010_pdsp[451];

    auto g_z_0_z_0_y_xx_0_z = buffer_1010_pdsp[452];

    auto g_z_0_z_0_y_xy_0_x = buffer_1010_pdsp[453];

    auto g_z_0_z_0_y_xy_0_y = buffer_1010_pdsp[454];

    auto g_z_0_z_0_y_xy_0_z = buffer_1010_pdsp[455];

    auto g_z_0_z_0_y_xz_0_x = buffer_1010_pdsp[456];

    auto g_z_0_z_0_y_xz_0_y = buffer_1010_pdsp[457];

    auto g_z_0_z_0_y_xz_0_z = buffer_1010_pdsp[458];

    auto g_z_0_z_0_y_yy_0_x = buffer_1010_pdsp[459];

    auto g_z_0_z_0_y_yy_0_y = buffer_1010_pdsp[460];

    auto g_z_0_z_0_y_yy_0_z = buffer_1010_pdsp[461];

    auto g_z_0_z_0_y_yz_0_x = buffer_1010_pdsp[462];

    auto g_z_0_z_0_y_yz_0_y = buffer_1010_pdsp[463];

    auto g_z_0_z_0_y_yz_0_z = buffer_1010_pdsp[464];

    auto g_z_0_z_0_y_zz_0_x = buffer_1010_pdsp[465];

    auto g_z_0_z_0_y_zz_0_y = buffer_1010_pdsp[466];

    auto g_z_0_z_0_y_zz_0_z = buffer_1010_pdsp[467];

    auto g_z_0_z_0_z_xx_0_x = buffer_1010_pdsp[468];

    auto g_z_0_z_0_z_xx_0_y = buffer_1010_pdsp[469];

    auto g_z_0_z_0_z_xx_0_z = buffer_1010_pdsp[470];

    auto g_z_0_z_0_z_xy_0_x = buffer_1010_pdsp[471];

    auto g_z_0_z_0_z_xy_0_y = buffer_1010_pdsp[472];

    auto g_z_0_z_0_z_xy_0_z = buffer_1010_pdsp[473];

    auto g_z_0_z_0_z_xz_0_x = buffer_1010_pdsp[474];

    auto g_z_0_z_0_z_xz_0_y = buffer_1010_pdsp[475];

    auto g_z_0_z_0_z_xz_0_z = buffer_1010_pdsp[476];

    auto g_z_0_z_0_z_yy_0_x = buffer_1010_pdsp[477];

    auto g_z_0_z_0_z_yy_0_y = buffer_1010_pdsp[478];

    auto g_z_0_z_0_z_yy_0_z = buffer_1010_pdsp[479];

    auto g_z_0_z_0_z_yz_0_x = buffer_1010_pdsp[480];

    auto g_z_0_z_0_z_yz_0_y = buffer_1010_pdsp[481];

    auto g_z_0_z_0_z_yz_0_z = buffer_1010_pdsp[482];

    auto g_z_0_z_0_z_zz_0_x = buffer_1010_pdsp[483];

    auto g_z_0_z_0_z_zz_0_y = buffer_1010_pdsp[484];

    auto g_z_0_z_0_z_zz_0_z = buffer_1010_pdsp[485];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_x_0_x_0_x_xx_0_x, g_x_0_x_0_x_xx_0_y, g_x_0_x_0_x_xx_0_z, g_xx_xx_x_x, g_xx_xx_x_y, g_xx_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_xx_0_x[i] = -2.0 * g_0_xx_x_x[i] * c_exps[i] + 4.0 * g_xx_xx_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xx_0_y[i] = -2.0 * g_0_xx_x_y[i] * c_exps[i] + 4.0 * g_xx_xx_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xx_0_z[i] = -2.0 * g_0_xx_x_z[i] * c_exps[i] + 4.0 * g_xx_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_x_0_x_0_x_xy_0_x, g_x_0_x_0_x_xy_0_y, g_x_0_x_0_x_xy_0_z, g_xx_xy_x_x, g_xx_xy_x_y, g_xx_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_xy_0_x[i] = -2.0 * g_0_xy_x_x[i] * c_exps[i] + 4.0 * g_xx_xy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xy_0_y[i] = -2.0 * g_0_xy_x_y[i] * c_exps[i] + 4.0 * g_xx_xy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xy_0_z[i] = -2.0 * g_0_xy_x_z[i] * c_exps[i] + 4.0 * g_xx_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_x_0_x_0_x_xz_0_x, g_x_0_x_0_x_xz_0_y, g_x_0_x_0_x_xz_0_z, g_xx_xz_x_x, g_xx_xz_x_y, g_xx_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_xz_0_x[i] = -2.0 * g_0_xz_x_x[i] * c_exps[i] + 4.0 * g_xx_xz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xz_0_y[i] = -2.0 * g_0_xz_x_y[i] * c_exps[i] + 4.0 * g_xx_xz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_xz_0_z[i] = -2.0 * g_0_xz_x_z[i] * c_exps[i] + 4.0 * g_xx_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_x_0_x_0_x_yy_0_x, g_x_0_x_0_x_yy_0_y, g_x_0_x_0_x_yy_0_z, g_xx_yy_x_x, g_xx_yy_x_y, g_xx_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_yy_0_x[i] = -2.0 * g_0_yy_x_x[i] * c_exps[i] + 4.0 * g_xx_yy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yy_0_y[i] = -2.0 * g_0_yy_x_y[i] * c_exps[i] + 4.0 * g_xx_yy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yy_0_z[i] = -2.0 * g_0_yy_x_z[i] * c_exps[i] + 4.0 * g_xx_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_x_0_x_0_x_yz_0_x, g_x_0_x_0_x_yz_0_y, g_x_0_x_0_x_yz_0_z, g_xx_yz_x_x, g_xx_yz_x_y, g_xx_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_yz_0_x[i] = -2.0 * g_0_yz_x_x[i] * c_exps[i] + 4.0 * g_xx_yz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yz_0_y[i] = -2.0 * g_0_yz_x_y[i] * c_exps[i] + 4.0 * g_xx_yz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_yz_0_z[i] = -2.0 * g_0_yz_x_z[i] * c_exps[i] + 4.0 * g_xx_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_x_0_x_0_x_zz_0_x, g_x_0_x_0_x_zz_0_y, g_x_0_x_0_x_zz_0_z, g_xx_zz_x_x, g_xx_zz_x_y, g_xx_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_zz_0_x[i] = -2.0 * g_0_zz_x_x[i] * c_exps[i] + 4.0 * g_xx_zz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_zz_0_y[i] = -2.0 * g_0_zz_x_y[i] * c_exps[i] + 4.0 * g_xx_zz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_zz_0_z[i] = -2.0 * g_0_zz_x_z[i] * c_exps[i] + 4.0 * g_xx_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_0_x_0_y_xx_0_x, g_x_0_x_0_y_xx_0_y, g_x_0_x_0_y_xx_0_z, g_xy_xx_x_x, g_xy_xx_x_y, g_xy_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_xx_0_x[i] = 4.0 * g_xy_xx_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xx_0_y[i] = 4.0 * g_xy_xx_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xx_0_z[i] = 4.0 * g_xy_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_0_x_0_y_xy_0_x, g_x_0_x_0_y_xy_0_y, g_x_0_x_0_y_xy_0_z, g_xy_xy_x_x, g_xy_xy_x_y, g_xy_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_xy_0_x[i] = 4.0 * g_xy_xy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xy_0_y[i] = 4.0 * g_xy_xy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xy_0_z[i] = 4.0 * g_xy_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_0_x_0_y_xz_0_x, g_x_0_x_0_y_xz_0_y, g_x_0_x_0_y_xz_0_z, g_xy_xz_x_x, g_xy_xz_x_y, g_xy_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_xz_0_x[i] = 4.0 * g_xy_xz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xz_0_y[i] = 4.0 * g_xy_xz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_xz_0_z[i] = 4.0 * g_xy_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_0_x_0_y_yy_0_x, g_x_0_x_0_y_yy_0_y, g_x_0_x_0_y_yy_0_z, g_xy_yy_x_x, g_xy_yy_x_y, g_xy_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_yy_0_x[i] = 4.0 * g_xy_yy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yy_0_y[i] = 4.0 * g_xy_yy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yy_0_z[i] = 4.0 * g_xy_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_x_0_x_0_y_yz_0_x, g_x_0_x_0_y_yz_0_y, g_x_0_x_0_y_yz_0_z, g_xy_yz_x_x, g_xy_yz_x_y, g_xy_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_yz_0_x[i] = 4.0 * g_xy_yz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yz_0_y[i] = 4.0 * g_xy_yz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_yz_0_z[i] = 4.0 * g_xy_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_x_0_x_0_y_zz_0_x, g_x_0_x_0_y_zz_0_y, g_x_0_x_0_y_zz_0_z, g_xy_zz_x_x, g_xy_zz_x_y, g_xy_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_zz_0_x[i] = 4.0 * g_xy_zz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_zz_0_y[i] = 4.0 * g_xy_zz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_zz_0_z[i] = 4.0 * g_xy_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_0_x_0_z_xx_0_x, g_x_0_x_0_z_xx_0_y, g_x_0_x_0_z_xx_0_z, g_xz_xx_x_x, g_xz_xx_x_y, g_xz_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_xx_0_x[i] = 4.0 * g_xz_xx_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xx_0_y[i] = 4.0 * g_xz_xx_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xx_0_z[i] = 4.0 * g_xz_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_0_x_0_z_xy_0_x, g_x_0_x_0_z_xy_0_y, g_x_0_x_0_z_xy_0_z, g_xz_xy_x_x, g_xz_xy_x_y, g_xz_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_xy_0_x[i] = 4.0 * g_xz_xy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xy_0_y[i] = 4.0 * g_xz_xy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xy_0_z[i] = 4.0 * g_xz_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_0_x_0_z_xz_0_x, g_x_0_x_0_z_xz_0_y, g_x_0_x_0_z_xz_0_z, g_xz_xz_x_x, g_xz_xz_x_y, g_xz_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_xz_0_x[i] = 4.0 * g_xz_xz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xz_0_y[i] = 4.0 * g_xz_xz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_xz_0_z[i] = 4.0 * g_xz_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_0_x_0_z_yy_0_x, g_x_0_x_0_z_yy_0_y, g_x_0_x_0_z_yy_0_z, g_xz_yy_x_x, g_xz_yy_x_y, g_xz_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_yy_0_x[i] = 4.0 * g_xz_yy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yy_0_y[i] = 4.0 * g_xz_yy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yy_0_z[i] = 4.0 * g_xz_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_0_x_0_z_yz_0_x, g_x_0_x_0_z_yz_0_y, g_x_0_x_0_z_yz_0_z, g_xz_yz_x_x, g_xz_yz_x_y, g_xz_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_yz_0_x[i] = 4.0 * g_xz_yz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yz_0_y[i] = 4.0 * g_xz_yz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_yz_0_z[i] = 4.0 * g_xz_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_0_x_0_z_zz_0_x, g_x_0_x_0_z_zz_0_y, g_x_0_x_0_z_zz_0_z, g_xz_zz_x_x, g_xz_zz_x_y, g_xz_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_zz_0_x[i] = 4.0 * g_xz_zz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_zz_0_y[i] = 4.0 * g_xz_zz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_zz_0_z[i] = 4.0 * g_xz_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_x_0_y_0_x_xx_0_x, g_x_0_y_0_x_xx_0_y, g_x_0_y_0_x_xx_0_z, g_xx_xx_y_x, g_xx_xx_y_y, g_xx_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_xx_0_x[i] = -2.0 * g_0_xx_y_x[i] * c_exps[i] + 4.0 * g_xx_xx_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xx_0_y[i] = -2.0 * g_0_xx_y_y[i] * c_exps[i] + 4.0 * g_xx_xx_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xx_0_z[i] = -2.0 * g_0_xx_y_z[i] * c_exps[i] + 4.0 * g_xx_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_x_0_y_0_x_xy_0_x, g_x_0_y_0_x_xy_0_y, g_x_0_y_0_x_xy_0_z, g_xx_xy_y_x, g_xx_xy_y_y, g_xx_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_xy_0_x[i] = -2.0 * g_0_xy_y_x[i] * c_exps[i] + 4.0 * g_xx_xy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xy_0_y[i] = -2.0 * g_0_xy_y_y[i] * c_exps[i] + 4.0 * g_xx_xy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xy_0_z[i] = -2.0 * g_0_xy_y_z[i] * c_exps[i] + 4.0 * g_xx_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_x_0_y_0_x_xz_0_x, g_x_0_y_0_x_xz_0_y, g_x_0_y_0_x_xz_0_z, g_xx_xz_y_x, g_xx_xz_y_y, g_xx_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_xz_0_x[i] = -2.0 * g_0_xz_y_x[i] * c_exps[i] + 4.0 * g_xx_xz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xz_0_y[i] = -2.0 * g_0_xz_y_y[i] * c_exps[i] + 4.0 * g_xx_xz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_xz_0_z[i] = -2.0 * g_0_xz_y_z[i] * c_exps[i] + 4.0 * g_xx_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_x_0_y_0_x_yy_0_x, g_x_0_y_0_x_yy_0_y, g_x_0_y_0_x_yy_0_z, g_xx_yy_y_x, g_xx_yy_y_y, g_xx_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_yy_0_x[i] = -2.0 * g_0_yy_y_x[i] * c_exps[i] + 4.0 * g_xx_yy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yy_0_y[i] = -2.0 * g_0_yy_y_y[i] * c_exps[i] + 4.0 * g_xx_yy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yy_0_z[i] = -2.0 * g_0_yy_y_z[i] * c_exps[i] + 4.0 * g_xx_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_x_0_y_0_x_yz_0_x, g_x_0_y_0_x_yz_0_y, g_x_0_y_0_x_yz_0_z, g_xx_yz_y_x, g_xx_yz_y_y, g_xx_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_yz_0_x[i] = -2.0 * g_0_yz_y_x[i] * c_exps[i] + 4.0 * g_xx_yz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yz_0_y[i] = -2.0 * g_0_yz_y_y[i] * c_exps[i] + 4.0 * g_xx_yz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_yz_0_z[i] = -2.0 * g_0_yz_y_z[i] * c_exps[i] + 4.0 * g_xx_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_x_0_y_0_x_zz_0_x, g_x_0_y_0_x_zz_0_y, g_x_0_y_0_x_zz_0_z, g_xx_zz_y_x, g_xx_zz_y_y, g_xx_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_zz_0_x[i] = -2.0 * g_0_zz_y_x[i] * c_exps[i] + 4.0 * g_xx_zz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_zz_0_y[i] = -2.0 * g_0_zz_y_y[i] * c_exps[i] + 4.0 * g_xx_zz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_zz_0_z[i] = -2.0 * g_0_zz_y_z[i] * c_exps[i] + 4.0 * g_xx_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_0_y_0_y_xx_0_x, g_x_0_y_0_y_xx_0_y, g_x_0_y_0_y_xx_0_z, g_xy_xx_y_x, g_xy_xx_y_y, g_xy_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_xx_0_x[i] = 4.0 * g_xy_xx_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xx_0_y[i] = 4.0 * g_xy_xx_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xx_0_z[i] = 4.0 * g_xy_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_0_y_0_y_xy_0_x, g_x_0_y_0_y_xy_0_y, g_x_0_y_0_y_xy_0_z, g_xy_xy_y_x, g_xy_xy_y_y, g_xy_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_xy_0_x[i] = 4.0 * g_xy_xy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xy_0_y[i] = 4.0 * g_xy_xy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xy_0_z[i] = 4.0 * g_xy_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_0_y_0_y_xz_0_x, g_x_0_y_0_y_xz_0_y, g_x_0_y_0_y_xz_0_z, g_xy_xz_y_x, g_xy_xz_y_y, g_xy_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_xz_0_x[i] = 4.0 * g_xy_xz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xz_0_y[i] = 4.0 * g_xy_xz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_xz_0_z[i] = 4.0 * g_xy_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_x_0_y_0_y_yy_0_x, g_x_0_y_0_y_yy_0_y, g_x_0_y_0_y_yy_0_z, g_xy_yy_y_x, g_xy_yy_y_y, g_xy_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_yy_0_x[i] = 4.0 * g_xy_yy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yy_0_y[i] = 4.0 * g_xy_yy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yy_0_z[i] = 4.0 * g_xy_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_x_0_y_0_y_yz_0_x, g_x_0_y_0_y_yz_0_y, g_x_0_y_0_y_yz_0_z, g_xy_yz_y_x, g_xy_yz_y_y, g_xy_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_yz_0_x[i] = 4.0 * g_xy_yz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yz_0_y[i] = 4.0 * g_xy_yz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_yz_0_z[i] = 4.0 * g_xy_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_x_0_y_0_y_zz_0_x, g_x_0_y_0_y_zz_0_y, g_x_0_y_0_y_zz_0_z, g_xy_zz_y_x, g_xy_zz_y_y, g_xy_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_zz_0_x[i] = 4.0 * g_xy_zz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_zz_0_y[i] = 4.0 * g_xy_zz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_zz_0_z[i] = 4.0 * g_xy_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_x_0_y_0_z_xx_0_x, g_x_0_y_0_z_xx_0_y, g_x_0_y_0_z_xx_0_z, g_xz_xx_y_x, g_xz_xx_y_y, g_xz_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_xx_0_x[i] = 4.0 * g_xz_xx_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xx_0_y[i] = 4.0 * g_xz_xx_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xx_0_z[i] = 4.0 * g_xz_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_x_0_y_0_z_xy_0_x, g_x_0_y_0_z_xy_0_y, g_x_0_y_0_z_xy_0_z, g_xz_xy_y_x, g_xz_xy_y_y, g_xz_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_xy_0_x[i] = 4.0 * g_xz_xy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xy_0_y[i] = 4.0 * g_xz_xy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xy_0_z[i] = 4.0 * g_xz_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_x_0_y_0_z_xz_0_x, g_x_0_y_0_z_xz_0_y, g_x_0_y_0_z_xz_0_z, g_xz_xz_y_x, g_xz_xz_y_y, g_xz_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_xz_0_x[i] = 4.0 * g_xz_xz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xz_0_y[i] = 4.0 * g_xz_xz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_xz_0_z[i] = 4.0 * g_xz_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_x_0_y_0_z_yy_0_x, g_x_0_y_0_z_yy_0_y, g_x_0_y_0_z_yy_0_z, g_xz_yy_y_x, g_xz_yy_y_y, g_xz_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_yy_0_x[i] = 4.0 * g_xz_yy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yy_0_y[i] = 4.0 * g_xz_yy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yy_0_z[i] = 4.0 * g_xz_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_x_0_y_0_z_yz_0_x, g_x_0_y_0_z_yz_0_y, g_x_0_y_0_z_yz_0_z, g_xz_yz_y_x, g_xz_yz_y_y, g_xz_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_yz_0_x[i] = 4.0 * g_xz_yz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yz_0_y[i] = 4.0 * g_xz_yz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_yz_0_z[i] = 4.0 * g_xz_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_x_0_y_0_z_zz_0_x, g_x_0_y_0_z_zz_0_y, g_x_0_y_0_z_zz_0_z, g_xz_zz_y_x, g_xz_zz_y_y, g_xz_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_zz_0_x[i] = 4.0 * g_xz_zz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_zz_0_y[i] = 4.0 * g_xz_zz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_zz_0_z[i] = 4.0 * g_xz_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_x_0_z_0_x_xx_0_x, g_x_0_z_0_x_xx_0_y, g_x_0_z_0_x_xx_0_z, g_xx_xx_z_x, g_xx_xx_z_y, g_xx_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_xx_0_x[i] = -2.0 * g_0_xx_z_x[i] * c_exps[i] + 4.0 * g_xx_xx_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xx_0_y[i] = -2.0 * g_0_xx_z_y[i] * c_exps[i] + 4.0 * g_xx_xx_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xx_0_z[i] = -2.0 * g_0_xx_z_z[i] * c_exps[i] + 4.0 * g_xx_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_x_0_z_0_x_xy_0_x, g_x_0_z_0_x_xy_0_y, g_x_0_z_0_x_xy_0_z, g_xx_xy_z_x, g_xx_xy_z_y, g_xx_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_xy_0_x[i] = -2.0 * g_0_xy_z_x[i] * c_exps[i] + 4.0 * g_xx_xy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xy_0_y[i] = -2.0 * g_0_xy_z_y[i] * c_exps[i] + 4.0 * g_xx_xy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xy_0_z[i] = -2.0 * g_0_xy_z_z[i] * c_exps[i] + 4.0 * g_xx_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_x_0_z_0_x_xz_0_x, g_x_0_z_0_x_xz_0_y, g_x_0_z_0_x_xz_0_z, g_xx_xz_z_x, g_xx_xz_z_y, g_xx_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_xz_0_x[i] = -2.0 * g_0_xz_z_x[i] * c_exps[i] + 4.0 * g_xx_xz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xz_0_y[i] = -2.0 * g_0_xz_z_y[i] * c_exps[i] + 4.0 * g_xx_xz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_xz_0_z[i] = -2.0 * g_0_xz_z_z[i] * c_exps[i] + 4.0 * g_xx_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_x_0_z_0_x_yy_0_x, g_x_0_z_0_x_yy_0_y, g_x_0_z_0_x_yy_0_z, g_xx_yy_z_x, g_xx_yy_z_y, g_xx_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_yy_0_x[i] = -2.0 * g_0_yy_z_x[i] * c_exps[i] + 4.0 * g_xx_yy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yy_0_y[i] = -2.0 * g_0_yy_z_y[i] * c_exps[i] + 4.0 * g_xx_yy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yy_0_z[i] = -2.0 * g_0_yy_z_z[i] * c_exps[i] + 4.0 * g_xx_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_x_0_z_0_x_yz_0_x, g_x_0_z_0_x_yz_0_y, g_x_0_z_0_x_yz_0_z, g_xx_yz_z_x, g_xx_yz_z_y, g_xx_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_yz_0_x[i] = -2.0 * g_0_yz_z_x[i] * c_exps[i] + 4.0 * g_xx_yz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yz_0_y[i] = -2.0 * g_0_yz_z_y[i] * c_exps[i] + 4.0 * g_xx_yz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_yz_0_z[i] = -2.0 * g_0_yz_z_z[i] * c_exps[i] + 4.0 * g_xx_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_x_0_z_0_x_zz_0_x, g_x_0_z_0_x_zz_0_y, g_x_0_z_0_x_zz_0_z, g_xx_zz_z_x, g_xx_zz_z_y, g_xx_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_zz_0_x[i] = -2.0 * g_0_zz_z_x[i] * c_exps[i] + 4.0 * g_xx_zz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_zz_0_y[i] = -2.0 * g_0_zz_z_y[i] * c_exps[i] + 4.0 * g_xx_zz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_zz_0_z[i] = -2.0 * g_0_zz_z_z[i] * c_exps[i] + 4.0 * g_xx_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_x_0_z_0_y_xx_0_x, g_x_0_z_0_y_xx_0_y, g_x_0_z_0_y_xx_0_z, g_xy_xx_z_x, g_xy_xx_z_y, g_xy_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_xx_0_x[i] = 4.0 * g_xy_xx_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xx_0_y[i] = 4.0 * g_xy_xx_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xx_0_z[i] = 4.0 * g_xy_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_x_0_z_0_y_xy_0_x, g_x_0_z_0_y_xy_0_y, g_x_0_z_0_y_xy_0_z, g_xy_xy_z_x, g_xy_xy_z_y, g_xy_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_xy_0_x[i] = 4.0 * g_xy_xy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xy_0_y[i] = 4.0 * g_xy_xy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xy_0_z[i] = 4.0 * g_xy_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_x_0_z_0_y_xz_0_x, g_x_0_z_0_y_xz_0_y, g_x_0_z_0_y_xz_0_z, g_xy_xz_z_x, g_xy_xz_z_y, g_xy_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_xz_0_x[i] = 4.0 * g_xy_xz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xz_0_y[i] = 4.0 * g_xy_xz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_xz_0_z[i] = 4.0 * g_xy_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_x_0_z_0_y_yy_0_x, g_x_0_z_0_y_yy_0_y, g_x_0_z_0_y_yy_0_z, g_xy_yy_z_x, g_xy_yy_z_y, g_xy_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_yy_0_x[i] = 4.0 * g_xy_yy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yy_0_y[i] = 4.0 * g_xy_yy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yy_0_z[i] = 4.0 * g_xy_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_x_0_z_0_y_yz_0_x, g_x_0_z_0_y_yz_0_y, g_x_0_z_0_y_yz_0_z, g_xy_yz_z_x, g_xy_yz_z_y, g_xy_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_yz_0_x[i] = 4.0 * g_xy_yz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yz_0_y[i] = 4.0 * g_xy_yz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_yz_0_z[i] = 4.0 * g_xy_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_x_0_z_0_y_zz_0_x, g_x_0_z_0_y_zz_0_y, g_x_0_z_0_y_zz_0_z, g_xy_zz_z_x, g_xy_zz_z_y, g_xy_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_zz_0_x[i] = 4.0 * g_xy_zz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_zz_0_y[i] = 4.0 * g_xy_zz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_zz_0_z[i] = 4.0 * g_xy_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_x_0_z_0_z_xx_0_x, g_x_0_z_0_z_xx_0_y, g_x_0_z_0_z_xx_0_z, g_xz_xx_z_x, g_xz_xx_z_y, g_xz_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_xx_0_x[i] = 4.0 * g_xz_xx_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xx_0_y[i] = 4.0 * g_xz_xx_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xx_0_z[i] = 4.0 * g_xz_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_x_0_z_0_z_xy_0_x, g_x_0_z_0_z_xy_0_y, g_x_0_z_0_z_xy_0_z, g_xz_xy_z_x, g_xz_xy_z_y, g_xz_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_xy_0_x[i] = 4.0 * g_xz_xy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xy_0_y[i] = 4.0 * g_xz_xy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xy_0_z[i] = 4.0 * g_xz_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_x_0_z_0_z_xz_0_x, g_x_0_z_0_z_xz_0_y, g_x_0_z_0_z_xz_0_z, g_xz_xz_z_x, g_xz_xz_z_y, g_xz_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_xz_0_x[i] = 4.0 * g_xz_xz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xz_0_y[i] = 4.0 * g_xz_xz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_xz_0_z[i] = 4.0 * g_xz_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_x_0_z_0_z_yy_0_x, g_x_0_z_0_z_yy_0_y, g_x_0_z_0_z_yy_0_z, g_xz_yy_z_x, g_xz_yy_z_y, g_xz_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_yy_0_x[i] = 4.0 * g_xz_yy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yy_0_y[i] = 4.0 * g_xz_yy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yy_0_z[i] = 4.0 * g_xz_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_x_0_z_0_z_yz_0_x, g_x_0_z_0_z_yz_0_y, g_x_0_z_0_z_yz_0_z, g_xz_yz_z_x, g_xz_yz_z_y, g_xz_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_yz_0_x[i] = 4.0 * g_xz_yz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yz_0_y[i] = 4.0 * g_xz_yz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_yz_0_z[i] = 4.0 * g_xz_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_x_0_z_0_z_zz_0_x, g_x_0_z_0_z_zz_0_y, g_x_0_z_0_z_zz_0_z, g_xz_zz_z_x, g_xz_zz_z_y, g_xz_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_zz_0_x[i] = 4.0 * g_xz_zz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_zz_0_y[i] = 4.0 * g_xz_zz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_zz_0_z[i] = 4.0 * g_xz_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_xy_xx_x_x, g_xy_xx_x_y, g_xy_xx_x_z, g_y_0_x_0_x_xx_0_x, g_y_0_x_0_x_xx_0_y, g_y_0_x_0_x_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_xx_0_x[i] = 4.0 * g_xy_xx_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xx_0_y[i] = 4.0 * g_xy_xx_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xx_0_z[i] = 4.0 * g_xy_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_xy_xy_x_x, g_xy_xy_x_y, g_xy_xy_x_z, g_y_0_x_0_x_xy_0_x, g_y_0_x_0_x_xy_0_y, g_y_0_x_0_x_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_xy_0_x[i] = 4.0 * g_xy_xy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xy_0_y[i] = 4.0 * g_xy_xy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xy_0_z[i] = 4.0 * g_xy_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_xy_xz_x_x, g_xy_xz_x_y, g_xy_xz_x_z, g_y_0_x_0_x_xz_0_x, g_y_0_x_0_x_xz_0_y, g_y_0_x_0_x_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_xz_0_x[i] = 4.0 * g_xy_xz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xz_0_y[i] = 4.0 * g_xy_xz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_xz_0_z[i] = 4.0 * g_xy_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_xy_yy_x_x, g_xy_yy_x_y, g_xy_yy_x_z, g_y_0_x_0_x_yy_0_x, g_y_0_x_0_x_yy_0_y, g_y_0_x_0_x_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_yy_0_x[i] = 4.0 * g_xy_yy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yy_0_y[i] = 4.0 * g_xy_yy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yy_0_z[i] = 4.0 * g_xy_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_xy_yz_x_x, g_xy_yz_x_y, g_xy_yz_x_z, g_y_0_x_0_x_yz_0_x, g_y_0_x_0_x_yz_0_y, g_y_0_x_0_x_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_yz_0_x[i] = 4.0 * g_xy_yz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yz_0_y[i] = 4.0 * g_xy_yz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_yz_0_z[i] = 4.0 * g_xy_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_xy_zz_x_x, g_xy_zz_x_y, g_xy_zz_x_z, g_y_0_x_0_x_zz_0_x, g_y_0_x_0_x_zz_0_y, g_y_0_x_0_x_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_zz_0_x[i] = 4.0 * g_xy_zz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_zz_0_y[i] = 4.0 * g_xy_zz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_zz_0_z[i] = 4.0 * g_xy_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_y_0_x_0_y_xx_0_x, g_y_0_x_0_y_xx_0_y, g_y_0_x_0_y_xx_0_z, g_yy_xx_x_x, g_yy_xx_x_y, g_yy_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_xx_0_x[i] = -2.0 * g_0_xx_x_x[i] * c_exps[i] + 4.0 * g_yy_xx_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xx_0_y[i] = -2.0 * g_0_xx_x_y[i] * c_exps[i] + 4.0 * g_yy_xx_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xx_0_z[i] = -2.0 * g_0_xx_x_z[i] * c_exps[i] + 4.0 * g_yy_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_y_0_x_0_y_xy_0_x, g_y_0_x_0_y_xy_0_y, g_y_0_x_0_y_xy_0_z, g_yy_xy_x_x, g_yy_xy_x_y, g_yy_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_xy_0_x[i] = -2.0 * g_0_xy_x_x[i] * c_exps[i] + 4.0 * g_yy_xy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xy_0_y[i] = -2.0 * g_0_xy_x_y[i] * c_exps[i] + 4.0 * g_yy_xy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xy_0_z[i] = -2.0 * g_0_xy_x_z[i] * c_exps[i] + 4.0 * g_yy_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_y_0_x_0_y_xz_0_x, g_y_0_x_0_y_xz_0_y, g_y_0_x_0_y_xz_0_z, g_yy_xz_x_x, g_yy_xz_x_y, g_yy_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_xz_0_x[i] = -2.0 * g_0_xz_x_x[i] * c_exps[i] + 4.0 * g_yy_xz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xz_0_y[i] = -2.0 * g_0_xz_x_y[i] * c_exps[i] + 4.0 * g_yy_xz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_xz_0_z[i] = -2.0 * g_0_xz_x_z[i] * c_exps[i] + 4.0 * g_yy_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_y_0_x_0_y_yy_0_x, g_y_0_x_0_y_yy_0_y, g_y_0_x_0_y_yy_0_z, g_yy_yy_x_x, g_yy_yy_x_y, g_yy_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_yy_0_x[i] = -2.0 * g_0_yy_x_x[i] * c_exps[i] + 4.0 * g_yy_yy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yy_0_y[i] = -2.0 * g_0_yy_x_y[i] * c_exps[i] + 4.0 * g_yy_yy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yy_0_z[i] = -2.0 * g_0_yy_x_z[i] * c_exps[i] + 4.0 * g_yy_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_y_0_x_0_y_yz_0_x, g_y_0_x_0_y_yz_0_y, g_y_0_x_0_y_yz_0_z, g_yy_yz_x_x, g_yy_yz_x_y, g_yy_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_yz_0_x[i] = -2.0 * g_0_yz_x_x[i] * c_exps[i] + 4.0 * g_yy_yz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yz_0_y[i] = -2.0 * g_0_yz_x_y[i] * c_exps[i] + 4.0 * g_yy_yz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_yz_0_z[i] = -2.0 * g_0_yz_x_z[i] * c_exps[i] + 4.0 * g_yy_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_y_0_x_0_y_zz_0_x, g_y_0_x_0_y_zz_0_y, g_y_0_x_0_y_zz_0_z, g_yy_zz_x_x, g_yy_zz_x_y, g_yy_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_zz_0_x[i] = -2.0 * g_0_zz_x_x[i] * c_exps[i] + 4.0 * g_yy_zz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_zz_0_y[i] = -2.0 * g_0_zz_x_y[i] * c_exps[i] + 4.0 * g_yy_zz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_zz_0_z[i] = -2.0 * g_0_zz_x_z[i] * c_exps[i] + 4.0 * g_yy_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_y_0_x_0_z_xx_0_x, g_y_0_x_0_z_xx_0_y, g_y_0_x_0_z_xx_0_z, g_yz_xx_x_x, g_yz_xx_x_y, g_yz_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_xx_0_x[i] = 4.0 * g_yz_xx_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xx_0_y[i] = 4.0 * g_yz_xx_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xx_0_z[i] = 4.0 * g_yz_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_y_0_x_0_z_xy_0_x, g_y_0_x_0_z_xy_0_y, g_y_0_x_0_z_xy_0_z, g_yz_xy_x_x, g_yz_xy_x_y, g_yz_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_xy_0_x[i] = 4.0 * g_yz_xy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xy_0_y[i] = 4.0 * g_yz_xy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xy_0_z[i] = 4.0 * g_yz_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_y_0_x_0_z_xz_0_x, g_y_0_x_0_z_xz_0_y, g_y_0_x_0_z_xz_0_z, g_yz_xz_x_x, g_yz_xz_x_y, g_yz_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_xz_0_x[i] = 4.0 * g_yz_xz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xz_0_y[i] = 4.0 * g_yz_xz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_xz_0_z[i] = 4.0 * g_yz_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_y_0_x_0_z_yy_0_x, g_y_0_x_0_z_yy_0_y, g_y_0_x_0_z_yy_0_z, g_yz_yy_x_x, g_yz_yy_x_y, g_yz_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_yy_0_x[i] = 4.0 * g_yz_yy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yy_0_y[i] = 4.0 * g_yz_yy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yy_0_z[i] = 4.0 * g_yz_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_y_0_x_0_z_yz_0_x, g_y_0_x_0_z_yz_0_y, g_y_0_x_0_z_yz_0_z, g_yz_yz_x_x, g_yz_yz_x_y, g_yz_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_yz_0_x[i] = 4.0 * g_yz_yz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yz_0_y[i] = 4.0 * g_yz_yz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_yz_0_z[i] = 4.0 * g_yz_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_y_0_x_0_z_zz_0_x, g_y_0_x_0_z_zz_0_y, g_y_0_x_0_z_zz_0_z, g_yz_zz_x_x, g_yz_zz_x_y, g_yz_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_zz_0_x[i] = 4.0 * g_yz_zz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_zz_0_y[i] = 4.0 * g_yz_zz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_zz_0_z[i] = 4.0 * g_yz_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_xy_xx_y_x, g_xy_xx_y_y, g_xy_xx_y_z, g_y_0_y_0_x_xx_0_x, g_y_0_y_0_x_xx_0_y, g_y_0_y_0_x_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_xx_0_x[i] = 4.0 * g_xy_xx_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xx_0_y[i] = 4.0 * g_xy_xx_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xx_0_z[i] = 4.0 * g_xy_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_xy_xy_y_x, g_xy_xy_y_y, g_xy_xy_y_z, g_y_0_y_0_x_xy_0_x, g_y_0_y_0_x_xy_0_y, g_y_0_y_0_x_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_xy_0_x[i] = 4.0 * g_xy_xy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xy_0_y[i] = 4.0 * g_xy_xy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xy_0_z[i] = 4.0 * g_xy_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_xy_xz_y_x, g_xy_xz_y_y, g_xy_xz_y_z, g_y_0_y_0_x_xz_0_x, g_y_0_y_0_x_xz_0_y, g_y_0_y_0_x_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_xz_0_x[i] = 4.0 * g_xy_xz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xz_0_y[i] = 4.0 * g_xy_xz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_xz_0_z[i] = 4.0 * g_xy_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_xy_yy_y_x, g_xy_yy_y_y, g_xy_yy_y_z, g_y_0_y_0_x_yy_0_x, g_y_0_y_0_x_yy_0_y, g_y_0_y_0_x_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_yy_0_x[i] = 4.0 * g_xy_yy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yy_0_y[i] = 4.0 * g_xy_yy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yy_0_z[i] = 4.0 * g_xy_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_xy_yz_y_x, g_xy_yz_y_y, g_xy_yz_y_z, g_y_0_y_0_x_yz_0_x, g_y_0_y_0_x_yz_0_y, g_y_0_y_0_x_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_yz_0_x[i] = 4.0 * g_xy_yz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yz_0_y[i] = 4.0 * g_xy_yz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_yz_0_z[i] = 4.0 * g_xy_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_xy_zz_y_x, g_xy_zz_y_y, g_xy_zz_y_z, g_y_0_y_0_x_zz_0_x, g_y_0_y_0_x_zz_0_y, g_y_0_y_0_x_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_zz_0_x[i] = 4.0 * g_xy_zz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_zz_0_y[i] = 4.0 * g_xy_zz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_zz_0_z[i] = 4.0 * g_xy_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_y_0_y_0_y_xx_0_x, g_y_0_y_0_y_xx_0_y, g_y_0_y_0_y_xx_0_z, g_yy_xx_y_x, g_yy_xx_y_y, g_yy_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_xx_0_x[i] = -2.0 * g_0_xx_y_x[i] * c_exps[i] + 4.0 * g_yy_xx_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xx_0_y[i] = -2.0 * g_0_xx_y_y[i] * c_exps[i] + 4.0 * g_yy_xx_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xx_0_z[i] = -2.0 * g_0_xx_y_z[i] * c_exps[i] + 4.0 * g_yy_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_y_0_y_0_y_xy_0_x, g_y_0_y_0_y_xy_0_y, g_y_0_y_0_y_xy_0_z, g_yy_xy_y_x, g_yy_xy_y_y, g_yy_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_xy_0_x[i] = -2.0 * g_0_xy_y_x[i] * c_exps[i] + 4.0 * g_yy_xy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xy_0_y[i] = -2.0 * g_0_xy_y_y[i] * c_exps[i] + 4.0 * g_yy_xy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xy_0_z[i] = -2.0 * g_0_xy_y_z[i] * c_exps[i] + 4.0 * g_yy_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_y_0_y_0_y_xz_0_x, g_y_0_y_0_y_xz_0_y, g_y_0_y_0_y_xz_0_z, g_yy_xz_y_x, g_yy_xz_y_y, g_yy_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_xz_0_x[i] = -2.0 * g_0_xz_y_x[i] * c_exps[i] + 4.0 * g_yy_xz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xz_0_y[i] = -2.0 * g_0_xz_y_y[i] * c_exps[i] + 4.0 * g_yy_xz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_xz_0_z[i] = -2.0 * g_0_xz_y_z[i] * c_exps[i] + 4.0 * g_yy_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_y_0_y_0_y_yy_0_x, g_y_0_y_0_y_yy_0_y, g_y_0_y_0_y_yy_0_z, g_yy_yy_y_x, g_yy_yy_y_y, g_yy_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_yy_0_x[i] = -2.0 * g_0_yy_y_x[i] * c_exps[i] + 4.0 * g_yy_yy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yy_0_y[i] = -2.0 * g_0_yy_y_y[i] * c_exps[i] + 4.0 * g_yy_yy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yy_0_z[i] = -2.0 * g_0_yy_y_z[i] * c_exps[i] + 4.0 * g_yy_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_y_0_y_0_y_yz_0_x, g_y_0_y_0_y_yz_0_y, g_y_0_y_0_y_yz_0_z, g_yy_yz_y_x, g_yy_yz_y_y, g_yy_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_yz_0_x[i] = -2.0 * g_0_yz_y_x[i] * c_exps[i] + 4.0 * g_yy_yz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yz_0_y[i] = -2.0 * g_0_yz_y_y[i] * c_exps[i] + 4.0 * g_yy_yz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_yz_0_z[i] = -2.0 * g_0_yz_y_z[i] * c_exps[i] + 4.0 * g_yy_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_y_0_y_0_y_zz_0_x, g_y_0_y_0_y_zz_0_y, g_y_0_y_0_y_zz_0_z, g_yy_zz_y_x, g_yy_zz_y_y, g_yy_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_zz_0_x[i] = -2.0 * g_0_zz_y_x[i] * c_exps[i] + 4.0 * g_yy_zz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_zz_0_y[i] = -2.0 * g_0_zz_y_y[i] * c_exps[i] + 4.0 * g_yy_zz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_zz_0_z[i] = -2.0 * g_0_zz_y_z[i] * c_exps[i] + 4.0 * g_yy_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_y_0_y_0_z_xx_0_x, g_y_0_y_0_z_xx_0_y, g_y_0_y_0_z_xx_0_z, g_yz_xx_y_x, g_yz_xx_y_y, g_yz_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_xx_0_x[i] = 4.0 * g_yz_xx_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xx_0_y[i] = 4.0 * g_yz_xx_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xx_0_z[i] = 4.0 * g_yz_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_y_0_y_0_z_xy_0_x, g_y_0_y_0_z_xy_0_y, g_y_0_y_0_z_xy_0_z, g_yz_xy_y_x, g_yz_xy_y_y, g_yz_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_xy_0_x[i] = 4.0 * g_yz_xy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xy_0_y[i] = 4.0 * g_yz_xy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xy_0_z[i] = 4.0 * g_yz_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_y_0_y_0_z_xz_0_x, g_y_0_y_0_z_xz_0_y, g_y_0_y_0_z_xz_0_z, g_yz_xz_y_x, g_yz_xz_y_y, g_yz_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_xz_0_x[i] = 4.0 * g_yz_xz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xz_0_y[i] = 4.0 * g_yz_xz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_xz_0_z[i] = 4.0 * g_yz_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_y_0_y_0_z_yy_0_x, g_y_0_y_0_z_yy_0_y, g_y_0_y_0_z_yy_0_z, g_yz_yy_y_x, g_yz_yy_y_y, g_yz_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_yy_0_x[i] = 4.0 * g_yz_yy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yy_0_y[i] = 4.0 * g_yz_yy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yy_0_z[i] = 4.0 * g_yz_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_y_0_y_0_z_yz_0_x, g_y_0_y_0_z_yz_0_y, g_y_0_y_0_z_yz_0_z, g_yz_yz_y_x, g_yz_yz_y_y, g_yz_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_yz_0_x[i] = 4.0 * g_yz_yz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yz_0_y[i] = 4.0 * g_yz_yz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_yz_0_z[i] = 4.0 * g_yz_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_y_0_y_0_z_zz_0_x, g_y_0_y_0_z_zz_0_y, g_y_0_y_0_z_zz_0_z, g_yz_zz_y_x, g_yz_zz_y_y, g_yz_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_zz_0_x[i] = 4.0 * g_yz_zz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_zz_0_y[i] = 4.0 * g_yz_zz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_zz_0_z[i] = 4.0 * g_yz_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_xy_xx_z_x, g_xy_xx_z_y, g_xy_xx_z_z, g_y_0_z_0_x_xx_0_x, g_y_0_z_0_x_xx_0_y, g_y_0_z_0_x_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_xx_0_x[i] = 4.0 * g_xy_xx_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xx_0_y[i] = 4.0 * g_xy_xx_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xx_0_z[i] = 4.0 * g_xy_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_xy_xy_z_x, g_xy_xy_z_y, g_xy_xy_z_z, g_y_0_z_0_x_xy_0_x, g_y_0_z_0_x_xy_0_y, g_y_0_z_0_x_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_xy_0_x[i] = 4.0 * g_xy_xy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xy_0_y[i] = 4.0 * g_xy_xy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xy_0_z[i] = 4.0 * g_xy_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_xy_xz_z_x, g_xy_xz_z_y, g_xy_xz_z_z, g_y_0_z_0_x_xz_0_x, g_y_0_z_0_x_xz_0_y, g_y_0_z_0_x_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_xz_0_x[i] = 4.0 * g_xy_xz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xz_0_y[i] = 4.0 * g_xy_xz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_xz_0_z[i] = 4.0 * g_xy_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_xy_yy_z_x, g_xy_yy_z_y, g_xy_yy_z_z, g_y_0_z_0_x_yy_0_x, g_y_0_z_0_x_yy_0_y, g_y_0_z_0_x_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_yy_0_x[i] = 4.0 * g_xy_yy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yy_0_y[i] = 4.0 * g_xy_yy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yy_0_z[i] = 4.0 * g_xy_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_xy_yz_z_x, g_xy_yz_z_y, g_xy_yz_z_z, g_y_0_z_0_x_yz_0_x, g_y_0_z_0_x_yz_0_y, g_y_0_z_0_x_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_yz_0_x[i] = 4.0 * g_xy_yz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yz_0_y[i] = 4.0 * g_xy_yz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_yz_0_z[i] = 4.0 * g_xy_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_xy_zz_z_x, g_xy_zz_z_y, g_xy_zz_z_z, g_y_0_z_0_x_zz_0_x, g_y_0_z_0_x_zz_0_y, g_y_0_z_0_x_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_zz_0_x[i] = 4.0 * g_xy_zz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_zz_0_y[i] = 4.0 * g_xy_zz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_zz_0_z[i] = 4.0 * g_xy_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_y_0_z_0_y_xx_0_x, g_y_0_z_0_y_xx_0_y, g_y_0_z_0_y_xx_0_z, g_yy_xx_z_x, g_yy_xx_z_y, g_yy_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_xx_0_x[i] = -2.0 * g_0_xx_z_x[i] * c_exps[i] + 4.0 * g_yy_xx_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xx_0_y[i] = -2.0 * g_0_xx_z_y[i] * c_exps[i] + 4.0 * g_yy_xx_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xx_0_z[i] = -2.0 * g_0_xx_z_z[i] * c_exps[i] + 4.0 * g_yy_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_y_0_z_0_y_xy_0_x, g_y_0_z_0_y_xy_0_y, g_y_0_z_0_y_xy_0_z, g_yy_xy_z_x, g_yy_xy_z_y, g_yy_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_xy_0_x[i] = -2.0 * g_0_xy_z_x[i] * c_exps[i] + 4.0 * g_yy_xy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xy_0_y[i] = -2.0 * g_0_xy_z_y[i] * c_exps[i] + 4.0 * g_yy_xy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xy_0_z[i] = -2.0 * g_0_xy_z_z[i] * c_exps[i] + 4.0 * g_yy_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_y_0_z_0_y_xz_0_x, g_y_0_z_0_y_xz_0_y, g_y_0_z_0_y_xz_0_z, g_yy_xz_z_x, g_yy_xz_z_y, g_yy_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_xz_0_x[i] = -2.0 * g_0_xz_z_x[i] * c_exps[i] + 4.0 * g_yy_xz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xz_0_y[i] = -2.0 * g_0_xz_z_y[i] * c_exps[i] + 4.0 * g_yy_xz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_xz_0_z[i] = -2.0 * g_0_xz_z_z[i] * c_exps[i] + 4.0 * g_yy_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_y_0_z_0_y_yy_0_x, g_y_0_z_0_y_yy_0_y, g_y_0_z_0_y_yy_0_z, g_yy_yy_z_x, g_yy_yy_z_y, g_yy_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_yy_0_x[i] = -2.0 * g_0_yy_z_x[i] * c_exps[i] + 4.0 * g_yy_yy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yy_0_y[i] = -2.0 * g_0_yy_z_y[i] * c_exps[i] + 4.0 * g_yy_yy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yy_0_z[i] = -2.0 * g_0_yy_z_z[i] * c_exps[i] + 4.0 * g_yy_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_y_0_z_0_y_yz_0_x, g_y_0_z_0_y_yz_0_y, g_y_0_z_0_y_yz_0_z, g_yy_yz_z_x, g_yy_yz_z_y, g_yy_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_yz_0_x[i] = -2.0 * g_0_yz_z_x[i] * c_exps[i] + 4.0 * g_yy_yz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yz_0_y[i] = -2.0 * g_0_yz_z_y[i] * c_exps[i] + 4.0 * g_yy_yz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_yz_0_z[i] = -2.0 * g_0_yz_z_z[i] * c_exps[i] + 4.0 * g_yy_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_y_0_z_0_y_zz_0_x, g_y_0_z_0_y_zz_0_y, g_y_0_z_0_y_zz_0_z, g_yy_zz_z_x, g_yy_zz_z_y, g_yy_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_zz_0_x[i] = -2.0 * g_0_zz_z_x[i] * c_exps[i] + 4.0 * g_yy_zz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_zz_0_y[i] = -2.0 * g_0_zz_z_y[i] * c_exps[i] + 4.0 * g_yy_zz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_zz_0_z[i] = -2.0 * g_0_zz_z_z[i] * c_exps[i] + 4.0 * g_yy_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_y_0_z_0_z_xx_0_x, g_y_0_z_0_z_xx_0_y, g_y_0_z_0_z_xx_0_z, g_yz_xx_z_x, g_yz_xx_z_y, g_yz_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_xx_0_x[i] = 4.0 * g_yz_xx_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xx_0_y[i] = 4.0 * g_yz_xx_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xx_0_z[i] = 4.0 * g_yz_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_y_0_z_0_z_xy_0_x, g_y_0_z_0_z_xy_0_y, g_y_0_z_0_z_xy_0_z, g_yz_xy_z_x, g_yz_xy_z_y, g_yz_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_xy_0_x[i] = 4.0 * g_yz_xy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xy_0_y[i] = 4.0 * g_yz_xy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xy_0_z[i] = 4.0 * g_yz_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_y_0_z_0_z_xz_0_x, g_y_0_z_0_z_xz_0_y, g_y_0_z_0_z_xz_0_z, g_yz_xz_z_x, g_yz_xz_z_y, g_yz_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_xz_0_x[i] = 4.0 * g_yz_xz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xz_0_y[i] = 4.0 * g_yz_xz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_xz_0_z[i] = 4.0 * g_yz_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_y_0_z_0_z_yy_0_x, g_y_0_z_0_z_yy_0_y, g_y_0_z_0_z_yy_0_z, g_yz_yy_z_x, g_yz_yy_z_y, g_yz_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_yy_0_x[i] = 4.0 * g_yz_yy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yy_0_y[i] = 4.0 * g_yz_yy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yy_0_z[i] = 4.0 * g_yz_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_y_0_z_0_z_yz_0_x, g_y_0_z_0_z_yz_0_y, g_y_0_z_0_z_yz_0_z, g_yz_yz_z_x, g_yz_yz_z_y, g_yz_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_yz_0_x[i] = 4.0 * g_yz_yz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yz_0_y[i] = 4.0 * g_yz_yz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_yz_0_z[i] = 4.0 * g_yz_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_y_0_z_0_z_zz_0_x, g_y_0_z_0_z_zz_0_y, g_y_0_z_0_z_zz_0_z, g_yz_zz_z_x, g_yz_zz_z_y, g_yz_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_zz_0_x[i] = 4.0 * g_yz_zz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_zz_0_y[i] = 4.0 * g_yz_zz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_zz_0_z[i] = 4.0 * g_yz_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (324-327)

    #pragma omp simd aligned(g_xz_xx_x_x, g_xz_xx_x_y, g_xz_xx_x_z, g_z_0_x_0_x_xx_0_x, g_z_0_x_0_x_xx_0_y, g_z_0_x_0_x_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_xx_0_x[i] = 4.0 * g_xz_xx_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xx_0_y[i] = 4.0 * g_xz_xx_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xx_0_z[i] = 4.0 * g_xz_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (327-330)

    #pragma omp simd aligned(g_xz_xy_x_x, g_xz_xy_x_y, g_xz_xy_x_z, g_z_0_x_0_x_xy_0_x, g_z_0_x_0_x_xy_0_y, g_z_0_x_0_x_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_xy_0_x[i] = 4.0 * g_xz_xy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xy_0_y[i] = 4.0 * g_xz_xy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xy_0_z[i] = 4.0 * g_xz_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (330-333)

    #pragma omp simd aligned(g_xz_xz_x_x, g_xz_xz_x_y, g_xz_xz_x_z, g_z_0_x_0_x_xz_0_x, g_z_0_x_0_x_xz_0_y, g_z_0_x_0_x_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_xz_0_x[i] = 4.0 * g_xz_xz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xz_0_y[i] = 4.0 * g_xz_xz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_xz_0_z[i] = 4.0 * g_xz_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (333-336)

    #pragma omp simd aligned(g_xz_yy_x_x, g_xz_yy_x_y, g_xz_yy_x_z, g_z_0_x_0_x_yy_0_x, g_z_0_x_0_x_yy_0_y, g_z_0_x_0_x_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_yy_0_x[i] = 4.0 * g_xz_yy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yy_0_y[i] = 4.0 * g_xz_yy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yy_0_z[i] = 4.0 * g_xz_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (336-339)

    #pragma omp simd aligned(g_xz_yz_x_x, g_xz_yz_x_y, g_xz_yz_x_z, g_z_0_x_0_x_yz_0_x, g_z_0_x_0_x_yz_0_y, g_z_0_x_0_x_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_yz_0_x[i] = 4.0 * g_xz_yz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yz_0_y[i] = 4.0 * g_xz_yz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_yz_0_z[i] = 4.0 * g_xz_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (339-342)

    #pragma omp simd aligned(g_xz_zz_x_x, g_xz_zz_x_y, g_xz_zz_x_z, g_z_0_x_0_x_zz_0_x, g_z_0_x_0_x_zz_0_y, g_z_0_x_0_x_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_zz_0_x[i] = 4.0 * g_xz_zz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_zz_0_y[i] = 4.0 * g_xz_zz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_zz_0_z[i] = 4.0 * g_xz_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (342-345)

    #pragma omp simd aligned(g_yz_xx_x_x, g_yz_xx_x_y, g_yz_xx_x_z, g_z_0_x_0_y_xx_0_x, g_z_0_x_0_y_xx_0_y, g_z_0_x_0_y_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_xx_0_x[i] = 4.0 * g_yz_xx_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xx_0_y[i] = 4.0 * g_yz_xx_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xx_0_z[i] = 4.0 * g_yz_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (345-348)

    #pragma omp simd aligned(g_yz_xy_x_x, g_yz_xy_x_y, g_yz_xy_x_z, g_z_0_x_0_y_xy_0_x, g_z_0_x_0_y_xy_0_y, g_z_0_x_0_y_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_xy_0_x[i] = 4.0 * g_yz_xy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xy_0_y[i] = 4.0 * g_yz_xy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xy_0_z[i] = 4.0 * g_yz_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (348-351)

    #pragma omp simd aligned(g_yz_xz_x_x, g_yz_xz_x_y, g_yz_xz_x_z, g_z_0_x_0_y_xz_0_x, g_z_0_x_0_y_xz_0_y, g_z_0_x_0_y_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_xz_0_x[i] = 4.0 * g_yz_xz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xz_0_y[i] = 4.0 * g_yz_xz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_xz_0_z[i] = 4.0 * g_yz_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (351-354)

    #pragma omp simd aligned(g_yz_yy_x_x, g_yz_yy_x_y, g_yz_yy_x_z, g_z_0_x_0_y_yy_0_x, g_z_0_x_0_y_yy_0_y, g_z_0_x_0_y_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_yy_0_x[i] = 4.0 * g_yz_yy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yy_0_y[i] = 4.0 * g_yz_yy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yy_0_z[i] = 4.0 * g_yz_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (354-357)

    #pragma omp simd aligned(g_yz_yz_x_x, g_yz_yz_x_y, g_yz_yz_x_z, g_z_0_x_0_y_yz_0_x, g_z_0_x_0_y_yz_0_y, g_z_0_x_0_y_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_yz_0_x[i] = 4.0 * g_yz_yz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yz_0_y[i] = 4.0 * g_yz_yz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_yz_0_z[i] = 4.0 * g_yz_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (357-360)

    #pragma omp simd aligned(g_yz_zz_x_x, g_yz_zz_x_y, g_yz_zz_x_z, g_z_0_x_0_y_zz_0_x, g_z_0_x_0_y_zz_0_y, g_z_0_x_0_y_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_zz_0_x[i] = 4.0 * g_yz_zz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_zz_0_y[i] = 4.0 * g_yz_zz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_zz_0_z[i] = 4.0 * g_yz_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (360-363)

    #pragma omp simd aligned(g_0_xx_x_x, g_0_xx_x_y, g_0_xx_x_z, g_z_0_x_0_z_xx_0_x, g_z_0_x_0_z_xx_0_y, g_z_0_x_0_z_xx_0_z, g_zz_xx_x_x, g_zz_xx_x_y, g_zz_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_xx_0_x[i] = -2.0 * g_0_xx_x_x[i] * c_exps[i] + 4.0 * g_zz_xx_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xx_0_y[i] = -2.0 * g_0_xx_x_y[i] * c_exps[i] + 4.0 * g_zz_xx_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xx_0_z[i] = -2.0 * g_0_xx_x_z[i] * c_exps[i] + 4.0 * g_zz_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (363-366)

    #pragma omp simd aligned(g_0_xy_x_x, g_0_xy_x_y, g_0_xy_x_z, g_z_0_x_0_z_xy_0_x, g_z_0_x_0_z_xy_0_y, g_z_0_x_0_z_xy_0_z, g_zz_xy_x_x, g_zz_xy_x_y, g_zz_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_xy_0_x[i] = -2.0 * g_0_xy_x_x[i] * c_exps[i] + 4.0 * g_zz_xy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xy_0_y[i] = -2.0 * g_0_xy_x_y[i] * c_exps[i] + 4.0 * g_zz_xy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xy_0_z[i] = -2.0 * g_0_xy_x_z[i] * c_exps[i] + 4.0 * g_zz_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (366-369)

    #pragma omp simd aligned(g_0_xz_x_x, g_0_xz_x_y, g_0_xz_x_z, g_z_0_x_0_z_xz_0_x, g_z_0_x_0_z_xz_0_y, g_z_0_x_0_z_xz_0_z, g_zz_xz_x_x, g_zz_xz_x_y, g_zz_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_xz_0_x[i] = -2.0 * g_0_xz_x_x[i] * c_exps[i] + 4.0 * g_zz_xz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xz_0_y[i] = -2.0 * g_0_xz_x_y[i] * c_exps[i] + 4.0 * g_zz_xz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_xz_0_z[i] = -2.0 * g_0_xz_x_z[i] * c_exps[i] + 4.0 * g_zz_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (369-372)

    #pragma omp simd aligned(g_0_yy_x_x, g_0_yy_x_y, g_0_yy_x_z, g_z_0_x_0_z_yy_0_x, g_z_0_x_0_z_yy_0_y, g_z_0_x_0_z_yy_0_z, g_zz_yy_x_x, g_zz_yy_x_y, g_zz_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_yy_0_x[i] = -2.0 * g_0_yy_x_x[i] * c_exps[i] + 4.0 * g_zz_yy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yy_0_y[i] = -2.0 * g_0_yy_x_y[i] * c_exps[i] + 4.0 * g_zz_yy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yy_0_z[i] = -2.0 * g_0_yy_x_z[i] * c_exps[i] + 4.0 * g_zz_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (372-375)

    #pragma omp simd aligned(g_0_yz_x_x, g_0_yz_x_y, g_0_yz_x_z, g_z_0_x_0_z_yz_0_x, g_z_0_x_0_z_yz_0_y, g_z_0_x_0_z_yz_0_z, g_zz_yz_x_x, g_zz_yz_x_y, g_zz_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_yz_0_x[i] = -2.0 * g_0_yz_x_x[i] * c_exps[i] + 4.0 * g_zz_yz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yz_0_y[i] = -2.0 * g_0_yz_x_y[i] * c_exps[i] + 4.0 * g_zz_yz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_yz_0_z[i] = -2.0 * g_0_yz_x_z[i] * c_exps[i] + 4.0 * g_zz_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (375-378)

    #pragma omp simd aligned(g_0_zz_x_x, g_0_zz_x_y, g_0_zz_x_z, g_z_0_x_0_z_zz_0_x, g_z_0_x_0_z_zz_0_y, g_z_0_x_0_z_zz_0_z, g_zz_zz_x_x, g_zz_zz_x_y, g_zz_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_zz_0_x[i] = -2.0 * g_0_zz_x_x[i] * c_exps[i] + 4.0 * g_zz_zz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_zz_0_y[i] = -2.0 * g_0_zz_x_y[i] * c_exps[i] + 4.0 * g_zz_zz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_zz_0_z[i] = -2.0 * g_0_zz_x_z[i] * c_exps[i] + 4.0 * g_zz_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (378-381)

    #pragma omp simd aligned(g_xz_xx_y_x, g_xz_xx_y_y, g_xz_xx_y_z, g_z_0_y_0_x_xx_0_x, g_z_0_y_0_x_xx_0_y, g_z_0_y_0_x_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_xx_0_x[i] = 4.0 * g_xz_xx_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xx_0_y[i] = 4.0 * g_xz_xx_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xx_0_z[i] = 4.0 * g_xz_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (381-384)

    #pragma omp simd aligned(g_xz_xy_y_x, g_xz_xy_y_y, g_xz_xy_y_z, g_z_0_y_0_x_xy_0_x, g_z_0_y_0_x_xy_0_y, g_z_0_y_0_x_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_xy_0_x[i] = 4.0 * g_xz_xy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xy_0_y[i] = 4.0 * g_xz_xy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xy_0_z[i] = 4.0 * g_xz_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (384-387)

    #pragma omp simd aligned(g_xz_xz_y_x, g_xz_xz_y_y, g_xz_xz_y_z, g_z_0_y_0_x_xz_0_x, g_z_0_y_0_x_xz_0_y, g_z_0_y_0_x_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_xz_0_x[i] = 4.0 * g_xz_xz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xz_0_y[i] = 4.0 * g_xz_xz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_xz_0_z[i] = 4.0 * g_xz_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (387-390)

    #pragma omp simd aligned(g_xz_yy_y_x, g_xz_yy_y_y, g_xz_yy_y_z, g_z_0_y_0_x_yy_0_x, g_z_0_y_0_x_yy_0_y, g_z_0_y_0_x_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_yy_0_x[i] = 4.0 * g_xz_yy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yy_0_y[i] = 4.0 * g_xz_yy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yy_0_z[i] = 4.0 * g_xz_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (390-393)

    #pragma omp simd aligned(g_xz_yz_y_x, g_xz_yz_y_y, g_xz_yz_y_z, g_z_0_y_0_x_yz_0_x, g_z_0_y_0_x_yz_0_y, g_z_0_y_0_x_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_yz_0_x[i] = 4.0 * g_xz_yz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yz_0_y[i] = 4.0 * g_xz_yz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_yz_0_z[i] = 4.0 * g_xz_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (393-396)

    #pragma omp simd aligned(g_xz_zz_y_x, g_xz_zz_y_y, g_xz_zz_y_z, g_z_0_y_0_x_zz_0_x, g_z_0_y_0_x_zz_0_y, g_z_0_y_0_x_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_zz_0_x[i] = 4.0 * g_xz_zz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_zz_0_y[i] = 4.0 * g_xz_zz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_zz_0_z[i] = 4.0 * g_xz_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (396-399)

    #pragma omp simd aligned(g_yz_xx_y_x, g_yz_xx_y_y, g_yz_xx_y_z, g_z_0_y_0_y_xx_0_x, g_z_0_y_0_y_xx_0_y, g_z_0_y_0_y_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_xx_0_x[i] = 4.0 * g_yz_xx_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xx_0_y[i] = 4.0 * g_yz_xx_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xx_0_z[i] = 4.0 * g_yz_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (399-402)

    #pragma omp simd aligned(g_yz_xy_y_x, g_yz_xy_y_y, g_yz_xy_y_z, g_z_0_y_0_y_xy_0_x, g_z_0_y_0_y_xy_0_y, g_z_0_y_0_y_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_xy_0_x[i] = 4.0 * g_yz_xy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xy_0_y[i] = 4.0 * g_yz_xy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xy_0_z[i] = 4.0 * g_yz_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (402-405)

    #pragma omp simd aligned(g_yz_xz_y_x, g_yz_xz_y_y, g_yz_xz_y_z, g_z_0_y_0_y_xz_0_x, g_z_0_y_0_y_xz_0_y, g_z_0_y_0_y_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_xz_0_x[i] = 4.0 * g_yz_xz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xz_0_y[i] = 4.0 * g_yz_xz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_xz_0_z[i] = 4.0 * g_yz_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (405-408)

    #pragma omp simd aligned(g_yz_yy_y_x, g_yz_yy_y_y, g_yz_yy_y_z, g_z_0_y_0_y_yy_0_x, g_z_0_y_0_y_yy_0_y, g_z_0_y_0_y_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_yy_0_x[i] = 4.0 * g_yz_yy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yy_0_y[i] = 4.0 * g_yz_yy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yy_0_z[i] = 4.0 * g_yz_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (408-411)

    #pragma omp simd aligned(g_yz_yz_y_x, g_yz_yz_y_y, g_yz_yz_y_z, g_z_0_y_0_y_yz_0_x, g_z_0_y_0_y_yz_0_y, g_z_0_y_0_y_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_yz_0_x[i] = 4.0 * g_yz_yz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yz_0_y[i] = 4.0 * g_yz_yz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_yz_0_z[i] = 4.0 * g_yz_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (411-414)

    #pragma omp simd aligned(g_yz_zz_y_x, g_yz_zz_y_y, g_yz_zz_y_z, g_z_0_y_0_y_zz_0_x, g_z_0_y_0_y_zz_0_y, g_z_0_y_0_y_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_zz_0_x[i] = 4.0 * g_yz_zz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_zz_0_y[i] = 4.0 * g_yz_zz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_zz_0_z[i] = 4.0 * g_yz_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (414-417)

    #pragma omp simd aligned(g_0_xx_y_x, g_0_xx_y_y, g_0_xx_y_z, g_z_0_y_0_z_xx_0_x, g_z_0_y_0_z_xx_0_y, g_z_0_y_0_z_xx_0_z, g_zz_xx_y_x, g_zz_xx_y_y, g_zz_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_xx_0_x[i] = -2.0 * g_0_xx_y_x[i] * c_exps[i] + 4.0 * g_zz_xx_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xx_0_y[i] = -2.0 * g_0_xx_y_y[i] * c_exps[i] + 4.0 * g_zz_xx_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xx_0_z[i] = -2.0 * g_0_xx_y_z[i] * c_exps[i] + 4.0 * g_zz_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (417-420)

    #pragma omp simd aligned(g_0_xy_y_x, g_0_xy_y_y, g_0_xy_y_z, g_z_0_y_0_z_xy_0_x, g_z_0_y_0_z_xy_0_y, g_z_0_y_0_z_xy_0_z, g_zz_xy_y_x, g_zz_xy_y_y, g_zz_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_xy_0_x[i] = -2.0 * g_0_xy_y_x[i] * c_exps[i] + 4.0 * g_zz_xy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xy_0_y[i] = -2.0 * g_0_xy_y_y[i] * c_exps[i] + 4.0 * g_zz_xy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xy_0_z[i] = -2.0 * g_0_xy_y_z[i] * c_exps[i] + 4.0 * g_zz_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (420-423)

    #pragma omp simd aligned(g_0_xz_y_x, g_0_xz_y_y, g_0_xz_y_z, g_z_0_y_0_z_xz_0_x, g_z_0_y_0_z_xz_0_y, g_z_0_y_0_z_xz_0_z, g_zz_xz_y_x, g_zz_xz_y_y, g_zz_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_xz_0_x[i] = -2.0 * g_0_xz_y_x[i] * c_exps[i] + 4.0 * g_zz_xz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xz_0_y[i] = -2.0 * g_0_xz_y_y[i] * c_exps[i] + 4.0 * g_zz_xz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_xz_0_z[i] = -2.0 * g_0_xz_y_z[i] * c_exps[i] + 4.0 * g_zz_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (423-426)

    #pragma omp simd aligned(g_0_yy_y_x, g_0_yy_y_y, g_0_yy_y_z, g_z_0_y_0_z_yy_0_x, g_z_0_y_0_z_yy_0_y, g_z_0_y_0_z_yy_0_z, g_zz_yy_y_x, g_zz_yy_y_y, g_zz_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_yy_0_x[i] = -2.0 * g_0_yy_y_x[i] * c_exps[i] + 4.0 * g_zz_yy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yy_0_y[i] = -2.0 * g_0_yy_y_y[i] * c_exps[i] + 4.0 * g_zz_yy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yy_0_z[i] = -2.0 * g_0_yy_y_z[i] * c_exps[i] + 4.0 * g_zz_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (426-429)

    #pragma omp simd aligned(g_0_yz_y_x, g_0_yz_y_y, g_0_yz_y_z, g_z_0_y_0_z_yz_0_x, g_z_0_y_0_z_yz_0_y, g_z_0_y_0_z_yz_0_z, g_zz_yz_y_x, g_zz_yz_y_y, g_zz_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_yz_0_x[i] = -2.0 * g_0_yz_y_x[i] * c_exps[i] + 4.0 * g_zz_yz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yz_0_y[i] = -2.0 * g_0_yz_y_y[i] * c_exps[i] + 4.0 * g_zz_yz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_yz_0_z[i] = -2.0 * g_0_yz_y_z[i] * c_exps[i] + 4.0 * g_zz_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (429-432)

    #pragma omp simd aligned(g_0_zz_y_x, g_0_zz_y_y, g_0_zz_y_z, g_z_0_y_0_z_zz_0_x, g_z_0_y_0_z_zz_0_y, g_z_0_y_0_z_zz_0_z, g_zz_zz_y_x, g_zz_zz_y_y, g_zz_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_zz_0_x[i] = -2.0 * g_0_zz_y_x[i] * c_exps[i] + 4.0 * g_zz_zz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_zz_0_y[i] = -2.0 * g_0_zz_y_y[i] * c_exps[i] + 4.0 * g_zz_zz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_zz_0_z[i] = -2.0 * g_0_zz_y_z[i] * c_exps[i] + 4.0 * g_zz_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (432-435)

    #pragma omp simd aligned(g_xz_xx_z_x, g_xz_xx_z_y, g_xz_xx_z_z, g_z_0_z_0_x_xx_0_x, g_z_0_z_0_x_xx_0_y, g_z_0_z_0_x_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_xx_0_x[i] = 4.0 * g_xz_xx_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xx_0_y[i] = 4.0 * g_xz_xx_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xx_0_z[i] = 4.0 * g_xz_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (435-438)

    #pragma omp simd aligned(g_xz_xy_z_x, g_xz_xy_z_y, g_xz_xy_z_z, g_z_0_z_0_x_xy_0_x, g_z_0_z_0_x_xy_0_y, g_z_0_z_0_x_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_xy_0_x[i] = 4.0 * g_xz_xy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xy_0_y[i] = 4.0 * g_xz_xy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xy_0_z[i] = 4.0 * g_xz_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (438-441)

    #pragma omp simd aligned(g_xz_xz_z_x, g_xz_xz_z_y, g_xz_xz_z_z, g_z_0_z_0_x_xz_0_x, g_z_0_z_0_x_xz_0_y, g_z_0_z_0_x_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_xz_0_x[i] = 4.0 * g_xz_xz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xz_0_y[i] = 4.0 * g_xz_xz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_xz_0_z[i] = 4.0 * g_xz_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (441-444)

    #pragma omp simd aligned(g_xz_yy_z_x, g_xz_yy_z_y, g_xz_yy_z_z, g_z_0_z_0_x_yy_0_x, g_z_0_z_0_x_yy_0_y, g_z_0_z_0_x_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_yy_0_x[i] = 4.0 * g_xz_yy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yy_0_y[i] = 4.0 * g_xz_yy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yy_0_z[i] = 4.0 * g_xz_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (444-447)

    #pragma omp simd aligned(g_xz_yz_z_x, g_xz_yz_z_y, g_xz_yz_z_z, g_z_0_z_0_x_yz_0_x, g_z_0_z_0_x_yz_0_y, g_z_0_z_0_x_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_yz_0_x[i] = 4.0 * g_xz_yz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yz_0_y[i] = 4.0 * g_xz_yz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_yz_0_z[i] = 4.0 * g_xz_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (447-450)

    #pragma omp simd aligned(g_xz_zz_z_x, g_xz_zz_z_y, g_xz_zz_z_z, g_z_0_z_0_x_zz_0_x, g_z_0_z_0_x_zz_0_y, g_z_0_z_0_x_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_zz_0_x[i] = 4.0 * g_xz_zz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_zz_0_y[i] = 4.0 * g_xz_zz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_zz_0_z[i] = 4.0 * g_xz_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (450-453)

    #pragma omp simd aligned(g_yz_xx_z_x, g_yz_xx_z_y, g_yz_xx_z_z, g_z_0_z_0_y_xx_0_x, g_z_0_z_0_y_xx_0_y, g_z_0_z_0_y_xx_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_xx_0_x[i] = 4.0 * g_yz_xx_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xx_0_y[i] = 4.0 * g_yz_xx_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xx_0_z[i] = 4.0 * g_yz_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (453-456)

    #pragma omp simd aligned(g_yz_xy_z_x, g_yz_xy_z_y, g_yz_xy_z_z, g_z_0_z_0_y_xy_0_x, g_z_0_z_0_y_xy_0_y, g_z_0_z_0_y_xy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_xy_0_x[i] = 4.0 * g_yz_xy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xy_0_y[i] = 4.0 * g_yz_xy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xy_0_z[i] = 4.0 * g_yz_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (456-459)

    #pragma omp simd aligned(g_yz_xz_z_x, g_yz_xz_z_y, g_yz_xz_z_z, g_z_0_z_0_y_xz_0_x, g_z_0_z_0_y_xz_0_y, g_z_0_z_0_y_xz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_xz_0_x[i] = 4.0 * g_yz_xz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xz_0_y[i] = 4.0 * g_yz_xz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_xz_0_z[i] = 4.0 * g_yz_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (459-462)

    #pragma omp simd aligned(g_yz_yy_z_x, g_yz_yy_z_y, g_yz_yy_z_z, g_z_0_z_0_y_yy_0_x, g_z_0_z_0_y_yy_0_y, g_z_0_z_0_y_yy_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_yy_0_x[i] = 4.0 * g_yz_yy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yy_0_y[i] = 4.0 * g_yz_yy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yy_0_z[i] = 4.0 * g_yz_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (462-465)

    #pragma omp simd aligned(g_yz_yz_z_x, g_yz_yz_z_y, g_yz_yz_z_z, g_z_0_z_0_y_yz_0_x, g_z_0_z_0_y_yz_0_y, g_z_0_z_0_y_yz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_yz_0_x[i] = 4.0 * g_yz_yz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yz_0_y[i] = 4.0 * g_yz_yz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_yz_0_z[i] = 4.0 * g_yz_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (465-468)

    #pragma omp simd aligned(g_yz_zz_z_x, g_yz_zz_z_y, g_yz_zz_z_z, g_z_0_z_0_y_zz_0_x, g_z_0_z_0_y_zz_0_y, g_z_0_z_0_y_zz_0_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_zz_0_x[i] = 4.0 * g_yz_zz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_zz_0_y[i] = 4.0 * g_yz_zz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_zz_0_z[i] = 4.0 * g_yz_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (468-471)

    #pragma omp simd aligned(g_0_xx_z_x, g_0_xx_z_y, g_0_xx_z_z, g_z_0_z_0_z_xx_0_x, g_z_0_z_0_z_xx_0_y, g_z_0_z_0_z_xx_0_z, g_zz_xx_z_x, g_zz_xx_z_y, g_zz_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_xx_0_x[i] = -2.0 * g_0_xx_z_x[i] * c_exps[i] + 4.0 * g_zz_xx_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xx_0_y[i] = -2.0 * g_0_xx_z_y[i] * c_exps[i] + 4.0 * g_zz_xx_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xx_0_z[i] = -2.0 * g_0_xx_z_z[i] * c_exps[i] + 4.0 * g_zz_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (471-474)

    #pragma omp simd aligned(g_0_xy_z_x, g_0_xy_z_y, g_0_xy_z_z, g_z_0_z_0_z_xy_0_x, g_z_0_z_0_z_xy_0_y, g_z_0_z_0_z_xy_0_z, g_zz_xy_z_x, g_zz_xy_z_y, g_zz_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_xy_0_x[i] = -2.0 * g_0_xy_z_x[i] * c_exps[i] + 4.0 * g_zz_xy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xy_0_y[i] = -2.0 * g_0_xy_z_y[i] * c_exps[i] + 4.0 * g_zz_xy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xy_0_z[i] = -2.0 * g_0_xy_z_z[i] * c_exps[i] + 4.0 * g_zz_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (474-477)

    #pragma omp simd aligned(g_0_xz_z_x, g_0_xz_z_y, g_0_xz_z_z, g_z_0_z_0_z_xz_0_x, g_z_0_z_0_z_xz_0_y, g_z_0_z_0_z_xz_0_z, g_zz_xz_z_x, g_zz_xz_z_y, g_zz_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_xz_0_x[i] = -2.0 * g_0_xz_z_x[i] * c_exps[i] + 4.0 * g_zz_xz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xz_0_y[i] = -2.0 * g_0_xz_z_y[i] * c_exps[i] + 4.0 * g_zz_xz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_xz_0_z[i] = -2.0 * g_0_xz_z_z[i] * c_exps[i] + 4.0 * g_zz_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (477-480)

    #pragma omp simd aligned(g_0_yy_z_x, g_0_yy_z_y, g_0_yy_z_z, g_z_0_z_0_z_yy_0_x, g_z_0_z_0_z_yy_0_y, g_z_0_z_0_z_yy_0_z, g_zz_yy_z_x, g_zz_yy_z_y, g_zz_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_yy_0_x[i] = -2.0 * g_0_yy_z_x[i] * c_exps[i] + 4.0 * g_zz_yy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yy_0_y[i] = -2.0 * g_0_yy_z_y[i] * c_exps[i] + 4.0 * g_zz_yy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yy_0_z[i] = -2.0 * g_0_yy_z_z[i] * c_exps[i] + 4.0 * g_zz_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (480-483)

    #pragma omp simd aligned(g_0_yz_z_x, g_0_yz_z_y, g_0_yz_z_z, g_z_0_z_0_z_yz_0_x, g_z_0_z_0_z_yz_0_y, g_z_0_z_0_z_yz_0_z, g_zz_yz_z_x, g_zz_yz_z_y, g_zz_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_yz_0_x[i] = -2.0 * g_0_yz_z_x[i] * c_exps[i] + 4.0 * g_zz_yz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yz_0_y[i] = -2.0 * g_0_yz_z_y[i] * c_exps[i] + 4.0 * g_zz_yz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_yz_0_z[i] = -2.0 * g_0_yz_z_z[i] * c_exps[i] + 4.0 * g_zz_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (483-486)

    #pragma omp simd aligned(g_0_zz_z_x, g_0_zz_z_y, g_0_zz_z_z, g_z_0_z_0_z_zz_0_x, g_z_0_z_0_z_zz_0_y, g_z_0_z_0_z_zz_0_z, g_zz_zz_z_x, g_zz_zz_z_y, g_zz_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_zz_0_x[i] = -2.0 * g_0_zz_z_x[i] * c_exps[i] + 4.0 * g_zz_zz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_zz_0_y[i] = -2.0 * g_0_zz_z_y[i] * c_exps[i] + 4.0 * g_zz_zz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_zz_0_z[i] = -2.0 * g_0_zz_z_z[i] * c_exps[i] + 4.0 * g_zz_zz_z_z[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace

