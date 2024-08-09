#include "GeomDeriv1010OfScalarForSDPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_sdpp_0(CSimdArray<double>& buffer_1010_sdpp,
                     const CSimdArray<double>& buffer_pdsp,
                     const CSimdArray<double>& buffer_pddp,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_sdpp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pdsp

    auto g_x_xx_0_x = buffer_pdsp[0];

    auto g_x_xx_0_y = buffer_pdsp[1];

    auto g_x_xx_0_z = buffer_pdsp[2];

    auto g_x_xy_0_x = buffer_pdsp[3];

    auto g_x_xy_0_y = buffer_pdsp[4];

    auto g_x_xy_0_z = buffer_pdsp[5];

    auto g_x_xz_0_x = buffer_pdsp[6];

    auto g_x_xz_0_y = buffer_pdsp[7];

    auto g_x_xz_0_z = buffer_pdsp[8];

    auto g_x_yy_0_x = buffer_pdsp[9];

    auto g_x_yy_0_y = buffer_pdsp[10];

    auto g_x_yy_0_z = buffer_pdsp[11];

    auto g_x_yz_0_x = buffer_pdsp[12];

    auto g_x_yz_0_y = buffer_pdsp[13];

    auto g_x_yz_0_z = buffer_pdsp[14];

    auto g_x_zz_0_x = buffer_pdsp[15];

    auto g_x_zz_0_y = buffer_pdsp[16];

    auto g_x_zz_0_z = buffer_pdsp[17];

    auto g_y_xx_0_x = buffer_pdsp[18];

    auto g_y_xx_0_y = buffer_pdsp[19];

    auto g_y_xx_0_z = buffer_pdsp[20];

    auto g_y_xy_0_x = buffer_pdsp[21];

    auto g_y_xy_0_y = buffer_pdsp[22];

    auto g_y_xy_0_z = buffer_pdsp[23];

    auto g_y_xz_0_x = buffer_pdsp[24];

    auto g_y_xz_0_y = buffer_pdsp[25];

    auto g_y_xz_0_z = buffer_pdsp[26];

    auto g_y_yy_0_x = buffer_pdsp[27];

    auto g_y_yy_0_y = buffer_pdsp[28];

    auto g_y_yy_0_z = buffer_pdsp[29];

    auto g_y_yz_0_x = buffer_pdsp[30];

    auto g_y_yz_0_y = buffer_pdsp[31];

    auto g_y_yz_0_z = buffer_pdsp[32];

    auto g_y_zz_0_x = buffer_pdsp[33];

    auto g_y_zz_0_y = buffer_pdsp[34];

    auto g_y_zz_0_z = buffer_pdsp[35];

    auto g_z_xx_0_x = buffer_pdsp[36];

    auto g_z_xx_0_y = buffer_pdsp[37];

    auto g_z_xx_0_z = buffer_pdsp[38];

    auto g_z_xy_0_x = buffer_pdsp[39];

    auto g_z_xy_0_y = buffer_pdsp[40];

    auto g_z_xy_0_z = buffer_pdsp[41];

    auto g_z_xz_0_x = buffer_pdsp[42];

    auto g_z_xz_0_y = buffer_pdsp[43];

    auto g_z_xz_0_z = buffer_pdsp[44];

    auto g_z_yy_0_x = buffer_pdsp[45];

    auto g_z_yy_0_y = buffer_pdsp[46];

    auto g_z_yy_0_z = buffer_pdsp[47];

    auto g_z_yz_0_x = buffer_pdsp[48];

    auto g_z_yz_0_y = buffer_pdsp[49];

    auto g_z_yz_0_z = buffer_pdsp[50];

    auto g_z_zz_0_x = buffer_pdsp[51];

    auto g_z_zz_0_y = buffer_pdsp[52];

    auto g_z_zz_0_z = buffer_pdsp[53];

    /// Set up components of auxilary buffer : buffer_pddp

    auto g_x_xx_xx_x = buffer_pddp[0];

    auto g_x_xx_xx_y = buffer_pddp[1];

    auto g_x_xx_xx_z = buffer_pddp[2];

    auto g_x_xx_xy_x = buffer_pddp[3];

    auto g_x_xx_xy_y = buffer_pddp[4];

    auto g_x_xx_xy_z = buffer_pddp[5];

    auto g_x_xx_xz_x = buffer_pddp[6];

    auto g_x_xx_xz_y = buffer_pddp[7];

    auto g_x_xx_xz_z = buffer_pddp[8];

    auto g_x_xx_yy_x = buffer_pddp[9];

    auto g_x_xx_yy_y = buffer_pddp[10];

    auto g_x_xx_yy_z = buffer_pddp[11];

    auto g_x_xx_yz_x = buffer_pddp[12];

    auto g_x_xx_yz_y = buffer_pddp[13];

    auto g_x_xx_yz_z = buffer_pddp[14];

    auto g_x_xx_zz_x = buffer_pddp[15];

    auto g_x_xx_zz_y = buffer_pddp[16];

    auto g_x_xx_zz_z = buffer_pddp[17];

    auto g_x_xy_xx_x = buffer_pddp[18];

    auto g_x_xy_xx_y = buffer_pddp[19];

    auto g_x_xy_xx_z = buffer_pddp[20];

    auto g_x_xy_xy_x = buffer_pddp[21];

    auto g_x_xy_xy_y = buffer_pddp[22];

    auto g_x_xy_xy_z = buffer_pddp[23];

    auto g_x_xy_xz_x = buffer_pddp[24];

    auto g_x_xy_xz_y = buffer_pddp[25];

    auto g_x_xy_xz_z = buffer_pddp[26];

    auto g_x_xy_yy_x = buffer_pddp[27];

    auto g_x_xy_yy_y = buffer_pddp[28];

    auto g_x_xy_yy_z = buffer_pddp[29];

    auto g_x_xy_yz_x = buffer_pddp[30];

    auto g_x_xy_yz_y = buffer_pddp[31];

    auto g_x_xy_yz_z = buffer_pddp[32];

    auto g_x_xy_zz_x = buffer_pddp[33];

    auto g_x_xy_zz_y = buffer_pddp[34];

    auto g_x_xy_zz_z = buffer_pddp[35];

    auto g_x_xz_xx_x = buffer_pddp[36];

    auto g_x_xz_xx_y = buffer_pddp[37];

    auto g_x_xz_xx_z = buffer_pddp[38];

    auto g_x_xz_xy_x = buffer_pddp[39];

    auto g_x_xz_xy_y = buffer_pddp[40];

    auto g_x_xz_xy_z = buffer_pddp[41];

    auto g_x_xz_xz_x = buffer_pddp[42];

    auto g_x_xz_xz_y = buffer_pddp[43];

    auto g_x_xz_xz_z = buffer_pddp[44];

    auto g_x_xz_yy_x = buffer_pddp[45];

    auto g_x_xz_yy_y = buffer_pddp[46];

    auto g_x_xz_yy_z = buffer_pddp[47];

    auto g_x_xz_yz_x = buffer_pddp[48];

    auto g_x_xz_yz_y = buffer_pddp[49];

    auto g_x_xz_yz_z = buffer_pddp[50];

    auto g_x_xz_zz_x = buffer_pddp[51];

    auto g_x_xz_zz_y = buffer_pddp[52];

    auto g_x_xz_zz_z = buffer_pddp[53];

    auto g_x_yy_xx_x = buffer_pddp[54];

    auto g_x_yy_xx_y = buffer_pddp[55];

    auto g_x_yy_xx_z = buffer_pddp[56];

    auto g_x_yy_xy_x = buffer_pddp[57];

    auto g_x_yy_xy_y = buffer_pddp[58];

    auto g_x_yy_xy_z = buffer_pddp[59];

    auto g_x_yy_xz_x = buffer_pddp[60];

    auto g_x_yy_xz_y = buffer_pddp[61];

    auto g_x_yy_xz_z = buffer_pddp[62];

    auto g_x_yy_yy_x = buffer_pddp[63];

    auto g_x_yy_yy_y = buffer_pddp[64];

    auto g_x_yy_yy_z = buffer_pddp[65];

    auto g_x_yy_yz_x = buffer_pddp[66];

    auto g_x_yy_yz_y = buffer_pddp[67];

    auto g_x_yy_yz_z = buffer_pddp[68];

    auto g_x_yy_zz_x = buffer_pddp[69];

    auto g_x_yy_zz_y = buffer_pddp[70];

    auto g_x_yy_zz_z = buffer_pddp[71];

    auto g_x_yz_xx_x = buffer_pddp[72];

    auto g_x_yz_xx_y = buffer_pddp[73];

    auto g_x_yz_xx_z = buffer_pddp[74];

    auto g_x_yz_xy_x = buffer_pddp[75];

    auto g_x_yz_xy_y = buffer_pddp[76];

    auto g_x_yz_xy_z = buffer_pddp[77];

    auto g_x_yz_xz_x = buffer_pddp[78];

    auto g_x_yz_xz_y = buffer_pddp[79];

    auto g_x_yz_xz_z = buffer_pddp[80];

    auto g_x_yz_yy_x = buffer_pddp[81];

    auto g_x_yz_yy_y = buffer_pddp[82];

    auto g_x_yz_yy_z = buffer_pddp[83];

    auto g_x_yz_yz_x = buffer_pddp[84];

    auto g_x_yz_yz_y = buffer_pddp[85];

    auto g_x_yz_yz_z = buffer_pddp[86];

    auto g_x_yz_zz_x = buffer_pddp[87];

    auto g_x_yz_zz_y = buffer_pddp[88];

    auto g_x_yz_zz_z = buffer_pddp[89];

    auto g_x_zz_xx_x = buffer_pddp[90];

    auto g_x_zz_xx_y = buffer_pddp[91];

    auto g_x_zz_xx_z = buffer_pddp[92];

    auto g_x_zz_xy_x = buffer_pddp[93];

    auto g_x_zz_xy_y = buffer_pddp[94];

    auto g_x_zz_xy_z = buffer_pddp[95];

    auto g_x_zz_xz_x = buffer_pddp[96];

    auto g_x_zz_xz_y = buffer_pddp[97];

    auto g_x_zz_xz_z = buffer_pddp[98];

    auto g_x_zz_yy_x = buffer_pddp[99];

    auto g_x_zz_yy_y = buffer_pddp[100];

    auto g_x_zz_yy_z = buffer_pddp[101];

    auto g_x_zz_yz_x = buffer_pddp[102];

    auto g_x_zz_yz_y = buffer_pddp[103];

    auto g_x_zz_yz_z = buffer_pddp[104];

    auto g_x_zz_zz_x = buffer_pddp[105];

    auto g_x_zz_zz_y = buffer_pddp[106];

    auto g_x_zz_zz_z = buffer_pddp[107];

    auto g_y_xx_xx_x = buffer_pddp[108];

    auto g_y_xx_xx_y = buffer_pddp[109];

    auto g_y_xx_xx_z = buffer_pddp[110];

    auto g_y_xx_xy_x = buffer_pddp[111];

    auto g_y_xx_xy_y = buffer_pddp[112];

    auto g_y_xx_xy_z = buffer_pddp[113];

    auto g_y_xx_xz_x = buffer_pddp[114];

    auto g_y_xx_xz_y = buffer_pddp[115];

    auto g_y_xx_xz_z = buffer_pddp[116];

    auto g_y_xx_yy_x = buffer_pddp[117];

    auto g_y_xx_yy_y = buffer_pddp[118];

    auto g_y_xx_yy_z = buffer_pddp[119];

    auto g_y_xx_yz_x = buffer_pddp[120];

    auto g_y_xx_yz_y = buffer_pddp[121];

    auto g_y_xx_yz_z = buffer_pddp[122];

    auto g_y_xx_zz_x = buffer_pddp[123];

    auto g_y_xx_zz_y = buffer_pddp[124];

    auto g_y_xx_zz_z = buffer_pddp[125];

    auto g_y_xy_xx_x = buffer_pddp[126];

    auto g_y_xy_xx_y = buffer_pddp[127];

    auto g_y_xy_xx_z = buffer_pddp[128];

    auto g_y_xy_xy_x = buffer_pddp[129];

    auto g_y_xy_xy_y = buffer_pddp[130];

    auto g_y_xy_xy_z = buffer_pddp[131];

    auto g_y_xy_xz_x = buffer_pddp[132];

    auto g_y_xy_xz_y = buffer_pddp[133];

    auto g_y_xy_xz_z = buffer_pddp[134];

    auto g_y_xy_yy_x = buffer_pddp[135];

    auto g_y_xy_yy_y = buffer_pddp[136];

    auto g_y_xy_yy_z = buffer_pddp[137];

    auto g_y_xy_yz_x = buffer_pddp[138];

    auto g_y_xy_yz_y = buffer_pddp[139];

    auto g_y_xy_yz_z = buffer_pddp[140];

    auto g_y_xy_zz_x = buffer_pddp[141];

    auto g_y_xy_zz_y = buffer_pddp[142];

    auto g_y_xy_zz_z = buffer_pddp[143];

    auto g_y_xz_xx_x = buffer_pddp[144];

    auto g_y_xz_xx_y = buffer_pddp[145];

    auto g_y_xz_xx_z = buffer_pddp[146];

    auto g_y_xz_xy_x = buffer_pddp[147];

    auto g_y_xz_xy_y = buffer_pddp[148];

    auto g_y_xz_xy_z = buffer_pddp[149];

    auto g_y_xz_xz_x = buffer_pddp[150];

    auto g_y_xz_xz_y = buffer_pddp[151];

    auto g_y_xz_xz_z = buffer_pddp[152];

    auto g_y_xz_yy_x = buffer_pddp[153];

    auto g_y_xz_yy_y = buffer_pddp[154];

    auto g_y_xz_yy_z = buffer_pddp[155];

    auto g_y_xz_yz_x = buffer_pddp[156];

    auto g_y_xz_yz_y = buffer_pddp[157];

    auto g_y_xz_yz_z = buffer_pddp[158];

    auto g_y_xz_zz_x = buffer_pddp[159];

    auto g_y_xz_zz_y = buffer_pddp[160];

    auto g_y_xz_zz_z = buffer_pddp[161];

    auto g_y_yy_xx_x = buffer_pddp[162];

    auto g_y_yy_xx_y = buffer_pddp[163];

    auto g_y_yy_xx_z = buffer_pddp[164];

    auto g_y_yy_xy_x = buffer_pddp[165];

    auto g_y_yy_xy_y = buffer_pddp[166];

    auto g_y_yy_xy_z = buffer_pddp[167];

    auto g_y_yy_xz_x = buffer_pddp[168];

    auto g_y_yy_xz_y = buffer_pddp[169];

    auto g_y_yy_xz_z = buffer_pddp[170];

    auto g_y_yy_yy_x = buffer_pddp[171];

    auto g_y_yy_yy_y = buffer_pddp[172];

    auto g_y_yy_yy_z = buffer_pddp[173];

    auto g_y_yy_yz_x = buffer_pddp[174];

    auto g_y_yy_yz_y = buffer_pddp[175];

    auto g_y_yy_yz_z = buffer_pddp[176];

    auto g_y_yy_zz_x = buffer_pddp[177];

    auto g_y_yy_zz_y = buffer_pddp[178];

    auto g_y_yy_zz_z = buffer_pddp[179];

    auto g_y_yz_xx_x = buffer_pddp[180];

    auto g_y_yz_xx_y = buffer_pddp[181];

    auto g_y_yz_xx_z = buffer_pddp[182];

    auto g_y_yz_xy_x = buffer_pddp[183];

    auto g_y_yz_xy_y = buffer_pddp[184];

    auto g_y_yz_xy_z = buffer_pddp[185];

    auto g_y_yz_xz_x = buffer_pddp[186];

    auto g_y_yz_xz_y = buffer_pddp[187];

    auto g_y_yz_xz_z = buffer_pddp[188];

    auto g_y_yz_yy_x = buffer_pddp[189];

    auto g_y_yz_yy_y = buffer_pddp[190];

    auto g_y_yz_yy_z = buffer_pddp[191];

    auto g_y_yz_yz_x = buffer_pddp[192];

    auto g_y_yz_yz_y = buffer_pddp[193];

    auto g_y_yz_yz_z = buffer_pddp[194];

    auto g_y_yz_zz_x = buffer_pddp[195];

    auto g_y_yz_zz_y = buffer_pddp[196];

    auto g_y_yz_zz_z = buffer_pddp[197];

    auto g_y_zz_xx_x = buffer_pddp[198];

    auto g_y_zz_xx_y = buffer_pddp[199];

    auto g_y_zz_xx_z = buffer_pddp[200];

    auto g_y_zz_xy_x = buffer_pddp[201];

    auto g_y_zz_xy_y = buffer_pddp[202];

    auto g_y_zz_xy_z = buffer_pddp[203];

    auto g_y_zz_xz_x = buffer_pddp[204];

    auto g_y_zz_xz_y = buffer_pddp[205];

    auto g_y_zz_xz_z = buffer_pddp[206];

    auto g_y_zz_yy_x = buffer_pddp[207];

    auto g_y_zz_yy_y = buffer_pddp[208];

    auto g_y_zz_yy_z = buffer_pddp[209];

    auto g_y_zz_yz_x = buffer_pddp[210];

    auto g_y_zz_yz_y = buffer_pddp[211];

    auto g_y_zz_yz_z = buffer_pddp[212];

    auto g_y_zz_zz_x = buffer_pddp[213];

    auto g_y_zz_zz_y = buffer_pddp[214];

    auto g_y_zz_zz_z = buffer_pddp[215];

    auto g_z_xx_xx_x = buffer_pddp[216];

    auto g_z_xx_xx_y = buffer_pddp[217];

    auto g_z_xx_xx_z = buffer_pddp[218];

    auto g_z_xx_xy_x = buffer_pddp[219];

    auto g_z_xx_xy_y = buffer_pddp[220];

    auto g_z_xx_xy_z = buffer_pddp[221];

    auto g_z_xx_xz_x = buffer_pddp[222];

    auto g_z_xx_xz_y = buffer_pddp[223];

    auto g_z_xx_xz_z = buffer_pddp[224];

    auto g_z_xx_yy_x = buffer_pddp[225];

    auto g_z_xx_yy_y = buffer_pddp[226];

    auto g_z_xx_yy_z = buffer_pddp[227];

    auto g_z_xx_yz_x = buffer_pddp[228];

    auto g_z_xx_yz_y = buffer_pddp[229];

    auto g_z_xx_yz_z = buffer_pddp[230];

    auto g_z_xx_zz_x = buffer_pddp[231];

    auto g_z_xx_zz_y = buffer_pddp[232];

    auto g_z_xx_zz_z = buffer_pddp[233];

    auto g_z_xy_xx_x = buffer_pddp[234];

    auto g_z_xy_xx_y = buffer_pddp[235];

    auto g_z_xy_xx_z = buffer_pddp[236];

    auto g_z_xy_xy_x = buffer_pddp[237];

    auto g_z_xy_xy_y = buffer_pddp[238];

    auto g_z_xy_xy_z = buffer_pddp[239];

    auto g_z_xy_xz_x = buffer_pddp[240];

    auto g_z_xy_xz_y = buffer_pddp[241];

    auto g_z_xy_xz_z = buffer_pddp[242];

    auto g_z_xy_yy_x = buffer_pddp[243];

    auto g_z_xy_yy_y = buffer_pddp[244];

    auto g_z_xy_yy_z = buffer_pddp[245];

    auto g_z_xy_yz_x = buffer_pddp[246];

    auto g_z_xy_yz_y = buffer_pddp[247];

    auto g_z_xy_yz_z = buffer_pddp[248];

    auto g_z_xy_zz_x = buffer_pddp[249];

    auto g_z_xy_zz_y = buffer_pddp[250];

    auto g_z_xy_zz_z = buffer_pddp[251];

    auto g_z_xz_xx_x = buffer_pddp[252];

    auto g_z_xz_xx_y = buffer_pddp[253];

    auto g_z_xz_xx_z = buffer_pddp[254];

    auto g_z_xz_xy_x = buffer_pddp[255];

    auto g_z_xz_xy_y = buffer_pddp[256];

    auto g_z_xz_xy_z = buffer_pddp[257];

    auto g_z_xz_xz_x = buffer_pddp[258];

    auto g_z_xz_xz_y = buffer_pddp[259];

    auto g_z_xz_xz_z = buffer_pddp[260];

    auto g_z_xz_yy_x = buffer_pddp[261];

    auto g_z_xz_yy_y = buffer_pddp[262];

    auto g_z_xz_yy_z = buffer_pddp[263];

    auto g_z_xz_yz_x = buffer_pddp[264];

    auto g_z_xz_yz_y = buffer_pddp[265];

    auto g_z_xz_yz_z = buffer_pddp[266];

    auto g_z_xz_zz_x = buffer_pddp[267];

    auto g_z_xz_zz_y = buffer_pddp[268];

    auto g_z_xz_zz_z = buffer_pddp[269];

    auto g_z_yy_xx_x = buffer_pddp[270];

    auto g_z_yy_xx_y = buffer_pddp[271];

    auto g_z_yy_xx_z = buffer_pddp[272];

    auto g_z_yy_xy_x = buffer_pddp[273];

    auto g_z_yy_xy_y = buffer_pddp[274];

    auto g_z_yy_xy_z = buffer_pddp[275];

    auto g_z_yy_xz_x = buffer_pddp[276];

    auto g_z_yy_xz_y = buffer_pddp[277];

    auto g_z_yy_xz_z = buffer_pddp[278];

    auto g_z_yy_yy_x = buffer_pddp[279];

    auto g_z_yy_yy_y = buffer_pddp[280];

    auto g_z_yy_yy_z = buffer_pddp[281];

    auto g_z_yy_yz_x = buffer_pddp[282];

    auto g_z_yy_yz_y = buffer_pddp[283];

    auto g_z_yy_yz_z = buffer_pddp[284];

    auto g_z_yy_zz_x = buffer_pddp[285];

    auto g_z_yy_zz_y = buffer_pddp[286];

    auto g_z_yy_zz_z = buffer_pddp[287];

    auto g_z_yz_xx_x = buffer_pddp[288];

    auto g_z_yz_xx_y = buffer_pddp[289];

    auto g_z_yz_xx_z = buffer_pddp[290];

    auto g_z_yz_xy_x = buffer_pddp[291];

    auto g_z_yz_xy_y = buffer_pddp[292];

    auto g_z_yz_xy_z = buffer_pddp[293];

    auto g_z_yz_xz_x = buffer_pddp[294];

    auto g_z_yz_xz_y = buffer_pddp[295];

    auto g_z_yz_xz_z = buffer_pddp[296];

    auto g_z_yz_yy_x = buffer_pddp[297];

    auto g_z_yz_yy_y = buffer_pddp[298];

    auto g_z_yz_yy_z = buffer_pddp[299];

    auto g_z_yz_yz_x = buffer_pddp[300];

    auto g_z_yz_yz_y = buffer_pddp[301];

    auto g_z_yz_yz_z = buffer_pddp[302];

    auto g_z_yz_zz_x = buffer_pddp[303];

    auto g_z_yz_zz_y = buffer_pddp[304];

    auto g_z_yz_zz_z = buffer_pddp[305];

    auto g_z_zz_xx_x = buffer_pddp[306];

    auto g_z_zz_xx_y = buffer_pddp[307];

    auto g_z_zz_xx_z = buffer_pddp[308];

    auto g_z_zz_xy_x = buffer_pddp[309];

    auto g_z_zz_xy_y = buffer_pddp[310];

    auto g_z_zz_xy_z = buffer_pddp[311];

    auto g_z_zz_xz_x = buffer_pddp[312];

    auto g_z_zz_xz_y = buffer_pddp[313];

    auto g_z_zz_xz_z = buffer_pddp[314];

    auto g_z_zz_yy_x = buffer_pddp[315];

    auto g_z_zz_yy_y = buffer_pddp[316];

    auto g_z_zz_yy_z = buffer_pddp[317];

    auto g_z_zz_yz_x = buffer_pddp[318];

    auto g_z_zz_yz_y = buffer_pddp[319];

    auto g_z_zz_yz_z = buffer_pddp[320];

    auto g_z_zz_zz_x = buffer_pddp[321];

    auto g_z_zz_zz_y = buffer_pddp[322];

    auto g_z_zz_zz_z = buffer_pddp[323];

    /// Set up components of integrals buffer : buffer_1010_sdpp

    auto g_x_0_x_0_0_xx_x_x = buffer_1010_sdpp[0];

    auto g_x_0_x_0_0_xx_x_y = buffer_1010_sdpp[1];

    auto g_x_0_x_0_0_xx_x_z = buffer_1010_sdpp[2];

    auto g_x_0_x_0_0_xx_y_x = buffer_1010_sdpp[3];

    auto g_x_0_x_0_0_xx_y_y = buffer_1010_sdpp[4];

    auto g_x_0_x_0_0_xx_y_z = buffer_1010_sdpp[5];

    auto g_x_0_x_0_0_xx_z_x = buffer_1010_sdpp[6];

    auto g_x_0_x_0_0_xx_z_y = buffer_1010_sdpp[7];

    auto g_x_0_x_0_0_xx_z_z = buffer_1010_sdpp[8];

    auto g_x_0_x_0_0_xy_x_x = buffer_1010_sdpp[9];

    auto g_x_0_x_0_0_xy_x_y = buffer_1010_sdpp[10];

    auto g_x_0_x_0_0_xy_x_z = buffer_1010_sdpp[11];

    auto g_x_0_x_0_0_xy_y_x = buffer_1010_sdpp[12];

    auto g_x_0_x_0_0_xy_y_y = buffer_1010_sdpp[13];

    auto g_x_0_x_0_0_xy_y_z = buffer_1010_sdpp[14];

    auto g_x_0_x_0_0_xy_z_x = buffer_1010_sdpp[15];

    auto g_x_0_x_0_0_xy_z_y = buffer_1010_sdpp[16];

    auto g_x_0_x_0_0_xy_z_z = buffer_1010_sdpp[17];

    auto g_x_0_x_0_0_xz_x_x = buffer_1010_sdpp[18];

    auto g_x_0_x_0_0_xz_x_y = buffer_1010_sdpp[19];

    auto g_x_0_x_0_0_xz_x_z = buffer_1010_sdpp[20];

    auto g_x_0_x_0_0_xz_y_x = buffer_1010_sdpp[21];

    auto g_x_0_x_0_0_xz_y_y = buffer_1010_sdpp[22];

    auto g_x_0_x_0_0_xz_y_z = buffer_1010_sdpp[23];

    auto g_x_0_x_0_0_xz_z_x = buffer_1010_sdpp[24];

    auto g_x_0_x_0_0_xz_z_y = buffer_1010_sdpp[25];

    auto g_x_0_x_0_0_xz_z_z = buffer_1010_sdpp[26];

    auto g_x_0_x_0_0_yy_x_x = buffer_1010_sdpp[27];

    auto g_x_0_x_0_0_yy_x_y = buffer_1010_sdpp[28];

    auto g_x_0_x_0_0_yy_x_z = buffer_1010_sdpp[29];

    auto g_x_0_x_0_0_yy_y_x = buffer_1010_sdpp[30];

    auto g_x_0_x_0_0_yy_y_y = buffer_1010_sdpp[31];

    auto g_x_0_x_0_0_yy_y_z = buffer_1010_sdpp[32];

    auto g_x_0_x_0_0_yy_z_x = buffer_1010_sdpp[33];

    auto g_x_0_x_0_0_yy_z_y = buffer_1010_sdpp[34];

    auto g_x_0_x_0_0_yy_z_z = buffer_1010_sdpp[35];

    auto g_x_0_x_0_0_yz_x_x = buffer_1010_sdpp[36];

    auto g_x_0_x_0_0_yz_x_y = buffer_1010_sdpp[37];

    auto g_x_0_x_0_0_yz_x_z = buffer_1010_sdpp[38];

    auto g_x_0_x_0_0_yz_y_x = buffer_1010_sdpp[39];

    auto g_x_0_x_0_0_yz_y_y = buffer_1010_sdpp[40];

    auto g_x_0_x_0_0_yz_y_z = buffer_1010_sdpp[41];

    auto g_x_0_x_0_0_yz_z_x = buffer_1010_sdpp[42];

    auto g_x_0_x_0_0_yz_z_y = buffer_1010_sdpp[43];

    auto g_x_0_x_0_0_yz_z_z = buffer_1010_sdpp[44];

    auto g_x_0_x_0_0_zz_x_x = buffer_1010_sdpp[45];

    auto g_x_0_x_0_0_zz_x_y = buffer_1010_sdpp[46];

    auto g_x_0_x_0_0_zz_x_z = buffer_1010_sdpp[47];

    auto g_x_0_x_0_0_zz_y_x = buffer_1010_sdpp[48];

    auto g_x_0_x_0_0_zz_y_y = buffer_1010_sdpp[49];

    auto g_x_0_x_0_0_zz_y_z = buffer_1010_sdpp[50];

    auto g_x_0_x_0_0_zz_z_x = buffer_1010_sdpp[51];

    auto g_x_0_x_0_0_zz_z_y = buffer_1010_sdpp[52];

    auto g_x_0_x_0_0_zz_z_z = buffer_1010_sdpp[53];

    auto g_x_0_y_0_0_xx_x_x = buffer_1010_sdpp[54];

    auto g_x_0_y_0_0_xx_x_y = buffer_1010_sdpp[55];

    auto g_x_0_y_0_0_xx_x_z = buffer_1010_sdpp[56];

    auto g_x_0_y_0_0_xx_y_x = buffer_1010_sdpp[57];

    auto g_x_0_y_0_0_xx_y_y = buffer_1010_sdpp[58];

    auto g_x_0_y_0_0_xx_y_z = buffer_1010_sdpp[59];

    auto g_x_0_y_0_0_xx_z_x = buffer_1010_sdpp[60];

    auto g_x_0_y_0_0_xx_z_y = buffer_1010_sdpp[61];

    auto g_x_0_y_0_0_xx_z_z = buffer_1010_sdpp[62];

    auto g_x_0_y_0_0_xy_x_x = buffer_1010_sdpp[63];

    auto g_x_0_y_0_0_xy_x_y = buffer_1010_sdpp[64];

    auto g_x_0_y_0_0_xy_x_z = buffer_1010_sdpp[65];

    auto g_x_0_y_0_0_xy_y_x = buffer_1010_sdpp[66];

    auto g_x_0_y_0_0_xy_y_y = buffer_1010_sdpp[67];

    auto g_x_0_y_0_0_xy_y_z = buffer_1010_sdpp[68];

    auto g_x_0_y_0_0_xy_z_x = buffer_1010_sdpp[69];

    auto g_x_0_y_0_0_xy_z_y = buffer_1010_sdpp[70];

    auto g_x_0_y_0_0_xy_z_z = buffer_1010_sdpp[71];

    auto g_x_0_y_0_0_xz_x_x = buffer_1010_sdpp[72];

    auto g_x_0_y_0_0_xz_x_y = buffer_1010_sdpp[73];

    auto g_x_0_y_0_0_xz_x_z = buffer_1010_sdpp[74];

    auto g_x_0_y_0_0_xz_y_x = buffer_1010_sdpp[75];

    auto g_x_0_y_0_0_xz_y_y = buffer_1010_sdpp[76];

    auto g_x_0_y_0_0_xz_y_z = buffer_1010_sdpp[77];

    auto g_x_0_y_0_0_xz_z_x = buffer_1010_sdpp[78];

    auto g_x_0_y_0_0_xz_z_y = buffer_1010_sdpp[79];

    auto g_x_0_y_0_0_xz_z_z = buffer_1010_sdpp[80];

    auto g_x_0_y_0_0_yy_x_x = buffer_1010_sdpp[81];

    auto g_x_0_y_0_0_yy_x_y = buffer_1010_sdpp[82];

    auto g_x_0_y_0_0_yy_x_z = buffer_1010_sdpp[83];

    auto g_x_0_y_0_0_yy_y_x = buffer_1010_sdpp[84];

    auto g_x_0_y_0_0_yy_y_y = buffer_1010_sdpp[85];

    auto g_x_0_y_0_0_yy_y_z = buffer_1010_sdpp[86];

    auto g_x_0_y_0_0_yy_z_x = buffer_1010_sdpp[87];

    auto g_x_0_y_0_0_yy_z_y = buffer_1010_sdpp[88];

    auto g_x_0_y_0_0_yy_z_z = buffer_1010_sdpp[89];

    auto g_x_0_y_0_0_yz_x_x = buffer_1010_sdpp[90];

    auto g_x_0_y_0_0_yz_x_y = buffer_1010_sdpp[91];

    auto g_x_0_y_0_0_yz_x_z = buffer_1010_sdpp[92];

    auto g_x_0_y_0_0_yz_y_x = buffer_1010_sdpp[93];

    auto g_x_0_y_0_0_yz_y_y = buffer_1010_sdpp[94];

    auto g_x_0_y_0_0_yz_y_z = buffer_1010_sdpp[95];

    auto g_x_0_y_0_0_yz_z_x = buffer_1010_sdpp[96];

    auto g_x_0_y_0_0_yz_z_y = buffer_1010_sdpp[97];

    auto g_x_0_y_0_0_yz_z_z = buffer_1010_sdpp[98];

    auto g_x_0_y_0_0_zz_x_x = buffer_1010_sdpp[99];

    auto g_x_0_y_0_0_zz_x_y = buffer_1010_sdpp[100];

    auto g_x_0_y_0_0_zz_x_z = buffer_1010_sdpp[101];

    auto g_x_0_y_0_0_zz_y_x = buffer_1010_sdpp[102];

    auto g_x_0_y_0_0_zz_y_y = buffer_1010_sdpp[103];

    auto g_x_0_y_0_0_zz_y_z = buffer_1010_sdpp[104];

    auto g_x_0_y_0_0_zz_z_x = buffer_1010_sdpp[105];

    auto g_x_0_y_0_0_zz_z_y = buffer_1010_sdpp[106];

    auto g_x_0_y_0_0_zz_z_z = buffer_1010_sdpp[107];

    auto g_x_0_z_0_0_xx_x_x = buffer_1010_sdpp[108];

    auto g_x_0_z_0_0_xx_x_y = buffer_1010_sdpp[109];

    auto g_x_0_z_0_0_xx_x_z = buffer_1010_sdpp[110];

    auto g_x_0_z_0_0_xx_y_x = buffer_1010_sdpp[111];

    auto g_x_0_z_0_0_xx_y_y = buffer_1010_sdpp[112];

    auto g_x_0_z_0_0_xx_y_z = buffer_1010_sdpp[113];

    auto g_x_0_z_0_0_xx_z_x = buffer_1010_sdpp[114];

    auto g_x_0_z_0_0_xx_z_y = buffer_1010_sdpp[115];

    auto g_x_0_z_0_0_xx_z_z = buffer_1010_sdpp[116];

    auto g_x_0_z_0_0_xy_x_x = buffer_1010_sdpp[117];

    auto g_x_0_z_0_0_xy_x_y = buffer_1010_sdpp[118];

    auto g_x_0_z_0_0_xy_x_z = buffer_1010_sdpp[119];

    auto g_x_0_z_0_0_xy_y_x = buffer_1010_sdpp[120];

    auto g_x_0_z_0_0_xy_y_y = buffer_1010_sdpp[121];

    auto g_x_0_z_0_0_xy_y_z = buffer_1010_sdpp[122];

    auto g_x_0_z_0_0_xy_z_x = buffer_1010_sdpp[123];

    auto g_x_0_z_0_0_xy_z_y = buffer_1010_sdpp[124];

    auto g_x_0_z_0_0_xy_z_z = buffer_1010_sdpp[125];

    auto g_x_0_z_0_0_xz_x_x = buffer_1010_sdpp[126];

    auto g_x_0_z_0_0_xz_x_y = buffer_1010_sdpp[127];

    auto g_x_0_z_0_0_xz_x_z = buffer_1010_sdpp[128];

    auto g_x_0_z_0_0_xz_y_x = buffer_1010_sdpp[129];

    auto g_x_0_z_0_0_xz_y_y = buffer_1010_sdpp[130];

    auto g_x_0_z_0_0_xz_y_z = buffer_1010_sdpp[131];

    auto g_x_0_z_0_0_xz_z_x = buffer_1010_sdpp[132];

    auto g_x_0_z_0_0_xz_z_y = buffer_1010_sdpp[133];

    auto g_x_0_z_0_0_xz_z_z = buffer_1010_sdpp[134];

    auto g_x_0_z_0_0_yy_x_x = buffer_1010_sdpp[135];

    auto g_x_0_z_0_0_yy_x_y = buffer_1010_sdpp[136];

    auto g_x_0_z_0_0_yy_x_z = buffer_1010_sdpp[137];

    auto g_x_0_z_0_0_yy_y_x = buffer_1010_sdpp[138];

    auto g_x_0_z_0_0_yy_y_y = buffer_1010_sdpp[139];

    auto g_x_0_z_0_0_yy_y_z = buffer_1010_sdpp[140];

    auto g_x_0_z_0_0_yy_z_x = buffer_1010_sdpp[141];

    auto g_x_0_z_0_0_yy_z_y = buffer_1010_sdpp[142];

    auto g_x_0_z_0_0_yy_z_z = buffer_1010_sdpp[143];

    auto g_x_0_z_0_0_yz_x_x = buffer_1010_sdpp[144];

    auto g_x_0_z_0_0_yz_x_y = buffer_1010_sdpp[145];

    auto g_x_0_z_0_0_yz_x_z = buffer_1010_sdpp[146];

    auto g_x_0_z_0_0_yz_y_x = buffer_1010_sdpp[147];

    auto g_x_0_z_0_0_yz_y_y = buffer_1010_sdpp[148];

    auto g_x_0_z_0_0_yz_y_z = buffer_1010_sdpp[149];

    auto g_x_0_z_0_0_yz_z_x = buffer_1010_sdpp[150];

    auto g_x_0_z_0_0_yz_z_y = buffer_1010_sdpp[151];

    auto g_x_0_z_0_0_yz_z_z = buffer_1010_sdpp[152];

    auto g_x_0_z_0_0_zz_x_x = buffer_1010_sdpp[153];

    auto g_x_0_z_0_0_zz_x_y = buffer_1010_sdpp[154];

    auto g_x_0_z_0_0_zz_x_z = buffer_1010_sdpp[155];

    auto g_x_0_z_0_0_zz_y_x = buffer_1010_sdpp[156];

    auto g_x_0_z_0_0_zz_y_y = buffer_1010_sdpp[157];

    auto g_x_0_z_0_0_zz_y_z = buffer_1010_sdpp[158];

    auto g_x_0_z_0_0_zz_z_x = buffer_1010_sdpp[159];

    auto g_x_0_z_0_0_zz_z_y = buffer_1010_sdpp[160];

    auto g_x_0_z_0_0_zz_z_z = buffer_1010_sdpp[161];

    auto g_y_0_x_0_0_xx_x_x = buffer_1010_sdpp[162];

    auto g_y_0_x_0_0_xx_x_y = buffer_1010_sdpp[163];

    auto g_y_0_x_0_0_xx_x_z = buffer_1010_sdpp[164];

    auto g_y_0_x_0_0_xx_y_x = buffer_1010_sdpp[165];

    auto g_y_0_x_0_0_xx_y_y = buffer_1010_sdpp[166];

    auto g_y_0_x_0_0_xx_y_z = buffer_1010_sdpp[167];

    auto g_y_0_x_0_0_xx_z_x = buffer_1010_sdpp[168];

    auto g_y_0_x_0_0_xx_z_y = buffer_1010_sdpp[169];

    auto g_y_0_x_0_0_xx_z_z = buffer_1010_sdpp[170];

    auto g_y_0_x_0_0_xy_x_x = buffer_1010_sdpp[171];

    auto g_y_0_x_0_0_xy_x_y = buffer_1010_sdpp[172];

    auto g_y_0_x_0_0_xy_x_z = buffer_1010_sdpp[173];

    auto g_y_0_x_0_0_xy_y_x = buffer_1010_sdpp[174];

    auto g_y_0_x_0_0_xy_y_y = buffer_1010_sdpp[175];

    auto g_y_0_x_0_0_xy_y_z = buffer_1010_sdpp[176];

    auto g_y_0_x_0_0_xy_z_x = buffer_1010_sdpp[177];

    auto g_y_0_x_0_0_xy_z_y = buffer_1010_sdpp[178];

    auto g_y_0_x_0_0_xy_z_z = buffer_1010_sdpp[179];

    auto g_y_0_x_0_0_xz_x_x = buffer_1010_sdpp[180];

    auto g_y_0_x_0_0_xz_x_y = buffer_1010_sdpp[181];

    auto g_y_0_x_0_0_xz_x_z = buffer_1010_sdpp[182];

    auto g_y_0_x_0_0_xz_y_x = buffer_1010_sdpp[183];

    auto g_y_0_x_0_0_xz_y_y = buffer_1010_sdpp[184];

    auto g_y_0_x_0_0_xz_y_z = buffer_1010_sdpp[185];

    auto g_y_0_x_0_0_xz_z_x = buffer_1010_sdpp[186];

    auto g_y_0_x_0_0_xz_z_y = buffer_1010_sdpp[187];

    auto g_y_0_x_0_0_xz_z_z = buffer_1010_sdpp[188];

    auto g_y_0_x_0_0_yy_x_x = buffer_1010_sdpp[189];

    auto g_y_0_x_0_0_yy_x_y = buffer_1010_sdpp[190];

    auto g_y_0_x_0_0_yy_x_z = buffer_1010_sdpp[191];

    auto g_y_0_x_0_0_yy_y_x = buffer_1010_sdpp[192];

    auto g_y_0_x_0_0_yy_y_y = buffer_1010_sdpp[193];

    auto g_y_0_x_0_0_yy_y_z = buffer_1010_sdpp[194];

    auto g_y_0_x_0_0_yy_z_x = buffer_1010_sdpp[195];

    auto g_y_0_x_0_0_yy_z_y = buffer_1010_sdpp[196];

    auto g_y_0_x_0_0_yy_z_z = buffer_1010_sdpp[197];

    auto g_y_0_x_0_0_yz_x_x = buffer_1010_sdpp[198];

    auto g_y_0_x_0_0_yz_x_y = buffer_1010_sdpp[199];

    auto g_y_0_x_0_0_yz_x_z = buffer_1010_sdpp[200];

    auto g_y_0_x_0_0_yz_y_x = buffer_1010_sdpp[201];

    auto g_y_0_x_0_0_yz_y_y = buffer_1010_sdpp[202];

    auto g_y_0_x_0_0_yz_y_z = buffer_1010_sdpp[203];

    auto g_y_0_x_0_0_yz_z_x = buffer_1010_sdpp[204];

    auto g_y_0_x_0_0_yz_z_y = buffer_1010_sdpp[205];

    auto g_y_0_x_0_0_yz_z_z = buffer_1010_sdpp[206];

    auto g_y_0_x_0_0_zz_x_x = buffer_1010_sdpp[207];

    auto g_y_0_x_0_0_zz_x_y = buffer_1010_sdpp[208];

    auto g_y_0_x_0_0_zz_x_z = buffer_1010_sdpp[209];

    auto g_y_0_x_0_0_zz_y_x = buffer_1010_sdpp[210];

    auto g_y_0_x_0_0_zz_y_y = buffer_1010_sdpp[211];

    auto g_y_0_x_0_0_zz_y_z = buffer_1010_sdpp[212];

    auto g_y_0_x_0_0_zz_z_x = buffer_1010_sdpp[213];

    auto g_y_0_x_0_0_zz_z_y = buffer_1010_sdpp[214];

    auto g_y_0_x_0_0_zz_z_z = buffer_1010_sdpp[215];

    auto g_y_0_y_0_0_xx_x_x = buffer_1010_sdpp[216];

    auto g_y_0_y_0_0_xx_x_y = buffer_1010_sdpp[217];

    auto g_y_0_y_0_0_xx_x_z = buffer_1010_sdpp[218];

    auto g_y_0_y_0_0_xx_y_x = buffer_1010_sdpp[219];

    auto g_y_0_y_0_0_xx_y_y = buffer_1010_sdpp[220];

    auto g_y_0_y_0_0_xx_y_z = buffer_1010_sdpp[221];

    auto g_y_0_y_0_0_xx_z_x = buffer_1010_sdpp[222];

    auto g_y_0_y_0_0_xx_z_y = buffer_1010_sdpp[223];

    auto g_y_0_y_0_0_xx_z_z = buffer_1010_sdpp[224];

    auto g_y_0_y_0_0_xy_x_x = buffer_1010_sdpp[225];

    auto g_y_0_y_0_0_xy_x_y = buffer_1010_sdpp[226];

    auto g_y_0_y_0_0_xy_x_z = buffer_1010_sdpp[227];

    auto g_y_0_y_0_0_xy_y_x = buffer_1010_sdpp[228];

    auto g_y_0_y_0_0_xy_y_y = buffer_1010_sdpp[229];

    auto g_y_0_y_0_0_xy_y_z = buffer_1010_sdpp[230];

    auto g_y_0_y_0_0_xy_z_x = buffer_1010_sdpp[231];

    auto g_y_0_y_0_0_xy_z_y = buffer_1010_sdpp[232];

    auto g_y_0_y_0_0_xy_z_z = buffer_1010_sdpp[233];

    auto g_y_0_y_0_0_xz_x_x = buffer_1010_sdpp[234];

    auto g_y_0_y_0_0_xz_x_y = buffer_1010_sdpp[235];

    auto g_y_0_y_0_0_xz_x_z = buffer_1010_sdpp[236];

    auto g_y_0_y_0_0_xz_y_x = buffer_1010_sdpp[237];

    auto g_y_0_y_0_0_xz_y_y = buffer_1010_sdpp[238];

    auto g_y_0_y_0_0_xz_y_z = buffer_1010_sdpp[239];

    auto g_y_0_y_0_0_xz_z_x = buffer_1010_sdpp[240];

    auto g_y_0_y_0_0_xz_z_y = buffer_1010_sdpp[241];

    auto g_y_0_y_0_0_xz_z_z = buffer_1010_sdpp[242];

    auto g_y_0_y_0_0_yy_x_x = buffer_1010_sdpp[243];

    auto g_y_0_y_0_0_yy_x_y = buffer_1010_sdpp[244];

    auto g_y_0_y_0_0_yy_x_z = buffer_1010_sdpp[245];

    auto g_y_0_y_0_0_yy_y_x = buffer_1010_sdpp[246];

    auto g_y_0_y_0_0_yy_y_y = buffer_1010_sdpp[247];

    auto g_y_0_y_0_0_yy_y_z = buffer_1010_sdpp[248];

    auto g_y_0_y_0_0_yy_z_x = buffer_1010_sdpp[249];

    auto g_y_0_y_0_0_yy_z_y = buffer_1010_sdpp[250];

    auto g_y_0_y_0_0_yy_z_z = buffer_1010_sdpp[251];

    auto g_y_0_y_0_0_yz_x_x = buffer_1010_sdpp[252];

    auto g_y_0_y_0_0_yz_x_y = buffer_1010_sdpp[253];

    auto g_y_0_y_0_0_yz_x_z = buffer_1010_sdpp[254];

    auto g_y_0_y_0_0_yz_y_x = buffer_1010_sdpp[255];

    auto g_y_0_y_0_0_yz_y_y = buffer_1010_sdpp[256];

    auto g_y_0_y_0_0_yz_y_z = buffer_1010_sdpp[257];

    auto g_y_0_y_0_0_yz_z_x = buffer_1010_sdpp[258];

    auto g_y_0_y_0_0_yz_z_y = buffer_1010_sdpp[259];

    auto g_y_0_y_0_0_yz_z_z = buffer_1010_sdpp[260];

    auto g_y_0_y_0_0_zz_x_x = buffer_1010_sdpp[261];

    auto g_y_0_y_0_0_zz_x_y = buffer_1010_sdpp[262];

    auto g_y_0_y_0_0_zz_x_z = buffer_1010_sdpp[263];

    auto g_y_0_y_0_0_zz_y_x = buffer_1010_sdpp[264];

    auto g_y_0_y_0_0_zz_y_y = buffer_1010_sdpp[265];

    auto g_y_0_y_0_0_zz_y_z = buffer_1010_sdpp[266];

    auto g_y_0_y_0_0_zz_z_x = buffer_1010_sdpp[267];

    auto g_y_0_y_0_0_zz_z_y = buffer_1010_sdpp[268];

    auto g_y_0_y_0_0_zz_z_z = buffer_1010_sdpp[269];

    auto g_y_0_z_0_0_xx_x_x = buffer_1010_sdpp[270];

    auto g_y_0_z_0_0_xx_x_y = buffer_1010_sdpp[271];

    auto g_y_0_z_0_0_xx_x_z = buffer_1010_sdpp[272];

    auto g_y_0_z_0_0_xx_y_x = buffer_1010_sdpp[273];

    auto g_y_0_z_0_0_xx_y_y = buffer_1010_sdpp[274];

    auto g_y_0_z_0_0_xx_y_z = buffer_1010_sdpp[275];

    auto g_y_0_z_0_0_xx_z_x = buffer_1010_sdpp[276];

    auto g_y_0_z_0_0_xx_z_y = buffer_1010_sdpp[277];

    auto g_y_0_z_0_0_xx_z_z = buffer_1010_sdpp[278];

    auto g_y_0_z_0_0_xy_x_x = buffer_1010_sdpp[279];

    auto g_y_0_z_0_0_xy_x_y = buffer_1010_sdpp[280];

    auto g_y_0_z_0_0_xy_x_z = buffer_1010_sdpp[281];

    auto g_y_0_z_0_0_xy_y_x = buffer_1010_sdpp[282];

    auto g_y_0_z_0_0_xy_y_y = buffer_1010_sdpp[283];

    auto g_y_0_z_0_0_xy_y_z = buffer_1010_sdpp[284];

    auto g_y_0_z_0_0_xy_z_x = buffer_1010_sdpp[285];

    auto g_y_0_z_0_0_xy_z_y = buffer_1010_sdpp[286];

    auto g_y_0_z_0_0_xy_z_z = buffer_1010_sdpp[287];

    auto g_y_0_z_0_0_xz_x_x = buffer_1010_sdpp[288];

    auto g_y_0_z_0_0_xz_x_y = buffer_1010_sdpp[289];

    auto g_y_0_z_0_0_xz_x_z = buffer_1010_sdpp[290];

    auto g_y_0_z_0_0_xz_y_x = buffer_1010_sdpp[291];

    auto g_y_0_z_0_0_xz_y_y = buffer_1010_sdpp[292];

    auto g_y_0_z_0_0_xz_y_z = buffer_1010_sdpp[293];

    auto g_y_0_z_0_0_xz_z_x = buffer_1010_sdpp[294];

    auto g_y_0_z_0_0_xz_z_y = buffer_1010_sdpp[295];

    auto g_y_0_z_0_0_xz_z_z = buffer_1010_sdpp[296];

    auto g_y_0_z_0_0_yy_x_x = buffer_1010_sdpp[297];

    auto g_y_0_z_0_0_yy_x_y = buffer_1010_sdpp[298];

    auto g_y_0_z_0_0_yy_x_z = buffer_1010_sdpp[299];

    auto g_y_0_z_0_0_yy_y_x = buffer_1010_sdpp[300];

    auto g_y_0_z_0_0_yy_y_y = buffer_1010_sdpp[301];

    auto g_y_0_z_0_0_yy_y_z = buffer_1010_sdpp[302];

    auto g_y_0_z_0_0_yy_z_x = buffer_1010_sdpp[303];

    auto g_y_0_z_0_0_yy_z_y = buffer_1010_sdpp[304];

    auto g_y_0_z_0_0_yy_z_z = buffer_1010_sdpp[305];

    auto g_y_0_z_0_0_yz_x_x = buffer_1010_sdpp[306];

    auto g_y_0_z_0_0_yz_x_y = buffer_1010_sdpp[307];

    auto g_y_0_z_0_0_yz_x_z = buffer_1010_sdpp[308];

    auto g_y_0_z_0_0_yz_y_x = buffer_1010_sdpp[309];

    auto g_y_0_z_0_0_yz_y_y = buffer_1010_sdpp[310];

    auto g_y_0_z_0_0_yz_y_z = buffer_1010_sdpp[311];

    auto g_y_0_z_0_0_yz_z_x = buffer_1010_sdpp[312];

    auto g_y_0_z_0_0_yz_z_y = buffer_1010_sdpp[313];

    auto g_y_0_z_0_0_yz_z_z = buffer_1010_sdpp[314];

    auto g_y_0_z_0_0_zz_x_x = buffer_1010_sdpp[315];

    auto g_y_0_z_0_0_zz_x_y = buffer_1010_sdpp[316];

    auto g_y_0_z_0_0_zz_x_z = buffer_1010_sdpp[317];

    auto g_y_0_z_0_0_zz_y_x = buffer_1010_sdpp[318];

    auto g_y_0_z_0_0_zz_y_y = buffer_1010_sdpp[319];

    auto g_y_0_z_0_0_zz_y_z = buffer_1010_sdpp[320];

    auto g_y_0_z_0_0_zz_z_x = buffer_1010_sdpp[321];

    auto g_y_0_z_0_0_zz_z_y = buffer_1010_sdpp[322];

    auto g_y_0_z_0_0_zz_z_z = buffer_1010_sdpp[323];

    auto g_z_0_x_0_0_xx_x_x = buffer_1010_sdpp[324];

    auto g_z_0_x_0_0_xx_x_y = buffer_1010_sdpp[325];

    auto g_z_0_x_0_0_xx_x_z = buffer_1010_sdpp[326];

    auto g_z_0_x_0_0_xx_y_x = buffer_1010_sdpp[327];

    auto g_z_0_x_0_0_xx_y_y = buffer_1010_sdpp[328];

    auto g_z_0_x_0_0_xx_y_z = buffer_1010_sdpp[329];

    auto g_z_0_x_0_0_xx_z_x = buffer_1010_sdpp[330];

    auto g_z_0_x_0_0_xx_z_y = buffer_1010_sdpp[331];

    auto g_z_0_x_0_0_xx_z_z = buffer_1010_sdpp[332];

    auto g_z_0_x_0_0_xy_x_x = buffer_1010_sdpp[333];

    auto g_z_0_x_0_0_xy_x_y = buffer_1010_sdpp[334];

    auto g_z_0_x_0_0_xy_x_z = buffer_1010_sdpp[335];

    auto g_z_0_x_0_0_xy_y_x = buffer_1010_sdpp[336];

    auto g_z_0_x_0_0_xy_y_y = buffer_1010_sdpp[337];

    auto g_z_0_x_0_0_xy_y_z = buffer_1010_sdpp[338];

    auto g_z_0_x_0_0_xy_z_x = buffer_1010_sdpp[339];

    auto g_z_0_x_0_0_xy_z_y = buffer_1010_sdpp[340];

    auto g_z_0_x_0_0_xy_z_z = buffer_1010_sdpp[341];

    auto g_z_0_x_0_0_xz_x_x = buffer_1010_sdpp[342];

    auto g_z_0_x_0_0_xz_x_y = buffer_1010_sdpp[343];

    auto g_z_0_x_0_0_xz_x_z = buffer_1010_sdpp[344];

    auto g_z_0_x_0_0_xz_y_x = buffer_1010_sdpp[345];

    auto g_z_0_x_0_0_xz_y_y = buffer_1010_sdpp[346];

    auto g_z_0_x_0_0_xz_y_z = buffer_1010_sdpp[347];

    auto g_z_0_x_0_0_xz_z_x = buffer_1010_sdpp[348];

    auto g_z_0_x_0_0_xz_z_y = buffer_1010_sdpp[349];

    auto g_z_0_x_0_0_xz_z_z = buffer_1010_sdpp[350];

    auto g_z_0_x_0_0_yy_x_x = buffer_1010_sdpp[351];

    auto g_z_0_x_0_0_yy_x_y = buffer_1010_sdpp[352];

    auto g_z_0_x_0_0_yy_x_z = buffer_1010_sdpp[353];

    auto g_z_0_x_0_0_yy_y_x = buffer_1010_sdpp[354];

    auto g_z_0_x_0_0_yy_y_y = buffer_1010_sdpp[355];

    auto g_z_0_x_0_0_yy_y_z = buffer_1010_sdpp[356];

    auto g_z_0_x_0_0_yy_z_x = buffer_1010_sdpp[357];

    auto g_z_0_x_0_0_yy_z_y = buffer_1010_sdpp[358];

    auto g_z_0_x_0_0_yy_z_z = buffer_1010_sdpp[359];

    auto g_z_0_x_0_0_yz_x_x = buffer_1010_sdpp[360];

    auto g_z_0_x_0_0_yz_x_y = buffer_1010_sdpp[361];

    auto g_z_0_x_0_0_yz_x_z = buffer_1010_sdpp[362];

    auto g_z_0_x_0_0_yz_y_x = buffer_1010_sdpp[363];

    auto g_z_0_x_0_0_yz_y_y = buffer_1010_sdpp[364];

    auto g_z_0_x_0_0_yz_y_z = buffer_1010_sdpp[365];

    auto g_z_0_x_0_0_yz_z_x = buffer_1010_sdpp[366];

    auto g_z_0_x_0_0_yz_z_y = buffer_1010_sdpp[367];

    auto g_z_0_x_0_0_yz_z_z = buffer_1010_sdpp[368];

    auto g_z_0_x_0_0_zz_x_x = buffer_1010_sdpp[369];

    auto g_z_0_x_0_0_zz_x_y = buffer_1010_sdpp[370];

    auto g_z_0_x_0_0_zz_x_z = buffer_1010_sdpp[371];

    auto g_z_0_x_0_0_zz_y_x = buffer_1010_sdpp[372];

    auto g_z_0_x_0_0_zz_y_y = buffer_1010_sdpp[373];

    auto g_z_0_x_0_0_zz_y_z = buffer_1010_sdpp[374];

    auto g_z_0_x_0_0_zz_z_x = buffer_1010_sdpp[375];

    auto g_z_0_x_0_0_zz_z_y = buffer_1010_sdpp[376];

    auto g_z_0_x_0_0_zz_z_z = buffer_1010_sdpp[377];

    auto g_z_0_y_0_0_xx_x_x = buffer_1010_sdpp[378];

    auto g_z_0_y_0_0_xx_x_y = buffer_1010_sdpp[379];

    auto g_z_0_y_0_0_xx_x_z = buffer_1010_sdpp[380];

    auto g_z_0_y_0_0_xx_y_x = buffer_1010_sdpp[381];

    auto g_z_0_y_0_0_xx_y_y = buffer_1010_sdpp[382];

    auto g_z_0_y_0_0_xx_y_z = buffer_1010_sdpp[383];

    auto g_z_0_y_0_0_xx_z_x = buffer_1010_sdpp[384];

    auto g_z_0_y_0_0_xx_z_y = buffer_1010_sdpp[385];

    auto g_z_0_y_0_0_xx_z_z = buffer_1010_sdpp[386];

    auto g_z_0_y_0_0_xy_x_x = buffer_1010_sdpp[387];

    auto g_z_0_y_0_0_xy_x_y = buffer_1010_sdpp[388];

    auto g_z_0_y_0_0_xy_x_z = buffer_1010_sdpp[389];

    auto g_z_0_y_0_0_xy_y_x = buffer_1010_sdpp[390];

    auto g_z_0_y_0_0_xy_y_y = buffer_1010_sdpp[391];

    auto g_z_0_y_0_0_xy_y_z = buffer_1010_sdpp[392];

    auto g_z_0_y_0_0_xy_z_x = buffer_1010_sdpp[393];

    auto g_z_0_y_0_0_xy_z_y = buffer_1010_sdpp[394];

    auto g_z_0_y_0_0_xy_z_z = buffer_1010_sdpp[395];

    auto g_z_0_y_0_0_xz_x_x = buffer_1010_sdpp[396];

    auto g_z_0_y_0_0_xz_x_y = buffer_1010_sdpp[397];

    auto g_z_0_y_0_0_xz_x_z = buffer_1010_sdpp[398];

    auto g_z_0_y_0_0_xz_y_x = buffer_1010_sdpp[399];

    auto g_z_0_y_0_0_xz_y_y = buffer_1010_sdpp[400];

    auto g_z_0_y_0_0_xz_y_z = buffer_1010_sdpp[401];

    auto g_z_0_y_0_0_xz_z_x = buffer_1010_sdpp[402];

    auto g_z_0_y_0_0_xz_z_y = buffer_1010_sdpp[403];

    auto g_z_0_y_0_0_xz_z_z = buffer_1010_sdpp[404];

    auto g_z_0_y_0_0_yy_x_x = buffer_1010_sdpp[405];

    auto g_z_0_y_0_0_yy_x_y = buffer_1010_sdpp[406];

    auto g_z_0_y_0_0_yy_x_z = buffer_1010_sdpp[407];

    auto g_z_0_y_0_0_yy_y_x = buffer_1010_sdpp[408];

    auto g_z_0_y_0_0_yy_y_y = buffer_1010_sdpp[409];

    auto g_z_0_y_0_0_yy_y_z = buffer_1010_sdpp[410];

    auto g_z_0_y_0_0_yy_z_x = buffer_1010_sdpp[411];

    auto g_z_0_y_0_0_yy_z_y = buffer_1010_sdpp[412];

    auto g_z_0_y_0_0_yy_z_z = buffer_1010_sdpp[413];

    auto g_z_0_y_0_0_yz_x_x = buffer_1010_sdpp[414];

    auto g_z_0_y_0_0_yz_x_y = buffer_1010_sdpp[415];

    auto g_z_0_y_0_0_yz_x_z = buffer_1010_sdpp[416];

    auto g_z_0_y_0_0_yz_y_x = buffer_1010_sdpp[417];

    auto g_z_0_y_0_0_yz_y_y = buffer_1010_sdpp[418];

    auto g_z_0_y_0_0_yz_y_z = buffer_1010_sdpp[419];

    auto g_z_0_y_0_0_yz_z_x = buffer_1010_sdpp[420];

    auto g_z_0_y_0_0_yz_z_y = buffer_1010_sdpp[421];

    auto g_z_0_y_0_0_yz_z_z = buffer_1010_sdpp[422];

    auto g_z_0_y_0_0_zz_x_x = buffer_1010_sdpp[423];

    auto g_z_0_y_0_0_zz_x_y = buffer_1010_sdpp[424];

    auto g_z_0_y_0_0_zz_x_z = buffer_1010_sdpp[425];

    auto g_z_0_y_0_0_zz_y_x = buffer_1010_sdpp[426];

    auto g_z_0_y_0_0_zz_y_y = buffer_1010_sdpp[427];

    auto g_z_0_y_0_0_zz_y_z = buffer_1010_sdpp[428];

    auto g_z_0_y_0_0_zz_z_x = buffer_1010_sdpp[429];

    auto g_z_0_y_0_0_zz_z_y = buffer_1010_sdpp[430];

    auto g_z_0_y_0_0_zz_z_z = buffer_1010_sdpp[431];

    auto g_z_0_z_0_0_xx_x_x = buffer_1010_sdpp[432];

    auto g_z_0_z_0_0_xx_x_y = buffer_1010_sdpp[433];

    auto g_z_0_z_0_0_xx_x_z = buffer_1010_sdpp[434];

    auto g_z_0_z_0_0_xx_y_x = buffer_1010_sdpp[435];

    auto g_z_0_z_0_0_xx_y_y = buffer_1010_sdpp[436];

    auto g_z_0_z_0_0_xx_y_z = buffer_1010_sdpp[437];

    auto g_z_0_z_0_0_xx_z_x = buffer_1010_sdpp[438];

    auto g_z_0_z_0_0_xx_z_y = buffer_1010_sdpp[439];

    auto g_z_0_z_0_0_xx_z_z = buffer_1010_sdpp[440];

    auto g_z_0_z_0_0_xy_x_x = buffer_1010_sdpp[441];

    auto g_z_0_z_0_0_xy_x_y = buffer_1010_sdpp[442];

    auto g_z_0_z_0_0_xy_x_z = buffer_1010_sdpp[443];

    auto g_z_0_z_0_0_xy_y_x = buffer_1010_sdpp[444];

    auto g_z_0_z_0_0_xy_y_y = buffer_1010_sdpp[445];

    auto g_z_0_z_0_0_xy_y_z = buffer_1010_sdpp[446];

    auto g_z_0_z_0_0_xy_z_x = buffer_1010_sdpp[447];

    auto g_z_0_z_0_0_xy_z_y = buffer_1010_sdpp[448];

    auto g_z_0_z_0_0_xy_z_z = buffer_1010_sdpp[449];

    auto g_z_0_z_0_0_xz_x_x = buffer_1010_sdpp[450];

    auto g_z_0_z_0_0_xz_x_y = buffer_1010_sdpp[451];

    auto g_z_0_z_0_0_xz_x_z = buffer_1010_sdpp[452];

    auto g_z_0_z_0_0_xz_y_x = buffer_1010_sdpp[453];

    auto g_z_0_z_0_0_xz_y_y = buffer_1010_sdpp[454];

    auto g_z_0_z_0_0_xz_y_z = buffer_1010_sdpp[455];

    auto g_z_0_z_0_0_xz_z_x = buffer_1010_sdpp[456];

    auto g_z_0_z_0_0_xz_z_y = buffer_1010_sdpp[457];

    auto g_z_0_z_0_0_xz_z_z = buffer_1010_sdpp[458];

    auto g_z_0_z_0_0_yy_x_x = buffer_1010_sdpp[459];

    auto g_z_0_z_0_0_yy_x_y = buffer_1010_sdpp[460];

    auto g_z_0_z_0_0_yy_x_z = buffer_1010_sdpp[461];

    auto g_z_0_z_0_0_yy_y_x = buffer_1010_sdpp[462];

    auto g_z_0_z_0_0_yy_y_y = buffer_1010_sdpp[463];

    auto g_z_0_z_0_0_yy_y_z = buffer_1010_sdpp[464];

    auto g_z_0_z_0_0_yy_z_x = buffer_1010_sdpp[465];

    auto g_z_0_z_0_0_yy_z_y = buffer_1010_sdpp[466];

    auto g_z_0_z_0_0_yy_z_z = buffer_1010_sdpp[467];

    auto g_z_0_z_0_0_yz_x_x = buffer_1010_sdpp[468];

    auto g_z_0_z_0_0_yz_x_y = buffer_1010_sdpp[469];

    auto g_z_0_z_0_0_yz_x_z = buffer_1010_sdpp[470];

    auto g_z_0_z_0_0_yz_y_x = buffer_1010_sdpp[471];

    auto g_z_0_z_0_0_yz_y_y = buffer_1010_sdpp[472];

    auto g_z_0_z_0_0_yz_y_z = buffer_1010_sdpp[473];

    auto g_z_0_z_0_0_yz_z_x = buffer_1010_sdpp[474];

    auto g_z_0_z_0_0_yz_z_y = buffer_1010_sdpp[475];

    auto g_z_0_z_0_0_yz_z_z = buffer_1010_sdpp[476];

    auto g_z_0_z_0_0_zz_x_x = buffer_1010_sdpp[477];

    auto g_z_0_z_0_0_zz_x_y = buffer_1010_sdpp[478];

    auto g_z_0_z_0_0_zz_x_z = buffer_1010_sdpp[479];

    auto g_z_0_z_0_0_zz_y_x = buffer_1010_sdpp[480];

    auto g_z_0_z_0_0_zz_y_y = buffer_1010_sdpp[481];

    auto g_z_0_z_0_0_zz_y_z = buffer_1010_sdpp[482];

    auto g_z_0_z_0_0_zz_z_x = buffer_1010_sdpp[483];

    auto g_z_0_z_0_0_zz_z_y = buffer_1010_sdpp[484];

    auto g_z_0_z_0_0_zz_z_z = buffer_1010_sdpp[485];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_x_0_0_xx_x_x, g_x_0_x_0_0_xx_x_y, g_x_0_x_0_0_xx_x_z, g_x_xx_0_x, g_x_xx_0_y, g_x_xx_0_z, g_x_xx_xx_x, g_x_xx_xx_y, g_x_xx_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xx_x_x[i] = -2.0 * g_x_xx_0_x[i] * a_exp + 4.0 * g_x_xx_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_x_y[i] = -2.0 * g_x_xx_0_y[i] * a_exp + 4.0 * g_x_xx_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_x_z[i] = -2.0 * g_x_xx_0_z[i] * a_exp + 4.0 * g_x_xx_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_x_0_0_xx_y_x, g_x_0_x_0_0_xx_y_y, g_x_0_x_0_0_xx_y_z, g_x_xx_xy_x, g_x_xx_xy_y, g_x_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xx_y_x[i] = 4.0 * g_x_xx_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_y_y[i] = 4.0 * g_x_xx_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_y_z[i] = 4.0 * g_x_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_x_0_0_xx_z_x, g_x_0_x_0_0_xx_z_y, g_x_0_x_0_0_xx_z_z, g_x_xx_xz_x, g_x_xx_xz_y, g_x_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xx_z_x[i] = 4.0 * g_x_xx_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_z_y[i] = 4.0 * g_x_xx_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_z_z[i] = 4.0 * g_x_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_0_x_0_0_xy_x_x, g_x_0_x_0_0_xy_x_y, g_x_0_x_0_0_xy_x_z, g_x_xy_0_x, g_x_xy_0_y, g_x_xy_0_z, g_x_xy_xx_x, g_x_xy_xx_y, g_x_xy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xy_x_x[i] = -2.0 * g_x_xy_0_x[i] * a_exp + 4.0 * g_x_xy_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_x_y[i] = -2.0 * g_x_xy_0_y[i] * a_exp + 4.0 * g_x_xy_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_x_z[i] = -2.0 * g_x_xy_0_z[i] * a_exp + 4.0 * g_x_xy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_0_x_0_0_xy_y_x, g_x_0_x_0_0_xy_y_y, g_x_0_x_0_0_xy_y_z, g_x_xy_xy_x, g_x_xy_xy_y, g_x_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xy_y_x[i] = 4.0 * g_x_xy_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_y_y[i] = 4.0 * g_x_xy_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_y_z[i] = 4.0 * g_x_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_0_x_0_0_xy_z_x, g_x_0_x_0_0_xy_z_y, g_x_0_x_0_0_xy_z_z, g_x_xy_xz_x, g_x_xy_xz_y, g_x_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xy_z_x[i] = 4.0 * g_x_xy_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_z_y[i] = 4.0 * g_x_xy_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_z_z[i] = 4.0 * g_x_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_0_x_0_0_xz_x_x, g_x_0_x_0_0_xz_x_y, g_x_0_x_0_0_xz_x_z, g_x_xz_0_x, g_x_xz_0_y, g_x_xz_0_z, g_x_xz_xx_x, g_x_xz_xx_y, g_x_xz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xz_x_x[i] = -2.0 * g_x_xz_0_x[i] * a_exp + 4.0 * g_x_xz_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_x_y[i] = -2.0 * g_x_xz_0_y[i] * a_exp + 4.0 * g_x_xz_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_x_z[i] = -2.0 * g_x_xz_0_z[i] * a_exp + 4.0 * g_x_xz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_0_x_0_0_xz_y_x, g_x_0_x_0_0_xz_y_y, g_x_0_x_0_0_xz_y_z, g_x_xz_xy_x, g_x_xz_xy_y, g_x_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xz_y_x[i] = 4.0 * g_x_xz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_y_y[i] = 4.0 * g_x_xz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_y_z[i] = 4.0 * g_x_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_0_x_0_0_xz_z_x, g_x_0_x_0_0_xz_z_y, g_x_0_x_0_0_xz_z_z, g_x_xz_xz_x, g_x_xz_xz_y, g_x_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xz_z_x[i] = 4.0 * g_x_xz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_z_y[i] = 4.0 * g_x_xz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_z_z[i] = 4.0 * g_x_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_0_x_0_0_yy_x_x, g_x_0_x_0_0_yy_x_y, g_x_0_x_0_0_yy_x_z, g_x_yy_0_x, g_x_yy_0_y, g_x_yy_0_z, g_x_yy_xx_x, g_x_yy_xx_y, g_x_yy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yy_x_x[i] = -2.0 * g_x_yy_0_x[i] * a_exp + 4.0 * g_x_yy_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_x_y[i] = -2.0 * g_x_yy_0_y[i] * a_exp + 4.0 * g_x_yy_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_x_z[i] = -2.0 * g_x_yy_0_z[i] * a_exp + 4.0 * g_x_yy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_x_0_x_0_0_yy_y_x, g_x_0_x_0_0_yy_y_y, g_x_0_x_0_0_yy_y_z, g_x_yy_xy_x, g_x_yy_xy_y, g_x_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yy_y_x[i] = 4.0 * g_x_yy_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_y_y[i] = 4.0 * g_x_yy_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_y_z[i] = 4.0 * g_x_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_x_0_x_0_0_yy_z_x, g_x_0_x_0_0_yy_z_y, g_x_0_x_0_0_yy_z_z, g_x_yy_xz_x, g_x_yy_xz_y, g_x_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yy_z_x[i] = 4.0 * g_x_yy_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_z_y[i] = 4.0 * g_x_yy_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_z_z[i] = 4.0 * g_x_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_0_x_0_0_yz_x_x, g_x_0_x_0_0_yz_x_y, g_x_0_x_0_0_yz_x_z, g_x_yz_0_x, g_x_yz_0_y, g_x_yz_0_z, g_x_yz_xx_x, g_x_yz_xx_y, g_x_yz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yz_x_x[i] = -2.0 * g_x_yz_0_x[i] * a_exp + 4.0 * g_x_yz_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_x_y[i] = -2.0 * g_x_yz_0_y[i] * a_exp + 4.0 * g_x_yz_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_x_z[i] = -2.0 * g_x_yz_0_z[i] * a_exp + 4.0 * g_x_yz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_0_x_0_0_yz_y_x, g_x_0_x_0_0_yz_y_y, g_x_0_x_0_0_yz_y_z, g_x_yz_xy_x, g_x_yz_xy_y, g_x_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yz_y_x[i] = 4.0 * g_x_yz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_y_y[i] = 4.0 * g_x_yz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_y_z[i] = 4.0 * g_x_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_0_x_0_0_yz_z_x, g_x_0_x_0_0_yz_z_y, g_x_0_x_0_0_yz_z_z, g_x_yz_xz_x, g_x_yz_xz_y, g_x_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yz_z_x[i] = 4.0 * g_x_yz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_z_y[i] = 4.0 * g_x_yz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_z_z[i] = 4.0 * g_x_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_0_x_0_0_zz_x_x, g_x_0_x_0_0_zz_x_y, g_x_0_x_0_0_zz_x_z, g_x_zz_0_x, g_x_zz_0_y, g_x_zz_0_z, g_x_zz_xx_x, g_x_zz_xx_y, g_x_zz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_zz_x_x[i] = -2.0 * g_x_zz_0_x[i] * a_exp + 4.0 * g_x_zz_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_x_y[i] = -2.0 * g_x_zz_0_y[i] * a_exp + 4.0 * g_x_zz_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_x_z[i] = -2.0 * g_x_zz_0_z[i] * a_exp + 4.0 * g_x_zz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_0_x_0_0_zz_y_x, g_x_0_x_0_0_zz_y_y, g_x_0_x_0_0_zz_y_z, g_x_zz_xy_x, g_x_zz_xy_y, g_x_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_zz_y_x[i] = 4.0 * g_x_zz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_y_y[i] = 4.0 * g_x_zz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_y_z[i] = 4.0 * g_x_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_0_x_0_0_zz_z_x, g_x_0_x_0_0_zz_z_y, g_x_0_x_0_0_zz_z_z, g_x_zz_xz_x, g_x_zz_xz_y, g_x_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_zz_z_x[i] = 4.0 * g_x_zz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_z_y[i] = 4.0 * g_x_zz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_z_z[i] = 4.0 * g_x_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_x_0_y_0_0_xx_x_x, g_x_0_y_0_0_xx_x_y, g_x_0_y_0_0_xx_x_z, g_x_xx_xy_x, g_x_xx_xy_y, g_x_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xx_x_x[i] = 4.0 * g_x_xx_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_x_y[i] = 4.0 * g_x_xx_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_x_z[i] = 4.0 * g_x_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_x_0_y_0_0_xx_y_x, g_x_0_y_0_0_xx_y_y, g_x_0_y_0_0_xx_y_z, g_x_xx_0_x, g_x_xx_0_y, g_x_xx_0_z, g_x_xx_yy_x, g_x_xx_yy_y, g_x_xx_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xx_y_x[i] = -2.0 * g_x_xx_0_x[i] * a_exp + 4.0 * g_x_xx_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_y_y[i] = -2.0 * g_x_xx_0_y[i] * a_exp + 4.0 * g_x_xx_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_y_z[i] = -2.0 * g_x_xx_0_z[i] * a_exp + 4.0 * g_x_xx_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_x_0_y_0_0_xx_z_x, g_x_0_y_0_0_xx_z_y, g_x_0_y_0_0_xx_z_z, g_x_xx_yz_x, g_x_xx_yz_y, g_x_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xx_z_x[i] = 4.0 * g_x_xx_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_z_y[i] = 4.0 * g_x_xx_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_z_z[i] = 4.0 * g_x_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_x_0_y_0_0_xy_x_x, g_x_0_y_0_0_xy_x_y, g_x_0_y_0_0_xy_x_z, g_x_xy_xy_x, g_x_xy_xy_y, g_x_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xy_x_x[i] = 4.0 * g_x_xy_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_x_y[i] = 4.0 * g_x_xy_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_x_z[i] = 4.0 * g_x_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_x_0_y_0_0_xy_y_x, g_x_0_y_0_0_xy_y_y, g_x_0_y_0_0_xy_y_z, g_x_xy_0_x, g_x_xy_0_y, g_x_xy_0_z, g_x_xy_yy_x, g_x_xy_yy_y, g_x_xy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xy_y_x[i] = -2.0 * g_x_xy_0_x[i] * a_exp + 4.0 * g_x_xy_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_y_y[i] = -2.0 * g_x_xy_0_y[i] * a_exp + 4.0 * g_x_xy_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_y_z[i] = -2.0 * g_x_xy_0_z[i] * a_exp + 4.0 * g_x_xy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_x_0_y_0_0_xy_z_x, g_x_0_y_0_0_xy_z_y, g_x_0_y_0_0_xy_z_z, g_x_xy_yz_x, g_x_xy_yz_y, g_x_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xy_z_x[i] = 4.0 * g_x_xy_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_z_y[i] = 4.0 * g_x_xy_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_z_z[i] = 4.0 * g_x_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_0_y_0_0_xz_x_x, g_x_0_y_0_0_xz_x_y, g_x_0_y_0_0_xz_x_z, g_x_xz_xy_x, g_x_xz_xy_y, g_x_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xz_x_x[i] = 4.0 * g_x_xz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_x_y[i] = 4.0 * g_x_xz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_x_z[i] = 4.0 * g_x_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_0_y_0_0_xz_y_x, g_x_0_y_0_0_xz_y_y, g_x_0_y_0_0_xz_y_z, g_x_xz_0_x, g_x_xz_0_y, g_x_xz_0_z, g_x_xz_yy_x, g_x_xz_yy_y, g_x_xz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xz_y_x[i] = -2.0 * g_x_xz_0_x[i] * a_exp + 4.0 * g_x_xz_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_y_y[i] = -2.0 * g_x_xz_0_y[i] * a_exp + 4.0 * g_x_xz_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_y_z[i] = -2.0 * g_x_xz_0_z[i] * a_exp + 4.0 * g_x_xz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_0_y_0_0_xz_z_x, g_x_0_y_0_0_xz_z_y, g_x_0_y_0_0_xz_z_z, g_x_xz_yz_x, g_x_xz_yz_y, g_x_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xz_z_x[i] = 4.0 * g_x_xz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_z_y[i] = 4.0 * g_x_xz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_z_z[i] = 4.0 * g_x_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_x_0_y_0_0_yy_x_x, g_x_0_y_0_0_yy_x_y, g_x_0_y_0_0_yy_x_z, g_x_yy_xy_x, g_x_yy_xy_y, g_x_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yy_x_x[i] = 4.0 * g_x_yy_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_x_y[i] = 4.0 * g_x_yy_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_x_z[i] = 4.0 * g_x_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_x_0_y_0_0_yy_y_x, g_x_0_y_0_0_yy_y_y, g_x_0_y_0_0_yy_y_z, g_x_yy_0_x, g_x_yy_0_y, g_x_yy_0_z, g_x_yy_yy_x, g_x_yy_yy_y, g_x_yy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yy_y_x[i] = -2.0 * g_x_yy_0_x[i] * a_exp + 4.0 * g_x_yy_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_y_y[i] = -2.0 * g_x_yy_0_y[i] * a_exp + 4.0 * g_x_yy_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_y_z[i] = -2.0 * g_x_yy_0_z[i] * a_exp + 4.0 * g_x_yy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_x_0_y_0_0_yy_z_x, g_x_0_y_0_0_yy_z_y, g_x_0_y_0_0_yy_z_z, g_x_yy_yz_x, g_x_yy_yz_y, g_x_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yy_z_x[i] = 4.0 * g_x_yy_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_z_y[i] = 4.0 * g_x_yy_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_z_z[i] = 4.0 * g_x_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_x_0_y_0_0_yz_x_x, g_x_0_y_0_0_yz_x_y, g_x_0_y_0_0_yz_x_z, g_x_yz_xy_x, g_x_yz_xy_y, g_x_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yz_x_x[i] = 4.0 * g_x_yz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_x_y[i] = 4.0 * g_x_yz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_x_z[i] = 4.0 * g_x_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_x_0_y_0_0_yz_y_x, g_x_0_y_0_0_yz_y_y, g_x_0_y_0_0_yz_y_z, g_x_yz_0_x, g_x_yz_0_y, g_x_yz_0_z, g_x_yz_yy_x, g_x_yz_yy_y, g_x_yz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yz_y_x[i] = -2.0 * g_x_yz_0_x[i] * a_exp + 4.0 * g_x_yz_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_y_y[i] = -2.0 * g_x_yz_0_y[i] * a_exp + 4.0 * g_x_yz_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_y_z[i] = -2.0 * g_x_yz_0_z[i] * a_exp + 4.0 * g_x_yz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_x_0_y_0_0_yz_z_x, g_x_0_y_0_0_yz_z_y, g_x_0_y_0_0_yz_z_z, g_x_yz_yz_x, g_x_yz_yz_y, g_x_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yz_z_x[i] = 4.0 * g_x_yz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_z_y[i] = 4.0 * g_x_yz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_z_z[i] = 4.0 * g_x_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_x_0_y_0_0_zz_x_x, g_x_0_y_0_0_zz_x_y, g_x_0_y_0_0_zz_x_z, g_x_zz_xy_x, g_x_zz_xy_y, g_x_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_zz_x_x[i] = 4.0 * g_x_zz_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_x_y[i] = 4.0 * g_x_zz_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_x_z[i] = 4.0 * g_x_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_x_0_y_0_0_zz_y_x, g_x_0_y_0_0_zz_y_y, g_x_0_y_0_0_zz_y_z, g_x_zz_0_x, g_x_zz_0_y, g_x_zz_0_z, g_x_zz_yy_x, g_x_zz_yy_y, g_x_zz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_zz_y_x[i] = -2.0 * g_x_zz_0_x[i] * a_exp + 4.0 * g_x_zz_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_y_y[i] = -2.0 * g_x_zz_0_y[i] * a_exp + 4.0 * g_x_zz_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_y_z[i] = -2.0 * g_x_zz_0_z[i] * a_exp + 4.0 * g_x_zz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_x_0_y_0_0_zz_z_x, g_x_0_y_0_0_zz_z_y, g_x_0_y_0_0_zz_z_z, g_x_zz_yz_x, g_x_zz_yz_y, g_x_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_zz_z_x[i] = 4.0 * g_x_zz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_z_y[i] = 4.0 * g_x_zz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_z_z[i] = 4.0 * g_x_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_x_0_z_0_0_xx_x_x, g_x_0_z_0_0_xx_x_y, g_x_0_z_0_0_xx_x_z, g_x_xx_xz_x, g_x_xx_xz_y, g_x_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xx_x_x[i] = 4.0 * g_x_xx_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_x_y[i] = 4.0 * g_x_xx_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_x_z[i] = 4.0 * g_x_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_x_0_z_0_0_xx_y_x, g_x_0_z_0_0_xx_y_y, g_x_0_z_0_0_xx_y_z, g_x_xx_yz_x, g_x_xx_yz_y, g_x_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xx_y_x[i] = 4.0 * g_x_xx_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_y_y[i] = 4.0 * g_x_xx_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_y_z[i] = 4.0 * g_x_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_x_0_z_0_0_xx_z_x, g_x_0_z_0_0_xx_z_y, g_x_0_z_0_0_xx_z_z, g_x_xx_0_x, g_x_xx_0_y, g_x_xx_0_z, g_x_xx_zz_x, g_x_xx_zz_y, g_x_xx_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xx_z_x[i] = -2.0 * g_x_xx_0_x[i] * a_exp + 4.0 * g_x_xx_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_z_y[i] = -2.0 * g_x_xx_0_y[i] * a_exp + 4.0 * g_x_xx_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_z_z[i] = -2.0 * g_x_xx_0_z[i] * a_exp + 4.0 * g_x_xx_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_x_0_z_0_0_xy_x_x, g_x_0_z_0_0_xy_x_y, g_x_0_z_0_0_xy_x_z, g_x_xy_xz_x, g_x_xy_xz_y, g_x_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xy_x_x[i] = 4.0 * g_x_xy_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_x_y[i] = 4.0 * g_x_xy_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_x_z[i] = 4.0 * g_x_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_x_0_z_0_0_xy_y_x, g_x_0_z_0_0_xy_y_y, g_x_0_z_0_0_xy_y_z, g_x_xy_yz_x, g_x_xy_yz_y, g_x_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xy_y_x[i] = 4.0 * g_x_xy_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_y_y[i] = 4.0 * g_x_xy_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_y_z[i] = 4.0 * g_x_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_x_0_z_0_0_xy_z_x, g_x_0_z_0_0_xy_z_y, g_x_0_z_0_0_xy_z_z, g_x_xy_0_x, g_x_xy_0_y, g_x_xy_0_z, g_x_xy_zz_x, g_x_xy_zz_y, g_x_xy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xy_z_x[i] = -2.0 * g_x_xy_0_x[i] * a_exp + 4.0 * g_x_xy_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_z_y[i] = -2.0 * g_x_xy_0_y[i] * a_exp + 4.0 * g_x_xy_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_z_z[i] = -2.0 * g_x_xy_0_z[i] * a_exp + 4.0 * g_x_xy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_x_0_z_0_0_xz_x_x, g_x_0_z_0_0_xz_x_y, g_x_0_z_0_0_xz_x_z, g_x_xz_xz_x, g_x_xz_xz_y, g_x_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xz_x_x[i] = 4.0 * g_x_xz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_x_y[i] = 4.0 * g_x_xz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_x_z[i] = 4.0 * g_x_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_x_0_z_0_0_xz_y_x, g_x_0_z_0_0_xz_y_y, g_x_0_z_0_0_xz_y_z, g_x_xz_yz_x, g_x_xz_yz_y, g_x_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xz_y_x[i] = 4.0 * g_x_xz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_y_y[i] = 4.0 * g_x_xz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_y_z[i] = 4.0 * g_x_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_x_0_z_0_0_xz_z_x, g_x_0_z_0_0_xz_z_y, g_x_0_z_0_0_xz_z_z, g_x_xz_0_x, g_x_xz_0_y, g_x_xz_0_z, g_x_xz_zz_x, g_x_xz_zz_y, g_x_xz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xz_z_x[i] = -2.0 * g_x_xz_0_x[i] * a_exp + 4.0 * g_x_xz_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_z_y[i] = -2.0 * g_x_xz_0_y[i] * a_exp + 4.0 * g_x_xz_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_z_z[i] = -2.0 * g_x_xz_0_z[i] * a_exp + 4.0 * g_x_xz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_x_0_z_0_0_yy_x_x, g_x_0_z_0_0_yy_x_y, g_x_0_z_0_0_yy_x_z, g_x_yy_xz_x, g_x_yy_xz_y, g_x_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yy_x_x[i] = 4.0 * g_x_yy_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_x_y[i] = 4.0 * g_x_yy_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_x_z[i] = 4.0 * g_x_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_x_0_z_0_0_yy_y_x, g_x_0_z_0_0_yy_y_y, g_x_0_z_0_0_yy_y_z, g_x_yy_yz_x, g_x_yy_yz_y, g_x_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yy_y_x[i] = 4.0 * g_x_yy_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_y_y[i] = 4.0 * g_x_yy_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_y_z[i] = 4.0 * g_x_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_x_0_z_0_0_yy_z_x, g_x_0_z_0_0_yy_z_y, g_x_0_z_0_0_yy_z_z, g_x_yy_0_x, g_x_yy_0_y, g_x_yy_0_z, g_x_yy_zz_x, g_x_yy_zz_y, g_x_yy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yy_z_x[i] = -2.0 * g_x_yy_0_x[i] * a_exp + 4.0 * g_x_yy_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_z_y[i] = -2.0 * g_x_yy_0_y[i] * a_exp + 4.0 * g_x_yy_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_z_z[i] = -2.0 * g_x_yy_0_z[i] * a_exp + 4.0 * g_x_yy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_x_0_z_0_0_yz_x_x, g_x_0_z_0_0_yz_x_y, g_x_0_z_0_0_yz_x_z, g_x_yz_xz_x, g_x_yz_xz_y, g_x_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yz_x_x[i] = 4.0 * g_x_yz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_x_y[i] = 4.0 * g_x_yz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_x_z[i] = 4.0 * g_x_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_x_0_z_0_0_yz_y_x, g_x_0_z_0_0_yz_y_y, g_x_0_z_0_0_yz_y_z, g_x_yz_yz_x, g_x_yz_yz_y, g_x_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yz_y_x[i] = 4.0 * g_x_yz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_y_y[i] = 4.0 * g_x_yz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_y_z[i] = 4.0 * g_x_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_x_0_z_0_0_yz_z_x, g_x_0_z_0_0_yz_z_y, g_x_0_z_0_0_yz_z_z, g_x_yz_0_x, g_x_yz_0_y, g_x_yz_0_z, g_x_yz_zz_x, g_x_yz_zz_y, g_x_yz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yz_z_x[i] = -2.0 * g_x_yz_0_x[i] * a_exp + 4.0 * g_x_yz_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_z_y[i] = -2.0 * g_x_yz_0_y[i] * a_exp + 4.0 * g_x_yz_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_z_z[i] = -2.0 * g_x_yz_0_z[i] * a_exp + 4.0 * g_x_yz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_x_0_z_0_0_zz_x_x, g_x_0_z_0_0_zz_x_y, g_x_0_z_0_0_zz_x_z, g_x_zz_xz_x, g_x_zz_xz_y, g_x_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_zz_x_x[i] = 4.0 * g_x_zz_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_x_y[i] = 4.0 * g_x_zz_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_x_z[i] = 4.0 * g_x_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_x_0_z_0_0_zz_y_x, g_x_0_z_0_0_zz_y_y, g_x_0_z_0_0_zz_y_z, g_x_zz_yz_x, g_x_zz_yz_y, g_x_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_zz_y_x[i] = 4.0 * g_x_zz_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_y_y[i] = 4.0 * g_x_zz_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_y_z[i] = 4.0 * g_x_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_x_0_z_0_0_zz_z_x, g_x_0_z_0_0_zz_z_y, g_x_0_z_0_0_zz_z_z, g_x_zz_0_x, g_x_zz_0_y, g_x_zz_0_z, g_x_zz_zz_x, g_x_zz_zz_y, g_x_zz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_zz_z_x[i] = -2.0 * g_x_zz_0_x[i] * a_exp + 4.0 * g_x_zz_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_z_y[i] = -2.0 * g_x_zz_0_y[i] * a_exp + 4.0 * g_x_zz_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_z_z[i] = -2.0 * g_x_zz_0_z[i] * a_exp + 4.0 * g_x_zz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_y_0_x_0_0_xx_x_x, g_y_0_x_0_0_xx_x_y, g_y_0_x_0_0_xx_x_z, g_y_xx_0_x, g_y_xx_0_y, g_y_xx_0_z, g_y_xx_xx_x, g_y_xx_xx_y, g_y_xx_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xx_x_x[i] = -2.0 * g_y_xx_0_x[i] * a_exp + 4.0 * g_y_xx_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_x_y[i] = -2.0 * g_y_xx_0_y[i] * a_exp + 4.0 * g_y_xx_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_x_z[i] = -2.0 * g_y_xx_0_z[i] * a_exp + 4.0 * g_y_xx_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_y_0_x_0_0_xx_y_x, g_y_0_x_0_0_xx_y_y, g_y_0_x_0_0_xx_y_z, g_y_xx_xy_x, g_y_xx_xy_y, g_y_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xx_y_x[i] = 4.0 * g_y_xx_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_y_y[i] = 4.0 * g_y_xx_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_y_z[i] = 4.0 * g_y_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_y_0_x_0_0_xx_z_x, g_y_0_x_0_0_xx_z_y, g_y_0_x_0_0_xx_z_z, g_y_xx_xz_x, g_y_xx_xz_y, g_y_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xx_z_x[i] = 4.0 * g_y_xx_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_z_y[i] = 4.0 * g_y_xx_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_z_z[i] = 4.0 * g_y_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_y_0_x_0_0_xy_x_x, g_y_0_x_0_0_xy_x_y, g_y_0_x_0_0_xy_x_z, g_y_xy_0_x, g_y_xy_0_y, g_y_xy_0_z, g_y_xy_xx_x, g_y_xy_xx_y, g_y_xy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xy_x_x[i] = -2.0 * g_y_xy_0_x[i] * a_exp + 4.0 * g_y_xy_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_x_y[i] = -2.0 * g_y_xy_0_y[i] * a_exp + 4.0 * g_y_xy_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_x_z[i] = -2.0 * g_y_xy_0_z[i] * a_exp + 4.0 * g_y_xy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_y_0_x_0_0_xy_y_x, g_y_0_x_0_0_xy_y_y, g_y_0_x_0_0_xy_y_z, g_y_xy_xy_x, g_y_xy_xy_y, g_y_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xy_y_x[i] = 4.0 * g_y_xy_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_y_y[i] = 4.0 * g_y_xy_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_y_z[i] = 4.0 * g_y_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_y_0_x_0_0_xy_z_x, g_y_0_x_0_0_xy_z_y, g_y_0_x_0_0_xy_z_z, g_y_xy_xz_x, g_y_xy_xz_y, g_y_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xy_z_x[i] = 4.0 * g_y_xy_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_z_y[i] = 4.0 * g_y_xy_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_z_z[i] = 4.0 * g_y_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_y_0_x_0_0_xz_x_x, g_y_0_x_0_0_xz_x_y, g_y_0_x_0_0_xz_x_z, g_y_xz_0_x, g_y_xz_0_y, g_y_xz_0_z, g_y_xz_xx_x, g_y_xz_xx_y, g_y_xz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xz_x_x[i] = -2.0 * g_y_xz_0_x[i] * a_exp + 4.0 * g_y_xz_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_x_y[i] = -2.0 * g_y_xz_0_y[i] * a_exp + 4.0 * g_y_xz_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_x_z[i] = -2.0 * g_y_xz_0_z[i] * a_exp + 4.0 * g_y_xz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_y_0_x_0_0_xz_y_x, g_y_0_x_0_0_xz_y_y, g_y_0_x_0_0_xz_y_z, g_y_xz_xy_x, g_y_xz_xy_y, g_y_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xz_y_x[i] = 4.0 * g_y_xz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_y_y[i] = 4.0 * g_y_xz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_y_z[i] = 4.0 * g_y_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_y_0_x_0_0_xz_z_x, g_y_0_x_0_0_xz_z_y, g_y_0_x_0_0_xz_z_z, g_y_xz_xz_x, g_y_xz_xz_y, g_y_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xz_z_x[i] = 4.0 * g_y_xz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_z_y[i] = 4.0 * g_y_xz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_z_z[i] = 4.0 * g_y_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_y_0_x_0_0_yy_x_x, g_y_0_x_0_0_yy_x_y, g_y_0_x_0_0_yy_x_z, g_y_yy_0_x, g_y_yy_0_y, g_y_yy_0_z, g_y_yy_xx_x, g_y_yy_xx_y, g_y_yy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yy_x_x[i] = -2.0 * g_y_yy_0_x[i] * a_exp + 4.0 * g_y_yy_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_x_y[i] = -2.0 * g_y_yy_0_y[i] * a_exp + 4.0 * g_y_yy_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_x_z[i] = -2.0 * g_y_yy_0_z[i] * a_exp + 4.0 * g_y_yy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_y_0_x_0_0_yy_y_x, g_y_0_x_0_0_yy_y_y, g_y_0_x_0_0_yy_y_z, g_y_yy_xy_x, g_y_yy_xy_y, g_y_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yy_y_x[i] = 4.0 * g_y_yy_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_y_y[i] = 4.0 * g_y_yy_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_y_z[i] = 4.0 * g_y_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_y_0_x_0_0_yy_z_x, g_y_0_x_0_0_yy_z_y, g_y_0_x_0_0_yy_z_z, g_y_yy_xz_x, g_y_yy_xz_y, g_y_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yy_z_x[i] = 4.0 * g_y_yy_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_z_y[i] = 4.0 * g_y_yy_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_z_z[i] = 4.0 * g_y_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_y_0_x_0_0_yz_x_x, g_y_0_x_0_0_yz_x_y, g_y_0_x_0_0_yz_x_z, g_y_yz_0_x, g_y_yz_0_y, g_y_yz_0_z, g_y_yz_xx_x, g_y_yz_xx_y, g_y_yz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yz_x_x[i] = -2.0 * g_y_yz_0_x[i] * a_exp + 4.0 * g_y_yz_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_x_y[i] = -2.0 * g_y_yz_0_y[i] * a_exp + 4.0 * g_y_yz_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_x_z[i] = -2.0 * g_y_yz_0_z[i] * a_exp + 4.0 * g_y_yz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_y_0_x_0_0_yz_y_x, g_y_0_x_0_0_yz_y_y, g_y_0_x_0_0_yz_y_z, g_y_yz_xy_x, g_y_yz_xy_y, g_y_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yz_y_x[i] = 4.0 * g_y_yz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_y_y[i] = 4.0 * g_y_yz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_y_z[i] = 4.0 * g_y_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_y_0_x_0_0_yz_z_x, g_y_0_x_0_0_yz_z_y, g_y_0_x_0_0_yz_z_z, g_y_yz_xz_x, g_y_yz_xz_y, g_y_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yz_z_x[i] = 4.0 * g_y_yz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_z_y[i] = 4.0 * g_y_yz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_z_z[i] = 4.0 * g_y_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_y_0_x_0_0_zz_x_x, g_y_0_x_0_0_zz_x_y, g_y_0_x_0_0_zz_x_z, g_y_zz_0_x, g_y_zz_0_y, g_y_zz_0_z, g_y_zz_xx_x, g_y_zz_xx_y, g_y_zz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_zz_x_x[i] = -2.0 * g_y_zz_0_x[i] * a_exp + 4.0 * g_y_zz_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_x_y[i] = -2.0 * g_y_zz_0_y[i] * a_exp + 4.0 * g_y_zz_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_x_z[i] = -2.0 * g_y_zz_0_z[i] * a_exp + 4.0 * g_y_zz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_y_0_x_0_0_zz_y_x, g_y_0_x_0_0_zz_y_y, g_y_0_x_0_0_zz_y_z, g_y_zz_xy_x, g_y_zz_xy_y, g_y_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_zz_y_x[i] = 4.0 * g_y_zz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_y_y[i] = 4.0 * g_y_zz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_y_z[i] = 4.0 * g_y_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_y_0_x_0_0_zz_z_x, g_y_0_x_0_0_zz_z_y, g_y_0_x_0_0_zz_z_z, g_y_zz_xz_x, g_y_zz_xz_y, g_y_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_zz_z_x[i] = 4.0 * g_y_zz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_z_y[i] = 4.0 * g_y_zz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_z_z[i] = 4.0 * g_y_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_y_0_y_0_0_xx_x_x, g_y_0_y_0_0_xx_x_y, g_y_0_y_0_0_xx_x_z, g_y_xx_xy_x, g_y_xx_xy_y, g_y_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xx_x_x[i] = 4.0 * g_y_xx_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_x_y[i] = 4.0 * g_y_xx_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_x_z[i] = 4.0 * g_y_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_y_0_y_0_0_xx_y_x, g_y_0_y_0_0_xx_y_y, g_y_0_y_0_0_xx_y_z, g_y_xx_0_x, g_y_xx_0_y, g_y_xx_0_z, g_y_xx_yy_x, g_y_xx_yy_y, g_y_xx_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xx_y_x[i] = -2.0 * g_y_xx_0_x[i] * a_exp + 4.0 * g_y_xx_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_y_y[i] = -2.0 * g_y_xx_0_y[i] * a_exp + 4.0 * g_y_xx_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_y_z[i] = -2.0 * g_y_xx_0_z[i] * a_exp + 4.0 * g_y_xx_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_y_0_y_0_0_xx_z_x, g_y_0_y_0_0_xx_z_y, g_y_0_y_0_0_xx_z_z, g_y_xx_yz_x, g_y_xx_yz_y, g_y_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xx_z_x[i] = 4.0 * g_y_xx_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_z_y[i] = 4.0 * g_y_xx_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_z_z[i] = 4.0 * g_y_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_y_0_y_0_0_xy_x_x, g_y_0_y_0_0_xy_x_y, g_y_0_y_0_0_xy_x_z, g_y_xy_xy_x, g_y_xy_xy_y, g_y_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xy_x_x[i] = 4.0 * g_y_xy_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_x_y[i] = 4.0 * g_y_xy_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_x_z[i] = 4.0 * g_y_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_y_0_y_0_0_xy_y_x, g_y_0_y_0_0_xy_y_y, g_y_0_y_0_0_xy_y_z, g_y_xy_0_x, g_y_xy_0_y, g_y_xy_0_z, g_y_xy_yy_x, g_y_xy_yy_y, g_y_xy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xy_y_x[i] = -2.0 * g_y_xy_0_x[i] * a_exp + 4.0 * g_y_xy_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_y_y[i] = -2.0 * g_y_xy_0_y[i] * a_exp + 4.0 * g_y_xy_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_y_z[i] = -2.0 * g_y_xy_0_z[i] * a_exp + 4.0 * g_y_xy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_y_0_y_0_0_xy_z_x, g_y_0_y_0_0_xy_z_y, g_y_0_y_0_0_xy_z_z, g_y_xy_yz_x, g_y_xy_yz_y, g_y_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xy_z_x[i] = 4.0 * g_y_xy_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_z_y[i] = 4.0 * g_y_xy_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_z_z[i] = 4.0 * g_y_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_y_0_y_0_0_xz_x_x, g_y_0_y_0_0_xz_x_y, g_y_0_y_0_0_xz_x_z, g_y_xz_xy_x, g_y_xz_xy_y, g_y_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xz_x_x[i] = 4.0 * g_y_xz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_x_y[i] = 4.0 * g_y_xz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_x_z[i] = 4.0 * g_y_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_y_0_y_0_0_xz_y_x, g_y_0_y_0_0_xz_y_y, g_y_0_y_0_0_xz_y_z, g_y_xz_0_x, g_y_xz_0_y, g_y_xz_0_z, g_y_xz_yy_x, g_y_xz_yy_y, g_y_xz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xz_y_x[i] = -2.0 * g_y_xz_0_x[i] * a_exp + 4.0 * g_y_xz_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_y_y[i] = -2.0 * g_y_xz_0_y[i] * a_exp + 4.0 * g_y_xz_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_y_z[i] = -2.0 * g_y_xz_0_z[i] * a_exp + 4.0 * g_y_xz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_y_0_y_0_0_xz_z_x, g_y_0_y_0_0_xz_z_y, g_y_0_y_0_0_xz_z_z, g_y_xz_yz_x, g_y_xz_yz_y, g_y_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xz_z_x[i] = 4.0 * g_y_xz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_z_y[i] = 4.0 * g_y_xz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_z_z[i] = 4.0 * g_y_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_y_0_y_0_0_yy_x_x, g_y_0_y_0_0_yy_x_y, g_y_0_y_0_0_yy_x_z, g_y_yy_xy_x, g_y_yy_xy_y, g_y_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yy_x_x[i] = 4.0 * g_y_yy_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_x_y[i] = 4.0 * g_y_yy_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_x_z[i] = 4.0 * g_y_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_y_0_y_0_0_yy_y_x, g_y_0_y_0_0_yy_y_y, g_y_0_y_0_0_yy_y_z, g_y_yy_0_x, g_y_yy_0_y, g_y_yy_0_z, g_y_yy_yy_x, g_y_yy_yy_y, g_y_yy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yy_y_x[i] = -2.0 * g_y_yy_0_x[i] * a_exp + 4.0 * g_y_yy_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_y_y[i] = -2.0 * g_y_yy_0_y[i] * a_exp + 4.0 * g_y_yy_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_y_z[i] = -2.0 * g_y_yy_0_z[i] * a_exp + 4.0 * g_y_yy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_y_0_y_0_0_yy_z_x, g_y_0_y_0_0_yy_z_y, g_y_0_y_0_0_yy_z_z, g_y_yy_yz_x, g_y_yy_yz_y, g_y_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yy_z_x[i] = 4.0 * g_y_yy_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_z_y[i] = 4.0 * g_y_yy_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_z_z[i] = 4.0 * g_y_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_y_0_y_0_0_yz_x_x, g_y_0_y_0_0_yz_x_y, g_y_0_y_0_0_yz_x_z, g_y_yz_xy_x, g_y_yz_xy_y, g_y_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yz_x_x[i] = 4.0 * g_y_yz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_x_y[i] = 4.0 * g_y_yz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_x_z[i] = 4.0 * g_y_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_y_0_y_0_0_yz_y_x, g_y_0_y_0_0_yz_y_y, g_y_0_y_0_0_yz_y_z, g_y_yz_0_x, g_y_yz_0_y, g_y_yz_0_z, g_y_yz_yy_x, g_y_yz_yy_y, g_y_yz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yz_y_x[i] = -2.0 * g_y_yz_0_x[i] * a_exp + 4.0 * g_y_yz_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_y_y[i] = -2.0 * g_y_yz_0_y[i] * a_exp + 4.0 * g_y_yz_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_y_z[i] = -2.0 * g_y_yz_0_z[i] * a_exp + 4.0 * g_y_yz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_y_0_y_0_0_yz_z_x, g_y_0_y_0_0_yz_z_y, g_y_0_y_0_0_yz_z_z, g_y_yz_yz_x, g_y_yz_yz_y, g_y_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yz_z_x[i] = 4.0 * g_y_yz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_z_y[i] = 4.0 * g_y_yz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_z_z[i] = 4.0 * g_y_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_y_0_y_0_0_zz_x_x, g_y_0_y_0_0_zz_x_y, g_y_0_y_0_0_zz_x_z, g_y_zz_xy_x, g_y_zz_xy_y, g_y_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_zz_x_x[i] = 4.0 * g_y_zz_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_x_y[i] = 4.0 * g_y_zz_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_x_z[i] = 4.0 * g_y_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_y_0_y_0_0_zz_y_x, g_y_0_y_0_0_zz_y_y, g_y_0_y_0_0_zz_y_z, g_y_zz_0_x, g_y_zz_0_y, g_y_zz_0_z, g_y_zz_yy_x, g_y_zz_yy_y, g_y_zz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_zz_y_x[i] = -2.0 * g_y_zz_0_x[i] * a_exp + 4.0 * g_y_zz_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_y_y[i] = -2.0 * g_y_zz_0_y[i] * a_exp + 4.0 * g_y_zz_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_y_z[i] = -2.0 * g_y_zz_0_z[i] * a_exp + 4.0 * g_y_zz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_y_0_y_0_0_zz_z_x, g_y_0_y_0_0_zz_z_y, g_y_0_y_0_0_zz_z_z, g_y_zz_yz_x, g_y_zz_yz_y, g_y_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_zz_z_x[i] = 4.0 * g_y_zz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_z_y[i] = 4.0 * g_y_zz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_z_z[i] = 4.0 * g_y_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_y_0_z_0_0_xx_x_x, g_y_0_z_0_0_xx_x_y, g_y_0_z_0_0_xx_x_z, g_y_xx_xz_x, g_y_xx_xz_y, g_y_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xx_x_x[i] = 4.0 * g_y_xx_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_x_y[i] = 4.0 * g_y_xx_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_x_z[i] = 4.0 * g_y_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_y_0_z_0_0_xx_y_x, g_y_0_z_0_0_xx_y_y, g_y_0_z_0_0_xx_y_z, g_y_xx_yz_x, g_y_xx_yz_y, g_y_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xx_y_x[i] = 4.0 * g_y_xx_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_y_y[i] = 4.0 * g_y_xx_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_y_z[i] = 4.0 * g_y_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_y_0_z_0_0_xx_z_x, g_y_0_z_0_0_xx_z_y, g_y_0_z_0_0_xx_z_z, g_y_xx_0_x, g_y_xx_0_y, g_y_xx_0_z, g_y_xx_zz_x, g_y_xx_zz_y, g_y_xx_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xx_z_x[i] = -2.0 * g_y_xx_0_x[i] * a_exp + 4.0 * g_y_xx_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_z_y[i] = -2.0 * g_y_xx_0_y[i] * a_exp + 4.0 * g_y_xx_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_z_z[i] = -2.0 * g_y_xx_0_z[i] * a_exp + 4.0 * g_y_xx_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_y_0_z_0_0_xy_x_x, g_y_0_z_0_0_xy_x_y, g_y_0_z_0_0_xy_x_z, g_y_xy_xz_x, g_y_xy_xz_y, g_y_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xy_x_x[i] = 4.0 * g_y_xy_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_x_y[i] = 4.0 * g_y_xy_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_x_z[i] = 4.0 * g_y_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_y_0_z_0_0_xy_y_x, g_y_0_z_0_0_xy_y_y, g_y_0_z_0_0_xy_y_z, g_y_xy_yz_x, g_y_xy_yz_y, g_y_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xy_y_x[i] = 4.0 * g_y_xy_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_y_y[i] = 4.0 * g_y_xy_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_y_z[i] = 4.0 * g_y_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_y_0_z_0_0_xy_z_x, g_y_0_z_0_0_xy_z_y, g_y_0_z_0_0_xy_z_z, g_y_xy_0_x, g_y_xy_0_y, g_y_xy_0_z, g_y_xy_zz_x, g_y_xy_zz_y, g_y_xy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xy_z_x[i] = -2.0 * g_y_xy_0_x[i] * a_exp + 4.0 * g_y_xy_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_z_y[i] = -2.0 * g_y_xy_0_y[i] * a_exp + 4.0 * g_y_xy_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_z_z[i] = -2.0 * g_y_xy_0_z[i] * a_exp + 4.0 * g_y_xy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_y_0_z_0_0_xz_x_x, g_y_0_z_0_0_xz_x_y, g_y_0_z_0_0_xz_x_z, g_y_xz_xz_x, g_y_xz_xz_y, g_y_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xz_x_x[i] = 4.0 * g_y_xz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_x_y[i] = 4.0 * g_y_xz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_x_z[i] = 4.0 * g_y_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_y_0_z_0_0_xz_y_x, g_y_0_z_0_0_xz_y_y, g_y_0_z_0_0_xz_y_z, g_y_xz_yz_x, g_y_xz_yz_y, g_y_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xz_y_x[i] = 4.0 * g_y_xz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_y_y[i] = 4.0 * g_y_xz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_y_z[i] = 4.0 * g_y_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_y_0_z_0_0_xz_z_x, g_y_0_z_0_0_xz_z_y, g_y_0_z_0_0_xz_z_z, g_y_xz_0_x, g_y_xz_0_y, g_y_xz_0_z, g_y_xz_zz_x, g_y_xz_zz_y, g_y_xz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xz_z_x[i] = -2.0 * g_y_xz_0_x[i] * a_exp + 4.0 * g_y_xz_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_z_y[i] = -2.0 * g_y_xz_0_y[i] * a_exp + 4.0 * g_y_xz_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_z_z[i] = -2.0 * g_y_xz_0_z[i] * a_exp + 4.0 * g_y_xz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_y_0_z_0_0_yy_x_x, g_y_0_z_0_0_yy_x_y, g_y_0_z_0_0_yy_x_z, g_y_yy_xz_x, g_y_yy_xz_y, g_y_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yy_x_x[i] = 4.0 * g_y_yy_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_x_y[i] = 4.0 * g_y_yy_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_x_z[i] = 4.0 * g_y_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_y_0_z_0_0_yy_y_x, g_y_0_z_0_0_yy_y_y, g_y_0_z_0_0_yy_y_z, g_y_yy_yz_x, g_y_yy_yz_y, g_y_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yy_y_x[i] = 4.0 * g_y_yy_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_y_y[i] = 4.0 * g_y_yy_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_y_z[i] = 4.0 * g_y_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_y_0_z_0_0_yy_z_x, g_y_0_z_0_0_yy_z_y, g_y_0_z_0_0_yy_z_z, g_y_yy_0_x, g_y_yy_0_y, g_y_yy_0_z, g_y_yy_zz_x, g_y_yy_zz_y, g_y_yy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yy_z_x[i] = -2.0 * g_y_yy_0_x[i] * a_exp + 4.0 * g_y_yy_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_z_y[i] = -2.0 * g_y_yy_0_y[i] * a_exp + 4.0 * g_y_yy_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_z_z[i] = -2.0 * g_y_yy_0_z[i] * a_exp + 4.0 * g_y_yy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_y_0_z_0_0_yz_x_x, g_y_0_z_0_0_yz_x_y, g_y_0_z_0_0_yz_x_z, g_y_yz_xz_x, g_y_yz_xz_y, g_y_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yz_x_x[i] = 4.0 * g_y_yz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_x_y[i] = 4.0 * g_y_yz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_x_z[i] = 4.0 * g_y_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_y_0_z_0_0_yz_y_x, g_y_0_z_0_0_yz_y_y, g_y_0_z_0_0_yz_y_z, g_y_yz_yz_x, g_y_yz_yz_y, g_y_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yz_y_x[i] = 4.0 * g_y_yz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_y_y[i] = 4.0 * g_y_yz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_y_z[i] = 4.0 * g_y_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_y_0_z_0_0_yz_z_x, g_y_0_z_0_0_yz_z_y, g_y_0_z_0_0_yz_z_z, g_y_yz_0_x, g_y_yz_0_y, g_y_yz_0_z, g_y_yz_zz_x, g_y_yz_zz_y, g_y_yz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yz_z_x[i] = -2.0 * g_y_yz_0_x[i] * a_exp + 4.0 * g_y_yz_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_z_y[i] = -2.0 * g_y_yz_0_y[i] * a_exp + 4.0 * g_y_yz_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_z_z[i] = -2.0 * g_y_yz_0_z[i] * a_exp + 4.0 * g_y_yz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_y_0_z_0_0_zz_x_x, g_y_0_z_0_0_zz_x_y, g_y_0_z_0_0_zz_x_z, g_y_zz_xz_x, g_y_zz_xz_y, g_y_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_zz_x_x[i] = 4.0 * g_y_zz_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_x_y[i] = 4.0 * g_y_zz_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_x_z[i] = 4.0 * g_y_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_y_0_z_0_0_zz_y_x, g_y_0_z_0_0_zz_y_y, g_y_0_z_0_0_zz_y_z, g_y_zz_yz_x, g_y_zz_yz_y, g_y_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_zz_y_x[i] = 4.0 * g_y_zz_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_y_y[i] = 4.0 * g_y_zz_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_y_z[i] = 4.0 * g_y_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_y_0_z_0_0_zz_z_x, g_y_0_z_0_0_zz_z_y, g_y_0_z_0_0_zz_z_z, g_y_zz_0_x, g_y_zz_0_y, g_y_zz_0_z, g_y_zz_zz_x, g_y_zz_zz_y, g_y_zz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_zz_z_x[i] = -2.0 * g_y_zz_0_x[i] * a_exp + 4.0 * g_y_zz_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_z_y[i] = -2.0 * g_y_zz_0_y[i] * a_exp + 4.0 * g_y_zz_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_z_z[i] = -2.0 * g_y_zz_0_z[i] * a_exp + 4.0 * g_y_zz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (324-327)

    #pragma omp simd aligned(g_z_0_x_0_0_xx_x_x, g_z_0_x_0_0_xx_x_y, g_z_0_x_0_0_xx_x_z, g_z_xx_0_x, g_z_xx_0_y, g_z_xx_0_z, g_z_xx_xx_x, g_z_xx_xx_y, g_z_xx_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xx_x_x[i] = -2.0 * g_z_xx_0_x[i] * a_exp + 4.0 * g_z_xx_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_x_y[i] = -2.0 * g_z_xx_0_y[i] * a_exp + 4.0 * g_z_xx_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_x_z[i] = -2.0 * g_z_xx_0_z[i] * a_exp + 4.0 * g_z_xx_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (327-330)

    #pragma omp simd aligned(g_z_0_x_0_0_xx_y_x, g_z_0_x_0_0_xx_y_y, g_z_0_x_0_0_xx_y_z, g_z_xx_xy_x, g_z_xx_xy_y, g_z_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xx_y_x[i] = 4.0 * g_z_xx_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_y_y[i] = 4.0 * g_z_xx_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_y_z[i] = 4.0 * g_z_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (330-333)

    #pragma omp simd aligned(g_z_0_x_0_0_xx_z_x, g_z_0_x_0_0_xx_z_y, g_z_0_x_0_0_xx_z_z, g_z_xx_xz_x, g_z_xx_xz_y, g_z_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xx_z_x[i] = 4.0 * g_z_xx_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_z_y[i] = 4.0 * g_z_xx_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_z_z[i] = 4.0 * g_z_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (333-336)

    #pragma omp simd aligned(g_z_0_x_0_0_xy_x_x, g_z_0_x_0_0_xy_x_y, g_z_0_x_0_0_xy_x_z, g_z_xy_0_x, g_z_xy_0_y, g_z_xy_0_z, g_z_xy_xx_x, g_z_xy_xx_y, g_z_xy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xy_x_x[i] = -2.0 * g_z_xy_0_x[i] * a_exp + 4.0 * g_z_xy_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_x_y[i] = -2.0 * g_z_xy_0_y[i] * a_exp + 4.0 * g_z_xy_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_x_z[i] = -2.0 * g_z_xy_0_z[i] * a_exp + 4.0 * g_z_xy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (336-339)

    #pragma omp simd aligned(g_z_0_x_0_0_xy_y_x, g_z_0_x_0_0_xy_y_y, g_z_0_x_0_0_xy_y_z, g_z_xy_xy_x, g_z_xy_xy_y, g_z_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xy_y_x[i] = 4.0 * g_z_xy_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_y_y[i] = 4.0 * g_z_xy_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_y_z[i] = 4.0 * g_z_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (339-342)

    #pragma omp simd aligned(g_z_0_x_0_0_xy_z_x, g_z_0_x_0_0_xy_z_y, g_z_0_x_0_0_xy_z_z, g_z_xy_xz_x, g_z_xy_xz_y, g_z_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xy_z_x[i] = 4.0 * g_z_xy_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_z_y[i] = 4.0 * g_z_xy_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_z_z[i] = 4.0 * g_z_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (342-345)

    #pragma omp simd aligned(g_z_0_x_0_0_xz_x_x, g_z_0_x_0_0_xz_x_y, g_z_0_x_0_0_xz_x_z, g_z_xz_0_x, g_z_xz_0_y, g_z_xz_0_z, g_z_xz_xx_x, g_z_xz_xx_y, g_z_xz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xz_x_x[i] = -2.0 * g_z_xz_0_x[i] * a_exp + 4.0 * g_z_xz_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_x_y[i] = -2.0 * g_z_xz_0_y[i] * a_exp + 4.0 * g_z_xz_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_x_z[i] = -2.0 * g_z_xz_0_z[i] * a_exp + 4.0 * g_z_xz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (345-348)

    #pragma omp simd aligned(g_z_0_x_0_0_xz_y_x, g_z_0_x_0_0_xz_y_y, g_z_0_x_0_0_xz_y_z, g_z_xz_xy_x, g_z_xz_xy_y, g_z_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xz_y_x[i] = 4.0 * g_z_xz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_y_y[i] = 4.0 * g_z_xz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_y_z[i] = 4.0 * g_z_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (348-351)

    #pragma omp simd aligned(g_z_0_x_0_0_xz_z_x, g_z_0_x_0_0_xz_z_y, g_z_0_x_0_0_xz_z_z, g_z_xz_xz_x, g_z_xz_xz_y, g_z_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xz_z_x[i] = 4.0 * g_z_xz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_z_y[i] = 4.0 * g_z_xz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_z_z[i] = 4.0 * g_z_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (351-354)

    #pragma omp simd aligned(g_z_0_x_0_0_yy_x_x, g_z_0_x_0_0_yy_x_y, g_z_0_x_0_0_yy_x_z, g_z_yy_0_x, g_z_yy_0_y, g_z_yy_0_z, g_z_yy_xx_x, g_z_yy_xx_y, g_z_yy_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yy_x_x[i] = -2.0 * g_z_yy_0_x[i] * a_exp + 4.0 * g_z_yy_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_x_y[i] = -2.0 * g_z_yy_0_y[i] * a_exp + 4.0 * g_z_yy_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_x_z[i] = -2.0 * g_z_yy_0_z[i] * a_exp + 4.0 * g_z_yy_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (354-357)

    #pragma omp simd aligned(g_z_0_x_0_0_yy_y_x, g_z_0_x_0_0_yy_y_y, g_z_0_x_0_0_yy_y_z, g_z_yy_xy_x, g_z_yy_xy_y, g_z_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yy_y_x[i] = 4.0 * g_z_yy_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_y_y[i] = 4.0 * g_z_yy_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_y_z[i] = 4.0 * g_z_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (357-360)

    #pragma omp simd aligned(g_z_0_x_0_0_yy_z_x, g_z_0_x_0_0_yy_z_y, g_z_0_x_0_0_yy_z_z, g_z_yy_xz_x, g_z_yy_xz_y, g_z_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yy_z_x[i] = 4.0 * g_z_yy_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_z_y[i] = 4.0 * g_z_yy_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_z_z[i] = 4.0 * g_z_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (360-363)

    #pragma omp simd aligned(g_z_0_x_0_0_yz_x_x, g_z_0_x_0_0_yz_x_y, g_z_0_x_0_0_yz_x_z, g_z_yz_0_x, g_z_yz_0_y, g_z_yz_0_z, g_z_yz_xx_x, g_z_yz_xx_y, g_z_yz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yz_x_x[i] = -2.0 * g_z_yz_0_x[i] * a_exp + 4.0 * g_z_yz_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_x_y[i] = -2.0 * g_z_yz_0_y[i] * a_exp + 4.0 * g_z_yz_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_x_z[i] = -2.0 * g_z_yz_0_z[i] * a_exp + 4.0 * g_z_yz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (363-366)

    #pragma omp simd aligned(g_z_0_x_0_0_yz_y_x, g_z_0_x_0_0_yz_y_y, g_z_0_x_0_0_yz_y_z, g_z_yz_xy_x, g_z_yz_xy_y, g_z_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yz_y_x[i] = 4.0 * g_z_yz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_y_y[i] = 4.0 * g_z_yz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_y_z[i] = 4.0 * g_z_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (366-369)

    #pragma omp simd aligned(g_z_0_x_0_0_yz_z_x, g_z_0_x_0_0_yz_z_y, g_z_0_x_0_0_yz_z_z, g_z_yz_xz_x, g_z_yz_xz_y, g_z_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yz_z_x[i] = 4.0 * g_z_yz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_z_y[i] = 4.0 * g_z_yz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_z_z[i] = 4.0 * g_z_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (369-372)

    #pragma omp simd aligned(g_z_0_x_0_0_zz_x_x, g_z_0_x_0_0_zz_x_y, g_z_0_x_0_0_zz_x_z, g_z_zz_0_x, g_z_zz_0_y, g_z_zz_0_z, g_z_zz_xx_x, g_z_zz_xx_y, g_z_zz_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_zz_x_x[i] = -2.0 * g_z_zz_0_x[i] * a_exp + 4.0 * g_z_zz_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_x_y[i] = -2.0 * g_z_zz_0_y[i] * a_exp + 4.0 * g_z_zz_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_x_z[i] = -2.0 * g_z_zz_0_z[i] * a_exp + 4.0 * g_z_zz_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (372-375)

    #pragma omp simd aligned(g_z_0_x_0_0_zz_y_x, g_z_0_x_0_0_zz_y_y, g_z_0_x_0_0_zz_y_z, g_z_zz_xy_x, g_z_zz_xy_y, g_z_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_zz_y_x[i] = 4.0 * g_z_zz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_y_y[i] = 4.0 * g_z_zz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_y_z[i] = 4.0 * g_z_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (375-378)

    #pragma omp simd aligned(g_z_0_x_0_0_zz_z_x, g_z_0_x_0_0_zz_z_y, g_z_0_x_0_0_zz_z_z, g_z_zz_xz_x, g_z_zz_xz_y, g_z_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_zz_z_x[i] = 4.0 * g_z_zz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_z_y[i] = 4.0 * g_z_zz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_z_z[i] = 4.0 * g_z_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (378-381)

    #pragma omp simd aligned(g_z_0_y_0_0_xx_x_x, g_z_0_y_0_0_xx_x_y, g_z_0_y_0_0_xx_x_z, g_z_xx_xy_x, g_z_xx_xy_y, g_z_xx_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xx_x_x[i] = 4.0 * g_z_xx_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_x_y[i] = 4.0 * g_z_xx_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_x_z[i] = 4.0 * g_z_xx_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (381-384)

    #pragma omp simd aligned(g_z_0_y_0_0_xx_y_x, g_z_0_y_0_0_xx_y_y, g_z_0_y_0_0_xx_y_z, g_z_xx_0_x, g_z_xx_0_y, g_z_xx_0_z, g_z_xx_yy_x, g_z_xx_yy_y, g_z_xx_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xx_y_x[i] = -2.0 * g_z_xx_0_x[i] * a_exp + 4.0 * g_z_xx_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_y_y[i] = -2.0 * g_z_xx_0_y[i] * a_exp + 4.0 * g_z_xx_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_y_z[i] = -2.0 * g_z_xx_0_z[i] * a_exp + 4.0 * g_z_xx_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (384-387)

    #pragma omp simd aligned(g_z_0_y_0_0_xx_z_x, g_z_0_y_0_0_xx_z_y, g_z_0_y_0_0_xx_z_z, g_z_xx_yz_x, g_z_xx_yz_y, g_z_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xx_z_x[i] = 4.0 * g_z_xx_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_z_y[i] = 4.0 * g_z_xx_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_z_z[i] = 4.0 * g_z_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (387-390)

    #pragma omp simd aligned(g_z_0_y_0_0_xy_x_x, g_z_0_y_0_0_xy_x_y, g_z_0_y_0_0_xy_x_z, g_z_xy_xy_x, g_z_xy_xy_y, g_z_xy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xy_x_x[i] = 4.0 * g_z_xy_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_x_y[i] = 4.0 * g_z_xy_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_x_z[i] = 4.0 * g_z_xy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (390-393)

    #pragma omp simd aligned(g_z_0_y_0_0_xy_y_x, g_z_0_y_0_0_xy_y_y, g_z_0_y_0_0_xy_y_z, g_z_xy_0_x, g_z_xy_0_y, g_z_xy_0_z, g_z_xy_yy_x, g_z_xy_yy_y, g_z_xy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xy_y_x[i] = -2.0 * g_z_xy_0_x[i] * a_exp + 4.0 * g_z_xy_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_y_y[i] = -2.0 * g_z_xy_0_y[i] * a_exp + 4.0 * g_z_xy_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_y_z[i] = -2.0 * g_z_xy_0_z[i] * a_exp + 4.0 * g_z_xy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (393-396)

    #pragma omp simd aligned(g_z_0_y_0_0_xy_z_x, g_z_0_y_0_0_xy_z_y, g_z_0_y_0_0_xy_z_z, g_z_xy_yz_x, g_z_xy_yz_y, g_z_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xy_z_x[i] = 4.0 * g_z_xy_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_z_y[i] = 4.0 * g_z_xy_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_z_z[i] = 4.0 * g_z_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (396-399)

    #pragma omp simd aligned(g_z_0_y_0_0_xz_x_x, g_z_0_y_0_0_xz_x_y, g_z_0_y_0_0_xz_x_z, g_z_xz_xy_x, g_z_xz_xy_y, g_z_xz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xz_x_x[i] = 4.0 * g_z_xz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_x_y[i] = 4.0 * g_z_xz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_x_z[i] = 4.0 * g_z_xz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (399-402)

    #pragma omp simd aligned(g_z_0_y_0_0_xz_y_x, g_z_0_y_0_0_xz_y_y, g_z_0_y_0_0_xz_y_z, g_z_xz_0_x, g_z_xz_0_y, g_z_xz_0_z, g_z_xz_yy_x, g_z_xz_yy_y, g_z_xz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xz_y_x[i] = -2.0 * g_z_xz_0_x[i] * a_exp + 4.0 * g_z_xz_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_y_y[i] = -2.0 * g_z_xz_0_y[i] * a_exp + 4.0 * g_z_xz_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_y_z[i] = -2.0 * g_z_xz_0_z[i] * a_exp + 4.0 * g_z_xz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (402-405)

    #pragma omp simd aligned(g_z_0_y_0_0_xz_z_x, g_z_0_y_0_0_xz_z_y, g_z_0_y_0_0_xz_z_z, g_z_xz_yz_x, g_z_xz_yz_y, g_z_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xz_z_x[i] = 4.0 * g_z_xz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_z_y[i] = 4.0 * g_z_xz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_z_z[i] = 4.0 * g_z_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (405-408)

    #pragma omp simd aligned(g_z_0_y_0_0_yy_x_x, g_z_0_y_0_0_yy_x_y, g_z_0_y_0_0_yy_x_z, g_z_yy_xy_x, g_z_yy_xy_y, g_z_yy_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yy_x_x[i] = 4.0 * g_z_yy_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_x_y[i] = 4.0 * g_z_yy_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_x_z[i] = 4.0 * g_z_yy_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (408-411)

    #pragma omp simd aligned(g_z_0_y_0_0_yy_y_x, g_z_0_y_0_0_yy_y_y, g_z_0_y_0_0_yy_y_z, g_z_yy_0_x, g_z_yy_0_y, g_z_yy_0_z, g_z_yy_yy_x, g_z_yy_yy_y, g_z_yy_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yy_y_x[i] = -2.0 * g_z_yy_0_x[i] * a_exp + 4.0 * g_z_yy_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_y_y[i] = -2.0 * g_z_yy_0_y[i] * a_exp + 4.0 * g_z_yy_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_y_z[i] = -2.0 * g_z_yy_0_z[i] * a_exp + 4.0 * g_z_yy_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (411-414)

    #pragma omp simd aligned(g_z_0_y_0_0_yy_z_x, g_z_0_y_0_0_yy_z_y, g_z_0_y_0_0_yy_z_z, g_z_yy_yz_x, g_z_yy_yz_y, g_z_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yy_z_x[i] = 4.0 * g_z_yy_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_z_y[i] = 4.0 * g_z_yy_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_z_z[i] = 4.0 * g_z_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (414-417)

    #pragma omp simd aligned(g_z_0_y_0_0_yz_x_x, g_z_0_y_0_0_yz_x_y, g_z_0_y_0_0_yz_x_z, g_z_yz_xy_x, g_z_yz_xy_y, g_z_yz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yz_x_x[i] = 4.0 * g_z_yz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_x_y[i] = 4.0 * g_z_yz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_x_z[i] = 4.0 * g_z_yz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (417-420)

    #pragma omp simd aligned(g_z_0_y_0_0_yz_y_x, g_z_0_y_0_0_yz_y_y, g_z_0_y_0_0_yz_y_z, g_z_yz_0_x, g_z_yz_0_y, g_z_yz_0_z, g_z_yz_yy_x, g_z_yz_yy_y, g_z_yz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yz_y_x[i] = -2.0 * g_z_yz_0_x[i] * a_exp + 4.0 * g_z_yz_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_y_y[i] = -2.0 * g_z_yz_0_y[i] * a_exp + 4.0 * g_z_yz_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_y_z[i] = -2.0 * g_z_yz_0_z[i] * a_exp + 4.0 * g_z_yz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (420-423)

    #pragma omp simd aligned(g_z_0_y_0_0_yz_z_x, g_z_0_y_0_0_yz_z_y, g_z_0_y_0_0_yz_z_z, g_z_yz_yz_x, g_z_yz_yz_y, g_z_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yz_z_x[i] = 4.0 * g_z_yz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_z_y[i] = 4.0 * g_z_yz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_z_z[i] = 4.0 * g_z_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (423-426)

    #pragma omp simd aligned(g_z_0_y_0_0_zz_x_x, g_z_0_y_0_0_zz_x_y, g_z_0_y_0_0_zz_x_z, g_z_zz_xy_x, g_z_zz_xy_y, g_z_zz_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_zz_x_x[i] = 4.0 * g_z_zz_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_x_y[i] = 4.0 * g_z_zz_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_x_z[i] = 4.0 * g_z_zz_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (426-429)

    #pragma omp simd aligned(g_z_0_y_0_0_zz_y_x, g_z_0_y_0_0_zz_y_y, g_z_0_y_0_0_zz_y_z, g_z_zz_0_x, g_z_zz_0_y, g_z_zz_0_z, g_z_zz_yy_x, g_z_zz_yy_y, g_z_zz_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_zz_y_x[i] = -2.0 * g_z_zz_0_x[i] * a_exp + 4.0 * g_z_zz_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_y_y[i] = -2.0 * g_z_zz_0_y[i] * a_exp + 4.0 * g_z_zz_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_y_z[i] = -2.0 * g_z_zz_0_z[i] * a_exp + 4.0 * g_z_zz_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (429-432)

    #pragma omp simd aligned(g_z_0_y_0_0_zz_z_x, g_z_0_y_0_0_zz_z_y, g_z_0_y_0_0_zz_z_z, g_z_zz_yz_x, g_z_zz_yz_y, g_z_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_zz_z_x[i] = 4.0 * g_z_zz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_z_y[i] = 4.0 * g_z_zz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_z_z[i] = 4.0 * g_z_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (432-435)

    #pragma omp simd aligned(g_z_0_z_0_0_xx_x_x, g_z_0_z_0_0_xx_x_y, g_z_0_z_0_0_xx_x_z, g_z_xx_xz_x, g_z_xx_xz_y, g_z_xx_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xx_x_x[i] = 4.0 * g_z_xx_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_x_y[i] = 4.0 * g_z_xx_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_x_z[i] = 4.0 * g_z_xx_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (435-438)

    #pragma omp simd aligned(g_z_0_z_0_0_xx_y_x, g_z_0_z_0_0_xx_y_y, g_z_0_z_0_0_xx_y_z, g_z_xx_yz_x, g_z_xx_yz_y, g_z_xx_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xx_y_x[i] = 4.0 * g_z_xx_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_y_y[i] = 4.0 * g_z_xx_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_y_z[i] = 4.0 * g_z_xx_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (438-441)

    #pragma omp simd aligned(g_z_0_z_0_0_xx_z_x, g_z_0_z_0_0_xx_z_y, g_z_0_z_0_0_xx_z_z, g_z_xx_0_x, g_z_xx_0_y, g_z_xx_0_z, g_z_xx_zz_x, g_z_xx_zz_y, g_z_xx_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xx_z_x[i] = -2.0 * g_z_xx_0_x[i] * a_exp + 4.0 * g_z_xx_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_z_y[i] = -2.0 * g_z_xx_0_y[i] * a_exp + 4.0 * g_z_xx_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_z_z[i] = -2.0 * g_z_xx_0_z[i] * a_exp + 4.0 * g_z_xx_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (441-444)

    #pragma omp simd aligned(g_z_0_z_0_0_xy_x_x, g_z_0_z_0_0_xy_x_y, g_z_0_z_0_0_xy_x_z, g_z_xy_xz_x, g_z_xy_xz_y, g_z_xy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xy_x_x[i] = 4.0 * g_z_xy_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_x_y[i] = 4.0 * g_z_xy_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_x_z[i] = 4.0 * g_z_xy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (444-447)

    #pragma omp simd aligned(g_z_0_z_0_0_xy_y_x, g_z_0_z_0_0_xy_y_y, g_z_0_z_0_0_xy_y_z, g_z_xy_yz_x, g_z_xy_yz_y, g_z_xy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xy_y_x[i] = 4.0 * g_z_xy_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_y_y[i] = 4.0 * g_z_xy_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_y_z[i] = 4.0 * g_z_xy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (447-450)

    #pragma omp simd aligned(g_z_0_z_0_0_xy_z_x, g_z_0_z_0_0_xy_z_y, g_z_0_z_0_0_xy_z_z, g_z_xy_0_x, g_z_xy_0_y, g_z_xy_0_z, g_z_xy_zz_x, g_z_xy_zz_y, g_z_xy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xy_z_x[i] = -2.0 * g_z_xy_0_x[i] * a_exp + 4.0 * g_z_xy_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_z_y[i] = -2.0 * g_z_xy_0_y[i] * a_exp + 4.0 * g_z_xy_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_z_z[i] = -2.0 * g_z_xy_0_z[i] * a_exp + 4.0 * g_z_xy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (450-453)

    #pragma omp simd aligned(g_z_0_z_0_0_xz_x_x, g_z_0_z_0_0_xz_x_y, g_z_0_z_0_0_xz_x_z, g_z_xz_xz_x, g_z_xz_xz_y, g_z_xz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xz_x_x[i] = 4.0 * g_z_xz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_x_y[i] = 4.0 * g_z_xz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_x_z[i] = 4.0 * g_z_xz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (453-456)

    #pragma omp simd aligned(g_z_0_z_0_0_xz_y_x, g_z_0_z_0_0_xz_y_y, g_z_0_z_0_0_xz_y_z, g_z_xz_yz_x, g_z_xz_yz_y, g_z_xz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xz_y_x[i] = 4.0 * g_z_xz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_y_y[i] = 4.0 * g_z_xz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_y_z[i] = 4.0 * g_z_xz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (456-459)

    #pragma omp simd aligned(g_z_0_z_0_0_xz_z_x, g_z_0_z_0_0_xz_z_y, g_z_0_z_0_0_xz_z_z, g_z_xz_0_x, g_z_xz_0_y, g_z_xz_0_z, g_z_xz_zz_x, g_z_xz_zz_y, g_z_xz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xz_z_x[i] = -2.0 * g_z_xz_0_x[i] * a_exp + 4.0 * g_z_xz_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_z_y[i] = -2.0 * g_z_xz_0_y[i] * a_exp + 4.0 * g_z_xz_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_z_z[i] = -2.0 * g_z_xz_0_z[i] * a_exp + 4.0 * g_z_xz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (459-462)

    #pragma omp simd aligned(g_z_0_z_0_0_yy_x_x, g_z_0_z_0_0_yy_x_y, g_z_0_z_0_0_yy_x_z, g_z_yy_xz_x, g_z_yy_xz_y, g_z_yy_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yy_x_x[i] = 4.0 * g_z_yy_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_x_y[i] = 4.0 * g_z_yy_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_x_z[i] = 4.0 * g_z_yy_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (462-465)

    #pragma omp simd aligned(g_z_0_z_0_0_yy_y_x, g_z_0_z_0_0_yy_y_y, g_z_0_z_0_0_yy_y_z, g_z_yy_yz_x, g_z_yy_yz_y, g_z_yy_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yy_y_x[i] = 4.0 * g_z_yy_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_y_y[i] = 4.0 * g_z_yy_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_y_z[i] = 4.0 * g_z_yy_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (465-468)

    #pragma omp simd aligned(g_z_0_z_0_0_yy_z_x, g_z_0_z_0_0_yy_z_y, g_z_0_z_0_0_yy_z_z, g_z_yy_0_x, g_z_yy_0_y, g_z_yy_0_z, g_z_yy_zz_x, g_z_yy_zz_y, g_z_yy_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yy_z_x[i] = -2.0 * g_z_yy_0_x[i] * a_exp + 4.0 * g_z_yy_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_z_y[i] = -2.0 * g_z_yy_0_y[i] * a_exp + 4.0 * g_z_yy_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_z_z[i] = -2.0 * g_z_yy_0_z[i] * a_exp + 4.0 * g_z_yy_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (468-471)

    #pragma omp simd aligned(g_z_0_z_0_0_yz_x_x, g_z_0_z_0_0_yz_x_y, g_z_0_z_0_0_yz_x_z, g_z_yz_xz_x, g_z_yz_xz_y, g_z_yz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yz_x_x[i] = 4.0 * g_z_yz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_x_y[i] = 4.0 * g_z_yz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_x_z[i] = 4.0 * g_z_yz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (471-474)

    #pragma omp simd aligned(g_z_0_z_0_0_yz_y_x, g_z_0_z_0_0_yz_y_y, g_z_0_z_0_0_yz_y_z, g_z_yz_yz_x, g_z_yz_yz_y, g_z_yz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yz_y_x[i] = 4.0 * g_z_yz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_y_y[i] = 4.0 * g_z_yz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_y_z[i] = 4.0 * g_z_yz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (474-477)

    #pragma omp simd aligned(g_z_0_z_0_0_yz_z_x, g_z_0_z_0_0_yz_z_y, g_z_0_z_0_0_yz_z_z, g_z_yz_0_x, g_z_yz_0_y, g_z_yz_0_z, g_z_yz_zz_x, g_z_yz_zz_y, g_z_yz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yz_z_x[i] = -2.0 * g_z_yz_0_x[i] * a_exp + 4.0 * g_z_yz_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_z_y[i] = -2.0 * g_z_yz_0_y[i] * a_exp + 4.0 * g_z_yz_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_z_z[i] = -2.0 * g_z_yz_0_z[i] * a_exp + 4.0 * g_z_yz_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (477-480)

    #pragma omp simd aligned(g_z_0_z_0_0_zz_x_x, g_z_0_z_0_0_zz_x_y, g_z_0_z_0_0_zz_x_z, g_z_zz_xz_x, g_z_zz_xz_y, g_z_zz_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_zz_x_x[i] = 4.0 * g_z_zz_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_x_y[i] = 4.0 * g_z_zz_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_x_z[i] = 4.0 * g_z_zz_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (480-483)

    #pragma omp simd aligned(g_z_0_z_0_0_zz_y_x, g_z_0_z_0_0_zz_y_y, g_z_0_z_0_0_zz_y_z, g_z_zz_yz_x, g_z_zz_yz_y, g_z_zz_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_zz_y_x[i] = 4.0 * g_z_zz_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_y_y[i] = 4.0 * g_z_zz_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_y_z[i] = 4.0 * g_z_zz_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (483-486)

    #pragma omp simd aligned(g_z_0_z_0_0_zz_z_x, g_z_0_z_0_0_zz_z_y, g_z_0_z_0_0_zz_z_z, g_z_zz_0_x, g_z_zz_0_y, g_z_zz_0_z, g_z_zz_zz_x, g_z_zz_zz_y, g_z_zz_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_zz_z_x[i] = -2.0 * g_z_zz_0_x[i] * a_exp + 4.0 * g_z_zz_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_z_y[i] = -2.0 * g_z_zz_0_y[i] * a_exp + 4.0 * g_z_zz_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_z_z[i] = -2.0 * g_z_zz_0_z[i] * a_exp + 4.0 * g_z_zz_zz_z[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace

