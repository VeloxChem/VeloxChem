#include "GeomDeriv1010OfScalarForSDSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_sdsd_0(CSimdArray<double>& buffer_1010_sdsd,
                     const CSimdArray<double>& buffer_pdpd,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_sdsd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pdpd

    auto g_x_xx_x_xx = buffer_pdpd[0];

    auto g_x_xx_x_xy = buffer_pdpd[1];

    auto g_x_xx_x_xz = buffer_pdpd[2];

    auto g_x_xx_x_yy = buffer_pdpd[3];

    auto g_x_xx_x_yz = buffer_pdpd[4];

    auto g_x_xx_x_zz = buffer_pdpd[5];

    auto g_x_xx_y_xx = buffer_pdpd[6];

    auto g_x_xx_y_xy = buffer_pdpd[7];

    auto g_x_xx_y_xz = buffer_pdpd[8];

    auto g_x_xx_y_yy = buffer_pdpd[9];

    auto g_x_xx_y_yz = buffer_pdpd[10];

    auto g_x_xx_y_zz = buffer_pdpd[11];

    auto g_x_xx_z_xx = buffer_pdpd[12];

    auto g_x_xx_z_xy = buffer_pdpd[13];

    auto g_x_xx_z_xz = buffer_pdpd[14];

    auto g_x_xx_z_yy = buffer_pdpd[15];

    auto g_x_xx_z_yz = buffer_pdpd[16];

    auto g_x_xx_z_zz = buffer_pdpd[17];

    auto g_x_xy_x_xx = buffer_pdpd[18];

    auto g_x_xy_x_xy = buffer_pdpd[19];

    auto g_x_xy_x_xz = buffer_pdpd[20];

    auto g_x_xy_x_yy = buffer_pdpd[21];

    auto g_x_xy_x_yz = buffer_pdpd[22];

    auto g_x_xy_x_zz = buffer_pdpd[23];

    auto g_x_xy_y_xx = buffer_pdpd[24];

    auto g_x_xy_y_xy = buffer_pdpd[25];

    auto g_x_xy_y_xz = buffer_pdpd[26];

    auto g_x_xy_y_yy = buffer_pdpd[27];

    auto g_x_xy_y_yz = buffer_pdpd[28];

    auto g_x_xy_y_zz = buffer_pdpd[29];

    auto g_x_xy_z_xx = buffer_pdpd[30];

    auto g_x_xy_z_xy = buffer_pdpd[31];

    auto g_x_xy_z_xz = buffer_pdpd[32];

    auto g_x_xy_z_yy = buffer_pdpd[33];

    auto g_x_xy_z_yz = buffer_pdpd[34];

    auto g_x_xy_z_zz = buffer_pdpd[35];

    auto g_x_xz_x_xx = buffer_pdpd[36];

    auto g_x_xz_x_xy = buffer_pdpd[37];

    auto g_x_xz_x_xz = buffer_pdpd[38];

    auto g_x_xz_x_yy = buffer_pdpd[39];

    auto g_x_xz_x_yz = buffer_pdpd[40];

    auto g_x_xz_x_zz = buffer_pdpd[41];

    auto g_x_xz_y_xx = buffer_pdpd[42];

    auto g_x_xz_y_xy = buffer_pdpd[43];

    auto g_x_xz_y_xz = buffer_pdpd[44];

    auto g_x_xz_y_yy = buffer_pdpd[45];

    auto g_x_xz_y_yz = buffer_pdpd[46];

    auto g_x_xz_y_zz = buffer_pdpd[47];

    auto g_x_xz_z_xx = buffer_pdpd[48];

    auto g_x_xz_z_xy = buffer_pdpd[49];

    auto g_x_xz_z_xz = buffer_pdpd[50];

    auto g_x_xz_z_yy = buffer_pdpd[51];

    auto g_x_xz_z_yz = buffer_pdpd[52];

    auto g_x_xz_z_zz = buffer_pdpd[53];

    auto g_x_yy_x_xx = buffer_pdpd[54];

    auto g_x_yy_x_xy = buffer_pdpd[55];

    auto g_x_yy_x_xz = buffer_pdpd[56];

    auto g_x_yy_x_yy = buffer_pdpd[57];

    auto g_x_yy_x_yz = buffer_pdpd[58];

    auto g_x_yy_x_zz = buffer_pdpd[59];

    auto g_x_yy_y_xx = buffer_pdpd[60];

    auto g_x_yy_y_xy = buffer_pdpd[61];

    auto g_x_yy_y_xz = buffer_pdpd[62];

    auto g_x_yy_y_yy = buffer_pdpd[63];

    auto g_x_yy_y_yz = buffer_pdpd[64];

    auto g_x_yy_y_zz = buffer_pdpd[65];

    auto g_x_yy_z_xx = buffer_pdpd[66];

    auto g_x_yy_z_xy = buffer_pdpd[67];

    auto g_x_yy_z_xz = buffer_pdpd[68];

    auto g_x_yy_z_yy = buffer_pdpd[69];

    auto g_x_yy_z_yz = buffer_pdpd[70];

    auto g_x_yy_z_zz = buffer_pdpd[71];

    auto g_x_yz_x_xx = buffer_pdpd[72];

    auto g_x_yz_x_xy = buffer_pdpd[73];

    auto g_x_yz_x_xz = buffer_pdpd[74];

    auto g_x_yz_x_yy = buffer_pdpd[75];

    auto g_x_yz_x_yz = buffer_pdpd[76];

    auto g_x_yz_x_zz = buffer_pdpd[77];

    auto g_x_yz_y_xx = buffer_pdpd[78];

    auto g_x_yz_y_xy = buffer_pdpd[79];

    auto g_x_yz_y_xz = buffer_pdpd[80];

    auto g_x_yz_y_yy = buffer_pdpd[81];

    auto g_x_yz_y_yz = buffer_pdpd[82];

    auto g_x_yz_y_zz = buffer_pdpd[83];

    auto g_x_yz_z_xx = buffer_pdpd[84];

    auto g_x_yz_z_xy = buffer_pdpd[85];

    auto g_x_yz_z_xz = buffer_pdpd[86];

    auto g_x_yz_z_yy = buffer_pdpd[87];

    auto g_x_yz_z_yz = buffer_pdpd[88];

    auto g_x_yz_z_zz = buffer_pdpd[89];

    auto g_x_zz_x_xx = buffer_pdpd[90];

    auto g_x_zz_x_xy = buffer_pdpd[91];

    auto g_x_zz_x_xz = buffer_pdpd[92];

    auto g_x_zz_x_yy = buffer_pdpd[93];

    auto g_x_zz_x_yz = buffer_pdpd[94];

    auto g_x_zz_x_zz = buffer_pdpd[95];

    auto g_x_zz_y_xx = buffer_pdpd[96];

    auto g_x_zz_y_xy = buffer_pdpd[97];

    auto g_x_zz_y_xz = buffer_pdpd[98];

    auto g_x_zz_y_yy = buffer_pdpd[99];

    auto g_x_zz_y_yz = buffer_pdpd[100];

    auto g_x_zz_y_zz = buffer_pdpd[101];

    auto g_x_zz_z_xx = buffer_pdpd[102];

    auto g_x_zz_z_xy = buffer_pdpd[103];

    auto g_x_zz_z_xz = buffer_pdpd[104];

    auto g_x_zz_z_yy = buffer_pdpd[105];

    auto g_x_zz_z_yz = buffer_pdpd[106];

    auto g_x_zz_z_zz = buffer_pdpd[107];

    auto g_y_xx_x_xx = buffer_pdpd[108];

    auto g_y_xx_x_xy = buffer_pdpd[109];

    auto g_y_xx_x_xz = buffer_pdpd[110];

    auto g_y_xx_x_yy = buffer_pdpd[111];

    auto g_y_xx_x_yz = buffer_pdpd[112];

    auto g_y_xx_x_zz = buffer_pdpd[113];

    auto g_y_xx_y_xx = buffer_pdpd[114];

    auto g_y_xx_y_xy = buffer_pdpd[115];

    auto g_y_xx_y_xz = buffer_pdpd[116];

    auto g_y_xx_y_yy = buffer_pdpd[117];

    auto g_y_xx_y_yz = buffer_pdpd[118];

    auto g_y_xx_y_zz = buffer_pdpd[119];

    auto g_y_xx_z_xx = buffer_pdpd[120];

    auto g_y_xx_z_xy = buffer_pdpd[121];

    auto g_y_xx_z_xz = buffer_pdpd[122];

    auto g_y_xx_z_yy = buffer_pdpd[123];

    auto g_y_xx_z_yz = buffer_pdpd[124];

    auto g_y_xx_z_zz = buffer_pdpd[125];

    auto g_y_xy_x_xx = buffer_pdpd[126];

    auto g_y_xy_x_xy = buffer_pdpd[127];

    auto g_y_xy_x_xz = buffer_pdpd[128];

    auto g_y_xy_x_yy = buffer_pdpd[129];

    auto g_y_xy_x_yz = buffer_pdpd[130];

    auto g_y_xy_x_zz = buffer_pdpd[131];

    auto g_y_xy_y_xx = buffer_pdpd[132];

    auto g_y_xy_y_xy = buffer_pdpd[133];

    auto g_y_xy_y_xz = buffer_pdpd[134];

    auto g_y_xy_y_yy = buffer_pdpd[135];

    auto g_y_xy_y_yz = buffer_pdpd[136];

    auto g_y_xy_y_zz = buffer_pdpd[137];

    auto g_y_xy_z_xx = buffer_pdpd[138];

    auto g_y_xy_z_xy = buffer_pdpd[139];

    auto g_y_xy_z_xz = buffer_pdpd[140];

    auto g_y_xy_z_yy = buffer_pdpd[141];

    auto g_y_xy_z_yz = buffer_pdpd[142];

    auto g_y_xy_z_zz = buffer_pdpd[143];

    auto g_y_xz_x_xx = buffer_pdpd[144];

    auto g_y_xz_x_xy = buffer_pdpd[145];

    auto g_y_xz_x_xz = buffer_pdpd[146];

    auto g_y_xz_x_yy = buffer_pdpd[147];

    auto g_y_xz_x_yz = buffer_pdpd[148];

    auto g_y_xz_x_zz = buffer_pdpd[149];

    auto g_y_xz_y_xx = buffer_pdpd[150];

    auto g_y_xz_y_xy = buffer_pdpd[151];

    auto g_y_xz_y_xz = buffer_pdpd[152];

    auto g_y_xz_y_yy = buffer_pdpd[153];

    auto g_y_xz_y_yz = buffer_pdpd[154];

    auto g_y_xz_y_zz = buffer_pdpd[155];

    auto g_y_xz_z_xx = buffer_pdpd[156];

    auto g_y_xz_z_xy = buffer_pdpd[157];

    auto g_y_xz_z_xz = buffer_pdpd[158];

    auto g_y_xz_z_yy = buffer_pdpd[159];

    auto g_y_xz_z_yz = buffer_pdpd[160];

    auto g_y_xz_z_zz = buffer_pdpd[161];

    auto g_y_yy_x_xx = buffer_pdpd[162];

    auto g_y_yy_x_xy = buffer_pdpd[163];

    auto g_y_yy_x_xz = buffer_pdpd[164];

    auto g_y_yy_x_yy = buffer_pdpd[165];

    auto g_y_yy_x_yz = buffer_pdpd[166];

    auto g_y_yy_x_zz = buffer_pdpd[167];

    auto g_y_yy_y_xx = buffer_pdpd[168];

    auto g_y_yy_y_xy = buffer_pdpd[169];

    auto g_y_yy_y_xz = buffer_pdpd[170];

    auto g_y_yy_y_yy = buffer_pdpd[171];

    auto g_y_yy_y_yz = buffer_pdpd[172];

    auto g_y_yy_y_zz = buffer_pdpd[173];

    auto g_y_yy_z_xx = buffer_pdpd[174];

    auto g_y_yy_z_xy = buffer_pdpd[175];

    auto g_y_yy_z_xz = buffer_pdpd[176];

    auto g_y_yy_z_yy = buffer_pdpd[177];

    auto g_y_yy_z_yz = buffer_pdpd[178];

    auto g_y_yy_z_zz = buffer_pdpd[179];

    auto g_y_yz_x_xx = buffer_pdpd[180];

    auto g_y_yz_x_xy = buffer_pdpd[181];

    auto g_y_yz_x_xz = buffer_pdpd[182];

    auto g_y_yz_x_yy = buffer_pdpd[183];

    auto g_y_yz_x_yz = buffer_pdpd[184];

    auto g_y_yz_x_zz = buffer_pdpd[185];

    auto g_y_yz_y_xx = buffer_pdpd[186];

    auto g_y_yz_y_xy = buffer_pdpd[187];

    auto g_y_yz_y_xz = buffer_pdpd[188];

    auto g_y_yz_y_yy = buffer_pdpd[189];

    auto g_y_yz_y_yz = buffer_pdpd[190];

    auto g_y_yz_y_zz = buffer_pdpd[191];

    auto g_y_yz_z_xx = buffer_pdpd[192];

    auto g_y_yz_z_xy = buffer_pdpd[193];

    auto g_y_yz_z_xz = buffer_pdpd[194];

    auto g_y_yz_z_yy = buffer_pdpd[195];

    auto g_y_yz_z_yz = buffer_pdpd[196];

    auto g_y_yz_z_zz = buffer_pdpd[197];

    auto g_y_zz_x_xx = buffer_pdpd[198];

    auto g_y_zz_x_xy = buffer_pdpd[199];

    auto g_y_zz_x_xz = buffer_pdpd[200];

    auto g_y_zz_x_yy = buffer_pdpd[201];

    auto g_y_zz_x_yz = buffer_pdpd[202];

    auto g_y_zz_x_zz = buffer_pdpd[203];

    auto g_y_zz_y_xx = buffer_pdpd[204];

    auto g_y_zz_y_xy = buffer_pdpd[205];

    auto g_y_zz_y_xz = buffer_pdpd[206];

    auto g_y_zz_y_yy = buffer_pdpd[207];

    auto g_y_zz_y_yz = buffer_pdpd[208];

    auto g_y_zz_y_zz = buffer_pdpd[209];

    auto g_y_zz_z_xx = buffer_pdpd[210];

    auto g_y_zz_z_xy = buffer_pdpd[211];

    auto g_y_zz_z_xz = buffer_pdpd[212];

    auto g_y_zz_z_yy = buffer_pdpd[213];

    auto g_y_zz_z_yz = buffer_pdpd[214];

    auto g_y_zz_z_zz = buffer_pdpd[215];

    auto g_z_xx_x_xx = buffer_pdpd[216];

    auto g_z_xx_x_xy = buffer_pdpd[217];

    auto g_z_xx_x_xz = buffer_pdpd[218];

    auto g_z_xx_x_yy = buffer_pdpd[219];

    auto g_z_xx_x_yz = buffer_pdpd[220];

    auto g_z_xx_x_zz = buffer_pdpd[221];

    auto g_z_xx_y_xx = buffer_pdpd[222];

    auto g_z_xx_y_xy = buffer_pdpd[223];

    auto g_z_xx_y_xz = buffer_pdpd[224];

    auto g_z_xx_y_yy = buffer_pdpd[225];

    auto g_z_xx_y_yz = buffer_pdpd[226];

    auto g_z_xx_y_zz = buffer_pdpd[227];

    auto g_z_xx_z_xx = buffer_pdpd[228];

    auto g_z_xx_z_xy = buffer_pdpd[229];

    auto g_z_xx_z_xz = buffer_pdpd[230];

    auto g_z_xx_z_yy = buffer_pdpd[231];

    auto g_z_xx_z_yz = buffer_pdpd[232];

    auto g_z_xx_z_zz = buffer_pdpd[233];

    auto g_z_xy_x_xx = buffer_pdpd[234];

    auto g_z_xy_x_xy = buffer_pdpd[235];

    auto g_z_xy_x_xz = buffer_pdpd[236];

    auto g_z_xy_x_yy = buffer_pdpd[237];

    auto g_z_xy_x_yz = buffer_pdpd[238];

    auto g_z_xy_x_zz = buffer_pdpd[239];

    auto g_z_xy_y_xx = buffer_pdpd[240];

    auto g_z_xy_y_xy = buffer_pdpd[241];

    auto g_z_xy_y_xz = buffer_pdpd[242];

    auto g_z_xy_y_yy = buffer_pdpd[243];

    auto g_z_xy_y_yz = buffer_pdpd[244];

    auto g_z_xy_y_zz = buffer_pdpd[245];

    auto g_z_xy_z_xx = buffer_pdpd[246];

    auto g_z_xy_z_xy = buffer_pdpd[247];

    auto g_z_xy_z_xz = buffer_pdpd[248];

    auto g_z_xy_z_yy = buffer_pdpd[249];

    auto g_z_xy_z_yz = buffer_pdpd[250];

    auto g_z_xy_z_zz = buffer_pdpd[251];

    auto g_z_xz_x_xx = buffer_pdpd[252];

    auto g_z_xz_x_xy = buffer_pdpd[253];

    auto g_z_xz_x_xz = buffer_pdpd[254];

    auto g_z_xz_x_yy = buffer_pdpd[255];

    auto g_z_xz_x_yz = buffer_pdpd[256];

    auto g_z_xz_x_zz = buffer_pdpd[257];

    auto g_z_xz_y_xx = buffer_pdpd[258];

    auto g_z_xz_y_xy = buffer_pdpd[259];

    auto g_z_xz_y_xz = buffer_pdpd[260];

    auto g_z_xz_y_yy = buffer_pdpd[261];

    auto g_z_xz_y_yz = buffer_pdpd[262];

    auto g_z_xz_y_zz = buffer_pdpd[263];

    auto g_z_xz_z_xx = buffer_pdpd[264];

    auto g_z_xz_z_xy = buffer_pdpd[265];

    auto g_z_xz_z_xz = buffer_pdpd[266];

    auto g_z_xz_z_yy = buffer_pdpd[267];

    auto g_z_xz_z_yz = buffer_pdpd[268];

    auto g_z_xz_z_zz = buffer_pdpd[269];

    auto g_z_yy_x_xx = buffer_pdpd[270];

    auto g_z_yy_x_xy = buffer_pdpd[271];

    auto g_z_yy_x_xz = buffer_pdpd[272];

    auto g_z_yy_x_yy = buffer_pdpd[273];

    auto g_z_yy_x_yz = buffer_pdpd[274];

    auto g_z_yy_x_zz = buffer_pdpd[275];

    auto g_z_yy_y_xx = buffer_pdpd[276];

    auto g_z_yy_y_xy = buffer_pdpd[277];

    auto g_z_yy_y_xz = buffer_pdpd[278];

    auto g_z_yy_y_yy = buffer_pdpd[279];

    auto g_z_yy_y_yz = buffer_pdpd[280];

    auto g_z_yy_y_zz = buffer_pdpd[281];

    auto g_z_yy_z_xx = buffer_pdpd[282];

    auto g_z_yy_z_xy = buffer_pdpd[283];

    auto g_z_yy_z_xz = buffer_pdpd[284];

    auto g_z_yy_z_yy = buffer_pdpd[285];

    auto g_z_yy_z_yz = buffer_pdpd[286];

    auto g_z_yy_z_zz = buffer_pdpd[287];

    auto g_z_yz_x_xx = buffer_pdpd[288];

    auto g_z_yz_x_xy = buffer_pdpd[289];

    auto g_z_yz_x_xz = buffer_pdpd[290];

    auto g_z_yz_x_yy = buffer_pdpd[291];

    auto g_z_yz_x_yz = buffer_pdpd[292];

    auto g_z_yz_x_zz = buffer_pdpd[293];

    auto g_z_yz_y_xx = buffer_pdpd[294];

    auto g_z_yz_y_xy = buffer_pdpd[295];

    auto g_z_yz_y_xz = buffer_pdpd[296];

    auto g_z_yz_y_yy = buffer_pdpd[297];

    auto g_z_yz_y_yz = buffer_pdpd[298];

    auto g_z_yz_y_zz = buffer_pdpd[299];

    auto g_z_yz_z_xx = buffer_pdpd[300];

    auto g_z_yz_z_xy = buffer_pdpd[301];

    auto g_z_yz_z_xz = buffer_pdpd[302];

    auto g_z_yz_z_yy = buffer_pdpd[303];

    auto g_z_yz_z_yz = buffer_pdpd[304];

    auto g_z_yz_z_zz = buffer_pdpd[305];

    auto g_z_zz_x_xx = buffer_pdpd[306];

    auto g_z_zz_x_xy = buffer_pdpd[307];

    auto g_z_zz_x_xz = buffer_pdpd[308];

    auto g_z_zz_x_yy = buffer_pdpd[309];

    auto g_z_zz_x_yz = buffer_pdpd[310];

    auto g_z_zz_x_zz = buffer_pdpd[311];

    auto g_z_zz_y_xx = buffer_pdpd[312];

    auto g_z_zz_y_xy = buffer_pdpd[313];

    auto g_z_zz_y_xz = buffer_pdpd[314];

    auto g_z_zz_y_yy = buffer_pdpd[315];

    auto g_z_zz_y_yz = buffer_pdpd[316];

    auto g_z_zz_y_zz = buffer_pdpd[317];

    auto g_z_zz_z_xx = buffer_pdpd[318];

    auto g_z_zz_z_xy = buffer_pdpd[319];

    auto g_z_zz_z_xz = buffer_pdpd[320];

    auto g_z_zz_z_yy = buffer_pdpd[321];

    auto g_z_zz_z_yz = buffer_pdpd[322];

    auto g_z_zz_z_zz = buffer_pdpd[323];

    /// Set up components of integrals buffer : buffer_1010_sdsd

    auto g_x_0_x_0_0_xx_0_xx = buffer_1010_sdsd[0];

    auto g_x_0_x_0_0_xx_0_xy = buffer_1010_sdsd[1];

    auto g_x_0_x_0_0_xx_0_xz = buffer_1010_sdsd[2];

    auto g_x_0_x_0_0_xx_0_yy = buffer_1010_sdsd[3];

    auto g_x_0_x_0_0_xx_0_yz = buffer_1010_sdsd[4];

    auto g_x_0_x_0_0_xx_0_zz = buffer_1010_sdsd[5];

    auto g_x_0_x_0_0_xy_0_xx = buffer_1010_sdsd[6];

    auto g_x_0_x_0_0_xy_0_xy = buffer_1010_sdsd[7];

    auto g_x_0_x_0_0_xy_0_xz = buffer_1010_sdsd[8];

    auto g_x_0_x_0_0_xy_0_yy = buffer_1010_sdsd[9];

    auto g_x_0_x_0_0_xy_0_yz = buffer_1010_sdsd[10];

    auto g_x_0_x_0_0_xy_0_zz = buffer_1010_sdsd[11];

    auto g_x_0_x_0_0_xz_0_xx = buffer_1010_sdsd[12];

    auto g_x_0_x_0_0_xz_0_xy = buffer_1010_sdsd[13];

    auto g_x_0_x_0_0_xz_0_xz = buffer_1010_sdsd[14];

    auto g_x_0_x_0_0_xz_0_yy = buffer_1010_sdsd[15];

    auto g_x_0_x_0_0_xz_0_yz = buffer_1010_sdsd[16];

    auto g_x_0_x_0_0_xz_0_zz = buffer_1010_sdsd[17];

    auto g_x_0_x_0_0_yy_0_xx = buffer_1010_sdsd[18];

    auto g_x_0_x_0_0_yy_0_xy = buffer_1010_sdsd[19];

    auto g_x_0_x_0_0_yy_0_xz = buffer_1010_sdsd[20];

    auto g_x_0_x_0_0_yy_0_yy = buffer_1010_sdsd[21];

    auto g_x_0_x_0_0_yy_0_yz = buffer_1010_sdsd[22];

    auto g_x_0_x_0_0_yy_0_zz = buffer_1010_sdsd[23];

    auto g_x_0_x_0_0_yz_0_xx = buffer_1010_sdsd[24];

    auto g_x_0_x_0_0_yz_0_xy = buffer_1010_sdsd[25];

    auto g_x_0_x_0_0_yz_0_xz = buffer_1010_sdsd[26];

    auto g_x_0_x_0_0_yz_0_yy = buffer_1010_sdsd[27];

    auto g_x_0_x_0_0_yz_0_yz = buffer_1010_sdsd[28];

    auto g_x_0_x_0_0_yz_0_zz = buffer_1010_sdsd[29];

    auto g_x_0_x_0_0_zz_0_xx = buffer_1010_sdsd[30];

    auto g_x_0_x_0_0_zz_0_xy = buffer_1010_sdsd[31];

    auto g_x_0_x_0_0_zz_0_xz = buffer_1010_sdsd[32];

    auto g_x_0_x_0_0_zz_0_yy = buffer_1010_sdsd[33];

    auto g_x_0_x_0_0_zz_0_yz = buffer_1010_sdsd[34];

    auto g_x_0_x_0_0_zz_0_zz = buffer_1010_sdsd[35];

    auto g_x_0_y_0_0_xx_0_xx = buffer_1010_sdsd[36];

    auto g_x_0_y_0_0_xx_0_xy = buffer_1010_sdsd[37];

    auto g_x_0_y_0_0_xx_0_xz = buffer_1010_sdsd[38];

    auto g_x_0_y_0_0_xx_0_yy = buffer_1010_sdsd[39];

    auto g_x_0_y_0_0_xx_0_yz = buffer_1010_sdsd[40];

    auto g_x_0_y_0_0_xx_0_zz = buffer_1010_sdsd[41];

    auto g_x_0_y_0_0_xy_0_xx = buffer_1010_sdsd[42];

    auto g_x_0_y_0_0_xy_0_xy = buffer_1010_sdsd[43];

    auto g_x_0_y_0_0_xy_0_xz = buffer_1010_sdsd[44];

    auto g_x_0_y_0_0_xy_0_yy = buffer_1010_sdsd[45];

    auto g_x_0_y_0_0_xy_0_yz = buffer_1010_sdsd[46];

    auto g_x_0_y_0_0_xy_0_zz = buffer_1010_sdsd[47];

    auto g_x_0_y_0_0_xz_0_xx = buffer_1010_sdsd[48];

    auto g_x_0_y_0_0_xz_0_xy = buffer_1010_sdsd[49];

    auto g_x_0_y_0_0_xz_0_xz = buffer_1010_sdsd[50];

    auto g_x_0_y_0_0_xz_0_yy = buffer_1010_sdsd[51];

    auto g_x_0_y_0_0_xz_0_yz = buffer_1010_sdsd[52];

    auto g_x_0_y_0_0_xz_0_zz = buffer_1010_sdsd[53];

    auto g_x_0_y_0_0_yy_0_xx = buffer_1010_sdsd[54];

    auto g_x_0_y_0_0_yy_0_xy = buffer_1010_sdsd[55];

    auto g_x_0_y_0_0_yy_0_xz = buffer_1010_sdsd[56];

    auto g_x_0_y_0_0_yy_0_yy = buffer_1010_sdsd[57];

    auto g_x_0_y_0_0_yy_0_yz = buffer_1010_sdsd[58];

    auto g_x_0_y_0_0_yy_0_zz = buffer_1010_sdsd[59];

    auto g_x_0_y_0_0_yz_0_xx = buffer_1010_sdsd[60];

    auto g_x_0_y_0_0_yz_0_xy = buffer_1010_sdsd[61];

    auto g_x_0_y_0_0_yz_0_xz = buffer_1010_sdsd[62];

    auto g_x_0_y_0_0_yz_0_yy = buffer_1010_sdsd[63];

    auto g_x_0_y_0_0_yz_0_yz = buffer_1010_sdsd[64];

    auto g_x_0_y_0_0_yz_0_zz = buffer_1010_sdsd[65];

    auto g_x_0_y_0_0_zz_0_xx = buffer_1010_sdsd[66];

    auto g_x_0_y_0_0_zz_0_xy = buffer_1010_sdsd[67];

    auto g_x_0_y_0_0_zz_0_xz = buffer_1010_sdsd[68];

    auto g_x_0_y_0_0_zz_0_yy = buffer_1010_sdsd[69];

    auto g_x_0_y_0_0_zz_0_yz = buffer_1010_sdsd[70];

    auto g_x_0_y_0_0_zz_0_zz = buffer_1010_sdsd[71];

    auto g_x_0_z_0_0_xx_0_xx = buffer_1010_sdsd[72];

    auto g_x_0_z_0_0_xx_0_xy = buffer_1010_sdsd[73];

    auto g_x_0_z_0_0_xx_0_xz = buffer_1010_sdsd[74];

    auto g_x_0_z_0_0_xx_0_yy = buffer_1010_sdsd[75];

    auto g_x_0_z_0_0_xx_0_yz = buffer_1010_sdsd[76];

    auto g_x_0_z_0_0_xx_0_zz = buffer_1010_sdsd[77];

    auto g_x_0_z_0_0_xy_0_xx = buffer_1010_sdsd[78];

    auto g_x_0_z_0_0_xy_0_xy = buffer_1010_sdsd[79];

    auto g_x_0_z_0_0_xy_0_xz = buffer_1010_sdsd[80];

    auto g_x_0_z_0_0_xy_0_yy = buffer_1010_sdsd[81];

    auto g_x_0_z_0_0_xy_0_yz = buffer_1010_sdsd[82];

    auto g_x_0_z_0_0_xy_0_zz = buffer_1010_sdsd[83];

    auto g_x_0_z_0_0_xz_0_xx = buffer_1010_sdsd[84];

    auto g_x_0_z_0_0_xz_0_xy = buffer_1010_sdsd[85];

    auto g_x_0_z_0_0_xz_0_xz = buffer_1010_sdsd[86];

    auto g_x_0_z_0_0_xz_0_yy = buffer_1010_sdsd[87];

    auto g_x_0_z_0_0_xz_0_yz = buffer_1010_sdsd[88];

    auto g_x_0_z_0_0_xz_0_zz = buffer_1010_sdsd[89];

    auto g_x_0_z_0_0_yy_0_xx = buffer_1010_sdsd[90];

    auto g_x_0_z_0_0_yy_0_xy = buffer_1010_sdsd[91];

    auto g_x_0_z_0_0_yy_0_xz = buffer_1010_sdsd[92];

    auto g_x_0_z_0_0_yy_0_yy = buffer_1010_sdsd[93];

    auto g_x_0_z_0_0_yy_0_yz = buffer_1010_sdsd[94];

    auto g_x_0_z_0_0_yy_0_zz = buffer_1010_sdsd[95];

    auto g_x_0_z_0_0_yz_0_xx = buffer_1010_sdsd[96];

    auto g_x_0_z_0_0_yz_0_xy = buffer_1010_sdsd[97];

    auto g_x_0_z_0_0_yz_0_xz = buffer_1010_sdsd[98];

    auto g_x_0_z_0_0_yz_0_yy = buffer_1010_sdsd[99];

    auto g_x_0_z_0_0_yz_0_yz = buffer_1010_sdsd[100];

    auto g_x_0_z_0_0_yz_0_zz = buffer_1010_sdsd[101];

    auto g_x_0_z_0_0_zz_0_xx = buffer_1010_sdsd[102];

    auto g_x_0_z_0_0_zz_0_xy = buffer_1010_sdsd[103];

    auto g_x_0_z_0_0_zz_0_xz = buffer_1010_sdsd[104];

    auto g_x_0_z_0_0_zz_0_yy = buffer_1010_sdsd[105];

    auto g_x_0_z_0_0_zz_0_yz = buffer_1010_sdsd[106];

    auto g_x_0_z_0_0_zz_0_zz = buffer_1010_sdsd[107];

    auto g_y_0_x_0_0_xx_0_xx = buffer_1010_sdsd[108];

    auto g_y_0_x_0_0_xx_0_xy = buffer_1010_sdsd[109];

    auto g_y_0_x_0_0_xx_0_xz = buffer_1010_sdsd[110];

    auto g_y_0_x_0_0_xx_0_yy = buffer_1010_sdsd[111];

    auto g_y_0_x_0_0_xx_0_yz = buffer_1010_sdsd[112];

    auto g_y_0_x_0_0_xx_0_zz = buffer_1010_sdsd[113];

    auto g_y_0_x_0_0_xy_0_xx = buffer_1010_sdsd[114];

    auto g_y_0_x_0_0_xy_0_xy = buffer_1010_sdsd[115];

    auto g_y_0_x_0_0_xy_0_xz = buffer_1010_sdsd[116];

    auto g_y_0_x_0_0_xy_0_yy = buffer_1010_sdsd[117];

    auto g_y_0_x_0_0_xy_0_yz = buffer_1010_sdsd[118];

    auto g_y_0_x_0_0_xy_0_zz = buffer_1010_sdsd[119];

    auto g_y_0_x_0_0_xz_0_xx = buffer_1010_sdsd[120];

    auto g_y_0_x_0_0_xz_0_xy = buffer_1010_sdsd[121];

    auto g_y_0_x_0_0_xz_0_xz = buffer_1010_sdsd[122];

    auto g_y_0_x_0_0_xz_0_yy = buffer_1010_sdsd[123];

    auto g_y_0_x_0_0_xz_0_yz = buffer_1010_sdsd[124];

    auto g_y_0_x_0_0_xz_0_zz = buffer_1010_sdsd[125];

    auto g_y_0_x_0_0_yy_0_xx = buffer_1010_sdsd[126];

    auto g_y_0_x_0_0_yy_0_xy = buffer_1010_sdsd[127];

    auto g_y_0_x_0_0_yy_0_xz = buffer_1010_sdsd[128];

    auto g_y_0_x_0_0_yy_0_yy = buffer_1010_sdsd[129];

    auto g_y_0_x_0_0_yy_0_yz = buffer_1010_sdsd[130];

    auto g_y_0_x_0_0_yy_0_zz = buffer_1010_sdsd[131];

    auto g_y_0_x_0_0_yz_0_xx = buffer_1010_sdsd[132];

    auto g_y_0_x_0_0_yz_0_xy = buffer_1010_sdsd[133];

    auto g_y_0_x_0_0_yz_0_xz = buffer_1010_sdsd[134];

    auto g_y_0_x_0_0_yz_0_yy = buffer_1010_sdsd[135];

    auto g_y_0_x_0_0_yz_0_yz = buffer_1010_sdsd[136];

    auto g_y_0_x_0_0_yz_0_zz = buffer_1010_sdsd[137];

    auto g_y_0_x_0_0_zz_0_xx = buffer_1010_sdsd[138];

    auto g_y_0_x_0_0_zz_0_xy = buffer_1010_sdsd[139];

    auto g_y_0_x_0_0_zz_0_xz = buffer_1010_sdsd[140];

    auto g_y_0_x_0_0_zz_0_yy = buffer_1010_sdsd[141];

    auto g_y_0_x_0_0_zz_0_yz = buffer_1010_sdsd[142];

    auto g_y_0_x_0_0_zz_0_zz = buffer_1010_sdsd[143];

    auto g_y_0_y_0_0_xx_0_xx = buffer_1010_sdsd[144];

    auto g_y_0_y_0_0_xx_0_xy = buffer_1010_sdsd[145];

    auto g_y_0_y_0_0_xx_0_xz = buffer_1010_sdsd[146];

    auto g_y_0_y_0_0_xx_0_yy = buffer_1010_sdsd[147];

    auto g_y_0_y_0_0_xx_0_yz = buffer_1010_sdsd[148];

    auto g_y_0_y_0_0_xx_0_zz = buffer_1010_sdsd[149];

    auto g_y_0_y_0_0_xy_0_xx = buffer_1010_sdsd[150];

    auto g_y_0_y_0_0_xy_0_xy = buffer_1010_sdsd[151];

    auto g_y_0_y_0_0_xy_0_xz = buffer_1010_sdsd[152];

    auto g_y_0_y_0_0_xy_0_yy = buffer_1010_sdsd[153];

    auto g_y_0_y_0_0_xy_0_yz = buffer_1010_sdsd[154];

    auto g_y_0_y_0_0_xy_0_zz = buffer_1010_sdsd[155];

    auto g_y_0_y_0_0_xz_0_xx = buffer_1010_sdsd[156];

    auto g_y_0_y_0_0_xz_0_xy = buffer_1010_sdsd[157];

    auto g_y_0_y_0_0_xz_0_xz = buffer_1010_sdsd[158];

    auto g_y_0_y_0_0_xz_0_yy = buffer_1010_sdsd[159];

    auto g_y_0_y_0_0_xz_0_yz = buffer_1010_sdsd[160];

    auto g_y_0_y_0_0_xz_0_zz = buffer_1010_sdsd[161];

    auto g_y_0_y_0_0_yy_0_xx = buffer_1010_sdsd[162];

    auto g_y_0_y_0_0_yy_0_xy = buffer_1010_sdsd[163];

    auto g_y_0_y_0_0_yy_0_xz = buffer_1010_sdsd[164];

    auto g_y_0_y_0_0_yy_0_yy = buffer_1010_sdsd[165];

    auto g_y_0_y_0_0_yy_0_yz = buffer_1010_sdsd[166];

    auto g_y_0_y_0_0_yy_0_zz = buffer_1010_sdsd[167];

    auto g_y_0_y_0_0_yz_0_xx = buffer_1010_sdsd[168];

    auto g_y_0_y_0_0_yz_0_xy = buffer_1010_sdsd[169];

    auto g_y_0_y_0_0_yz_0_xz = buffer_1010_sdsd[170];

    auto g_y_0_y_0_0_yz_0_yy = buffer_1010_sdsd[171];

    auto g_y_0_y_0_0_yz_0_yz = buffer_1010_sdsd[172];

    auto g_y_0_y_0_0_yz_0_zz = buffer_1010_sdsd[173];

    auto g_y_0_y_0_0_zz_0_xx = buffer_1010_sdsd[174];

    auto g_y_0_y_0_0_zz_0_xy = buffer_1010_sdsd[175];

    auto g_y_0_y_0_0_zz_0_xz = buffer_1010_sdsd[176];

    auto g_y_0_y_0_0_zz_0_yy = buffer_1010_sdsd[177];

    auto g_y_0_y_0_0_zz_0_yz = buffer_1010_sdsd[178];

    auto g_y_0_y_0_0_zz_0_zz = buffer_1010_sdsd[179];

    auto g_y_0_z_0_0_xx_0_xx = buffer_1010_sdsd[180];

    auto g_y_0_z_0_0_xx_0_xy = buffer_1010_sdsd[181];

    auto g_y_0_z_0_0_xx_0_xz = buffer_1010_sdsd[182];

    auto g_y_0_z_0_0_xx_0_yy = buffer_1010_sdsd[183];

    auto g_y_0_z_0_0_xx_0_yz = buffer_1010_sdsd[184];

    auto g_y_0_z_0_0_xx_0_zz = buffer_1010_sdsd[185];

    auto g_y_0_z_0_0_xy_0_xx = buffer_1010_sdsd[186];

    auto g_y_0_z_0_0_xy_0_xy = buffer_1010_sdsd[187];

    auto g_y_0_z_0_0_xy_0_xz = buffer_1010_sdsd[188];

    auto g_y_0_z_0_0_xy_0_yy = buffer_1010_sdsd[189];

    auto g_y_0_z_0_0_xy_0_yz = buffer_1010_sdsd[190];

    auto g_y_0_z_0_0_xy_0_zz = buffer_1010_sdsd[191];

    auto g_y_0_z_0_0_xz_0_xx = buffer_1010_sdsd[192];

    auto g_y_0_z_0_0_xz_0_xy = buffer_1010_sdsd[193];

    auto g_y_0_z_0_0_xz_0_xz = buffer_1010_sdsd[194];

    auto g_y_0_z_0_0_xz_0_yy = buffer_1010_sdsd[195];

    auto g_y_0_z_0_0_xz_0_yz = buffer_1010_sdsd[196];

    auto g_y_0_z_0_0_xz_0_zz = buffer_1010_sdsd[197];

    auto g_y_0_z_0_0_yy_0_xx = buffer_1010_sdsd[198];

    auto g_y_0_z_0_0_yy_0_xy = buffer_1010_sdsd[199];

    auto g_y_0_z_0_0_yy_0_xz = buffer_1010_sdsd[200];

    auto g_y_0_z_0_0_yy_0_yy = buffer_1010_sdsd[201];

    auto g_y_0_z_0_0_yy_0_yz = buffer_1010_sdsd[202];

    auto g_y_0_z_0_0_yy_0_zz = buffer_1010_sdsd[203];

    auto g_y_0_z_0_0_yz_0_xx = buffer_1010_sdsd[204];

    auto g_y_0_z_0_0_yz_0_xy = buffer_1010_sdsd[205];

    auto g_y_0_z_0_0_yz_0_xz = buffer_1010_sdsd[206];

    auto g_y_0_z_0_0_yz_0_yy = buffer_1010_sdsd[207];

    auto g_y_0_z_0_0_yz_0_yz = buffer_1010_sdsd[208];

    auto g_y_0_z_0_0_yz_0_zz = buffer_1010_sdsd[209];

    auto g_y_0_z_0_0_zz_0_xx = buffer_1010_sdsd[210];

    auto g_y_0_z_0_0_zz_0_xy = buffer_1010_sdsd[211];

    auto g_y_0_z_0_0_zz_0_xz = buffer_1010_sdsd[212];

    auto g_y_0_z_0_0_zz_0_yy = buffer_1010_sdsd[213];

    auto g_y_0_z_0_0_zz_0_yz = buffer_1010_sdsd[214];

    auto g_y_0_z_0_0_zz_0_zz = buffer_1010_sdsd[215];

    auto g_z_0_x_0_0_xx_0_xx = buffer_1010_sdsd[216];

    auto g_z_0_x_0_0_xx_0_xy = buffer_1010_sdsd[217];

    auto g_z_0_x_0_0_xx_0_xz = buffer_1010_sdsd[218];

    auto g_z_0_x_0_0_xx_0_yy = buffer_1010_sdsd[219];

    auto g_z_0_x_0_0_xx_0_yz = buffer_1010_sdsd[220];

    auto g_z_0_x_0_0_xx_0_zz = buffer_1010_sdsd[221];

    auto g_z_0_x_0_0_xy_0_xx = buffer_1010_sdsd[222];

    auto g_z_0_x_0_0_xy_0_xy = buffer_1010_sdsd[223];

    auto g_z_0_x_0_0_xy_0_xz = buffer_1010_sdsd[224];

    auto g_z_0_x_0_0_xy_0_yy = buffer_1010_sdsd[225];

    auto g_z_0_x_0_0_xy_0_yz = buffer_1010_sdsd[226];

    auto g_z_0_x_0_0_xy_0_zz = buffer_1010_sdsd[227];

    auto g_z_0_x_0_0_xz_0_xx = buffer_1010_sdsd[228];

    auto g_z_0_x_0_0_xz_0_xy = buffer_1010_sdsd[229];

    auto g_z_0_x_0_0_xz_0_xz = buffer_1010_sdsd[230];

    auto g_z_0_x_0_0_xz_0_yy = buffer_1010_sdsd[231];

    auto g_z_0_x_0_0_xz_0_yz = buffer_1010_sdsd[232];

    auto g_z_0_x_0_0_xz_0_zz = buffer_1010_sdsd[233];

    auto g_z_0_x_0_0_yy_0_xx = buffer_1010_sdsd[234];

    auto g_z_0_x_0_0_yy_0_xy = buffer_1010_sdsd[235];

    auto g_z_0_x_0_0_yy_0_xz = buffer_1010_sdsd[236];

    auto g_z_0_x_0_0_yy_0_yy = buffer_1010_sdsd[237];

    auto g_z_0_x_0_0_yy_0_yz = buffer_1010_sdsd[238];

    auto g_z_0_x_0_0_yy_0_zz = buffer_1010_sdsd[239];

    auto g_z_0_x_0_0_yz_0_xx = buffer_1010_sdsd[240];

    auto g_z_0_x_0_0_yz_0_xy = buffer_1010_sdsd[241];

    auto g_z_0_x_0_0_yz_0_xz = buffer_1010_sdsd[242];

    auto g_z_0_x_0_0_yz_0_yy = buffer_1010_sdsd[243];

    auto g_z_0_x_0_0_yz_0_yz = buffer_1010_sdsd[244];

    auto g_z_0_x_0_0_yz_0_zz = buffer_1010_sdsd[245];

    auto g_z_0_x_0_0_zz_0_xx = buffer_1010_sdsd[246];

    auto g_z_0_x_0_0_zz_0_xy = buffer_1010_sdsd[247];

    auto g_z_0_x_0_0_zz_0_xz = buffer_1010_sdsd[248];

    auto g_z_0_x_0_0_zz_0_yy = buffer_1010_sdsd[249];

    auto g_z_0_x_0_0_zz_0_yz = buffer_1010_sdsd[250];

    auto g_z_0_x_0_0_zz_0_zz = buffer_1010_sdsd[251];

    auto g_z_0_y_0_0_xx_0_xx = buffer_1010_sdsd[252];

    auto g_z_0_y_0_0_xx_0_xy = buffer_1010_sdsd[253];

    auto g_z_0_y_0_0_xx_0_xz = buffer_1010_sdsd[254];

    auto g_z_0_y_0_0_xx_0_yy = buffer_1010_sdsd[255];

    auto g_z_0_y_0_0_xx_0_yz = buffer_1010_sdsd[256];

    auto g_z_0_y_0_0_xx_0_zz = buffer_1010_sdsd[257];

    auto g_z_0_y_0_0_xy_0_xx = buffer_1010_sdsd[258];

    auto g_z_0_y_0_0_xy_0_xy = buffer_1010_sdsd[259];

    auto g_z_0_y_0_0_xy_0_xz = buffer_1010_sdsd[260];

    auto g_z_0_y_0_0_xy_0_yy = buffer_1010_sdsd[261];

    auto g_z_0_y_0_0_xy_0_yz = buffer_1010_sdsd[262];

    auto g_z_0_y_0_0_xy_0_zz = buffer_1010_sdsd[263];

    auto g_z_0_y_0_0_xz_0_xx = buffer_1010_sdsd[264];

    auto g_z_0_y_0_0_xz_0_xy = buffer_1010_sdsd[265];

    auto g_z_0_y_0_0_xz_0_xz = buffer_1010_sdsd[266];

    auto g_z_0_y_0_0_xz_0_yy = buffer_1010_sdsd[267];

    auto g_z_0_y_0_0_xz_0_yz = buffer_1010_sdsd[268];

    auto g_z_0_y_0_0_xz_0_zz = buffer_1010_sdsd[269];

    auto g_z_0_y_0_0_yy_0_xx = buffer_1010_sdsd[270];

    auto g_z_0_y_0_0_yy_0_xy = buffer_1010_sdsd[271];

    auto g_z_0_y_0_0_yy_0_xz = buffer_1010_sdsd[272];

    auto g_z_0_y_0_0_yy_0_yy = buffer_1010_sdsd[273];

    auto g_z_0_y_0_0_yy_0_yz = buffer_1010_sdsd[274];

    auto g_z_0_y_0_0_yy_0_zz = buffer_1010_sdsd[275];

    auto g_z_0_y_0_0_yz_0_xx = buffer_1010_sdsd[276];

    auto g_z_0_y_0_0_yz_0_xy = buffer_1010_sdsd[277];

    auto g_z_0_y_0_0_yz_0_xz = buffer_1010_sdsd[278];

    auto g_z_0_y_0_0_yz_0_yy = buffer_1010_sdsd[279];

    auto g_z_0_y_0_0_yz_0_yz = buffer_1010_sdsd[280];

    auto g_z_0_y_0_0_yz_0_zz = buffer_1010_sdsd[281];

    auto g_z_0_y_0_0_zz_0_xx = buffer_1010_sdsd[282];

    auto g_z_0_y_0_0_zz_0_xy = buffer_1010_sdsd[283];

    auto g_z_0_y_0_0_zz_0_xz = buffer_1010_sdsd[284];

    auto g_z_0_y_0_0_zz_0_yy = buffer_1010_sdsd[285];

    auto g_z_0_y_0_0_zz_0_yz = buffer_1010_sdsd[286];

    auto g_z_0_y_0_0_zz_0_zz = buffer_1010_sdsd[287];

    auto g_z_0_z_0_0_xx_0_xx = buffer_1010_sdsd[288];

    auto g_z_0_z_0_0_xx_0_xy = buffer_1010_sdsd[289];

    auto g_z_0_z_0_0_xx_0_xz = buffer_1010_sdsd[290];

    auto g_z_0_z_0_0_xx_0_yy = buffer_1010_sdsd[291];

    auto g_z_0_z_0_0_xx_0_yz = buffer_1010_sdsd[292];

    auto g_z_0_z_0_0_xx_0_zz = buffer_1010_sdsd[293];

    auto g_z_0_z_0_0_xy_0_xx = buffer_1010_sdsd[294];

    auto g_z_0_z_0_0_xy_0_xy = buffer_1010_sdsd[295];

    auto g_z_0_z_0_0_xy_0_xz = buffer_1010_sdsd[296];

    auto g_z_0_z_0_0_xy_0_yy = buffer_1010_sdsd[297];

    auto g_z_0_z_0_0_xy_0_yz = buffer_1010_sdsd[298];

    auto g_z_0_z_0_0_xy_0_zz = buffer_1010_sdsd[299];

    auto g_z_0_z_0_0_xz_0_xx = buffer_1010_sdsd[300];

    auto g_z_0_z_0_0_xz_0_xy = buffer_1010_sdsd[301];

    auto g_z_0_z_0_0_xz_0_xz = buffer_1010_sdsd[302];

    auto g_z_0_z_0_0_xz_0_yy = buffer_1010_sdsd[303];

    auto g_z_0_z_0_0_xz_0_yz = buffer_1010_sdsd[304];

    auto g_z_0_z_0_0_xz_0_zz = buffer_1010_sdsd[305];

    auto g_z_0_z_0_0_yy_0_xx = buffer_1010_sdsd[306];

    auto g_z_0_z_0_0_yy_0_xy = buffer_1010_sdsd[307];

    auto g_z_0_z_0_0_yy_0_xz = buffer_1010_sdsd[308];

    auto g_z_0_z_0_0_yy_0_yy = buffer_1010_sdsd[309];

    auto g_z_0_z_0_0_yy_0_yz = buffer_1010_sdsd[310];

    auto g_z_0_z_0_0_yy_0_zz = buffer_1010_sdsd[311];

    auto g_z_0_z_0_0_yz_0_xx = buffer_1010_sdsd[312];

    auto g_z_0_z_0_0_yz_0_xy = buffer_1010_sdsd[313];

    auto g_z_0_z_0_0_yz_0_xz = buffer_1010_sdsd[314];

    auto g_z_0_z_0_0_yz_0_yy = buffer_1010_sdsd[315];

    auto g_z_0_z_0_0_yz_0_yz = buffer_1010_sdsd[316];

    auto g_z_0_z_0_0_yz_0_zz = buffer_1010_sdsd[317];

    auto g_z_0_z_0_0_zz_0_xx = buffer_1010_sdsd[318];

    auto g_z_0_z_0_0_zz_0_xy = buffer_1010_sdsd[319];

    auto g_z_0_z_0_0_zz_0_xz = buffer_1010_sdsd[320];

    auto g_z_0_z_0_0_zz_0_yy = buffer_1010_sdsd[321];

    auto g_z_0_z_0_0_zz_0_yz = buffer_1010_sdsd[322];

    auto g_z_0_z_0_0_zz_0_zz = buffer_1010_sdsd[323];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_x_0_0_xx_0_xx, g_x_0_x_0_0_xx_0_xy, g_x_0_x_0_0_xx_0_xz, g_x_0_x_0_0_xx_0_yy, g_x_0_x_0_0_xx_0_yz, g_x_0_x_0_0_xx_0_zz, g_x_xx_x_xx, g_x_xx_x_xy, g_x_xx_x_xz, g_x_xx_x_yy, g_x_xx_x_yz, g_x_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xx_0_xx[i] = 4.0 * g_x_xx_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_0_xy[i] = 4.0 * g_x_xx_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_0_xz[i] = 4.0 * g_x_xx_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_0_yy[i] = 4.0 * g_x_xx_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_0_yz[i] = 4.0 * g_x_xx_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_0_zz[i] = 4.0 * g_x_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_x_0_0_xy_0_xx, g_x_0_x_0_0_xy_0_xy, g_x_0_x_0_0_xy_0_xz, g_x_0_x_0_0_xy_0_yy, g_x_0_x_0_0_xy_0_yz, g_x_0_x_0_0_xy_0_zz, g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xy_0_xx[i] = 4.0 * g_x_xy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_0_xy[i] = 4.0 * g_x_xy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_0_xz[i] = 4.0 * g_x_xy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_0_yy[i] = 4.0 * g_x_xy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_0_yz[i] = 4.0 * g_x_xy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_0_zz[i] = 4.0 * g_x_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_x_0_0_xz_0_xx, g_x_0_x_0_0_xz_0_xy, g_x_0_x_0_0_xz_0_xz, g_x_0_x_0_0_xz_0_yy, g_x_0_x_0_0_xz_0_yz, g_x_0_x_0_0_xz_0_zz, g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xz_0_xx[i] = 4.0 * g_x_xz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_0_xy[i] = 4.0 * g_x_xz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_0_xz[i] = 4.0 * g_x_xz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_0_yy[i] = 4.0 * g_x_xz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_0_yz[i] = 4.0 * g_x_xz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_0_zz[i] = 4.0 * g_x_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_x_0_0_yy_0_xx, g_x_0_x_0_0_yy_0_xy, g_x_0_x_0_0_yy_0_xz, g_x_0_x_0_0_yy_0_yy, g_x_0_x_0_0_yy_0_yz, g_x_0_x_0_0_yy_0_zz, g_x_yy_x_xx, g_x_yy_x_xy, g_x_yy_x_xz, g_x_yy_x_yy, g_x_yy_x_yz, g_x_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yy_0_xx[i] = 4.0 * g_x_yy_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_0_xy[i] = 4.0 * g_x_yy_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_0_xz[i] = 4.0 * g_x_yy_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_0_yy[i] = 4.0 * g_x_yy_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_0_yz[i] = 4.0 * g_x_yy_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_0_zz[i] = 4.0 * g_x_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_x_0_0_yz_0_xx, g_x_0_x_0_0_yz_0_xy, g_x_0_x_0_0_yz_0_xz, g_x_0_x_0_0_yz_0_yy, g_x_0_x_0_0_yz_0_yz, g_x_0_x_0_0_yz_0_zz, g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yz_0_xx[i] = 4.0 * g_x_yz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_0_xy[i] = 4.0 * g_x_yz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_0_xz[i] = 4.0 * g_x_yz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_0_yy[i] = 4.0 * g_x_yz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_0_yz[i] = 4.0 * g_x_yz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_0_zz[i] = 4.0 * g_x_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_x_0_0_zz_0_xx, g_x_0_x_0_0_zz_0_xy, g_x_0_x_0_0_zz_0_xz, g_x_0_x_0_0_zz_0_yy, g_x_0_x_0_0_zz_0_yz, g_x_0_x_0_0_zz_0_zz, g_x_zz_x_xx, g_x_zz_x_xy, g_x_zz_x_xz, g_x_zz_x_yy, g_x_zz_x_yz, g_x_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_zz_0_xx[i] = 4.0 * g_x_zz_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_0_xy[i] = 4.0 * g_x_zz_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_0_xz[i] = 4.0 * g_x_zz_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_0_yy[i] = 4.0 * g_x_zz_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_0_yz[i] = 4.0 * g_x_zz_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_0_zz[i] = 4.0 * g_x_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_y_0_0_xx_0_xx, g_x_0_y_0_0_xx_0_xy, g_x_0_y_0_0_xx_0_xz, g_x_0_y_0_0_xx_0_yy, g_x_0_y_0_0_xx_0_yz, g_x_0_y_0_0_xx_0_zz, g_x_xx_y_xx, g_x_xx_y_xy, g_x_xx_y_xz, g_x_xx_y_yy, g_x_xx_y_yz, g_x_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xx_0_xx[i] = 4.0 * g_x_xx_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_0_xy[i] = 4.0 * g_x_xx_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_0_xz[i] = 4.0 * g_x_xx_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_0_yy[i] = 4.0 * g_x_xx_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_0_yz[i] = 4.0 * g_x_xx_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_0_zz[i] = 4.0 * g_x_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_y_0_0_xy_0_xx, g_x_0_y_0_0_xy_0_xy, g_x_0_y_0_0_xy_0_xz, g_x_0_y_0_0_xy_0_yy, g_x_0_y_0_0_xy_0_yz, g_x_0_y_0_0_xy_0_zz, g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xy_0_xx[i] = 4.0 * g_x_xy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_0_xy[i] = 4.0 * g_x_xy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_0_xz[i] = 4.0 * g_x_xy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_0_yy[i] = 4.0 * g_x_xy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_0_yz[i] = 4.0 * g_x_xy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_0_zz[i] = 4.0 * g_x_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_y_0_0_xz_0_xx, g_x_0_y_0_0_xz_0_xy, g_x_0_y_0_0_xz_0_xz, g_x_0_y_0_0_xz_0_yy, g_x_0_y_0_0_xz_0_yz, g_x_0_y_0_0_xz_0_zz, g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xz_0_xx[i] = 4.0 * g_x_xz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_0_xy[i] = 4.0 * g_x_xz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_0_xz[i] = 4.0 * g_x_xz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_0_yy[i] = 4.0 * g_x_xz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_0_yz[i] = 4.0 * g_x_xz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_0_zz[i] = 4.0 * g_x_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_0_y_0_0_yy_0_xx, g_x_0_y_0_0_yy_0_xy, g_x_0_y_0_0_yy_0_xz, g_x_0_y_0_0_yy_0_yy, g_x_0_y_0_0_yy_0_yz, g_x_0_y_0_0_yy_0_zz, g_x_yy_y_xx, g_x_yy_y_xy, g_x_yy_y_xz, g_x_yy_y_yy, g_x_yy_y_yz, g_x_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yy_0_xx[i] = 4.0 * g_x_yy_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_0_xy[i] = 4.0 * g_x_yy_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_0_xz[i] = 4.0 * g_x_yy_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_0_yy[i] = 4.0 * g_x_yy_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_0_yz[i] = 4.0 * g_x_yy_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_0_zz[i] = 4.0 * g_x_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_0_y_0_0_yz_0_xx, g_x_0_y_0_0_yz_0_xy, g_x_0_y_0_0_yz_0_xz, g_x_0_y_0_0_yz_0_yy, g_x_0_y_0_0_yz_0_yz, g_x_0_y_0_0_yz_0_zz, g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yz_0_xx[i] = 4.0 * g_x_yz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_0_xy[i] = 4.0 * g_x_yz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_0_xz[i] = 4.0 * g_x_yz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_0_yy[i] = 4.0 * g_x_yz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_0_yz[i] = 4.0 * g_x_yz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_0_zz[i] = 4.0 * g_x_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_0_y_0_0_zz_0_xx, g_x_0_y_0_0_zz_0_xy, g_x_0_y_0_0_zz_0_xz, g_x_0_y_0_0_zz_0_yy, g_x_0_y_0_0_zz_0_yz, g_x_0_y_0_0_zz_0_zz, g_x_zz_y_xx, g_x_zz_y_xy, g_x_zz_y_xz, g_x_zz_y_yy, g_x_zz_y_yz, g_x_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_zz_0_xx[i] = 4.0 * g_x_zz_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_0_xy[i] = 4.0 * g_x_zz_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_0_xz[i] = 4.0 * g_x_zz_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_0_yy[i] = 4.0 * g_x_zz_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_0_yz[i] = 4.0 * g_x_zz_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_0_zz[i] = 4.0 * g_x_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_z_0_0_xx_0_xx, g_x_0_z_0_0_xx_0_xy, g_x_0_z_0_0_xx_0_xz, g_x_0_z_0_0_xx_0_yy, g_x_0_z_0_0_xx_0_yz, g_x_0_z_0_0_xx_0_zz, g_x_xx_z_xx, g_x_xx_z_xy, g_x_xx_z_xz, g_x_xx_z_yy, g_x_xx_z_yz, g_x_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xx_0_xx[i] = 4.0 * g_x_xx_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_0_xy[i] = 4.0 * g_x_xx_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_0_xz[i] = 4.0 * g_x_xx_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_0_yy[i] = 4.0 * g_x_xx_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_0_yz[i] = 4.0 * g_x_xx_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_0_zz[i] = 4.0 * g_x_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_z_0_0_xy_0_xx, g_x_0_z_0_0_xy_0_xy, g_x_0_z_0_0_xy_0_xz, g_x_0_z_0_0_xy_0_yy, g_x_0_z_0_0_xy_0_yz, g_x_0_z_0_0_xy_0_zz, g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xy_0_xx[i] = 4.0 * g_x_xy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_0_xy[i] = 4.0 * g_x_xy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_0_xz[i] = 4.0 * g_x_xy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_0_yy[i] = 4.0 * g_x_xy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_0_yz[i] = 4.0 * g_x_xy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_0_zz[i] = 4.0 * g_x_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_z_0_0_xz_0_xx, g_x_0_z_0_0_xz_0_xy, g_x_0_z_0_0_xz_0_xz, g_x_0_z_0_0_xz_0_yy, g_x_0_z_0_0_xz_0_yz, g_x_0_z_0_0_xz_0_zz, g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xz_0_xx[i] = 4.0 * g_x_xz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_0_xy[i] = 4.0 * g_x_xz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_0_xz[i] = 4.0 * g_x_xz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_0_yy[i] = 4.0 * g_x_xz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_0_yz[i] = 4.0 * g_x_xz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_0_zz[i] = 4.0 * g_x_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_z_0_0_yy_0_xx, g_x_0_z_0_0_yy_0_xy, g_x_0_z_0_0_yy_0_xz, g_x_0_z_0_0_yy_0_yy, g_x_0_z_0_0_yy_0_yz, g_x_0_z_0_0_yy_0_zz, g_x_yy_z_xx, g_x_yy_z_xy, g_x_yy_z_xz, g_x_yy_z_yy, g_x_yy_z_yz, g_x_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yy_0_xx[i] = 4.0 * g_x_yy_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_0_xy[i] = 4.0 * g_x_yy_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_0_xz[i] = 4.0 * g_x_yy_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_0_yy[i] = 4.0 * g_x_yy_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_0_yz[i] = 4.0 * g_x_yy_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_0_zz[i] = 4.0 * g_x_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_z_0_0_yz_0_xx, g_x_0_z_0_0_yz_0_xy, g_x_0_z_0_0_yz_0_xz, g_x_0_z_0_0_yz_0_yy, g_x_0_z_0_0_yz_0_yz, g_x_0_z_0_0_yz_0_zz, g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yz_0_xx[i] = 4.0 * g_x_yz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_0_xy[i] = 4.0 * g_x_yz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_0_xz[i] = 4.0 * g_x_yz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_0_yy[i] = 4.0 * g_x_yz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_0_yz[i] = 4.0 * g_x_yz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_0_zz[i] = 4.0 * g_x_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_z_0_0_zz_0_xx, g_x_0_z_0_0_zz_0_xy, g_x_0_z_0_0_zz_0_xz, g_x_0_z_0_0_zz_0_yy, g_x_0_z_0_0_zz_0_yz, g_x_0_z_0_0_zz_0_zz, g_x_zz_z_xx, g_x_zz_z_xy, g_x_zz_z_xz, g_x_zz_z_yy, g_x_zz_z_yz, g_x_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_zz_0_xx[i] = 4.0 * g_x_zz_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_0_xy[i] = 4.0 * g_x_zz_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_0_xz[i] = 4.0 * g_x_zz_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_0_yy[i] = 4.0 * g_x_zz_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_0_yz[i] = 4.0 * g_x_zz_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_0_zz[i] = 4.0 * g_x_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_y_0_x_0_0_xx_0_xx, g_y_0_x_0_0_xx_0_xy, g_y_0_x_0_0_xx_0_xz, g_y_0_x_0_0_xx_0_yy, g_y_0_x_0_0_xx_0_yz, g_y_0_x_0_0_xx_0_zz, g_y_xx_x_xx, g_y_xx_x_xy, g_y_xx_x_xz, g_y_xx_x_yy, g_y_xx_x_yz, g_y_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xx_0_xx[i] = 4.0 * g_y_xx_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_0_xy[i] = 4.0 * g_y_xx_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_0_xz[i] = 4.0 * g_y_xx_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_0_yy[i] = 4.0 * g_y_xx_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_0_yz[i] = 4.0 * g_y_xx_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_0_zz[i] = 4.0 * g_y_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_y_0_x_0_0_xy_0_xx, g_y_0_x_0_0_xy_0_xy, g_y_0_x_0_0_xy_0_xz, g_y_0_x_0_0_xy_0_yy, g_y_0_x_0_0_xy_0_yz, g_y_0_x_0_0_xy_0_zz, g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xy_0_xx[i] = 4.0 * g_y_xy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_0_xy[i] = 4.0 * g_y_xy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_0_xz[i] = 4.0 * g_y_xy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_0_yy[i] = 4.0 * g_y_xy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_0_yz[i] = 4.0 * g_y_xy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_0_zz[i] = 4.0 * g_y_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_y_0_x_0_0_xz_0_xx, g_y_0_x_0_0_xz_0_xy, g_y_0_x_0_0_xz_0_xz, g_y_0_x_0_0_xz_0_yy, g_y_0_x_0_0_xz_0_yz, g_y_0_x_0_0_xz_0_zz, g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xz_0_xx[i] = 4.0 * g_y_xz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_0_xy[i] = 4.0 * g_y_xz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_0_xz[i] = 4.0 * g_y_xz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_0_yy[i] = 4.0 * g_y_xz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_0_yz[i] = 4.0 * g_y_xz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_0_zz[i] = 4.0 * g_y_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_y_0_x_0_0_yy_0_xx, g_y_0_x_0_0_yy_0_xy, g_y_0_x_0_0_yy_0_xz, g_y_0_x_0_0_yy_0_yy, g_y_0_x_0_0_yy_0_yz, g_y_0_x_0_0_yy_0_zz, g_y_yy_x_xx, g_y_yy_x_xy, g_y_yy_x_xz, g_y_yy_x_yy, g_y_yy_x_yz, g_y_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yy_0_xx[i] = 4.0 * g_y_yy_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_0_xy[i] = 4.0 * g_y_yy_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_0_xz[i] = 4.0 * g_y_yy_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_0_yy[i] = 4.0 * g_y_yy_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_0_yz[i] = 4.0 * g_y_yy_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_0_zz[i] = 4.0 * g_y_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_y_0_x_0_0_yz_0_xx, g_y_0_x_0_0_yz_0_xy, g_y_0_x_0_0_yz_0_xz, g_y_0_x_0_0_yz_0_yy, g_y_0_x_0_0_yz_0_yz, g_y_0_x_0_0_yz_0_zz, g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yz_0_xx[i] = 4.0 * g_y_yz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_0_xy[i] = 4.0 * g_y_yz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_0_xz[i] = 4.0 * g_y_yz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_0_yy[i] = 4.0 * g_y_yz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_0_yz[i] = 4.0 * g_y_yz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_0_zz[i] = 4.0 * g_y_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_y_0_x_0_0_zz_0_xx, g_y_0_x_0_0_zz_0_xy, g_y_0_x_0_0_zz_0_xz, g_y_0_x_0_0_zz_0_yy, g_y_0_x_0_0_zz_0_yz, g_y_0_x_0_0_zz_0_zz, g_y_zz_x_xx, g_y_zz_x_xy, g_y_zz_x_xz, g_y_zz_x_yy, g_y_zz_x_yz, g_y_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_zz_0_xx[i] = 4.0 * g_y_zz_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_0_xy[i] = 4.0 * g_y_zz_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_0_xz[i] = 4.0 * g_y_zz_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_0_yy[i] = 4.0 * g_y_zz_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_0_yz[i] = 4.0 * g_y_zz_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_0_zz[i] = 4.0 * g_y_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_y_0_y_0_0_xx_0_xx, g_y_0_y_0_0_xx_0_xy, g_y_0_y_0_0_xx_0_xz, g_y_0_y_0_0_xx_0_yy, g_y_0_y_0_0_xx_0_yz, g_y_0_y_0_0_xx_0_zz, g_y_xx_y_xx, g_y_xx_y_xy, g_y_xx_y_xz, g_y_xx_y_yy, g_y_xx_y_yz, g_y_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xx_0_xx[i] = 4.0 * g_y_xx_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_0_xy[i] = 4.0 * g_y_xx_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_0_xz[i] = 4.0 * g_y_xx_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_0_yy[i] = 4.0 * g_y_xx_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_0_yz[i] = 4.0 * g_y_xx_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_0_zz[i] = 4.0 * g_y_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_y_0_y_0_0_xy_0_xx, g_y_0_y_0_0_xy_0_xy, g_y_0_y_0_0_xy_0_xz, g_y_0_y_0_0_xy_0_yy, g_y_0_y_0_0_xy_0_yz, g_y_0_y_0_0_xy_0_zz, g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xy_0_xx[i] = 4.0 * g_y_xy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_0_xy[i] = 4.0 * g_y_xy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_0_xz[i] = 4.0 * g_y_xy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_0_yy[i] = 4.0 * g_y_xy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_0_yz[i] = 4.0 * g_y_xy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_0_zz[i] = 4.0 * g_y_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_y_0_y_0_0_xz_0_xx, g_y_0_y_0_0_xz_0_xy, g_y_0_y_0_0_xz_0_xz, g_y_0_y_0_0_xz_0_yy, g_y_0_y_0_0_xz_0_yz, g_y_0_y_0_0_xz_0_zz, g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xz_0_xx[i] = 4.0 * g_y_xz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_0_xy[i] = 4.0 * g_y_xz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_0_xz[i] = 4.0 * g_y_xz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_0_yy[i] = 4.0 * g_y_xz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_0_yz[i] = 4.0 * g_y_xz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_0_zz[i] = 4.0 * g_y_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_y_0_y_0_0_yy_0_xx, g_y_0_y_0_0_yy_0_xy, g_y_0_y_0_0_yy_0_xz, g_y_0_y_0_0_yy_0_yy, g_y_0_y_0_0_yy_0_yz, g_y_0_y_0_0_yy_0_zz, g_y_yy_y_xx, g_y_yy_y_xy, g_y_yy_y_xz, g_y_yy_y_yy, g_y_yy_y_yz, g_y_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yy_0_xx[i] = 4.0 * g_y_yy_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_0_xy[i] = 4.0 * g_y_yy_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_0_xz[i] = 4.0 * g_y_yy_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_0_yy[i] = 4.0 * g_y_yy_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_0_yz[i] = 4.0 * g_y_yy_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_0_zz[i] = 4.0 * g_y_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_y_0_y_0_0_yz_0_xx, g_y_0_y_0_0_yz_0_xy, g_y_0_y_0_0_yz_0_xz, g_y_0_y_0_0_yz_0_yy, g_y_0_y_0_0_yz_0_yz, g_y_0_y_0_0_yz_0_zz, g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yz_0_xx[i] = 4.0 * g_y_yz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_0_xy[i] = 4.0 * g_y_yz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_0_xz[i] = 4.0 * g_y_yz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_0_yy[i] = 4.0 * g_y_yz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_0_yz[i] = 4.0 * g_y_yz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_0_zz[i] = 4.0 * g_y_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_y_0_y_0_0_zz_0_xx, g_y_0_y_0_0_zz_0_xy, g_y_0_y_0_0_zz_0_xz, g_y_0_y_0_0_zz_0_yy, g_y_0_y_0_0_zz_0_yz, g_y_0_y_0_0_zz_0_zz, g_y_zz_y_xx, g_y_zz_y_xy, g_y_zz_y_xz, g_y_zz_y_yy, g_y_zz_y_yz, g_y_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_zz_0_xx[i] = 4.0 * g_y_zz_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_0_xy[i] = 4.0 * g_y_zz_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_0_xz[i] = 4.0 * g_y_zz_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_0_yy[i] = 4.0 * g_y_zz_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_0_yz[i] = 4.0 * g_y_zz_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_0_zz[i] = 4.0 * g_y_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_y_0_z_0_0_xx_0_xx, g_y_0_z_0_0_xx_0_xy, g_y_0_z_0_0_xx_0_xz, g_y_0_z_0_0_xx_0_yy, g_y_0_z_0_0_xx_0_yz, g_y_0_z_0_0_xx_0_zz, g_y_xx_z_xx, g_y_xx_z_xy, g_y_xx_z_xz, g_y_xx_z_yy, g_y_xx_z_yz, g_y_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xx_0_xx[i] = 4.0 * g_y_xx_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_0_xy[i] = 4.0 * g_y_xx_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_0_xz[i] = 4.0 * g_y_xx_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_0_yy[i] = 4.0 * g_y_xx_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_0_yz[i] = 4.0 * g_y_xx_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_0_zz[i] = 4.0 * g_y_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_y_0_z_0_0_xy_0_xx, g_y_0_z_0_0_xy_0_xy, g_y_0_z_0_0_xy_0_xz, g_y_0_z_0_0_xy_0_yy, g_y_0_z_0_0_xy_0_yz, g_y_0_z_0_0_xy_0_zz, g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xy_0_xx[i] = 4.0 * g_y_xy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_0_xy[i] = 4.0 * g_y_xy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_0_xz[i] = 4.0 * g_y_xy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_0_yy[i] = 4.0 * g_y_xy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_0_yz[i] = 4.0 * g_y_xy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_0_zz[i] = 4.0 * g_y_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_y_0_z_0_0_xz_0_xx, g_y_0_z_0_0_xz_0_xy, g_y_0_z_0_0_xz_0_xz, g_y_0_z_0_0_xz_0_yy, g_y_0_z_0_0_xz_0_yz, g_y_0_z_0_0_xz_0_zz, g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xz_0_xx[i] = 4.0 * g_y_xz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_0_xy[i] = 4.0 * g_y_xz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_0_xz[i] = 4.0 * g_y_xz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_0_yy[i] = 4.0 * g_y_xz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_0_yz[i] = 4.0 * g_y_xz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_0_zz[i] = 4.0 * g_y_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_y_0_z_0_0_yy_0_xx, g_y_0_z_0_0_yy_0_xy, g_y_0_z_0_0_yy_0_xz, g_y_0_z_0_0_yy_0_yy, g_y_0_z_0_0_yy_0_yz, g_y_0_z_0_0_yy_0_zz, g_y_yy_z_xx, g_y_yy_z_xy, g_y_yy_z_xz, g_y_yy_z_yy, g_y_yy_z_yz, g_y_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yy_0_xx[i] = 4.0 * g_y_yy_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_0_xy[i] = 4.0 * g_y_yy_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_0_xz[i] = 4.0 * g_y_yy_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_0_yy[i] = 4.0 * g_y_yy_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_0_yz[i] = 4.0 * g_y_yy_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_0_zz[i] = 4.0 * g_y_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_y_0_z_0_0_yz_0_xx, g_y_0_z_0_0_yz_0_xy, g_y_0_z_0_0_yz_0_xz, g_y_0_z_0_0_yz_0_yy, g_y_0_z_0_0_yz_0_yz, g_y_0_z_0_0_yz_0_zz, g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yz_0_xx[i] = 4.0 * g_y_yz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_0_xy[i] = 4.0 * g_y_yz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_0_xz[i] = 4.0 * g_y_yz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_0_yy[i] = 4.0 * g_y_yz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_0_yz[i] = 4.0 * g_y_yz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_0_zz[i] = 4.0 * g_y_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_y_0_z_0_0_zz_0_xx, g_y_0_z_0_0_zz_0_xy, g_y_0_z_0_0_zz_0_xz, g_y_0_z_0_0_zz_0_yy, g_y_0_z_0_0_zz_0_yz, g_y_0_z_0_0_zz_0_zz, g_y_zz_z_xx, g_y_zz_z_xy, g_y_zz_z_xz, g_y_zz_z_yy, g_y_zz_z_yz, g_y_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_zz_0_xx[i] = 4.0 * g_y_zz_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_0_xy[i] = 4.0 * g_y_zz_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_0_xz[i] = 4.0 * g_y_zz_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_0_yy[i] = 4.0 * g_y_zz_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_0_yz[i] = 4.0 * g_y_zz_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_0_zz[i] = 4.0 * g_y_zz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_z_0_x_0_0_xx_0_xx, g_z_0_x_0_0_xx_0_xy, g_z_0_x_0_0_xx_0_xz, g_z_0_x_0_0_xx_0_yy, g_z_0_x_0_0_xx_0_yz, g_z_0_x_0_0_xx_0_zz, g_z_xx_x_xx, g_z_xx_x_xy, g_z_xx_x_xz, g_z_xx_x_yy, g_z_xx_x_yz, g_z_xx_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xx_0_xx[i] = 4.0 * g_z_xx_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_0_xy[i] = 4.0 * g_z_xx_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_0_xz[i] = 4.0 * g_z_xx_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_0_yy[i] = 4.0 * g_z_xx_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_0_yz[i] = 4.0 * g_z_xx_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_0_zz[i] = 4.0 * g_z_xx_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_z_0_x_0_0_xy_0_xx, g_z_0_x_0_0_xy_0_xy, g_z_0_x_0_0_xy_0_xz, g_z_0_x_0_0_xy_0_yy, g_z_0_x_0_0_xy_0_yz, g_z_0_x_0_0_xy_0_zz, g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xy_0_xx[i] = 4.0 * g_z_xy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_0_xy[i] = 4.0 * g_z_xy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_0_xz[i] = 4.0 * g_z_xy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_0_yy[i] = 4.0 * g_z_xy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_0_yz[i] = 4.0 * g_z_xy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_0_zz[i] = 4.0 * g_z_xy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_z_0_x_0_0_xz_0_xx, g_z_0_x_0_0_xz_0_xy, g_z_0_x_0_0_xz_0_xz, g_z_0_x_0_0_xz_0_yy, g_z_0_x_0_0_xz_0_yz, g_z_0_x_0_0_xz_0_zz, g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xz_0_xx[i] = 4.0 * g_z_xz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_0_xy[i] = 4.0 * g_z_xz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_0_xz[i] = 4.0 * g_z_xz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_0_yy[i] = 4.0 * g_z_xz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_0_yz[i] = 4.0 * g_z_xz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_0_zz[i] = 4.0 * g_z_xz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_z_0_x_0_0_yy_0_xx, g_z_0_x_0_0_yy_0_xy, g_z_0_x_0_0_yy_0_xz, g_z_0_x_0_0_yy_0_yy, g_z_0_x_0_0_yy_0_yz, g_z_0_x_0_0_yy_0_zz, g_z_yy_x_xx, g_z_yy_x_xy, g_z_yy_x_xz, g_z_yy_x_yy, g_z_yy_x_yz, g_z_yy_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yy_0_xx[i] = 4.0 * g_z_yy_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_0_xy[i] = 4.0 * g_z_yy_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_0_xz[i] = 4.0 * g_z_yy_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_0_yy[i] = 4.0 * g_z_yy_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_0_yz[i] = 4.0 * g_z_yy_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_0_zz[i] = 4.0 * g_z_yy_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_z_0_x_0_0_yz_0_xx, g_z_0_x_0_0_yz_0_xy, g_z_0_x_0_0_yz_0_xz, g_z_0_x_0_0_yz_0_yy, g_z_0_x_0_0_yz_0_yz, g_z_0_x_0_0_yz_0_zz, g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yz_0_xx[i] = 4.0 * g_z_yz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_0_xy[i] = 4.0 * g_z_yz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_0_xz[i] = 4.0 * g_z_yz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_0_yy[i] = 4.0 * g_z_yz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_0_yz[i] = 4.0 * g_z_yz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_0_zz[i] = 4.0 * g_z_yz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_z_0_x_0_0_zz_0_xx, g_z_0_x_0_0_zz_0_xy, g_z_0_x_0_0_zz_0_xz, g_z_0_x_0_0_zz_0_yy, g_z_0_x_0_0_zz_0_yz, g_z_0_x_0_0_zz_0_zz, g_z_zz_x_xx, g_z_zz_x_xy, g_z_zz_x_xz, g_z_zz_x_yy, g_z_zz_x_yz, g_z_zz_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_zz_0_xx[i] = 4.0 * g_z_zz_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_0_xy[i] = 4.0 * g_z_zz_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_0_xz[i] = 4.0 * g_z_zz_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_0_yy[i] = 4.0 * g_z_zz_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_0_yz[i] = 4.0 * g_z_zz_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_0_zz[i] = 4.0 * g_z_zz_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_z_0_y_0_0_xx_0_xx, g_z_0_y_0_0_xx_0_xy, g_z_0_y_0_0_xx_0_xz, g_z_0_y_0_0_xx_0_yy, g_z_0_y_0_0_xx_0_yz, g_z_0_y_0_0_xx_0_zz, g_z_xx_y_xx, g_z_xx_y_xy, g_z_xx_y_xz, g_z_xx_y_yy, g_z_xx_y_yz, g_z_xx_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xx_0_xx[i] = 4.0 * g_z_xx_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_0_xy[i] = 4.0 * g_z_xx_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_0_xz[i] = 4.0 * g_z_xx_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_0_yy[i] = 4.0 * g_z_xx_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_0_yz[i] = 4.0 * g_z_xx_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_0_zz[i] = 4.0 * g_z_xx_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_z_0_y_0_0_xy_0_xx, g_z_0_y_0_0_xy_0_xy, g_z_0_y_0_0_xy_0_xz, g_z_0_y_0_0_xy_0_yy, g_z_0_y_0_0_xy_0_yz, g_z_0_y_0_0_xy_0_zz, g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xy_0_xx[i] = 4.0 * g_z_xy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_0_xy[i] = 4.0 * g_z_xy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_0_xz[i] = 4.0 * g_z_xy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_0_yy[i] = 4.0 * g_z_xy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_0_yz[i] = 4.0 * g_z_xy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_0_zz[i] = 4.0 * g_z_xy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_z_0_y_0_0_xz_0_xx, g_z_0_y_0_0_xz_0_xy, g_z_0_y_0_0_xz_0_xz, g_z_0_y_0_0_xz_0_yy, g_z_0_y_0_0_xz_0_yz, g_z_0_y_0_0_xz_0_zz, g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xz_0_xx[i] = 4.0 * g_z_xz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_0_xy[i] = 4.0 * g_z_xz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_0_xz[i] = 4.0 * g_z_xz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_0_yy[i] = 4.0 * g_z_xz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_0_yz[i] = 4.0 * g_z_xz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_0_zz[i] = 4.0 * g_z_xz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_z_0_y_0_0_yy_0_xx, g_z_0_y_0_0_yy_0_xy, g_z_0_y_0_0_yy_0_xz, g_z_0_y_0_0_yy_0_yy, g_z_0_y_0_0_yy_0_yz, g_z_0_y_0_0_yy_0_zz, g_z_yy_y_xx, g_z_yy_y_xy, g_z_yy_y_xz, g_z_yy_y_yy, g_z_yy_y_yz, g_z_yy_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yy_0_xx[i] = 4.0 * g_z_yy_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_0_xy[i] = 4.0 * g_z_yy_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_0_xz[i] = 4.0 * g_z_yy_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_0_yy[i] = 4.0 * g_z_yy_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_0_yz[i] = 4.0 * g_z_yy_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_0_zz[i] = 4.0 * g_z_yy_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_z_0_y_0_0_yz_0_xx, g_z_0_y_0_0_yz_0_xy, g_z_0_y_0_0_yz_0_xz, g_z_0_y_0_0_yz_0_yy, g_z_0_y_0_0_yz_0_yz, g_z_0_y_0_0_yz_0_zz, g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yz_0_xx[i] = 4.0 * g_z_yz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_0_xy[i] = 4.0 * g_z_yz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_0_xz[i] = 4.0 * g_z_yz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_0_yy[i] = 4.0 * g_z_yz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_0_yz[i] = 4.0 * g_z_yz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_0_zz[i] = 4.0 * g_z_yz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_z_0_y_0_0_zz_0_xx, g_z_0_y_0_0_zz_0_xy, g_z_0_y_0_0_zz_0_xz, g_z_0_y_0_0_zz_0_yy, g_z_0_y_0_0_zz_0_yz, g_z_0_y_0_0_zz_0_zz, g_z_zz_y_xx, g_z_zz_y_xy, g_z_zz_y_xz, g_z_zz_y_yy, g_z_zz_y_yz, g_z_zz_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_zz_0_xx[i] = 4.0 * g_z_zz_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_0_xy[i] = 4.0 * g_z_zz_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_0_xz[i] = 4.0 * g_z_zz_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_0_yy[i] = 4.0 * g_z_zz_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_0_yz[i] = 4.0 * g_z_zz_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_0_zz[i] = 4.0 * g_z_zz_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_z_0_z_0_0_xx_0_xx, g_z_0_z_0_0_xx_0_xy, g_z_0_z_0_0_xx_0_xz, g_z_0_z_0_0_xx_0_yy, g_z_0_z_0_0_xx_0_yz, g_z_0_z_0_0_xx_0_zz, g_z_xx_z_xx, g_z_xx_z_xy, g_z_xx_z_xz, g_z_xx_z_yy, g_z_xx_z_yz, g_z_xx_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xx_0_xx[i] = 4.0 * g_z_xx_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_0_xy[i] = 4.0 * g_z_xx_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_0_xz[i] = 4.0 * g_z_xx_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_0_yy[i] = 4.0 * g_z_xx_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_0_yz[i] = 4.0 * g_z_xx_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_0_zz[i] = 4.0 * g_z_xx_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_z_0_z_0_0_xy_0_xx, g_z_0_z_0_0_xy_0_xy, g_z_0_z_0_0_xy_0_xz, g_z_0_z_0_0_xy_0_yy, g_z_0_z_0_0_xy_0_yz, g_z_0_z_0_0_xy_0_zz, g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xy_0_xx[i] = 4.0 * g_z_xy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_0_xy[i] = 4.0 * g_z_xy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_0_xz[i] = 4.0 * g_z_xy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_0_yy[i] = 4.0 * g_z_xy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_0_yz[i] = 4.0 * g_z_xy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_0_zz[i] = 4.0 * g_z_xy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_z_0_z_0_0_xz_0_xx, g_z_0_z_0_0_xz_0_xy, g_z_0_z_0_0_xz_0_xz, g_z_0_z_0_0_xz_0_yy, g_z_0_z_0_0_xz_0_yz, g_z_0_z_0_0_xz_0_zz, g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xz_0_xx[i] = 4.0 * g_z_xz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_0_xy[i] = 4.0 * g_z_xz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_0_xz[i] = 4.0 * g_z_xz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_0_yy[i] = 4.0 * g_z_xz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_0_yz[i] = 4.0 * g_z_xz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_0_zz[i] = 4.0 * g_z_xz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_z_0_z_0_0_yy_0_xx, g_z_0_z_0_0_yy_0_xy, g_z_0_z_0_0_yy_0_xz, g_z_0_z_0_0_yy_0_yy, g_z_0_z_0_0_yy_0_yz, g_z_0_z_0_0_yy_0_zz, g_z_yy_z_xx, g_z_yy_z_xy, g_z_yy_z_xz, g_z_yy_z_yy, g_z_yy_z_yz, g_z_yy_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yy_0_xx[i] = 4.0 * g_z_yy_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_0_xy[i] = 4.0 * g_z_yy_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_0_xz[i] = 4.0 * g_z_yy_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_0_yy[i] = 4.0 * g_z_yy_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_0_yz[i] = 4.0 * g_z_yy_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_0_zz[i] = 4.0 * g_z_yy_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_z_0_z_0_0_yz_0_xx, g_z_0_z_0_0_yz_0_xy, g_z_0_z_0_0_yz_0_xz, g_z_0_z_0_0_yz_0_yy, g_z_0_z_0_0_yz_0_yz, g_z_0_z_0_0_yz_0_zz, g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yz_0_xx[i] = 4.0 * g_z_yz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_0_xy[i] = 4.0 * g_z_yz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_0_xz[i] = 4.0 * g_z_yz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_0_yy[i] = 4.0 * g_z_yz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_0_yz[i] = 4.0 * g_z_yz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_0_zz[i] = 4.0 * g_z_yz_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_z_0_z_0_0_zz_0_xx, g_z_0_z_0_0_zz_0_xy, g_z_0_z_0_0_zz_0_xz, g_z_0_z_0_0_zz_0_yy, g_z_0_z_0_0_zz_0_yz, g_z_0_z_0_0_zz_0_zz, g_z_zz_z_xx, g_z_zz_z_xy, g_z_zz_z_xz, g_z_zz_z_yy, g_z_zz_z_yz, g_z_zz_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_zz_0_xx[i] = 4.0 * g_z_zz_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_0_xy[i] = 4.0 * g_z_zz_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_0_xz[i] = 4.0 * g_z_zz_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_0_yy[i] = 4.0 * g_z_zz_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_0_yz[i] = 4.0 * g_z_zz_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_0_zz[i] = 4.0 * g_z_zz_z_zz[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace

