#include "GeomDeriv1100OfScalarForSDDD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_sddd_0(CSimdArray<double>& buffer_1100_sddd,
                     const CSimdArray<double>& buffer_ppdd,
                     const CSimdArray<double>& buffer_pfdd,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_sddd.number_of_columns();

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

    /// Set up components of auxilary buffer : buffer_pfdd

    auto g_x_xxx_xx_xx = buffer_pfdd[0];

    auto g_x_xxx_xx_xy = buffer_pfdd[1];

    auto g_x_xxx_xx_xz = buffer_pfdd[2];

    auto g_x_xxx_xx_yy = buffer_pfdd[3];

    auto g_x_xxx_xx_yz = buffer_pfdd[4];

    auto g_x_xxx_xx_zz = buffer_pfdd[5];

    auto g_x_xxx_xy_xx = buffer_pfdd[6];

    auto g_x_xxx_xy_xy = buffer_pfdd[7];

    auto g_x_xxx_xy_xz = buffer_pfdd[8];

    auto g_x_xxx_xy_yy = buffer_pfdd[9];

    auto g_x_xxx_xy_yz = buffer_pfdd[10];

    auto g_x_xxx_xy_zz = buffer_pfdd[11];

    auto g_x_xxx_xz_xx = buffer_pfdd[12];

    auto g_x_xxx_xz_xy = buffer_pfdd[13];

    auto g_x_xxx_xz_xz = buffer_pfdd[14];

    auto g_x_xxx_xz_yy = buffer_pfdd[15];

    auto g_x_xxx_xz_yz = buffer_pfdd[16];

    auto g_x_xxx_xz_zz = buffer_pfdd[17];

    auto g_x_xxx_yy_xx = buffer_pfdd[18];

    auto g_x_xxx_yy_xy = buffer_pfdd[19];

    auto g_x_xxx_yy_xz = buffer_pfdd[20];

    auto g_x_xxx_yy_yy = buffer_pfdd[21];

    auto g_x_xxx_yy_yz = buffer_pfdd[22];

    auto g_x_xxx_yy_zz = buffer_pfdd[23];

    auto g_x_xxx_yz_xx = buffer_pfdd[24];

    auto g_x_xxx_yz_xy = buffer_pfdd[25];

    auto g_x_xxx_yz_xz = buffer_pfdd[26];

    auto g_x_xxx_yz_yy = buffer_pfdd[27];

    auto g_x_xxx_yz_yz = buffer_pfdd[28];

    auto g_x_xxx_yz_zz = buffer_pfdd[29];

    auto g_x_xxx_zz_xx = buffer_pfdd[30];

    auto g_x_xxx_zz_xy = buffer_pfdd[31];

    auto g_x_xxx_zz_xz = buffer_pfdd[32];

    auto g_x_xxx_zz_yy = buffer_pfdd[33];

    auto g_x_xxx_zz_yz = buffer_pfdd[34];

    auto g_x_xxx_zz_zz = buffer_pfdd[35];

    auto g_x_xxy_xx_xx = buffer_pfdd[36];

    auto g_x_xxy_xx_xy = buffer_pfdd[37];

    auto g_x_xxy_xx_xz = buffer_pfdd[38];

    auto g_x_xxy_xx_yy = buffer_pfdd[39];

    auto g_x_xxy_xx_yz = buffer_pfdd[40];

    auto g_x_xxy_xx_zz = buffer_pfdd[41];

    auto g_x_xxy_xy_xx = buffer_pfdd[42];

    auto g_x_xxy_xy_xy = buffer_pfdd[43];

    auto g_x_xxy_xy_xz = buffer_pfdd[44];

    auto g_x_xxy_xy_yy = buffer_pfdd[45];

    auto g_x_xxy_xy_yz = buffer_pfdd[46];

    auto g_x_xxy_xy_zz = buffer_pfdd[47];

    auto g_x_xxy_xz_xx = buffer_pfdd[48];

    auto g_x_xxy_xz_xy = buffer_pfdd[49];

    auto g_x_xxy_xz_xz = buffer_pfdd[50];

    auto g_x_xxy_xz_yy = buffer_pfdd[51];

    auto g_x_xxy_xz_yz = buffer_pfdd[52];

    auto g_x_xxy_xz_zz = buffer_pfdd[53];

    auto g_x_xxy_yy_xx = buffer_pfdd[54];

    auto g_x_xxy_yy_xy = buffer_pfdd[55];

    auto g_x_xxy_yy_xz = buffer_pfdd[56];

    auto g_x_xxy_yy_yy = buffer_pfdd[57];

    auto g_x_xxy_yy_yz = buffer_pfdd[58];

    auto g_x_xxy_yy_zz = buffer_pfdd[59];

    auto g_x_xxy_yz_xx = buffer_pfdd[60];

    auto g_x_xxy_yz_xy = buffer_pfdd[61];

    auto g_x_xxy_yz_xz = buffer_pfdd[62];

    auto g_x_xxy_yz_yy = buffer_pfdd[63];

    auto g_x_xxy_yz_yz = buffer_pfdd[64];

    auto g_x_xxy_yz_zz = buffer_pfdd[65];

    auto g_x_xxy_zz_xx = buffer_pfdd[66];

    auto g_x_xxy_zz_xy = buffer_pfdd[67];

    auto g_x_xxy_zz_xz = buffer_pfdd[68];

    auto g_x_xxy_zz_yy = buffer_pfdd[69];

    auto g_x_xxy_zz_yz = buffer_pfdd[70];

    auto g_x_xxy_zz_zz = buffer_pfdd[71];

    auto g_x_xxz_xx_xx = buffer_pfdd[72];

    auto g_x_xxz_xx_xy = buffer_pfdd[73];

    auto g_x_xxz_xx_xz = buffer_pfdd[74];

    auto g_x_xxz_xx_yy = buffer_pfdd[75];

    auto g_x_xxz_xx_yz = buffer_pfdd[76];

    auto g_x_xxz_xx_zz = buffer_pfdd[77];

    auto g_x_xxz_xy_xx = buffer_pfdd[78];

    auto g_x_xxz_xy_xy = buffer_pfdd[79];

    auto g_x_xxz_xy_xz = buffer_pfdd[80];

    auto g_x_xxz_xy_yy = buffer_pfdd[81];

    auto g_x_xxz_xy_yz = buffer_pfdd[82];

    auto g_x_xxz_xy_zz = buffer_pfdd[83];

    auto g_x_xxz_xz_xx = buffer_pfdd[84];

    auto g_x_xxz_xz_xy = buffer_pfdd[85];

    auto g_x_xxz_xz_xz = buffer_pfdd[86];

    auto g_x_xxz_xz_yy = buffer_pfdd[87];

    auto g_x_xxz_xz_yz = buffer_pfdd[88];

    auto g_x_xxz_xz_zz = buffer_pfdd[89];

    auto g_x_xxz_yy_xx = buffer_pfdd[90];

    auto g_x_xxz_yy_xy = buffer_pfdd[91];

    auto g_x_xxz_yy_xz = buffer_pfdd[92];

    auto g_x_xxz_yy_yy = buffer_pfdd[93];

    auto g_x_xxz_yy_yz = buffer_pfdd[94];

    auto g_x_xxz_yy_zz = buffer_pfdd[95];

    auto g_x_xxz_yz_xx = buffer_pfdd[96];

    auto g_x_xxz_yz_xy = buffer_pfdd[97];

    auto g_x_xxz_yz_xz = buffer_pfdd[98];

    auto g_x_xxz_yz_yy = buffer_pfdd[99];

    auto g_x_xxz_yz_yz = buffer_pfdd[100];

    auto g_x_xxz_yz_zz = buffer_pfdd[101];

    auto g_x_xxz_zz_xx = buffer_pfdd[102];

    auto g_x_xxz_zz_xy = buffer_pfdd[103];

    auto g_x_xxz_zz_xz = buffer_pfdd[104];

    auto g_x_xxz_zz_yy = buffer_pfdd[105];

    auto g_x_xxz_zz_yz = buffer_pfdd[106];

    auto g_x_xxz_zz_zz = buffer_pfdd[107];

    auto g_x_xyy_xx_xx = buffer_pfdd[108];

    auto g_x_xyy_xx_xy = buffer_pfdd[109];

    auto g_x_xyy_xx_xz = buffer_pfdd[110];

    auto g_x_xyy_xx_yy = buffer_pfdd[111];

    auto g_x_xyy_xx_yz = buffer_pfdd[112];

    auto g_x_xyy_xx_zz = buffer_pfdd[113];

    auto g_x_xyy_xy_xx = buffer_pfdd[114];

    auto g_x_xyy_xy_xy = buffer_pfdd[115];

    auto g_x_xyy_xy_xz = buffer_pfdd[116];

    auto g_x_xyy_xy_yy = buffer_pfdd[117];

    auto g_x_xyy_xy_yz = buffer_pfdd[118];

    auto g_x_xyy_xy_zz = buffer_pfdd[119];

    auto g_x_xyy_xz_xx = buffer_pfdd[120];

    auto g_x_xyy_xz_xy = buffer_pfdd[121];

    auto g_x_xyy_xz_xz = buffer_pfdd[122];

    auto g_x_xyy_xz_yy = buffer_pfdd[123];

    auto g_x_xyy_xz_yz = buffer_pfdd[124];

    auto g_x_xyy_xz_zz = buffer_pfdd[125];

    auto g_x_xyy_yy_xx = buffer_pfdd[126];

    auto g_x_xyy_yy_xy = buffer_pfdd[127];

    auto g_x_xyy_yy_xz = buffer_pfdd[128];

    auto g_x_xyy_yy_yy = buffer_pfdd[129];

    auto g_x_xyy_yy_yz = buffer_pfdd[130];

    auto g_x_xyy_yy_zz = buffer_pfdd[131];

    auto g_x_xyy_yz_xx = buffer_pfdd[132];

    auto g_x_xyy_yz_xy = buffer_pfdd[133];

    auto g_x_xyy_yz_xz = buffer_pfdd[134];

    auto g_x_xyy_yz_yy = buffer_pfdd[135];

    auto g_x_xyy_yz_yz = buffer_pfdd[136];

    auto g_x_xyy_yz_zz = buffer_pfdd[137];

    auto g_x_xyy_zz_xx = buffer_pfdd[138];

    auto g_x_xyy_zz_xy = buffer_pfdd[139];

    auto g_x_xyy_zz_xz = buffer_pfdd[140];

    auto g_x_xyy_zz_yy = buffer_pfdd[141];

    auto g_x_xyy_zz_yz = buffer_pfdd[142];

    auto g_x_xyy_zz_zz = buffer_pfdd[143];

    auto g_x_xyz_xx_xx = buffer_pfdd[144];

    auto g_x_xyz_xx_xy = buffer_pfdd[145];

    auto g_x_xyz_xx_xz = buffer_pfdd[146];

    auto g_x_xyz_xx_yy = buffer_pfdd[147];

    auto g_x_xyz_xx_yz = buffer_pfdd[148];

    auto g_x_xyz_xx_zz = buffer_pfdd[149];

    auto g_x_xyz_xy_xx = buffer_pfdd[150];

    auto g_x_xyz_xy_xy = buffer_pfdd[151];

    auto g_x_xyz_xy_xz = buffer_pfdd[152];

    auto g_x_xyz_xy_yy = buffer_pfdd[153];

    auto g_x_xyz_xy_yz = buffer_pfdd[154];

    auto g_x_xyz_xy_zz = buffer_pfdd[155];

    auto g_x_xyz_xz_xx = buffer_pfdd[156];

    auto g_x_xyz_xz_xy = buffer_pfdd[157];

    auto g_x_xyz_xz_xz = buffer_pfdd[158];

    auto g_x_xyz_xz_yy = buffer_pfdd[159];

    auto g_x_xyz_xz_yz = buffer_pfdd[160];

    auto g_x_xyz_xz_zz = buffer_pfdd[161];

    auto g_x_xyz_yy_xx = buffer_pfdd[162];

    auto g_x_xyz_yy_xy = buffer_pfdd[163];

    auto g_x_xyz_yy_xz = buffer_pfdd[164];

    auto g_x_xyz_yy_yy = buffer_pfdd[165];

    auto g_x_xyz_yy_yz = buffer_pfdd[166];

    auto g_x_xyz_yy_zz = buffer_pfdd[167];

    auto g_x_xyz_yz_xx = buffer_pfdd[168];

    auto g_x_xyz_yz_xy = buffer_pfdd[169];

    auto g_x_xyz_yz_xz = buffer_pfdd[170];

    auto g_x_xyz_yz_yy = buffer_pfdd[171];

    auto g_x_xyz_yz_yz = buffer_pfdd[172];

    auto g_x_xyz_yz_zz = buffer_pfdd[173];

    auto g_x_xyz_zz_xx = buffer_pfdd[174];

    auto g_x_xyz_zz_xy = buffer_pfdd[175];

    auto g_x_xyz_zz_xz = buffer_pfdd[176];

    auto g_x_xyz_zz_yy = buffer_pfdd[177];

    auto g_x_xyz_zz_yz = buffer_pfdd[178];

    auto g_x_xyz_zz_zz = buffer_pfdd[179];

    auto g_x_xzz_xx_xx = buffer_pfdd[180];

    auto g_x_xzz_xx_xy = buffer_pfdd[181];

    auto g_x_xzz_xx_xz = buffer_pfdd[182];

    auto g_x_xzz_xx_yy = buffer_pfdd[183];

    auto g_x_xzz_xx_yz = buffer_pfdd[184];

    auto g_x_xzz_xx_zz = buffer_pfdd[185];

    auto g_x_xzz_xy_xx = buffer_pfdd[186];

    auto g_x_xzz_xy_xy = buffer_pfdd[187];

    auto g_x_xzz_xy_xz = buffer_pfdd[188];

    auto g_x_xzz_xy_yy = buffer_pfdd[189];

    auto g_x_xzz_xy_yz = buffer_pfdd[190];

    auto g_x_xzz_xy_zz = buffer_pfdd[191];

    auto g_x_xzz_xz_xx = buffer_pfdd[192];

    auto g_x_xzz_xz_xy = buffer_pfdd[193];

    auto g_x_xzz_xz_xz = buffer_pfdd[194];

    auto g_x_xzz_xz_yy = buffer_pfdd[195];

    auto g_x_xzz_xz_yz = buffer_pfdd[196];

    auto g_x_xzz_xz_zz = buffer_pfdd[197];

    auto g_x_xzz_yy_xx = buffer_pfdd[198];

    auto g_x_xzz_yy_xy = buffer_pfdd[199];

    auto g_x_xzz_yy_xz = buffer_pfdd[200];

    auto g_x_xzz_yy_yy = buffer_pfdd[201];

    auto g_x_xzz_yy_yz = buffer_pfdd[202];

    auto g_x_xzz_yy_zz = buffer_pfdd[203];

    auto g_x_xzz_yz_xx = buffer_pfdd[204];

    auto g_x_xzz_yz_xy = buffer_pfdd[205];

    auto g_x_xzz_yz_xz = buffer_pfdd[206];

    auto g_x_xzz_yz_yy = buffer_pfdd[207];

    auto g_x_xzz_yz_yz = buffer_pfdd[208];

    auto g_x_xzz_yz_zz = buffer_pfdd[209];

    auto g_x_xzz_zz_xx = buffer_pfdd[210];

    auto g_x_xzz_zz_xy = buffer_pfdd[211];

    auto g_x_xzz_zz_xz = buffer_pfdd[212];

    auto g_x_xzz_zz_yy = buffer_pfdd[213];

    auto g_x_xzz_zz_yz = buffer_pfdd[214];

    auto g_x_xzz_zz_zz = buffer_pfdd[215];

    auto g_x_yyy_xx_xx = buffer_pfdd[216];

    auto g_x_yyy_xx_xy = buffer_pfdd[217];

    auto g_x_yyy_xx_xz = buffer_pfdd[218];

    auto g_x_yyy_xx_yy = buffer_pfdd[219];

    auto g_x_yyy_xx_yz = buffer_pfdd[220];

    auto g_x_yyy_xx_zz = buffer_pfdd[221];

    auto g_x_yyy_xy_xx = buffer_pfdd[222];

    auto g_x_yyy_xy_xy = buffer_pfdd[223];

    auto g_x_yyy_xy_xz = buffer_pfdd[224];

    auto g_x_yyy_xy_yy = buffer_pfdd[225];

    auto g_x_yyy_xy_yz = buffer_pfdd[226];

    auto g_x_yyy_xy_zz = buffer_pfdd[227];

    auto g_x_yyy_xz_xx = buffer_pfdd[228];

    auto g_x_yyy_xz_xy = buffer_pfdd[229];

    auto g_x_yyy_xz_xz = buffer_pfdd[230];

    auto g_x_yyy_xz_yy = buffer_pfdd[231];

    auto g_x_yyy_xz_yz = buffer_pfdd[232];

    auto g_x_yyy_xz_zz = buffer_pfdd[233];

    auto g_x_yyy_yy_xx = buffer_pfdd[234];

    auto g_x_yyy_yy_xy = buffer_pfdd[235];

    auto g_x_yyy_yy_xz = buffer_pfdd[236];

    auto g_x_yyy_yy_yy = buffer_pfdd[237];

    auto g_x_yyy_yy_yz = buffer_pfdd[238];

    auto g_x_yyy_yy_zz = buffer_pfdd[239];

    auto g_x_yyy_yz_xx = buffer_pfdd[240];

    auto g_x_yyy_yz_xy = buffer_pfdd[241];

    auto g_x_yyy_yz_xz = buffer_pfdd[242];

    auto g_x_yyy_yz_yy = buffer_pfdd[243];

    auto g_x_yyy_yz_yz = buffer_pfdd[244];

    auto g_x_yyy_yz_zz = buffer_pfdd[245];

    auto g_x_yyy_zz_xx = buffer_pfdd[246];

    auto g_x_yyy_zz_xy = buffer_pfdd[247];

    auto g_x_yyy_zz_xz = buffer_pfdd[248];

    auto g_x_yyy_zz_yy = buffer_pfdd[249];

    auto g_x_yyy_zz_yz = buffer_pfdd[250];

    auto g_x_yyy_zz_zz = buffer_pfdd[251];

    auto g_x_yyz_xx_xx = buffer_pfdd[252];

    auto g_x_yyz_xx_xy = buffer_pfdd[253];

    auto g_x_yyz_xx_xz = buffer_pfdd[254];

    auto g_x_yyz_xx_yy = buffer_pfdd[255];

    auto g_x_yyz_xx_yz = buffer_pfdd[256];

    auto g_x_yyz_xx_zz = buffer_pfdd[257];

    auto g_x_yyz_xy_xx = buffer_pfdd[258];

    auto g_x_yyz_xy_xy = buffer_pfdd[259];

    auto g_x_yyz_xy_xz = buffer_pfdd[260];

    auto g_x_yyz_xy_yy = buffer_pfdd[261];

    auto g_x_yyz_xy_yz = buffer_pfdd[262];

    auto g_x_yyz_xy_zz = buffer_pfdd[263];

    auto g_x_yyz_xz_xx = buffer_pfdd[264];

    auto g_x_yyz_xz_xy = buffer_pfdd[265];

    auto g_x_yyz_xz_xz = buffer_pfdd[266];

    auto g_x_yyz_xz_yy = buffer_pfdd[267];

    auto g_x_yyz_xz_yz = buffer_pfdd[268];

    auto g_x_yyz_xz_zz = buffer_pfdd[269];

    auto g_x_yyz_yy_xx = buffer_pfdd[270];

    auto g_x_yyz_yy_xy = buffer_pfdd[271];

    auto g_x_yyz_yy_xz = buffer_pfdd[272];

    auto g_x_yyz_yy_yy = buffer_pfdd[273];

    auto g_x_yyz_yy_yz = buffer_pfdd[274];

    auto g_x_yyz_yy_zz = buffer_pfdd[275];

    auto g_x_yyz_yz_xx = buffer_pfdd[276];

    auto g_x_yyz_yz_xy = buffer_pfdd[277];

    auto g_x_yyz_yz_xz = buffer_pfdd[278];

    auto g_x_yyz_yz_yy = buffer_pfdd[279];

    auto g_x_yyz_yz_yz = buffer_pfdd[280];

    auto g_x_yyz_yz_zz = buffer_pfdd[281];

    auto g_x_yyz_zz_xx = buffer_pfdd[282];

    auto g_x_yyz_zz_xy = buffer_pfdd[283];

    auto g_x_yyz_zz_xz = buffer_pfdd[284];

    auto g_x_yyz_zz_yy = buffer_pfdd[285];

    auto g_x_yyz_zz_yz = buffer_pfdd[286];

    auto g_x_yyz_zz_zz = buffer_pfdd[287];

    auto g_x_yzz_xx_xx = buffer_pfdd[288];

    auto g_x_yzz_xx_xy = buffer_pfdd[289];

    auto g_x_yzz_xx_xz = buffer_pfdd[290];

    auto g_x_yzz_xx_yy = buffer_pfdd[291];

    auto g_x_yzz_xx_yz = buffer_pfdd[292];

    auto g_x_yzz_xx_zz = buffer_pfdd[293];

    auto g_x_yzz_xy_xx = buffer_pfdd[294];

    auto g_x_yzz_xy_xy = buffer_pfdd[295];

    auto g_x_yzz_xy_xz = buffer_pfdd[296];

    auto g_x_yzz_xy_yy = buffer_pfdd[297];

    auto g_x_yzz_xy_yz = buffer_pfdd[298];

    auto g_x_yzz_xy_zz = buffer_pfdd[299];

    auto g_x_yzz_xz_xx = buffer_pfdd[300];

    auto g_x_yzz_xz_xy = buffer_pfdd[301];

    auto g_x_yzz_xz_xz = buffer_pfdd[302];

    auto g_x_yzz_xz_yy = buffer_pfdd[303];

    auto g_x_yzz_xz_yz = buffer_pfdd[304];

    auto g_x_yzz_xz_zz = buffer_pfdd[305];

    auto g_x_yzz_yy_xx = buffer_pfdd[306];

    auto g_x_yzz_yy_xy = buffer_pfdd[307];

    auto g_x_yzz_yy_xz = buffer_pfdd[308];

    auto g_x_yzz_yy_yy = buffer_pfdd[309];

    auto g_x_yzz_yy_yz = buffer_pfdd[310];

    auto g_x_yzz_yy_zz = buffer_pfdd[311];

    auto g_x_yzz_yz_xx = buffer_pfdd[312];

    auto g_x_yzz_yz_xy = buffer_pfdd[313];

    auto g_x_yzz_yz_xz = buffer_pfdd[314];

    auto g_x_yzz_yz_yy = buffer_pfdd[315];

    auto g_x_yzz_yz_yz = buffer_pfdd[316];

    auto g_x_yzz_yz_zz = buffer_pfdd[317];

    auto g_x_yzz_zz_xx = buffer_pfdd[318];

    auto g_x_yzz_zz_xy = buffer_pfdd[319];

    auto g_x_yzz_zz_xz = buffer_pfdd[320];

    auto g_x_yzz_zz_yy = buffer_pfdd[321];

    auto g_x_yzz_zz_yz = buffer_pfdd[322];

    auto g_x_yzz_zz_zz = buffer_pfdd[323];

    auto g_x_zzz_xx_xx = buffer_pfdd[324];

    auto g_x_zzz_xx_xy = buffer_pfdd[325];

    auto g_x_zzz_xx_xz = buffer_pfdd[326];

    auto g_x_zzz_xx_yy = buffer_pfdd[327];

    auto g_x_zzz_xx_yz = buffer_pfdd[328];

    auto g_x_zzz_xx_zz = buffer_pfdd[329];

    auto g_x_zzz_xy_xx = buffer_pfdd[330];

    auto g_x_zzz_xy_xy = buffer_pfdd[331];

    auto g_x_zzz_xy_xz = buffer_pfdd[332];

    auto g_x_zzz_xy_yy = buffer_pfdd[333];

    auto g_x_zzz_xy_yz = buffer_pfdd[334];

    auto g_x_zzz_xy_zz = buffer_pfdd[335];

    auto g_x_zzz_xz_xx = buffer_pfdd[336];

    auto g_x_zzz_xz_xy = buffer_pfdd[337];

    auto g_x_zzz_xz_xz = buffer_pfdd[338];

    auto g_x_zzz_xz_yy = buffer_pfdd[339];

    auto g_x_zzz_xz_yz = buffer_pfdd[340];

    auto g_x_zzz_xz_zz = buffer_pfdd[341];

    auto g_x_zzz_yy_xx = buffer_pfdd[342];

    auto g_x_zzz_yy_xy = buffer_pfdd[343];

    auto g_x_zzz_yy_xz = buffer_pfdd[344];

    auto g_x_zzz_yy_yy = buffer_pfdd[345];

    auto g_x_zzz_yy_yz = buffer_pfdd[346];

    auto g_x_zzz_yy_zz = buffer_pfdd[347];

    auto g_x_zzz_yz_xx = buffer_pfdd[348];

    auto g_x_zzz_yz_xy = buffer_pfdd[349];

    auto g_x_zzz_yz_xz = buffer_pfdd[350];

    auto g_x_zzz_yz_yy = buffer_pfdd[351];

    auto g_x_zzz_yz_yz = buffer_pfdd[352];

    auto g_x_zzz_yz_zz = buffer_pfdd[353];

    auto g_x_zzz_zz_xx = buffer_pfdd[354];

    auto g_x_zzz_zz_xy = buffer_pfdd[355];

    auto g_x_zzz_zz_xz = buffer_pfdd[356];

    auto g_x_zzz_zz_yy = buffer_pfdd[357];

    auto g_x_zzz_zz_yz = buffer_pfdd[358];

    auto g_x_zzz_zz_zz = buffer_pfdd[359];

    auto g_y_xxx_xx_xx = buffer_pfdd[360];

    auto g_y_xxx_xx_xy = buffer_pfdd[361];

    auto g_y_xxx_xx_xz = buffer_pfdd[362];

    auto g_y_xxx_xx_yy = buffer_pfdd[363];

    auto g_y_xxx_xx_yz = buffer_pfdd[364];

    auto g_y_xxx_xx_zz = buffer_pfdd[365];

    auto g_y_xxx_xy_xx = buffer_pfdd[366];

    auto g_y_xxx_xy_xy = buffer_pfdd[367];

    auto g_y_xxx_xy_xz = buffer_pfdd[368];

    auto g_y_xxx_xy_yy = buffer_pfdd[369];

    auto g_y_xxx_xy_yz = buffer_pfdd[370];

    auto g_y_xxx_xy_zz = buffer_pfdd[371];

    auto g_y_xxx_xz_xx = buffer_pfdd[372];

    auto g_y_xxx_xz_xy = buffer_pfdd[373];

    auto g_y_xxx_xz_xz = buffer_pfdd[374];

    auto g_y_xxx_xz_yy = buffer_pfdd[375];

    auto g_y_xxx_xz_yz = buffer_pfdd[376];

    auto g_y_xxx_xz_zz = buffer_pfdd[377];

    auto g_y_xxx_yy_xx = buffer_pfdd[378];

    auto g_y_xxx_yy_xy = buffer_pfdd[379];

    auto g_y_xxx_yy_xz = buffer_pfdd[380];

    auto g_y_xxx_yy_yy = buffer_pfdd[381];

    auto g_y_xxx_yy_yz = buffer_pfdd[382];

    auto g_y_xxx_yy_zz = buffer_pfdd[383];

    auto g_y_xxx_yz_xx = buffer_pfdd[384];

    auto g_y_xxx_yz_xy = buffer_pfdd[385];

    auto g_y_xxx_yz_xz = buffer_pfdd[386];

    auto g_y_xxx_yz_yy = buffer_pfdd[387];

    auto g_y_xxx_yz_yz = buffer_pfdd[388];

    auto g_y_xxx_yz_zz = buffer_pfdd[389];

    auto g_y_xxx_zz_xx = buffer_pfdd[390];

    auto g_y_xxx_zz_xy = buffer_pfdd[391];

    auto g_y_xxx_zz_xz = buffer_pfdd[392];

    auto g_y_xxx_zz_yy = buffer_pfdd[393];

    auto g_y_xxx_zz_yz = buffer_pfdd[394];

    auto g_y_xxx_zz_zz = buffer_pfdd[395];

    auto g_y_xxy_xx_xx = buffer_pfdd[396];

    auto g_y_xxy_xx_xy = buffer_pfdd[397];

    auto g_y_xxy_xx_xz = buffer_pfdd[398];

    auto g_y_xxy_xx_yy = buffer_pfdd[399];

    auto g_y_xxy_xx_yz = buffer_pfdd[400];

    auto g_y_xxy_xx_zz = buffer_pfdd[401];

    auto g_y_xxy_xy_xx = buffer_pfdd[402];

    auto g_y_xxy_xy_xy = buffer_pfdd[403];

    auto g_y_xxy_xy_xz = buffer_pfdd[404];

    auto g_y_xxy_xy_yy = buffer_pfdd[405];

    auto g_y_xxy_xy_yz = buffer_pfdd[406];

    auto g_y_xxy_xy_zz = buffer_pfdd[407];

    auto g_y_xxy_xz_xx = buffer_pfdd[408];

    auto g_y_xxy_xz_xy = buffer_pfdd[409];

    auto g_y_xxy_xz_xz = buffer_pfdd[410];

    auto g_y_xxy_xz_yy = buffer_pfdd[411];

    auto g_y_xxy_xz_yz = buffer_pfdd[412];

    auto g_y_xxy_xz_zz = buffer_pfdd[413];

    auto g_y_xxy_yy_xx = buffer_pfdd[414];

    auto g_y_xxy_yy_xy = buffer_pfdd[415];

    auto g_y_xxy_yy_xz = buffer_pfdd[416];

    auto g_y_xxy_yy_yy = buffer_pfdd[417];

    auto g_y_xxy_yy_yz = buffer_pfdd[418];

    auto g_y_xxy_yy_zz = buffer_pfdd[419];

    auto g_y_xxy_yz_xx = buffer_pfdd[420];

    auto g_y_xxy_yz_xy = buffer_pfdd[421];

    auto g_y_xxy_yz_xz = buffer_pfdd[422];

    auto g_y_xxy_yz_yy = buffer_pfdd[423];

    auto g_y_xxy_yz_yz = buffer_pfdd[424];

    auto g_y_xxy_yz_zz = buffer_pfdd[425];

    auto g_y_xxy_zz_xx = buffer_pfdd[426];

    auto g_y_xxy_zz_xy = buffer_pfdd[427];

    auto g_y_xxy_zz_xz = buffer_pfdd[428];

    auto g_y_xxy_zz_yy = buffer_pfdd[429];

    auto g_y_xxy_zz_yz = buffer_pfdd[430];

    auto g_y_xxy_zz_zz = buffer_pfdd[431];

    auto g_y_xxz_xx_xx = buffer_pfdd[432];

    auto g_y_xxz_xx_xy = buffer_pfdd[433];

    auto g_y_xxz_xx_xz = buffer_pfdd[434];

    auto g_y_xxz_xx_yy = buffer_pfdd[435];

    auto g_y_xxz_xx_yz = buffer_pfdd[436];

    auto g_y_xxz_xx_zz = buffer_pfdd[437];

    auto g_y_xxz_xy_xx = buffer_pfdd[438];

    auto g_y_xxz_xy_xy = buffer_pfdd[439];

    auto g_y_xxz_xy_xz = buffer_pfdd[440];

    auto g_y_xxz_xy_yy = buffer_pfdd[441];

    auto g_y_xxz_xy_yz = buffer_pfdd[442];

    auto g_y_xxz_xy_zz = buffer_pfdd[443];

    auto g_y_xxz_xz_xx = buffer_pfdd[444];

    auto g_y_xxz_xz_xy = buffer_pfdd[445];

    auto g_y_xxz_xz_xz = buffer_pfdd[446];

    auto g_y_xxz_xz_yy = buffer_pfdd[447];

    auto g_y_xxz_xz_yz = buffer_pfdd[448];

    auto g_y_xxz_xz_zz = buffer_pfdd[449];

    auto g_y_xxz_yy_xx = buffer_pfdd[450];

    auto g_y_xxz_yy_xy = buffer_pfdd[451];

    auto g_y_xxz_yy_xz = buffer_pfdd[452];

    auto g_y_xxz_yy_yy = buffer_pfdd[453];

    auto g_y_xxz_yy_yz = buffer_pfdd[454];

    auto g_y_xxz_yy_zz = buffer_pfdd[455];

    auto g_y_xxz_yz_xx = buffer_pfdd[456];

    auto g_y_xxz_yz_xy = buffer_pfdd[457];

    auto g_y_xxz_yz_xz = buffer_pfdd[458];

    auto g_y_xxz_yz_yy = buffer_pfdd[459];

    auto g_y_xxz_yz_yz = buffer_pfdd[460];

    auto g_y_xxz_yz_zz = buffer_pfdd[461];

    auto g_y_xxz_zz_xx = buffer_pfdd[462];

    auto g_y_xxz_zz_xy = buffer_pfdd[463];

    auto g_y_xxz_zz_xz = buffer_pfdd[464];

    auto g_y_xxz_zz_yy = buffer_pfdd[465];

    auto g_y_xxz_zz_yz = buffer_pfdd[466];

    auto g_y_xxz_zz_zz = buffer_pfdd[467];

    auto g_y_xyy_xx_xx = buffer_pfdd[468];

    auto g_y_xyy_xx_xy = buffer_pfdd[469];

    auto g_y_xyy_xx_xz = buffer_pfdd[470];

    auto g_y_xyy_xx_yy = buffer_pfdd[471];

    auto g_y_xyy_xx_yz = buffer_pfdd[472];

    auto g_y_xyy_xx_zz = buffer_pfdd[473];

    auto g_y_xyy_xy_xx = buffer_pfdd[474];

    auto g_y_xyy_xy_xy = buffer_pfdd[475];

    auto g_y_xyy_xy_xz = buffer_pfdd[476];

    auto g_y_xyy_xy_yy = buffer_pfdd[477];

    auto g_y_xyy_xy_yz = buffer_pfdd[478];

    auto g_y_xyy_xy_zz = buffer_pfdd[479];

    auto g_y_xyy_xz_xx = buffer_pfdd[480];

    auto g_y_xyy_xz_xy = buffer_pfdd[481];

    auto g_y_xyy_xz_xz = buffer_pfdd[482];

    auto g_y_xyy_xz_yy = buffer_pfdd[483];

    auto g_y_xyy_xz_yz = buffer_pfdd[484];

    auto g_y_xyy_xz_zz = buffer_pfdd[485];

    auto g_y_xyy_yy_xx = buffer_pfdd[486];

    auto g_y_xyy_yy_xy = buffer_pfdd[487];

    auto g_y_xyy_yy_xz = buffer_pfdd[488];

    auto g_y_xyy_yy_yy = buffer_pfdd[489];

    auto g_y_xyy_yy_yz = buffer_pfdd[490];

    auto g_y_xyy_yy_zz = buffer_pfdd[491];

    auto g_y_xyy_yz_xx = buffer_pfdd[492];

    auto g_y_xyy_yz_xy = buffer_pfdd[493];

    auto g_y_xyy_yz_xz = buffer_pfdd[494];

    auto g_y_xyy_yz_yy = buffer_pfdd[495];

    auto g_y_xyy_yz_yz = buffer_pfdd[496];

    auto g_y_xyy_yz_zz = buffer_pfdd[497];

    auto g_y_xyy_zz_xx = buffer_pfdd[498];

    auto g_y_xyy_zz_xy = buffer_pfdd[499];

    auto g_y_xyy_zz_xz = buffer_pfdd[500];

    auto g_y_xyy_zz_yy = buffer_pfdd[501];

    auto g_y_xyy_zz_yz = buffer_pfdd[502];

    auto g_y_xyy_zz_zz = buffer_pfdd[503];

    auto g_y_xyz_xx_xx = buffer_pfdd[504];

    auto g_y_xyz_xx_xy = buffer_pfdd[505];

    auto g_y_xyz_xx_xz = buffer_pfdd[506];

    auto g_y_xyz_xx_yy = buffer_pfdd[507];

    auto g_y_xyz_xx_yz = buffer_pfdd[508];

    auto g_y_xyz_xx_zz = buffer_pfdd[509];

    auto g_y_xyz_xy_xx = buffer_pfdd[510];

    auto g_y_xyz_xy_xy = buffer_pfdd[511];

    auto g_y_xyz_xy_xz = buffer_pfdd[512];

    auto g_y_xyz_xy_yy = buffer_pfdd[513];

    auto g_y_xyz_xy_yz = buffer_pfdd[514];

    auto g_y_xyz_xy_zz = buffer_pfdd[515];

    auto g_y_xyz_xz_xx = buffer_pfdd[516];

    auto g_y_xyz_xz_xy = buffer_pfdd[517];

    auto g_y_xyz_xz_xz = buffer_pfdd[518];

    auto g_y_xyz_xz_yy = buffer_pfdd[519];

    auto g_y_xyz_xz_yz = buffer_pfdd[520];

    auto g_y_xyz_xz_zz = buffer_pfdd[521];

    auto g_y_xyz_yy_xx = buffer_pfdd[522];

    auto g_y_xyz_yy_xy = buffer_pfdd[523];

    auto g_y_xyz_yy_xz = buffer_pfdd[524];

    auto g_y_xyz_yy_yy = buffer_pfdd[525];

    auto g_y_xyz_yy_yz = buffer_pfdd[526];

    auto g_y_xyz_yy_zz = buffer_pfdd[527];

    auto g_y_xyz_yz_xx = buffer_pfdd[528];

    auto g_y_xyz_yz_xy = buffer_pfdd[529];

    auto g_y_xyz_yz_xz = buffer_pfdd[530];

    auto g_y_xyz_yz_yy = buffer_pfdd[531];

    auto g_y_xyz_yz_yz = buffer_pfdd[532];

    auto g_y_xyz_yz_zz = buffer_pfdd[533];

    auto g_y_xyz_zz_xx = buffer_pfdd[534];

    auto g_y_xyz_zz_xy = buffer_pfdd[535];

    auto g_y_xyz_zz_xz = buffer_pfdd[536];

    auto g_y_xyz_zz_yy = buffer_pfdd[537];

    auto g_y_xyz_zz_yz = buffer_pfdd[538];

    auto g_y_xyz_zz_zz = buffer_pfdd[539];

    auto g_y_xzz_xx_xx = buffer_pfdd[540];

    auto g_y_xzz_xx_xy = buffer_pfdd[541];

    auto g_y_xzz_xx_xz = buffer_pfdd[542];

    auto g_y_xzz_xx_yy = buffer_pfdd[543];

    auto g_y_xzz_xx_yz = buffer_pfdd[544];

    auto g_y_xzz_xx_zz = buffer_pfdd[545];

    auto g_y_xzz_xy_xx = buffer_pfdd[546];

    auto g_y_xzz_xy_xy = buffer_pfdd[547];

    auto g_y_xzz_xy_xz = buffer_pfdd[548];

    auto g_y_xzz_xy_yy = buffer_pfdd[549];

    auto g_y_xzz_xy_yz = buffer_pfdd[550];

    auto g_y_xzz_xy_zz = buffer_pfdd[551];

    auto g_y_xzz_xz_xx = buffer_pfdd[552];

    auto g_y_xzz_xz_xy = buffer_pfdd[553];

    auto g_y_xzz_xz_xz = buffer_pfdd[554];

    auto g_y_xzz_xz_yy = buffer_pfdd[555];

    auto g_y_xzz_xz_yz = buffer_pfdd[556];

    auto g_y_xzz_xz_zz = buffer_pfdd[557];

    auto g_y_xzz_yy_xx = buffer_pfdd[558];

    auto g_y_xzz_yy_xy = buffer_pfdd[559];

    auto g_y_xzz_yy_xz = buffer_pfdd[560];

    auto g_y_xzz_yy_yy = buffer_pfdd[561];

    auto g_y_xzz_yy_yz = buffer_pfdd[562];

    auto g_y_xzz_yy_zz = buffer_pfdd[563];

    auto g_y_xzz_yz_xx = buffer_pfdd[564];

    auto g_y_xzz_yz_xy = buffer_pfdd[565];

    auto g_y_xzz_yz_xz = buffer_pfdd[566];

    auto g_y_xzz_yz_yy = buffer_pfdd[567];

    auto g_y_xzz_yz_yz = buffer_pfdd[568];

    auto g_y_xzz_yz_zz = buffer_pfdd[569];

    auto g_y_xzz_zz_xx = buffer_pfdd[570];

    auto g_y_xzz_zz_xy = buffer_pfdd[571];

    auto g_y_xzz_zz_xz = buffer_pfdd[572];

    auto g_y_xzz_zz_yy = buffer_pfdd[573];

    auto g_y_xzz_zz_yz = buffer_pfdd[574];

    auto g_y_xzz_zz_zz = buffer_pfdd[575];

    auto g_y_yyy_xx_xx = buffer_pfdd[576];

    auto g_y_yyy_xx_xy = buffer_pfdd[577];

    auto g_y_yyy_xx_xz = buffer_pfdd[578];

    auto g_y_yyy_xx_yy = buffer_pfdd[579];

    auto g_y_yyy_xx_yz = buffer_pfdd[580];

    auto g_y_yyy_xx_zz = buffer_pfdd[581];

    auto g_y_yyy_xy_xx = buffer_pfdd[582];

    auto g_y_yyy_xy_xy = buffer_pfdd[583];

    auto g_y_yyy_xy_xz = buffer_pfdd[584];

    auto g_y_yyy_xy_yy = buffer_pfdd[585];

    auto g_y_yyy_xy_yz = buffer_pfdd[586];

    auto g_y_yyy_xy_zz = buffer_pfdd[587];

    auto g_y_yyy_xz_xx = buffer_pfdd[588];

    auto g_y_yyy_xz_xy = buffer_pfdd[589];

    auto g_y_yyy_xz_xz = buffer_pfdd[590];

    auto g_y_yyy_xz_yy = buffer_pfdd[591];

    auto g_y_yyy_xz_yz = buffer_pfdd[592];

    auto g_y_yyy_xz_zz = buffer_pfdd[593];

    auto g_y_yyy_yy_xx = buffer_pfdd[594];

    auto g_y_yyy_yy_xy = buffer_pfdd[595];

    auto g_y_yyy_yy_xz = buffer_pfdd[596];

    auto g_y_yyy_yy_yy = buffer_pfdd[597];

    auto g_y_yyy_yy_yz = buffer_pfdd[598];

    auto g_y_yyy_yy_zz = buffer_pfdd[599];

    auto g_y_yyy_yz_xx = buffer_pfdd[600];

    auto g_y_yyy_yz_xy = buffer_pfdd[601];

    auto g_y_yyy_yz_xz = buffer_pfdd[602];

    auto g_y_yyy_yz_yy = buffer_pfdd[603];

    auto g_y_yyy_yz_yz = buffer_pfdd[604];

    auto g_y_yyy_yz_zz = buffer_pfdd[605];

    auto g_y_yyy_zz_xx = buffer_pfdd[606];

    auto g_y_yyy_zz_xy = buffer_pfdd[607];

    auto g_y_yyy_zz_xz = buffer_pfdd[608];

    auto g_y_yyy_zz_yy = buffer_pfdd[609];

    auto g_y_yyy_zz_yz = buffer_pfdd[610];

    auto g_y_yyy_zz_zz = buffer_pfdd[611];

    auto g_y_yyz_xx_xx = buffer_pfdd[612];

    auto g_y_yyz_xx_xy = buffer_pfdd[613];

    auto g_y_yyz_xx_xz = buffer_pfdd[614];

    auto g_y_yyz_xx_yy = buffer_pfdd[615];

    auto g_y_yyz_xx_yz = buffer_pfdd[616];

    auto g_y_yyz_xx_zz = buffer_pfdd[617];

    auto g_y_yyz_xy_xx = buffer_pfdd[618];

    auto g_y_yyz_xy_xy = buffer_pfdd[619];

    auto g_y_yyz_xy_xz = buffer_pfdd[620];

    auto g_y_yyz_xy_yy = buffer_pfdd[621];

    auto g_y_yyz_xy_yz = buffer_pfdd[622];

    auto g_y_yyz_xy_zz = buffer_pfdd[623];

    auto g_y_yyz_xz_xx = buffer_pfdd[624];

    auto g_y_yyz_xz_xy = buffer_pfdd[625];

    auto g_y_yyz_xz_xz = buffer_pfdd[626];

    auto g_y_yyz_xz_yy = buffer_pfdd[627];

    auto g_y_yyz_xz_yz = buffer_pfdd[628];

    auto g_y_yyz_xz_zz = buffer_pfdd[629];

    auto g_y_yyz_yy_xx = buffer_pfdd[630];

    auto g_y_yyz_yy_xy = buffer_pfdd[631];

    auto g_y_yyz_yy_xz = buffer_pfdd[632];

    auto g_y_yyz_yy_yy = buffer_pfdd[633];

    auto g_y_yyz_yy_yz = buffer_pfdd[634];

    auto g_y_yyz_yy_zz = buffer_pfdd[635];

    auto g_y_yyz_yz_xx = buffer_pfdd[636];

    auto g_y_yyz_yz_xy = buffer_pfdd[637];

    auto g_y_yyz_yz_xz = buffer_pfdd[638];

    auto g_y_yyz_yz_yy = buffer_pfdd[639];

    auto g_y_yyz_yz_yz = buffer_pfdd[640];

    auto g_y_yyz_yz_zz = buffer_pfdd[641];

    auto g_y_yyz_zz_xx = buffer_pfdd[642];

    auto g_y_yyz_zz_xy = buffer_pfdd[643];

    auto g_y_yyz_zz_xz = buffer_pfdd[644];

    auto g_y_yyz_zz_yy = buffer_pfdd[645];

    auto g_y_yyz_zz_yz = buffer_pfdd[646];

    auto g_y_yyz_zz_zz = buffer_pfdd[647];

    auto g_y_yzz_xx_xx = buffer_pfdd[648];

    auto g_y_yzz_xx_xy = buffer_pfdd[649];

    auto g_y_yzz_xx_xz = buffer_pfdd[650];

    auto g_y_yzz_xx_yy = buffer_pfdd[651];

    auto g_y_yzz_xx_yz = buffer_pfdd[652];

    auto g_y_yzz_xx_zz = buffer_pfdd[653];

    auto g_y_yzz_xy_xx = buffer_pfdd[654];

    auto g_y_yzz_xy_xy = buffer_pfdd[655];

    auto g_y_yzz_xy_xz = buffer_pfdd[656];

    auto g_y_yzz_xy_yy = buffer_pfdd[657];

    auto g_y_yzz_xy_yz = buffer_pfdd[658];

    auto g_y_yzz_xy_zz = buffer_pfdd[659];

    auto g_y_yzz_xz_xx = buffer_pfdd[660];

    auto g_y_yzz_xz_xy = buffer_pfdd[661];

    auto g_y_yzz_xz_xz = buffer_pfdd[662];

    auto g_y_yzz_xz_yy = buffer_pfdd[663];

    auto g_y_yzz_xz_yz = buffer_pfdd[664];

    auto g_y_yzz_xz_zz = buffer_pfdd[665];

    auto g_y_yzz_yy_xx = buffer_pfdd[666];

    auto g_y_yzz_yy_xy = buffer_pfdd[667];

    auto g_y_yzz_yy_xz = buffer_pfdd[668];

    auto g_y_yzz_yy_yy = buffer_pfdd[669];

    auto g_y_yzz_yy_yz = buffer_pfdd[670];

    auto g_y_yzz_yy_zz = buffer_pfdd[671];

    auto g_y_yzz_yz_xx = buffer_pfdd[672];

    auto g_y_yzz_yz_xy = buffer_pfdd[673];

    auto g_y_yzz_yz_xz = buffer_pfdd[674];

    auto g_y_yzz_yz_yy = buffer_pfdd[675];

    auto g_y_yzz_yz_yz = buffer_pfdd[676];

    auto g_y_yzz_yz_zz = buffer_pfdd[677];

    auto g_y_yzz_zz_xx = buffer_pfdd[678];

    auto g_y_yzz_zz_xy = buffer_pfdd[679];

    auto g_y_yzz_zz_xz = buffer_pfdd[680];

    auto g_y_yzz_zz_yy = buffer_pfdd[681];

    auto g_y_yzz_zz_yz = buffer_pfdd[682];

    auto g_y_yzz_zz_zz = buffer_pfdd[683];

    auto g_y_zzz_xx_xx = buffer_pfdd[684];

    auto g_y_zzz_xx_xy = buffer_pfdd[685];

    auto g_y_zzz_xx_xz = buffer_pfdd[686];

    auto g_y_zzz_xx_yy = buffer_pfdd[687];

    auto g_y_zzz_xx_yz = buffer_pfdd[688];

    auto g_y_zzz_xx_zz = buffer_pfdd[689];

    auto g_y_zzz_xy_xx = buffer_pfdd[690];

    auto g_y_zzz_xy_xy = buffer_pfdd[691];

    auto g_y_zzz_xy_xz = buffer_pfdd[692];

    auto g_y_zzz_xy_yy = buffer_pfdd[693];

    auto g_y_zzz_xy_yz = buffer_pfdd[694];

    auto g_y_zzz_xy_zz = buffer_pfdd[695];

    auto g_y_zzz_xz_xx = buffer_pfdd[696];

    auto g_y_zzz_xz_xy = buffer_pfdd[697];

    auto g_y_zzz_xz_xz = buffer_pfdd[698];

    auto g_y_zzz_xz_yy = buffer_pfdd[699];

    auto g_y_zzz_xz_yz = buffer_pfdd[700];

    auto g_y_zzz_xz_zz = buffer_pfdd[701];

    auto g_y_zzz_yy_xx = buffer_pfdd[702];

    auto g_y_zzz_yy_xy = buffer_pfdd[703];

    auto g_y_zzz_yy_xz = buffer_pfdd[704];

    auto g_y_zzz_yy_yy = buffer_pfdd[705];

    auto g_y_zzz_yy_yz = buffer_pfdd[706];

    auto g_y_zzz_yy_zz = buffer_pfdd[707];

    auto g_y_zzz_yz_xx = buffer_pfdd[708];

    auto g_y_zzz_yz_xy = buffer_pfdd[709];

    auto g_y_zzz_yz_xz = buffer_pfdd[710];

    auto g_y_zzz_yz_yy = buffer_pfdd[711];

    auto g_y_zzz_yz_yz = buffer_pfdd[712];

    auto g_y_zzz_yz_zz = buffer_pfdd[713];

    auto g_y_zzz_zz_xx = buffer_pfdd[714];

    auto g_y_zzz_zz_xy = buffer_pfdd[715];

    auto g_y_zzz_zz_xz = buffer_pfdd[716];

    auto g_y_zzz_zz_yy = buffer_pfdd[717];

    auto g_y_zzz_zz_yz = buffer_pfdd[718];

    auto g_y_zzz_zz_zz = buffer_pfdd[719];

    auto g_z_xxx_xx_xx = buffer_pfdd[720];

    auto g_z_xxx_xx_xy = buffer_pfdd[721];

    auto g_z_xxx_xx_xz = buffer_pfdd[722];

    auto g_z_xxx_xx_yy = buffer_pfdd[723];

    auto g_z_xxx_xx_yz = buffer_pfdd[724];

    auto g_z_xxx_xx_zz = buffer_pfdd[725];

    auto g_z_xxx_xy_xx = buffer_pfdd[726];

    auto g_z_xxx_xy_xy = buffer_pfdd[727];

    auto g_z_xxx_xy_xz = buffer_pfdd[728];

    auto g_z_xxx_xy_yy = buffer_pfdd[729];

    auto g_z_xxx_xy_yz = buffer_pfdd[730];

    auto g_z_xxx_xy_zz = buffer_pfdd[731];

    auto g_z_xxx_xz_xx = buffer_pfdd[732];

    auto g_z_xxx_xz_xy = buffer_pfdd[733];

    auto g_z_xxx_xz_xz = buffer_pfdd[734];

    auto g_z_xxx_xz_yy = buffer_pfdd[735];

    auto g_z_xxx_xz_yz = buffer_pfdd[736];

    auto g_z_xxx_xz_zz = buffer_pfdd[737];

    auto g_z_xxx_yy_xx = buffer_pfdd[738];

    auto g_z_xxx_yy_xy = buffer_pfdd[739];

    auto g_z_xxx_yy_xz = buffer_pfdd[740];

    auto g_z_xxx_yy_yy = buffer_pfdd[741];

    auto g_z_xxx_yy_yz = buffer_pfdd[742];

    auto g_z_xxx_yy_zz = buffer_pfdd[743];

    auto g_z_xxx_yz_xx = buffer_pfdd[744];

    auto g_z_xxx_yz_xy = buffer_pfdd[745];

    auto g_z_xxx_yz_xz = buffer_pfdd[746];

    auto g_z_xxx_yz_yy = buffer_pfdd[747];

    auto g_z_xxx_yz_yz = buffer_pfdd[748];

    auto g_z_xxx_yz_zz = buffer_pfdd[749];

    auto g_z_xxx_zz_xx = buffer_pfdd[750];

    auto g_z_xxx_zz_xy = buffer_pfdd[751];

    auto g_z_xxx_zz_xz = buffer_pfdd[752];

    auto g_z_xxx_zz_yy = buffer_pfdd[753];

    auto g_z_xxx_zz_yz = buffer_pfdd[754];

    auto g_z_xxx_zz_zz = buffer_pfdd[755];

    auto g_z_xxy_xx_xx = buffer_pfdd[756];

    auto g_z_xxy_xx_xy = buffer_pfdd[757];

    auto g_z_xxy_xx_xz = buffer_pfdd[758];

    auto g_z_xxy_xx_yy = buffer_pfdd[759];

    auto g_z_xxy_xx_yz = buffer_pfdd[760];

    auto g_z_xxy_xx_zz = buffer_pfdd[761];

    auto g_z_xxy_xy_xx = buffer_pfdd[762];

    auto g_z_xxy_xy_xy = buffer_pfdd[763];

    auto g_z_xxy_xy_xz = buffer_pfdd[764];

    auto g_z_xxy_xy_yy = buffer_pfdd[765];

    auto g_z_xxy_xy_yz = buffer_pfdd[766];

    auto g_z_xxy_xy_zz = buffer_pfdd[767];

    auto g_z_xxy_xz_xx = buffer_pfdd[768];

    auto g_z_xxy_xz_xy = buffer_pfdd[769];

    auto g_z_xxy_xz_xz = buffer_pfdd[770];

    auto g_z_xxy_xz_yy = buffer_pfdd[771];

    auto g_z_xxy_xz_yz = buffer_pfdd[772];

    auto g_z_xxy_xz_zz = buffer_pfdd[773];

    auto g_z_xxy_yy_xx = buffer_pfdd[774];

    auto g_z_xxy_yy_xy = buffer_pfdd[775];

    auto g_z_xxy_yy_xz = buffer_pfdd[776];

    auto g_z_xxy_yy_yy = buffer_pfdd[777];

    auto g_z_xxy_yy_yz = buffer_pfdd[778];

    auto g_z_xxy_yy_zz = buffer_pfdd[779];

    auto g_z_xxy_yz_xx = buffer_pfdd[780];

    auto g_z_xxy_yz_xy = buffer_pfdd[781];

    auto g_z_xxy_yz_xz = buffer_pfdd[782];

    auto g_z_xxy_yz_yy = buffer_pfdd[783];

    auto g_z_xxy_yz_yz = buffer_pfdd[784];

    auto g_z_xxy_yz_zz = buffer_pfdd[785];

    auto g_z_xxy_zz_xx = buffer_pfdd[786];

    auto g_z_xxy_zz_xy = buffer_pfdd[787];

    auto g_z_xxy_zz_xz = buffer_pfdd[788];

    auto g_z_xxy_zz_yy = buffer_pfdd[789];

    auto g_z_xxy_zz_yz = buffer_pfdd[790];

    auto g_z_xxy_zz_zz = buffer_pfdd[791];

    auto g_z_xxz_xx_xx = buffer_pfdd[792];

    auto g_z_xxz_xx_xy = buffer_pfdd[793];

    auto g_z_xxz_xx_xz = buffer_pfdd[794];

    auto g_z_xxz_xx_yy = buffer_pfdd[795];

    auto g_z_xxz_xx_yz = buffer_pfdd[796];

    auto g_z_xxz_xx_zz = buffer_pfdd[797];

    auto g_z_xxz_xy_xx = buffer_pfdd[798];

    auto g_z_xxz_xy_xy = buffer_pfdd[799];

    auto g_z_xxz_xy_xz = buffer_pfdd[800];

    auto g_z_xxz_xy_yy = buffer_pfdd[801];

    auto g_z_xxz_xy_yz = buffer_pfdd[802];

    auto g_z_xxz_xy_zz = buffer_pfdd[803];

    auto g_z_xxz_xz_xx = buffer_pfdd[804];

    auto g_z_xxz_xz_xy = buffer_pfdd[805];

    auto g_z_xxz_xz_xz = buffer_pfdd[806];

    auto g_z_xxz_xz_yy = buffer_pfdd[807];

    auto g_z_xxz_xz_yz = buffer_pfdd[808];

    auto g_z_xxz_xz_zz = buffer_pfdd[809];

    auto g_z_xxz_yy_xx = buffer_pfdd[810];

    auto g_z_xxz_yy_xy = buffer_pfdd[811];

    auto g_z_xxz_yy_xz = buffer_pfdd[812];

    auto g_z_xxz_yy_yy = buffer_pfdd[813];

    auto g_z_xxz_yy_yz = buffer_pfdd[814];

    auto g_z_xxz_yy_zz = buffer_pfdd[815];

    auto g_z_xxz_yz_xx = buffer_pfdd[816];

    auto g_z_xxz_yz_xy = buffer_pfdd[817];

    auto g_z_xxz_yz_xz = buffer_pfdd[818];

    auto g_z_xxz_yz_yy = buffer_pfdd[819];

    auto g_z_xxz_yz_yz = buffer_pfdd[820];

    auto g_z_xxz_yz_zz = buffer_pfdd[821];

    auto g_z_xxz_zz_xx = buffer_pfdd[822];

    auto g_z_xxz_zz_xy = buffer_pfdd[823];

    auto g_z_xxz_zz_xz = buffer_pfdd[824];

    auto g_z_xxz_zz_yy = buffer_pfdd[825];

    auto g_z_xxz_zz_yz = buffer_pfdd[826];

    auto g_z_xxz_zz_zz = buffer_pfdd[827];

    auto g_z_xyy_xx_xx = buffer_pfdd[828];

    auto g_z_xyy_xx_xy = buffer_pfdd[829];

    auto g_z_xyy_xx_xz = buffer_pfdd[830];

    auto g_z_xyy_xx_yy = buffer_pfdd[831];

    auto g_z_xyy_xx_yz = buffer_pfdd[832];

    auto g_z_xyy_xx_zz = buffer_pfdd[833];

    auto g_z_xyy_xy_xx = buffer_pfdd[834];

    auto g_z_xyy_xy_xy = buffer_pfdd[835];

    auto g_z_xyy_xy_xz = buffer_pfdd[836];

    auto g_z_xyy_xy_yy = buffer_pfdd[837];

    auto g_z_xyy_xy_yz = buffer_pfdd[838];

    auto g_z_xyy_xy_zz = buffer_pfdd[839];

    auto g_z_xyy_xz_xx = buffer_pfdd[840];

    auto g_z_xyy_xz_xy = buffer_pfdd[841];

    auto g_z_xyy_xz_xz = buffer_pfdd[842];

    auto g_z_xyy_xz_yy = buffer_pfdd[843];

    auto g_z_xyy_xz_yz = buffer_pfdd[844];

    auto g_z_xyy_xz_zz = buffer_pfdd[845];

    auto g_z_xyy_yy_xx = buffer_pfdd[846];

    auto g_z_xyy_yy_xy = buffer_pfdd[847];

    auto g_z_xyy_yy_xz = buffer_pfdd[848];

    auto g_z_xyy_yy_yy = buffer_pfdd[849];

    auto g_z_xyy_yy_yz = buffer_pfdd[850];

    auto g_z_xyy_yy_zz = buffer_pfdd[851];

    auto g_z_xyy_yz_xx = buffer_pfdd[852];

    auto g_z_xyy_yz_xy = buffer_pfdd[853];

    auto g_z_xyy_yz_xz = buffer_pfdd[854];

    auto g_z_xyy_yz_yy = buffer_pfdd[855];

    auto g_z_xyy_yz_yz = buffer_pfdd[856];

    auto g_z_xyy_yz_zz = buffer_pfdd[857];

    auto g_z_xyy_zz_xx = buffer_pfdd[858];

    auto g_z_xyy_zz_xy = buffer_pfdd[859];

    auto g_z_xyy_zz_xz = buffer_pfdd[860];

    auto g_z_xyy_zz_yy = buffer_pfdd[861];

    auto g_z_xyy_zz_yz = buffer_pfdd[862];

    auto g_z_xyy_zz_zz = buffer_pfdd[863];

    auto g_z_xyz_xx_xx = buffer_pfdd[864];

    auto g_z_xyz_xx_xy = buffer_pfdd[865];

    auto g_z_xyz_xx_xz = buffer_pfdd[866];

    auto g_z_xyz_xx_yy = buffer_pfdd[867];

    auto g_z_xyz_xx_yz = buffer_pfdd[868];

    auto g_z_xyz_xx_zz = buffer_pfdd[869];

    auto g_z_xyz_xy_xx = buffer_pfdd[870];

    auto g_z_xyz_xy_xy = buffer_pfdd[871];

    auto g_z_xyz_xy_xz = buffer_pfdd[872];

    auto g_z_xyz_xy_yy = buffer_pfdd[873];

    auto g_z_xyz_xy_yz = buffer_pfdd[874];

    auto g_z_xyz_xy_zz = buffer_pfdd[875];

    auto g_z_xyz_xz_xx = buffer_pfdd[876];

    auto g_z_xyz_xz_xy = buffer_pfdd[877];

    auto g_z_xyz_xz_xz = buffer_pfdd[878];

    auto g_z_xyz_xz_yy = buffer_pfdd[879];

    auto g_z_xyz_xz_yz = buffer_pfdd[880];

    auto g_z_xyz_xz_zz = buffer_pfdd[881];

    auto g_z_xyz_yy_xx = buffer_pfdd[882];

    auto g_z_xyz_yy_xy = buffer_pfdd[883];

    auto g_z_xyz_yy_xz = buffer_pfdd[884];

    auto g_z_xyz_yy_yy = buffer_pfdd[885];

    auto g_z_xyz_yy_yz = buffer_pfdd[886];

    auto g_z_xyz_yy_zz = buffer_pfdd[887];

    auto g_z_xyz_yz_xx = buffer_pfdd[888];

    auto g_z_xyz_yz_xy = buffer_pfdd[889];

    auto g_z_xyz_yz_xz = buffer_pfdd[890];

    auto g_z_xyz_yz_yy = buffer_pfdd[891];

    auto g_z_xyz_yz_yz = buffer_pfdd[892];

    auto g_z_xyz_yz_zz = buffer_pfdd[893];

    auto g_z_xyz_zz_xx = buffer_pfdd[894];

    auto g_z_xyz_zz_xy = buffer_pfdd[895];

    auto g_z_xyz_zz_xz = buffer_pfdd[896];

    auto g_z_xyz_zz_yy = buffer_pfdd[897];

    auto g_z_xyz_zz_yz = buffer_pfdd[898];

    auto g_z_xyz_zz_zz = buffer_pfdd[899];

    auto g_z_xzz_xx_xx = buffer_pfdd[900];

    auto g_z_xzz_xx_xy = buffer_pfdd[901];

    auto g_z_xzz_xx_xz = buffer_pfdd[902];

    auto g_z_xzz_xx_yy = buffer_pfdd[903];

    auto g_z_xzz_xx_yz = buffer_pfdd[904];

    auto g_z_xzz_xx_zz = buffer_pfdd[905];

    auto g_z_xzz_xy_xx = buffer_pfdd[906];

    auto g_z_xzz_xy_xy = buffer_pfdd[907];

    auto g_z_xzz_xy_xz = buffer_pfdd[908];

    auto g_z_xzz_xy_yy = buffer_pfdd[909];

    auto g_z_xzz_xy_yz = buffer_pfdd[910];

    auto g_z_xzz_xy_zz = buffer_pfdd[911];

    auto g_z_xzz_xz_xx = buffer_pfdd[912];

    auto g_z_xzz_xz_xy = buffer_pfdd[913];

    auto g_z_xzz_xz_xz = buffer_pfdd[914];

    auto g_z_xzz_xz_yy = buffer_pfdd[915];

    auto g_z_xzz_xz_yz = buffer_pfdd[916];

    auto g_z_xzz_xz_zz = buffer_pfdd[917];

    auto g_z_xzz_yy_xx = buffer_pfdd[918];

    auto g_z_xzz_yy_xy = buffer_pfdd[919];

    auto g_z_xzz_yy_xz = buffer_pfdd[920];

    auto g_z_xzz_yy_yy = buffer_pfdd[921];

    auto g_z_xzz_yy_yz = buffer_pfdd[922];

    auto g_z_xzz_yy_zz = buffer_pfdd[923];

    auto g_z_xzz_yz_xx = buffer_pfdd[924];

    auto g_z_xzz_yz_xy = buffer_pfdd[925];

    auto g_z_xzz_yz_xz = buffer_pfdd[926];

    auto g_z_xzz_yz_yy = buffer_pfdd[927];

    auto g_z_xzz_yz_yz = buffer_pfdd[928];

    auto g_z_xzz_yz_zz = buffer_pfdd[929];

    auto g_z_xzz_zz_xx = buffer_pfdd[930];

    auto g_z_xzz_zz_xy = buffer_pfdd[931];

    auto g_z_xzz_zz_xz = buffer_pfdd[932];

    auto g_z_xzz_zz_yy = buffer_pfdd[933];

    auto g_z_xzz_zz_yz = buffer_pfdd[934];

    auto g_z_xzz_zz_zz = buffer_pfdd[935];

    auto g_z_yyy_xx_xx = buffer_pfdd[936];

    auto g_z_yyy_xx_xy = buffer_pfdd[937];

    auto g_z_yyy_xx_xz = buffer_pfdd[938];

    auto g_z_yyy_xx_yy = buffer_pfdd[939];

    auto g_z_yyy_xx_yz = buffer_pfdd[940];

    auto g_z_yyy_xx_zz = buffer_pfdd[941];

    auto g_z_yyy_xy_xx = buffer_pfdd[942];

    auto g_z_yyy_xy_xy = buffer_pfdd[943];

    auto g_z_yyy_xy_xz = buffer_pfdd[944];

    auto g_z_yyy_xy_yy = buffer_pfdd[945];

    auto g_z_yyy_xy_yz = buffer_pfdd[946];

    auto g_z_yyy_xy_zz = buffer_pfdd[947];

    auto g_z_yyy_xz_xx = buffer_pfdd[948];

    auto g_z_yyy_xz_xy = buffer_pfdd[949];

    auto g_z_yyy_xz_xz = buffer_pfdd[950];

    auto g_z_yyy_xz_yy = buffer_pfdd[951];

    auto g_z_yyy_xz_yz = buffer_pfdd[952];

    auto g_z_yyy_xz_zz = buffer_pfdd[953];

    auto g_z_yyy_yy_xx = buffer_pfdd[954];

    auto g_z_yyy_yy_xy = buffer_pfdd[955];

    auto g_z_yyy_yy_xz = buffer_pfdd[956];

    auto g_z_yyy_yy_yy = buffer_pfdd[957];

    auto g_z_yyy_yy_yz = buffer_pfdd[958];

    auto g_z_yyy_yy_zz = buffer_pfdd[959];

    auto g_z_yyy_yz_xx = buffer_pfdd[960];

    auto g_z_yyy_yz_xy = buffer_pfdd[961];

    auto g_z_yyy_yz_xz = buffer_pfdd[962];

    auto g_z_yyy_yz_yy = buffer_pfdd[963];

    auto g_z_yyy_yz_yz = buffer_pfdd[964];

    auto g_z_yyy_yz_zz = buffer_pfdd[965];

    auto g_z_yyy_zz_xx = buffer_pfdd[966];

    auto g_z_yyy_zz_xy = buffer_pfdd[967];

    auto g_z_yyy_zz_xz = buffer_pfdd[968];

    auto g_z_yyy_zz_yy = buffer_pfdd[969];

    auto g_z_yyy_zz_yz = buffer_pfdd[970];

    auto g_z_yyy_zz_zz = buffer_pfdd[971];

    auto g_z_yyz_xx_xx = buffer_pfdd[972];

    auto g_z_yyz_xx_xy = buffer_pfdd[973];

    auto g_z_yyz_xx_xz = buffer_pfdd[974];

    auto g_z_yyz_xx_yy = buffer_pfdd[975];

    auto g_z_yyz_xx_yz = buffer_pfdd[976];

    auto g_z_yyz_xx_zz = buffer_pfdd[977];

    auto g_z_yyz_xy_xx = buffer_pfdd[978];

    auto g_z_yyz_xy_xy = buffer_pfdd[979];

    auto g_z_yyz_xy_xz = buffer_pfdd[980];

    auto g_z_yyz_xy_yy = buffer_pfdd[981];

    auto g_z_yyz_xy_yz = buffer_pfdd[982];

    auto g_z_yyz_xy_zz = buffer_pfdd[983];

    auto g_z_yyz_xz_xx = buffer_pfdd[984];

    auto g_z_yyz_xz_xy = buffer_pfdd[985];

    auto g_z_yyz_xz_xz = buffer_pfdd[986];

    auto g_z_yyz_xz_yy = buffer_pfdd[987];

    auto g_z_yyz_xz_yz = buffer_pfdd[988];

    auto g_z_yyz_xz_zz = buffer_pfdd[989];

    auto g_z_yyz_yy_xx = buffer_pfdd[990];

    auto g_z_yyz_yy_xy = buffer_pfdd[991];

    auto g_z_yyz_yy_xz = buffer_pfdd[992];

    auto g_z_yyz_yy_yy = buffer_pfdd[993];

    auto g_z_yyz_yy_yz = buffer_pfdd[994];

    auto g_z_yyz_yy_zz = buffer_pfdd[995];

    auto g_z_yyz_yz_xx = buffer_pfdd[996];

    auto g_z_yyz_yz_xy = buffer_pfdd[997];

    auto g_z_yyz_yz_xz = buffer_pfdd[998];

    auto g_z_yyz_yz_yy = buffer_pfdd[999];

    auto g_z_yyz_yz_yz = buffer_pfdd[1000];

    auto g_z_yyz_yz_zz = buffer_pfdd[1001];

    auto g_z_yyz_zz_xx = buffer_pfdd[1002];

    auto g_z_yyz_zz_xy = buffer_pfdd[1003];

    auto g_z_yyz_zz_xz = buffer_pfdd[1004];

    auto g_z_yyz_zz_yy = buffer_pfdd[1005];

    auto g_z_yyz_zz_yz = buffer_pfdd[1006];

    auto g_z_yyz_zz_zz = buffer_pfdd[1007];

    auto g_z_yzz_xx_xx = buffer_pfdd[1008];

    auto g_z_yzz_xx_xy = buffer_pfdd[1009];

    auto g_z_yzz_xx_xz = buffer_pfdd[1010];

    auto g_z_yzz_xx_yy = buffer_pfdd[1011];

    auto g_z_yzz_xx_yz = buffer_pfdd[1012];

    auto g_z_yzz_xx_zz = buffer_pfdd[1013];

    auto g_z_yzz_xy_xx = buffer_pfdd[1014];

    auto g_z_yzz_xy_xy = buffer_pfdd[1015];

    auto g_z_yzz_xy_xz = buffer_pfdd[1016];

    auto g_z_yzz_xy_yy = buffer_pfdd[1017];

    auto g_z_yzz_xy_yz = buffer_pfdd[1018];

    auto g_z_yzz_xy_zz = buffer_pfdd[1019];

    auto g_z_yzz_xz_xx = buffer_pfdd[1020];

    auto g_z_yzz_xz_xy = buffer_pfdd[1021];

    auto g_z_yzz_xz_xz = buffer_pfdd[1022];

    auto g_z_yzz_xz_yy = buffer_pfdd[1023];

    auto g_z_yzz_xz_yz = buffer_pfdd[1024];

    auto g_z_yzz_xz_zz = buffer_pfdd[1025];

    auto g_z_yzz_yy_xx = buffer_pfdd[1026];

    auto g_z_yzz_yy_xy = buffer_pfdd[1027];

    auto g_z_yzz_yy_xz = buffer_pfdd[1028];

    auto g_z_yzz_yy_yy = buffer_pfdd[1029];

    auto g_z_yzz_yy_yz = buffer_pfdd[1030];

    auto g_z_yzz_yy_zz = buffer_pfdd[1031];

    auto g_z_yzz_yz_xx = buffer_pfdd[1032];

    auto g_z_yzz_yz_xy = buffer_pfdd[1033];

    auto g_z_yzz_yz_xz = buffer_pfdd[1034];

    auto g_z_yzz_yz_yy = buffer_pfdd[1035];

    auto g_z_yzz_yz_yz = buffer_pfdd[1036];

    auto g_z_yzz_yz_zz = buffer_pfdd[1037];

    auto g_z_yzz_zz_xx = buffer_pfdd[1038];

    auto g_z_yzz_zz_xy = buffer_pfdd[1039];

    auto g_z_yzz_zz_xz = buffer_pfdd[1040];

    auto g_z_yzz_zz_yy = buffer_pfdd[1041];

    auto g_z_yzz_zz_yz = buffer_pfdd[1042];

    auto g_z_yzz_zz_zz = buffer_pfdd[1043];

    auto g_z_zzz_xx_xx = buffer_pfdd[1044];

    auto g_z_zzz_xx_xy = buffer_pfdd[1045];

    auto g_z_zzz_xx_xz = buffer_pfdd[1046];

    auto g_z_zzz_xx_yy = buffer_pfdd[1047];

    auto g_z_zzz_xx_yz = buffer_pfdd[1048];

    auto g_z_zzz_xx_zz = buffer_pfdd[1049];

    auto g_z_zzz_xy_xx = buffer_pfdd[1050];

    auto g_z_zzz_xy_xy = buffer_pfdd[1051];

    auto g_z_zzz_xy_xz = buffer_pfdd[1052];

    auto g_z_zzz_xy_yy = buffer_pfdd[1053];

    auto g_z_zzz_xy_yz = buffer_pfdd[1054];

    auto g_z_zzz_xy_zz = buffer_pfdd[1055];

    auto g_z_zzz_xz_xx = buffer_pfdd[1056];

    auto g_z_zzz_xz_xy = buffer_pfdd[1057];

    auto g_z_zzz_xz_xz = buffer_pfdd[1058];

    auto g_z_zzz_xz_yy = buffer_pfdd[1059];

    auto g_z_zzz_xz_yz = buffer_pfdd[1060];

    auto g_z_zzz_xz_zz = buffer_pfdd[1061];

    auto g_z_zzz_yy_xx = buffer_pfdd[1062];

    auto g_z_zzz_yy_xy = buffer_pfdd[1063];

    auto g_z_zzz_yy_xz = buffer_pfdd[1064];

    auto g_z_zzz_yy_yy = buffer_pfdd[1065];

    auto g_z_zzz_yy_yz = buffer_pfdd[1066];

    auto g_z_zzz_yy_zz = buffer_pfdd[1067];

    auto g_z_zzz_yz_xx = buffer_pfdd[1068];

    auto g_z_zzz_yz_xy = buffer_pfdd[1069];

    auto g_z_zzz_yz_xz = buffer_pfdd[1070];

    auto g_z_zzz_yz_yy = buffer_pfdd[1071];

    auto g_z_zzz_yz_yz = buffer_pfdd[1072];

    auto g_z_zzz_yz_zz = buffer_pfdd[1073];

    auto g_z_zzz_zz_xx = buffer_pfdd[1074];

    auto g_z_zzz_zz_xy = buffer_pfdd[1075];

    auto g_z_zzz_zz_xz = buffer_pfdd[1076];

    auto g_z_zzz_zz_yy = buffer_pfdd[1077];

    auto g_z_zzz_zz_yz = buffer_pfdd[1078];

    auto g_z_zzz_zz_zz = buffer_pfdd[1079];

    /// Set up components of integrals buffer : buffer_1100_sddd

    auto g_x_x_0_0_0_xx_xx_xx = buffer_1100_sddd[0];

    auto g_x_x_0_0_0_xx_xx_xy = buffer_1100_sddd[1];

    auto g_x_x_0_0_0_xx_xx_xz = buffer_1100_sddd[2];

    auto g_x_x_0_0_0_xx_xx_yy = buffer_1100_sddd[3];

    auto g_x_x_0_0_0_xx_xx_yz = buffer_1100_sddd[4];

    auto g_x_x_0_0_0_xx_xx_zz = buffer_1100_sddd[5];

    auto g_x_x_0_0_0_xx_xy_xx = buffer_1100_sddd[6];

    auto g_x_x_0_0_0_xx_xy_xy = buffer_1100_sddd[7];

    auto g_x_x_0_0_0_xx_xy_xz = buffer_1100_sddd[8];

    auto g_x_x_0_0_0_xx_xy_yy = buffer_1100_sddd[9];

    auto g_x_x_0_0_0_xx_xy_yz = buffer_1100_sddd[10];

    auto g_x_x_0_0_0_xx_xy_zz = buffer_1100_sddd[11];

    auto g_x_x_0_0_0_xx_xz_xx = buffer_1100_sddd[12];

    auto g_x_x_0_0_0_xx_xz_xy = buffer_1100_sddd[13];

    auto g_x_x_0_0_0_xx_xz_xz = buffer_1100_sddd[14];

    auto g_x_x_0_0_0_xx_xz_yy = buffer_1100_sddd[15];

    auto g_x_x_0_0_0_xx_xz_yz = buffer_1100_sddd[16];

    auto g_x_x_0_0_0_xx_xz_zz = buffer_1100_sddd[17];

    auto g_x_x_0_0_0_xx_yy_xx = buffer_1100_sddd[18];

    auto g_x_x_0_0_0_xx_yy_xy = buffer_1100_sddd[19];

    auto g_x_x_0_0_0_xx_yy_xz = buffer_1100_sddd[20];

    auto g_x_x_0_0_0_xx_yy_yy = buffer_1100_sddd[21];

    auto g_x_x_0_0_0_xx_yy_yz = buffer_1100_sddd[22];

    auto g_x_x_0_0_0_xx_yy_zz = buffer_1100_sddd[23];

    auto g_x_x_0_0_0_xx_yz_xx = buffer_1100_sddd[24];

    auto g_x_x_0_0_0_xx_yz_xy = buffer_1100_sddd[25];

    auto g_x_x_0_0_0_xx_yz_xz = buffer_1100_sddd[26];

    auto g_x_x_0_0_0_xx_yz_yy = buffer_1100_sddd[27];

    auto g_x_x_0_0_0_xx_yz_yz = buffer_1100_sddd[28];

    auto g_x_x_0_0_0_xx_yz_zz = buffer_1100_sddd[29];

    auto g_x_x_0_0_0_xx_zz_xx = buffer_1100_sddd[30];

    auto g_x_x_0_0_0_xx_zz_xy = buffer_1100_sddd[31];

    auto g_x_x_0_0_0_xx_zz_xz = buffer_1100_sddd[32];

    auto g_x_x_0_0_0_xx_zz_yy = buffer_1100_sddd[33];

    auto g_x_x_0_0_0_xx_zz_yz = buffer_1100_sddd[34];

    auto g_x_x_0_0_0_xx_zz_zz = buffer_1100_sddd[35];

    auto g_x_x_0_0_0_xy_xx_xx = buffer_1100_sddd[36];

    auto g_x_x_0_0_0_xy_xx_xy = buffer_1100_sddd[37];

    auto g_x_x_0_0_0_xy_xx_xz = buffer_1100_sddd[38];

    auto g_x_x_0_0_0_xy_xx_yy = buffer_1100_sddd[39];

    auto g_x_x_0_0_0_xy_xx_yz = buffer_1100_sddd[40];

    auto g_x_x_0_0_0_xy_xx_zz = buffer_1100_sddd[41];

    auto g_x_x_0_0_0_xy_xy_xx = buffer_1100_sddd[42];

    auto g_x_x_0_0_0_xy_xy_xy = buffer_1100_sddd[43];

    auto g_x_x_0_0_0_xy_xy_xz = buffer_1100_sddd[44];

    auto g_x_x_0_0_0_xy_xy_yy = buffer_1100_sddd[45];

    auto g_x_x_0_0_0_xy_xy_yz = buffer_1100_sddd[46];

    auto g_x_x_0_0_0_xy_xy_zz = buffer_1100_sddd[47];

    auto g_x_x_0_0_0_xy_xz_xx = buffer_1100_sddd[48];

    auto g_x_x_0_0_0_xy_xz_xy = buffer_1100_sddd[49];

    auto g_x_x_0_0_0_xy_xz_xz = buffer_1100_sddd[50];

    auto g_x_x_0_0_0_xy_xz_yy = buffer_1100_sddd[51];

    auto g_x_x_0_0_0_xy_xz_yz = buffer_1100_sddd[52];

    auto g_x_x_0_0_0_xy_xz_zz = buffer_1100_sddd[53];

    auto g_x_x_0_0_0_xy_yy_xx = buffer_1100_sddd[54];

    auto g_x_x_0_0_0_xy_yy_xy = buffer_1100_sddd[55];

    auto g_x_x_0_0_0_xy_yy_xz = buffer_1100_sddd[56];

    auto g_x_x_0_0_0_xy_yy_yy = buffer_1100_sddd[57];

    auto g_x_x_0_0_0_xy_yy_yz = buffer_1100_sddd[58];

    auto g_x_x_0_0_0_xy_yy_zz = buffer_1100_sddd[59];

    auto g_x_x_0_0_0_xy_yz_xx = buffer_1100_sddd[60];

    auto g_x_x_0_0_0_xy_yz_xy = buffer_1100_sddd[61];

    auto g_x_x_0_0_0_xy_yz_xz = buffer_1100_sddd[62];

    auto g_x_x_0_0_0_xy_yz_yy = buffer_1100_sddd[63];

    auto g_x_x_0_0_0_xy_yz_yz = buffer_1100_sddd[64];

    auto g_x_x_0_0_0_xy_yz_zz = buffer_1100_sddd[65];

    auto g_x_x_0_0_0_xy_zz_xx = buffer_1100_sddd[66];

    auto g_x_x_0_0_0_xy_zz_xy = buffer_1100_sddd[67];

    auto g_x_x_0_0_0_xy_zz_xz = buffer_1100_sddd[68];

    auto g_x_x_0_0_0_xy_zz_yy = buffer_1100_sddd[69];

    auto g_x_x_0_0_0_xy_zz_yz = buffer_1100_sddd[70];

    auto g_x_x_0_0_0_xy_zz_zz = buffer_1100_sddd[71];

    auto g_x_x_0_0_0_xz_xx_xx = buffer_1100_sddd[72];

    auto g_x_x_0_0_0_xz_xx_xy = buffer_1100_sddd[73];

    auto g_x_x_0_0_0_xz_xx_xz = buffer_1100_sddd[74];

    auto g_x_x_0_0_0_xz_xx_yy = buffer_1100_sddd[75];

    auto g_x_x_0_0_0_xz_xx_yz = buffer_1100_sddd[76];

    auto g_x_x_0_0_0_xz_xx_zz = buffer_1100_sddd[77];

    auto g_x_x_0_0_0_xz_xy_xx = buffer_1100_sddd[78];

    auto g_x_x_0_0_0_xz_xy_xy = buffer_1100_sddd[79];

    auto g_x_x_0_0_0_xz_xy_xz = buffer_1100_sddd[80];

    auto g_x_x_0_0_0_xz_xy_yy = buffer_1100_sddd[81];

    auto g_x_x_0_0_0_xz_xy_yz = buffer_1100_sddd[82];

    auto g_x_x_0_0_0_xz_xy_zz = buffer_1100_sddd[83];

    auto g_x_x_0_0_0_xz_xz_xx = buffer_1100_sddd[84];

    auto g_x_x_0_0_0_xz_xz_xy = buffer_1100_sddd[85];

    auto g_x_x_0_0_0_xz_xz_xz = buffer_1100_sddd[86];

    auto g_x_x_0_0_0_xz_xz_yy = buffer_1100_sddd[87];

    auto g_x_x_0_0_0_xz_xz_yz = buffer_1100_sddd[88];

    auto g_x_x_0_0_0_xz_xz_zz = buffer_1100_sddd[89];

    auto g_x_x_0_0_0_xz_yy_xx = buffer_1100_sddd[90];

    auto g_x_x_0_0_0_xz_yy_xy = buffer_1100_sddd[91];

    auto g_x_x_0_0_0_xz_yy_xz = buffer_1100_sddd[92];

    auto g_x_x_0_0_0_xz_yy_yy = buffer_1100_sddd[93];

    auto g_x_x_0_0_0_xz_yy_yz = buffer_1100_sddd[94];

    auto g_x_x_0_0_0_xz_yy_zz = buffer_1100_sddd[95];

    auto g_x_x_0_0_0_xz_yz_xx = buffer_1100_sddd[96];

    auto g_x_x_0_0_0_xz_yz_xy = buffer_1100_sddd[97];

    auto g_x_x_0_0_0_xz_yz_xz = buffer_1100_sddd[98];

    auto g_x_x_0_0_0_xz_yz_yy = buffer_1100_sddd[99];

    auto g_x_x_0_0_0_xz_yz_yz = buffer_1100_sddd[100];

    auto g_x_x_0_0_0_xz_yz_zz = buffer_1100_sddd[101];

    auto g_x_x_0_0_0_xz_zz_xx = buffer_1100_sddd[102];

    auto g_x_x_0_0_0_xz_zz_xy = buffer_1100_sddd[103];

    auto g_x_x_0_0_0_xz_zz_xz = buffer_1100_sddd[104];

    auto g_x_x_0_0_0_xz_zz_yy = buffer_1100_sddd[105];

    auto g_x_x_0_0_0_xz_zz_yz = buffer_1100_sddd[106];

    auto g_x_x_0_0_0_xz_zz_zz = buffer_1100_sddd[107];

    auto g_x_x_0_0_0_yy_xx_xx = buffer_1100_sddd[108];

    auto g_x_x_0_0_0_yy_xx_xy = buffer_1100_sddd[109];

    auto g_x_x_0_0_0_yy_xx_xz = buffer_1100_sddd[110];

    auto g_x_x_0_0_0_yy_xx_yy = buffer_1100_sddd[111];

    auto g_x_x_0_0_0_yy_xx_yz = buffer_1100_sddd[112];

    auto g_x_x_0_0_0_yy_xx_zz = buffer_1100_sddd[113];

    auto g_x_x_0_0_0_yy_xy_xx = buffer_1100_sddd[114];

    auto g_x_x_0_0_0_yy_xy_xy = buffer_1100_sddd[115];

    auto g_x_x_0_0_0_yy_xy_xz = buffer_1100_sddd[116];

    auto g_x_x_0_0_0_yy_xy_yy = buffer_1100_sddd[117];

    auto g_x_x_0_0_0_yy_xy_yz = buffer_1100_sddd[118];

    auto g_x_x_0_0_0_yy_xy_zz = buffer_1100_sddd[119];

    auto g_x_x_0_0_0_yy_xz_xx = buffer_1100_sddd[120];

    auto g_x_x_0_0_0_yy_xz_xy = buffer_1100_sddd[121];

    auto g_x_x_0_0_0_yy_xz_xz = buffer_1100_sddd[122];

    auto g_x_x_0_0_0_yy_xz_yy = buffer_1100_sddd[123];

    auto g_x_x_0_0_0_yy_xz_yz = buffer_1100_sddd[124];

    auto g_x_x_0_0_0_yy_xz_zz = buffer_1100_sddd[125];

    auto g_x_x_0_0_0_yy_yy_xx = buffer_1100_sddd[126];

    auto g_x_x_0_0_0_yy_yy_xy = buffer_1100_sddd[127];

    auto g_x_x_0_0_0_yy_yy_xz = buffer_1100_sddd[128];

    auto g_x_x_0_0_0_yy_yy_yy = buffer_1100_sddd[129];

    auto g_x_x_0_0_0_yy_yy_yz = buffer_1100_sddd[130];

    auto g_x_x_0_0_0_yy_yy_zz = buffer_1100_sddd[131];

    auto g_x_x_0_0_0_yy_yz_xx = buffer_1100_sddd[132];

    auto g_x_x_0_0_0_yy_yz_xy = buffer_1100_sddd[133];

    auto g_x_x_0_0_0_yy_yz_xz = buffer_1100_sddd[134];

    auto g_x_x_0_0_0_yy_yz_yy = buffer_1100_sddd[135];

    auto g_x_x_0_0_0_yy_yz_yz = buffer_1100_sddd[136];

    auto g_x_x_0_0_0_yy_yz_zz = buffer_1100_sddd[137];

    auto g_x_x_0_0_0_yy_zz_xx = buffer_1100_sddd[138];

    auto g_x_x_0_0_0_yy_zz_xy = buffer_1100_sddd[139];

    auto g_x_x_0_0_0_yy_zz_xz = buffer_1100_sddd[140];

    auto g_x_x_0_0_0_yy_zz_yy = buffer_1100_sddd[141];

    auto g_x_x_0_0_0_yy_zz_yz = buffer_1100_sddd[142];

    auto g_x_x_0_0_0_yy_zz_zz = buffer_1100_sddd[143];

    auto g_x_x_0_0_0_yz_xx_xx = buffer_1100_sddd[144];

    auto g_x_x_0_0_0_yz_xx_xy = buffer_1100_sddd[145];

    auto g_x_x_0_0_0_yz_xx_xz = buffer_1100_sddd[146];

    auto g_x_x_0_0_0_yz_xx_yy = buffer_1100_sddd[147];

    auto g_x_x_0_0_0_yz_xx_yz = buffer_1100_sddd[148];

    auto g_x_x_0_0_0_yz_xx_zz = buffer_1100_sddd[149];

    auto g_x_x_0_0_0_yz_xy_xx = buffer_1100_sddd[150];

    auto g_x_x_0_0_0_yz_xy_xy = buffer_1100_sddd[151];

    auto g_x_x_0_0_0_yz_xy_xz = buffer_1100_sddd[152];

    auto g_x_x_0_0_0_yz_xy_yy = buffer_1100_sddd[153];

    auto g_x_x_0_0_0_yz_xy_yz = buffer_1100_sddd[154];

    auto g_x_x_0_0_0_yz_xy_zz = buffer_1100_sddd[155];

    auto g_x_x_0_0_0_yz_xz_xx = buffer_1100_sddd[156];

    auto g_x_x_0_0_0_yz_xz_xy = buffer_1100_sddd[157];

    auto g_x_x_0_0_0_yz_xz_xz = buffer_1100_sddd[158];

    auto g_x_x_0_0_0_yz_xz_yy = buffer_1100_sddd[159];

    auto g_x_x_0_0_0_yz_xz_yz = buffer_1100_sddd[160];

    auto g_x_x_0_0_0_yz_xz_zz = buffer_1100_sddd[161];

    auto g_x_x_0_0_0_yz_yy_xx = buffer_1100_sddd[162];

    auto g_x_x_0_0_0_yz_yy_xy = buffer_1100_sddd[163];

    auto g_x_x_0_0_0_yz_yy_xz = buffer_1100_sddd[164];

    auto g_x_x_0_0_0_yz_yy_yy = buffer_1100_sddd[165];

    auto g_x_x_0_0_0_yz_yy_yz = buffer_1100_sddd[166];

    auto g_x_x_0_0_0_yz_yy_zz = buffer_1100_sddd[167];

    auto g_x_x_0_0_0_yz_yz_xx = buffer_1100_sddd[168];

    auto g_x_x_0_0_0_yz_yz_xy = buffer_1100_sddd[169];

    auto g_x_x_0_0_0_yz_yz_xz = buffer_1100_sddd[170];

    auto g_x_x_0_0_0_yz_yz_yy = buffer_1100_sddd[171];

    auto g_x_x_0_0_0_yz_yz_yz = buffer_1100_sddd[172];

    auto g_x_x_0_0_0_yz_yz_zz = buffer_1100_sddd[173];

    auto g_x_x_0_0_0_yz_zz_xx = buffer_1100_sddd[174];

    auto g_x_x_0_0_0_yz_zz_xy = buffer_1100_sddd[175];

    auto g_x_x_0_0_0_yz_zz_xz = buffer_1100_sddd[176];

    auto g_x_x_0_0_0_yz_zz_yy = buffer_1100_sddd[177];

    auto g_x_x_0_0_0_yz_zz_yz = buffer_1100_sddd[178];

    auto g_x_x_0_0_0_yz_zz_zz = buffer_1100_sddd[179];

    auto g_x_x_0_0_0_zz_xx_xx = buffer_1100_sddd[180];

    auto g_x_x_0_0_0_zz_xx_xy = buffer_1100_sddd[181];

    auto g_x_x_0_0_0_zz_xx_xz = buffer_1100_sddd[182];

    auto g_x_x_0_0_0_zz_xx_yy = buffer_1100_sddd[183];

    auto g_x_x_0_0_0_zz_xx_yz = buffer_1100_sddd[184];

    auto g_x_x_0_0_0_zz_xx_zz = buffer_1100_sddd[185];

    auto g_x_x_0_0_0_zz_xy_xx = buffer_1100_sddd[186];

    auto g_x_x_0_0_0_zz_xy_xy = buffer_1100_sddd[187];

    auto g_x_x_0_0_0_zz_xy_xz = buffer_1100_sddd[188];

    auto g_x_x_0_0_0_zz_xy_yy = buffer_1100_sddd[189];

    auto g_x_x_0_0_0_zz_xy_yz = buffer_1100_sddd[190];

    auto g_x_x_0_0_0_zz_xy_zz = buffer_1100_sddd[191];

    auto g_x_x_0_0_0_zz_xz_xx = buffer_1100_sddd[192];

    auto g_x_x_0_0_0_zz_xz_xy = buffer_1100_sddd[193];

    auto g_x_x_0_0_0_zz_xz_xz = buffer_1100_sddd[194];

    auto g_x_x_0_0_0_zz_xz_yy = buffer_1100_sddd[195];

    auto g_x_x_0_0_0_zz_xz_yz = buffer_1100_sddd[196];

    auto g_x_x_0_0_0_zz_xz_zz = buffer_1100_sddd[197];

    auto g_x_x_0_0_0_zz_yy_xx = buffer_1100_sddd[198];

    auto g_x_x_0_0_0_zz_yy_xy = buffer_1100_sddd[199];

    auto g_x_x_0_0_0_zz_yy_xz = buffer_1100_sddd[200];

    auto g_x_x_0_0_0_zz_yy_yy = buffer_1100_sddd[201];

    auto g_x_x_0_0_0_zz_yy_yz = buffer_1100_sddd[202];

    auto g_x_x_0_0_0_zz_yy_zz = buffer_1100_sddd[203];

    auto g_x_x_0_0_0_zz_yz_xx = buffer_1100_sddd[204];

    auto g_x_x_0_0_0_zz_yz_xy = buffer_1100_sddd[205];

    auto g_x_x_0_0_0_zz_yz_xz = buffer_1100_sddd[206];

    auto g_x_x_0_0_0_zz_yz_yy = buffer_1100_sddd[207];

    auto g_x_x_0_0_0_zz_yz_yz = buffer_1100_sddd[208];

    auto g_x_x_0_0_0_zz_yz_zz = buffer_1100_sddd[209];

    auto g_x_x_0_0_0_zz_zz_xx = buffer_1100_sddd[210];

    auto g_x_x_0_0_0_zz_zz_xy = buffer_1100_sddd[211];

    auto g_x_x_0_0_0_zz_zz_xz = buffer_1100_sddd[212];

    auto g_x_x_0_0_0_zz_zz_yy = buffer_1100_sddd[213];

    auto g_x_x_0_0_0_zz_zz_yz = buffer_1100_sddd[214];

    auto g_x_x_0_0_0_zz_zz_zz = buffer_1100_sddd[215];

    auto g_x_y_0_0_0_xx_xx_xx = buffer_1100_sddd[216];

    auto g_x_y_0_0_0_xx_xx_xy = buffer_1100_sddd[217];

    auto g_x_y_0_0_0_xx_xx_xz = buffer_1100_sddd[218];

    auto g_x_y_0_0_0_xx_xx_yy = buffer_1100_sddd[219];

    auto g_x_y_0_0_0_xx_xx_yz = buffer_1100_sddd[220];

    auto g_x_y_0_0_0_xx_xx_zz = buffer_1100_sddd[221];

    auto g_x_y_0_0_0_xx_xy_xx = buffer_1100_sddd[222];

    auto g_x_y_0_0_0_xx_xy_xy = buffer_1100_sddd[223];

    auto g_x_y_0_0_0_xx_xy_xz = buffer_1100_sddd[224];

    auto g_x_y_0_0_0_xx_xy_yy = buffer_1100_sddd[225];

    auto g_x_y_0_0_0_xx_xy_yz = buffer_1100_sddd[226];

    auto g_x_y_0_0_0_xx_xy_zz = buffer_1100_sddd[227];

    auto g_x_y_0_0_0_xx_xz_xx = buffer_1100_sddd[228];

    auto g_x_y_0_0_0_xx_xz_xy = buffer_1100_sddd[229];

    auto g_x_y_0_0_0_xx_xz_xz = buffer_1100_sddd[230];

    auto g_x_y_0_0_0_xx_xz_yy = buffer_1100_sddd[231];

    auto g_x_y_0_0_0_xx_xz_yz = buffer_1100_sddd[232];

    auto g_x_y_0_0_0_xx_xz_zz = buffer_1100_sddd[233];

    auto g_x_y_0_0_0_xx_yy_xx = buffer_1100_sddd[234];

    auto g_x_y_0_0_0_xx_yy_xy = buffer_1100_sddd[235];

    auto g_x_y_0_0_0_xx_yy_xz = buffer_1100_sddd[236];

    auto g_x_y_0_0_0_xx_yy_yy = buffer_1100_sddd[237];

    auto g_x_y_0_0_0_xx_yy_yz = buffer_1100_sddd[238];

    auto g_x_y_0_0_0_xx_yy_zz = buffer_1100_sddd[239];

    auto g_x_y_0_0_0_xx_yz_xx = buffer_1100_sddd[240];

    auto g_x_y_0_0_0_xx_yz_xy = buffer_1100_sddd[241];

    auto g_x_y_0_0_0_xx_yz_xz = buffer_1100_sddd[242];

    auto g_x_y_0_0_0_xx_yz_yy = buffer_1100_sddd[243];

    auto g_x_y_0_0_0_xx_yz_yz = buffer_1100_sddd[244];

    auto g_x_y_0_0_0_xx_yz_zz = buffer_1100_sddd[245];

    auto g_x_y_0_0_0_xx_zz_xx = buffer_1100_sddd[246];

    auto g_x_y_0_0_0_xx_zz_xy = buffer_1100_sddd[247];

    auto g_x_y_0_0_0_xx_zz_xz = buffer_1100_sddd[248];

    auto g_x_y_0_0_0_xx_zz_yy = buffer_1100_sddd[249];

    auto g_x_y_0_0_0_xx_zz_yz = buffer_1100_sddd[250];

    auto g_x_y_0_0_0_xx_zz_zz = buffer_1100_sddd[251];

    auto g_x_y_0_0_0_xy_xx_xx = buffer_1100_sddd[252];

    auto g_x_y_0_0_0_xy_xx_xy = buffer_1100_sddd[253];

    auto g_x_y_0_0_0_xy_xx_xz = buffer_1100_sddd[254];

    auto g_x_y_0_0_0_xy_xx_yy = buffer_1100_sddd[255];

    auto g_x_y_0_0_0_xy_xx_yz = buffer_1100_sddd[256];

    auto g_x_y_0_0_0_xy_xx_zz = buffer_1100_sddd[257];

    auto g_x_y_0_0_0_xy_xy_xx = buffer_1100_sddd[258];

    auto g_x_y_0_0_0_xy_xy_xy = buffer_1100_sddd[259];

    auto g_x_y_0_0_0_xy_xy_xz = buffer_1100_sddd[260];

    auto g_x_y_0_0_0_xy_xy_yy = buffer_1100_sddd[261];

    auto g_x_y_0_0_0_xy_xy_yz = buffer_1100_sddd[262];

    auto g_x_y_0_0_0_xy_xy_zz = buffer_1100_sddd[263];

    auto g_x_y_0_0_0_xy_xz_xx = buffer_1100_sddd[264];

    auto g_x_y_0_0_0_xy_xz_xy = buffer_1100_sddd[265];

    auto g_x_y_0_0_0_xy_xz_xz = buffer_1100_sddd[266];

    auto g_x_y_0_0_0_xy_xz_yy = buffer_1100_sddd[267];

    auto g_x_y_0_0_0_xy_xz_yz = buffer_1100_sddd[268];

    auto g_x_y_0_0_0_xy_xz_zz = buffer_1100_sddd[269];

    auto g_x_y_0_0_0_xy_yy_xx = buffer_1100_sddd[270];

    auto g_x_y_0_0_0_xy_yy_xy = buffer_1100_sddd[271];

    auto g_x_y_0_0_0_xy_yy_xz = buffer_1100_sddd[272];

    auto g_x_y_0_0_0_xy_yy_yy = buffer_1100_sddd[273];

    auto g_x_y_0_0_0_xy_yy_yz = buffer_1100_sddd[274];

    auto g_x_y_0_0_0_xy_yy_zz = buffer_1100_sddd[275];

    auto g_x_y_0_0_0_xy_yz_xx = buffer_1100_sddd[276];

    auto g_x_y_0_0_0_xy_yz_xy = buffer_1100_sddd[277];

    auto g_x_y_0_0_0_xy_yz_xz = buffer_1100_sddd[278];

    auto g_x_y_0_0_0_xy_yz_yy = buffer_1100_sddd[279];

    auto g_x_y_0_0_0_xy_yz_yz = buffer_1100_sddd[280];

    auto g_x_y_0_0_0_xy_yz_zz = buffer_1100_sddd[281];

    auto g_x_y_0_0_0_xy_zz_xx = buffer_1100_sddd[282];

    auto g_x_y_0_0_0_xy_zz_xy = buffer_1100_sddd[283];

    auto g_x_y_0_0_0_xy_zz_xz = buffer_1100_sddd[284];

    auto g_x_y_0_0_0_xy_zz_yy = buffer_1100_sddd[285];

    auto g_x_y_0_0_0_xy_zz_yz = buffer_1100_sddd[286];

    auto g_x_y_0_0_0_xy_zz_zz = buffer_1100_sddd[287];

    auto g_x_y_0_0_0_xz_xx_xx = buffer_1100_sddd[288];

    auto g_x_y_0_0_0_xz_xx_xy = buffer_1100_sddd[289];

    auto g_x_y_0_0_0_xz_xx_xz = buffer_1100_sddd[290];

    auto g_x_y_0_0_0_xz_xx_yy = buffer_1100_sddd[291];

    auto g_x_y_0_0_0_xz_xx_yz = buffer_1100_sddd[292];

    auto g_x_y_0_0_0_xz_xx_zz = buffer_1100_sddd[293];

    auto g_x_y_0_0_0_xz_xy_xx = buffer_1100_sddd[294];

    auto g_x_y_0_0_0_xz_xy_xy = buffer_1100_sddd[295];

    auto g_x_y_0_0_0_xz_xy_xz = buffer_1100_sddd[296];

    auto g_x_y_0_0_0_xz_xy_yy = buffer_1100_sddd[297];

    auto g_x_y_0_0_0_xz_xy_yz = buffer_1100_sddd[298];

    auto g_x_y_0_0_0_xz_xy_zz = buffer_1100_sddd[299];

    auto g_x_y_0_0_0_xz_xz_xx = buffer_1100_sddd[300];

    auto g_x_y_0_0_0_xz_xz_xy = buffer_1100_sddd[301];

    auto g_x_y_0_0_0_xz_xz_xz = buffer_1100_sddd[302];

    auto g_x_y_0_0_0_xz_xz_yy = buffer_1100_sddd[303];

    auto g_x_y_0_0_0_xz_xz_yz = buffer_1100_sddd[304];

    auto g_x_y_0_0_0_xz_xz_zz = buffer_1100_sddd[305];

    auto g_x_y_0_0_0_xz_yy_xx = buffer_1100_sddd[306];

    auto g_x_y_0_0_0_xz_yy_xy = buffer_1100_sddd[307];

    auto g_x_y_0_0_0_xz_yy_xz = buffer_1100_sddd[308];

    auto g_x_y_0_0_0_xz_yy_yy = buffer_1100_sddd[309];

    auto g_x_y_0_0_0_xz_yy_yz = buffer_1100_sddd[310];

    auto g_x_y_0_0_0_xz_yy_zz = buffer_1100_sddd[311];

    auto g_x_y_0_0_0_xz_yz_xx = buffer_1100_sddd[312];

    auto g_x_y_0_0_0_xz_yz_xy = buffer_1100_sddd[313];

    auto g_x_y_0_0_0_xz_yz_xz = buffer_1100_sddd[314];

    auto g_x_y_0_0_0_xz_yz_yy = buffer_1100_sddd[315];

    auto g_x_y_0_0_0_xz_yz_yz = buffer_1100_sddd[316];

    auto g_x_y_0_0_0_xz_yz_zz = buffer_1100_sddd[317];

    auto g_x_y_0_0_0_xz_zz_xx = buffer_1100_sddd[318];

    auto g_x_y_0_0_0_xz_zz_xy = buffer_1100_sddd[319];

    auto g_x_y_0_0_0_xz_zz_xz = buffer_1100_sddd[320];

    auto g_x_y_0_0_0_xz_zz_yy = buffer_1100_sddd[321];

    auto g_x_y_0_0_0_xz_zz_yz = buffer_1100_sddd[322];

    auto g_x_y_0_0_0_xz_zz_zz = buffer_1100_sddd[323];

    auto g_x_y_0_0_0_yy_xx_xx = buffer_1100_sddd[324];

    auto g_x_y_0_0_0_yy_xx_xy = buffer_1100_sddd[325];

    auto g_x_y_0_0_0_yy_xx_xz = buffer_1100_sddd[326];

    auto g_x_y_0_0_0_yy_xx_yy = buffer_1100_sddd[327];

    auto g_x_y_0_0_0_yy_xx_yz = buffer_1100_sddd[328];

    auto g_x_y_0_0_0_yy_xx_zz = buffer_1100_sddd[329];

    auto g_x_y_0_0_0_yy_xy_xx = buffer_1100_sddd[330];

    auto g_x_y_0_0_0_yy_xy_xy = buffer_1100_sddd[331];

    auto g_x_y_0_0_0_yy_xy_xz = buffer_1100_sddd[332];

    auto g_x_y_0_0_0_yy_xy_yy = buffer_1100_sddd[333];

    auto g_x_y_0_0_0_yy_xy_yz = buffer_1100_sddd[334];

    auto g_x_y_0_0_0_yy_xy_zz = buffer_1100_sddd[335];

    auto g_x_y_0_0_0_yy_xz_xx = buffer_1100_sddd[336];

    auto g_x_y_0_0_0_yy_xz_xy = buffer_1100_sddd[337];

    auto g_x_y_0_0_0_yy_xz_xz = buffer_1100_sddd[338];

    auto g_x_y_0_0_0_yy_xz_yy = buffer_1100_sddd[339];

    auto g_x_y_0_0_0_yy_xz_yz = buffer_1100_sddd[340];

    auto g_x_y_0_0_0_yy_xz_zz = buffer_1100_sddd[341];

    auto g_x_y_0_0_0_yy_yy_xx = buffer_1100_sddd[342];

    auto g_x_y_0_0_0_yy_yy_xy = buffer_1100_sddd[343];

    auto g_x_y_0_0_0_yy_yy_xz = buffer_1100_sddd[344];

    auto g_x_y_0_0_0_yy_yy_yy = buffer_1100_sddd[345];

    auto g_x_y_0_0_0_yy_yy_yz = buffer_1100_sddd[346];

    auto g_x_y_0_0_0_yy_yy_zz = buffer_1100_sddd[347];

    auto g_x_y_0_0_0_yy_yz_xx = buffer_1100_sddd[348];

    auto g_x_y_0_0_0_yy_yz_xy = buffer_1100_sddd[349];

    auto g_x_y_0_0_0_yy_yz_xz = buffer_1100_sddd[350];

    auto g_x_y_0_0_0_yy_yz_yy = buffer_1100_sddd[351];

    auto g_x_y_0_0_0_yy_yz_yz = buffer_1100_sddd[352];

    auto g_x_y_0_0_0_yy_yz_zz = buffer_1100_sddd[353];

    auto g_x_y_0_0_0_yy_zz_xx = buffer_1100_sddd[354];

    auto g_x_y_0_0_0_yy_zz_xy = buffer_1100_sddd[355];

    auto g_x_y_0_0_0_yy_zz_xz = buffer_1100_sddd[356];

    auto g_x_y_0_0_0_yy_zz_yy = buffer_1100_sddd[357];

    auto g_x_y_0_0_0_yy_zz_yz = buffer_1100_sddd[358];

    auto g_x_y_0_0_0_yy_zz_zz = buffer_1100_sddd[359];

    auto g_x_y_0_0_0_yz_xx_xx = buffer_1100_sddd[360];

    auto g_x_y_0_0_0_yz_xx_xy = buffer_1100_sddd[361];

    auto g_x_y_0_0_0_yz_xx_xz = buffer_1100_sddd[362];

    auto g_x_y_0_0_0_yz_xx_yy = buffer_1100_sddd[363];

    auto g_x_y_0_0_0_yz_xx_yz = buffer_1100_sddd[364];

    auto g_x_y_0_0_0_yz_xx_zz = buffer_1100_sddd[365];

    auto g_x_y_0_0_0_yz_xy_xx = buffer_1100_sddd[366];

    auto g_x_y_0_0_0_yz_xy_xy = buffer_1100_sddd[367];

    auto g_x_y_0_0_0_yz_xy_xz = buffer_1100_sddd[368];

    auto g_x_y_0_0_0_yz_xy_yy = buffer_1100_sddd[369];

    auto g_x_y_0_0_0_yz_xy_yz = buffer_1100_sddd[370];

    auto g_x_y_0_0_0_yz_xy_zz = buffer_1100_sddd[371];

    auto g_x_y_0_0_0_yz_xz_xx = buffer_1100_sddd[372];

    auto g_x_y_0_0_0_yz_xz_xy = buffer_1100_sddd[373];

    auto g_x_y_0_0_0_yz_xz_xz = buffer_1100_sddd[374];

    auto g_x_y_0_0_0_yz_xz_yy = buffer_1100_sddd[375];

    auto g_x_y_0_0_0_yz_xz_yz = buffer_1100_sddd[376];

    auto g_x_y_0_0_0_yz_xz_zz = buffer_1100_sddd[377];

    auto g_x_y_0_0_0_yz_yy_xx = buffer_1100_sddd[378];

    auto g_x_y_0_0_0_yz_yy_xy = buffer_1100_sddd[379];

    auto g_x_y_0_0_0_yz_yy_xz = buffer_1100_sddd[380];

    auto g_x_y_0_0_0_yz_yy_yy = buffer_1100_sddd[381];

    auto g_x_y_0_0_0_yz_yy_yz = buffer_1100_sddd[382];

    auto g_x_y_0_0_0_yz_yy_zz = buffer_1100_sddd[383];

    auto g_x_y_0_0_0_yz_yz_xx = buffer_1100_sddd[384];

    auto g_x_y_0_0_0_yz_yz_xy = buffer_1100_sddd[385];

    auto g_x_y_0_0_0_yz_yz_xz = buffer_1100_sddd[386];

    auto g_x_y_0_0_0_yz_yz_yy = buffer_1100_sddd[387];

    auto g_x_y_0_0_0_yz_yz_yz = buffer_1100_sddd[388];

    auto g_x_y_0_0_0_yz_yz_zz = buffer_1100_sddd[389];

    auto g_x_y_0_0_0_yz_zz_xx = buffer_1100_sddd[390];

    auto g_x_y_0_0_0_yz_zz_xy = buffer_1100_sddd[391];

    auto g_x_y_0_0_0_yz_zz_xz = buffer_1100_sddd[392];

    auto g_x_y_0_0_0_yz_zz_yy = buffer_1100_sddd[393];

    auto g_x_y_0_0_0_yz_zz_yz = buffer_1100_sddd[394];

    auto g_x_y_0_0_0_yz_zz_zz = buffer_1100_sddd[395];

    auto g_x_y_0_0_0_zz_xx_xx = buffer_1100_sddd[396];

    auto g_x_y_0_0_0_zz_xx_xy = buffer_1100_sddd[397];

    auto g_x_y_0_0_0_zz_xx_xz = buffer_1100_sddd[398];

    auto g_x_y_0_0_0_zz_xx_yy = buffer_1100_sddd[399];

    auto g_x_y_0_0_0_zz_xx_yz = buffer_1100_sddd[400];

    auto g_x_y_0_0_0_zz_xx_zz = buffer_1100_sddd[401];

    auto g_x_y_0_0_0_zz_xy_xx = buffer_1100_sddd[402];

    auto g_x_y_0_0_0_zz_xy_xy = buffer_1100_sddd[403];

    auto g_x_y_0_0_0_zz_xy_xz = buffer_1100_sddd[404];

    auto g_x_y_0_0_0_zz_xy_yy = buffer_1100_sddd[405];

    auto g_x_y_0_0_0_zz_xy_yz = buffer_1100_sddd[406];

    auto g_x_y_0_0_0_zz_xy_zz = buffer_1100_sddd[407];

    auto g_x_y_0_0_0_zz_xz_xx = buffer_1100_sddd[408];

    auto g_x_y_0_0_0_zz_xz_xy = buffer_1100_sddd[409];

    auto g_x_y_0_0_0_zz_xz_xz = buffer_1100_sddd[410];

    auto g_x_y_0_0_0_zz_xz_yy = buffer_1100_sddd[411];

    auto g_x_y_0_0_0_zz_xz_yz = buffer_1100_sddd[412];

    auto g_x_y_0_0_0_zz_xz_zz = buffer_1100_sddd[413];

    auto g_x_y_0_0_0_zz_yy_xx = buffer_1100_sddd[414];

    auto g_x_y_0_0_0_zz_yy_xy = buffer_1100_sddd[415];

    auto g_x_y_0_0_0_zz_yy_xz = buffer_1100_sddd[416];

    auto g_x_y_0_0_0_zz_yy_yy = buffer_1100_sddd[417];

    auto g_x_y_0_0_0_zz_yy_yz = buffer_1100_sddd[418];

    auto g_x_y_0_0_0_zz_yy_zz = buffer_1100_sddd[419];

    auto g_x_y_0_0_0_zz_yz_xx = buffer_1100_sddd[420];

    auto g_x_y_0_0_0_zz_yz_xy = buffer_1100_sddd[421];

    auto g_x_y_0_0_0_zz_yz_xz = buffer_1100_sddd[422];

    auto g_x_y_0_0_0_zz_yz_yy = buffer_1100_sddd[423];

    auto g_x_y_0_0_0_zz_yz_yz = buffer_1100_sddd[424];

    auto g_x_y_0_0_0_zz_yz_zz = buffer_1100_sddd[425];

    auto g_x_y_0_0_0_zz_zz_xx = buffer_1100_sddd[426];

    auto g_x_y_0_0_0_zz_zz_xy = buffer_1100_sddd[427];

    auto g_x_y_0_0_0_zz_zz_xz = buffer_1100_sddd[428];

    auto g_x_y_0_0_0_zz_zz_yy = buffer_1100_sddd[429];

    auto g_x_y_0_0_0_zz_zz_yz = buffer_1100_sddd[430];

    auto g_x_y_0_0_0_zz_zz_zz = buffer_1100_sddd[431];

    auto g_x_z_0_0_0_xx_xx_xx = buffer_1100_sddd[432];

    auto g_x_z_0_0_0_xx_xx_xy = buffer_1100_sddd[433];

    auto g_x_z_0_0_0_xx_xx_xz = buffer_1100_sddd[434];

    auto g_x_z_0_0_0_xx_xx_yy = buffer_1100_sddd[435];

    auto g_x_z_0_0_0_xx_xx_yz = buffer_1100_sddd[436];

    auto g_x_z_0_0_0_xx_xx_zz = buffer_1100_sddd[437];

    auto g_x_z_0_0_0_xx_xy_xx = buffer_1100_sddd[438];

    auto g_x_z_0_0_0_xx_xy_xy = buffer_1100_sddd[439];

    auto g_x_z_0_0_0_xx_xy_xz = buffer_1100_sddd[440];

    auto g_x_z_0_0_0_xx_xy_yy = buffer_1100_sddd[441];

    auto g_x_z_0_0_0_xx_xy_yz = buffer_1100_sddd[442];

    auto g_x_z_0_0_0_xx_xy_zz = buffer_1100_sddd[443];

    auto g_x_z_0_0_0_xx_xz_xx = buffer_1100_sddd[444];

    auto g_x_z_0_0_0_xx_xz_xy = buffer_1100_sddd[445];

    auto g_x_z_0_0_0_xx_xz_xz = buffer_1100_sddd[446];

    auto g_x_z_0_0_0_xx_xz_yy = buffer_1100_sddd[447];

    auto g_x_z_0_0_0_xx_xz_yz = buffer_1100_sddd[448];

    auto g_x_z_0_0_0_xx_xz_zz = buffer_1100_sddd[449];

    auto g_x_z_0_0_0_xx_yy_xx = buffer_1100_sddd[450];

    auto g_x_z_0_0_0_xx_yy_xy = buffer_1100_sddd[451];

    auto g_x_z_0_0_0_xx_yy_xz = buffer_1100_sddd[452];

    auto g_x_z_0_0_0_xx_yy_yy = buffer_1100_sddd[453];

    auto g_x_z_0_0_0_xx_yy_yz = buffer_1100_sddd[454];

    auto g_x_z_0_0_0_xx_yy_zz = buffer_1100_sddd[455];

    auto g_x_z_0_0_0_xx_yz_xx = buffer_1100_sddd[456];

    auto g_x_z_0_0_0_xx_yz_xy = buffer_1100_sddd[457];

    auto g_x_z_0_0_0_xx_yz_xz = buffer_1100_sddd[458];

    auto g_x_z_0_0_0_xx_yz_yy = buffer_1100_sddd[459];

    auto g_x_z_0_0_0_xx_yz_yz = buffer_1100_sddd[460];

    auto g_x_z_0_0_0_xx_yz_zz = buffer_1100_sddd[461];

    auto g_x_z_0_0_0_xx_zz_xx = buffer_1100_sddd[462];

    auto g_x_z_0_0_0_xx_zz_xy = buffer_1100_sddd[463];

    auto g_x_z_0_0_0_xx_zz_xz = buffer_1100_sddd[464];

    auto g_x_z_0_0_0_xx_zz_yy = buffer_1100_sddd[465];

    auto g_x_z_0_0_0_xx_zz_yz = buffer_1100_sddd[466];

    auto g_x_z_0_0_0_xx_zz_zz = buffer_1100_sddd[467];

    auto g_x_z_0_0_0_xy_xx_xx = buffer_1100_sddd[468];

    auto g_x_z_0_0_0_xy_xx_xy = buffer_1100_sddd[469];

    auto g_x_z_0_0_0_xy_xx_xz = buffer_1100_sddd[470];

    auto g_x_z_0_0_0_xy_xx_yy = buffer_1100_sddd[471];

    auto g_x_z_0_0_0_xy_xx_yz = buffer_1100_sddd[472];

    auto g_x_z_0_0_0_xy_xx_zz = buffer_1100_sddd[473];

    auto g_x_z_0_0_0_xy_xy_xx = buffer_1100_sddd[474];

    auto g_x_z_0_0_0_xy_xy_xy = buffer_1100_sddd[475];

    auto g_x_z_0_0_0_xy_xy_xz = buffer_1100_sddd[476];

    auto g_x_z_0_0_0_xy_xy_yy = buffer_1100_sddd[477];

    auto g_x_z_0_0_0_xy_xy_yz = buffer_1100_sddd[478];

    auto g_x_z_0_0_0_xy_xy_zz = buffer_1100_sddd[479];

    auto g_x_z_0_0_0_xy_xz_xx = buffer_1100_sddd[480];

    auto g_x_z_0_0_0_xy_xz_xy = buffer_1100_sddd[481];

    auto g_x_z_0_0_0_xy_xz_xz = buffer_1100_sddd[482];

    auto g_x_z_0_0_0_xy_xz_yy = buffer_1100_sddd[483];

    auto g_x_z_0_0_0_xy_xz_yz = buffer_1100_sddd[484];

    auto g_x_z_0_0_0_xy_xz_zz = buffer_1100_sddd[485];

    auto g_x_z_0_0_0_xy_yy_xx = buffer_1100_sddd[486];

    auto g_x_z_0_0_0_xy_yy_xy = buffer_1100_sddd[487];

    auto g_x_z_0_0_0_xy_yy_xz = buffer_1100_sddd[488];

    auto g_x_z_0_0_0_xy_yy_yy = buffer_1100_sddd[489];

    auto g_x_z_0_0_0_xy_yy_yz = buffer_1100_sddd[490];

    auto g_x_z_0_0_0_xy_yy_zz = buffer_1100_sddd[491];

    auto g_x_z_0_0_0_xy_yz_xx = buffer_1100_sddd[492];

    auto g_x_z_0_0_0_xy_yz_xy = buffer_1100_sddd[493];

    auto g_x_z_0_0_0_xy_yz_xz = buffer_1100_sddd[494];

    auto g_x_z_0_0_0_xy_yz_yy = buffer_1100_sddd[495];

    auto g_x_z_0_0_0_xy_yz_yz = buffer_1100_sddd[496];

    auto g_x_z_0_0_0_xy_yz_zz = buffer_1100_sddd[497];

    auto g_x_z_0_0_0_xy_zz_xx = buffer_1100_sddd[498];

    auto g_x_z_0_0_0_xy_zz_xy = buffer_1100_sddd[499];

    auto g_x_z_0_0_0_xy_zz_xz = buffer_1100_sddd[500];

    auto g_x_z_0_0_0_xy_zz_yy = buffer_1100_sddd[501];

    auto g_x_z_0_0_0_xy_zz_yz = buffer_1100_sddd[502];

    auto g_x_z_0_0_0_xy_zz_zz = buffer_1100_sddd[503];

    auto g_x_z_0_0_0_xz_xx_xx = buffer_1100_sddd[504];

    auto g_x_z_0_0_0_xz_xx_xy = buffer_1100_sddd[505];

    auto g_x_z_0_0_0_xz_xx_xz = buffer_1100_sddd[506];

    auto g_x_z_0_0_0_xz_xx_yy = buffer_1100_sddd[507];

    auto g_x_z_0_0_0_xz_xx_yz = buffer_1100_sddd[508];

    auto g_x_z_0_0_0_xz_xx_zz = buffer_1100_sddd[509];

    auto g_x_z_0_0_0_xz_xy_xx = buffer_1100_sddd[510];

    auto g_x_z_0_0_0_xz_xy_xy = buffer_1100_sddd[511];

    auto g_x_z_0_0_0_xz_xy_xz = buffer_1100_sddd[512];

    auto g_x_z_0_0_0_xz_xy_yy = buffer_1100_sddd[513];

    auto g_x_z_0_0_0_xz_xy_yz = buffer_1100_sddd[514];

    auto g_x_z_0_0_0_xz_xy_zz = buffer_1100_sddd[515];

    auto g_x_z_0_0_0_xz_xz_xx = buffer_1100_sddd[516];

    auto g_x_z_0_0_0_xz_xz_xy = buffer_1100_sddd[517];

    auto g_x_z_0_0_0_xz_xz_xz = buffer_1100_sddd[518];

    auto g_x_z_0_0_0_xz_xz_yy = buffer_1100_sddd[519];

    auto g_x_z_0_0_0_xz_xz_yz = buffer_1100_sddd[520];

    auto g_x_z_0_0_0_xz_xz_zz = buffer_1100_sddd[521];

    auto g_x_z_0_0_0_xz_yy_xx = buffer_1100_sddd[522];

    auto g_x_z_0_0_0_xz_yy_xy = buffer_1100_sddd[523];

    auto g_x_z_0_0_0_xz_yy_xz = buffer_1100_sddd[524];

    auto g_x_z_0_0_0_xz_yy_yy = buffer_1100_sddd[525];

    auto g_x_z_0_0_0_xz_yy_yz = buffer_1100_sddd[526];

    auto g_x_z_0_0_0_xz_yy_zz = buffer_1100_sddd[527];

    auto g_x_z_0_0_0_xz_yz_xx = buffer_1100_sddd[528];

    auto g_x_z_0_0_0_xz_yz_xy = buffer_1100_sddd[529];

    auto g_x_z_0_0_0_xz_yz_xz = buffer_1100_sddd[530];

    auto g_x_z_0_0_0_xz_yz_yy = buffer_1100_sddd[531];

    auto g_x_z_0_0_0_xz_yz_yz = buffer_1100_sddd[532];

    auto g_x_z_0_0_0_xz_yz_zz = buffer_1100_sddd[533];

    auto g_x_z_0_0_0_xz_zz_xx = buffer_1100_sddd[534];

    auto g_x_z_0_0_0_xz_zz_xy = buffer_1100_sddd[535];

    auto g_x_z_0_0_0_xz_zz_xz = buffer_1100_sddd[536];

    auto g_x_z_0_0_0_xz_zz_yy = buffer_1100_sddd[537];

    auto g_x_z_0_0_0_xz_zz_yz = buffer_1100_sddd[538];

    auto g_x_z_0_0_0_xz_zz_zz = buffer_1100_sddd[539];

    auto g_x_z_0_0_0_yy_xx_xx = buffer_1100_sddd[540];

    auto g_x_z_0_0_0_yy_xx_xy = buffer_1100_sddd[541];

    auto g_x_z_0_0_0_yy_xx_xz = buffer_1100_sddd[542];

    auto g_x_z_0_0_0_yy_xx_yy = buffer_1100_sddd[543];

    auto g_x_z_0_0_0_yy_xx_yz = buffer_1100_sddd[544];

    auto g_x_z_0_0_0_yy_xx_zz = buffer_1100_sddd[545];

    auto g_x_z_0_0_0_yy_xy_xx = buffer_1100_sddd[546];

    auto g_x_z_0_0_0_yy_xy_xy = buffer_1100_sddd[547];

    auto g_x_z_0_0_0_yy_xy_xz = buffer_1100_sddd[548];

    auto g_x_z_0_0_0_yy_xy_yy = buffer_1100_sddd[549];

    auto g_x_z_0_0_0_yy_xy_yz = buffer_1100_sddd[550];

    auto g_x_z_0_0_0_yy_xy_zz = buffer_1100_sddd[551];

    auto g_x_z_0_0_0_yy_xz_xx = buffer_1100_sddd[552];

    auto g_x_z_0_0_0_yy_xz_xy = buffer_1100_sddd[553];

    auto g_x_z_0_0_0_yy_xz_xz = buffer_1100_sddd[554];

    auto g_x_z_0_0_0_yy_xz_yy = buffer_1100_sddd[555];

    auto g_x_z_0_0_0_yy_xz_yz = buffer_1100_sddd[556];

    auto g_x_z_0_0_0_yy_xz_zz = buffer_1100_sddd[557];

    auto g_x_z_0_0_0_yy_yy_xx = buffer_1100_sddd[558];

    auto g_x_z_0_0_0_yy_yy_xy = buffer_1100_sddd[559];

    auto g_x_z_0_0_0_yy_yy_xz = buffer_1100_sddd[560];

    auto g_x_z_0_0_0_yy_yy_yy = buffer_1100_sddd[561];

    auto g_x_z_0_0_0_yy_yy_yz = buffer_1100_sddd[562];

    auto g_x_z_0_0_0_yy_yy_zz = buffer_1100_sddd[563];

    auto g_x_z_0_0_0_yy_yz_xx = buffer_1100_sddd[564];

    auto g_x_z_0_0_0_yy_yz_xy = buffer_1100_sddd[565];

    auto g_x_z_0_0_0_yy_yz_xz = buffer_1100_sddd[566];

    auto g_x_z_0_0_0_yy_yz_yy = buffer_1100_sddd[567];

    auto g_x_z_0_0_0_yy_yz_yz = buffer_1100_sddd[568];

    auto g_x_z_0_0_0_yy_yz_zz = buffer_1100_sddd[569];

    auto g_x_z_0_0_0_yy_zz_xx = buffer_1100_sddd[570];

    auto g_x_z_0_0_0_yy_zz_xy = buffer_1100_sddd[571];

    auto g_x_z_0_0_0_yy_zz_xz = buffer_1100_sddd[572];

    auto g_x_z_0_0_0_yy_zz_yy = buffer_1100_sddd[573];

    auto g_x_z_0_0_0_yy_zz_yz = buffer_1100_sddd[574];

    auto g_x_z_0_0_0_yy_zz_zz = buffer_1100_sddd[575];

    auto g_x_z_0_0_0_yz_xx_xx = buffer_1100_sddd[576];

    auto g_x_z_0_0_0_yz_xx_xy = buffer_1100_sddd[577];

    auto g_x_z_0_0_0_yz_xx_xz = buffer_1100_sddd[578];

    auto g_x_z_0_0_0_yz_xx_yy = buffer_1100_sddd[579];

    auto g_x_z_0_0_0_yz_xx_yz = buffer_1100_sddd[580];

    auto g_x_z_0_0_0_yz_xx_zz = buffer_1100_sddd[581];

    auto g_x_z_0_0_0_yz_xy_xx = buffer_1100_sddd[582];

    auto g_x_z_0_0_0_yz_xy_xy = buffer_1100_sddd[583];

    auto g_x_z_0_0_0_yz_xy_xz = buffer_1100_sddd[584];

    auto g_x_z_0_0_0_yz_xy_yy = buffer_1100_sddd[585];

    auto g_x_z_0_0_0_yz_xy_yz = buffer_1100_sddd[586];

    auto g_x_z_0_0_0_yz_xy_zz = buffer_1100_sddd[587];

    auto g_x_z_0_0_0_yz_xz_xx = buffer_1100_sddd[588];

    auto g_x_z_0_0_0_yz_xz_xy = buffer_1100_sddd[589];

    auto g_x_z_0_0_0_yz_xz_xz = buffer_1100_sddd[590];

    auto g_x_z_0_0_0_yz_xz_yy = buffer_1100_sddd[591];

    auto g_x_z_0_0_0_yz_xz_yz = buffer_1100_sddd[592];

    auto g_x_z_0_0_0_yz_xz_zz = buffer_1100_sddd[593];

    auto g_x_z_0_0_0_yz_yy_xx = buffer_1100_sddd[594];

    auto g_x_z_0_0_0_yz_yy_xy = buffer_1100_sddd[595];

    auto g_x_z_0_0_0_yz_yy_xz = buffer_1100_sddd[596];

    auto g_x_z_0_0_0_yz_yy_yy = buffer_1100_sddd[597];

    auto g_x_z_0_0_0_yz_yy_yz = buffer_1100_sddd[598];

    auto g_x_z_0_0_0_yz_yy_zz = buffer_1100_sddd[599];

    auto g_x_z_0_0_0_yz_yz_xx = buffer_1100_sddd[600];

    auto g_x_z_0_0_0_yz_yz_xy = buffer_1100_sddd[601];

    auto g_x_z_0_0_0_yz_yz_xz = buffer_1100_sddd[602];

    auto g_x_z_0_0_0_yz_yz_yy = buffer_1100_sddd[603];

    auto g_x_z_0_0_0_yz_yz_yz = buffer_1100_sddd[604];

    auto g_x_z_0_0_0_yz_yz_zz = buffer_1100_sddd[605];

    auto g_x_z_0_0_0_yz_zz_xx = buffer_1100_sddd[606];

    auto g_x_z_0_0_0_yz_zz_xy = buffer_1100_sddd[607];

    auto g_x_z_0_0_0_yz_zz_xz = buffer_1100_sddd[608];

    auto g_x_z_0_0_0_yz_zz_yy = buffer_1100_sddd[609];

    auto g_x_z_0_0_0_yz_zz_yz = buffer_1100_sddd[610];

    auto g_x_z_0_0_0_yz_zz_zz = buffer_1100_sddd[611];

    auto g_x_z_0_0_0_zz_xx_xx = buffer_1100_sddd[612];

    auto g_x_z_0_0_0_zz_xx_xy = buffer_1100_sddd[613];

    auto g_x_z_0_0_0_zz_xx_xz = buffer_1100_sddd[614];

    auto g_x_z_0_0_0_zz_xx_yy = buffer_1100_sddd[615];

    auto g_x_z_0_0_0_zz_xx_yz = buffer_1100_sddd[616];

    auto g_x_z_0_0_0_zz_xx_zz = buffer_1100_sddd[617];

    auto g_x_z_0_0_0_zz_xy_xx = buffer_1100_sddd[618];

    auto g_x_z_0_0_0_zz_xy_xy = buffer_1100_sddd[619];

    auto g_x_z_0_0_0_zz_xy_xz = buffer_1100_sddd[620];

    auto g_x_z_0_0_0_zz_xy_yy = buffer_1100_sddd[621];

    auto g_x_z_0_0_0_zz_xy_yz = buffer_1100_sddd[622];

    auto g_x_z_0_0_0_zz_xy_zz = buffer_1100_sddd[623];

    auto g_x_z_0_0_0_zz_xz_xx = buffer_1100_sddd[624];

    auto g_x_z_0_0_0_zz_xz_xy = buffer_1100_sddd[625];

    auto g_x_z_0_0_0_zz_xz_xz = buffer_1100_sddd[626];

    auto g_x_z_0_0_0_zz_xz_yy = buffer_1100_sddd[627];

    auto g_x_z_0_0_0_zz_xz_yz = buffer_1100_sddd[628];

    auto g_x_z_0_0_0_zz_xz_zz = buffer_1100_sddd[629];

    auto g_x_z_0_0_0_zz_yy_xx = buffer_1100_sddd[630];

    auto g_x_z_0_0_0_zz_yy_xy = buffer_1100_sddd[631];

    auto g_x_z_0_0_0_zz_yy_xz = buffer_1100_sddd[632];

    auto g_x_z_0_0_0_zz_yy_yy = buffer_1100_sddd[633];

    auto g_x_z_0_0_0_zz_yy_yz = buffer_1100_sddd[634];

    auto g_x_z_0_0_0_zz_yy_zz = buffer_1100_sddd[635];

    auto g_x_z_0_0_0_zz_yz_xx = buffer_1100_sddd[636];

    auto g_x_z_0_0_0_zz_yz_xy = buffer_1100_sddd[637];

    auto g_x_z_0_0_0_zz_yz_xz = buffer_1100_sddd[638];

    auto g_x_z_0_0_0_zz_yz_yy = buffer_1100_sddd[639];

    auto g_x_z_0_0_0_zz_yz_yz = buffer_1100_sddd[640];

    auto g_x_z_0_0_0_zz_yz_zz = buffer_1100_sddd[641];

    auto g_x_z_0_0_0_zz_zz_xx = buffer_1100_sddd[642];

    auto g_x_z_0_0_0_zz_zz_xy = buffer_1100_sddd[643];

    auto g_x_z_0_0_0_zz_zz_xz = buffer_1100_sddd[644];

    auto g_x_z_0_0_0_zz_zz_yy = buffer_1100_sddd[645];

    auto g_x_z_0_0_0_zz_zz_yz = buffer_1100_sddd[646];

    auto g_x_z_0_0_0_zz_zz_zz = buffer_1100_sddd[647];

    auto g_y_x_0_0_0_xx_xx_xx = buffer_1100_sddd[648];

    auto g_y_x_0_0_0_xx_xx_xy = buffer_1100_sddd[649];

    auto g_y_x_0_0_0_xx_xx_xz = buffer_1100_sddd[650];

    auto g_y_x_0_0_0_xx_xx_yy = buffer_1100_sddd[651];

    auto g_y_x_0_0_0_xx_xx_yz = buffer_1100_sddd[652];

    auto g_y_x_0_0_0_xx_xx_zz = buffer_1100_sddd[653];

    auto g_y_x_0_0_0_xx_xy_xx = buffer_1100_sddd[654];

    auto g_y_x_0_0_0_xx_xy_xy = buffer_1100_sddd[655];

    auto g_y_x_0_0_0_xx_xy_xz = buffer_1100_sddd[656];

    auto g_y_x_0_0_0_xx_xy_yy = buffer_1100_sddd[657];

    auto g_y_x_0_0_0_xx_xy_yz = buffer_1100_sddd[658];

    auto g_y_x_0_0_0_xx_xy_zz = buffer_1100_sddd[659];

    auto g_y_x_0_0_0_xx_xz_xx = buffer_1100_sddd[660];

    auto g_y_x_0_0_0_xx_xz_xy = buffer_1100_sddd[661];

    auto g_y_x_0_0_0_xx_xz_xz = buffer_1100_sddd[662];

    auto g_y_x_0_0_0_xx_xz_yy = buffer_1100_sddd[663];

    auto g_y_x_0_0_0_xx_xz_yz = buffer_1100_sddd[664];

    auto g_y_x_0_0_0_xx_xz_zz = buffer_1100_sddd[665];

    auto g_y_x_0_0_0_xx_yy_xx = buffer_1100_sddd[666];

    auto g_y_x_0_0_0_xx_yy_xy = buffer_1100_sddd[667];

    auto g_y_x_0_0_0_xx_yy_xz = buffer_1100_sddd[668];

    auto g_y_x_0_0_0_xx_yy_yy = buffer_1100_sddd[669];

    auto g_y_x_0_0_0_xx_yy_yz = buffer_1100_sddd[670];

    auto g_y_x_0_0_0_xx_yy_zz = buffer_1100_sddd[671];

    auto g_y_x_0_0_0_xx_yz_xx = buffer_1100_sddd[672];

    auto g_y_x_0_0_0_xx_yz_xy = buffer_1100_sddd[673];

    auto g_y_x_0_0_0_xx_yz_xz = buffer_1100_sddd[674];

    auto g_y_x_0_0_0_xx_yz_yy = buffer_1100_sddd[675];

    auto g_y_x_0_0_0_xx_yz_yz = buffer_1100_sddd[676];

    auto g_y_x_0_0_0_xx_yz_zz = buffer_1100_sddd[677];

    auto g_y_x_0_0_0_xx_zz_xx = buffer_1100_sddd[678];

    auto g_y_x_0_0_0_xx_zz_xy = buffer_1100_sddd[679];

    auto g_y_x_0_0_0_xx_zz_xz = buffer_1100_sddd[680];

    auto g_y_x_0_0_0_xx_zz_yy = buffer_1100_sddd[681];

    auto g_y_x_0_0_0_xx_zz_yz = buffer_1100_sddd[682];

    auto g_y_x_0_0_0_xx_zz_zz = buffer_1100_sddd[683];

    auto g_y_x_0_0_0_xy_xx_xx = buffer_1100_sddd[684];

    auto g_y_x_0_0_0_xy_xx_xy = buffer_1100_sddd[685];

    auto g_y_x_0_0_0_xy_xx_xz = buffer_1100_sddd[686];

    auto g_y_x_0_0_0_xy_xx_yy = buffer_1100_sddd[687];

    auto g_y_x_0_0_0_xy_xx_yz = buffer_1100_sddd[688];

    auto g_y_x_0_0_0_xy_xx_zz = buffer_1100_sddd[689];

    auto g_y_x_0_0_0_xy_xy_xx = buffer_1100_sddd[690];

    auto g_y_x_0_0_0_xy_xy_xy = buffer_1100_sddd[691];

    auto g_y_x_0_0_0_xy_xy_xz = buffer_1100_sddd[692];

    auto g_y_x_0_0_0_xy_xy_yy = buffer_1100_sddd[693];

    auto g_y_x_0_0_0_xy_xy_yz = buffer_1100_sddd[694];

    auto g_y_x_0_0_0_xy_xy_zz = buffer_1100_sddd[695];

    auto g_y_x_0_0_0_xy_xz_xx = buffer_1100_sddd[696];

    auto g_y_x_0_0_0_xy_xz_xy = buffer_1100_sddd[697];

    auto g_y_x_0_0_0_xy_xz_xz = buffer_1100_sddd[698];

    auto g_y_x_0_0_0_xy_xz_yy = buffer_1100_sddd[699];

    auto g_y_x_0_0_0_xy_xz_yz = buffer_1100_sddd[700];

    auto g_y_x_0_0_0_xy_xz_zz = buffer_1100_sddd[701];

    auto g_y_x_0_0_0_xy_yy_xx = buffer_1100_sddd[702];

    auto g_y_x_0_0_0_xy_yy_xy = buffer_1100_sddd[703];

    auto g_y_x_0_0_0_xy_yy_xz = buffer_1100_sddd[704];

    auto g_y_x_0_0_0_xy_yy_yy = buffer_1100_sddd[705];

    auto g_y_x_0_0_0_xy_yy_yz = buffer_1100_sddd[706];

    auto g_y_x_0_0_0_xy_yy_zz = buffer_1100_sddd[707];

    auto g_y_x_0_0_0_xy_yz_xx = buffer_1100_sddd[708];

    auto g_y_x_0_0_0_xy_yz_xy = buffer_1100_sddd[709];

    auto g_y_x_0_0_0_xy_yz_xz = buffer_1100_sddd[710];

    auto g_y_x_0_0_0_xy_yz_yy = buffer_1100_sddd[711];

    auto g_y_x_0_0_0_xy_yz_yz = buffer_1100_sddd[712];

    auto g_y_x_0_0_0_xy_yz_zz = buffer_1100_sddd[713];

    auto g_y_x_0_0_0_xy_zz_xx = buffer_1100_sddd[714];

    auto g_y_x_0_0_0_xy_zz_xy = buffer_1100_sddd[715];

    auto g_y_x_0_0_0_xy_zz_xz = buffer_1100_sddd[716];

    auto g_y_x_0_0_0_xy_zz_yy = buffer_1100_sddd[717];

    auto g_y_x_0_0_0_xy_zz_yz = buffer_1100_sddd[718];

    auto g_y_x_0_0_0_xy_zz_zz = buffer_1100_sddd[719];

    auto g_y_x_0_0_0_xz_xx_xx = buffer_1100_sddd[720];

    auto g_y_x_0_0_0_xz_xx_xy = buffer_1100_sddd[721];

    auto g_y_x_0_0_0_xz_xx_xz = buffer_1100_sddd[722];

    auto g_y_x_0_0_0_xz_xx_yy = buffer_1100_sddd[723];

    auto g_y_x_0_0_0_xz_xx_yz = buffer_1100_sddd[724];

    auto g_y_x_0_0_0_xz_xx_zz = buffer_1100_sddd[725];

    auto g_y_x_0_0_0_xz_xy_xx = buffer_1100_sddd[726];

    auto g_y_x_0_0_0_xz_xy_xy = buffer_1100_sddd[727];

    auto g_y_x_0_0_0_xz_xy_xz = buffer_1100_sddd[728];

    auto g_y_x_0_0_0_xz_xy_yy = buffer_1100_sddd[729];

    auto g_y_x_0_0_0_xz_xy_yz = buffer_1100_sddd[730];

    auto g_y_x_0_0_0_xz_xy_zz = buffer_1100_sddd[731];

    auto g_y_x_0_0_0_xz_xz_xx = buffer_1100_sddd[732];

    auto g_y_x_0_0_0_xz_xz_xy = buffer_1100_sddd[733];

    auto g_y_x_0_0_0_xz_xz_xz = buffer_1100_sddd[734];

    auto g_y_x_0_0_0_xz_xz_yy = buffer_1100_sddd[735];

    auto g_y_x_0_0_0_xz_xz_yz = buffer_1100_sddd[736];

    auto g_y_x_0_0_0_xz_xz_zz = buffer_1100_sddd[737];

    auto g_y_x_0_0_0_xz_yy_xx = buffer_1100_sddd[738];

    auto g_y_x_0_0_0_xz_yy_xy = buffer_1100_sddd[739];

    auto g_y_x_0_0_0_xz_yy_xz = buffer_1100_sddd[740];

    auto g_y_x_0_0_0_xz_yy_yy = buffer_1100_sddd[741];

    auto g_y_x_0_0_0_xz_yy_yz = buffer_1100_sddd[742];

    auto g_y_x_0_0_0_xz_yy_zz = buffer_1100_sddd[743];

    auto g_y_x_0_0_0_xz_yz_xx = buffer_1100_sddd[744];

    auto g_y_x_0_0_0_xz_yz_xy = buffer_1100_sddd[745];

    auto g_y_x_0_0_0_xz_yz_xz = buffer_1100_sddd[746];

    auto g_y_x_0_0_0_xz_yz_yy = buffer_1100_sddd[747];

    auto g_y_x_0_0_0_xz_yz_yz = buffer_1100_sddd[748];

    auto g_y_x_0_0_0_xz_yz_zz = buffer_1100_sddd[749];

    auto g_y_x_0_0_0_xz_zz_xx = buffer_1100_sddd[750];

    auto g_y_x_0_0_0_xz_zz_xy = buffer_1100_sddd[751];

    auto g_y_x_0_0_0_xz_zz_xz = buffer_1100_sddd[752];

    auto g_y_x_0_0_0_xz_zz_yy = buffer_1100_sddd[753];

    auto g_y_x_0_0_0_xz_zz_yz = buffer_1100_sddd[754];

    auto g_y_x_0_0_0_xz_zz_zz = buffer_1100_sddd[755];

    auto g_y_x_0_0_0_yy_xx_xx = buffer_1100_sddd[756];

    auto g_y_x_0_0_0_yy_xx_xy = buffer_1100_sddd[757];

    auto g_y_x_0_0_0_yy_xx_xz = buffer_1100_sddd[758];

    auto g_y_x_0_0_0_yy_xx_yy = buffer_1100_sddd[759];

    auto g_y_x_0_0_0_yy_xx_yz = buffer_1100_sddd[760];

    auto g_y_x_0_0_0_yy_xx_zz = buffer_1100_sddd[761];

    auto g_y_x_0_0_0_yy_xy_xx = buffer_1100_sddd[762];

    auto g_y_x_0_0_0_yy_xy_xy = buffer_1100_sddd[763];

    auto g_y_x_0_0_0_yy_xy_xz = buffer_1100_sddd[764];

    auto g_y_x_0_0_0_yy_xy_yy = buffer_1100_sddd[765];

    auto g_y_x_0_0_0_yy_xy_yz = buffer_1100_sddd[766];

    auto g_y_x_0_0_0_yy_xy_zz = buffer_1100_sddd[767];

    auto g_y_x_0_0_0_yy_xz_xx = buffer_1100_sddd[768];

    auto g_y_x_0_0_0_yy_xz_xy = buffer_1100_sddd[769];

    auto g_y_x_0_0_0_yy_xz_xz = buffer_1100_sddd[770];

    auto g_y_x_0_0_0_yy_xz_yy = buffer_1100_sddd[771];

    auto g_y_x_0_0_0_yy_xz_yz = buffer_1100_sddd[772];

    auto g_y_x_0_0_0_yy_xz_zz = buffer_1100_sddd[773];

    auto g_y_x_0_0_0_yy_yy_xx = buffer_1100_sddd[774];

    auto g_y_x_0_0_0_yy_yy_xy = buffer_1100_sddd[775];

    auto g_y_x_0_0_0_yy_yy_xz = buffer_1100_sddd[776];

    auto g_y_x_0_0_0_yy_yy_yy = buffer_1100_sddd[777];

    auto g_y_x_0_0_0_yy_yy_yz = buffer_1100_sddd[778];

    auto g_y_x_0_0_0_yy_yy_zz = buffer_1100_sddd[779];

    auto g_y_x_0_0_0_yy_yz_xx = buffer_1100_sddd[780];

    auto g_y_x_0_0_0_yy_yz_xy = buffer_1100_sddd[781];

    auto g_y_x_0_0_0_yy_yz_xz = buffer_1100_sddd[782];

    auto g_y_x_0_0_0_yy_yz_yy = buffer_1100_sddd[783];

    auto g_y_x_0_0_0_yy_yz_yz = buffer_1100_sddd[784];

    auto g_y_x_0_0_0_yy_yz_zz = buffer_1100_sddd[785];

    auto g_y_x_0_0_0_yy_zz_xx = buffer_1100_sddd[786];

    auto g_y_x_0_0_0_yy_zz_xy = buffer_1100_sddd[787];

    auto g_y_x_0_0_0_yy_zz_xz = buffer_1100_sddd[788];

    auto g_y_x_0_0_0_yy_zz_yy = buffer_1100_sddd[789];

    auto g_y_x_0_0_0_yy_zz_yz = buffer_1100_sddd[790];

    auto g_y_x_0_0_0_yy_zz_zz = buffer_1100_sddd[791];

    auto g_y_x_0_0_0_yz_xx_xx = buffer_1100_sddd[792];

    auto g_y_x_0_0_0_yz_xx_xy = buffer_1100_sddd[793];

    auto g_y_x_0_0_0_yz_xx_xz = buffer_1100_sddd[794];

    auto g_y_x_0_0_0_yz_xx_yy = buffer_1100_sddd[795];

    auto g_y_x_0_0_0_yz_xx_yz = buffer_1100_sddd[796];

    auto g_y_x_0_0_0_yz_xx_zz = buffer_1100_sddd[797];

    auto g_y_x_0_0_0_yz_xy_xx = buffer_1100_sddd[798];

    auto g_y_x_0_0_0_yz_xy_xy = buffer_1100_sddd[799];

    auto g_y_x_0_0_0_yz_xy_xz = buffer_1100_sddd[800];

    auto g_y_x_0_0_0_yz_xy_yy = buffer_1100_sddd[801];

    auto g_y_x_0_0_0_yz_xy_yz = buffer_1100_sddd[802];

    auto g_y_x_0_0_0_yz_xy_zz = buffer_1100_sddd[803];

    auto g_y_x_0_0_0_yz_xz_xx = buffer_1100_sddd[804];

    auto g_y_x_0_0_0_yz_xz_xy = buffer_1100_sddd[805];

    auto g_y_x_0_0_0_yz_xz_xz = buffer_1100_sddd[806];

    auto g_y_x_0_0_0_yz_xz_yy = buffer_1100_sddd[807];

    auto g_y_x_0_0_0_yz_xz_yz = buffer_1100_sddd[808];

    auto g_y_x_0_0_0_yz_xz_zz = buffer_1100_sddd[809];

    auto g_y_x_0_0_0_yz_yy_xx = buffer_1100_sddd[810];

    auto g_y_x_0_0_0_yz_yy_xy = buffer_1100_sddd[811];

    auto g_y_x_0_0_0_yz_yy_xz = buffer_1100_sddd[812];

    auto g_y_x_0_0_0_yz_yy_yy = buffer_1100_sddd[813];

    auto g_y_x_0_0_0_yz_yy_yz = buffer_1100_sddd[814];

    auto g_y_x_0_0_0_yz_yy_zz = buffer_1100_sddd[815];

    auto g_y_x_0_0_0_yz_yz_xx = buffer_1100_sddd[816];

    auto g_y_x_0_0_0_yz_yz_xy = buffer_1100_sddd[817];

    auto g_y_x_0_0_0_yz_yz_xz = buffer_1100_sddd[818];

    auto g_y_x_0_0_0_yz_yz_yy = buffer_1100_sddd[819];

    auto g_y_x_0_0_0_yz_yz_yz = buffer_1100_sddd[820];

    auto g_y_x_0_0_0_yz_yz_zz = buffer_1100_sddd[821];

    auto g_y_x_0_0_0_yz_zz_xx = buffer_1100_sddd[822];

    auto g_y_x_0_0_0_yz_zz_xy = buffer_1100_sddd[823];

    auto g_y_x_0_0_0_yz_zz_xz = buffer_1100_sddd[824];

    auto g_y_x_0_0_0_yz_zz_yy = buffer_1100_sddd[825];

    auto g_y_x_0_0_0_yz_zz_yz = buffer_1100_sddd[826];

    auto g_y_x_0_0_0_yz_zz_zz = buffer_1100_sddd[827];

    auto g_y_x_0_0_0_zz_xx_xx = buffer_1100_sddd[828];

    auto g_y_x_0_0_0_zz_xx_xy = buffer_1100_sddd[829];

    auto g_y_x_0_0_0_zz_xx_xz = buffer_1100_sddd[830];

    auto g_y_x_0_0_0_zz_xx_yy = buffer_1100_sddd[831];

    auto g_y_x_0_0_0_zz_xx_yz = buffer_1100_sddd[832];

    auto g_y_x_0_0_0_zz_xx_zz = buffer_1100_sddd[833];

    auto g_y_x_0_0_0_zz_xy_xx = buffer_1100_sddd[834];

    auto g_y_x_0_0_0_zz_xy_xy = buffer_1100_sddd[835];

    auto g_y_x_0_0_0_zz_xy_xz = buffer_1100_sddd[836];

    auto g_y_x_0_0_0_zz_xy_yy = buffer_1100_sddd[837];

    auto g_y_x_0_0_0_zz_xy_yz = buffer_1100_sddd[838];

    auto g_y_x_0_0_0_zz_xy_zz = buffer_1100_sddd[839];

    auto g_y_x_0_0_0_zz_xz_xx = buffer_1100_sddd[840];

    auto g_y_x_0_0_0_zz_xz_xy = buffer_1100_sddd[841];

    auto g_y_x_0_0_0_zz_xz_xz = buffer_1100_sddd[842];

    auto g_y_x_0_0_0_zz_xz_yy = buffer_1100_sddd[843];

    auto g_y_x_0_0_0_zz_xz_yz = buffer_1100_sddd[844];

    auto g_y_x_0_0_0_zz_xz_zz = buffer_1100_sddd[845];

    auto g_y_x_0_0_0_zz_yy_xx = buffer_1100_sddd[846];

    auto g_y_x_0_0_0_zz_yy_xy = buffer_1100_sddd[847];

    auto g_y_x_0_0_0_zz_yy_xz = buffer_1100_sddd[848];

    auto g_y_x_0_0_0_zz_yy_yy = buffer_1100_sddd[849];

    auto g_y_x_0_0_0_zz_yy_yz = buffer_1100_sddd[850];

    auto g_y_x_0_0_0_zz_yy_zz = buffer_1100_sddd[851];

    auto g_y_x_0_0_0_zz_yz_xx = buffer_1100_sddd[852];

    auto g_y_x_0_0_0_zz_yz_xy = buffer_1100_sddd[853];

    auto g_y_x_0_0_0_zz_yz_xz = buffer_1100_sddd[854];

    auto g_y_x_0_0_0_zz_yz_yy = buffer_1100_sddd[855];

    auto g_y_x_0_0_0_zz_yz_yz = buffer_1100_sddd[856];

    auto g_y_x_0_0_0_zz_yz_zz = buffer_1100_sddd[857];

    auto g_y_x_0_0_0_zz_zz_xx = buffer_1100_sddd[858];

    auto g_y_x_0_0_0_zz_zz_xy = buffer_1100_sddd[859];

    auto g_y_x_0_0_0_zz_zz_xz = buffer_1100_sddd[860];

    auto g_y_x_0_0_0_zz_zz_yy = buffer_1100_sddd[861];

    auto g_y_x_0_0_0_zz_zz_yz = buffer_1100_sddd[862];

    auto g_y_x_0_0_0_zz_zz_zz = buffer_1100_sddd[863];

    auto g_y_y_0_0_0_xx_xx_xx = buffer_1100_sddd[864];

    auto g_y_y_0_0_0_xx_xx_xy = buffer_1100_sddd[865];

    auto g_y_y_0_0_0_xx_xx_xz = buffer_1100_sddd[866];

    auto g_y_y_0_0_0_xx_xx_yy = buffer_1100_sddd[867];

    auto g_y_y_0_0_0_xx_xx_yz = buffer_1100_sddd[868];

    auto g_y_y_0_0_0_xx_xx_zz = buffer_1100_sddd[869];

    auto g_y_y_0_0_0_xx_xy_xx = buffer_1100_sddd[870];

    auto g_y_y_0_0_0_xx_xy_xy = buffer_1100_sddd[871];

    auto g_y_y_0_0_0_xx_xy_xz = buffer_1100_sddd[872];

    auto g_y_y_0_0_0_xx_xy_yy = buffer_1100_sddd[873];

    auto g_y_y_0_0_0_xx_xy_yz = buffer_1100_sddd[874];

    auto g_y_y_0_0_0_xx_xy_zz = buffer_1100_sddd[875];

    auto g_y_y_0_0_0_xx_xz_xx = buffer_1100_sddd[876];

    auto g_y_y_0_0_0_xx_xz_xy = buffer_1100_sddd[877];

    auto g_y_y_0_0_0_xx_xz_xz = buffer_1100_sddd[878];

    auto g_y_y_0_0_0_xx_xz_yy = buffer_1100_sddd[879];

    auto g_y_y_0_0_0_xx_xz_yz = buffer_1100_sddd[880];

    auto g_y_y_0_0_0_xx_xz_zz = buffer_1100_sddd[881];

    auto g_y_y_0_0_0_xx_yy_xx = buffer_1100_sddd[882];

    auto g_y_y_0_0_0_xx_yy_xy = buffer_1100_sddd[883];

    auto g_y_y_0_0_0_xx_yy_xz = buffer_1100_sddd[884];

    auto g_y_y_0_0_0_xx_yy_yy = buffer_1100_sddd[885];

    auto g_y_y_0_0_0_xx_yy_yz = buffer_1100_sddd[886];

    auto g_y_y_0_0_0_xx_yy_zz = buffer_1100_sddd[887];

    auto g_y_y_0_0_0_xx_yz_xx = buffer_1100_sddd[888];

    auto g_y_y_0_0_0_xx_yz_xy = buffer_1100_sddd[889];

    auto g_y_y_0_0_0_xx_yz_xz = buffer_1100_sddd[890];

    auto g_y_y_0_0_0_xx_yz_yy = buffer_1100_sddd[891];

    auto g_y_y_0_0_0_xx_yz_yz = buffer_1100_sddd[892];

    auto g_y_y_0_0_0_xx_yz_zz = buffer_1100_sddd[893];

    auto g_y_y_0_0_0_xx_zz_xx = buffer_1100_sddd[894];

    auto g_y_y_0_0_0_xx_zz_xy = buffer_1100_sddd[895];

    auto g_y_y_0_0_0_xx_zz_xz = buffer_1100_sddd[896];

    auto g_y_y_0_0_0_xx_zz_yy = buffer_1100_sddd[897];

    auto g_y_y_0_0_0_xx_zz_yz = buffer_1100_sddd[898];

    auto g_y_y_0_0_0_xx_zz_zz = buffer_1100_sddd[899];

    auto g_y_y_0_0_0_xy_xx_xx = buffer_1100_sddd[900];

    auto g_y_y_0_0_0_xy_xx_xy = buffer_1100_sddd[901];

    auto g_y_y_0_0_0_xy_xx_xz = buffer_1100_sddd[902];

    auto g_y_y_0_0_0_xy_xx_yy = buffer_1100_sddd[903];

    auto g_y_y_0_0_0_xy_xx_yz = buffer_1100_sddd[904];

    auto g_y_y_0_0_0_xy_xx_zz = buffer_1100_sddd[905];

    auto g_y_y_0_0_0_xy_xy_xx = buffer_1100_sddd[906];

    auto g_y_y_0_0_0_xy_xy_xy = buffer_1100_sddd[907];

    auto g_y_y_0_0_0_xy_xy_xz = buffer_1100_sddd[908];

    auto g_y_y_0_0_0_xy_xy_yy = buffer_1100_sddd[909];

    auto g_y_y_0_0_0_xy_xy_yz = buffer_1100_sddd[910];

    auto g_y_y_0_0_0_xy_xy_zz = buffer_1100_sddd[911];

    auto g_y_y_0_0_0_xy_xz_xx = buffer_1100_sddd[912];

    auto g_y_y_0_0_0_xy_xz_xy = buffer_1100_sddd[913];

    auto g_y_y_0_0_0_xy_xz_xz = buffer_1100_sddd[914];

    auto g_y_y_0_0_0_xy_xz_yy = buffer_1100_sddd[915];

    auto g_y_y_0_0_0_xy_xz_yz = buffer_1100_sddd[916];

    auto g_y_y_0_0_0_xy_xz_zz = buffer_1100_sddd[917];

    auto g_y_y_0_0_0_xy_yy_xx = buffer_1100_sddd[918];

    auto g_y_y_0_0_0_xy_yy_xy = buffer_1100_sddd[919];

    auto g_y_y_0_0_0_xy_yy_xz = buffer_1100_sddd[920];

    auto g_y_y_0_0_0_xy_yy_yy = buffer_1100_sddd[921];

    auto g_y_y_0_0_0_xy_yy_yz = buffer_1100_sddd[922];

    auto g_y_y_0_0_0_xy_yy_zz = buffer_1100_sddd[923];

    auto g_y_y_0_0_0_xy_yz_xx = buffer_1100_sddd[924];

    auto g_y_y_0_0_0_xy_yz_xy = buffer_1100_sddd[925];

    auto g_y_y_0_0_0_xy_yz_xz = buffer_1100_sddd[926];

    auto g_y_y_0_0_0_xy_yz_yy = buffer_1100_sddd[927];

    auto g_y_y_0_0_0_xy_yz_yz = buffer_1100_sddd[928];

    auto g_y_y_0_0_0_xy_yz_zz = buffer_1100_sddd[929];

    auto g_y_y_0_0_0_xy_zz_xx = buffer_1100_sddd[930];

    auto g_y_y_0_0_0_xy_zz_xy = buffer_1100_sddd[931];

    auto g_y_y_0_0_0_xy_zz_xz = buffer_1100_sddd[932];

    auto g_y_y_0_0_0_xy_zz_yy = buffer_1100_sddd[933];

    auto g_y_y_0_0_0_xy_zz_yz = buffer_1100_sddd[934];

    auto g_y_y_0_0_0_xy_zz_zz = buffer_1100_sddd[935];

    auto g_y_y_0_0_0_xz_xx_xx = buffer_1100_sddd[936];

    auto g_y_y_0_0_0_xz_xx_xy = buffer_1100_sddd[937];

    auto g_y_y_0_0_0_xz_xx_xz = buffer_1100_sddd[938];

    auto g_y_y_0_0_0_xz_xx_yy = buffer_1100_sddd[939];

    auto g_y_y_0_0_0_xz_xx_yz = buffer_1100_sddd[940];

    auto g_y_y_0_0_0_xz_xx_zz = buffer_1100_sddd[941];

    auto g_y_y_0_0_0_xz_xy_xx = buffer_1100_sddd[942];

    auto g_y_y_0_0_0_xz_xy_xy = buffer_1100_sddd[943];

    auto g_y_y_0_0_0_xz_xy_xz = buffer_1100_sddd[944];

    auto g_y_y_0_0_0_xz_xy_yy = buffer_1100_sddd[945];

    auto g_y_y_0_0_0_xz_xy_yz = buffer_1100_sddd[946];

    auto g_y_y_0_0_0_xz_xy_zz = buffer_1100_sddd[947];

    auto g_y_y_0_0_0_xz_xz_xx = buffer_1100_sddd[948];

    auto g_y_y_0_0_0_xz_xz_xy = buffer_1100_sddd[949];

    auto g_y_y_0_0_0_xz_xz_xz = buffer_1100_sddd[950];

    auto g_y_y_0_0_0_xz_xz_yy = buffer_1100_sddd[951];

    auto g_y_y_0_0_0_xz_xz_yz = buffer_1100_sddd[952];

    auto g_y_y_0_0_0_xz_xz_zz = buffer_1100_sddd[953];

    auto g_y_y_0_0_0_xz_yy_xx = buffer_1100_sddd[954];

    auto g_y_y_0_0_0_xz_yy_xy = buffer_1100_sddd[955];

    auto g_y_y_0_0_0_xz_yy_xz = buffer_1100_sddd[956];

    auto g_y_y_0_0_0_xz_yy_yy = buffer_1100_sddd[957];

    auto g_y_y_0_0_0_xz_yy_yz = buffer_1100_sddd[958];

    auto g_y_y_0_0_0_xz_yy_zz = buffer_1100_sddd[959];

    auto g_y_y_0_0_0_xz_yz_xx = buffer_1100_sddd[960];

    auto g_y_y_0_0_0_xz_yz_xy = buffer_1100_sddd[961];

    auto g_y_y_0_0_0_xz_yz_xz = buffer_1100_sddd[962];

    auto g_y_y_0_0_0_xz_yz_yy = buffer_1100_sddd[963];

    auto g_y_y_0_0_0_xz_yz_yz = buffer_1100_sddd[964];

    auto g_y_y_0_0_0_xz_yz_zz = buffer_1100_sddd[965];

    auto g_y_y_0_0_0_xz_zz_xx = buffer_1100_sddd[966];

    auto g_y_y_0_0_0_xz_zz_xy = buffer_1100_sddd[967];

    auto g_y_y_0_0_0_xz_zz_xz = buffer_1100_sddd[968];

    auto g_y_y_0_0_0_xz_zz_yy = buffer_1100_sddd[969];

    auto g_y_y_0_0_0_xz_zz_yz = buffer_1100_sddd[970];

    auto g_y_y_0_0_0_xz_zz_zz = buffer_1100_sddd[971];

    auto g_y_y_0_0_0_yy_xx_xx = buffer_1100_sddd[972];

    auto g_y_y_0_0_0_yy_xx_xy = buffer_1100_sddd[973];

    auto g_y_y_0_0_0_yy_xx_xz = buffer_1100_sddd[974];

    auto g_y_y_0_0_0_yy_xx_yy = buffer_1100_sddd[975];

    auto g_y_y_0_0_0_yy_xx_yz = buffer_1100_sddd[976];

    auto g_y_y_0_0_0_yy_xx_zz = buffer_1100_sddd[977];

    auto g_y_y_0_0_0_yy_xy_xx = buffer_1100_sddd[978];

    auto g_y_y_0_0_0_yy_xy_xy = buffer_1100_sddd[979];

    auto g_y_y_0_0_0_yy_xy_xz = buffer_1100_sddd[980];

    auto g_y_y_0_0_0_yy_xy_yy = buffer_1100_sddd[981];

    auto g_y_y_0_0_0_yy_xy_yz = buffer_1100_sddd[982];

    auto g_y_y_0_0_0_yy_xy_zz = buffer_1100_sddd[983];

    auto g_y_y_0_0_0_yy_xz_xx = buffer_1100_sddd[984];

    auto g_y_y_0_0_0_yy_xz_xy = buffer_1100_sddd[985];

    auto g_y_y_0_0_0_yy_xz_xz = buffer_1100_sddd[986];

    auto g_y_y_0_0_0_yy_xz_yy = buffer_1100_sddd[987];

    auto g_y_y_0_0_0_yy_xz_yz = buffer_1100_sddd[988];

    auto g_y_y_0_0_0_yy_xz_zz = buffer_1100_sddd[989];

    auto g_y_y_0_0_0_yy_yy_xx = buffer_1100_sddd[990];

    auto g_y_y_0_0_0_yy_yy_xy = buffer_1100_sddd[991];

    auto g_y_y_0_0_0_yy_yy_xz = buffer_1100_sddd[992];

    auto g_y_y_0_0_0_yy_yy_yy = buffer_1100_sddd[993];

    auto g_y_y_0_0_0_yy_yy_yz = buffer_1100_sddd[994];

    auto g_y_y_0_0_0_yy_yy_zz = buffer_1100_sddd[995];

    auto g_y_y_0_0_0_yy_yz_xx = buffer_1100_sddd[996];

    auto g_y_y_0_0_0_yy_yz_xy = buffer_1100_sddd[997];

    auto g_y_y_0_0_0_yy_yz_xz = buffer_1100_sddd[998];

    auto g_y_y_0_0_0_yy_yz_yy = buffer_1100_sddd[999];

    auto g_y_y_0_0_0_yy_yz_yz = buffer_1100_sddd[1000];

    auto g_y_y_0_0_0_yy_yz_zz = buffer_1100_sddd[1001];

    auto g_y_y_0_0_0_yy_zz_xx = buffer_1100_sddd[1002];

    auto g_y_y_0_0_0_yy_zz_xy = buffer_1100_sddd[1003];

    auto g_y_y_0_0_0_yy_zz_xz = buffer_1100_sddd[1004];

    auto g_y_y_0_0_0_yy_zz_yy = buffer_1100_sddd[1005];

    auto g_y_y_0_0_0_yy_zz_yz = buffer_1100_sddd[1006];

    auto g_y_y_0_0_0_yy_zz_zz = buffer_1100_sddd[1007];

    auto g_y_y_0_0_0_yz_xx_xx = buffer_1100_sddd[1008];

    auto g_y_y_0_0_0_yz_xx_xy = buffer_1100_sddd[1009];

    auto g_y_y_0_0_0_yz_xx_xz = buffer_1100_sddd[1010];

    auto g_y_y_0_0_0_yz_xx_yy = buffer_1100_sddd[1011];

    auto g_y_y_0_0_0_yz_xx_yz = buffer_1100_sddd[1012];

    auto g_y_y_0_0_0_yz_xx_zz = buffer_1100_sddd[1013];

    auto g_y_y_0_0_0_yz_xy_xx = buffer_1100_sddd[1014];

    auto g_y_y_0_0_0_yz_xy_xy = buffer_1100_sddd[1015];

    auto g_y_y_0_0_0_yz_xy_xz = buffer_1100_sddd[1016];

    auto g_y_y_0_0_0_yz_xy_yy = buffer_1100_sddd[1017];

    auto g_y_y_0_0_0_yz_xy_yz = buffer_1100_sddd[1018];

    auto g_y_y_0_0_0_yz_xy_zz = buffer_1100_sddd[1019];

    auto g_y_y_0_0_0_yz_xz_xx = buffer_1100_sddd[1020];

    auto g_y_y_0_0_0_yz_xz_xy = buffer_1100_sddd[1021];

    auto g_y_y_0_0_0_yz_xz_xz = buffer_1100_sddd[1022];

    auto g_y_y_0_0_0_yz_xz_yy = buffer_1100_sddd[1023];

    auto g_y_y_0_0_0_yz_xz_yz = buffer_1100_sddd[1024];

    auto g_y_y_0_0_0_yz_xz_zz = buffer_1100_sddd[1025];

    auto g_y_y_0_0_0_yz_yy_xx = buffer_1100_sddd[1026];

    auto g_y_y_0_0_0_yz_yy_xy = buffer_1100_sddd[1027];

    auto g_y_y_0_0_0_yz_yy_xz = buffer_1100_sddd[1028];

    auto g_y_y_0_0_0_yz_yy_yy = buffer_1100_sddd[1029];

    auto g_y_y_0_0_0_yz_yy_yz = buffer_1100_sddd[1030];

    auto g_y_y_0_0_0_yz_yy_zz = buffer_1100_sddd[1031];

    auto g_y_y_0_0_0_yz_yz_xx = buffer_1100_sddd[1032];

    auto g_y_y_0_0_0_yz_yz_xy = buffer_1100_sddd[1033];

    auto g_y_y_0_0_0_yz_yz_xz = buffer_1100_sddd[1034];

    auto g_y_y_0_0_0_yz_yz_yy = buffer_1100_sddd[1035];

    auto g_y_y_0_0_0_yz_yz_yz = buffer_1100_sddd[1036];

    auto g_y_y_0_0_0_yz_yz_zz = buffer_1100_sddd[1037];

    auto g_y_y_0_0_0_yz_zz_xx = buffer_1100_sddd[1038];

    auto g_y_y_0_0_0_yz_zz_xy = buffer_1100_sddd[1039];

    auto g_y_y_0_0_0_yz_zz_xz = buffer_1100_sddd[1040];

    auto g_y_y_0_0_0_yz_zz_yy = buffer_1100_sddd[1041];

    auto g_y_y_0_0_0_yz_zz_yz = buffer_1100_sddd[1042];

    auto g_y_y_0_0_0_yz_zz_zz = buffer_1100_sddd[1043];

    auto g_y_y_0_0_0_zz_xx_xx = buffer_1100_sddd[1044];

    auto g_y_y_0_0_0_zz_xx_xy = buffer_1100_sddd[1045];

    auto g_y_y_0_0_0_zz_xx_xz = buffer_1100_sddd[1046];

    auto g_y_y_0_0_0_zz_xx_yy = buffer_1100_sddd[1047];

    auto g_y_y_0_0_0_zz_xx_yz = buffer_1100_sddd[1048];

    auto g_y_y_0_0_0_zz_xx_zz = buffer_1100_sddd[1049];

    auto g_y_y_0_0_0_zz_xy_xx = buffer_1100_sddd[1050];

    auto g_y_y_0_0_0_zz_xy_xy = buffer_1100_sddd[1051];

    auto g_y_y_0_0_0_zz_xy_xz = buffer_1100_sddd[1052];

    auto g_y_y_0_0_0_zz_xy_yy = buffer_1100_sddd[1053];

    auto g_y_y_0_0_0_zz_xy_yz = buffer_1100_sddd[1054];

    auto g_y_y_0_0_0_zz_xy_zz = buffer_1100_sddd[1055];

    auto g_y_y_0_0_0_zz_xz_xx = buffer_1100_sddd[1056];

    auto g_y_y_0_0_0_zz_xz_xy = buffer_1100_sddd[1057];

    auto g_y_y_0_0_0_zz_xz_xz = buffer_1100_sddd[1058];

    auto g_y_y_0_0_0_zz_xz_yy = buffer_1100_sddd[1059];

    auto g_y_y_0_0_0_zz_xz_yz = buffer_1100_sddd[1060];

    auto g_y_y_0_0_0_zz_xz_zz = buffer_1100_sddd[1061];

    auto g_y_y_0_0_0_zz_yy_xx = buffer_1100_sddd[1062];

    auto g_y_y_0_0_0_zz_yy_xy = buffer_1100_sddd[1063];

    auto g_y_y_0_0_0_zz_yy_xz = buffer_1100_sddd[1064];

    auto g_y_y_0_0_0_zz_yy_yy = buffer_1100_sddd[1065];

    auto g_y_y_0_0_0_zz_yy_yz = buffer_1100_sddd[1066];

    auto g_y_y_0_0_0_zz_yy_zz = buffer_1100_sddd[1067];

    auto g_y_y_0_0_0_zz_yz_xx = buffer_1100_sddd[1068];

    auto g_y_y_0_0_0_zz_yz_xy = buffer_1100_sddd[1069];

    auto g_y_y_0_0_0_zz_yz_xz = buffer_1100_sddd[1070];

    auto g_y_y_0_0_0_zz_yz_yy = buffer_1100_sddd[1071];

    auto g_y_y_0_0_0_zz_yz_yz = buffer_1100_sddd[1072];

    auto g_y_y_0_0_0_zz_yz_zz = buffer_1100_sddd[1073];

    auto g_y_y_0_0_0_zz_zz_xx = buffer_1100_sddd[1074];

    auto g_y_y_0_0_0_zz_zz_xy = buffer_1100_sddd[1075];

    auto g_y_y_0_0_0_zz_zz_xz = buffer_1100_sddd[1076];

    auto g_y_y_0_0_0_zz_zz_yy = buffer_1100_sddd[1077];

    auto g_y_y_0_0_0_zz_zz_yz = buffer_1100_sddd[1078];

    auto g_y_y_0_0_0_zz_zz_zz = buffer_1100_sddd[1079];

    auto g_y_z_0_0_0_xx_xx_xx = buffer_1100_sddd[1080];

    auto g_y_z_0_0_0_xx_xx_xy = buffer_1100_sddd[1081];

    auto g_y_z_0_0_0_xx_xx_xz = buffer_1100_sddd[1082];

    auto g_y_z_0_0_0_xx_xx_yy = buffer_1100_sddd[1083];

    auto g_y_z_0_0_0_xx_xx_yz = buffer_1100_sddd[1084];

    auto g_y_z_0_0_0_xx_xx_zz = buffer_1100_sddd[1085];

    auto g_y_z_0_0_0_xx_xy_xx = buffer_1100_sddd[1086];

    auto g_y_z_0_0_0_xx_xy_xy = buffer_1100_sddd[1087];

    auto g_y_z_0_0_0_xx_xy_xz = buffer_1100_sddd[1088];

    auto g_y_z_0_0_0_xx_xy_yy = buffer_1100_sddd[1089];

    auto g_y_z_0_0_0_xx_xy_yz = buffer_1100_sddd[1090];

    auto g_y_z_0_0_0_xx_xy_zz = buffer_1100_sddd[1091];

    auto g_y_z_0_0_0_xx_xz_xx = buffer_1100_sddd[1092];

    auto g_y_z_0_0_0_xx_xz_xy = buffer_1100_sddd[1093];

    auto g_y_z_0_0_0_xx_xz_xz = buffer_1100_sddd[1094];

    auto g_y_z_0_0_0_xx_xz_yy = buffer_1100_sddd[1095];

    auto g_y_z_0_0_0_xx_xz_yz = buffer_1100_sddd[1096];

    auto g_y_z_0_0_0_xx_xz_zz = buffer_1100_sddd[1097];

    auto g_y_z_0_0_0_xx_yy_xx = buffer_1100_sddd[1098];

    auto g_y_z_0_0_0_xx_yy_xy = buffer_1100_sddd[1099];

    auto g_y_z_0_0_0_xx_yy_xz = buffer_1100_sddd[1100];

    auto g_y_z_0_0_0_xx_yy_yy = buffer_1100_sddd[1101];

    auto g_y_z_0_0_0_xx_yy_yz = buffer_1100_sddd[1102];

    auto g_y_z_0_0_0_xx_yy_zz = buffer_1100_sddd[1103];

    auto g_y_z_0_0_0_xx_yz_xx = buffer_1100_sddd[1104];

    auto g_y_z_0_0_0_xx_yz_xy = buffer_1100_sddd[1105];

    auto g_y_z_0_0_0_xx_yz_xz = buffer_1100_sddd[1106];

    auto g_y_z_0_0_0_xx_yz_yy = buffer_1100_sddd[1107];

    auto g_y_z_0_0_0_xx_yz_yz = buffer_1100_sddd[1108];

    auto g_y_z_0_0_0_xx_yz_zz = buffer_1100_sddd[1109];

    auto g_y_z_0_0_0_xx_zz_xx = buffer_1100_sddd[1110];

    auto g_y_z_0_0_0_xx_zz_xy = buffer_1100_sddd[1111];

    auto g_y_z_0_0_0_xx_zz_xz = buffer_1100_sddd[1112];

    auto g_y_z_0_0_0_xx_zz_yy = buffer_1100_sddd[1113];

    auto g_y_z_0_0_0_xx_zz_yz = buffer_1100_sddd[1114];

    auto g_y_z_0_0_0_xx_zz_zz = buffer_1100_sddd[1115];

    auto g_y_z_0_0_0_xy_xx_xx = buffer_1100_sddd[1116];

    auto g_y_z_0_0_0_xy_xx_xy = buffer_1100_sddd[1117];

    auto g_y_z_0_0_0_xy_xx_xz = buffer_1100_sddd[1118];

    auto g_y_z_0_0_0_xy_xx_yy = buffer_1100_sddd[1119];

    auto g_y_z_0_0_0_xy_xx_yz = buffer_1100_sddd[1120];

    auto g_y_z_0_0_0_xy_xx_zz = buffer_1100_sddd[1121];

    auto g_y_z_0_0_0_xy_xy_xx = buffer_1100_sddd[1122];

    auto g_y_z_0_0_0_xy_xy_xy = buffer_1100_sddd[1123];

    auto g_y_z_0_0_0_xy_xy_xz = buffer_1100_sddd[1124];

    auto g_y_z_0_0_0_xy_xy_yy = buffer_1100_sddd[1125];

    auto g_y_z_0_0_0_xy_xy_yz = buffer_1100_sddd[1126];

    auto g_y_z_0_0_0_xy_xy_zz = buffer_1100_sddd[1127];

    auto g_y_z_0_0_0_xy_xz_xx = buffer_1100_sddd[1128];

    auto g_y_z_0_0_0_xy_xz_xy = buffer_1100_sddd[1129];

    auto g_y_z_0_0_0_xy_xz_xz = buffer_1100_sddd[1130];

    auto g_y_z_0_0_0_xy_xz_yy = buffer_1100_sddd[1131];

    auto g_y_z_0_0_0_xy_xz_yz = buffer_1100_sddd[1132];

    auto g_y_z_0_0_0_xy_xz_zz = buffer_1100_sddd[1133];

    auto g_y_z_0_0_0_xy_yy_xx = buffer_1100_sddd[1134];

    auto g_y_z_0_0_0_xy_yy_xy = buffer_1100_sddd[1135];

    auto g_y_z_0_0_0_xy_yy_xz = buffer_1100_sddd[1136];

    auto g_y_z_0_0_0_xy_yy_yy = buffer_1100_sddd[1137];

    auto g_y_z_0_0_0_xy_yy_yz = buffer_1100_sddd[1138];

    auto g_y_z_0_0_0_xy_yy_zz = buffer_1100_sddd[1139];

    auto g_y_z_0_0_0_xy_yz_xx = buffer_1100_sddd[1140];

    auto g_y_z_0_0_0_xy_yz_xy = buffer_1100_sddd[1141];

    auto g_y_z_0_0_0_xy_yz_xz = buffer_1100_sddd[1142];

    auto g_y_z_0_0_0_xy_yz_yy = buffer_1100_sddd[1143];

    auto g_y_z_0_0_0_xy_yz_yz = buffer_1100_sddd[1144];

    auto g_y_z_0_0_0_xy_yz_zz = buffer_1100_sddd[1145];

    auto g_y_z_0_0_0_xy_zz_xx = buffer_1100_sddd[1146];

    auto g_y_z_0_0_0_xy_zz_xy = buffer_1100_sddd[1147];

    auto g_y_z_0_0_0_xy_zz_xz = buffer_1100_sddd[1148];

    auto g_y_z_0_0_0_xy_zz_yy = buffer_1100_sddd[1149];

    auto g_y_z_0_0_0_xy_zz_yz = buffer_1100_sddd[1150];

    auto g_y_z_0_0_0_xy_zz_zz = buffer_1100_sddd[1151];

    auto g_y_z_0_0_0_xz_xx_xx = buffer_1100_sddd[1152];

    auto g_y_z_0_0_0_xz_xx_xy = buffer_1100_sddd[1153];

    auto g_y_z_0_0_0_xz_xx_xz = buffer_1100_sddd[1154];

    auto g_y_z_0_0_0_xz_xx_yy = buffer_1100_sddd[1155];

    auto g_y_z_0_0_0_xz_xx_yz = buffer_1100_sddd[1156];

    auto g_y_z_0_0_0_xz_xx_zz = buffer_1100_sddd[1157];

    auto g_y_z_0_0_0_xz_xy_xx = buffer_1100_sddd[1158];

    auto g_y_z_0_0_0_xz_xy_xy = buffer_1100_sddd[1159];

    auto g_y_z_0_0_0_xz_xy_xz = buffer_1100_sddd[1160];

    auto g_y_z_0_0_0_xz_xy_yy = buffer_1100_sddd[1161];

    auto g_y_z_0_0_0_xz_xy_yz = buffer_1100_sddd[1162];

    auto g_y_z_0_0_0_xz_xy_zz = buffer_1100_sddd[1163];

    auto g_y_z_0_0_0_xz_xz_xx = buffer_1100_sddd[1164];

    auto g_y_z_0_0_0_xz_xz_xy = buffer_1100_sddd[1165];

    auto g_y_z_0_0_0_xz_xz_xz = buffer_1100_sddd[1166];

    auto g_y_z_0_0_0_xz_xz_yy = buffer_1100_sddd[1167];

    auto g_y_z_0_0_0_xz_xz_yz = buffer_1100_sddd[1168];

    auto g_y_z_0_0_0_xz_xz_zz = buffer_1100_sddd[1169];

    auto g_y_z_0_0_0_xz_yy_xx = buffer_1100_sddd[1170];

    auto g_y_z_0_0_0_xz_yy_xy = buffer_1100_sddd[1171];

    auto g_y_z_0_0_0_xz_yy_xz = buffer_1100_sddd[1172];

    auto g_y_z_0_0_0_xz_yy_yy = buffer_1100_sddd[1173];

    auto g_y_z_0_0_0_xz_yy_yz = buffer_1100_sddd[1174];

    auto g_y_z_0_0_0_xz_yy_zz = buffer_1100_sddd[1175];

    auto g_y_z_0_0_0_xz_yz_xx = buffer_1100_sddd[1176];

    auto g_y_z_0_0_0_xz_yz_xy = buffer_1100_sddd[1177];

    auto g_y_z_0_0_0_xz_yz_xz = buffer_1100_sddd[1178];

    auto g_y_z_0_0_0_xz_yz_yy = buffer_1100_sddd[1179];

    auto g_y_z_0_0_0_xz_yz_yz = buffer_1100_sddd[1180];

    auto g_y_z_0_0_0_xz_yz_zz = buffer_1100_sddd[1181];

    auto g_y_z_0_0_0_xz_zz_xx = buffer_1100_sddd[1182];

    auto g_y_z_0_0_0_xz_zz_xy = buffer_1100_sddd[1183];

    auto g_y_z_0_0_0_xz_zz_xz = buffer_1100_sddd[1184];

    auto g_y_z_0_0_0_xz_zz_yy = buffer_1100_sddd[1185];

    auto g_y_z_0_0_0_xz_zz_yz = buffer_1100_sddd[1186];

    auto g_y_z_0_0_0_xz_zz_zz = buffer_1100_sddd[1187];

    auto g_y_z_0_0_0_yy_xx_xx = buffer_1100_sddd[1188];

    auto g_y_z_0_0_0_yy_xx_xy = buffer_1100_sddd[1189];

    auto g_y_z_0_0_0_yy_xx_xz = buffer_1100_sddd[1190];

    auto g_y_z_0_0_0_yy_xx_yy = buffer_1100_sddd[1191];

    auto g_y_z_0_0_0_yy_xx_yz = buffer_1100_sddd[1192];

    auto g_y_z_0_0_0_yy_xx_zz = buffer_1100_sddd[1193];

    auto g_y_z_0_0_0_yy_xy_xx = buffer_1100_sddd[1194];

    auto g_y_z_0_0_0_yy_xy_xy = buffer_1100_sddd[1195];

    auto g_y_z_0_0_0_yy_xy_xz = buffer_1100_sddd[1196];

    auto g_y_z_0_0_0_yy_xy_yy = buffer_1100_sddd[1197];

    auto g_y_z_0_0_0_yy_xy_yz = buffer_1100_sddd[1198];

    auto g_y_z_0_0_0_yy_xy_zz = buffer_1100_sddd[1199];

    auto g_y_z_0_0_0_yy_xz_xx = buffer_1100_sddd[1200];

    auto g_y_z_0_0_0_yy_xz_xy = buffer_1100_sddd[1201];

    auto g_y_z_0_0_0_yy_xz_xz = buffer_1100_sddd[1202];

    auto g_y_z_0_0_0_yy_xz_yy = buffer_1100_sddd[1203];

    auto g_y_z_0_0_0_yy_xz_yz = buffer_1100_sddd[1204];

    auto g_y_z_0_0_0_yy_xz_zz = buffer_1100_sddd[1205];

    auto g_y_z_0_0_0_yy_yy_xx = buffer_1100_sddd[1206];

    auto g_y_z_0_0_0_yy_yy_xy = buffer_1100_sddd[1207];

    auto g_y_z_0_0_0_yy_yy_xz = buffer_1100_sddd[1208];

    auto g_y_z_0_0_0_yy_yy_yy = buffer_1100_sddd[1209];

    auto g_y_z_0_0_0_yy_yy_yz = buffer_1100_sddd[1210];

    auto g_y_z_0_0_0_yy_yy_zz = buffer_1100_sddd[1211];

    auto g_y_z_0_0_0_yy_yz_xx = buffer_1100_sddd[1212];

    auto g_y_z_0_0_0_yy_yz_xy = buffer_1100_sddd[1213];

    auto g_y_z_0_0_0_yy_yz_xz = buffer_1100_sddd[1214];

    auto g_y_z_0_0_0_yy_yz_yy = buffer_1100_sddd[1215];

    auto g_y_z_0_0_0_yy_yz_yz = buffer_1100_sddd[1216];

    auto g_y_z_0_0_0_yy_yz_zz = buffer_1100_sddd[1217];

    auto g_y_z_0_0_0_yy_zz_xx = buffer_1100_sddd[1218];

    auto g_y_z_0_0_0_yy_zz_xy = buffer_1100_sddd[1219];

    auto g_y_z_0_0_0_yy_zz_xz = buffer_1100_sddd[1220];

    auto g_y_z_0_0_0_yy_zz_yy = buffer_1100_sddd[1221];

    auto g_y_z_0_0_0_yy_zz_yz = buffer_1100_sddd[1222];

    auto g_y_z_0_0_0_yy_zz_zz = buffer_1100_sddd[1223];

    auto g_y_z_0_0_0_yz_xx_xx = buffer_1100_sddd[1224];

    auto g_y_z_0_0_0_yz_xx_xy = buffer_1100_sddd[1225];

    auto g_y_z_0_0_0_yz_xx_xz = buffer_1100_sddd[1226];

    auto g_y_z_0_0_0_yz_xx_yy = buffer_1100_sddd[1227];

    auto g_y_z_0_0_0_yz_xx_yz = buffer_1100_sddd[1228];

    auto g_y_z_0_0_0_yz_xx_zz = buffer_1100_sddd[1229];

    auto g_y_z_0_0_0_yz_xy_xx = buffer_1100_sddd[1230];

    auto g_y_z_0_0_0_yz_xy_xy = buffer_1100_sddd[1231];

    auto g_y_z_0_0_0_yz_xy_xz = buffer_1100_sddd[1232];

    auto g_y_z_0_0_0_yz_xy_yy = buffer_1100_sddd[1233];

    auto g_y_z_0_0_0_yz_xy_yz = buffer_1100_sddd[1234];

    auto g_y_z_0_0_0_yz_xy_zz = buffer_1100_sddd[1235];

    auto g_y_z_0_0_0_yz_xz_xx = buffer_1100_sddd[1236];

    auto g_y_z_0_0_0_yz_xz_xy = buffer_1100_sddd[1237];

    auto g_y_z_0_0_0_yz_xz_xz = buffer_1100_sddd[1238];

    auto g_y_z_0_0_0_yz_xz_yy = buffer_1100_sddd[1239];

    auto g_y_z_0_0_0_yz_xz_yz = buffer_1100_sddd[1240];

    auto g_y_z_0_0_0_yz_xz_zz = buffer_1100_sddd[1241];

    auto g_y_z_0_0_0_yz_yy_xx = buffer_1100_sddd[1242];

    auto g_y_z_0_0_0_yz_yy_xy = buffer_1100_sddd[1243];

    auto g_y_z_0_0_0_yz_yy_xz = buffer_1100_sddd[1244];

    auto g_y_z_0_0_0_yz_yy_yy = buffer_1100_sddd[1245];

    auto g_y_z_0_0_0_yz_yy_yz = buffer_1100_sddd[1246];

    auto g_y_z_0_0_0_yz_yy_zz = buffer_1100_sddd[1247];

    auto g_y_z_0_0_0_yz_yz_xx = buffer_1100_sddd[1248];

    auto g_y_z_0_0_0_yz_yz_xy = buffer_1100_sddd[1249];

    auto g_y_z_0_0_0_yz_yz_xz = buffer_1100_sddd[1250];

    auto g_y_z_0_0_0_yz_yz_yy = buffer_1100_sddd[1251];

    auto g_y_z_0_0_0_yz_yz_yz = buffer_1100_sddd[1252];

    auto g_y_z_0_0_0_yz_yz_zz = buffer_1100_sddd[1253];

    auto g_y_z_0_0_0_yz_zz_xx = buffer_1100_sddd[1254];

    auto g_y_z_0_0_0_yz_zz_xy = buffer_1100_sddd[1255];

    auto g_y_z_0_0_0_yz_zz_xz = buffer_1100_sddd[1256];

    auto g_y_z_0_0_0_yz_zz_yy = buffer_1100_sddd[1257];

    auto g_y_z_0_0_0_yz_zz_yz = buffer_1100_sddd[1258];

    auto g_y_z_0_0_0_yz_zz_zz = buffer_1100_sddd[1259];

    auto g_y_z_0_0_0_zz_xx_xx = buffer_1100_sddd[1260];

    auto g_y_z_0_0_0_zz_xx_xy = buffer_1100_sddd[1261];

    auto g_y_z_0_0_0_zz_xx_xz = buffer_1100_sddd[1262];

    auto g_y_z_0_0_0_zz_xx_yy = buffer_1100_sddd[1263];

    auto g_y_z_0_0_0_zz_xx_yz = buffer_1100_sddd[1264];

    auto g_y_z_0_0_0_zz_xx_zz = buffer_1100_sddd[1265];

    auto g_y_z_0_0_0_zz_xy_xx = buffer_1100_sddd[1266];

    auto g_y_z_0_0_0_zz_xy_xy = buffer_1100_sddd[1267];

    auto g_y_z_0_0_0_zz_xy_xz = buffer_1100_sddd[1268];

    auto g_y_z_0_0_0_zz_xy_yy = buffer_1100_sddd[1269];

    auto g_y_z_0_0_0_zz_xy_yz = buffer_1100_sddd[1270];

    auto g_y_z_0_0_0_zz_xy_zz = buffer_1100_sddd[1271];

    auto g_y_z_0_0_0_zz_xz_xx = buffer_1100_sddd[1272];

    auto g_y_z_0_0_0_zz_xz_xy = buffer_1100_sddd[1273];

    auto g_y_z_0_0_0_zz_xz_xz = buffer_1100_sddd[1274];

    auto g_y_z_0_0_0_zz_xz_yy = buffer_1100_sddd[1275];

    auto g_y_z_0_0_0_zz_xz_yz = buffer_1100_sddd[1276];

    auto g_y_z_0_0_0_zz_xz_zz = buffer_1100_sddd[1277];

    auto g_y_z_0_0_0_zz_yy_xx = buffer_1100_sddd[1278];

    auto g_y_z_0_0_0_zz_yy_xy = buffer_1100_sddd[1279];

    auto g_y_z_0_0_0_zz_yy_xz = buffer_1100_sddd[1280];

    auto g_y_z_0_0_0_zz_yy_yy = buffer_1100_sddd[1281];

    auto g_y_z_0_0_0_zz_yy_yz = buffer_1100_sddd[1282];

    auto g_y_z_0_0_0_zz_yy_zz = buffer_1100_sddd[1283];

    auto g_y_z_0_0_0_zz_yz_xx = buffer_1100_sddd[1284];

    auto g_y_z_0_0_0_zz_yz_xy = buffer_1100_sddd[1285];

    auto g_y_z_0_0_0_zz_yz_xz = buffer_1100_sddd[1286];

    auto g_y_z_0_0_0_zz_yz_yy = buffer_1100_sddd[1287];

    auto g_y_z_0_0_0_zz_yz_yz = buffer_1100_sddd[1288];

    auto g_y_z_0_0_0_zz_yz_zz = buffer_1100_sddd[1289];

    auto g_y_z_0_0_0_zz_zz_xx = buffer_1100_sddd[1290];

    auto g_y_z_0_0_0_zz_zz_xy = buffer_1100_sddd[1291];

    auto g_y_z_0_0_0_zz_zz_xz = buffer_1100_sddd[1292];

    auto g_y_z_0_0_0_zz_zz_yy = buffer_1100_sddd[1293];

    auto g_y_z_0_0_0_zz_zz_yz = buffer_1100_sddd[1294];

    auto g_y_z_0_0_0_zz_zz_zz = buffer_1100_sddd[1295];

    auto g_z_x_0_0_0_xx_xx_xx = buffer_1100_sddd[1296];

    auto g_z_x_0_0_0_xx_xx_xy = buffer_1100_sddd[1297];

    auto g_z_x_0_0_0_xx_xx_xz = buffer_1100_sddd[1298];

    auto g_z_x_0_0_0_xx_xx_yy = buffer_1100_sddd[1299];

    auto g_z_x_0_0_0_xx_xx_yz = buffer_1100_sddd[1300];

    auto g_z_x_0_0_0_xx_xx_zz = buffer_1100_sddd[1301];

    auto g_z_x_0_0_0_xx_xy_xx = buffer_1100_sddd[1302];

    auto g_z_x_0_0_0_xx_xy_xy = buffer_1100_sddd[1303];

    auto g_z_x_0_0_0_xx_xy_xz = buffer_1100_sddd[1304];

    auto g_z_x_0_0_0_xx_xy_yy = buffer_1100_sddd[1305];

    auto g_z_x_0_0_0_xx_xy_yz = buffer_1100_sddd[1306];

    auto g_z_x_0_0_0_xx_xy_zz = buffer_1100_sddd[1307];

    auto g_z_x_0_0_0_xx_xz_xx = buffer_1100_sddd[1308];

    auto g_z_x_0_0_0_xx_xz_xy = buffer_1100_sddd[1309];

    auto g_z_x_0_0_0_xx_xz_xz = buffer_1100_sddd[1310];

    auto g_z_x_0_0_0_xx_xz_yy = buffer_1100_sddd[1311];

    auto g_z_x_0_0_0_xx_xz_yz = buffer_1100_sddd[1312];

    auto g_z_x_0_0_0_xx_xz_zz = buffer_1100_sddd[1313];

    auto g_z_x_0_0_0_xx_yy_xx = buffer_1100_sddd[1314];

    auto g_z_x_0_0_0_xx_yy_xy = buffer_1100_sddd[1315];

    auto g_z_x_0_0_0_xx_yy_xz = buffer_1100_sddd[1316];

    auto g_z_x_0_0_0_xx_yy_yy = buffer_1100_sddd[1317];

    auto g_z_x_0_0_0_xx_yy_yz = buffer_1100_sddd[1318];

    auto g_z_x_0_0_0_xx_yy_zz = buffer_1100_sddd[1319];

    auto g_z_x_0_0_0_xx_yz_xx = buffer_1100_sddd[1320];

    auto g_z_x_0_0_0_xx_yz_xy = buffer_1100_sddd[1321];

    auto g_z_x_0_0_0_xx_yz_xz = buffer_1100_sddd[1322];

    auto g_z_x_0_0_0_xx_yz_yy = buffer_1100_sddd[1323];

    auto g_z_x_0_0_0_xx_yz_yz = buffer_1100_sddd[1324];

    auto g_z_x_0_0_0_xx_yz_zz = buffer_1100_sddd[1325];

    auto g_z_x_0_0_0_xx_zz_xx = buffer_1100_sddd[1326];

    auto g_z_x_0_0_0_xx_zz_xy = buffer_1100_sddd[1327];

    auto g_z_x_0_0_0_xx_zz_xz = buffer_1100_sddd[1328];

    auto g_z_x_0_0_0_xx_zz_yy = buffer_1100_sddd[1329];

    auto g_z_x_0_0_0_xx_zz_yz = buffer_1100_sddd[1330];

    auto g_z_x_0_0_0_xx_zz_zz = buffer_1100_sddd[1331];

    auto g_z_x_0_0_0_xy_xx_xx = buffer_1100_sddd[1332];

    auto g_z_x_0_0_0_xy_xx_xy = buffer_1100_sddd[1333];

    auto g_z_x_0_0_0_xy_xx_xz = buffer_1100_sddd[1334];

    auto g_z_x_0_0_0_xy_xx_yy = buffer_1100_sddd[1335];

    auto g_z_x_0_0_0_xy_xx_yz = buffer_1100_sddd[1336];

    auto g_z_x_0_0_0_xy_xx_zz = buffer_1100_sddd[1337];

    auto g_z_x_0_0_0_xy_xy_xx = buffer_1100_sddd[1338];

    auto g_z_x_0_0_0_xy_xy_xy = buffer_1100_sddd[1339];

    auto g_z_x_0_0_0_xy_xy_xz = buffer_1100_sddd[1340];

    auto g_z_x_0_0_0_xy_xy_yy = buffer_1100_sddd[1341];

    auto g_z_x_0_0_0_xy_xy_yz = buffer_1100_sddd[1342];

    auto g_z_x_0_0_0_xy_xy_zz = buffer_1100_sddd[1343];

    auto g_z_x_0_0_0_xy_xz_xx = buffer_1100_sddd[1344];

    auto g_z_x_0_0_0_xy_xz_xy = buffer_1100_sddd[1345];

    auto g_z_x_0_0_0_xy_xz_xz = buffer_1100_sddd[1346];

    auto g_z_x_0_0_0_xy_xz_yy = buffer_1100_sddd[1347];

    auto g_z_x_0_0_0_xy_xz_yz = buffer_1100_sddd[1348];

    auto g_z_x_0_0_0_xy_xz_zz = buffer_1100_sddd[1349];

    auto g_z_x_0_0_0_xy_yy_xx = buffer_1100_sddd[1350];

    auto g_z_x_0_0_0_xy_yy_xy = buffer_1100_sddd[1351];

    auto g_z_x_0_0_0_xy_yy_xz = buffer_1100_sddd[1352];

    auto g_z_x_0_0_0_xy_yy_yy = buffer_1100_sddd[1353];

    auto g_z_x_0_0_0_xy_yy_yz = buffer_1100_sddd[1354];

    auto g_z_x_0_0_0_xy_yy_zz = buffer_1100_sddd[1355];

    auto g_z_x_0_0_0_xy_yz_xx = buffer_1100_sddd[1356];

    auto g_z_x_0_0_0_xy_yz_xy = buffer_1100_sddd[1357];

    auto g_z_x_0_0_0_xy_yz_xz = buffer_1100_sddd[1358];

    auto g_z_x_0_0_0_xy_yz_yy = buffer_1100_sddd[1359];

    auto g_z_x_0_0_0_xy_yz_yz = buffer_1100_sddd[1360];

    auto g_z_x_0_0_0_xy_yz_zz = buffer_1100_sddd[1361];

    auto g_z_x_0_0_0_xy_zz_xx = buffer_1100_sddd[1362];

    auto g_z_x_0_0_0_xy_zz_xy = buffer_1100_sddd[1363];

    auto g_z_x_0_0_0_xy_zz_xz = buffer_1100_sddd[1364];

    auto g_z_x_0_0_0_xy_zz_yy = buffer_1100_sddd[1365];

    auto g_z_x_0_0_0_xy_zz_yz = buffer_1100_sddd[1366];

    auto g_z_x_0_0_0_xy_zz_zz = buffer_1100_sddd[1367];

    auto g_z_x_0_0_0_xz_xx_xx = buffer_1100_sddd[1368];

    auto g_z_x_0_0_0_xz_xx_xy = buffer_1100_sddd[1369];

    auto g_z_x_0_0_0_xz_xx_xz = buffer_1100_sddd[1370];

    auto g_z_x_0_0_0_xz_xx_yy = buffer_1100_sddd[1371];

    auto g_z_x_0_0_0_xz_xx_yz = buffer_1100_sddd[1372];

    auto g_z_x_0_0_0_xz_xx_zz = buffer_1100_sddd[1373];

    auto g_z_x_0_0_0_xz_xy_xx = buffer_1100_sddd[1374];

    auto g_z_x_0_0_0_xz_xy_xy = buffer_1100_sddd[1375];

    auto g_z_x_0_0_0_xz_xy_xz = buffer_1100_sddd[1376];

    auto g_z_x_0_0_0_xz_xy_yy = buffer_1100_sddd[1377];

    auto g_z_x_0_0_0_xz_xy_yz = buffer_1100_sddd[1378];

    auto g_z_x_0_0_0_xz_xy_zz = buffer_1100_sddd[1379];

    auto g_z_x_0_0_0_xz_xz_xx = buffer_1100_sddd[1380];

    auto g_z_x_0_0_0_xz_xz_xy = buffer_1100_sddd[1381];

    auto g_z_x_0_0_0_xz_xz_xz = buffer_1100_sddd[1382];

    auto g_z_x_0_0_0_xz_xz_yy = buffer_1100_sddd[1383];

    auto g_z_x_0_0_0_xz_xz_yz = buffer_1100_sddd[1384];

    auto g_z_x_0_0_0_xz_xz_zz = buffer_1100_sddd[1385];

    auto g_z_x_0_0_0_xz_yy_xx = buffer_1100_sddd[1386];

    auto g_z_x_0_0_0_xz_yy_xy = buffer_1100_sddd[1387];

    auto g_z_x_0_0_0_xz_yy_xz = buffer_1100_sddd[1388];

    auto g_z_x_0_0_0_xz_yy_yy = buffer_1100_sddd[1389];

    auto g_z_x_0_0_0_xz_yy_yz = buffer_1100_sddd[1390];

    auto g_z_x_0_0_0_xz_yy_zz = buffer_1100_sddd[1391];

    auto g_z_x_0_0_0_xz_yz_xx = buffer_1100_sddd[1392];

    auto g_z_x_0_0_0_xz_yz_xy = buffer_1100_sddd[1393];

    auto g_z_x_0_0_0_xz_yz_xz = buffer_1100_sddd[1394];

    auto g_z_x_0_0_0_xz_yz_yy = buffer_1100_sddd[1395];

    auto g_z_x_0_0_0_xz_yz_yz = buffer_1100_sddd[1396];

    auto g_z_x_0_0_0_xz_yz_zz = buffer_1100_sddd[1397];

    auto g_z_x_0_0_0_xz_zz_xx = buffer_1100_sddd[1398];

    auto g_z_x_0_0_0_xz_zz_xy = buffer_1100_sddd[1399];

    auto g_z_x_0_0_0_xz_zz_xz = buffer_1100_sddd[1400];

    auto g_z_x_0_0_0_xz_zz_yy = buffer_1100_sddd[1401];

    auto g_z_x_0_0_0_xz_zz_yz = buffer_1100_sddd[1402];

    auto g_z_x_0_0_0_xz_zz_zz = buffer_1100_sddd[1403];

    auto g_z_x_0_0_0_yy_xx_xx = buffer_1100_sddd[1404];

    auto g_z_x_0_0_0_yy_xx_xy = buffer_1100_sddd[1405];

    auto g_z_x_0_0_0_yy_xx_xz = buffer_1100_sddd[1406];

    auto g_z_x_0_0_0_yy_xx_yy = buffer_1100_sddd[1407];

    auto g_z_x_0_0_0_yy_xx_yz = buffer_1100_sddd[1408];

    auto g_z_x_0_0_0_yy_xx_zz = buffer_1100_sddd[1409];

    auto g_z_x_0_0_0_yy_xy_xx = buffer_1100_sddd[1410];

    auto g_z_x_0_0_0_yy_xy_xy = buffer_1100_sddd[1411];

    auto g_z_x_0_0_0_yy_xy_xz = buffer_1100_sddd[1412];

    auto g_z_x_0_0_0_yy_xy_yy = buffer_1100_sddd[1413];

    auto g_z_x_0_0_0_yy_xy_yz = buffer_1100_sddd[1414];

    auto g_z_x_0_0_0_yy_xy_zz = buffer_1100_sddd[1415];

    auto g_z_x_0_0_0_yy_xz_xx = buffer_1100_sddd[1416];

    auto g_z_x_0_0_0_yy_xz_xy = buffer_1100_sddd[1417];

    auto g_z_x_0_0_0_yy_xz_xz = buffer_1100_sddd[1418];

    auto g_z_x_0_0_0_yy_xz_yy = buffer_1100_sddd[1419];

    auto g_z_x_0_0_0_yy_xz_yz = buffer_1100_sddd[1420];

    auto g_z_x_0_0_0_yy_xz_zz = buffer_1100_sddd[1421];

    auto g_z_x_0_0_0_yy_yy_xx = buffer_1100_sddd[1422];

    auto g_z_x_0_0_0_yy_yy_xy = buffer_1100_sddd[1423];

    auto g_z_x_0_0_0_yy_yy_xz = buffer_1100_sddd[1424];

    auto g_z_x_0_0_0_yy_yy_yy = buffer_1100_sddd[1425];

    auto g_z_x_0_0_0_yy_yy_yz = buffer_1100_sddd[1426];

    auto g_z_x_0_0_0_yy_yy_zz = buffer_1100_sddd[1427];

    auto g_z_x_0_0_0_yy_yz_xx = buffer_1100_sddd[1428];

    auto g_z_x_0_0_0_yy_yz_xy = buffer_1100_sddd[1429];

    auto g_z_x_0_0_0_yy_yz_xz = buffer_1100_sddd[1430];

    auto g_z_x_0_0_0_yy_yz_yy = buffer_1100_sddd[1431];

    auto g_z_x_0_0_0_yy_yz_yz = buffer_1100_sddd[1432];

    auto g_z_x_0_0_0_yy_yz_zz = buffer_1100_sddd[1433];

    auto g_z_x_0_0_0_yy_zz_xx = buffer_1100_sddd[1434];

    auto g_z_x_0_0_0_yy_zz_xy = buffer_1100_sddd[1435];

    auto g_z_x_0_0_0_yy_zz_xz = buffer_1100_sddd[1436];

    auto g_z_x_0_0_0_yy_zz_yy = buffer_1100_sddd[1437];

    auto g_z_x_0_0_0_yy_zz_yz = buffer_1100_sddd[1438];

    auto g_z_x_0_0_0_yy_zz_zz = buffer_1100_sddd[1439];

    auto g_z_x_0_0_0_yz_xx_xx = buffer_1100_sddd[1440];

    auto g_z_x_0_0_0_yz_xx_xy = buffer_1100_sddd[1441];

    auto g_z_x_0_0_0_yz_xx_xz = buffer_1100_sddd[1442];

    auto g_z_x_0_0_0_yz_xx_yy = buffer_1100_sddd[1443];

    auto g_z_x_0_0_0_yz_xx_yz = buffer_1100_sddd[1444];

    auto g_z_x_0_0_0_yz_xx_zz = buffer_1100_sddd[1445];

    auto g_z_x_0_0_0_yz_xy_xx = buffer_1100_sddd[1446];

    auto g_z_x_0_0_0_yz_xy_xy = buffer_1100_sddd[1447];

    auto g_z_x_0_0_0_yz_xy_xz = buffer_1100_sddd[1448];

    auto g_z_x_0_0_0_yz_xy_yy = buffer_1100_sddd[1449];

    auto g_z_x_0_0_0_yz_xy_yz = buffer_1100_sddd[1450];

    auto g_z_x_0_0_0_yz_xy_zz = buffer_1100_sddd[1451];

    auto g_z_x_0_0_0_yz_xz_xx = buffer_1100_sddd[1452];

    auto g_z_x_0_0_0_yz_xz_xy = buffer_1100_sddd[1453];

    auto g_z_x_0_0_0_yz_xz_xz = buffer_1100_sddd[1454];

    auto g_z_x_0_0_0_yz_xz_yy = buffer_1100_sddd[1455];

    auto g_z_x_0_0_0_yz_xz_yz = buffer_1100_sddd[1456];

    auto g_z_x_0_0_0_yz_xz_zz = buffer_1100_sddd[1457];

    auto g_z_x_0_0_0_yz_yy_xx = buffer_1100_sddd[1458];

    auto g_z_x_0_0_0_yz_yy_xy = buffer_1100_sddd[1459];

    auto g_z_x_0_0_0_yz_yy_xz = buffer_1100_sddd[1460];

    auto g_z_x_0_0_0_yz_yy_yy = buffer_1100_sddd[1461];

    auto g_z_x_0_0_0_yz_yy_yz = buffer_1100_sddd[1462];

    auto g_z_x_0_0_0_yz_yy_zz = buffer_1100_sddd[1463];

    auto g_z_x_0_0_0_yz_yz_xx = buffer_1100_sddd[1464];

    auto g_z_x_0_0_0_yz_yz_xy = buffer_1100_sddd[1465];

    auto g_z_x_0_0_0_yz_yz_xz = buffer_1100_sddd[1466];

    auto g_z_x_0_0_0_yz_yz_yy = buffer_1100_sddd[1467];

    auto g_z_x_0_0_0_yz_yz_yz = buffer_1100_sddd[1468];

    auto g_z_x_0_0_0_yz_yz_zz = buffer_1100_sddd[1469];

    auto g_z_x_0_0_0_yz_zz_xx = buffer_1100_sddd[1470];

    auto g_z_x_0_0_0_yz_zz_xy = buffer_1100_sddd[1471];

    auto g_z_x_0_0_0_yz_zz_xz = buffer_1100_sddd[1472];

    auto g_z_x_0_0_0_yz_zz_yy = buffer_1100_sddd[1473];

    auto g_z_x_0_0_0_yz_zz_yz = buffer_1100_sddd[1474];

    auto g_z_x_0_0_0_yz_zz_zz = buffer_1100_sddd[1475];

    auto g_z_x_0_0_0_zz_xx_xx = buffer_1100_sddd[1476];

    auto g_z_x_0_0_0_zz_xx_xy = buffer_1100_sddd[1477];

    auto g_z_x_0_0_0_zz_xx_xz = buffer_1100_sddd[1478];

    auto g_z_x_0_0_0_zz_xx_yy = buffer_1100_sddd[1479];

    auto g_z_x_0_0_0_zz_xx_yz = buffer_1100_sddd[1480];

    auto g_z_x_0_0_0_zz_xx_zz = buffer_1100_sddd[1481];

    auto g_z_x_0_0_0_zz_xy_xx = buffer_1100_sddd[1482];

    auto g_z_x_0_0_0_zz_xy_xy = buffer_1100_sddd[1483];

    auto g_z_x_0_0_0_zz_xy_xz = buffer_1100_sddd[1484];

    auto g_z_x_0_0_0_zz_xy_yy = buffer_1100_sddd[1485];

    auto g_z_x_0_0_0_zz_xy_yz = buffer_1100_sddd[1486];

    auto g_z_x_0_0_0_zz_xy_zz = buffer_1100_sddd[1487];

    auto g_z_x_0_0_0_zz_xz_xx = buffer_1100_sddd[1488];

    auto g_z_x_0_0_0_zz_xz_xy = buffer_1100_sddd[1489];

    auto g_z_x_0_0_0_zz_xz_xz = buffer_1100_sddd[1490];

    auto g_z_x_0_0_0_zz_xz_yy = buffer_1100_sddd[1491];

    auto g_z_x_0_0_0_zz_xz_yz = buffer_1100_sddd[1492];

    auto g_z_x_0_0_0_zz_xz_zz = buffer_1100_sddd[1493];

    auto g_z_x_0_0_0_zz_yy_xx = buffer_1100_sddd[1494];

    auto g_z_x_0_0_0_zz_yy_xy = buffer_1100_sddd[1495];

    auto g_z_x_0_0_0_zz_yy_xz = buffer_1100_sddd[1496];

    auto g_z_x_0_0_0_zz_yy_yy = buffer_1100_sddd[1497];

    auto g_z_x_0_0_0_zz_yy_yz = buffer_1100_sddd[1498];

    auto g_z_x_0_0_0_zz_yy_zz = buffer_1100_sddd[1499];

    auto g_z_x_0_0_0_zz_yz_xx = buffer_1100_sddd[1500];

    auto g_z_x_0_0_0_zz_yz_xy = buffer_1100_sddd[1501];

    auto g_z_x_0_0_0_zz_yz_xz = buffer_1100_sddd[1502];

    auto g_z_x_0_0_0_zz_yz_yy = buffer_1100_sddd[1503];

    auto g_z_x_0_0_0_zz_yz_yz = buffer_1100_sddd[1504];

    auto g_z_x_0_0_0_zz_yz_zz = buffer_1100_sddd[1505];

    auto g_z_x_0_0_0_zz_zz_xx = buffer_1100_sddd[1506];

    auto g_z_x_0_0_0_zz_zz_xy = buffer_1100_sddd[1507];

    auto g_z_x_0_0_0_zz_zz_xz = buffer_1100_sddd[1508];

    auto g_z_x_0_0_0_zz_zz_yy = buffer_1100_sddd[1509];

    auto g_z_x_0_0_0_zz_zz_yz = buffer_1100_sddd[1510];

    auto g_z_x_0_0_0_zz_zz_zz = buffer_1100_sddd[1511];

    auto g_z_y_0_0_0_xx_xx_xx = buffer_1100_sddd[1512];

    auto g_z_y_0_0_0_xx_xx_xy = buffer_1100_sddd[1513];

    auto g_z_y_0_0_0_xx_xx_xz = buffer_1100_sddd[1514];

    auto g_z_y_0_0_0_xx_xx_yy = buffer_1100_sddd[1515];

    auto g_z_y_0_0_0_xx_xx_yz = buffer_1100_sddd[1516];

    auto g_z_y_0_0_0_xx_xx_zz = buffer_1100_sddd[1517];

    auto g_z_y_0_0_0_xx_xy_xx = buffer_1100_sddd[1518];

    auto g_z_y_0_0_0_xx_xy_xy = buffer_1100_sddd[1519];

    auto g_z_y_0_0_0_xx_xy_xz = buffer_1100_sddd[1520];

    auto g_z_y_0_0_0_xx_xy_yy = buffer_1100_sddd[1521];

    auto g_z_y_0_0_0_xx_xy_yz = buffer_1100_sddd[1522];

    auto g_z_y_0_0_0_xx_xy_zz = buffer_1100_sddd[1523];

    auto g_z_y_0_0_0_xx_xz_xx = buffer_1100_sddd[1524];

    auto g_z_y_0_0_0_xx_xz_xy = buffer_1100_sddd[1525];

    auto g_z_y_0_0_0_xx_xz_xz = buffer_1100_sddd[1526];

    auto g_z_y_0_0_0_xx_xz_yy = buffer_1100_sddd[1527];

    auto g_z_y_0_0_0_xx_xz_yz = buffer_1100_sddd[1528];

    auto g_z_y_0_0_0_xx_xz_zz = buffer_1100_sddd[1529];

    auto g_z_y_0_0_0_xx_yy_xx = buffer_1100_sddd[1530];

    auto g_z_y_0_0_0_xx_yy_xy = buffer_1100_sddd[1531];

    auto g_z_y_0_0_0_xx_yy_xz = buffer_1100_sddd[1532];

    auto g_z_y_0_0_0_xx_yy_yy = buffer_1100_sddd[1533];

    auto g_z_y_0_0_0_xx_yy_yz = buffer_1100_sddd[1534];

    auto g_z_y_0_0_0_xx_yy_zz = buffer_1100_sddd[1535];

    auto g_z_y_0_0_0_xx_yz_xx = buffer_1100_sddd[1536];

    auto g_z_y_0_0_0_xx_yz_xy = buffer_1100_sddd[1537];

    auto g_z_y_0_0_0_xx_yz_xz = buffer_1100_sddd[1538];

    auto g_z_y_0_0_0_xx_yz_yy = buffer_1100_sddd[1539];

    auto g_z_y_0_0_0_xx_yz_yz = buffer_1100_sddd[1540];

    auto g_z_y_0_0_0_xx_yz_zz = buffer_1100_sddd[1541];

    auto g_z_y_0_0_0_xx_zz_xx = buffer_1100_sddd[1542];

    auto g_z_y_0_0_0_xx_zz_xy = buffer_1100_sddd[1543];

    auto g_z_y_0_0_0_xx_zz_xz = buffer_1100_sddd[1544];

    auto g_z_y_0_0_0_xx_zz_yy = buffer_1100_sddd[1545];

    auto g_z_y_0_0_0_xx_zz_yz = buffer_1100_sddd[1546];

    auto g_z_y_0_0_0_xx_zz_zz = buffer_1100_sddd[1547];

    auto g_z_y_0_0_0_xy_xx_xx = buffer_1100_sddd[1548];

    auto g_z_y_0_0_0_xy_xx_xy = buffer_1100_sddd[1549];

    auto g_z_y_0_0_0_xy_xx_xz = buffer_1100_sddd[1550];

    auto g_z_y_0_0_0_xy_xx_yy = buffer_1100_sddd[1551];

    auto g_z_y_0_0_0_xy_xx_yz = buffer_1100_sddd[1552];

    auto g_z_y_0_0_0_xy_xx_zz = buffer_1100_sddd[1553];

    auto g_z_y_0_0_0_xy_xy_xx = buffer_1100_sddd[1554];

    auto g_z_y_0_0_0_xy_xy_xy = buffer_1100_sddd[1555];

    auto g_z_y_0_0_0_xy_xy_xz = buffer_1100_sddd[1556];

    auto g_z_y_0_0_0_xy_xy_yy = buffer_1100_sddd[1557];

    auto g_z_y_0_0_0_xy_xy_yz = buffer_1100_sddd[1558];

    auto g_z_y_0_0_0_xy_xy_zz = buffer_1100_sddd[1559];

    auto g_z_y_0_0_0_xy_xz_xx = buffer_1100_sddd[1560];

    auto g_z_y_0_0_0_xy_xz_xy = buffer_1100_sddd[1561];

    auto g_z_y_0_0_0_xy_xz_xz = buffer_1100_sddd[1562];

    auto g_z_y_0_0_0_xy_xz_yy = buffer_1100_sddd[1563];

    auto g_z_y_0_0_0_xy_xz_yz = buffer_1100_sddd[1564];

    auto g_z_y_0_0_0_xy_xz_zz = buffer_1100_sddd[1565];

    auto g_z_y_0_0_0_xy_yy_xx = buffer_1100_sddd[1566];

    auto g_z_y_0_0_0_xy_yy_xy = buffer_1100_sddd[1567];

    auto g_z_y_0_0_0_xy_yy_xz = buffer_1100_sddd[1568];

    auto g_z_y_0_0_0_xy_yy_yy = buffer_1100_sddd[1569];

    auto g_z_y_0_0_0_xy_yy_yz = buffer_1100_sddd[1570];

    auto g_z_y_0_0_0_xy_yy_zz = buffer_1100_sddd[1571];

    auto g_z_y_0_0_0_xy_yz_xx = buffer_1100_sddd[1572];

    auto g_z_y_0_0_0_xy_yz_xy = buffer_1100_sddd[1573];

    auto g_z_y_0_0_0_xy_yz_xz = buffer_1100_sddd[1574];

    auto g_z_y_0_0_0_xy_yz_yy = buffer_1100_sddd[1575];

    auto g_z_y_0_0_0_xy_yz_yz = buffer_1100_sddd[1576];

    auto g_z_y_0_0_0_xy_yz_zz = buffer_1100_sddd[1577];

    auto g_z_y_0_0_0_xy_zz_xx = buffer_1100_sddd[1578];

    auto g_z_y_0_0_0_xy_zz_xy = buffer_1100_sddd[1579];

    auto g_z_y_0_0_0_xy_zz_xz = buffer_1100_sddd[1580];

    auto g_z_y_0_0_0_xy_zz_yy = buffer_1100_sddd[1581];

    auto g_z_y_0_0_0_xy_zz_yz = buffer_1100_sddd[1582];

    auto g_z_y_0_0_0_xy_zz_zz = buffer_1100_sddd[1583];

    auto g_z_y_0_0_0_xz_xx_xx = buffer_1100_sddd[1584];

    auto g_z_y_0_0_0_xz_xx_xy = buffer_1100_sddd[1585];

    auto g_z_y_0_0_0_xz_xx_xz = buffer_1100_sddd[1586];

    auto g_z_y_0_0_0_xz_xx_yy = buffer_1100_sddd[1587];

    auto g_z_y_0_0_0_xz_xx_yz = buffer_1100_sddd[1588];

    auto g_z_y_0_0_0_xz_xx_zz = buffer_1100_sddd[1589];

    auto g_z_y_0_0_0_xz_xy_xx = buffer_1100_sddd[1590];

    auto g_z_y_0_0_0_xz_xy_xy = buffer_1100_sddd[1591];

    auto g_z_y_0_0_0_xz_xy_xz = buffer_1100_sddd[1592];

    auto g_z_y_0_0_0_xz_xy_yy = buffer_1100_sddd[1593];

    auto g_z_y_0_0_0_xz_xy_yz = buffer_1100_sddd[1594];

    auto g_z_y_0_0_0_xz_xy_zz = buffer_1100_sddd[1595];

    auto g_z_y_0_0_0_xz_xz_xx = buffer_1100_sddd[1596];

    auto g_z_y_0_0_0_xz_xz_xy = buffer_1100_sddd[1597];

    auto g_z_y_0_0_0_xz_xz_xz = buffer_1100_sddd[1598];

    auto g_z_y_0_0_0_xz_xz_yy = buffer_1100_sddd[1599];

    auto g_z_y_0_0_0_xz_xz_yz = buffer_1100_sddd[1600];

    auto g_z_y_0_0_0_xz_xz_zz = buffer_1100_sddd[1601];

    auto g_z_y_0_0_0_xz_yy_xx = buffer_1100_sddd[1602];

    auto g_z_y_0_0_0_xz_yy_xy = buffer_1100_sddd[1603];

    auto g_z_y_0_0_0_xz_yy_xz = buffer_1100_sddd[1604];

    auto g_z_y_0_0_0_xz_yy_yy = buffer_1100_sddd[1605];

    auto g_z_y_0_0_0_xz_yy_yz = buffer_1100_sddd[1606];

    auto g_z_y_0_0_0_xz_yy_zz = buffer_1100_sddd[1607];

    auto g_z_y_0_0_0_xz_yz_xx = buffer_1100_sddd[1608];

    auto g_z_y_0_0_0_xz_yz_xy = buffer_1100_sddd[1609];

    auto g_z_y_0_0_0_xz_yz_xz = buffer_1100_sddd[1610];

    auto g_z_y_0_0_0_xz_yz_yy = buffer_1100_sddd[1611];

    auto g_z_y_0_0_0_xz_yz_yz = buffer_1100_sddd[1612];

    auto g_z_y_0_0_0_xz_yz_zz = buffer_1100_sddd[1613];

    auto g_z_y_0_0_0_xz_zz_xx = buffer_1100_sddd[1614];

    auto g_z_y_0_0_0_xz_zz_xy = buffer_1100_sddd[1615];

    auto g_z_y_0_0_0_xz_zz_xz = buffer_1100_sddd[1616];

    auto g_z_y_0_0_0_xz_zz_yy = buffer_1100_sddd[1617];

    auto g_z_y_0_0_0_xz_zz_yz = buffer_1100_sddd[1618];

    auto g_z_y_0_0_0_xz_zz_zz = buffer_1100_sddd[1619];

    auto g_z_y_0_0_0_yy_xx_xx = buffer_1100_sddd[1620];

    auto g_z_y_0_0_0_yy_xx_xy = buffer_1100_sddd[1621];

    auto g_z_y_0_0_0_yy_xx_xz = buffer_1100_sddd[1622];

    auto g_z_y_0_0_0_yy_xx_yy = buffer_1100_sddd[1623];

    auto g_z_y_0_0_0_yy_xx_yz = buffer_1100_sddd[1624];

    auto g_z_y_0_0_0_yy_xx_zz = buffer_1100_sddd[1625];

    auto g_z_y_0_0_0_yy_xy_xx = buffer_1100_sddd[1626];

    auto g_z_y_0_0_0_yy_xy_xy = buffer_1100_sddd[1627];

    auto g_z_y_0_0_0_yy_xy_xz = buffer_1100_sddd[1628];

    auto g_z_y_0_0_0_yy_xy_yy = buffer_1100_sddd[1629];

    auto g_z_y_0_0_0_yy_xy_yz = buffer_1100_sddd[1630];

    auto g_z_y_0_0_0_yy_xy_zz = buffer_1100_sddd[1631];

    auto g_z_y_0_0_0_yy_xz_xx = buffer_1100_sddd[1632];

    auto g_z_y_0_0_0_yy_xz_xy = buffer_1100_sddd[1633];

    auto g_z_y_0_0_0_yy_xz_xz = buffer_1100_sddd[1634];

    auto g_z_y_0_0_0_yy_xz_yy = buffer_1100_sddd[1635];

    auto g_z_y_0_0_0_yy_xz_yz = buffer_1100_sddd[1636];

    auto g_z_y_0_0_0_yy_xz_zz = buffer_1100_sddd[1637];

    auto g_z_y_0_0_0_yy_yy_xx = buffer_1100_sddd[1638];

    auto g_z_y_0_0_0_yy_yy_xy = buffer_1100_sddd[1639];

    auto g_z_y_0_0_0_yy_yy_xz = buffer_1100_sddd[1640];

    auto g_z_y_0_0_0_yy_yy_yy = buffer_1100_sddd[1641];

    auto g_z_y_0_0_0_yy_yy_yz = buffer_1100_sddd[1642];

    auto g_z_y_0_0_0_yy_yy_zz = buffer_1100_sddd[1643];

    auto g_z_y_0_0_0_yy_yz_xx = buffer_1100_sddd[1644];

    auto g_z_y_0_0_0_yy_yz_xy = buffer_1100_sddd[1645];

    auto g_z_y_0_0_0_yy_yz_xz = buffer_1100_sddd[1646];

    auto g_z_y_0_0_0_yy_yz_yy = buffer_1100_sddd[1647];

    auto g_z_y_0_0_0_yy_yz_yz = buffer_1100_sddd[1648];

    auto g_z_y_0_0_0_yy_yz_zz = buffer_1100_sddd[1649];

    auto g_z_y_0_0_0_yy_zz_xx = buffer_1100_sddd[1650];

    auto g_z_y_0_0_0_yy_zz_xy = buffer_1100_sddd[1651];

    auto g_z_y_0_0_0_yy_zz_xz = buffer_1100_sddd[1652];

    auto g_z_y_0_0_0_yy_zz_yy = buffer_1100_sddd[1653];

    auto g_z_y_0_0_0_yy_zz_yz = buffer_1100_sddd[1654];

    auto g_z_y_0_0_0_yy_zz_zz = buffer_1100_sddd[1655];

    auto g_z_y_0_0_0_yz_xx_xx = buffer_1100_sddd[1656];

    auto g_z_y_0_0_0_yz_xx_xy = buffer_1100_sddd[1657];

    auto g_z_y_0_0_0_yz_xx_xz = buffer_1100_sddd[1658];

    auto g_z_y_0_0_0_yz_xx_yy = buffer_1100_sddd[1659];

    auto g_z_y_0_0_0_yz_xx_yz = buffer_1100_sddd[1660];

    auto g_z_y_0_0_0_yz_xx_zz = buffer_1100_sddd[1661];

    auto g_z_y_0_0_0_yz_xy_xx = buffer_1100_sddd[1662];

    auto g_z_y_0_0_0_yz_xy_xy = buffer_1100_sddd[1663];

    auto g_z_y_0_0_0_yz_xy_xz = buffer_1100_sddd[1664];

    auto g_z_y_0_0_0_yz_xy_yy = buffer_1100_sddd[1665];

    auto g_z_y_0_0_0_yz_xy_yz = buffer_1100_sddd[1666];

    auto g_z_y_0_0_0_yz_xy_zz = buffer_1100_sddd[1667];

    auto g_z_y_0_0_0_yz_xz_xx = buffer_1100_sddd[1668];

    auto g_z_y_0_0_0_yz_xz_xy = buffer_1100_sddd[1669];

    auto g_z_y_0_0_0_yz_xz_xz = buffer_1100_sddd[1670];

    auto g_z_y_0_0_0_yz_xz_yy = buffer_1100_sddd[1671];

    auto g_z_y_0_0_0_yz_xz_yz = buffer_1100_sddd[1672];

    auto g_z_y_0_0_0_yz_xz_zz = buffer_1100_sddd[1673];

    auto g_z_y_0_0_0_yz_yy_xx = buffer_1100_sddd[1674];

    auto g_z_y_0_0_0_yz_yy_xy = buffer_1100_sddd[1675];

    auto g_z_y_0_0_0_yz_yy_xz = buffer_1100_sddd[1676];

    auto g_z_y_0_0_0_yz_yy_yy = buffer_1100_sddd[1677];

    auto g_z_y_0_0_0_yz_yy_yz = buffer_1100_sddd[1678];

    auto g_z_y_0_0_0_yz_yy_zz = buffer_1100_sddd[1679];

    auto g_z_y_0_0_0_yz_yz_xx = buffer_1100_sddd[1680];

    auto g_z_y_0_0_0_yz_yz_xy = buffer_1100_sddd[1681];

    auto g_z_y_0_0_0_yz_yz_xz = buffer_1100_sddd[1682];

    auto g_z_y_0_0_0_yz_yz_yy = buffer_1100_sddd[1683];

    auto g_z_y_0_0_0_yz_yz_yz = buffer_1100_sddd[1684];

    auto g_z_y_0_0_0_yz_yz_zz = buffer_1100_sddd[1685];

    auto g_z_y_0_0_0_yz_zz_xx = buffer_1100_sddd[1686];

    auto g_z_y_0_0_0_yz_zz_xy = buffer_1100_sddd[1687];

    auto g_z_y_0_0_0_yz_zz_xz = buffer_1100_sddd[1688];

    auto g_z_y_0_0_0_yz_zz_yy = buffer_1100_sddd[1689];

    auto g_z_y_0_0_0_yz_zz_yz = buffer_1100_sddd[1690];

    auto g_z_y_0_0_0_yz_zz_zz = buffer_1100_sddd[1691];

    auto g_z_y_0_0_0_zz_xx_xx = buffer_1100_sddd[1692];

    auto g_z_y_0_0_0_zz_xx_xy = buffer_1100_sddd[1693];

    auto g_z_y_0_0_0_zz_xx_xz = buffer_1100_sddd[1694];

    auto g_z_y_0_0_0_zz_xx_yy = buffer_1100_sddd[1695];

    auto g_z_y_0_0_0_zz_xx_yz = buffer_1100_sddd[1696];

    auto g_z_y_0_0_0_zz_xx_zz = buffer_1100_sddd[1697];

    auto g_z_y_0_0_0_zz_xy_xx = buffer_1100_sddd[1698];

    auto g_z_y_0_0_0_zz_xy_xy = buffer_1100_sddd[1699];

    auto g_z_y_0_0_0_zz_xy_xz = buffer_1100_sddd[1700];

    auto g_z_y_0_0_0_zz_xy_yy = buffer_1100_sddd[1701];

    auto g_z_y_0_0_0_zz_xy_yz = buffer_1100_sddd[1702];

    auto g_z_y_0_0_0_zz_xy_zz = buffer_1100_sddd[1703];

    auto g_z_y_0_0_0_zz_xz_xx = buffer_1100_sddd[1704];

    auto g_z_y_0_0_0_zz_xz_xy = buffer_1100_sddd[1705];

    auto g_z_y_0_0_0_zz_xz_xz = buffer_1100_sddd[1706];

    auto g_z_y_0_0_0_zz_xz_yy = buffer_1100_sddd[1707];

    auto g_z_y_0_0_0_zz_xz_yz = buffer_1100_sddd[1708];

    auto g_z_y_0_0_0_zz_xz_zz = buffer_1100_sddd[1709];

    auto g_z_y_0_0_0_zz_yy_xx = buffer_1100_sddd[1710];

    auto g_z_y_0_0_0_zz_yy_xy = buffer_1100_sddd[1711];

    auto g_z_y_0_0_0_zz_yy_xz = buffer_1100_sddd[1712];

    auto g_z_y_0_0_0_zz_yy_yy = buffer_1100_sddd[1713];

    auto g_z_y_0_0_0_zz_yy_yz = buffer_1100_sddd[1714];

    auto g_z_y_0_0_0_zz_yy_zz = buffer_1100_sddd[1715];

    auto g_z_y_0_0_0_zz_yz_xx = buffer_1100_sddd[1716];

    auto g_z_y_0_0_0_zz_yz_xy = buffer_1100_sddd[1717];

    auto g_z_y_0_0_0_zz_yz_xz = buffer_1100_sddd[1718];

    auto g_z_y_0_0_0_zz_yz_yy = buffer_1100_sddd[1719];

    auto g_z_y_0_0_0_zz_yz_yz = buffer_1100_sddd[1720];

    auto g_z_y_0_0_0_zz_yz_zz = buffer_1100_sddd[1721];

    auto g_z_y_0_0_0_zz_zz_xx = buffer_1100_sddd[1722];

    auto g_z_y_0_0_0_zz_zz_xy = buffer_1100_sddd[1723];

    auto g_z_y_0_0_0_zz_zz_xz = buffer_1100_sddd[1724];

    auto g_z_y_0_0_0_zz_zz_yy = buffer_1100_sddd[1725];

    auto g_z_y_0_0_0_zz_zz_yz = buffer_1100_sddd[1726];

    auto g_z_y_0_0_0_zz_zz_zz = buffer_1100_sddd[1727];

    auto g_z_z_0_0_0_xx_xx_xx = buffer_1100_sddd[1728];

    auto g_z_z_0_0_0_xx_xx_xy = buffer_1100_sddd[1729];

    auto g_z_z_0_0_0_xx_xx_xz = buffer_1100_sddd[1730];

    auto g_z_z_0_0_0_xx_xx_yy = buffer_1100_sddd[1731];

    auto g_z_z_0_0_0_xx_xx_yz = buffer_1100_sddd[1732];

    auto g_z_z_0_0_0_xx_xx_zz = buffer_1100_sddd[1733];

    auto g_z_z_0_0_0_xx_xy_xx = buffer_1100_sddd[1734];

    auto g_z_z_0_0_0_xx_xy_xy = buffer_1100_sddd[1735];

    auto g_z_z_0_0_0_xx_xy_xz = buffer_1100_sddd[1736];

    auto g_z_z_0_0_0_xx_xy_yy = buffer_1100_sddd[1737];

    auto g_z_z_0_0_0_xx_xy_yz = buffer_1100_sddd[1738];

    auto g_z_z_0_0_0_xx_xy_zz = buffer_1100_sddd[1739];

    auto g_z_z_0_0_0_xx_xz_xx = buffer_1100_sddd[1740];

    auto g_z_z_0_0_0_xx_xz_xy = buffer_1100_sddd[1741];

    auto g_z_z_0_0_0_xx_xz_xz = buffer_1100_sddd[1742];

    auto g_z_z_0_0_0_xx_xz_yy = buffer_1100_sddd[1743];

    auto g_z_z_0_0_0_xx_xz_yz = buffer_1100_sddd[1744];

    auto g_z_z_0_0_0_xx_xz_zz = buffer_1100_sddd[1745];

    auto g_z_z_0_0_0_xx_yy_xx = buffer_1100_sddd[1746];

    auto g_z_z_0_0_0_xx_yy_xy = buffer_1100_sddd[1747];

    auto g_z_z_0_0_0_xx_yy_xz = buffer_1100_sddd[1748];

    auto g_z_z_0_0_0_xx_yy_yy = buffer_1100_sddd[1749];

    auto g_z_z_0_0_0_xx_yy_yz = buffer_1100_sddd[1750];

    auto g_z_z_0_0_0_xx_yy_zz = buffer_1100_sddd[1751];

    auto g_z_z_0_0_0_xx_yz_xx = buffer_1100_sddd[1752];

    auto g_z_z_0_0_0_xx_yz_xy = buffer_1100_sddd[1753];

    auto g_z_z_0_0_0_xx_yz_xz = buffer_1100_sddd[1754];

    auto g_z_z_0_0_0_xx_yz_yy = buffer_1100_sddd[1755];

    auto g_z_z_0_0_0_xx_yz_yz = buffer_1100_sddd[1756];

    auto g_z_z_0_0_0_xx_yz_zz = buffer_1100_sddd[1757];

    auto g_z_z_0_0_0_xx_zz_xx = buffer_1100_sddd[1758];

    auto g_z_z_0_0_0_xx_zz_xy = buffer_1100_sddd[1759];

    auto g_z_z_0_0_0_xx_zz_xz = buffer_1100_sddd[1760];

    auto g_z_z_0_0_0_xx_zz_yy = buffer_1100_sddd[1761];

    auto g_z_z_0_0_0_xx_zz_yz = buffer_1100_sddd[1762];

    auto g_z_z_0_0_0_xx_zz_zz = buffer_1100_sddd[1763];

    auto g_z_z_0_0_0_xy_xx_xx = buffer_1100_sddd[1764];

    auto g_z_z_0_0_0_xy_xx_xy = buffer_1100_sddd[1765];

    auto g_z_z_0_0_0_xy_xx_xz = buffer_1100_sddd[1766];

    auto g_z_z_0_0_0_xy_xx_yy = buffer_1100_sddd[1767];

    auto g_z_z_0_0_0_xy_xx_yz = buffer_1100_sddd[1768];

    auto g_z_z_0_0_0_xy_xx_zz = buffer_1100_sddd[1769];

    auto g_z_z_0_0_0_xy_xy_xx = buffer_1100_sddd[1770];

    auto g_z_z_0_0_0_xy_xy_xy = buffer_1100_sddd[1771];

    auto g_z_z_0_0_0_xy_xy_xz = buffer_1100_sddd[1772];

    auto g_z_z_0_0_0_xy_xy_yy = buffer_1100_sddd[1773];

    auto g_z_z_0_0_0_xy_xy_yz = buffer_1100_sddd[1774];

    auto g_z_z_0_0_0_xy_xy_zz = buffer_1100_sddd[1775];

    auto g_z_z_0_0_0_xy_xz_xx = buffer_1100_sddd[1776];

    auto g_z_z_0_0_0_xy_xz_xy = buffer_1100_sddd[1777];

    auto g_z_z_0_0_0_xy_xz_xz = buffer_1100_sddd[1778];

    auto g_z_z_0_0_0_xy_xz_yy = buffer_1100_sddd[1779];

    auto g_z_z_0_0_0_xy_xz_yz = buffer_1100_sddd[1780];

    auto g_z_z_0_0_0_xy_xz_zz = buffer_1100_sddd[1781];

    auto g_z_z_0_0_0_xy_yy_xx = buffer_1100_sddd[1782];

    auto g_z_z_0_0_0_xy_yy_xy = buffer_1100_sddd[1783];

    auto g_z_z_0_0_0_xy_yy_xz = buffer_1100_sddd[1784];

    auto g_z_z_0_0_0_xy_yy_yy = buffer_1100_sddd[1785];

    auto g_z_z_0_0_0_xy_yy_yz = buffer_1100_sddd[1786];

    auto g_z_z_0_0_0_xy_yy_zz = buffer_1100_sddd[1787];

    auto g_z_z_0_0_0_xy_yz_xx = buffer_1100_sddd[1788];

    auto g_z_z_0_0_0_xy_yz_xy = buffer_1100_sddd[1789];

    auto g_z_z_0_0_0_xy_yz_xz = buffer_1100_sddd[1790];

    auto g_z_z_0_0_0_xy_yz_yy = buffer_1100_sddd[1791];

    auto g_z_z_0_0_0_xy_yz_yz = buffer_1100_sddd[1792];

    auto g_z_z_0_0_0_xy_yz_zz = buffer_1100_sddd[1793];

    auto g_z_z_0_0_0_xy_zz_xx = buffer_1100_sddd[1794];

    auto g_z_z_0_0_0_xy_zz_xy = buffer_1100_sddd[1795];

    auto g_z_z_0_0_0_xy_zz_xz = buffer_1100_sddd[1796];

    auto g_z_z_0_0_0_xy_zz_yy = buffer_1100_sddd[1797];

    auto g_z_z_0_0_0_xy_zz_yz = buffer_1100_sddd[1798];

    auto g_z_z_0_0_0_xy_zz_zz = buffer_1100_sddd[1799];

    auto g_z_z_0_0_0_xz_xx_xx = buffer_1100_sddd[1800];

    auto g_z_z_0_0_0_xz_xx_xy = buffer_1100_sddd[1801];

    auto g_z_z_0_0_0_xz_xx_xz = buffer_1100_sddd[1802];

    auto g_z_z_0_0_0_xz_xx_yy = buffer_1100_sddd[1803];

    auto g_z_z_0_0_0_xz_xx_yz = buffer_1100_sddd[1804];

    auto g_z_z_0_0_0_xz_xx_zz = buffer_1100_sddd[1805];

    auto g_z_z_0_0_0_xz_xy_xx = buffer_1100_sddd[1806];

    auto g_z_z_0_0_0_xz_xy_xy = buffer_1100_sddd[1807];

    auto g_z_z_0_0_0_xz_xy_xz = buffer_1100_sddd[1808];

    auto g_z_z_0_0_0_xz_xy_yy = buffer_1100_sddd[1809];

    auto g_z_z_0_0_0_xz_xy_yz = buffer_1100_sddd[1810];

    auto g_z_z_0_0_0_xz_xy_zz = buffer_1100_sddd[1811];

    auto g_z_z_0_0_0_xz_xz_xx = buffer_1100_sddd[1812];

    auto g_z_z_0_0_0_xz_xz_xy = buffer_1100_sddd[1813];

    auto g_z_z_0_0_0_xz_xz_xz = buffer_1100_sddd[1814];

    auto g_z_z_0_0_0_xz_xz_yy = buffer_1100_sddd[1815];

    auto g_z_z_0_0_0_xz_xz_yz = buffer_1100_sddd[1816];

    auto g_z_z_0_0_0_xz_xz_zz = buffer_1100_sddd[1817];

    auto g_z_z_0_0_0_xz_yy_xx = buffer_1100_sddd[1818];

    auto g_z_z_0_0_0_xz_yy_xy = buffer_1100_sddd[1819];

    auto g_z_z_0_0_0_xz_yy_xz = buffer_1100_sddd[1820];

    auto g_z_z_0_0_0_xz_yy_yy = buffer_1100_sddd[1821];

    auto g_z_z_0_0_0_xz_yy_yz = buffer_1100_sddd[1822];

    auto g_z_z_0_0_0_xz_yy_zz = buffer_1100_sddd[1823];

    auto g_z_z_0_0_0_xz_yz_xx = buffer_1100_sddd[1824];

    auto g_z_z_0_0_0_xz_yz_xy = buffer_1100_sddd[1825];

    auto g_z_z_0_0_0_xz_yz_xz = buffer_1100_sddd[1826];

    auto g_z_z_0_0_0_xz_yz_yy = buffer_1100_sddd[1827];

    auto g_z_z_0_0_0_xz_yz_yz = buffer_1100_sddd[1828];

    auto g_z_z_0_0_0_xz_yz_zz = buffer_1100_sddd[1829];

    auto g_z_z_0_0_0_xz_zz_xx = buffer_1100_sddd[1830];

    auto g_z_z_0_0_0_xz_zz_xy = buffer_1100_sddd[1831];

    auto g_z_z_0_0_0_xz_zz_xz = buffer_1100_sddd[1832];

    auto g_z_z_0_0_0_xz_zz_yy = buffer_1100_sddd[1833];

    auto g_z_z_0_0_0_xz_zz_yz = buffer_1100_sddd[1834];

    auto g_z_z_0_0_0_xz_zz_zz = buffer_1100_sddd[1835];

    auto g_z_z_0_0_0_yy_xx_xx = buffer_1100_sddd[1836];

    auto g_z_z_0_0_0_yy_xx_xy = buffer_1100_sddd[1837];

    auto g_z_z_0_0_0_yy_xx_xz = buffer_1100_sddd[1838];

    auto g_z_z_0_0_0_yy_xx_yy = buffer_1100_sddd[1839];

    auto g_z_z_0_0_0_yy_xx_yz = buffer_1100_sddd[1840];

    auto g_z_z_0_0_0_yy_xx_zz = buffer_1100_sddd[1841];

    auto g_z_z_0_0_0_yy_xy_xx = buffer_1100_sddd[1842];

    auto g_z_z_0_0_0_yy_xy_xy = buffer_1100_sddd[1843];

    auto g_z_z_0_0_0_yy_xy_xz = buffer_1100_sddd[1844];

    auto g_z_z_0_0_0_yy_xy_yy = buffer_1100_sddd[1845];

    auto g_z_z_0_0_0_yy_xy_yz = buffer_1100_sddd[1846];

    auto g_z_z_0_0_0_yy_xy_zz = buffer_1100_sddd[1847];

    auto g_z_z_0_0_0_yy_xz_xx = buffer_1100_sddd[1848];

    auto g_z_z_0_0_0_yy_xz_xy = buffer_1100_sddd[1849];

    auto g_z_z_0_0_0_yy_xz_xz = buffer_1100_sddd[1850];

    auto g_z_z_0_0_0_yy_xz_yy = buffer_1100_sddd[1851];

    auto g_z_z_0_0_0_yy_xz_yz = buffer_1100_sddd[1852];

    auto g_z_z_0_0_0_yy_xz_zz = buffer_1100_sddd[1853];

    auto g_z_z_0_0_0_yy_yy_xx = buffer_1100_sddd[1854];

    auto g_z_z_0_0_0_yy_yy_xy = buffer_1100_sddd[1855];

    auto g_z_z_0_0_0_yy_yy_xz = buffer_1100_sddd[1856];

    auto g_z_z_0_0_0_yy_yy_yy = buffer_1100_sddd[1857];

    auto g_z_z_0_0_0_yy_yy_yz = buffer_1100_sddd[1858];

    auto g_z_z_0_0_0_yy_yy_zz = buffer_1100_sddd[1859];

    auto g_z_z_0_0_0_yy_yz_xx = buffer_1100_sddd[1860];

    auto g_z_z_0_0_0_yy_yz_xy = buffer_1100_sddd[1861];

    auto g_z_z_0_0_0_yy_yz_xz = buffer_1100_sddd[1862];

    auto g_z_z_0_0_0_yy_yz_yy = buffer_1100_sddd[1863];

    auto g_z_z_0_0_0_yy_yz_yz = buffer_1100_sddd[1864];

    auto g_z_z_0_0_0_yy_yz_zz = buffer_1100_sddd[1865];

    auto g_z_z_0_0_0_yy_zz_xx = buffer_1100_sddd[1866];

    auto g_z_z_0_0_0_yy_zz_xy = buffer_1100_sddd[1867];

    auto g_z_z_0_0_0_yy_zz_xz = buffer_1100_sddd[1868];

    auto g_z_z_0_0_0_yy_zz_yy = buffer_1100_sddd[1869];

    auto g_z_z_0_0_0_yy_zz_yz = buffer_1100_sddd[1870];

    auto g_z_z_0_0_0_yy_zz_zz = buffer_1100_sddd[1871];

    auto g_z_z_0_0_0_yz_xx_xx = buffer_1100_sddd[1872];

    auto g_z_z_0_0_0_yz_xx_xy = buffer_1100_sddd[1873];

    auto g_z_z_0_0_0_yz_xx_xz = buffer_1100_sddd[1874];

    auto g_z_z_0_0_0_yz_xx_yy = buffer_1100_sddd[1875];

    auto g_z_z_0_0_0_yz_xx_yz = buffer_1100_sddd[1876];

    auto g_z_z_0_0_0_yz_xx_zz = buffer_1100_sddd[1877];

    auto g_z_z_0_0_0_yz_xy_xx = buffer_1100_sddd[1878];

    auto g_z_z_0_0_0_yz_xy_xy = buffer_1100_sddd[1879];

    auto g_z_z_0_0_0_yz_xy_xz = buffer_1100_sddd[1880];

    auto g_z_z_0_0_0_yz_xy_yy = buffer_1100_sddd[1881];

    auto g_z_z_0_0_0_yz_xy_yz = buffer_1100_sddd[1882];

    auto g_z_z_0_0_0_yz_xy_zz = buffer_1100_sddd[1883];

    auto g_z_z_0_0_0_yz_xz_xx = buffer_1100_sddd[1884];

    auto g_z_z_0_0_0_yz_xz_xy = buffer_1100_sddd[1885];

    auto g_z_z_0_0_0_yz_xz_xz = buffer_1100_sddd[1886];

    auto g_z_z_0_0_0_yz_xz_yy = buffer_1100_sddd[1887];

    auto g_z_z_0_0_0_yz_xz_yz = buffer_1100_sddd[1888];

    auto g_z_z_0_0_0_yz_xz_zz = buffer_1100_sddd[1889];

    auto g_z_z_0_0_0_yz_yy_xx = buffer_1100_sddd[1890];

    auto g_z_z_0_0_0_yz_yy_xy = buffer_1100_sddd[1891];

    auto g_z_z_0_0_0_yz_yy_xz = buffer_1100_sddd[1892];

    auto g_z_z_0_0_0_yz_yy_yy = buffer_1100_sddd[1893];

    auto g_z_z_0_0_0_yz_yy_yz = buffer_1100_sddd[1894];

    auto g_z_z_0_0_0_yz_yy_zz = buffer_1100_sddd[1895];

    auto g_z_z_0_0_0_yz_yz_xx = buffer_1100_sddd[1896];

    auto g_z_z_0_0_0_yz_yz_xy = buffer_1100_sddd[1897];

    auto g_z_z_0_0_0_yz_yz_xz = buffer_1100_sddd[1898];

    auto g_z_z_0_0_0_yz_yz_yy = buffer_1100_sddd[1899];

    auto g_z_z_0_0_0_yz_yz_yz = buffer_1100_sddd[1900];

    auto g_z_z_0_0_0_yz_yz_zz = buffer_1100_sddd[1901];

    auto g_z_z_0_0_0_yz_zz_xx = buffer_1100_sddd[1902];

    auto g_z_z_0_0_0_yz_zz_xy = buffer_1100_sddd[1903];

    auto g_z_z_0_0_0_yz_zz_xz = buffer_1100_sddd[1904];

    auto g_z_z_0_0_0_yz_zz_yy = buffer_1100_sddd[1905];

    auto g_z_z_0_0_0_yz_zz_yz = buffer_1100_sddd[1906];

    auto g_z_z_0_0_0_yz_zz_zz = buffer_1100_sddd[1907];

    auto g_z_z_0_0_0_zz_xx_xx = buffer_1100_sddd[1908];

    auto g_z_z_0_0_0_zz_xx_xy = buffer_1100_sddd[1909];

    auto g_z_z_0_0_0_zz_xx_xz = buffer_1100_sddd[1910];

    auto g_z_z_0_0_0_zz_xx_yy = buffer_1100_sddd[1911];

    auto g_z_z_0_0_0_zz_xx_yz = buffer_1100_sddd[1912];

    auto g_z_z_0_0_0_zz_xx_zz = buffer_1100_sddd[1913];

    auto g_z_z_0_0_0_zz_xy_xx = buffer_1100_sddd[1914];

    auto g_z_z_0_0_0_zz_xy_xy = buffer_1100_sddd[1915];

    auto g_z_z_0_0_0_zz_xy_xz = buffer_1100_sddd[1916];

    auto g_z_z_0_0_0_zz_xy_yy = buffer_1100_sddd[1917];

    auto g_z_z_0_0_0_zz_xy_yz = buffer_1100_sddd[1918];

    auto g_z_z_0_0_0_zz_xy_zz = buffer_1100_sddd[1919];

    auto g_z_z_0_0_0_zz_xz_xx = buffer_1100_sddd[1920];

    auto g_z_z_0_0_0_zz_xz_xy = buffer_1100_sddd[1921];

    auto g_z_z_0_0_0_zz_xz_xz = buffer_1100_sddd[1922];

    auto g_z_z_0_0_0_zz_xz_yy = buffer_1100_sddd[1923];

    auto g_z_z_0_0_0_zz_xz_yz = buffer_1100_sddd[1924];

    auto g_z_z_0_0_0_zz_xz_zz = buffer_1100_sddd[1925];

    auto g_z_z_0_0_0_zz_yy_xx = buffer_1100_sddd[1926];

    auto g_z_z_0_0_0_zz_yy_xy = buffer_1100_sddd[1927];

    auto g_z_z_0_0_0_zz_yy_xz = buffer_1100_sddd[1928];

    auto g_z_z_0_0_0_zz_yy_yy = buffer_1100_sddd[1929];

    auto g_z_z_0_0_0_zz_yy_yz = buffer_1100_sddd[1930];

    auto g_z_z_0_0_0_zz_yy_zz = buffer_1100_sddd[1931];

    auto g_z_z_0_0_0_zz_yz_xx = buffer_1100_sddd[1932];

    auto g_z_z_0_0_0_zz_yz_xy = buffer_1100_sddd[1933];

    auto g_z_z_0_0_0_zz_yz_xz = buffer_1100_sddd[1934];

    auto g_z_z_0_0_0_zz_yz_yy = buffer_1100_sddd[1935];

    auto g_z_z_0_0_0_zz_yz_yz = buffer_1100_sddd[1936];

    auto g_z_z_0_0_0_zz_yz_zz = buffer_1100_sddd[1937];

    auto g_z_z_0_0_0_zz_zz_xx = buffer_1100_sddd[1938];

    auto g_z_z_0_0_0_zz_zz_xy = buffer_1100_sddd[1939];

    auto g_z_z_0_0_0_zz_zz_xz = buffer_1100_sddd[1940];

    auto g_z_z_0_0_0_zz_zz_yy = buffer_1100_sddd[1941];

    auto g_z_z_0_0_0_zz_zz_yz = buffer_1100_sddd[1942];

    auto g_z_z_0_0_0_zz_zz_zz = buffer_1100_sddd[1943];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_x_0_0_0_xx_xx_xx, g_x_x_0_0_0_xx_xx_xy, g_x_x_0_0_0_xx_xx_xz, g_x_x_0_0_0_xx_xx_yy, g_x_x_0_0_0_xx_xx_yz, g_x_x_0_0_0_xx_xx_zz, g_x_x_xx_xx, g_x_x_xx_xy, g_x_x_xx_xz, g_x_x_xx_yy, g_x_x_xx_yz, g_x_x_xx_zz, g_x_xxx_xx_xx, g_x_xxx_xx_xy, g_x_xxx_xx_xz, g_x_xxx_xx_yy, g_x_xxx_xx_yz, g_x_xxx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xx_xx_xx[i] = -4.0 * g_x_x_xx_xx[i] * a_exp + 4.0 * g_x_xxx_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_xx_xy[i] = -4.0 * g_x_x_xx_xy[i] * a_exp + 4.0 * g_x_xxx_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_xx_xz[i] = -4.0 * g_x_x_xx_xz[i] * a_exp + 4.0 * g_x_xxx_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_xx_yy[i] = -4.0 * g_x_x_xx_yy[i] * a_exp + 4.0 * g_x_xxx_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_xx_yz[i] = -4.0 * g_x_x_xx_yz[i] * a_exp + 4.0 * g_x_xxx_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_xx_zz[i] = -4.0 * g_x_x_xx_zz[i] * a_exp + 4.0 * g_x_xxx_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_x_0_0_0_xx_xy_xx, g_x_x_0_0_0_xx_xy_xy, g_x_x_0_0_0_xx_xy_xz, g_x_x_0_0_0_xx_xy_yy, g_x_x_0_0_0_xx_xy_yz, g_x_x_0_0_0_xx_xy_zz, g_x_x_xy_xx, g_x_x_xy_xy, g_x_x_xy_xz, g_x_x_xy_yy, g_x_x_xy_yz, g_x_x_xy_zz, g_x_xxx_xy_xx, g_x_xxx_xy_xy, g_x_xxx_xy_xz, g_x_xxx_xy_yy, g_x_xxx_xy_yz, g_x_xxx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xx_xy_xx[i] = -4.0 * g_x_x_xy_xx[i] * a_exp + 4.0 * g_x_xxx_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_xy_xy[i] = -4.0 * g_x_x_xy_xy[i] * a_exp + 4.0 * g_x_xxx_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_xy_xz[i] = -4.0 * g_x_x_xy_xz[i] * a_exp + 4.0 * g_x_xxx_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_xy_yy[i] = -4.0 * g_x_x_xy_yy[i] * a_exp + 4.0 * g_x_xxx_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_xy_yz[i] = -4.0 * g_x_x_xy_yz[i] * a_exp + 4.0 * g_x_xxx_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_xy_zz[i] = -4.0 * g_x_x_xy_zz[i] * a_exp + 4.0 * g_x_xxx_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_x_0_0_0_xx_xz_xx, g_x_x_0_0_0_xx_xz_xy, g_x_x_0_0_0_xx_xz_xz, g_x_x_0_0_0_xx_xz_yy, g_x_x_0_0_0_xx_xz_yz, g_x_x_0_0_0_xx_xz_zz, g_x_x_xz_xx, g_x_x_xz_xy, g_x_x_xz_xz, g_x_x_xz_yy, g_x_x_xz_yz, g_x_x_xz_zz, g_x_xxx_xz_xx, g_x_xxx_xz_xy, g_x_xxx_xz_xz, g_x_xxx_xz_yy, g_x_xxx_xz_yz, g_x_xxx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xx_xz_xx[i] = -4.0 * g_x_x_xz_xx[i] * a_exp + 4.0 * g_x_xxx_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_xz_xy[i] = -4.0 * g_x_x_xz_xy[i] * a_exp + 4.0 * g_x_xxx_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_xz_xz[i] = -4.0 * g_x_x_xz_xz[i] * a_exp + 4.0 * g_x_xxx_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_xz_yy[i] = -4.0 * g_x_x_xz_yy[i] * a_exp + 4.0 * g_x_xxx_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_xz_yz[i] = -4.0 * g_x_x_xz_yz[i] * a_exp + 4.0 * g_x_xxx_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_xz_zz[i] = -4.0 * g_x_x_xz_zz[i] * a_exp + 4.0 * g_x_xxx_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_x_0_0_0_xx_yy_xx, g_x_x_0_0_0_xx_yy_xy, g_x_x_0_0_0_xx_yy_xz, g_x_x_0_0_0_xx_yy_yy, g_x_x_0_0_0_xx_yy_yz, g_x_x_0_0_0_xx_yy_zz, g_x_x_yy_xx, g_x_x_yy_xy, g_x_x_yy_xz, g_x_x_yy_yy, g_x_x_yy_yz, g_x_x_yy_zz, g_x_xxx_yy_xx, g_x_xxx_yy_xy, g_x_xxx_yy_xz, g_x_xxx_yy_yy, g_x_xxx_yy_yz, g_x_xxx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xx_yy_xx[i] = -4.0 * g_x_x_yy_xx[i] * a_exp + 4.0 * g_x_xxx_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_yy_xy[i] = -4.0 * g_x_x_yy_xy[i] * a_exp + 4.0 * g_x_xxx_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_yy_xz[i] = -4.0 * g_x_x_yy_xz[i] * a_exp + 4.0 * g_x_xxx_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_yy_yy[i] = -4.0 * g_x_x_yy_yy[i] * a_exp + 4.0 * g_x_xxx_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_yy_yz[i] = -4.0 * g_x_x_yy_yz[i] * a_exp + 4.0 * g_x_xxx_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_yy_zz[i] = -4.0 * g_x_x_yy_zz[i] * a_exp + 4.0 * g_x_xxx_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_x_0_0_0_xx_yz_xx, g_x_x_0_0_0_xx_yz_xy, g_x_x_0_0_0_xx_yz_xz, g_x_x_0_0_0_xx_yz_yy, g_x_x_0_0_0_xx_yz_yz, g_x_x_0_0_0_xx_yz_zz, g_x_x_yz_xx, g_x_x_yz_xy, g_x_x_yz_xz, g_x_x_yz_yy, g_x_x_yz_yz, g_x_x_yz_zz, g_x_xxx_yz_xx, g_x_xxx_yz_xy, g_x_xxx_yz_xz, g_x_xxx_yz_yy, g_x_xxx_yz_yz, g_x_xxx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xx_yz_xx[i] = -4.0 * g_x_x_yz_xx[i] * a_exp + 4.0 * g_x_xxx_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_yz_xy[i] = -4.0 * g_x_x_yz_xy[i] * a_exp + 4.0 * g_x_xxx_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_yz_xz[i] = -4.0 * g_x_x_yz_xz[i] * a_exp + 4.0 * g_x_xxx_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_yz_yy[i] = -4.0 * g_x_x_yz_yy[i] * a_exp + 4.0 * g_x_xxx_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_yz_yz[i] = -4.0 * g_x_x_yz_yz[i] * a_exp + 4.0 * g_x_xxx_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_yz_zz[i] = -4.0 * g_x_x_yz_zz[i] * a_exp + 4.0 * g_x_xxx_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_x_0_0_0_xx_zz_xx, g_x_x_0_0_0_xx_zz_xy, g_x_x_0_0_0_xx_zz_xz, g_x_x_0_0_0_xx_zz_yy, g_x_x_0_0_0_xx_zz_yz, g_x_x_0_0_0_xx_zz_zz, g_x_x_zz_xx, g_x_x_zz_xy, g_x_x_zz_xz, g_x_x_zz_yy, g_x_x_zz_yz, g_x_x_zz_zz, g_x_xxx_zz_xx, g_x_xxx_zz_xy, g_x_xxx_zz_xz, g_x_xxx_zz_yy, g_x_xxx_zz_yz, g_x_xxx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xx_zz_xx[i] = -4.0 * g_x_x_zz_xx[i] * a_exp + 4.0 * g_x_xxx_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_zz_xy[i] = -4.0 * g_x_x_zz_xy[i] * a_exp + 4.0 * g_x_xxx_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_zz_xz[i] = -4.0 * g_x_x_zz_xz[i] * a_exp + 4.0 * g_x_xxx_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_zz_yy[i] = -4.0 * g_x_x_zz_yy[i] * a_exp + 4.0 * g_x_xxx_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_zz_yz[i] = -4.0 * g_x_x_zz_yz[i] * a_exp + 4.0 * g_x_xxx_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_zz_zz[i] = -4.0 * g_x_x_zz_zz[i] * a_exp + 4.0 * g_x_xxx_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_x_0_0_0_xy_xx_xx, g_x_x_0_0_0_xy_xx_xy, g_x_x_0_0_0_xy_xx_xz, g_x_x_0_0_0_xy_xx_yy, g_x_x_0_0_0_xy_xx_yz, g_x_x_0_0_0_xy_xx_zz, g_x_xxy_xx_xx, g_x_xxy_xx_xy, g_x_xxy_xx_xz, g_x_xxy_xx_yy, g_x_xxy_xx_yz, g_x_xxy_xx_zz, g_x_y_xx_xx, g_x_y_xx_xy, g_x_y_xx_xz, g_x_y_xx_yy, g_x_y_xx_yz, g_x_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xy_xx_xx[i] = -2.0 * g_x_y_xx_xx[i] * a_exp + 4.0 * g_x_xxy_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_xx_xy[i] = -2.0 * g_x_y_xx_xy[i] * a_exp + 4.0 * g_x_xxy_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_xx_xz[i] = -2.0 * g_x_y_xx_xz[i] * a_exp + 4.0 * g_x_xxy_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_xx_yy[i] = -2.0 * g_x_y_xx_yy[i] * a_exp + 4.0 * g_x_xxy_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_xx_yz[i] = -2.0 * g_x_y_xx_yz[i] * a_exp + 4.0 * g_x_xxy_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_xx_zz[i] = -2.0 * g_x_y_xx_zz[i] * a_exp + 4.0 * g_x_xxy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_x_0_0_0_xy_xy_xx, g_x_x_0_0_0_xy_xy_xy, g_x_x_0_0_0_xy_xy_xz, g_x_x_0_0_0_xy_xy_yy, g_x_x_0_0_0_xy_xy_yz, g_x_x_0_0_0_xy_xy_zz, g_x_xxy_xy_xx, g_x_xxy_xy_xy, g_x_xxy_xy_xz, g_x_xxy_xy_yy, g_x_xxy_xy_yz, g_x_xxy_xy_zz, g_x_y_xy_xx, g_x_y_xy_xy, g_x_y_xy_xz, g_x_y_xy_yy, g_x_y_xy_yz, g_x_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xy_xy_xx[i] = -2.0 * g_x_y_xy_xx[i] * a_exp + 4.0 * g_x_xxy_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_xy_xy[i] = -2.0 * g_x_y_xy_xy[i] * a_exp + 4.0 * g_x_xxy_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_xy_xz[i] = -2.0 * g_x_y_xy_xz[i] * a_exp + 4.0 * g_x_xxy_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_xy_yy[i] = -2.0 * g_x_y_xy_yy[i] * a_exp + 4.0 * g_x_xxy_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_xy_yz[i] = -2.0 * g_x_y_xy_yz[i] * a_exp + 4.0 * g_x_xxy_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_xy_zz[i] = -2.0 * g_x_y_xy_zz[i] * a_exp + 4.0 * g_x_xxy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_x_0_0_0_xy_xz_xx, g_x_x_0_0_0_xy_xz_xy, g_x_x_0_0_0_xy_xz_xz, g_x_x_0_0_0_xy_xz_yy, g_x_x_0_0_0_xy_xz_yz, g_x_x_0_0_0_xy_xz_zz, g_x_xxy_xz_xx, g_x_xxy_xz_xy, g_x_xxy_xz_xz, g_x_xxy_xz_yy, g_x_xxy_xz_yz, g_x_xxy_xz_zz, g_x_y_xz_xx, g_x_y_xz_xy, g_x_y_xz_xz, g_x_y_xz_yy, g_x_y_xz_yz, g_x_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xy_xz_xx[i] = -2.0 * g_x_y_xz_xx[i] * a_exp + 4.0 * g_x_xxy_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_xz_xy[i] = -2.0 * g_x_y_xz_xy[i] * a_exp + 4.0 * g_x_xxy_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_xz_xz[i] = -2.0 * g_x_y_xz_xz[i] * a_exp + 4.0 * g_x_xxy_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_xz_yy[i] = -2.0 * g_x_y_xz_yy[i] * a_exp + 4.0 * g_x_xxy_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_xz_yz[i] = -2.0 * g_x_y_xz_yz[i] * a_exp + 4.0 * g_x_xxy_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_xz_zz[i] = -2.0 * g_x_y_xz_zz[i] * a_exp + 4.0 * g_x_xxy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_x_0_0_0_xy_yy_xx, g_x_x_0_0_0_xy_yy_xy, g_x_x_0_0_0_xy_yy_xz, g_x_x_0_0_0_xy_yy_yy, g_x_x_0_0_0_xy_yy_yz, g_x_x_0_0_0_xy_yy_zz, g_x_xxy_yy_xx, g_x_xxy_yy_xy, g_x_xxy_yy_xz, g_x_xxy_yy_yy, g_x_xxy_yy_yz, g_x_xxy_yy_zz, g_x_y_yy_xx, g_x_y_yy_xy, g_x_y_yy_xz, g_x_y_yy_yy, g_x_y_yy_yz, g_x_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xy_yy_xx[i] = -2.0 * g_x_y_yy_xx[i] * a_exp + 4.0 * g_x_xxy_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_yy_xy[i] = -2.0 * g_x_y_yy_xy[i] * a_exp + 4.0 * g_x_xxy_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_yy_xz[i] = -2.0 * g_x_y_yy_xz[i] * a_exp + 4.0 * g_x_xxy_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_yy_yy[i] = -2.0 * g_x_y_yy_yy[i] * a_exp + 4.0 * g_x_xxy_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_yy_yz[i] = -2.0 * g_x_y_yy_yz[i] * a_exp + 4.0 * g_x_xxy_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_yy_zz[i] = -2.0 * g_x_y_yy_zz[i] * a_exp + 4.0 * g_x_xxy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_x_0_0_0_xy_yz_xx, g_x_x_0_0_0_xy_yz_xy, g_x_x_0_0_0_xy_yz_xz, g_x_x_0_0_0_xy_yz_yy, g_x_x_0_0_0_xy_yz_yz, g_x_x_0_0_0_xy_yz_zz, g_x_xxy_yz_xx, g_x_xxy_yz_xy, g_x_xxy_yz_xz, g_x_xxy_yz_yy, g_x_xxy_yz_yz, g_x_xxy_yz_zz, g_x_y_yz_xx, g_x_y_yz_xy, g_x_y_yz_xz, g_x_y_yz_yy, g_x_y_yz_yz, g_x_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xy_yz_xx[i] = -2.0 * g_x_y_yz_xx[i] * a_exp + 4.0 * g_x_xxy_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_yz_xy[i] = -2.0 * g_x_y_yz_xy[i] * a_exp + 4.0 * g_x_xxy_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_yz_xz[i] = -2.0 * g_x_y_yz_xz[i] * a_exp + 4.0 * g_x_xxy_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_yz_yy[i] = -2.0 * g_x_y_yz_yy[i] * a_exp + 4.0 * g_x_xxy_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_yz_yz[i] = -2.0 * g_x_y_yz_yz[i] * a_exp + 4.0 * g_x_xxy_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_yz_zz[i] = -2.0 * g_x_y_yz_zz[i] * a_exp + 4.0 * g_x_xxy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_x_0_0_0_xy_zz_xx, g_x_x_0_0_0_xy_zz_xy, g_x_x_0_0_0_xy_zz_xz, g_x_x_0_0_0_xy_zz_yy, g_x_x_0_0_0_xy_zz_yz, g_x_x_0_0_0_xy_zz_zz, g_x_xxy_zz_xx, g_x_xxy_zz_xy, g_x_xxy_zz_xz, g_x_xxy_zz_yy, g_x_xxy_zz_yz, g_x_xxy_zz_zz, g_x_y_zz_xx, g_x_y_zz_xy, g_x_y_zz_xz, g_x_y_zz_yy, g_x_y_zz_yz, g_x_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xy_zz_xx[i] = -2.0 * g_x_y_zz_xx[i] * a_exp + 4.0 * g_x_xxy_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_zz_xy[i] = -2.0 * g_x_y_zz_xy[i] * a_exp + 4.0 * g_x_xxy_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_zz_xz[i] = -2.0 * g_x_y_zz_xz[i] * a_exp + 4.0 * g_x_xxy_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_zz_yy[i] = -2.0 * g_x_y_zz_yy[i] * a_exp + 4.0 * g_x_xxy_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_zz_yz[i] = -2.0 * g_x_y_zz_yz[i] * a_exp + 4.0 * g_x_xxy_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_zz_zz[i] = -2.0 * g_x_y_zz_zz[i] * a_exp + 4.0 * g_x_xxy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_x_0_0_0_xz_xx_xx, g_x_x_0_0_0_xz_xx_xy, g_x_x_0_0_0_xz_xx_xz, g_x_x_0_0_0_xz_xx_yy, g_x_x_0_0_0_xz_xx_yz, g_x_x_0_0_0_xz_xx_zz, g_x_xxz_xx_xx, g_x_xxz_xx_xy, g_x_xxz_xx_xz, g_x_xxz_xx_yy, g_x_xxz_xx_yz, g_x_xxz_xx_zz, g_x_z_xx_xx, g_x_z_xx_xy, g_x_z_xx_xz, g_x_z_xx_yy, g_x_z_xx_yz, g_x_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xz_xx_xx[i] = -2.0 * g_x_z_xx_xx[i] * a_exp + 4.0 * g_x_xxz_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_xx_xy[i] = -2.0 * g_x_z_xx_xy[i] * a_exp + 4.0 * g_x_xxz_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_xx_xz[i] = -2.0 * g_x_z_xx_xz[i] * a_exp + 4.0 * g_x_xxz_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_xx_yy[i] = -2.0 * g_x_z_xx_yy[i] * a_exp + 4.0 * g_x_xxz_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_xx_yz[i] = -2.0 * g_x_z_xx_yz[i] * a_exp + 4.0 * g_x_xxz_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_xx_zz[i] = -2.0 * g_x_z_xx_zz[i] * a_exp + 4.0 * g_x_xxz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_x_0_0_0_xz_xy_xx, g_x_x_0_0_0_xz_xy_xy, g_x_x_0_0_0_xz_xy_xz, g_x_x_0_0_0_xz_xy_yy, g_x_x_0_0_0_xz_xy_yz, g_x_x_0_0_0_xz_xy_zz, g_x_xxz_xy_xx, g_x_xxz_xy_xy, g_x_xxz_xy_xz, g_x_xxz_xy_yy, g_x_xxz_xy_yz, g_x_xxz_xy_zz, g_x_z_xy_xx, g_x_z_xy_xy, g_x_z_xy_xz, g_x_z_xy_yy, g_x_z_xy_yz, g_x_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xz_xy_xx[i] = -2.0 * g_x_z_xy_xx[i] * a_exp + 4.0 * g_x_xxz_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_xy_xy[i] = -2.0 * g_x_z_xy_xy[i] * a_exp + 4.0 * g_x_xxz_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_xy_xz[i] = -2.0 * g_x_z_xy_xz[i] * a_exp + 4.0 * g_x_xxz_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_xy_yy[i] = -2.0 * g_x_z_xy_yy[i] * a_exp + 4.0 * g_x_xxz_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_xy_yz[i] = -2.0 * g_x_z_xy_yz[i] * a_exp + 4.0 * g_x_xxz_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_xy_zz[i] = -2.0 * g_x_z_xy_zz[i] * a_exp + 4.0 * g_x_xxz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_x_0_0_0_xz_xz_xx, g_x_x_0_0_0_xz_xz_xy, g_x_x_0_0_0_xz_xz_xz, g_x_x_0_0_0_xz_xz_yy, g_x_x_0_0_0_xz_xz_yz, g_x_x_0_0_0_xz_xz_zz, g_x_xxz_xz_xx, g_x_xxz_xz_xy, g_x_xxz_xz_xz, g_x_xxz_xz_yy, g_x_xxz_xz_yz, g_x_xxz_xz_zz, g_x_z_xz_xx, g_x_z_xz_xy, g_x_z_xz_xz, g_x_z_xz_yy, g_x_z_xz_yz, g_x_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xz_xz_xx[i] = -2.0 * g_x_z_xz_xx[i] * a_exp + 4.0 * g_x_xxz_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_xz_xy[i] = -2.0 * g_x_z_xz_xy[i] * a_exp + 4.0 * g_x_xxz_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_xz_xz[i] = -2.0 * g_x_z_xz_xz[i] * a_exp + 4.0 * g_x_xxz_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_xz_yy[i] = -2.0 * g_x_z_xz_yy[i] * a_exp + 4.0 * g_x_xxz_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_xz_yz[i] = -2.0 * g_x_z_xz_yz[i] * a_exp + 4.0 * g_x_xxz_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_xz_zz[i] = -2.0 * g_x_z_xz_zz[i] * a_exp + 4.0 * g_x_xxz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_x_0_0_0_xz_yy_xx, g_x_x_0_0_0_xz_yy_xy, g_x_x_0_0_0_xz_yy_xz, g_x_x_0_0_0_xz_yy_yy, g_x_x_0_0_0_xz_yy_yz, g_x_x_0_0_0_xz_yy_zz, g_x_xxz_yy_xx, g_x_xxz_yy_xy, g_x_xxz_yy_xz, g_x_xxz_yy_yy, g_x_xxz_yy_yz, g_x_xxz_yy_zz, g_x_z_yy_xx, g_x_z_yy_xy, g_x_z_yy_xz, g_x_z_yy_yy, g_x_z_yy_yz, g_x_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xz_yy_xx[i] = -2.0 * g_x_z_yy_xx[i] * a_exp + 4.0 * g_x_xxz_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_yy_xy[i] = -2.0 * g_x_z_yy_xy[i] * a_exp + 4.0 * g_x_xxz_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_yy_xz[i] = -2.0 * g_x_z_yy_xz[i] * a_exp + 4.0 * g_x_xxz_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_yy_yy[i] = -2.0 * g_x_z_yy_yy[i] * a_exp + 4.0 * g_x_xxz_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_yy_yz[i] = -2.0 * g_x_z_yy_yz[i] * a_exp + 4.0 * g_x_xxz_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_yy_zz[i] = -2.0 * g_x_z_yy_zz[i] * a_exp + 4.0 * g_x_xxz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_x_0_0_0_xz_yz_xx, g_x_x_0_0_0_xz_yz_xy, g_x_x_0_0_0_xz_yz_xz, g_x_x_0_0_0_xz_yz_yy, g_x_x_0_0_0_xz_yz_yz, g_x_x_0_0_0_xz_yz_zz, g_x_xxz_yz_xx, g_x_xxz_yz_xy, g_x_xxz_yz_xz, g_x_xxz_yz_yy, g_x_xxz_yz_yz, g_x_xxz_yz_zz, g_x_z_yz_xx, g_x_z_yz_xy, g_x_z_yz_xz, g_x_z_yz_yy, g_x_z_yz_yz, g_x_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xz_yz_xx[i] = -2.0 * g_x_z_yz_xx[i] * a_exp + 4.0 * g_x_xxz_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_yz_xy[i] = -2.0 * g_x_z_yz_xy[i] * a_exp + 4.0 * g_x_xxz_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_yz_xz[i] = -2.0 * g_x_z_yz_xz[i] * a_exp + 4.0 * g_x_xxz_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_yz_yy[i] = -2.0 * g_x_z_yz_yy[i] * a_exp + 4.0 * g_x_xxz_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_yz_yz[i] = -2.0 * g_x_z_yz_yz[i] * a_exp + 4.0 * g_x_xxz_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_yz_zz[i] = -2.0 * g_x_z_yz_zz[i] * a_exp + 4.0 * g_x_xxz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_x_0_0_0_xz_zz_xx, g_x_x_0_0_0_xz_zz_xy, g_x_x_0_0_0_xz_zz_xz, g_x_x_0_0_0_xz_zz_yy, g_x_x_0_0_0_xz_zz_yz, g_x_x_0_0_0_xz_zz_zz, g_x_xxz_zz_xx, g_x_xxz_zz_xy, g_x_xxz_zz_xz, g_x_xxz_zz_yy, g_x_xxz_zz_yz, g_x_xxz_zz_zz, g_x_z_zz_xx, g_x_z_zz_xy, g_x_z_zz_xz, g_x_z_zz_yy, g_x_z_zz_yz, g_x_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xz_zz_xx[i] = -2.0 * g_x_z_zz_xx[i] * a_exp + 4.0 * g_x_xxz_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_zz_xy[i] = -2.0 * g_x_z_zz_xy[i] * a_exp + 4.0 * g_x_xxz_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_zz_xz[i] = -2.0 * g_x_z_zz_xz[i] * a_exp + 4.0 * g_x_xxz_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_zz_yy[i] = -2.0 * g_x_z_zz_yy[i] * a_exp + 4.0 * g_x_xxz_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_zz_yz[i] = -2.0 * g_x_z_zz_yz[i] * a_exp + 4.0 * g_x_xxz_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_zz_zz[i] = -2.0 * g_x_z_zz_zz[i] * a_exp + 4.0 * g_x_xxz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_x_0_0_0_yy_xx_xx, g_x_x_0_0_0_yy_xx_xy, g_x_x_0_0_0_yy_xx_xz, g_x_x_0_0_0_yy_xx_yy, g_x_x_0_0_0_yy_xx_yz, g_x_x_0_0_0_yy_xx_zz, g_x_xyy_xx_xx, g_x_xyy_xx_xy, g_x_xyy_xx_xz, g_x_xyy_xx_yy, g_x_xyy_xx_yz, g_x_xyy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yy_xx_xx[i] = 4.0 * g_x_xyy_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_xx_xy[i] = 4.0 * g_x_xyy_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_xx_xz[i] = 4.0 * g_x_xyy_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_xx_yy[i] = 4.0 * g_x_xyy_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_xx_yz[i] = 4.0 * g_x_xyy_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_xx_zz[i] = 4.0 * g_x_xyy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_x_0_0_0_yy_xy_xx, g_x_x_0_0_0_yy_xy_xy, g_x_x_0_0_0_yy_xy_xz, g_x_x_0_0_0_yy_xy_yy, g_x_x_0_0_0_yy_xy_yz, g_x_x_0_0_0_yy_xy_zz, g_x_xyy_xy_xx, g_x_xyy_xy_xy, g_x_xyy_xy_xz, g_x_xyy_xy_yy, g_x_xyy_xy_yz, g_x_xyy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yy_xy_xx[i] = 4.0 * g_x_xyy_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_xy_xy[i] = 4.0 * g_x_xyy_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_xy_xz[i] = 4.0 * g_x_xyy_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_xy_yy[i] = 4.0 * g_x_xyy_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_xy_yz[i] = 4.0 * g_x_xyy_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_xy_zz[i] = 4.0 * g_x_xyy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_x_0_0_0_yy_xz_xx, g_x_x_0_0_0_yy_xz_xy, g_x_x_0_0_0_yy_xz_xz, g_x_x_0_0_0_yy_xz_yy, g_x_x_0_0_0_yy_xz_yz, g_x_x_0_0_0_yy_xz_zz, g_x_xyy_xz_xx, g_x_xyy_xz_xy, g_x_xyy_xz_xz, g_x_xyy_xz_yy, g_x_xyy_xz_yz, g_x_xyy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yy_xz_xx[i] = 4.0 * g_x_xyy_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_xz_xy[i] = 4.0 * g_x_xyy_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_xz_xz[i] = 4.0 * g_x_xyy_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_xz_yy[i] = 4.0 * g_x_xyy_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_xz_yz[i] = 4.0 * g_x_xyy_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_xz_zz[i] = 4.0 * g_x_xyy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_x_0_0_0_yy_yy_xx, g_x_x_0_0_0_yy_yy_xy, g_x_x_0_0_0_yy_yy_xz, g_x_x_0_0_0_yy_yy_yy, g_x_x_0_0_0_yy_yy_yz, g_x_x_0_0_0_yy_yy_zz, g_x_xyy_yy_xx, g_x_xyy_yy_xy, g_x_xyy_yy_xz, g_x_xyy_yy_yy, g_x_xyy_yy_yz, g_x_xyy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yy_yy_xx[i] = 4.0 * g_x_xyy_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_yy_xy[i] = 4.0 * g_x_xyy_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_yy_xz[i] = 4.0 * g_x_xyy_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_yy_yy[i] = 4.0 * g_x_xyy_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_yy_yz[i] = 4.0 * g_x_xyy_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_yy_zz[i] = 4.0 * g_x_xyy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_x_0_0_0_yy_yz_xx, g_x_x_0_0_0_yy_yz_xy, g_x_x_0_0_0_yy_yz_xz, g_x_x_0_0_0_yy_yz_yy, g_x_x_0_0_0_yy_yz_yz, g_x_x_0_0_0_yy_yz_zz, g_x_xyy_yz_xx, g_x_xyy_yz_xy, g_x_xyy_yz_xz, g_x_xyy_yz_yy, g_x_xyy_yz_yz, g_x_xyy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yy_yz_xx[i] = 4.0 * g_x_xyy_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_yz_xy[i] = 4.0 * g_x_xyy_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_yz_xz[i] = 4.0 * g_x_xyy_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_yz_yy[i] = 4.0 * g_x_xyy_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_yz_yz[i] = 4.0 * g_x_xyy_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_yz_zz[i] = 4.0 * g_x_xyy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_x_0_0_0_yy_zz_xx, g_x_x_0_0_0_yy_zz_xy, g_x_x_0_0_0_yy_zz_xz, g_x_x_0_0_0_yy_zz_yy, g_x_x_0_0_0_yy_zz_yz, g_x_x_0_0_0_yy_zz_zz, g_x_xyy_zz_xx, g_x_xyy_zz_xy, g_x_xyy_zz_xz, g_x_xyy_zz_yy, g_x_xyy_zz_yz, g_x_xyy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yy_zz_xx[i] = 4.0 * g_x_xyy_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_zz_xy[i] = 4.0 * g_x_xyy_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_zz_xz[i] = 4.0 * g_x_xyy_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_zz_yy[i] = 4.0 * g_x_xyy_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_zz_yz[i] = 4.0 * g_x_xyy_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_zz_zz[i] = 4.0 * g_x_xyy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_x_0_0_0_yz_xx_xx, g_x_x_0_0_0_yz_xx_xy, g_x_x_0_0_0_yz_xx_xz, g_x_x_0_0_0_yz_xx_yy, g_x_x_0_0_0_yz_xx_yz, g_x_x_0_0_0_yz_xx_zz, g_x_xyz_xx_xx, g_x_xyz_xx_xy, g_x_xyz_xx_xz, g_x_xyz_xx_yy, g_x_xyz_xx_yz, g_x_xyz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yz_xx_xx[i] = 4.0 * g_x_xyz_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_xx_xy[i] = 4.0 * g_x_xyz_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_xx_xz[i] = 4.0 * g_x_xyz_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_xx_yy[i] = 4.0 * g_x_xyz_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_xx_yz[i] = 4.0 * g_x_xyz_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_xx_zz[i] = 4.0 * g_x_xyz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_x_0_0_0_yz_xy_xx, g_x_x_0_0_0_yz_xy_xy, g_x_x_0_0_0_yz_xy_xz, g_x_x_0_0_0_yz_xy_yy, g_x_x_0_0_0_yz_xy_yz, g_x_x_0_0_0_yz_xy_zz, g_x_xyz_xy_xx, g_x_xyz_xy_xy, g_x_xyz_xy_xz, g_x_xyz_xy_yy, g_x_xyz_xy_yz, g_x_xyz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yz_xy_xx[i] = 4.0 * g_x_xyz_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_xy_xy[i] = 4.0 * g_x_xyz_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_xy_xz[i] = 4.0 * g_x_xyz_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_xy_yy[i] = 4.0 * g_x_xyz_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_xy_yz[i] = 4.0 * g_x_xyz_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_xy_zz[i] = 4.0 * g_x_xyz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_x_0_0_0_yz_xz_xx, g_x_x_0_0_0_yz_xz_xy, g_x_x_0_0_0_yz_xz_xz, g_x_x_0_0_0_yz_xz_yy, g_x_x_0_0_0_yz_xz_yz, g_x_x_0_0_0_yz_xz_zz, g_x_xyz_xz_xx, g_x_xyz_xz_xy, g_x_xyz_xz_xz, g_x_xyz_xz_yy, g_x_xyz_xz_yz, g_x_xyz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yz_xz_xx[i] = 4.0 * g_x_xyz_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_xz_xy[i] = 4.0 * g_x_xyz_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_xz_xz[i] = 4.0 * g_x_xyz_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_xz_yy[i] = 4.0 * g_x_xyz_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_xz_yz[i] = 4.0 * g_x_xyz_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_xz_zz[i] = 4.0 * g_x_xyz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_x_0_0_0_yz_yy_xx, g_x_x_0_0_0_yz_yy_xy, g_x_x_0_0_0_yz_yy_xz, g_x_x_0_0_0_yz_yy_yy, g_x_x_0_0_0_yz_yy_yz, g_x_x_0_0_0_yz_yy_zz, g_x_xyz_yy_xx, g_x_xyz_yy_xy, g_x_xyz_yy_xz, g_x_xyz_yy_yy, g_x_xyz_yy_yz, g_x_xyz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yz_yy_xx[i] = 4.0 * g_x_xyz_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_yy_xy[i] = 4.0 * g_x_xyz_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_yy_xz[i] = 4.0 * g_x_xyz_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_yy_yy[i] = 4.0 * g_x_xyz_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_yy_yz[i] = 4.0 * g_x_xyz_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_yy_zz[i] = 4.0 * g_x_xyz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_x_0_0_0_yz_yz_xx, g_x_x_0_0_0_yz_yz_xy, g_x_x_0_0_0_yz_yz_xz, g_x_x_0_0_0_yz_yz_yy, g_x_x_0_0_0_yz_yz_yz, g_x_x_0_0_0_yz_yz_zz, g_x_xyz_yz_xx, g_x_xyz_yz_xy, g_x_xyz_yz_xz, g_x_xyz_yz_yy, g_x_xyz_yz_yz, g_x_xyz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yz_yz_xx[i] = 4.0 * g_x_xyz_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_yz_xy[i] = 4.0 * g_x_xyz_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_yz_xz[i] = 4.0 * g_x_xyz_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_yz_yy[i] = 4.0 * g_x_xyz_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_yz_yz[i] = 4.0 * g_x_xyz_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_yz_zz[i] = 4.0 * g_x_xyz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_x_0_0_0_yz_zz_xx, g_x_x_0_0_0_yz_zz_xy, g_x_x_0_0_0_yz_zz_xz, g_x_x_0_0_0_yz_zz_yy, g_x_x_0_0_0_yz_zz_yz, g_x_x_0_0_0_yz_zz_zz, g_x_xyz_zz_xx, g_x_xyz_zz_xy, g_x_xyz_zz_xz, g_x_xyz_zz_yy, g_x_xyz_zz_yz, g_x_xyz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yz_zz_xx[i] = 4.0 * g_x_xyz_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_zz_xy[i] = 4.0 * g_x_xyz_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_zz_xz[i] = 4.0 * g_x_xyz_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_zz_yy[i] = 4.0 * g_x_xyz_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_zz_yz[i] = 4.0 * g_x_xyz_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_zz_zz[i] = 4.0 * g_x_xyz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_x_0_0_0_zz_xx_xx, g_x_x_0_0_0_zz_xx_xy, g_x_x_0_0_0_zz_xx_xz, g_x_x_0_0_0_zz_xx_yy, g_x_x_0_0_0_zz_xx_yz, g_x_x_0_0_0_zz_xx_zz, g_x_xzz_xx_xx, g_x_xzz_xx_xy, g_x_xzz_xx_xz, g_x_xzz_xx_yy, g_x_xzz_xx_yz, g_x_xzz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_zz_xx_xx[i] = 4.0 * g_x_xzz_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_xx_xy[i] = 4.0 * g_x_xzz_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_xx_xz[i] = 4.0 * g_x_xzz_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_xx_yy[i] = 4.0 * g_x_xzz_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_xx_yz[i] = 4.0 * g_x_xzz_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_xx_zz[i] = 4.0 * g_x_xzz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_x_0_0_0_zz_xy_xx, g_x_x_0_0_0_zz_xy_xy, g_x_x_0_0_0_zz_xy_xz, g_x_x_0_0_0_zz_xy_yy, g_x_x_0_0_0_zz_xy_yz, g_x_x_0_0_0_zz_xy_zz, g_x_xzz_xy_xx, g_x_xzz_xy_xy, g_x_xzz_xy_xz, g_x_xzz_xy_yy, g_x_xzz_xy_yz, g_x_xzz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_zz_xy_xx[i] = 4.0 * g_x_xzz_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_xy_xy[i] = 4.0 * g_x_xzz_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_xy_xz[i] = 4.0 * g_x_xzz_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_xy_yy[i] = 4.0 * g_x_xzz_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_xy_yz[i] = 4.0 * g_x_xzz_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_xy_zz[i] = 4.0 * g_x_xzz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_x_0_0_0_zz_xz_xx, g_x_x_0_0_0_zz_xz_xy, g_x_x_0_0_0_zz_xz_xz, g_x_x_0_0_0_zz_xz_yy, g_x_x_0_0_0_zz_xz_yz, g_x_x_0_0_0_zz_xz_zz, g_x_xzz_xz_xx, g_x_xzz_xz_xy, g_x_xzz_xz_xz, g_x_xzz_xz_yy, g_x_xzz_xz_yz, g_x_xzz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_zz_xz_xx[i] = 4.0 * g_x_xzz_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_xz_xy[i] = 4.0 * g_x_xzz_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_xz_xz[i] = 4.0 * g_x_xzz_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_xz_yy[i] = 4.0 * g_x_xzz_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_xz_yz[i] = 4.0 * g_x_xzz_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_xz_zz[i] = 4.0 * g_x_xzz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_x_0_0_0_zz_yy_xx, g_x_x_0_0_0_zz_yy_xy, g_x_x_0_0_0_zz_yy_xz, g_x_x_0_0_0_zz_yy_yy, g_x_x_0_0_0_zz_yy_yz, g_x_x_0_0_0_zz_yy_zz, g_x_xzz_yy_xx, g_x_xzz_yy_xy, g_x_xzz_yy_xz, g_x_xzz_yy_yy, g_x_xzz_yy_yz, g_x_xzz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_zz_yy_xx[i] = 4.0 * g_x_xzz_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_yy_xy[i] = 4.0 * g_x_xzz_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_yy_xz[i] = 4.0 * g_x_xzz_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_yy_yy[i] = 4.0 * g_x_xzz_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_yy_yz[i] = 4.0 * g_x_xzz_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_yy_zz[i] = 4.0 * g_x_xzz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_x_0_0_0_zz_yz_xx, g_x_x_0_0_0_zz_yz_xy, g_x_x_0_0_0_zz_yz_xz, g_x_x_0_0_0_zz_yz_yy, g_x_x_0_0_0_zz_yz_yz, g_x_x_0_0_0_zz_yz_zz, g_x_xzz_yz_xx, g_x_xzz_yz_xy, g_x_xzz_yz_xz, g_x_xzz_yz_yy, g_x_xzz_yz_yz, g_x_xzz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_zz_yz_xx[i] = 4.0 * g_x_xzz_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_yz_xy[i] = 4.0 * g_x_xzz_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_yz_xz[i] = 4.0 * g_x_xzz_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_yz_yy[i] = 4.0 * g_x_xzz_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_yz_yz[i] = 4.0 * g_x_xzz_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_yz_zz[i] = 4.0 * g_x_xzz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_x_0_0_0_zz_zz_xx, g_x_x_0_0_0_zz_zz_xy, g_x_x_0_0_0_zz_zz_xz, g_x_x_0_0_0_zz_zz_yy, g_x_x_0_0_0_zz_zz_yz, g_x_x_0_0_0_zz_zz_zz, g_x_xzz_zz_xx, g_x_xzz_zz_xy, g_x_xzz_zz_xz, g_x_xzz_zz_yy, g_x_xzz_zz_yz, g_x_xzz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_zz_zz_xx[i] = 4.0 * g_x_xzz_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_zz_xy[i] = 4.0 * g_x_xzz_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_zz_xz[i] = 4.0 * g_x_xzz_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_zz_yy[i] = 4.0 * g_x_xzz_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_zz_yz[i] = 4.0 * g_x_xzz_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_zz_zz[i] = 4.0 * g_x_xzz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_xxy_xx_xx, g_x_xxy_xx_xy, g_x_xxy_xx_xz, g_x_xxy_xx_yy, g_x_xxy_xx_yz, g_x_xxy_xx_zz, g_x_y_0_0_0_xx_xx_xx, g_x_y_0_0_0_xx_xx_xy, g_x_y_0_0_0_xx_xx_xz, g_x_y_0_0_0_xx_xx_yy, g_x_y_0_0_0_xx_xx_yz, g_x_y_0_0_0_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xx_xx_xx[i] = 4.0 * g_x_xxy_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_xx_xy[i] = 4.0 * g_x_xxy_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_xx_xz[i] = 4.0 * g_x_xxy_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_xx_yy[i] = 4.0 * g_x_xxy_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_xx_yz[i] = 4.0 * g_x_xxy_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_xx_zz[i] = 4.0 * g_x_xxy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_xxy_xy_xx, g_x_xxy_xy_xy, g_x_xxy_xy_xz, g_x_xxy_xy_yy, g_x_xxy_xy_yz, g_x_xxy_xy_zz, g_x_y_0_0_0_xx_xy_xx, g_x_y_0_0_0_xx_xy_xy, g_x_y_0_0_0_xx_xy_xz, g_x_y_0_0_0_xx_xy_yy, g_x_y_0_0_0_xx_xy_yz, g_x_y_0_0_0_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xx_xy_xx[i] = 4.0 * g_x_xxy_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_xy_xy[i] = 4.0 * g_x_xxy_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_xy_xz[i] = 4.0 * g_x_xxy_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_xy_yy[i] = 4.0 * g_x_xxy_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_xy_yz[i] = 4.0 * g_x_xxy_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_xy_zz[i] = 4.0 * g_x_xxy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_xxy_xz_xx, g_x_xxy_xz_xy, g_x_xxy_xz_xz, g_x_xxy_xz_yy, g_x_xxy_xz_yz, g_x_xxy_xz_zz, g_x_y_0_0_0_xx_xz_xx, g_x_y_0_0_0_xx_xz_xy, g_x_y_0_0_0_xx_xz_xz, g_x_y_0_0_0_xx_xz_yy, g_x_y_0_0_0_xx_xz_yz, g_x_y_0_0_0_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xx_xz_xx[i] = 4.0 * g_x_xxy_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_xz_xy[i] = 4.0 * g_x_xxy_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_xz_xz[i] = 4.0 * g_x_xxy_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_xz_yy[i] = 4.0 * g_x_xxy_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_xz_yz[i] = 4.0 * g_x_xxy_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_xz_zz[i] = 4.0 * g_x_xxy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_xxy_yy_xx, g_x_xxy_yy_xy, g_x_xxy_yy_xz, g_x_xxy_yy_yy, g_x_xxy_yy_yz, g_x_xxy_yy_zz, g_x_y_0_0_0_xx_yy_xx, g_x_y_0_0_0_xx_yy_xy, g_x_y_0_0_0_xx_yy_xz, g_x_y_0_0_0_xx_yy_yy, g_x_y_0_0_0_xx_yy_yz, g_x_y_0_0_0_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xx_yy_xx[i] = 4.0 * g_x_xxy_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_yy_xy[i] = 4.0 * g_x_xxy_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_yy_xz[i] = 4.0 * g_x_xxy_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_yy_yy[i] = 4.0 * g_x_xxy_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_yy_yz[i] = 4.0 * g_x_xxy_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_yy_zz[i] = 4.0 * g_x_xxy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_xxy_yz_xx, g_x_xxy_yz_xy, g_x_xxy_yz_xz, g_x_xxy_yz_yy, g_x_xxy_yz_yz, g_x_xxy_yz_zz, g_x_y_0_0_0_xx_yz_xx, g_x_y_0_0_0_xx_yz_xy, g_x_y_0_0_0_xx_yz_xz, g_x_y_0_0_0_xx_yz_yy, g_x_y_0_0_0_xx_yz_yz, g_x_y_0_0_0_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xx_yz_xx[i] = 4.0 * g_x_xxy_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_yz_xy[i] = 4.0 * g_x_xxy_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_yz_xz[i] = 4.0 * g_x_xxy_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_yz_yy[i] = 4.0 * g_x_xxy_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_yz_yz[i] = 4.0 * g_x_xxy_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_yz_zz[i] = 4.0 * g_x_xxy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_xxy_zz_xx, g_x_xxy_zz_xy, g_x_xxy_zz_xz, g_x_xxy_zz_yy, g_x_xxy_zz_yz, g_x_xxy_zz_zz, g_x_y_0_0_0_xx_zz_xx, g_x_y_0_0_0_xx_zz_xy, g_x_y_0_0_0_xx_zz_xz, g_x_y_0_0_0_xx_zz_yy, g_x_y_0_0_0_xx_zz_yz, g_x_y_0_0_0_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xx_zz_xx[i] = 4.0 * g_x_xxy_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_zz_xy[i] = 4.0 * g_x_xxy_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_zz_xz[i] = 4.0 * g_x_xxy_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_zz_yy[i] = 4.0 * g_x_xxy_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_zz_yz[i] = 4.0 * g_x_xxy_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_zz_zz[i] = 4.0 * g_x_xxy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_x_xx_xx, g_x_x_xx_xy, g_x_x_xx_xz, g_x_x_xx_yy, g_x_x_xx_yz, g_x_x_xx_zz, g_x_xyy_xx_xx, g_x_xyy_xx_xy, g_x_xyy_xx_xz, g_x_xyy_xx_yy, g_x_xyy_xx_yz, g_x_xyy_xx_zz, g_x_y_0_0_0_xy_xx_xx, g_x_y_0_0_0_xy_xx_xy, g_x_y_0_0_0_xy_xx_xz, g_x_y_0_0_0_xy_xx_yy, g_x_y_0_0_0_xy_xx_yz, g_x_y_0_0_0_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xy_xx_xx[i] = -2.0 * g_x_x_xx_xx[i] * a_exp + 4.0 * g_x_xyy_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_xx_xy[i] = -2.0 * g_x_x_xx_xy[i] * a_exp + 4.0 * g_x_xyy_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_xx_xz[i] = -2.0 * g_x_x_xx_xz[i] * a_exp + 4.0 * g_x_xyy_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_xx_yy[i] = -2.0 * g_x_x_xx_yy[i] * a_exp + 4.0 * g_x_xyy_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_xx_yz[i] = -2.0 * g_x_x_xx_yz[i] * a_exp + 4.0 * g_x_xyy_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_xx_zz[i] = -2.0 * g_x_x_xx_zz[i] * a_exp + 4.0 * g_x_xyy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_x_xy_xx, g_x_x_xy_xy, g_x_x_xy_xz, g_x_x_xy_yy, g_x_x_xy_yz, g_x_x_xy_zz, g_x_xyy_xy_xx, g_x_xyy_xy_xy, g_x_xyy_xy_xz, g_x_xyy_xy_yy, g_x_xyy_xy_yz, g_x_xyy_xy_zz, g_x_y_0_0_0_xy_xy_xx, g_x_y_0_0_0_xy_xy_xy, g_x_y_0_0_0_xy_xy_xz, g_x_y_0_0_0_xy_xy_yy, g_x_y_0_0_0_xy_xy_yz, g_x_y_0_0_0_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xy_xy_xx[i] = -2.0 * g_x_x_xy_xx[i] * a_exp + 4.0 * g_x_xyy_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_xy_xy[i] = -2.0 * g_x_x_xy_xy[i] * a_exp + 4.0 * g_x_xyy_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_xy_xz[i] = -2.0 * g_x_x_xy_xz[i] * a_exp + 4.0 * g_x_xyy_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_xy_yy[i] = -2.0 * g_x_x_xy_yy[i] * a_exp + 4.0 * g_x_xyy_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_xy_yz[i] = -2.0 * g_x_x_xy_yz[i] * a_exp + 4.0 * g_x_xyy_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_xy_zz[i] = -2.0 * g_x_x_xy_zz[i] * a_exp + 4.0 * g_x_xyy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_x_xz_xx, g_x_x_xz_xy, g_x_x_xz_xz, g_x_x_xz_yy, g_x_x_xz_yz, g_x_x_xz_zz, g_x_xyy_xz_xx, g_x_xyy_xz_xy, g_x_xyy_xz_xz, g_x_xyy_xz_yy, g_x_xyy_xz_yz, g_x_xyy_xz_zz, g_x_y_0_0_0_xy_xz_xx, g_x_y_0_0_0_xy_xz_xy, g_x_y_0_0_0_xy_xz_xz, g_x_y_0_0_0_xy_xz_yy, g_x_y_0_0_0_xy_xz_yz, g_x_y_0_0_0_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xy_xz_xx[i] = -2.0 * g_x_x_xz_xx[i] * a_exp + 4.0 * g_x_xyy_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_xz_xy[i] = -2.0 * g_x_x_xz_xy[i] * a_exp + 4.0 * g_x_xyy_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_xz_xz[i] = -2.0 * g_x_x_xz_xz[i] * a_exp + 4.0 * g_x_xyy_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_xz_yy[i] = -2.0 * g_x_x_xz_yy[i] * a_exp + 4.0 * g_x_xyy_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_xz_yz[i] = -2.0 * g_x_x_xz_yz[i] * a_exp + 4.0 * g_x_xyy_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_xz_zz[i] = -2.0 * g_x_x_xz_zz[i] * a_exp + 4.0 * g_x_xyy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_x_yy_xx, g_x_x_yy_xy, g_x_x_yy_xz, g_x_x_yy_yy, g_x_x_yy_yz, g_x_x_yy_zz, g_x_xyy_yy_xx, g_x_xyy_yy_xy, g_x_xyy_yy_xz, g_x_xyy_yy_yy, g_x_xyy_yy_yz, g_x_xyy_yy_zz, g_x_y_0_0_0_xy_yy_xx, g_x_y_0_0_0_xy_yy_xy, g_x_y_0_0_0_xy_yy_xz, g_x_y_0_0_0_xy_yy_yy, g_x_y_0_0_0_xy_yy_yz, g_x_y_0_0_0_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xy_yy_xx[i] = -2.0 * g_x_x_yy_xx[i] * a_exp + 4.0 * g_x_xyy_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_yy_xy[i] = -2.0 * g_x_x_yy_xy[i] * a_exp + 4.0 * g_x_xyy_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_yy_xz[i] = -2.0 * g_x_x_yy_xz[i] * a_exp + 4.0 * g_x_xyy_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_yy_yy[i] = -2.0 * g_x_x_yy_yy[i] * a_exp + 4.0 * g_x_xyy_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_yy_yz[i] = -2.0 * g_x_x_yy_yz[i] * a_exp + 4.0 * g_x_xyy_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_yy_zz[i] = -2.0 * g_x_x_yy_zz[i] * a_exp + 4.0 * g_x_xyy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_x_yz_xx, g_x_x_yz_xy, g_x_x_yz_xz, g_x_x_yz_yy, g_x_x_yz_yz, g_x_x_yz_zz, g_x_xyy_yz_xx, g_x_xyy_yz_xy, g_x_xyy_yz_xz, g_x_xyy_yz_yy, g_x_xyy_yz_yz, g_x_xyy_yz_zz, g_x_y_0_0_0_xy_yz_xx, g_x_y_0_0_0_xy_yz_xy, g_x_y_0_0_0_xy_yz_xz, g_x_y_0_0_0_xy_yz_yy, g_x_y_0_0_0_xy_yz_yz, g_x_y_0_0_0_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xy_yz_xx[i] = -2.0 * g_x_x_yz_xx[i] * a_exp + 4.0 * g_x_xyy_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_yz_xy[i] = -2.0 * g_x_x_yz_xy[i] * a_exp + 4.0 * g_x_xyy_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_yz_xz[i] = -2.0 * g_x_x_yz_xz[i] * a_exp + 4.0 * g_x_xyy_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_yz_yy[i] = -2.0 * g_x_x_yz_yy[i] * a_exp + 4.0 * g_x_xyy_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_yz_yz[i] = -2.0 * g_x_x_yz_yz[i] * a_exp + 4.0 * g_x_xyy_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_yz_zz[i] = -2.0 * g_x_x_yz_zz[i] * a_exp + 4.0 * g_x_xyy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_x_zz_xx, g_x_x_zz_xy, g_x_x_zz_xz, g_x_x_zz_yy, g_x_x_zz_yz, g_x_x_zz_zz, g_x_xyy_zz_xx, g_x_xyy_zz_xy, g_x_xyy_zz_xz, g_x_xyy_zz_yy, g_x_xyy_zz_yz, g_x_xyy_zz_zz, g_x_y_0_0_0_xy_zz_xx, g_x_y_0_0_0_xy_zz_xy, g_x_y_0_0_0_xy_zz_xz, g_x_y_0_0_0_xy_zz_yy, g_x_y_0_0_0_xy_zz_yz, g_x_y_0_0_0_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xy_zz_xx[i] = -2.0 * g_x_x_zz_xx[i] * a_exp + 4.0 * g_x_xyy_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_zz_xy[i] = -2.0 * g_x_x_zz_xy[i] * a_exp + 4.0 * g_x_xyy_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_zz_xz[i] = -2.0 * g_x_x_zz_xz[i] * a_exp + 4.0 * g_x_xyy_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_zz_yy[i] = -2.0 * g_x_x_zz_yy[i] * a_exp + 4.0 * g_x_xyy_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_zz_yz[i] = -2.0 * g_x_x_zz_yz[i] * a_exp + 4.0 * g_x_xyy_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_zz_zz[i] = -2.0 * g_x_x_zz_zz[i] * a_exp + 4.0 * g_x_xyy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_xyz_xx_xx, g_x_xyz_xx_xy, g_x_xyz_xx_xz, g_x_xyz_xx_yy, g_x_xyz_xx_yz, g_x_xyz_xx_zz, g_x_y_0_0_0_xz_xx_xx, g_x_y_0_0_0_xz_xx_xy, g_x_y_0_0_0_xz_xx_xz, g_x_y_0_0_0_xz_xx_yy, g_x_y_0_0_0_xz_xx_yz, g_x_y_0_0_0_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xz_xx_xx[i] = 4.0 * g_x_xyz_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_xx_xy[i] = 4.0 * g_x_xyz_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_xx_xz[i] = 4.0 * g_x_xyz_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_xx_yy[i] = 4.0 * g_x_xyz_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_xx_yz[i] = 4.0 * g_x_xyz_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_xx_zz[i] = 4.0 * g_x_xyz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_xyz_xy_xx, g_x_xyz_xy_xy, g_x_xyz_xy_xz, g_x_xyz_xy_yy, g_x_xyz_xy_yz, g_x_xyz_xy_zz, g_x_y_0_0_0_xz_xy_xx, g_x_y_0_0_0_xz_xy_xy, g_x_y_0_0_0_xz_xy_xz, g_x_y_0_0_0_xz_xy_yy, g_x_y_0_0_0_xz_xy_yz, g_x_y_0_0_0_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xz_xy_xx[i] = 4.0 * g_x_xyz_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_xy_xy[i] = 4.0 * g_x_xyz_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_xy_xz[i] = 4.0 * g_x_xyz_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_xy_yy[i] = 4.0 * g_x_xyz_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_xy_yz[i] = 4.0 * g_x_xyz_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_xy_zz[i] = 4.0 * g_x_xyz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_xyz_xz_xx, g_x_xyz_xz_xy, g_x_xyz_xz_xz, g_x_xyz_xz_yy, g_x_xyz_xz_yz, g_x_xyz_xz_zz, g_x_y_0_0_0_xz_xz_xx, g_x_y_0_0_0_xz_xz_xy, g_x_y_0_0_0_xz_xz_xz, g_x_y_0_0_0_xz_xz_yy, g_x_y_0_0_0_xz_xz_yz, g_x_y_0_0_0_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xz_xz_xx[i] = 4.0 * g_x_xyz_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_xz_xy[i] = 4.0 * g_x_xyz_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_xz_xz[i] = 4.0 * g_x_xyz_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_xz_yy[i] = 4.0 * g_x_xyz_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_xz_yz[i] = 4.0 * g_x_xyz_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_xz_zz[i] = 4.0 * g_x_xyz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_xyz_yy_xx, g_x_xyz_yy_xy, g_x_xyz_yy_xz, g_x_xyz_yy_yy, g_x_xyz_yy_yz, g_x_xyz_yy_zz, g_x_y_0_0_0_xz_yy_xx, g_x_y_0_0_0_xz_yy_xy, g_x_y_0_0_0_xz_yy_xz, g_x_y_0_0_0_xz_yy_yy, g_x_y_0_0_0_xz_yy_yz, g_x_y_0_0_0_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xz_yy_xx[i] = 4.0 * g_x_xyz_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_yy_xy[i] = 4.0 * g_x_xyz_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_yy_xz[i] = 4.0 * g_x_xyz_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_yy_yy[i] = 4.0 * g_x_xyz_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_yy_yz[i] = 4.0 * g_x_xyz_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_yy_zz[i] = 4.0 * g_x_xyz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_xyz_yz_xx, g_x_xyz_yz_xy, g_x_xyz_yz_xz, g_x_xyz_yz_yy, g_x_xyz_yz_yz, g_x_xyz_yz_zz, g_x_y_0_0_0_xz_yz_xx, g_x_y_0_0_0_xz_yz_xy, g_x_y_0_0_0_xz_yz_xz, g_x_y_0_0_0_xz_yz_yy, g_x_y_0_0_0_xz_yz_yz, g_x_y_0_0_0_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xz_yz_xx[i] = 4.0 * g_x_xyz_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_yz_xy[i] = 4.0 * g_x_xyz_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_yz_xz[i] = 4.0 * g_x_xyz_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_yz_yy[i] = 4.0 * g_x_xyz_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_yz_yz[i] = 4.0 * g_x_xyz_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_yz_zz[i] = 4.0 * g_x_xyz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_xyz_zz_xx, g_x_xyz_zz_xy, g_x_xyz_zz_xz, g_x_xyz_zz_yy, g_x_xyz_zz_yz, g_x_xyz_zz_zz, g_x_y_0_0_0_xz_zz_xx, g_x_y_0_0_0_xz_zz_xy, g_x_y_0_0_0_xz_zz_xz, g_x_y_0_0_0_xz_zz_yy, g_x_y_0_0_0_xz_zz_yz, g_x_y_0_0_0_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xz_zz_xx[i] = 4.0 * g_x_xyz_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_zz_xy[i] = 4.0 * g_x_xyz_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_zz_xz[i] = 4.0 * g_x_xyz_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_zz_yy[i] = 4.0 * g_x_xyz_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_zz_yz[i] = 4.0 * g_x_xyz_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_zz_zz[i] = 4.0 * g_x_xyz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_x_y_0_0_0_yy_xx_xx, g_x_y_0_0_0_yy_xx_xy, g_x_y_0_0_0_yy_xx_xz, g_x_y_0_0_0_yy_xx_yy, g_x_y_0_0_0_yy_xx_yz, g_x_y_0_0_0_yy_xx_zz, g_x_y_xx_xx, g_x_y_xx_xy, g_x_y_xx_xz, g_x_y_xx_yy, g_x_y_xx_yz, g_x_y_xx_zz, g_x_yyy_xx_xx, g_x_yyy_xx_xy, g_x_yyy_xx_xz, g_x_yyy_xx_yy, g_x_yyy_xx_yz, g_x_yyy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yy_xx_xx[i] = -4.0 * g_x_y_xx_xx[i] * a_exp + 4.0 * g_x_yyy_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_xx_xy[i] = -4.0 * g_x_y_xx_xy[i] * a_exp + 4.0 * g_x_yyy_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_xx_xz[i] = -4.0 * g_x_y_xx_xz[i] * a_exp + 4.0 * g_x_yyy_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_xx_yy[i] = -4.0 * g_x_y_xx_yy[i] * a_exp + 4.0 * g_x_yyy_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_xx_yz[i] = -4.0 * g_x_y_xx_yz[i] * a_exp + 4.0 * g_x_yyy_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_xx_zz[i] = -4.0 * g_x_y_xx_zz[i] * a_exp + 4.0 * g_x_yyy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_x_y_0_0_0_yy_xy_xx, g_x_y_0_0_0_yy_xy_xy, g_x_y_0_0_0_yy_xy_xz, g_x_y_0_0_0_yy_xy_yy, g_x_y_0_0_0_yy_xy_yz, g_x_y_0_0_0_yy_xy_zz, g_x_y_xy_xx, g_x_y_xy_xy, g_x_y_xy_xz, g_x_y_xy_yy, g_x_y_xy_yz, g_x_y_xy_zz, g_x_yyy_xy_xx, g_x_yyy_xy_xy, g_x_yyy_xy_xz, g_x_yyy_xy_yy, g_x_yyy_xy_yz, g_x_yyy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yy_xy_xx[i] = -4.0 * g_x_y_xy_xx[i] * a_exp + 4.0 * g_x_yyy_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_xy_xy[i] = -4.0 * g_x_y_xy_xy[i] * a_exp + 4.0 * g_x_yyy_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_xy_xz[i] = -4.0 * g_x_y_xy_xz[i] * a_exp + 4.0 * g_x_yyy_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_xy_yy[i] = -4.0 * g_x_y_xy_yy[i] * a_exp + 4.0 * g_x_yyy_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_xy_yz[i] = -4.0 * g_x_y_xy_yz[i] * a_exp + 4.0 * g_x_yyy_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_xy_zz[i] = -4.0 * g_x_y_xy_zz[i] * a_exp + 4.0 * g_x_yyy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_x_y_0_0_0_yy_xz_xx, g_x_y_0_0_0_yy_xz_xy, g_x_y_0_0_0_yy_xz_xz, g_x_y_0_0_0_yy_xz_yy, g_x_y_0_0_0_yy_xz_yz, g_x_y_0_0_0_yy_xz_zz, g_x_y_xz_xx, g_x_y_xz_xy, g_x_y_xz_xz, g_x_y_xz_yy, g_x_y_xz_yz, g_x_y_xz_zz, g_x_yyy_xz_xx, g_x_yyy_xz_xy, g_x_yyy_xz_xz, g_x_yyy_xz_yy, g_x_yyy_xz_yz, g_x_yyy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yy_xz_xx[i] = -4.0 * g_x_y_xz_xx[i] * a_exp + 4.0 * g_x_yyy_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_xz_xy[i] = -4.0 * g_x_y_xz_xy[i] * a_exp + 4.0 * g_x_yyy_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_xz_xz[i] = -4.0 * g_x_y_xz_xz[i] * a_exp + 4.0 * g_x_yyy_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_xz_yy[i] = -4.0 * g_x_y_xz_yy[i] * a_exp + 4.0 * g_x_yyy_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_xz_yz[i] = -4.0 * g_x_y_xz_yz[i] * a_exp + 4.0 * g_x_yyy_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_xz_zz[i] = -4.0 * g_x_y_xz_zz[i] * a_exp + 4.0 * g_x_yyy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_x_y_0_0_0_yy_yy_xx, g_x_y_0_0_0_yy_yy_xy, g_x_y_0_0_0_yy_yy_xz, g_x_y_0_0_0_yy_yy_yy, g_x_y_0_0_0_yy_yy_yz, g_x_y_0_0_0_yy_yy_zz, g_x_y_yy_xx, g_x_y_yy_xy, g_x_y_yy_xz, g_x_y_yy_yy, g_x_y_yy_yz, g_x_y_yy_zz, g_x_yyy_yy_xx, g_x_yyy_yy_xy, g_x_yyy_yy_xz, g_x_yyy_yy_yy, g_x_yyy_yy_yz, g_x_yyy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yy_yy_xx[i] = -4.0 * g_x_y_yy_xx[i] * a_exp + 4.0 * g_x_yyy_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_yy_xy[i] = -4.0 * g_x_y_yy_xy[i] * a_exp + 4.0 * g_x_yyy_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_yy_xz[i] = -4.0 * g_x_y_yy_xz[i] * a_exp + 4.0 * g_x_yyy_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_yy_yy[i] = -4.0 * g_x_y_yy_yy[i] * a_exp + 4.0 * g_x_yyy_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_yy_yz[i] = -4.0 * g_x_y_yy_yz[i] * a_exp + 4.0 * g_x_yyy_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_yy_zz[i] = -4.0 * g_x_y_yy_zz[i] * a_exp + 4.0 * g_x_yyy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_x_y_0_0_0_yy_yz_xx, g_x_y_0_0_0_yy_yz_xy, g_x_y_0_0_0_yy_yz_xz, g_x_y_0_0_0_yy_yz_yy, g_x_y_0_0_0_yy_yz_yz, g_x_y_0_0_0_yy_yz_zz, g_x_y_yz_xx, g_x_y_yz_xy, g_x_y_yz_xz, g_x_y_yz_yy, g_x_y_yz_yz, g_x_y_yz_zz, g_x_yyy_yz_xx, g_x_yyy_yz_xy, g_x_yyy_yz_xz, g_x_yyy_yz_yy, g_x_yyy_yz_yz, g_x_yyy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yy_yz_xx[i] = -4.0 * g_x_y_yz_xx[i] * a_exp + 4.0 * g_x_yyy_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_yz_xy[i] = -4.0 * g_x_y_yz_xy[i] * a_exp + 4.0 * g_x_yyy_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_yz_xz[i] = -4.0 * g_x_y_yz_xz[i] * a_exp + 4.0 * g_x_yyy_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_yz_yy[i] = -4.0 * g_x_y_yz_yy[i] * a_exp + 4.0 * g_x_yyy_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_yz_yz[i] = -4.0 * g_x_y_yz_yz[i] * a_exp + 4.0 * g_x_yyy_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_yz_zz[i] = -4.0 * g_x_y_yz_zz[i] * a_exp + 4.0 * g_x_yyy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_x_y_0_0_0_yy_zz_xx, g_x_y_0_0_0_yy_zz_xy, g_x_y_0_0_0_yy_zz_xz, g_x_y_0_0_0_yy_zz_yy, g_x_y_0_0_0_yy_zz_yz, g_x_y_0_0_0_yy_zz_zz, g_x_y_zz_xx, g_x_y_zz_xy, g_x_y_zz_xz, g_x_y_zz_yy, g_x_y_zz_yz, g_x_y_zz_zz, g_x_yyy_zz_xx, g_x_yyy_zz_xy, g_x_yyy_zz_xz, g_x_yyy_zz_yy, g_x_yyy_zz_yz, g_x_yyy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yy_zz_xx[i] = -4.0 * g_x_y_zz_xx[i] * a_exp + 4.0 * g_x_yyy_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_zz_xy[i] = -4.0 * g_x_y_zz_xy[i] * a_exp + 4.0 * g_x_yyy_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_zz_xz[i] = -4.0 * g_x_y_zz_xz[i] * a_exp + 4.0 * g_x_yyy_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_zz_yy[i] = -4.0 * g_x_y_zz_yy[i] * a_exp + 4.0 * g_x_yyy_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_zz_yz[i] = -4.0 * g_x_y_zz_yz[i] * a_exp + 4.0 * g_x_yyy_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_zz_zz[i] = -4.0 * g_x_y_zz_zz[i] * a_exp + 4.0 * g_x_yyy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_x_y_0_0_0_yz_xx_xx, g_x_y_0_0_0_yz_xx_xy, g_x_y_0_0_0_yz_xx_xz, g_x_y_0_0_0_yz_xx_yy, g_x_y_0_0_0_yz_xx_yz, g_x_y_0_0_0_yz_xx_zz, g_x_yyz_xx_xx, g_x_yyz_xx_xy, g_x_yyz_xx_xz, g_x_yyz_xx_yy, g_x_yyz_xx_yz, g_x_yyz_xx_zz, g_x_z_xx_xx, g_x_z_xx_xy, g_x_z_xx_xz, g_x_z_xx_yy, g_x_z_xx_yz, g_x_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yz_xx_xx[i] = -2.0 * g_x_z_xx_xx[i] * a_exp + 4.0 * g_x_yyz_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_xx_xy[i] = -2.0 * g_x_z_xx_xy[i] * a_exp + 4.0 * g_x_yyz_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_xx_xz[i] = -2.0 * g_x_z_xx_xz[i] * a_exp + 4.0 * g_x_yyz_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_xx_yy[i] = -2.0 * g_x_z_xx_yy[i] * a_exp + 4.0 * g_x_yyz_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_xx_yz[i] = -2.0 * g_x_z_xx_yz[i] * a_exp + 4.0 * g_x_yyz_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_xx_zz[i] = -2.0 * g_x_z_xx_zz[i] * a_exp + 4.0 * g_x_yyz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_x_y_0_0_0_yz_xy_xx, g_x_y_0_0_0_yz_xy_xy, g_x_y_0_0_0_yz_xy_xz, g_x_y_0_0_0_yz_xy_yy, g_x_y_0_0_0_yz_xy_yz, g_x_y_0_0_0_yz_xy_zz, g_x_yyz_xy_xx, g_x_yyz_xy_xy, g_x_yyz_xy_xz, g_x_yyz_xy_yy, g_x_yyz_xy_yz, g_x_yyz_xy_zz, g_x_z_xy_xx, g_x_z_xy_xy, g_x_z_xy_xz, g_x_z_xy_yy, g_x_z_xy_yz, g_x_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yz_xy_xx[i] = -2.0 * g_x_z_xy_xx[i] * a_exp + 4.0 * g_x_yyz_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_xy_xy[i] = -2.0 * g_x_z_xy_xy[i] * a_exp + 4.0 * g_x_yyz_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_xy_xz[i] = -2.0 * g_x_z_xy_xz[i] * a_exp + 4.0 * g_x_yyz_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_xy_yy[i] = -2.0 * g_x_z_xy_yy[i] * a_exp + 4.0 * g_x_yyz_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_xy_yz[i] = -2.0 * g_x_z_xy_yz[i] * a_exp + 4.0 * g_x_yyz_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_xy_zz[i] = -2.0 * g_x_z_xy_zz[i] * a_exp + 4.0 * g_x_yyz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_x_y_0_0_0_yz_xz_xx, g_x_y_0_0_0_yz_xz_xy, g_x_y_0_0_0_yz_xz_xz, g_x_y_0_0_0_yz_xz_yy, g_x_y_0_0_0_yz_xz_yz, g_x_y_0_0_0_yz_xz_zz, g_x_yyz_xz_xx, g_x_yyz_xz_xy, g_x_yyz_xz_xz, g_x_yyz_xz_yy, g_x_yyz_xz_yz, g_x_yyz_xz_zz, g_x_z_xz_xx, g_x_z_xz_xy, g_x_z_xz_xz, g_x_z_xz_yy, g_x_z_xz_yz, g_x_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yz_xz_xx[i] = -2.0 * g_x_z_xz_xx[i] * a_exp + 4.0 * g_x_yyz_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_xz_xy[i] = -2.0 * g_x_z_xz_xy[i] * a_exp + 4.0 * g_x_yyz_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_xz_xz[i] = -2.0 * g_x_z_xz_xz[i] * a_exp + 4.0 * g_x_yyz_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_xz_yy[i] = -2.0 * g_x_z_xz_yy[i] * a_exp + 4.0 * g_x_yyz_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_xz_yz[i] = -2.0 * g_x_z_xz_yz[i] * a_exp + 4.0 * g_x_yyz_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_xz_zz[i] = -2.0 * g_x_z_xz_zz[i] * a_exp + 4.0 * g_x_yyz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_x_y_0_0_0_yz_yy_xx, g_x_y_0_0_0_yz_yy_xy, g_x_y_0_0_0_yz_yy_xz, g_x_y_0_0_0_yz_yy_yy, g_x_y_0_0_0_yz_yy_yz, g_x_y_0_0_0_yz_yy_zz, g_x_yyz_yy_xx, g_x_yyz_yy_xy, g_x_yyz_yy_xz, g_x_yyz_yy_yy, g_x_yyz_yy_yz, g_x_yyz_yy_zz, g_x_z_yy_xx, g_x_z_yy_xy, g_x_z_yy_xz, g_x_z_yy_yy, g_x_z_yy_yz, g_x_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yz_yy_xx[i] = -2.0 * g_x_z_yy_xx[i] * a_exp + 4.0 * g_x_yyz_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_yy_xy[i] = -2.0 * g_x_z_yy_xy[i] * a_exp + 4.0 * g_x_yyz_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_yy_xz[i] = -2.0 * g_x_z_yy_xz[i] * a_exp + 4.0 * g_x_yyz_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_yy_yy[i] = -2.0 * g_x_z_yy_yy[i] * a_exp + 4.0 * g_x_yyz_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_yy_yz[i] = -2.0 * g_x_z_yy_yz[i] * a_exp + 4.0 * g_x_yyz_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_yy_zz[i] = -2.0 * g_x_z_yy_zz[i] * a_exp + 4.0 * g_x_yyz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_x_y_0_0_0_yz_yz_xx, g_x_y_0_0_0_yz_yz_xy, g_x_y_0_0_0_yz_yz_xz, g_x_y_0_0_0_yz_yz_yy, g_x_y_0_0_0_yz_yz_yz, g_x_y_0_0_0_yz_yz_zz, g_x_yyz_yz_xx, g_x_yyz_yz_xy, g_x_yyz_yz_xz, g_x_yyz_yz_yy, g_x_yyz_yz_yz, g_x_yyz_yz_zz, g_x_z_yz_xx, g_x_z_yz_xy, g_x_z_yz_xz, g_x_z_yz_yy, g_x_z_yz_yz, g_x_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yz_yz_xx[i] = -2.0 * g_x_z_yz_xx[i] * a_exp + 4.0 * g_x_yyz_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_yz_xy[i] = -2.0 * g_x_z_yz_xy[i] * a_exp + 4.0 * g_x_yyz_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_yz_xz[i] = -2.0 * g_x_z_yz_xz[i] * a_exp + 4.0 * g_x_yyz_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_yz_yy[i] = -2.0 * g_x_z_yz_yy[i] * a_exp + 4.0 * g_x_yyz_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_yz_yz[i] = -2.0 * g_x_z_yz_yz[i] * a_exp + 4.0 * g_x_yyz_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_yz_zz[i] = -2.0 * g_x_z_yz_zz[i] * a_exp + 4.0 * g_x_yyz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_x_y_0_0_0_yz_zz_xx, g_x_y_0_0_0_yz_zz_xy, g_x_y_0_0_0_yz_zz_xz, g_x_y_0_0_0_yz_zz_yy, g_x_y_0_0_0_yz_zz_yz, g_x_y_0_0_0_yz_zz_zz, g_x_yyz_zz_xx, g_x_yyz_zz_xy, g_x_yyz_zz_xz, g_x_yyz_zz_yy, g_x_yyz_zz_yz, g_x_yyz_zz_zz, g_x_z_zz_xx, g_x_z_zz_xy, g_x_z_zz_xz, g_x_z_zz_yy, g_x_z_zz_yz, g_x_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yz_zz_xx[i] = -2.0 * g_x_z_zz_xx[i] * a_exp + 4.0 * g_x_yyz_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_zz_xy[i] = -2.0 * g_x_z_zz_xy[i] * a_exp + 4.0 * g_x_yyz_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_zz_xz[i] = -2.0 * g_x_z_zz_xz[i] * a_exp + 4.0 * g_x_yyz_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_zz_yy[i] = -2.0 * g_x_z_zz_yy[i] * a_exp + 4.0 * g_x_yyz_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_zz_yz[i] = -2.0 * g_x_z_zz_yz[i] * a_exp + 4.0 * g_x_yyz_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_zz_zz[i] = -2.0 * g_x_z_zz_zz[i] * a_exp + 4.0 * g_x_yyz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_x_y_0_0_0_zz_xx_xx, g_x_y_0_0_0_zz_xx_xy, g_x_y_0_0_0_zz_xx_xz, g_x_y_0_0_0_zz_xx_yy, g_x_y_0_0_0_zz_xx_yz, g_x_y_0_0_0_zz_xx_zz, g_x_yzz_xx_xx, g_x_yzz_xx_xy, g_x_yzz_xx_xz, g_x_yzz_xx_yy, g_x_yzz_xx_yz, g_x_yzz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_zz_xx_xx[i] = 4.0 * g_x_yzz_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_xx_xy[i] = 4.0 * g_x_yzz_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_xx_xz[i] = 4.0 * g_x_yzz_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_xx_yy[i] = 4.0 * g_x_yzz_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_xx_yz[i] = 4.0 * g_x_yzz_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_xx_zz[i] = 4.0 * g_x_yzz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_x_y_0_0_0_zz_xy_xx, g_x_y_0_0_0_zz_xy_xy, g_x_y_0_0_0_zz_xy_xz, g_x_y_0_0_0_zz_xy_yy, g_x_y_0_0_0_zz_xy_yz, g_x_y_0_0_0_zz_xy_zz, g_x_yzz_xy_xx, g_x_yzz_xy_xy, g_x_yzz_xy_xz, g_x_yzz_xy_yy, g_x_yzz_xy_yz, g_x_yzz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_zz_xy_xx[i] = 4.0 * g_x_yzz_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_xy_xy[i] = 4.0 * g_x_yzz_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_xy_xz[i] = 4.0 * g_x_yzz_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_xy_yy[i] = 4.0 * g_x_yzz_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_xy_yz[i] = 4.0 * g_x_yzz_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_xy_zz[i] = 4.0 * g_x_yzz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_x_y_0_0_0_zz_xz_xx, g_x_y_0_0_0_zz_xz_xy, g_x_y_0_0_0_zz_xz_xz, g_x_y_0_0_0_zz_xz_yy, g_x_y_0_0_0_zz_xz_yz, g_x_y_0_0_0_zz_xz_zz, g_x_yzz_xz_xx, g_x_yzz_xz_xy, g_x_yzz_xz_xz, g_x_yzz_xz_yy, g_x_yzz_xz_yz, g_x_yzz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_zz_xz_xx[i] = 4.0 * g_x_yzz_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_xz_xy[i] = 4.0 * g_x_yzz_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_xz_xz[i] = 4.0 * g_x_yzz_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_xz_yy[i] = 4.0 * g_x_yzz_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_xz_yz[i] = 4.0 * g_x_yzz_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_xz_zz[i] = 4.0 * g_x_yzz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_x_y_0_0_0_zz_yy_xx, g_x_y_0_0_0_zz_yy_xy, g_x_y_0_0_0_zz_yy_xz, g_x_y_0_0_0_zz_yy_yy, g_x_y_0_0_0_zz_yy_yz, g_x_y_0_0_0_zz_yy_zz, g_x_yzz_yy_xx, g_x_yzz_yy_xy, g_x_yzz_yy_xz, g_x_yzz_yy_yy, g_x_yzz_yy_yz, g_x_yzz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_zz_yy_xx[i] = 4.0 * g_x_yzz_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_yy_xy[i] = 4.0 * g_x_yzz_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_yy_xz[i] = 4.0 * g_x_yzz_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_yy_yy[i] = 4.0 * g_x_yzz_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_yy_yz[i] = 4.0 * g_x_yzz_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_yy_zz[i] = 4.0 * g_x_yzz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_x_y_0_0_0_zz_yz_xx, g_x_y_0_0_0_zz_yz_xy, g_x_y_0_0_0_zz_yz_xz, g_x_y_0_0_0_zz_yz_yy, g_x_y_0_0_0_zz_yz_yz, g_x_y_0_0_0_zz_yz_zz, g_x_yzz_yz_xx, g_x_yzz_yz_xy, g_x_yzz_yz_xz, g_x_yzz_yz_yy, g_x_yzz_yz_yz, g_x_yzz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_zz_yz_xx[i] = 4.0 * g_x_yzz_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_yz_xy[i] = 4.0 * g_x_yzz_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_yz_xz[i] = 4.0 * g_x_yzz_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_yz_yy[i] = 4.0 * g_x_yzz_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_yz_yz[i] = 4.0 * g_x_yzz_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_yz_zz[i] = 4.0 * g_x_yzz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_x_y_0_0_0_zz_zz_xx, g_x_y_0_0_0_zz_zz_xy, g_x_y_0_0_0_zz_zz_xz, g_x_y_0_0_0_zz_zz_yy, g_x_y_0_0_0_zz_zz_yz, g_x_y_0_0_0_zz_zz_zz, g_x_yzz_zz_xx, g_x_yzz_zz_xy, g_x_yzz_zz_xz, g_x_yzz_zz_yy, g_x_yzz_zz_yz, g_x_yzz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_zz_zz_xx[i] = 4.0 * g_x_yzz_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_zz_xy[i] = 4.0 * g_x_yzz_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_zz_xz[i] = 4.0 * g_x_yzz_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_zz_yy[i] = 4.0 * g_x_yzz_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_zz_yz[i] = 4.0 * g_x_yzz_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_zz_zz[i] = 4.0 * g_x_yzz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_x_xxz_xx_xx, g_x_xxz_xx_xy, g_x_xxz_xx_xz, g_x_xxz_xx_yy, g_x_xxz_xx_yz, g_x_xxz_xx_zz, g_x_z_0_0_0_xx_xx_xx, g_x_z_0_0_0_xx_xx_xy, g_x_z_0_0_0_xx_xx_xz, g_x_z_0_0_0_xx_xx_yy, g_x_z_0_0_0_xx_xx_yz, g_x_z_0_0_0_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xx_xx_xx[i] = 4.0 * g_x_xxz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_xx_xy[i] = 4.0 * g_x_xxz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_xx_xz[i] = 4.0 * g_x_xxz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_xx_yy[i] = 4.0 * g_x_xxz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_xx_yz[i] = 4.0 * g_x_xxz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_xx_zz[i] = 4.0 * g_x_xxz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_x_xxz_xy_xx, g_x_xxz_xy_xy, g_x_xxz_xy_xz, g_x_xxz_xy_yy, g_x_xxz_xy_yz, g_x_xxz_xy_zz, g_x_z_0_0_0_xx_xy_xx, g_x_z_0_0_0_xx_xy_xy, g_x_z_0_0_0_xx_xy_xz, g_x_z_0_0_0_xx_xy_yy, g_x_z_0_0_0_xx_xy_yz, g_x_z_0_0_0_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xx_xy_xx[i] = 4.0 * g_x_xxz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_xy_xy[i] = 4.0 * g_x_xxz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_xy_xz[i] = 4.0 * g_x_xxz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_xy_yy[i] = 4.0 * g_x_xxz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_xy_yz[i] = 4.0 * g_x_xxz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_xy_zz[i] = 4.0 * g_x_xxz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_x_xxz_xz_xx, g_x_xxz_xz_xy, g_x_xxz_xz_xz, g_x_xxz_xz_yy, g_x_xxz_xz_yz, g_x_xxz_xz_zz, g_x_z_0_0_0_xx_xz_xx, g_x_z_0_0_0_xx_xz_xy, g_x_z_0_0_0_xx_xz_xz, g_x_z_0_0_0_xx_xz_yy, g_x_z_0_0_0_xx_xz_yz, g_x_z_0_0_0_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xx_xz_xx[i] = 4.0 * g_x_xxz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_xz_xy[i] = 4.0 * g_x_xxz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_xz_xz[i] = 4.0 * g_x_xxz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_xz_yy[i] = 4.0 * g_x_xxz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_xz_yz[i] = 4.0 * g_x_xxz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_xz_zz[i] = 4.0 * g_x_xxz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_x_xxz_yy_xx, g_x_xxz_yy_xy, g_x_xxz_yy_xz, g_x_xxz_yy_yy, g_x_xxz_yy_yz, g_x_xxz_yy_zz, g_x_z_0_0_0_xx_yy_xx, g_x_z_0_0_0_xx_yy_xy, g_x_z_0_0_0_xx_yy_xz, g_x_z_0_0_0_xx_yy_yy, g_x_z_0_0_0_xx_yy_yz, g_x_z_0_0_0_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xx_yy_xx[i] = 4.0 * g_x_xxz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_yy_xy[i] = 4.0 * g_x_xxz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_yy_xz[i] = 4.0 * g_x_xxz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_yy_yy[i] = 4.0 * g_x_xxz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_yy_yz[i] = 4.0 * g_x_xxz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_yy_zz[i] = 4.0 * g_x_xxz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_x_xxz_yz_xx, g_x_xxz_yz_xy, g_x_xxz_yz_xz, g_x_xxz_yz_yy, g_x_xxz_yz_yz, g_x_xxz_yz_zz, g_x_z_0_0_0_xx_yz_xx, g_x_z_0_0_0_xx_yz_xy, g_x_z_0_0_0_xx_yz_xz, g_x_z_0_0_0_xx_yz_yy, g_x_z_0_0_0_xx_yz_yz, g_x_z_0_0_0_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xx_yz_xx[i] = 4.0 * g_x_xxz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_yz_xy[i] = 4.0 * g_x_xxz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_yz_xz[i] = 4.0 * g_x_xxz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_yz_yy[i] = 4.0 * g_x_xxz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_yz_yz[i] = 4.0 * g_x_xxz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_yz_zz[i] = 4.0 * g_x_xxz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_x_xxz_zz_xx, g_x_xxz_zz_xy, g_x_xxz_zz_xz, g_x_xxz_zz_yy, g_x_xxz_zz_yz, g_x_xxz_zz_zz, g_x_z_0_0_0_xx_zz_xx, g_x_z_0_0_0_xx_zz_xy, g_x_z_0_0_0_xx_zz_xz, g_x_z_0_0_0_xx_zz_yy, g_x_z_0_0_0_xx_zz_yz, g_x_z_0_0_0_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xx_zz_xx[i] = 4.0 * g_x_xxz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_zz_xy[i] = 4.0 * g_x_xxz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_zz_xz[i] = 4.0 * g_x_xxz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_zz_yy[i] = 4.0 * g_x_xxz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_zz_yz[i] = 4.0 * g_x_xxz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_zz_zz[i] = 4.0 * g_x_xxz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_x_xyz_xx_xx, g_x_xyz_xx_xy, g_x_xyz_xx_xz, g_x_xyz_xx_yy, g_x_xyz_xx_yz, g_x_xyz_xx_zz, g_x_z_0_0_0_xy_xx_xx, g_x_z_0_0_0_xy_xx_xy, g_x_z_0_0_0_xy_xx_xz, g_x_z_0_0_0_xy_xx_yy, g_x_z_0_0_0_xy_xx_yz, g_x_z_0_0_0_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xy_xx_xx[i] = 4.0 * g_x_xyz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_xx_xy[i] = 4.0 * g_x_xyz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_xx_xz[i] = 4.0 * g_x_xyz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_xx_yy[i] = 4.0 * g_x_xyz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_xx_yz[i] = 4.0 * g_x_xyz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_xx_zz[i] = 4.0 * g_x_xyz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_x_xyz_xy_xx, g_x_xyz_xy_xy, g_x_xyz_xy_xz, g_x_xyz_xy_yy, g_x_xyz_xy_yz, g_x_xyz_xy_zz, g_x_z_0_0_0_xy_xy_xx, g_x_z_0_0_0_xy_xy_xy, g_x_z_0_0_0_xy_xy_xz, g_x_z_0_0_0_xy_xy_yy, g_x_z_0_0_0_xy_xy_yz, g_x_z_0_0_0_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xy_xy_xx[i] = 4.0 * g_x_xyz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_xy_xy[i] = 4.0 * g_x_xyz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_xy_xz[i] = 4.0 * g_x_xyz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_xy_yy[i] = 4.0 * g_x_xyz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_xy_yz[i] = 4.0 * g_x_xyz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_xy_zz[i] = 4.0 * g_x_xyz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_x_xyz_xz_xx, g_x_xyz_xz_xy, g_x_xyz_xz_xz, g_x_xyz_xz_yy, g_x_xyz_xz_yz, g_x_xyz_xz_zz, g_x_z_0_0_0_xy_xz_xx, g_x_z_0_0_0_xy_xz_xy, g_x_z_0_0_0_xy_xz_xz, g_x_z_0_0_0_xy_xz_yy, g_x_z_0_0_0_xy_xz_yz, g_x_z_0_0_0_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xy_xz_xx[i] = 4.0 * g_x_xyz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_xz_xy[i] = 4.0 * g_x_xyz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_xz_xz[i] = 4.0 * g_x_xyz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_xz_yy[i] = 4.0 * g_x_xyz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_xz_yz[i] = 4.0 * g_x_xyz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_xz_zz[i] = 4.0 * g_x_xyz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_x_xyz_yy_xx, g_x_xyz_yy_xy, g_x_xyz_yy_xz, g_x_xyz_yy_yy, g_x_xyz_yy_yz, g_x_xyz_yy_zz, g_x_z_0_0_0_xy_yy_xx, g_x_z_0_0_0_xy_yy_xy, g_x_z_0_0_0_xy_yy_xz, g_x_z_0_0_0_xy_yy_yy, g_x_z_0_0_0_xy_yy_yz, g_x_z_0_0_0_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xy_yy_xx[i] = 4.0 * g_x_xyz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_yy_xy[i] = 4.0 * g_x_xyz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_yy_xz[i] = 4.0 * g_x_xyz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_yy_yy[i] = 4.0 * g_x_xyz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_yy_yz[i] = 4.0 * g_x_xyz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_yy_zz[i] = 4.0 * g_x_xyz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_x_xyz_yz_xx, g_x_xyz_yz_xy, g_x_xyz_yz_xz, g_x_xyz_yz_yy, g_x_xyz_yz_yz, g_x_xyz_yz_zz, g_x_z_0_0_0_xy_yz_xx, g_x_z_0_0_0_xy_yz_xy, g_x_z_0_0_0_xy_yz_xz, g_x_z_0_0_0_xy_yz_yy, g_x_z_0_0_0_xy_yz_yz, g_x_z_0_0_0_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xy_yz_xx[i] = 4.0 * g_x_xyz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_yz_xy[i] = 4.0 * g_x_xyz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_yz_xz[i] = 4.0 * g_x_xyz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_yz_yy[i] = 4.0 * g_x_xyz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_yz_yz[i] = 4.0 * g_x_xyz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_yz_zz[i] = 4.0 * g_x_xyz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_x_xyz_zz_xx, g_x_xyz_zz_xy, g_x_xyz_zz_xz, g_x_xyz_zz_yy, g_x_xyz_zz_yz, g_x_xyz_zz_zz, g_x_z_0_0_0_xy_zz_xx, g_x_z_0_0_0_xy_zz_xy, g_x_z_0_0_0_xy_zz_xz, g_x_z_0_0_0_xy_zz_yy, g_x_z_0_0_0_xy_zz_yz, g_x_z_0_0_0_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xy_zz_xx[i] = 4.0 * g_x_xyz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_zz_xy[i] = 4.0 * g_x_xyz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_zz_xz[i] = 4.0 * g_x_xyz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_zz_yy[i] = 4.0 * g_x_xyz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_zz_yz[i] = 4.0 * g_x_xyz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_zz_zz[i] = 4.0 * g_x_xyz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_x_x_xx_xx, g_x_x_xx_xy, g_x_x_xx_xz, g_x_x_xx_yy, g_x_x_xx_yz, g_x_x_xx_zz, g_x_xzz_xx_xx, g_x_xzz_xx_xy, g_x_xzz_xx_xz, g_x_xzz_xx_yy, g_x_xzz_xx_yz, g_x_xzz_xx_zz, g_x_z_0_0_0_xz_xx_xx, g_x_z_0_0_0_xz_xx_xy, g_x_z_0_0_0_xz_xx_xz, g_x_z_0_0_0_xz_xx_yy, g_x_z_0_0_0_xz_xx_yz, g_x_z_0_0_0_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xz_xx_xx[i] = -2.0 * g_x_x_xx_xx[i] * a_exp + 4.0 * g_x_xzz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_xx_xy[i] = -2.0 * g_x_x_xx_xy[i] * a_exp + 4.0 * g_x_xzz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_xx_xz[i] = -2.0 * g_x_x_xx_xz[i] * a_exp + 4.0 * g_x_xzz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_xx_yy[i] = -2.0 * g_x_x_xx_yy[i] * a_exp + 4.0 * g_x_xzz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_xx_yz[i] = -2.0 * g_x_x_xx_yz[i] * a_exp + 4.0 * g_x_xzz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_xx_zz[i] = -2.0 * g_x_x_xx_zz[i] * a_exp + 4.0 * g_x_xzz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_x_x_xy_xx, g_x_x_xy_xy, g_x_x_xy_xz, g_x_x_xy_yy, g_x_x_xy_yz, g_x_x_xy_zz, g_x_xzz_xy_xx, g_x_xzz_xy_xy, g_x_xzz_xy_xz, g_x_xzz_xy_yy, g_x_xzz_xy_yz, g_x_xzz_xy_zz, g_x_z_0_0_0_xz_xy_xx, g_x_z_0_0_0_xz_xy_xy, g_x_z_0_0_0_xz_xy_xz, g_x_z_0_0_0_xz_xy_yy, g_x_z_0_0_0_xz_xy_yz, g_x_z_0_0_0_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xz_xy_xx[i] = -2.0 * g_x_x_xy_xx[i] * a_exp + 4.0 * g_x_xzz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_xy_xy[i] = -2.0 * g_x_x_xy_xy[i] * a_exp + 4.0 * g_x_xzz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_xy_xz[i] = -2.0 * g_x_x_xy_xz[i] * a_exp + 4.0 * g_x_xzz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_xy_yy[i] = -2.0 * g_x_x_xy_yy[i] * a_exp + 4.0 * g_x_xzz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_xy_yz[i] = -2.0 * g_x_x_xy_yz[i] * a_exp + 4.0 * g_x_xzz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_xy_zz[i] = -2.0 * g_x_x_xy_zz[i] * a_exp + 4.0 * g_x_xzz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_x_x_xz_xx, g_x_x_xz_xy, g_x_x_xz_xz, g_x_x_xz_yy, g_x_x_xz_yz, g_x_x_xz_zz, g_x_xzz_xz_xx, g_x_xzz_xz_xy, g_x_xzz_xz_xz, g_x_xzz_xz_yy, g_x_xzz_xz_yz, g_x_xzz_xz_zz, g_x_z_0_0_0_xz_xz_xx, g_x_z_0_0_0_xz_xz_xy, g_x_z_0_0_0_xz_xz_xz, g_x_z_0_0_0_xz_xz_yy, g_x_z_0_0_0_xz_xz_yz, g_x_z_0_0_0_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xz_xz_xx[i] = -2.0 * g_x_x_xz_xx[i] * a_exp + 4.0 * g_x_xzz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_xz_xy[i] = -2.0 * g_x_x_xz_xy[i] * a_exp + 4.0 * g_x_xzz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_xz_xz[i] = -2.0 * g_x_x_xz_xz[i] * a_exp + 4.0 * g_x_xzz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_xz_yy[i] = -2.0 * g_x_x_xz_yy[i] * a_exp + 4.0 * g_x_xzz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_xz_yz[i] = -2.0 * g_x_x_xz_yz[i] * a_exp + 4.0 * g_x_xzz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_xz_zz[i] = -2.0 * g_x_x_xz_zz[i] * a_exp + 4.0 * g_x_xzz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_x_x_yy_xx, g_x_x_yy_xy, g_x_x_yy_xz, g_x_x_yy_yy, g_x_x_yy_yz, g_x_x_yy_zz, g_x_xzz_yy_xx, g_x_xzz_yy_xy, g_x_xzz_yy_xz, g_x_xzz_yy_yy, g_x_xzz_yy_yz, g_x_xzz_yy_zz, g_x_z_0_0_0_xz_yy_xx, g_x_z_0_0_0_xz_yy_xy, g_x_z_0_0_0_xz_yy_xz, g_x_z_0_0_0_xz_yy_yy, g_x_z_0_0_0_xz_yy_yz, g_x_z_0_0_0_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xz_yy_xx[i] = -2.0 * g_x_x_yy_xx[i] * a_exp + 4.0 * g_x_xzz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_yy_xy[i] = -2.0 * g_x_x_yy_xy[i] * a_exp + 4.0 * g_x_xzz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_yy_xz[i] = -2.0 * g_x_x_yy_xz[i] * a_exp + 4.0 * g_x_xzz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_yy_yy[i] = -2.0 * g_x_x_yy_yy[i] * a_exp + 4.0 * g_x_xzz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_yy_yz[i] = -2.0 * g_x_x_yy_yz[i] * a_exp + 4.0 * g_x_xzz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_yy_zz[i] = -2.0 * g_x_x_yy_zz[i] * a_exp + 4.0 * g_x_xzz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_x_x_yz_xx, g_x_x_yz_xy, g_x_x_yz_xz, g_x_x_yz_yy, g_x_x_yz_yz, g_x_x_yz_zz, g_x_xzz_yz_xx, g_x_xzz_yz_xy, g_x_xzz_yz_xz, g_x_xzz_yz_yy, g_x_xzz_yz_yz, g_x_xzz_yz_zz, g_x_z_0_0_0_xz_yz_xx, g_x_z_0_0_0_xz_yz_xy, g_x_z_0_0_0_xz_yz_xz, g_x_z_0_0_0_xz_yz_yy, g_x_z_0_0_0_xz_yz_yz, g_x_z_0_0_0_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xz_yz_xx[i] = -2.0 * g_x_x_yz_xx[i] * a_exp + 4.0 * g_x_xzz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_yz_xy[i] = -2.0 * g_x_x_yz_xy[i] * a_exp + 4.0 * g_x_xzz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_yz_xz[i] = -2.0 * g_x_x_yz_xz[i] * a_exp + 4.0 * g_x_xzz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_yz_yy[i] = -2.0 * g_x_x_yz_yy[i] * a_exp + 4.0 * g_x_xzz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_yz_yz[i] = -2.0 * g_x_x_yz_yz[i] * a_exp + 4.0 * g_x_xzz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_yz_zz[i] = -2.0 * g_x_x_yz_zz[i] * a_exp + 4.0 * g_x_xzz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_x_x_zz_xx, g_x_x_zz_xy, g_x_x_zz_xz, g_x_x_zz_yy, g_x_x_zz_yz, g_x_x_zz_zz, g_x_xzz_zz_xx, g_x_xzz_zz_xy, g_x_xzz_zz_xz, g_x_xzz_zz_yy, g_x_xzz_zz_yz, g_x_xzz_zz_zz, g_x_z_0_0_0_xz_zz_xx, g_x_z_0_0_0_xz_zz_xy, g_x_z_0_0_0_xz_zz_xz, g_x_z_0_0_0_xz_zz_yy, g_x_z_0_0_0_xz_zz_yz, g_x_z_0_0_0_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xz_zz_xx[i] = -2.0 * g_x_x_zz_xx[i] * a_exp + 4.0 * g_x_xzz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_zz_xy[i] = -2.0 * g_x_x_zz_xy[i] * a_exp + 4.0 * g_x_xzz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_zz_xz[i] = -2.0 * g_x_x_zz_xz[i] * a_exp + 4.0 * g_x_xzz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_zz_yy[i] = -2.0 * g_x_x_zz_yy[i] * a_exp + 4.0 * g_x_xzz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_zz_yz[i] = -2.0 * g_x_x_zz_yz[i] * a_exp + 4.0 * g_x_xzz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_zz_zz[i] = -2.0 * g_x_x_zz_zz[i] * a_exp + 4.0 * g_x_xzz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_x_yyz_xx_xx, g_x_yyz_xx_xy, g_x_yyz_xx_xz, g_x_yyz_xx_yy, g_x_yyz_xx_yz, g_x_yyz_xx_zz, g_x_z_0_0_0_yy_xx_xx, g_x_z_0_0_0_yy_xx_xy, g_x_z_0_0_0_yy_xx_xz, g_x_z_0_0_0_yy_xx_yy, g_x_z_0_0_0_yy_xx_yz, g_x_z_0_0_0_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yy_xx_xx[i] = 4.0 * g_x_yyz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_xx_xy[i] = 4.0 * g_x_yyz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_xx_xz[i] = 4.0 * g_x_yyz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_xx_yy[i] = 4.0 * g_x_yyz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_xx_yz[i] = 4.0 * g_x_yyz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_xx_zz[i] = 4.0 * g_x_yyz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_x_yyz_xy_xx, g_x_yyz_xy_xy, g_x_yyz_xy_xz, g_x_yyz_xy_yy, g_x_yyz_xy_yz, g_x_yyz_xy_zz, g_x_z_0_0_0_yy_xy_xx, g_x_z_0_0_0_yy_xy_xy, g_x_z_0_0_0_yy_xy_xz, g_x_z_0_0_0_yy_xy_yy, g_x_z_0_0_0_yy_xy_yz, g_x_z_0_0_0_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yy_xy_xx[i] = 4.0 * g_x_yyz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_xy_xy[i] = 4.0 * g_x_yyz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_xy_xz[i] = 4.0 * g_x_yyz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_xy_yy[i] = 4.0 * g_x_yyz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_xy_yz[i] = 4.0 * g_x_yyz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_xy_zz[i] = 4.0 * g_x_yyz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_x_yyz_xz_xx, g_x_yyz_xz_xy, g_x_yyz_xz_xz, g_x_yyz_xz_yy, g_x_yyz_xz_yz, g_x_yyz_xz_zz, g_x_z_0_0_0_yy_xz_xx, g_x_z_0_0_0_yy_xz_xy, g_x_z_0_0_0_yy_xz_xz, g_x_z_0_0_0_yy_xz_yy, g_x_z_0_0_0_yy_xz_yz, g_x_z_0_0_0_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yy_xz_xx[i] = 4.0 * g_x_yyz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_xz_xy[i] = 4.0 * g_x_yyz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_xz_xz[i] = 4.0 * g_x_yyz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_xz_yy[i] = 4.0 * g_x_yyz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_xz_yz[i] = 4.0 * g_x_yyz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_xz_zz[i] = 4.0 * g_x_yyz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_x_yyz_yy_xx, g_x_yyz_yy_xy, g_x_yyz_yy_xz, g_x_yyz_yy_yy, g_x_yyz_yy_yz, g_x_yyz_yy_zz, g_x_z_0_0_0_yy_yy_xx, g_x_z_0_0_0_yy_yy_xy, g_x_z_0_0_0_yy_yy_xz, g_x_z_0_0_0_yy_yy_yy, g_x_z_0_0_0_yy_yy_yz, g_x_z_0_0_0_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yy_yy_xx[i] = 4.0 * g_x_yyz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_yy_xy[i] = 4.0 * g_x_yyz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_yy_xz[i] = 4.0 * g_x_yyz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_yy_yy[i] = 4.0 * g_x_yyz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_yy_yz[i] = 4.0 * g_x_yyz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_yy_zz[i] = 4.0 * g_x_yyz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_x_yyz_yz_xx, g_x_yyz_yz_xy, g_x_yyz_yz_xz, g_x_yyz_yz_yy, g_x_yyz_yz_yz, g_x_yyz_yz_zz, g_x_z_0_0_0_yy_yz_xx, g_x_z_0_0_0_yy_yz_xy, g_x_z_0_0_0_yy_yz_xz, g_x_z_0_0_0_yy_yz_yy, g_x_z_0_0_0_yy_yz_yz, g_x_z_0_0_0_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yy_yz_xx[i] = 4.0 * g_x_yyz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_yz_xy[i] = 4.0 * g_x_yyz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_yz_xz[i] = 4.0 * g_x_yyz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_yz_yy[i] = 4.0 * g_x_yyz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_yz_yz[i] = 4.0 * g_x_yyz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_yz_zz[i] = 4.0 * g_x_yyz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_x_yyz_zz_xx, g_x_yyz_zz_xy, g_x_yyz_zz_xz, g_x_yyz_zz_yy, g_x_yyz_zz_yz, g_x_yyz_zz_zz, g_x_z_0_0_0_yy_zz_xx, g_x_z_0_0_0_yy_zz_xy, g_x_z_0_0_0_yy_zz_xz, g_x_z_0_0_0_yy_zz_yy, g_x_z_0_0_0_yy_zz_yz, g_x_z_0_0_0_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yy_zz_xx[i] = 4.0 * g_x_yyz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_zz_xy[i] = 4.0 * g_x_yyz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_zz_xz[i] = 4.0 * g_x_yyz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_zz_yy[i] = 4.0 * g_x_yyz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_zz_yz[i] = 4.0 * g_x_yyz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_zz_zz[i] = 4.0 * g_x_yyz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_x_y_xx_xx, g_x_y_xx_xy, g_x_y_xx_xz, g_x_y_xx_yy, g_x_y_xx_yz, g_x_y_xx_zz, g_x_yzz_xx_xx, g_x_yzz_xx_xy, g_x_yzz_xx_xz, g_x_yzz_xx_yy, g_x_yzz_xx_yz, g_x_yzz_xx_zz, g_x_z_0_0_0_yz_xx_xx, g_x_z_0_0_0_yz_xx_xy, g_x_z_0_0_0_yz_xx_xz, g_x_z_0_0_0_yz_xx_yy, g_x_z_0_0_0_yz_xx_yz, g_x_z_0_0_0_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yz_xx_xx[i] = -2.0 * g_x_y_xx_xx[i] * a_exp + 4.0 * g_x_yzz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_xx_xy[i] = -2.0 * g_x_y_xx_xy[i] * a_exp + 4.0 * g_x_yzz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_xx_xz[i] = -2.0 * g_x_y_xx_xz[i] * a_exp + 4.0 * g_x_yzz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_xx_yy[i] = -2.0 * g_x_y_xx_yy[i] * a_exp + 4.0 * g_x_yzz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_xx_yz[i] = -2.0 * g_x_y_xx_yz[i] * a_exp + 4.0 * g_x_yzz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_xx_zz[i] = -2.0 * g_x_y_xx_zz[i] * a_exp + 4.0 * g_x_yzz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_x_y_xy_xx, g_x_y_xy_xy, g_x_y_xy_xz, g_x_y_xy_yy, g_x_y_xy_yz, g_x_y_xy_zz, g_x_yzz_xy_xx, g_x_yzz_xy_xy, g_x_yzz_xy_xz, g_x_yzz_xy_yy, g_x_yzz_xy_yz, g_x_yzz_xy_zz, g_x_z_0_0_0_yz_xy_xx, g_x_z_0_0_0_yz_xy_xy, g_x_z_0_0_0_yz_xy_xz, g_x_z_0_0_0_yz_xy_yy, g_x_z_0_0_0_yz_xy_yz, g_x_z_0_0_0_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yz_xy_xx[i] = -2.0 * g_x_y_xy_xx[i] * a_exp + 4.0 * g_x_yzz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_xy_xy[i] = -2.0 * g_x_y_xy_xy[i] * a_exp + 4.0 * g_x_yzz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_xy_xz[i] = -2.0 * g_x_y_xy_xz[i] * a_exp + 4.0 * g_x_yzz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_xy_yy[i] = -2.0 * g_x_y_xy_yy[i] * a_exp + 4.0 * g_x_yzz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_xy_yz[i] = -2.0 * g_x_y_xy_yz[i] * a_exp + 4.0 * g_x_yzz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_xy_zz[i] = -2.0 * g_x_y_xy_zz[i] * a_exp + 4.0 * g_x_yzz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_x_y_xz_xx, g_x_y_xz_xy, g_x_y_xz_xz, g_x_y_xz_yy, g_x_y_xz_yz, g_x_y_xz_zz, g_x_yzz_xz_xx, g_x_yzz_xz_xy, g_x_yzz_xz_xz, g_x_yzz_xz_yy, g_x_yzz_xz_yz, g_x_yzz_xz_zz, g_x_z_0_0_0_yz_xz_xx, g_x_z_0_0_0_yz_xz_xy, g_x_z_0_0_0_yz_xz_xz, g_x_z_0_0_0_yz_xz_yy, g_x_z_0_0_0_yz_xz_yz, g_x_z_0_0_0_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yz_xz_xx[i] = -2.0 * g_x_y_xz_xx[i] * a_exp + 4.0 * g_x_yzz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_xz_xy[i] = -2.0 * g_x_y_xz_xy[i] * a_exp + 4.0 * g_x_yzz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_xz_xz[i] = -2.0 * g_x_y_xz_xz[i] * a_exp + 4.0 * g_x_yzz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_xz_yy[i] = -2.0 * g_x_y_xz_yy[i] * a_exp + 4.0 * g_x_yzz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_xz_yz[i] = -2.0 * g_x_y_xz_yz[i] * a_exp + 4.0 * g_x_yzz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_xz_zz[i] = -2.0 * g_x_y_xz_zz[i] * a_exp + 4.0 * g_x_yzz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_x_y_yy_xx, g_x_y_yy_xy, g_x_y_yy_xz, g_x_y_yy_yy, g_x_y_yy_yz, g_x_y_yy_zz, g_x_yzz_yy_xx, g_x_yzz_yy_xy, g_x_yzz_yy_xz, g_x_yzz_yy_yy, g_x_yzz_yy_yz, g_x_yzz_yy_zz, g_x_z_0_0_0_yz_yy_xx, g_x_z_0_0_0_yz_yy_xy, g_x_z_0_0_0_yz_yy_xz, g_x_z_0_0_0_yz_yy_yy, g_x_z_0_0_0_yz_yy_yz, g_x_z_0_0_0_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yz_yy_xx[i] = -2.0 * g_x_y_yy_xx[i] * a_exp + 4.0 * g_x_yzz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_yy_xy[i] = -2.0 * g_x_y_yy_xy[i] * a_exp + 4.0 * g_x_yzz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_yy_xz[i] = -2.0 * g_x_y_yy_xz[i] * a_exp + 4.0 * g_x_yzz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_yy_yy[i] = -2.0 * g_x_y_yy_yy[i] * a_exp + 4.0 * g_x_yzz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_yy_yz[i] = -2.0 * g_x_y_yy_yz[i] * a_exp + 4.0 * g_x_yzz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_yy_zz[i] = -2.0 * g_x_y_yy_zz[i] * a_exp + 4.0 * g_x_yzz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_x_y_yz_xx, g_x_y_yz_xy, g_x_y_yz_xz, g_x_y_yz_yy, g_x_y_yz_yz, g_x_y_yz_zz, g_x_yzz_yz_xx, g_x_yzz_yz_xy, g_x_yzz_yz_xz, g_x_yzz_yz_yy, g_x_yzz_yz_yz, g_x_yzz_yz_zz, g_x_z_0_0_0_yz_yz_xx, g_x_z_0_0_0_yz_yz_xy, g_x_z_0_0_0_yz_yz_xz, g_x_z_0_0_0_yz_yz_yy, g_x_z_0_0_0_yz_yz_yz, g_x_z_0_0_0_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yz_yz_xx[i] = -2.0 * g_x_y_yz_xx[i] * a_exp + 4.0 * g_x_yzz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_yz_xy[i] = -2.0 * g_x_y_yz_xy[i] * a_exp + 4.0 * g_x_yzz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_yz_xz[i] = -2.0 * g_x_y_yz_xz[i] * a_exp + 4.0 * g_x_yzz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_yz_yy[i] = -2.0 * g_x_y_yz_yy[i] * a_exp + 4.0 * g_x_yzz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_yz_yz[i] = -2.0 * g_x_y_yz_yz[i] * a_exp + 4.0 * g_x_yzz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_yz_zz[i] = -2.0 * g_x_y_yz_zz[i] * a_exp + 4.0 * g_x_yzz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_x_y_zz_xx, g_x_y_zz_xy, g_x_y_zz_xz, g_x_y_zz_yy, g_x_y_zz_yz, g_x_y_zz_zz, g_x_yzz_zz_xx, g_x_yzz_zz_xy, g_x_yzz_zz_xz, g_x_yzz_zz_yy, g_x_yzz_zz_yz, g_x_yzz_zz_zz, g_x_z_0_0_0_yz_zz_xx, g_x_z_0_0_0_yz_zz_xy, g_x_z_0_0_0_yz_zz_xz, g_x_z_0_0_0_yz_zz_yy, g_x_z_0_0_0_yz_zz_yz, g_x_z_0_0_0_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yz_zz_xx[i] = -2.0 * g_x_y_zz_xx[i] * a_exp + 4.0 * g_x_yzz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_zz_xy[i] = -2.0 * g_x_y_zz_xy[i] * a_exp + 4.0 * g_x_yzz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_zz_xz[i] = -2.0 * g_x_y_zz_xz[i] * a_exp + 4.0 * g_x_yzz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_zz_yy[i] = -2.0 * g_x_y_zz_yy[i] * a_exp + 4.0 * g_x_yzz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_zz_yz[i] = -2.0 * g_x_y_zz_yz[i] * a_exp + 4.0 * g_x_yzz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_zz_zz[i] = -2.0 * g_x_y_zz_zz[i] * a_exp + 4.0 * g_x_yzz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_x_z_0_0_0_zz_xx_xx, g_x_z_0_0_0_zz_xx_xy, g_x_z_0_0_0_zz_xx_xz, g_x_z_0_0_0_zz_xx_yy, g_x_z_0_0_0_zz_xx_yz, g_x_z_0_0_0_zz_xx_zz, g_x_z_xx_xx, g_x_z_xx_xy, g_x_z_xx_xz, g_x_z_xx_yy, g_x_z_xx_yz, g_x_z_xx_zz, g_x_zzz_xx_xx, g_x_zzz_xx_xy, g_x_zzz_xx_xz, g_x_zzz_xx_yy, g_x_zzz_xx_yz, g_x_zzz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_zz_xx_xx[i] = -4.0 * g_x_z_xx_xx[i] * a_exp + 4.0 * g_x_zzz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_xx_xy[i] = -4.0 * g_x_z_xx_xy[i] * a_exp + 4.0 * g_x_zzz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_xx_xz[i] = -4.0 * g_x_z_xx_xz[i] * a_exp + 4.0 * g_x_zzz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_xx_yy[i] = -4.0 * g_x_z_xx_yy[i] * a_exp + 4.0 * g_x_zzz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_xx_yz[i] = -4.0 * g_x_z_xx_yz[i] * a_exp + 4.0 * g_x_zzz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_xx_zz[i] = -4.0 * g_x_z_xx_zz[i] * a_exp + 4.0 * g_x_zzz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_x_z_0_0_0_zz_xy_xx, g_x_z_0_0_0_zz_xy_xy, g_x_z_0_0_0_zz_xy_xz, g_x_z_0_0_0_zz_xy_yy, g_x_z_0_0_0_zz_xy_yz, g_x_z_0_0_0_zz_xy_zz, g_x_z_xy_xx, g_x_z_xy_xy, g_x_z_xy_xz, g_x_z_xy_yy, g_x_z_xy_yz, g_x_z_xy_zz, g_x_zzz_xy_xx, g_x_zzz_xy_xy, g_x_zzz_xy_xz, g_x_zzz_xy_yy, g_x_zzz_xy_yz, g_x_zzz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_zz_xy_xx[i] = -4.0 * g_x_z_xy_xx[i] * a_exp + 4.0 * g_x_zzz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_xy_xy[i] = -4.0 * g_x_z_xy_xy[i] * a_exp + 4.0 * g_x_zzz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_xy_xz[i] = -4.0 * g_x_z_xy_xz[i] * a_exp + 4.0 * g_x_zzz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_xy_yy[i] = -4.0 * g_x_z_xy_yy[i] * a_exp + 4.0 * g_x_zzz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_xy_yz[i] = -4.0 * g_x_z_xy_yz[i] * a_exp + 4.0 * g_x_zzz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_xy_zz[i] = -4.0 * g_x_z_xy_zz[i] * a_exp + 4.0 * g_x_zzz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_x_z_0_0_0_zz_xz_xx, g_x_z_0_0_0_zz_xz_xy, g_x_z_0_0_0_zz_xz_xz, g_x_z_0_0_0_zz_xz_yy, g_x_z_0_0_0_zz_xz_yz, g_x_z_0_0_0_zz_xz_zz, g_x_z_xz_xx, g_x_z_xz_xy, g_x_z_xz_xz, g_x_z_xz_yy, g_x_z_xz_yz, g_x_z_xz_zz, g_x_zzz_xz_xx, g_x_zzz_xz_xy, g_x_zzz_xz_xz, g_x_zzz_xz_yy, g_x_zzz_xz_yz, g_x_zzz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_zz_xz_xx[i] = -4.0 * g_x_z_xz_xx[i] * a_exp + 4.0 * g_x_zzz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_xz_xy[i] = -4.0 * g_x_z_xz_xy[i] * a_exp + 4.0 * g_x_zzz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_xz_xz[i] = -4.0 * g_x_z_xz_xz[i] * a_exp + 4.0 * g_x_zzz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_xz_yy[i] = -4.0 * g_x_z_xz_yy[i] * a_exp + 4.0 * g_x_zzz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_xz_yz[i] = -4.0 * g_x_z_xz_yz[i] * a_exp + 4.0 * g_x_zzz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_xz_zz[i] = -4.0 * g_x_z_xz_zz[i] * a_exp + 4.0 * g_x_zzz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_x_z_0_0_0_zz_yy_xx, g_x_z_0_0_0_zz_yy_xy, g_x_z_0_0_0_zz_yy_xz, g_x_z_0_0_0_zz_yy_yy, g_x_z_0_0_0_zz_yy_yz, g_x_z_0_0_0_zz_yy_zz, g_x_z_yy_xx, g_x_z_yy_xy, g_x_z_yy_xz, g_x_z_yy_yy, g_x_z_yy_yz, g_x_z_yy_zz, g_x_zzz_yy_xx, g_x_zzz_yy_xy, g_x_zzz_yy_xz, g_x_zzz_yy_yy, g_x_zzz_yy_yz, g_x_zzz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_zz_yy_xx[i] = -4.0 * g_x_z_yy_xx[i] * a_exp + 4.0 * g_x_zzz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_yy_xy[i] = -4.0 * g_x_z_yy_xy[i] * a_exp + 4.0 * g_x_zzz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_yy_xz[i] = -4.0 * g_x_z_yy_xz[i] * a_exp + 4.0 * g_x_zzz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_yy_yy[i] = -4.0 * g_x_z_yy_yy[i] * a_exp + 4.0 * g_x_zzz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_yy_yz[i] = -4.0 * g_x_z_yy_yz[i] * a_exp + 4.0 * g_x_zzz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_yy_zz[i] = -4.0 * g_x_z_yy_zz[i] * a_exp + 4.0 * g_x_zzz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_x_z_0_0_0_zz_yz_xx, g_x_z_0_0_0_zz_yz_xy, g_x_z_0_0_0_zz_yz_xz, g_x_z_0_0_0_zz_yz_yy, g_x_z_0_0_0_zz_yz_yz, g_x_z_0_0_0_zz_yz_zz, g_x_z_yz_xx, g_x_z_yz_xy, g_x_z_yz_xz, g_x_z_yz_yy, g_x_z_yz_yz, g_x_z_yz_zz, g_x_zzz_yz_xx, g_x_zzz_yz_xy, g_x_zzz_yz_xz, g_x_zzz_yz_yy, g_x_zzz_yz_yz, g_x_zzz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_zz_yz_xx[i] = -4.0 * g_x_z_yz_xx[i] * a_exp + 4.0 * g_x_zzz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_yz_xy[i] = -4.0 * g_x_z_yz_xy[i] * a_exp + 4.0 * g_x_zzz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_yz_xz[i] = -4.0 * g_x_z_yz_xz[i] * a_exp + 4.0 * g_x_zzz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_yz_yy[i] = -4.0 * g_x_z_yz_yy[i] * a_exp + 4.0 * g_x_zzz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_yz_yz[i] = -4.0 * g_x_z_yz_yz[i] * a_exp + 4.0 * g_x_zzz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_yz_zz[i] = -4.0 * g_x_z_yz_zz[i] * a_exp + 4.0 * g_x_zzz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_x_z_0_0_0_zz_zz_xx, g_x_z_0_0_0_zz_zz_xy, g_x_z_0_0_0_zz_zz_xz, g_x_z_0_0_0_zz_zz_yy, g_x_z_0_0_0_zz_zz_yz, g_x_z_0_0_0_zz_zz_zz, g_x_z_zz_xx, g_x_z_zz_xy, g_x_z_zz_xz, g_x_z_zz_yy, g_x_z_zz_yz, g_x_z_zz_zz, g_x_zzz_zz_xx, g_x_zzz_zz_xy, g_x_zzz_zz_xz, g_x_zzz_zz_yy, g_x_zzz_zz_yz, g_x_zzz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_zz_zz_xx[i] = -4.0 * g_x_z_zz_xx[i] * a_exp + 4.0 * g_x_zzz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_zz_xy[i] = -4.0 * g_x_z_zz_xy[i] * a_exp + 4.0 * g_x_zzz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_zz_xz[i] = -4.0 * g_x_z_zz_xz[i] * a_exp + 4.0 * g_x_zzz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_zz_yy[i] = -4.0 * g_x_z_zz_yy[i] * a_exp + 4.0 * g_x_zzz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_zz_yz[i] = -4.0 * g_x_z_zz_yz[i] * a_exp + 4.0 * g_x_zzz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_zz_zz[i] = -4.0 * g_x_z_zz_zz[i] * a_exp + 4.0 * g_x_zzz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_y_x_0_0_0_xx_xx_xx, g_y_x_0_0_0_xx_xx_xy, g_y_x_0_0_0_xx_xx_xz, g_y_x_0_0_0_xx_xx_yy, g_y_x_0_0_0_xx_xx_yz, g_y_x_0_0_0_xx_xx_zz, g_y_x_xx_xx, g_y_x_xx_xy, g_y_x_xx_xz, g_y_x_xx_yy, g_y_x_xx_yz, g_y_x_xx_zz, g_y_xxx_xx_xx, g_y_xxx_xx_xy, g_y_xxx_xx_xz, g_y_xxx_xx_yy, g_y_xxx_xx_yz, g_y_xxx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xx_xx_xx[i] = -4.0 * g_y_x_xx_xx[i] * a_exp + 4.0 * g_y_xxx_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_xx_xy[i] = -4.0 * g_y_x_xx_xy[i] * a_exp + 4.0 * g_y_xxx_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_xx_xz[i] = -4.0 * g_y_x_xx_xz[i] * a_exp + 4.0 * g_y_xxx_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_xx_yy[i] = -4.0 * g_y_x_xx_yy[i] * a_exp + 4.0 * g_y_xxx_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_xx_yz[i] = -4.0 * g_y_x_xx_yz[i] * a_exp + 4.0 * g_y_xxx_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_xx_zz[i] = -4.0 * g_y_x_xx_zz[i] * a_exp + 4.0 * g_y_xxx_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_y_x_0_0_0_xx_xy_xx, g_y_x_0_0_0_xx_xy_xy, g_y_x_0_0_0_xx_xy_xz, g_y_x_0_0_0_xx_xy_yy, g_y_x_0_0_0_xx_xy_yz, g_y_x_0_0_0_xx_xy_zz, g_y_x_xy_xx, g_y_x_xy_xy, g_y_x_xy_xz, g_y_x_xy_yy, g_y_x_xy_yz, g_y_x_xy_zz, g_y_xxx_xy_xx, g_y_xxx_xy_xy, g_y_xxx_xy_xz, g_y_xxx_xy_yy, g_y_xxx_xy_yz, g_y_xxx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xx_xy_xx[i] = -4.0 * g_y_x_xy_xx[i] * a_exp + 4.0 * g_y_xxx_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_xy_xy[i] = -4.0 * g_y_x_xy_xy[i] * a_exp + 4.0 * g_y_xxx_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_xy_xz[i] = -4.0 * g_y_x_xy_xz[i] * a_exp + 4.0 * g_y_xxx_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_xy_yy[i] = -4.0 * g_y_x_xy_yy[i] * a_exp + 4.0 * g_y_xxx_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_xy_yz[i] = -4.0 * g_y_x_xy_yz[i] * a_exp + 4.0 * g_y_xxx_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_xy_zz[i] = -4.0 * g_y_x_xy_zz[i] * a_exp + 4.0 * g_y_xxx_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_y_x_0_0_0_xx_xz_xx, g_y_x_0_0_0_xx_xz_xy, g_y_x_0_0_0_xx_xz_xz, g_y_x_0_0_0_xx_xz_yy, g_y_x_0_0_0_xx_xz_yz, g_y_x_0_0_0_xx_xz_zz, g_y_x_xz_xx, g_y_x_xz_xy, g_y_x_xz_xz, g_y_x_xz_yy, g_y_x_xz_yz, g_y_x_xz_zz, g_y_xxx_xz_xx, g_y_xxx_xz_xy, g_y_xxx_xz_xz, g_y_xxx_xz_yy, g_y_xxx_xz_yz, g_y_xxx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xx_xz_xx[i] = -4.0 * g_y_x_xz_xx[i] * a_exp + 4.0 * g_y_xxx_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_xz_xy[i] = -4.0 * g_y_x_xz_xy[i] * a_exp + 4.0 * g_y_xxx_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_xz_xz[i] = -4.0 * g_y_x_xz_xz[i] * a_exp + 4.0 * g_y_xxx_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_xz_yy[i] = -4.0 * g_y_x_xz_yy[i] * a_exp + 4.0 * g_y_xxx_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_xz_yz[i] = -4.0 * g_y_x_xz_yz[i] * a_exp + 4.0 * g_y_xxx_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_xz_zz[i] = -4.0 * g_y_x_xz_zz[i] * a_exp + 4.0 * g_y_xxx_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_y_x_0_0_0_xx_yy_xx, g_y_x_0_0_0_xx_yy_xy, g_y_x_0_0_0_xx_yy_xz, g_y_x_0_0_0_xx_yy_yy, g_y_x_0_0_0_xx_yy_yz, g_y_x_0_0_0_xx_yy_zz, g_y_x_yy_xx, g_y_x_yy_xy, g_y_x_yy_xz, g_y_x_yy_yy, g_y_x_yy_yz, g_y_x_yy_zz, g_y_xxx_yy_xx, g_y_xxx_yy_xy, g_y_xxx_yy_xz, g_y_xxx_yy_yy, g_y_xxx_yy_yz, g_y_xxx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xx_yy_xx[i] = -4.0 * g_y_x_yy_xx[i] * a_exp + 4.0 * g_y_xxx_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_yy_xy[i] = -4.0 * g_y_x_yy_xy[i] * a_exp + 4.0 * g_y_xxx_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_yy_xz[i] = -4.0 * g_y_x_yy_xz[i] * a_exp + 4.0 * g_y_xxx_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_yy_yy[i] = -4.0 * g_y_x_yy_yy[i] * a_exp + 4.0 * g_y_xxx_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_yy_yz[i] = -4.0 * g_y_x_yy_yz[i] * a_exp + 4.0 * g_y_xxx_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_yy_zz[i] = -4.0 * g_y_x_yy_zz[i] * a_exp + 4.0 * g_y_xxx_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_y_x_0_0_0_xx_yz_xx, g_y_x_0_0_0_xx_yz_xy, g_y_x_0_0_0_xx_yz_xz, g_y_x_0_0_0_xx_yz_yy, g_y_x_0_0_0_xx_yz_yz, g_y_x_0_0_0_xx_yz_zz, g_y_x_yz_xx, g_y_x_yz_xy, g_y_x_yz_xz, g_y_x_yz_yy, g_y_x_yz_yz, g_y_x_yz_zz, g_y_xxx_yz_xx, g_y_xxx_yz_xy, g_y_xxx_yz_xz, g_y_xxx_yz_yy, g_y_xxx_yz_yz, g_y_xxx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xx_yz_xx[i] = -4.0 * g_y_x_yz_xx[i] * a_exp + 4.0 * g_y_xxx_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_yz_xy[i] = -4.0 * g_y_x_yz_xy[i] * a_exp + 4.0 * g_y_xxx_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_yz_xz[i] = -4.0 * g_y_x_yz_xz[i] * a_exp + 4.0 * g_y_xxx_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_yz_yy[i] = -4.0 * g_y_x_yz_yy[i] * a_exp + 4.0 * g_y_xxx_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_yz_yz[i] = -4.0 * g_y_x_yz_yz[i] * a_exp + 4.0 * g_y_xxx_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_yz_zz[i] = -4.0 * g_y_x_yz_zz[i] * a_exp + 4.0 * g_y_xxx_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_y_x_0_0_0_xx_zz_xx, g_y_x_0_0_0_xx_zz_xy, g_y_x_0_0_0_xx_zz_xz, g_y_x_0_0_0_xx_zz_yy, g_y_x_0_0_0_xx_zz_yz, g_y_x_0_0_0_xx_zz_zz, g_y_x_zz_xx, g_y_x_zz_xy, g_y_x_zz_xz, g_y_x_zz_yy, g_y_x_zz_yz, g_y_x_zz_zz, g_y_xxx_zz_xx, g_y_xxx_zz_xy, g_y_xxx_zz_xz, g_y_xxx_zz_yy, g_y_xxx_zz_yz, g_y_xxx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xx_zz_xx[i] = -4.0 * g_y_x_zz_xx[i] * a_exp + 4.0 * g_y_xxx_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_zz_xy[i] = -4.0 * g_y_x_zz_xy[i] * a_exp + 4.0 * g_y_xxx_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_zz_xz[i] = -4.0 * g_y_x_zz_xz[i] * a_exp + 4.0 * g_y_xxx_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_zz_yy[i] = -4.0 * g_y_x_zz_yy[i] * a_exp + 4.0 * g_y_xxx_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_zz_yz[i] = -4.0 * g_y_x_zz_yz[i] * a_exp + 4.0 * g_y_xxx_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_zz_zz[i] = -4.0 * g_y_x_zz_zz[i] * a_exp + 4.0 * g_y_xxx_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_y_x_0_0_0_xy_xx_xx, g_y_x_0_0_0_xy_xx_xy, g_y_x_0_0_0_xy_xx_xz, g_y_x_0_0_0_xy_xx_yy, g_y_x_0_0_0_xy_xx_yz, g_y_x_0_0_0_xy_xx_zz, g_y_xxy_xx_xx, g_y_xxy_xx_xy, g_y_xxy_xx_xz, g_y_xxy_xx_yy, g_y_xxy_xx_yz, g_y_xxy_xx_zz, g_y_y_xx_xx, g_y_y_xx_xy, g_y_y_xx_xz, g_y_y_xx_yy, g_y_y_xx_yz, g_y_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xy_xx_xx[i] = -2.0 * g_y_y_xx_xx[i] * a_exp + 4.0 * g_y_xxy_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_xx_xy[i] = -2.0 * g_y_y_xx_xy[i] * a_exp + 4.0 * g_y_xxy_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_xx_xz[i] = -2.0 * g_y_y_xx_xz[i] * a_exp + 4.0 * g_y_xxy_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_xx_yy[i] = -2.0 * g_y_y_xx_yy[i] * a_exp + 4.0 * g_y_xxy_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_xx_yz[i] = -2.0 * g_y_y_xx_yz[i] * a_exp + 4.0 * g_y_xxy_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_xx_zz[i] = -2.0 * g_y_y_xx_zz[i] * a_exp + 4.0 * g_y_xxy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_y_x_0_0_0_xy_xy_xx, g_y_x_0_0_0_xy_xy_xy, g_y_x_0_0_0_xy_xy_xz, g_y_x_0_0_0_xy_xy_yy, g_y_x_0_0_0_xy_xy_yz, g_y_x_0_0_0_xy_xy_zz, g_y_xxy_xy_xx, g_y_xxy_xy_xy, g_y_xxy_xy_xz, g_y_xxy_xy_yy, g_y_xxy_xy_yz, g_y_xxy_xy_zz, g_y_y_xy_xx, g_y_y_xy_xy, g_y_y_xy_xz, g_y_y_xy_yy, g_y_y_xy_yz, g_y_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xy_xy_xx[i] = -2.0 * g_y_y_xy_xx[i] * a_exp + 4.0 * g_y_xxy_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_xy_xy[i] = -2.0 * g_y_y_xy_xy[i] * a_exp + 4.0 * g_y_xxy_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_xy_xz[i] = -2.0 * g_y_y_xy_xz[i] * a_exp + 4.0 * g_y_xxy_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_xy_yy[i] = -2.0 * g_y_y_xy_yy[i] * a_exp + 4.0 * g_y_xxy_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_xy_yz[i] = -2.0 * g_y_y_xy_yz[i] * a_exp + 4.0 * g_y_xxy_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_xy_zz[i] = -2.0 * g_y_y_xy_zz[i] * a_exp + 4.0 * g_y_xxy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_y_x_0_0_0_xy_xz_xx, g_y_x_0_0_0_xy_xz_xy, g_y_x_0_0_0_xy_xz_xz, g_y_x_0_0_0_xy_xz_yy, g_y_x_0_0_0_xy_xz_yz, g_y_x_0_0_0_xy_xz_zz, g_y_xxy_xz_xx, g_y_xxy_xz_xy, g_y_xxy_xz_xz, g_y_xxy_xz_yy, g_y_xxy_xz_yz, g_y_xxy_xz_zz, g_y_y_xz_xx, g_y_y_xz_xy, g_y_y_xz_xz, g_y_y_xz_yy, g_y_y_xz_yz, g_y_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xy_xz_xx[i] = -2.0 * g_y_y_xz_xx[i] * a_exp + 4.0 * g_y_xxy_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_xz_xy[i] = -2.0 * g_y_y_xz_xy[i] * a_exp + 4.0 * g_y_xxy_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_xz_xz[i] = -2.0 * g_y_y_xz_xz[i] * a_exp + 4.0 * g_y_xxy_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_xz_yy[i] = -2.0 * g_y_y_xz_yy[i] * a_exp + 4.0 * g_y_xxy_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_xz_yz[i] = -2.0 * g_y_y_xz_yz[i] * a_exp + 4.0 * g_y_xxy_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_xz_zz[i] = -2.0 * g_y_y_xz_zz[i] * a_exp + 4.0 * g_y_xxy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_y_x_0_0_0_xy_yy_xx, g_y_x_0_0_0_xy_yy_xy, g_y_x_0_0_0_xy_yy_xz, g_y_x_0_0_0_xy_yy_yy, g_y_x_0_0_0_xy_yy_yz, g_y_x_0_0_0_xy_yy_zz, g_y_xxy_yy_xx, g_y_xxy_yy_xy, g_y_xxy_yy_xz, g_y_xxy_yy_yy, g_y_xxy_yy_yz, g_y_xxy_yy_zz, g_y_y_yy_xx, g_y_y_yy_xy, g_y_y_yy_xz, g_y_y_yy_yy, g_y_y_yy_yz, g_y_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xy_yy_xx[i] = -2.0 * g_y_y_yy_xx[i] * a_exp + 4.0 * g_y_xxy_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_yy_xy[i] = -2.0 * g_y_y_yy_xy[i] * a_exp + 4.0 * g_y_xxy_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_yy_xz[i] = -2.0 * g_y_y_yy_xz[i] * a_exp + 4.0 * g_y_xxy_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_yy_yy[i] = -2.0 * g_y_y_yy_yy[i] * a_exp + 4.0 * g_y_xxy_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_yy_yz[i] = -2.0 * g_y_y_yy_yz[i] * a_exp + 4.0 * g_y_xxy_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_yy_zz[i] = -2.0 * g_y_y_yy_zz[i] * a_exp + 4.0 * g_y_xxy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_y_x_0_0_0_xy_yz_xx, g_y_x_0_0_0_xy_yz_xy, g_y_x_0_0_0_xy_yz_xz, g_y_x_0_0_0_xy_yz_yy, g_y_x_0_0_0_xy_yz_yz, g_y_x_0_0_0_xy_yz_zz, g_y_xxy_yz_xx, g_y_xxy_yz_xy, g_y_xxy_yz_xz, g_y_xxy_yz_yy, g_y_xxy_yz_yz, g_y_xxy_yz_zz, g_y_y_yz_xx, g_y_y_yz_xy, g_y_y_yz_xz, g_y_y_yz_yy, g_y_y_yz_yz, g_y_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xy_yz_xx[i] = -2.0 * g_y_y_yz_xx[i] * a_exp + 4.0 * g_y_xxy_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_yz_xy[i] = -2.0 * g_y_y_yz_xy[i] * a_exp + 4.0 * g_y_xxy_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_yz_xz[i] = -2.0 * g_y_y_yz_xz[i] * a_exp + 4.0 * g_y_xxy_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_yz_yy[i] = -2.0 * g_y_y_yz_yy[i] * a_exp + 4.0 * g_y_xxy_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_yz_yz[i] = -2.0 * g_y_y_yz_yz[i] * a_exp + 4.0 * g_y_xxy_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_yz_zz[i] = -2.0 * g_y_y_yz_zz[i] * a_exp + 4.0 * g_y_xxy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_y_x_0_0_0_xy_zz_xx, g_y_x_0_0_0_xy_zz_xy, g_y_x_0_0_0_xy_zz_xz, g_y_x_0_0_0_xy_zz_yy, g_y_x_0_0_0_xy_zz_yz, g_y_x_0_0_0_xy_zz_zz, g_y_xxy_zz_xx, g_y_xxy_zz_xy, g_y_xxy_zz_xz, g_y_xxy_zz_yy, g_y_xxy_zz_yz, g_y_xxy_zz_zz, g_y_y_zz_xx, g_y_y_zz_xy, g_y_y_zz_xz, g_y_y_zz_yy, g_y_y_zz_yz, g_y_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xy_zz_xx[i] = -2.0 * g_y_y_zz_xx[i] * a_exp + 4.0 * g_y_xxy_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_zz_xy[i] = -2.0 * g_y_y_zz_xy[i] * a_exp + 4.0 * g_y_xxy_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_zz_xz[i] = -2.0 * g_y_y_zz_xz[i] * a_exp + 4.0 * g_y_xxy_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_zz_yy[i] = -2.0 * g_y_y_zz_yy[i] * a_exp + 4.0 * g_y_xxy_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_zz_yz[i] = -2.0 * g_y_y_zz_yz[i] * a_exp + 4.0 * g_y_xxy_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_zz_zz[i] = -2.0 * g_y_y_zz_zz[i] * a_exp + 4.0 * g_y_xxy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_y_x_0_0_0_xz_xx_xx, g_y_x_0_0_0_xz_xx_xy, g_y_x_0_0_0_xz_xx_xz, g_y_x_0_0_0_xz_xx_yy, g_y_x_0_0_0_xz_xx_yz, g_y_x_0_0_0_xz_xx_zz, g_y_xxz_xx_xx, g_y_xxz_xx_xy, g_y_xxz_xx_xz, g_y_xxz_xx_yy, g_y_xxz_xx_yz, g_y_xxz_xx_zz, g_y_z_xx_xx, g_y_z_xx_xy, g_y_z_xx_xz, g_y_z_xx_yy, g_y_z_xx_yz, g_y_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xz_xx_xx[i] = -2.0 * g_y_z_xx_xx[i] * a_exp + 4.0 * g_y_xxz_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_xx_xy[i] = -2.0 * g_y_z_xx_xy[i] * a_exp + 4.0 * g_y_xxz_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_xx_xz[i] = -2.0 * g_y_z_xx_xz[i] * a_exp + 4.0 * g_y_xxz_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_xx_yy[i] = -2.0 * g_y_z_xx_yy[i] * a_exp + 4.0 * g_y_xxz_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_xx_yz[i] = -2.0 * g_y_z_xx_yz[i] * a_exp + 4.0 * g_y_xxz_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_xx_zz[i] = -2.0 * g_y_z_xx_zz[i] * a_exp + 4.0 * g_y_xxz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_y_x_0_0_0_xz_xy_xx, g_y_x_0_0_0_xz_xy_xy, g_y_x_0_0_0_xz_xy_xz, g_y_x_0_0_0_xz_xy_yy, g_y_x_0_0_0_xz_xy_yz, g_y_x_0_0_0_xz_xy_zz, g_y_xxz_xy_xx, g_y_xxz_xy_xy, g_y_xxz_xy_xz, g_y_xxz_xy_yy, g_y_xxz_xy_yz, g_y_xxz_xy_zz, g_y_z_xy_xx, g_y_z_xy_xy, g_y_z_xy_xz, g_y_z_xy_yy, g_y_z_xy_yz, g_y_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xz_xy_xx[i] = -2.0 * g_y_z_xy_xx[i] * a_exp + 4.0 * g_y_xxz_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_xy_xy[i] = -2.0 * g_y_z_xy_xy[i] * a_exp + 4.0 * g_y_xxz_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_xy_xz[i] = -2.0 * g_y_z_xy_xz[i] * a_exp + 4.0 * g_y_xxz_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_xy_yy[i] = -2.0 * g_y_z_xy_yy[i] * a_exp + 4.0 * g_y_xxz_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_xy_yz[i] = -2.0 * g_y_z_xy_yz[i] * a_exp + 4.0 * g_y_xxz_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_xy_zz[i] = -2.0 * g_y_z_xy_zz[i] * a_exp + 4.0 * g_y_xxz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_y_x_0_0_0_xz_xz_xx, g_y_x_0_0_0_xz_xz_xy, g_y_x_0_0_0_xz_xz_xz, g_y_x_0_0_0_xz_xz_yy, g_y_x_0_0_0_xz_xz_yz, g_y_x_0_0_0_xz_xz_zz, g_y_xxz_xz_xx, g_y_xxz_xz_xy, g_y_xxz_xz_xz, g_y_xxz_xz_yy, g_y_xxz_xz_yz, g_y_xxz_xz_zz, g_y_z_xz_xx, g_y_z_xz_xy, g_y_z_xz_xz, g_y_z_xz_yy, g_y_z_xz_yz, g_y_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xz_xz_xx[i] = -2.0 * g_y_z_xz_xx[i] * a_exp + 4.0 * g_y_xxz_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_xz_xy[i] = -2.0 * g_y_z_xz_xy[i] * a_exp + 4.0 * g_y_xxz_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_xz_xz[i] = -2.0 * g_y_z_xz_xz[i] * a_exp + 4.0 * g_y_xxz_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_xz_yy[i] = -2.0 * g_y_z_xz_yy[i] * a_exp + 4.0 * g_y_xxz_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_xz_yz[i] = -2.0 * g_y_z_xz_yz[i] * a_exp + 4.0 * g_y_xxz_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_xz_zz[i] = -2.0 * g_y_z_xz_zz[i] * a_exp + 4.0 * g_y_xxz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_y_x_0_0_0_xz_yy_xx, g_y_x_0_0_0_xz_yy_xy, g_y_x_0_0_0_xz_yy_xz, g_y_x_0_0_0_xz_yy_yy, g_y_x_0_0_0_xz_yy_yz, g_y_x_0_0_0_xz_yy_zz, g_y_xxz_yy_xx, g_y_xxz_yy_xy, g_y_xxz_yy_xz, g_y_xxz_yy_yy, g_y_xxz_yy_yz, g_y_xxz_yy_zz, g_y_z_yy_xx, g_y_z_yy_xy, g_y_z_yy_xz, g_y_z_yy_yy, g_y_z_yy_yz, g_y_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xz_yy_xx[i] = -2.0 * g_y_z_yy_xx[i] * a_exp + 4.0 * g_y_xxz_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_yy_xy[i] = -2.0 * g_y_z_yy_xy[i] * a_exp + 4.0 * g_y_xxz_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_yy_xz[i] = -2.0 * g_y_z_yy_xz[i] * a_exp + 4.0 * g_y_xxz_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_yy_yy[i] = -2.0 * g_y_z_yy_yy[i] * a_exp + 4.0 * g_y_xxz_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_yy_yz[i] = -2.0 * g_y_z_yy_yz[i] * a_exp + 4.0 * g_y_xxz_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_yy_zz[i] = -2.0 * g_y_z_yy_zz[i] * a_exp + 4.0 * g_y_xxz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_y_x_0_0_0_xz_yz_xx, g_y_x_0_0_0_xz_yz_xy, g_y_x_0_0_0_xz_yz_xz, g_y_x_0_0_0_xz_yz_yy, g_y_x_0_0_0_xz_yz_yz, g_y_x_0_0_0_xz_yz_zz, g_y_xxz_yz_xx, g_y_xxz_yz_xy, g_y_xxz_yz_xz, g_y_xxz_yz_yy, g_y_xxz_yz_yz, g_y_xxz_yz_zz, g_y_z_yz_xx, g_y_z_yz_xy, g_y_z_yz_xz, g_y_z_yz_yy, g_y_z_yz_yz, g_y_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xz_yz_xx[i] = -2.0 * g_y_z_yz_xx[i] * a_exp + 4.0 * g_y_xxz_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_yz_xy[i] = -2.0 * g_y_z_yz_xy[i] * a_exp + 4.0 * g_y_xxz_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_yz_xz[i] = -2.0 * g_y_z_yz_xz[i] * a_exp + 4.0 * g_y_xxz_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_yz_yy[i] = -2.0 * g_y_z_yz_yy[i] * a_exp + 4.0 * g_y_xxz_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_yz_yz[i] = -2.0 * g_y_z_yz_yz[i] * a_exp + 4.0 * g_y_xxz_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_yz_zz[i] = -2.0 * g_y_z_yz_zz[i] * a_exp + 4.0 * g_y_xxz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_y_x_0_0_0_xz_zz_xx, g_y_x_0_0_0_xz_zz_xy, g_y_x_0_0_0_xz_zz_xz, g_y_x_0_0_0_xz_zz_yy, g_y_x_0_0_0_xz_zz_yz, g_y_x_0_0_0_xz_zz_zz, g_y_xxz_zz_xx, g_y_xxz_zz_xy, g_y_xxz_zz_xz, g_y_xxz_zz_yy, g_y_xxz_zz_yz, g_y_xxz_zz_zz, g_y_z_zz_xx, g_y_z_zz_xy, g_y_z_zz_xz, g_y_z_zz_yy, g_y_z_zz_yz, g_y_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xz_zz_xx[i] = -2.0 * g_y_z_zz_xx[i] * a_exp + 4.0 * g_y_xxz_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_zz_xy[i] = -2.0 * g_y_z_zz_xy[i] * a_exp + 4.0 * g_y_xxz_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_zz_xz[i] = -2.0 * g_y_z_zz_xz[i] * a_exp + 4.0 * g_y_xxz_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_zz_yy[i] = -2.0 * g_y_z_zz_yy[i] * a_exp + 4.0 * g_y_xxz_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_zz_yz[i] = -2.0 * g_y_z_zz_yz[i] * a_exp + 4.0 * g_y_xxz_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_zz_zz[i] = -2.0 * g_y_z_zz_zz[i] * a_exp + 4.0 * g_y_xxz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_y_x_0_0_0_yy_xx_xx, g_y_x_0_0_0_yy_xx_xy, g_y_x_0_0_0_yy_xx_xz, g_y_x_0_0_0_yy_xx_yy, g_y_x_0_0_0_yy_xx_yz, g_y_x_0_0_0_yy_xx_zz, g_y_xyy_xx_xx, g_y_xyy_xx_xy, g_y_xyy_xx_xz, g_y_xyy_xx_yy, g_y_xyy_xx_yz, g_y_xyy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yy_xx_xx[i] = 4.0 * g_y_xyy_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_xx_xy[i] = 4.0 * g_y_xyy_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_xx_xz[i] = 4.0 * g_y_xyy_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_xx_yy[i] = 4.0 * g_y_xyy_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_xx_yz[i] = 4.0 * g_y_xyy_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_xx_zz[i] = 4.0 * g_y_xyy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_y_x_0_0_0_yy_xy_xx, g_y_x_0_0_0_yy_xy_xy, g_y_x_0_0_0_yy_xy_xz, g_y_x_0_0_0_yy_xy_yy, g_y_x_0_0_0_yy_xy_yz, g_y_x_0_0_0_yy_xy_zz, g_y_xyy_xy_xx, g_y_xyy_xy_xy, g_y_xyy_xy_xz, g_y_xyy_xy_yy, g_y_xyy_xy_yz, g_y_xyy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yy_xy_xx[i] = 4.0 * g_y_xyy_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_xy_xy[i] = 4.0 * g_y_xyy_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_xy_xz[i] = 4.0 * g_y_xyy_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_xy_yy[i] = 4.0 * g_y_xyy_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_xy_yz[i] = 4.0 * g_y_xyy_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_xy_zz[i] = 4.0 * g_y_xyy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_y_x_0_0_0_yy_xz_xx, g_y_x_0_0_0_yy_xz_xy, g_y_x_0_0_0_yy_xz_xz, g_y_x_0_0_0_yy_xz_yy, g_y_x_0_0_0_yy_xz_yz, g_y_x_0_0_0_yy_xz_zz, g_y_xyy_xz_xx, g_y_xyy_xz_xy, g_y_xyy_xz_xz, g_y_xyy_xz_yy, g_y_xyy_xz_yz, g_y_xyy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yy_xz_xx[i] = 4.0 * g_y_xyy_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_xz_xy[i] = 4.0 * g_y_xyy_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_xz_xz[i] = 4.0 * g_y_xyy_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_xz_yy[i] = 4.0 * g_y_xyy_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_xz_yz[i] = 4.0 * g_y_xyy_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_xz_zz[i] = 4.0 * g_y_xyy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_y_x_0_0_0_yy_yy_xx, g_y_x_0_0_0_yy_yy_xy, g_y_x_0_0_0_yy_yy_xz, g_y_x_0_0_0_yy_yy_yy, g_y_x_0_0_0_yy_yy_yz, g_y_x_0_0_0_yy_yy_zz, g_y_xyy_yy_xx, g_y_xyy_yy_xy, g_y_xyy_yy_xz, g_y_xyy_yy_yy, g_y_xyy_yy_yz, g_y_xyy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yy_yy_xx[i] = 4.0 * g_y_xyy_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_yy_xy[i] = 4.0 * g_y_xyy_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_yy_xz[i] = 4.0 * g_y_xyy_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_yy_yy[i] = 4.0 * g_y_xyy_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_yy_yz[i] = 4.0 * g_y_xyy_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_yy_zz[i] = 4.0 * g_y_xyy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_y_x_0_0_0_yy_yz_xx, g_y_x_0_0_0_yy_yz_xy, g_y_x_0_0_0_yy_yz_xz, g_y_x_0_0_0_yy_yz_yy, g_y_x_0_0_0_yy_yz_yz, g_y_x_0_0_0_yy_yz_zz, g_y_xyy_yz_xx, g_y_xyy_yz_xy, g_y_xyy_yz_xz, g_y_xyy_yz_yy, g_y_xyy_yz_yz, g_y_xyy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yy_yz_xx[i] = 4.0 * g_y_xyy_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_yz_xy[i] = 4.0 * g_y_xyy_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_yz_xz[i] = 4.0 * g_y_xyy_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_yz_yy[i] = 4.0 * g_y_xyy_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_yz_yz[i] = 4.0 * g_y_xyy_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_yz_zz[i] = 4.0 * g_y_xyy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_y_x_0_0_0_yy_zz_xx, g_y_x_0_0_0_yy_zz_xy, g_y_x_0_0_0_yy_zz_xz, g_y_x_0_0_0_yy_zz_yy, g_y_x_0_0_0_yy_zz_yz, g_y_x_0_0_0_yy_zz_zz, g_y_xyy_zz_xx, g_y_xyy_zz_xy, g_y_xyy_zz_xz, g_y_xyy_zz_yy, g_y_xyy_zz_yz, g_y_xyy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yy_zz_xx[i] = 4.0 * g_y_xyy_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_zz_xy[i] = 4.0 * g_y_xyy_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_zz_xz[i] = 4.0 * g_y_xyy_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_zz_yy[i] = 4.0 * g_y_xyy_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_zz_yz[i] = 4.0 * g_y_xyy_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_zz_zz[i] = 4.0 * g_y_xyy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_y_x_0_0_0_yz_xx_xx, g_y_x_0_0_0_yz_xx_xy, g_y_x_0_0_0_yz_xx_xz, g_y_x_0_0_0_yz_xx_yy, g_y_x_0_0_0_yz_xx_yz, g_y_x_0_0_0_yz_xx_zz, g_y_xyz_xx_xx, g_y_xyz_xx_xy, g_y_xyz_xx_xz, g_y_xyz_xx_yy, g_y_xyz_xx_yz, g_y_xyz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yz_xx_xx[i] = 4.0 * g_y_xyz_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_xx_xy[i] = 4.0 * g_y_xyz_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_xx_xz[i] = 4.0 * g_y_xyz_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_xx_yy[i] = 4.0 * g_y_xyz_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_xx_yz[i] = 4.0 * g_y_xyz_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_xx_zz[i] = 4.0 * g_y_xyz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_y_x_0_0_0_yz_xy_xx, g_y_x_0_0_0_yz_xy_xy, g_y_x_0_0_0_yz_xy_xz, g_y_x_0_0_0_yz_xy_yy, g_y_x_0_0_0_yz_xy_yz, g_y_x_0_0_0_yz_xy_zz, g_y_xyz_xy_xx, g_y_xyz_xy_xy, g_y_xyz_xy_xz, g_y_xyz_xy_yy, g_y_xyz_xy_yz, g_y_xyz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yz_xy_xx[i] = 4.0 * g_y_xyz_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_xy_xy[i] = 4.0 * g_y_xyz_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_xy_xz[i] = 4.0 * g_y_xyz_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_xy_yy[i] = 4.0 * g_y_xyz_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_xy_yz[i] = 4.0 * g_y_xyz_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_xy_zz[i] = 4.0 * g_y_xyz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_y_x_0_0_0_yz_xz_xx, g_y_x_0_0_0_yz_xz_xy, g_y_x_0_0_0_yz_xz_xz, g_y_x_0_0_0_yz_xz_yy, g_y_x_0_0_0_yz_xz_yz, g_y_x_0_0_0_yz_xz_zz, g_y_xyz_xz_xx, g_y_xyz_xz_xy, g_y_xyz_xz_xz, g_y_xyz_xz_yy, g_y_xyz_xz_yz, g_y_xyz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yz_xz_xx[i] = 4.0 * g_y_xyz_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_xz_xy[i] = 4.0 * g_y_xyz_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_xz_xz[i] = 4.0 * g_y_xyz_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_xz_yy[i] = 4.0 * g_y_xyz_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_xz_yz[i] = 4.0 * g_y_xyz_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_xz_zz[i] = 4.0 * g_y_xyz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_y_x_0_0_0_yz_yy_xx, g_y_x_0_0_0_yz_yy_xy, g_y_x_0_0_0_yz_yy_xz, g_y_x_0_0_0_yz_yy_yy, g_y_x_0_0_0_yz_yy_yz, g_y_x_0_0_0_yz_yy_zz, g_y_xyz_yy_xx, g_y_xyz_yy_xy, g_y_xyz_yy_xz, g_y_xyz_yy_yy, g_y_xyz_yy_yz, g_y_xyz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yz_yy_xx[i] = 4.0 * g_y_xyz_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_yy_xy[i] = 4.0 * g_y_xyz_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_yy_xz[i] = 4.0 * g_y_xyz_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_yy_yy[i] = 4.0 * g_y_xyz_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_yy_yz[i] = 4.0 * g_y_xyz_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_yy_zz[i] = 4.0 * g_y_xyz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_y_x_0_0_0_yz_yz_xx, g_y_x_0_0_0_yz_yz_xy, g_y_x_0_0_0_yz_yz_xz, g_y_x_0_0_0_yz_yz_yy, g_y_x_0_0_0_yz_yz_yz, g_y_x_0_0_0_yz_yz_zz, g_y_xyz_yz_xx, g_y_xyz_yz_xy, g_y_xyz_yz_xz, g_y_xyz_yz_yy, g_y_xyz_yz_yz, g_y_xyz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yz_yz_xx[i] = 4.0 * g_y_xyz_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_yz_xy[i] = 4.0 * g_y_xyz_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_yz_xz[i] = 4.0 * g_y_xyz_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_yz_yy[i] = 4.0 * g_y_xyz_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_yz_yz[i] = 4.0 * g_y_xyz_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_yz_zz[i] = 4.0 * g_y_xyz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_y_x_0_0_0_yz_zz_xx, g_y_x_0_0_0_yz_zz_xy, g_y_x_0_0_0_yz_zz_xz, g_y_x_0_0_0_yz_zz_yy, g_y_x_0_0_0_yz_zz_yz, g_y_x_0_0_0_yz_zz_zz, g_y_xyz_zz_xx, g_y_xyz_zz_xy, g_y_xyz_zz_xz, g_y_xyz_zz_yy, g_y_xyz_zz_yz, g_y_xyz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yz_zz_xx[i] = 4.0 * g_y_xyz_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_zz_xy[i] = 4.0 * g_y_xyz_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_zz_xz[i] = 4.0 * g_y_xyz_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_zz_yy[i] = 4.0 * g_y_xyz_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_zz_yz[i] = 4.0 * g_y_xyz_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_zz_zz[i] = 4.0 * g_y_xyz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_y_x_0_0_0_zz_xx_xx, g_y_x_0_0_0_zz_xx_xy, g_y_x_0_0_0_zz_xx_xz, g_y_x_0_0_0_zz_xx_yy, g_y_x_0_0_0_zz_xx_yz, g_y_x_0_0_0_zz_xx_zz, g_y_xzz_xx_xx, g_y_xzz_xx_xy, g_y_xzz_xx_xz, g_y_xzz_xx_yy, g_y_xzz_xx_yz, g_y_xzz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_zz_xx_xx[i] = 4.0 * g_y_xzz_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_xx_xy[i] = 4.0 * g_y_xzz_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_xx_xz[i] = 4.0 * g_y_xzz_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_xx_yy[i] = 4.0 * g_y_xzz_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_xx_yz[i] = 4.0 * g_y_xzz_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_xx_zz[i] = 4.0 * g_y_xzz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_y_x_0_0_0_zz_xy_xx, g_y_x_0_0_0_zz_xy_xy, g_y_x_0_0_0_zz_xy_xz, g_y_x_0_0_0_zz_xy_yy, g_y_x_0_0_0_zz_xy_yz, g_y_x_0_0_0_zz_xy_zz, g_y_xzz_xy_xx, g_y_xzz_xy_xy, g_y_xzz_xy_xz, g_y_xzz_xy_yy, g_y_xzz_xy_yz, g_y_xzz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_zz_xy_xx[i] = 4.0 * g_y_xzz_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_xy_xy[i] = 4.0 * g_y_xzz_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_xy_xz[i] = 4.0 * g_y_xzz_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_xy_yy[i] = 4.0 * g_y_xzz_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_xy_yz[i] = 4.0 * g_y_xzz_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_xy_zz[i] = 4.0 * g_y_xzz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_y_x_0_0_0_zz_xz_xx, g_y_x_0_0_0_zz_xz_xy, g_y_x_0_0_0_zz_xz_xz, g_y_x_0_0_0_zz_xz_yy, g_y_x_0_0_0_zz_xz_yz, g_y_x_0_0_0_zz_xz_zz, g_y_xzz_xz_xx, g_y_xzz_xz_xy, g_y_xzz_xz_xz, g_y_xzz_xz_yy, g_y_xzz_xz_yz, g_y_xzz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_zz_xz_xx[i] = 4.0 * g_y_xzz_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_xz_xy[i] = 4.0 * g_y_xzz_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_xz_xz[i] = 4.0 * g_y_xzz_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_xz_yy[i] = 4.0 * g_y_xzz_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_xz_yz[i] = 4.0 * g_y_xzz_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_xz_zz[i] = 4.0 * g_y_xzz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_y_x_0_0_0_zz_yy_xx, g_y_x_0_0_0_zz_yy_xy, g_y_x_0_0_0_zz_yy_xz, g_y_x_0_0_0_zz_yy_yy, g_y_x_0_0_0_zz_yy_yz, g_y_x_0_0_0_zz_yy_zz, g_y_xzz_yy_xx, g_y_xzz_yy_xy, g_y_xzz_yy_xz, g_y_xzz_yy_yy, g_y_xzz_yy_yz, g_y_xzz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_zz_yy_xx[i] = 4.0 * g_y_xzz_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_yy_xy[i] = 4.0 * g_y_xzz_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_yy_xz[i] = 4.0 * g_y_xzz_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_yy_yy[i] = 4.0 * g_y_xzz_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_yy_yz[i] = 4.0 * g_y_xzz_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_yy_zz[i] = 4.0 * g_y_xzz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_y_x_0_0_0_zz_yz_xx, g_y_x_0_0_0_zz_yz_xy, g_y_x_0_0_0_zz_yz_xz, g_y_x_0_0_0_zz_yz_yy, g_y_x_0_0_0_zz_yz_yz, g_y_x_0_0_0_zz_yz_zz, g_y_xzz_yz_xx, g_y_xzz_yz_xy, g_y_xzz_yz_xz, g_y_xzz_yz_yy, g_y_xzz_yz_yz, g_y_xzz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_zz_yz_xx[i] = 4.0 * g_y_xzz_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_yz_xy[i] = 4.0 * g_y_xzz_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_yz_xz[i] = 4.0 * g_y_xzz_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_yz_yy[i] = 4.0 * g_y_xzz_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_yz_yz[i] = 4.0 * g_y_xzz_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_yz_zz[i] = 4.0 * g_y_xzz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_y_x_0_0_0_zz_zz_xx, g_y_x_0_0_0_zz_zz_xy, g_y_x_0_0_0_zz_zz_xz, g_y_x_0_0_0_zz_zz_yy, g_y_x_0_0_0_zz_zz_yz, g_y_x_0_0_0_zz_zz_zz, g_y_xzz_zz_xx, g_y_xzz_zz_xy, g_y_xzz_zz_xz, g_y_xzz_zz_yy, g_y_xzz_zz_yz, g_y_xzz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_zz_zz_xx[i] = 4.0 * g_y_xzz_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_zz_xy[i] = 4.0 * g_y_xzz_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_zz_xz[i] = 4.0 * g_y_xzz_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_zz_yy[i] = 4.0 * g_y_xzz_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_zz_yz[i] = 4.0 * g_y_xzz_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_zz_zz[i] = 4.0 * g_y_xzz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_y_xxy_xx_xx, g_y_xxy_xx_xy, g_y_xxy_xx_xz, g_y_xxy_xx_yy, g_y_xxy_xx_yz, g_y_xxy_xx_zz, g_y_y_0_0_0_xx_xx_xx, g_y_y_0_0_0_xx_xx_xy, g_y_y_0_0_0_xx_xx_xz, g_y_y_0_0_0_xx_xx_yy, g_y_y_0_0_0_xx_xx_yz, g_y_y_0_0_0_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xx_xx_xx[i] = 4.0 * g_y_xxy_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_xx_xy[i] = 4.0 * g_y_xxy_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_xx_xz[i] = 4.0 * g_y_xxy_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_xx_yy[i] = 4.0 * g_y_xxy_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_xx_yz[i] = 4.0 * g_y_xxy_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_xx_zz[i] = 4.0 * g_y_xxy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_y_xxy_xy_xx, g_y_xxy_xy_xy, g_y_xxy_xy_xz, g_y_xxy_xy_yy, g_y_xxy_xy_yz, g_y_xxy_xy_zz, g_y_y_0_0_0_xx_xy_xx, g_y_y_0_0_0_xx_xy_xy, g_y_y_0_0_0_xx_xy_xz, g_y_y_0_0_0_xx_xy_yy, g_y_y_0_0_0_xx_xy_yz, g_y_y_0_0_0_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xx_xy_xx[i] = 4.0 * g_y_xxy_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_xy_xy[i] = 4.0 * g_y_xxy_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_xy_xz[i] = 4.0 * g_y_xxy_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_xy_yy[i] = 4.0 * g_y_xxy_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_xy_yz[i] = 4.0 * g_y_xxy_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_xy_zz[i] = 4.0 * g_y_xxy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_y_xxy_xz_xx, g_y_xxy_xz_xy, g_y_xxy_xz_xz, g_y_xxy_xz_yy, g_y_xxy_xz_yz, g_y_xxy_xz_zz, g_y_y_0_0_0_xx_xz_xx, g_y_y_0_0_0_xx_xz_xy, g_y_y_0_0_0_xx_xz_xz, g_y_y_0_0_0_xx_xz_yy, g_y_y_0_0_0_xx_xz_yz, g_y_y_0_0_0_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xx_xz_xx[i] = 4.0 * g_y_xxy_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_xz_xy[i] = 4.0 * g_y_xxy_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_xz_xz[i] = 4.0 * g_y_xxy_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_xz_yy[i] = 4.0 * g_y_xxy_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_xz_yz[i] = 4.0 * g_y_xxy_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_xz_zz[i] = 4.0 * g_y_xxy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_y_xxy_yy_xx, g_y_xxy_yy_xy, g_y_xxy_yy_xz, g_y_xxy_yy_yy, g_y_xxy_yy_yz, g_y_xxy_yy_zz, g_y_y_0_0_0_xx_yy_xx, g_y_y_0_0_0_xx_yy_xy, g_y_y_0_0_0_xx_yy_xz, g_y_y_0_0_0_xx_yy_yy, g_y_y_0_0_0_xx_yy_yz, g_y_y_0_0_0_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xx_yy_xx[i] = 4.0 * g_y_xxy_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_yy_xy[i] = 4.0 * g_y_xxy_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_yy_xz[i] = 4.0 * g_y_xxy_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_yy_yy[i] = 4.0 * g_y_xxy_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_yy_yz[i] = 4.0 * g_y_xxy_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_yy_zz[i] = 4.0 * g_y_xxy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_y_xxy_yz_xx, g_y_xxy_yz_xy, g_y_xxy_yz_xz, g_y_xxy_yz_yy, g_y_xxy_yz_yz, g_y_xxy_yz_zz, g_y_y_0_0_0_xx_yz_xx, g_y_y_0_0_0_xx_yz_xy, g_y_y_0_0_0_xx_yz_xz, g_y_y_0_0_0_xx_yz_yy, g_y_y_0_0_0_xx_yz_yz, g_y_y_0_0_0_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xx_yz_xx[i] = 4.0 * g_y_xxy_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_yz_xy[i] = 4.0 * g_y_xxy_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_yz_xz[i] = 4.0 * g_y_xxy_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_yz_yy[i] = 4.0 * g_y_xxy_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_yz_yz[i] = 4.0 * g_y_xxy_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_yz_zz[i] = 4.0 * g_y_xxy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_y_xxy_zz_xx, g_y_xxy_zz_xy, g_y_xxy_zz_xz, g_y_xxy_zz_yy, g_y_xxy_zz_yz, g_y_xxy_zz_zz, g_y_y_0_0_0_xx_zz_xx, g_y_y_0_0_0_xx_zz_xy, g_y_y_0_0_0_xx_zz_xz, g_y_y_0_0_0_xx_zz_yy, g_y_y_0_0_0_xx_zz_yz, g_y_y_0_0_0_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xx_zz_xx[i] = 4.0 * g_y_xxy_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_zz_xy[i] = 4.0 * g_y_xxy_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_zz_xz[i] = 4.0 * g_y_xxy_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_zz_yy[i] = 4.0 * g_y_xxy_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_zz_yz[i] = 4.0 * g_y_xxy_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_zz_zz[i] = 4.0 * g_y_xxy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_y_x_xx_xx, g_y_x_xx_xy, g_y_x_xx_xz, g_y_x_xx_yy, g_y_x_xx_yz, g_y_x_xx_zz, g_y_xyy_xx_xx, g_y_xyy_xx_xy, g_y_xyy_xx_xz, g_y_xyy_xx_yy, g_y_xyy_xx_yz, g_y_xyy_xx_zz, g_y_y_0_0_0_xy_xx_xx, g_y_y_0_0_0_xy_xx_xy, g_y_y_0_0_0_xy_xx_xz, g_y_y_0_0_0_xy_xx_yy, g_y_y_0_0_0_xy_xx_yz, g_y_y_0_0_0_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xy_xx_xx[i] = -2.0 * g_y_x_xx_xx[i] * a_exp + 4.0 * g_y_xyy_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_xx_xy[i] = -2.0 * g_y_x_xx_xy[i] * a_exp + 4.0 * g_y_xyy_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_xx_xz[i] = -2.0 * g_y_x_xx_xz[i] * a_exp + 4.0 * g_y_xyy_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_xx_yy[i] = -2.0 * g_y_x_xx_yy[i] * a_exp + 4.0 * g_y_xyy_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_xx_yz[i] = -2.0 * g_y_x_xx_yz[i] * a_exp + 4.0 * g_y_xyy_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_xx_zz[i] = -2.0 * g_y_x_xx_zz[i] * a_exp + 4.0 * g_y_xyy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_y_x_xy_xx, g_y_x_xy_xy, g_y_x_xy_xz, g_y_x_xy_yy, g_y_x_xy_yz, g_y_x_xy_zz, g_y_xyy_xy_xx, g_y_xyy_xy_xy, g_y_xyy_xy_xz, g_y_xyy_xy_yy, g_y_xyy_xy_yz, g_y_xyy_xy_zz, g_y_y_0_0_0_xy_xy_xx, g_y_y_0_0_0_xy_xy_xy, g_y_y_0_0_0_xy_xy_xz, g_y_y_0_0_0_xy_xy_yy, g_y_y_0_0_0_xy_xy_yz, g_y_y_0_0_0_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xy_xy_xx[i] = -2.0 * g_y_x_xy_xx[i] * a_exp + 4.0 * g_y_xyy_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_xy_xy[i] = -2.0 * g_y_x_xy_xy[i] * a_exp + 4.0 * g_y_xyy_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_xy_xz[i] = -2.0 * g_y_x_xy_xz[i] * a_exp + 4.0 * g_y_xyy_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_xy_yy[i] = -2.0 * g_y_x_xy_yy[i] * a_exp + 4.0 * g_y_xyy_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_xy_yz[i] = -2.0 * g_y_x_xy_yz[i] * a_exp + 4.0 * g_y_xyy_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_xy_zz[i] = -2.0 * g_y_x_xy_zz[i] * a_exp + 4.0 * g_y_xyy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_y_x_xz_xx, g_y_x_xz_xy, g_y_x_xz_xz, g_y_x_xz_yy, g_y_x_xz_yz, g_y_x_xz_zz, g_y_xyy_xz_xx, g_y_xyy_xz_xy, g_y_xyy_xz_xz, g_y_xyy_xz_yy, g_y_xyy_xz_yz, g_y_xyy_xz_zz, g_y_y_0_0_0_xy_xz_xx, g_y_y_0_0_0_xy_xz_xy, g_y_y_0_0_0_xy_xz_xz, g_y_y_0_0_0_xy_xz_yy, g_y_y_0_0_0_xy_xz_yz, g_y_y_0_0_0_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xy_xz_xx[i] = -2.0 * g_y_x_xz_xx[i] * a_exp + 4.0 * g_y_xyy_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_xz_xy[i] = -2.0 * g_y_x_xz_xy[i] * a_exp + 4.0 * g_y_xyy_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_xz_xz[i] = -2.0 * g_y_x_xz_xz[i] * a_exp + 4.0 * g_y_xyy_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_xz_yy[i] = -2.0 * g_y_x_xz_yy[i] * a_exp + 4.0 * g_y_xyy_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_xz_yz[i] = -2.0 * g_y_x_xz_yz[i] * a_exp + 4.0 * g_y_xyy_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_xz_zz[i] = -2.0 * g_y_x_xz_zz[i] * a_exp + 4.0 * g_y_xyy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_y_x_yy_xx, g_y_x_yy_xy, g_y_x_yy_xz, g_y_x_yy_yy, g_y_x_yy_yz, g_y_x_yy_zz, g_y_xyy_yy_xx, g_y_xyy_yy_xy, g_y_xyy_yy_xz, g_y_xyy_yy_yy, g_y_xyy_yy_yz, g_y_xyy_yy_zz, g_y_y_0_0_0_xy_yy_xx, g_y_y_0_0_0_xy_yy_xy, g_y_y_0_0_0_xy_yy_xz, g_y_y_0_0_0_xy_yy_yy, g_y_y_0_0_0_xy_yy_yz, g_y_y_0_0_0_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xy_yy_xx[i] = -2.0 * g_y_x_yy_xx[i] * a_exp + 4.0 * g_y_xyy_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_yy_xy[i] = -2.0 * g_y_x_yy_xy[i] * a_exp + 4.0 * g_y_xyy_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_yy_xz[i] = -2.0 * g_y_x_yy_xz[i] * a_exp + 4.0 * g_y_xyy_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_yy_yy[i] = -2.0 * g_y_x_yy_yy[i] * a_exp + 4.0 * g_y_xyy_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_yy_yz[i] = -2.0 * g_y_x_yy_yz[i] * a_exp + 4.0 * g_y_xyy_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_yy_zz[i] = -2.0 * g_y_x_yy_zz[i] * a_exp + 4.0 * g_y_xyy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_y_x_yz_xx, g_y_x_yz_xy, g_y_x_yz_xz, g_y_x_yz_yy, g_y_x_yz_yz, g_y_x_yz_zz, g_y_xyy_yz_xx, g_y_xyy_yz_xy, g_y_xyy_yz_xz, g_y_xyy_yz_yy, g_y_xyy_yz_yz, g_y_xyy_yz_zz, g_y_y_0_0_0_xy_yz_xx, g_y_y_0_0_0_xy_yz_xy, g_y_y_0_0_0_xy_yz_xz, g_y_y_0_0_0_xy_yz_yy, g_y_y_0_0_0_xy_yz_yz, g_y_y_0_0_0_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xy_yz_xx[i] = -2.0 * g_y_x_yz_xx[i] * a_exp + 4.0 * g_y_xyy_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_yz_xy[i] = -2.0 * g_y_x_yz_xy[i] * a_exp + 4.0 * g_y_xyy_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_yz_xz[i] = -2.0 * g_y_x_yz_xz[i] * a_exp + 4.0 * g_y_xyy_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_yz_yy[i] = -2.0 * g_y_x_yz_yy[i] * a_exp + 4.0 * g_y_xyy_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_yz_yz[i] = -2.0 * g_y_x_yz_yz[i] * a_exp + 4.0 * g_y_xyy_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_yz_zz[i] = -2.0 * g_y_x_yz_zz[i] * a_exp + 4.0 * g_y_xyy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_y_x_zz_xx, g_y_x_zz_xy, g_y_x_zz_xz, g_y_x_zz_yy, g_y_x_zz_yz, g_y_x_zz_zz, g_y_xyy_zz_xx, g_y_xyy_zz_xy, g_y_xyy_zz_xz, g_y_xyy_zz_yy, g_y_xyy_zz_yz, g_y_xyy_zz_zz, g_y_y_0_0_0_xy_zz_xx, g_y_y_0_0_0_xy_zz_xy, g_y_y_0_0_0_xy_zz_xz, g_y_y_0_0_0_xy_zz_yy, g_y_y_0_0_0_xy_zz_yz, g_y_y_0_0_0_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xy_zz_xx[i] = -2.0 * g_y_x_zz_xx[i] * a_exp + 4.0 * g_y_xyy_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_zz_xy[i] = -2.0 * g_y_x_zz_xy[i] * a_exp + 4.0 * g_y_xyy_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_zz_xz[i] = -2.0 * g_y_x_zz_xz[i] * a_exp + 4.0 * g_y_xyy_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_zz_yy[i] = -2.0 * g_y_x_zz_yy[i] * a_exp + 4.0 * g_y_xyy_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_zz_yz[i] = -2.0 * g_y_x_zz_yz[i] * a_exp + 4.0 * g_y_xyy_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_zz_zz[i] = -2.0 * g_y_x_zz_zz[i] * a_exp + 4.0 * g_y_xyy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_y_xyz_xx_xx, g_y_xyz_xx_xy, g_y_xyz_xx_xz, g_y_xyz_xx_yy, g_y_xyz_xx_yz, g_y_xyz_xx_zz, g_y_y_0_0_0_xz_xx_xx, g_y_y_0_0_0_xz_xx_xy, g_y_y_0_0_0_xz_xx_xz, g_y_y_0_0_0_xz_xx_yy, g_y_y_0_0_0_xz_xx_yz, g_y_y_0_0_0_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xz_xx_xx[i] = 4.0 * g_y_xyz_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_xx_xy[i] = 4.0 * g_y_xyz_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_xx_xz[i] = 4.0 * g_y_xyz_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_xx_yy[i] = 4.0 * g_y_xyz_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_xx_yz[i] = 4.0 * g_y_xyz_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_xx_zz[i] = 4.0 * g_y_xyz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_y_xyz_xy_xx, g_y_xyz_xy_xy, g_y_xyz_xy_xz, g_y_xyz_xy_yy, g_y_xyz_xy_yz, g_y_xyz_xy_zz, g_y_y_0_0_0_xz_xy_xx, g_y_y_0_0_0_xz_xy_xy, g_y_y_0_0_0_xz_xy_xz, g_y_y_0_0_0_xz_xy_yy, g_y_y_0_0_0_xz_xy_yz, g_y_y_0_0_0_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xz_xy_xx[i] = 4.0 * g_y_xyz_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_xy_xy[i] = 4.0 * g_y_xyz_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_xy_xz[i] = 4.0 * g_y_xyz_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_xy_yy[i] = 4.0 * g_y_xyz_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_xy_yz[i] = 4.0 * g_y_xyz_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_xy_zz[i] = 4.0 * g_y_xyz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_y_xyz_xz_xx, g_y_xyz_xz_xy, g_y_xyz_xz_xz, g_y_xyz_xz_yy, g_y_xyz_xz_yz, g_y_xyz_xz_zz, g_y_y_0_0_0_xz_xz_xx, g_y_y_0_0_0_xz_xz_xy, g_y_y_0_0_0_xz_xz_xz, g_y_y_0_0_0_xz_xz_yy, g_y_y_0_0_0_xz_xz_yz, g_y_y_0_0_0_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xz_xz_xx[i] = 4.0 * g_y_xyz_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_xz_xy[i] = 4.0 * g_y_xyz_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_xz_xz[i] = 4.0 * g_y_xyz_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_xz_yy[i] = 4.0 * g_y_xyz_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_xz_yz[i] = 4.0 * g_y_xyz_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_xz_zz[i] = 4.0 * g_y_xyz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_y_xyz_yy_xx, g_y_xyz_yy_xy, g_y_xyz_yy_xz, g_y_xyz_yy_yy, g_y_xyz_yy_yz, g_y_xyz_yy_zz, g_y_y_0_0_0_xz_yy_xx, g_y_y_0_0_0_xz_yy_xy, g_y_y_0_0_0_xz_yy_xz, g_y_y_0_0_0_xz_yy_yy, g_y_y_0_0_0_xz_yy_yz, g_y_y_0_0_0_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xz_yy_xx[i] = 4.0 * g_y_xyz_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_yy_xy[i] = 4.0 * g_y_xyz_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_yy_xz[i] = 4.0 * g_y_xyz_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_yy_yy[i] = 4.0 * g_y_xyz_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_yy_yz[i] = 4.0 * g_y_xyz_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_yy_zz[i] = 4.0 * g_y_xyz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_y_xyz_yz_xx, g_y_xyz_yz_xy, g_y_xyz_yz_xz, g_y_xyz_yz_yy, g_y_xyz_yz_yz, g_y_xyz_yz_zz, g_y_y_0_0_0_xz_yz_xx, g_y_y_0_0_0_xz_yz_xy, g_y_y_0_0_0_xz_yz_xz, g_y_y_0_0_0_xz_yz_yy, g_y_y_0_0_0_xz_yz_yz, g_y_y_0_0_0_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xz_yz_xx[i] = 4.0 * g_y_xyz_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_yz_xy[i] = 4.0 * g_y_xyz_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_yz_xz[i] = 4.0 * g_y_xyz_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_yz_yy[i] = 4.0 * g_y_xyz_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_yz_yz[i] = 4.0 * g_y_xyz_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_yz_zz[i] = 4.0 * g_y_xyz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_y_xyz_zz_xx, g_y_xyz_zz_xy, g_y_xyz_zz_xz, g_y_xyz_zz_yy, g_y_xyz_zz_yz, g_y_xyz_zz_zz, g_y_y_0_0_0_xz_zz_xx, g_y_y_0_0_0_xz_zz_xy, g_y_y_0_0_0_xz_zz_xz, g_y_y_0_0_0_xz_zz_yy, g_y_y_0_0_0_xz_zz_yz, g_y_y_0_0_0_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xz_zz_xx[i] = 4.0 * g_y_xyz_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_zz_xy[i] = 4.0 * g_y_xyz_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_zz_xz[i] = 4.0 * g_y_xyz_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_zz_yy[i] = 4.0 * g_y_xyz_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_zz_yz[i] = 4.0 * g_y_xyz_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_zz_zz[i] = 4.0 * g_y_xyz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (972-978)

    #pragma omp simd aligned(g_y_y_0_0_0_yy_xx_xx, g_y_y_0_0_0_yy_xx_xy, g_y_y_0_0_0_yy_xx_xz, g_y_y_0_0_0_yy_xx_yy, g_y_y_0_0_0_yy_xx_yz, g_y_y_0_0_0_yy_xx_zz, g_y_y_xx_xx, g_y_y_xx_xy, g_y_y_xx_xz, g_y_y_xx_yy, g_y_y_xx_yz, g_y_y_xx_zz, g_y_yyy_xx_xx, g_y_yyy_xx_xy, g_y_yyy_xx_xz, g_y_yyy_xx_yy, g_y_yyy_xx_yz, g_y_yyy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yy_xx_xx[i] = -4.0 * g_y_y_xx_xx[i] * a_exp + 4.0 * g_y_yyy_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_xx_xy[i] = -4.0 * g_y_y_xx_xy[i] * a_exp + 4.0 * g_y_yyy_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_xx_xz[i] = -4.0 * g_y_y_xx_xz[i] * a_exp + 4.0 * g_y_yyy_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_xx_yy[i] = -4.0 * g_y_y_xx_yy[i] * a_exp + 4.0 * g_y_yyy_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_xx_yz[i] = -4.0 * g_y_y_xx_yz[i] * a_exp + 4.0 * g_y_yyy_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_xx_zz[i] = -4.0 * g_y_y_xx_zz[i] * a_exp + 4.0 * g_y_yyy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (978-984)

    #pragma omp simd aligned(g_y_y_0_0_0_yy_xy_xx, g_y_y_0_0_0_yy_xy_xy, g_y_y_0_0_0_yy_xy_xz, g_y_y_0_0_0_yy_xy_yy, g_y_y_0_0_0_yy_xy_yz, g_y_y_0_0_0_yy_xy_zz, g_y_y_xy_xx, g_y_y_xy_xy, g_y_y_xy_xz, g_y_y_xy_yy, g_y_y_xy_yz, g_y_y_xy_zz, g_y_yyy_xy_xx, g_y_yyy_xy_xy, g_y_yyy_xy_xz, g_y_yyy_xy_yy, g_y_yyy_xy_yz, g_y_yyy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yy_xy_xx[i] = -4.0 * g_y_y_xy_xx[i] * a_exp + 4.0 * g_y_yyy_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_xy_xy[i] = -4.0 * g_y_y_xy_xy[i] * a_exp + 4.0 * g_y_yyy_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_xy_xz[i] = -4.0 * g_y_y_xy_xz[i] * a_exp + 4.0 * g_y_yyy_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_xy_yy[i] = -4.0 * g_y_y_xy_yy[i] * a_exp + 4.0 * g_y_yyy_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_xy_yz[i] = -4.0 * g_y_y_xy_yz[i] * a_exp + 4.0 * g_y_yyy_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_xy_zz[i] = -4.0 * g_y_y_xy_zz[i] * a_exp + 4.0 * g_y_yyy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (984-990)

    #pragma omp simd aligned(g_y_y_0_0_0_yy_xz_xx, g_y_y_0_0_0_yy_xz_xy, g_y_y_0_0_0_yy_xz_xz, g_y_y_0_0_0_yy_xz_yy, g_y_y_0_0_0_yy_xz_yz, g_y_y_0_0_0_yy_xz_zz, g_y_y_xz_xx, g_y_y_xz_xy, g_y_y_xz_xz, g_y_y_xz_yy, g_y_y_xz_yz, g_y_y_xz_zz, g_y_yyy_xz_xx, g_y_yyy_xz_xy, g_y_yyy_xz_xz, g_y_yyy_xz_yy, g_y_yyy_xz_yz, g_y_yyy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yy_xz_xx[i] = -4.0 * g_y_y_xz_xx[i] * a_exp + 4.0 * g_y_yyy_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_xz_xy[i] = -4.0 * g_y_y_xz_xy[i] * a_exp + 4.0 * g_y_yyy_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_xz_xz[i] = -4.0 * g_y_y_xz_xz[i] * a_exp + 4.0 * g_y_yyy_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_xz_yy[i] = -4.0 * g_y_y_xz_yy[i] * a_exp + 4.0 * g_y_yyy_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_xz_yz[i] = -4.0 * g_y_y_xz_yz[i] * a_exp + 4.0 * g_y_yyy_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_xz_zz[i] = -4.0 * g_y_y_xz_zz[i] * a_exp + 4.0 * g_y_yyy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (990-996)

    #pragma omp simd aligned(g_y_y_0_0_0_yy_yy_xx, g_y_y_0_0_0_yy_yy_xy, g_y_y_0_0_0_yy_yy_xz, g_y_y_0_0_0_yy_yy_yy, g_y_y_0_0_0_yy_yy_yz, g_y_y_0_0_0_yy_yy_zz, g_y_y_yy_xx, g_y_y_yy_xy, g_y_y_yy_xz, g_y_y_yy_yy, g_y_y_yy_yz, g_y_y_yy_zz, g_y_yyy_yy_xx, g_y_yyy_yy_xy, g_y_yyy_yy_xz, g_y_yyy_yy_yy, g_y_yyy_yy_yz, g_y_yyy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yy_yy_xx[i] = -4.0 * g_y_y_yy_xx[i] * a_exp + 4.0 * g_y_yyy_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_yy_xy[i] = -4.0 * g_y_y_yy_xy[i] * a_exp + 4.0 * g_y_yyy_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_yy_xz[i] = -4.0 * g_y_y_yy_xz[i] * a_exp + 4.0 * g_y_yyy_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_yy_yy[i] = -4.0 * g_y_y_yy_yy[i] * a_exp + 4.0 * g_y_yyy_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_yy_yz[i] = -4.0 * g_y_y_yy_yz[i] * a_exp + 4.0 * g_y_yyy_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_yy_zz[i] = -4.0 * g_y_y_yy_zz[i] * a_exp + 4.0 * g_y_yyy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (996-1002)

    #pragma omp simd aligned(g_y_y_0_0_0_yy_yz_xx, g_y_y_0_0_0_yy_yz_xy, g_y_y_0_0_0_yy_yz_xz, g_y_y_0_0_0_yy_yz_yy, g_y_y_0_0_0_yy_yz_yz, g_y_y_0_0_0_yy_yz_zz, g_y_y_yz_xx, g_y_y_yz_xy, g_y_y_yz_xz, g_y_y_yz_yy, g_y_y_yz_yz, g_y_y_yz_zz, g_y_yyy_yz_xx, g_y_yyy_yz_xy, g_y_yyy_yz_xz, g_y_yyy_yz_yy, g_y_yyy_yz_yz, g_y_yyy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yy_yz_xx[i] = -4.0 * g_y_y_yz_xx[i] * a_exp + 4.0 * g_y_yyy_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_yz_xy[i] = -4.0 * g_y_y_yz_xy[i] * a_exp + 4.0 * g_y_yyy_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_yz_xz[i] = -4.0 * g_y_y_yz_xz[i] * a_exp + 4.0 * g_y_yyy_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_yz_yy[i] = -4.0 * g_y_y_yz_yy[i] * a_exp + 4.0 * g_y_yyy_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_yz_yz[i] = -4.0 * g_y_y_yz_yz[i] * a_exp + 4.0 * g_y_yyy_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_yz_zz[i] = -4.0 * g_y_y_yz_zz[i] * a_exp + 4.0 * g_y_yyy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1002-1008)

    #pragma omp simd aligned(g_y_y_0_0_0_yy_zz_xx, g_y_y_0_0_0_yy_zz_xy, g_y_y_0_0_0_yy_zz_xz, g_y_y_0_0_0_yy_zz_yy, g_y_y_0_0_0_yy_zz_yz, g_y_y_0_0_0_yy_zz_zz, g_y_y_zz_xx, g_y_y_zz_xy, g_y_y_zz_xz, g_y_y_zz_yy, g_y_y_zz_yz, g_y_y_zz_zz, g_y_yyy_zz_xx, g_y_yyy_zz_xy, g_y_yyy_zz_xz, g_y_yyy_zz_yy, g_y_yyy_zz_yz, g_y_yyy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yy_zz_xx[i] = -4.0 * g_y_y_zz_xx[i] * a_exp + 4.0 * g_y_yyy_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_zz_xy[i] = -4.0 * g_y_y_zz_xy[i] * a_exp + 4.0 * g_y_yyy_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_zz_xz[i] = -4.0 * g_y_y_zz_xz[i] * a_exp + 4.0 * g_y_yyy_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_zz_yy[i] = -4.0 * g_y_y_zz_yy[i] * a_exp + 4.0 * g_y_yyy_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_zz_yz[i] = -4.0 * g_y_y_zz_yz[i] * a_exp + 4.0 * g_y_yyy_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_zz_zz[i] = -4.0 * g_y_y_zz_zz[i] * a_exp + 4.0 * g_y_yyy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1008-1014)

    #pragma omp simd aligned(g_y_y_0_0_0_yz_xx_xx, g_y_y_0_0_0_yz_xx_xy, g_y_y_0_0_0_yz_xx_xz, g_y_y_0_0_0_yz_xx_yy, g_y_y_0_0_0_yz_xx_yz, g_y_y_0_0_0_yz_xx_zz, g_y_yyz_xx_xx, g_y_yyz_xx_xy, g_y_yyz_xx_xz, g_y_yyz_xx_yy, g_y_yyz_xx_yz, g_y_yyz_xx_zz, g_y_z_xx_xx, g_y_z_xx_xy, g_y_z_xx_xz, g_y_z_xx_yy, g_y_z_xx_yz, g_y_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yz_xx_xx[i] = -2.0 * g_y_z_xx_xx[i] * a_exp + 4.0 * g_y_yyz_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_xx_xy[i] = -2.0 * g_y_z_xx_xy[i] * a_exp + 4.0 * g_y_yyz_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_xx_xz[i] = -2.0 * g_y_z_xx_xz[i] * a_exp + 4.0 * g_y_yyz_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_xx_yy[i] = -2.0 * g_y_z_xx_yy[i] * a_exp + 4.0 * g_y_yyz_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_xx_yz[i] = -2.0 * g_y_z_xx_yz[i] * a_exp + 4.0 * g_y_yyz_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_xx_zz[i] = -2.0 * g_y_z_xx_zz[i] * a_exp + 4.0 * g_y_yyz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1014-1020)

    #pragma omp simd aligned(g_y_y_0_0_0_yz_xy_xx, g_y_y_0_0_0_yz_xy_xy, g_y_y_0_0_0_yz_xy_xz, g_y_y_0_0_0_yz_xy_yy, g_y_y_0_0_0_yz_xy_yz, g_y_y_0_0_0_yz_xy_zz, g_y_yyz_xy_xx, g_y_yyz_xy_xy, g_y_yyz_xy_xz, g_y_yyz_xy_yy, g_y_yyz_xy_yz, g_y_yyz_xy_zz, g_y_z_xy_xx, g_y_z_xy_xy, g_y_z_xy_xz, g_y_z_xy_yy, g_y_z_xy_yz, g_y_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yz_xy_xx[i] = -2.0 * g_y_z_xy_xx[i] * a_exp + 4.0 * g_y_yyz_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_xy_xy[i] = -2.0 * g_y_z_xy_xy[i] * a_exp + 4.0 * g_y_yyz_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_xy_xz[i] = -2.0 * g_y_z_xy_xz[i] * a_exp + 4.0 * g_y_yyz_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_xy_yy[i] = -2.0 * g_y_z_xy_yy[i] * a_exp + 4.0 * g_y_yyz_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_xy_yz[i] = -2.0 * g_y_z_xy_yz[i] * a_exp + 4.0 * g_y_yyz_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_xy_zz[i] = -2.0 * g_y_z_xy_zz[i] * a_exp + 4.0 * g_y_yyz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1020-1026)

    #pragma omp simd aligned(g_y_y_0_0_0_yz_xz_xx, g_y_y_0_0_0_yz_xz_xy, g_y_y_0_0_0_yz_xz_xz, g_y_y_0_0_0_yz_xz_yy, g_y_y_0_0_0_yz_xz_yz, g_y_y_0_0_0_yz_xz_zz, g_y_yyz_xz_xx, g_y_yyz_xz_xy, g_y_yyz_xz_xz, g_y_yyz_xz_yy, g_y_yyz_xz_yz, g_y_yyz_xz_zz, g_y_z_xz_xx, g_y_z_xz_xy, g_y_z_xz_xz, g_y_z_xz_yy, g_y_z_xz_yz, g_y_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yz_xz_xx[i] = -2.0 * g_y_z_xz_xx[i] * a_exp + 4.0 * g_y_yyz_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_xz_xy[i] = -2.0 * g_y_z_xz_xy[i] * a_exp + 4.0 * g_y_yyz_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_xz_xz[i] = -2.0 * g_y_z_xz_xz[i] * a_exp + 4.0 * g_y_yyz_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_xz_yy[i] = -2.0 * g_y_z_xz_yy[i] * a_exp + 4.0 * g_y_yyz_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_xz_yz[i] = -2.0 * g_y_z_xz_yz[i] * a_exp + 4.0 * g_y_yyz_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_xz_zz[i] = -2.0 * g_y_z_xz_zz[i] * a_exp + 4.0 * g_y_yyz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1026-1032)

    #pragma omp simd aligned(g_y_y_0_0_0_yz_yy_xx, g_y_y_0_0_0_yz_yy_xy, g_y_y_0_0_0_yz_yy_xz, g_y_y_0_0_0_yz_yy_yy, g_y_y_0_0_0_yz_yy_yz, g_y_y_0_0_0_yz_yy_zz, g_y_yyz_yy_xx, g_y_yyz_yy_xy, g_y_yyz_yy_xz, g_y_yyz_yy_yy, g_y_yyz_yy_yz, g_y_yyz_yy_zz, g_y_z_yy_xx, g_y_z_yy_xy, g_y_z_yy_xz, g_y_z_yy_yy, g_y_z_yy_yz, g_y_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yz_yy_xx[i] = -2.0 * g_y_z_yy_xx[i] * a_exp + 4.0 * g_y_yyz_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_yy_xy[i] = -2.0 * g_y_z_yy_xy[i] * a_exp + 4.0 * g_y_yyz_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_yy_xz[i] = -2.0 * g_y_z_yy_xz[i] * a_exp + 4.0 * g_y_yyz_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_yy_yy[i] = -2.0 * g_y_z_yy_yy[i] * a_exp + 4.0 * g_y_yyz_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_yy_yz[i] = -2.0 * g_y_z_yy_yz[i] * a_exp + 4.0 * g_y_yyz_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_yy_zz[i] = -2.0 * g_y_z_yy_zz[i] * a_exp + 4.0 * g_y_yyz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1032-1038)

    #pragma omp simd aligned(g_y_y_0_0_0_yz_yz_xx, g_y_y_0_0_0_yz_yz_xy, g_y_y_0_0_0_yz_yz_xz, g_y_y_0_0_0_yz_yz_yy, g_y_y_0_0_0_yz_yz_yz, g_y_y_0_0_0_yz_yz_zz, g_y_yyz_yz_xx, g_y_yyz_yz_xy, g_y_yyz_yz_xz, g_y_yyz_yz_yy, g_y_yyz_yz_yz, g_y_yyz_yz_zz, g_y_z_yz_xx, g_y_z_yz_xy, g_y_z_yz_xz, g_y_z_yz_yy, g_y_z_yz_yz, g_y_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yz_yz_xx[i] = -2.0 * g_y_z_yz_xx[i] * a_exp + 4.0 * g_y_yyz_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_yz_xy[i] = -2.0 * g_y_z_yz_xy[i] * a_exp + 4.0 * g_y_yyz_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_yz_xz[i] = -2.0 * g_y_z_yz_xz[i] * a_exp + 4.0 * g_y_yyz_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_yz_yy[i] = -2.0 * g_y_z_yz_yy[i] * a_exp + 4.0 * g_y_yyz_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_yz_yz[i] = -2.0 * g_y_z_yz_yz[i] * a_exp + 4.0 * g_y_yyz_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_yz_zz[i] = -2.0 * g_y_z_yz_zz[i] * a_exp + 4.0 * g_y_yyz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1038-1044)

    #pragma omp simd aligned(g_y_y_0_0_0_yz_zz_xx, g_y_y_0_0_0_yz_zz_xy, g_y_y_0_0_0_yz_zz_xz, g_y_y_0_0_0_yz_zz_yy, g_y_y_0_0_0_yz_zz_yz, g_y_y_0_0_0_yz_zz_zz, g_y_yyz_zz_xx, g_y_yyz_zz_xy, g_y_yyz_zz_xz, g_y_yyz_zz_yy, g_y_yyz_zz_yz, g_y_yyz_zz_zz, g_y_z_zz_xx, g_y_z_zz_xy, g_y_z_zz_xz, g_y_z_zz_yy, g_y_z_zz_yz, g_y_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yz_zz_xx[i] = -2.0 * g_y_z_zz_xx[i] * a_exp + 4.0 * g_y_yyz_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_zz_xy[i] = -2.0 * g_y_z_zz_xy[i] * a_exp + 4.0 * g_y_yyz_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_zz_xz[i] = -2.0 * g_y_z_zz_xz[i] * a_exp + 4.0 * g_y_yyz_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_zz_yy[i] = -2.0 * g_y_z_zz_yy[i] * a_exp + 4.0 * g_y_yyz_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_zz_yz[i] = -2.0 * g_y_z_zz_yz[i] * a_exp + 4.0 * g_y_yyz_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_zz_zz[i] = -2.0 * g_y_z_zz_zz[i] * a_exp + 4.0 * g_y_yyz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1044-1050)

    #pragma omp simd aligned(g_y_y_0_0_0_zz_xx_xx, g_y_y_0_0_0_zz_xx_xy, g_y_y_0_0_0_zz_xx_xz, g_y_y_0_0_0_zz_xx_yy, g_y_y_0_0_0_zz_xx_yz, g_y_y_0_0_0_zz_xx_zz, g_y_yzz_xx_xx, g_y_yzz_xx_xy, g_y_yzz_xx_xz, g_y_yzz_xx_yy, g_y_yzz_xx_yz, g_y_yzz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_zz_xx_xx[i] = 4.0 * g_y_yzz_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_xx_xy[i] = 4.0 * g_y_yzz_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_xx_xz[i] = 4.0 * g_y_yzz_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_xx_yy[i] = 4.0 * g_y_yzz_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_xx_yz[i] = 4.0 * g_y_yzz_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_xx_zz[i] = 4.0 * g_y_yzz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1050-1056)

    #pragma omp simd aligned(g_y_y_0_0_0_zz_xy_xx, g_y_y_0_0_0_zz_xy_xy, g_y_y_0_0_0_zz_xy_xz, g_y_y_0_0_0_zz_xy_yy, g_y_y_0_0_0_zz_xy_yz, g_y_y_0_0_0_zz_xy_zz, g_y_yzz_xy_xx, g_y_yzz_xy_xy, g_y_yzz_xy_xz, g_y_yzz_xy_yy, g_y_yzz_xy_yz, g_y_yzz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_zz_xy_xx[i] = 4.0 * g_y_yzz_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_xy_xy[i] = 4.0 * g_y_yzz_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_xy_xz[i] = 4.0 * g_y_yzz_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_xy_yy[i] = 4.0 * g_y_yzz_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_xy_yz[i] = 4.0 * g_y_yzz_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_xy_zz[i] = 4.0 * g_y_yzz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1056-1062)

    #pragma omp simd aligned(g_y_y_0_0_0_zz_xz_xx, g_y_y_0_0_0_zz_xz_xy, g_y_y_0_0_0_zz_xz_xz, g_y_y_0_0_0_zz_xz_yy, g_y_y_0_0_0_zz_xz_yz, g_y_y_0_0_0_zz_xz_zz, g_y_yzz_xz_xx, g_y_yzz_xz_xy, g_y_yzz_xz_xz, g_y_yzz_xz_yy, g_y_yzz_xz_yz, g_y_yzz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_zz_xz_xx[i] = 4.0 * g_y_yzz_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_xz_xy[i] = 4.0 * g_y_yzz_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_xz_xz[i] = 4.0 * g_y_yzz_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_xz_yy[i] = 4.0 * g_y_yzz_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_xz_yz[i] = 4.0 * g_y_yzz_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_xz_zz[i] = 4.0 * g_y_yzz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1062-1068)

    #pragma omp simd aligned(g_y_y_0_0_0_zz_yy_xx, g_y_y_0_0_0_zz_yy_xy, g_y_y_0_0_0_zz_yy_xz, g_y_y_0_0_0_zz_yy_yy, g_y_y_0_0_0_zz_yy_yz, g_y_y_0_0_0_zz_yy_zz, g_y_yzz_yy_xx, g_y_yzz_yy_xy, g_y_yzz_yy_xz, g_y_yzz_yy_yy, g_y_yzz_yy_yz, g_y_yzz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_zz_yy_xx[i] = 4.0 * g_y_yzz_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_yy_xy[i] = 4.0 * g_y_yzz_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_yy_xz[i] = 4.0 * g_y_yzz_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_yy_yy[i] = 4.0 * g_y_yzz_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_yy_yz[i] = 4.0 * g_y_yzz_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_yy_zz[i] = 4.0 * g_y_yzz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1068-1074)

    #pragma omp simd aligned(g_y_y_0_0_0_zz_yz_xx, g_y_y_0_0_0_zz_yz_xy, g_y_y_0_0_0_zz_yz_xz, g_y_y_0_0_0_zz_yz_yy, g_y_y_0_0_0_zz_yz_yz, g_y_y_0_0_0_zz_yz_zz, g_y_yzz_yz_xx, g_y_yzz_yz_xy, g_y_yzz_yz_xz, g_y_yzz_yz_yy, g_y_yzz_yz_yz, g_y_yzz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_zz_yz_xx[i] = 4.0 * g_y_yzz_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_yz_xy[i] = 4.0 * g_y_yzz_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_yz_xz[i] = 4.0 * g_y_yzz_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_yz_yy[i] = 4.0 * g_y_yzz_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_yz_yz[i] = 4.0 * g_y_yzz_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_yz_zz[i] = 4.0 * g_y_yzz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1074-1080)

    #pragma omp simd aligned(g_y_y_0_0_0_zz_zz_xx, g_y_y_0_0_0_zz_zz_xy, g_y_y_0_0_0_zz_zz_xz, g_y_y_0_0_0_zz_zz_yy, g_y_y_0_0_0_zz_zz_yz, g_y_y_0_0_0_zz_zz_zz, g_y_yzz_zz_xx, g_y_yzz_zz_xy, g_y_yzz_zz_xz, g_y_yzz_zz_yy, g_y_yzz_zz_yz, g_y_yzz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_zz_zz_xx[i] = 4.0 * g_y_yzz_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_zz_xy[i] = 4.0 * g_y_yzz_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_zz_xz[i] = 4.0 * g_y_yzz_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_zz_yy[i] = 4.0 * g_y_yzz_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_zz_yz[i] = 4.0 * g_y_yzz_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_zz_zz[i] = 4.0 * g_y_yzz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1080-1086)

    #pragma omp simd aligned(g_y_xxz_xx_xx, g_y_xxz_xx_xy, g_y_xxz_xx_xz, g_y_xxz_xx_yy, g_y_xxz_xx_yz, g_y_xxz_xx_zz, g_y_z_0_0_0_xx_xx_xx, g_y_z_0_0_0_xx_xx_xy, g_y_z_0_0_0_xx_xx_xz, g_y_z_0_0_0_xx_xx_yy, g_y_z_0_0_0_xx_xx_yz, g_y_z_0_0_0_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xx_xx_xx[i] = 4.0 * g_y_xxz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_xx_xy[i] = 4.0 * g_y_xxz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_xx_xz[i] = 4.0 * g_y_xxz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_xx_yy[i] = 4.0 * g_y_xxz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_xx_yz[i] = 4.0 * g_y_xxz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_xx_zz[i] = 4.0 * g_y_xxz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1086-1092)

    #pragma omp simd aligned(g_y_xxz_xy_xx, g_y_xxz_xy_xy, g_y_xxz_xy_xz, g_y_xxz_xy_yy, g_y_xxz_xy_yz, g_y_xxz_xy_zz, g_y_z_0_0_0_xx_xy_xx, g_y_z_0_0_0_xx_xy_xy, g_y_z_0_0_0_xx_xy_xz, g_y_z_0_0_0_xx_xy_yy, g_y_z_0_0_0_xx_xy_yz, g_y_z_0_0_0_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xx_xy_xx[i] = 4.0 * g_y_xxz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_xy_xy[i] = 4.0 * g_y_xxz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_xy_xz[i] = 4.0 * g_y_xxz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_xy_yy[i] = 4.0 * g_y_xxz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_xy_yz[i] = 4.0 * g_y_xxz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_xy_zz[i] = 4.0 * g_y_xxz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1092-1098)

    #pragma omp simd aligned(g_y_xxz_xz_xx, g_y_xxz_xz_xy, g_y_xxz_xz_xz, g_y_xxz_xz_yy, g_y_xxz_xz_yz, g_y_xxz_xz_zz, g_y_z_0_0_0_xx_xz_xx, g_y_z_0_0_0_xx_xz_xy, g_y_z_0_0_0_xx_xz_xz, g_y_z_0_0_0_xx_xz_yy, g_y_z_0_0_0_xx_xz_yz, g_y_z_0_0_0_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xx_xz_xx[i] = 4.0 * g_y_xxz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_xz_xy[i] = 4.0 * g_y_xxz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_xz_xz[i] = 4.0 * g_y_xxz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_xz_yy[i] = 4.0 * g_y_xxz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_xz_yz[i] = 4.0 * g_y_xxz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_xz_zz[i] = 4.0 * g_y_xxz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1098-1104)

    #pragma omp simd aligned(g_y_xxz_yy_xx, g_y_xxz_yy_xy, g_y_xxz_yy_xz, g_y_xxz_yy_yy, g_y_xxz_yy_yz, g_y_xxz_yy_zz, g_y_z_0_0_0_xx_yy_xx, g_y_z_0_0_0_xx_yy_xy, g_y_z_0_0_0_xx_yy_xz, g_y_z_0_0_0_xx_yy_yy, g_y_z_0_0_0_xx_yy_yz, g_y_z_0_0_0_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xx_yy_xx[i] = 4.0 * g_y_xxz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_yy_xy[i] = 4.0 * g_y_xxz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_yy_xz[i] = 4.0 * g_y_xxz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_yy_yy[i] = 4.0 * g_y_xxz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_yy_yz[i] = 4.0 * g_y_xxz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_yy_zz[i] = 4.0 * g_y_xxz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1104-1110)

    #pragma omp simd aligned(g_y_xxz_yz_xx, g_y_xxz_yz_xy, g_y_xxz_yz_xz, g_y_xxz_yz_yy, g_y_xxz_yz_yz, g_y_xxz_yz_zz, g_y_z_0_0_0_xx_yz_xx, g_y_z_0_0_0_xx_yz_xy, g_y_z_0_0_0_xx_yz_xz, g_y_z_0_0_0_xx_yz_yy, g_y_z_0_0_0_xx_yz_yz, g_y_z_0_0_0_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xx_yz_xx[i] = 4.0 * g_y_xxz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_yz_xy[i] = 4.0 * g_y_xxz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_yz_xz[i] = 4.0 * g_y_xxz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_yz_yy[i] = 4.0 * g_y_xxz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_yz_yz[i] = 4.0 * g_y_xxz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_yz_zz[i] = 4.0 * g_y_xxz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1110-1116)

    #pragma omp simd aligned(g_y_xxz_zz_xx, g_y_xxz_zz_xy, g_y_xxz_zz_xz, g_y_xxz_zz_yy, g_y_xxz_zz_yz, g_y_xxz_zz_zz, g_y_z_0_0_0_xx_zz_xx, g_y_z_0_0_0_xx_zz_xy, g_y_z_0_0_0_xx_zz_xz, g_y_z_0_0_0_xx_zz_yy, g_y_z_0_0_0_xx_zz_yz, g_y_z_0_0_0_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xx_zz_xx[i] = 4.0 * g_y_xxz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_zz_xy[i] = 4.0 * g_y_xxz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_zz_xz[i] = 4.0 * g_y_xxz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_zz_yy[i] = 4.0 * g_y_xxz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_zz_yz[i] = 4.0 * g_y_xxz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_zz_zz[i] = 4.0 * g_y_xxz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1116-1122)

    #pragma omp simd aligned(g_y_xyz_xx_xx, g_y_xyz_xx_xy, g_y_xyz_xx_xz, g_y_xyz_xx_yy, g_y_xyz_xx_yz, g_y_xyz_xx_zz, g_y_z_0_0_0_xy_xx_xx, g_y_z_0_0_0_xy_xx_xy, g_y_z_0_0_0_xy_xx_xz, g_y_z_0_0_0_xy_xx_yy, g_y_z_0_0_0_xy_xx_yz, g_y_z_0_0_0_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xy_xx_xx[i] = 4.0 * g_y_xyz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_xx_xy[i] = 4.0 * g_y_xyz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_xx_xz[i] = 4.0 * g_y_xyz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_xx_yy[i] = 4.0 * g_y_xyz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_xx_yz[i] = 4.0 * g_y_xyz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_xx_zz[i] = 4.0 * g_y_xyz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1122-1128)

    #pragma omp simd aligned(g_y_xyz_xy_xx, g_y_xyz_xy_xy, g_y_xyz_xy_xz, g_y_xyz_xy_yy, g_y_xyz_xy_yz, g_y_xyz_xy_zz, g_y_z_0_0_0_xy_xy_xx, g_y_z_0_0_0_xy_xy_xy, g_y_z_0_0_0_xy_xy_xz, g_y_z_0_0_0_xy_xy_yy, g_y_z_0_0_0_xy_xy_yz, g_y_z_0_0_0_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xy_xy_xx[i] = 4.0 * g_y_xyz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_xy_xy[i] = 4.0 * g_y_xyz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_xy_xz[i] = 4.0 * g_y_xyz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_xy_yy[i] = 4.0 * g_y_xyz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_xy_yz[i] = 4.0 * g_y_xyz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_xy_zz[i] = 4.0 * g_y_xyz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1128-1134)

    #pragma omp simd aligned(g_y_xyz_xz_xx, g_y_xyz_xz_xy, g_y_xyz_xz_xz, g_y_xyz_xz_yy, g_y_xyz_xz_yz, g_y_xyz_xz_zz, g_y_z_0_0_0_xy_xz_xx, g_y_z_0_0_0_xy_xz_xy, g_y_z_0_0_0_xy_xz_xz, g_y_z_0_0_0_xy_xz_yy, g_y_z_0_0_0_xy_xz_yz, g_y_z_0_0_0_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xy_xz_xx[i] = 4.0 * g_y_xyz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_xz_xy[i] = 4.0 * g_y_xyz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_xz_xz[i] = 4.0 * g_y_xyz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_xz_yy[i] = 4.0 * g_y_xyz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_xz_yz[i] = 4.0 * g_y_xyz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_xz_zz[i] = 4.0 * g_y_xyz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1134-1140)

    #pragma omp simd aligned(g_y_xyz_yy_xx, g_y_xyz_yy_xy, g_y_xyz_yy_xz, g_y_xyz_yy_yy, g_y_xyz_yy_yz, g_y_xyz_yy_zz, g_y_z_0_0_0_xy_yy_xx, g_y_z_0_0_0_xy_yy_xy, g_y_z_0_0_0_xy_yy_xz, g_y_z_0_0_0_xy_yy_yy, g_y_z_0_0_0_xy_yy_yz, g_y_z_0_0_0_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xy_yy_xx[i] = 4.0 * g_y_xyz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_yy_xy[i] = 4.0 * g_y_xyz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_yy_xz[i] = 4.0 * g_y_xyz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_yy_yy[i] = 4.0 * g_y_xyz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_yy_yz[i] = 4.0 * g_y_xyz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_yy_zz[i] = 4.0 * g_y_xyz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1140-1146)

    #pragma omp simd aligned(g_y_xyz_yz_xx, g_y_xyz_yz_xy, g_y_xyz_yz_xz, g_y_xyz_yz_yy, g_y_xyz_yz_yz, g_y_xyz_yz_zz, g_y_z_0_0_0_xy_yz_xx, g_y_z_0_0_0_xy_yz_xy, g_y_z_0_0_0_xy_yz_xz, g_y_z_0_0_0_xy_yz_yy, g_y_z_0_0_0_xy_yz_yz, g_y_z_0_0_0_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xy_yz_xx[i] = 4.0 * g_y_xyz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_yz_xy[i] = 4.0 * g_y_xyz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_yz_xz[i] = 4.0 * g_y_xyz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_yz_yy[i] = 4.0 * g_y_xyz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_yz_yz[i] = 4.0 * g_y_xyz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_yz_zz[i] = 4.0 * g_y_xyz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1146-1152)

    #pragma omp simd aligned(g_y_xyz_zz_xx, g_y_xyz_zz_xy, g_y_xyz_zz_xz, g_y_xyz_zz_yy, g_y_xyz_zz_yz, g_y_xyz_zz_zz, g_y_z_0_0_0_xy_zz_xx, g_y_z_0_0_0_xy_zz_xy, g_y_z_0_0_0_xy_zz_xz, g_y_z_0_0_0_xy_zz_yy, g_y_z_0_0_0_xy_zz_yz, g_y_z_0_0_0_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xy_zz_xx[i] = 4.0 * g_y_xyz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_zz_xy[i] = 4.0 * g_y_xyz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_zz_xz[i] = 4.0 * g_y_xyz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_zz_yy[i] = 4.0 * g_y_xyz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_zz_yz[i] = 4.0 * g_y_xyz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_zz_zz[i] = 4.0 * g_y_xyz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1152-1158)

    #pragma omp simd aligned(g_y_x_xx_xx, g_y_x_xx_xy, g_y_x_xx_xz, g_y_x_xx_yy, g_y_x_xx_yz, g_y_x_xx_zz, g_y_xzz_xx_xx, g_y_xzz_xx_xy, g_y_xzz_xx_xz, g_y_xzz_xx_yy, g_y_xzz_xx_yz, g_y_xzz_xx_zz, g_y_z_0_0_0_xz_xx_xx, g_y_z_0_0_0_xz_xx_xy, g_y_z_0_0_0_xz_xx_xz, g_y_z_0_0_0_xz_xx_yy, g_y_z_0_0_0_xz_xx_yz, g_y_z_0_0_0_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xz_xx_xx[i] = -2.0 * g_y_x_xx_xx[i] * a_exp + 4.0 * g_y_xzz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_xx_xy[i] = -2.0 * g_y_x_xx_xy[i] * a_exp + 4.0 * g_y_xzz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_xx_xz[i] = -2.0 * g_y_x_xx_xz[i] * a_exp + 4.0 * g_y_xzz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_xx_yy[i] = -2.0 * g_y_x_xx_yy[i] * a_exp + 4.0 * g_y_xzz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_xx_yz[i] = -2.0 * g_y_x_xx_yz[i] * a_exp + 4.0 * g_y_xzz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_xx_zz[i] = -2.0 * g_y_x_xx_zz[i] * a_exp + 4.0 * g_y_xzz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1158-1164)

    #pragma omp simd aligned(g_y_x_xy_xx, g_y_x_xy_xy, g_y_x_xy_xz, g_y_x_xy_yy, g_y_x_xy_yz, g_y_x_xy_zz, g_y_xzz_xy_xx, g_y_xzz_xy_xy, g_y_xzz_xy_xz, g_y_xzz_xy_yy, g_y_xzz_xy_yz, g_y_xzz_xy_zz, g_y_z_0_0_0_xz_xy_xx, g_y_z_0_0_0_xz_xy_xy, g_y_z_0_0_0_xz_xy_xz, g_y_z_0_0_0_xz_xy_yy, g_y_z_0_0_0_xz_xy_yz, g_y_z_0_0_0_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xz_xy_xx[i] = -2.0 * g_y_x_xy_xx[i] * a_exp + 4.0 * g_y_xzz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_xy_xy[i] = -2.0 * g_y_x_xy_xy[i] * a_exp + 4.0 * g_y_xzz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_xy_xz[i] = -2.0 * g_y_x_xy_xz[i] * a_exp + 4.0 * g_y_xzz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_xy_yy[i] = -2.0 * g_y_x_xy_yy[i] * a_exp + 4.0 * g_y_xzz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_xy_yz[i] = -2.0 * g_y_x_xy_yz[i] * a_exp + 4.0 * g_y_xzz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_xy_zz[i] = -2.0 * g_y_x_xy_zz[i] * a_exp + 4.0 * g_y_xzz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1164-1170)

    #pragma omp simd aligned(g_y_x_xz_xx, g_y_x_xz_xy, g_y_x_xz_xz, g_y_x_xz_yy, g_y_x_xz_yz, g_y_x_xz_zz, g_y_xzz_xz_xx, g_y_xzz_xz_xy, g_y_xzz_xz_xz, g_y_xzz_xz_yy, g_y_xzz_xz_yz, g_y_xzz_xz_zz, g_y_z_0_0_0_xz_xz_xx, g_y_z_0_0_0_xz_xz_xy, g_y_z_0_0_0_xz_xz_xz, g_y_z_0_0_0_xz_xz_yy, g_y_z_0_0_0_xz_xz_yz, g_y_z_0_0_0_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xz_xz_xx[i] = -2.0 * g_y_x_xz_xx[i] * a_exp + 4.0 * g_y_xzz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_xz_xy[i] = -2.0 * g_y_x_xz_xy[i] * a_exp + 4.0 * g_y_xzz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_xz_xz[i] = -2.0 * g_y_x_xz_xz[i] * a_exp + 4.0 * g_y_xzz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_xz_yy[i] = -2.0 * g_y_x_xz_yy[i] * a_exp + 4.0 * g_y_xzz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_xz_yz[i] = -2.0 * g_y_x_xz_yz[i] * a_exp + 4.0 * g_y_xzz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_xz_zz[i] = -2.0 * g_y_x_xz_zz[i] * a_exp + 4.0 * g_y_xzz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1170-1176)

    #pragma omp simd aligned(g_y_x_yy_xx, g_y_x_yy_xy, g_y_x_yy_xz, g_y_x_yy_yy, g_y_x_yy_yz, g_y_x_yy_zz, g_y_xzz_yy_xx, g_y_xzz_yy_xy, g_y_xzz_yy_xz, g_y_xzz_yy_yy, g_y_xzz_yy_yz, g_y_xzz_yy_zz, g_y_z_0_0_0_xz_yy_xx, g_y_z_0_0_0_xz_yy_xy, g_y_z_0_0_0_xz_yy_xz, g_y_z_0_0_0_xz_yy_yy, g_y_z_0_0_0_xz_yy_yz, g_y_z_0_0_0_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xz_yy_xx[i] = -2.0 * g_y_x_yy_xx[i] * a_exp + 4.0 * g_y_xzz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_yy_xy[i] = -2.0 * g_y_x_yy_xy[i] * a_exp + 4.0 * g_y_xzz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_yy_xz[i] = -2.0 * g_y_x_yy_xz[i] * a_exp + 4.0 * g_y_xzz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_yy_yy[i] = -2.0 * g_y_x_yy_yy[i] * a_exp + 4.0 * g_y_xzz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_yy_yz[i] = -2.0 * g_y_x_yy_yz[i] * a_exp + 4.0 * g_y_xzz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_yy_zz[i] = -2.0 * g_y_x_yy_zz[i] * a_exp + 4.0 * g_y_xzz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1176-1182)

    #pragma omp simd aligned(g_y_x_yz_xx, g_y_x_yz_xy, g_y_x_yz_xz, g_y_x_yz_yy, g_y_x_yz_yz, g_y_x_yz_zz, g_y_xzz_yz_xx, g_y_xzz_yz_xy, g_y_xzz_yz_xz, g_y_xzz_yz_yy, g_y_xzz_yz_yz, g_y_xzz_yz_zz, g_y_z_0_0_0_xz_yz_xx, g_y_z_0_0_0_xz_yz_xy, g_y_z_0_0_0_xz_yz_xz, g_y_z_0_0_0_xz_yz_yy, g_y_z_0_0_0_xz_yz_yz, g_y_z_0_0_0_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xz_yz_xx[i] = -2.0 * g_y_x_yz_xx[i] * a_exp + 4.0 * g_y_xzz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_yz_xy[i] = -2.0 * g_y_x_yz_xy[i] * a_exp + 4.0 * g_y_xzz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_yz_xz[i] = -2.0 * g_y_x_yz_xz[i] * a_exp + 4.0 * g_y_xzz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_yz_yy[i] = -2.0 * g_y_x_yz_yy[i] * a_exp + 4.0 * g_y_xzz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_yz_yz[i] = -2.0 * g_y_x_yz_yz[i] * a_exp + 4.0 * g_y_xzz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_yz_zz[i] = -2.0 * g_y_x_yz_zz[i] * a_exp + 4.0 * g_y_xzz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1182-1188)

    #pragma omp simd aligned(g_y_x_zz_xx, g_y_x_zz_xy, g_y_x_zz_xz, g_y_x_zz_yy, g_y_x_zz_yz, g_y_x_zz_zz, g_y_xzz_zz_xx, g_y_xzz_zz_xy, g_y_xzz_zz_xz, g_y_xzz_zz_yy, g_y_xzz_zz_yz, g_y_xzz_zz_zz, g_y_z_0_0_0_xz_zz_xx, g_y_z_0_0_0_xz_zz_xy, g_y_z_0_0_0_xz_zz_xz, g_y_z_0_0_0_xz_zz_yy, g_y_z_0_0_0_xz_zz_yz, g_y_z_0_0_0_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xz_zz_xx[i] = -2.0 * g_y_x_zz_xx[i] * a_exp + 4.0 * g_y_xzz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_zz_xy[i] = -2.0 * g_y_x_zz_xy[i] * a_exp + 4.0 * g_y_xzz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_zz_xz[i] = -2.0 * g_y_x_zz_xz[i] * a_exp + 4.0 * g_y_xzz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_zz_yy[i] = -2.0 * g_y_x_zz_yy[i] * a_exp + 4.0 * g_y_xzz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_zz_yz[i] = -2.0 * g_y_x_zz_yz[i] * a_exp + 4.0 * g_y_xzz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_zz_zz[i] = -2.0 * g_y_x_zz_zz[i] * a_exp + 4.0 * g_y_xzz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1188-1194)

    #pragma omp simd aligned(g_y_yyz_xx_xx, g_y_yyz_xx_xy, g_y_yyz_xx_xz, g_y_yyz_xx_yy, g_y_yyz_xx_yz, g_y_yyz_xx_zz, g_y_z_0_0_0_yy_xx_xx, g_y_z_0_0_0_yy_xx_xy, g_y_z_0_0_0_yy_xx_xz, g_y_z_0_0_0_yy_xx_yy, g_y_z_0_0_0_yy_xx_yz, g_y_z_0_0_0_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yy_xx_xx[i] = 4.0 * g_y_yyz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_xx_xy[i] = 4.0 * g_y_yyz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_xx_xz[i] = 4.0 * g_y_yyz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_xx_yy[i] = 4.0 * g_y_yyz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_xx_yz[i] = 4.0 * g_y_yyz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_xx_zz[i] = 4.0 * g_y_yyz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1194-1200)

    #pragma omp simd aligned(g_y_yyz_xy_xx, g_y_yyz_xy_xy, g_y_yyz_xy_xz, g_y_yyz_xy_yy, g_y_yyz_xy_yz, g_y_yyz_xy_zz, g_y_z_0_0_0_yy_xy_xx, g_y_z_0_0_0_yy_xy_xy, g_y_z_0_0_0_yy_xy_xz, g_y_z_0_0_0_yy_xy_yy, g_y_z_0_0_0_yy_xy_yz, g_y_z_0_0_0_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yy_xy_xx[i] = 4.0 * g_y_yyz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_xy_xy[i] = 4.0 * g_y_yyz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_xy_xz[i] = 4.0 * g_y_yyz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_xy_yy[i] = 4.0 * g_y_yyz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_xy_yz[i] = 4.0 * g_y_yyz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_xy_zz[i] = 4.0 * g_y_yyz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1200-1206)

    #pragma omp simd aligned(g_y_yyz_xz_xx, g_y_yyz_xz_xy, g_y_yyz_xz_xz, g_y_yyz_xz_yy, g_y_yyz_xz_yz, g_y_yyz_xz_zz, g_y_z_0_0_0_yy_xz_xx, g_y_z_0_0_0_yy_xz_xy, g_y_z_0_0_0_yy_xz_xz, g_y_z_0_0_0_yy_xz_yy, g_y_z_0_0_0_yy_xz_yz, g_y_z_0_0_0_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yy_xz_xx[i] = 4.0 * g_y_yyz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_xz_xy[i] = 4.0 * g_y_yyz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_xz_xz[i] = 4.0 * g_y_yyz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_xz_yy[i] = 4.0 * g_y_yyz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_xz_yz[i] = 4.0 * g_y_yyz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_xz_zz[i] = 4.0 * g_y_yyz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1206-1212)

    #pragma omp simd aligned(g_y_yyz_yy_xx, g_y_yyz_yy_xy, g_y_yyz_yy_xz, g_y_yyz_yy_yy, g_y_yyz_yy_yz, g_y_yyz_yy_zz, g_y_z_0_0_0_yy_yy_xx, g_y_z_0_0_0_yy_yy_xy, g_y_z_0_0_0_yy_yy_xz, g_y_z_0_0_0_yy_yy_yy, g_y_z_0_0_0_yy_yy_yz, g_y_z_0_0_0_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yy_yy_xx[i] = 4.0 * g_y_yyz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_yy_xy[i] = 4.0 * g_y_yyz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_yy_xz[i] = 4.0 * g_y_yyz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_yy_yy[i] = 4.0 * g_y_yyz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_yy_yz[i] = 4.0 * g_y_yyz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_yy_zz[i] = 4.0 * g_y_yyz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1212-1218)

    #pragma omp simd aligned(g_y_yyz_yz_xx, g_y_yyz_yz_xy, g_y_yyz_yz_xz, g_y_yyz_yz_yy, g_y_yyz_yz_yz, g_y_yyz_yz_zz, g_y_z_0_0_0_yy_yz_xx, g_y_z_0_0_0_yy_yz_xy, g_y_z_0_0_0_yy_yz_xz, g_y_z_0_0_0_yy_yz_yy, g_y_z_0_0_0_yy_yz_yz, g_y_z_0_0_0_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yy_yz_xx[i] = 4.0 * g_y_yyz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_yz_xy[i] = 4.0 * g_y_yyz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_yz_xz[i] = 4.0 * g_y_yyz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_yz_yy[i] = 4.0 * g_y_yyz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_yz_yz[i] = 4.0 * g_y_yyz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_yz_zz[i] = 4.0 * g_y_yyz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1218-1224)

    #pragma omp simd aligned(g_y_yyz_zz_xx, g_y_yyz_zz_xy, g_y_yyz_zz_xz, g_y_yyz_zz_yy, g_y_yyz_zz_yz, g_y_yyz_zz_zz, g_y_z_0_0_0_yy_zz_xx, g_y_z_0_0_0_yy_zz_xy, g_y_z_0_0_0_yy_zz_xz, g_y_z_0_0_0_yy_zz_yy, g_y_z_0_0_0_yy_zz_yz, g_y_z_0_0_0_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yy_zz_xx[i] = 4.0 * g_y_yyz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_zz_xy[i] = 4.0 * g_y_yyz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_zz_xz[i] = 4.0 * g_y_yyz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_zz_yy[i] = 4.0 * g_y_yyz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_zz_yz[i] = 4.0 * g_y_yyz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_zz_zz[i] = 4.0 * g_y_yyz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1224-1230)

    #pragma omp simd aligned(g_y_y_xx_xx, g_y_y_xx_xy, g_y_y_xx_xz, g_y_y_xx_yy, g_y_y_xx_yz, g_y_y_xx_zz, g_y_yzz_xx_xx, g_y_yzz_xx_xy, g_y_yzz_xx_xz, g_y_yzz_xx_yy, g_y_yzz_xx_yz, g_y_yzz_xx_zz, g_y_z_0_0_0_yz_xx_xx, g_y_z_0_0_0_yz_xx_xy, g_y_z_0_0_0_yz_xx_xz, g_y_z_0_0_0_yz_xx_yy, g_y_z_0_0_0_yz_xx_yz, g_y_z_0_0_0_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yz_xx_xx[i] = -2.0 * g_y_y_xx_xx[i] * a_exp + 4.0 * g_y_yzz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_xx_xy[i] = -2.0 * g_y_y_xx_xy[i] * a_exp + 4.0 * g_y_yzz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_xx_xz[i] = -2.0 * g_y_y_xx_xz[i] * a_exp + 4.0 * g_y_yzz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_xx_yy[i] = -2.0 * g_y_y_xx_yy[i] * a_exp + 4.0 * g_y_yzz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_xx_yz[i] = -2.0 * g_y_y_xx_yz[i] * a_exp + 4.0 * g_y_yzz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_xx_zz[i] = -2.0 * g_y_y_xx_zz[i] * a_exp + 4.0 * g_y_yzz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1230-1236)

    #pragma omp simd aligned(g_y_y_xy_xx, g_y_y_xy_xy, g_y_y_xy_xz, g_y_y_xy_yy, g_y_y_xy_yz, g_y_y_xy_zz, g_y_yzz_xy_xx, g_y_yzz_xy_xy, g_y_yzz_xy_xz, g_y_yzz_xy_yy, g_y_yzz_xy_yz, g_y_yzz_xy_zz, g_y_z_0_0_0_yz_xy_xx, g_y_z_0_0_0_yz_xy_xy, g_y_z_0_0_0_yz_xy_xz, g_y_z_0_0_0_yz_xy_yy, g_y_z_0_0_0_yz_xy_yz, g_y_z_0_0_0_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yz_xy_xx[i] = -2.0 * g_y_y_xy_xx[i] * a_exp + 4.0 * g_y_yzz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_xy_xy[i] = -2.0 * g_y_y_xy_xy[i] * a_exp + 4.0 * g_y_yzz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_xy_xz[i] = -2.0 * g_y_y_xy_xz[i] * a_exp + 4.0 * g_y_yzz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_xy_yy[i] = -2.0 * g_y_y_xy_yy[i] * a_exp + 4.0 * g_y_yzz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_xy_yz[i] = -2.0 * g_y_y_xy_yz[i] * a_exp + 4.0 * g_y_yzz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_xy_zz[i] = -2.0 * g_y_y_xy_zz[i] * a_exp + 4.0 * g_y_yzz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1236-1242)

    #pragma omp simd aligned(g_y_y_xz_xx, g_y_y_xz_xy, g_y_y_xz_xz, g_y_y_xz_yy, g_y_y_xz_yz, g_y_y_xz_zz, g_y_yzz_xz_xx, g_y_yzz_xz_xy, g_y_yzz_xz_xz, g_y_yzz_xz_yy, g_y_yzz_xz_yz, g_y_yzz_xz_zz, g_y_z_0_0_0_yz_xz_xx, g_y_z_0_0_0_yz_xz_xy, g_y_z_0_0_0_yz_xz_xz, g_y_z_0_0_0_yz_xz_yy, g_y_z_0_0_0_yz_xz_yz, g_y_z_0_0_0_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yz_xz_xx[i] = -2.0 * g_y_y_xz_xx[i] * a_exp + 4.0 * g_y_yzz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_xz_xy[i] = -2.0 * g_y_y_xz_xy[i] * a_exp + 4.0 * g_y_yzz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_xz_xz[i] = -2.0 * g_y_y_xz_xz[i] * a_exp + 4.0 * g_y_yzz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_xz_yy[i] = -2.0 * g_y_y_xz_yy[i] * a_exp + 4.0 * g_y_yzz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_xz_yz[i] = -2.0 * g_y_y_xz_yz[i] * a_exp + 4.0 * g_y_yzz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_xz_zz[i] = -2.0 * g_y_y_xz_zz[i] * a_exp + 4.0 * g_y_yzz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1242-1248)

    #pragma omp simd aligned(g_y_y_yy_xx, g_y_y_yy_xy, g_y_y_yy_xz, g_y_y_yy_yy, g_y_y_yy_yz, g_y_y_yy_zz, g_y_yzz_yy_xx, g_y_yzz_yy_xy, g_y_yzz_yy_xz, g_y_yzz_yy_yy, g_y_yzz_yy_yz, g_y_yzz_yy_zz, g_y_z_0_0_0_yz_yy_xx, g_y_z_0_0_0_yz_yy_xy, g_y_z_0_0_0_yz_yy_xz, g_y_z_0_0_0_yz_yy_yy, g_y_z_0_0_0_yz_yy_yz, g_y_z_0_0_0_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yz_yy_xx[i] = -2.0 * g_y_y_yy_xx[i] * a_exp + 4.0 * g_y_yzz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_yy_xy[i] = -2.0 * g_y_y_yy_xy[i] * a_exp + 4.0 * g_y_yzz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_yy_xz[i] = -2.0 * g_y_y_yy_xz[i] * a_exp + 4.0 * g_y_yzz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_yy_yy[i] = -2.0 * g_y_y_yy_yy[i] * a_exp + 4.0 * g_y_yzz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_yy_yz[i] = -2.0 * g_y_y_yy_yz[i] * a_exp + 4.0 * g_y_yzz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_yy_zz[i] = -2.0 * g_y_y_yy_zz[i] * a_exp + 4.0 * g_y_yzz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1248-1254)

    #pragma omp simd aligned(g_y_y_yz_xx, g_y_y_yz_xy, g_y_y_yz_xz, g_y_y_yz_yy, g_y_y_yz_yz, g_y_y_yz_zz, g_y_yzz_yz_xx, g_y_yzz_yz_xy, g_y_yzz_yz_xz, g_y_yzz_yz_yy, g_y_yzz_yz_yz, g_y_yzz_yz_zz, g_y_z_0_0_0_yz_yz_xx, g_y_z_0_0_0_yz_yz_xy, g_y_z_0_0_0_yz_yz_xz, g_y_z_0_0_0_yz_yz_yy, g_y_z_0_0_0_yz_yz_yz, g_y_z_0_0_0_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yz_yz_xx[i] = -2.0 * g_y_y_yz_xx[i] * a_exp + 4.0 * g_y_yzz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_yz_xy[i] = -2.0 * g_y_y_yz_xy[i] * a_exp + 4.0 * g_y_yzz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_yz_xz[i] = -2.0 * g_y_y_yz_xz[i] * a_exp + 4.0 * g_y_yzz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_yz_yy[i] = -2.0 * g_y_y_yz_yy[i] * a_exp + 4.0 * g_y_yzz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_yz_yz[i] = -2.0 * g_y_y_yz_yz[i] * a_exp + 4.0 * g_y_yzz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_yz_zz[i] = -2.0 * g_y_y_yz_zz[i] * a_exp + 4.0 * g_y_yzz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1254-1260)

    #pragma omp simd aligned(g_y_y_zz_xx, g_y_y_zz_xy, g_y_y_zz_xz, g_y_y_zz_yy, g_y_y_zz_yz, g_y_y_zz_zz, g_y_yzz_zz_xx, g_y_yzz_zz_xy, g_y_yzz_zz_xz, g_y_yzz_zz_yy, g_y_yzz_zz_yz, g_y_yzz_zz_zz, g_y_z_0_0_0_yz_zz_xx, g_y_z_0_0_0_yz_zz_xy, g_y_z_0_0_0_yz_zz_xz, g_y_z_0_0_0_yz_zz_yy, g_y_z_0_0_0_yz_zz_yz, g_y_z_0_0_0_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yz_zz_xx[i] = -2.0 * g_y_y_zz_xx[i] * a_exp + 4.0 * g_y_yzz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_zz_xy[i] = -2.0 * g_y_y_zz_xy[i] * a_exp + 4.0 * g_y_yzz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_zz_xz[i] = -2.0 * g_y_y_zz_xz[i] * a_exp + 4.0 * g_y_yzz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_zz_yy[i] = -2.0 * g_y_y_zz_yy[i] * a_exp + 4.0 * g_y_yzz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_zz_yz[i] = -2.0 * g_y_y_zz_yz[i] * a_exp + 4.0 * g_y_yzz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_zz_zz[i] = -2.0 * g_y_y_zz_zz[i] * a_exp + 4.0 * g_y_yzz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1260-1266)

    #pragma omp simd aligned(g_y_z_0_0_0_zz_xx_xx, g_y_z_0_0_0_zz_xx_xy, g_y_z_0_0_0_zz_xx_xz, g_y_z_0_0_0_zz_xx_yy, g_y_z_0_0_0_zz_xx_yz, g_y_z_0_0_0_zz_xx_zz, g_y_z_xx_xx, g_y_z_xx_xy, g_y_z_xx_xz, g_y_z_xx_yy, g_y_z_xx_yz, g_y_z_xx_zz, g_y_zzz_xx_xx, g_y_zzz_xx_xy, g_y_zzz_xx_xz, g_y_zzz_xx_yy, g_y_zzz_xx_yz, g_y_zzz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_zz_xx_xx[i] = -4.0 * g_y_z_xx_xx[i] * a_exp + 4.0 * g_y_zzz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_xx_xy[i] = -4.0 * g_y_z_xx_xy[i] * a_exp + 4.0 * g_y_zzz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_xx_xz[i] = -4.0 * g_y_z_xx_xz[i] * a_exp + 4.0 * g_y_zzz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_xx_yy[i] = -4.0 * g_y_z_xx_yy[i] * a_exp + 4.0 * g_y_zzz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_xx_yz[i] = -4.0 * g_y_z_xx_yz[i] * a_exp + 4.0 * g_y_zzz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_xx_zz[i] = -4.0 * g_y_z_xx_zz[i] * a_exp + 4.0 * g_y_zzz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1266-1272)

    #pragma omp simd aligned(g_y_z_0_0_0_zz_xy_xx, g_y_z_0_0_0_zz_xy_xy, g_y_z_0_0_0_zz_xy_xz, g_y_z_0_0_0_zz_xy_yy, g_y_z_0_0_0_zz_xy_yz, g_y_z_0_0_0_zz_xy_zz, g_y_z_xy_xx, g_y_z_xy_xy, g_y_z_xy_xz, g_y_z_xy_yy, g_y_z_xy_yz, g_y_z_xy_zz, g_y_zzz_xy_xx, g_y_zzz_xy_xy, g_y_zzz_xy_xz, g_y_zzz_xy_yy, g_y_zzz_xy_yz, g_y_zzz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_zz_xy_xx[i] = -4.0 * g_y_z_xy_xx[i] * a_exp + 4.0 * g_y_zzz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_xy_xy[i] = -4.0 * g_y_z_xy_xy[i] * a_exp + 4.0 * g_y_zzz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_xy_xz[i] = -4.0 * g_y_z_xy_xz[i] * a_exp + 4.0 * g_y_zzz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_xy_yy[i] = -4.0 * g_y_z_xy_yy[i] * a_exp + 4.0 * g_y_zzz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_xy_yz[i] = -4.0 * g_y_z_xy_yz[i] * a_exp + 4.0 * g_y_zzz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_xy_zz[i] = -4.0 * g_y_z_xy_zz[i] * a_exp + 4.0 * g_y_zzz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1272-1278)

    #pragma omp simd aligned(g_y_z_0_0_0_zz_xz_xx, g_y_z_0_0_0_zz_xz_xy, g_y_z_0_0_0_zz_xz_xz, g_y_z_0_0_0_zz_xz_yy, g_y_z_0_0_0_zz_xz_yz, g_y_z_0_0_0_zz_xz_zz, g_y_z_xz_xx, g_y_z_xz_xy, g_y_z_xz_xz, g_y_z_xz_yy, g_y_z_xz_yz, g_y_z_xz_zz, g_y_zzz_xz_xx, g_y_zzz_xz_xy, g_y_zzz_xz_xz, g_y_zzz_xz_yy, g_y_zzz_xz_yz, g_y_zzz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_zz_xz_xx[i] = -4.0 * g_y_z_xz_xx[i] * a_exp + 4.0 * g_y_zzz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_xz_xy[i] = -4.0 * g_y_z_xz_xy[i] * a_exp + 4.0 * g_y_zzz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_xz_xz[i] = -4.0 * g_y_z_xz_xz[i] * a_exp + 4.0 * g_y_zzz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_xz_yy[i] = -4.0 * g_y_z_xz_yy[i] * a_exp + 4.0 * g_y_zzz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_xz_yz[i] = -4.0 * g_y_z_xz_yz[i] * a_exp + 4.0 * g_y_zzz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_xz_zz[i] = -4.0 * g_y_z_xz_zz[i] * a_exp + 4.0 * g_y_zzz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1278-1284)

    #pragma omp simd aligned(g_y_z_0_0_0_zz_yy_xx, g_y_z_0_0_0_zz_yy_xy, g_y_z_0_0_0_zz_yy_xz, g_y_z_0_0_0_zz_yy_yy, g_y_z_0_0_0_zz_yy_yz, g_y_z_0_0_0_zz_yy_zz, g_y_z_yy_xx, g_y_z_yy_xy, g_y_z_yy_xz, g_y_z_yy_yy, g_y_z_yy_yz, g_y_z_yy_zz, g_y_zzz_yy_xx, g_y_zzz_yy_xy, g_y_zzz_yy_xz, g_y_zzz_yy_yy, g_y_zzz_yy_yz, g_y_zzz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_zz_yy_xx[i] = -4.0 * g_y_z_yy_xx[i] * a_exp + 4.0 * g_y_zzz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_yy_xy[i] = -4.0 * g_y_z_yy_xy[i] * a_exp + 4.0 * g_y_zzz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_yy_xz[i] = -4.0 * g_y_z_yy_xz[i] * a_exp + 4.0 * g_y_zzz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_yy_yy[i] = -4.0 * g_y_z_yy_yy[i] * a_exp + 4.0 * g_y_zzz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_yy_yz[i] = -4.0 * g_y_z_yy_yz[i] * a_exp + 4.0 * g_y_zzz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_yy_zz[i] = -4.0 * g_y_z_yy_zz[i] * a_exp + 4.0 * g_y_zzz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1284-1290)

    #pragma omp simd aligned(g_y_z_0_0_0_zz_yz_xx, g_y_z_0_0_0_zz_yz_xy, g_y_z_0_0_0_zz_yz_xz, g_y_z_0_0_0_zz_yz_yy, g_y_z_0_0_0_zz_yz_yz, g_y_z_0_0_0_zz_yz_zz, g_y_z_yz_xx, g_y_z_yz_xy, g_y_z_yz_xz, g_y_z_yz_yy, g_y_z_yz_yz, g_y_z_yz_zz, g_y_zzz_yz_xx, g_y_zzz_yz_xy, g_y_zzz_yz_xz, g_y_zzz_yz_yy, g_y_zzz_yz_yz, g_y_zzz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_zz_yz_xx[i] = -4.0 * g_y_z_yz_xx[i] * a_exp + 4.0 * g_y_zzz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_yz_xy[i] = -4.0 * g_y_z_yz_xy[i] * a_exp + 4.0 * g_y_zzz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_yz_xz[i] = -4.0 * g_y_z_yz_xz[i] * a_exp + 4.0 * g_y_zzz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_yz_yy[i] = -4.0 * g_y_z_yz_yy[i] * a_exp + 4.0 * g_y_zzz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_yz_yz[i] = -4.0 * g_y_z_yz_yz[i] * a_exp + 4.0 * g_y_zzz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_yz_zz[i] = -4.0 * g_y_z_yz_zz[i] * a_exp + 4.0 * g_y_zzz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1290-1296)

    #pragma omp simd aligned(g_y_z_0_0_0_zz_zz_xx, g_y_z_0_0_0_zz_zz_xy, g_y_z_0_0_0_zz_zz_xz, g_y_z_0_0_0_zz_zz_yy, g_y_z_0_0_0_zz_zz_yz, g_y_z_0_0_0_zz_zz_zz, g_y_z_zz_xx, g_y_z_zz_xy, g_y_z_zz_xz, g_y_z_zz_yy, g_y_z_zz_yz, g_y_z_zz_zz, g_y_zzz_zz_xx, g_y_zzz_zz_xy, g_y_zzz_zz_xz, g_y_zzz_zz_yy, g_y_zzz_zz_yz, g_y_zzz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_zz_zz_xx[i] = -4.0 * g_y_z_zz_xx[i] * a_exp + 4.0 * g_y_zzz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_zz_xy[i] = -4.0 * g_y_z_zz_xy[i] * a_exp + 4.0 * g_y_zzz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_zz_xz[i] = -4.0 * g_y_z_zz_xz[i] * a_exp + 4.0 * g_y_zzz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_zz_yy[i] = -4.0 * g_y_z_zz_yy[i] * a_exp + 4.0 * g_y_zzz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_zz_yz[i] = -4.0 * g_y_z_zz_yz[i] * a_exp + 4.0 * g_y_zzz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_zz_zz[i] = -4.0 * g_y_z_zz_zz[i] * a_exp + 4.0 * g_y_zzz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1296-1302)

    #pragma omp simd aligned(g_z_x_0_0_0_xx_xx_xx, g_z_x_0_0_0_xx_xx_xy, g_z_x_0_0_0_xx_xx_xz, g_z_x_0_0_0_xx_xx_yy, g_z_x_0_0_0_xx_xx_yz, g_z_x_0_0_0_xx_xx_zz, g_z_x_xx_xx, g_z_x_xx_xy, g_z_x_xx_xz, g_z_x_xx_yy, g_z_x_xx_yz, g_z_x_xx_zz, g_z_xxx_xx_xx, g_z_xxx_xx_xy, g_z_xxx_xx_xz, g_z_xxx_xx_yy, g_z_xxx_xx_yz, g_z_xxx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xx_xx_xx[i] = -4.0 * g_z_x_xx_xx[i] * a_exp + 4.0 * g_z_xxx_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_xx_xy[i] = -4.0 * g_z_x_xx_xy[i] * a_exp + 4.0 * g_z_xxx_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_xx_xz[i] = -4.0 * g_z_x_xx_xz[i] * a_exp + 4.0 * g_z_xxx_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_xx_yy[i] = -4.0 * g_z_x_xx_yy[i] * a_exp + 4.0 * g_z_xxx_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_xx_yz[i] = -4.0 * g_z_x_xx_yz[i] * a_exp + 4.0 * g_z_xxx_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_xx_zz[i] = -4.0 * g_z_x_xx_zz[i] * a_exp + 4.0 * g_z_xxx_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1302-1308)

    #pragma omp simd aligned(g_z_x_0_0_0_xx_xy_xx, g_z_x_0_0_0_xx_xy_xy, g_z_x_0_0_0_xx_xy_xz, g_z_x_0_0_0_xx_xy_yy, g_z_x_0_0_0_xx_xy_yz, g_z_x_0_0_0_xx_xy_zz, g_z_x_xy_xx, g_z_x_xy_xy, g_z_x_xy_xz, g_z_x_xy_yy, g_z_x_xy_yz, g_z_x_xy_zz, g_z_xxx_xy_xx, g_z_xxx_xy_xy, g_z_xxx_xy_xz, g_z_xxx_xy_yy, g_z_xxx_xy_yz, g_z_xxx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xx_xy_xx[i] = -4.0 * g_z_x_xy_xx[i] * a_exp + 4.0 * g_z_xxx_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_xy_xy[i] = -4.0 * g_z_x_xy_xy[i] * a_exp + 4.0 * g_z_xxx_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_xy_xz[i] = -4.0 * g_z_x_xy_xz[i] * a_exp + 4.0 * g_z_xxx_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_xy_yy[i] = -4.0 * g_z_x_xy_yy[i] * a_exp + 4.0 * g_z_xxx_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_xy_yz[i] = -4.0 * g_z_x_xy_yz[i] * a_exp + 4.0 * g_z_xxx_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_xy_zz[i] = -4.0 * g_z_x_xy_zz[i] * a_exp + 4.0 * g_z_xxx_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1308-1314)

    #pragma omp simd aligned(g_z_x_0_0_0_xx_xz_xx, g_z_x_0_0_0_xx_xz_xy, g_z_x_0_0_0_xx_xz_xz, g_z_x_0_0_0_xx_xz_yy, g_z_x_0_0_0_xx_xz_yz, g_z_x_0_0_0_xx_xz_zz, g_z_x_xz_xx, g_z_x_xz_xy, g_z_x_xz_xz, g_z_x_xz_yy, g_z_x_xz_yz, g_z_x_xz_zz, g_z_xxx_xz_xx, g_z_xxx_xz_xy, g_z_xxx_xz_xz, g_z_xxx_xz_yy, g_z_xxx_xz_yz, g_z_xxx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xx_xz_xx[i] = -4.0 * g_z_x_xz_xx[i] * a_exp + 4.0 * g_z_xxx_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_xz_xy[i] = -4.0 * g_z_x_xz_xy[i] * a_exp + 4.0 * g_z_xxx_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_xz_xz[i] = -4.0 * g_z_x_xz_xz[i] * a_exp + 4.0 * g_z_xxx_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_xz_yy[i] = -4.0 * g_z_x_xz_yy[i] * a_exp + 4.0 * g_z_xxx_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_xz_yz[i] = -4.0 * g_z_x_xz_yz[i] * a_exp + 4.0 * g_z_xxx_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_xz_zz[i] = -4.0 * g_z_x_xz_zz[i] * a_exp + 4.0 * g_z_xxx_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1314-1320)

    #pragma omp simd aligned(g_z_x_0_0_0_xx_yy_xx, g_z_x_0_0_0_xx_yy_xy, g_z_x_0_0_0_xx_yy_xz, g_z_x_0_0_0_xx_yy_yy, g_z_x_0_0_0_xx_yy_yz, g_z_x_0_0_0_xx_yy_zz, g_z_x_yy_xx, g_z_x_yy_xy, g_z_x_yy_xz, g_z_x_yy_yy, g_z_x_yy_yz, g_z_x_yy_zz, g_z_xxx_yy_xx, g_z_xxx_yy_xy, g_z_xxx_yy_xz, g_z_xxx_yy_yy, g_z_xxx_yy_yz, g_z_xxx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xx_yy_xx[i] = -4.0 * g_z_x_yy_xx[i] * a_exp + 4.0 * g_z_xxx_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_yy_xy[i] = -4.0 * g_z_x_yy_xy[i] * a_exp + 4.0 * g_z_xxx_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_yy_xz[i] = -4.0 * g_z_x_yy_xz[i] * a_exp + 4.0 * g_z_xxx_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_yy_yy[i] = -4.0 * g_z_x_yy_yy[i] * a_exp + 4.0 * g_z_xxx_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_yy_yz[i] = -4.0 * g_z_x_yy_yz[i] * a_exp + 4.0 * g_z_xxx_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_yy_zz[i] = -4.0 * g_z_x_yy_zz[i] * a_exp + 4.0 * g_z_xxx_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1320-1326)

    #pragma omp simd aligned(g_z_x_0_0_0_xx_yz_xx, g_z_x_0_0_0_xx_yz_xy, g_z_x_0_0_0_xx_yz_xz, g_z_x_0_0_0_xx_yz_yy, g_z_x_0_0_0_xx_yz_yz, g_z_x_0_0_0_xx_yz_zz, g_z_x_yz_xx, g_z_x_yz_xy, g_z_x_yz_xz, g_z_x_yz_yy, g_z_x_yz_yz, g_z_x_yz_zz, g_z_xxx_yz_xx, g_z_xxx_yz_xy, g_z_xxx_yz_xz, g_z_xxx_yz_yy, g_z_xxx_yz_yz, g_z_xxx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xx_yz_xx[i] = -4.0 * g_z_x_yz_xx[i] * a_exp + 4.0 * g_z_xxx_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_yz_xy[i] = -4.0 * g_z_x_yz_xy[i] * a_exp + 4.0 * g_z_xxx_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_yz_xz[i] = -4.0 * g_z_x_yz_xz[i] * a_exp + 4.0 * g_z_xxx_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_yz_yy[i] = -4.0 * g_z_x_yz_yy[i] * a_exp + 4.0 * g_z_xxx_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_yz_yz[i] = -4.0 * g_z_x_yz_yz[i] * a_exp + 4.0 * g_z_xxx_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_yz_zz[i] = -4.0 * g_z_x_yz_zz[i] * a_exp + 4.0 * g_z_xxx_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1326-1332)

    #pragma omp simd aligned(g_z_x_0_0_0_xx_zz_xx, g_z_x_0_0_0_xx_zz_xy, g_z_x_0_0_0_xx_zz_xz, g_z_x_0_0_0_xx_zz_yy, g_z_x_0_0_0_xx_zz_yz, g_z_x_0_0_0_xx_zz_zz, g_z_x_zz_xx, g_z_x_zz_xy, g_z_x_zz_xz, g_z_x_zz_yy, g_z_x_zz_yz, g_z_x_zz_zz, g_z_xxx_zz_xx, g_z_xxx_zz_xy, g_z_xxx_zz_xz, g_z_xxx_zz_yy, g_z_xxx_zz_yz, g_z_xxx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xx_zz_xx[i] = -4.0 * g_z_x_zz_xx[i] * a_exp + 4.0 * g_z_xxx_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_zz_xy[i] = -4.0 * g_z_x_zz_xy[i] * a_exp + 4.0 * g_z_xxx_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_zz_xz[i] = -4.0 * g_z_x_zz_xz[i] * a_exp + 4.0 * g_z_xxx_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_zz_yy[i] = -4.0 * g_z_x_zz_yy[i] * a_exp + 4.0 * g_z_xxx_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_zz_yz[i] = -4.0 * g_z_x_zz_yz[i] * a_exp + 4.0 * g_z_xxx_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_zz_zz[i] = -4.0 * g_z_x_zz_zz[i] * a_exp + 4.0 * g_z_xxx_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1332-1338)

    #pragma omp simd aligned(g_z_x_0_0_0_xy_xx_xx, g_z_x_0_0_0_xy_xx_xy, g_z_x_0_0_0_xy_xx_xz, g_z_x_0_0_0_xy_xx_yy, g_z_x_0_0_0_xy_xx_yz, g_z_x_0_0_0_xy_xx_zz, g_z_xxy_xx_xx, g_z_xxy_xx_xy, g_z_xxy_xx_xz, g_z_xxy_xx_yy, g_z_xxy_xx_yz, g_z_xxy_xx_zz, g_z_y_xx_xx, g_z_y_xx_xy, g_z_y_xx_xz, g_z_y_xx_yy, g_z_y_xx_yz, g_z_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xy_xx_xx[i] = -2.0 * g_z_y_xx_xx[i] * a_exp + 4.0 * g_z_xxy_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_xx_xy[i] = -2.0 * g_z_y_xx_xy[i] * a_exp + 4.0 * g_z_xxy_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_xx_xz[i] = -2.0 * g_z_y_xx_xz[i] * a_exp + 4.0 * g_z_xxy_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_xx_yy[i] = -2.0 * g_z_y_xx_yy[i] * a_exp + 4.0 * g_z_xxy_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_xx_yz[i] = -2.0 * g_z_y_xx_yz[i] * a_exp + 4.0 * g_z_xxy_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_xx_zz[i] = -2.0 * g_z_y_xx_zz[i] * a_exp + 4.0 * g_z_xxy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1338-1344)

    #pragma omp simd aligned(g_z_x_0_0_0_xy_xy_xx, g_z_x_0_0_0_xy_xy_xy, g_z_x_0_0_0_xy_xy_xz, g_z_x_0_0_0_xy_xy_yy, g_z_x_0_0_0_xy_xy_yz, g_z_x_0_0_0_xy_xy_zz, g_z_xxy_xy_xx, g_z_xxy_xy_xy, g_z_xxy_xy_xz, g_z_xxy_xy_yy, g_z_xxy_xy_yz, g_z_xxy_xy_zz, g_z_y_xy_xx, g_z_y_xy_xy, g_z_y_xy_xz, g_z_y_xy_yy, g_z_y_xy_yz, g_z_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xy_xy_xx[i] = -2.0 * g_z_y_xy_xx[i] * a_exp + 4.0 * g_z_xxy_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_xy_xy[i] = -2.0 * g_z_y_xy_xy[i] * a_exp + 4.0 * g_z_xxy_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_xy_xz[i] = -2.0 * g_z_y_xy_xz[i] * a_exp + 4.0 * g_z_xxy_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_xy_yy[i] = -2.0 * g_z_y_xy_yy[i] * a_exp + 4.0 * g_z_xxy_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_xy_yz[i] = -2.0 * g_z_y_xy_yz[i] * a_exp + 4.0 * g_z_xxy_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_xy_zz[i] = -2.0 * g_z_y_xy_zz[i] * a_exp + 4.0 * g_z_xxy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1344-1350)

    #pragma omp simd aligned(g_z_x_0_0_0_xy_xz_xx, g_z_x_0_0_0_xy_xz_xy, g_z_x_0_0_0_xy_xz_xz, g_z_x_0_0_0_xy_xz_yy, g_z_x_0_0_0_xy_xz_yz, g_z_x_0_0_0_xy_xz_zz, g_z_xxy_xz_xx, g_z_xxy_xz_xy, g_z_xxy_xz_xz, g_z_xxy_xz_yy, g_z_xxy_xz_yz, g_z_xxy_xz_zz, g_z_y_xz_xx, g_z_y_xz_xy, g_z_y_xz_xz, g_z_y_xz_yy, g_z_y_xz_yz, g_z_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xy_xz_xx[i] = -2.0 * g_z_y_xz_xx[i] * a_exp + 4.0 * g_z_xxy_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_xz_xy[i] = -2.0 * g_z_y_xz_xy[i] * a_exp + 4.0 * g_z_xxy_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_xz_xz[i] = -2.0 * g_z_y_xz_xz[i] * a_exp + 4.0 * g_z_xxy_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_xz_yy[i] = -2.0 * g_z_y_xz_yy[i] * a_exp + 4.0 * g_z_xxy_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_xz_yz[i] = -2.0 * g_z_y_xz_yz[i] * a_exp + 4.0 * g_z_xxy_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_xz_zz[i] = -2.0 * g_z_y_xz_zz[i] * a_exp + 4.0 * g_z_xxy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1350-1356)

    #pragma omp simd aligned(g_z_x_0_0_0_xy_yy_xx, g_z_x_0_0_0_xy_yy_xy, g_z_x_0_0_0_xy_yy_xz, g_z_x_0_0_0_xy_yy_yy, g_z_x_0_0_0_xy_yy_yz, g_z_x_0_0_0_xy_yy_zz, g_z_xxy_yy_xx, g_z_xxy_yy_xy, g_z_xxy_yy_xz, g_z_xxy_yy_yy, g_z_xxy_yy_yz, g_z_xxy_yy_zz, g_z_y_yy_xx, g_z_y_yy_xy, g_z_y_yy_xz, g_z_y_yy_yy, g_z_y_yy_yz, g_z_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xy_yy_xx[i] = -2.0 * g_z_y_yy_xx[i] * a_exp + 4.0 * g_z_xxy_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_yy_xy[i] = -2.0 * g_z_y_yy_xy[i] * a_exp + 4.0 * g_z_xxy_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_yy_xz[i] = -2.0 * g_z_y_yy_xz[i] * a_exp + 4.0 * g_z_xxy_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_yy_yy[i] = -2.0 * g_z_y_yy_yy[i] * a_exp + 4.0 * g_z_xxy_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_yy_yz[i] = -2.0 * g_z_y_yy_yz[i] * a_exp + 4.0 * g_z_xxy_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_yy_zz[i] = -2.0 * g_z_y_yy_zz[i] * a_exp + 4.0 * g_z_xxy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1356-1362)

    #pragma omp simd aligned(g_z_x_0_0_0_xy_yz_xx, g_z_x_0_0_0_xy_yz_xy, g_z_x_0_0_0_xy_yz_xz, g_z_x_0_0_0_xy_yz_yy, g_z_x_0_0_0_xy_yz_yz, g_z_x_0_0_0_xy_yz_zz, g_z_xxy_yz_xx, g_z_xxy_yz_xy, g_z_xxy_yz_xz, g_z_xxy_yz_yy, g_z_xxy_yz_yz, g_z_xxy_yz_zz, g_z_y_yz_xx, g_z_y_yz_xy, g_z_y_yz_xz, g_z_y_yz_yy, g_z_y_yz_yz, g_z_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xy_yz_xx[i] = -2.0 * g_z_y_yz_xx[i] * a_exp + 4.0 * g_z_xxy_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_yz_xy[i] = -2.0 * g_z_y_yz_xy[i] * a_exp + 4.0 * g_z_xxy_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_yz_xz[i] = -2.0 * g_z_y_yz_xz[i] * a_exp + 4.0 * g_z_xxy_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_yz_yy[i] = -2.0 * g_z_y_yz_yy[i] * a_exp + 4.0 * g_z_xxy_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_yz_yz[i] = -2.0 * g_z_y_yz_yz[i] * a_exp + 4.0 * g_z_xxy_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_yz_zz[i] = -2.0 * g_z_y_yz_zz[i] * a_exp + 4.0 * g_z_xxy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1362-1368)

    #pragma omp simd aligned(g_z_x_0_0_0_xy_zz_xx, g_z_x_0_0_0_xy_zz_xy, g_z_x_0_0_0_xy_zz_xz, g_z_x_0_0_0_xy_zz_yy, g_z_x_0_0_0_xy_zz_yz, g_z_x_0_0_0_xy_zz_zz, g_z_xxy_zz_xx, g_z_xxy_zz_xy, g_z_xxy_zz_xz, g_z_xxy_zz_yy, g_z_xxy_zz_yz, g_z_xxy_zz_zz, g_z_y_zz_xx, g_z_y_zz_xy, g_z_y_zz_xz, g_z_y_zz_yy, g_z_y_zz_yz, g_z_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xy_zz_xx[i] = -2.0 * g_z_y_zz_xx[i] * a_exp + 4.0 * g_z_xxy_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_zz_xy[i] = -2.0 * g_z_y_zz_xy[i] * a_exp + 4.0 * g_z_xxy_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_zz_xz[i] = -2.0 * g_z_y_zz_xz[i] * a_exp + 4.0 * g_z_xxy_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_zz_yy[i] = -2.0 * g_z_y_zz_yy[i] * a_exp + 4.0 * g_z_xxy_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_zz_yz[i] = -2.0 * g_z_y_zz_yz[i] * a_exp + 4.0 * g_z_xxy_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_zz_zz[i] = -2.0 * g_z_y_zz_zz[i] * a_exp + 4.0 * g_z_xxy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1368-1374)

    #pragma omp simd aligned(g_z_x_0_0_0_xz_xx_xx, g_z_x_0_0_0_xz_xx_xy, g_z_x_0_0_0_xz_xx_xz, g_z_x_0_0_0_xz_xx_yy, g_z_x_0_0_0_xz_xx_yz, g_z_x_0_0_0_xz_xx_zz, g_z_xxz_xx_xx, g_z_xxz_xx_xy, g_z_xxz_xx_xz, g_z_xxz_xx_yy, g_z_xxz_xx_yz, g_z_xxz_xx_zz, g_z_z_xx_xx, g_z_z_xx_xy, g_z_z_xx_xz, g_z_z_xx_yy, g_z_z_xx_yz, g_z_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xz_xx_xx[i] = -2.0 * g_z_z_xx_xx[i] * a_exp + 4.0 * g_z_xxz_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_xx_xy[i] = -2.0 * g_z_z_xx_xy[i] * a_exp + 4.0 * g_z_xxz_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_xx_xz[i] = -2.0 * g_z_z_xx_xz[i] * a_exp + 4.0 * g_z_xxz_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_xx_yy[i] = -2.0 * g_z_z_xx_yy[i] * a_exp + 4.0 * g_z_xxz_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_xx_yz[i] = -2.0 * g_z_z_xx_yz[i] * a_exp + 4.0 * g_z_xxz_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_xx_zz[i] = -2.0 * g_z_z_xx_zz[i] * a_exp + 4.0 * g_z_xxz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1374-1380)

    #pragma omp simd aligned(g_z_x_0_0_0_xz_xy_xx, g_z_x_0_0_0_xz_xy_xy, g_z_x_0_0_0_xz_xy_xz, g_z_x_0_0_0_xz_xy_yy, g_z_x_0_0_0_xz_xy_yz, g_z_x_0_0_0_xz_xy_zz, g_z_xxz_xy_xx, g_z_xxz_xy_xy, g_z_xxz_xy_xz, g_z_xxz_xy_yy, g_z_xxz_xy_yz, g_z_xxz_xy_zz, g_z_z_xy_xx, g_z_z_xy_xy, g_z_z_xy_xz, g_z_z_xy_yy, g_z_z_xy_yz, g_z_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xz_xy_xx[i] = -2.0 * g_z_z_xy_xx[i] * a_exp + 4.0 * g_z_xxz_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_xy_xy[i] = -2.0 * g_z_z_xy_xy[i] * a_exp + 4.0 * g_z_xxz_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_xy_xz[i] = -2.0 * g_z_z_xy_xz[i] * a_exp + 4.0 * g_z_xxz_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_xy_yy[i] = -2.0 * g_z_z_xy_yy[i] * a_exp + 4.0 * g_z_xxz_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_xy_yz[i] = -2.0 * g_z_z_xy_yz[i] * a_exp + 4.0 * g_z_xxz_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_xy_zz[i] = -2.0 * g_z_z_xy_zz[i] * a_exp + 4.0 * g_z_xxz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1380-1386)

    #pragma omp simd aligned(g_z_x_0_0_0_xz_xz_xx, g_z_x_0_0_0_xz_xz_xy, g_z_x_0_0_0_xz_xz_xz, g_z_x_0_0_0_xz_xz_yy, g_z_x_0_0_0_xz_xz_yz, g_z_x_0_0_0_xz_xz_zz, g_z_xxz_xz_xx, g_z_xxz_xz_xy, g_z_xxz_xz_xz, g_z_xxz_xz_yy, g_z_xxz_xz_yz, g_z_xxz_xz_zz, g_z_z_xz_xx, g_z_z_xz_xy, g_z_z_xz_xz, g_z_z_xz_yy, g_z_z_xz_yz, g_z_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xz_xz_xx[i] = -2.0 * g_z_z_xz_xx[i] * a_exp + 4.0 * g_z_xxz_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_xz_xy[i] = -2.0 * g_z_z_xz_xy[i] * a_exp + 4.0 * g_z_xxz_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_xz_xz[i] = -2.0 * g_z_z_xz_xz[i] * a_exp + 4.0 * g_z_xxz_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_xz_yy[i] = -2.0 * g_z_z_xz_yy[i] * a_exp + 4.0 * g_z_xxz_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_xz_yz[i] = -2.0 * g_z_z_xz_yz[i] * a_exp + 4.0 * g_z_xxz_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_xz_zz[i] = -2.0 * g_z_z_xz_zz[i] * a_exp + 4.0 * g_z_xxz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1386-1392)

    #pragma omp simd aligned(g_z_x_0_0_0_xz_yy_xx, g_z_x_0_0_0_xz_yy_xy, g_z_x_0_0_0_xz_yy_xz, g_z_x_0_0_0_xz_yy_yy, g_z_x_0_0_0_xz_yy_yz, g_z_x_0_0_0_xz_yy_zz, g_z_xxz_yy_xx, g_z_xxz_yy_xy, g_z_xxz_yy_xz, g_z_xxz_yy_yy, g_z_xxz_yy_yz, g_z_xxz_yy_zz, g_z_z_yy_xx, g_z_z_yy_xy, g_z_z_yy_xz, g_z_z_yy_yy, g_z_z_yy_yz, g_z_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xz_yy_xx[i] = -2.0 * g_z_z_yy_xx[i] * a_exp + 4.0 * g_z_xxz_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_yy_xy[i] = -2.0 * g_z_z_yy_xy[i] * a_exp + 4.0 * g_z_xxz_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_yy_xz[i] = -2.0 * g_z_z_yy_xz[i] * a_exp + 4.0 * g_z_xxz_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_yy_yy[i] = -2.0 * g_z_z_yy_yy[i] * a_exp + 4.0 * g_z_xxz_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_yy_yz[i] = -2.0 * g_z_z_yy_yz[i] * a_exp + 4.0 * g_z_xxz_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_yy_zz[i] = -2.0 * g_z_z_yy_zz[i] * a_exp + 4.0 * g_z_xxz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1392-1398)

    #pragma omp simd aligned(g_z_x_0_0_0_xz_yz_xx, g_z_x_0_0_0_xz_yz_xy, g_z_x_0_0_0_xz_yz_xz, g_z_x_0_0_0_xz_yz_yy, g_z_x_0_0_0_xz_yz_yz, g_z_x_0_0_0_xz_yz_zz, g_z_xxz_yz_xx, g_z_xxz_yz_xy, g_z_xxz_yz_xz, g_z_xxz_yz_yy, g_z_xxz_yz_yz, g_z_xxz_yz_zz, g_z_z_yz_xx, g_z_z_yz_xy, g_z_z_yz_xz, g_z_z_yz_yy, g_z_z_yz_yz, g_z_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xz_yz_xx[i] = -2.0 * g_z_z_yz_xx[i] * a_exp + 4.0 * g_z_xxz_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_yz_xy[i] = -2.0 * g_z_z_yz_xy[i] * a_exp + 4.0 * g_z_xxz_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_yz_xz[i] = -2.0 * g_z_z_yz_xz[i] * a_exp + 4.0 * g_z_xxz_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_yz_yy[i] = -2.0 * g_z_z_yz_yy[i] * a_exp + 4.0 * g_z_xxz_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_yz_yz[i] = -2.0 * g_z_z_yz_yz[i] * a_exp + 4.0 * g_z_xxz_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_yz_zz[i] = -2.0 * g_z_z_yz_zz[i] * a_exp + 4.0 * g_z_xxz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1398-1404)

    #pragma omp simd aligned(g_z_x_0_0_0_xz_zz_xx, g_z_x_0_0_0_xz_zz_xy, g_z_x_0_0_0_xz_zz_xz, g_z_x_0_0_0_xz_zz_yy, g_z_x_0_0_0_xz_zz_yz, g_z_x_0_0_0_xz_zz_zz, g_z_xxz_zz_xx, g_z_xxz_zz_xy, g_z_xxz_zz_xz, g_z_xxz_zz_yy, g_z_xxz_zz_yz, g_z_xxz_zz_zz, g_z_z_zz_xx, g_z_z_zz_xy, g_z_z_zz_xz, g_z_z_zz_yy, g_z_z_zz_yz, g_z_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xz_zz_xx[i] = -2.0 * g_z_z_zz_xx[i] * a_exp + 4.0 * g_z_xxz_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_zz_xy[i] = -2.0 * g_z_z_zz_xy[i] * a_exp + 4.0 * g_z_xxz_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_zz_xz[i] = -2.0 * g_z_z_zz_xz[i] * a_exp + 4.0 * g_z_xxz_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_zz_yy[i] = -2.0 * g_z_z_zz_yy[i] * a_exp + 4.0 * g_z_xxz_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_zz_yz[i] = -2.0 * g_z_z_zz_yz[i] * a_exp + 4.0 * g_z_xxz_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_zz_zz[i] = -2.0 * g_z_z_zz_zz[i] * a_exp + 4.0 * g_z_xxz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1404-1410)

    #pragma omp simd aligned(g_z_x_0_0_0_yy_xx_xx, g_z_x_0_0_0_yy_xx_xy, g_z_x_0_0_0_yy_xx_xz, g_z_x_0_0_0_yy_xx_yy, g_z_x_0_0_0_yy_xx_yz, g_z_x_0_0_0_yy_xx_zz, g_z_xyy_xx_xx, g_z_xyy_xx_xy, g_z_xyy_xx_xz, g_z_xyy_xx_yy, g_z_xyy_xx_yz, g_z_xyy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yy_xx_xx[i] = 4.0 * g_z_xyy_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_xx_xy[i] = 4.0 * g_z_xyy_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_xx_xz[i] = 4.0 * g_z_xyy_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_xx_yy[i] = 4.0 * g_z_xyy_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_xx_yz[i] = 4.0 * g_z_xyy_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_xx_zz[i] = 4.0 * g_z_xyy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1410-1416)

    #pragma omp simd aligned(g_z_x_0_0_0_yy_xy_xx, g_z_x_0_0_0_yy_xy_xy, g_z_x_0_0_0_yy_xy_xz, g_z_x_0_0_0_yy_xy_yy, g_z_x_0_0_0_yy_xy_yz, g_z_x_0_0_0_yy_xy_zz, g_z_xyy_xy_xx, g_z_xyy_xy_xy, g_z_xyy_xy_xz, g_z_xyy_xy_yy, g_z_xyy_xy_yz, g_z_xyy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yy_xy_xx[i] = 4.0 * g_z_xyy_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_xy_xy[i] = 4.0 * g_z_xyy_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_xy_xz[i] = 4.0 * g_z_xyy_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_xy_yy[i] = 4.0 * g_z_xyy_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_xy_yz[i] = 4.0 * g_z_xyy_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_xy_zz[i] = 4.0 * g_z_xyy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1416-1422)

    #pragma omp simd aligned(g_z_x_0_0_0_yy_xz_xx, g_z_x_0_0_0_yy_xz_xy, g_z_x_0_0_0_yy_xz_xz, g_z_x_0_0_0_yy_xz_yy, g_z_x_0_0_0_yy_xz_yz, g_z_x_0_0_0_yy_xz_zz, g_z_xyy_xz_xx, g_z_xyy_xz_xy, g_z_xyy_xz_xz, g_z_xyy_xz_yy, g_z_xyy_xz_yz, g_z_xyy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yy_xz_xx[i] = 4.0 * g_z_xyy_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_xz_xy[i] = 4.0 * g_z_xyy_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_xz_xz[i] = 4.0 * g_z_xyy_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_xz_yy[i] = 4.0 * g_z_xyy_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_xz_yz[i] = 4.0 * g_z_xyy_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_xz_zz[i] = 4.0 * g_z_xyy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1422-1428)

    #pragma omp simd aligned(g_z_x_0_0_0_yy_yy_xx, g_z_x_0_0_0_yy_yy_xy, g_z_x_0_0_0_yy_yy_xz, g_z_x_0_0_0_yy_yy_yy, g_z_x_0_0_0_yy_yy_yz, g_z_x_0_0_0_yy_yy_zz, g_z_xyy_yy_xx, g_z_xyy_yy_xy, g_z_xyy_yy_xz, g_z_xyy_yy_yy, g_z_xyy_yy_yz, g_z_xyy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yy_yy_xx[i] = 4.0 * g_z_xyy_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_yy_xy[i] = 4.0 * g_z_xyy_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_yy_xz[i] = 4.0 * g_z_xyy_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_yy_yy[i] = 4.0 * g_z_xyy_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_yy_yz[i] = 4.0 * g_z_xyy_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_yy_zz[i] = 4.0 * g_z_xyy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1428-1434)

    #pragma omp simd aligned(g_z_x_0_0_0_yy_yz_xx, g_z_x_0_0_0_yy_yz_xy, g_z_x_0_0_0_yy_yz_xz, g_z_x_0_0_0_yy_yz_yy, g_z_x_0_0_0_yy_yz_yz, g_z_x_0_0_0_yy_yz_zz, g_z_xyy_yz_xx, g_z_xyy_yz_xy, g_z_xyy_yz_xz, g_z_xyy_yz_yy, g_z_xyy_yz_yz, g_z_xyy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yy_yz_xx[i] = 4.0 * g_z_xyy_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_yz_xy[i] = 4.0 * g_z_xyy_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_yz_xz[i] = 4.0 * g_z_xyy_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_yz_yy[i] = 4.0 * g_z_xyy_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_yz_yz[i] = 4.0 * g_z_xyy_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_yz_zz[i] = 4.0 * g_z_xyy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1434-1440)

    #pragma omp simd aligned(g_z_x_0_0_0_yy_zz_xx, g_z_x_0_0_0_yy_zz_xy, g_z_x_0_0_0_yy_zz_xz, g_z_x_0_0_0_yy_zz_yy, g_z_x_0_0_0_yy_zz_yz, g_z_x_0_0_0_yy_zz_zz, g_z_xyy_zz_xx, g_z_xyy_zz_xy, g_z_xyy_zz_xz, g_z_xyy_zz_yy, g_z_xyy_zz_yz, g_z_xyy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yy_zz_xx[i] = 4.0 * g_z_xyy_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_zz_xy[i] = 4.0 * g_z_xyy_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_zz_xz[i] = 4.0 * g_z_xyy_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_zz_yy[i] = 4.0 * g_z_xyy_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_zz_yz[i] = 4.0 * g_z_xyy_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_zz_zz[i] = 4.0 * g_z_xyy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1440-1446)

    #pragma omp simd aligned(g_z_x_0_0_0_yz_xx_xx, g_z_x_0_0_0_yz_xx_xy, g_z_x_0_0_0_yz_xx_xz, g_z_x_0_0_0_yz_xx_yy, g_z_x_0_0_0_yz_xx_yz, g_z_x_0_0_0_yz_xx_zz, g_z_xyz_xx_xx, g_z_xyz_xx_xy, g_z_xyz_xx_xz, g_z_xyz_xx_yy, g_z_xyz_xx_yz, g_z_xyz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yz_xx_xx[i] = 4.0 * g_z_xyz_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_xx_xy[i] = 4.0 * g_z_xyz_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_xx_xz[i] = 4.0 * g_z_xyz_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_xx_yy[i] = 4.0 * g_z_xyz_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_xx_yz[i] = 4.0 * g_z_xyz_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_xx_zz[i] = 4.0 * g_z_xyz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1446-1452)

    #pragma omp simd aligned(g_z_x_0_0_0_yz_xy_xx, g_z_x_0_0_0_yz_xy_xy, g_z_x_0_0_0_yz_xy_xz, g_z_x_0_0_0_yz_xy_yy, g_z_x_0_0_0_yz_xy_yz, g_z_x_0_0_0_yz_xy_zz, g_z_xyz_xy_xx, g_z_xyz_xy_xy, g_z_xyz_xy_xz, g_z_xyz_xy_yy, g_z_xyz_xy_yz, g_z_xyz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yz_xy_xx[i] = 4.0 * g_z_xyz_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_xy_xy[i] = 4.0 * g_z_xyz_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_xy_xz[i] = 4.0 * g_z_xyz_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_xy_yy[i] = 4.0 * g_z_xyz_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_xy_yz[i] = 4.0 * g_z_xyz_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_xy_zz[i] = 4.0 * g_z_xyz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1452-1458)

    #pragma omp simd aligned(g_z_x_0_0_0_yz_xz_xx, g_z_x_0_0_0_yz_xz_xy, g_z_x_0_0_0_yz_xz_xz, g_z_x_0_0_0_yz_xz_yy, g_z_x_0_0_0_yz_xz_yz, g_z_x_0_0_0_yz_xz_zz, g_z_xyz_xz_xx, g_z_xyz_xz_xy, g_z_xyz_xz_xz, g_z_xyz_xz_yy, g_z_xyz_xz_yz, g_z_xyz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yz_xz_xx[i] = 4.0 * g_z_xyz_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_xz_xy[i] = 4.0 * g_z_xyz_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_xz_xz[i] = 4.0 * g_z_xyz_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_xz_yy[i] = 4.0 * g_z_xyz_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_xz_yz[i] = 4.0 * g_z_xyz_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_xz_zz[i] = 4.0 * g_z_xyz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1458-1464)

    #pragma omp simd aligned(g_z_x_0_0_0_yz_yy_xx, g_z_x_0_0_0_yz_yy_xy, g_z_x_0_0_0_yz_yy_xz, g_z_x_0_0_0_yz_yy_yy, g_z_x_0_0_0_yz_yy_yz, g_z_x_0_0_0_yz_yy_zz, g_z_xyz_yy_xx, g_z_xyz_yy_xy, g_z_xyz_yy_xz, g_z_xyz_yy_yy, g_z_xyz_yy_yz, g_z_xyz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yz_yy_xx[i] = 4.0 * g_z_xyz_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_yy_xy[i] = 4.0 * g_z_xyz_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_yy_xz[i] = 4.0 * g_z_xyz_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_yy_yy[i] = 4.0 * g_z_xyz_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_yy_yz[i] = 4.0 * g_z_xyz_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_yy_zz[i] = 4.0 * g_z_xyz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1464-1470)

    #pragma omp simd aligned(g_z_x_0_0_0_yz_yz_xx, g_z_x_0_0_0_yz_yz_xy, g_z_x_0_0_0_yz_yz_xz, g_z_x_0_0_0_yz_yz_yy, g_z_x_0_0_0_yz_yz_yz, g_z_x_0_0_0_yz_yz_zz, g_z_xyz_yz_xx, g_z_xyz_yz_xy, g_z_xyz_yz_xz, g_z_xyz_yz_yy, g_z_xyz_yz_yz, g_z_xyz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yz_yz_xx[i] = 4.0 * g_z_xyz_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_yz_xy[i] = 4.0 * g_z_xyz_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_yz_xz[i] = 4.0 * g_z_xyz_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_yz_yy[i] = 4.0 * g_z_xyz_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_yz_yz[i] = 4.0 * g_z_xyz_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_yz_zz[i] = 4.0 * g_z_xyz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1470-1476)

    #pragma omp simd aligned(g_z_x_0_0_0_yz_zz_xx, g_z_x_0_0_0_yz_zz_xy, g_z_x_0_0_0_yz_zz_xz, g_z_x_0_0_0_yz_zz_yy, g_z_x_0_0_0_yz_zz_yz, g_z_x_0_0_0_yz_zz_zz, g_z_xyz_zz_xx, g_z_xyz_zz_xy, g_z_xyz_zz_xz, g_z_xyz_zz_yy, g_z_xyz_zz_yz, g_z_xyz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yz_zz_xx[i] = 4.0 * g_z_xyz_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_zz_xy[i] = 4.0 * g_z_xyz_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_zz_xz[i] = 4.0 * g_z_xyz_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_zz_yy[i] = 4.0 * g_z_xyz_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_zz_yz[i] = 4.0 * g_z_xyz_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_zz_zz[i] = 4.0 * g_z_xyz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1476-1482)

    #pragma omp simd aligned(g_z_x_0_0_0_zz_xx_xx, g_z_x_0_0_0_zz_xx_xy, g_z_x_0_0_0_zz_xx_xz, g_z_x_0_0_0_zz_xx_yy, g_z_x_0_0_0_zz_xx_yz, g_z_x_0_0_0_zz_xx_zz, g_z_xzz_xx_xx, g_z_xzz_xx_xy, g_z_xzz_xx_xz, g_z_xzz_xx_yy, g_z_xzz_xx_yz, g_z_xzz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_zz_xx_xx[i] = 4.0 * g_z_xzz_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_xx_xy[i] = 4.0 * g_z_xzz_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_xx_xz[i] = 4.0 * g_z_xzz_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_xx_yy[i] = 4.0 * g_z_xzz_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_xx_yz[i] = 4.0 * g_z_xzz_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_xx_zz[i] = 4.0 * g_z_xzz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1482-1488)

    #pragma omp simd aligned(g_z_x_0_0_0_zz_xy_xx, g_z_x_0_0_0_zz_xy_xy, g_z_x_0_0_0_zz_xy_xz, g_z_x_0_0_0_zz_xy_yy, g_z_x_0_0_0_zz_xy_yz, g_z_x_0_0_0_zz_xy_zz, g_z_xzz_xy_xx, g_z_xzz_xy_xy, g_z_xzz_xy_xz, g_z_xzz_xy_yy, g_z_xzz_xy_yz, g_z_xzz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_zz_xy_xx[i] = 4.0 * g_z_xzz_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_xy_xy[i] = 4.0 * g_z_xzz_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_xy_xz[i] = 4.0 * g_z_xzz_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_xy_yy[i] = 4.0 * g_z_xzz_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_xy_yz[i] = 4.0 * g_z_xzz_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_xy_zz[i] = 4.0 * g_z_xzz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1488-1494)

    #pragma omp simd aligned(g_z_x_0_0_0_zz_xz_xx, g_z_x_0_0_0_zz_xz_xy, g_z_x_0_0_0_zz_xz_xz, g_z_x_0_0_0_zz_xz_yy, g_z_x_0_0_0_zz_xz_yz, g_z_x_0_0_0_zz_xz_zz, g_z_xzz_xz_xx, g_z_xzz_xz_xy, g_z_xzz_xz_xz, g_z_xzz_xz_yy, g_z_xzz_xz_yz, g_z_xzz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_zz_xz_xx[i] = 4.0 * g_z_xzz_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_xz_xy[i] = 4.0 * g_z_xzz_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_xz_xz[i] = 4.0 * g_z_xzz_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_xz_yy[i] = 4.0 * g_z_xzz_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_xz_yz[i] = 4.0 * g_z_xzz_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_xz_zz[i] = 4.0 * g_z_xzz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1494-1500)

    #pragma omp simd aligned(g_z_x_0_0_0_zz_yy_xx, g_z_x_0_0_0_zz_yy_xy, g_z_x_0_0_0_zz_yy_xz, g_z_x_0_0_0_zz_yy_yy, g_z_x_0_0_0_zz_yy_yz, g_z_x_0_0_0_zz_yy_zz, g_z_xzz_yy_xx, g_z_xzz_yy_xy, g_z_xzz_yy_xz, g_z_xzz_yy_yy, g_z_xzz_yy_yz, g_z_xzz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_zz_yy_xx[i] = 4.0 * g_z_xzz_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_yy_xy[i] = 4.0 * g_z_xzz_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_yy_xz[i] = 4.0 * g_z_xzz_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_yy_yy[i] = 4.0 * g_z_xzz_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_yy_yz[i] = 4.0 * g_z_xzz_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_yy_zz[i] = 4.0 * g_z_xzz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1500-1506)

    #pragma omp simd aligned(g_z_x_0_0_0_zz_yz_xx, g_z_x_0_0_0_zz_yz_xy, g_z_x_0_0_0_zz_yz_xz, g_z_x_0_0_0_zz_yz_yy, g_z_x_0_0_0_zz_yz_yz, g_z_x_0_0_0_zz_yz_zz, g_z_xzz_yz_xx, g_z_xzz_yz_xy, g_z_xzz_yz_xz, g_z_xzz_yz_yy, g_z_xzz_yz_yz, g_z_xzz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_zz_yz_xx[i] = 4.0 * g_z_xzz_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_yz_xy[i] = 4.0 * g_z_xzz_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_yz_xz[i] = 4.0 * g_z_xzz_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_yz_yy[i] = 4.0 * g_z_xzz_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_yz_yz[i] = 4.0 * g_z_xzz_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_yz_zz[i] = 4.0 * g_z_xzz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1506-1512)

    #pragma omp simd aligned(g_z_x_0_0_0_zz_zz_xx, g_z_x_0_0_0_zz_zz_xy, g_z_x_0_0_0_zz_zz_xz, g_z_x_0_0_0_zz_zz_yy, g_z_x_0_0_0_zz_zz_yz, g_z_x_0_0_0_zz_zz_zz, g_z_xzz_zz_xx, g_z_xzz_zz_xy, g_z_xzz_zz_xz, g_z_xzz_zz_yy, g_z_xzz_zz_yz, g_z_xzz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_zz_zz_xx[i] = 4.0 * g_z_xzz_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_zz_xy[i] = 4.0 * g_z_xzz_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_zz_xz[i] = 4.0 * g_z_xzz_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_zz_yy[i] = 4.0 * g_z_xzz_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_zz_yz[i] = 4.0 * g_z_xzz_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_zz_zz[i] = 4.0 * g_z_xzz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1512-1518)

    #pragma omp simd aligned(g_z_xxy_xx_xx, g_z_xxy_xx_xy, g_z_xxy_xx_xz, g_z_xxy_xx_yy, g_z_xxy_xx_yz, g_z_xxy_xx_zz, g_z_y_0_0_0_xx_xx_xx, g_z_y_0_0_0_xx_xx_xy, g_z_y_0_0_0_xx_xx_xz, g_z_y_0_0_0_xx_xx_yy, g_z_y_0_0_0_xx_xx_yz, g_z_y_0_0_0_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xx_xx_xx[i] = 4.0 * g_z_xxy_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_xx_xy[i] = 4.0 * g_z_xxy_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_xx_xz[i] = 4.0 * g_z_xxy_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_xx_yy[i] = 4.0 * g_z_xxy_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_xx_yz[i] = 4.0 * g_z_xxy_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_xx_zz[i] = 4.0 * g_z_xxy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1518-1524)

    #pragma omp simd aligned(g_z_xxy_xy_xx, g_z_xxy_xy_xy, g_z_xxy_xy_xz, g_z_xxy_xy_yy, g_z_xxy_xy_yz, g_z_xxy_xy_zz, g_z_y_0_0_0_xx_xy_xx, g_z_y_0_0_0_xx_xy_xy, g_z_y_0_0_0_xx_xy_xz, g_z_y_0_0_0_xx_xy_yy, g_z_y_0_0_0_xx_xy_yz, g_z_y_0_0_0_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xx_xy_xx[i] = 4.0 * g_z_xxy_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_xy_xy[i] = 4.0 * g_z_xxy_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_xy_xz[i] = 4.0 * g_z_xxy_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_xy_yy[i] = 4.0 * g_z_xxy_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_xy_yz[i] = 4.0 * g_z_xxy_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_xy_zz[i] = 4.0 * g_z_xxy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1524-1530)

    #pragma omp simd aligned(g_z_xxy_xz_xx, g_z_xxy_xz_xy, g_z_xxy_xz_xz, g_z_xxy_xz_yy, g_z_xxy_xz_yz, g_z_xxy_xz_zz, g_z_y_0_0_0_xx_xz_xx, g_z_y_0_0_0_xx_xz_xy, g_z_y_0_0_0_xx_xz_xz, g_z_y_0_0_0_xx_xz_yy, g_z_y_0_0_0_xx_xz_yz, g_z_y_0_0_0_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xx_xz_xx[i] = 4.0 * g_z_xxy_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_xz_xy[i] = 4.0 * g_z_xxy_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_xz_xz[i] = 4.0 * g_z_xxy_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_xz_yy[i] = 4.0 * g_z_xxy_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_xz_yz[i] = 4.0 * g_z_xxy_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_xz_zz[i] = 4.0 * g_z_xxy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1530-1536)

    #pragma omp simd aligned(g_z_xxy_yy_xx, g_z_xxy_yy_xy, g_z_xxy_yy_xz, g_z_xxy_yy_yy, g_z_xxy_yy_yz, g_z_xxy_yy_zz, g_z_y_0_0_0_xx_yy_xx, g_z_y_0_0_0_xx_yy_xy, g_z_y_0_0_0_xx_yy_xz, g_z_y_0_0_0_xx_yy_yy, g_z_y_0_0_0_xx_yy_yz, g_z_y_0_0_0_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xx_yy_xx[i] = 4.0 * g_z_xxy_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_yy_xy[i] = 4.0 * g_z_xxy_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_yy_xz[i] = 4.0 * g_z_xxy_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_yy_yy[i] = 4.0 * g_z_xxy_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_yy_yz[i] = 4.0 * g_z_xxy_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_yy_zz[i] = 4.0 * g_z_xxy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1536-1542)

    #pragma omp simd aligned(g_z_xxy_yz_xx, g_z_xxy_yz_xy, g_z_xxy_yz_xz, g_z_xxy_yz_yy, g_z_xxy_yz_yz, g_z_xxy_yz_zz, g_z_y_0_0_0_xx_yz_xx, g_z_y_0_0_0_xx_yz_xy, g_z_y_0_0_0_xx_yz_xz, g_z_y_0_0_0_xx_yz_yy, g_z_y_0_0_0_xx_yz_yz, g_z_y_0_0_0_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xx_yz_xx[i] = 4.0 * g_z_xxy_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_yz_xy[i] = 4.0 * g_z_xxy_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_yz_xz[i] = 4.0 * g_z_xxy_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_yz_yy[i] = 4.0 * g_z_xxy_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_yz_yz[i] = 4.0 * g_z_xxy_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_yz_zz[i] = 4.0 * g_z_xxy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1542-1548)

    #pragma omp simd aligned(g_z_xxy_zz_xx, g_z_xxy_zz_xy, g_z_xxy_zz_xz, g_z_xxy_zz_yy, g_z_xxy_zz_yz, g_z_xxy_zz_zz, g_z_y_0_0_0_xx_zz_xx, g_z_y_0_0_0_xx_zz_xy, g_z_y_0_0_0_xx_zz_xz, g_z_y_0_0_0_xx_zz_yy, g_z_y_0_0_0_xx_zz_yz, g_z_y_0_0_0_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xx_zz_xx[i] = 4.0 * g_z_xxy_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_zz_xy[i] = 4.0 * g_z_xxy_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_zz_xz[i] = 4.0 * g_z_xxy_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_zz_yy[i] = 4.0 * g_z_xxy_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_zz_yz[i] = 4.0 * g_z_xxy_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_zz_zz[i] = 4.0 * g_z_xxy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1548-1554)

    #pragma omp simd aligned(g_z_x_xx_xx, g_z_x_xx_xy, g_z_x_xx_xz, g_z_x_xx_yy, g_z_x_xx_yz, g_z_x_xx_zz, g_z_xyy_xx_xx, g_z_xyy_xx_xy, g_z_xyy_xx_xz, g_z_xyy_xx_yy, g_z_xyy_xx_yz, g_z_xyy_xx_zz, g_z_y_0_0_0_xy_xx_xx, g_z_y_0_0_0_xy_xx_xy, g_z_y_0_0_0_xy_xx_xz, g_z_y_0_0_0_xy_xx_yy, g_z_y_0_0_0_xy_xx_yz, g_z_y_0_0_0_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xy_xx_xx[i] = -2.0 * g_z_x_xx_xx[i] * a_exp + 4.0 * g_z_xyy_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_xx_xy[i] = -2.0 * g_z_x_xx_xy[i] * a_exp + 4.0 * g_z_xyy_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_xx_xz[i] = -2.0 * g_z_x_xx_xz[i] * a_exp + 4.0 * g_z_xyy_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_xx_yy[i] = -2.0 * g_z_x_xx_yy[i] * a_exp + 4.0 * g_z_xyy_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_xx_yz[i] = -2.0 * g_z_x_xx_yz[i] * a_exp + 4.0 * g_z_xyy_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_xx_zz[i] = -2.0 * g_z_x_xx_zz[i] * a_exp + 4.0 * g_z_xyy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1554-1560)

    #pragma omp simd aligned(g_z_x_xy_xx, g_z_x_xy_xy, g_z_x_xy_xz, g_z_x_xy_yy, g_z_x_xy_yz, g_z_x_xy_zz, g_z_xyy_xy_xx, g_z_xyy_xy_xy, g_z_xyy_xy_xz, g_z_xyy_xy_yy, g_z_xyy_xy_yz, g_z_xyy_xy_zz, g_z_y_0_0_0_xy_xy_xx, g_z_y_0_0_0_xy_xy_xy, g_z_y_0_0_0_xy_xy_xz, g_z_y_0_0_0_xy_xy_yy, g_z_y_0_0_0_xy_xy_yz, g_z_y_0_0_0_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xy_xy_xx[i] = -2.0 * g_z_x_xy_xx[i] * a_exp + 4.0 * g_z_xyy_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_xy_xy[i] = -2.0 * g_z_x_xy_xy[i] * a_exp + 4.0 * g_z_xyy_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_xy_xz[i] = -2.0 * g_z_x_xy_xz[i] * a_exp + 4.0 * g_z_xyy_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_xy_yy[i] = -2.0 * g_z_x_xy_yy[i] * a_exp + 4.0 * g_z_xyy_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_xy_yz[i] = -2.0 * g_z_x_xy_yz[i] * a_exp + 4.0 * g_z_xyy_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_xy_zz[i] = -2.0 * g_z_x_xy_zz[i] * a_exp + 4.0 * g_z_xyy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1560-1566)

    #pragma omp simd aligned(g_z_x_xz_xx, g_z_x_xz_xy, g_z_x_xz_xz, g_z_x_xz_yy, g_z_x_xz_yz, g_z_x_xz_zz, g_z_xyy_xz_xx, g_z_xyy_xz_xy, g_z_xyy_xz_xz, g_z_xyy_xz_yy, g_z_xyy_xz_yz, g_z_xyy_xz_zz, g_z_y_0_0_0_xy_xz_xx, g_z_y_0_0_0_xy_xz_xy, g_z_y_0_0_0_xy_xz_xz, g_z_y_0_0_0_xy_xz_yy, g_z_y_0_0_0_xy_xz_yz, g_z_y_0_0_0_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xy_xz_xx[i] = -2.0 * g_z_x_xz_xx[i] * a_exp + 4.0 * g_z_xyy_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_xz_xy[i] = -2.0 * g_z_x_xz_xy[i] * a_exp + 4.0 * g_z_xyy_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_xz_xz[i] = -2.0 * g_z_x_xz_xz[i] * a_exp + 4.0 * g_z_xyy_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_xz_yy[i] = -2.0 * g_z_x_xz_yy[i] * a_exp + 4.0 * g_z_xyy_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_xz_yz[i] = -2.0 * g_z_x_xz_yz[i] * a_exp + 4.0 * g_z_xyy_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_xz_zz[i] = -2.0 * g_z_x_xz_zz[i] * a_exp + 4.0 * g_z_xyy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1566-1572)

    #pragma omp simd aligned(g_z_x_yy_xx, g_z_x_yy_xy, g_z_x_yy_xz, g_z_x_yy_yy, g_z_x_yy_yz, g_z_x_yy_zz, g_z_xyy_yy_xx, g_z_xyy_yy_xy, g_z_xyy_yy_xz, g_z_xyy_yy_yy, g_z_xyy_yy_yz, g_z_xyy_yy_zz, g_z_y_0_0_0_xy_yy_xx, g_z_y_0_0_0_xy_yy_xy, g_z_y_0_0_0_xy_yy_xz, g_z_y_0_0_0_xy_yy_yy, g_z_y_0_0_0_xy_yy_yz, g_z_y_0_0_0_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xy_yy_xx[i] = -2.0 * g_z_x_yy_xx[i] * a_exp + 4.0 * g_z_xyy_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_yy_xy[i] = -2.0 * g_z_x_yy_xy[i] * a_exp + 4.0 * g_z_xyy_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_yy_xz[i] = -2.0 * g_z_x_yy_xz[i] * a_exp + 4.0 * g_z_xyy_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_yy_yy[i] = -2.0 * g_z_x_yy_yy[i] * a_exp + 4.0 * g_z_xyy_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_yy_yz[i] = -2.0 * g_z_x_yy_yz[i] * a_exp + 4.0 * g_z_xyy_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_yy_zz[i] = -2.0 * g_z_x_yy_zz[i] * a_exp + 4.0 * g_z_xyy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1572-1578)

    #pragma omp simd aligned(g_z_x_yz_xx, g_z_x_yz_xy, g_z_x_yz_xz, g_z_x_yz_yy, g_z_x_yz_yz, g_z_x_yz_zz, g_z_xyy_yz_xx, g_z_xyy_yz_xy, g_z_xyy_yz_xz, g_z_xyy_yz_yy, g_z_xyy_yz_yz, g_z_xyy_yz_zz, g_z_y_0_0_0_xy_yz_xx, g_z_y_0_0_0_xy_yz_xy, g_z_y_0_0_0_xy_yz_xz, g_z_y_0_0_0_xy_yz_yy, g_z_y_0_0_0_xy_yz_yz, g_z_y_0_0_0_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xy_yz_xx[i] = -2.0 * g_z_x_yz_xx[i] * a_exp + 4.0 * g_z_xyy_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_yz_xy[i] = -2.0 * g_z_x_yz_xy[i] * a_exp + 4.0 * g_z_xyy_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_yz_xz[i] = -2.0 * g_z_x_yz_xz[i] * a_exp + 4.0 * g_z_xyy_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_yz_yy[i] = -2.0 * g_z_x_yz_yy[i] * a_exp + 4.0 * g_z_xyy_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_yz_yz[i] = -2.0 * g_z_x_yz_yz[i] * a_exp + 4.0 * g_z_xyy_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_yz_zz[i] = -2.0 * g_z_x_yz_zz[i] * a_exp + 4.0 * g_z_xyy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1578-1584)

    #pragma omp simd aligned(g_z_x_zz_xx, g_z_x_zz_xy, g_z_x_zz_xz, g_z_x_zz_yy, g_z_x_zz_yz, g_z_x_zz_zz, g_z_xyy_zz_xx, g_z_xyy_zz_xy, g_z_xyy_zz_xz, g_z_xyy_zz_yy, g_z_xyy_zz_yz, g_z_xyy_zz_zz, g_z_y_0_0_0_xy_zz_xx, g_z_y_0_0_0_xy_zz_xy, g_z_y_0_0_0_xy_zz_xz, g_z_y_0_0_0_xy_zz_yy, g_z_y_0_0_0_xy_zz_yz, g_z_y_0_0_0_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xy_zz_xx[i] = -2.0 * g_z_x_zz_xx[i] * a_exp + 4.0 * g_z_xyy_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_zz_xy[i] = -2.0 * g_z_x_zz_xy[i] * a_exp + 4.0 * g_z_xyy_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_zz_xz[i] = -2.0 * g_z_x_zz_xz[i] * a_exp + 4.0 * g_z_xyy_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_zz_yy[i] = -2.0 * g_z_x_zz_yy[i] * a_exp + 4.0 * g_z_xyy_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_zz_yz[i] = -2.0 * g_z_x_zz_yz[i] * a_exp + 4.0 * g_z_xyy_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_zz_zz[i] = -2.0 * g_z_x_zz_zz[i] * a_exp + 4.0 * g_z_xyy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1584-1590)

    #pragma omp simd aligned(g_z_xyz_xx_xx, g_z_xyz_xx_xy, g_z_xyz_xx_xz, g_z_xyz_xx_yy, g_z_xyz_xx_yz, g_z_xyz_xx_zz, g_z_y_0_0_0_xz_xx_xx, g_z_y_0_0_0_xz_xx_xy, g_z_y_0_0_0_xz_xx_xz, g_z_y_0_0_0_xz_xx_yy, g_z_y_0_0_0_xz_xx_yz, g_z_y_0_0_0_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xz_xx_xx[i] = 4.0 * g_z_xyz_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_xx_xy[i] = 4.0 * g_z_xyz_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_xx_xz[i] = 4.0 * g_z_xyz_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_xx_yy[i] = 4.0 * g_z_xyz_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_xx_yz[i] = 4.0 * g_z_xyz_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_xx_zz[i] = 4.0 * g_z_xyz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1590-1596)

    #pragma omp simd aligned(g_z_xyz_xy_xx, g_z_xyz_xy_xy, g_z_xyz_xy_xz, g_z_xyz_xy_yy, g_z_xyz_xy_yz, g_z_xyz_xy_zz, g_z_y_0_0_0_xz_xy_xx, g_z_y_0_0_0_xz_xy_xy, g_z_y_0_0_0_xz_xy_xz, g_z_y_0_0_0_xz_xy_yy, g_z_y_0_0_0_xz_xy_yz, g_z_y_0_0_0_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xz_xy_xx[i] = 4.0 * g_z_xyz_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_xy_xy[i] = 4.0 * g_z_xyz_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_xy_xz[i] = 4.0 * g_z_xyz_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_xy_yy[i] = 4.0 * g_z_xyz_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_xy_yz[i] = 4.0 * g_z_xyz_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_xy_zz[i] = 4.0 * g_z_xyz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1596-1602)

    #pragma omp simd aligned(g_z_xyz_xz_xx, g_z_xyz_xz_xy, g_z_xyz_xz_xz, g_z_xyz_xz_yy, g_z_xyz_xz_yz, g_z_xyz_xz_zz, g_z_y_0_0_0_xz_xz_xx, g_z_y_0_0_0_xz_xz_xy, g_z_y_0_0_0_xz_xz_xz, g_z_y_0_0_0_xz_xz_yy, g_z_y_0_0_0_xz_xz_yz, g_z_y_0_0_0_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xz_xz_xx[i] = 4.0 * g_z_xyz_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_xz_xy[i] = 4.0 * g_z_xyz_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_xz_xz[i] = 4.0 * g_z_xyz_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_xz_yy[i] = 4.0 * g_z_xyz_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_xz_yz[i] = 4.0 * g_z_xyz_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_xz_zz[i] = 4.0 * g_z_xyz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1602-1608)

    #pragma omp simd aligned(g_z_xyz_yy_xx, g_z_xyz_yy_xy, g_z_xyz_yy_xz, g_z_xyz_yy_yy, g_z_xyz_yy_yz, g_z_xyz_yy_zz, g_z_y_0_0_0_xz_yy_xx, g_z_y_0_0_0_xz_yy_xy, g_z_y_0_0_0_xz_yy_xz, g_z_y_0_0_0_xz_yy_yy, g_z_y_0_0_0_xz_yy_yz, g_z_y_0_0_0_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xz_yy_xx[i] = 4.0 * g_z_xyz_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_yy_xy[i] = 4.0 * g_z_xyz_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_yy_xz[i] = 4.0 * g_z_xyz_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_yy_yy[i] = 4.0 * g_z_xyz_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_yy_yz[i] = 4.0 * g_z_xyz_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_yy_zz[i] = 4.0 * g_z_xyz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1608-1614)

    #pragma omp simd aligned(g_z_xyz_yz_xx, g_z_xyz_yz_xy, g_z_xyz_yz_xz, g_z_xyz_yz_yy, g_z_xyz_yz_yz, g_z_xyz_yz_zz, g_z_y_0_0_0_xz_yz_xx, g_z_y_0_0_0_xz_yz_xy, g_z_y_0_0_0_xz_yz_xz, g_z_y_0_0_0_xz_yz_yy, g_z_y_0_0_0_xz_yz_yz, g_z_y_0_0_0_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xz_yz_xx[i] = 4.0 * g_z_xyz_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_yz_xy[i] = 4.0 * g_z_xyz_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_yz_xz[i] = 4.0 * g_z_xyz_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_yz_yy[i] = 4.0 * g_z_xyz_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_yz_yz[i] = 4.0 * g_z_xyz_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_yz_zz[i] = 4.0 * g_z_xyz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1614-1620)

    #pragma omp simd aligned(g_z_xyz_zz_xx, g_z_xyz_zz_xy, g_z_xyz_zz_xz, g_z_xyz_zz_yy, g_z_xyz_zz_yz, g_z_xyz_zz_zz, g_z_y_0_0_0_xz_zz_xx, g_z_y_0_0_0_xz_zz_xy, g_z_y_0_0_0_xz_zz_xz, g_z_y_0_0_0_xz_zz_yy, g_z_y_0_0_0_xz_zz_yz, g_z_y_0_0_0_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xz_zz_xx[i] = 4.0 * g_z_xyz_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_zz_xy[i] = 4.0 * g_z_xyz_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_zz_xz[i] = 4.0 * g_z_xyz_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_zz_yy[i] = 4.0 * g_z_xyz_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_zz_yz[i] = 4.0 * g_z_xyz_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_zz_zz[i] = 4.0 * g_z_xyz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1620-1626)

    #pragma omp simd aligned(g_z_y_0_0_0_yy_xx_xx, g_z_y_0_0_0_yy_xx_xy, g_z_y_0_0_0_yy_xx_xz, g_z_y_0_0_0_yy_xx_yy, g_z_y_0_0_0_yy_xx_yz, g_z_y_0_0_0_yy_xx_zz, g_z_y_xx_xx, g_z_y_xx_xy, g_z_y_xx_xz, g_z_y_xx_yy, g_z_y_xx_yz, g_z_y_xx_zz, g_z_yyy_xx_xx, g_z_yyy_xx_xy, g_z_yyy_xx_xz, g_z_yyy_xx_yy, g_z_yyy_xx_yz, g_z_yyy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yy_xx_xx[i] = -4.0 * g_z_y_xx_xx[i] * a_exp + 4.0 * g_z_yyy_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_xx_xy[i] = -4.0 * g_z_y_xx_xy[i] * a_exp + 4.0 * g_z_yyy_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_xx_xz[i] = -4.0 * g_z_y_xx_xz[i] * a_exp + 4.0 * g_z_yyy_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_xx_yy[i] = -4.0 * g_z_y_xx_yy[i] * a_exp + 4.0 * g_z_yyy_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_xx_yz[i] = -4.0 * g_z_y_xx_yz[i] * a_exp + 4.0 * g_z_yyy_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_xx_zz[i] = -4.0 * g_z_y_xx_zz[i] * a_exp + 4.0 * g_z_yyy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1626-1632)

    #pragma omp simd aligned(g_z_y_0_0_0_yy_xy_xx, g_z_y_0_0_0_yy_xy_xy, g_z_y_0_0_0_yy_xy_xz, g_z_y_0_0_0_yy_xy_yy, g_z_y_0_0_0_yy_xy_yz, g_z_y_0_0_0_yy_xy_zz, g_z_y_xy_xx, g_z_y_xy_xy, g_z_y_xy_xz, g_z_y_xy_yy, g_z_y_xy_yz, g_z_y_xy_zz, g_z_yyy_xy_xx, g_z_yyy_xy_xy, g_z_yyy_xy_xz, g_z_yyy_xy_yy, g_z_yyy_xy_yz, g_z_yyy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yy_xy_xx[i] = -4.0 * g_z_y_xy_xx[i] * a_exp + 4.0 * g_z_yyy_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_xy_xy[i] = -4.0 * g_z_y_xy_xy[i] * a_exp + 4.0 * g_z_yyy_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_xy_xz[i] = -4.0 * g_z_y_xy_xz[i] * a_exp + 4.0 * g_z_yyy_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_xy_yy[i] = -4.0 * g_z_y_xy_yy[i] * a_exp + 4.0 * g_z_yyy_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_xy_yz[i] = -4.0 * g_z_y_xy_yz[i] * a_exp + 4.0 * g_z_yyy_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_xy_zz[i] = -4.0 * g_z_y_xy_zz[i] * a_exp + 4.0 * g_z_yyy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1632-1638)

    #pragma omp simd aligned(g_z_y_0_0_0_yy_xz_xx, g_z_y_0_0_0_yy_xz_xy, g_z_y_0_0_0_yy_xz_xz, g_z_y_0_0_0_yy_xz_yy, g_z_y_0_0_0_yy_xz_yz, g_z_y_0_0_0_yy_xz_zz, g_z_y_xz_xx, g_z_y_xz_xy, g_z_y_xz_xz, g_z_y_xz_yy, g_z_y_xz_yz, g_z_y_xz_zz, g_z_yyy_xz_xx, g_z_yyy_xz_xy, g_z_yyy_xz_xz, g_z_yyy_xz_yy, g_z_yyy_xz_yz, g_z_yyy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yy_xz_xx[i] = -4.0 * g_z_y_xz_xx[i] * a_exp + 4.0 * g_z_yyy_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_xz_xy[i] = -4.0 * g_z_y_xz_xy[i] * a_exp + 4.0 * g_z_yyy_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_xz_xz[i] = -4.0 * g_z_y_xz_xz[i] * a_exp + 4.0 * g_z_yyy_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_xz_yy[i] = -4.0 * g_z_y_xz_yy[i] * a_exp + 4.0 * g_z_yyy_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_xz_yz[i] = -4.0 * g_z_y_xz_yz[i] * a_exp + 4.0 * g_z_yyy_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_xz_zz[i] = -4.0 * g_z_y_xz_zz[i] * a_exp + 4.0 * g_z_yyy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1638-1644)

    #pragma omp simd aligned(g_z_y_0_0_0_yy_yy_xx, g_z_y_0_0_0_yy_yy_xy, g_z_y_0_0_0_yy_yy_xz, g_z_y_0_0_0_yy_yy_yy, g_z_y_0_0_0_yy_yy_yz, g_z_y_0_0_0_yy_yy_zz, g_z_y_yy_xx, g_z_y_yy_xy, g_z_y_yy_xz, g_z_y_yy_yy, g_z_y_yy_yz, g_z_y_yy_zz, g_z_yyy_yy_xx, g_z_yyy_yy_xy, g_z_yyy_yy_xz, g_z_yyy_yy_yy, g_z_yyy_yy_yz, g_z_yyy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yy_yy_xx[i] = -4.0 * g_z_y_yy_xx[i] * a_exp + 4.0 * g_z_yyy_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_yy_xy[i] = -4.0 * g_z_y_yy_xy[i] * a_exp + 4.0 * g_z_yyy_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_yy_xz[i] = -4.0 * g_z_y_yy_xz[i] * a_exp + 4.0 * g_z_yyy_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_yy_yy[i] = -4.0 * g_z_y_yy_yy[i] * a_exp + 4.0 * g_z_yyy_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_yy_yz[i] = -4.0 * g_z_y_yy_yz[i] * a_exp + 4.0 * g_z_yyy_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_yy_zz[i] = -4.0 * g_z_y_yy_zz[i] * a_exp + 4.0 * g_z_yyy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1644-1650)

    #pragma omp simd aligned(g_z_y_0_0_0_yy_yz_xx, g_z_y_0_0_0_yy_yz_xy, g_z_y_0_0_0_yy_yz_xz, g_z_y_0_0_0_yy_yz_yy, g_z_y_0_0_0_yy_yz_yz, g_z_y_0_0_0_yy_yz_zz, g_z_y_yz_xx, g_z_y_yz_xy, g_z_y_yz_xz, g_z_y_yz_yy, g_z_y_yz_yz, g_z_y_yz_zz, g_z_yyy_yz_xx, g_z_yyy_yz_xy, g_z_yyy_yz_xz, g_z_yyy_yz_yy, g_z_yyy_yz_yz, g_z_yyy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yy_yz_xx[i] = -4.0 * g_z_y_yz_xx[i] * a_exp + 4.0 * g_z_yyy_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_yz_xy[i] = -4.0 * g_z_y_yz_xy[i] * a_exp + 4.0 * g_z_yyy_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_yz_xz[i] = -4.0 * g_z_y_yz_xz[i] * a_exp + 4.0 * g_z_yyy_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_yz_yy[i] = -4.0 * g_z_y_yz_yy[i] * a_exp + 4.0 * g_z_yyy_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_yz_yz[i] = -4.0 * g_z_y_yz_yz[i] * a_exp + 4.0 * g_z_yyy_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_yz_zz[i] = -4.0 * g_z_y_yz_zz[i] * a_exp + 4.0 * g_z_yyy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1650-1656)

    #pragma omp simd aligned(g_z_y_0_0_0_yy_zz_xx, g_z_y_0_0_0_yy_zz_xy, g_z_y_0_0_0_yy_zz_xz, g_z_y_0_0_0_yy_zz_yy, g_z_y_0_0_0_yy_zz_yz, g_z_y_0_0_0_yy_zz_zz, g_z_y_zz_xx, g_z_y_zz_xy, g_z_y_zz_xz, g_z_y_zz_yy, g_z_y_zz_yz, g_z_y_zz_zz, g_z_yyy_zz_xx, g_z_yyy_zz_xy, g_z_yyy_zz_xz, g_z_yyy_zz_yy, g_z_yyy_zz_yz, g_z_yyy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yy_zz_xx[i] = -4.0 * g_z_y_zz_xx[i] * a_exp + 4.0 * g_z_yyy_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_zz_xy[i] = -4.0 * g_z_y_zz_xy[i] * a_exp + 4.0 * g_z_yyy_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_zz_xz[i] = -4.0 * g_z_y_zz_xz[i] * a_exp + 4.0 * g_z_yyy_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_zz_yy[i] = -4.0 * g_z_y_zz_yy[i] * a_exp + 4.0 * g_z_yyy_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_zz_yz[i] = -4.0 * g_z_y_zz_yz[i] * a_exp + 4.0 * g_z_yyy_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_zz_zz[i] = -4.0 * g_z_y_zz_zz[i] * a_exp + 4.0 * g_z_yyy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1656-1662)

    #pragma omp simd aligned(g_z_y_0_0_0_yz_xx_xx, g_z_y_0_0_0_yz_xx_xy, g_z_y_0_0_0_yz_xx_xz, g_z_y_0_0_0_yz_xx_yy, g_z_y_0_0_0_yz_xx_yz, g_z_y_0_0_0_yz_xx_zz, g_z_yyz_xx_xx, g_z_yyz_xx_xy, g_z_yyz_xx_xz, g_z_yyz_xx_yy, g_z_yyz_xx_yz, g_z_yyz_xx_zz, g_z_z_xx_xx, g_z_z_xx_xy, g_z_z_xx_xz, g_z_z_xx_yy, g_z_z_xx_yz, g_z_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yz_xx_xx[i] = -2.0 * g_z_z_xx_xx[i] * a_exp + 4.0 * g_z_yyz_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_xx_xy[i] = -2.0 * g_z_z_xx_xy[i] * a_exp + 4.0 * g_z_yyz_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_xx_xz[i] = -2.0 * g_z_z_xx_xz[i] * a_exp + 4.0 * g_z_yyz_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_xx_yy[i] = -2.0 * g_z_z_xx_yy[i] * a_exp + 4.0 * g_z_yyz_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_xx_yz[i] = -2.0 * g_z_z_xx_yz[i] * a_exp + 4.0 * g_z_yyz_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_xx_zz[i] = -2.0 * g_z_z_xx_zz[i] * a_exp + 4.0 * g_z_yyz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1662-1668)

    #pragma omp simd aligned(g_z_y_0_0_0_yz_xy_xx, g_z_y_0_0_0_yz_xy_xy, g_z_y_0_0_0_yz_xy_xz, g_z_y_0_0_0_yz_xy_yy, g_z_y_0_0_0_yz_xy_yz, g_z_y_0_0_0_yz_xy_zz, g_z_yyz_xy_xx, g_z_yyz_xy_xy, g_z_yyz_xy_xz, g_z_yyz_xy_yy, g_z_yyz_xy_yz, g_z_yyz_xy_zz, g_z_z_xy_xx, g_z_z_xy_xy, g_z_z_xy_xz, g_z_z_xy_yy, g_z_z_xy_yz, g_z_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yz_xy_xx[i] = -2.0 * g_z_z_xy_xx[i] * a_exp + 4.0 * g_z_yyz_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_xy_xy[i] = -2.0 * g_z_z_xy_xy[i] * a_exp + 4.0 * g_z_yyz_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_xy_xz[i] = -2.0 * g_z_z_xy_xz[i] * a_exp + 4.0 * g_z_yyz_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_xy_yy[i] = -2.0 * g_z_z_xy_yy[i] * a_exp + 4.0 * g_z_yyz_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_xy_yz[i] = -2.0 * g_z_z_xy_yz[i] * a_exp + 4.0 * g_z_yyz_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_xy_zz[i] = -2.0 * g_z_z_xy_zz[i] * a_exp + 4.0 * g_z_yyz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1668-1674)

    #pragma omp simd aligned(g_z_y_0_0_0_yz_xz_xx, g_z_y_0_0_0_yz_xz_xy, g_z_y_0_0_0_yz_xz_xz, g_z_y_0_0_0_yz_xz_yy, g_z_y_0_0_0_yz_xz_yz, g_z_y_0_0_0_yz_xz_zz, g_z_yyz_xz_xx, g_z_yyz_xz_xy, g_z_yyz_xz_xz, g_z_yyz_xz_yy, g_z_yyz_xz_yz, g_z_yyz_xz_zz, g_z_z_xz_xx, g_z_z_xz_xy, g_z_z_xz_xz, g_z_z_xz_yy, g_z_z_xz_yz, g_z_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yz_xz_xx[i] = -2.0 * g_z_z_xz_xx[i] * a_exp + 4.0 * g_z_yyz_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_xz_xy[i] = -2.0 * g_z_z_xz_xy[i] * a_exp + 4.0 * g_z_yyz_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_xz_xz[i] = -2.0 * g_z_z_xz_xz[i] * a_exp + 4.0 * g_z_yyz_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_xz_yy[i] = -2.0 * g_z_z_xz_yy[i] * a_exp + 4.0 * g_z_yyz_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_xz_yz[i] = -2.0 * g_z_z_xz_yz[i] * a_exp + 4.0 * g_z_yyz_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_xz_zz[i] = -2.0 * g_z_z_xz_zz[i] * a_exp + 4.0 * g_z_yyz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1674-1680)

    #pragma omp simd aligned(g_z_y_0_0_0_yz_yy_xx, g_z_y_0_0_0_yz_yy_xy, g_z_y_0_0_0_yz_yy_xz, g_z_y_0_0_0_yz_yy_yy, g_z_y_0_0_0_yz_yy_yz, g_z_y_0_0_0_yz_yy_zz, g_z_yyz_yy_xx, g_z_yyz_yy_xy, g_z_yyz_yy_xz, g_z_yyz_yy_yy, g_z_yyz_yy_yz, g_z_yyz_yy_zz, g_z_z_yy_xx, g_z_z_yy_xy, g_z_z_yy_xz, g_z_z_yy_yy, g_z_z_yy_yz, g_z_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yz_yy_xx[i] = -2.0 * g_z_z_yy_xx[i] * a_exp + 4.0 * g_z_yyz_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_yy_xy[i] = -2.0 * g_z_z_yy_xy[i] * a_exp + 4.0 * g_z_yyz_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_yy_xz[i] = -2.0 * g_z_z_yy_xz[i] * a_exp + 4.0 * g_z_yyz_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_yy_yy[i] = -2.0 * g_z_z_yy_yy[i] * a_exp + 4.0 * g_z_yyz_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_yy_yz[i] = -2.0 * g_z_z_yy_yz[i] * a_exp + 4.0 * g_z_yyz_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_yy_zz[i] = -2.0 * g_z_z_yy_zz[i] * a_exp + 4.0 * g_z_yyz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1680-1686)

    #pragma omp simd aligned(g_z_y_0_0_0_yz_yz_xx, g_z_y_0_0_0_yz_yz_xy, g_z_y_0_0_0_yz_yz_xz, g_z_y_0_0_0_yz_yz_yy, g_z_y_0_0_0_yz_yz_yz, g_z_y_0_0_0_yz_yz_zz, g_z_yyz_yz_xx, g_z_yyz_yz_xy, g_z_yyz_yz_xz, g_z_yyz_yz_yy, g_z_yyz_yz_yz, g_z_yyz_yz_zz, g_z_z_yz_xx, g_z_z_yz_xy, g_z_z_yz_xz, g_z_z_yz_yy, g_z_z_yz_yz, g_z_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yz_yz_xx[i] = -2.0 * g_z_z_yz_xx[i] * a_exp + 4.0 * g_z_yyz_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_yz_xy[i] = -2.0 * g_z_z_yz_xy[i] * a_exp + 4.0 * g_z_yyz_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_yz_xz[i] = -2.0 * g_z_z_yz_xz[i] * a_exp + 4.0 * g_z_yyz_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_yz_yy[i] = -2.0 * g_z_z_yz_yy[i] * a_exp + 4.0 * g_z_yyz_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_yz_yz[i] = -2.0 * g_z_z_yz_yz[i] * a_exp + 4.0 * g_z_yyz_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_yz_zz[i] = -2.0 * g_z_z_yz_zz[i] * a_exp + 4.0 * g_z_yyz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1686-1692)

    #pragma omp simd aligned(g_z_y_0_0_0_yz_zz_xx, g_z_y_0_0_0_yz_zz_xy, g_z_y_0_0_0_yz_zz_xz, g_z_y_0_0_0_yz_zz_yy, g_z_y_0_0_0_yz_zz_yz, g_z_y_0_0_0_yz_zz_zz, g_z_yyz_zz_xx, g_z_yyz_zz_xy, g_z_yyz_zz_xz, g_z_yyz_zz_yy, g_z_yyz_zz_yz, g_z_yyz_zz_zz, g_z_z_zz_xx, g_z_z_zz_xy, g_z_z_zz_xz, g_z_z_zz_yy, g_z_z_zz_yz, g_z_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yz_zz_xx[i] = -2.0 * g_z_z_zz_xx[i] * a_exp + 4.0 * g_z_yyz_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_zz_xy[i] = -2.0 * g_z_z_zz_xy[i] * a_exp + 4.0 * g_z_yyz_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_zz_xz[i] = -2.0 * g_z_z_zz_xz[i] * a_exp + 4.0 * g_z_yyz_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_zz_yy[i] = -2.0 * g_z_z_zz_yy[i] * a_exp + 4.0 * g_z_yyz_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_zz_yz[i] = -2.0 * g_z_z_zz_yz[i] * a_exp + 4.0 * g_z_yyz_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_zz_zz[i] = -2.0 * g_z_z_zz_zz[i] * a_exp + 4.0 * g_z_yyz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1692-1698)

    #pragma omp simd aligned(g_z_y_0_0_0_zz_xx_xx, g_z_y_0_0_0_zz_xx_xy, g_z_y_0_0_0_zz_xx_xz, g_z_y_0_0_0_zz_xx_yy, g_z_y_0_0_0_zz_xx_yz, g_z_y_0_0_0_zz_xx_zz, g_z_yzz_xx_xx, g_z_yzz_xx_xy, g_z_yzz_xx_xz, g_z_yzz_xx_yy, g_z_yzz_xx_yz, g_z_yzz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_zz_xx_xx[i] = 4.0 * g_z_yzz_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_xx_xy[i] = 4.0 * g_z_yzz_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_xx_xz[i] = 4.0 * g_z_yzz_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_xx_yy[i] = 4.0 * g_z_yzz_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_xx_yz[i] = 4.0 * g_z_yzz_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_xx_zz[i] = 4.0 * g_z_yzz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1698-1704)

    #pragma omp simd aligned(g_z_y_0_0_0_zz_xy_xx, g_z_y_0_0_0_zz_xy_xy, g_z_y_0_0_0_zz_xy_xz, g_z_y_0_0_0_zz_xy_yy, g_z_y_0_0_0_zz_xy_yz, g_z_y_0_0_0_zz_xy_zz, g_z_yzz_xy_xx, g_z_yzz_xy_xy, g_z_yzz_xy_xz, g_z_yzz_xy_yy, g_z_yzz_xy_yz, g_z_yzz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_zz_xy_xx[i] = 4.0 * g_z_yzz_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_xy_xy[i] = 4.0 * g_z_yzz_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_xy_xz[i] = 4.0 * g_z_yzz_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_xy_yy[i] = 4.0 * g_z_yzz_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_xy_yz[i] = 4.0 * g_z_yzz_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_xy_zz[i] = 4.0 * g_z_yzz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1704-1710)

    #pragma omp simd aligned(g_z_y_0_0_0_zz_xz_xx, g_z_y_0_0_0_zz_xz_xy, g_z_y_0_0_0_zz_xz_xz, g_z_y_0_0_0_zz_xz_yy, g_z_y_0_0_0_zz_xz_yz, g_z_y_0_0_0_zz_xz_zz, g_z_yzz_xz_xx, g_z_yzz_xz_xy, g_z_yzz_xz_xz, g_z_yzz_xz_yy, g_z_yzz_xz_yz, g_z_yzz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_zz_xz_xx[i] = 4.0 * g_z_yzz_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_xz_xy[i] = 4.0 * g_z_yzz_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_xz_xz[i] = 4.0 * g_z_yzz_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_xz_yy[i] = 4.0 * g_z_yzz_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_xz_yz[i] = 4.0 * g_z_yzz_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_xz_zz[i] = 4.0 * g_z_yzz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1710-1716)

    #pragma omp simd aligned(g_z_y_0_0_0_zz_yy_xx, g_z_y_0_0_0_zz_yy_xy, g_z_y_0_0_0_zz_yy_xz, g_z_y_0_0_0_zz_yy_yy, g_z_y_0_0_0_zz_yy_yz, g_z_y_0_0_0_zz_yy_zz, g_z_yzz_yy_xx, g_z_yzz_yy_xy, g_z_yzz_yy_xz, g_z_yzz_yy_yy, g_z_yzz_yy_yz, g_z_yzz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_zz_yy_xx[i] = 4.0 * g_z_yzz_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_yy_xy[i] = 4.0 * g_z_yzz_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_yy_xz[i] = 4.0 * g_z_yzz_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_yy_yy[i] = 4.0 * g_z_yzz_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_yy_yz[i] = 4.0 * g_z_yzz_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_yy_zz[i] = 4.0 * g_z_yzz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1716-1722)

    #pragma omp simd aligned(g_z_y_0_0_0_zz_yz_xx, g_z_y_0_0_0_zz_yz_xy, g_z_y_0_0_0_zz_yz_xz, g_z_y_0_0_0_zz_yz_yy, g_z_y_0_0_0_zz_yz_yz, g_z_y_0_0_0_zz_yz_zz, g_z_yzz_yz_xx, g_z_yzz_yz_xy, g_z_yzz_yz_xz, g_z_yzz_yz_yy, g_z_yzz_yz_yz, g_z_yzz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_zz_yz_xx[i] = 4.0 * g_z_yzz_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_yz_xy[i] = 4.0 * g_z_yzz_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_yz_xz[i] = 4.0 * g_z_yzz_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_yz_yy[i] = 4.0 * g_z_yzz_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_yz_yz[i] = 4.0 * g_z_yzz_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_yz_zz[i] = 4.0 * g_z_yzz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1722-1728)

    #pragma omp simd aligned(g_z_y_0_0_0_zz_zz_xx, g_z_y_0_0_0_zz_zz_xy, g_z_y_0_0_0_zz_zz_xz, g_z_y_0_0_0_zz_zz_yy, g_z_y_0_0_0_zz_zz_yz, g_z_y_0_0_0_zz_zz_zz, g_z_yzz_zz_xx, g_z_yzz_zz_xy, g_z_yzz_zz_xz, g_z_yzz_zz_yy, g_z_yzz_zz_yz, g_z_yzz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_zz_zz_xx[i] = 4.0 * g_z_yzz_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_zz_xy[i] = 4.0 * g_z_yzz_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_zz_xz[i] = 4.0 * g_z_yzz_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_zz_yy[i] = 4.0 * g_z_yzz_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_zz_yz[i] = 4.0 * g_z_yzz_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_zz_zz[i] = 4.0 * g_z_yzz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1728-1734)

    #pragma omp simd aligned(g_z_xxz_xx_xx, g_z_xxz_xx_xy, g_z_xxz_xx_xz, g_z_xxz_xx_yy, g_z_xxz_xx_yz, g_z_xxz_xx_zz, g_z_z_0_0_0_xx_xx_xx, g_z_z_0_0_0_xx_xx_xy, g_z_z_0_0_0_xx_xx_xz, g_z_z_0_0_0_xx_xx_yy, g_z_z_0_0_0_xx_xx_yz, g_z_z_0_0_0_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xx_xx_xx[i] = 4.0 * g_z_xxz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_xx_xy[i] = 4.0 * g_z_xxz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_xx_xz[i] = 4.0 * g_z_xxz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_xx_yy[i] = 4.0 * g_z_xxz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_xx_yz[i] = 4.0 * g_z_xxz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_xx_zz[i] = 4.0 * g_z_xxz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1734-1740)

    #pragma omp simd aligned(g_z_xxz_xy_xx, g_z_xxz_xy_xy, g_z_xxz_xy_xz, g_z_xxz_xy_yy, g_z_xxz_xy_yz, g_z_xxz_xy_zz, g_z_z_0_0_0_xx_xy_xx, g_z_z_0_0_0_xx_xy_xy, g_z_z_0_0_0_xx_xy_xz, g_z_z_0_0_0_xx_xy_yy, g_z_z_0_0_0_xx_xy_yz, g_z_z_0_0_0_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xx_xy_xx[i] = 4.0 * g_z_xxz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_xy_xy[i] = 4.0 * g_z_xxz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_xy_xz[i] = 4.0 * g_z_xxz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_xy_yy[i] = 4.0 * g_z_xxz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_xy_yz[i] = 4.0 * g_z_xxz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_xy_zz[i] = 4.0 * g_z_xxz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1740-1746)

    #pragma omp simd aligned(g_z_xxz_xz_xx, g_z_xxz_xz_xy, g_z_xxz_xz_xz, g_z_xxz_xz_yy, g_z_xxz_xz_yz, g_z_xxz_xz_zz, g_z_z_0_0_0_xx_xz_xx, g_z_z_0_0_0_xx_xz_xy, g_z_z_0_0_0_xx_xz_xz, g_z_z_0_0_0_xx_xz_yy, g_z_z_0_0_0_xx_xz_yz, g_z_z_0_0_0_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xx_xz_xx[i] = 4.0 * g_z_xxz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_xz_xy[i] = 4.0 * g_z_xxz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_xz_xz[i] = 4.0 * g_z_xxz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_xz_yy[i] = 4.0 * g_z_xxz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_xz_yz[i] = 4.0 * g_z_xxz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_xz_zz[i] = 4.0 * g_z_xxz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1746-1752)

    #pragma omp simd aligned(g_z_xxz_yy_xx, g_z_xxz_yy_xy, g_z_xxz_yy_xz, g_z_xxz_yy_yy, g_z_xxz_yy_yz, g_z_xxz_yy_zz, g_z_z_0_0_0_xx_yy_xx, g_z_z_0_0_0_xx_yy_xy, g_z_z_0_0_0_xx_yy_xz, g_z_z_0_0_0_xx_yy_yy, g_z_z_0_0_0_xx_yy_yz, g_z_z_0_0_0_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xx_yy_xx[i] = 4.0 * g_z_xxz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_yy_xy[i] = 4.0 * g_z_xxz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_yy_xz[i] = 4.0 * g_z_xxz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_yy_yy[i] = 4.0 * g_z_xxz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_yy_yz[i] = 4.0 * g_z_xxz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_yy_zz[i] = 4.0 * g_z_xxz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1752-1758)

    #pragma omp simd aligned(g_z_xxz_yz_xx, g_z_xxz_yz_xy, g_z_xxz_yz_xz, g_z_xxz_yz_yy, g_z_xxz_yz_yz, g_z_xxz_yz_zz, g_z_z_0_0_0_xx_yz_xx, g_z_z_0_0_0_xx_yz_xy, g_z_z_0_0_0_xx_yz_xz, g_z_z_0_0_0_xx_yz_yy, g_z_z_0_0_0_xx_yz_yz, g_z_z_0_0_0_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xx_yz_xx[i] = 4.0 * g_z_xxz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_yz_xy[i] = 4.0 * g_z_xxz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_yz_xz[i] = 4.0 * g_z_xxz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_yz_yy[i] = 4.0 * g_z_xxz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_yz_yz[i] = 4.0 * g_z_xxz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_yz_zz[i] = 4.0 * g_z_xxz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1758-1764)

    #pragma omp simd aligned(g_z_xxz_zz_xx, g_z_xxz_zz_xy, g_z_xxz_zz_xz, g_z_xxz_zz_yy, g_z_xxz_zz_yz, g_z_xxz_zz_zz, g_z_z_0_0_0_xx_zz_xx, g_z_z_0_0_0_xx_zz_xy, g_z_z_0_0_0_xx_zz_xz, g_z_z_0_0_0_xx_zz_yy, g_z_z_0_0_0_xx_zz_yz, g_z_z_0_0_0_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xx_zz_xx[i] = 4.0 * g_z_xxz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_zz_xy[i] = 4.0 * g_z_xxz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_zz_xz[i] = 4.0 * g_z_xxz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_zz_yy[i] = 4.0 * g_z_xxz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_zz_yz[i] = 4.0 * g_z_xxz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_zz_zz[i] = 4.0 * g_z_xxz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1764-1770)

    #pragma omp simd aligned(g_z_xyz_xx_xx, g_z_xyz_xx_xy, g_z_xyz_xx_xz, g_z_xyz_xx_yy, g_z_xyz_xx_yz, g_z_xyz_xx_zz, g_z_z_0_0_0_xy_xx_xx, g_z_z_0_0_0_xy_xx_xy, g_z_z_0_0_0_xy_xx_xz, g_z_z_0_0_0_xy_xx_yy, g_z_z_0_0_0_xy_xx_yz, g_z_z_0_0_0_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xy_xx_xx[i] = 4.0 * g_z_xyz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_xx_xy[i] = 4.0 * g_z_xyz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_xx_xz[i] = 4.0 * g_z_xyz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_xx_yy[i] = 4.0 * g_z_xyz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_xx_yz[i] = 4.0 * g_z_xyz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_xx_zz[i] = 4.0 * g_z_xyz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1770-1776)

    #pragma omp simd aligned(g_z_xyz_xy_xx, g_z_xyz_xy_xy, g_z_xyz_xy_xz, g_z_xyz_xy_yy, g_z_xyz_xy_yz, g_z_xyz_xy_zz, g_z_z_0_0_0_xy_xy_xx, g_z_z_0_0_0_xy_xy_xy, g_z_z_0_0_0_xy_xy_xz, g_z_z_0_0_0_xy_xy_yy, g_z_z_0_0_0_xy_xy_yz, g_z_z_0_0_0_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xy_xy_xx[i] = 4.0 * g_z_xyz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_xy_xy[i] = 4.0 * g_z_xyz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_xy_xz[i] = 4.0 * g_z_xyz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_xy_yy[i] = 4.0 * g_z_xyz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_xy_yz[i] = 4.0 * g_z_xyz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_xy_zz[i] = 4.0 * g_z_xyz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1776-1782)

    #pragma omp simd aligned(g_z_xyz_xz_xx, g_z_xyz_xz_xy, g_z_xyz_xz_xz, g_z_xyz_xz_yy, g_z_xyz_xz_yz, g_z_xyz_xz_zz, g_z_z_0_0_0_xy_xz_xx, g_z_z_0_0_0_xy_xz_xy, g_z_z_0_0_0_xy_xz_xz, g_z_z_0_0_0_xy_xz_yy, g_z_z_0_0_0_xy_xz_yz, g_z_z_0_0_0_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xy_xz_xx[i] = 4.0 * g_z_xyz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_xz_xy[i] = 4.0 * g_z_xyz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_xz_xz[i] = 4.0 * g_z_xyz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_xz_yy[i] = 4.0 * g_z_xyz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_xz_yz[i] = 4.0 * g_z_xyz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_xz_zz[i] = 4.0 * g_z_xyz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1782-1788)

    #pragma omp simd aligned(g_z_xyz_yy_xx, g_z_xyz_yy_xy, g_z_xyz_yy_xz, g_z_xyz_yy_yy, g_z_xyz_yy_yz, g_z_xyz_yy_zz, g_z_z_0_0_0_xy_yy_xx, g_z_z_0_0_0_xy_yy_xy, g_z_z_0_0_0_xy_yy_xz, g_z_z_0_0_0_xy_yy_yy, g_z_z_0_0_0_xy_yy_yz, g_z_z_0_0_0_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xy_yy_xx[i] = 4.0 * g_z_xyz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_yy_xy[i] = 4.0 * g_z_xyz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_yy_xz[i] = 4.0 * g_z_xyz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_yy_yy[i] = 4.0 * g_z_xyz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_yy_yz[i] = 4.0 * g_z_xyz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_yy_zz[i] = 4.0 * g_z_xyz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1788-1794)

    #pragma omp simd aligned(g_z_xyz_yz_xx, g_z_xyz_yz_xy, g_z_xyz_yz_xz, g_z_xyz_yz_yy, g_z_xyz_yz_yz, g_z_xyz_yz_zz, g_z_z_0_0_0_xy_yz_xx, g_z_z_0_0_0_xy_yz_xy, g_z_z_0_0_0_xy_yz_xz, g_z_z_0_0_0_xy_yz_yy, g_z_z_0_0_0_xy_yz_yz, g_z_z_0_0_0_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xy_yz_xx[i] = 4.0 * g_z_xyz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_yz_xy[i] = 4.0 * g_z_xyz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_yz_xz[i] = 4.0 * g_z_xyz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_yz_yy[i] = 4.0 * g_z_xyz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_yz_yz[i] = 4.0 * g_z_xyz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_yz_zz[i] = 4.0 * g_z_xyz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1794-1800)

    #pragma omp simd aligned(g_z_xyz_zz_xx, g_z_xyz_zz_xy, g_z_xyz_zz_xz, g_z_xyz_zz_yy, g_z_xyz_zz_yz, g_z_xyz_zz_zz, g_z_z_0_0_0_xy_zz_xx, g_z_z_0_0_0_xy_zz_xy, g_z_z_0_0_0_xy_zz_xz, g_z_z_0_0_0_xy_zz_yy, g_z_z_0_0_0_xy_zz_yz, g_z_z_0_0_0_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xy_zz_xx[i] = 4.0 * g_z_xyz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_zz_xy[i] = 4.0 * g_z_xyz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_zz_xz[i] = 4.0 * g_z_xyz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_zz_yy[i] = 4.0 * g_z_xyz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_zz_yz[i] = 4.0 * g_z_xyz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_zz_zz[i] = 4.0 * g_z_xyz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1800-1806)

    #pragma omp simd aligned(g_z_x_xx_xx, g_z_x_xx_xy, g_z_x_xx_xz, g_z_x_xx_yy, g_z_x_xx_yz, g_z_x_xx_zz, g_z_xzz_xx_xx, g_z_xzz_xx_xy, g_z_xzz_xx_xz, g_z_xzz_xx_yy, g_z_xzz_xx_yz, g_z_xzz_xx_zz, g_z_z_0_0_0_xz_xx_xx, g_z_z_0_0_0_xz_xx_xy, g_z_z_0_0_0_xz_xx_xz, g_z_z_0_0_0_xz_xx_yy, g_z_z_0_0_0_xz_xx_yz, g_z_z_0_0_0_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xz_xx_xx[i] = -2.0 * g_z_x_xx_xx[i] * a_exp + 4.0 * g_z_xzz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_xx_xy[i] = -2.0 * g_z_x_xx_xy[i] * a_exp + 4.0 * g_z_xzz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_xx_xz[i] = -2.0 * g_z_x_xx_xz[i] * a_exp + 4.0 * g_z_xzz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_xx_yy[i] = -2.0 * g_z_x_xx_yy[i] * a_exp + 4.0 * g_z_xzz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_xx_yz[i] = -2.0 * g_z_x_xx_yz[i] * a_exp + 4.0 * g_z_xzz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_xx_zz[i] = -2.0 * g_z_x_xx_zz[i] * a_exp + 4.0 * g_z_xzz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1806-1812)

    #pragma omp simd aligned(g_z_x_xy_xx, g_z_x_xy_xy, g_z_x_xy_xz, g_z_x_xy_yy, g_z_x_xy_yz, g_z_x_xy_zz, g_z_xzz_xy_xx, g_z_xzz_xy_xy, g_z_xzz_xy_xz, g_z_xzz_xy_yy, g_z_xzz_xy_yz, g_z_xzz_xy_zz, g_z_z_0_0_0_xz_xy_xx, g_z_z_0_0_0_xz_xy_xy, g_z_z_0_0_0_xz_xy_xz, g_z_z_0_0_0_xz_xy_yy, g_z_z_0_0_0_xz_xy_yz, g_z_z_0_0_0_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xz_xy_xx[i] = -2.0 * g_z_x_xy_xx[i] * a_exp + 4.0 * g_z_xzz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_xy_xy[i] = -2.0 * g_z_x_xy_xy[i] * a_exp + 4.0 * g_z_xzz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_xy_xz[i] = -2.0 * g_z_x_xy_xz[i] * a_exp + 4.0 * g_z_xzz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_xy_yy[i] = -2.0 * g_z_x_xy_yy[i] * a_exp + 4.0 * g_z_xzz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_xy_yz[i] = -2.0 * g_z_x_xy_yz[i] * a_exp + 4.0 * g_z_xzz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_xy_zz[i] = -2.0 * g_z_x_xy_zz[i] * a_exp + 4.0 * g_z_xzz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1812-1818)

    #pragma omp simd aligned(g_z_x_xz_xx, g_z_x_xz_xy, g_z_x_xz_xz, g_z_x_xz_yy, g_z_x_xz_yz, g_z_x_xz_zz, g_z_xzz_xz_xx, g_z_xzz_xz_xy, g_z_xzz_xz_xz, g_z_xzz_xz_yy, g_z_xzz_xz_yz, g_z_xzz_xz_zz, g_z_z_0_0_0_xz_xz_xx, g_z_z_0_0_0_xz_xz_xy, g_z_z_0_0_0_xz_xz_xz, g_z_z_0_0_0_xz_xz_yy, g_z_z_0_0_0_xz_xz_yz, g_z_z_0_0_0_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xz_xz_xx[i] = -2.0 * g_z_x_xz_xx[i] * a_exp + 4.0 * g_z_xzz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_xz_xy[i] = -2.0 * g_z_x_xz_xy[i] * a_exp + 4.0 * g_z_xzz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_xz_xz[i] = -2.0 * g_z_x_xz_xz[i] * a_exp + 4.0 * g_z_xzz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_xz_yy[i] = -2.0 * g_z_x_xz_yy[i] * a_exp + 4.0 * g_z_xzz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_xz_yz[i] = -2.0 * g_z_x_xz_yz[i] * a_exp + 4.0 * g_z_xzz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_xz_zz[i] = -2.0 * g_z_x_xz_zz[i] * a_exp + 4.0 * g_z_xzz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1818-1824)

    #pragma omp simd aligned(g_z_x_yy_xx, g_z_x_yy_xy, g_z_x_yy_xz, g_z_x_yy_yy, g_z_x_yy_yz, g_z_x_yy_zz, g_z_xzz_yy_xx, g_z_xzz_yy_xy, g_z_xzz_yy_xz, g_z_xzz_yy_yy, g_z_xzz_yy_yz, g_z_xzz_yy_zz, g_z_z_0_0_0_xz_yy_xx, g_z_z_0_0_0_xz_yy_xy, g_z_z_0_0_0_xz_yy_xz, g_z_z_0_0_0_xz_yy_yy, g_z_z_0_0_0_xz_yy_yz, g_z_z_0_0_0_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xz_yy_xx[i] = -2.0 * g_z_x_yy_xx[i] * a_exp + 4.0 * g_z_xzz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_yy_xy[i] = -2.0 * g_z_x_yy_xy[i] * a_exp + 4.0 * g_z_xzz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_yy_xz[i] = -2.0 * g_z_x_yy_xz[i] * a_exp + 4.0 * g_z_xzz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_yy_yy[i] = -2.0 * g_z_x_yy_yy[i] * a_exp + 4.0 * g_z_xzz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_yy_yz[i] = -2.0 * g_z_x_yy_yz[i] * a_exp + 4.0 * g_z_xzz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_yy_zz[i] = -2.0 * g_z_x_yy_zz[i] * a_exp + 4.0 * g_z_xzz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1824-1830)

    #pragma omp simd aligned(g_z_x_yz_xx, g_z_x_yz_xy, g_z_x_yz_xz, g_z_x_yz_yy, g_z_x_yz_yz, g_z_x_yz_zz, g_z_xzz_yz_xx, g_z_xzz_yz_xy, g_z_xzz_yz_xz, g_z_xzz_yz_yy, g_z_xzz_yz_yz, g_z_xzz_yz_zz, g_z_z_0_0_0_xz_yz_xx, g_z_z_0_0_0_xz_yz_xy, g_z_z_0_0_0_xz_yz_xz, g_z_z_0_0_0_xz_yz_yy, g_z_z_0_0_0_xz_yz_yz, g_z_z_0_0_0_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xz_yz_xx[i] = -2.0 * g_z_x_yz_xx[i] * a_exp + 4.0 * g_z_xzz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_yz_xy[i] = -2.0 * g_z_x_yz_xy[i] * a_exp + 4.0 * g_z_xzz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_yz_xz[i] = -2.0 * g_z_x_yz_xz[i] * a_exp + 4.0 * g_z_xzz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_yz_yy[i] = -2.0 * g_z_x_yz_yy[i] * a_exp + 4.0 * g_z_xzz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_yz_yz[i] = -2.0 * g_z_x_yz_yz[i] * a_exp + 4.0 * g_z_xzz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_yz_zz[i] = -2.0 * g_z_x_yz_zz[i] * a_exp + 4.0 * g_z_xzz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1830-1836)

    #pragma omp simd aligned(g_z_x_zz_xx, g_z_x_zz_xy, g_z_x_zz_xz, g_z_x_zz_yy, g_z_x_zz_yz, g_z_x_zz_zz, g_z_xzz_zz_xx, g_z_xzz_zz_xy, g_z_xzz_zz_xz, g_z_xzz_zz_yy, g_z_xzz_zz_yz, g_z_xzz_zz_zz, g_z_z_0_0_0_xz_zz_xx, g_z_z_0_0_0_xz_zz_xy, g_z_z_0_0_0_xz_zz_xz, g_z_z_0_0_0_xz_zz_yy, g_z_z_0_0_0_xz_zz_yz, g_z_z_0_0_0_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xz_zz_xx[i] = -2.0 * g_z_x_zz_xx[i] * a_exp + 4.0 * g_z_xzz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_zz_xy[i] = -2.0 * g_z_x_zz_xy[i] * a_exp + 4.0 * g_z_xzz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_zz_xz[i] = -2.0 * g_z_x_zz_xz[i] * a_exp + 4.0 * g_z_xzz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_zz_yy[i] = -2.0 * g_z_x_zz_yy[i] * a_exp + 4.0 * g_z_xzz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_zz_yz[i] = -2.0 * g_z_x_zz_yz[i] * a_exp + 4.0 * g_z_xzz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_zz_zz[i] = -2.0 * g_z_x_zz_zz[i] * a_exp + 4.0 * g_z_xzz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1836-1842)

    #pragma omp simd aligned(g_z_yyz_xx_xx, g_z_yyz_xx_xy, g_z_yyz_xx_xz, g_z_yyz_xx_yy, g_z_yyz_xx_yz, g_z_yyz_xx_zz, g_z_z_0_0_0_yy_xx_xx, g_z_z_0_0_0_yy_xx_xy, g_z_z_0_0_0_yy_xx_xz, g_z_z_0_0_0_yy_xx_yy, g_z_z_0_0_0_yy_xx_yz, g_z_z_0_0_0_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yy_xx_xx[i] = 4.0 * g_z_yyz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_xx_xy[i] = 4.0 * g_z_yyz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_xx_xz[i] = 4.0 * g_z_yyz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_xx_yy[i] = 4.0 * g_z_yyz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_xx_yz[i] = 4.0 * g_z_yyz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_xx_zz[i] = 4.0 * g_z_yyz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1842-1848)

    #pragma omp simd aligned(g_z_yyz_xy_xx, g_z_yyz_xy_xy, g_z_yyz_xy_xz, g_z_yyz_xy_yy, g_z_yyz_xy_yz, g_z_yyz_xy_zz, g_z_z_0_0_0_yy_xy_xx, g_z_z_0_0_0_yy_xy_xy, g_z_z_0_0_0_yy_xy_xz, g_z_z_0_0_0_yy_xy_yy, g_z_z_0_0_0_yy_xy_yz, g_z_z_0_0_0_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yy_xy_xx[i] = 4.0 * g_z_yyz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_xy_xy[i] = 4.0 * g_z_yyz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_xy_xz[i] = 4.0 * g_z_yyz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_xy_yy[i] = 4.0 * g_z_yyz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_xy_yz[i] = 4.0 * g_z_yyz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_xy_zz[i] = 4.0 * g_z_yyz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1848-1854)

    #pragma omp simd aligned(g_z_yyz_xz_xx, g_z_yyz_xz_xy, g_z_yyz_xz_xz, g_z_yyz_xz_yy, g_z_yyz_xz_yz, g_z_yyz_xz_zz, g_z_z_0_0_0_yy_xz_xx, g_z_z_0_0_0_yy_xz_xy, g_z_z_0_0_0_yy_xz_xz, g_z_z_0_0_0_yy_xz_yy, g_z_z_0_0_0_yy_xz_yz, g_z_z_0_0_0_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yy_xz_xx[i] = 4.0 * g_z_yyz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_xz_xy[i] = 4.0 * g_z_yyz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_xz_xz[i] = 4.0 * g_z_yyz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_xz_yy[i] = 4.0 * g_z_yyz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_xz_yz[i] = 4.0 * g_z_yyz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_xz_zz[i] = 4.0 * g_z_yyz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1854-1860)

    #pragma omp simd aligned(g_z_yyz_yy_xx, g_z_yyz_yy_xy, g_z_yyz_yy_xz, g_z_yyz_yy_yy, g_z_yyz_yy_yz, g_z_yyz_yy_zz, g_z_z_0_0_0_yy_yy_xx, g_z_z_0_0_0_yy_yy_xy, g_z_z_0_0_0_yy_yy_xz, g_z_z_0_0_0_yy_yy_yy, g_z_z_0_0_0_yy_yy_yz, g_z_z_0_0_0_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yy_yy_xx[i] = 4.0 * g_z_yyz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_yy_xy[i] = 4.0 * g_z_yyz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_yy_xz[i] = 4.0 * g_z_yyz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_yy_yy[i] = 4.0 * g_z_yyz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_yy_yz[i] = 4.0 * g_z_yyz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_yy_zz[i] = 4.0 * g_z_yyz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1860-1866)

    #pragma omp simd aligned(g_z_yyz_yz_xx, g_z_yyz_yz_xy, g_z_yyz_yz_xz, g_z_yyz_yz_yy, g_z_yyz_yz_yz, g_z_yyz_yz_zz, g_z_z_0_0_0_yy_yz_xx, g_z_z_0_0_0_yy_yz_xy, g_z_z_0_0_0_yy_yz_xz, g_z_z_0_0_0_yy_yz_yy, g_z_z_0_0_0_yy_yz_yz, g_z_z_0_0_0_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yy_yz_xx[i] = 4.0 * g_z_yyz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_yz_xy[i] = 4.0 * g_z_yyz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_yz_xz[i] = 4.0 * g_z_yyz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_yz_yy[i] = 4.0 * g_z_yyz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_yz_yz[i] = 4.0 * g_z_yyz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_yz_zz[i] = 4.0 * g_z_yyz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1866-1872)

    #pragma omp simd aligned(g_z_yyz_zz_xx, g_z_yyz_zz_xy, g_z_yyz_zz_xz, g_z_yyz_zz_yy, g_z_yyz_zz_yz, g_z_yyz_zz_zz, g_z_z_0_0_0_yy_zz_xx, g_z_z_0_0_0_yy_zz_xy, g_z_z_0_0_0_yy_zz_xz, g_z_z_0_0_0_yy_zz_yy, g_z_z_0_0_0_yy_zz_yz, g_z_z_0_0_0_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yy_zz_xx[i] = 4.0 * g_z_yyz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_zz_xy[i] = 4.0 * g_z_yyz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_zz_xz[i] = 4.0 * g_z_yyz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_zz_yy[i] = 4.0 * g_z_yyz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_zz_yz[i] = 4.0 * g_z_yyz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_zz_zz[i] = 4.0 * g_z_yyz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1872-1878)

    #pragma omp simd aligned(g_z_y_xx_xx, g_z_y_xx_xy, g_z_y_xx_xz, g_z_y_xx_yy, g_z_y_xx_yz, g_z_y_xx_zz, g_z_yzz_xx_xx, g_z_yzz_xx_xy, g_z_yzz_xx_xz, g_z_yzz_xx_yy, g_z_yzz_xx_yz, g_z_yzz_xx_zz, g_z_z_0_0_0_yz_xx_xx, g_z_z_0_0_0_yz_xx_xy, g_z_z_0_0_0_yz_xx_xz, g_z_z_0_0_0_yz_xx_yy, g_z_z_0_0_0_yz_xx_yz, g_z_z_0_0_0_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yz_xx_xx[i] = -2.0 * g_z_y_xx_xx[i] * a_exp + 4.0 * g_z_yzz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_xx_xy[i] = -2.0 * g_z_y_xx_xy[i] * a_exp + 4.0 * g_z_yzz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_xx_xz[i] = -2.0 * g_z_y_xx_xz[i] * a_exp + 4.0 * g_z_yzz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_xx_yy[i] = -2.0 * g_z_y_xx_yy[i] * a_exp + 4.0 * g_z_yzz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_xx_yz[i] = -2.0 * g_z_y_xx_yz[i] * a_exp + 4.0 * g_z_yzz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_xx_zz[i] = -2.0 * g_z_y_xx_zz[i] * a_exp + 4.0 * g_z_yzz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1878-1884)

    #pragma omp simd aligned(g_z_y_xy_xx, g_z_y_xy_xy, g_z_y_xy_xz, g_z_y_xy_yy, g_z_y_xy_yz, g_z_y_xy_zz, g_z_yzz_xy_xx, g_z_yzz_xy_xy, g_z_yzz_xy_xz, g_z_yzz_xy_yy, g_z_yzz_xy_yz, g_z_yzz_xy_zz, g_z_z_0_0_0_yz_xy_xx, g_z_z_0_0_0_yz_xy_xy, g_z_z_0_0_0_yz_xy_xz, g_z_z_0_0_0_yz_xy_yy, g_z_z_0_0_0_yz_xy_yz, g_z_z_0_0_0_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yz_xy_xx[i] = -2.0 * g_z_y_xy_xx[i] * a_exp + 4.0 * g_z_yzz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_xy_xy[i] = -2.0 * g_z_y_xy_xy[i] * a_exp + 4.0 * g_z_yzz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_xy_xz[i] = -2.0 * g_z_y_xy_xz[i] * a_exp + 4.0 * g_z_yzz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_xy_yy[i] = -2.0 * g_z_y_xy_yy[i] * a_exp + 4.0 * g_z_yzz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_xy_yz[i] = -2.0 * g_z_y_xy_yz[i] * a_exp + 4.0 * g_z_yzz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_xy_zz[i] = -2.0 * g_z_y_xy_zz[i] * a_exp + 4.0 * g_z_yzz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1884-1890)

    #pragma omp simd aligned(g_z_y_xz_xx, g_z_y_xz_xy, g_z_y_xz_xz, g_z_y_xz_yy, g_z_y_xz_yz, g_z_y_xz_zz, g_z_yzz_xz_xx, g_z_yzz_xz_xy, g_z_yzz_xz_xz, g_z_yzz_xz_yy, g_z_yzz_xz_yz, g_z_yzz_xz_zz, g_z_z_0_0_0_yz_xz_xx, g_z_z_0_0_0_yz_xz_xy, g_z_z_0_0_0_yz_xz_xz, g_z_z_0_0_0_yz_xz_yy, g_z_z_0_0_0_yz_xz_yz, g_z_z_0_0_0_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yz_xz_xx[i] = -2.0 * g_z_y_xz_xx[i] * a_exp + 4.0 * g_z_yzz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_xz_xy[i] = -2.0 * g_z_y_xz_xy[i] * a_exp + 4.0 * g_z_yzz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_xz_xz[i] = -2.0 * g_z_y_xz_xz[i] * a_exp + 4.0 * g_z_yzz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_xz_yy[i] = -2.0 * g_z_y_xz_yy[i] * a_exp + 4.0 * g_z_yzz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_xz_yz[i] = -2.0 * g_z_y_xz_yz[i] * a_exp + 4.0 * g_z_yzz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_xz_zz[i] = -2.0 * g_z_y_xz_zz[i] * a_exp + 4.0 * g_z_yzz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1890-1896)

    #pragma omp simd aligned(g_z_y_yy_xx, g_z_y_yy_xy, g_z_y_yy_xz, g_z_y_yy_yy, g_z_y_yy_yz, g_z_y_yy_zz, g_z_yzz_yy_xx, g_z_yzz_yy_xy, g_z_yzz_yy_xz, g_z_yzz_yy_yy, g_z_yzz_yy_yz, g_z_yzz_yy_zz, g_z_z_0_0_0_yz_yy_xx, g_z_z_0_0_0_yz_yy_xy, g_z_z_0_0_0_yz_yy_xz, g_z_z_0_0_0_yz_yy_yy, g_z_z_0_0_0_yz_yy_yz, g_z_z_0_0_0_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yz_yy_xx[i] = -2.0 * g_z_y_yy_xx[i] * a_exp + 4.0 * g_z_yzz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_yy_xy[i] = -2.0 * g_z_y_yy_xy[i] * a_exp + 4.0 * g_z_yzz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_yy_xz[i] = -2.0 * g_z_y_yy_xz[i] * a_exp + 4.0 * g_z_yzz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_yy_yy[i] = -2.0 * g_z_y_yy_yy[i] * a_exp + 4.0 * g_z_yzz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_yy_yz[i] = -2.0 * g_z_y_yy_yz[i] * a_exp + 4.0 * g_z_yzz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_yy_zz[i] = -2.0 * g_z_y_yy_zz[i] * a_exp + 4.0 * g_z_yzz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1896-1902)

    #pragma omp simd aligned(g_z_y_yz_xx, g_z_y_yz_xy, g_z_y_yz_xz, g_z_y_yz_yy, g_z_y_yz_yz, g_z_y_yz_zz, g_z_yzz_yz_xx, g_z_yzz_yz_xy, g_z_yzz_yz_xz, g_z_yzz_yz_yy, g_z_yzz_yz_yz, g_z_yzz_yz_zz, g_z_z_0_0_0_yz_yz_xx, g_z_z_0_0_0_yz_yz_xy, g_z_z_0_0_0_yz_yz_xz, g_z_z_0_0_0_yz_yz_yy, g_z_z_0_0_0_yz_yz_yz, g_z_z_0_0_0_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yz_yz_xx[i] = -2.0 * g_z_y_yz_xx[i] * a_exp + 4.0 * g_z_yzz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_yz_xy[i] = -2.0 * g_z_y_yz_xy[i] * a_exp + 4.0 * g_z_yzz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_yz_xz[i] = -2.0 * g_z_y_yz_xz[i] * a_exp + 4.0 * g_z_yzz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_yz_yy[i] = -2.0 * g_z_y_yz_yy[i] * a_exp + 4.0 * g_z_yzz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_yz_yz[i] = -2.0 * g_z_y_yz_yz[i] * a_exp + 4.0 * g_z_yzz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_yz_zz[i] = -2.0 * g_z_y_yz_zz[i] * a_exp + 4.0 * g_z_yzz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1902-1908)

    #pragma omp simd aligned(g_z_y_zz_xx, g_z_y_zz_xy, g_z_y_zz_xz, g_z_y_zz_yy, g_z_y_zz_yz, g_z_y_zz_zz, g_z_yzz_zz_xx, g_z_yzz_zz_xy, g_z_yzz_zz_xz, g_z_yzz_zz_yy, g_z_yzz_zz_yz, g_z_yzz_zz_zz, g_z_z_0_0_0_yz_zz_xx, g_z_z_0_0_0_yz_zz_xy, g_z_z_0_0_0_yz_zz_xz, g_z_z_0_0_0_yz_zz_yy, g_z_z_0_0_0_yz_zz_yz, g_z_z_0_0_0_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yz_zz_xx[i] = -2.0 * g_z_y_zz_xx[i] * a_exp + 4.0 * g_z_yzz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_zz_xy[i] = -2.0 * g_z_y_zz_xy[i] * a_exp + 4.0 * g_z_yzz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_zz_xz[i] = -2.0 * g_z_y_zz_xz[i] * a_exp + 4.0 * g_z_yzz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_zz_yy[i] = -2.0 * g_z_y_zz_yy[i] * a_exp + 4.0 * g_z_yzz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_zz_yz[i] = -2.0 * g_z_y_zz_yz[i] * a_exp + 4.0 * g_z_yzz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_zz_zz[i] = -2.0 * g_z_y_zz_zz[i] * a_exp + 4.0 * g_z_yzz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1908-1914)

    #pragma omp simd aligned(g_z_z_0_0_0_zz_xx_xx, g_z_z_0_0_0_zz_xx_xy, g_z_z_0_0_0_zz_xx_xz, g_z_z_0_0_0_zz_xx_yy, g_z_z_0_0_0_zz_xx_yz, g_z_z_0_0_0_zz_xx_zz, g_z_z_xx_xx, g_z_z_xx_xy, g_z_z_xx_xz, g_z_z_xx_yy, g_z_z_xx_yz, g_z_z_xx_zz, g_z_zzz_xx_xx, g_z_zzz_xx_xy, g_z_zzz_xx_xz, g_z_zzz_xx_yy, g_z_zzz_xx_yz, g_z_zzz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_zz_xx_xx[i] = -4.0 * g_z_z_xx_xx[i] * a_exp + 4.0 * g_z_zzz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_xx_xy[i] = -4.0 * g_z_z_xx_xy[i] * a_exp + 4.0 * g_z_zzz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_xx_xz[i] = -4.0 * g_z_z_xx_xz[i] * a_exp + 4.0 * g_z_zzz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_xx_yy[i] = -4.0 * g_z_z_xx_yy[i] * a_exp + 4.0 * g_z_zzz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_xx_yz[i] = -4.0 * g_z_z_xx_yz[i] * a_exp + 4.0 * g_z_zzz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_xx_zz[i] = -4.0 * g_z_z_xx_zz[i] * a_exp + 4.0 * g_z_zzz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1914-1920)

    #pragma omp simd aligned(g_z_z_0_0_0_zz_xy_xx, g_z_z_0_0_0_zz_xy_xy, g_z_z_0_0_0_zz_xy_xz, g_z_z_0_0_0_zz_xy_yy, g_z_z_0_0_0_zz_xy_yz, g_z_z_0_0_0_zz_xy_zz, g_z_z_xy_xx, g_z_z_xy_xy, g_z_z_xy_xz, g_z_z_xy_yy, g_z_z_xy_yz, g_z_z_xy_zz, g_z_zzz_xy_xx, g_z_zzz_xy_xy, g_z_zzz_xy_xz, g_z_zzz_xy_yy, g_z_zzz_xy_yz, g_z_zzz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_zz_xy_xx[i] = -4.0 * g_z_z_xy_xx[i] * a_exp + 4.0 * g_z_zzz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_xy_xy[i] = -4.0 * g_z_z_xy_xy[i] * a_exp + 4.0 * g_z_zzz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_xy_xz[i] = -4.0 * g_z_z_xy_xz[i] * a_exp + 4.0 * g_z_zzz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_xy_yy[i] = -4.0 * g_z_z_xy_yy[i] * a_exp + 4.0 * g_z_zzz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_xy_yz[i] = -4.0 * g_z_z_xy_yz[i] * a_exp + 4.0 * g_z_zzz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_xy_zz[i] = -4.0 * g_z_z_xy_zz[i] * a_exp + 4.0 * g_z_zzz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1920-1926)

    #pragma omp simd aligned(g_z_z_0_0_0_zz_xz_xx, g_z_z_0_0_0_zz_xz_xy, g_z_z_0_0_0_zz_xz_xz, g_z_z_0_0_0_zz_xz_yy, g_z_z_0_0_0_zz_xz_yz, g_z_z_0_0_0_zz_xz_zz, g_z_z_xz_xx, g_z_z_xz_xy, g_z_z_xz_xz, g_z_z_xz_yy, g_z_z_xz_yz, g_z_z_xz_zz, g_z_zzz_xz_xx, g_z_zzz_xz_xy, g_z_zzz_xz_xz, g_z_zzz_xz_yy, g_z_zzz_xz_yz, g_z_zzz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_zz_xz_xx[i] = -4.0 * g_z_z_xz_xx[i] * a_exp + 4.0 * g_z_zzz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_xz_xy[i] = -4.0 * g_z_z_xz_xy[i] * a_exp + 4.0 * g_z_zzz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_xz_xz[i] = -4.0 * g_z_z_xz_xz[i] * a_exp + 4.0 * g_z_zzz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_xz_yy[i] = -4.0 * g_z_z_xz_yy[i] * a_exp + 4.0 * g_z_zzz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_xz_yz[i] = -4.0 * g_z_z_xz_yz[i] * a_exp + 4.0 * g_z_zzz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_xz_zz[i] = -4.0 * g_z_z_xz_zz[i] * a_exp + 4.0 * g_z_zzz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1926-1932)

    #pragma omp simd aligned(g_z_z_0_0_0_zz_yy_xx, g_z_z_0_0_0_zz_yy_xy, g_z_z_0_0_0_zz_yy_xz, g_z_z_0_0_0_zz_yy_yy, g_z_z_0_0_0_zz_yy_yz, g_z_z_0_0_0_zz_yy_zz, g_z_z_yy_xx, g_z_z_yy_xy, g_z_z_yy_xz, g_z_z_yy_yy, g_z_z_yy_yz, g_z_z_yy_zz, g_z_zzz_yy_xx, g_z_zzz_yy_xy, g_z_zzz_yy_xz, g_z_zzz_yy_yy, g_z_zzz_yy_yz, g_z_zzz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_zz_yy_xx[i] = -4.0 * g_z_z_yy_xx[i] * a_exp + 4.0 * g_z_zzz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_yy_xy[i] = -4.0 * g_z_z_yy_xy[i] * a_exp + 4.0 * g_z_zzz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_yy_xz[i] = -4.0 * g_z_z_yy_xz[i] * a_exp + 4.0 * g_z_zzz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_yy_yy[i] = -4.0 * g_z_z_yy_yy[i] * a_exp + 4.0 * g_z_zzz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_yy_yz[i] = -4.0 * g_z_z_yy_yz[i] * a_exp + 4.0 * g_z_zzz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_yy_zz[i] = -4.0 * g_z_z_yy_zz[i] * a_exp + 4.0 * g_z_zzz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1932-1938)

    #pragma omp simd aligned(g_z_z_0_0_0_zz_yz_xx, g_z_z_0_0_0_zz_yz_xy, g_z_z_0_0_0_zz_yz_xz, g_z_z_0_0_0_zz_yz_yy, g_z_z_0_0_0_zz_yz_yz, g_z_z_0_0_0_zz_yz_zz, g_z_z_yz_xx, g_z_z_yz_xy, g_z_z_yz_xz, g_z_z_yz_yy, g_z_z_yz_yz, g_z_z_yz_zz, g_z_zzz_yz_xx, g_z_zzz_yz_xy, g_z_zzz_yz_xz, g_z_zzz_yz_yy, g_z_zzz_yz_yz, g_z_zzz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_zz_yz_xx[i] = -4.0 * g_z_z_yz_xx[i] * a_exp + 4.0 * g_z_zzz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_yz_xy[i] = -4.0 * g_z_z_yz_xy[i] * a_exp + 4.0 * g_z_zzz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_yz_xz[i] = -4.0 * g_z_z_yz_xz[i] * a_exp + 4.0 * g_z_zzz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_yz_yy[i] = -4.0 * g_z_z_yz_yy[i] * a_exp + 4.0 * g_z_zzz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_yz_yz[i] = -4.0 * g_z_z_yz_yz[i] * a_exp + 4.0 * g_z_zzz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_yz_zz[i] = -4.0 * g_z_z_yz_zz[i] * a_exp + 4.0 * g_z_zzz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1938-1944)

    #pragma omp simd aligned(g_z_z_0_0_0_zz_zz_xx, g_z_z_0_0_0_zz_zz_xy, g_z_z_0_0_0_zz_zz_xz, g_z_z_0_0_0_zz_zz_yy, g_z_z_0_0_0_zz_zz_yz, g_z_z_0_0_0_zz_zz_zz, g_z_z_zz_xx, g_z_z_zz_xy, g_z_z_zz_xz, g_z_z_zz_yy, g_z_z_zz_yz, g_z_z_zz_zz, g_z_zzz_zz_xx, g_z_zzz_zz_xy, g_z_zzz_zz_xz, g_z_zzz_zz_yy, g_z_zzz_zz_yz, g_z_zzz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_zz_zz_xx[i] = -4.0 * g_z_z_zz_xx[i] * a_exp + 4.0 * g_z_zzz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_zz_xy[i] = -4.0 * g_z_z_zz_xy[i] * a_exp + 4.0 * g_z_zzz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_zz_xz[i] = -4.0 * g_z_z_zz_xz[i] * a_exp + 4.0 * g_z_zzz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_zz_yy[i] = -4.0 * g_z_z_zz_yy[i] * a_exp + 4.0 * g_z_zzz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_zz_yz[i] = -4.0 * g_z_z_zz_yz[i] * a_exp + 4.0 * g_z_zzz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_zz_zz[i] = -4.0 * g_z_z_zz_zz[i] * a_exp + 4.0 * g_z_zzz_zz_zz[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

