#include "GeomDeriv1010OfScalarForSPPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_sppd_0(CSimdArray<double>& buffer_1010_sppd,
                     const CSimdArray<double>& buffer_ppsd,
                     const CSimdArray<double>& buffer_ppdd,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_sppd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_ppsd

    auto g_x_x_0_xx = buffer_ppsd[0];

    auto g_x_x_0_xy = buffer_ppsd[1];

    auto g_x_x_0_xz = buffer_ppsd[2];

    auto g_x_x_0_yy = buffer_ppsd[3];

    auto g_x_x_0_yz = buffer_ppsd[4];

    auto g_x_x_0_zz = buffer_ppsd[5];

    auto g_x_y_0_xx = buffer_ppsd[6];

    auto g_x_y_0_xy = buffer_ppsd[7];

    auto g_x_y_0_xz = buffer_ppsd[8];

    auto g_x_y_0_yy = buffer_ppsd[9];

    auto g_x_y_0_yz = buffer_ppsd[10];

    auto g_x_y_0_zz = buffer_ppsd[11];

    auto g_x_z_0_xx = buffer_ppsd[12];

    auto g_x_z_0_xy = buffer_ppsd[13];

    auto g_x_z_0_xz = buffer_ppsd[14];

    auto g_x_z_0_yy = buffer_ppsd[15];

    auto g_x_z_0_yz = buffer_ppsd[16];

    auto g_x_z_0_zz = buffer_ppsd[17];

    auto g_y_x_0_xx = buffer_ppsd[18];

    auto g_y_x_0_xy = buffer_ppsd[19];

    auto g_y_x_0_xz = buffer_ppsd[20];

    auto g_y_x_0_yy = buffer_ppsd[21];

    auto g_y_x_0_yz = buffer_ppsd[22];

    auto g_y_x_0_zz = buffer_ppsd[23];

    auto g_y_y_0_xx = buffer_ppsd[24];

    auto g_y_y_0_xy = buffer_ppsd[25];

    auto g_y_y_0_xz = buffer_ppsd[26];

    auto g_y_y_0_yy = buffer_ppsd[27];

    auto g_y_y_0_yz = buffer_ppsd[28];

    auto g_y_y_0_zz = buffer_ppsd[29];

    auto g_y_z_0_xx = buffer_ppsd[30];

    auto g_y_z_0_xy = buffer_ppsd[31];

    auto g_y_z_0_xz = buffer_ppsd[32];

    auto g_y_z_0_yy = buffer_ppsd[33];

    auto g_y_z_0_yz = buffer_ppsd[34];

    auto g_y_z_0_zz = buffer_ppsd[35];

    auto g_z_x_0_xx = buffer_ppsd[36];

    auto g_z_x_0_xy = buffer_ppsd[37];

    auto g_z_x_0_xz = buffer_ppsd[38];

    auto g_z_x_0_yy = buffer_ppsd[39];

    auto g_z_x_0_yz = buffer_ppsd[40];

    auto g_z_x_0_zz = buffer_ppsd[41];

    auto g_z_y_0_xx = buffer_ppsd[42];

    auto g_z_y_0_xy = buffer_ppsd[43];

    auto g_z_y_0_xz = buffer_ppsd[44];

    auto g_z_y_0_yy = buffer_ppsd[45];

    auto g_z_y_0_yz = buffer_ppsd[46];

    auto g_z_y_0_zz = buffer_ppsd[47];

    auto g_z_z_0_xx = buffer_ppsd[48];

    auto g_z_z_0_xy = buffer_ppsd[49];

    auto g_z_z_0_xz = buffer_ppsd[50];

    auto g_z_z_0_yy = buffer_ppsd[51];

    auto g_z_z_0_yz = buffer_ppsd[52];

    auto g_z_z_0_zz = buffer_ppsd[53];

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

    /// Set up components of integrals buffer : buffer_1010_sppd

    auto g_x_0_x_0_0_x_x_xx = buffer_1010_sppd[0];

    auto g_x_0_x_0_0_x_x_xy = buffer_1010_sppd[1];

    auto g_x_0_x_0_0_x_x_xz = buffer_1010_sppd[2];

    auto g_x_0_x_0_0_x_x_yy = buffer_1010_sppd[3];

    auto g_x_0_x_0_0_x_x_yz = buffer_1010_sppd[4];

    auto g_x_0_x_0_0_x_x_zz = buffer_1010_sppd[5];

    auto g_x_0_x_0_0_x_y_xx = buffer_1010_sppd[6];

    auto g_x_0_x_0_0_x_y_xy = buffer_1010_sppd[7];

    auto g_x_0_x_0_0_x_y_xz = buffer_1010_sppd[8];

    auto g_x_0_x_0_0_x_y_yy = buffer_1010_sppd[9];

    auto g_x_0_x_0_0_x_y_yz = buffer_1010_sppd[10];

    auto g_x_0_x_0_0_x_y_zz = buffer_1010_sppd[11];

    auto g_x_0_x_0_0_x_z_xx = buffer_1010_sppd[12];

    auto g_x_0_x_0_0_x_z_xy = buffer_1010_sppd[13];

    auto g_x_0_x_0_0_x_z_xz = buffer_1010_sppd[14];

    auto g_x_0_x_0_0_x_z_yy = buffer_1010_sppd[15];

    auto g_x_0_x_0_0_x_z_yz = buffer_1010_sppd[16];

    auto g_x_0_x_0_0_x_z_zz = buffer_1010_sppd[17];

    auto g_x_0_x_0_0_y_x_xx = buffer_1010_sppd[18];

    auto g_x_0_x_0_0_y_x_xy = buffer_1010_sppd[19];

    auto g_x_0_x_0_0_y_x_xz = buffer_1010_sppd[20];

    auto g_x_0_x_0_0_y_x_yy = buffer_1010_sppd[21];

    auto g_x_0_x_0_0_y_x_yz = buffer_1010_sppd[22];

    auto g_x_0_x_0_0_y_x_zz = buffer_1010_sppd[23];

    auto g_x_0_x_0_0_y_y_xx = buffer_1010_sppd[24];

    auto g_x_0_x_0_0_y_y_xy = buffer_1010_sppd[25];

    auto g_x_0_x_0_0_y_y_xz = buffer_1010_sppd[26];

    auto g_x_0_x_0_0_y_y_yy = buffer_1010_sppd[27];

    auto g_x_0_x_0_0_y_y_yz = buffer_1010_sppd[28];

    auto g_x_0_x_0_0_y_y_zz = buffer_1010_sppd[29];

    auto g_x_0_x_0_0_y_z_xx = buffer_1010_sppd[30];

    auto g_x_0_x_0_0_y_z_xy = buffer_1010_sppd[31];

    auto g_x_0_x_0_0_y_z_xz = buffer_1010_sppd[32];

    auto g_x_0_x_0_0_y_z_yy = buffer_1010_sppd[33];

    auto g_x_0_x_0_0_y_z_yz = buffer_1010_sppd[34];

    auto g_x_0_x_0_0_y_z_zz = buffer_1010_sppd[35];

    auto g_x_0_x_0_0_z_x_xx = buffer_1010_sppd[36];

    auto g_x_0_x_0_0_z_x_xy = buffer_1010_sppd[37];

    auto g_x_0_x_0_0_z_x_xz = buffer_1010_sppd[38];

    auto g_x_0_x_0_0_z_x_yy = buffer_1010_sppd[39];

    auto g_x_0_x_0_0_z_x_yz = buffer_1010_sppd[40];

    auto g_x_0_x_0_0_z_x_zz = buffer_1010_sppd[41];

    auto g_x_0_x_0_0_z_y_xx = buffer_1010_sppd[42];

    auto g_x_0_x_0_0_z_y_xy = buffer_1010_sppd[43];

    auto g_x_0_x_0_0_z_y_xz = buffer_1010_sppd[44];

    auto g_x_0_x_0_0_z_y_yy = buffer_1010_sppd[45];

    auto g_x_0_x_0_0_z_y_yz = buffer_1010_sppd[46];

    auto g_x_0_x_0_0_z_y_zz = buffer_1010_sppd[47];

    auto g_x_0_x_0_0_z_z_xx = buffer_1010_sppd[48];

    auto g_x_0_x_0_0_z_z_xy = buffer_1010_sppd[49];

    auto g_x_0_x_0_0_z_z_xz = buffer_1010_sppd[50];

    auto g_x_0_x_0_0_z_z_yy = buffer_1010_sppd[51];

    auto g_x_0_x_0_0_z_z_yz = buffer_1010_sppd[52];

    auto g_x_0_x_0_0_z_z_zz = buffer_1010_sppd[53];

    auto g_x_0_y_0_0_x_x_xx = buffer_1010_sppd[54];

    auto g_x_0_y_0_0_x_x_xy = buffer_1010_sppd[55];

    auto g_x_0_y_0_0_x_x_xz = buffer_1010_sppd[56];

    auto g_x_0_y_0_0_x_x_yy = buffer_1010_sppd[57];

    auto g_x_0_y_0_0_x_x_yz = buffer_1010_sppd[58];

    auto g_x_0_y_0_0_x_x_zz = buffer_1010_sppd[59];

    auto g_x_0_y_0_0_x_y_xx = buffer_1010_sppd[60];

    auto g_x_0_y_0_0_x_y_xy = buffer_1010_sppd[61];

    auto g_x_0_y_0_0_x_y_xz = buffer_1010_sppd[62];

    auto g_x_0_y_0_0_x_y_yy = buffer_1010_sppd[63];

    auto g_x_0_y_0_0_x_y_yz = buffer_1010_sppd[64];

    auto g_x_0_y_0_0_x_y_zz = buffer_1010_sppd[65];

    auto g_x_0_y_0_0_x_z_xx = buffer_1010_sppd[66];

    auto g_x_0_y_0_0_x_z_xy = buffer_1010_sppd[67];

    auto g_x_0_y_0_0_x_z_xz = buffer_1010_sppd[68];

    auto g_x_0_y_0_0_x_z_yy = buffer_1010_sppd[69];

    auto g_x_0_y_0_0_x_z_yz = buffer_1010_sppd[70];

    auto g_x_0_y_0_0_x_z_zz = buffer_1010_sppd[71];

    auto g_x_0_y_0_0_y_x_xx = buffer_1010_sppd[72];

    auto g_x_0_y_0_0_y_x_xy = buffer_1010_sppd[73];

    auto g_x_0_y_0_0_y_x_xz = buffer_1010_sppd[74];

    auto g_x_0_y_0_0_y_x_yy = buffer_1010_sppd[75];

    auto g_x_0_y_0_0_y_x_yz = buffer_1010_sppd[76];

    auto g_x_0_y_0_0_y_x_zz = buffer_1010_sppd[77];

    auto g_x_0_y_0_0_y_y_xx = buffer_1010_sppd[78];

    auto g_x_0_y_0_0_y_y_xy = buffer_1010_sppd[79];

    auto g_x_0_y_0_0_y_y_xz = buffer_1010_sppd[80];

    auto g_x_0_y_0_0_y_y_yy = buffer_1010_sppd[81];

    auto g_x_0_y_0_0_y_y_yz = buffer_1010_sppd[82];

    auto g_x_0_y_0_0_y_y_zz = buffer_1010_sppd[83];

    auto g_x_0_y_0_0_y_z_xx = buffer_1010_sppd[84];

    auto g_x_0_y_0_0_y_z_xy = buffer_1010_sppd[85];

    auto g_x_0_y_0_0_y_z_xz = buffer_1010_sppd[86];

    auto g_x_0_y_0_0_y_z_yy = buffer_1010_sppd[87];

    auto g_x_0_y_0_0_y_z_yz = buffer_1010_sppd[88];

    auto g_x_0_y_0_0_y_z_zz = buffer_1010_sppd[89];

    auto g_x_0_y_0_0_z_x_xx = buffer_1010_sppd[90];

    auto g_x_0_y_0_0_z_x_xy = buffer_1010_sppd[91];

    auto g_x_0_y_0_0_z_x_xz = buffer_1010_sppd[92];

    auto g_x_0_y_0_0_z_x_yy = buffer_1010_sppd[93];

    auto g_x_0_y_0_0_z_x_yz = buffer_1010_sppd[94];

    auto g_x_0_y_0_0_z_x_zz = buffer_1010_sppd[95];

    auto g_x_0_y_0_0_z_y_xx = buffer_1010_sppd[96];

    auto g_x_0_y_0_0_z_y_xy = buffer_1010_sppd[97];

    auto g_x_0_y_0_0_z_y_xz = buffer_1010_sppd[98];

    auto g_x_0_y_0_0_z_y_yy = buffer_1010_sppd[99];

    auto g_x_0_y_0_0_z_y_yz = buffer_1010_sppd[100];

    auto g_x_0_y_0_0_z_y_zz = buffer_1010_sppd[101];

    auto g_x_0_y_0_0_z_z_xx = buffer_1010_sppd[102];

    auto g_x_0_y_0_0_z_z_xy = buffer_1010_sppd[103];

    auto g_x_0_y_0_0_z_z_xz = buffer_1010_sppd[104];

    auto g_x_0_y_0_0_z_z_yy = buffer_1010_sppd[105];

    auto g_x_0_y_0_0_z_z_yz = buffer_1010_sppd[106];

    auto g_x_0_y_0_0_z_z_zz = buffer_1010_sppd[107];

    auto g_x_0_z_0_0_x_x_xx = buffer_1010_sppd[108];

    auto g_x_0_z_0_0_x_x_xy = buffer_1010_sppd[109];

    auto g_x_0_z_0_0_x_x_xz = buffer_1010_sppd[110];

    auto g_x_0_z_0_0_x_x_yy = buffer_1010_sppd[111];

    auto g_x_0_z_0_0_x_x_yz = buffer_1010_sppd[112];

    auto g_x_0_z_0_0_x_x_zz = buffer_1010_sppd[113];

    auto g_x_0_z_0_0_x_y_xx = buffer_1010_sppd[114];

    auto g_x_0_z_0_0_x_y_xy = buffer_1010_sppd[115];

    auto g_x_0_z_0_0_x_y_xz = buffer_1010_sppd[116];

    auto g_x_0_z_0_0_x_y_yy = buffer_1010_sppd[117];

    auto g_x_0_z_0_0_x_y_yz = buffer_1010_sppd[118];

    auto g_x_0_z_0_0_x_y_zz = buffer_1010_sppd[119];

    auto g_x_0_z_0_0_x_z_xx = buffer_1010_sppd[120];

    auto g_x_0_z_0_0_x_z_xy = buffer_1010_sppd[121];

    auto g_x_0_z_0_0_x_z_xz = buffer_1010_sppd[122];

    auto g_x_0_z_0_0_x_z_yy = buffer_1010_sppd[123];

    auto g_x_0_z_0_0_x_z_yz = buffer_1010_sppd[124];

    auto g_x_0_z_0_0_x_z_zz = buffer_1010_sppd[125];

    auto g_x_0_z_0_0_y_x_xx = buffer_1010_sppd[126];

    auto g_x_0_z_0_0_y_x_xy = buffer_1010_sppd[127];

    auto g_x_0_z_0_0_y_x_xz = buffer_1010_sppd[128];

    auto g_x_0_z_0_0_y_x_yy = buffer_1010_sppd[129];

    auto g_x_0_z_0_0_y_x_yz = buffer_1010_sppd[130];

    auto g_x_0_z_0_0_y_x_zz = buffer_1010_sppd[131];

    auto g_x_0_z_0_0_y_y_xx = buffer_1010_sppd[132];

    auto g_x_0_z_0_0_y_y_xy = buffer_1010_sppd[133];

    auto g_x_0_z_0_0_y_y_xz = buffer_1010_sppd[134];

    auto g_x_0_z_0_0_y_y_yy = buffer_1010_sppd[135];

    auto g_x_0_z_0_0_y_y_yz = buffer_1010_sppd[136];

    auto g_x_0_z_0_0_y_y_zz = buffer_1010_sppd[137];

    auto g_x_0_z_0_0_y_z_xx = buffer_1010_sppd[138];

    auto g_x_0_z_0_0_y_z_xy = buffer_1010_sppd[139];

    auto g_x_0_z_0_0_y_z_xz = buffer_1010_sppd[140];

    auto g_x_0_z_0_0_y_z_yy = buffer_1010_sppd[141];

    auto g_x_0_z_0_0_y_z_yz = buffer_1010_sppd[142];

    auto g_x_0_z_0_0_y_z_zz = buffer_1010_sppd[143];

    auto g_x_0_z_0_0_z_x_xx = buffer_1010_sppd[144];

    auto g_x_0_z_0_0_z_x_xy = buffer_1010_sppd[145];

    auto g_x_0_z_0_0_z_x_xz = buffer_1010_sppd[146];

    auto g_x_0_z_0_0_z_x_yy = buffer_1010_sppd[147];

    auto g_x_0_z_0_0_z_x_yz = buffer_1010_sppd[148];

    auto g_x_0_z_0_0_z_x_zz = buffer_1010_sppd[149];

    auto g_x_0_z_0_0_z_y_xx = buffer_1010_sppd[150];

    auto g_x_0_z_0_0_z_y_xy = buffer_1010_sppd[151];

    auto g_x_0_z_0_0_z_y_xz = buffer_1010_sppd[152];

    auto g_x_0_z_0_0_z_y_yy = buffer_1010_sppd[153];

    auto g_x_0_z_0_0_z_y_yz = buffer_1010_sppd[154];

    auto g_x_0_z_0_0_z_y_zz = buffer_1010_sppd[155];

    auto g_x_0_z_0_0_z_z_xx = buffer_1010_sppd[156];

    auto g_x_0_z_0_0_z_z_xy = buffer_1010_sppd[157];

    auto g_x_0_z_0_0_z_z_xz = buffer_1010_sppd[158];

    auto g_x_0_z_0_0_z_z_yy = buffer_1010_sppd[159];

    auto g_x_0_z_0_0_z_z_yz = buffer_1010_sppd[160];

    auto g_x_0_z_0_0_z_z_zz = buffer_1010_sppd[161];

    auto g_y_0_x_0_0_x_x_xx = buffer_1010_sppd[162];

    auto g_y_0_x_0_0_x_x_xy = buffer_1010_sppd[163];

    auto g_y_0_x_0_0_x_x_xz = buffer_1010_sppd[164];

    auto g_y_0_x_0_0_x_x_yy = buffer_1010_sppd[165];

    auto g_y_0_x_0_0_x_x_yz = buffer_1010_sppd[166];

    auto g_y_0_x_0_0_x_x_zz = buffer_1010_sppd[167];

    auto g_y_0_x_0_0_x_y_xx = buffer_1010_sppd[168];

    auto g_y_0_x_0_0_x_y_xy = buffer_1010_sppd[169];

    auto g_y_0_x_0_0_x_y_xz = buffer_1010_sppd[170];

    auto g_y_0_x_0_0_x_y_yy = buffer_1010_sppd[171];

    auto g_y_0_x_0_0_x_y_yz = buffer_1010_sppd[172];

    auto g_y_0_x_0_0_x_y_zz = buffer_1010_sppd[173];

    auto g_y_0_x_0_0_x_z_xx = buffer_1010_sppd[174];

    auto g_y_0_x_0_0_x_z_xy = buffer_1010_sppd[175];

    auto g_y_0_x_0_0_x_z_xz = buffer_1010_sppd[176];

    auto g_y_0_x_0_0_x_z_yy = buffer_1010_sppd[177];

    auto g_y_0_x_0_0_x_z_yz = buffer_1010_sppd[178];

    auto g_y_0_x_0_0_x_z_zz = buffer_1010_sppd[179];

    auto g_y_0_x_0_0_y_x_xx = buffer_1010_sppd[180];

    auto g_y_0_x_0_0_y_x_xy = buffer_1010_sppd[181];

    auto g_y_0_x_0_0_y_x_xz = buffer_1010_sppd[182];

    auto g_y_0_x_0_0_y_x_yy = buffer_1010_sppd[183];

    auto g_y_0_x_0_0_y_x_yz = buffer_1010_sppd[184];

    auto g_y_0_x_0_0_y_x_zz = buffer_1010_sppd[185];

    auto g_y_0_x_0_0_y_y_xx = buffer_1010_sppd[186];

    auto g_y_0_x_0_0_y_y_xy = buffer_1010_sppd[187];

    auto g_y_0_x_0_0_y_y_xz = buffer_1010_sppd[188];

    auto g_y_0_x_0_0_y_y_yy = buffer_1010_sppd[189];

    auto g_y_0_x_0_0_y_y_yz = buffer_1010_sppd[190];

    auto g_y_0_x_0_0_y_y_zz = buffer_1010_sppd[191];

    auto g_y_0_x_0_0_y_z_xx = buffer_1010_sppd[192];

    auto g_y_0_x_0_0_y_z_xy = buffer_1010_sppd[193];

    auto g_y_0_x_0_0_y_z_xz = buffer_1010_sppd[194];

    auto g_y_0_x_0_0_y_z_yy = buffer_1010_sppd[195];

    auto g_y_0_x_0_0_y_z_yz = buffer_1010_sppd[196];

    auto g_y_0_x_0_0_y_z_zz = buffer_1010_sppd[197];

    auto g_y_0_x_0_0_z_x_xx = buffer_1010_sppd[198];

    auto g_y_0_x_0_0_z_x_xy = buffer_1010_sppd[199];

    auto g_y_0_x_0_0_z_x_xz = buffer_1010_sppd[200];

    auto g_y_0_x_0_0_z_x_yy = buffer_1010_sppd[201];

    auto g_y_0_x_0_0_z_x_yz = buffer_1010_sppd[202];

    auto g_y_0_x_0_0_z_x_zz = buffer_1010_sppd[203];

    auto g_y_0_x_0_0_z_y_xx = buffer_1010_sppd[204];

    auto g_y_0_x_0_0_z_y_xy = buffer_1010_sppd[205];

    auto g_y_0_x_0_0_z_y_xz = buffer_1010_sppd[206];

    auto g_y_0_x_0_0_z_y_yy = buffer_1010_sppd[207];

    auto g_y_0_x_0_0_z_y_yz = buffer_1010_sppd[208];

    auto g_y_0_x_0_0_z_y_zz = buffer_1010_sppd[209];

    auto g_y_0_x_0_0_z_z_xx = buffer_1010_sppd[210];

    auto g_y_0_x_0_0_z_z_xy = buffer_1010_sppd[211];

    auto g_y_0_x_0_0_z_z_xz = buffer_1010_sppd[212];

    auto g_y_0_x_0_0_z_z_yy = buffer_1010_sppd[213];

    auto g_y_0_x_0_0_z_z_yz = buffer_1010_sppd[214];

    auto g_y_0_x_0_0_z_z_zz = buffer_1010_sppd[215];

    auto g_y_0_y_0_0_x_x_xx = buffer_1010_sppd[216];

    auto g_y_0_y_0_0_x_x_xy = buffer_1010_sppd[217];

    auto g_y_0_y_0_0_x_x_xz = buffer_1010_sppd[218];

    auto g_y_0_y_0_0_x_x_yy = buffer_1010_sppd[219];

    auto g_y_0_y_0_0_x_x_yz = buffer_1010_sppd[220];

    auto g_y_0_y_0_0_x_x_zz = buffer_1010_sppd[221];

    auto g_y_0_y_0_0_x_y_xx = buffer_1010_sppd[222];

    auto g_y_0_y_0_0_x_y_xy = buffer_1010_sppd[223];

    auto g_y_0_y_0_0_x_y_xz = buffer_1010_sppd[224];

    auto g_y_0_y_0_0_x_y_yy = buffer_1010_sppd[225];

    auto g_y_0_y_0_0_x_y_yz = buffer_1010_sppd[226];

    auto g_y_0_y_0_0_x_y_zz = buffer_1010_sppd[227];

    auto g_y_0_y_0_0_x_z_xx = buffer_1010_sppd[228];

    auto g_y_0_y_0_0_x_z_xy = buffer_1010_sppd[229];

    auto g_y_0_y_0_0_x_z_xz = buffer_1010_sppd[230];

    auto g_y_0_y_0_0_x_z_yy = buffer_1010_sppd[231];

    auto g_y_0_y_0_0_x_z_yz = buffer_1010_sppd[232];

    auto g_y_0_y_0_0_x_z_zz = buffer_1010_sppd[233];

    auto g_y_0_y_0_0_y_x_xx = buffer_1010_sppd[234];

    auto g_y_0_y_0_0_y_x_xy = buffer_1010_sppd[235];

    auto g_y_0_y_0_0_y_x_xz = buffer_1010_sppd[236];

    auto g_y_0_y_0_0_y_x_yy = buffer_1010_sppd[237];

    auto g_y_0_y_0_0_y_x_yz = buffer_1010_sppd[238];

    auto g_y_0_y_0_0_y_x_zz = buffer_1010_sppd[239];

    auto g_y_0_y_0_0_y_y_xx = buffer_1010_sppd[240];

    auto g_y_0_y_0_0_y_y_xy = buffer_1010_sppd[241];

    auto g_y_0_y_0_0_y_y_xz = buffer_1010_sppd[242];

    auto g_y_0_y_0_0_y_y_yy = buffer_1010_sppd[243];

    auto g_y_0_y_0_0_y_y_yz = buffer_1010_sppd[244];

    auto g_y_0_y_0_0_y_y_zz = buffer_1010_sppd[245];

    auto g_y_0_y_0_0_y_z_xx = buffer_1010_sppd[246];

    auto g_y_0_y_0_0_y_z_xy = buffer_1010_sppd[247];

    auto g_y_0_y_0_0_y_z_xz = buffer_1010_sppd[248];

    auto g_y_0_y_0_0_y_z_yy = buffer_1010_sppd[249];

    auto g_y_0_y_0_0_y_z_yz = buffer_1010_sppd[250];

    auto g_y_0_y_0_0_y_z_zz = buffer_1010_sppd[251];

    auto g_y_0_y_0_0_z_x_xx = buffer_1010_sppd[252];

    auto g_y_0_y_0_0_z_x_xy = buffer_1010_sppd[253];

    auto g_y_0_y_0_0_z_x_xz = buffer_1010_sppd[254];

    auto g_y_0_y_0_0_z_x_yy = buffer_1010_sppd[255];

    auto g_y_0_y_0_0_z_x_yz = buffer_1010_sppd[256];

    auto g_y_0_y_0_0_z_x_zz = buffer_1010_sppd[257];

    auto g_y_0_y_0_0_z_y_xx = buffer_1010_sppd[258];

    auto g_y_0_y_0_0_z_y_xy = buffer_1010_sppd[259];

    auto g_y_0_y_0_0_z_y_xz = buffer_1010_sppd[260];

    auto g_y_0_y_0_0_z_y_yy = buffer_1010_sppd[261];

    auto g_y_0_y_0_0_z_y_yz = buffer_1010_sppd[262];

    auto g_y_0_y_0_0_z_y_zz = buffer_1010_sppd[263];

    auto g_y_0_y_0_0_z_z_xx = buffer_1010_sppd[264];

    auto g_y_0_y_0_0_z_z_xy = buffer_1010_sppd[265];

    auto g_y_0_y_0_0_z_z_xz = buffer_1010_sppd[266];

    auto g_y_0_y_0_0_z_z_yy = buffer_1010_sppd[267];

    auto g_y_0_y_0_0_z_z_yz = buffer_1010_sppd[268];

    auto g_y_0_y_0_0_z_z_zz = buffer_1010_sppd[269];

    auto g_y_0_z_0_0_x_x_xx = buffer_1010_sppd[270];

    auto g_y_0_z_0_0_x_x_xy = buffer_1010_sppd[271];

    auto g_y_0_z_0_0_x_x_xz = buffer_1010_sppd[272];

    auto g_y_0_z_0_0_x_x_yy = buffer_1010_sppd[273];

    auto g_y_0_z_0_0_x_x_yz = buffer_1010_sppd[274];

    auto g_y_0_z_0_0_x_x_zz = buffer_1010_sppd[275];

    auto g_y_0_z_0_0_x_y_xx = buffer_1010_sppd[276];

    auto g_y_0_z_0_0_x_y_xy = buffer_1010_sppd[277];

    auto g_y_0_z_0_0_x_y_xz = buffer_1010_sppd[278];

    auto g_y_0_z_0_0_x_y_yy = buffer_1010_sppd[279];

    auto g_y_0_z_0_0_x_y_yz = buffer_1010_sppd[280];

    auto g_y_0_z_0_0_x_y_zz = buffer_1010_sppd[281];

    auto g_y_0_z_0_0_x_z_xx = buffer_1010_sppd[282];

    auto g_y_0_z_0_0_x_z_xy = buffer_1010_sppd[283];

    auto g_y_0_z_0_0_x_z_xz = buffer_1010_sppd[284];

    auto g_y_0_z_0_0_x_z_yy = buffer_1010_sppd[285];

    auto g_y_0_z_0_0_x_z_yz = buffer_1010_sppd[286];

    auto g_y_0_z_0_0_x_z_zz = buffer_1010_sppd[287];

    auto g_y_0_z_0_0_y_x_xx = buffer_1010_sppd[288];

    auto g_y_0_z_0_0_y_x_xy = buffer_1010_sppd[289];

    auto g_y_0_z_0_0_y_x_xz = buffer_1010_sppd[290];

    auto g_y_0_z_0_0_y_x_yy = buffer_1010_sppd[291];

    auto g_y_0_z_0_0_y_x_yz = buffer_1010_sppd[292];

    auto g_y_0_z_0_0_y_x_zz = buffer_1010_sppd[293];

    auto g_y_0_z_0_0_y_y_xx = buffer_1010_sppd[294];

    auto g_y_0_z_0_0_y_y_xy = buffer_1010_sppd[295];

    auto g_y_0_z_0_0_y_y_xz = buffer_1010_sppd[296];

    auto g_y_0_z_0_0_y_y_yy = buffer_1010_sppd[297];

    auto g_y_0_z_0_0_y_y_yz = buffer_1010_sppd[298];

    auto g_y_0_z_0_0_y_y_zz = buffer_1010_sppd[299];

    auto g_y_0_z_0_0_y_z_xx = buffer_1010_sppd[300];

    auto g_y_0_z_0_0_y_z_xy = buffer_1010_sppd[301];

    auto g_y_0_z_0_0_y_z_xz = buffer_1010_sppd[302];

    auto g_y_0_z_0_0_y_z_yy = buffer_1010_sppd[303];

    auto g_y_0_z_0_0_y_z_yz = buffer_1010_sppd[304];

    auto g_y_0_z_0_0_y_z_zz = buffer_1010_sppd[305];

    auto g_y_0_z_0_0_z_x_xx = buffer_1010_sppd[306];

    auto g_y_0_z_0_0_z_x_xy = buffer_1010_sppd[307];

    auto g_y_0_z_0_0_z_x_xz = buffer_1010_sppd[308];

    auto g_y_0_z_0_0_z_x_yy = buffer_1010_sppd[309];

    auto g_y_0_z_0_0_z_x_yz = buffer_1010_sppd[310];

    auto g_y_0_z_0_0_z_x_zz = buffer_1010_sppd[311];

    auto g_y_0_z_0_0_z_y_xx = buffer_1010_sppd[312];

    auto g_y_0_z_0_0_z_y_xy = buffer_1010_sppd[313];

    auto g_y_0_z_0_0_z_y_xz = buffer_1010_sppd[314];

    auto g_y_0_z_0_0_z_y_yy = buffer_1010_sppd[315];

    auto g_y_0_z_0_0_z_y_yz = buffer_1010_sppd[316];

    auto g_y_0_z_0_0_z_y_zz = buffer_1010_sppd[317];

    auto g_y_0_z_0_0_z_z_xx = buffer_1010_sppd[318];

    auto g_y_0_z_0_0_z_z_xy = buffer_1010_sppd[319];

    auto g_y_0_z_0_0_z_z_xz = buffer_1010_sppd[320];

    auto g_y_0_z_0_0_z_z_yy = buffer_1010_sppd[321];

    auto g_y_0_z_0_0_z_z_yz = buffer_1010_sppd[322];

    auto g_y_0_z_0_0_z_z_zz = buffer_1010_sppd[323];

    auto g_z_0_x_0_0_x_x_xx = buffer_1010_sppd[324];

    auto g_z_0_x_0_0_x_x_xy = buffer_1010_sppd[325];

    auto g_z_0_x_0_0_x_x_xz = buffer_1010_sppd[326];

    auto g_z_0_x_0_0_x_x_yy = buffer_1010_sppd[327];

    auto g_z_0_x_0_0_x_x_yz = buffer_1010_sppd[328];

    auto g_z_0_x_0_0_x_x_zz = buffer_1010_sppd[329];

    auto g_z_0_x_0_0_x_y_xx = buffer_1010_sppd[330];

    auto g_z_0_x_0_0_x_y_xy = buffer_1010_sppd[331];

    auto g_z_0_x_0_0_x_y_xz = buffer_1010_sppd[332];

    auto g_z_0_x_0_0_x_y_yy = buffer_1010_sppd[333];

    auto g_z_0_x_0_0_x_y_yz = buffer_1010_sppd[334];

    auto g_z_0_x_0_0_x_y_zz = buffer_1010_sppd[335];

    auto g_z_0_x_0_0_x_z_xx = buffer_1010_sppd[336];

    auto g_z_0_x_0_0_x_z_xy = buffer_1010_sppd[337];

    auto g_z_0_x_0_0_x_z_xz = buffer_1010_sppd[338];

    auto g_z_0_x_0_0_x_z_yy = buffer_1010_sppd[339];

    auto g_z_0_x_0_0_x_z_yz = buffer_1010_sppd[340];

    auto g_z_0_x_0_0_x_z_zz = buffer_1010_sppd[341];

    auto g_z_0_x_0_0_y_x_xx = buffer_1010_sppd[342];

    auto g_z_0_x_0_0_y_x_xy = buffer_1010_sppd[343];

    auto g_z_0_x_0_0_y_x_xz = buffer_1010_sppd[344];

    auto g_z_0_x_0_0_y_x_yy = buffer_1010_sppd[345];

    auto g_z_0_x_0_0_y_x_yz = buffer_1010_sppd[346];

    auto g_z_0_x_0_0_y_x_zz = buffer_1010_sppd[347];

    auto g_z_0_x_0_0_y_y_xx = buffer_1010_sppd[348];

    auto g_z_0_x_0_0_y_y_xy = buffer_1010_sppd[349];

    auto g_z_0_x_0_0_y_y_xz = buffer_1010_sppd[350];

    auto g_z_0_x_0_0_y_y_yy = buffer_1010_sppd[351];

    auto g_z_0_x_0_0_y_y_yz = buffer_1010_sppd[352];

    auto g_z_0_x_0_0_y_y_zz = buffer_1010_sppd[353];

    auto g_z_0_x_0_0_y_z_xx = buffer_1010_sppd[354];

    auto g_z_0_x_0_0_y_z_xy = buffer_1010_sppd[355];

    auto g_z_0_x_0_0_y_z_xz = buffer_1010_sppd[356];

    auto g_z_0_x_0_0_y_z_yy = buffer_1010_sppd[357];

    auto g_z_0_x_0_0_y_z_yz = buffer_1010_sppd[358];

    auto g_z_0_x_0_0_y_z_zz = buffer_1010_sppd[359];

    auto g_z_0_x_0_0_z_x_xx = buffer_1010_sppd[360];

    auto g_z_0_x_0_0_z_x_xy = buffer_1010_sppd[361];

    auto g_z_0_x_0_0_z_x_xz = buffer_1010_sppd[362];

    auto g_z_0_x_0_0_z_x_yy = buffer_1010_sppd[363];

    auto g_z_0_x_0_0_z_x_yz = buffer_1010_sppd[364];

    auto g_z_0_x_0_0_z_x_zz = buffer_1010_sppd[365];

    auto g_z_0_x_0_0_z_y_xx = buffer_1010_sppd[366];

    auto g_z_0_x_0_0_z_y_xy = buffer_1010_sppd[367];

    auto g_z_0_x_0_0_z_y_xz = buffer_1010_sppd[368];

    auto g_z_0_x_0_0_z_y_yy = buffer_1010_sppd[369];

    auto g_z_0_x_0_0_z_y_yz = buffer_1010_sppd[370];

    auto g_z_0_x_0_0_z_y_zz = buffer_1010_sppd[371];

    auto g_z_0_x_0_0_z_z_xx = buffer_1010_sppd[372];

    auto g_z_0_x_0_0_z_z_xy = buffer_1010_sppd[373];

    auto g_z_0_x_0_0_z_z_xz = buffer_1010_sppd[374];

    auto g_z_0_x_0_0_z_z_yy = buffer_1010_sppd[375];

    auto g_z_0_x_0_0_z_z_yz = buffer_1010_sppd[376];

    auto g_z_0_x_0_0_z_z_zz = buffer_1010_sppd[377];

    auto g_z_0_y_0_0_x_x_xx = buffer_1010_sppd[378];

    auto g_z_0_y_0_0_x_x_xy = buffer_1010_sppd[379];

    auto g_z_0_y_0_0_x_x_xz = buffer_1010_sppd[380];

    auto g_z_0_y_0_0_x_x_yy = buffer_1010_sppd[381];

    auto g_z_0_y_0_0_x_x_yz = buffer_1010_sppd[382];

    auto g_z_0_y_0_0_x_x_zz = buffer_1010_sppd[383];

    auto g_z_0_y_0_0_x_y_xx = buffer_1010_sppd[384];

    auto g_z_0_y_0_0_x_y_xy = buffer_1010_sppd[385];

    auto g_z_0_y_0_0_x_y_xz = buffer_1010_sppd[386];

    auto g_z_0_y_0_0_x_y_yy = buffer_1010_sppd[387];

    auto g_z_0_y_0_0_x_y_yz = buffer_1010_sppd[388];

    auto g_z_0_y_0_0_x_y_zz = buffer_1010_sppd[389];

    auto g_z_0_y_0_0_x_z_xx = buffer_1010_sppd[390];

    auto g_z_0_y_0_0_x_z_xy = buffer_1010_sppd[391];

    auto g_z_0_y_0_0_x_z_xz = buffer_1010_sppd[392];

    auto g_z_0_y_0_0_x_z_yy = buffer_1010_sppd[393];

    auto g_z_0_y_0_0_x_z_yz = buffer_1010_sppd[394];

    auto g_z_0_y_0_0_x_z_zz = buffer_1010_sppd[395];

    auto g_z_0_y_0_0_y_x_xx = buffer_1010_sppd[396];

    auto g_z_0_y_0_0_y_x_xy = buffer_1010_sppd[397];

    auto g_z_0_y_0_0_y_x_xz = buffer_1010_sppd[398];

    auto g_z_0_y_0_0_y_x_yy = buffer_1010_sppd[399];

    auto g_z_0_y_0_0_y_x_yz = buffer_1010_sppd[400];

    auto g_z_0_y_0_0_y_x_zz = buffer_1010_sppd[401];

    auto g_z_0_y_0_0_y_y_xx = buffer_1010_sppd[402];

    auto g_z_0_y_0_0_y_y_xy = buffer_1010_sppd[403];

    auto g_z_0_y_0_0_y_y_xz = buffer_1010_sppd[404];

    auto g_z_0_y_0_0_y_y_yy = buffer_1010_sppd[405];

    auto g_z_0_y_0_0_y_y_yz = buffer_1010_sppd[406];

    auto g_z_0_y_0_0_y_y_zz = buffer_1010_sppd[407];

    auto g_z_0_y_0_0_y_z_xx = buffer_1010_sppd[408];

    auto g_z_0_y_0_0_y_z_xy = buffer_1010_sppd[409];

    auto g_z_0_y_0_0_y_z_xz = buffer_1010_sppd[410];

    auto g_z_0_y_0_0_y_z_yy = buffer_1010_sppd[411];

    auto g_z_0_y_0_0_y_z_yz = buffer_1010_sppd[412];

    auto g_z_0_y_0_0_y_z_zz = buffer_1010_sppd[413];

    auto g_z_0_y_0_0_z_x_xx = buffer_1010_sppd[414];

    auto g_z_0_y_0_0_z_x_xy = buffer_1010_sppd[415];

    auto g_z_0_y_0_0_z_x_xz = buffer_1010_sppd[416];

    auto g_z_0_y_0_0_z_x_yy = buffer_1010_sppd[417];

    auto g_z_0_y_0_0_z_x_yz = buffer_1010_sppd[418];

    auto g_z_0_y_0_0_z_x_zz = buffer_1010_sppd[419];

    auto g_z_0_y_0_0_z_y_xx = buffer_1010_sppd[420];

    auto g_z_0_y_0_0_z_y_xy = buffer_1010_sppd[421];

    auto g_z_0_y_0_0_z_y_xz = buffer_1010_sppd[422];

    auto g_z_0_y_0_0_z_y_yy = buffer_1010_sppd[423];

    auto g_z_0_y_0_0_z_y_yz = buffer_1010_sppd[424];

    auto g_z_0_y_0_0_z_y_zz = buffer_1010_sppd[425];

    auto g_z_0_y_0_0_z_z_xx = buffer_1010_sppd[426];

    auto g_z_0_y_0_0_z_z_xy = buffer_1010_sppd[427];

    auto g_z_0_y_0_0_z_z_xz = buffer_1010_sppd[428];

    auto g_z_0_y_0_0_z_z_yy = buffer_1010_sppd[429];

    auto g_z_0_y_0_0_z_z_yz = buffer_1010_sppd[430];

    auto g_z_0_y_0_0_z_z_zz = buffer_1010_sppd[431];

    auto g_z_0_z_0_0_x_x_xx = buffer_1010_sppd[432];

    auto g_z_0_z_0_0_x_x_xy = buffer_1010_sppd[433];

    auto g_z_0_z_0_0_x_x_xz = buffer_1010_sppd[434];

    auto g_z_0_z_0_0_x_x_yy = buffer_1010_sppd[435];

    auto g_z_0_z_0_0_x_x_yz = buffer_1010_sppd[436];

    auto g_z_0_z_0_0_x_x_zz = buffer_1010_sppd[437];

    auto g_z_0_z_0_0_x_y_xx = buffer_1010_sppd[438];

    auto g_z_0_z_0_0_x_y_xy = buffer_1010_sppd[439];

    auto g_z_0_z_0_0_x_y_xz = buffer_1010_sppd[440];

    auto g_z_0_z_0_0_x_y_yy = buffer_1010_sppd[441];

    auto g_z_0_z_0_0_x_y_yz = buffer_1010_sppd[442];

    auto g_z_0_z_0_0_x_y_zz = buffer_1010_sppd[443];

    auto g_z_0_z_0_0_x_z_xx = buffer_1010_sppd[444];

    auto g_z_0_z_0_0_x_z_xy = buffer_1010_sppd[445];

    auto g_z_0_z_0_0_x_z_xz = buffer_1010_sppd[446];

    auto g_z_0_z_0_0_x_z_yy = buffer_1010_sppd[447];

    auto g_z_0_z_0_0_x_z_yz = buffer_1010_sppd[448];

    auto g_z_0_z_0_0_x_z_zz = buffer_1010_sppd[449];

    auto g_z_0_z_0_0_y_x_xx = buffer_1010_sppd[450];

    auto g_z_0_z_0_0_y_x_xy = buffer_1010_sppd[451];

    auto g_z_0_z_0_0_y_x_xz = buffer_1010_sppd[452];

    auto g_z_0_z_0_0_y_x_yy = buffer_1010_sppd[453];

    auto g_z_0_z_0_0_y_x_yz = buffer_1010_sppd[454];

    auto g_z_0_z_0_0_y_x_zz = buffer_1010_sppd[455];

    auto g_z_0_z_0_0_y_y_xx = buffer_1010_sppd[456];

    auto g_z_0_z_0_0_y_y_xy = buffer_1010_sppd[457];

    auto g_z_0_z_0_0_y_y_xz = buffer_1010_sppd[458];

    auto g_z_0_z_0_0_y_y_yy = buffer_1010_sppd[459];

    auto g_z_0_z_0_0_y_y_yz = buffer_1010_sppd[460];

    auto g_z_0_z_0_0_y_y_zz = buffer_1010_sppd[461];

    auto g_z_0_z_0_0_y_z_xx = buffer_1010_sppd[462];

    auto g_z_0_z_0_0_y_z_xy = buffer_1010_sppd[463];

    auto g_z_0_z_0_0_y_z_xz = buffer_1010_sppd[464];

    auto g_z_0_z_0_0_y_z_yy = buffer_1010_sppd[465];

    auto g_z_0_z_0_0_y_z_yz = buffer_1010_sppd[466];

    auto g_z_0_z_0_0_y_z_zz = buffer_1010_sppd[467];

    auto g_z_0_z_0_0_z_x_xx = buffer_1010_sppd[468];

    auto g_z_0_z_0_0_z_x_xy = buffer_1010_sppd[469];

    auto g_z_0_z_0_0_z_x_xz = buffer_1010_sppd[470];

    auto g_z_0_z_0_0_z_x_yy = buffer_1010_sppd[471];

    auto g_z_0_z_0_0_z_x_yz = buffer_1010_sppd[472];

    auto g_z_0_z_0_0_z_x_zz = buffer_1010_sppd[473];

    auto g_z_0_z_0_0_z_y_xx = buffer_1010_sppd[474];

    auto g_z_0_z_0_0_z_y_xy = buffer_1010_sppd[475];

    auto g_z_0_z_0_0_z_y_xz = buffer_1010_sppd[476];

    auto g_z_0_z_0_0_z_y_yy = buffer_1010_sppd[477];

    auto g_z_0_z_0_0_z_y_yz = buffer_1010_sppd[478];

    auto g_z_0_z_0_0_z_y_zz = buffer_1010_sppd[479];

    auto g_z_0_z_0_0_z_z_xx = buffer_1010_sppd[480];

    auto g_z_0_z_0_0_z_z_xy = buffer_1010_sppd[481];

    auto g_z_0_z_0_0_z_z_xz = buffer_1010_sppd[482];

    auto g_z_0_z_0_0_z_z_yy = buffer_1010_sppd[483];

    auto g_z_0_z_0_0_z_z_yz = buffer_1010_sppd[484];

    auto g_z_0_z_0_0_z_z_zz = buffer_1010_sppd[485];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_x_0_0_x_x_xx, g_x_0_x_0_0_x_x_xy, g_x_0_x_0_0_x_x_xz, g_x_0_x_0_0_x_x_yy, g_x_0_x_0_0_x_x_yz, g_x_0_x_0_0_x_x_zz, g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_x_x_xx_xx, g_x_x_xx_xy, g_x_x_xx_xz, g_x_x_xx_yy, g_x_x_xx_yz, g_x_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_x_x_xx[i] = -2.0 * g_x_x_0_xx[i] * a_exp + 4.0 * g_x_x_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_x_xy[i] = -2.0 * g_x_x_0_xy[i] * a_exp + 4.0 * g_x_x_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_x_xz[i] = -2.0 * g_x_x_0_xz[i] * a_exp + 4.0 * g_x_x_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_x_yy[i] = -2.0 * g_x_x_0_yy[i] * a_exp + 4.0 * g_x_x_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_x_yz[i] = -2.0 * g_x_x_0_yz[i] * a_exp + 4.0 * g_x_x_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_x_zz[i] = -2.0 * g_x_x_0_zz[i] * a_exp + 4.0 * g_x_x_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_x_0_0_x_y_xx, g_x_0_x_0_0_x_y_xy, g_x_0_x_0_0_x_y_xz, g_x_0_x_0_0_x_y_yy, g_x_0_x_0_0_x_y_yz, g_x_0_x_0_0_x_y_zz, g_x_x_xy_xx, g_x_x_xy_xy, g_x_x_xy_xz, g_x_x_xy_yy, g_x_x_xy_yz, g_x_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_x_y_xx[i] = 4.0 * g_x_x_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_y_xy[i] = 4.0 * g_x_x_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_y_xz[i] = 4.0 * g_x_x_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_y_yy[i] = 4.0 * g_x_x_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_y_yz[i] = 4.0 * g_x_x_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_y_zz[i] = 4.0 * g_x_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_x_0_0_x_z_xx, g_x_0_x_0_0_x_z_xy, g_x_0_x_0_0_x_z_xz, g_x_0_x_0_0_x_z_yy, g_x_0_x_0_0_x_z_yz, g_x_0_x_0_0_x_z_zz, g_x_x_xz_xx, g_x_x_xz_xy, g_x_x_xz_xz, g_x_x_xz_yy, g_x_x_xz_yz, g_x_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_x_z_xx[i] = 4.0 * g_x_x_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_z_xy[i] = 4.0 * g_x_x_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_z_xz[i] = 4.0 * g_x_x_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_z_yy[i] = 4.0 * g_x_x_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_z_yz[i] = 4.0 * g_x_x_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_z_zz[i] = 4.0 * g_x_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_x_0_0_y_x_xx, g_x_0_x_0_0_y_x_xy, g_x_0_x_0_0_y_x_xz, g_x_0_x_0_0_y_x_yy, g_x_0_x_0_0_y_x_yz, g_x_0_x_0_0_y_x_zz, g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_x_y_xx_xx, g_x_y_xx_xy, g_x_y_xx_xz, g_x_y_xx_yy, g_x_y_xx_yz, g_x_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_y_x_xx[i] = -2.0 * g_x_y_0_xx[i] * a_exp + 4.0 * g_x_y_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_x_xy[i] = -2.0 * g_x_y_0_xy[i] * a_exp + 4.0 * g_x_y_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_x_xz[i] = -2.0 * g_x_y_0_xz[i] * a_exp + 4.0 * g_x_y_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_x_yy[i] = -2.0 * g_x_y_0_yy[i] * a_exp + 4.0 * g_x_y_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_x_yz[i] = -2.0 * g_x_y_0_yz[i] * a_exp + 4.0 * g_x_y_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_x_zz[i] = -2.0 * g_x_y_0_zz[i] * a_exp + 4.0 * g_x_y_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_x_0_0_y_y_xx, g_x_0_x_0_0_y_y_xy, g_x_0_x_0_0_y_y_xz, g_x_0_x_0_0_y_y_yy, g_x_0_x_0_0_y_y_yz, g_x_0_x_0_0_y_y_zz, g_x_y_xy_xx, g_x_y_xy_xy, g_x_y_xy_xz, g_x_y_xy_yy, g_x_y_xy_yz, g_x_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_y_y_xx[i] = 4.0 * g_x_y_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_y_xy[i] = 4.0 * g_x_y_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_y_xz[i] = 4.0 * g_x_y_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_y_yy[i] = 4.0 * g_x_y_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_y_yz[i] = 4.0 * g_x_y_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_y_zz[i] = 4.0 * g_x_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_x_0_0_y_z_xx, g_x_0_x_0_0_y_z_xy, g_x_0_x_0_0_y_z_xz, g_x_0_x_0_0_y_z_yy, g_x_0_x_0_0_y_z_yz, g_x_0_x_0_0_y_z_zz, g_x_y_xz_xx, g_x_y_xz_xy, g_x_y_xz_xz, g_x_y_xz_yy, g_x_y_xz_yz, g_x_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_y_z_xx[i] = 4.0 * g_x_y_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_z_xy[i] = 4.0 * g_x_y_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_z_xz[i] = 4.0 * g_x_y_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_z_yy[i] = 4.0 * g_x_y_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_z_yz[i] = 4.0 * g_x_y_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_z_zz[i] = 4.0 * g_x_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_x_0_0_z_x_xx, g_x_0_x_0_0_z_x_xy, g_x_0_x_0_0_z_x_xz, g_x_0_x_0_0_z_x_yy, g_x_0_x_0_0_z_x_yz, g_x_0_x_0_0_z_x_zz, g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_x_z_xx_xx, g_x_z_xx_xy, g_x_z_xx_xz, g_x_z_xx_yy, g_x_z_xx_yz, g_x_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_z_x_xx[i] = -2.0 * g_x_z_0_xx[i] * a_exp + 4.0 * g_x_z_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_x_xy[i] = -2.0 * g_x_z_0_xy[i] * a_exp + 4.0 * g_x_z_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_x_xz[i] = -2.0 * g_x_z_0_xz[i] * a_exp + 4.0 * g_x_z_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_x_yy[i] = -2.0 * g_x_z_0_yy[i] * a_exp + 4.0 * g_x_z_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_x_yz[i] = -2.0 * g_x_z_0_yz[i] * a_exp + 4.0 * g_x_z_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_x_zz[i] = -2.0 * g_x_z_0_zz[i] * a_exp + 4.0 * g_x_z_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_x_0_0_z_y_xx, g_x_0_x_0_0_z_y_xy, g_x_0_x_0_0_z_y_xz, g_x_0_x_0_0_z_y_yy, g_x_0_x_0_0_z_y_yz, g_x_0_x_0_0_z_y_zz, g_x_z_xy_xx, g_x_z_xy_xy, g_x_z_xy_xz, g_x_z_xy_yy, g_x_z_xy_yz, g_x_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_z_y_xx[i] = 4.0 * g_x_z_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_y_xy[i] = 4.0 * g_x_z_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_y_xz[i] = 4.0 * g_x_z_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_y_yy[i] = 4.0 * g_x_z_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_y_yz[i] = 4.0 * g_x_z_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_y_zz[i] = 4.0 * g_x_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_x_0_0_z_z_xx, g_x_0_x_0_0_z_z_xy, g_x_0_x_0_0_z_z_xz, g_x_0_x_0_0_z_z_yy, g_x_0_x_0_0_z_z_yz, g_x_0_x_0_0_z_z_zz, g_x_z_xz_xx, g_x_z_xz_xy, g_x_z_xz_xz, g_x_z_xz_yy, g_x_z_xz_yz, g_x_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_z_z_xx[i] = 4.0 * g_x_z_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_z_xy[i] = 4.0 * g_x_z_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_z_xz[i] = 4.0 * g_x_z_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_z_yy[i] = 4.0 * g_x_z_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_z_yz[i] = 4.0 * g_x_z_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_z_zz[i] = 4.0 * g_x_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_0_y_0_0_x_x_xx, g_x_0_y_0_0_x_x_xy, g_x_0_y_0_0_x_x_xz, g_x_0_y_0_0_x_x_yy, g_x_0_y_0_0_x_x_yz, g_x_0_y_0_0_x_x_zz, g_x_x_xy_xx, g_x_x_xy_xy, g_x_x_xy_xz, g_x_x_xy_yy, g_x_x_xy_yz, g_x_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_x_x_xx[i] = 4.0 * g_x_x_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_x_xy[i] = 4.0 * g_x_x_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_x_xz[i] = 4.0 * g_x_x_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_x_yy[i] = 4.0 * g_x_x_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_x_yz[i] = 4.0 * g_x_x_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_x_zz[i] = 4.0 * g_x_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_0_y_0_0_x_y_xx, g_x_0_y_0_0_x_y_xy, g_x_0_y_0_0_x_y_xz, g_x_0_y_0_0_x_y_yy, g_x_0_y_0_0_x_y_yz, g_x_0_y_0_0_x_y_zz, g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_x_x_yy_xx, g_x_x_yy_xy, g_x_x_yy_xz, g_x_x_yy_yy, g_x_x_yy_yz, g_x_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_x_y_xx[i] = -2.0 * g_x_x_0_xx[i] * a_exp + 4.0 * g_x_x_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_y_xy[i] = -2.0 * g_x_x_0_xy[i] * a_exp + 4.0 * g_x_x_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_y_xz[i] = -2.0 * g_x_x_0_xz[i] * a_exp + 4.0 * g_x_x_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_y_yy[i] = -2.0 * g_x_x_0_yy[i] * a_exp + 4.0 * g_x_x_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_y_yz[i] = -2.0 * g_x_x_0_yz[i] * a_exp + 4.0 * g_x_x_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_y_zz[i] = -2.0 * g_x_x_0_zz[i] * a_exp + 4.0 * g_x_x_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_0_y_0_0_x_z_xx, g_x_0_y_0_0_x_z_xy, g_x_0_y_0_0_x_z_xz, g_x_0_y_0_0_x_z_yy, g_x_0_y_0_0_x_z_yz, g_x_0_y_0_0_x_z_zz, g_x_x_yz_xx, g_x_x_yz_xy, g_x_x_yz_xz, g_x_x_yz_yy, g_x_x_yz_yz, g_x_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_x_z_xx[i] = 4.0 * g_x_x_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_z_xy[i] = 4.0 * g_x_x_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_z_xz[i] = 4.0 * g_x_x_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_z_yy[i] = 4.0 * g_x_x_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_z_yz[i] = 4.0 * g_x_x_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_z_zz[i] = 4.0 * g_x_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_y_0_0_y_x_xx, g_x_0_y_0_0_y_x_xy, g_x_0_y_0_0_y_x_xz, g_x_0_y_0_0_y_x_yy, g_x_0_y_0_0_y_x_yz, g_x_0_y_0_0_y_x_zz, g_x_y_xy_xx, g_x_y_xy_xy, g_x_y_xy_xz, g_x_y_xy_yy, g_x_y_xy_yz, g_x_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_y_x_xx[i] = 4.0 * g_x_y_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_x_xy[i] = 4.0 * g_x_y_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_x_xz[i] = 4.0 * g_x_y_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_x_yy[i] = 4.0 * g_x_y_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_x_yz[i] = 4.0 * g_x_y_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_x_zz[i] = 4.0 * g_x_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_y_0_0_y_y_xx, g_x_0_y_0_0_y_y_xy, g_x_0_y_0_0_y_y_xz, g_x_0_y_0_0_y_y_yy, g_x_0_y_0_0_y_y_yz, g_x_0_y_0_0_y_y_zz, g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_x_y_yy_xx, g_x_y_yy_xy, g_x_y_yy_xz, g_x_y_yy_yy, g_x_y_yy_yz, g_x_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_y_y_xx[i] = -2.0 * g_x_y_0_xx[i] * a_exp + 4.0 * g_x_y_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_y_xy[i] = -2.0 * g_x_y_0_xy[i] * a_exp + 4.0 * g_x_y_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_y_xz[i] = -2.0 * g_x_y_0_xz[i] * a_exp + 4.0 * g_x_y_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_y_yy[i] = -2.0 * g_x_y_0_yy[i] * a_exp + 4.0 * g_x_y_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_y_yz[i] = -2.0 * g_x_y_0_yz[i] * a_exp + 4.0 * g_x_y_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_y_zz[i] = -2.0 * g_x_y_0_zz[i] * a_exp + 4.0 * g_x_y_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_y_0_0_y_z_xx, g_x_0_y_0_0_y_z_xy, g_x_0_y_0_0_y_z_xz, g_x_0_y_0_0_y_z_yy, g_x_0_y_0_0_y_z_yz, g_x_0_y_0_0_y_z_zz, g_x_y_yz_xx, g_x_y_yz_xy, g_x_y_yz_xz, g_x_y_yz_yy, g_x_y_yz_yz, g_x_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_y_z_xx[i] = 4.0 * g_x_y_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_z_xy[i] = 4.0 * g_x_y_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_z_xz[i] = 4.0 * g_x_y_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_z_yy[i] = 4.0 * g_x_y_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_z_yz[i] = 4.0 * g_x_y_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_z_zz[i] = 4.0 * g_x_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_y_0_0_z_x_xx, g_x_0_y_0_0_z_x_xy, g_x_0_y_0_0_z_x_xz, g_x_0_y_0_0_z_x_yy, g_x_0_y_0_0_z_x_yz, g_x_0_y_0_0_z_x_zz, g_x_z_xy_xx, g_x_z_xy_xy, g_x_z_xy_xz, g_x_z_xy_yy, g_x_z_xy_yz, g_x_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_z_x_xx[i] = 4.0 * g_x_z_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_x_xy[i] = 4.0 * g_x_z_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_x_xz[i] = 4.0 * g_x_z_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_x_yy[i] = 4.0 * g_x_z_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_x_yz[i] = 4.0 * g_x_z_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_x_zz[i] = 4.0 * g_x_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_y_0_0_z_y_xx, g_x_0_y_0_0_z_y_xy, g_x_0_y_0_0_z_y_xz, g_x_0_y_0_0_z_y_yy, g_x_0_y_0_0_z_y_yz, g_x_0_y_0_0_z_y_zz, g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_x_z_yy_xx, g_x_z_yy_xy, g_x_z_yy_xz, g_x_z_yy_yy, g_x_z_yy_yz, g_x_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_z_y_xx[i] = -2.0 * g_x_z_0_xx[i] * a_exp + 4.0 * g_x_z_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_y_xy[i] = -2.0 * g_x_z_0_xy[i] * a_exp + 4.0 * g_x_z_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_y_xz[i] = -2.0 * g_x_z_0_xz[i] * a_exp + 4.0 * g_x_z_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_y_yy[i] = -2.0 * g_x_z_0_yy[i] * a_exp + 4.0 * g_x_z_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_y_yz[i] = -2.0 * g_x_z_0_yz[i] * a_exp + 4.0 * g_x_z_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_y_zz[i] = -2.0 * g_x_z_0_zz[i] * a_exp + 4.0 * g_x_z_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_y_0_0_z_z_xx, g_x_0_y_0_0_z_z_xy, g_x_0_y_0_0_z_z_xz, g_x_0_y_0_0_z_z_yy, g_x_0_y_0_0_z_z_yz, g_x_0_y_0_0_z_z_zz, g_x_z_yz_xx, g_x_z_yz_xy, g_x_z_yz_xz, g_x_z_yz_yy, g_x_z_yz_yz, g_x_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_z_z_xx[i] = 4.0 * g_x_z_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_z_xy[i] = 4.0 * g_x_z_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_z_xz[i] = 4.0 * g_x_z_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_z_yy[i] = 4.0 * g_x_z_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_z_yz[i] = 4.0 * g_x_z_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_z_zz[i] = 4.0 * g_x_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_0_z_0_0_x_x_xx, g_x_0_z_0_0_x_x_xy, g_x_0_z_0_0_x_x_xz, g_x_0_z_0_0_x_x_yy, g_x_0_z_0_0_x_x_yz, g_x_0_z_0_0_x_x_zz, g_x_x_xz_xx, g_x_x_xz_xy, g_x_x_xz_xz, g_x_x_xz_yy, g_x_x_xz_yz, g_x_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_x_x_xx[i] = 4.0 * g_x_x_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_x_xy[i] = 4.0 * g_x_x_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_x_xz[i] = 4.0 * g_x_x_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_x_yy[i] = 4.0 * g_x_x_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_x_yz[i] = 4.0 * g_x_x_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_x_zz[i] = 4.0 * g_x_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_0_z_0_0_x_y_xx, g_x_0_z_0_0_x_y_xy, g_x_0_z_0_0_x_y_xz, g_x_0_z_0_0_x_y_yy, g_x_0_z_0_0_x_y_yz, g_x_0_z_0_0_x_y_zz, g_x_x_yz_xx, g_x_x_yz_xy, g_x_x_yz_xz, g_x_x_yz_yy, g_x_x_yz_yz, g_x_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_x_y_xx[i] = 4.0 * g_x_x_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_y_xy[i] = 4.0 * g_x_x_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_y_xz[i] = 4.0 * g_x_x_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_y_yy[i] = 4.0 * g_x_x_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_y_yz[i] = 4.0 * g_x_x_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_y_zz[i] = 4.0 * g_x_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_0_z_0_0_x_z_xx, g_x_0_z_0_0_x_z_xy, g_x_0_z_0_0_x_z_xz, g_x_0_z_0_0_x_z_yy, g_x_0_z_0_0_x_z_yz, g_x_0_z_0_0_x_z_zz, g_x_x_0_xx, g_x_x_0_xy, g_x_x_0_xz, g_x_x_0_yy, g_x_x_0_yz, g_x_x_0_zz, g_x_x_zz_xx, g_x_x_zz_xy, g_x_x_zz_xz, g_x_x_zz_yy, g_x_x_zz_yz, g_x_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_x_z_xx[i] = -2.0 * g_x_x_0_xx[i] * a_exp + 4.0 * g_x_x_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_z_xy[i] = -2.0 * g_x_x_0_xy[i] * a_exp + 4.0 * g_x_x_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_z_xz[i] = -2.0 * g_x_x_0_xz[i] * a_exp + 4.0 * g_x_x_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_z_yy[i] = -2.0 * g_x_x_0_yy[i] * a_exp + 4.0 * g_x_x_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_z_yz[i] = -2.0 * g_x_x_0_yz[i] * a_exp + 4.0 * g_x_x_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_z_zz[i] = -2.0 * g_x_x_0_zz[i] * a_exp + 4.0 * g_x_x_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_0_z_0_0_y_x_xx, g_x_0_z_0_0_y_x_xy, g_x_0_z_0_0_y_x_xz, g_x_0_z_0_0_y_x_yy, g_x_0_z_0_0_y_x_yz, g_x_0_z_0_0_y_x_zz, g_x_y_xz_xx, g_x_y_xz_xy, g_x_y_xz_xz, g_x_y_xz_yy, g_x_y_xz_yz, g_x_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_y_x_xx[i] = 4.0 * g_x_y_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_x_xy[i] = 4.0 * g_x_y_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_x_xz[i] = 4.0 * g_x_y_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_x_yy[i] = 4.0 * g_x_y_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_x_yz[i] = 4.0 * g_x_y_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_x_zz[i] = 4.0 * g_x_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_0_z_0_0_y_y_xx, g_x_0_z_0_0_y_y_xy, g_x_0_z_0_0_y_y_xz, g_x_0_z_0_0_y_y_yy, g_x_0_z_0_0_y_y_yz, g_x_0_z_0_0_y_y_zz, g_x_y_yz_xx, g_x_y_yz_xy, g_x_y_yz_xz, g_x_y_yz_yy, g_x_y_yz_yz, g_x_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_y_y_xx[i] = 4.0 * g_x_y_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_y_xy[i] = 4.0 * g_x_y_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_y_xz[i] = 4.0 * g_x_y_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_y_yy[i] = 4.0 * g_x_y_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_y_yz[i] = 4.0 * g_x_y_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_y_zz[i] = 4.0 * g_x_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_0_z_0_0_y_z_xx, g_x_0_z_0_0_y_z_xy, g_x_0_z_0_0_y_z_xz, g_x_0_z_0_0_y_z_yy, g_x_0_z_0_0_y_z_yz, g_x_0_z_0_0_y_z_zz, g_x_y_0_xx, g_x_y_0_xy, g_x_y_0_xz, g_x_y_0_yy, g_x_y_0_yz, g_x_y_0_zz, g_x_y_zz_xx, g_x_y_zz_xy, g_x_y_zz_xz, g_x_y_zz_yy, g_x_y_zz_yz, g_x_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_y_z_xx[i] = -2.0 * g_x_y_0_xx[i] * a_exp + 4.0 * g_x_y_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_z_xy[i] = -2.0 * g_x_y_0_xy[i] * a_exp + 4.0 * g_x_y_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_z_xz[i] = -2.0 * g_x_y_0_xz[i] * a_exp + 4.0 * g_x_y_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_z_yy[i] = -2.0 * g_x_y_0_yy[i] * a_exp + 4.0 * g_x_y_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_z_yz[i] = -2.0 * g_x_y_0_yz[i] * a_exp + 4.0 * g_x_y_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_z_zz[i] = -2.0 * g_x_y_0_zz[i] * a_exp + 4.0 * g_x_y_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_z_0_0_z_x_xx, g_x_0_z_0_0_z_x_xy, g_x_0_z_0_0_z_x_xz, g_x_0_z_0_0_z_x_yy, g_x_0_z_0_0_z_x_yz, g_x_0_z_0_0_z_x_zz, g_x_z_xz_xx, g_x_z_xz_xy, g_x_z_xz_xz, g_x_z_xz_yy, g_x_z_xz_yz, g_x_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_z_x_xx[i] = 4.0 * g_x_z_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_x_xy[i] = 4.0 * g_x_z_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_x_xz[i] = 4.0 * g_x_z_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_x_yy[i] = 4.0 * g_x_z_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_x_yz[i] = 4.0 * g_x_z_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_x_zz[i] = 4.0 * g_x_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_z_0_0_z_y_xx, g_x_0_z_0_0_z_y_xy, g_x_0_z_0_0_z_y_xz, g_x_0_z_0_0_z_y_yy, g_x_0_z_0_0_z_y_yz, g_x_0_z_0_0_z_y_zz, g_x_z_yz_xx, g_x_z_yz_xy, g_x_z_yz_xz, g_x_z_yz_yy, g_x_z_yz_yz, g_x_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_z_y_xx[i] = 4.0 * g_x_z_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_y_xy[i] = 4.0 * g_x_z_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_y_xz[i] = 4.0 * g_x_z_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_y_yy[i] = 4.0 * g_x_z_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_y_yz[i] = 4.0 * g_x_z_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_y_zz[i] = 4.0 * g_x_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_z_0_0_z_z_xx, g_x_0_z_0_0_z_z_xy, g_x_0_z_0_0_z_z_xz, g_x_0_z_0_0_z_z_yy, g_x_0_z_0_0_z_z_yz, g_x_0_z_0_0_z_z_zz, g_x_z_0_xx, g_x_z_0_xy, g_x_z_0_xz, g_x_z_0_yy, g_x_z_0_yz, g_x_z_0_zz, g_x_z_zz_xx, g_x_z_zz_xy, g_x_z_zz_xz, g_x_z_zz_yy, g_x_z_zz_yz, g_x_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_z_z_xx[i] = -2.0 * g_x_z_0_xx[i] * a_exp + 4.0 * g_x_z_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_z_xy[i] = -2.0 * g_x_z_0_xy[i] * a_exp + 4.0 * g_x_z_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_z_xz[i] = -2.0 * g_x_z_0_xz[i] * a_exp + 4.0 * g_x_z_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_z_yy[i] = -2.0 * g_x_z_0_yy[i] * a_exp + 4.0 * g_x_z_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_z_yz[i] = -2.0 * g_x_z_0_yz[i] * a_exp + 4.0 * g_x_z_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_z_zz[i] = -2.0 * g_x_z_0_zz[i] * a_exp + 4.0 * g_x_z_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_y_0_x_0_0_x_x_xx, g_y_0_x_0_0_x_x_xy, g_y_0_x_0_0_x_x_xz, g_y_0_x_0_0_x_x_yy, g_y_0_x_0_0_x_x_yz, g_y_0_x_0_0_x_x_zz, g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_y_x_xx_xx, g_y_x_xx_xy, g_y_x_xx_xz, g_y_x_xx_yy, g_y_x_xx_yz, g_y_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_x_x_xx[i] = -2.0 * g_y_x_0_xx[i] * a_exp + 4.0 * g_y_x_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_x_xy[i] = -2.0 * g_y_x_0_xy[i] * a_exp + 4.0 * g_y_x_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_x_xz[i] = -2.0 * g_y_x_0_xz[i] * a_exp + 4.0 * g_y_x_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_x_yy[i] = -2.0 * g_y_x_0_yy[i] * a_exp + 4.0 * g_y_x_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_x_yz[i] = -2.0 * g_y_x_0_yz[i] * a_exp + 4.0 * g_y_x_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_x_zz[i] = -2.0 * g_y_x_0_zz[i] * a_exp + 4.0 * g_y_x_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_y_0_x_0_0_x_y_xx, g_y_0_x_0_0_x_y_xy, g_y_0_x_0_0_x_y_xz, g_y_0_x_0_0_x_y_yy, g_y_0_x_0_0_x_y_yz, g_y_0_x_0_0_x_y_zz, g_y_x_xy_xx, g_y_x_xy_xy, g_y_x_xy_xz, g_y_x_xy_yy, g_y_x_xy_yz, g_y_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_x_y_xx[i] = 4.0 * g_y_x_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_y_xy[i] = 4.0 * g_y_x_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_y_xz[i] = 4.0 * g_y_x_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_y_yy[i] = 4.0 * g_y_x_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_y_yz[i] = 4.0 * g_y_x_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_y_zz[i] = 4.0 * g_y_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_y_0_x_0_0_x_z_xx, g_y_0_x_0_0_x_z_xy, g_y_0_x_0_0_x_z_xz, g_y_0_x_0_0_x_z_yy, g_y_0_x_0_0_x_z_yz, g_y_0_x_0_0_x_z_zz, g_y_x_xz_xx, g_y_x_xz_xy, g_y_x_xz_xz, g_y_x_xz_yy, g_y_x_xz_yz, g_y_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_x_z_xx[i] = 4.0 * g_y_x_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_z_xy[i] = 4.0 * g_y_x_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_z_xz[i] = 4.0 * g_y_x_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_z_yy[i] = 4.0 * g_y_x_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_z_yz[i] = 4.0 * g_y_x_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_z_zz[i] = 4.0 * g_y_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_y_0_x_0_0_y_x_xx, g_y_0_x_0_0_y_x_xy, g_y_0_x_0_0_y_x_xz, g_y_0_x_0_0_y_x_yy, g_y_0_x_0_0_y_x_yz, g_y_0_x_0_0_y_x_zz, g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz, g_y_y_xx_xx, g_y_y_xx_xy, g_y_y_xx_xz, g_y_y_xx_yy, g_y_y_xx_yz, g_y_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_y_x_xx[i] = -2.0 * g_y_y_0_xx[i] * a_exp + 4.0 * g_y_y_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_x_xy[i] = -2.0 * g_y_y_0_xy[i] * a_exp + 4.0 * g_y_y_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_x_xz[i] = -2.0 * g_y_y_0_xz[i] * a_exp + 4.0 * g_y_y_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_x_yy[i] = -2.0 * g_y_y_0_yy[i] * a_exp + 4.0 * g_y_y_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_x_yz[i] = -2.0 * g_y_y_0_yz[i] * a_exp + 4.0 * g_y_y_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_x_zz[i] = -2.0 * g_y_y_0_zz[i] * a_exp + 4.0 * g_y_y_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_y_0_x_0_0_y_y_xx, g_y_0_x_0_0_y_y_xy, g_y_0_x_0_0_y_y_xz, g_y_0_x_0_0_y_y_yy, g_y_0_x_0_0_y_y_yz, g_y_0_x_0_0_y_y_zz, g_y_y_xy_xx, g_y_y_xy_xy, g_y_y_xy_xz, g_y_y_xy_yy, g_y_y_xy_yz, g_y_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_y_y_xx[i] = 4.0 * g_y_y_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_y_xy[i] = 4.0 * g_y_y_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_y_xz[i] = 4.0 * g_y_y_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_y_yy[i] = 4.0 * g_y_y_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_y_yz[i] = 4.0 * g_y_y_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_y_zz[i] = 4.0 * g_y_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_y_0_x_0_0_y_z_xx, g_y_0_x_0_0_y_z_xy, g_y_0_x_0_0_y_z_xz, g_y_0_x_0_0_y_z_yy, g_y_0_x_0_0_y_z_yz, g_y_0_x_0_0_y_z_zz, g_y_y_xz_xx, g_y_y_xz_xy, g_y_y_xz_xz, g_y_y_xz_yy, g_y_y_xz_yz, g_y_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_y_z_xx[i] = 4.0 * g_y_y_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_z_xy[i] = 4.0 * g_y_y_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_z_xz[i] = 4.0 * g_y_y_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_z_yy[i] = 4.0 * g_y_y_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_z_yz[i] = 4.0 * g_y_y_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_z_zz[i] = 4.0 * g_y_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_y_0_x_0_0_z_x_xx, g_y_0_x_0_0_z_x_xy, g_y_0_x_0_0_z_x_xz, g_y_0_x_0_0_z_x_yy, g_y_0_x_0_0_z_x_yz, g_y_0_x_0_0_z_x_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz, g_y_z_xx_xx, g_y_z_xx_xy, g_y_z_xx_xz, g_y_z_xx_yy, g_y_z_xx_yz, g_y_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_z_x_xx[i] = -2.0 * g_y_z_0_xx[i] * a_exp + 4.0 * g_y_z_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_x_xy[i] = -2.0 * g_y_z_0_xy[i] * a_exp + 4.0 * g_y_z_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_x_xz[i] = -2.0 * g_y_z_0_xz[i] * a_exp + 4.0 * g_y_z_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_x_yy[i] = -2.0 * g_y_z_0_yy[i] * a_exp + 4.0 * g_y_z_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_x_yz[i] = -2.0 * g_y_z_0_yz[i] * a_exp + 4.0 * g_y_z_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_x_zz[i] = -2.0 * g_y_z_0_zz[i] * a_exp + 4.0 * g_y_z_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_y_0_x_0_0_z_y_xx, g_y_0_x_0_0_z_y_xy, g_y_0_x_0_0_z_y_xz, g_y_0_x_0_0_z_y_yy, g_y_0_x_0_0_z_y_yz, g_y_0_x_0_0_z_y_zz, g_y_z_xy_xx, g_y_z_xy_xy, g_y_z_xy_xz, g_y_z_xy_yy, g_y_z_xy_yz, g_y_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_z_y_xx[i] = 4.0 * g_y_z_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_y_xy[i] = 4.0 * g_y_z_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_y_xz[i] = 4.0 * g_y_z_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_y_yy[i] = 4.0 * g_y_z_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_y_yz[i] = 4.0 * g_y_z_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_y_zz[i] = 4.0 * g_y_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_y_0_x_0_0_z_z_xx, g_y_0_x_0_0_z_z_xy, g_y_0_x_0_0_z_z_xz, g_y_0_x_0_0_z_z_yy, g_y_0_x_0_0_z_z_yz, g_y_0_x_0_0_z_z_zz, g_y_z_xz_xx, g_y_z_xz_xy, g_y_z_xz_xz, g_y_z_xz_yy, g_y_z_xz_yz, g_y_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_z_z_xx[i] = 4.0 * g_y_z_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_z_xy[i] = 4.0 * g_y_z_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_z_xz[i] = 4.0 * g_y_z_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_z_yy[i] = 4.0 * g_y_z_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_z_yz[i] = 4.0 * g_y_z_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_z_zz[i] = 4.0 * g_y_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_y_0_y_0_0_x_x_xx, g_y_0_y_0_0_x_x_xy, g_y_0_y_0_0_x_x_xz, g_y_0_y_0_0_x_x_yy, g_y_0_y_0_0_x_x_yz, g_y_0_y_0_0_x_x_zz, g_y_x_xy_xx, g_y_x_xy_xy, g_y_x_xy_xz, g_y_x_xy_yy, g_y_x_xy_yz, g_y_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_x_x_xx[i] = 4.0 * g_y_x_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_x_xy[i] = 4.0 * g_y_x_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_x_xz[i] = 4.0 * g_y_x_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_x_yy[i] = 4.0 * g_y_x_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_x_yz[i] = 4.0 * g_y_x_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_x_zz[i] = 4.0 * g_y_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_y_0_y_0_0_x_y_xx, g_y_0_y_0_0_x_y_xy, g_y_0_y_0_0_x_y_xz, g_y_0_y_0_0_x_y_yy, g_y_0_y_0_0_x_y_yz, g_y_0_y_0_0_x_y_zz, g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_y_x_yy_xx, g_y_x_yy_xy, g_y_x_yy_xz, g_y_x_yy_yy, g_y_x_yy_yz, g_y_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_x_y_xx[i] = -2.0 * g_y_x_0_xx[i] * a_exp + 4.0 * g_y_x_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_y_xy[i] = -2.0 * g_y_x_0_xy[i] * a_exp + 4.0 * g_y_x_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_y_xz[i] = -2.0 * g_y_x_0_xz[i] * a_exp + 4.0 * g_y_x_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_y_yy[i] = -2.0 * g_y_x_0_yy[i] * a_exp + 4.0 * g_y_x_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_y_yz[i] = -2.0 * g_y_x_0_yz[i] * a_exp + 4.0 * g_y_x_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_y_zz[i] = -2.0 * g_y_x_0_zz[i] * a_exp + 4.0 * g_y_x_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_y_0_y_0_0_x_z_xx, g_y_0_y_0_0_x_z_xy, g_y_0_y_0_0_x_z_xz, g_y_0_y_0_0_x_z_yy, g_y_0_y_0_0_x_z_yz, g_y_0_y_0_0_x_z_zz, g_y_x_yz_xx, g_y_x_yz_xy, g_y_x_yz_xz, g_y_x_yz_yy, g_y_x_yz_yz, g_y_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_x_z_xx[i] = 4.0 * g_y_x_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_z_xy[i] = 4.0 * g_y_x_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_z_xz[i] = 4.0 * g_y_x_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_z_yy[i] = 4.0 * g_y_x_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_z_yz[i] = 4.0 * g_y_x_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_z_zz[i] = 4.0 * g_y_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_y_0_y_0_0_y_x_xx, g_y_0_y_0_0_y_x_xy, g_y_0_y_0_0_y_x_xz, g_y_0_y_0_0_y_x_yy, g_y_0_y_0_0_y_x_yz, g_y_0_y_0_0_y_x_zz, g_y_y_xy_xx, g_y_y_xy_xy, g_y_y_xy_xz, g_y_y_xy_yy, g_y_y_xy_yz, g_y_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_y_x_xx[i] = 4.0 * g_y_y_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_x_xy[i] = 4.0 * g_y_y_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_x_xz[i] = 4.0 * g_y_y_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_x_yy[i] = 4.0 * g_y_y_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_x_yz[i] = 4.0 * g_y_y_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_x_zz[i] = 4.0 * g_y_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_y_0_y_0_0_y_y_xx, g_y_0_y_0_0_y_y_xy, g_y_0_y_0_0_y_y_xz, g_y_0_y_0_0_y_y_yy, g_y_0_y_0_0_y_y_yz, g_y_0_y_0_0_y_y_zz, g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz, g_y_y_yy_xx, g_y_y_yy_xy, g_y_y_yy_xz, g_y_y_yy_yy, g_y_y_yy_yz, g_y_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_y_y_xx[i] = -2.0 * g_y_y_0_xx[i] * a_exp + 4.0 * g_y_y_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_y_xy[i] = -2.0 * g_y_y_0_xy[i] * a_exp + 4.0 * g_y_y_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_y_xz[i] = -2.0 * g_y_y_0_xz[i] * a_exp + 4.0 * g_y_y_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_y_yy[i] = -2.0 * g_y_y_0_yy[i] * a_exp + 4.0 * g_y_y_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_y_yz[i] = -2.0 * g_y_y_0_yz[i] * a_exp + 4.0 * g_y_y_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_y_zz[i] = -2.0 * g_y_y_0_zz[i] * a_exp + 4.0 * g_y_y_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_y_0_y_0_0_y_z_xx, g_y_0_y_0_0_y_z_xy, g_y_0_y_0_0_y_z_xz, g_y_0_y_0_0_y_z_yy, g_y_0_y_0_0_y_z_yz, g_y_0_y_0_0_y_z_zz, g_y_y_yz_xx, g_y_y_yz_xy, g_y_y_yz_xz, g_y_y_yz_yy, g_y_y_yz_yz, g_y_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_y_z_xx[i] = 4.0 * g_y_y_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_z_xy[i] = 4.0 * g_y_y_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_z_xz[i] = 4.0 * g_y_y_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_z_yy[i] = 4.0 * g_y_y_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_z_yz[i] = 4.0 * g_y_y_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_z_zz[i] = 4.0 * g_y_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_y_0_y_0_0_z_x_xx, g_y_0_y_0_0_z_x_xy, g_y_0_y_0_0_z_x_xz, g_y_0_y_0_0_z_x_yy, g_y_0_y_0_0_z_x_yz, g_y_0_y_0_0_z_x_zz, g_y_z_xy_xx, g_y_z_xy_xy, g_y_z_xy_xz, g_y_z_xy_yy, g_y_z_xy_yz, g_y_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_z_x_xx[i] = 4.0 * g_y_z_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_x_xy[i] = 4.0 * g_y_z_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_x_xz[i] = 4.0 * g_y_z_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_x_yy[i] = 4.0 * g_y_z_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_x_yz[i] = 4.0 * g_y_z_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_x_zz[i] = 4.0 * g_y_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_y_0_y_0_0_z_y_xx, g_y_0_y_0_0_z_y_xy, g_y_0_y_0_0_z_y_xz, g_y_0_y_0_0_z_y_yy, g_y_0_y_0_0_z_y_yz, g_y_0_y_0_0_z_y_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz, g_y_z_yy_xx, g_y_z_yy_xy, g_y_z_yy_xz, g_y_z_yy_yy, g_y_z_yy_yz, g_y_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_z_y_xx[i] = -2.0 * g_y_z_0_xx[i] * a_exp + 4.0 * g_y_z_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_y_xy[i] = -2.0 * g_y_z_0_xy[i] * a_exp + 4.0 * g_y_z_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_y_xz[i] = -2.0 * g_y_z_0_xz[i] * a_exp + 4.0 * g_y_z_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_y_yy[i] = -2.0 * g_y_z_0_yy[i] * a_exp + 4.0 * g_y_z_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_y_yz[i] = -2.0 * g_y_z_0_yz[i] * a_exp + 4.0 * g_y_z_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_y_zz[i] = -2.0 * g_y_z_0_zz[i] * a_exp + 4.0 * g_y_z_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_y_0_y_0_0_z_z_xx, g_y_0_y_0_0_z_z_xy, g_y_0_y_0_0_z_z_xz, g_y_0_y_0_0_z_z_yy, g_y_0_y_0_0_z_z_yz, g_y_0_y_0_0_z_z_zz, g_y_z_yz_xx, g_y_z_yz_xy, g_y_z_yz_xz, g_y_z_yz_yy, g_y_z_yz_yz, g_y_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_z_z_xx[i] = 4.0 * g_y_z_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_z_xy[i] = 4.0 * g_y_z_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_z_xz[i] = 4.0 * g_y_z_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_z_yy[i] = 4.0 * g_y_z_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_z_yz[i] = 4.0 * g_y_z_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_z_zz[i] = 4.0 * g_y_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_y_0_z_0_0_x_x_xx, g_y_0_z_0_0_x_x_xy, g_y_0_z_0_0_x_x_xz, g_y_0_z_0_0_x_x_yy, g_y_0_z_0_0_x_x_yz, g_y_0_z_0_0_x_x_zz, g_y_x_xz_xx, g_y_x_xz_xy, g_y_x_xz_xz, g_y_x_xz_yy, g_y_x_xz_yz, g_y_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_x_x_xx[i] = 4.0 * g_y_x_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_x_xy[i] = 4.0 * g_y_x_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_x_xz[i] = 4.0 * g_y_x_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_x_yy[i] = 4.0 * g_y_x_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_x_yz[i] = 4.0 * g_y_x_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_x_zz[i] = 4.0 * g_y_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_y_0_z_0_0_x_y_xx, g_y_0_z_0_0_x_y_xy, g_y_0_z_0_0_x_y_xz, g_y_0_z_0_0_x_y_yy, g_y_0_z_0_0_x_y_yz, g_y_0_z_0_0_x_y_zz, g_y_x_yz_xx, g_y_x_yz_xy, g_y_x_yz_xz, g_y_x_yz_yy, g_y_x_yz_yz, g_y_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_x_y_xx[i] = 4.0 * g_y_x_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_y_xy[i] = 4.0 * g_y_x_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_y_xz[i] = 4.0 * g_y_x_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_y_yy[i] = 4.0 * g_y_x_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_y_yz[i] = 4.0 * g_y_x_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_y_zz[i] = 4.0 * g_y_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_y_0_z_0_0_x_z_xx, g_y_0_z_0_0_x_z_xy, g_y_0_z_0_0_x_z_xz, g_y_0_z_0_0_x_z_yy, g_y_0_z_0_0_x_z_yz, g_y_0_z_0_0_x_z_zz, g_y_x_0_xx, g_y_x_0_xy, g_y_x_0_xz, g_y_x_0_yy, g_y_x_0_yz, g_y_x_0_zz, g_y_x_zz_xx, g_y_x_zz_xy, g_y_x_zz_xz, g_y_x_zz_yy, g_y_x_zz_yz, g_y_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_x_z_xx[i] = -2.0 * g_y_x_0_xx[i] * a_exp + 4.0 * g_y_x_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_z_xy[i] = -2.0 * g_y_x_0_xy[i] * a_exp + 4.0 * g_y_x_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_z_xz[i] = -2.0 * g_y_x_0_xz[i] * a_exp + 4.0 * g_y_x_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_z_yy[i] = -2.0 * g_y_x_0_yy[i] * a_exp + 4.0 * g_y_x_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_z_yz[i] = -2.0 * g_y_x_0_yz[i] * a_exp + 4.0 * g_y_x_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_z_zz[i] = -2.0 * g_y_x_0_zz[i] * a_exp + 4.0 * g_y_x_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_y_0_z_0_0_y_x_xx, g_y_0_z_0_0_y_x_xy, g_y_0_z_0_0_y_x_xz, g_y_0_z_0_0_y_x_yy, g_y_0_z_0_0_y_x_yz, g_y_0_z_0_0_y_x_zz, g_y_y_xz_xx, g_y_y_xz_xy, g_y_y_xz_xz, g_y_y_xz_yy, g_y_y_xz_yz, g_y_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_y_x_xx[i] = 4.0 * g_y_y_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_x_xy[i] = 4.0 * g_y_y_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_x_xz[i] = 4.0 * g_y_y_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_x_yy[i] = 4.0 * g_y_y_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_x_yz[i] = 4.0 * g_y_y_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_x_zz[i] = 4.0 * g_y_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_y_0_z_0_0_y_y_xx, g_y_0_z_0_0_y_y_xy, g_y_0_z_0_0_y_y_xz, g_y_0_z_0_0_y_y_yy, g_y_0_z_0_0_y_y_yz, g_y_0_z_0_0_y_y_zz, g_y_y_yz_xx, g_y_y_yz_xy, g_y_y_yz_xz, g_y_y_yz_yy, g_y_y_yz_yz, g_y_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_y_y_xx[i] = 4.0 * g_y_y_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_y_xy[i] = 4.0 * g_y_y_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_y_xz[i] = 4.0 * g_y_y_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_y_yy[i] = 4.0 * g_y_y_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_y_yz[i] = 4.0 * g_y_y_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_y_zz[i] = 4.0 * g_y_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_y_0_z_0_0_y_z_xx, g_y_0_z_0_0_y_z_xy, g_y_0_z_0_0_y_z_xz, g_y_0_z_0_0_y_z_yy, g_y_0_z_0_0_y_z_yz, g_y_0_z_0_0_y_z_zz, g_y_y_0_xx, g_y_y_0_xy, g_y_y_0_xz, g_y_y_0_yy, g_y_y_0_yz, g_y_y_0_zz, g_y_y_zz_xx, g_y_y_zz_xy, g_y_y_zz_xz, g_y_y_zz_yy, g_y_y_zz_yz, g_y_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_y_z_xx[i] = -2.0 * g_y_y_0_xx[i] * a_exp + 4.0 * g_y_y_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_z_xy[i] = -2.0 * g_y_y_0_xy[i] * a_exp + 4.0 * g_y_y_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_z_xz[i] = -2.0 * g_y_y_0_xz[i] * a_exp + 4.0 * g_y_y_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_z_yy[i] = -2.0 * g_y_y_0_yy[i] * a_exp + 4.0 * g_y_y_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_z_yz[i] = -2.0 * g_y_y_0_yz[i] * a_exp + 4.0 * g_y_y_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_z_zz[i] = -2.0 * g_y_y_0_zz[i] * a_exp + 4.0 * g_y_y_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_y_0_z_0_0_z_x_xx, g_y_0_z_0_0_z_x_xy, g_y_0_z_0_0_z_x_xz, g_y_0_z_0_0_z_x_yy, g_y_0_z_0_0_z_x_yz, g_y_0_z_0_0_z_x_zz, g_y_z_xz_xx, g_y_z_xz_xy, g_y_z_xz_xz, g_y_z_xz_yy, g_y_z_xz_yz, g_y_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_z_x_xx[i] = 4.0 * g_y_z_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_x_xy[i] = 4.0 * g_y_z_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_x_xz[i] = 4.0 * g_y_z_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_x_yy[i] = 4.0 * g_y_z_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_x_yz[i] = 4.0 * g_y_z_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_x_zz[i] = 4.0 * g_y_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_y_0_z_0_0_z_y_xx, g_y_0_z_0_0_z_y_xy, g_y_0_z_0_0_z_y_xz, g_y_0_z_0_0_z_y_yy, g_y_0_z_0_0_z_y_yz, g_y_0_z_0_0_z_y_zz, g_y_z_yz_xx, g_y_z_yz_xy, g_y_z_yz_xz, g_y_z_yz_yy, g_y_z_yz_yz, g_y_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_z_y_xx[i] = 4.0 * g_y_z_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_y_xy[i] = 4.0 * g_y_z_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_y_xz[i] = 4.0 * g_y_z_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_y_yy[i] = 4.0 * g_y_z_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_y_yz[i] = 4.0 * g_y_z_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_y_zz[i] = 4.0 * g_y_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_y_0_z_0_0_z_z_xx, g_y_0_z_0_0_z_z_xy, g_y_0_z_0_0_z_z_xz, g_y_0_z_0_0_z_z_yy, g_y_0_z_0_0_z_z_yz, g_y_0_z_0_0_z_z_zz, g_y_z_0_xx, g_y_z_0_xy, g_y_z_0_xz, g_y_z_0_yy, g_y_z_0_yz, g_y_z_0_zz, g_y_z_zz_xx, g_y_z_zz_xy, g_y_z_zz_xz, g_y_z_zz_yy, g_y_z_zz_yz, g_y_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_z_z_xx[i] = -2.0 * g_y_z_0_xx[i] * a_exp + 4.0 * g_y_z_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_z_xy[i] = -2.0 * g_y_z_0_xy[i] * a_exp + 4.0 * g_y_z_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_z_xz[i] = -2.0 * g_y_z_0_xz[i] * a_exp + 4.0 * g_y_z_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_z_yy[i] = -2.0 * g_y_z_0_yy[i] * a_exp + 4.0 * g_y_z_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_z_yz[i] = -2.0 * g_y_z_0_yz[i] * a_exp + 4.0 * g_y_z_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_z_zz[i] = -2.0 * g_y_z_0_zz[i] * a_exp + 4.0 * g_y_z_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_z_0_x_0_0_x_x_xx, g_z_0_x_0_0_x_x_xy, g_z_0_x_0_0_x_x_xz, g_z_0_x_0_0_x_x_yy, g_z_0_x_0_0_x_x_yz, g_z_0_x_0_0_x_x_zz, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz, g_z_x_xx_xx, g_z_x_xx_xy, g_z_x_xx_xz, g_z_x_xx_yy, g_z_x_xx_yz, g_z_x_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_x_x_xx[i] = -2.0 * g_z_x_0_xx[i] * a_exp + 4.0 * g_z_x_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_x_xy[i] = -2.0 * g_z_x_0_xy[i] * a_exp + 4.0 * g_z_x_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_x_xz[i] = -2.0 * g_z_x_0_xz[i] * a_exp + 4.0 * g_z_x_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_x_yy[i] = -2.0 * g_z_x_0_yy[i] * a_exp + 4.0 * g_z_x_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_x_yz[i] = -2.0 * g_z_x_0_yz[i] * a_exp + 4.0 * g_z_x_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_x_zz[i] = -2.0 * g_z_x_0_zz[i] * a_exp + 4.0 * g_z_x_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_z_0_x_0_0_x_y_xx, g_z_0_x_0_0_x_y_xy, g_z_0_x_0_0_x_y_xz, g_z_0_x_0_0_x_y_yy, g_z_0_x_0_0_x_y_yz, g_z_0_x_0_0_x_y_zz, g_z_x_xy_xx, g_z_x_xy_xy, g_z_x_xy_xz, g_z_x_xy_yy, g_z_x_xy_yz, g_z_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_x_y_xx[i] = 4.0 * g_z_x_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_y_xy[i] = 4.0 * g_z_x_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_y_xz[i] = 4.0 * g_z_x_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_y_yy[i] = 4.0 * g_z_x_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_y_yz[i] = 4.0 * g_z_x_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_y_zz[i] = 4.0 * g_z_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_z_0_x_0_0_x_z_xx, g_z_0_x_0_0_x_z_xy, g_z_0_x_0_0_x_z_xz, g_z_0_x_0_0_x_z_yy, g_z_0_x_0_0_x_z_yz, g_z_0_x_0_0_x_z_zz, g_z_x_xz_xx, g_z_x_xz_xy, g_z_x_xz_xz, g_z_x_xz_yy, g_z_x_xz_yz, g_z_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_x_z_xx[i] = 4.0 * g_z_x_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_z_xy[i] = 4.0 * g_z_x_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_z_xz[i] = 4.0 * g_z_x_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_z_yy[i] = 4.0 * g_z_x_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_z_yz[i] = 4.0 * g_z_x_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_z_zz[i] = 4.0 * g_z_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_z_0_x_0_0_y_x_xx, g_z_0_x_0_0_y_x_xy, g_z_0_x_0_0_y_x_xz, g_z_0_x_0_0_y_x_yy, g_z_0_x_0_0_y_x_yz, g_z_0_x_0_0_y_x_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz, g_z_y_xx_xx, g_z_y_xx_xy, g_z_y_xx_xz, g_z_y_xx_yy, g_z_y_xx_yz, g_z_y_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_y_x_xx[i] = -2.0 * g_z_y_0_xx[i] * a_exp + 4.0 * g_z_y_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_x_xy[i] = -2.0 * g_z_y_0_xy[i] * a_exp + 4.0 * g_z_y_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_x_xz[i] = -2.0 * g_z_y_0_xz[i] * a_exp + 4.0 * g_z_y_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_x_yy[i] = -2.0 * g_z_y_0_yy[i] * a_exp + 4.0 * g_z_y_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_x_yz[i] = -2.0 * g_z_y_0_yz[i] * a_exp + 4.0 * g_z_y_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_x_zz[i] = -2.0 * g_z_y_0_zz[i] * a_exp + 4.0 * g_z_y_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_z_0_x_0_0_y_y_xx, g_z_0_x_0_0_y_y_xy, g_z_0_x_0_0_y_y_xz, g_z_0_x_0_0_y_y_yy, g_z_0_x_0_0_y_y_yz, g_z_0_x_0_0_y_y_zz, g_z_y_xy_xx, g_z_y_xy_xy, g_z_y_xy_xz, g_z_y_xy_yy, g_z_y_xy_yz, g_z_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_y_y_xx[i] = 4.0 * g_z_y_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_y_xy[i] = 4.0 * g_z_y_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_y_xz[i] = 4.0 * g_z_y_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_y_yy[i] = 4.0 * g_z_y_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_y_yz[i] = 4.0 * g_z_y_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_y_zz[i] = 4.0 * g_z_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_z_0_x_0_0_y_z_xx, g_z_0_x_0_0_y_z_xy, g_z_0_x_0_0_y_z_xz, g_z_0_x_0_0_y_z_yy, g_z_0_x_0_0_y_z_yz, g_z_0_x_0_0_y_z_zz, g_z_y_xz_xx, g_z_y_xz_xy, g_z_y_xz_xz, g_z_y_xz_yy, g_z_y_xz_yz, g_z_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_y_z_xx[i] = 4.0 * g_z_y_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_z_xy[i] = 4.0 * g_z_y_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_z_xz[i] = 4.0 * g_z_y_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_z_yy[i] = 4.0 * g_z_y_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_z_yz[i] = 4.0 * g_z_y_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_z_zz[i] = 4.0 * g_z_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_z_0_x_0_0_z_x_xx, g_z_0_x_0_0_z_x_xy, g_z_0_x_0_0_z_x_xz, g_z_0_x_0_0_z_x_yy, g_z_0_x_0_0_z_x_yz, g_z_0_x_0_0_z_x_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz, g_z_z_xx_xx, g_z_z_xx_xy, g_z_z_xx_xz, g_z_z_xx_yy, g_z_z_xx_yz, g_z_z_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_z_x_xx[i] = -2.0 * g_z_z_0_xx[i] * a_exp + 4.0 * g_z_z_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_x_xy[i] = -2.0 * g_z_z_0_xy[i] * a_exp + 4.0 * g_z_z_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_x_xz[i] = -2.0 * g_z_z_0_xz[i] * a_exp + 4.0 * g_z_z_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_x_yy[i] = -2.0 * g_z_z_0_yy[i] * a_exp + 4.0 * g_z_z_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_x_yz[i] = -2.0 * g_z_z_0_yz[i] * a_exp + 4.0 * g_z_z_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_x_zz[i] = -2.0 * g_z_z_0_zz[i] * a_exp + 4.0 * g_z_z_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_z_0_x_0_0_z_y_xx, g_z_0_x_0_0_z_y_xy, g_z_0_x_0_0_z_y_xz, g_z_0_x_0_0_z_y_yy, g_z_0_x_0_0_z_y_yz, g_z_0_x_0_0_z_y_zz, g_z_z_xy_xx, g_z_z_xy_xy, g_z_z_xy_xz, g_z_z_xy_yy, g_z_z_xy_yz, g_z_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_z_y_xx[i] = 4.0 * g_z_z_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_y_xy[i] = 4.0 * g_z_z_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_y_xz[i] = 4.0 * g_z_z_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_y_yy[i] = 4.0 * g_z_z_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_y_yz[i] = 4.0 * g_z_z_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_y_zz[i] = 4.0 * g_z_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_z_0_x_0_0_z_z_xx, g_z_0_x_0_0_z_z_xy, g_z_0_x_0_0_z_z_xz, g_z_0_x_0_0_z_z_yy, g_z_0_x_0_0_z_z_yz, g_z_0_x_0_0_z_z_zz, g_z_z_xz_xx, g_z_z_xz_xy, g_z_z_xz_xz, g_z_z_xz_yy, g_z_z_xz_yz, g_z_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_z_z_xx[i] = 4.0 * g_z_z_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_z_xy[i] = 4.0 * g_z_z_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_z_xz[i] = 4.0 * g_z_z_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_z_yy[i] = 4.0 * g_z_z_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_z_yz[i] = 4.0 * g_z_z_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_z_zz[i] = 4.0 * g_z_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_z_0_y_0_0_x_x_xx, g_z_0_y_0_0_x_x_xy, g_z_0_y_0_0_x_x_xz, g_z_0_y_0_0_x_x_yy, g_z_0_y_0_0_x_x_yz, g_z_0_y_0_0_x_x_zz, g_z_x_xy_xx, g_z_x_xy_xy, g_z_x_xy_xz, g_z_x_xy_yy, g_z_x_xy_yz, g_z_x_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_x_x_xx[i] = 4.0 * g_z_x_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_x_xy[i] = 4.0 * g_z_x_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_x_xz[i] = 4.0 * g_z_x_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_x_yy[i] = 4.0 * g_z_x_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_x_yz[i] = 4.0 * g_z_x_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_x_zz[i] = 4.0 * g_z_x_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_z_0_y_0_0_x_y_xx, g_z_0_y_0_0_x_y_xy, g_z_0_y_0_0_x_y_xz, g_z_0_y_0_0_x_y_yy, g_z_0_y_0_0_x_y_yz, g_z_0_y_0_0_x_y_zz, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz, g_z_x_yy_xx, g_z_x_yy_xy, g_z_x_yy_xz, g_z_x_yy_yy, g_z_x_yy_yz, g_z_x_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_x_y_xx[i] = -2.0 * g_z_x_0_xx[i] * a_exp + 4.0 * g_z_x_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_y_xy[i] = -2.0 * g_z_x_0_xy[i] * a_exp + 4.0 * g_z_x_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_y_xz[i] = -2.0 * g_z_x_0_xz[i] * a_exp + 4.0 * g_z_x_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_y_yy[i] = -2.0 * g_z_x_0_yy[i] * a_exp + 4.0 * g_z_x_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_y_yz[i] = -2.0 * g_z_x_0_yz[i] * a_exp + 4.0 * g_z_x_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_y_zz[i] = -2.0 * g_z_x_0_zz[i] * a_exp + 4.0 * g_z_x_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_z_0_y_0_0_x_z_xx, g_z_0_y_0_0_x_z_xy, g_z_0_y_0_0_x_z_xz, g_z_0_y_0_0_x_z_yy, g_z_0_y_0_0_x_z_yz, g_z_0_y_0_0_x_z_zz, g_z_x_yz_xx, g_z_x_yz_xy, g_z_x_yz_xz, g_z_x_yz_yy, g_z_x_yz_yz, g_z_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_x_z_xx[i] = 4.0 * g_z_x_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_z_xy[i] = 4.0 * g_z_x_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_z_xz[i] = 4.0 * g_z_x_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_z_yy[i] = 4.0 * g_z_x_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_z_yz[i] = 4.0 * g_z_x_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_z_zz[i] = 4.0 * g_z_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_z_0_y_0_0_y_x_xx, g_z_0_y_0_0_y_x_xy, g_z_0_y_0_0_y_x_xz, g_z_0_y_0_0_y_x_yy, g_z_0_y_0_0_y_x_yz, g_z_0_y_0_0_y_x_zz, g_z_y_xy_xx, g_z_y_xy_xy, g_z_y_xy_xz, g_z_y_xy_yy, g_z_y_xy_yz, g_z_y_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_y_x_xx[i] = 4.0 * g_z_y_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_x_xy[i] = 4.0 * g_z_y_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_x_xz[i] = 4.0 * g_z_y_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_x_yy[i] = 4.0 * g_z_y_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_x_yz[i] = 4.0 * g_z_y_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_x_zz[i] = 4.0 * g_z_y_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_z_0_y_0_0_y_y_xx, g_z_0_y_0_0_y_y_xy, g_z_0_y_0_0_y_y_xz, g_z_0_y_0_0_y_y_yy, g_z_0_y_0_0_y_y_yz, g_z_0_y_0_0_y_y_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz, g_z_y_yy_xx, g_z_y_yy_xy, g_z_y_yy_xz, g_z_y_yy_yy, g_z_y_yy_yz, g_z_y_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_y_y_xx[i] = -2.0 * g_z_y_0_xx[i] * a_exp + 4.0 * g_z_y_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_y_xy[i] = -2.0 * g_z_y_0_xy[i] * a_exp + 4.0 * g_z_y_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_y_xz[i] = -2.0 * g_z_y_0_xz[i] * a_exp + 4.0 * g_z_y_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_y_yy[i] = -2.0 * g_z_y_0_yy[i] * a_exp + 4.0 * g_z_y_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_y_yz[i] = -2.0 * g_z_y_0_yz[i] * a_exp + 4.0 * g_z_y_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_y_zz[i] = -2.0 * g_z_y_0_zz[i] * a_exp + 4.0 * g_z_y_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_z_0_y_0_0_y_z_xx, g_z_0_y_0_0_y_z_xy, g_z_0_y_0_0_y_z_xz, g_z_0_y_0_0_y_z_yy, g_z_0_y_0_0_y_z_yz, g_z_0_y_0_0_y_z_zz, g_z_y_yz_xx, g_z_y_yz_xy, g_z_y_yz_xz, g_z_y_yz_yy, g_z_y_yz_yz, g_z_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_y_z_xx[i] = 4.0 * g_z_y_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_z_xy[i] = 4.0 * g_z_y_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_z_xz[i] = 4.0 * g_z_y_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_z_yy[i] = 4.0 * g_z_y_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_z_yz[i] = 4.0 * g_z_y_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_z_zz[i] = 4.0 * g_z_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_z_0_y_0_0_z_x_xx, g_z_0_y_0_0_z_x_xy, g_z_0_y_0_0_z_x_xz, g_z_0_y_0_0_z_x_yy, g_z_0_y_0_0_z_x_yz, g_z_0_y_0_0_z_x_zz, g_z_z_xy_xx, g_z_z_xy_xy, g_z_z_xy_xz, g_z_z_xy_yy, g_z_z_xy_yz, g_z_z_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_z_x_xx[i] = 4.0 * g_z_z_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_x_xy[i] = 4.0 * g_z_z_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_x_xz[i] = 4.0 * g_z_z_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_x_yy[i] = 4.0 * g_z_z_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_x_yz[i] = 4.0 * g_z_z_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_x_zz[i] = 4.0 * g_z_z_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_z_0_y_0_0_z_y_xx, g_z_0_y_0_0_z_y_xy, g_z_0_y_0_0_z_y_xz, g_z_0_y_0_0_z_y_yy, g_z_0_y_0_0_z_y_yz, g_z_0_y_0_0_z_y_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz, g_z_z_yy_xx, g_z_z_yy_xy, g_z_z_yy_xz, g_z_z_yy_yy, g_z_z_yy_yz, g_z_z_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_z_y_xx[i] = -2.0 * g_z_z_0_xx[i] * a_exp + 4.0 * g_z_z_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_y_xy[i] = -2.0 * g_z_z_0_xy[i] * a_exp + 4.0 * g_z_z_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_y_xz[i] = -2.0 * g_z_z_0_xz[i] * a_exp + 4.0 * g_z_z_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_y_yy[i] = -2.0 * g_z_z_0_yy[i] * a_exp + 4.0 * g_z_z_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_y_yz[i] = -2.0 * g_z_z_0_yz[i] * a_exp + 4.0 * g_z_z_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_y_zz[i] = -2.0 * g_z_z_0_zz[i] * a_exp + 4.0 * g_z_z_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_z_0_y_0_0_z_z_xx, g_z_0_y_0_0_z_z_xy, g_z_0_y_0_0_z_z_xz, g_z_0_y_0_0_z_z_yy, g_z_0_y_0_0_z_z_yz, g_z_0_y_0_0_z_z_zz, g_z_z_yz_xx, g_z_z_yz_xy, g_z_z_yz_xz, g_z_z_yz_yy, g_z_z_yz_yz, g_z_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_z_z_xx[i] = 4.0 * g_z_z_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_z_xy[i] = 4.0 * g_z_z_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_z_xz[i] = 4.0 * g_z_z_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_z_yy[i] = 4.0 * g_z_z_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_z_yz[i] = 4.0 * g_z_z_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_z_zz[i] = 4.0 * g_z_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_z_0_z_0_0_x_x_xx, g_z_0_z_0_0_x_x_xy, g_z_0_z_0_0_x_x_xz, g_z_0_z_0_0_x_x_yy, g_z_0_z_0_0_x_x_yz, g_z_0_z_0_0_x_x_zz, g_z_x_xz_xx, g_z_x_xz_xy, g_z_x_xz_xz, g_z_x_xz_yy, g_z_x_xz_yz, g_z_x_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_x_x_xx[i] = 4.0 * g_z_x_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_x_xy[i] = 4.0 * g_z_x_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_x_xz[i] = 4.0 * g_z_x_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_x_yy[i] = 4.0 * g_z_x_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_x_yz[i] = 4.0 * g_z_x_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_x_zz[i] = 4.0 * g_z_x_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_z_0_z_0_0_x_y_xx, g_z_0_z_0_0_x_y_xy, g_z_0_z_0_0_x_y_xz, g_z_0_z_0_0_x_y_yy, g_z_0_z_0_0_x_y_yz, g_z_0_z_0_0_x_y_zz, g_z_x_yz_xx, g_z_x_yz_xy, g_z_x_yz_xz, g_z_x_yz_yy, g_z_x_yz_yz, g_z_x_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_x_y_xx[i] = 4.0 * g_z_x_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_y_xy[i] = 4.0 * g_z_x_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_y_xz[i] = 4.0 * g_z_x_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_y_yy[i] = 4.0 * g_z_x_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_y_yz[i] = 4.0 * g_z_x_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_y_zz[i] = 4.0 * g_z_x_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_z_0_z_0_0_x_z_xx, g_z_0_z_0_0_x_z_xy, g_z_0_z_0_0_x_z_xz, g_z_0_z_0_0_x_z_yy, g_z_0_z_0_0_x_z_yz, g_z_0_z_0_0_x_z_zz, g_z_x_0_xx, g_z_x_0_xy, g_z_x_0_xz, g_z_x_0_yy, g_z_x_0_yz, g_z_x_0_zz, g_z_x_zz_xx, g_z_x_zz_xy, g_z_x_zz_xz, g_z_x_zz_yy, g_z_x_zz_yz, g_z_x_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_x_z_xx[i] = -2.0 * g_z_x_0_xx[i] * a_exp + 4.0 * g_z_x_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_z_xy[i] = -2.0 * g_z_x_0_xy[i] * a_exp + 4.0 * g_z_x_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_z_xz[i] = -2.0 * g_z_x_0_xz[i] * a_exp + 4.0 * g_z_x_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_z_yy[i] = -2.0 * g_z_x_0_yy[i] * a_exp + 4.0 * g_z_x_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_z_yz[i] = -2.0 * g_z_x_0_yz[i] * a_exp + 4.0 * g_z_x_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_z_zz[i] = -2.0 * g_z_x_0_zz[i] * a_exp + 4.0 * g_z_x_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_z_0_z_0_0_y_x_xx, g_z_0_z_0_0_y_x_xy, g_z_0_z_0_0_y_x_xz, g_z_0_z_0_0_y_x_yy, g_z_0_z_0_0_y_x_yz, g_z_0_z_0_0_y_x_zz, g_z_y_xz_xx, g_z_y_xz_xy, g_z_y_xz_xz, g_z_y_xz_yy, g_z_y_xz_yz, g_z_y_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_y_x_xx[i] = 4.0 * g_z_y_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_x_xy[i] = 4.0 * g_z_y_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_x_xz[i] = 4.0 * g_z_y_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_x_yy[i] = 4.0 * g_z_y_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_x_yz[i] = 4.0 * g_z_y_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_x_zz[i] = 4.0 * g_z_y_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_z_0_z_0_0_y_y_xx, g_z_0_z_0_0_y_y_xy, g_z_0_z_0_0_y_y_xz, g_z_0_z_0_0_y_y_yy, g_z_0_z_0_0_y_y_yz, g_z_0_z_0_0_y_y_zz, g_z_y_yz_xx, g_z_y_yz_xy, g_z_y_yz_xz, g_z_y_yz_yy, g_z_y_yz_yz, g_z_y_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_y_y_xx[i] = 4.0 * g_z_y_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_y_xy[i] = 4.0 * g_z_y_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_y_xz[i] = 4.0 * g_z_y_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_y_yy[i] = 4.0 * g_z_y_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_y_yz[i] = 4.0 * g_z_y_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_y_zz[i] = 4.0 * g_z_y_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_z_0_z_0_0_y_z_xx, g_z_0_z_0_0_y_z_xy, g_z_0_z_0_0_y_z_xz, g_z_0_z_0_0_y_z_yy, g_z_0_z_0_0_y_z_yz, g_z_0_z_0_0_y_z_zz, g_z_y_0_xx, g_z_y_0_xy, g_z_y_0_xz, g_z_y_0_yy, g_z_y_0_yz, g_z_y_0_zz, g_z_y_zz_xx, g_z_y_zz_xy, g_z_y_zz_xz, g_z_y_zz_yy, g_z_y_zz_yz, g_z_y_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_y_z_xx[i] = -2.0 * g_z_y_0_xx[i] * a_exp + 4.0 * g_z_y_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_z_xy[i] = -2.0 * g_z_y_0_xy[i] * a_exp + 4.0 * g_z_y_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_z_xz[i] = -2.0 * g_z_y_0_xz[i] * a_exp + 4.0 * g_z_y_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_z_yy[i] = -2.0 * g_z_y_0_yy[i] * a_exp + 4.0 * g_z_y_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_z_yz[i] = -2.0 * g_z_y_0_yz[i] * a_exp + 4.0 * g_z_y_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_z_zz[i] = -2.0 * g_z_y_0_zz[i] * a_exp + 4.0 * g_z_y_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_z_0_z_0_0_z_x_xx, g_z_0_z_0_0_z_x_xy, g_z_0_z_0_0_z_x_xz, g_z_0_z_0_0_z_x_yy, g_z_0_z_0_0_z_x_yz, g_z_0_z_0_0_z_x_zz, g_z_z_xz_xx, g_z_z_xz_xy, g_z_z_xz_xz, g_z_z_xz_yy, g_z_z_xz_yz, g_z_z_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_z_x_xx[i] = 4.0 * g_z_z_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_x_xy[i] = 4.0 * g_z_z_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_x_xz[i] = 4.0 * g_z_z_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_x_yy[i] = 4.0 * g_z_z_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_x_yz[i] = 4.0 * g_z_z_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_x_zz[i] = 4.0 * g_z_z_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_z_0_z_0_0_z_y_xx, g_z_0_z_0_0_z_y_xy, g_z_0_z_0_0_z_y_xz, g_z_0_z_0_0_z_y_yy, g_z_0_z_0_0_z_y_yz, g_z_0_z_0_0_z_y_zz, g_z_z_yz_xx, g_z_z_yz_xy, g_z_z_yz_xz, g_z_z_yz_yy, g_z_z_yz_yz, g_z_z_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_z_y_xx[i] = 4.0 * g_z_z_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_y_xy[i] = 4.0 * g_z_z_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_y_xz[i] = 4.0 * g_z_z_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_y_yy[i] = 4.0 * g_z_z_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_y_yz[i] = 4.0 * g_z_z_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_y_zz[i] = 4.0 * g_z_z_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_z_0_z_0_0_z_z_xx, g_z_0_z_0_0_z_z_xy, g_z_0_z_0_0_z_z_xz, g_z_0_z_0_0_z_z_yy, g_z_0_z_0_0_z_z_yz, g_z_0_z_0_0_z_z_zz, g_z_z_0_xx, g_z_z_0_xy, g_z_z_0_xz, g_z_z_0_yy, g_z_z_0_yz, g_z_z_0_zz, g_z_z_zz_xx, g_z_z_zz_xy, g_z_z_zz_xz, g_z_z_zz_yy, g_z_z_zz_yz, g_z_z_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_z_z_xx[i] = -2.0 * g_z_z_0_xx[i] * a_exp + 4.0 * g_z_z_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_z_xy[i] = -2.0 * g_z_z_0_xy[i] * a_exp + 4.0 * g_z_z_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_z_xz[i] = -2.0 * g_z_z_0_xz[i] * a_exp + 4.0 * g_z_z_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_z_yy[i] = -2.0 * g_z_z_0_yy[i] * a_exp + 4.0 * g_z_z_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_z_yz[i] = -2.0 * g_z_z_0_yz[i] * a_exp + 4.0 * g_z_z_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_z_zz[i] = -2.0 * g_z_z_0_zz[i] * a_exp + 4.0 * g_z_z_zz_zz[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace

