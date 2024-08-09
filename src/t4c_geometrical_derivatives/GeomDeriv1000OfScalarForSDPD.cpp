#include "GeomDeriv1000OfScalarForSDPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_sdpd_0(CSimdArray<double>& buffer_1000_sdpd,
                     const CSimdArray<double>& buffer_pdpd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_sdpd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_sdpd

    auto g_x_0_0_0_0_xx_x_xx = buffer_1000_sdpd[0];

    auto g_x_0_0_0_0_xx_x_xy = buffer_1000_sdpd[1];

    auto g_x_0_0_0_0_xx_x_xz = buffer_1000_sdpd[2];

    auto g_x_0_0_0_0_xx_x_yy = buffer_1000_sdpd[3];

    auto g_x_0_0_0_0_xx_x_yz = buffer_1000_sdpd[4];

    auto g_x_0_0_0_0_xx_x_zz = buffer_1000_sdpd[5];

    auto g_x_0_0_0_0_xx_y_xx = buffer_1000_sdpd[6];

    auto g_x_0_0_0_0_xx_y_xy = buffer_1000_sdpd[7];

    auto g_x_0_0_0_0_xx_y_xz = buffer_1000_sdpd[8];

    auto g_x_0_0_0_0_xx_y_yy = buffer_1000_sdpd[9];

    auto g_x_0_0_0_0_xx_y_yz = buffer_1000_sdpd[10];

    auto g_x_0_0_0_0_xx_y_zz = buffer_1000_sdpd[11];

    auto g_x_0_0_0_0_xx_z_xx = buffer_1000_sdpd[12];

    auto g_x_0_0_0_0_xx_z_xy = buffer_1000_sdpd[13];

    auto g_x_0_0_0_0_xx_z_xz = buffer_1000_sdpd[14];

    auto g_x_0_0_0_0_xx_z_yy = buffer_1000_sdpd[15];

    auto g_x_0_0_0_0_xx_z_yz = buffer_1000_sdpd[16];

    auto g_x_0_0_0_0_xx_z_zz = buffer_1000_sdpd[17];

    auto g_x_0_0_0_0_xy_x_xx = buffer_1000_sdpd[18];

    auto g_x_0_0_0_0_xy_x_xy = buffer_1000_sdpd[19];

    auto g_x_0_0_0_0_xy_x_xz = buffer_1000_sdpd[20];

    auto g_x_0_0_0_0_xy_x_yy = buffer_1000_sdpd[21];

    auto g_x_0_0_0_0_xy_x_yz = buffer_1000_sdpd[22];

    auto g_x_0_0_0_0_xy_x_zz = buffer_1000_sdpd[23];

    auto g_x_0_0_0_0_xy_y_xx = buffer_1000_sdpd[24];

    auto g_x_0_0_0_0_xy_y_xy = buffer_1000_sdpd[25];

    auto g_x_0_0_0_0_xy_y_xz = buffer_1000_sdpd[26];

    auto g_x_0_0_0_0_xy_y_yy = buffer_1000_sdpd[27];

    auto g_x_0_0_0_0_xy_y_yz = buffer_1000_sdpd[28];

    auto g_x_0_0_0_0_xy_y_zz = buffer_1000_sdpd[29];

    auto g_x_0_0_0_0_xy_z_xx = buffer_1000_sdpd[30];

    auto g_x_0_0_0_0_xy_z_xy = buffer_1000_sdpd[31];

    auto g_x_0_0_0_0_xy_z_xz = buffer_1000_sdpd[32];

    auto g_x_0_0_0_0_xy_z_yy = buffer_1000_sdpd[33];

    auto g_x_0_0_0_0_xy_z_yz = buffer_1000_sdpd[34];

    auto g_x_0_0_0_0_xy_z_zz = buffer_1000_sdpd[35];

    auto g_x_0_0_0_0_xz_x_xx = buffer_1000_sdpd[36];

    auto g_x_0_0_0_0_xz_x_xy = buffer_1000_sdpd[37];

    auto g_x_0_0_0_0_xz_x_xz = buffer_1000_sdpd[38];

    auto g_x_0_0_0_0_xz_x_yy = buffer_1000_sdpd[39];

    auto g_x_0_0_0_0_xz_x_yz = buffer_1000_sdpd[40];

    auto g_x_0_0_0_0_xz_x_zz = buffer_1000_sdpd[41];

    auto g_x_0_0_0_0_xz_y_xx = buffer_1000_sdpd[42];

    auto g_x_0_0_0_0_xz_y_xy = buffer_1000_sdpd[43];

    auto g_x_0_0_0_0_xz_y_xz = buffer_1000_sdpd[44];

    auto g_x_0_0_0_0_xz_y_yy = buffer_1000_sdpd[45];

    auto g_x_0_0_0_0_xz_y_yz = buffer_1000_sdpd[46];

    auto g_x_0_0_0_0_xz_y_zz = buffer_1000_sdpd[47];

    auto g_x_0_0_0_0_xz_z_xx = buffer_1000_sdpd[48];

    auto g_x_0_0_0_0_xz_z_xy = buffer_1000_sdpd[49];

    auto g_x_0_0_0_0_xz_z_xz = buffer_1000_sdpd[50];

    auto g_x_0_0_0_0_xz_z_yy = buffer_1000_sdpd[51];

    auto g_x_0_0_0_0_xz_z_yz = buffer_1000_sdpd[52];

    auto g_x_0_0_0_0_xz_z_zz = buffer_1000_sdpd[53];

    auto g_x_0_0_0_0_yy_x_xx = buffer_1000_sdpd[54];

    auto g_x_0_0_0_0_yy_x_xy = buffer_1000_sdpd[55];

    auto g_x_0_0_0_0_yy_x_xz = buffer_1000_sdpd[56];

    auto g_x_0_0_0_0_yy_x_yy = buffer_1000_sdpd[57];

    auto g_x_0_0_0_0_yy_x_yz = buffer_1000_sdpd[58];

    auto g_x_0_0_0_0_yy_x_zz = buffer_1000_sdpd[59];

    auto g_x_0_0_0_0_yy_y_xx = buffer_1000_sdpd[60];

    auto g_x_0_0_0_0_yy_y_xy = buffer_1000_sdpd[61];

    auto g_x_0_0_0_0_yy_y_xz = buffer_1000_sdpd[62];

    auto g_x_0_0_0_0_yy_y_yy = buffer_1000_sdpd[63];

    auto g_x_0_0_0_0_yy_y_yz = buffer_1000_sdpd[64];

    auto g_x_0_0_0_0_yy_y_zz = buffer_1000_sdpd[65];

    auto g_x_0_0_0_0_yy_z_xx = buffer_1000_sdpd[66];

    auto g_x_0_0_0_0_yy_z_xy = buffer_1000_sdpd[67];

    auto g_x_0_0_0_0_yy_z_xz = buffer_1000_sdpd[68];

    auto g_x_0_0_0_0_yy_z_yy = buffer_1000_sdpd[69];

    auto g_x_0_0_0_0_yy_z_yz = buffer_1000_sdpd[70];

    auto g_x_0_0_0_0_yy_z_zz = buffer_1000_sdpd[71];

    auto g_x_0_0_0_0_yz_x_xx = buffer_1000_sdpd[72];

    auto g_x_0_0_0_0_yz_x_xy = buffer_1000_sdpd[73];

    auto g_x_0_0_0_0_yz_x_xz = buffer_1000_sdpd[74];

    auto g_x_0_0_0_0_yz_x_yy = buffer_1000_sdpd[75];

    auto g_x_0_0_0_0_yz_x_yz = buffer_1000_sdpd[76];

    auto g_x_0_0_0_0_yz_x_zz = buffer_1000_sdpd[77];

    auto g_x_0_0_0_0_yz_y_xx = buffer_1000_sdpd[78];

    auto g_x_0_0_0_0_yz_y_xy = buffer_1000_sdpd[79];

    auto g_x_0_0_0_0_yz_y_xz = buffer_1000_sdpd[80];

    auto g_x_0_0_0_0_yz_y_yy = buffer_1000_sdpd[81];

    auto g_x_0_0_0_0_yz_y_yz = buffer_1000_sdpd[82];

    auto g_x_0_0_0_0_yz_y_zz = buffer_1000_sdpd[83];

    auto g_x_0_0_0_0_yz_z_xx = buffer_1000_sdpd[84];

    auto g_x_0_0_0_0_yz_z_xy = buffer_1000_sdpd[85];

    auto g_x_0_0_0_0_yz_z_xz = buffer_1000_sdpd[86];

    auto g_x_0_0_0_0_yz_z_yy = buffer_1000_sdpd[87];

    auto g_x_0_0_0_0_yz_z_yz = buffer_1000_sdpd[88];

    auto g_x_0_0_0_0_yz_z_zz = buffer_1000_sdpd[89];

    auto g_x_0_0_0_0_zz_x_xx = buffer_1000_sdpd[90];

    auto g_x_0_0_0_0_zz_x_xy = buffer_1000_sdpd[91];

    auto g_x_0_0_0_0_zz_x_xz = buffer_1000_sdpd[92];

    auto g_x_0_0_0_0_zz_x_yy = buffer_1000_sdpd[93];

    auto g_x_0_0_0_0_zz_x_yz = buffer_1000_sdpd[94];

    auto g_x_0_0_0_0_zz_x_zz = buffer_1000_sdpd[95];

    auto g_x_0_0_0_0_zz_y_xx = buffer_1000_sdpd[96];

    auto g_x_0_0_0_0_zz_y_xy = buffer_1000_sdpd[97];

    auto g_x_0_0_0_0_zz_y_xz = buffer_1000_sdpd[98];

    auto g_x_0_0_0_0_zz_y_yy = buffer_1000_sdpd[99];

    auto g_x_0_0_0_0_zz_y_yz = buffer_1000_sdpd[100];

    auto g_x_0_0_0_0_zz_y_zz = buffer_1000_sdpd[101];

    auto g_x_0_0_0_0_zz_z_xx = buffer_1000_sdpd[102];

    auto g_x_0_0_0_0_zz_z_xy = buffer_1000_sdpd[103];

    auto g_x_0_0_0_0_zz_z_xz = buffer_1000_sdpd[104];

    auto g_x_0_0_0_0_zz_z_yy = buffer_1000_sdpd[105];

    auto g_x_0_0_0_0_zz_z_yz = buffer_1000_sdpd[106];

    auto g_x_0_0_0_0_zz_z_zz = buffer_1000_sdpd[107];

    auto g_y_0_0_0_0_xx_x_xx = buffer_1000_sdpd[108];

    auto g_y_0_0_0_0_xx_x_xy = buffer_1000_sdpd[109];

    auto g_y_0_0_0_0_xx_x_xz = buffer_1000_sdpd[110];

    auto g_y_0_0_0_0_xx_x_yy = buffer_1000_sdpd[111];

    auto g_y_0_0_0_0_xx_x_yz = buffer_1000_sdpd[112];

    auto g_y_0_0_0_0_xx_x_zz = buffer_1000_sdpd[113];

    auto g_y_0_0_0_0_xx_y_xx = buffer_1000_sdpd[114];

    auto g_y_0_0_0_0_xx_y_xy = buffer_1000_sdpd[115];

    auto g_y_0_0_0_0_xx_y_xz = buffer_1000_sdpd[116];

    auto g_y_0_0_0_0_xx_y_yy = buffer_1000_sdpd[117];

    auto g_y_0_0_0_0_xx_y_yz = buffer_1000_sdpd[118];

    auto g_y_0_0_0_0_xx_y_zz = buffer_1000_sdpd[119];

    auto g_y_0_0_0_0_xx_z_xx = buffer_1000_sdpd[120];

    auto g_y_0_0_0_0_xx_z_xy = buffer_1000_sdpd[121];

    auto g_y_0_0_0_0_xx_z_xz = buffer_1000_sdpd[122];

    auto g_y_0_0_0_0_xx_z_yy = buffer_1000_sdpd[123];

    auto g_y_0_0_0_0_xx_z_yz = buffer_1000_sdpd[124];

    auto g_y_0_0_0_0_xx_z_zz = buffer_1000_sdpd[125];

    auto g_y_0_0_0_0_xy_x_xx = buffer_1000_sdpd[126];

    auto g_y_0_0_0_0_xy_x_xy = buffer_1000_sdpd[127];

    auto g_y_0_0_0_0_xy_x_xz = buffer_1000_sdpd[128];

    auto g_y_0_0_0_0_xy_x_yy = buffer_1000_sdpd[129];

    auto g_y_0_0_0_0_xy_x_yz = buffer_1000_sdpd[130];

    auto g_y_0_0_0_0_xy_x_zz = buffer_1000_sdpd[131];

    auto g_y_0_0_0_0_xy_y_xx = buffer_1000_sdpd[132];

    auto g_y_0_0_0_0_xy_y_xy = buffer_1000_sdpd[133];

    auto g_y_0_0_0_0_xy_y_xz = buffer_1000_sdpd[134];

    auto g_y_0_0_0_0_xy_y_yy = buffer_1000_sdpd[135];

    auto g_y_0_0_0_0_xy_y_yz = buffer_1000_sdpd[136];

    auto g_y_0_0_0_0_xy_y_zz = buffer_1000_sdpd[137];

    auto g_y_0_0_0_0_xy_z_xx = buffer_1000_sdpd[138];

    auto g_y_0_0_0_0_xy_z_xy = buffer_1000_sdpd[139];

    auto g_y_0_0_0_0_xy_z_xz = buffer_1000_sdpd[140];

    auto g_y_0_0_0_0_xy_z_yy = buffer_1000_sdpd[141];

    auto g_y_0_0_0_0_xy_z_yz = buffer_1000_sdpd[142];

    auto g_y_0_0_0_0_xy_z_zz = buffer_1000_sdpd[143];

    auto g_y_0_0_0_0_xz_x_xx = buffer_1000_sdpd[144];

    auto g_y_0_0_0_0_xz_x_xy = buffer_1000_sdpd[145];

    auto g_y_0_0_0_0_xz_x_xz = buffer_1000_sdpd[146];

    auto g_y_0_0_0_0_xz_x_yy = buffer_1000_sdpd[147];

    auto g_y_0_0_0_0_xz_x_yz = buffer_1000_sdpd[148];

    auto g_y_0_0_0_0_xz_x_zz = buffer_1000_sdpd[149];

    auto g_y_0_0_0_0_xz_y_xx = buffer_1000_sdpd[150];

    auto g_y_0_0_0_0_xz_y_xy = buffer_1000_sdpd[151];

    auto g_y_0_0_0_0_xz_y_xz = buffer_1000_sdpd[152];

    auto g_y_0_0_0_0_xz_y_yy = buffer_1000_sdpd[153];

    auto g_y_0_0_0_0_xz_y_yz = buffer_1000_sdpd[154];

    auto g_y_0_0_0_0_xz_y_zz = buffer_1000_sdpd[155];

    auto g_y_0_0_0_0_xz_z_xx = buffer_1000_sdpd[156];

    auto g_y_0_0_0_0_xz_z_xy = buffer_1000_sdpd[157];

    auto g_y_0_0_0_0_xz_z_xz = buffer_1000_sdpd[158];

    auto g_y_0_0_0_0_xz_z_yy = buffer_1000_sdpd[159];

    auto g_y_0_0_0_0_xz_z_yz = buffer_1000_sdpd[160];

    auto g_y_0_0_0_0_xz_z_zz = buffer_1000_sdpd[161];

    auto g_y_0_0_0_0_yy_x_xx = buffer_1000_sdpd[162];

    auto g_y_0_0_0_0_yy_x_xy = buffer_1000_sdpd[163];

    auto g_y_0_0_0_0_yy_x_xz = buffer_1000_sdpd[164];

    auto g_y_0_0_0_0_yy_x_yy = buffer_1000_sdpd[165];

    auto g_y_0_0_0_0_yy_x_yz = buffer_1000_sdpd[166];

    auto g_y_0_0_0_0_yy_x_zz = buffer_1000_sdpd[167];

    auto g_y_0_0_0_0_yy_y_xx = buffer_1000_sdpd[168];

    auto g_y_0_0_0_0_yy_y_xy = buffer_1000_sdpd[169];

    auto g_y_0_0_0_0_yy_y_xz = buffer_1000_sdpd[170];

    auto g_y_0_0_0_0_yy_y_yy = buffer_1000_sdpd[171];

    auto g_y_0_0_0_0_yy_y_yz = buffer_1000_sdpd[172];

    auto g_y_0_0_0_0_yy_y_zz = buffer_1000_sdpd[173];

    auto g_y_0_0_0_0_yy_z_xx = buffer_1000_sdpd[174];

    auto g_y_0_0_0_0_yy_z_xy = buffer_1000_sdpd[175];

    auto g_y_0_0_0_0_yy_z_xz = buffer_1000_sdpd[176];

    auto g_y_0_0_0_0_yy_z_yy = buffer_1000_sdpd[177];

    auto g_y_0_0_0_0_yy_z_yz = buffer_1000_sdpd[178];

    auto g_y_0_0_0_0_yy_z_zz = buffer_1000_sdpd[179];

    auto g_y_0_0_0_0_yz_x_xx = buffer_1000_sdpd[180];

    auto g_y_0_0_0_0_yz_x_xy = buffer_1000_sdpd[181];

    auto g_y_0_0_0_0_yz_x_xz = buffer_1000_sdpd[182];

    auto g_y_0_0_0_0_yz_x_yy = buffer_1000_sdpd[183];

    auto g_y_0_0_0_0_yz_x_yz = buffer_1000_sdpd[184];

    auto g_y_0_0_0_0_yz_x_zz = buffer_1000_sdpd[185];

    auto g_y_0_0_0_0_yz_y_xx = buffer_1000_sdpd[186];

    auto g_y_0_0_0_0_yz_y_xy = buffer_1000_sdpd[187];

    auto g_y_0_0_0_0_yz_y_xz = buffer_1000_sdpd[188];

    auto g_y_0_0_0_0_yz_y_yy = buffer_1000_sdpd[189];

    auto g_y_0_0_0_0_yz_y_yz = buffer_1000_sdpd[190];

    auto g_y_0_0_0_0_yz_y_zz = buffer_1000_sdpd[191];

    auto g_y_0_0_0_0_yz_z_xx = buffer_1000_sdpd[192];

    auto g_y_0_0_0_0_yz_z_xy = buffer_1000_sdpd[193];

    auto g_y_0_0_0_0_yz_z_xz = buffer_1000_sdpd[194];

    auto g_y_0_0_0_0_yz_z_yy = buffer_1000_sdpd[195];

    auto g_y_0_0_0_0_yz_z_yz = buffer_1000_sdpd[196];

    auto g_y_0_0_0_0_yz_z_zz = buffer_1000_sdpd[197];

    auto g_y_0_0_0_0_zz_x_xx = buffer_1000_sdpd[198];

    auto g_y_0_0_0_0_zz_x_xy = buffer_1000_sdpd[199];

    auto g_y_0_0_0_0_zz_x_xz = buffer_1000_sdpd[200];

    auto g_y_0_0_0_0_zz_x_yy = buffer_1000_sdpd[201];

    auto g_y_0_0_0_0_zz_x_yz = buffer_1000_sdpd[202];

    auto g_y_0_0_0_0_zz_x_zz = buffer_1000_sdpd[203];

    auto g_y_0_0_0_0_zz_y_xx = buffer_1000_sdpd[204];

    auto g_y_0_0_0_0_zz_y_xy = buffer_1000_sdpd[205];

    auto g_y_0_0_0_0_zz_y_xz = buffer_1000_sdpd[206];

    auto g_y_0_0_0_0_zz_y_yy = buffer_1000_sdpd[207];

    auto g_y_0_0_0_0_zz_y_yz = buffer_1000_sdpd[208];

    auto g_y_0_0_0_0_zz_y_zz = buffer_1000_sdpd[209];

    auto g_y_0_0_0_0_zz_z_xx = buffer_1000_sdpd[210];

    auto g_y_0_0_0_0_zz_z_xy = buffer_1000_sdpd[211];

    auto g_y_0_0_0_0_zz_z_xz = buffer_1000_sdpd[212];

    auto g_y_0_0_0_0_zz_z_yy = buffer_1000_sdpd[213];

    auto g_y_0_0_0_0_zz_z_yz = buffer_1000_sdpd[214];

    auto g_y_0_0_0_0_zz_z_zz = buffer_1000_sdpd[215];

    auto g_z_0_0_0_0_xx_x_xx = buffer_1000_sdpd[216];

    auto g_z_0_0_0_0_xx_x_xy = buffer_1000_sdpd[217];

    auto g_z_0_0_0_0_xx_x_xz = buffer_1000_sdpd[218];

    auto g_z_0_0_0_0_xx_x_yy = buffer_1000_sdpd[219];

    auto g_z_0_0_0_0_xx_x_yz = buffer_1000_sdpd[220];

    auto g_z_0_0_0_0_xx_x_zz = buffer_1000_sdpd[221];

    auto g_z_0_0_0_0_xx_y_xx = buffer_1000_sdpd[222];

    auto g_z_0_0_0_0_xx_y_xy = buffer_1000_sdpd[223];

    auto g_z_0_0_0_0_xx_y_xz = buffer_1000_sdpd[224];

    auto g_z_0_0_0_0_xx_y_yy = buffer_1000_sdpd[225];

    auto g_z_0_0_0_0_xx_y_yz = buffer_1000_sdpd[226];

    auto g_z_0_0_0_0_xx_y_zz = buffer_1000_sdpd[227];

    auto g_z_0_0_0_0_xx_z_xx = buffer_1000_sdpd[228];

    auto g_z_0_0_0_0_xx_z_xy = buffer_1000_sdpd[229];

    auto g_z_0_0_0_0_xx_z_xz = buffer_1000_sdpd[230];

    auto g_z_0_0_0_0_xx_z_yy = buffer_1000_sdpd[231];

    auto g_z_0_0_0_0_xx_z_yz = buffer_1000_sdpd[232];

    auto g_z_0_0_0_0_xx_z_zz = buffer_1000_sdpd[233];

    auto g_z_0_0_0_0_xy_x_xx = buffer_1000_sdpd[234];

    auto g_z_0_0_0_0_xy_x_xy = buffer_1000_sdpd[235];

    auto g_z_0_0_0_0_xy_x_xz = buffer_1000_sdpd[236];

    auto g_z_0_0_0_0_xy_x_yy = buffer_1000_sdpd[237];

    auto g_z_0_0_0_0_xy_x_yz = buffer_1000_sdpd[238];

    auto g_z_0_0_0_0_xy_x_zz = buffer_1000_sdpd[239];

    auto g_z_0_0_0_0_xy_y_xx = buffer_1000_sdpd[240];

    auto g_z_0_0_0_0_xy_y_xy = buffer_1000_sdpd[241];

    auto g_z_0_0_0_0_xy_y_xz = buffer_1000_sdpd[242];

    auto g_z_0_0_0_0_xy_y_yy = buffer_1000_sdpd[243];

    auto g_z_0_0_0_0_xy_y_yz = buffer_1000_sdpd[244];

    auto g_z_0_0_0_0_xy_y_zz = buffer_1000_sdpd[245];

    auto g_z_0_0_0_0_xy_z_xx = buffer_1000_sdpd[246];

    auto g_z_0_0_0_0_xy_z_xy = buffer_1000_sdpd[247];

    auto g_z_0_0_0_0_xy_z_xz = buffer_1000_sdpd[248];

    auto g_z_0_0_0_0_xy_z_yy = buffer_1000_sdpd[249];

    auto g_z_0_0_0_0_xy_z_yz = buffer_1000_sdpd[250];

    auto g_z_0_0_0_0_xy_z_zz = buffer_1000_sdpd[251];

    auto g_z_0_0_0_0_xz_x_xx = buffer_1000_sdpd[252];

    auto g_z_0_0_0_0_xz_x_xy = buffer_1000_sdpd[253];

    auto g_z_0_0_0_0_xz_x_xz = buffer_1000_sdpd[254];

    auto g_z_0_0_0_0_xz_x_yy = buffer_1000_sdpd[255];

    auto g_z_0_0_0_0_xz_x_yz = buffer_1000_sdpd[256];

    auto g_z_0_0_0_0_xz_x_zz = buffer_1000_sdpd[257];

    auto g_z_0_0_0_0_xz_y_xx = buffer_1000_sdpd[258];

    auto g_z_0_0_0_0_xz_y_xy = buffer_1000_sdpd[259];

    auto g_z_0_0_0_0_xz_y_xz = buffer_1000_sdpd[260];

    auto g_z_0_0_0_0_xz_y_yy = buffer_1000_sdpd[261];

    auto g_z_0_0_0_0_xz_y_yz = buffer_1000_sdpd[262];

    auto g_z_0_0_0_0_xz_y_zz = buffer_1000_sdpd[263];

    auto g_z_0_0_0_0_xz_z_xx = buffer_1000_sdpd[264];

    auto g_z_0_0_0_0_xz_z_xy = buffer_1000_sdpd[265];

    auto g_z_0_0_0_0_xz_z_xz = buffer_1000_sdpd[266];

    auto g_z_0_0_0_0_xz_z_yy = buffer_1000_sdpd[267];

    auto g_z_0_0_0_0_xz_z_yz = buffer_1000_sdpd[268];

    auto g_z_0_0_0_0_xz_z_zz = buffer_1000_sdpd[269];

    auto g_z_0_0_0_0_yy_x_xx = buffer_1000_sdpd[270];

    auto g_z_0_0_0_0_yy_x_xy = buffer_1000_sdpd[271];

    auto g_z_0_0_0_0_yy_x_xz = buffer_1000_sdpd[272];

    auto g_z_0_0_0_0_yy_x_yy = buffer_1000_sdpd[273];

    auto g_z_0_0_0_0_yy_x_yz = buffer_1000_sdpd[274];

    auto g_z_0_0_0_0_yy_x_zz = buffer_1000_sdpd[275];

    auto g_z_0_0_0_0_yy_y_xx = buffer_1000_sdpd[276];

    auto g_z_0_0_0_0_yy_y_xy = buffer_1000_sdpd[277];

    auto g_z_0_0_0_0_yy_y_xz = buffer_1000_sdpd[278];

    auto g_z_0_0_0_0_yy_y_yy = buffer_1000_sdpd[279];

    auto g_z_0_0_0_0_yy_y_yz = buffer_1000_sdpd[280];

    auto g_z_0_0_0_0_yy_y_zz = buffer_1000_sdpd[281];

    auto g_z_0_0_0_0_yy_z_xx = buffer_1000_sdpd[282];

    auto g_z_0_0_0_0_yy_z_xy = buffer_1000_sdpd[283];

    auto g_z_0_0_0_0_yy_z_xz = buffer_1000_sdpd[284];

    auto g_z_0_0_0_0_yy_z_yy = buffer_1000_sdpd[285];

    auto g_z_0_0_0_0_yy_z_yz = buffer_1000_sdpd[286];

    auto g_z_0_0_0_0_yy_z_zz = buffer_1000_sdpd[287];

    auto g_z_0_0_0_0_yz_x_xx = buffer_1000_sdpd[288];

    auto g_z_0_0_0_0_yz_x_xy = buffer_1000_sdpd[289];

    auto g_z_0_0_0_0_yz_x_xz = buffer_1000_sdpd[290];

    auto g_z_0_0_0_0_yz_x_yy = buffer_1000_sdpd[291];

    auto g_z_0_0_0_0_yz_x_yz = buffer_1000_sdpd[292];

    auto g_z_0_0_0_0_yz_x_zz = buffer_1000_sdpd[293];

    auto g_z_0_0_0_0_yz_y_xx = buffer_1000_sdpd[294];

    auto g_z_0_0_0_0_yz_y_xy = buffer_1000_sdpd[295];

    auto g_z_0_0_0_0_yz_y_xz = buffer_1000_sdpd[296];

    auto g_z_0_0_0_0_yz_y_yy = buffer_1000_sdpd[297];

    auto g_z_0_0_0_0_yz_y_yz = buffer_1000_sdpd[298];

    auto g_z_0_0_0_0_yz_y_zz = buffer_1000_sdpd[299];

    auto g_z_0_0_0_0_yz_z_xx = buffer_1000_sdpd[300];

    auto g_z_0_0_0_0_yz_z_xy = buffer_1000_sdpd[301];

    auto g_z_0_0_0_0_yz_z_xz = buffer_1000_sdpd[302];

    auto g_z_0_0_0_0_yz_z_yy = buffer_1000_sdpd[303];

    auto g_z_0_0_0_0_yz_z_yz = buffer_1000_sdpd[304];

    auto g_z_0_0_0_0_yz_z_zz = buffer_1000_sdpd[305];

    auto g_z_0_0_0_0_zz_x_xx = buffer_1000_sdpd[306];

    auto g_z_0_0_0_0_zz_x_xy = buffer_1000_sdpd[307];

    auto g_z_0_0_0_0_zz_x_xz = buffer_1000_sdpd[308];

    auto g_z_0_0_0_0_zz_x_yy = buffer_1000_sdpd[309];

    auto g_z_0_0_0_0_zz_x_yz = buffer_1000_sdpd[310];

    auto g_z_0_0_0_0_zz_x_zz = buffer_1000_sdpd[311];

    auto g_z_0_0_0_0_zz_y_xx = buffer_1000_sdpd[312];

    auto g_z_0_0_0_0_zz_y_xy = buffer_1000_sdpd[313];

    auto g_z_0_0_0_0_zz_y_xz = buffer_1000_sdpd[314];

    auto g_z_0_0_0_0_zz_y_yy = buffer_1000_sdpd[315];

    auto g_z_0_0_0_0_zz_y_yz = buffer_1000_sdpd[316];

    auto g_z_0_0_0_0_zz_y_zz = buffer_1000_sdpd[317];

    auto g_z_0_0_0_0_zz_z_xx = buffer_1000_sdpd[318];

    auto g_z_0_0_0_0_zz_z_xy = buffer_1000_sdpd[319];

    auto g_z_0_0_0_0_zz_z_xz = buffer_1000_sdpd[320];

    auto g_z_0_0_0_0_zz_z_yy = buffer_1000_sdpd[321];

    auto g_z_0_0_0_0_zz_z_yz = buffer_1000_sdpd[322];

    auto g_z_0_0_0_0_zz_z_zz = buffer_1000_sdpd[323];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_0_0_0_xx_x_xx, g_x_0_0_0_0_xx_x_xy, g_x_0_0_0_0_xx_x_xz, g_x_0_0_0_0_xx_x_yy, g_x_0_0_0_0_xx_x_yz, g_x_0_0_0_0_xx_x_zz, g_x_xx_x_xx, g_x_xx_x_xy, g_x_xx_x_xz, g_x_xx_x_yy, g_x_xx_x_yz, g_x_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xx_x_xx[i] = 2.0 * g_x_xx_x_xx[i] * a_exp;

        g_x_0_0_0_0_xx_x_xy[i] = 2.0 * g_x_xx_x_xy[i] * a_exp;

        g_x_0_0_0_0_xx_x_xz[i] = 2.0 * g_x_xx_x_xz[i] * a_exp;

        g_x_0_0_0_0_xx_x_yy[i] = 2.0 * g_x_xx_x_yy[i] * a_exp;

        g_x_0_0_0_0_xx_x_yz[i] = 2.0 * g_x_xx_x_yz[i] * a_exp;

        g_x_0_0_0_0_xx_x_zz[i] = 2.0 * g_x_xx_x_zz[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_0_0_0_xx_y_xx, g_x_0_0_0_0_xx_y_xy, g_x_0_0_0_0_xx_y_xz, g_x_0_0_0_0_xx_y_yy, g_x_0_0_0_0_xx_y_yz, g_x_0_0_0_0_xx_y_zz, g_x_xx_y_xx, g_x_xx_y_xy, g_x_xx_y_xz, g_x_xx_y_yy, g_x_xx_y_yz, g_x_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xx_y_xx[i] = 2.0 * g_x_xx_y_xx[i] * a_exp;

        g_x_0_0_0_0_xx_y_xy[i] = 2.0 * g_x_xx_y_xy[i] * a_exp;

        g_x_0_0_0_0_xx_y_xz[i] = 2.0 * g_x_xx_y_xz[i] * a_exp;

        g_x_0_0_0_0_xx_y_yy[i] = 2.0 * g_x_xx_y_yy[i] * a_exp;

        g_x_0_0_0_0_xx_y_yz[i] = 2.0 * g_x_xx_y_yz[i] * a_exp;

        g_x_0_0_0_0_xx_y_zz[i] = 2.0 * g_x_xx_y_zz[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_0_0_0_xx_z_xx, g_x_0_0_0_0_xx_z_xy, g_x_0_0_0_0_xx_z_xz, g_x_0_0_0_0_xx_z_yy, g_x_0_0_0_0_xx_z_yz, g_x_0_0_0_0_xx_z_zz, g_x_xx_z_xx, g_x_xx_z_xy, g_x_xx_z_xz, g_x_xx_z_yy, g_x_xx_z_yz, g_x_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xx_z_xx[i] = 2.0 * g_x_xx_z_xx[i] * a_exp;

        g_x_0_0_0_0_xx_z_xy[i] = 2.0 * g_x_xx_z_xy[i] * a_exp;

        g_x_0_0_0_0_xx_z_xz[i] = 2.0 * g_x_xx_z_xz[i] * a_exp;

        g_x_0_0_0_0_xx_z_yy[i] = 2.0 * g_x_xx_z_yy[i] * a_exp;

        g_x_0_0_0_0_xx_z_yz[i] = 2.0 * g_x_xx_z_yz[i] * a_exp;

        g_x_0_0_0_0_xx_z_zz[i] = 2.0 * g_x_xx_z_zz[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_0_0_0_xy_x_xx, g_x_0_0_0_0_xy_x_xy, g_x_0_0_0_0_xy_x_xz, g_x_0_0_0_0_xy_x_yy, g_x_0_0_0_0_xy_x_yz, g_x_0_0_0_0_xy_x_zz, g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xy_x_xx[i] = 2.0 * g_x_xy_x_xx[i] * a_exp;

        g_x_0_0_0_0_xy_x_xy[i] = 2.0 * g_x_xy_x_xy[i] * a_exp;

        g_x_0_0_0_0_xy_x_xz[i] = 2.0 * g_x_xy_x_xz[i] * a_exp;

        g_x_0_0_0_0_xy_x_yy[i] = 2.0 * g_x_xy_x_yy[i] * a_exp;

        g_x_0_0_0_0_xy_x_yz[i] = 2.0 * g_x_xy_x_yz[i] * a_exp;

        g_x_0_0_0_0_xy_x_zz[i] = 2.0 * g_x_xy_x_zz[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_0_0_0_xy_y_xx, g_x_0_0_0_0_xy_y_xy, g_x_0_0_0_0_xy_y_xz, g_x_0_0_0_0_xy_y_yy, g_x_0_0_0_0_xy_y_yz, g_x_0_0_0_0_xy_y_zz, g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xy_y_xx[i] = 2.0 * g_x_xy_y_xx[i] * a_exp;

        g_x_0_0_0_0_xy_y_xy[i] = 2.0 * g_x_xy_y_xy[i] * a_exp;

        g_x_0_0_0_0_xy_y_xz[i] = 2.0 * g_x_xy_y_xz[i] * a_exp;

        g_x_0_0_0_0_xy_y_yy[i] = 2.0 * g_x_xy_y_yy[i] * a_exp;

        g_x_0_0_0_0_xy_y_yz[i] = 2.0 * g_x_xy_y_yz[i] * a_exp;

        g_x_0_0_0_0_xy_y_zz[i] = 2.0 * g_x_xy_y_zz[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_0_0_0_xy_z_xx, g_x_0_0_0_0_xy_z_xy, g_x_0_0_0_0_xy_z_xz, g_x_0_0_0_0_xy_z_yy, g_x_0_0_0_0_xy_z_yz, g_x_0_0_0_0_xy_z_zz, g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xy_z_xx[i] = 2.0 * g_x_xy_z_xx[i] * a_exp;

        g_x_0_0_0_0_xy_z_xy[i] = 2.0 * g_x_xy_z_xy[i] * a_exp;

        g_x_0_0_0_0_xy_z_xz[i] = 2.0 * g_x_xy_z_xz[i] * a_exp;

        g_x_0_0_0_0_xy_z_yy[i] = 2.0 * g_x_xy_z_yy[i] * a_exp;

        g_x_0_0_0_0_xy_z_yz[i] = 2.0 * g_x_xy_z_yz[i] * a_exp;

        g_x_0_0_0_0_xy_z_zz[i] = 2.0 * g_x_xy_z_zz[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_0_0_0_xz_x_xx, g_x_0_0_0_0_xz_x_xy, g_x_0_0_0_0_xz_x_xz, g_x_0_0_0_0_xz_x_yy, g_x_0_0_0_0_xz_x_yz, g_x_0_0_0_0_xz_x_zz, g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xz_x_xx[i] = 2.0 * g_x_xz_x_xx[i] * a_exp;

        g_x_0_0_0_0_xz_x_xy[i] = 2.0 * g_x_xz_x_xy[i] * a_exp;

        g_x_0_0_0_0_xz_x_xz[i] = 2.0 * g_x_xz_x_xz[i] * a_exp;

        g_x_0_0_0_0_xz_x_yy[i] = 2.0 * g_x_xz_x_yy[i] * a_exp;

        g_x_0_0_0_0_xz_x_yz[i] = 2.0 * g_x_xz_x_yz[i] * a_exp;

        g_x_0_0_0_0_xz_x_zz[i] = 2.0 * g_x_xz_x_zz[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_0_0_0_xz_y_xx, g_x_0_0_0_0_xz_y_xy, g_x_0_0_0_0_xz_y_xz, g_x_0_0_0_0_xz_y_yy, g_x_0_0_0_0_xz_y_yz, g_x_0_0_0_0_xz_y_zz, g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xz_y_xx[i] = 2.0 * g_x_xz_y_xx[i] * a_exp;

        g_x_0_0_0_0_xz_y_xy[i] = 2.0 * g_x_xz_y_xy[i] * a_exp;

        g_x_0_0_0_0_xz_y_xz[i] = 2.0 * g_x_xz_y_xz[i] * a_exp;

        g_x_0_0_0_0_xz_y_yy[i] = 2.0 * g_x_xz_y_yy[i] * a_exp;

        g_x_0_0_0_0_xz_y_yz[i] = 2.0 * g_x_xz_y_yz[i] * a_exp;

        g_x_0_0_0_0_xz_y_zz[i] = 2.0 * g_x_xz_y_zz[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_0_0_0_xz_z_xx, g_x_0_0_0_0_xz_z_xy, g_x_0_0_0_0_xz_z_xz, g_x_0_0_0_0_xz_z_yy, g_x_0_0_0_0_xz_z_yz, g_x_0_0_0_0_xz_z_zz, g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xz_z_xx[i] = 2.0 * g_x_xz_z_xx[i] * a_exp;

        g_x_0_0_0_0_xz_z_xy[i] = 2.0 * g_x_xz_z_xy[i] * a_exp;

        g_x_0_0_0_0_xz_z_xz[i] = 2.0 * g_x_xz_z_xz[i] * a_exp;

        g_x_0_0_0_0_xz_z_yy[i] = 2.0 * g_x_xz_z_yy[i] * a_exp;

        g_x_0_0_0_0_xz_z_yz[i] = 2.0 * g_x_xz_z_yz[i] * a_exp;

        g_x_0_0_0_0_xz_z_zz[i] = 2.0 * g_x_xz_z_zz[i] * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_0_0_0_0_yy_x_xx, g_x_0_0_0_0_yy_x_xy, g_x_0_0_0_0_yy_x_xz, g_x_0_0_0_0_yy_x_yy, g_x_0_0_0_0_yy_x_yz, g_x_0_0_0_0_yy_x_zz, g_x_yy_x_xx, g_x_yy_x_xy, g_x_yy_x_xz, g_x_yy_x_yy, g_x_yy_x_yz, g_x_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yy_x_xx[i] = 2.0 * g_x_yy_x_xx[i] * a_exp;

        g_x_0_0_0_0_yy_x_xy[i] = 2.0 * g_x_yy_x_xy[i] * a_exp;

        g_x_0_0_0_0_yy_x_xz[i] = 2.0 * g_x_yy_x_xz[i] * a_exp;

        g_x_0_0_0_0_yy_x_yy[i] = 2.0 * g_x_yy_x_yy[i] * a_exp;

        g_x_0_0_0_0_yy_x_yz[i] = 2.0 * g_x_yy_x_yz[i] * a_exp;

        g_x_0_0_0_0_yy_x_zz[i] = 2.0 * g_x_yy_x_zz[i] * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_0_0_0_0_yy_y_xx, g_x_0_0_0_0_yy_y_xy, g_x_0_0_0_0_yy_y_xz, g_x_0_0_0_0_yy_y_yy, g_x_0_0_0_0_yy_y_yz, g_x_0_0_0_0_yy_y_zz, g_x_yy_y_xx, g_x_yy_y_xy, g_x_yy_y_xz, g_x_yy_y_yy, g_x_yy_y_yz, g_x_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yy_y_xx[i] = 2.0 * g_x_yy_y_xx[i] * a_exp;

        g_x_0_0_0_0_yy_y_xy[i] = 2.0 * g_x_yy_y_xy[i] * a_exp;

        g_x_0_0_0_0_yy_y_xz[i] = 2.0 * g_x_yy_y_xz[i] * a_exp;

        g_x_0_0_0_0_yy_y_yy[i] = 2.0 * g_x_yy_y_yy[i] * a_exp;

        g_x_0_0_0_0_yy_y_yz[i] = 2.0 * g_x_yy_y_yz[i] * a_exp;

        g_x_0_0_0_0_yy_y_zz[i] = 2.0 * g_x_yy_y_zz[i] * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_0_0_0_0_yy_z_xx, g_x_0_0_0_0_yy_z_xy, g_x_0_0_0_0_yy_z_xz, g_x_0_0_0_0_yy_z_yy, g_x_0_0_0_0_yy_z_yz, g_x_0_0_0_0_yy_z_zz, g_x_yy_z_xx, g_x_yy_z_xy, g_x_yy_z_xz, g_x_yy_z_yy, g_x_yy_z_yz, g_x_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yy_z_xx[i] = 2.0 * g_x_yy_z_xx[i] * a_exp;

        g_x_0_0_0_0_yy_z_xy[i] = 2.0 * g_x_yy_z_xy[i] * a_exp;

        g_x_0_0_0_0_yy_z_xz[i] = 2.0 * g_x_yy_z_xz[i] * a_exp;

        g_x_0_0_0_0_yy_z_yy[i] = 2.0 * g_x_yy_z_yy[i] * a_exp;

        g_x_0_0_0_0_yy_z_yz[i] = 2.0 * g_x_yy_z_yz[i] * a_exp;

        g_x_0_0_0_0_yy_z_zz[i] = 2.0 * g_x_yy_z_zz[i] * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_0_0_0_yz_x_xx, g_x_0_0_0_0_yz_x_xy, g_x_0_0_0_0_yz_x_xz, g_x_0_0_0_0_yz_x_yy, g_x_0_0_0_0_yz_x_yz, g_x_0_0_0_0_yz_x_zz, g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yz_x_xx[i] = 2.0 * g_x_yz_x_xx[i] * a_exp;

        g_x_0_0_0_0_yz_x_xy[i] = 2.0 * g_x_yz_x_xy[i] * a_exp;

        g_x_0_0_0_0_yz_x_xz[i] = 2.0 * g_x_yz_x_xz[i] * a_exp;

        g_x_0_0_0_0_yz_x_yy[i] = 2.0 * g_x_yz_x_yy[i] * a_exp;

        g_x_0_0_0_0_yz_x_yz[i] = 2.0 * g_x_yz_x_yz[i] * a_exp;

        g_x_0_0_0_0_yz_x_zz[i] = 2.0 * g_x_yz_x_zz[i] * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_0_0_0_yz_y_xx, g_x_0_0_0_0_yz_y_xy, g_x_0_0_0_0_yz_y_xz, g_x_0_0_0_0_yz_y_yy, g_x_0_0_0_0_yz_y_yz, g_x_0_0_0_0_yz_y_zz, g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yz_y_xx[i] = 2.0 * g_x_yz_y_xx[i] * a_exp;

        g_x_0_0_0_0_yz_y_xy[i] = 2.0 * g_x_yz_y_xy[i] * a_exp;

        g_x_0_0_0_0_yz_y_xz[i] = 2.0 * g_x_yz_y_xz[i] * a_exp;

        g_x_0_0_0_0_yz_y_yy[i] = 2.0 * g_x_yz_y_yy[i] * a_exp;

        g_x_0_0_0_0_yz_y_yz[i] = 2.0 * g_x_yz_y_yz[i] * a_exp;

        g_x_0_0_0_0_yz_y_zz[i] = 2.0 * g_x_yz_y_zz[i] * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_0_0_0_yz_z_xx, g_x_0_0_0_0_yz_z_xy, g_x_0_0_0_0_yz_z_xz, g_x_0_0_0_0_yz_z_yy, g_x_0_0_0_0_yz_z_yz, g_x_0_0_0_0_yz_z_zz, g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yz_z_xx[i] = 2.0 * g_x_yz_z_xx[i] * a_exp;

        g_x_0_0_0_0_yz_z_xy[i] = 2.0 * g_x_yz_z_xy[i] * a_exp;

        g_x_0_0_0_0_yz_z_xz[i] = 2.0 * g_x_yz_z_xz[i] * a_exp;

        g_x_0_0_0_0_yz_z_yy[i] = 2.0 * g_x_yz_z_yy[i] * a_exp;

        g_x_0_0_0_0_yz_z_yz[i] = 2.0 * g_x_yz_z_yz[i] * a_exp;

        g_x_0_0_0_0_yz_z_zz[i] = 2.0 * g_x_yz_z_zz[i] * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_0_0_0_zz_x_xx, g_x_0_0_0_0_zz_x_xy, g_x_0_0_0_0_zz_x_xz, g_x_0_0_0_0_zz_x_yy, g_x_0_0_0_0_zz_x_yz, g_x_0_0_0_0_zz_x_zz, g_x_zz_x_xx, g_x_zz_x_xy, g_x_zz_x_xz, g_x_zz_x_yy, g_x_zz_x_yz, g_x_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_zz_x_xx[i] = 2.0 * g_x_zz_x_xx[i] * a_exp;

        g_x_0_0_0_0_zz_x_xy[i] = 2.0 * g_x_zz_x_xy[i] * a_exp;

        g_x_0_0_0_0_zz_x_xz[i] = 2.0 * g_x_zz_x_xz[i] * a_exp;

        g_x_0_0_0_0_zz_x_yy[i] = 2.0 * g_x_zz_x_yy[i] * a_exp;

        g_x_0_0_0_0_zz_x_yz[i] = 2.0 * g_x_zz_x_yz[i] * a_exp;

        g_x_0_0_0_0_zz_x_zz[i] = 2.0 * g_x_zz_x_zz[i] * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_0_0_0_zz_y_xx, g_x_0_0_0_0_zz_y_xy, g_x_0_0_0_0_zz_y_xz, g_x_0_0_0_0_zz_y_yy, g_x_0_0_0_0_zz_y_yz, g_x_0_0_0_0_zz_y_zz, g_x_zz_y_xx, g_x_zz_y_xy, g_x_zz_y_xz, g_x_zz_y_yy, g_x_zz_y_yz, g_x_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_zz_y_xx[i] = 2.0 * g_x_zz_y_xx[i] * a_exp;

        g_x_0_0_0_0_zz_y_xy[i] = 2.0 * g_x_zz_y_xy[i] * a_exp;

        g_x_0_0_0_0_zz_y_xz[i] = 2.0 * g_x_zz_y_xz[i] * a_exp;

        g_x_0_0_0_0_zz_y_yy[i] = 2.0 * g_x_zz_y_yy[i] * a_exp;

        g_x_0_0_0_0_zz_y_yz[i] = 2.0 * g_x_zz_y_yz[i] * a_exp;

        g_x_0_0_0_0_zz_y_zz[i] = 2.0 * g_x_zz_y_zz[i] * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_0_0_0_zz_z_xx, g_x_0_0_0_0_zz_z_xy, g_x_0_0_0_0_zz_z_xz, g_x_0_0_0_0_zz_z_yy, g_x_0_0_0_0_zz_z_yz, g_x_0_0_0_0_zz_z_zz, g_x_zz_z_xx, g_x_zz_z_xy, g_x_zz_z_xz, g_x_zz_z_yy, g_x_zz_z_yz, g_x_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_zz_z_xx[i] = 2.0 * g_x_zz_z_xx[i] * a_exp;

        g_x_0_0_0_0_zz_z_xy[i] = 2.0 * g_x_zz_z_xy[i] * a_exp;

        g_x_0_0_0_0_zz_z_xz[i] = 2.0 * g_x_zz_z_xz[i] * a_exp;

        g_x_0_0_0_0_zz_z_yy[i] = 2.0 * g_x_zz_z_yy[i] * a_exp;

        g_x_0_0_0_0_zz_z_yz[i] = 2.0 * g_x_zz_z_yz[i] * a_exp;

        g_x_0_0_0_0_zz_z_zz[i] = 2.0 * g_x_zz_z_zz[i] * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_y_0_0_0_0_xx_x_xx, g_y_0_0_0_0_xx_x_xy, g_y_0_0_0_0_xx_x_xz, g_y_0_0_0_0_xx_x_yy, g_y_0_0_0_0_xx_x_yz, g_y_0_0_0_0_xx_x_zz, g_y_xx_x_xx, g_y_xx_x_xy, g_y_xx_x_xz, g_y_xx_x_yy, g_y_xx_x_yz, g_y_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xx_x_xx[i] = 2.0 * g_y_xx_x_xx[i] * a_exp;

        g_y_0_0_0_0_xx_x_xy[i] = 2.0 * g_y_xx_x_xy[i] * a_exp;

        g_y_0_0_0_0_xx_x_xz[i] = 2.0 * g_y_xx_x_xz[i] * a_exp;

        g_y_0_0_0_0_xx_x_yy[i] = 2.0 * g_y_xx_x_yy[i] * a_exp;

        g_y_0_0_0_0_xx_x_yz[i] = 2.0 * g_y_xx_x_yz[i] * a_exp;

        g_y_0_0_0_0_xx_x_zz[i] = 2.0 * g_y_xx_x_zz[i] * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_y_0_0_0_0_xx_y_xx, g_y_0_0_0_0_xx_y_xy, g_y_0_0_0_0_xx_y_xz, g_y_0_0_0_0_xx_y_yy, g_y_0_0_0_0_xx_y_yz, g_y_0_0_0_0_xx_y_zz, g_y_xx_y_xx, g_y_xx_y_xy, g_y_xx_y_xz, g_y_xx_y_yy, g_y_xx_y_yz, g_y_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xx_y_xx[i] = 2.0 * g_y_xx_y_xx[i] * a_exp;

        g_y_0_0_0_0_xx_y_xy[i] = 2.0 * g_y_xx_y_xy[i] * a_exp;

        g_y_0_0_0_0_xx_y_xz[i] = 2.0 * g_y_xx_y_xz[i] * a_exp;

        g_y_0_0_0_0_xx_y_yy[i] = 2.0 * g_y_xx_y_yy[i] * a_exp;

        g_y_0_0_0_0_xx_y_yz[i] = 2.0 * g_y_xx_y_yz[i] * a_exp;

        g_y_0_0_0_0_xx_y_zz[i] = 2.0 * g_y_xx_y_zz[i] * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_y_0_0_0_0_xx_z_xx, g_y_0_0_0_0_xx_z_xy, g_y_0_0_0_0_xx_z_xz, g_y_0_0_0_0_xx_z_yy, g_y_0_0_0_0_xx_z_yz, g_y_0_0_0_0_xx_z_zz, g_y_xx_z_xx, g_y_xx_z_xy, g_y_xx_z_xz, g_y_xx_z_yy, g_y_xx_z_yz, g_y_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xx_z_xx[i] = 2.0 * g_y_xx_z_xx[i] * a_exp;

        g_y_0_0_0_0_xx_z_xy[i] = 2.0 * g_y_xx_z_xy[i] * a_exp;

        g_y_0_0_0_0_xx_z_xz[i] = 2.0 * g_y_xx_z_xz[i] * a_exp;

        g_y_0_0_0_0_xx_z_yy[i] = 2.0 * g_y_xx_z_yy[i] * a_exp;

        g_y_0_0_0_0_xx_z_yz[i] = 2.0 * g_y_xx_z_yz[i] * a_exp;

        g_y_0_0_0_0_xx_z_zz[i] = 2.0 * g_y_xx_z_zz[i] * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_y_0_0_0_0_xy_x_xx, g_y_0_0_0_0_xy_x_xy, g_y_0_0_0_0_xy_x_xz, g_y_0_0_0_0_xy_x_yy, g_y_0_0_0_0_xy_x_yz, g_y_0_0_0_0_xy_x_zz, g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xy_x_xx[i] = 2.0 * g_y_xy_x_xx[i] * a_exp;

        g_y_0_0_0_0_xy_x_xy[i] = 2.0 * g_y_xy_x_xy[i] * a_exp;

        g_y_0_0_0_0_xy_x_xz[i] = 2.0 * g_y_xy_x_xz[i] * a_exp;

        g_y_0_0_0_0_xy_x_yy[i] = 2.0 * g_y_xy_x_yy[i] * a_exp;

        g_y_0_0_0_0_xy_x_yz[i] = 2.0 * g_y_xy_x_yz[i] * a_exp;

        g_y_0_0_0_0_xy_x_zz[i] = 2.0 * g_y_xy_x_zz[i] * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_y_0_0_0_0_xy_y_xx, g_y_0_0_0_0_xy_y_xy, g_y_0_0_0_0_xy_y_xz, g_y_0_0_0_0_xy_y_yy, g_y_0_0_0_0_xy_y_yz, g_y_0_0_0_0_xy_y_zz, g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xy_y_xx[i] = 2.0 * g_y_xy_y_xx[i] * a_exp;

        g_y_0_0_0_0_xy_y_xy[i] = 2.0 * g_y_xy_y_xy[i] * a_exp;

        g_y_0_0_0_0_xy_y_xz[i] = 2.0 * g_y_xy_y_xz[i] * a_exp;

        g_y_0_0_0_0_xy_y_yy[i] = 2.0 * g_y_xy_y_yy[i] * a_exp;

        g_y_0_0_0_0_xy_y_yz[i] = 2.0 * g_y_xy_y_yz[i] * a_exp;

        g_y_0_0_0_0_xy_y_zz[i] = 2.0 * g_y_xy_y_zz[i] * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_y_0_0_0_0_xy_z_xx, g_y_0_0_0_0_xy_z_xy, g_y_0_0_0_0_xy_z_xz, g_y_0_0_0_0_xy_z_yy, g_y_0_0_0_0_xy_z_yz, g_y_0_0_0_0_xy_z_zz, g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xy_z_xx[i] = 2.0 * g_y_xy_z_xx[i] * a_exp;

        g_y_0_0_0_0_xy_z_xy[i] = 2.0 * g_y_xy_z_xy[i] * a_exp;

        g_y_0_0_0_0_xy_z_xz[i] = 2.0 * g_y_xy_z_xz[i] * a_exp;

        g_y_0_0_0_0_xy_z_yy[i] = 2.0 * g_y_xy_z_yy[i] * a_exp;

        g_y_0_0_0_0_xy_z_yz[i] = 2.0 * g_y_xy_z_yz[i] * a_exp;

        g_y_0_0_0_0_xy_z_zz[i] = 2.0 * g_y_xy_z_zz[i] * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_y_0_0_0_0_xz_x_xx, g_y_0_0_0_0_xz_x_xy, g_y_0_0_0_0_xz_x_xz, g_y_0_0_0_0_xz_x_yy, g_y_0_0_0_0_xz_x_yz, g_y_0_0_0_0_xz_x_zz, g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xz_x_xx[i] = 2.0 * g_y_xz_x_xx[i] * a_exp;

        g_y_0_0_0_0_xz_x_xy[i] = 2.0 * g_y_xz_x_xy[i] * a_exp;

        g_y_0_0_0_0_xz_x_xz[i] = 2.0 * g_y_xz_x_xz[i] * a_exp;

        g_y_0_0_0_0_xz_x_yy[i] = 2.0 * g_y_xz_x_yy[i] * a_exp;

        g_y_0_0_0_0_xz_x_yz[i] = 2.0 * g_y_xz_x_yz[i] * a_exp;

        g_y_0_0_0_0_xz_x_zz[i] = 2.0 * g_y_xz_x_zz[i] * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_y_0_0_0_0_xz_y_xx, g_y_0_0_0_0_xz_y_xy, g_y_0_0_0_0_xz_y_xz, g_y_0_0_0_0_xz_y_yy, g_y_0_0_0_0_xz_y_yz, g_y_0_0_0_0_xz_y_zz, g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xz_y_xx[i] = 2.0 * g_y_xz_y_xx[i] * a_exp;

        g_y_0_0_0_0_xz_y_xy[i] = 2.0 * g_y_xz_y_xy[i] * a_exp;

        g_y_0_0_0_0_xz_y_xz[i] = 2.0 * g_y_xz_y_xz[i] * a_exp;

        g_y_0_0_0_0_xz_y_yy[i] = 2.0 * g_y_xz_y_yy[i] * a_exp;

        g_y_0_0_0_0_xz_y_yz[i] = 2.0 * g_y_xz_y_yz[i] * a_exp;

        g_y_0_0_0_0_xz_y_zz[i] = 2.0 * g_y_xz_y_zz[i] * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_y_0_0_0_0_xz_z_xx, g_y_0_0_0_0_xz_z_xy, g_y_0_0_0_0_xz_z_xz, g_y_0_0_0_0_xz_z_yy, g_y_0_0_0_0_xz_z_yz, g_y_0_0_0_0_xz_z_zz, g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xz_z_xx[i] = 2.0 * g_y_xz_z_xx[i] * a_exp;

        g_y_0_0_0_0_xz_z_xy[i] = 2.0 * g_y_xz_z_xy[i] * a_exp;

        g_y_0_0_0_0_xz_z_xz[i] = 2.0 * g_y_xz_z_xz[i] * a_exp;

        g_y_0_0_0_0_xz_z_yy[i] = 2.0 * g_y_xz_z_yy[i] * a_exp;

        g_y_0_0_0_0_xz_z_yz[i] = 2.0 * g_y_xz_z_yz[i] * a_exp;

        g_y_0_0_0_0_xz_z_zz[i] = 2.0 * g_y_xz_z_zz[i] * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_y_0_0_0_0_yy_x_xx, g_y_0_0_0_0_yy_x_xy, g_y_0_0_0_0_yy_x_xz, g_y_0_0_0_0_yy_x_yy, g_y_0_0_0_0_yy_x_yz, g_y_0_0_0_0_yy_x_zz, g_y_yy_x_xx, g_y_yy_x_xy, g_y_yy_x_xz, g_y_yy_x_yy, g_y_yy_x_yz, g_y_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yy_x_xx[i] = 2.0 * g_y_yy_x_xx[i] * a_exp;

        g_y_0_0_0_0_yy_x_xy[i] = 2.0 * g_y_yy_x_xy[i] * a_exp;

        g_y_0_0_0_0_yy_x_xz[i] = 2.0 * g_y_yy_x_xz[i] * a_exp;

        g_y_0_0_0_0_yy_x_yy[i] = 2.0 * g_y_yy_x_yy[i] * a_exp;

        g_y_0_0_0_0_yy_x_yz[i] = 2.0 * g_y_yy_x_yz[i] * a_exp;

        g_y_0_0_0_0_yy_x_zz[i] = 2.0 * g_y_yy_x_zz[i] * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_y_0_0_0_0_yy_y_xx, g_y_0_0_0_0_yy_y_xy, g_y_0_0_0_0_yy_y_xz, g_y_0_0_0_0_yy_y_yy, g_y_0_0_0_0_yy_y_yz, g_y_0_0_0_0_yy_y_zz, g_y_yy_y_xx, g_y_yy_y_xy, g_y_yy_y_xz, g_y_yy_y_yy, g_y_yy_y_yz, g_y_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yy_y_xx[i] = 2.0 * g_y_yy_y_xx[i] * a_exp;

        g_y_0_0_0_0_yy_y_xy[i] = 2.0 * g_y_yy_y_xy[i] * a_exp;

        g_y_0_0_0_0_yy_y_xz[i] = 2.0 * g_y_yy_y_xz[i] * a_exp;

        g_y_0_0_0_0_yy_y_yy[i] = 2.0 * g_y_yy_y_yy[i] * a_exp;

        g_y_0_0_0_0_yy_y_yz[i] = 2.0 * g_y_yy_y_yz[i] * a_exp;

        g_y_0_0_0_0_yy_y_zz[i] = 2.0 * g_y_yy_y_zz[i] * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_y_0_0_0_0_yy_z_xx, g_y_0_0_0_0_yy_z_xy, g_y_0_0_0_0_yy_z_xz, g_y_0_0_0_0_yy_z_yy, g_y_0_0_0_0_yy_z_yz, g_y_0_0_0_0_yy_z_zz, g_y_yy_z_xx, g_y_yy_z_xy, g_y_yy_z_xz, g_y_yy_z_yy, g_y_yy_z_yz, g_y_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yy_z_xx[i] = 2.0 * g_y_yy_z_xx[i] * a_exp;

        g_y_0_0_0_0_yy_z_xy[i] = 2.0 * g_y_yy_z_xy[i] * a_exp;

        g_y_0_0_0_0_yy_z_xz[i] = 2.0 * g_y_yy_z_xz[i] * a_exp;

        g_y_0_0_0_0_yy_z_yy[i] = 2.0 * g_y_yy_z_yy[i] * a_exp;

        g_y_0_0_0_0_yy_z_yz[i] = 2.0 * g_y_yy_z_yz[i] * a_exp;

        g_y_0_0_0_0_yy_z_zz[i] = 2.0 * g_y_yy_z_zz[i] * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_y_0_0_0_0_yz_x_xx, g_y_0_0_0_0_yz_x_xy, g_y_0_0_0_0_yz_x_xz, g_y_0_0_0_0_yz_x_yy, g_y_0_0_0_0_yz_x_yz, g_y_0_0_0_0_yz_x_zz, g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yz_x_xx[i] = 2.0 * g_y_yz_x_xx[i] * a_exp;

        g_y_0_0_0_0_yz_x_xy[i] = 2.0 * g_y_yz_x_xy[i] * a_exp;

        g_y_0_0_0_0_yz_x_xz[i] = 2.0 * g_y_yz_x_xz[i] * a_exp;

        g_y_0_0_0_0_yz_x_yy[i] = 2.0 * g_y_yz_x_yy[i] * a_exp;

        g_y_0_0_0_0_yz_x_yz[i] = 2.0 * g_y_yz_x_yz[i] * a_exp;

        g_y_0_0_0_0_yz_x_zz[i] = 2.0 * g_y_yz_x_zz[i] * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_y_0_0_0_0_yz_y_xx, g_y_0_0_0_0_yz_y_xy, g_y_0_0_0_0_yz_y_xz, g_y_0_0_0_0_yz_y_yy, g_y_0_0_0_0_yz_y_yz, g_y_0_0_0_0_yz_y_zz, g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yz_y_xx[i] = 2.0 * g_y_yz_y_xx[i] * a_exp;

        g_y_0_0_0_0_yz_y_xy[i] = 2.0 * g_y_yz_y_xy[i] * a_exp;

        g_y_0_0_0_0_yz_y_xz[i] = 2.0 * g_y_yz_y_xz[i] * a_exp;

        g_y_0_0_0_0_yz_y_yy[i] = 2.0 * g_y_yz_y_yy[i] * a_exp;

        g_y_0_0_0_0_yz_y_yz[i] = 2.0 * g_y_yz_y_yz[i] * a_exp;

        g_y_0_0_0_0_yz_y_zz[i] = 2.0 * g_y_yz_y_zz[i] * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_y_0_0_0_0_yz_z_xx, g_y_0_0_0_0_yz_z_xy, g_y_0_0_0_0_yz_z_xz, g_y_0_0_0_0_yz_z_yy, g_y_0_0_0_0_yz_z_yz, g_y_0_0_0_0_yz_z_zz, g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yz_z_xx[i] = 2.0 * g_y_yz_z_xx[i] * a_exp;

        g_y_0_0_0_0_yz_z_xy[i] = 2.0 * g_y_yz_z_xy[i] * a_exp;

        g_y_0_0_0_0_yz_z_xz[i] = 2.0 * g_y_yz_z_xz[i] * a_exp;

        g_y_0_0_0_0_yz_z_yy[i] = 2.0 * g_y_yz_z_yy[i] * a_exp;

        g_y_0_0_0_0_yz_z_yz[i] = 2.0 * g_y_yz_z_yz[i] * a_exp;

        g_y_0_0_0_0_yz_z_zz[i] = 2.0 * g_y_yz_z_zz[i] * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_y_0_0_0_0_zz_x_xx, g_y_0_0_0_0_zz_x_xy, g_y_0_0_0_0_zz_x_xz, g_y_0_0_0_0_zz_x_yy, g_y_0_0_0_0_zz_x_yz, g_y_0_0_0_0_zz_x_zz, g_y_zz_x_xx, g_y_zz_x_xy, g_y_zz_x_xz, g_y_zz_x_yy, g_y_zz_x_yz, g_y_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_zz_x_xx[i] = 2.0 * g_y_zz_x_xx[i] * a_exp;

        g_y_0_0_0_0_zz_x_xy[i] = 2.0 * g_y_zz_x_xy[i] * a_exp;

        g_y_0_0_0_0_zz_x_xz[i] = 2.0 * g_y_zz_x_xz[i] * a_exp;

        g_y_0_0_0_0_zz_x_yy[i] = 2.0 * g_y_zz_x_yy[i] * a_exp;

        g_y_0_0_0_0_zz_x_yz[i] = 2.0 * g_y_zz_x_yz[i] * a_exp;

        g_y_0_0_0_0_zz_x_zz[i] = 2.0 * g_y_zz_x_zz[i] * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_y_0_0_0_0_zz_y_xx, g_y_0_0_0_0_zz_y_xy, g_y_0_0_0_0_zz_y_xz, g_y_0_0_0_0_zz_y_yy, g_y_0_0_0_0_zz_y_yz, g_y_0_0_0_0_zz_y_zz, g_y_zz_y_xx, g_y_zz_y_xy, g_y_zz_y_xz, g_y_zz_y_yy, g_y_zz_y_yz, g_y_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_zz_y_xx[i] = 2.0 * g_y_zz_y_xx[i] * a_exp;

        g_y_0_0_0_0_zz_y_xy[i] = 2.0 * g_y_zz_y_xy[i] * a_exp;

        g_y_0_0_0_0_zz_y_xz[i] = 2.0 * g_y_zz_y_xz[i] * a_exp;

        g_y_0_0_0_0_zz_y_yy[i] = 2.0 * g_y_zz_y_yy[i] * a_exp;

        g_y_0_0_0_0_zz_y_yz[i] = 2.0 * g_y_zz_y_yz[i] * a_exp;

        g_y_0_0_0_0_zz_y_zz[i] = 2.0 * g_y_zz_y_zz[i] * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_y_0_0_0_0_zz_z_xx, g_y_0_0_0_0_zz_z_xy, g_y_0_0_0_0_zz_z_xz, g_y_0_0_0_0_zz_z_yy, g_y_0_0_0_0_zz_z_yz, g_y_0_0_0_0_zz_z_zz, g_y_zz_z_xx, g_y_zz_z_xy, g_y_zz_z_xz, g_y_zz_z_yy, g_y_zz_z_yz, g_y_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_zz_z_xx[i] = 2.0 * g_y_zz_z_xx[i] * a_exp;

        g_y_0_0_0_0_zz_z_xy[i] = 2.0 * g_y_zz_z_xy[i] * a_exp;

        g_y_0_0_0_0_zz_z_xz[i] = 2.0 * g_y_zz_z_xz[i] * a_exp;

        g_y_0_0_0_0_zz_z_yy[i] = 2.0 * g_y_zz_z_yy[i] * a_exp;

        g_y_0_0_0_0_zz_z_yz[i] = 2.0 * g_y_zz_z_yz[i] * a_exp;

        g_y_0_0_0_0_zz_z_zz[i] = 2.0 * g_y_zz_z_zz[i] * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_z_0_0_0_0_xx_x_xx, g_z_0_0_0_0_xx_x_xy, g_z_0_0_0_0_xx_x_xz, g_z_0_0_0_0_xx_x_yy, g_z_0_0_0_0_xx_x_yz, g_z_0_0_0_0_xx_x_zz, g_z_xx_x_xx, g_z_xx_x_xy, g_z_xx_x_xz, g_z_xx_x_yy, g_z_xx_x_yz, g_z_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xx_x_xx[i] = 2.0 * g_z_xx_x_xx[i] * a_exp;

        g_z_0_0_0_0_xx_x_xy[i] = 2.0 * g_z_xx_x_xy[i] * a_exp;

        g_z_0_0_0_0_xx_x_xz[i] = 2.0 * g_z_xx_x_xz[i] * a_exp;

        g_z_0_0_0_0_xx_x_yy[i] = 2.0 * g_z_xx_x_yy[i] * a_exp;

        g_z_0_0_0_0_xx_x_yz[i] = 2.0 * g_z_xx_x_yz[i] * a_exp;

        g_z_0_0_0_0_xx_x_zz[i] = 2.0 * g_z_xx_x_zz[i] * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_z_0_0_0_0_xx_y_xx, g_z_0_0_0_0_xx_y_xy, g_z_0_0_0_0_xx_y_xz, g_z_0_0_0_0_xx_y_yy, g_z_0_0_0_0_xx_y_yz, g_z_0_0_0_0_xx_y_zz, g_z_xx_y_xx, g_z_xx_y_xy, g_z_xx_y_xz, g_z_xx_y_yy, g_z_xx_y_yz, g_z_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xx_y_xx[i] = 2.0 * g_z_xx_y_xx[i] * a_exp;

        g_z_0_0_0_0_xx_y_xy[i] = 2.0 * g_z_xx_y_xy[i] * a_exp;

        g_z_0_0_0_0_xx_y_xz[i] = 2.0 * g_z_xx_y_xz[i] * a_exp;

        g_z_0_0_0_0_xx_y_yy[i] = 2.0 * g_z_xx_y_yy[i] * a_exp;

        g_z_0_0_0_0_xx_y_yz[i] = 2.0 * g_z_xx_y_yz[i] * a_exp;

        g_z_0_0_0_0_xx_y_zz[i] = 2.0 * g_z_xx_y_zz[i] * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_z_0_0_0_0_xx_z_xx, g_z_0_0_0_0_xx_z_xy, g_z_0_0_0_0_xx_z_xz, g_z_0_0_0_0_xx_z_yy, g_z_0_0_0_0_xx_z_yz, g_z_0_0_0_0_xx_z_zz, g_z_xx_z_xx, g_z_xx_z_xy, g_z_xx_z_xz, g_z_xx_z_yy, g_z_xx_z_yz, g_z_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xx_z_xx[i] = 2.0 * g_z_xx_z_xx[i] * a_exp;

        g_z_0_0_0_0_xx_z_xy[i] = 2.0 * g_z_xx_z_xy[i] * a_exp;

        g_z_0_0_0_0_xx_z_xz[i] = 2.0 * g_z_xx_z_xz[i] * a_exp;

        g_z_0_0_0_0_xx_z_yy[i] = 2.0 * g_z_xx_z_yy[i] * a_exp;

        g_z_0_0_0_0_xx_z_yz[i] = 2.0 * g_z_xx_z_yz[i] * a_exp;

        g_z_0_0_0_0_xx_z_zz[i] = 2.0 * g_z_xx_z_zz[i] * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_z_0_0_0_0_xy_x_xx, g_z_0_0_0_0_xy_x_xy, g_z_0_0_0_0_xy_x_xz, g_z_0_0_0_0_xy_x_yy, g_z_0_0_0_0_xy_x_yz, g_z_0_0_0_0_xy_x_zz, g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xy_x_xx[i] = 2.0 * g_z_xy_x_xx[i] * a_exp;

        g_z_0_0_0_0_xy_x_xy[i] = 2.0 * g_z_xy_x_xy[i] * a_exp;

        g_z_0_0_0_0_xy_x_xz[i] = 2.0 * g_z_xy_x_xz[i] * a_exp;

        g_z_0_0_0_0_xy_x_yy[i] = 2.0 * g_z_xy_x_yy[i] * a_exp;

        g_z_0_0_0_0_xy_x_yz[i] = 2.0 * g_z_xy_x_yz[i] * a_exp;

        g_z_0_0_0_0_xy_x_zz[i] = 2.0 * g_z_xy_x_zz[i] * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_z_0_0_0_0_xy_y_xx, g_z_0_0_0_0_xy_y_xy, g_z_0_0_0_0_xy_y_xz, g_z_0_0_0_0_xy_y_yy, g_z_0_0_0_0_xy_y_yz, g_z_0_0_0_0_xy_y_zz, g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xy_y_xx[i] = 2.0 * g_z_xy_y_xx[i] * a_exp;

        g_z_0_0_0_0_xy_y_xy[i] = 2.0 * g_z_xy_y_xy[i] * a_exp;

        g_z_0_0_0_0_xy_y_xz[i] = 2.0 * g_z_xy_y_xz[i] * a_exp;

        g_z_0_0_0_0_xy_y_yy[i] = 2.0 * g_z_xy_y_yy[i] * a_exp;

        g_z_0_0_0_0_xy_y_yz[i] = 2.0 * g_z_xy_y_yz[i] * a_exp;

        g_z_0_0_0_0_xy_y_zz[i] = 2.0 * g_z_xy_y_zz[i] * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_z_0_0_0_0_xy_z_xx, g_z_0_0_0_0_xy_z_xy, g_z_0_0_0_0_xy_z_xz, g_z_0_0_0_0_xy_z_yy, g_z_0_0_0_0_xy_z_yz, g_z_0_0_0_0_xy_z_zz, g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xy_z_xx[i] = 2.0 * g_z_xy_z_xx[i] * a_exp;

        g_z_0_0_0_0_xy_z_xy[i] = 2.0 * g_z_xy_z_xy[i] * a_exp;

        g_z_0_0_0_0_xy_z_xz[i] = 2.0 * g_z_xy_z_xz[i] * a_exp;

        g_z_0_0_0_0_xy_z_yy[i] = 2.0 * g_z_xy_z_yy[i] * a_exp;

        g_z_0_0_0_0_xy_z_yz[i] = 2.0 * g_z_xy_z_yz[i] * a_exp;

        g_z_0_0_0_0_xy_z_zz[i] = 2.0 * g_z_xy_z_zz[i] * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_z_0_0_0_0_xz_x_xx, g_z_0_0_0_0_xz_x_xy, g_z_0_0_0_0_xz_x_xz, g_z_0_0_0_0_xz_x_yy, g_z_0_0_0_0_xz_x_yz, g_z_0_0_0_0_xz_x_zz, g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xz_x_xx[i] = 2.0 * g_z_xz_x_xx[i] * a_exp;

        g_z_0_0_0_0_xz_x_xy[i] = 2.0 * g_z_xz_x_xy[i] * a_exp;

        g_z_0_0_0_0_xz_x_xz[i] = 2.0 * g_z_xz_x_xz[i] * a_exp;

        g_z_0_0_0_0_xz_x_yy[i] = 2.0 * g_z_xz_x_yy[i] * a_exp;

        g_z_0_0_0_0_xz_x_yz[i] = 2.0 * g_z_xz_x_yz[i] * a_exp;

        g_z_0_0_0_0_xz_x_zz[i] = 2.0 * g_z_xz_x_zz[i] * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_z_0_0_0_0_xz_y_xx, g_z_0_0_0_0_xz_y_xy, g_z_0_0_0_0_xz_y_xz, g_z_0_0_0_0_xz_y_yy, g_z_0_0_0_0_xz_y_yz, g_z_0_0_0_0_xz_y_zz, g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xz_y_xx[i] = 2.0 * g_z_xz_y_xx[i] * a_exp;

        g_z_0_0_0_0_xz_y_xy[i] = 2.0 * g_z_xz_y_xy[i] * a_exp;

        g_z_0_0_0_0_xz_y_xz[i] = 2.0 * g_z_xz_y_xz[i] * a_exp;

        g_z_0_0_0_0_xz_y_yy[i] = 2.0 * g_z_xz_y_yy[i] * a_exp;

        g_z_0_0_0_0_xz_y_yz[i] = 2.0 * g_z_xz_y_yz[i] * a_exp;

        g_z_0_0_0_0_xz_y_zz[i] = 2.0 * g_z_xz_y_zz[i] * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_z_0_0_0_0_xz_z_xx, g_z_0_0_0_0_xz_z_xy, g_z_0_0_0_0_xz_z_xz, g_z_0_0_0_0_xz_z_yy, g_z_0_0_0_0_xz_z_yz, g_z_0_0_0_0_xz_z_zz, g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xz_z_xx[i] = 2.0 * g_z_xz_z_xx[i] * a_exp;

        g_z_0_0_0_0_xz_z_xy[i] = 2.0 * g_z_xz_z_xy[i] * a_exp;

        g_z_0_0_0_0_xz_z_xz[i] = 2.0 * g_z_xz_z_xz[i] * a_exp;

        g_z_0_0_0_0_xz_z_yy[i] = 2.0 * g_z_xz_z_yy[i] * a_exp;

        g_z_0_0_0_0_xz_z_yz[i] = 2.0 * g_z_xz_z_yz[i] * a_exp;

        g_z_0_0_0_0_xz_z_zz[i] = 2.0 * g_z_xz_z_zz[i] * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_z_0_0_0_0_yy_x_xx, g_z_0_0_0_0_yy_x_xy, g_z_0_0_0_0_yy_x_xz, g_z_0_0_0_0_yy_x_yy, g_z_0_0_0_0_yy_x_yz, g_z_0_0_0_0_yy_x_zz, g_z_yy_x_xx, g_z_yy_x_xy, g_z_yy_x_xz, g_z_yy_x_yy, g_z_yy_x_yz, g_z_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yy_x_xx[i] = 2.0 * g_z_yy_x_xx[i] * a_exp;

        g_z_0_0_0_0_yy_x_xy[i] = 2.0 * g_z_yy_x_xy[i] * a_exp;

        g_z_0_0_0_0_yy_x_xz[i] = 2.0 * g_z_yy_x_xz[i] * a_exp;

        g_z_0_0_0_0_yy_x_yy[i] = 2.0 * g_z_yy_x_yy[i] * a_exp;

        g_z_0_0_0_0_yy_x_yz[i] = 2.0 * g_z_yy_x_yz[i] * a_exp;

        g_z_0_0_0_0_yy_x_zz[i] = 2.0 * g_z_yy_x_zz[i] * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_z_0_0_0_0_yy_y_xx, g_z_0_0_0_0_yy_y_xy, g_z_0_0_0_0_yy_y_xz, g_z_0_0_0_0_yy_y_yy, g_z_0_0_0_0_yy_y_yz, g_z_0_0_0_0_yy_y_zz, g_z_yy_y_xx, g_z_yy_y_xy, g_z_yy_y_xz, g_z_yy_y_yy, g_z_yy_y_yz, g_z_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yy_y_xx[i] = 2.0 * g_z_yy_y_xx[i] * a_exp;

        g_z_0_0_0_0_yy_y_xy[i] = 2.0 * g_z_yy_y_xy[i] * a_exp;

        g_z_0_0_0_0_yy_y_xz[i] = 2.0 * g_z_yy_y_xz[i] * a_exp;

        g_z_0_0_0_0_yy_y_yy[i] = 2.0 * g_z_yy_y_yy[i] * a_exp;

        g_z_0_0_0_0_yy_y_yz[i] = 2.0 * g_z_yy_y_yz[i] * a_exp;

        g_z_0_0_0_0_yy_y_zz[i] = 2.0 * g_z_yy_y_zz[i] * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_z_0_0_0_0_yy_z_xx, g_z_0_0_0_0_yy_z_xy, g_z_0_0_0_0_yy_z_xz, g_z_0_0_0_0_yy_z_yy, g_z_0_0_0_0_yy_z_yz, g_z_0_0_0_0_yy_z_zz, g_z_yy_z_xx, g_z_yy_z_xy, g_z_yy_z_xz, g_z_yy_z_yy, g_z_yy_z_yz, g_z_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yy_z_xx[i] = 2.0 * g_z_yy_z_xx[i] * a_exp;

        g_z_0_0_0_0_yy_z_xy[i] = 2.0 * g_z_yy_z_xy[i] * a_exp;

        g_z_0_0_0_0_yy_z_xz[i] = 2.0 * g_z_yy_z_xz[i] * a_exp;

        g_z_0_0_0_0_yy_z_yy[i] = 2.0 * g_z_yy_z_yy[i] * a_exp;

        g_z_0_0_0_0_yy_z_yz[i] = 2.0 * g_z_yy_z_yz[i] * a_exp;

        g_z_0_0_0_0_yy_z_zz[i] = 2.0 * g_z_yy_z_zz[i] * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_z_0_0_0_0_yz_x_xx, g_z_0_0_0_0_yz_x_xy, g_z_0_0_0_0_yz_x_xz, g_z_0_0_0_0_yz_x_yy, g_z_0_0_0_0_yz_x_yz, g_z_0_0_0_0_yz_x_zz, g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yz_x_xx[i] = 2.0 * g_z_yz_x_xx[i] * a_exp;

        g_z_0_0_0_0_yz_x_xy[i] = 2.0 * g_z_yz_x_xy[i] * a_exp;

        g_z_0_0_0_0_yz_x_xz[i] = 2.0 * g_z_yz_x_xz[i] * a_exp;

        g_z_0_0_0_0_yz_x_yy[i] = 2.0 * g_z_yz_x_yy[i] * a_exp;

        g_z_0_0_0_0_yz_x_yz[i] = 2.0 * g_z_yz_x_yz[i] * a_exp;

        g_z_0_0_0_0_yz_x_zz[i] = 2.0 * g_z_yz_x_zz[i] * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_z_0_0_0_0_yz_y_xx, g_z_0_0_0_0_yz_y_xy, g_z_0_0_0_0_yz_y_xz, g_z_0_0_0_0_yz_y_yy, g_z_0_0_0_0_yz_y_yz, g_z_0_0_0_0_yz_y_zz, g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yz_y_xx[i] = 2.0 * g_z_yz_y_xx[i] * a_exp;

        g_z_0_0_0_0_yz_y_xy[i] = 2.0 * g_z_yz_y_xy[i] * a_exp;

        g_z_0_0_0_0_yz_y_xz[i] = 2.0 * g_z_yz_y_xz[i] * a_exp;

        g_z_0_0_0_0_yz_y_yy[i] = 2.0 * g_z_yz_y_yy[i] * a_exp;

        g_z_0_0_0_0_yz_y_yz[i] = 2.0 * g_z_yz_y_yz[i] * a_exp;

        g_z_0_0_0_0_yz_y_zz[i] = 2.0 * g_z_yz_y_zz[i] * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_z_0_0_0_0_yz_z_xx, g_z_0_0_0_0_yz_z_xy, g_z_0_0_0_0_yz_z_xz, g_z_0_0_0_0_yz_z_yy, g_z_0_0_0_0_yz_z_yz, g_z_0_0_0_0_yz_z_zz, g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yz_z_xx[i] = 2.0 * g_z_yz_z_xx[i] * a_exp;

        g_z_0_0_0_0_yz_z_xy[i] = 2.0 * g_z_yz_z_xy[i] * a_exp;

        g_z_0_0_0_0_yz_z_xz[i] = 2.0 * g_z_yz_z_xz[i] * a_exp;

        g_z_0_0_0_0_yz_z_yy[i] = 2.0 * g_z_yz_z_yy[i] * a_exp;

        g_z_0_0_0_0_yz_z_yz[i] = 2.0 * g_z_yz_z_yz[i] * a_exp;

        g_z_0_0_0_0_yz_z_zz[i] = 2.0 * g_z_yz_z_zz[i] * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_z_0_0_0_0_zz_x_xx, g_z_0_0_0_0_zz_x_xy, g_z_0_0_0_0_zz_x_xz, g_z_0_0_0_0_zz_x_yy, g_z_0_0_0_0_zz_x_yz, g_z_0_0_0_0_zz_x_zz, g_z_zz_x_xx, g_z_zz_x_xy, g_z_zz_x_xz, g_z_zz_x_yy, g_z_zz_x_yz, g_z_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_zz_x_xx[i] = 2.0 * g_z_zz_x_xx[i] * a_exp;

        g_z_0_0_0_0_zz_x_xy[i] = 2.0 * g_z_zz_x_xy[i] * a_exp;

        g_z_0_0_0_0_zz_x_xz[i] = 2.0 * g_z_zz_x_xz[i] * a_exp;

        g_z_0_0_0_0_zz_x_yy[i] = 2.0 * g_z_zz_x_yy[i] * a_exp;

        g_z_0_0_0_0_zz_x_yz[i] = 2.0 * g_z_zz_x_yz[i] * a_exp;

        g_z_0_0_0_0_zz_x_zz[i] = 2.0 * g_z_zz_x_zz[i] * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_z_0_0_0_0_zz_y_xx, g_z_0_0_0_0_zz_y_xy, g_z_0_0_0_0_zz_y_xz, g_z_0_0_0_0_zz_y_yy, g_z_0_0_0_0_zz_y_yz, g_z_0_0_0_0_zz_y_zz, g_z_zz_y_xx, g_z_zz_y_xy, g_z_zz_y_xz, g_z_zz_y_yy, g_z_zz_y_yz, g_z_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_zz_y_xx[i] = 2.0 * g_z_zz_y_xx[i] * a_exp;

        g_z_0_0_0_0_zz_y_xy[i] = 2.0 * g_z_zz_y_xy[i] * a_exp;

        g_z_0_0_0_0_zz_y_xz[i] = 2.0 * g_z_zz_y_xz[i] * a_exp;

        g_z_0_0_0_0_zz_y_yy[i] = 2.0 * g_z_zz_y_yy[i] * a_exp;

        g_z_0_0_0_0_zz_y_yz[i] = 2.0 * g_z_zz_y_yz[i] * a_exp;

        g_z_0_0_0_0_zz_y_zz[i] = 2.0 * g_z_zz_y_zz[i] * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_z_0_0_0_0_zz_z_xx, g_z_0_0_0_0_zz_z_xy, g_z_0_0_0_0_zz_z_xz, g_z_0_0_0_0_zz_z_yy, g_z_0_0_0_0_zz_z_yz, g_z_0_0_0_0_zz_z_zz, g_z_zz_z_xx, g_z_zz_z_xy, g_z_zz_z_xz, g_z_zz_z_yy, g_z_zz_z_yz, g_z_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_zz_z_xx[i] = 2.0 * g_z_zz_z_xx[i] * a_exp;

        g_z_0_0_0_0_zz_z_xy[i] = 2.0 * g_z_zz_z_xy[i] * a_exp;

        g_z_0_0_0_0_zz_z_xz[i] = 2.0 * g_z_zz_z_xz[i] * a_exp;

        g_z_0_0_0_0_zz_z_yy[i] = 2.0 * g_z_zz_z_yy[i] * a_exp;

        g_z_0_0_0_0_zz_z_yz[i] = 2.0 * g_z_zz_z_yz[i] * a_exp;

        g_z_0_0_0_0_zz_z_zz[i] = 2.0 * g_z_zz_z_zz[i] * a_exp;
    }
}

} // t4c_geom namespace

