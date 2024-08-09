#include "GeomDeriv2000OfScalarForSPPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_sppd_0(CSimdArray<double>& buffer_2000_sppd,
                     const CSimdArray<double>& buffer_sppd,
                     const CSimdArray<double>& buffer_dppd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_sppd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sppd

    auto g_0_x_x_xx = buffer_sppd[0];

    auto g_0_x_x_xy = buffer_sppd[1];

    auto g_0_x_x_xz = buffer_sppd[2];

    auto g_0_x_x_yy = buffer_sppd[3];

    auto g_0_x_x_yz = buffer_sppd[4];

    auto g_0_x_x_zz = buffer_sppd[5];

    auto g_0_x_y_xx = buffer_sppd[6];

    auto g_0_x_y_xy = buffer_sppd[7];

    auto g_0_x_y_xz = buffer_sppd[8];

    auto g_0_x_y_yy = buffer_sppd[9];

    auto g_0_x_y_yz = buffer_sppd[10];

    auto g_0_x_y_zz = buffer_sppd[11];

    auto g_0_x_z_xx = buffer_sppd[12];

    auto g_0_x_z_xy = buffer_sppd[13];

    auto g_0_x_z_xz = buffer_sppd[14];

    auto g_0_x_z_yy = buffer_sppd[15];

    auto g_0_x_z_yz = buffer_sppd[16];

    auto g_0_x_z_zz = buffer_sppd[17];

    auto g_0_y_x_xx = buffer_sppd[18];

    auto g_0_y_x_xy = buffer_sppd[19];

    auto g_0_y_x_xz = buffer_sppd[20];

    auto g_0_y_x_yy = buffer_sppd[21];

    auto g_0_y_x_yz = buffer_sppd[22];

    auto g_0_y_x_zz = buffer_sppd[23];

    auto g_0_y_y_xx = buffer_sppd[24];

    auto g_0_y_y_xy = buffer_sppd[25];

    auto g_0_y_y_xz = buffer_sppd[26];

    auto g_0_y_y_yy = buffer_sppd[27];

    auto g_0_y_y_yz = buffer_sppd[28];

    auto g_0_y_y_zz = buffer_sppd[29];

    auto g_0_y_z_xx = buffer_sppd[30];

    auto g_0_y_z_xy = buffer_sppd[31];

    auto g_0_y_z_xz = buffer_sppd[32];

    auto g_0_y_z_yy = buffer_sppd[33];

    auto g_0_y_z_yz = buffer_sppd[34];

    auto g_0_y_z_zz = buffer_sppd[35];

    auto g_0_z_x_xx = buffer_sppd[36];

    auto g_0_z_x_xy = buffer_sppd[37];

    auto g_0_z_x_xz = buffer_sppd[38];

    auto g_0_z_x_yy = buffer_sppd[39];

    auto g_0_z_x_yz = buffer_sppd[40];

    auto g_0_z_x_zz = buffer_sppd[41];

    auto g_0_z_y_xx = buffer_sppd[42];

    auto g_0_z_y_xy = buffer_sppd[43];

    auto g_0_z_y_xz = buffer_sppd[44];

    auto g_0_z_y_yy = buffer_sppd[45];

    auto g_0_z_y_yz = buffer_sppd[46];

    auto g_0_z_y_zz = buffer_sppd[47];

    auto g_0_z_z_xx = buffer_sppd[48];

    auto g_0_z_z_xy = buffer_sppd[49];

    auto g_0_z_z_xz = buffer_sppd[50];

    auto g_0_z_z_yy = buffer_sppd[51];

    auto g_0_z_z_yz = buffer_sppd[52];

    auto g_0_z_z_zz = buffer_sppd[53];

    /// Set up components of auxilary buffer : buffer_dppd

    auto g_xx_x_x_xx = buffer_dppd[0];

    auto g_xx_x_x_xy = buffer_dppd[1];

    auto g_xx_x_x_xz = buffer_dppd[2];

    auto g_xx_x_x_yy = buffer_dppd[3];

    auto g_xx_x_x_yz = buffer_dppd[4];

    auto g_xx_x_x_zz = buffer_dppd[5];

    auto g_xx_x_y_xx = buffer_dppd[6];

    auto g_xx_x_y_xy = buffer_dppd[7];

    auto g_xx_x_y_xz = buffer_dppd[8];

    auto g_xx_x_y_yy = buffer_dppd[9];

    auto g_xx_x_y_yz = buffer_dppd[10];

    auto g_xx_x_y_zz = buffer_dppd[11];

    auto g_xx_x_z_xx = buffer_dppd[12];

    auto g_xx_x_z_xy = buffer_dppd[13];

    auto g_xx_x_z_xz = buffer_dppd[14];

    auto g_xx_x_z_yy = buffer_dppd[15];

    auto g_xx_x_z_yz = buffer_dppd[16];

    auto g_xx_x_z_zz = buffer_dppd[17];

    auto g_xx_y_x_xx = buffer_dppd[18];

    auto g_xx_y_x_xy = buffer_dppd[19];

    auto g_xx_y_x_xz = buffer_dppd[20];

    auto g_xx_y_x_yy = buffer_dppd[21];

    auto g_xx_y_x_yz = buffer_dppd[22];

    auto g_xx_y_x_zz = buffer_dppd[23];

    auto g_xx_y_y_xx = buffer_dppd[24];

    auto g_xx_y_y_xy = buffer_dppd[25];

    auto g_xx_y_y_xz = buffer_dppd[26];

    auto g_xx_y_y_yy = buffer_dppd[27];

    auto g_xx_y_y_yz = buffer_dppd[28];

    auto g_xx_y_y_zz = buffer_dppd[29];

    auto g_xx_y_z_xx = buffer_dppd[30];

    auto g_xx_y_z_xy = buffer_dppd[31];

    auto g_xx_y_z_xz = buffer_dppd[32];

    auto g_xx_y_z_yy = buffer_dppd[33];

    auto g_xx_y_z_yz = buffer_dppd[34];

    auto g_xx_y_z_zz = buffer_dppd[35];

    auto g_xx_z_x_xx = buffer_dppd[36];

    auto g_xx_z_x_xy = buffer_dppd[37];

    auto g_xx_z_x_xz = buffer_dppd[38];

    auto g_xx_z_x_yy = buffer_dppd[39];

    auto g_xx_z_x_yz = buffer_dppd[40];

    auto g_xx_z_x_zz = buffer_dppd[41];

    auto g_xx_z_y_xx = buffer_dppd[42];

    auto g_xx_z_y_xy = buffer_dppd[43];

    auto g_xx_z_y_xz = buffer_dppd[44];

    auto g_xx_z_y_yy = buffer_dppd[45];

    auto g_xx_z_y_yz = buffer_dppd[46];

    auto g_xx_z_y_zz = buffer_dppd[47];

    auto g_xx_z_z_xx = buffer_dppd[48];

    auto g_xx_z_z_xy = buffer_dppd[49];

    auto g_xx_z_z_xz = buffer_dppd[50];

    auto g_xx_z_z_yy = buffer_dppd[51];

    auto g_xx_z_z_yz = buffer_dppd[52];

    auto g_xx_z_z_zz = buffer_dppd[53];

    auto g_xy_x_x_xx = buffer_dppd[54];

    auto g_xy_x_x_xy = buffer_dppd[55];

    auto g_xy_x_x_xz = buffer_dppd[56];

    auto g_xy_x_x_yy = buffer_dppd[57];

    auto g_xy_x_x_yz = buffer_dppd[58];

    auto g_xy_x_x_zz = buffer_dppd[59];

    auto g_xy_x_y_xx = buffer_dppd[60];

    auto g_xy_x_y_xy = buffer_dppd[61];

    auto g_xy_x_y_xz = buffer_dppd[62];

    auto g_xy_x_y_yy = buffer_dppd[63];

    auto g_xy_x_y_yz = buffer_dppd[64];

    auto g_xy_x_y_zz = buffer_dppd[65];

    auto g_xy_x_z_xx = buffer_dppd[66];

    auto g_xy_x_z_xy = buffer_dppd[67];

    auto g_xy_x_z_xz = buffer_dppd[68];

    auto g_xy_x_z_yy = buffer_dppd[69];

    auto g_xy_x_z_yz = buffer_dppd[70];

    auto g_xy_x_z_zz = buffer_dppd[71];

    auto g_xy_y_x_xx = buffer_dppd[72];

    auto g_xy_y_x_xy = buffer_dppd[73];

    auto g_xy_y_x_xz = buffer_dppd[74];

    auto g_xy_y_x_yy = buffer_dppd[75];

    auto g_xy_y_x_yz = buffer_dppd[76];

    auto g_xy_y_x_zz = buffer_dppd[77];

    auto g_xy_y_y_xx = buffer_dppd[78];

    auto g_xy_y_y_xy = buffer_dppd[79];

    auto g_xy_y_y_xz = buffer_dppd[80];

    auto g_xy_y_y_yy = buffer_dppd[81];

    auto g_xy_y_y_yz = buffer_dppd[82];

    auto g_xy_y_y_zz = buffer_dppd[83];

    auto g_xy_y_z_xx = buffer_dppd[84];

    auto g_xy_y_z_xy = buffer_dppd[85];

    auto g_xy_y_z_xz = buffer_dppd[86];

    auto g_xy_y_z_yy = buffer_dppd[87];

    auto g_xy_y_z_yz = buffer_dppd[88];

    auto g_xy_y_z_zz = buffer_dppd[89];

    auto g_xy_z_x_xx = buffer_dppd[90];

    auto g_xy_z_x_xy = buffer_dppd[91];

    auto g_xy_z_x_xz = buffer_dppd[92];

    auto g_xy_z_x_yy = buffer_dppd[93];

    auto g_xy_z_x_yz = buffer_dppd[94];

    auto g_xy_z_x_zz = buffer_dppd[95];

    auto g_xy_z_y_xx = buffer_dppd[96];

    auto g_xy_z_y_xy = buffer_dppd[97];

    auto g_xy_z_y_xz = buffer_dppd[98];

    auto g_xy_z_y_yy = buffer_dppd[99];

    auto g_xy_z_y_yz = buffer_dppd[100];

    auto g_xy_z_y_zz = buffer_dppd[101];

    auto g_xy_z_z_xx = buffer_dppd[102];

    auto g_xy_z_z_xy = buffer_dppd[103];

    auto g_xy_z_z_xz = buffer_dppd[104];

    auto g_xy_z_z_yy = buffer_dppd[105];

    auto g_xy_z_z_yz = buffer_dppd[106];

    auto g_xy_z_z_zz = buffer_dppd[107];

    auto g_xz_x_x_xx = buffer_dppd[108];

    auto g_xz_x_x_xy = buffer_dppd[109];

    auto g_xz_x_x_xz = buffer_dppd[110];

    auto g_xz_x_x_yy = buffer_dppd[111];

    auto g_xz_x_x_yz = buffer_dppd[112];

    auto g_xz_x_x_zz = buffer_dppd[113];

    auto g_xz_x_y_xx = buffer_dppd[114];

    auto g_xz_x_y_xy = buffer_dppd[115];

    auto g_xz_x_y_xz = buffer_dppd[116];

    auto g_xz_x_y_yy = buffer_dppd[117];

    auto g_xz_x_y_yz = buffer_dppd[118];

    auto g_xz_x_y_zz = buffer_dppd[119];

    auto g_xz_x_z_xx = buffer_dppd[120];

    auto g_xz_x_z_xy = buffer_dppd[121];

    auto g_xz_x_z_xz = buffer_dppd[122];

    auto g_xz_x_z_yy = buffer_dppd[123];

    auto g_xz_x_z_yz = buffer_dppd[124];

    auto g_xz_x_z_zz = buffer_dppd[125];

    auto g_xz_y_x_xx = buffer_dppd[126];

    auto g_xz_y_x_xy = buffer_dppd[127];

    auto g_xz_y_x_xz = buffer_dppd[128];

    auto g_xz_y_x_yy = buffer_dppd[129];

    auto g_xz_y_x_yz = buffer_dppd[130];

    auto g_xz_y_x_zz = buffer_dppd[131];

    auto g_xz_y_y_xx = buffer_dppd[132];

    auto g_xz_y_y_xy = buffer_dppd[133];

    auto g_xz_y_y_xz = buffer_dppd[134];

    auto g_xz_y_y_yy = buffer_dppd[135];

    auto g_xz_y_y_yz = buffer_dppd[136];

    auto g_xz_y_y_zz = buffer_dppd[137];

    auto g_xz_y_z_xx = buffer_dppd[138];

    auto g_xz_y_z_xy = buffer_dppd[139];

    auto g_xz_y_z_xz = buffer_dppd[140];

    auto g_xz_y_z_yy = buffer_dppd[141];

    auto g_xz_y_z_yz = buffer_dppd[142];

    auto g_xz_y_z_zz = buffer_dppd[143];

    auto g_xz_z_x_xx = buffer_dppd[144];

    auto g_xz_z_x_xy = buffer_dppd[145];

    auto g_xz_z_x_xz = buffer_dppd[146];

    auto g_xz_z_x_yy = buffer_dppd[147];

    auto g_xz_z_x_yz = buffer_dppd[148];

    auto g_xz_z_x_zz = buffer_dppd[149];

    auto g_xz_z_y_xx = buffer_dppd[150];

    auto g_xz_z_y_xy = buffer_dppd[151];

    auto g_xz_z_y_xz = buffer_dppd[152];

    auto g_xz_z_y_yy = buffer_dppd[153];

    auto g_xz_z_y_yz = buffer_dppd[154];

    auto g_xz_z_y_zz = buffer_dppd[155];

    auto g_xz_z_z_xx = buffer_dppd[156];

    auto g_xz_z_z_xy = buffer_dppd[157];

    auto g_xz_z_z_xz = buffer_dppd[158];

    auto g_xz_z_z_yy = buffer_dppd[159];

    auto g_xz_z_z_yz = buffer_dppd[160];

    auto g_xz_z_z_zz = buffer_dppd[161];

    auto g_yy_x_x_xx = buffer_dppd[162];

    auto g_yy_x_x_xy = buffer_dppd[163];

    auto g_yy_x_x_xz = buffer_dppd[164];

    auto g_yy_x_x_yy = buffer_dppd[165];

    auto g_yy_x_x_yz = buffer_dppd[166];

    auto g_yy_x_x_zz = buffer_dppd[167];

    auto g_yy_x_y_xx = buffer_dppd[168];

    auto g_yy_x_y_xy = buffer_dppd[169];

    auto g_yy_x_y_xz = buffer_dppd[170];

    auto g_yy_x_y_yy = buffer_dppd[171];

    auto g_yy_x_y_yz = buffer_dppd[172];

    auto g_yy_x_y_zz = buffer_dppd[173];

    auto g_yy_x_z_xx = buffer_dppd[174];

    auto g_yy_x_z_xy = buffer_dppd[175];

    auto g_yy_x_z_xz = buffer_dppd[176];

    auto g_yy_x_z_yy = buffer_dppd[177];

    auto g_yy_x_z_yz = buffer_dppd[178];

    auto g_yy_x_z_zz = buffer_dppd[179];

    auto g_yy_y_x_xx = buffer_dppd[180];

    auto g_yy_y_x_xy = buffer_dppd[181];

    auto g_yy_y_x_xz = buffer_dppd[182];

    auto g_yy_y_x_yy = buffer_dppd[183];

    auto g_yy_y_x_yz = buffer_dppd[184];

    auto g_yy_y_x_zz = buffer_dppd[185];

    auto g_yy_y_y_xx = buffer_dppd[186];

    auto g_yy_y_y_xy = buffer_dppd[187];

    auto g_yy_y_y_xz = buffer_dppd[188];

    auto g_yy_y_y_yy = buffer_dppd[189];

    auto g_yy_y_y_yz = buffer_dppd[190];

    auto g_yy_y_y_zz = buffer_dppd[191];

    auto g_yy_y_z_xx = buffer_dppd[192];

    auto g_yy_y_z_xy = buffer_dppd[193];

    auto g_yy_y_z_xz = buffer_dppd[194];

    auto g_yy_y_z_yy = buffer_dppd[195];

    auto g_yy_y_z_yz = buffer_dppd[196];

    auto g_yy_y_z_zz = buffer_dppd[197];

    auto g_yy_z_x_xx = buffer_dppd[198];

    auto g_yy_z_x_xy = buffer_dppd[199];

    auto g_yy_z_x_xz = buffer_dppd[200];

    auto g_yy_z_x_yy = buffer_dppd[201];

    auto g_yy_z_x_yz = buffer_dppd[202];

    auto g_yy_z_x_zz = buffer_dppd[203];

    auto g_yy_z_y_xx = buffer_dppd[204];

    auto g_yy_z_y_xy = buffer_dppd[205];

    auto g_yy_z_y_xz = buffer_dppd[206];

    auto g_yy_z_y_yy = buffer_dppd[207];

    auto g_yy_z_y_yz = buffer_dppd[208];

    auto g_yy_z_y_zz = buffer_dppd[209];

    auto g_yy_z_z_xx = buffer_dppd[210];

    auto g_yy_z_z_xy = buffer_dppd[211];

    auto g_yy_z_z_xz = buffer_dppd[212];

    auto g_yy_z_z_yy = buffer_dppd[213];

    auto g_yy_z_z_yz = buffer_dppd[214];

    auto g_yy_z_z_zz = buffer_dppd[215];

    auto g_yz_x_x_xx = buffer_dppd[216];

    auto g_yz_x_x_xy = buffer_dppd[217];

    auto g_yz_x_x_xz = buffer_dppd[218];

    auto g_yz_x_x_yy = buffer_dppd[219];

    auto g_yz_x_x_yz = buffer_dppd[220];

    auto g_yz_x_x_zz = buffer_dppd[221];

    auto g_yz_x_y_xx = buffer_dppd[222];

    auto g_yz_x_y_xy = buffer_dppd[223];

    auto g_yz_x_y_xz = buffer_dppd[224];

    auto g_yz_x_y_yy = buffer_dppd[225];

    auto g_yz_x_y_yz = buffer_dppd[226];

    auto g_yz_x_y_zz = buffer_dppd[227];

    auto g_yz_x_z_xx = buffer_dppd[228];

    auto g_yz_x_z_xy = buffer_dppd[229];

    auto g_yz_x_z_xz = buffer_dppd[230];

    auto g_yz_x_z_yy = buffer_dppd[231];

    auto g_yz_x_z_yz = buffer_dppd[232];

    auto g_yz_x_z_zz = buffer_dppd[233];

    auto g_yz_y_x_xx = buffer_dppd[234];

    auto g_yz_y_x_xy = buffer_dppd[235];

    auto g_yz_y_x_xz = buffer_dppd[236];

    auto g_yz_y_x_yy = buffer_dppd[237];

    auto g_yz_y_x_yz = buffer_dppd[238];

    auto g_yz_y_x_zz = buffer_dppd[239];

    auto g_yz_y_y_xx = buffer_dppd[240];

    auto g_yz_y_y_xy = buffer_dppd[241];

    auto g_yz_y_y_xz = buffer_dppd[242];

    auto g_yz_y_y_yy = buffer_dppd[243];

    auto g_yz_y_y_yz = buffer_dppd[244];

    auto g_yz_y_y_zz = buffer_dppd[245];

    auto g_yz_y_z_xx = buffer_dppd[246];

    auto g_yz_y_z_xy = buffer_dppd[247];

    auto g_yz_y_z_xz = buffer_dppd[248];

    auto g_yz_y_z_yy = buffer_dppd[249];

    auto g_yz_y_z_yz = buffer_dppd[250];

    auto g_yz_y_z_zz = buffer_dppd[251];

    auto g_yz_z_x_xx = buffer_dppd[252];

    auto g_yz_z_x_xy = buffer_dppd[253];

    auto g_yz_z_x_xz = buffer_dppd[254];

    auto g_yz_z_x_yy = buffer_dppd[255];

    auto g_yz_z_x_yz = buffer_dppd[256];

    auto g_yz_z_x_zz = buffer_dppd[257];

    auto g_yz_z_y_xx = buffer_dppd[258];

    auto g_yz_z_y_xy = buffer_dppd[259];

    auto g_yz_z_y_xz = buffer_dppd[260];

    auto g_yz_z_y_yy = buffer_dppd[261];

    auto g_yz_z_y_yz = buffer_dppd[262];

    auto g_yz_z_y_zz = buffer_dppd[263];

    auto g_yz_z_z_xx = buffer_dppd[264];

    auto g_yz_z_z_xy = buffer_dppd[265];

    auto g_yz_z_z_xz = buffer_dppd[266];

    auto g_yz_z_z_yy = buffer_dppd[267];

    auto g_yz_z_z_yz = buffer_dppd[268];

    auto g_yz_z_z_zz = buffer_dppd[269];

    auto g_zz_x_x_xx = buffer_dppd[270];

    auto g_zz_x_x_xy = buffer_dppd[271];

    auto g_zz_x_x_xz = buffer_dppd[272];

    auto g_zz_x_x_yy = buffer_dppd[273];

    auto g_zz_x_x_yz = buffer_dppd[274];

    auto g_zz_x_x_zz = buffer_dppd[275];

    auto g_zz_x_y_xx = buffer_dppd[276];

    auto g_zz_x_y_xy = buffer_dppd[277];

    auto g_zz_x_y_xz = buffer_dppd[278];

    auto g_zz_x_y_yy = buffer_dppd[279];

    auto g_zz_x_y_yz = buffer_dppd[280];

    auto g_zz_x_y_zz = buffer_dppd[281];

    auto g_zz_x_z_xx = buffer_dppd[282];

    auto g_zz_x_z_xy = buffer_dppd[283];

    auto g_zz_x_z_xz = buffer_dppd[284];

    auto g_zz_x_z_yy = buffer_dppd[285];

    auto g_zz_x_z_yz = buffer_dppd[286];

    auto g_zz_x_z_zz = buffer_dppd[287];

    auto g_zz_y_x_xx = buffer_dppd[288];

    auto g_zz_y_x_xy = buffer_dppd[289];

    auto g_zz_y_x_xz = buffer_dppd[290];

    auto g_zz_y_x_yy = buffer_dppd[291];

    auto g_zz_y_x_yz = buffer_dppd[292];

    auto g_zz_y_x_zz = buffer_dppd[293];

    auto g_zz_y_y_xx = buffer_dppd[294];

    auto g_zz_y_y_xy = buffer_dppd[295];

    auto g_zz_y_y_xz = buffer_dppd[296];

    auto g_zz_y_y_yy = buffer_dppd[297];

    auto g_zz_y_y_yz = buffer_dppd[298];

    auto g_zz_y_y_zz = buffer_dppd[299];

    auto g_zz_y_z_xx = buffer_dppd[300];

    auto g_zz_y_z_xy = buffer_dppd[301];

    auto g_zz_y_z_xz = buffer_dppd[302];

    auto g_zz_y_z_yy = buffer_dppd[303];

    auto g_zz_y_z_yz = buffer_dppd[304];

    auto g_zz_y_z_zz = buffer_dppd[305];

    auto g_zz_z_x_xx = buffer_dppd[306];

    auto g_zz_z_x_xy = buffer_dppd[307];

    auto g_zz_z_x_xz = buffer_dppd[308];

    auto g_zz_z_x_yy = buffer_dppd[309];

    auto g_zz_z_x_yz = buffer_dppd[310];

    auto g_zz_z_x_zz = buffer_dppd[311];

    auto g_zz_z_y_xx = buffer_dppd[312];

    auto g_zz_z_y_xy = buffer_dppd[313];

    auto g_zz_z_y_xz = buffer_dppd[314];

    auto g_zz_z_y_yy = buffer_dppd[315];

    auto g_zz_z_y_yz = buffer_dppd[316];

    auto g_zz_z_y_zz = buffer_dppd[317];

    auto g_zz_z_z_xx = buffer_dppd[318];

    auto g_zz_z_z_xy = buffer_dppd[319];

    auto g_zz_z_z_xz = buffer_dppd[320];

    auto g_zz_z_z_yy = buffer_dppd[321];

    auto g_zz_z_z_yz = buffer_dppd[322];

    auto g_zz_z_z_zz = buffer_dppd[323];

    /// Set up components of integrals buffer : buffer_2000_sppd

    auto g_xx_0_0_0_0_x_x_xx = buffer_2000_sppd[0];

    auto g_xx_0_0_0_0_x_x_xy = buffer_2000_sppd[1];

    auto g_xx_0_0_0_0_x_x_xz = buffer_2000_sppd[2];

    auto g_xx_0_0_0_0_x_x_yy = buffer_2000_sppd[3];

    auto g_xx_0_0_0_0_x_x_yz = buffer_2000_sppd[4];

    auto g_xx_0_0_0_0_x_x_zz = buffer_2000_sppd[5];

    auto g_xx_0_0_0_0_x_y_xx = buffer_2000_sppd[6];

    auto g_xx_0_0_0_0_x_y_xy = buffer_2000_sppd[7];

    auto g_xx_0_0_0_0_x_y_xz = buffer_2000_sppd[8];

    auto g_xx_0_0_0_0_x_y_yy = buffer_2000_sppd[9];

    auto g_xx_0_0_0_0_x_y_yz = buffer_2000_sppd[10];

    auto g_xx_0_0_0_0_x_y_zz = buffer_2000_sppd[11];

    auto g_xx_0_0_0_0_x_z_xx = buffer_2000_sppd[12];

    auto g_xx_0_0_0_0_x_z_xy = buffer_2000_sppd[13];

    auto g_xx_0_0_0_0_x_z_xz = buffer_2000_sppd[14];

    auto g_xx_0_0_0_0_x_z_yy = buffer_2000_sppd[15];

    auto g_xx_0_0_0_0_x_z_yz = buffer_2000_sppd[16];

    auto g_xx_0_0_0_0_x_z_zz = buffer_2000_sppd[17];

    auto g_xx_0_0_0_0_y_x_xx = buffer_2000_sppd[18];

    auto g_xx_0_0_0_0_y_x_xy = buffer_2000_sppd[19];

    auto g_xx_0_0_0_0_y_x_xz = buffer_2000_sppd[20];

    auto g_xx_0_0_0_0_y_x_yy = buffer_2000_sppd[21];

    auto g_xx_0_0_0_0_y_x_yz = buffer_2000_sppd[22];

    auto g_xx_0_0_0_0_y_x_zz = buffer_2000_sppd[23];

    auto g_xx_0_0_0_0_y_y_xx = buffer_2000_sppd[24];

    auto g_xx_0_0_0_0_y_y_xy = buffer_2000_sppd[25];

    auto g_xx_0_0_0_0_y_y_xz = buffer_2000_sppd[26];

    auto g_xx_0_0_0_0_y_y_yy = buffer_2000_sppd[27];

    auto g_xx_0_0_0_0_y_y_yz = buffer_2000_sppd[28];

    auto g_xx_0_0_0_0_y_y_zz = buffer_2000_sppd[29];

    auto g_xx_0_0_0_0_y_z_xx = buffer_2000_sppd[30];

    auto g_xx_0_0_0_0_y_z_xy = buffer_2000_sppd[31];

    auto g_xx_0_0_0_0_y_z_xz = buffer_2000_sppd[32];

    auto g_xx_0_0_0_0_y_z_yy = buffer_2000_sppd[33];

    auto g_xx_0_0_0_0_y_z_yz = buffer_2000_sppd[34];

    auto g_xx_0_0_0_0_y_z_zz = buffer_2000_sppd[35];

    auto g_xx_0_0_0_0_z_x_xx = buffer_2000_sppd[36];

    auto g_xx_0_0_0_0_z_x_xy = buffer_2000_sppd[37];

    auto g_xx_0_0_0_0_z_x_xz = buffer_2000_sppd[38];

    auto g_xx_0_0_0_0_z_x_yy = buffer_2000_sppd[39];

    auto g_xx_0_0_0_0_z_x_yz = buffer_2000_sppd[40];

    auto g_xx_0_0_0_0_z_x_zz = buffer_2000_sppd[41];

    auto g_xx_0_0_0_0_z_y_xx = buffer_2000_sppd[42];

    auto g_xx_0_0_0_0_z_y_xy = buffer_2000_sppd[43];

    auto g_xx_0_0_0_0_z_y_xz = buffer_2000_sppd[44];

    auto g_xx_0_0_0_0_z_y_yy = buffer_2000_sppd[45];

    auto g_xx_0_0_0_0_z_y_yz = buffer_2000_sppd[46];

    auto g_xx_0_0_0_0_z_y_zz = buffer_2000_sppd[47];

    auto g_xx_0_0_0_0_z_z_xx = buffer_2000_sppd[48];

    auto g_xx_0_0_0_0_z_z_xy = buffer_2000_sppd[49];

    auto g_xx_0_0_0_0_z_z_xz = buffer_2000_sppd[50];

    auto g_xx_0_0_0_0_z_z_yy = buffer_2000_sppd[51];

    auto g_xx_0_0_0_0_z_z_yz = buffer_2000_sppd[52];

    auto g_xx_0_0_0_0_z_z_zz = buffer_2000_sppd[53];

    auto g_xy_0_0_0_0_x_x_xx = buffer_2000_sppd[54];

    auto g_xy_0_0_0_0_x_x_xy = buffer_2000_sppd[55];

    auto g_xy_0_0_0_0_x_x_xz = buffer_2000_sppd[56];

    auto g_xy_0_0_0_0_x_x_yy = buffer_2000_sppd[57];

    auto g_xy_0_0_0_0_x_x_yz = buffer_2000_sppd[58];

    auto g_xy_0_0_0_0_x_x_zz = buffer_2000_sppd[59];

    auto g_xy_0_0_0_0_x_y_xx = buffer_2000_sppd[60];

    auto g_xy_0_0_0_0_x_y_xy = buffer_2000_sppd[61];

    auto g_xy_0_0_0_0_x_y_xz = buffer_2000_sppd[62];

    auto g_xy_0_0_0_0_x_y_yy = buffer_2000_sppd[63];

    auto g_xy_0_0_0_0_x_y_yz = buffer_2000_sppd[64];

    auto g_xy_0_0_0_0_x_y_zz = buffer_2000_sppd[65];

    auto g_xy_0_0_0_0_x_z_xx = buffer_2000_sppd[66];

    auto g_xy_0_0_0_0_x_z_xy = buffer_2000_sppd[67];

    auto g_xy_0_0_0_0_x_z_xz = buffer_2000_sppd[68];

    auto g_xy_0_0_0_0_x_z_yy = buffer_2000_sppd[69];

    auto g_xy_0_0_0_0_x_z_yz = buffer_2000_sppd[70];

    auto g_xy_0_0_0_0_x_z_zz = buffer_2000_sppd[71];

    auto g_xy_0_0_0_0_y_x_xx = buffer_2000_sppd[72];

    auto g_xy_0_0_0_0_y_x_xy = buffer_2000_sppd[73];

    auto g_xy_0_0_0_0_y_x_xz = buffer_2000_sppd[74];

    auto g_xy_0_0_0_0_y_x_yy = buffer_2000_sppd[75];

    auto g_xy_0_0_0_0_y_x_yz = buffer_2000_sppd[76];

    auto g_xy_0_0_0_0_y_x_zz = buffer_2000_sppd[77];

    auto g_xy_0_0_0_0_y_y_xx = buffer_2000_sppd[78];

    auto g_xy_0_0_0_0_y_y_xy = buffer_2000_sppd[79];

    auto g_xy_0_0_0_0_y_y_xz = buffer_2000_sppd[80];

    auto g_xy_0_0_0_0_y_y_yy = buffer_2000_sppd[81];

    auto g_xy_0_0_0_0_y_y_yz = buffer_2000_sppd[82];

    auto g_xy_0_0_0_0_y_y_zz = buffer_2000_sppd[83];

    auto g_xy_0_0_0_0_y_z_xx = buffer_2000_sppd[84];

    auto g_xy_0_0_0_0_y_z_xy = buffer_2000_sppd[85];

    auto g_xy_0_0_0_0_y_z_xz = buffer_2000_sppd[86];

    auto g_xy_0_0_0_0_y_z_yy = buffer_2000_sppd[87];

    auto g_xy_0_0_0_0_y_z_yz = buffer_2000_sppd[88];

    auto g_xy_0_0_0_0_y_z_zz = buffer_2000_sppd[89];

    auto g_xy_0_0_0_0_z_x_xx = buffer_2000_sppd[90];

    auto g_xy_0_0_0_0_z_x_xy = buffer_2000_sppd[91];

    auto g_xy_0_0_0_0_z_x_xz = buffer_2000_sppd[92];

    auto g_xy_0_0_0_0_z_x_yy = buffer_2000_sppd[93];

    auto g_xy_0_0_0_0_z_x_yz = buffer_2000_sppd[94];

    auto g_xy_0_0_0_0_z_x_zz = buffer_2000_sppd[95];

    auto g_xy_0_0_0_0_z_y_xx = buffer_2000_sppd[96];

    auto g_xy_0_0_0_0_z_y_xy = buffer_2000_sppd[97];

    auto g_xy_0_0_0_0_z_y_xz = buffer_2000_sppd[98];

    auto g_xy_0_0_0_0_z_y_yy = buffer_2000_sppd[99];

    auto g_xy_0_0_0_0_z_y_yz = buffer_2000_sppd[100];

    auto g_xy_0_0_0_0_z_y_zz = buffer_2000_sppd[101];

    auto g_xy_0_0_0_0_z_z_xx = buffer_2000_sppd[102];

    auto g_xy_0_0_0_0_z_z_xy = buffer_2000_sppd[103];

    auto g_xy_0_0_0_0_z_z_xz = buffer_2000_sppd[104];

    auto g_xy_0_0_0_0_z_z_yy = buffer_2000_sppd[105];

    auto g_xy_0_0_0_0_z_z_yz = buffer_2000_sppd[106];

    auto g_xy_0_0_0_0_z_z_zz = buffer_2000_sppd[107];

    auto g_xz_0_0_0_0_x_x_xx = buffer_2000_sppd[108];

    auto g_xz_0_0_0_0_x_x_xy = buffer_2000_sppd[109];

    auto g_xz_0_0_0_0_x_x_xz = buffer_2000_sppd[110];

    auto g_xz_0_0_0_0_x_x_yy = buffer_2000_sppd[111];

    auto g_xz_0_0_0_0_x_x_yz = buffer_2000_sppd[112];

    auto g_xz_0_0_0_0_x_x_zz = buffer_2000_sppd[113];

    auto g_xz_0_0_0_0_x_y_xx = buffer_2000_sppd[114];

    auto g_xz_0_0_0_0_x_y_xy = buffer_2000_sppd[115];

    auto g_xz_0_0_0_0_x_y_xz = buffer_2000_sppd[116];

    auto g_xz_0_0_0_0_x_y_yy = buffer_2000_sppd[117];

    auto g_xz_0_0_0_0_x_y_yz = buffer_2000_sppd[118];

    auto g_xz_0_0_0_0_x_y_zz = buffer_2000_sppd[119];

    auto g_xz_0_0_0_0_x_z_xx = buffer_2000_sppd[120];

    auto g_xz_0_0_0_0_x_z_xy = buffer_2000_sppd[121];

    auto g_xz_0_0_0_0_x_z_xz = buffer_2000_sppd[122];

    auto g_xz_0_0_0_0_x_z_yy = buffer_2000_sppd[123];

    auto g_xz_0_0_0_0_x_z_yz = buffer_2000_sppd[124];

    auto g_xz_0_0_0_0_x_z_zz = buffer_2000_sppd[125];

    auto g_xz_0_0_0_0_y_x_xx = buffer_2000_sppd[126];

    auto g_xz_0_0_0_0_y_x_xy = buffer_2000_sppd[127];

    auto g_xz_0_0_0_0_y_x_xz = buffer_2000_sppd[128];

    auto g_xz_0_0_0_0_y_x_yy = buffer_2000_sppd[129];

    auto g_xz_0_0_0_0_y_x_yz = buffer_2000_sppd[130];

    auto g_xz_0_0_0_0_y_x_zz = buffer_2000_sppd[131];

    auto g_xz_0_0_0_0_y_y_xx = buffer_2000_sppd[132];

    auto g_xz_0_0_0_0_y_y_xy = buffer_2000_sppd[133];

    auto g_xz_0_0_0_0_y_y_xz = buffer_2000_sppd[134];

    auto g_xz_0_0_0_0_y_y_yy = buffer_2000_sppd[135];

    auto g_xz_0_0_0_0_y_y_yz = buffer_2000_sppd[136];

    auto g_xz_0_0_0_0_y_y_zz = buffer_2000_sppd[137];

    auto g_xz_0_0_0_0_y_z_xx = buffer_2000_sppd[138];

    auto g_xz_0_0_0_0_y_z_xy = buffer_2000_sppd[139];

    auto g_xz_0_0_0_0_y_z_xz = buffer_2000_sppd[140];

    auto g_xz_0_0_0_0_y_z_yy = buffer_2000_sppd[141];

    auto g_xz_0_0_0_0_y_z_yz = buffer_2000_sppd[142];

    auto g_xz_0_0_0_0_y_z_zz = buffer_2000_sppd[143];

    auto g_xz_0_0_0_0_z_x_xx = buffer_2000_sppd[144];

    auto g_xz_0_0_0_0_z_x_xy = buffer_2000_sppd[145];

    auto g_xz_0_0_0_0_z_x_xz = buffer_2000_sppd[146];

    auto g_xz_0_0_0_0_z_x_yy = buffer_2000_sppd[147];

    auto g_xz_0_0_0_0_z_x_yz = buffer_2000_sppd[148];

    auto g_xz_0_0_0_0_z_x_zz = buffer_2000_sppd[149];

    auto g_xz_0_0_0_0_z_y_xx = buffer_2000_sppd[150];

    auto g_xz_0_0_0_0_z_y_xy = buffer_2000_sppd[151];

    auto g_xz_0_0_0_0_z_y_xz = buffer_2000_sppd[152];

    auto g_xz_0_0_0_0_z_y_yy = buffer_2000_sppd[153];

    auto g_xz_0_0_0_0_z_y_yz = buffer_2000_sppd[154];

    auto g_xz_0_0_0_0_z_y_zz = buffer_2000_sppd[155];

    auto g_xz_0_0_0_0_z_z_xx = buffer_2000_sppd[156];

    auto g_xz_0_0_0_0_z_z_xy = buffer_2000_sppd[157];

    auto g_xz_0_0_0_0_z_z_xz = buffer_2000_sppd[158];

    auto g_xz_0_0_0_0_z_z_yy = buffer_2000_sppd[159];

    auto g_xz_0_0_0_0_z_z_yz = buffer_2000_sppd[160];

    auto g_xz_0_0_0_0_z_z_zz = buffer_2000_sppd[161];

    auto g_yy_0_0_0_0_x_x_xx = buffer_2000_sppd[162];

    auto g_yy_0_0_0_0_x_x_xy = buffer_2000_sppd[163];

    auto g_yy_0_0_0_0_x_x_xz = buffer_2000_sppd[164];

    auto g_yy_0_0_0_0_x_x_yy = buffer_2000_sppd[165];

    auto g_yy_0_0_0_0_x_x_yz = buffer_2000_sppd[166];

    auto g_yy_0_0_0_0_x_x_zz = buffer_2000_sppd[167];

    auto g_yy_0_0_0_0_x_y_xx = buffer_2000_sppd[168];

    auto g_yy_0_0_0_0_x_y_xy = buffer_2000_sppd[169];

    auto g_yy_0_0_0_0_x_y_xz = buffer_2000_sppd[170];

    auto g_yy_0_0_0_0_x_y_yy = buffer_2000_sppd[171];

    auto g_yy_0_0_0_0_x_y_yz = buffer_2000_sppd[172];

    auto g_yy_0_0_0_0_x_y_zz = buffer_2000_sppd[173];

    auto g_yy_0_0_0_0_x_z_xx = buffer_2000_sppd[174];

    auto g_yy_0_0_0_0_x_z_xy = buffer_2000_sppd[175];

    auto g_yy_0_0_0_0_x_z_xz = buffer_2000_sppd[176];

    auto g_yy_0_0_0_0_x_z_yy = buffer_2000_sppd[177];

    auto g_yy_0_0_0_0_x_z_yz = buffer_2000_sppd[178];

    auto g_yy_0_0_0_0_x_z_zz = buffer_2000_sppd[179];

    auto g_yy_0_0_0_0_y_x_xx = buffer_2000_sppd[180];

    auto g_yy_0_0_0_0_y_x_xy = buffer_2000_sppd[181];

    auto g_yy_0_0_0_0_y_x_xz = buffer_2000_sppd[182];

    auto g_yy_0_0_0_0_y_x_yy = buffer_2000_sppd[183];

    auto g_yy_0_0_0_0_y_x_yz = buffer_2000_sppd[184];

    auto g_yy_0_0_0_0_y_x_zz = buffer_2000_sppd[185];

    auto g_yy_0_0_0_0_y_y_xx = buffer_2000_sppd[186];

    auto g_yy_0_0_0_0_y_y_xy = buffer_2000_sppd[187];

    auto g_yy_0_0_0_0_y_y_xz = buffer_2000_sppd[188];

    auto g_yy_0_0_0_0_y_y_yy = buffer_2000_sppd[189];

    auto g_yy_0_0_0_0_y_y_yz = buffer_2000_sppd[190];

    auto g_yy_0_0_0_0_y_y_zz = buffer_2000_sppd[191];

    auto g_yy_0_0_0_0_y_z_xx = buffer_2000_sppd[192];

    auto g_yy_0_0_0_0_y_z_xy = buffer_2000_sppd[193];

    auto g_yy_0_0_0_0_y_z_xz = buffer_2000_sppd[194];

    auto g_yy_0_0_0_0_y_z_yy = buffer_2000_sppd[195];

    auto g_yy_0_0_0_0_y_z_yz = buffer_2000_sppd[196];

    auto g_yy_0_0_0_0_y_z_zz = buffer_2000_sppd[197];

    auto g_yy_0_0_0_0_z_x_xx = buffer_2000_sppd[198];

    auto g_yy_0_0_0_0_z_x_xy = buffer_2000_sppd[199];

    auto g_yy_0_0_0_0_z_x_xz = buffer_2000_sppd[200];

    auto g_yy_0_0_0_0_z_x_yy = buffer_2000_sppd[201];

    auto g_yy_0_0_0_0_z_x_yz = buffer_2000_sppd[202];

    auto g_yy_0_0_0_0_z_x_zz = buffer_2000_sppd[203];

    auto g_yy_0_0_0_0_z_y_xx = buffer_2000_sppd[204];

    auto g_yy_0_0_0_0_z_y_xy = buffer_2000_sppd[205];

    auto g_yy_0_0_0_0_z_y_xz = buffer_2000_sppd[206];

    auto g_yy_0_0_0_0_z_y_yy = buffer_2000_sppd[207];

    auto g_yy_0_0_0_0_z_y_yz = buffer_2000_sppd[208];

    auto g_yy_0_0_0_0_z_y_zz = buffer_2000_sppd[209];

    auto g_yy_0_0_0_0_z_z_xx = buffer_2000_sppd[210];

    auto g_yy_0_0_0_0_z_z_xy = buffer_2000_sppd[211];

    auto g_yy_0_0_0_0_z_z_xz = buffer_2000_sppd[212];

    auto g_yy_0_0_0_0_z_z_yy = buffer_2000_sppd[213];

    auto g_yy_0_0_0_0_z_z_yz = buffer_2000_sppd[214];

    auto g_yy_0_0_0_0_z_z_zz = buffer_2000_sppd[215];

    auto g_yz_0_0_0_0_x_x_xx = buffer_2000_sppd[216];

    auto g_yz_0_0_0_0_x_x_xy = buffer_2000_sppd[217];

    auto g_yz_0_0_0_0_x_x_xz = buffer_2000_sppd[218];

    auto g_yz_0_0_0_0_x_x_yy = buffer_2000_sppd[219];

    auto g_yz_0_0_0_0_x_x_yz = buffer_2000_sppd[220];

    auto g_yz_0_0_0_0_x_x_zz = buffer_2000_sppd[221];

    auto g_yz_0_0_0_0_x_y_xx = buffer_2000_sppd[222];

    auto g_yz_0_0_0_0_x_y_xy = buffer_2000_sppd[223];

    auto g_yz_0_0_0_0_x_y_xz = buffer_2000_sppd[224];

    auto g_yz_0_0_0_0_x_y_yy = buffer_2000_sppd[225];

    auto g_yz_0_0_0_0_x_y_yz = buffer_2000_sppd[226];

    auto g_yz_0_0_0_0_x_y_zz = buffer_2000_sppd[227];

    auto g_yz_0_0_0_0_x_z_xx = buffer_2000_sppd[228];

    auto g_yz_0_0_0_0_x_z_xy = buffer_2000_sppd[229];

    auto g_yz_0_0_0_0_x_z_xz = buffer_2000_sppd[230];

    auto g_yz_0_0_0_0_x_z_yy = buffer_2000_sppd[231];

    auto g_yz_0_0_0_0_x_z_yz = buffer_2000_sppd[232];

    auto g_yz_0_0_0_0_x_z_zz = buffer_2000_sppd[233];

    auto g_yz_0_0_0_0_y_x_xx = buffer_2000_sppd[234];

    auto g_yz_0_0_0_0_y_x_xy = buffer_2000_sppd[235];

    auto g_yz_0_0_0_0_y_x_xz = buffer_2000_sppd[236];

    auto g_yz_0_0_0_0_y_x_yy = buffer_2000_sppd[237];

    auto g_yz_0_0_0_0_y_x_yz = buffer_2000_sppd[238];

    auto g_yz_0_0_0_0_y_x_zz = buffer_2000_sppd[239];

    auto g_yz_0_0_0_0_y_y_xx = buffer_2000_sppd[240];

    auto g_yz_0_0_0_0_y_y_xy = buffer_2000_sppd[241];

    auto g_yz_0_0_0_0_y_y_xz = buffer_2000_sppd[242];

    auto g_yz_0_0_0_0_y_y_yy = buffer_2000_sppd[243];

    auto g_yz_0_0_0_0_y_y_yz = buffer_2000_sppd[244];

    auto g_yz_0_0_0_0_y_y_zz = buffer_2000_sppd[245];

    auto g_yz_0_0_0_0_y_z_xx = buffer_2000_sppd[246];

    auto g_yz_0_0_0_0_y_z_xy = buffer_2000_sppd[247];

    auto g_yz_0_0_0_0_y_z_xz = buffer_2000_sppd[248];

    auto g_yz_0_0_0_0_y_z_yy = buffer_2000_sppd[249];

    auto g_yz_0_0_0_0_y_z_yz = buffer_2000_sppd[250];

    auto g_yz_0_0_0_0_y_z_zz = buffer_2000_sppd[251];

    auto g_yz_0_0_0_0_z_x_xx = buffer_2000_sppd[252];

    auto g_yz_0_0_0_0_z_x_xy = buffer_2000_sppd[253];

    auto g_yz_0_0_0_0_z_x_xz = buffer_2000_sppd[254];

    auto g_yz_0_0_0_0_z_x_yy = buffer_2000_sppd[255];

    auto g_yz_0_0_0_0_z_x_yz = buffer_2000_sppd[256];

    auto g_yz_0_0_0_0_z_x_zz = buffer_2000_sppd[257];

    auto g_yz_0_0_0_0_z_y_xx = buffer_2000_sppd[258];

    auto g_yz_0_0_0_0_z_y_xy = buffer_2000_sppd[259];

    auto g_yz_0_0_0_0_z_y_xz = buffer_2000_sppd[260];

    auto g_yz_0_0_0_0_z_y_yy = buffer_2000_sppd[261];

    auto g_yz_0_0_0_0_z_y_yz = buffer_2000_sppd[262];

    auto g_yz_0_0_0_0_z_y_zz = buffer_2000_sppd[263];

    auto g_yz_0_0_0_0_z_z_xx = buffer_2000_sppd[264];

    auto g_yz_0_0_0_0_z_z_xy = buffer_2000_sppd[265];

    auto g_yz_0_0_0_0_z_z_xz = buffer_2000_sppd[266];

    auto g_yz_0_0_0_0_z_z_yy = buffer_2000_sppd[267];

    auto g_yz_0_0_0_0_z_z_yz = buffer_2000_sppd[268];

    auto g_yz_0_0_0_0_z_z_zz = buffer_2000_sppd[269];

    auto g_zz_0_0_0_0_x_x_xx = buffer_2000_sppd[270];

    auto g_zz_0_0_0_0_x_x_xy = buffer_2000_sppd[271];

    auto g_zz_0_0_0_0_x_x_xz = buffer_2000_sppd[272];

    auto g_zz_0_0_0_0_x_x_yy = buffer_2000_sppd[273];

    auto g_zz_0_0_0_0_x_x_yz = buffer_2000_sppd[274];

    auto g_zz_0_0_0_0_x_x_zz = buffer_2000_sppd[275];

    auto g_zz_0_0_0_0_x_y_xx = buffer_2000_sppd[276];

    auto g_zz_0_0_0_0_x_y_xy = buffer_2000_sppd[277];

    auto g_zz_0_0_0_0_x_y_xz = buffer_2000_sppd[278];

    auto g_zz_0_0_0_0_x_y_yy = buffer_2000_sppd[279];

    auto g_zz_0_0_0_0_x_y_yz = buffer_2000_sppd[280];

    auto g_zz_0_0_0_0_x_y_zz = buffer_2000_sppd[281];

    auto g_zz_0_0_0_0_x_z_xx = buffer_2000_sppd[282];

    auto g_zz_0_0_0_0_x_z_xy = buffer_2000_sppd[283];

    auto g_zz_0_0_0_0_x_z_xz = buffer_2000_sppd[284];

    auto g_zz_0_0_0_0_x_z_yy = buffer_2000_sppd[285];

    auto g_zz_0_0_0_0_x_z_yz = buffer_2000_sppd[286];

    auto g_zz_0_0_0_0_x_z_zz = buffer_2000_sppd[287];

    auto g_zz_0_0_0_0_y_x_xx = buffer_2000_sppd[288];

    auto g_zz_0_0_0_0_y_x_xy = buffer_2000_sppd[289];

    auto g_zz_0_0_0_0_y_x_xz = buffer_2000_sppd[290];

    auto g_zz_0_0_0_0_y_x_yy = buffer_2000_sppd[291];

    auto g_zz_0_0_0_0_y_x_yz = buffer_2000_sppd[292];

    auto g_zz_0_0_0_0_y_x_zz = buffer_2000_sppd[293];

    auto g_zz_0_0_0_0_y_y_xx = buffer_2000_sppd[294];

    auto g_zz_0_0_0_0_y_y_xy = buffer_2000_sppd[295];

    auto g_zz_0_0_0_0_y_y_xz = buffer_2000_sppd[296];

    auto g_zz_0_0_0_0_y_y_yy = buffer_2000_sppd[297];

    auto g_zz_0_0_0_0_y_y_yz = buffer_2000_sppd[298];

    auto g_zz_0_0_0_0_y_y_zz = buffer_2000_sppd[299];

    auto g_zz_0_0_0_0_y_z_xx = buffer_2000_sppd[300];

    auto g_zz_0_0_0_0_y_z_xy = buffer_2000_sppd[301];

    auto g_zz_0_0_0_0_y_z_xz = buffer_2000_sppd[302];

    auto g_zz_0_0_0_0_y_z_yy = buffer_2000_sppd[303];

    auto g_zz_0_0_0_0_y_z_yz = buffer_2000_sppd[304];

    auto g_zz_0_0_0_0_y_z_zz = buffer_2000_sppd[305];

    auto g_zz_0_0_0_0_z_x_xx = buffer_2000_sppd[306];

    auto g_zz_0_0_0_0_z_x_xy = buffer_2000_sppd[307];

    auto g_zz_0_0_0_0_z_x_xz = buffer_2000_sppd[308];

    auto g_zz_0_0_0_0_z_x_yy = buffer_2000_sppd[309];

    auto g_zz_0_0_0_0_z_x_yz = buffer_2000_sppd[310];

    auto g_zz_0_0_0_0_z_x_zz = buffer_2000_sppd[311];

    auto g_zz_0_0_0_0_z_y_xx = buffer_2000_sppd[312];

    auto g_zz_0_0_0_0_z_y_xy = buffer_2000_sppd[313];

    auto g_zz_0_0_0_0_z_y_xz = buffer_2000_sppd[314];

    auto g_zz_0_0_0_0_z_y_yy = buffer_2000_sppd[315];

    auto g_zz_0_0_0_0_z_y_yz = buffer_2000_sppd[316];

    auto g_zz_0_0_0_0_z_y_zz = buffer_2000_sppd[317];

    auto g_zz_0_0_0_0_z_z_xx = buffer_2000_sppd[318];

    auto g_zz_0_0_0_0_z_z_xy = buffer_2000_sppd[319];

    auto g_zz_0_0_0_0_z_z_xz = buffer_2000_sppd[320];

    auto g_zz_0_0_0_0_z_z_yy = buffer_2000_sppd[321];

    auto g_zz_0_0_0_0_z_z_yz = buffer_2000_sppd[322];

    auto g_zz_0_0_0_0_z_z_zz = buffer_2000_sppd[323];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_xx_0_0_0_0_x_x_xx, g_xx_0_0_0_0_x_x_xy, g_xx_0_0_0_0_x_x_xz, g_xx_0_0_0_0_x_x_yy, g_xx_0_0_0_0_x_x_yz, g_xx_0_0_0_0_x_x_zz, g_xx_x_x_xx, g_xx_x_x_xy, g_xx_x_x_xz, g_xx_x_x_yy, g_xx_x_x_yz, g_xx_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_x_x_xx[i] = -2.0 * g_0_x_x_xx[i] * a_exp + 4.0 * g_xx_x_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_x_xy[i] = -2.0 * g_0_x_x_xy[i] * a_exp + 4.0 * g_xx_x_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_x_xz[i] = -2.0 * g_0_x_x_xz[i] * a_exp + 4.0 * g_xx_x_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_x_yy[i] = -2.0 * g_0_x_x_yy[i] * a_exp + 4.0 * g_xx_x_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_x_yz[i] = -2.0 * g_0_x_x_yz[i] * a_exp + 4.0 * g_xx_x_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_x_zz[i] = -2.0 * g_0_x_x_zz[i] * a_exp + 4.0 * g_xx_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_xx_0_0_0_0_x_y_xx, g_xx_0_0_0_0_x_y_xy, g_xx_0_0_0_0_x_y_xz, g_xx_0_0_0_0_x_y_yy, g_xx_0_0_0_0_x_y_yz, g_xx_0_0_0_0_x_y_zz, g_xx_x_y_xx, g_xx_x_y_xy, g_xx_x_y_xz, g_xx_x_y_yy, g_xx_x_y_yz, g_xx_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_x_y_xx[i] = -2.0 * g_0_x_y_xx[i] * a_exp + 4.0 * g_xx_x_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_y_xy[i] = -2.0 * g_0_x_y_xy[i] * a_exp + 4.0 * g_xx_x_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_y_xz[i] = -2.0 * g_0_x_y_xz[i] * a_exp + 4.0 * g_xx_x_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_y_yy[i] = -2.0 * g_0_x_y_yy[i] * a_exp + 4.0 * g_xx_x_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_y_yz[i] = -2.0 * g_0_x_y_yz[i] * a_exp + 4.0 * g_xx_x_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_y_zz[i] = -2.0 * g_0_x_y_zz[i] * a_exp + 4.0 * g_xx_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_xx_0_0_0_0_x_z_xx, g_xx_0_0_0_0_x_z_xy, g_xx_0_0_0_0_x_z_xz, g_xx_0_0_0_0_x_z_yy, g_xx_0_0_0_0_x_z_yz, g_xx_0_0_0_0_x_z_zz, g_xx_x_z_xx, g_xx_x_z_xy, g_xx_x_z_xz, g_xx_x_z_yy, g_xx_x_z_yz, g_xx_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_x_z_xx[i] = -2.0 * g_0_x_z_xx[i] * a_exp + 4.0 * g_xx_x_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_z_xy[i] = -2.0 * g_0_x_z_xy[i] * a_exp + 4.0 * g_xx_x_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_z_xz[i] = -2.0 * g_0_x_z_xz[i] * a_exp + 4.0 * g_xx_x_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_z_yy[i] = -2.0 * g_0_x_z_yy[i] * a_exp + 4.0 * g_xx_x_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_z_yz[i] = -2.0 * g_0_x_z_yz[i] * a_exp + 4.0 * g_xx_x_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_z_zz[i] = -2.0 * g_0_x_z_zz[i] * a_exp + 4.0 * g_xx_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_xx_0_0_0_0_y_x_xx, g_xx_0_0_0_0_y_x_xy, g_xx_0_0_0_0_y_x_xz, g_xx_0_0_0_0_y_x_yy, g_xx_0_0_0_0_y_x_yz, g_xx_0_0_0_0_y_x_zz, g_xx_y_x_xx, g_xx_y_x_xy, g_xx_y_x_xz, g_xx_y_x_yy, g_xx_y_x_yz, g_xx_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_y_x_xx[i] = -2.0 * g_0_y_x_xx[i] * a_exp + 4.0 * g_xx_y_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_x_xy[i] = -2.0 * g_0_y_x_xy[i] * a_exp + 4.0 * g_xx_y_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_x_xz[i] = -2.0 * g_0_y_x_xz[i] * a_exp + 4.0 * g_xx_y_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_x_yy[i] = -2.0 * g_0_y_x_yy[i] * a_exp + 4.0 * g_xx_y_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_x_yz[i] = -2.0 * g_0_y_x_yz[i] * a_exp + 4.0 * g_xx_y_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_x_zz[i] = -2.0 * g_0_y_x_zz[i] * a_exp + 4.0 * g_xx_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_xx_0_0_0_0_y_y_xx, g_xx_0_0_0_0_y_y_xy, g_xx_0_0_0_0_y_y_xz, g_xx_0_0_0_0_y_y_yy, g_xx_0_0_0_0_y_y_yz, g_xx_0_0_0_0_y_y_zz, g_xx_y_y_xx, g_xx_y_y_xy, g_xx_y_y_xz, g_xx_y_y_yy, g_xx_y_y_yz, g_xx_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_y_y_xx[i] = -2.0 * g_0_y_y_xx[i] * a_exp + 4.0 * g_xx_y_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_y_xy[i] = -2.0 * g_0_y_y_xy[i] * a_exp + 4.0 * g_xx_y_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_y_xz[i] = -2.0 * g_0_y_y_xz[i] * a_exp + 4.0 * g_xx_y_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_y_yy[i] = -2.0 * g_0_y_y_yy[i] * a_exp + 4.0 * g_xx_y_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_y_yz[i] = -2.0 * g_0_y_y_yz[i] * a_exp + 4.0 * g_xx_y_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_y_zz[i] = -2.0 * g_0_y_y_zz[i] * a_exp + 4.0 * g_xx_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_xx_0_0_0_0_y_z_xx, g_xx_0_0_0_0_y_z_xy, g_xx_0_0_0_0_y_z_xz, g_xx_0_0_0_0_y_z_yy, g_xx_0_0_0_0_y_z_yz, g_xx_0_0_0_0_y_z_zz, g_xx_y_z_xx, g_xx_y_z_xy, g_xx_y_z_xz, g_xx_y_z_yy, g_xx_y_z_yz, g_xx_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_y_z_xx[i] = -2.0 * g_0_y_z_xx[i] * a_exp + 4.0 * g_xx_y_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_z_xy[i] = -2.0 * g_0_y_z_xy[i] * a_exp + 4.0 * g_xx_y_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_z_xz[i] = -2.0 * g_0_y_z_xz[i] * a_exp + 4.0 * g_xx_y_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_z_yy[i] = -2.0 * g_0_y_z_yy[i] * a_exp + 4.0 * g_xx_y_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_z_yz[i] = -2.0 * g_0_y_z_yz[i] * a_exp + 4.0 * g_xx_y_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_z_zz[i] = -2.0 * g_0_y_z_zz[i] * a_exp + 4.0 * g_xx_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_xx_0_0_0_0_z_x_xx, g_xx_0_0_0_0_z_x_xy, g_xx_0_0_0_0_z_x_xz, g_xx_0_0_0_0_z_x_yy, g_xx_0_0_0_0_z_x_yz, g_xx_0_0_0_0_z_x_zz, g_xx_z_x_xx, g_xx_z_x_xy, g_xx_z_x_xz, g_xx_z_x_yy, g_xx_z_x_yz, g_xx_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_z_x_xx[i] = -2.0 * g_0_z_x_xx[i] * a_exp + 4.0 * g_xx_z_x_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_x_xy[i] = -2.0 * g_0_z_x_xy[i] * a_exp + 4.0 * g_xx_z_x_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_x_xz[i] = -2.0 * g_0_z_x_xz[i] * a_exp + 4.0 * g_xx_z_x_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_x_yy[i] = -2.0 * g_0_z_x_yy[i] * a_exp + 4.0 * g_xx_z_x_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_x_yz[i] = -2.0 * g_0_z_x_yz[i] * a_exp + 4.0 * g_xx_z_x_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_x_zz[i] = -2.0 * g_0_z_x_zz[i] * a_exp + 4.0 * g_xx_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_xx_0_0_0_0_z_y_xx, g_xx_0_0_0_0_z_y_xy, g_xx_0_0_0_0_z_y_xz, g_xx_0_0_0_0_z_y_yy, g_xx_0_0_0_0_z_y_yz, g_xx_0_0_0_0_z_y_zz, g_xx_z_y_xx, g_xx_z_y_xy, g_xx_z_y_xz, g_xx_z_y_yy, g_xx_z_y_yz, g_xx_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_z_y_xx[i] = -2.0 * g_0_z_y_xx[i] * a_exp + 4.0 * g_xx_z_y_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_y_xy[i] = -2.0 * g_0_z_y_xy[i] * a_exp + 4.0 * g_xx_z_y_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_y_xz[i] = -2.0 * g_0_z_y_xz[i] * a_exp + 4.0 * g_xx_z_y_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_y_yy[i] = -2.0 * g_0_z_y_yy[i] * a_exp + 4.0 * g_xx_z_y_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_y_yz[i] = -2.0 * g_0_z_y_yz[i] * a_exp + 4.0 * g_xx_z_y_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_y_zz[i] = -2.0 * g_0_z_y_zz[i] * a_exp + 4.0 * g_xx_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_xx_0_0_0_0_z_z_xx, g_xx_0_0_0_0_z_z_xy, g_xx_0_0_0_0_z_z_xz, g_xx_0_0_0_0_z_z_yy, g_xx_0_0_0_0_z_z_yz, g_xx_0_0_0_0_z_z_zz, g_xx_z_z_xx, g_xx_z_z_xy, g_xx_z_z_xz, g_xx_z_z_yy, g_xx_z_z_yz, g_xx_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_z_z_xx[i] = -2.0 * g_0_z_z_xx[i] * a_exp + 4.0 * g_xx_z_z_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_z_xy[i] = -2.0 * g_0_z_z_xy[i] * a_exp + 4.0 * g_xx_z_z_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_z_xz[i] = -2.0 * g_0_z_z_xz[i] * a_exp + 4.0 * g_xx_z_z_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_z_yy[i] = -2.0 * g_0_z_z_yy[i] * a_exp + 4.0 * g_xx_z_z_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_z_yz[i] = -2.0 * g_0_z_z_yz[i] * a_exp + 4.0 * g_xx_z_z_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_z_zz[i] = -2.0 * g_0_z_z_zz[i] * a_exp + 4.0 * g_xx_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_xy_0_0_0_0_x_x_xx, g_xy_0_0_0_0_x_x_xy, g_xy_0_0_0_0_x_x_xz, g_xy_0_0_0_0_x_x_yy, g_xy_0_0_0_0_x_x_yz, g_xy_0_0_0_0_x_x_zz, g_xy_x_x_xx, g_xy_x_x_xy, g_xy_x_x_xz, g_xy_x_x_yy, g_xy_x_x_yz, g_xy_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_x_x_xx[i] = 4.0 * g_xy_x_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_x_xy[i] = 4.0 * g_xy_x_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_x_xz[i] = 4.0 * g_xy_x_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_x_yy[i] = 4.0 * g_xy_x_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_x_yz[i] = 4.0 * g_xy_x_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_x_zz[i] = 4.0 * g_xy_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_xy_0_0_0_0_x_y_xx, g_xy_0_0_0_0_x_y_xy, g_xy_0_0_0_0_x_y_xz, g_xy_0_0_0_0_x_y_yy, g_xy_0_0_0_0_x_y_yz, g_xy_0_0_0_0_x_y_zz, g_xy_x_y_xx, g_xy_x_y_xy, g_xy_x_y_xz, g_xy_x_y_yy, g_xy_x_y_yz, g_xy_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_x_y_xx[i] = 4.0 * g_xy_x_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_y_xy[i] = 4.0 * g_xy_x_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_y_xz[i] = 4.0 * g_xy_x_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_y_yy[i] = 4.0 * g_xy_x_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_y_yz[i] = 4.0 * g_xy_x_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_y_zz[i] = 4.0 * g_xy_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_xy_0_0_0_0_x_z_xx, g_xy_0_0_0_0_x_z_xy, g_xy_0_0_0_0_x_z_xz, g_xy_0_0_0_0_x_z_yy, g_xy_0_0_0_0_x_z_yz, g_xy_0_0_0_0_x_z_zz, g_xy_x_z_xx, g_xy_x_z_xy, g_xy_x_z_xz, g_xy_x_z_yy, g_xy_x_z_yz, g_xy_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_x_z_xx[i] = 4.0 * g_xy_x_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_z_xy[i] = 4.0 * g_xy_x_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_z_xz[i] = 4.0 * g_xy_x_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_z_yy[i] = 4.0 * g_xy_x_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_z_yz[i] = 4.0 * g_xy_x_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_z_zz[i] = 4.0 * g_xy_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_xy_0_0_0_0_y_x_xx, g_xy_0_0_0_0_y_x_xy, g_xy_0_0_0_0_y_x_xz, g_xy_0_0_0_0_y_x_yy, g_xy_0_0_0_0_y_x_yz, g_xy_0_0_0_0_y_x_zz, g_xy_y_x_xx, g_xy_y_x_xy, g_xy_y_x_xz, g_xy_y_x_yy, g_xy_y_x_yz, g_xy_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_y_x_xx[i] = 4.0 * g_xy_y_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_x_xy[i] = 4.0 * g_xy_y_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_x_xz[i] = 4.0 * g_xy_y_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_x_yy[i] = 4.0 * g_xy_y_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_x_yz[i] = 4.0 * g_xy_y_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_x_zz[i] = 4.0 * g_xy_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_xy_0_0_0_0_y_y_xx, g_xy_0_0_0_0_y_y_xy, g_xy_0_0_0_0_y_y_xz, g_xy_0_0_0_0_y_y_yy, g_xy_0_0_0_0_y_y_yz, g_xy_0_0_0_0_y_y_zz, g_xy_y_y_xx, g_xy_y_y_xy, g_xy_y_y_xz, g_xy_y_y_yy, g_xy_y_y_yz, g_xy_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_y_y_xx[i] = 4.0 * g_xy_y_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_y_xy[i] = 4.0 * g_xy_y_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_y_xz[i] = 4.0 * g_xy_y_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_y_yy[i] = 4.0 * g_xy_y_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_y_yz[i] = 4.0 * g_xy_y_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_y_zz[i] = 4.0 * g_xy_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_xy_0_0_0_0_y_z_xx, g_xy_0_0_0_0_y_z_xy, g_xy_0_0_0_0_y_z_xz, g_xy_0_0_0_0_y_z_yy, g_xy_0_0_0_0_y_z_yz, g_xy_0_0_0_0_y_z_zz, g_xy_y_z_xx, g_xy_y_z_xy, g_xy_y_z_xz, g_xy_y_z_yy, g_xy_y_z_yz, g_xy_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_y_z_xx[i] = 4.0 * g_xy_y_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_z_xy[i] = 4.0 * g_xy_y_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_z_xz[i] = 4.0 * g_xy_y_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_z_yy[i] = 4.0 * g_xy_y_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_z_yz[i] = 4.0 * g_xy_y_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_z_zz[i] = 4.0 * g_xy_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_xy_0_0_0_0_z_x_xx, g_xy_0_0_0_0_z_x_xy, g_xy_0_0_0_0_z_x_xz, g_xy_0_0_0_0_z_x_yy, g_xy_0_0_0_0_z_x_yz, g_xy_0_0_0_0_z_x_zz, g_xy_z_x_xx, g_xy_z_x_xy, g_xy_z_x_xz, g_xy_z_x_yy, g_xy_z_x_yz, g_xy_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_z_x_xx[i] = 4.0 * g_xy_z_x_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_x_xy[i] = 4.0 * g_xy_z_x_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_x_xz[i] = 4.0 * g_xy_z_x_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_x_yy[i] = 4.0 * g_xy_z_x_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_x_yz[i] = 4.0 * g_xy_z_x_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_x_zz[i] = 4.0 * g_xy_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_xy_0_0_0_0_z_y_xx, g_xy_0_0_0_0_z_y_xy, g_xy_0_0_0_0_z_y_xz, g_xy_0_0_0_0_z_y_yy, g_xy_0_0_0_0_z_y_yz, g_xy_0_0_0_0_z_y_zz, g_xy_z_y_xx, g_xy_z_y_xy, g_xy_z_y_xz, g_xy_z_y_yy, g_xy_z_y_yz, g_xy_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_z_y_xx[i] = 4.0 * g_xy_z_y_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_y_xy[i] = 4.0 * g_xy_z_y_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_y_xz[i] = 4.0 * g_xy_z_y_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_y_yy[i] = 4.0 * g_xy_z_y_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_y_yz[i] = 4.0 * g_xy_z_y_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_y_zz[i] = 4.0 * g_xy_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_xy_0_0_0_0_z_z_xx, g_xy_0_0_0_0_z_z_xy, g_xy_0_0_0_0_z_z_xz, g_xy_0_0_0_0_z_z_yy, g_xy_0_0_0_0_z_z_yz, g_xy_0_0_0_0_z_z_zz, g_xy_z_z_xx, g_xy_z_z_xy, g_xy_z_z_xz, g_xy_z_z_yy, g_xy_z_z_yz, g_xy_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_z_z_xx[i] = 4.0 * g_xy_z_z_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_z_xy[i] = 4.0 * g_xy_z_z_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_z_xz[i] = 4.0 * g_xy_z_z_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_z_yy[i] = 4.0 * g_xy_z_z_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_z_yz[i] = 4.0 * g_xy_z_z_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_z_zz[i] = 4.0 * g_xy_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_xz_0_0_0_0_x_x_xx, g_xz_0_0_0_0_x_x_xy, g_xz_0_0_0_0_x_x_xz, g_xz_0_0_0_0_x_x_yy, g_xz_0_0_0_0_x_x_yz, g_xz_0_0_0_0_x_x_zz, g_xz_x_x_xx, g_xz_x_x_xy, g_xz_x_x_xz, g_xz_x_x_yy, g_xz_x_x_yz, g_xz_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_x_x_xx[i] = 4.0 * g_xz_x_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_x_xy[i] = 4.0 * g_xz_x_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_x_xz[i] = 4.0 * g_xz_x_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_x_yy[i] = 4.0 * g_xz_x_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_x_yz[i] = 4.0 * g_xz_x_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_x_zz[i] = 4.0 * g_xz_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_xz_0_0_0_0_x_y_xx, g_xz_0_0_0_0_x_y_xy, g_xz_0_0_0_0_x_y_xz, g_xz_0_0_0_0_x_y_yy, g_xz_0_0_0_0_x_y_yz, g_xz_0_0_0_0_x_y_zz, g_xz_x_y_xx, g_xz_x_y_xy, g_xz_x_y_xz, g_xz_x_y_yy, g_xz_x_y_yz, g_xz_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_x_y_xx[i] = 4.0 * g_xz_x_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_y_xy[i] = 4.0 * g_xz_x_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_y_xz[i] = 4.0 * g_xz_x_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_y_yy[i] = 4.0 * g_xz_x_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_y_yz[i] = 4.0 * g_xz_x_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_y_zz[i] = 4.0 * g_xz_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_xz_0_0_0_0_x_z_xx, g_xz_0_0_0_0_x_z_xy, g_xz_0_0_0_0_x_z_xz, g_xz_0_0_0_0_x_z_yy, g_xz_0_0_0_0_x_z_yz, g_xz_0_0_0_0_x_z_zz, g_xz_x_z_xx, g_xz_x_z_xy, g_xz_x_z_xz, g_xz_x_z_yy, g_xz_x_z_yz, g_xz_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_x_z_xx[i] = 4.0 * g_xz_x_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_z_xy[i] = 4.0 * g_xz_x_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_z_xz[i] = 4.0 * g_xz_x_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_z_yy[i] = 4.0 * g_xz_x_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_z_yz[i] = 4.0 * g_xz_x_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_z_zz[i] = 4.0 * g_xz_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_xz_0_0_0_0_y_x_xx, g_xz_0_0_0_0_y_x_xy, g_xz_0_0_0_0_y_x_xz, g_xz_0_0_0_0_y_x_yy, g_xz_0_0_0_0_y_x_yz, g_xz_0_0_0_0_y_x_zz, g_xz_y_x_xx, g_xz_y_x_xy, g_xz_y_x_xz, g_xz_y_x_yy, g_xz_y_x_yz, g_xz_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_y_x_xx[i] = 4.0 * g_xz_y_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_x_xy[i] = 4.0 * g_xz_y_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_x_xz[i] = 4.0 * g_xz_y_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_x_yy[i] = 4.0 * g_xz_y_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_x_yz[i] = 4.0 * g_xz_y_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_x_zz[i] = 4.0 * g_xz_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_xz_0_0_0_0_y_y_xx, g_xz_0_0_0_0_y_y_xy, g_xz_0_0_0_0_y_y_xz, g_xz_0_0_0_0_y_y_yy, g_xz_0_0_0_0_y_y_yz, g_xz_0_0_0_0_y_y_zz, g_xz_y_y_xx, g_xz_y_y_xy, g_xz_y_y_xz, g_xz_y_y_yy, g_xz_y_y_yz, g_xz_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_y_y_xx[i] = 4.0 * g_xz_y_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_y_xy[i] = 4.0 * g_xz_y_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_y_xz[i] = 4.0 * g_xz_y_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_y_yy[i] = 4.0 * g_xz_y_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_y_yz[i] = 4.0 * g_xz_y_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_y_zz[i] = 4.0 * g_xz_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_xz_0_0_0_0_y_z_xx, g_xz_0_0_0_0_y_z_xy, g_xz_0_0_0_0_y_z_xz, g_xz_0_0_0_0_y_z_yy, g_xz_0_0_0_0_y_z_yz, g_xz_0_0_0_0_y_z_zz, g_xz_y_z_xx, g_xz_y_z_xy, g_xz_y_z_xz, g_xz_y_z_yy, g_xz_y_z_yz, g_xz_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_y_z_xx[i] = 4.0 * g_xz_y_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_z_xy[i] = 4.0 * g_xz_y_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_z_xz[i] = 4.0 * g_xz_y_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_z_yy[i] = 4.0 * g_xz_y_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_z_yz[i] = 4.0 * g_xz_y_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_z_zz[i] = 4.0 * g_xz_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_xz_0_0_0_0_z_x_xx, g_xz_0_0_0_0_z_x_xy, g_xz_0_0_0_0_z_x_xz, g_xz_0_0_0_0_z_x_yy, g_xz_0_0_0_0_z_x_yz, g_xz_0_0_0_0_z_x_zz, g_xz_z_x_xx, g_xz_z_x_xy, g_xz_z_x_xz, g_xz_z_x_yy, g_xz_z_x_yz, g_xz_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_z_x_xx[i] = 4.0 * g_xz_z_x_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_x_xy[i] = 4.0 * g_xz_z_x_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_x_xz[i] = 4.0 * g_xz_z_x_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_x_yy[i] = 4.0 * g_xz_z_x_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_x_yz[i] = 4.0 * g_xz_z_x_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_x_zz[i] = 4.0 * g_xz_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_xz_0_0_0_0_z_y_xx, g_xz_0_0_0_0_z_y_xy, g_xz_0_0_0_0_z_y_xz, g_xz_0_0_0_0_z_y_yy, g_xz_0_0_0_0_z_y_yz, g_xz_0_0_0_0_z_y_zz, g_xz_z_y_xx, g_xz_z_y_xy, g_xz_z_y_xz, g_xz_z_y_yy, g_xz_z_y_yz, g_xz_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_z_y_xx[i] = 4.0 * g_xz_z_y_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_y_xy[i] = 4.0 * g_xz_z_y_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_y_xz[i] = 4.0 * g_xz_z_y_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_y_yy[i] = 4.0 * g_xz_z_y_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_y_yz[i] = 4.0 * g_xz_z_y_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_y_zz[i] = 4.0 * g_xz_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_xz_0_0_0_0_z_z_xx, g_xz_0_0_0_0_z_z_xy, g_xz_0_0_0_0_z_z_xz, g_xz_0_0_0_0_z_z_yy, g_xz_0_0_0_0_z_z_yz, g_xz_0_0_0_0_z_z_zz, g_xz_z_z_xx, g_xz_z_z_xy, g_xz_z_z_xz, g_xz_z_z_yy, g_xz_z_z_yz, g_xz_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_z_z_xx[i] = 4.0 * g_xz_z_z_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_z_xy[i] = 4.0 * g_xz_z_z_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_z_xz[i] = 4.0 * g_xz_z_z_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_z_yy[i] = 4.0 * g_xz_z_z_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_z_yz[i] = 4.0 * g_xz_z_z_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_z_zz[i] = 4.0 * g_xz_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_yy_0_0_0_0_x_x_xx, g_yy_0_0_0_0_x_x_xy, g_yy_0_0_0_0_x_x_xz, g_yy_0_0_0_0_x_x_yy, g_yy_0_0_0_0_x_x_yz, g_yy_0_0_0_0_x_x_zz, g_yy_x_x_xx, g_yy_x_x_xy, g_yy_x_x_xz, g_yy_x_x_yy, g_yy_x_x_yz, g_yy_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_x_x_xx[i] = -2.0 * g_0_x_x_xx[i] * a_exp + 4.0 * g_yy_x_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_x_xy[i] = -2.0 * g_0_x_x_xy[i] * a_exp + 4.0 * g_yy_x_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_x_xz[i] = -2.0 * g_0_x_x_xz[i] * a_exp + 4.0 * g_yy_x_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_x_yy[i] = -2.0 * g_0_x_x_yy[i] * a_exp + 4.0 * g_yy_x_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_x_yz[i] = -2.0 * g_0_x_x_yz[i] * a_exp + 4.0 * g_yy_x_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_x_zz[i] = -2.0 * g_0_x_x_zz[i] * a_exp + 4.0 * g_yy_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_yy_0_0_0_0_x_y_xx, g_yy_0_0_0_0_x_y_xy, g_yy_0_0_0_0_x_y_xz, g_yy_0_0_0_0_x_y_yy, g_yy_0_0_0_0_x_y_yz, g_yy_0_0_0_0_x_y_zz, g_yy_x_y_xx, g_yy_x_y_xy, g_yy_x_y_xz, g_yy_x_y_yy, g_yy_x_y_yz, g_yy_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_x_y_xx[i] = -2.0 * g_0_x_y_xx[i] * a_exp + 4.0 * g_yy_x_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_y_xy[i] = -2.0 * g_0_x_y_xy[i] * a_exp + 4.0 * g_yy_x_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_y_xz[i] = -2.0 * g_0_x_y_xz[i] * a_exp + 4.0 * g_yy_x_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_y_yy[i] = -2.0 * g_0_x_y_yy[i] * a_exp + 4.0 * g_yy_x_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_y_yz[i] = -2.0 * g_0_x_y_yz[i] * a_exp + 4.0 * g_yy_x_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_y_zz[i] = -2.0 * g_0_x_y_zz[i] * a_exp + 4.0 * g_yy_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_yy_0_0_0_0_x_z_xx, g_yy_0_0_0_0_x_z_xy, g_yy_0_0_0_0_x_z_xz, g_yy_0_0_0_0_x_z_yy, g_yy_0_0_0_0_x_z_yz, g_yy_0_0_0_0_x_z_zz, g_yy_x_z_xx, g_yy_x_z_xy, g_yy_x_z_xz, g_yy_x_z_yy, g_yy_x_z_yz, g_yy_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_x_z_xx[i] = -2.0 * g_0_x_z_xx[i] * a_exp + 4.0 * g_yy_x_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_z_xy[i] = -2.0 * g_0_x_z_xy[i] * a_exp + 4.0 * g_yy_x_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_z_xz[i] = -2.0 * g_0_x_z_xz[i] * a_exp + 4.0 * g_yy_x_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_z_yy[i] = -2.0 * g_0_x_z_yy[i] * a_exp + 4.0 * g_yy_x_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_z_yz[i] = -2.0 * g_0_x_z_yz[i] * a_exp + 4.0 * g_yy_x_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_z_zz[i] = -2.0 * g_0_x_z_zz[i] * a_exp + 4.0 * g_yy_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_yy_0_0_0_0_y_x_xx, g_yy_0_0_0_0_y_x_xy, g_yy_0_0_0_0_y_x_xz, g_yy_0_0_0_0_y_x_yy, g_yy_0_0_0_0_y_x_yz, g_yy_0_0_0_0_y_x_zz, g_yy_y_x_xx, g_yy_y_x_xy, g_yy_y_x_xz, g_yy_y_x_yy, g_yy_y_x_yz, g_yy_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_y_x_xx[i] = -2.0 * g_0_y_x_xx[i] * a_exp + 4.0 * g_yy_y_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_x_xy[i] = -2.0 * g_0_y_x_xy[i] * a_exp + 4.0 * g_yy_y_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_x_xz[i] = -2.0 * g_0_y_x_xz[i] * a_exp + 4.0 * g_yy_y_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_x_yy[i] = -2.0 * g_0_y_x_yy[i] * a_exp + 4.0 * g_yy_y_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_x_yz[i] = -2.0 * g_0_y_x_yz[i] * a_exp + 4.0 * g_yy_y_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_x_zz[i] = -2.0 * g_0_y_x_zz[i] * a_exp + 4.0 * g_yy_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_yy_0_0_0_0_y_y_xx, g_yy_0_0_0_0_y_y_xy, g_yy_0_0_0_0_y_y_xz, g_yy_0_0_0_0_y_y_yy, g_yy_0_0_0_0_y_y_yz, g_yy_0_0_0_0_y_y_zz, g_yy_y_y_xx, g_yy_y_y_xy, g_yy_y_y_xz, g_yy_y_y_yy, g_yy_y_y_yz, g_yy_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_y_y_xx[i] = -2.0 * g_0_y_y_xx[i] * a_exp + 4.0 * g_yy_y_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_y_xy[i] = -2.0 * g_0_y_y_xy[i] * a_exp + 4.0 * g_yy_y_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_y_xz[i] = -2.0 * g_0_y_y_xz[i] * a_exp + 4.0 * g_yy_y_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_y_yy[i] = -2.0 * g_0_y_y_yy[i] * a_exp + 4.0 * g_yy_y_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_y_yz[i] = -2.0 * g_0_y_y_yz[i] * a_exp + 4.0 * g_yy_y_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_y_zz[i] = -2.0 * g_0_y_y_zz[i] * a_exp + 4.0 * g_yy_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_yy_0_0_0_0_y_z_xx, g_yy_0_0_0_0_y_z_xy, g_yy_0_0_0_0_y_z_xz, g_yy_0_0_0_0_y_z_yy, g_yy_0_0_0_0_y_z_yz, g_yy_0_0_0_0_y_z_zz, g_yy_y_z_xx, g_yy_y_z_xy, g_yy_y_z_xz, g_yy_y_z_yy, g_yy_y_z_yz, g_yy_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_y_z_xx[i] = -2.0 * g_0_y_z_xx[i] * a_exp + 4.0 * g_yy_y_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_z_xy[i] = -2.0 * g_0_y_z_xy[i] * a_exp + 4.0 * g_yy_y_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_z_xz[i] = -2.0 * g_0_y_z_xz[i] * a_exp + 4.0 * g_yy_y_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_z_yy[i] = -2.0 * g_0_y_z_yy[i] * a_exp + 4.0 * g_yy_y_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_z_yz[i] = -2.0 * g_0_y_z_yz[i] * a_exp + 4.0 * g_yy_y_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_z_zz[i] = -2.0 * g_0_y_z_zz[i] * a_exp + 4.0 * g_yy_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_yy_0_0_0_0_z_x_xx, g_yy_0_0_0_0_z_x_xy, g_yy_0_0_0_0_z_x_xz, g_yy_0_0_0_0_z_x_yy, g_yy_0_0_0_0_z_x_yz, g_yy_0_0_0_0_z_x_zz, g_yy_z_x_xx, g_yy_z_x_xy, g_yy_z_x_xz, g_yy_z_x_yy, g_yy_z_x_yz, g_yy_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_z_x_xx[i] = -2.0 * g_0_z_x_xx[i] * a_exp + 4.0 * g_yy_z_x_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_x_xy[i] = -2.0 * g_0_z_x_xy[i] * a_exp + 4.0 * g_yy_z_x_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_x_xz[i] = -2.0 * g_0_z_x_xz[i] * a_exp + 4.0 * g_yy_z_x_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_x_yy[i] = -2.0 * g_0_z_x_yy[i] * a_exp + 4.0 * g_yy_z_x_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_x_yz[i] = -2.0 * g_0_z_x_yz[i] * a_exp + 4.0 * g_yy_z_x_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_x_zz[i] = -2.0 * g_0_z_x_zz[i] * a_exp + 4.0 * g_yy_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_yy_0_0_0_0_z_y_xx, g_yy_0_0_0_0_z_y_xy, g_yy_0_0_0_0_z_y_xz, g_yy_0_0_0_0_z_y_yy, g_yy_0_0_0_0_z_y_yz, g_yy_0_0_0_0_z_y_zz, g_yy_z_y_xx, g_yy_z_y_xy, g_yy_z_y_xz, g_yy_z_y_yy, g_yy_z_y_yz, g_yy_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_z_y_xx[i] = -2.0 * g_0_z_y_xx[i] * a_exp + 4.0 * g_yy_z_y_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_y_xy[i] = -2.0 * g_0_z_y_xy[i] * a_exp + 4.0 * g_yy_z_y_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_y_xz[i] = -2.0 * g_0_z_y_xz[i] * a_exp + 4.0 * g_yy_z_y_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_y_yy[i] = -2.0 * g_0_z_y_yy[i] * a_exp + 4.0 * g_yy_z_y_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_y_yz[i] = -2.0 * g_0_z_y_yz[i] * a_exp + 4.0 * g_yy_z_y_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_y_zz[i] = -2.0 * g_0_z_y_zz[i] * a_exp + 4.0 * g_yy_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_yy_0_0_0_0_z_z_xx, g_yy_0_0_0_0_z_z_xy, g_yy_0_0_0_0_z_z_xz, g_yy_0_0_0_0_z_z_yy, g_yy_0_0_0_0_z_z_yz, g_yy_0_0_0_0_z_z_zz, g_yy_z_z_xx, g_yy_z_z_xy, g_yy_z_z_xz, g_yy_z_z_yy, g_yy_z_z_yz, g_yy_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_z_z_xx[i] = -2.0 * g_0_z_z_xx[i] * a_exp + 4.0 * g_yy_z_z_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_z_xy[i] = -2.0 * g_0_z_z_xy[i] * a_exp + 4.0 * g_yy_z_z_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_z_xz[i] = -2.0 * g_0_z_z_xz[i] * a_exp + 4.0 * g_yy_z_z_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_z_yy[i] = -2.0 * g_0_z_z_yy[i] * a_exp + 4.0 * g_yy_z_z_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_z_yz[i] = -2.0 * g_0_z_z_yz[i] * a_exp + 4.0 * g_yy_z_z_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_z_zz[i] = -2.0 * g_0_z_z_zz[i] * a_exp + 4.0 * g_yy_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_yz_0_0_0_0_x_x_xx, g_yz_0_0_0_0_x_x_xy, g_yz_0_0_0_0_x_x_xz, g_yz_0_0_0_0_x_x_yy, g_yz_0_0_0_0_x_x_yz, g_yz_0_0_0_0_x_x_zz, g_yz_x_x_xx, g_yz_x_x_xy, g_yz_x_x_xz, g_yz_x_x_yy, g_yz_x_x_yz, g_yz_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_x_x_xx[i] = 4.0 * g_yz_x_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_x_xy[i] = 4.0 * g_yz_x_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_x_xz[i] = 4.0 * g_yz_x_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_x_yy[i] = 4.0 * g_yz_x_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_x_yz[i] = 4.0 * g_yz_x_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_x_zz[i] = 4.0 * g_yz_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_yz_0_0_0_0_x_y_xx, g_yz_0_0_0_0_x_y_xy, g_yz_0_0_0_0_x_y_xz, g_yz_0_0_0_0_x_y_yy, g_yz_0_0_0_0_x_y_yz, g_yz_0_0_0_0_x_y_zz, g_yz_x_y_xx, g_yz_x_y_xy, g_yz_x_y_xz, g_yz_x_y_yy, g_yz_x_y_yz, g_yz_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_x_y_xx[i] = 4.0 * g_yz_x_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_y_xy[i] = 4.0 * g_yz_x_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_y_xz[i] = 4.0 * g_yz_x_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_y_yy[i] = 4.0 * g_yz_x_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_y_yz[i] = 4.0 * g_yz_x_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_y_zz[i] = 4.0 * g_yz_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_yz_0_0_0_0_x_z_xx, g_yz_0_0_0_0_x_z_xy, g_yz_0_0_0_0_x_z_xz, g_yz_0_0_0_0_x_z_yy, g_yz_0_0_0_0_x_z_yz, g_yz_0_0_0_0_x_z_zz, g_yz_x_z_xx, g_yz_x_z_xy, g_yz_x_z_xz, g_yz_x_z_yy, g_yz_x_z_yz, g_yz_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_x_z_xx[i] = 4.0 * g_yz_x_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_z_xy[i] = 4.0 * g_yz_x_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_z_xz[i] = 4.0 * g_yz_x_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_z_yy[i] = 4.0 * g_yz_x_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_z_yz[i] = 4.0 * g_yz_x_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_z_zz[i] = 4.0 * g_yz_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_yz_0_0_0_0_y_x_xx, g_yz_0_0_0_0_y_x_xy, g_yz_0_0_0_0_y_x_xz, g_yz_0_0_0_0_y_x_yy, g_yz_0_0_0_0_y_x_yz, g_yz_0_0_0_0_y_x_zz, g_yz_y_x_xx, g_yz_y_x_xy, g_yz_y_x_xz, g_yz_y_x_yy, g_yz_y_x_yz, g_yz_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_y_x_xx[i] = 4.0 * g_yz_y_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_x_xy[i] = 4.0 * g_yz_y_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_x_xz[i] = 4.0 * g_yz_y_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_x_yy[i] = 4.0 * g_yz_y_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_x_yz[i] = 4.0 * g_yz_y_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_x_zz[i] = 4.0 * g_yz_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_yz_0_0_0_0_y_y_xx, g_yz_0_0_0_0_y_y_xy, g_yz_0_0_0_0_y_y_xz, g_yz_0_0_0_0_y_y_yy, g_yz_0_0_0_0_y_y_yz, g_yz_0_0_0_0_y_y_zz, g_yz_y_y_xx, g_yz_y_y_xy, g_yz_y_y_xz, g_yz_y_y_yy, g_yz_y_y_yz, g_yz_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_y_y_xx[i] = 4.0 * g_yz_y_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_y_xy[i] = 4.0 * g_yz_y_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_y_xz[i] = 4.0 * g_yz_y_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_y_yy[i] = 4.0 * g_yz_y_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_y_yz[i] = 4.0 * g_yz_y_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_y_zz[i] = 4.0 * g_yz_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_yz_0_0_0_0_y_z_xx, g_yz_0_0_0_0_y_z_xy, g_yz_0_0_0_0_y_z_xz, g_yz_0_0_0_0_y_z_yy, g_yz_0_0_0_0_y_z_yz, g_yz_0_0_0_0_y_z_zz, g_yz_y_z_xx, g_yz_y_z_xy, g_yz_y_z_xz, g_yz_y_z_yy, g_yz_y_z_yz, g_yz_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_y_z_xx[i] = 4.0 * g_yz_y_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_z_xy[i] = 4.0 * g_yz_y_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_z_xz[i] = 4.0 * g_yz_y_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_z_yy[i] = 4.0 * g_yz_y_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_z_yz[i] = 4.0 * g_yz_y_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_z_zz[i] = 4.0 * g_yz_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_yz_0_0_0_0_z_x_xx, g_yz_0_0_0_0_z_x_xy, g_yz_0_0_0_0_z_x_xz, g_yz_0_0_0_0_z_x_yy, g_yz_0_0_0_0_z_x_yz, g_yz_0_0_0_0_z_x_zz, g_yz_z_x_xx, g_yz_z_x_xy, g_yz_z_x_xz, g_yz_z_x_yy, g_yz_z_x_yz, g_yz_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_z_x_xx[i] = 4.0 * g_yz_z_x_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_x_xy[i] = 4.0 * g_yz_z_x_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_x_xz[i] = 4.0 * g_yz_z_x_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_x_yy[i] = 4.0 * g_yz_z_x_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_x_yz[i] = 4.0 * g_yz_z_x_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_x_zz[i] = 4.0 * g_yz_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_yz_0_0_0_0_z_y_xx, g_yz_0_0_0_0_z_y_xy, g_yz_0_0_0_0_z_y_xz, g_yz_0_0_0_0_z_y_yy, g_yz_0_0_0_0_z_y_yz, g_yz_0_0_0_0_z_y_zz, g_yz_z_y_xx, g_yz_z_y_xy, g_yz_z_y_xz, g_yz_z_y_yy, g_yz_z_y_yz, g_yz_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_z_y_xx[i] = 4.0 * g_yz_z_y_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_y_xy[i] = 4.0 * g_yz_z_y_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_y_xz[i] = 4.0 * g_yz_z_y_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_y_yy[i] = 4.0 * g_yz_z_y_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_y_yz[i] = 4.0 * g_yz_z_y_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_y_zz[i] = 4.0 * g_yz_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_yz_0_0_0_0_z_z_xx, g_yz_0_0_0_0_z_z_xy, g_yz_0_0_0_0_z_z_xz, g_yz_0_0_0_0_z_z_yy, g_yz_0_0_0_0_z_z_yz, g_yz_0_0_0_0_z_z_zz, g_yz_z_z_xx, g_yz_z_z_xy, g_yz_z_z_xz, g_yz_z_z_yy, g_yz_z_z_yz, g_yz_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_z_z_xx[i] = 4.0 * g_yz_z_z_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_z_xy[i] = 4.0 * g_yz_z_z_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_z_xz[i] = 4.0 * g_yz_z_z_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_z_yy[i] = 4.0 * g_yz_z_z_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_z_yz[i] = 4.0 * g_yz_z_z_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_z_zz[i] = 4.0 * g_yz_z_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_zz_0_0_0_0_x_x_xx, g_zz_0_0_0_0_x_x_xy, g_zz_0_0_0_0_x_x_xz, g_zz_0_0_0_0_x_x_yy, g_zz_0_0_0_0_x_x_yz, g_zz_0_0_0_0_x_x_zz, g_zz_x_x_xx, g_zz_x_x_xy, g_zz_x_x_xz, g_zz_x_x_yy, g_zz_x_x_yz, g_zz_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_x_x_xx[i] = -2.0 * g_0_x_x_xx[i] * a_exp + 4.0 * g_zz_x_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_x_xy[i] = -2.0 * g_0_x_x_xy[i] * a_exp + 4.0 * g_zz_x_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_x_xz[i] = -2.0 * g_0_x_x_xz[i] * a_exp + 4.0 * g_zz_x_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_x_yy[i] = -2.0 * g_0_x_x_yy[i] * a_exp + 4.0 * g_zz_x_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_x_yz[i] = -2.0 * g_0_x_x_yz[i] * a_exp + 4.0 * g_zz_x_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_x_zz[i] = -2.0 * g_0_x_x_zz[i] * a_exp + 4.0 * g_zz_x_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_zz_0_0_0_0_x_y_xx, g_zz_0_0_0_0_x_y_xy, g_zz_0_0_0_0_x_y_xz, g_zz_0_0_0_0_x_y_yy, g_zz_0_0_0_0_x_y_yz, g_zz_0_0_0_0_x_y_zz, g_zz_x_y_xx, g_zz_x_y_xy, g_zz_x_y_xz, g_zz_x_y_yy, g_zz_x_y_yz, g_zz_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_x_y_xx[i] = -2.0 * g_0_x_y_xx[i] * a_exp + 4.0 * g_zz_x_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_y_xy[i] = -2.0 * g_0_x_y_xy[i] * a_exp + 4.0 * g_zz_x_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_y_xz[i] = -2.0 * g_0_x_y_xz[i] * a_exp + 4.0 * g_zz_x_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_y_yy[i] = -2.0 * g_0_x_y_yy[i] * a_exp + 4.0 * g_zz_x_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_y_yz[i] = -2.0 * g_0_x_y_yz[i] * a_exp + 4.0 * g_zz_x_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_y_zz[i] = -2.0 * g_0_x_y_zz[i] * a_exp + 4.0 * g_zz_x_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_zz_0_0_0_0_x_z_xx, g_zz_0_0_0_0_x_z_xy, g_zz_0_0_0_0_x_z_xz, g_zz_0_0_0_0_x_z_yy, g_zz_0_0_0_0_x_z_yz, g_zz_0_0_0_0_x_z_zz, g_zz_x_z_xx, g_zz_x_z_xy, g_zz_x_z_xz, g_zz_x_z_yy, g_zz_x_z_yz, g_zz_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_x_z_xx[i] = -2.0 * g_0_x_z_xx[i] * a_exp + 4.0 * g_zz_x_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_z_xy[i] = -2.0 * g_0_x_z_xy[i] * a_exp + 4.0 * g_zz_x_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_z_xz[i] = -2.0 * g_0_x_z_xz[i] * a_exp + 4.0 * g_zz_x_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_z_yy[i] = -2.0 * g_0_x_z_yy[i] * a_exp + 4.0 * g_zz_x_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_z_yz[i] = -2.0 * g_0_x_z_yz[i] * a_exp + 4.0 * g_zz_x_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_z_zz[i] = -2.0 * g_0_x_z_zz[i] * a_exp + 4.0 * g_zz_x_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_zz_0_0_0_0_y_x_xx, g_zz_0_0_0_0_y_x_xy, g_zz_0_0_0_0_y_x_xz, g_zz_0_0_0_0_y_x_yy, g_zz_0_0_0_0_y_x_yz, g_zz_0_0_0_0_y_x_zz, g_zz_y_x_xx, g_zz_y_x_xy, g_zz_y_x_xz, g_zz_y_x_yy, g_zz_y_x_yz, g_zz_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_y_x_xx[i] = -2.0 * g_0_y_x_xx[i] * a_exp + 4.0 * g_zz_y_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_x_xy[i] = -2.0 * g_0_y_x_xy[i] * a_exp + 4.0 * g_zz_y_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_x_xz[i] = -2.0 * g_0_y_x_xz[i] * a_exp + 4.0 * g_zz_y_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_x_yy[i] = -2.0 * g_0_y_x_yy[i] * a_exp + 4.0 * g_zz_y_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_x_yz[i] = -2.0 * g_0_y_x_yz[i] * a_exp + 4.0 * g_zz_y_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_x_zz[i] = -2.0 * g_0_y_x_zz[i] * a_exp + 4.0 * g_zz_y_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_zz_0_0_0_0_y_y_xx, g_zz_0_0_0_0_y_y_xy, g_zz_0_0_0_0_y_y_xz, g_zz_0_0_0_0_y_y_yy, g_zz_0_0_0_0_y_y_yz, g_zz_0_0_0_0_y_y_zz, g_zz_y_y_xx, g_zz_y_y_xy, g_zz_y_y_xz, g_zz_y_y_yy, g_zz_y_y_yz, g_zz_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_y_y_xx[i] = -2.0 * g_0_y_y_xx[i] * a_exp + 4.0 * g_zz_y_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_y_xy[i] = -2.0 * g_0_y_y_xy[i] * a_exp + 4.0 * g_zz_y_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_y_xz[i] = -2.0 * g_0_y_y_xz[i] * a_exp + 4.0 * g_zz_y_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_y_yy[i] = -2.0 * g_0_y_y_yy[i] * a_exp + 4.0 * g_zz_y_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_y_yz[i] = -2.0 * g_0_y_y_yz[i] * a_exp + 4.0 * g_zz_y_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_y_zz[i] = -2.0 * g_0_y_y_zz[i] * a_exp + 4.0 * g_zz_y_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_zz_0_0_0_0_y_z_xx, g_zz_0_0_0_0_y_z_xy, g_zz_0_0_0_0_y_z_xz, g_zz_0_0_0_0_y_z_yy, g_zz_0_0_0_0_y_z_yz, g_zz_0_0_0_0_y_z_zz, g_zz_y_z_xx, g_zz_y_z_xy, g_zz_y_z_xz, g_zz_y_z_yy, g_zz_y_z_yz, g_zz_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_y_z_xx[i] = -2.0 * g_0_y_z_xx[i] * a_exp + 4.0 * g_zz_y_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_z_xy[i] = -2.0 * g_0_y_z_xy[i] * a_exp + 4.0 * g_zz_y_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_z_xz[i] = -2.0 * g_0_y_z_xz[i] * a_exp + 4.0 * g_zz_y_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_z_yy[i] = -2.0 * g_0_y_z_yy[i] * a_exp + 4.0 * g_zz_y_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_z_yz[i] = -2.0 * g_0_y_z_yz[i] * a_exp + 4.0 * g_zz_y_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_z_zz[i] = -2.0 * g_0_y_z_zz[i] * a_exp + 4.0 * g_zz_y_z_zz[i] * a_exp * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_zz_0_0_0_0_z_x_xx, g_zz_0_0_0_0_z_x_xy, g_zz_0_0_0_0_z_x_xz, g_zz_0_0_0_0_z_x_yy, g_zz_0_0_0_0_z_x_yz, g_zz_0_0_0_0_z_x_zz, g_zz_z_x_xx, g_zz_z_x_xy, g_zz_z_x_xz, g_zz_z_x_yy, g_zz_z_x_yz, g_zz_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_z_x_xx[i] = -2.0 * g_0_z_x_xx[i] * a_exp + 4.0 * g_zz_z_x_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_x_xy[i] = -2.0 * g_0_z_x_xy[i] * a_exp + 4.0 * g_zz_z_x_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_x_xz[i] = -2.0 * g_0_z_x_xz[i] * a_exp + 4.0 * g_zz_z_x_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_x_yy[i] = -2.0 * g_0_z_x_yy[i] * a_exp + 4.0 * g_zz_z_x_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_x_yz[i] = -2.0 * g_0_z_x_yz[i] * a_exp + 4.0 * g_zz_z_x_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_x_zz[i] = -2.0 * g_0_z_x_zz[i] * a_exp + 4.0 * g_zz_z_x_zz[i] * a_exp * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_zz_0_0_0_0_z_y_xx, g_zz_0_0_0_0_z_y_xy, g_zz_0_0_0_0_z_y_xz, g_zz_0_0_0_0_z_y_yy, g_zz_0_0_0_0_z_y_yz, g_zz_0_0_0_0_z_y_zz, g_zz_z_y_xx, g_zz_z_y_xy, g_zz_z_y_xz, g_zz_z_y_yy, g_zz_z_y_yz, g_zz_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_z_y_xx[i] = -2.0 * g_0_z_y_xx[i] * a_exp + 4.0 * g_zz_z_y_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_y_xy[i] = -2.0 * g_0_z_y_xy[i] * a_exp + 4.0 * g_zz_z_y_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_y_xz[i] = -2.0 * g_0_z_y_xz[i] * a_exp + 4.0 * g_zz_z_y_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_y_yy[i] = -2.0 * g_0_z_y_yy[i] * a_exp + 4.0 * g_zz_z_y_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_y_yz[i] = -2.0 * g_0_z_y_yz[i] * a_exp + 4.0 * g_zz_z_y_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_y_zz[i] = -2.0 * g_0_z_y_zz[i] * a_exp + 4.0 * g_zz_z_y_zz[i] * a_exp * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_zz_0_0_0_0_z_z_xx, g_zz_0_0_0_0_z_z_xy, g_zz_0_0_0_0_z_z_xz, g_zz_0_0_0_0_z_z_yy, g_zz_0_0_0_0_z_z_yz, g_zz_0_0_0_0_z_z_zz, g_zz_z_z_xx, g_zz_z_z_xy, g_zz_z_z_xz, g_zz_z_z_yy, g_zz_z_z_yz, g_zz_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_z_z_xx[i] = -2.0 * g_0_z_z_xx[i] * a_exp + 4.0 * g_zz_z_z_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_z_xy[i] = -2.0 * g_0_z_z_xy[i] * a_exp + 4.0 * g_zz_z_z_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_z_xz[i] = -2.0 * g_0_z_z_xz[i] * a_exp + 4.0 * g_zz_z_z_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_z_yy[i] = -2.0 * g_0_z_z_yy[i] * a_exp + 4.0 * g_zz_z_z_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_z_yz[i] = -2.0 * g_0_z_z_yz[i] * a_exp + 4.0 * g_zz_z_z_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_z_zz[i] = -2.0 * g_0_z_z_zz[i] * a_exp + 4.0 * g_zz_z_z_zz[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

