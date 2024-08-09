#include "GeomDeriv1000OfScalarForPPPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_pppd_0(CSimdArray<double>& buffer_1000_pppd,
                     const CSimdArray<double>& buffer_sppd,
                     const CSimdArray<double>& buffer_dppd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_pppd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_pppd

    auto g_x_0_0_0_x_x_x_xx = buffer_1000_pppd[0];

    auto g_x_0_0_0_x_x_x_xy = buffer_1000_pppd[1];

    auto g_x_0_0_0_x_x_x_xz = buffer_1000_pppd[2];

    auto g_x_0_0_0_x_x_x_yy = buffer_1000_pppd[3];

    auto g_x_0_0_0_x_x_x_yz = buffer_1000_pppd[4];

    auto g_x_0_0_0_x_x_x_zz = buffer_1000_pppd[5];

    auto g_x_0_0_0_x_x_y_xx = buffer_1000_pppd[6];

    auto g_x_0_0_0_x_x_y_xy = buffer_1000_pppd[7];

    auto g_x_0_0_0_x_x_y_xz = buffer_1000_pppd[8];

    auto g_x_0_0_0_x_x_y_yy = buffer_1000_pppd[9];

    auto g_x_0_0_0_x_x_y_yz = buffer_1000_pppd[10];

    auto g_x_0_0_0_x_x_y_zz = buffer_1000_pppd[11];

    auto g_x_0_0_0_x_x_z_xx = buffer_1000_pppd[12];

    auto g_x_0_0_0_x_x_z_xy = buffer_1000_pppd[13];

    auto g_x_0_0_0_x_x_z_xz = buffer_1000_pppd[14];

    auto g_x_0_0_0_x_x_z_yy = buffer_1000_pppd[15];

    auto g_x_0_0_0_x_x_z_yz = buffer_1000_pppd[16];

    auto g_x_0_0_0_x_x_z_zz = buffer_1000_pppd[17];

    auto g_x_0_0_0_x_y_x_xx = buffer_1000_pppd[18];

    auto g_x_0_0_0_x_y_x_xy = buffer_1000_pppd[19];

    auto g_x_0_0_0_x_y_x_xz = buffer_1000_pppd[20];

    auto g_x_0_0_0_x_y_x_yy = buffer_1000_pppd[21];

    auto g_x_0_0_0_x_y_x_yz = buffer_1000_pppd[22];

    auto g_x_0_0_0_x_y_x_zz = buffer_1000_pppd[23];

    auto g_x_0_0_0_x_y_y_xx = buffer_1000_pppd[24];

    auto g_x_0_0_0_x_y_y_xy = buffer_1000_pppd[25];

    auto g_x_0_0_0_x_y_y_xz = buffer_1000_pppd[26];

    auto g_x_0_0_0_x_y_y_yy = buffer_1000_pppd[27];

    auto g_x_0_0_0_x_y_y_yz = buffer_1000_pppd[28];

    auto g_x_0_0_0_x_y_y_zz = buffer_1000_pppd[29];

    auto g_x_0_0_0_x_y_z_xx = buffer_1000_pppd[30];

    auto g_x_0_0_0_x_y_z_xy = buffer_1000_pppd[31];

    auto g_x_0_0_0_x_y_z_xz = buffer_1000_pppd[32];

    auto g_x_0_0_0_x_y_z_yy = buffer_1000_pppd[33];

    auto g_x_0_0_0_x_y_z_yz = buffer_1000_pppd[34];

    auto g_x_0_0_0_x_y_z_zz = buffer_1000_pppd[35];

    auto g_x_0_0_0_x_z_x_xx = buffer_1000_pppd[36];

    auto g_x_0_0_0_x_z_x_xy = buffer_1000_pppd[37];

    auto g_x_0_0_0_x_z_x_xz = buffer_1000_pppd[38];

    auto g_x_0_0_0_x_z_x_yy = buffer_1000_pppd[39];

    auto g_x_0_0_0_x_z_x_yz = buffer_1000_pppd[40];

    auto g_x_0_0_0_x_z_x_zz = buffer_1000_pppd[41];

    auto g_x_0_0_0_x_z_y_xx = buffer_1000_pppd[42];

    auto g_x_0_0_0_x_z_y_xy = buffer_1000_pppd[43];

    auto g_x_0_0_0_x_z_y_xz = buffer_1000_pppd[44];

    auto g_x_0_0_0_x_z_y_yy = buffer_1000_pppd[45];

    auto g_x_0_0_0_x_z_y_yz = buffer_1000_pppd[46];

    auto g_x_0_0_0_x_z_y_zz = buffer_1000_pppd[47];

    auto g_x_0_0_0_x_z_z_xx = buffer_1000_pppd[48];

    auto g_x_0_0_0_x_z_z_xy = buffer_1000_pppd[49];

    auto g_x_0_0_0_x_z_z_xz = buffer_1000_pppd[50];

    auto g_x_0_0_0_x_z_z_yy = buffer_1000_pppd[51];

    auto g_x_0_0_0_x_z_z_yz = buffer_1000_pppd[52];

    auto g_x_0_0_0_x_z_z_zz = buffer_1000_pppd[53];

    auto g_x_0_0_0_y_x_x_xx = buffer_1000_pppd[54];

    auto g_x_0_0_0_y_x_x_xy = buffer_1000_pppd[55];

    auto g_x_0_0_0_y_x_x_xz = buffer_1000_pppd[56];

    auto g_x_0_0_0_y_x_x_yy = buffer_1000_pppd[57];

    auto g_x_0_0_0_y_x_x_yz = buffer_1000_pppd[58];

    auto g_x_0_0_0_y_x_x_zz = buffer_1000_pppd[59];

    auto g_x_0_0_0_y_x_y_xx = buffer_1000_pppd[60];

    auto g_x_0_0_0_y_x_y_xy = buffer_1000_pppd[61];

    auto g_x_0_0_0_y_x_y_xz = buffer_1000_pppd[62];

    auto g_x_0_0_0_y_x_y_yy = buffer_1000_pppd[63];

    auto g_x_0_0_0_y_x_y_yz = buffer_1000_pppd[64];

    auto g_x_0_0_0_y_x_y_zz = buffer_1000_pppd[65];

    auto g_x_0_0_0_y_x_z_xx = buffer_1000_pppd[66];

    auto g_x_0_0_0_y_x_z_xy = buffer_1000_pppd[67];

    auto g_x_0_0_0_y_x_z_xz = buffer_1000_pppd[68];

    auto g_x_0_0_0_y_x_z_yy = buffer_1000_pppd[69];

    auto g_x_0_0_0_y_x_z_yz = buffer_1000_pppd[70];

    auto g_x_0_0_0_y_x_z_zz = buffer_1000_pppd[71];

    auto g_x_0_0_0_y_y_x_xx = buffer_1000_pppd[72];

    auto g_x_0_0_0_y_y_x_xy = buffer_1000_pppd[73];

    auto g_x_0_0_0_y_y_x_xz = buffer_1000_pppd[74];

    auto g_x_0_0_0_y_y_x_yy = buffer_1000_pppd[75];

    auto g_x_0_0_0_y_y_x_yz = buffer_1000_pppd[76];

    auto g_x_0_0_0_y_y_x_zz = buffer_1000_pppd[77];

    auto g_x_0_0_0_y_y_y_xx = buffer_1000_pppd[78];

    auto g_x_0_0_0_y_y_y_xy = buffer_1000_pppd[79];

    auto g_x_0_0_0_y_y_y_xz = buffer_1000_pppd[80];

    auto g_x_0_0_0_y_y_y_yy = buffer_1000_pppd[81];

    auto g_x_0_0_0_y_y_y_yz = buffer_1000_pppd[82];

    auto g_x_0_0_0_y_y_y_zz = buffer_1000_pppd[83];

    auto g_x_0_0_0_y_y_z_xx = buffer_1000_pppd[84];

    auto g_x_0_0_0_y_y_z_xy = buffer_1000_pppd[85];

    auto g_x_0_0_0_y_y_z_xz = buffer_1000_pppd[86];

    auto g_x_0_0_0_y_y_z_yy = buffer_1000_pppd[87];

    auto g_x_0_0_0_y_y_z_yz = buffer_1000_pppd[88];

    auto g_x_0_0_0_y_y_z_zz = buffer_1000_pppd[89];

    auto g_x_0_0_0_y_z_x_xx = buffer_1000_pppd[90];

    auto g_x_0_0_0_y_z_x_xy = buffer_1000_pppd[91];

    auto g_x_0_0_0_y_z_x_xz = buffer_1000_pppd[92];

    auto g_x_0_0_0_y_z_x_yy = buffer_1000_pppd[93];

    auto g_x_0_0_0_y_z_x_yz = buffer_1000_pppd[94];

    auto g_x_0_0_0_y_z_x_zz = buffer_1000_pppd[95];

    auto g_x_0_0_0_y_z_y_xx = buffer_1000_pppd[96];

    auto g_x_0_0_0_y_z_y_xy = buffer_1000_pppd[97];

    auto g_x_0_0_0_y_z_y_xz = buffer_1000_pppd[98];

    auto g_x_0_0_0_y_z_y_yy = buffer_1000_pppd[99];

    auto g_x_0_0_0_y_z_y_yz = buffer_1000_pppd[100];

    auto g_x_0_0_0_y_z_y_zz = buffer_1000_pppd[101];

    auto g_x_0_0_0_y_z_z_xx = buffer_1000_pppd[102];

    auto g_x_0_0_0_y_z_z_xy = buffer_1000_pppd[103];

    auto g_x_0_0_0_y_z_z_xz = buffer_1000_pppd[104];

    auto g_x_0_0_0_y_z_z_yy = buffer_1000_pppd[105];

    auto g_x_0_0_0_y_z_z_yz = buffer_1000_pppd[106];

    auto g_x_0_0_0_y_z_z_zz = buffer_1000_pppd[107];

    auto g_x_0_0_0_z_x_x_xx = buffer_1000_pppd[108];

    auto g_x_0_0_0_z_x_x_xy = buffer_1000_pppd[109];

    auto g_x_0_0_0_z_x_x_xz = buffer_1000_pppd[110];

    auto g_x_0_0_0_z_x_x_yy = buffer_1000_pppd[111];

    auto g_x_0_0_0_z_x_x_yz = buffer_1000_pppd[112];

    auto g_x_0_0_0_z_x_x_zz = buffer_1000_pppd[113];

    auto g_x_0_0_0_z_x_y_xx = buffer_1000_pppd[114];

    auto g_x_0_0_0_z_x_y_xy = buffer_1000_pppd[115];

    auto g_x_0_0_0_z_x_y_xz = buffer_1000_pppd[116];

    auto g_x_0_0_0_z_x_y_yy = buffer_1000_pppd[117];

    auto g_x_0_0_0_z_x_y_yz = buffer_1000_pppd[118];

    auto g_x_0_0_0_z_x_y_zz = buffer_1000_pppd[119];

    auto g_x_0_0_0_z_x_z_xx = buffer_1000_pppd[120];

    auto g_x_0_0_0_z_x_z_xy = buffer_1000_pppd[121];

    auto g_x_0_0_0_z_x_z_xz = buffer_1000_pppd[122];

    auto g_x_0_0_0_z_x_z_yy = buffer_1000_pppd[123];

    auto g_x_0_0_0_z_x_z_yz = buffer_1000_pppd[124];

    auto g_x_0_0_0_z_x_z_zz = buffer_1000_pppd[125];

    auto g_x_0_0_0_z_y_x_xx = buffer_1000_pppd[126];

    auto g_x_0_0_0_z_y_x_xy = buffer_1000_pppd[127];

    auto g_x_0_0_0_z_y_x_xz = buffer_1000_pppd[128];

    auto g_x_0_0_0_z_y_x_yy = buffer_1000_pppd[129];

    auto g_x_0_0_0_z_y_x_yz = buffer_1000_pppd[130];

    auto g_x_0_0_0_z_y_x_zz = buffer_1000_pppd[131];

    auto g_x_0_0_0_z_y_y_xx = buffer_1000_pppd[132];

    auto g_x_0_0_0_z_y_y_xy = buffer_1000_pppd[133];

    auto g_x_0_0_0_z_y_y_xz = buffer_1000_pppd[134];

    auto g_x_0_0_0_z_y_y_yy = buffer_1000_pppd[135];

    auto g_x_0_0_0_z_y_y_yz = buffer_1000_pppd[136];

    auto g_x_0_0_0_z_y_y_zz = buffer_1000_pppd[137];

    auto g_x_0_0_0_z_y_z_xx = buffer_1000_pppd[138];

    auto g_x_0_0_0_z_y_z_xy = buffer_1000_pppd[139];

    auto g_x_0_0_0_z_y_z_xz = buffer_1000_pppd[140];

    auto g_x_0_0_0_z_y_z_yy = buffer_1000_pppd[141];

    auto g_x_0_0_0_z_y_z_yz = buffer_1000_pppd[142];

    auto g_x_0_0_0_z_y_z_zz = buffer_1000_pppd[143];

    auto g_x_0_0_0_z_z_x_xx = buffer_1000_pppd[144];

    auto g_x_0_0_0_z_z_x_xy = buffer_1000_pppd[145];

    auto g_x_0_0_0_z_z_x_xz = buffer_1000_pppd[146];

    auto g_x_0_0_0_z_z_x_yy = buffer_1000_pppd[147];

    auto g_x_0_0_0_z_z_x_yz = buffer_1000_pppd[148];

    auto g_x_0_0_0_z_z_x_zz = buffer_1000_pppd[149];

    auto g_x_0_0_0_z_z_y_xx = buffer_1000_pppd[150];

    auto g_x_0_0_0_z_z_y_xy = buffer_1000_pppd[151];

    auto g_x_0_0_0_z_z_y_xz = buffer_1000_pppd[152];

    auto g_x_0_0_0_z_z_y_yy = buffer_1000_pppd[153];

    auto g_x_0_0_0_z_z_y_yz = buffer_1000_pppd[154];

    auto g_x_0_0_0_z_z_y_zz = buffer_1000_pppd[155];

    auto g_x_0_0_0_z_z_z_xx = buffer_1000_pppd[156];

    auto g_x_0_0_0_z_z_z_xy = buffer_1000_pppd[157];

    auto g_x_0_0_0_z_z_z_xz = buffer_1000_pppd[158];

    auto g_x_0_0_0_z_z_z_yy = buffer_1000_pppd[159];

    auto g_x_0_0_0_z_z_z_yz = buffer_1000_pppd[160];

    auto g_x_0_0_0_z_z_z_zz = buffer_1000_pppd[161];

    auto g_y_0_0_0_x_x_x_xx = buffer_1000_pppd[162];

    auto g_y_0_0_0_x_x_x_xy = buffer_1000_pppd[163];

    auto g_y_0_0_0_x_x_x_xz = buffer_1000_pppd[164];

    auto g_y_0_0_0_x_x_x_yy = buffer_1000_pppd[165];

    auto g_y_0_0_0_x_x_x_yz = buffer_1000_pppd[166];

    auto g_y_0_0_0_x_x_x_zz = buffer_1000_pppd[167];

    auto g_y_0_0_0_x_x_y_xx = buffer_1000_pppd[168];

    auto g_y_0_0_0_x_x_y_xy = buffer_1000_pppd[169];

    auto g_y_0_0_0_x_x_y_xz = buffer_1000_pppd[170];

    auto g_y_0_0_0_x_x_y_yy = buffer_1000_pppd[171];

    auto g_y_0_0_0_x_x_y_yz = buffer_1000_pppd[172];

    auto g_y_0_0_0_x_x_y_zz = buffer_1000_pppd[173];

    auto g_y_0_0_0_x_x_z_xx = buffer_1000_pppd[174];

    auto g_y_0_0_0_x_x_z_xy = buffer_1000_pppd[175];

    auto g_y_0_0_0_x_x_z_xz = buffer_1000_pppd[176];

    auto g_y_0_0_0_x_x_z_yy = buffer_1000_pppd[177];

    auto g_y_0_0_0_x_x_z_yz = buffer_1000_pppd[178];

    auto g_y_0_0_0_x_x_z_zz = buffer_1000_pppd[179];

    auto g_y_0_0_0_x_y_x_xx = buffer_1000_pppd[180];

    auto g_y_0_0_0_x_y_x_xy = buffer_1000_pppd[181];

    auto g_y_0_0_0_x_y_x_xz = buffer_1000_pppd[182];

    auto g_y_0_0_0_x_y_x_yy = buffer_1000_pppd[183];

    auto g_y_0_0_0_x_y_x_yz = buffer_1000_pppd[184];

    auto g_y_0_0_0_x_y_x_zz = buffer_1000_pppd[185];

    auto g_y_0_0_0_x_y_y_xx = buffer_1000_pppd[186];

    auto g_y_0_0_0_x_y_y_xy = buffer_1000_pppd[187];

    auto g_y_0_0_0_x_y_y_xz = buffer_1000_pppd[188];

    auto g_y_0_0_0_x_y_y_yy = buffer_1000_pppd[189];

    auto g_y_0_0_0_x_y_y_yz = buffer_1000_pppd[190];

    auto g_y_0_0_0_x_y_y_zz = buffer_1000_pppd[191];

    auto g_y_0_0_0_x_y_z_xx = buffer_1000_pppd[192];

    auto g_y_0_0_0_x_y_z_xy = buffer_1000_pppd[193];

    auto g_y_0_0_0_x_y_z_xz = buffer_1000_pppd[194];

    auto g_y_0_0_0_x_y_z_yy = buffer_1000_pppd[195];

    auto g_y_0_0_0_x_y_z_yz = buffer_1000_pppd[196];

    auto g_y_0_0_0_x_y_z_zz = buffer_1000_pppd[197];

    auto g_y_0_0_0_x_z_x_xx = buffer_1000_pppd[198];

    auto g_y_0_0_0_x_z_x_xy = buffer_1000_pppd[199];

    auto g_y_0_0_0_x_z_x_xz = buffer_1000_pppd[200];

    auto g_y_0_0_0_x_z_x_yy = buffer_1000_pppd[201];

    auto g_y_0_0_0_x_z_x_yz = buffer_1000_pppd[202];

    auto g_y_0_0_0_x_z_x_zz = buffer_1000_pppd[203];

    auto g_y_0_0_0_x_z_y_xx = buffer_1000_pppd[204];

    auto g_y_0_0_0_x_z_y_xy = buffer_1000_pppd[205];

    auto g_y_0_0_0_x_z_y_xz = buffer_1000_pppd[206];

    auto g_y_0_0_0_x_z_y_yy = buffer_1000_pppd[207];

    auto g_y_0_0_0_x_z_y_yz = buffer_1000_pppd[208];

    auto g_y_0_0_0_x_z_y_zz = buffer_1000_pppd[209];

    auto g_y_0_0_0_x_z_z_xx = buffer_1000_pppd[210];

    auto g_y_0_0_0_x_z_z_xy = buffer_1000_pppd[211];

    auto g_y_0_0_0_x_z_z_xz = buffer_1000_pppd[212];

    auto g_y_0_0_0_x_z_z_yy = buffer_1000_pppd[213];

    auto g_y_0_0_0_x_z_z_yz = buffer_1000_pppd[214];

    auto g_y_0_0_0_x_z_z_zz = buffer_1000_pppd[215];

    auto g_y_0_0_0_y_x_x_xx = buffer_1000_pppd[216];

    auto g_y_0_0_0_y_x_x_xy = buffer_1000_pppd[217];

    auto g_y_0_0_0_y_x_x_xz = buffer_1000_pppd[218];

    auto g_y_0_0_0_y_x_x_yy = buffer_1000_pppd[219];

    auto g_y_0_0_0_y_x_x_yz = buffer_1000_pppd[220];

    auto g_y_0_0_0_y_x_x_zz = buffer_1000_pppd[221];

    auto g_y_0_0_0_y_x_y_xx = buffer_1000_pppd[222];

    auto g_y_0_0_0_y_x_y_xy = buffer_1000_pppd[223];

    auto g_y_0_0_0_y_x_y_xz = buffer_1000_pppd[224];

    auto g_y_0_0_0_y_x_y_yy = buffer_1000_pppd[225];

    auto g_y_0_0_0_y_x_y_yz = buffer_1000_pppd[226];

    auto g_y_0_0_0_y_x_y_zz = buffer_1000_pppd[227];

    auto g_y_0_0_0_y_x_z_xx = buffer_1000_pppd[228];

    auto g_y_0_0_0_y_x_z_xy = buffer_1000_pppd[229];

    auto g_y_0_0_0_y_x_z_xz = buffer_1000_pppd[230];

    auto g_y_0_0_0_y_x_z_yy = buffer_1000_pppd[231];

    auto g_y_0_0_0_y_x_z_yz = buffer_1000_pppd[232];

    auto g_y_0_0_0_y_x_z_zz = buffer_1000_pppd[233];

    auto g_y_0_0_0_y_y_x_xx = buffer_1000_pppd[234];

    auto g_y_0_0_0_y_y_x_xy = buffer_1000_pppd[235];

    auto g_y_0_0_0_y_y_x_xz = buffer_1000_pppd[236];

    auto g_y_0_0_0_y_y_x_yy = buffer_1000_pppd[237];

    auto g_y_0_0_0_y_y_x_yz = buffer_1000_pppd[238];

    auto g_y_0_0_0_y_y_x_zz = buffer_1000_pppd[239];

    auto g_y_0_0_0_y_y_y_xx = buffer_1000_pppd[240];

    auto g_y_0_0_0_y_y_y_xy = buffer_1000_pppd[241];

    auto g_y_0_0_0_y_y_y_xz = buffer_1000_pppd[242];

    auto g_y_0_0_0_y_y_y_yy = buffer_1000_pppd[243];

    auto g_y_0_0_0_y_y_y_yz = buffer_1000_pppd[244];

    auto g_y_0_0_0_y_y_y_zz = buffer_1000_pppd[245];

    auto g_y_0_0_0_y_y_z_xx = buffer_1000_pppd[246];

    auto g_y_0_0_0_y_y_z_xy = buffer_1000_pppd[247];

    auto g_y_0_0_0_y_y_z_xz = buffer_1000_pppd[248];

    auto g_y_0_0_0_y_y_z_yy = buffer_1000_pppd[249];

    auto g_y_0_0_0_y_y_z_yz = buffer_1000_pppd[250];

    auto g_y_0_0_0_y_y_z_zz = buffer_1000_pppd[251];

    auto g_y_0_0_0_y_z_x_xx = buffer_1000_pppd[252];

    auto g_y_0_0_0_y_z_x_xy = buffer_1000_pppd[253];

    auto g_y_0_0_0_y_z_x_xz = buffer_1000_pppd[254];

    auto g_y_0_0_0_y_z_x_yy = buffer_1000_pppd[255];

    auto g_y_0_0_0_y_z_x_yz = buffer_1000_pppd[256];

    auto g_y_0_0_0_y_z_x_zz = buffer_1000_pppd[257];

    auto g_y_0_0_0_y_z_y_xx = buffer_1000_pppd[258];

    auto g_y_0_0_0_y_z_y_xy = buffer_1000_pppd[259];

    auto g_y_0_0_0_y_z_y_xz = buffer_1000_pppd[260];

    auto g_y_0_0_0_y_z_y_yy = buffer_1000_pppd[261];

    auto g_y_0_0_0_y_z_y_yz = buffer_1000_pppd[262];

    auto g_y_0_0_0_y_z_y_zz = buffer_1000_pppd[263];

    auto g_y_0_0_0_y_z_z_xx = buffer_1000_pppd[264];

    auto g_y_0_0_0_y_z_z_xy = buffer_1000_pppd[265];

    auto g_y_0_0_0_y_z_z_xz = buffer_1000_pppd[266];

    auto g_y_0_0_0_y_z_z_yy = buffer_1000_pppd[267];

    auto g_y_0_0_0_y_z_z_yz = buffer_1000_pppd[268];

    auto g_y_0_0_0_y_z_z_zz = buffer_1000_pppd[269];

    auto g_y_0_0_0_z_x_x_xx = buffer_1000_pppd[270];

    auto g_y_0_0_0_z_x_x_xy = buffer_1000_pppd[271];

    auto g_y_0_0_0_z_x_x_xz = buffer_1000_pppd[272];

    auto g_y_0_0_0_z_x_x_yy = buffer_1000_pppd[273];

    auto g_y_0_0_0_z_x_x_yz = buffer_1000_pppd[274];

    auto g_y_0_0_0_z_x_x_zz = buffer_1000_pppd[275];

    auto g_y_0_0_0_z_x_y_xx = buffer_1000_pppd[276];

    auto g_y_0_0_0_z_x_y_xy = buffer_1000_pppd[277];

    auto g_y_0_0_0_z_x_y_xz = buffer_1000_pppd[278];

    auto g_y_0_0_0_z_x_y_yy = buffer_1000_pppd[279];

    auto g_y_0_0_0_z_x_y_yz = buffer_1000_pppd[280];

    auto g_y_0_0_0_z_x_y_zz = buffer_1000_pppd[281];

    auto g_y_0_0_0_z_x_z_xx = buffer_1000_pppd[282];

    auto g_y_0_0_0_z_x_z_xy = buffer_1000_pppd[283];

    auto g_y_0_0_0_z_x_z_xz = buffer_1000_pppd[284];

    auto g_y_0_0_0_z_x_z_yy = buffer_1000_pppd[285];

    auto g_y_0_0_0_z_x_z_yz = buffer_1000_pppd[286];

    auto g_y_0_0_0_z_x_z_zz = buffer_1000_pppd[287];

    auto g_y_0_0_0_z_y_x_xx = buffer_1000_pppd[288];

    auto g_y_0_0_0_z_y_x_xy = buffer_1000_pppd[289];

    auto g_y_0_0_0_z_y_x_xz = buffer_1000_pppd[290];

    auto g_y_0_0_0_z_y_x_yy = buffer_1000_pppd[291];

    auto g_y_0_0_0_z_y_x_yz = buffer_1000_pppd[292];

    auto g_y_0_0_0_z_y_x_zz = buffer_1000_pppd[293];

    auto g_y_0_0_0_z_y_y_xx = buffer_1000_pppd[294];

    auto g_y_0_0_0_z_y_y_xy = buffer_1000_pppd[295];

    auto g_y_0_0_0_z_y_y_xz = buffer_1000_pppd[296];

    auto g_y_0_0_0_z_y_y_yy = buffer_1000_pppd[297];

    auto g_y_0_0_0_z_y_y_yz = buffer_1000_pppd[298];

    auto g_y_0_0_0_z_y_y_zz = buffer_1000_pppd[299];

    auto g_y_0_0_0_z_y_z_xx = buffer_1000_pppd[300];

    auto g_y_0_0_0_z_y_z_xy = buffer_1000_pppd[301];

    auto g_y_0_0_0_z_y_z_xz = buffer_1000_pppd[302];

    auto g_y_0_0_0_z_y_z_yy = buffer_1000_pppd[303];

    auto g_y_0_0_0_z_y_z_yz = buffer_1000_pppd[304];

    auto g_y_0_0_0_z_y_z_zz = buffer_1000_pppd[305];

    auto g_y_0_0_0_z_z_x_xx = buffer_1000_pppd[306];

    auto g_y_0_0_0_z_z_x_xy = buffer_1000_pppd[307];

    auto g_y_0_0_0_z_z_x_xz = buffer_1000_pppd[308];

    auto g_y_0_0_0_z_z_x_yy = buffer_1000_pppd[309];

    auto g_y_0_0_0_z_z_x_yz = buffer_1000_pppd[310];

    auto g_y_0_0_0_z_z_x_zz = buffer_1000_pppd[311];

    auto g_y_0_0_0_z_z_y_xx = buffer_1000_pppd[312];

    auto g_y_0_0_0_z_z_y_xy = buffer_1000_pppd[313];

    auto g_y_0_0_0_z_z_y_xz = buffer_1000_pppd[314];

    auto g_y_0_0_0_z_z_y_yy = buffer_1000_pppd[315];

    auto g_y_0_0_0_z_z_y_yz = buffer_1000_pppd[316];

    auto g_y_0_0_0_z_z_y_zz = buffer_1000_pppd[317];

    auto g_y_0_0_0_z_z_z_xx = buffer_1000_pppd[318];

    auto g_y_0_0_0_z_z_z_xy = buffer_1000_pppd[319];

    auto g_y_0_0_0_z_z_z_xz = buffer_1000_pppd[320];

    auto g_y_0_0_0_z_z_z_yy = buffer_1000_pppd[321];

    auto g_y_0_0_0_z_z_z_yz = buffer_1000_pppd[322];

    auto g_y_0_0_0_z_z_z_zz = buffer_1000_pppd[323];

    auto g_z_0_0_0_x_x_x_xx = buffer_1000_pppd[324];

    auto g_z_0_0_0_x_x_x_xy = buffer_1000_pppd[325];

    auto g_z_0_0_0_x_x_x_xz = buffer_1000_pppd[326];

    auto g_z_0_0_0_x_x_x_yy = buffer_1000_pppd[327];

    auto g_z_0_0_0_x_x_x_yz = buffer_1000_pppd[328];

    auto g_z_0_0_0_x_x_x_zz = buffer_1000_pppd[329];

    auto g_z_0_0_0_x_x_y_xx = buffer_1000_pppd[330];

    auto g_z_0_0_0_x_x_y_xy = buffer_1000_pppd[331];

    auto g_z_0_0_0_x_x_y_xz = buffer_1000_pppd[332];

    auto g_z_0_0_0_x_x_y_yy = buffer_1000_pppd[333];

    auto g_z_0_0_0_x_x_y_yz = buffer_1000_pppd[334];

    auto g_z_0_0_0_x_x_y_zz = buffer_1000_pppd[335];

    auto g_z_0_0_0_x_x_z_xx = buffer_1000_pppd[336];

    auto g_z_0_0_0_x_x_z_xy = buffer_1000_pppd[337];

    auto g_z_0_0_0_x_x_z_xz = buffer_1000_pppd[338];

    auto g_z_0_0_0_x_x_z_yy = buffer_1000_pppd[339];

    auto g_z_0_0_0_x_x_z_yz = buffer_1000_pppd[340];

    auto g_z_0_0_0_x_x_z_zz = buffer_1000_pppd[341];

    auto g_z_0_0_0_x_y_x_xx = buffer_1000_pppd[342];

    auto g_z_0_0_0_x_y_x_xy = buffer_1000_pppd[343];

    auto g_z_0_0_0_x_y_x_xz = buffer_1000_pppd[344];

    auto g_z_0_0_0_x_y_x_yy = buffer_1000_pppd[345];

    auto g_z_0_0_0_x_y_x_yz = buffer_1000_pppd[346];

    auto g_z_0_0_0_x_y_x_zz = buffer_1000_pppd[347];

    auto g_z_0_0_0_x_y_y_xx = buffer_1000_pppd[348];

    auto g_z_0_0_0_x_y_y_xy = buffer_1000_pppd[349];

    auto g_z_0_0_0_x_y_y_xz = buffer_1000_pppd[350];

    auto g_z_0_0_0_x_y_y_yy = buffer_1000_pppd[351];

    auto g_z_0_0_0_x_y_y_yz = buffer_1000_pppd[352];

    auto g_z_0_0_0_x_y_y_zz = buffer_1000_pppd[353];

    auto g_z_0_0_0_x_y_z_xx = buffer_1000_pppd[354];

    auto g_z_0_0_0_x_y_z_xy = buffer_1000_pppd[355];

    auto g_z_0_0_0_x_y_z_xz = buffer_1000_pppd[356];

    auto g_z_0_0_0_x_y_z_yy = buffer_1000_pppd[357];

    auto g_z_0_0_0_x_y_z_yz = buffer_1000_pppd[358];

    auto g_z_0_0_0_x_y_z_zz = buffer_1000_pppd[359];

    auto g_z_0_0_0_x_z_x_xx = buffer_1000_pppd[360];

    auto g_z_0_0_0_x_z_x_xy = buffer_1000_pppd[361];

    auto g_z_0_0_0_x_z_x_xz = buffer_1000_pppd[362];

    auto g_z_0_0_0_x_z_x_yy = buffer_1000_pppd[363];

    auto g_z_0_0_0_x_z_x_yz = buffer_1000_pppd[364];

    auto g_z_0_0_0_x_z_x_zz = buffer_1000_pppd[365];

    auto g_z_0_0_0_x_z_y_xx = buffer_1000_pppd[366];

    auto g_z_0_0_0_x_z_y_xy = buffer_1000_pppd[367];

    auto g_z_0_0_0_x_z_y_xz = buffer_1000_pppd[368];

    auto g_z_0_0_0_x_z_y_yy = buffer_1000_pppd[369];

    auto g_z_0_0_0_x_z_y_yz = buffer_1000_pppd[370];

    auto g_z_0_0_0_x_z_y_zz = buffer_1000_pppd[371];

    auto g_z_0_0_0_x_z_z_xx = buffer_1000_pppd[372];

    auto g_z_0_0_0_x_z_z_xy = buffer_1000_pppd[373];

    auto g_z_0_0_0_x_z_z_xz = buffer_1000_pppd[374];

    auto g_z_0_0_0_x_z_z_yy = buffer_1000_pppd[375];

    auto g_z_0_0_0_x_z_z_yz = buffer_1000_pppd[376];

    auto g_z_0_0_0_x_z_z_zz = buffer_1000_pppd[377];

    auto g_z_0_0_0_y_x_x_xx = buffer_1000_pppd[378];

    auto g_z_0_0_0_y_x_x_xy = buffer_1000_pppd[379];

    auto g_z_0_0_0_y_x_x_xz = buffer_1000_pppd[380];

    auto g_z_0_0_0_y_x_x_yy = buffer_1000_pppd[381];

    auto g_z_0_0_0_y_x_x_yz = buffer_1000_pppd[382];

    auto g_z_0_0_0_y_x_x_zz = buffer_1000_pppd[383];

    auto g_z_0_0_0_y_x_y_xx = buffer_1000_pppd[384];

    auto g_z_0_0_0_y_x_y_xy = buffer_1000_pppd[385];

    auto g_z_0_0_0_y_x_y_xz = buffer_1000_pppd[386];

    auto g_z_0_0_0_y_x_y_yy = buffer_1000_pppd[387];

    auto g_z_0_0_0_y_x_y_yz = buffer_1000_pppd[388];

    auto g_z_0_0_0_y_x_y_zz = buffer_1000_pppd[389];

    auto g_z_0_0_0_y_x_z_xx = buffer_1000_pppd[390];

    auto g_z_0_0_0_y_x_z_xy = buffer_1000_pppd[391];

    auto g_z_0_0_0_y_x_z_xz = buffer_1000_pppd[392];

    auto g_z_0_0_0_y_x_z_yy = buffer_1000_pppd[393];

    auto g_z_0_0_0_y_x_z_yz = buffer_1000_pppd[394];

    auto g_z_0_0_0_y_x_z_zz = buffer_1000_pppd[395];

    auto g_z_0_0_0_y_y_x_xx = buffer_1000_pppd[396];

    auto g_z_0_0_0_y_y_x_xy = buffer_1000_pppd[397];

    auto g_z_0_0_0_y_y_x_xz = buffer_1000_pppd[398];

    auto g_z_0_0_0_y_y_x_yy = buffer_1000_pppd[399];

    auto g_z_0_0_0_y_y_x_yz = buffer_1000_pppd[400];

    auto g_z_0_0_0_y_y_x_zz = buffer_1000_pppd[401];

    auto g_z_0_0_0_y_y_y_xx = buffer_1000_pppd[402];

    auto g_z_0_0_0_y_y_y_xy = buffer_1000_pppd[403];

    auto g_z_0_0_0_y_y_y_xz = buffer_1000_pppd[404];

    auto g_z_0_0_0_y_y_y_yy = buffer_1000_pppd[405];

    auto g_z_0_0_0_y_y_y_yz = buffer_1000_pppd[406];

    auto g_z_0_0_0_y_y_y_zz = buffer_1000_pppd[407];

    auto g_z_0_0_0_y_y_z_xx = buffer_1000_pppd[408];

    auto g_z_0_0_0_y_y_z_xy = buffer_1000_pppd[409];

    auto g_z_0_0_0_y_y_z_xz = buffer_1000_pppd[410];

    auto g_z_0_0_0_y_y_z_yy = buffer_1000_pppd[411];

    auto g_z_0_0_0_y_y_z_yz = buffer_1000_pppd[412];

    auto g_z_0_0_0_y_y_z_zz = buffer_1000_pppd[413];

    auto g_z_0_0_0_y_z_x_xx = buffer_1000_pppd[414];

    auto g_z_0_0_0_y_z_x_xy = buffer_1000_pppd[415];

    auto g_z_0_0_0_y_z_x_xz = buffer_1000_pppd[416];

    auto g_z_0_0_0_y_z_x_yy = buffer_1000_pppd[417];

    auto g_z_0_0_0_y_z_x_yz = buffer_1000_pppd[418];

    auto g_z_0_0_0_y_z_x_zz = buffer_1000_pppd[419];

    auto g_z_0_0_0_y_z_y_xx = buffer_1000_pppd[420];

    auto g_z_0_0_0_y_z_y_xy = buffer_1000_pppd[421];

    auto g_z_0_0_0_y_z_y_xz = buffer_1000_pppd[422];

    auto g_z_0_0_0_y_z_y_yy = buffer_1000_pppd[423];

    auto g_z_0_0_0_y_z_y_yz = buffer_1000_pppd[424];

    auto g_z_0_0_0_y_z_y_zz = buffer_1000_pppd[425];

    auto g_z_0_0_0_y_z_z_xx = buffer_1000_pppd[426];

    auto g_z_0_0_0_y_z_z_xy = buffer_1000_pppd[427];

    auto g_z_0_0_0_y_z_z_xz = buffer_1000_pppd[428];

    auto g_z_0_0_0_y_z_z_yy = buffer_1000_pppd[429];

    auto g_z_0_0_0_y_z_z_yz = buffer_1000_pppd[430];

    auto g_z_0_0_0_y_z_z_zz = buffer_1000_pppd[431];

    auto g_z_0_0_0_z_x_x_xx = buffer_1000_pppd[432];

    auto g_z_0_0_0_z_x_x_xy = buffer_1000_pppd[433];

    auto g_z_0_0_0_z_x_x_xz = buffer_1000_pppd[434];

    auto g_z_0_0_0_z_x_x_yy = buffer_1000_pppd[435];

    auto g_z_0_0_0_z_x_x_yz = buffer_1000_pppd[436];

    auto g_z_0_0_0_z_x_x_zz = buffer_1000_pppd[437];

    auto g_z_0_0_0_z_x_y_xx = buffer_1000_pppd[438];

    auto g_z_0_0_0_z_x_y_xy = buffer_1000_pppd[439];

    auto g_z_0_0_0_z_x_y_xz = buffer_1000_pppd[440];

    auto g_z_0_0_0_z_x_y_yy = buffer_1000_pppd[441];

    auto g_z_0_0_0_z_x_y_yz = buffer_1000_pppd[442];

    auto g_z_0_0_0_z_x_y_zz = buffer_1000_pppd[443];

    auto g_z_0_0_0_z_x_z_xx = buffer_1000_pppd[444];

    auto g_z_0_0_0_z_x_z_xy = buffer_1000_pppd[445];

    auto g_z_0_0_0_z_x_z_xz = buffer_1000_pppd[446];

    auto g_z_0_0_0_z_x_z_yy = buffer_1000_pppd[447];

    auto g_z_0_0_0_z_x_z_yz = buffer_1000_pppd[448];

    auto g_z_0_0_0_z_x_z_zz = buffer_1000_pppd[449];

    auto g_z_0_0_0_z_y_x_xx = buffer_1000_pppd[450];

    auto g_z_0_0_0_z_y_x_xy = buffer_1000_pppd[451];

    auto g_z_0_0_0_z_y_x_xz = buffer_1000_pppd[452];

    auto g_z_0_0_0_z_y_x_yy = buffer_1000_pppd[453];

    auto g_z_0_0_0_z_y_x_yz = buffer_1000_pppd[454];

    auto g_z_0_0_0_z_y_x_zz = buffer_1000_pppd[455];

    auto g_z_0_0_0_z_y_y_xx = buffer_1000_pppd[456];

    auto g_z_0_0_0_z_y_y_xy = buffer_1000_pppd[457];

    auto g_z_0_0_0_z_y_y_xz = buffer_1000_pppd[458];

    auto g_z_0_0_0_z_y_y_yy = buffer_1000_pppd[459];

    auto g_z_0_0_0_z_y_y_yz = buffer_1000_pppd[460];

    auto g_z_0_0_0_z_y_y_zz = buffer_1000_pppd[461];

    auto g_z_0_0_0_z_y_z_xx = buffer_1000_pppd[462];

    auto g_z_0_0_0_z_y_z_xy = buffer_1000_pppd[463];

    auto g_z_0_0_0_z_y_z_xz = buffer_1000_pppd[464];

    auto g_z_0_0_0_z_y_z_yy = buffer_1000_pppd[465];

    auto g_z_0_0_0_z_y_z_yz = buffer_1000_pppd[466];

    auto g_z_0_0_0_z_y_z_zz = buffer_1000_pppd[467];

    auto g_z_0_0_0_z_z_x_xx = buffer_1000_pppd[468];

    auto g_z_0_0_0_z_z_x_xy = buffer_1000_pppd[469];

    auto g_z_0_0_0_z_z_x_xz = buffer_1000_pppd[470];

    auto g_z_0_0_0_z_z_x_yy = buffer_1000_pppd[471];

    auto g_z_0_0_0_z_z_x_yz = buffer_1000_pppd[472];

    auto g_z_0_0_0_z_z_x_zz = buffer_1000_pppd[473];

    auto g_z_0_0_0_z_z_y_xx = buffer_1000_pppd[474];

    auto g_z_0_0_0_z_z_y_xy = buffer_1000_pppd[475];

    auto g_z_0_0_0_z_z_y_xz = buffer_1000_pppd[476];

    auto g_z_0_0_0_z_z_y_yy = buffer_1000_pppd[477];

    auto g_z_0_0_0_z_z_y_yz = buffer_1000_pppd[478];

    auto g_z_0_0_0_z_z_y_zz = buffer_1000_pppd[479];

    auto g_z_0_0_0_z_z_z_xx = buffer_1000_pppd[480];

    auto g_z_0_0_0_z_z_z_xy = buffer_1000_pppd[481];

    auto g_z_0_0_0_z_z_z_xz = buffer_1000_pppd[482];

    auto g_z_0_0_0_z_z_z_yy = buffer_1000_pppd[483];

    auto g_z_0_0_0_z_z_z_yz = buffer_1000_pppd[484];

    auto g_z_0_0_0_z_z_z_zz = buffer_1000_pppd[485];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_x_0_0_0_x_x_x_xx, g_x_0_0_0_x_x_x_xy, g_x_0_0_0_x_x_x_xz, g_x_0_0_0_x_x_x_yy, g_x_0_0_0_x_x_x_yz, g_x_0_0_0_x_x_x_zz, g_xx_x_x_xx, g_xx_x_x_xy, g_xx_x_x_xz, g_xx_x_x_yy, g_xx_x_x_yz, g_xx_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_x_x_xx[i] = -g_0_x_x_xx[i] + 2.0 * g_xx_x_x_xx[i] * a_exp;

        g_x_0_0_0_x_x_x_xy[i] = -g_0_x_x_xy[i] + 2.0 * g_xx_x_x_xy[i] * a_exp;

        g_x_0_0_0_x_x_x_xz[i] = -g_0_x_x_xz[i] + 2.0 * g_xx_x_x_xz[i] * a_exp;

        g_x_0_0_0_x_x_x_yy[i] = -g_0_x_x_yy[i] + 2.0 * g_xx_x_x_yy[i] * a_exp;

        g_x_0_0_0_x_x_x_yz[i] = -g_0_x_x_yz[i] + 2.0 * g_xx_x_x_yz[i] * a_exp;

        g_x_0_0_0_x_x_x_zz[i] = -g_0_x_x_zz[i] + 2.0 * g_xx_x_x_zz[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_x_0_0_0_x_x_y_xx, g_x_0_0_0_x_x_y_xy, g_x_0_0_0_x_x_y_xz, g_x_0_0_0_x_x_y_yy, g_x_0_0_0_x_x_y_yz, g_x_0_0_0_x_x_y_zz, g_xx_x_y_xx, g_xx_x_y_xy, g_xx_x_y_xz, g_xx_x_y_yy, g_xx_x_y_yz, g_xx_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_x_y_xx[i] = -g_0_x_y_xx[i] + 2.0 * g_xx_x_y_xx[i] * a_exp;

        g_x_0_0_0_x_x_y_xy[i] = -g_0_x_y_xy[i] + 2.0 * g_xx_x_y_xy[i] * a_exp;

        g_x_0_0_0_x_x_y_xz[i] = -g_0_x_y_xz[i] + 2.0 * g_xx_x_y_xz[i] * a_exp;

        g_x_0_0_0_x_x_y_yy[i] = -g_0_x_y_yy[i] + 2.0 * g_xx_x_y_yy[i] * a_exp;

        g_x_0_0_0_x_x_y_yz[i] = -g_0_x_y_yz[i] + 2.0 * g_xx_x_y_yz[i] * a_exp;

        g_x_0_0_0_x_x_y_zz[i] = -g_0_x_y_zz[i] + 2.0 * g_xx_x_y_zz[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_x_0_0_0_x_x_z_xx, g_x_0_0_0_x_x_z_xy, g_x_0_0_0_x_x_z_xz, g_x_0_0_0_x_x_z_yy, g_x_0_0_0_x_x_z_yz, g_x_0_0_0_x_x_z_zz, g_xx_x_z_xx, g_xx_x_z_xy, g_xx_x_z_xz, g_xx_x_z_yy, g_xx_x_z_yz, g_xx_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_x_z_xx[i] = -g_0_x_z_xx[i] + 2.0 * g_xx_x_z_xx[i] * a_exp;

        g_x_0_0_0_x_x_z_xy[i] = -g_0_x_z_xy[i] + 2.0 * g_xx_x_z_xy[i] * a_exp;

        g_x_0_0_0_x_x_z_xz[i] = -g_0_x_z_xz[i] + 2.0 * g_xx_x_z_xz[i] * a_exp;

        g_x_0_0_0_x_x_z_yy[i] = -g_0_x_z_yy[i] + 2.0 * g_xx_x_z_yy[i] * a_exp;

        g_x_0_0_0_x_x_z_yz[i] = -g_0_x_z_yz[i] + 2.0 * g_xx_x_z_yz[i] * a_exp;

        g_x_0_0_0_x_x_z_zz[i] = -g_0_x_z_zz[i] + 2.0 * g_xx_x_z_zz[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_x_0_0_0_x_y_x_xx, g_x_0_0_0_x_y_x_xy, g_x_0_0_0_x_y_x_xz, g_x_0_0_0_x_y_x_yy, g_x_0_0_0_x_y_x_yz, g_x_0_0_0_x_y_x_zz, g_xx_y_x_xx, g_xx_y_x_xy, g_xx_y_x_xz, g_xx_y_x_yy, g_xx_y_x_yz, g_xx_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_y_x_xx[i] = -g_0_y_x_xx[i] + 2.0 * g_xx_y_x_xx[i] * a_exp;

        g_x_0_0_0_x_y_x_xy[i] = -g_0_y_x_xy[i] + 2.0 * g_xx_y_x_xy[i] * a_exp;

        g_x_0_0_0_x_y_x_xz[i] = -g_0_y_x_xz[i] + 2.0 * g_xx_y_x_xz[i] * a_exp;

        g_x_0_0_0_x_y_x_yy[i] = -g_0_y_x_yy[i] + 2.0 * g_xx_y_x_yy[i] * a_exp;

        g_x_0_0_0_x_y_x_yz[i] = -g_0_y_x_yz[i] + 2.0 * g_xx_y_x_yz[i] * a_exp;

        g_x_0_0_0_x_y_x_zz[i] = -g_0_y_x_zz[i] + 2.0 * g_xx_y_x_zz[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_x_0_0_0_x_y_y_xx, g_x_0_0_0_x_y_y_xy, g_x_0_0_0_x_y_y_xz, g_x_0_0_0_x_y_y_yy, g_x_0_0_0_x_y_y_yz, g_x_0_0_0_x_y_y_zz, g_xx_y_y_xx, g_xx_y_y_xy, g_xx_y_y_xz, g_xx_y_y_yy, g_xx_y_y_yz, g_xx_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_y_y_xx[i] = -g_0_y_y_xx[i] + 2.0 * g_xx_y_y_xx[i] * a_exp;

        g_x_0_0_0_x_y_y_xy[i] = -g_0_y_y_xy[i] + 2.0 * g_xx_y_y_xy[i] * a_exp;

        g_x_0_0_0_x_y_y_xz[i] = -g_0_y_y_xz[i] + 2.0 * g_xx_y_y_xz[i] * a_exp;

        g_x_0_0_0_x_y_y_yy[i] = -g_0_y_y_yy[i] + 2.0 * g_xx_y_y_yy[i] * a_exp;

        g_x_0_0_0_x_y_y_yz[i] = -g_0_y_y_yz[i] + 2.0 * g_xx_y_y_yz[i] * a_exp;

        g_x_0_0_0_x_y_y_zz[i] = -g_0_y_y_zz[i] + 2.0 * g_xx_y_y_zz[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_x_0_0_0_x_y_z_xx, g_x_0_0_0_x_y_z_xy, g_x_0_0_0_x_y_z_xz, g_x_0_0_0_x_y_z_yy, g_x_0_0_0_x_y_z_yz, g_x_0_0_0_x_y_z_zz, g_xx_y_z_xx, g_xx_y_z_xy, g_xx_y_z_xz, g_xx_y_z_yy, g_xx_y_z_yz, g_xx_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_y_z_xx[i] = -g_0_y_z_xx[i] + 2.0 * g_xx_y_z_xx[i] * a_exp;

        g_x_0_0_0_x_y_z_xy[i] = -g_0_y_z_xy[i] + 2.0 * g_xx_y_z_xy[i] * a_exp;

        g_x_0_0_0_x_y_z_xz[i] = -g_0_y_z_xz[i] + 2.0 * g_xx_y_z_xz[i] * a_exp;

        g_x_0_0_0_x_y_z_yy[i] = -g_0_y_z_yy[i] + 2.0 * g_xx_y_z_yy[i] * a_exp;

        g_x_0_0_0_x_y_z_yz[i] = -g_0_y_z_yz[i] + 2.0 * g_xx_y_z_yz[i] * a_exp;

        g_x_0_0_0_x_y_z_zz[i] = -g_0_y_z_zz[i] + 2.0 * g_xx_y_z_zz[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_x_0_0_0_x_z_x_xx, g_x_0_0_0_x_z_x_xy, g_x_0_0_0_x_z_x_xz, g_x_0_0_0_x_z_x_yy, g_x_0_0_0_x_z_x_yz, g_x_0_0_0_x_z_x_zz, g_xx_z_x_xx, g_xx_z_x_xy, g_xx_z_x_xz, g_xx_z_x_yy, g_xx_z_x_yz, g_xx_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_z_x_xx[i] = -g_0_z_x_xx[i] + 2.0 * g_xx_z_x_xx[i] * a_exp;

        g_x_0_0_0_x_z_x_xy[i] = -g_0_z_x_xy[i] + 2.0 * g_xx_z_x_xy[i] * a_exp;

        g_x_0_0_0_x_z_x_xz[i] = -g_0_z_x_xz[i] + 2.0 * g_xx_z_x_xz[i] * a_exp;

        g_x_0_0_0_x_z_x_yy[i] = -g_0_z_x_yy[i] + 2.0 * g_xx_z_x_yy[i] * a_exp;

        g_x_0_0_0_x_z_x_yz[i] = -g_0_z_x_yz[i] + 2.0 * g_xx_z_x_yz[i] * a_exp;

        g_x_0_0_0_x_z_x_zz[i] = -g_0_z_x_zz[i] + 2.0 * g_xx_z_x_zz[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_x_0_0_0_x_z_y_xx, g_x_0_0_0_x_z_y_xy, g_x_0_0_0_x_z_y_xz, g_x_0_0_0_x_z_y_yy, g_x_0_0_0_x_z_y_yz, g_x_0_0_0_x_z_y_zz, g_xx_z_y_xx, g_xx_z_y_xy, g_xx_z_y_xz, g_xx_z_y_yy, g_xx_z_y_yz, g_xx_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_z_y_xx[i] = -g_0_z_y_xx[i] + 2.0 * g_xx_z_y_xx[i] * a_exp;

        g_x_0_0_0_x_z_y_xy[i] = -g_0_z_y_xy[i] + 2.0 * g_xx_z_y_xy[i] * a_exp;

        g_x_0_0_0_x_z_y_xz[i] = -g_0_z_y_xz[i] + 2.0 * g_xx_z_y_xz[i] * a_exp;

        g_x_0_0_0_x_z_y_yy[i] = -g_0_z_y_yy[i] + 2.0 * g_xx_z_y_yy[i] * a_exp;

        g_x_0_0_0_x_z_y_yz[i] = -g_0_z_y_yz[i] + 2.0 * g_xx_z_y_yz[i] * a_exp;

        g_x_0_0_0_x_z_y_zz[i] = -g_0_z_y_zz[i] + 2.0 * g_xx_z_y_zz[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_x_0_0_0_x_z_z_xx, g_x_0_0_0_x_z_z_xy, g_x_0_0_0_x_z_z_xz, g_x_0_0_0_x_z_z_yy, g_x_0_0_0_x_z_z_yz, g_x_0_0_0_x_z_z_zz, g_xx_z_z_xx, g_xx_z_z_xy, g_xx_z_z_xz, g_xx_z_z_yy, g_xx_z_z_yz, g_xx_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_z_z_xx[i] = -g_0_z_z_xx[i] + 2.0 * g_xx_z_z_xx[i] * a_exp;

        g_x_0_0_0_x_z_z_xy[i] = -g_0_z_z_xy[i] + 2.0 * g_xx_z_z_xy[i] * a_exp;

        g_x_0_0_0_x_z_z_xz[i] = -g_0_z_z_xz[i] + 2.0 * g_xx_z_z_xz[i] * a_exp;

        g_x_0_0_0_x_z_z_yy[i] = -g_0_z_z_yy[i] + 2.0 * g_xx_z_z_yy[i] * a_exp;

        g_x_0_0_0_x_z_z_yz[i] = -g_0_z_z_yz[i] + 2.0 * g_xx_z_z_yz[i] * a_exp;

        g_x_0_0_0_x_z_z_zz[i] = -g_0_z_z_zz[i] + 2.0 * g_xx_z_z_zz[i] * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_0_0_0_y_x_x_xx, g_x_0_0_0_y_x_x_xy, g_x_0_0_0_y_x_x_xz, g_x_0_0_0_y_x_x_yy, g_x_0_0_0_y_x_x_yz, g_x_0_0_0_y_x_x_zz, g_xy_x_x_xx, g_xy_x_x_xy, g_xy_x_x_xz, g_xy_x_x_yy, g_xy_x_x_yz, g_xy_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_x_x_xx[i] = 2.0 * g_xy_x_x_xx[i] * a_exp;

        g_x_0_0_0_y_x_x_xy[i] = 2.0 * g_xy_x_x_xy[i] * a_exp;

        g_x_0_0_0_y_x_x_xz[i] = 2.0 * g_xy_x_x_xz[i] * a_exp;

        g_x_0_0_0_y_x_x_yy[i] = 2.0 * g_xy_x_x_yy[i] * a_exp;

        g_x_0_0_0_y_x_x_yz[i] = 2.0 * g_xy_x_x_yz[i] * a_exp;

        g_x_0_0_0_y_x_x_zz[i] = 2.0 * g_xy_x_x_zz[i] * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_0_0_0_y_x_y_xx, g_x_0_0_0_y_x_y_xy, g_x_0_0_0_y_x_y_xz, g_x_0_0_0_y_x_y_yy, g_x_0_0_0_y_x_y_yz, g_x_0_0_0_y_x_y_zz, g_xy_x_y_xx, g_xy_x_y_xy, g_xy_x_y_xz, g_xy_x_y_yy, g_xy_x_y_yz, g_xy_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_x_y_xx[i] = 2.0 * g_xy_x_y_xx[i] * a_exp;

        g_x_0_0_0_y_x_y_xy[i] = 2.0 * g_xy_x_y_xy[i] * a_exp;

        g_x_0_0_0_y_x_y_xz[i] = 2.0 * g_xy_x_y_xz[i] * a_exp;

        g_x_0_0_0_y_x_y_yy[i] = 2.0 * g_xy_x_y_yy[i] * a_exp;

        g_x_0_0_0_y_x_y_yz[i] = 2.0 * g_xy_x_y_yz[i] * a_exp;

        g_x_0_0_0_y_x_y_zz[i] = 2.0 * g_xy_x_y_zz[i] * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_0_0_0_y_x_z_xx, g_x_0_0_0_y_x_z_xy, g_x_0_0_0_y_x_z_xz, g_x_0_0_0_y_x_z_yy, g_x_0_0_0_y_x_z_yz, g_x_0_0_0_y_x_z_zz, g_xy_x_z_xx, g_xy_x_z_xy, g_xy_x_z_xz, g_xy_x_z_yy, g_xy_x_z_yz, g_xy_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_x_z_xx[i] = 2.0 * g_xy_x_z_xx[i] * a_exp;

        g_x_0_0_0_y_x_z_xy[i] = 2.0 * g_xy_x_z_xy[i] * a_exp;

        g_x_0_0_0_y_x_z_xz[i] = 2.0 * g_xy_x_z_xz[i] * a_exp;

        g_x_0_0_0_y_x_z_yy[i] = 2.0 * g_xy_x_z_yy[i] * a_exp;

        g_x_0_0_0_y_x_z_yz[i] = 2.0 * g_xy_x_z_yz[i] * a_exp;

        g_x_0_0_0_y_x_z_zz[i] = 2.0 * g_xy_x_z_zz[i] * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_0_0_y_y_x_xx, g_x_0_0_0_y_y_x_xy, g_x_0_0_0_y_y_x_xz, g_x_0_0_0_y_y_x_yy, g_x_0_0_0_y_y_x_yz, g_x_0_0_0_y_y_x_zz, g_xy_y_x_xx, g_xy_y_x_xy, g_xy_y_x_xz, g_xy_y_x_yy, g_xy_y_x_yz, g_xy_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_y_x_xx[i] = 2.0 * g_xy_y_x_xx[i] * a_exp;

        g_x_0_0_0_y_y_x_xy[i] = 2.0 * g_xy_y_x_xy[i] * a_exp;

        g_x_0_0_0_y_y_x_xz[i] = 2.0 * g_xy_y_x_xz[i] * a_exp;

        g_x_0_0_0_y_y_x_yy[i] = 2.0 * g_xy_y_x_yy[i] * a_exp;

        g_x_0_0_0_y_y_x_yz[i] = 2.0 * g_xy_y_x_yz[i] * a_exp;

        g_x_0_0_0_y_y_x_zz[i] = 2.0 * g_xy_y_x_zz[i] * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_0_0_y_y_y_xx, g_x_0_0_0_y_y_y_xy, g_x_0_0_0_y_y_y_xz, g_x_0_0_0_y_y_y_yy, g_x_0_0_0_y_y_y_yz, g_x_0_0_0_y_y_y_zz, g_xy_y_y_xx, g_xy_y_y_xy, g_xy_y_y_xz, g_xy_y_y_yy, g_xy_y_y_yz, g_xy_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_y_y_xx[i] = 2.0 * g_xy_y_y_xx[i] * a_exp;

        g_x_0_0_0_y_y_y_xy[i] = 2.0 * g_xy_y_y_xy[i] * a_exp;

        g_x_0_0_0_y_y_y_xz[i] = 2.0 * g_xy_y_y_xz[i] * a_exp;

        g_x_0_0_0_y_y_y_yy[i] = 2.0 * g_xy_y_y_yy[i] * a_exp;

        g_x_0_0_0_y_y_y_yz[i] = 2.0 * g_xy_y_y_yz[i] * a_exp;

        g_x_0_0_0_y_y_y_zz[i] = 2.0 * g_xy_y_y_zz[i] * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_0_0_y_y_z_xx, g_x_0_0_0_y_y_z_xy, g_x_0_0_0_y_y_z_xz, g_x_0_0_0_y_y_z_yy, g_x_0_0_0_y_y_z_yz, g_x_0_0_0_y_y_z_zz, g_xy_y_z_xx, g_xy_y_z_xy, g_xy_y_z_xz, g_xy_y_z_yy, g_xy_y_z_yz, g_xy_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_y_z_xx[i] = 2.0 * g_xy_y_z_xx[i] * a_exp;

        g_x_0_0_0_y_y_z_xy[i] = 2.0 * g_xy_y_z_xy[i] * a_exp;

        g_x_0_0_0_y_y_z_xz[i] = 2.0 * g_xy_y_z_xz[i] * a_exp;

        g_x_0_0_0_y_y_z_yy[i] = 2.0 * g_xy_y_z_yy[i] * a_exp;

        g_x_0_0_0_y_y_z_yz[i] = 2.0 * g_xy_y_z_yz[i] * a_exp;

        g_x_0_0_0_y_y_z_zz[i] = 2.0 * g_xy_y_z_zz[i] * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_0_0_y_z_x_xx, g_x_0_0_0_y_z_x_xy, g_x_0_0_0_y_z_x_xz, g_x_0_0_0_y_z_x_yy, g_x_0_0_0_y_z_x_yz, g_x_0_0_0_y_z_x_zz, g_xy_z_x_xx, g_xy_z_x_xy, g_xy_z_x_xz, g_xy_z_x_yy, g_xy_z_x_yz, g_xy_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_z_x_xx[i] = 2.0 * g_xy_z_x_xx[i] * a_exp;

        g_x_0_0_0_y_z_x_xy[i] = 2.0 * g_xy_z_x_xy[i] * a_exp;

        g_x_0_0_0_y_z_x_xz[i] = 2.0 * g_xy_z_x_xz[i] * a_exp;

        g_x_0_0_0_y_z_x_yy[i] = 2.0 * g_xy_z_x_yy[i] * a_exp;

        g_x_0_0_0_y_z_x_yz[i] = 2.0 * g_xy_z_x_yz[i] * a_exp;

        g_x_0_0_0_y_z_x_zz[i] = 2.0 * g_xy_z_x_zz[i] * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_0_0_y_z_y_xx, g_x_0_0_0_y_z_y_xy, g_x_0_0_0_y_z_y_xz, g_x_0_0_0_y_z_y_yy, g_x_0_0_0_y_z_y_yz, g_x_0_0_0_y_z_y_zz, g_xy_z_y_xx, g_xy_z_y_xy, g_xy_z_y_xz, g_xy_z_y_yy, g_xy_z_y_yz, g_xy_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_z_y_xx[i] = 2.0 * g_xy_z_y_xx[i] * a_exp;

        g_x_0_0_0_y_z_y_xy[i] = 2.0 * g_xy_z_y_xy[i] * a_exp;

        g_x_0_0_0_y_z_y_xz[i] = 2.0 * g_xy_z_y_xz[i] * a_exp;

        g_x_0_0_0_y_z_y_yy[i] = 2.0 * g_xy_z_y_yy[i] * a_exp;

        g_x_0_0_0_y_z_y_yz[i] = 2.0 * g_xy_z_y_yz[i] * a_exp;

        g_x_0_0_0_y_z_y_zz[i] = 2.0 * g_xy_z_y_zz[i] * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_0_0_y_z_z_xx, g_x_0_0_0_y_z_z_xy, g_x_0_0_0_y_z_z_xz, g_x_0_0_0_y_z_z_yy, g_x_0_0_0_y_z_z_yz, g_x_0_0_0_y_z_z_zz, g_xy_z_z_xx, g_xy_z_z_xy, g_xy_z_z_xz, g_xy_z_z_yy, g_xy_z_z_yz, g_xy_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_z_z_xx[i] = 2.0 * g_xy_z_z_xx[i] * a_exp;

        g_x_0_0_0_y_z_z_xy[i] = 2.0 * g_xy_z_z_xy[i] * a_exp;

        g_x_0_0_0_y_z_z_xz[i] = 2.0 * g_xy_z_z_xz[i] * a_exp;

        g_x_0_0_0_y_z_z_yy[i] = 2.0 * g_xy_z_z_yy[i] * a_exp;

        g_x_0_0_0_y_z_z_yz[i] = 2.0 * g_xy_z_z_yz[i] * a_exp;

        g_x_0_0_0_y_z_z_zz[i] = 2.0 * g_xy_z_z_zz[i] * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_0_0_0_z_x_x_xx, g_x_0_0_0_z_x_x_xy, g_x_0_0_0_z_x_x_xz, g_x_0_0_0_z_x_x_yy, g_x_0_0_0_z_x_x_yz, g_x_0_0_0_z_x_x_zz, g_xz_x_x_xx, g_xz_x_x_xy, g_xz_x_x_xz, g_xz_x_x_yy, g_xz_x_x_yz, g_xz_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_x_x_xx[i] = 2.0 * g_xz_x_x_xx[i] * a_exp;

        g_x_0_0_0_z_x_x_xy[i] = 2.0 * g_xz_x_x_xy[i] * a_exp;

        g_x_0_0_0_z_x_x_xz[i] = 2.0 * g_xz_x_x_xz[i] * a_exp;

        g_x_0_0_0_z_x_x_yy[i] = 2.0 * g_xz_x_x_yy[i] * a_exp;

        g_x_0_0_0_z_x_x_yz[i] = 2.0 * g_xz_x_x_yz[i] * a_exp;

        g_x_0_0_0_z_x_x_zz[i] = 2.0 * g_xz_x_x_zz[i] * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_0_0_0_z_x_y_xx, g_x_0_0_0_z_x_y_xy, g_x_0_0_0_z_x_y_xz, g_x_0_0_0_z_x_y_yy, g_x_0_0_0_z_x_y_yz, g_x_0_0_0_z_x_y_zz, g_xz_x_y_xx, g_xz_x_y_xy, g_xz_x_y_xz, g_xz_x_y_yy, g_xz_x_y_yz, g_xz_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_x_y_xx[i] = 2.0 * g_xz_x_y_xx[i] * a_exp;

        g_x_0_0_0_z_x_y_xy[i] = 2.0 * g_xz_x_y_xy[i] * a_exp;

        g_x_0_0_0_z_x_y_xz[i] = 2.0 * g_xz_x_y_xz[i] * a_exp;

        g_x_0_0_0_z_x_y_yy[i] = 2.0 * g_xz_x_y_yy[i] * a_exp;

        g_x_0_0_0_z_x_y_yz[i] = 2.0 * g_xz_x_y_yz[i] * a_exp;

        g_x_0_0_0_z_x_y_zz[i] = 2.0 * g_xz_x_y_zz[i] * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_0_0_0_z_x_z_xx, g_x_0_0_0_z_x_z_xy, g_x_0_0_0_z_x_z_xz, g_x_0_0_0_z_x_z_yy, g_x_0_0_0_z_x_z_yz, g_x_0_0_0_z_x_z_zz, g_xz_x_z_xx, g_xz_x_z_xy, g_xz_x_z_xz, g_xz_x_z_yy, g_xz_x_z_yz, g_xz_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_x_z_xx[i] = 2.0 * g_xz_x_z_xx[i] * a_exp;

        g_x_0_0_0_z_x_z_xy[i] = 2.0 * g_xz_x_z_xy[i] * a_exp;

        g_x_0_0_0_z_x_z_xz[i] = 2.0 * g_xz_x_z_xz[i] * a_exp;

        g_x_0_0_0_z_x_z_yy[i] = 2.0 * g_xz_x_z_yy[i] * a_exp;

        g_x_0_0_0_z_x_z_yz[i] = 2.0 * g_xz_x_z_yz[i] * a_exp;

        g_x_0_0_0_z_x_z_zz[i] = 2.0 * g_xz_x_z_zz[i] * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_0_0_0_z_y_x_xx, g_x_0_0_0_z_y_x_xy, g_x_0_0_0_z_y_x_xz, g_x_0_0_0_z_y_x_yy, g_x_0_0_0_z_y_x_yz, g_x_0_0_0_z_y_x_zz, g_xz_y_x_xx, g_xz_y_x_xy, g_xz_y_x_xz, g_xz_y_x_yy, g_xz_y_x_yz, g_xz_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_y_x_xx[i] = 2.0 * g_xz_y_x_xx[i] * a_exp;

        g_x_0_0_0_z_y_x_xy[i] = 2.0 * g_xz_y_x_xy[i] * a_exp;

        g_x_0_0_0_z_y_x_xz[i] = 2.0 * g_xz_y_x_xz[i] * a_exp;

        g_x_0_0_0_z_y_x_yy[i] = 2.0 * g_xz_y_x_yy[i] * a_exp;

        g_x_0_0_0_z_y_x_yz[i] = 2.0 * g_xz_y_x_yz[i] * a_exp;

        g_x_0_0_0_z_y_x_zz[i] = 2.0 * g_xz_y_x_zz[i] * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_0_0_0_z_y_y_xx, g_x_0_0_0_z_y_y_xy, g_x_0_0_0_z_y_y_xz, g_x_0_0_0_z_y_y_yy, g_x_0_0_0_z_y_y_yz, g_x_0_0_0_z_y_y_zz, g_xz_y_y_xx, g_xz_y_y_xy, g_xz_y_y_xz, g_xz_y_y_yy, g_xz_y_y_yz, g_xz_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_y_y_xx[i] = 2.0 * g_xz_y_y_xx[i] * a_exp;

        g_x_0_0_0_z_y_y_xy[i] = 2.0 * g_xz_y_y_xy[i] * a_exp;

        g_x_0_0_0_z_y_y_xz[i] = 2.0 * g_xz_y_y_xz[i] * a_exp;

        g_x_0_0_0_z_y_y_yy[i] = 2.0 * g_xz_y_y_yy[i] * a_exp;

        g_x_0_0_0_z_y_y_yz[i] = 2.0 * g_xz_y_y_yz[i] * a_exp;

        g_x_0_0_0_z_y_y_zz[i] = 2.0 * g_xz_y_y_zz[i] * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_0_0_0_z_y_z_xx, g_x_0_0_0_z_y_z_xy, g_x_0_0_0_z_y_z_xz, g_x_0_0_0_z_y_z_yy, g_x_0_0_0_z_y_z_yz, g_x_0_0_0_z_y_z_zz, g_xz_y_z_xx, g_xz_y_z_xy, g_xz_y_z_xz, g_xz_y_z_yy, g_xz_y_z_yz, g_xz_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_y_z_xx[i] = 2.0 * g_xz_y_z_xx[i] * a_exp;

        g_x_0_0_0_z_y_z_xy[i] = 2.0 * g_xz_y_z_xy[i] * a_exp;

        g_x_0_0_0_z_y_z_xz[i] = 2.0 * g_xz_y_z_xz[i] * a_exp;

        g_x_0_0_0_z_y_z_yy[i] = 2.0 * g_xz_y_z_yy[i] * a_exp;

        g_x_0_0_0_z_y_z_yz[i] = 2.0 * g_xz_y_z_yz[i] * a_exp;

        g_x_0_0_0_z_y_z_zz[i] = 2.0 * g_xz_y_z_zz[i] * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_0_0_z_z_x_xx, g_x_0_0_0_z_z_x_xy, g_x_0_0_0_z_z_x_xz, g_x_0_0_0_z_z_x_yy, g_x_0_0_0_z_z_x_yz, g_x_0_0_0_z_z_x_zz, g_xz_z_x_xx, g_xz_z_x_xy, g_xz_z_x_xz, g_xz_z_x_yy, g_xz_z_x_yz, g_xz_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_z_x_xx[i] = 2.0 * g_xz_z_x_xx[i] * a_exp;

        g_x_0_0_0_z_z_x_xy[i] = 2.0 * g_xz_z_x_xy[i] * a_exp;

        g_x_0_0_0_z_z_x_xz[i] = 2.0 * g_xz_z_x_xz[i] * a_exp;

        g_x_0_0_0_z_z_x_yy[i] = 2.0 * g_xz_z_x_yy[i] * a_exp;

        g_x_0_0_0_z_z_x_yz[i] = 2.0 * g_xz_z_x_yz[i] * a_exp;

        g_x_0_0_0_z_z_x_zz[i] = 2.0 * g_xz_z_x_zz[i] * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_0_0_z_z_y_xx, g_x_0_0_0_z_z_y_xy, g_x_0_0_0_z_z_y_xz, g_x_0_0_0_z_z_y_yy, g_x_0_0_0_z_z_y_yz, g_x_0_0_0_z_z_y_zz, g_xz_z_y_xx, g_xz_z_y_xy, g_xz_z_y_xz, g_xz_z_y_yy, g_xz_z_y_yz, g_xz_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_z_y_xx[i] = 2.0 * g_xz_z_y_xx[i] * a_exp;

        g_x_0_0_0_z_z_y_xy[i] = 2.0 * g_xz_z_y_xy[i] * a_exp;

        g_x_0_0_0_z_z_y_xz[i] = 2.0 * g_xz_z_y_xz[i] * a_exp;

        g_x_0_0_0_z_z_y_yy[i] = 2.0 * g_xz_z_y_yy[i] * a_exp;

        g_x_0_0_0_z_z_y_yz[i] = 2.0 * g_xz_z_y_yz[i] * a_exp;

        g_x_0_0_0_z_z_y_zz[i] = 2.0 * g_xz_z_y_zz[i] * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_0_0_z_z_z_xx, g_x_0_0_0_z_z_z_xy, g_x_0_0_0_z_z_z_xz, g_x_0_0_0_z_z_z_yy, g_x_0_0_0_z_z_z_yz, g_x_0_0_0_z_z_z_zz, g_xz_z_z_xx, g_xz_z_z_xy, g_xz_z_z_xz, g_xz_z_z_yy, g_xz_z_z_yz, g_xz_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_z_z_xx[i] = 2.0 * g_xz_z_z_xx[i] * a_exp;

        g_x_0_0_0_z_z_z_xy[i] = 2.0 * g_xz_z_z_xy[i] * a_exp;

        g_x_0_0_0_z_z_z_xz[i] = 2.0 * g_xz_z_z_xz[i] * a_exp;

        g_x_0_0_0_z_z_z_yy[i] = 2.0 * g_xz_z_z_yy[i] * a_exp;

        g_x_0_0_0_z_z_z_yz[i] = 2.0 * g_xz_z_z_yz[i] * a_exp;

        g_x_0_0_0_z_z_z_zz[i] = 2.0 * g_xz_z_z_zz[i] * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_xy_x_x_xx, g_xy_x_x_xy, g_xy_x_x_xz, g_xy_x_x_yy, g_xy_x_x_yz, g_xy_x_x_zz, g_y_0_0_0_x_x_x_xx, g_y_0_0_0_x_x_x_xy, g_y_0_0_0_x_x_x_xz, g_y_0_0_0_x_x_x_yy, g_y_0_0_0_x_x_x_yz, g_y_0_0_0_x_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_x_x_xx[i] = 2.0 * g_xy_x_x_xx[i] * a_exp;

        g_y_0_0_0_x_x_x_xy[i] = 2.0 * g_xy_x_x_xy[i] * a_exp;

        g_y_0_0_0_x_x_x_xz[i] = 2.0 * g_xy_x_x_xz[i] * a_exp;

        g_y_0_0_0_x_x_x_yy[i] = 2.0 * g_xy_x_x_yy[i] * a_exp;

        g_y_0_0_0_x_x_x_yz[i] = 2.0 * g_xy_x_x_yz[i] * a_exp;

        g_y_0_0_0_x_x_x_zz[i] = 2.0 * g_xy_x_x_zz[i] * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_xy_x_y_xx, g_xy_x_y_xy, g_xy_x_y_xz, g_xy_x_y_yy, g_xy_x_y_yz, g_xy_x_y_zz, g_y_0_0_0_x_x_y_xx, g_y_0_0_0_x_x_y_xy, g_y_0_0_0_x_x_y_xz, g_y_0_0_0_x_x_y_yy, g_y_0_0_0_x_x_y_yz, g_y_0_0_0_x_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_x_y_xx[i] = 2.0 * g_xy_x_y_xx[i] * a_exp;

        g_y_0_0_0_x_x_y_xy[i] = 2.0 * g_xy_x_y_xy[i] * a_exp;

        g_y_0_0_0_x_x_y_xz[i] = 2.0 * g_xy_x_y_xz[i] * a_exp;

        g_y_0_0_0_x_x_y_yy[i] = 2.0 * g_xy_x_y_yy[i] * a_exp;

        g_y_0_0_0_x_x_y_yz[i] = 2.0 * g_xy_x_y_yz[i] * a_exp;

        g_y_0_0_0_x_x_y_zz[i] = 2.0 * g_xy_x_y_zz[i] * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_xy_x_z_xx, g_xy_x_z_xy, g_xy_x_z_xz, g_xy_x_z_yy, g_xy_x_z_yz, g_xy_x_z_zz, g_y_0_0_0_x_x_z_xx, g_y_0_0_0_x_x_z_xy, g_y_0_0_0_x_x_z_xz, g_y_0_0_0_x_x_z_yy, g_y_0_0_0_x_x_z_yz, g_y_0_0_0_x_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_x_z_xx[i] = 2.0 * g_xy_x_z_xx[i] * a_exp;

        g_y_0_0_0_x_x_z_xy[i] = 2.0 * g_xy_x_z_xy[i] * a_exp;

        g_y_0_0_0_x_x_z_xz[i] = 2.0 * g_xy_x_z_xz[i] * a_exp;

        g_y_0_0_0_x_x_z_yy[i] = 2.0 * g_xy_x_z_yy[i] * a_exp;

        g_y_0_0_0_x_x_z_yz[i] = 2.0 * g_xy_x_z_yz[i] * a_exp;

        g_y_0_0_0_x_x_z_zz[i] = 2.0 * g_xy_x_z_zz[i] * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_xy_y_x_xx, g_xy_y_x_xy, g_xy_y_x_xz, g_xy_y_x_yy, g_xy_y_x_yz, g_xy_y_x_zz, g_y_0_0_0_x_y_x_xx, g_y_0_0_0_x_y_x_xy, g_y_0_0_0_x_y_x_xz, g_y_0_0_0_x_y_x_yy, g_y_0_0_0_x_y_x_yz, g_y_0_0_0_x_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_y_x_xx[i] = 2.0 * g_xy_y_x_xx[i] * a_exp;

        g_y_0_0_0_x_y_x_xy[i] = 2.0 * g_xy_y_x_xy[i] * a_exp;

        g_y_0_0_0_x_y_x_xz[i] = 2.0 * g_xy_y_x_xz[i] * a_exp;

        g_y_0_0_0_x_y_x_yy[i] = 2.0 * g_xy_y_x_yy[i] * a_exp;

        g_y_0_0_0_x_y_x_yz[i] = 2.0 * g_xy_y_x_yz[i] * a_exp;

        g_y_0_0_0_x_y_x_zz[i] = 2.0 * g_xy_y_x_zz[i] * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_xy_y_y_xx, g_xy_y_y_xy, g_xy_y_y_xz, g_xy_y_y_yy, g_xy_y_y_yz, g_xy_y_y_zz, g_y_0_0_0_x_y_y_xx, g_y_0_0_0_x_y_y_xy, g_y_0_0_0_x_y_y_xz, g_y_0_0_0_x_y_y_yy, g_y_0_0_0_x_y_y_yz, g_y_0_0_0_x_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_y_y_xx[i] = 2.0 * g_xy_y_y_xx[i] * a_exp;

        g_y_0_0_0_x_y_y_xy[i] = 2.0 * g_xy_y_y_xy[i] * a_exp;

        g_y_0_0_0_x_y_y_xz[i] = 2.0 * g_xy_y_y_xz[i] * a_exp;

        g_y_0_0_0_x_y_y_yy[i] = 2.0 * g_xy_y_y_yy[i] * a_exp;

        g_y_0_0_0_x_y_y_yz[i] = 2.0 * g_xy_y_y_yz[i] * a_exp;

        g_y_0_0_0_x_y_y_zz[i] = 2.0 * g_xy_y_y_zz[i] * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_xy_y_z_xx, g_xy_y_z_xy, g_xy_y_z_xz, g_xy_y_z_yy, g_xy_y_z_yz, g_xy_y_z_zz, g_y_0_0_0_x_y_z_xx, g_y_0_0_0_x_y_z_xy, g_y_0_0_0_x_y_z_xz, g_y_0_0_0_x_y_z_yy, g_y_0_0_0_x_y_z_yz, g_y_0_0_0_x_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_y_z_xx[i] = 2.0 * g_xy_y_z_xx[i] * a_exp;

        g_y_0_0_0_x_y_z_xy[i] = 2.0 * g_xy_y_z_xy[i] * a_exp;

        g_y_0_0_0_x_y_z_xz[i] = 2.0 * g_xy_y_z_xz[i] * a_exp;

        g_y_0_0_0_x_y_z_yy[i] = 2.0 * g_xy_y_z_yy[i] * a_exp;

        g_y_0_0_0_x_y_z_yz[i] = 2.0 * g_xy_y_z_yz[i] * a_exp;

        g_y_0_0_0_x_y_z_zz[i] = 2.0 * g_xy_y_z_zz[i] * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_xy_z_x_xx, g_xy_z_x_xy, g_xy_z_x_xz, g_xy_z_x_yy, g_xy_z_x_yz, g_xy_z_x_zz, g_y_0_0_0_x_z_x_xx, g_y_0_0_0_x_z_x_xy, g_y_0_0_0_x_z_x_xz, g_y_0_0_0_x_z_x_yy, g_y_0_0_0_x_z_x_yz, g_y_0_0_0_x_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_z_x_xx[i] = 2.0 * g_xy_z_x_xx[i] * a_exp;

        g_y_0_0_0_x_z_x_xy[i] = 2.0 * g_xy_z_x_xy[i] * a_exp;

        g_y_0_0_0_x_z_x_xz[i] = 2.0 * g_xy_z_x_xz[i] * a_exp;

        g_y_0_0_0_x_z_x_yy[i] = 2.0 * g_xy_z_x_yy[i] * a_exp;

        g_y_0_0_0_x_z_x_yz[i] = 2.0 * g_xy_z_x_yz[i] * a_exp;

        g_y_0_0_0_x_z_x_zz[i] = 2.0 * g_xy_z_x_zz[i] * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_xy_z_y_xx, g_xy_z_y_xy, g_xy_z_y_xz, g_xy_z_y_yy, g_xy_z_y_yz, g_xy_z_y_zz, g_y_0_0_0_x_z_y_xx, g_y_0_0_0_x_z_y_xy, g_y_0_0_0_x_z_y_xz, g_y_0_0_0_x_z_y_yy, g_y_0_0_0_x_z_y_yz, g_y_0_0_0_x_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_z_y_xx[i] = 2.0 * g_xy_z_y_xx[i] * a_exp;

        g_y_0_0_0_x_z_y_xy[i] = 2.0 * g_xy_z_y_xy[i] * a_exp;

        g_y_0_0_0_x_z_y_xz[i] = 2.0 * g_xy_z_y_xz[i] * a_exp;

        g_y_0_0_0_x_z_y_yy[i] = 2.0 * g_xy_z_y_yy[i] * a_exp;

        g_y_0_0_0_x_z_y_yz[i] = 2.0 * g_xy_z_y_yz[i] * a_exp;

        g_y_0_0_0_x_z_y_zz[i] = 2.0 * g_xy_z_y_zz[i] * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_xy_z_z_xx, g_xy_z_z_xy, g_xy_z_z_xz, g_xy_z_z_yy, g_xy_z_z_yz, g_xy_z_z_zz, g_y_0_0_0_x_z_z_xx, g_y_0_0_0_x_z_z_xy, g_y_0_0_0_x_z_z_xz, g_y_0_0_0_x_z_z_yy, g_y_0_0_0_x_z_z_yz, g_y_0_0_0_x_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_z_z_xx[i] = 2.0 * g_xy_z_z_xx[i] * a_exp;

        g_y_0_0_0_x_z_z_xy[i] = 2.0 * g_xy_z_z_xy[i] * a_exp;

        g_y_0_0_0_x_z_z_xz[i] = 2.0 * g_xy_z_z_xz[i] * a_exp;

        g_y_0_0_0_x_z_z_yy[i] = 2.0 * g_xy_z_z_yy[i] * a_exp;

        g_y_0_0_0_x_z_z_yz[i] = 2.0 * g_xy_z_z_yz[i] * a_exp;

        g_y_0_0_0_x_z_z_zz[i] = 2.0 * g_xy_z_z_zz[i] * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_y_0_0_0_y_x_x_xx, g_y_0_0_0_y_x_x_xy, g_y_0_0_0_y_x_x_xz, g_y_0_0_0_y_x_x_yy, g_y_0_0_0_y_x_x_yz, g_y_0_0_0_y_x_x_zz, g_yy_x_x_xx, g_yy_x_x_xy, g_yy_x_x_xz, g_yy_x_x_yy, g_yy_x_x_yz, g_yy_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_x_x_xx[i] = -g_0_x_x_xx[i] + 2.0 * g_yy_x_x_xx[i] * a_exp;

        g_y_0_0_0_y_x_x_xy[i] = -g_0_x_x_xy[i] + 2.0 * g_yy_x_x_xy[i] * a_exp;

        g_y_0_0_0_y_x_x_xz[i] = -g_0_x_x_xz[i] + 2.0 * g_yy_x_x_xz[i] * a_exp;

        g_y_0_0_0_y_x_x_yy[i] = -g_0_x_x_yy[i] + 2.0 * g_yy_x_x_yy[i] * a_exp;

        g_y_0_0_0_y_x_x_yz[i] = -g_0_x_x_yz[i] + 2.0 * g_yy_x_x_yz[i] * a_exp;

        g_y_0_0_0_y_x_x_zz[i] = -g_0_x_x_zz[i] + 2.0 * g_yy_x_x_zz[i] * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_y_0_0_0_y_x_y_xx, g_y_0_0_0_y_x_y_xy, g_y_0_0_0_y_x_y_xz, g_y_0_0_0_y_x_y_yy, g_y_0_0_0_y_x_y_yz, g_y_0_0_0_y_x_y_zz, g_yy_x_y_xx, g_yy_x_y_xy, g_yy_x_y_xz, g_yy_x_y_yy, g_yy_x_y_yz, g_yy_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_x_y_xx[i] = -g_0_x_y_xx[i] + 2.0 * g_yy_x_y_xx[i] * a_exp;

        g_y_0_0_0_y_x_y_xy[i] = -g_0_x_y_xy[i] + 2.0 * g_yy_x_y_xy[i] * a_exp;

        g_y_0_0_0_y_x_y_xz[i] = -g_0_x_y_xz[i] + 2.0 * g_yy_x_y_xz[i] * a_exp;

        g_y_0_0_0_y_x_y_yy[i] = -g_0_x_y_yy[i] + 2.0 * g_yy_x_y_yy[i] * a_exp;

        g_y_0_0_0_y_x_y_yz[i] = -g_0_x_y_yz[i] + 2.0 * g_yy_x_y_yz[i] * a_exp;

        g_y_0_0_0_y_x_y_zz[i] = -g_0_x_y_zz[i] + 2.0 * g_yy_x_y_zz[i] * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_y_0_0_0_y_x_z_xx, g_y_0_0_0_y_x_z_xy, g_y_0_0_0_y_x_z_xz, g_y_0_0_0_y_x_z_yy, g_y_0_0_0_y_x_z_yz, g_y_0_0_0_y_x_z_zz, g_yy_x_z_xx, g_yy_x_z_xy, g_yy_x_z_xz, g_yy_x_z_yy, g_yy_x_z_yz, g_yy_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_x_z_xx[i] = -g_0_x_z_xx[i] + 2.0 * g_yy_x_z_xx[i] * a_exp;

        g_y_0_0_0_y_x_z_xy[i] = -g_0_x_z_xy[i] + 2.0 * g_yy_x_z_xy[i] * a_exp;

        g_y_0_0_0_y_x_z_xz[i] = -g_0_x_z_xz[i] + 2.0 * g_yy_x_z_xz[i] * a_exp;

        g_y_0_0_0_y_x_z_yy[i] = -g_0_x_z_yy[i] + 2.0 * g_yy_x_z_yy[i] * a_exp;

        g_y_0_0_0_y_x_z_yz[i] = -g_0_x_z_yz[i] + 2.0 * g_yy_x_z_yz[i] * a_exp;

        g_y_0_0_0_y_x_z_zz[i] = -g_0_x_z_zz[i] + 2.0 * g_yy_x_z_zz[i] * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_y_0_0_0_y_y_x_xx, g_y_0_0_0_y_y_x_xy, g_y_0_0_0_y_y_x_xz, g_y_0_0_0_y_y_x_yy, g_y_0_0_0_y_y_x_yz, g_y_0_0_0_y_y_x_zz, g_yy_y_x_xx, g_yy_y_x_xy, g_yy_y_x_xz, g_yy_y_x_yy, g_yy_y_x_yz, g_yy_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_y_x_xx[i] = -g_0_y_x_xx[i] + 2.0 * g_yy_y_x_xx[i] * a_exp;

        g_y_0_0_0_y_y_x_xy[i] = -g_0_y_x_xy[i] + 2.0 * g_yy_y_x_xy[i] * a_exp;

        g_y_0_0_0_y_y_x_xz[i] = -g_0_y_x_xz[i] + 2.0 * g_yy_y_x_xz[i] * a_exp;

        g_y_0_0_0_y_y_x_yy[i] = -g_0_y_x_yy[i] + 2.0 * g_yy_y_x_yy[i] * a_exp;

        g_y_0_0_0_y_y_x_yz[i] = -g_0_y_x_yz[i] + 2.0 * g_yy_y_x_yz[i] * a_exp;

        g_y_0_0_0_y_y_x_zz[i] = -g_0_y_x_zz[i] + 2.0 * g_yy_y_x_zz[i] * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_y_0_0_0_y_y_y_xx, g_y_0_0_0_y_y_y_xy, g_y_0_0_0_y_y_y_xz, g_y_0_0_0_y_y_y_yy, g_y_0_0_0_y_y_y_yz, g_y_0_0_0_y_y_y_zz, g_yy_y_y_xx, g_yy_y_y_xy, g_yy_y_y_xz, g_yy_y_y_yy, g_yy_y_y_yz, g_yy_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_y_y_xx[i] = -g_0_y_y_xx[i] + 2.0 * g_yy_y_y_xx[i] * a_exp;

        g_y_0_0_0_y_y_y_xy[i] = -g_0_y_y_xy[i] + 2.0 * g_yy_y_y_xy[i] * a_exp;

        g_y_0_0_0_y_y_y_xz[i] = -g_0_y_y_xz[i] + 2.0 * g_yy_y_y_xz[i] * a_exp;

        g_y_0_0_0_y_y_y_yy[i] = -g_0_y_y_yy[i] + 2.0 * g_yy_y_y_yy[i] * a_exp;

        g_y_0_0_0_y_y_y_yz[i] = -g_0_y_y_yz[i] + 2.0 * g_yy_y_y_yz[i] * a_exp;

        g_y_0_0_0_y_y_y_zz[i] = -g_0_y_y_zz[i] + 2.0 * g_yy_y_y_zz[i] * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_y_0_0_0_y_y_z_xx, g_y_0_0_0_y_y_z_xy, g_y_0_0_0_y_y_z_xz, g_y_0_0_0_y_y_z_yy, g_y_0_0_0_y_y_z_yz, g_y_0_0_0_y_y_z_zz, g_yy_y_z_xx, g_yy_y_z_xy, g_yy_y_z_xz, g_yy_y_z_yy, g_yy_y_z_yz, g_yy_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_y_z_xx[i] = -g_0_y_z_xx[i] + 2.0 * g_yy_y_z_xx[i] * a_exp;

        g_y_0_0_0_y_y_z_xy[i] = -g_0_y_z_xy[i] + 2.0 * g_yy_y_z_xy[i] * a_exp;

        g_y_0_0_0_y_y_z_xz[i] = -g_0_y_z_xz[i] + 2.0 * g_yy_y_z_xz[i] * a_exp;

        g_y_0_0_0_y_y_z_yy[i] = -g_0_y_z_yy[i] + 2.0 * g_yy_y_z_yy[i] * a_exp;

        g_y_0_0_0_y_y_z_yz[i] = -g_0_y_z_yz[i] + 2.0 * g_yy_y_z_yz[i] * a_exp;

        g_y_0_0_0_y_y_z_zz[i] = -g_0_y_z_zz[i] + 2.0 * g_yy_y_z_zz[i] * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_y_0_0_0_y_z_x_xx, g_y_0_0_0_y_z_x_xy, g_y_0_0_0_y_z_x_xz, g_y_0_0_0_y_z_x_yy, g_y_0_0_0_y_z_x_yz, g_y_0_0_0_y_z_x_zz, g_yy_z_x_xx, g_yy_z_x_xy, g_yy_z_x_xz, g_yy_z_x_yy, g_yy_z_x_yz, g_yy_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_z_x_xx[i] = -g_0_z_x_xx[i] + 2.0 * g_yy_z_x_xx[i] * a_exp;

        g_y_0_0_0_y_z_x_xy[i] = -g_0_z_x_xy[i] + 2.0 * g_yy_z_x_xy[i] * a_exp;

        g_y_0_0_0_y_z_x_xz[i] = -g_0_z_x_xz[i] + 2.0 * g_yy_z_x_xz[i] * a_exp;

        g_y_0_0_0_y_z_x_yy[i] = -g_0_z_x_yy[i] + 2.0 * g_yy_z_x_yy[i] * a_exp;

        g_y_0_0_0_y_z_x_yz[i] = -g_0_z_x_yz[i] + 2.0 * g_yy_z_x_yz[i] * a_exp;

        g_y_0_0_0_y_z_x_zz[i] = -g_0_z_x_zz[i] + 2.0 * g_yy_z_x_zz[i] * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_y_0_0_0_y_z_y_xx, g_y_0_0_0_y_z_y_xy, g_y_0_0_0_y_z_y_xz, g_y_0_0_0_y_z_y_yy, g_y_0_0_0_y_z_y_yz, g_y_0_0_0_y_z_y_zz, g_yy_z_y_xx, g_yy_z_y_xy, g_yy_z_y_xz, g_yy_z_y_yy, g_yy_z_y_yz, g_yy_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_z_y_xx[i] = -g_0_z_y_xx[i] + 2.0 * g_yy_z_y_xx[i] * a_exp;

        g_y_0_0_0_y_z_y_xy[i] = -g_0_z_y_xy[i] + 2.0 * g_yy_z_y_xy[i] * a_exp;

        g_y_0_0_0_y_z_y_xz[i] = -g_0_z_y_xz[i] + 2.0 * g_yy_z_y_xz[i] * a_exp;

        g_y_0_0_0_y_z_y_yy[i] = -g_0_z_y_yy[i] + 2.0 * g_yy_z_y_yy[i] * a_exp;

        g_y_0_0_0_y_z_y_yz[i] = -g_0_z_y_yz[i] + 2.0 * g_yy_z_y_yz[i] * a_exp;

        g_y_0_0_0_y_z_y_zz[i] = -g_0_z_y_zz[i] + 2.0 * g_yy_z_y_zz[i] * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_y_0_0_0_y_z_z_xx, g_y_0_0_0_y_z_z_xy, g_y_0_0_0_y_z_z_xz, g_y_0_0_0_y_z_z_yy, g_y_0_0_0_y_z_z_yz, g_y_0_0_0_y_z_z_zz, g_yy_z_z_xx, g_yy_z_z_xy, g_yy_z_z_xz, g_yy_z_z_yy, g_yy_z_z_yz, g_yy_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_z_z_xx[i] = -g_0_z_z_xx[i] + 2.0 * g_yy_z_z_xx[i] * a_exp;

        g_y_0_0_0_y_z_z_xy[i] = -g_0_z_z_xy[i] + 2.0 * g_yy_z_z_xy[i] * a_exp;

        g_y_0_0_0_y_z_z_xz[i] = -g_0_z_z_xz[i] + 2.0 * g_yy_z_z_xz[i] * a_exp;

        g_y_0_0_0_y_z_z_yy[i] = -g_0_z_z_yy[i] + 2.0 * g_yy_z_z_yy[i] * a_exp;

        g_y_0_0_0_y_z_z_yz[i] = -g_0_z_z_yz[i] + 2.0 * g_yy_z_z_yz[i] * a_exp;

        g_y_0_0_0_y_z_z_zz[i] = -g_0_z_z_zz[i] + 2.0 * g_yy_z_z_zz[i] * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_y_0_0_0_z_x_x_xx, g_y_0_0_0_z_x_x_xy, g_y_0_0_0_z_x_x_xz, g_y_0_0_0_z_x_x_yy, g_y_0_0_0_z_x_x_yz, g_y_0_0_0_z_x_x_zz, g_yz_x_x_xx, g_yz_x_x_xy, g_yz_x_x_xz, g_yz_x_x_yy, g_yz_x_x_yz, g_yz_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_x_x_xx[i] = 2.0 * g_yz_x_x_xx[i] * a_exp;

        g_y_0_0_0_z_x_x_xy[i] = 2.0 * g_yz_x_x_xy[i] * a_exp;

        g_y_0_0_0_z_x_x_xz[i] = 2.0 * g_yz_x_x_xz[i] * a_exp;

        g_y_0_0_0_z_x_x_yy[i] = 2.0 * g_yz_x_x_yy[i] * a_exp;

        g_y_0_0_0_z_x_x_yz[i] = 2.0 * g_yz_x_x_yz[i] * a_exp;

        g_y_0_0_0_z_x_x_zz[i] = 2.0 * g_yz_x_x_zz[i] * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_y_0_0_0_z_x_y_xx, g_y_0_0_0_z_x_y_xy, g_y_0_0_0_z_x_y_xz, g_y_0_0_0_z_x_y_yy, g_y_0_0_0_z_x_y_yz, g_y_0_0_0_z_x_y_zz, g_yz_x_y_xx, g_yz_x_y_xy, g_yz_x_y_xz, g_yz_x_y_yy, g_yz_x_y_yz, g_yz_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_x_y_xx[i] = 2.0 * g_yz_x_y_xx[i] * a_exp;

        g_y_0_0_0_z_x_y_xy[i] = 2.0 * g_yz_x_y_xy[i] * a_exp;

        g_y_0_0_0_z_x_y_xz[i] = 2.0 * g_yz_x_y_xz[i] * a_exp;

        g_y_0_0_0_z_x_y_yy[i] = 2.0 * g_yz_x_y_yy[i] * a_exp;

        g_y_0_0_0_z_x_y_yz[i] = 2.0 * g_yz_x_y_yz[i] * a_exp;

        g_y_0_0_0_z_x_y_zz[i] = 2.0 * g_yz_x_y_zz[i] * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_y_0_0_0_z_x_z_xx, g_y_0_0_0_z_x_z_xy, g_y_0_0_0_z_x_z_xz, g_y_0_0_0_z_x_z_yy, g_y_0_0_0_z_x_z_yz, g_y_0_0_0_z_x_z_zz, g_yz_x_z_xx, g_yz_x_z_xy, g_yz_x_z_xz, g_yz_x_z_yy, g_yz_x_z_yz, g_yz_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_x_z_xx[i] = 2.0 * g_yz_x_z_xx[i] * a_exp;

        g_y_0_0_0_z_x_z_xy[i] = 2.0 * g_yz_x_z_xy[i] * a_exp;

        g_y_0_0_0_z_x_z_xz[i] = 2.0 * g_yz_x_z_xz[i] * a_exp;

        g_y_0_0_0_z_x_z_yy[i] = 2.0 * g_yz_x_z_yy[i] * a_exp;

        g_y_0_0_0_z_x_z_yz[i] = 2.0 * g_yz_x_z_yz[i] * a_exp;

        g_y_0_0_0_z_x_z_zz[i] = 2.0 * g_yz_x_z_zz[i] * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_y_0_0_0_z_y_x_xx, g_y_0_0_0_z_y_x_xy, g_y_0_0_0_z_y_x_xz, g_y_0_0_0_z_y_x_yy, g_y_0_0_0_z_y_x_yz, g_y_0_0_0_z_y_x_zz, g_yz_y_x_xx, g_yz_y_x_xy, g_yz_y_x_xz, g_yz_y_x_yy, g_yz_y_x_yz, g_yz_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_y_x_xx[i] = 2.0 * g_yz_y_x_xx[i] * a_exp;

        g_y_0_0_0_z_y_x_xy[i] = 2.0 * g_yz_y_x_xy[i] * a_exp;

        g_y_0_0_0_z_y_x_xz[i] = 2.0 * g_yz_y_x_xz[i] * a_exp;

        g_y_0_0_0_z_y_x_yy[i] = 2.0 * g_yz_y_x_yy[i] * a_exp;

        g_y_0_0_0_z_y_x_yz[i] = 2.0 * g_yz_y_x_yz[i] * a_exp;

        g_y_0_0_0_z_y_x_zz[i] = 2.0 * g_yz_y_x_zz[i] * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_y_0_0_0_z_y_y_xx, g_y_0_0_0_z_y_y_xy, g_y_0_0_0_z_y_y_xz, g_y_0_0_0_z_y_y_yy, g_y_0_0_0_z_y_y_yz, g_y_0_0_0_z_y_y_zz, g_yz_y_y_xx, g_yz_y_y_xy, g_yz_y_y_xz, g_yz_y_y_yy, g_yz_y_y_yz, g_yz_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_y_y_xx[i] = 2.0 * g_yz_y_y_xx[i] * a_exp;

        g_y_0_0_0_z_y_y_xy[i] = 2.0 * g_yz_y_y_xy[i] * a_exp;

        g_y_0_0_0_z_y_y_xz[i] = 2.0 * g_yz_y_y_xz[i] * a_exp;

        g_y_0_0_0_z_y_y_yy[i] = 2.0 * g_yz_y_y_yy[i] * a_exp;

        g_y_0_0_0_z_y_y_yz[i] = 2.0 * g_yz_y_y_yz[i] * a_exp;

        g_y_0_0_0_z_y_y_zz[i] = 2.0 * g_yz_y_y_zz[i] * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_y_0_0_0_z_y_z_xx, g_y_0_0_0_z_y_z_xy, g_y_0_0_0_z_y_z_xz, g_y_0_0_0_z_y_z_yy, g_y_0_0_0_z_y_z_yz, g_y_0_0_0_z_y_z_zz, g_yz_y_z_xx, g_yz_y_z_xy, g_yz_y_z_xz, g_yz_y_z_yy, g_yz_y_z_yz, g_yz_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_y_z_xx[i] = 2.0 * g_yz_y_z_xx[i] * a_exp;

        g_y_0_0_0_z_y_z_xy[i] = 2.0 * g_yz_y_z_xy[i] * a_exp;

        g_y_0_0_0_z_y_z_xz[i] = 2.0 * g_yz_y_z_xz[i] * a_exp;

        g_y_0_0_0_z_y_z_yy[i] = 2.0 * g_yz_y_z_yy[i] * a_exp;

        g_y_0_0_0_z_y_z_yz[i] = 2.0 * g_yz_y_z_yz[i] * a_exp;

        g_y_0_0_0_z_y_z_zz[i] = 2.0 * g_yz_y_z_zz[i] * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_y_0_0_0_z_z_x_xx, g_y_0_0_0_z_z_x_xy, g_y_0_0_0_z_z_x_xz, g_y_0_0_0_z_z_x_yy, g_y_0_0_0_z_z_x_yz, g_y_0_0_0_z_z_x_zz, g_yz_z_x_xx, g_yz_z_x_xy, g_yz_z_x_xz, g_yz_z_x_yy, g_yz_z_x_yz, g_yz_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_z_x_xx[i] = 2.0 * g_yz_z_x_xx[i] * a_exp;

        g_y_0_0_0_z_z_x_xy[i] = 2.0 * g_yz_z_x_xy[i] * a_exp;

        g_y_0_0_0_z_z_x_xz[i] = 2.0 * g_yz_z_x_xz[i] * a_exp;

        g_y_0_0_0_z_z_x_yy[i] = 2.0 * g_yz_z_x_yy[i] * a_exp;

        g_y_0_0_0_z_z_x_yz[i] = 2.0 * g_yz_z_x_yz[i] * a_exp;

        g_y_0_0_0_z_z_x_zz[i] = 2.0 * g_yz_z_x_zz[i] * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_y_0_0_0_z_z_y_xx, g_y_0_0_0_z_z_y_xy, g_y_0_0_0_z_z_y_xz, g_y_0_0_0_z_z_y_yy, g_y_0_0_0_z_z_y_yz, g_y_0_0_0_z_z_y_zz, g_yz_z_y_xx, g_yz_z_y_xy, g_yz_z_y_xz, g_yz_z_y_yy, g_yz_z_y_yz, g_yz_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_z_y_xx[i] = 2.0 * g_yz_z_y_xx[i] * a_exp;

        g_y_0_0_0_z_z_y_xy[i] = 2.0 * g_yz_z_y_xy[i] * a_exp;

        g_y_0_0_0_z_z_y_xz[i] = 2.0 * g_yz_z_y_xz[i] * a_exp;

        g_y_0_0_0_z_z_y_yy[i] = 2.0 * g_yz_z_y_yy[i] * a_exp;

        g_y_0_0_0_z_z_y_yz[i] = 2.0 * g_yz_z_y_yz[i] * a_exp;

        g_y_0_0_0_z_z_y_zz[i] = 2.0 * g_yz_z_y_zz[i] * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_y_0_0_0_z_z_z_xx, g_y_0_0_0_z_z_z_xy, g_y_0_0_0_z_z_z_xz, g_y_0_0_0_z_z_z_yy, g_y_0_0_0_z_z_z_yz, g_y_0_0_0_z_z_z_zz, g_yz_z_z_xx, g_yz_z_z_xy, g_yz_z_z_xz, g_yz_z_z_yy, g_yz_z_z_yz, g_yz_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_z_z_xx[i] = 2.0 * g_yz_z_z_xx[i] * a_exp;

        g_y_0_0_0_z_z_z_xy[i] = 2.0 * g_yz_z_z_xy[i] * a_exp;

        g_y_0_0_0_z_z_z_xz[i] = 2.0 * g_yz_z_z_xz[i] * a_exp;

        g_y_0_0_0_z_z_z_yy[i] = 2.0 * g_yz_z_z_yy[i] * a_exp;

        g_y_0_0_0_z_z_z_yz[i] = 2.0 * g_yz_z_z_yz[i] * a_exp;

        g_y_0_0_0_z_z_z_zz[i] = 2.0 * g_yz_z_z_zz[i] * a_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_xz_x_x_xx, g_xz_x_x_xy, g_xz_x_x_xz, g_xz_x_x_yy, g_xz_x_x_yz, g_xz_x_x_zz, g_z_0_0_0_x_x_x_xx, g_z_0_0_0_x_x_x_xy, g_z_0_0_0_x_x_x_xz, g_z_0_0_0_x_x_x_yy, g_z_0_0_0_x_x_x_yz, g_z_0_0_0_x_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_x_x_xx[i] = 2.0 * g_xz_x_x_xx[i] * a_exp;

        g_z_0_0_0_x_x_x_xy[i] = 2.0 * g_xz_x_x_xy[i] * a_exp;

        g_z_0_0_0_x_x_x_xz[i] = 2.0 * g_xz_x_x_xz[i] * a_exp;

        g_z_0_0_0_x_x_x_yy[i] = 2.0 * g_xz_x_x_yy[i] * a_exp;

        g_z_0_0_0_x_x_x_yz[i] = 2.0 * g_xz_x_x_yz[i] * a_exp;

        g_z_0_0_0_x_x_x_zz[i] = 2.0 * g_xz_x_x_zz[i] * a_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_xz_x_y_xx, g_xz_x_y_xy, g_xz_x_y_xz, g_xz_x_y_yy, g_xz_x_y_yz, g_xz_x_y_zz, g_z_0_0_0_x_x_y_xx, g_z_0_0_0_x_x_y_xy, g_z_0_0_0_x_x_y_xz, g_z_0_0_0_x_x_y_yy, g_z_0_0_0_x_x_y_yz, g_z_0_0_0_x_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_x_y_xx[i] = 2.0 * g_xz_x_y_xx[i] * a_exp;

        g_z_0_0_0_x_x_y_xy[i] = 2.0 * g_xz_x_y_xy[i] * a_exp;

        g_z_0_0_0_x_x_y_xz[i] = 2.0 * g_xz_x_y_xz[i] * a_exp;

        g_z_0_0_0_x_x_y_yy[i] = 2.0 * g_xz_x_y_yy[i] * a_exp;

        g_z_0_0_0_x_x_y_yz[i] = 2.0 * g_xz_x_y_yz[i] * a_exp;

        g_z_0_0_0_x_x_y_zz[i] = 2.0 * g_xz_x_y_zz[i] * a_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_xz_x_z_xx, g_xz_x_z_xy, g_xz_x_z_xz, g_xz_x_z_yy, g_xz_x_z_yz, g_xz_x_z_zz, g_z_0_0_0_x_x_z_xx, g_z_0_0_0_x_x_z_xy, g_z_0_0_0_x_x_z_xz, g_z_0_0_0_x_x_z_yy, g_z_0_0_0_x_x_z_yz, g_z_0_0_0_x_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_x_z_xx[i] = 2.0 * g_xz_x_z_xx[i] * a_exp;

        g_z_0_0_0_x_x_z_xy[i] = 2.0 * g_xz_x_z_xy[i] * a_exp;

        g_z_0_0_0_x_x_z_xz[i] = 2.0 * g_xz_x_z_xz[i] * a_exp;

        g_z_0_0_0_x_x_z_yy[i] = 2.0 * g_xz_x_z_yy[i] * a_exp;

        g_z_0_0_0_x_x_z_yz[i] = 2.0 * g_xz_x_z_yz[i] * a_exp;

        g_z_0_0_0_x_x_z_zz[i] = 2.0 * g_xz_x_z_zz[i] * a_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_xz_y_x_xx, g_xz_y_x_xy, g_xz_y_x_xz, g_xz_y_x_yy, g_xz_y_x_yz, g_xz_y_x_zz, g_z_0_0_0_x_y_x_xx, g_z_0_0_0_x_y_x_xy, g_z_0_0_0_x_y_x_xz, g_z_0_0_0_x_y_x_yy, g_z_0_0_0_x_y_x_yz, g_z_0_0_0_x_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_y_x_xx[i] = 2.0 * g_xz_y_x_xx[i] * a_exp;

        g_z_0_0_0_x_y_x_xy[i] = 2.0 * g_xz_y_x_xy[i] * a_exp;

        g_z_0_0_0_x_y_x_xz[i] = 2.0 * g_xz_y_x_xz[i] * a_exp;

        g_z_0_0_0_x_y_x_yy[i] = 2.0 * g_xz_y_x_yy[i] * a_exp;

        g_z_0_0_0_x_y_x_yz[i] = 2.0 * g_xz_y_x_yz[i] * a_exp;

        g_z_0_0_0_x_y_x_zz[i] = 2.0 * g_xz_y_x_zz[i] * a_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_xz_y_y_xx, g_xz_y_y_xy, g_xz_y_y_xz, g_xz_y_y_yy, g_xz_y_y_yz, g_xz_y_y_zz, g_z_0_0_0_x_y_y_xx, g_z_0_0_0_x_y_y_xy, g_z_0_0_0_x_y_y_xz, g_z_0_0_0_x_y_y_yy, g_z_0_0_0_x_y_y_yz, g_z_0_0_0_x_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_y_y_xx[i] = 2.0 * g_xz_y_y_xx[i] * a_exp;

        g_z_0_0_0_x_y_y_xy[i] = 2.0 * g_xz_y_y_xy[i] * a_exp;

        g_z_0_0_0_x_y_y_xz[i] = 2.0 * g_xz_y_y_xz[i] * a_exp;

        g_z_0_0_0_x_y_y_yy[i] = 2.0 * g_xz_y_y_yy[i] * a_exp;

        g_z_0_0_0_x_y_y_yz[i] = 2.0 * g_xz_y_y_yz[i] * a_exp;

        g_z_0_0_0_x_y_y_zz[i] = 2.0 * g_xz_y_y_zz[i] * a_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_xz_y_z_xx, g_xz_y_z_xy, g_xz_y_z_xz, g_xz_y_z_yy, g_xz_y_z_yz, g_xz_y_z_zz, g_z_0_0_0_x_y_z_xx, g_z_0_0_0_x_y_z_xy, g_z_0_0_0_x_y_z_xz, g_z_0_0_0_x_y_z_yy, g_z_0_0_0_x_y_z_yz, g_z_0_0_0_x_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_y_z_xx[i] = 2.0 * g_xz_y_z_xx[i] * a_exp;

        g_z_0_0_0_x_y_z_xy[i] = 2.0 * g_xz_y_z_xy[i] * a_exp;

        g_z_0_0_0_x_y_z_xz[i] = 2.0 * g_xz_y_z_xz[i] * a_exp;

        g_z_0_0_0_x_y_z_yy[i] = 2.0 * g_xz_y_z_yy[i] * a_exp;

        g_z_0_0_0_x_y_z_yz[i] = 2.0 * g_xz_y_z_yz[i] * a_exp;

        g_z_0_0_0_x_y_z_zz[i] = 2.0 * g_xz_y_z_zz[i] * a_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_xz_z_x_xx, g_xz_z_x_xy, g_xz_z_x_xz, g_xz_z_x_yy, g_xz_z_x_yz, g_xz_z_x_zz, g_z_0_0_0_x_z_x_xx, g_z_0_0_0_x_z_x_xy, g_z_0_0_0_x_z_x_xz, g_z_0_0_0_x_z_x_yy, g_z_0_0_0_x_z_x_yz, g_z_0_0_0_x_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_z_x_xx[i] = 2.0 * g_xz_z_x_xx[i] * a_exp;

        g_z_0_0_0_x_z_x_xy[i] = 2.0 * g_xz_z_x_xy[i] * a_exp;

        g_z_0_0_0_x_z_x_xz[i] = 2.0 * g_xz_z_x_xz[i] * a_exp;

        g_z_0_0_0_x_z_x_yy[i] = 2.0 * g_xz_z_x_yy[i] * a_exp;

        g_z_0_0_0_x_z_x_yz[i] = 2.0 * g_xz_z_x_yz[i] * a_exp;

        g_z_0_0_0_x_z_x_zz[i] = 2.0 * g_xz_z_x_zz[i] * a_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_xz_z_y_xx, g_xz_z_y_xy, g_xz_z_y_xz, g_xz_z_y_yy, g_xz_z_y_yz, g_xz_z_y_zz, g_z_0_0_0_x_z_y_xx, g_z_0_0_0_x_z_y_xy, g_z_0_0_0_x_z_y_xz, g_z_0_0_0_x_z_y_yy, g_z_0_0_0_x_z_y_yz, g_z_0_0_0_x_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_z_y_xx[i] = 2.0 * g_xz_z_y_xx[i] * a_exp;

        g_z_0_0_0_x_z_y_xy[i] = 2.0 * g_xz_z_y_xy[i] * a_exp;

        g_z_0_0_0_x_z_y_xz[i] = 2.0 * g_xz_z_y_xz[i] * a_exp;

        g_z_0_0_0_x_z_y_yy[i] = 2.0 * g_xz_z_y_yy[i] * a_exp;

        g_z_0_0_0_x_z_y_yz[i] = 2.0 * g_xz_z_y_yz[i] * a_exp;

        g_z_0_0_0_x_z_y_zz[i] = 2.0 * g_xz_z_y_zz[i] * a_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_xz_z_z_xx, g_xz_z_z_xy, g_xz_z_z_xz, g_xz_z_z_yy, g_xz_z_z_yz, g_xz_z_z_zz, g_z_0_0_0_x_z_z_xx, g_z_0_0_0_x_z_z_xy, g_z_0_0_0_x_z_z_xz, g_z_0_0_0_x_z_z_yy, g_z_0_0_0_x_z_z_yz, g_z_0_0_0_x_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_z_z_xx[i] = 2.0 * g_xz_z_z_xx[i] * a_exp;

        g_z_0_0_0_x_z_z_xy[i] = 2.0 * g_xz_z_z_xy[i] * a_exp;

        g_z_0_0_0_x_z_z_xz[i] = 2.0 * g_xz_z_z_xz[i] * a_exp;

        g_z_0_0_0_x_z_z_yy[i] = 2.0 * g_xz_z_z_yy[i] * a_exp;

        g_z_0_0_0_x_z_z_yz[i] = 2.0 * g_xz_z_z_yz[i] * a_exp;

        g_z_0_0_0_x_z_z_zz[i] = 2.0 * g_xz_z_z_zz[i] * a_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_yz_x_x_xx, g_yz_x_x_xy, g_yz_x_x_xz, g_yz_x_x_yy, g_yz_x_x_yz, g_yz_x_x_zz, g_z_0_0_0_y_x_x_xx, g_z_0_0_0_y_x_x_xy, g_z_0_0_0_y_x_x_xz, g_z_0_0_0_y_x_x_yy, g_z_0_0_0_y_x_x_yz, g_z_0_0_0_y_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_x_x_xx[i] = 2.0 * g_yz_x_x_xx[i] * a_exp;

        g_z_0_0_0_y_x_x_xy[i] = 2.0 * g_yz_x_x_xy[i] * a_exp;

        g_z_0_0_0_y_x_x_xz[i] = 2.0 * g_yz_x_x_xz[i] * a_exp;

        g_z_0_0_0_y_x_x_yy[i] = 2.0 * g_yz_x_x_yy[i] * a_exp;

        g_z_0_0_0_y_x_x_yz[i] = 2.0 * g_yz_x_x_yz[i] * a_exp;

        g_z_0_0_0_y_x_x_zz[i] = 2.0 * g_yz_x_x_zz[i] * a_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_yz_x_y_xx, g_yz_x_y_xy, g_yz_x_y_xz, g_yz_x_y_yy, g_yz_x_y_yz, g_yz_x_y_zz, g_z_0_0_0_y_x_y_xx, g_z_0_0_0_y_x_y_xy, g_z_0_0_0_y_x_y_xz, g_z_0_0_0_y_x_y_yy, g_z_0_0_0_y_x_y_yz, g_z_0_0_0_y_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_x_y_xx[i] = 2.0 * g_yz_x_y_xx[i] * a_exp;

        g_z_0_0_0_y_x_y_xy[i] = 2.0 * g_yz_x_y_xy[i] * a_exp;

        g_z_0_0_0_y_x_y_xz[i] = 2.0 * g_yz_x_y_xz[i] * a_exp;

        g_z_0_0_0_y_x_y_yy[i] = 2.0 * g_yz_x_y_yy[i] * a_exp;

        g_z_0_0_0_y_x_y_yz[i] = 2.0 * g_yz_x_y_yz[i] * a_exp;

        g_z_0_0_0_y_x_y_zz[i] = 2.0 * g_yz_x_y_zz[i] * a_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_yz_x_z_xx, g_yz_x_z_xy, g_yz_x_z_xz, g_yz_x_z_yy, g_yz_x_z_yz, g_yz_x_z_zz, g_z_0_0_0_y_x_z_xx, g_z_0_0_0_y_x_z_xy, g_z_0_0_0_y_x_z_xz, g_z_0_0_0_y_x_z_yy, g_z_0_0_0_y_x_z_yz, g_z_0_0_0_y_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_x_z_xx[i] = 2.0 * g_yz_x_z_xx[i] * a_exp;

        g_z_0_0_0_y_x_z_xy[i] = 2.0 * g_yz_x_z_xy[i] * a_exp;

        g_z_0_0_0_y_x_z_xz[i] = 2.0 * g_yz_x_z_xz[i] * a_exp;

        g_z_0_0_0_y_x_z_yy[i] = 2.0 * g_yz_x_z_yy[i] * a_exp;

        g_z_0_0_0_y_x_z_yz[i] = 2.0 * g_yz_x_z_yz[i] * a_exp;

        g_z_0_0_0_y_x_z_zz[i] = 2.0 * g_yz_x_z_zz[i] * a_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_yz_y_x_xx, g_yz_y_x_xy, g_yz_y_x_xz, g_yz_y_x_yy, g_yz_y_x_yz, g_yz_y_x_zz, g_z_0_0_0_y_y_x_xx, g_z_0_0_0_y_y_x_xy, g_z_0_0_0_y_y_x_xz, g_z_0_0_0_y_y_x_yy, g_z_0_0_0_y_y_x_yz, g_z_0_0_0_y_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_y_x_xx[i] = 2.0 * g_yz_y_x_xx[i] * a_exp;

        g_z_0_0_0_y_y_x_xy[i] = 2.0 * g_yz_y_x_xy[i] * a_exp;

        g_z_0_0_0_y_y_x_xz[i] = 2.0 * g_yz_y_x_xz[i] * a_exp;

        g_z_0_0_0_y_y_x_yy[i] = 2.0 * g_yz_y_x_yy[i] * a_exp;

        g_z_0_0_0_y_y_x_yz[i] = 2.0 * g_yz_y_x_yz[i] * a_exp;

        g_z_0_0_0_y_y_x_zz[i] = 2.0 * g_yz_y_x_zz[i] * a_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_yz_y_y_xx, g_yz_y_y_xy, g_yz_y_y_xz, g_yz_y_y_yy, g_yz_y_y_yz, g_yz_y_y_zz, g_z_0_0_0_y_y_y_xx, g_z_0_0_0_y_y_y_xy, g_z_0_0_0_y_y_y_xz, g_z_0_0_0_y_y_y_yy, g_z_0_0_0_y_y_y_yz, g_z_0_0_0_y_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_y_y_xx[i] = 2.0 * g_yz_y_y_xx[i] * a_exp;

        g_z_0_0_0_y_y_y_xy[i] = 2.0 * g_yz_y_y_xy[i] * a_exp;

        g_z_0_0_0_y_y_y_xz[i] = 2.0 * g_yz_y_y_xz[i] * a_exp;

        g_z_0_0_0_y_y_y_yy[i] = 2.0 * g_yz_y_y_yy[i] * a_exp;

        g_z_0_0_0_y_y_y_yz[i] = 2.0 * g_yz_y_y_yz[i] * a_exp;

        g_z_0_0_0_y_y_y_zz[i] = 2.0 * g_yz_y_y_zz[i] * a_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_yz_y_z_xx, g_yz_y_z_xy, g_yz_y_z_xz, g_yz_y_z_yy, g_yz_y_z_yz, g_yz_y_z_zz, g_z_0_0_0_y_y_z_xx, g_z_0_0_0_y_y_z_xy, g_z_0_0_0_y_y_z_xz, g_z_0_0_0_y_y_z_yy, g_z_0_0_0_y_y_z_yz, g_z_0_0_0_y_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_y_z_xx[i] = 2.0 * g_yz_y_z_xx[i] * a_exp;

        g_z_0_0_0_y_y_z_xy[i] = 2.0 * g_yz_y_z_xy[i] * a_exp;

        g_z_0_0_0_y_y_z_xz[i] = 2.0 * g_yz_y_z_xz[i] * a_exp;

        g_z_0_0_0_y_y_z_yy[i] = 2.0 * g_yz_y_z_yy[i] * a_exp;

        g_z_0_0_0_y_y_z_yz[i] = 2.0 * g_yz_y_z_yz[i] * a_exp;

        g_z_0_0_0_y_y_z_zz[i] = 2.0 * g_yz_y_z_zz[i] * a_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_yz_z_x_xx, g_yz_z_x_xy, g_yz_z_x_xz, g_yz_z_x_yy, g_yz_z_x_yz, g_yz_z_x_zz, g_z_0_0_0_y_z_x_xx, g_z_0_0_0_y_z_x_xy, g_z_0_0_0_y_z_x_xz, g_z_0_0_0_y_z_x_yy, g_z_0_0_0_y_z_x_yz, g_z_0_0_0_y_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_z_x_xx[i] = 2.0 * g_yz_z_x_xx[i] * a_exp;

        g_z_0_0_0_y_z_x_xy[i] = 2.0 * g_yz_z_x_xy[i] * a_exp;

        g_z_0_0_0_y_z_x_xz[i] = 2.0 * g_yz_z_x_xz[i] * a_exp;

        g_z_0_0_0_y_z_x_yy[i] = 2.0 * g_yz_z_x_yy[i] * a_exp;

        g_z_0_0_0_y_z_x_yz[i] = 2.0 * g_yz_z_x_yz[i] * a_exp;

        g_z_0_0_0_y_z_x_zz[i] = 2.0 * g_yz_z_x_zz[i] * a_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_yz_z_y_xx, g_yz_z_y_xy, g_yz_z_y_xz, g_yz_z_y_yy, g_yz_z_y_yz, g_yz_z_y_zz, g_z_0_0_0_y_z_y_xx, g_z_0_0_0_y_z_y_xy, g_z_0_0_0_y_z_y_xz, g_z_0_0_0_y_z_y_yy, g_z_0_0_0_y_z_y_yz, g_z_0_0_0_y_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_z_y_xx[i] = 2.0 * g_yz_z_y_xx[i] * a_exp;

        g_z_0_0_0_y_z_y_xy[i] = 2.0 * g_yz_z_y_xy[i] * a_exp;

        g_z_0_0_0_y_z_y_xz[i] = 2.0 * g_yz_z_y_xz[i] * a_exp;

        g_z_0_0_0_y_z_y_yy[i] = 2.0 * g_yz_z_y_yy[i] * a_exp;

        g_z_0_0_0_y_z_y_yz[i] = 2.0 * g_yz_z_y_yz[i] * a_exp;

        g_z_0_0_0_y_z_y_zz[i] = 2.0 * g_yz_z_y_zz[i] * a_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_yz_z_z_xx, g_yz_z_z_xy, g_yz_z_z_xz, g_yz_z_z_yy, g_yz_z_z_yz, g_yz_z_z_zz, g_z_0_0_0_y_z_z_xx, g_z_0_0_0_y_z_z_xy, g_z_0_0_0_y_z_z_xz, g_z_0_0_0_y_z_z_yy, g_z_0_0_0_y_z_z_yz, g_z_0_0_0_y_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_z_z_xx[i] = 2.0 * g_yz_z_z_xx[i] * a_exp;

        g_z_0_0_0_y_z_z_xy[i] = 2.0 * g_yz_z_z_xy[i] * a_exp;

        g_z_0_0_0_y_z_z_xz[i] = 2.0 * g_yz_z_z_xz[i] * a_exp;

        g_z_0_0_0_y_z_z_yy[i] = 2.0 * g_yz_z_z_yy[i] * a_exp;

        g_z_0_0_0_y_z_z_yz[i] = 2.0 * g_yz_z_z_yz[i] * a_exp;

        g_z_0_0_0_y_z_z_zz[i] = 2.0 * g_yz_z_z_zz[i] * a_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_z_0_0_0_z_x_x_xx, g_z_0_0_0_z_x_x_xy, g_z_0_0_0_z_x_x_xz, g_z_0_0_0_z_x_x_yy, g_z_0_0_0_z_x_x_yz, g_z_0_0_0_z_x_x_zz, g_zz_x_x_xx, g_zz_x_x_xy, g_zz_x_x_xz, g_zz_x_x_yy, g_zz_x_x_yz, g_zz_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_x_x_xx[i] = -g_0_x_x_xx[i] + 2.0 * g_zz_x_x_xx[i] * a_exp;

        g_z_0_0_0_z_x_x_xy[i] = -g_0_x_x_xy[i] + 2.0 * g_zz_x_x_xy[i] * a_exp;

        g_z_0_0_0_z_x_x_xz[i] = -g_0_x_x_xz[i] + 2.0 * g_zz_x_x_xz[i] * a_exp;

        g_z_0_0_0_z_x_x_yy[i] = -g_0_x_x_yy[i] + 2.0 * g_zz_x_x_yy[i] * a_exp;

        g_z_0_0_0_z_x_x_yz[i] = -g_0_x_x_yz[i] + 2.0 * g_zz_x_x_yz[i] * a_exp;

        g_z_0_0_0_z_x_x_zz[i] = -g_0_x_x_zz[i] + 2.0 * g_zz_x_x_zz[i] * a_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_z_0_0_0_z_x_y_xx, g_z_0_0_0_z_x_y_xy, g_z_0_0_0_z_x_y_xz, g_z_0_0_0_z_x_y_yy, g_z_0_0_0_z_x_y_yz, g_z_0_0_0_z_x_y_zz, g_zz_x_y_xx, g_zz_x_y_xy, g_zz_x_y_xz, g_zz_x_y_yy, g_zz_x_y_yz, g_zz_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_x_y_xx[i] = -g_0_x_y_xx[i] + 2.0 * g_zz_x_y_xx[i] * a_exp;

        g_z_0_0_0_z_x_y_xy[i] = -g_0_x_y_xy[i] + 2.0 * g_zz_x_y_xy[i] * a_exp;

        g_z_0_0_0_z_x_y_xz[i] = -g_0_x_y_xz[i] + 2.0 * g_zz_x_y_xz[i] * a_exp;

        g_z_0_0_0_z_x_y_yy[i] = -g_0_x_y_yy[i] + 2.0 * g_zz_x_y_yy[i] * a_exp;

        g_z_0_0_0_z_x_y_yz[i] = -g_0_x_y_yz[i] + 2.0 * g_zz_x_y_yz[i] * a_exp;

        g_z_0_0_0_z_x_y_zz[i] = -g_0_x_y_zz[i] + 2.0 * g_zz_x_y_zz[i] * a_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_z_0_0_0_z_x_z_xx, g_z_0_0_0_z_x_z_xy, g_z_0_0_0_z_x_z_xz, g_z_0_0_0_z_x_z_yy, g_z_0_0_0_z_x_z_yz, g_z_0_0_0_z_x_z_zz, g_zz_x_z_xx, g_zz_x_z_xy, g_zz_x_z_xz, g_zz_x_z_yy, g_zz_x_z_yz, g_zz_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_x_z_xx[i] = -g_0_x_z_xx[i] + 2.0 * g_zz_x_z_xx[i] * a_exp;

        g_z_0_0_0_z_x_z_xy[i] = -g_0_x_z_xy[i] + 2.0 * g_zz_x_z_xy[i] * a_exp;

        g_z_0_0_0_z_x_z_xz[i] = -g_0_x_z_xz[i] + 2.0 * g_zz_x_z_xz[i] * a_exp;

        g_z_0_0_0_z_x_z_yy[i] = -g_0_x_z_yy[i] + 2.0 * g_zz_x_z_yy[i] * a_exp;

        g_z_0_0_0_z_x_z_yz[i] = -g_0_x_z_yz[i] + 2.0 * g_zz_x_z_yz[i] * a_exp;

        g_z_0_0_0_z_x_z_zz[i] = -g_0_x_z_zz[i] + 2.0 * g_zz_x_z_zz[i] * a_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_z_0_0_0_z_y_x_xx, g_z_0_0_0_z_y_x_xy, g_z_0_0_0_z_y_x_xz, g_z_0_0_0_z_y_x_yy, g_z_0_0_0_z_y_x_yz, g_z_0_0_0_z_y_x_zz, g_zz_y_x_xx, g_zz_y_x_xy, g_zz_y_x_xz, g_zz_y_x_yy, g_zz_y_x_yz, g_zz_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_y_x_xx[i] = -g_0_y_x_xx[i] + 2.0 * g_zz_y_x_xx[i] * a_exp;

        g_z_0_0_0_z_y_x_xy[i] = -g_0_y_x_xy[i] + 2.0 * g_zz_y_x_xy[i] * a_exp;

        g_z_0_0_0_z_y_x_xz[i] = -g_0_y_x_xz[i] + 2.0 * g_zz_y_x_xz[i] * a_exp;

        g_z_0_0_0_z_y_x_yy[i] = -g_0_y_x_yy[i] + 2.0 * g_zz_y_x_yy[i] * a_exp;

        g_z_0_0_0_z_y_x_yz[i] = -g_0_y_x_yz[i] + 2.0 * g_zz_y_x_yz[i] * a_exp;

        g_z_0_0_0_z_y_x_zz[i] = -g_0_y_x_zz[i] + 2.0 * g_zz_y_x_zz[i] * a_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_z_0_0_0_z_y_y_xx, g_z_0_0_0_z_y_y_xy, g_z_0_0_0_z_y_y_xz, g_z_0_0_0_z_y_y_yy, g_z_0_0_0_z_y_y_yz, g_z_0_0_0_z_y_y_zz, g_zz_y_y_xx, g_zz_y_y_xy, g_zz_y_y_xz, g_zz_y_y_yy, g_zz_y_y_yz, g_zz_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_y_y_xx[i] = -g_0_y_y_xx[i] + 2.0 * g_zz_y_y_xx[i] * a_exp;

        g_z_0_0_0_z_y_y_xy[i] = -g_0_y_y_xy[i] + 2.0 * g_zz_y_y_xy[i] * a_exp;

        g_z_0_0_0_z_y_y_xz[i] = -g_0_y_y_xz[i] + 2.0 * g_zz_y_y_xz[i] * a_exp;

        g_z_0_0_0_z_y_y_yy[i] = -g_0_y_y_yy[i] + 2.0 * g_zz_y_y_yy[i] * a_exp;

        g_z_0_0_0_z_y_y_yz[i] = -g_0_y_y_yz[i] + 2.0 * g_zz_y_y_yz[i] * a_exp;

        g_z_0_0_0_z_y_y_zz[i] = -g_0_y_y_zz[i] + 2.0 * g_zz_y_y_zz[i] * a_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_z_0_0_0_z_y_z_xx, g_z_0_0_0_z_y_z_xy, g_z_0_0_0_z_y_z_xz, g_z_0_0_0_z_y_z_yy, g_z_0_0_0_z_y_z_yz, g_z_0_0_0_z_y_z_zz, g_zz_y_z_xx, g_zz_y_z_xy, g_zz_y_z_xz, g_zz_y_z_yy, g_zz_y_z_yz, g_zz_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_y_z_xx[i] = -g_0_y_z_xx[i] + 2.0 * g_zz_y_z_xx[i] * a_exp;

        g_z_0_0_0_z_y_z_xy[i] = -g_0_y_z_xy[i] + 2.0 * g_zz_y_z_xy[i] * a_exp;

        g_z_0_0_0_z_y_z_xz[i] = -g_0_y_z_xz[i] + 2.0 * g_zz_y_z_xz[i] * a_exp;

        g_z_0_0_0_z_y_z_yy[i] = -g_0_y_z_yy[i] + 2.0 * g_zz_y_z_yy[i] * a_exp;

        g_z_0_0_0_z_y_z_yz[i] = -g_0_y_z_yz[i] + 2.0 * g_zz_y_z_yz[i] * a_exp;

        g_z_0_0_0_z_y_z_zz[i] = -g_0_y_z_zz[i] + 2.0 * g_zz_y_z_zz[i] * a_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_z_0_0_0_z_z_x_xx, g_z_0_0_0_z_z_x_xy, g_z_0_0_0_z_z_x_xz, g_z_0_0_0_z_z_x_yy, g_z_0_0_0_z_z_x_yz, g_z_0_0_0_z_z_x_zz, g_zz_z_x_xx, g_zz_z_x_xy, g_zz_z_x_xz, g_zz_z_x_yy, g_zz_z_x_yz, g_zz_z_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_z_x_xx[i] = -g_0_z_x_xx[i] + 2.0 * g_zz_z_x_xx[i] * a_exp;

        g_z_0_0_0_z_z_x_xy[i] = -g_0_z_x_xy[i] + 2.0 * g_zz_z_x_xy[i] * a_exp;

        g_z_0_0_0_z_z_x_xz[i] = -g_0_z_x_xz[i] + 2.0 * g_zz_z_x_xz[i] * a_exp;

        g_z_0_0_0_z_z_x_yy[i] = -g_0_z_x_yy[i] + 2.0 * g_zz_z_x_yy[i] * a_exp;

        g_z_0_0_0_z_z_x_yz[i] = -g_0_z_x_yz[i] + 2.0 * g_zz_z_x_yz[i] * a_exp;

        g_z_0_0_0_z_z_x_zz[i] = -g_0_z_x_zz[i] + 2.0 * g_zz_z_x_zz[i] * a_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_z_0_0_0_z_z_y_xx, g_z_0_0_0_z_z_y_xy, g_z_0_0_0_z_z_y_xz, g_z_0_0_0_z_z_y_yy, g_z_0_0_0_z_z_y_yz, g_z_0_0_0_z_z_y_zz, g_zz_z_y_xx, g_zz_z_y_xy, g_zz_z_y_xz, g_zz_z_y_yy, g_zz_z_y_yz, g_zz_z_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_z_y_xx[i] = -g_0_z_y_xx[i] + 2.0 * g_zz_z_y_xx[i] * a_exp;

        g_z_0_0_0_z_z_y_xy[i] = -g_0_z_y_xy[i] + 2.0 * g_zz_z_y_xy[i] * a_exp;

        g_z_0_0_0_z_z_y_xz[i] = -g_0_z_y_xz[i] + 2.0 * g_zz_z_y_xz[i] * a_exp;

        g_z_0_0_0_z_z_y_yy[i] = -g_0_z_y_yy[i] + 2.0 * g_zz_z_y_yy[i] * a_exp;

        g_z_0_0_0_z_z_y_yz[i] = -g_0_z_y_yz[i] + 2.0 * g_zz_z_y_yz[i] * a_exp;

        g_z_0_0_0_z_z_y_zz[i] = -g_0_z_y_zz[i] + 2.0 * g_zz_z_y_zz[i] * a_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_z_0_0_0_z_z_z_xx, g_z_0_0_0_z_z_z_xy, g_z_0_0_0_z_z_z_xz, g_z_0_0_0_z_z_z_yy, g_z_0_0_0_z_z_z_yz, g_z_0_0_0_z_z_z_zz, g_zz_z_z_xx, g_zz_z_z_xy, g_zz_z_z_xz, g_zz_z_z_yy, g_zz_z_z_yz, g_zz_z_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_z_z_xx[i] = -g_0_z_z_xx[i] + 2.0 * g_zz_z_z_xx[i] * a_exp;

        g_z_0_0_0_z_z_z_xy[i] = -g_0_z_z_xy[i] + 2.0 * g_zz_z_z_xy[i] * a_exp;

        g_z_0_0_0_z_z_z_xz[i] = -g_0_z_z_xz[i] + 2.0 * g_zz_z_z_xz[i] * a_exp;

        g_z_0_0_0_z_z_z_yy[i] = -g_0_z_z_yy[i] + 2.0 * g_zz_z_z_yy[i] * a_exp;

        g_z_0_0_0_z_z_z_yz[i] = -g_0_z_z_yz[i] + 2.0 * g_zz_z_z_yz[i] * a_exp;

        g_z_0_0_0_z_z_z_zz[i] = -g_0_z_z_zz[i] + 2.0 * g_zz_z_z_zz[i] * a_exp;
    }
}

} // t4c_geom namespace

