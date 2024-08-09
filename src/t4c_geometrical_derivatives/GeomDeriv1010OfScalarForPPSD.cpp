#include "GeomDeriv1010OfScalarForPPSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_ppsd_0(CSimdArray<double>& buffer_1010_ppsd,
                     const CSimdArray<double>& buffer_sppd,
                     const CSimdArray<double>& buffer_dppd,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_ppsd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1010_ppsd

    auto g_x_0_x_0_x_x_0_xx = buffer_1010_ppsd[0];

    auto g_x_0_x_0_x_x_0_xy = buffer_1010_ppsd[1];

    auto g_x_0_x_0_x_x_0_xz = buffer_1010_ppsd[2];

    auto g_x_0_x_0_x_x_0_yy = buffer_1010_ppsd[3];

    auto g_x_0_x_0_x_x_0_yz = buffer_1010_ppsd[4];

    auto g_x_0_x_0_x_x_0_zz = buffer_1010_ppsd[5];

    auto g_x_0_x_0_x_y_0_xx = buffer_1010_ppsd[6];

    auto g_x_0_x_0_x_y_0_xy = buffer_1010_ppsd[7];

    auto g_x_0_x_0_x_y_0_xz = buffer_1010_ppsd[8];

    auto g_x_0_x_0_x_y_0_yy = buffer_1010_ppsd[9];

    auto g_x_0_x_0_x_y_0_yz = buffer_1010_ppsd[10];

    auto g_x_0_x_0_x_y_0_zz = buffer_1010_ppsd[11];

    auto g_x_0_x_0_x_z_0_xx = buffer_1010_ppsd[12];

    auto g_x_0_x_0_x_z_0_xy = buffer_1010_ppsd[13];

    auto g_x_0_x_0_x_z_0_xz = buffer_1010_ppsd[14];

    auto g_x_0_x_0_x_z_0_yy = buffer_1010_ppsd[15];

    auto g_x_0_x_0_x_z_0_yz = buffer_1010_ppsd[16];

    auto g_x_0_x_0_x_z_0_zz = buffer_1010_ppsd[17];

    auto g_x_0_x_0_y_x_0_xx = buffer_1010_ppsd[18];

    auto g_x_0_x_0_y_x_0_xy = buffer_1010_ppsd[19];

    auto g_x_0_x_0_y_x_0_xz = buffer_1010_ppsd[20];

    auto g_x_0_x_0_y_x_0_yy = buffer_1010_ppsd[21];

    auto g_x_0_x_0_y_x_0_yz = buffer_1010_ppsd[22];

    auto g_x_0_x_0_y_x_0_zz = buffer_1010_ppsd[23];

    auto g_x_0_x_0_y_y_0_xx = buffer_1010_ppsd[24];

    auto g_x_0_x_0_y_y_0_xy = buffer_1010_ppsd[25];

    auto g_x_0_x_0_y_y_0_xz = buffer_1010_ppsd[26];

    auto g_x_0_x_0_y_y_0_yy = buffer_1010_ppsd[27];

    auto g_x_0_x_0_y_y_0_yz = buffer_1010_ppsd[28];

    auto g_x_0_x_0_y_y_0_zz = buffer_1010_ppsd[29];

    auto g_x_0_x_0_y_z_0_xx = buffer_1010_ppsd[30];

    auto g_x_0_x_0_y_z_0_xy = buffer_1010_ppsd[31];

    auto g_x_0_x_0_y_z_0_xz = buffer_1010_ppsd[32];

    auto g_x_0_x_0_y_z_0_yy = buffer_1010_ppsd[33];

    auto g_x_0_x_0_y_z_0_yz = buffer_1010_ppsd[34];

    auto g_x_0_x_0_y_z_0_zz = buffer_1010_ppsd[35];

    auto g_x_0_x_0_z_x_0_xx = buffer_1010_ppsd[36];

    auto g_x_0_x_0_z_x_0_xy = buffer_1010_ppsd[37];

    auto g_x_0_x_0_z_x_0_xz = buffer_1010_ppsd[38];

    auto g_x_0_x_0_z_x_0_yy = buffer_1010_ppsd[39];

    auto g_x_0_x_0_z_x_0_yz = buffer_1010_ppsd[40];

    auto g_x_0_x_0_z_x_0_zz = buffer_1010_ppsd[41];

    auto g_x_0_x_0_z_y_0_xx = buffer_1010_ppsd[42];

    auto g_x_0_x_0_z_y_0_xy = buffer_1010_ppsd[43];

    auto g_x_0_x_0_z_y_0_xz = buffer_1010_ppsd[44];

    auto g_x_0_x_0_z_y_0_yy = buffer_1010_ppsd[45];

    auto g_x_0_x_0_z_y_0_yz = buffer_1010_ppsd[46];

    auto g_x_0_x_0_z_y_0_zz = buffer_1010_ppsd[47];

    auto g_x_0_x_0_z_z_0_xx = buffer_1010_ppsd[48];

    auto g_x_0_x_0_z_z_0_xy = buffer_1010_ppsd[49];

    auto g_x_0_x_0_z_z_0_xz = buffer_1010_ppsd[50];

    auto g_x_0_x_0_z_z_0_yy = buffer_1010_ppsd[51];

    auto g_x_0_x_0_z_z_0_yz = buffer_1010_ppsd[52];

    auto g_x_0_x_0_z_z_0_zz = buffer_1010_ppsd[53];

    auto g_x_0_y_0_x_x_0_xx = buffer_1010_ppsd[54];

    auto g_x_0_y_0_x_x_0_xy = buffer_1010_ppsd[55];

    auto g_x_0_y_0_x_x_0_xz = buffer_1010_ppsd[56];

    auto g_x_0_y_0_x_x_0_yy = buffer_1010_ppsd[57];

    auto g_x_0_y_0_x_x_0_yz = buffer_1010_ppsd[58];

    auto g_x_0_y_0_x_x_0_zz = buffer_1010_ppsd[59];

    auto g_x_0_y_0_x_y_0_xx = buffer_1010_ppsd[60];

    auto g_x_0_y_0_x_y_0_xy = buffer_1010_ppsd[61];

    auto g_x_0_y_0_x_y_0_xz = buffer_1010_ppsd[62];

    auto g_x_0_y_0_x_y_0_yy = buffer_1010_ppsd[63];

    auto g_x_0_y_0_x_y_0_yz = buffer_1010_ppsd[64];

    auto g_x_0_y_0_x_y_0_zz = buffer_1010_ppsd[65];

    auto g_x_0_y_0_x_z_0_xx = buffer_1010_ppsd[66];

    auto g_x_0_y_0_x_z_0_xy = buffer_1010_ppsd[67];

    auto g_x_0_y_0_x_z_0_xz = buffer_1010_ppsd[68];

    auto g_x_0_y_0_x_z_0_yy = buffer_1010_ppsd[69];

    auto g_x_0_y_0_x_z_0_yz = buffer_1010_ppsd[70];

    auto g_x_0_y_0_x_z_0_zz = buffer_1010_ppsd[71];

    auto g_x_0_y_0_y_x_0_xx = buffer_1010_ppsd[72];

    auto g_x_0_y_0_y_x_0_xy = buffer_1010_ppsd[73];

    auto g_x_0_y_0_y_x_0_xz = buffer_1010_ppsd[74];

    auto g_x_0_y_0_y_x_0_yy = buffer_1010_ppsd[75];

    auto g_x_0_y_0_y_x_0_yz = buffer_1010_ppsd[76];

    auto g_x_0_y_0_y_x_0_zz = buffer_1010_ppsd[77];

    auto g_x_0_y_0_y_y_0_xx = buffer_1010_ppsd[78];

    auto g_x_0_y_0_y_y_0_xy = buffer_1010_ppsd[79];

    auto g_x_0_y_0_y_y_0_xz = buffer_1010_ppsd[80];

    auto g_x_0_y_0_y_y_0_yy = buffer_1010_ppsd[81];

    auto g_x_0_y_0_y_y_0_yz = buffer_1010_ppsd[82];

    auto g_x_0_y_0_y_y_0_zz = buffer_1010_ppsd[83];

    auto g_x_0_y_0_y_z_0_xx = buffer_1010_ppsd[84];

    auto g_x_0_y_0_y_z_0_xy = buffer_1010_ppsd[85];

    auto g_x_0_y_0_y_z_0_xz = buffer_1010_ppsd[86];

    auto g_x_0_y_0_y_z_0_yy = buffer_1010_ppsd[87];

    auto g_x_0_y_0_y_z_0_yz = buffer_1010_ppsd[88];

    auto g_x_0_y_0_y_z_0_zz = buffer_1010_ppsd[89];

    auto g_x_0_y_0_z_x_0_xx = buffer_1010_ppsd[90];

    auto g_x_0_y_0_z_x_0_xy = buffer_1010_ppsd[91];

    auto g_x_0_y_0_z_x_0_xz = buffer_1010_ppsd[92];

    auto g_x_0_y_0_z_x_0_yy = buffer_1010_ppsd[93];

    auto g_x_0_y_0_z_x_0_yz = buffer_1010_ppsd[94];

    auto g_x_0_y_0_z_x_0_zz = buffer_1010_ppsd[95];

    auto g_x_0_y_0_z_y_0_xx = buffer_1010_ppsd[96];

    auto g_x_0_y_0_z_y_0_xy = buffer_1010_ppsd[97];

    auto g_x_0_y_0_z_y_0_xz = buffer_1010_ppsd[98];

    auto g_x_0_y_0_z_y_0_yy = buffer_1010_ppsd[99];

    auto g_x_0_y_0_z_y_0_yz = buffer_1010_ppsd[100];

    auto g_x_0_y_0_z_y_0_zz = buffer_1010_ppsd[101];

    auto g_x_0_y_0_z_z_0_xx = buffer_1010_ppsd[102];

    auto g_x_0_y_0_z_z_0_xy = buffer_1010_ppsd[103];

    auto g_x_0_y_0_z_z_0_xz = buffer_1010_ppsd[104];

    auto g_x_0_y_0_z_z_0_yy = buffer_1010_ppsd[105];

    auto g_x_0_y_0_z_z_0_yz = buffer_1010_ppsd[106];

    auto g_x_0_y_0_z_z_0_zz = buffer_1010_ppsd[107];

    auto g_x_0_z_0_x_x_0_xx = buffer_1010_ppsd[108];

    auto g_x_0_z_0_x_x_0_xy = buffer_1010_ppsd[109];

    auto g_x_0_z_0_x_x_0_xz = buffer_1010_ppsd[110];

    auto g_x_0_z_0_x_x_0_yy = buffer_1010_ppsd[111];

    auto g_x_0_z_0_x_x_0_yz = buffer_1010_ppsd[112];

    auto g_x_0_z_0_x_x_0_zz = buffer_1010_ppsd[113];

    auto g_x_0_z_0_x_y_0_xx = buffer_1010_ppsd[114];

    auto g_x_0_z_0_x_y_0_xy = buffer_1010_ppsd[115];

    auto g_x_0_z_0_x_y_0_xz = buffer_1010_ppsd[116];

    auto g_x_0_z_0_x_y_0_yy = buffer_1010_ppsd[117];

    auto g_x_0_z_0_x_y_0_yz = buffer_1010_ppsd[118];

    auto g_x_0_z_0_x_y_0_zz = buffer_1010_ppsd[119];

    auto g_x_0_z_0_x_z_0_xx = buffer_1010_ppsd[120];

    auto g_x_0_z_0_x_z_0_xy = buffer_1010_ppsd[121];

    auto g_x_0_z_0_x_z_0_xz = buffer_1010_ppsd[122];

    auto g_x_0_z_0_x_z_0_yy = buffer_1010_ppsd[123];

    auto g_x_0_z_0_x_z_0_yz = buffer_1010_ppsd[124];

    auto g_x_0_z_0_x_z_0_zz = buffer_1010_ppsd[125];

    auto g_x_0_z_0_y_x_0_xx = buffer_1010_ppsd[126];

    auto g_x_0_z_0_y_x_0_xy = buffer_1010_ppsd[127];

    auto g_x_0_z_0_y_x_0_xz = buffer_1010_ppsd[128];

    auto g_x_0_z_0_y_x_0_yy = buffer_1010_ppsd[129];

    auto g_x_0_z_0_y_x_0_yz = buffer_1010_ppsd[130];

    auto g_x_0_z_0_y_x_0_zz = buffer_1010_ppsd[131];

    auto g_x_0_z_0_y_y_0_xx = buffer_1010_ppsd[132];

    auto g_x_0_z_0_y_y_0_xy = buffer_1010_ppsd[133];

    auto g_x_0_z_0_y_y_0_xz = buffer_1010_ppsd[134];

    auto g_x_0_z_0_y_y_0_yy = buffer_1010_ppsd[135];

    auto g_x_0_z_0_y_y_0_yz = buffer_1010_ppsd[136];

    auto g_x_0_z_0_y_y_0_zz = buffer_1010_ppsd[137];

    auto g_x_0_z_0_y_z_0_xx = buffer_1010_ppsd[138];

    auto g_x_0_z_0_y_z_0_xy = buffer_1010_ppsd[139];

    auto g_x_0_z_0_y_z_0_xz = buffer_1010_ppsd[140];

    auto g_x_0_z_0_y_z_0_yy = buffer_1010_ppsd[141];

    auto g_x_0_z_0_y_z_0_yz = buffer_1010_ppsd[142];

    auto g_x_0_z_0_y_z_0_zz = buffer_1010_ppsd[143];

    auto g_x_0_z_0_z_x_0_xx = buffer_1010_ppsd[144];

    auto g_x_0_z_0_z_x_0_xy = buffer_1010_ppsd[145];

    auto g_x_0_z_0_z_x_0_xz = buffer_1010_ppsd[146];

    auto g_x_0_z_0_z_x_0_yy = buffer_1010_ppsd[147];

    auto g_x_0_z_0_z_x_0_yz = buffer_1010_ppsd[148];

    auto g_x_0_z_0_z_x_0_zz = buffer_1010_ppsd[149];

    auto g_x_0_z_0_z_y_0_xx = buffer_1010_ppsd[150];

    auto g_x_0_z_0_z_y_0_xy = buffer_1010_ppsd[151];

    auto g_x_0_z_0_z_y_0_xz = buffer_1010_ppsd[152];

    auto g_x_0_z_0_z_y_0_yy = buffer_1010_ppsd[153];

    auto g_x_0_z_0_z_y_0_yz = buffer_1010_ppsd[154];

    auto g_x_0_z_0_z_y_0_zz = buffer_1010_ppsd[155];

    auto g_x_0_z_0_z_z_0_xx = buffer_1010_ppsd[156];

    auto g_x_0_z_0_z_z_0_xy = buffer_1010_ppsd[157];

    auto g_x_0_z_0_z_z_0_xz = buffer_1010_ppsd[158];

    auto g_x_0_z_0_z_z_0_yy = buffer_1010_ppsd[159];

    auto g_x_0_z_0_z_z_0_yz = buffer_1010_ppsd[160];

    auto g_x_0_z_0_z_z_0_zz = buffer_1010_ppsd[161];

    auto g_y_0_x_0_x_x_0_xx = buffer_1010_ppsd[162];

    auto g_y_0_x_0_x_x_0_xy = buffer_1010_ppsd[163];

    auto g_y_0_x_0_x_x_0_xz = buffer_1010_ppsd[164];

    auto g_y_0_x_0_x_x_0_yy = buffer_1010_ppsd[165];

    auto g_y_0_x_0_x_x_0_yz = buffer_1010_ppsd[166];

    auto g_y_0_x_0_x_x_0_zz = buffer_1010_ppsd[167];

    auto g_y_0_x_0_x_y_0_xx = buffer_1010_ppsd[168];

    auto g_y_0_x_0_x_y_0_xy = buffer_1010_ppsd[169];

    auto g_y_0_x_0_x_y_0_xz = buffer_1010_ppsd[170];

    auto g_y_0_x_0_x_y_0_yy = buffer_1010_ppsd[171];

    auto g_y_0_x_0_x_y_0_yz = buffer_1010_ppsd[172];

    auto g_y_0_x_0_x_y_0_zz = buffer_1010_ppsd[173];

    auto g_y_0_x_0_x_z_0_xx = buffer_1010_ppsd[174];

    auto g_y_0_x_0_x_z_0_xy = buffer_1010_ppsd[175];

    auto g_y_0_x_0_x_z_0_xz = buffer_1010_ppsd[176];

    auto g_y_0_x_0_x_z_0_yy = buffer_1010_ppsd[177];

    auto g_y_0_x_0_x_z_0_yz = buffer_1010_ppsd[178];

    auto g_y_0_x_0_x_z_0_zz = buffer_1010_ppsd[179];

    auto g_y_0_x_0_y_x_0_xx = buffer_1010_ppsd[180];

    auto g_y_0_x_0_y_x_0_xy = buffer_1010_ppsd[181];

    auto g_y_0_x_0_y_x_0_xz = buffer_1010_ppsd[182];

    auto g_y_0_x_0_y_x_0_yy = buffer_1010_ppsd[183];

    auto g_y_0_x_0_y_x_0_yz = buffer_1010_ppsd[184];

    auto g_y_0_x_0_y_x_0_zz = buffer_1010_ppsd[185];

    auto g_y_0_x_0_y_y_0_xx = buffer_1010_ppsd[186];

    auto g_y_0_x_0_y_y_0_xy = buffer_1010_ppsd[187];

    auto g_y_0_x_0_y_y_0_xz = buffer_1010_ppsd[188];

    auto g_y_0_x_0_y_y_0_yy = buffer_1010_ppsd[189];

    auto g_y_0_x_0_y_y_0_yz = buffer_1010_ppsd[190];

    auto g_y_0_x_0_y_y_0_zz = buffer_1010_ppsd[191];

    auto g_y_0_x_0_y_z_0_xx = buffer_1010_ppsd[192];

    auto g_y_0_x_0_y_z_0_xy = buffer_1010_ppsd[193];

    auto g_y_0_x_0_y_z_0_xz = buffer_1010_ppsd[194];

    auto g_y_0_x_0_y_z_0_yy = buffer_1010_ppsd[195];

    auto g_y_0_x_0_y_z_0_yz = buffer_1010_ppsd[196];

    auto g_y_0_x_0_y_z_0_zz = buffer_1010_ppsd[197];

    auto g_y_0_x_0_z_x_0_xx = buffer_1010_ppsd[198];

    auto g_y_0_x_0_z_x_0_xy = buffer_1010_ppsd[199];

    auto g_y_0_x_0_z_x_0_xz = buffer_1010_ppsd[200];

    auto g_y_0_x_0_z_x_0_yy = buffer_1010_ppsd[201];

    auto g_y_0_x_0_z_x_0_yz = buffer_1010_ppsd[202];

    auto g_y_0_x_0_z_x_0_zz = buffer_1010_ppsd[203];

    auto g_y_0_x_0_z_y_0_xx = buffer_1010_ppsd[204];

    auto g_y_0_x_0_z_y_0_xy = buffer_1010_ppsd[205];

    auto g_y_0_x_0_z_y_0_xz = buffer_1010_ppsd[206];

    auto g_y_0_x_0_z_y_0_yy = buffer_1010_ppsd[207];

    auto g_y_0_x_0_z_y_0_yz = buffer_1010_ppsd[208];

    auto g_y_0_x_0_z_y_0_zz = buffer_1010_ppsd[209];

    auto g_y_0_x_0_z_z_0_xx = buffer_1010_ppsd[210];

    auto g_y_0_x_0_z_z_0_xy = buffer_1010_ppsd[211];

    auto g_y_0_x_0_z_z_0_xz = buffer_1010_ppsd[212];

    auto g_y_0_x_0_z_z_0_yy = buffer_1010_ppsd[213];

    auto g_y_0_x_0_z_z_0_yz = buffer_1010_ppsd[214];

    auto g_y_0_x_0_z_z_0_zz = buffer_1010_ppsd[215];

    auto g_y_0_y_0_x_x_0_xx = buffer_1010_ppsd[216];

    auto g_y_0_y_0_x_x_0_xy = buffer_1010_ppsd[217];

    auto g_y_0_y_0_x_x_0_xz = buffer_1010_ppsd[218];

    auto g_y_0_y_0_x_x_0_yy = buffer_1010_ppsd[219];

    auto g_y_0_y_0_x_x_0_yz = buffer_1010_ppsd[220];

    auto g_y_0_y_0_x_x_0_zz = buffer_1010_ppsd[221];

    auto g_y_0_y_0_x_y_0_xx = buffer_1010_ppsd[222];

    auto g_y_0_y_0_x_y_0_xy = buffer_1010_ppsd[223];

    auto g_y_0_y_0_x_y_0_xz = buffer_1010_ppsd[224];

    auto g_y_0_y_0_x_y_0_yy = buffer_1010_ppsd[225];

    auto g_y_0_y_0_x_y_0_yz = buffer_1010_ppsd[226];

    auto g_y_0_y_0_x_y_0_zz = buffer_1010_ppsd[227];

    auto g_y_0_y_0_x_z_0_xx = buffer_1010_ppsd[228];

    auto g_y_0_y_0_x_z_0_xy = buffer_1010_ppsd[229];

    auto g_y_0_y_0_x_z_0_xz = buffer_1010_ppsd[230];

    auto g_y_0_y_0_x_z_0_yy = buffer_1010_ppsd[231];

    auto g_y_0_y_0_x_z_0_yz = buffer_1010_ppsd[232];

    auto g_y_0_y_0_x_z_0_zz = buffer_1010_ppsd[233];

    auto g_y_0_y_0_y_x_0_xx = buffer_1010_ppsd[234];

    auto g_y_0_y_0_y_x_0_xy = buffer_1010_ppsd[235];

    auto g_y_0_y_0_y_x_0_xz = buffer_1010_ppsd[236];

    auto g_y_0_y_0_y_x_0_yy = buffer_1010_ppsd[237];

    auto g_y_0_y_0_y_x_0_yz = buffer_1010_ppsd[238];

    auto g_y_0_y_0_y_x_0_zz = buffer_1010_ppsd[239];

    auto g_y_0_y_0_y_y_0_xx = buffer_1010_ppsd[240];

    auto g_y_0_y_0_y_y_0_xy = buffer_1010_ppsd[241];

    auto g_y_0_y_0_y_y_0_xz = buffer_1010_ppsd[242];

    auto g_y_0_y_0_y_y_0_yy = buffer_1010_ppsd[243];

    auto g_y_0_y_0_y_y_0_yz = buffer_1010_ppsd[244];

    auto g_y_0_y_0_y_y_0_zz = buffer_1010_ppsd[245];

    auto g_y_0_y_0_y_z_0_xx = buffer_1010_ppsd[246];

    auto g_y_0_y_0_y_z_0_xy = buffer_1010_ppsd[247];

    auto g_y_0_y_0_y_z_0_xz = buffer_1010_ppsd[248];

    auto g_y_0_y_0_y_z_0_yy = buffer_1010_ppsd[249];

    auto g_y_0_y_0_y_z_0_yz = buffer_1010_ppsd[250];

    auto g_y_0_y_0_y_z_0_zz = buffer_1010_ppsd[251];

    auto g_y_0_y_0_z_x_0_xx = buffer_1010_ppsd[252];

    auto g_y_0_y_0_z_x_0_xy = buffer_1010_ppsd[253];

    auto g_y_0_y_0_z_x_0_xz = buffer_1010_ppsd[254];

    auto g_y_0_y_0_z_x_0_yy = buffer_1010_ppsd[255];

    auto g_y_0_y_0_z_x_0_yz = buffer_1010_ppsd[256];

    auto g_y_0_y_0_z_x_0_zz = buffer_1010_ppsd[257];

    auto g_y_0_y_0_z_y_0_xx = buffer_1010_ppsd[258];

    auto g_y_0_y_0_z_y_0_xy = buffer_1010_ppsd[259];

    auto g_y_0_y_0_z_y_0_xz = buffer_1010_ppsd[260];

    auto g_y_0_y_0_z_y_0_yy = buffer_1010_ppsd[261];

    auto g_y_0_y_0_z_y_0_yz = buffer_1010_ppsd[262];

    auto g_y_0_y_0_z_y_0_zz = buffer_1010_ppsd[263];

    auto g_y_0_y_0_z_z_0_xx = buffer_1010_ppsd[264];

    auto g_y_0_y_0_z_z_0_xy = buffer_1010_ppsd[265];

    auto g_y_0_y_0_z_z_0_xz = buffer_1010_ppsd[266];

    auto g_y_0_y_0_z_z_0_yy = buffer_1010_ppsd[267];

    auto g_y_0_y_0_z_z_0_yz = buffer_1010_ppsd[268];

    auto g_y_0_y_0_z_z_0_zz = buffer_1010_ppsd[269];

    auto g_y_0_z_0_x_x_0_xx = buffer_1010_ppsd[270];

    auto g_y_0_z_0_x_x_0_xy = buffer_1010_ppsd[271];

    auto g_y_0_z_0_x_x_0_xz = buffer_1010_ppsd[272];

    auto g_y_0_z_0_x_x_0_yy = buffer_1010_ppsd[273];

    auto g_y_0_z_0_x_x_0_yz = buffer_1010_ppsd[274];

    auto g_y_0_z_0_x_x_0_zz = buffer_1010_ppsd[275];

    auto g_y_0_z_0_x_y_0_xx = buffer_1010_ppsd[276];

    auto g_y_0_z_0_x_y_0_xy = buffer_1010_ppsd[277];

    auto g_y_0_z_0_x_y_0_xz = buffer_1010_ppsd[278];

    auto g_y_0_z_0_x_y_0_yy = buffer_1010_ppsd[279];

    auto g_y_0_z_0_x_y_0_yz = buffer_1010_ppsd[280];

    auto g_y_0_z_0_x_y_0_zz = buffer_1010_ppsd[281];

    auto g_y_0_z_0_x_z_0_xx = buffer_1010_ppsd[282];

    auto g_y_0_z_0_x_z_0_xy = buffer_1010_ppsd[283];

    auto g_y_0_z_0_x_z_0_xz = buffer_1010_ppsd[284];

    auto g_y_0_z_0_x_z_0_yy = buffer_1010_ppsd[285];

    auto g_y_0_z_0_x_z_0_yz = buffer_1010_ppsd[286];

    auto g_y_0_z_0_x_z_0_zz = buffer_1010_ppsd[287];

    auto g_y_0_z_0_y_x_0_xx = buffer_1010_ppsd[288];

    auto g_y_0_z_0_y_x_0_xy = buffer_1010_ppsd[289];

    auto g_y_0_z_0_y_x_0_xz = buffer_1010_ppsd[290];

    auto g_y_0_z_0_y_x_0_yy = buffer_1010_ppsd[291];

    auto g_y_0_z_0_y_x_0_yz = buffer_1010_ppsd[292];

    auto g_y_0_z_0_y_x_0_zz = buffer_1010_ppsd[293];

    auto g_y_0_z_0_y_y_0_xx = buffer_1010_ppsd[294];

    auto g_y_0_z_0_y_y_0_xy = buffer_1010_ppsd[295];

    auto g_y_0_z_0_y_y_0_xz = buffer_1010_ppsd[296];

    auto g_y_0_z_0_y_y_0_yy = buffer_1010_ppsd[297];

    auto g_y_0_z_0_y_y_0_yz = buffer_1010_ppsd[298];

    auto g_y_0_z_0_y_y_0_zz = buffer_1010_ppsd[299];

    auto g_y_0_z_0_y_z_0_xx = buffer_1010_ppsd[300];

    auto g_y_0_z_0_y_z_0_xy = buffer_1010_ppsd[301];

    auto g_y_0_z_0_y_z_0_xz = buffer_1010_ppsd[302];

    auto g_y_0_z_0_y_z_0_yy = buffer_1010_ppsd[303];

    auto g_y_0_z_0_y_z_0_yz = buffer_1010_ppsd[304];

    auto g_y_0_z_0_y_z_0_zz = buffer_1010_ppsd[305];

    auto g_y_0_z_0_z_x_0_xx = buffer_1010_ppsd[306];

    auto g_y_0_z_0_z_x_0_xy = buffer_1010_ppsd[307];

    auto g_y_0_z_0_z_x_0_xz = buffer_1010_ppsd[308];

    auto g_y_0_z_0_z_x_0_yy = buffer_1010_ppsd[309];

    auto g_y_0_z_0_z_x_0_yz = buffer_1010_ppsd[310];

    auto g_y_0_z_0_z_x_0_zz = buffer_1010_ppsd[311];

    auto g_y_0_z_0_z_y_0_xx = buffer_1010_ppsd[312];

    auto g_y_0_z_0_z_y_0_xy = buffer_1010_ppsd[313];

    auto g_y_0_z_0_z_y_0_xz = buffer_1010_ppsd[314];

    auto g_y_0_z_0_z_y_0_yy = buffer_1010_ppsd[315];

    auto g_y_0_z_0_z_y_0_yz = buffer_1010_ppsd[316];

    auto g_y_0_z_0_z_y_0_zz = buffer_1010_ppsd[317];

    auto g_y_0_z_0_z_z_0_xx = buffer_1010_ppsd[318];

    auto g_y_0_z_0_z_z_0_xy = buffer_1010_ppsd[319];

    auto g_y_0_z_0_z_z_0_xz = buffer_1010_ppsd[320];

    auto g_y_0_z_0_z_z_0_yy = buffer_1010_ppsd[321];

    auto g_y_0_z_0_z_z_0_yz = buffer_1010_ppsd[322];

    auto g_y_0_z_0_z_z_0_zz = buffer_1010_ppsd[323];

    auto g_z_0_x_0_x_x_0_xx = buffer_1010_ppsd[324];

    auto g_z_0_x_0_x_x_0_xy = buffer_1010_ppsd[325];

    auto g_z_0_x_0_x_x_0_xz = buffer_1010_ppsd[326];

    auto g_z_0_x_0_x_x_0_yy = buffer_1010_ppsd[327];

    auto g_z_0_x_0_x_x_0_yz = buffer_1010_ppsd[328];

    auto g_z_0_x_0_x_x_0_zz = buffer_1010_ppsd[329];

    auto g_z_0_x_0_x_y_0_xx = buffer_1010_ppsd[330];

    auto g_z_0_x_0_x_y_0_xy = buffer_1010_ppsd[331];

    auto g_z_0_x_0_x_y_0_xz = buffer_1010_ppsd[332];

    auto g_z_0_x_0_x_y_0_yy = buffer_1010_ppsd[333];

    auto g_z_0_x_0_x_y_0_yz = buffer_1010_ppsd[334];

    auto g_z_0_x_0_x_y_0_zz = buffer_1010_ppsd[335];

    auto g_z_0_x_0_x_z_0_xx = buffer_1010_ppsd[336];

    auto g_z_0_x_0_x_z_0_xy = buffer_1010_ppsd[337];

    auto g_z_0_x_0_x_z_0_xz = buffer_1010_ppsd[338];

    auto g_z_0_x_0_x_z_0_yy = buffer_1010_ppsd[339];

    auto g_z_0_x_0_x_z_0_yz = buffer_1010_ppsd[340];

    auto g_z_0_x_0_x_z_0_zz = buffer_1010_ppsd[341];

    auto g_z_0_x_0_y_x_0_xx = buffer_1010_ppsd[342];

    auto g_z_0_x_0_y_x_0_xy = buffer_1010_ppsd[343];

    auto g_z_0_x_0_y_x_0_xz = buffer_1010_ppsd[344];

    auto g_z_0_x_0_y_x_0_yy = buffer_1010_ppsd[345];

    auto g_z_0_x_0_y_x_0_yz = buffer_1010_ppsd[346];

    auto g_z_0_x_0_y_x_0_zz = buffer_1010_ppsd[347];

    auto g_z_0_x_0_y_y_0_xx = buffer_1010_ppsd[348];

    auto g_z_0_x_0_y_y_0_xy = buffer_1010_ppsd[349];

    auto g_z_0_x_0_y_y_0_xz = buffer_1010_ppsd[350];

    auto g_z_0_x_0_y_y_0_yy = buffer_1010_ppsd[351];

    auto g_z_0_x_0_y_y_0_yz = buffer_1010_ppsd[352];

    auto g_z_0_x_0_y_y_0_zz = buffer_1010_ppsd[353];

    auto g_z_0_x_0_y_z_0_xx = buffer_1010_ppsd[354];

    auto g_z_0_x_0_y_z_0_xy = buffer_1010_ppsd[355];

    auto g_z_0_x_0_y_z_0_xz = buffer_1010_ppsd[356];

    auto g_z_0_x_0_y_z_0_yy = buffer_1010_ppsd[357];

    auto g_z_0_x_0_y_z_0_yz = buffer_1010_ppsd[358];

    auto g_z_0_x_0_y_z_0_zz = buffer_1010_ppsd[359];

    auto g_z_0_x_0_z_x_0_xx = buffer_1010_ppsd[360];

    auto g_z_0_x_0_z_x_0_xy = buffer_1010_ppsd[361];

    auto g_z_0_x_0_z_x_0_xz = buffer_1010_ppsd[362];

    auto g_z_0_x_0_z_x_0_yy = buffer_1010_ppsd[363];

    auto g_z_0_x_0_z_x_0_yz = buffer_1010_ppsd[364];

    auto g_z_0_x_0_z_x_0_zz = buffer_1010_ppsd[365];

    auto g_z_0_x_0_z_y_0_xx = buffer_1010_ppsd[366];

    auto g_z_0_x_0_z_y_0_xy = buffer_1010_ppsd[367];

    auto g_z_0_x_0_z_y_0_xz = buffer_1010_ppsd[368];

    auto g_z_0_x_0_z_y_0_yy = buffer_1010_ppsd[369];

    auto g_z_0_x_0_z_y_0_yz = buffer_1010_ppsd[370];

    auto g_z_0_x_0_z_y_0_zz = buffer_1010_ppsd[371];

    auto g_z_0_x_0_z_z_0_xx = buffer_1010_ppsd[372];

    auto g_z_0_x_0_z_z_0_xy = buffer_1010_ppsd[373];

    auto g_z_0_x_0_z_z_0_xz = buffer_1010_ppsd[374];

    auto g_z_0_x_0_z_z_0_yy = buffer_1010_ppsd[375];

    auto g_z_0_x_0_z_z_0_yz = buffer_1010_ppsd[376];

    auto g_z_0_x_0_z_z_0_zz = buffer_1010_ppsd[377];

    auto g_z_0_y_0_x_x_0_xx = buffer_1010_ppsd[378];

    auto g_z_0_y_0_x_x_0_xy = buffer_1010_ppsd[379];

    auto g_z_0_y_0_x_x_0_xz = buffer_1010_ppsd[380];

    auto g_z_0_y_0_x_x_0_yy = buffer_1010_ppsd[381];

    auto g_z_0_y_0_x_x_0_yz = buffer_1010_ppsd[382];

    auto g_z_0_y_0_x_x_0_zz = buffer_1010_ppsd[383];

    auto g_z_0_y_0_x_y_0_xx = buffer_1010_ppsd[384];

    auto g_z_0_y_0_x_y_0_xy = buffer_1010_ppsd[385];

    auto g_z_0_y_0_x_y_0_xz = buffer_1010_ppsd[386];

    auto g_z_0_y_0_x_y_0_yy = buffer_1010_ppsd[387];

    auto g_z_0_y_0_x_y_0_yz = buffer_1010_ppsd[388];

    auto g_z_0_y_0_x_y_0_zz = buffer_1010_ppsd[389];

    auto g_z_0_y_0_x_z_0_xx = buffer_1010_ppsd[390];

    auto g_z_0_y_0_x_z_0_xy = buffer_1010_ppsd[391];

    auto g_z_0_y_0_x_z_0_xz = buffer_1010_ppsd[392];

    auto g_z_0_y_0_x_z_0_yy = buffer_1010_ppsd[393];

    auto g_z_0_y_0_x_z_0_yz = buffer_1010_ppsd[394];

    auto g_z_0_y_0_x_z_0_zz = buffer_1010_ppsd[395];

    auto g_z_0_y_0_y_x_0_xx = buffer_1010_ppsd[396];

    auto g_z_0_y_0_y_x_0_xy = buffer_1010_ppsd[397];

    auto g_z_0_y_0_y_x_0_xz = buffer_1010_ppsd[398];

    auto g_z_0_y_0_y_x_0_yy = buffer_1010_ppsd[399];

    auto g_z_0_y_0_y_x_0_yz = buffer_1010_ppsd[400];

    auto g_z_0_y_0_y_x_0_zz = buffer_1010_ppsd[401];

    auto g_z_0_y_0_y_y_0_xx = buffer_1010_ppsd[402];

    auto g_z_0_y_0_y_y_0_xy = buffer_1010_ppsd[403];

    auto g_z_0_y_0_y_y_0_xz = buffer_1010_ppsd[404];

    auto g_z_0_y_0_y_y_0_yy = buffer_1010_ppsd[405];

    auto g_z_0_y_0_y_y_0_yz = buffer_1010_ppsd[406];

    auto g_z_0_y_0_y_y_0_zz = buffer_1010_ppsd[407];

    auto g_z_0_y_0_y_z_0_xx = buffer_1010_ppsd[408];

    auto g_z_0_y_0_y_z_0_xy = buffer_1010_ppsd[409];

    auto g_z_0_y_0_y_z_0_xz = buffer_1010_ppsd[410];

    auto g_z_0_y_0_y_z_0_yy = buffer_1010_ppsd[411];

    auto g_z_0_y_0_y_z_0_yz = buffer_1010_ppsd[412];

    auto g_z_0_y_0_y_z_0_zz = buffer_1010_ppsd[413];

    auto g_z_0_y_0_z_x_0_xx = buffer_1010_ppsd[414];

    auto g_z_0_y_0_z_x_0_xy = buffer_1010_ppsd[415];

    auto g_z_0_y_0_z_x_0_xz = buffer_1010_ppsd[416];

    auto g_z_0_y_0_z_x_0_yy = buffer_1010_ppsd[417];

    auto g_z_0_y_0_z_x_0_yz = buffer_1010_ppsd[418];

    auto g_z_0_y_0_z_x_0_zz = buffer_1010_ppsd[419];

    auto g_z_0_y_0_z_y_0_xx = buffer_1010_ppsd[420];

    auto g_z_0_y_0_z_y_0_xy = buffer_1010_ppsd[421];

    auto g_z_0_y_0_z_y_0_xz = buffer_1010_ppsd[422];

    auto g_z_0_y_0_z_y_0_yy = buffer_1010_ppsd[423];

    auto g_z_0_y_0_z_y_0_yz = buffer_1010_ppsd[424];

    auto g_z_0_y_0_z_y_0_zz = buffer_1010_ppsd[425];

    auto g_z_0_y_0_z_z_0_xx = buffer_1010_ppsd[426];

    auto g_z_0_y_0_z_z_0_xy = buffer_1010_ppsd[427];

    auto g_z_0_y_0_z_z_0_xz = buffer_1010_ppsd[428];

    auto g_z_0_y_0_z_z_0_yy = buffer_1010_ppsd[429];

    auto g_z_0_y_0_z_z_0_yz = buffer_1010_ppsd[430];

    auto g_z_0_y_0_z_z_0_zz = buffer_1010_ppsd[431];

    auto g_z_0_z_0_x_x_0_xx = buffer_1010_ppsd[432];

    auto g_z_0_z_0_x_x_0_xy = buffer_1010_ppsd[433];

    auto g_z_0_z_0_x_x_0_xz = buffer_1010_ppsd[434];

    auto g_z_0_z_0_x_x_0_yy = buffer_1010_ppsd[435];

    auto g_z_0_z_0_x_x_0_yz = buffer_1010_ppsd[436];

    auto g_z_0_z_0_x_x_0_zz = buffer_1010_ppsd[437];

    auto g_z_0_z_0_x_y_0_xx = buffer_1010_ppsd[438];

    auto g_z_0_z_0_x_y_0_xy = buffer_1010_ppsd[439];

    auto g_z_0_z_0_x_y_0_xz = buffer_1010_ppsd[440];

    auto g_z_0_z_0_x_y_0_yy = buffer_1010_ppsd[441];

    auto g_z_0_z_0_x_y_0_yz = buffer_1010_ppsd[442];

    auto g_z_0_z_0_x_y_0_zz = buffer_1010_ppsd[443];

    auto g_z_0_z_0_x_z_0_xx = buffer_1010_ppsd[444];

    auto g_z_0_z_0_x_z_0_xy = buffer_1010_ppsd[445];

    auto g_z_0_z_0_x_z_0_xz = buffer_1010_ppsd[446];

    auto g_z_0_z_0_x_z_0_yy = buffer_1010_ppsd[447];

    auto g_z_0_z_0_x_z_0_yz = buffer_1010_ppsd[448];

    auto g_z_0_z_0_x_z_0_zz = buffer_1010_ppsd[449];

    auto g_z_0_z_0_y_x_0_xx = buffer_1010_ppsd[450];

    auto g_z_0_z_0_y_x_0_xy = buffer_1010_ppsd[451];

    auto g_z_0_z_0_y_x_0_xz = buffer_1010_ppsd[452];

    auto g_z_0_z_0_y_x_0_yy = buffer_1010_ppsd[453];

    auto g_z_0_z_0_y_x_0_yz = buffer_1010_ppsd[454];

    auto g_z_0_z_0_y_x_0_zz = buffer_1010_ppsd[455];

    auto g_z_0_z_0_y_y_0_xx = buffer_1010_ppsd[456];

    auto g_z_0_z_0_y_y_0_xy = buffer_1010_ppsd[457];

    auto g_z_0_z_0_y_y_0_xz = buffer_1010_ppsd[458];

    auto g_z_0_z_0_y_y_0_yy = buffer_1010_ppsd[459];

    auto g_z_0_z_0_y_y_0_yz = buffer_1010_ppsd[460];

    auto g_z_0_z_0_y_y_0_zz = buffer_1010_ppsd[461];

    auto g_z_0_z_0_y_z_0_xx = buffer_1010_ppsd[462];

    auto g_z_0_z_0_y_z_0_xy = buffer_1010_ppsd[463];

    auto g_z_0_z_0_y_z_0_xz = buffer_1010_ppsd[464];

    auto g_z_0_z_0_y_z_0_yy = buffer_1010_ppsd[465];

    auto g_z_0_z_0_y_z_0_yz = buffer_1010_ppsd[466];

    auto g_z_0_z_0_y_z_0_zz = buffer_1010_ppsd[467];

    auto g_z_0_z_0_z_x_0_xx = buffer_1010_ppsd[468];

    auto g_z_0_z_0_z_x_0_xy = buffer_1010_ppsd[469];

    auto g_z_0_z_0_z_x_0_xz = buffer_1010_ppsd[470];

    auto g_z_0_z_0_z_x_0_yy = buffer_1010_ppsd[471];

    auto g_z_0_z_0_z_x_0_yz = buffer_1010_ppsd[472];

    auto g_z_0_z_0_z_x_0_zz = buffer_1010_ppsd[473];

    auto g_z_0_z_0_z_y_0_xx = buffer_1010_ppsd[474];

    auto g_z_0_z_0_z_y_0_xy = buffer_1010_ppsd[475];

    auto g_z_0_z_0_z_y_0_xz = buffer_1010_ppsd[476];

    auto g_z_0_z_0_z_y_0_yy = buffer_1010_ppsd[477];

    auto g_z_0_z_0_z_y_0_yz = buffer_1010_ppsd[478];

    auto g_z_0_z_0_z_y_0_zz = buffer_1010_ppsd[479];

    auto g_z_0_z_0_z_z_0_xx = buffer_1010_ppsd[480];

    auto g_z_0_z_0_z_z_0_xy = buffer_1010_ppsd[481];

    auto g_z_0_z_0_z_z_0_xz = buffer_1010_ppsd[482];

    auto g_z_0_z_0_z_z_0_yy = buffer_1010_ppsd[483];

    auto g_z_0_z_0_z_z_0_yz = buffer_1010_ppsd[484];

    auto g_z_0_z_0_z_z_0_zz = buffer_1010_ppsd[485];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_x_0_x_0_x_x_0_xx, g_x_0_x_0_x_x_0_xy, g_x_0_x_0_x_x_0_xz, g_x_0_x_0_x_x_0_yy, g_x_0_x_0_x_x_0_yz, g_x_0_x_0_x_x_0_zz, g_xx_x_x_xx, g_xx_x_x_xy, g_xx_x_x_xz, g_xx_x_x_yy, g_xx_x_x_yz, g_xx_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_x_0_xx[i] = -2.0 * g_0_x_x_xx[i] * c_exps[i] + 4.0 * g_xx_x_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_0_xy[i] = -2.0 * g_0_x_x_xy[i] * c_exps[i] + 4.0 * g_xx_x_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_0_xz[i] = -2.0 * g_0_x_x_xz[i] * c_exps[i] + 4.0 * g_xx_x_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_0_yy[i] = -2.0 * g_0_x_x_yy[i] * c_exps[i] + 4.0 * g_xx_x_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_0_yz[i] = -2.0 * g_0_x_x_yz[i] * c_exps[i] + 4.0 * g_xx_x_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_x_0_zz[i] = -2.0 * g_0_x_x_zz[i] * c_exps[i] + 4.0 * g_xx_x_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_x_0_x_0_x_y_0_xx, g_x_0_x_0_x_y_0_xy, g_x_0_x_0_x_y_0_xz, g_x_0_x_0_x_y_0_yy, g_x_0_x_0_x_y_0_yz, g_x_0_x_0_x_y_0_zz, g_xx_y_x_xx, g_xx_y_x_xy, g_xx_y_x_xz, g_xx_y_x_yy, g_xx_y_x_yz, g_xx_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_y_0_xx[i] = -2.0 * g_0_y_x_xx[i] * c_exps[i] + 4.0 * g_xx_y_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_0_xy[i] = -2.0 * g_0_y_x_xy[i] * c_exps[i] + 4.0 * g_xx_y_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_0_xz[i] = -2.0 * g_0_y_x_xz[i] * c_exps[i] + 4.0 * g_xx_y_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_0_yy[i] = -2.0 * g_0_y_x_yy[i] * c_exps[i] + 4.0 * g_xx_y_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_0_yz[i] = -2.0 * g_0_y_x_yz[i] * c_exps[i] + 4.0 * g_xx_y_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_y_0_zz[i] = -2.0 * g_0_y_x_zz[i] * c_exps[i] + 4.0 * g_xx_y_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_x_0_x_0_x_z_0_xx, g_x_0_x_0_x_z_0_xy, g_x_0_x_0_x_z_0_xz, g_x_0_x_0_x_z_0_yy, g_x_0_x_0_x_z_0_yz, g_x_0_x_0_x_z_0_zz, g_xx_z_x_xx, g_xx_z_x_xy, g_xx_z_x_xz, g_xx_z_x_yy, g_xx_z_x_yz, g_xx_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_x_z_0_xx[i] = -2.0 * g_0_z_x_xx[i] * c_exps[i] + 4.0 * g_xx_z_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_0_xy[i] = -2.0 * g_0_z_x_xy[i] * c_exps[i] + 4.0 * g_xx_z_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_0_xz[i] = -2.0 * g_0_z_x_xz[i] * c_exps[i] + 4.0 * g_xx_z_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_0_yy[i] = -2.0 * g_0_z_x_yy[i] * c_exps[i] + 4.0 * g_xx_z_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_0_yz[i] = -2.0 * g_0_z_x_yz[i] * c_exps[i] + 4.0 * g_xx_z_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_x_z_0_zz[i] = -2.0 * g_0_z_x_zz[i] * c_exps[i] + 4.0 * g_xx_z_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_x_0_y_x_0_xx, g_x_0_x_0_y_x_0_xy, g_x_0_x_0_y_x_0_xz, g_x_0_x_0_y_x_0_yy, g_x_0_x_0_y_x_0_yz, g_x_0_x_0_y_x_0_zz, g_xy_x_x_xx, g_xy_x_x_xy, g_xy_x_x_xz, g_xy_x_x_yy, g_xy_x_x_yz, g_xy_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_x_0_xx[i] = 4.0 * g_xy_x_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_0_xy[i] = 4.0 * g_xy_x_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_0_xz[i] = 4.0 * g_xy_x_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_0_yy[i] = 4.0 * g_xy_x_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_0_yz[i] = 4.0 * g_xy_x_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_x_0_zz[i] = 4.0 * g_xy_x_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_x_0_y_y_0_xx, g_x_0_x_0_y_y_0_xy, g_x_0_x_0_y_y_0_xz, g_x_0_x_0_y_y_0_yy, g_x_0_x_0_y_y_0_yz, g_x_0_x_0_y_y_0_zz, g_xy_y_x_xx, g_xy_y_x_xy, g_xy_y_x_xz, g_xy_y_x_yy, g_xy_y_x_yz, g_xy_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_y_0_xx[i] = 4.0 * g_xy_y_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_0_xy[i] = 4.0 * g_xy_y_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_0_xz[i] = 4.0 * g_xy_y_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_0_yy[i] = 4.0 * g_xy_y_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_0_yz[i] = 4.0 * g_xy_y_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_y_0_zz[i] = 4.0 * g_xy_y_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_x_0_y_z_0_xx, g_x_0_x_0_y_z_0_xy, g_x_0_x_0_y_z_0_xz, g_x_0_x_0_y_z_0_yy, g_x_0_x_0_y_z_0_yz, g_x_0_x_0_y_z_0_zz, g_xy_z_x_xx, g_xy_z_x_xy, g_xy_z_x_xz, g_xy_z_x_yy, g_xy_z_x_yz, g_xy_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_y_z_0_xx[i] = 4.0 * g_xy_z_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_0_xy[i] = 4.0 * g_xy_z_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_0_xz[i] = 4.0 * g_xy_z_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_0_yy[i] = 4.0 * g_xy_z_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_0_yz[i] = 4.0 * g_xy_z_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_y_z_0_zz[i] = 4.0 * g_xy_z_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_x_0_z_x_0_xx, g_x_0_x_0_z_x_0_xy, g_x_0_x_0_z_x_0_xz, g_x_0_x_0_z_x_0_yy, g_x_0_x_0_z_x_0_yz, g_x_0_x_0_z_x_0_zz, g_xz_x_x_xx, g_xz_x_x_xy, g_xz_x_x_xz, g_xz_x_x_yy, g_xz_x_x_yz, g_xz_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_x_0_xx[i] = 4.0 * g_xz_x_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_0_xy[i] = 4.0 * g_xz_x_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_0_xz[i] = 4.0 * g_xz_x_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_0_yy[i] = 4.0 * g_xz_x_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_0_yz[i] = 4.0 * g_xz_x_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_x_0_zz[i] = 4.0 * g_xz_x_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_x_0_z_y_0_xx, g_x_0_x_0_z_y_0_xy, g_x_0_x_0_z_y_0_xz, g_x_0_x_0_z_y_0_yy, g_x_0_x_0_z_y_0_yz, g_x_0_x_0_z_y_0_zz, g_xz_y_x_xx, g_xz_y_x_xy, g_xz_y_x_xz, g_xz_y_x_yy, g_xz_y_x_yz, g_xz_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_y_0_xx[i] = 4.0 * g_xz_y_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_0_xy[i] = 4.0 * g_xz_y_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_0_xz[i] = 4.0 * g_xz_y_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_0_yy[i] = 4.0 * g_xz_y_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_0_yz[i] = 4.0 * g_xz_y_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_y_0_zz[i] = 4.0 * g_xz_y_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_x_0_z_z_0_xx, g_x_0_x_0_z_z_0_xy, g_x_0_x_0_z_z_0_xz, g_x_0_x_0_z_z_0_yy, g_x_0_x_0_z_z_0_yz, g_x_0_x_0_z_z_0_zz, g_xz_z_x_xx, g_xz_z_x_xy, g_xz_z_x_xz, g_xz_z_x_yy, g_xz_z_x_yz, g_xz_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_z_z_0_xx[i] = 4.0 * g_xz_z_x_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_0_xy[i] = 4.0 * g_xz_z_x_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_0_xz[i] = 4.0 * g_xz_z_x_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_0_yy[i] = 4.0 * g_xz_z_x_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_0_yz[i] = 4.0 * g_xz_z_x_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_z_z_0_zz[i] = 4.0 * g_xz_z_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_x_0_y_0_x_x_0_xx, g_x_0_y_0_x_x_0_xy, g_x_0_y_0_x_x_0_xz, g_x_0_y_0_x_x_0_yy, g_x_0_y_0_x_x_0_yz, g_x_0_y_0_x_x_0_zz, g_xx_x_y_xx, g_xx_x_y_xy, g_xx_x_y_xz, g_xx_x_y_yy, g_xx_x_y_yz, g_xx_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_x_0_xx[i] = -2.0 * g_0_x_y_xx[i] * c_exps[i] + 4.0 * g_xx_x_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_0_xy[i] = -2.0 * g_0_x_y_xy[i] * c_exps[i] + 4.0 * g_xx_x_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_0_xz[i] = -2.0 * g_0_x_y_xz[i] * c_exps[i] + 4.0 * g_xx_x_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_0_yy[i] = -2.0 * g_0_x_y_yy[i] * c_exps[i] + 4.0 * g_xx_x_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_0_yz[i] = -2.0 * g_0_x_y_yz[i] * c_exps[i] + 4.0 * g_xx_x_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_x_0_zz[i] = -2.0 * g_0_x_y_zz[i] * c_exps[i] + 4.0 * g_xx_x_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_x_0_y_0_x_y_0_xx, g_x_0_y_0_x_y_0_xy, g_x_0_y_0_x_y_0_xz, g_x_0_y_0_x_y_0_yy, g_x_0_y_0_x_y_0_yz, g_x_0_y_0_x_y_0_zz, g_xx_y_y_xx, g_xx_y_y_xy, g_xx_y_y_xz, g_xx_y_y_yy, g_xx_y_y_yz, g_xx_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_y_0_xx[i] = -2.0 * g_0_y_y_xx[i] * c_exps[i] + 4.0 * g_xx_y_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_0_xy[i] = -2.0 * g_0_y_y_xy[i] * c_exps[i] + 4.0 * g_xx_y_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_0_xz[i] = -2.0 * g_0_y_y_xz[i] * c_exps[i] + 4.0 * g_xx_y_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_0_yy[i] = -2.0 * g_0_y_y_yy[i] * c_exps[i] + 4.0 * g_xx_y_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_0_yz[i] = -2.0 * g_0_y_y_yz[i] * c_exps[i] + 4.0 * g_xx_y_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_y_0_zz[i] = -2.0 * g_0_y_y_zz[i] * c_exps[i] + 4.0 * g_xx_y_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_x_0_y_0_x_z_0_xx, g_x_0_y_0_x_z_0_xy, g_x_0_y_0_x_z_0_xz, g_x_0_y_0_x_z_0_yy, g_x_0_y_0_x_z_0_yz, g_x_0_y_0_x_z_0_zz, g_xx_z_y_xx, g_xx_z_y_xy, g_xx_z_y_xz, g_xx_z_y_yy, g_xx_z_y_yz, g_xx_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_x_z_0_xx[i] = -2.0 * g_0_z_y_xx[i] * c_exps[i] + 4.0 * g_xx_z_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_0_xy[i] = -2.0 * g_0_z_y_xy[i] * c_exps[i] + 4.0 * g_xx_z_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_0_xz[i] = -2.0 * g_0_z_y_xz[i] * c_exps[i] + 4.0 * g_xx_z_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_0_yy[i] = -2.0 * g_0_z_y_yy[i] * c_exps[i] + 4.0 * g_xx_z_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_0_yz[i] = -2.0 * g_0_z_y_yz[i] * c_exps[i] + 4.0 * g_xx_z_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_x_z_0_zz[i] = -2.0 * g_0_z_y_zz[i] * c_exps[i] + 4.0 * g_xx_z_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_y_0_y_x_0_xx, g_x_0_y_0_y_x_0_xy, g_x_0_y_0_y_x_0_xz, g_x_0_y_0_y_x_0_yy, g_x_0_y_0_y_x_0_yz, g_x_0_y_0_y_x_0_zz, g_xy_x_y_xx, g_xy_x_y_xy, g_xy_x_y_xz, g_xy_x_y_yy, g_xy_x_y_yz, g_xy_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_x_0_xx[i] = 4.0 * g_xy_x_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_0_xy[i] = 4.0 * g_xy_x_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_0_xz[i] = 4.0 * g_xy_x_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_0_yy[i] = 4.0 * g_xy_x_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_0_yz[i] = 4.0 * g_xy_x_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_x_0_zz[i] = 4.0 * g_xy_x_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_y_0_y_y_0_xx, g_x_0_y_0_y_y_0_xy, g_x_0_y_0_y_y_0_xz, g_x_0_y_0_y_y_0_yy, g_x_0_y_0_y_y_0_yz, g_x_0_y_0_y_y_0_zz, g_xy_y_y_xx, g_xy_y_y_xy, g_xy_y_y_xz, g_xy_y_y_yy, g_xy_y_y_yz, g_xy_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_y_0_xx[i] = 4.0 * g_xy_y_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_0_xy[i] = 4.0 * g_xy_y_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_0_xz[i] = 4.0 * g_xy_y_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_0_yy[i] = 4.0 * g_xy_y_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_0_yz[i] = 4.0 * g_xy_y_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_y_0_zz[i] = 4.0 * g_xy_y_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_y_0_y_z_0_xx, g_x_0_y_0_y_z_0_xy, g_x_0_y_0_y_z_0_xz, g_x_0_y_0_y_z_0_yy, g_x_0_y_0_y_z_0_yz, g_x_0_y_0_y_z_0_zz, g_xy_z_y_xx, g_xy_z_y_xy, g_xy_z_y_xz, g_xy_z_y_yy, g_xy_z_y_yz, g_xy_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_y_z_0_xx[i] = 4.0 * g_xy_z_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_0_xy[i] = 4.0 * g_xy_z_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_0_xz[i] = 4.0 * g_xy_z_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_0_yy[i] = 4.0 * g_xy_z_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_0_yz[i] = 4.0 * g_xy_z_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_y_z_0_zz[i] = 4.0 * g_xy_z_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_y_0_z_x_0_xx, g_x_0_y_0_z_x_0_xy, g_x_0_y_0_z_x_0_xz, g_x_0_y_0_z_x_0_yy, g_x_0_y_0_z_x_0_yz, g_x_0_y_0_z_x_0_zz, g_xz_x_y_xx, g_xz_x_y_xy, g_xz_x_y_xz, g_xz_x_y_yy, g_xz_x_y_yz, g_xz_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_x_0_xx[i] = 4.0 * g_xz_x_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_0_xy[i] = 4.0 * g_xz_x_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_0_xz[i] = 4.0 * g_xz_x_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_0_yy[i] = 4.0 * g_xz_x_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_0_yz[i] = 4.0 * g_xz_x_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_x_0_zz[i] = 4.0 * g_xz_x_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_y_0_z_y_0_xx, g_x_0_y_0_z_y_0_xy, g_x_0_y_0_z_y_0_xz, g_x_0_y_0_z_y_0_yy, g_x_0_y_0_z_y_0_yz, g_x_0_y_0_z_y_0_zz, g_xz_y_y_xx, g_xz_y_y_xy, g_xz_y_y_xz, g_xz_y_y_yy, g_xz_y_y_yz, g_xz_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_y_0_xx[i] = 4.0 * g_xz_y_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_0_xy[i] = 4.0 * g_xz_y_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_0_xz[i] = 4.0 * g_xz_y_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_0_yy[i] = 4.0 * g_xz_y_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_0_yz[i] = 4.0 * g_xz_y_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_y_0_zz[i] = 4.0 * g_xz_y_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_y_0_z_z_0_xx, g_x_0_y_0_z_z_0_xy, g_x_0_y_0_z_z_0_xz, g_x_0_y_0_z_z_0_yy, g_x_0_y_0_z_z_0_yz, g_x_0_y_0_z_z_0_zz, g_xz_z_y_xx, g_xz_z_y_xy, g_xz_z_y_xz, g_xz_z_y_yy, g_xz_z_y_yz, g_xz_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_z_z_0_xx[i] = 4.0 * g_xz_z_y_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_0_xy[i] = 4.0 * g_xz_z_y_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_0_xz[i] = 4.0 * g_xz_z_y_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_0_yy[i] = 4.0 * g_xz_z_y_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_0_yz[i] = 4.0 * g_xz_z_y_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_z_z_0_zz[i] = 4.0 * g_xz_z_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_x_0_z_0_x_x_0_xx, g_x_0_z_0_x_x_0_xy, g_x_0_z_0_x_x_0_xz, g_x_0_z_0_x_x_0_yy, g_x_0_z_0_x_x_0_yz, g_x_0_z_0_x_x_0_zz, g_xx_x_z_xx, g_xx_x_z_xy, g_xx_x_z_xz, g_xx_x_z_yy, g_xx_x_z_yz, g_xx_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_x_0_xx[i] = -2.0 * g_0_x_z_xx[i] * c_exps[i] + 4.0 * g_xx_x_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_0_xy[i] = -2.0 * g_0_x_z_xy[i] * c_exps[i] + 4.0 * g_xx_x_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_0_xz[i] = -2.0 * g_0_x_z_xz[i] * c_exps[i] + 4.0 * g_xx_x_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_0_yy[i] = -2.0 * g_0_x_z_yy[i] * c_exps[i] + 4.0 * g_xx_x_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_0_yz[i] = -2.0 * g_0_x_z_yz[i] * c_exps[i] + 4.0 * g_xx_x_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_x_0_zz[i] = -2.0 * g_0_x_z_zz[i] * c_exps[i] + 4.0 * g_xx_x_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_x_0_z_0_x_y_0_xx, g_x_0_z_0_x_y_0_xy, g_x_0_z_0_x_y_0_xz, g_x_0_z_0_x_y_0_yy, g_x_0_z_0_x_y_0_yz, g_x_0_z_0_x_y_0_zz, g_xx_y_z_xx, g_xx_y_z_xy, g_xx_y_z_xz, g_xx_y_z_yy, g_xx_y_z_yz, g_xx_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_y_0_xx[i] = -2.0 * g_0_y_z_xx[i] * c_exps[i] + 4.0 * g_xx_y_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_0_xy[i] = -2.0 * g_0_y_z_xy[i] * c_exps[i] + 4.0 * g_xx_y_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_0_xz[i] = -2.0 * g_0_y_z_xz[i] * c_exps[i] + 4.0 * g_xx_y_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_0_yy[i] = -2.0 * g_0_y_z_yy[i] * c_exps[i] + 4.0 * g_xx_y_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_0_yz[i] = -2.0 * g_0_y_z_yz[i] * c_exps[i] + 4.0 * g_xx_y_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_y_0_zz[i] = -2.0 * g_0_y_z_zz[i] * c_exps[i] + 4.0 * g_xx_y_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_x_0_z_0_x_z_0_xx, g_x_0_z_0_x_z_0_xy, g_x_0_z_0_x_z_0_xz, g_x_0_z_0_x_z_0_yy, g_x_0_z_0_x_z_0_yz, g_x_0_z_0_x_z_0_zz, g_xx_z_z_xx, g_xx_z_z_xy, g_xx_z_z_xz, g_xx_z_z_yy, g_xx_z_z_yz, g_xx_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_x_z_0_xx[i] = -2.0 * g_0_z_z_xx[i] * c_exps[i] + 4.0 * g_xx_z_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_0_xy[i] = -2.0 * g_0_z_z_xy[i] * c_exps[i] + 4.0 * g_xx_z_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_0_xz[i] = -2.0 * g_0_z_z_xz[i] * c_exps[i] + 4.0 * g_xx_z_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_0_yy[i] = -2.0 * g_0_z_z_yy[i] * c_exps[i] + 4.0 * g_xx_z_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_0_yz[i] = -2.0 * g_0_z_z_yz[i] * c_exps[i] + 4.0 * g_xx_z_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_x_z_0_zz[i] = -2.0 * g_0_z_z_zz[i] * c_exps[i] + 4.0 * g_xx_z_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_0_z_0_y_x_0_xx, g_x_0_z_0_y_x_0_xy, g_x_0_z_0_y_x_0_xz, g_x_0_z_0_y_x_0_yy, g_x_0_z_0_y_x_0_yz, g_x_0_z_0_y_x_0_zz, g_xy_x_z_xx, g_xy_x_z_xy, g_xy_x_z_xz, g_xy_x_z_yy, g_xy_x_z_yz, g_xy_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_x_0_xx[i] = 4.0 * g_xy_x_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_0_xy[i] = 4.0 * g_xy_x_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_0_xz[i] = 4.0 * g_xy_x_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_0_yy[i] = 4.0 * g_xy_x_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_0_yz[i] = 4.0 * g_xy_x_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_x_0_zz[i] = 4.0 * g_xy_x_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_0_z_0_y_y_0_xx, g_x_0_z_0_y_y_0_xy, g_x_0_z_0_y_y_0_xz, g_x_0_z_0_y_y_0_yy, g_x_0_z_0_y_y_0_yz, g_x_0_z_0_y_y_0_zz, g_xy_y_z_xx, g_xy_y_z_xy, g_xy_y_z_xz, g_xy_y_z_yy, g_xy_y_z_yz, g_xy_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_y_0_xx[i] = 4.0 * g_xy_y_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_0_xy[i] = 4.0 * g_xy_y_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_0_xz[i] = 4.0 * g_xy_y_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_0_yy[i] = 4.0 * g_xy_y_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_0_yz[i] = 4.0 * g_xy_y_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_y_0_zz[i] = 4.0 * g_xy_y_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_0_z_0_y_z_0_xx, g_x_0_z_0_y_z_0_xy, g_x_0_z_0_y_z_0_xz, g_x_0_z_0_y_z_0_yy, g_x_0_z_0_y_z_0_yz, g_x_0_z_0_y_z_0_zz, g_xy_z_z_xx, g_xy_z_z_xy, g_xy_z_z_xz, g_xy_z_z_yy, g_xy_z_z_yz, g_xy_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_y_z_0_xx[i] = 4.0 * g_xy_z_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_0_xy[i] = 4.0 * g_xy_z_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_0_xz[i] = 4.0 * g_xy_z_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_0_yy[i] = 4.0 * g_xy_z_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_0_yz[i] = 4.0 * g_xy_z_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_y_z_0_zz[i] = 4.0 * g_xy_z_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_z_0_z_x_0_xx, g_x_0_z_0_z_x_0_xy, g_x_0_z_0_z_x_0_xz, g_x_0_z_0_z_x_0_yy, g_x_0_z_0_z_x_0_yz, g_x_0_z_0_z_x_0_zz, g_xz_x_z_xx, g_xz_x_z_xy, g_xz_x_z_xz, g_xz_x_z_yy, g_xz_x_z_yz, g_xz_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_x_0_xx[i] = 4.0 * g_xz_x_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_0_xy[i] = 4.0 * g_xz_x_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_0_xz[i] = 4.0 * g_xz_x_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_0_yy[i] = 4.0 * g_xz_x_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_0_yz[i] = 4.0 * g_xz_x_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_x_0_zz[i] = 4.0 * g_xz_x_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_z_0_z_y_0_xx, g_x_0_z_0_z_y_0_xy, g_x_0_z_0_z_y_0_xz, g_x_0_z_0_z_y_0_yy, g_x_0_z_0_z_y_0_yz, g_x_0_z_0_z_y_0_zz, g_xz_y_z_xx, g_xz_y_z_xy, g_xz_y_z_xz, g_xz_y_z_yy, g_xz_y_z_yz, g_xz_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_y_0_xx[i] = 4.0 * g_xz_y_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_0_xy[i] = 4.0 * g_xz_y_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_0_xz[i] = 4.0 * g_xz_y_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_0_yy[i] = 4.0 * g_xz_y_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_0_yz[i] = 4.0 * g_xz_y_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_y_0_zz[i] = 4.0 * g_xz_y_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_z_0_z_z_0_xx, g_x_0_z_0_z_z_0_xy, g_x_0_z_0_z_z_0_xz, g_x_0_z_0_z_z_0_yy, g_x_0_z_0_z_z_0_yz, g_x_0_z_0_z_z_0_zz, g_xz_z_z_xx, g_xz_z_z_xy, g_xz_z_z_xz, g_xz_z_z_yy, g_xz_z_z_yz, g_xz_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_z_z_0_xx[i] = 4.0 * g_xz_z_z_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_0_xy[i] = 4.0 * g_xz_z_z_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_0_xz[i] = 4.0 * g_xz_z_z_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_0_yy[i] = 4.0 * g_xz_z_z_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_0_yz[i] = 4.0 * g_xz_z_z_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_z_z_0_zz[i] = 4.0 * g_xz_z_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_xy_x_x_xx, g_xy_x_x_xy, g_xy_x_x_xz, g_xy_x_x_yy, g_xy_x_x_yz, g_xy_x_x_zz, g_y_0_x_0_x_x_0_xx, g_y_0_x_0_x_x_0_xy, g_y_0_x_0_x_x_0_xz, g_y_0_x_0_x_x_0_yy, g_y_0_x_0_x_x_0_yz, g_y_0_x_0_x_x_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_x_0_xx[i] = 4.0 * g_xy_x_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_0_xy[i] = 4.0 * g_xy_x_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_0_xz[i] = 4.0 * g_xy_x_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_0_yy[i] = 4.0 * g_xy_x_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_0_yz[i] = 4.0 * g_xy_x_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_x_0_zz[i] = 4.0 * g_xy_x_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_xy_y_x_xx, g_xy_y_x_xy, g_xy_y_x_xz, g_xy_y_x_yy, g_xy_y_x_yz, g_xy_y_x_zz, g_y_0_x_0_x_y_0_xx, g_y_0_x_0_x_y_0_xy, g_y_0_x_0_x_y_0_xz, g_y_0_x_0_x_y_0_yy, g_y_0_x_0_x_y_0_yz, g_y_0_x_0_x_y_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_y_0_xx[i] = 4.0 * g_xy_y_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_0_xy[i] = 4.0 * g_xy_y_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_0_xz[i] = 4.0 * g_xy_y_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_0_yy[i] = 4.0 * g_xy_y_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_0_yz[i] = 4.0 * g_xy_y_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_y_0_zz[i] = 4.0 * g_xy_y_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_xy_z_x_xx, g_xy_z_x_xy, g_xy_z_x_xz, g_xy_z_x_yy, g_xy_z_x_yz, g_xy_z_x_zz, g_y_0_x_0_x_z_0_xx, g_y_0_x_0_x_z_0_xy, g_y_0_x_0_x_z_0_xz, g_y_0_x_0_x_z_0_yy, g_y_0_x_0_x_z_0_yz, g_y_0_x_0_x_z_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_x_z_0_xx[i] = 4.0 * g_xy_z_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_0_xy[i] = 4.0 * g_xy_z_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_0_xz[i] = 4.0 * g_xy_z_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_0_yy[i] = 4.0 * g_xy_z_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_0_yz[i] = 4.0 * g_xy_z_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_x_z_0_zz[i] = 4.0 * g_xy_z_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_y_0_x_0_y_x_0_xx, g_y_0_x_0_y_x_0_xy, g_y_0_x_0_y_x_0_xz, g_y_0_x_0_y_x_0_yy, g_y_0_x_0_y_x_0_yz, g_y_0_x_0_y_x_0_zz, g_yy_x_x_xx, g_yy_x_x_xy, g_yy_x_x_xz, g_yy_x_x_yy, g_yy_x_x_yz, g_yy_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_x_0_xx[i] = -2.0 * g_0_x_x_xx[i] * c_exps[i] + 4.0 * g_yy_x_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_0_xy[i] = -2.0 * g_0_x_x_xy[i] * c_exps[i] + 4.0 * g_yy_x_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_0_xz[i] = -2.0 * g_0_x_x_xz[i] * c_exps[i] + 4.0 * g_yy_x_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_0_yy[i] = -2.0 * g_0_x_x_yy[i] * c_exps[i] + 4.0 * g_yy_x_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_0_yz[i] = -2.0 * g_0_x_x_yz[i] * c_exps[i] + 4.0 * g_yy_x_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_x_0_zz[i] = -2.0 * g_0_x_x_zz[i] * c_exps[i] + 4.0 * g_yy_x_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_y_0_x_0_y_y_0_xx, g_y_0_x_0_y_y_0_xy, g_y_0_x_0_y_y_0_xz, g_y_0_x_0_y_y_0_yy, g_y_0_x_0_y_y_0_yz, g_y_0_x_0_y_y_0_zz, g_yy_y_x_xx, g_yy_y_x_xy, g_yy_y_x_xz, g_yy_y_x_yy, g_yy_y_x_yz, g_yy_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_y_0_xx[i] = -2.0 * g_0_y_x_xx[i] * c_exps[i] + 4.0 * g_yy_y_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_0_xy[i] = -2.0 * g_0_y_x_xy[i] * c_exps[i] + 4.0 * g_yy_y_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_0_xz[i] = -2.0 * g_0_y_x_xz[i] * c_exps[i] + 4.0 * g_yy_y_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_0_yy[i] = -2.0 * g_0_y_x_yy[i] * c_exps[i] + 4.0 * g_yy_y_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_0_yz[i] = -2.0 * g_0_y_x_yz[i] * c_exps[i] + 4.0 * g_yy_y_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_y_0_zz[i] = -2.0 * g_0_y_x_zz[i] * c_exps[i] + 4.0 * g_yy_y_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_y_0_x_0_y_z_0_xx, g_y_0_x_0_y_z_0_xy, g_y_0_x_0_y_z_0_xz, g_y_0_x_0_y_z_0_yy, g_y_0_x_0_y_z_0_yz, g_y_0_x_0_y_z_0_zz, g_yy_z_x_xx, g_yy_z_x_xy, g_yy_z_x_xz, g_yy_z_x_yy, g_yy_z_x_yz, g_yy_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_y_z_0_xx[i] = -2.0 * g_0_z_x_xx[i] * c_exps[i] + 4.0 * g_yy_z_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_0_xy[i] = -2.0 * g_0_z_x_xy[i] * c_exps[i] + 4.0 * g_yy_z_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_0_xz[i] = -2.0 * g_0_z_x_xz[i] * c_exps[i] + 4.0 * g_yy_z_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_0_yy[i] = -2.0 * g_0_z_x_yy[i] * c_exps[i] + 4.0 * g_yy_z_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_0_yz[i] = -2.0 * g_0_z_x_yz[i] * c_exps[i] + 4.0 * g_yy_z_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_y_z_0_zz[i] = -2.0 * g_0_z_x_zz[i] * c_exps[i] + 4.0 * g_yy_z_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_y_0_x_0_z_x_0_xx, g_y_0_x_0_z_x_0_xy, g_y_0_x_0_z_x_0_xz, g_y_0_x_0_z_x_0_yy, g_y_0_x_0_z_x_0_yz, g_y_0_x_0_z_x_0_zz, g_yz_x_x_xx, g_yz_x_x_xy, g_yz_x_x_xz, g_yz_x_x_yy, g_yz_x_x_yz, g_yz_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_x_0_xx[i] = 4.0 * g_yz_x_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_0_xy[i] = 4.0 * g_yz_x_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_0_xz[i] = 4.0 * g_yz_x_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_0_yy[i] = 4.0 * g_yz_x_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_0_yz[i] = 4.0 * g_yz_x_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_x_0_zz[i] = 4.0 * g_yz_x_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_y_0_x_0_z_y_0_xx, g_y_0_x_0_z_y_0_xy, g_y_0_x_0_z_y_0_xz, g_y_0_x_0_z_y_0_yy, g_y_0_x_0_z_y_0_yz, g_y_0_x_0_z_y_0_zz, g_yz_y_x_xx, g_yz_y_x_xy, g_yz_y_x_xz, g_yz_y_x_yy, g_yz_y_x_yz, g_yz_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_y_0_xx[i] = 4.0 * g_yz_y_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_0_xy[i] = 4.0 * g_yz_y_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_0_xz[i] = 4.0 * g_yz_y_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_0_yy[i] = 4.0 * g_yz_y_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_0_yz[i] = 4.0 * g_yz_y_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_y_0_zz[i] = 4.0 * g_yz_y_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_y_0_x_0_z_z_0_xx, g_y_0_x_0_z_z_0_xy, g_y_0_x_0_z_z_0_xz, g_y_0_x_0_z_z_0_yy, g_y_0_x_0_z_z_0_yz, g_y_0_x_0_z_z_0_zz, g_yz_z_x_xx, g_yz_z_x_xy, g_yz_z_x_xz, g_yz_z_x_yy, g_yz_z_x_yz, g_yz_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_z_z_0_xx[i] = 4.0 * g_yz_z_x_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_0_xy[i] = 4.0 * g_yz_z_x_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_0_xz[i] = 4.0 * g_yz_z_x_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_0_yy[i] = 4.0 * g_yz_z_x_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_0_yz[i] = 4.0 * g_yz_z_x_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_z_z_0_zz[i] = 4.0 * g_yz_z_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_xy_x_y_xx, g_xy_x_y_xy, g_xy_x_y_xz, g_xy_x_y_yy, g_xy_x_y_yz, g_xy_x_y_zz, g_y_0_y_0_x_x_0_xx, g_y_0_y_0_x_x_0_xy, g_y_0_y_0_x_x_0_xz, g_y_0_y_0_x_x_0_yy, g_y_0_y_0_x_x_0_yz, g_y_0_y_0_x_x_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_x_0_xx[i] = 4.0 * g_xy_x_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_0_xy[i] = 4.0 * g_xy_x_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_0_xz[i] = 4.0 * g_xy_x_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_0_yy[i] = 4.0 * g_xy_x_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_0_yz[i] = 4.0 * g_xy_x_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_x_0_zz[i] = 4.0 * g_xy_x_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_xy_y_y_xx, g_xy_y_y_xy, g_xy_y_y_xz, g_xy_y_y_yy, g_xy_y_y_yz, g_xy_y_y_zz, g_y_0_y_0_x_y_0_xx, g_y_0_y_0_x_y_0_xy, g_y_0_y_0_x_y_0_xz, g_y_0_y_0_x_y_0_yy, g_y_0_y_0_x_y_0_yz, g_y_0_y_0_x_y_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_y_0_xx[i] = 4.0 * g_xy_y_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_0_xy[i] = 4.0 * g_xy_y_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_0_xz[i] = 4.0 * g_xy_y_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_0_yy[i] = 4.0 * g_xy_y_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_0_yz[i] = 4.0 * g_xy_y_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_y_0_zz[i] = 4.0 * g_xy_y_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_xy_z_y_xx, g_xy_z_y_xy, g_xy_z_y_xz, g_xy_z_y_yy, g_xy_z_y_yz, g_xy_z_y_zz, g_y_0_y_0_x_z_0_xx, g_y_0_y_0_x_z_0_xy, g_y_0_y_0_x_z_0_xz, g_y_0_y_0_x_z_0_yy, g_y_0_y_0_x_z_0_yz, g_y_0_y_0_x_z_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_x_z_0_xx[i] = 4.0 * g_xy_z_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_0_xy[i] = 4.0 * g_xy_z_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_0_xz[i] = 4.0 * g_xy_z_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_0_yy[i] = 4.0 * g_xy_z_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_0_yz[i] = 4.0 * g_xy_z_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_x_z_0_zz[i] = 4.0 * g_xy_z_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_y_0_y_0_y_x_0_xx, g_y_0_y_0_y_x_0_xy, g_y_0_y_0_y_x_0_xz, g_y_0_y_0_y_x_0_yy, g_y_0_y_0_y_x_0_yz, g_y_0_y_0_y_x_0_zz, g_yy_x_y_xx, g_yy_x_y_xy, g_yy_x_y_xz, g_yy_x_y_yy, g_yy_x_y_yz, g_yy_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_x_0_xx[i] = -2.0 * g_0_x_y_xx[i] * c_exps[i] + 4.0 * g_yy_x_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_0_xy[i] = -2.0 * g_0_x_y_xy[i] * c_exps[i] + 4.0 * g_yy_x_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_0_xz[i] = -2.0 * g_0_x_y_xz[i] * c_exps[i] + 4.0 * g_yy_x_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_0_yy[i] = -2.0 * g_0_x_y_yy[i] * c_exps[i] + 4.0 * g_yy_x_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_0_yz[i] = -2.0 * g_0_x_y_yz[i] * c_exps[i] + 4.0 * g_yy_x_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_x_0_zz[i] = -2.0 * g_0_x_y_zz[i] * c_exps[i] + 4.0 * g_yy_x_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_y_0_y_0_y_y_0_xx, g_y_0_y_0_y_y_0_xy, g_y_0_y_0_y_y_0_xz, g_y_0_y_0_y_y_0_yy, g_y_0_y_0_y_y_0_yz, g_y_0_y_0_y_y_0_zz, g_yy_y_y_xx, g_yy_y_y_xy, g_yy_y_y_xz, g_yy_y_y_yy, g_yy_y_y_yz, g_yy_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_y_0_xx[i] = -2.0 * g_0_y_y_xx[i] * c_exps[i] + 4.0 * g_yy_y_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_0_xy[i] = -2.0 * g_0_y_y_xy[i] * c_exps[i] + 4.0 * g_yy_y_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_0_xz[i] = -2.0 * g_0_y_y_xz[i] * c_exps[i] + 4.0 * g_yy_y_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_0_yy[i] = -2.0 * g_0_y_y_yy[i] * c_exps[i] + 4.0 * g_yy_y_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_0_yz[i] = -2.0 * g_0_y_y_yz[i] * c_exps[i] + 4.0 * g_yy_y_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_y_0_zz[i] = -2.0 * g_0_y_y_zz[i] * c_exps[i] + 4.0 * g_yy_y_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_y_0_y_0_y_z_0_xx, g_y_0_y_0_y_z_0_xy, g_y_0_y_0_y_z_0_xz, g_y_0_y_0_y_z_0_yy, g_y_0_y_0_y_z_0_yz, g_y_0_y_0_y_z_0_zz, g_yy_z_y_xx, g_yy_z_y_xy, g_yy_z_y_xz, g_yy_z_y_yy, g_yy_z_y_yz, g_yy_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_y_z_0_xx[i] = -2.0 * g_0_z_y_xx[i] * c_exps[i] + 4.0 * g_yy_z_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_0_xy[i] = -2.0 * g_0_z_y_xy[i] * c_exps[i] + 4.0 * g_yy_z_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_0_xz[i] = -2.0 * g_0_z_y_xz[i] * c_exps[i] + 4.0 * g_yy_z_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_0_yy[i] = -2.0 * g_0_z_y_yy[i] * c_exps[i] + 4.0 * g_yy_z_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_0_yz[i] = -2.0 * g_0_z_y_yz[i] * c_exps[i] + 4.0 * g_yy_z_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_y_z_0_zz[i] = -2.0 * g_0_z_y_zz[i] * c_exps[i] + 4.0 * g_yy_z_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_y_0_y_0_z_x_0_xx, g_y_0_y_0_z_x_0_xy, g_y_0_y_0_z_x_0_xz, g_y_0_y_0_z_x_0_yy, g_y_0_y_0_z_x_0_yz, g_y_0_y_0_z_x_0_zz, g_yz_x_y_xx, g_yz_x_y_xy, g_yz_x_y_xz, g_yz_x_y_yy, g_yz_x_y_yz, g_yz_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_x_0_xx[i] = 4.0 * g_yz_x_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_0_xy[i] = 4.0 * g_yz_x_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_0_xz[i] = 4.0 * g_yz_x_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_0_yy[i] = 4.0 * g_yz_x_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_0_yz[i] = 4.0 * g_yz_x_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_x_0_zz[i] = 4.0 * g_yz_x_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_y_0_y_0_z_y_0_xx, g_y_0_y_0_z_y_0_xy, g_y_0_y_0_z_y_0_xz, g_y_0_y_0_z_y_0_yy, g_y_0_y_0_z_y_0_yz, g_y_0_y_0_z_y_0_zz, g_yz_y_y_xx, g_yz_y_y_xy, g_yz_y_y_xz, g_yz_y_y_yy, g_yz_y_y_yz, g_yz_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_y_0_xx[i] = 4.0 * g_yz_y_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_0_xy[i] = 4.0 * g_yz_y_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_0_xz[i] = 4.0 * g_yz_y_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_0_yy[i] = 4.0 * g_yz_y_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_0_yz[i] = 4.0 * g_yz_y_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_y_0_zz[i] = 4.0 * g_yz_y_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_y_0_y_0_z_z_0_xx, g_y_0_y_0_z_z_0_xy, g_y_0_y_0_z_z_0_xz, g_y_0_y_0_z_z_0_yy, g_y_0_y_0_z_z_0_yz, g_y_0_y_0_z_z_0_zz, g_yz_z_y_xx, g_yz_z_y_xy, g_yz_z_y_xz, g_yz_z_y_yy, g_yz_z_y_yz, g_yz_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_z_z_0_xx[i] = 4.0 * g_yz_z_y_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_0_xy[i] = 4.0 * g_yz_z_y_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_0_xz[i] = 4.0 * g_yz_z_y_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_0_yy[i] = 4.0 * g_yz_z_y_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_0_yz[i] = 4.0 * g_yz_z_y_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_z_z_0_zz[i] = 4.0 * g_yz_z_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_xy_x_z_xx, g_xy_x_z_xy, g_xy_x_z_xz, g_xy_x_z_yy, g_xy_x_z_yz, g_xy_x_z_zz, g_y_0_z_0_x_x_0_xx, g_y_0_z_0_x_x_0_xy, g_y_0_z_0_x_x_0_xz, g_y_0_z_0_x_x_0_yy, g_y_0_z_0_x_x_0_yz, g_y_0_z_0_x_x_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_x_0_xx[i] = 4.0 * g_xy_x_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_0_xy[i] = 4.0 * g_xy_x_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_0_xz[i] = 4.0 * g_xy_x_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_0_yy[i] = 4.0 * g_xy_x_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_0_yz[i] = 4.0 * g_xy_x_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_x_0_zz[i] = 4.0 * g_xy_x_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_xy_y_z_xx, g_xy_y_z_xy, g_xy_y_z_xz, g_xy_y_z_yy, g_xy_y_z_yz, g_xy_y_z_zz, g_y_0_z_0_x_y_0_xx, g_y_0_z_0_x_y_0_xy, g_y_0_z_0_x_y_0_xz, g_y_0_z_0_x_y_0_yy, g_y_0_z_0_x_y_0_yz, g_y_0_z_0_x_y_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_y_0_xx[i] = 4.0 * g_xy_y_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_0_xy[i] = 4.0 * g_xy_y_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_0_xz[i] = 4.0 * g_xy_y_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_0_yy[i] = 4.0 * g_xy_y_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_0_yz[i] = 4.0 * g_xy_y_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_y_0_zz[i] = 4.0 * g_xy_y_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_xy_z_z_xx, g_xy_z_z_xy, g_xy_z_z_xz, g_xy_z_z_yy, g_xy_z_z_yz, g_xy_z_z_zz, g_y_0_z_0_x_z_0_xx, g_y_0_z_0_x_z_0_xy, g_y_0_z_0_x_z_0_xz, g_y_0_z_0_x_z_0_yy, g_y_0_z_0_x_z_0_yz, g_y_0_z_0_x_z_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_x_z_0_xx[i] = 4.0 * g_xy_z_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_0_xy[i] = 4.0 * g_xy_z_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_0_xz[i] = 4.0 * g_xy_z_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_0_yy[i] = 4.0 * g_xy_z_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_0_yz[i] = 4.0 * g_xy_z_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_x_z_0_zz[i] = 4.0 * g_xy_z_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_y_0_z_0_y_x_0_xx, g_y_0_z_0_y_x_0_xy, g_y_0_z_0_y_x_0_xz, g_y_0_z_0_y_x_0_yy, g_y_0_z_0_y_x_0_yz, g_y_0_z_0_y_x_0_zz, g_yy_x_z_xx, g_yy_x_z_xy, g_yy_x_z_xz, g_yy_x_z_yy, g_yy_x_z_yz, g_yy_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_x_0_xx[i] = -2.0 * g_0_x_z_xx[i] * c_exps[i] + 4.0 * g_yy_x_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_0_xy[i] = -2.0 * g_0_x_z_xy[i] * c_exps[i] + 4.0 * g_yy_x_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_0_xz[i] = -2.0 * g_0_x_z_xz[i] * c_exps[i] + 4.0 * g_yy_x_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_0_yy[i] = -2.0 * g_0_x_z_yy[i] * c_exps[i] + 4.0 * g_yy_x_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_0_yz[i] = -2.0 * g_0_x_z_yz[i] * c_exps[i] + 4.0 * g_yy_x_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_x_0_zz[i] = -2.0 * g_0_x_z_zz[i] * c_exps[i] + 4.0 * g_yy_x_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_y_0_z_0_y_y_0_xx, g_y_0_z_0_y_y_0_xy, g_y_0_z_0_y_y_0_xz, g_y_0_z_0_y_y_0_yy, g_y_0_z_0_y_y_0_yz, g_y_0_z_0_y_y_0_zz, g_yy_y_z_xx, g_yy_y_z_xy, g_yy_y_z_xz, g_yy_y_z_yy, g_yy_y_z_yz, g_yy_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_y_0_xx[i] = -2.0 * g_0_y_z_xx[i] * c_exps[i] + 4.0 * g_yy_y_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_0_xy[i] = -2.0 * g_0_y_z_xy[i] * c_exps[i] + 4.0 * g_yy_y_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_0_xz[i] = -2.0 * g_0_y_z_xz[i] * c_exps[i] + 4.0 * g_yy_y_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_0_yy[i] = -2.0 * g_0_y_z_yy[i] * c_exps[i] + 4.0 * g_yy_y_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_0_yz[i] = -2.0 * g_0_y_z_yz[i] * c_exps[i] + 4.0 * g_yy_y_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_y_0_zz[i] = -2.0 * g_0_y_z_zz[i] * c_exps[i] + 4.0 * g_yy_y_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_y_0_z_0_y_z_0_xx, g_y_0_z_0_y_z_0_xy, g_y_0_z_0_y_z_0_xz, g_y_0_z_0_y_z_0_yy, g_y_0_z_0_y_z_0_yz, g_y_0_z_0_y_z_0_zz, g_yy_z_z_xx, g_yy_z_z_xy, g_yy_z_z_xz, g_yy_z_z_yy, g_yy_z_z_yz, g_yy_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_y_z_0_xx[i] = -2.0 * g_0_z_z_xx[i] * c_exps[i] + 4.0 * g_yy_z_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_0_xy[i] = -2.0 * g_0_z_z_xy[i] * c_exps[i] + 4.0 * g_yy_z_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_0_xz[i] = -2.0 * g_0_z_z_xz[i] * c_exps[i] + 4.0 * g_yy_z_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_0_yy[i] = -2.0 * g_0_z_z_yy[i] * c_exps[i] + 4.0 * g_yy_z_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_0_yz[i] = -2.0 * g_0_z_z_yz[i] * c_exps[i] + 4.0 * g_yy_z_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_y_z_0_zz[i] = -2.0 * g_0_z_z_zz[i] * c_exps[i] + 4.0 * g_yy_z_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_y_0_z_0_z_x_0_xx, g_y_0_z_0_z_x_0_xy, g_y_0_z_0_z_x_0_xz, g_y_0_z_0_z_x_0_yy, g_y_0_z_0_z_x_0_yz, g_y_0_z_0_z_x_0_zz, g_yz_x_z_xx, g_yz_x_z_xy, g_yz_x_z_xz, g_yz_x_z_yy, g_yz_x_z_yz, g_yz_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_x_0_xx[i] = 4.0 * g_yz_x_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_0_xy[i] = 4.0 * g_yz_x_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_0_xz[i] = 4.0 * g_yz_x_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_0_yy[i] = 4.0 * g_yz_x_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_0_yz[i] = 4.0 * g_yz_x_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_x_0_zz[i] = 4.0 * g_yz_x_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_y_0_z_0_z_y_0_xx, g_y_0_z_0_z_y_0_xy, g_y_0_z_0_z_y_0_xz, g_y_0_z_0_z_y_0_yy, g_y_0_z_0_z_y_0_yz, g_y_0_z_0_z_y_0_zz, g_yz_y_z_xx, g_yz_y_z_xy, g_yz_y_z_xz, g_yz_y_z_yy, g_yz_y_z_yz, g_yz_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_y_0_xx[i] = 4.0 * g_yz_y_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_0_xy[i] = 4.0 * g_yz_y_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_0_xz[i] = 4.0 * g_yz_y_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_0_yy[i] = 4.0 * g_yz_y_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_0_yz[i] = 4.0 * g_yz_y_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_y_0_zz[i] = 4.0 * g_yz_y_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_y_0_z_0_z_z_0_xx, g_y_0_z_0_z_z_0_xy, g_y_0_z_0_z_z_0_xz, g_y_0_z_0_z_z_0_yy, g_y_0_z_0_z_z_0_yz, g_y_0_z_0_z_z_0_zz, g_yz_z_z_xx, g_yz_z_z_xy, g_yz_z_z_xz, g_yz_z_z_yy, g_yz_z_z_yz, g_yz_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_z_z_0_xx[i] = 4.0 * g_yz_z_z_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_0_xy[i] = 4.0 * g_yz_z_z_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_0_xz[i] = 4.0 * g_yz_z_z_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_0_yy[i] = 4.0 * g_yz_z_z_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_0_yz[i] = 4.0 * g_yz_z_z_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_z_z_0_zz[i] = 4.0 * g_yz_z_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_xz_x_x_xx, g_xz_x_x_xy, g_xz_x_x_xz, g_xz_x_x_yy, g_xz_x_x_yz, g_xz_x_x_zz, g_z_0_x_0_x_x_0_xx, g_z_0_x_0_x_x_0_xy, g_z_0_x_0_x_x_0_xz, g_z_0_x_0_x_x_0_yy, g_z_0_x_0_x_x_0_yz, g_z_0_x_0_x_x_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_x_0_xx[i] = 4.0 * g_xz_x_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_0_xy[i] = 4.0 * g_xz_x_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_0_xz[i] = 4.0 * g_xz_x_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_0_yy[i] = 4.0 * g_xz_x_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_0_yz[i] = 4.0 * g_xz_x_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_x_0_zz[i] = 4.0 * g_xz_x_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_xz_y_x_xx, g_xz_y_x_xy, g_xz_y_x_xz, g_xz_y_x_yy, g_xz_y_x_yz, g_xz_y_x_zz, g_z_0_x_0_x_y_0_xx, g_z_0_x_0_x_y_0_xy, g_z_0_x_0_x_y_0_xz, g_z_0_x_0_x_y_0_yy, g_z_0_x_0_x_y_0_yz, g_z_0_x_0_x_y_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_y_0_xx[i] = 4.0 * g_xz_y_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_0_xy[i] = 4.0 * g_xz_y_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_0_xz[i] = 4.0 * g_xz_y_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_0_yy[i] = 4.0 * g_xz_y_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_0_yz[i] = 4.0 * g_xz_y_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_y_0_zz[i] = 4.0 * g_xz_y_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_xz_z_x_xx, g_xz_z_x_xy, g_xz_z_x_xz, g_xz_z_x_yy, g_xz_z_x_yz, g_xz_z_x_zz, g_z_0_x_0_x_z_0_xx, g_z_0_x_0_x_z_0_xy, g_z_0_x_0_x_z_0_xz, g_z_0_x_0_x_z_0_yy, g_z_0_x_0_x_z_0_yz, g_z_0_x_0_x_z_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_x_z_0_xx[i] = 4.0 * g_xz_z_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_0_xy[i] = 4.0 * g_xz_z_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_0_xz[i] = 4.0 * g_xz_z_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_0_yy[i] = 4.0 * g_xz_z_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_0_yz[i] = 4.0 * g_xz_z_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_x_z_0_zz[i] = 4.0 * g_xz_z_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_yz_x_x_xx, g_yz_x_x_xy, g_yz_x_x_xz, g_yz_x_x_yy, g_yz_x_x_yz, g_yz_x_x_zz, g_z_0_x_0_y_x_0_xx, g_z_0_x_0_y_x_0_xy, g_z_0_x_0_y_x_0_xz, g_z_0_x_0_y_x_0_yy, g_z_0_x_0_y_x_0_yz, g_z_0_x_0_y_x_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_x_0_xx[i] = 4.0 * g_yz_x_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_0_xy[i] = 4.0 * g_yz_x_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_0_xz[i] = 4.0 * g_yz_x_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_0_yy[i] = 4.0 * g_yz_x_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_0_yz[i] = 4.0 * g_yz_x_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_x_0_zz[i] = 4.0 * g_yz_x_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_yz_y_x_xx, g_yz_y_x_xy, g_yz_y_x_xz, g_yz_y_x_yy, g_yz_y_x_yz, g_yz_y_x_zz, g_z_0_x_0_y_y_0_xx, g_z_0_x_0_y_y_0_xy, g_z_0_x_0_y_y_0_xz, g_z_0_x_0_y_y_0_yy, g_z_0_x_0_y_y_0_yz, g_z_0_x_0_y_y_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_y_0_xx[i] = 4.0 * g_yz_y_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_0_xy[i] = 4.0 * g_yz_y_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_0_xz[i] = 4.0 * g_yz_y_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_0_yy[i] = 4.0 * g_yz_y_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_0_yz[i] = 4.0 * g_yz_y_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_y_0_zz[i] = 4.0 * g_yz_y_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_yz_z_x_xx, g_yz_z_x_xy, g_yz_z_x_xz, g_yz_z_x_yy, g_yz_z_x_yz, g_yz_z_x_zz, g_z_0_x_0_y_z_0_xx, g_z_0_x_0_y_z_0_xy, g_z_0_x_0_y_z_0_xz, g_z_0_x_0_y_z_0_yy, g_z_0_x_0_y_z_0_yz, g_z_0_x_0_y_z_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_y_z_0_xx[i] = 4.0 * g_yz_z_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_0_xy[i] = 4.0 * g_yz_z_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_0_xz[i] = 4.0 * g_yz_z_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_0_yy[i] = 4.0 * g_yz_z_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_0_yz[i] = 4.0 * g_yz_z_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_y_z_0_zz[i] = 4.0 * g_yz_z_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_z_0_x_0_z_x_0_xx, g_z_0_x_0_z_x_0_xy, g_z_0_x_0_z_x_0_xz, g_z_0_x_0_z_x_0_yy, g_z_0_x_0_z_x_0_yz, g_z_0_x_0_z_x_0_zz, g_zz_x_x_xx, g_zz_x_x_xy, g_zz_x_x_xz, g_zz_x_x_yy, g_zz_x_x_yz, g_zz_x_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_x_0_xx[i] = -2.0 * g_0_x_x_xx[i] * c_exps[i] + 4.0 * g_zz_x_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_0_xy[i] = -2.0 * g_0_x_x_xy[i] * c_exps[i] + 4.0 * g_zz_x_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_0_xz[i] = -2.0 * g_0_x_x_xz[i] * c_exps[i] + 4.0 * g_zz_x_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_0_yy[i] = -2.0 * g_0_x_x_yy[i] * c_exps[i] + 4.0 * g_zz_x_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_0_yz[i] = -2.0 * g_0_x_x_yz[i] * c_exps[i] + 4.0 * g_zz_x_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_x_0_zz[i] = -2.0 * g_0_x_x_zz[i] * c_exps[i] + 4.0 * g_zz_x_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz, g_z_0_x_0_z_y_0_xx, g_z_0_x_0_z_y_0_xy, g_z_0_x_0_z_y_0_xz, g_z_0_x_0_z_y_0_yy, g_z_0_x_0_z_y_0_yz, g_z_0_x_0_z_y_0_zz, g_zz_y_x_xx, g_zz_y_x_xy, g_zz_y_x_xz, g_zz_y_x_yy, g_zz_y_x_yz, g_zz_y_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_y_0_xx[i] = -2.0 * g_0_y_x_xx[i] * c_exps[i] + 4.0 * g_zz_y_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_0_xy[i] = -2.0 * g_0_y_x_xy[i] * c_exps[i] + 4.0 * g_zz_y_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_0_xz[i] = -2.0 * g_0_y_x_xz[i] * c_exps[i] + 4.0 * g_zz_y_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_0_yy[i] = -2.0 * g_0_y_x_yy[i] * c_exps[i] + 4.0 * g_zz_y_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_0_yz[i] = -2.0 * g_0_y_x_yz[i] * c_exps[i] + 4.0 * g_zz_y_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_y_0_zz[i] = -2.0 * g_0_y_x_zz[i] * c_exps[i] + 4.0 * g_zz_y_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz, g_z_0_x_0_z_z_0_xx, g_z_0_x_0_z_z_0_xy, g_z_0_x_0_z_z_0_xz, g_z_0_x_0_z_z_0_yy, g_z_0_x_0_z_z_0_yz, g_z_0_x_0_z_z_0_zz, g_zz_z_x_xx, g_zz_z_x_xy, g_zz_z_x_xz, g_zz_z_x_yy, g_zz_z_x_yz, g_zz_z_x_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_z_z_0_xx[i] = -2.0 * g_0_z_x_xx[i] * c_exps[i] + 4.0 * g_zz_z_x_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_0_xy[i] = -2.0 * g_0_z_x_xy[i] * c_exps[i] + 4.0 * g_zz_z_x_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_0_xz[i] = -2.0 * g_0_z_x_xz[i] * c_exps[i] + 4.0 * g_zz_z_x_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_0_yy[i] = -2.0 * g_0_z_x_yy[i] * c_exps[i] + 4.0 * g_zz_z_x_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_0_yz[i] = -2.0 * g_0_z_x_yz[i] * c_exps[i] + 4.0 * g_zz_z_x_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_z_z_0_zz[i] = -2.0 * g_0_z_x_zz[i] * c_exps[i] + 4.0 * g_zz_z_x_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_xz_x_y_xx, g_xz_x_y_xy, g_xz_x_y_xz, g_xz_x_y_yy, g_xz_x_y_yz, g_xz_x_y_zz, g_z_0_y_0_x_x_0_xx, g_z_0_y_0_x_x_0_xy, g_z_0_y_0_x_x_0_xz, g_z_0_y_0_x_x_0_yy, g_z_0_y_0_x_x_0_yz, g_z_0_y_0_x_x_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_x_0_xx[i] = 4.0 * g_xz_x_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_0_xy[i] = 4.0 * g_xz_x_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_0_xz[i] = 4.0 * g_xz_x_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_0_yy[i] = 4.0 * g_xz_x_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_0_yz[i] = 4.0 * g_xz_x_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_x_0_zz[i] = 4.0 * g_xz_x_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_xz_y_y_xx, g_xz_y_y_xy, g_xz_y_y_xz, g_xz_y_y_yy, g_xz_y_y_yz, g_xz_y_y_zz, g_z_0_y_0_x_y_0_xx, g_z_0_y_0_x_y_0_xy, g_z_0_y_0_x_y_0_xz, g_z_0_y_0_x_y_0_yy, g_z_0_y_0_x_y_0_yz, g_z_0_y_0_x_y_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_y_0_xx[i] = 4.0 * g_xz_y_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_0_xy[i] = 4.0 * g_xz_y_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_0_xz[i] = 4.0 * g_xz_y_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_0_yy[i] = 4.0 * g_xz_y_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_0_yz[i] = 4.0 * g_xz_y_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_y_0_zz[i] = 4.0 * g_xz_y_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_xz_z_y_xx, g_xz_z_y_xy, g_xz_z_y_xz, g_xz_z_y_yy, g_xz_z_y_yz, g_xz_z_y_zz, g_z_0_y_0_x_z_0_xx, g_z_0_y_0_x_z_0_xy, g_z_0_y_0_x_z_0_xz, g_z_0_y_0_x_z_0_yy, g_z_0_y_0_x_z_0_yz, g_z_0_y_0_x_z_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_x_z_0_xx[i] = 4.0 * g_xz_z_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_0_xy[i] = 4.0 * g_xz_z_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_0_xz[i] = 4.0 * g_xz_z_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_0_yy[i] = 4.0 * g_xz_z_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_0_yz[i] = 4.0 * g_xz_z_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_x_z_0_zz[i] = 4.0 * g_xz_z_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_yz_x_y_xx, g_yz_x_y_xy, g_yz_x_y_xz, g_yz_x_y_yy, g_yz_x_y_yz, g_yz_x_y_zz, g_z_0_y_0_y_x_0_xx, g_z_0_y_0_y_x_0_xy, g_z_0_y_0_y_x_0_xz, g_z_0_y_0_y_x_0_yy, g_z_0_y_0_y_x_0_yz, g_z_0_y_0_y_x_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_x_0_xx[i] = 4.0 * g_yz_x_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_0_xy[i] = 4.0 * g_yz_x_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_0_xz[i] = 4.0 * g_yz_x_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_0_yy[i] = 4.0 * g_yz_x_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_0_yz[i] = 4.0 * g_yz_x_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_x_0_zz[i] = 4.0 * g_yz_x_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_yz_y_y_xx, g_yz_y_y_xy, g_yz_y_y_xz, g_yz_y_y_yy, g_yz_y_y_yz, g_yz_y_y_zz, g_z_0_y_0_y_y_0_xx, g_z_0_y_0_y_y_0_xy, g_z_0_y_0_y_y_0_xz, g_z_0_y_0_y_y_0_yy, g_z_0_y_0_y_y_0_yz, g_z_0_y_0_y_y_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_y_0_xx[i] = 4.0 * g_yz_y_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_0_xy[i] = 4.0 * g_yz_y_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_0_xz[i] = 4.0 * g_yz_y_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_0_yy[i] = 4.0 * g_yz_y_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_0_yz[i] = 4.0 * g_yz_y_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_y_0_zz[i] = 4.0 * g_yz_y_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_yz_z_y_xx, g_yz_z_y_xy, g_yz_z_y_xz, g_yz_z_y_yy, g_yz_z_y_yz, g_yz_z_y_zz, g_z_0_y_0_y_z_0_xx, g_z_0_y_0_y_z_0_xy, g_z_0_y_0_y_z_0_xz, g_z_0_y_0_y_z_0_yy, g_z_0_y_0_y_z_0_yz, g_z_0_y_0_y_z_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_y_z_0_xx[i] = 4.0 * g_yz_z_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_0_xy[i] = 4.0 * g_yz_z_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_0_xz[i] = 4.0 * g_yz_z_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_0_yy[i] = 4.0 * g_yz_z_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_0_yz[i] = 4.0 * g_yz_z_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_y_z_0_zz[i] = 4.0 * g_yz_z_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz, g_z_0_y_0_z_x_0_xx, g_z_0_y_0_z_x_0_xy, g_z_0_y_0_z_x_0_xz, g_z_0_y_0_z_x_0_yy, g_z_0_y_0_z_x_0_yz, g_z_0_y_0_z_x_0_zz, g_zz_x_y_xx, g_zz_x_y_xy, g_zz_x_y_xz, g_zz_x_y_yy, g_zz_x_y_yz, g_zz_x_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_x_0_xx[i] = -2.0 * g_0_x_y_xx[i] * c_exps[i] + 4.0 * g_zz_x_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_0_xy[i] = -2.0 * g_0_x_y_xy[i] * c_exps[i] + 4.0 * g_zz_x_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_0_xz[i] = -2.0 * g_0_x_y_xz[i] * c_exps[i] + 4.0 * g_zz_x_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_0_yy[i] = -2.0 * g_0_x_y_yy[i] * c_exps[i] + 4.0 * g_zz_x_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_0_yz[i] = -2.0 * g_0_x_y_yz[i] * c_exps[i] + 4.0 * g_zz_x_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_x_0_zz[i] = -2.0 * g_0_x_y_zz[i] * c_exps[i] + 4.0 * g_zz_x_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_z_0_y_0_z_y_0_xx, g_z_0_y_0_z_y_0_xy, g_z_0_y_0_z_y_0_xz, g_z_0_y_0_z_y_0_yy, g_z_0_y_0_z_y_0_yz, g_z_0_y_0_z_y_0_zz, g_zz_y_y_xx, g_zz_y_y_xy, g_zz_y_y_xz, g_zz_y_y_yy, g_zz_y_y_yz, g_zz_y_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_y_0_xx[i] = -2.0 * g_0_y_y_xx[i] * c_exps[i] + 4.0 * g_zz_y_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_0_xy[i] = -2.0 * g_0_y_y_xy[i] * c_exps[i] + 4.0 * g_zz_y_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_0_xz[i] = -2.0 * g_0_y_y_xz[i] * c_exps[i] + 4.0 * g_zz_y_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_0_yy[i] = -2.0 * g_0_y_y_yy[i] * c_exps[i] + 4.0 * g_zz_y_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_0_yz[i] = -2.0 * g_0_y_y_yz[i] * c_exps[i] + 4.0 * g_zz_y_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_y_0_zz[i] = -2.0 * g_0_y_y_zz[i] * c_exps[i] + 4.0 * g_zz_y_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz, g_z_0_y_0_z_z_0_xx, g_z_0_y_0_z_z_0_xy, g_z_0_y_0_z_z_0_xz, g_z_0_y_0_z_z_0_yy, g_z_0_y_0_z_z_0_yz, g_z_0_y_0_z_z_0_zz, g_zz_z_y_xx, g_zz_z_y_xy, g_zz_z_y_xz, g_zz_z_y_yy, g_zz_z_y_yz, g_zz_z_y_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_z_z_0_xx[i] = -2.0 * g_0_z_y_xx[i] * c_exps[i] + 4.0 * g_zz_z_y_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_0_xy[i] = -2.0 * g_0_z_y_xy[i] * c_exps[i] + 4.0 * g_zz_z_y_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_0_xz[i] = -2.0 * g_0_z_y_xz[i] * c_exps[i] + 4.0 * g_zz_z_y_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_0_yy[i] = -2.0 * g_0_z_y_yy[i] * c_exps[i] + 4.0 * g_zz_z_y_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_0_yz[i] = -2.0 * g_0_z_y_yz[i] * c_exps[i] + 4.0 * g_zz_z_y_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_z_z_0_zz[i] = -2.0 * g_0_z_y_zz[i] * c_exps[i] + 4.0 * g_zz_z_y_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_xz_x_z_xx, g_xz_x_z_xy, g_xz_x_z_xz, g_xz_x_z_yy, g_xz_x_z_yz, g_xz_x_z_zz, g_z_0_z_0_x_x_0_xx, g_z_0_z_0_x_x_0_xy, g_z_0_z_0_x_x_0_xz, g_z_0_z_0_x_x_0_yy, g_z_0_z_0_x_x_0_yz, g_z_0_z_0_x_x_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_x_0_xx[i] = 4.0 * g_xz_x_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_0_xy[i] = 4.0 * g_xz_x_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_0_xz[i] = 4.0 * g_xz_x_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_0_yy[i] = 4.0 * g_xz_x_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_0_yz[i] = 4.0 * g_xz_x_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_x_0_zz[i] = 4.0 * g_xz_x_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_xz_y_z_xx, g_xz_y_z_xy, g_xz_y_z_xz, g_xz_y_z_yy, g_xz_y_z_yz, g_xz_y_z_zz, g_z_0_z_0_x_y_0_xx, g_z_0_z_0_x_y_0_xy, g_z_0_z_0_x_y_0_xz, g_z_0_z_0_x_y_0_yy, g_z_0_z_0_x_y_0_yz, g_z_0_z_0_x_y_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_y_0_xx[i] = 4.0 * g_xz_y_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_0_xy[i] = 4.0 * g_xz_y_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_0_xz[i] = 4.0 * g_xz_y_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_0_yy[i] = 4.0 * g_xz_y_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_0_yz[i] = 4.0 * g_xz_y_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_y_0_zz[i] = 4.0 * g_xz_y_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_xz_z_z_xx, g_xz_z_z_xy, g_xz_z_z_xz, g_xz_z_z_yy, g_xz_z_z_yz, g_xz_z_z_zz, g_z_0_z_0_x_z_0_xx, g_z_0_z_0_x_z_0_xy, g_z_0_z_0_x_z_0_xz, g_z_0_z_0_x_z_0_yy, g_z_0_z_0_x_z_0_yz, g_z_0_z_0_x_z_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_x_z_0_xx[i] = 4.0 * g_xz_z_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_0_xy[i] = 4.0 * g_xz_z_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_0_xz[i] = 4.0 * g_xz_z_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_0_yy[i] = 4.0 * g_xz_z_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_0_yz[i] = 4.0 * g_xz_z_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_x_z_0_zz[i] = 4.0 * g_xz_z_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_yz_x_z_xx, g_yz_x_z_xy, g_yz_x_z_xz, g_yz_x_z_yy, g_yz_x_z_yz, g_yz_x_z_zz, g_z_0_z_0_y_x_0_xx, g_z_0_z_0_y_x_0_xy, g_z_0_z_0_y_x_0_xz, g_z_0_z_0_y_x_0_yy, g_z_0_z_0_y_x_0_yz, g_z_0_z_0_y_x_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_x_0_xx[i] = 4.0 * g_yz_x_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_0_xy[i] = 4.0 * g_yz_x_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_0_xz[i] = 4.0 * g_yz_x_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_0_yy[i] = 4.0 * g_yz_x_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_0_yz[i] = 4.0 * g_yz_x_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_x_0_zz[i] = 4.0 * g_yz_x_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_yz_y_z_xx, g_yz_y_z_xy, g_yz_y_z_xz, g_yz_y_z_yy, g_yz_y_z_yz, g_yz_y_z_zz, g_z_0_z_0_y_y_0_xx, g_z_0_z_0_y_y_0_xy, g_z_0_z_0_y_y_0_xz, g_z_0_z_0_y_y_0_yy, g_z_0_z_0_y_y_0_yz, g_z_0_z_0_y_y_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_y_0_xx[i] = 4.0 * g_yz_y_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_0_xy[i] = 4.0 * g_yz_y_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_0_xz[i] = 4.0 * g_yz_y_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_0_yy[i] = 4.0 * g_yz_y_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_0_yz[i] = 4.0 * g_yz_y_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_y_0_zz[i] = 4.0 * g_yz_y_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_yz_z_z_xx, g_yz_z_z_xy, g_yz_z_z_xz, g_yz_z_z_yy, g_yz_z_z_yz, g_yz_z_z_zz, g_z_0_z_0_y_z_0_xx, g_z_0_z_0_y_z_0_xy, g_z_0_z_0_y_z_0_xz, g_z_0_z_0_y_z_0_yy, g_z_0_z_0_y_z_0_yz, g_z_0_z_0_y_z_0_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_y_z_0_xx[i] = 4.0 * g_yz_z_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_0_xy[i] = 4.0 * g_yz_z_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_0_xz[i] = 4.0 * g_yz_z_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_0_yy[i] = 4.0 * g_yz_z_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_0_yz[i] = 4.0 * g_yz_z_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_y_z_0_zz[i] = 4.0 * g_yz_z_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz, g_z_0_z_0_z_x_0_xx, g_z_0_z_0_z_x_0_xy, g_z_0_z_0_z_x_0_xz, g_z_0_z_0_z_x_0_yy, g_z_0_z_0_z_x_0_yz, g_z_0_z_0_z_x_0_zz, g_zz_x_z_xx, g_zz_x_z_xy, g_zz_x_z_xz, g_zz_x_z_yy, g_zz_x_z_yz, g_zz_x_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_x_0_xx[i] = -2.0 * g_0_x_z_xx[i] * c_exps[i] + 4.0 * g_zz_x_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_0_xy[i] = -2.0 * g_0_x_z_xy[i] * c_exps[i] + 4.0 * g_zz_x_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_0_xz[i] = -2.0 * g_0_x_z_xz[i] * c_exps[i] + 4.0 * g_zz_x_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_0_yy[i] = -2.0 * g_0_x_z_yy[i] * c_exps[i] + 4.0 * g_zz_x_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_0_yz[i] = -2.0 * g_0_x_z_yz[i] * c_exps[i] + 4.0 * g_zz_x_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_x_0_zz[i] = -2.0 * g_0_x_z_zz[i] * c_exps[i] + 4.0 * g_zz_x_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz, g_z_0_z_0_z_y_0_xx, g_z_0_z_0_z_y_0_xy, g_z_0_z_0_z_y_0_xz, g_z_0_z_0_z_y_0_yy, g_z_0_z_0_z_y_0_yz, g_z_0_z_0_z_y_0_zz, g_zz_y_z_xx, g_zz_y_z_xy, g_zz_y_z_xz, g_zz_y_z_yy, g_zz_y_z_yz, g_zz_y_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_y_0_xx[i] = -2.0 * g_0_y_z_xx[i] * c_exps[i] + 4.0 * g_zz_y_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_0_xy[i] = -2.0 * g_0_y_z_xy[i] * c_exps[i] + 4.0 * g_zz_y_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_0_xz[i] = -2.0 * g_0_y_z_xz[i] * c_exps[i] + 4.0 * g_zz_y_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_0_yy[i] = -2.0 * g_0_y_z_yy[i] * c_exps[i] + 4.0 * g_zz_y_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_0_yz[i] = -2.0 * g_0_y_z_yz[i] * c_exps[i] + 4.0 * g_zz_y_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_y_0_zz[i] = -2.0 * g_0_y_z_zz[i] * c_exps[i] + 4.0 * g_zz_y_z_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_z_0_z_0_z_z_0_xx, g_z_0_z_0_z_z_0_xy, g_z_0_z_0_z_z_0_xz, g_z_0_z_0_z_z_0_yy, g_z_0_z_0_z_z_0_yz, g_z_0_z_0_z_z_0_zz, g_zz_z_z_xx, g_zz_z_z_xy, g_zz_z_z_xz, g_zz_z_z_yy, g_zz_z_z_yz, g_zz_z_z_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_z_z_0_xx[i] = -2.0 * g_0_z_z_xx[i] * c_exps[i] + 4.0 * g_zz_z_z_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_0_xy[i] = -2.0 * g_0_z_z_xy[i] * c_exps[i] + 4.0 * g_zz_z_z_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_0_xz[i] = -2.0 * g_0_z_z_xz[i] * c_exps[i] + 4.0 * g_zz_z_z_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_0_yy[i] = -2.0 * g_0_z_z_yy[i] * c_exps[i] + 4.0 * g_zz_z_z_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_0_yz[i] = -2.0 * g_0_z_z_yz[i] * c_exps[i] + 4.0 * g_zz_z_z_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_z_z_0_zz[i] = -2.0 * g_0_z_z_zz[i] * c_exps[i] + 4.0 * g_zz_z_z_zz[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace

