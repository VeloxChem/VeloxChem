#include "GeomDeriv1100OfScalarForSPPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_sppd_0(CSimdArray<double>& buffer_1100_sppd,
                     const CSimdArray<double>& buffer_pspd,
                     const CSimdArray<double>& buffer_pdpd,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_sppd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pspd

    auto g_x_0_x_xx = buffer_pspd[0];

    auto g_x_0_x_xy = buffer_pspd[1];

    auto g_x_0_x_xz = buffer_pspd[2];

    auto g_x_0_x_yy = buffer_pspd[3];

    auto g_x_0_x_yz = buffer_pspd[4];

    auto g_x_0_x_zz = buffer_pspd[5];

    auto g_x_0_y_xx = buffer_pspd[6];

    auto g_x_0_y_xy = buffer_pspd[7];

    auto g_x_0_y_xz = buffer_pspd[8];

    auto g_x_0_y_yy = buffer_pspd[9];

    auto g_x_0_y_yz = buffer_pspd[10];

    auto g_x_0_y_zz = buffer_pspd[11];

    auto g_x_0_z_xx = buffer_pspd[12];

    auto g_x_0_z_xy = buffer_pspd[13];

    auto g_x_0_z_xz = buffer_pspd[14];

    auto g_x_0_z_yy = buffer_pspd[15];

    auto g_x_0_z_yz = buffer_pspd[16];

    auto g_x_0_z_zz = buffer_pspd[17];

    auto g_y_0_x_xx = buffer_pspd[18];

    auto g_y_0_x_xy = buffer_pspd[19];

    auto g_y_0_x_xz = buffer_pspd[20];

    auto g_y_0_x_yy = buffer_pspd[21];

    auto g_y_0_x_yz = buffer_pspd[22];

    auto g_y_0_x_zz = buffer_pspd[23];

    auto g_y_0_y_xx = buffer_pspd[24];

    auto g_y_0_y_xy = buffer_pspd[25];

    auto g_y_0_y_xz = buffer_pspd[26];

    auto g_y_0_y_yy = buffer_pspd[27];

    auto g_y_0_y_yz = buffer_pspd[28];

    auto g_y_0_y_zz = buffer_pspd[29];

    auto g_y_0_z_xx = buffer_pspd[30];

    auto g_y_0_z_xy = buffer_pspd[31];

    auto g_y_0_z_xz = buffer_pspd[32];

    auto g_y_0_z_yy = buffer_pspd[33];

    auto g_y_0_z_yz = buffer_pspd[34];

    auto g_y_0_z_zz = buffer_pspd[35];

    auto g_z_0_x_xx = buffer_pspd[36];

    auto g_z_0_x_xy = buffer_pspd[37];

    auto g_z_0_x_xz = buffer_pspd[38];

    auto g_z_0_x_yy = buffer_pspd[39];

    auto g_z_0_x_yz = buffer_pspd[40];

    auto g_z_0_x_zz = buffer_pspd[41];

    auto g_z_0_y_xx = buffer_pspd[42];

    auto g_z_0_y_xy = buffer_pspd[43];

    auto g_z_0_y_xz = buffer_pspd[44];

    auto g_z_0_y_yy = buffer_pspd[45];

    auto g_z_0_y_yz = buffer_pspd[46];

    auto g_z_0_y_zz = buffer_pspd[47];

    auto g_z_0_z_xx = buffer_pspd[48];

    auto g_z_0_z_xy = buffer_pspd[49];

    auto g_z_0_z_xz = buffer_pspd[50];

    auto g_z_0_z_yy = buffer_pspd[51];

    auto g_z_0_z_yz = buffer_pspd[52];

    auto g_z_0_z_zz = buffer_pspd[53];

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

    /// Set up components of integrals buffer : buffer_1100_sppd

    auto g_x_x_0_0_0_x_x_xx = buffer_1100_sppd[0];

    auto g_x_x_0_0_0_x_x_xy = buffer_1100_sppd[1];

    auto g_x_x_0_0_0_x_x_xz = buffer_1100_sppd[2];

    auto g_x_x_0_0_0_x_x_yy = buffer_1100_sppd[3];

    auto g_x_x_0_0_0_x_x_yz = buffer_1100_sppd[4];

    auto g_x_x_0_0_0_x_x_zz = buffer_1100_sppd[5];

    auto g_x_x_0_0_0_x_y_xx = buffer_1100_sppd[6];

    auto g_x_x_0_0_0_x_y_xy = buffer_1100_sppd[7];

    auto g_x_x_0_0_0_x_y_xz = buffer_1100_sppd[8];

    auto g_x_x_0_0_0_x_y_yy = buffer_1100_sppd[9];

    auto g_x_x_0_0_0_x_y_yz = buffer_1100_sppd[10];

    auto g_x_x_0_0_0_x_y_zz = buffer_1100_sppd[11];

    auto g_x_x_0_0_0_x_z_xx = buffer_1100_sppd[12];

    auto g_x_x_0_0_0_x_z_xy = buffer_1100_sppd[13];

    auto g_x_x_0_0_0_x_z_xz = buffer_1100_sppd[14];

    auto g_x_x_0_0_0_x_z_yy = buffer_1100_sppd[15];

    auto g_x_x_0_0_0_x_z_yz = buffer_1100_sppd[16];

    auto g_x_x_0_0_0_x_z_zz = buffer_1100_sppd[17];

    auto g_x_x_0_0_0_y_x_xx = buffer_1100_sppd[18];

    auto g_x_x_0_0_0_y_x_xy = buffer_1100_sppd[19];

    auto g_x_x_0_0_0_y_x_xz = buffer_1100_sppd[20];

    auto g_x_x_0_0_0_y_x_yy = buffer_1100_sppd[21];

    auto g_x_x_0_0_0_y_x_yz = buffer_1100_sppd[22];

    auto g_x_x_0_0_0_y_x_zz = buffer_1100_sppd[23];

    auto g_x_x_0_0_0_y_y_xx = buffer_1100_sppd[24];

    auto g_x_x_0_0_0_y_y_xy = buffer_1100_sppd[25];

    auto g_x_x_0_0_0_y_y_xz = buffer_1100_sppd[26];

    auto g_x_x_0_0_0_y_y_yy = buffer_1100_sppd[27];

    auto g_x_x_0_0_0_y_y_yz = buffer_1100_sppd[28];

    auto g_x_x_0_0_0_y_y_zz = buffer_1100_sppd[29];

    auto g_x_x_0_0_0_y_z_xx = buffer_1100_sppd[30];

    auto g_x_x_0_0_0_y_z_xy = buffer_1100_sppd[31];

    auto g_x_x_0_0_0_y_z_xz = buffer_1100_sppd[32];

    auto g_x_x_0_0_0_y_z_yy = buffer_1100_sppd[33];

    auto g_x_x_0_0_0_y_z_yz = buffer_1100_sppd[34];

    auto g_x_x_0_0_0_y_z_zz = buffer_1100_sppd[35];

    auto g_x_x_0_0_0_z_x_xx = buffer_1100_sppd[36];

    auto g_x_x_0_0_0_z_x_xy = buffer_1100_sppd[37];

    auto g_x_x_0_0_0_z_x_xz = buffer_1100_sppd[38];

    auto g_x_x_0_0_0_z_x_yy = buffer_1100_sppd[39];

    auto g_x_x_0_0_0_z_x_yz = buffer_1100_sppd[40];

    auto g_x_x_0_0_0_z_x_zz = buffer_1100_sppd[41];

    auto g_x_x_0_0_0_z_y_xx = buffer_1100_sppd[42];

    auto g_x_x_0_0_0_z_y_xy = buffer_1100_sppd[43];

    auto g_x_x_0_0_0_z_y_xz = buffer_1100_sppd[44];

    auto g_x_x_0_0_0_z_y_yy = buffer_1100_sppd[45];

    auto g_x_x_0_0_0_z_y_yz = buffer_1100_sppd[46];

    auto g_x_x_0_0_0_z_y_zz = buffer_1100_sppd[47];

    auto g_x_x_0_0_0_z_z_xx = buffer_1100_sppd[48];

    auto g_x_x_0_0_0_z_z_xy = buffer_1100_sppd[49];

    auto g_x_x_0_0_0_z_z_xz = buffer_1100_sppd[50];

    auto g_x_x_0_0_0_z_z_yy = buffer_1100_sppd[51];

    auto g_x_x_0_0_0_z_z_yz = buffer_1100_sppd[52];

    auto g_x_x_0_0_0_z_z_zz = buffer_1100_sppd[53];

    auto g_x_y_0_0_0_x_x_xx = buffer_1100_sppd[54];

    auto g_x_y_0_0_0_x_x_xy = buffer_1100_sppd[55];

    auto g_x_y_0_0_0_x_x_xz = buffer_1100_sppd[56];

    auto g_x_y_0_0_0_x_x_yy = buffer_1100_sppd[57];

    auto g_x_y_0_0_0_x_x_yz = buffer_1100_sppd[58];

    auto g_x_y_0_0_0_x_x_zz = buffer_1100_sppd[59];

    auto g_x_y_0_0_0_x_y_xx = buffer_1100_sppd[60];

    auto g_x_y_0_0_0_x_y_xy = buffer_1100_sppd[61];

    auto g_x_y_0_0_0_x_y_xz = buffer_1100_sppd[62];

    auto g_x_y_0_0_0_x_y_yy = buffer_1100_sppd[63];

    auto g_x_y_0_0_0_x_y_yz = buffer_1100_sppd[64];

    auto g_x_y_0_0_0_x_y_zz = buffer_1100_sppd[65];

    auto g_x_y_0_0_0_x_z_xx = buffer_1100_sppd[66];

    auto g_x_y_0_0_0_x_z_xy = buffer_1100_sppd[67];

    auto g_x_y_0_0_0_x_z_xz = buffer_1100_sppd[68];

    auto g_x_y_0_0_0_x_z_yy = buffer_1100_sppd[69];

    auto g_x_y_0_0_0_x_z_yz = buffer_1100_sppd[70];

    auto g_x_y_0_0_0_x_z_zz = buffer_1100_sppd[71];

    auto g_x_y_0_0_0_y_x_xx = buffer_1100_sppd[72];

    auto g_x_y_0_0_0_y_x_xy = buffer_1100_sppd[73];

    auto g_x_y_0_0_0_y_x_xz = buffer_1100_sppd[74];

    auto g_x_y_0_0_0_y_x_yy = buffer_1100_sppd[75];

    auto g_x_y_0_0_0_y_x_yz = buffer_1100_sppd[76];

    auto g_x_y_0_0_0_y_x_zz = buffer_1100_sppd[77];

    auto g_x_y_0_0_0_y_y_xx = buffer_1100_sppd[78];

    auto g_x_y_0_0_0_y_y_xy = buffer_1100_sppd[79];

    auto g_x_y_0_0_0_y_y_xz = buffer_1100_sppd[80];

    auto g_x_y_0_0_0_y_y_yy = buffer_1100_sppd[81];

    auto g_x_y_0_0_0_y_y_yz = buffer_1100_sppd[82];

    auto g_x_y_0_0_0_y_y_zz = buffer_1100_sppd[83];

    auto g_x_y_0_0_0_y_z_xx = buffer_1100_sppd[84];

    auto g_x_y_0_0_0_y_z_xy = buffer_1100_sppd[85];

    auto g_x_y_0_0_0_y_z_xz = buffer_1100_sppd[86];

    auto g_x_y_0_0_0_y_z_yy = buffer_1100_sppd[87];

    auto g_x_y_0_0_0_y_z_yz = buffer_1100_sppd[88];

    auto g_x_y_0_0_0_y_z_zz = buffer_1100_sppd[89];

    auto g_x_y_0_0_0_z_x_xx = buffer_1100_sppd[90];

    auto g_x_y_0_0_0_z_x_xy = buffer_1100_sppd[91];

    auto g_x_y_0_0_0_z_x_xz = buffer_1100_sppd[92];

    auto g_x_y_0_0_0_z_x_yy = buffer_1100_sppd[93];

    auto g_x_y_0_0_0_z_x_yz = buffer_1100_sppd[94];

    auto g_x_y_0_0_0_z_x_zz = buffer_1100_sppd[95];

    auto g_x_y_0_0_0_z_y_xx = buffer_1100_sppd[96];

    auto g_x_y_0_0_0_z_y_xy = buffer_1100_sppd[97];

    auto g_x_y_0_0_0_z_y_xz = buffer_1100_sppd[98];

    auto g_x_y_0_0_0_z_y_yy = buffer_1100_sppd[99];

    auto g_x_y_0_0_0_z_y_yz = buffer_1100_sppd[100];

    auto g_x_y_0_0_0_z_y_zz = buffer_1100_sppd[101];

    auto g_x_y_0_0_0_z_z_xx = buffer_1100_sppd[102];

    auto g_x_y_0_0_0_z_z_xy = buffer_1100_sppd[103];

    auto g_x_y_0_0_0_z_z_xz = buffer_1100_sppd[104];

    auto g_x_y_0_0_0_z_z_yy = buffer_1100_sppd[105];

    auto g_x_y_0_0_0_z_z_yz = buffer_1100_sppd[106];

    auto g_x_y_0_0_0_z_z_zz = buffer_1100_sppd[107];

    auto g_x_z_0_0_0_x_x_xx = buffer_1100_sppd[108];

    auto g_x_z_0_0_0_x_x_xy = buffer_1100_sppd[109];

    auto g_x_z_0_0_0_x_x_xz = buffer_1100_sppd[110];

    auto g_x_z_0_0_0_x_x_yy = buffer_1100_sppd[111];

    auto g_x_z_0_0_0_x_x_yz = buffer_1100_sppd[112];

    auto g_x_z_0_0_0_x_x_zz = buffer_1100_sppd[113];

    auto g_x_z_0_0_0_x_y_xx = buffer_1100_sppd[114];

    auto g_x_z_0_0_0_x_y_xy = buffer_1100_sppd[115];

    auto g_x_z_0_0_0_x_y_xz = buffer_1100_sppd[116];

    auto g_x_z_0_0_0_x_y_yy = buffer_1100_sppd[117];

    auto g_x_z_0_0_0_x_y_yz = buffer_1100_sppd[118];

    auto g_x_z_0_0_0_x_y_zz = buffer_1100_sppd[119];

    auto g_x_z_0_0_0_x_z_xx = buffer_1100_sppd[120];

    auto g_x_z_0_0_0_x_z_xy = buffer_1100_sppd[121];

    auto g_x_z_0_0_0_x_z_xz = buffer_1100_sppd[122];

    auto g_x_z_0_0_0_x_z_yy = buffer_1100_sppd[123];

    auto g_x_z_0_0_0_x_z_yz = buffer_1100_sppd[124];

    auto g_x_z_0_0_0_x_z_zz = buffer_1100_sppd[125];

    auto g_x_z_0_0_0_y_x_xx = buffer_1100_sppd[126];

    auto g_x_z_0_0_0_y_x_xy = buffer_1100_sppd[127];

    auto g_x_z_0_0_0_y_x_xz = buffer_1100_sppd[128];

    auto g_x_z_0_0_0_y_x_yy = buffer_1100_sppd[129];

    auto g_x_z_0_0_0_y_x_yz = buffer_1100_sppd[130];

    auto g_x_z_0_0_0_y_x_zz = buffer_1100_sppd[131];

    auto g_x_z_0_0_0_y_y_xx = buffer_1100_sppd[132];

    auto g_x_z_0_0_0_y_y_xy = buffer_1100_sppd[133];

    auto g_x_z_0_0_0_y_y_xz = buffer_1100_sppd[134];

    auto g_x_z_0_0_0_y_y_yy = buffer_1100_sppd[135];

    auto g_x_z_0_0_0_y_y_yz = buffer_1100_sppd[136];

    auto g_x_z_0_0_0_y_y_zz = buffer_1100_sppd[137];

    auto g_x_z_0_0_0_y_z_xx = buffer_1100_sppd[138];

    auto g_x_z_0_0_0_y_z_xy = buffer_1100_sppd[139];

    auto g_x_z_0_0_0_y_z_xz = buffer_1100_sppd[140];

    auto g_x_z_0_0_0_y_z_yy = buffer_1100_sppd[141];

    auto g_x_z_0_0_0_y_z_yz = buffer_1100_sppd[142];

    auto g_x_z_0_0_0_y_z_zz = buffer_1100_sppd[143];

    auto g_x_z_0_0_0_z_x_xx = buffer_1100_sppd[144];

    auto g_x_z_0_0_0_z_x_xy = buffer_1100_sppd[145];

    auto g_x_z_0_0_0_z_x_xz = buffer_1100_sppd[146];

    auto g_x_z_0_0_0_z_x_yy = buffer_1100_sppd[147];

    auto g_x_z_0_0_0_z_x_yz = buffer_1100_sppd[148];

    auto g_x_z_0_0_0_z_x_zz = buffer_1100_sppd[149];

    auto g_x_z_0_0_0_z_y_xx = buffer_1100_sppd[150];

    auto g_x_z_0_0_0_z_y_xy = buffer_1100_sppd[151];

    auto g_x_z_0_0_0_z_y_xz = buffer_1100_sppd[152];

    auto g_x_z_0_0_0_z_y_yy = buffer_1100_sppd[153];

    auto g_x_z_0_0_0_z_y_yz = buffer_1100_sppd[154];

    auto g_x_z_0_0_0_z_y_zz = buffer_1100_sppd[155];

    auto g_x_z_0_0_0_z_z_xx = buffer_1100_sppd[156];

    auto g_x_z_0_0_0_z_z_xy = buffer_1100_sppd[157];

    auto g_x_z_0_0_0_z_z_xz = buffer_1100_sppd[158];

    auto g_x_z_0_0_0_z_z_yy = buffer_1100_sppd[159];

    auto g_x_z_0_0_0_z_z_yz = buffer_1100_sppd[160];

    auto g_x_z_0_0_0_z_z_zz = buffer_1100_sppd[161];

    auto g_y_x_0_0_0_x_x_xx = buffer_1100_sppd[162];

    auto g_y_x_0_0_0_x_x_xy = buffer_1100_sppd[163];

    auto g_y_x_0_0_0_x_x_xz = buffer_1100_sppd[164];

    auto g_y_x_0_0_0_x_x_yy = buffer_1100_sppd[165];

    auto g_y_x_0_0_0_x_x_yz = buffer_1100_sppd[166];

    auto g_y_x_0_0_0_x_x_zz = buffer_1100_sppd[167];

    auto g_y_x_0_0_0_x_y_xx = buffer_1100_sppd[168];

    auto g_y_x_0_0_0_x_y_xy = buffer_1100_sppd[169];

    auto g_y_x_0_0_0_x_y_xz = buffer_1100_sppd[170];

    auto g_y_x_0_0_0_x_y_yy = buffer_1100_sppd[171];

    auto g_y_x_0_0_0_x_y_yz = buffer_1100_sppd[172];

    auto g_y_x_0_0_0_x_y_zz = buffer_1100_sppd[173];

    auto g_y_x_0_0_0_x_z_xx = buffer_1100_sppd[174];

    auto g_y_x_0_0_0_x_z_xy = buffer_1100_sppd[175];

    auto g_y_x_0_0_0_x_z_xz = buffer_1100_sppd[176];

    auto g_y_x_0_0_0_x_z_yy = buffer_1100_sppd[177];

    auto g_y_x_0_0_0_x_z_yz = buffer_1100_sppd[178];

    auto g_y_x_0_0_0_x_z_zz = buffer_1100_sppd[179];

    auto g_y_x_0_0_0_y_x_xx = buffer_1100_sppd[180];

    auto g_y_x_0_0_0_y_x_xy = buffer_1100_sppd[181];

    auto g_y_x_0_0_0_y_x_xz = buffer_1100_sppd[182];

    auto g_y_x_0_0_0_y_x_yy = buffer_1100_sppd[183];

    auto g_y_x_0_0_0_y_x_yz = buffer_1100_sppd[184];

    auto g_y_x_0_0_0_y_x_zz = buffer_1100_sppd[185];

    auto g_y_x_0_0_0_y_y_xx = buffer_1100_sppd[186];

    auto g_y_x_0_0_0_y_y_xy = buffer_1100_sppd[187];

    auto g_y_x_0_0_0_y_y_xz = buffer_1100_sppd[188];

    auto g_y_x_0_0_0_y_y_yy = buffer_1100_sppd[189];

    auto g_y_x_0_0_0_y_y_yz = buffer_1100_sppd[190];

    auto g_y_x_0_0_0_y_y_zz = buffer_1100_sppd[191];

    auto g_y_x_0_0_0_y_z_xx = buffer_1100_sppd[192];

    auto g_y_x_0_0_0_y_z_xy = buffer_1100_sppd[193];

    auto g_y_x_0_0_0_y_z_xz = buffer_1100_sppd[194];

    auto g_y_x_0_0_0_y_z_yy = buffer_1100_sppd[195];

    auto g_y_x_0_0_0_y_z_yz = buffer_1100_sppd[196];

    auto g_y_x_0_0_0_y_z_zz = buffer_1100_sppd[197];

    auto g_y_x_0_0_0_z_x_xx = buffer_1100_sppd[198];

    auto g_y_x_0_0_0_z_x_xy = buffer_1100_sppd[199];

    auto g_y_x_0_0_0_z_x_xz = buffer_1100_sppd[200];

    auto g_y_x_0_0_0_z_x_yy = buffer_1100_sppd[201];

    auto g_y_x_0_0_0_z_x_yz = buffer_1100_sppd[202];

    auto g_y_x_0_0_0_z_x_zz = buffer_1100_sppd[203];

    auto g_y_x_0_0_0_z_y_xx = buffer_1100_sppd[204];

    auto g_y_x_0_0_0_z_y_xy = buffer_1100_sppd[205];

    auto g_y_x_0_0_0_z_y_xz = buffer_1100_sppd[206];

    auto g_y_x_0_0_0_z_y_yy = buffer_1100_sppd[207];

    auto g_y_x_0_0_0_z_y_yz = buffer_1100_sppd[208];

    auto g_y_x_0_0_0_z_y_zz = buffer_1100_sppd[209];

    auto g_y_x_0_0_0_z_z_xx = buffer_1100_sppd[210];

    auto g_y_x_0_0_0_z_z_xy = buffer_1100_sppd[211];

    auto g_y_x_0_0_0_z_z_xz = buffer_1100_sppd[212];

    auto g_y_x_0_0_0_z_z_yy = buffer_1100_sppd[213];

    auto g_y_x_0_0_0_z_z_yz = buffer_1100_sppd[214];

    auto g_y_x_0_0_0_z_z_zz = buffer_1100_sppd[215];

    auto g_y_y_0_0_0_x_x_xx = buffer_1100_sppd[216];

    auto g_y_y_0_0_0_x_x_xy = buffer_1100_sppd[217];

    auto g_y_y_0_0_0_x_x_xz = buffer_1100_sppd[218];

    auto g_y_y_0_0_0_x_x_yy = buffer_1100_sppd[219];

    auto g_y_y_0_0_0_x_x_yz = buffer_1100_sppd[220];

    auto g_y_y_0_0_0_x_x_zz = buffer_1100_sppd[221];

    auto g_y_y_0_0_0_x_y_xx = buffer_1100_sppd[222];

    auto g_y_y_0_0_0_x_y_xy = buffer_1100_sppd[223];

    auto g_y_y_0_0_0_x_y_xz = buffer_1100_sppd[224];

    auto g_y_y_0_0_0_x_y_yy = buffer_1100_sppd[225];

    auto g_y_y_0_0_0_x_y_yz = buffer_1100_sppd[226];

    auto g_y_y_0_0_0_x_y_zz = buffer_1100_sppd[227];

    auto g_y_y_0_0_0_x_z_xx = buffer_1100_sppd[228];

    auto g_y_y_0_0_0_x_z_xy = buffer_1100_sppd[229];

    auto g_y_y_0_0_0_x_z_xz = buffer_1100_sppd[230];

    auto g_y_y_0_0_0_x_z_yy = buffer_1100_sppd[231];

    auto g_y_y_0_0_0_x_z_yz = buffer_1100_sppd[232];

    auto g_y_y_0_0_0_x_z_zz = buffer_1100_sppd[233];

    auto g_y_y_0_0_0_y_x_xx = buffer_1100_sppd[234];

    auto g_y_y_0_0_0_y_x_xy = buffer_1100_sppd[235];

    auto g_y_y_0_0_0_y_x_xz = buffer_1100_sppd[236];

    auto g_y_y_0_0_0_y_x_yy = buffer_1100_sppd[237];

    auto g_y_y_0_0_0_y_x_yz = buffer_1100_sppd[238];

    auto g_y_y_0_0_0_y_x_zz = buffer_1100_sppd[239];

    auto g_y_y_0_0_0_y_y_xx = buffer_1100_sppd[240];

    auto g_y_y_0_0_0_y_y_xy = buffer_1100_sppd[241];

    auto g_y_y_0_0_0_y_y_xz = buffer_1100_sppd[242];

    auto g_y_y_0_0_0_y_y_yy = buffer_1100_sppd[243];

    auto g_y_y_0_0_0_y_y_yz = buffer_1100_sppd[244];

    auto g_y_y_0_0_0_y_y_zz = buffer_1100_sppd[245];

    auto g_y_y_0_0_0_y_z_xx = buffer_1100_sppd[246];

    auto g_y_y_0_0_0_y_z_xy = buffer_1100_sppd[247];

    auto g_y_y_0_0_0_y_z_xz = buffer_1100_sppd[248];

    auto g_y_y_0_0_0_y_z_yy = buffer_1100_sppd[249];

    auto g_y_y_0_0_0_y_z_yz = buffer_1100_sppd[250];

    auto g_y_y_0_0_0_y_z_zz = buffer_1100_sppd[251];

    auto g_y_y_0_0_0_z_x_xx = buffer_1100_sppd[252];

    auto g_y_y_0_0_0_z_x_xy = buffer_1100_sppd[253];

    auto g_y_y_0_0_0_z_x_xz = buffer_1100_sppd[254];

    auto g_y_y_0_0_0_z_x_yy = buffer_1100_sppd[255];

    auto g_y_y_0_0_0_z_x_yz = buffer_1100_sppd[256];

    auto g_y_y_0_0_0_z_x_zz = buffer_1100_sppd[257];

    auto g_y_y_0_0_0_z_y_xx = buffer_1100_sppd[258];

    auto g_y_y_0_0_0_z_y_xy = buffer_1100_sppd[259];

    auto g_y_y_0_0_0_z_y_xz = buffer_1100_sppd[260];

    auto g_y_y_0_0_0_z_y_yy = buffer_1100_sppd[261];

    auto g_y_y_0_0_0_z_y_yz = buffer_1100_sppd[262];

    auto g_y_y_0_0_0_z_y_zz = buffer_1100_sppd[263];

    auto g_y_y_0_0_0_z_z_xx = buffer_1100_sppd[264];

    auto g_y_y_0_0_0_z_z_xy = buffer_1100_sppd[265];

    auto g_y_y_0_0_0_z_z_xz = buffer_1100_sppd[266];

    auto g_y_y_0_0_0_z_z_yy = buffer_1100_sppd[267];

    auto g_y_y_0_0_0_z_z_yz = buffer_1100_sppd[268];

    auto g_y_y_0_0_0_z_z_zz = buffer_1100_sppd[269];

    auto g_y_z_0_0_0_x_x_xx = buffer_1100_sppd[270];

    auto g_y_z_0_0_0_x_x_xy = buffer_1100_sppd[271];

    auto g_y_z_0_0_0_x_x_xz = buffer_1100_sppd[272];

    auto g_y_z_0_0_0_x_x_yy = buffer_1100_sppd[273];

    auto g_y_z_0_0_0_x_x_yz = buffer_1100_sppd[274];

    auto g_y_z_0_0_0_x_x_zz = buffer_1100_sppd[275];

    auto g_y_z_0_0_0_x_y_xx = buffer_1100_sppd[276];

    auto g_y_z_0_0_0_x_y_xy = buffer_1100_sppd[277];

    auto g_y_z_0_0_0_x_y_xz = buffer_1100_sppd[278];

    auto g_y_z_0_0_0_x_y_yy = buffer_1100_sppd[279];

    auto g_y_z_0_0_0_x_y_yz = buffer_1100_sppd[280];

    auto g_y_z_0_0_0_x_y_zz = buffer_1100_sppd[281];

    auto g_y_z_0_0_0_x_z_xx = buffer_1100_sppd[282];

    auto g_y_z_0_0_0_x_z_xy = buffer_1100_sppd[283];

    auto g_y_z_0_0_0_x_z_xz = buffer_1100_sppd[284];

    auto g_y_z_0_0_0_x_z_yy = buffer_1100_sppd[285];

    auto g_y_z_0_0_0_x_z_yz = buffer_1100_sppd[286];

    auto g_y_z_0_0_0_x_z_zz = buffer_1100_sppd[287];

    auto g_y_z_0_0_0_y_x_xx = buffer_1100_sppd[288];

    auto g_y_z_0_0_0_y_x_xy = buffer_1100_sppd[289];

    auto g_y_z_0_0_0_y_x_xz = buffer_1100_sppd[290];

    auto g_y_z_0_0_0_y_x_yy = buffer_1100_sppd[291];

    auto g_y_z_0_0_0_y_x_yz = buffer_1100_sppd[292];

    auto g_y_z_0_0_0_y_x_zz = buffer_1100_sppd[293];

    auto g_y_z_0_0_0_y_y_xx = buffer_1100_sppd[294];

    auto g_y_z_0_0_0_y_y_xy = buffer_1100_sppd[295];

    auto g_y_z_0_0_0_y_y_xz = buffer_1100_sppd[296];

    auto g_y_z_0_0_0_y_y_yy = buffer_1100_sppd[297];

    auto g_y_z_0_0_0_y_y_yz = buffer_1100_sppd[298];

    auto g_y_z_0_0_0_y_y_zz = buffer_1100_sppd[299];

    auto g_y_z_0_0_0_y_z_xx = buffer_1100_sppd[300];

    auto g_y_z_0_0_0_y_z_xy = buffer_1100_sppd[301];

    auto g_y_z_0_0_0_y_z_xz = buffer_1100_sppd[302];

    auto g_y_z_0_0_0_y_z_yy = buffer_1100_sppd[303];

    auto g_y_z_0_0_0_y_z_yz = buffer_1100_sppd[304];

    auto g_y_z_0_0_0_y_z_zz = buffer_1100_sppd[305];

    auto g_y_z_0_0_0_z_x_xx = buffer_1100_sppd[306];

    auto g_y_z_0_0_0_z_x_xy = buffer_1100_sppd[307];

    auto g_y_z_0_0_0_z_x_xz = buffer_1100_sppd[308];

    auto g_y_z_0_0_0_z_x_yy = buffer_1100_sppd[309];

    auto g_y_z_0_0_0_z_x_yz = buffer_1100_sppd[310];

    auto g_y_z_0_0_0_z_x_zz = buffer_1100_sppd[311];

    auto g_y_z_0_0_0_z_y_xx = buffer_1100_sppd[312];

    auto g_y_z_0_0_0_z_y_xy = buffer_1100_sppd[313];

    auto g_y_z_0_0_0_z_y_xz = buffer_1100_sppd[314];

    auto g_y_z_0_0_0_z_y_yy = buffer_1100_sppd[315];

    auto g_y_z_0_0_0_z_y_yz = buffer_1100_sppd[316];

    auto g_y_z_0_0_0_z_y_zz = buffer_1100_sppd[317];

    auto g_y_z_0_0_0_z_z_xx = buffer_1100_sppd[318];

    auto g_y_z_0_0_0_z_z_xy = buffer_1100_sppd[319];

    auto g_y_z_0_0_0_z_z_xz = buffer_1100_sppd[320];

    auto g_y_z_0_0_0_z_z_yy = buffer_1100_sppd[321];

    auto g_y_z_0_0_0_z_z_yz = buffer_1100_sppd[322];

    auto g_y_z_0_0_0_z_z_zz = buffer_1100_sppd[323];

    auto g_z_x_0_0_0_x_x_xx = buffer_1100_sppd[324];

    auto g_z_x_0_0_0_x_x_xy = buffer_1100_sppd[325];

    auto g_z_x_0_0_0_x_x_xz = buffer_1100_sppd[326];

    auto g_z_x_0_0_0_x_x_yy = buffer_1100_sppd[327];

    auto g_z_x_0_0_0_x_x_yz = buffer_1100_sppd[328];

    auto g_z_x_0_0_0_x_x_zz = buffer_1100_sppd[329];

    auto g_z_x_0_0_0_x_y_xx = buffer_1100_sppd[330];

    auto g_z_x_0_0_0_x_y_xy = buffer_1100_sppd[331];

    auto g_z_x_0_0_0_x_y_xz = buffer_1100_sppd[332];

    auto g_z_x_0_0_0_x_y_yy = buffer_1100_sppd[333];

    auto g_z_x_0_0_0_x_y_yz = buffer_1100_sppd[334];

    auto g_z_x_0_0_0_x_y_zz = buffer_1100_sppd[335];

    auto g_z_x_0_0_0_x_z_xx = buffer_1100_sppd[336];

    auto g_z_x_0_0_0_x_z_xy = buffer_1100_sppd[337];

    auto g_z_x_0_0_0_x_z_xz = buffer_1100_sppd[338];

    auto g_z_x_0_0_0_x_z_yy = buffer_1100_sppd[339];

    auto g_z_x_0_0_0_x_z_yz = buffer_1100_sppd[340];

    auto g_z_x_0_0_0_x_z_zz = buffer_1100_sppd[341];

    auto g_z_x_0_0_0_y_x_xx = buffer_1100_sppd[342];

    auto g_z_x_0_0_0_y_x_xy = buffer_1100_sppd[343];

    auto g_z_x_0_0_0_y_x_xz = buffer_1100_sppd[344];

    auto g_z_x_0_0_0_y_x_yy = buffer_1100_sppd[345];

    auto g_z_x_0_0_0_y_x_yz = buffer_1100_sppd[346];

    auto g_z_x_0_0_0_y_x_zz = buffer_1100_sppd[347];

    auto g_z_x_0_0_0_y_y_xx = buffer_1100_sppd[348];

    auto g_z_x_0_0_0_y_y_xy = buffer_1100_sppd[349];

    auto g_z_x_0_0_0_y_y_xz = buffer_1100_sppd[350];

    auto g_z_x_0_0_0_y_y_yy = buffer_1100_sppd[351];

    auto g_z_x_0_0_0_y_y_yz = buffer_1100_sppd[352];

    auto g_z_x_0_0_0_y_y_zz = buffer_1100_sppd[353];

    auto g_z_x_0_0_0_y_z_xx = buffer_1100_sppd[354];

    auto g_z_x_0_0_0_y_z_xy = buffer_1100_sppd[355];

    auto g_z_x_0_0_0_y_z_xz = buffer_1100_sppd[356];

    auto g_z_x_0_0_0_y_z_yy = buffer_1100_sppd[357];

    auto g_z_x_0_0_0_y_z_yz = buffer_1100_sppd[358];

    auto g_z_x_0_0_0_y_z_zz = buffer_1100_sppd[359];

    auto g_z_x_0_0_0_z_x_xx = buffer_1100_sppd[360];

    auto g_z_x_0_0_0_z_x_xy = buffer_1100_sppd[361];

    auto g_z_x_0_0_0_z_x_xz = buffer_1100_sppd[362];

    auto g_z_x_0_0_0_z_x_yy = buffer_1100_sppd[363];

    auto g_z_x_0_0_0_z_x_yz = buffer_1100_sppd[364];

    auto g_z_x_0_0_0_z_x_zz = buffer_1100_sppd[365];

    auto g_z_x_0_0_0_z_y_xx = buffer_1100_sppd[366];

    auto g_z_x_0_0_0_z_y_xy = buffer_1100_sppd[367];

    auto g_z_x_0_0_0_z_y_xz = buffer_1100_sppd[368];

    auto g_z_x_0_0_0_z_y_yy = buffer_1100_sppd[369];

    auto g_z_x_0_0_0_z_y_yz = buffer_1100_sppd[370];

    auto g_z_x_0_0_0_z_y_zz = buffer_1100_sppd[371];

    auto g_z_x_0_0_0_z_z_xx = buffer_1100_sppd[372];

    auto g_z_x_0_0_0_z_z_xy = buffer_1100_sppd[373];

    auto g_z_x_0_0_0_z_z_xz = buffer_1100_sppd[374];

    auto g_z_x_0_0_0_z_z_yy = buffer_1100_sppd[375];

    auto g_z_x_0_0_0_z_z_yz = buffer_1100_sppd[376];

    auto g_z_x_0_0_0_z_z_zz = buffer_1100_sppd[377];

    auto g_z_y_0_0_0_x_x_xx = buffer_1100_sppd[378];

    auto g_z_y_0_0_0_x_x_xy = buffer_1100_sppd[379];

    auto g_z_y_0_0_0_x_x_xz = buffer_1100_sppd[380];

    auto g_z_y_0_0_0_x_x_yy = buffer_1100_sppd[381];

    auto g_z_y_0_0_0_x_x_yz = buffer_1100_sppd[382];

    auto g_z_y_0_0_0_x_x_zz = buffer_1100_sppd[383];

    auto g_z_y_0_0_0_x_y_xx = buffer_1100_sppd[384];

    auto g_z_y_0_0_0_x_y_xy = buffer_1100_sppd[385];

    auto g_z_y_0_0_0_x_y_xz = buffer_1100_sppd[386];

    auto g_z_y_0_0_0_x_y_yy = buffer_1100_sppd[387];

    auto g_z_y_0_0_0_x_y_yz = buffer_1100_sppd[388];

    auto g_z_y_0_0_0_x_y_zz = buffer_1100_sppd[389];

    auto g_z_y_0_0_0_x_z_xx = buffer_1100_sppd[390];

    auto g_z_y_0_0_0_x_z_xy = buffer_1100_sppd[391];

    auto g_z_y_0_0_0_x_z_xz = buffer_1100_sppd[392];

    auto g_z_y_0_0_0_x_z_yy = buffer_1100_sppd[393];

    auto g_z_y_0_0_0_x_z_yz = buffer_1100_sppd[394];

    auto g_z_y_0_0_0_x_z_zz = buffer_1100_sppd[395];

    auto g_z_y_0_0_0_y_x_xx = buffer_1100_sppd[396];

    auto g_z_y_0_0_0_y_x_xy = buffer_1100_sppd[397];

    auto g_z_y_0_0_0_y_x_xz = buffer_1100_sppd[398];

    auto g_z_y_0_0_0_y_x_yy = buffer_1100_sppd[399];

    auto g_z_y_0_0_0_y_x_yz = buffer_1100_sppd[400];

    auto g_z_y_0_0_0_y_x_zz = buffer_1100_sppd[401];

    auto g_z_y_0_0_0_y_y_xx = buffer_1100_sppd[402];

    auto g_z_y_0_0_0_y_y_xy = buffer_1100_sppd[403];

    auto g_z_y_0_0_0_y_y_xz = buffer_1100_sppd[404];

    auto g_z_y_0_0_0_y_y_yy = buffer_1100_sppd[405];

    auto g_z_y_0_0_0_y_y_yz = buffer_1100_sppd[406];

    auto g_z_y_0_0_0_y_y_zz = buffer_1100_sppd[407];

    auto g_z_y_0_0_0_y_z_xx = buffer_1100_sppd[408];

    auto g_z_y_0_0_0_y_z_xy = buffer_1100_sppd[409];

    auto g_z_y_0_0_0_y_z_xz = buffer_1100_sppd[410];

    auto g_z_y_0_0_0_y_z_yy = buffer_1100_sppd[411];

    auto g_z_y_0_0_0_y_z_yz = buffer_1100_sppd[412];

    auto g_z_y_0_0_0_y_z_zz = buffer_1100_sppd[413];

    auto g_z_y_0_0_0_z_x_xx = buffer_1100_sppd[414];

    auto g_z_y_0_0_0_z_x_xy = buffer_1100_sppd[415];

    auto g_z_y_0_0_0_z_x_xz = buffer_1100_sppd[416];

    auto g_z_y_0_0_0_z_x_yy = buffer_1100_sppd[417];

    auto g_z_y_0_0_0_z_x_yz = buffer_1100_sppd[418];

    auto g_z_y_0_0_0_z_x_zz = buffer_1100_sppd[419];

    auto g_z_y_0_0_0_z_y_xx = buffer_1100_sppd[420];

    auto g_z_y_0_0_0_z_y_xy = buffer_1100_sppd[421];

    auto g_z_y_0_0_0_z_y_xz = buffer_1100_sppd[422];

    auto g_z_y_0_0_0_z_y_yy = buffer_1100_sppd[423];

    auto g_z_y_0_0_0_z_y_yz = buffer_1100_sppd[424];

    auto g_z_y_0_0_0_z_y_zz = buffer_1100_sppd[425];

    auto g_z_y_0_0_0_z_z_xx = buffer_1100_sppd[426];

    auto g_z_y_0_0_0_z_z_xy = buffer_1100_sppd[427];

    auto g_z_y_0_0_0_z_z_xz = buffer_1100_sppd[428];

    auto g_z_y_0_0_0_z_z_yy = buffer_1100_sppd[429];

    auto g_z_y_0_0_0_z_z_yz = buffer_1100_sppd[430];

    auto g_z_y_0_0_0_z_z_zz = buffer_1100_sppd[431];

    auto g_z_z_0_0_0_x_x_xx = buffer_1100_sppd[432];

    auto g_z_z_0_0_0_x_x_xy = buffer_1100_sppd[433];

    auto g_z_z_0_0_0_x_x_xz = buffer_1100_sppd[434];

    auto g_z_z_0_0_0_x_x_yy = buffer_1100_sppd[435];

    auto g_z_z_0_0_0_x_x_yz = buffer_1100_sppd[436];

    auto g_z_z_0_0_0_x_x_zz = buffer_1100_sppd[437];

    auto g_z_z_0_0_0_x_y_xx = buffer_1100_sppd[438];

    auto g_z_z_0_0_0_x_y_xy = buffer_1100_sppd[439];

    auto g_z_z_0_0_0_x_y_xz = buffer_1100_sppd[440];

    auto g_z_z_0_0_0_x_y_yy = buffer_1100_sppd[441];

    auto g_z_z_0_0_0_x_y_yz = buffer_1100_sppd[442];

    auto g_z_z_0_0_0_x_y_zz = buffer_1100_sppd[443];

    auto g_z_z_0_0_0_x_z_xx = buffer_1100_sppd[444];

    auto g_z_z_0_0_0_x_z_xy = buffer_1100_sppd[445];

    auto g_z_z_0_0_0_x_z_xz = buffer_1100_sppd[446];

    auto g_z_z_0_0_0_x_z_yy = buffer_1100_sppd[447];

    auto g_z_z_0_0_0_x_z_yz = buffer_1100_sppd[448];

    auto g_z_z_0_0_0_x_z_zz = buffer_1100_sppd[449];

    auto g_z_z_0_0_0_y_x_xx = buffer_1100_sppd[450];

    auto g_z_z_0_0_0_y_x_xy = buffer_1100_sppd[451];

    auto g_z_z_0_0_0_y_x_xz = buffer_1100_sppd[452];

    auto g_z_z_0_0_0_y_x_yy = buffer_1100_sppd[453];

    auto g_z_z_0_0_0_y_x_yz = buffer_1100_sppd[454];

    auto g_z_z_0_0_0_y_x_zz = buffer_1100_sppd[455];

    auto g_z_z_0_0_0_y_y_xx = buffer_1100_sppd[456];

    auto g_z_z_0_0_0_y_y_xy = buffer_1100_sppd[457];

    auto g_z_z_0_0_0_y_y_xz = buffer_1100_sppd[458];

    auto g_z_z_0_0_0_y_y_yy = buffer_1100_sppd[459];

    auto g_z_z_0_0_0_y_y_yz = buffer_1100_sppd[460];

    auto g_z_z_0_0_0_y_y_zz = buffer_1100_sppd[461];

    auto g_z_z_0_0_0_y_z_xx = buffer_1100_sppd[462];

    auto g_z_z_0_0_0_y_z_xy = buffer_1100_sppd[463];

    auto g_z_z_0_0_0_y_z_xz = buffer_1100_sppd[464];

    auto g_z_z_0_0_0_y_z_yy = buffer_1100_sppd[465];

    auto g_z_z_0_0_0_y_z_yz = buffer_1100_sppd[466];

    auto g_z_z_0_0_0_y_z_zz = buffer_1100_sppd[467];

    auto g_z_z_0_0_0_z_x_xx = buffer_1100_sppd[468];

    auto g_z_z_0_0_0_z_x_xy = buffer_1100_sppd[469];

    auto g_z_z_0_0_0_z_x_xz = buffer_1100_sppd[470];

    auto g_z_z_0_0_0_z_x_yy = buffer_1100_sppd[471];

    auto g_z_z_0_0_0_z_x_yz = buffer_1100_sppd[472];

    auto g_z_z_0_0_0_z_x_zz = buffer_1100_sppd[473];

    auto g_z_z_0_0_0_z_y_xx = buffer_1100_sppd[474];

    auto g_z_z_0_0_0_z_y_xy = buffer_1100_sppd[475];

    auto g_z_z_0_0_0_z_y_xz = buffer_1100_sppd[476];

    auto g_z_z_0_0_0_z_y_yy = buffer_1100_sppd[477];

    auto g_z_z_0_0_0_z_y_yz = buffer_1100_sppd[478];

    auto g_z_z_0_0_0_z_y_zz = buffer_1100_sppd[479];

    auto g_z_z_0_0_0_z_z_xx = buffer_1100_sppd[480];

    auto g_z_z_0_0_0_z_z_xy = buffer_1100_sppd[481];

    auto g_z_z_0_0_0_z_z_xz = buffer_1100_sppd[482];

    auto g_z_z_0_0_0_z_z_yy = buffer_1100_sppd[483];

    auto g_z_z_0_0_0_z_z_yz = buffer_1100_sppd[484];

    auto g_z_z_0_0_0_z_z_zz = buffer_1100_sppd[485];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_x_xx, g_x_0_x_xy, g_x_0_x_xz, g_x_0_x_yy, g_x_0_x_yz, g_x_0_x_zz, g_x_x_0_0_0_x_x_xx, g_x_x_0_0_0_x_x_xy, g_x_x_0_0_0_x_x_xz, g_x_x_0_0_0_x_x_yy, g_x_x_0_0_0_x_x_yz, g_x_x_0_0_0_x_x_zz, g_x_xx_x_xx, g_x_xx_x_xy, g_x_xx_x_xz, g_x_xx_x_yy, g_x_xx_x_yz, g_x_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_x_x_xx[i] = -2.0 * g_x_0_x_xx[i] * a_exp + 4.0 * g_x_xx_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_x_xy[i] = -2.0 * g_x_0_x_xy[i] * a_exp + 4.0 * g_x_xx_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_x_xz[i] = -2.0 * g_x_0_x_xz[i] * a_exp + 4.0 * g_x_xx_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_x_yy[i] = -2.0 * g_x_0_x_yy[i] * a_exp + 4.0 * g_x_xx_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_x_yz[i] = -2.0 * g_x_0_x_yz[i] * a_exp + 4.0 * g_x_xx_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_x_zz[i] = -2.0 * g_x_0_x_zz[i] * a_exp + 4.0 * g_x_xx_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_y_xx, g_x_0_y_xy, g_x_0_y_xz, g_x_0_y_yy, g_x_0_y_yz, g_x_0_y_zz, g_x_x_0_0_0_x_y_xx, g_x_x_0_0_0_x_y_xy, g_x_x_0_0_0_x_y_xz, g_x_x_0_0_0_x_y_yy, g_x_x_0_0_0_x_y_yz, g_x_x_0_0_0_x_y_zz, g_x_xx_y_xx, g_x_xx_y_xy, g_x_xx_y_xz, g_x_xx_y_yy, g_x_xx_y_yz, g_x_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_x_y_xx[i] = -2.0 * g_x_0_y_xx[i] * a_exp + 4.0 * g_x_xx_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_y_xy[i] = -2.0 * g_x_0_y_xy[i] * a_exp + 4.0 * g_x_xx_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_y_xz[i] = -2.0 * g_x_0_y_xz[i] * a_exp + 4.0 * g_x_xx_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_y_yy[i] = -2.0 * g_x_0_y_yy[i] * a_exp + 4.0 * g_x_xx_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_y_yz[i] = -2.0 * g_x_0_y_yz[i] * a_exp + 4.0 * g_x_xx_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_y_zz[i] = -2.0 * g_x_0_y_zz[i] * a_exp + 4.0 * g_x_xx_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_z_xx, g_x_0_z_xy, g_x_0_z_xz, g_x_0_z_yy, g_x_0_z_yz, g_x_0_z_zz, g_x_x_0_0_0_x_z_xx, g_x_x_0_0_0_x_z_xy, g_x_x_0_0_0_x_z_xz, g_x_x_0_0_0_x_z_yy, g_x_x_0_0_0_x_z_yz, g_x_x_0_0_0_x_z_zz, g_x_xx_z_xx, g_x_xx_z_xy, g_x_xx_z_xz, g_x_xx_z_yy, g_x_xx_z_yz, g_x_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_x_z_xx[i] = -2.0 * g_x_0_z_xx[i] * a_exp + 4.0 * g_x_xx_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_z_xy[i] = -2.0 * g_x_0_z_xy[i] * a_exp + 4.0 * g_x_xx_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_z_xz[i] = -2.0 * g_x_0_z_xz[i] * a_exp + 4.0 * g_x_xx_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_z_yy[i] = -2.0 * g_x_0_z_yy[i] * a_exp + 4.0 * g_x_xx_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_z_yz[i] = -2.0 * g_x_0_z_yz[i] * a_exp + 4.0 * g_x_xx_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_z_zz[i] = -2.0 * g_x_0_z_zz[i] * a_exp + 4.0 * g_x_xx_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_x_0_0_0_y_x_xx, g_x_x_0_0_0_y_x_xy, g_x_x_0_0_0_y_x_xz, g_x_x_0_0_0_y_x_yy, g_x_x_0_0_0_y_x_yz, g_x_x_0_0_0_y_x_zz, g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_y_x_xx[i] = 4.0 * g_x_xy_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_x_xy[i] = 4.0 * g_x_xy_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_x_xz[i] = 4.0 * g_x_xy_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_x_yy[i] = 4.0 * g_x_xy_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_x_yz[i] = 4.0 * g_x_xy_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_x_zz[i] = 4.0 * g_x_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_x_0_0_0_y_y_xx, g_x_x_0_0_0_y_y_xy, g_x_x_0_0_0_y_y_xz, g_x_x_0_0_0_y_y_yy, g_x_x_0_0_0_y_y_yz, g_x_x_0_0_0_y_y_zz, g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_y_y_xx[i] = 4.0 * g_x_xy_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_y_xy[i] = 4.0 * g_x_xy_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_y_xz[i] = 4.0 * g_x_xy_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_y_yy[i] = 4.0 * g_x_xy_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_y_yz[i] = 4.0 * g_x_xy_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_y_zz[i] = 4.0 * g_x_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_x_0_0_0_y_z_xx, g_x_x_0_0_0_y_z_xy, g_x_x_0_0_0_y_z_xz, g_x_x_0_0_0_y_z_yy, g_x_x_0_0_0_y_z_yz, g_x_x_0_0_0_y_z_zz, g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_y_z_xx[i] = 4.0 * g_x_xy_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_z_xy[i] = 4.0 * g_x_xy_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_z_xz[i] = 4.0 * g_x_xy_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_z_yy[i] = 4.0 * g_x_xy_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_z_yz[i] = 4.0 * g_x_xy_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_z_zz[i] = 4.0 * g_x_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_x_0_0_0_z_x_xx, g_x_x_0_0_0_z_x_xy, g_x_x_0_0_0_z_x_xz, g_x_x_0_0_0_z_x_yy, g_x_x_0_0_0_z_x_yz, g_x_x_0_0_0_z_x_zz, g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_z_x_xx[i] = 4.0 * g_x_xz_x_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_x_xy[i] = 4.0 * g_x_xz_x_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_x_xz[i] = 4.0 * g_x_xz_x_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_x_yy[i] = 4.0 * g_x_xz_x_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_x_yz[i] = 4.0 * g_x_xz_x_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_x_zz[i] = 4.0 * g_x_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_x_0_0_0_z_y_xx, g_x_x_0_0_0_z_y_xy, g_x_x_0_0_0_z_y_xz, g_x_x_0_0_0_z_y_yy, g_x_x_0_0_0_z_y_yz, g_x_x_0_0_0_z_y_zz, g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_z_y_xx[i] = 4.0 * g_x_xz_y_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_y_xy[i] = 4.0 * g_x_xz_y_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_y_xz[i] = 4.0 * g_x_xz_y_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_y_yy[i] = 4.0 * g_x_xz_y_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_y_yz[i] = 4.0 * g_x_xz_y_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_y_zz[i] = 4.0 * g_x_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_x_0_0_0_z_z_xx, g_x_x_0_0_0_z_z_xy, g_x_x_0_0_0_z_z_xz, g_x_x_0_0_0_z_z_yy, g_x_x_0_0_0_z_z_yz, g_x_x_0_0_0_z_z_zz, g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_z_z_xx[i] = 4.0 * g_x_xz_z_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_z_xy[i] = 4.0 * g_x_xz_z_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_z_xz[i] = 4.0 * g_x_xz_z_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_z_yy[i] = 4.0 * g_x_xz_z_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_z_yz[i] = 4.0 * g_x_xz_z_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_z_zz[i] = 4.0 * g_x_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_xy_x_xx, g_x_xy_x_xy, g_x_xy_x_xz, g_x_xy_x_yy, g_x_xy_x_yz, g_x_xy_x_zz, g_x_y_0_0_0_x_x_xx, g_x_y_0_0_0_x_x_xy, g_x_y_0_0_0_x_x_xz, g_x_y_0_0_0_x_x_yy, g_x_y_0_0_0_x_x_yz, g_x_y_0_0_0_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_x_x_xx[i] = 4.0 * g_x_xy_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_x_xy[i] = 4.0 * g_x_xy_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_x_xz[i] = 4.0 * g_x_xy_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_x_yy[i] = 4.0 * g_x_xy_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_x_yz[i] = 4.0 * g_x_xy_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_x_zz[i] = 4.0 * g_x_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_xy_y_xx, g_x_xy_y_xy, g_x_xy_y_xz, g_x_xy_y_yy, g_x_xy_y_yz, g_x_xy_y_zz, g_x_y_0_0_0_x_y_xx, g_x_y_0_0_0_x_y_xy, g_x_y_0_0_0_x_y_xz, g_x_y_0_0_0_x_y_yy, g_x_y_0_0_0_x_y_yz, g_x_y_0_0_0_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_x_y_xx[i] = 4.0 * g_x_xy_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_y_xy[i] = 4.0 * g_x_xy_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_y_xz[i] = 4.0 * g_x_xy_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_y_yy[i] = 4.0 * g_x_xy_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_y_yz[i] = 4.0 * g_x_xy_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_y_zz[i] = 4.0 * g_x_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_xy_z_xx, g_x_xy_z_xy, g_x_xy_z_xz, g_x_xy_z_yy, g_x_xy_z_yz, g_x_xy_z_zz, g_x_y_0_0_0_x_z_xx, g_x_y_0_0_0_x_z_xy, g_x_y_0_0_0_x_z_xz, g_x_y_0_0_0_x_z_yy, g_x_y_0_0_0_x_z_yz, g_x_y_0_0_0_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_x_z_xx[i] = 4.0 * g_x_xy_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_z_xy[i] = 4.0 * g_x_xy_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_z_xz[i] = 4.0 * g_x_xy_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_z_yy[i] = 4.0 * g_x_xy_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_z_yz[i] = 4.0 * g_x_xy_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_z_zz[i] = 4.0 * g_x_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_x_xx, g_x_0_x_xy, g_x_0_x_xz, g_x_0_x_yy, g_x_0_x_yz, g_x_0_x_zz, g_x_y_0_0_0_y_x_xx, g_x_y_0_0_0_y_x_xy, g_x_y_0_0_0_y_x_xz, g_x_y_0_0_0_y_x_yy, g_x_y_0_0_0_y_x_yz, g_x_y_0_0_0_y_x_zz, g_x_yy_x_xx, g_x_yy_x_xy, g_x_yy_x_xz, g_x_yy_x_yy, g_x_yy_x_yz, g_x_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_y_x_xx[i] = -2.0 * g_x_0_x_xx[i] * a_exp + 4.0 * g_x_yy_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_x_xy[i] = -2.0 * g_x_0_x_xy[i] * a_exp + 4.0 * g_x_yy_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_x_xz[i] = -2.0 * g_x_0_x_xz[i] * a_exp + 4.0 * g_x_yy_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_x_yy[i] = -2.0 * g_x_0_x_yy[i] * a_exp + 4.0 * g_x_yy_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_x_yz[i] = -2.0 * g_x_0_x_yz[i] * a_exp + 4.0 * g_x_yy_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_x_zz[i] = -2.0 * g_x_0_x_zz[i] * a_exp + 4.0 * g_x_yy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_y_xx, g_x_0_y_xy, g_x_0_y_xz, g_x_0_y_yy, g_x_0_y_yz, g_x_0_y_zz, g_x_y_0_0_0_y_y_xx, g_x_y_0_0_0_y_y_xy, g_x_y_0_0_0_y_y_xz, g_x_y_0_0_0_y_y_yy, g_x_y_0_0_0_y_y_yz, g_x_y_0_0_0_y_y_zz, g_x_yy_y_xx, g_x_yy_y_xy, g_x_yy_y_xz, g_x_yy_y_yy, g_x_yy_y_yz, g_x_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_y_y_xx[i] = -2.0 * g_x_0_y_xx[i] * a_exp + 4.0 * g_x_yy_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_y_xy[i] = -2.0 * g_x_0_y_xy[i] * a_exp + 4.0 * g_x_yy_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_y_xz[i] = -2.0 * g_x_0_y_xz[i] * a_exp + 4.0 * g_x_yy_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_y_yy[i] = -2.0 * g_x_0_y_yy[i] * a_exp + 4.0 * g_x_yy_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_y_yz[i] = -2.0 * g_x_0_y_yz[i] * a_exp + 4.0 * g_x_yy_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_y_zz[i] = -2.0 * g_x_0_y_zz[i] * a_exp + 4.0 * g_x_yy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_z_xx, g_x_0_z_xy, g_x_0_z_xz, g_x_0_z_yy, g_x_0_z_yz, g_x_0_z_zz, g_x_y_0_0_0_y_z_xx, g_x_y_0_0_0_y_z_xy, g_x_y_0_0_0_y_z_xz, g_x_y_0_0_0_y_z_yy, g_x_y_0_0_0_y_z_yz, g_x_y_0_0_0_y_z_zz, g_x_yy_z_xx, g_x_yy_z_xy, g_x_yy_z_xz, g_x_yy_z_yy, g_x_yy_z_yz, g_x_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_y_z_xx[i] = -2.0 * g_x_0_z_xx[i] * a_exp + 4.0 * g_x_yy_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_z_xy[i] = -2.0 * g_x_0_z_xy[i] * a_exp + 4.0 * g_x_yy_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_z_xz[i] = -2.0 * g_x_0_z_xz[i] * a_exp + 4.0 * g_x_yy_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_z_yy[i] = -2.0 * g_x_0_z_yy[i] * a_exp + 4.0 * g_x_yy_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_z_yz[i] = -2.0 * g_x_0_z_yz[i] * a_exp + 4.0 * g_x_yy_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_z_zz[i] = -2.0 * g_x_0_z_zz[i] * a_exp + 4.0 * g_x_yy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_y_0_0_0_z_x_xx, g_x_y_0_0_0_z_x_xy, g_x_y_0_0_0_z_x_xz, g_x_y_0_0_0_z_x_yy, g_x_y_0_0_0_z_x_yz, g_x_y_0_0_0_z_x_zz, g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_z_x_xx[i] = 4.0 * g_x_yz_x_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_x_xy[i] = 4.0 * g_x_yz_x_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_x_xz[i] = 4.0 * g_x_yz_x_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_x_yy[i] = 4.0 * g_x_yz_x_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_x_yz[i] = 4.0 * g_x_yz_x_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_x_zz[i] = 4.0 * g_x_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_y_0_0_0_z_y_xx, g_x_y_0_0_0_z_y_xy, g_x_y_0_0_0_z_y_xz, g_x_y_0_0_0_z_y_yy, g_x_y_0_0_0_z_y_yz, g_x_y_0_0_0_z_y_zz, g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_z_y_xx[i] = 4.0 * g_x_yz_y_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_y_xy[i] = 4.0 * g_x_yz_y_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_y_xz[i] = 4.0 * g_x_yz_y_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_y_yy[i] = 4.0 * g_x_yz_y_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_y_yz[i] = 4.0 * g_x_yz_y_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_y_zz[i] = 4.0 * g_x_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_y_0_0_0_z_z_xx, g_x_y_0_0_0_z_z_xy, g_x_y_0_0_0_z_z_xz, g_x_y_0_0_0_z_z_yy, g_x_y_0_0_0_z_z_yz, g_x_y_0_0_0_z_z_zz, g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_z_z_xx[i] = 4.0 * g_x_yz_z_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_z_xy[i] = 4.0 * g_x_yz_z_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_z_xz[i] = 4.0 * g_x_yz_z_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_z_yy[i] = 4.0 * g_x_yz_z_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_z_yz[i] = 4.0 * g_x_yz_z_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_z_zz[i] = 4.0 * g_x_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_xz_x_xx, g_x_xz_x_xy, g_x_xz_x_xz, g_x_xz_x_yy, g_x_xz_x_yz, g_x_xz_x_zz, g_x_z_0_0_0_x_x_xx, g_x_z_0_0_0_x_x_xy, g_x_z_0_0_0_x_x_xz, g_x_z_0_0_0_x_x_yy, g_x_z_0_0_0_x_x_yz, g_x_z_0_0_0_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_x_x_xx[i] = 4.0 * g_x_xz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_x_xy[i] = 4.0 * g_x_xz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_x_xz[i] = 4.0 * g_x_xz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_x_yy[i] = 4.0 * g_x_xz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_x_yz[i] = 4.0 * g_x_xz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_x_zz[i] = 4.0 * g_x_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_xz_y_xx, g_x_xz_y_xy, g_x_xz_y_xz, g_x_xz_y_yy, g_x_xz_y_yz, g_x_xz_y_zz, g_x_z_0_0_0_x_y_xx, g_x_z_0_0_0_x_y_xy, g_x_z_0_0_0_x_y_xz, g_x_z_0_0_0_x_y_yy, g_x_z_0_0_0_x_y_yz, g_x_z_0_0_0_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_x_y_xx[i] = 4.0 * g_x_xz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_y_xy[i] = 4.0 * g_x_xz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_y_xz[i] = 4.0 * g_x_xz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_y_yy[i] = 4.0 * g_x_xz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_y_yz[i] = 4.0 * g_x_xz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_y_zz[i] = 4.0 * g_x_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_xz_z_xx, g_x_xz_z_xy, g_x_xz_z_xz, g_x_xz_z_yy, g_x_xz_z_yz, g_x_xz_z_zz, g_x_z_0_0_0_x_z_xx, g_x_z_0_0_0_x_z_xy, g_x_z_0_0_0_x_z_xz, g_x_z_0_0_0_x_z_yy, g_x_z_0_0_0_x_z_yz, g_x_z_0_0_0_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_x_z_xx[i] = 4.0 * g_x_xz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_z_xy[i] = 4.0 * g_x_xz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_z_xz[i] = 4.0 * g_x_xz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_z_yy[i] = 4.0 * g_x_xz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_z_yz[i] = 4.0 * g_x_xz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_z_zz[i] = 4.0 * g_x_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_yz_x_xx, g_x_yz_x_xy, g_x_yz_x_xz, g_x_yz_x_yy, g_x_yz_x_yz, g_x_yz_x_zz, g_x_z_0_0_0_y_x_xx, g_x_z_0_0_0_y_x_xy, g_x_z_0_0_0_y_x_xz, g_x_z_0_0_0_y_x_yy, g_x_z_0_0_0_y_x_yz, g_x_z_0_0_0_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_y_x_xx[i] = 4.0 * g_x_yz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_x_xy[i] = 4.0 * g_x_yz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_x_xz[i] = 4.0 * g_x_yz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_x_yy[i] = 4.0 * g_x_yz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_x_yz[i] = 4.0 * g_x_yz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_x_zz[i] = 4.0 * g_x_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_yz_y_xx, g_x_yz_y_xy, g_x_yz_y_xz, g_x_yz_y_yy, g_x_yz_y_yz, g_x_yz_y_zz, g_x_z_0_0_0_y_y_xx, g_x_z_0_0_0_y_y_xy, g_x_z_0_0_0_y_y_xz, g_x_z_0_0_0_y_y_yy, g_x_z_0_0_0_y_y_yz, g_x_z_0_0_0_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_y_y_xx[i] = 4.0 * g_x_yz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_y_xy[i] = 4.0 * g_x_yz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_y_xz[i] = 4.0 * g_x_yz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_y_yy[i] = 4.0 * g_x_yz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_y_yz[i] = 4.0 * g_x_yz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_y_zz[i] = 4.0 * g_x_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_yz_z_xx, g_x_yz_z_xy, g_x_yz_z_xz, g_x_yz_z_yy, g_x_yz_z_yz, g_x_yz_z_zz, g_x_z_0_0_0_y_z_xx, g_x_z_0_0_0_y_z_xy, g_x_z_0_0_0_y_z_xz, g_x_z_0_0_0_y_z_yy, g_x_z_0_0_0_y_z_yz, g_x_z_0_0_0_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_y_z_xx[i] = 4.0 * g_x_yz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_z_xy[i] = 4.0 * g_x_yz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_z_xz[i] = 4.0 * g_x_yz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_z_yy[i] = 4.0 * g_x_yz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_z_yz[i] = 4.0 * g_x_yz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_z_zz[i] = 4.0 * g_x_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_x_xx, g_x_0_x_xy, g_x_0_x_xz, g_x_0_x_yy, g_x_0_x_yz, g_x_0_x_zz, g_x_z_0_0_0_z_x_xx, g_x_z_0_0_0_z_x_xy, g_x_z_0_0_0_z_x_xz, g_x_z_0_0_0_z_x_yy, g_x_z_0_0_0_z_x_yz, g_x_z_0_0_0_z_x_zz, g_x_zz_x_xx, g_x_zz_x_xy, g_x_zz_x_xz, g_x_zz_x_yy, g_x_zz_x_yz, g_x_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_z_x_xx[i] = -2.0 * g_x_0_x_xx[i] * a_exp + 4.0 * g_x_zz_x_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_x_xy[i] = -2.0 * g_x_0_x_xy[i] * a_exp + 4.0 * g_x_zz_x_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_x_xz[i] = -2.0 * g_x_0_x_xz[i] * a_exp + 4.0 * g_x_zz_x_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_x_yy[i] = -2.0 * g_x_0_x_yy[i] * a_exp + 4.0 * g_x_zz_x_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_x_yz[i] = -2.0 * g_x_0_x_yz[i] * a_exp + 4.0 * g_x_zz_x_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_x_zz[i] = -2.0 * g_x_0_x_zz[i] * a_exp + 4.0 * g_x_zz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_y_xx, g_x_0_y_xy, g_x_0_y_xz, g_x_0_y_yy, g_x_0_y_yz, g_x_0_y_zz, g_x_z_0_0_0_z_y_xx, g_x_z_0_0_0_z_y_xy, g_x_z_0_0_0_z_y_xz, g_x_z_0_0_0_z_y_yy, g_x_z_0_0_0_z_y_yz, g_x_z_0_0_0_z_y_zz, g_x_zz_y_xx, g_x_zz_y_xy, g_x_zz_y_xz, g_x_zz_y_yy, g_x_zz_y_yz, g_x_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_z_y_xx[i] = -2.0 * g_x_0_y_xx[i] * a_exp + 4.0 * g_x_zz_y_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_y_xy[i] = -2.0 * g_x_0_y_xy[i] * a_exp + 4.0 * g_x_zz_y_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_y_xz[i] = -2.0 * g_x_0_y_xz[i] * a_exp + 4.0 * g_x_zz_y_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_y_yy[i] = -2.0 * g_x_0_y_yy[i] * a_exp + 4.0 * g_x_zz_y_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_y_yz[i] = -2.0 * g_x_0_y_yz[i] * a_exp + 4.0 * g_x_zz_y_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_y_zz[i] = -2.0 * g_x_0_y_zz[i] * a_exp + 4.0 * g_x_zz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_z_xx, g_x_0_z_xy, g_x_0_z_xz, g_x_0_z_yy, g_x_0_z_yz, g_x_0_z_zz, g_x_z_0_0_0_z_z_xx, g_x_z_0_0_0_z_z_xy, g_x_z_0_0_0_z_z_xz, g_x_z_0_0_0_z_z_yy, g_x_z_0_0_0_z_z_yz, g_x_z_0_0_0_z_z_zz, g_x_zz_z_xx, g_x_zz_z_xy, g_x_zz_z_xz, g_x_zz_z_yy, g_x_zz_z_yz, g_x_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_z_z_xx[i] = -2.0 * g_x_0_z_xx[i] * a_exp + 4.0 * g_x_zz_z_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_z_xy[i] = -2.0 * g_x_0_z_xy[i] * a_exp + 4.0 * g_x_zz_z_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_z_xz[i] = -2.0 * g_x_0_z_xz[i] * a_exp + 4.0 * g_x_zz_z_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_z_yy[i] = -2.0 * g_x_0_z_yy[i] * a_exp + 4.0 * g_x_zz_z_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_z_yz[i] = -2.0 * g_x_0_z_yz[i] * a_exp + 4.0 * g_x_zz_z_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_z_zz[i] = -2.0 * g_x_0_z_zz[i] * a_exp + 4.0 * g_x_zz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_y_0_x_xx, g_y_0_x_xy, g_y_0_x_xz, g_y_0_x_yy, g_y_0_x_yz, g_y_0_x_zz, g_y_x_0_0_0_x_x_xx, g_y_x_0_0_0_x_x_xy, g_y_x_0_0_0_x_x_xz, g_y_x_0_0_0_x_x_yy, g_y_x_0_0_0_x_x_yz, g_y_x_0_0_0_x_x_zz, g_y_xx_x_xx, g_y_xx_x_xy, g_y_xx_x_xz, g_y_xx_x_yy, g_y_xx_x_yz, g_y_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_x_x_xx[i] = -2.0 * g_y_0_x_xx[i] * a_exp + 4.0 * g_y_xx_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_x_xy[i] = -2.0 * g_y_0_x_xy[i] * a_exp + 4.0 * g_y_xx_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_x_xz[i] = -2.0 * g_y_0_x_xz[i] * a_exp + 4.0 * g_y_xx_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_x_yy[i] = -2.0 * g_y_0_x_yy[i] * a_exp + 4.0 * g_y_xx_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_x_yz[i] = -2.0 * g_y_0_x_yz[i] * a_exp + 4.0 * g_y_xx_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_x_zz[i] = -2.0 * g_y_0_x_zz[i] * a_exp + 4.0 * g_y_xx_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_y_0_y_xx, g_y_0_y_xy, g_y_0_y_xz, g_y_0_y_yy, g_y_0_y_yz, g_y_0_y_zz, g_y_x_0_0_0_x_y_xx, g_y_x_0_0_0_x_y_xy, g_y_x_0_0_0_x_y_xz, g_y_x_0_0_0_x_y_yy, g_y_x_0_0_0_x_y_yz, g_y_x_0_0_0_x_y_zz, g_y_xx_y_xx, g_y_xx_y_xy, g_y_xx_y_xz, g_y_xx_y_yy, g_y_xx_y_yz, g_y_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_x_y_xx[i] = -2.0 * g_y_0_y_xx[i] * a_exp + 4.0 * g_y_xx_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_y_xy[i] = -2.0 * g_y_0_y_xy[i] * a_exp + 4.0 * g_y_xx_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_y_xz[i] = -2.0 * g_y_0_y_xz[i] * a_exp + 4.0 * g_y_xx_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_y_yy[i] = -2.0 * g_y_0_y_yy[i] * a_exp + 4.0 * g_y_xx_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_y_yz[i] = -2.0 * g_y_0_y_yz[i] * a_exp + 4.0 * g_y_xx_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_y_zz[i] = -2.0 * g_y_0_y_zz[i] * a_exp + 4.0 * g_y_xx_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_y_0_z_xx, g_y_0_z_xy, g_y_0_z_xz, g_y_0_z_yy, g_y_0_z_yz, g_y_0_z_zz, g_y_x_0_0_0_x_z_xx, g_y_x_0_0_0_x_z_xy, g_y_x_0_0_0_x_z_xz, g_y_x_0_0_0_x_z_yy, g_y_x_0_0_0_x_z_yz, g_y_x_0_0_0_x_z_zz, g_y_xx_z_xx, g_y_xx_z_xy, g_y_xx_z_xz, g_y_xx_z_yy, g_y_xx_z_yz, g_y_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_x_z_xx[i] = -2.0 * g_y_0_z_xx[i] * a_exp + 4.0 * g_y_xx_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_z_xy[i] = -2.0 * g_y_0_z_xy[i] * a_exp + 4.0 * g_y_xx_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_z_xz[i] = -2.0 * g_y_0_z_xz[i] * a_exp + 4.0 * g_y_xx_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_z_yy[i] = -2.0 * g_y_0_z_yy[i] * a_exp + 4.0 * g_y_xx_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_z_yz[i] = -2.0 * g_y_0_z_yz[i] * a_exp + 4.0 * g_y_xx_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_z_zz[i] = -2.0 * g_y_0_z_zz[i] * a_exp + 4.0 * g_y_xx_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_y_x_0_0_0_y_x_xx, g_y_x_0_0_0_y_x_xy, g_y_x_0_0_0_y_x_xz, g_y_x_0_0_0_y_x_yy, g_y_x_0_0_0_y_x_yz, g_y_x_0_0_0_y_x_zz, g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_y_x_xx[i] = 4.0 * g_y_xy_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_x_xy[i] = 4.0 * g_y_xy_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_x_xz[i] = 4.0 * g_y_xy_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_x_yy[i] = 4.0 * g_y_xy_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_x_yz[i] = 4.0 * g_y_xy_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_x_zz[i] = 4.0 * g_y_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_y_x_0_0_0_y_y_xx, g_y_x_0_0_0_y_y_xy, g_y_x_0_0_0_y_y_xz, g_y_x_0_0_0_y_y_yy, g_y_x_0_0_0_y_y_yz, g_y_x_0_0_0_y_y_zz, g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_y_y_xx[i] = 4.0 * g_y_xy_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_y_xy[i] = 4.0 * g_y_xy_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_y_xz[i] = 4.0 * g_y_xy_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_y_yy[i] = 4.0 * g_y_xy_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_y_yz[i] = 4.0 * g_y_xy_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_y_zz[i] = 4.0 * g_y_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_y_x_0_0_0_y_z_xx, g_y_x_0_0_0_y_z_xy, g_y_x_0_0_0_y_z_xz, g_y_x_0_0_0_y_z_yy, g_y_x_0_0_0_y_z_yz, g_y_x_0_0_0_y_z_zz, g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_y_z_xx[i] = 4.0 * g_y_xy_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_z_xy[i] = 4.0 * g_y_xy_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_z_xz[i] = 4.0 * g_y_xy_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_z_yy[i] = 4.0 * g_y_xy_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_z_yz[i] = 4.0 * g_y_xy_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_z_zz[i] = 4.0 * g_y_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_y_x_0_0_0_z_x_xx, g_y_x_0_0_0_z_x_xy, g_y_x_0_0_0_z_x_xz, g_y_x_0_0_0_z_x_yy, g_y_x_0_0_0_z_x_yz, g_y_x_0_0_0_z_x_zz, g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_z_x_xx[i] = 4.0 * g_y_xz_x_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_x_xy[i] = 4.0 * g_y_xz_x_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_x_xz[i] = 4.0 * g_y_xz_x_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_x_yy[i] = 4.0 * g_y_xz_x_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_x_yz[i] = 4.0 * g_y_xz_x_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_x_zz[i] = 4.0 * g_y_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_y_x_0_0_0_z_y_xx, g_y_x_0_0_0_z_y_xy, g_y_x_0_0_0_z_y_xz, g_y_x_0_0_0_z_y_yy, g_y_x_0_0_0_z_y_yz, g_y_x_0_0_0_z_y_zz, g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_z_y_xx[i] = 4.0 * g_y_xz_y_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_y_xy[i] = 4.0 * g_y_xz_y_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_y_xz[i] = 4.0 * g_y_xz_y_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_y_yy[i] = 4.0 * g_y_xz_y_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_y_yz[i] = 4.0 * g_y_xz_y_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_y_zz[i] = 4.0 * g_y_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_y_x_0_0_0_z_z_xx, g_y_x_0_0_0_z_z_xy, g_y_x_0_0_0_z_z_xz, g_y_x_0_0_0_z_z_yy, g_y_x_0_0_0_z_z_yz, g_y_x_0_0_0_z_z_zz, g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_z_z_xx[i] = 4.0 * g_y_xz_z_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_z_xy[i] = 4.0 * g_y_xz_z_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_z_xz[i] = 4.0 * g_y_xz_z_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_z_yy[i] = 4.0 * g_y_xz_z_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_z_yz[i] = 4.0 * g_y_xz_z_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_z_zz[i] = 4.0 * g_y_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_y_xy_x_xx, g_y_xy_x_xy, g_y_xy_x_xz, g_y_xy_x_yy, g_y_xy_x_yz, g_y_xy_x_zz, g_y_y_0_0_0_x_x_xx, g_y_y_0_0_0_x_x_xy, g_y_y_0_0_0_x_x_xz, g_y_y_0_0_0_x_x_yy, g_y_y_0_0_0_x_x_yz, g_y_y_0_0_0_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_x_x_xx[i] = 4.0 * g_y_xy_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_x_xy[i] = 4.0 * g_y_xy_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_x_xz[i] = 4.0 * g_y_xy_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_x_yy[i] = 4.0 * g_y_xy_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_x_yz[i] = 4.0 * g_y_xy_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_x_zz[i] = 4.0 * g_y_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_y_xy_y_xx, g_y_xy_y_xy, g_y_xy_y_xz, g_y_xy_y_yy, g_y_xy_y_yz, g_y_xy_y_zz, g_y_y_0_0_0_x_y_xx, g_y_y_0_0_0_x_y_xy, g_y_y_0_0_0_x_y_xz, g_y_y_0_0_0_x_y_yy, g_y_y_0_0_0_x_y_yz, g_y_y_0_0_0_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_x_y_xx[i] = 4.0 * g_y_xy_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_y_xy[i] = 4.0 * g_y_xy_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_y_xz[i] = 4.0 * g_y_xy_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_y_yy[i] = 4.0 * g_y_xy_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_y_yz[i] = 4.0 * g_y_xy_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_y_zz[i] = 4.0 * g_y_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_y_xy_z_xx, g_y_xy_z_xy, g_y_xy_z_xz, g_y_xy_z_yy, g_y_xy_z_yz, g_y_xy_z_zz, g_y_y_0_0_0_x_z_xx, g_y_y_0_0_0_x_z_xy, g_y_y_0_0_0_x_z_xz, g_y_y_0_0_0_x_z_yy, g_y_y_0_0_0_x_z_yz, g_y_y_0_0_0_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_x_z_xx[i] = 4.0 * g_y_xy_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_z_xy[i] = 4.0 * g_y_xy_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_z_xz[i] = 4.0 * g_y_xy_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_z_yy[i] = 4.0 * g_y_xy_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_z_yz[i] = 4.0 * g_y_xy_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_z_zz[i] = 4.0 * g_y_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_y_0_x_xx, g_y_0_x_xy, g_y_0_x_xz, g_y_0_x_yy, g_y_0_x_yz, g_y_0_x_zz, g_y_y_0_0_0_y_x_xx, g_y_y_0_0_0_y_x_xy, g_y_y_0_0_0_y_x_xz, g_y_y_0_0_0_y_x_yy, g_y_y_0_0_0_y_x_yz, g_y_y_0_0_0_y_x_zz, g_y_yy_x_xx, g_y_yy_x_xy, g_y_yy_x_xz, g_y_yy_x_yy, g_y_yy_x_yz, g_y_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_y_x_xx[i] = -2.0 * g_y_0_x_xx[i] * a_exp + 4.0 * g_y_yy_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_x_xy[i] = -2.0 * g_y_0_x_xy[i] * a_exp + 4.0 * g_y_yy_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_x_xz[i] = -2.0 * g_y_0_x_xz[i] * a_exp + 4.0 * g_y_yy_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_x_yy[i] = -2.0 * g_y_0_x_yy[i] * a_exp + 4.0 * g_y_yy_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_x_yz[i] = -2.0 * g_y_0_x_yz[i] * a_exp + 4.0 * g_y_yy_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_x_zz[i] = -2.0 * g_y_0_x_zz[i] * a_exp + 4.0 * g_y_yy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_y_0_y_xx, g_y_0_y_xy, g_y_0_y_xz, g_y_0_y_yy, g_y_0_y_yz, g_y_0_y_zz, g_y_y_0_0_0_y_y_xx, g_y_y_0_0_0_y_y_xy, g_y_y_0_0_0_y_y_xz, g_y_y_0_0_0_y_y_yy, g_y_y_0_0_0_y_y_yz, g_y_y_0_0_0_y_y_zz, g_y_yy_y_xx, g_y_yy_y_xy, g_y_yy_y_xz, g_y_yy_y_yy, g_y_yy_y_yz, g_y_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_y_y_xx[i] = -2.0 * g_y_0_y_xx[i] * a_exp + 4.0 * g_y_yy_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_y_xy[i] = -2.0 * g_y_0_y_xy[i] * a_exp + 4.0 * g_y_yy_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_y_xz[i] = -2.0 * g_y_0_y_xz[i] * a_exp + 4.0 * g_y_yy_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_y_yy[i] = -2.0 * g_y_0_y_yy[i] * a_exp + 4.0 * g_y_yy_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_y_yz[i] = -2.0 * g_y_0_y_yz[i] * a_exp + 4.0 * g_y_yy_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_y_zz[i] = -2.0 * g_y_0_y_zz[i] * a_exp + 4.0 * g_y_yy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_y_0_z_xx, g_y_0_z_xy, g_y_0_z_xz, g_y_0_z_yy, g_y_0_z_yz, g_y_0_z_zz, g_y_y_0_0_0_y_z_xx, g_y_y_0_0_0_y_z_xy, g_y_y_0_0_0_y_z_xz, g_y_y_0_0_0_y_z_yy, g_y_y_0_0_0_y_z_yz, g_y_y_0_0_0_y_z_zz, g_y_yy_z_xx, g_y_yy_z_xy, g_y_yy_z_xz, g_y_yy_z_yy, g_y_yy_z_yz, g_y_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_y_z_xx[i] = -2.0 * g_y_0_z_xx[i] * a_exp + 4.0 * g_y_yy_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_z_xy[i] = -2.0 * g_y_0_z_xy[i] * a_exp + 4.0 * g_y_yy_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_z_xz[i] = -2.0 * g_y_0_z_xz[i] * a_exp + 4.0 * g_y_yy_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_z_yy[i] = -2.0 * g_y_0_z_yy[i] * a_exp + 4.0 * g_y_yy_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_z_yz[i] = -2.0 * g_y_0_z_yz[i] * a_exp + 4.0 * g_y_yy_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_z_zz[i] = -2.0 * g_y_0_z_zz[i] * a_exp + 4.0 * g_y_yy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_y_y_0_0_0_z_x_xx, g_y_y_0_0_0_z_x_xy, g_y_y_0_0_0_z_x_xz, g_y_y_0_0_0_z_x_yy, g_y_y_0_0_0_z_x_yz, g_y_y_0_0_0_z_x_zz, g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_z_x_xx[i] = 4.0 * g_y_yz_x_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_x_xy[i] = 4.0 * g_y_yz_x_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_x_xz[i] = 4.0 * g_y_yz_x_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_x_yy[i] = 4.0 * g_y_yz_x_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_x_yz[i] = 4.0 * g_y_yz_x_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_x_zz[i] = 4.0 * g_y_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_y_y_0_0_0_z_y_xx, g_y_y_0_0_0_z_y_xy, g_y_y_0_0_0_z_y_xz, g_y_y_0_0_0_z_y_yy, g_y_y_0_0_0_z_y_yz, g_y_y_0_0_0_z_y_zz, g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_z_y_xx[i] = 4.0 * g_y_yz_y_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_y_xy[i] = 4.0 * g_y_yz_y_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_y_xz[i] = 4.0 * g_y_yz_y_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_y_yy[i] = 4.0 * g_y_yz_y_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_y_yz[i] = 4.0 * g_y_yz_y_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_y_zz[i] = 4.0 * g_y_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_y_y_0_0_0_z_z_xx, g_y_y_0_0_0_z_z_xy, g_y_y_0_0_0_z_z_xz, g_y_y_0_0_0_z_z_yy, g_y_y_0_0_0_z_z_yz, g_y_y_0_0_0_z_z_zz, g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_z_z_xx[i] = 4.0 * g_y_yz_z_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_z_xy[i] = 4.0 * g_y_yz_z_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_z_xz[i] = 4.0 * g_y_yz_z_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_z_yy[i] = 4.0 * g_y_yz_z_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_z_yz[i] = 4.0 * g_y_yz_z_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_z_zz[i] = 4.0 * g_y_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_y_xz_x_xx, g_y_xz_x_xy, g_y_xz_x_xz, g_y_xz_x_yy, g_y_xz_x_yz, g_y_xz_x_zz, g_y_z_0_0_0_x_x_xx, g_y_z_0_0_0_x_x_xy, g_y_z_0_0_0_x_x_xz, g_y_z_0_0_0_x_x_yy, g_y_z_0_0_0_x_x_yz, g_y_z_0_0_0_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_x_x_xx[i] = 4.0 * g_y_xz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_x_xy[i] = 4.0 * g_y_xz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_x_xz[i] = 4.0 * g_y_xz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_x_yy[i] = 4.0 * g_y_xz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_x_yz[i] = 4.0 * g_y_xz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_x_zz[i] = 4.0 * g_y_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_y_xz_y_xx, g_y_xz_y_xy, g_y_xz_y_xz, g_y_xz_y_yy, g_y_xz_y_yz, g_y_xz_y_zz, g_y_z_0_0_0_x_y_xx, g_y_z_0_0_0_x_y_xy, g_y_z_0_0_0_x_y_xz, g_y_z_0_0_0_x_y_yy, g_y_z_0_0_0_x_y_yz, g_y_z_0_0_0_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_x_y_xx[i] = 4.0 * g_y_xz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_y_xy[i] = 4.0 * g_y_xz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_y_xz[i] = 4.0 * g_y_xz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_y_yy[i] = 4.0 * g_y_xz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_y_yz[i] = 4.0 * g_y_xz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_y_zz[i] = 4.0 * g_y_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_y_xz_z_xx, g_y_xz_z_xy, g_y_xz_z_xz, g_y_xz_z_yy, g_y_xz_z_yz, g_y_xz_z_zz, g_y_z_0_0_0_x_z_xx, g_y_z_0_0_0_x_z_xy, g_y_z_0_0_0_x_z_xz, g_y_z_0_0_0_x_z_yy, g_y_z_0_0_0_x_z_yz, g_y_z_0_0_0_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_x_z_xx[i] = 4.0 * g_y_xz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_z_xy[i] = 4.0 * g_y_xz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_z_xz[i] = 4.0 * g_y_xz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_z_yy[i] = 4.0 * g_y_xz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_z_yz[i] = 4.0 * g_y_xz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_z_zz[i] = 4.0 * g_y_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_y_yz_x_xx, g_y_yz_x_xy, g_y_yz_x_xz, g_y_yz_x_yy, g_y_yz_x_yz, g_y_yz_x_zz, g_y_z_0_0_0_y_x_xx, g_y_z_0_0_0_y_x_xy, g_y_z_0_0_0_y_x_xz, g_y_z_0_0_0_y_x_yy, g_y_z_0_0_0_y_x_yz, g_y_z_0_0_0_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_y_x_xx[i] = 4.0 * g_y_yz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_x_xy[i] = 4.0 * g_y_yz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_x_xz[i] = 4.0 * g_y_yz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_x_yy[i] = 4.0 * g_y_yz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_x_yz[i] = 4.0 * g_y_yz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_x_zz[i] = 4.0 * g_y_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_y_yz_y_xx, g_y_yz_y_xy, g_y_yz_y_xz, g_y_yz_y_yy, g_y_yz_y_yz, g_y_yz_y_zz, g_y_z_0_0_0_y_y_xx, g_y_z_0_0_0_y_y_xy, g_y_z_0_0_0_y_y_xz, g_y_z_0_0_0_y_y_yy, g_y_z_0_0_0_y_y_yz, g_y_z_0_0_0_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_y_y_xx[i] = 4.0 * g_y_yz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_y_xy[i] = 4.0 * g_y_yz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_y_xz[i] = 4.0 * g_y_yz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_y_yy[i] = 4.0 * g_y_yz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_y_yz[i] = 4.0 * g_y_yz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_y_zz[i] = 4.0 * g_y_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_y_yz_z_xx, g_y_yz_z_xy, g_y_yz_z_xz, g_y_yz_z_yy, g_y_yz_z_yz, g_y_yz_z_zz, g_y_z_0_0_0_y_z_xx, g_y_z_0_0_0_y_z_xy, g_y_z_0_0_0_y_z_xz, g_y_z_0_0_0_y_z_yy, g_y_z_0_0_0_y_z_yz, g_y_z_0_0_0_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_y_z_xx[i] = 4.0 * g_y_yz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_z_xy[i] = 4.0 * g_y_yz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_z_xz[i] = 4.0 * g_y_yz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_z_yy[i] = 4.0 * g_y_yz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_z_yz[i] = 4.0 * g_y_yz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_z_zz[i] = 4.0 * g_y_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_y_0_x_xx, g_y_0_x_xy, g_y_0_x_xz, g_y_0_x_yy, g_y_0_x_yz, g_y_0_x_zz, g_y_z_0_0_0_z_x_xx, g_y_z_0_0_0_z_x_xy, g_y_z_0_0_0_z_x_xz, g_y_z_0_0_0_z_x_yy, g_y_z_0_0_0_z_x_yz, g_y_z_0_0_0_z_x_zz, g_y_zz_x_xx, g_y_zz_x_xy, g_y_zz_x_xz, g_y_zz_x_yy, g_y_zz_x_yz, g_y_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_z_x_xx[i] = -2.0 * g_y_0_x_xx[i] * a_exp + 4.0 * g_y_zz_x_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_x_xy[i] = -2.0 * g_y_0_x_xy[i] * a_exp + 4.0 * g_y_zz_x_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_x_xz[i] = -2.0 * g_y_0_x_xz[i] * a_exp + 4.0 * g_y_zz_x_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_x_yy[i] = -2.0 * g_y_0_x_yy[i] * a_exp + 4.0 * g_y_zz_x_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_x_yz[i] = -2.0 * g_y_0_x_yz[i] * a_exp + 4.0 * g_y_zz_x_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_x_zz[i] = -2.0 * g_y_0_x_zz[i] * a_exp + 4.0 * g_y_zz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_y_0_y_xx, g_y_0_y_xy, g_y_0_y_xz, g_y_0_y_yy, g_y_0_y_yz, g_y_0_y_zz, g_y_z_0_0_0_z_y_xx, g_y_z_0_0_0_z_y_xy, g_y_z_0_0_0_z_y_xz, g_y_z_0_0_0_z_y_yy, g_y_z_0_0_0_z_y_yz, g_y_z_0_0_0_z_y_zz, g_y_zz_y_xx, g_y_zz_y_xy, g_y_zz_y_xz, g_y_zz_y_yy, g_y_zz_y_yz, g_y_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_z_y_xx[i] = -2.0 * g_y_0_y_xx[i] * a_exp + 4.0 * g_y_zz_y_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_y_xy[i] = -2.0 * g_y_0_y_xy[i] * a_exp + 4.0 * g_y_zz_y_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_y_xz[i] = -2.0 * g_y_0_y_xz[i] * a_exp + 4.0 * g_y_zz_y_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_y_yy[i] = -2.0 * g_y_0_y_yy[i] * a_exp + 4.0 * g_y_zz_y_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_y_yz[i] = -2.0 * g_y_0_y_yz[i] * a_exp + 4.0 * g_y_zz_y_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_y_zz[i] = -2.0 * g_y_0_y_zz[i] * a_exp + 4.0 * g_y_zz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_y_0_z_xx, g_y_0_z_xy, g_y_0_z_xz, g_y_0_z_yy, g_y_0_z_yz, g_y_0_z_zz, g_y_z_0_0_0_z_z_xx, g_y_z_0_0_0_z_z_xy, g_y_z_0_0_0_z_z_xz, g_y_z_0_0_0_z_z_yy, g_y_z_0_0_0_z_z_yz, g_y_z_0_0_0_z_z_zz, g_y_zz_z_xx, g_y_zz_z_xy, g_y_zz_z_xz, g_y_zz_z_yy, g_y_zz_z_yz, g_y_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_z_z_xx[i] = -2.0 * g_y_0_z_xx[i] * a_exp + 4.0 * g_y_zz_z_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_z_xy[i] = -2.0 * g_y_0_z_xy[i] * a_exp + 4.0 * g_y_zz_z_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_z_xz[i] = -2.0 * g_y_0_z_xz[i] * a_exp + 4.0 * g_y_zz_z_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_z_yy[i] = -2.0 * g_y_0_z_yy[i] * a_exp + 4.0 * g_y_zz_z_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_z_yz[i] = -2.0 * g_y_0_z_yz[i] * a_exp + 4.0 * g_y_zz_z_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_z_zz[i] = -2.0 * g_y_0_z_zz[i] * a_exp + 4.0 * g_y_zz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_z_0_x_xx, g_z_0_x_xy, g_z_0_x_xz, g_z_0_x_yy, g_z_0_x_yz, g_z_0_x_zz, g_z_x_0_0_0_x_x_xx, g_z_x_0_0_0_x_x_xy, g_z_x_0_0_0_x_x_xz, g_z_x_0_0_0_x_x_yy, g_z_x_0_0_0_x_x_yz, g_z_x_0_0_0_x_x_zz, g_z_xx_x_xx, g_z_xx_x_xy, g_z_xx_x_xz, g_z_xx_x_yy, g_z_xx_x_yz, g_z_xx_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_x_x_xx[i] = -2.0 * g_z_0_x_xx[i] * a_exp + 4.0 * g_z_xx_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_x_xy[i] = -2.0 * g_z_0_x_xy[i] * a_exp + 4.0 * g_z_xx_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_x_xz[i] = -2.0 * g_z_0_x_xz[i] * a_exp + 4.0 * g_z_xx_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_x_yy[i] = -2.0 * g_z_0_x_yy[i] * a_exp + 4.0 * g_z_xx_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_x_yz[i] = -2.0 * g_z_0_x_yz[i] * a_exp + 4.0 * g_z_xx_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_x_zz[i] = -2.0 * g_z_0_x_zz[i] * a_exp + 4.0 * g_z_xx_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_z_0_y_xx, g_z_0_y_xy, g_z_0_y_xz, g_z_0_y_yy, g_z_0_y_yz, g_z_0_y_zz, g_z_x_0_0_0_x_y_xx, g_z_x_0_0_0_x_y_xy, g_z_x_0_0_0_x_y_xz, g_z_x_0_0_0_x_y_yy, g_z_x_0_0_0_x_y_yz, g_z_x_0_0_0_x_y_zz, g_z_xx_y_xx, g_z_xx_y_xy, g_z_xx_y_xz, g_z_xx_y_yy, g_z_xx_y_yz, g_z_xx_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_x_y_xx[i] = -2.0 * g_z_0_y_xx[i] * a_exp + 4.0 * g_z_xx_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_y_xy[i] = -2.0 * g_z_0_y_xy[i] * a_exp + 4.0 * g_z_xx_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_y_xz[i] = -2.0 * g_z_0_y_xz[i] * a_exp + 4.0 * g_z_xx_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_y_yy[i] = -2.0 * g_z_0_y_yy[i] * a_exp + 4.0 * g_z_xx_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_y_yz[i] = -2.0 * g_z_0_y_yz[i] * a_exp + 4.0 * g_z_xx_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_y_zz[i] = -2.0 * g_z_0_y_zz[i] * a_exp + 4.0 * g_z_xx_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_z_0_z_xx, g_z_0_z_xy, g_z_0_z_xz, g_z_0_z_yy, g_z_0_z_yz, g_z_0_z_zz, g_z_x_0_0_0_x_z_xx, g_z_x_0_0_0_x_z_xy, g_z_x_0_0_0_x_z_xz, g_z_x_0_0_0_x_z_yy, g_z_x_0_0_0_x_z_yz, g_z_x_0_0_0_x_z_zz, g_z_xx_z_xx, g_z_xx_z_xy, g_z_xx_z_xz, g_z_xx_z_yy, g_z_xx_z_yz, g_z_xx_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_x_z_xx[i] = -2.0 * g_z_0_z_xx[i] * a_exp + 4.0 * g_z_xx_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_z_xy[i] = -2.0 * g_z_0_z_xy[i] * a_exp + 4.0 * g_z_xx_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_z_xz[i] = -2.0 * g_z_0_z_xz[i] * a_exp + 4.0 * g_z_xx_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_z_yy[i] = -2.0 * g_z_0_z_yy[i] * a_exp + 4.0 * g_z_xx_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_z_yz[i] = -2.0 * g_z_0_z_yz[i] * a_exp + 4.0 * g_z_xx_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_z_zz[i] = -2.0 * g_z_0_z_zz[i] * a_exp + 4.0 * g_z_xx_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_z_x_0_0_0_y_x_xx, g_z_x_0_0_0_y_x_xy, g_z_x_0_0_0_y_x_xz, g_z_x_0_0_0_y_x_yy, g_z_x_0_0_0_y_x_yz, g_z_x_0_0_0_y_x_zz, g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_y_x_xx[i] = 4.0 * g_z_xy_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_x_xy[i] = 4.0 * g_z_xy_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_x_xz[i] = 4.0 * g_z_xy_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_x_yy[i] = 4.0 * g_z_xy_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_x_yz[i] = 4.0 * g_z_xy_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_x_zz[i] = 4.0 * g_z_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_z_x_0_0_0_y_y_xx, g_z_x_0_0_0_y_y_xy, g_z_x_0_0_0_y_y_xz, g_z_x_0_0_0_y_y_yy, g_z_x_0_0_0_y_y_yz, g_z_x_0_0_0_y_y_zz, g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_y_y_xx[i] = 4.0 * g_z_xy_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_y_xy[i] = 4.0 * g_z_xy_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_y_xz[i] = 4.0 * g_z_xy_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_y_yy[i] = 4.0 * g_z_xy_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_y_yz[i] = 4.0 * g_z_xy_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_y_zz[i] = 4.0 * g_z_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_z_x_0_0_0_y_z_xx, g_z_x_0_0_0_y_z_xy, g_z_x_0_0_0_y_z_xz, g_z_x_0_0_0_y_z_yy, g_z_x_0_0_0_y_z_yz, g_z_x_0_0_0_y_z_zz, g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_y_z_xx[i] = 4.0 * g_z_xy_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_z_xy[i] = 4.0 * g_z_xy_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_z_xz[i] = 4.0 * g_z_xy_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_z_yy[i] = 4.0 * g_z_xy_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_z_yz[i] = 4.0 * g_z_xy_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_z_zz[i] = 4.0 * g_z_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_z_x_0_0_0_z_x_xx, g_z_x_0_0_0_z_x_xy, g_z_x_0_0_0_z_x_xz, g_z_x_0_0_0_z_x_yy, g_z_x_0_0_0_z_x_yz, g_z_x_0_0_0_z_x_zz, g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_z_x_xx[i] = 4.0 * g_z_xz_x_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_x_xy[i] = 4.0 * g_z_xz_x_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_x_xz[i] = 4.0 * g_z_xz_x_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_x_yy[i] = 4.0 * g_z_xz_x_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_x_yz[i] = 4.0 * g_z_xz_x_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_x_zz[i] = 4.0 * g_z_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_z_x_0_0_0_z_y_xx, g_z_x_0_0_0_z_y_xy, g_z_x_0_0_0_z_y_xz, g_z_x_0_0_0_z_y_yy, g_z_x_0_0_0_z_y_yz, g_z_x_0_0_0_z_y_zz, g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_z_y_xx[i] = 4.0 * g_z_xz_y_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_y_xy[i] = 4.0 * g_z_xz_y_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_y_xz[i] = 4.0 * g_z_xz_y_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_y_yy[i] = 4.0 * g_z_xz_y_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_y_yz[i] = 4.0 * g_z_xz_y_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_y_zz[i] = 4.0 * g_z_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_z_x_0_0_0_z_z_xx, g_z_x_0_0_0_z_z_xy, g_z_x_0_0_0_z_z_xz, g_z_x_0_0_0_z_z_yy, g_z_x_0_0_0_z_z_yz, g_z_x_0_0_0_z_z_zz, g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_z_z_xx[i] = 4.0 * g_z_xz_z_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_z_xy[i] = 4.0 * g_z_xz_z_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_z_xz[i] = 4.0 * g_z_xz_z_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_z_yy[i] = 4.0 * g_z_xz_z_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_z_yz[i] = 4.0 * g_z_xz_z_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_z_zz[i] = 4.0 * g_z_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_z_xy_x_xx, g_z_xy_x_xy, g_z_xy_x_xz, g_z_xy_x_yy, g_z_xy_x_yz, g_z_xy_x_zz, g_z_y_0_0_0_x_x_xx, g_z_y_0_0_0_x_x_xy, g_z_y_0_0_0_x_x_xz, g_z_y_0_0_0_x_x_yy, g_z_y_0_0_0_x_x_yz, g_z_y_0_0_0_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_x_x_xx[i] = 4.0 * g_z_xy_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_x_xy[i] = 4.0 * g_z_xy_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_x_xz[i] = 4.0 * g_z_xy_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_x_yy[i] = 4.0 * g_z_xy_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_x_yz[i] = 4.0 * g_z_xy_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_x_zz[i] = 4.0 * g_z_xy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_z_xy_y_xx, g_z_xy_y_xy, g_z_xy_y_xz, g_z_xy_y_yy, g_z_xy_y_yz, g_z_xy_y_zz, g_z_y_0_0_0_x_y_xx, g_z_y_0_0_0_x_y_xy, g_z_y_0_0_0_x_y_xz, g_z_y_0_0_0_x_y_yy, g_z_y_0_0_0_x_y_yz, g_z_y_0_0_0_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_x_y_xx[i] = 4.0 * g_z_xy_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_y_xy[i] = 4.0 * g_z_xy_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_y_xz[i] = 4.0 * g_z_xy_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_y_yy[i] = 4.0 * g_z_xy_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_y_yz[i] = 4.0 * g_z_xy_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_y_zz[i] = 4.0 * g_z_xy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_z_xy_z_xx, g_z_xy_z_xy, g_z_xy_z_xz, g_z_xy_z_yy, g_z_xy_z_yz, g_z_xy_z_zz, g_z_y_0_0_0_x_z_xx, g_z_y_0_0_0_x_z_xy, g_z_y_0_0_0_x_z_xz, g_z_y_0_0_0_x_z_yy, g_z_y_0_0_0_x_z_yz, g_z_y_0_0_0_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_x_z_xx[i] = 4.0 * g_z_xy_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_z_xy[i] = 4.0 * g_z_xy_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_z_xz[i] = 4.0 * g_z_xy_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_z_yy[i] = 4.0 * g_z_xy_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_z_yz[i] = 4.0 * g_z_xy_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_z_zz[i] = 4.0 * g_z_xy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_z_0_x_xx, g_z_0_x_xy, g_z_0_x_xz, g_z_0_x_yy, g_z_0_x_yz, g_z_0_x_zz, g_z_y_0_0_0_y_x_xx, g_z_y_0_0_0_y_x_xy, g_z_y_0_0_0_y_x_xz, g_z_y_0_0_0_y_x_yy, g_z_y_0_0_0_y_x_yz, g_z_y_0_0_0_y_x_zz, g_z_yy_x_xx, g_z_yy_x_xy, g_z_yy_x_xz, g_z_yy_x_yy, g_z_yy_x_yz, g_z_yy_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_y_x_xx[i] = -2.0 * g_z_0_x_xx[i] * a_exp + 4.0 * g_z_yy_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_x_xy[i] = -2.0 * g_z_0_x_xy[i] * a_exp + 4.0 * g_z_yy_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_x_xz[i] = -2.0 * g_z_0_x_xz[i] * a_exp + 4.0 * g_z_yy_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_x_yy[i] = -2.0 * g_z_0_x_yy[i] * a_exp + 4.0 * g_z_yy_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_x_yz[i] = -2.0 * g_z_0_x_yz[i] * a_exp + 4.0 * g_z_yy_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_x_zz[i] = -2.0 * g_z_0_x_zz[i] * a_exp + 4.0 * g_z_yy_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_z_0_y_xx, g_z_0_y_xy, g_z_0_y_xz, g_z_0_y_yy, g_z_0_y_yz, g_z_0_y_zz, g_z_y_0_0_0_y_y_xx, g_z_y_0_0_0_y_y_xy, g_z_y_0_0_0_y_y_xz, g_z_y_0_0_0_y_y_yy, g_z_y_0_0_0_y_y_yz, g_z_y_0_0_0_y_y_zz, g_z_yy_y_xx, g_z_yy_y_xy, g_z_yy_y_xz, g_z_yy_y_yy, g_z_yy_y_yz, g_z_yy_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_y_y_xx[i] = -2.0 * g_z_0_y_xx[i] * a_exp + 4.0 * g_z_yy_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_y_xy[i] = -2.0 * g_z_0_y_xy[i] * a_exp + 4.0 * g_z_yy_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_y_xz[i] = -2.0 * g_z_0_y_xz[i] * a_exp + 4.0 * g_z_yy_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_y_yy[i] = -2.0 * g_z_0_y_yy[i] * a_exp + 4.0 * g_z_yy_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_y_yz[i] = -2.0 * g_z_0_y_yz[i] * a_exp + 4.0 * g_z_yy_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_y_zz[i] = -2.0 * g_z_0_y_zz[i] * a_exp + 4.0 * g_z_yy_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_z_0_z_xx, g_z_0_z_xy, g_z_0_z_xz, g_z_0_z_yy, g_z_0_z_yz, g_z_0_z_zz, g_z_y_0_0_0_y_z_xx, g_z_y_0_0_0_y_z_xy, g_z_y_0_0_0_y_z_xz, g_z_y_0_0_0_y_z_yy, g_z_y_0_0_0_y_z_yz, g_z_y_0_0_0_y_z_zz, g_z_yy_z_xx, g_z_yy_z_xy, g_z_yy_z_xz, g_z_yy_z_yy, g_z_yy_z_yz, g_z_yy_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_y_z_xx[i] = -2.0 * g_z_0_z_xx[i] * a_exp + 4.0 * g_z_yy_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_z_xy[i] = -2.0 * g_z_0_z_xy[i] * a_exp + 4.0 * g_z_yy_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_z_xz[i] = -2.0 * g_z_0_z_xz[i] * a_exp + 4.0 * g_z_yy_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_z_yy[i] = -2.0 * g_z_0_z_yy[i] * a_exp + 4.0 * g_z_yy_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_z_yz[i] = -2.0 * g_z_0_z_yz[i] * a_exp + 4.0 * g_z_yy_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_z_zz[i] = -2.0 * g_z_0_z_zz[i] * a_exp + 4.0 * g_z_yy_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_z_y_0_0_0_z_x_xx, g_z_y_0_0_0_z_x_xy, g_z_y_0_0_0_z_x_xz, g_z_y_0_0_0_z_x_yy, g_z_y_0_0_0_z_x_yz, g_z_y_0_0_0_z_x_zz, g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_z_x_xx[i] = 4.0 * g_z_yz_x_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_x_xy[i] = 4.0 * g_z_yz_x_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_x_xz[i] = 4.0 * g_z_yz_x_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_x_yy[i] = 4.0 * g_z_yz_x_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_x_yz[i] = 4.0 * g_z_yz_x_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_x_zz[i] = 4.0 * g_z_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_z_y_0_0_0_z_y_xx, g_z_y_0_0_0_z_y_xy, g_z_y_0_0_0_z_y_xz, g_z_y_0_0_0_z_y_yy, g_z_y_0_0_0_z_y_yz, g_z_y_0_0_0_z_y_zz, g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_z_y_xx[i] = 4.0 * g_z_yz_y_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_y_xy[i] = 4.0 * g_z_yz_y_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_y_xz[i] = 4.0 * g_z_yz_y_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_y_yy[i] = 4.0 * g_z_yz_y_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_y_yz[i] = 4.0 * g_z_yz_y_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_y_zz[i] = 4.0 * g_z_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_z_y_0_0_0_z_z_xx, g_z_y_0_0_0_z_z_xy, g_z_y_0_0_0_z_z_xz, g_z_y_0_0_0_z_z_yy, g_z_y_0_0_0_z_z_yz, g_z_y_0_0_0_z_z_zz, g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_z_z_xx[i] = 4.0 * g_z_yz_z_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_z_xy[i] = 4.0 * g_z_yz_z_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_z_xz[i] = 4.0 * g_z_yz_z_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_z_yy[i] = 4.0 * g_z_yz_z_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_z_yz[i] = 4.0 * g_z_yz_z_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_z_zz[i] = 4.0 * g_z_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_z_xz_x_xx, g_z_xz_x_xy, g_z_xz_x_xz, g_z_xz_x_yy, g_z_xz_x_yz, g_z_xz_x_zz, g_z_z_0_0_0_x_x_xx, g_z_z_0_0_0_x_x_xy, g_z_z_0_0_0_x_x_xz, g_z_z_0_0_0_x_x_yy, g_z_z_0_0_0_x_x_yz, g_z_z_0_0_0_x_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_x_x_xx[i] = 4.0 * g_z_xz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_x_xy[i] = 4.0 * g_z_xz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_x_xz[i] = 4.0 * g_z_xz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_x_yy[i] = 4.0 * g_z_xz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_x_yz[i] = 4.0 * g_z_xz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_x_zz[i] = 4.0 * g_z_xz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_z_xz_y_xx, g_z_xz_y_xy, g_z_xz_y_xz, g_z_xz_y_yy, g_z_xz_y_yz, g_z_xz_y_zz, g_z_z_0_0_0_x_y_xx, g_z_z_0_0_0_x_y_xy, g_z_z_0_0_0_x_y_xz, g_z_z_0_0_0_x_y_yy, g_z_z_0_0_0_x_y_yz, g_z_z_0_0_0_x_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_x_y_xx[i] = 4.0 * g_z_xz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_y_xy[i] = 4.0 * g_z_xz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_y_xz[i] = 4.0 * g_z_xz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_y_yy[i] = 4.0 * g_z_xz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_y_yz[i] = 4.0 * g_z_xz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_y_zz[i] = 4.0 * g_z_xz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_z_xz_z_xx, g_z_xz_z_xy, g_z_xz_z_xz, g_z_xz_z_yy, g_z_xz_z_yz, g_z_xz_z_zz, g_z_z_0_0_0_x_z_xx, g_z_z_0_0_0_x_z_xy, g_z_z_0_0_0_x_z_xz, g_z_z_0_0_0_x_z_yy, g_z_z_0_0_0_x_z_yz, g_z_z_0_0_0_x_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_x_z_xx[i] = 4.0 * g_z_xz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_z_xy[i] = 4.0 * g_z_xz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_z_xz[i] = 4.0 * g_z_xz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_z_yy[i] = 4.0 * g_z_xz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_z_yz[i] = 4.0 * g_z_xz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_z_zz[i] = 4.0 * g_z_xz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_z_yz_x_xx, g_z_yz_x_xy, g_z_yz_x_xz, g_z_yz_x_yy, g_z_yz_x_yz, g_z_yz_x_zz, g_z_z_0_0_0_y_x_xx, g_z_z_0_0_0_y_x_xy, g_z_z_0_0_0_y_x_xz, g_z_z_0_0_0_y_x_yy, g_z_z_0_0_0_y_x_yz, g_z_z_0_0_0_y_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_y_x_xx[i] = 4.0 * g_z_yz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_x_xy[i] = 4.0 * g_z_yz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_x_xz[i] = 4.0 * g_z_yz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_x_yy[i] = 4.0 * g_z_yz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_x_yz[i] = 4.0 * g_z_yz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_x_zz[i] = 4.0 * g_z_yz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_z_yz_y_xx, g_z_yz_y_xy, g_z_yz_y_xz, g_z_yz_y_yy, g_z_yz_y_yz, g_z_yz_y_zz, g_z_z_0_0_0_y_y_xx, g_z_z_0_0_0_y_y_xy, g_z_z_0_0_0_y_y_xz, g_z_z_0_0_0_y_y_yy, g_z_z_0_0_0_y_y_yz, g_z_z_0_0_0_y_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_y_y_xx[i] = 4.0 * g_z_yz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_y_xy[i] = 4.0 * g_z_yz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_y_xz[i] = 4.0 * g_z_yz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_y_yy[i] = 4.0 * g_z_yz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_y_yz[i] = 4.0 * g_z_yz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_y_zz[i] = 4.0 * g_z_yz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_z_yz_z_xx, g_z_yz_z_xy, g_z_yz_z_xz, g_z_yz_z_yy, g_z_yz_z_yz, g_z_yz_z_zz, g_z_z_0_0_0_y_z_xx, g_z_z_0_0_0_y_z_xy, g_z_z_0_0_0_y_z_xz, g_z_z_0_0_0_y_z_yy, g_z_z_0_0_0_y_z_yz, g_z_z_0_0_0_y_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_y_z_xx[i] = 4.0 * g_z_yz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_z_xy[i] = 4.0 * g_z_yz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_z_xz[i] = 4.0 * g_z_yz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_z_yy[i] = 4.0 * g_z_yz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_z_yz[i] = 4.0 * g_z_yz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_z_zz[i] = 4.0 * g_z_yz_z_zz[i] * a_exp * b_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_z_0_x_xx, g_z_0_x_xy, g_z_0_x_xz, g_z_0_x_yy, g_z_0_x_yz, g_z_0_x_zz, g_z_z_0_0_0_z_x_xx, g_z_z_0_0_0_z_x_xy, g_z_z_0_0_0_z_x_xz, g_z_z_0_0_0_z_x_yy, g_z_z_0_0_0_z_x_yz, g_z_z_0_0_0_z_x_zz, g_z_zz_x_xx, g_z_zz_x_xy, g_z_zz_x_xz, g_z_zz_x_yy, g_z_zz_x_yz, g_z_zz_x_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_z_x_xx[i] = -2.0 * g_z_0_x_xx[i] * a_exp + 4.0 * g_z_zz_x_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_x_xy[i] = -2.0 * g_z_0_x_xy[i] * a_exp + 4.0 * g_z_zz_x_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_x_xz[i] = -2.0 * g_z_0_x_xz[i] * a_exp + 4.0 * g_z_zz_x_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_x_yy[i] = -2.0 * g_z_0_x_yy[i] * a_exp + 4.0 * g_z_zz_x_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_x_yz[i] = -2.0 * g_z_0_x_yz[i] * a_exp + 4.0 * g_z_zz_x_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_x_zz[i] = -2.0 * g_z_0_x_zz[i] * a_exp + 4.0 * g_z_zz_x_zz[i] * a_exp * b_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_z_0_y_xx, g_z_0_y_xy, g_z_0_y_xz, g_z_0_y_yy, g_z_0_y_yz, g_z_0_y_zz, g_z_z_0_0_0_z_y_xx, g_z_z_0_0_0_z_y_xy, g_z_z_0_0_0_z_y_xz, g_z_z_0_0_0_z_y_yy, g_z_z_0_0_0_z_y_yz, g_z_z_0_0_0_z_y_zz, g_z_zz_y_xx, g_z_zz_y_xy, g_z_zz_y_xz, g_z_zz_y_yy, g_z_zz_y_yz, g_z_zz_y_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_z_y_xx[i] = -2.0 * g_z_0_y_xx[i] * a_exp + 4.0 * g_z_zz_y_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_y_xy[i] = -2.0 * g_z_0_y_xy[i] * a_exp + 4.0 * g_z_zz_y_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_y_xz[i] = -2.0 * g_z_0_y_xz[i] * a_exp + 4.0 * g_z_zz_y_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_y_yy[i] = -2.0 * g_z_0_y_yy[i] * a_exp + 4.0 * g_z_zz_y_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_y_yz[i] = -2.0 * g_z_0_y_yz[i] * a_exp + 4.0 * g_z_zz_y_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_y_zz[i] = -2.0 * g_z_0_y_zz[i] * a_exp + 4.0 * g_z_zz_y_zz[i] * a_exp * b_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_z_0_z_xx, g_z_0_z_xy, g_z_0_z_xz, g_z_0_z_yy, g_z_0_z_yz, g_z_0_z_zz, g_z_z_0_0_0_z_z_xx, g_z_z_0_0_0_z_z_xy, g_z_z_0_0_0_z_z_xz, g_z_z_0_0_0_z_z_yy, g_z_z_0_0_0_z_z_yz, g_z_z_0_0_0_z_z_zz, g_z_zz_z_xx, g_z_zz_z_xy, g_z_zz_z_xz, g_z_zz_z_yy, g_z_zz_z_yz, g_z_zz_z_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_z_z_xx[i] = -2.0 * g_z_0_z_xx[i] * a_exp + 4.0 * g_z_zz_z_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_z_xy[i] = -2.0 * g_z_0_z_xy[i] * a_exp + 4.0 * g_z_zz_z_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_z_xz[i] = -2.0 * g_z_0_z_xz[i] * a_exp + 4.0 * g_z_zz_z_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_z_yy[i] = -2.0 * g_z_0_z_yy[i] * a_exp + 4.0 * g_z_zz_z_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_z_yz[i] = -2.0 * g_z_0_z_yz[i] * a_exp + 4.0 * g_z_zz_z_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_z_zz[i] = -2.0 * g_z_0_z_zz[i] * a_exp + 4.0 * g_z_zz_z_zz[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

