#include "GeomDeriv1010OfScalarForSDPD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_sdpd_0(CSimdArray<double>& buffer_1010_sdpd,
                     const CSimdArray<double>& buffer_pdsd,
                     const CSimdArray<double>& buffer_pddd,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_sdpd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pdsd

    auto g_x_xx_0_xx = buffer_pdsd[0];

    auto g_x_xx_0_xy = buffer_pdsd[1];

    auto g_x_xx_0_xz = buffer_pdsd[2];

    auto g_x_xx_0_yy = buffer_pdsd[3];

    auto g_x_xx_0_yz = buffer_pdsd[4];

    auto g_x_xx_0_zz = buffer_pdsd[5];

    auto g_x_xy_0_xx = buffer_pdsd[6];

    auto g_x_xy_0_xy = buffer_pdsd[7];

    auto g_x_xy_0_xz = buffer_pdsd[8];

    auto g_x_xy_0_yy = buffer_pdsd[9];

    auto g_x_xy_0_yz = buffer_pdsd[10];

    auto g_x_xy_0_zz = buffer_pdsd[11];

    auto g_x_xz_0_xx = buffer_pdsd[12];

    auto g_x_xz_0_xy = buffer_pdsd[13];

    auto g_x_xz_0_xz = buffer_pdsd[14];

    auto g_x_xz_0_yy = buffer_pdsd[15];

    auto g_x_xz_0_yz = buffer_pdsd[16];

    auto g_x_xz_0_zz = buffer_pdsd[17];

    auto g_x_yy_0_xx = buffer_pdsd[18];

    auto g_x_yy_0_xy = buffer_pdsd[19];

    auto g_x_yy_0_xz = buffer_pdsd[20];

    auto g_x_yy_0_yy = buffer_pdsd[21];

    auto g_x_yy_0_yz = buffer_pdsd[22];

    auto g_x_yy_0_zz = buffer_pdsd[23];

    auto g_x_yz_0_xx = buffer_pdsd[24];

    auto g_x_yz_0_xy = buffer_pdsd[25];

    auto g_x_yz_0_xz = buffer_pdsd[26];

    auto g_x_yz_0_yy = buffer_pdsd[27];

    auto g_x_yz_0_yz = buffer_pdsd[28];

    auto g_x_yz_0_zz = buffer_pdsd[29];

    auto g_x_zz_0_xx = buffer_pdsd[30];

    auto g_x_zz_0_xy = buffer_pdsd[31];

    auto g_x_zz_0_xz = buffer_pdsd[32];

    auto g_x_zz_0_yy = buffer_pdsd[33];

    auto g_x_zz_0_yz = buffer_pdsd[34];

    auto g_x_zz_0_zz = buffer_pdsd[35];

    auto g_y_xx_0_xx = buffer_pdsd[36];

    auto g_y_xx_0_xy = buffer_pdsd[37];

    auto g_y_xx_0_xz = buffer_pdsd[38];

    auto g_y_xx_0_yy = buffer_pdsd[39];

    auto g_y_xx_0_yz = buffer_pdsd[40];

    auto g_y_xx_0_zz = buffer_pdsd[41];

    auto g_y_xy_0_xx = buffer_pdsd[42];

    auto g_y_xy_0_xy = buffer_pdsd[43];

    auto g_y_xy_0_xz = buffer_pdsd[44];

    auto g_y_xy_0_yy = buffer_pdsd[45];

    auto g_y_xy_0_yz = buffer_pdsd[46];

    auto g_y_xy_0_zz = buffer_pdsd[47];

    auto g_y_xz_0_xx = buffer_pdsd[48];

    auto g_y_xz_0_xy = buffer_pdsd[49];

    auto g_y_xz_0_xz = buffer_pdsd[50];

    auto g_y_xz_0_yy = buffer_pdsd[51];

    auto g_y_xz_0_yz = buffer_pdsd[52];

    auto g_y_xz_0_zz = buffer_pdsd[53];

    auto g_y_yy_0_xx = buffer_pdsd[54];

    auto g_y_yy_0_xy = buffer_pdsd[55];

    auto g_y_yy_0_xz = buffer_pdsd[56];

    auto g_y_yy_0_yy = buffer_pdsd[57];

    auto g_y_yy_0_yz = buffer_pdsd[58];

    auto g_y_yy_0_zz = buffer_pdsd[59];

    auto g_y_yz_0_xx = buffer_pdsd[60];

    auto g_y_yz_0_xy = buffer_pdsd[61];

    auto g_y_yz_0_xz = buffer_pdsd[62];

    auto g_y_yz_0_yy = buffer_pdsd[63];

    auto g_y_yz_0_yz = buffer_pdsd[64];

    auto g_y_yz_0_zz = buffer_pdsd[65];

    auto g_y_zz_0_xx = buffer_pdsd[66];

    auto g_y_zz_0_xy = buffer_pdsd[67];

    auto g_y_zz_0_xz = buffer_pdsd[68];

    auto g_y_zz_0_yy = buffer_pdsd[69];

    auto g_y_zz_0_yz = buffer_pdsd[70];

    auto g_y_zz_0_zz = buffer_pdsd[71];

    auto g_z_xx_0_xx = buffer_pdsd[72];

    auto g_z_xx_0_xy = buffer_pdsd[73];

    auto g_z_xx_0_xz = buffer_pdsd[74];

    auto g_z_xx_0_yy = buffer_pdsd[75];

    auto g_z_xx_0_yz = buffer_pdsd[76];

    auto g_z_xx_0_zz = buffer_pdsd[77];

    auto g_z_xy_0_xx = buffer_pdsd[78];

    auto g_z_xy_0_xy = buffer_pdsd[79];

    auto g_z_xy_0_xz = buffer_pdsd[80];

    auto g_z_xy_0_yy = buffer_pdsd[81];

    auto g_z_xy_0_yz = buffer_pdsd[82];

    auto g_z_xy_0_zz = buffer_pdsd[83];

    auto g_z_xz_0_xx = buffer_pdsd[84];

    auto g_z_xz_0_xy = buffer_pdsd[85];

    auto g_z_xz_0_xz = buffer_pdsd[86];

    auto g_z_xz_0_yy = buffer_pdsd[87];

    auto g_z_xz_0_yz = buffer_pdsd[88];

    auto g_z_xz_0_zz = buffer_pdsd[89];

    auto g_z_yy_0_xx = buffer_pdsd[90];

    auto g_z_yy_0_xy = buffer_pdsd[91];

    auto g_z_yy_0_xz = buffer_pdsd[92];

    auto g_z_yy_0_yy = buffer_pdsd[93];

    auto g_z_yy_0_yz = buffer_pdsd[94];

    auto g_z_yy_0_zz = buffer_pdsd[95];

    auto g_z_yz_0_xx = buffer_pdsd[96];

    auto g_z_yz_0_xy = buffer_pdsd[97];

    auto g_z_yz_0_xz = buffer_pdsd[98];

    auto g_z_yz_0_yy = buffer_pdsd[99];

    auto g_z_yz_0_yz = buffer_pdsd[100];

    auto g_z_yz_0_zz = buffer_pdsd[101];

    auto g_z_zz_0_xx = buffer_pdsd[102];

    auto g_z_zz_0_xy = buffer_pdsd[103];

    auto g_z_zz_0_xz = buffer_pdsd[104];

    auto g_z_zz_0_yy = buffer_pdsd[105];

    auto g_z_zz_0_yz = buffer_pdsd[106];

    auto g_z_zz_0_zz = buffer_pdsd[107];

    /// Set up components of auxilary buffer : buffer_pddd

    auto g_x_xx_xx_xx = buffer_pddd[0];

    auto g_x_xx_xx_xy = buffer_pddd[1];

    auto g_x_xx_xx_xz = buffer_pddd[2];

    auto g_x_xx_xx_yy = buffer_pddd[3];

    auto g_x_xx_xx_yz = buffer_pddd[4];

    auto g_x_xx_xx_zz = buffer_pddd[5];

    auto g_x_xx_xy_xx = buffer_pddd[6];

    auto g_x_xx_xy_xy = buffer_pddd[7];

    auto g_x_xx_xy_xz = buffer_pddd[8];

    auto g_x_xx_xy_yy = buffer_pddd[9];

    auto g_x_xx_xy_yz = buffer_pddd[10];

    auto g_x_xx_xy_zz = buffer_pddd[11];

    auto g_x_xx_xz_xx = buffer_pddd[12];

    auto g_x_xx_xz_xy = buffer_pddd[13];

    auto g_x_xx_xz_xz = buffer_pddd[14];

    auto g_x_xx_xz_yy = buffer_pddd[15];

    auto g_x_xx_xz_yz = buffer_pddd[16];

    auto g_x_xx_xz_zz = buffer_pddd[17];

    auto g_x_xx_yy_xx = buffer_pddd[18];

    auto g_x_xx_yy_xy = buffer_pddd[19];

    auto g_x_xx_yy_xz = buffer_pddd[20];

    auto g_x_xx_yy_yy = buffer_pddd[21];

    auto g_x_xx_yy_yz = buffer_pddd[22];

    auto g_x_xx_yy_zz = buffer_pddd[23];

    auto g_x_xx_yz_xx = buffer_pddd[24];

    auto g_x_xx_yz_xy = buffer_pddd[25];

    auto g_x_xx_yz_xz = buffer_pddd[26];

    auto g_x_xx_yz_yy = buffer_pddd[27];

    auto g_x_xx_yz_yz = buffer_pddd[28];

    auto g_x_xx_yz_zz = buffer_pddd[29];

    auto g_x_xx_zz_xx = buffer_pddd[30];

    auto g_x_xx_zz_xy = buffer_pddd[31];

    auto g_x_xx_zz_xz = buffer_pddd[32];

    auto g_x_xx_zz_yy = buffer_pddd[33];

    auto g_x_xx_zz_yz = buffer_pddd[34];

    auto g_x_xx_zz_zz = buffer_pddd[35];

    auto g_x_xy_xx_xx = buffer_pddd[36];

    auto g_x_xy_xx_xy = buffer_pddd[37];

    auto g_x_xy_xx_xz = buffer_pddd[38];

    auto g_x_xy_xx_yy = buffer_pddd[39];

    auto g_x_xy_xx_yz = buffer_pddd[40];

    auto g_x_xy_xx_zz = buffer_pddd[41];

    auto g_x_xy_xy_xx = buffer_pddd[42];

    auto g_x_xy_xy_xy = buffer_pddd[43];

    auto g_x_xy_xy_xz = buffer_pddd[44];

    auto g_x_xy_xy_yy = buffer_pddd[45];

    auto g_x_xy_xy_yz = buffer_pddd[46];

    auto g_x_xy_xy_zz = buffer_pddd[47];

    auto g_x_xy_xz_xx = buffer_pddd[48];

    auto g_x_xy_xz_xy = buffer_pddd[49];

    auto g_x_xy_xz_xz = buffer_pddd[50];

    auto g_x_xy_xz_yy = buffer_pddd[51];

    auto g_x_xy_xz_yz = buffer_pddd[52];

    auto g_x_xy_xz_zz = buffer_pddd[53];

    auto g_x_xy_yy_xx = buffer_pddd[54];

    auto g_x_xy_yy_xy = buffer_pddd[55];

    auto g_x_xy_yy_xz = buffer_pddd[56];

    auto g_x_xy_yy_yy = buffer_pddd[57];

    auto g_x_xy_yy_yz = buffer_pddd[58];

    auto g_x_xy_yy_zz = buffer_pddd[59];

    auto g_x_xy_yz_xx = buffer_pddd[60];

    auto g_x_xy_yz_xy = buffer_pddd[61];

    auto g_x_xy_yz_xz = buffer_pddd[62];

    auto g_x_xy_yz_yy = buffer_pddd[63];

    auto g_x_xy_yz_yz = buffer_pddd[64];

    auto g_x_xy_yz_zz = buffer_pddd[65];

    auto g_x_xy_zz_xx = buffer_pddd[66];

    auto g_x_xy_zz_xy = buffer_pddd[67];

    auto g_x_xy_zz_xz = buffer_pddd[68];

    auto g_x_xy_zz_yy = buffer_pddd[69];

    auto g_x_xy_zz_yz = buffer_pddd[70];

    auto g_x_xy_zz_zz = buffer_pddd[71];

    auto g_x_xz_xx_xx = buffer_pddd[72];

    auto g_x_xz_xx_xy = buffer_pddd[73];

    auto g_x_xz_xx_xz = buffer_pddd[74];

    auto g_x_xz_xx_yy = buffer_pddd[75];

    auto g_x_xz_xx_yz = buffer_pddd[76];

    auto g_x_xz_xx_zz = buffer_pddd[77];

    auto g_x_xz_xy_xx = buffer_pddd[78];

    auto g_x_xz_xy_xy = buffer_pddd[79];

    auto g_x_xz_xy_xz = buffer_pddd[80];

    auto g_x_xz_xy_yy = buffer_pddd[81];

    auto g_x_xz_xy_yz = buffer_pddd[82];

    auto g_x_xz_xy_zz = buffer_pddd[83];

    auto g_x_xz_xz_xx = buffer_pddd[84];

    auto g_x_xz_xz_xy = buffer_pddd[85];

    auto g_x_xz_xz_xz = buffer_pddd[86];

    auto g_x_xz_xz_yy = buffer_pddd[87];

    auto g_x_xz_xz_yz = buffer_pddd[88];

    auto g_x_xz_xz_zz = buffer_pddd[89];

    auto g_x_xz_yy_xx = buffer_pddd[90];

    auto g_x_xz_yy_xy = buffer_pddd[91];

    auto g_x_xz_yy_xz = buffer_pddd[92];

    auto g_x_xz_yy_yy = buffer_pddd[93];

    auto g_x_xz_yy_yz = buffer_pddd[94];

    auto g_x_xz_yy_zz = buffer_pddd[95];

    auto g_x_xz_yz_xx = buffer_pddd[96];

    auto g_x_xz_yz_xy = buffer_pddd[97];

    auto g_x_xz_yz_xz = buffer_pddd[98];

    auto g_x_xz_yz_yy = buffer_pddd[99];

    auto g_x_xz_yz_yz = buffer_pddd[100];

    auto g_x_xz_yz_zz = buffer_pddd[101];

    auto g_x_xz_zz_xx = buffer_pddd[102];

    auto g_x_xz_zz_xy = buffer_pddd[103];

    auto g_x_xz_zz_xz = buffer_pddd[104];

    auto g_x_xz_zz_yy = buffer_pddd[105];

    auto g_x_xz_zz_yz = buffer_pddd[106];

    auto g_x_xz_zz_zz = buffer_pddd[107];

    auto g_x_yy_xx_xx = buffer_pddd[108];

    auto g_x_yy_xx_xy = buffer_pddd[109];

    auto g_x_yy_xx_xz = buffer_pddd[110];

    auto g_x_yy_xx_yy = buffer_pddd[111];

    auto g_x_yy_xx_yz = buffer_pddd[112];

    auto g_x_yy_xx_zz = buffer_pddd[113];

    auto g_x_yy_xy_xx = buffer_pddd[114];

    auto g_x_yy_xy_xy = buffer_pddd[115];

    auto g_x_yy_xy_xz = buffer_pddd[116];

    auto g_x_yy_xy_yy = buffer_pddd[117];

    auto g_x_yy_xy_yz = buffer_pddd[118];

    auto g_x_yy_xy_zz = buffer_pddd[119];

    auto g_x_yy_xz_xx = buffer_pddd[120];

    auto g_x_yy_xz_xy = buffer_pddd[121];

    auto g_x_yy_xz_xz = buffer_pddd[122];

    auto g_x_yy_xz_yy = buffer_pddd[123];

    auto g_x_yy_xz_yz = buffer_pddd[124];

    auto g_x_yy_xz_zz = buffer_pddd[125];

    auto g_x_yy_yy_xx = buffer_pddd[126];

    auto g_x_yy_yy_xy = buffer_pddd[127];

    auto g_x_yy_yy_xz = buffer_pddd[128];

    auto g_x_yy_yy_yy = buffer_pddd[129];

    auto g_x_yy_yy_yz = buffer_pddd[130];

    auto g_x_yy_yy_zz = buffer_pddd[131];

    auto g_x_yy_yz_xx = buffer_pddd[132];

    auto g_x_yy_yz_xy = buffer_pddd[133];

    auto g_x_yy_yz_xz = buffer_pddd[134];

    auto g_x_yy_yz_yy = buffer_pddd[135];

    auto g_x_yy_yz_yz = buffer_pddd[136];

    auto g_x_yy_yz_zz = buffer_pddd[137];

    auto g_x_yy_zz_xx = buffer_pddd[138];

    auto g_x_yy_zz_xy = buffer_pddd[139];

    auto g_x_yy_zz_xz = buffer_pddd[140];

    auto g_x_yy_zz_yy = buffer_pddd[141];

    auto g_x_yy_zz_yz = buffer_pddd[142];

    auto g_x_yy_zz_zz = buffer_pddd[143];

    auto g_x_yz_xx_xx = buffer_pddd[144];

    auto g_x_yz_xx_xy = buffer_pddd[145];

    auto g_x_yz_xx_xz = buffer_pddd[146];

    auto g_x_yz_xx_yy = buffer_pddd[147];

    auto g_x_yz_xx_yz = buffer_pddd[148];

    auto g_x_yz_xx_zz = buffer_pddd[149];

    auto g_x_yz_xy_xx = buffer_pddd[150];

    auto g_x_yz_xy_xy = buffer_pddd[151];

    auto g_x_yz_xy_xz = buffer_pddd[152];

    auto g_x_yz_xy_yy = buffer_pddd[153];

    auto g_x_yz_xy_yz = buffer_pddd[154];

    auto g_x_yz_xy_zz = buffer_pddd[155];

    auto g_x_yz_xz_xx = buffer_pddd[156];

    auto g_x_yz_xz_xy = buffer_pddd[157];

    auto g_x_yz_xz_xz = buffer_pddd[158];

    auto g_x_yz_xz_yy = buffer_pddd[159];

    auto g_x_yz_xz_yz = buffer_pddd[160];

    auto g_x_yz_xz_zz = buffer_pddd[161];

    auto g_x_yz_yy_xx = buffer_pddd[162];

    auto g_x_yz_yy_xy = buffer_pddd[163];

    auto g_x_yz_yy_xz = buffer_pddd[164];

    auto g_x_yz_yy_yy = buffer_pddd[165];

    auto g_x_yz_yy_yz = buffer_pddd[166];

    auto g_x_yz_yy_zz = buffer_pddd[167];

    auto g_x_yz_yz_xx = buffer_pddd[168];

    auto g_x_yz_yz_xy = buffer_pddd[169];

    auto g_x_yz_yz_xz = buffer_pddd[170];

    auto g_x_yz_yz_yy = buffer_pddd[171];

    auto g_x_yz_yz_yz = buffer_pddd[172];

    auto g_x_yz_yz_zz = buffer_pddd[173];

    auto g_x_yz_zz_xx = buffer_pddd[174];

    auto g_x_yz_zz_xy = buffer_pddd[175];

    auto g_x_yz_zz_xz = buffer_pddd[176];

    auto g_x_yz_zz_yy = buffer_pddd[177];

    auto g_x_yz_zz_yz = buffer_pddd[178];

    auto g_x_yz_zz_zz = buffer_pddd[179];

    auto g_x_zz_xx_xx = buffer_pddd[180];

    auto g_x_zz_xx_xy = buffer_pddd[181];

    auto g_x_zz_xx_xz = buffer_pddd[182];

    auto g_x_zz_xx_yy = buffer_pddd[183];

    auto g_x_zz_xx_yz = buffer_pddd[184];

    auto g_x_zz_xx_zz = buffer_pddd[185];

    auto g_x_zz_xy_xx = buffer_pddd[186];

    auto g_x_zz_xy_xy = buffer_pddd[187];

    auto g_x_zz_xy_xz = buffer_pddd[188];

    auto g_x_zz_xy_yy = buffer_pddd[189];

    auto g_x_zz_xy_yz = buffer_pddd[190];

    auto g_x_zz_xy_zz = buffer_pddd[191];

    auto g_x_zz_xz_xx = buffer_pddd[192];

    auto g_x_zz_xz_xy = buffer_pddd[193];

    auto g_x_zz_xz_xz = buffer_pddd[194];

    auto g_x_zz_xz_yy = buffer_pddd[195];

    auto g_x_zz_xz_yz = buffer_pddd[196];

    auto g_x_zz_xz_zz = buffer_pddd[197];

    auto g_x_zz_yy_xx = buffer_pddd[198];

    auto g_x_zz_yy_xy = buffer_pddd[199];

    auto g_x_zz_yy_xz = buffer_pddd[200];

    auto g_x_zz_yy_yy = buffer_pddd[201];

    auto g_x_zz_yy_yz = buffer_pddd[202];

    auto g_x_zz_yy_zz = buffer_pddd[203];

    auto g_x_zz_yz_xx = buffer_pddd[204];

    auto g_x_zz_yz_xy = buffer_pddd[205];

    auto g_x_zz_yz_xz = buffer_pddd[206];

    auto g_x_zz_yz_yy = buffer_pddd[207];

    auto g_x_zz_yz_yz = buffer_pddd[208];

    auto g_x_zz_yz_zz = buffer_pddd[209];

    auto g_x_zz_zz_xx = buffer_pddd[210];

    auto g_x_zz_zz_xy = buffer_pddd[211];

    auto g_x_zz_zz_xz = buffer_pddd[212];

    auto g_x_zz_zz_yy = buffer_pddd[213];

    auto g_x_zz_zz_yz = buffer_pddd[214];

    auto g_x_zz_zz_zz = buffer_pddd[215];

    auto g_y_xx_xx_xx = buffer_pddd[216];

    auto g_y_xx_xx_xy = buffer_pddd[217];

    auto g_y_xx_xx_xz = buffer_pddd[218];

    auto g_y_xx_xx_yy = buffer_pddd[219];

    auto g_y_xx_xx_yz = buffer_pddd[220];

    auto g_y_xx_xx_zz = buffer_pddd[221];

    auto g_y_xx_xy_xx = buffer_pddd[222];

    auto g_y_xx_xy_xy = buffer_pddd[223];

    auto g_y_xx_xy_xz = buffer_pddd[224];

    auto g_y_xx_xy_yy = buffer_pddd[225];

    auto g_y_xx_xy_yz = buffer_pddd[226];

    auto g_y_xx_xy_zz = buffer_pddd[227];

    auto g_y_xx_xz_xx = buffer_pddd[228];

    auto g_y_xx_xz_xy = buffer_pddd[229];

    auto g_y_xx_xz_xz = buffer_pddd[230];

    auto g_y_xx_xz_yy = buffer_pddd[231];

    auto g_y_xx_xz_yz = buffer_pddd[232];

    auto g_y_xx_xz_zz = buffer_pddd[233];

    auto g_y_xx_yy_xx = buffer_pddd[234];

    auto g_y_xx_yy_xy = buffer_pddd[235];

    auto g_y_xx_yy_xz = buffer_pddd[236];

    auto g_y_xx_yy_yy = buffer_pddd[237];

    auto g_y_xx_yy_yz = buffer_pddd[238];

    auto g_y_xx_yy_zz = buffer_pddd[239];

    auto g_y_xx_yz_xx = buffer_pddd[240];

    auto g_y_xx_yz_xy = buffer_pddd[241];

    auto g_y_xx_yz_xz = buffer_pddd[242];

    auto g_y_xx_yz_yy = buffer_pddd[243];

    auto g_y_xx_yz_yz = buffer_pddd[244];

    auto g_y_xx_yz_zz = buffer_pddd[245];

    auto g_y_xx_zz_xx = buffer_pddd[246];

    auto g_y_xx_zz_xy = buffer_pddd[247];

    auto g_y_xx_zz_xz = buffer_pddd[248];

    auto g_y_xx_zz_yy = buffer_pddd[249];

    auto g_y_xx_zz_yz = buffer_pddd[250];

    auto g_y_xx_zz_zz = buffer_pddd[251];

    auto g_y_xy_xx_xx = buffer_pddd[252];

    auto g_y_xy_xx_xy = buffer_pddd[253];

    auto g_y_xy_xx_xz = buffer_pddd[254];

    auto g_y_xy_xx_yy = buffer_pddd[255];

    auto g_y_xy_xx_yz = buffer_pddd[256];

    auto g_y_xy_xx_zz = buffer_pddd[257];

    auto g_y_xy_xy_xx = buffer_pddd[258];

    auto g_y_xy_xy_xy = buffer_pddd[259];

    auto g_y_xy_xy_xz = buffer_pddd[260];

    auto g_y_xy_xy_yy = buffer_pddd[261];

    auto g_y_xy_xy_yz = buffer_pddd[262];

    auto g_y_xy_xy_zz = buffer_pddd[263];

    auto g_y_xy_xz_xx = buffer_pddd[264];

    auto g_y_xy_xz_xy = buffer_pddd[265];

    auto g_y_xy_xz_xz = buffer_pddd[266];

    auto g_y_xy_xz_yy = buffer_pddd[267];

    auto g_y_xy_xz_yz = buffer_pddd[268];

    auto g_y_xy_xz_zz = buffer_pddd[269];

    auto g_y_xy_yy_xx = buffer_pddd[270];

    auto g_y_xy_yy_xy = buffer_pddd[271];

    auto g_y_xy_yy_xz = buffer_pddd[272];

    auto g_y_xy_yy_yy = buffer_pddd[273];

    auto g_y_xy_yy_yz = buffer_pddd[274];

    auto g_y_xy_yy_zz = buffer_pddd[275];

    auto g_y_xy_yz_xx = buffer_pddd[276];

    auto g_y_xy_yz_xy = buffer_pddd[277];

    auto g_y_xy_yz_xz = buffer_pddd[278];

    auto g_y_xy_yz_yy = buffer_pddd[279];

    auto g_y_xy_yz_yz = buffer_pddd[280];

    auto g_y_xy_yz_zz = buffer_pddd[281];

    auto g_y_xy_zz_xx = buffer_pddd[282];

    auto g_y_xy_zz_xy = buffer_pddd[283];

    auto g_y_xy_zz_xz = buffer_pddd[284];

    auto g_y_xy_zz_yy = buffer_pddd[285];

    auto g_y_xy_zz_yz = buffer_pddd[286];

    auto g_y_xy_zz_zz = buffer_pddd[287];

    auto g_y_xz_xx_xx = buffer_pddd[288];

    auto g_y_xz_xx_xy = buffer_pddd[289];

    auto g_y_xz_xx_xz = buffer_pddd[290];

    auto g_y_xz_xx_yy = buffer_pddd[291];

    auto g_y_xz_xx_yz = buffer_pddd[292];

    auto g_y_xz_xx_zz = buffer_pddd[293];

    auto g_y_xz_xy_xx = buffer_pddd[294];

    auto g_y_xz_xy_xy = buffer_pddd[295];

    auto g_y_xz_xy_xz = buffer_pddd[296];

    auto g_y_xz_xy_yy = buffer_pddd[297];

    auto g_y_xz_xy_yz = buffer_pddd[298];

    auto g_y_xz_xy_zz = buffer_pddd[299];

    auto g_y_xz_xz_xx = buffer_pddd[300];

    auto g_y_xz_xz_xy = buffer_pddd[301];

    auto g_y_xz_xz_xz = buffer_pddd[302];

    auto g_y_xz_xz_yy = buffer_pddd[303];

    auto g_y_xz_xz_yz = buffer_pddd[304];

    auto g_y_xz_xz_zz = buffer_pddd[305];

    auto g_y_xz_yy_xx = buffer_pddd[306];

    auto g_y_xz_yy_xy = buffer_pddd[307];

    auto g_y_xz_yy_xz = buffer_pddd[308];

    auto g_y_xz_yy_yy = buffer_pddd[309];

    auto g_y_xz_yy_yz = buffer_pddd[310];

    auto g_y_xz_yy_zz = buffer_pddd[311];

    auto g_y_xz_yz_xx = buffer_pddd[312];

    auto g_y_xz_yz_xy = buffer_pddd[313];

    auto g_y_xz_yz_xz = buffer_pddd[314];

    auto g_y_xz_yz_yy = buffer_pddd[315];

    auto g_y_xz_yz_yz = buffer_pddd[316];

    auto g_y_xz_yz_zz = buffer_pddd[317];

    auto g_y_xz_zz_xx = buffer_pddd[318];

    auto g_y_xz_zz_xy = buffer_pddd[319];

    auto g_y_xz_zz_xz = buffer_pddd[320];

    auto g_y_xz_zz_yy = buffer_pddd[321];

    auto g_y_xz_zz_yz = buffer_pddd[322];

    auto g_y_xz_zz_zz = buffer_pddd[323];

    auto g_y_yy_xx_xx = buffer_pddd[324];

    auto g_y_yy_xx_xy = buffer_pddd[325];

    auto g_y_yy_xx_xz = buffer_pddd[326];

    auto g_y_yy_xx_yy = buffer_pddd[327];

    auto g_y_yy_xx_yz = buffer_pddd[328];

    auto g_y_yy_xx_zz = buffer_pddd[329];

    auto g_y_yy_xy_xx = buffer_pddd[330];

    auto g_y_yy_xy_xy = buffer_pddd[331];

    auto g_y_yy_xy_xz = buffer_pddd[332];

    auto g_y_yy_xy_yy = buffer_pddd[333];

    auto g_y_yy_xy_yz = buffer_pddd[334];

    auto g_y_yy_xy_zz = buffer_pddd[335];

    auto g_y_yy_xz_xx = buffer_pddd[336];

    auto g_y_yy_xz_xy = buffer_pddd[337];

    auto g_y_yy_xz_xz = buffer_pddd[338];

    auto g_y_yy_xz_yy = buffer_pddd[339];

    auto g_y_yy_xz_yz = buffer_pddd[340];

    auto g_y_yy_xz_zz = buffer_pddd[341];

    auto g_y_yy_yy_xx = buffer_pddd[342];

    auto g_y_yy_yy_xy = buffer_pddd[343];

    auto g_y_yy_yy_xz = buffer_pddd[344];

    auto g_y_yy_yy_yy = buffer_pddd[345];

    auto g_y_yy_yy_yz = buffer_pddd[346];

    auto g_y_yy_yy_zz = buffer_pddd[347];

    auto g_y_yy_yz_xx = buffer_pddd[348];

    auto g_y_yy_yz_xy = buffer_pddd[349];

    auto g_y_yy_yz_xz = buffer_pddd[350];

    auto g_y_yy_yz_yy = buffer_pddd[351];

    auto g_y_yy_yz_yz = buffer_pddd[352];

    auto g_y_yy_yz_zz = buffer_pddd[353];

    auto g_y_yy_zz_xx = buffer_pddd[354];

    auto g_y_yy_zz_xy = buffer_pddd[355];

    auto g_y_yy_zz_xz = buffer_pddd[356];

    auto g_y_yy_zz_yy = buffer_pddd[357];

    auto g_y_yy_zz_yz = buffer_pddd[358];

    auto g_y_yy_zz_zz = buffer_pddd[359];

    auto g_y_yz_xx_xx = buffer_pddd[360];

    auto g_y_yz_xx_xy = buffer_pddd[361];

    auto g_y_yz_xx_xz = buffer_pddd[362];

    auto g_y_yz_xx_yy = buffer_pddd[363];

    auto g_y_yz_xx_yz = buffer_pddd[364];

    auto g_y_yz_xx_zz = buffer_pddd[365];

    auto g_y_yz_xy_xx = buffer_pddd[366];

    auto g_y_yz_xy_xy = buffer_pddd[367];

    auto g_y_yz_xy_xz = buffer_pddd[368];

    auto g_y_yz_xy_yy = buffer_pddd[369];

    auto g_y_yz_xy_yz = buffer_pddd[370];

    auto g_y_yz_xy_zz = buffer_pddd[371];

    auto g_y_yz_xz_xx = buffer_pddd[372];

    auto g_y_yz_xz_xy = buffer_pddd[373];

    auto g_y_yz_xz_xz = buffer_pddd[374];

    auto g_y_yz_xz_yy = buffer_pddd[375];

    auto g_y_yz_xz_yz = buffer_pddd[376];

    auto g_y_yz_xz_zz = buffer_pddd[377];

    auto g_y_yz_yy_xx = buffer_pddd[378];

    auto g_y_yz_yy_xy = buffer_pddd[379];

    auto g_y_yz_yy_xz = buffer_pddd[380];

    auto g_y_yz_yy_yy = buffer_pddd[381];

    auto g_y_yz_yy_yz = buffer_pddd[382];

    auto g_y_yz_yy_zz = buffer_pddd[383];

    auto g_y_yz_yz_xx = buffer_pddd[384];

    auto g_y_yz_yz_xy = buffer_pddd[385];

    auto g_y_yz_yz_xz = buffer_pddd[386];

    auto g_y_yz_yz_yy = buffer_pddd[387];

    auto g_y_yz_yz_yz = buffer_pddd[388];

    auto g_y_yz_yz_zz = buffer_pddd[389];

    auto g_y_yz_zz_xx = buffer_pddd[390];

    auto g_y_yz_zz_xy = buffer_pddd[391];

    auto g_y_yz_zz_xz = buffer_pddd[392];

    auto g_y_yz_zz_yy = buffer_pddd[393];

    auto g_y_yz_zz_yz = buffer_pddd[394];

    auto g_y_yz_zz_zz = buffer_pddd[395];

    auto g_y_zz_xx_xx = buffer_pddd[396];

    auto g_y_zz_xx_xy = buffer_pddd[397];

    auto g_y_zz_xx_xz = buffer_pddd[398];

    auto g_y_zz_xx_yy = buffer_pddd[399];

    auto g_y_zz_xx_yz = buffer_pddd[400];

    auto g_y_zz_xx_zz = buffer_pddd[401];

    auto g_y_zz_xy_xx = buffer_pddd[402];

    auto g_y_zz_xy_xy = buffer_pddd[403];

    auto g_y_zz_xy_xz = buffer_pddd[404];

    auto g_y_zz_xy_yy = buffer_pddd[405];

    auto g_y_zz_xy_yz = buffer_pddd[406];

    auto g_y_zz_xy_zz = buffer_pddd[407];

    auto g_y_zz_xz_xx = buffer_pddd[408];

    auto g_y_zz_xz_xy = buffer_pddd[409];

    auto g_y_zz_xz_xz = buffer_pddd[410];

    auto g_y_zz_xz_yy = buffer_pddd[411];

    auto g_y_zz_xz_yz = buffer_pddd[412];

    auto g_y_zz_xz_zz = buffer_pddd[413];

    auto g_y_zz_yy_xx = buffer_pddd[414];

    auto g_y_zz_yy_xy = buffer_pddd[415];

    auto g_y_zz_yy_xz = buffer_pddd[416];

    auto g_y_zz_yy_yy = buffer_pddd[417];

    auto g_y_zz_yy_yz = buffer_pddd[418];

    auto g_y_zz_yy_zz = buffer_pddd[419];

    auto g_y_zz_yz_xx = buffer_pddd[420];

    auto g_y_zz_yz_xy = buffer_pddd[421];

    auto g_y_zz_yz_xz = buffer_pddd[422];

    auto g_y_zz_yz_yy = buffer_pddd[423];

    auto g_y_zz_yz_yz = buffer_pddd[424];

    auto g_y_zz_yz_zz = buffer_pddd[425];

    auto g_y_zz_zz_xx = buffer_pddd[426];

    auto g_y_zz_zz_xy = buffer_pddd[427];

    auto g_y_zz_zz_xz = buffer_pddd[428];

    auto g_y_zz_zz_yy = buffer_pddd[429];

    auto g_y_zz_zz_yz = buffer_pddd[430];

    auto g_y_zz_zz_zz = buffer_pddd[431];

    auto g_z_xx_xx_xx = buffer_pddd[432];

    auto g_z_xx_xx_xy = buffer_pddd[433];

    auto g_z_xx_xx_xz = buffer_pddd[434];

    auto g_z_xx_xx_yy = buffer_pddd[435];

    auto g_z_xx_xx_yz = buffer_pddd[436];

    auto g_z_xx_xx_zz = buffer_pddd[437];

    auto g_z_xx_xy_xx = buffer_pddd[438];

    auto g_z_xx_xy_xy = buffer_pddd[439];

    auto g_z_xx_xy_xz = buffer_pddd[440];

    auto g_z_xx_xy_yy = buffer_pddd[441];

    auto g_z_xx_xy_yz = buffer_pddd[442];

    auto g_z_xx_xy_zz = buffer_pddd[443];

    auto g_z_xx_xz_xx = buffer_pddd[444];

    auto g_z_xx_xz_xy = buffer_pddd[445];

    auto g_z_xx_xz_xz = buffer_pddd[446];

    auto g_z_xx_xz_yy = buffer_pddd[447];

    auto g_z_xx_xz_yz = buffer_pddd[448];

    auto g_z_xx_xz_zz = buffer_pddd[449];

    auto g_z_xx_yy_xx = buffer_pddd[450];

    auto g_z_xx_yy_xy = buffer_pddd[451];

    auto g_z_xx_yy_xz = buffer_pddd[452];

    auto g_z_xx_yy_yy = buffer_pddd[453];

    auto g_z_xx_yy_yz = buffer_pddd[454];

    auto g_z_xx_yy_zz = buffer_pddd[455];

    auto g_z_xx_yz_xx = buffer_pddd[456];

    auto g_z_xx_yz_xy = buffer_pddd[457];

    auto g_z_xx_yz_xz = buffer_pddd[458];

    auto g_z_xx_yz_yy = buffer_pddd[459];

    auto g_z_xx_yz_yz = buffer_pddd[460];

    auto g_z_xx_yz_zz = buffer_pddd[461];

    auto g_z_xx_zz_xx = buffer_pddd[462];

    auto g_z_xx_zz_xy = buffer_pddd[463];

    auto g_z_xx_zz_xz = buffer_pddd[464];

    auto g_z_xx_zz_yy = buffer_pddd[465];

    auto g_z_xx_zz_yz = buffer_pddd[466];

    auto g_z_xx_zz_zz = buffer_pddd[467];

    auto g_z_xy_xx_xx = buffer_pddd[468];

    auto g_z_xy_xx_xy = buffer_pddd[469];

    auto g_z_xy_xx_xz = buffer_pddd[470];

    auto g_z_xy_xx_yy = buffer_pddd[471];

    auto g_z_xy_xx_yz = buffer_pddd[472];

    auto g_z_xy_xx_zz = buffer_pddd[473];

    auto g_z_xy_xy_xx = buffer_pddd[474];

    auto g_z_xy_xy_xy = buffer_pddd[475];

    auto g_z_xy_xy_xz = buffer_pddd[476];

    auto g_z_xy_xy_yy = buffer_pddd[477];

    auto g_z_xy_xy_yz = buffer_pddd[478];

    auto g_z_xy_xy_zz = buffer_pddd[479];

    auto g_z_xy_xz_xx = buffer_pddd[480];

    auto g_z_xy_xz_xy = buffer_pddd[481];

    auto g_z_xy_xz_xz = buffer_pddd[482];

    auto g_z_xy_xz_yy = buffer_pddd[483];

    auto g_z_xy_xz_yz = buffer_pddd[484];

    auto g_z_xy_xz_zz = buffer_pddd[485];

    auto g_z_xy_yy_xx = buffer_pddd[486];

    auto g_z_xy_yy_xy = buffer_pddd[487];

    auto g_z_xy_yy_xz = buffer_pddd[488];

    auto g_z_xy_yy_yy = buffer_pddd[489];

    auto g_z_xy_yy_yz = buffer_pddd[490];

    auto g_z_xy_yy_zz = buffer_pddd[491];

    auto g_z_xy_yz_xx = buffer_pddd[492];

    auto g_z_xy_yz_xy = buffer_pddd[493];

    auto g_z_xy_yz_xz = buffer_pddd[494];

    auto g_z_xy_yz_yy = buffer_pddd[495];

    auto g_z_xy_yz_yz = buffer_pddd[496];

    auto g_z_xy_yz_zz = buffer_pddd[497];

    auto g_z_xy_zz_xx = buffer_pddd[498];

    auto g_z_xy_zz_xy = buffer_pddd[499];

    auto g_z_xy_zz_xz = buffer_pddd[500];

    auto g_z_xy_zz_yy = buffer_pddd[501];

    auto g_z_xy_zz_yz = buffer_pddd[502];

    auto g_z_xy_zz_zz = buffer_pddd[503];

    auto g_z_xz_xx_xx = buffer_pddd[504];

    auto g_z_xz_xx_xy = buffer_pddd[505];

    auto g_z_xz_xx_xz = buffer_pddd[506];

    auto g_z_xz_xx_yy = buffer_pddd[507];

    auto g_z_xz_xx_yz = buffer_pddd[508];

    auto g_z_xz_xx_zz = buffer_pddd[509];

    auto g_z_xz_xy_xx = buffer_pddd[510];

    auto g_z_xz_xy_xy = buffer_pddd[511];

    auto g_z_xz_xy_xz = buffer_pddd[512];

    auto g_z_xz_xy_yy = buffer_pddd[513];

    auto g_z_xz_xy_yz = buffer_pddd[514];

    auto g_z_xz_xy_zz = buffer_pddd[515];

    auto g_z_xz_xz_xx = buffer_pddd[516];

    auto g_z_xz_xz_xy = buffer_pddd[517];

    auto g_z_xz_xz_xz = buffer_pddd[518];

    auto g_z_xz_xz_yy = buffer_pddd[519];

    auto g_z_xz_xz_yz = buffer_pddd[520];

    auto g_z_xz_xz_zz = buffer_pddd[521];

    auto g_z_xz_yy_xx = buffer_pddd[522];

    auto g_z_xz_yy_xy = buffer_pddd[523];

    auto g_z_xz_yy_xz = buffer_pddd[524];

    auto g_z_xz_yy_yy = buffer_pddd[525];

    auto g_z_xz_yy_yz = buffer_pddd[526];

    auto g_z_xz_yy_zz = buffer_pddd[527];

    auto g_z_xz_yz_xx = buffer_pddd[528];

    auto g_z_xz_yz_xy = buffer_pddd[529];

    auto g_z_xz_yz_xz = buffer_pddd[530];

    auto g_z_xz_yz_yy = buffer_pddd[531];

    auto g_z_xz_yz_yz = buffer_pddd[532];

    auto g_z_xz_yz_zz = buffer_pddd[533];

    auto g_z_xz_zz_xx = buffer_pddd[534];

    auto g_z_xz_zz_xy = buffer_pddd[535];

    auto g_z_xz_zz_xz = buffer_pddd[536];

    auto g_z_xz_zz_yy = buffer_pddd[537];

    auto g_z_xz_zz_yz = buffer_pddd[538];

    auto g_z_xz_zz_zz = buffer_pddd[539];

    auto g_z_yy_xx_xx = buffer_pddd[540];

    auto g_z_yy_xx_xy = buffer_pddd[541];

    auto g_z_yy_xx_xz = buffer_pddd[542];

    auto g_z_yy_xx_yy = buffer_pddd[543];

    auto g_z_yy_xx_yz = buffer_pddd[544];

    auto g_z_yy_xx_zz = buffer_pddd[545];

    auto g_z_yy_xy_xx = buffer_pddd[546];

    auto g_z_yy_xy_xy = buffer_pddd[547];

    auto g_z_yy_xy_xz = buffer_pddd[548];

    auto g_z_yy_xy_yy = buffer_pddd[549];

    auto g_z_yy_xy_yz = buffer_pddd[550];

    auto g_z_yy_xy_zz = buffer_pddd[551];

    auto g_z_yy_xz_xx = buffer_pddd[552];

    auto g_z_yy_xz_xy = buffer_pddd[553];

    auto g_z_yy_xz_xz = buffer_pddd[554];

    auto g_z_yy_xz_yy = buffer_pddd[555];

    auto g_z_yy_xz_yz = buffer_pddd[556];

    auto g_z_yy_xz_zz = buffer_pddd[557];

    auto g_z_yy_yy_xx = buffer_pddd[558];

    auto g_z_yy_yy_xy = buffer_pddd[559];

    auto g_z_yy_yy_xz = buffer_pddd[560];

    auto g_z_yy_yy_yy = buffer_pddd[561];

    auto g_z_yy_yy_yz = buffer_pddd[562];

    auto g_z_yy_yy_zz = buffer_pddd[563];

    auto g_z_yy_yz_xx = buffer_pddd[564];

    auto g_z_yy_yz_xy = buffer_pddd[565];

    auto g_z_yy_yz_xz = buffer_pddd[566];

    auto g_z_yy_yz_yy = buffer_pddd[567];

    auto g_z_yy_yz_yz = buffer_pddd[568];

    auto g_z_yy_yz_zz = buffer_pddd[569];

    auto g_z_yy_zz_xx = buffer_pddd[570];

    auto g_z_yy_zz_xy = buffer_pddd[571];

    auto g_z_yy_zz_xz = buffer_pddd[572];

    auto g_z_yy_zz_yy = buffer_pddd[573];

    auto g_z_yy_zz_yz = buffer_pddd[574];

    auto g_z_yy_zz_zz = buffer_pddd[575];

    auto g_z_yz_xx_xx = buffer_pddd[576];

    auto g_z_yz_xx_xy = buffer_pddd[577];

    auto g_z_yz_xx_xz = buffer_pddd[578];

    auto g_z_yz_xx_yy = buffer_pddd[579];

    auto g_z_yz_xx_yz = buffer_pddd[580];

    auto g_z_yz_xx_zz = buffer_pddd[581];

    auto g_z_yz_xy_xx = buffer_pddd[582];

    auto g_z_yz_xy_xy = buffer_pddd[583];

    auto g_z_yz_xy_xz = buffer_pddd[584];

    auto g_z_yz_xy_yy = buffer_pddd[585];

    auto g_z_yz_xy_yz = buffer_pddd[586];

    auto g_z_yz_xy_zz = buffer_pddd[587];

    auto g_z_yz_xz_xx = buffer_pddd[588];

    auto g_z_yz_xz_xy = buffer_pddd[589];

    auto g_z_yz_xz_xz = buffer_pddd[590];

    auto g_z_yz_xz_yy = buffer_pddd[591];

    auto g_z_yz_xz_yz = buffer_pddd[592];

    auto g_z_yz_xz_zz = buffer_pddd[593];

    auto g_z_yz_yy_xx = buffer_pddd[594];

    auto g_z_yz_yy_xy = buffer_pddd[595];

    auto g_z_yz_yy_xz = buffer_pddd[596];

    auto g_z_yz_yy_yy = buffer_pddd[597];

    auto g_z_yz_yy_yz = buffer_pddd[598];

    auto g_z_yz_yy_zz = buffer_pddd[599];

    auto g_z_yz_yz_xx = buffer_pddd[600];

    auto g_z_yz_yz_xy = buffer_pddd[601];

    auto g_z_yz_yz_xz = buffer_pddd[602];

    auto g_z_yz_yz_yy = buffer_pddd[603];

    auto g_z_yz_yz_yz = buffer_pddd[604];

    auto g_z_yz_yz_zz = buffer_pddd[605];

    auto g_z_yz_zz_xx = buffer_pddd[606];

    auto g_z_yz_zz_xy = buffer_pddd[607];

    auto g_z_yz_zz_xz = buffer_pddd[608];

    auto g_z_yz_zz_yy = buffer_pddd[609];

    auto g_z_yz_zz_yz = buffer_pddd[610];

    auto g_z_yz_zz_zz = buffer_pddd[611];

    auto g_z_zz_xx_xx = buffer_pddd[612];

    auto g_z_zz_xx_xy = buffer_pddd[613];

    auto g_z_zz_xx_xz = buffer_pddd[614];

    auto g_z_zz_xx_yy = buffer_pddd[615];

    auto g_z_zz_xx_yz = buffer_pddd[616];

    auto g_z_zz_xx_zz = buffer_pddd[617];

    auto g_z_zz_xy_xx = buffer_pddd[618];

    auto g_z_zz_xy_xy = buffer_pddd[619];

    auto g_z_zz_xy_xz = buffer_pddd[620];

    auto g_z_zz_xy_yy = buffer_pddd[621];

    auto g_z_zz_xy_yz = buffer_pddd[622];

    auto g_z_zz_xy_zz = buffer_pddd[623];

    auto g_z_zz_xz_xx = buffer_pddd[624];

    auto g_z_zz_xz_xy = buffer_pddd[625];

    auto g_z_zz_xz_xz = buffer_pddd[626];

    auto g_z_zz_xz_yy = buffer_pddd[627];

    auto g_z_zz_xz_yz = buffer_pddd[628];

    auto g_z_zz_xz_zz = buffer_pddd[629];

    auto g_z_zz_yy_xx = buffer_pddd[630];

    auto g_z_zz_yy_xy = buffer_pddd[631];

    auto g_z_zz_yy_xz = buffer_pddd[632];

    auto g_z_zz_yy_yy = buffer_pddd[633];

    auto g_z_zz_yy_yz = buffer_pddd[634];

    auto g_z_zz_yy_zz = buffer_pddd[635];

    auto g_z_zz_yz_xx = buffer_pddd[636];

    auto g_z_zz_yz_xy = buffer_pddd[637];

    auto g_z_zz_yz_xz = buffer_pddd[638];

    auto g_z_zz_yz_yy = buffer_pddd[639];

    auto g_z_zz_yz_yz = buffer_pddd[640];

    auto g_z_zz_yz_zz = buffer_pddd[641];

    auto g_z_zz_zz_xx = buffer_pddd[642];

    auto g_z_zz_zz_xy = buffer_pddd[643];

    auto g_z_zz_zz_xz = buffer_pddd[644];

    auto g_z_zz_zz_yy = buffer_pddd[645];

    auto g_z_zz_zz_yz = buffer_pddd[646];

    auto g_z_zz_zz_zz = buffer_pddd[647];

    /// Set up components of integrals buffer : buffer_1010_sdpd

    auto g_x_0_x_0_0_xx_x_xx = buffer_1010_sdpd[0];

    auto g_x_0_x_0_0_xx_x_xy = buffer_1010_sdpd[1];

    auto g_x_0_x_0_0_xx_x_xz = buffer_1010_sdpd[2];

    auto g_x_0_x_0_0_xx_x_yy = buffer_1010_sdpd[3];

    auto g_x_0_x_0_0_xx_x_yz = buffer_1010_sdpd[4];

    auto g_x_0_x_0_0_xx_x_zz = buffer_1010_sdpd[5];

    auto g_x_0_x_0_0_xx_y_xx = buffer_1010_sdpd[6];

    auto g_x_0_x_0_0_xx_y_xy = buffer_1010_sdpd[7];

    auto g_x_0_x_0_0_xx_y_xz = buffer_1010_sdpd[8];

    auto g_x_0_x_0_0_xx_y_yy = buffer_1010_sdpd[9];

    auto g_x_0_x_0_0_xx_y_yz = buffer_1010_sdpd[10];

    auto g_x_0_x_0_0_xx_y_zz = buffer_1010_sdpd[11];

    auto g_x_0_x_0_0_xx_z_xx = buffer_1010_sdpd[12];

    auto g_x_0_x_0_0_xx_z_xy = buffer_1010_sdpd[13];

    auto g_x_0_x_0_0_xx_z_xz = buffer_1010_sdpd[14];

    auto g_x_0_x_0_0_xx_z_yy = buffer_1010_sdpd[15];

    auto g_x_0_x_0_0_xx_z_yz = buffer_1010_sdpd[16];

    auto g_x_0_x_0_0_xx_z_zz = buffer_1010_sdpd[17];

    auto g_x_0_x_0_0_xy_x_xx = buffer_1010_sdpd[18];

    auto g_x_0_x_0_0_xy_x_xy = buffer_1010_sdpd[19];

    auto g_x_0_x_0_0_xy_x_xz = buffer_1010_sdpd[20];

    auto g_x_0_x_0_0_xy_x_yy = buffer_1010_sdpd[21];

    auto g_x_0_x_0_0_xy_x_yz = buffer_1010_sdpd[22];

    auto g_x_0_x_0_0_xy_x_zz = buffer_1010_sdpd[23];

    auto g_x_0_x_0_0_xy_y_xx = buffer_1010_sdpd[24];

    auto g_x_0_x_0_0_xy_y_xy = buffer_1010_sdpd[25];

    auto g_x_0_x_0_0_xy_y_xz = buffer_1010_sdpd[26];

    auto g_x_0_x_0_0_xy_y_yy = buffer_1010_sdpd[27];

    auto g_x_0_x_0_0_xy_y_yz = buffer_1010_sdpd[28];

    auto g_x_0_x_0_0_xy_y_zz = buffer_1010_sdpd[29];

    auto g_x_0_x_0_0_xy_z_xx = buffer_1010_sdpd[30];

    auto g_x_0_x_0_0_xy_z_xy = buffer_1010_sdpd[31];

    auto g_x_0_x_0_0_xy_z_xz = buffer_1010_sdpd[32];

    auto g_x_0_x_0_0_xy_z_yy = buffer_1010_sdpd[33];

    auto g_x_0_x_0_0_xy_z_yz = buffer_1010_sdpd[34];

    auto g_x_0_x_0_0_xy_z_zz = buffer_1010_sdpd[35];

    auto g_x_0_x_0_0_xz_x_xx = buffer_1010_sdpd[36];

    auto g_x_0_x_0_0_xz_x_xy = buffer_1010_sdpd[37];

    auto g_x_0_x_0_0_xz_x_xz = buffer_1010_sdpd[38];

    auto g_x_0_x_0_0_xz_x_yy = buffer_1010_sdpd[39];

    auto g_x_0_x_0_0_xz_x_yz = buffer_1010_sdpd[40];

    auto g_x_0_x_0_0_xz_x_zz = buffer_1010_sdpd[41];

    auto g_x_0_x_0_0_xz_y_xx = buffer_1010_sdpd[42];

    auto g_x_0_x_0_0_xz_y_xy = buffer_1010_sdpd[43];

    auto g_x_0_x_0_0_xz_y_xz = buffer_1010_sdpd[44];

    auto g_x_0_x_0_0_xz_y_yy = buffer_1010_sdpd[45];

    auto g_x_0_x_0_0_xz_y_yz = buffer_1010_sdpd[46];

    auto g_x_0_x_0_0_xz_y_zz = buffer_1010_sdpd[47];

    auto g_x_0_x_0_0_xz_z_xx = buffer_1010_sdpd[48];

    auto g_x_0_x_0_0_xz_z_xy = buffer_1010_sdpd[49];

    auto g_x_0_x_0_0_xz_z_xz = buffer_1010_sdpd[50];

    auto g_x_0_x_0_0_xz_z_yy = buffer_1010_sdpd[51];

    auto g_x_0_x_0_0_xz_z_yz = buffer_1010_sdpd[52];

    auto g_x_0_x_0_0_xz_z_zz = buffer_1010_sdpd[53];

    auto g_x_0_x_0_0_yy_x_xx = buffer_1010_sdpd[54];

    auto g_x_0_x_0_0_yy_x_xy = buffer_1010_sdpd[55];

    auto g_x_0_x_0_0_yy_x_xz = buffer_1010_sdpd[56];

    auto g_x_0_x_0_0_yy_x_yy = buffer_1010_sdpd[57];

    auto g_x_0_x_0_0_yy_x_yz = buffer_1010_sdpd[58];

    auto g_x_0_x_0_0_yy_x_zz = buffer_1010_sdpd[59];

    auto g_x_0_x_0_0_yy_y_xx = buffer_1010_sdpd[60];

    auto g_x_0_x_0_0_yy_y_xy = buffer_1010_sdpd[61];

    auto g_x_0_x_0_0_yy_y_xz = buffer_1010_sdpd[62];

    auto g_x_0_x_0_0_yy_y_yy = buffer_1010_sdpd[63];

    auto g_x_0_x_0_0_yy_y_yz = buffer_1010_sdpd[64];

    auto g_x_0_x_0_0_yy_y_zz = buffer_1010_sdpd[65];

    auto g_x_0_x_0_0_yy_z_xx = buffer_1010_sdpd[66];

    auto g_x_0_x_0_0_yy_z_xy = buffer_1010_sdpd[67];

    auto g_x_0_x_0_0_yy_z_xz = buffer_1010_sdpd[68];

    auto g_x_0_x_0_0_yy_z_yy = buffer_1010_sdpd[69];

    auto g_x_0_x_0_0_yy_z_yz = buffer_1010_sdpd[70];

    auto g_x_0_x_0_0_yy_z_zz = buffer_1010_sdpd[71];

    auto g_x_0_x_0_0_yz_x_xx = buffer_1010_sdpd[72];

    auto g_x_0_x_0_0_yz_x_xy = buffer_1010_sdpd[73];

    auto g_x_0_x_0_0_yz_x_xz = buffer_1010_sdpd[74];

    auto g_x_0_x_0_0_yz_x_yy = buffer_1010_sdpd[75];

    auto g_x_0_x_0_0_yz_x_yz = buffer_1010_sdpd[76];

    auto g_x_0_x_0_0_yz_x_zz = buffer_1010_sdpd[77];

    auto g_x_0_x_0_0_yz_y_xx = buffer_1010_sdpd[78];

    auto g_x_0_x_0_0_yz_y_xy = buffer_1010_sdpd[79];

    auto g_x_0_x_0_0_yz_y_xz = buffer_1010_sdpd[80];

    auto g_x_0_x_0_0_yz_y_yy = buffer_1010_sdpd[81];

    auto g_x_0_x_0_0_yz_y_yz = buffer_1010_sdpd[82];

    auto g_x_0_x_0_0_yz_y_zz = buffer_1010_sdpd[83];

    auto g_x_0_x_0_0_yz_z_xx = buffer_1010_sdpd[84];

    auto g_x_0_x_0_0_yz_z_xy = buffer_1010_sdpd[85];

    auto g_x_0_x_0_0_yz_z_xz = buffer_1010_sdpd[86];

    auto g_x_0_x_0_0_yz_z_yy = buffer_1010_sdpd[87];

    auto g_x_0_x_0_0_yz_z_yz = buffer_1010_sdpd[88];

    auto g_x_0_x_0_0_yz_z_zz = buffer_1010_sdpd[89];

    auto g_x_0_x_0_0_zz_x_xx = buffer_1010_sdpd[90];

    auto g_x_0_x_0_0_zz_x_xy = buffer_1010_sdpd[91];

    auto g_x_0_x_0_0_zz_x_xz = buffer_1010_sdpd[92];

    auto g_x_0_x_0_0_zz_x_yy = buffer_1010_sdpd[93];

    auto g_x_0_x_0_0_zz_x_yz = buffer_1010_sdpd[94];

    auto g_x_0_x_0_0_zz_x_zz = buffer_1010_sdpd[95];

    auto g_x_0_x_0_0_zz_y_xx = buffer_1010_sdpd[96];

    auto g_x_0_x_0_0_zz_y_xy = buffer_1010_sdpd[97];

    auto g_x_0_x_0_0_zz_y_xz = buffer_1010_sdpd[98];

    auto g_x_0_x_0_0_zz_y_yy = buffer_1010_sdpd[99];

    auto g_x_0_x_0_0_zz_y_yz = buffer_1010_sdpd[100];

    auto g_x_0_x_0_0_zz_y_zz = buffer_1010_sdpd[101];

    auto g_x_0_x_0_0_zz_z_xx = buffer_1010_sdpd[102];

    auto g_x_0_x_0_0_zz_z_xy = buffer_1010_sdpd[103];

    auto g_x_0_x_0_0_zz_z_xz = buffer_1010_sdpd[104];

    auto g_x_0_x_0_0_zz_z_yy = buffer_1010_sdpd[105];

    auto g_x_0_x_0_0_zz_z_yz = buffer_1010_sdpd[106];

    auto g_x_0_x_0_0_zz_z_zz = buffer_1010_sdpd[107];

    auto g_x_0_y_0_0_xx_x_xx = buffer_1010_sdpd[108];

    auto g_x_0_y_0_0_xx_x_xy = buffer_1010_sdpd[109];

    auto g_x_0_y_0_0_xx_x_xz = buffer_1010_sdpd[110];

    auto g_x_0_y_0_0_xx_x_yy = buffer_1010_sdpd[111];

    auto g_x_0_y_0_0_xx_x_yz = buffer_1010_sdpd[112];

    auto g_x_0_y_0_0_xx_x_zz = buffer_1010_sdpd[113];

    auto g_x_0_y_0_0_xx_y_xx = buffer_1010_sdpd[114];

    auto g_x_0_y_0_0_xx_y_xy = buffer_1010_sdpd[115];

    auto g_x_0_y_0_0_xx_y_xz = buffer_1010_sdpd[116];

    auto g_x_0_y_0_0_xx_y_yy = buffer_1010_sdpd[117];

    auto g_x_0_y_0_0_xx_y_yz = buffer_1010_sdpd[118];

    auto g_x_0_y_0_0_xx_y_zz = buffer_1010_sdpd[119];

    auto g_x_0_y_0_0_xx_z_xx = buffer_1010_sdpd[120];

    auto g_x_0_y_0_0_xx_z_xy = buffer_1010_sdpd[121];

    auto g_x_0_y_0_0_xx_z_xz = buffer_1010_sdpd[122];

    auto g_x_0_y_0_0_xx_z_yy = buffer_1010_sdpd[123];

    auto g_x_0_y_0_0_xx_z_yz = buffer_1010_sdpd[124];

    auto g_x_0_y_0_0_xx_z_zz = buffer_1010_sdpd[125];

    auto g_x_0_y_0_0_xy_x_xx = buffer_1010_sdpd[126];

    auto g_x_0_y_0_0_xy_x_xy = buffer_1010_sdpd[127];

    auto g_x_0_y_0_0_xy_x_xz = buffer_1010_sdpd[128];

    auto g_x_0_y_0_0_xy_x_yy = buffer_1010_sdpd[129];

    auto g_x_0_y_0_0_xy_x_yz = buffer_1010_sdpd[130];

    auto g_x_0_y_0_0_xy_x_zz = buffer_1010_sdpd[131];

    auto g_x_0_y_0_0_xy_y_xx = buffer_1010_sdpd[132];

    auto g_x_0_y_0_0_xy_y_xy = buffer_1010_sdpd[133];

    auto g_x_0_y_0_0_xy_y_xz = buffer_1010_sdpd[134];

    auto g_x_0_y_0_0_xy_y_yy = buffer_1010_sdpd[135];

    auto g_x_0_y_0_0_xy_y_yz = buffer_1010_sdpd[136];

    auto g_x_0_y_0_0_xy_y_zz = buffer_1010_sdpd[137];

    auto g_x_0_y_0_0_xy_z_xx = buffer_1010_sdpd[138];

    auto g_x_0_y_0_0_xy_z_xy = buffer_1010_sdpd[139];

    auto g_x_0_y_0_0_xy_z_xz = buffer_1010_sdpd[140];

    auto g_x_0_y_0_0_xy_z_yy = buffer_1010_sdpd[141];

    auto g_x_0_y_0_0_xy_z_yz = buffer_1010_sdpd[142];

    auto g_x_0_y_0_0_xy_z_zz = buffer_1010_sdpd[143];

    auto g_x_0_y_0_0_xz_x_xx = buffer_1010_sdpd[144];

    auto g_x_0_y_0_0_xz_x_xy = buffer_1010_sdpd[145];

    auto g_x_0_y_0_0_xz_x_xz = buffer_1010_sdpd[146];

    auto g_x_0_y_0_0_xz_x_yy = buffer_1010_sdpd[147];

    auto g_x_0_y_0_0_xz_x_yz = buffer_1010_sdpd[148];

    auto g_x_0_y_0_0_xz_x_zz = buffer_1010_sdpd[149];

    auto g_x_0_y_0_0_xz_y_xx = buffer_1010_sdpd[150];

    auto g_x_0_y_0_0_xz_y_xy = buffer_1010_sdpd[151];

    auto g_x_0_y_0_0_xz_y_xz = buffer_1010_sdpd[152];

    auto g_x_0_y_0_0_xz_y_yy = buffer_1010_sdpd[153];

    auto g_x_0_y_0_0_xz_y_yz = buffer_1010_sdpd[154];

    auto g_x_0_y_0_0_xz_y_zz = buffer_1010_sdpd[155];

    auto g_x_0_y_0_0_xz_z_xx = buffer_1010_sdpd[156];

    auto g_x_0_y_0_0_xz_z_xy = buffer_1010_sdpd[157];

    auto g_x_0_y_0_0_xz_z_xz = buffer_1010_sdpd[158];

    auto g_x_0_y_0_0_xz_z_yy = buffer_1010_sdpd[159];

    auto g_x_0_y_0_0_xz_z_yz = buffer_1010_sdpd[160];

    auto g_x_0_y_0_0_xz_z_zz = buffer_1010_sdpd[161];

    auto g_x_0_y_0_0_yy_x_xx = buffer_1010_sdpd[162];

    auto g_x_0_y_0_0_yy_x_xy = buffer_1010_sdpd[163];

    auto g_x_0_y_0_0_yy_x_xz = buffer_1010_sdpd[164];

    auto g_x_0_y_0_0_yy_x_yy = buffer_1010_sdpd[165];

    auto g_x_0_y_0_0_yy_x_yz = buffer_1010_sdpd[166];

    auto g_x_0_y_0_0_yy_x_zz = buffer_1010_sdpd[167];

    auto g_x_0_y_0_0_yy_y_xx = buffer_1010_sdpd[168];

    auto g_x_0_y_0_0_yy_y_xy = buffer_1010_sdpd[169];

    auto g_x_0_y_0_0_yy_y_xz = buffer_1010_sdpd[170];

    auto g_x_0_y_0_0_yy_y_yy = buffer_1010_sdpd[171];

    auto g_x_0_y_0_0_yy_y_yz = buffer_1010_sdpd[172];

    auto g_x_0_y_0_0_yy_y_zz = buffer_1010_sdpd[173];

    auto g_x_0_y_0_0_yy_z_xx = buffer_1010_sdpd[174];

    auto g_x_0_y_0_0_yy_z_xy = buffer_1010_sdpd[175];

    auto g_x_0_y_0_0_yy_z_xz = buffer_1010_sdpd[176];

    auto g_x_0_y_0_0_yy_z_yy = buffer_1010_sdpd[177];

    auto g_x_0_y_0_0_yy_z_yz = buffer_1010_sdpd[178];

    auto g_x_0_y_0_0_yy_z_zz = buffer_1010_sdpd[179];

    auto g_x_0_y_0_0_yz_x_xx = buffer_1010_sdpd[180];

    auto g_x_0_y_0_0_yz_x_xy = buffer_1010_sdpd[181];

    auto g_x_0_y_0_0_yz_x_xz = buffer_1010_sdpd[182];

    auto g_x_0_y_0_0_yz_x_yy = buffer_1010_sdpd[183];

    auto g_x_0_y_0_0_yz_x_yz = buffer_1010_sdpd[184];

    auto g_x_0_y_0_0_yz_x_zz = buffer_1010_sdpd[185];

    auto g_x_0_y_0_0_yz_y_xx = buffer_1010_sdpd[186];

    auto g_x_0_y_0_0_yz_y_xy = buffer_1010_sdpd[187];

    auto g_x_0_y_0_0_yz_y_xz = buffer_1010_sdpd[188];

    auto g_x_0_y_0_0_yz_y_yy = buffer_1010_sdpd[189];

    auto g_x_0_y_0_0_yz_y_yz = buffer_1010_sdpd[190];

    auto g_x_0_y_0_0_yz_y_zz = buffer_1010_sdpd[191];

    auto g_x_0_y_0_0_yz_z_xx = buffer_1010_sdpd[192];

    auto g_x_0_y_0_0_yz_z_xy = buffer_1010_sdpd[193];

    auto g_x_0_y_0_0_yz_z_xz = buffer_1010_sdpd[194];

    auto g_x_0_y_0_0_yz_z_yy = buffer_1010_sdpd[195];

    auto g_x_0_y_0_0_yz_z_yz = buffer_1010_sdpd[196];

    auto g_x_0_y_0_0_yz_z_zz = buffer_1010_sdpd[197];

    auto g_x_0_y_0_0_zz_x_xx = buffer_1010_sdpd[198];

    auto g_x_0_y_0_0_zz_x_xy = buffer_1010_sdpd[199];

    auto g_x_0_y_0_0_zz_x_xz = buffer_1010_sdpd[200];

    auto g_x_0_y_0_0_zz_x_yy = buffer_1010_sdpd[201];

    auto g_x_0_y_0_0_zz_x_yz = buffer_1010_sdpd[202];

    auto g_x_0_y_0_0_zz_x_zz = buffer_1010_sdpd[203];

    auto g_x_0_y_0_0_zz_y_xx = buffer_1010_sdpd[204];

    auto g_x_0_y_0_0_zz_y_xy = buffer_1010_sdpd[205];

    auto g_x_0_y_0_0_zz_y_xz = buffer_1010_sdpd[206];

    auto g_x_0_y_0_0_zz_y_yy = buffer_1010_sdpd[207];

    auto g_x_0_y_0_0_zz_y_yz = buffer_1010_sdpd[208];

    auto g_x_0_y_0_0_zz_y_zz = buffer_1010_sdpd[209];

    auto g_x_0_y_0_0_zz_z_xx = buffer_1010_sdpd[210];

    auto g_x_0_y_0_0_zz_z_xy = buffer_1010_sdpd[211];

    auto g_x_0_y_0_0_zz_z_xz = buffer_1010_sdpd[212];

    auto g_x_0_y_0_0_zz_z_yy = buffer_1010_sdpd[213];

    auto g_x_0_y_0_0_zz_z_yz = buffer_1010_sdpd[214];

    auto g_x_0_y_0_0_zz_z_zz = buffer_1010_sdpd[215];

    auto g_x_0_z_0_0_xx_x_xx = buffer_1010_sdpd[216];

    auto g_x_0_z_0_0_xx_x_xy = buffer_1010_sdpd[217];

    auto g_x_0_z_0_0_xx_x_xz = buffer_1010_sdpd[218];

    auto g_x_0_z_0_0_xx_x_yy = buffer_1010_sdpd[219];

    auto g_x_0_z_0_0_xx_x_yz = buffer_1010_sdpd[220];

    auto g_x_0_z_0_0_xx_x_zz = buffer_1010_sdpd[221];

    auto g_x_0_z_0_0_xx_y_xx = buffer_1010_sdpd[222];

    auto g_x_0_z_0_0_xx_y_xy = buffer_1010_sdpd[223];

    auto g_x_0_z_0_0_xx_y_xz = buffer_1010_sdpd[224];

    auto g_x_0_z_0_0_xx_y_yy = buffer_1010_sdpd[225];

    auto g_x_0_z_0_0_xx_y_yz = buffer_1010_sdpd[226];

    auto g_x_0_z_0_0_xx_y_zz = buffer_1010_sdpd[227];

    auto g_x_0_z_0_0_xx_z_xx = buffer_1010_sdpd[228];

    auto g_x_0_z_0_0_xx_z_xy = buffer_1010_sdpd[229];

    auto g_x_0_z_0_0_xx_z_xz = buffer_1010_sdpd[230];

    auto g_x_0_z_0_0_xx_z_yy = buffer_1010_sdpd[231];

    auto g_x_0_z_0_0_xx_z_yz = buffer_1010_sdpd[232];

    auto g_x_0_z_0_0_xx_z_zz = buffer_1010_sdpd[233];

    auto g_x_0_z_0_0_xy_x_xx = buffer_1010_sdpd[234];

    auto g_x_0_z_0_0_xy_x_xy = buffer_1010_sdpd[235];

    auto g_x_0_z_0_0_xy_x_xz = buffer_1010_sdpd[236];

    auto g_x_0_z_0_0_xy_x_yy = buffer_1010_sdpd[237];

    auto g_x_0_z_0_0_xy_x_yz = buffer_1010_sdpd[238];

    auto g_x_0_z_0_0_xy_x_zz = buffer_1010_sdpd[239];

    auto g_x_0_z_0_0_xy_y_xx = buffer_1010_sdpd[240];

    auto g_x_0_z_0_0_xy_y_xy = buffer_1010_sdpd[241];

    auto g_x_0_z_0_0_xy_y_xz = buffer_1010_sdpd[242];

    auto g_x_0_z_0_0_xy_y_yy = buffer_1010_sdpd[243];

    auto g_x_0_z_0_0_xy_y_yz = buffer_1010_sdpd[244];

    auto g_x_0_z_0_0_xy_y_zz = buffer_1010_sdpd[245];

    auto g_x_0_z_0_0_xy_z_xx = buffer_1010_sdpd[246];

    auto g_x_0_z_0_0_xy_z_xy = buffer_1010_sdpd[247];

    auto g_x_0_z_0_0_xy_z_xz = buffer_1010_sdpd[248];

    auto g_x_0_z_0_0_xy_z_yy = buffer_1010_sdpd[249];

    auto g_x_0_z_0_0_xy_z_yz = buffer_1010_sdpd[250];

    auto g_x_0_z_0_0_xy_z_zz = buffer_1010_sdpd[251];

    auto g_x_0_z_0_0_xz_x_xx = buffer_1010_sdpd[252];

    auto g_x_0_z_0_0_xz_x_xy = buffer_1010_sdpd[253];

    auto g_x_0_z_0_0_xz_x_xz = buffer_1010_sdpd[254];

    auto g_x_0_z_0_0_xz_x_yy = buffer_1010_sdpd[255];

    auto g_x_0_z_0_0_xz_x_yz = buffer_1010_sdpd[256];

    auto g_x_0_z_0_0_xz_x_zz = buffer_1010_sdpd[257];

    auto g_x_0_z_0_0_xz_y_xx = buffer_1010_sdpd[258];

    auto g_x_0_z_0_0_xz_y_xy = buffer_1010_sdpd[259];

    auto g_x_0_z_0_0_xz_y_xz = buffer_1010_sdpd[260];

    auto g_x_0_z_0_0_xz_y_yy = buffer_1010_sdpd[261];

    auto g_x_0_z_0_0_xz_y_yz = buffer_1010_sdpd[262];

    auto g_x_0_z_0_0_xz_y_zz = buffer_1010_sdpd[263];

    auto g_x_0_z_0_0_xz_z_xx = buffer_1010_sdpd[264];

    auto g_x_0_z_0_0_xz_z_xy = buffer_1010_sdpd[265];

    auto g_x_0_z_0_0_xz_z_xz = buffer_1010_sdpd[266];

    auto g_x_0_z_0_0_xz_z_yy = buffer_1010_sdpd[267];

    auto g_x_0_z_0_0_xz_z_yz = buffer_1010_sdpd[268];

    auto g_x_0_z_0_0_xz_z_zz = buffer_1010_sdpd[269];

    auto g_x_0_z_0_0_yy_x_xx = buffer_1010_sdpd[270];

    auto g_x_0_z_0_0_yy_x_xy = buffer_1010_sdpd[271];

    auto g_x_0_z_0_0_yy_x_xz = buffer_1010_sdpd[272];

    auto g_x_0_z_0_0_yy_x_yy = buffer_1010_sdpd[273];

    auto g_x_0_z_0_0_yy_x_yz = buffer_1010_sdpd[274];

    auto g_x_0_z_0_0_yy_x_zz = buffer_1010_sdpd[275];

    auto g_x_0_z_0_0_yy_y_xx = buffer_1010_sdpd[276];

    auto g_x_0_z_0_0_yy_y_xy = buffer_1010_sdpd[277];

    auto g_x_0_z_0_0_yy_y_xz = buffer_1010_sdpd[278];

    auto g_x_0_z_0_0_yy_y_yy = buffer_1010_sdpd[279];

    auto g_x_0_z_0_0_yy_y_yz = buffer_1010_sdpd[280];

    auto g_x_0_z_0_0_yy_y_zz = buffer_1010_sdpd[281];

    auto g_x_0_z_0_0_yy_z_xx = buffer_1010_sdpd[282];

    auto g_x_0_z_0_0_yy_z_xy = buffer_1010_sdpd[283];

    auto g_x_0_z_0_0_yy_z_xz = buffer_1010_sdpd[284];

    auto g_x_0_z_0_0_yy_z_yy = buffer_1010_sdpd[285];

    auto g_x_0_z_0_0_yy_z_yz = buffer_1010_sdpd[286];

    auto g_x_0_z_0_0_yy_z_zz = buffer_1010_sdpd[287];

    auto g_x_0_z_0_0_yz_x_xx = buffer_1010_sdpd[288];

    auto g_x_0_z_0_0_yz_x_xy = buffer_1010_sdpd[289];

    auto g_x_0_z_0_0_yz_x_xz = buffer_1010_sdpd[290];

    auto g_x_0_z_0_0_yz_x_yy = buffer_1010_sdpd[291];

    auto g_x_0_z_0_0_yz_x_yz = buffer_1010_sdpd[292];

    auto g_x_0_z_0_0_yz_x_zz = buffer_1010_sdpd[293];

    auto g_x_0_z_0_0_yz_y_xx = buffer_1010_sdpd[294];

    auto g_x_0_z_0_0_yz_y_xy = buffer_1010_sdpd[295];

    auto g_x_0_z_0_0_yz_y_xz = buffer_1010_sdpd[296];

    auto g_x_0_z_0_0_yz_y_yy = buffer_1010_sdpd[297];

    auto g_x_0_z_0_0_yz_y_yz = buffer_1010_sdpd[298];

    auto g_x_0_z_0_0_yz_y_zz = buffer_1010_sdpd[299];

    auto g_x_0_z_0_0_yz_z_xx = buffer_1010_sdpd[300];

    auto g_x_0_z_0_0_yz_z_xy = buffer_1010_sdpd[301];

    auto g_x_0_z_0_0_yz_z_xz = buffer_1010_sdpd[302];

    auto g_x_0_z_0_0_yz_z_yy = buffer_1010_sdpd[303];

    auto g_x_0_z_0_0_yz_z_yz = buffer_1010_sdpd[304];

    auto g_x_0_z_0_0_yz_z_zz = buffer_1010_sdpd[305];

    auto g_x_0_z_0_0_zz_x_xx = buffer_1010_sdpd[306];

    auto g_x_0_z_0_0_zz_x_xy = buffer_1010_sdpd[307];

    auto g_x_0_z_0_0_zz_x_xz = buffer_1010_sdpd[308];

    auto g_x_0_z_0_0_zz_x_yy = buffer_1010_sdpd[309];

    auto g_x_0_z_0_0_zz_x_yz = buffer_1010_sdpd[310];

    auto g_x_0_z_0_0_zz_x_zz = buffer_1010_sdpd[311];

    auto g_x_0_z_0_0_zz_y_xx = buffer_1010_sdpd[312];

    auto g_x_0_z_0_0_zz_y_xy = buffer_1010_sdpd[313];

    auto g_x_0_z_0_0_zz_y_xz = buffer_1010_sdpd[314];

    auto g_x_0_z_0_0_zz_y_yy = buffer_1010_sdpd[315];

    auto g_x_0_z_0_0_zz_y_yz = buffer_1010_sdpd[316];

    auto g_x_0_z_0_0_zz_y_zz = buffer_1010_sdpd[317];

    auto g_x_0_z_0_0_zz_z_xx = buffer_1010_sdpd[318];

    auto g_x_0_z_0_0_zz_z_xy = buffer_1010_sdpd[319];

    auto g_x_0_z_0_0_zz_z_xz = buffer_1010_sdpd[320];

    auto g_x_0_z_0_0_zz_z_yy = buffer_1010_sdpd[321];

    auto g_x_0_z_0_0_zz_z_yz = buffer_1010_sdpd[322];

    auto g_x_0_z_0_0_zz_z_zz = buffer_1010_sdpd[323];

    auto g_y_0_x_0_0_xx_x_xx = buffer_1010_sdpd[324];

    auto g_y_0_x_0_0_xx_x_xy = buffer_1010_sdpd[325];

    auto g_y_0_x_0_0_xx_x_xz = buffer_1010_sdpd[326];

    auto g_y_0_x_0_0_xx_x_yy = buffer_1010_sdpd[327];

    auto g_y_0_x_0_0_xx_x_yz = buffer_1010_sdpd[328];

    auto g_y_0_x_0_0_xx_x_zz = buffer_1010_sdpd[329];

    auto g_y_0_x_0_0_xx_y_xx = buffer_1010_sdpd[330];

    auto g_y_0_x_0_0_xx_y_xy = buffer_1010_sdpd[331];

    auto g_y_0_x_0_0_xx_y_xz = buffer_1010_sdpd[332];

    auto g_y_0_x_0_0_xx_y_yy = buffer_1010_sdpd[333];

    auto g_y_0_x_0_0_xx_y_yz = buffer_1010_sdpd[334];

    auto g_y_0_x_0_0_xx_y_zz = buffer_1010_sdpd[335];

    auto g_y_0_x_0_0_xx_z_xx = buffer_1010_sdpd[336];

    auto g_y_0_x_0_0_xx_z_xy = buffer_1010_sdpd[337];

    auto g_y_0_x_0_0_xx_z_xz = buffer_1010_sdpd[338];

    auto g_y_0_x_0_0_xx_z_yy = buffer_1010_sdpd[339];

    auto g_y_0_x_0_0_xx_z_yz = buffer_1010_sdpd[340];

    auto g_y_0_x_0_0_xx_z_zz = buffer_1010_sdpd[341];

    auto g_y_0_x_0_0_xy_x_xx = buffer_1010_sdpd[342];

    auto g_y_0_x_0_0_xy_x_xy = buffer_1010_sdpd[343];

    auto g_y_0_x_0_0_xy_x_xz = buffer_1010_sdpd[344];

    auto g_y_0_x_0_0_xy_x_yy = buffer_1010_sdpd[345];

    auto g_y_0_x_0_0_xy_x_yz = buffer_1010_sdpd[346];

    auto g_y_0_x_0_0_xy_x_zz = buffer_1010_sdpd[347];

    auto g_y_0_x_0_0_xy_y_xx = buffer_1010_sdpd[348];

    auto g_y_0_x_0_0_xy_y_xy = buffer_1010_sdpd[349];

    auto g_y_0_x_0_0_xy_y_xz = buffer_1010_sdpd[350];

    auto g_y_0_x_0_0_xy_y_yy = buffer_1010_sdpd[351];

    auto g_y_0_x_0_0_xy_y_yz = buffer_1010_sdpd[352];

    auto g_y_0_x_0_0_xy_y_zz = buffer_1010_sdpd[353];

    auto g_y_0_x_0_0_xy_z_xx = buffer_1010_sdpd[354];

    auto g_y_0_x_0_0_xy_z_xy = buffer_1010_sdpd[355];

    auto g_y_0_x_0_0_xy_z_xz = buffer_1010_sdpd[356];

    auto g_y_0_x_0_0_xy_z_yy = buffer_1010_sdpd[357];

    auto g_y_0_x_0_0_xy_z_yz = buffer_1010_sdpd[358];

    auto g_y_0_x_0_0_xy_z_zz = buffer_1010_sdpd[359];

    auto g_y_0_x_0_0_xz_x_xx = buffer_1010_sdpd[360];

    auto g_y_0_x_0_0_xz_x_xy = buffer_1010_sdpd[361];

    auto g_y_0_x_0_0_xz_x_xz = buffer_1010_sdpd[362];

    auto g_y_0_x_0_0_xz_x_yy = buffer_1010_sdpd[363];

    auto g_y_0_x_0_0_xz_x_yz = buffer_1010_sdpd[364];

    auto g_y_0_x_0_0_xz_x_zz = buffer_1010_sdpd[365];

    auto g_y_0_x_0_0_xz_y_xx = buffer_1010_sdpd[366];

    auto g_y_0_x_0_0_xz_y_xy = buffer_1010_sdpd[367];

    auto g_y_0_x_0_0_xz_y_xz = buffer_1010_sdpd[368];

    auto g_y_0_x_0_0_xz_y_yy = buffer_1010_sdpd[369];

    auto g_y_0_x_0_0_xz_y_yz = buffer_1010_sdpd[370];

    auto g_y_0_x_0_0_xz_y_zz = buffer_1010_sdpd[371];

    auto g_y_0_x_0_0_xz_z_xx = buffer_1010_sdpd[372];

    auto g_y_0_x_0_0_xz_z_xy = buffer_1010_sdpd[373];

    auto g_y_0_x_0_0_xz_z_xz = buffer_1010_sdpd[374];

    auto g_y_0_x_0_0_xz_z_yy = buffer_1010_sdpd[375];

    auto g_y_0_x_0_0_xz_z_yz = buffer_1010_sdpd[376];

    auto g_y_0_x_0_0_xz_z_zz = buffer_1010_sdpd[377];

    auto g_y_0_x_0_0_yy_x_xx = buffer_1010_sdpd[378];

    auto g_y_0_x_0_0_yy_x_xy = buffer_1010_sdpd[379];

    auto g_y_0_x_0_0_yy_x_xz = buffer_1010_sdpd[380];

    auto g_y_0_x_0_0_yy_x_yy = buffer_1010_sdpd[381];

    auto g_y_0_x_0_0_yy_x_yz = buffer_1010_sdpd[382];

    auto g_y_0_x_0_0_yy_x_zz = buffer_1010_sdpd[383];

    auto g_y_0_x_0_0_yy_y_xx = buffer_1010_sdpd[384];

    auto g_y_0_x_0_0_yy_y_xy = buffer_1010_sdpd[385];

    auto g_y_0_x_0_0_yy_y_xz = buffer_1010_sdpd[386];

    auto g_y_0_x_0_0_yy_y_yy = buffer_1010_sdpd[387];

    auto g_y_0_x_0_0_yy_y_yz = buffer_1010_sdpd[388];

    auto g_y_0_x_0_0_yy_y_zz = buffer_1010_sdpd[389];

    auto g_y_0_x_0_0_yy_z_xx = buffer_1010_sdpd[390];

    auto g_y_0_x_0_0_yy_z_xy = buffer_1010_sdpd[391];

    auto g_y_0_x_0_0_yy_z_xz = buffer_1010_sdpd[392];

    auto g_y_0_x_0_0_yy_z_yy = buffer_1010_sdpd[393];

    auto g_y_0_x_0_0_yy_z_yz = buffer_1010_sdpd[394];

    auto g_y_0_x_0_0_yy_z_zz = buffer_1010_sdpd[395];

    auto g_y_0_x_0_0_yz_x_xx = buffer_1010_sdpd[396];

    auto g_y_0_x_0_0_yz_x_xy = buffer_1010_sdpd[397];

    auto g_y_0_x_0_0_yz_x_xz = buffer_1010_sdpd[398];

    auto g_y_0_x_0_0_yz_x_yy = buffer_1010_sdpd[399];

    auto g_y_0_x_0_0_yz_x_yz = buffer_1010_sdpd[400];

    auto g_y_0_x_0_0_yz_x_zz = buffer_1010_sdpd[401];

    auto g_y_0_x_0_0_yz_y_xx = buffer_1010_sdpd[402];

    auto g_y_0_x_0_0_yz_y_xy = buffer_1010_sdpd[403];

    auto g_y_0_x_0_0_yz_y_xz = buffer_1010_sdpd[404];

    auto g_y_0_x_0_0_yz_y_yy = buffer_1010_sdpd[405];

    auto g_y_0_x_0_0_yz_y_yz = buffer_1010_sdpd[406];

    auto g_y_0_x_0_0_yz_y_zz = buffer_1010_sdpd[407];

    auto g_y_0_x_0_0_yz_z_xx = buffer_1010_sdpd[408];

    auto g_y_0_x_0_0_yz_z_xy = buffer_1010_sdpd[409];

    auto g_y_0_x_0_0_yz_z_xz = buffer_1010_sdpd[410];

    auto g_y_0_x_0_0_yz_z_yy = buffer_1010_sdpd[411];

    auto g_y_0_x_0_0_yz_z_yz = buffer_1010_sdpd[412];

    auto g_y_0_x_0_0_yz_z_zz = buffer_1010_sdpd[413];

    auto g_y_0_x_0_0_zz_x_xx = buffer_1010_sdpd[414];

    auto g_y_0_x_0_0_zz_x_xy = buffer_1010_sdpd[415];

    auto g_y_0_x_0_0_zz_x_xz = buffer_1010_sdpd[416];

    auto g_y_0_x_0_0_zz_x_yy = buffer_1010_sdpd[417];

    auto g_y_0_x_0_0_zz_x_yz = buffer_1010_sdpd[418];

    auto g_y_0_x_0_0_zz_x_zz = buffer_1010_sdpd[419];

    auto g_y_0_x_0_0_zz_y_xx = buffer_1010_sdpd[420];

    auto g_y_0_x_0_0_zz_y_xy = buffer_1010_sdpd[421];

    auto g_y_0_x_0_0_zz_y_xz = buffer_1010_sdpd[422];

    auto g_y_0_x_0_0_zz_y_yy = buffer_1010_sdpd[423];

    auto g_y_0_x_0_0_zz_y_yz = buffer_1010_sdpd[424];

    auto g_y_0_x_0_0_zz_y_zz = buffer_1010_sdpd[425];

    auto g_y_0_x_0_0_zz_z_xx = buffer_1010_sdpd[426];

    auto g_y_0_x_0_0_zz_z_xy = buffer_1010_sdpd[427];

    auto g_y_0_x_0_0_zz_z_xz = buffer_1010_sdpd[428];

    auto g_y_0_x_0_0_zz_z_yy = buffer_1010_sdpd[429];

    auto g_y_0_x_0_0_zz_z_yz = buffer_1010_sdpd[430];

    auto g_y_0_x_0_0_zz_z_zz = buffer_1010_sdpd[431];

    auto g_y_0_y_0_0_xx_x_xx = buffer_1010_sdpd[432];

    auto g_y_0_y_0_0_xx_x_xy = buffer_1010_sdpd[433];

    auto g_y_0_y_0_0_xx_x_xz = buffer_1010_sdpd[434];

    auto g_y_0_y_0_0_xx_x_yy = buffer_1010_sdpd[435];

    auto g_y_0_y_0_0_xx_x_yz = buffer_1010_sdpd[436];

    auto g_y_0_y_0_0_xx_x_zz = buffer_1010_sdpd[437];

    auto g_y_0_y_0_0_xx_y_xx = buffer_1010_sdpd[438];

    auto g_y_0_y_0_0_xx_y_xy = buffer_1010_sdpd[439];

    auto g_y_0_y_0_0_xx_y_xz = buffer_1010_sdpd[440];

    auto g_y_0_y_0_0_xx_y_yy = buffer_1010_sdpd[441];

    auto g_y_0_y_0_0_xx_y_yz = buffer_1010_sdpd[442];

    auto g_y_0_y_0_0_xx_y_zz = buffer_1010_sdpd[443];

    auto g_y_0_y_0_0_xx_z_xx = buffer_1010_sdpd[444];

    auto g_y_0_y_0_0_xx_z_xy = buffer_1010_sdpd[445];

    auto g_y_0_y_0_0_xx_z_xz = buffer_1010_sdpd[446];

    auto g_y_0_y_0_0_xx_z_yy = buffer_1010_sdpd[447];

    auto g_y_0_y_0_0_xx_z_yz = buffer_1010_sdpd[448];

    auto g_y_0_y_0_0_xx_z_zz = buffer_1010_sdpd[449];

    auto g_y_0_y_0_0_xy_x_xx = buffer_1010_sdpd[450];

    auto g_y_0_y_0_0_xy_x_xy = buffer_1010_sdpd[451];

    auto g_y_0_y_0_0_xy_x_xz = buffer_1010_sdpd[452];

    auto g_y_0_y_0_0_xy_x_yy = buffer_1010_sdpd[453];

    auto g_y_0_y_0_0_xy_x_yz = buffer_1010_sdpd[454];

    auto g_y_0_y_0_0_xy_x_zz = buffer_1010_sdpd[455];

    auto g_y_0_y_0_0_xy_y_xx = buffer_1010_sdpd[456];

    auto g_y_0_y_0_0_xy_y_xy = buffer_1010_sdpd[457];

    auto g_y_0_y_0_0_xy_y_xz = buffer_1010_sdpd[458];

    auto g_y_0_y_0_0_xy_y_yy = buffer_1010_sdpd[459];

    auto g_y_0_y_0_0_xy_y_yz = buffer_1010_sdpd[460];

    auto g_y_0_y_0_0_xy_y_zz = buffer_1010_sdpd[461];

    auto g_y_0_y_0_0_xy_z_xx = buffer_1010_sdpd[462];

    auto g_y_0_y_0_0_xy_z_xy = buffer_1010_sdpd[463];

    auto g_y_0_y_0_0_xy_z_xz = buffer_1010_sdpd[464];

    auto g_y_0_y_0_0_xy_z_yy = buffer_1010_sdpd[465];

    auto g_y_0_y_0_0_xy_z_yz = buffer_1010_sdpd[466];

    auto g_y_0_y_0_0_xy_z_zz = buffer_1010_sdpd[467];

    auto g_y_0_y_0_0_xz_x_xx = buffer_1010_sdpd[468];

    auto g_y_0_y_0_0_xz_x_xy = buffer_1010_sdpd[469];

    auto g_y_0_y_0_0_xz_x_xz = buffer_1010_sdpd[470];

    auto g_y_0_y_0_0_xz_x_yy = buffer_1010_sdpd[471];

    auto g_y_0_y_0_0_xz_x_yz = buffer_1010_sdpd[472];

    auto g_y_0_y_0_0_xz_x_zz = buffer_1010_sdpd[473];

    auto g_y_0_y_0_0_xz_y_xx = buffer_1010_sdpd[474];

    auto g_y_0_y_0_0_xz_y_xy = buffer_1010_sdpd[475];

    auto g_y_0_y_0_0_xz_y_xz = buffer_1010_sdpd[476];

    auto g_y_0_y_0_0_xz_y_yy = buffer_1010_sdpd[477];

    auto g_y_0_y_0_0_xz_y_yz = buffer_1010_sdpd[478];

    auto g_y_0_y_0_0_xz_y_zz = buffer_1010_sdpd[479];

    auto g_y_0_y_0_0_xz_z_xx = buffer_1010_sdpd[480];

    auto g_y_0_y_0_0_xz_z_xy = buffer_1010_sdpd[481];

    auto g_y_0_y_0_0_xz_z_xz = buffer_1010_sdpd[482];

    auto g_y_0_y_0_0_xz_z_yy = buffer_1010_sdpd[483];

    auto g_y_0_y_0_0_xz_z_yz = buffer_1010_sdpd[484];

    auto g_y_0_y_0_0_xz_z_zz = buffer_1010_sdpd[485];

    auto g_y_0_y_0_0_yy_x_xx = buffer_1010_sdpd[486];

    auto g_y_0_y_0_0_yy_x_xy = buffer_1010_sdpd[487];

    auto g_y_0_y_0_0_yy_x_xz = buffer_1010_sdpd[488];

    auto g_y_0_y_0_0_yy_x_yy = buffer_1010_sdpd[489];

    auto g_y_0_y_0_0_yy_x_yz = buffer_1010_sdpd[490];

    auto g_y_0_y_0_0_yy_x_zz = buffer_1010_sdpd[491];

    auto g_y_0_y_0_0_yy_y_xx = buffer_1010_sdpd[492];

    auto g_y_0_y_0_0_yy_y_xy = buffer_1010_sdpd[493];

    auto g_y_0_y_0_0_yy_y_xz = buffer_1010_sdpd[494];

    auto g_y_0_y_0_0_yy_y_yy = buffer_1010_sdpd[495];

    auto g_y_0_y_0_0_yy_y_yz = buffer_1010_sdpd[496];

    auto g_y_0_y_0_0_yy_y_zz = buffer_1010_sdpd[497];

    auto g_y_0_y_0_0_yy_z_xx = buffer_1010_sdpd[498];

    auto g_y_0_y_0_0_yy_z_xy = buffer_1010_sdpd[499];

    auto g_y_0_y_0_0_yy_z_xz = buffer_1010_sdpd[500];

    auto g_y_0_y_0_0_yy_z_yy = buffer_1010_sdpd[501];

    auto g_y_0_y_0_0_yy_z_yz = buffer_1010_sdpd[502];

    auto g_y_0_y_0_0_yy_z_zz = buffer_1010_sdpd[503];

    auto g_y_0_y_0_0_yz_x_xx = buffer_1010_sdpd[504];

    auto g_y_0_y_0_0_yz_x_xy = buffer_1010_sdpd[505];

    auto g_y_0_y_0_0_yz_x_xz = buffer_1010_sdpd[506];

    auto g_y_0_y_0_0_yz_x_yy = buffer_1010_sdpd[507];

    auto g_y_0_y_0_0_yz_x_yz = buffer_1010_sdpd[508];

    auto g_y_0_y_0_0_yz_x_zz = buffer_1010_sdpd[509];

    auto g_y_0_y_0_0_yz_y_xx = buffer_1010_sdpd[510];

    auto g_y_0_y_0_0_yz_y_xy = buffer_1010_sdpd[511];

    auto g_y_0_y_0_0_yz_y_xz = buffer_1010_sdpd[512];

    auto g_y_0_y_0_0_yz_y_yy = buffer_1010_sdpd[513];

    auto g_y_0_y_0_0_yz_y_yz = buffer_1010_sdpd[514];

    auto g_y_0_y_0_0_yz_y_zz = buffer_1010_sdpd[515];

    auto g_y_0_y_0_0_yz_z_xx = buffer_1010_sdpd[516];

    auto g_y_0_y_0_0_yz_z_xy = buffer_1010_sdpd[517];

    auto g_y_0_y_0_0_yz_z_xz = buffer_1010_sdpd[518];

    auto g_y_0_y_0_0_yz_z_yy = buffer_1010_sdpd[519];

    auto g_y_0_y_0_0_yz_z_yz = buffer_1010_sdpd[520];

    auto g_y_0_y_0_0_yz_z_zz = buffer_1010_sdpd[521];

    auto g_y_0_y_0_0_zz_x_xx = buffer_1010_sdpd[522];

    auto g_y_0_y_0_0_zz_x_xy = buffer_1010_sdpd[523];

    auto g_y_0_y_0_0_zz_x_xz = buffer_1010_sdpd[524];

    auto g_y_0_y_0_0_zz_x_yy = buffer_1010_sdpd[525];

    auto g_y_0_y_0_0_zz_x_yz = buffer_1010_sdpd[526];

    auto g_y_0_y_0_0_zz_x_zz = buffer_1010_sdpd[527];

    auto g_y_0_y_0_0_zz_y_xx = buffer_1010_sdpd[528];

    auto g_y_0_y_0_0_zz_y_xy = buffer_1010_sdpd[529];

    auto g_y_0_y_0_0_zz_y_xz = buffer_1010_sdpd[530];

    auto g_y_0_y_0_0_zz_y_yy = buffer_1010_sdpd[531];

    auto g_y_0_y_0_0_zz_y_yz = buffer_1010_sdpd[532];

    auto g_y_0_y_0_0_zz_y_zz = buffer_1010_sdpd[533];

    auto g_y_0_y_0_0_zz_z_xx = buffer_1010_sdpd[534];

    auto g_y_0_y_0_0_zz_z_xy = buffer_1010_sdpd[535];

    auto g_y_0_y_0_0_zz_z_xz = buffer_1010_sdpd[536];

    auto g_y_0_y_0_0_zz_z_yy = buffer_1010_sdpd[537];

    auto g_y_0_y_0_0_zz_z_yz = buffer_1010_sdpd[538];

    auto g_y_0_y_0_0_zz_z_zz = buffer_1010_sdpd[539];

    auto g_y_0_z_0_0_xx_x_xx = buffer_1010_sdpd[540];

    auto g_y_0_z_0_0_xx_x_xy = buffer_1010_sdpd[541];

    auto g_y_0_z_0_0_xx_x_xz = buffer_1010_sdpd[542];

    auto g_y_0_z_0_0_xx_x_yy = buffer_1010_sdpd[543];

    auto g_y_0_z_0_0_xx_x_yz = buffer_1010_sdpd[544];

    auto g_y_0_z_0_0_xx_x_zz = buffer_1010_sdpd[545];

    auto g_y_0_z_0_0_xx_y_xx = buffer_1010_sdpd[546];

    auto g_y_0_z_0_0_xx_y_xy = buffer_1010_sdpd[547];

    auto g_y_0_z_0_0_xx_y_xz = buffer_1010_sdpd[548];

    auto g_y_0_z_0_0_xx_y_yy = buffer_1010_sdpd[549];

    auto g_y_0_z_0_0_xx_y_yz = buffer_1010_sdpd[550];

    auto g_y_0_z_0_0_xx_y_zz = buffer_1010_sdpd[551];

    auto g_y_0_z_0_0_xx_z_xx = buffer_1010_sdpd[552];

    auto g_y_0_z_0_0_xx_z_xy = buffer_1010_sdpd[553];

    auto g_y_0_z_0_0_xx_z_xz = buffer_1010_sdpd[554];

    auto g_y_0_z_0_0_xx_z_yy = buffer_1010_sdpd[555];

    auto g_y_0_z_0_0_xx_z_yz = buffer_1010_sdpd[556];

    auto g_y_0_z_0_0_xx_z_zz = buffer_1010_sdpd[557];

    auto g_y_0_z_0_0_xy_x_xx = buffer_1010_sdpd[558];

    auto g_y_0_z_0_0_xy_x_xy = buffer_1010_sdpd[559];

    auto g_y_0_z_0_0_xy_x_xz = buffer_1010_sdpd[560];

    auto g_y_0_z_0_0_xy_x_yy = buffer_1010_sdpd[561];

    auto g_y_0_z_0_0_xy_x_yz = buffer_1010_sdpd[562];

    auto g_y_0_z_0_0_xy_x_zz = buffer_1010_sdpd[563];

    auto g_y_0_z_0_0_xy_y_xx = buffer_1010_sdpd[564];

    auto g_y_0_z_0_0_xy_y_xy = buffer_1010_sdpd[565];

    auto g_y_0_z_0_0_xy_y_xz = buffer_1010_sdpd[566];

    auto g_y_0_z_0_0_xy_y_yy = buffer_1010_sdpd[567];

    auto g_y_0_z_0_0_xy_y_yz = buffer_1010_sdpd[568];

    auto g_y_0_z_0_0_xy_y_zz = buffer_1010_sdpd[569];

    auto g_y_0_z_0_0_xy_z_xx = buffer_1010_sdpd[570];

    auto g_y_0_z_0_0_xy_z_xy = buffer_1010_sdpd[571];

    auto g_y_0_z_0_0_xy_z_xz = buffer_1010_sdpd[572];

    auto g_y_0_z_0_0_xy_z_yy = buffer_1010_sdpd[573];

    auto g_y_0_z_0_0_xy_z_yz = buffer_1010_sdpd[574];

    auto g_y_0_z_0_0_xy_z_zz = buffer_1010_sdpd[575];

    auto g_y_0_z_0_0_xz_x_xx = buffer_1010_sdpd[576];

    auto g_y_0_z_0_0_xz_x_xy = buffer_1010_sdpd[577];

    auto g_y_0_z_0_0_xz_x_xz = buffer_1010_sdpd[578];

    auto g_y_0_z_0_0_xz_x_yy = buffer_1010_sdpd[579];

    auto g_y_0_z_0_0_xz_x_yz = buffer_1010_sdpd[580];

    auto g_y_0_z_0_0_xz_x_zz = buffer_1010_sdpd[581];

    auto g_y_0_z_0_0_xz_y_xx = buffer_1010_sdpd[582];

    auto g_y_0_z_0_0_xz_y_xy = buffer_1010_sdpd[583];

    auto g_y_0_z_0_0_xz_y_xz = buffer_1010_sdpd[584];

    auto g_y_0_z_0_0_xz_y_yy = buffer_1010_sdpd[585];

    auto g_y_0_z_0_0_xz_y_yz = buffer_1010_sdpd[586];

    auto g_y_0_z_0_0_xz_y_zz = buffer_1010_sdpd[587];

    auto g_y_0_z_0_0_xz_z_xx = buffer_1010_sdpd[588];

    auto g_y_0_z_0_0_xz_z_xy = buffer_1010_sdpd[589];

    auto g_y_0_z_0_0_xz_z_xz = buffer_1010_sdpd[590];

    auto g_y_0_z_0_0_xz_z_yy = buffer_1010_sdpd[591];

    auto g_y_0_z_0_0_xz_z_yz = buffer_1010_sdpd[592];

    auto g_y_0_z_0_0_xz_z_zz = buffer_1010_sdpd[593];

    auto g_y_0_z_0_0_yy_x_xx = buffer_1010_sdpd[594];

    auto g_y_0_z_0_0_yy_x_xy = buffer_1010_sdpd[595];

    auto g_y_0_z_0_0_yy_x_xz = buffer_1010_sdpd[596];

    auto g_y_0_z_0_0_yy_x_yy = buffer_1010_sdpd[597];

    auto g_y_0_z_0_0_yy_x_yz = buffer_1010_sdpd[598];

    auto g_y_0_z_0_0_yy_x_zz = buffer_1010_sdpd[599];

    auto g_y_0_z_0_0_yy_y_xx = buffer_1010_sdpd[600];

    auto g_y_0_z_0_0_yy_y_xy = buffer_1010_sdpd[601];

    auto g_y_0_z_0_0_yy_y_xz = buffer_1010_sdpd[602];

    auto g_y_0_z_0_0_yy_y_yy = buffer_1010_sdpd[603];

    auto g_y_0_z_0_0_yy_y_yz = buffer_1010_sdpd[604];

    auto g_y_0_z_0_0_yy_y_zz = buffer_1010_sdpd[605];

    auto g_y_0_z_0_0_yy_z_xx = buffer_1010_sdpd[606];

    auto g_y_0_z_0_0_yy_z_xy = buffer_1010_sdpd[607];

    auto g_y_0_z_0_0_yy_z_xz = buffer_1010_sdpd[608];

    auto g_y_0_z_0_0_yy_z_yy = buffer_1010_sdpd[609];

    auto g_y_0_z_0_0_yy_z_yz = buffer_1010_sdpd[610];

    auto g_y_0_z_0_0_yy_z_zz = buffer_1010_sdpd[611];

    auto g_y_0_z_0_0_yz_x_xx = buffer_1010_sdpd[612];

    auto g_y_0_z_0_0_yz_x_xy = buffer_1010_sdpd[613];

    auto g_y_0_z_0_0_yz_x_xz = buffer_1010_sdpd[614];

    auto g_y_0_z_0_0_yz_x_yy = buffer_1010_sdpd[615];

    auto g_y_0_z_0_0_yz_x_yz = buffer_1010_sdpd[616];

    auto g_y_0_z_0_0_yz_x_zz = buffer_1010_sdpd[617];

    auto g_y_0_z_0_0_yz_y_xx = buffer_1010_sdpd[618];

    auto g_y_0_z_0_0_yz_y_xy = buffer_1010_sdpd[619];

    auto g_y_0_z_0_0_yz_y_xz = buffer_1010_sdpd[620];

    auto g_y_0_z_0_0_yz_y_yy = buffer_1010_sdpd[621];

    auto g_y_0_z_0_0_yz_y_yz = buffer_1010_sdpd[622];

    auto g_y_0_z_0_0_yz_y_zz = buffer_1010_sdpd[623];

    auto g_y_0_z_0_0_yz_z_xx = buffer_1010_sdpd[624];

    auto g_y_0_z_0_0_yz_z_xy = buffer_1010_sdpd[625];

    auto g_y_0_z_0_0_yz_z_xz = buffer_1010_sdpd[626];

    auto g_y_0_z_0_0_yz_z_yy = buffer_1010_sdpd[627];

    auto g_y_0_z_0_0_yz_z_yz = buffer_1010_sdpd[628];

    auto g_y_0_z_0_0_yz_z_zz = buffer_1010_sdpd[629];

    auto g_y_0_z_0_0_zz_x_xx = buffer_1010_sdpd[630];

    auto g_y_0_z_0_0_zz_x_xy = buffer_1010_sdpd[631];

    auto g_y_0_z_0_0_zz_x_xz = buffer_1010_sdpd[632];

    auto g_y_0_z_0_0_zz_x_yy = buffer_1010_sdpd[633];

    auto g_y_0_z_0_0_zz_x_yz = buffer_1010_sdpd[634];

    auto g_y_0_z_0_0_zz_x_zz = buffer_1010_sdpd[635];

    auto g_y_0_z_0_0_zz_y_xx = buffer_1010_sdpd[636];

    auto g_y_0_z_0_0_zz_y_xy = buffer_1010_sdpd[637];

    auto g_y_0_z_0_0_zz_y_xz = buffer_1010_sdpd[638];

    auto g_y_0_z_0_0_zz_y_yy = buffer_1010_sdpd[639];

    auto g_y_0_z_0_0_zz_y_yz = buffer_1010_sdpd[640];

    auto g_y_0_z_0_0_zz_y_zz = buffer_1010_sdpd[641];

    auto g_y_0_z_0_0_zz_z_xx = buffer_1010_sdpd[642];

    auto g_y_0_z_0_0_zz_z_xy = buffer_1010_sdpd[643];

    auto g_y_0_z_0_0_zz_z_xz = buffer_1010_sdpd[644];

    auto g_y_0_z_0_0_zz_z_yy = buffer_1010_sdpd[645];

    auto g_y_0_z_0_0_zz_z_yz = buffer_1010_sdpd[646];

    auto g_y_0_z_0_0_zz_z_zz = buffer_1010_sdpd[647];

    auto g_z_0_x_0_0_xx_x_xx = buffer_1010_sdpd[648];

    auto g_z_0_x_0_0_xx_x_xy = buffer_1010_sdpd[649];

    auto g_z_0_x_0_0_xx_x_xz = buffer_1010_sdpd[650];

    auto g_z_0_x_0_0_xx_x_yy = buffer_1010_sdpd[651];

    auto g_z_0_x_0_0_xx_x_yz = buffer_1010_sdpd[652];

    auto g_z_0_x_0_0_xx_x_zz = buffer_1010_sdpd[653];

    auto g_z_0_x_0_0_xx_y_xx = buffer_1010_sdpd[654];

    auto g_z_0_x_0_0_xx_y_xy = buffer_1010_sdpd[655];

    auto g_z_0_x_0_0_xx_y_xz = buffer_1010_sdpd[656];

    auto g_z_0_x_0_0_xx_y_yy = buffer_1010_sdpd[657];

    auto g_z_0_x_0_0_xx_y_yz = buffer_1010_sdpd[658];

    auto g_z_0_x_0_0_xx_y_zz = buffer_1010_sdpd[659];

    auto g_z_0_x_0_0_xx_z_xx = buffer_1010_sdpd[660];

    auto g_z_0_x_0_0_xx_z_xy = buffer_1010_sdpd[661];

    auto g_z_0_x_0_0_xx_z_xz = buffer_1010_sdpd[662];

    auto g_z_0_x_0_0_xx_z_yy = buffer_1010_sdpd[663];

    auto g_z_0_x_0_0_xx_z_yz = buffer_1010_sdpd[664];

    auto g_z_0_x_0_0_xx_z_zz = buffer_1010_sdpd[665];

    auto g_z_0_x_0_0_xy_x_xx = buffer_1010_sdpd[666];

    auto g_z_0_x_0_0_xy_x_xy = buffer_1010_sdpd[667];

    auto g_z_0_x_0_0_xy_x_xz = buffer_1010_sdpd[668];

    auto g_z_0_x_0_0_xy_x_yy = buffer_1010_sdpd[669];

    auto g_z_0_x_0_0_xy_x_yz = buffer_1010_sdpd[670];

    auto g_z_0_x_0_0_xy_x_zz = buffer_1010_sdpd[671];

    auto g_z_0_x_0_0_xy_y_xx = buffer_1010_sdpd[672];

    auto g_z_0_x_0_0_xy_y_xy = buffer_1010_sdpd[673];

    auto g_z_0_x_0_0_xy_y_xz = buffer_1010_sdpd[674];

    auto g_z_0_x_0_0_xy_y_yy = buffer_1010_sdpd[675];

    auto g_z_0_x_0_0_xy_y_yz = buffer_1010_sdpd[676];

    auto g_z_0_x_0_0_xy_y_zz = buffer_1010_sdpd[677];

    auto g_z_0_x_0_0_xy_z_xx = buffer_1010_sdpd[678];

    auto g_z_0_x_0_0_xy_z_xy = buffer_1010_sdpd[679];

    auto g_z_0_x_0_0_xy_z_xz = buffer_1010_sdpd[680];

    auto g_z_0_x_0_0_xy_z_yy = buffer_1010_sdpd[681];

    auto g_z_0_x_0_0_xy_z_yz = buffer_1010_sdpd[682];

    auto g_z_0_x_0_0_xy_z_zz = buffer_1010_sdpd[683];

    auto g_z_0_x_0_0_xz_x_xx = buffer_1010_sdpd[684];

    auto g_z_0_x_0_0_xz_x_xy = buffer_1010_sdpd[685];

    auto g_z_0_x_0_0_xz_x_xz = buffer_1010_sdpd[686];

    auto g_z_0_x_0_0_xz_x_yy = buffer_1010_sdpd[687];

    auto g_z_0_x_0_0_xz_x_yz = buffer_1010_sdpd[688];

    auto g_z_0_x_0_0_xz_x_zz = buffer_1010_sdpd[689];

    auto g_z_0_x_0_0_xz_y_xx = buffer_1010_sdpd[690];

    auto g_z_0_x_0_0_xz_y_xy = buffer_1010_sdpd[691];

    auto g_z_0_x_0_0_xz_y_xz = buffer_1010_sdpd[692];

    auto g_z_0_x_0_0_xz_y_yy = buffer_1010_sdpd[693];

    auto g_z_0_x_0_0_xz_y_yz = buffer_1010_sdpd[694];

    auto g_z_0_x_0_0_xz_y_zz = buffer_1010_sdpd[695];

    auto g_z_0_x_0_0_xz_z_xx = buffer_1010_sdpd[696];

    auto g_z_0_x_0_0_xz_z_xy = buffer_1010_sdpd[697];

    auto g_z_0_x_0_0_xz_z_xz = buffer_1010_sdpd[698];

    auto g_z_0_x_0_0_xz_z_yy = buffer_1010_sdpd[699];

    auto g_z_0_x_0_0_xz_z_yz = buffer_1010_sdpd[700];

    auto g_z_0_x_0_0_xz_z_zz = buffer_1010_sdpd[701];

    auto g_z_0_x_0_0_yy_x_xx = buffer_1010_sdpd[702];

    auto g_z_0_x_0_0_yy_x_xy = buffer_1010_sdpd[703];

    auto g_z_0_x_0_0_yy_x_xz = buffer_1010_sdpd[704];

    auto g_z_0_x_0_0_yy_x_yy = buffer_1010_sdpd[705];

    auto g_z_0_x_0_0_yy_x_yz = buffer_1010_sdpd[706];

    auto g_z_0_x_0_0_yy_x_zz = buffer_1010_sdpd[707];

    auto g_z_0_x_0_0_yy_y_xx = buffer_1010_sdpd[708];

    auto g_z_0_x_0_0_yy_y_xy = buffer_1010_sdpd[709];

    auto g_z_0_x_0_0_yy_y_xz = buffer_1010_sdpd[710];

    auto g_z_0_x_0_0_yy_y_yy = buffer_1010_sdpd[711];

    auto g_z_0_x_0_0_yy_y_yz = buffer_1010_sdpd[712];

    auto g_z_0_x_0_0_yy_y_zz = buffer_1010_sdpd[713];

    auto g_z_0_x_0_0_yy_z_xx = buffer_1010_sdpd[714];

    auto g_z_0_x_0_0_yy_z_xy = buffer_1010_sdpd[715];

    auto g_z_0_x_0_0_yy_z_xz = buffer_1010_sdpd[716];

    auto g_z_0_x_0_0_yy_z_yy = buffer_1010_sdpd[717];

    auto g_z_0_x_0_0_yy_z_yz = buffer_1010_sdpd[718];

    auto g_z_0_x_0_0_yy_z_zz = buffer_1010_sdpd[719];

    auto g_z_0_x_0_0_yz_x_xx = buffer_1010_sdpd[720];

    auto g_z_0_x_0_0_yz_x_xy = buffer_1010_sdpd[721];

    auto g_z_0_x_0_0_yz_x_xz = buffer_1010_sdpd[722];

    auto g_z_0_x_0_0_yz_x_yy = buffer_1010_sdpd[723];

    auto g_z_0_x_0_0_yz_x_yz = buffer_1010_sdpd[724];

    auto g_z_0_x_0_0_yz_x_zz = buffer_1010_sdpd[725];

    auto g_z_0_x_0_0_yz_y_xx = buffer_1010_sdpd[726];

    auto g_z_0_x_0_0_yz_y_xy = buffer_1010_sdpd[727];

    auto g_z_0_x_0_0_yz_y_xz = buffer_1010_sdpd[728];

    auto g_z_0_x_0_0_yz_y_yy = buffer_1010_sdpd[729];

    auto g_z_0_x_0_0_yz_y_yz = buffer_1010_sdpd[730];

    auto g_z_0_x_0_0_yz_y_zz = buffer_1010_sdpd[731];

    auto g_z_0_x_0_0_yz_z_xx = buffer_1010_sdpd[732];

    auto g_z_0_x_0_0_yz_z_xy = buffer_1010_sdpd[733];

    auto g_z_0_x_0_0_yz_z_xz = buffer_1010_sdpd[734];

    auto g_z_0_x_0_0_yz_z_yy = buffer_1010_sdpd[735];

    auto g_z_0_x_0_0_yz_z_yz = buffer_1010_sdpd[736];

    auto g_z_0_x_0_0_yz_z_zz = buffer_1010_sdpd[737];

    auto g_z_0_x_0_0_zz_x_xx = buffer_1010_sdpd[738];

    auto g_z_0_x_0_0_zz_x_xy = buffer_1010_sdpd[739];

    auto g_z_0_x_0_0_zz_x_xz = buffer_1010_sdpd[740];

    auto g_z_0_x_0_0_zz_x_yy = buffer_1010_sdpd[741];

    auto g_z_0_x_0_0_zz_x_yz = buffer_1010_sdpd[742];

    auto g_z_0_x_0_0_zz_x_zz = buffer_1010_sdpd[743];

    auto g_z_0_x_0_0_zz_y_xx = buffer_1010_sdpd[744];

    auto g_z_0_x_0_0_zz_y_xy = buffer_1010_sdpd[745];

    auto g_z_0_x_0_0_zz_y_xz = buffer_1010_sdpd[746];

    auto g_z_0_x_0_0_zz_y_yy = buffer_1010_sdpd[747];

    auto g_z_0_x_0_0_zz_y_yz = buffer_1010_sdpd[748];

    auto g_z_0_x_0_0_zz_y_zz = buffer_1010_sdpd[749];

    auto g_z_0_x_0_0_zz_z_xx = buffer_1010_sdpd[750];

    auto g_z_0_x_0_0_zz_z_xy = buffer_1010_sdpd[751];

    auto g_z_0_x_0_0_zz_z_xz = buffer_1010_sdpd[752];

    auto g_z_0_x_0_0_zz_z_yy = buffer_1010_sdpd[753];

    auto g_z_0_x_0_0_zz_z_yz = buffer_1010_sdpd[754];

    auto g_z_0_x_0_0_zz_z_zz = buffer_1010_sdpd[755];

    auto g_z_0_y_0_0_xx_x_xx = buffer_1010_sdpd[756];

    auto g_z_0_y_0_0_xx_x_xy = buffer_1010_sdpd[757];

    auto g_z_0_y_0_0_xx_x_xz = buffer_1010_sdpd[758];

    auto g_z_0_y_0_0_xx_x_yy = buffer_1010_sdpd[759];

    auto g_z_0_y_0_0_xx_x_yz = buffer_1010_sdpd[760];

    auto g_z_0_y_0_0_xx_x_zz = buffer_1010_sdpd[761];

    auto g_z_0_y_0_0_xx_y_xx = buffer_1010_sdpd[762];

    auto g_z_0_y_0_0_xx_y_xy = buffer_1010_sdpd[763];

    auto g_z_0_y_0_0_xx_y_xz = buffer_1010_sdpd[764];

    auto g_z_0_y_0_0_xx_y_yy = buffer_1010_sdpd[765];

    auto g_z_0_y_0_0_xx_y_yz = buffer_1010_sdpd[766];

    auto g_z_0_y_0_0_xx_y_zz = buffer_1010_sdpd[767];

    auto g_z_0_y_0_0_xx_z_xx = buffer_1010_sdpd[768];

    auto g_z_0_y_0_0_xx_z_xy = buffer_1010_sdpd[769];

    auto g_z_0_y_0_0_xx_z_xz = buffer_1010_sdpd[770];

    auto g_z_0_y_0_0_xx_z_yy = buffer_1010_sdpd[771];

    auto g_z_0_y_0_0_xx_z_yz = buffer_1010_sdpd[772];

    auto g_z_0_y_0_0_xx_z_zz = buffer_1010_sdpd[773];

    auto g_z_0_y_0_0_xy_x_xx = buffer_1010_sdpd[774];

    auto g_z_0_y_0_0_xy_x_xy = buffer_1010_sdpd[775];

    auto g_z_0_y_0_0_xy_x_xz = buffer_1010_sdpd[776];

    auto g_z_0_y_0_0_xy_x_yy = buffer_1010_sdpd[777];

    auto g_z_0_y_0_0_xy_x_yz = buffer_1010_sdpd[778];

    auto g_z_0_y_0_0_xy_x_zz = buffer_1010_sdpd[779];

    auto g_z_0_y_0_0_xy_y_xx = buffer_1010_sdpd[780];

    auto g_z_0_y_0_0_xy_y_xy = buffer_1010_sdpd[781];

    auto g_z_0_y_0_0_xy_y_xz = buffer_1010_sdpd[782];

    auto g_z_0_y_0_0_xy_y_yy = buffer_1010_sdpd[783];

    auto g_z_0_y_0_0_xy_y_yz = buffer_1010_sdpd[784];

    auto g_z_0_y_0_0_xy_y_zz = buffer_1010_sdpd[785];

    auto g_z_0_y_0_0_xy_z_xx = buffer_1010_sdpd[786];

    auto g_z_0_y_0_0_xy_z_xy = buffer_1010_sdpd[787];

    auto g_z_0_y_0_0_xy_z_xz = buffer_1010_sdpd[788];

    auto g_z_0_y_0_0_xy_z_yy = buffer_1010_sdpd[789];

    auto g_z_0_y_0_0_xy_z_yz = buffer_1010_sdpd[790];

    auto g_z_0_y_0_0_xy_z_zz = buffer_1010_sdpd[791];

    auto g_z_0_y_0_0_xz_x_xx = buffer_1010_sdpd[792];

    auto g_z_0_y_0_0_xz_x_xy = buffer_1010_sdpd[793];

    auto g_z_0_y_0_0_xz_x_xz = buffer_1010_sdpd[794];

    auto g_z_0_y_0_0_xz_x_yy = buffer_1010_sdpd[795];

    auto g_z_0_y_0_0_xz_x_yz = buffer_1010_sdpd[796];

    auto g_z_0_y_0_0_xz_x_zz = buffer_1010_sdpd[797];

    auto g_z_0_y_0_0_xz_y_xx = buffer_1010_sdpd[798];

    auto g_z_0_y_0_0_xz_y_xy = buffer_1010_sdpd[799];

    auto g_z_0_y_0_0_xz_y_xz = buffer_1010_sdpd[800];

    auto g_z_0_y_0_0_xz_y_yy = buffer_1010_sdpd[801];

    auto g_z_0_y_0_0_xz_y_yz = buffer_1010_sdpd[802];

    auto g_z_0_y_0_0_xz_y_zz = buffer_1010_sdpd[803];

    auto g_z_0_y_0_0_xz_z_xx = buffer_1010_sdpd[804];

    auto g_z_0_y_0_0_xz_z_xy = buffer_1010_sdpd[805];

    auto g_z_0_y_0_0_xz_z_xz = buffer_1010_sdpd[806];

    auto g_z_0_y_0_0_xz_z_yy = buffer_1010_sdpd[807];

    auto g_z_0_y_0_0_xz_z_yz = buffer_1010_sdpd[808];

    auto g_z_0_y_0_0_xz_z_zz = buffer_1010_sdpd[809];

    auto g_z_0_y_0_0_yy_x_xx = buffer_1010_sdpd[810];

    auto g_z_0_y_0_0_yy_x_xy = buffer_1010_sdpd[811];

    auto g_z_0_y_0_0_yy_x_xz = buffer_1010_sdpd[812];

    auto g_z_0_y_0_0_yy_x_yy = buffer_1010_sdpd[813];

    auto g_z_0_y_0_0_yy_x_yz = buffer_1010_sdpd[814];

    auto g_z_0_y_0_0_yy_x_zz = buffer_1010_sdpd[815];

    auto g_z_0_y_0_0_yy_y_xx = buffer_1010_sdpd[816];

    auto g_z_0_y_0_0_yy_y_xy = buffer_1010_sdpd[817];

    auto g_z_0_y_0_0_yy_y_xz = buffer_1010_sdpd[818];

    auto g_z_0_y_0_0_yy_y_yy = buffer_1010_sdpd[819];

    auto g_z_0_y_0_0_yy_y_yz = buffer_1010_sdpd[820];

    auto g_z_0_y_0_0_yy_y_zz = buffer_1010_sdpd[821];

    auto g_z_0_y_0_0_yy_z_xx = buffer_1010_sdpd[822];

    auto g_z_0_y_0_0_yy_z_xy = buffer_1010_sdpd[823];

    auto g_z_0_y_0_0_yy_z_xz = buffer_1010_sdpd[824];

    auto g_z_0_y_0_0_yy_z_yy = buffer_1010_sdpd[825];

    auto g_z_0_y_0_0_yy_z_yz = buffer_1010_sdpd[826];

    auto g_z_0_y_0_0_yy_z_zz = buffer_1010_sdpd[827];

    auto g_z_0_y_0_0_yz_x_xx = buffer_1010_sdpd[828];

    auto g_z_0_y_0_0_yz_x_xy = buffer_1010_sdpd[829];

    auto g_z_0_y_0_0_yz_x_xz = buffer_1010_sdpd[830];

    auto g_z_0_y_0_0_yz_x_yy = buffer_1010_sdpd[831];

    auto g_z_0_y_0_0_yz_x_yz = buffer_1010_sdpd[832];

    auto g_z_0_y_0_0_yz_x_zz = buffer_1010_sdpd[833];

    auto g_z_0_y_0_0_yz_y_xx = buffer_1010_sdpd[834];

    auto g_z_0_y_0_0_yz_y_xy = buffer_1010_sdpd[835];

    auto g_z_0_y_0_0_yz_y_xz = buffer_1010_sdpd[836];

    auto g_z_0_y_0_0_yz_y_yy = buffer_1010_sdpd[837];

    auto g_z_0_y_0_0_yz_y_yz = buffer_1010_sdpd[838];

    auto g_z_0_y_0_0_yz_y_zz = buffer_1010_sdpd[839];

    auto g_z_0_y_0_0_yz_z_xx = buffer_1010_sdpd[840];

    auto g_z_0_y_0_0_yz_z_xy = buffer_1010_sdpd[841];

    auto g_z_0_y_0_0_yz_z_xz = buffer_1010_sdpd[842];

    auto g_z_0_y_0_0_yz_z_yy = buffer_1010_sdpd[843];

    auto g_z_0_y_0_0_yz_z_yz = buffer_1010_sdpd[844];

    auto g_z_0_y_0_0_yz_z_zz = buffer_1010_sdpd[845];

    auto g_z_0_y_0_0_zz_x_xx = buffer_1010_sdpd[846];

    auto g_z_0_y_0_0_zz_x_xy = buffer_1010_sdpd[847];

    auto g_z_0_y_0_0_zz_x_xz = buffer_1010_sdpd[848];

    auto g_z_0_y_0_0_zz_x_yy = buffer_1010_sdpd[849];

    auto g_z_0_y_0_0_zz_x_yz = buffer_1010_sdpd[850];

    auto g_z_0_y_0_0_zz_x_zz = buffer_1010_sdpd[851];

    auto g_z_0_y_0_0_zz_y_xx = buffer_1010_sdpd[852];

    auto g_z_0_y_0_0_zz_y_xy = buffer_1010_sdpd[853];

    auto g_z_0_y_0_0_zz_y_xz = buffer_1010_sdpd[854];

    auto g_z_0_y_0_0_zz_y_yy = buffer_1010_sdpd[855];

    auto g_z_0_y_0_0_zz_y_yz = buffer_1010_sdpd[856];

    auto g_z_0_y_0_0_zz_y_zz = buffer_1010_sdpd[857];

    auto g_z_0_y_0_0_zz_z_xx = buffer_1010_sdpd[858];

    auto g_z_0_y_0_0_zz_z_xy = buffer_1010_sdpd[859];

    auto g_z_0_y_0_0_zz_z_xz = buffer_1010_sdpd[860];

    auto g_z_0_y_0_0_zz_z_yy = buffer_1010_sdpd[861];

    auto g_z_0_y_0_0_zz_z_yz = buffer_1010_sdpd[862];

    auto g_z_0_y_0_0_zz_z_zz = buffer_1010_sdpd[863];

    auto g_z_0_z_0_0_xx_x_xx = buffer_1010_sdpd[864];

    auto g_z_0_z_0_0_xx_x_xy = buffer_1010_sdpd[865];

    auto g_z_0_z_0_0_xx_x_xz = buffer_1010_sdpd[866];

    auto g_z_0_z_0_0_xx_x_yy = buffer_1010_sdpd[867];

    auto g_z_0_z_0_0_xx_x_yz = buffer_1010_sdpd[868];

    auto g_z_0_z_0_0_xx_x_zz = buffer_1010_sdpd[869];

    auto g_z_0_z_0_0_xx_y_xx = buffer_1010_sdpd[870];

    auto g_z_0_z_0_0_xx_y_xy = buffer_1010_sdpd[871];

    auto g_z_0_z_0_0_xx_y_xz = buffer_1010_sdpd[872];

    auto g_z_0_z_0_0_xx_y_yy = buffer_1010_sdpd[873];

    auto g_z_0_z_0_0_xx_y_yz = buffer_1010_sdpd[874];

    auto g_z_0_z_0_0_xx_y_zz = buffer_1010_sdpd[875];

    auto g_z_0_z_0_0_xx_z_xx = buffer_1010_sdpd[876];

    auto g_z_0_z_0_0_xx_z_xy = buffer_1010_sdpd[877];

    auto g_z_0_z_0_0_xx_z_xz = buffer_1010_sdpd[878];

    auto g_z_0_z_0_0_xx_z_yy = buffer_1010_sdpd[879];

    auto g_z_0_z_0_0_xx_z_yz = buffer_1010_sdpd[880];

    auto g_z_0_z_0_0_xx_z_zz = buffer_1010_sdpd[881];

    auto g_z_0_z_0_0_xy_x_xx = buffer_1010_sdpd[882];

    auto g_z_0_z_0_0_xy_x_xy = buffer_1010_sdpd[883];

    auto g_z_0_z_0_0_xy_x_xz = buffer_1010_sdpd[884];

    auto g_z_0_z_0_0_xy_x_yy = buffer_1010_sdpd[885];

    auto g_z_0_z_0_0_xy_x_yz = buffer_1010_sdpd[886];

    auto g_z_0_z_0_0_xy_x_zz = buffer_1010_sdpd[887];

    auto g_z_0_z_0_0_xy_y_xx = buffer_1010_sdpd[888];

    auto g_z_0_z_0_0_xy_y_xy = buffer_1010_sdpd[889];

    auto g_z_0_z_0_0_xy_y_xz = buffer_1010_sdpd[890];

    auto g_z_0_z_0_0_xy_y_yy = buffer_1010_sdpd[891];

    auto g_z_0_z_0_0_xy_y_yz = buffer_1010_sdpd[892];

    auto g_z_0_z_0_0_xy_y_zz = buffer_1010_sdpd[893];

    auto g_z_0_z_0_0_xy_z_xx = buffer_1010_sdpd[894];

    auto g_z_0_z_0_0_xy_z_xy = buffer_1010_sdpd[895];

    auto g_z_0_z_0_0_xy_z_xz = buffer_1010_sdpd[896];

    auto g_z_0_z_0_0_xy_z_yy = buffer_1010_sdpd[897];

    auto g_z_0_z_0_0_xy_z_yz = buffer_1010_sdpd[898];

    auto g_z_0_z_0_0_xy_z_zz = buffer_1010_sdpd[899];

    auto g_z_0_z_0_0_xz_x_xx = buffer_1010_sdpd[900];

    auto g_z_0_z_0_0_xz_x_xy = buffer_1010_sdpd[901];

    auto g_z_0_z_0_0_xz_x_xz = buffer_1010_sdpd[902];

    auto g_z_0_z_0_0_xz_x_yy = buffer_1010_sdpd[903];

    auto g_z_0_z_0_0_xz_x_yz = buffer_1010_sdpd[904];

    auto g_z_0_z_0_0_xz_x_zz = buffer_1010_sdpd[905];

    auto g_z_0_z_0_0_xz_y_xx = buffer_1010_sdpd[906];

    auto g_z_0_z_0_0_xz_y_xy = buffer_1010_sdpd[907];

    auto g_z_0_z_0_0_xz_y_xz = buffer_1010_sdpd[908];

    auto g_z_0_z_0_0_xz_y_yy = buffer_1010_sdpd[909];

    auto g_z_0_z_0_0_xz_y_yz = buffer_1010_sdpd[910];

    auto g_z_0_z_0_0_xz_y_zz = buffer_1010_sdpd[911];

    auto g_z_0_z_0_0_xz_z_xx = buffer_1010_sdpd[912];

    auto g_z_0_z_0_0_xz_z_xy = buffer_1010_sdpd[913];

    auto g_z_0_z_0_0_xz_z_xz = buffer_1010_sdpd[914];

    auto g_z_0_z_0_0_xz_z_yy = buffer_1010_sdpd[915];

    auto g_z_0_z_0_0_xz_z_yz = buffer_1010_sdpd[916];

    auto g_z_0_z_0_0_xz_z_zz = buffer_1010_sdpd[917];

    auto g_z_0_z_0_0_yy_x_xx = buffer_1010_sdpd[918];

    auto g_z_0_z_0_0_yy_x_xy = buffer_1010_sdpd[919];

    auto g_z_0_z_0_0_yy_x_xz = buffer_1010_sdpd[920];

    auto g_z_0_z_0_0_yy_x_yy = buffer_1010_sdpd[921];

    auto g_z_0_z_0_0_yy_x_yz = buffer_1010_sdpd[922];

    auto g_z_0_z_0_0_yy_x_zz = buffer_1010_sdpd[923];

    auto g_z_0_z_0_0_yy_y_xx = buffer_1010_sdpd[924];

    auto g_z_0_z_0_0_yy_y_xy = buffer_1010_sdpd[925];

    auto g_z_0_z_0_0_yy_y_xz = buffer_1010_sdpd[926];

    auto g_z_0_z_0_0_yy_y_yy = buffer_1010_sdpd[927];

    auto g_z_0_z_0_0_yy_y_yz = buffer_1010_sdpd[928];

    auto g_z_0_z_0_0_yy_y_zz = buffer_1010_sdpd[929];

    auto g_z_0_z_0_0_yy_z_xx = buffer_1010_sdpd[930];

    auto g_z_0_z_0_0_yy_z_xy = buffer_1010_sdpd[931];

    auto g_z_0_z_0_0_yy_z_xz = buffer_1010_sdpd[932];

    auto g_z_0_z_0_0_yy_z_yy = buffer_1010_sdpd[933];

    auto g_z_0_z_0_0_yy_z_yz = buffer_1010_sdpd[934];

    auto g_z_0_z_0_0_yy_z_zz = buffer_1010_sdpd[935];

    auto g_z_0_z_0_0_yz_x_xx = buffer_1010_sdpd[936];

    auto g_z_0_z_0_0_yz_x_xy = buffer_1010_sdpd[937];

    auto g_z_0_z_0_0_yz_x_xz = buffer_1010_sdpd[938];

    auto g_z_0_z_0_0_yz_x_yy = buffer_1010_sdpd[939];

    auto g_z_0_z_0_0_yz_x_yz = buffer_1010_sdpd[940];

    auto g_z_0_z_0_0_yz_x_zz = buffer_1010_sdpd[941];

    auto g_z_0_z_0_0_yz_y_xx = buffer_1010_sdpd[942];

    auto g_z_0_z_0_0_yz_y_xy = buffer_1010_sdpd[943];

    auto g_z_0_z_0_0_yz_y_xz = buffer_1010_sdpd[944];

    auto g_z_0_z_0_0_yz_y_yy = buffer_1010_sdpd[945];

    auto g_z_0_z_0_0_yz_y_yz = buffer_1010_sdpd[946];

    auto g_z_0_z_0_0_yz_y_zz = buffer_1010_sdpd[947];

    auto g_z_0_z_0_0_yz_z_xx = buffer_1010_sdpd[948];

    auto g_z_0_z_0_0_yz_z_xy = buffer_1010_sdpd[949];

    auto g_z_0_z_0_0_yz_z_xz = buffer_1010_sdpd[950];

    auto g_z_0_z_0_0_yz_z_yy = buffer_1010_sdpd[951];

    auto g_z_0_z_0_0_yz_z_yz = buffer_1010_sdpd[952];

    auto g_z_0_z_0_0_yz_z_zz = buffer_1010_sdpd[953];

    auto g_z_0_z_0_0_zz_x_xx = buffer_1010_sdpd[954];

    auto g_z_0_z_0_0_zz_x_xy = buffer_1010_sdpd[955];

    auto g_z_0_z_0_0_zz_x_xz = buffer_1010_sdpd[956];

    auto g_z_0_z_0_0_zz_x_yy = buffer_1010_sdpd[957];

    auto g_z_0_z_0_0_zz_x_yz = buffer_1010_sdpd[958];

    auto g_z_0_z_0_0_zz_x_zz = buffer_1010_sdpd[959];

    auto g_z_0_z_0_0_zz_y_xx = buffer_1010_sdpd[960];

    auto g_z_0_z_0_0_zz_y_xy = buffer_1010_sdpd[961];

    auto g_z_0_z_0_0_zz_y_xz = buffer_1010_sdpd[962];

    auto g_z_0_z_0_0_zz_y_yy = buffer_1010_sdpd[963];

    auto g_z_0_z_0_0_zz_y_yz = buffer_1010_sdpd[964];

    auto g_z_0_z_0_0_zz_y_zz = buffer_1010_sdpd[965];

    auto g_z_0_z_0_0_zz_z_xx = buffer_1010_sdpd[966];

    auto g_z_0_z_0_0_zz_z_xy = buffer_1010_sdpd[967];

    auto g_z_0_z_0_0_zz_z_xz = buffer_1010_sdpd[968];

    auto g_z_0_z_0_0_zz_z_yy = buffer_1010_sdpd[969];

    auto g_z_0_z_0_0_zz_z_yz = buffer_1010_sdpd[970];

    auto g_z_0_z_0_0_zz_z_zz = buffer_1010_sdpd[971];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_x_0_0_xx_x_xx, g_x_0_x_0_0_xx_x_xy, g_x_0_x_0_0_xx_x_xz, g_x_0_x_0_0_xx_x_yy, g_x_0_x_0_0_xx_x_yz, g_x_0_x_0_0_xx_x_zz, g_x_xx_0_xx, g_x_xx_0_xy, g_x_xx_0_xz, g_x_xx_0_yy, g_x_xx_0_yz, g_x_xx_0_zz, g_x_xx_xx_xx, g_x_xx_xx_xy, g_x_xx_xx_xz, g_x_xx_xx_yy, g_x_xx_xx_yz, g_x_xx_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xx_x_xx[i] = -2.0 * g_x_xx_0_xx[i] * a_exp + 4.0 * g_x_xx_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_x_xy[i] = -2.0 * g_x_xx_0_xy[i] * a_exp + 4.0 * g_x_xx_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_x_xz[i] = -2.0 * g_x_xx_0_xz[i] * a_exp + 4.0 * g_x_xx_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_x_yy[i] = -2.0 * g_x_xx_0_yy[i] * a_exp + 4.0 * g_x_xx_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_x_yz[i] = -2.0 * g_x_xx_0_yz[i] * a_exp + 4.0 * g_x_xx_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_x_zz[i] = -2.0 * g_x_xx_0_zz[i] * a_exp + 4.0 * g_x_xx_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_x_0_0_xx_y_xx, g_x_0_x_0_0_xx_y_xy, g_x_0_x_0_0_xx_y_xz, g_x_0_x_0_0_xx_y_yy, g_x_0_x_0_0_xx_y_yz, g_x_0_x_0_0_xx_y_zz, g_x_xx_xy_xx, g_x_xx_xy_xy, g_x_xx_xy_xz, g_x_xx_xy_yy, g_x_xx_xy_yz, g_x_xx_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xx_y_xx[i] = 4.0 * g_x_xx_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_y_xy[i] = 4.0 * g_x_xx_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_y_xz[i] = 4.0 * g_x_xx_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_y_yy[i] = 4.0 * g_x_xx_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_y_yz[i] = 4.0 * g_x_xx_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_y_zz[i] = 4.0 * g_x_xx_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_x_0_0_xx_z_xx, g_x_0_x_0_0_xx_z_xy, g_x_0_x_0_0_xx_z_xz, g_x_0_x_0_0_xx_z_yy, g_x_0_x_0_0_xx_z_yz, g_x_0_x_0_0_xx_z_zz, g_x_xx_xz_xx, g_x_xx_xz_xy, g_x_xx_xz_xz, g_x_xx_xz_yy, g_x_xx_xz_yz, g_x_xx_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xx_z_xx[i] = 4.0 * g_x_xx_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_z_xy[i] = 4.0 * g_x_xx_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_z_xz[i] = 4.0 * g_x_xx_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_z_yy[i] = 4.0 * g_x_xx_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_z_yz[i] = 4.0 * g_x_xx_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_z_zz[i] = 4.0 * g_x_xx_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_x_0_0_xy_x_xx, g_x_0_x_0_0_xy_x_xy, g_x_0_x_0_0_xy_x_xz, g_x_0_x_0_0_xy_x_yy, g_x_0_x_0_0_xy_x_yz, g_x_0_x_0_0_xy_x_zz, g_x_xy_0_xx, g_x_xy_0_xy, g_x_xy_0_xz, g_x_xy_0_yy, g_x_xy_0_yz, g_x_xy_0_zz, g_x_xy_xx_xx, g_x_xy_xx_xy, g_x_xy_xx_xz, g_x_xy_xx_yy, g_x_xy_xx_yz, g_x_xy_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xy_x_xx[i] = -2.0 * g_x_xy_0_xx[i] * a_exp + 4.0 * g_x_xy_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_x_xy[i] = -2.0 * g_x_xy_0_xy[i] * a_exp + 4.0 * g_x_xy_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_x_xz[i] = -2.0 * g_x_xy_0_xz[i] * a_exp + 4.0 * g_x_xy_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_x_yy[i] = -2.0 * g_x_xy_0_yy[i] * a_exp + 4.0 * g_x_xy_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_x_yz[i] = -2.0 * g_x_xy_0_yz[i] * a_exp + 4.0 * g_x_xy_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_x_zz[i] = -2.0 * g_x_xy_0_zz[i] * a_exp + 4.0 * g_x_xy_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_x_0_0_xy_y_xx, g_x_0_x_0_0_xy_y_xy, g_x_0_x_0_0_xy_y_xz, g_x_0_x_0_0_xy_y_yy, g_x_0_x_0_0_xy_y_yz, g_x_0_x_0_0_xy_y_zz, g_x_xy_xy_xx, g_x_xy_xy_xy, g_x_xy_xy_xz, g_x_xy_xy_yy, g_x_xy_xy_yz, g_x_xy_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xy_y_xx[i] = 4.0 * g_x_xy_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_y_xy[i] = 4.0 * g_x_xy_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_y_xz[i] = 4.0 * g_x_xy_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_y_yy[i] = 4.0 * g_x_xy_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_y_yz[i] = 4.0 * g_x_xy_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_y_zz[i] = 4.0 * g_x_xy_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_x_0_0_xy_z_xx, g_x_0_x_0_0_xy_z_xy, g_x_0_x_0_0_xy_z_xz, g_x_0_x_0_0_xy_z_yy, g_x_0_x_0_0_xy_z_yz, g_x_0_x_0_0_xy_z_zz, g_x_xy_xz_xx, g_x_xy_xz_xy, g_x_xy_xz_xz, g_x_xy_xz_yy, g_x_xy_xz_yz, g_x_xy_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xy_z_xx[i] = 4.0 * g_x_xy_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_z_xy[i] = 4.0 * g_x_xy_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_z_xz[i] = 4.0 * g_x_xy_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_z_yy[i] = 4.0 * g_x_xy_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_z_yz[i] = 4.0 * g_x_xy_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_z_zz[i] = 4.0 * g_x_xy_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_x_0_0_xz_x_xx, g_x_0_x_0_0_xz_x_xy, g_x_0_x_0_0_xz_x_xz, g_x_0_x_0_0_xz_x_yy, g_x_0_x_0_0_xz_x_yz, g_x_0_x_0_0_xz_x_zz, g_x_xz_0_xx, g_x_xz_0_xy, g_x_xz_0_xz, g_x_xz_0_yy, g_x_xz_0_yz, g_x_xz_0_zz, g_x_xz_xx_xx, g_x_xz_xx_xy, g_x_xz_xx_xz, g_x_xz_xx_yy, g_x_xz_xx_yz, g_x_xz_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xz_x_xx[i] = -2.0 * g_x_xz_0_xx[i] * a_exp + 4.0 * g_x_xz_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_x_xy[i] = -2.0 * g_x_xz_0_xy[i] * a_exp + 4.0 * g_x_xz_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_x_xz[i] = -2.0 * g_x_xz_0_xz[i] * a_exp + 4.0 * g_x_xz_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_x_yy[i] = -2.0 * g_x_xz_0_yy[i] * a_exp + 4.0 * g_x_xz_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_x_yz[i] = -2.0 * g_x_xz_0_yz[i] * a_exp + 4.0 * g_x_xz_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_x_zz[i] = -2.0 * g_x_xz_0_zz[i] * a_exp + 4.0 * g_x_xz_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_x_0_0_xz_y_xx, g_x_0_x_0_0_xz_y_xy, g_x_0_x_0_0_xz_y_xz, g_x_0_x_0_0_xz_y_yy, g_x_0_x_0_0_xz_y_yz, g_x_0_x_0_0_xz_y_zz, g_x_xz_xy_xx, g_x_xz_xy_xy, g_x_xz_xy_xz, g_x_xz_xy_yy, g_x_xz_xy_yz, g_x_xz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xz_y_xx[i] = 4.0 * g_x_xz_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_y_xy[i] = 4.0 * g_x_xz_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_y_xz[i] = 4.0 * g_x_xz_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_y_yy[i] = 4.0 * g_x_xz_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_y_yz[i] = 4.0 * g_x_xz_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_y_zz[i] = 4.0 * g_x_xz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_x_0_0_xz_z_xx, g_x_0_x_0_0_xz_z_xy, g_x_0_x_0_0_xz_z_xz, g_x_0_x_0_0_xz_z_yy, g_x_0_x_0_0_xz_z_yz, g_x_0_x_0_0_xz_z_zz, g_x_xz_xz_xx, g_x_xz_xz_xy, g_x_xz_xz_xz, g_x_xz_xz_yy, g_x_xz_xz_yz, g_x_xz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xz_z_xx[i] = 4.0 * g_x_xz_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_z_xy[i] = 4.0 * g_x_xz_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_z_xz[i] = 4.0 * g_x_xz_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_z_yy[i] = 4.0 * g_x_xz_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_z_yz[i] = 4.0 * g_x_xz_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_z_zz[i] = 4.0 * g_x_xz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_0_x_0_0_yy_x_xx, g_x_0_x_0_0_yy_x_xy, g_x_0_x_0_0_yy_x_xz, g_x_0_x_0_0_yy_x_yy, g_x_0_x_0_0_yy_x_yz, g_x_0_x_0_0_yy_x_zz, g_x_yy_0_xx, g_x_yy_0_xy, g_x_yy_0_xz, g_x_yy_0_yy, g_x_yy_0_yz, g_x_yy_0_zz, g_x_yy_xx_xx, g_x_yy_xx_xy, g_x_yy_xx_xz, g_x_yy_xx_yy, g_x_yy_xx_yz, g_x_yy_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yy_x_xx[i] = -2.0 * g_x_yy_0_xx[i] * a_exp + 4.0 * g_x_yy_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_x_xy[i] = -2.0 * g_x_yy_0_xy[i] * a_exp + 4.0 * g_x_yy_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_x_xz[i] = -2.0 * g_x_yy_0_xz[i] * a_exp + 4.0 * g_x_yy_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_x_yy[i] = -2.0 * g_x_yy_0_yy[i] * a_exp + 4.0 * g_x_yy_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_x_yz[i] = -2.0 * g_x_yy_0_yz[i] * a_exp + 4.0 * g_x_yy_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_x_zz[i] = -2.0 * g_x_yy_0_zz[i] * a_exp + 4.0 * g_x_yy_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_0_x_0_0_yy_y_xx, g_x_0_x_0_0_yy_y_xy, g_x_0_x_0_0_yy_y_xz, g_x_0_x_0_0_yy_y_yy, g_x_0_x_0_0_yy_y_yz, g_x_0_x_0_0_yy_y_zz, g_x_yy_xy_xx, g_x_yy_xy_xy, g_x_yy_xy_xz, g_x_yy_xy_yy, g_x_yy_xy_yz, g_x_yy_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yy_y_xx[i] = 4.0 * g_x_yy_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_y_xy[i] = 4.0 * g_x_yy_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_y_xz[i] = 4.0 * g_x_yy_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_y_yy[i] = 4.0 * g_x_yy_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_y_yz[i] = 4.0 * g_x_yy_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_y_zz[i] = 4.0 * g_x_yy_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_0_x_0_0_yy_z_xx, g_x_0_x_0_0_yy_z_xy, g_x_0_x_0_0_yy_z_xz, g_x_0_x_0_0_yy_z_yy, g_x_0_x_0_0_yy_z_yz, g_x_0_x_0_0_yy_z_zz, g_x_yy_xz_xx, g_x_yy_xz_xy, g_x_yy_xz_xz, g_x_yy_xz_yy, g_x_yy_xz_yz, g_x_yy_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yy_z_xx[i] = 4.0 * g_x_yy_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_z_xy[i] = 4.0 * g_x_yy_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_z_xz[i] = 4.0 * g_x_yy_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_z_yy[i] = 4.0 * g_x_yy_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_z_yz[i] = 4.0 * g_x_yy_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_z_zz[i] = 4.0 * g_x_yy_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_x_0_0_yz_x_xx, g_x_0_x_0_0_yz_x_xy, g_x_0_x_0_0_yz_x_xz, g_x_0_x_0_0_yz_x_yy, g_x_0_x_0_0_yz_x_yz, g_x_0_x_0_0_yz_x_zz, g_x_yz_0_xx, g_x_yz_0_xy, g_x_yz_0_xz, g_x_yz_0_yy, g_x_yz_0_yz, g_x_yz_0_zz, g_x_yz_xx_xx, g_x_yz_xx_xy, g_x_yz_xx_xz, g_x_yz_xx_yy, g_x_yz_xx_yz, g_x_yz_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yz_x_xx[i] = -2.0 * g_x_yz_0_xx[i] * a_exp + 4.0 * g_x_yz_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_x_xy[i] = -2.0 * g_x_yz_0_xy[i] * a_exp + 4.0 * g_x_yz_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_x_xz[i] = -2.0 * g_x_yz_0_xz[i] * a_exp + 4.0 * g_x_yz_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_x_yy[i] = -2.0 * g_x_yz_0_yy[i] * a_exp + 4.0 * g_x_yz_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_x_yz[i] = -2.0 * g_x_yz_0_yz[i] * a_exp + 4.0 * g_x_yz_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_x_zz[i] = -2.0 * g_x_yz_0_zz[i] * a_exp + 4.0 * g_x_yz_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_x_0_0_yz_y_xx, g_x_0_x_0_0_yz_y_xy, g_x_0_x_0_0_yz_y_xz, g_x_0_x_0_0_yz_y_yy, g_x_0_x_0_0_yz_y_yz, g_x_0_x_0_0_yz_y_zz, g_x_yz_xy_xx, g_x_yz_xy_xy, g_x_yz_xy_xz, g_x_yz_xy_yy, g_x_yz_xy_yz, g_x_yz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yz_y_xx[i] = 4.0 * g_x_yz_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_y_xy[i] = 4.0 * g_x_yz_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_y_xz[i] = 4.0 * g_x_yz_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_y_yy[i] = 4.0 * g_x_yz_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_y_yz[i] = 4.0 * g_x_yz_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_y_zz[i] = 4.0 * g_x_yz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_x_0_0_yz_z_xx, g_x_0_x_0_0_yz_z_xy, g_x_0_x_0_0_yz_z_xz, g_x_0_x_0_0_yz_z_yy, g_x_0_x_0_0_yz_z_yz, g_x_0_x_0_0_yz_z_zz, g_x_yz_xz_xx, g_x_yz_xz_xy, g_x_yz_xz_xz, g_x_yz_xz_yy, g_x_yz_xz_yz, g_x_yz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yz_z_xx[i] = 4.0 * g_x_yz_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_z_xy[i] = 4.0 * g_x_yz_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_z_xz[i] = 4.0 * g_x_yz_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_z_yy[i] = 4.0 * g_x_yz_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_z_yz[i] = 4.0 * g_x_yz_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_z_zz[i] = 4.0 * g_x_yz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_x_0_0_zz_x_xx, g_x_0_x_0_0_zz_x_xy, g_x_0_x_0_0_zz_x_xz, g_x_0_x_0_0_zz_x_yy, g_x_0_x_0_0_zz_x_yz, g_x_0_x_0_0_zz_x_zz, g_x_zz_0_xx, g_x_zz_0_xy, g_x_zz_0_xz, g_x_zz_0_yy, g_x_zz_0_yz, g_x_zz_0_zz, g_x_zz_xx_xx, g_x_zz_xx_xy, g_x_zz_xx_xz, g_x_zz_xx_yy, g_x_zz_xx_yz, g_x_zz_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_zz_x_xx[i] = -2.0 * g_x_zz_0_xx[i] * a_exp + 4.0 * g_x_zz_xx_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_x_xy[i] = -2.0 * g_x_zz_0_xy[i] * a_exp + 4.0 * g_x_zz_xx_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_x_xz[i] = -2.0 * g_x_zz_0_xz[i] * a_exp + 4.0 * g_x_zz_xx_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_x_yy[i] = -2.0 * g_x_zz_0_yy[i] * a_exp + 4.0 * g_x_zz_xx_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_x_yz[i] = -2.0 * g_x_zz_0_yz[i] * a_exp + 4.0 * g_x_zz_xx_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_x_zz[i] = -2.0 * g_x_zz_0_zz[i] * a_exp + 4.0 * g_x_zz_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_x_0_0_zz_y_xx, g_x_0_x_0_0_zz_y_xy, g_x_0_x_0_0_zz_y_xz, g_x_0_x_0_0_zz_y_yy, g_x_0_x_0_0_zz_y_yz, g_x_0_x_0_0_zz_y_zz, g_x_zz_xy_xx, g_x_zz_xy_xy, g_x_zz_xy_xz, g_x_zz_xy_yy, g_x_zz_xy_yz, g_x_zz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_zz_y_xx[i] = 4.0 * g_x_zz_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_y_xy[i] = 4.0 * g_x_zz_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_y_xz[i] = 4.0 * g_x_zz_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_y_yy[i] = 4.0 * g_x_zz_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_y_yz[i] = 4.0 * g_x_zz_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_y_zz[i] = 4.0 * g_x_zz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_x_0_0_zz_z_xx, g_x_0_x_0_0_zz_z_xy, g_x_0_x_0_0_zz_z_xz, g_x_0_x_0_0_zz_z_yy, g_x_0_x_0_0_zz_z_yz, g_x_0_x_0_0_zz_z_zz, g_x_zz_xz_xx, g_x_zz_xz_xy, g_x_zz_xz_xz, g_x_zz_xz_yy, g_x_zz_xz_yz, g_x_zz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_zz_z_xx[i] = 4.0 * g_x_zz_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_z_xy[i] = 4.0 * g_x_zz_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_z_xz[i] = 4.0 * g_x_zz_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_z_yy[i] = 4.0 * g_x_zz_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_z_yz[i] = 4.0 * g_x_zz_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_z_zz[i] = 4.0 * g_x_zz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_0_y_0_0_xx_x_xx, g_x_0_y_0_0_xx_x_xy, g_x_0_y_0_0_xx_x_xz, g_x_0_y_0_0_xx_x_yy, g_x_0_y_0_0_xx_x_yz, g_x_0_y_0_0_xx_x_zz, g_x_xx_xy_xx, g_x_xx_xy_xy, g_x_xx_xy_xz, g_x_xx_xy_yy, g_x_xx_xy_yz, g_x_xx_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xx_x_xx[i] = 4.0 * g_x_xx_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_x_xy[i] = 4.0 * g_x_xx_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_x_xz[i] = 4.0 * g_x_xx_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_x_yy[i] = 4.0 * g_x_xx_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_x_yz[i] = 4.0 * g_x_xx_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_x_zz[i] = 4.0 * g_x_xx_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_0_y_0_0_xx_y_xx, g_x_0_y_0_0_xx_y_xy, g_x_0_y_0_0_xx_y_xz, g_x_0_y_0_0_xx_y_yy, g_x_0_y_0_0_xx_y_yz, g_x_0_y_0_0_xx_y_zz, g_x_xx_0_xx, g_x_xx_0_xy, g_x_xx_0_xz, g_x_xx_0_yy, g_x_xx_0_yz, g_x_xx_0_zz, g_x_xx_yy_xx, g_x_xx_yy_xy, g_x_xx_yy_xz, g_x_xx_yy_yy, g_x_xx_yy_yz, g_x_xx_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xx_y_xx[i] = -2.0 * g_x_xx_0_xx[i] * a_exp + 4.0 * g_x_xx_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_y_xy[i] = -2.0 * g_x_xx_0_xy[i] * a_exp + 4.0 * g_x_xx_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_y_xz[i] = -2.0 * g_x_xx_0_xz[i] * a_exp + 4.0 * g_x_xx_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_y_yy[i] = -2.0 * g_x_xx_0_yy[i] * a_exp + 4.0 * g_x_xx_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_y_yz[i] = -2.0 * g_x_xx_0_yz[i] * a_exp + 4.0 * g_x_xx_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_y_zz[i] = -2.0 * g_x_xx_0_zz[i] * a_exp + 4.0 * g_x_xx_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_0_y_0_0_xx_z_xx, g_x_0_y_0_0_xx_z_xy, g_x_0_y_0_0_xx_z_xz, g_x_0_y_0_0_xx_z_yy, g_x_0_y_0_0_xx_z_yz, g_x_0_y_0_0_xx_z_zz, g_x_xx_yz_xx, g_x_xx_yz_xy, g_x_xx_yz_xz, g_x_xx_yz_yy, g_x_xx_yz_yz, g_x_xx_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xx_z_xx[i] = 4.0 * g_x_xx_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_z_xy[i] = 4.0 * g_x_xx_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_z_xz[i] = 4.0 * g_x_xx_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_z_yy[i] = 4.0 * g_x_xx_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_z_yz[i] = 4.0 * g_x_xx_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_z_zz[i] = 4.0 * g_x_xx_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_0_y_0_0_xy_x_xx, g_x_0_y_0_0_xy_x_xy, g_x_0_y_0_0_xy_x_xz, g_x_0_y_0_0_xy_x_yy, g_x_0_y_0_0_xy_x_yz, g_x_0_y_0_0_xy_x_zz, g_x_xy_xy_xx, g_x_xy_xy_xy, g_x_xy_xy_xz, g_x_xy_xy_yy, g_x_xy_xy_yz, g_x_xy_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xy_x_xx[i] = 4.0 * g_x_xy_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_x_xy[i] = 4.0 * g_x_xy_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_x_xz[i] = 4.0 * g_x_xy_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_x_yy[i] = 4.0 * g_x_xy_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_x_yz[i] = 4.0 * g_x_xy_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_x_zz[i] = 4.0 * g_x_xy_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_0_y_0_0_xy_y_xx, g_x_0_y_0_0_xy_y_xy, g_x_0_y_0_0_xy_y_xz, g_x_0_y_0_0_xy_y_yy, g_x_0_y_0_0_xy_y_yz, g_x_0_y_0_0_xy_y_zz, g_x_xy_0_xx, g_x_xy_0_xy, g_x_xy_0_xz, g_x_xy_0_yy, g_x_xy_0_yz, g_x_xy_0_zz, g_x_xy_yy_xx, g_x_xy_yy_xy, g_x_xy_yy_xz, g_x_xy_yy_yy, g_x_xy_yy_yz, g_x_xy_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xy_y_xx[i] = -2.0 * g_x_xy_0_xx[i] * a_exp + 4.0 * g_x_xy_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_y_xy[i] = -2.0 * g_x_xy_0_xy[i] * a_exp + 4.0 * g_x_xy_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_y_xz[i] = -2.0 * g_x_xy_0_xz[i] * a_exp + 4.0 * g_x_xy_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_y_yy[i] = -2.0 * g_x_xy_0_yy[i] * a_exp + 4.0 * g_x_xy_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_y_yz[i] = -2.0 * g_x_xy_0_yz[i] * a_exp + 4.0 * g_x_xy_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_y_zz[i] = -2.0 * g_x_xy_0_zz[i] * a_exp + 4.0 * g_x_xy_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_0_y_0_0_xy_z_xx, g_x_0_y_0_0_xy_z_xy, g_x_0_y_0_0_xy_z_xz, g_x_0_y_0_0_xy_z_yy, g_x_0_y_0_0_xy_z_yz, g_x_0_y_0_0_xy_z_zz, g_x_xy_yz_xx, g_x_xy_yz_xy, g_x_xy_yz_xz, g_x_xy_yz_yy, g_x_xy_yz_yz, g_x_xy_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xy_z_xx[i] = 4.0 * g_x_xy_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_z_xy[i] = 4.0 * g_x_xy_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_z_xz[i] = 4.0 * g_x_xy_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_z_yy[i] = 4.0 * g_x_xy_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_z_yz[i] = 4.0 * g_x_xy_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_z_zz[i] = 4.0 * g_x_xy_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_y_0_0_xz_x_xx, g_x_0_y_0_0_xz_x_xy, g_x_0_y_0_0_xz_x_xz, g_x_0_y_0_0_xz_x_yy, g_x_0_y_0_0_xz_x_yz, g_x_0_y_0_0_xz_x_zz, g_x_xz_xy_xx, g_x_xz_xy_xy, g_x_xz_xy_xz, g_x_xz_xy_yy, g_x_xz_xy_yz, g_x_xz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xz_x_xx[i] = 4.0 * g_x_xz_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_x_xy[i] = 4.0 * g_x_xz_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_x_xz[i] = 4.0 * g_x_xz_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_x_yy[i] = 4.0 * g_x_xz_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_x_yz[i] = 4.0 * g_x_xz_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_x_zz[i] = 4.0 * g_x_xz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_y_0_0_xz_y_xx, g_x_0_y_0_0_xz_y_xy, g_x_0_y_0_0_xz_y_xz, g_x_0_y_0_0_xz_y_yy, g_x_0_y_0_0_xz_y_yz, g_x_0_y_0_0_xz_y_zz, g_x_xz_0_xx, g_x_xz_0_xy, g_x_xz_0_xz, g_x_xz_0_yy, g_x_xz_0_yz, g_x_xz_0_zz, g_x_xz_yy_xx, g_x_xz_yy_xy, g_x_xz_yy_xz, g_x_xz_yy_yy, g_x_xz_yy_yz, g_x_xz_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xz_y_xx[i] = -2.0 * g_x_xz_0_xx[i] * a_exp + 4.0 * g_x_xz_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_y_xy[i] = -2.0 * g_x_xz_0_xy[i] * a_exp + 4.0 * g_x_xz_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_y_xz[i] = -2.0 * g_x_xz_0_xz[i] * a_exp + 4.0 * g_x_xz_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_y_yy[i] = -2.0 * g_x_xz_0_yy[i] * a_exp + 4.0 * g_x_xz_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_y_yz[i] = -2.0 * g_x_xz_0_yz[i] * a_exp + 4.0 * g_x_xz_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_y_zz[i] = -2.0 * g_x_xz_0_zz[i] * a_exp + 4.0 * g_x_xz_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_y_0_0_xz_z_xx, g_x_0_y_0_0_xz_z_xy, g_x_0_y_0_0_xz_z_xz, g_x_0_y_0_0_xz_z_yy, g_x_0_y_0_0_xz_z_yz, g_x_0_y_0_0_xz_z_zz, g_x_xz_yz_xx, g_x_xz_yz_xy, g_x_xz_yz_xz, g_x_xz_yz_yy, g_x_xz_yz_yz, g_x_xz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xz_z_xx[i] = 4.0 * g_x_xz_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_z_xy[i] = 4.0 * g_x_xz_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_z_xz[i] = 4.0 * g_x_xz_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_z_yy[i] = 4.0 * g_x_xz_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_z_yz[i] = 4.0 * g_x_xz_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_z_zz[i] = 4.0 * g_x_xz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_0_y_0_0_yy_x_xx, g_x_0_y_0_0_yy_x_xy, g_x_0_y_0_0_yy_x_xz, g_x_0_y_0_0_yy_x_yy, g_x_0_y_0_0_yy_x_yz, g_x_0_y_0_0_yy_x_zz, g_x_yy_xy_xx, g_x_yy_xy_xy, g_x_yy_xy_xz, g_x_yy_xy_yy, g_x_yy_xy_yz, g_x_yy_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yy_x_xx[i] = 4.0 * g_x_yy_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_x_xy[i] = 4.0 * g_x_yy_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_x_xz[i] = 4.0 * g_x_yy_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_x_yy[i] = 4.0 * g_x_yy_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_x_yz[i] = 4.0 * g_x_yy_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_x_zz[i] = 4.0 * g_x_yy_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_0_y_0_0_yy_y_xx, g_x_0_y_0_0_yy_y_xy, g_x_0_y_0_0_yy_y_xz, g_x_0_y_0_0_yy_y_yy, g_x_0_y_0_0_yy_y_yz, g_x_0_y_0_0_yy_y_zz, g_x_yy_0_xx, g_x_yy_0_xy, g_x_yy_0_xz, g_x_yy_0_yy, g_x_yy_0_yz, g_x_yy_0_zz, g_x_yy_yy_xx, g_x_yy_yy_xy, g_x_yy_yy_xz, g_x_yy_yy_yy, g_x_yy_yy_yz, g_x_yy_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yy_y_xx[i] = -2.0 * g_x_yy_0_xx[i] * a_exp + 4.0 * g_x_yy_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_y_xy[i] = -2.0 * g_x_yy_0_xy[i] * a_exp + 4.0 * g_x_yy_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_y_xz[i] = -2.0 * g_x_yy_0_xz[i] * a_exp + 4.0 * g_x_yy_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_y_yy[i] = -2.0 * g_x_yy_0_yy[i] * a_exp + 4.0 * g_x_yy_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_y_yz[i] = -2.0 * g_x_yy_0_yz[i] * a_exp + 4.0 * g_x_yy_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_y_zz[i] = -2.0 * g_x_yy_0_zz[i] * a_exp + 4.0 * g_x_yy_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_0_y_0_0_yy_z_xx, g_x_0_y_0_0_yy_z_xy, g_x_0_y_0_0_yy_z_xz, g_x_0_y_0_0_yy_z_yy, g_x_0_y_0_0_yy_z_yz, g_x_0_y_0_0_yy_z_zz, g_x_yy_yz_xx, g_x_yy_yz_xy, g_x_yy_yz_xz, g_x_yy_yz_yy, g_x_yy_yz_yz, g_x_yy_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yy_z_xx[i] = 4.0 * g_x_yy_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_z_xy[i] = 4.0 * g_x_yy_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_z_xz[i] = 4.0 * g_x_yy_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_z_yy[i] = 4.0 * g_x_yy_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_z_yz[i] = 4.0 * g_x_yy_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_z_zz[i] = 4.0 * g_x_yy_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_0_y_0_0_yz_x_xx, g_x_0_y_0_0_yz_x_xy, g_x_0_y_0_0_yz_x_xz, g_x_0_y_0_0_yz_x_yy, g_x_0_y_0_0_yz_x_yz, g_x_0_y_0_0_yz_x_zz, g_x_yz_xy_xx, g_x_yz_xy_xy, g_x_yz_xy_xz, g_x_yz_xy_yy, g_x_yz_xy_yz, g_x_yz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yz_x_xx[i] = 4.0 * g_x_yz_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_x_xy[i] = 4.0 * g_x_yz_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_x_xz[i] = 4.0 * g_x_yz_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_x_yy[i] = 4.0 * g_x_yz_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_x_yz[i] = 4.0 * g_x_yz_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_x_zz[i] = 4.0 * g_x_yz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_0_y_0_0_yz_y_xx, g_x_0_y_0_0_yz_y_xy, g_x_0_y_0_0_yz_y_xz, g_x_0_y_0_0_yz_y_yy, g_x_0_y_0_0_yz_y_yz, g_x_0_y_0_0_yz_y_zz, g_x_yz_0_xx, g_x_yz_0_xy, g_x_yz_0_xz, g_x_yz_0_yy, g_x_yz_0_yz, g_x_yz_0_zz, g_x_yz_yy_xx, g_x_yz_yy_xy, g_x_yz_yy_xz, g_x_yz_yy_yy, g_x_yz_yy_yz, g_x_yz_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yz_y_xx[i] = -2.0 * g_x_yz_0_xx[i] * a_exp + 4.0 * g_x_yz_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_y_xy[i] = -2.0 * g_x_yz_0_xy[i] * a_exp + 4.0 * g_x_yz_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_y_xz[i] = -2.0 * g_x_yz_0_xz[i] * a_exp + 4.0 * g_x_yz_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_y_yy[i] = -2.0 * g_x_yz_0_yy[i] * a_exp + 4.0 * g_x_yz_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_y_yz[i] = -2.0 * g_x_yz_0_yz[i] * a_exp + 4.0 * g_x_yz_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_y_zz[i] = -2.0 * g_x_yz_0_zz[i] * a_exp + 4.0 * g_x_yz_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_0_y_0_0_yz_z_xx, g_x_0_y_0_0_yz_z_xy, g_x_0_y_0_0_yz_z_xz, g_x_0_y_0_0_yz_z_yy, g_x_0_y_0_0_yz_z_yz, g_x_0_y_0_0_yz_z_zz, g_x_yz_yz_xx, g_x_yz_yz_xy, g_x_yz_yz_xz, g_x_yz_yz_yy, g_x_yz_yz_yz, g_x_yz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yz_z_xx[i] = 4.0 * g_x_yz_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_z_xy[i] = 4.0 * g_x_yz_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_z_xz[i] = 4.0 * g_x_yz_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_z_yy[i] = 4.0 * g_x_yz_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_z_yz[i] = 4.0 * g_x_yz_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_z_zz[i] = 4.0 * g_x_yz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_0_y_0_0_zz_x_xx, g_x_0_y_0_0_zz_x_xy, g_x_0_y_0_0_zz_x_xz, g_x_0_y_0_0_zz_x_yy, g_x_0_y_0_0_zz_x_yz, g_x_0_y_0_0_zz_x_zz, g_x_zz_xy_xx, g_x_zz_xy_xy, g_x_zz_xy_xz, g_x_zz_xy_yy, g_x_zz_xy_yz, g_x_zz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_zz_x_xx[i] = 4.0 * g_x_zz_xy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_x_xy[i] = 4.0 * g_x_zz_xy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_x_xz[i] = 4.0 * g_x_zz_xy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_x_yy[i] = 4.0 * g_x_zz_xy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_x_yz[i] = 4.0 * g_x_zz_xy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_x_zz[i] = 4.0 * g_x_zz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_0_y_0_0_zz_y_xx, g_x_0_y_0_0_zz_y_xy, g_x_0_y_0_0_zz_y_xz, g_x_0_y_0_0_zz_y_yy, g_x_0_y_0_0_zz_y_yz, g_x_0_y_0_0_zz_y_zz, g_x_zz_0_xx, g_x_zz_0_xy, g_x_zz_0_xz, g_x_zz_0_yy, g_x_zz_0_yz, g_x_zz_0_zz, g_x_zz_yy_xx, g_x_zz_yy_xy, g_x_zz_yy_xz, g_x_zz_yy_yy, g_x_zz_yy_yz, g_x_zz_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_zz_y_xx[i] = -2.0 * g_x_zz_0_xx[i] * a_exp + 4.0 * g_x_zz_yy_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_y_xy[i] = -2.0 * g_x_zz_0_xy[i] * a_exp + 4.0 * g_x_zz_yy_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_y_xz[i] = -2.0 * g_x_zz_0_xz[i] * a_exp + 4.0 * g_x_zz_yy_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_y_yy[i] = -2.0 * g_x_zz_0_yy[i] * a_exp + 4.0 * g_x_zz_yy_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_y_yz[i] = -2.0 * g_x_zz_0_yz[i] * a_exp + 4.0 * g_x_zz_yy_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_y_zz[i] = -2.0 * g_x_zz_0_zz[i] * a_exp + 4.0 * g_x_zz_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_0_y_0_0_zz_z_xx, g_x_0_y_0_0_zz_z_xy, g_x_0_y_0_0_zz_z_xz, g_x_0_y_0_0_zz_z_yy, g_x_0_y_0_0_zz_z_yz, g_x_0_y_0_0_zz_z_zz, g_x_zz_yz_xx, g_x_zz_yz_xy, g_x_zz_yz_xz, g_x_zz_yz_yy, g_x_zz_yz_yz, g_x_zz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_zz_z_xx[i] = 4.0 * g_x_zz_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_z_xy[i] = 4.0 * g_x_zz_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_z_xz[i] = 4.0 * g_x_zz_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_z_yy[i] = 4.0 * g_x_zz_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_z_yz[i] = 4.0 * g_x_zz_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_z_zz[i] = 4.0 * g_x_zz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_0_z_0_0_xx_x_xx, g_x_0_z_0_0_xx_x_xy, g_x_0_z_0_0_xx_x_xz, g_x_0_z_0_0_xx_x_yy, g_x_0_z_0_0_xx_x_yz, g_x_0_z_0_0_xx_x_zz, g_x_xx_xz_xx, g_x_xx_xz_xy, g_x_xx_xz_xz, g_x_xx_xz_yy, g_x_xx_xz_yz, g_x_xx_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xx_x_xx[i] = 4.0 * g_x_xx_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_x_xy[i] = 4.0 * g_x_xx_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_x_xz[i] = 4.0 * g_x_xx_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_x_yy[i] = 4.0 * g_x_xx_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_x_yz[i] = 4.0 * g_x_xx_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_x_zz[i] = 4.0 * g_x_xx_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_0_z_0_0_xx_y_xx, g_x_0_z_0_0_xx_y_xy, g_x_0_z_0_0_xx_y_xz, g_x_0_z_0_0_xx_y_yy, g_x_0_z_0_0_xx_y_yz, g_x_0_z_0_0_xx_y_zz, g_x_xx_yz_xx, g_x_xx_yz_xy, g_x_xx_yz_xz, g_x_xx_yz_yy, g_x_xx_yz_yz, g_x_xx_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xx_y_xx[i] = 4.0 * g_x_xx_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_y_xy[i] = 4.0 * g_x_xx_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_y_xz[i] = 4.0 * g_x_xx_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_y_yy[i] = 4.0 * g_x_xx_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_y_yz[i] = 4.0 * g_x_xx_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_y_zz[i] = 4.0 * g_x_xx_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_0_z_0_0_xx_z_xx, g_x_0_z_0_0_xx_z_xy, g_x_0_z_0_0_xx_z_xz, g_x_0_z_0_0_xx_z_yy, g_x_0_z_0_0_xx_z_yz, g_x_0_z_0_0_xx_z_zz, g_x_xx_0_xx, g_x_xx_0_xy, g_x_xx_0_xz, g_x_xx_0_yy, g_x_xx_0_yz, g_x_xx_0_zz, g_x_xx_zz_xx, g_x_xx_zz_xy, g_x_xx_zz_xz, g_x_xx_zz_yy, g_x_xx_zz_yz, g_x_xx_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xx_z_xx[i] = -2.0 * g_x_xx_0_xx[i] * a_exp + 4.0 * g_x_xx_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_z_xy[i] = -2.0 * g_x_xx_0_xy[i] * a_exp + 4.0 * g_x_xx_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_z_xz[i] = -2.0 * g_x_xx_0_xz[i] * a_exp + 4.0 * g_x_xx_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_z_yy[i] = -2.0 * g_x_xx_0_yy[i] * a_exp + 4.0 * g_x_xx_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_z_yz[i] = -2.0 * g_x_xx_0_yz[i] * a_exp + 4.0 * g_x_xx_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_z_zz[i] = -2.0 * g_x_xx_0_zz[i] * a_exp + 4.0 * g_x_xx_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_0_z_0_0_xy_x_xx, g_x_0_z_0_0_xy_x_xy, g_x_0_z_0_0_xy_x_xz, g_x_0_z_0_0_xy_x_yy, g_x_0_z_0_0_xy_x_yz, g_x_0_z_0_0_xy_x_zz, g_x_xy_xz_xx, g_x_xy_xz_xy, g_x_xy_xz_xz, g_x_xy_xz_yy, g_x_xy_xz_yz, g_x_xy_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xy_x_xx[i] = 4.0 * g_x_xy_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_x_xy[i] = 4.0 * g_x_xy_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_x_xz[i] = 4.0 * g_x_xy_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_x_yy[i] = 4.0 * g_x_xy_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_x_yz[i] = 4.0 * g_x_xy_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_x_zz[i] = 4.0 * g_x_xy_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_0_z_0_0_xy_y_xx, g_x_0_z_0_0_xy_y_xy, g_x_0_z_0_0_xy_y_xz, g_x_0_z_0_0_xy_y_yy, g_x_0_z_0_0_xy_y_yz, g_x_0_z_0_0_xy_y_zz, g_x_xy_yz_xx, g_x_xy_yz_xy, g_x_xy_yz_xz, g_x_xy_yz_yy, g_x_xy_yz_yz, g_x_xy_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xy_y_xx[i] = 4.0 * g_x_xy_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_y_xy[i] = 4.0 * g_x_xy_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_y_xz[i] = 4.0 * g_x_xy_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_y_yy[i] = 4.0 * g_x_xy_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_y_yz[i] = 4.0 * g_x_xy_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_y_zz[i] = 4.0 * g_x_xy_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_0_z_0_0_xy_z_xx, g_x_0_z_0_0_xy_z_xy, g_x_0_z_0_0_xy_z_xz, g_x_0_z_0_0_xy_z_yy, g_x_0_z_0_0_xy_z_yz, g_x_0_z_0_0_xy_z_zz, g_x_xy_0_xx, g_x_xy_0_xy, g_x_xy_0_xz, g_x_xy_0_yy, g_x_xy_0_yz, g_x_xy_0_zz, g_x_xy_zz_xx, g_x_xy_zz_xy, g_x_xy_zz_xz, g_x_xy_zz_yy, g_x_xy_zz_yz, g_x_xy_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xy_z_xx[i] = -2.0 * g_x_xy_0_xx[i] * a_exp + 4.0 * g_x_xy_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_z_xy[i] = -2.0 * g_x_xy_0_xy[i] * a_exp + 4.0 * g_x_xy_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_z_xz[i] = -2.0 * g_x_xy_0_xz[i] * a_exp + 4.0 * g_x_xy_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_z_yy[i] = -2.0 * g_x_xy_0_yy[i] * a_exp + 4.0 * g_x_xy_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_z_yz[i] = -2.0 * g_x_xy_0_yz[i] * a_exp + 4.0 * g_x_xy_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_z_zz[i] = -2.0 * g_x_xy_0_zz[i] * a_exp + 4.0 * g_x_xy_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_0_z_0_0_xz_x_xx, g_x_0_z_0_0_xz_x_xy, g_x_0_z_0_0_xz_x_xz, g_x_0_z_0_0_xz_x_yy, g_x_0_z_0_0_xz_x_yz, g_x_0_z_0_0_xz_x_zz, g_x_xz_xz_xx, g_x_xz_xz_xy, g_x_xz_xz_xz, g_x_xz_xz_yy, g_x_xz_xz_yz, g_x_xz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xz_x_xx[i] = 4.0 * g_x_xz_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_x_xy[i] = 4.0 * g_x_xz_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_x_xz[i] = 4.0 * g_x_xz_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_x_yy[i] = 4.0 * g_x_xz_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_x_yz[i] = 4.0 * g_x_xz_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_x_zz[i] = 4.0 * g_x_xz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_0_z_0_0_xz_y_xx, g_x_0_z_0_0_xz_y_xy, g_x_0_z_0_0_xz_y_xz, g_x_0_z_0_0_xz_y_yy, g_x_0_z_0_0_xz_y_yz, g_x_0_z_0_0_xz_y_zz, g_x_xz_yz_xx, g_x_xz_yz_xy, g_x_xz_yz_xz, g_x_xz_yz_yy, g_x_xz_yz_yz, g_x_xz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xz_y_xx[i] = 4.0 * g_x_xz_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_y_xy[i] = 4.0 * g_x_xz_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_y_xz[i] = 4.0 * g_x_xz_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_y_yy[i] = 4.0 * g_x_xz_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_y_yz[i] = 4.0 * g_x_xz_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_y_zz[i] = 4.0 * g_x_xz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_0_z_0_0_xz_z_xx, g_x_0_z_0_0_xz_z_xy, g_x_0_z_0_0_xz_z_xz, g_x_0_z_0_0_xz_z_yy, g_x_0_z_0_0_xz_z_yz, g_x_0_z_0_0_xz_z_zz, g_x_xz_0_xx, g_x_xz_0_xy, g_x_xz_0_xz, g_x_xz_0_yy, g_x_xz_0_yz, g_x_xz_0_zz, g_x_xz_zz_xx, g_x_xz_zz_xy, g_x_xz_zz_xz, g_x_xz_zz_yy, g_x_xz_zz_yz, g_x_xz_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xz_z_xx[i] = -2.0 * g_x_xz_0_xx[i] * a_exp + 4.0 * g_x_xz_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_z_xy[i] = -2.0 * g_x_xz_0_xy[i] * a_exp + 4.0 * g_x_xz_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_z_xz[i] = -2.0 * g_x_xz_0_xz[i] * a_exp + 4.0 * g_x_xz_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_z_yy[i] = -2.0 * g_x_xz_0_yy[i] * a_exp + 4.0 * g_x_xz_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_z_yz[i] = -2.0 * g_x_xz_0_yz[i] * a_exp + 4.0 * g_x_xz_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_z_zz[i] = -2.0 * g_x_xz_0_zz[i] * a_exp + 4.0 * g_x_xz_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_0_z_0_0_yy_x_xx, g_x_0_z_0_0_yy_x_xy, g_x_0_z_0_0_yy_x_xz, g_x_0_z_0_0_yy_x_yy, g_x_0_z_0_0_yy_x_yz, g_x_0_z_0_0_yy_x_zz, g_x_yy_xz_xx, g_x_yy_xz_xy, g_x_yy_xz_xz, g_x_yy_xz_yy, g_x_yy_xz_yz, g_x_yy_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yy_x_xx[i] = 4.0 * g_x_yy_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_x_xy[i] = 4.0 * g_x_yy_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_x_xz[i] = 4.0 * g_x_yy_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_x_yy[i] = 4.0 * g_x_yy_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_x_yz[i] = 4.0 * g_x_yy_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_x_zz[i] = 4.0 * g_x_yy_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_0_z_0_0_yy_y_xx, g_x_0_z_0_0_yy_y_xy, g_x_0_z_0_0_yy_y_xz, g_x_0_z_0_0_yy_y_yy, g_x_0_z_0_0_yy_y_yz, g_x_0_z_0_0_yy_y_zz, g_x_yy_yz_xx, g_x_yy_yz_xy, g_x_yy_yz_xz, g_x_yy_yz_yy, g_x_yy_yz_yz, g_x_yy_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yy_y_xx[i] = 4.0 * g_x_yy_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_y_xy[i] = 4.0 * g_x_yy_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_y_xz[i] = 4.0 * g_x_yy_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_y_yy[i] = 4.0 * g_x_yy_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_y_yz[i] = 4.0 * g_x_yy_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_y_zz[i] = 4.0 * g_x_yy_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_0_z_0_0_yy_z_xx, g_x_0_z_0_0_yy_z_xy, g_x_0_z_0_0_yy_z_xz, g_x_0_z_0_0_yy_z_yy, g_x_0_z_0_0_yy_z_yz, g_x_0_z_0_0_yy_z_zz, g_x_yy_0_xx, g_x_yy_0_xy, g_x_yy_0_xz, g_x_yy_0_yy, g_x_yy_0_yz, g_x_yy_0_zz, g_x_yy_zz_xx, g_x_yy_zz_xy, g_x_yy_zz_xz, g_x_yy_zz_yy, g_x_yy_zz_yz, g_x_yy_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yy_z_xx[i] = -2.0 * g_x_yy_0_xx[i] * a_exp + 4.0 * g_x_yy_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_z_xy[i] = -2.0 * g_x_yy_0_xy[i] * a_exp + 4.0 * g_x_yy_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_z_xz[i] = -2.0 * g_x_yy_0_xz[i] * a_exp + 4.0 * g_x_yy_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_z_yy[i] = -2.0 * g_x_yy_0_yy[i] * a_exp + 4.0 * g_x_yy_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_z_yz[i] = -2.0 * g_x_yy_0_yz[i] * a_exp + 4.0 * g_x_yy_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_z_zz[i] = -2.0 * g_x_yy_0_zz[i] * a_exp + 4.0 * g_x_yy_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_0_z_0_0_yz_x_xx, g_x_0_z_0_0_yz_x_xy, g_x_0_z_0_0_yz_x_xz, g_x_0_z_0_0_yz_x_yy, g_x_0_z_0_0_yz_x_yz, g_x_0_z_0_0_yz_x_zz, g_x_yz_xz_xx, g_x_yz_xz_xy, g_x_yz_xz_xz, g_x_yz_xz_yy, g_x_yz_xz_yz, g_x_yz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yz_x_xx[i] = 4.0 * g_x_yz_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_x_xy[i] = 4.0 * g_x_yz_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_x_xz[i] = 4.0 * g_x_yz_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_x_yy[i] = 4.0 * g_x_yz_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_x_yz[i] = 4.0 * g_x_yz_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_x_zz[i] = 4.0 * g_x_yz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_0_z_0_0_yz_y_xx, g_x_0_z_0_0_yz_y_xy, g_x_0_z_0_0_yz_y_xz, g_x_0_z_0_0_yz_y_yy, g_x_0_z_0_0_yz_y_yz, g_x_0_z_0_0_yz_y_zz, g_x_yz_yz_xx, g_x_yz_yz_xy, g_x_yz_yz_xz, g_x_yz_yz_yy, g_x_yz_yz_yz, g_x_yz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yz_y_xx[i] = 4.0 * g_x_yz_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_y_xy[i] = 4.0 * g_x_yz_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_y_xz[i] = 4.0 * g_x_yz_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_y_yy[i] = 4.0 * g_x_yz_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_y_yz[i] = 4.0 * g_x_yz_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_y_zz[i] = 4.0 * g_x_yz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_0_z_0_0_yz_z_xx, g_x_0_z_0_0_yz_z_xy, g_x_0_z_0_0_yz_z_xz, g_x_0_z_0_0_yz_z_yy, g_x_0_z_0_0_yz_z_yz, g_x_0_z_0_0_yz_z_zz, g_x_yz_0_xx, g_x_yz_0_xy, g_x_yz_0_xz, g_x_yz_0_yy, g_x_yz_0_yz, g_x_yz_0_zz, g_x_yz_zz_xx, g_x_yz_zz_xy, g_x_yz_zz_xz, g_x_yz_zz_yy, g_x_yz_zz_yz, g_x_yz_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yz_z_xx[i] = -2.0 * g_x_yz_0_xx[i] * a_exp + 4.0 * g_x_yz_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_z_xy[i] = -2.0 * g_x_yz_0_xy[i] * a_exp + 4.0 * g_x_yz_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_z_xz[i] = -2.0 * g_x_yz_0_xz[i] * a_exp + 4.0 * g_x_yz_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_z_yy[i] = -2.0 * g_x_yz_0_yy[i] * a_exp + 4.0 * g_x_yz_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_z_yz[i] = -2.0 * g_x_yz_0_yz[i] * a_exp + 4.0 * g_x_yz_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_z_zz[i] = -2.0 * g_x_yz_0_zz[i] * a_exp + 4.0 * g_x_yz_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_0_z_0_0_zz_x_xx, g_x_0_z_0_0_zz_x_xy, g_x_0_z_0_0_zz_x_xz, g_x_0_z_0_0_zz_x_yy, g_x_0_z_0_0_zz_x_yz, g_x_0_z_0_0_zz_x_zz, g_x_zz_xz_xx, g_x_zz_xz_xy, g_x_zz_xz_xz, g_x_zz_xz_yy, g_x_zz_xz_yz, g_x_zz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_zz_x_xx[i] = 4.0 * g_x_zz_xz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_x_xy[i] = 4.0 * g_x_zz_xz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_x_xz[i] = 4.0 * g_x_zz_xz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_x_yy[i] = 4.0 * g_x_zz_xz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_x_yz[i] = 4.0 * g_x_zz_xz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_x_zz[i] = 4.0 * g_x_zz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_0_z_0_0_zz_y_xx, g_x_0_z_0_0_zz_y_xy, g_x_0_z_0_0_zz_y_xz, g_x_0_z_0_0_zz_y_yy, g_x_0_z_0_0_zz_y_yz, g_x_0_z_0_0_zz_y_zz, g_x_zz_yz_xx, g_x_zz_yz_xy, g_x_zz_yz_xz, g_x_zz_yz_yy, g_x_zz_yz_yz, g_x_zz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_zz_y_xx[i] = 4.0 * g_x_zz_yz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_y_xy[i] = 4.0 * g_x_zz_yz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_y_xz[i] = 4.0 * g_x_zz_yz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_y_yy[i] = 4.0 * g_x_zz_yz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_y_yz[i] = 4.0 * g_x_zz_yz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_y_zz[i] = 4.0 * g_x_zz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_0_z_0_0_zz_z_xx, g_x_0_z_0_0_zz_z_xy, g_x_0_z_0_0_zz_z_xz, g_x_0_z_0_0_zz_z_yy, g_x_0_z_0_0_zz_z_yz, g_x_0_z_0_0_zz_z_zz, g_x_zz_0_xx, g_x_zz_0_xy, g_x_zz_0_xz, g_x_zz_0_yy, g_x_zz_0_yz, g_x_zz_0_zz, g_x_zz_zz_xx, g_x_zz_zz_xy, g_x_zz_zz_xz, g_x_zz_zz_yy, g_x_zz_zz_yz, g_x_zz_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_zz_z_xx[i] = -2.0 * g_x_zz_0_xx[i] * a_exp + 4.0 * g_x_zz_zz_xx[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_z_xy[i] = -2.0 * g_x_zz_0_xy[i] * a_exp + 4.0 * g_x_zz_zz_xy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_z_xz[i] = -2.0 * g_x_zz_0_xz[i] * a_exp + 4.0 * g_x_zz_zz_xz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_z_yy[i] = -2.0 * g_x_zz_0_yy[i] * a_exp + 4.0 * g_x_zz_zz_yy[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_z_yz[i] = -2.0 * g_x_zz_0_yz[i] * a_exp + 4.0 * g_x_zz_zz_yz[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_z_zz[i] = -2.0 * g_x_zz_0_zz[i] * a_exp + 4.0 * g_x_zz_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_y_0_x_0_0_xx_x_xx, g_y_0_x_0_0_xx_x_xy, g_y_0_x_0_0_xx_x_xz, g_y_0_x_0_0_xx_x_yy, g_y_0_x_0_0_xx_x_yz, g_y_0_x_0_0_xx_x_zz, g_y_xx_0_xx, g_y_xx_0_xy, g_y_xx_0_xz, g_y_xx_0_yy, g_y_xx_0_yz, g_y_xx_0_zz, g_y_xx_xx_xx, g_y_xx_xx_xy, g_y_xx_xx_xz, g_y_xx_xx_yy, g_y_xx_xx_yz, g_y_xx_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xx_x_xx[i] = -2.0 * g_y_xx_0_xx[i] * a_exp + 4.0 * g_y_xx_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_x_xy[i] = -2.0 * g_y_xx_0_xy[i] * a_exp + 4.0 * g_y_xx_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_x_xz[i] = -2.0 * g_y_xx_0_xz[i] * a_exp + 4.0 * g_y_xx_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_x_yy[i] = -2.0 * g_y_xx_0_yy[i] * a_exp + 4.0 * g_y_xx_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_x_yz[i] = -2.0 * g_y_xx_0_yz[i] * a_exp + 4.0 * g_y_xx_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_x_zz[i] = -2.0 * g_y_xx_0_zz[i] * a_exp + 4.0 * g_y_xx_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_y_0_x_0_0_xx_y_xx, g_y_0_x_0_0_xx_y_xy, g_y_0_x_0_0_xx_y_xz, g_y_0_x_0_0_xx_y_yy, g_y_0_x_0_0_xx_y_yz, g_y_0_x_0_0_xx_y_zz, g_y_xx_xy_xx, g_y_xx_xy_xy, g_y_xx_xy_xz, g_y_xx_xy_yy, g_y_xx_xy_yz, g_y_xx_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xx_y_xx[i] = 4.0 * g_y_xx_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_y_xy[i] = 4.0 * g_y_xx_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_y_xz[i] = 4.0 * g_y_xx_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_y_yy[i] = 4.0 * g_y_xx_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_y_yz[i] = 4.0 * g_y_xx_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_y_zz[i] = 4.0 * g_y_xx_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_y_0_x_0_0_xx_z_xx, g_y_0_x_0_0_xx_z_xy, g_y_0_x_0_0_xx_z_xz, g_y_0_x_0_0_xx_z_yy, g_y_0_x_0_0_xx_z_yz, g_y_0_x_0_0_xx_z_zz, g_y_xx_xz_xx, g_y_xx_xz_xy, g_y_xx_xz_xz, g_y_xx_xz_yy, g_y_xx_xz_yz, g_y_xx_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xx_z_xx[i] = 4.0 * g_y_xx_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_z_xy[i] = 4.0 * g_y_xx_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_z_xz[i] = 4.0 * g_y_xx_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_z_yy[i] = 4.0 * g_y_xx_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_z_yz[i] = 4.0 * g_y_xx_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_z_zz[i] = 4.0 * g_y_xx_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_y_0_x_0_0_xy_x_xx, g_y_0_x_0_0_xy_x_xy, g_y_0_x_0_0_xy_x_xz, g_y_0_x_0_0_xy_x_yy, g_y_0_x_0_0_xy_x_yz, g_y_0_x_0_0_xy_x_zz, g_y_xy_0_xx, g_y_xy_0_xy, g_y_xy_0_xz, g_y_xy_0_yy, g_y_xy_0_yz, g_y_xy_0_zz, g_y_xy_xx_xx, g_y_xy_xx_xy, g_y_xy_xx_xz, g_y_xy_xx_yy, g_y_xy_xx_yz, g_y_xy_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xy_x_xx[i] = -2.0 * g_y_xy_0_xx[i] * a_exp + 4.0 * g_y_xy_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_x_xy[i] = -2.0 * g_y_xy_0_xy[i] * a_exp + 4.0 * g_y_xy_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_x_xz[i] = -2.0 * g_y_xy_0_xz[i] * a_exp + 4.0 * g_y_xy_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_x_yy[i] = -2.0 * g_y_xy_0_yy[i] * a_exp + 4.0 * g_y_xy_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_x_yz[i] = -2.0 * g_y_xy_0_yz[i] * a_exp + 4.0 * g_y_xy_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_x_zz[i] = -2.0 * g_y_xy_0_zz[i] * a_exp + 4.0 * g_y_xy_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_y_0_x_0_0_xy_y_xx, g_y_0_x_0_0_xy_y_xy, g_y_0_x_0_0_xy_y_xz, g_y_0_x_0_0_xy_y_yy, g_y_0_x_0_0_xy_y_yz, g_y_0_x_0_0_xy_y_zz, g_y_xy_xy_xx, g_y_xy_xy_xy, g_y_xy_xy_xz, g_y_xy_xy_yy, g_y_xy_xy_yz, g_y_xy_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xy_y_xx[i] = 4.0 * g_y_xy_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_y_xy[i] = 4.0 * g_y_xy_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_y_xz[i] = 4.0 * g_y_xy_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_y_yy[i] = 4.0 * g_y_xy_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_y_yz[i] = 4.0 * g_y_xy_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_y_zz[i] = 4.0 * g_y_xy_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_y_0_x_0_0_xy_z_xx, g_y_0_x_0_0_xy_z_xy, g_y_0_x_0_0_xy_z_xz, g_y_0_x_0_0_xy_z_yy, g_y_0_x_0_0_xy_z_yz, g_y_0_x_0_0_xy_z_zz, g_y_xy_xz_xx, g_y_xy_xz_xy, g_y_xy_xz_xz, g_y_xy_xz_yy, g_y_xy_xz_yz, g_y_xy_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xy_z_xx[i] = 4.0 * g_y_xy_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_z_xy[i] = 4.0 * g_y_xy_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_z_xz[i] = 4.0 * g_y_xy_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_z_yy[i] = 4.0 * g_y_xy_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_z_yz[i] = 4.0 * g_y_xy_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_z_zz[i] = 4.0 * g_y_xy_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_y_0_x_0_0_xz_x_xx, g_y_0_x_0_0_xz_x_xy, g_y_0_x_0_0_xz_x_xz, g_y_0_x_0_0_xz_x_yy, g_y_0_x_0_0_xz_x_yz, g_y_0_x_0_0_xz_x_zz, g_y_xz_0_xx, g_y_xz_0_xy, g_y_xz_0_xz, g_y_xz_0_yy, g_y_xz_0_yz, g_y_xz_0_zz, g_y_xz_xx_xx, g_y_xz_xx_xy, g_y_xz_xx_xz, g_y_xz_xx_yy, g_y_xz_xx_yz, g_y_xz_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xz_x_xx[i] = -2.0 * g_y_xz_0_xx[i] * a_exp + 4.0 * g_y_xz_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_x_xy[i] = -2.0 * g_y_xz_0_xy[i] * a_exp + 4.0 * g_y_xz_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_x_xz[i] = -2.0 * g_y_xz_0_xz[i] * a_exp + 4.0 * g_y_xz_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_x_yy[i] = -2.0 * g_y_xz_0_yy[i] * a_exp + 4.0 * g_y_xz_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_x_yz[i] = -2.0 * g_y_xz_0_yz[i] * a_exp + 4.0 * g_y_xz_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_x_zz[i] = -2.0 * g_y_xz_0_zz[i] * a_exp + 4.0 * g_y_xz_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_y_0_x_0_0_xz_y_xx, g_y_0_x_0_0_xz_y_xy, g_y_0_x_0_0_xz_y_xz, g_y_0_x_0_0_xz_y_yy, g_y_0_x_0_0_xz_y_yz, g_y_0_x_0_0_xz_y_zz, g_y_xz_xy_xx, g_y_xz_xy_xy, g_y_xz_xy_xz, g_y_xz_xy_yy, g_y_xz_xy_yz, g_y_xz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xz_y_xx[i] = 4.0 * g_y_xz_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_y_xy[i] = 4.0 * g_y_xz_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_y_xz[i] = 4.0 * g_y_xz_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_y_yy[i] = 4.0 * g_y_xz_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_y_yz[i] = 4.0 * g_y_xz_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_y_zz[i] = 4.0 * g_y_xz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_y_0_x_0_0_xz_z_xx, g_y_0_x_0_0_xz_z_xy, g_y_0_x_0_0_xz_z_xz, g_y_0_x_0_0_xz_z_yy, g_y_0_x_0_0_xz_z_yz, g_y_0_x_0_0_xz_z_zz, g_y_xz_xz_xx, g_y_xz_xz_xy, g_y_xz_xz_xz, g_y_xz_xz_yy, g_y_xz_xz_yz, g_y_xz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xz_z_xx[i] = 4.0 * g_y_xz_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_z_xy[i] = 4.0 * g_y_xz_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_z_xz[i] = 4.0 * g_y_xz_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_z_yy[i] = 4.0 * g_y_xz_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_z_yz[i] = 4.0 * g_y_xz_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_z_zz[i] = 4.0 * g_y_xz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_y_0_x_0_0_yy_x_xx, g_y_0_x_0_0_yy_x_xy, g_y_0_x_0_0_yy_x_xz, g_y_0_x_0_0_yy_x_yy, g_y_0_x_0_0_yy_x_yz, g_y_0_x_0_0_yy_x_zz, g_y_yy_0_xx, g_y_yy_0_xy, g_y_yy_0_xz, g_y_yy_0_yy, g_y_yy_0_yz, g_y_yy_0_zz, g_y_yy_xx_xx, g_y_yy_xx_xy, g_y_yy_xx_xz, g_y_yy_xx_yy, g_y_yy_xx_yz, g_y_yy_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yy_x_xx[i] = -2.0 * g_y_yy_0_xx[i] * a_exp + 4.0 * g_y_yy_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_x_xy[i] = -2.0 * g_y_yy_0_xy[i] * a_exp + 4.0 * g_y_yy_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_x_xz[i] = -2.0 * g_y_yy_0_xz[i] * a_exp + 4.0 * g_y_yy_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_x_yy[i] = -2.0 * g_y_yy_0_yy[i] * a_exp + 4.0 * g_y_yy_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_x_yz[i] = -2.0 * g_y_yy_0_yz[i] * a_exp + 4.0 * g_y_yy_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_x_zz[i] = -2.0 * g_y_yy_0_zz[i] * a_exp + 4.0 * g_y_yy_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_y_0_x_0_0_yy_y_xx, g_y_0_x_0_0_yy_y_xy, g_y_0_x_0_0_yy_y_xz, g_y_0_x_0_0_yy_y_yy, g_y_0_x_0_0_yy_y_yz, g_y_0_x_0_0_yy_y_zz, g_y_yy_xy_xx, g_y_yy_xy_xy, g_y_yy_xy_xz, g_y_yy_xy_yy, g_y_yy_xy_yz, g_y_yy_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yy_y_xx[i] = 4.0 * g_y_yy_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_y_xy[i] = 4.0 * g_y_yy_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_y_xz[i] = 4.0 * g_y_yy_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_y_yy[i] = 4.0 * g_y_yy_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_y_yz[i] = 4.0 * g_y_yy_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_y_zz[i] = 4.0 * g_y_yy_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_y_0_x_0_0_yy_z_xx, g_y_0_x_0_0_yy_z_xy, g_y_0_x_0_0_yy_z_xz, g_y_0_x_0_0_yy_z_yy, g_y_0_x_0_0_yy_z_yz, g_y_0_x_0_0_yy_z_zz, g_y_yy_xz_xx, g_y_yy_xz_xy, g_y_yy_xz_xz, g_y_yy_xz_yy, g_y_yy_xz_yz, g_y_yy_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yy_z_xx[i] = 4.0 * g_y_yy_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_z_xy[i] = 4.0 * g_y_yy_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_z_xz[i] = 4.0 * g_y_yy_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_z_yy[i] = 4.0 * g_y_yy_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_z_yz[i] = 4.0 * g_y_yy_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_z_zz[i] = 4.0 * g_y_yy_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_y_0_x_0_0_yz_x_xx, g_y_0_x_0_0_yz_x_xy, g_y_0_x_0_0_yz_x_xz, g_y_0_x_0_0_yz_x_yy, g_y_0_x_0_0_yz_x_yz, g_y_0_x_0_0_yz_x_zz, g_y_yz_0_xx, g_y_yz_0_xy, g_y_yz_0_xz, g_y_yz_0_yy, g_y_yz_0_yz, g_y_yz_0_zz, g_y_yz_xx_xx, g_y_yz_xx_xy, g_y_yz_xx_xz, g_y_yz_xx_yy, g_y_yz_xx_yz, g_y_yz_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yz_x_xx[i] = -2.0 * g_y_yz_0_xx[i] * a_exp + 4.0 * g_y_yz_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_x_xy[i] = -2.0 * g_y_yz_0_xy[i] * a_exp + 4.0 * g_y_yz_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_x_xz[i] = -2.0 * g_y_yz_0_xz[i] * a_exp + 4.0 * g_y_yz_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_x_yy[i] = -2.0 * g_y_yz_0_yy[i] * a_exp + 4.0 * g_y_yz_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_x_yz[i] = -2.0 * g_y_yz_0_yz[i] * a_exp + 4.0 * g_y_yz_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_x_zz[i] = -2.0 * g_y_yz_0_zz[i] * a_exp + 4.0 * g_y_yz_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_y_0_x_0_0_yz_y_xx, g_y_0_x_0_0_yz_y_xy, g_y_0_x_0_0_yz_y_xz, g_y_0_x_0_0_yz_y_yy, g_y_0_x_0_0_yz_y_yz, g_y_0_x_0_0_yz_y_zz, g_y_yz_xy_xx, g_y_yz_xy_xy, g_y_yz_xy_xz, g_y_yz_xy_yy, g_y_yz_xy_yz, g_y_yz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yz_y_xx[i] = 4.0 * g_y_yz_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_y_xy[i] = 4.0 * g_y_yz_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_y_xz[i] = 4.0 * g_y_yz_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_y_yy[i] = 4.0 * g_y_yz_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_y_yz[i] = 4.0 * g_y_yz_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_y_zz[i] = 4.0 * g_y_yz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_y_0_x_0_0_yz_z_xx, g_y_0_x_0_0_yz_z_xy, g_y_0_x_0_0_yz_z_xz, g_y_0_x_0_0_yz_z_yy, g_y_0_x_0_0_yz_z_yz, g_y_0_x_0_0_yz_z_zz, g_y_yz_xz_xx, g_y_yz_xz_xy, g_y_yz_xz_xz, g_y_yz_xz_yy, g_y_yz_xz_yz, g_y_yz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yz_z_xx[i] = 4.0 * g_y_yz_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_z_xy[i] = 4.0 * g_y_yz_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_z_xz[i] = 4.0 * g_y_yz_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_z_yy[i] = 4.0 * g_y_yz_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_z_yz[i] = 4.0 * g_y_yz_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_z_zz[i] = 4.0 * g_y_yz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_y_0_x_0_0_zz_x_xx, g_y_0_x_0_0_zz_x_xy, g_y_0_x_0_0_zz_x_xz, g_y_0_x_0_0_zz_x_yy, g_y_0_x_0_0_zz_x_yz, g_y_0_x_0_0_zz_x_zz, g_y_zz_0_xx, g_y_zz_0_xy, g_y_zz_0_xz, g_y_zz_0_yy, g_y_zz_0_yz, g_y_zz_0_zz, g_y_zz_xx_xx, g_y_zz_xx_xy, g_y_zz_xx_xz, g_y_zz_xx_yy, g_y_zz_xx_yz, g_y_zz_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_zz_x_xx[i] = -2.0 * g_y_zz_0_xx[i] * a_exp + 4.0 * g_y_zz_xx_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_x_xy[i] = -2.0 * g_y_zz_0_xy[i] * a_exp + 4.0 * g_y_zz_xx_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_x_xz[i] = -2.0 * g_y_zz_0_xz[i] * a_exp + 4.0 * g_y_zz_xx_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_x_yy[i] = -2.0 * g_y_zz_0_yy[i] * a_exp + 4.0 * g_y_zz_xx_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_x_yz[i] = -2.0 * g_y_zz_0_yz[i] * a_exp + 4.0 * g_y_zz_xx_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_x_zz[i] = -2.0 * g_y_zz_0_zz[i] * a_exp + 4.0 * g_y_zz_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_y_0_x_0_0_zz_y_xx, g_y_0_x_0_0_zz_y_xy, g_y_0_x_0_0_zz_y_xz, g_y_0_x_0_0_zz_y_yy, g_y_0_x_0_0_zz_y_yz, g_y_0_x_0_0_zz_y_zz, g_y_zz_xy_xx, g_y_zz_xy_xy, g_y_zz_xy_xz, g_y_zz_xy_yy, g_y_zz_xy_yz, g_y_zz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_zz_y_xx[i] = 4.0 * g_y_zz_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_y_xy[i] = 4.0 * g_y_zz_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_y_xz[i] = 4.0 * g_y_zz_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_y_yy[i] = 4.0 * g_y_zz_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_y_yz[i] = 4.0 * g_y_zz_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_y_zz[i] = 4.0 * g_y_zz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_y_0_x_0_0_zz_z_xx, g_y_0_x_0_0_zz_z_xy, g_y_0_x_0_0_zz_z_xz, g_y_0_x_0_0_zz_z_yy, g_y_0_x_0_0_zz_z_yz, g_y_0_x_0_0_zz_z_zz, g_y_zz_xz_xx, g_y_zz_xz_xy, g_y_zz_xz_xz, g_y_zz_xz_yy, g_y_zz_xz_yz, g_y_zz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_zz_z_xx[i] = 4.0 * g_y_zz_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_z_xy[i] = 4.0 * g_y_zz_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_z_xz[i] = 4.0 * g_y_zz_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_z_yy[i] = 4.0 * g_y_zz_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_z_yz[i] = 4.0 * g_y_zz_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_z_zz[i] = 4.0 * g_y_zz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_y_0_y_0_0_xx_x_xx, g_y_0_y_0_0_xx_x_xy, g_y_0_y_0_0_xx_x_xz, g_y_0_y_0_0_xx_x_yy, g_y_0_y_0_0_xx_x_yz, g_y_0_y_0_0_xx_x_zz, g_y_xx_xy_xx, g_y_xx_xy_xy, g_y_xx_xy_xz, g_y_xx_xy_yy, g_y_xx_xy_yz, g_y_xx_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xx_x_xx[i] = 4.0 * g_y_xx_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_x_xy[i] = 4.0 * g_y_xx_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_x_xz[i] = 4.0 * g_y_xx_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_x_yy[i] = 4.0 * g_y_xx_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_x_yz[i] = 4.0 * g_y_xx_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_x_zz[i] = 4.0 * g_y_xx_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_y_0_y_0_0_xx_y_xx, g_y_0_y_0_0_xx_y_xy, g_y_0_y_0_0_xx_y_xz, g_y_0_y_0_0_xx_y_yy, g_y_0_y_0_0_xx_y_yz, g_y_0_y_0_0_xx_y_zz, g_y_xx_0_xx, g_y_xx_0_xy, g_y_xx_0_xz, g_y_xx_0_yy, g_y_xx_0_yz, g_y_xx_0_zz, g_y_xx_yy_xx, g_y_xx_yy_xy, g_y_xx_yy_xz, g_y_xx_yy_yy, g_y_xx_yy_yz, g_y_xx_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xx_y_xx[i] = -2.0 * g_y_xx_0_xx[i] * a_exp + 4.0 * g_y_xx_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_y_xy[i] = -2.0 * g_y_xx_0_xy[i] * a_exp + 4.0 * g_y_xx_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_y_xz[i] = -2.0 * g_y_xx_0_xz[i] * a_exp + 4.0 * g_y_xx_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_y_yy[i] = -2.0 * g_y_xx_0_yy[i] * a_exp + 4.0 * g_y_xx_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_y_yz[i] = -2.0 * g_y_xx_0_yz[i] * a_exp + 4.0 * g_y_xx_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_y_zz[i] = -2.0 * g_y_xx_0_zz[i] * a_exp + 4.0 * g_y_xx_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_y_0_y_0_0_xx_z_xx, g_y_0_y_0_0_xx_z_xy, g_y_0_y_0_0_xx_z_xz, g_y_0_y_0_0_xx_z_yy, g_y_0_y_0_0_xx_z_yz, g_y_0_y_0_0_xx_z_zz, g_y_xx_yz_xx, g_y_xx_yz_xy, g_y_xx_yz_xz, g_y_xx_yz_yy, g_y_xx_yz_yz, g_y_xx_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xx_z_xx[i] = 4.0 * g_y_xx_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_z_xy[i] = 4.0 * g_y_xx_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_z_xz[i] = 4.0 * g_y_xx_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_z_yy[i] = 4.0 * g_y_xx_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_z_yz[i] = 4.0 * g_y_xx_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_z_zz[i] = 4.0 * g_y_xx_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_y_0_y_0_0_xy_x_xx, g_y_0_y_0_0_xy_x_xy, g_y_0_y_0_0_xy_x_xz, g_y_0_y_0_0_xy_x_yy, g_y_0_y_0_0_xy_x_yz, g_y_0_y_0_0_xy_x_zz, g_y_xy_xy_xx, g_y_xy_xy_xy, g_y_xy_xy_xz, g_y_xy_xy_yy, g_y_xy_xy_yz, g_y_xy_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xy_x_xx[i] = 4.0 * g_y_xy_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_x_xy[i] = 4.0 * g_y_xy_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_x_xz[i] = 4.0 * g_y_xy_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_x_yy[i] = 4.0 * g_y_xy_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_x_yz[i] = 4.0 * g_y_xy_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_x_zz[i] = 4.0 * g_y_xy_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_y_0_y_0_0_xy_y_xx, g_y_0_y_0_0_xy_y_xy, g_y_0_y_0_0_xy_y_xz, g_y_0_y_0_0_xy_y_yy, g_y_0_y_0_0_xy_y_yz, g_y_0_y_0_0_xy_y_zz, g_y_xy_0_xx, g_y_xy_0_xy, g_y_xy_0_xz, g_y_xy_0_yy, g_y_xy_0_yz, g_y_xy_0_zz, g_y_xy_yy_xx, g_y_xy_yy_xy, g_y_xy_yy_xz, g_y_xy_yy_yy, g_y_xy_yy_yz, g_y_xy_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xy_y_xx[i] = -2.0 * g_y_xy_0_xx[i] * a_exp + 4.0 * g_y_xy_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_y_xy[i] = -2.0 * g_y_xy_0_xy[i] * a_exp + 4.0 * g_y_xy_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_y_xz[i] = -2.0 * g_y_xy_0_xz[i] * a_exp + 4.0 * g_y_xy_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_y_yy[i] = -2.0 * g_y_xy_0_yy[i] * a_exp + 4.0 * g_y_xy_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_y_yz[i] = -2.0 * g_y_xy_0_yz[i] * a_exp + 4.0 * g_y_xy_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_y_zz[i] = -2.0 * g_y_xy_0_zz[i] * a_exp + 4.0 * g_y_xy_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_y_0_y_0_0_xy_z_xx, g_y_0_y_0_0_xy_z_xy, g_y_0_y_0_0_xy_z_xz, g_y_0_y_0_0_xy_z_yy, g_y_0_y_0_0_xy_z_yz, g_y_0_y_0_0_xy_z_zz, g_y_xy_yz_xx, g_y_xy_yz_xy, g_y_xy_yz_xz, g_y_xy_yz_yy, g_y_xy_yz_yz, g_y_xy_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xy_z_xx[i] = 4.0 * g_y_xy_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_z_xy[i] = 4.0 * g_y_xy_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_z_xz[i] = 4.0 * g_y_xy_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_z_yy[i] = 4.0 * g_y_xy_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_z_yz[i] = 4.0 * g_y_xy_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_z_zz[i] = 4.0 * g_y_xy_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_y_0_y_0_0_xz_x_xx, g_y_0_y_0_0_xz_x_xy, g_y_0_y_0_0_xz_x_xz, g_y_0_y_0_0_xz_x_yy, g_y_0_y_0_0_xz_x_yz, g_y_0_y_0_0_xz_x_zz, g_y_xz_xy_xx, g_y_xz_xy_xy, g_y_xz_xy_xz, g_y_xz_xy_yy, g_y_xz_xy_yz, g_y_xz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xz_x_xx[i] = 4.0 * g_y_xz_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_x_xy[i] = 4.0 * g_y_xz_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_x_xz[i] = 4.0 * g_y_xz_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_x_yy[i] = 4.0 * g_y_xz_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_x_yz[i] = 4.0 * g_y_xz_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_x_zz[i] = 4.0 * g_y_xz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_y_0_y_0_0_xz_y_xx, g_y_0_y_0_0_xz_y_xy, g_y_0_y_0_0_xz_y_xz, g_y_0_y_0_0_xz_y_yy, g_y_0_y_0_0_xz_y_yz, g_y_0_y_0_0_xz_y_zz, g_y_xz_0_xx, g_y_xz_0_xy, g_y_xz_0_xz, g_y_xz_0_yy, g_y_xz_0_yz, g_y_xz_0_zz, g_y_xz_yy_xx, g_y_xz_yy_xy, g_y_xz_yy_xz, g_y_xz_yy_yy, g_y_xz_yy_yz, g_y_xz_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xz_y_xx[i] = -2.0 * g_y_xz_0_xx[i] * a_exp + 4.0 * g_y_xz_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_y_xy[i] = -2.0 * g_y_xz_0_xy[i] * a_exp + 4.0 * g_y_xz_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_y_xz[i] = -2.0 * g_y_xz_0_xz[i] * a_exp + 4.0 * g_y_xz_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_y_yy[i] = -2.0 * g_y_xz_0_yy[i] * a_exp + 4.0 * g_y_xz_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_y_yz[i] = -2.0 * g_y_xz_0_yz[i] * a_exp + 4.0 * g_y_xz_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_y_zz[i] = -2.0 * g_y_xz_0_zz[i] * a_exp + 4.0 * g_y_xz_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_y_0_y_0_0_xz_z_xx, g_y_0_y_0_0_xz_z_xy, g_y_0_y_0_0_xz_z_xz, g_y_0_y_0_0_xz_z_yy, g_y_0_y_0_0_xz_z_yz, g_y_0_y_0_0_xz_z_zz, g_y_xz_yz_xx, g_y_xz_yz_xy, g_y_xz_yz_xz, g_y_xz_yz_yy, g_y_xz_yz_yz, g_y_xz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xz_z_xx[i] = 4.0 * g_y_xz_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_z_xy[i] = 4.0 * g_y_xz_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_z_xz[i] = 4.0 * g_y_xz_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_z_yy[i] = 4.0 * g_y_xz_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_z_yz[i] = 4.0 * g_y_xz_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_z_zz[i] = 4.0 * g_y_xz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_y_0_y_0_0_yy_x_xx, g_y_0_y_0_0_yy_x_xy, g_y_0_y_0_0_yy_x_xz, g_y_0_y_0_0_yy_x_yy, g_y_0_y_0_0_yy_x_yz, g_y_0_y_0_0_yy_x_zz, g_y_yy_xy_xx, g_y_yy_xy_xy, g_y_yy_xy_xz, g_y_yy_xy_yy, g_y_yy_xy_yz, g_y_yy_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yy_x_xx[i] = 4.0 * g_y_yy_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_x_xy[i] = 4.0 * g_y_yy_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_x_xz[i] = 4.0 * g_y_yy_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_x_yy[i] = 4.0 * g_y_yy_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_x_yz[i] = 4.0 * g_y_yy_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_x_zz[i] = 4.0 * g_y_yy_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_y_0_y_0_0_yy_y_xx, g_y_0_y_0_0_yy_y_xy, g_y_0_y_0_0_yy_y_xz, g_y_0_y_0_0_yy_y_yy, g_y_0_y_0_0_yy_y_yz, g_y_0_y_0_0_yy_y_zz, g_y_yy_0_xx, g_y_yy_0_xy, g_y_yy_0_xz, g_y_yy_0_yy, g_y_yy_0_yz, g_y_yy_0_zz, g_y_yy_yy_xx, g_y_yy_yy_xy, g_y_yy_yy_xz, g_y_yy_yy_yy, g_y_yy_yy_yz, g_y_yy_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yy_y_xx[i] = -2.0 * g_y_yy_0_xx[i] * a_exp + 4.0 * g_y_yy_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_y_xy[i] = -2.0 * g_y_yy_0_xy[i] * a_exp + 4.0 * g_y_yy_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_y_xz[i] = -2.0 * g_y_yy_0_xz[i] * a_exp + 4.0 * g_y_yy_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_y_yy[i] = -2.0 * g_y_yy_0_yy[i] * a_exp + 4.0 * g_y_yy_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_y_yz[i] = -2.0 * g_y_yy_0_yz[i] * a_exp + 4.0 * g_y_yy_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_y_zz[i] = -2.0 * g_y_yy_0_zz[i] * a_exp + 4.0 * g_y_yy_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_y_0_y_0_0_yy_z_xx, g_y_0_y_0_0_yy_z_xy, g_y_0_y_0_0_yy_z_xz, g_y_0_y_0_0_yy_z_yy, g_y_0_y_0_0_yy_z_yz, g_y_0_y_0_0_yy_z_zz, g_y_yy_yz_xx, g_y_yy_yz_xy, g_y_yy_yz_xz, g_y_yy_yz_yy, g_y_yy_yz_yz, g_y_yy_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yy_z_xx[i] = 4.0 * g_y_yy_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_z_xy[i] = 4.0 * g_y_yy_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_z_xz[i] = 4.0 * g_y_yy_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_z_yy[i] = 4.0 * g_y_yy_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_z_yz[i] = 4.0 * g_y_yy_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_z_zz[i] = 4.0 * g_y_yy_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_y_0_y_0_0_yz_x_xx, g_y_0_y_0_0_yz_x_xy, g_y_0_y_0_0_yz_x_xz, g_y_0_y_0_0_yz_x_yy, g_y_0_y_0_0_yz_x_yz, g_y_0_y_0_0_yz_x_zz, g_y_yz_xy_xx, g_y_yz_xy_xy, g_y_yz_xy_xz, g_y_yz_xy_yy, g_y_yz_xy_yz, g_y_yz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yz_x_xx[i] = 4.0 * g_y_yz_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_x_xy[i] = 4.0 * g_y_yz_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_x_xz[i] = 4.0 * g_y_yz_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_x_yy[i] = 4.0 * g_y_yz_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_x_yz[i] = 4.0 * g_y_yz_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_x_zz[i] = 4.0 * g_y_yz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_y_0_y_0_0_yz_y_xx, g_y_0_y_0_0_yz_y_xy, g_y_0_y_0_0_yz_y_xz, g_y_0_y_0_0_yz_y_yy, g_y_0_y_0_0_yz_y_yz, g_y_0_y_0_0_yz_y_zz, g_y_yz_0_xx, g_y_yz_0_xy, g_y_yz_0_xz, g_y_yz_0_yy, g_y_yz_0_yz, g_y_yz_0_zz, g_y_yz_yy_xx, g_y_yz_yy_xy, g_y_yz_yy_xz, g_y_yz_yy_yy, g_y_yz_yy_yz, g_y_yz_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yz_y_xx[i] = -2.0 * g_y_yz_0_xx[i] * a_exp + 4.0 * g_y_yz_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_y_xy[i] = -2.0 * g_y_yz_0_xy[i] * a_exp + 4.0 * g_y_yz_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_y_xz[i] = -2.0 * g_y_yz_0_xz[i] * a_exp + 4.0 * g_y_yz_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_y_yy[i] = -2.0 * g_y_yz_0_yy[i] * a_exp + 4.0 * g_y_yz_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_y_yz[i] = -2.0 * g_y_yz_0_yz[i] * a_exp + 4.0 * g_y_yz_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_y_zz[i] = -2.0 * g_y_yz_0_zz[i] * a_exp + 4.0 * g_y_yz_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_y_0_y_0_0_yz_z_xx, g_y_0_y_0_0_yz_z_xy, g_y_0_y_0_0_yz_z_xz, g_y_0_y_0_0_yz_z_yy, g_y_0_y_0_0_yz_z_yz, g_y_0_y_0_0_yz_z_zz, g_y_yz_yz_xx, g_y_yz_yz_xy, g_y_yz_yz_xz, g_y_yz_yz_yy, g_y_yz_yz_yz, g_y_yz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yz_z_xx[i] = 4.0 * g_y_yz_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_z_xy[i] = 4.0 * g_y_yz_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_z_xz[i] = 4.0 * g_y_yz_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_z_yy[i] = 4.0 * g_y_yz_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_z_yz[i] = 4.0 * g_y_yz_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_z_zz[i] = 4.0 * g_y_yz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_y_0_y_0_0_zz_x_xx, g_y_0_y_0_0_zz_x_xy, g_y_0_y_0_0_zz_x_xz, g_y_0_y_0_0_zz_x_yy, g_y_0_y_0_0_zz_x_yz, g_y_0_y_0_0_zz_x_zz, g_y_zz_xy_xx, g_y_zz_xy_xy, g_y_zz_xy_xz, g_y_zz_xy_yy, g_y_zz_xy_yz, g_y_zz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_zz_x_xx[i] = 4.0 * g_y_zz_xy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_x_xy[i] = 4.0 * g_y_zz_xy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_x_xz[i] = 4.0 * g_y_zz_xy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_x_yy[i] = 4.0 * g_y_zz_xy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_x_yz[i] = 4.0 * g_y_zz_xy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_x_zz[i] = 4.0 * g_y_zz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_y_0_y_0_0_zz_y_xx, g_y_0_y_0_0_zz_y_xy, g_y_0_y_0_0_zz_y_xz, g_y_0_y_0_0_zz_y_yy, g_y_0_y_0_0_zz_y_yz, g_y_0_y_0_0_zz_y_zz, g_y_zz_0_xx, g_y_zz_0_xy, g_y_zz_0_xz, g_y_zz_0_yy, g_y_zz_0_yz, g_y_zz_0_zz, g_y_zz_yy_xx, g_y_zz_yy_xy, g_y_zz_yy_xz, g_y_zz_yy_yy, g_y_zz_yy_yz, g_y_zz_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_zz_y_xx[i] = -2.0 * g_y_zz_0_xx[i] * a_exp + 4.0 * g_y_zz_yy_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_y_xy[i] = -2.0 * g_y_zz_0_xy[i] * a_exp + 4.0 * g_y_zz_yy_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_y_xz[i] = -2.0 * g_y_zz_0_xz[i] * a_exp + 4.0 * g_y_zz_yy_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_y_yy[i] = -2.0 * g_y_zz_0_yy[i] * a_exp + 4.0 * g_y_zz_yy_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_y_yz[i] = -2.0 * g_y_zz_0_yz[i] * a_exp + 4.0 * g_y_zz_yy_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_y_zz[i] = -2.0 * g_y_zz_0_zz[i] * a_exp + 4.0 * g_y_zz_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_y_0_y_0_0_zz_z_xx, g_y_0_y_0_0_zz_z_xy, g_y_0_y_0_0_zz_z_xz, g_y_0_y_0_0_zz_z_yy, g_y_0_y_0_0_zz_z_yz, g_y_0_y_0_0_zz_z_zz, g_y_zz_yz_xx, g_y_zz_yz_xy, g_y_zz_yz_xz, g_y_zz_yz_yy, g_y_zz_yz_yz, g_y_zz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_zz_z_xx[i] = 4.0 * g_y_zz_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_z_xy[i] = 4.0 * g_y_zz_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_z_xz[i] = 4.0 * g_y_zz_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_z_yy[i] = 4.0 * g_y_zz_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_z_yz[i] = 4.0 * g_y_zz_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_z_zz[i] = 4.0 * g_y_zz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_y_0_z_0_0_xx_x_xx, g_y_0_z_0_0_xx_x_xy, g_y_0_z_0_0_xx_x_xz, g_y_0_z_0_0_xx_x_yy, g_y_0_z_0_0_xx_x_yz, g_y_0_z_0_0_xx_x_zz, g_y_xx_xz_xx, g_y_xx_xz_xy, g_y_xx_xz_xz, g_y_xx_xz_yy, g_y_xx_xz_yz, g_y_xx_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xx_x_xx[i] = 4.0 * g_y_xx_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_x_xy[i] = 4.0 * g_y_xx_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_x_xz[i] = 4.0 * g_y_xx_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_x_yy[i] = 4.0 * g_y_xx_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_x_yz[i] = 4.0 * g_y_xx_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_x_zz[i] = 4.0 * g_y_xx_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_y_0_z_0_0_xx_y_xx, g_y_0_z_0_0_xx_y_xy, g_y_0_z_0_0_xx_y_xz, g_y_0_z_0_0_xx_y_yy, g_y_0_z_0_0_xx_y_yz, g_y_0_z_0_0_xx_y_zz, g_y_xx_yz_xx, g_y_xx_yz_xy, g_y_xx_yz_xz, g_y_xx_yz_yy, g_y_xx_yz_yz, g_y_xx_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xx_y_xx[i] = 4.0 * g_y_xx_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_y_xy[i] = 4.0 * g_y_xx_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_y_xz[i] = 4.0 * g_y_xx_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_y_yy[i] = 4.0 * g_y_xx_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_y_yz[i] = 4.0 * g_y_xx_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_y_zz[i] = 4.0 * g_y_xx_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_y_0_z_0_0_xx_z_xx, g_y_0_z_0_0_xx_z_xy, g_y_0_z_0_0_xx_z_xz, g_y_0_z_0_0_xx_z_yy, g_y_0_z_0_0_xx_z_yz, g_y_0_z_0_0_xx_z_zz, g_y_xx_0_xx, g_y_xx_0_xy, g_y_xx_0_xz, g_y_xx_0_yy, g_y_xx_0_yz, g_y_xx_0_zz, g_y_xx_zz_xx, g_y_xx_zz_xy, g_y_xx_zz_xz, g_y_xx_zz_yy, g_y_xx_zz_yz, g_y_xx_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xx_z_xx[i] = -2.0 * g_y_xx_0_xx[i] * a_exp + 4.0 * g_y_xx_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_z_xy[i] = -2.0 * g_y_xx_0_xy[i] * a_exp + 4.0 * g_y_xx_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_z_xz[i] = -2.0 * g_y_xx_0_xz[i] * a_exp + 4.0 * g_y_xx_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_z_yy[i] = -2.0 * g_y_xx_0_yy[i] * a_exp + 4.0 * g_y_xx_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_z_yz[i] = -2.0 * g_y_xx_0_yz[i] * a_exp + 4.0 * g_y_xx_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_z_zz[i] = -2.0 * g_y_xx_0_zz[i] * a_exp + 4.0 * g_y_xx_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_y_0_z_0_0_xy_x_xx, g_y_0_z_0_0_xy_x_xy, g_y_0_z_0_0_xy_x_xz, g_y_0_z_0_0_xy_x_yy, g_y_0_z_0_0_xy_x_yz, g_y_0_z_0_0_xy_x_zz, g_y_xy_xz_xx, g_y_xy_xz_xy, g_y_xy_xz_xz, g_y_xy_xz_yy, g_y_xy_xz_yz, g_y_xy_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xy_x_xx[i] = 4.0 * g_y_xy_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_x_xy[i] = 4.0 * g_y_xy_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_x_xz[i] = 4.0 * g_y_xy_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_x_yy[i] = 4.0 * g_y_xy_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_x_yz[i] = 4.0 * g_y_xy_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_x_zz[i] = 4.0 * g_y_xy_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_y_0_z_0_0_xy_y_xx, g_y_0_z_0_0_xy_y_xy, g_y_0_z_0_0_xy_y_xz, g_y_0_z_0_0_xy_y_yy, g_y_0_z_0_0_xy_y_yz, g_y_0_z_0_0_xy_y_zz, g_y_xy_yz_xx, g_y_xy_yz_xy, g_y_xy_yz_xz, g_y_xy_yz_yy, g_y_xy_yz_yz, g_y_xy_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xy_y_xx[i] = 4.0 * g_y_xy_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_y_xy[i] = 4.0 * g_y_xy_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_y_xz[i] = 4.0 * g_y_xy_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_y_yy[i] = 4.0 * g_y_xy_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_y_yz[i] = 4.0 * g_y_xy_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_y_zz[i] = 4.0 * g_y_xy_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_y_0_z_0_0_xy_z_xx, g_y_0_z_0_0_xy_z_xy, g_y_0_z_0_0_xy_z_xz, g_y_0_z_0_0_xy_z_yy, g_y_0_z_0_0_xy_z_yz, g_y_0_z_0_0_xy_z_zz, g_y_xy_0_xx, g_y_xy_0_xy, g_y_xy_0_xz, g_y_xy_0_yy, g_y_xy_0_yz, g_y_xy_0_zz, g_y_xy_zz_xx, g_y_xy_zz_xy, g_y_xy_zz_xz, g_y_xy_zz_yy, g_y_xy_zz_yz, g_y_xy_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xy_z_xx[i] = -2.0 * g_y_xy_0_xx[i] * a_exp + 4.0 * g_y_xy_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_z_xy[i] = -2.0 * g_y_xy_0_xy[i] * a_exp + 4.0 * g_y_xy_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_z_xz[i] = -2.0 * g_y_xy_0_xz[i] * a_exp + 4.0 * g_y_xy_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_z_yy[i] = -2.0 * g_y_xy_0_yy[i] * a_exp + 4.0 * g_y_xy_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_z_yz[i] = -2.0 * g_y_xy_0_yz[i] * a_exp + 4.0 * g_y_xy_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_z_zz[i] = -2.0 * g_y_xy_0_zz[i] * a_exp + 4.0 * g_y_xy_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_y_0_z_0_0_xz_x_xx, g_y_0_z_0_0_xz_x_xy, g_y_0_z_0_0_xz_x_xz, g_y_0_z_0_0_xz_x_yy, g_y_0_z_0_0_xz_x_yz, g_y_0_z_0_0_xz_x_zz, g_y_xz_xz_xx, g_y_xz_xz_xy, g_y_xz_xz_xz, g_y_xz_xz_yy, g_y_xz_xz_yz, g_y_xz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xz_x_xx[i] = 4.0 * g_y_xz_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_x_xy[i] = 4.0 * g_y_xz_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_x_xz[i] = 4.0 * g_y_xz_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_x_yy[i] = 4.0 * g_y_xz_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_x_yz[i] = 4.0 * g_y_xz_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_x_zz[i] = 4.0 * g_y_xz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_y_0_z_0_0_xz_y_xx, g_y_0_z_0_0_xz_y_xy, g_y_0_z_0_0_xz_y_xz, g_y_0_z_0_0_xz_y_yy, g_y_0_z_0_0_xz_y_yz, g_y_0_z_0_0_xz_y_zz, g_y_xz_yz_xx, g_y_xz_yz_xy, g_y_xz_yz_xz, g_y_xz_yz_yy, g_y_xz_yz_yz, g_y_xz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xz_y_xx[i] = 4.0 * g_y_xz_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_y_xy[i] = 4.0 * g_y_xz_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_y_xz[i] = 4.0 * g_y_xz_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_y_yy[i] = 4.0 * g_y_xz_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_y_yz[i] = 4.0 * g_y_xz_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_y_zz[i] = 4.0 * g_y_xz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_y_0_z_0_0_xz_z_xx, g_y_0_z_0_0_xz_z_xy, g_y_0_z_0_0_xz_z_xz, g_y_0_z_0_0_xz_z_yy, g_y_0_z_0_0_xz_z_yz, g_y_0_z_0_0_xz_z_zz, g_y_xz_0_xx, g_y_xz_0_xy, g_y_xz_0_xz, g_y_xz_0_yy, g_y_xz_0_yz, g_y_xz_0_zz, g_y_xz_zz_xx, g_y_xz_zz_xy, g_y_xz_zz_xz, g_y_xz_zz_yy, g_y_xz_zz_yz, g_y_xz_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xz_z_xx[i] = -2.0 * g_y_xz_0_xx[i] * a_exp + 4.0 * g_y_xz_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_z_xy[i] = -2.0 * g_y_xz_0_xy[i] * a_exp + 4.0 * g_y_xz_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_z_xz[i] = -2.0 * g_y_xz_0_xz[i] * a_exp + 4.0 * g_y_xz_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_z_yy[i] = -2.0 * g_y_xz_0_yy[i] * a_exp + 4.0 * g_y_xz_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_z_yz[i] = -2.0 * g_y_xz_0_yz[i] * a_exp + 4.0 * g_y_xz_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_z_zz[i] = -2.0 * g_y_xz_0_zz[i] * a_exp + 4.0 * g_y_xz_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_y_0_z_0_0_yy_x_xx, g_y_0_z_0_0_yy_x_xy, g_y_0_z_0_0_yy_x_xz, g_y_0_z_0_0_yy_x_yy, g_y_0_z_0_0_yy_x_yz, g_y_0_z_0_0_yy_x_zz, g_y_yy_xz_xx, g_y_yy_xz_xy, g_y_yy_xz_xz, g_y_yy_xz_yy, g_y_yy_xz_yz, g_y_yy_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yy_x_xx[i] = 4.0 * g_y_yy_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_x_xy[i] = 4.0 * g_y_yy_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_x_xz[i] = 4.0 * g_y_yy_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_x_yy[i] = 4.0 * g_y_yy_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_x_yz[i] = 4.0 * g_y_yy_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_x_zz[i] = 4.0 * g_y_yy_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_y_0_z_0_0_yy_y_xx, g_y_0_z_0_0_yy_y_xy, g_y_0_z_0_0_yy_y_xz, g_y_0_z_0_0_yy_y_yy, g_y_0_z_0_0_yy_y_yz, g_y_0_z_0_0_yy_y_zz, g_y_yy_yz_xx, g_y_yy_yz_xy, g_y_yy_yz_xz, g_y_yy_yz_yy, g_y_yy_yz_yz, g_y_yy_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yy_y_xx[i] = 4.0 * g_y_yy_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_y_xy[i] = 4.0 * g_y_yy_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_y_xz[i] = 4.0 * g_y_yy_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_y_yy[i] = 4.0 * g_y_yy_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_y_yz[i] = 4.0 * g_y_yy_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_y_zz[i] = 4.0 * g_y_yy_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_y_0_z_0_0_yy_z_xx, g_y_0_z_0_0_yy_z_xy, g_y_0_z_0_0_yy_z_xz, g_y_0_z_0_0_yy_z_yy, g_y_0_z_0_0_yy_z_yz, g_y_0_z_0_0_yy_z_zz, g_y_yy_0_xx, g_y_yy_0_xy, g_y_yy_0_xz, g_y_yy_0_yy, g_y_yy_0_yz, g_y_yy_0_zz, g_y_yy_zz_xx, g_y_yy_zz_xy, g_y_yy_zz_xz, g_y_yy_zz_yy, g_y_yy_zz_yz, g_y_yy_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yy_z_xx[i] = -2.0 * g_y_yy_0_xx[i] * a_exp + 4.0 * g_y_yy_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_z_xy[i] = -2.0 * g_y_yy_0_xy[i] * a_exp + 4.0 * g_y_yy_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_z_xz[i] = -2.0 * g_y_yy_0_xz[i] * a_exp + 4.0 * g_y_yy_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_z_yy[i] = -2.0 * g_y_yy_0_yy[i] * a_exp + 4.0 * g_y_yy_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_z_yz[i] = -2.0 * g_y_yy_0_yz[i] * a_exp + 4.0 * g_y_yy_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_z_zz[i] = -2.0 * g_y_yy_0_zz[i] * a_exp + 4.0 * g_y_yy_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_y_0_z_0_0_yz_x_xx, g_y_0_z_0_0_yz_x_xy, g_y_0_z_0_0_yz_x_xz, g_y_0_z_0_0_yz_x_yy, g_y_0_z_0_0_yz_x_yz, g_y_0_z_0_0_yz_x_zz, g_y_yz_xz_xx, g_y_yz_xz_xy, g_y_yz_xz_xz, g_y_yz_xz_yy, g_y_yz_xz_yz, g_y_yz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yz_x_xx[i] = 4.0 * g_y_yz_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_x_xy[i] = 4.0 * g_y_yz_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_x_xz[i] = 4.0 * g_y_yz_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_x_yy[i] = 4.0 * g_y_yz_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_x_yz[i] = 4.0 * g_y_yz_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_x_zz[i] = 4.0 * g_y_yz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_y_0_z_0_0_yz_y_xx, g_y_0_z_0_0_yz_y_xy, g_y_0_z_0_0_yz_y_xz, g_y_0_z_0_0_yz_y_yy, g_y_0_z_0_0_yz_y_yz, g_y_0_z_0_0_yz_y_zz, g_y_yz_yz_xx, g_y_yz_yz_xy, g_y_yz_yz_xz, g_y_yz_yz_yy, g_y_yz_yz_yz, g_y_yz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yz_y_xx[i] = 4.0 * g_y_yz_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_y_xy[i] = 4.0 * g_y_yz_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_y_xz[i] = 4.0 * g_y_yz_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_y_yy[i] = 4.0 * g_y_yz_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_y_yz[i] = 4.0 * g_y_yz_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_y_zz[i] = 4.0 * g_y_yz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_y_0_z_0_0_yz_z_xx, g_y_0_z_0_0_yz_z_xy, g_y_0_z_0_0_yz_z_xz, g_y_0_z_0_0_yz_z_yy, g_y_0_z_0_0_yz_z_yz, g_y_0_z_0_0_yz_z_zz, g_y_yz_0_xx, g_y_yz_0_xy, g_y_yz_0_xz, g_y_yz_0_yy, g_y_yz_0_yz, g_y_yz_0_zz, g_y_yz_zz_xx, g_y_yz_zz_xy, g_y_yz_zz_xz, g_y_yz_zz_yy, g_y_yz_zz_yz, g_y_yz_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yz_z_xx[i] = -2.0 * g_y_yz_0_xx[i] * a_exp + 4.0 * g_y_yz_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_z_xy[i] = -2.0 * g_y_yz_0_xy[i] * a_exp + 4.0 * g_y_yz_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_z_xz[i] = -2.0 * g_y_yz_0_xz[i] * a_exp + 4.0 * g_y_yz_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_z_yy[i] = -2.0 * g_y_yz_0_yy[i] * a_exp + 4.0 * g_y_yz_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_z_yz[i] = -2.0 * g_y_yz_0_yz[i] * a_exp + 4.0 * g_y_yz_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_z_zz[i] = -2.0 * g_y_yz_0_zz[i] * a_exp + 4.0 * g_y_yz_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_y_0_z_0_0_zz_x_xx, g_y_0_z_0_0_zz_x_xy, g_y_0_z_0_0_zz_x_xz, g_y_0_z_0_0_zz_x_yy, g_y_0_z_0_0_zz_x_yz, g_y_0_z_0_0_zz_x_zz, g_y_zz_xz_xx, g_y_zz_xz_xy, g_y_zz_xz_xz, g_y_zz_xz_yy, g_y_zz_xz_yz, g_y_zz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_zz_x_xx[i] = 4.0 * g_y_zz_xz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_x_xy[i] = 4.0 * g_y_zz_xz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_x_xz[i] = 4.0 * g_y_zz_xz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_x_yy[i] = 4.0 * g_y_zz_xz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_x_yz[i] = 4.0 * g_y_zz_xz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_x_zz[i] = 4.0 * g_y_zz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_y_0_z_0_0_zz_y_xx, g_y_0_z_0_0_zz_y_xy, g_y_0_z_0_0_zz_y_xz, g_y_0_z_0_0_zz_y_yy, g_y_0_z_0_0_zz_y_yz, g_y_0_z_0_0_zz_y_zz, g_y_zz_yz_xx, g_y_zz_yz_xy, g_y_zz_yz_xz, g_y_zz_yz_yy, g_y_zz_yz_yz, g_y_zz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_zz_y_xx[i] = 4.0 * g_y_zz_yz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_y_xy[i] = 4.0 * g_y_zz_yz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_y_xz[i] = 4.0 * g_y_zz_yz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_y_yy[i] = 4.0 * g_y_zz_yz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_y_yz[i] = 4.0 * g_y_zz_yz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_y_zz[i] = 4.0 * g_y_zz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_y_0_z_0_0_zz_z_xx, g_y_0_z_0_0_zz_z_xy, g_y_0_z_0_0_zz_z_xz, g_y_0_z_0_0_zz_z_yy, g_y_0_z_0_0_zz_z_yz, g_y_0_z_0_0_zz_z_zz, g_y_zz_0_xx, g_y_zz_0_xy, g_y_zz_0_xz, g_y_zz_0_yy, g_y_zz_0_yz, g_y_zz_0_zz, g_y_zz_zz_xx, g_y_zz_zz_xy, g_y_zz_zz_xz, g_y_zz_zz_yy, g_y_zz_zz_yz, g_y_zz_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_zz_z_xx[i] = -2.0 * g_y_zz_0_xx[i] * a_exp + 4.0 * g_y_zz_zz_xx[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_z_xy[i] = -2.0 * g_y_zz_0_xy[i] * a_exp + 4.0 * g_y_zz_zz_xy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_z_xz[i] = -2.0 * g_y_zz_0_xz[i] * a_exp + 4.0 * g_y_zz_zz_xz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_z_yy[i] = -2.0 * g_y_zz_0_yy[i] * a_exp + 4.0 * g_y_zz_zz_yy[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_z_yz[i] = -2.0 * g_y_zz_0_yz[i] * a_exp + 4.0 * g_y_zz_zz_yz[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_z_zz[i] = -2.0 * g_y_zz_0_zz[i] * a_exp + 4.0 * g_y_zz_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_z_0_x_0_0_xx_x_xx, g_z_0_x_0_0_xx_x_xy, g_z_0_x_0_0_xx_x_xz, g_z_0_x_0_0_xx_x_yy, g_z_0_x_0_0_xx_x_yz, g_z_0_x_0_0_xx_x_zz, g_z_xx_0_xx, g_z_xx_0_xy, g_z_xx_0_xz, g_z_xx_0_yy, g_z_xx_0_yz, g_z_xx_0_zz, g_z_xx_xx_xx, g_z_xx_xx_xy, g_z_xx_xx_xz, g_z_xx_xx_yy, g_z_xx_xx_yz, g_z_xx_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xx_x_xx[i] = -2.0 * g_z_xx_0_xx[i] * a_exp + 4.0 * g_z_xx_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_x_xy[i] = -2.0 * g_z_xx_0_xy[i] * a_exp + 4.0 * g_z_xx_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_x_xz[i] = -2.0 * g_z_xx_0_xz[i] * a_exp + 4.0 * g_z_xx_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_x_yy[i] = -2.0 * g_z_xx_0_yy[i] * a_exp + 4.0 * g_z_xx_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_x_yz[i] = -2.0 * g_z_xx_0_yz[i] * a_exp + 4.0 * g_z_xx_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_x_zz[i] = -2.0 * g_z_xx_0_zz[i] * a_exp + 4.0 * g_z_xx_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_z_0_x_0_0_xx_y_xx, g_z_0_x_0_0_xx_y_xy, g_z_0_x_0_0_xx_y_xz, g_z_0_x_0_0_xx_y_yy, g_z_0_x_0_0_xx_y_yz, g_z_0_x_0_0_xx_y_zz, g_z_xx_xy_xx, g_z_xx_xy_xy, g_z_xx_xy_xz, g_z_xx_xy_yy, g_z_xx_xy_yz, g_z_xx_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xx_y_xx[i] = 4.0 * g_z_xx_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_y_xy[i] = 4.0 * g_z_xx_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_y_xz[i] = 4.0 * g_z_xx_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_y_yy[i] = 4.0 * g_z_xx_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_y_yz[i] = 4.0 * g_z_xx_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_y_zz[i] = 4.0 * g_z_xx_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_z_0_x_0_0_xx_z_xx, g_z_0_x_0_0_xx_z_xy, g_z_0_x_0_0_xx_z_xz, g_z_0_x_0_0_xx_z_yy, g_z_0_x_0_0_xx_z_yz, g_z_0_x_0_0_xx_z_zz, g_z_xx_xz_xx, g_z_xx_xz_xy, g_z_xx_xz_xz, g_z_xx_xz_yy, g_z_xx_xz_yz, g_z_xx_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xx_z_xx[i] = 4.0 * g_z_xx_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_z_xy[i] = 4.0 * g_z_xx_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_z_xz[i] = 4.0 * g_z_xx_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_z_yy[i] = 4.0 * g_z_xx_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_z_yz[i] = 4.0 * g_z_xx_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_z_zz[i] = 4.0 * g_z_xx_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_z_0_x_0_0_xy_x_xx, g_z_0_x_0_0_xy_x_xy, g_z_0_x_0_0_xy_x_xz, g_z_0_x_0_0_xy_x_yy, g_z_0_x_0_0_xy_x_yz, g_z_0_x_0_0_xy_x_zz, g_z_xy_0_xx, g_z_xy_0_xy, g_z_xy_0_xz, g_z_xy_0_yy, g_z_xy_0_yz, g_z_xy_0_zz, g_z_xy_xx_xx, g_z_xy_xx_xy, g_z_xy_xx_xz, g_z_xy_xx_yy, g_z_xy_xx_yz, g_z_xy_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xy_x_xx[i] = -2.0 * g_z_xy_0_xx[i] * a_exp + 4.0 * g_z_xy_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_x_xy[i] = -2.0 * g_z_xy_0_xy[i] * a_exp + 4.0 * g_z_xy_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_x_xz[i] = -2.0 * g_z_xy_0_xz[i] * a_exp + 4.0 * g_z_xy_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_x_yy[i] = -2.0 * g_z_xy_0_yy[i] * a_exp + 4.0 * g_z_xy_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_x_yz[i] = -2.0 * g_z_xy_0_yz[i] * a_exp + 4.0 * g_z_xy_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_x_zz[i] = -2.0 * g_z_xy_0_zz[i] * a_exp + 4.0 * g_z_xy_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_z_0_x_0_0_xy_y_xx, g_z_0_x_0_0_xy_y_xy, g_z_0_x_0_0_xy_y_xz, g_z_0_x_0_0_xy_y_yy, g_z_0_x_0_0_xy_y_yz, g_z_0_x_0_0_xy_y_zz, g_z_xy_xy_xx, g_z_xy_xy_xy, g_z_xy_xy_xz, g_z_xy_xy_yy, g_z_xy_xy_yz, g_z_xy_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xy_y_xx[i] = 4.0 * g_z_xy_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_y_xy[i] = 4.0 * g_z_xy_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_y_xz[i] = 4.0 * g_z_xy_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_y_yy[i] = 4.0 * g_z_xy_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_y_yz[i] = 4.0 * g_z_xy_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_y_zz[i] = 4.0 * g_z_xy_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_z_0_x_0_0_xy_z_xx, g_z_0_x_0_0_xy_z_xy, g_z_0_x_0_0_xy_z_xz, g_z_0_x_0_0_xy_z_yy, g_z_0_x_0_0_xy_z_yz, g_z_0_x_0_0_xy_z_zz, g_z_xy_xz_xx, g_z_xy_xz_xy, g_z_xy_xz_xz, g_z_xy_xz_yy, g_z_xy_xz_yz, g_z_xy_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xy_z_xx[i] = 4.0 * g_z_xy_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_z_xy[i] = 4.0 * g_z_xy_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_z_xz[i] = 4.0 * g_z_xy_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_z_yy[i] = 4.0 * g_z_xy_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_z_yz[i] = 4.0 * g_z_xy_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_z_zz[i] = 4.0 * g_z_xy_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_z_0_x_0_0_xz_x_xx, g_z_0_x_0_0_xz_x_xy, g_z_0_x_0_0_xz_x_xz, g_z_0_x_0_0_xz_x_yy, g_z_0_x_0_0_xz_x_yz, g_z_0_x_0_0_xz_x_zz, g_z_xz_0_xx, g_z_xz_0_xy, g_z_xz_0_xz, g_z_xz_0_yy, g_z_xz_0_yz, g_z_xz_0_zz, g_z_xz_xx_xx, g_z_xz_xx_xy, g_z_xz_xx_xz, g_z_xz_xx_yy, g_z_xz_xx_yz, g_z_xz_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xz_x_xx[i] = -2.0 * g_z_xz_0_xx[i] * a_exp + 4.0 * g_z_xz_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_x_xy[i] = -2.0 * g_z_xz_0_xy[i] * a_exp + 4.0 * g_z_xz_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_x_xz[i] = -2.0 * g_z_xz_0_xz[i] * a_exp + 4.0 * g_z_xz_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_x_yy[i] = -2.0 * g_z_xz_0_yy[i] * a_exp + 4.0 * g_z_xz_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_x_yz[i] = -2.0 * g_z_xz_0_yz[i] * a_exp + 4.0 * g_z_xz_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_x_zz[i] = -2.0 * g_z_xz_0_zz[i] * a_exp + 4.0 * g_z_xz_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_z_0_x_0_0_xz_y_xx, g_z_0_x_0_0_xz_y_xy, g_z_0_x_0_0_xz_y_xz, g_z_0_x_0_0_xz_y_yy, g_z_0_x_0_0_xz_y_yz, g_z_0_x_0_0_xz_y_zz, g_z_xz_xy_xx, g_z_xz_xy_xy, g_z_xz_xy_xz, g_z_xz_xy_yy, g_z_xz_xy_yz, g_z_xz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xz_y_xx[i] = 4.0 * g_z_xz_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_y_xy[i] = 4.0 * g_z_xz_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_y_xz[i] = 4.0 * g_z_xz_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_y_yy[i] = 4.0 * g_z_xz_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_y_yz[i] = 4.0 * g_z_xz_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_y_zz[i] = 4.0 * g_z_xz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_z_0_x_0_0_xz_z_xx, g_z_0_x_0_0_xz_z_xy, g_z_0_x_0_0_xz_z_xz, g_z_0_x_0_0_xz_z_yy, g_z_0_x_0_0_xz_z_yz, g_z_0_x_0_0_xz_z_zz, g_z_xz_xz_xx, g_z_xz_xz_xy, g_z_xz_xz_xz, g_z_xz_xz_yy, g_z_xz_xz_yz, g_z_xz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xz_z_xx[i] = 4.0 * g_z_xz_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_z_xy[i] = 4.0 * g_z_xz_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_z_xz[i] = 4.0 * g_z_xz_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_z_yy[i] = 4.0 * g_z_xz_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_z_yz[i] = 4.0 * g_z_xz_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_z_zz[i] = 4.0 * g_z_xz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_z_0_x_0_0_yy_x_xx, g_z_0_x_0_0_yy_x_xy, g_z_0_x_0_0_yy_x_xz, g_z_0_x_0_0_yy_x_yy, g_z_0_x_0_0_yy_x_yz, g_z_0_x_0_0_yy_x_zz, g_z_yy_0_xx, g_z_yy_0_xy, g_z_yy_0_xz, g_z_yy_0_yy, g_z_yy_0_yz, g_z_yy_0_zz, g_z_yy_xx_xx, g_z_yy_xx_xy, g_z_yy_xx_xz, g_z_yy_xx_yy, g_z_yy_xx_yz, g_z_yy_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yy_x_xx[i] = -2.0 * g_z_yy_0_xx[i] * a_exp + 4.0 * g_z_yy_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_x_xy[i] = -2.0 * g_z_yy_0_xy[i] * a_exp + 4.0 * g_z_yy_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_x_xz[i] = -2.0 * g_z_yy_0_xz[i] * a_exp + 4.0 * g_z_yy_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_x_yy[i] = -2.0 * g_z_yy_0_yy[i] * a_exp + 4.0 * g_z_yy_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_x_yz[i] = -2.0 * g_z_yy_0_yz[i] * a_exp + 4.0 * g_z_yy_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_x_zz[i] = -2.0 * g_z_yy_0_zz[i] * a_exp + 4.0 * g_z_yy_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_z_0_x_0_0_yy_y_xx, g_z_0_x_0_0_yy_y_xy, g_z_0_x_0_0_yy_y_xz, g_z_0_x_0_0_yy_y_yy, g_z_0_x_0_0_yy_y_yz, g_z_0_x_0_0_yy_y_zz, g_z_yy_xy_xx, g_z_yy_xy_xy, g_z_yy_xy_xz, g_z_yy_xy_yy, g_z_yy_xy_yz, g_z_yy_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yy_y_xx[i] = 4.0 * g_z_yy_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_y_xy[i] = 4.0 * g_z_yy_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_y_xz[i] = 4.0 * g_z_yy_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_y_yy[i] = 4.0 * g_z_yy_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_y_yz[i] = 4.0 * g_z_yy_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_y_zz[i] = 4.0 * g_z_yy_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_z_0_x_0_0_yy_z_xx, g_z_0_x_0_0_yy_z_xy, g_z_0_x_0_0_yy_z_xz, g_z_0_x_0_0_yy_z_yy, g_z_0_x_0_0_yy_z_yz, g_z_0_x_0_0_yy_z_zz, g_z_yy_xz_xx, g_z_yy_xz_xy, g_z_yy_xz_xz, g_z_yy_xz_yy, g_z_yy_xz_yz, g_z_yy_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yy_z_xx[i] = 4.0 * g_z_yy_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_z_xy[i] = 4.0 * g_z_yy_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_z_xz[i] = 4.0 * g_z_yy_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_z_yy[i] = 4.0 * g_z_yy_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_z_yz[i] = 4.0 * g_z_yy_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_z_zz[i] = 4.0 * g_z_yy_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_z_0_x_0_0_yz_x_xx, g_z_0_x_0_0_yz_x_xy, g_z_0_x_0_0_yz_x_xz, g_z_0_x_0_0_yz_x_yy, g_z_0_x_0_0_yz_x_yz, g_z_0_x_0_0_yz_x_zz, g_z_yz_0_xx, g_z_yz_0_xy, g_z_yz_0_xz, g_z_yz_0_yy, g_z_yz_0_yz, g_z_yz_0_zz, g_z_yz_xx_xx, g_z_yz_xx_xy, g_z_yz_xx_xz, g_z_yz_xx_yy, g_z_yz_xx_yz, g_z_yz_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yz_x_xx[i] = -2.0 * g_z_yz_0_xx[i] * a_exp + 4.0 * g_z_yz_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_x_xy[i] = -2.0 * g_z_yz_0_xy[i] * a_exp + 4.0 * g_z_yz_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_x_xz[i] = -2.0 * g_z_yz_0_xz[i] * a_exp + 4.0 * g_z_yz_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_x_yy[i] = -2.0 * g_z_yz_0_yy[i] * a_exp + 4.0 * g_z_yz_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_x_yz[i] = -2.0 * g_z_yz_0_yz[i] * a_exp + 4.0 * g_z_yz_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_x_zz[i] = -2.0 * g_z_yz_0_zz[i] * a_exp + 4.0 * g_z_yz_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_z_0_x_0_0_yz_y_xx, g_z_0_x_0_0_yz_y_xy, g_z_0_x_0_0_yz_y_xz, g_z_0_x_0_0_yz_y_yy, g_z_0_x_0_0_yz_y_yz, g_z_0_x_0_0_yz_y_zz, g_z_yz_xy_xx, g_z_yz_xy_xy, g_z_yz_xy_xz, g_z_yz_xy_yy, g_z_yz_xy_yz, g_z_yz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yz_y_xx[i] = 4.0 * g_z_yz_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_y_xy[i] = 4.0 * g_z_yz_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_y_xz[i] = 4.0 * g_z_yz_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_y_yy[i] = 4.0 * g_z_yz_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_y_yz[i] = 4.0 * g_z_yz_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_y_zz[i] = 4.0 * g_z_yz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_z_0_x_0_0_yz_z_xx, g_z_0_x_0_0_yz_z_xy, g_z_0_x_0_0_yz_z_xz, g_z_0_x_0_0_yz_z_yy, g_z_0_x_0_0_yz_z_yz, g_z_0_x_0_0_yz_z_zz, g_z_yz_xz_xx, g_z_yz_xz_xy, g_z_yz_xz_xz, g_z_yz_xz_yy, g_z_yz_xz_yz, g_z_yz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yz_z_xx[i] = 4.0 * g_z_yz_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_z_xy[i] = 4.0 * g_z_yz_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_z_xz[i] = 4.0 * g_z_yz_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_z_yy[i] = 4.0 * g_z_yz_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_z_yz[i] = 4.0 * g_z_yz_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_z_zz[i] = 4.0 * g_z_yz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_z_0_x_0_0_zz_x_xx, g_z_0_x_0_0_zz_x_xy, g_z_0_x_0_0_zz_x_xz, g_z_0_x_0_0_zz_x_yy, g_z_0_x_0_0_zz_x_yz, g_z_0_x_0_0_zz_x_zz, g_z_zz_0_xx, g_z_zz_0_xy, g_z_zz_0_xz, g_z_zz_0_yy, g_z_zz_0_yz, g_z_zz_0_zz, g_z_zz_xx_xx, g_z_zz_xx_xy, g_z_zz_xx_xz, g_z_zz_xx_yy, g_z_zz_xx_yz, g_z_zz_xx_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_zz_x_xx[i] = -2.0 * g_z_zz_0_xx[i] * a_exp + 4.0 * g_z_zz_xx_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_x_xy[i] = -2.0 * g_z_zz_0_xy[i] * a_exp + 4.0 * g_z_zz_xx_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_x_xz[i] = -2.0 * g_z_zz_0_xz[i] * a_exp + 4.0 * g_z_zz_xx_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_x_yy[i] = -2.0 * g_z_zz_0_yy[i] * a_exp + 4.0 * g_z_zz_xx_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_x_yz[i] = -2.0 * g_z_zz_0_yz[i] * a_exp + 4.0 * g_z_zz_xx_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_x_zz[i] = -2.0 * g_z_zz_0_zz[i] * a_exp + 4.0 * g_z_zz_xx_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_z_0_x_0_0_zz_y_xx, g_z_0_x_0_0_zz_y_xy, g_z_0_x_0_0_zz_y_xz, g_z_0_x_0_0_zz_y_yy, g_z_0_x_0_0_zz_y_yz, g_z_0_x_0_0_zz_y_zz, g_z_zz_xy_xx, g_z_zz_xy_xy, g_z_zz_xy_xz, g_z_zz_xy_yy, g_z_zz_xy_yz, g_z_zz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_zz_y_xx[i] = 4.0 * g_z_zz_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_y_xy[i] = 4.0 * g_z_zz_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_y_xz[i] = 4.0 * g_z_zz_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_y_yy[i] = 4.0 * g_z_zz_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_y_yz[i] = 4.0 * g_z_zz_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_y_zz[i] = 4.0 * g_z_zz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_z_0_x_0_0_zz_z_xx, g_z_0_x_0_0_zz_z_xy, g_z_0_x_0_0_zz_z_xz, g_z_0_x_0_0_zz_z_yy, g_z_0_x_0_0_zz_z_yz, g_z_0_x_0_0_zz_z_zz, g_z_zz_xz_xx, g_z_zz_xz_xy, g_z_zz_xz_xz, g_z_zz_xz_yy, g_z_zz_xz_yz, g_z_zz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_zz_z_xx[i] = 4.0 * g_z_zz_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_z_xy[i] = 4.0 * g_z_zz_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_z_xz[i] = 4.0 * g_z_zz_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_z_yy[i] = 4.0 * g_z_zz_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_z_yz[i] = 4.0 * g_z_zz_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_z_zz[i] = 4.0 * g_z_zz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_z_0_y_0_0_xx_x_xx, g_z_0_y_0_0_xx_x_xy, g_z_0_y_0_0_xx_x_xz, g_z_0_y_0_0_xx_x_yy, g_z_0_y_0_0_xx_x_yz, g_z_0_y_0_0_xx_x_zz, g_z_xx_xy_xx, g_z_xx_xy_xy, g_z_xx_xy_xz, g_z_xx_xy_yy, g_z_xx_xy_yz, g_z_xx_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xx_x_xx[i] = 4.0 * g_z_xx_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_x_xy[i] = 4.0 * g_z_xx_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_x_xz[i] = 4.0 * g_z_xx_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_x_yy[i] = 4.0 * g_z_xx_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_x_yz[i] = 4.0 * g_z_xx_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_x_zz[i] = 4.0 * g_z_xx_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_z_0_y_0_0_xx_y_xx, g_z_0_y_0_0_xx_y_xy, g_z_0_y_0_0_xx_y_xz, g_z_0_y_0_0_xx_y_yy, g_z_0_y_0_0_xx_y_yz, g_z_0_y_0_0_xx_y_zz, g_z_xx_0_xx, g_z_xx_0_xy, g_z_xx_0_xz, g_z_xx_0_yy, g_z_xx_0_yz, g_z_xx_0_zz, g_z_xx_yy_xx, g_z_xx_yy_xy, g_z_xx_yy_xz, g_z_xx_yy_yy, g_z_xx_yy_yz, g_z_xx_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xx_y_xx[i] = -2.0 * g_z_xx_0_xx[i] * a_exp + 4.0 * g_z_xx_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_y_xy[i] = -2.0 * g_z_xx_0_xy[i] * a_exp + 4.0 * g_z_xx_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_y_xz[i] = -2.0 * g_z_xx_0_xz[i] * a_exp + 4.0 * g_z_xx_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_y_yy[i] = -2.0 * g_z_xx_0_yy[i] * a_exp + 4.0 * g_z_xx_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_y_yz[i] = -2.0 * g_z_xx_0_yz[i] * a_exp + 4.0 * g_z_xx_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_y_zz[i] = -2.0 * g_z_xx_0_zz[i] * a_exp + 4.0 * g_z_xx_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_z_0_y_0_0_xx_z_xx, g_z_0_y_0_0_xx_z_xy, g_z_0_y_0_0_xx_z_xz, g_z_0_y_0_0_xx_z_yy, g_z_0_y_0_0_xx_z_yz, g_z_0_y_0_0_xx_z_zz, g_z_xx_yz_xx, g_z_xx_yz_xy, g_z_xx_yz_xz, g_z_xx_yz_yy, g_z_xx_yz_yz, g_z_xx_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xx_z_xx[i] = 4.0 * g_z_xx_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_z_xy[i] = 4.0 * g_z_xx_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_z_xz[i] = 4.0 * g_z_xx_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_z_yy[i] = 4.0 * g_z_xx_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_z_yz[i] = 4.0 * g_z_xx_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_z_zz[i] = 4.0 * g_z_xx_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_z_0_y_0_0_xy_x_xx, g_z_0_y_0_0_xy_x_xy, g_z_0_y_0_0_xy_x_xz, g_z_0_y_0_0_xy_x_yy, g_z_0_y_0_0_xy_x_yz, g_z_0_y_0_0_xy_x_zz, g_z_xy_xy_xx, g_z_xy_xy_xy, g_z_xy_xy_xz, g_z_xy_xy_yy, g_z_xy_xy_yz, g_z_xy_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xy_x_xx[i] = 4.0 * g_z_xy_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_x_xy[i] = 4.0 * g_z_xy_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_x_xz[i] = 4.0 * g_z_xy_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_x_yy[i] = 4.0 * g_z_xy_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_x_yz[i] = 4.0 * g_z_xy_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_x_zz[i] = 4.0 * g_z_xy_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_z_0_y_0_0_xy_y_xx, g_z_0_y_0_0_xy_y_xy, g_z_0_y_0_0_xy_y_xz, g_z_0_y_0_0_xy_y_yy, g_z_0_y_0_0_xy_y_yz, g_z_0_y_0_0_xy_y_zz, g_z_xy_0_xx, g_z_xy_0_xy, g_z_xy_0_xz, g_z_xy_0_yy, g_z_xy_0_yz, g_z_xy_0_zz, g_z_xy_yy_xx, g_z_xy_yy_xy, g_z_xy_yy_xz, g_z_xy_yy_yy, g_z_xy_yy_yz, g_z_xy_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xy_y_xx[i] = -2.0 * g_z_xy_0_xx[i] * a_exp + 4.0 * g_z_xy_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_y_xy[i] = -2.0 * g_z_xy_0_xy[i] * a_exp + 4.0 * g_z_xy_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_y_xz[i] = -2.0 * g_z_xy_0_xz[i] * a_exp + 4.0 * g_z_xy_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_y_yy[i] = -2.0 * g_z_xy_0_yy[i] * a_exp + 4.0 * g_z_xy_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_y_yz[i] = -2.0 * g_z_xy_0_yz[i] * a_exp + 4.0 * g_z_xy_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_y_zz[i] = -2.0 * g_z_xy_0_zz[i] * a_exp + 4.0 * g_z_xy_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_z_0_y_0_0_xy_z_xx, g_z_0_y_0_0_xy_z_xy, g_z_0_y_0_0_xy_z_xz, g_z_0_y_0_0_xy_z_yy, g_z_0_y_0_0_xy_z_yz, g_z_0_y_0_0_xy_z_zz, g_z_xy_yz_xx, g_z_xy_yz_xy, g_z_xy_yz_xz, g_z_xy_yz_yy, g_z_xy_yz_yz, g_z_xy_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xy_z_xx[i] = 4.0 * g_z_xy_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_z_xy[i] = 4.0 * g_z_xy_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_z_xz[i] = 4.0 * g_z_xy_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_z_yy[i] = 4.0 * g_z_xy_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_z_yz[i] = 4.0 * g_z_xy_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_z_zz[i] = 4.0 * g_z_xy_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_z_0_y_0_0_xz_x_xx, g_z_0_y_0_0_xz_x_xy, g_z_0_y_0_0_xz_x_xz, g_z_0_y_0_0_xz_x_yy, g_z_0_y_0_0_xz_x_yz, g_z_0_y_0_0_xz_x_zz, g_z_xz_xy_xx, g_z_xz_xy_xy, g_z_xz_xy_xz, g_z_xz_xy_yy, g_z_xz_xy_yz, g_z_xz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xz_x_xx[i] = 4.0 * g_z_xz_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_x_xy[i] = 4.0 * g_z_xz_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_x_xz[i] = 4.0 * g_z_xz_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_x_yy[i] = 4.0 * g_z_xz_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_x_yz[i] = 4.0 * g_z_xz_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_x_zz[i] = 4.0 * g_z_xz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_z_0_y_0_0_xz_y_xx, g_z_0_y_0_0_xz_y_xy, g_z_0_y_0_0_xz_y_xz, g_z_0_y_0_0_xz_y_yy, g_z_0_y_0_0_xz_y_yz, g_z_0_y_0_0_xz_y_zz, g_z_xz_0_xx, g_z_xz_0_xy, g_z_xz_0_xz, g_z_xz_0_yy, g_z_xz_0_yz, g_z_xz_0_zz, g_z_xz_yy_xx, g_z_xz_yy_xy, g_z_xz_yy_xz, g_z_xz_yy_yy, g_z_xz_yy_yz, g_z_xz_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xz_y_xx[i] = -2.0 * g_z_xz_0_xx[i] * a_exp + 4.0 * g_z_xz_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_y_xy[i] = -2.0 * g_z_xz_0_xy[i] * a_exp + 4.0 * g_z_xz_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_y_xz[i] = -2.0 * g_z_xz_0_xz[i] * a_exp + 4.0 * g_z_xz_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_y_yy[i] = -2.0 * g_z_xz_0_yy[i] * a_exp + 4.0 * g_z_xz_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_y_yz[i] = -2.0 * g_z_xz_0_yz[i] * a_exp + 4.0 * g_z_xz_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_y_zz[i] = -2.0 * g_z_xz_0_zz[i] * a_exp + 4.0 * g_z_xz_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_z_0_y_0_0_xz_z_xx, g_z_0_y_0_0_xz_z_xy, g_z_0_y_0_0_xz_z_xz, g_z_0_y_0_0_xz_z_yy, g_z_0_y_0_0_xz_z_yz, g_z_0_y_0_0_xz_z_zz, g_z_xz_yz_xx, g_z_xz_yz_xy, g_z_xz_yz_xz, g_z_xz_yz_yy, g_z_xz_yz_yz, g_z_xz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xz_z_xx[i] = 4.0 * g_z_xz_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_z_xy[i] = 4.0 * g_z_xz_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_z_xz[i] = 4.0 * g_z_xz_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_z_yy[i] = 4.0 * g_z_xz_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_z_yz[i] = 4.0 * g_z_xz_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_z_zz[i] = 4.0 * g_z_xz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_z_0_y_0_0_yy_x_xx, g_z_0_y_0_0_yy_x_xy, g_z_0_y_0_0_yy_x_xz, g_z_0_y_0_0_yy_x_yy, g_z_0_y_0_0_yy_x_yz, g_z_0_y_0_0_yy_x_zz, g_z_yy_xy_xx, g_z_yy_xy_xy, g_z_yy_xy_xz, g_z_yy_xy_yy, g_z_yy_xy_yz, g_z_yy_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yy_x_xx[i] = 4.0 * g_z_yy_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_x_xy[i] = 4.0 * g_z_yy_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_x_xz[i] = 4.0 * g_z_yy_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_x_yy[i] = 4.0 * g_z_yy_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_x_yz[i] = 4.0 * g_z_yy_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_x_zz[i] = 4.0 * g_z_yy_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_z_0_y_0_0_yy_y_xx, g_z_0_y_0_0_yy_y_xy, g_z_0_y_0_0_yy_y_xz, g_z_0_y_0_0_yy_y_yy, g_z_0_y_0_0_yy_y_yz, g_z_0_y_0_0_yy_y_zz, g_z_yy_0_xx, g_z_yy_0_xy, g_z_yy_0_xz, g_z_yy_0_yy, g_z_yy_0_yz, g_z_yy_0_zz, g_z_yy_yy_xx, g_z_yy_yy_xy, g_z_yy_yy_xz, g_z_yy_yy_yy, g_z_yy_yy_yz, g_z_yy_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yy_y_xx[i] = -2.0 * g_z_yy_0_xx[i] * a_exp + 4.0 * g_z_yy_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_y_xy[i] = -2.0 * g_z_yy_0_xy[i] * a_exp + 4.0 * g_z_yy_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_y_xz[i] = -2.0 * g_z_yy_0_xz[i] * a_exp + 4.0 * g_z_yy_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_y_yy[i] = -2.0 * g_z_yy_0_yy[i] * a_exp + 4.0 * g_z_yy_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_y_yz[i] = -2.0 * g_z_yy_0_yz[i] * a_exp + 4.0 * g_z_yy_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_y_zz[i] = -2.0 * g_z_yy_0_zz[i] * a_exp + 4.0 * g_z_yy_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_z_0_y_0_0_yy_z_xx, g_z_0_y_0_0_yy_z_xy, g_z_0_y_0_0_yy_z_xz, g_z_0_y_0_0_yy_z_yy, g_z_0_y_0_0_yy_z_yz, g_z_0_y_0_0_yy_z_zz, g_z_yy_yz_xx, g_z_yy_yz_xy, g_z_yy_yz_xz, g_z_yy_yz_yy, g_z_yy_yz_yz, g_z_yy_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yy_z_xx[i] = 4.0 * g_z_yy_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_z_xy[i] = 4.0 * g_z_yy_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_z_xz[i] = 4.0 * g_z_yy_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_z_yy[i] = 4.0 * g_z_yy_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_z_yz[i] = 4.0 * g_z_yy_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_z_zz[i] = 4.0 * g_z_yy_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_z_0_y_0_0_yz_x_xx, g_z_0_y_0_0_yz_x_xy, g_z_0_y_0_0_yz_x_xz, g_z_0_y_0_0_yz_x_yy, g_z_0_y_0_0_yz_x_yz, g_z_0_y_0_0_yz_x_zz, g_z_yz_xy_xx, g_z_yz_xy_xy, g_z_yz_xy_xz, g_z_yz_xy_yy, g_z_yz_xy_yz, g_z_yz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yz_x_xx[i] = 4.0 * g_z_yz_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_x_xy[i] = 4.0 * g_z_yz_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_x_xz[i] = 4.0 * g_z_yz_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_x_yy[i] = 4.0 * g_z_yz_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_x_yz[i] = 4.0 * g_z_yz_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_x_zz[i] = 4.0 * g_z_yz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_z_0_y_0_0_yz_y_xx, g_z_0_y_0_0_yz_y_xy, g_z_0_y_0_0_yz_y_xz, g_z_0_y_0_0_yz_y_yy, g_z_0_y_0_0_yz_y_yz, g_z_0_y_0_0_yz_y_zz, g_z_yz_0_xx, g_z_yz_0_xy, g_z_yz_0_xz, g_z_yz_0_yy, g_z_yz_0_yz, g_z_yz_0_zz, g_z_yz_yy_xx, g_z_yz_yy_xy, g_z_yz_yy_xz, g_z_yz_yy_yy, g_z_yz_yy_yz, g_z_yz_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yz_y_xx[i] = -2.0 * g_z_yz_0_xx[i] * a_exp + 4.0 * g_z_yz_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_y_xy[i] = -2.0 * g_z_yz_0_xy[i] * a_exp + 4.0 * g_z_yz_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_y_xz[i] = -2.0 * g_z_yz_0_xz[i] * a_exp + 4.0 * g_z_yz_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_y_yy[i] = -2.0 * g_z_yz_0_yy[i] * a_exp + 4.0 * g_z_yz_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_y_yz[i] = -2.0 * g_z_yz_0_yz[i] * a_exp + 4.0 * g_z_yz_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_y_zz[i] = -2.0 * g_z_yz_0_zz[i] * a_exp + 4.0 * g_z_yz_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_z_0_y_0_0_yz_z_xx, g_z_0_y_0_0_yz_z_xy, g_z_0_y_0_0_yz_z_xz, g_z_0_y_0_0_yz_z_yy, g_z_0_y_0_0_yz_z_yz, g_z_0_y_0_0_yz_z_zz, g_z_yz_yz_xx, g_z_yz_yz_xy, g_z_yz_yz_xz, g_z_yz_yz_yy, g_z_yz_yz_yz, g_z_yz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yz_z_xx[i] = 4.0 * g_z_yz_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_z_xy[i] = 4.0 * g_z_yz_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_z_xz[i] = 4.0 * g_z_yz_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_z_yy[i] = 4.0 * g_z_yz_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_z_yz[i] = 4.0 * g_z_yz_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_z_zz[i] = 4.0 * g_z_yz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_z_0_y_0_0_zz_x_xx, g_z_0_y_0_0_zz_x_xy, g_z_0_y_0_0_zz_x_xz, g_z_0_y_0_0_zz_x_yy, g_z_0_y_0_0_zz_x_yz, g_z_0_y_0_0_zz_x_zz, g_z_zz_xy_xx, g_z_zz_xy_xy, g_z_zz_xy_xz, g_z_zz_xy_yy, g_z_zz_xy_yz, g_z_zz_xy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_zz_x_xx[i] = 4.0 * g_z_zz_xy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_x_xy[i] = 4.0 * g_z_zz_xy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_x_xz[i] = 4.0 * g_z_zz_xy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_x_yy[i] = 4.0 * g_z_zz_xy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_x_yz[i] = 4.0 * g_z_zz_xy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_x_zz[i] = 4.0 * g_z_zz_xy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_z_0_y_0_0_zz_y_xx, g_z_0_y_0_0_zz_y_xy, g_z_0_y_0_0_zz_y_xz, g_z_0_y_0_0_zz_y_yy, g_z_0_y_0_0_zz_y_yz, g_z_0_y_0_0_zz_y_zz, g_z_zz_0_xx, g_z_zz_0_xy, g_z_zz_0_xz, g_z_zz_0_yy, g_z_zz_0_yz, g_z_zz_0_zz, g_z_zz_yy_xx, g_z_zz_yy_xy, g_z_zz_yy_xz, g_z_zz_yy_yy, g_z_zz_yy_yz, g_z_zz_yy_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_zz_y_xx[i] = -2.0 * g_z_zz_0_xx[i] * a_exp + 4.0 * g_z_zz_yy_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_y_xy[i] = -2.0 * g_z_zz_0_xy[i] * a_exp + 4.0 * g_z_zz_yy_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_y_xz[i] = -2.0 * g_z_zz_0_xz[i] * a_exp + 4.0 * g_z_zz_yy_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_y_yy[i] = -2.0 * g_z_zz_0_yy[i] * a_exp + 4.0 * g_z_zz_yy_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_y_yz[i] = -2.0 * g_z_zz_0_yz[i] * a_exp + 4.0 * g_z_zz_yy_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_y_zz[i] = -2.0 * g_z_zz_0_zz[i] * a_exp + 4.0 * g_z_zz_yy_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_z_0_y_0_0_zz_z_xx, g_z_0_y_0_0_zz_z_xy, g_z_0_y_0_0_zz_z_xz, g_z_0_y_0_0_zz_z_yy, g_z_0_y_0_0_zz_z_yz, g_z_0_y_0_0_zz_z_zz, g_z_zz_yz_xx, g_z_zz_yz_xy, g_z_zz_yz_xz, g_z_zz_yz_yy, g_z_zz_yz_yz, g_z_zz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_zz_z_xx[i] = 4.0 * g_z_zz_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_z_xy[i] = 4.0 * g_z_zz_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_z_xz[i] = 4.0 * g_z_zz_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_z_yy[i] = 4.0 * g_z_zz_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_z_yz[i] = 4.0 * g_z_zz_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_z_zz[i] = 4.0 * g_z_zz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_z_0_z_0_0_xx_x_xx, g_z_0_z_0_0_xx_x_xy, g_z_0_z_0_0_xx_x_xz, g_z_0_z_0_0_xx_x_yy, g_z_0_z_0_0_xx_x_yz, g_z_0_z_0_0_xx_x_zz, g_z_xx_xz_xx, g_z_xx_xz_xy, g_z_xx_xz_xz, g_z_xx_xz_yy, g_z_xx_xz_yz, g_z_xx_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xx_x_xx[i] = 4.0 * g_z_xx_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_x_xy[i] = 4.0 * g_z_xx_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_x_xz[i] = 4.0 * g_z_xx_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_x_yy[i] = 4.0 * g_z_xx_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_x_yz[i] = 4.0 * g_z_xx_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_x_zz[i] = 4.0 * g_z_xx_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_z_0_z_0_0_xx_y_xx, g_z_0_z_0_0_xx_y_xy, g_z_0_z_0_0_xx_y_xz, g_z_0_z_0_0_xx_y_yy, g_z_0_z_0_0_xx_y_yz, g_z_0_z_0_0_xx_y_zz, g_z_xx_yz_xx, g_z_xx_yz_xy, g_z_xx_yz_xz, g_z_xx_yz_yy, g_z_xx_yz_yz, g_z_xx_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xx_y_xx[i] = 4.0 * g_z_xx_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_y_xy[i] = 4.0 * g_z_xx_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_y_xz[i] = 4.0 * g_z_xx_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_y_yy[i] = 4.0 * g_z_xx_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_y_yz[i] = 4.0 * g_z_xx_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_y_zz[i] = 4.0 * g_z_xx_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_z_0_z_0_0_xx_z_xx, g_z_0_z_0_0_xx_z_xy, g_z_0_z_0_0_xx_z_xz, g_z_0_z_0_0_xx_z_yy, g_z_0_z_0_0_xx_z_yz, g_z_0_z_0_0_xx_z_zz, g_z_xx_0_xx, g_z_xx_0_xy, g_z_xx_0_xz, g_z_xx_0_yy, g_z_xx_0_yz, g_z_xx_0_zz, g_z_xx_zz_xx, g_z_xx_zz_xy, g_z_xx_zz_xz, g_z_xx_zz_yy, g_z_xx_zz_yz, g_z_xx_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xx_z_xx[i] = -2.0 * g_z_xx_0_xx[i] * a_exp + 4.0 * g_z_xx_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_z_xy[i] = -2.0 * g_z_xx_0_xy[i] * a_exp + 4.0 * g_z_xx_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_z_xz[i] = -2.0 * g_z_xx_0_xz[i] * a_exp + 4.0 * g_z_xx_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_z_yy[i] = -2.0 * g_z_xx_0_yy[i] * a_exp + 4.0 * g_z_xx_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_z_yz[i] = -2.0 * g_z_xx_0_yz[i] * a_exp + 4.0 * g_z_xx_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_z_zz[i] = -2.0 * g_z_xx_0_zz[i] * a_exp + 4.0 * g_z_xx_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_z_0_z_0_0_xy_x_xx, g_z_0_z_0_0_xy_x_xy, g_z_0_z_0_0_xy_x_xz, g_z_0_z_0_0_xy_x_yy, g_z_0_z_0_0_xy_x_yz, g_z_0_z_0_0_xy_x_zz, g_z_xy_xz_xx, g_z_xy_xz_xy, g_z_xy_xz_xz, g_z_xy_xz_yy, g_z_xy_xz_yz, g_z_xy_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xy_x_xx[i] = 4.0 * g_z_xy_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_x_xy[i] = 4.0 * g_z_xy_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_x_xz[i] = 4.0 * g_z_xy_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_x_yy[i] = 4.0 * g_z_xy_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_x_yz[i] = 4.0 * g_z_xy_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_x_zz[i] = 4.0 * g_z_xy_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_z_0_z_0_0_xy_y_xx, g_z_0_z_0_0_xy_y_xy, g_z_0_z_0_0_xy_y_xz, g_z_0_z_0_0_xy_y_yy, g_z_0_z_0_0_xy_y_yz, g_z_0_z_0_0_xy_y_zz, g_z_xy_yz_xx, g_z_xy_yz_xy, g_z_xy_yz_xz, g_z_xy_yz_yy, g_z_xy_yz_yz, g_z_xy_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xy_y_xx[i] = 4.0 * g_z_xy_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_y_xy[i] = 4.0 * g_z_xy_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_y_xz[i] = 4.0 * g_z_xy_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_y_yy[i] = 4.0 * g_z_xy_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_y_yz[i] = 4.0 * g_z_xy_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_y_zz[i] = 4.0 * g_z_xy_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_z_0_z_0_0_xy_z_xx, g_z_0_z_0_0_xy_z_xy, g_z_0_z_0_0_xy_z_xz, g_z_0_z_0_0_xy_z_yy, g_z_0_z_0_0_xy_z_yz, g_z_0_z_0_0_xy_z_zz, g_z_xy_0_xx, g_z_xy_0_xy, g_z_xy_0_xz, g_z_xy_0_yy, g_z_xy_0_yz, g_z_xy_0_zz, g_z_xy_zz_xx, g_z_xy_zz_xy, g_z_xy_zz_xz, g_z_xy_zz_yy, g_z_xy_zz_yz, g_z_xy_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xy_z_xx[i] = -2.0 * g_z_xy_0_xx[i] * a_exp + 4.0 * g_z_xy_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_z_xy[i] = -2.0 * g_z_xy_0_xy[i] * a_exp + 4.0 * g_z_xy_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_z_xz[i] = -2.0 * g_z_xy_0_xz[i] * a_exp + 4.0 * g_z_xy_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_z_yy[i] = -2.0 * g_z_xy_0_yy[i] * a_exp + 4.0 * g_z_xy_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_z_yz[i] = -2.0 * g_z_xy_0_yz[i] * a_exp + 4.0 * g_z_xy_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_z_zz[i] = -2.0 * g_z_xy_0_zz[i] * a_exp + 4.0 * g_z_xy_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_z_0_z_0_0_xz_x_xx, g_z_0_z_0_0_xz_x_xy, g_z_0_z_0_0_xz_x_xz, g_z_0_z_0_0_xz_x_yy, g_z_0_z_0_0_xz_x_yz, g_z_0_z_0_0_xz_x_zz, g_z_xz_xz_xx, g_z_xz_xz_xy, g_z_xz_xz_xz, g_z_xz_xz_yy, g_z_xz_xz_yz, g_z_xz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xz_x_xx[i] = 4.0 * g_z_xz_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_x_xy[i] = 4.0 * g_z_xz_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_x_xz[i] = 4.0 * g_z_xz_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_x_yy[i] = 4.0 * g_z_xz_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_x_yz[i] = 4.0 * g_z_xz_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_x_zz[i] = 4.0 * g_z_xz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_z_0_z_0_0_xz_y_xx, g_z_0_z_0_0_xz_y_xy, g_z_0_z_0_0_xz_y_xz, g_z_0_z_0_0_xz_y_yy, g_z_0_z_0_0_xz_y_yz, g_z_0_z_0_0_xz_y_zz, g_z_xz_yz_xx, g_z_xz_yz_xy, g_z_xz_yz_xz, g_z_xz_yz_yy, g_z_xz_yz_yz, g_z_xz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xz_y_xx[i] = 4.0 * g_z_xz_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_y_xy[i] = 4.0 * g_z_xz_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_y_xz[i] = 4.0 * g_z_xz_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_y_yy[i] = 4.0 * g_z_xz_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_y_yz[i] = 4.0 * g_z_xz_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_y_zz[i] = 4.0 * g_z_xz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_z_0_z_0_0_xz_z_xx, g_z_0_z_0_0_xz_z_xy, g_z_0_z_0_0_xz_z_xz, g_z_0_z_0_0_xz_z_yy, g_z_0_z_0_0_xz_z_yz, g_z_0_z_0_0_xz_z_zz, g_z_xz_0_xx, g_z_xz_0_xy, g_z_xz_0_xz, g_z_xz_0_yy, g_z_xz_0_yz, g_z_xz_0_zz, g_z_xz_zz_xx, g_z_xz_zz_xy, g_z_xz_zz_xz, g_z_xz_zz_yy, g_z_xz_zz_yz, g_z_xz_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xz_z_xx[i] = -2.0 * g_z_xz_0_xx[i] * a_exp + 4.0 * g_z_xz_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_z_xy[i] = -2.0 * g_z_xz_0_xy[i] * a_exp + 4.0 * g_z_xz_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_z_xz[i] = -2.0 * g_z_xz_0_xz[i] * a_exp + 4.0 * g_z_xz_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_z_yy[i] = -2.0 * g_z_xz_0_yy[i] * a_exp + 4.0 * g_z_xz_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_z_yz[i] = -2.0 * g_z_xz_0_yz[i] * a_exp + 4.0 * g_z_xz_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_z_zz[i] = -2.0 * g_z_xz_0_zz[i] * a_exp + 4.0 * g_z_xz_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_z_0_z_0_0_yy_x_xx, g_z_0_z_0_0_yy_x_xy, g_z_0_z_0_0_yy_x_xz, g_z_0_z_0_0_yy_x_yy, g_z_0_z_0_0_yy_x_yz, g_z_0_z_0_0_yy_x_zz, g_z_yy_xz_xx, g_z_yy_xz_xy, g_z_yy_xz_xz, g_z_yy_xz_yy, g_z_yy_xz_yz, g_z_yy_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yy_x_xx[i] = 4.0 * g_z_yy_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_x_xy[i] = 4.0 * g_z_yy_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_x_xz[i] = 4.0 * g_z_yy_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_x_yy[i] = 4.0 * g_z_yy_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_x_yz[i] = 4.0 * g_z_yy_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_x_zz[i] = 4.0 * g_z_yy_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_z_0_z_0_0_yy_y_xx, g_z_0_z_0_0_yy_y_xy, g_z_0_z_0_0_yy_y_xz, g_z_0_z_0_0_yy_y_yy, g_z_0_z_0_0_yy_y_yz, g_z_0_z_0_0_yy_y_zz, g_z_yy_yz_xx, g_z_yy_yz_xy, g_z_yy_yz_xz, g_z_yy_yz_yy, g_z_yy_yz_yz, g_z_yy_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yy_y_xx[i] = 4.0 * g_z_yy_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_y_xy[i] = 4.0 * g_z_yy_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_y_xz[i] = 4.0 * g_z_yy_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_y_yy[i] = 4.0 * g_z_yy_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_y_yz[i] = 4.0 * g_z_yy_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_y_zz[i] = 4.0 * g_z_yy_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_z_0_z_0_0_yy_z_xx, g_z_0_z_0_0_yy_z_xy, g_z_0_z_0_0_yy_z_xz, g_z_0_z_0_0_yy_z_yy, g_z_0_z_0_0_yy_z_yz, g_z_0_z_0_0_yy_z_zz, g_z_yy_0_xx, g_z_yy_0_xy, g_z_yy_0_xz, g_z_yy_0_yy, g_z_yy_0_yz, g_z_yy_0_zz, g_z_yy_zz_xx, g_z_yy_zz_xy, g_z_yy_zz_xz, g_z_yy_zz_yy, g_z_yy_zz_yz, g_z_yy_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yy_z_xx[i] = -2.0 * g_z_yy_0_xx[i] * a_exp + 4.0 * g_z_yy_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_z_xy[i] = -2.0 * g_z_yy_0_xy[i] * a_exp + 4.0 * g_z_yy_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_z_xz[i] = -2.0 * g_z_yy_0_xz[i] * a_exp + 4.0 * g_z_yy_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_z_yy[i] = -2.0 * g_z_yy_0_yy[i] * a_exp + 4.0 * g_z_yy_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_z_yz[i] = -2.0 * g_z_yy_0_yz[i] * a_exp + 4.0 * g_z_yy_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_z_zz[i] = -2.0 * g_z_yy_0_zz[i] * a_exp + 4.0 * g_z_yy_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_z_0_z_0_0_yz_x_xx, g_z_0_z_0_0_yz_x_xy, g_z_0_z_0_0_yz_x_xz, g_z_0_z_0_0_yz_x_yy, g_z_0_z_0_0_yz_x_yz, g_z_0_z_0_0_yz_x_zz, g_z_yz_xz_xx, g_z_yz_xz_xy, g_z_yz_xz_xz, g_z_yz_xz_yy, g_z_yz_xz_yz, g_z_yz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yz_x_xx[i] = 4.0 * g_z_yz_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_x_xy[i] = 4.0 * g_z_yz_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_x_xz[i] = 4.0 * g_z_yz_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_x_yy[i] = 4.0 * g_z_yz_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_x_yz[i] = 4.0 * g_z_yz_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_x_zz[i] = 4.0 * g_z_yz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_z_0_z_0_0_yz_y_xx, g_z_0_z_0_0_yz_y_xy, g_z_0_z_0_0_yz_y_xz, g_z_0_z_0_0_yz_y_yy, g_z_0_z_0_0_yz_y_yz, g_z_0_z_0_0_yz_y_zz, g_z_yz_yz_xx, g_z_yz_yz_xy, g_z_yz_yz_xz, g_z_yz_yz_yy, g_z_yz_yz_yz, g_z_yz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yz_y_xx[i] = 4.0 * g_z_yz_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_y_xy[i] = 4.0 * g_z_yz_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_y_xz[i] = 4.0 * g_z_yz_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_y_yy[i] = 4.0 * g_z_yz_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_y_yz[i] = 4.0 * g_z_yz_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_y_zz[i] = 4.0 * g_z_yz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_z_0_z_0_0_yz_z_xx, g_z_0_z_0_0_yz_z_xy, g_z_0_z_0_0_yz_z_xz, g_z_0_z_0_0_yz_z_yy, g_z_0_z_0_0_yz_z_yz, g_z_0_z_0_0_yz_z_zz, g_z_yz_0_xx, g_z_yz_0_xy, g_z_yz_0_xz, g_z_yz_0_yy, g_z_yz_0_yz, g_z_yz_0_zz, g_z_yz_zz_xx, g_z_yz_zz_xy, g_z_yz_zz_xz, g_z_yz_zz_yy, g_z_yz_zz_yz, g_z_yz_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yz_z_xx[i] = -2.0 * g_z_yz_0_xx[i] * a_exp + 4.0 * g_z_yz_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_z_xy[i] = -2.0 * g_z_yz_0_xy[i] * a_exp + 4.0 * g_z_yz_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_z_xz[i] = -2.0 * g_z_yz_0_xz[i] * a_exp + 4.0 * g_z_yz_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_z_yy[i] = -2.0 * g_z_yz_0_yy[i] * a_exp + 4.0 * g_z_yz_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_z_yz[i] = -2.0 * g_z_yz_0_yz[i] * a_exp + 4.0 * g_z_yz_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_z_zz[i] = -2.0 * g_z_yz_0_zz[i] * a_exp + 4.0 * g_z_yz_zz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_z_0_z_0_0_zz_x_xx, g_z_0_z_0_0_zz_x_xy, g_z_0_z_0_0_zz_x_xz, g_z_0_z_0_0_zz_x_yy, g_z_0_z_0_0_zz_x_yz, g_z_0_z_0_0_zz_x_zz, g_z_zz_xz_xx, g_z_zz_xz_xy, g_z_zz_xz_xz, g_z_zz_xz_yy, g_z_zz_xz_yz, g_z_zz_xz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_zz_x_xx[i] = 4.0 * g_z_zz_xz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_x_xy[i] = 4.0 * g_z_zz_xz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_x_xz[i] = 4.0 * g_z_zz_xz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_x_yy[i] = 4.0 * g_z_zz_xz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_x_yz[i] = 4.0 * g_z_zz_xz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_x_zz[i] = 4.0 * g_z_zz_xz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_z_0_z_0_0_zz_y_xx, g_z_0_z_0_0_zz_y_xy, g_z_0_z_0_0_zz_y_xz, g_z_0_z_0_0_zz_y_yy, g_z_0_z_0_0_zz_y_yz, g_z_0_z_0_0_zz_y_zz, g_z_zz_yz_xx, g_z_zz_yz_xy, g_z_zz_yz_xz, g_z_zz_yz_yy, g_z_zz_yz_yz, g_z_zz_yz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_zz_y_xx[i] = 4.0 * g_z_zz_yz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_y_xy[i] = 4.0 * g_z_zz_yz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_y_xz[i] = 4.0 * g_z_zz_yz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_y_yy[i] = 4.0 * g_z_zz_yz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_y_yz[i] = 4.0 * g_z_zz_yz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_y_zz[i] = 4.0 * g_z_zz_yz_zz[i] * a_exp * c_exps[i];
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_z_0_z_0_0_zz_z_xx, g_z_0_z_0_0_zz_z_xy, g_z_0_z_0_0_zz_z_xz, g_z_0_z_0_0_zz_z_yy, g_z_0_z_0_0_zz_z_yz, g_z_0_z_0_0_zz_z_zz, g_z_zz_0_xx, g_z_zz_0_xy, g_z_zz_0_xz, g_z_zz_0_yy, g_z_zz_0_yz, g_z_zz_0_zz, g_z_zz_zz_xx, g_z_zz_zz_xy, g_z_zz_zz_xz, g_z_zz_zz_yy, g_z_zz_zz_yz, g_z_zz_zz_zz, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_zz_z_xx[i] = -2.0 * g_z_zz_0_xx[i] * a_exp + 4.0 * g_z_zz_zz_xx[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_z_xy[i] = -2.0 * g_z_zz_0_xy[i] * a_exp + 4.0 * g_z_zz_zz_xy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_z_xz[i] = -2.0 * g_z_zz_0_xz[i] * a_exp + 4.0 * g_z_zz_zz_xz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_z_yy[i] = -2.0 * g_z_zz_0_yy[i] * a_exp + 4.0 * g_z_zz_zz_yy[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_z_yz[i] = -2.0 * g_z_zz_0_yz[i] * a_exp + 4.0 * g_z_zz_zz_yz[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_z_zz[i] = -2.0 * g_z_zz_0_zz[i] * a_exp + 4.0 * g_z_zz_zz_zz[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace

