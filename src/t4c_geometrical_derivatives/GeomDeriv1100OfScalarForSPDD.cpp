#include "GeomDeriv1100OfScalarForSPDD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_spdd_0(CSimdArray<double>& buffer_1100_spdd,
                     const CSimdArray<double>& buffer_psdd,
                     const CSimdArray<double>& buffer_pddd,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_spdd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_psdd

    auto g_x_0_xx_xx = buffer_psdd[0];

    auto g_x_0_xx_xy = buffer_psdd[1];

    auto g_x_0_xx_xz = buffer_psdd[2];

    auto g_x_0_xx_yy = buffer_psdd[3];

    auto g_x_0_xx_yz = buffer_psdd[4];

    auto g_x_0_xx_zz = buffer_psdd[5];

    auto g_x_0_xy_xx = buffer_psdd[6];

    auto g_x_0_xy_xy = buffer_psdd[7];

    auto g_x_0_xy_xz = buffer_psdd[8];

    auto g_x_0_xy_yy = buffer_psdd[9];

    auto g_x_0_xy_yz = buffer_psdd[10];

    auto g_x_0_xy_zz = buffer_psdd[11];

    auto g_x_0_xz_xx = buffer_psdd[12];

    auto g_x_0_xz_xy = buffer_psdd[13];

    auto g_x_0_xz_xz = buffer_psdd[14];

    auto g_x_0_xz_yy = buffer_psdd[15];

    auto g_x_0_xz_yz = buffer_psdd[16];

    auto g_x_0_xz_zz = buffer_psdd[17];

    auto g_x_0_yy_xx = buffer_psdd[18];

    auto g_x_0_yy_xy = buffer_psdd[19];

    auto g_x_0_yy_xz = buffer_psdd[20];

    auto g_x_0_yy_yy = buffer_psdd[21];

    auto g_x_0_yy_yz = buffer_psdd[22];

    auto g_x_0_yy_zz = buffer_psdd[23];

    auto g_x_0_yz_xx = buffer_psdd[24];

    auto g_x_0_yz_xy = buffer_psdd[25];

    auto g_x_0_yz_xz = buffer_psdd[26];

    auto g_x_0_yz_yy = buffer_psdd[27];

    auto g_x_0_yz_yz = buffer_psdd[28];

    auto g_x_0_yz_zz = buffer_psdd[29];

    auto g_x_0_zz_xx = buffer_psdd[30];

    auto g_x_0_zz_xy = buffer_psdd[31];

    auto g_x_0_zz_xz = buffer_psdd[32];

    auto g_x_0_zz_yy = buffer_psdd[33];

    auto g_x_0_zz_yz = buffer_psdd[34];

    auto g_x_0_zz_zz = buffer_psdd[35];

    auto g_y_0_xx_xx = buffer_psdd[36];

    auto g_y_0_xx_xy = buffer_psdd[37];

    auto g_y_0_xx_xz = buffer_psdd[38];

    auto g_y_0_xx_yy = buffer_psdd[39];

    auto g_y_0_xx_yz = buffer_psdd[40];

    auto g_y_0_xx_zz = buffer_psdd[41];

    auto g_y_0_xy_xx = buffer_psdd[42];

    auto g_y_0_xy_xy = buffer_psdd[43];

    auto g_y_0_xy_xz = buffer_psdd[44];

    auto g_y_0_xy_yy = buffer_psdd[45];

    auto g_y_0_xy_yz = buffer_psdd[46];

    auto g_y_0_xy_zz = buffer_psdd[47];

    auto g_y_0_xz_xx = buffer_psdd[48];

    auto g_y_0_xz_xy = buffer_psdd[49];

    auto g_y_0_xz_xz = buffer_psdd[50];

    auto g_y_0_xz_yy = buffer_psdd[51];

    auto g_y_0_xz_yz = buffer_psdd[52];

    auto g_y_0_xz_zz = buffer_psdd[53];

    auto g_y_0_yy_xx = buffer_psdd[54];

    auto g_y_0_yy_xy = buffer_psdd[55];

    auto g_y_0_yy_xz = buffer_psdd[56];

    auto g_y_0_yy_yy = buffer_psdd[57];

    auto g_y_0_yy_yz = buffer_psdd[58];

    auto g_y_0_yy_zz = buffer_psdd[59];

    auto g_y_0_yz_xx = buffer_psdd[60];

    auto g_y_0_yz_xy = buffer_psdd[61];

    auto g_y_0_yz_xz = buffer_psdd[62];

    auto g_y_0_yz_yy = buffer_psdd[63];

    auto g_y_0_yz_yz = buffer_psdd[64];

    auto g_y_0_yz_zz = buffer_psdd[65];

    auto g_y_0_zz_xx = buffer_psdd[66];

    auto g_y_0_zz_xy = buffer_psdd[67];

    auto g_y_0_zz_xz = buffer_psdd[68];

    auto g_y_0_zz_yy = buffer_psdd[69];

    auto g_y_0_zz_yz = buffer_psdd[70];

    auto g_y_0_zz_zz = buffer_psdd[71];

    auto g_z_0_xx_xx = buffer_psdd[72];

    auto g_z_0_xx_xy = buffer_psdd[73];

    auto g_z_0_xx_xz = buffer_psdd[74];

    auto g_z_0_xx_yy = buffer_psdd[75];

    auto g_z_0_xx_yz = buffer_psdd[76];

    auto g_z_0_xx_zz = buffer_psdd[77];

    auto g_z_0_xy_xx = buffer_psdd[78];

    auto g_z_0_xy_xy = buffer_psdd[79];

    auto g_z_0_xy_xz = buffer_psdd[80];

    auto g_z_0_xy_yy = buffer_psdd[81];

    auto g_z_0_xy_yz = buffer_psdd[82];

    auto g_z_0_xy_zz = buffer_psdd[83];

    auto g_z_0_xz_xx = buffer_psdd[84];

    auto g_z_0_xz_xy = buffer_psdd[85];

    auto g_z_0_xz_xz = buffer_psdd[86];

    auto g_z_0_xz_yy = buffer_psdd[87];

    auto g_z_0_xz_yz = buffer_psdd[88];

    auto g_z_0_xz_zz = buffer_psdd[89];

    auto g_z_0_yy_xx = buffer_psdd[90];

    auto g_z_0_yy_xy = buffer_psdd[91];

    auto g_z_0_yy_xz = buffer_psdd[92];

    auto g_z_0_yy_yy = buffer_psdd[93];

    auto g_z_0_yy_yz = buffer_psdd[94];

    auto g_z_0_yy_zz = buffer_psdd[95];

    auto g_z_0_yz_xx = buffer_psdd[96];

    auto g_z_0_yz_xy = buffer_psdd[97];

    auto g_z_0_yz_xz = buffer_psdd[98];

    auto g_z_0_yz_yy = buffer_psdd[99];

    auto g_z_0_yz_yz = buffer_psdd[100];

    auto g_z_0_yz_zz = buffer_psdd[101];

    auto g_z_0_zz_xx = buffer_psdd[102];

    auto g_z_0_zz_xy = buffer_psdd[103];

    auto g_z_0_zz_xz = buffer_psdd[104];

    auto g_z_0_zz_yy = buffer_psdd[105];

    auto g_z_0_zz_yz = buffer_psdd[106];

    auto g_z_0_zz_zz = buffer_psdd[107];

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

    /// Set up components of integrals buffer : buffer_1100_spdd

    auto g_x_x_0_0_0_x_xx_xx = buffer_1100_spdd[0];

    auto g_x_x_0_0_0_x_xx_xy = buffer_1100_spdd[1];

    auto g_x_x_0_0_0_x_xx_xz = buffer_1100_spdd[2];

    auto g_x_x_0_0_0_x_xx_yy = buffer_1100_spdd[3];

    auto g_x_x_0_0_0_x_xx_yz = buffer_1100_spdd[4];

    auto g_x_x_0_0_0_x_xx_zz = buffer_1100_spdd[5];

    auto g_x_x_0_0_0_x_xy_xx = buffer_1100_spdd[6];

    auto g_x_x_0_0_0_x_xy_xy = buffer_1100_spdd[7];

    auto g_x_x_0_0_0_x_xy_xz = buffer_1100_spdd[8];

    auto g_x_x_0_0_0_x_xy_yy = buffer_1100_spdd[9];

    auto g_x_x_0_0_0_x_xy_yz = buffer_1100_spdd[10];

    auto g_x_x_0_0_0_x_xy_zz = buffer_1100_spdd[11];

    auto g_x_x_0_0_0_x_xz_xx = buffer_1100_spdd[12];

    auto g_x_x_0_0_0_x_xz_xy = buffer_1100_spdd[13];

    auto g_x_x_0_0_0_x_xz_xz = buffer_1100_spdd[14];

    auto g_x_x_0_0_0_x_xz_yy = buffer_1100_spdd[15];

    auto g_x_x_0_0_0_x_xz_yz = buffer_1100_spdd[16];

    auto g_x_x_0_0_0_x_xz_zz = buffer_1100_spdd[17];

    auto g_x_x_0_0_0_x_yy_xx = buffer_1100_spdd[18];

    auto g_x_x_0_0_0_x_yy_xy = buffer_1100_spdd[19];

    auto g_x_x_0_0_0_x_yy_xz = buffer_1100_spdd[20];

    auto g_x_x_0_0_0_x_yy_yy = buffer_1100_spdd[21];

    auto g_x_x_0_0_0_x_yy_yz = buffer_1100_spdd[22];

    auto g_x_x_0_0_0_x_yy_zz = buffer_1100_spdd[23];

    auto g_x_x_0_0_0_x_yz_xx = buffer_1100_spdd[24];

    auto g_x_x_0_0_0_x_yz_xy = buffer_1100_spdd[25];

    auto g_x_x_0_0_0_x_yz_xz = buffer_1100_spdd[26];

    auto g_x_x_0_0_0_x_yz_yy = buffer_1100_spdd[27];

    auto g_x_x_0_0_0_x_yz_yz = buffer_1100_spdd[28];

    auto g_x_x_0_0_0_x_yz_zz = buffer_1100_spdd[29];

    auto g_x_x_0_0_0_x_zz_xx = buffer_1100_spdd[30];

    auto g_x_x_0_0_0_x_zz_xy = buffer_1100_spdd[31];

    auto g_x_x_0_0_0_x_zz_xz = buffer_1100_spdd[32];

    auto g_x_x_0_0_0_x_zz_yy = buffer_1100_spdd[33];

    auto g_x_x_0_0_0_x_zz_yz = buffer_1100_spdd[34];

    auto g_x_x_0_0_0_x_zz_zz = buffer_1100_spdd[35];

    auto g_x_x_0_0_0_y_xx_xx = buffer_1100_spdd[36];

    auto g_x_x_0_0_0_y_xx_xy = buffer_1100_spdd[37];

    auto g_x_x_0_0_0_y_xx_xz = buffer_1100_spdd[38];

    auto g_x_x_0_0_0_y_xx_yy = buffer_1100_spdd[39];

    auto g_x_x_0_0_0_y_xx_yz = buffer_1100_spdd[40];

    auto g_x_x_0_0_0_y_xx_zz = buffer_1100_spdd[41];

    auto g_x_x_0_0_0_y_xy_xx = buffer_1100_spdd[42];

    auto g_x_x_0_0_0_y_xy_xy = buffer_1100_spdd[43];

    auto g_x_x_0_0_0_y_xy_xz = buffer_1100_spdd[44];

    auto g_x_x_0_0_0_y_xy_yy = buffer_1100_spdd[45];

    auto g_x_x_0_0_0_y_xy_yz = buffer_1100_spdd[46];

    auto g_x_x_0_0_0_y_xy_zz = buffer_1100_spdd[47];

    auto g_x_x_0_0_0_y_xz_xx = buffer_1100_spdd[48];

    auto g_x_x_0_0_0_y_xz_xy = buffer_1100_spdd[49];

    auto g_x_x_0_0_0_y_xz_xz = buffer_1100_spdd[50];

    auto g_x_x_0_0_0_y_xz_yy = buffer_1100_spdd[51];

    auto g_x_x_0_0_0_y_xz_yz = buffer_1100_spdd[52];

    auto g_x_x_0_0_0_y_xz_zz = buffer_1100_spdd[53];

    auto g_x_x_0_0_0_y_yy_xx = buffer_1100_spdd[54];

    auto g_x_x_0_0_0_y_yy_xy = buffer_1100_spdd[55];

    auto g_x_x_0_0_0_y_yy_xz = buffer_1100_spdd[56];

    auto g_x_x_0_0_0_y_yy_yy = buffer_1100_spdd[57];

    auto g_x_x_0_0_0_y_yy_yz = buffer_1100_spdd[58];

    auto g_x_x_0_0_0_y_yy_zz = buffer_1100_spdd[59];

    auto g_x_x_0_0_0_y_yz_xx = buffer_1100_spdd[60];

    auto g_x_x_0_0_0_y_yz_xy = buffer_1100_spdd[61];

    auto g_x_x_0_0_0_y_yz_xz = buffer_1100_spdd[62];

    auto g_x_x_0_0_0_y_yz_yy = buffer_1100_spdd[63];

    auto g_x_x_0_0_0_y_yz_yz = buffer_1100_spdd[64];

    auto g_x_x_0_0_0_y_yz_zz = buffer_1100_spdd[65];

    auto g_x_x_0_0_0_y_zz_xx = buffer_1100_spdd[66];

    auto g_x_x_0_0_0_y_zz_xy = buffer_1100_spdd[67];

    auto g_x_x_0_0_0_y_zz_xz = buffer_1100_spdd[68];

    auto g_x_x_0_0_0_y_zz_yy = buffer_1100_spdd[69];

    auto g_x_x_0_0_0_y_zz_yz = buffer_1100_spdd[70];

    auto g_x_x_0_0_0_y_zz_zz = buffer_1100_spdd[71];

    auto g_x_x_0_0_0_z_xx_xx = buffer_1100_spdd[72];

    auto g_x_x_0_0_0_z_xx_xy = buffer_1100_spdd[73];

    auto g_x_x_0_0_0_z_xx_xz = buffer_1100_spdd[74];

    auto g_x_x_0_0_0_z_xx_yy = buffer_1100_spdd[75];

    auto g_x_x_0_0_0_z_xx_yz = buffer_1100_spdd[76];

    auto g_x_x_0_0_0_z_xx_zz = buffer_1100_spdd[77];

    auto g_x_x_0_0_0_z_xy_xx = buffer_1100_spdd[78];

    auto g_x_x_0_0_0_z_xy_xy = buffer_1100_spdd[79];

    auto g_x_x_0_0_0_z_xy_xz = buffer_1100_spdd[80];

    auto g_x_x_0_0_0_z_xy_yy = buffer_1100_spdd[81];

    auto g_x_x_0_0_0_z_xy_yz = buffer_1100_spdd[82];

    auto g_x_x_0_0_0_z_xy_zz = buffer_1100_spdd[83];

    auto g_x_x_0_0_0_z_xz_xx = buffer_1100_spdd[84];

    auto g_x_x_0_0_0_z_xz_xy = buffer_1100_spdd[85];

    auto g_x_x_0_0_0_z_xz_xz = buffer_1100_spdd[86];

    auto g_x_x_0_0_0_z_xz_yy = buffer_1100_spdd[87];

    auto g_x_x_0_0_0_z_xz_yz = buffer_1100_spdd[88];

    auto g_x_x_0_0_0_z_xz_zz = buffer_1100_spdd[89];

    auto g_x_x_0_0_0_z_yy_xx = buffer_1100_spdd[90];

    auto g_x_x_0_0_0_z_yy_xy = buffer_1100_spdd[91];

    auto g_x_x_0_0_0_z_yy_xz = buffer_1100_spdd[92];

    auto g_x_x_0_0_0_z_yy_yy = buffer_1100_spdd[93];

    auto g_x_x_0_0_0_z_yy_yz = buffer_1100_spdd[94];

    auto g_x_x_0_0_0_z_yy_zz = buffer_1100_spdd[95];

    auto g_x_x_0_0_0_z_yz_xx = buffer_1100_spdd[96];

    auto g_x_x_0_0_0_z_yz_xy = buffer_1100_spdd[97];

    auto g_x_x_0_0_0_z_yz_xz = buffer_1100_spdd[98];

    auto g_x_x_0_0_0_z_yz_yy = buffer_1100_spdd[99];

    auto g_x_x_0_0_0_z_yz_yz = buffer_1100_spdd[100];

    auto g_x_x_0_0_0_z_yz_zz = buffer_1100_spdd[101];

    auto g_x_x_0_0_0_z_zz_xx = buffer_1100_spdd[102];

    auto g_x_x_0_0_0_z_zz_xy = buffer_1100_spdd[103];

    auto g_x_x_0_0_0_z_zz_xz = buffer_1100_spdd[104];

    auto g_x_x_0_0_0_z_zz_yy = buffer_1100_spdd[105];

    auto g_x_x_0_0_0_z_zz_yz = buffer_1100_spdd[106];

    auto g_x_x_0_0_0_z_zz_zz = buffer_1100_spdd[107];

    auto g_x_y_0_0_0_x_xx_xx = buffer_1100_spdd[108];

    auto g_x_y_0_0_0_x_xx_xy = buffer_1100_spdd[109];

    auto g_x_y_0_0_0_x_xx_xz = buffer_1100_spdd[110];

    auto g_x_y_0_0_0_x_xx_yy = buffer_1100_spdd[111];

    auto g_x_y_0_0_0_x_xx_yz = buffer_1100_spdd[112];

    auto g_x_y_0_0_0_x_xx_zz = buffer_1100_spdd[113];

    auto g_x_y_0_0_0_x_xy_xx = buffer_1100_spdd[114];

    auto g_x_y_0_0_0_x_xy_xy = buffer_1100_spdd[115];

    auto g_x_y_0_0_0_x_xy_xz = buffer_1100_spdd[116];

    auto g_x_y_0_0_0_x_xy_yy = buffer_1100_spdd[117];

    auto g_x_y_0_0_0_x_xy_yz = buffer_1100_spdd[118];

    auto g_x_y_0_0_0_x_xy_zz = buffer_1100_spdd[119];

    auto g_x_y_0_0_0_x_xz_xx = buffer_1100_spdd[120];

    auto g_x_y_0_0_0_x_xz_xy = buffer_1100_spdd[121];

    auto g_x_y_0_0_0_x_xz_xz = buffer_1100_spdd[122];

    auto g_x_y_0_0_0_x_xz_yy = buffer_1100_spdd[123];

    auto g_x_y_0_0_0_x_xz_yz = buffer_1100_spdd[124];

    auto g_x_y_0_0_0_x_xz_zz = buffer_1100_spdd[125];

    auto g_x_y_0_0_0_x_yy_xx = buffer_1100_spdd[126];

    auto g_x_y_0_0_0_x_yy_xy = buffer_1100_spdd[127];

    auto g_x_y_0_0_0_x_yy_xz = buffer_1100_spdd[128];

    auto g_x_y_0_0_0_x_yy_yy = buffer_1100_spdd[129];

    auto g_x_y_0_0_0_x_yy_yz = buffer_1100_spdd[130];

    auto g_x_y_0_0_0_x_yy_zz = buffer_1100_spdd[131];

    auto g_x_y_0_0_0_x_yz_xx = buffer_1100_spdd[132];

    auto g_x_y_0_0_0_x_yz_xy = buffer_1100_spdd[133];

    auto g_x_y_0_0_0_x_yz_xz = buffer_1100_spdd[134];

    auto g_x_y_0_0_0_x_yz_yy = buffer_1100_spdd[135];

    auto g_x_y_0_0_0_x_yz_yz = buffer_1100_spdd[136];

    auto g_x_y_0_0_0_x_yz_zz = buffer_1100_spdd[137];

    auto g_x_y_0_0_0_x_zz_xx = buffer_1100_spdd[138];

    auto g_x_y_0_0_0_x_zz_xy = buffer_1100_spdd[139];

    auto g_x_y_0_0_0_x_zz_xz = buffer_1100_spdd[140];

    auto g_x_y_0_0_0_x_zz_yy = buffer_1100_spdd[141];

    auto g_x_y_0_0_0_x_zz_yz = buffer_1100_spdd[142];

    auto g_x_y_0_0_0_x_zz_zz = buffer_1100_spdd[143];

    auto g_x_y_0_0_0_y_xx_xx = buffer_1100_spdd[144];

    auto g_x_y_0_0_0_y_xx_xy = buffer_1100_spdd[145];

    auto g_x_y_0_0_0_y_xx_xz = buffer_1100_spdd[146];

    auto g_x_y_0_0_0_y_xx_yy = buffer_1100_spdd[147];

    auto g_x_y_0_0_0_y_xx_yz = buffer_1100_spdd[148];

    auto g_x_y_0_0_0_y_xx_zz = buffer_1100_spdd[149];

    auto g_x_y_0_0_0_y_xy_xx = buffer_1100_spdd[150];

    auto g_x_y_0_0_0_y_xy_xy = buffer_1100_spdd[151];

    auto g_x_y_0_0_0_y_xy_xz = buffer_1100_spdd[152];

    auto g_x_y_0_0_0_y_xy_yy = buffer_1100_spdd[153];

    auto g_x_y_0_0_0_y_xy_yz = buffer_1100_spdd[154];

    auto g_x_y_0_0_0_y_xy_zz = buffer_1100_spdd[155];

    auto g_x_y_0_0_0_y_xz_xx = buffer_1100_spdd[156];

    auto g_x_y_0_0_0_y_xz_xy = buffer_1100_spdd[157];

    auto g_x_y_0_0_0_y_xz_xz = buffer_1100_spdd[158];

    auto g_x_y_0_0_0_y_xz_yy = buffer_1100_spdd[159];

    auto g_x_y_0_0_0_y_xz_yz = buffer_1100_spdd[160];

    auto g_x_y_0_0_0_y_xz_zz = buffer_1100_spdd[161];

    auto g_x_y_0_0_0_y_yy_xx = buffer_1100_spdd[162];

    auto g_x_y_0_0_0_y_yy_xy = buffer_1100_spdd[163];

    auto g_x_y_0_0_0_y_yy_xz = buffer_1100_spdd[164];

    auto g_x_y_0_0_0_y_yy_yy = buffer_1100_spdd[165];

    auto g_x_y_0_0_0_y_yy_yz = buffer_1100_spdd[166];

    auto g_x_y_0_0_0_y_yy_zz = buffer_1100_spdd[167];

    auto g_x_y_0_0_0_y_yz_xx = buffer_1100_spdd[168];

    auto g_x_y_0_0_0_y_yz_xy = buffer_1100_spdd[169];

    auto g_x_y_0_0_0_y_yz_xz = buffer_1100_spdd[170];

    auto g_x_y_0_0_0_y_yz_yy = buffer_1100_spdd[171];

    auto g_x_y_0_0_0_y_yz_yz = buffer_1100_spdd[172];

    auto g_x_y_0_0_0_y_yz_zz = buffer_1100_spdd[173];

    auto g_x_y_0_0_0_y_zz_xx = buffer_1100_spdd[174];

    auto g_x_y_0_0_0_y_zz_xy = buffer_1100_spdd[175];

    auto g_x_y_0_0_0_y_zz_xz = buffer_1100_spdd[176];

    auto g_x_y_0_0_0_y_zz_yy = buffer_1100_spdd[177];

    auto g_x_y_0_0_0_y_zz_yz = buffer_1100_spdd[178];

    auto g_x_y_0_0_0_y_zz_zz = buffer_1100_spdd[179];

    auto g_x_y_0_0_0_z_xx_xx = buffer_1100_spdd[180];

    auto g_x_y_0_0_0_z_xx_xy = buffer_1100_spdd[181];

    auto g_x_y_0_0_0_z_xx_xz = buffer_1100_spdd[182];

    auto g_x_y_0_0_0_z_xx_yy = buffer_1100_spdd[183];

    auto g_x_y_0_0_0_z_xx_yz = buffer_1100_spdd[184];

    auto g_x_y_0_0_0_z_xx_zz = buffer_1100_spdd[185];

    auto g_x_y_0_0_0_z_xy_xx = buffer_1100_spdd[186];

    auto g_x_y_0_0_0_z_xy_xy = buffer_1100_spdd[187];

    auto g_x_y_0_0_0_z_xy_xz = buffer_1100_spdd[188];

    auto g_x_y_0_0_0_z_xy_yy = buffer_1100_spdd[189];

    auto g_x_y_0_0_0_z_xy_yz = buffer_1100_spdd[190];

    auto g_x_y_0_0_0_z_xy_zz = buffer_1100_spdd[191];

    auto g_x_y_0_0_0_z_xz_xx = buffer_1100_spdd[192];

    auto g_x_y_0_0_0_z_xz_xy = buffer_1100_spdd[193];

    auto g_x_y_0_0_0_z_xz_xz = buffer_1100_spdd[194];

    auto g_x_y_0_0_0_z_xz_yy = buffer_1100_spdd[195];

    auto g_x_y_0_0_0_z_xz_yz = buffer_1100_spdd[196];

    auto g_x_y_0_0_0_z_xz_zz = buffer_1100_spdd[197];

    auto g_x_y_0_0_0_z_yy_xx = buffer_1100_spdd[198];

    auto g_x_y_0_0_0_z_yy_xy = buffer_1100_spdd[199];

    auto g_x_y_0_0_0_z_yy_xz = buffer_1100_spdd[200];

    auto g_x_y_0_0_0_z_yy_yy = buffer_1100_spdd[201];

    auto g_x_y_0_0_0_z_yy_yz = buffer_1100_spdd[202];

    auto g_x_y_0_0_0_z_yy_zz = buffer_1100_spdd[203];

    auto g_x_y_0_0_0_z_yz_xx = buffer_1100_spdd[204];

    auto g_x_y_0_0_0_z_yz_xy = buffer_1100_spdd[205];

    auto g_x_y_0_0_0_z_yz_xz = buffer_1100_spdd[206];

    auto g_x_y_0_0_0_z_yz_yy = buffer_1100_spdd[207];

    auto g_x_y_0_0_0_z_yz_yz = buffer_1100_spdd[208];

    auto g_x_y_0_0_0_z_yz_zz = buffer_1100_spdd[209];

    auto g_x_y_0_0_0_z_zz_xx = buffer_1100_spdd[210];

    auto g_x_y_0_0_0_z_zz_xy = buffer_1100_spdd[211];

    auto g_x_y_0_0_0_z_zz_xz = buffer_1100_spdd[212];

    auto g_x_y_0_0_0_z_zz_yy = buffer_1100_spdd[213];

    auto g_x_y_0_0_0_z_zz_yz = buffer_1100_spdd[214];

    auto g_x_y_0_0_0_z_zz_zz = buffer_1100_spdd[215];

    auto g_x_z_0_0_0_x_xx_xx = buffer_1100_spdd[216];

    auto g_x_z_0_0_0_x_xx_xy = buffer_1100_spdd[217];

    auto g_x_z_0_0_0_x_xx_xz = buffer_1100_spdd[218];

    auto g_x_z_0_0_0_x_xx_yy = buffer_1100_spdd[219];

    auto g_x_z_0_0_0_x_xx_yz = buffer_1100_spdd[220];

    auto g_x_z_0_0_0_x_xx_zz = buffer_1100_spdd[221];

    auto g_x_z_0_0_0_x_xy_xx = buffer_1100_spdd[222];

    auto g_x_z_0_0_0_x_xy_xy = buffer_1100_spdd[223];

    auto g_x_z_0_0_0_x_xy_xz = buffer_1100_spdd[224];

    auto g_x_z_0_0_0_x_xy_yy = buffer_1100_spdd[225];

    auto g_x_z_0_0_0_x_xy_yz = buffer_1100_spdd[226];

    auto g_x_z_0_0_0_x_xy_zz = buffer_1100_spdd[227];

    auto g_x_z_0_0_0_x_xz_xx = buffer_1100_spdd[228];

    auto g_x_z_0_0_0_x_xz_xy = buffer_1100_spdd[229];

    auto g_x_z_0_0_0_x_xz_xz = buffer_1100_spdd[230];

    auto g_x_z_0_0_0_x_xz_yy = buffer_1100_spdd[231];

    auto g_x_z_0_0_0_x_xz_yz = buffer_1100_spdd[232];

    auto g_x_z_0_0_0_x_xz_zz = buffer_1100_spdd[233];

    auto g_x_z_0_0_0_x_yy_xx = buffer_1100_spdd[234];

    auto g_x_z_0_0_0_x_yy_xy = buffer_1100_spdd[235];

    auto g_x_z_0_0_0_x_yy_xz = buffer_1100_spdd[236];

    auto g_x_z_0_0_0_x_yy_yy = buffer_1100_spdd[237];

    auto g_x_z_0_0_0_x_yy_yz = buffer_1100_spdd[238];

    auto g_x_z_0_0_0_x_yy_zz = buffer_1100_spdd[239];

    auto g_x_z_0_0_0_x_yz_xx = buffer_1100_spdd[240];

    auto g_x_z_0_0_0_x_yz_xy = buffer_1100_spdd[241];

    auto g_x_z_0_0_0_x_yz_xz = buffer_1100_spdd[242];

    auto g_x_z_0_0_0_x_yz_yy = buffer_1100_spdd[243];

    auto g_x_z_0_0_0_x_yz_yz = buffer_1100_spdd[244];

    auto g_x_z_0_0_0_x_yz_zz = buffer_1100_spdd[245];

    auto g_x_z_0_0_0_x_zz_xx = buffer_1100_spdd[246];

    auto g_x_z_0_0_0_x_zz_xy = buffer_1100_spdd[247];

    auto g_x_z_0_0_0_x_zz_xz = buffer_1100_spdd[248];

    auto g_x_z_0_0_0_x_zz_yy = buffer_1100_spdd[249];

    auto g_x_z_0_0_0_x_zz_yz = buffer_1100_spdd[250];

    auto g_x_z_0_0_0_x_zz_zz = buffer_1100_spdd[251];

    auto g_x_z_0_0_0_y_xx_xx = buffer_1100_spdd[252];

    auto g_x_z_0_0_0_y_xx_xy = buffer_1100_spdd[253];

    auto g_x_z_0_0_0_y_xx_xz = buffer_1100_spdd[254];

    auto g_x_z_0_0_0_y_xx_yy = buffer_1100_spdd[255];

    auto g_x_z_0_0_0_y_xx_yz = buffer_1100_spdd[256];

    auto g_x_z_0_0_0_y_xx_zz = buffer_1100_spdd[257];

    auto g_x_z_0_0_0_y_xy_xx = buffer_1100_spdd[258];

    auto g_x_z_0_0_0_y_xy_xy = buffer_1100_spdd[259];

    auto g_x_z_0_0_0_y_xy_xz = buffer_1100_spdd[260];

    auto g_x_z_0_0_0_y_xy_yy = buffer_1100_spdd[261];

    auto g_x_z_0_0_0_y_xy_yz = buffer_1100_spdd[262];

    auto g_x_z_0_0_0_y_xy_zz = buffer_1100_spdd[263];

    auto g_x_z_0_0_0_y_xz_xx = buffer_1100_spdd[264];

    auto g_x_z_0_0_0_y_xz_xy = buffer_1100_spdd[265];

    auto g_x_z_0_0_0_y_xz_xz = buffer_1100_spdd[266];

    auto g_x_z_0_0_0_y_xz_yy = buffer_1100_spdd[267];

    auto g_x_z_0_0_0_y_xz_yz = buffer_1100_spdd[268];

    auto g_x_z_0_0_0_y_xz_zz = buffer_1100_spdd[269];

    auto g_x_z_0_0_0_y_yy_xx = buffer_1100_spdd[270];

    auto g_x_z_0_0_0_y_yy_xy = buffer_1100_spdd[271];

    auto g_x_z_0_0_0_y_yy_xz = buffer_1100_spdd[272];

    auto g_x_z_0_0_0_y_yy_yy = buffer_1100_spdd[273];

    auto g_x_z_0_0_0_y_yy_yz = buffer_1100_spdd[274];

    auto g_x_z_0_0_0_y_yy_zz = buffer_1100_spdd[275];

    auto g_x_z_0_0_0_y_yz_xx = buffer_1100_spdd[276];

    auto g_x_z_0_0_0_y_yz_xy = buffer_1100_spdd[277];

    auto g_x_z_0_0_0_y_yz_xz = buffer_1100_spdd[278];

    auto g_x_z_0_0_0_y_yz_yy = buffer_1100_spdd[279];

    auto g_x_z_0_0_0_y_yz_yz = buffer_1100_spdd[280];

    auto g_x_z_0_0_0_y_yz_zz = buffer_1100_spdd[281];

    auto g_x_z_0_0_0_y_zz_xx = buffer_1100_spdd[282];

    auto g_x_z_0_0_0_y_zz_xy = buffer_1100_spdd[283];

    auto g_x_z_0_0_0_y_zz_xz = buffer_1100_spdd[284];

    auto g_x_z_0_0_0_y_zz_yy = buffer_1100_spdd[285];

    auto g_x_z_0_0_0_y_zz_yz = buffer_1100_spdd[286];

    auto g_x_z_0_0_0_y_zz_zz = buffer_1100_spdd[287];

    auto g_x_z_0_0_0_z_xx_xx = buffer_1100_spdd[288];

    auto g_x_z_0_0_0_z_xx_xy = buffer_1100_spdd[289];

    auto g_x_z_0_0_0_z_xx_xz = buffer_1100_spdd[290];

    auto g_x_z_0_0_0_z_xx_yy = buffer_1100_spdd[291];

    auto g_x_z_0_0_0_z_xx_yz = buffer_1100_spdd[292];

    auto g_x_z_0_0_0_z_xx_zz = buffer_1100_spdd[293];

    auto g_x_z_0_0_0_z_xy_xx = buffer_1100_spdd[294];

    auto g_x_z_0_0_0_z_xy_xy = buffer_1100_spdd[295];

    auto g_x_z_0_0_0_z_xy_xz = buffer_1100_spdd[296];

    auto g_x_z_0_0_0_z_xy_yy = buffer_1100_spdd[297];

    auto g_x_z_0_0_0_z_xy_yz = buffer_1100_spdd[298];

    auto g_x_z_0_0_0_z_xy_zz = buffer_1100_spdd[299];

    auto g_x_z_0_0_0_z_xz_xx = buffer_1100_spdd[300];

    auto g_x_z_0_0_0_z_xz_xy = buffer_1100_spdd[301];

    auto g_x_z_0_0_0_z_xz_xz = buffer_1100_spdd[302];

    auto g_x_z_0_0_0_z_xz_yy = buffer_1100_spdd[303];

    auto g_x_z_0_0_0_z_xz_yz = buffer_1100_spdd[304];

    auto g_x_z_0_0_0_z_xz_zz = buffer_1100_spdd[305];

    auto g_x_z_0_0_0_z_yy_xx = buffer_1100_spdd[306];

    auto g_x_z_0_0_0_z_yy_xy = buffer_1100_spdd[307];

    auto g_x_z_0_0_0_z_yy_xz = buffer_1100_spdd[308];

    auto g_x_z_0_0_0_z_yy_yy = buffer_1100_spdd[309];

    auto g_x_z_0_0_0_z_yy_yz = buffer_1100_spdd[310];

    auto g_x_z_0_0_0_z_yy_zz = buffer_1100_spdd[311];

    auto g_x_z_0_0_0_z_yz_xx = buffer_1100_spdd[312];

    auto g_x_z_0_0_0_z_yz_xy = buffer_1100_spdd[313];

    auto g_x_z_0_0_0_z_yz_xz = buffer_1100_spdd[314];

    auto g_x_z_0_0_0_z_yz_yy = buffer_1100_spdd[315];

    auto g_x_z_0_0_0_z_yz_yz = buffer_1100_spdd[316];

    auto g_x_z_0_0_0_z_yz_zz = buffer_1100_spdd[317];

    auto g_x_z_0_0_0_z_zz_xx = buffer_1100_spdd[318];

    auto g_x_z_0_0_0_z_zz_xy = buffer_1100_spdd[319];

    auto g_x_z_0_0_0_z_zz_xz = buffer_1100_spdd[320];

    auto g_x_z_0_0_0_z_zz_yy = buffer_1100_spdd[321];

    auto g_x_z_0_0_0_z_zz_yz = buffer_1100_spdd[322];

    auto g_x_z_0_0_0_z_zz_zz = buffer_1100_spdd[323];

    auto g_y_x_0_0_0_x_xx_xx = buffer_1100_spdd[324];

    auto g_y_x_0_0_0_x_xx_xy = buffer_1100_spdd[325];

    auto g_y_x_0_0_0_x_xx_xz = buffer_1100_spdd[326];

    auto g_y_x_0_0_0_x_xx_yy = buffer_1100_spdd[327];

    auto g_y_x_0_0_0_x_xx_yz = buffer_1100_spdd[328];

    auto g_y_x_0_0_0_x_xx_zz = buffer_1100_spdd[329];

    auto g_y_x_0_0_0_x_xy_xx = buffer_1100_spdd[330];

    auto g_y_x_0_0_0_x_xy_xy = buffer_1100_spdd[331];

    auto g_y_x_0_0_0_x_xy_xz = buffer_1100_spdd[332];

    auto g_y_x_0_0_0_x_xy_yy = buffer_1100_spdd[333];

    auto g_y_x_0_0_0_x_xy_yz = buffer_1100_spdd[334];

    auto g_y_x_0_0_0_x_xy_zz = buffer_1100_spdd[335];

    auto g_y_x_0_0_0_x_xz_xx = buffer_1100_spdd[336];

    auto g_y_x_0_0_0_x_xz_xy = buffer_1100_spdd[337];

    auto g_y_x_0_0_0_x_xz_xz = buffer_1100_spdd[338];

    auto g_y_x_0_0_0_x_xz_yy = buffer_1100_spdd[339];

    auto g_y_x_0_0_0_x_xz_yz = buffer_1100_spdd[340];

    auto g_y_x_0_0_0_x_xz_zz = buffer_1100_spdd[341];

    auto g_y_x_0_0_0_x_yy_xx = buffer_1100_spdd[342];

    auto g_y_x_0_0_0_x_yy_xy = buffer_1100_spdd[343];

    auto g_y_x_0_0_0_x_yy_xz = buffer_1100_spdd[344];

    auto g_y_x_0_0_0_x_yy_yy = buffer_1100_spdd[345];

    auto g_y_x_0_0_0_x_yy_yz = buffer_1100_spdd[346];

    auto g_y_x_0_0_0_x_yy_zz = buffer_1100_spdd[347];

    auto g_y_x_0_0_0_x_yz_xx = buffer_1100_spdd[348];

    auto g_y_x_0_0_0_x_yz_xy = buffer_1100_spdd[349];

    auto g_y_x_0_0_0_x_yz_xz = buffer_1100_spdd[350];

    auto g_y_x_0_0_0_x_yz_yy = buffer_1100_spdd[351];

    auto g_y_x_0_0_0_x_yz_yz = buffer_1100_spdd[352];

    auto g_y_x_0_0_0_x_yz_zz = buffer_1100_spdd[353];

    auto g_y_x_0_0_0_x_zz_xx = buffer_1100_spdd[354];

    auto g_y_x_0_0_0_x_zz_xy = buffer_1100_spdd[355];

    auto g_y_x_0_0_0_x_zz_xz = buffer_1100_spdd[356];

    auto g_y_x_0_0_0_x_zz_yy = buffer_1100_spdd[357];

    auto g_y_x_0_0_0_x_zz_yz = buffer_1100_spdd[358];

    auto g_y_x_0_0_0_x_zz_zz = buffer_1100_spdd[359];

    auto g_y_x_0_0_0_y_xx_xx = buffer_1100_spdd[360];

    auto g_y_x_0_0_0_y_xx_xy = buffer_1100_spdd[361];

    auto g_y_x_0_0_0_y_xx_xz = buffer_1100_spdd[362];

    auto g_y_x_0_0_0_y_xx_yy = buffer_1100_spdd[363];

    auto g_y_x_0_0_0_y_xx_yz = buffer_1100_spdd[364];

    auto g_y_x_0_0_0_y_xx_zz = buffer_1100_spdd[365];

    auto g_y_x_0_0_0_y_xy_xx = buffer_1100_spdd[366];

    auto g_y_x_0_0_0_y_xy_xy = buffer_1100_spdd[367];

    auto g_y_x_0_0_0_y_xy_xz = buffer_1100_spdd[368];

    auto g_y_x_0_0_0_y_xy_yy = buffer_1100_spdd[369];

    auto g_y_x_0_0_0_y_xy_yz = buffer_1100_spdd[370];

    auto g_y_x_0_0_0_y_xy_zz = buffer_1100_spdd[371];

    auto g_y_x_0_0_0_y_xz_xx = buffer_1100_spdd[372];

    auto g_y_x_0_0_0_y_xz_xy = buffer_1100_spdd[373];

    auto g_y_x_0_0_0_y_xz_xz = buffer_1100_spdd[374];

    auto g_y_x_0_0_0_y_xz_yy = buffer_1100_spdd[375];

    auto g_y_x_0_0_0_y_xz_yz = buffer_1100_spdd[376];

    auto g_y_x_0_0_0_y_xz_zz = buffer_1100_spdd[377];

    auto g_y_x_0_0_0_y_yy_xx = buffer_1100_spdd[378];

    auto g_y_x_0_0_0_y_yy_xy = buffer_1100_spdd[379];

    auto g_y_x_0_0_0_y_yy_xz = buffer_1100_spdd[380];

    auto g_y_x_0_0_0_y_yy_yy = buffer_1100_spdd[381];

    auto g_y_x_0_0_0_y_yy_yz = buffer_1100_spdd[382];

    auto g_y_x_0_0_0_y_yy_zz = buffer_1100_spdd[383];

    auto g_y_x_0_0_0_y_yz_xx = buffer_1100_spdd[384];

    auto g_y_x_0_0_0_y_yz_xy = buffer_1100_spdd[385];

    auto g_y_x_0_0_0_y_yz_xz = buffer_1100_spdd[386];

    auto g_y_x_0_0_0_y_yz_yy = buffer_1100_spdd[387];

    auto g_y_x_0_0_0_y_yz_yz = buffer_1100_spdd[388];

    auto g_y_x_0_0_0_y_yz_zz = buffer_1100_spdd[389];

    auto g_y_x_0_0_0_y_zz_xx = buffer_1100_spdd[390];

    auto g_y_x_0_0_0_y_zz_xy = buffer_1100_spdd[391];

    auto g_y_x_0_0_0_y_zz_xz = buffer_1100_spdd[392];

    auto g_y_x_0_0_0_y_zz_yy = buffer_1100_spdd[393];

    auto g_y_x_0_0_0_y_zz_yz = buffer_1100_spdd[394];

    auto g_y_x_0_0_0_y_zz_zz = buffer_1100_spdd[395];

    auto g_y_x_0_0_0_z_xx_xx = buffer_1100_spdd[396];

    auto g_y_x_0_0_0_z_xx_xy = buffer_1100_spdd[397];

    auto g_y_x_0_0_0_z_xx_xz = buffer_1100_spdd[398];

    auto g_y_x_0_0_0_z_xx_yy = buffer_1100_spdd[399];

    auto g_y_x_0_0_0_z_xx_yz = buffer_1100_spdd[400];

    auto g_y_x_0_0_0_z_xx_zz = buffer_1100_spdd[401];

    auto g_y_x_0_0_0_z_xy_xx = buffer_1100_spdd[402];

    auto g_y_x_0_0_0_z_xy_xy = buffer_1100_spdd[403];

    auto g_y_x_0_0_0_z_xy_xz = buffer_1100_spdd[404];

    auto g_y_x_0_0_0_z_xy_yy = buffer_1100_spdd[405];

    auto g_y_x_0_0_0_z_xy_yz = buffer_1100_spdd[406];

    auto g_y_x_0_0_0_z_xy_zz = buffer_1100_spdd[407];

    auto g_y_x_0_0_0_z_xz_xx = buffer_1100_spdd[408];

    auto g_y_x_0_0_0_z_xz_xy = buffer_1100_spdd[409];

    auto g_y_x_0_0_0_z_xz_xz = buffer_1100_spdd[410];

    auto g_y_x_0_0_0_z_xz_yy = buffer_1100_spdd[411];

    auto g_y_x_0_0_0_z_xz_yz = buffer_1100_spdd[412];

    auto g_y_x_0_0_0_z_xz_zz = buffer_1100_spdd[413];

    auto g_y_x_0_0_0_z_yy_xx = buffer_1100_spdd[414];

    auto g_y_x_0_0_0_z_yy_xy = buffer_1100_spdd[415];

    auto g_y_x_0_0_0_z_yy_xz = buffer_1100_spdd[416];

    auto g_y_x_0_0_0_z_yy_yy = buffer_1100_spdd[417];

    auto g_y_x_0_0_0_z_yy_yz = buffer_1100_spdd[418];

    auto g_y_x_0_0_0_z_yy_zz = buffer_1100_spdd[419];

    auto g_y_x_0_0_0_z_yz_xx = buffer_1100_spdd[420];

    auto g_y_x_0_0_0_z_yz_xy = buffer_1100_spdd[421];

    auto g_y_x_0_0_0_z_yz_xz = buffer_1100_spdd[422];

    auto g_y_x_0_0_0_z_yz_yy = buffer_1100_spdd[423];

    auto g_y_x_0_0_0_z_yz_yz = buffer_1100_spdd[424];

    auto g_y_x_0_0_0_z_yz_zz = buffer_1100_spdd[425];

    auto g_y_x_0_0_0_z_zz_xx = buffer_1100_spdd[426];

    auto g_y_x_0_0_0_z_zz_xy = buffer_1100_spdd[427];

    auto g_y_x_0_0_0_z_zz_xz = buffer_1100_spdd[428];

    auto g_y_x_0_0_0_z_zz_yy = buffer_1100_spdd[429];

    auto g_y_x_0_0_0_z_zz_yz = buffer_1100_spdd[430];

    auto g_y_x_0_0_0_z_zz_zz = buffer_1100_spdd[431];

    auto g_y_y_0_0_0_x_xx_xx = buffer_1100_spdd[432];

    auto g_y_y_0_0_0_x_xx_xy = buffer_1100_spdd[433];

    auto g_y_y_0_0_0_x_xx_xz = buffer_1100_spdd[434];

    auto g_y_y_0_0_0_x_xx_yy = buffer_1100_spdd[435];

    auto g_y_y_0_0_0_x_xx_yz = buffer_1100_spdd[436];

    auto g_y_y_0_0_0_x_xx_zz = buffer_1100_spdd[437];

    auto g_y_y_0_0_0_x_xy_xx = buffer_1100_spdd[438];

    auto g_y_y_0_0_0_x_xy_xy = buffer_1100_spdd[439];

    auto g_y_y_0_0_0_x_xy_xz = buffer_1100_spdd[440];

    auto g_y_y_0_0_0_x_xy_yy = buffer_1100_spdd[441];

    auto g_y_y_0_0_0_x_xy_yz = buffer_1100_spdd[442];

    auto g_y_y_0_0_0_x_xy_zz = buffer_1100_spdd[443];

    auto g_y_y_0_0_0_x_xz_xx = buffer_1100_spdd[444];

    auto g_y_y_0_0_0_x_xz_xy = buffer_1100_spdd[445];

    auto g_y_y_0_0_0_x_xz_xz = buffer_1100_spdd[446];

    auto g_y_y_0_0_0_x_xz_yy = buffer_1100_spdd[447];

    auto g_y_y_0_0_0_x_xz_yz = buffer_1100_spdd[448];

    auto g_y_y_0_0_0_x_xz_zz = buffer_1100_spdd[449];

    auto g_y_y_0_0_0_x_yy_xx = buffer_1100_spdd[450];

    auto g_y_y_0_0_0_x_yy_xy = buffer_1100_spdd[451];

    auto g_y_y_0_0_0_x_yy_xz = buffer_1100_spdd[452];

    auto g_y_y_0_0_0_x_yy_yy = buffer_1100_spdd[453];

    auto g_y_y_0_0_0_x_yy_yz = buffer_1100_spdd[454];

    auto g_y_y_0_0_0_x_yy_zz = buffer_1100_spdd[455];

    auto g_y_y_0_0_0_x_yz_xx = buffer_1100_spdd[456];

    auto g_y_y_0_0_0_x_yz_xy = buffer_1100_spdd[457];

    auto g_y_y_0_0_0_x_yz_xz = buffer_1100_spdd[458];

    auto g_y_y_0_0_0_x_yz_yy = buffer_1100_spdd[459];

    auto g_y_y_0_0_0_x_yz_yz = buffer_1100_spdd[460];

    auto g_y_y_0_0_0_x_yz_zz = buffer_1100_spdd[461];

    auto g_y_y_0_0_0_x_zz_xx = buffer_1100_spdd[462];

    auto g_y_y_0_0_0_x_zz_xy = buffer_1100_spdd[463];

    auto g_y_y_0_0_0_x_zz_xz = buffer_1100_spdd[464];

    auto g_y_y_0_0_0_x_zz_yy = buffer_1100_spdd[465];

    auto g_y_y_0_0_0_x_zz_yz = buffer_1100_spdd[466];

    auto g_y_y_0_0_0_x_zz_zz = buffer_1100_spdd[467];

    auto g_y_y_0_0_0_y_xx_xx = buffer_1100_spdd[468];

    auto g_y_y_0_0_0_y_xx_xy = buffer_1100_spdd[469];

    auto g_y_y_0_0_0_y_xx_xz = buffer_1100_spdd[470];

    auto g_y_y_0_0_0_y_xx_yy = buffer_1100_spdd[471];

    auto g_y_y_0_0_0_y_xx_yz = buffer_1100_spdd[472];

    auto g_y_y_0_0_0_y_xx_zz = buffer_1100_spdd[473];

    auto g_y_y_0_0_0_y_xy_xx = buffer_1100_spdd[474];

    auto g_y_y_0_0_0_y_xy_xy = buffer_1100_spdd[475];

    auto g_y_y_0_0_0_y_xy_xz = buffer_1100_spdd[476];

    auto g_y_y_0_0_0_y_xy_yy = buffer_1100_spdd[477];

    auto g_y_y_0_0_0_y_xy_yz = buffer_1100_spdd[478];

    auto g_y_y_0_0_0_y_xy_zz = buffer_1100_spdd[479];

    auto g_y_y_0_0_0_y_xz_xx = buffer_1100_spdd[480];

    auto g_y_y_0_0_0_y_xz_xy = buffer_1100_spdd[481];

    auto g_y_y_0_0_0_y_xz_xz = buffer_1100_spdd[482];

    auto g_y_y_0_0_0_y_xz_yy = buffer_1100_spdd[483];

    auto g_y_y_0_0_0_y_xz_yz = buffer_1100_spdd[484];

    auto g_y_y_0_0_0_y_xz_zz = buffer_1100_spdd[485];

    auto g_y_y_0_0_0_y_yy_xx = buffer_1100_spdd[486];

    auto g_y_y_0_0_0_y_yy_xy = buffer_1100_spdd[487];

    auto g_y_y_0_0_0_y_yy_xz = buffer_1100_spdd[488];

    auto g_y_y_0_0_0_y_yy_yy = buffer_1100_spdd[489];

    auto g_y_y_0_0_0_y_yy_yz = buffer_1100_spdd[490];

    auto g_y_y_0_0_0_y_yy_zz = buffer_1100_spdd[491];

    auto g_y_y_0_0_0_y_yz_xx = buffer_1100_spdd[492];

    auto g_y_y_0_0_0_y_yz_xy = buffer_1100_spdd[493];

    auto g_y_y_0_0_0_y_yz_xz = buffer_1100_spdd[494];

    auto g_y_y_0_0_0_y_yz_yy = buffer_1100_spdd[495];

    auto g_y_y_0_0_0_y_yz_yz = buffer_1100_spdd[496];

    auto g_y_y_0_0_0_y_yz_zz = buffer_1100_spdd[497];

    auto g_y_y_0_0_0_y_zz_xx = buffer_1100_spdd[498];

    auto g_y_y_0_0_0_y_zz_xy = buffer_1100_spdd[499];

    auto g_y_y_0_0_0_y_zz_xz = buffer_1100_spdd[500];

    auto g_y_y_0_0_0_y_zz_yy = buffer_1100_spdd[501];

    auto g_y_y_0_0_0_y_zz_yz = buffer_1100_spdd[502];

    auto g_y_y_0_0_0_y_zz_zz = buffer_1100_spdd[503];

    auto g_y_y_0_0_0_z_xx_xx = buffer_1100_spdd[504];

    auto g_y_y_0_0_0_z_xx_xy = buffer_1100_spdd[505];

    auto g_y_y_0_0_0_z_xx_xz = buffer_1100_spdd[506];

    auto g_y_y_0_0_0_z_xx_yy = buffer_1100_spdd[507];

    auto g_y_y_0_0_0_z_xx_yz = buffer_1100_spdd[508];

    auto g_y_y_0_0_0_z_xx_zz = buffer_1100_spdd[509];

    auto g_y_y_0_0_0_z_xy_xx = buffer_1100_spdd[510];

    auto g_y_y_0_0_0_z_xy_xy = buffer_1100_spdd[511];

    auto g_y_y_0_0_0_z_xy_xz = buffer_1100_spdd[512];

    auto g_y_y_0_0_0_z_xy_yy = buffer_1100_spdd[513];

    auto g_y_y_0_0_0_z_xy_yz = buffer_1100_spdd[514];

    auto g_y_y_0_0_0_z_xy_zz = buffer_1100_spdd[515];

    auto g_y_y_0_0_0_z_xz_xx = buffer_1100_spdd[516];

    auto g_y_y_0_0_0_z_xz_xy = buffer_1100_spdd[517];

    auto g_y_y_0_0_0_z_xz_xz = buffer_1100_spdd[518];

    auto g_y_y_0_0_0_z_xz_yy = buffer_1100_spdd[519];

    auto g_y_y_0_0_0_z_xz_yz = buffer_1100_spdd[520];

    auto g_y_y_0_0_0_z_xz_zz = buffer_1100_spdd[521];

    auto g_y_y_0_0_0_z_yy_xx = buffer_1100_spdd[522];

    auto g_y_y_0_0_0_z_yy_xy = buffer_1100_spdd[523];

    auto g_y_y_0_0_0_z_yy_xz = buffer_1100_spdd[524];

    auto g_y_y_0_0_0_z_yy_yy = buffer_1100_spdd[525];

    auto g_y_y_0_0_0_z_yy_yz = buffer_1100_spdd[526];

    auto g_y_y_0_0_0_z_yy_zz = buffer_1100_spdd[527];

    auto g_y_y_0_0_0_z_yz_xx = buffer_1100_spdd[528];

    auto g_y_y_0_0_0_z_yz_xy = buffer_1100_spdd[529];

    auto g_y_y_0_0_0_z_yz_xz = buffer_1100_spdd[530];

    auto g_y_y_0_0_0_z_yz_yy = buffer_1100_spdd[531];

    auto g_y_y_0_0_0_z_yz_yz = buffer_1100_spdd[532];

    auto g_y_y_0_0_0_z_yz_zz = buffer_1100_spdd[533];

    auto g_y_y_0_0_0_z_zz_xx = buffer_1100_spdd[534];

    auto g_y_y_0_0_0_z_zz_xy = buffer_1100_spdd[535];

    auto g_y_y_0_0_0_z_zz_xz = buffer_1100_spdd[536];

    auto g_y_y_0_0_0_z_zz_yy = buffer_1100_spdd[537];

    auto g_y_y_0_0_0_z_zz_yz = buffer_1100_spdd[538];

    auto g_y_y_0_0_0_z_zz_zz = buffer_1100_spdd[539];

    auto g_y_z_0_0_0_x_xx_xx = buffer_1100_spdd[540];

    auto g_y_z_0_0_0_x_xx_xy = buffer_1100_spdd[541];

    auto g_y_z_0_0_0_x_xx_xz = buffer_1100_spdd[542];

    auto g_y_z_0_0_0_x_xx_yy = buffer_1100_spdd[543];

    auto g_y_z_0_0_0_x_xx_yz = buffer_1100_spdd[544];

    auto g_y_z_0_0_0_x_xx_zz = buffer_1100_spdd[545];

    auto g_y_z_0_0_0_x_xy_xx = buffer_1100_spdd[546];

    auto g_y_z_0_0_0_x_xy_xy = buffer_1100_spdd[547];

    auto g_y_z_0_0_0_x_xy_xz = buffer_1100_spdd[548];

    auto g_y_z_0_0_0_x_xy_yy = buffer_1100_spdd[549];

    auto g_y_z_0_0_0_x_xy_yz = buffer_1100_spdd[550];

    auto g_y_z_0_0_0_x_xy_zz = buffer_1100_spdd[551];

    auto g_y_z_0_0_0_x_xz_xx = buffer_1100_spdd[552];

    auto g_y_z_0_0_0_x_xz_xy = buffer_1100_spdd[553];

    auto g_y_z_0_0_0_x_xz_xz = buffer_1100_spdd[554];

    auto g_y_z_0_0_0_x_xz_yy = buffer_1100_spdd[555];

    auto g_y_z_0_0_0_x_xz_yz = buffer_1100_spdd[556];

    auto g_y_z_0_0_0_x_xz_zz = buffer_1100_spdd[557];

    auto g_y_z_0_0_0_x_yy_xx = buffer_1100_spdd[558];

    auto g_y_z_0_0_0_x_yy_xy = buffer_1100_spdd[559];

    auto g_y_z_0_0_0_x_yy_xz = buffer_1100_spdd[560];

    auto g_y_z_0_0_0_x_yy_yy = buffer_1100_spdd[561];

    auto g_y_z_0_0_0_x_yy_yz = buffer_1100_spdd[562];

    auto g_y_z_0_0_0_x_yy_zz = buffer_1100_spdd[563];

    auto g_y_z_0_0_0_x_yz_xx = buffer_1100_spdd[564];

    auto g_y_z_0_0_0_x_yz_xy = buffer_1100_spdd[565];

    auto g_y_z_0_0_0_x_yz_xz = buffer_1100_spdd[566];

    auto g_y_z_0_0_0_x_yz_yy = buffer_1100_spdd[567];

    auto g_y_z_0_0_0_x_yz_yz = buffer_1100_spdd[568];

    auto g_y_z_0_0_0_x_yz_zz = buffer_1100_spdd[569];

    auto g_y_z_0_0_0_x_zz_xx = buffer_1100_spdd[570];

    auto g_y_z_0_0_0_x_zz_xy = buffer_1100_spdd[571];

    auto g_y_z_0_0_0_x_zz_xz = buffer_1100_spdd[572];

    auto g_y_z_0_0_0_x_zz_yy = buffer_1100_spdd[573];

    auto g_y_z_0_0_0_x_zz_yz = buffer_1100_spdd[574];

    auto g_y_z_0_0_0_x_zz_zz = buffer_1100_spdd[575];

    auto g_y_z_0_0_0_y_xx_xx = buffer_1100_spdd[576];

    auto g_y_z_0_0_0_y_xx_xy = buffer_1100_spdd[577];

    auto g_y_z_0_0_0_y_xx_xz = buffer_1100_spdd[578];

    auto g_y_z_0_0_0_y_xx_yy = buffer_1100_spdd[579];

    auto g_y_z_0_0_0_y_xx_yz = buffer_1100_spdd[580];

    auto g_y_z_0_0_0_y_xx_zz = buffer_1100_spdd[581];

    auto g_y_z_0_0_0_y_xy_xx = buffer_1100_spdd[582];

    auto g_y_z_0_0_0_y_xy_xy = buffer_1100_spdd[583];

    auto g_y_z_0_0_0_y_xy_xz = buffer_1100_spdd[584];

    auto g_y_z_0_0_0_y_xy_yy = buffer_1100_spdd[585];

    auto g_y_z_0_0_0_y_xy_yz = buffer_1100_spdd[586];

    auto g_y_z_0_0_0_y_xy_zz = buffer_1100_spdd[587];

    auto g_y_z_0_0_0_y_xz_xx = buffer_1100_spdd[588];

    auto g_y_z_0_0_0_y_xz_xy = buffer_1100_spdd[589];

    auto g_y_z_0_0_0_y_xz_xz = buffer_1100_spdd[590];

    auto g_y_z_0_0_0_y_xz_yy = buffer_1100_spdd[591];

    auto g_y_z_0_0_0_y_xz_yz = buffer_1100_spdd[592];

    auto g_y_z_0_0_0_y_xz_zz = buffer_1100_spdd[593];

    auto g_y_z_0_0_0_y_yy_xx = buffer_1100_spdd[594];

    auto g_y_z_0_0_0_y_yy_xy = buffer_1100_spdd[595];

    auto g_y_z_0_0_0_y_yy_xz = buffer_1100_spdd[596];

    auto g_y_z_0_0_0_y_yy_yy = buffer_1100_spdd[597];

    auto g_y_z_0_0_0_y_yy_yz = buffer_1100_spdd[598];

    auto g_y_z_0_0_0_y_yy_zz = buffer_1100_spdd[599];

    auto g_y_z_0_0_0_y_yz_xx = buffer_1100_spdd[600];

    auto g_y_z_0_0_0_y_yz_xy = buffer_1100_spdd[601];

    auto g_y_z_0_0_0_y_yz_xz = buffer_1100_spdd[602];

    auto g_y_z_0_0_0_y_yz_yy = buffer_1100_spdd[603];

    auto g_y_z_0_0_0_y_yz_yz = buffer_1100_spdd[604];

    auto g_y_z_0_0_0_y_yz_zz = buffer_1100_spdd[605];

    auto g_y_z_0_0_0_y_zz_xx = buffer_1100_spdd[606];

    auto g_y_z_0_0_0_y_zz_xy = buffer_1100_spdd[607];

    auto g_y_z_0_0_0_y_zz_xz = buffer_1100_spdd[608];

    auto g_y_z_0_0_0_y_zz_yy = buffer_1100_spdd[609];

    auto g_y_z_0_0_0_y_zz_yz = buffer_1100_spdd[610];

    auto g_y_z_0_0_0_y_zz_zz = buffer_1100_spdd[611];

    auto g_y_z_0_0_0_z_xx_xx = buffer_1100_spdd[612];

    auto g_y_z_0_0_0_z_xx_xy = buffer_1100_spdd[613];

    auto g_y_z_0_0_0_z_xx_xz = buffer_1100_spdd[614];

    auto g_y_z_0_0_0_z_xx_yy = buffer_1100_spdd[615];

    auto g_y_z_0_0_0_z_xx_yz = buffer_1100_spdd[616];

    auto g_y_z_0_0_0_z_xx_zz = buffer_1100_spdd[617];

    auto g_y_z_0_0_0_z_xy_xx = buffer_1100_spdd[618];

    auto g_y_z_0_0_0_z_xy_xy = buffer_1100_spdd[619];

    auto g_y_z_0_0_0_z_xy_xz = buffer_1100_spdd[620];

    auto g_y_z_0_0_0_z_xy_yy = buffer_1100_spdd[621];

    auto g_y_z_0_0_0_z_xy_yz = buffer_1100_spdd[622];

    auto g_y_z_0_0_0_z_xy_zz = buffer_1100_spdd[623];

    auto g_y_z_0_0_0_z_xz_xx = buffer_1100_spdd[624];

    auto g_y_z_0_0_0_z_xz_xy = buffer_1100_spdd[625];

    auto g_y_z_0_0_0_z_xz_xz = buffer_1100_spdd[626];

    auto g_y_z_0_0_0_z_xz_yy = buffer_1100_spdd[627];

    auto g_y_z_0_0_0_z_xz_yz = buffer_1100_spdd[628];

    auto g_y_z_0_0_0_z_xz_zz = buffer_1100_spdd[629];

    auto g_y_z_0_0_0_z_yy_xx = buffer_1100_spdd[630];

    auto g_y_z_0_0_0_z_yy_xy = buffer_1100_spdd[631];

    auto g_y_z_0_0_0_z_yy_xz = buffer_1100_spdd[632];

    auto g_y_z_0_0_0_z_yy_yy = buffer_1100_spdd[633];

    auto g_y_z_0_0_0_z_yy_yz = buffer_1100_spdd[634];

    auto g_y_z_0_0_0_z_yy_zz = buffer_1100_spdd[635];

    auto g_y_z_0_0_0_z_yz_xx = buffer_1100_spdd[636];

    auto g_y_z_0_0_0_z_yz_xy = buffer_1100_spdd[637];

    auto g_y_z_0_0_0_z_yz_xz = buffer_1100_spdd[638];

    auto g_y_z_0_0_0_z_yz_yy = buffer_1100_spdd[639];

    auto g_y_z_0_0_0_z_yz_yz = buffer_1100_spdd[640];

    auto g_y_z_0_0_0_z_yz_zz = buffer_1100_spdd[641];

    auto g_y_z_0_0_0_z_zz_xx = buffer_1100_spdd[642];

    auto g_y_z_0_0_0_z_zz_xy = buffer_1100_spdd[643];

    auto g_y_z_0_0_0_z_zz_xz = buffer_1100_spdd[644];

    auto g_y_z_0_0_0_z_zz_yy = buffer_1100_spdd[645];

    auto g_y_z_0_0_0_z_zz_yz = buffer_1100_spdd[646];

    auto g_y_z_0_0_0_z_zz_zz = buffer_1100_spdd[647];

    auto g_z_x_0_0_0_x_xx_xx = buffer_1100_spdd[648];

    auto g_z_x_0_0_0_x_xx_xy = buffer_1100_spdd[649];

    auto g_z_x_0_0_0_x_xx_xz = buffer_1100_spdd[650];

    auto g_z_x_0_0_0_x_xx_yy = buffer_1100_spdd[651];

    auto g_z_x_0_0_0_x_xx_yz = buffer_1100_spdd[652];

    auto g_z_x_0_0_0_x_xx_zz = buffer_1100_spdd[653];

    auto g_z_x_0_0_0_x_xy_xx = buffer_1100_spdd[654];

    auto g_z_x_0_0_0_x_xy_xy = buffer_1100_spdd[655];

    auto g_z_x_0_0_0_x_xy_xz = buffer_1100_spdd[656];

    auto g_z_x_0_0_0_x_xy_yy = buffer_1100_spdd[657];

    auto g_z_x_0_0_0_x_xy_yz = buffer_1100_spdd[658];

    auto g_z_x_0_0_0_x_xy_zz = buffer_1100_spdd[659];

    auto g_z_x_0_0_0_x_xz_xx = buffer_1100_spdd[660];

    auto g_z_x_0_0_0_x_xz_xy = buffer_1100_spdd[661];

    auto g_z_x_0_0_0_x_xz_xz = buffer_1100_spdd[662];

    auto g_z_x_0_0_0_x_xz_yy = buffer_1100_spdd[663];

    auto g_z_x_0_0_0_x_xz_yz = buffer_1100_spdd[664];

    auto g_z_x_0_0_0_x_xz_zz = buffer_1100_spdd[665];

    auto g_z_x_0_0_0_x_yy_xx = buffer_1100_spdd[666];

    auto g_z_x_0_0_0_x_yy_xy = buffer_1100_spdd[667];

    auto g_z_x_0_0_0_x_yy_xz = buffer_1100_spdd[668];

    auto g_z_x_0_0_0_x_yy_yy = buffer_1100_spdd[669];

    auto g_z_x_0_0_0_x_yy_yz = buffer_1100_spdd[670];

    auto g_z_x_0_0_0_x_yy_zz = buffer_1100_spdd[671];

    auto g_z_x_0_0_0_x_yz_xx = buffer_1100_spdd[672];

    auto g_z_x_0_0_0_x_yz_xy = buffer_1100_spdd[673];

    auto g_z_x_0_0_0_x_yz_xz = buffer_1100_spdd[674];

    auto g_z_x_0_0_0_x_yz_yy = buffer_1100_spdd[675];

    auto g_z_x_0_0_0_x_yz_yz = buffer_1100_spdd[676];

    auto g_z_x_0_0_0_x_yz_zz = buffer_1100_spdd[677];

    auto g_z_x_0_0_0_x_zz_xx = buffer_1100_spdd[678];

    auto g_z_x_0_0_0_x_zz_xy = buffer_1100_spdd[679];

    auto g_z_x_0_0_0_x_zz_xz = buffer_1100_spdd[680];

    auto g_z_x_0_0_0_x_zz_yy = buffer_1100_spdd[681];

    auto g_z_x_0_0_0_x_zz_yz = buffer_1100_spdd[682];

    auto g_z_x_0_0_0_x_zz_zz = buffer_1100_spdd[683];

    auto g_z_x_0_0_0_y_xx_xx = buffer_1100_spdd[684];

    auto g_z_x_0_0_0_y_xx_xy = buffer_1100_spdd[685];

    auto g_z_x_0_0_0_y_xx_xz = buffer_1100_spdd[686];

    auto g_z_x_0_0_0_y_xx_yy = buffer_1100_spdd[687];

    auto g_z_x_0_0_0_y_xx_yz = buffer_1100_spdd[688];

    auto g_z_x_0_0_0_y_xx_zz = buffer_1100_spdd[689];

    auto g_z_x_0_0_0_y_xy_xx = buffer_1100_spdd[690];

    auto g_z_x_0_0_0_y_xy_xy = buffer_1100_spdd[691];

    auto g_z_x_0_0_0_y_xy_xz = buffer_1100_spdd[692];

    auto g_z_x_0_0_0_y_xy_yy = buffer_1100_spdd[693];

    auto g_z_x_0_0_0_y_xy_yz = buffer_1100_spdd[694];

    auto g_z_x_0_0_0_y_xy_zz = buffer_1100_spdd[695];

    auto g_z_x_0_0_0_y_xz_xx = buffer_1100_spdd[696];

    auto g_z_x_0_0_0_y_xz_xy = buffer_1100_spdd[697];

    auto g_z_x_0_0_0_y_xz_xz = buffer_1100_spdd[698];

    auto g_z_x_0_0_0_y_xz_yy = buffer_1100_spdd[699];

    auto g_z_x_0_0_0_y_xz_yz = buffer_1100_spdd[700];

    auto g_z_x_0_0_0_y_xz_zz = buffer_1100_spdd[701];

    auto g_z_x_0_0_0_y_yy_xx = buffer_1100_spdd[702];

    auto g_z_x_0_0_0_y_yy_xy = buffer_1100_spdd[703];

    auto g_z_x_0_0_0_y_yy_xz = buffer_1100_spdd[704];

    auto g_z_x_0_0_0_y_yy_yy = buffer_1100_spdd[705];

    auto g_z_x_0_0_0_y_yy_yz = buffer_1100_spdd[706];

    auto g_z_x_0_0_0_y_yy_zz = buffer_1100_spdd[707];

    auto g_z_x_0_0_0_y_yz_xx = buffer_1100_spdd[708];

    auto g_z_x_0_0_0_y_yz_xy = buffer_1100_spdd[709];

    auto g_z_x_0_0_0_y_yz_xz = buffer_1100_spdd[710];

    auto g_z_x_0_0_0_y_yz_yy = buffer_1100_spdd[711];

    auto g_z_x_0_0_0_y_yz_yz = buffer_1100_spdd[712];

    auto g_z_x_0_0_0_y_yz_zz = buffer_1100_spdd[713];

    auto g_z_x_0_0_0_y_zz_xx = buffer_1100_spdd[714];

    auto g_z_x_0_0_0_y_zz_xy = buffer_1100_spdd[715];

    auto g_z_x_0_0_0_y_zz_xz = buffer_1100_spdd[716];

    auto g_z_x_0_0_0_y_zz_yy = buffer_1100_spdd[717];

    auto g_z_x_0_0_0_y_zz_yz = buffer_1100_spdd[718];

    auto g_z_x_0_0_0_y_zz_zz = buffer_1100_spdd[719];

    auto g_z_x_0_0_0_z_xx_xx = buffer_1100_spdd[720];

    auto g_z_x_0_0_0_z_xx_xy = buffer_1100_spdd[721];

    auto g_z_x_0_0_0_z_xx_xz = buffer_1100_spdd[722];

    auto g_z_x_0_0_0_z_xx_yy = buffer_1100_spdd[723];

    auto g_z_x_0_0_0_z_xx_yz = buffer_1100_spdd[724];

    auto g_z_x_0_0_0_z_xx_zz = buffer_1100_spdd[725];

    auto g_z_x_0_0_0_z_xy_xx = buffer_1100_spdd[726];

    auto g_z_x_0_0_0_z_xy_xy = buffer_1100_spdd[727];

    auto g_z_x_0_0_0_z_xy_xz = buffer_1100_spdd[728];

    auto g_z_x_0_0_0_z_xy_yy = buffer_1100_spdd[729];

    auto g_z_x_0_0_0_z_xy_yz = buffer_1100_spdd[730];

    auto g_z_x_0_0_0_z_xy_zz = buffer_1100_spdd[731];

    auto g_z_x_0_0_0_z_xz_xx = buffer_1100_spdd[732];

    auto g_z_x_0_0_0_z_xz_xy = buffer_1100_spdd[733];

    auto g_z_x_0_0_0_z_xz_xz = buffer_1100_spdd[734];

    auto g_z_x_0_0_0_z_xz_yy = buffer_1100_spdd[735];

    auto g_z_x_0_0_0_z_xz_yz = buffer_1100_spdd[736];

    auto g_z_x_0_0_0_z_xz_zz = buffer_1100_spdd[737];

    auto g_z_x_0_0_0_z_yy_xx = buffer_1100_spdd[738];

    auto g_z_x_0_0_0_z_yy_xy = buffer_1100_spdd[739];

    auto g_z_x_0_0_0_z_yy_xz = buffer_1100_spdd[740];

    auto g_z_x_0_0_0_z_yy_yy = buffer_1100_spdd[741];

    auto g_z_x_0_0_0_z_yy_yz = buffer_1100_spdd[742];

    auto g_z_x_0_0_0_z_yy_zz = buffer_1100_spdd[743];

    auto g_z_x_0_0_0_z_yz_xx = buffer_1100_spdd[744];

    auto g_z_x_0_0_0_z_yz_xy = buffer_1100_spdd[745];

    auto g_z_x_0_0_0_z_yz_xz = buffer_1100_spdd[746];

    auto g_z_x_0_0_0_z_yz_yy = buffer_1100_spdd[747];

    auto g_z_x_0_0_0_z_yz_yz = buffer_1100_spdd[748];

    auto g_z_x_0_0_0_z_yz_zz = buffer_1100_spdd[749];

    auto g_z_x_0_0_0_z_zz_xx = buffer_1100_spdd[750];

    auto g_z_x_0_0_0_z_zz_xy = buffer_1100_spdd[751];

    auto g_z_x_0_0_0_z_zz_xz = buffer_1100_spdd[752];

    auto g_z_x_0_0_0_z_zz_yy = buffer_1100_spdd[753];

    auto g_z_x_0_0_0_z_zz_yz = buffer_1100_spdd[754];

    auto g_z_x_0_0_0_z_zz_zz = buffer_1100_spdd[755];

    auto g_z_y_0_0_0_x_xx_xx = buffer_1100_spdd[756];

    auto g_z_y_0_0_0_x_xx_xy = buffer_1100_spdd[757];

    auto g_z_y_0_0_0_x_xx_xz = buffer_1100_spdd[758];

    auto g_z_y_0_0_0_x_xx_yy = buffer_1100_spdd[759];

    auto g_z_y_0_0_0_x_xx_yz = buffer_1100_spdd[760];

    auto g_z_y_0_0_0_x_xx_zz = buffer_1100_spdd[761];

    auto g_z_y_0_0_0_x_xy_xx = buffer_1100_spdd[762];

    auto g_z_y_0_0_0_x_xy_xy = buffer_1100_spdd[763];

    auto g_z_y_0_0_0_x_xy_xz = buffer_1100_spdd[764];

    auto g_z_y_0_0_0_x_xy_yy = buffer_1100_spdd[765];

    auto g_z_y_0_0_0_x_xy_yz = buffer_1100_spdd[766];

    auto g_z_y_0_0_0_x_xy_zz = buffer_1100_spdd[767];

    auto g_z_y_0_0_0_x_xz_xx = buffer_1100_spdd[768];

    auto g_z_y_0_0_0_x_xz_xy = buffer_1100_spdd[769];

    auto g_z_y_0_0_0_x_xz_xz = buffer_1100_spdd[770];

    auto g_z_y_0_0_0_x_xz_yy = buffer_1100_spdd[771];

    auto g_z_y_0_0_0_x_xz_yz = buffer_1100_spdd[772];

    auto g_z_y_0_0_0_x_xz_zz = buffer_1100_spdd[773];

    auto g_z_y_0_0_0_x_yy_xx = buffer_1100_spdd[774];

    auto g_z_y_0_0_0_x_yy_xy = buffer_1100_spdd[775];

    auto g_z_y_0_0_0_x_yy_xz = buffer_1100_spdd[776];

    auto g_z_y_0_0_0_x_yy_yy = buffer_1100_spdd[777];

    auto g_z_y_0_0_0_x_yy_yz = buffer_1100_spdd[778];

    auto g_z_y_0_0_0_x_yy_zz = buffer_1100_spdd[779];

    auto g_z_y_0_0_0_x_yz_xx = buffer_1100_spdd[780];

    auto g_z_y_0_0_0_x_yz_xy = buffer_1100_spdd[781];

    auto g_z_y_0_0_0_x_yz_xz = buffer_1100_spdd[782];

    auto g_z_y_0_0_0_x_yz_yy = buffer_1100_spdd[783];

    auto g_z_y_0_0_0_x_yz_yz = buffer_1100_spdd[784];

    auto g_z_y_0_0_0_x_yz_zz = buffer_1100_spdd[785];

    auto g_z_y_0_0_0_x_zz_xx = buffer_1100_spdd[786];

    auto g_z_y_0_0_0_x_zz_xy = buffer_1100_spdd[787];

    auto g_z_y_0_0_0_x_zz_xz = buffer_1100_spdd[788];

    auto g_z_y_0_0_0_x_zz_yy = buffer_1100_spdd[789];

    auto g_z_y_0_0_0_x_zz_yz = buffer_1100_spdd[790];

    auto g_z_y_0_0_0_x_zz_zz = buffer_1100_spdd[791];

    auto g_z_y_0_0_0_y_xx_xx = buffer_1100_spdd[792];

    auto g_z_y_0_0_0_y_xx_xy = buffer_1100_spdd[793];

    auto g_z_y_0_0_0_y_xx_xz = buffer_1100_spdd[794];

    auto g_z_y_0_0_0_y_xx_yy = buffer_1100_spdd[795];

    auto g_z_y_0_0_0_y_xx_yz = buffer_1100_spdd[796];

    auto g_z_y_0_0_0_y_xx_zz = buffer_1100_spdd[797];

    auto g_z_y_0_0_0_y_xy_xx = buffer_1100_spdd[798];

    auto g_z_y_0_0_0_y_xy_xy = buffer_1100_spdd[799];

    auto g_z_y_0_0_0_y_xy_xz = buffer_1100_spdd[800];

    auto g_z_y_0_0_0_y_xy_yy = buffer_1100_spdd[801];

    auto g_z_y_0_0_0_y_xy_yz = buffer_1100_spdd[802];

    auto g_z_y_0_0_0_y_xy_zz = buffer_1100_spdd[803];

    auto g_z_y_0_0_0_y_xz_xx = buffer_1100_spdd[804];

    auto g_z_y_0_0_0_y_xz_xy = buffer_1100_spdd[805];

    auto g_z_y_0_0_0_y_xz_xz = buffer_1100_spdd[806];

    auto g_z_y_0_0_0_y_xz_yy = buffer_1100_spdd[807];

    auto g_z_y_0_0_0_y_xz_yz = buffer_1100_spdd[808];

    auto g_z_y_0_0_0_y_xz_zz = buffer_1100_spdd[809];

    auto g_z_y_0_0_0_y_yy_xx = buffer_1100_spdd[810];

    auto g_z_y_0_0_0_y_yy_xy = buffer_1100_spdd[811];

    auto g_z_y_0_0_0_y_yy_xz = buffer_1100_spdd[812];

    auto g_z_y_0_0_0_y_yy_yy = buffer_1100_spdd[813];

    auto g_z_y_0_0_0_y_yy_yz = buffer_1100_spdd[814];

    auto g_z_y_0_0_0_y_yy_zz = buffer_1100_spdd[815];

    auto g_z_y_0_0_0_y_yz_xx = buffer_1100_spdd[816];

    auto g_z_y_0_0_0_y_yz_xy = buffer_1100_spdd[817];

    auto g_z_y_0_0_0_y_yz_xz = buffer_1100_spdd[818];

    auto g_z_y_0_0_0_y_yz_yy = buffer_1100_spdd[819];

    auto g_z_y_0_0_0_y_yz_yz = buffer_1100_spdd[820];

    auto g_z_y_0_0_0_y_yz_zz = buffer_1100_spdd[821];

    auto g_z_y_0_0_0_y_zz_xx = buffer_1100_spdd[822];

    auto g_z_y_0_0_0_y_zz_xy = buffer_1100_spdd[823];

    auto g_z_y_0_0_0_y_zz_xz = buffer_1100_spdd[824];

    auto g_z_y_0_0_0_y_zz_yy = buffer_1100_spdd[825];

    auto g_z_y_0_0_0_y_zz_yz = buffer_1100_spdd[826];

    auto g_z_y_0_0_0_y_zz_zz = buffer_1100_spdd[827];

    auto g_z_y_0_0_0_z_xx_xx = buffer_1100_spdd[828];

    auto g_z_y_0_0_0_z_xx_xy = buffer_1100_spdd[829];

    auto g_z_y_0_0_0_z_xx_xz = buffer_1100_spdd[830];

    auto g_z_y_0_0_0_z_xx_yy = buffer_1100_spdd[831];

    auto g_z_y_0_0_0_z_xx_yz = buffer_1100_spdd[832];

    auto g_z_y_0_0_0_z_xx_zz = buffer_1100_spdd[833];

    auto g_z_y_0_0_0_z_xy_xx = buffer_1100_spdd[834];

    auto g_z_y_0_0_0_z_xy_xy = buffer_1100_spdd[835];

    auto g_z_y_0_0_0_z_xy_xz = buffer_1100_spdd[836];

    auto g_z_y_0_0_0_z_xy_yy = buffer_1100_spdd[837];

    auto g_z_y_0_0_0_z_xy_yz = buffer_1100_spdd[838];

    auto g_z_y_0_0_0_z_xy_zz = buffer_1100_spdd[839];

    auto g_z_y_0_0_0_z_xz_xx = buffer_1100_spdd[840];

    auto g_z_y_0_0_0_z_xz_xy = buffer_1100_spdd[841];

    auto g_z_y_0_0_0_z_xz_xz = buffer_1100_spdd[842];

    auto g_z_y_0_0_0_z_xz_yy = buffer_1100_spdd[843];

    auto g_z_y_0_0_0_z_xz_yz = buffer_1100_spdd[844];

    auto g_z_y_0_0_0_z_xz_zz = buffer_1100_spdd[845];

    auto g_z_y_0_0_0_z_yy_xx = buffer_1100_spdd[846];

    auto g_z_y_0_0_0_z_yy_xy = buffer_1100_spdd[847];

    auto g_z_y_0_0_0_z_yy_xz = buffer_1100_spdd[848];

    auto g_z_y_0_0_0_z_yy_yy = buffer_1100_spdd[849];

    auto g_z_y_0_0_0_z_yy_yz = buffer_1100_spdd[850];

    auto g_z_y_0_0_0_z_yy_zz = buffer_1100_spdd[851];

    auto g_z_y_0_0_0_z_yz_xx = buffer_1100_spdd[852];

    auto g_z_y_0_0_0_z_yz_xy = buffer_1100_spdd[853];

    auto g_z_y_0_0_0_z_yz_xz = buffer_1100_spdd[854];

    auto g_z_y_0_0_0_z_yz_yy = buffer_1100_spdd[855];

    auto g_z_y_0_0_0_z_yz_yz = buffer_1100_spdd[856];

    auto g_z_y_0_0_0_z_yz_zz = buffer_1100_spdd[857];

    auto g_z_y_0_0_0_z_zz_xx = buffer_1100_spdd[858];

    auto g_z_y_0_0_0_z_zz_xy = buffer_1100_spdd[859];

    auto g_z_y_0_0_0_z_zz_xz = buffer_1100_spdd[860];

    auto g_z_y_0_0_0_z_zz_yy = buffer_1100_spdd[861];

    auto g_z_y_0_0_0_z_zz_yz = buffer_1100_spdd[862];

    auto g_z_y_0_0_0_z_zz_zz = buffer_1100_spdd[863];

    auto g_z_z_0_0_0_x_xx_xx = buffer_1100_spdd[864];

    auto g_z_z_0_0_0_x_xx_xy = buffer_1100_spdd[865];

    auto g_z_z_0_0_0_x_xx_xz = buffer_1100_spdd[866];

    auto g_z_z_0_0_0_x_xx_yy = buffer_1100_spdd[867];

    auto g_z_z_0_0_0_x_xx_yz = buffer_1100_spdd[868];

    auto g_z_z_0_0_0_x_xx_zz = buffer_1100_spdd[869];

    auto g_z_z_0_0_0_x_xy_xx = buffer_1100_spdd[870];

    auto g_z_z_0_0_0_x_xy_xy = buffer_1100_spdd[871];

    auto g_z_z_0_0_0_x_xy_xz = buffer_1100_spdd[872];

    auto g_z_z_0_0_0_x_xy_yy = buffer_1100_spdd[873];

    auto g_z_z_0_0_0_x_xy_yz = buffer_1100_spdd[874];

    auto g_z_z_0_0_0_x_xy_zz = buffer_1100_spdd[875];

    auto g_z_z_0_0_0_x_xz_xx = buffer_1100_spdd[876];

    auto g_z_z_0_0_0_x_xz_xy = buffer_1100_spdd[877];

    auto g_z_z_0_0_0_x_xz_xz = buffer_1100_spdd[878];

    auto g_z_z_0_0_0_x_xz_yy = buffer_1100_spdd[879];

    auto g_z_z_0_0_0_x_xz_yz = buffer_1100_spdd[880];

    auto g_z_z_0_0_0_x_xz_zz = buffer_1100_spdd[881];

    auto g_z_z_0_0_0_x_yy_xx = buffer_1100_spdd[882];

    auto g_z_z_0_0_0_x_yy_xy = buffer_1100_spdd[883];

    auto g_z_z_0_0_0_x_yy_xz = buffer_1100_spdd[884];

    auto g_z_z_0_0_0_x_yy_yy = buffer_1100_spdd[885];

    auto g_z_z_0_0_0_x_yy_yz = buffer_1100_spdd[886];

    auto g_z_z_0_0_0_x_yy_zz = buffer_1100_spdd[887];

    auto g_z_z_0_0_0_x_yz_xx = buffer_1100_spdd[888];

    auto g_z_z_0_0_0_x_yz_xy = buffer_1100_spdd[889];

    auto g_z_z_0_0_0_x_yz_xz = buffer_1100_spdd[890];

    auto g_z_z_0_0_0_x_yz_yy = buffer_1100_spdd[891];

    auto g_z_z_0_0_0_x_yz_yz = buffer_1100_spdd[892];

    auto g_z_z_0_0_0_x_yz_zz = buffer_1100_spdd[893];

    auto g_z_z_0_0_0_x_zz_xx = buffer_1100_spdd[894];

    auto g_z_z_0_0_0_x_zz_xy = buffer_1100_spdd[895];

    auto g_z_z_0_0_0_x_zz_xz = buffer_1100_spdd[896];

    auto g_z_z_0_0_0_x_zz_yy = buffer_1100_spdd[897];

    auto g_z_z_0_0_0_x_zz_yz = buffer_1100_spdd[898];

    auto g_z_z_0_0_0_x_zz_zz = buffer_1100_spdd[899];

    auto g_z_z_0_0_0_y_xx_xx = buffer_1100_spdd[900];

    auto g_z_z_0_0_0_y_xx_xy = buffer_1100_spdd[901];

    auto g_z_z_0_0_0_y_xx_xz = buffer_1100_spdd[902];

    auto g_z_z_0_0_0_y_xx_yy = buffer_1100_spdd[903];

    auto g_z_z_0_0_0_y_xx_yz = buffer_1100_spdd[904];

    auto g_z_z_0_0_0_y_xx_zz = buffer_1100_spdd[905];

    auto g_z_z_0_0_0_y_xy_xx = buffer_1100_spdd[906];

    auto g_z_z_0_0_0_y_xy_xy = buffer_1100_spdd[907];

    auto g_z_z_0_0_0_y_xy_xz = buffer_1100_spdd[908];

    auto g_z_z_0_0_0_y_xy_yy = buffer_1100_spdd[909];

    auto g_z_z_0_0_0_y_xy_yz = buffer_1100_spdd[910];

    auto g_z_z_0_0_0_y_xy_zz = buffer_1100_spdd[911];

    auto g_z_z_0_0_0_y_xz_xx = buffer_1100_spdd[912];

    auto g_z_z_0_0_0_y_xz_xy = buffer_1100_spdd[913];

    auto g_z_z_0_0_0_y_xz_xz = buffer_1100_spdd[914];

    auto g_z_z_0_0_0_y_xz_yy = buffer_1100_spdd[915];

    auto g_z_z_0_0_0_y_xz_yz = buffer_1100_spdd[916];

    auto g_z_z_0_0_0_y_xz_zz = buffer_1100_spdd[917];

    auto g_z_z_0_0_0_y_yy_xx = buffer_1100_spdd[918];

    auto g_z_z_0_0_0_y_yy_xy = buffer_1100_spdd[919];

    auto g_z_z_0_0_0_y_yy_xz = buffer_1100_spdd[920];

    auto g_z_z_0_0_0_y_yy_yy = buffer_1100_spdd[921];

    auto g_z_z_0_0_0_y_yy_yz = buffer_1100_spdd[922];

    auto g_z_z_0_0_0_y_yy_zz = buffer_1100_spdd[923];

    auto g_z_z_0_0_0_y_yz_xx = buffer_1100_spdd[924];

    auto g_z_z_0_0_0_y_yz_xy = buffer_1100_spdd[925];

    auto g_z_z_0_0_0_y_yz_xz = buffer_1100_spdd[926];

    auto g_z_z_0_0_0_y_yz_yy = buffer_1100_spdd[927];

    auto g_z_z_0_0_0_y_yz_yz = buffer_1100_spdd[928];

    auto g_z_z_0_0_0_y_yz_zz = buffer_1100_spdd[929];

    auto g_z_z_0_0_0_y_zz_xx = buffer_1100_spdd[930];

    auto g_z_z_0_0_0_y_zz_xy = buffer_1100_spdd[931];

    auto g_z_z_0_0_0_y_zz_xz = buffer_1100_spdd[932];

    auto g_z_z_0_0_0_y_zz_yy = buffer_1100_spdd[933];

    auto g_z_z_0_0_0_y_zz_yz = buffer_1100_spdd[934];

    auto g_z_z_0_0_0_y_zz_zz = buffer_1100_spdd[935];

    auto g_z_z_0_0_0_z_xx_xx = buffer_1100_spdd[936];

    auto g_z_z_0_0_0_z_xx_xy = buffer_1100_spdd[937];

    auto g_z_z_0_0_0_z_xx_xz = buffer_1100_spdd[938];

    auto g_z_z_0_0_0_z_xx_yy = buffer_1100_spdd[939];

    auto g_z_z_0_0_0_z_xx_yz = buffer_1100_spdd[940];

    auto g_z_z_0_0_0_z_xx_zz = buffer_1100_spdd[941];

    auto g_z_z_0_0_0_z_xy_xx = buffer_1100_spdd[942];

    auto g_z_z_0_0_0_z_xy_xy = buffer_1100_spdd[943];

    auto g_z_z_0_0_0_z_xy_xz = buffer_1100_spdd[944];

    auto g_z_z_0_0_0_z_xy_yy = buffer_1100_spdd[945];

    auto g_z_z_0_0_0_z_xy_yz = buffer_1100_spdd[946];

    auto g_z_z_0_0_0_z_xy_zz = buffer_1100_spdd[947];

    auto g_z_z_0_0_0_z_xz_xx = buffer_1100_spdd[948];

    auto g_z_z_0_0_0_z_xz_xy = buffer_1100_spdd[949];

    auto g_z_z_0_0_0_z_xz_xz = buffer_1100_spdd[950];

    auto g_z_z_0_0_0_z_xz_yy = buffer_1100_spdd[951];

    auto g_z_z_0_0_0_z_xz_yz = buffer_1100_spdd[952];

    auto g_z_z_0_0_0_z_xz_zz = buffer_1100_spdd[953];

    auto g_z_z_0_0_0_z_yy_xx = buffer_1100_spdd[954];

    auto g_z_z_0_0_0_z_yy_xy = buffer_1100_spdd[955];

    auto g_z_z_0_0_0_z_yy_xz = buffer_1100_spdd[956];

    auto g_z_z_0_0_0_z_yy_yy = buffer_1100_spdd[957];

    auto g_z_z_0_0_0_z_yy_yz = buffer_1100_spdd[958];

    auto g_z_z_0_0_0_z_yy_zz = buffer_1100_spdd[959];

    auto g_z_z_0_0_0_z_yz_xx = buffer_1100_spdd[960];

    auto g_z_z_0_0_0_z_yz_xy = buffer_1100_spdd[961];

    auto g_z_z_0_0_0_z_yz_xz = buffer_1100_spdd[962];

    auto g_z_z_0_0_0_z_yz_yy = buffer_1100_spdd[963];

    auto g_z_z_0_0_0_z_yz_yz = buffer_1100_spdd[964];

    auto g_z_z_0_0_0_z_yz_zz = buffer_1100_spdd[965];

    auto g_z_z_0_0_0_z_zz_xx = buffer_1100_spdd[966];

    auto g_z_z_0_0_0_z_zz_xy = buffer_1100_spdd[967];

    auto g_z_z_0_0_0_z_zz_xz = buffer_1100_spdd[968];

    auto g_z_z_0_0_0_z_zz_yy = buffer_1100_spdd[969];

    auto g_z_z_0_0_0_z_zz_yz = buffer_1100_spdd[970];

    auto g_z_z_0_0_0_z_zz_zz = buffer_1100_spdd[971];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_xx_xx, g_x_0_xx_xy, g_x_0_xx_xz, g_x_0_xx_yy, g_x_0_xx_yz, g_x_0_xx_zz, g_x_x_0_0_0_x_xx_xx, g_x_x_0_0_0_x_xx_xy, g_x_x_0_0_0_x_xx_xz, g_x_x_0_0_0_x_xx_yy, g_x_x_0_0_0_x_xx_yz, g_x_x_0_0_0_x_xx_zz, g_x_xx_xx_xx, g_x_xx_xx_xy, g_x_xx_xx_xz, g_x_xx_xx_yy, g_x_xx_xx_yz, g_x_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_x_xx_xx[i] = -2.0 * g_x_0_xx_xx[i] * a_exp + 4.0 * g_x_xx_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_xx_xy[i] = -2.0 * g_x_0_xx_xy[i] * a_exp + 4.0 * g_x_xx_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_xx_xz[i] = -2.0 * g_x_0_xx_xz[i] * a_exp + 4.0 * g_x_xx_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_xx_yy[i] = -2.0 * g_x_0_xx_yy[i] * a_exp + 4.0 * g_x_xx_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_xx_yz[i] = -2.0 * g_x_0_xx_yz[i] * a_exp + 4.0 * g_x_xx_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_xx_zz[i] = -2.0 * g_x_0_xx_zz[i] * a_exp + 4.0 * g_x_xx_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_xy_xx, g_x_0_xy_xy, g_x_0_xy_xz, g_x_0_xy_yy, g_x_0_xy_yz, g_x_0_xy_zz, g_x_x_0_0_0_x_xy_xx, g_x_x_0_0_0_x_xy_xy, g_x_x_0_0_0_x_xy_xz, g_x_x_0_0_0_x_xy_yy, g_x_x_0_0_0_x_xy_yz, g_x_x_0_0_0_x_xy_zz, g_x_xx_xy_xx, g_x_xx_xy_xy, g_x_xx_xy_xz, g_x_xx_xy_yy, g_x_xx_xy_yz, g_x_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_x_xy_xx[i] = -2.0 * g_x_0_xy_xx[i] * a_exp + 4.0 * g_x_xx_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_xy_xy[i] = -2.0 * g_x_0_xy_xy[i] * a_exp + 4.0 * g_x_xx_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_xy_xz[i] = -2.0 * g_x_0_xy_xz[i] * a_exp + 4.0 * g_x_xx_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_xy_yy[i] = -2.0 * g_x_0_xy_yy[i] * a_exp + 4.0 * g_x_xx_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_xy_yz[i] = -2.0 * g_x_0_xy_yz[i] * a_exp + 4.0 * g_x_xx_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_xy_zz[i] = -2.0 * g_x_0_xy_zz[i] * a_exp + 4.0 * g_x_xx_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_xz_xx, g_x_0_xz_xy, g_x_0_xz_xz, g_x_0_xz_yy, g_x_0_xz_yz, g_x_0_xz_zz, g_x_x_0_0_0_x_xz_xx, g_x_x_0_0_0_x_xz_xy, g_x_x_0_0_0_x_xz_xz, g_x_x_0_0_0_x_xz_yy, g_x_x_0_0_0_x_xz_yz, g_x_x_0_0_0_x_xz_zz, g_x_xx_xz_xx, g_x_xx_xz_xy, g_x_xx_xz_xz, g_x_xx_xz_yy, g_x_xx_xz_yz, g_x_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_x_xz_xx[i] = -2.0 * g_x_0_xz_xx[i] * a_exp + 4.0 * g_x_xx_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_xz_xy[i] = -2.0 * g_x_0_xz_xy[i] * a_exp + 4.0 * g_x_xx_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_xz_xz[i] = -2.0 * g_x_0_xz_xz[i] * a_exp + 4.0 * g_x_xx_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_xz_yy[i] = -2.0 * g_x_0_xz_yy[i] * a_exp + 4.0 * g_x_xx_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_xz_yz[i] = -2.0 * g_x_0_xz_yz[i] * a_exp + 4.0 * g_x_xx_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_xz_zz[i] = -2.0 * g_x_0_xz_zz[i] * a_exp + 4.0 * g_x_xx_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_yy_xx, g_x_0_yy_xy, g_x_0_yy_xz, g_x_0_yy_yy, g_x_0_yy_yz, g_x_0_yy_zz, g_x_x_0_0_0_x_yy_xx, g_x_x_0_0_0_x_yy_xy, g_x_x_0_0_0_x_yy_xz, g_x_x_0_0_0_x_yy_yy, g_x_x_0_0_0_x_yy_yz, g_x_x_0_0_0_x_yy_zz, g_x_xx_yy_xx, g_x_xx_yy_xy, g_x_xx_yy_xz, g_x_xx_yy_yy, g_x_xx_yy_yz, g_x_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_x_yy_xx[i] = -2.0 * g_x_0_yy_xx[i] * a_exp + 4.0 * g_x_xx_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_yy_xy[i] = -2.0 * g_x_0_yy_xy[i] * a_exp + 4.0 * g_x_xx_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_yy_xz[i] = -2.0 * g_x_0_yy_xz[i] * a_exp + 4.0 * g_x_xx_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_yy_yy[i] = -2.0 * g_x_0_yy_yy[i] * a_exp + 4.0 * g_x_xx_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_yy_yz[i] = -2.0 * g_x_0_yy_yz[i] * a_exp + 4.0 * g_x_xx_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_yy_zz[i] = -2.0 * g_x_0_yy_zz[i] * a_exp + 4.0 * g_x_xx_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_yz_xx, g_x_0_yz_xy, g_x_0_yz_xz, g_x_0_yz_yy, g_x_0_yz_yz, g_x_0_yz_zz, g_x_x_0_0_0_x_yz_xx, g_x_x_0_0_0_x_yz_xy, g_x_x_0_0_0_x_yz_xz, g_x_x_0_0_0_x_yz_yy, g_x_x_0_0_0_x_yz_yz, g_x_x_0_0_0_x_yz_zz, g_x_xx_yz_xx, g_x_xx_yz_xy, g_x_xx_yz_xz, g_x_xx_yz_yy, g_x_xx_yz_yz, g_x_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_x_yz_xx[i] = -2.0 * g_x_0_yz_xx[i] * a_exp + 4.0 * g_x_xx_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_yz_xy[i] = -2.0 * g_x_0_yz_xy[i] * a_exp + 4.0 * g_x_xx_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_yz_xz[i] = -2.0 * g_x_0_yz_xz[i] * a_exp + 4.0 * g_x_xx_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_yz_yy[i] = -2.0 * g_x_0_yz_yy[i] * a_exp + 4.0 * g_x_xx_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_yz_yz[i] = -2.0 * g_x_0_yz_yz[i] * a_exp + 4.0 * g_x_xx_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_yz_zz[i] = -2.0 * g_x_0_yz_zz[i] * a_exp + 4.0 * g_x_xx_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_zz_xx, g_x_0_zz_xy, g_x_0_zz_xz, g_x_0_zz_yy, g_x_0_zz_yz, g_x_0_zz_zz, g_x_x_0_0_0_x_zz_xx, g_x_x_0_0_0_x_zz_xy, g_x_x_0_0_0_x_zz_xz, g_x_x_0_0_0_x_zz_yy, g_x_x_0_0_0_x_zz_yz, g_x_x_0_0_0_x_zz_zz, g_x_xx_zz_xx, g_x_xx_zz_xy, g_x_xx_zz_xz, g_x_xx_zz_yy, g_x_xx_zz_yz, g_x_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_x_zz_xx[i] = -2.0 * g_x_0_zz_xx[i] * a_exp + 4.0 * g_x_xx_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_zz_xy[i] = -2.0 * g_x_0_zz_xy[i] * a_exp + 4.0 * g_x_xx_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_zz_xz[i] = -2.0 * g_x_0_zz_xz[i] * a_exp + 4.0 * g_x_xx_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_zz_yy[i] = -2.0 * g_x_0_zz_yy[i] * a_exp + 4.0 * g_x_xx_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_zz_yz[i] = -2.0 * g_x_0_zz_yz[i] * a_exp + 4.0 * g_x_xx_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_zz_zz[i] = -2.0 * g_x_0_zz_zz[i] * a_exp + 4.0 * g_x_xx_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_x_0_0_0_y_xx_xx, g_x_x_0_0_0_y_xx_xy, g_x_x_0_0_0_y_xx_xz, g_x_x_0_0_0_y_xx_yy, g_x_x_0_0_0_y_xx_yz, g_x_x_0_0_0_y_xx_zz, g_x_xy_xx_xx, g_x_xy_xx_xy, g_x_xy_xx_xz, g_x_xy_xx_yy, g_x_xy_xx_yz, g_x_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_y_xx_xx[i] = 4.0 * g_x_xy_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_xx_xy[i] = 4.0 * g_x_xy_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_xx_xz[i] = 4.0 * g_x_xy_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_xx_yy[i] = 4.0 * g_x_xy_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_xx_yz[i] = 4.0 * g_x_xy_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_xx_zz[i] = 4.0 * g_x_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_x_0_0_0_y_xy_xx, g_x_x_0_0_0_y_xy_xy, g_x_x_0_0_0_y_xy_xz, g_x_x_0_0_0_y_xy_yy, g_x_x_0_0_0_y_xy_yz, g_x_x_0_0_0_y_xy_zz, g_x_xy_xy_xx, g_x_xy_xy_xy, g_x_xy_xy_xz, g_x_xy_xy_yy, g_x_xy_xy_yz, g_x_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_y_xy_xx[i] = 4.0 * g_x_xy_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_xy_xy[i] = 4.0 * g_x_xy_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_xy_xz[i] = 4.0 * g_x_xy_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_xy_yy[i] = 4.0 * g_x_xy_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_xy_yz[i] = 4.0 * g_x_xy_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_xy_zz[i] = 4.0 * g_x_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_x_0_0_0_y_xz_xx, g_x_x_0_0_0_y_xz_xy, g_x_x_0_0_0_y_xz_xz, g_x_x_0_0_0_y_xz_yy, g_x_x_0_0_0_y_xz_yz, g_x_x_0_0_0_y_xz_zz, g_x_xy_xz_xx, g_x_xy_xz_xy, g_x_xy_xz_xz, g_x_xy_xz_yy, g_x_xy_xz_yz, g_x_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_y_xz_xx[i] = 4.0 * g_x_xy_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_xz_xy[i] = 4.0 * g_x_xy_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_xz_xz[i] = 4.0 * g_x_xy_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_xz_yy[i] = 4.0 * g_x_xy_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_xz_yz[i] = 4.0 * g_x_xy_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_xz_zz[i] = 4.0 * g_x_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_x_0_0_0_y_yy_xx, g_x_x_0_0_0_y_yy_xy, g_x_x_0_0_0_y_yy_xz, g_x_x_0_0_0_y_yy_yy, g_x_x_0_0_0_y_yy_yz, g_x_x_0_0_0_y_yy_zz, g_x_xy_yy_xx, g_x_xy_yy_xy, g_x_xy_yy_xz, g_x_xy_yy_yy, g_x_xy_yy_yz, g_x_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_y_yy_xx[i] = 4.0 * g_x_xy_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_yy_xy[i] = 4.0 * g_x_xy_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_yy_xz[i] = 4.0 * g_x_xy_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_yy_yy[i] = 4.0 * g_x_xy_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_yy_yz[i] = 4.0 * g_x_xy_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_yy_zz[i] = 4.0 * g_x_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_x_0_0_0_y_yz_xx, g_x_x_0_0_0_y_yz_xy, g_x_x_0_0_0_y_yz_xz, g_x_x_0_0_0_y_yz_yy, g_x_x_0_0_0_y_yz_yz, g_x_x_0_0_0_y_yz_zz, g_x_xy_yz_xx, g_x_xy_yz_xy, g_x_xy_yz_xz, g_x_xy_yz_yy, g_x_xy_yz_yz, g_x_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_y_yz_xx[i] = 4.0 * g_x_xy_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_yz_xy[i] = 4.0 * g_x_xy_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_yz_xz[i] = 4.0 * g_x_xy_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_yz_yy[i] = 4.0 * g_x_xy_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_yz_yz[i] = 4.0 * g_x_xy_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_yz_zz[i] = 4.0 * g_x_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_x_0_0_0_y_zz_xx, g_x_x_0_0_0_y_zz_xy, g_x_x_0_0_0_y_zz_xz, g_x_x_0_0_0_y_zz_yy, g_x_x_0_0_0_y_zz_yz, g_x_x_0_0_0_y_zz_zz, g_x_xy_zz_xx, g_x_xy_zz_xy, g_x_xy_zz_xz, g_x_xy_zz_yy, g_x_xy_zz_yz, g_x_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_y_zz_xx[i] = 4.0 * g_x_xy_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_zz_xy[i] = 4.0 * g_x_xy_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_zz_xz[i] = 4.0 * g_x_xy_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_zz_yy[i] = 4.0 * g_x_xy_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_zz_yz[i] = 4.0 * g_x_xy_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_zz_zz[i] = 4.0 * g_x_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_x_0_0_0_z_xx_xx, g_x_x_0_0_0_z_xx_xy, g_x_x_0_0_0_z_xx_xz, g_x_x_0_0_0_z_xx_yy, g_x_x_0_0_0_z_xx_yz, g_x_x_0_0_0_z_xx_zz, g_x_xz_xx_xx, g_x_xz_xx_xy, g_x_xz_xx_xz, g_x_xz_xx_yy, g_x_xz_xx_yz, g_x_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_z_xx_xx[i] = 4.0 * g_x_xz_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_xx_xy[i] = 4.0 * g_x_xz_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_xx_xz[i] = 4.0 * g_x_xz_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_xx_yy[i] = 4.0 * g_x_xz_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_xx_yz[i] = 4.0 * g_x_xz_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_xx_zz[i] = 4.0 * g_x_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_x_0_0_0_z_xy_xx, g_x_x_0_0_0_z_xy_xy, g_x_x_0_0_0_z_xy_xz, g_x_x_0_0_0_z_xy_yy, g_x_x_0_0_0_z_xy_yz, g_x_x_0_0_0_z_xy_zz, g_x_xz_xy_xx, g_x_xz_xy_xy, g_x_xz_xy_xz, g_x_xz_xy_yy, g_x_xz_xy_yz, g_x_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_z_xy_xx[i] = 4.0 * g_x_xz_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_xy_xy[i] = 4.0 * g_x_xz_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_xy_xz[i] = 4.0 * g_x_xz_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_xy_yy[i] = 4.0 * g_x_xz_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_xy_yz[i] = 4.0 * g_x_xz_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_xy_zz[i] = 4.0 * g_x_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_x_0_0_0_z_xz_xx, g_x_x_0_0_0_z_xz_xy, g_x_x_0_0_0_z_xz_xz, g_x_x_0_0_0_z_xz_yy, g_x_x_0_0_0_z_xz_yz, g_x_x_0_0_0_z_xz_zz, g_x_xz_xz_xx, g_x_xz_xz_xy, g_x_xz_xz_xz, g_x_xz_xz_yy, g_x_xz_xz_yz, g_x_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_z_xz_xx[i] = 4.0 * g_x_xz_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_xz_xy[i] = 4.0 * g_x_xz_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_xz_xz[i] = 4.0 * g_x_xz_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_xz_yy[i] = 4.0 * g_x_xz_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_xz_yz[i] = 4.0 * g_x_xz_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_xz_zz[i] = 4.0 * g_x_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_x_0_0_0_z_yy_xx, g_x_x_0_0_0_z_yy_xy, g_x_x_0_0_0_z_yy_xz, g_x_x_0_0_0_z_yy_yy, g_x_x_0_0_0_z_yy_yz, g_x_x_0_0_0_z_yy_zz, g_x_xz_yy_xx, g_x_xz_yy_xy, g_x_xz_yy_xz, g_x_xz_yy_yy, g_x_xz_yy_yz, g_x_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_z_yy_xx[i] = 4.0 * g_x_xz_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_yy_xy[i] = 4.0 * g_x_xz_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_yy_xz[i] = 4.0 * g_x_xz_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_yy_yy[i] = 4.0 * g_x_xz_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_yy_yz[i] = 4.0 * g_x_xz_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_yy_zz[i] = 4.0 * g_x_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_x_0_0_0_z_yz_xx, g_x_x_0_0_0_z_yz_xy, g_x_x_0_0_0_z_yz_xz, g_x_x_0_0_0_z_yz_yy, g_x_x_0_0_0_z_yz_yz, g_x_x_0_0_0_z_yz_zz, g_x_xz_yz_xx, g_x_xz_yz_xy, g_x_xz_yz_xz, g_x_xz_yz_yy, g_x_xz_yz_yz, g_x_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_z_yz_xx[i] = 4.0 * g_x_xz_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_yz_xy[i] = 4.0 * g_x_xz_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_yz_xz[i] = 4.0 * g_x_xz_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_yz_yy[i] = 4.0 * g_x_xz_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_yz_yz[i] = 4.0 * g_x_xz_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_yz_zz[i] = 4.0 * g_x_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_x_0_0_0_z_zz_xx, g_x_x_0_0_0_z_zz_xy, g_x_x_0_0_0_z_zz_xz, g_x_x_0_0_0_z_zz_yy, g_x_x_0_0_0_z_zz_yz, g_x_x_0_0_0_z_zz_zz, g_x_xz_zz_xx, g_x_xz_zz_xy, g_x_xz_zz_xz, g_x_xz_zz_yy, g_x_xz_zz_yz, g_x_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_z_zz_xx[i] = 4.0 * g_x_xz_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_zz_xy[i] = 4.0 * g_x_xz_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_zz_xz[i] = 4.0 * g_x_xz_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_zz_yy[i] = 4.0 * g_x_xz_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_zz_yz[i] = 4.0 * g_x_xz_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_zz_zz[i] = 4.0 * g_x_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_xy_xx_xx, g_x_xy_xx_xy, g_x_xy_xx_xz, g_x_xy_xx_yy, g_x_xy_xx_yz, g_x_xy_xx_zz, g_x_y_0_0_0_x_xx_xx, g_x_y_0_0_0_x_xx_xy, g_x_y_0_0_0_x_xx_xz, g_x_y_0_0_0_x_xx_yy, g_x_y_0_0_0_x_xx_yz, g_x_y_0_0_0_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_x_xx_xx[i] = 4.0 * g_x_xy_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_xx_xy[i] = 4.0 * g_x_xy_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_xx_xz[i] = 4.0 * g_x_xy_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_xx_yy[i] = 4.0 * g_x_xy_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_xx_yz[i] = 4.0 * g_x_xy_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_xx_zz[i] = 4.0 * g_x_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_xy_xy_xx, g_x_xy_xy_xy, g_x_xy_xy_xz, g_x_xy_xy_yy, g_x_xy_xy_yz, g_x_xy_xy_zz, g_x_y_0_0_0_x_xy_xx, g_x_y_0_0_0_x_xy_xy, g_x_y_0_0_0_x_xy_xz, g_x_y_0_0_0_x_xy_yy, g_x_y_0_0_0_x_xy_yz, g_x_y_0_0_0_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_x_xy_xx[i] = 4.0 * g_x_xy_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_xy_xy[i] = 4.0 * g_x_xy_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_xy_xz[i] = 4.0 * g_x_xy_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_xy_yy[i] = 4.0 * g_x_xy_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_xy_yz[i] = 4.0 * g_x_xy_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_xy_zz[i] = 4.0 * g_x_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_xy_xz_xx, g_x_xy_xz_xy, g_x_xy_xz_xz, g_x_xy_xz_yy, g_x_xy_xz_yz, g_x_xy_xz_zz, g_x_y_0_0_0_x_xz_xx, g_x_y_0_0_0_x_xz_xy, g_x_y_0_0_0_x_xz_xz, g_x_y_0_0_0_x_xz_yy, g_x_y_0_0_0_x_xz_yz, g_x_y_0_0_0_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_x_xz_xx[i] = 4.0 * g_x_xy_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_xz_xy[i] = 4.0 * g_x_xy_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_xz_xz[i] = 4.0 * g_x_xy_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_xz_yy[i] = 4.0 * g_x_xy_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_xz_yz[i] = 4.0 * g_x_xy_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_xz_zz[i] = 4.0 * g_x_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_xy_yy_xx, g_x_xy_yy_xy, g_x_xy_yy_xz, g_x_xy_yy_yy, g_x_xy_yy_yz, g_x_xy_yy_zz, g_x_y_0_0_0_x_yy_xx, g_x_y_0_0_0_x_yy_xy, g_x_y_0_0_0_x_yy_xz, g_x_y_0_0_0_x_yy_yy, g_x_y_0_0_0_x_yy_yz, g_x_y_0_0_0_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_x_yy_xx[i] = 4.0 * g_x_xy_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_yy_xy[i] = 4.0 * g_x_xy_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_yy_xz[i] = 4.0 * g_x_xy_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_yy_yy[i] = 4.0 * g_x_xy_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_yy_yz[i] = 4.0 * g_x_xy_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_yy_zz[i] = 4.0 * g_x_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_xy_yz_xx, g_x_xy_yz_xy, g_x_xy_yz_xz, g_x_xy_yz_yy, g_x_xy_yz_yz, g_x_xy_yz_zz, g_x_y_0_0_0_x_yz_xx, g_x_y_0_0_0_x_yz_xy, g_x_y_0_0_0_x_yz_xz, g_x_y_0_0_0_x_yz_yy, g_x_y_0_0_0_x_yz_yz, g_x_y_0_0_0_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_x_yz_xx[i] = 4.0 * g_x_xy_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_yz_xy[i] = 4.0 * g_x_xy_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_yz_xz[i] = 4.0 * g_x_xy_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_yz_yy[i] = 4.0 * g_x_xy_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_yz_yz[i] = 4.0 * g_x_xy_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_yz_zz[i] = 4.0 * g_x_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_xy_zz_xx, g_x_xy_zz_xy, g_x_xy_zz_xz, g_x_xy_zz_yy, g_x_xy_zz_yz, g_x_xy_zz_zz, g_x_y_0_0_0_x_zz_xx, g_x_y_0_0_0_x_zz_xy, g_x_y_0_0_0_x_zz_xz, g_x_y_0_0_0_x_zz_yy, g_x_y_0_0_0_x_zz_yz, g_x_y_0_0_0_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_x_zz_xx[i] = 4.0 * g_x_xy_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_zz_xy[i] = 4.0 * g_x_xy_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_zz_xz[i] = 4.0 * g_x_xy_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_zz_yy[i] = 4.0 * g_x_xy_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_zz_yz[i] = 4.0 * g_x_xy_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_zz_zz[i] = 4.0 * g_x_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_xx_xx, g_x_0_xx_xy, g_x_0_xx_xz, g_x_0_xx_yy, g_x_0_xx_yz, g_x_0_xx_zz, g_x_y_0_0_0_y_xx_xx, g_x_y_0_0_0_y_xx_xy, g_x_y_0_0_0_y_xx_xz, g_x_y_0_0_0_y_xx_yy, g_x_y_0_0_0_y_xx_yz, g_x_y_0_0_0_y_xx_zz, g_x_yy_xx_xx, g_x_yy_xx_xy, g_x_yy_xx_xz, g_x_yy_xx_yy, g_x_yy_xx_yz, g_x_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_y_xx_xx[i] = -2.0 * g_x_0_xx_xx[i] * a_exp + 4.0 * g_x_yy_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_xx_xy[i] = -2.0 * g_x_0_xx_xy[i] * a_exp + 4.0 * g_x_yy_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_xx_xz[i] = -2.0 * g_x_0_xx_xz[i] * a_exp + 4.0 * g_x_yy_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_xx_yy[i] = -2.0 * g_x_0_xx_yy[i] * a_exp + 4.0 * g_x_yy_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_xx_yz[i] = -2.0 * g_x_0_xx_yz[i] * a_exp + 4.0 * g_x_yy_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_xx_zz[i] = -2.0 * g_x_0_xx_zz[i] * a_exp + 4.0 * g_x_yy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_xy_xx, g_x_0_xy_xy, g_x_0_xy_xz, g_x_0_xy_yy, g_x_0_xy_yz, g_x_0_xy_zz, g_x_y_0_0_0_y_xy_xx, g_x_y_0_0_0_y_xy_xy, g_x_y_0_0_0_y_xy_xz, g_x_y_0_0_0_y_xy_yy, g_x_y_0_0_0_y_xy_yz, g_x_y_0_0_0_y_xy_zz, g_x_yy_xy_xx, g_x_yy_xy_xy, g_x_yy_xy_xz, g_x_yy_xy_yy, g_x_yy_xy_yz, g_x_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_y_xy_xx[i] = -2.0 * g_x_0_xy_xx[i] * a_exp + 4.0 * g_x_yy_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_xy_xy[i] = -2.0 * g_x_0_xy_xy[i] * a_exp + 4.0 * g_x_yy_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_xy_xz[i] = -2.0 * g_x_0_xy_xz[i] * a_exp + 4.0 * g_x_yy_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_xy_yy[i] = -2.0 * g_x_0_xy_yy[i] * a_exp + 4.0 * g_x_yy_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_xy_yz[i] = -2.0 * g_x_0_xy_yz[i] * a_exp + 4.0 * g_x_yy_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_xy_zz[i] = -2.0 * g_x_0_xy_zz[i] * a_exp + 4.0 * g_x_yy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_xz_xx, g_x_0_xz_xy, g_x_0_xz_xz, g_x_0_xz_yy, g_x_0_xz_yz, g_x_0_xz_zz, g_x_y_0_0_0_y_xz_xx, g_x_y_0_0_0_y_xz_xy, g_x_y_0_0_0_y_xz_xz, g_x_y_0_0_0_y_xz_yy, g_x_y_0_0_0_y_xz_yz, g_x_y_0_0_0_y_xz_zz, g_x_yy_xz_xx, g_x_yy_xz_xy, g_x_yy_xz_xz, g_x_yy_xz_yy, g_x_yy_xz_yz, g_x_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_y_xz_xx[i] = -2.0 * g_x_0_xz_xx[i] * a_exp + 4.0 * g_x_yy_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_xz_xy[i] = -2.0 * g_x_0_xz_xy[i] * a_exp + 4.0 * g_x_yy_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_xz_xz[i] = -2.0 * g_x_0_xz_xz[i] * a_exp + 4.0 * g_x_yy_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_xz_yy[i] = -2.0 * g_x_0_xz_yy[i] * a_exp + 4.0 * g_x_yy_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_xz_yz[i] = -2.0 * g_x_0_xz_yz[i] * a_exp + 4.0 * g_x_yy_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_xz_zz[i] = -2.0 * g_x_0_xz_zz[i] * a_exp + 4.0 * g_x_yy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_0_yy_xx, g_x_0_yy_xy, g_x_0_yy_xz, g_x_0_yy_yy, g_x_0_yy_yz, g_x_0_yy_zz, g_x_y_0_0_0_y_yy_xx, g_x_y_0_0_0_y_yy_xy, g_x_y_0_0_0_y_yy_xz, g_x_y_0_0_0_y_yy_yy, g_x_y_0_0_0_y_yy_yz, g_x_y_0_0_0_y_yy_zz, g_x_yy_yy_xx, g_x_yy_yy_xy, g_x_yy_yy_xz, g_x_yy_yy_yy, g_x_yy_yy_yz, g_x_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_y_yy_xx[i] = -2.0 * g_x_0_yy_xx[i] * a_exp + 4.0 * g_x_yy_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_yy_xy[i] = -2.0 * g_x_0_yy_xy[i] * a_exp + 4.0 * g_x_yy_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_yy_xz[i] = -2.0 * g_x_0_yy_xz[i] * a_exp + 4.0 * g_x_yy_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_yy_yy[i] = -2.0 * g_x_0_yy_yy[i] * a_exp + 4.0 * g_x_yy_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_yy_yz[i] = -2.0 * g_x_0_yy_yz[i] * a_exp + 4.0 * g_x_yy_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_yy_zz[i] = -2.0 * g_x_0_yy_zz[i] * a_exp + 4.0 * g_x_yy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_0_yz_xx, g_x_0_yz_xy, g_x_0_yz_xz, g_x_0_yz_yy, g_x_0_yz_yz, g_x_0_yz_zz, g_x_y_0_0_0_y_yz_xx, g_x_y_0_0_0_y_yz_xy, g_x_y_0_0_0_y_yz_xz, g_x_y_0_0_0_y_yz_yy, g_x_y_0_0_0_y_yz_yz, g_x_y_0_0_0_y_yz_zz, g_x_yy_yz_xx, g_x_yy_yz_xy, g_x_yy_yz_xz, g_x_yy_yz_yy, g_x_yy_yz_yz, g_x_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_y_yz_xx[i] = -2.0 * g_x_0_yz_xx[i] * a_exp + 4.0 * g_x_yy_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_yz_xy[i] = -2.0 * g_x_0_yz_xy[i] * a_exp + 4.0 * g_x_yy_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_yz_xz[i] = -2.0 * g_x_0_yz_xz[i] * a_exp + 4.0 * g_x_yy_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_yz_yy[i] = -2.0 * g_x_0_yz_yy[i] * a_exp + 4.0 * g_x_yy_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_yz_yz[i] = -2.0 * g_x_0_yz_yz[i] * a_exp + 4.0 * g_x_yy_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_yz_zz[i] = -2.0 * g_x_0_yz_zz[i] * a_exp + 4.0 * g_x_yy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_0_zz_xx, g_x_0_zz_xy, g_x_0_zz_xz, g_x_0_zz_yy, g_x_0_zz_yz, g_x_0_zz_zz, g_x_y_0_0_0_y_zz_xx, g_x_y_0_0_0_y_zz_xy, g_x_y_0_0_0_y_zz_xz, g_x_y_0_0_0_y_zz_yy, g_x_y_0_0_0_y_zz_yz, g_x_y_0_0_0_y_zz_zz, g_x_yy_zz_xx, g_x_yy_zz_xy, g_x_yy_zz_xz, g_x_yy_zz_yy, g_x_yy_zz_yz, g_x_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_y_zz_xx[i] = -2.0 * g_x_0_zz_xx[i] * a_exp + 4.0 * g_x_yy_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_zz_xy[i] = -2.0 * g_x_0_zz_xy[i] * a_exp + 4.0 * g_x_yy_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_zz_xz[i] = -2.0 * g_x_0_zz_xz[i] * a_exp + 4.0 * g_x_yy_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_zz_yy[i] = -2.0 * g_x_0_zz_yy[i] * a_exp + 4.0 * g_x_yy_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_zz_yz[i] = -2.0 * g_x_0_zz_yz[i] * a_exp + 4.0 * g_x_yy_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_zz_zz[i] = -2.0 * g_x_0_zz_zz[i] * a_exp + 4.0 * g_x_yy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_y_0_0_0_z_xx_xx, g_x_y_0_0_0_z_xx_xy, g_x_y_0_0_0_z_xx_xz, g_x_y_0_0_0_z_xx_yy, g_x_y_0_0_0_z_xx_yz, g_x_y_0_0_0_z_xx_zz, g_x_yz_xx_xx, g_x_yz_xx_xy, g_x_yz_xx_xz, g_x_yz_xx_yy, g_x_yz_xx_yz, g_x_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_z_xx_xx[i] = 4.0 * g_x_yz_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_xx_xy[i] = 4.0 * g_x_yz_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_xx_xz[i] = 4.0 * g_x_yz_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_xx_yy[i] = 4.0 * g_x_yz_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_xx_yz[i] = 4.0 * g_x_yz_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_xx_zz[i] = 4.0 * g_x_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_y_0_0_0_z_xy_xx, g_x_y_0_0_0_z_xy_xy, g_x_y_0_0_0_z_xy_xz, g_x_y_0_0_0_z_xy_yy, g_x_y_0_0_0_z_xy_yz, g_x_y_0_0_0_z_xy_zz, g_x_yz_xy_xx, g_x_yz_xy_xy, g_x_yz_xy_xz, g_x_yz_xy_yy, g_x_yz_xy_yz, g_x_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_z_xy_xx[i] = 4.0 * g_x_yz_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_xy_xy[i] = 4.0 * g_x_yz_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_xy_xz[i] = 4.0 * g_x_yz_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_xy_yy[i] = 4.0 * g_x_yz_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_xy_yz[i] = 4.0 * g_x_yz_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_xy_zz[i] = 4.0 * g_x_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_y_0_0_0_z_xz_xx, g_x_y_0_0_0_z_xz_xy, g_x_y_0_0_0_z_xz_xz, g_x_y_0_0_0_z_xz_yy, g_x_y_0_0_0_z_xz_yz, g_x_y_0_0_0_z_xz_zz, g_x_yz_xz_xx, g_x_yz_xz_xy, g_x_yz_xz_xz, g_x_yz_xz_yy, g_x_yz_xz_yz, g_x_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_z_xz_xx[i] = 4.0 * g_x_yz_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_xz_xy[i] = 4.0 * g_x_yz_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_xz_xz[i] = 4.0 * g_x_yz_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_xz_yy[i] = 4.0 * g_x_yz_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_xz_yz[i] = 4.0 * g_x_yz_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_xz_zz[i] = 4.0 * g_x_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_y_0_0_0_z_yy_xx, g_x_y_0_0_0_z_yy_xy, g_x_y_0_0_0_z_yy_xz, g_x_y_0_0_0_z_yy_yy, g_x_y_0_0_0_z_yy_yz, g_x_y_0_0_0_z_yy_zz, g_x_yz_yy_xx, g_x_yz_yy_xy, g_x_yz_yy_xz, g_x_yz_yy_yy, g_x_yz_yy_yz, g_x_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_z_yy_xx[i] = 4.0 * g_x_yz_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_yy_xy[i] = 4.0 * g_x_yz_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_yy_xz[i] = 4.0 * g_x_yz_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_yy_yy[i] = 4.0 * g_x_yz_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_yy_yz[i] = 4.0 * g_x_yz_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_yy_zz[i] = 4.0 * g_x_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_y_0_0_0_z_yz_xx, g_x_y_0_0_0_z_yz_xy, g_x_y_0_0_0_z_yz_xz, g_x_y_0_0_0_z_yz_yy, g_x_y_0_0_0_z_yz_yz, g_x_y_0_0_0_z_yz_zz, g_x_yz_yz_xx, g_x_yz_yz_xy, g_x_yz_yz_xz, g_x_yz_yz_yy, g_x_yz_yz_yz, g_x_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_z_yz_xx[i] = 4.0 * g_x_yz_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_yz_xy[i] = 4.0 * g_x_yz_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_yz_xz[i] = 4.0 * g_x_yz_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_yz_yy[i] = 4.0 * g_x_yz_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_yz_yz[i] = 4.0 * g_x_yz_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_yz_zz[i] = 4.0 * g_x_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_y_0_0_0_z_zz_xx, g_x_y_0_0_0_z_zz_xy, g_x_y_0_0_0_z_zz_xz, g_x_y_0_0_0_z_zz_yy, g_x_y_0_0_0_z_zz_yz, g_x_y_0_0_0_z_zz_zz, g_x_yz_zz_xx, g_x_yz_zz_xy, g_x_yz_zz_xz, g_x_yz_zz_yy, g_x_yz_zz_yz, g_x_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_z_zz_xx[i] = 4.0 * g_x_yz_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_zz_xy[i] = 4.0 * g_x_yz_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_zz_xz[i] = 4.0 * g_x_yz_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_zz_yy[i] = 4.0 * g_x_yz_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_zz_yz[i] = 4.0 * g_x_yz_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_zz_zz[i] = 4.0 * g_x_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_xz_xx_xx, g_x_xz_xx_xy, g_x_xz_xx_xz, g_x_xz_xx_yy, g_x_xz_xx_yz, g_x_xz_xx_zz, g_x_z_0_0_0_x_xx_xx, g_x_z_0_0_0_x_xx_xy, g_x_z_0_0_0_x_xx_xz, g_x_z_0_0_0_x_xx_yy, g_x_z_0_0_0_x_xx_yz, g_x_z_0_0_0_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_x_xx_xx[i] = 4.0 * g_x_xz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_xx_xy[i] = 4.0 * g_x_xz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_xx_xz[i] = 4.0 * g_x_xz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_xx_yy[i] = 4.0 * g_x_xz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_xx_yz[i] = 4.0 * g_x_xz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_xx_zz[i] = 4.0 * g_x_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_xz_xy_xx, g_x_xz_xy_xy, g_x_xz_xy_xz, g_x_xz_xy_yy, g_x_xz_xy_yz, g_x_xz_xy_zz, g_x_z_0_0_0_x_xy_xx, g_x_z_0_0_0_x_xy_xy, g_x_z_0_0_0_x_xy_xz, g_x_z_0_0_0_x_xy_yy, g_x_z_0_0_0_x_xy_yz, g_x_z_0_0_0_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_x_xy_xx[i] = 4.0 * g_x_xz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_xy_xy[i] = 4.0 * g_x_xz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_xy_xz[i] = 4.0 * g_x_xz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_xy_yy[i] = 4.0 * g_x_xz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_xy_yz[i] = 4.0 * g_x_xz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_xy_zz[i] = 4.0 * g_x_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_xz_xz_xx, g_x_xz_xz_xy, g_x_xz_xz_xz, g_x_xz_xz_yy, g_x_xz_xz_yz, g_x_xz_xz_zz, g_x_z_0_0_0_x_xz_xx, g_x_z_0_0_0_x_xz_xy, g_x_z_0_0_0_x_xz_xz, g_x_z_0_0_0_x_xz_yy, g_x_z_0_0_0_x_xz_yz, g_x_z_0_0_0_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_x_xz_xx[i] = 4.0 * g_x_xz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_xz_xy[i] = 4.0 * g_x_xz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_xz_xz[i] = 4.0 * g_x_xz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_xz_yy[i] = 4.0 * g_x_xz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_xz_yz[i] = 4.0 * g_x_xz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_xz_zz[i] = 4.0 * g_x_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_xz_yy_xx, g_x_xz_yy_xy, g_x_xz_yy_xz, g_x_xz_yy_yy, g_x_xz_yy_yz, g_x_xz_yy_zz, g_x_z_0_0_0_x_yy_xx, g_x_z_0_0_0_x_yy_xy, g_x_z_0_0_0_x_yy_xz, g_x_z_0_0_0_x_yy_yy, g_x_z_0_0_0_x_yy_yz, g_x_z_0_0_0_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_x_yy_xx[i] = 4.0 * g_x_xz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_yy_xy[i] = 4.0 * g_x_xz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_yy_xz[i] = 4.0 * g_x_xz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_yy_yy[i] = 4.0 * g_x_xz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_yy_yz[i] = 4.0 * g_x_xz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_yy_zz[i] = 4.0 * g_x_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_xz_yz_xx, g_x_xz_yz_xy, g_x_xz_yz_xz, g_x_xz_yz_yy, g_x_xz_yz_yz, g_x_xz_yz_zz, g_x_z_0_0_0_x_yz_xx, g_x_z_0_0_0_x_yz_xy, g_x_z_0_0_0_x_yz_xz, g_x_z_0_0_0_x_yz_yy, g_x_z_0_0_0_x_yz_yz, g_x_z_0_0_0_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_x_yz_xx[i] = 4.0 * g_x_xz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_yz_xy[i] = 4.0 * g_x_xz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_yz_xz[i] = 4.0 * g_x_xz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_yz_yy[i] = 4.0 * g_x_xz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_yz_yz[i] = 4.0 * g_x_xz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_yz_zz[i] = 4.0 * g_x_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_xz_zz_xx, g_x_xz_zz_xy, g_x_xz_zz_xz, g_x_xz_zz_yy, g_x_xz_zz_yz, g_x_xz_zz_zz, g_x_z_0_0_0_x_zz_xx, g_x_z_0_0_0_x_zz_xy, g_x_z_0_0_0_x_zz_xz, g_x_z_0_0_0_x_zz_yy, g_x_z_0_0_0_x_zz_yz, g_x_z_0_0_0_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_x_zz_xx[i] = 4.0 * g_x_xz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_zz_xy[i] = 4.0 * g_x_xz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_zz_xz[i] = 4.0 * g_x_xz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_zz_yy[i] = 4.0 * g_x_xz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_zz_yz[i] = 4.0 * g_x_xz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_zz_zz[i] = 4.0 * g_x_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_yz_xx_xx, g_x_yz_xx_xy, g_x_yz_xx_xz, g_x_yz_xx_yy, g_x_yz_xx_yz, g_x_yz_xx_zz, g_x_z_0_0_0_y_xx_xx, g_x_z_0_0_0_y_xx_xy, g_x_z_0_0_0_y_xx_xz, g_x_z_0_0_0_y_xx_yy, g_x_z_0_0_0_y_xx_yz, g_x_z_0_0_0_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_y_xx_xx[i] = 4.0 * g_x_yz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_xx_xy[i] = 4.0 * g_x_yz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_xx_xz[i] = 4.0 * g_x_yz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_xx_yy[i] = 4.0 * g_x_yz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_xx_yz[i] = 4.0 * g_x_yz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_xx_zz[i] = 4.0 * g_x_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_yz_xy_xx, g_x_yz_xy_xy, g_x_yz_xy_xz, g_x_yz_xy_yy, g_x_yz_xy_yz, g_x_yz_xy_zz, g_x_z_0_0_0_y_xy_xx, g_x_z_0_0_0_y_xy_xy, g_x_z_0_0_0_y_xy_xz, g_x_z_0_0_0_y_xy_yy, g_x_z_0_0_0_y_xy_yz, g_x_z_0_0_0_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_y_xy_xx[i] = 4.0 * g_x_yz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_xy_xy[i] = 4.0 * g_x_yz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_xy_xz[i] = 4.0 * g_x_yz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_xy_yy[i] = 4.0 * g_x_yz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_xy_yz[i] = 4.0 * g_x_yz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_xy_zz[i] = 4.0 * g_x_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_yz_xz_xx, g_x_yz_xz_xy, g_x_yz_xz_xz, g_x_yz_xz_yy, g_x_yz_xz_yz, g_x_yz_xz_zz, g_x_z_0_0_0_y_xz_xx, g_x_z_0_0_0_y_xz_xy, g_x_z_0_0_0_y_xz_xz, g_x_z_0_0_0_y_xz_yy, g_x_z_0_0_0_y_xz_yz, g_x_z_0_0_0_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_y_xz_xx[i] = 4.0 * g_x_yz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_xz_xy[i] = 4.0 * g_x_yz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_xz_xz[i] = 4.0 * g_x_yz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_xz_yy[i] = 4.0 * g_x_yz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_xz_yz[i] = 4.0 * g_x_yz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_xz_zz[i] = 4.0 * g_x_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_yz_yy_xx, g_x_yz_yy_xy, g_x_yz_yy_xz, g_x_yz_yy_yy, g_x_yz_yy_yz, g_x_yz_yy_zz, g_x_z_0_0_0_y_yy_xx, g_x_z_0_0_0_y_yy_xy, g_x_z_0_0_0_y_yy_xz, g_x_z_0_0_0_y_yy_yy, g_x_z_0_0_0_y_yy_yz, g_x_z_0_0_0_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_y_yy_xx[i] = 4.0 * g_x_yz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_yy_xy[i] = 4.0 * g_x_yz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_yy_xz[i] = 4.0 * g_x_yz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_yy_yy[i] = 4.0 * g_x_yz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_yy_yz[i] = 4.0 * g_x_yz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_yy_zz[i] = 4.0 * g_x_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_yz_yz_xx, g_x_yz_yz_xy, g_x_yz_yz_xz, g_x_yz_yz_yy, g_x_yz_yz_yz, g_x_yz_yz_zz, g_x_z_0_0_0_y_yz_xx, g_x_z_0_0_0_y_yz_xy, g_x_z_0_0_0_y_yz_xz, g_x_z_0_0_0_y_yz_yy, g_x_z_0_0_0_y_yz_yz, g_x_z_0_0_0_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_y_yz_xx[i] = 4.0 * g_x_yz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_yz_xy[i] = 4.0 * g_x_yz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_yz_xz[i] = 4.0 * g_x_yz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_yz_yy[i] = 4.0 * g_x_yz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_yz_yz[i] = 4.0 * g_x_yz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_yz_zz[i] = 4.0 * g_x_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_yz_zz_xx, g_x_yz_zz_xy, g_x_yz_zz_xz, g_x_yz_zz_yy, g_x_yz_zz_yz, g_x_yz_zz_zz, g_x_z_0_0_0_y_zz_xx, g_x_z_0_0_0_y_zz_xy, g_x_z_0_0_0_y_zz_xz, g_x_z_0_0_0_y_zz_yy, g_x_z_0_0_0_y_zz_yz, g_x_z_0_0_0_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_y_zz_xx[i] = 4.0 * g_x_yz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_zz_xy[i] = 4.0 * g_x_yz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_zz_xz[i] = 4.0 * g_x_yz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_zz_yy[i] = 4.0 * g_x_yz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_zz_yz[i] = 4.0 * g_x_yz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_zz_zz[i] = 4.0 * g_x_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_0_xx_xx, g_x_0_xx_xy, g_x_0_xx_xz, g_x_0_xx_yy, g_x_0_xx_yz, g_x_0_xx_zz, g_x_z_0_0_0_z_xx_xx, g_x_z_0_0_0_z_xx_xy, g_x_z_0_0_0_z_xx_xz, g_x_z_0_0_0_z_xx_yy, g_x_z_0_0_0_z_xx_yz, g_x_z_0_0_0_z_xx_zz, g_x_zz_xx_xx, g_x_zz_xx_xy, g_x_zz_xx_xz, g_x_zz_xx_yy, g_x_zz_xx_yz, g_x_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_z_xx_xx[i] = -2.0 * g_x_0_xx_xx[i] * a_exp + 4.0 * g_x_zz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_xx_xy[i] = -2.0 * g_x_0_xx_xy[i] * a_exp + 4.0 * g_x_zz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_xx_xz[i] = -2.0 * g_x_0_xx_xz[i] * a_exp + 4.0 * g_x_zz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_xx_yy[i] = -2.0 * g_x_0_xx_yy[i] * a_exp + 4.0 * g_x_zz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_xx_yz[i] = -2.0 * g_x_0_xx_yz[i] * a_exp + 4.0 * g_x_zz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_xx_zz[i] = -2.0 * g_x_0_xx_zz[i] * a_exp + 4.0 * g_x_zz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_0_xy_xx, g_x_0_xy_xy, g_x_0_xy_xz, g_x_0_xy_yy, g_x_0_xy_yz, g_x_0_xy_zz, g_x_z_0_0_0_z_xy_xx, g_x_z_0_0_0_z_xy_xy, g_x_z_0_0_0_z_xy_xz, g_x_z_0_0_0_z_xy_yy, g_x_z_0_0_0_z_xy_yz, g_x_z_0_0_0_z_xy_zz, g_x_zz_xy_xx, g_x_zz_xy_xy, g_x_zz_xy_xz, g_x_zz_xy_yy, g_x_zz_xy_yz, g_x_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_z_xy_xx[i] = -2.0 * g_x_0_xy_xx[i] * a_exp + 4.0 * g_x_zz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_xy_xy[i] = -2.0 * g_x_0_xy_xy[i] * a_exp + 4.0 * g_x_zz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_xy_xz[i] = -2.0 * g_x_0_xy_xz[i] * a_exp + 4.0 * g_x_zz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_xy_yy[i] = -2.0 * g_x_0_xy_yy[i] * a_exp + 4.0 * g_x_zz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_xy_yz[i] = -2.0 * g_x_0_xy_yz[i] * a_exp + 4.0 * g_x_zz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_xy_zz[i] = -2.0 * g_x_0_xy_zz[i] * a_exp + 4.0 * g_x_zz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_0_xz_xx, g_x_0_xz_xy, g_x_0_xz_xz, g_x_0_xz_yy, g_x_0_xz_yz, g_x_0_xz_zz, g_x_z_0_0_0_z_xz_xx, g_x_z_0_0_0_z_xz_xy, g_x_z_0_0_0_z_xz_xz, g_x_z_0_0_0_z_xz_yy, g_x_z_0_0_0_z_xz_yz, g_x_z_0_0_0_z_xz_zz, g_x_zz_xz_xx, g_x_zz_xz_xy, g_x_zz_xz_xz, g_x_zz_xz_yy, g_x_zz_xz_yz, g_x_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_z_xz_xx[i] = -2.0 * g_x_0_xz_xx[i] * a_exp + 4.0 * g_x_zz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_xz_xy[i] = -2.0 * g_x_0_xz_xy[i] * a_exp + 4.0 * g_x_zz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_xz_xz[i] = -2.0 * g_x_0_xz_xz[i] * a_exp + 4.0 * g_x_zz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_xz_yy[i] = -2.0 * g_x_0_xz_yy[i] * a_exp + 4.0 * g_x_zz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_xz_yz[i] = -2.0 * g_x_0_xz_yz[i] * a_exp + 4.0 * g_x_zz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_xz_zz[i] = -2.0 * g_x_0_xz_zz[i] * a_exp + 4.0 * g_x_zz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_0_yy_xx, g_x_0_yy_xy, g_x_0_yy_xz, g_x_0_yy_yy, g_x_0_yy_yz, g_x_0_yy_zz, g_x_z_0_0_0_z_yy_xx, g_x_z_0_0_0_z_yy_xy, g_x_z_0_0_0_z_yy_xz, g_x_z_0_0_0_z_yy_yy, g_x_z_0_0_0_z_yy_yz, g_x_z_0_0_0_z_yy_zz, g_x_zz_yy_xx, g_x_zz_yy_xy, g_x_zz_yy_xz, g_x_zz_yy_yy, g_x_zz_yy_yz, g_x_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_z_yy_xx[i] = -2.0 * g_x_0_yy_xx[i] * a_exp + 4.0 * g_x_zz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_yy_xy[i] = -2.0 * g_x_0_yy_xy[i] * a_exp + 4.0 * g_x_zz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_yy_xz[i] = -2.0 * g_x_0_yy_xz[i] * a_exp + 4.0 * g_x_zz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_yy_yy[i] = -2.0 * g_x_0_yy_yy[i] * a_exp + 4.0 * g_x_zz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_yy_yz[i] = -2.0 * g_x_0_yy_yz[i] * a_exp + 4.0 * g_x_zz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_yy_zz[i] = -2.0 * g_x_0_yy_zz[i] * a_exp + 4.0 * g_x_zz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_0_yz_xx, g_x_0_yz_xy, g_x_0_yz_xz, g_x_0_yz_yy, g_x_0_yz_yz, g_x_0_yz_zz, g_x_z_0_0_0_z_yz_xx, g_x_z_0_0_0_z_yz_xy, g_x_z_0_0_0_z_yz_xz, g_x_z_0_0_0_z_yz_yy, g_x_z_0_0_0_z_yz_yz, g_x_z_0_0_0_z_yz_zz, g_x_zz_yz_xx, g_x_zz_yz_xy, g_x_zz_yz_xz, g_x_zz_yz_yy, g_x_zz_yz_yz, g_x_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_z_yz_xx[i] = -2.0 * g_x_0_yz_xx[i] * a_exp + 4.0 * g_x_zz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_yz_xy[i] = -2.0 * g_x_0_yz_xy[i] * a_exp + 4.0 * g_x_zz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_yz_xz[i] = -2.0 * g_x_0_yz_xz[i] * a_exp + 4.0 * g_x_zz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_yz_yy[i] = -2.0 * g_x_0_yz_yy[i] * a_exp + 4.0 * g_x_zz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_yz_yz[i] = -2.0 * g_x_0_yz_yz[i] * a_exp + 4.0 * g_x_zz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_yz_zz[i] = -2.0 * g_x_0_yz_zz[i] * a_exp + 4.0 * g_x_zz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_0_zz_xx, g_x_0_zz_xy, g_x_0_zz_xz, g_x_0_zz_yy, g_x_0_zz_yz, g_x_0_zz_zz, g_x_z_0_0_0_z_zz_xx, g_x_z_0_0_0_z_zz_xy, g_x_z_0_0_0_z_zz_xz, g_x_z_0_0_0_z_zz_yy, g_x_z_0_0_0_z_zz_yz, g_x_z_0_0_0_z_zz_zz, g_x_zz_zz_xx, g_x_zz_zz_xy, g_x_zz_zz_xz, g_x_zz_zz_yy, g_x_zz_zz_yz, g_x_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_z_zz_xx[i] = -2.0 * g_x_0_zz_xx[i] * a_exp + 4.0 * g_x_zz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_zz_xy[i] = -2.0 * g_x_0_zz_xy[i] * a_exp + 4.0 * g_x_zz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_zz_xz[i] = -2.0 * g_x_0_zz_xz[i] * a_exp + 4.0 * g_x_zz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_zz_yy[i] = -2.0 * g_x_0_zz_yy[i] * a_exp + 4.0 * g_x_zz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_zz_yz[i] = -2.0 * g_x_0_zz_yz[i] * a_exp + 4.0 * g_x_zz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_zz_zz[i] = -2.0 * g_x_0_zz_zz[i] * a_exp + 4.0 * g_x_zz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_y_0_xx_xx, g_y_0_xx_xy, g_y_0_xx_xz, g_y_0_xx_yy, g_y_0_xx_yz, g_y_0_xx_zz, g_y_x_0_0_0_x_xx_xx, g_y_x_0_0_0_x_xx_xy, g_y_x_0_0_0_x_xx_xz, g_y_x_0_0_0_x_xx_yy, g_y_x_0_0_0_x_xx_yz, g_y_x_0_0_0_x_xx_zz, g_y_xx_xx_xx, g_y_xx_xx_xy, g_y_xx_xx_xz, g_y_xx_xx_yy, g_y_xx_xx_yz, g_y_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_x_xx_xx[i] = -2.0 * g_y_0_xx_xx[i] * a_exp + 4.0 * g_y_xx_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_xx_xy[i] = -2.0 * g_y_0_xx_xy[i] * a_exp + 4.0 * g_y_xx_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_xx_xz[i] = -2.0 * g_y_0_xx_xz[i] * a_exp + 4.0 * g_y_xx_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_xx_yy[i] = -2.0 * g_y_0_xx_yy[i] * a_exp + 4.0 * g_y_xx_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_xx_yz[i] = -2.0 * g_y_0_xx_yz[i] * a_exp + 4.0 * g_y_xx_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_xx_zz[i] = -2.0 * g_y_0_xx_zz[i] * a_exp + 4.0 * g_y_xx_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_y_0_xy_xx, g_y_0_xy_xy, g_y_0_xy_xz, g_y_0_xy_yy, g_y_0_xy_yz, g_y_0_xy_zz, g_y_x_0_0_0_x_xy_xx, g_y_x_0_0_0_x_xy_xy, g_y_x_0_0_0_x_xy_xz, g_y_x_0_0_0_x_xy_yy, g_y_x_0_0_0_x_xy_yz, g_y_x_0_0_0_x_xy_zz, g_y_xx_xy_xx, g_y_xx_xy_xy, g_y_xx_xy_xz, g_y_xx_xy_yy, g_y_xx_xy_yz, g_y_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_x_xy_xx[i] = -2.0 * g_y_0_xy_xx[i] * a_exp + 4.0 * g_y_xx_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_xy_xy[i] = -2.0 * g_y_0_xy_xy[i] * a_exp + 4.0 * g_y_xx_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_xy_xz[i] = -2.0 * g_y_0_xy_xz[i] * a_exp + 4.0 * g_y_xx_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_xy_yy[i] = -2.0 * g_y_0_xy_yy[i] * a_exp + 4.0 * g_y_xx_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_xy_yz[i] = -2.0 * g_y_0_xy_yz[i] * a_exp + 4.0 * g_y_xx_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_xy_zz[i] = -2.0 * g_y_0_xy_zz[i] * a_exp + 4.0 * g_y_xx_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_y_0_xz_xx, g_y_0_xz_xy, g_y_0_xz_xz, g_y_0_xz_yy, g_y_0_xz_yz, g_y_0_xz_zz, g_y_x_0_0_0_x_xz_xx, g_y_x_0_0_0_x_xz_xy, g_y_x_0_0_0_x_xz_xz, g_y_x_0_0_0_x_xz_yy, g_y_x_0_0_0_x_xz_yz, g_y_x_0_0_0_x_xz_zz, g_y_xx_xz_xx, g_y_xx_xz_xy, g_y_xx_xz_xz, g_y_xx_xz_yy, g_y_xx_xz_yz, g_y_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_x_xz_xx[i] = -2.0 * g_y_0_xz_xx[i] * a_exp + 4.0 * g_y_xx_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_xz_xy[i] = -2.0 * g_y_0_xz_xy[i] * a_exp + 4.0 * g_y_xx_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_xz_xz[i] = -2.0 * g_y_0_xz_xz[i] * a_exp + 4.0 * g_y_xx_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_xz_yy[i] = -2.0 * g_y_0_xz_yy[i] * a_exp + 4.0 * g_y_xx_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_xz_yz[i] = -2.0 * g_y_0_xz_yz[i] * a_exp + 4.0 * g_y_xx_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_xz_zz[i] = -2.0 * g_y_0_xz_zz[i] * a_exp + 4.0 * g_y_xx_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_y_0_yy_xx, g_y_0_yy_xy, g_y_0_yy_xz, g_y_0_yy_yy, g_y_0_yy_yz, g_y_0_yy_zz, g_y_x_0_0_0_x_yy_xx, g_y_x_0_0_0_x_yy_xy, g_y_x_0_0_0_x_yy_xz, g_y_x_0_0_0_x_yy_yy, g_y_x_0_0_0_x_yy_yz, g_y_x_0_0_0_x_yy_zz, g_y_xx_yy_xx, g_y_xx_yy_xy, g_y_xx_yy_xz, g_y_xx_yy_yy, g_y_xx_yy_yz, g_y_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_x_yy_xx[i] = -2.0 * g_y_0_yy_xx[i] * a_exp + 4.0 * g_y_xx_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_yy_xy[i] = -2.0 * g_y_0_yy_xy[i] * a_exp + 4.0 * g_y_xx_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_yy_xz[i] = -2.0 * g_y_0_yy_xz[i] * a_exp + 4.0 * g_y_xx_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_yy_yy[i] = -2.0 * g_y_0_yy_yy[i] * a_exp + 4.0 * g_y_xx_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_yy_yz[i] = -2.0 * g_y_0_yy_yz[i] * a_exp + 4.0 * g_y_xx_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_yy_zz[i] = -2.0 * g_y_0_yy_zz[i] * a_exp + 4.0 * g_y_xx_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_y_0_yz_xx, g_y_0_yz_xy, g_y_0_yz_xz, g_y_0_yz_yy, g_y_0_yz_yz, g_y_0_yz_zz, g_y_x_0_0_0_x_yz_xx, g_y_x_0_0_0_x_yz_xy, g_y_x_0_0_0_x_yz_xz, g_y_x_0_0_0_x_yz_yy, g_y_x_0_0_0_x_yz_yz, g_y_x_0_0_0_x_yz_zz, g_y_xx_yz_xx, g_y_xx_yz_xy, g_y_xx_yz_xz, g_y_xx_yz_yy, g_y_xx_yz_yz, g_y_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_x_yz_xx[i] = -2.0 * g_y_0_yz_xx[i] * a_exp + 4.0 * g_y_xx_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_yz_xy[i] = -2.0 * g_y_0_yz_xy[i] * a_exp + 4.0 * g_y_xx_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_yz_xz[i] = -2.0 * g_y_0_yz_xz[i] * a_exp + 4.0 * g_y_xx_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_yz_yy[i] = -2.0 * g_y_0_yz_yy[i] * a_exp + 4.0 * g_y_xx_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_yz_yz[i] = -2.0 * g_y_0_yz_yz[i] * a_exp + 4.0 * g_y_xx_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_yz_zz[i] = -2.0 * g_y_0_yz_zz[i] * a_exp + 4.0 * g_y_xx_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_y_0_zz_xx, g_y_0_zz_xy, g_y_0_zz_xz, g_y_0_zz_yy, g_y_0_zz_yz, g_y_0_zz_zz, g_y_x_0_0_0_x_zz_xx, g_y_x_0_0_0_x_zz_xy, g_y_x_0_0_0_x_zz_xz, g_y_x_0_0_0_x_zz_yy, g_y_x_0_0_0_x_zz_yz, g_y_x_0_0_0_x_zz_zz, g_y_xx_zz_xx, g_y_xx_zz_xy, g_y_xx_zz_xz, g_y_xx_zz_yy, g_y_xx_zz_yz, g_y_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_x_zz_xx[i] = -2.0 * g_y_0_zz_xx[i] * a_exp + 4.0 * g_y_xx_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_zz_xy[i] = -2.0 * g_y_0_zz_xy[i] * a_exp + 4.0 * g_y_xx_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_zz_xz[i] = -2.0 * g_y_0_zz_xz[i] * a_exp + 4.0 * g_y_xx_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_zz_yy[i] = -2.0 * g_y_0_zz_yy[i] * a_exp + 4.0 * g_y_xx_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_zz_yz[i] = -2.0 * g_y_0_zz_yz[i] * a_exp + 4.0 * g_y_xx_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_zz_zz[i] = -2.0 * g_y_0_zz_zz[i] * a_exp + 4.0 * g_y_xx_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_y_x_0_0_0_y_xx_xx, g_y_x_0_0_0_y_xx_xy, g_y_x_0_0_0_y_xx_xz, g_y_x_0_0_0_y_xx_yy, g_y_x_0_0_0_y_xx_yz, g_y_x_0_0_0_y_xx_zz, g_y_xy_xx_xx, g_y_xy_xx_xy, g_y_xy_xx_xz, g_y_xy_xx_yy, g_y_xy_xx_yz, g_y_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_y_xx_xx[i] = 4.0 * g_y_xy_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_xx_xy[i] = 4.0 * g_y_xy_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_xx_xz[i] = 4.0 * g_y_xy_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_xx_yy[i] = 4.0 * g_y_xy_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_xx_yz[i] = 4.0 * g_y_xy_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_xx_zz[i] = 4.0 * g_y_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_y_x_0_0_0_y_xy_xx, g_y_x_0_0_0_y_xy_xy, g_y_x_0_0_0_y_xy_xz, g_y_x_0_0_0_y_xy_yy, g_y_x_0_0_0_y_xy_yz, g_y_x_0_0_0_y_xy_zz, g_y_xy_xy_xx, g_y_xy_xy_xy, g_y_xy_xy_xz, g_y_xy_xy_yy, g_y_xy_xy_yz, g_y_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_y_xy_xx[i] = 4.0 * g_y_xy_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_xy_xy[i] = 4.0 * g_y_xy_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_xy_xz[i] = 4.0 * g_y_xy_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_xy_yy[i] = 4.0 * g_y_xy_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_xy_yz[i] = 4.0 * g_y_xy_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_xy_zz[i] = 4.0 * g_y_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_y_x_0_0_0_y_xz_xx, g_y_x_0_0_0_y_xz_xy, g_y_x_0_0_0_y_xz_xz, g_y_x_0_0_0_y_xz_yy, g_y_x_0_0_0_y_xz_yz, g_y_x_0_0_0_y_xz_zz, g_y_xy_xz_xx, g_y_xy_xz_xy, g_y_xy_xz_xz, g_y_xy_xz_yy, g_y_xy_xz_yz, g_y_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_y_xz_xx[i] = 4.0 * g_y_xy_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_xz_xy[i] = 4.0 * g_y_xy_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_xz_xz[i] = 4.0 * g_y_xy_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_xz_yy[i] = 4.0 * g_y_xy_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_xz_yz[i] = 4.0 * g_y_xy_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_xz_zz[i] = 4.0 * g_y_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_y_x_0_0_0_y_yy_xx, g_y_x_0_0_0_y_yy_xy, g_y_x_0_0_0_y_yy_xz, g_y_x_0_0_0_y_yy_yy, g_y_x_0_0_0_y_yy_yz, g_y_x_0_0_0_y_yy_zz, g_y_xy_yy_xx, g_y_xy_yy_xy, g_y_xy_yy_xz, g_y_xy_yy_yy, g_y_xy_yy_yz, g_y_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_y_yy_xx[i] = 4.0 * g_y_xy_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_yy_xy[i] = 4.0 * g_y_xy_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_yy_xz[i] = 4.0 * g_y_xy_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_yy_yy[i] = 4.0 * g_y_xy_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_yy_yz[i] = 4.0 * g_y_xy_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_yy_zz[i] = 4.0 * g_y_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_y_x_0_0_0_y_yz_xx, g_y_x_0_0_0_y_yz_xy, g_y_x_0_0_0_y_yz_xz, g_y_x_0_0_0_y_yz_yy, g_y_x_0_0_0_y_yz_yz, g_y_x_0_0_0_y_yz_zz, g_y_xy_yz_xx, g_y_xy_yz_xy, g_y_xy_yz_xz, g_y_xy_yz_yy, g_y_xy_yz_yz, g_y_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_y_yz_xx[i] = 4.0 * g_y_xy_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_yz_xy[i] = 4.0 * g_y_xy_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_yz_xz[i] = 4.0 * g_y_xy_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_yz_yy[i] = 4.0 * g_y_xy_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_yz_yz[i] = 4.0 * g_y_xy_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_yz_zz[i] = 4.0 * g_y_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_y_x_0_0_0_y_zz_xx, g_y_x_0_0_0_y_zz_xy, g_y_x_0_0_0_y_zz_xz, g_y_x_0_0_0_y_zz_yy, g_y_x_0_0_0_y_zz_yz, g_y_x_0_0_0_y_zz_zz, g_y_xy_zz_xx, g_y_xy_zz_xy, g_y_xy_zz_xz, g_y_xy_zz_yy, g_y_xy_zz_yz, g_y_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_y_zz_xx[i] = 4.0 * g_y_xy_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_zz_xy[i] = 4.0 * g_y_xy_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_zz_xz[i] = 4.0 * g_y_xy_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_zz_yy[i] = 4.0 * g_y_xy_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_zz_yz[i] = 4.0 * g_y_xy_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_zz_zz[i] = 4.0 * g_y_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_y_x_0_0_0_z_xx_xx, g_y_x_0_0_0_z_xx_xy, g_y_x_0_0_0_z_xx_xz, g_y_x_0_0_0_z_xx_yy, g_y_x_0_0_0_z_xx_yz, g_y_x_0_0_0_z_xx_zz, g_y_xz_xx_xx, g_y_xz_xx_xy, g_y_xz_xx_xz, g_y_xz_xx_yy, g_y_xz_xx_yz, g_y_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_z_xx_xx[i] = 4.0 * g_y_xz_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_xx_xy[i] = 4.0 * g_y_xz_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_xx_xz[i] = 4.0 * g_y_xz_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_xx_yy[i] = 4.0 * g_y_xz_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_xx_yz[i] = 4.0 * g_y_xz_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_xx_zz[i] = 4.0 * g_y_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_y_x_0_0_0_z_xy_xx, g_y_x_0_0_0_z_xy_xy, g_y_x_0_0_0_z_xy_xz, g_y_x_0_0_0_z_xy_yy, g_y_x_0_0_0_z_xy_yz, g_y_x_0_0_0_z_xy_zz, g_y_xz_xy_xx, g_y_xz_xy_xy, g_y_xz_xy_xz, g_y_xz_xy_yy, g_y_xz_xy_yz, g_y_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_z_xy_xx[i] = 4.0 * g_y_xz_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_xy_xy[i] = 4.0 * g_y_xz_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_xy_xz[i] = 4.0 * g_y_xz_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_xy_yy[i] = 4.0 * g_y_xz_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_xy_yz[i] = 4.0 * g_y_xz_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_xy_zz[i] = 4.0 * g_y_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_y_x_0_0_0_z_xz_xx, g_y_x_0_0_0_z_xz_xy, g_y_x_0_0_0_z_xz_xz, g_y_x_0_0_0_z_xz_yy, g_y_x_0_0_0_z_xz_yz, g_y_x_0_0_0_z_xz_zz, g_y_xz_xz_xx, g_y_xz_xz_xy, g_y_xz_xz_xz, g_y_xz_xz_yy, g_y_xz_xz_yz, g_y_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_z_xz_xx[i] = 4.0 * g_y_xz_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_xz_xy[i] = 4.0 * g_y_xz_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_xz_xz[i] = 4.0 * g_y_xz_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_xz_yy[i] = 4.0 * g_y_xz_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_xz_yz[i] = 4.0 * g_y_xz_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_xz_zz[i] = 4.0 * g_y_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_y_x_0_0_0_z_yy_xx, g_y_x_0_0_0_z_yy_xy, g_y_x_0_0_0_z_yy_xz, g_y_x_0_0_0_z_yy_yy, g_y_x_0_0_0_z_yy_yz, g_y_x_0_0_0_z_yy_zz, g_y_xz_yy_xx, g_y_xz_yy_xy, g_y_xz_yy_xz, g_y_xz_yy_yy, g_y_xz_yy_yz, g_y_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_z_yy_xx[i] = 4.0 * g_y_xz_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_yy_xy[i] = 4.0 * g_y_xz_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_yy_xz[i] = 4.0 * g_y_xz_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_yy_yy[i] = 4.0 * g_y_xz_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_yy_yz[i] = 4.0 * g_y_xz_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_yy_zz[i] = 4.0 * g_y_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_y_x_0_0_0_z_yz_xx, g_y_x_0_0_0_z_yz_xy, g_y_x_0_0_0_z_yz_xz, g_y_x_0_0_0_z_yz_yy, g_y_x_0_0_0_z_yz_yz, g_y_x_0_0_0_z_yz_zz, g_y_xz_yz_xx, g_y_xz_yz_xy, g_y_xz_yz_xz, g_y_xz_yz_yy, g_y_xz_yz_yz, g_y_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_z_yz_xx[i] = 4.0 * g_y_xz_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_yz_xy[i] = 4.0 * g_y_xz_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_yz_xz[i] = 4.0 * g_y_xz_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_yz_yy[i] = 4.0 * g_y_xz_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_yz_yz[i] = 4.0 * g_y_xz_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_yz_zz[i] = 4.0 * g_y_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_y_x_0_0_0_z_zz_xx, g_y_x_0_0_0_z_zz_xy, g_y_x_0_0_0_z_zz_xz, g_y_x_0_0_0_z_zz_yy, g_y_x_0_0_0_z_zz_yz, g_y_x_0_0_0_z_zz_zz, g_y_xz_zz_xx, g_y_xz_zz_xy, g_y_xz_zz_xz, g_y_xz_zz_yy, g_y_xz_zz_yz, g_y_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_z_zz_xx[i] = 4.0 * g_y_xz_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_zz_xy[i] = 4.0 * g_y_xz_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_zz_xz[i] = 4.0 * g_y_xz_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_zz_yy[i] = 4.0 * g_y_xz_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_zz_yz[i] = 4.0 * g_y_xz_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_zz_zz[i] = 4.0 * g_y_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_y_xy_xx_xx, g_y_xy_xx_xy, g_y_xy_xx_xz, g_y_xy_xx_yy, g_y_xy_xx_yz, g_y_xy_xx_zz, g_y_y_0_0_0_x_xx_xx, g_y_y_0_0_0_x_xx_xy, g_y_y_0_0_0_x_xx_xz, g_y_y_0_0_0_x_xx_yy, g_y_y_0_0_0_x_xx_yz, g_y_y_0_0_0_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_x_xx_xx[i] = 4.0 * g_y_xy_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_xx_xy[i] = 4.0 * g_y_xy_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_xx_xz[i] = 4.0 * g_y_xy_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_xx_yy[i] = 4.0 * g_y_xy_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_xx_yz[i] = 4.0 * g_y_xy_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_xx_zz[i] = 4.0 * g_y_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_y_xy_xy_xx, g_y_xy_xy_xy, g_y_xy_xy_xz, g_y_xy_xy_yy, g_y_xy_xy_yz, g_y_xy_xy_zz, g_y_y_0_0_0_x_xy_xx, g_y_y_0_0_0_x_xy_xy, g_y_y_0_0_0_x_xy_xz, g_y_y_0_0_0_x_xy_yy, g_y_y_0_0_0_x_xy_yz, g_y_y_0_0_0_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_x_xy_xx[i] = 4.0 * g_y_xy_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_xy_xy[i] = 4.0 * g_y_xy_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_xy_xz[i] = 4.0 * g_y_xy_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_xy_yy[i] = 4.0 * g_y_xy_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_xy_yz[i] = 4.0 * g_y_xy_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_xy_zz[i] = 4.0 * g_y_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_y_xy_xz_xx, g_y_xy_xz_xy, g_y_xy_xz_xz, g_y_xy_xz_yy, g_y_xy_xz_yz, g_y_xy_xz_zz, g_y_y_0_0_0_x_xz_xx, g_y_y_0_0_0_x_xz_xy, g_y_y_0_0_0_x_xz_xz, g_y_y_0_0_0_x_xz_yy, g_y_y_0_0_0_x_xz_yz, g_y_y_0_0_0_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_x_xz_xx[i] = 4.0 * g_y_xy_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_xz_xy[i] = 4.0 * g_y_xy_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_xz_xz[i] = 4.0 * g_y_xy_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_xz_yy[i] = 4.0 * g_y_xy_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_xz_yz[i] = 4.0 * g_y_xy_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_xz_zz[i] = 4.0 * g_y_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_y_xy_yy_xx, g_y_xy_yy_xy, g_y_xy_yy_xz, g_y_xy_yy_yy, g_y_xy_yy_yz, g_y_xy_yy_zz, g_y_y_0_0_0_x_yy_xx, g_y_y_0_0_0_x_yy_xy, g_y_y_0_0_0_x_yy_xz, g_y_y_0_0_0_x_yy_yy, g_y_y_0_0_0_x_yy_yz, g_y_y_0_0_0_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_x_yy_xx[i] = 4.0 * g_y_xy_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_yy_xy[i] = 4.0 * g_y_xy_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_yy_xz[i] = 4.0 * g_y_xy_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_yy_yy[i] = 4.0 * g_y_xy_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_yy_yz[i] = 4.0 * g_y_xy_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_yy_zz[i] = 4.0 * g_y_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_y_xy_yz_xx, g_y_xy_yz_xy, g_y_xy_yz_xz, g_y_xy_yz_yy, g_y_xy_yz_yz, g_y_xy_yz_zz, g_y_y_0_0_0_x_yz_xx, g_y_y_0_0_0_x_yz_xy, g_y_y_0_0_0_x_yz_xz, g_y_y_0_0_0_x_yz_yy, g_y_y_0_0_0_x_yz_yz, g_y_y_0_0_0_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_x_yz_xx[i] = 4.0 * g_y_xy_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_yz_xy[i] = 4.0 * g_y_xy_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_yz_xz[i] = 4.0 * g_y_xy_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_yz_yy[i] = 4.0 * g_y_xy_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_yz_yz[i] = 4.0 * g_y_xy_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_yz_zz[i] = 4.0 * g_y_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_y_xy_zz_xx, g_y_xy_zz_xy, g_y_xy_zz_xz, g_y_xy_zz_yy, g_y_xy_zz_yz, g_y_xy_zz_zz, g_y_y_0_0_0_x_zz_xx, g_y_y_0_0_0_x_zz_xy, g_y_y_0_0_0_x_zz_xz, g_y_y_0_0_0_x_zz_yy, g_y_y_0_0_0_x_zz_yz, g_y_y_0_0_0_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_x_zz_xx[i] = 4.0 * g_y_xy_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_zz_xy[i] = 4.0 * g_y_xy_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_zz_xz[i] = 4.0 * g_y_xy_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_zz_yy[i] = 4.0 * g_y_xy_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_zz_yz[i] = 4.0 * g_y_xy_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_zz_zz[i] = 4.0 * g_y_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_y_0_xx_xx, g_y_0_xx_xy, g_y_0_xx_xz, g_y_0_xx_yy, g_y_0_xx_yz, g_y_0_xx_zz, g_y_y_0_0_0_y_xx_xx, g_y_y_0_0_0_y_xx_xy, g_y_y_0_0_0_y_xx_xz, g_y_y_0_0_0_y_xx_yy, g_y_y_0_0_0_y_xx_yz, g_y_y_0_0_0_y_xx_zz, g_y_yy_xx_xx, g_y_yy_xx_xy, g_y_yy_xx_xz, g_y_yy_xx_yy, g_y_yy_xx_yz, g_y_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_y_xx_xx[i] = -2.0 * g_y_0_xx_xx[i] * a_exp + 4.0 * g_y_yy_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_xx_xy[i] = -2.0 * g_y_0_xx_xy[i] * a_exp + 4.0 * g_y_yy_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_xx_xz[i] = -2.0 * g_y_0_xx_xz[i] * a_exp + 4.0 * g_y_yy_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_xx_yy[i] = -2.0 * g_y_0_xx_yy[i] * a_exp + 4.0 * g_y_yy_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_xx_yz[i] = -2.0 * g_y_0_xx_yz[i] * a_exp + 4.0 * g_y_yy_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_xx_zz[i] = -2.0 * g_y_0_xx_zz[i] * a_exp + 4.0 * g_y_yy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_y_0_xy_xx, g_y_0_xy_xy, g_y_0_xy_xz, g_y_0_xy_yy, g_y_0_xy_yz, g_y_0_xy_zz, g_y_y_0_0_0_y_xy_xx, g_y_y_0_0_0_y_xy_xy, g_y_y_0_0_0_y_xy_xz, g_y_y_0_0_0_y_xy_yy, g_y_y_0_0_0_y_xy_yz, g_y_y_0_0_0_y_xy_zz, g_y_yy_xy_xx, g_y_yy_xy_xy, g_y_yy_xy_xz, g_y_yy_xy_yy, g_y_yy_xy_yz, g_y_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_y_xy_xx[i] = -2.0 * g_y_0_xy_xx[i] * a_exp + 4.0 * g_y_yy_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_xy_xy[i] = -2.0 * g_y_0_xy_xy[i] * a_exp + 4.0 * g_y_yy_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_xy_xz[i] = -2.0 * g_y_0_xy_xz[i] * a_exp + 4.0 * g_y_yy_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_xy_yy[i] = -2.0 * g_y_0_xy_yy[i] * a_exp + 4.0 * g_y_yy_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_xy_yz[i] = -2.0 * g_y_0_xy_yz[i] * a_exp + 4.0 * g_y_yy_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_xy_zz[i] = -2.0 * g_y_0_xy_zz[i] * a_exp + 4.0 * g_y_yy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_y_0_xz_xx, g_y_0_xz_xy, g_y_0_xz_xz, g_y_0_xz_yy, g_y_0_xz_yz, g_y_0_xz_zz, g_y_y_0_0_0_y_xz_xx, g_y_y_0_0_0_y_xz_xy, g_y_y_0_0_0_y_xz_xz, g_y_y_0_0_0_y_xz_yy, g_y_y_0_0_0_y_xz_yz, g_y_y_0_0_0_y_xz_zz, g_y_yy_xz_xx, g_y_yy_xz_xy, g_y_yy_xz_xz, g_y_yy_xz_yy, g_y_yy_xz_yz, g_y_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_y_xz_xx[i] = -2.0 * g_y_0_xz_xx[i] * a_exp + 4.0 * g_y_yy_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_xz_xy[i] = -2.0 * g_y_0_xz_xy[i] * a_exp + 4.0 * g_y_yy_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_xz_xz[i] = -2.0 * g_y_0_xz_xz[i] * a_exp + 4.0 * g_y_yy_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_xz_yy[i] = -2.0 * g_y_0_xz_yy[i] * a_exp + 4.0 * g_y_yy_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_xz_yz[i] = -2.0 * g_y_0_xz_yz[i] * a_exp + 4.0 * g_y_yy_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_xz_zz[i] = -2.0 * g_y_0_xz_zz[i] * a_exp + 4.0 * g_y_yy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_y_0_yy_xx, g_y_0_yy_xy, g_y_0_yy_xz, g_y_0_yy_yy, g_y_0_yy_yz, g_y_0_yy_zz, g_y_y_0_0_0_y_yy_xx, g_y_y_0_0_0_y_yy_xy, g_y_y_0_0_0_y_yy_xz, g_y_y_0_0_0_y_yy_yy, g_y_y_0_0_0_y_yy_yz, g_y_y_0_0_0_y_yy_zz, g_y_yy_yy_xx, g_y_yy_yy_xy, g_y_yy_yy_xz, g_y_yy_yy_yy, g_y_yy_yy_yz, g_y_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_y_yy_xx[i] = -2.0 * g_y_0_yy_xx[i] * a_exp + 4.0 * g_y_yy_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_yy_xy[i] = -2.0 * g_y_0_yy_xy[i] * a_exp + 4.0 * g_y_yy_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_yy_xz[i] = -2.0 * g_y_0_yy_xz[i] * a_exp + 4.0 * g_y_yy_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_yy_yy[i] = -2.0 * g_y_0_yy_yy[i] * a_exp + 4.0 * g_y_yy_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_yy_yz[i] = -2.0 * g_y_0_yy_yz[i] * a_exp + 4.0 * g_y_yy_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_yy_zz[i] = -2.0 * g_y_0_yy_zz[i] * a_exp + 4.0 * g_y_yy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_y_0_yz_xx, g_y_0_yz_xy, g_y_0_yz_xz, g_y_0_yz_yy, g_y_0_yz_yz, g_y_0_yz_zz, g_y_y_0_0_0_y_yz_xx, g_y_y_0_0_0_y_yz_xy, g_y_y_0_0_0_y_yz_xz, g_y_y_0_0_0_y_yz_yy, g_y_y_0_0_0_y_yz_yz, g_y_y_0_0_0_y_yz_zz, g_y_yy_yz_xx, g_y_yy_yz_xy, g_y_yy_yz_xz, g_y_yy_yz_yy, g_y_yy_yz_yz, g_y_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_y_yz_xx[i] = -2.0 * g_y_0_yz_xx[i] * a_exp + 4.0 * g_y_yy_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_yz_xy[i] = -2.0 * g_y_0_yz_xy[i] * a_exp + 4.0 * g_y_yy_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_yz_xz[i] = -2.0 * g_y_0_yz_xz[i] * a_exp + 4.0 * g_y_yy_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_yz_yy[i] = -2.0 * g_y_0_yz_yy[i] * a_exp + 4.0 * g_y_yy_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_yz_yz[i] = -2.0 * g_y_0_yz_yz[i] * a_exp + 4.0 * g_y_yy_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_yz_zz[i] = -2.0 * g_y_0_yz_zz[i] * a_exp + 4.0 * g_y_yy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_y_0_zz_xx, g_y_0_zz_xy, g_y_0_zz_xz, g_y_0_zz_yy, g_y_0_zz_yz, g_y_0_zz_zz, g_y_y_0_0_0_y_zz_xx, g_y_y_0_0_0_y_zz_xy, g_y_y_0_0_0_y_zz_xz, g_y_y_0_0_0_y_zz_yy, g_y_y_0_0_0_y_zz_yz, g_y_y_0_0_0_y_zz_zz, g_y_yy_zz_xx, g_y_yy_zz_xy, g_y_yy_zz_xz, g_y_yy_zz_yy, g_y_yy_zz_yz, g_y_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_y_zz_xx[i] = -2.0 * g_y_0_zz_xx[i] * a_exp + 4.0 * g_y_yy_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_zz_xy[i] = -2.0 * g_y_0_zz_xy[i] * a_exp + 4.0 * g_y_yy_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_zz_xz[i] = -2.0 * g_y_0_zz_xz[i] * a_exp + 4.0 * g_y_yy_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_zz_yy[i] = -2.0 * g_y_0_zz_yy[i] * a_exp + 4.0 * g_y_yy_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_zz_yz[i] = -2.0 * g_y_0_zz_yz[i] * a_exp + 4.0 * g_y_yy_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_zz_zz[i] = -2.0 * g_y_0_zz_zz[i] * a_exp + 4.0 * g_y_yy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_y_y_0_0_0_z_xx_xx, g_y_y_0_0_0_z_xx_xy, g_y_y_0_0_0_z_xx_xz, g_y_y_0_0_0_z_xx_yy, g_y_y_0_0_0_z_xx_yz, g_y_y_0_0_0_z_xx_zz, g_y_yz_xx_xx, g_y_yz_xx_xy, g_y_yz_xx_xz, g_y_yz_xx_yy, g_y_yz_xx_yz, g_y_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_z_xx_xx[i] = 4.0 * g_y_yz_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_xx_xy[i] = 4.0 * g_y_yz_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_xx_xz[i] = 4.0 * g_y_yz_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_xx_yy[i] = 4.0 * g_y_yz_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_xx_yz[i] = 4.0 * g_y_yz_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_xx_zz[i] = 4.0 * g_y_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_y_y_0_0_0_z_xy_xx, g_y_y_0_0_0_z_xy_xy, g_y_y_0_0_0_z_xy_xz, g_y_y_0_0_0_z_xy_yy, g_y_y_0_0_0_z_xy_yz, g_y_y_0_0_0_z_xy_zz, g_y_yz_xy_xx, g_y_yz_xy_xy, g_y_yz_xy_xz, g_y_yz_xy_yy, g_y_yz_xy_yz, g_y_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_z_xy_xx[i] = 4.0 * g_y_yz_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_xy_xy[i] = 4.0 * g_y_yz_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_xy_xz[i] = 4.0 * g_y_yz_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_xy_yy[i] = 4.0 * g_y_yz_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_xy_yz[i] = 4.0 * g_y_yz_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_xy_zz[i] = 4.0 * g_y_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_y_y_0_0_0_z_xz_xx, g_y_y_0_0_0_z_xz_xy, g_y_y_0_0_0_z_xz_xz, g_y_y_0_0_0_z_xz_yy, g_y_y_0_0_0_z_xz_yz, g_y_y_0_0_0_z_xz_zz, g_y_yz_xz_xx, g_y_yz_xz_xy, g_y_yz_xz_xz, g_y_yz_xz_yy, g_y_yz_xz_yz, g_y_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_z_xz_xx[i] = 4.0 * g_y_yz_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_xz_xy[i] = 4.0 * g_y_yz_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_xz_xz[i] = 4.0 * g_y_yz_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_xz_yy[i] = 4.0 * g_y_yz_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_xz_yz[i] = 4.0 * g_y_yz_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_xz_zz[i] = 4.0 * g_y_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_y_y_0_0_0_z_yy_xx, g_y_y_0_0_0_z_yy_xy, g_y_y_0_0_0_z_yy_xz, g_y_y_0_0_0_z_yy_yy, g_y_y_0_0_0_z_yy_yz, g_y_y_0_0_0_z_yy_zz, g_y_yz_yy_xx, g_y_yz_yy_xy, g_y_yz_yy_xz, g_y_yz_yy_yy, g_y_yz_yy_yz, g_y_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_z_yy_xx[i] = 4.0 * g_y_yz_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_yy_xy[i] = 4.0 * g_y_yz_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_yy_xz[i] = 4.0 * g_y_yz_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_yy_yy[i] = 4.0 * g_y_yz_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_yy_yz[i] = 4.0 * g_y_yz_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_yy_zz[i] = 4.0 * g_y_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_y_y_0_0_0_z_yz_xx, g_y_y_0_0_0_z_yz_xy, g_y_y_0_0_0_z_yz_xz, g_y_y_0_0_0_z_yz_yy, g_y_y_0_0_0_z_yz_yz, g_y_y_0_0_0_z_yz_zz, g_y_yz_yz_xx, g_y_yz_yz_xy, g_y_yz_yz_xz, g_y_yz_yz_yy, g_y_yz_yz_yz, g_y_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_z_yz_xx[i] = 4.0 * g_y_yz_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_yz_xy[i] = 4.0 * g_y_yz_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_yz_xz[i] = 4.0 * g_y_yz_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_yz_yy[i] = 4.0 * g_y_yz_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_yz_yz[i] = 4.0 * g_y_yz_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_yz_zz[i] = 4.0 * g_y_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_y_y_0_0_0_z_zz_xx, g_y_y_0_0_0_z_zz_xy, g_y_y_0_0_0_z_zz_xz, g_y_y_0_0_0_z_zz_yy, g_y_y_0_0_0_z_zz_yz, g_y_y_0_0_0_z_zz_zz, g_y_yz_zz_xx, g_y_yz_zz_xy, g_y_yz_zz_xz, g_y_yz_zz_yy, g_y_yz_zz_yz, g_y_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_z_zz_xx[i] = 4.0 * g_y_yz_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_zz_xy[i] = 4.0 * g_y_yz_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_zz_xz[i] = 4.0 * g_y_yz_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_zz_yy[i] = 4.0 * g_y_yz_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_zz_yz[i] = 4.0 * g_y_yz_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_zz_zz[i] = 4.0 * g_y_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_y_xz_xx_xx, g_y_xz_xx_xy, g_y_xz_xx_xz, g_y_xz_xx_yy, g_y_xz_xx_yz, g_y_xz_xx_zz, g_y_z_0_0_0_x_xx_xx, g_y_z_0_0_0_x_xx_xy, g_y_z_0_0_0_x_xx_xz, g_y_z_0_0_0_x_xx_yy, g_y_z_0_0_0_x_xx_yz, g_y_z_0_0_0_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_x_xx_xx[i] = 4.0 * g_y_xz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_xx_xy[i] = 4.0 * g_y_xz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_xx_xz[i] = 4.0 * g_y_xz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_xx_yy[i] = 4.0 * g_y_xz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_xx_yz[i] = 4.0 * g_y_xz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_xx_zz[i] = 4.0 * g_y_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_y_xz_xy_xx, g_y_xz_xy_xy, g_y_xz_xy_xz, g_y_xz_xy_yy, g_y_xz_xy_yz, g_y_xz_xy_zz, g_y_z_0_0_0_x_xy_xx, g_y_z_0_0_0_x_xy_xy, g_y_z_0_0_0_x_xy_xz, g_y_z_0_0_0_x_xy_yy, g_y_z_0_0_0_x_xy_yz, g_y_z_0_0_0_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_x_xy_xx[i] = 4.0 * g_y_xz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_xy_xy[i] = 4.0 * g_y_xz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_xy_xz[i] = 4.0 * g_y_xz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_xy_yy[i] = 4.0 * g_y_xz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_xy_yz[i] = 4.0 * g_y_xz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_xy_zz[i] = 4.0 * g_y_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_y_xz_xz_xx, g_y_xz_xz_xy, g_y_xz_xz_xz, g_y_xz_xz_yy, g_y_xz_xz_yz, g_y_xz_xz_zz, g_y_z_0_0_0_x_xz_xx, g_y_z_0_0_0_x_xz_xy, g_y_z_0_0_0_x_xz_xz, g_y_z_0_0_0_x_xz_yy, g_y_z_0_0_0_x_xz_yz, g_y_z_0_0_0_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_x_xz_xx[i] = 4.0 * g_y_xz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_xz_xy[i] = 4.0 * g_y_xz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_xz_xz[i] = 4.0 * g_y_xz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_xz_yy[i] = 4.0 * g_y_xz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_xz_yz[i] = 4.0 * g_y_xz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_xz_zz[i] = 4.0 * g_y_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_y_xz_yy_xx, g_y_xz_yy_xy, g_y_xz_yy_xz, g_y_xz_yy_yy, g_y_xz_yy_yz, g_y_xz_yy_zz, g_y_z_0_0_0_x_yy_xx, g_y_z_0_0_0_x_yy_xy, g_y_z_0_0_0_x_yy_xz, g_y_z_0_0_0_x_yy_yy, g_y_z_0_0_0_x_yy_yz, g_y_z_0_0_0_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_x_yy_xx[i] = 4.0 * g_y_xz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_yy_xy[i] = 4.0 * g_y_xz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_yy_xz[i] = 4.0 * g_y_xz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_yy_yy[i] = 4.0 * g_y_xz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_yy_yz[i] = 4.0 * g_y_xz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_yy_zz[i] = 4.0 * g_y_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_y_xz_yz_xx, g_y_xz_yz_xy, g_y_xz_yz_xz, g_y_xz_yz_yy, g_y_xz_yz_yz, g_y_xz_yz_zz, g_y_z_0_0_0_x_yz_xx, g_y_z_0_0_0_x_yz_xy, g_y_z_0_0_0_x_yz_xz, g_y_z_0_0_0_x_yz_yy, g_y_z_0_0_0_x_yz_yz, g_y_z_0_0_0_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_x_yz_xx[i] = 4.0 * g_y_xz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_yz_xy[i] = 4.0 * g_y_xz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_yz_xz[i] = 4.0 * g_y_xz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_yz_yy[i] = 4.0 * g_y_xz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_yz_yz[i] = 4.0 * g_y_xz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_yz_zz[i] = 4.0 * g_y_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_y_xz_zz_xx, g_y_xz_zz_xy, g_y_xz_zz_xz, g_y_xz_zz_yy, g_y_xz_zz_yz, g_y_xz_zz_zz, g_y_z_0_0_0_x_zz_xx, g_y_z_0_0_0_x_zz_xy, g_y_z_0_0_0_x_zz_xz, g_y_z_0_0_0_x_zz_yy, g_y_z_0_0_0_x_zz_yz, g_y_z_0_0_0_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_x_zz_xx[i] = 4.0 * g_y_xz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_zz_xy[i] = 4.0 * g_y_xz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_zz_xz[i] = 4.0 * g_y_xz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_zz_yy[i] = 4.0 * g_y_xz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_zz_yz[i] = 4.0 * g_y_xz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_zz_zz[i] = 4.0 * g_y_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_y_yz_xx_xx, g_y_yz_xx_xy, g_y_yz_xx_xz, g_y_yz_xx_yy, g_y_yz_xx_yz, g_y_yz_xx_zz, g_y_z_0_0_0_y_xx_xx, g_y_z_0_0_0_y_xx_xy, g_y_z_0_0_0_y_xx_xz, g_y_z_0_0_0_y_xx_yy, g_y_z_0_0_0_y_xx_yz, g_y_z_0_0_0_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_y_xx_xx[i] = 4.0 * g_y_yz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_xx_xy[i] = 4.0 * g_y_yz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_xx_xz[i] = 4.0 * g_y_yz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_xx_yy[i] = 4.0 * g_y_yz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_xx_yz[i] = 4.0 * g_y_yz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_xx_zz[i] = 4.0 * g_y_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_y_yz_xy_xx, g_y_yz_xy_xy, g_y_yz_xy_xz, g_y_yz_xy_yy, g_y_yz_xy_yz, g_y_yz_xy_zz, g_y_z_0_0_0_y_xy_xx, g_y_z_0_0_0_y_xy_xy, g_y_z_0_0_0_y_xy_xz, g_y_z_0_0_0_y_xy_yy, g_y_z_0_0_0_y_xy_yz, g_y_z_0_0_0_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_y_xy_xx[i] = 4.0 * g_y_yz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_xy_xy[i] = 4.0 * g_y_yz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_xy_xz[i] = 4.0 * g_y_yz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_xy_yy[i] = 4.0 * g_y_yz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_xy_yz[i] = 4.0 * g_y_yz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_xy_zz[i] = 4.0 * g_y_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_y_yz_xz_xx, g_y_yz_xz_xy, g_y_yz_xz_xz, g_y_yz_xz_yy, g_y_yz_xz_yz, g_y_yz_xz_zz, g_y_z_0_0_0_y_xz_xx, g_y_z_0_0_0_y_xz_xy, g_y_z_0_0_0_y_xz_xz, g_y_z_0_0_0_y_xz_yy, g_y_z_0_0_0_y_xz_yz, g_y_z_0_0_0_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_y_xz_xx[i] = 4.0 * g_y_yz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_xz_xy[i] = 4.0 * g_y_yz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_xz_xz[i] = 4.0 * g_y_yz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_xz_yy[i] = 4.0 * g_y_yz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_xz_yz[i] = 4.0 * g_y_yz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_xz_zz[i] = 4.0 * g_y_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_y_yz_yy_xx, g_y_yz_yy_xy, g_y_yz_yy_xz, g_y_yz_yy_yy, g_y_yz_yy_yz, g_y_yz_yy_zz, g_y_z_0_0_0_y_yy_xx, g_y_z_0_0_0_y_yy_xy, g_y_z_0_0_0_y_yy_xz, g_y_z_0_0_0_y_yy_yy, g_y_z_0_0_0_y_yy_yz, g_y_z_0_0_0_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_y_yy_xx[i] = 4.0 * g_y_yz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_yy_xy[i] = 4.0 * g_y_yz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_yy_xz[i] = 4.0 * g_y_yz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_yy_yy[i] = 4.0 * g_y_yz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_yy_yz[i] = 4.0 * g_y_yz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_yy_zz[i] = 4.0 * g_y_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_y_yz_yz_xx, g_y_yz_yz_xy, g_y_yz_yz_xz, g_y_yz_yz_yy, g_y_yz_yz_yz, g_y_yz_yz_zz, g_y_z_0_0_0_y_yz_xx, g_y_z_0_0_0_y_yz_xy, g_y_z_0_0_0_y_yz_xz, g_y_z_0_0_0_y_yz_yy, g_y_z_0_0_0_y_yz_yz, g_y_z_0_0_0_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_y_yz_xx[i] = 4.0 * g_y_yz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_yz_xy[i] = 4.0 * g_y_yz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_yz_xz[i] = 4.0 * g_y_yz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_yz_yy[i] = 4.0 * g_y_yz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_yz_yz[i] = 4.0 * g_y_yz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_yz_zz[i] = 4.0 * g_y_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_y_yz_zz_xx, g_y_yz_zz_xy, g_y_yz_zz_xz, g_y_yz_zz_yy, g_y_yz_zz_yz, g_y_yz_zz_zz, g_y_z_0_0_0_y_zz_xx, g_y_z_0_0_0_y_zz_xy, g_y_z_0_0_0_y_zz_xz, g_y_z_0_0_0_y_zz_yy, g_y_z_0_0_0_y_zz_yz, g_y_z_0_0_0_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_y_zz_xx[i] = 4.0 * g_y_yz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_zz_xy[i] = 4.0 * g_y_yz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_zz_xz[i] = 4.0 * g_y_yz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_zz_yy[i] = 4.0 * g_y_yz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_zz_yz[i] = 4.0 * g_y_yz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_zz_zz[i] = 4.0 * g_y_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_y_0_xx_xx, g_y_0_xx_xy, g_y_0_xx_xz, g_y_0_xx_yy, g_y_0_xx_yz, g_y_0_xx_zz, g_y_z_0_0_0_z_xx_xx, g_y_z_0_0_0_z_xx_xy, g_y_z_0_0_0_z_xx_xz, g_y_z_0_0_0_z_xx_yy, g_y_z_0_0_0_z_xx_yz, g_y_z_0_0_0_z_xx_zz, g_y_zz_xx_xx, g_y_zz_xx_xy, g_y_zz_xx_xz, g_y_zz_xx_yy, g_y_zz_xx_yz, g_y_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_z_xx_xx[i] = -2.0 * g_y_0_xx_xx[i] * a_exp + 4.0 * g_y_zz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_xx_xy[i] = -2.0 * g_y_0_xx_xy[i] * a_exp + 4.0 * g_y_zz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_xx_xz[i] = -2.0 * g_y_0_xx_xz[i] * a_exp + 4.0 * g_y_zz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_xx_yy[i] = -2.0 * g_y_0_xx_yy[i] * a_exp + 4.0 * g_y_zz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_xx_yz[i] = -2.0 * g_y_0_xx_yz[i] * a_exp + 4.0 * g_y_zz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_xx_zz[i] = -2.0 * g_y_0_xx_zz[i] * a_exp + 4.0 * g_y_zz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_y_0_xy_xx, g_y_0_xy_xy, g_y_0_xy_xz, g_y_0_xy_yy, g_y_0_xy_yz, g_y_0_xy_zz, g_y_z_0_0_0_z_xy_xx, g_y_z_0_0_0_z_xy_xy, g_y_z_0_0_0_z_xy_xz, g_y_z_0_0_0_z_xy_yy, g_y_z_0_0_0_z_xy_yz, g_y_z_0_0_0_z_xy_zz, g_y_zz_xy_xx, g_y_zz_xy_xy, g_y_zz_xy_xz, g_y_zz_xy_yy, g_y_zz_xy_yz, g_y_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_z_xy_xx[i] = -2.0 * g_y_0_xy_xx[i] * a_exp + 4.0 * g_y_zz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_xy_xy[i] = -2.0 * g_y_0_xy_xy[i] * a_exp + 4.0 * g_y_zz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_xy_xz[i] = -2.0 * g_y_0_xy_xz[i] * a_exp + 4.0 * g_y_zz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_xy_yy[i] = -2.0 * g_y_0_xy_yy[i] * a_exp + 4.0 * g_y_zz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_xy_yz[i] = -2.0 * g_y_0_xy_yz[i] * a_exp + 4.0 * g_y_zz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_xy_zz[i] = -2.0 * g_y_0_xy_zz[i] * a_exp + 4.0 * g_y_zz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_y_0_xz_xx, g_y_0_xz_xy, g_y_0_xz_xz, g_y_0_xz_yy, g_y_0_xz_yz, g_y_0_xz_zz, g_y_z_0_0_0_z_xz_xx, g_y_z_0_0_0_z_xz_xy, g_y_z_0_0_0_z_xz_xz, g_y_z_0_0_0_z_xz_yy, g_y_z_0_0_0_z_xz_yz, g_y_z_0_0_0_z_xz_zz, g_y_zz_xz_xx, g_y_zz_xz_xy, g_y_zz_xz_xz, g_y_zz_xz_yy, g_y_zz_xz_yz, g_y_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_z_xz_xx[i] = -2.0 * g_y_0_xz_xx[i] * a_exp + 4.0 * g_y_zz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_xz_xy[i] = -2.0 * g_y_0_xz_xy[i] * a_exp + 4.0 * g_y_zz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_xz_xz[i] = -2.0 * g_y_0_xz_xz[i] * a_exp + 4.0 * g_y_zz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_xz_yy[i] = -2.0 * g_y_0_xz_yy[i] * a_exp + 4.0 * g_y_zz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_xz_yz[i] = -2.0 * g_y_0_xz_yz[i] * a_exp + 4.0 * g_y_zz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_xz_zz[i] = -2.0 * g_y_0_xz_zz[i] * a_exp + 4.0 * g_y_zz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_y_0_yy_xx, g_y_0_yy_xy, g_y_0_yy_xz, g_y_0_yy_yy, g_y_0_yy_yz, g_y_0_yy_zz, g_y_z_0_0_0_z_yy_xx, g_y_z_0_0_0_z_yy_xy, g_y_z_0_0_0_z_yy_xz, g_y_z_0_0_0_z_yy_yy, g_y_z_0_0_0_z_yy_yz, g_y_z_0_0_0_z_yy_zz, g_y_zz_yy_xx, g_y_zz_yy_xy, g_y_zz_yy_xz, g_y_zz_yy_yy, g_y_zz_yy_yz, g_y_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_z_yy_xx[i] = -2.0 * g_y_0_yy_xx[i] * a_exp + 4.0 * g_y_zz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_yy_xy[i] = -2.0 * g_y_0_yy_xy[i] * a_exp + 4.0 * g_y_zz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_yy_xz[i] = -2.0 * g_y_0_yy_xz[i] * a_exp + 4.0 * g_y_zz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_yy_yy[i] = -2.0 * g_y_0_yy_yy[i] * a_exp + 4.0 * g_y_zz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_yy_yz[i] = -2.0 * g_y_0_yy_yz[i] * a_exp + 4.0 * g_y_zz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_yy_zz[i] = -2.0 * g_y_0_yy_zz[i] * a_exp + 4.0 * g_y_zz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_y_0_yz_xx, g_y_0_yz_xy, g_y_0_yz_xz, g_y_0_yz_yy, g_y_0_yz_yz, g_y_0_yz_zz, g_y_z_0_0_0_z_yz_xx, g_y_z_0_0_0_z_yz_xy, g_y_z_0_0_0_z_yz_xz, g_y_z_0_0_0_z_yz_yy, g_y_z_0_0_0_z_yz_yz, g_y_z_0_0_0_z_yz_zz, g_y_zz_yz_xx, g_y_zz_yz_xy, g_y_zz_yz_xz, g_y_zz_yz_yy, g_y_zz_yz_yz, g_y_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_z_yz_xx[i] = -2.0 * g_y_0_yz_xx[i] * a_exp + 4.0 * g_y_zz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_yz_xy[i] = -2.0 * g_y_0_yz_xy[i] * a_exp + 4.0 * g_y_zz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_yz_xz[i] = -2.0 * g_y_0_yz_xz[i] * a_exp + 4.0 * g_y_zz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_yz_yy[i] = -2.0 * g_y_0_yz_yy[i] * a_exp + 4.0 * g_y_zz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_yz_yz[i] = -2.0 * g_y_0_yz_yz[i] * a_exp + 4.0 * g_y_zz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_yz_zz[i] = -2.0 * g_y_0_yz_zz[i] * a_exp + 4.0 * g_y_zz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_y_0_zz_xx, g_y_0_zz_xy, g_y_0_zz_xz, g_y_0_zz_yy, g_y_0_zz_yz, g_y_0_zz_zz, g_y_z_0_0_0_z_zz_xx, g_y_z_0_0_0_z_zz_xy, g_y_z_0_0_0_z_zz_xz, g_y_z_0_0_0_z_zz_yy, g_y_z_0_0_0_z_zz_yz, g_y_z_0_0_0_z_zz_zz, g_y_zz_zz_xx, g_y_zz_zz_xy, g_y_zz_zz_xz, g_y_zz_zz_yy, g_y_zz_zz_yz, g_y_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_z_zz_xx[i] = -2.0 * g_y_0_zz_xx[i] * a_exp + 4.0 * g_y_zz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_zz_xy[i] = -2.0 * g_y_0_zz_xy[i] * a_exp + 4.0 * g_y_zz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_zz_xz[i] = -2.0 * g_y_0_zz_xz[i] * a_exp + 4.0 * g_y_zz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_zz_yy[i] = -2.0 * g_y_0_zz_yy[i] * a_exp + 4.0 * g_y_zz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_zz_yz[i] = -2.0 * g_y_0_zz_yz[i] * a_exp + 4.0 * g_y_zz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_zz_zz[i] = -2.0 * g_y_0_zz_zz[i] * a_exp + 4.0 * g_y_zz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_z_0_xx_xx, g_z_0_xx_xy, g_z_0_xx_xz, g_z_0_xx_yy, g_z_0_xx_yz, g_z_0_xx_zz, g_z_x_0_0_0_x_xx_xx, g_z_x_0_0_0_x_xx_xy, g_z_x_0_0_0_x_xx_xz, g_z_x_0_0_0_x_xx_yy, g_z_x_0_0_0_x_xx_yz, g_z_x_0_0_0_x_xx_zz, g_z_xx_xx_xx, g_z_xx_xx_xy, g_z_xx_xx_xz, g_z_xx_xx_yy, g_z_xx_xx_yz, g_z_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_x_xx_xx[i] = -2.0 * g_z_0_xx_xx[i] * a_exp + 4.0 * g_z_xx_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_xx_xy[i] = -2.0 * g_z_0_xx_xy[i] * a_exp + 4.0 * g_z_xx_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_xx_xz[i] = -2.0 * g_z_0_xx_xz[i] * a_exp + 4.0 * g_z_xx_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_xx_yy[i] = -2.0 * g_z_0_xx_yy[i] * a_exp + 4.0 * g_z_xx_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_xx_yz[i] = -2.0 * g_z_0_xx_yz[i] * a_exp + 4.0 * g_z_xx_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_xx_zz[i] = -2.0 * g_z_0_xx_zz[i] * a_exp + 4.0 * g_z_xx_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_z_0_xy_xx, g_z_0_xy_xy, g_z_0_xy_xz, g_z_0_xy_yy, g_z_0_xy_yz, g_z_0_xy_zz, g_z_x_0_0_0_x_xy_xx, g_z_x_0_0_0_x_xy_xy, g_z_x_0_0_0_x_xy_xz, g_z_x_0_0_0_x_xy_yy, g_z_x_0_0_0_x_xy_yz, g_z_x_0_0_0_x_xy_zz, g_z_xx_xy_xx, g_z_xx_xy_xy, g_z_xx_xy_xz, g_z_xx_xy_yy, g_z_xx_xy_yz, g_z_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_x_xy_xx[i] = -2.0 * g_z_0_xy_xx[i] * a_exp + 4.0 * g_z_xx_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_xy_xy[i] = -2.0 * g_z_0_xy_xy[i] * a_exp + 4.0 * g_z_xx_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_xy_xz[i] = -2.0 * g_z_0_xy_xz[i] * a_exp + 4.0 * g_z_xx_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_xy_yy[i] = -2.0 * g_z_0_xy_yy[i] * a_exp + 4.0 * g_z_xx_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_xy_yz[i] = -2.0 * g_z_0_xy_yz[i] * a_exp + 4.0 * g_z_xx_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_xy_zz[i] = -2.0 * g_z_0_xy_zz[i] * a_exp + 4.0 * g_z_xx_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_z_0_xz_xx, g_z_0_xz_xy, g_z_0_xz_xz, g_z_0_xz_yy, g_z_0_xz_yz, g_z_0_xz_zz, g_z_x_0_0_0_x_xz_xx, g_z_x_0_0_0_x_xz_xy, g_z_x_0_0_0_x_xz_xz, g_z_x_0_0_0_x_xz_yy, g_z_x_0_0_0_x_xz_yz, g_z_x_0_0_0_x_xz_zz, g_z_xx_xz_xx, g_z_xx_xz_xy, g_z_xx_xz_xz, g_z_xx_xz_yy, g_z_xx_xz_yz, g_z_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_x_xz_xx[i] = -2.0 * g_z_0_xz_xx[i] * a_exp + 4.0 * g_z_xx_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_xz_xy[i] = -2.0 * g_z_0_xz_xy[i] * a_exp + 4.0 * g_z_xx_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_xz_xz[i] = -2.0 * g_z_0_xz_xz[i] * a_exp + 4.0 * g_z_xx_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_xz_yy[i] = -2.0 * g_z_0_xz_yy[i] * a_exp + 4.0 * g_z_xx_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_xz_yz[i] = -2.0 * g_z_0_xz_yz[i] * a_exp + 4.0 * g_z_xx_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_xz_zz[i] = -2.0 * g_z_0_xz_zz[i] * a_exp + 4.0 * g_z_xx_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_z_0_yy_xx, g_z_0_yy_xy, g_z_0_yy_xz, g_z_0_yy_yy, g_z_0_yy_yz, g_z_0_yy_zz, g_z_x_0_0_0_x_yy_xx, g_z_x_0_0_0_x_yy_xy, g_z_x_0_0_0_x_yy_xz, g_z_x_0_0_0_x_yy_yy, g_z_x_0_0_0_x_yy_yz, g_z_x_0_0_0_x_yy_zz, g_z_xx_yy_xx, g_z_xx_yy_xy, g_z_xx_yy_xz, g_z_xx_yy_yy, g_z_xx_yy_yz, g_z_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_x_yy_xx[i] = -2.0 * g_z_0_yy_xx[i] * a_exp + 4.0 * g_z_xx_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_yy_xy[i] = -2.0 * g_z_0_yy_xy[i] * a_exp + 4.0 * g_z_xx_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_yy_xz[i] = -2.0 * g_z_0_yy_xz[i] * a_exp + 4.0 * g_z_xx_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_yy_yy[i] = -2.0 * g_z_0_yy_yy[i] * a_exp + 4.0 * g_z_xx_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_yy_yz[i] = -2.0 * g_z_0_yy_yz[i] * a_exp + 4.0 * g_z_xx_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_yy_zz[i] = -2.0 * g_z_0_yy_zz[i] * a_exp + 4.0 * g_z_xx_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_z_0_yz_xx, g_z_0_yz_xy, g_z_0_yz_xz, g_z_0_yz_yy, g_z_0_yz_yz, g_z_0_yz_zz, g_z_x_0_0_0_x_yz_xx, g_z_x_0_0_0_x_yz_xy, g_z_x_0_0_0_x_yz_xz, g_z_x_0_0_0_x_yz_yy, g_z_x_0_0_0_x_yz_yz, g_z_x_0_0_0_x_yz_zz, g_z_xx_yz_xx, g_z_xx_yz_xy, g_z_xx_yz_xz, g_z_xx_yz_yy, g_z_xx_yz_yz, g_z_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_x_yz_xx[i] = -2.0 * g_z_0_yz_xx[i] * a_exp + 4.0 * g_z_xx_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_yz_xy[i] = -2.0 * g_z_0_yz_xy[i] * a_exp + 4.0 * g_z_xx_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_yz_xz[i] = -2.0 * g_z_0_yz_xz[i] * a_exp + 4.0 * g_z_xx_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_yz_yy[i] = -2.0 * g_z_0_yz_yy[i] * a_exp + 4.0 * g_z_xx_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_yz_yz[i] = -2.0 * g_z_0_yz_yz[i] * a_exp + 4.0 * g_z_xx_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_yz_zz[i] = -2.0 * g_z_0_yz_zz[i] * a_exp + 4.0 * g_z_xx_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_z_0_zz_xx, g_z_0_zz_xy, g_z_0_zz_xz, g_z_0_zz_yy, g_z_0_zz_yz, g_z_0_zz_zz, g_z_x_0_0_0_x_zz_xx, g_z_x_0_0_0_x_zz_xy, g_z_x_0_0_0_x_zz_xz, g_z_x_0_0_0_x_zz_yy, g_z_x_0_0_0_x_zz_yz, g_z_x_0_0_0_x_zz_zz, g_z_xx_zz_xx, g_z_xx_zz_xy, g_z_xx_zz_xz, g_z_xx_zz_yy, g_z_xx_zz_yz, g_z_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_x_zz_xx[i] = -2.0 * g_z_0_zz_xx[i] * a_exp + 4.0 * g_z_xx_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_zz_xy[i] = -2.0 * g_z_0_zz_xy[i] * a_exp + 4.0 * g_z_xx_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_zz_xz[i] = -2.0 * g_z_0_zz_xz[i] * a_exp + 4.0 * g_z_xx_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_zz_yy[i] = -2.0 * g_z_0_zz_yy[i] * a_exp + 4.0 * g_z_xx_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_zz_yz[i] = -2.0 * g_z_0_zz_yz[i] * a_exp + 4.0 * g_z_xx_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_zz_zz[i] = -2.0 * g_z_0_zz_zz[i] * a_exp + 4.0 * g_z_xx_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_z_x_0_0_0_y_xx_xx, g_z_x_0_0_0_y_xx_xy, g_z_x_0_0_0_y_xx_xz, g_z_x_0_0_0_y_xx_yy, g_z_x_0_0_0_y_xx_yz, g_z_x_0_0_0_y_xx_zz, g_z_xy_xx_xx, g_z_xy_xx_xy, g_z_xy_xx_xz, g_z_xy_xx_yy, g_z_xy_xx_yz, g_z_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_y_xx_xx[i] = 4.0 * g_z_xy_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_xx_xy[i] = 4.0 * g_z_xy_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_xx_xz[i] = 4.0 * g_z_xy_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_xx_yy[i] = 4.0 * g_z_xy_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_xx_yz[i] = 4.0 * g_z_xy_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_xx_zz[i] = 4.0 * g_z_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_z_x_0_0_0_y_xy_xx, g_z_x_0_0_0_y_xy_xy, g_z_x_0_0_0_y_xy_xz, g_z_x_0_0_0_y_xy_yy, g_z_x_0_0_0_y_xy_yz, g_z_x_0_0_0_y_xy_zz, g_z_xy_xy_xx, g_z_xy_xy_xy, g_z_xy_xy_xz, g_z_xy_xy_yy, g_z_xy_xy_yz, g_z_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_y_xy_xx[i] = 4.0 * g_z_xy_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_xy_xy[i] = 4.0 * g_z_xy_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_xy_xz[i] = 4.0 * g_z_xy_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_xy_yy[i] = 4.0 * g_z_xy_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_xy_yz[i] = 4.0 * g_z_xy_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_xy_zz[i] = 4.0 * g_z_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_z_x_0_0_0_y_xz_xx, g_z_x_0_0_0_y_xz_xy, g_z_x_0_0_0_y_xz_xz, g_z_x_0_0_0_y_xz_yy, g_z_x_0_0_0_y_xz_yz, g_z_x_0_0_0_y_xz_zz, g_z_xy_xz_xx, g_z_xy_xz_xy, g_z_xy_xz_xz, g_z_xy_xz_yy, g_z_xy_xz_yz, g_z_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_y_xz_xx[i] = 4.0 * g_z_xy_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_xz_xy[i] = 4.0 * g_z_xy_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_xz_xz[i] = 4.0 * g_z_xy_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_xz_yy[i] = 4.0 * g_z_xy_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_xz_yz[i] = 4.0 * g_z_xy_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_xz_zz[i] = 4.0 * g_z_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_z_x_0_0_0_y_yy_xx, g_z_x_0_0_0_y_yy_xy, g_z_x_0_0_0_y_yy_xz, g_z_x_0_0_0_y_yy_yy, g_z_x_0_0_0_y_yy_yz, g_z_x_0_0_0_y_yy_zz, g_z_xy_yy_xx, g_z_xy_yy_xy, g_z_xy_yy_xz, g_z_xy_yy_yy, g_z_xy_yy_yz, g_z_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_y_yy_xx[i] = 4.0 * g_z_xy_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_yy_xy[i] = 4.0 * g_z_xy_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_yy_xz[i] = 4.0 * g_z_xy_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_yy_yy[i] = 4.0 * g_z_xy_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_yy_yz[i] = 4.0 * g_z_xy_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_yy_zz[i] = 4.0 * g_z_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_z_x_0_0_0_y_yz_xx, g_z_x_0_0_0_y_yz_xy, g_z_x_0_0_0_y_yz_xz, g_z_x_0_0_0_y_yz_yy, g_z_x_0_0_0_y_yz_yz, g_z_x_0_0_0_y_yz_zz, g_z_xy_yz_xx, g_z_xy_yz_xy, g_z_xy_yz_xz, g_z_xy_yz_yy, g_z_xy_yz_yz, g_z_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_y_yz_xx[i] = 4.0 * g_z_xy_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_yz_xy[i] = 4.0 * g_z_xy_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_yz_xz[i] = 4.0 * g_z_xy_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_yz_yy[i] = 4.0 * g_z_xy_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_yz_yz[i] = 4.0 * g_z_xy_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_yz_zz[i] = 4.0 * g_z_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_z_x_0_0_0_y_zz_xx, g_z_x_0_0_0_y_zz_xy, g_z_x_0_0_0_y_zz_xz, g_z_x_0_0_0_y_zz_yy, g_z_x_0_0_0_y_zz_yz, g_z_x_0_0_0_y_zz_zz, g_z_xy_zz_xx, g_z_xy_zz_xy, g_z_xy_zz_xz, g_z_xy_zz_yy, g_z_xy_zz_yz, g_z_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_y_zz_xx[i] = 4.0 * g_z_xy_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_zz_xy[i] = 4.0 * g_z_xy_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_zz_xz[i] = 4.0 * g_z_xy_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_zz_yy[i] = 4.0 * g_z_xy_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_zz_yz[i] = 4.0 * g_z_xy_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_zz_zz[i] = 4.0 * g_z_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_z_x_0_0_0_z_xx_xx, g_z_x_0_0_0_z_xx_xy, g_z_x_0_0_0_z_xx_xz, g_z_x_0_0_0_z_xx_yy, g_z_x_0_0_0_z_xx_yz, g_z_x_0_0_0_z_xx_zz, g_z_xz_xx_xx, g_z_xz_xx_xy, g_z_xz_xx_xz, g_z_xz_xx_yy, g_z_xz_xx_yz, g_z_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_z_xx_xx[i] = 4.0 * g_z_xz_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_xx_xy[i] = 4.0 * g_z_xz_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_xx_xz[i] = 4.0 * g_z_xz_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_xx_yy[i] = 4.0 * g_z_xz_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_xx_yz[i] = 4.0 * g_z_xz_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_xx_zz[i] = 4.0 * g_z_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_z_x_0_0_0_z_xy_xx, g_z_x_0_0_0_z_xy_xy, g_z_x_0_0_0_z_xy_xz, g_z_x_0_0_0_z_xy_yy, g_z_x_0_0_0_z_xy_yz, g_z_x_0_0_0_z_xy_zz, g_z_xz_xy_xx, g_z_xz_xy_xy, g_z_xz_xy_xz, g_z_xz_xy_yy, g_z_xz_xy_yz, g_z_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_z_xy_xx[i] = 4.0 * g_z_xz_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_xy_xy[i] = 4.0 * g_z_xz_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_xy_xz[i] = 4.0 * g_z_xz_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_xy_yy[i] = 4.0 * g_z_xz_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_xy_yz[i] = 4.0 * g_z_xz_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_xy_zz[i] = 4.0 * g_z_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_z_x_0_0_0_z_xz_xx, g_z_x_0_0_0_z_xz_xy, g_z_x_0_0_0_z_xz_xz, g_z_x_0_0_0_z_xz_yy, g_z_x_0_0_0_z_xz_yz, g_z_x_0_0_0_z_xz_zz, g_z_xz_xz_xx, g_z_xz_xz_xy, g_z_xz_xz_xz, g_z_xz_xz_yy, g_z_xz_xz_yz, g_z_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_z_xz_xx[i] = 4.0 * g_z_xz_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_xz_xy[i] = 4.0 * g_z_xz_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_xz_xz[i] = 4.0 * g_z_xz_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_xz_yy[i] = 4.0 * g_z_xz_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_xz_yz[i] = 4.0 * g_z_xz_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_xz_zz[i] = 4.0 * g_z_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_z_x_0_0_0_z_yy_xx, g_z_x_0_0_0_z_yy_xy, g_z_x_0_0_0_z_yy_xz, g_z_x_0_0_0_z_yy_yy, g_z_x_0_0_0_z_yy_yz, g_z_x_0_0_0_z_yy_zz, g_z_xz_yy_xx, g_z_xz_yy_xy, g_z_xz_yy_xz, g_z_xz_yy_yy, g_z_xz_yy_yz, g_z_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_z_yy_xx[i] = 4.0 * g_z_xz_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_yy_xy[i] = 4.0 * g_z_xz_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_yy_xz[i] = 4.0 * g_z_xz_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_yy_yy[i] = 4.0 * g_z_xz_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_yy_yz[i] = 4.0 * g_z_xz_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_yy_zz[i] = 4.0 * g_z_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_z_x_0_0_0_z_yz_xx, g_z_x_0_0_0_z_yz_xy, g_z_x_0_0_0_z_yz_xz, g_z_x_0_0_0_z_yz_yy, g_z_x_0_0_0_z_yz_yz, g_z_x_0_0_0_z_yz_zz, g_z_xz_yz_xx, g_z_xz_yz_xy, g_z_xz_yz_xz, g_z_xz_yz_yy, g_z_xz_yz_yz, g_z_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_z_yz_xx[i] = 4.0 * g_z_xz_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_yz_xy[i] = 4.0 * g_z_xz_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_yz_xz[i] = 4.0 * g_z_xz_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_yz_yy[i] = 4.0 * g_z_xz_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_yz_yz[i] = 4.0 * g_z_xz_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_yz_zz[i] = 4.0 * g_z_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_z_x_0_0_0_z_zz_xx, g_z_x_0_0_0_z_zz_xy, g_z_x_0_0_0_z_zz_xz, g_z_x_0_0_0_z_zz_yy, g_z_x_0_0_0_z_zz_yz, g_z_x_0_0_0_z_zz_zz, g_z_xz_zz_xx, g_z_xz_zz_xy, g_z_xz_zz_xz, g_z_xz_zz_yy, g_z_xz_zz_yz, g_z_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_z_zz_xx[i] = 4.0 * g_z_xz_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_zz_xy[i] = 4.0 * g_z_xz_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_zz_xz[i] = 4.0 * g_z_xz_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_zz_yy[i] = 4.0 * g_z_xz_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_zz_yz[i] = 4.0 * g_z_xz_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_zz_zz[i] = 4.0 * g_z_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_z_xy_xx_xx, g_z_xy_xx_xy, g_z_xy_xx_xz, g_z_xy_xx_yy, g_z_xy_xx_yz, g_z_xy_xx_zz, g_z_y_0_0_0_x_xx_xx, g_z_y_0_0_0_x_xx_xy, g_z_y_0_0_0_x_xx_xz, g_z_y_0_0_0_x_xx_yy, g_z_y_0_0_0_x_xx_yz, g_z_y_0_0_0_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_x_xx_xx[i] = 4.0 * g_z_xy_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_xx_xy[i] = 4.0 * g_z_xy_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_xx_xz[i] = 4.0 * g_z_xy_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_xx_yy[i] = 4.0 * g_z_xy_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_xx_yz[i] = 4.0 * g_z_xy_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_xx_zz[i] = 4.0 * g_z_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_z_xy_xy_xx, g_z_xy_xy_xy, g_z_xy_xy_xz, g_z_xy_xy_yy, g_z_xy_xy_yz, g_z_xy_xy_zz, g_z_y_0_0_0_x_xy_xx, g_z_y_0_0_0_x_xy_xy, g_z_y_0_0_0_x_xy_xz, g_z_y_0_0_0_x_xy_yy, g_z_y_0_0_0_x_xy_yz, g_z_y_0_0_0_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_x_xy_xx[i] = 4.0 * g_z_xy_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_xy_xy[i] = 4.0 * g_z_xy_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_xy_xz[i] = 4.0 * g_z_xy_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_xy_yy[i] = 4.0 * g_z_xy_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_xy_yz[i] = 4.0 * g_z_xy_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_xy_zz[i] = 4.0 * g_z_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_z_xy_xz_xx, g_z_xy_xz_xy, g_z_xy_xz_xz, g_z_xy_xz_yy, g_z_xy_xz_yz, g_z_xy_xz_zz, g_z_y_0_0_0_x_xz_xx, g_z_y_0_0_0_x_xz_xy, g_z_y_0_0_0_x_xz_xz, g_z_y_0_0_0_x_xz_yy, g_z_y_0_0_0_x_xz_yz, g_z_y_0_0_0_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_x_xz_xx[i] = 4.0 * g_z_xy_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_xz_xy[i] = 4.0 * g_z_xy_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_xz_xz[i] = 4.0 * g_z_xy_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_xz_yy[i] = 4.0 * g_z_xy_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_xz_yz[i] = 4.0 * g_z_xy_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_xz_zz[i] = 4.0 * g_z_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_z_xy_yy_xx, g_z_xy_yy_xy, g_z_xy_yy_xz, g_z_xy_yy_yy, g_z_xy_yy_yz, g_z_xy_yy_zz, g_z_y_0_0_0_x_yy_xx, g_z_y_0_0_0_x_yy_xy, g_z_y_0_0_0_x_yy_xz, g_z_y_0_0_0_x_yy_yy, g_z_y_0_0_0_x_yy_yz, g_z_y_0_0_0_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_x_yy_xx[i] = 4.0 * g_z_xy_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_yy_xy[i] = 4.0 * g_z_xy_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_yy_xz[i] = 4.0 * g_z_xy_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_yy_yy[i] = 4.0 * g_z_xy_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_yy_yz[i] = 4.0 * g_z_xy_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_yy_zz[i] = 4.0 * g_z_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_z_xy_yz_xx, g_z_xy_yz_xy, g_z_xy_yz_xz, g_z_xy_yz_yy, g_z_xy_yz_yz, g_z_xy_yz_zz, g_z_y_0_0_0_x_yz_xx, g_z_y_0_0_0_x_yz_xy, g_z_y_0_0_0_x_yz_xz, g_z_y_0_0_0_x_yz_yy, g_z_y_0_0_0_x_yz_yz, g_z_y_0_0_0_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_x_yz_xx[i] = 4.0 * g_z_xy_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_yz_xy[i] = 4.0 * g_z_xy_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_yz_xz[i] = 4.0 * g_z_xy_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_yz_yy[i] = 4.0 * g_z_xy_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_yz_yz[i] = 4.0 * g_z_xy_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_yz_zz[i] = 4.0 * g_z_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_z_xy_zz_xx, g_z_xy_zz_xy, g_z_xy_zz_xz, g_z_xy_zz_yy, g_z_xy_zz_yz, g_z_xy_zz_zz, g_z_y_0_0_0_x_zz_xx, g_z_y_0_0_0_x_zz_xy, g_z_y_0_0_0_x_zz_xz, g_z_y_0_0_0_x_zz_yy, g_z_y_0_0_0_x_zz_yz, g_z_y_0_0_0_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_x_zz_xx[i] = 4.0 * g_z_xy_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_zz_xy[i] = 4.0 * g_z_xy_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_zz_xz[i] = 4.0 * g_z_xy_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_zz_yy[i] = 4.0 * g_z_xy_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_zz_yz[i] = 4.0 * g_z_xy_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_zz_zz[i] = 4.0 * g_z_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_z_0_xx_xx, g_z_0_xx_xy, g_z_0_xx_xz, g_z_0_xx_yy, g_z_0_xx_yz, g_z_0_xx_zz, g_z_y_0_0_0_y_xx_xx, g_z_y_0_0_0_y_xx_xy, g_z_y_0_0_0_y_xx_xz, g_z_y_0_0_0_y_xx_yy, g_z_y_0_0_0_y_xx_yz, g_z_y_0_0_0_y_xx_zz, g_z_yy_xx_xx, g_z_yy_xx_xy, g_z_yy_xx_xz, g_z_yy_xx_yy, g_z_yy_xx_yz, g_z_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_y_xx_xx[i] = -2.0 * g_z_0_xx_xx[i] * a_exp + 4.0 * g_z_yy_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_xx_xy[i] = -2.0 * g_z_0_xx_xy[i] * a_exp + 4.0 * g_z_yy_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_xx_xz[i] = -2.0 * g_z_0_xx_xz[i] * a_exp + 4.0 * g_z_yy_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_xx_yy[i] = -2.0 * g_z_0_xx_yy[i] * a_exp + 4.0 * g_z_yy_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_xx_yz[i] = -2.0 * g_z_0_xx_yz[i] * a_exp + 4.0 * g_z_yy_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_xx_zz[i] = -2.0 * g_z_0_xx_zz[i] * a_exp + 4.0 * g_z_yy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_z_0_xy_xx, g_z_0_xy_xy, g_z_0_xy_xz, g_z_0_xy_yy, g_z_0_xy_yz, g_z_0_xy_zz, g_z_y_0_0_0_y_xy_xx, g_z_y_0_0_0_y_xy_xy, g_z_y_0_0_0_y_xy_xz, g_z_y_0_0_0_y_xy_yy, g_z_y_0_0_0_y_xy_yz, g_z_y_0_0_0_y_xy_zz, g_z_yy_xy_xx, g_z_yy_xy_xy, g_z_yy_xy_xz, g_z_yy_xy_yy, g_z_yy_xy_yz, g_z_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_y_xy_xx[i] = -2.0 * g_z_0_xy_xx[i] * a_exp + 4.0 * g_z_yy_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_xy_xy[i] = -2.0 * g_z_0_xy_xy[i] * a_exp + 4.0 * g_z_yy_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_xy_xz[i] = -2.0 * g_z_0_xy_xz[i] * a_exp + 4.0 * g_z_yy_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_xy_yy[i] = -2.0 * g_z_0_xy_yy[i] * a_exp + 4.0 * g_z_yy_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_xy_yz[i] = -2.0 * g_z_0_xy_yz[i] * a_exp + 4.0 * g_z_yy_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_xy_zz[i] = -2.0 * g_z_0_xy_zz[i] * a_exp + 4.0 * g_z_yy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_z_0_xz_xx, g_z_0_xz_xy, g_z_0_xz_xz, g_z_0_xz_yy, g_z_0_xz_yz, g_z_0_xz_zz, g_z_y_0_0_0_y_xz_xx, g_z_y_0_0_0_y_xz_xy, g_z_y_0_0_0_y_xz_xz, g_z_y_0_0_0_y_xz_yy, g_z_y_0_0_0_y_xz_yz, g_z_y_0_0_0_y_xz_zz, g_z_yy_xz_xx, g_z_yy_xz_xy, g_z_yy_xz_xz, g_z_yy_xz_yy, g_z_yy_xz_yz, g_z_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_y_xz_xx[i] = -2.0 * g_z_0_xz_xx[i] * a_exp + 4.0 * g_z_yy_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_xz_xy[i] = -2.0 * g_z_0_xz_xy[i] * a_exp + 4.0 * g_z_yy_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_xz_xz[i] = -2.0 * g_z_0_xz_xz[i] * a_exp + 4.0 * g_z_yy_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_xz_yy[i] = -2.0 * g_z_0_xz_yy[i] * a_exp + 4.0 * g_z_yy_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_xz_yz[i] = -2.0 * g_z_0_xz_yz[i] * a_exp + 4.0 * g_z_yy_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_xz_zz[i] = -2.0 * g_z_0_xz_zz[i] * a_exp + 4.0 * g_z_yy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_z_0_yy_xx, g_z_0_yy_xy, g_z_0_yy_xz, g_z_0_yy_yy, g_z_0_yy_yz, g_z_0_yy_zz, g_z_y_0_0_0_y_yy_xx, g_z_y_0_0_0_y_yy_xy, g_z_y_0_0_0_y_yy_xz, g_z_y_0_0_0_y_yy_yy, g_z_y_0_0_0_y_yy_yz, g_z_y_0_0_0_y_yy_zz, g_z_yy_yy_xx, g_z_yy_yy_xy, g_z_yy_yy_xz, g_z_yy_yy_yy, g_z_yy_yy_yz, g_z_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_y_yy_xx[i] = -2.0 * g_z_0_yy_xx[i] * a_exp + 4.0 * g_z_yy_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_yy_xy[i] = -2.0 * g_z_0_yy_xy[i] * a_exp + 4.0 * g_z_yy_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_yy_xz[i] = -2.0 * g_z_0_yy_xz[i] * a_exp + 4.0 * g_z_yy_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_yy_yy[i] = -2.0 * g_z_0_yy_yy[i] * a_exp + 4.0 * g_z_yy_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_yy_yz[i] = -2.0 * g_z_0_yy_yz[i] * a_exp + 4.0 * g_z_yy_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_yy_zz[i] = -2.0 * g_z_0_yy_zz[i] * a_exp + 4.0 * g_z_yy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_z_0_yz_xx, g_z_0_yz_xy, g_z_0_yz_xz, g_z_0_yz_yy, g_z_0_yz_yz, g_z_0_yz_zz, g_z_y_0_0_0_y_yz_xx, g_z_y_0_0_0_y_yz_xy, g_z_y_0_0_0_y_yz_xz, g_z_y_0_0_0_y_yz_yy, g_z_y_0_0_0_y_yz_yz, g_z_y_0_0_0_y_yz_zz, g_z_yy_yz_xx, g_z_yy_yz_xy, g_z_yy_yz_xz, g_z_yy_yz_yy, g_z_yy_yz_yz, g_z_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_y_yz_xx[i] = -2.0 * g_z_0_yz_xx[i] * a_exp + 4.0 * g_z_yy_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_yz_xy[i] = -2.0 * g_z_0_yz_xy[i] * a_exp + 4.0 * g_z_yy_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_yz_xz[i] = -2.0 * g_z_0_yz_xz[i] * a_exp + 4.0 * g_z_yy_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_yz_yy[i] = -2.0 * g_z_0_yz_yy[i] * a_exp + 4.0 * g_z_yy_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_yz_yz[i] = -2.0 * g_z_0_yz_yz[i] * a_exp + 4.0 * g_z_yy_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_yz_zz[i] = -2.0 * g_z_0_yz_zz[i] * a_exp + 4.0 * g_z_yy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_z_0_zz_xx, g_z_0_zz_xy, g_z_0_zz_xz, g_z_0_zz_yy, g_z_0_zz_yz, g_z_0_zz_zz, g_z_y_0_0_0_y_zz_xx, g_z_y_0_0_0_y_zz_xy, g_z_y_0_0_0_y_zz_xz, g_z_y_0_0_0_y_zz_yy, g_z_y_0_0_0_y_zz_yz, g_z_y_0_0_0_y_zz_zz, g_z_yy_zz_xx, g_z_yy_zz_xy, g_z_yy_zz_xz, g_z_yy_zz_yy, g_z_yy_zz_yz, g_z_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_y_zz_xx[i] = -2.0 * g_z_0_zz_xx[i] * a_exp + 4.0 * g_z_yy_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_zz_xy[i] = -2.0 * g_z_0_zz_xy[i] * a_exp + 4.0 * g_z_yy_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_zz_xz[i] = -2.0 * g_z_0_zz_xz[i] * a_exp + 4.0 * g_z_yy_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_zz_yy[i] = -2.0 * g_z_0_zz_yy[i] * a_exp + 4.0 * g_z_yy_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_zz_yz[i] = -2.0 * g_z_0_zz_yz[i] * a_exp + 4.0 * g_z_yy_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_zz_zz[i] = -2.0 * g_z_0_zz_zz[i] * a_exp + 4.0 * g_z_yy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_z_y_0_0_0_z_xx_xx, g_z_y_0_0_0_z_xx_xy, g_z_y_0_0_0_z_xx_xz, g_z_y_0_0_0_z_xx_yy, g_z_y_0_0_0_z_xx_yz, g_z_y_0_0_0_z_xx_zz, g_z_yz_xx_xx, g_z_yz_xx_xy, g_z_yz_xx_xz, g_z_yz_xx_yy, g_z_yz_xx_yz, g_z_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_z_xx_xx[i] = 4.0 * g_z_yz_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_xx_xy[i] = 4.0 * g_z_yz_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_xx_xz[i] = 4.0 * g_z_yz_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_xx_yy[i] = 4.0 * g_z_yz_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_xx_yz[i] = 4.0 * g_z_yz_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_xx_zz[i] = 4.0 * g_z_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_z_y_0_0_0_z_xy_xx, g_z_y_0_0_0_z_xy_xy, g_z_y_0_0_0_z_xy_xz, g_z_y_0_0_0_z_xy_yy, g_z_y_0_0_0_z_xy_yz, g_z_y_0_0_0_z_xy_zz, g_z_yz_xy_xx, g_z_yz_xy_xy, g_z_yz_xy_xz, g_z_yz_xy_yy, g_z_yz_xy_yz, g_z_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_z_xy_xx[i] = 4.0 * g_z_yz_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_xy_xy[i] = 4.0 * g_z_yz_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_xy_xz[i] = 4.0 * g_z_yz_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_xy_yy[i] = 4.0 * g_z_yz_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_xy_yz[i] = 4.0 * g_z_yz_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_xy_zz[i] = 4.0 * g_z_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_z_y_0_0_0_z_xz_xx, g_z_y_0_0_0_z_xz_xy, g_z_y_0_0_0_z_xz_xz, g_z_y_0_0_0_z_xz_yy, g_z_y_0_0_0_z_xz_yz, g_z_y_0_0_0_z_xz_zz, g_z_yz_xz_xx, g_z_yz_xz_xy, g_z_yz_xz_xz, g_z_yz_xz_yy, g_z_yz_xz_yz, g_z_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_z_xz_xx[i] = 4.0 * g_z_yz_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_xz_xy[i] = 4.0 * g_z_yz_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_xz_xz[i] = 4.0 * g_z_yz_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_xz_yy[i] = 4.0 * g_z_yz_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_xz_yz[i] = 4.0 * g_z_yz_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_xz_zz[i] = 4.0 * g_z_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_z_y_0_0_0_z_yy_xx, g_z_y_0_0_0_z_yy_xy, g_z_y_0_0_0_z_yy_xz, g_z_y_0_0_0_z_yy_yy, g_z_y_0_0_0_z_yy_yz, g_z_y_0_0_0_z_yy_zz, g_z_yz_yy_xx, g_z_yz_yy_xy, g_z_yz_yy_xz, g_z_yz_yy_yy, g_z_yz_yy_yz, g_z_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_z_yy_xx[i] = 4.0 * g_z_yz_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_yy_xy[i] = 4.0 * g_z_yz_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_yy_xz[i] = 4.0 * g_z_yz_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_yy_yy[i] = 4.0 * g_z_yz_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_yy_yz[i] = 4.0 * g_z_yz_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_yy_zz[i] = 4.0 * g_z_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_z_y_0_0_0_z_yz_xx, g_z_y_0_0_0_z_yz_xy, g_z_y_0_0_0_z_yz_xz, g_z_y_0_0_0_z_yz_yy, g_z_y_0_0_0_z_yz_yz, g_z_y_0_0_0_z_yz_zz, g_z_yz_yz_xx, g_z_yz_yz_xy, g_z_yz_yz_xz, g_z_yz_yz_yy, g_z_yz_yz_yz, g_z_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_z_yz_xx[i] = 4.0 * g_z_yz_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_yz_xy[i] = 4.0 * g_z_yz_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_yz_xz[i] = 4.0 * g_z_yz_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_yz_yy[i] = 4.0 * g_z_yz_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_yz_yz[i] = 4.0 * g_z_yz_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_yz_zz[i] = 4.0 * g_z_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_z_y_0_0_0_z_zz_xx, g_z_y_0_0_0_z_zz_xy, g_z_y_0_0_0_z_zz_xz, g_z_y_0_0_0_z_zz_yy, g_z_y_0_0_0_z_zz_yz, g_z_y_0_0_0_z_zz_zz, g_z_yz_zz_xx, g_z_yz_zz_xy, g_z_yz_zz_xz, g_z_yz_zz_yy, g_z_yz_zz_yz, g_z_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_z_zz_xx[i] = 4.0 * g_z_yz_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_zz_xy[i] = 4.0 * g_z_yz_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_zz_xz[i] = 4.0 * g_z_yz_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_zz_yy[i] = 4.0 * g_z_yz_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_zz_yz[i] = 4.0 * g_z_yz_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_zz_zz[i] = 4.0 * g_z_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_z_xz_xx_xx, g_z_xz_xx_xy, g_z_xz_xx_xz, g_z_xz_xx_yy, g_z_xz_xx_yz, g_z_xz_xx_zz, g_z_z_0_0_0_x_xx_xx, g_z_z_0_0_0_x_xx_xy, g_z_z_0_0_0_x_xx_xz, g_z_z_0_0_0_x_xx_yy, g_z_z_0_0_0_x_xx_yz, g_z_z_0_0_0_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_x_xx_xx[i] = 4.0 * g_z_xz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_xx_xy[i] = 4.0 * g_z_xz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_xx_xz[i] = 4.0 * g_z_xz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_xx_yy[i] = 4.0 * g_z_xz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_xx_yz[i] = 4.0 * g_z_xz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_xx_zz[i] = 4.0 * g_z_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_z_xz_xy_xx, g_z_xz_xy_xy, g_z_xz_xy_xz, g_z_xz_xy_yy, g_z_xz_xy_yz, g_z_xz_xy_zz, g_z_z_0_0_0_x_xy_xx, g_z_z_0_0_0_x_xy_xy, g_z_z_0_0_0_x_xy_xz, g_z_z_0_0_0_x_xy_yy, g_z_z_0_0_0_x_xy_yz, g_z_z_0_0_0_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_x_xy_xx[i] = 4.0 * g_z_xz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_xy_xy[i] = 4.0 * g_z_xz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_xy_xz[i] = 4.0 * g_z_xz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_xy_yy[i] = 4.0 * g_z_xz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_xy_yz[i] = 4.0 * g_z_xz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_xy_zz[i] = 4.0 * g_z_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_z_xz_xz_xx, g_z_xz_xz_xy, g_z_xz_xz_xz, g_z_xz_xz_yy, g_z_xz_xz_yz, g_z_xz_xz_zz, g_z_z_0_0_0_x_xz_xx, g_z_z_0_0_0_x_xz_xy, g_z_z_0_0_0_x_xz_xz, g_z_z_0_0_0_x_xz_yy, g_z_z_0_0_0_x_xz_yz, g_z_z_0_0_0_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_x_xz_xx[i] = 4.0 * g_z_xz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_xz_xy[i] = 4.0 * g_z_xz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_xz_xz[i] = 4.0 * g_z_xz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_xz_yy[i] = 4.0 * g_z_xz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_xz_yz[i] = 4.0 * g_z_xz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_xz_zz[i] = 4.0 * g_z_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_z_xz_yy_xx, g_z_xz_yy_xy, g_z_xz_yy_xz, g_z_xz_yy_yy, g_z_xz_yy_yz, g_z_xz_yy_zz, g_z_z_0_0_0_x_yy_xx, g_z_z_0_0_0_x_yy_xy, g_z_z_0_0_0_x_yy_xz, g_z_z_0_0_0_x_yy_yy, g_z_z_0_0_0_x_yy_yz, g_z_z_0_0_0_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_x_yy_xx[i] = 4.0 * g_z_xz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_yy_xy[i] = 4.0 * g_z_xz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_yy_xz[i] = 4.0 * g_z_xz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_yy_yy[i] = 4.0 * g_z_xz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_yy_yz[i] = 4.0 * g_z_xz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_yy_zz[i] = 4.0 * g_z_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_z_xz_yz_xx, g_z_xz_yz_xy, g_z_xz_yz_xz, g_z_xz_yz_yy, g_z_xz_yz_yz, g_z_xz_yz_zz, g_z_z_0_0_0_x_yz_xx, g_z_z_0_0_0_x_yz_xy, g_z_z_0_0_0_x_yz_xz, g_z_z_0_0_0_x_yz_yy, g_z_z_0_0_0_x_yz_yz, g_z_z_0_0_0_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_x_yz_xx[i] = 4.0 * g_z_xz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_yz_xy[i] = 4.0 * g_z_xz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_yz_xz[i] = 4.0 * g_z_xz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_yz_yy[i] = 4.0 * g_z_xz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_yz_yz[i] = 4.0 * g_z_xz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_yz_zz[i] = 4.0 * g_z_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_z_xz_zz_xx, g_z_xz_zz_xy, g_z_xz_zz_xz, g_z_xz_zz_yy, g_z_xz_zz_yz, g_z_xz_zz_zz, g_z_z_0_0_0_x_zz_xx, g_z_z_0_0_0_x_zz_xy, g_z_z_0_0_0_x_zz_xz, g_z_z_0_0_0_x_zz_yy, g_z_z_0_0_0_x_zz_yz, g_z_z_0_0_0_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_x_zz_xx[i] = 4.0 * g_z_xz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_zz_xy[i] = 4.0 * g_z_xz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_zz_xz[i] = 4.0 * g_z_xz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_zz_yy[i] = 4.0 * g_z_xz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_zz_yz[i] = 4.0 * g_z_xz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_zz_zz[i] = 4.0 * g_z_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_z_yz_xx_xx, g_z_yz_xx_xy, g_z_yz_xx_xz, g_z_yz_xx_yy, g_z_yz_xx_yz, g_z_yz_xx_zz, g_z_z_0_0_0_y_xx_xx, g_z_z_0_0_0_y_xx_xy, g_z_z_0_0_0_y_xx_xz, g_z_z_0_0_0_y_xx_yy, g_z_z_0_0_0_y_xx_yz, g_z_z_0_0_0_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_y_xx_xx[i] = 4.0 * g_z_yz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_xx_xy[i] = 4.0 * g_z_yz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_xx_xz[i] = 4.0 * g_z_yz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_xx_yy[i] = 4.0 * g_z_yz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_xx_yz[i] = 4.0 * g_z_yz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_xx_zz[i] = 4.0 * g_z_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_z_yz_xy_xx, g_z_yz_xy_xy, g_z_yz_xy_xz, g_z_yz_xy_yy, g_z_yz_xy_yz, g_z_yz_xy_zz, g_z_z_0_0_0_y_xy_xx, g_z_z_0_0_0_y_xy_xy, g_z_z_0_0_0_y_xy_xz, g_z_z_0_0_0_y_xy_yy, g_z_z_0_0_0_y_xy_yz, g_z_z_0_0_0_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_y_xy_xx[i] = 4.0 * g_z_yz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_xy_xy[i] = 4.0 * g_z_yz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_xy_xz[i] = 4.0 * g_z_yz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_xy_yy[i] = 4.0 * g_z_yz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_xy_yz[i] = 4.0 * g_z_yz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_xy_zz[i] = 4.0 * g_z_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_z_yz_xz_xx, g_z_yz_xz_xy, g_z_yz_xz_xz, g_z_yz_xz_yy, g_z_yz_xz_yz, g_z_yz_xz_zz, g_z_z_0_0_0_y_xz_xx, g_z_z_0_0_0_y_xz_xy, g_z_z_0_0_0_y_xz_xz, g_z_z_0_0_0_y_xz_yy, g_z_z_0_0_0_y_xz_yz, g_z_z_0_0_0_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_y_xz_xx[i] = 4.0 * g_z_yz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_xz_xy[i] = 4.0 * g_z_yz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_xz_xz[i] = 4.0 * g_z_yz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_xz_yy[i] = 4.0 * g_z_yz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_xz_yz[i] = 4.0 * g_z_yz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_xz_zz[i] = 4.0 * g_z_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_z_yz_yy_xx, g_z_yz_yy_xy, g_z_yz_yy_xz, g_z_yz_yy_yy, g_z_yz_yy_yz, g_z_yz_yy_zz, g_z_z_0_0_0_y_yy_xx, g_z_z_0_0_0_y_yy_xy, g_z_z_0_0_0_y_yy_xz, g_z_z_0_0_0_y_yy_yy, g_z_z_0_0_0_y_yy_yz, g_z_z_0_0_0_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_y_yy_xx[i] = 4.0 * g_z_yz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_yy_xy[i] = 4.0 * g_z_yz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_yy_xz[i] = 4.0 * g_z_yz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_yy_yy[i] = 4.0 * g_z_yz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_yy_yz[i] = 4.0 * g_z_yz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_yy_zz[i] = 4.0 * g_z_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_z_yz_yz_xx, g_z_yz_yz_xy, g_z_yz_yz_xz, g_z_yz_yz_yy, g_z_yz_yz_yz, g_z_yz_yz_zz, g_z_z_0_0_0_y_yz_xx, g_z_z_0_0_0_y_yz_xy, g_z_z_0_0_0_y_yz_xz, g_z_z_0_0_0_y_yz_yy, g_z_z_0_0_0_y_yz_yz, g_z_z_0_0_0_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_y_yz_xx[i] = 4.0 * g_z_yz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_yz_xy[i] = 4.0 * g_z_yz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_yz_xz[i] = 4.0 * g_z_yz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_yz_yy[i] = 4.0 * g_z_yz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_yz_yz[i] = 4.0 * g_z_yz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_yz_zz[i] = 4.0 * g_z_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_z_yz_zz_xx, g_z_yz_zz_xy, g_z_yz_zz_xz, g_z_yz_zz_yy, g_z_yz_zz_yz, g_z_yz_zz_zz, g_z_z_0_0_0_y_zz_xx, g_z_z_0_0_0_y_zz_xy, g_z_z_0_0_0_y_zz_xz, g_z_z_0_0_0_y_zz_yy, g_z_z_0_0_0_y_zz_yz, g_z_z_0_0_0_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_y_zz_xx[i] = 4.0 * g_z_yz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_zz_xy[i] = 4.0 * g_z_yz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_zz_xz[i] = 4.0 * g_z_yz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_zz_yy[i] = 4.0 * g_z_yz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_zz_yz[i] = 4.0 * g_z_yz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_zz_zz[i] = 4.0 * g_z_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_z_0_xx_xx, g_z_0_xx_xy, g_z_0_xx_xz, g_z_0_xx_yy, g_z_0_xx_yz, g_z_0_xx_zz, g_z_z_0_0_0_z_xx_xx, g_z_z_0_0_0_z_xx_xy, g_z_z_0_0_0_z_xx_xz, g_z_z_0_0_0_z_xx_yy, g_z_z_0_0_0_z_xx_yz, g_z_z_0_0_0_z_xx_zz, g_z_zz_xx_xx, g_z_zz_xx_xy, g_z_zz_xx_xz, g_z_zz_xx_yy, g_z_zz_xx_yz, g_z_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_z_xx_xx[i] = -2.0 * g_z_0_xx_xx[i] * a_exp + 4.0 * g_z_zz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_xx_xy[i] = -2.0 * g_z_0_xx_xy[i] * a_exp + 4.0 * g_z_zz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_xx_xz[i] = -2.0 * g_z_0_xx_xz[i] * a_exp + 4.0 * g_z_zz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_xx_yy[i] = -2.0 * g_z_0_xx_yy[i] * a_exp + 4.0 * g_z_zz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_xx_yz[i] = -2.0 * g_z_0_xx_yz[i] * a_exp + 4.0 * g_z_zz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_xx_zz[i] = -2.0 * g_z_0_xx_zz[i] * a_exp + 4.0 * g_z_zz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_z_0_xy_xx, g_z_0_xy_xy, g_z_0_xy_xz, g_z_0_xy_yy, g_z_0_xy_yz, g_z_0_xy_zz, g_z_z_0_0_0_z_xy_xx, g_z_z_0_0_0_z_xy_xy, g_z_z_0_0_0_z_xy_xz, g_z_z_0_0_0_z_xy_yy, g_z_z_0_0_0_z_xy_yz, g_z_z_0_0_0_z_xy_zz, g_z_zz_xy_xx, g_z_zz_xy_xy, g_z_zz_xy_xz, g_z_zz_xy_yy, g_z_zz_xy_yz, g_z_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_z_xy_xx[i] = -2.0 * g_z_0_xy_xx[i] * a_exp + 4.0 * g_z_zz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_xy_xy[i] = -2.0 * g_z_0_xy_xy[i] * a_exp + 4.0 * g_z_zz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_xy_xz[i] = -2.0 * g_z_0_xy_xz[i] * a_exp + 4.0 * g_z_zz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_xy_yy[i] = -2.0 * g_z_0_xy_yy[i] * a_exp + 4.0 * g_z_zz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_xy_yz[i] = -2.0 * g_z_0_xy_yz[i] * a_exp + 4.0 * g_z_zz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_xy_zz[i] = -2.0 * g_z_0_xy_zz[i] * a_exp + 4.0 * g_z_zz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_z_0_xz_xx, g_z_0_xz_xy, g_z_0_xz_xz, g_z_0_xz_yy, g_z_0_xz_yz, g_z_0_xz_zz, g_z_z_0_0_0_z_xz_xx, g_z_z_0_0_0_z_xz_xy, g_z_z_0_0_0_z_xz_xz, g_z_z_0_0_0_z_xz_yy, g_z_z_0_0_0_z_xz_yz, g_z_z_0_0_0_z_xz_zz, g_z_zz_xz_xx, g_z_zz_xz_xy, g_z_zz_xz_xz, g_z_zz_xz_yy, g_z_zz_xz_yz, g_z_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_z_xz_xx[i] = -2.0 * g_z_0_xz_xx[i] * a_exp + 4.0 * g_z_zz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_xz_xy[i] = -2.0 * g_z_0_xz_xy[i] * a_exp + 4.0 * g_z_zz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_xz_xz[i] = -2.0 * g_z_0_xz_xz[i] * a_exp + 4.0 * g_z_zz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_xz_yy[i] = -2.0 * g_z_0_xz_yy[i] * a_exp + 4.0 * g_z_zz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_xz_yz[i] = -2.0 * g_z_0_xz_yz[i] * a_exp + 4.0 * g_z_zz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_xz_zz[i] = -2.0 * g_z_0_xz_zz[i] * a_exp + 4.0 * g_z_zz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_z_0_yy_xx, g_z_0_yy_xy, g_z_0_yy_xz, g_z_0_yy_yy, g_z_0_yy_yz, g_z_0_yy_zz, g_z_z_0_0_0_z_yy_xx, g_z_z_0_0_0_z_yy_xy, g_z_z_0_0_0_z_yy_xz, g_z_z_0_0_0_z_yy_yy, g_z_z_0_0_0_z_yy_yz, g_z_z_0_0_0_z_yy_zz, g_z_zz_yy_xx, g_z_zz_yy_xy, g_z_zz_yy_xz, g_z_zz_yy_yy, g_z_zz_yy_yz, g_z_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_z_yy_xx[i] = -2.0 * g_z_0_yy_xx[i] * a_exp + 4.0 * g_z_zz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_yy_xy[i] = -2.0 * g_z_0_yy_xy[i] * a_exp + 4.0 * g_z_zz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_yy_xz[i] = -2.0 * g_z_0_yy_xz[i] * a_exp + 4.0 * g_z_zz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_yy_yy[i] = -2.0 * g_z_0_yy_yy[i] * a_exp + 4.0 * g_z_zz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_yy_yz[i] = -2.0 * g_z_0_yy_yz[i] * a_exp + 4.0 * g_z_zz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_yy_zz[i] = -2.0 * g_z_0_yy_zz[i] * a_exp + 4.0 * g_z_zz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_z_0_yz_xx, g_z_0_yz_xy, g_z_0_yz_xz, g_z_0_yz_yy, g_z_0_yz_yz, g_z_0_yz_zz, g_z_z_0_0_0_z_yz_xx, g_z_z_0_0_0_z_yz_xy, g_z_z_0_0_0_z_yz_xz, g_z_z_0_0_0_z_yz_yy, g_z_z_0_0_0_z_yz_yz, g_z_z_0_0_0_z_yz_zz, g_z_zz_yz_xx, g_z_zz_yz_xy, g_z_zz_yz_xz, g_z_zz_yz_yy, g_z_zz_yz_yz, g_z_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_z_yz_xx[i] = -2.0 * g_z_0_yz_xx[i] * a_exp + 4.0 * g_z_zz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_yz_xy[i] = -2.0 * g_z_0_yz_xy[i] * a_exp + 4.0 * g_z_zz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_yz_xz[i] = -2.0 * g_z_0_yz_xz[i] * a_exp + 4.0 * g_z_zz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_yz_yy[i] = -2.0 * g_z_0_yz_yy[i] * a_exp + 4.0 * g_z_zz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_yz_yz[i] = -2.0 * g_z_0_yz_yz[i] * a_exp + 4.0 * g_z_zz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_yz_zz[i] = -2.0 * g_z_0_yz_zz[i] * a_exp + 4.0 * g_z_zz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_z_0_zz_xx, g_z_0_zz_xy, g_z_0_zz_xz, g_z_0_zz_yy, g_z_0_zz_yz, g_z_0_zz_zz, g_z_z_0_0_0_z_zz_xx, g_z_z_0_0_0_z_zz_xy, g_z_z_0_0_0_z_zz_xz, g_z_z_0_0_0_z_zz_yy, g_z_z_0_0_0_z_zz_yz, g_z_z_0_0_0_z_zz_zz, g_z_zz_zz_xx, g_z_zz_zz_xy, g_z_zz_zz_xz, g_z_zz_zz_yy, g_z_zz_zz_yz, g_z_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_z_zz_xx[i] = -2.0 * g_z_0_zz_xx[i] * a_exp + 4.0 * g_z_zz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_zz_xy[i] = -2.0 * g_z_0_zz_xy[i] * a_exp + 4.0 * g_z_zz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_zz_xz[i] = -2.0 * g_z_0_zz_xz[i] * a_exp + 4.0 * g_z_zz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_zz_yy[i] = -2.0 * g_z_0_zz_yy[i] * a_exp + 4.0 * g_z_zz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_zz_yz[i] = -2.0 * g_z_0_zz_yz[i] * a_exp + 4.0 * g_z_zz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_zz_zz[i] = -2.0 * g_z_0_zz_zz[i] * a_exp + 4.0 * g_z_zz_zz_zz[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

