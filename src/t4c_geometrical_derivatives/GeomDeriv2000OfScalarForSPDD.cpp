#include "GeomDeriv2000OfScalarForSPDD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_spdd_0(CSimdArray<double>& buffer_2000_spdd,
                     const CSimdArray<double>& buffer_spdd,
                     const CSimdArray<double>& buffer_dpdd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_spdd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_spdd

    auto g_0_x_xx_xx = buffer_spdd[0];

    auto g_0_x_xx_xy = buffer_spdd[1];

    auto g_0_x_xx_xz = buffer_spdd[2];

    auto g_0_x_xx_yy = buffer_spdd[3];

    auto g_0_x_xx_yz = buffer_spdd[4];

    auto g_0_x_xx_zz = buffer_spdd[5];

    auto g_0_x_xy_xx = buffer_spdd[6];

    auto g_0_x_xy_xy = buffer_spdd[7];

    auto g_0_x_xy_xz = buffer_spdd[8];

    auto g_0_x_xy_yy = buffer_spdd[9];

    auto g_0_x_xy_yz = buffer_spdd[10];

    auto g_0_x_xy_zz = buffer_spdd[11];

    auto g_0_x_xz_xx = buffer_spdd[12];

    auto g_0_x_xz_xy = buffer_spdd[13];

    auto g_0_x_xz_xz = buffer_spdd[14];

    auto g_0_x_xz_yy = buffer_spdd[15];

    auto g_0_x_xz_yz = buffer_spdd[16];

    auto g_0_x_xz_zz = buffer_spdd[17];

    auto g_0_x_yy_xx = buffer_spdd[18];

    auto g_0_x_yy_xy = buffer_spdd[19];

    auto g_0_x_yy_xz = buffer_spdd[20];

    auto g_0_x_yy_yy = buffer_spdd[21];

    auto g_0_x_yy_yz = buffer_spdd[22];

    auto g_0_x_yy_zz = buffer_spdd[23];

    auto g_0_x_yz_xx = buffer_spdd[24];

    auto g_0_x_yz_xy = buffer_spdd[25];

    auto g_0_x_yz_xz = buffer_spdd[26];

    auto g_0_x_yz_yy = buffer_spdd[27];

    auto g_0_x_yz_yz = buffer_spdd[28];

    auto g_0_x_yz_zz = buffer_spdd[29];

    auto g_0_x_zz_xx = buffer_spdd[30];

    auto g_0_x_zz_xy = buffer_spdd[31];

    auto g_0_x_zz_xz = buffer_spdd[32];

    auto g_0_x_zz_yy = buffer_spdd[33];

    auto g_0_x_zz_yz = buffer_spdd[34];

    auto g_0_x_zz_zz = buffer_spdd[35];

    auto g_0_y_xx_xx = buffer_spdd[36];

    auto g_0_y_xx_xy = buffer_spdd[37];

    auto g_0_y_xx_xz = buffer_spdd[38];

    auto g_0_y_xx_yy = buffer_spdd[39];

    auto g_0_y_xx_yz = buffer_spdd[40];

    auto g_0_y_xx_zz = buffer_spdd[41];

    auto g_0_y_xy_xx = buffer_spdd[42];

    auto g_0_y_xy_xy = buffer_spdd[43];

    auto g_0_y_xy_xz = buffer_spdd[44];

    auto g_0_y_xy_yy = buffer_spdd[45];

    auto g_0_y_xy_yz = buffer_spdd[46];

    auto g_0_y_xy_zz = buffer_spdd[47];

    auto g_0_y_xz_xx = buffer_spdd[48];

    auto g_0_y_xz_xy = buffer_spdd[49];

    auto g_0_y_xz_xz = buffer_spdd[50];

    auto g_0_y_xz_yy = buffer_spdd[51];

    auto g_0_y_xz_yz = buffer_spdd[52];

    auto g_0_y_xz_zz = buffer_spdd[53];

    auto g_0_y_yy_xx = buffer_spdd[54];

    auto g_0_y_yy_xy = buffer_spdd[55];

    auto g_0_y_yy_xz = buffer_spdd[56];

    auto g_0_y_yy_yy = buffer_spdd[57];

    auto g_0_y_yy_yz = buffer_spdd[58];

    auto g_0_y_yy_zz = buffer_spdd[59];

    auto g_0_y_yz_xx = buffer_spdd[60];

    auto g_0_y_yz_xy = buffer_spdd[61];

    auto g_0_y_yz_xz = buffer_spdd[62];

    auto g_0_y_yz_yy = buffer_spdd[63];

    auto g_0_y_yz_yz = buffer_spdd[64];

    auto g_0_y_yz_zz = buffer_spdd[65];

    auto g_0_y_zz_xx = buffer_spdd[66];

    auto g_0_y_zz_xy = buffer_spdd[67];

    auto g_0_y_zz_xz = buffer_spdd[68];

    auto g_0_y_zz_yy = buffer_spdd[69];

    auto g_0_y_zz_yz = buffer_spdd[70];

    auto g_0_y_zz_zz = buffer_spdd[71];

    auto g_0_z_xx_xx = buffer_spdd[72];

    auto g_0_z_xx_xy = buffer_spdd[73];

    auto g_0_z_xx_xz = buffer_spdd[74];

    auto g_0_z_xx_yy = buffer_spdd[75];

    auto g_0_z_xx_yz = buffer_spdd[76];

    auto g_0_z_xx_zz = buffer_spdd[77];

    auto g_0_z_xy_xx = buffer_spdd[78];

    auto g_0_z_xy_xy = buffer_spdd[79];

    auto g_0_z_xy_xz = buffer_spdd[80];

    auto g_0_z_xy_yy = buffer_spdd[81];

    auto g_0_z_xy_yz = buffer_spdd[82];

    auto g_0_z_xy_zz = buffer_spdd[83];

    auto g_0_z_xz_xx = buffer_spdd[84];

    auto g_0_z_xz_xy = buffer_spdd[85];

    auto g_0_z_xz_xz = buffer_spdd[86];

    auto g_0_z_xz_yy = buffer_spdd[87];

    auto g_0_z_xz_yz = buffer_spdd[88];

    auto g_0_z_xz_zz = buffer_spdd[89];

    auto g_0_z_yy_xx = buffer_spdd[90];

    auto g_0_z_yy_xy = buffer_spdd[91];

    auto g_0_z_yy_xz = buffer_spdd[92];

    auto g_0_z_yy_yy = buffer_spdd[93];

    auto g_0_z_yy_yz = buffer_spdd[94];

    auto g_0_z_yy_zz = buffer_spdd[95];

    auto g_0_z_yz_xx = buffer_spdd[96];

    auto g_0_z_yz_xy = buffer_spdd[97];

    auto g_0_z_yz_xz = buffer_spdd[98];

    auto g_0_z_yz_yy = buffer_spdd[99];

    auto g_0_z_yz_yz = buffer_spdd[100];

    auto g_0_z_yz_zz = buffer_spdd[101];

    auto g_0_z_zz_xx = buffer_spdd[102];

    auto g_0_z_zz_xy = buffer_spdd[103];

    auto g_0_z_zz_xz = buffer_spdd[104];

    auto g_0_z_zz_yy = buffer_spdd[105];

    auto g_0_z_zz_yz = buffer_spdd[106];

    auto g_0_z_zz_zz = buffer_spdd[107];

    /// Set up components of auxilary buffer : buffer_dpdd

    auto g_xx_x_xx_xx = buffer_dpdd[0];

    auto g_xx_x_xx_xy = buffer_dpdd[1];

    auto g_xx_x_xx_xz = buffer_dpdd[2];

    auto g_xx_x_xx_yy = buffer_dpdd[3];

    auto g_xx_x_xx_yz = buffer_dpdd[4];

    auto g_xx_x_xx_zz = buffer_dpdd[5];

    auto g_xx_x_xy_xx = buffer_dpdd[6];

    auto g_xx_x_xy_xy = buffer_dpdd[7];

    auto g_xx_x_xy_xz = buffer_dpdd[8];

    auto g_xx_x_xy_yy = buffer_dpdd[9];

    auto g_xx_x_xy_yz = buffer_dpdd[10];

    auto g_xx_x_xy_zz = buffer_dpdd[11];

    auto g_xx_x_xz_xx = buffer_dpdd[12];

    auto g_xx_x_xz_xy = buffer_dpdd[13];

    auto g_xx_x_xz_xz = buffer_dpdd[14];

    auto g_xx_x_xz_yy = buffer_dpdd[15];

    auto g_xx_x_xz_yz = buffer_dpdd[16];

    auto g_xx_x_xz_zz = buffer_dpdd[17];

    auto g_xx_x_yy_xx = buffer_dpdd[18];

    auto g_xx_x_yy_xy = buffer_dpdd[19];

    auto g_xx_x_yy_xz = buffer_dpdd[20];

    auto g_xx_x_yy_yy = buffer_dpdd[21];

    auto g_xx_x_yy_yz = buffer_dpdd[22];

    auto g_xx_x_yy_zz = buffer_dpdd[23];

    auto g_xx_x_yz_xx = buffer_dpdd[24];

    auto g_xx_x_yz_xy = buffer_dpdd[25];

    auto g_xx_x_yz_xz = buffer_dpdd[26];

    auto g_xx_x_yz_yy = buffer_dpdd[27];

    auto g_xx_x_yz_yz = buffer_dpdd[28];

    auto g_xx_x_yz_zz = buffer_dpdd[29];

    auto g_xx_x_zz_xx = buffer_dpdd[30];

    auto g_xx_x_zz_xy = buffer_dpdd[31];

    auto g_xx_x_zz_xz = buffer_dpdd[32];

    auto g_xx_x_zz_yy = buffer_dpdd[33];

    auto g_xx_x_zz_yz = buffer_dpdd[34];

    auto g_xx_x_zz_zz = buffer_dpdd[35];

    auto g_xx_y_xx_xx = buffer_dpdd[36];

    auto g_xx_y_xx_xy = buffer_dpdd[37];

    auto g_xx_y_xx_xz = buffer_dpdd[38];

    auto g_xx_y_xx_yy = buffer_dpdd[39];

    auto g_xx_y_xx_yz = buffer_dpdd[40];

    auto g_xx_y_xx_zz = buffer_dpdd[41];

    auto g_xx_y_xy_xx = buffer_dpdd[42];

    auto g_xx_y_xy_xy = buffer_dpdd[43];

    auto g_xx_y_xy_xz = buffer_dpdd[44];

    auto g_xx_y_xy_yy = buffer_dpdd[45];

    auto g_xx_y_xy_yz = buffer_dpdd[46];

    auto g_xx_y_xy_zz = buffer_dpdd[47];

    auto g_xx_y_xz_xx = buffer_dpdd[48];

    auto g_xx_y_xz_xy = buffer_dpdd[49];

    auto g_xx_y_xz_xz = buffer_dpdd[50];

    auto g_xx_y_xz_yy = buffer_dpdd[51];

    auto g_xx_y_xz_yz = buffer_dpdd[52];

    auto g_xx_y_xz_zz = buffer_dpdd[53];

    auto g_xx_y_yy_xx = buffer_dpdd[54];

    auto g_xx_y_yy_xy = buffer_dpdd[55];

    auto g_xx_y_yy_xz = buffer_dpdd[56];

    auto g_xx_y_yy_yy = buffer_dpdd[57];

    auto g_xx_y_yy_yz = buffer_dpdd[58];

    auto g_xx_y_yy_zz = buffer_dpdd[59];

    auto g_xx_y_yz_xx = buffer_dpdd[60];

    auto g_xx_y_yz_xy = buffer_dpdd[61];

    auto g_xx_y_yz_xz = buffer_dpdd[62];

    auto g_xx_y_yz_yy = buffer_dpdd[63];

    auto g_xx_y_yz_yz = buffer_dpdd[64];

    auto g_xx_y_yz_zz = buffer_dpdd[65];

    auto g_xx_y_zz_xx = buffer_dpdd[66];

    auto g_xx_y_zz_xy = buffer_dpdd[67];

    auto g_xx_y_zz_xz = buffer_dpdd[68];

    auto g_xx_y_zz_yy = buffer_dpdd[69];

    auto g_xx_y_zz_yz = buffer_dpdd[70];

    auto g_xx_y_zz_zz = buffer_dpdd[71];

    auto g_xx_z_xx_xx = buffer_dpdd[72];

    auto g_xx_z_xx_xy = buffer_dpdd[73];

    auto g_xx_z_xx_xz = buffer_dpdd[74];

    auto g_xx_z_xx_yy = buffer_dpdd[75];

    auto g_xx_z_xx_yz = buffer_dpdd[76];

    auto g_xx_z_xx_zz = buffer_dpdd[77];

    auto g_xx_z_xy_xx = buffer_dpdd[78];

    auto g_xx_z_xy_xy = buffer_dpdd[79];

    auto g_xx_z_xy_xz = buffer_dpdd[80];

    auto g_xx_z_xy_yy = buffer_dpdd[81];

    auto g_xx_z_xy_yz = buffer_dpdd[82];

    auto g_xx_z_xy_zz = buffer_dpdd[83];

    auto g_xx_z_xz_xx = buffer_dpdd[84];

    auto g_xx_z_xz_xy = buffer_dpdd[85];

    auto g_xx_z_xz_xz = buffer_dpdd[86];

    auto g_xx_z_xz_yy = buffer_dpdd[87];

    auto g_xx_z_xz_yz = buffer_dpdd[88];

    auto g_xx_z_xz_zz = buffer_dpdd[89];

    auto g_xx_z_yy_xx = buffer_dpdd[90];

    auto g_xx_z_yy_xy = buffer_dpdd[91];

    auto g_xx_z_yy_xz = buffer_dpdd[92];

    auto g_xx_z_yy_yy = buffer_dpdd[93];

    auto g_xx_z_yy_yz = buffer_dpdd[94];

    auto g_xx_z_yy_zz = buffer_dpdd[95];

    auto g_xx_z_yz_xx = buffer_dpdd[96];

    auto g_xx_z_yz_xy = buffer_dpdd[97];

    auto g_xx_z_yz_xz = buffer_dpdd[98];

    auto g_xx_z_yz_yy = buffer_dpdd[99];

    auto g_xx_z_yz_yz = buffer_dpdd[100];

    auto g_xx_z_yz_zz = buffer_dpdd[101];

    auto g_xx_z_zz_xx = buffer_dpdd[102];

    auto g_xx_z_zz_xy = buffer_dpdd[103];

    auto g_xx_z_zz_xz = buffer_dpdd[104];

    auto g_xx_z_zz_yy = buffer_dpdd[105];

    auto g_xx_z_zz_yz = buffer_dpdd[106];

    auto g_xx_z_zz_zz = buffer_dpdd[107];

    auto g_xy_x_xx_xx = buffer_dpdd[108];

    auto g_xy_x_xx_xy = buffer_dpdd[109];

    auto g_xy_x_xx_xz = buffer_dpdd[110];

    auto g_xy_x_xx_yy = buffer_dpdd[111];

    auto g_xy_x_xx_yz = buffer_dpdd[112];

    auto g_xy_x_xx_zz = buffer_dpdd[113];

    auto g_xy_x_xy_xx = buffer_dpdd[114];

    auto g_xy_x_xy_xy = buffer_dpdd[115];

    auto g_xy_x_xy_xz = buffer_dpdd[116];

    auto g_xy_x_xy_yy = buffer_dpdd[117];

    auto g_xy_x_xy_yz = buffer_dpdd[118];

    auto g_xy_x_xy_zz = buffer_dpdd[119];

    auto g_xy_x_xz_xx = buffer_dpdd[120];

    auto g_xy_x_xz_xy = buffer_dpdd[121];

    auto g_xy_x_xz_xz = buffer_dpdd[122];

    auto g_xy_x_xz_yy = buffer_dpdd[123];

    auto g_xy_x_xz_yz = buffer_dpdd[124];

    auto g_xy_x_xz_zz = buffer_dpdd[125];

    auto g_xy_x_yy_xx = buffer_dpdd[126];

    auto g_xy_x_yy_xy = buffer_dpdd[127];

    auto g_xy_x_yy_xz = buffer_dpdd[128];

    auto g_xy_x_yy_yy = buffer_dpdd[129];

    auto g_xy_x_yy_yz = buffer_dpdd[130];

    auto g_xy_x_yy_zz = buffer_dpdd[131];

    auto g_xy_x_yz_xx = buffer_dpdd[132];

    auto g_xy_x_yz_xy = buffer_dpdd[133];

    auto g_xy_x_yz_xz = buffer_dpdd[134];

    auto g_xy_x_yz_yy = buffer_dpdd[135];

    auto g_xy_x_yz_yz = buffer_dpdd[136];

    auto g_xy_x_yz_zz = buffer_dpdd[137];

    auto g_xy_x_zz_xx = buffer_dpdd[138];

    auto g_xy_x_zz_xy = buffer_dpdd[139];

    auto g_xy_x_zz_xz = buffer_dpdd[140];

    auto g_xy_x_zz_yy = buffer_dpdd[141];

    auto g_xy_x_zz_yz = buffer_dpdd[142];

    auto g_xy_x_zz_zz = buffer_dpdd[143];

    auto g_xy_y_xx_xx = buffer_dpdd[144];

    auto g_xy_y_xx_xy = buffer_dpdd[145];

    auto g_xy_y_xx_xz = buffer_dpdd[146];

    auto g_xy_y_xx_yy = buffer_dpdd[147];

    auto g_xy_y_xx_yz = buffer_dpdd[148];

    auto g_xy_y_xx_zz = buffer_dpdd[149];

    auto g_xy_y_xy_xx = buffer_dpdd[150];

    auto g_xy_y_xy_xy = buffer_dpdd[151];

    auto g_xy_y_xy_xz = buffer_dpdd[152];

    auto g_xy_y_xy_yy = buffer_dpdd[153];

    auto g_xy_y_xy_yz = buffer_dpdd[154];

    auto g_xy_y_xy_zz = buffer_dpdd[155];

    auto g_xy_y_xz_xx = buffer_dpdd[156];

    auto g_xy_y_xz_xy = buffer_dpdd[157];

    auto g_xy_y_xz_xz = buffer_dpdd[158];

    auto g_xy_y_xz_yy = buffer_dpdd[159];

    auto g_xy_y_xz_yz = buffer_dpdd[160];

    auto g_xy_y_xz_zz = buffer_dpdd[161];

    auto g_xy_y_yy_xx = buffer_dpdd[162];

    auto g_xy_y_yy_xy = buffer_dpdd[163];

    auto g_xy_y_yy_xz = buffer_dpdd[164];

    auto g_xy_y_yy_yy = buffer_dpdd[165];

    auto g_xy_y_yy_yz = buffer_dpdd[166];

    auto g_xy_y_yy_zz = buffer_dpdd[167];

    auto g_xy_y_yz_xx = buffer_dpdd[168];

    auto g_xy_y_yz_xy = buffer_dpdd[169];

    auto g_xy_y_yz_xz = buffer_dpdd[170];

    auto g_xy_y_yz_yy = buffer_dpdd[171];

    auto g_xy_y_yz_yz = buffer_dpdd[172];

    auto g_xy_y_yz_zz = buffer_dpdd[173];

    auto g_xy_y_zz_xx = buffer_dpdd[174];

    auto g_xy_y_zz_xy = buffer_dpdd[175];

    auto g_xy_y_zz_xz = buffer_dpdd[176];

    auto g_xy_y_zz_yy = buffer_dpdd[177];

    auto g_xy_y_zz_yz = buffer_dpdd[178];

    auto g_xy_y_zz_zz = buffer_dpdd[179];

    auto g_xy_z_xx_xx = buffer_dpdd[180];

    auto g_xy_z_xx_xy = buffer_dpdd[181];

    auto g_xy_z_xx_xz = buffer_dpdd[182];

    auto g_xy_z_xx_yy = buffer_dpdd[183];

    auto g_xy_z_xx_yz = buffer_dpdd[184];

    auto g_xy_z_xx_zz = buffer_dpdd[185];

    auto g_xy_z_xy_xx = buffer_dpdd[186];

    auto g_xy_z_xy_xy = buffer_dpdd[187];

    auto g_xy_z_xy_xz = buffer_dpdd[188];

    auto g_xy_z_xy_yy = buffer_dpdd[189];

    auto g_xy_z_xy_yz = buffer_dpdd[190];

    auto g_xy_z_xy_zz = buffer_dpdd[191];

    auto g_xy_z_xz_xx = buffer_dpdd[192];

    auto g_xy_z_xz_xy = buffer_dpdd[193];

    auto g_xy_z_xz_xz = buffer_dpdd[194];

    auto g_xy_z_xz_yy = buffer_dpdd[195];

    auto g_xy_z_xz_yz = buffer_dpdd[196];

    auto g_xy_z_xz_zz = buffer_dpdd[197];

    auto g_xy_z_yy_xx = buffer_dpdd[198];

    auto g_xy_z_yy_xy = buffer_dpdd[199];

    auto g_xy_z_yy_xz = buffer_dpdd[200];

    auto g_xy_z_yy_yy = buffer_dpdd[201];

    auto g_xy_z_yy_yz = buffer_dpdd[202];

    auto g_xy_z_yy_zz = buffer_dpdd[203];

    auto g_xy_z_yz_xx = buffer_dpdd[204];

    auto g_xy_z_yz_xy = buffer_dpdd[205];

    auto g_xy_z_yz_xz = buffer_dpdd[206];

    auto g_xy_z_yz_yy = buffer_dpdd[207];

    auto g_xy_z_yz_yz = buffer_dpdd[208];

    auto g_xy_z_yz_zz = buffer_dpdd[209];

    auto g_xy_z_zz_xx = buffer_dpdd[210];

    auto g_xy_z_zz_xy = buffer_dpdd[211];

    auto g_xy_z_zz_xz = buffer_dpdd[212];

    auto g_xy_z_zz_yy = buffer_dpdd[213];

    auto g_xy_z_zz_yz = buffer_dpdd[214];

    auto g_xy_z_zz_zz = buffer_dpdd[215];

    auto g_xz_x_xx_xx = buffer_dpdd[216];

    auto g_xz_x_xx_xy = buffer_dpdd[217];

    auto g_xz_x_xx_xz = buffer_dpdd[218];

    auto g_xz_x_xx_yy = buffer_dpdd[219];

    auto g_xz_x_xx_yz = buffer_dpdd[220];

    auto g_xz_x_xx_zz = buffer_dpdd[221];

    auto g_xz_x_xy_xx = buffer_dpdd[222];

    auto g_xz_x_xy_xy = buffer_dpdd[223];

    auto g_xz_x_xy_xz = buffer_dpdd[224];

    auto g_xz_x_xy_yy = buffer_dpdd[225];

    auto g_xz_x_xy_yz = buffer_dpdd[226];

    auto g_xz_x_xy_zz = buffer_dpdd[227];

    auto g_xz_x_xz_xx = buffer_dpdd[228];

    auto g_xz_x_xz_xy = buffer_dpdd[229];

    auto g_xz_x_xz_xz = buffer_dpdd[230];

    auto g_xz_x_xz_yy = buffer_dpdd[231];

    auto g_xz_x_xz_yz = buffer_dpdd[232];

    auto g_xz_x_xz_zz = buffer_dpdd[233];

    auto g_xz_x_yy_xx = buffer_dpdd[234];

    auto g_xz_x_yy_xy = buffer_dpdd[235];

    auto g_xz_x_yy_xz = buffer_dpdd[236];

    auto g_xz_x_yy_yy = buffer_dpdd[237];

    auto g_xz_x_yy_yz = buffer_dpdd[238];

    auto g_xz_x_yy_zz = buffer_dpdd[239];

    auto g_xz_x_yz_xx = buffer_dpdd[240];

    auto g_xz_x_yz_xy = buffer_dpdd[241];

    auto g_xz_x_yz_xz = buffer_dpdd[242];

    auto g_xz_x_yz_yy = buffer_dpdd[243];

    auto g_xz_x_yz_yz = buffer_dpdd[244];

    auto g_xz_x_yz_zz = buffer_dpdd[245];

    auto g_xz_x_zz_xx = buffer_dpdd[246];

    auto g_xz_x_zz_xy = buffer_dpdd[247];

    auto g_xz_x_zz_xz = buffer_dpdd[248];

    auto g_xz_x_zz_yy = buffer_dpdd[249];

    auto g_xz_x_zz_yz = buffer_dpdd[250];

    auto g_xz_x_zz_zz = buffer_dpdd[251];

    auto g_xz_y_xx_xx = buffer_dpdd[252];

    auto g_xz_y_xx_xy = buffer_dpdd[253];

    auto g_xz_y_xx_xz = buffer_dpdd[254];

    auto g_xz_y_xx_yy = buffer_dpdd[255];

    auto g_xz_y_xx_yz = buffer_dpdd[256];

    auto g_xz_y_xx_zz = buffer_dpdd[257];

    auto g_xz_y_xy_xx = buffer_dpdd[258];

    auto g_xz_y_xy_xy = buffer_dpdd[259];

    auto g_xz_y_xy_xz = buffer_dpdd[260];

    auto g_xz_y_xy_yy = buffer_dpdd[261];

    auto g_xz_y_xy_yz = buffer_dpdd[262];

    auto g_xz_y_xy_zz = buffer_dpdd[263];

    auto g_xz_y_xz_xx = buffer_dpdd[264];

    auto g_xz_y_xz_xy = buffer_dpdd[265];

    auto g_xz_y_xz_xz = buffer_dpdd[266];

    auto g_xz_y_xz_yy = buffer_dpdd[267];

    auto g_xz_y_xz_yz = buffer_dpdd[268];

    auto g_xz_y_xz_zz = buffer_dpdd[269];

    auto g_xz_y_yy_xx = buffer_dpdd[270];

    auto g_xz_y_yy_xy = buffer_dpdd[271];

    auto g_xz_y_yy_xz = buffer_dpdd[272];

    auto g_xz_y_yy_yy = buffer_dpdd[273];

    auto g_xz_y_yy_yz = buffer_dpdd[274];

    auto g_xz_y_yy_zz = buffer_dpdd[275];

    auto g_xz_y_yz_xx = buffer_dpdd[276];

    auto g_xz_y_yz_xy = buffer_dpdd[277];

    auto g_xz_y_yz_xz = buffer_dpdd[278];

    auto g_xz_y_yz_yy = buffer_dpdd[279];

    auto g_xz_y_yz_yz = buffer_dpdd[280];

    auto g_xz_y_yz_zz = buffer_dpdd[281];

    auto g_xz_y_zz_xx = buffer_dpdd[282];

    auto g_xz_y_zz_xy = buffer_dpdd[283];

    auto g_xz_y_zz_xz = buffer_dpdd[284];

    auto g_xz_y_zz_yy = buffer_dpdd[285];

    auto g_xz_y_zz_yz = buffer_dpdd[286];

    auto g_xz_y_zz_zz = buffer_dpdd[287];

    auto g_xz_z_xx_xx = buffer_dpdd[288];

    auto g_xz_z_xx_xy = buffer_dpdd[289];

    auto g_xz_z_xx_xz = buffer_dpdd[290];

    auto g_xz_z_xx_yy = buffer_dpdd[291];

    auto g_xz_z_xx_yz = buffer_dpdd[292];

    auto g_xz_z_xx_zz = buffer_dpdd[293];

    auto g_xz_z_xy_xx = buffer_dpdd[294];

    auto g_xz_z_xy_xy = buffer_dpdd[295];

    auto g_xz_z_xy_xz = buffer_dpdd[296];

    auto g_xz_z_xy_yy = buffer_dpdd[297];

    auto g_xz_z_xy_yz = buffer_dpdd[298];

    auto g_xz_z_xy_zz = buffer_dpdd[299];

    auto g_xz_z_xz_xx = buffer_dpdd[300];

    auto g_xz_z_xz_xy = buffer_dpdd[301];

    auto g_xz_z_xz_xz = buffer_dpdd[302];

    auto g_xz_z_xz_yy = buffer_dpdd[303];

    auto g_xz_z_xz_yz = buffer_dpdd[304];

    auto g_xz_z_xz_zz = buffer_dpdd[305];

    auto g_xz_z_yy_xx = buffer_dpdd[306];

    auto g_xz_z_yy_xy = buffer_dpdd[307];

    auto g_xz_z_yy_xz = buffer_dpdd[308];

    auto g_xz_z_yy_yy = buffer_dpdd[309];

    auto g_xz_z_yy_yz = buffer_dpdd[310];

    auto g_xz_z_yy_zz = buffer_dpdd[311];

    auto g_xz_z_yz_xx = buffer_dpdd[312];

    auto g_xz_z_yz_xy = buffer_dpdd[313];

    auto g_xz_z_yz_xz = buffer_dpdd[314];

    auto g_xz_z_yz_yy = buffer_dpdd[315];

    auto g_xz_z_yz_yz = buffer_dpdd[316];

    auto g_xz_z_yz_zz = buffer_dpdd[317];

    auto g_xz_z_zz_xx = buffer_dpdd[318];

    auto g_xz_z_zz_xy = buffer_dpdd[319];

    auto g_xz_z_zz_xz = buffer_dpdd[320];

    auto g_xz_z_zz_yy = buffer_dpdd[321];

    auto g_xz_z_zz_yz = buffer_dpdd[322];

    auto g_xz_z_zz_zz = buffer_dpdd[323];

    auto g_yy_x_xx_xx = buffer_dpdd[324];

    auto g_yy_x_xx_xy = buffer_dpdd[325];

    auto g_yy_x_xx_xz = buffer_dpdd[326];

    auto g_yy_x_xx_yy = buffer_dpdd[327];

    auto g_yy_x_xx_yz = buffer_dpdd[328];

    auto g_yy_x_xx_zz = buffer_dpdd[329];

    auto g_yy_x_xy_xx = buffer_dpdd[330];

    auto g_yy_x_xy_xy = buffer_dpdd[331];

    auto g_yy_x_xy_xz = buffer_dpdd[332];

    auto g_yy_x_xy_yy = buffer_dpdd[333];

    auto g_yy_x_xy_yz = buffer_dpdd[334];

    auto g_yy_x_xy_zz = buffer_dpdd[335];

    auto g_yy_x_xz_xx = buffer_dpdd[336];

    auto g_yy_x_xz_xy = buffer_dpdd[337];

    auto g_yy_x_xz_xz = buffer_dpdd[338];

    auto g_yy_x_xz_yy = buffer_dpdd[339];

    auto g_yy_x_xz_yz = buffer_dpdd[340];

    auto g_yy_x_xz_zz = buffer_dpdd[341];

    auto g_yy_x_yy_xx = buffer_dpdd[342];

    auto g_yy_x_yy_xy = buffer_dpdd[343];

    auto g_yy_x_yy_xz = buffer_dpdd[344];

    auto g_yy_x_yy_yy = buffer_dpdd[345];

    auto g_yy_x_yy_yz = buffer_dpdd[346];

    auto g_yy_x_yy_zz = buffer_dpdd[347];

    auto g_yy_x_yz_xx = buffer_dpdd[348];

    auto g_yy_x_yz_xy = buffer_dpdd[349];

    auto g_yy_x_yz_xz = buffer_dpdd[350];

    auto g_yy_x_yz_yy = buffer_dpdd[351];

    auto g_yy_x_yz_yz = buffer_dpdd[352];

    auto g_yy_x_yz_zz = buffer_dpdd[353];

    auto g_yy_x_zz_xx = buffer_dpdd[354];

    auto g_yy_x_zz_xy = buffer_dpdd[355];

    auto g_yy_x_zz_xz = buffer_dpdd[356];

    auto g_yy_x_zz_yy = buffer_dpdd[357];

    auto g_yy_x_zz_yz = buffer_dpdd[358];

    auto g_yy_x_zz_zz = buffer_dpdd[359];

    auto g_yy_y_xx_xx = buffer_dpdd[360];

    auto g_yy_y_xx_xy = buffer_dpdd[361];

    auto g_yy_y_xx_xz = buffer_dpdd[362];

    auto g_yy_y_xx_yy = buffer_dpdd[363];

    auto g_yy_y_xx_yz = buffer_dpdd[364];

    auto g_yy_y_xx_zz = buffer_dpdd[365];

    auto g_yy_y_xy_xx = buffer_dpdd[366];

    auto g_yy_y_xy_xy = buffer_dpdd[367];

    auto g_yy_y_xy_xz = buffer_dpdd[368];

    auto g_yy_y_xy_yy = buffer_dpdd[369];

    auto g_yy_y_xy_yz = buffer_dpdd[370];

    auto g_yy_y_xy_zz = buffer_dpdd[371];

    auto g_yy_y_xz_xx = buffer_dpdd[372];

    auto g_yy_y_xz_xy = buffer_dpdd[373];

    auto g_yy_y_xz_xz = buffer_dpdd[374];

    auto g_yy_y_xz_yy = buffer_dpdd[375];

    auto g_yy_y_xz_yz = buffer_dpdd[376];

    auto g_yy_y_xz_zz = buffer_dpdd[377];

    auto g_yy_y_yy_xx = buffer_dpdd[378];

    auto g_yy_y_yy_xy = buffer_dpdd[379];

    auto g_yy_y_yy_xz = buffer_dpdd[380];

    auto g_yy_y_yy_yy = buffer_dpdd[381];

    auto g_yy_y_yy_yz = buffer_dpdd[382];

    auto g_yy_y_yy_zz = buffer_dpdd[383];

    auto g_yy_y_yz_xx = buffer_dpdd[384];

    auto g_yy_y_yz_xy = buffer_dpdd[385];

    auto g_yy_y_yz_xz = buffer_dpdd[386];

    auto g_yy_y_yz_yy = buffer_dpdd[387];

    auto g_yy_y_yz_yz = buffer_dpdd[388];

    auto g_yy_y_yz_zz = buffer_dpdd[389];

    auto g_yy_y_zz_xx = buffer_dpdd[390];

    auto g_yy_y_zz_xy = buffer_dpdd[391];

    auto g_yy_y_zz_xz = buffer_dpdd[392];

    auto g_yy_y_zz_yy = buffer_dpdd[393];

    auto g_yy_y_zz_yz = buffer_dpdd[394];

    auto g_yy_y_zz_zz = buffer_dpdd[395];

    auto g_yy_z_xx_xx = buffer_dpdd[396];

    auto g_yy_z_xx_xy = buffer_dpdd[397];

    auto g_yy_z_xx_xz = buffer_dpdd[398];

    auto g_yy_z_xx_yy = buffer_dpdd[399];

    auto g_yy_z_xx_yz = buffer_dpdd[400];

    auto g_yy_z_xx_zz = buffer_dpdd[401];

    auto g_yy_z_xy_xx = buffer_dpdd[402];

    auto g_yy_z_xy_xy = buffer_dpdd[403];

    auto g_yy_z_xy_xz = buffer_dpdd[404];

    auto g_yy_z_xy_yy = buffer_dpdd[405];

    auto g_yy_z_xy_yz = buffer_dpdd[406];

    auto g_yy_z_xy_zz = buffer_dpdd[407];

    auto g_yy_z_xz_xx = buffer_dpdd[408];

    auto g_yy_z_xz_xy = buffer_dpdd[409];

    auto g_yy_z_xz_xz = buffer_dpdd[410];

    auto g_yy_z_xz_yy = buffer_dpdd[411];

    auto g_yy_z_xz_yz = buffer_dpdd[412];

    auto g_yy_z_xz_zz = buffer_dpdd[413];

    auto g_yy_z_yy_xx = buffer_dpdd[414];

    auto g_yy_z_yy_xy = buffer_dpdd[415];

    auto g_yy_z_yy_xz = buffer_dpdd[416];

    auto g_yy_z_yy_yy = buffer_dpdd[417];

    auto g_yy_z_yy_yz = buffer_dpdd[418];

    auto g_yy_z_yy_zz = buffer_dpdd[419];

    auto g_yy_z_yz_xx = buffer_dpdd[420];

    auto g_yy_z_yz_xy = buffer_dpdd[421];

    auto g_yy_z_yz_xz = buffer_dpdd[422];

    auto g_yy_z_yz_yy = buffer_dpdd[423];

    auto g_yy_z_yz_yz = buffer_dpdd[424];

    auto g_yy_z_yz_zz = buffer_dpdd[425];

    auto g_yy_z_zz_xx = buffer_dpdd[426];

    auto g_yy_z_zz_xy = buffer_dpdd[427];

    auto g_yy_z_zz_xz = buffer_dpdd[428];

    auto g_yy_z_zz_yy = buffer_dpdd[429];

    auto g_yy_z_zz_yz = buffer_dpdd[430];

    auto g_yy_z_zz_zz = buffer_dpdd[431];

    auto g_yz_x_xx_xx = buffer_dpdd[432];

    auto g_yz_x_xx_xy = buffer_dpdd[433];

    auto g_yz_x_xx_xz = buffer_dpdd[434];

    auto g_yz_x_xx_yy = buffer_dpdd[435];

    auto g_yz_x_xx_yz = buffer_dpdd[436];

    auto g_yz_x_xx_zz = buffer_dpdd[437];

    auto g_yz_x_xy_xx = buffer_dpdd[438];

    auto g_yz_x_xy_xy = buffer_dpdd[439];

    auto g_yz_x_xy_xz = buffer_dpdd[440];

    auto g_yz_x_xy_yy = buffer_dpdd[441];

    auto g_yz_x_xy_yz = buffer_dpdd[442];

    auto g_yz_x_xy_zz = buffer_dpdd[443];

    auto g_yz_x_xz_xx = buffer_dpdd[444];

    auto g_yz_x_xz_xy = buffer_dpdd[445];

    auto g_yz_x_xz_xz = buffer_dpdd[446];

    auto g_yz_x_xz_yy = buffer_dpdd[447];

    auto g_yz_x_xz_yz = buffer_dpdd[448];

    auto g_yz_x_xz_zz = buffer_dpdd[449];

    auto g_yz_x_yy_xx = buffer_dpdd[450];

    auto g_yz_x_yy_xy = buffer_dpdd[451];

    auto g_yz_x_yy_xz = buffer_dpdd[452];

    auto g_yz_x_yy_yy = buffer_dpdd[453];

    auto g_yz_x_yy_yz = buffer_dpdd[454];

    auto g_yz_x_yy_zz = buffer_dpdd[455];

    auto g_yz_x_yz_xx = buffer_dpdd[456];

    auto g_yz_x_yz_xy = buffer_dpdd[457];

    auto g_yz_x_yz_xz = buffer_dpdd[458];

    auto g_yz_x_yz_yy = buffer_dpdd[459];

    auto g_yz_x_yz_yz = buffer_dpdd[460];

    auto g_yz_x_yz_zz = buffer_dpdd[461];

    auto g_yz_x_zz_xx = buffer_dpdd[462];

    auto g_yz_x_zz_xy = buffer_dpdd[463];

    auto g_yz_x_zz_xz = buffer_dpdd[464];

    auto g_yz_x_zz_yy = buffer_dpdd[465];

    auto g_yz_x_zz_yz = buffer_dpdd[466];

    auto g_yz_x_zz_zz = buffer_dpdd[467];

    auto g_yz_y_xx_xx = buffer_dpdd[468];

    auto g_yz_y_xx_xy = buffer_dpdd[469];

    auto g_yz_y_xx_xz = buffer_dpdd[470];

    auto g_yz_y_xx_yy = buffer_dpdd[471];

    auto g_yz_y_xx_yz = buffer_dpdd[472];

    auto g_yz_y_xx_zz = buffer_dpdd[473];

    auto g_yz_y_xy_xx = buffer_dpdd[474];

    auto g_yz_y_xy_xy = buffer_dpdd[475];

    auto g_yz_y_xy_xz = buffer_dpdd[476];

    auto g_yz_y_xy_yy = buffer_dpdd[477];

    auto g_yz_y_xy_yz = buffer_dpdd[478];

    auto g_yz_y_xy_zz = buffer_dpdd[479];

    auto g_yz_y_xz_xx = buffer_dpdd[480];

    auto g_yz_y_xz_xy = buffer_dpdd[481];

    auto g_yz_y_xz_xz = buffer_dpdd[482];

    auto g_yz_y_xz_yy = buffer_dpdd[483];

    auto g_yz_y_xz_yz = buffer_dpdd[484];

    auto g_yz_y_xz_zz = buffer_dpdd[485];

    auto g_yz_y_yy_xx = buffer_dpdd[486];

    auto g_yz_y_yy_xy = buffer_dpdd[487];

    auto g_yz_y_yy_xz = buffer_dpdd[488];

    auto g_yz_y_yy_yy = buffer_dpdd[489];

    auto g_yz_y_yy_yz = buffer_dpdd[490];

    auto g_yz_y_yy_zz = buffer_dpdd[491];

    auto g_yz_y_yz_xx = buffer_dpdd[492];

    auto g_yz_y_yz_xy = buffer_dpdd[493];

    auto g_yz_y_yz_xz = buffer_dpdd[494];

    auto g_yz_y_yz_yy = buffer_dpdd[495];

    auto g_yz_y_yz_yz = buffer_dpdd[496];

    auto g_yz_y_yz_zz = buffer_dpdd[497];

    auto g_yz_y_zz_xx = buffer_dpdd[498];

    auto g_yz_y_zz_xy = buffer_dpdd[499];

    auto g_yz_y_zz_xz = buffer_dpdd[500];

    auto g_yz_y_zz_yy = buffer_dpdd[501];

    auto g_yz_y_zz_yz = buffer_dpdd[502];

    auto g_yz_y_zz_zz = buffer_dpdd[503];

    auto g_yz_z_xx_xx = buffer_dpdd[504];

    auto g_yz_z_xx_xy = buffer_dpdd[505];

    auto g_yz_z_xx_xz = buffer_dpdd[506];

    auto g_yz_z_xx_yy = buffer_dpdd[507];

    auto g_yz_z_xx_yz = buffer_dpdd[508];

    auto g_yz_z_xx_zz = buffer_dpdd[509];

    auto g_yz_z_xy_xx = buffer_dpdd[510];

    auto g_yz_z_xy_xy = buffer_dpdd[511];

    auto g_yz_z_xy_xz = buffer_dpdd[512];

    auto g_yz_z_xy_yy = buffer_dpdd[513];

    auto g_yz_z_xy_yz = buffer_dpdd[514];

    auto g_yz_z_xy_zz = buffer_dpdd[515];

    auto g_yz_z_xz_xx = buffer_dpdd[516];

    auto g_yz_z_xz_xy = buffer_dpdd[517];

    auto g_yz_z_xz_xz = buffer_dpdd[518];

    auto g_yz_z_xz_yy = buffer_dpdd[519];

    auto g_yz_z_xz_yz = buffer_dpdd[520];

    auto g_yz_z_xz_zz = buffer_dpdd[521];

    auto g_yz_z_yy_xx = buffer_dpdd[522];

    auto g_yz_z_yy_xy = buffer_dpdd[523];

    auto g_yz_z_yy_xz = buffer_dpdd[524];

    auto g_yz_z_yy_yy = buffer_dpdd[525];

    auto g_yz_z_yy_yz = buffer_dpdd[526];

    auto g_yz_z_yy_zz = buffer_dpdd[527];

    auto g_yz_z_yz_xx = buffer_dpdd[528];

    auto g_yz_z_yz_xy = buffer_dpdd[529];

    auto g_yz_z_yz_xz = buffer_dpdd[530];

    auto g_yz_z_yz_yy = buffer_dpdd[531];

    auto g_yz_z_yz_yz = buffer_dpdd[532];

    auto g_yz_z_yz_zz = buffer_dpdd[533];

    auto g_yz_z_zz_xx = buffer_dpdd[534];

    auto g_yz_z_zz_xy = buffer_dpdd[535];

    auto g_yz_z_zz_xz = buffer_dpdd[536];

    auto g_yz_z_zz_yy = buffer_dpdd[537];

    auto g_yz_z_zz_yz = buffer_dpdd[538];

    auto g_yz_z_zz_zz = buffer_dpdd[539];

    auto g_zz_x_xx_xx = buffer_dpdd[540];

    auto g_zz_x_xx_xy = buffer_dpdd[541];

    auto g_zz_x_xx_xz = buffer_dpdd[542];

    auto g_zz_x_xx_yy = buffer_dpdd[543];

    auto g_zz_x_xx_yz = buffer_dpdd[544];

    auto g_zz_x_xx_zz = buffer_dpdd[545];

    auto g_zz_x_xy_xx = buffer_dpdd[546];

    auto g_zz_x_xy_xy = buffer_dpdd[547];

    auto g_zz_x_xy_xz = buffer_dpdd[548];

    auto g_zz_x_xy_yy = buffer_dpdd[549];

    auto g_zz_x_xy_yz = buffer_dpdd[550];

    auto g_zz_x_xy_zz = buffer_dpdd[551];

    auto g_zz_x_xz_xx = buffer_dpdd[552];

    auto g_zz_x_xz_xy = buffer_dpdd[553];

    auto g_zz_x_xz_xz = buffer_dpdd[554];

    auto g_zz_x_xz_yy = buffer_dpdd[555];

    auto g_zz_x_xz_yz = buffer_dpdd[556];

    auto g_zz_x_xz_zz = buffer_dpdd[557];

    auto g_zz_x_yy_xx = buffer_dpdd[558];

    auto g_zz_x_yy_xy = buffer_dpdd[559];

    auto g_zz_x_yy_xz = buffer_dpdd[560];

    auto g_zz_x_yy_yy = buffer_dpdd[561];

    auto g_zz_x_yy_yz = buffer_dpdd[562];

    auto g_zz_x_yy_zz = buffer_dpdd[563];

    auto g_zz_x_yz_xx = buffer_dpdd[564];

    auto g_zz_x_yz_xy = buffer_dpdd[565];

    auto g_zz_x_yz_xz = buffer_dpdd[566];

    auto g_zz_x_yz_yy = buffer_dpdd[567];

    auto g_zz_x_yz_yz = buffer_dpdd[568];

    auto g_zz_x_yz_zz = buffer_dpdd[569];

    auto g_zz_x_zz_xx = buffer_dpdd[570];

    auto g_zz_x_zz_xy = buffer_dpdd[571];

    auto g_zz_x_zz_xz = buffer_dpdd[572];

    auto g_zz_x_zz_yy = buffer_dpdd[573];

    auto g_zz_x_zz_yz = buffer_dpdd[574];

    auto g_zz_x_zz_zz = buffer_dpdd[575];

    auto g_zz_y_xx_xx = buffer_dpdd[576];

    auto g_zz_y_xx_xy = buffer_dpdd[577];

    auto g_zz_y_xx_xz = buffer_dpdd[578];

    auto g_zz_y_xx_yy = buffer_dpdd[579];

    auto g_zz_y_xx_yz = buffer_dpdd[580];

    auto g_zz_y_xx_zz = buffer_dpdd[581];

    auto g_zz_y_xy_xx = buffer_dpdd[582];

    auto g_zz_y_xy_xy = buffer_dpdd[583];

    auto g_zz_y_xy_xz = buffer_dpdd[584];

    auto g_zz_y_xy_yy = buffer_dpdd[585];

    auto g_zz_y_xy_yz = buffer_dpdd[586];

    auto g_zz_y_xy_zz = buffer_dpdd[587];

    auto g_zz_y_xz_xx = buffer_dpdd[588];

    auto g_zz_y_xz_xy = buffer_dpdd[589];

    auto g_zz_y_xz_xz = buffer_dpdd[590];

    auto g_zz_y_xz_yy = buffer_dpdd[591];

    auto g_zz_y_xz_yz = buffer_dpdd[592];

    auto g_zz_y_xz_zz = buffer_dpdd[593];

    auto g_zz_y_yy_xx = buffer_dpdd[594];

    auto g_zz_y_yy_xy = buffer_dpdd[595];

    auto g_zz_y_yy_xz = buffer_dpdd[596];

    auto g_zz_y_yy_yy = buffer_dpdd[597];

    auto g_zz_y_yy_yz = buffer_dpdd[598];

    auto g_zz_y_yy_zz = buffer_dpdd[599];

    auto g_zz_y_yz_xx = buffer_dpdd[600];

    auto g_zz_y_yz_xy = buffer_dpdd[601];

    auto g_zz_y_yz_xz = buffer_dpdd[602];

    auto g_zz_y_yz_yy = buffer_dpdd[603];

    auto g_zz_y_yz_yz = buffer_dpdd[604];

    auto g_zz_y_yz_zz = buffer_dpdd[605];

    auto g_zz_y_zz_xx = buffer_dpdd[606];

    auto g_zz_y_zz_xy = buffer_dpdd[607];

    auto g_zz_y_zz_xz = buffer_dpdd[608];

    auto g_zz_y_zz_yy = buffer_dpdd[609];

    auto g_zz_y_zz_yz = buffer_dpdd[610];

    auto g_zz_y_zz_zz = buffer_dpdd[611];

    auto g_zz_z_xx_xx = buffer_dpdd[612];

    auto g_zz_z_xx_xy = buffer_dpdd[613];

    auto g_zz_z_xx_xz = buffer_dpdd[614];

    auto g_zz_z_xx_yy = buffer_dpdd[615];

    auto g_zz_z_xx_yz = buffer_dpdd[616];

    auto g_zz_z_xx_zz = buffer_dpdd[617];

    auto g_zz_z_xy_xx = buffer_dpdd[618];

    auto g_zz_z_xy_xy = buffer_dpdd[619];

    auto g_zz_z_xy_xz = buffer_dpdd[620];

    auto g_zz_z_xy_yy = buffer_dpdd[621];

    auto g_zz_z_xy_yz = buffer_dpdd[622];

    auto g_zz_z_xy_zz = buffer_dpdd[623];

    auto g_zz_z_xz_xx = buffer_dpdd[624];

    auto g_zz_z_xz_xy = buffer_dpdd[625];

    auto g_zz_z_xz_xz = buffer_dpdd[626];

    auto g_zz_z_xz_yy = buffer_dpdd[627];

    auto g_zz_z_xz_yz = buffer_dpdd[628];

    auto g_zz_z_xz_zz = buffer_dpdd[629];

    auto g_zz_z_yy_xx = buffer_dpdd[630];

    auto g_zz_z_yy_xy = buffer_dpdd[631];

    auto g_zz_z_yy_xz = buffer_dpdd[632];

    auto g_zz_z_yy_yy = buffer_dpdd[633];

    auto g_zz_z_yy_yz = buffer_dpdd[634];

    auto g_zz_z_yy_zz = buffer_dpdd[635];

    auto g_zz_z_yz_xx = buffer_dpdd[636];

    auto g_zz_z_yz_xy = buffer_dpdd[637];

    auto g_zz_z_yz_xz = buffer_dpdd[638];

    auto g_zz_z_yz_yy = buffer_dpdd[639];

    auto g_zz_z_yz_yz = buffer_dpdd[640];

    auto g_zz_z_yz_zz = buffer_dpdd[641];

    auto g_zz_z_zz_xx = buffer_dpdd[642];

    auto g_zz_z_zz_xy = buffer_dpdd[643];

    auto g_zz_z_zz_xz = buffer_dpdd[644];

    auto g_zz_z_zz_yy = buffer_dpdd[645];

    auto g_zz_z_zz_yz = buffer_dpdd[646];

    auto g_zz_z_zz_zz = buffer_dpdd[647];

    /// Set up components of integrals buffer : buffer_2000_spdd

    auto g_xx_0_0_0_0_x_xx_xx = buffer_2000_spdd[0];

    auto g_xx_0_0_0_0_x_xx_xy = buffer_2000_spdd[1];

    auto g_xx_0_0_0_0_x_xx_xz = buffer_2000_spdd[2];

    auto g_xx_0_0_0_0_x_xx_yy = buffer_2000_spdd[3];

    auto g_xx_0_0_0_0_x_xx_yz = buffer_2000_spdd[4];

    auto g_xx_0_0_0_0_x_xx_zz = buffer_2000_spdd[5];

    auto g_xx_0_0_0_0_x_xy_xx = buffer_2000_spdd[6];

    auto g_xx_0_0_0_0_x_xy_xy = buffer_2000_spdd[7];

    auto g_xx_0_0_0_0_x_xy_xz = buffer_2000_spdd[8];

    auto g_xx_0_0_0_0_x_xy_yy = buffer_2000_spdd[9];

    auto g_xx_0_0_0_0_x_xy_yz = buffer_2000_spdd[10];

    auto g_xx_0_0_0_0_x_xy_zz = buffer_2000_spdd[11];

    auto g_xx_0_0_0_0_x_xz_xx = buffer_2000_spdd[12];

    auto g_xx_0_0_0_0_x_xz_xy = buffer_2000_spdd[13];

    auto g_xx_0_0_0_0_x_xz_xz = buffer_2000_spdd[14];

    auto g_xx_0_0_0_0_x_xz_yy = buffer_2000_spdd[15];

    auto g_xx_0_0_0_0_x_xz_yz = buffer_2000_spdd[16];

    auto g_xx_0_0_0_0_x_xz_zz = buffer_2000_spdd[17];

    auto g_xx_0_0_0_0_x_yy_xx = buffer_2000_spdd[18];

    auto g_xx_0_0_0_0_x_yy_xy = buffer_2000_spdd[19];

    auto g_xx_0_0_0_0_x_yy_xz = buffer_2000_spdd[20];

    auto g_xx_0_0_0_0_x_yy_yy = buffer_2000_spdd[21];

    auto g_xx_0_0_0_0_x_yy_yz = buffer_2000_spdd[22];

    auto g_xx_0_0_0_0_x_yy_zz = buffer_2000_spdd[23];

    auto g_xx_0_0_0_0_x_yz_xx = buffer_2000_spdd[24];

    auto g_xx_0_0_0_0_x_yz_xy = buffer_2000_spdd[25];

    auto g_xx_0_0_0_0_x_yz_xz = buffer_2000_spdd[26];

    auto g_xx_0_0_0_0_x_yz_yy = buffer_2000_spdd[27];

    auto g_xx_0_0_0_0_x_yz_yz = buffer_2000_spdd[28];

    auto g_xx_0_0_0_0_x_yz_zz = buffer_2000_spdd[29];

    auto g_xx_0_0_0_0_x_zz_xx = buffer_2000_spdd[30];

    auto g_xx_0_0_0_0_x_zz_xy = buffer_2000_spdd[31];

    auto g_xx_0_0_0_0_x_zz_xz = buffer_2000_spdd[32];

    auto g_xx_0_0_0_0_x_zz_yy = buffer_2000_spdd[33];

    auto g_xx_0_0_0_0_x_zz_yz = buffer_2000_spdd[34];

    auto g_xx_0_0_0_0_x_zz_zz = buffer_2000_spdd[35];

    auto g_xx_0_0_0_0_y_xx_xx = buffer_2000_spdd[36];

    auto g_xx_0_0_0_0_y_xx_xy = buffer_2000_spdd[37];

    auto g_xx_0_0_0_0_y_xx_xz = buffer_2000_spdd[38];

    auto g_xx_0_0_0_0_y_xx_yy = buffer_2000_spdd[39];

    auto g_xx_0_0_0_0_y_xx_yz = buffer_2000_spdd[40];

    auto g_xx_0_0_0_0_y_xx_zz = buffer_2000_spdd[41];

    auto g_xx_0_0_0_0_y_xy_xx = buffer_2000_spdd[42];

    auto g_xx_0_0_0_0_y_xy_xy = buffer_2000_spdd[43];

    auto g_xx_0_0_0_0_y_xy_xz = buffer_2000_spdd[44];

    auto g_xx_0_0_0_0_y_xy_yy = buffer_2000_spdd[45];

    auto g_xx_0_0_0_0_y_xy_yz = buffer_2000_spdd[46];

    auto g_xx_0_0_0_0_y_xy_zz = buffer_2000_spdd[47];

    auto g_xx_0_0_0_0_y_xz_xx = buffer_2000_spdd[48];

    auto g_xx_0_0_0_0_y_xz_xy = buffer_2000_spdd[49];

    auto g_xx_0_0_0_0_y_xz_xz = buffer_2000_spdd[50];

    auto g_xx_0_0_0_0_y_xz_yy = buffer_2000_spdd[51];

    auto g_xx_0_0_0_0_y_xz_yz = buffer_2000_spdd[52];

    auto g_xx_0_0_0_0_y_xz_zz = buffer_2000_spdd[53];

    auto g_xx_0_0_0_0_y_yy_xx = buffer_2000_spdd[54];

    auto g_xx_0_0_0_0_y_yy_xy = buffer_2000_spdd[55];

    auto g_xx_0_0_0_0_y_yy_xz = buffer_2000_spdd[56];

    auto g_xx_0_0_0_0_y_yy_yy = buffer_2000_spdd[57];

    auto g_xx_0_0_0_0_y_yy_yz = buffer_2000_spdd[58];

    auto g_xx_0_0_0_0_y_yy_zz = buffer_2000_spdd[59];

    auto g_xx_0_0_0_0_y_yz_xx = buffer_2000_spdd[60];

    auto g_xx_0_0_0_0_y_yz_xy = buffer_2000_spdd[61];

    auto g_xx_0_0_0_0_y_yz_xz = buffer_2000_spdd[62];

    auto g_xx_0_0_0_0_y_yz_yy = buffer_2000_spdd[63];

    auto g_xx_0_0_0_0_y_yz_yz = buffer_2000_spdd[64];

    auto g_xx_0_0_0_0_y_yz_zz = buffer_2000_spdd[65];

    auto g_xx_0_0_0_0_y_zz_xx = buffer_2000_spdd[66];

    auto g_xx_0_0_0_0_y_zz_xy = buffer_2000_spdd[67];

    auto g_xx_0_0_0_0_y_zz_xz = buffer_2000_spdd[68];

    auto g_xx_0_0_0_0_y_zz_yy = buffer_2000_spdd[69];

    auto g_xx_0_0_0_0_y_zz_yz = buffer_2000_spdd[70];

    auto g_xx_0_0_0_0_y_zz_zz = buffer_2000_spdd[71];

    auto g_xx_0_0_0_0_z_xx_xx = buffer_2000_spdd[72];

    auto g_xx_0_0_0_0_z_xx_xy = buffer_2000_spdd[73];

    auto g_xx_0_0_0_0_z_xx_xz = buffer_2000_spdd[74];

    auto g_xx_0_0_0_0_z_xx_yy = buffer_2000_spdd[75];

    auto g_xx_0_0_0_0_z_xx_yz = buffer_2000_spdd[76];

    auto g_xx_0_0_0_0_z_xx_zz = buffer_2000_spdd[77];

    auto g_xx_0_0_0_0_z_xy_xx = buffer_2000_spdd[78];

    auto g_xx_0_0_0_0_z_xy_xy = buffer_2000_spdd[79];

    auto g_xx_0_0_0_0_z_xy_xz = buffer_2000_spdd[80];

    auto g_xx_0_0_0_0_z_xy_yy = buffer_2000_spdd[81];

    auto g_xx_0_0_0_0_z_xy_yz = buffer_2000_spdd[82];

    auto g_xx_0_0_0_0_z_xy_zz = buffer_2000_spdd[83];

    auto g_xx_0_0_0_0_z_xz_xx = buffer_2000_spdd[84];

    auto g_xx_0_0_0_0_z_xz_xy = buffer_2000_spdd[85];

    auto g_xx_0_0_0_0_z_xz_xz = buffer_2000_spdd[86];

    auto g_xx_0_0_0_0_z_xz_yy = buffer_2000_spdd[87];

    auto g_xx_0_0_0_0_z_xz_yz = buffer_2000_spdd[88];

    auto g_xx_0_0_0_0_z_xz_zz = buffer_2000_spdd[89];

    auto g_xx_0_0_0_0_z_yy_xx = buffer_2000_spdd[90];

    auto g_xx_0_0_0_0_z_yy_xy = buffer_2000_spdd[91];

    auto g_xx_0_0_0_0_z_yy_xz = buffer_2000_spdd[92];

    auto g_xx_0_0_0_0_z_yy_yy = buffer_2000_spdd[93];

    auto g_xx_0_0_0_0_z_yy_yz = buffer_2000_spdd[94];

    auto g_xx_0_0_0_0_z_yy_zz = buffer_2000_spdd[95];

    auto g_xx_0_0_0_0_z_yz_xx = buffer_2000_spdd[96];

    auto g_xx_0_0_0_0_z_yz_xy = buffer_2000_spdd[97];

    auto g_xx_0_0_0_0_z_yz_xz = buffer_2000_spdd[98];

    auto g_xx_0_0_0_0_z_yz_yy = buffer_2000_spdd[99];

    auto g_xx_0_0_0_0_z_yz_yz = buffer_2000_spdd[100];

    auto g_xx_0_0_0_0_z_yz_zz = buffer_2000_spdd[101];

    auto g_xx_0_0_0_0_z_zz_xx = buffer_2000_spdd[102];

    auto g_xx_0_0_0_0_z_zz_xy = buffer_2000_spdd[103];

    auto g_xx_0_0_0_0_z_zz_xz = buffer_2000_spdd[104];

    auto g_xx_0_0_0_0_z_zz_yy = buffer_2000_spdd[105];

    auto g_xx_0_0_0_0_z_zz_yz = buffer_2000_spdd[106];

    auto g_xx_0_0_0_0_z_zz_zz = buffer_2000_spdd[107];

    auto g_xy_0_0_0_0_x_xx_xx = buffer_2000_spdd[108];

    auto g_xy_0_0_0_0_x_xx_xy = buffer_2000_spdd[109];

    auto g_xy_0_0_0_0_x_xx_xz = buffer_2000_spdd[110];

    auto g_xy_0_0_0_0_x_xx_yy = buffer_2000_spdd[111];

    auto g_xy_0_0_0_0_x_xx_yz = buffer_2000_spdd[112];

    auto g_xy_0_0_0_0_x_xx_zz = buffer_2000_spdd[113];

    auto g_xy_0_0_0_0_x_xy_xx = buffer_2000_spdd[114];

    auto g_xy_0_0_0_0_x_xy_xy = buffer_2000_spdd[115];

    auto g_xy_0_0_0_0_x_xy_xz = buffer_2000_spdd[116];

    auto g_xy_0_0_0_0_x_xy_yy = buffer_2000_spdd[117];

    auto g_xy_0_0_0_0_x_xy_yz = buffer_2000_spdd[118];

    auto g_xy_0_0_0_0_x_xy_zz = buffer_2000_spdd[119];

    auto g_xy_0_0_0_0_x_xz_xx = buffer_2000_spdd[120];

    auto g_xy_0_0_0_0_x_xz_xy = buffer_2000_spdd[121];

    auto g_xy_0_0_0_0_x_xz_xz = buffer_2000_spdd[122];

    auto g_xy_0_0_0_0_x_xz_yy = buffer_2000_spdd[123];

    auto g_xy_0_0_0_0_x_xz_yz = buffer_2000_spdd[124];

    auto g_xy_0_0_0_0_x_xz_zz = buffer_2000_spdd[125];

    auto g_xy_0_0_0_0_x_yy_xx = buffer_2000_spdd[126];

    auto g_xy_0_0_0_0_x_yy_xy = buffer_2000_spdd[127];

    auto g_xy_0_0_0_0_x_yy_xz = buffer_2000_spdd[128];

    auto g_xy_0_0_0_0_x_yy_yy = buffer_2000_spdd[129];

    auto g_xy_0_0_0_0_x_yy_yz = buffer_2000_spdd[130];

    auto g_xy_0_0_0_0_x_yy_zz = buffer_2000_spdd[131];

    auto g_xy_0_0_0_0_x_yz_xx = buffer_2000_spdd[132];

    auto g_xy_0_0_0_0_x_yz_xy = buffer_2000_spdd[133];

    auto g_xy_0_0_0_0_x_yz_xz = buffer_2000_spdd[134];

    auto g_xy_0_0_0_0_x_yz_yy = buffer_2000_spdd[135];

    auto g_xy_0_0_0_0_x_yz_yz = buffer_2000_spdd[136];

    auto g_xy_0_0_0_0_x_yz_zz = buffer_2000_spdd[137];

    auto g_xy_0_0_0_0_x_zz_xx = buffer_2000_spdd[138];

    auto g_xy_0_0_0_0_x_zz_xy = buffer_2000_spdd[139];

    auto g_xy_0_0_0_0_x_zz_xz = buffer_2000_spdd[140];

    auto g_xy_0_0_0_0_x_zz_yy = buffer_2000_spdd[141];

    auto g_xy_0_0_0_0_x_zz_yz = buffer_2000_spdd[142];

    auto g_xy_0_0_0_0_x_zz_zz = buffer_2000_spdd[143];

    auto g_xy_0_0_0_0_y_xx_xx = buffer_2000_spdd[144];

    auto g_xy_0_0_0_0_y_xx_xy = buffer_2000_spdd[145];

    auto g_xy_0_0_0_0_y_xx_xz = buffer_2000_spdd[146];

    auto g_xy_0_0_0_0_y_xx_yy = buffer_2000_spdd[147];

    auto g_xy_0_0_0_0_y_xx_yz = buffer_2000_spdd[148];

    auto g_xy_0_0_0_0_y_xx_zz = buffer_2000_spdd[149];

    auto g_xy_0_0_0_0_y_xy_xx = buffer_2000_spdd[150];

    auto g_xy_0_0_0_0_y_xy_xy = buffer_2000_spdd[151];

    auto g_xy_0_0_0_0_y_xy_xz = buffer_2000_spdd[152];

    auto g_xy_0_0_0_0_y_xy_yy = buffer_2000_spdd[153];

    auto g_xy_0_0_0_0_y_xy_yz = buffer_2000_spdd[154];

    auto g_xy_0_0_0_0_y_xy_zz = buffer_2000_spdd[155];

    auto g_xy_0_0_0_0_y_xz_xx = buffer_2000_spdd[156];

    auto g_xy_0_0_0_0_y_xz_xy = buffer_2000_spdd[157];

    auto g_xy_0_0_0_0_y_xz_xz = buffer_2000_spdd[158];

    auto g_xy_0_0_0_0_y_xz_yy = buffer_2000_spdd[159];

    auto g_xy_0_0_0_0_y_xz_yz = buffer_2000_spdd[160];

    auto g_xy_0_0_0_0_y_xz_zz = buffer_2000_spdd[161];

    auto g_xy_0_0_0_0_y_yy_xx = buffer_2000_spdd[162];

    auto g_xy_0_0_0_0_y_yy_xy = buffer_2000_spdd[163];

    auto g_xy_0_0_0_0_y_yy_xz = buffer_2000_spdd[164];

    auto g_xy_0_0_0_0_y_yy_yy = buffer_2000_spdd[165];

    auto g_xy_0_0_0_0_y_yy_yz = buffer_2000_spdd[166];

    auto g_xy_0_0_0_0_y_yy_zz = buffer_2000_spdd[167];

    auto g_xy_0_0_0_0_y_yz_xx = buffer_2000_spdd[168];

    auto g_xy_0_0_0_0_y_yz_xy = buffer_2000_spdd[169];

    auto g_xy_0_0_0_0_y_yz_xz = buffer_2000_spdd[170];

    auto g_xy_0_0_0_0_y_yz_yy = buffer_2000_spdd[171];

    auto g_xy_0_0_0_0_y_yz_yz = buffer_2000_spdd[172];

    auto g_xy_0_0_0_0_y_yz_zz = buffer_2000_spdd[173];

    auto g_xy_0_0_0_0_y_zz_xx = buffer_2000_spdd[174];

    auto g_xy_0_0_0_0_y_zz_xy = buffer_2000_spdd[175];

    auto g_xy_0_0_0_0_y_zz_xz = buffer_2000_spdd[176];

    auto g_xy_0_0_0_0_y_zz_yy = buffer_2000_spdd[177];

    auto g_xy_0_0_0_0_y_zz_yz = buffer_2000_spdd[178];

    auto g_xy_0_0_0_0_y_zz_zz = buffer_2000_spdd[179];

    auto g_xy_0_0_0_0_z_xx_xx = buffer_2000_spdd[180];

    auto g_xy_0_0_0_0_z_xx_xy = buffer_2000_spdd[181];

    auto g_xy_0_0_0_0_z_xx_xz = buffer_2000_spdd[182];

    auto g_xy_0_0_0_0_z_xx_yy = buffer_2000_spdd[183];

    auto g_xy_0_0_0_0_z_xx_yz = buffer_2000_spdd[184];

    auto g_xy_0_0_0_0_z_xx_zz = buffer_2000_spdd[185];

    auto g_xy_0_0_0_0_z_xy_xx = buffer_2000_spdd[186];

    auto g_xy_0_0_0_0_z_xy_xy = buffer_2000_spdd[187];

    auto g_xy_0_0_0_0_z_xy_xz = buffer_2000_spdd[188];

    auto g_xy_0_0_0_0_z_xy_yy = buffer_2000_spdd[189];

    auto g_xy_0_0_0_0_z_xy_yz = buffer_2000_spdd[190];

    auto g_xy_0_0_0_0_z_xy_zz = buffer_2000_spdd[191];

    auto g_xy_0_0_0_0_z_xz_xx = buffer_2000_spdd[192];

    auto g_xy_0_0_0_0_z_xz_xy = buffer_2000_spdd[193];

    auto g_xy_0_0_0_0_z_xz_xz = buffer_2000_spdd[194];

    auto g_xy_0_0_0_0_z_xz_yy = buffer_2000_spdd[195];

    auto g_xy_0_0_0_0_z_xz_yz = buffer_2000_spdd[196];

    auto g_xy_0_0_0_0_z_xz_zz = buffer_2000_spdd[197];

    auto g_xy_0_0_0_0_z_yy_xx = buffer_2000_spdd[198];

    auto g_xy_0_0_0_0_z_yy_xy = buffer_2000_spdd[199];

    auto g_xy_0_0_0_0_z_yy_xz = buffer_2000_spdd[200];

    auto g_xy_0_0_0_0_z_yy_yy = buffer_2000_spdd[201];

    auto g_xy_0_0_0_0_z_yy_yz = buffer_2000_spdd[202];

    auto g_xy_0_0_0_0_z_yy_zz = buffer_2000_spdd[203];

    auto g_xy_0_0_0_0_z_yz_xx = buffer_2000_spdd[204];

    auto g_xy_0_0_0_0_z_yz_xy = buffer_2000_spdd[205];

    auto g_xy_0_0_0_0_z_yz_xz = buffer_2000_spdd[206];

    auto g_xy_0_0_0_0_z_yz_yy = buffer_2000_spdd[207];

    auto g_xy_0_0_0_0_z_yz_yz = buffer_2000_spdd[208];

    auto g_xy_0_0_0_0_z_yz_zz = buffer_2000_spdd[209];

    auto g_xy_0_0_0_0_z_zz_xx = buffer_2000_spdd[210];

    auto g_xy_0_0_0_0_z_zz_xy = buffer_2000_spdd[211];

    auto g_xy_0_0_0_0_z_zz_xz = buffer_2000_spdd[212];

    auto g_xy_0_0_0_0_z_zz_yy = buffer_2000_spdd[213];

    auto g_xy_0_0_0_0_z_zz_yz = buffer_2000_spdd[214];

    auto g_xy_0_0_0_0_z_zz_zz = buffer_2000_spdd[215];

    auto g_xz_0_0_0_0_x_xx_xx = buffer_2000_spdd[216];

    auto g_xz_0_0_0_0_x_xx_xy = buffer_2000_spdd[217];

    auto g_xz_0_0_0_0_x_xx_xz = buffer_2000_spdd[218];

    auto g_xz_0_0_0_0_x_xx_yy = buffer_2000_spdd[219];

    auto g_xz_0_0_0_0_x_xx_yz = buffer_2000_spdd[220];

    auto g_xz_0_0_0_0_x_xx_zz = buffer_2000_spdd[221];

    auto g_xz_0_0_0_0_x_xy_xx = buffer_2000_spdd[222];

    auto g_xz_0_0_0_0_x_xy_xy = buffer_2000_spdd[223];

    auto g_xz_0_0_0_0_x_xy_xz = buffer_2000_spdd[224];

    auto g_xz_0_0_0_0_x_xy_yy = buffer_2000_spdd[225];

    auto g_xz_0_0_0_0_x_xy_yz = buffer_2000_spdd[226];

    auto g_xz_0_0_0_0_x_xy_zz = buffer_2000_spdd[227];

    auto g_xz_0_0_0_0_x_xz_xx = buffer_2000_spdd[228];

    auto g_xz_0_0_0_0_x_xz_xy = buffer_2000_spdd[229];

    auto g_xz_0_0_0_0_x_xz_xz = buffer_2000_spdd[230];

    auto g_xz_0_0_0_0_x_xz_yy = buffer_2000_spdd[231];

    auto g_xz_0_0_0_0_x_xz_yz = buffer_2000_spdd[232];

    auto g_xz_0_0_0_0_x_xz_zz = buffer_2000_spdd[233];

    auto g_xz_0_0_0_0_x_yy_xx = buffer_2000_spdd[234];

    auto g_xz_0_0_0_0_x_yy_xy = buffer_2000_spdd[235];

    auto g_xz_0_0_0_0_x_yy_xz = buffer_2000_spdd[236];

    auto g_xz_0_0_0_0_x_yy_yy = buffer_2000_spdd[237];

    auto g_xz_0_0_0_0_x_yy_yz = buffer_2000_spdd[238];

    auto g_xz_0_0_0_0_x_yy_zz = buffer_2000_spdd[239];

    auto g_xz_0_0_0_0_x_yz_xx = buffer_2000_spdd[240];

    auto g_xz_0_0_0_0_x_yz_xy = buffer_2000_spdd[241];

    auto g_xz_0_0_0_0_x_yz_xz = buffer_2000_spdd[242];

    auto g_xz_0_0_0_0_x_yz_yy = buffer_2000_spdd[243];

    auto g_xz_0_0_0_0_x_yz_yz = buffer_2000_spdd[244];

    auto g_xz_0_0_0_0_x_yz_zz = buffer_2000_spdd[245];

    auto g_xz_0_0_0_0_x_zz_xx = buffer_2000_spdd[246];

    auto g_xz_0_0_0_0_x_zz_xy = buffer_2000_spdd[247];

    auto g_xz_0_0_0_0_x_zz_xz = buffer_2000_spdd[248];

    auto g_xz_0_0_0_0_x_zz_yy = buffer_2000_spdd[249];

    auto g_xz_0_0_0_0_x_zz_yz = buffer_2000_spdd[250];

    auto g_xz_0_0_0_0_x_zz_zz = buffer_2000_spdd[251];

    auto g_xz_0_0_0_0_y_xx_xx = buffer_2000_spdd[252];

    auto g_xz_0_0_0_0_y_xx_xy = buffer_2000_spdd[253];

    auto g_xz_0_0_0_0_y_xx_xz = buffer_2000_spdd[254];

    auto g_xz_0_0_0_0_y_xx_yy = buffer_2000_spdd[255];

    auto g_xz_0_0_0_0_y_xx_yz = buffer_2000_spdd[256];

    auto g_xz_0_0_0_0_y_xx_zz = buffer_2000_spdd[257];

    auto g_xz_0_0_0_0_y_xy_xx = buffer_2000_spdd[258];

    auto g_xz_0_0_0_0_y_xy_xy = buffer_2000_spdd[259];

    auto g_xz_0_0_0_0_y_xy_xz = buffer_2000_spdd[260];

    auto g_xz_0_0_0_0_y_xy_yy = buffer_2000_spdd[261];

    auto g_xz_0_0_0_0_y_xy_yz = buffer_2000_spdd[262];

    auto g_xz_0_0_0_0_y_xy_zz = buffer_2000_spdd[263];

    auto g_xz_0_0_0_0_y_xz_xx = buffer_2000_spdd[264];

    auto g_xz_0_0_0_0_y_xz_xy = buffer_2000_spdd[265];

    auto g_xz_0_0_0_0_y_xz_xz = buffer_2000_spdd[266];

    auto g_xz_0_0_0_0_y_xz_yy = buffer_2000_spdd[267];

    auto g_xz_0_0_0_0_y_xz_yz = buffer_2000_spdd[268];

    auto g_xz_0_0_0_0_y_xz_zz = buffer_2000_spdd[269];

    auto g_xz_0_0_0_0_y_yy_xx = buffer_2000_spdd[270];

    auto g_xz_0_0_0_0_y_yy_xy = buffer_2000_spdd[271];

    auto g_xz_0_0_0_0_y_yy_xz = buffer_2000_spdd[272];

    auto g_xz_0_0_0_0_y_yy_yy = buffer_2000_spdd[273];

    auto g_xz_0_0_0_0_y_yy_yz = buffer_2000_spdd[274];

    auto g_xz_0_0_0_0_y_yy_zz = buffer_2000_spdd[275];

    auto g_xz_0_0_0_0_y_yz_xx = buffer_2000_spdd[276];

    auto g_xz_0_0_0_0_y_yz_xy = buffer_2000_spdd[277];

    auto g_xz_0_0_0_0_y_yz_xz = buffer_2000_spdd[278];

    auto g_xz_0_0_0_0_y_yz_yy = buffer_2000_spdd[279];

    auto g_xz_0_0_0_0_y_yz_yz = buffer_2000_spdd[280];

    auto g_xz_0_0_0_0_y_yz_zz = buffer_2000_spdd[281];

    auto g_xz_0_0_0_0_y_zz_xx = buffer_2000_spdd[282];

    auto g_xz_0_0_0_0_y_zz_xy = buffer_2000_spdd[283];

    auto g_xz_0_0_0_0_y_zz_xz = buffer_2000_spdd[284];

    auto g_xz_0_0_0_0_y_zz_yy = buffer_2000_spdd[285];

    auto g_xz_0_0_0_0_y_zz_yz = buffer_2000_spdd[286];

    auto g_xz_0_0_0_0_y_zz_zz = buffer_2000_spdd[287];

    auto g_xz_0_0_0_0_z_xx_xx = buffer_2000_spdd[288];

    auto g_xz_0_0_0_0_z_xx_xy = buffer_2000_spdd[289];

    auto g_xz_0_0_0_0_z_xx_xz = buffer_2000_spdd[290];

    auto g_xz_0_0_0_0_z_xx_yy = buffer_2000_spdd[291];

    auto g_xz_0_0_0_0_z_xx_yz = buffer_2000_spdd[292];

    auto g_xz_0_0_0_0_z_xx_zz = buffer_2000_spdd[293];

    auto g_xz_0_0_0_0_z_xy_xx = buffer_2000_spdd[294];

    auto g_xz_0_0_0_0_z_xy_xy = buffer_2000_spdd[295];

    auto g_xz_0_0_0_0_z_xy_xz = buffer_2000_spdd[296];

    auto g_xz_0_0_0_0_z_xy_yy = buffer_2000_spdd[297];

    auto g_xz_0_0_0_0_z_xy_yz = buffer_2000_spdd[298];

    auto g_xz_0_0_0_0_z_xy_zz = buffer_2000_spdd[299];

    auto g_xz_0_0_0_0_z_xz_xx = buffer_2000_spdd[300];

    auto g_xz_0_0_0_0_z_xz_xy = buffer_2000_spdd[301];

    auto g_xz_0_0_0_0_z_xz_xz = buffer_2000_spdd[302];

    auto g_xz_0_0_0_0_z_xz_yy = buffer_2000_spdd[303];

    auto g_xz_0_0_0_0_z_xz_yz = buffer_2000_spdd[304];

    auto g_xz_0_0_0_0_z_xz_zz = buffer_2000_spdd[305];

    auto g_xz_0_0_0_0_z_yy_xx = buffer_2000_spdd[306];

    auto g_xz_0_0_0_0_z_yy_xy = buffer_2000_spdd[307];

    auto g_xz_0_0_0_0_z_yy_xz = buffer_2000_spdd[308];

    auto g_xz_0_0_0_0_z_yy_yy = buffer_2000_spdd[309];

    auto g_xz_0_0_0_0_z_yy_yz = buffer_2000_spdd[310];

    auto g_xz_0_0_0_0_z_yy_zz = buffer_2000_spdd[311];

    auto g_xz_0_0_0_0_z_yz_xx = buffer_2000_spdd[312];

    auto g_xz_0_0_0_0_z_yz_xy = buffer_2000_spdd[313];

    auto g_xz_0_0_0_0_z_yz_xz = buffer_2000_spdd[314];

    auto g_xz_0_0_0_0_z_yz_yy = buffer_2000_spdd[315];

    auto g_xz_0_0_0_0_z_yz_yz = buffer_2000_spdd[316];

    auto g_xz_0_0_0_0_z_yz_zz = buffer_2000_spdd[317];

    auto g_xz_0_0_0_0_z_zz_xx = buffer_2000_spdd[318];

    auto g_xz_0_0_0_0_z_zz_xy = buffer_2000_spdd[319];

    auto g_xz_0_0_0_0_z_zz_xz = buffer_2000_spdd[320];

    auto g_xz_0_0_0_0_z_zz_yy = buffer_2000_spdd[321];

    auto g_xz_0_0_0_0_z_zz_yz = buffer_2000_spdd[322];

    auto g_xz_0_0_0_0_z_zz_zz = buffer_2000_spdd[323];

    auto g_yy_0_0_0_0_x_xx_xx = buffer_2000_spdd[324];

    auto g_yy_0_0_0_0_x_xx_xy = buffer_2000_spdd[325];

    auto g_yy_0_0_0_0_x_xx_xz = buffer_2000_spdd[326];

    auto g_yy_0_0_0_0_x_xx_yy = buffer_2000_spdd[327];

    auto g_yy_0_0_0_0_x_xx_yz = buffer_2000_spdd[328];

    auto g_yy_0_0_0_0_x_xx_zz = buffer_2000_spdd[329];

    auto g_yy_0_0_0_0_x_xy_xx = buffer_2000_spdd[330];

    auto g_yy_0_0_0_0_x_xy_xy = buffer_2000_spdd[331];

    auto g_yy_0_0_0_0_x_xy_xz = buffer_2000_spdd[332];

    auto g_yy_0_0_0_0_x_xy_yy = buffer_2000_spdd[333];

    auto g_yy_0_0_0_0_x_xy_yz = buffer_2000_spdd[334];

    auto g_yy_0_0_0_0_x_xy_zz = buffer_2000_spdd[335];

    auto g_yy_0_0_0_0_x_xz_xx = buffer_2000_spdd[336];

    auto g_yy_0_0_0_0_x_xz_xy = buffer_2000_spdd[337];

    auto g_yy_0_0_0_0_x_xz_xz = buffer_2000_spdd[338];

    auto g_yy_0_0_0_0_x_xz_yy = buffer_2000_spdd[339];

    auto g_yy_0_0_0_0_x_xz_yz = buffer_2000_spdd[340];

    auto g_yy_0_0_0_0_x_xz_zz = buffer_2000_spdd[341];

    auto g_yy_0_0_0_0_x_yy_xx = buffer_2000_spdd[342];

    auto g_yy_0_0_0_0_x_yy_xy = buffer_2000_spdd[343];

    auto g_yy_0_0_0_0_x_yy_xz = buffer_2000_spdd[344];

    auto g_yy_0_0_0_0_x_yy_yy = buffer_2000_spdd[345];

    auto g_yy_0_0_0_0_x_yy_yz = buffer_2000_spdd[346];

    auto g_yy_0_0_0_0_x_yy_zz = buffer_2000_spdd[347];

    auto g_yy_0_0_0_0_x_yz_xx = buffer_2000_spdd[348];

    auto g_yy_0_0_0_0_x_yz_xy = buffer_2000_spdd[349];

    auto g_yy_0_0_0_0_x_yz_xz = buffer_2000_spdd[350];

    auto g_yy_0_0_0_0_x_yz_yy = buffer_2000_spdd[351];

    auto g_yy_0_0_0_0_x_yz_yz = buffer_2000_spdd[352];

    auto g_yy_0_0_0_0_x_yz_zz = buffer_2000_spdd[353];

    auto g_yy_0_0_0_0_x_zz_xx = buffer_2000_spdd[354];

    auto g_yy_0_0_0_0_x_zz_xy = buffer_2000_spdd[355];

    auto g_yy_0_0_0_0_x_zz_xz = buffer_2000_spdd[356];

    auto g_yy_0_0_0_0_x_zz_yy = buffer_2000_spdd[357];

    auto g_yy_0_0_0_0_x_zz_yz = buffer_2000_spdd[358];

    auto g_yy_0_0_0_0_x_zz_zz = buffer_2000_spdd[359];

    auto g_yy_0_0_0_0_y_xx_xx = buffer_2000_spdd[360];

    auto g_yy_0_0_0_0_y_xx_xy = buffer_2000_spdd[361];

    auto g_yy_0_0_0_0_y_xx_xz = buffer_2000_spdd[362];

    auto g_yy_0_0_0_0_y_xx_yy = buffer_2000_spdd[363];

    auto g_yy_0_0_0_0_y_xx_yz = buffer_2000_spdd[364];

    auto g_yy_0_0_0_0_y_xx_zz = buffer_2000_spdd[365];

    auto g_yy_0_0_0_0_y_xy_xx = buffer_2000_spdd[366];

    auto g_yy_0_0_0_0_y_xy_xy = buffer_2000_spdd[367];

    auto g_yy_0_0_0_0_y_xy_xz = buffer_2000_spdd[368];

    auto g_yy_0_0_0_0_y_xy_yy = buffer_2000_spdd[369];

    auto g_yy_0_0_0_0_y_xy_yz = buffer_2000_spdd[370];

    auto g_yy_0_0_0_0_y_xy_zz = buffer_2000_spdd[371];

    auto g_yy_0_0_0_0_y_xz_xx = buffer_2000_spdd[372];

    auto g_yy_0_0_0_0_y_xz_xy = buffer_2000_spdd[373];

    auto g_yy_0_0_0_0_y_xz_xz = buffer_2000_spdd[374];

    auto g_yy_0_0_0_0_y_xz_yy = buffer_2000_spdd[375];

    auto g_yy_0_0_0_0_y_xz_yz = buffer_2000_spdd[376];

    auto g_yy_0_0_0_0_y_xz_zz = buffer_2000_spdd[377];

    auto g_yy_0_0_0_0_y_yy_xx = buffer_2000_spdd[378];

    auto g_yy_0_0_0_0_y_yy_xy = buffer_2000_spdd[379];

    auto g_yy_0_0_0_0_y_yy_xz = buffer_2000_spdd[380];

    auto g_yy_0_0_0_0_y_yy_yy = buffer_2000_spdd[381];

    auto g_yy_0_0_0_0_y_yy_yz = buffer_2000_spdd[382];

    auto g_yy_0_0_0_0_y_yy_zz = buffer_2000_spdd[383];

    auto g_yy_0_0_0_0_y_yz_xx = buffer_2000_spdd[384];

    auto g_yy_0_0_0_0_y_yz_xy = buffer_2000_spdd[385];

    auto g_yy_0_0_0_0_y_yz_xz = buffer_2000_spdd[386];

    auto g_yy_0_0_0_0_y_yz_yy = buffer_2000_spdd[387];

    auto g_yy_0_0_0_0_y_yz_yz = buffer_2000_spdd[388];

    auto g_yy_0_0_0_0_y_yz_zz = buffer_2000_spdd[389];

    auto g_yy_0_0_0_0_y_zz_xx = buffer_2000_spdd[390];

    auto g_yy_0_0_0_0_y_zz_xy = buffer_2000_spdd[391];

    auto g_yy_0_0_0_0_y_zz_xz = buffer_2000_spdd[392];

    auto g_yy_0_0_0_0_y_zz_yy = buffer_2000_spdd[393];

    auto g_yy_0_0_0_0_y_zz_yz = buffer_2000_spdd[394];

    auto g_yy_0_0_0_0_y_zz_zz = buffer_2000_spdd[395];

    auto g_yy_0_0_0_0_z_xx_xx = buffer_2000_spdd[396];

    auto g_yy_0_0_0_0_z_xx_xy = buffer_2000_spdd[397];

    auto g_yy_0_0_0_0_z_xx_xz = buffer_2000_spdd[398];

    auto g_yy_0_0_0_0_z_xx_yy = buffer_2000_spdd[399];

    auto g_yy_0_0_0_0_z_xx_yz = buffer_2000_spdd[400];

    auto g_yy_0_0_0_0_z_xx_zz = buffer_2000_spdd[401];

    auto g_yy_0_0_0_0_z_xy_xx = buffer_2000_spdd[402];

    auto g_yy_0_0_0_0_z_xy_xy = buffer_2000_spdd[403];

    auto g_yy_0_0_0_0_z_xy_xz = buffer_2000_spdd[404];

    auto g_yy_0_0_0_0_z_xy_yy = buffer_2000_spdd[405];

    auto g_yy_0_0_0_0_z_xy_yz = buffer_2000_spdd[406];

    auto g_yy_0_0_0_0_z_xy_zz = buffer_2000_spdd[407];

    auto g_yy_0_0_0_0_z_xz_xx = buffer_2000_spdd[408];

    auto g_yy_0_0_0_0_z_xz_xy = buffer_2000_spdd[409];

    auto g_yy_0_0_0_0_z_xz_xz = buffer_2000_spdd[410];

    auto g_yy_0_0_0_0_z_xz_yy = buffer_2000_spdd[411];

    auto g_yy_0_0_0_0_z_xz_yz = buffer_2000_spdd[412];

    auto g_yy_0_0_0_0_z_xz_zz = buffer_2000_spdd[413];

    auto g_yy_0_0_0_0_z_yy_xx = buffer_2000_spdd[414];

    auto g_yy_0_0_0_0_z_yy_xy = buffer_2000_spdd[415];

    auto g_yy_0_0_0_0_z_yy_xz = buffer_2000_spdd[416];

    auto g_yy_0_0_0_0_z_yy_yy = buffer_2000_spdd[417];

    auto g_yy_0_0_0_0_z_yy_yz = buffer_2000_spdd[418];

    auto g_yy_0_0_0_0_z_yy_zz = buffer_2000_spdd[419];

    auto g_yy_0_0_0_0_z_yz_xx = buffer_2000_spdd[420];

    auto g_yy_0_0_0_0_z_yz_xy = buffer_2000_spdd[421];

    auto g_yy_0_0_0_0_z_yz_xz = buffer_2000_spdd[422];

    auto g_yy_0_0_0_0_z_yz_yy = buffer_2000_spdd[423];

    auto g_yy_0_0_0_0_z_yz_yz = buffer_2000_spdd[424];

    auto g_yy_0_0_0_0_z_yz_zz = buffer_2000_spdd[425];

    auto g_yy_0_0_0_0_z_zz_xx = buffer_2000_spdd[426];

    auto g_yy_0_0_0_0_z_zz_xy = buffer_2000_spdd[427];

    auto g_yy_0_0_0_0_z_zz_xz = buffer_2000_spdd[428];

    auto g_yy_0_0_0_0_z_zz_yy = buffer_2000_spdd[429];

    auto g_yy_0_0_0_0_z_zz_yz = buffer_2000_spdd[430];

    auto g_yy_0_0_0_0_z_zz_zz = buffer_2000_spdd[431];

    auto g_yz_0_0_0_0_x_xx_xx = buffer_2000_spdd[432];

    auto g_yz_0_0_0_0_x_xx_xy = buffer_2000_spdd[433];

    auto g_yz_0_0_0_0_x_xx_xz = buffer_2000_spdd[434];

    auto g_yz_0_0_0_0_x_xx_yy = buffer_2000_spdd[435];

    auto g_yz_0_0_0_0_x_xx_yz = buffer_2000_spdd[436];

    auto g_yz_0_0_0_0_x_xx_zz = buffer_2000_spdd[437];

    auto g_yz_0_0_0_0_x_xy_xx = buffer_2000_spdd[438];

    auto g_yz_0_0_0_0_x_xy_xy = buffer_2000_spdd[439];

    auto g_yz_0_0_0_0_x_xy_xz = buffer_2000_spdd[440];

    auto g_yz_0_0_0_0_x_xy_yy = buffer_2000_spdd[441];

    auto g_yz_0_0_0_0_x_xy_yz = buffer_2000_spdd[442];

    auto g_yz_0_0_0_0_x_xy_zz = buffer_2000_spdd[443];

    auto g_yz_0_0_0_0_x_xz_xx = buffer_2000_spdd[444];

    auto g_yz_0_0_0_0_x_xz_xy = buffer_2000_spdd[445];

    auto g_yz_0_0_0_0_x_xz_xz = buffer_2000_spdd[446];

    auto g_yz_0_0_0_0_x_xz_yy = buffer_2000_spdd[447];

    auto g_yz_0_0_0_0_x_xz_yz = buffer_2000_spdd[448];

    auto g_yz_0_0_0_0_x_xz_zz = buffer_2000_spdd[449];

    auto g_yz_0_0_0_0_x_yy_xx = buffer_2000_spdd[450];

    auto g_yz_0_0_0_0_x_yy_xy = buffer_2000_spdd[451];

    auto g_yz_0_0_0_0_x_yy_xz = buffer_2000_spdd[452];

    auto g_yz_0_0_0_0_x_yy_yy = buffer_2000_spdd[453];

    auto g_yz_0_0_0_0_x_yy_yz = buffer_2000_spdd[454];

    auto g_yz_0_0_0_0_x_yy_zz = buffer_2000_spdd[455];

    auto g_yz_0_0_0_0_x_yz_xx = buffer_2000_spdd[456];

    auto g_yz_0_0_0_0_x_yz_xy = buffer_2000_spdd[457];

    auto g_yz_0_0_0_0_x_yz_xz = buffer_2000_spdd[458];

    auto g_yz_0_0_0_0_x_yz_yy = buffer_2000_spdd[459];

    auto g_yz_0_0_0_0_x_yz_yz = buffer_2000_spdd[460];

    auto g_yz_0_0_0_0_x_yz_zz = buffer_2000_spdd[461];

    auto g_yz_0_0_0_0_x_zz_xx = buffer_2000_spdd[462];

    auto g_yz_0_0_0_0_x_zz_xy = buffer_2000_spdd[463];

    auto g_yz_0_0_0_0_x_zz_xz = buffer_2000_spdd[464];

    auto g_yz_0_0_0_0_x_zz_yy = buffer_2000_spdd[465];

    auto g_yz_0_0_0_0_x_zz_yz = buffer_2000_spdd[466];

    auto g_yz_0_0_0_0_x_zz_zz = buffer_2000_spdd[467];

    auto g_yz_0_0_0_0_y_xx_xx = buffer_2000_spdd[468];

    auto g_yz_0_0_0_0_y_xx_xy = buffer_2000_spdd[469];

    auto g_yz_0_0_0_0_y_xx_xz = buffer_2000_spdd[470];

    auto g_yz_0_0_0_0_y_xx_yy = buffer_2000_spdd[471];

    auto g_yz_0_0_0_0_y_xx_yz = buffer_2000_spdd[472];

    auto g_yz_0_0_0_0_y_xx_zz = buffer_2000_spdd[473];

    auto g_yz_0_0_0_0_y_xy_xx = buffer_2000_spdd[474];

    auto g_yz_0_0_0_0_y_xy_xy = buffer_2000_spdd[475];

    auto g_yz_0_0_0_0_y_xy_xz = buffer_2000_spdd[476];

    auto g_yz_0_0_0_0_y_xy_yy = buffer_2000_spdd[477];

    auto g_yz_0_0_0_0_y_xy_yz = buffer_2000_spdd[478];

    auto g_yz_0_0_0_0_y_xy_zz = buffer_2000_spdd[479];

    auto g_yz_0_0_0_0_y_xz_xx = buffer_2000_spdd[480];

    auto g_yz_0_0_0_0_y_xz_xy = buffer_2000_spdd[481];

    auto g_yz_0_0_0_0_y_xz_xz = buffer_2000_spdd[482];

    auto g_yz_0_0_0_0_y_xz_yy = buffer_2000_spdd[483];

    auto g_yz_0_0_0_0_y_xz_yz = buffer_2000_spdd[484];

    auto g_yz_0_0_0_0_y_xz_zz = buffer_2000_spdd[485];

    auto g_yz_0_0_0_0_y_yy_xx = buffer_2000_spdd[486];

    auto g_yz_0_0_0_0_y_yy_xy = buffer_2000_spdd[487];

    auto g_yz_0_0_0_0_y_yy_xz = buffer_2000_spdd[488];

    auto g_yz_0_0_0_0_y_yy_yy = buffer_2000_spdd[489];

    auto g_yz_0_0_0_0_y_yy_yz = buffer_2000_spdd[490];

    auto g_yz_0_0_0_0_y_yy_zz = buffer_2000_spdd[491];

    auto g_yz_0_0_0_0_y_yz_xx = buffer_2000_spdd[492];

    auto g_yz_0_0_0_0_y_yz_xy = buffer_2000_spdd[493];

    auto g_yz_0_0_0_0_y_yz_xz = buffer_2000_spdd[494];

    auto g_yz_0_0_0_0_y_yz_yy = buffer_2000_spdd[495];

    auto g_yz_0_0_0_0_y_yz_yz = buffer_2000_spdd[496];

    auto g_yz_0_0_0_0_y_yz_zz = buffer_2000_spdd[497];

    auto g_yz_0_0_0_0_y_zz_xx = buffer_2000_spdd[498];

    auto g_yz_0_0_0_0_y_zz_xy = buffer_2000_spdd[499];

    auto g_yz_0_0_0_0_y_zz_xz = buffer_2000_spdd[500];

    auto g_yz_0_0_0_0_y_zz_yy = buffer_2000_spdd[501];

    auto g_yz_0_0_0_0_y_zz_yz = buffer_2000_spdd[502];

    auto g_yz_0_0_0_0_y_zz_zz = buffer_2000_spdd[503];

    auto g_yz_0_0_0_0_z_xx_xx = buffer_2000_spdd[504];

    auto g_yz_0_0_0_0_z_xx_xy = buffer_2000_spdd[505];

    auto g_yz_0_0_0_0_z_xx_xz = buffer_2000_spdd[506];

    auto g_yz_0_0_0_0_z_xx_yy = buffer_2000_spdd[507];

    auto g_yz_0_0_0_0_z_xx_yz = buffer_2000_spdd[508];

    auto g_yz_0_0_0_0_z_xx_zz = buffer_2000_spdd[509];

    auto g_yz_0_0_0_0_z_xy_xx = buffer_2000_spdd[510];

    auto g_yz_0_0_0_0_z_xy_xy = buffer_2000_spdd[511];

    auto g_yz_0_0_0_0_z_xy_xz = buffer_2000_spdd[512];

    auto g_yz_0_0_0_0_z_xy_yy = buffer_2000_spdd[513];

    auto g_yz_0_0_0_0_z_xy_yz = buffer_2000_spdd[514];

    auto g_yz_0_0_0_0_z_xy_zz = buffer_2000_spdd[515];

    auto g_yz_0_0_0_0_z_xz_xx = buffer_2000_spdd[516];

    auto g_yz_0_0_0_0_z_xz_xy = buffer_2000_spdd[517];

    auto g_yz_0_0_0_0_z_xz_xz = buffer_2000_spdd[518];

    auto g_yz_0_0_0_0_z_xz_yy = buffer_2000_spdd[519];

    auto g_yz_0_0_0_0_z_xz_yz = buffer_2000_spdd[520];

    auto g_yz_0_0_0_0_z_xz_zz = buffer_2000_spdd[521];

    auto g_yz_0_0_0_0_z_yy_xx = buffer_2000_spdd[522];

    auto g_yz_0_0_0_0_z_yy_xy = buffer_2000_spdd[523];

    auto g_yz_0_0_0_0_z_yy_xz = buffer_2000_spdd[524];

    auto g_yz_0_0_0_0_z_yy_yy = buffer_2000_spdd[525];

    auto g_yz_0_0_0_0_z_yy_yz = buffer_2000_spdd[526];

    auto g_yz_0_0_0_0_z_yy_zz = buffer_2000_spdd[527];

    auto g_yz_0_0_0_0_z_yz_xx = buffer_2000_spdd[528];

    auto g_yz_0_0_0_0_z_yz_xy = buffer_2000_spdd[529];

    auto g_yz_0_0_0_0_z_yz_xz = buffer_2000_spdd[530];

    auto g_yz_0_0_0_0_z_yz_yy = buffer_2000_spdd[531];

    auto g_yz_0_0_0_0_z_yz_yz = buffer_2000_spdd[532];

    auto g_yz_0_0_0_0_z_yz_zz = buffer_2000_spdd[533];

    auto g_yz_0_0_0_0_z_zz_xx = buffer_2000_spdd[534];

    auto g_yz_0_0_0_0_z_zz_xy = buffer_2000_spdd[535];

    auto g_yz_0_0_0_0_z_zz_xz = buffer_2000_spdd[536];

    auto g_yz_0_0_0_0_z_zz_yy = buffer_2000_spdd[537];

    auto g_yz_0_0_0_0_z_zz_yz = buffer_2000_spdd[538];

    auto g_yz_0_0_0_0_z_zz_zz = buffer_2000_spdd[539];

    auto g_zz_0_0_0_0_x_xx_xx = buffer_2000_spdd[540];

    auto g_zz_0_0_0_0_x_xx_xy = buffer_2000_spdd[541];

    auto g_zz_0_0_0_0_x_xx_xz = buffer_2000_spdd[542];

    auto g_zz_0_0_0_0_x_xx_yy = buffer_2000_spdd[543];

    auto g_zz_0_0_0_0_x_xx_yz = buffer_2000_spdd[544];

    auto g_zz_0_0_0_0_x_xx_zz = buffer_2000_spdd[545];

    auto g_zz_0_0_0_0_x_xy_xx = buffer_2000_spdd[546];

    auto g_zz_0_0_0_0_x_xy_xy = buffer_2000_spdd[547];

    auto g_zz_0_0_0_0_x_xy_xz = buffer_2000_spdd[548];

    auto g_zz_0_0_0_0_x_xy_yy = buffer_2000_spdd[549];

    auto g_zz_0_0_0_0_x_xy_yz = buffer_2000_spdd[550];

    auto g_zz_0_0_0_0_x_xy_zz = buffer_2000_spdd[551];

    auto g_zz_0_0_0_0_x_xz_xx = buffer_2000_spdd[552];

    auto g_zz_0_0_0_0_x_xz_xy = buffer_2000_spdd[553];

    auto g_zz_0_0_0_0_x_xz_xz = buffer_2000_spdd[554];

    auto g_zz_0_0_0_0_x_xz_yy = buffer_2000_spdd[555];

    auto g_zz_0_0_0_0_x_xz_yz = buffer_2000_spdd[556];

    auto g_zz_0_0_0_0_x_xz_zz = buffer_2000_spdd[557];

    auto g_zz_0_0_0_0_x_yy_xx = buffer_2000_spdd[558];

    auto g_zz_0_0_0_0_x_yy_xy = buffer_2000_spdd[559];

    auto g_zz_0_0_0_0_x_yy_xz = buffer_2000_spdd[560];

    auto g_zz_0_0_0_0_x_yy_yy = buffer_2000_spdd[561];

    auto g_zz_0_0_0_0_x_yy_yz = buffer_2000_spdd[562];

    auto g_zz_0_0_0_0_x_yy_zz = buffer_2000_spdd[563];

    auto g_zz_0_0_0_0_x_yz_xx = buffer_2000_spdd[564];

    auto g_zz_0_0_0_0_x_yz_xy = buffer_2000_spdd[565];

    auto g_zz_0_0_0_0_x_yz_xz = buffer_2000_spdd[566];

    auto g_zz_0_0_0_0_x_yz_yy = buffer_2000_spdd[567];

    auto g_zz_0_0_0_0_x_yz_yz = buffer_2000_spdd[568];

    auto g_zz_0_0_0_0_x_yz_zz = buffer_2000_spdd[569];

    auto g_zz_0_0_0_0_x_zz_xx = buffer_2000_spdd[570];

    auto g_zz_0_0_0_0_x_zz_xy = buffer_2000_spdd[571];

    auto g_zz_0_0_0_0_x_zz_xz = buffer_2000_spdd[572];

    auto g_zz_0_0_0_0_x_zz_yy = buffer_2000_spdd[573];

    auto g_zz_0_0_0_0_x_zz_yz = buffer_2000_spdd[574];

    auto g_zz_0_0_0_0_x_zz_zz = buffer_2000_spdd[575];

    auto g_zz_0_0_0_0_y_xx_xx = buffer_2000_spdd[576];

    auto g_zz_0_0_0_0_y_xx_xy = buffer_2000_spdd[577];

    auto g_zz_0_0_0_0_y_xx_xz = buffer_2000_spdd[578];

    auto g_zz_0_0_0_0_y_xx_yy = buffer_2000_spdd[579];

    auto g_zz_0_0_0_0_y_xx_yz = buffer_2000_spdd[580];

    auto g_zz_0_0_0_0_y_xx_zz = buffer_2000_spdd[581];

    auto g_zz_0_0_0_0_y_xy_xx = buffer_2000_spdd[582];

    auto g_zz_0_0_0_0_y_xy_xy = buffer_2000_spdd[583];

    auto g_zz_0_0_0_0_y_xy_xz = buffer_2000_spdd[584];

    auto g_zz_0_0_0_0_y_xy_yy = buffer_2000_spdd[585];

    auto g_zz_0_0_0_0_y_xy_yz = buffer_2000_spdd[586];

    auto g_zz_0_0_0_0_y_xy_zz = buffer_2000_spdd[587];

    auto g_zz_0_0_0_0_y_xz_xx = buffer_2000_spdd[588];

    auto g_zz_0_0_0_0_y_xz_xy = buffer_2000_spdd[589];

    auto g_zz_0_0_0_0_y_xz_xz = buffer_2000_spdd[590];

    auto g_zz_0_0_0_0_y_xz_yy = buffer_2000_spdd[591];

    auto g_zz_0_0_0_0_y_xz_yz = buffer_2000_spdd[592];

    auto g_zz_0_0_0_0_y_xz_zz = buffer_2000_spdd[593];

    auto g_zz_0_0_0_0_y_yy_xx = buffer_2000_spdd[594];

    auto g_zz_0_0_0_0_y_yy_xy = buffer_2000_spdd[595];

    auto g_zz_0_0_0_0_y_yy_xz = buffer_2000_spdd[596];

    auto g_zz_0_0_0_0_y_yy_yy = buffer_2000_spdd[597];

    auto g_zz_0_0_0_0_y_yy_yz = buffer_2000_spdd[598];

    auto g_zz_0_0_0_0_y_yy_zz = buffer_2000_spdd[599];

    auto g_zz_0_0_0_0_y_yz_xx = buffer_2000_spdd[600];

    auto g_zz_0_0_0_0_y_yz_xy = buffer_2000_spdd[601];

    auto g_zz_0_0_0_0_y_yz_xz = buffer_2000_spdd[602];

    auto g_zz_0_0_0_0_y_yz_yy = buffer_2000_spdd[603];

    auto g_zz_0_0_0_0_y_yz_yz = buffer_2000_spdd[604];

    auto g_zz_0_0_0_0_y_yz_zz = buffer_2000_spdd[605];

    auto g_zz_0_0_0_0_y_zz_xx = buffer_2000_spdd[606];

    auto g_zz_0_0_0_0_y_zz_xy = buffer_2000_spdd[607];

    auto g_zz_0_0_0_0_y_zz_xz = buffer_2000_spdd[608];

    auto g_zz_0_0_0_0_y_zz_yy = buffer_2000_spdd[609];

    auto g_zz_0_0_0_0_y_zz_yz = buffer_2000_spdd[610];

    auto g_zz_0_0_0_0_y_zz_zz = buffer_2000_spdd[611];

    auto g_zz_0_0_0_0_z_xx_xx = buffer_2000_spdd[612];

    auto g_zz_0_0_0_0_z_xx_xy = buffer_2000_spdd[613];

    auto g_zz_0_0_0_0_z_xx_xz = buffer_2000_spdd[614];

    auto g_zz_0_0_0_0_z_xx_yy = buffer_2000_spdd[615];

    auto g_zz_0_0_0_0_z_xx_yz = buffer_2000_spdd[616];

    auto g_zz_0_0_0_0_z_xx_zz = buffer_2000_spdd[617];

    auto g_zz_0_0_0_0_z_xy_xx = buffer_2000_spdd[618];

    auto g_zz_0_0_0_0_z_xy_xy = buffer_2000_spdd[619];

    auto g_zz_0_0_0_0_z_xy_xz = buffer_2000_spdd[620];

    auto g_zz_0_0_0_0_z_xy_yy = buffer_2000_spdd[621];

    auto g_zz_0_0_0_0_z_xy_yz = buffer_2000_spdd[622];

    auto g_zz_0_0_0_0_z_xy_zz = buffer_2000_spdd[623];

    auto g_zz_0_0_0_0_z_xz_xx = buffer_2000_spdd[624];

    auto g_zz_0_0_0_0_z_xz_xy = buffer_2000_spdd[625];

    auto g_zz_0_0_0_0_z_xz_xz = buffer_2000_spdd[626];

    auto g_zz_0_0_0_0_z_xz_yy = buffer_2000_spdd[627];

    auto g_zz_0_0_0_0_z_xz_yz = buffer_2000_spdd[628];

    auto g_zz_0_0_0_0_z_xz_zz = buffer_2000_spdd[629];

    auto g_zz_0_0_0_0_z_yy_xx = buffer_2000_spdd[630];

    auto g_zz_0_0_0_0_z_yy_xy = buffer_2000_spdd[631];

    auto g_zz_0_0_0_0_z_yy_xz = buffer_2000_spdd[632];

    auto g_zz_0_0_0_0_z_yy_yy = buffer_2000_spdd[633];

    auto g_zz_0_0_0_0_z_yy_yz = buffer_2000_spdd[634];

    auto g_zz_0_0_0_0_z_yy_zz = buffer_2000_spdd[635];

    auto g_zz_0_0_0_0_z_yz_xx = buffer_2000_spdd[636];

    auto g_zz_0_0_0_0_z_yz_xy = buffer_2000_spdd[637];

    auto g_zz_0_0_0_0_z_yz_xz = buffer_2000_spdd[638];

    auto g_zz_0_0_0_0_z_yz_yy = buffer_2000_spdd[639];

    auto g_zz_0_0_0_0_z_yz_yz = buffer_2000_spdd[640];

    auto g_zz_0_0_0_0_z_yz_zz = buffer_2000_spdd[641];

    auto g_zz_0_0_0_0_z_zz_xx = buffer_2000_spdd[642];

    auto g_zz_0_0_0_0_z_zz_xy = buffer_2000_spdd[643];

    auto g_zz_0_0_0_0_z_zz_xz = buffer_2000_spdd[644];

    auto g_zz_0_0_0_0_z_zz_yy = buffer_2000_spdd[645];

    auto g_zz_0_0_0_0_z_zz_yz = buffer_2000_spdd[646];

    auto g_zz_0_0_0_0_z_zz_zz = buffer_2000_spdd[647];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_x_xx_xx, g_0_x_xx_xy, g_0_x_xx_xz, g_0_x_xx_yy, g_0_x_xx_yz, g_0_x_xx_zz, g_xx_0_0_0_0_x_xx_xx, g_xx_0_0_0_0_x_xx_xy, g_xx_0_0_0_0_x_xx_xz, g_xx_0_0_0_0_x_xx_yy, g_xx_0_0_0_0_x_xx_yz, g_xx_0_0_0_0_x_xx_zz, g_xx_x_xx_xx, g_xx_x_xx_xy, g_xx_x_xx_xz, g_xx_x_xx_yy, g_xx_x_xx_yz, g_xx_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_x_xx_xx[i] = -2.0 * g_0_x_xx_xx[i] * a_exp + 4.0 * g_xx_x_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_xx_xy[i] = -2.0 * g_0_x_xx_xy[i] * a_exp + 4.0 * g_xx_x_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_xx_xz[i] = -2.0 * g_0_x_xx_xz[i] * a_exp + 4.0 * g_xx_x_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_xx_yy[i] = -2.0 * g_0_x_xx_yy[i] * a_exp + 4.0 * g_xx_x_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_xx_yz[i] = -2.0 * g_0_x_xx_yz[i] * a_exp + 4.0 * g_xx_x_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_xx_zz[i] = -2.0 * g_0_x_xx_zz[i] * a_exp + 4.0 * g_xx_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_x_xy_xx, g_0_x_xy_xy, g_0_x_xy_xz, g_0_x_xy_yy, g_0_x_xy_yz, g_0_x_xy_zz, g_xx_0_0_0_0_x_xy_xx, g_xx_0_0_0_0_x_xy_xy, g_xx_0_0_0_0_x_xy_xz, g_xx_0_0_0_0_x_xy_yy, g_xx_0_0_0_0_x_xy_yz, g_xx_0_0_0_0_x_xy_zz, g_xx_x_xy_xx, g_xx_x_xy_xy, g_xx_x_xy_xz, g_xx_x_xy_yy, g_xx_x_xy_yz, g_xx_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_x_xy_xx[i] = -2.0 * g_0_x_xy_xx[i] * a_exp + 4.0 * g_xx_x_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_xy_xy[i] = -2.0 * g_0_x_xy_xy[i] * a_exp + 4.0 * g_xx_x_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_xy_xz[i] = -2.0 * g_0_x_xy_xz[i] * a_exp + 4.0 * g_xx_x_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_xy_yy[i] = -2.0 * g_0_x_xy_yy[i] * a_exp + 4.0 * g_xx_x_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_xy_yz[i] = -2.0 * g_0_x_xy_yz[i] * a_exp + 4.0 * g_xx_x_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_xy_zz[i] = -2.0 * g_0_x_xy_zz[i] * a_exp + 4.0 * g_xx_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_x_xz_xx, g_0_x_xz_xy, g_0_x_xz_xz, g_0_x_xz_yy, g_0_x_xz_yz, g_0_x_xz_zz, g_xx_0_0_0_0_x_xz_xx, g_xx_0_0_0_0_x_xz_xy, g_xx_0_0_0_0_x_xz_xz, g_xx_0_0_0_0_x_xz_yy, g_xx_0_0_0_0_x_xz_yz, g_xx_0_0_0_0_x_xz_zz, g_xx_x_xz_xx, g_xx_x_xz_xy, g_xx_x_xz_xz, g_xx_x_xz_yy, g_xx_x_xz_yz, g_xx_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_x_xz_xx[i] = -2.0 * g_0_x_xz_xx[i] * a_exp + 4.0 * g_xx_x_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_xz_xy[i] = -2.0 * g_0_x_xz_xy[i] * a_exp + 4.0 * g_xx_x_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_xz_xz[i] = -2.0 * g_0_x_xz_xz[i] * a_exp + 4.0 * g_xx_x_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_xz_yy[i] = -2.0 * g_0_x_xz_yy[i] * a_exp + 4.0 * g_xx_x_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_xz_yz[i] = -2.0 * g_0_x_xz_yz[i] * a_exp + 4.0 * g_xx_x_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_xz_zz[i] = -2.0 * g_0_x_xz_zz[i] * a_exp + 4.0 * g_xx_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_x_yy_xx, g_0_x_yy_xy, g_0_x_yy_xz, g_0_x_yy_yy, g_0_x_yy_yz, g_0_x_yy_zz, g_xx_0_0_0_0_x_yy_xx, g_xx_0_0_0_0_x_yy_xy, g_xx_0_0_0_0_x_yy_xz, g_xx_0_0_0_0_x_yy_yy, g_xx_0_0_0_0_x_yy_yz, g_xx_0_0_0_0_x_yy_zz, g_xx_x_yy_xx, g_xx_x_yy_xy, g_xx_x_yy_xz, g_xx_x_yy_yy, g_xx_x_yy_yz, g_xx_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_x_yy_xx[i] = -2.0 * g_0_x_yy_xx[i] * a_exp + 4.0 * g_xx_x_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_yy_xy[i] = -2.0 * g_0_x_yy_xy[i] * a_exp + 4.0 * g_xx_x_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_yy_xz[i] = -2.0 * g_0_x_yy_xz[i] * a_exp + 4.0 * g_xx_x_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_yy_yy[i] = -2.0 * g_0_x_yy_yy[i] * a_exp + 4.0 * g_xx_x_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_yy_yz[i] = -2.0 * g_0_x_yy_yz[i] * a_exp + 4.0 * g_xx_x_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_yy_zz[i] = -2.0 * g_0_x_yy_zz[i] * a_exp + 4.0 * g_xx_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_x_yz_xx, g_0_x_yz_xy, g_0_x_yz_xz, g_0_x_yz_yy, g_0_x_yz_yz, g_0_x_yz_zz, g_xx_0_0_0_0_x_yz_xx, g_xx_0_0_0_0_x_yz_xy, g_xx_0_0_0_0_x_yz_xz, g_xx_0_0_0_0_x_yz_yy, g_xx_0_0_0_0_x_yz_yz, g_xx_0_0_0_0_x_yz_zz, g_xx_x_yz_xx, g_xx_x_yz_xy, g_xx_x_yz_xz, g_xx_x_yz_yy, g_xx_x_yz_yz, g_xx_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_x_yz_xx[i] = -2.0 * g_0_x_yz_xx[i] * a_exp + 4.0 * g_xx_x_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_yz_xy[i] = -2.0 * g_0_x_yz_xy[i] * a_exp + 4.0 * g_xx_x_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_yz_xz[i] = -2.0 * g_0_x_yz_xz[i] * a_exp + 4.0 * g_xx_x_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_yz_yy[i] = -2.0 * g_0_x_yz_yy[i] * a_exp + 4.0 * g_xx_x_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_yz_yz[i] = -2.0 * g_0_x_yz_yz[i] * a_exp + 4.0 * g_xx_x_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_yz_zz[i] = -2.0 * g_0_x_yz_zz[i] * a_exp + 4.0 * g_xx_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_x_zz_xx, g_0_x_zz_xy, g_0_x_zz_xz, g_0_x_zz_yy, g_0_x_zz_yz, g_0_x_zz_zz, g_xx_0_0_0_0_x_zz_xx, g_xx_0_0_0_0_x_zz_xy, g_xx_0_0_0_0_x_zz_xz, g_xx_0_0_0_0_x_zz_yy, g_xx_0_0_0_0_x_zz_yz, g_xx_0_0_0_0_x_zz_zz, g_xx_x_zz_xx, g_xx_x_zz_xy, g_xx_x_zz_xz, g_xx_x_zz_yy, g_xx_x_zz_yz, g_xx_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_x_zz_xx[i] = -2.0 * g_0_x_zz_xx[i] * a_exp + 4.0 * g_xx_x_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_zz_xy[i] = -2.0 * g_0_x_zz_xy[i] * a_exp + 4.0 * g_xx_x_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_zz_xz[i] = -2.0 * g_0_x_zz_xz[i] * a_exp + 4.0 * g_xx_x_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_zz_yy[i] = -2.0 * g_0_x_zz_yy[i] * a_exp + 4.0 * g_xx_x_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_zz_yz[i] = -2.0 * g_0_x_zz_yz[i] * a_exp + 4.0 * g_xx_x_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_zz_zz[i] = -2.0 * g_0_x_zz_zz[i] * a_exp + 4.0 * g_xx_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_0_y_xx_xx, g_0_y_xx_xy, g_0_y_xx_xz, g_0_y_xx_yy, g_0_y_xx_yz, g_0_y_xx_zz, g_xx_0_0_0_0_y_xx_xx, g_xx_0_0_0_0_y_xx_xy, g_xx_0_0_0_0_y_xx_xz, g_xx_0_0_0_0_y_xx_yy, g_xx_0_0_0_0_y_xx_yz, g_xx_0_0_0_0_y_xx_zz, g_xx_y_xx_xx, g_xx_y_xx_xy, g_xx_y_xx_xz, g_xx_y_xx_yy, g_xx_y_xx_yz, g_xx_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_y_xx_xx[i] = -2.0 * g_0_y_xx_xx[i] * a_exp + 4.0 * g_xx_y_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_xx_xy[i] = -2.0 * g_0_y_xx_xy[i] * a_exp + 4.0 * g_xx_y_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_xx_xz[i] = -2.0 * g_0_y_xx_xz[i] * a_exp + 4.0 * g_xx_y_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_xx_yy[i] = -2.0 * g_0_y_xx_yy[i] * a_exp + 4.0 * g_xx_y_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_xx_yz[i] = -2.0 * g_0_y_xx_yz[i] * a_exp + 4.0 * g_xx_y_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_xx_zz[i] = -2.0 * g_0_y_xx_zz[i] * a_exp + 4.0 * g_xx_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_0_y_xy_xx, g_0_y_xy_xy, g_0_y_xy_xz, g_0_y_xy_yy, g_0_y_xy_yz, g_0_y_xy_zz, g_xx_0_0_0_0_y_xy_xx, g_xx_0_0_0_0_y_xy_xy, g_xx_0_0_0_0_y_xy_xz, g_xx_0_0_0_0_y_xy_yy, g_xx_0_0_0_0_y_xy_yz, g_xx_0_0_0_0_y_xy_zz, g_xx_y_xy_xx, g_xx_y_xy_xy, g_xx_y_xy_xz, g_xx_y_xy_yy, g_xx_y_xy_yz, g_xx_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_y_xy_xx[i] = -2.0 * g_0_y_xy_xx[i] * a_exp + 4.0 * g_xx_y_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_xy_xy[i] = -2.0 * g_0_y_xy_xy[i] * a_exp + 4.0 * g_xx_y_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_xy_xz[i] = -2.0 * g_0_y_xy_xz[i] * a_exp + 4.0 * g_xx_y_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_xy_yy[i] = -2.0 * g_0_y_xy_yy[i] * a_exp + 4.0 * g_xx_y_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_xy_yz[i] = -2.0 * g_0_y_xy_yz[i] * a_exp + 4.0 * g_xx_y_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_xy_zz[i] = -2.0 * g_0_y_xy_zz[i] * a_exp + 4.0 * g_xx_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_0_y_xz_xx, g_0_y_xz_xy, g_0_y_xz_xz, g_0_y_xz_yy, g_0_y_xz_yz, g_0_y_xz_zz, g_xx_0_0_0_0_y_xz_xx, g_xx_0_0_0_0_y_xz_xy, g_xx_0_0_0_0_y_xz_xz, g_xx_0_0_0_0_y_xz_yy, g_xx_0_0_0_0_y_xz_yz, g_xx_0_0_0_0_y_xz_zz, g_xx_y_xz_xx, g_xx_y_xz_xy, g_xx_y_xz_xz, g_xx_y_xz_yy, g_xx_y_xz_yz, g_xx_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_y_xz_xx[i] = -2.0 * g_0_y_xz_xx[i] * a_exp + 4.0 * g_xx_y_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_xz_xy[i] = -2.0 * g_0_y_xz_xy[i] * a_exp + 4.0 * g_xx_y_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_xz_xz[i] = -2.0 * g_0_y_xz_xz[i] * a_exp + 4.0 * g_xx_y_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_xz_yy[i] = -2.0 * g_0_y_xz_yy[i] * a_exp + 4.0 * g_xx_y_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_xz_yz[i] = -2.0 * g_0_y_xz_yz[i] * a_exp + 4.0 * g_xx_y_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_xz_zz[i] = -2.0 * g_0_y_xz_zz[i] * a_exp + 4.0 * g_xx_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_0_y_yy_xx, g_0_y_yy_xy, g_0_y_yy_xz, g_0_y_yy_yy, g_0_y_yy_yz, g_0_y_yy_zz, g_xx_0_0_0_0_y_yy_xx, g_xx_0_0_0_0_y_yy_xy, g_xx_0_0_0_0_y_yy_xz, g_xx_0_0_0_0_y_yy_yy, g_xx_0_0_0_0_y_yy_yz, g_xx_0_0_0_0_y_yy_zz, g_xx_y_yy_xx, g_xx_y_yy_xy, g_xx_y_yy_xz, g_xx_y_yy_yy, g_xx_y_yy_yz, g_xx_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_y_yy_xx[i] = -2.0 * g_0_y_yy_xx[i] * a_exp + 4.0 * g_xx_y_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_yy_xy[i] = -2.0 * g_0_y_yy_xy[i] * a_exp + 4.0 * g_xx_y_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_yy_xz[i] = -2.0 * g_0_y_yy_xz[i] * a_exp + 4.0 * g_xx_y_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_yy_yy[i] = -2.0 * g_0_y_yy_yy[i] * a_exp + 4.0 * g_xx_y_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_yy_yz[i] = -2.0 * g_0_y_yy_yz[i] * a_exp + 4.0 * g_xx_y_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_yy_zz[i] = -2.0 * g_0_y_yy_zz[i] * a_exp + 4.0 * g_xx_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_0_y_yz_xx, g_0_y_yz_xy, g_0_y_yz_xz, g_0_y_yz_yy, g_0_y_yz_yz, g_0_y_yz_zz, g_xx_0_0_0_0_y_yz_xx, g_xx_0_0_0_0_y_yz_xy, g_xx_0_0_0_0_y_yz_xz, g_xx_0_0_0_0_y_yz_yy, g_xx_0_0_0_0_y_yz_yz, g_xx_0_0_0_0_y_yz_zz, g_xx_y_yz_xx, g_xx_y_yz_xy, g_xx_y_yz_xz, g_xx_y_yz_yy, g_xx_y_yz_yz, g_xx_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_y_yz_xx[i] = -2.0 * g_0_y_yz_xx[i] * a_exp + 4.0 * g_xx_y_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_yz_xy[i] = -2.0 * g_0_y_yz_xy[i] * a_exp + 4.0 * g_xx_y_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_yz_xz[i] = -2.0 * g_0_y_yz_xz[i] * a_exp + 4.0 * g_xx_y_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_yz_yy[i] = -2.0 * g_0_y_yz_yy[i] * a_exp + 4.0 * g_xx_y_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_yz_yz[i] = -2.0 * g_0_y_yz_yz[i] * a_exp + 4.0 * g_xx_y_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_yz_zz[i] = -2.0 * g_0_y_yz_zz[i] * a_exp + 4.0 * g_xx_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_0_y_zz_xx, g_0_y_zz_xy, g_0_y_zz_xz, g_0_y_zz_yy, g_0_y_zz_yz, g_0_y_zz_zz, g_xx_0_0_0_0_y_zz_xx, g_xx_0_0_0_0_y_zz_xy, g_xx_0_0_0_0_y_zz_xz, g_xx_0_0_0_0_y_zz_yy, g_xx_0_0_0_0_y_zz_yz, g_xx_0_0_0_0_y_zz_zz, g_xx_y_zz_xx, g_xx_y_zz_xy, g_xx_y_zz_xz, g_xx_y_zz_yy, g_xx_y_zz_yz, g_xx_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_y_zz_xx[i] = -2.0 * g_0_y_zz_xx[i] * a_exp + 4.0 * g_xx_y_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_zz_xy[i] = -2.0 * g_0_y_zz_xy[i] * a_exp + 4.0 * g_xx_y_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_zz_xz[i] = -2.0 * g_0_y_zz_xz[i] * a_exp + 4.0 * g_xx_y_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_zz_yy[i] = -2.0 * g_0_y_zz_yy[i] * a_exp + 4.0 * g_xx_y_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_zz_yz[i] = -2.0 * g_0_y_zz_yz[i] * a_exp + 4.0 * g_xx_y_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_zz_zz[i] = -2.0 * g_0_y_zz_zz[i] * a_exp + 4.0 * g_xx_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_0_z_xx_xx, g_0_z_xx_xy, g_0_z_xx_xz, g_0_z_xx_yy, g_0_z_xx_yz, g_0_z_xx_zz, g_xx_0_0_0_0_z_xx_xx, g_xx_0_0_0_0_z_xx_xy, g_xx_0_0_0_0_z_xx_xz, g_xx_0_0_0_0_z_xx_yy, g_xx_0_0_0_0_z_xx_yz, g_xx_0_0_0_0_z_xx_zz, g_xx_z_xx_xx, g_xx_z_xx_xy, g_xx_z_xx_xz, g_xx_z_xx_yy, g_xx_z_xx_yz, g_xx_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_z_xx_xx[i] = -2.0 * g_0_z_xx_xx[i] * a_exp + 4.0 * g_xx_z_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_xx_xy[i] = -2.0 * g_0_z_xx_xy[i] * a_exp + 4.0 * g_xx_z_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_xx_xz[i] = -2.0 * g_0_z_xx_xz[i] * a_exp + 4.0 * g_xx_z_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_xx_yy[i] = -2.0 * g_0_z_xx_yy[i] * a_exp + 4.0 * g_xx_z_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_xx_yz[i] = -2.0 * g_0_z_xx_yz[i] * a_exp + 4.0 * g_xx_z_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_xx_zz[i] = -2.0 * g_0_z_xx_zz[i] * a_exp + 4.0 * g_xx_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_0_z_xy_xx, g_0_z_xy_xy, g_0_z_xy_xz, g_0_z_xy_yy, g_0_z_xy_yz, g_0_z_xy_zz, g_xx_0_0_0_0_z_xy_xx, g_xx_0_0_0_0_z_xy_xy, g_xx_0_0_0_0_z_xy_xz, g_xx_0_0_0_0_z_xy_yy, g_xx_0_0_0_0_z_xy_yz, g_xx_0_0_0_0_z_xy_zz, g_xx_z_xy_xx, g_xx_z_xy_xy, g_xx_z_xy_xz, g_xx_z_xy_yy, g_xx_z_xy_yz, g_xx_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_z_xy_xx[i] = -2.0 * g_0_z_xy_xx[i] * a_exp + 4.0 * g_xx_z_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_xy_xy[i] = -2.0 * g_0_z_xy_xy[i] * a_exp + 4.0 * g_xx_z_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_xy_xz[i] = -2.0 * g_0_z_xy_xz[i] * a_exp + 4.0 * g_xx_z_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_xy_yy[i] = -2.0 * g_0_z_xy_yy[i] * a_exp + 4.0 * g_xx_z_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_xy_yz[i] = -2.0 * g_0_z_xy_yz[i] * a_exp + 4.0 * g_xx_z_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_xy_zz[i] = -2.0 * g_0_z_xy_zz[i] * a_exp + 4.0 * g_xx_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_0_z_xz_xx, g_0_z_xz_xy, g_0_z_xz_xz, g_0_z_xz_yy, g_0_z_xz_yz, g_0_z_xz_zz, g_xx_0_0_0_0_z_xz_xx, g_xx_0_0_0_0_z_xz_xy, g_xx_0_0_0_0_z_xz_xz, g_xx_0_0_0_0_z_xz_yy, g_xx_0_0_0_0_z_xz_yz, g_xx_0_0_0_0_z_xz_zz, g_xx_z_xz_xx, g_xx_z_xz_xy, g_xx_z_xz_xz, g_xx_z_xz_yy, g_xx_z_xz_yz, g_xx_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_z_xz_xx[i] = -2.0 * g_0_z_xz_xx[i] * a_exp + 4.0 * g_xx_z_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_xz_xy[i] = -2.0 * g_0_z_xz_xy[i] * a_exp + 4.0 * g_xx_z_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_xz_xz[i] = -2.0 * g_0_z_xz_xz[i] * a_exp + 4.0 * g_xx_z_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_xz_yy[i] = -2.0 * g_0_z_xz_yy[i] * a_exp + 4.0 * g_xx_z_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_xz_yz[i] = -2.0 * g_0_z_xz_yz[i] * a_exp + 4.0 * g_xx_z_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_xz_zz[i] = -2.0 * g_0_z_xz_zz[i] * a_exp + 4.0 * g_xx_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_0_z_yy_xx, g_0_z_yy_xy, g_0_z_yy_xz, g_0_z_yy_yy, g_0_z_yy_yz, g_0_z_yy_zz, g_xx_0_0_0_0_z_yy_xx, g_xx_0_0_0_0_z_yy_xy, g_xx_0_0_0_0_z_yy_xz, g_xx_0_0_0_0_z_yy_yy, g_xx_0_0_0_0_z_yy_yz, g_xx_0_0_0_0_z_yy_zz, g_xx_z_yy_xx, g_xx_z_yy_xy, g_xx_z_yy_xz, g_xx_z_yy_yy, g_xx_z_yy_yz, g_xx_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_z_yy_xx[i] = -2.0 * g_0_z_yy_xx[i] * a_exp + 4.0 * g_xx_z_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_yy_xy[i] = -2.0 * g_0_z_yy_xy[i] * a_exp + 4.0 * g_xx_z_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_yy_xz[i] = -2.0 * g_0_z_yy_xz[i] * a_exp + 4.0 * g_xx_z_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_yy_yy[i] = -2.0 * g_0_z_yy_yy[i] * a_exp + 4.0 * g_xx_z_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_yy_yz[i] = -2.0 * g_0_z_yy_yz[i] * a_exp + 4.0 * g_xx_z_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_yy_zz[i] = -2.0 * g_0_z_yy_zz[i] * a_exp + 4.0 * g_xx_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_0_z_yz_xx, g_0_z_yz_xy, g_0_z_yz_xz, g_0_z_yz_yy, g_0_z_yz_yz, g_0_z_yz_zz, g_xx_0_0_0_0_z_yz_xx, g_xx_0_0_0_0_z_yz_xy, g_xx_0_0_0_0_z_yz_xz, g_xx_0_0_0_0_z_yz_yy, g_xx_0_0_0_0_z_yz_yz, g_xx_0_0_0_0_z_yz_zz, g_xx_z_yz_xx, g_xx_z_yz_xy, g_xx_z_yz_xz, g_xx_z_yz_yy, g_xx_z_yz_yz, g_xx_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_z_yz_xx[i] = -2.0 * g_0_z_yz_xx[i] * a_exp + 4.0 * g_xx_z_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_yz_xy[i] = -2.0 * g_0_z_yz_xy[i] * a_exp + 4.0 * g_xx_z_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_yz_xz[i] = -2.0 * g_0_z_yz_xz[i] * a_exp + 4.0 * g_xx_z_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_yz_yy[i] = -2.0 * g_0_z_yz_yy[i] * a_exp + 4.0 * g_xx_z_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_yz_yz[i] = -2.0 * g_0_z_yz_yz[i] * a_exp + 4.0 * g_xx_z_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_yz_zz[i] = -2.0 * g_0_z_yz_zz[i] * a_exp + 4.0 * g_xx_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_0_z_zz_xx, g_0_z_zz_xy, g_0_z_zz_xz, g_0_z_zz_yy, g_0_z_zz_yz, g_0_z_zz_zz, g_xx_0_0_0_0_z_zz_xx, g_xx_0_0_0_0_z_zz_xy, g_xx_0_0_0_0_z_zz_xz, g_xx_0_0_0_0_z_zz_yy, g_xx_0_0_0_0_z_zz_yz, g_xx_0_0_0_0_z_zz_zz, g_xx_z_zz_xx, g_xx_z_zz_xy, g_xx_z_zz_xz, g_xx_z_zz_yy, g_xx_z_zz_yz, g_xx_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_z_zz_xx[i] = -2.0 * g_0_z_zz_xx[i] * a_exp + 4.0 * g_xx_z_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_zz_xy[i] = -2.0 * g_0_z_zz_xy[i] * a_exp + 4.0 * g_xx_z_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_zz_xz[i] = -2.0 * g_0_z_zz_xz[i] * a_exp + 4.0 * g_xx_z_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_zz_yy[i] = -2.0 * g_0_z_zz_yy[i] * a_exp + 4.0 * g_xx_z_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_zz_yz[i] = -2.0 * g_0_z_zz_yz[i] * a_exp + 4.0 * g_xx_z_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_zz_zz[i] = -2.0 * g_0_z_zz_zz[i] * a_exp + 4.0 * g_xx_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_xy_0_0_0_0_x_xx_xx, g_xy_0_0_0_0_x_xx_xy, g_xy_0_0_0_0_x_xx_xz, g_xy_0_0_0_0_x_xx_yy, g_xy_0_0_0_0_x_xx_yz, g_xy_0_0_0_0_x_xx_zz, g_xy_x_xx_xx, g_xy_x_xx_xy, g_xy_x_xx_xz, g_xy_x_xx_yy, g_xy_x_xx_yz, g_xy_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_x_xx_xx[i] = 4.0 * g_xy_x_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_xx_xy[i] = 4.0 * g_xy_x_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_xx_xz[i] = 4.0 * g_xy_x_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_xx_yy[i] = 4.0 * g_xy_x_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_xx_yz[i] = 4.0 * g_xy_x_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_xx_zz[i] = 4.0 * g_xy_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_xy_0_0_0_0_x_xy_xx, g_xy_0_0_0_0_x_xy_xy, g_xy_0_0_0_0_x_xy_xz, g_xy_0_0_0_0_x_xy_yy, g_xy_0_0_0_0_x_xy_yz, g_xy_0_0_0_0_x_xy_zz, g_xy_x_xy_xx, g_xy_x_xy_xy, g_xy_x_xy_xz, g_xy_x_xy_yy, g_xy_x_xy_yz, g_xy_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_x_xy_xx[i] = 4.0 * g_xy_x_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_xy_xy[i] = 4.0 * g_xy_x_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_xy_xz[i] = 4.0 * g_xy_x_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_xy_yy[i] = 4.0 * g_xy_x_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_xy_yz[i] = 4.0 * g_xy_x_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_xy_zz[i] = 4.0 * g_xy_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_xy_0_0_0_0_x_xz_xx, g_xy_0_0_0_0_x_xz_xy, g_xy_0_0_0_0_x_xz_xz, g_xy_0_0_0_0_x_xz_yy, g_xy_0_0_0_0_x_xz_yz, g_xy_0_0_0_0_x_xz_zz, g_xy_x_xz_xx, g_xy_x_xz_xy, g_xy_x_xz_xz, g_xy_x_xz_yy, g_xy_x_xz_yz, g_xy_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_x_xz_xx[i] = 4.0 * g_xy_x_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_xz_xy[i] = 4.0 * g_xy_x_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_xz_xz[i] = 4.0 * g_xy_x_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_xz_yy[i] = 4.0 * g_xy_x_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_xz_yz[i] = 4.0 * g_xy_x_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_xz_zz[i] = 4.0 * g_xy_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_xy_0_0_0_0_x_yy_xx, g_xy_0_0_0_0_x_yy_xy, g_xy_0_0_0_0_x_yy_xz, g_xy_0_0_0_0_x_yy_yy, g_xy_0_0_0_0_x_yy_yz, g_xy_0_0_0_0_x_yy_zz, g_xy_x_yy_xx, g_xy_x_yy_xy, g_xy_x_yy_xz, g_xy_x_yy_yy, g_xy_x_yy_yz, g_xy_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_x_yy_xx[i] = 4.0 * g_xy_x_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_yy_xy[i] = 4.0 * g_xy_x_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_yy_xz[i] = 4.0 * g_xy_x_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_yy_yy[i] = 4.0 * g_xy_x_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_yy_yz[i] = 4.0 * g_xy_x_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_yy_zz[i] = 4.0 * g_xy_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_xy_0_0_0_0_x_yz_xx, g_xy_0_0_0_0_x_yz_xy, g_xy_0_0_0_0_x_yz_xz, g_xy_0_0_0_0_x_yz_yy, g_xy_0_0_0_0_x_yz_yz, g_xy_0_0_0_0_x_yz_zz, g_xy_x_yz_xx, g_xy_x_yz_xy, g_xy_x_yz_xz, g_xy_x_yz_yy, g_xy_x_yz_yz, g_xy_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_x_yz_xx[i] = 4.0 * g_xy_x_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_yz_xy[i] = 4.0 * g_xy_x_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_yz_xz[i] = 4.0 * g_xy_x_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_yz_yy[i] = 4.0 * g_xy_x_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_yz_yz[i] = 4.0 * g_xy_x_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_yz_zz[i] = 4.0 * g_xy_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_xy_0_0_0_0_x_zz_xx, g_xy_0_0_0_0_x_zz_xy, g_xy_0_0_0_0_x_zz_xz, g_xy_0_0_0_0_x_zz_yy, g_xy_0_0_0_0_x_zz_yz, g_xy_0_0_0_0_x_zz_zz, g_xy_x_zz_xx, g_xy_x_zz_xy, g_xy_x_zz_xz, g_xy_x_zz_yy, g_xy_x_zz_yz, g_xy_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_x_zz_xx[i] = 4.0 * g_xy_x_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_zz_xy[i] = 4.0 * g_xy_x_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_zz_xz[i] = 4.0 * g_xy_x_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_zz_yy[i] = 4.0 * g_xy_x_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_zz_yz[i] = 4.0 * g_xy_x_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_zz_zz[i] = 4.0 * g_xy_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_xy_0_0_0_0_y_xx_xx, g_xy_0_0_0_0_y_xx_xy, g_xy_0_0_0_0_y_xx_xz, g_xy_0_0_0_0_y_xx_yy, g_xy_0_0_0_0_y_xx_yz, g_xy_0_0_0_0_y_xx_zz, g_xy_y_xx_xx, g_xy_y_xx_xy, g_xy_y_xx_xz, g_xy_y_xx_yy, g_xy_y_xx_yz, g_xy_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_y_xx_xx[i] = 4.0 * g_xy_y_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_xx_xy[i] = 4.0 * g_xy_y_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_xx_xz[i] = 4.0 * g_xy_y_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_xx_yy[i] = 4.0 * g_xy_y_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_xx_yz[i] = 4.0 * g_xy_y_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_xx_zz[i] = 4.0 * g_xy_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_xy_0_0_0_0_y_xy_xx, g_xy_0_0_0_0_y_xy_xy, g_xy_0_0_0_0_y_xy_xz, g_xy_0_0_0_0_y_xy_yy, g_xy_0_0_0_0_y_xy_yz, g_xy_0_0_0_0_y_xy_zz, g_xy_y_xy_xx, g_xy_y_xy_xy, g_xy_y_xy_xz, g_xy_y_xy_yy, g_xy_y_xy_yz, g_xy_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_y_xy_xx[i] = 4.0 * g_xy_y_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_xy_xy[i] = 4.0 * g_xy_y_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_xy_xz[i] = 4.0 * g_xy_y_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_xy_yy[i] = 4.0 * g_xy_y_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_xy_yz[i] = 4.0 * g_xy_y_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_xy_zz[i] = 4.0 * g_xy_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_xy_0_0_0_0_y_xz_xx, g_xy_0_0_0_0_y_xz_xy, g_xy_0_0_0_0_y_xz_xz, g_xy_0_0_0_0_y_xz_yy, g_xy_0_0_0_0_y_xz_yz, g_xy_0_0_0_0_y_xz_zz, g_xy_y_xz_xx, g_xy_y_xz_xy, g_xy_y_xz_xz, g_xy_y_xz_yy, g_xy_y_xz_yz, g_xy_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_y_xz_xx[i] = 4.0 * g_xy_y_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_xz_xy[i] = 4.0 * g_xy_y_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_xz_xz[i] = 4.0 * g_xy_y_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_xz_yy[i] = 4.0 * g_xy_y_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_xz_yz[i] = 4.0 * g_xy_y_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_xz_zz[i] = 4.0 * g_xy_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_xy_0_0_0_0_y_yy_xx, g_xy_0_0_0_0_y_yy_xy, g_xy_0_0_0_0_y_yy_xz, g_xy_0_0_0_0_y_yy_yy, g_xy_0_0_0_0_y_yy_yz, g_xy_0_0_0_0_y_yy_zz, g_xy_y_yy_xx, g_xy_y_yy_xy, g_xy_y_yy_xz, g_xy_y_yy_yy, g_xy_y_yy_yz, g_xy_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_y_yy_xx[i] = 4.0 * g_xy_y_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_yy_xy[i] = 4.0 * g_xy_y_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_yy_xz[i] = 4.0 * g_xy_y_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_yy_yy[i] = 4.0 * g_xy_y_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_yy_yz[i] = 4.0 * g_xy_y_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_yy_zz[i] = 4.0 * g_xy_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_xy_0_0_0_0_y_yz_xx, g_xy_0_0_0_0_y_yz_xy, g_xy_0_0_0_0_y_yz_xz, g_xy_0_0_0_0_y_yz_yy, g_xy_0_0_0_0_y_yz_yz, g_xy_0_0_0_0_y_yz_zz, g_xy_y_yz_xx, g_xy_y_yz_xy, g_xy_y_yz_xz, g_xy_y_yz_yy, g_xy_y_yz_yz, g_xy_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_y_yz_xx[i] = 4.0 * g_xy_y_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_yz_xy[i] = 4.0 * g_xy_y_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_yz_xz[i] = 4.0 * g_xy_y_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_yz_yy[i] = 4.0 * g_xy_y_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_yz_yz[i] = 4.0 * g_xy_y_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_yz_zz[i] = 4.0 * g_xy_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_xy_0_0_0_0_y_zz_xx, g_xy_0_0_0_0_y_zz_xy, g_xy_0_0_0_0_y_zz_xz, g_xy_0_0_0_0_y_zz_yy, g_xy_0_0_0_0_y_zz_yz, g_xy_0_0_0_0_y_zz_zz, g_xy_y_zz_xx, g_xy_y_zz_xy, g_xy_y_zz_xz, g_xy_y_zz_yy, g_xy_y_zz_yz, g_xy_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_y_zz_xx[i] = 4.0 * g_xy_y_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_zz_xy[i] = 4.0 * g_xy_y_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_zz_xz[i] = 4.0 * g_xy_y_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_zz_yy[i] = 4.0 * g_xy_y_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_zz_yz[i] = 4.0 * g_xy_y_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_zz_zz[i] = 4.0 * g_xy_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_xy_0_0_0_0_z_xx_xx, g_xy_0_0_0_0_z_xx_xy, g_xy_0_0_0_0_z_xx_xz, g_xy_0_0_0_0_z_xx_yy, g_xy_0_0_0_0_z_xx_yz, g_xy_0_0_0_0_z_xx_zz, g_xy_z_xx_xx, g_xy_z_xx_xy, g_xy_z_xx_xz, g_xy_z_xx_yy, g_xy_z_xx_yz, g_xy_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_z_xx_xx[i] = 4.0 * g_xy_z_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_xx_xy[i] = 4.0 * g_xy_z_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_xx_xz[i] = 4.0 * g_xy_z_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_xx_yy[i] = 4.0 * g_xy_z_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_xx_yz[i] = 4.0 * g_xy_z_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_xx_zz[i] = 4.0 * g_xy_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_xy_0_0_0_0_z_xy_xx, g_xy_0_0_0_0_z_xy_xy, g_xy_0_0_0_0_z_xy_xz, g_xy_0_0_0_0_z_xy_yy, g_xy_0_0_0_0_z_xy_yz, g_xy_0_0_0_0_z_xy_zz, g_xy_z_xy_xx, g_xy_z_xy_xy, g_xy_z_xy_xz, g_xy_z_xy_yy, g_xy_z_xy_yz, g_xy_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_z_xy_xx[i] = 4.0 * g_xy_z_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_xy_xy[i] = 4.0 * g_xy_z_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_xy_xz[i] = 4.0 * g_xy_z_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_xy_yy[i] = 4.0 * g_xy_z_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_xy_yz[i] = 4.0 * g_xy_z_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_xy_zz[i] = 4.0 * g_xy_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_xy_0_0_0_0_z_xz_xx, g_xy_0_0_0_0_z_xz_xy, g_xy_0_0_0_0_z_xz_xz, g_xy_0_0_0_0_z_xz_yy, g_xy_0_0_0_0_z_xz_yz, g_xy_0_0_0_0_z_xz_zz, g_xy_z_xz_xx, g_xy_z_xz_xy, g_xy_z_xz_xz, g_xy_z_xz_yy, g_xy_z_xz_yz, g_xy_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_z_xz_xx[i] = 4.0 * g_xy_z_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_xz_xy[i] = 4.0 * g_xy_z_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_xz_xz[i] = 4.0 * g_xy_z_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_xz_yy[i] = 4.0 * g_xy_z_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_xz_yz[i] = 4.0 * g_xy_z_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_xz_zz[i] = 4.0 * g_xy_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_xy_0_0_0_0_z_yy_xx, g_xy_0_0_0_0_z_yy_xy, g_xy_0_0_0_0_z_yy_xz, g_xy_0_0_0_0_z_yy_yy, g_xy_0_0_0_0_z_yy_yz, g_xy_0_0_0_0_z_yy_zz, g_xy_z_yy_xx, g_xy_z_yy_xy, g_xy_z_yy_xz, g_xy_z_yy_yy, g_xy_z_yy_yz, g_xy_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_z_yy_xx[i] = 4.0 * g_xy_z_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_yy_xy[i] = 4.0 * g_xy_z_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_yy_xz[i] = 4.0 * g_xy_z_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_yy_yy[i] = 4.0 * g_xy_z_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_yy_yz[i] = 4.0 * g_xy_z_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_yy_zz[i] = 4.0 * g_xy_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_xy_0_0_0_0_z_yz_xx, g_xy_0_0_0_0_z_yz_xy, g_xy_0_0_0_0_z_yz_xz, g_xy_0_0_0_0_z_yz_yy, g_xy_0_0_0_0_z_yz_yz, g_xy_0_0_0_0_z_yz_zz, g_xy_z_yz_xx, g_xy_z_yz_xy, g_xy_z_yz_xz, g_xy_z_yz_yy, g_xy_z_yz_yz, g_xy_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_z_yz_xx[i] = 4.0 * g_xy_z_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_yz_xy[i] = 4.0 * g_xy_z_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_yz_xz[i] = 4.0 * g_xy_z_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_yz_yy[i] = 4.0 * g_xy_z_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_yz_yz[i] = 4.0 * g_xy_z_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_yz_zz[i] = 4.0 * g_xy_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_xy_0_0_0_0_z_zz_xx, g_xy_0_0_0_0_z_zz_xy, g_xy_0_0_0_0_z_zz_xz, g_xy_0_0_0_0_z_zz_yy, g_xy_0_0_0_0_z_zz_yz, g_xy_0_0_0_0_z_zz_zz, g_xy_z_zz_xx, g_xy_z_zz_xy, g_xy_z_zz_xz, g_xy_z_zz_yy, g_xy_z_zz_yz, g_xy_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_z_zz_xx[i] = 4.0 * g_xy_z_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_zz_xy[i] = 4.0 * g_xy_z_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_zz_xz[i] = 4.0 * g_xy_z_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_zz_yy[i] = 4.0 * g_xy_z_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_zz_yz[i] = 4.0 * g_xy_z_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_zz_zz[i] = 4.0 * g_xy_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_xz_0_0_0_0_x_xx_xx, g_xz_0_0_0_0_x_xx_xy, g_xz_0_0_0_0_x_xx_xz, g_xz_0_0_0_0_x_xx_yy, g_xz_0_0_0_0_x_xx_yz, g_xz_0_0_0_0_x_xx_zz, g_xz_x_xx_xx, g_xz_x_xx_xy, g_xz_x_xx_xz, g_xz_x_xx_yy, g_xz_x_xx_yz, g_xz_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_x_xx_xx[i] = 4.0 * g_xz_x_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_xx_xy[i] = 4.0 * g_xz_x_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_xx_xz[i] = 4.0 * g_xz_x_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_xx_yy[i] = 4.0 * g_xz_x_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_xx_yz[i] = 4.0 * g_xz_x_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_xx_zz[i] = 4.0 * g_xz_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_xz_0_0_0_0_x_xy_xx, g_xz_0_0_0_0_x_xy_xy, g_xz_0_0_0_0_x_xy_xz, g_xz_0_0_0_0_x_xy_yy, g_xz_0_0_0_0_x_xy_yz, g_xz_0_0_0_0_x_xy_zz, g_xz_x_xy_xx, g_xz_x_xy_xy, g_xz_x_xy_xz, g_xz_x_xy_yy, g_xz_x_xy_yz, g_xz_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_x_xy_xx[i] = 4.0 * g_xz_x_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_xy_xy[i] = 4.0 * g_xz_x_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_xy_xz[i] = 4.0 * g_xz_x_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_xy_yy[i] = 4.0 * g_xz_x_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_xy_yz[i] = 4.0 * g_xz_x_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_xy_zz[i] = 4.0 * g_xz_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_xz_0_0_0_0_x_xz_xx, g_xz_0_0_0_0_x_xz_xy, g_xz_0_0_0_0_x_xz_xz, g_xz_0_0_0_0_x_xz_yy, g_xz_0_0_0_0_x_xz_yz, g_xz_0_0_0_0_x_xz_zz, g_xz_x_xz_xx, g_xz_x_xz_xy, g_xz_x_xz_xz, g_xz_x_xz_yy, g_xz_x_xz_yz, g_xz_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_x_xz_xx[i] = 4.0 * g_xz_x_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_xz_xy[i] = 4.0 * g_xz_x_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_xz_xz[i] = 4.0 * g_xz_x_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_xz_yy[i] = 4.0 * g_xz_x_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_xz_yz[i] = 4.0 * g_xz_x_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_xz_zz[i] = 4.0 * g_xz_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_xz_0_0_0_0_x_yy_xx, g_xz_0_0_0_0_x_yy_xy, g_xz_0_0_0_0_x_yy_xz, g_xz_0_0_0_0_x_yy_yy, g_xz_0_0_0_0_x_yy_yz, g_xz_0_0_0_0_x_yy_zz, g_xz_x_yy_xx, g_xz_x_yy_xy, g_xz_x_yy_xz, g_xz_x_yy_yy, g_xz_x_yy_yz, g_xz_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_x_yy_xx[i] = 4.0 * g_xz_x_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_yy_xy[i] = 4.0 * g_xz_x_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_yy_xz[i] = 4.0 * g_xz_x_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_yy_yy[i] = 4.0 * g_xz_x_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_yy_yz[i] = 4.0 * g_xz_x_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_yy_zz[i] = 4.0 * g_xz_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_xz_0_0_0_0_x_yz_xx, g_xz_0_0_0_0_x_yz_xy, g_xz_0_0_0_0_x_yz_xz, g_xz_0_0_0_0_x_yz_yy, g_xz_0_0_0_0_x_yz_yz, g_xz_0_0_0_0_x_yz_zz, g_xz_x_yz_xx, g_xz_x_yz_xy, g_xz_x_yz_xz, g_xz_x_yz_yy, g_xz_x_yz_yz, g_xz_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_x_yz_xx[i] = 4.0 * g_xz_x_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_yz_xy[i] = 4.0 * g_xz_x_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_yz_xz[i] = 4.0 * g_xz_x_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_yz_yy[i] = 4.0 * g_xz_x_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_yz_yz[i] = 4.0 * g_xz_x_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_yz_zz[i] = 4.0 * g_xz_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_xz_0_0_0_0_x_zz_xx, g_xz_0_0_0_0_x_zz_xy, g_xz_0_0_0_0_x_zz_xz, g_xz_0_0_0_0_x_zz_yy, g_xz_0_0_0_0_x_zz_yz, g_xz_0_0_0_0_x_zz_zz, g_xz_x_zz_xx, g_xz_x_zz_xy, g_xz_x_zz_xz, g_xz_x_zz_yy, g_xz_x_zz_yz, g_xz_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_x_zz_xx[i] = 4.0 * g_xz_x_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_zz_xy[i] = 4.0 * g_xz_x_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_zz_xz[i] = 4.0 * g_xz_x_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_zz_yy[i] = 4.0 * g_xz_x_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_zz_yz[i] = 4.0 * g_xz_x_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_zz_zz[i] = 4.0 * g_xz_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_xz_0_0_0_0_y_xx_xx, g_xz_0_0_0_0_y_xx_xy, g_xz_0_0_0_0_y_xx_xz, g_xz_0_0_0_0_y_xx_yy, g_xz_0_0_0_0_y_xx_yz, g_xz_0_0_0_0_y_xx_zz, g_xz_y_xx_xx, g_xz_y_xx_xy, g_xz_y_xx_xz, g_xz_y_xx_yy, g_xz_y_xx_yz, g_xz_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_y_xx_xx[i] = 4.0 * g_xz_y_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_xx_xy[i] = 4.0 * g_xz_y_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_xx_xz[i] = 4.0 * g_xz_y_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_xx_yy[i] = 4.0 * g_xz_y_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_xx_yz[i] = 4.0 * g_xz_y_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_xx_zz[i] = 4.0 * g_xz_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_xz_0_0_0_0_y_xy_xx, g_xz_0_0_0_0_y_xy_xy, g_xz_0_0_0_0_y_xy_xz, g_xz_0_0_0_0_y_xy_yy, g_xz_0_0_0_0_y_xy_yz, g_xz_0_0_0_0_y_xy_zz, g_xz_y_xy_xx, g_xz_y_xy_xy, g_xz_y_xy_xz, g_xz_y_xy_yy, g_xz_y_xy_yz, g_xz_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_y_xy_xx[i] = 4.0 * g_xz_y_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_xy_xy[i] = 4.0 * g_xz_y_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_xy_xz[i] = 4.0 * g_xz_y_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_xy_yy[i] = 4.0 * g_xz_y_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_xy_yz[i] = 4.0 * g_xz_y_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_xy_zz[i] = 4.0 * g_xz_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_xz_0_0_0_0_y_xz_xx, g_xz_0_0_0_0_y_xz_xy, g_xz_0_0_0_0_y_xz_xz, g_xz_0_0_0_0_y_xz_yy, g_xz_0_0_0_0_y_xz_yz, g_xz_0_0_0_0_y_xz_zz, g_xz_y_xz_xx, g_xz_y_xz_xy, g_xz_y_xz_xz, g_xz_y_xz_yy, g_xz_y_xz_yz, g_xz_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_y_xz_xx[i] = 4.0 * g_xz_y_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_xz_xy[i] = 4.0 * g_xz_y_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_xz_xz[i] = 4.0 * g_xz_y_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_xz_yy[i] = 4.0 * g_xz_y_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_xz_yz[i] = 4.0 * g_xz_y_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_xz_zz[i] = 4.0 * g_xz_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_xz_0_0_0_0_y_yy_xx, g_xz_0_0_0_0_y_yy_xy, g_xz_0_0_0_0_y_yy_xz, g_xz_0_0_0_0_y_yy_yy, g_xz_0_0_0_0_y_yy_yz, g_xz_0_0_0_0_y_yy_zz, g_xz_y_yy_xx, g_xz_y_yy_xy, g_xz_y_yy_xz, g_xz_y_yy_yy, g_xz_y_yy_yz, g_xz_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_y_yy_xx[i] = 4.0 * g_xz_y_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_yy_xy[i] = 4.0 * g_xz_y_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_yy_xz[i] = 4.0 * g_xz_y_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_yy_yy[i] = 4.0 * g_xz_y_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_yy_yz[i] = 4.0 * g_xz_y_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_yy_zz[i] = 4.0 * g_xz_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_xz_0_0_0_0_y_yz_xx, g_xz_0_0_0_0_y_yz_xy, g_xz_0_0_0_0_y_yz_xz, g_xz_0_0_0_0_y_yz_yy, g_xz_0_0_0_0_y_yz_yz, g_xz_0_0_0_0_y_yz_zz, g_xz_y_yz_xx, g_xz_y_yz_xy, g_xz_y_yz_xz, g_xz_y_yz_yy, g_xz_y_yz_yz, g_xz_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_y_yz_xx[i] = 4.0 * g_xz_y_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_yz_xy[i] = 4.0 * g_xz_y_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_yz_xz[i] = 4.0 * g_xz_y_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_yz_yy[i] = 4.0 * g_xz_y_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_yz_yz[i] = 4.0 * g_xz_y_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_yz_zz[i] = 4.0 * g_xz_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_xz_0_0_0_0_y_zz_xx, g_xz_0_0_0_0_y_zz_xy, g_xz_0_0_0_0_y_zz_xz, g_xz_0_0_0_0_y_zz_yy, g_xz_0_0_0_0_y_zz_yz, g_xz_0_0_0_0_y_zz_zz, g_xz_y_zz_xx, g_xz_y_zz_xy, g_xz_y_zz_xz, g_xz_y_zz_yy, g_xz_y_zz_yz, g_xz_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_y_zz_xx[i] = 4.0 * g_xz_y_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_zz_xy[i] = 4.0 * g_xz_y_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_zz_xz[i] = 4.0 * g_xz_y_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_zz_yy[i] = 4.0 * g_xz_y_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_zz_yz[i] = 4.0 * g_xz_y_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_zz_zz[i] = 4.0 * g_xz_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_xz_0_0_0_0_z_xx_xx, g_xz_0_0_0_0_z_xx_xy, g_xz_0_0_0_0_z_xx_xz, g_xz_0_0_0_0_z_xx_yy, g_xz_0_0_0_0_z_xx_yz, g_xz_0_0_0_0_z_xx_zz, g_xz_z_xx_xx, g_xz_z_xx_xy, g_xz_z_xx_xz, g_xz_z_xx_yy, g_xz_z_xx_yz, g_xz_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_z_xx_xx[i] = 4.0 * g_xz_z_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_xx_xy[i] = 4.0 * g_xz_z_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_xx_xz[i] = 4.0 * g_xz_z_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_xx_yy[i] = 4.0 * g_xz_z_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_xx_yz[i] = 4.0 * g_xz_z_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_xx_zz[i] = 4.0 * g_xz_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_xz_0_0_0_0_z_xy_xx, g_xz_0_0_0_0_z_xy_xy, g_xz_0_0_0_0_z_xy_xz, g_xz_0_0_0_0_z_xy_yy, g_xz_0_0_0_0_z_xy_yz, g_xz_0_0_0_0_z_xy_zz, g_xz_z_xy_xx, g_xz_z_xy_xy, g_xz_z_xy_xz, g_xz_z_xy_yy, g_xz_z_xy_yz, g_xz_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_z_xy_xx[i] = 4.0 * g_xz_z_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_xy_xy[i] = 4.0 * g_xz_z_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_xy_xz[i] = 4.0 * g_xz_z_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_xy_yy[i] = 4.0 * g_xz_z_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_xy_yz[i] = 4.0 * g_xz_z_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_xy_zz[i] = 4.0 * g_xz_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_xz_0_0_0_0_z_xz_xx, g_xz_0_0_0_0_z_xz_xy, g_xz_0_0_0_0_z_xz_xz, g_xz_0_0_0_0_z_xz_yy, g_xz_0_0_0_0_z_xz_yz, g_xz_0_0_0_0_z_xz_zz, g_xz_z_xz_xx, g_xz_z_xz_xy, g_xz_z_xz_xz, g_xz_z_xz_yy, g_xz_z_xz_yz, g_xz_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_z_xz_xx[i] = 4.0 * g_xz_z_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_xz_xy[i] = 4.0 * g_xz_z_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_xz_xz[i] = 4.0 * g_xz_z_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_xz_yy[i] = 4.0 * g_xz_z_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_xz_yz[i] = 4.0 * g_xz_z_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_xz_zz[i] = 4.0 * g_xz_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_xz_0_0_0_0_z_yy_xx, g_xz_0_0_0_0_z_yy_xy, g_xz_0_0_0_0_z_yy_xz, g_xz_0_0_0_0_z_yy_yy, g_xz_0_0_0_0_z_yy_yz, g_xz_0_0_0_0_z_yy_zz, g_xz_z_yy_xx, g_xz_z_yy_xy, g_xz_z_yy_xz, g_xz_z_yy_yy, g_xz_z_yy_yz, g_xz_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_z_yy_xx[i] = 4.0 * g_xz_z_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_yy_xy[i] = 4.0 * g_xz_z_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_yy_xz[i] = 4.0 * g_xz_z_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_yy_yy[i] = 4.0 * g_xz_z_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_yy_yz[i] = 4.0 * g_xz_z_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_yy_zz[i] = 4.0 * g_xz_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_xz_0_0_0_0_z_yz_xx, g_xz_0_0_0_0_z_yz_xy, g_xz_0_0_0_0_z_yz_xz, g_xz_0_0_0_0_z_yz_yy, g_xz_0_0_0_0_z_yz_yz, g_xz_0_0_0_0_z_yz_zz, g_xz_z_yz_xx, g_xz_z_yz_xy, g_xz_z_yz_xz, g_xz_z_yz_yy, g_xz_z_yz_yz, g_xz_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_z_yz_xx[i] = 4.0 * g_xz_z_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_yz_xy[i] = 4.0 * g_xz_z_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_yz_xz[i] = 4.0 * g_xz_z_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_yz_yy[i] = 4.0 * g_xz_z_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_yz_yz[i] = 4.0 * g_xz_z_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_yz_zz[i] = 4.0 * g_xz_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_xz_0_0_0_0_z_zz_xx, g_xz_0_0_0_0_z_zz_xy, g_xz_0_0_0_0_z_zz_xz, g_xz_0_0_0_0_z_zz_yy, g_xz_0_0_0_0_z_zz_yz, g_xz_0_0_0_0_z_zz_zz, g_xz_z_zz_xx, g_xz_z_zz_xy, g_xz_z_zz_xz, g_xz_z_zz_yy, g_xz_z_zz_yz, g_xz_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_z_zz_xx[i] = 4.0 * g_xz_z_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_zz_xy[i] = 4.0 * g_xz_z_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_zz_xz[i] = 4.0 * g_xz_z_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_zz_yy[i] = 4.0 * g_xz_z_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_zz_yz[i] = 4.0 * g_xz_z_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_zz_zz[i] = 4.0 * g_xz_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_0_x_xx_xx, g_0_x_xx_xy, g_0_x_xx_xz, g_0_x_xx_yy, g_0_x_xx_yz, g_0_x_xx_zz, g_yy_0_0_0_0_x_xx_xx, g_yy_0_0_0_0_x_xx_xy, g_yy_0_0_0_0_x_xx_xz, g_yy_0_0_0_0_x_xx_yy, g_yy_0_0_0_0_x_xx_yz, g_yy_0_0_0_0_x_xx_zz, g_yy_x_xx_xx, g_yy_x_xx_xy, g_yy_x_xx_xz, g_yy_x_xx_yy, g_yy_x_xx_yz, g_yy_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_x_xx_xx[i] = -2.0 * g_0_x_xx_xx[i] * a_exp + 4.0 * g_yy_x_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_xx_xy[i] = -2.0 * g_0_x_xx_xy[i] * a_exp + 4.0 * g_yy_x_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_xx_xz[i] = -2.0 * g_0_x_xx_xz[i] * a_exp + 4.0 * g_yy_x_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_xx_yy[i] = -2.0 * g_0_x_xx_yy[i] * a_exp + 4.0 * g_yy_x_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_xx_yz[i] = -2.0 * g_0_x_xx_yz[i] * a_exp + 4.0 * g_yy_x_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_xx_zz[i] = -2.0 * g_0_x_xx_zz[i] * a_exp + 4.0 * g_yy_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_0_x_xy_xx, g_0_x_xy_xy, g_0_x_xy_xz, g_0_x_xy_yy, g_0_x_xy_yz, g_0_x_xy_zz, g_yy_0_0_0_0_x_xy_xx, g_yy_0_0_0_0_x_xy_xy, g_yy_0_0_0_0_x_xy_xz, g_yy_0_0_0_0_x_xy_yy, g_yy_0_0_0_0_x_xy_yz, g_yy_0_0_0_0_x_xy_zz, g_yy_x_xy_xx, g_yy_x_xy_xy, g_yy_x_xy_xz, g_yy_x_xy_yy, g_yy_x_xy_yz, g_yy_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_x_xy_xx[i] = -2.0 * g_0_x_xy_xx[i] * a_exp + 4.0 * g_yy_x_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_xy_xy[i] = -2.0 * g_0_x_xy_xy[i] * a_exp + 4.0 * g_yy_x_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_xy_xz[i] = -2.0 * g_0_x_xy_xz[i] * a_exp + 4.0 * g_yy_x_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_xy_yy[i] = -2.0 * g_0_x_xy_yy[i] * a_exp + 4.0 * g_yy_x_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_xy_yz[i] = -2.0 * g_0_x_xy_yz[i] * a_exp + 4.0 * g_yy_x_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_xy_zz[i] = -2.0 * g_0_x_xy_zz[i] * a_exp + 4.0 * g_yy_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_0_x_xz_xx, g_0_x_xz_xy, g_0_x_xz_xz, g_0_x_xz_yy, g_0_x_xz_yz, g_0_x_xz_zz, g_yy_0_0_0_0_x_xz_xx, g_yy_0_0_0_0_x_xz_xy, g_yy_0_0_0_0_x_xz_xz, g_yy_0_0_0_0_x_xz_yy, g_yy_0_0_0_0_x_xz_yz, g_yy_0_0_0_0_x_xz_zz, g_yy_x_xz_xx, g_yy_x_xz_xy, g_yy_x_xz_xz, g_yy_x_xz_yy, g_yy_x_xz_yz, g_yy_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_x_xz_xx[i] = -2.0 * g_0_x_xz_xx[i] * a_exp + 4.0 * g_yy_x_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_xz_xy[i] = -2.0 * g_0_x_xz_xy[i] * a_exp + 4.0 * g_yy_x_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_xz_xz[i] = -2.0 * g_0_x_xz_xz[i] * a_exp + 4.0 * g_yy_x_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_xz_yy[i] = -2.0 * g_0_x_xz_yy[i] * a_exp + 4.0 * g_yy_x_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_xz_yz[i] = -2.0 * g_0_x_xz_yz[i] * a_exp + 4.0 * g_yy_x_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_xz_zz[i] = -2.0 * g_0_x_xz_zz[i] * a_exp + 4.0 * g_yy_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_0_x_yy_xx, g_0_x_yy_xy, g_0_x_yy_xz, g_0_x_yy_yy, g_0_x_yy_yz, g_0_x_yy_zz, g_yy_0_0_0_0_x_yy_xx, g_yy_0_0_0_0_x_yy_xy, g_yy_0_0_0_0_x_yy_xz, g_yy_0_0_0_0_x_yy_yy, g_yy_0_0_0_0_x_yy_yz, g_yy_0_0_0_0_x_yy_zz, g_yy_x_yy_xx, g_yy_x_yy_xy, g_yy_x_yy_xz, g_yy_x_yy_yy, g_yy_x_yy_yz, g_yy_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_x_yy_xx[i] = -2.0 * g_0_x_yy_xx[i] * a_exp + 4.0 * g_yy_x_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_yy_xy[i] = -2.0 * g_0_x_yy_xy[i] * a_exp + 4.0 * g_yy_x_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_yy_xz[i] = -2.0 * g_0_x_yy_xz[i] * a_exp + 4.0 * g_yy_x_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_yy_yy[i] = -2.0 * g_0_x_yy_yy[i] * a_exp + 4.0 * g_yy_x_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_yy_yz[i] = -2.0 * g_0_x_yy_yz[i] * a_exp + 4.0 * g_yy_x_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_yy_zz[i] = -2.0 * g_0_x_yy_zz[i] * a_exp + 4.0 * g_yy_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_0_x_yz_xx, g_0_x_yz_xy, g_0_x_yz_xz, g_0_x_yz_yy, g_0_x_yz_yz, g_0_x_yz_zz, g_yy_0_0_0_0_x_yz_xx, g_yy_0_0_0_0_x_yz_xy, g_yy_0_0_0_0_x_yz_xz, g_yy_0_0_0_0_x_yz_yy, g_yy_0_0_0_0_x_yz_yz, g_yy_0_0_0_0_x_yz_zz, g_yy_x_yz_xx, g_yy_x_yz_xy, g_yy_x_yz_xz, g_yy_x_yz_yy, g_yy_x_yz_yz, g_yy_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_x_yz_xx[i] = -2.0 * g_0_x_yz_xx[i] * a_exp + 4.0 * g_yy_x_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_yz_xy[i] = -2.0 * g_0_x_yz_xy[i] * a_exp + 4.0 * g_yy_x_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_yz_xz[i] = -2.0 * g_0_x_yz_xz[i] * a_exp + 4.0 * g_yy_x_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_yz_yy[i] = -2.0 * g_0_x_yz_yy[i] * a_exp + 4.0 * g_yy_x_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_yz_yz[i] = -2.0 * g_0_x_yz_yz[i] * a_exp + 4.0 * g_yy_x_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_yz_zz[i] = -2.0 * g_0_x_yz_zz[i] * a_exp + 4.0 * g_yy_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_0_x_zz_xx, g_0_x_zz_xy, g_0_x_zz_xz, g_0_x_zz_yy, g_0_x_zz_yz, g_0_x_zz_zz, g_yy_0_0_0_0_x_zz_xx, g_yy_0_0_0_0_x_zz_xy, g_yy_0_0_0_0_x_zz_xz, g_yy_0_0_0_0_x_zz_yy, g_yy_0_0_0_0_x_zz_yz, g_yy_0_0_0_0_x_zz_zz, g_yy_x_zz_xx, g_yy_x_zz_xy, g_yy_x_zz_xz, g_yy_x_zz_yy, g_yy_x_zz_yz, g_yy_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_x_zz_xx[i] = -2.0 * g_0_x_zz_xx[i] * a_exp + 4.0 * g_yy_x_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_zz_xy[i] = -2.0 * g_0_x_zz_xy[i] * a_exp + 4.0 * g_yy_x_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_zz_xz[i] = -2.0 * g_0_x_zz_xz[i] * a_exp + 4.0 * g_yy_x_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_zz_yy[i] = -2.0 * g_0_x_zz_yy[i] * a_exp + 4.0 * g_yy_x_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_zz_yz[i] = -2.0 * g_0_x_zz_yz[i] * a_exp + 4.0 * g_yy_x_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_zz_zz[i] = -2.0 * g_0_x_zz_zz[i] * a_exp + 4.0 * g_yy_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_0_y_xx_xx, g_0_y_xx_xy, g_0_y_xx_xz, g_0_y_xx_yy, g_0_y_xx_yz, g_0_y_xx_zz, g_yy_0_0_0_0_y_xx_xx, g_yy_0_0_0_0_y_xx_xy, g_yy_0_0_0_0_y_xx_xz, g_yy_0_0_0_0_y_xx_yy, g_yy_0_0_0_0_y_xx_yz, g_yy_0_0_0_0_y_xx_zz, g_yy_y_xx_xx, g_yy_y_xx_xy, g_yy_y_xx_xz, g_yy_y_xx_yy, g_yy_y_xx_yz, g_yy_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_y_xx_xx[i] = -2.0 * g_0_y_xx_xx[i] * a_exp + 4.0 * g_yy_y_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_xx_xy[i] = -2.0 * g_0_y_xx_xy[i] * a_exp + 4.0 * g_yy_y_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_xx_xz[i] = -2.0 * g_0_y_xx_xz[i] * a_exp + 4.0 * g_yy_y_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_xx_yy[i] = -2.0 * g_0_y_xx_yy[i] * a_exp + 4.0 * g_yy_y_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_xx_yz[i] = -2.0 * g_0_y_xx_yz[i] * a_exp + 4.0 * g_yy_y_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_xx_zz[i] = -2.0 * g_0_y_xx_zz[i] * a_exp + 4.0 * g_yy_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_0_y_xy_xx, g_0_y_xy_xy, g_0_y_xy_xz, g_0_y_xy_yy, g_0_y_xy_yz, g_0_y_xy_zz, g_yy_0_0_0_0_y_xy_xx, g_yy_0_0_0_0_y_xy_xy, g_yy_0_0_0_0_y_xy_xz, g_yy_0_0_0_0_y_xy_yy, g_yy_0_0_0_0_y_xy_yz, g_yy_0_0_0_0_y_xy_zz, g_yy_y_xy_xx, g_yy_y_xy_xy, g_yy_y_xy_xz, g_yy_y_xy_yy, g_yy_y_xy_yz, g_yy_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_y_xy_xx[i] = -2.0 * g_0_y_xy_xx[i] * a_exp + 4.0 * g_yy_y_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_xy_xy[i] = -2.0 * g_0_y_xy_xy[i] * a_exp + 4.0 * g_yy_y_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_xy_xz[i] = -2.0 * g_0_y_xy_xz[i] * a_exp + 4.0 * g_yy_y_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_xy_yy[i] = -2.0 * g_0_y_xy_yy[i] * a_exp + 4.0 * g_yy_y_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_xy_yz[i] = -2.0 * g_0_y_xy_yz[i] * a_exp + 4.0 * g_yy_y_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_xy_zz[i] = -2.0 * g_0_y_xy_zz[i] * a_exp + 4.0 * g_yy_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_0_y_xz_xx, g_0_y_xz_xy, g_0_y_xz_xz, g_0_y_xz_yy, g_0_y_xz_yz, g_0_y_xz_zz, g_yy_0_0_0_0_y_xz_xx, g_yy_0_0_0_0_y_xz_xy, g_yy_0_0_0_0_y_xz_xz, g_yy_0_0_0_0_y_xz_yy, g_yy_0_0_0_0_y_xz_yz, g_yy_0_0_0_0_y_xz_zz, g_yy_y_xz_xx, g_yy_y_xz_xy, g_yy_y_xz_xz, g_yy_y_xz_yy, g_yy_y_xz_yz, g_yy_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_y_xz_xx[i] = -2.0 * g_0_y_xz_xx[i] * a_exp + 4.0 * g_yy_y_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_xz_xy[i] = -2.0 * g_0_y_xz_xy[i] * a_exp + 4.0 * g_yy_y_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_xz_xz[i] = -2.0 * g_0_y_xz_xz[i] * a_exp + 4.0 * g_yy_y_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_xz_yy[i] = -2.0 * g_0_y_xz_yy[i] * a_exp + 4.0 * g_yy_y_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_xz_yz[i] = -2.0 * g_0_y_xz_yz[i] * a_exp + 4.0 * g_yy_y_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_xz_zz[i] = -2.0 * g_0_y_xz_zz[i] * a_exp + 4.0 * g_yy_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_0_y_yy_xx, g_0_y_yy_xy, g_0_y_yy_xz, g_0_y_yy_yy, g_0_y_yy_yz, g_0_y_yy_zz, g_yy_0_0_0_0_y_yy_xx, g_yy_0_0_0_0_y_yy_xy, g_yy_0_0_0_0_y_yy_xz, g_yy_0_0_0_0_y_yy_yy, g_yy_0_0_0_0_y_yy_yz, g_yy_0_0_0_0_y_yy_zz, g_yy_y_yy_xx, g_yy_y_yy_xy, g_yy_y_yy_xz, g_yy_y_yy_yy, g_yy_y_yy_yz, g_yy_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_y_yy_xx[i] = -2.0 * g_0_y_yy_xx[i] * a_exp + 4.0 * g_yy_y_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_yy_xy[i] = -2.0 * g_0_y_yy_xy[i] * a_exp + 4.0 * g_yy_y_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_yy_xz[i] = -2.0 * g_0_y_yy_xz[i] * a_exp + 4.0 * g_yy_y_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_yy_yy[i] = -2.0 * g_0_y_yy_yy[i] * a_exp + 4.0 * g_yy_y_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_yy_yz[i] = -2.0 * g_0_y_yy_yz[i] * a_exp + 4.0 * g_yy_y_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_yy_zz[i] = -2.0 * g_0_y_yy_zz[i] * a_exp + 4.0 * g_yy_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_0_y_yz_xx, g_0_y_yz_xy, g_0_y_yz_xz, g_0_y_yz_yy, g_0_y_yz_yz, g_0_y_yz_zz, g_yy_0_0_0_0_y_yz_xx, g_yy_0_0_0_0_y_yz_xy, g_yy_0_0_0_0_y_yz_xz, g_yy_0_0_0_0_y_yz_yy, g_yy_0_0_0_0_y_yz_yz, g_yy_0_0_0_0_y_yz_zz, g_yy_y_yz_xx, g_yy_y_yz_xy, g_yy_y_yz_xz, g_yy_y_yz_yy, g_yy_y_yz_yz, g_yy_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_y_yz_xx[i] = -2.0 * g_0_y_yz_xx[i] * a_exp + 4.0 * g_yy_y_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_yz_xy[i] = -2.0 * g_0_y_yz_xy[i] * a_exp + 4.0 * g_yy_y_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_yz_xz[i] = -2.0 * g_0_y_yz_xz[i] * a_exp + 4.0 * g_yy_y_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_yz_yy[i] = -2.0 * g_0_y_yz_yy[i] * a_exp + 4.0 * g_yy_y_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_yz_yz[i] = -2.0 * g_0_y_yz_yz[i] * a_exp + 4.0 * g_yy_y_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_yz_zz[i] = -2.0 * g_0_y_yz_zz[i] * a_exp + 4.0 * g_yy_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_0_y_zz_xx, g_0_y_zz_xy, g_0_y_zz_xz, g_0_y_zz_yy, g_0_y_zz_yz, g_0_y_zz_zz, g_yy_0_0_0_0_y_zz_xx, g_yy_0_0_0_0_y_zz_xy, g_yy_0_0_0_0_y_zz_xz, g_yy_0_0_0_0_y_zz_yy, g_yy_0_0_0_0_y_zz_yz, g_yy_0_0_0_0_y_zz_zz, g_yy_y_zz_xx, g_yy_y_zz_xy, g_yy_y_zz_xz, g_yy_y_zz_yy, g_yy_y_zz_yz, g_yy_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_y_zz_xx[i] = -2.0 * g_0_y_zz_xx[i] * a_exp + 4.0 * g_yy_y_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_zz_xy[i] = -2.0 * g_0_y_zz_xy[i] * a_exp + 4.0 * g_yy_y_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_zz_xz[i] = -2.0 * g_0_y_zz_xz[i] * a_exp + 4.0 * g_yy_y_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_zz_yy[i] = -2.0 * g_0_y_zz_yy[i] * a_exp + 4.0 * g_yy_y_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_zz_yz[i] = -2.0 * g_0_y_zz_yz[i] * a_exp + 4.0 * g_yy_y_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_zz_zz[i] = -2.0 * g_0_y_zz_zz[i] * a_exp + 4.0 * g_yy_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_0_z_xx_xx, g_0_z_xx_xy, g_0_z_xx_xz, g_0_z_xx_yy, g_0_z_xx_yz, g_0_z_xx_zz, g_yy_0_0_0_0_z_xx_xx, g_yy_0_0_0_0_z_xx_xy, g_yy_0_0_0_0_z_xx_xz, g_yy_0_0_0_0_z_xx_yy, g_yy_0_0_0_0_z_xx_yz, g_yy_0_0_0_0_z_xx_zz, g_yy_z_xx_xx, g_yy_z_xx_xy, g_yy_z_xx_xz, g_yy_z_xx_yy, g_yy_z_xx_yz, g_yy_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_z_xx_xx[i] = -2.0 * g_0_z_xx_xx[i] * a_exp + 4.0 * g_yy_z_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_xx_xy[i] = -2.0 * g_0_z_xx_xy[i] * a_exp + 4.0 * g_yy_z_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_xx_xz[i] = -2.0 * g_0_z_xx_xz[i] * a_exp + 4.0 * g_yy_z_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_xx_yy[i] = -2.0 * g_0_z_xx_yy[i] * a_exp + 4.0 * g_yy_z_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_xx_yz[i] = -2.0 * g_0_z_xx_yz[i] * a_exp + 4.0 * g_yy_z_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_xx_zz[i] = -2.0 * g_0_z_xx_zz[i] * a_exp + 4.0 * g_yy_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_0_z_xy_xx, g_0_z_xy_xy, g_0_z_xy_xz, g_0_z_xy_yy, g_0_z_xy_yz, g_0_z_xy_zz, g_yy_0_0_0_0_z_xy_xx, g_yy_0_0_0_0_z_xy_xy, g_yy_0_0_0_0_z_xy_xz, g_yy_0_0_0_0_z_xy_yy, g_yy_0_0_0_0_z_xy_yz, g_yy_0_0_0_0_z_xy_zz, g_yy_z_xy_xx, g_yy_z_xy_xy, g_yy_z_xy_xz, g_yy_z_xy_yy, g_yy_z_xy_yz, g_yy_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_z_xy_xx[i] = -2.0 * g_0_z_xy_xx[i] * a_exp + 4.0 * g_yy_z_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_xy_xy[i] = -2.0 * g_0_z_xy_xy[i] * a_exp + 4.0 * g_yy_z_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_xy_xz[i] = -2.0 * g_0_z_xy_xz[i] * a_exp + 4.0 * g_yy_z_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_xy_yy[i] = -2.0 * g_0_z_xy_yy[i] * a_exp + 4.0 * g_yy_z_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_xy_yz[i] = -2.0 * g_0_z_xy_yz[i] * a_exp + 4.0 * g_yy_z_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_xy_zz[i] = -2.0 * g_0_z_xy_zz[i] * a_exp + 4.0 * g_yy_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_0_z_xz_xx, g_0_z_xz_xy, g_0_z_xz_xz, g_0_z_xz_yy, g_0_z_xz_yz, g_0_z_xz_zz, g_yy_0_0_0_0_z_xz_xx, g_yy_0_0_0_0_z_xz_xy, g_yy_0_0_0_0_z_xz_xz, g_yy_0_0_0_0_z_xz_yy, g_yy_0_0_0_0_z_xz_yz, g_yy_0_0_0_0_z_xz_zz, g_yy_z_xz_xx, g_yy_z_xz_xy, g_yy_z_xz_xz, g_yy_z_xz_yy, g_yy_z_xz_yz, g_yy_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_z_xz_xx[i] = -2.0 * g_0_z_xz_xx[i] * a_exp + 4.0 * g_yy_z_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_xz_xy[i] = -2.0 * g_0_z_xz_xy[i] * a_exp + 4.0 * g_yy_z_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_xz_xz[i] = -2.0 * g_0_z_xz_xz[i] * a_exp + 4.0 * g_yy_z_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_xz_yy[i] = -2.0 * g_0_z_xz_yy[i] * a_exp + 4.0 * g_yy_z_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_xz_yz[i] = -2.0 * g_0_z_xz_yz[i] * a_exp + 4.0 * g_yy_z_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_xz_zz[i] = -2.0 * g_0_z_xz_zz[i] * a_exp + 4.0 * g_yy_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_0_z_yy_xx, g_0_z_yy_xy, g_0_z_yy_xz, g_0_z_yy_yy, g_0_z_yy_yz, g_0_z_yy_zz, g_yy_0_0_0_0_z_yy_xx, g_yy_0_0_0_0_z_yy_xy, g_yy_0_0_0_0_z_yy_xz, g_yy_0_0_0_0_z_yy_yy, g_yy_0_0_0_0_z_yy_yz, g_yy_0_0_0_0_z_yy_zz, g_yy_z_yy_xx, g_yy_z_yy_xy, g_yy_z_yy_xz, g_yy_z_yy_yy, g_yy_z_yy_yz, g_yy_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_z_yy_xx[i] = -2.0 * g_0_z_yy_xx[i] * a_exp + 4.0 * g_yy_z_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_yy_xy[i] = -2.0 * g_0_z_yy_xy[i] * a_exp + 4.0 * g_yy_z_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_yy_xz[i] = -2.0 * g_0_z_yy_xz[i] * a_exp + 4.0 * g_yy_z_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_yy_yy[i] = -2.0 * g_0_z_yy_yy[i] * a_exp + 4.0 * g_yy_z_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_yy_yz[i] = -2.0 * g_0_z_yy_yz[i] * a_exp + 4.0 * g_yy_z_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_yy_zz[i] = -2.0 * g_0_z_yy_zz[i] * a_exp + 4.0 * g_yy_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_0_z_yz_xx, g_0_z_yz_xy, g_0_z_yz_xz, g_0_z_yz_yy, g_0_z_yz_yz, g_0_z_yz_zz, g_yy_0_0_0_0_z_yz_xx, g_yy_0_0_0_0_z_yz_xy, g_yy_0_0_0_0_z_yz_xz, g_yy_0_0_0_0_z_yz_yy, g_yy_0_0_0_0_z_yz_yz, g_yy_0_0_0_0_z_yz_zz, g_yy_z_yz_xx, g_yy_z_yz_xy, g_yy_z_yz_xz, g_yy_z_yz_yy, g_yy_z_yz_yz, g_yy_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_z_yz_xx[i] = -2.0 * g_0_z_yz_xx[i] * a_exp + 4.0 * g_yy_z_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_yz_xy[i] = -2.0 * g_0_z_yz_xy[i] * a_exp + 4.0 * g_yy_z_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_yz_xz[i] = -2.0 * g_0_z_yz_xz[i] * a_exp + 4.0 * g_yy_z_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_yz_yy[i] = -2.0 * g_0_z_yz_yy[i] * a_exp + 4.0 * g_yy_z_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_yz_yz[i] = -2.0 * g_0_z_yz_yz[i] * a_exp + 4.0 * g_yy_z_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_yz_zz[i] = -2.0 * g_0_z_yz_zz[i] * a_exp + 4.0 * g_yy_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_0_z_zz_xx, g_0_z_zz_xy, g_0_z_zz_xz, g_0_z_zz_yy, g_0_z_zz_yz, g_0_z_zz_zz, g_yy_0_0_0_0_z_zz_xx, g_yy_0_0_0_0_z_zz_xy, g_yy_0_0_0_0_z_zz_xz, g_yy_0_0_0_0_z_zz_yy, g_yy_0_0_0_0_z_zz_yz, g_yy_0_0_0_0_z_zz_zz, g_yy_z_zz_xx, g_yy_z_zz_xy, g_yy_z_zz_xz, g_yy_z_zz_yy, g_yy_z_zz_yz, g_yy_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_z_zz_xx[i] = -2.0 * g_0_z_zz_xx[i] * a_exp + 4.0 * g_yy_z_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_zz_xy[i] = -2.0 * g_0_z_zz_xy[i] * a_exp + 4.0 * g_yy_z_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_zz_xz[i] = -2.0 * g_0_z_zz_xz[i] * a_exp + 4.0 * g_yy_z_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_zz_yy[i] = -2.0 * g_0_z_zz_yy[i] * a_exp + 4.0 * g_yy_z_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_zz_yz[i] = -2.0 * g_0_z_zz_yz[i] * a_exp + 4.0 * g_yy_z_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_zz_zz[i] = -2.0 * g_0_z_zz_zz[i] * a_exp + 4.0 * g_yy_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_yz_0_0_0_0_x_xx_xx, g_yz_0_0_0_0_x_xx_xy, g_yz_0_0_0_0_x_xx_xz, g_yz_0_0_0_0_x_xx_yy, g_yz_0_0_0_0_x_xx_yz, g_yz_0_0_0_0_x_xx_zz, g_yz_x_xx_xx, g_yz_x_xx_xy, g_yz_x_xx_xz, g_yz_x_xx_yy, g_yz_x_xx_yz, g_yz_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_x_xx_xx[i] = 4.0 * g_yz_x_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_xx_xy[i] = 4.0 * g_yz_x_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_xx_xz[i] = 4.0 * g_yz_x_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_xx_yy[i] = 4.0 * g_yz_x_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_xx_yz[i] = 4.0 * g_yz_x_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_xx_zz[i] = 4.0 * g_yz_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_yz_0_0_0_0_x_xy_xx, g_yz_0_0_0_0_x_xy_xy, g_yz_0_0_0_0_x_xy_xz, g_yz_0_0_0_0_x_xy_yy, g_yz_0_0_0_0_x_xy_yz, g_yz_0_0_0_0_x_xy_zz, g_yz_x_xy_xx, g_yz_x_xy_xy, g_yz_x_xy_xz, g_yz_x_xy_yy, g_yz_x_xy_yz, g_yz_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_x_xy_xx[i] = 4.0 * g_yz_x_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_xy_xy[i] = 4.0 * g_yz_x_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_xy_xz[i] = 4.0 * g_yz_x_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_xy_yy[i] = 4.0 * g_yz_x_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_xy_yz[i] = 4.0 * g_yz_x_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_xy_zz[i] = 4.0 * g_yz_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_yz_0_0_0_0_x_xz_xx, g_yz_0_0_0_0_x_xz_xy, g_yz_0_0_0_0_x_xz_xz, g_yz_0_0_0_0_x_xz_yy, g_yz_0_0_0_0_x_xz_yz, g_yz_0_0_0_0_x_xz_zz, g_yz_x_xz_xx, g_yz_x_xz_xy, g_yz_x_xz_xz, g_yz_x_xz_yy, g_yz_x_xz_yz, g_yz_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_x_xz_xx[i] = 4.0 * g_yz_x_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_xz_xy[i] = 4.0 * g_yz_x_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_xz_xz[i] = 4.0 * g_yz_x_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_xz_yy[i] = 4.0 * g_yz_x_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_xz_yz[i] = 4.0 * g_yz_x_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_xz_zz[i] = 4.0 * g_yz_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_yz_0_0_0_0_x_yy_xx, g_yz_0_0_0_0_x_yy_xy, g_yz_0_0_0_0_x_yy_xz, g_yz_0_0_0_0_x_yy_yy, g_yz_0_0_0_0_x_yy_yz, g_yz_0_0_0_0_x_yy_zz, g_yz_x_yy_xx, g_yz_x_yy_xy, g_yz_x_yy_xz, g_yz_x_yy_yy, g_yz_x_yy_yz, g_yz_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_x_yy_xx[i] = 4.0 * g_yz_x_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_yy_xy[i] = 4.0 * g_yz_x_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_yy_xz[i] = 4.0 * g_yz_x_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_yy_yy[i] = 4.0 * g_yz_x_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_yy_yz[i] = 4.0 * g_yz_x_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_yy_zz[i] = 4.0 * g_yz_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_yz_0_0_0_0_x_yz_xx, g_yz_0_0_0_0_x_yz_xy, g_yz_0_0_0_0_x_yz_xz, g_yz_0_0_0_0_x_yz_yy, g_yz_0_0_0_0_x_yz_yz, g_yz_0_0_0_0_x_yz_zz, g_yz_x_yz_xx, g_yz_x_yz_xy, g_yz_x_yz_xz, g_yz_x_yz_yy, g_yz_x_yz_yz, g_yz_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_x_yz_xx[i] = 4.0 * g_yz_x_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_yz_xy[i] = 4.0 * g_yz_x_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_yz_xz[i] = 4.0 * g_yz_x_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_yz_yy[i] = 4.0 * g_yz_x_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_yz_yz[i] = 4.0 * g_yz_x_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_yz_zz[i] = 4.0 * g_yz_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_yz_0_0_0_0_x_zz_xx, g_yz_0_0_0_0_x_zz_xy, g_yz_0_0_0_0_x_zz_xz, g_yz_0_0_0_0_x_zz_yy, g_yz_0_0_0_0_x_zz_yz, g_yz_0_0_0_0_x_zz_zz, g_yz_x_zz_xx, g_yz_x_zz_xy, g_yz_x_zz_xz, g_yz_x_zz_yy, g_yz_x_zz_yz, g_yz_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_x_zz_xx[i] = 4.0 * g_yz_x_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_zz_xy[i] = 4.0 * g_yz_x_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_zz_xz[i] = 4.0 * g_yz_x_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_zz_yy[i] = 4.0 * g_yz_x_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_zz_yz[i] = 4.0 * g_yz_x_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_zz_zz[i] = 4.0 * g_yz_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_yz_0_0_0_0_y_xx_xx, g_yz_0_0_0_0_y_xx_xy, g_yz_0_0_0_0_y_xx_xz, g_yz_0_0_0_0_y_xx_yy, g_yz_0_0_0_0_y_xx_yz, g_yz_0_0_0_0_y_xx_zz, g_yz_y_xx_xx, g_yz_y_xx_xy, g_yz_y_xx_xz, g_yz_y_xx_yy, g_yz_y_xx_yz, g_yz_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_y_xx_xx[i] = 4.0 * g_yz_y_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_xx_xy[i] = 4.0 * g_yz_y_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_xx_xz[i] = 4.0 * g_yz_y_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_xx_yy[i] = 4.0 * g_yz_y_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_xx_yz[i] = 4.0 * g_yz_y_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_xx_zz[i] = 4.0 * g_yz_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_yz_0_0_0_0_y_xy_xx, g_yz_0_0_0_0_y_xy_xy, g_yz_0_0_0_0_y_xy_xz, g_yz_0_0_0_0_y_xy_yy, g_yz_0_0_0_0_y_xy_yz, g_yz_0_0_0_0_y_xy_zz, g_yz_y_xy_xx, g_yz_y_xy_xy, g_yz_y_xy_xz, g_yz_y_xy_yy, g_yz_y_xy_yz, g_yz_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_y_xy_xx[i] = 4.0 * g_yz_y_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_xy_xy[i] = 4.0 * g_yz_y_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_xy_xz[i] = 4.0 * g_yz_y_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_xy_yy[i] = 4.0 * g_yz_y_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_xy_yz[i] = 4.0 * g_yz_y_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_xy_zz[i] = 4.0 * g_yz_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_yz_0_0_0_0_y_xz_xx, g_yz_0_0_0_0_y_xz_xy, g_yz_0_0_0_0_y_xz_xz, g_yz_0_0_0_0_y_xz_yy, g_yz_0_0_0_0_y_xz_yz, g_yz_0_0_0_0_y_xz_zz, g_yz_y_xz_xx, g_yz_y_xz_xy, g_yz_y_xz_xz, g_yz_y_xz_yy, g_yz_y_xz_yz, g_yz_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_y_xz_xx[i] = 4.0 * g_yz_y_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_xz_xy[i] = 4.0 * g_yz_y_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_xz_xz[i] = 4.0 * g_yz_y_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_xz_yy[i] = 4.0 * g_yz_y_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_xz_yz[i] = 4.0 * g_yz_y_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_xz_zz[i] = 4.0 * g_yz_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_yz_0_0_0_0_y_yy_xx, g_yz_0_0_0_0_y_yy_xy, g_yz_0_0_0_0_y_yy_xz, g_yz_0_0_0_0_y_yy_yy, g_yz_0_0_0_0_y_yy_yz, g_yz_0_0_0_0_y_yy_zz, g_yz_y_yy_xx, g_yz_y_yy_xy, g_yz_y_yy_xz, g_yz_y_yy_yy, g_yz_y_yy_yz, g_yz_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_y_yy_xx[i] = 4.0 * g_yz_y_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_yy_xy[i] = 4.0 * g_yz_y_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_yy_xz[i] = 4.0 * g_yz_y_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_yy_yy[i] = 4.0 * g_yz_y_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_yy_yz[i] = 4.0 * g_yz_y_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_yy_zz[i] = 4.0 * g_yz_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_yz_0_0_0_0_y_yz_xx, g_yz_0_0_0_0_y_yz_xy, g_yz_0_0_0_0_y_yz_xz, g_yz_0_0_0_0_y_yz_yy, g_yz_0_0_0_0_y_yz_yz, g_yz_0_0_0_0_y_yz_zz, g_yz_y_yz_xx, g_yz_y_yz_xy, g_yz_y_yz_xz, g_yz_y_yz_yy, g_yz_y_yz_yz, g_yz_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_y_yz_xx[i] = 4.0 * g_yz_y_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_yz_xy[i] = 4.0 * g_yz_y_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_yz_xz[i] = 4.0 * g_yz_y_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_yz_yy[i] = 4.0 * g_yz_y_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_yz_yz[i] = 4.0 * g_yz_y_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_yz_zz[i] = 4.0 * g_yz_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_yz_0_0_0_0_y_zz_xx, g_yz_0_0_0_0_y_zz_xy, g_yz_0_0_0_0_y_zz_xz, g_yz_0_0_0_0_y_zz_yy, g_yz_0_0_0_0_y_zz_yz, g_yz_0_0_0_0_y_zz_zz, g_yz_y_zz_xx, g_yz_y_zz_xy, g_yz_y_zz_xz, g_yz_y_zz_yy, g_yz_y_zz_yz, g_yz_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_y_zz_xx[i] = 4.0 * g_yz_y_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_zz_xy[i] = 4.0 * g_yz_y_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_zz_xz[i] = 4.0 * g_yz_y_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_zz_yy[i] = 4.0 * g_yz_y_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_zz_yz[i] = 4.0 * g_yz_y_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_zz_zz[i] = 4.0 * g_yz_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_yz_0_0_0_0_z_xx_xx, g_yz_0_0_0_0_z_xx_xy, g_yz_0_0_0_0_z_xx_xz, g_yz_0_0_0_0_z_xx_yy, g_yz_0_0_0_0_z_xx_yz, g_yz_0_0_0_0_z_xx_zz, g_yz_z_xx_xx, g_yz_z_xx_xy, g_yz_z_xx_xz, g_yz_z_xx_yy, g_yz_z_xx_yz, g_yz_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_z_xx_xx[i] = 4.0 * g_yz_z_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_xx_xy[i] = 4.0 * g_yz_z_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_xx_xz[i] = 4.0 * g_yz_z_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_xx_yy[i] = 4.0 * g_yz_z_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_xx_yz[i] = 4.0 * g_yz_z_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_xx_zz[i] = 4.0 * g_yz_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_yz_0_0_0_0_z_xy_xx, g_yz_0_0_0_0_z_xy_xy, g_yz_0_0_0_0_z_xy_xz, g_yz_0_0_0_0_z_xy_yy, g_yz_0_0_0_0_z_xy_yz, g_yz_0_0_0_0_z_xy_zz, g_yz_z_xy_xx, g_yz_z_xy_xy, g_yz_z_xy_xz, g_yz_z_xy_yy, g_yz_z_xy_yz, g_yz_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_z_xy_xx[i] = 4.0 * g_yz_z_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_xy_xy[i] = 4.0 * g_yz_z_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_xy_xz[i] = 4.0 * g_yz_z_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_xy_yy[i] = 4.0 * g_yz_z_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_xy_yz[i] = 4.0 * g_yz_z_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_xy_zz[i] = 4.0 * g_yz_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_yz_0_0_0_0_z_xz_xx, g_yz_0_0_0_0_z_xz_xy, g_yz_0_0_0_0_z_xz_xz, g_yz_0_0_0_0_z_xz_yy, g_yz_0_0_0_0_z_xz_yz, g_yz_0_0_0_0_z_xz_zz, g_yz_z_xz_xx, g_yz_z_xz_xy, g_yz_z_xz_xz, g_yz_z_xz_yy, g_yz_z_xz_yz, g_yz_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_z_xz_xx[i] = 4.0 * g_yz_z_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_xz_xy[i] = 4.0 * g_yz_z_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_xz_xz[i] = 4.0 * g_yz_z_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_xz_yy[i] = 4.0 * g_yz_z_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_xz_yz[i] = 4.0 * g_yz_z_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_xz_zz[i] = 4.0 * g_yz_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_yz_0_0_0_0_z_yy_xx, g_yz_0_0_0_0_z_yy_xy, g_yz_0_0_0_0_z_yy_xz, g_yz_0_0_0_0_z_yy_yy, g_yz_0_0_0_0_z_yy_yz, g_yz_0_0_0_0_z_yy_zz, g_yz_z_yy_xx, g_yz_z_yy_xy, g_yz_z_yy_xz, g_yz_z_yy_yy, g_yz_z_yy_yz, g_yz_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_z_yy_xx[i] = 4.0 * g_yz_z_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_yy_xy[i] = 4.0 * g_yz_z_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_yy_xz[i] = 4.0 * g_yz_z_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_yy_yy[i] = 4.0 * g_yz_z_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_yy_yz[i] = 4.0 * g_yz_z_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_yy_zz[i] = 4.0 * g_yz_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_yz_0_0_0_0_z_yz_xx, g_yz_0_0_0_0_z_yz_xy, g_yz_0_0_0_0_z_yz_xz, g_yz_0_0_0_0_z_yz_yy, g_yz_0_0_0_0_z_yz_yz, g_yz_0_0_0_0_z_yz_zz, g_yz_z_yz_xx, g_yz_z_yz_xy, g_yz_z_yz_xz, g_yz_z_yz_yy, g_yz_z_yz_yz, g_yz_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_z_yz_xx[i] = 4.0 * g_yz_z_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_yz_xy[i] = 4.0 * g_yz_z_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_yz_xz[i] = 4.0 * g_yz_z_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_yz_yy[i] = 4.0 * g_yz_z_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_yz_yz[i] = 4.0 * g_yz_z_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_yz_zz[i] = 4.0 * g_yz_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_yz_0_0_0_0_z_zz_xx, g_yz_0_0_0_0_z_zz_xy, g_yz_0_0_0_0_z_zz_xz, g_yz_0_0_0_0_z_zz_yy, g_yz_0_0_0_0_z_zz_yz, g_yz_0_0_0_0_z_zz_zz, g_yz_z_zz_xx, g_yz_z_zz_xy, g_yz_z_zz_xz, g_yz_z_zz_yy, g_yz_z_zz_yz, g_yz_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_z_zz_xx[i] = 4.0 * g_yz_z_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_zz_xy[i] = 4.0 * g_yz_z_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_zz_xz[i] = 4.0 * g_yz_z_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_zz_yy[i] = 4.0 * g_yz_z_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_zz_yz[i] = 4.0 * g_yz_z_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_zz_zz[i] = 4.0 * g_yz_z_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_0_x_xx_xx, g_0_x_xx_xy, g_0_x_xx_xz, g_0_x_xx_yy, g_0_x_xx_yz, g_0_x_xx_zz, g_zz_0_0_0_0_x_xx_xx, g_zz_0_0_0_0_x_xx_xy, g_zz_0_0_0_0_x_xx_xz, g_zz_0_0_0_0_x_xx_yy, g_zz_0_0_0_0_x_xx_yz, g_zz_0_0_0_0_x_xx_zz, g_zz_x_xx_xx, g_zz_x_xx_xy, g_zz_x_xx_xz, g_zz_x_xx_yy, g_zz_x_xx_yz, g_zz_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_x_xx_xx[i] = -2.0 * g_0_x_xx_xx[i] * a_exp + 4.0 * g_zz_x_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_xx_xy[i] = -2.0 * g_0_x_xx_xy[i] * a_exp + 4.0 * g_zz_x_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_xx_xz[i] = -2.0 * g_0_x_xx_xz[i] * a_exp + 4.0 * g_zz_x_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_xx_yy[i] = -2.0 * g_0_x_xx_yy[i] * a_exp + 4.0 * g_zz_x_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_xx_yz[i] = -2.0 * g_0_x_xx_yz[i] * a_exp + 4.0 * g_zz_x_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_xx_zz[i] = -2.0 * g_0_x_xx_zz[i] * a_exp + 4.0 * g_zz_x_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_0_x_xy_xx, g_0_x_xy_xy, g_0_x_xy_xz, g_0_x_xy_yy, g_0_x_xy_yz, g_0_x_xy_zz, g_zz_0_0_0_0_x_xy_xx, g_zz_0_0_0_0_x_xy_xy, g_zz_0_0_0_0_x_xy_xz, g_zz_0_0_0_0_x_xy_yy, g_zz_0_0_0_0_x_xy_yz, g_zz_0_0_0_0_x_xy_zz, g_zz_x_xy_xx, g_zz_x_xy_xy, g_zz_x_xy_xz, g_zz_x_xy_yy, g_zz_x_xy_yz, g_zz_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_x_xy_xx[i] = -2.0 * g_0_x_xy_xx[i] * a_exp + 4.0 * g_zz_x_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_xy_xy[i] = -2.0 * g_0_x_xy_xy[i] * a_exp + 4.0 * g_zz_x_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_xy_xz[i] = -2.0 * g_0_x_xy_xz[i] * a_exp + 4.0 * g_zz_x_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_xy_yy[i] = -2.0 * g_0_x_xy_yy[i] * a_exp + 4.0 * g_zz_x_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_xy_yz[i] = -2.0 * g_0_x_xy_yz[i] * a_exp + 4.0 * g_zz_x_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_xy_zz[i] = -2.0 * g_0_x_xy_zz[i] * a_exp + 4.0 * g_zz_x_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_0_x_xz_xx, g_0_x_xz_xy, g_0_x_xz_xz, g_0_x_xz_yy, g_0_x_xz_yz, g_0_x_xz_zz, g_zz_0_0_0_0_x_xz_xx, g_zz_0_0_0_0_x_xz_xy, g_zz_0_0_0_0_x_xz_xz, g_zz_0_0_0_0_x_xz_yy, g_zz_0_0_0_0_x_xz_yz, g_zz_0_0_0_0_x_xz_zz, g_zz_x_xz_xx, g_zz_x_xz_xy, g_zz_x_xz_xz, g_zz_x_xz_yy, g_zz_x_xz_yz, g_zz_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_x_xz_xx[i] = -2.0 * g_0_x_xz_xx[i] * a_exp + 4.0 * g_zz_x_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_xz_xy[i] = -2.0 * g_0_x_xz_xy[i] * a_exp + 4.0 * g_zz_x_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_xz_xz[i] = -2.0 * g_0_x_xz_xz[i] * a_exp + 4.0 * g_zz_x_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_xz_yy[i] = -2.0 * g_0_x_xz_yy[i] * a_exp + 4.0 * g_zz_x_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_xz_yz[i] = -2.0 * g_0_x_xz_yz[i] * a_exp + 4.0 * g_zz_x_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_xz_zz[i] = -2.0 * g_0_x_xz_zz[i] * a_exp + 4.0 * g_zz_x_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_0_x_yy_xx, g_0_x_yy_xy, g_0_x_yy_xz, g_0_x_yy_yy, g_0_x_yy_yz, g_0_x_yy_zz, g_zz_0_0_0_0_x_yy_xx, g_zz_0_0_0_0_x_yy_xy, g_zz_0_0_0_0_x_yy_xz, g_zz_0_0_0_0_x_yy_yy, g_zz_0_0_0_0_x_yy_yz, g_zz_0_0_0_0_x_yy_zz, g_zz_x_yy_xx, g_zz_x_yy_xy, g_zz_x_yy_xz, g_zz_x_yy_yy, g_zz_x_yy_yz, g_zz_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_x_yy_xx[i] = -2.0 * g_0_x_yy_xx[i] * a_exp + 4.0 * g_zz_x_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_yy_xy[i] = -2.0 * g_0_x_yy_xy[i] * a_exp + 4.0 * g_zz_x_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_yy_xz[i] = -2.0 * g_0_x_yy_xz[i] * a_exp + 4.0 * g_zz_x_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_yy_yy[i] = -2.0 * g_0_x_yy_yy[i] * a_exp + 4.0 * g_zz_x_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_yy_yz[i] = -2.0 * g_0_x_yy_yz[i] * a_exp + 4.0 * g_zz_x_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_yy_zz[i] = -2.0 * g_0_x_yy_zz[i] * a_exp + 4.0 * g_zz_x_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_0_x_yz_xx, g_0_x_yz_xy, g_0_x_yz_xz, g_0_x_yz_yy, g_0_x_yz_yz, g_0_x_yz_zz, g_zz_0_0_0_0_x_yz_xx, g_zz_0_0_0_0_x_yz_xy, g_zz_0_0_0_0_x_yz_xz, g_zz_0_0_0_0_x_yz_yy, g_zz_0_0_0_0_x_yz_yz, g_zz_0_0_0_0_x_yz_zz, g_zz_x_yz_xx, g_zz_x_yz_xy, g_zz_x_yz_xz, g_zz_x_yz_yy, g_zz_x_yz_yz, g_zz_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_x_yz_xx[i] = -2.0 * g_0_x_yz_xx[i] * a_exp + 4.0 * g_zz_x_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_yz_xy[i] = -2.0 * g_0_x_yz_xy[i] * a_exp + 4.0 * g_zz_x_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_yz_xz[i] = -2.0 * g_0_x_yz_xz[i] * a_exp + 4.0 * g_zz_x_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_yz_yy[i] = -2.0 * g_0_x_yz_yy[i] * a_exp + 4.0 * g_zz_x_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_yz_yz[i] = -2.0 * g_0_x_yz_yz[i] * a_exp + 4.0 * g_zz_x_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_yz_zz[i] = -2.0 * g_0_x_yz_zz[i] * a_exp + 4.0 * g_zz_x_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_0_x_zz_xx, g_0_x_zz_xy, g_0_x_zz_xz, g_0_x_zz_yy, g_0_x_zz_yz, g_0_x_zz_zz, g_zz_0_0_0_0_x_zz_xx, g_zz_0_0_0_0_x_zz_xy, g_zz_0_0_0_0_x_zz_xz, g_zz_0_0_0_0_x_zz_yy, g_zz_0_0_0_0_x_zz_yz, g_zz_0_0_0_0_x_zz_zz, g_zz_x_zz_xx, g_zz_x_zz_xy, g_zz_x_zz_xz, g_zz_x_zz_yy, g_zz_x_zz_yz, g_zz_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_x_zz_xx[i] = -2.0 * g_0_x_zz_xx[i] * a_exp + 4.0 * g_zz_x_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_zz_xy[i] = -2.0 * g_0_x_zz_xy[i] * a_exp + 4.0 * g_zz_x_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_zz_xz[i] = -2.0 * g_0_x_zz_xz[i] * a_exp + 4.0 * g_zz_x_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_zz_yy[i] = -2.0 * g_0_x_zz_yy[i] * a_exp + 4.0 * g_zz_x_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_zz_yz[i] = -2.0 * g_0_x_zz_yz[i] * a_exp + 4.0 * g_zz_x_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_zz_zz[i] = -2.0 * g_0_x_zz_zz[i] * a_exp + 4.0 * g_zz_x_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_0_y_xx_xx, g_0_y_xx_xy, g_0_y_xx_xz, g_0_y_xx_yy, g_0_y_xx_yz, g_0_y_xx_zz, g_zz_0_0_0_0_y_xx_xx, g_zz_0_0_0_0_y_xx_xy, g_zz_0_0_0_0_y_xx_xz, g_zz_0_0_0_0_y_xx_yy, g_zz_0_0_0_0_y_xx_yz, g_zz_0_0_0_0_y_xx_zz, g_zz_y_xx_xx, g_zz_y_xx_xy, g_zz_y_xx_xz, g_zz_y_xx_yy, g_zz_y_xx_yz, g_zz_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_y_xx_xx[i] = -2.0 * g_0_y_xx_xx[i] * a_exp + 4.0 * g_zz_y_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_xx_xy[i] = -2.0 * g_0_y_xx_xy[i] * a_exp + 4.0 * g_zz_y_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_xx_xz[i] = -2.0 * g_0_y_xx_xz[i] * a_exp + 4.0 * g_zz_y_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_xx_yy[i] = -2.0 * g_0_y_xx_yy[i] * a_exp + 4.0 * g_zz_y_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_xx_yz[i] = -2.0 * g_0_y_xx_yz[i] * a_exp + 4.0 * g_zz_y_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_xx_zz[i] = -2.0 * g_0_y_xx_zz[i] * a_exp + 4.0 * g_zz_y_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_0_y_xy_xx, g_0_y_xy_xy, g_0_y_xy_xz, g_0_y_xy_yy, g_0_y_xy_yz, g_0_y_xy_zz, g_zz_0_0_0_0_y_xy_xx, g_zz_0_0_0_0_y_xy_xy, g_zz_0_0_0_0_y_xy_xz, g_zz_0_0_0_0_y_xy_yy, g_zz_0_0_0_0_y_xy_yz, g_zz_0_0_0_0_y_xy_zz, g_zz_y_xy_xx, g_zz_y_xy_xy, g_zz_y_xy_xz, g_zz_y_xy_yy, g_zz_y_xy_yz, g_zz_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_y_xy_xx[i] = -2.0 * g_0_y_xy_xx[i] * a_exp + 4.0 * g_zz_y_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_xy_xy[i] = -2.0 * g_0_y_xy_xy[i] * a_exp + 4.0 * g_zz_y_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_xy_xz[i] = -2.0 * g_0_y_xy_xz[i] * a_exp + 4.0 * g_zz_y_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_xy_yy[i] = -2.0 * g_0_y_xy_yy[i] * a_exp + 4.0 * g_zz_y_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_xy_yz[i] = -2.0 * g_0_y_xy_yz[i] * a_exp + 4.0 * g_zz_y_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_xy_zz[i] = -2.0 * g_0_y_xy_zz[i] * a_exp + 4.0 * g_zz_y_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_0_y_xz_xx, g_0_y_xz_xy, g_0_y_xz_xz, g_0_y_xz_yy, g_0_y_xz_yz, g_0_y_xz_zz, g_zz_0_0_0_0_y_xz_xx, g_zz_0_0_0_0_y_xz_xy, g_zz_0_0_0_0_y_xz_xz, g_zz_0_0_0_0_y_xz_yy, g_zz_0_0_0_0_y_xz_yz, g_zz_0_0_0_0_y_xz_zz, g_zz_y_xz_xx, g_zz_y_xz_xy, g_zz_y_xz_xz, g_zz_y_xz_yy, g_zz_y_xz_yz, g_zz_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_y_xz_xx[i] = -2.0 * g_0_y_xz_xx[i] * a_exp + 4.0 * g_zz_y_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_xz_xy[i] = -2.0 * g_0_y_xz_xy[i] * a_exp + 4.0 * g_zz_y_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_xz_xz[i] = -2.0 * g_0_y_xz_xz[i] * a_exp + 4.0 * g_zz_y_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_xz_yy[i] = -2.0 * g_0_y_xz_yy[i] * a_exp + 4.0 * g_zz_y_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_xz_yz[i] = -2.0 * g_0_y_xz_yz[i] * a_exp + 4.0 * g_zz_y_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_xz_zz[i] = -2.0 * g_0_y_xz_zz[i] * a_exp + 4.0 * g_zz_y_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_0_y_yy_xx, g_0_y_yy_xy, g_0_y_yy_xz, g_0_y_yy_yy, g_0_y_yy_yz, g_0_y_yy_zz, g_zz_0_0_0_0_y_yy_xx, g_zz_0_0_0_0_y_yy_xy, g_zz_0_0_0_0_y_yy_xz, g_zz_0_0_0_0_y_yy_yy, g_zz_0_0_0_0_y_yy_yz, g_zz_0_0_0_0_y_yy_zz, g_zz_y_yy_xx, g_zz_y_yy_xy, g_zz_y_yy_xz, g_zz_y_yy_yy, g_zz_y_yy_yz, g_zz_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_y_yy_xx[i] = -2.0 * g_0_y_yy_xx[i] * a_exp + 4.0 * g_zz_y_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_yy_xy[i] = -2.0 * g_0_y_yy_xy[i] * a_exp + 4.0 * g_zz_y_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_yy_xz[i] = -2.0 * g_0_y_yy_xz[i] * a_exp + 4.0 * g_zz_y_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_yy_yy[i] = -2.0 * g_0_y_yy_yy[i] * a_exp + 4.0 * g_zz_y_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_yy_yz[i] = -2.0 * g_0_y_yy_yz[i] * a_exp + 4.0 * g_zz_y_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_yy_zz[i] = -2.0 * g_0_y_yy_zz[i] * a_exp + 4.0 * g_zz_y_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_0_y_yz_xx, g_0_y_yz_xy, g_0_y_yz_xz, g_0_y_yz_yy, g_0_y_yz_yz, g_0_y_yz_zz, g_zz_0_0_0_0_y_yz_xx, g_zz_0_0_0_0_y_yz_xy, g_zz_0_0_0_0_y_yz_xz, g_zz_0_0_0_0_y_yz_yy, g_zz_0_0_0_0_y_yz_yz, g_zz_0_0_0_0_y_yz_zz, g_zz_y_yz_xx, g_zz_y_yz_xy, g_zz_y_yz_xz, g_zz_y_yz_yy, g_zz_y_yz_yz, g_zz_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_y_yz_xx[i] = -2.0 * g_0_y_yz_xx[i] * a_exp + 4.0 * g_zz_y_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_yz_xy[i] = -2.0 * g_0_y_yz_xy[i] * a_exp + 4.0 * g_zz_y_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_yz_xz[i] = -2.0 * g_0_y_yz_xz[i] * a_exp + 4.0 * g_zz_y_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_yz_yy[i] = -2.0 * g_0_y_yz_yy[i] * a_exp + 4.0 * g_zz_y_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_yz_yz[i] = -2.0 * g_0_y_yz_yz[i] * a_exp + 4.0 * g_zz_y_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_yz_zz[i] = -2.0 * g_0_y_yz_zz[i] * a_exp + 4.0 * g_zz_y_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_0_y_zz_xx, g_0_y_zz_xy, g_0_y_zz_xz, g_0_y_zz_yy, g_0_y_zz_yz, g_0_y_zz_zz, g_zz_0_0_0_0_y_zz_xx, g_zz_0_0_0_0_y_zz_xy, g_zz_0_0_0_0_y_zz_xz, g_zz_0_0_0_0_y_zz_yy, g_zz_0_0_0_0_y_zz_yz, g_zz_0_0_0_0_y_zz_zz, g_zz_y_zz_xx, g_zz_y_zz_xy, g_zz_y_zz_xz, g_zz_y_zz_yy, g_zz_y_zz_yz, g_zz_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_y_zz_xx[i] = -2.0 * g_0_y_zz_xx[i] * a_exp + 4.0 * g_zz_y_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_zz_xy[i] = -2.0 * g_0_y_zz_xy[i] * a_exp + 4.0 * g_zz_y_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_zz_xz[i] = -2.0 * g_0_y_zz_xz[i] * a_exp + 4.0 * g_zz_y_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_zz_yy[i] = -2.0 * g_0_y_zz_yy[i] * a_exp + 4.0 * g_zz_y_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_zz_yz[i] = -2.0 * g_0_y_zz_yz[i] * a_exp + 4.0 * g_zz_y_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_zz_zz[i] = -2.0 * g_0_y_zz_zz[i] * a_exp + 4.0 * g_zz_y_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_0_z_xx_xx, g_0_z_xx_xy, g_0_z_xx_xz, g_0_z_xx_yy, g_0_z_xx_yz, g_0_z_xx_zz, g_zz_0_0_0_0_z_xx_xx, g_zz_0_0_0_0_z_xx_xy, g_zz_0_0_0_0_z_xx_xz, g_zz_0_0_0_0_z_xx_yy, g_zz_0_0_0_0_z_xx_yz, g_zz_0_0_0_0_z_xx_zz, g_zz_z_xx_xx, g_zz_z_xx_xy, g_zz_z_xx_xz, g_zz_z_xx_yy, g_zz_z_xx_yz, g_zz_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_z_xx_xx[i] = -2.0 * g_0_z_xx_xx[i] * a_exp + 4.0 * g_zz_z_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_xx_xy[i] = -2.0 * g_0_z_xx_xy[i] * a_exp + 4.0 * g_zz_z_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_xx_xz[i] = -2.0 * g_0_z_xx_xz[i] * a_exp + 4.0 * g_zz_z_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_xx_yy[i] = -2.0 * g_0_z_xx_yy[i] * a_exp + 4.0 * g_zz_z_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_xx_yz[i] = -2.0 * g_0_z_xx_yz[i] * a_exp + 4.0 * g_zz_z_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_xx_zz[i] = -2.0 * g_0_z_xx_zz[i] * a_exp + 4.0 * g_zz_z_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_0_z_xy_xx, g_0_z_xy_xy, g_0_z_xy_xz, g_0_z_xy_yy, g_0_z_xy_yz, g_0_z_xy_zz, g_zz_0_0_0_0_z_xy_xx, g_zz_0_0_0_0_z_xy_xy, g_zz_0_0_0_0_z_xy_xz, g_zz_0_0_0_0_z_xy_yy, g_zz_0_0_0_0_z_xy_yz, g_zz_0_0_0_0_z_xy_zz, g_zz_z_xy_xx, g_zz_z_xy_xy, g_zz_z_xy_xz, g_zz_z_xy_yy, g_zz_z_xy_yz, g_zz_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_z_xy_xx[i] = -2.0 * g_0_z_xy_xx[i] * a_exp + 4.0 * g_zz_z_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_xy_xy[i] = -2.0 * g_0_z_xy_xy[i] * a_exp + 4.0 * g_zz_z_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_xy_xz[i] = -2.0 * g_0_z_xy_xz[i] * a_exp + 4.0 * g_zz_z_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_xy_yy[i] = -2.0 * g_0_z_xy_yy[i] * a_exp + 4.0 * g_zz_z_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_xy_yz[i] = -2.0 * g_0_z_xy_yz[i] * a_exp + 4.0 * g_zz_z_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_xy_zz[i] = -2.0 * g_0_z_xy_zz[i] * a_exp + 4.0 * g_zz_z_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_0_z_xz_xx, g_0_z_xz_xy, g_0_z_xz_xz, g_0_z_xz_yy, g_0_z_xz_yz, g_0_z_xz_zz, g_zz_0_0_0_0_z_xz_xx, g_zz_0_0_0_0_z_xz_xy, g_zz_0_0_0_0_z_xz_xz, g_zz_0_0_0_0_z_xz_yy, g_zz_0_0_0_0_z_xz_yz, g_zz_0_0_0_0_z_xz_zz, g_zz_z_xz_xx, g_zz_z_xz_xy, g_zz_z_xz_xz, g_zz_z_xz_yy, g_zz_z_xz_yz, g_zz_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_z_xz_xx[i] = -2.0 * g_0_z_xz_xx[i] * a_exp + 4.0 * g_zz_z_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_xz_xy[i] = -2.0 * g_0_z_xz_xy[i] * a_exp + 4.0 * g_zz_z_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_xz_xz[i] = -2.0 * g_0_z_xz_xz[i] * a_exp + 4.0 * g_zz_z_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_xz_yy[i] = -2.0 * g_0_z_xz_yy[i] * a_exp + 4.0 * g_zz_z_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_xz_yz[i] = -2.0 * g_0_z_xz_yz[i] * a_exp + 4.0 * g_zz_z_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_xz_zz[i] = -2.0 * g_0_z_xz_zz[i] * a_exp + 4.0 * g_zz_z_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_0_z_yy_xx, g_0_z_yy_xy, g_0_z_yy_xz, g_0_z_yy_yy, g_0_z_yy_yz, g_0_z_yy_zz, g_zz_0_0_0_0_z_yy_xx, g_zz_0_0_0_0_z_yy_xy, g_zz_0_0_0_0_z_yy_xz, g_zz_0_0_0_0_z_yy_yy, g_zz_0_0_0_0_z_yy_yz, g_zz_0_0_0_0_z_yy_zz, g_zz_z_yy_xx, g_zz_z_yy_xy, g_zz_z_yy_xz, g_zz_z_yy_yy, g_zz_z_yy_yz, g_zz_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_z_yy_xx[i] = -2.0 * g_0_z_yy_xx[i] * a_exp + 4.0 * g_zz_z_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_yy_xy[i] = -2.0 * g_0_z_yy_xy[i] * a_exp + 4.0 * g_zz_z_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_yy_xz[i] = -2.0 * g_0_z_yy_xz[i] * a_exp + 4.0 * g_zz_z_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_yy_yy[i] = -2.0 * g_0_z_yy_yy[i] * a_exp + 4.0 * g_zz_z_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_yy_yz[i] = -2.0 * g_0_z_yy_yz[i] * a_exp + 4.0 * g_zz_z_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_yy_zz[i] = -2.0 * g_0_z_yy_zz[i] * a_exp + 4.0 * g_zz_z_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_0_z_yz_xx, g_0_z_yz_xy, g_0_z_yz_xz, g_0_z_yz_yy, g_0_z_yz_yz, g_0_z_yz_zz, g_zz_0_0_0_0_z_yz_xx, g_zz_0_0_0_0_z_yz_xy, g_zz_0_0_0_0_z_yz_xz, g_zz_0_0_0_0_z_yz_yy, g_zz_0_0_0_0_z_yz_yz, g_zz_0_0_0_0_z_yz_zz, g_zz_z_yz_xx, g_zz_z_yz_xy, g_zz_z_yz_xz, g_zz_z_yz_yy, g_zz_z_yz_yz, g_zz_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_z_yz_xx[i] = -2.0 * g_0_z_yz_xx[i] * a_exp + 4.0 * g_zz_z_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_yz_xy[i] = -2.0 * g_0_z_yz_xy[i] * a_exp + 4.0 * g_zz_z_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_yz_xz[i] = -2.0 * g_0_z_yz_xz[i] * a_exp + 4.0 * g_zz_z_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_yz_yy[i] = -2.0 * g_0_z_yz_yy[i] * a_exp + 4.0 * g_zz_z_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_yz_yz[i] = -2.0 * g_0_z_yz_yz[i] * a_exp + 4.0 * g_zz_z_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_yz_zz[i] = -2.0 * g_0_z_yz_zz[i] * a_exp + 4.0 * g_zz_z_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_0_z_zz_xx, g_0_z_zz_xy, g_0_z_zz_xz, g_0_z_zz_yy, g_0_z_zz_yz, g_0_z_zz_zz, g_zz_0_0_0_0_z_zz_xx, g_zz_0_0_0_0_z_zz_xy, g_zz_0_0_0_0_z_zz_xz, g_zz_0_0_0_0_z_zz_yy, g_zz_0_0_0_0_z_zz_yz, g_zz_0_0_0_0_z_zz_zz, g_zz_z_zz_xx, g_zz_z_zz_xy, g_zz_z_zz_xz, g_zz_z_zz_yy, g_zz_z_zz_yz, g_zz_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_z_zz_xx[i] = -2.0 * g_0_z_zz_xx[i] * a_exp + 4.0 * g_zz_z_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_zz_xy[i] = -2.0 * g_0_z_zz_xy[i] * a_exp + 4.0 * g_zz_z_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_zz_xz[i] = -2.0 * g_0_z_zz_xz[i] * a_exp + 4.0 * g_zz_z_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_zz_yy[i] = -2.0 * g_0_z_zz_yy[i] * a_exp + 4.0 * g_zz_z_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_zz_yz[i] = -2.0 * g_0_z_zz_yz[i] * a_exp + 4.0 * g_zz_z_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_zz_zz[i] = -2.0 * g_0_z_zz_zz[i] * a_exp + 4.0 * g_zz_z_zz_zz[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

