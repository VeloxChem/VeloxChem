#include "GeomDeriv1000OfScalarForPPDD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_ppdd_0(CSimdArray<double>& buffer_1000_ppdd,
                     const CSimdArray<double>& buffer_spdd,
                     const CSimdArray<double>& buffer_dpdd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_ppdd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_ppdd

    auto g_x_0_0_0_x_x_xx_xx = buffer_1000_ppdd[0];

    auto g_x_0_0_0_x_x_xx_xy = buffer_1000_ppdd[1];

    auto g_x_0_0_0_x_x_xx_xz = buffer_1000_ppdd[2];

    auto g_x_0_0_0_x_x_xx_yy = buffer_1000_ppdd[3];

    auto g_x_0_0_0_x_x_xx_yz = buffer_1000_ppdd[4];

    auto g_x_0_0_0_x_x_xx_zz = buffer_1000_ppdd[5];

    auto g_x_0_0_0_x_x_xy_xx = buffer_1000_ppdd[6];

    auto g_x_0_0_0_x_x_xy_xy = buffer_1000_ppdd[7];

    auto g_x_0_0_0_x_x_xy_xz = buffer_1000_ppdd[8];

    auto g_x_0_0_0_x_x_xy_yy = buffer_1000_ppdd[9];

    auto g_x_0_0_0_x_x_xy_yz = buffer_1000_ppdd[10];

    auto g_x_0_0_0_x_x_xy_zz = buffer_1000_ppdd[11];

    auto g_x_0_0_0_x_x_xz_xx = buffer_1000_ppdd[12];

    auto g_x_0_0_0_x_x_xz_xy = buffer_1000_ppdd[13];

    auto g_x_0_0_0_x_x_xz_xz = buffer_1000_ppdd[14];

    auto g_x_0_0_0_x_x_xz_yy = buffer_1000_ppdd[15];

    auto g_x_0_0_0_x_x_xz_yz = buffer_1000_ppdd[16];

    auto g_x_0_0_0_x_x_xz_zz = buffer_1000_ppdd[17];

    auto g_x_0_0_0_x_x_yy_xx = buffer_1000_ppdd[18];

    auto g_x_0_0_0_x_x_yy_xy = buffer_1000_ppdd[19];

    auto g_x_0_0_0_x_x_yy_xz = buffer_1000_ppdd[20];

    auto g_x_0_0_0_x_x_yy_yy = buffer_1000_ppdd[21];

    auto g_x_0_0_0_x_x_yy_yz = buffer_1000_ppdd[22];

    auto g_x_0_0_0_x_x_yy_zz = buffer_1000_ppdd[23];

    auto g_x_0_0_0_x_x_yz_xx = buffer_1000_ppdd[24];

    auto g_x_0_0_0_x_x_yz_xy = buffer_1000_ppdd[25];

    auto g_x_0_0_0_x_x_yz_xz = buffer_1000_ppdd[26];

    auto g_x_0_0_0_x_x_yz_yy = buffer_1000_ppdd[27];

    auto g_x_0_0_0_x_x_yz_yz = buffer_1000_ppdd[28];

    auto g_x_0_0_0_x_x_yz_zz = buffer_1000_ppdd[29];

    auto g_x_0_0_0_x_x_zz_xx = buffer_1000_ppdd[30];

    auto g_x_0_0_0_x_x_zz_xy = buffer_1000_ppdd[31];

    auto g_x_0_0_0_x_x_zz_xz = buffer_1000_ppdd[32];

    auto g_x_0_0_0_x_x_zz_yy = buffer_1000_ppdd[33];

    auto g_x_0_0_0_x_x_zz_yz = buffer_1000_ppdd[34];

    auto g_x_0_0_0_x_x_zz_zz = buffer_1000_ppdd[35];

    auto g_x_0_0_0_x_y_xx_xx = buffer_1000_ppdd[36];

    auto g_x_0_0_0_x_y_xx_xy = buffer_1000_ppdd[37];

    auto g_x_0_0_0_x_y_xx_xz = buffer_1000_ppdd[38];

    auto g_x_0_0_0_x_y_xx_yy = buffer_1000_ppdd[39];

    auto g_x_0_0_0_x_y_xx_yz = buffer_1000_ppdd[40];

    auto g_x_0_0_0_x_y_xx_zz = buffer_1000_ppdd[41];

    auto g_x_0_0_0_x_y_xy_xx = buffer_1000_ppdd[42];

    auto g_x_0_0_0_x_y_xy_xy = buffer_1000_ppdd[43];

    auto g_x_0_0_0_x_y_xy_xz = buffer_1000_ppdd[44];

    auto g_x_0_0_0_x_y_xy_yy = buffer_1000_ppdd[45];

    auto g_x_0_0_0_x_y_xy_yz = buffer_1000_ppdd[46];

    auto g_x_0_0_0_x_y_xy_zz = buffer_1000_ppdd[47];

    auto g_x_0_0_0_x_y_xz_xx = buffer_1000_ppdd[48];

    auto g_x_0_0_0_x_y_xz_xy = buffer_1000_ppdd[49];

    auto g_x_0_0_0_x_y_xz_xz = buffer_1000_ppdd[50];

    auto g_x_0_0_0_x_y_xz_yy = buffer_1000_ppdd[51];

    auto g_x_0_0_0_x_y_xz_yz = buffer_1000_ppdd[52];

    auto g_x_0_0_0_x_y_xz_zz = buffer_1000_ppdd[53];

    auto g_x_0_0_0_x_y_yy_xx = buffer_1000_ppdd[54];

    auto g_x_0_0_0_x_y_yy_xy = buffer_1000_ppdd[55];

    auto g_x_0_0_0_x_y_yy_xz = buffer_1000_ppdd[56];

    auto g_x_0_0_0_x_y_yy_yy = buffer_1000_ppdd[57];

    auto g_x_0_0_0_x_y_yy_yz = buffer_1000_ppdd[58];

    auto g_x_0_0_0_x_y_yy_zz = buffer_1000_ppdd[59];

    auto g_x_0_0_0_x_y_yz_xx = buffer_1000_ppdd[60];

    auto g_x_0_0_0_x_y_yz_xy = buffer_1000_ppdd[61];

    auto g_x_0_0_0_x_y_yz_xz = buffer_1000_ppdd[62];

    auto g_x_0_0_0_x_y_yz_yy = buffer_1000_ppdd[63];

    auto g_x_0_0_0_x_y_yz_yz = buffer_1000_ppdd[64];

    auto g_x_0_0_0_x_y_yz_zz = buffer_1000_ppdd[65];

    auto g_x_0_0_0_x_y_zz_xx = buffer_1000_ppdd[66];

    auto g_x_0_0_0_x_y_zz_xy = buffer_1000_ppdd[67];

    auto g_x_0_0_0_x_y_zz_xz = buffer_1000_ppdd[68];

    auto g_x_0_0_0_x_y_zz_yy = buffer_1000_ppdd[69];

    auto g_x_0_0_0_x_y_zz_yz = buffer_1000_ppdd[70];

    auto g_x_0_0_0_x_y_zz_zz = buffer_1000_ppdd[71];

    auto g_x_0_0_0_x_z_xx_xx = buffer_1000_ppdd[72];

    auto g_x_0_0_0_x_z_xx_xy = buffer_1000_ppdd[73];

    auto g_x_0_0_0_x_z_xx_xz = buffer_1000_ppdd[74];

    auto g_x_0_0_0_x_z_xx_yy = buffer_1000_ppdd[75];

    auto g_x_0_0_0_x_z_xx_yz = buffer_1000_ppdd[76];

    auto g_x_0_0_0_x_z_xx_zz = buffer_1000_ppdd[77];

    auto g_x_0_0_0_x_z_xy_xx = buffer_1000_ppdd[78];

    auto g_x_0_0_0_x_z_xy_xy = buffer_1000_ppdd[79];

    auto g_x_0_0_0_x_z_xy_xz = buffer_1000_ppdd[80];

    auto g_x_0_0_0_x_z_xy_yy = buffer_1000_ppdd[81];

    auto g_x_0_0_0_x_z_xy_yz = buffer_1000_ppdd[82];

    auto g_x_0_0_0_x_z_xy_zz = buffer_1000_ppdd[83];

    auto g_x_0_0_0_x_z_xz_xx = buffer_1000_ppdd[84];

    auto g_x_0_0_0_x_z_xz_xy = buffer_1000_ppdd[85];

    auto g_x_0_0_0_x_z_xz_xz = buffer_1000_ppdd[86];

    auto g_x_0_0_0_x_z_xz_yy = buffer_1000_ppdd[87];

    auto g_x_0_0_0_x_z_xz_yz = buffer_1000_ppdd[88];

    auto g_x_0_0_0_x_z_xz_zz = buffer_1000_ppdd[89];

    auto g_x_0_0_0_x_z_yy_xx = buffer_1000_ppdd[90];

    auto g_x_0_0_0_x_z_yy_xy = buffer_1000_ppdd[91];

    auto g_x_0_0_0_x_z_yy_xz = buffer_1000_ppdd[92];

    auto g_x_0_0_0_x_z_yy_yy = buffer_1000_ppdd[93];

    auto g_x_0_0_0_x_z_yy_yz = buffer_1000_ppdd[94];

    auto g_x_0_0_0_x_z_yy_zz = buffer_1000_ppdd[95];

    auto g_x_0_0_0_x_z_yz_xx = buffer_1000_ppdd[96];

    auto g_x_0_0_0_x_z_yz_xy = buffer_1000_ppdd[97];

    auto g_x_0_0_0_x_z_yz_xz = buffer_1000_ppdd[98];

    auto g_x_0_0_0_x_z_yz_yy = buffer_1000_ppdd[99];

    auto g_x_0_0_0_x_z_yz_yz = buffer_1000_ppdd[100];

    auto g_x_0_0_0_x_z_yz_zz = buffer_1000_ppdd[101];

    auto g_x_0_0_0_x_z_zz_xx = buffer_1000_ppdd[102];

    auto g_x_0_0_0_x_z_zz_xy = buffer_1000_ppdd[103];

    auto g_x_0_0_0_x_z_zz_xz = buffer_1000_ppdd[104];

    auto g_x_0_0_0_x_z_zz_yy = buffer_1000_ppdd[105];

    auto g_x_0_0_0_x_z_zz_yz = buffer_1000_ppdd[106];

    auto g_x_0_0_0_x_z_zz_zz = buffer_1000_ppdd[107];

    auto g_x_0_0_0_y_x_xx_xx = buffer_1000_ppdd[108];

    auto g_x_0_0_0_y_x_xx_xy = buffer_1000_ppdd[109];

    auto g_x_0_0_0_y_x_xx_xz = buffer_1000_ppdd[110];

    auto g_x_0_0_0_y_x_xx_yy = buffer_1000_ppdd[111];

    auto g_x_0_0_0_y_x_xx_yz = buffer_1000_ppdd[112];

    auto g_x_0_0_0_y_x_xx_zz = buffer_1000_ppdd[113];

    auto g_x_0_0_0_y_x_xy_xx = buffer_1000_ppdd[114];

    auto g_x_0_0_0_y_x_xy_xy = buffer_1000_ppdd[115];

    auto g_x_0_0_0_y_x_xy_xz = buffer_1000_ppdd[116];

    auto g_x_0_0_0_y_x_xy_yy = buffer_1000_ppdd[117];

    auto g_x_0_0_0_y_x_xy_yz = buffer_1000_ppdd[118];

    auto g_x_0_0_0_y_x_xy_zz = buffer_1000_ppdd[119];

    auto g_x_0_0_0_y_x_xz_xx = buffer_1000_ppdd[120];

    auto g_x_0_0_0_y_x_xz_xy = buffer_1000_ppdd[121];

    auto g_x_0_0_0_y_x_xz_xz = buffer_1000_ppdd[122];

    auto g_x_0_0_0_y_x_xz_yy = buffer_1000_ppdd[123];

    auto g_x_0_0_0_y_x_xz_yz = buffer_1000_ppdd[124];

    auto g_x_0_0_0_y_x_xz_zz = buffer_1000_ppdd[125];

    auto g_x_0_0_0_y_x_yy_xx = buffer_1000_ppdd[126];

    auto g_x_0_0_0_y_x_yy_xy = buffer_1000_ppdd[127];

    auto g_x_0_0_0_y_x_yy_xz = buffer_1000_ppdd[128];

    auto g_x_0_0_0_y_x_yy_yy = buffer_1000_ppdd[129];

    auto g_x_0_0_0_y_x_yy_yz = buffer_1000_ppdd[130];

    auto g_x_0_0_0_y_x_yy_zz = buffer_1000_ppdd[131];

    auto g_x_0_0_0_y_x_yz_xx = buffer_1000_ppdd[132];

    auto g_x_0_0_0_y_x_yz_xy = buffer_1000_ppdd[133];

    auto g_x_0_0_0_y_x_yz_xz = buffer_1000_ppdd[134];

    auto g_x_0_0_0_y_x_yz_yy = buffer_1000_ppdd[135];

    auto g_x_0_0_0_y_x_yz_yz = buffer_1000_ppdd[136];

    auto g_x_0_0_0_y_x_yz_zz = buffer_1000_ppdd[137];

    auto g_x_0_0_0_y_x_zz_xx = buffer_1000_ppdd[138];

    auto g_x_0_0_0_y_x_zz_xy = buffer_1000_ppdd[139];

    auto g_x_0_0_0_y_x_zz_xz = buffer_1000_ppdd[140];

    auto g_x_0_0_0_y_x_zz_yy = buffer_1000_ppdd[141];

    auto g_x_0_0_0_y_x_zz_yz = buffer_1000_ppdd[142];

    auto g_x_0_0_0_y_x_zz_zz = buffer_1000_ppdd[143];

    auto g_x_0_0_0_y_y_xx_xx = buffer_1000_ppdd[144];

    auto g_x_0_0_0_y_y_xx_xy = buffer_1000_ppdd[145];

    auto g_x_0_0_0_y_y_xx_xz = buffer_1000_ppdd[146];

    auto g_x_0_0_0_y_y_xx_yy = buffer_1000_ppdd[147];

    auto g_x_0_0_0_y_y_xx_yz = buffer_1000_ppdd[148];

    auto g_x_0_0_0_y_y_xx_zz = buffer_1000_ppdd[149];

    auto g_x_0_0_0_y_y_xy_xx = buffer_1000_ppdd[150];

    auto g_x_0_0_0_y_y_xy_xy = buffer_1000_ppdd[151];

    auto g_x_0_0_0_y_y_xy_xz = buffer_1000_ppdd[152];

    auto g_x_0_0_0_y_y_xy_yy = buffer_1000_ppdd[153];

    auto g_x_0_0_0_y_y_xy_yz = buffer_1000_ppdd[154];

    auto g_x_0_0_0_y_y_xy_zz = buffer_1000_ppdd[155];

    auto g_x_0_0_0_y_y_xz_xx = buffer_1000_ppdd[156];

    auto g_x_0_0_0_y_y_xz_xy = buffer_1000_ppdd[157];

    auto g_x_0_0_0_y_y_xz_xz = buffer_1000_ppdd[158];

    auto g_x_0_0_0_y_y_xz_yy = buffer_1000_ppdd[159];

    auto g_x_0_0_0_y_y_xz_yz = buffer_1000_ppdd[160];

    auto g_x_0_0_0_y_y_xz_zz = buffer_1000_ppdd[161];

    auto g_x_0_0_0_y_y_yy_xx = buffer_1000_ppdd[162];

    auto g_x_0_0_0_y_y_yy_xy = buffer_1000_ppdd[163];

    auto g_x_0_0_0_y_y_yy_xz = buffer_1000_ppdd[164];

    auto g_x_0_0_0_y_y_yy_yy = buffer_1000_ppdd[165];

    auto g_x_0_0_0_y_y_yy_yz = buffer_1000_ppdd[166];

    auto g_x_0_0_0_y_y_yy_zz = buffer_1000_ppdd[167];

    auto g_x_0_0_0_y_y_yz_xx = buffer_1000_ppdd[168];

    auto g_x_0_0_0_y_y_yz_xy = buffer_1000_ppdd[169];

    auto g_x_0_0_0_y_y_yz_xz = buffer_1000_ppdd[170];

    auto g_x_0_0_0_y_y_yz_yy = buffer_1000_ppdd[171];

    auto g_x_0_0_0_y_y_yz_yz = buffer_1000_ppdd[172];

    auto g_x_0_0_0_y_y_yz_zz = buffer_1000_ppdd[173];

    auto g_x_0_0_0_y_y_zz_xx = buffer_1000_ppdd[174];

    auto g_x_0_0_0_y_y_zz_xy = buffer_1000_ppdd[175];

    auto g_x_0_0_0_y_y_zz_xz = buffer_1000_ppdd[176];

    auto g_x_0_0_0_y_y_zz_yy = buffer_1000_ppdd[177];

    auto g_x_0_0_0_y_y_zz_yz = buffer_1000_ppdd[178];

    auto g_x_0_0_0_y_y_zz_zz = buffer_1000_ppdd[179];

    auto g_x_0_0_0_y_z_xx_xx = buffer_1000_ppdd[180];

    auto g_x_0_0_0_y_z_xx_xy = buffer_1000_ppdd[181];

    auto g_x_0_0_0_y_z_xx_xz = buffer_1000_ppdd[182];

    auto g_x_0_0_0_y_z_xx_yy = buffer_1000_ppdd[183];

    auto g_x_0_0_0_y_z_xx_yz = buffer_1000_ppdd[184];

    auto g_x_0_0_0_y_z_xx_zz = buffer_1000_ppdd[185];

    auto g_x_0_0_0_y_z_xy_xx = buffer_1000_ppdd[186];

    auto g_x_0_0_0_y_z_xy_xy = buffer_1000_ppdd[187];

    auto g_x_0_0_0_y_z_xy_xz = buffer_1000_ppdd[188];

    auto g_x_0_0_0_y_z_xy_yy = buffer_1000_ppdd[189];

    auto g_x_0_0_0_y_z_xy_yz = buffer_1000_ppdd[190];

    auto g_x_0_0_0_y_z_xy_zz = buffer_1000_ppdd[191];

    auto g_x_0_0_0_y_z_xz_xx = buffer_1000_ppdd[192];

    auto g_x_0_0_0_y_z_xz_xy = buffer_1000_ppdd[193];

    auto g_x_0_0_0_y_z_xz_xz = buffer_1000_ppdd[194];

    auto g_x_0_0_0_y_z_xz_yy = buffer_1000_ppdd[195];

    auto g_x_0_0_0_y_z_xz_yz = buffer_1000_ppdd[196];

    auto g_x_0_0_0_y_z_xz_zz = buffer_1000_ppdd[197];

    auto g_x_0_0_0_y_z_yy_xx = buffer_1000_ppdd[198];

    auto g_x_0_0_0_y_z_yy_xy = buffer_1000_ppdd[199];

    auto g_x_0_0_0_y_z_yy_xz = buffer_1000_ppdd[200];

    auto g_x_0_0_0_y_z_yy_yy = buffer_1000_ppdd[201];

    auto g_x_0_0_0_y_z_yy_yz = buffer_1000_ppdd[202];

    auto g_x_0_0_0_y_z_yy_zz = buffer_1000_ppdd[203];

    auto g_x_0_0_0_y_z_yz_xx = buffer_1000_ppdd[204];

    auto g_x_0_0_0_y_z_yz_xy = buffer_1000_ppdd[205];

    auto g_x_0_0_0_y_z_yz_xz = buffer_1000_ppdd[206];

    auto g_x_0_0_0_y_z_yz_yy = buffer_1000_ppdd[207];

    auto g_x_0_0_0_y_z_yz_yz = buffer_1000_ppdd[208];

    auto g_x_0_0_0_y_z_yz_zz = buffer_1000_ppdd[209];

    auto g_x_0_0_0_y_z_zz_xx = buffer_1000_ppdd[210];

    auto g_x_0_0_0_y_z_zz_xy = buffer_1000_ppdd[211];

    auto g_x_0_0_0_y_z_zz_xz = buffer_1000_ppdd[212];

    auto g_x_0_0_0_y_z_zz_yy = buffer_1000_ppdd[213];

    auto g_x_0_0_0_y_z_zz_yz = buffer_1000_ppdd[214];

    auto g_x_0_0_0_y_z_zz_zz = buffer_1000_ppdd[215];

    auto g_x_0_0_0_z_x_xx_xx = buffer_1000_ppdd[216];

    auto g_x_0_0_0_z_x_xx_xy = buffer_1000_ppdd[217];

    auto g_x_0_0_0_z_x_xx_xz = buffer_1000_ppdd[218];

    auto g_x_0_0_0_z_x_xx_yy = buffer_1000_ppdd[219];

    auto g_x_0_0_0_z_x_xx_yz = buffer_1000_ppdd[220];

    auto g_x_0_0_0_z_x_xx_zz = buffer_1000_ppdd[221];

    auto g_x_0_0_0_z_x_xy_xx = buffer_1000_ppdd[222];

    auto g_x_0_0_0_z_x_xy_xy = buffer_1000_ppdd[223];

    auto g_x_0_0_0_z_x_xy_xz = buffer_1000_ppdd[224];

    auto g_x_0_0_0_z_x_xy_yy = buffer_1000_ppdd[225];

    auto g_x_0_0_0_z_x_xy_yz = buffer_1000_ppdd[226];

    auto g_x_0_0_0_z_x_xy_zz = buffer_1000_ppdd[227];

    auto g_x_0_0_0_z_x_xz_xx = buffer_1000_ppdd[228];

    auto g_x_0_0_0_z_x_xz_xy = buffer_1000_ppdd[229];

    auto g_x_0_0_0_z_x_xz_xz = buffer_1000_ppdd[230];

    auto g_x_0_0_0_z_x_xz_yy = buffer_1000_ppdd[231];

    auto g_x_0_0_0_z_x_xz_yz = buffer_1000_ppdd[232];

    auto g_x_0_0_0_z_x_xz_zz = buffer_1000_ppdd[233];

    auto g_x_0_0_0_z_x_yy_xx = buffer_1000_ppdd[234];

    auto g_x_0_0_0_z_x_yy_xy = buffer_1000_ppdd[235];

    auto g_x_0_0_0_z_x_yy_xz = buffer_1000_ppdd[236];

    auto g_x_0_0_0_z_x_yy_yy = buffer_1000_ppdd[237];

    auto g_x_0_0_0_z_x_yy_yz = buffer_1000_ppdd[238];

    auto g_x_0_0_0_z_x_yy_zz = buffer_1000_ppdd[239];

    auto g_x_0_0_0_z_x_yz_xx = buffer_1000_ppdd[240];

    auto g_x_0_0_0_z_x_yz_xy = buffer_1000_ppdd[241];

    auto g_x_0_0_0_z_x_yz_xz = buffer_1000_ppdd[242];

    auto g_x_0_0_0_z_x_yz_yy = buffer_1000_ppdd[243];

    auto g_x_0_0_0_z_x_yz_yz = buffer_1000_ppdd[244];

    auto g_x_0_0_0_z_x_yz_zz = buffer_1000_ppdd[245];

    auto g_x_0_0_0_z_x_zz_xx = buffer_1000_ppdd[246];

    auto g_x_0_0_0_z_x_zz_xy = buffer_1000_ppdd[247];

    auto g_x_0_0_0_z_x_zz_xz = buffer_1000_ppdd[248];

    auto g_x_0_0_0_z_x_zz_yy = buffer_1000_ppdd[249];

    auto g_x_0_0_0_z_x_zz_yz = buffer_1000_ppdd[250];

    auto g_x_0_0_0_z_x_zz_zz = buffer_1000_ppdd[251];

    auto g_x_0_0_0_z_y_xx_xx = buffer_1000_ppdd[252];

    auto g_x_0_0_0_z_y_xx_xy = buffer_1000_ppdd[253];

    auto g_x_0_0_0_z_y_xx_xz = buffer_1000_ppdd[254];

    auto g_x_0_0_0_z_y_xx_yy = buffer_1000_ppdd[255];

    auto g_x_0_0_0_z_y_xx_yz = buffer_1000_ppdd[256];

    auto g_x_0_0_0_z_y_xx_zz = buffer_1000_ppdd[257];

    auto g_x_0_0_0_z_y_xy_xx = buffer_1000_ppdd[258];

    auto g_x_0_0_0_z_y_xy_xy = buffer_1000_ppdd[259];

    auto g_x_0_0_0_z_y_xy_xz = buffer_1000_ppdd[260];

    auto g_x_0_0_0_z_y_xy_yy = buffer_1000_ppdd[261];

    auto g_x_0_0_0_z_y_xy_yz = buffer_1000_ppdd[262];

    auto g_x_0_0_0_z_y_xy_zz = buffer_1000_ppdd[263];

    auto g_x_0_0_0_z_y_xz_xx = buffer_1000_ppdd[264];

    auto g_x_0_0_0_z_y_xz_xy = buffer_1000_ppdd[265];

    auto g_x_0_0_0_z_y_xz_xz = buffer_1000_ppdd[266];

    auto g_x_0_0_0_z_y_xz_yy = buffer_1000_ppdd[267];

    auto g_x_0_0_0_z_y_xz_yz = buffer_1000_ppdd[268];

    auto g_x_0_0_0_z_y_xz_zz = buffer_1000_ppdd[269];

    auto g_x_0_0_0_z_y_yy_xx = buffer_1000_ppdd[270];

    auto g_x_0_0_0_z_y_yy_xy = buffer_1000_ppdd[271];

    auto g_x_0_0_0_z_y_yy_xz = buffer_1000_ppdd[272];

    auto g_x_0_0_0_z_y_yy_yy = buffer_1000_ppdd[273];

    auto g_x_0_0_0_z_y_yy_yz = buffer_1000_ppdd[274];

    auto g_x_0_0_0_z_y_yy_zz = buffer_1000_ppdd[275];

    auto g_x_0_0_0_z_y_yz_xx = buffer_1000_ppdd[276];

    auto g_x_0_0_0_z_y_yz_xy = buffer_1000_ppdd[277];

    auto g_x_0_0_0_z_y_yz_xz = buffer_1000_ppdd[278];

    auto g_x_0_0_0_z_y_yz_yy = buffer_1000_ppdd[279];

    auto g_x_0_0_0_z_y_yz_yz = buffer_1000_ppdd[280];

    auto g_x_0_0_0_z_y_yz_zz = buffer_1000_ppdd[281];

    auto g_x_0_0_0_z_y_zz_xx = buffer_1000_ppdd[282];

    auto g_x_0_0_0_z_y_zz_xy = buffer_1000_ppdd[283];

    auto g_x_0_0_0_z_y_zz_xz = buffer_1000_ppdd[284];

    auto g_x_0_0_0_z_y_zz_yy = buffer_1000_ppdd[285];

    auto g_x_0_0_0_z_y_zz_yz = buffer_1000_ppdd[286];

    auto g_x_0_0_0_z_y_zz_zz = buffer_1000_ppdd[287];

    auto g_x_0_0_0_z_z_xx_xx = buffer_1000_ppdd[288];

    auto g_x_0_0_0_z_z_xx_xy = buffer_1000_ppdd[289];

    auto g_x_0_0_0_z_z_xx_xz = buffer_1000_ppdd[290];

    auto g_x_0_0_0_z_z_xx_yy = buffer_1000_ppdd[291];

    auto g_x_0_0_0_z_z_xx_yz = buffer_1000_ppdd[292];

    auto g_x_0_0_0_z_z_xx_zz = buffer_1000_ppdd[293];

    auto g_x_0_0_0_z_z_xy_xx = buffer_1000_ppdd[294];

    auto g_x_0_0_0_z_z_xy_xy = buffer_1000_ppdd[295];

    auto g_x_0_0_0_z_z_xy_xz = buffer_1000_ppdd[296];

    auto g_x_0_0_0_z_z_xy_yy = buffer_1000_ppdd[297];

    auto g_x_0_0_0_z_z_xy_yz = buffer_1000_ppdd[298];

    auto g_x_0_0_0_z_z_xy_zz = buffer_1000_ppdd[299];

    auto g_x_0_0_0_z_z_xz_xx = buffer_1000_ppdd[300];

    auto g_x_0_0_0_z_z_xz_xy = buffer_1000_ppdd[301];

    auto g_x_0_0_0_z_z_xz_xz = buffer_1000_ppdd[302];

    auto g_x_0_0_0_z_z_xz_yy = buffer_1000_ppdd[303];

    auto g_x_0_0_0_z_z_xz_yz = buffer_1000_ppdd[304];

    auto g_x_0_0_0_z_z_xz_zz = buffer_1000_ppdd[305];

    auto g_x_0_0_0_z_z_yy_xx = buffer_1000_ppdd[306];

    auto g_x_0_0_0_z_z_yy_xy = buffer_1000_ppdd[307];

    auto g_x_0_0_0_z_z_yy_xz = buffer_1000_ppdd[308];

    auto g_x_0_0_0_z_z_yy_yy = buffer_1000_ppdd[309];

    auto g_x_0_0_0_z_z_yy_yz = buffer_1000_ppdd[310];

    auto g_x_0_0_0_z_z_yy_zz = buffer_1000_ppdd[311];

    auto g_x_0_0_0_z_z_yz_xx = buffer_1000_ppdd[312];

    auto g_x_0_0_0_z_z_yz_xy = buffer_1000_ppdd[313];

    auto g_x_0_0_0_z_z_yz_xz = buffer_1000_ppdd[314];

    auto g_x_0_0_0_z_z_yz_yy = buffer_1000_ppdd[315];

    auto g_x_0_0_0_z_z_yz_yz = buffer_1000_ppdd[316];

    auto g_x_0_0_0_z_z_yz_zz = buffer_1000_ppdd[317];

    auto g_x_0_0_0_z_z_zz_xx = buffer_1000_ppdd[318];

    auto g_x_0_0_0_z_z_zz_xy = buffer_1000_ppdd[319];

    auto g_x_0_0_0_z_z_zz_xz = buffer_1000_ppdd[320];

    auto g_x_0_0_0_z_z_zz_yy = buffer_1000_ppdd[321];

    auto g_x_0_0_0_z_z_zz_yz = buffer_1000_ppdd[322];

    auto g_x_0_0_0_z_z_zz_zz = buffer_1000_ppdd[323];

    auto g_y_0_0_0_x_x_xx_xx = buffer_1000_ppdd[324];

    auto g_y_0_0_0_x_x_xx_xy = buffer_1000_ppdd[325];

    auto g_y_0_0_0_x_x_xx_xz = buffer_1000_ppdd[326];

    auto g_y_0_0_0_x_x_xx_yy = buffer_1000_ppdd[327];

    auto g_y_0_0_0_x_x_xx_yz = buffer_1000_ppdd[328];

    auto g_y_0_0_0_x_x_xx_zz = buffer_1000_ppdd[329];

    auto g_y_0_0_0_x_x_xy_xx = buffer_1000_ppdd[330];

    auto g_y_0_0_0_x_x_xy_xy = buffer_1000_ppdd[331];

    auto g_y_0_0_0_x_x_xy_xz = buffer_1000_ppdd[332];

    auto g_y_0_0_0_x_x_xy_yy = buffer_1000_ppdd[333];

    auto g_y_0_0_0_x_x_xy_yz = buffer_1000_ppdd[334];

    auto g_y_0_0_0_x_x_xy_zz = buffer_1000_ppdd[335];

    auto g_y_0_0_0_x_x_xz_xx = buffer_1000_ppdd[336];

    auto g_y_0_0_0_x_x_xz_xy = buffer_1000_ppdd[337];

    auto g_y_0_0_0_x_x_xz_xz = buffer_1000_ppdd[338];

    auto g_y_0_0_0_x_x_xz_yy = buffer_1000_ppdd[339];

    auto g_y_0_0_0_x_x_xz_yz = buffer_1000_ppdd[340];

    auto g_y_0_0_0_x_x_xz_zz = buffer_1000_ppdd[341];

    auto g_y_0_0_0_x_x_yy_xx = buffer_1000_ppdd[342];

    auto g_y_0_0_0_x_x_yy_xy = buffer_1000_ppdd[343];

    auto g_y_0_0_0_x_x_yy_xz = buffer_1000_ppdd[344];

    auto g_y_0_0_0_x_x_yy_yy = buffer_1000_ppdd[345];

    auto g_y_0_0_0_x_x_yy_yz = buffer_1000_ppdd[346];

    auto g_y_0_0_0_x_x_yy_zz = buffer_1000_ppdd[347];

    auto g_y_0_0_0_x_x_yz_xx = buffer_1000_ppdd[348];

    auto g_y_0_0_0_x_x_yz_xy = buffer_1000_ppdd[349];

    auto g_y_0_0_0_x_x_yz_xz = buffer_1000_ppdd[350];

    auto g_y_0_0_0_x_x_yz_yy = buffer_1000_ppdd[351];

    auto g_y_0_0_0_x_x_yz_yz = buffer_1000_ppdd[352];

    auto g_y_0_0_0_x_x_yz_zz = buffer_1000_ppdd[353];

    auto g_y_0_0_0_x_x_zz_xx = buffer_1000_ppdd[354];

    auto g_y_0_0_0_x_x_zz_xy = buffer_1000_ppdd[355];

    auto g_y_0_0_0_x_x_zz_xz = buffer_1000_ppdd[356];

    auto g_y_0_0_0_x_x_zz_yy = buffer_1000_ppdd[357];

    auto g_y_0_0_0_x_x_zz_yz = buffer_1000_ppdd[358];

    auto g_y_0_0_0_x_x_zz_zz = buffer_1000_ppdd[359];

    auto g_y_0_0_0_x_y_xx_xx = buffer_1000_ppdd[360];

    auto g_y_0_0_0_x_y_xx_xy = buffer_1000_ppdd[361];

    auto g_y_0_0_0_x_y_xx_xz = buffer_1000_ppdd[362];

    auto g_y_0_0_0_x_y_xx_yy = buffer_1000_ppdd[363];

    auto g_y_0_0_0_x_y_xx_yz = buffer_1000_ppdd[364];

    auto g_y_0_0_0_x_y_xx_zz = buffer_1000_ppdd[365];

    auto g_y_0_0_0_x_y_xy_xx = buffer_1000_ppdd[366];

    auto g_y_0_0_0_x_y_xy_xy = buffer_1000_ppdd[367];

    auto g_y_0_0_0_x_y_xy_xz = buffer_1000_ppdd[368];

    auto g_y_0_0_0_x_y_xy_yy = buffer_1000_ppdd[369];

    auto g_y_0_0_0_x_y_xy_yz = buffer_1000_ppdd[370];

    auto g_y_0_0_0_x_y_xy_zz = buffer_1000_ppdd[371];

    auto g_y_0_0_0_x_y_xz_xx = buffer_1000_ppdd[372];

    auto g_y_0_0_0_x_y_xz_xy = buffer_1000_ppdd[373];

    auto g_y_0_0_0_x_y_xz_xz = buffer_1000_ppdd[374];

    auto g_y_0_0_0_x_y_xz_yy = buffer_1000_ppdd[375];

    auto g_y_0_0_0_x_y_xz_yz = buffer_1000_ppdd[376];

    auto g_y_0_0_0_x_y_xz_zz = buffer_1000_ppdd[377];

    auto g_y_0_0_0_x_y_yy_xx = buffer_1000_ppdd[378];

    auto g_y_0_0_0_x_y_yy_xy = buffer_1000_ppdd[379];

    auto g_y_0_0_0_x_y_yy_xz = buffer_1000_ppdd[380];

    auto g_y_0_0_0_x_y_yy_yy = buffer_1000_ppdd[381];

    auto g_y_0_0_0_x_y_yy_yz = buffer_1000_ppdd[382];

    auto g_y_0_0_0_x_y_yy_zz = buffer_1000_ppdd[383];

    auto g_y_0_0_0_x_y_yz_xx = buffer_1000_ppdd[384];

    auto g_y_0_0_0_x_y_yz_xy = buffer_1000_ppdd[385];

    auto g_y_0_0_0_x_y_yz_xz = buffer_1000_ppdd[386];

    auto g_y_0_0_0_x_y_yz_yy = buffer_1000_ppdd[387];

    auto g_y_0_0_0_x_y_yz_yz = buffer_1000_ppdd[388];

    auto g_y_0_0_0_x_y_yz_zz = buffer_1000_ppdd[389];

    auto g_y_0_0_0_x_y_zz_xx = buffer_1000_ppdd[390];

    auto g_y_0_0_0_x_y_zz_xy = buffer_1000_ppdd[391];

    auto g_y_0_0_0_x_y_zz_xz = buffer_1000_ppdd[392];

    auto g_y_0_0_0_x_y_zz_yy = buffer_1000_ppdd[393];

    auto g_y_0_0_0_x_y_zz_yz = buffer_1000_ppdd[394];

    auto g_y_0_0_0_x_y_zz_zz = buffer_1000_ppdd[395];

    auto g_y_0_0_0_x_z_xx_xx = buffer_1000_ppdd[396];

    auto g_y_0_0_0_x_z_xx_xy = buffer_1000_ppdd[397];

    auto g_y_0_0_0_x_z_xx_xz = buffer_1000_ppdd[398];

    auto g_y_0_0_0_x_z_xx_yy = buffer_1000_ppdd[399];

    auto g_y_0_0_0_x_z_xx_yz = buffer_1000_ppdd[400];

    auto g_y_0_0_0_x_z_xx_zz = buffer_1000_ppdd[401];

    auto g_y_0_0_0_x_z_xy_xx = buffer_1000_ppdd[402];

    auto g_y_0_0_0_x_z_xy_xy = buffer_1000_ppdd[403];

    auto g_y_0_0_0_x_z_xy_xz = buffer_1000_ppdd[404];

    auto g_y_0_0_0_x_z_xy_yy = buffer_1000_ppdd[405];

    auto g_y_0_0_0_x_z_xy_yz = buffer_1000_ppdd[406];

    auto g_y_0_0_0_x_z_xy_zz = buffer_1000_ppdd[407];

    auto g_y_0_0_0_x_z_xz_xx = buffer_1000_ppdd[408];

    auto g_y_0_0_0_x_z_xz_xy = buffer_1000_ppdd[409];

    auto g_y_0_0_0_x_z_xz_xz = buffer_1000_ppdd[410];

    auto g_y_0_0_0_x_z_xz_yy = buffer_1000_ppdd[411];

    auto g_y_0_0_0_x_z_xz_yz = buffer_1000_ppdd[412];

    auto g_y_0_0_0_x_z_xz_zz = buffer_1000_ppdd[413];

    auto g_y_0_0_0_x_z_yy_xx = buffer_1000_ppdd[414];

    auto g_y_0_0_0_x_z_yy_xy = buffer_1000_ppdd[415];

    auto g_y_0_0_0_x_z_yy_xz = buffer_1000_ppdd[416];

    auto g_y_0_0_0_x_z_yy_yy = buffer_1000_ppdd[417];

    auto g_y_0_0_0_x_z_yy_yz = buffer_1000_ppdd[418];

    auto g_y_0_0_0_x_z_yy_zz = buffer_1000_ppdd[419];

    auto g_y_0_0_0_x_z_yz_xx = buffer_1000_ppdd[420];

    auto g_y_0_0_0_x_z_yz_xy = buffer_1000_ppdd[421];

    auto g_y_0_0_0_x_z_yz_xz = buffer_1000_ppdd[422];

    auto g_y_0_0_0_x_z_yz_yy = buffer_1000_ppdd[423];

    auto g_y_0_0_0_x_z_yz_yz = buffer_1000_ppdd[424];

    auto g_y_0_0_0_x_z_yz_zz = buffer_1000_ppdd[425];

    auto g_y_0_0_0_x_z_zz_xx = buffer_1000_ppdd[426];

    auto g_y_0_0_0_x_z_zz_xy = buffer_1000_ppdd[427];

    auto g_y_0_0_0_x_z_zz_xz = buffer_1000_ppdd[428];

    auto g_y_0_0_0_x_z_zz_yy = buffer_1000_ppdd[429];

    auto g_y_0_0_0_x_z_zz_yz = buffer_1000_ppdd[430];

    auto g_y_0_0_0_x_z_zz_zz = buffer_1000_ppdd[431];

    auto g_y_0_0_0_y_x_xx_xx = buffer_1000_ppdd[432];

    auto g_y_0_0_0_y_x_xx_xy = buffer_1000_ppdd[433];

    auto g_y_0_0_0_y_x_xx_xz = buffer_1000_ppdd[434];

    auto g_y_0_0_0_y_x_xx_yy = buffer_1000_ppdd[435];

    auto g_y_0_0_0_y_x_xx_yz = buffer_1000_ppdd[436];

    auto g_y_0_0_0_y_x_xx_zz = buffer_1000_ppdd[437];

    auto g_y_0_0_0_y_x_xy_xx = buffer_1000_ppdd[438];

    auto g_y_0_0_0_y_x_xy_xy = buffer_1000_ppdd[439];

    auto g_y_0_0_0_y_x_xy_xz = buffer_1000_ppdd[440];

    auto g_y_0_0_0_y_x_xy_yy = buffer_1000_ppdd[441];

    auto g_y_0_0_0_y_x_xy_yz = buffer_1000_ppdd[442];

    auto g_y_0_0_0_y_x_xy_zz = buffer_1000_ppdd[443];

    auto g_y_0_0_0_y_x_xz_xx = buffer_1000_ppdd[444];

    auto g_y_0_0_0_y_x_xz_xy = buffer_1000_ppdd[445];

    auto g_y_0_0_0_y_x_xz_xz = buffer_1000_ppdd[446];

    auto g_y_0_0_0_y_x_xz_yy = buffer_1000_ppdd[447];

    auto g_y_0_0_0_y_x_xz_yz = buffer_1000_ppdd[448];

    auto g_y_0_0_0_y_x_xz_zz = buffer_1000_ppdd[449];

    auto g_y_0_0_0_y_x_yy_xx = buffer_1000_ppdd[450];

    auto g_y_0_0_0_y_x_yy_xy = buffer_1000_ppdd[451];

    auto g_y_0_0_0_y_x_yy_xz = buffer_1000_ppdd[452];

    auto g_y_0_0_0_y_x_yy_yy = buffer_1000_ppdd[453];

    auto g_y_0_0_0_y_x_yy_yz = buffer_1000_ppdd[454];

    auto g_y_0_0_0_y_x_yy_zz = buffer_1000_ppdd[455];

    auto g_y_0_0_0_y_x_yz_xx = buffer_1000_ppdd[456];

    auto g_y_0_0_0_y_x_yz_xy = buffer_1000_ppdd[457];

    auto g_y_0_0_0_y_x_yz_xz = buffer_1000_ppdd[458];

    auto g_y_0_0_0_y_x_yz_yy = buffer_1000_ppdd[459];

    auto g_y_0_0_0_y_x_yz_yz = buffer_1000_ppdd[460];

    auto g_y_0_0_0_y_x_yz_zz = buffer_1000_ppdd[461];

    auto g_y_0_0_0_y_x_zz_xx = buffer_1000_ppdd[462];

    auto g_y_0_0_0_y_x_zz_xy = buffer_1000_ppdd[463];

    auto g_y_0_0_0_y_x_zz_xz = buffer_1000_ppdd[464];

    auto g_y_0_0_0_y_x_zz_yy = buffer_1000_ppdd[465];

    auto g_y_0_0_0_y_x_zz_yz = buffer_1000_ppdd[466];

    auto g_y_0_0_0_y_x_zz_zz = buffer_1000_ppdd[467];

    auto g_y_0_0_0_y_y_xx_xx = buffer_1000_ppdd[468];

    auto g_y_0_0_0_y_y_xx_xy = buffer_1000_ppdd[469];

    auto g_y_0_0_0_y_y_xx_xz = buffer_1000_ppdd[470];

    auto g_y_0_0_0_y_y_xx_yy = buffer_1000_ppdd[471];

    auto g_y_0_0_0_y_y_xx_yz = buffer_1000_ppdd[472];

    auto g_y_0_0_0_y_y_xx_zz = buffer_1000_ppdd[473];

    auto g_y_0_0_0_y_y_xy_xx = buffer_1000_ppdd[474];

    auto g_y_0_0_0_y_y_xy_xy = buffer_1000_ppdd[475];

    auto g_y_0_0_0_y_y_xy_xz = buffer_1000_ppdd[476];

    auto g_y_0_0_0_y_y_xy_yy = buffer_1000_ppdd[477];

    auto g_y_0_0_0_y_y_xy_yz = buffer_1000_ppdd[478];

    auto g_y_0_0_0_y_y_xy_zz = buffer_1000_ppdd[479];

    auto g_y_0_0_0_y_y_xz_xx = buffer_1000_ppdd[480];

    auto g_y_0_0_0_y_y_xz_xy = buffer_1000_ppdd[481];

    auto g_y_0_0_0_y_y_xz_xz = buffer_1000_ppdd[482];

    auto g_y_0_0_0_y_y_xz_yy = buffer_1000_ppdd[483];

    auto g_y_0_0_0_y_y_xz_yz = buffer_1000_ppdd[484];

    auto g_y_0_0_0_y_y_xz_zz = buffer_1000_ppdd[485];

    auto g_y_0_0_0_y_y_yy_xx = buffer_1000_ppdd[486];

    auto g_y_0_0_0_y_y_yy_xy = buffer_1000_ppdd[487];

    auto g_y_0_0_0_y_y_yy_xz = buffer_1000_ppdd[488];

    auto g_y_0_0_0_y_y_yy_yy = buffer_1000_ppdd[489];

    auto g_y_0_0_0_y_y_yy_yz = buffer_1000_ppdd[490];

    auto g_y_0_0_0_y_y_yy_zz = buffer_1000_ppdd[491];

    auto g_y_0_0_0_y_y_yz_xx = buffer_1000_ppdd[492];

    auto g_y_0_0_0_y_y_yz_xy = buffer_1000_ppdd[493];

    auto g_y_0_0_0_y_y_yz_xz = buffer_1000_ppdd[494];

    auto g_y_0_0_0_y_y_yz_yy = buffer_1000_ppdd[495];

    auto g_y_0_0_0_y_y_yz_yz = buffer_1000_ppdd[496];

    auto g_y_0_0_0_y_y_yz_zz = buffer_1000_ppdd[497];

    auto g_y_0_0_0_y_y_zz_xx = buffer_1000_ppdd[498];

    auto g_y_0_0_0_y_y_zz_xy = buffer_1000_ppdd[499];

    auto g_y_0_0_0_y_y_zz_xz = buffer_1000_ppdd[500];

    auto g_y_0_0_0_y_y_zz_yy = buffer_1000_ppdd[501];

    auto g_y_0_0_0_y_y_zz_yz = buffer_1000_ppdd[502];

    auto g_y_0_0_0_y_y_zz_zz = buffer_1000_ppdd[503];

    auto g_y_0_0_0_y_z_xx_xx = buffer_1000_ppdd[504];

    auto g_y_0_0_0_y_z_xx_xy = buffer_1000_ppdd[505];

    auto g_y_0_0_0_y_z_xx_xz = buffer_1000_ppdd[506];

    auto g_y_0_0_0_y_z_xx_yy = buffer_1000_ppdd[507];

    auto g_y_0_0_0_y_z_xx_yz = buffer_1000_ppdd[508];

    auto g_y_0_0_0_y_z_xx_zz = buffer_1000_ppdd[509];

    auto g_y_0_0_0_y_z_xy_xx = buffer_1000_ppdd[510];

    auto g_y_0_0_0_y_z_xy_xy = buffer_1000_ppdd[511];

    auto g_y_0_0_0_y_z_xy_xz = buffer_1000_ppdd[512];

    auto g_y_0_0_0_y_z_xy_yy = buffer_1000_ppdd[513];

    auto g_y_0_0_0_y_z_xy_yz = buffer_1000_ppdd[514];

    auto g_y_0_0_0_y_z_xy_zz = buffer_1000_ppdd[515];

    auto g_y_0_0_0_y_z_xz_xx = buffer_1000_ppdd[516];

    auto g_y_0_0_0_y_z_xz_xy = buffer_1000_ppdd[517];

    auto g_y_0_0_0_y_z_xz_xz = buffer_1000_ppdd[518];

    auto g_y_0_0_0_y_z_xz_yy = buffer_1000_ppdd[519];

    auto g_y_0_0_0_y_z_xz_yz = buffer_1000_ppdd[520];

    auto g_y_0_0_0_y_z_xz_zz = buffer_1000_ppdd[521];

    auto g_y_0_0_0_y_z_yy_xx = buffer_1000_ppdd[522];

    auto g_y_0_0_0_y_z_yy_xy = buffer_1000_ppdd[523];

    auto g_y_0_0_0_y_z_yy_xz = buffer_1000_ppdd[524];

    auto g_y_0_0_0_y_z_yy_yy = buffer_1000_ppdd[525];

    auto g_y_0_0_0_y_z_yy_yz = buffer_1000_ppdd[526];

    auto g_y_0_0_0_y_z_yy_zz = buffer_1000_ppdd[527];

    auto g_y_0_0_0_y_z_yz_xx = buffer_1000_ppdd[528];

    auto g_y_0_0_0_y_z_yz_xy = buffer_1000_ppdd[529];

    auto g_y_0_0_0_y_z_yz_xz = buffer_1000_ppdd[530];

    auto g_y_0_0_0_y_z_yz_yy = buffer_1000_ppdd[531];

    auto g_y_0_0_0_y_z_yz_yz = buffer_1000_ppdd[532];

    auto g_y_0_0_0_y_z_yz_zz = buffer_1000_ppdd[533];

    auto g_y_0_0_0_y_z_zz_xx = buffer_1000_ppdd[534];

    auto g_y_0_0_0_y_z_zz_xy = buffer_1000_ppdd[535];

    auto g_y_0_0_0_y_z_zz_xz = buffer_1000_ppdd[536];

    auto g_y_0_0_0_y_z_zz_yy = buffer_1000_ppdd[537];

    auto g_y_0_0_0_y_z_zz_yz = buffer_1000_ppdd[538];

    auto g_y_0_0_0_y_z_zz_zz = buffer_1000_ppdd[539];

    auto g_y_0_0_0_z_x_xx_xx = buffer_1000_ppdd[540];

    auto g_y_0_0_0_z_x_xx_xy = buffer_1000_ppdd[541];

    auto g_y_0_0_0_z_x_xx_xz = buffer_1000_ppdd[542];

    auto g_y_0_0_0_z_x_xx_yy = buffer_1000_ppdd[543];

    auto g_y_0_0_0_z_x_xx_yz = buffer_1000_ppdd[544];

    auto g_y_0_0_0_z_x_xx_zz = buffer_1000_ppdd[545];

    auto g_y_0_0_0_z_x_xy_xx = buffer_1000_ppdd[546];

    auto g_y_0_0_0_z_x_xy_xy = buffer_1000_ppdd[547];

    auto g_y_0_0_0_z_x_xy_xz = buffer_1000_ppdd[548];

    auto g_y_0_0_0_z_x_xy_yy = buffer_1000_ppdd[549];

    auto g_y_0_0_0_z_x_xy_yz = buffer_1000_ppdd[550];

    auto g_y_0_0_0_z_x_xy_zz = buffer_1000_ppdd[551];

    auto g_y_0_0_0_z_x_xz_xx = buffer_1000_ppdd[552];

    auto g_y_0_0_0_z_x_xz_xy = buffer_1000_ppdd[553];

    auto g_y_0_0_0_z_x_xz_xz = buffer_1000_ppdd[554];

    auto g_y_0_0_0_z_x_xz_yy = buffer_1000_ppdd[555];

    auto g_y_0_0_0_z_x_xz_yz = buffer_1000_ppdd[556];

    auto g_y_0_0_0_z_x_xz_zz = buffer_1000_ppdd[557];

    auto g_y_0_0_0_z_x_yy_xx = buffer_1000_ppdd[558];

    auto g_y_0_0_0_z_x_yy_xy = buffer_1000_ppdd[559];

    auto g_y_0_0_0_z_x_yy_xz = buffer_1000_ppdd[560];

    auto g_y_0_0_0_z_x_yy_yy = buffer_1000_ppdd[561];

    auto g_y_0_0_0_z_x_yy_yz = buffer_1000_ppdd[562];

    auto g_y_0_0_0_z_x_yy_zz = buffer_1000_ppdd[563];

    auto g_y_0_0_0_z_x_yz_xx = buffer_1000_ppdd[564];

    auto g_y_0_0_0_z_x_yz_xy = buffer_1000_ppdd[565];

    auto g_y_0_0_0_z_x_yz_xz = buffer_1000_ppdd[566];

    auto g_y_0_0_0_z_x_yz_yy = buffer_1000_ppdd[567];

    auto g_y_0_0_0_z_x_yz_yz = buffer_1000_ppdd[568];

    auto g_y_0_0_0_z_x_yz_zz = buffer_1000_ppdd[569];

    auto g_y_0_0_0_z_x_zz_xx = buffer_1000_ppdd[570];

    auto g_y_0_0_0_z_x_zz_xy = buffer_1000_ppdd[571];

    auto g_y_0_0_0_z_x_zz_xz = buffer_1000_ppdd[572];

    auto g_y_0_0_0_z_x_zz_yy = buffer_1000_ppdd[573];

    auto g_y_0_0_0_z_x_zz_yz = buffer_1000_ppdd[574];

    auto g_y_0_0_0_z_x_zz_zz = buffer_1000_ppdd[575];

    auto g_y_0_0_0_z_y_xx_xx = buffer_1000_ppdd[576];

    auto g_y_0_0_0_z_y_xx_xy = buffer_1000_ppdd[577];

    auto g_y_0_0_0_z_y_xx_xz = buffer_1000_ppdd[578];

    auto g_y_0_0_0_z_y_xx_yy = buffer_1000_ppdd[579];

    auto g_y_0_0_0_z_y_xx_yz = buffer_1000_ppdd[580];

    auto g_y_0_0_0_z_y_xx_zz = buffer_1000_ppdd[581];

    auto g_y_0_0_0_z_y_xy_xx = buffer_1000_ppdd[582];

    auto g_y_0_0_0_z_y_xy_xy = buffer_1000_ppdd[583];

    auto g_y_0_0_0_z_y_xy_xz = buffer_1000_ppdd[584];

    auto g_y_0_0_0_z_y_xy_yy = buffer_1000_ppdd[585];

    auto g_y_0_0_0_z_y_xy_yz = buffer_1000_ppdd[586];

    auto g_y_0_0_0_z_y_xy_zz = buffer_1000_ppdd[587];

    auto g_y_0_0_0_z_y_xz_xx = buffer_1000_ppdd[588];

    auto g_y_0_0_0_z_y_xz_xy = buffer_1000_ppdd[589];

    auto g_y_0_0_0_z_y_xz_xz = buffer_1000_ppdd[590];

    auto g_y_0_0_0_z_y_xz_yy = buffer_1000_ppdd[591];

    auto g_y_0_0_0_z_y_xz_yz = buffer_1000_ppdd[592];

    auto g_y_0_0_0_z_y_xz_zz = buffer_1000_ppdd[593];

    auto g_y_0_0_0_z_y_yy_xx = buffer_1000_ppdd[594];

    auto g_y_0_0_0_z_y_yy_xy = buffer_1000_ppdd[595];

    auto g_y_0_0_0_z_y_yy_xz = buffer_1000_ppdd[596];

    auto g_y_0_0_0_z_y_yy_yy = buffer_1000_ppdd[597];

    auto g_y_0_0_0_z_y_yy_yz = buffer_1000_ppdd[598];

    auto g_y_0_0_0_z_y_yy_zz = buffer_1000_ppdd[599];

    auto g_y_0_0_0_z_y_yz_xx = buffer_1000_ppdd[600];

    auto g_y_0_0_0_z_y_yz_xy = buffer_1000_ppdd[601];

    auto g_y_0_0_0_z_y_yz_xz = buffer_1000_ppdd[602];

    auto g_y_0_0_0_z_y_yz_yy = buffer_1000_ppdd[603];

    auto g_y_0_0_0_z_y_yz_yz = buffer_1000_ppdd[604];

    auto g_y_0_0_0_z_y_yz_zz = buffer_1000_ppdd[605];

    auto g_y_0_0_0_z_y_zz_xx = buffer_1000_ppdd[606];

    auto g_y_0_0_0_z_y_zz_xy = buffer_1000_ppdd[607];

    auto g_y_0_0_0_z_y_zz_xz = buffer_1000_ppdd[608];

    auto g_y_0_0_0_z_y_zz_yy = buffer_1000_ppdd[609];

    auto g_y_0_0_0_z_y_zz_yz = buffer_1000_ppdd[610];

    auto g_y_0_0_0_z_y_zz_zz = buffer_1000_ppdd[611];

    auto g_y_0_0_0_z_z_xx_xx = buffer_1000_ppdd[612];

    auto g_y_0_0_0_z_z_xx_xy = buffer_1000_ppdd[613];

    auto g_y_0_0_0_z_z_xx_xz = buffer_1000_ppdd[614];

    auto g_y_0_0_0_z_z_xx_yy = buffer_1000_ppdd[615];

    auto g_y_0_0_0_z_z_xx_yz = buffer_1000_ppdd[616];

    auto g_y_0_0_0_z_z_xx_zz = buffer_1000_ppdd[617];

    auto g_y_0_0_0_z_z_xy_xx = buffer_1000_ppdd[618];

    auto g_y_0_0_0_z_z_xy_xy = buffer_1000_ppdd[619];

    auto g_y_0_0_0_z_z_xy_xz = buffer_1000_ppdd[620];

    auto g_y_0_0_0_z_z_xy_yy = buffer_1000_ppdd[621];

    auto g_y_0_0_0_z_z_xy_yz = buffer_1000_ppdd[622];

    auto g_y_0_0_0_z_z_xy_zz = buffer_1000_ppdd[623];

    auto g_y_0_0_0_z_z_xz_xx = buffer_1000_ppdd[624];

    auto g_y_0_0_0_z_z_xz_xy = buffer_1000_ppdd[625];

    auto g_y_0_0_0_z_z_xz_xz = buffer_1000_ppdd[626];

    auto g_y_0_0_0_z_z_xz_yy = buffer_1000_ppdd[627];

    auto g_y_0_0_0_z_z_xz_yz = buffer_1000_ppdd[628];

    auto g_y_0_0_0_z_z_xz_zz = buffer_1000_ppdd[629];

    auto g_y_0_0_0_z_z_yy_xx = buffer_1000_ppdd[630];

    auto g_y_0_0_0_z_z_yy_xy = buffer_1000_ppdd[631];

    auto g_y_0_0_0_z_z_yy_xz = buffer_1000_ppdd[632];

    auto g_y_0_0_0_z_z_yy_yy = buffer_1000_ppdd[633];

    auto g_y_0_0_0_z_z_yy_yz = buffer_1000_ppdd[634];

    auto g_y_0_0_0_z_z_yy_zz = buffer_1000_ppdd[635];

    auto g_y_0_0_0_z_z_yz_xx = buffer_1000_ppdd[636];

    auto g_y_0_0_0_z_z_yz_xy = buffer_1000_ppdd[637];

    auto g_y_0_0_0_z_z_yz_xz = buffer_1000_ppdd[638];

    auto g_y_0_0_0_z_z_yz_yy = buffer_1000_ppdd[639];

    auto g_y_0_0_0_z_z_yz_yz = buffer_1000_ppdd[640];

    auto g_y_0_0_0_z_z_yz_zz = buffer_1000_ppdd[641];

    auto g_y_0_0_0_z_z_zz_xx = buffer_1000_ppdd[642];

    auto g_y_0_0_0_z_z_zz_xy = buffer_1000_ppdd[643];

    auto g_y_0_0_0_z_z_zz_xz = buffer_1000_ppdd[644];

    auto g_y_0_0_0_z_z_zz_yy = buffer_1000_ppdd[645];

    auto g_y_0_0_0_z_z_zz_yz = buffer_1000_ppdd[646];

    auto g_y_0_0_0_z_z_zz_zz = buffer_1000_ppdd[647];

    auto g_z_0_0_0_x_x_xx_xx = buffer_1000_ppdd[648];

    auto g_z_0_0_0_x_x_xx_xy = buffer_1000_ppdd[649];

    auto g_z_0_0_0_x_x_xx_xz = buffer_1000_ppdd[650];

    auto g_z_0_0_0_x_x_xx_yy = buffer_1000_ppdd[651];

    auto g_z_0_0_0_x_x_xx_yz = buffer_1000_ppdd[652];

    auto g_z_0_0_0_x_x_xx_zz = buffer_1000_ppdd[653];

    auto g_z_0_0_0_x_x_xy_xx = buffer_1000_ppdd[654];

    auto g_z_0_0_0_x_x_xy_xy = buffer_1000_ppdd[655];

    auto g_z_0_0_0_x_x_xy_xz = buffer_1000_ppdd[656];

    auto g_z_0_0_0_x_x_xy_yy = buffer_1000_ppdd[657];

    auto g_z_0_0_0_x_x_xy_yz = buffer_1000_ppdd[658];

    auto g_z_0_0_0_x_x_xy_zz = buffer_1000_ppdd[659];

    auto g_z_0_0_0_x_x_xz_xx = buffer_1000_ppdd[660];

    auto g_z_0_0_0_x_x_xz_xy = buffer_1000_ppdd[661];

    auto g_z_0_0_0_x_x_xz_xz = buffer_1000_ppdd[662];

    auto g_z_0_0_0_x_x_xz_yy = buffer_1000_ppdd[663];

    auto g_z_0_0_0_x_x_xz_yz = buffer_1000_ppdd[664];

    auto g_z_0_0_0_x_x_xz_zz = buffer_1000_ppdd[665];

    auto g_z_0_0_0_x_x_yy_xx = buffer_1000_ppdd[666];

    auto g_z_0_0_0_x_x_yy_xy = buffer_1000_ppdd[667];

    auto g_z_0_0_0_x_x_yy_xz = buffer_1000_ppdd[668];

    auto g_z_0_0_0_x_x_yy_yy = buffer_1000_ppdd[669];

    auto g_z_0_0_0_x_x_yy_yz = buffer_1000_ppdd[670];

    auto g_z_0_0_0_x_x_yy_zz = buffer_1000_ppdd[671];

    auto g_z_0_0_0_x_x_yz_xx = buffer_1000_ppdd[672];

    auto g_z_0_0_0_x_x_yz_xy = buffer_1000_ppdd[673];

    auto g_z_0_0_0_x_x_yz_xz = buffer_1000_ppdd[674];

    auto g_z_0_0_0_x_x_yz_yy = buffer_1000_ppdd[675];

    auto g_z_0_0_0_x_x_yz_yz = buffer_1000_ppdd[676];

    auto g_z_0_0_0_x_x_yz_zz = buffer_1000_ppdd[677];

    auto g_z_0_0_0_x_x_zz_xx = buffer_1000_ppdd[678];

    auto g_z_0_0_0_x_x_zz_xy = buffer_1000_ppdd[679];

    auto g_z_0_0_0_x_x_zz_xz = buffer_1000_ppdd[680];

    auto g_z_0_0_0_x_x_zz_yy = buffer_1000_ppdd[681];

    auto g_z_0_0_0_x_x_zz_yz = buffer_1000_ppdd[682];

    auto g_z_0_0_0_x_x_zz_zz = buffer_1000_ppdd[683];

    auto g_z_0_0_0_x_y_xx_xx = buffer_1000_ppdd[684];

    auto g_z_0_0_0_x_y_xx_xy = buffer_1000_ppdd[685];

    auto g_z_0_0_0_x_y_xx_xz = buffer_1000_ppdd[686];

    auto g_z_0_0_0_x_y_xx_yy = buffer_1000_ppdd[687];

    auto g_z_0_0_0_x_y_xx_yz = buffer_1000_ppdd[688];

    auto g_z_0_0_0_x_y_xx_zz = buffer_1000_ppdd[689];

    auto g_z_0_0_0_x_y_xy_xx = buffer_1000_ppdd[690];

    auto g_z_0_0_0_x_y_xy_xy = buffer_1000_ppdd[691];

    auto g_z_0_0_0_x_y_xy_xz = buffer_1000_ppdd[692];

    auto g_z_0_0_0_x_y_xy_yy = buffer_1000_ppdd[693];

    auto g_z_0_0_0_x_y_xy_yz = buffer_1000_ppdd[694];

    auto g_z_0_0_0_x_y_xy_zz = buffer_1000_ppdd[695];

    auto g_z_0_0_0_x_y_xz_xx = buffer_1000_ppdd[696];

    auto g_z_0_0_0_x_y_xz_xy = buffer_1000_ppdd[697];

    auto g_z_0_0_0_x_y_xz_xz = buffer_1000_ppdd[698];

    auto g_z_0_0_0_x_y_xz_yy = buffer_1000_ppdd[699];

    auto g_z_0_0_0_x_y_xz_yz = buffer_1000_ppdd[700];

    auto g_z_0_0_0_x_y_xz_zz = buffer_1000_ppdd[701];

    auto g_z_0_0_0_x_y_yy_xx = buffer_1000_ppdd[702];

    auto g_z_0_0_0_x_y_yy_xy = buffer_1000_ppdd[703];

    auto g_z_0_0_0_x_y_yy_xz = buffer_1000_ppdd[704];

    auto g_z_0_0_0_x_y_yy_yy = buffer_1000_ppdd[705];

    auto g_z_0_0_0_x_y_yy_yz = buffer_1000_ppdd[706];

    auto g_z_0_0_0_x_y_yy_zz = buffer_1000_ppdd[707];

    auto g_z_0_0_0_x_y_yz_xx = buffer_1000_ppdd[708];

    auto g_z_0_0_0_x_y_yz_xy = buffer_1000_ppdd[709];

    auto g_z_0_0_0_x_y_yz_xz = buffer_1000_ppdd[710];

    auto g_z_0_0_0_x_y_yz_yy = buffer_1000_ppdd[711];

    auto g_z_0_0_0_x_y_yz_yz = buffer_1000_ppdd[712];

    auto g_z_0_0_0_x_y_yz_zz = buffer_1000_ppdd[713];

    auto g_z_0_0_0_x_y_zz_xx = buffer_1000_ppdd[714];

    auto g_z_0_0_0_x_y_zz_xy = buffer_1000_ppdd[715];

    auto g_z_0_0_0_x_y_zz_xz = buffer_1000_ppdd[716];

    auto g_z_0_0_0_x_y_zz_yy = buffer_1000_ppdd[717];

    auto g_z_0_0_0_x_y_zz_yz = buffer_1000_ppdd[718];

    auto g_z_0_0_0_x_y_zz_zz = buffer_1000_ppdd[719];

    auto g_z_0_0_0_x_z_xx_xx = buffer_1000_ppdd[720];

    auto g_z_0_0_0_x_z_xx_xy = buffer_1000_ppdd[721];

    auto g_z_0_0_0_x_z_xx_xz = buffer_1000_ppdd[722];

    auto g_z_0_0_0_x_z_xx_yy = buffer_1000_ppdd[723];

    auto g_z_0_0_0_x_z_xx_yz = buffer_1000_ppdd[724];

    auto g_z_0_0_0_x_z_xx_zz = buffer_1000_ppdd[725];

    auto g_z_0_0_0_x_z_xy_xx = buffer_1000_ppdd[726];

    auto g_z_0_0_0_x_z_xy_xy = buffer_1000_ppdd[727];

    auto g_z_0_0_0_x_z_xy_xz = buffer_1000_ppdd[728];

    auto g_z_0_0_0_x_z_xy_yy = buffer_1000_ppdd[729];

    auto g_z_0_0_0_x_z_xy_yz = buffer_1000_ppdd[730];

    auto g_z_0_0_0_x_z_xy_zz = buffer_1000_ppdd[731];

    auto g_z_0_0_0_x_z_xz_xx = buffer_1000_ppdd[732];

    auto g_z_0_0_0_x_z_xz_xy = buffer_1000_ppdd[733];

    auto g_z_0_0_0_x_z_xz_xz = buffer_1000_ppdd[734];

    auto g_z_0_0_0_x_z_xz_yy = buffer_1000_ppdd[735];

    auto g_z_0_0_0_x_z_xz_yz = buffer_1000_ppdd[736];

    auto g_z_0_0_0_x_z_xz_zz = buffer_1000_ppdd[737];

    auto g_z_0_0_0_x_z_yy_xx = buffer_1000_ppdd[738];

    auto g_z_0_0_0_x_z_yy_xy = buffer_1000_ppdd[739];

    auto g_z_0_0_0_x_z_yy_xz = buffer_1000_ppdd[740];

    auto g_z_0_0_0_x_z_yy_yy = buffer_1000_ppdd[741];

    auto g_z_0_0_0_x_z_yy_yz = buffer_1000_ppdd[742];

    auto g_z_0_0_0_x_z_yy_zz = buffer_1000_ppdd[743];

    auto g_z_0_0_0_x_z_yz_xx = buffer_1000_ppdd[744];

    auto g_z_0_0_0_x_z_yz_xy = buffer_1000_ppdd[745];

    auto g_z_0_0_0_x_z_yz_xz = buffer_1000_ppdd[746];

    auto g_z_0_0_0_x_z_yz_yy = buffer_1000_ppdd[747];

    auto g_z_0_0_0_x_z_yz_yz = buffer_1000_ppdd[748];

    auto g_z_0_0_0_x_z_yz_zz = buffer_1000_ppdd[749];

    auto g_z_0_0_0_x_z_zz_xx = buffer_1000_ppdd[750];

    auto g_z_0_0_0_x_z_zz_xy = buffer_1000_ppdd[751];

    auto g_z_0_0_0_x_z_zz_xz = buffer_1000_ppdd[752];

    auto g_z_0_0_0_x_z_zz_yy = buffer_1000_ppdd[753];

    auto g_z_0_0_0_x_z_zz_yz = buffer_1000_ppdd[754];

    auto g_z_0_0_0_x_z_zz_zz = buffer_1000_ppdd[755];

    auto g_z_0_0_0_y_x_xx_xx = buffer_1000_ppdd[756];

    auto g_z_0_0_0_y_x_xx_xy = buffer_1000_ppdd[757];

    auto g_z_0_0_0_y_x_xx_xz = buffer_1000_ppdd[758];

    auto g_z_0_0_0_y_x_xx_yy = buffer_1000_ppdd[759];

    auto g_z_0_0_0_y_x_xx_yz = buffer_1000_ppdd[760];

    auto g_z_0_0_0_y_x_xx_zz = buffer_1000_ppdd[761];

    auto g_z_0_0_0_y_x_xy_xx = buffer_1000_ppdd[762];

    auto g_z_0_0_0_y_x_xy_xy = buffer_1000_ppdd[763];

    auto g_z_0_0_0_y_x_xy_xz = buffer_1000_ppdd[764];

    auto g_z_0_0_0_y_x_xy_yy = buffer_1000_ppdd[765];

    auto g_z_0_0_0_y_x_xy_yz = buffer_1000_ppdd[766];

    auto g_z_0_0_0_y_x_xy_zz = buffer_1000_ppdd[767];

    auto g_z_0_0_0_y_x_xz_xx = buffer_1000_ppdd[768];

    auto g_z_0_0_0_y_x_xz_xy = buffer_1000_ppdd[769];

    auto g_z_0_0_0_y_x_xz_xz = buffer_1000_ppdd[770];

    auto g_z_0_0_0_y_x_xz_yy = buffer_1000_ppdd[771];

    auto g_z_0_0_0_y_x_xz_yz = buffer_1000_ppdd[772];

    auto g_z_0_0_0_y_x_xz_zz = buffer_1000_ppdd[773];

    auto g_z_0_0_0_y_x_yy_xx = buffer_1000_ppdd[774];

    auto g_z_0_0_0_y_x_yy_xy = buffer_1000_ppdd[775];

    auto g_z_0_0_0_y_x_yy_xz = buffer_1000_ppdd[776];

    auto g_z_0_0_0_y_x_yy_yy = buffer_1000_ppdd[777];

    auto g_z_0_0_0_y_x_yy_yz = buffer_1000_ppdd[778];

    auto g_z_0_0_0_y_x_yy_zz = buffer_1000_ppdd[779];

    auto g_z_0_0_0_y_x_yz_xx = buffer_1000_ppdd[780];

    auto g_z_0_0_0_y_x_yz_xy = buffer_1000_ppdd[781];

    auto g_z_0_0_0_y_x_yz_xz = buffer_1000_ppdd[782];

    auto g_z_0_0_0_y_x_yz_yy = buffer_1000_ppdd[783];

    auto g_z_0_0_0_y_x_yz_yz = buffer_1000_ppdd[784];

    auto g_z_0_0_0_y_x_yz_zz = buffer_1000_ppdd[785];

    auto g_z_0_0_0_y_x_zz_xx = buffer_1000_ppdd[786];

    auto g_z_0_0_0_y_x_zz_xy = buffer_1000_ppdd[787];

    auto g_z_0_0_0_y_x_zz_xz = buffer_1000_ppdd[788];

    auto g_z_0_0_0_y_x_zz_yy = buffer_1000_ppdd[789];

    auto g_z_0_0_0_y_x_zz_yz = buffer_1000_ppdd[790];

    auto g_z_0_0_0_y_x_zz_zz = buffer_1000_ppdd[791];

    auto g_z_0_0_0_y_y_xx_xx = buffer_1000_ppdd[792];

    auto g_z_0_0_0_y_y_xx_xy = buffer_1000_ppdd[793];

    auto g_z_0_0_0_y_y_xx_xz = buffer_1000_ppdd[794];

    auto g_z_0_0_0_y_y_xx_yy = buffer_1000_ppdd[795];

    auto g_z_0_0_0_y_y_xx_yz = buffer_1000_ppdd[796];

    auto g_z_0_0_0_y_y_xx_zz = buffer_1000_ppdd[797];

    auto g_z_0_0_0_y_y_xy_xx = buffer_1000_ppdd[798];

    auto g_z_0_0_0_y_y_xy_xy = buffer_1000_ppdd[799];

    auto g_z_0_0_0_y_y_xy_xz = buffer_1000_ppdd[800];

    auto g_z_0_0_0_y_y_xy_yy = buffer_1000_ppdd[801];

    auto g_z_0_0_0_y_y_xy_yz = buffer_1000_ppdd[802];

    auto g_z_0_0_0_y_y_xy_zz = buffer_1000_ppdd[803];

    auto g_z_0_0_0_y_y_xz_xx = buffer_1000_ppdd[804];

    auto g_z_0_0_0_y_y_xz_xy = buffer_1000_ppdd[805];

    auto g_z_0_0_0_y_y_xz_xz = buffer_1000_ppdd[806];

    auto g_z_0_0_0_y_y_xz_yy = buffer_1000_ppdd[807];

    auto g_z_0_0_0_y_y_xz_yz = buffer_1000_ppdd[808];

    auto g_z_0_0_0_y_y_xz_zz = buffer_1000_ppdd[809];

    auto g_z_0_0_0_y_y_yy_xx = buffer_1000_ppdd[810];

    auto g_z_0_0_0_y_y_yy_xy = buffer_1000_ppdd[811];

    auto g_z_0_0_0_y_y_yy_xz = buffer_1000_ppdd[812];

    auto g_z_0_0_0_y_y_yy_yy = buffer_1000_ppdd[813];

    auto g_z_0_0_0_y_y_yy_yz = buffer_1000_ppdd[814];

    auto g_z_0_0_0_y_y_yy_zz = buffer_1000_ppdd[815];

    auto g_z_0_0_0_y_y_yz_xx = buffer_1000_ppdd[816];

    auto g_z_0_0_0_y_y_yz_xy = buffer_1000_ppdd[817];

    auto g_z_0_0_0_y_y_yz_xz = buffer_1000_ppdd[818];

    auto g_z_0_0_0_y_y_yz_yy = buffer_1000_ppdd[819];

    auto g_z_0_0_0_y_y_yz_yz = buffer_1000_ppdd[820];

    auto g_z_0_0_0_y_y_yz_zz = buffer_1000_ppdd[821];

    auto g_z_0_0_0_y_y_zz_xx = buffer_1000_ppdd[822];

    auto g_z_0_0_0_y_y_zz_xy = buffer_1000_ppdd[823];

    auto g_z_0_0_0_y_y_zz_xz = buffer_1000_ppdd[824];

    auto g_z_0_0_0_y_y_zz_yy = buffer_1000_ppdd[825];

    auto g_z_0_0_0_y_y_zz_yz = buffer_1000_ppdd[826];

    auto g_z_0_0_0_y_y_zz_zz = buffer_1000_ppdd[827];

    auto g_z_0_0_0_y_z_xx_xx = buffer_1000_ppdd[828];

    auto g_z_0_0_0_y_z_xx_xy = buffer_1000_ppdd[829];

    auto g_z_0_0_0_y_z_xx_xz = buffer_1000_ppdd[830];

    auto g_z_0_0_0_y_z_xx_yy = buffer_1000_ppdd[831];

    auto g_z_0_0_0_y_z_xx_yz = buffer_1000_ppdd[832];

    auto g_z_0_0_0_y_z_xx_zz = buffer_1000_ppdd[833];

    auto g_z_0_0_0_y_z_xy_xx = buffer_1000_ppdd[834];

    auto g_z_0_0_0_y_z_xy_xy = buffer_1000_ppdd[835];

    auto g_z_0_0_0_y_z_xy_xz = buffer_1000_ppdd[836];

    auto g_z_0_0_0_y_z_xy_yy = buffer_1000_ppdd[837];

    auto g_z_0_0_0_y_z_xy_yz = buffer_1000_ppdd[838];

    auto g_z_0_0_0_y_z_xy_zz = buffer_1000_ppdd[839];

    auto g_z_0_0_0_y_z_xz_xx = buffer_1000_ppdd[840];

    auto g_z_0_0_0_y_z_xz_xy = buffer_1000_ppdd[841];

    auto g_z_0_0_0_y_z_xz_xz = buffer_1000_ppdd[842];

    auto g_z_0_0_0_y_z_xz_yy = buffer_1000_ppdd[843];

    auto g_z_0_0_0_y_z_xz_yz = buffer_1000_ppdd[844];

    auto g_z_0_0_0_y_z_xz_zz = buffer_1000_ppdd[845];

    auto g_z_0_0_0_y_z_yy_xx = buffer_1000_ppdd[846];

    auto g_z_0_0_0_y_z_yy_xy = buffer_1000_ppdd[847];

    auto g_z_0_0_0_y_z_yy_xz = buffer_1000_ppdd[848];

    auto g_z_0_0_0_y_z_yy_yy = buffer_1000_ppdd[849];

    auto g_z_0_0_0_y_z_yy_yz = buffer_1000_ppdd[850];

    auto g_z_0_0_0_y_z_yy_zz = buffer_1000_ppdd[851];

    auto g_z_0_0_0_y_z_yz_xx = buffer_1000_ppdd[852];

    auto g_z_0_0_0_y_z_yz_xy = buffer_1000_ppdd[853];

    auto g_z_0_0_0_y_z_yz_xz = buffer_1000_ppdd[854];

    auto g_z_0_0_0_y_z_yz_yy = buffer_1000_ppdd[855];

    auto g_z_0_0_0_y_z_yz_yz = buffer_1000_ppdd[856];

    auto g_z_0_0_0_y_z_yz_zz = buffer_1000_ppdd[857];

    auto g_z_0_0_0_y_z_zz_xx = buffer_1000_ppdd[858];

    auto g_z_0_0_0_y_z_zz_xy = buffer_1000_ppdd[859];

    auto g_z_0_0_0_y_z_zz_xz = buffer_1000_ppdd[860];

    auto g_z_0_0_0_y_z_zz_yy = buffer_1000_ppdd[861];

    auto g_z_0_0_0_y_z_zz_yz = buffer_1000_ppdd[862];

    auto g_z_0_0_0_y_z_zz_zz = buffer_1000_ppdd[863];

    auto g_z_0_0_0_z_x_xx_xx = buffer_1000_ppdd[864];

    auto g_z_0_0_0_z_x_xx_xy = buffer_1000_ppdd[865];

    auto g_z_0_0_0_z_x_xx_xz = buffer_1000_ppdd[866];

    auto g_z_0_0_0_z_x_xx_yy = buffer_1000_ppdd[867];

    auto g_z_0_0_0_z_x_xx_yz = buffer_1000_ppdd[868];

    auto g_z_0_0_0_z_x_xx_zz = buffer_1000_ppdd[869];

    auto g_z_0_0_0_z_x_xy_xx = buffer_1000_ppdd[870];

    auto g_z_0_0_0_z_x_xy_xy = buffer_1000_ppdd[871];

    auto g_z_0_0_0_z_x_xy_xz = buffer_1000_ppdd[872];

    auto g_z_0_0_0_z_x_xy_yy = buffer_1000_ppdd[873];

    auto g_z_0_0_0_z_x_xy_yz = buffer_1000_ppdd[874];

    auto g_z_0_0_0_z_x_xy_zz = buffer_1000_ppdd[875];

    auto g_z_0_0_0_z_x_xz_xx = buffer_1000_ppdd[876];

    auto g_z_0_0_0_z_x_xz_xy = buffer_1000_ppdd[877];

    auto g_z_0_0_0_z_x_xz_xz = buffer_1000_ppdd[878];

    auto g_z_0_0_0_z_x_xz_yy = buffer_1000_ppdd[879];

    auto g_z_0_0_0_z_x_xz_yz = buffer_1000_ppdd[880];

    auto g_z_0_0_0_z_x_xz_zz = buffer_1000_ppdd[881];

    auto g_z_0_0_0_z_x_yy_xx = buffer_1000_ppdd[882];

    auto g_z_0_0_0_z_x_yy_xy = buffer_1000_ppdd[883];

    auto g_z_0_0_0_z_x_yy_xz = buffer_1000_ppdd[884];

    auto g_z_0_0_0_z_x_yy_yy = buffer_1000_ppdd[885];

    auto g_z_0_0_0_z_x_yy_yz = buffer_1000_ppdd[886];

    auto g_z_0_0_0_z_x_yy_zz = buffer_1000_ppdd[887];

    auto g_z_0_0_0_z_x_yz_xx = buffer_1000_ppdd[888];

    auto g_z_0_0_0_z_x_yz_xy = buffer_1000_ppdd[889];

    auto g_z_0_0_0_z_x_yz_xz = buffer_1000_ppdd[890];

    auto g_z_0_0_0_z_x_yz_yy = buffer_1000_ppdd[891];

    auto g_z_0_0_0_z_x_yz_yz = buffer_1000_ppdd[892];

    auto g_z_0_0_0_z_x_yz_zz = buffer_1000_ppdd[893];

    auto g_z_0_0_0_z_x_zz_xx = buffer_1000_ppdd[894];

    auto g_z_0_0_0_z_x_zz_xy = buffer_1000_ppdd[895];

    auto g_z_0_0_0_z_x_zz_xz = buffer_1000_ppdd[896];

    auto g_z_0_0_0_z_x_zz_yy = buffer_1000_ppdd[897];

    auto g_z_0_0_0_z_x_zz_yz = buffer_1000_ppdd[898];

    auto g_z_0_0_0_z_x_zz_zz = buffer_1000_ppdd[899];

    auto g_z_0_0_0_z_y_xx_xx = buffer_1000_ppdd[900];

    auto g_z_0_0_0_z_y_xx_xy = buffer_1000_ppdd[901];

    auto g_z_0_0_0_z_y_xx_xz = buffer_1000_ppdd[902];

    auto g_z_0_0_0_z_y_xx_yy = buffer_1000_ppdd[903];

    auto g_z_0_0_0_z_y_xx_yz = buffer_1000_ppdd[904];

    auto g_z_0_0_0_z_y_xx_zz = buffer_1000_ppdd[905];

    auto g_z_0_0_0_z_y_xy_xx = buffer_1000_ppdd[906];

    auto g_z_0_0_0_z_y_xy_xy = buffer_1000_ppdd[907];

    auto g_z_0_0_0_z_y_xy_xz = buffer_1000_ppdd[908];

    auto g_z_0_0_0_z_y_xy_yy = buffer_1000_ppdd[909];

    auto g_z_0_0_0_z_y_xy_yz = buffer_1000_ppdd[910];

    auto g_z_0_0_0_z_y_xy_zz = buffer_1000_ppdd[911];

    auto g_z_0_0_0_z_y_xz_xx = buffer_1000_ppdd[912];

    auto g_z_0_0_0_z_y_xz_xy = buffer_1000_ppdd[913];

    auto g_z_0_0_0_z_y_xz_xz = buffer_1000_ppdd[914];

    auto g_z_0_0_0_z_y_xz_yy = buffer_1000_ppdd[915];

    auto g_z_0_0_0_z_y_xz_yz = buffer_1000_ppdd[916];

    auto g_z_0_0_0_z_y_xz_zz = buffer_1000_ppdd[917];

    auto g_z_0_0_0_z_y_yy_xx = buffer_1000_ppdd[918];

    auto g_z_0_0_0_z_y_yy_xy = buffer_1000_ppdd[919];

    auto g_z_0_0_0_z_y_yy_xz = buffer_1000_ppdd[920];

    auto g_z_0_0_0_z_y_yy_yy = buffer_1000_ppdd[921];

    auto g_z_0_0_0_z_y_yy_yz = buffer_1000_ppdd[922];

    auto g_z_0_0_0_z_y_yy_zz = buffer_1000_ppdd[923];

    auto g_z_0_0_0_z_y_yz_xx = buffer_1000_ppdd[924];

    auto g_z_0_0_0_z_y_yz_xy = buffer_1000_ppdd[925];

    auto g_z_0_0_0_z_y_yz_xz = buffer_1000_ppdd[926];

    auto g_z_0_0_0_z_y_yz_yy = buffer_1000_ppdd[927];

    auto g_z_0_0_0_z_y_yz_yz = buffer_1000_ppdd[928];

    auto g_z_0_0_0_z_y_yz_zz = buffer_1000_ppdd[929];

    auto g_z_0_0_0_z_y_zz_xx = buffer_1000_ppdd[930];

    auto g_z_0_0_0_z_y_zz_xy = buffer_1000_ppdd[931];

    auto g_z_0_0_0_z_y_zz_xz = buffer_1000_ppdd[932];

    auto g_z_0_0_0_z_y_zz_yy = buffer_1000_ppdd[933];

    auto g_z_0_0_0_z_y_zz_yz = buffer_1000_ppdd[934];

    auto g_z_0_0_0_z_y_zz_zz = buffer_1000_ppdd[935];

    auto g_z_0_0_0_z_z_xx_xx = buffer_1000_ppdd[936];

    auto g_z_0_0_0_z_z_xx_xy = buffer_1000_ppdd[937];

    auto g_z_0_0_0_z_z_xx_xz = buffer_1000_ppdd[938];

    auto g_z_0_0_0_z_z_xx_yy = buffer_1000_ppdd[939];

    auto g_z_0_0_0_z_z_xx_yz = buffer_1000_ppdd[940];

    auto g_z_0_0_0_z_z_xx_zz = buffer_1000_ppdd[941];

    auto g_z_0_0_0_z_z_xy_xx = buffer_1000_ppdd[942];

    auto g_z_0_0_0_z_z_xy_xy = buffer_1000_ppdd[943];

    auto g_z_0_0_0_z_z_xy_xz = buffer_1000_ppdd[944];

    auto g_z_0_0_0_z_z_xy_yy = buffer_1000_ppdd[945];

    auto g_z_0_0_0_z_z_xy_yz = buffer_1000_ppdd[946];

    auto g_z_0_0_0_z_z_xy_zz = buffer_1000_ppdd[947];

    auto g_z_0_0_0_z_z_xz_xx = buffer_1000_ppdd[948];

    auto g_z_0_0_0_z_z_xz_xy = buffer_1000_ppdd[949];

    auto g_z_0_0_0_z_z_xz_xz = buffer_1000_ppdd[950];

    auto g_z_0_0_0_z_z_xz_yy = buffer_1000_ppdd[951];

    auto g_z_0_0_0_z_z_xz_yz = buffer_1000_ppdd[952];

    auto g_z_0_0_0_z_z_xz_zz = buffer_1000_ppdd[953];

    auto g_z_0_0_0_z_z_yy_xx = buffer_1000_ppdd[954];

    auto g_z_0_0_0_z_z_yy_xy = buffer_1000_ppdd[955];

    auto g_z_0_0_0_z_z_yy_xz = buffer_1000_ppdd[956];

    auto g_z_0_0_0_z_z_yy_yy = buffer_1000_ppdd[957];

    auto g_z_0_0_0_z_z_yy_yz = buffer_1000_ppdd[958];

    auto g_z_0_0_0_z_z_yy_zz = buffer_1000_ppdd[959];

    auto g_z_0_0_0_z_z_yz_xx = buffer_1000_ppdd[960];

    auto g_z_0_0_0_z_z_yz_xy = buffer_1000_ppdd[961];

    auto g_z_0_0_0_z_z_yz_xz = buffer_1000_ppdd[962];

    auto g_z_0_0_0_z_z_yz_yy = buffer_1000_ppdd[963];

    auto g_z_0_0_0_z_z_yz_yz = buffer_1000_ppdd[964];

    auto g_z_0_0_0_z_z_yz_zz = buffer_1000_ppdd[965];

    auto g_z_0_0_0_z_z_zz_xx = buffer_1000_ppdd[966];

    auto g_z_0_0_0_z_z_zz_xy = buffer_1000_ppdd[967];

    auto g_z_0_0_0_z_z_zz_xz = buffer_1000_ppdd[968];

    auto g_z_0_0_0_z_z_zz_yy = buffer_1000_ppdd[969];

    auto g_z_0_0_0_z_z_zz_yz = buffer_1000_ppdd[970];

    auto g_z_0_0_0_z_z_zz_zz = buffer_1000_ppdd[971];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_x_xx_xx, g_0_x_xx_xy, g_0_x_xx_xz, g_0_x_xx_yy, g_0_x_xx_yz, g_0_x_xx_zz, g_x_0_0_0_x_x_xx_xx, g_x_0_0_0_x_x_xx_xy, g_x_0_0_0_x_x_xx_xz, g_x_0_0_0_x_x_xx_yy, g_x_0_0_0_x_x_xx_yz, g_x_0_0_0_x_x_xx_zz, g_xx_x_xx_xx, g_xx_x_xx_xy, g_xx_x_xx_xz, g_xx_x_xx_yy, g_xx_x_xx_yz, g_xx_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_x_xx_xx[i] = -g_0_x_xx_xx[i] + 2.0 * g_xx_x_xx_xx[i] * a_exp;

        g_x_0_0_0_x_x_xx_xy[i] = -g_0_x_xx_xy[i] + 2.0 * g_xx_x_xx_xy[i] * a_exp;

        g_x_0_0_0_x_x_xx_xz[i] = -g_0_x_xx_xz[i] + 2.0 * g_xx_x_xx_xz[i] * a_exp;

        g_x_0_0_0_x_x_xx_yy[i] = -g_0_x_xx_yy[i] + 2.0 * g_xx_x_xx_yy[i] * a_exp;

        g_x_0_0_0_x_x_xx_yz[i] = -g_0_x_xx_yz[i] + 2.0 * g_xx_x_xx_yz[i] * a_exp;

        g_x_0_0_0_x_x_xx_zz[i] = -g_0_x_xx_zz[i] + 2.0 * g_xx_x_xx_zz[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_x_xy_xx, g_0_x_xy_xy, g_0_x_xy_xz, g_0_x_xy_yy, g_0_x_xy_yz, g_0_x_xy_zz, g_x_0_0_0_x_x_xy_xx, g_x_0_0_0_x_x_xy_xy, g_x_0_0_0_x_x_xy_xz, g_x_0_0_0_x_x_xy_yy, g_x_0_0_0_x_x_xy_yz, g_x_0_0_0_x_x_xy_zz, g_xx_x_xy_xx, g_xx_x_xy_xy, g_xx_x_xy_xz, g_xx_x_xy_yy, g_xx_x_xy_yz, g_xx_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_x_xy_xx[i] = -g_0_x_xy_xx[i] + 2.0 * g_xx_x_xy_xx[i] * a_exp;

        g_x_0_0_0_x_x_xy_xy[i] = -g_0_x_xy_xy[i] + 2.0 * g_xx_x_xy_xy[i] * a_exp;

        g_x_0_0_0_x_x_xy_xz[i] = -g_0_x_xy_xz[i] + 2.0 * g_xx_x_xy_xz[i] * a_exp;

        g_x_0_0_0_x_x_xy_yy[i] = -g_0_x_xy_yy[i] + 2.0 * g_xx_x_xy_yy[i] * a_exp;

        g_x_0_0_0_x_x_xy_yz[i] = -g_0_x_xy_yz[i] + 2.0 * g_xx_x_xy_yz[i] * a_exp;

        g_x_0_0_0_x_x_xy_zz[i] = -g_0_x_xy_zz[i] + 2.0 * g_xx_x_xy_zz[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_x_xz_xx, g_0_x_xz_xy, g_0_x_xz_xz, g_0_x_xz_yy, g_0_x_xz_yz, g_0_x_xz_zz, g_x_0_0_0_x_x_xz_xx, g_x_0_0_0_x_x_xz_xy, g_x_0_0_0_x_x_xz_xz, g_x_0_0_0_x_x_xz_yy, g_x_0_0_0_x_x_xz_yz, g_x_0_0_0_x_x_xz_zz, g_xx_x_xz_xx, g_xx_x_xz_xy, g_xx_x_xz_xz, g_xx_x_xz_yy, g_xx_x_xz_yz, g_xx_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_x_xz_xx[i] = -g_0_x_xz_xx[i] + 2.0 * g_xx_x_xz_xx[i] * a_exp;

        g_x_0_0_0_x_x_xz_xy[i] = -g_0_x_xz_xy[i] + 2.0 * g_xx_x_xz_xy[i] * a_exp;

        g_x_0_0_0_x_x_xz_xz[i] = -g_0_x_xz_xz[i] + 2.0 * g_xx_x_xz_xz[i] * a_exp;

        g_x_0_0_0_x_x_xz_yy[i] = -g_0_x_xz_yy[i] + 2.0 * g_xx_x_xz_yy[i] * a_exp;

        g_x_0_0_0_x_x_xz_yz[i] = -g_0_x_xz_yz[i] + 2.0 * g_xx_x_xz_yz[i] * a_exp;

        g_x_0_0_0_x_x_xz_zz[i] = -g_0_x_xz_zz[i] + 2.0 * g_xx_x_xz_zz[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_x_yy_xx, g_0_x_yy_xy, g_0_x_yy_xz, g_0_x_yy_yy, g_0_x_yy_yz, g_0_x_yy_zz, g_x_0_0_0_x_x_yy_xx, g_x_0_0_0_x_x_yy_xy, g_x_0_0_0_x_x_yy_xz, g_x_0_0_0_x_x_yy_yy, g_x_0_0_0_x_x_yy_yz, g_x_0_0_0_x_x_yy_zz, g_xx_x_yy_xx, g_xx_x_yy_xy, g_xx_x_yy_xz, g_xx_x_yy_yy, g_xx_x_yy_yz, g_xx_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_x_yy_xx[i] = -g_0_x_yy_xx[i] + 2.0 * g_xx_x_yy_xx[i] * a_exp;

        g_x_0_0_0_x_x_yy_xy[i] = -g_0_x_yy_xy[i] + 2.0 * g_xx_x_yy_xy[i] * a_exp;

        g_x_0_0_0_x_x_yy_xz[i] = -g_0_x_yy_xz[i] + 2.0 * g_xx_x_yy_xz[i] * a_exp;

        g_x_0_0_0_x_x_yy_yy[i] = -g_0_x_yy_yy[i] + 2.0 * g_xx_x_yy_yy[i] * a_exp;

        g_x_0_0_0_x_x_yy_yz[i] = -g_0_x_yy_yz[i] + 2.0 * g_xx_x_yy_yz[i] * a_exp;

        g_x_0_0_0_x_x_yy_zz[i] = -g_0_x_yy_zz[i] + 2.0 * g_xx_x_yy_zz[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_x_yz_xx, g_0_x_yz_xy, g_0_x_yz_xz, g_0_x_yz_yy, g_0_x_yz_yz, g_0_x_yz_zz, g_x_0_0_0_x_x_yz_xx, g_x_0_0_0_x_x_yz_xy, g_x_0_0_0_x_x_yz_xz, g_x_0_0_0_x_x_yz_yy, g_x_0_0_0_x_x_yz_yz, g_x_0_0_0_x_x_yz_zz, g_xx_x_yz_xx, g_xx_x_yz_xy, g_xx_x_yz_xz, g_xx_x_yz_yy, g_xx_x_yz_yz, g_xx_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_x_yz_xx[i] = -g_0_x_yz_xx[i] + 2.0 * g_xx_x_yz_xx[i] * a_exp;

        g_x_0_0_0_x_x_yz_xy[i] = -g_0_x_yz_xy[i] + 2.0 * g_xx_x_yz_xy[i] * a_exp;

        g_x_0_0_0_x_x_yz_xz[i] = -g_0_x_yz_xz[i] + 2.0 * g_xx_x_yz_xz[i] * a_exp;

        g_x_0_0_0_x_x_yz_yy[i] = -g_0_x_yz_yy[i] + 2.0 * g_xx_x_yz_yy[i] * a_exp;

        g_x_0_0_0_x_x_yz_yz[i] = -g_0_x_yz_yz[i] + 2.0 * g_xx_x_yz_yz[i] * a_exp;

        g_x_0_0_0_x_x_yz_zz[i] = -g_0_x_yz_zz[i] + 2.0 * g_xx_x_yz_zz[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_x_zz_xx, g_0_x_zz_xy, g_0_x_zz_xz, g_0_x_zz_yy, g_0_x_zz_yz, g_0_x_zz_zz, g_x_0_0_0_x_x_zz_xx, g_x_0_0_0_x_x_zz_xy, g_x_0_0_0_x_x_zz_xz, g_x_0_0_0_x_x_zz_yy, g_x_0_0_0_x_x_zz_yz, g_x_0_0_0_x_x_zz_zz, g_xx_x_zz_xx, g_xx_x_zz_xy, g_xx_x_zz_xz, g_xx_x_zz_yy, g_xx_x_zz_yz, g_xx_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_x_zz_xx[i] = -g_0_x_zz_xx[i] + 2.0 * g_xx_x_zz_xx[i] * a_exp;

        g_x_0_0_0_x_x_zz_xy[i] = -g_0_x_zz_xy[i] + 2.0 * g_xx_x_zz_xy[i] * a_exp;

        g_x_0_0_0_x_x_zz_xz[i] = -g_0_x_zz_xz[i] + 2.0 * g_xx_x_zz_xz[i] * a_exp;

        g_x_0_0_0_x_x_zz_yy[i] = -g_0_x_zz_yy[i] + 2.0 * g_xx_x_zz_yy[i] * a_exp;

        g_x_0_0_0_x_x_zz_yz[i] = -g_0_x_zz_yz[i] + 2.0 * g_xx_x_zz_yz[i] * a_exp;

        g_x_0_0_0_x_x_zz_zz[i] = -g_0_x_zz_zz[i] + 2.0 * g_xx_x_zz_zz[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_0_y_xx_xx, g_0_y_xx_xy, g_0_y_xx_xz, g_0_y_xx_yy, g_0_y_xx_yz, g_0_y_xx_zz, g_x_0_0_0_x_y_xx_xx, g_x_0_0_0_x_y_xx_xy, g_x_0_0_0_x_y_xx_xz, g_x_0_0_0_x_y_xx_yy, g_x_0_0_0_x_y_xx_yz, g_x_0_0_0_x_y_xx_zz, g_xx_y_xx_xx, g_xx_y_xx_xy, g_xx_y_xx_xz, g_xx_y_xx_yy, g_xx_y_xx_yz, g_xx_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_y_xx_xx[i] = -g_0_y_xx_xx[i] + 2.0 * g_xx_y_xx_xx[i] * a_exp;

        g_x_0_0_0_x_y_xx_xy[i] = -g_0_y_xx_xy[i] + 2.0 * g_xx_y_xx_xy[i] * a_exp;

        g_x_0_0_0_x_y_xx_xz[i] = -g_0_y_xx_xz[i] + 2.0 * g_xx_y_xx_xz[i] * a_exp;

        g_x_0_0_0_x_y_xx_yy[i] = -g_0_y_xx_yy[i] + 2.0 * g_xx_y_xx_yy[i] * a_exp;

        g_x_0_0_0_x_y_xx_yz[i] = -g_0_y_xx_yz[i] + 2.0 * g_xx_y_xx_yz[i] * a_exp;

        g_x_0_0_0_x_y_xx_zz[i] = -g_0_y_xx_zz[i] + 2.0 * g_xx_y_xx_zz[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_0_y_xy_xx, g_0_y_xy_xy, g_0_y_xy_xz, g_0_y_xy_yy, g_0_y_xy_yz, g_0_y_xy_zz, g_x_0_0_0_x_y_xy_xx, g_x_0_0_0_x_y_xy_xy, g_x_0_0_0_x_y_xy_xz, g_x_0_0_0_x_y_xy_yy, g_x_0_0_0_x_y_xy_yz, g_x_0_0_0_x_y_xy_zz, g_xx_y_xy_xx, g_xx_y_xy_xy, g_xx_y_xy_xz, g_xx_y_xy_yy, g_xx_y_xy_yz, g_xx_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_y_xy_xx[i] = -g_0_y_xy_xx[i] + 2.0 * g_xx_y_xy_xx[i] * a_exp;

        g_x_0_0_0_x_y_xy_xy[i] = -g_0_y_xy_xy[i] + 2.0 * g_xx_y_xy_xy[i] * a_exp;

        g_x_0_0_0_x_y_xy_xz[i] = -g_0_y_xy_xz[i] + 2.0 * g_xx_y_xy_xz[i] * a_exp;

        g_x_0_0_0_x_y_xy_yy[i] = -g_0_y_xy_yy[i] + 2.0 * g_xx_y_xy_yy[i] * a_exp;

        g_x_0_0_0_x_y_xy_yz[i] = -g_0_y_xy_yz[i] + 2.0 * g_xx_y_xy_yz[i] * a_exp;

        g_x_0_0_0_x_y_xy_zz[i] = -g_0_y_xy_zz[i] + 2.0 * g_xx_y_xy_zz[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_0_y_xz_xx, g_0_y_xz_xy, g_0_y_xz_xz, g_0_y_xz_yy, g_0_y_xz_yz, g_0_y_xz_zz, g_x_0_0_0_x_y_xz_xx, g_x_0_0_0_x_y_xz_xy, g_x_0_0_0_x_y_xz_xz, g_x_0_0_0_x_y_xz_yy, g_x_0_0_0_x_y_xz_yz, g_x_0_0_0_x_y_xz_zz, g_xx_y_xz_xx, g_xx_y_xz_xy, g_xx_y_xz_xz, g_xx_y_xz_yy, g_xx_y_xz_yz, g_xx_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_y_xz_xx[i] = -g_0_y_xz_xx[i] + 2.0 * g_xx_y_xz_xx[i] * a_exp;

        g_x_0_0_0_x_y_xz_xy[i] = -g_0_y_xz_xy[i] + 2.0 * g_xx_y_xz_xy[i] * a_exp;

        g_x_0_0_0_x_y_xz_xz[i] = -g_0_y_xz_xz[i] + 2.0 * g_xx_y_xz_xz[i] * a_exp;

        g_x_0_0_0_x_y_xz_yy[i] = -g_0_y_xz_yy[i] + 2.0 * g_xx_y_xz_yy[i] * a_exp;

        g_x_0_0_0_x_y_xz_yz[i] = -g_0_y_xz_yz[i] + 2.0 * g_xx_y_xz_yz[i] * a_exp;

        g_x_0_0_0_x_y_xz_zz[i] = -g_0_y_xz_zz[i] + 2.0 * g_xx_y_xz_zz[i] * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_0_y_yy_xx, g_0_y_yy_xy, g_0_y_yy_xz, g_0_y_yy_yy, g_0_y_yy_yz, g_0_y_yy_zz, g_x_0_0_0_x_y_yy_xx, g_x_0_0_0_x_y_yy_xy, g_x_0_0_0_x_y_yy_xz, g_x_0_0_0_x_y_yy_yy, g_x_0_0_0_x_y_yy_yz, g_x_0_0_0_x_y_yy_zz, g_xx_y_yy_xx, g_xx_y_yy_xy, g_xx_y_yy_xz, g_xx_y_yy_yy, g_xx_y_yy_yz, g_xx_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_y_yy_xx[i] = -g_0_y_yy_xx[i] + 2.0 * g_xx_y_yy_xx[i] * a_exp;

        g_x_0_0_0_x_y_yy_xy[i] = -g_0_y_yy_xy[i] + 2.0 * g_xx_y_yy_xy[i] * a_exp;

        g_x_0_0_0_x_y_yy_xz[i] = -g_0_y_yy_xz[i] + 2.0 * g_xx_y_yy_xz[i] * a_exp;

        g_x_0_0_0_x_y_yy_yy[i] = -g_0_y_yy_yy[i] + 2.0 * g_xx_y_yy_yy[i] * a_exp;

        g_x_0_0_0_x_y_yy_yz[i] = -g_0_y_yy_yz[i] + 2.0 * g_xx_y_yy_yz[i] * a_exp;

        g_x_0_0_0_x_y_yy_zz[i] = -g_0_y_yy_zz[i] + 2.0 * g_xx_y_yy_zz[i] * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_0_y_yz_xx, g_0_y_yz_xy, g_0_y_yz_xz, g_0_y_yz_yy, g_0_y_yz_yz, g_0_y_yz_zz, g_x_0_0_0_x_y_yz_xx, g_x_0_0_0_x_y_yz_xy, g_x_0_0_0_x_y_yz_xz, g_x_0_0_0_x_y_yz_yy, g_x_0_0_0_x_y_yz_yz, g_x_0_0_0_x_y_yz_zz, g_xx_y_yz_xx, g_xx_y_yz_xy, g_xx_y_yz_xz, g_xx_y_yz_yy, g_xx_y_yz_yz, g_xx_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_y_yz_xx[i] = -g_0_y_yz_xx[i] + 2.0 * g_xx_y_yz_xx[i] * a_exp;

        g_x_0_0_0_x_y_yz_xy[i] = -g_0_y_yz_xy[i] + 2.0 * g_xx_y_yz_xy[i] * a_exp;

        g_x_0_0_0_x_y_yz_xz[i] = -g_0_y_yz_xz[i] + 2.0 * g_xx_y_yz_xz[i] * a_exp;

        g_x_0_0_0_x_y_yz_yy[i] = -g_0_y_yz_yy[i] + 2.0 * g_xx_y_yz_yy[i] * a_exp;

        g_x_0_0_0_x_y_yz_yz[i] = -g_0_y_yz_yz[i] + 2.0 * g_xx_y_yz_yz[i] * a_exp;

        g_x_0_0_0_x_y_yz_zz[i] = -g_0_y_yz_zz[i] + 2.0 * g_xx_y_yz_zz[i] * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_0_y_zz_xx, g_0_y_zz_xy, g_0_y_zz_xz, g_0_y_zz_yy, g_0_y_zz_yz, g_0_y_zz_zz, g_x_0_0_0_x_y_zz_xx, g_x_0_0_0_x_y_zz_xy, g_x_0_0_0_x_y_zz_xz, g_x_0_0_0_x_y_zz_yy, g_x_0_0_0_x_y_zz_yz, g_x_0_0_0_x_y_zz_zz, g_xx_y_zz_xx, g_xx_y_zz_xy, g_xx_y_zz_xz, g_xx_y_zz_yy, g_xx_y_zz_yz, g_xx_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_y_zz_xx[i] = -g_0_y_zz_xx[i] + 2.0 * g_xx_y_zz_xx[i] * a_exp;

        g_x_0_0_0_x_y_zz_xy[i] = -g_0_y_zz_xy[i] + 2.0 * g_xx_y_zz_xy[i] * a_exp;

        g_x_0_0_0_x_y_zz_xz[i] = -g_0_y_zz_xz[i] + 2.0 * g_xx_y_zz_xz[i] * a_exp;

        g_x_0_0_0_x_y_zz_yy[i] = -g_0_y_zz_yy[i] + 2.0 * g_xx_y_zz_yy[i] * a_exp;

        g_x_0_0_0_x_y_zz_yz[i] = -g_0_y_zz_yz[i] + 2.0 * g_xx_y_zz_yz[i] * a_exp;

        g_x_0_0_0_x_y_zz_zz[i] = -g_0_y_zz_zz[i] + 2.0 * g_xx_y_zz_zz[i] * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_0_z_xx_xx, g_0_z_xx_xy, g_0_z_xx_xz, g_0_z_xx_yy, g_0_z_xx_yz, g_0_z_xx_zz, g_x_0_0_0_x_z_xx_xx, g_x_0_0_0_x_z_xx_xy, g_x_0_0_0_x_z_xx_xz, g_x_0_0_0_x_z_xx_yy, g_x_0_0_0_x_z_xx_yz, g_x_0_0_0_x_z_xx_zz, g_xx_z_xx_xx, g_xx_z_xx_xy, g_xx_z_xx_xz, g_xx_z_xx_yy, g_xx_z_xx_yz, g_xx_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_z_xx_xx[i] = -g_0_z_xx_xx[i] + 2.0 * g_xx_z_xx_xx[i] * a_exp;

        g_x_0_0_0_x_z_xx_xy[i] = -g_0_z_xx_xy[i] + 2.0 * g_xx_z_xx_xy[i] * a_exp;

        g_x_0_0_0_x_z_xx_xz[i] = -g_0_z_xx_xz[i] + 2.0 * g_xx_z_xx_xz[i] * a_exp;

        g_x_0_0_0_x_z_xx_yy[i] = -g_0_z_xx_yy[i] + 2.0 * g_xx_z_xx_yy[i] * a_exp;

        g_x_0_0_0_x_z_xx_yz[i] = -g_0_z_xx_yz[i] + 2.0 * g_xx_z_xx_yz[i] * a_exp;

        g_x_0_0_0_x_z_xx_zz[i] = -g_0_z_xx_zz[i] + 2.0 * g_xx_z_xx_zz[i] * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_0_z_xy_xx, g_0_z_xy_xy, g_0_z_xy_xz, g_0_z_xy_yy, g_0_z_xy_yz, g_0_z_xy_zz, g_x_0_0_0_x_z_xy_xx, g_x_0_0_0_x_z_xy_xy, g_x_0_0_0_x_z_xy_xz, g_x_0_0_0_x_z_xy_yy, g_x_0_0_0_x_z_xy_yz, g_x_0_0_0_x_z_xy_zz, g_xx_z_xy_xx, g_xx_z_xy_xy, g_xx_z_xy_xz, g_xx_z_xy_yy, g_xx_z_xy_yz, g_xx_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_z_xy_xx[i] = -g_0_z_xy_xx[i] + 2.0 * g_xx_z_xy_xx[i] * a_exp;

        g_x_0_0_0_x_z_xy_xy[i] = -g_0_z_xy_xy[i] + 2.0 * g_xx_z_xy_xy[i] * a_exp;

        g_x_0_0_0_x_z_xy_xz[i] = -g_0_z_xy_xz[i] + 2.0 * g_xx_z_xy_xz[i] * a_exp;

        g_x_0_0_0_x_z_xy_yy[i] = -g_0_z_xy_yy[i] + 2.0 * g_xx_z_xy_yy[i] * a_exp;

        g_x_0_0_0_x_z_xy_yz[i] = -g_0_z_xy_yz[i] + 2.0 * g_xx_z_xy_yz[i] * a_exp;

        g_x_0_0_0_x_z_xy_zz[i] = -g_0_z_xy_zz[i] + 2.0 * g_xx_z_xy_zz[i] * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_0_z_xz_xx, g_0_z_xz_xy, g_0_z_xz_xz, g_0_z_xz_yy, g_0_z_xz_yz, g_0_z_xz_zz, g_x_0_0_0_x_z_xz_xx, g_x_0_0_0_x_z_xz_xy, g_x_0_0_0_x_z_xz_xz, g_x_0_0_0_x_z_xz_yy, g_x_0_0_0_x_z_xz_yz, g_x_0_0_0_x_z_xz_zz, g_xx_z_xz_xx, g_xx_z_xz_xy, g_xx_z_xz_xz, g_xx_z_xz_yy, g_xx_z_xz_yz, g_xx_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_z_xz_xx[i] = -g_0_z_xz_xx[i] + 2.0 * g_xx_z_xz_xx[i] * a_exp;

        g_x_0_0_0_x_z_xz_xy[i] = -g_0_z_xz_xy[i] + 2.0 * g_xx_z_xz_xy[i] * a_exp;

        g_x_0_0_0_x_z_xz_xz[i] = -g_0_z_xz_xz[i] + 2.0 * g_xx_z_xz_xz[i] * a_exp;

        g_x_0_0_0_x_z_xz_yy[i] = -g_0_z_xz_yy[i] + 2.0 * g_xx_z_xz_yy[i] * a_exp;

        g_x_0_0_0_x_z_xz_yz[i] = -g_0_z_xz_yz[i] + 2.0 * g_xx_z_xz_yz[i] * a_exp;

        g_x_0_0_0_x_z_xz_zz[i] = -g_0_z_xz_zz[i] + 2.0 * g_xx_z_xz_zz[i] * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_0_z_yy_xx, g_0_z_yy_xy, g_0_z_yy_xz, g_0_z_yy_yy, g_0_z_yy_yz, g_0_z_yy_zz, g_x_0_0_0_x_z_yy_xx, g_x_0_0_0_x_z_yy_xy, g_x_0_0_0_x_z_yy_xz, g_x_0_0_0_x_z_yy_yy, g_x_0_0_0_x_z_yy_yz, g_x_0_0_0_x_z_yy_zz, g_xx_z_yy_xx, g_xx_z_yy_xy, g_xx_z_yy_xz, g_xx_z_yy_yy, g_xx_z_yy_yz, g_xx_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_z_yy_xx[i] = -g_0_z_yy_xx[i] + 2.0 * g_xx_z_yy_xx[i] * a_exp;

        g_x_0_0_0_x_z_yy_xy[i] = -g_0_z_yy_xy[i] + 2.0 * g_xx_z_yy_xy[i] * a_exp;

        g_x_0_0_0_x_z_yy_xz[i] = -g_0_z_yy_xz[i] + 2.0 * g_xx_z_yy_xz[i] * a_exp;

        g_x_0_0_0_x_z_yy_yy[i] = -g_0_z_yy_yy[i] + 2.0 * g_xx_z_yy_yy[i] * a_exp;

        g_x_0_0_0_x_z_yy_yz[i] = -g_0_z_yy_yz[i] + 2.0 * g_xx_z_yy_yz[i] * a_exp;

        g_x_0_0_0_x_z_yy_zz[i] = -g_0_z_yy_zz[i] + 2.0 * g_xx_z_yy_zz[i] * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_0_z_yz_xx, g_0_z_yz_xy, g_0_z_yz_xz, g_0_z_yz_yy, g_0_z_yz_yz, g_0_z_yz_zz, g_x_0_0_0_x_z_yz_xx, g_x_0_0_0_x_z_yz_xy, g_x_0_0_0_x_z_yz_xz, g_x_0_0_0_x_z_yz_yy, g_x_0_0_0_x_z_yz_yz, g_x_0_0_0_x_z_yz_zz, g_xx_z_yz_xx, g_xx_z_yz_xy, g_xx_z_yz_xz, g_xx_z_yz_yy, g_xx_z_yz_yz, g_xx_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_z_yz_xx[i] = -g_0_z_yz_xx[i] + 2.0 * g_xx_z_yz_xx[i] * a_exp;

        g_x_0_0_0_x_z_yz_xy[i] = -g_0_z_yz_xy[i] + 2.0 * g_xx_z_yz_xy[i] * a_exp;

        g_x_0_0_0_x_z_yz_xz[i] = -g_0_z_yz_xz[i] + 2.0 * g_xx_z_yz_xz[i] * a_exp;

        g_x_0_0_0_x_z_yz_yy[i] = -g_0_z_yz_yy[i] + 2.0 * g_xx_z_yz_yy[i] * a_exp;

        g_x_0_0_0_x_z_yz_yz[i] = -g_0_z_yz_yz[i] + 2.0 * g_xx_z_yz_yz[i] * a_exp;

        g_x_0_0_0_x_z_yz_zz[i] = -g_0_z_yz_zz[i] + 2.0 * g_xx_z_yz_zz[i] * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_0_z_zz_xx, g_0_z_zz_xy, g_0_z_zz_xz, g_0_z_zz_yy, g_0_z_zz_yz, g_0_z_zz_zz, g_x_0_0_0_x_z_zz_xx, g_x_0_0_0_x_z_zz_xy, g_x_0_0_0_x_z_zz_xz, g_x_0_0_0_x_z_zz_yy, g_x_0_0_0_x_z_zz_yz, g_x_0_0_0_x_z_zz_zz, g_xx_z_zz_xx, g_xx_z_zz_xy, g_xx_z_zz_xz, g_xx_z_zz_yy, g_xx_z_zz_yz, g_xx_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_z_zz_xx[i] = -g_0_z_zz_xx[i] + 2.0 * g_xx_z_zz_xx[i] * a_exp;

        g_x_0_0_0_x_z_zz_xy[i] = -g_0_z_zz_xy[i] + 2.0 * g_xx_z_zz_xy[i] * a_exp;

        g_x_0_0_0_x_z_zz_xz[i] = -g_0_z_zz_xz[i] + 2.0 * g_xx_z_zz_xz[i] * a_exp;

        g_x_0_0_0_x_z_zz_yy[i] = -g_0_z_zz_yy[i] + 2.0 * g_xx_z_zz_yy[i] * a_exp;

        g_x_0_0_0_x_z_zz_yz[i] = -g_0_z_zz_yz[i] + 2.0 * g_xx_z_zz_yz[i] * a_exp;

        g_x_0_0_0_x_z_zz_zz[i] = -g_0_z_zz_zz[i] + 2.0 * g_xx_z_zz_zz[i] * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_0_0_0_y_x_xx_xx, g_x_0_0_0_y_x_xx_xy, g_x_0_0_0_y_x_xx_xz, g_x_0_0_0_y_x_xx_yy, g_x_0_0_0_y_x_xx_yz, g_x_0_0_0_y_x_xx_zz, g_xy_x_xx_xx, g_xy_x_xx_xy, g_xy_x_xx_xz, g_xy_x_xx_yy, g_xy_x_xx_yz, g_xy_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_x_xx_xx[i] = 2.0 * g_xy_x_xx_xx[i] * a_exp;

        g_x_0_0_0_y_x_xx_xy[i] = 2.0 * g_xy_x_xx_xy[i] * a_exp;

        g_x_0_0_0_y_x_xx_xz[i] = 2.0 * g_xy_x_xx_xz[i] * a_exp;

        g_x_0_0_0_y_x_xx_yy[i] = 2.0 * g_xy_x_xx_yy[i] * a_exp;

        g_x_0_0_0_y_x_xx_yz[i] = 2.0 * g_xy_x_xx_yz[i] * a_exp;

        g_x_0_0_0_y_x_xx_zz[i] = 2.0 * g_xy_x_xx_zz[i] * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_0_0_0_y_x_xy_xx, g_x_0_0_0_y_x_xy_xy, g_x_0_0_0_y_x_xy_xz, g_x_0_0_0_y_x_xy_yy, g_x_0_0_0_y_x_xy_yz, g_x_0_0_0_y_x_xy_zz, g_xy_x_xy_xx, g_xy_x_xy_xy, g_xy_x_xy_xz, g_xy_x_xy_yy, g_xy_x_xy_yz, g_xy_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_x_xy_xx[i] = 2.0 * g_xy_x_xy_xx[i] * a_exp;

        g_x_0_0_0_y_x_xy_xy[i] = 2.0 * g_xy_x_xy_xy[i] * a_exp;

        g_x_0_0_0_y_x_xy_xz[i] = 2.0 * g_xy_x_xy_xz[i] * a_exp;

        g_x_0_0_0_y_x_xy_yy[i] = 2.0 * g_xy_x_xy_yy[i] * a_exp;

        g_x_0_0_0_y_x_xy_yz[i] = 2.0 * g_xy_x_xy_yz[i] * a_exp;

        g_x_0_0_0_y_x_xy_zz[i] = 2.0 * g_xy_x_xy_zz[i] * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_0_0_0_y_x_xz_xx, g_x_0_0_0_y_x_xz_xy, g_x_0_0_0_y_x_xz_xz, g_x_0_0_0_y_x_xz_yy, g_x_0_0_0_y_x_xz_yz, g_x_0_0_0_y_x_xz_zz, g_xy_x_xz_xx, g_xy_x_xz_xy, g_xy_x_xz_xz, g_xy_x_xz_yy, g_xy_x_xz_yz, g_xy_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_x_xz_xx[i] = 2.0 * g_xy_x_xz_xx[i] * a_exp;

        g_x_0_0_0_y_x_xz_xy[i] = 2.0 * g_xy_x_xz_xy[i] * a_exp;

        g_x_0_0_0_y_x_xz_xz[i] = 2.0 * g_xy_x_xz_xz[i] * a_exp;

        g_x_0_0_0_y_x_xz_yy[i] = 2.0 * g_xy_x_xz_yy[i] * a_exp;

        g_x_0_0_0_y_x_xz_yz[i] = 2.0 * g_xy_x_xz_yz[i] * a_exp;

        g_x_0_0_0_y_x_xz_zz[i] = 2.0 * g_xy_x_xz_zz[i] * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_0_0_0_y_x_yy_xx, g_x_0_0_0_y_x_yy_xy, g_x_0_0_0_y_x_yy_xz, g_x_0_0_0_y_x_yy_yy, g_x_0_0_0_y_x_yy_yz, g_x_0_0_0_y_x_yy_zz, g_xy_x_yy_xx, g_xy_x_yy_xy, g_xy_x_yy_xz, g_xy_x_yy_yy, g_xy_x_yy_yz, g_xy_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_x_yy_xx[i] = 2.0 * g_xy_x_yy_xx[i] * a_exp;

        g_x_0_0_0_y_x_yy_xy[i] = 2.0 * g_xy_x_yy_xy[i] * a_exp;

        g_x_0_0_0_y_x_yy_xz[i] = 2.0 * g_xy_x_yy_xz[i] * a_exp;

        g_x_0_0_0_y_x_yy_yy[i] = 2.0 * g_xy_x_yy_yy[i] * a_exp;

        g_x_0_0_0_y_x_yy_yz[i] = 2.0 * g_xy_x_yy_yz[i] * a_exp;

        g_x_0_0_0_y_x_yy_zz[i] = 2.0 * g_xy_x_yy_zz[i] * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_0_0_0_y_x_yz_xx, g_x_0_0_0_y_x_yz_xy, g_x_0_0_0_y_x_yz_xz, g_x_0_0_0_y_x_yz_yy, g_x_0_0_0_y_x_yz_yz, g_x_0_0_0_y_x_yz_zz, g_xy_x_yz_xx, g_xy_x_yz_xy, g_xy_x_yz_xz, g_xy_x_yz_yy, g_xy_x_yz_yz, g_xy_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_x_yz_xx[i] = 2.0 * g_xy_x_yz_xx[i] * a_exp;

        g_x_0_0_0_y_x_yz_xy[i] = 2.0 * g_xy_x_yz_xy[i] * a_exp;

        g_x_0_0_0_y_x_yz_xz[i] = 2.0 * g_xy_x_yz_xz[i] * a_exp;

        g_x_0_0_0_y_x_yz_yy[i] = 2.0 * g_xy_x_yz_yy[i] * a_exp;

        g_x_0_0_0_y_x_yz_yz[i] = 2.0 * g_xy_x_yz_yz[i] * a_exp;

        g_x_0_0_0_y_x_yz_zz[i] = 2.0 * g_xy_x_yz_zz[i] * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_0_0_0_y_x_zz_xx, g_x_0_0_0_y_x_zz_xy, g_x_0_0_0_y_x_zz_xz, g_x_0_0_0_y_x_zz_yy, g_x_0_0_0_y_x_zz_yz, g_x_0_0_0_y_x_zz_zz, g_xy_x_zz_xx, g_xy_x_zz_xy, g_xy_x_zz_xz, g_xy_x_zz_yy, g_xy_x_zz_yz, g_xy_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_x_zz_xx[i] = 2.0 * g_xy_x_zz_xx[i] * a_exp;

        g_x_0_0_0_y_x_zz_xy[i] = 2.0 * g_xy_x_zz_xy[i] * a_exp;

        g_x_0_0_0_y_x_zz_xz[i] = 2.0 * g_xy_x_zz_xz[i] * a_exp;

        g_x_0_0_0_y_x_zz_yy[i] = 2.0 * g_xy_x_zz_yy[i] * a_exp;

        g_x_0_0_0_y_x_zz_yz[i] = 2.0 * g_xy_x_zz_yz[i] * a_exp;

        g_x_0_0_0_y_x_zz_zz[i] = 2.0 * g_xy_x_zz_zz[i] * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_0_0_0_y_y_xx_xx, g_x_0_0_0_y_y_xx_xy, g_x_0_0_0_y_y_xx_xz, g_x_0_0_0_y_y_xx_yy, g_x_0_0_0_y_y_xx_yz, g_x_0_0_0_y_y_xx_zz, g_xy_y_xx_xx, g_xy_y_xx_xy, g_xy_y_xx_xz, g_xy_y_xx_yy, g_xy_y_xx_yz, g_xy_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_y_xx_xx[i] = 2.0 * g_xy_y_xx_xx[i] * a_exp;

        g_x_0_0_0_y_y_xx_xy[i] = 2.0 * g_xy_y_xx_xy[i] * a_exp;

        g_x_0_0_0_y_y_xx_xz[i] = 2.0 * g_xy_y_xx_xz[i] * a_exp;

        g_x_0_0_0_y_y_xx_yy[i] = 2.0 * g_xy_y_xx_yy[i] * a_exp;

        g_x_0_0_0_y_y_xx_yz[i] = 2.0 * g_xy_y_xx_yz[i] * a_exp;

        g_x_0_0_0_y_y_xx_zz[i] = 2.0 * g_xy_y_xx_zz[i] * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_0_0_0_y_y_xy_xx, g_x_0_0_0_y_y_xy_xy, g_x_0_0_0_y_y_xy_xz, g_x_0_0_0_y_y_xy_yy, g_x_0_0_0_y_y_xy_yz, g_x_0_0_0_y_y_xy_zz, g_xy_y_xy_xx, g_xy_y_xy_xy, g_xy_y_xy_xz, g_xy_y_xy_yy, g_xy_y_xy_yz, g_xy_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_y_xy_xx[i] = 2.0 * g_xy_y_xy_xx[i] * a_exp;

        g_x_0_0_0_y_y_xy_xy[i] = 2.0 * g_xy_y_xy_xy[i] * a_exp;

        g_x_0_0_0_y_y_xy_xz[i] = 2.0 * g_xy_y_xy_xz[i] * a_exp;

        g_x_0_0_0_y_y_xy_yy[i] = 2.0 * g_xy_y_xy_yy[i] * a_exp;

        g_x_0_0_0_y_y_xy_yz[i] = 2.0 * g_xy_y_xy_yz[i] * a_exp;

        g_x_0_0_0_y_y_xy_zz[i] = 2.0 * g_xy_y_xy_zz[i] * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_0_0_0_y_y_xz_xx, g_x_0_0_0_y_y_xz_xy, g_x_0_0_0_y_y_xz_xz, g_x_0_0_0_y_y_xz_yy, g_x_0_0_0_y_y_xz_yz, g_x_0_0_0_y_y_xz_zz, g_xy_y_xz_xx, g_xy_y_xz_xy, g_xy_y_xz_xz, g_xy_y_xz_yy, g_xy_y_xz_yz, g_xy_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_y_xz_xx[i] = 2.0 * g_xy_y_xz_xx[i] * a_exp;

        g_x_0_0_0_y_y_xz_xy[i] = 2.0 * g_xy_y_xz_xy[i] * a_exp;

        g_x_0_0_0_y_y_xz_xz[i] = 2.0 * g_xy_y_xz_xz[i] * a_exp;

        g_x_0_0_0_y_y_xz_yy[i] = 2.0 * g_xy_y_xz_yy[i] * a_exp;

        g_x_0_0_0_y_y_xz_yz[i] = 2.0 * g_xy_y_xz_yz[i] * a_exp;

        g_x_0_0_0_y_y_xz_zz[i] = 2.0 * g_xy_y_xz_zz[i] * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_0_0_0_y_y_yy_xx, g_x_0_0_0_y_y_yy_xy, g_x_0_0_0_y_y_yy_xz, g_x_0_0_0_y_y_yy_yy, g_x_0_0_0_y_y_yy_yz, g_x_0_0_0_y_y_yy_zz, g_xy_y_yy_xx, g_xy_y_yy_xy, g_xy_y_yy_xz, g_xy_y_yy_yy, g_xy_y_yy_yz, g_xy_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_y_yy_xx[i] = 2.0 * g_xy_y_yy_xx[i] * a_exp;

        g_x_0_0_0_y_y_yy_xy[i] = 2.0 * g_xy_y_yy_xy[i] * a_exp;

        g_x_0_0_0_y_y_yy_xz[i] = 2.0 * g_xy_y_yy_xz[i] * a_exp;

        g_x_0_0_0_y_y_yy_yy[i] = 2.0 * g_xy_y_yy_yy[i] * a_exp;

        g_x_0_0_0_y_y_yy_yz[i] = 2.0 * g_xy_y_yy_yz[i] * a_exp;

        g_x_0_0_0_y_y_yy_zz[i] = 2.0 * g_xy_y_yy_zz[i] * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_0_0_0_y_y_yz_xx, g_x_0_0_0_y_y_yz_xy, g_x_0_0_0_y_y_yz_xz, g_x_0_0_0_y_y_yz_yy, g_x_0_0_0_y_y_yz_yz, g_x_0_0_0_y_y_yz_zz, g_xy_y_yz_xx, g_xy_y_yz_xy, g_xy_y_yz_xz, g_xy_y_yz_yy, g_xy_y_yz_yz, g_xy_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_y_yz_xx[i] = 2.0 * g_xy_y_yz_xx[i] * a_exp;

        g_x_0_0_0_y_y_yz_xy[i] = 2.0 * g_xy_y_yz_xy[i] * a_exp;

        g_x_0_0_0_y_y_yz_xz[i] = 2.0 * g_xy_y_yz_xz[i] * a_exp;

        g_x_0_0_0_y_y_yz_yy[i] = 2.0 * g_xy_y_yz_yy[i] * a_exp;

        g_x_0_0_0_y_y_yz_yz[i] = 2.0 * g_xy_y_yz_yz[i] * a_exp;

        g_x_0_0_0_y_y_yz_zz[i] = 2.0 * g_xy_y_yz_zz[i] * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_0_0_0_y_y_zz_xx, g_x_0_0_0_y_y_zz_xy, g_x_0_0_0_y_y_zz_xz, g_x_0_0_0_y_y_zz_yy, g_x_0_0_0_y_y_zz_yz, g_x_0_0_0_y_y_zz_zz, g_xy_y_zz_xx, g_xy_y_zz_xy, g_xy_y_zz_xz, g_xy_y_zz_yy, g_xy_y_zz_yz, g_xy_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_y_zz_xx[i] = 2.0 * g_xy_y_zz_xx[i] * a_exp;

        g_x_0_0_0_y_y_zz_xy[i] = 2.0 * g_xy_y_zz_xy[i] * a_exp;

        g_x_0_0_0_y_y_zz_xz[i] = 2.0 * g_xy_y_zz_xz[i] * a_exp;

        g_x_0_0_0_y_y_zz_yy[i] = 2.0 * g_xy_y_zz_yy[i] * a_exp;

        g_x_0_0_0_y_y_zz_yz[i] = 2.0 * g_xy_y_zz_yz[i] * a_exp;

        g_x_0_0_0_y_y_zz_zz[i] = 2.0 * g_xy_y_zz_zz[i] * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_0_0_0_y_z_xx_xx, g_x_0_0_0_y_z_xx_xy, g_x_0_0_0_y_z_xx_xz, g_x_0_0_0_y_z_xx_yy, g_x_0_0_0_y_z_xx_yz, g_x_0_0_0_y_z_xx_zz, g_xy_z_xx_xx, g_xy_z_xx_xy, g_xy_z_xx_xz, g_xy_z_xx_yy, g_xy_z_xx_yz, g_xy_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_z_xx_xx[i] = 2.0 * g_xy_z_xx_xx[i] * a_exp;

        g_x_0_0_0_y_z_xx_xy[i] = 2.0 * g_xy_z_xx_xy[i] * a_exp;

        g_x_0_0_0_y_z_xx_xz[i] = 2.0 * g_xy_z_xx_xz[i] * a_exp;

        g_x_0_0_0_y_z_xx_yy[i] = 2.0 * g_xy_z_xx_yy[i] * a_exp;

        g_x_0_0_0_y_z_xx_yz[i] = 2.0 * g_xy_z_xx_yz[i] * a_exp;

        g_x_0_0_0_y_z_xx_zz[i] = 2.0 * g_xy_z_xx_zz[i] * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_0_0_0_y_z_xy_xx, g_x_0_0_0_y_z_xy_xy, g_x_0_0_0_y_z_xy_xz, g_x_0_0_0_y_z_xy_yy, g_x_0_0_0_y_z_xy_yz, g_x_0_0_0_y_z_xy_zz, g_xy_z_xy_xx, g_xy_z_xy_xy, g_xy_z_xy_xz, g_xy_z_xy_yy, g_xy_z_xy_yz, g_xy_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_z_xy_xx[i] = 2.0 * g_xy_z_xy_xx[i] * a_exp;

        g_x_0_0_0_y_z_xy_xy[i] = 2.0 * g_xy_z_xy_xy[i] * a_exp;

        g_x_0_0_0_y_z_xy_xz[i] = 2.0 * g_xy_z_xy_xz[i] * a_exp;

        g_x_0_0_0_y_z_xy_yy[i] = 2.0 * g_xy_z_xy_yy[i] * a_exp;

        g_x_0_0_0_y_z_xy_yz[i] = 2.0 * g_xy_z_xy_yz[i] * a_exp;

        g_x_0_0_0_y_z_xy_zz[i] = 2.0 * g_xy_z_xy_zz[i] * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_0_0_0_y_z_xz_xx, g_x_0_0_0_y_z_xz_xy, g_x_0_0_0_y_z_xz_xz, g_x_0_0_0_y_z_xz_yy, g_x_0_0_0_y_z_xz_yz, g_x_0_0_0_y_z_xz_zz, g_xy_z_xz_xx, g_xy_z_xz_xy, g_xy_z_xz_xz, g_xy_z_xz_yy, g_xy_z_xz_yz, g_xy_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_z_xz_xx[i] = 2.0 * g_xy_z_xz_xx[i] * a_exp;

        g_x_0_0_0_y_z_xz_xy[i] = 2.0 * g_xy_z_xz_xy[i] * a_exp;

        g_x_0_0_0_y_z_xz_xz[i] = 2.0 * g_xy_z_xz_xz[i] * a_exp;

        g_x_0_0_0_y_z_xz_yy[i] = 2.0 * g_xy_z_xz_yy[i] * a_exp;

        g_x_0_0_0_y_z_xz_yz[i] = 2.0 * g_xy_z_xz_yz[i] * a_exp;

        g_x_0_0_0_y_z_xz_zz[i] = 2.0 * g_xy_z_xz_zz[i] * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_0_0_0_y_z_yy_xx, g_x_0_0_0_y_z_yy_xy, g_x_0_0_0_y_z_yy_xz, g_x_0_0_0_y_z_yy_yy, g_x_0_0_0_y_z_yy_yz, g_x_0_0_0_y_z_yy_zz, g_xy_z_yy_xx, g_xy_z_yy_xy, g_xy_z_yy_xz, g_xy_z_yy_yy, g_xy_z_yy_yz, g_xy_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_z_yy_xx[i] = 2.0 * g_xy_z_yy_xx[i] * a_exp;

        g_x_0_0_0_y_z_yy_xy[i] = 2.0 * g_xy_z_yy_xy[i] * a_exp;

        g_x_0_0_0_y_z_yy_xz[i] = 2.0 * g_xy_z_yy_xz[i] * a_exp;

        g_x_0_0_0_y_z_yy_yy[i] = 2.0 * g_xy_z_yy_yy[i] * a_exp;

        g_x_0_0_0_y_z_yy_yz[i] = 2.0 * g_xy_z_yy_yz[i] * a_exp;

        g_x_0_0_0_y_z_yy_zz[i] = 2.0 * g_xy_z_yy_zz[i] * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_0_0_0_y_z_yz_xx, g_x_0_0_0_y_z_yz_xy, g_x_0_0_0_y_z_yz_xz, g_x_0_0_0_y_z_yz_yy, g_x_0_0_0_y_z_yz_yz, g_x_0_0_0_y_z_yz_zz, g_xy_z_yz_xx, g_xy_z_yz_xy, g_xy_z_yz_xz, g_xy_z_yz_yy, g_xy_z_yz_yz, g_xy_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_z_yz_xx[i] = 2.0 * g_xy_z_yz_xx[i] * a_exp;

        g_x_0_0_0_y_z_yz_xy[i] = 2.0 * g_xy_z_yz_xy[i] * a_exp;

        g_x_0_0_0_y_z_yz_xz[i] = 2.0 * g_xy_z_yz_xz[i] * a_exp;

        g_x_0_0_0_y_z_yz_yy[i] = 2.0 * g_xy_z_yz_yy[i] * a_exp;

        g_x_0_0_0_y_z_yz_yz[i] = 2.0 * g_xy_z_yz_yz[i] * a_exp;

        g_x_0_0_0_y_z_yz_zz[i] = 2.0 * g_xy_z_yz_zz[i] * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_0_0_0_y_z_zz_xx, g_x_0_0_0_y_z_zz_xy, g_x_0_0_0_y_z_zz_xz, g_x_0_0_0_y_z_zz_yy, g_x_0_0_0_y_z_zz_yz, g_x_0_0_0_y_z_zz_zz, g_xy_z_zz_xx, g_xy_z_zz_xy, g_xy_z_zz_xz, g_xy_z_zz_yy, g_xy_z_zz_yz, g_xy_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_z_zz_xx[i] = 2.0 * g_xy_z_zz_xx[i] * a_exp;

        g_x_0_0_0_y_z_zz_xy[i] = 2.0 * g_xy_z_zz_xy[i] * a_exp;

        g_x_0_0_0_y_z_zz_xz[i] = 2.0 * g_xy_z_zz_xz[i] * a_exp;

        g_x_0_0_0_y_z_zz_yy[i] = 2.0 * g_xy_z_zz_yy[i] * a_exp;

        g_x_0_0_0_y_z_zz_yz[i] = 2.0 * g_xy_z_zz_yz[i] * a_exp;

        g_x_0_0_0_y_z_zz_zz[i] = 2.0 * g_xy_z_zz_zz[i] * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_0_0_0_z_x_xx_xx, g_x_0_0_0_z_x_xx_xy, g_x_0_0_0_z_x_xx_xz, g_x_0_0_0_z_x_xx_yy, g_x_0_0_0_z_x_xx_yz, g_x_0_0_0_z_x_xx_zz, g_xz_x_xx_xx, g_xz_x_xx_xy, g_xz_x_xx_xz, g_xz_x_xx_yy, g_xz_x_xx_yz, g_xz_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_x_xx_xx[i] = 2.0 * g_xz_x_xx_xx[i] * a_exp;

        g_x_0_0_0_z_x_xx_xy[i] = 2.0 * g_xz_x_xx_xy[i] * a_exp;

        g_x_0_0_0_z_x_xx_xz[i] = 2.0 * g_xz_x_xx_xz[i] * a_exp;

        g_x_0_0_0_z_x_xx_yy[i] = 2.0 * g_xz_x_xx_yy[i] * a_exp;

        g_x_0_0_0_z_x_xx_yz[i] = 2.0 * g_xz_x_xx_yz[i] * a_exp;

        g_x_0_0_0_z_x_xx_zz[i] = 2.0 * g_xz_x_xx_zz[i] * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_0_0_0_z_x_xy_xx, g_x_0_0_0_z_x_xy_xy, g_x_0_0_0_z_x_xy_xz, g_x_0_0_0_z_x_xy_yy, g_x_0_0_0_z_x_xy_yz, g_x_0_0_0_z_x_xy_zz, g_xz_x_xy_xx, g_xz_x_xy_xy, g_xz_x_xy_xz, g_xz_x_xy_yy, g_xz_x_xy_yz, g_xz_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_x_xy_xx[i] = 2.0 * g_xz_x_xy_xx[i] * a_exp;

        g_x_0_0_0_z_x_xy_xy[i] = 2.0 * g_xz_x_xy_xy[i] * a_exp;

        g_x_0_0_0_z_x_xy_xz[i] = 2.0 * g_xz_x_xy_xz[i] * a_exp;

        g_x_0_0_0_z_x_xy_yy[i] = 2.0 * g_xz_x_xy_yy[i] * a_exp;

        g_x_0_0_0_z_x_xy_yz[i] = 2.0 * g_xz_x_xy_yz[i] * a_exp;

        g_x_0_0_0_z_x_xy_zz[i] = 2.0 * g_xz_x_xy_zz[i] * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_0_0_0_z_x_xz_xx, g_x_0_0_0_z_x_xz_xy, g_x_0_0_0_z_x_xz_xz, g_x_0_0_0_z_x_xz_yy, g_x_0_0_0_z_x_xz_yz, g_x_0_0_0_z_x_xz_zz, g_xz_x_xz_xx, g_xz_x_xz_xy, g_xz_x_xz_xz, g_xz_x_xz_yy, g_xz_x_xz_yz, g_xz_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_x_xz_xx[i] = 2.0 * g_xz_x_xz_xx[i] * a_exp;

        g_x_0_0_0_z_x_xz_xy[i] = 2.0 * g_xz_x_xz_xy[i] * a_exp;

        g_x_0_0_0_z_x_xz_xz[i] = 2.0 * g_xz_x_xz_xz[i] * a_exp;

        g_x_0_0_0_z_x_xz_yy[i] = 2.0 * g_xz_x_xz_yy[i] * a_exp;

        g_x_0_0_0_z_x_xz_yz[i] = 2.0 * g_xz_x_xz_yz[i] * a_exp;

        g_x_0_0_0_z_x_xz_zz[i] = 2.0 * g_xz_x_xz_zz[i] * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_0_0_0_z_x_yy_xx, g_x_0_0_0_z_x_yy_xy, g_x_0_0_0_z_x_yy_xz, g_x_0_0_0_z_x_yy_yy, g_x_0_0_0_z_x_yy_yz, g_x_0_0_0_z_x_yy_zz, g_xz_x_yy_xx, g_xz_x_yy_xy, g_xz_x_yy_xz, g_xz_x_yy_yy, g_xz_x_yy_yz, g_xz_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_x_yy_xx[i] = 2.0 * g_xz_x_yy_xx[i] * a_exp;

        g_x_0_0_0_z_x_yy_xy[i] = 2.0 * g_xz_x_yy_xy[i] * a_exp;

        g_x_0_0_0_z_x_yy_xz[i] = 2.0 * g_xz_x_yy_xz[i] * a_exp;

        g_x_0_0_0_z_x_yy_yy[i] = 2.0 * g_xz_x_yy_yy[i] * a_exp;

        g_x_0_0_0_z_x_yy_yz[i] = 2.0 * g_xz_x_yy_yz[i] * a_exp;

        g_x_0_0_0_z_x_yy_zz[i] = 2.0 * g_xz_x_yy_zz[i] * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_0_0_0_z_x_yz_xx, g_x_0_0_0_z_x_yz_xy, g_x_0_0_0_z_x_yz_xz, g_x_0_0_0_z_x_yz_yy, g_x_0_0_0_z_x_yz_yz, g_x_0_0_0_z_x_yz_zz, g_xz_x_yz_xx, g_xz_x_yz_xy, g_xz_x_yz_xz, g_xz_x_yz_yy, g_xz_x_yz_yz, g_xz_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_x_yz_xx[i] = 2.0 * g_xz_x_yz_xx[i] * a_exp;

        g_x_0_0_0_z_x_yz_xy[i] = 2.0 * g_xz_x_yz_xy[i] * a_exp;

        g_x_0_0_0_z_x_yz_xz[i] = 2.0 * g_xz_x_yz_xz[i] * a_exp;

        g_x_0_0_0_z_x_yz_yy[i] = 2.0 * g_xz_x_yz_yy[i] * a_exp;

        g_x_0_0_0_z_x_yz_yz[i] = 2.0 * g_xz_x_yz_yz[i] * a_exp;

        g_x_0_0_0_z_x_yz_zz[i] = 2.0 * g_xz_x_yz_zz[i] * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_0_0_0_z_x_zz_xx, g_x_0_0_0_z_x_zz_xy, g_x_0_0_0_z_x_zz_xz, g_x_0_0_0_z_x_zz_yy, g_x_0_0_0_z_x_zz_yz, g_x_0_0_0_z_x_zz_zz, g_xz_x_zz_xx, g_xz_x_zz_xy, g_xz_x_zz_xz, g_xz_x_zz_yy, g_xz_x_zz_yz, g_xz_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_x_zz_xx[i] = 2.0 * g_xz_x_zz_xx[i] * a_exp;

        g_x_0_0_0_z_x_zz_xy[i] = 2.0 * g_xz_x_zz_xy[i] * a_exp;

        g_x_0_0_0_z_x_zz_xz[i] = 2.0 * g_xz_x_zz_xz[i] * a_exp;

        g_x_0_0_0_z_x_zz_yy[i] = 2.0 * g_xz_x_zz_yy[i] * a_exp;

        g_x_0_0_0_z_x_zz_yz[i] = 2.0 * g_xz_x_zz_yz[i] * a_exp;

        g_x_0_0_0_z_x_zz_zz[i] = 2.0 * g_xz_x_zz_zz[i] * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_0_0_0_z_y_xx_xx, g_x_0_0_0_z_y_xx_xy, g_x_0_0_0_z_y_xx_xz, g_x_0_0_0_z_y_xx_yy, g_x_0_0_0_z_y_xx_yz, g_x_0_0_0_z_y_xx_zz, g_xz_y_xx_xx, g_xz_y_xx_xy, g_xz_y_xx_xz, g_xz_y_xx_yy, g_xz_y_xx_yz, g_xz_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_y_xx_xx[i] = 2.0 * g_xz_y_xx_xx[i] * a_exp;

        g_x_0_0_0_z_y_xx_xy[i] = 2.0 * g_xz_y_xx_xy[i] * a_exp;

        g_x_0_0_0_z_y_xx_xz[i] = 2.0 * g_xz_y_xx_xz[i] * a_exp;

        g_x_0_0_0_z_y_xx_yy[i] = 2.0 * g_xz_y_xx_yy[i] * a_exp;

        g_x_0_0_0_z_y_xx_yz[i] = 2.0 * g_xz_y_xx_yz[i] * a_exp;

        g_x_0_0_0_z_y_xx_zz[i] = 2.0 * g_xz_y_xx_zz[i] * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_0_0_0_z_y_xy_xx, g_x_0_0_0_z_y_xy_xy, g_x_0_0_0_z_y_xy_xz, g_x_0_0_0_z_y_xy_yy, g_x_0_0_0_z_y_xy_yz, g_x_0_0_0_z_y_xy_zz, g_xz_y_xy_xx, g_xz_y_xy_xy, g_xz_y_xy_xz, g_xz_y_xy_yy, g_xz_y_xy_yz, g_xz_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_y_xy_xx[i] = 2.0 * g_xz_y_xy_xx[i] * a_exp;

        g_x_0_0_0_z_y_xy_xy[i] = 2.0 * g_xz_y_xy_xy[i] * a_exp;

        g_x_0_0_0_z_y_xy_xz[i] = 2.0 * g_xz_y_xy_xz[i] * a_exp;

        g_x_0_0_0_z_y_xy_yy[i] = 2.0 * g_xz_y_xy_yy[i] * a_exp;

        g_x_0_0_0_z_y_xy_yz[i] = 2.0 * g_xz_y_xy_yz[i] * a_exp;

        g_x_0_0_0_z_y_xy_zz[i] = 2.0 * g_xz_y_xy_zz[i] * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_0_0_0_z_y_xz_xx, g_x_0_0_0_z_y_xz_xy, g_x_0_0_0_z_y_xz_xz, g_x_0_0_0_z_y_xz_yy, g_x_0_0_0_z_y_xz_yz, g_x_0_0_0_z_y_xz_zz, g_xz_y_xz_xx, g_xz_y_xz_xy, g_xz_y_xz_xz, g_xz_y_xz_yy, g_xz_y_xz_yz, g_xz_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_y_xz_xx[i] = 2.0 * g_xz_y_xz_xx[i] * a_exp;

        g_x_0_0_0_z_y_xz_xy[i] = 2.0 * g_xz_y_xz_xy[i] * a_exp;

        g_x_0_0_0_z_y_xz_xz[i] = 2.0 * g_xz_y_xz_xz[i] * a_exp;

        g_x_0_0_0_z_y_xz_yy[i] = 2.0 * g_xz_y_xz_yy[i] * a_exp;

        g_x_0_0_0_z_y_xz_yz[i] = 2.0 * g_xz_y_xz_yz[i] * a_exp;

        g_x_0_0_0_z_y_xz_zz[i] = 2.0 * g_xz_y_xz_zz[i] * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_0_0_0_z_y_yy_xx, g_x_0_0_0_z_y_yy_xy, g_x_0_0_0_z_y_yy_xz, g_x_0_0_0_z_y_yy_yy, g_x_0_0_0_z_y_yy_yz, g_x_0_0_0_z_y_yy_zz, g_xz_y_yy_xx, g_xz_y_yy_xy, g_xz_y_yy_xz, g_xz_y_yy_yy, g_xz_y_yy_yz, g_xz_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_y_yy_xx[i] = 2.0 * g_xz_y_yy_xx[i] * a_exp;

        g_x_0_0_0_z_y_yy_xy[i] = 2.0 * g_xz_y_yy_xy[i] * a_exp;

        g_x_0_0_0_z_y_yy_xz[i] = 2.0 * g_xz_y_yy_xz[i] * a_exp;

        g_x_0_0_0_z_y_yy_yy[i] = 2.0 * g_xz_y_yy_yy[i] * a_exp;

        g_x_0_0_0_z_y_yy_yz[i] = 2.0 * g_xz_y_yy_yz[i] * a_exp;

        g_x_0_0_0_z_y_yy_zz[i] = 2.0 * g_xz_y_yy_zz[i] * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_0_0_0_z_y_yz_xx, g_x_0_0_0_z_y_yz_xy, g_x_0_0_0_z_y_yz_xz, g_x_0_0_0_z_y_yz_yy, g_x_0_0_0_z_y_yz_yz, g_x_0_0_0_z_y_yz_zz, g_xz_y_yz_xx, g_xz_y_yz_xy, g_xz_y_yz_xz, g_xz_y_yz_yy, g_xz_y_yz_yz, g_xz_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_y_yz_xx[i] = 2.0 * g_xz_y_yz_xx[i] * a_exp;

        g_x_0_0_0_z_y_yz_xy[i] = 2.0 * g_xz_y_yz_xy[i] * a_exp;

        g_x_0_0_0_z_y_yz_xz[i] = 2.0 * g_xz_y_yz_xz[i] * a_exp;

        g_x_0_0_0_z_y_yz_yy[i] = 2.0 * g_xz_y_yz_yy[i] * a_exp;

        g_x_0_0_0_z_y_yz_yz[i] = 2.0 * g_xz_y_yz_yz[i] * a_exp;

        g_x_0_0_0_z_y_yz_zz[i] = 2.0 * g_xz_y_yz_zz[i] * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_0_0_0_z_y_zz_xx, g_x_0_0_0_z_y_zz_xy, g_x_0_0_0_z_y_zz_xz, g_x_0_0_0_z_y_zz_yy, g_x_0_0_0_z_y_zz_yz, g_x_0_0_0_z_y_zz_zz, g_xz_y_zz_xx, g_xz_y_zz_xy, g_xz_y_zz_xz, g_xz_y_zz_yy, g_xz_y_zz_yz, g_xz_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_y_zz_xx[i] = 2.0 * g_xz_y_zz_xx[i] * a_exp;

        g_x_0_0_0_z_y_zz_xy[i] = 2.0 * g_xz_y_zz_xy[i] * a_exp;

        g_x_0_0_0_z_y_zz_xz[i] = 2.0 * g_xz_y_zz_xz[i] * a_exp;

        g_x_0_0_0_z_y_zz_yy[i] = 2.0 * g_xz_y_zz_yy[i] * a_exp;

        g_x_0_0_0_z_y_zz_yz[i] = 2.0 * g_xz_y_zz_yz[i] * a_exp;

        g_x_0_0_0_z_y_zz_zz[i] = 2.0 * g_xz_y_zz_zz[i] * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_0_0_0_z_z_xx_xx, g_x_0_0_0_z_z_xx_xy, g_x_0_0_0_z_z_xx_xz, g_x_0_0_0_z_z_xx_yy, g_x_0_0_0_z_z_xx_yz, g_x_0_0_0_z_z_xx_zz, g_xz_z_xx_xx, g_xz_z_xx_xy, g_xz_z_xx_xz, g_xz_z_xx_yy, g_xz_z_xx_yz, g_xz_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_z_xx_xx[i] = 2.0 * g_xz_z_xx_xx[i] * a_exp;

        g_x_0_0_0_z_z_xx_xy[i] = 2.0 * g_xz_z_xx_xy[i] * a_exp;

        g_x_0_0_0_z_z_xx_xz[i] = 2.0 * g_xz_z_xx_xz[i] * a_exp;

        g_x_0_0_0_z_z_xx_yy[i] = 2.0 * g_xz_z_xx_yy[i] * a_exp;

        g_x_0_0_0_z_z_xx_yz[i] = 2.0 * g_xz_z_xx_yz[i] * a_exp;

        g_x_0_0_0_z_z_xx_zz[i] = 2.0 * g_xz_z_xx_zz[i] * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_0_0_0_z_z_xy_xx, g_x_0_0_0_z_z_xy_xy, g_x_0_0_0_z_z_xy_xz, g_x_0_0_0_z_z_xy_yy, g_x_0_0_0_z_z_xy_yz, g_x_0_0_0_z_z_xy_zz, g_xz_z_xy_xx, g_xz_z_xy_xy, g_xz_z_xy_xz, g_xz_z_xy_yy, g_xz_z_xy_yz, g_xz_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_z_xy_xx[i] = 2.0 * g_xz_z_xy_xx[i] * a_exp;

        g_x_0_0_0_z_z_xy_xy[i] = 2.0 * g_xz_z_xy_xy[i] * a_exp;

        g_x_0_0_0_z_z_xy_xz[i] = 2.0 * g_xz_z_xy_xz[i] * a_exp;

        g_x_0_0_0_z_z_xy_yy[i] = 2.0 * g_xz_z_xy_yy[i] * a_exp;

        g_x_0_0_0_z_z_xy_yz[i] = 2.0 * g_xz_z_xy_yz[i] * a_exp;

        g_x_0_0_0_z_z_xy_zz[i] = 2.0 * g_xz_z_xy_zz[i] * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_0_0_0_z_z_xz_xx, g_x_0_0_0_z_z_xz_xy, g_x_0_0_0_z_z_xz_xz, g_x_0_0_0_z_z_xz_yy, g_x_0_0_0_z_z_xz_yz, g_x_0_0_0_z_z_xz_zz, g_xz_z_xz_xx, g_xz_z_xz_xy, g_xz_z_xz_xz, g_xz_z_xz_yy, g_xz_z_xz_yz, g_xz_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_z_xz_xx[i] = 2.0 * g_xz_z_xz_xx[i] * a_exp;

        g_x_0_0_0_z_z_xz_xy[i] = 2.0 * g_xz_z_xz_xy[i] * a_exp;

        g_x_0_0_0_z_z_xz_xz[i] = 2.0 * g_xz_z_xz_xz[i] * a_exp;

        g_x_0_0_0_z_z_xz_yy[i] = 2.0 * g_xz_z_xz_yy[i] * a_exp;

        g_x_0_0_0_z_z_xz_yz[i] = 2.0 * g_xz_z_xz_yz[i] * a_exp;

        g_x_0_0_0_z_z_xz_zz[i] = 2.0 * g_xz_z_xz_zz[i] * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_0_0_0_z_z_yy_xx, g_x_0_0_0_z_z_yy_xy, g_x_0_0_0_z_z_yy_xz, g_x_0_0_0_z_z_yy_yy, g_x_0_0_0_z_z_yy_yz, g_x_0_0_0_z_z_yy_zz, g_xz_z_yy_xx, g_xz_z_yy_xy, g_xz_z_yy_xz, g_xz_z_yy_yy, g_xz_z_yy_yz, g_xz_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_z_yy_xx[i] = 2.0 * g_xz_z_yy_xx[i] * a_exp;

        g_x_0_0_0_z_z_yy_xy[i] = 2.0 * g_xz_z_yy_xy[i] * a_exp;

        g_x_0_0_0_z_z_yy_xz[i] = 2.0 * g_xz_z_yy_xz[i] * a_exp;

        g_x_0_0_0_z_z_yy_yy[i] = 2.0 * g_xz_z_yy_yy[i] * a_exp;

        g_x_0_0_0_z_z_yy_yz[i] = 2.0 * g_xz_z_yy_yz[i] * a_exp;

        g_x_0_0_0_z_z_yy_zz[i] = 2.0 * g_xz_z_yy_zz[i] * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_0_0_0_z_z_yz_xx, g_x_0_0_0_z_z_yz_xy, g_x_0_0_0_z_z_yz_xz, g_x_0_0_0_z_z_yz_yy, g_x_0_0_0_z_z_yz_yz, g_x_0_0_0_z_z_yz_zz, g_xz_z_yz_xx, g_xz_z_yz_xy, g_xz_z_yz_xz, g_xz_z_yz_yy, g_xz_z_yz_yz, g_xz_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_z_yz_xx[i] = 2.0 * g_xz_z_yz_xx[i] * a_exp;

        g_x_0_0_0_z_z_yz_xy[i] = 2.0 * g_xz_z_yz_xy[i] * a_exp;

        g_x_0_0_0_z_z_yz_xz[i] = 2.0 * g_xz_z_yz_xz[i] * a_exp;

        g_x_0_0_0_z_z_yz_yy[i] = 2.0 * g_xz_z_yz_yy[i] * a_exp;

        g_x_0_0_0_z_z_yz_yz[i] = 2.0 * g_xz_z_yz_yz[i] * a_exp;

        g_x_0_0_0_z_z_yz_zz[i] = 2.0 * g_xz_z_yz_zz[i] * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_0_0_0_z_z_zz_xx, g_x_0_0_0_z_z_zz_xy, g_x_0_0_0_z_z_zz_xz, g_x_0_0_0_z_z_zz_yy, g_x_0_0_0_z_z_zz_yz, g_x_0_0_0_z_z_zz_zz, g_xz_z_zz_xx, g_xz_z_zz_xy, g_xz_z_zz_xz, g_xz_z_zz_yy, g_xz_z_zz_yz, g_xz_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_z_zz_xx[i] = 2.0 * g_xz_z_zz_xx[i] * a_exp;

        g_x_0_0_0_z_z_zz_xy[i] = 2.0 * g_xz_z_zz_xy[i] * a_exp;

        g_x_0_0_0_z_z_zz_xz[i] = 2.0 * g_xz_z_zz_xz[i] * a_exp;

        g_x_0_0_0_z_z_zz_yy[i] = 2.0 * g_xz_z_zz_yy[i] * a_exp;

        g_x_0_0_0_z_z_zz_yz[i] = 2.0 * g_xz_z_zz_yz[i] * a_exp;

        g_x_0_0_0_z_z_zz_zz[i] = 2.0 * g_xz_z_zz_zz[i] * a_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_xy_x_xx_xx, g_xy_x_xx_xy, g_xy_x_xx_xz, g_xy_x_xx_yy, g_xy_x_xx_yz, g_xy_x_xx_zz, g_y_0_0_0_x_x_xx_xx, g_y_0_0_0_x_x_xx_xy, g_y_0_0_0_x_x_xx_xz, g_y_0_0_0_x_x_xx_yy, g_y_0_0_0_x_x_xx_yz, g_y_0_0_0_x_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_x_xx_xx[i] = 2.0 * g_xy_x_xx_xx[i] * a_exp;

        g_y_0_0_0_x_x_xx_xy[i] = 2.0 * g_xy_x_xx_xy[i] * a_exp;

        g_y_0_0_0_x_x_xx_xz[i] = 2.0 * g_xy_x_xx_xz[i] * a_exp;

        g_y_0_0_0_x_x_xx_yy[i] = 2.0 * g_xy_x_xx_yy[i] * a_exp;

        g_y_0_0_0_x_x_xx_yz[i] = 2.0 * g_xy_x_xx_yz[i] * a_exp;

        g_y_0_0_0_x_x_xx_zz[i] = 2.0 * g_xy_x_xx_zz[i] * a_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_xy_x_xy_xx, g_xy_x_xy_xy, g_xy_x_xy_xz, g_xy_x_xy_yy, g_xy_x_xy_yz, g_xy_x_xy_zz, g_y_0_0_0_x_x_xy_xx, g_y_0_0_0_x_x_xy_xy, g_y_0_0_0_x_x_xy_xz, g_y_0_0_0_x_x_xy_yy, g_y_0_0_0_x_x_xy_yz, g_y_0_0_0_x_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_x_xy_xx[i] = 2.0 * g_xy_x_xy_xx[i] * a_exp;

        g_y_0_0_0_x_x_xy_xy[i] = 2.0 * g_xy_x_xy_xy[i] * a_exp;

        g_y_0_0_0_x_x_xy_xz[i] = 2.0 * g_xy_x_xy_xz[i] * a_exp;

        g_y_0_0_0_x_x_xy_yy[i] = 2.0 * g_xy_x_xy_yy[i] * a_exp;

        g_y_0_0_0_x_x_xy_yz[i] = 2.0 * g_xy_x_xy_yz[i] * a_exp;

        g_y_0_0_0_x_x_xy_zz[i] = 2.0 * g_xy_x_xy_zz[i] * a_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_xy_x_xz_xx, g_xy_x_xz_xy, g_xy_x_xz_xz, g_xy_x_xz_yy, g_xy_x_xz_yz, g_xy_x_xz_zz, g_y_0_0_0_x_x_xz_xx, g_y_0_0_0_x_x_xz_xy, g_y_0_0_0_x_x_xz_xz, g_y_0_0_0_x_x_xz_yy, g_y_0_0_0_x_x_xz_yz, g_y_0_0_0_x_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_x_xz_xx[i] = 2.0 * g_xy_x_xz_xx[i] * a_exp;

        g_y_0_0_0_x_x_xz_xy[i] = 2.0 * g_xy_x_xz_xy[i] * a_exp;

        g_y_0_0_0_x_x_xz_xz[i] = 2.0 * g_xy_x_xz_xz[i] * a_exp;

        g_y_0_0_0_x_x_xz_yy[i] = 2.0 * g_xy_x_xz_yy[i] * a_exp;

        g_y_0_0_0_x_x_xz_yz[i] = 2.0 * g_xy_x_xz_yz[i] * a_exp;

        g_y_0_0_0_x_x_xz_zz[i] = 2.0 * g_xy_x_xz_zz[i] * a_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_xy_x_yy_xx, g_xy_x_yy_xy, g_xy_x_yy_xz, g_xy_x_yy_yy, g_xy_x_yy_yz, g_xy_x_yy_zz, g_y_0_0_0_x_x_yy_xx, g_y_0_0_0_x_x_yy_xy, g_y_0_0_0_x_x_yy_xz, g_y_0_0_0_x_x_yy_yy, g_y_0_0_0_x_x_yy_yz, g_y_0_0_0_x_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_x_yy_xx[i] = 2.0 * g_xy_x_yy_xx[i] * a_exp;

        g_y_0_0_0_x_x_yy_xy[i] = 2.0 * g_xy_x_yy_xy[i] * a_exp;

        g_y_0_0_0_x_x_yy_xz[i] = 2.0 * g_xy_x_yy_xz[i] * a_exp;

        g_y_0_0_0_x_x_yy_yy[i] = 2.0 * g_xy_x_yy_yy[i] * a_exp;

        g_y_0_0_0_x_x_yy_yz[i] = 2.0 * g_xy_x_yy_yz[i] * a_exp;

        g_y_0_0_0_x_x_yy_zz[i] = 2.0 * g_xy_x_yy_zz[i] * a_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_xy_x_yz_xx, g_xy_x_yz_xy, g_xy_x_yz_xz, g_xy_x_yz_yy, g_xy_x_yz_yz, g_xy_x_yz_zz, g_y_0_0_0_x_x_yz_xx, g_y_0_0_0_x_x_yz_xy, g_y_0_0_0_x_x_yz_xz, g_y_0_0_0_x_x_yz_yy, g_y_0_0_0_x_x_yz_yz, g_y_0_0_0_x_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_x_yz_xx[i] = 2.0 * g_xy_x_yz_xx[i] * a_exp;

        g_y_0_0_0_x_x_yz_xy[i] = 2.0 * g_xy_x_yz_xy[i] * a_exp;

        g_y_0_0_0_x_x_yz_xz[i] = 2.0 * g_xy_x_yz_xz[i] * a_exp;

        g_y_0_0_0_x_x_yz_yy[i] = 2.0 * g_xy_x_yz_yy[i] * a_exp;

        g_y_0_0_0_x_x_yz_yz[i] = 2.0 * g_xy_x_yz_yz[i] * a_exp;

        g_y_0_0_0_x_x_yz_zz[i] = 2.0 * g_xy_x_yz_zz[i] * a_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_xy_x_zz_xx, g_xy_x_zz_xy, g_xy_x_zz_xz, g_xy_x_zz_yy, g_xy_x_zz_yz, g_xy_x_zz_zz, g_y_0_0_0_x_x_zz_xx, g_y_0_0_0_x_x_zz_xy, g_y_0_0_0_x_x_zz_xz, g_y_0_0_0_x_x_zz_yy, g_y_0_0_0_x_x_zz_yz, g_y_0_0_0_x_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_x_zz_xx[i] = 2.0 * g_xy_x_zz_xx[i] * a_exp;

        g_y_0_0_0_x_x_zz_xy[i] = 2.0 * g_xy_x_zz_xy[i] * a_exp;

        g_y_0_0_0_x_x_zz_xz[i] = 2.0 * g_xy_x_zz_xz[i] * a_exp;

        g_y_0_0_0_x_x_zz_yy[i] = 2.0 * g_xy_x_zz_yy[i] * a_exp;

        g_y_0_0_0_x_x_zz_yz[i] = 2.0 * g_xy_x_zz_yz[i] * a_exp;

        g_y_0_0_0_x_x_zz_zz[i] = 2.0 * g_xy_x_zz_zz[i] * a_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_xy_y_xx_xx, g_xy_y_xx_xy, g_xy_y_xx_xz, g_xy_y_xx_yy, g_xy_y_xx_yz, g_xy_y_xx_zz, g_y_0_0_0_x_y_xx_xx, g_y_0_0_0_x_y_xx_xy, g_y_0_0_0_x_y_xx_xz, g_y_0_0_0_x_y_xx_yy, g_y_0_0_0_x_y_xx_yz, g_y_0_0_0_x_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_y_xx_xx[i] = 2.0 * g_xy_y_xx_xx[i] * a_exp;

        g_y_0_0_0_x_y_xx_xy[i] = 2.0 * g_xy_y_xx_xy[i] * a_exp;

        g_y_0_0_0_x_y_xx_xz[i] = 2.0 * g_xy_y_xx_xz[i] * a_exp;

        g_y_0_0_0_x_y_xx_yy[i] = 2.0 * g_xy_y_xx_yy[i] * a_exp;

        g_y_0_0_0_x_y_xx_yz[i] = 2.0 * g_xy_y_xx_yz[i] * a_exp;

        g_y_0_0_0_x_y_xx_zz[i] = 2.0 * g_xy_y_xx_zz[i] * a_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_xy_y_xy_xx, g_xy_y_xy_xy, g_xy_y_xy_xz, g_xy_y_xy_yy, g_xy_y_xy_yz, g_xy_y_xy_zz, g_y_0_0_0_x_y_xy_xx, g_y_0_0_0_x_y_xy_xy, g_y_0_0_0_x_y_xy_xz, g_y_0_0_0_x_y_xy_yy, g_y_0_0_0_x_y_xy_yz, g_y_0_0_0_x_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_y_xy_xx[i] = 2.0 * g_xy_y_xy_xx[i] * a_exp;

        g_y_0_0_0_x_y_xy_xy[i] = 2.0 * g_xy_y_xy_xy[i] * a_exp;

        g_y_0_0_0_x_y_xy_xz[i] = 2.0 * g_xy_y_xy_xz[i] * a_exp;

        g_y_0_0_0_x_y_xy_yy[i] = 2.0 * g_xy_y_xy_yy[i] * a_exp;

        g_y_0_0_0_x_y_xy_yz[i] = 2.0 * g_xy_y_xy_yz[i] * a_exp;

        g_y_0_0_0_x_y_xy_zz[i] = 2.0 * g_xy_y_xy_zz[i] * a_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_xy_y_xz_xx, g_xy_y_xz_xy, g_xy_y_xz_xz, g_xy_y_xz_yy, g_xy_y_xz_yz, g_xy_y_xz_zz, g_y_0_0_0_x_y_xz_xx, g_y_0_0_0_x_y_xz_xy, g_y_0_0_0_x_y_xz_xz, g_y_0_0_0_x_y_xz_yy, g_y_0_0_0_x_y_xz_yz, g_y_0_0_0_x_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_y_xz_xx[i] = 2.0 * g_xy_y_xz_xx[i] * a_exp;

        g_y_0_0_0_x_y_xz_xy[i] = 2.0 * g_xy_y_xz_xy[i] * a_exp;

        g_y_0_0_0_x_y_xz_xz[i] = 2.0 * g_xy_y_xz_xz[i] * a_exp;

        g_y_0_0_0_x_y_xz_yy[i] = 2.0 * g_xy_y_xz_yy[i] * a_exp;

        g_y_0_0_0_x_y_xz_yz[i] = 2.0 * g_xy_y_xz_yz[i] * a_exp;

        g_y_0_0_0_x_y_xz_zz[i] = 2.0 * g_xy_y_xz_zz[i] * a_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_xy_y_yy_xx, g_xy_y_yy_xy, g_xy_y_yy_xz, g_xy_y_yy_yy, g_xy_y_yy_yz, g_xy_y_yy_zz, g_y_0_0_0_x_y_yy_xx, g_y_0_0_0_x_y_yy_xy, g_y_0_0_0_x_y_yy_xz, g_y_0_0_0_x_y_yy_yy, g_y_0_0_0_x_y_yy_yz, g_y_0_0_0_x_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_y_yy_xx[i] = 2.0 * g_xy_y_yy_xx[i] * a_exp;

        g_y_0_0_0_x_y_yy_xy[i] = 2.0 * g_xy_y_yy_xy[i] * a_exp;

        g_y_0_0_0_x_y_yy_xz[i] = 2.0 * g_xy_y_yy_xz[i] * a_exp;

        g_y_0_0_0_x_y_yy_yy[i] = 2.0 * g_xy_y_yy_yy[i] * a_exp;

        g_y_0_0_0_x_y_yy_yz[i] = 2.0 * g_xy_y_yy_yz[i] * a_exp;

        g_y_0_0_0_x_y_yy_zz[i] = 2.0 * g_xy_y_yy_zz[i] * a_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_xy_y_yz_xx, g_xy_y_yz_xy, g_xy_y_yz_xz, g_xy_y_yz_yy, g_xy_y_yz_yz, g_xy_y_yz_zz, g_y_0_0_0_x_y_yz_xx, g_y_0_0_0_x_y_yz_xy, g_y_0_0_0_x_y_yz_xz, g_y_0_0_0_x_y_yz_yy, g_y_0_0_0_x_y_yz_yz, g_y_0_0_0_x_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_y_yz_xx[i] = 2.0 * g_xy_y_yz_xx[i] * a_exp;

        g_y_0_0_0_x_y_yz_xy[i] = 2.0 * g_xy_y_yz_xy[i] * a_exp;

        g_y_0_0_0_x_y_yz_xz[i] = 2.0 * g_xy_y_yz_xz[i] * a_exp;

        g_y_0_0_0_x_y_yz_yy[i] = 2.0 * g_xy_y_yz_yy[i] * a_exp;

        g_y_0_0_0_x_y_yz_yz[i] = 2.0 * g_xy_y_yz_yz[i] * a_exp;

        g_y_0_0_0_x_y_yz_zz[i] = 2.0 * g_xy_y_yz_zz[i] * a_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_xy_y_zz_xx, g_xy_y_zz_xy, g_xy_y_zz_xz, g_xy_y_zz_yy, g_xy_y_zz_yz, g_xy_y_zz_zz, g_y_0_0_0_x_y_zz_xx, g_y_0_0_0_x_y_zz_xy, g_y_0_0_0_x_y_zz_xz, g_y_0_0_0_x_y_zz_yy, g_y_0_0_0_x_y_zz_yz, g_y_0_0_0_x_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_y_zz_xx[i] = 2.0 * g_xy_y_zz_xx[i] * a_exp;

        g_y_0_0_0_x_y_zz_xy[i] = 2.0 * g_xy_y_zz_xy[i] * a_exp;

        g_y_0_0_0_x_y_zz_xz[i] = 2.0 * g_xy_y_zz_xz[i] * a_exp;

        g_y_0_0_0_x_y_zz_yy[i] = 2.0 * g_xy_y_zz_yy[i] * a_exp;

        g_y_0_0_0_x_y_zz_yz[i] = 2.0 * g_xy_y_zz_yz[i] * a_exp;

        g_y_0_0_0_x_y_zz_zz[i] = 2.0 * g_xy_y_zz_zz[i] * a_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_xy_z_xx_xx, g_xy_z_xx_xy, g_xy_z_xx_xz, g_xy_z_xx_yy, g_xy_z_xx_yz, g_xy_z_xx_zz, g_y_0_0_0_x_z_xx_xx, g_y_0_0_0_x_z_xx_xy, g_y_0_0_0_x_z_xx_xz, g_y_0_0_0_x_z_xx_yy, g_y_0_0_0_x_z_xx_yz, g_y_0_0_0_x_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_z_xx_xx[i] = 2.0 * g_xy_z_xx_xx[i] * a_exp;

        g_y_0_0_0_x_z_xx_xy[i] = 2.0 * g_xy_z_xx_xy[i] * a_exp;

        g_y_0_0_0_x_z_xx_xz[i] = 2.0 * g_xy_z_xx_xz[i] * a_exp;

        g_y_0_0_0_x_z_xx_yy[i] = 2.0 * g_xy_z_xx_yy[i] * a_exp;

        g_y_0_0_0_x_z_xx_yz[i] = 2.0 * g_xy_z_xx_yz[i] * a_exp;

        g_y_0_0_0_x_z_xx_zz[i] = 2.0 * g_xy_z_xx_zz[i] * a_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_xy_z_xy_xx, g_xy_z_xy_xy, g_xy_z_xy_xz, g_xy_z_xy_yy, g_xy_z_xy_yz, g_xy_z_xy_zz, g_y_0_0_0_x_z_xy_xx, g_y_0_0_0_x_z_xy_xy, g_y_0_0_0_x_z_xy_xz, g_y_0_0_0_x_z_xy_yy, g_y_0_0_0_x_z_xy_yz, g_y_0_0_0_x_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_z_xy_xx[i] = 2.0 * g_xy_z_xy_xx[i] * a_exp;

        g_y_0_0_0_x_z_xy_xy[i] = 2.0 * g_xy_z_xy_xy[i] * a_exp;

        g_y_0_0_0_x_z_xy_xz[i] = 2.0 * g_xy_z_xy_xz[i] * a_exp;

        g_y_0_0_0_x_z_xy_yy[i] = 2.0 * g_xy_z_xy_yy[i] * a_exp;

        g_y_0_0_0_x_z_xy_yz[i] = 2.0 * g_xy_z_xy_yz[i] * a_exp;

        g_y_0_0_0_x_z_xy_zz[i] = 2.0 * g_xy_z_xy_zz[i] * a_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_xy_z_xz_xx, g_xy_z_xz_xy, g_xy_z_xz_xz, g_xy_z_xz_yy, g_xy_z_xz_yz, g_xy_z_xz_zz, g_y_0_0_0_x_z_xz_xx, g_y_0_0_0_x_z_xz_xy, g_y_0_0_0_x_z_xz_xz, g_y_0_0_0_x_z_xz_yy, g_y_0_0_0_x_z_xz_yz, g_y_0_0_0_x_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_z_xz_xx[i] = 2.0 * g_xy_z_xz_xx[i] * a_exp;

        g_y_0_0_0_x_z_xz_xy[i] = 2.0 * g_xy_z_xz_xy[i] * a_exp;

        g_y_0_0_0_x_z_xz_xz[i] = 2.0 * g_xy_z_xz_xz[i] * a_exp;

        g_y_0_0_0_x_z_xz_yy[i] = 2.0 * g_xy_z_xz_yy[i] * a_exp;

        g_y_0_0_0_x_z_xz_yz[i] = 2.0 * g_xy_z_xz_yz[i] * a_exp;

        g_y_0_0_0_x_z_xz_zz[i] = 2.0 * g_xy_z_xz_zz[i] * a_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_xy_z_yy_xx, g_xy_z_yy_xy, g_xy_z_yy_xz, g_xy_z_yy_yy, g_xy_z_yy_yz, g_xy_z_yy_zz, g_y_0_0_0_x_z_yy_xx, g_y_0_0_0_x_z_yy_xy, g_y_0_0_0_x_z_yy_xz, g_y_0_0_0_x_z_yy_yy, g_y_0_0_0_x_z_yy_yz, g_y_0_0_0_x_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_z_yy_xx[i] = 2.0 * g_xy_z_yy_xx[i] * a_exp;

        g_y_0_0_0_x_z_yy_xy[i] = 2.0 * g_xy_z_yy_xy[i] * a_exp;

        g_y_0_0_0_x_z_yy_xz[i] = 2.0 * g_xy_z_yy_xz[i] * a_exp;

        g_y_0_0_0_x_z_yy_yy[i] = 2.0 * g_xy_z_yy_yy[i] * a_exp;

        g_y_0_0_0_x_z_yy_yz[i] = 2.0 * g_xy_z_yy_yz[i] * a_exp;

        g_y_0_0_0_x_z_yy_zz[i] = 2.0 * g_xy_z_yy_zz[i] * a_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_xy_z_yz_xx, g_xy_z_yz_xy, g_xy_z_yz_xz, g_xy_z_yz_yy, g_xy_z_yz_yz, g_xy_z_yz_zz, g_y_0_0_0_x_z_yz_xx, g_y_0_0_0_x_z_yz_xy, g_y_0_0_0_x_z_yz_xz, g_y_0_0_0_x_z_yz_yy, g_y_0_0_0_x_z_yz_yz, g_y_0_0_0_x_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_z_yz_xx[i] = 2.0 * g_xy_z_yz_xx[i] * a_exp;

        g_y_0_0_0_x_z_yz_xy[i] = 2.0 * g_xy_z_yz_xy[i] * a_exp;

        g_y_0_0_0_x_z_yz_xz[i] = 2.0 * g_xy_z_yz_xz[i] * a_exp;

        g_y_0_0_0_x_z_yz_yy[i] = 2.0 * g_xy_z_yz_yy[i] * a_exp;

        g_y_0_0_0_x_z_yz_yz[i] = 2.0 * g_xy_z_yz_yz[i] * a_exp;

        g_y_0_0_0_x_z_yz_zz[i] = 2.0 * g_xy_z_yz_zz[i] * a_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_xy_z_zz_xx, g_xy_z_zz_xy, g_xy_z_zz_xz, g_xy_z_zz_yy, g_xy_z_zz_yz, g_xy_z_zz_zz, g_y_0_0_0_x_z_zz_xx, g_y_0_0_0_x_z_zz_xy, g_y_0_0_0_x_z_zz_xz, g_y_0_0_0_x_z_zz_yy, g_y_0_0_0_x_z_zz_yz, g_y_0_0_0_x_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_z_zz_xx[i] = 2.0 * g_xy_z_zz_xx[i] * a_exp;

        g_y_0_0_0_x_z_zz_xy[i] = 2.0 * g_xy_z_zz_xy[i] * a_exp;

        g_y_0_0_0_x_z_zz_xz[i] = 2.0 * g_xy_z_zz_xz[i] * a_exp;

        g_y_0_0_0_x_z_zz_yy[i] = 2.0 * g_xy_z_zz_yy[i] * a_exp;

        g_y_0_0_0_x_z_zz_yz[i] = 2.0 * g_xy_z_zz_yz[i] * a_exp;

        g_y_0_0_0_x_z_zz_zz[i] = 2.0 * g_xy_z_zz_zz[i] * a_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_0_x_xx_xx, g_0_x_xx_xy, g_0_x_xx_xz, g_0_x_xx_yy, g_0_x_xx_yz, g_0_x_xx_zz, g_y_0_0_0_y_x_xx_xx, g_y_0_0_0_y_x_xx_xy, g_y_0_0_0_y_x_xx_xz, g_y_0_0_0_y_x_xx_yy, g_y_0_0_0_y_x_xx_yz, g_y_0_0_0_y_x_xx_zz, g_yy_x_xx_xx, g_yy_x_xx_xy, g_yy_x_xx_xz, g_yy_x_xx_yy, g_yy_x_xx_yz, g_yy_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_x_xx_xx[i] = -g_0_x_xx_xx[i] + 2.0 * g_yy_x_xx_xx[i] * a_exp;

        g_y_0_0_0_y_x_xx_xy[i] = -g_0_x_xx_xy[i] + 2.0 * g_yy_x_xx_xy[i] * a_exp;

        g_y_0_0_0_y_x_xx_xz[i] = -g_0_x_xx_xz[i] + 2.0 * g_yy_x_xx_xz[i] * a_exp;

        g_y_0_0_0_y_x_xx_yy[i] = -g_0_x_xx_yy[i] + 2.0 * g_yy_x_xx_yy[i] * a_exp;

        g_y_0_0_0_y_x_xx_yz[i] = -g_0_x_xx_yz[i] + 2.0 * g_yy_x_xx_yz[i] * a_exp;

        g_y_0_0_0_y_x_xx_zz[i] = -g_0_x_xx_zz[i] + 2.0 * g_yy_x_xx_zz[i] * a_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_0_x_xy_xx, g_0_x_xy_xy, g_0_x_xy_xz, g_0_x_xy_yy, g_0_x_xy_yz, g_0_x_xy_zz, g_y_0_0_0_y_x_xy_xx, g_y_0_0_0_y_x_xy_xy, g_y_0_0_0_y_x_xy_xz, g_y_0_0_0_y_x_xy_yy, g_y_0_0_0_y_x_xy_yz, g_y_0_0_0_y_x_xy_zz, g_yy_x_xy_xx, g_yy_x_xy_xy, g_yy_x_xy_xz, g_yy_x_xy_yy, g_yy_x_xy_yz, g_yy_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_x_xy_xx[i] = -g_0_x_xy_xx[i] + 2.0 * g_yy_x_xy_xx[i] * a_exp;

        g_y_0_0_0_y_x_xy_xy[i] = -g_0_x_xy_xy[i] + 2.0 * g_yy_x_xy_xy[i] * a_exp;

        g_y_0_0_0_y_x_xy_xz[i] = -g_0_x_xy_xz[i] + 2.0 * g_yy_x_xy_xz[i] * a_exp;

        g_y_0_0_0_y_x_xy_yy[i] = -g_0_x_xy_yy[i] + 2.0 * g_yy_x_xy_yy[i] * a_exp;

        g_y_0_0_0_y_x_xy_yz[i] = -g_0_x_xy_yz[i] + 2.0 * g_yy_x_xy_yz[i] * a_exp;

        g_y_0_0_0_y_x_xy_zz[i] = -g_0_x_xy_zz[i] + 2.0 * g_yy_x_xy_zz[i] * a_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_0_x_xz_xx, g_0_x_xz_xy, g_0_x_xz_xz, g_0_x_xz_yy, g_0_x_xz_yz, g_0_x_xz_zz, g_y_0_0_0_y_x_xz_xx, g_y_0_0_0_y_x_xz_xy, g_y_0_0_0_y_x_xz_xz, g_y_0_0_0_y_x_xz_yy, g_y_0_0_0_y_x_xz_yz, g_y_0_0_0_y_x_xz_zz, g_yy_x_xz_xx, g_yy_x_xz_xy, g_yy_x_xz_xz, g_yy_x_xz_yy, g_yy_x_xz_yz, g_yy_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_x_xz_xx[i] = -g_0_x_xz_xx[i] + 2.0 * g_yy_x_xz_xx[i] * a_exp;

        g_y_0_0_0_y_x_xz_xy[i] = -g_0_x_xz_xy[i] + 2.0 * g_yy_x_xz_xy[i] * a_exp;

        g_y_0_0_0_y_x_xz_xz[i] = -g_0_x_xz_xz[i] + 2.0 * g_yy_x_xz_xz[i] * a_exp;

        g_y_0_0_0_y_x_xz_yy[i] = -g_0_x_xz_yy[i] + 2.0 * g_yy_x_xz_yy[i] * a_exp;

        g_y_0_0_0_y_x_xz_yz[i] = -g_0_x_xz_yz[i] + 2.0 * g_yy_x_xz_yz[i] * a_exp;

        g_y_0_0_0_y_x_xz_zz[i] = -g_0_x_xz_zz[i] + 2.0 * g_yy_x_xz_zz[i] * a_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_0_x_yy_xx, g_0_x_yy_xy, g_0_x_yy_xz, g_0_x_yy_yy, g_0_x_yy_yz, g_0_x_yy_zz, g_y_0_0_0_y_x_yy_xx, g_y_0_0_0_y_x_yy_xy, g_y_0_0_0_y_x_yy_xz, g_y_0_0_0_y_x_yy_yy, g_y_0_0_0_y_x_yy_yz, g_y_0_0_0_y_x_yy_zz, g_yy_x_yy_xx, g_yy_x_yy_xy, g_yy_x_yy_xz, g_yy_x_yy_yy, g_yy_x_yy_yz, g_yy_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_x_yy_xx[i] = -g_0_x_yy_xx[i] + 2.0 * g_yy_x_yy_xx[i] * a_exp;

        g_y_0_0_0_y_x_yy_xy[i] = -g_0_x_yy_xy[i] + 2.0 * g_yy_x_yy_xy[i] * a_exp;

        g_y_0_0_0_y_x_yy_xz[i] = -g_0_x_yy_xz[i] + 2.0 * g_yy_x_yy_xz[i] * a_exp;

        g_y_0_0_0_y_x_yy_yy[i] = -g_0_x_yy_yy[i] + 2.0 * g_yy_x_yy_yy[i] * a_exp;

        g_y_0_0_0_y_x_yy_yz[i] = -g_0_x_yy_yz[i] + 2.0 * g_yy_x_yy_yz[i] * a_exp;

        g_y_0_0_0_y_x_yy_zz[i] = -g_0_x_yy_zz[i] + 2.0 * g_yy_x_yy_zz[i] * a_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_0_x_yz_xx, g_0_x_yz_xy, g_0_x_yz_xz, g_0_x_yz_yy, g_0_x_yz_yz, g_0_x_yz_zz, g_y_0_0_0_y_x_yz_xx, g_y_0_0_0_y_x_yz_xy, g_y_0_0_0_y_x_yz_xz, g_y_0_0_0_y_x_yz_yy, g_y_0_0_0_y_x_yz_yz, g_y_0_0_0_y_x_yz_zz, g_yy_x_yz_xx, g_yy_x_yz_xy, g_yy_x_yz_xz, g_yy_x_yz_yy, g_yy_x_yz_yz, g_yy_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_x_yz_xx[i] = -g_0_x_yz_xx[i] + 2.0 * g_yy_x_yz_xx[i] * a_exp;

        g_y_0_0_0_y_x_yz_xy[i] = -g_0_x_yz_xy[i] + 2.0 * g_yy_x_yz_xy[i] * a_exp;

        g_y_0_0_0_y_x_yz_xz[i] = -g_0_x_yz_xz[i] + 2.0 * g_yy_x_yz_xz[i] * a_exp;

        g_y_0_0_0_y_x_yz_yy[i] = -g_0_x_yz_yy[i] + 2.0 * g_yy_x_yz_yy[i] * a_exp;

        g_y_0_0_0_y_x_yz_yz[i] = -g_0_x_yz_yz[i] + 2.0 * g_yy_x_yz_yz[i] * a_exp;

        g_y_0_0_0_y_x_yz_zz[i] = -g_0_x_yz_zz[i] + 2.0 * g_yy_x_yz_zz[i] * a_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_0_x_zz_xx, g_0_x_zz_xy, g_0_x_zz_xz, g_0_x_zz_yy, g_0_x_zz_yz, g_0_x_zz_zz, g_y_0_0_0_y_x_zz_xx, g_y_0_0_0_y_x_zz_xy, g_y_0_0_0_y_x_zz_xz, g_y_0_0_0_y_x_zz_yy, g_y_0_0_0_y_x_zz_yz, g_y_0_0_0_y_x_zz_zz, g_yy_x_zz_xx, g_yy_x_zz_xy, g_yy_x_zz_xz, g_yy_x_zz_yy, g_yy_x_zz_yz, g_yy_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_x_zz_xx[i] = -g_0_x_zz_xx[i] + 2.0 * g_yy_x_zz_xx[i] * a_exp;

        g_y_0_0_0_y_x_zz_xy[i] = -g_0_x_zz_xy[i] + 2.0 * g_yy_x_zz_xy[i] * a_exp;

        g_y_0_0_0_y_x_zz_xz[i] = -g_0_x_zz_xz[i] + 2.0 * g_yy_x_zz_xz[i] * a_exp;

        g_y_0_0_0_y_x_zz_yy[i] = -g_0_x_zz_yy[i] + 2.0 * g_yy_x_zz_yy[i] * a_exp;

        g_y_0_0_0_y_x_zz_yz[i] = -g_0_x_zz_yz[i] + 2.0 * g_yy_x_zz_yz[i] * a_exp;

        g_y_0_0_0_y_x_zz_zz[i] = -g_0_x_zz_zz[i] + 2.0 * g_yy_x_zz_zz[i] * a_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_0_y_xx_xx, g_0_y_xx_xy, g_0_y_xx_xz, g_0_y_xx_yy, g_0_y_xx_yz, g_0_y_xx_zz, g_y_0_0_0_y_y_xx_xx, g_y_0_0_0_y_y_xx_xy, g_y_0_0_0_y_y_xx_xz, g_y_0_0_0_y_y_xx_yy, g_y_0_0_0_y_y_xx_yz, g_y_0_0_0_y_y_xx_zz, g_yy_y_xx_xx, g_yy_y_xx_xy, g_yy_y_xx_xz, g_yy_y_xx_yy, g_yy_y_xx_yz, g_yy_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_y_xx_xx[i] = -g_0_y_xx_xx[i] + 2.0 * g_yy_y_xx_xx[i] * a_exp;

        g_y_0_0_0_y_y_xx_xy[i] = -g_0_y_xx_xy[i] + 2.0 * g_yy_y_xx_xy[i] * a_exp;

        g_y_0_0_0_y_y_xx_xz[i] = -g_0_y_xx_xz[i] + 2.0 * g_yy_y_xx_xz[i] * a_exp;

        g_y_0_0_0_y_y_xx_yy[i] = -g_0_y_xx_yy[i] + 2.0 * g_yy_y_xx_yy[i] * a_exp;

        g_y_0_0_0_y_y_xx_yz[i] = -g_0_y_xx_yz[i] + 2.0 * g_yy_y_xx_yz[i] * a_exp;

        g_y_0_0_0_y_y_xx_zz[i] = -g_0_y_xx_zz[i] + 2.0 * g_yy_y_xx_zz[i] * a_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_0_y_xy_xx, g_0_y_xy_xy, g_0_y_xy_xz, g_0_y_xy_yy, g_0_y_xy_yz, g_0_y_xy_zz, g_y_0_0_0_y_y_xy_xx, g_y_0_0_0_y_y_xy_xy, g_y_0_0_0_y_y_xy_xz, g_y_0_0_0_y_y_xy_yy, g_y_0_0_0_y_y_xy_yz, g_y_0_0_0_y_y_xy_zz, g_yy_y_xy_xx, g_yy_y_xy_xy, g_yy_y_xy_xz, g_yy_y_xy_yy, g_yy_y_xy_yz, g_yy_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_y_xy_xx[i] = -g_0_y_xy_xx[i] + 2.0 * g_yy_y_xy_xx[i] * a_exp;

        g_y_0_0_0_y_y_xy_xy[i] = -g_0_y_xy_xy[i] + 2.0 * g_yy_y_xy_xy[i] * a_exp;

        g_y_0_0_0_y_y_xy_xz[i] = -g_0_y_xy_xz[i] + 2.0 * g_yy_y_xy_xz[i] * a_exp;

        g_y_0_0_0_y_y_xy_yy[i] = -g_0_y_xy_yy[i] + 2.0 * g_yy_y_xy_yy[i] * a_exp;

        g_y_0_0_0_y_y_xy_yz[i] = -g_0_y_xy_yz[i] + 2.0 * g_yy_y_xy_yz[i] * a_exp;

        g_y_0_0_0_y_y_xy_zz[i] = -g_0_y_xy_zz[i] + 2.0 * g_yy_y_xy_zz[i] * a_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_0_y_xz_xx, g_0_y_xz_xy, g_0_y_xz_xz, g_0_y_xz_yy, g_0_y_xz_yz, g_0_y_xz_zz, g_y_0_0_0_y_y_xz_xx, g_y_0_0_0_y_y_xz_xy, g_y_0_0_0_y_y_xz_xz, g_y_0_0_0_y_y_xz_yy, g_y_0_0_0_y_y_xz_yz, g_y_0_0_0_y_y_xz_zz, g_yy_y_xz_xx, g_yy_y_xz_xy, g_yy_y_xz_xz, g_yy_y_xz_yy, g_yy_y_xz_yz, g_yy_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_y_xz_xx[i] = -g_0_y_xz_xx[i] + 2.0 * g_yy_y_xz_xx[i] * a_exp;

        g_y_0_0_0_y_y_xz_xy[i] = -g_0_y_xz_xy[i] + 2.0 * g_yy_y_xz_xy[i] * a_exp;

        g_y_0_0_0_y_y_xz_xz[i] = -g_0_y_xz_xz[i] + 2.0 * g_yy_y_xz_xz[i] * a_exp;

        g_y_0_0_0_y_y_xz_yy[i] = -g_0_y_xz_yy[i] + 2.0 * g_yy_y_xz_yy[i] * a_exp;

        g_y_0_0_0_y_y_xz_yz[i] = -g_0_y_xz_yz[i] + 2.0 * g_yy_y_xz_yz[i] * a_exp;

        g_y_0_0_0_y_y_xz_zz[i] = -g_0_y_xz_zz[i] + 2.0 * g_yy_y_xz_zz[i] * a_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_0_y_yy_xx, g_0_y_yy_xy, g_0_y_yy_xz, g_0_y_yy_yy, g_0_y_yy_yz, g_0_y_yy_zz, g_y_0_0_0_y_y_yy_xx, g_y_0_0_0_y_y_yy_xy, g_y_0_0_0_y_y_yy_xz, g_y_0_0_0_y_y_yy_yy, g_y_0_0_0_y_y_yy_yz, g_y_0_0_0_y_y_yy_zz, g_yy_y_yy_xx, g_yy_y_yy_xy, g_yy_y_yy_xz, g_yy_y_yy_yy, g_yy_y_yy_yz, g_yy_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_y_yy_xx[i] = -g_0_y_yy_xx[i] + 2.0 * g_yy_y_yy_xx[i] * a_exp;

        g_y_0_0_0_y_y_yy_xy[i] = -g_0_y_yy_xy[i] + 2.0 * g_yy_y_yy_xy[i] * a_exp;

        g_y_0_0_0_y_y_yy_xz[i] = -g_0_y_yy_xz[i] + 2.0 * g_yy_y_yy_xz[i] * a_exp;

        g_y_0_0_0_y_y_yy_yy[i] = -g_0_y_yy_yy[i] + 2.0 * g_yy_y_yy_yy[i] * a_exp;

        g_y_0_0_0_y_y_yy_yz[i] = -g_0_y_yy_yz[i] + 2.0 * g_yy_y_yy_yz[i] * a_exp;

        g_y_0_0_0_y_y_yy_zz[i] = -g_0_y_yy_zz[i] + 2.0 * g_yy_y_yy_zz[i] * a_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_0_y_yz_xx, g_0_y_yz_xy, g_0_y_yz_xz, g_0_y_yz_yy, g_0_y_yz_yz, g_0_y_yz_zz, g_y_0_0_0_y_y_yz_xx, g_y_0_0_0_y_y_yz_xy, g_y_0_0_0_y_y_yz_xz, g_y_0_0_0_y_y_yz_yy, g_y_0_0_0_y_y_yz_yz, g_y_0_0_0_y_y_yz_zz, g_yy_y_yz_xx, g_yy_y_yz_xy, g_yy_y_yz_xz, g_yy_y_yz_yy, g_yy_y_yz_yz, g_yy_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_y_yz_xx[i] = -g_0_y_yz_xx[i] + 2.0 * g_yy_y_yz_xx[i] * a_exp;

        g_y_0_0_0_y_y_yz_xy[i] = -g_0_y_yz_xy[i] + 2.0 * g_yy_y_yz_xy[i] * a_exp;

        g_y_0_0_0_y_y_yz_xz[i] = -g_0_y_yz_xz[i] + 2.0 * g_yy_y_yz_xz[i] * a_exp;

        g_y_0_0_0_y_y_yz_yy[i] = -g_0_y_yz_yy[i] + 2.0 * g_yy_y_yz_yy[i] * a_exp;

        g_y_0_0_0_y_y_yz_yz[i] = -g_0_y_yz_yz[i] + 2.0 * g_yy_y_yz_yz[i] * a_exp;

        g_y_0_0_0_y_y_yz_zz[i] = -g_0_y_yz_zz[i] + 2.0 * g_yy_y_yz_zz[i] * a_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_0_y_zz_xx, g_0_y_zz_xy, g_0_y_zz_xz, g_0_y_zz_yy, g_0_y_zz_yz, g_0_y_zz_zz, g_y_0_0_0_y_y_zz_xx, g_y_0_0_0_y_y_zz_xy, g_y_0_0_0_y_y_zz_xz, g_y_0_0_0_y_y_zz_yy, g_y_0_0_0_y_y_zz_yz, g_y_0_0_0_y_y_zz_zz, g_yy_y_zz_xx, g_yy_y_zz_xy, g_yy_y_zz_xz, g_yy_y_zz_yy, g_yy_y_zz_yz, g_yy_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_y_zz_xx[i] = -g_0_y_zz_xx[i] + 2.0 * g_yy_y_zz_xx[i] * a_exp;

        g_y_0_0_0_y_y_zz_xy[i] = -g_0_y_zz_xy[i] + 2.0 * g_yy_y_zz_xy[i] * a_exp;

        g_y_0_0_0_y_y_zz_xz[i] = -g_0_y_zz_xz[i] + 2.0 * g_yy_y_zz_xz[i] * a_exp;

        g_y_0_0_0_y_y_zz_yy[i] = -g_0_y_zz_yy[i] + 2.0 * g_yy_y_zz_yy[i] * a_exp;

        g_y_0_0_0_y_y_zz_yz[i] = -g_0_y_zz_yz[i] + 2.0 * g_yy_y_zz_yz[i] * a_exp;

        g_y_0_0_0_y_y_zz_zz[i] = -g_0_y_zz_zz[i] + 2.0 * g_yy_y_zz_zz[i] * a_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_0_z_xx_xx, g_0_z_xx_xy, g_0_z_xx_xz, g_0_z_xx_yy, g_0_z_xx_yz, g_0_z_xx_zz, g_y_0_0_0_y_z_xx_xx, g_y_0_0_0_y_z_xx_xy, g_y_0_0_0_y_z_xx_xz, g_y_0_0_0_y_z_xx_yy, g_y_0_0_0_y_z_xx_yz, g_y_0_0_0_y_z_xx_zz, g_yy_z_xx_xx, g_yy_z_xx_xy, g_yy_z_xx_xz, g_yy_z_xx_yy, g_yy_z_xx_yz, g_yy_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_z_xx_xx[i] = -g_0_z_xx_xx[i] + 2.0 * g_yy_z_xx_xx[i] * a_exp;

        g_y_0_0_0_y_z_xx_xy[i] = -g_0_z_xx_xy[i] + 2.0 * g_yy_z_xx_xy[i] * a_exp;

        g_y_0_0_0_y_z_xx_xz[i] = -g_0_z_xx_xz[i] + 2.0 * g_yy_z_xx_xz[i] * a_exp;

        g_y_0_0_0_y_z_xx_yy[i] = -g_0_z_xx_yy[i] + 2.0 * g_yy_z_xx_yy[i] * a_exp;

        g_y_0_0_0_y_z_xx_yz[i] = -g_0_z_xx_yz[i] + 2.0 * g_yy_z_xx_yz[i] * a_exp;

        g_y_0_0_0_y_z_xx_zz[i] = -g_0_z_xx_zz[i] + 2.0 * g_yy_z_xx_zz[i] * a_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_0_z_xy_xx, g_0_z_xy_xy, g_0_z_xy_xz, g_0_z_xy_yy, g_0_z_xy_yz, g_0_z_xy_zz, g_y_0_0_0_y_z_xy_xx, g_y_0_0_0_y_z_xy_xy, g_y_0_0_0_y_z_xy_xz, g_y_0_0_0_y_z_xy_yy, g_y_0_0_0_y_z_xy_yz, g_y_0_0_0_y_z_xy_zz, g_yy_z_xy_xx, g_yy_z_xy_xy, g_yy_z_xy_xz, g_yy_z_xy_yy, g_yy_z_xy_yz, g_yy_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_z_xy_xx[i] = -g_0_z_xy_xx[i] + 2.0 * g_yy_z_xy_xx[i] * a_exp;

        g_y_0_0_0_y_z_xy_xy[i] = -g_0_z_xy_xy[i] + 2.0 * g_yy_z_xy_xy[i] * a_exp;

        g_y_0_0_0_y_z_xy_xz[i] = -g_0_z_xy_xz[i] + 2.0 * g_yy_z_xy_xz[i] * a_exp;

        g_y_0_0_0_y_z_xy_yy[i] = -g_0_z_xy_yy[i] + 2.0 * g_yy_z_xy_yy[i] * a_exp;

        g_y_0_0_0_y_z_xy_yz[i] = -g_0_z_xy_yz[i] + 2.0 * g_yy_z_xy_yz[i] * a_exp;

        g_y_0_0_0_y_z_xy_zz[i] = -g_0_z_xy_zz[i] + 2.0 * g_yy_z_xy_zz[i] * a_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_0_z_xz_xx, g_0_z_xz_xy, g_0_z_xz_xz, g_0_z_xz_yy, g_0_z_xz_yz, g_0_z_xz_zz, g_y_0_0_0_y_z_xz_xx, g_y_0_0_0_y_z_xz_xy, g_y_0_0_0_y_z_xz_xz, g_y_0_0_0_y_z_xz_yy, g_y_0_0_0_y_z_xz_yz, g_y_0_0_0_y_z_xz_zz, g_yy_z_xz_xx, g_yy_z_xz_xy, g_yy_z_xz_xz, g_yy_z_xz_yy, g_yy_z_xz_yz, g_yy_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_z_xz_xx[i] = -g_0_z_xz_xx[i] + 2.0 * g_yy_z_xz_xx[i] * a_exp;

        g_y_0_0_0_y_z_xz_xy[i] = -g_0_z_xz_xy[i] + 2.0 * g_yy_z_xz_xy[i] * a_exp;

        g_y_0_0_0_y_z_xz_xz[i] = -g_0_z_xz_xz[i] + 2.0 * g_yy_z_xz_xz[i] * a_exp;

        g_y_0_0_0_y_z_xz_yy[i] = -g_0_z_xz_yy[i] + 2.0 * g_yy_z_xz_yy[i] * a_exp;

        g_y_0_0_0_y_z_xz_yz[i] = -g_0_z_xz_yz[i] + 2.0 * g_yy_z_xz_yz[i] * a_exp;

        g_y_0_0_0_y_z_xz_zz[i] = -g_0_z_xz_zz[i] + 2.0 * g_yy_z_xz_zz[i] * a_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_0_z_yy_xx, g_0_z_yy_xy, g_0_z_yy_xz, g_0_z_yy_yy, g_0_z_yy_yz, g_0_z_yy_zz, g_y_0_0_0_y_z_yy_xx, g_y_0_0_0_y_z_yy_xy, g_y_0_0_0_y_z_yy_xz, g_y_0_0_0_y_z_yy_yy, g_y_0_0_0_y_z_yy_yz, g_y_0_0_0_y_z_yy_zz, g_yy_z_yy_xx, g_yy_z_yy_xy, g_yy_z_yy_xz, g_yy_z_yy_yy, g_yy_z_yy_yz, g_yy_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_z_yy_xx[i] = -g_0_z_yy_xx[i] + 2.0 * g_yy_z_yy_xx[i] * a_exp;

        g_y_0_0_0_y_z_yy_xy[i] = -g_0_z_yy_xy[i] + 2.0 * g_yy_z_yy_xy[i] * a_exp;

        g_y_0_0_0_y_z_yy_xz[i] = -g_0_z_yy_xz[i] + 2.0 * g_yy_z_yy_xz[i] * a_exp;

        g_y_0_0_0_y_z_yy_yy[i] = -g_0_z_yy_yy[i] + 2.0 * g_yy_z_yy_yy[i] * a_exp;

        g_y_0_0_0_y_z_yy_yz[i] = -g_0_z_yy_yz[i] + 2.0 * g_yy_z_yy_yz[i] * a_exp;

        g_y_0_0_0_y_z_yy_zz[i] = -g_0_z_yy_zz[i] + 2.0 * g_yy_z_yy_zz[i] * a_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_0_z_yz_xx, g_0_z_yz_xy, g_0_z_yz_xz, g_0_z_yz_yy, g_0_z_yz_yz, g_0_z_yz_zz, g_y_0_0_0_y_z_yz_xx, g_y_0_0_0_y_z_yz_xy, g_y_0_0_0_y_z_yz_xz, g_y_0_0_0_y_z_yz_yy, g_y_0_0_0_y_z_yz_yz, g_y_0_0_0_y_z_yz_zz, g_yy_z_yz_xx, g_yy_z_yz_xy, g_yy_z_yz_xz, g_yy_z_yz_yy, g_yy_z_yz_yz, g_yy_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_z_yz_xx[i] = -g_0_z_yz_xx[i] + 2.0 * g_yy_z_yz_xx[i] * a_exp;

        g_y_0_0_0_y_z_yz_xy[i] = -g_0_z_yz_xy[i] + 2.0 * g_yy_z_yz_xy[i] * a_exp;

        g_y_0_0_0_y_z_yz_xz[i] = -g_0_z_yz_xz[i] + 2.0 * g_yy_z_yz_xz[i] * a_exp;

        g_y_0_0_0_y_z_yz_yy[i] = -g_0_z_yz_yy[i] + 2.0 * g_yy_z_yz_yy[i] * a_exp;

        g_y_0_0_0_y_z_yz_yz[i] = -g_0_z_yz_yz[i] + 2.0 * g_yy_z_yz_yz[i] * a_exp;

        g_y_0_0_0_y_z_yz_zz[i] = -g_0_z_yz_zz[i] + 2.0 * g_yy_z_yz_zz[i] * a_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_0_z_zz_xx, g_0_z_zz_xy, g_0_z_zz_xz, g_0_z_zz_yy, g_0_z_zz_yz, g_0_z_zz_zz, g_y_0_0_0_y_z_zz_xx, g_y_0_0_0_y_z_zz_xy, g_y_0_0_0_y_z_zz_xz, g_y_0_0_0_y_z_zz_yy, g_y_0_0_0_y_z_zz_yz, g_y_0_0_0_y_z_zz_zz, g_yy_z_zz_xx, g_yy_z_zz_xy, g_yy_z_zz_xz, g_yy_z_zz_yy, g_yy_z_zz_yz, g_yy_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_z_zz_xx[i] = -g_0_z_zz_xx[i] + 2.0 * g_yy_z_zz_xx[i] * a_exp;

        g_y_0_0_0_y_z_zz_xy[i] = -g_0_z_zz_xy[i] + 2.0 * g_yy_z_zz_xy[i] * a_exp;

        g_y_0_0_0_y_z_zz_xz[i] = -g_0_z_zz_xz[i] + 2.0 * g_yy_z_zz_xz[i] * a_exp;

        g_y_0_0_0_y_z_zz_yy[i] = -g_0_z_zz_yy[i] + 2.0 * g_yy_z_zz_yy[i] * a_exp;

        g_y_0_0_0_y_z_zz_yz[i] = -g_0_z_zz_yz[i] + 2.0 * g_yy_z_zz_yz[i] * a_exp;

        g_y_0_0_0_y_z_zz_zz[i] = -g_0_z_zz_zz[i] + 2.0 * g_yy_z_zz_zz[i] * a_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_y_0_0_0_z_x_xx_xx, g_y_0_0_0_z_x_xx_xy, g_y_0_0_0_z_x_xx_xz, g_y_0_0_0_z_x_xx_yy, g_y_0_0_0_z_x_xx_yz, g_y_0_0_0_z_x_xx_zz, g_yz_x_xx_xx, g_yz_x_xx_xy, g_yz_x_xx_xz, g_yz_x_xx_yy, g_yz_x_xx_yz, g_yz_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_x_xx_xx[i] = 2.0 * g_yz_x_xx_xx[i] * a_exp;

        g_y_0_0_0_z_x_xx_xy[i] = 2.0 * g_yz_x_xx_xy[i] * a_exp;

        g_y_0_0_0_z_x_xx_xz[i] = 2.0 * g_yz_x_xx_xz[i] * a_exp;

        g_y_0_0_0_z_x_xx_yy[i] = 2.0 * g_yz_x_xx_yy[i] * a_exp;

        g_y_0_0_0_z_x_xx_yz[i] = 2.0 * g_yz_x_xx_yz[i] * a_exp;

        g_y_0_0_0_z_x_xx_zz[i] = 2.0 * g_yz_x_xx_zz[i] * a_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_y_0_0_0_z_x_xy_xx, g_y_0_0_0_z_x_xy_xy, g_y_0_0_0_z_x_xy_xz, g_y_0_0_0_z_x_xy_yy, g_y_0_0_0_z_x_xy_yz, g_y_0_0_0_z_x_xy_zz, g_yz_x_xy_xx, g_yz_x_xy_xy, g_yz_x_xy_xz, g_yz_x_xy_yy, g_yz_x_xy_yz, g_yz_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_x_xy_xx[i] = 2.0 * g_yz_x_xy_xx[i] * a_exp;

        g_y_0_0_0_z_x_xy_xy[i] = 2.0 * g_yz_x_xy_xy[i] * a_exp;

        g_y_0_0_0_z_x_xy_xz[i] = 2.0 * g_yz_x_xy_xz[i] * a_exp;

        g_y_0_0_0_z_x_xy_yy[i] = 2.0 * g_yz_x_xy_yy[i] * a_exp;

        g_y_0_0_0_z_x_xy_yz[i] = 2.0 * g_yz_x_xy_yz[i] * a_exp;

        g_y_0_0_0_z_x_xy_zz[i] = 2.0 * g_yz_x_xy_zz[i] * a_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_y_0_0_0_z_x_xz_xx, g_y_0_0_0_z_x_xz_xy, g_y_0_0_0_z_x_xz_xz, g_y_0_0_0_z_x_xz_yy, g_y_0_0_0_z_x_xz_yz, g_y_0_0_0_z_x_xz_zz, g_yz_x_xz_xx, g_yz_x_xz_xy, g_yz_x_xz_xz, g_yz_x_xz_yy, g_yz_x_xz_yz, g_yz_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_x_xz_xx[i] = 2.0 * g_yz_x_xz_xx[i] * a_exp;

        g_y_0_0_0_z_x_xz_xy[i] = 2.0 * g_yz_x_xz_xy[i] * a_exp;

        g_y_0_0_0_z_x_xz_xz[i] = 2.0 * g_yz_x_xz_xz[i] * a_exp;

        g_y_0_0_0_z_x_xz_yy[i] = 2.0 * g_yz_x_xz_yy[i] * a_exp;

        g_y_0_0_0_z_x_xz_yz[i] = 2.0 * g_yz_x_xz_yz[i] * a_exp;

        g_y_0_0_0_z_x_xz_zz[i] = 2.0 * g_yz_x_xz_zz[i] * a_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_y_0_0_0_z_x_yy_xx, g_y_0_0_0_z_x_yy_xy, g_y_0_0_0_z_x_yy_xz, g_y_0_0_0_z_x_yy_yy, g_y_0_0_0_z_x_yy_yz, g_y_0_0_0_z_x_yy_zz, g_yz_x_yy_xx, g_yz_x_yy_xy, g_yz_x_yy_xz, g_yz_x_yy_yy, g_yz_x_yy_yz, g_yz_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_x_yy_xx[i] = 2.0 * g_yz_x_yy_xx[i] * a_exp;

        g_y_0_0_0_z_x_yy_xy[i] = 2.0 * g_yz_x_yy_xy[i] * a_exp;

        g_y_0_0_0_z_x_yy_xz[i] = 2.0 * g_yz_x_yy_xz[i] * a_exp;

        g_y_0_0_0_z_x_yy_yy[i] = 2.0 * g_yz_x_yy_yy[i] * a_exp;

        g_y_0_0_0_z_x_yy_yz[i] = 2.0 * g_yz_x_yy_yz[i] * a_exp;

        g_y_0_0_0_z_x_yy_zz[i] = 2.0 * g_yz_x_yy_zz[i] * a_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_y_0_0_0_z_x_yz_xx, g_y_0_0_0_z_x_yz_xy, g_y_0_0_0_z_x_yz_xz, g_y_0_0_0_z_x_yz_yy, g_y_0_0_0_z_x_yz_yz, g_y_0_0_0_z_x_yz_zz, g_yz_x_yz_xx, g_yz_x_yz_xy, g_yz_x_yz_xz, g_yz_x_yz_yy, g_yz_x_yz_yz, g_yz_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_x_yz_xx[i] = 2.0 * g_yz_x_yz_xx[i] * a_exp;

        g_y_0_0_0_z_x_yz_xy[i] = 2.0 * g_yz_x_yz_xy[i] * a_exp;

        g_y_0_0_0_z_x_yz_xz[i] = 2.0 * g_yz_x_yz_xz[i] * a_exp;

        g_y_0_0_0_z_x_yz_yy[i] = 2.0 * g_yz_x_yz_yy[i] * a_exp;

        g_y_0_0_0_z_x_yz_yz[i] = 2.0 * g_yz_x_yz_yz[i] * a_exp;

        g_y_0_0_0_z_x_yz_zz[i] = 2.0 * g_yz_x_yz_zz[i] * a_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_y_0_0_0_z_x_zz_xx, g_y_0_0_0_z_x_zz_xy, g_y_0_0_0_z_x_zz_xz, g_y_0_0_0_z_x_zz_yy, g_y_0_0_0_z_x_zz_yz, g_y_0_0_0_z_x_zz_zz, g_yz_x_zz_xx, g_yz_x_zz_xy, g_yz_x_zz_xz, g_yz_x_zz_yy, g_yz_x_zz_yz, g_yz_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_x_zz_xx[i] = 2.0 * g_yz_x_zz_xx[i] * a_exp;

        g_y_0_0_0_z_x_zz_xy[i] = 2.0 * g_yz_x_zz_xy[i] * a_exp;

        g_y_0_0_0_z_x_zz_xz[i] = 2.0 * g_yz_x_zz_xz[i] * a_exp;

        g_y_0_0_0_z_x_zz_yy[i] = 2.0 * g_yz_x_zz_yy[i] * a_exp;

        g_y_0_0_0_z_x_zz_yz[i] = 2.0 * g_yz_x_zz_yz[i] * a_exp;

        g_y_0_0_0_z_x_zz_zz[i] = 2.0 * g_yz_x_zz_zz[i] * a_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_y_0_0_0_z_y_xx_xx, g_y_0_0_0_z_y_xx_xy, g_y_0_0_0_z_y_xx_xz, g_y_0_0_0_z_y_xx_yy, g_y_0_0_0_z_y_xx_yz, g_y_0_0_0_z_y_xx_zz, g_yz_y_xx_xx, g_yz_y_xx_xy, g_yz_y_xx_xz, g_yz_y_xx_yy, g_yz_y_xx_yz, g_yz_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_y_xx_xx[i] = 2.0 * g_yz_y_xx_xx[i] * a_exp;

        g_y_0_0_0_z_y_xx_xy[i] = 2.0 * g_yz_y_xx_xy[i] * a_exp;

        g_y_0_0_0_z_y_xx_xz[i] = 2.0 * g_yz_y_xx_xz[i] * a_exp;

        g_y_0_0_0_z_y_xx_yy[i] = 2.0 * g_yz_y_xx_yy[i] * a_exp;

        g_y_0_0_0_z_y_xx_yz[i] = 2.0 * g_yz_y_xx_yz[i] * a_exp;

        g_y_0_0_0_z_y_xx_zz[i] = 2.0 * g_yz_y_xx_zz[i] * a_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_y_0_0_0_z_y_xy_xx, g_y_0_0_0_z_y_xy_xy, g_y_0_0_0_z_y_xy_xz, g_y_0_0_0_z_y_xy_yy, g_y_0_0_0_z_y_xy_yz, g_y_0_0_0_z_y_xy_zz, g_yz_y_xy_xx, g_yz_y_xy_xy, g_yz_y_xy_xz, g_yz_y_xy_yy, g_yz_y_xy_yz, g_yz_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_y_xy_xx[i] = 2.0 * g_yz_y_xy_xx[i] * a_exp;

        g_y_0_0_0_z_y_xy_xy[i] = 2.0 * g_yz_y_xy_xy[i] * a_exp;

        g_y_0_0_0_z_y_xy_xz[i] = 2.0 * g_yz_y_xy_xz[i] * a_exp;

        g_y_0_0_0_z_y_xy_yy[i] = 2.0 * g_yz_y_xy_yy[i] * a_exp;

        g_y_0_0_0_z_y_xy_yz[i] = 2.0 * g_yz_y_xy_yz[i] * a_exp;

        g_y_0_0_0_z_y_xy_zz[i] = 2.0 * g_yz_y_xy_zz[i] * a_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_y_0_0_0_z_y_xz_xx, g_y_0_0_0_z_y_xz_xy, g_y_0_0_0_z_y_xz_xz, g_y_0_0_0_z_y_xz_yy, g_y_0_0_0_z_y_xz_yz, g_y_0_0_0_z_y_xz_zz, g_yz_y_xz_xx, g_yz_y_xz_xy, g_yz_y_xz_xz, g_yz_y_xz_yy, g_yz_y_xz_yz, g_yz_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_y_xz_xx[i] = 2.0 * g_yz_y_xz_xx[i] * a_exp;

        g_y_0_0_0_z_y_xz_xy[i] = 2.0 * g_yz_y_xz_xy[i] * a_exp;

        g_y_0_0_0_z_y_xz_xz[i] = 2.0 * g_yz_y_xz_xz[i] * a_exp;

        g_y_0_0_0_z_y_xz_yy[i] = 2.0 * g_yz_y_xz_yy[i] * a_exp;

        g_y_0_0_0_z_y_xz_yz[i] = 2.0 * g_yz_y_xz_yz[i] * a_exp;

        g_y_0_0_0_z_y_xz_zz[i] = 2.0 * g_yz_y_xz_zz[i] * a_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_y_0_0_0_z_y_yy_xx, g_y_0_0_0_z_y_yy_xy, g_y_0_0_0_z_y_yy_xz, g_y_0_0_0_z_y_yy_yy, g_y_0_0_0_z_y_yy_yz, g_y_0_0_0_z_y_yy_zz, g_yz_y_yy_xx, g_yz_y_yy_xy, g_yz_y_yy_xz, g_yz_y_yy_yy, g_yz_y_yy_yz, g_yz_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_y_yy_xx[i] = 2.0 * g_yz_y_yy_xx[i] * a_exp;

        g_y_0_0_0_z_y_yy_xy[i] = 2.0 * g_yz_y_yy_xy[i] * a_exp;

        g_y_0_0_0_z_y_yy_xz[i] = 2.0 * g_yz_y_yy_xz[i] * a_exp;

        g_y_0_0_0_z_y_yy_yy[i] = 2.0 * g_yz_y_yy_yy[i] * a_exp;

        g_y_0_0_0_z_y_yy_yz[i] = 2.0 * g_yz_y_yy_yz[i] * a_exp;

        g_y_0_0_0_z_y_yy_zz[i] = 2.0 * g_yz_y_yy_zz[i] * a_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_y_0_0_0_z_y_yz_xx, g_y_0_0_0_z_y_yz_xy, g_y_0_0_0_z_y_yz_xz, g_y_0_0_0_z_y_yz_yy, g_y_0_0_0_z_y_yz_yz, g_y_0_0_0_z_y_yz_zz, g_yz_y_yz_xx, g_yz_y_yz_xy, g_yz_y_yz_xz, g_yz_y_yz_yy, g_yz_y_yz_yz, g_yz_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_y_yz_xx[i] = 2.0 * g_yz_y_yz_xx[i] * a_exp;

        g_y_0_0_0_z_y_yz_xy[i] = 2.0 * g_yz_y_yz_xy[i] * a_exp;

        g_y_0_0_0_z_y_yz_xz[i] = 2.0 * g_yz_y_yz_xz[i] * a_exp;

        g_y_0_0_0_z_y_yz_yy[i] = 2.0 * g_yz_y_yz_yy[i] * a_exp;

        g_y_0_0_0_z_y_yz_yz[i] = 2.0 * g_yz_y_yz_yz[i] * a_exp;

        g_y_0_0_0_z_y_yz_zz[i] = 2.0 * g_yz_y_yz_zz[i] * a_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_y_0_0_0_z_y_zz_xx, g_y_0_0_0_z_y_zz_xy, g_y_0_0_0_z_y_zz_xz, g_y_0_0_0_z_y_zz_yy, g_y_0_0_0_z_y_zz_yz, g_y_0_0_0_z_y_zz_zz, g_yz_y_zz_xx, g_yz_y_zz_xy, g_yz_y_zz_xz, g_yz_y_zz_yy, g_yz_y_zz_yz, g_yz_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_y_zz_xx[i] = 2.0 * g_yz_y_zz_xx[i] * a_exp;

        g_y_0_0_0_z_y_zz_xy[i] = 2.0 * g_yz_y_zz_xy[i] * a_exp;

        g_y_0_0_0_z_y_zz_xz[i] = 2.0 * g_yz_y_zz_xz[i] * a_exp;

        g_y_0_0_0_z_y_zz_yy[i] = 2.0 * g_yz_y_zz_yy[i] * a_exp;

        g_y_0_0_0_z_y_zz_yz[i] = 2.0 * g_yz_y_zz_yz[i] * a_exp;

        g_y_0_0_0_z_y_zz_zz[i] = 2.0 * g_yz_y_zz_zz[i] * a_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_y_0_0_0_z_z_xx_xx, g_y_0_0_0_z_z_xx_xy, g_y_0_0_0_z_z_xx_xz, g_y_0_0_0_z_z_xx_yy, g_y_0_0_0_z_z_xx_yz, g_y_0_0_0_z_z_xx_zz, g_yz_z_xx_xx, g_yz_z_xx_xy, g_yz_z_xx_xz, g_yz_z_xx_yy, g_yz_z_xx_yz, g_yz_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_z_xx_xx[i] = 2.0 * g_yz_z_xx_xx[i] * a_exp;

        g_y_0_0_0_z_z_xx_xy[i] = 2.0 * g_yz_z_xx_xy[i] * a_exp;

        g_y_0_0_0_z_z_xx_xz[i] = 2.0 * g_yz_z_xx_xz[i] * a_exp;

        g_y_0_0_0_z_z_xx_yy[i] = 2.0 * g_yz_z_xx_yy[i] * a_exp;

        g_y_0_0_0_z_z_xx_yz[i] = 2.0 * g_yz_z_xx_yz[i] * a_exp;

        g_y_0_0_0_z_z_xx_zz[i] = 2.0 * g_yz_z_xx_zz[i] * a_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_y_0_0_0_z_z_xy_xx, g_y_0_0_0_z_z_xy_xy, g_y_0_0_0_z_z_xy_xz, g_y_0_0_0_z_z_xy_yy, g_y_0_0_0_z_z_xy_yz, g_y_0_0_0_z_z_xy_zz, g_yz_z_xy_xx, g_yz_z_xy_xy, g_yz_z_xy_xz, g_yz_z_xy_yy, g_yz_z_xy_yz, g_yz_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_z_xy_xx[i] = 2.0 * g_yz_z_xy_xx[i] * a_exp;

        g_y_0_0_0_z_z_xy_xy[i] = 2.0 * g_yz_z_xy_xy[i] * a_exp;

        g_y_0_0_0_z_z_xy_xz[i] = 2.0 * g_yz_z_xy_xz[i] * a_exp;

        g_y_0_0_0_z_z_xy_yy[i] = 2.0 * g_yz_z_xy_yy[i] * a_exp;

        g_y_0_0_0_z_z_xy_yz[i] = 2.0 * g_yz_z_xy_yz[i] * a_exp;

        g_y_0_0_0_z_z_xy_zz[i] = 2.0 * g_yz_z_xy_zz[i] * a_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_y_0_0_0_z_z_xz_xx, g_y_0_0_0_z_z_xz_xy, g_y_0_0_0_z_z_xz_xz, g_y_0_0_0_z_z_xz_yy, g_y_0_0_0_z_z_xz_yz, g_y_0_0_0_z_z_xz_zz, g_yz_z_xz_xx, g_yz_z_xz_xy, g_yz_z_xz_xz, g_yz_z_xz_yy, g_yz_z_xz_yz, g_yz_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_z_xz_xx[i] = 2.0 * g_yz_z_xz_xx[i] * a_exp;

        g_y_0_0_0_z_z_xz_xy[i] = 2.0 * g_yz_z_xz_xy[i] * a_exp;

        g_y_0_0_0_z_z_xz_xz[i] = 2.0 * g_yz_z_xz_xz[i] * a_exp;

        g_y_0_0_0_z_z_xz_yy[i] = 2.0 * g_yz_z_xz_yy[i] * a_exp;

        g_y_0_0_0_z_z_xz_yz[i] = 2.0 * g_yz_z_xz_yz[i] * a_exp;

        g_y_0_0_0_z_z_xz_zz[i] = 2.0 * g_yz_z_xz_zz[i] * a_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_y_0_0_0_z_z_yy_xx, g_y_0_0_0_z_z_yy_xy, g_y_0_0_0_z_z_yy_xz, g_y_0_0_0_z_z_yy_yy, g_y_0_0_0_z_z_yy_yz, g_y_0_0_0_z_z_yy_zz, g_yz_z_yy_xx, g_yz_z_yy_xy, g_yz_z_yy_xz, g_yz_z_yy_yy, g_yz_z_yy_yz, g_yz_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_z_yy_xx[i] = 2.0 * g_yz_z_yy_xx[i] * a_exp;

        g_y_0_0_0_z_z_yy_xy[i] = 2.0 * g_yz_z_yy_xy[i] * a_exp;

        g_y_0_0_0_z_z_yy_xz[i] = 2.0 * g_yz_z_yy_xz[i] * a_exp;

        g_y_0_0_0_z_z_yy_yy[i] = 2.0 * g_yz_z_yy_yy[i] * a_exp;

        g_y_0_0_0_z_z_yy_yz[i] = 2.0 * g_yz_z_yy_yz[i] * a_exp;

        g_y_0_0_0_z_z_yy_zz[i] = 2.0 * g_yz_z_yy_zz[i] * a_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_y_0_0_0_z_z_yz_xx, g_y_0_0_0_z_z_yz_xy, g_y_0_0_0_z_z_yz_xz, g_y_0_0_0_z_z_yz_yy, g_y_0_0_0_z_z_yz_yz, g_y_0_0_0_z_z_yz_zz, g_yz_z_yz_xx, g_yz_z_yz_xy, g_yz_z_yz_xz, g_yz_z_yz_yy, g_yz_z_yz_yz, g_yz_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_z_yz_xx[i] = 2.0 * g_yz_z_yz_xx[i] * a_exp;

        g_y_0_0_0_z_z_yz_xy[i] = 2.0 * g_yz_z_yz_xy[i] * a_exp;

        g_y_0_0_0_z_z_yz_xz[i] = 2.0 * g_yz_z_yz_xz[i] * a_exp;

        g_y_0_0_0_z_z_yz_yy[i] = 2.0 * g_yz_z_yz_yy[i] * a_exp;

        g_y_0_0_0_z_z_yz_yz[i] = 2.0 * g_yz_z_yz_yz[i] * a_exp;

        g_y_0_0_0_z_z_yz_zz[i] = 2.0 * g_yz_z_yz_zz[i] * a_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_y_0_0_0_z_z_zz_xx, g_y_0_0_0_z_z_zz_xy, g_y_0_0_0_z_z_zz_xz, g_y_0_0_0_z_z_zz_yy, g_y_0_0_0_z_z_zz_yz, g_y_0_0_0_z_z_zz_zz, g_yz_z_zz_xx, g_yz_z_zz_xy, g_yz_z_zz_xz, g_yz_z_zz_yy, g_yz_z_zz_yz, g_yz_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_z_zz_xx[i] = 2.0 * g_yz_z_zz_xx[i] * a_exp;

        g_y_0_0_0_z_z_zz_xy[i] = 2.0 * g_yz_z_zz_xy[i] * a_exp;

        g_y_0_0_0_z_z_zz_xz[i] = 2.0 * g_yz_z_zz_xz[i] * a_exp;

        g_y_0_0_0_z_z_zz_yy[i] = 2.0 * g_yz_z_zz_yy[i] * a_exp;

        g_y_0_0_0_z_z_zz_yz[i] = 2.0 * g_yz_z_zz_yz[i] * a_exp;

        g_y_0_0_0_z_z_zz_zz[i] = 2.0 * g_yz_z_zz_zz[i] * a_exp;
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_xz_x_xx_xx, g_xz_x_xx_xy, g_xz_x_xx_xz, g_xz_x_xx_yy, g_xz_x_xx_yz, g_xz_x_xx_zz, g_z_0_0_0_x_x_xx_xx, g_z_0_0_0_x_x_xx_xy, g_z_0_0_0_x_x_xx_xz, g_z_0_0_0_x_x_xx_yy, g_z_0_0_0_x_x_xx_yz, g_z_0_0_0_x_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_x_xx_xx[i] = 2.0 * g_xz_x_xx_xx[i] * a_exp;

        g_z_0_0_0_x_x_xx_xy[i] = 2.0 * g_xz_x_xx_xy[i] * a_exp;

        g_z_0_0_0_x_x_xx_xz[i] = 2.0 * g_xz_x_xx_xz[i] * a_exp;

        g_z_0_0_0_x_x_xx_yy[i] = 2.0 * g_xz_x_xx_yy[i] * a_exp;

        g_z_0_0_0_x_x_xx_yz[i] = 2.0 * g_xz_x_xx_yz[i] * a_exp;

        g_z_0_0_0_x_x_xx_zz[i] = 2.0 * g_xz_x_xx_zz[i] * a_exp;
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_xz_x_xy_xx, g_xz_x_xy_xy, g_xz_x_xy_xz, g_xz_x_xy_yy, g_xz_x_xy_yz, g_xz_x_xy_zz, g_z_0_0_0_x_x_xy_xx, g_z_0_0_0_x_x_xy_xy, g_z_0_0_0_x_x_xy_xz, g_z_0_0_0_x_x_xy_yy, g_z_0_0_0_x_x_xy_yz, g_z_0_0_0_x_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_x_xy_xx[i] = 2.0 * g_xz_x_xy_xx[i] * a_exp;

        g_z_0_0_0_x_x_xy_xy[i] = 2.0 * g_xz_x_xy_xy[i] * a_exp;

        g_z_0_0_0_x_x_xy_xz[i] = 2.0 * g_xz_x_xy_xz[i] * a_exp;

        g_z_0_0_0_x_x_xy_yy[i] = 2.0 * g_xz_x_xy_yy[i] * a_exp;

        g_z_0_0_0_x_x_xy_yz[i] = 2.0 * g_xz_x_xy_yz[i] * a_exp;

        g_z_0_0_0_x_x_xy_zz[i] = 2.0 * g_xz_x_xy_zz[i] * a_exp;
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_xz_x_xz_xx, g_xz_x_xz_xy, g_xz_x_xz_xz, g_xz_x_xz_yy, g_xz_x_xz_yz, g_xz_x_xz_zz, g_z_0_0_0_x_x_xz_xx, g_z_0_0_0_x_x_xz_xy, g_z_0_0_0_x_x_xz_xz, g_z_0_0_0_x_x_xz_yy, g_z_0_0_0_x_x_xz_yz, g_z_0_0_0_x_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_x_xz_xx[i] = 2.0 * g_xz_x_xz_xx[i] * a_exp;

        g_z_0_0_0_x_x_xz_xy[i] = 2.0 * g_xz_x_xz_xy[i] * a_exp;

        g_z_0_0_0_x_x_xz_xz[i] = 2.0 * g_xz_x_xz_xz[i] * a_exp;

        g_z_0_0_0_x_x_xz_yy[i] = 2.0 * g_xz_x_xz_yy[i] * a_exp;

        g_z_0_0_0_x_x_xz_yz[i] = 2.0 * g_xz_x_xz_yz[i] * a_exp;

        g_z_0_0_0_x_x_xz_zz[i] = 2.0 * g_xz_x_xz_zz[i] * a_exp;
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_xz_x_yy_xx, g_xz_x_yy_xy, g_xz_x_yy_xz, g_xz_x_yy_yy, g_xz_x_yy_yz, g_xz_x_yy_zz, g_z_0_0_0_x_x_yy_xx, g_z_0_0_0_x_x_yy_xy, g_z_0_0_0_x_x_yy_xz, g_z_0_0_0_x_x_yy_yy, g_z_0_0_0_x_x_yy_yz, g_z_0_0_0_x_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_x_yy_xx[i] = 2.0 * g_xz_x_yy_xx[i] * a_exp;

        g_z_0_0_0_x_x_yy_xy[i] = 2.0 * g_xz_x_yy_xy[i] * a_exp;

        g_z_0_0_0_x_x_yy_xz[i] = 2.0 * g_xz_x_yy_xz[i] * a_exp;

        g_z_0_0_0_x_x_yy_yy[i] = 2.0 * g_xz_x_yy_yy[i] * a_exp;

        g_z_0_0_0_x_x_yy_yz[i] = 2.0 * g_xz_x_yy_yz[i] * a_exp;

        g_z_0_0_0_x_x_yy_zz[i] = 2.0 * g_xz_x_yy_zz[i] * a_exp;
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_xz_x_yz_xx, g_xz_x_yz_xy, g_xz_x_yz_xz, g_xz_x_yz_yy, g_xz_x_yz_yz, g_xz_x_yz_zz, g_z_0_0_0_x_x_yz_xx, g_z_0_0_0_x_x_yz_xy, g_z_0_0_0_x_x_yz_xz, g_z_0_0_0_x_x_yz_yy, g_z_0_0_0_x_x_yz_yz, g_z_0_0_0_x_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_x_yz_xx[i] = 2.0 * g_xz_x_yz_xx[i] * a_exp;

        g_z_0_0_0_x_x_yz_xy[i] = 2.0 * g_xz_x_yz_xy[i] * a_exp;

        g_z_0_0_0_x_x_yz_xz[i] = 2.0 * g_xz_x_yz_xz[i] * a_exp;

        g_z_0_0_0_x_x_yz_yy[i] = 2.0 * g_xz_x_yz_yy[i] * a_exp;

        g_z_0_0_0_x_x_yz_yz[i] = 2.0 * g_xz_x_yz_yz[i] * a_exp;

        g_z_0_0_0_x_x_yz_zz[i] = 2.0 * g_xz_x_yz_zz[i] * a_exp;
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_xz_x_zz_xx, g_xz_x_zz_xy, g_xz_x_zz_xz, g_xz_x_zz_yy, g_xz_x_zz_yz, g_xz_x_zz_zz, g_z_0_0_0_x_x_zz_xx, g_z_0_0_0_x_x_zz_xy, g_z_0_0_0_x_x_zz_xz, g_z_0_0_0_x_x_zz_yy, g_z_0_0_0_x_x_zz_yz, g_z_0_0_0_x_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_x_zz_xx[i] = 2.0 * g_xz_x_zz_xx[i] * a_exp;

        g_z_0_0_0_x_x_zz_xy[i] = 2.0 * g_xz_x_zz_xy[i] * a_exp;

        g_z_0_0_0_x_x_zz_xz[i] = 2.0 * g_xz_x_zz_xz[i] * a_exp;

        g_z_0_0_0_x_x_zz_yy[i] = 2.0 * g_xz_x_zz_yy[i] * a_exp;

        g_z_0_0_0_x_x_zz_yz[i] = 2.0 * g_xz_x_zz_yz[i] * a_exp;

        g_z_0_0_0_x_x_zz_zz[i] = 2.0 * g_xz_x_zz_zz[i] * a_exp;
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_xz_y_xx_xx, g_xz_y_xx_xy, g_xz_y_xx_xz, g_xz_y_xx_yy, g_xz_y_xx_yz, g_xz_y_xx_zz, g_z_0_0_0_x_y_xx_xx, g_z_0_0_0_x_y_xx_xy, g_z_0_0_0_x_y_xx_xz, g_z_0_0_0_x_y_xx_yy, g_z_0_0_0_x_y_xx_yz, g_z_0_0_0_x_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_y_xx_xx[i] = 2.0 * g_xz_y_xx_xx[i] * a_exp;

        g_z_0_0_0_x_y_xx_xy[i] = 2.0 * g_xz_y_xx_xy[i] * a_exp;

        g_z_0_0_0_x_y_xx_xz[i] = 2.0 * g_xz_y_xx_xz[i] * a_exp;

        g_z_0_0_0_x_y_xx_yy[i] = 2.0 * g_xz_y_xx_yy[i] * a_exp;

        g_z_0_0_0_x_y_xx_yz[i] = 2.0 * g_xz_y_xx_yz[i] * a_exp;

        g_z_0_0_0_x_y_xx_zz[i] = 2.0 * g_xz_y_xx_zz[i] * a_exp;
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_xz_y_xy_xx, g_xz_y_xy_xy, g_xz_y_xy_xz, g_xz_y_xy_yy, g_xz_y_xy_yz, g_xz_y_xy_zz, g_z_0_0_0_x_y_xy_xx, g_z_0_0_0_x_y_xy_xy, g_z_0_0_0_x_y_xy_xz, g_z_0_0_0_x_y_xy_yy, g_z_0_0_0_x_y_xy_yz, g_z_0_0_0_x_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_y_xy_xx[i] = 2.0 * g_xz_y_xy_xx[i] * a_exp;

        g_z_0_0_0_x_y_xy_xy[i] = 2.0 * g_xz_y_xy_xy[i] * a_exp;

        g_z_0_0_0_x_y_xy_xz[i] = 2.0 * g_xz_y_xy_xz[i] * a_exp;

        g_z_0_0_0_x_y_xy_yy[i] = 2.0 * g_xz_y_xy_yy[i] * a_exp;

        g_z_0_0_0_x_y_xy_yz[i] = 2.0 * g_xz_y_xy_yz[i] * a_exp;

        g_z_0_0_0_x_y_xy_zz[i] = 2.0 * g_xz_y_xy_zz[i] * a_exp;
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_xz_y_xz_xx, g_xz_y_xz_xy, g_xz_y_xz_xz, g_xz_y_xz_yy, g_xz_y_xz_yz, g_xz_y_xz_zz, g_z_0_0_0_x_y_xz_xx, g_z_0_0_0_x_y_xz_xy, g_z_0_0_0_x_y_xz_xz, g_z_0_0_0_x_y_xz_yy, g_z_0_0_0_x_y_xz_yz, g_z_0_0_0_x_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_y_xz_xx[i] = 2.0 * g_xz_y_xz_xx[i] * a_exp;

        g_z_0_0_0_x_y_xz_xy[i] = 2.0 * g_xz_y_xz_xy[i] * a_exp;

        g_z_0_0_0_x_y_xz_xz[i] = 2.0 * g_xz_y_xz_xz[i] * a_exp;

        g_z_0_0_0_x_y_xz_yy[i] = 2.0 * g_xz_y_xz_yy[i] * a_exp;

        g_z_0_0_0_x_y_xz_yz[i] = 2.0 * g_xz_y_xz_yz[i] * a_exp;

        g_z_0_0_0_x_y_xz_zz[i] = 2.0 * g_xz_y_xz_zz[i] * a_exp;
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_xz_y_yy_xx, g_xz_y_yy_xy, g_xz_y_yy_xz, g_xz_y_yy_yy, g_xz_y_yy_yz, g_xz_y_yy_zz, g_z_0_0_0_x_y_yy_xx, g_z_0_0_0_x_y_yy_xy, g_z_0_0_0_x_y_yy_xz, g_z_0_0_0_x_y_yy_yy, g_z_0_0_0_x_y_yy_yz, g_z_0_0_0_x_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_y_yy_xx[i] = 2.0 * g_xz_y_yy_xx[i] * a_exp;

        g_z_0_0_0_x_y_yy_xy[i] = 2.0 * g_xz_y_yy_xy[i] * a_exp;

        g_z_0_0_0_x_y_yy_xz[i] = 2.0 * g_xz_y_yy_xz[i] * a_exp;

        g_z_0_0_0_x_y_yy_yy[i] = 2.0 * g_xz_y_yy_yy[i] * a_exp;

        g_z_0_0_0_x_y_yy_yz[i] = 2.0 * g_xz_y_yy_yz[i] * a_exp;

        g_z_0_0_0_x_y_yy_zz[i] = 2.0 * g_xz_y_yy_zz[i] * a_exp;
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_xz_y_yz_xx, g_xz_y_yz_xy, g_xz_y_yz_xz, g_xz_y_yz_yy, g_xz_y_yz_yz, g_xz_y_yz_zz, g_z_0_0_0_x_y_yz_xx, g_z_0_0_0_x_y_yz_xy, g_z_0_0_0_x_y_yz_xz, g_z_0_0_0_x_y_yz_yy, g_z_0_0_0_x_y_yz_yz, g_z_0_0_0_x_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_y_yz_xx[i] = 2.0 * g_xz_y_yz_xx[i] * a_exp;

        g_z_0_0_0_x_y_yz_xy[i] = 2.0 * g_xz_y_yz_xy[i] * a_exp;

        g_z_0_0_0_x_y_yz_xz[i] = 2.0 * g_xz_y_yz_xz[i] * a_exp;

        g_z_0_0_0_x_y_yz_yy[i] = 2.0 * g_xz_y_yz_yy[i] * a_exp;

        g_z_0_0_0_x_y_yz_yz[i] = 2.0 * g_xz_y_yz_yz[i] * a_exp;

        g_z_0_0_0_x_y_yz_zz[i] = 2.0 * g_xz_y_yz_zz[i] * a_exp;
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_xz_y_zz_xx, g_xz_y_zz_xy, g_xz_y_zz_xz, g_xz_y_zz_yy, g_xz_y_zz_yz, g_xz_y_zz_zz, g_z_0_0_0_x_y_zz_xx, g_z_0_0_0_x_y_zz_xy, g_z_0_0_0_x_y_zz_xz, g_z_0_0_0_x_y_zz_yy, g_z_0_0_0_x_y_zz_yz, g_z_0_0_0_x_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_y_zz_xx[i] = 2.0 * g_xz_y_zz_xx[i] * a_exp;

        g_z_0_0_0_x_y_zz_xy[i] = 2.0 * g_xz_y_zz_xy[i] * a_exp;

        g_z_0_0_0_x_y_zz_xz[i] = 2.0 * g_xz_y_zz_xz[i] * a_exp;

        g_z_0_0_0_x_y_zz_yy[i] = 2.0 * g_xz_y_zz_yy[i] * a_exp;

        g_z_0_0_0_x_y_zz_yz[i] = 2.0 * g_xz_y_zz_yz[i] * a_exp;

        g_z_0_0_0_x_y_zz_zz[i] = 2.0 * g_xz_y_zz_zz[i] * a_exp;
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_xz_z_xx_xx, g_xz_z_xx_xy, g_xz_z_xx_xz, g_xz_z_xx_yy, g_xz_z_xx_yz, g_xz_z_xx_zz, g_z_0_0_0_x_z_xx_xx, g_z_0_0_0_x_z_xx_xy, g_z_0_0_0_x_z_xx_xz, g_z_0_0_0_x_z_xx_yy, g_z_0_0_0_x_z_xx_yz, g_z_0_0_0_x_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_z_xx_xx[i] = 2.0 * g_xz_z_xx_xx[i] * a_exp;

        g_z_0_0_0_x_z_xx_xy[i] = 2.0 * g_xz_z_xx_xy[i] * a_exp;

        g_z_0_0_0_x_z_xx_xz[i] = 2.0 * g_xz_z_xx_xz[i] * a_exp;

        g_z_0_0_0_x_z_xx_yy[i] = 2.0 * g_xz_z_xx_yy[i] * a_exp;

        g_z_0_0_0_x_z_xx_yz[i] = 2.0 * g_xz_z_xx_yz[i] * a_exp;

        g_z_0_0_0_x_z_xx_zz[i] = 2.0 * g_xz_z_xx_zz[i] * a_exp;
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_xz_z_xy_xx, g_xz_z_xy_xy, g_xz_z_xy_xz, g_xz_z_xy_yy, g_xz_z_xy_yz, g_xz_z_xy_zz, g_z_0_0_0_x_z_xy_xx, g_z_0_0_0_x_z_xy_xy, g_z_0_0_0_x_z_xy_xz, g_z_0_0_0_x_z_xy_yy, g_z_0_0_0_x_z_xy_yz, g_z_0_0_0_x_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_z_xy_xx[i] = 2.0 * g_xz_z_xy_xx[i] * a_exp;

        g_z_0_0_0_x_z_xy_xy[i] = 2.0 * g_xz_z_xy_xy[i] * a_exp;

        g_z_0_0_0_x_z_xy_xz[i] = 2.0 * g_xz_z_xy_xz[i] * a_exp;

        g_z_0_0_0_x_z_xy_yy[i] = 2.0 * g_xz_z_xy_yy[i] * a_exp;

        g_z_0_0_0_x_z_xy_yz[i] = 2.0 * g_xz_z_xy_yz[i] * a_exp;

        g_z_0_0_0_x_z_xy_zz[i] = 2.0 * g_xz_z_xy_zz[i] * a_exp;
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_xz_z_xz_xx, g_xz_z_xz_xy, g_xz_z_xz_xz, g_xz_z_xz_yy, g_xz_z_xz_yz, g_xz_z_xz_zz, g_z_0_0_0_x_z_xz_xx, g_z_0_0_0_x_z_xz_xy, g_z_0_0_0_x_z_xz_xz, g_z_0_0_0_x_z_xz_yy, g_z_0_0_0_x_z_xz_yz, g_z_0_0_0_x_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_z_xz_xx[i] = 2.0 * g_xz_z_xz_xx[i] * a_exp;

        g_z_0_0_0_x_z_xz_xy[i] = 2.0 * g_xz_z_xz_xy[i] * a_exp;

        g_z_0_0_0_x_z_xz_xz[i] = 2.0 * g_xz_z_xz_xz[i] * a_exp;

        g_z_0_0_0_x_z_xz_yy[i] = 2.0 * g_xz_z_xz_yy[i] * a_exp;

        g_z_0_0_0_x_z_xz_yz[i] = 2.0 * g_xz_z_xz_yz[i] * a_exp;

        g_z_0_0_0_x_z_xz_zz[i] = 2.0 * g_xz_z_xz_zz[i] * a_exp;
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_xz_z_yy_xx, g_xz_z_yy_xy, g_xz_z_yy_xz, g_xz_z_yy_yy, g_xz_z_yy_yz, g_xz_z_yy_zz, g_z_0_0_0_x_z_yy_xx, g_z_0_0_0_x_z_yy_xy, g_z_0_0_0_x_z_yy_xz, g_z_0_0_0_x_z_yy_yy, g_z_0_0_0_x_z_yy_yz, g_z_0_0_0_x_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_z_yy_xx[i] = 2.0 * g_xz_z_yy_xx[i] * a_exp;

        g_z_0_0_0_x_z_yy_xy[i] = 2.0 * g_xz_z_yy_xy[i] * a_exp;

        g_z_0_0_0_x_z_yy_xz[i] = 2.0 * g_xz_z_yy_xz[i] * a_exp;

        g_z_0_0_0_x_z_yy_yy[i] = 2.0 * g_xz_z_yy_yy[i] * a_exp;

        g_z_0_0_0_x_z_yy_yz[i] = 2.0 * g_xz_z_yy_yz[i] * a_exp;

        g_z_0_0_0_x_z_yy_zz[i] = 2.0 * g_xz_z_yy_zz[i] * a_exp;
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_xz_z_yz_xx, g_xz_z_yz_xy, g_xz_z_yz_xz, g_xz_z_yz_yy, g_xz_z_yz_yz, g_xz_z_yz_zz, g_z_0_0_0_x_z_yz_xx, g_z_0_0_0_x_z_yz_xy, g_z_0_0_0_x_z_yz_xz, g_z_0_0_0_x_z_yz_yy, g_z_0_0_0_x_z_yz_yz, g_z_0_0_0_x_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_z_yz_xx[i] = 2.0 * g_xz_z_yz_xx[i] * a_exp;

        g_z_0_0_0_x_z_yz_xy[i] = 2.0 * g_xz_z_yz_xy[i] * a_exp;

        g_z_0_0_0_x_z_yz_xz[i] = 2.0 * g_xz_z_yz_xz[i] * a_exp;

        g_z_0_0_0_x_z_yz_yy[i] = 2.0 * g_xz_z_yz_yy[i] * a_exp;

        g_z_0_0_0_x_z_yz_yz[i] = 2.0 * g_xz_z_yz_yz[i] * a_exp;

        g_z_0_0_0_x_z_yz_zz[i] = 2.0 * g_xz_z_yz_zz[i] * a_exp;
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_xz_z_zz_xx, g_xz_z_zz_xy, g_xz_z_zz_xz, g_xz_z_zz_yy, g_xz_z_zz_yz, g_xz_z_zz_zz, g_z_0_0_0_x_z_zz_xx, g_z_0_0_0_x_z_zz_xy, g_z_0_0_0_x_z_zz_xz, g_z_0_0_0_x_z_zz_yy, g_z_0_0_0_x_z_zz_yz, g_z_0_0_0_x_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_z_zz_xx[i] = 2.0 * g_xz_z_zz_xx[i] * a_exp;

        g_z_0_0_0_x_z_zz_xy[i] = 2.0 * g_xz_z_zz_xy[i] * a_exp;

        g_z_0_0_0_x_z_zz_xz[i] = 2.0 * g_xz_z_zz_xz[i] * a_exp;

        g_z_0_0_0_x_z_zz_yy[i] = 2.0 * g_xz_z_zz_yy[i] * a_exp;

        g_z_0_0_0_x_z_zz_yz[i] = 2.0 * g_xz_z_zz_yz[i] * a_exp;

        g_z_0_0_0_x_z_zz_zz[i] = 2.0 * g_xz_z_zz_zz[i] * a_exp;
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_yz_x_xx_xx, g_yz_x_xx_xy, g_yz_x_xx_xz, g_yz_x_xx_yy, g_yz_x_xx_yz, g_yz_x_xx_zz, g_z_0_0_0_y_x_xx_xx, g_z_0_0_0_y_x_xx_xy, g_z_0_0_0_y_x_xx_xz, g_z_0_0_0_y_x_xx_yy, g_z_0_0_0_y_x_xx_yz, g_z_0_0_0_y_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_x_xx_xx[i] = 2.0 * g_yz_x_xx_xx[i] * a_exp;

        g_z_0_0_0_y_x_xx_xy[i] = 2.0 * g_yz_x_xx_xy[i] * a_exp;

        g_z_0_0_0_y_x_xx_xz[i] = 2.0 * g_yz_x_xx_xz[i] * a_exp;

        g_z_0_0_0_y_x_xx_yy[i] = 2.0 * g_yz_x_xx_yy[i] * a_exp;

        g_z_0_0_0_y_x_xx_yz[i] = 2.0 * g_yz_x_xx_yz[i] * a_exp;

        g_z_0_0_0_y_x_xx_zz[i] = 2.0 * g_yz_x_xx_zz[i] * a_exp;
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_yz_x_xy_xx, g_yz_x_xy_xy, g_yz_x_xy_xz, g_yz_x_xy_yy, g_yz_x_xy_yz, g_yz_x_xy_zz, g_z_0_0_0_y_x_xy_xx, g_z_0_0_0_y_x_xy_xy, g_z_0_0_0_y_x_xy_xz, g_z_0_0_0_y_x_xy_yy, g_z_0_0_0_y_x_xy_yz, g_z_0_0_0_y_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_x_xy_xx[i] = 2.0 * g_yz_x_xy_xx[i] * a_exp;

        g_z_0_0_0_y_x_xy_xy[i] = 2.0 * g_yz_x_xy_xy[i] * a_exp;

        g_z_0_0_0_y_x_xy_xz[i] = 2.0 * g_yz_x_xy_xz[i] * a_exp;

        g_z_0_0_0_y_x_xy_yy[i] = 2.0 * g_yz_x_xy_yy[i] * a_exp;

        g_z_0_0_0_y_x_xy_yz[i] = 2.0 * g_yz_x_xy_yz[i] * a_exp;

        g_z_0_0_0_y_x_xy_zz[i] = 2.0 * g_yz_x_xy_zz[i] * a_exp;
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_yz_x_xz_xx, g_yz_x_xz_xy, g_yz_x_xz_xz, g_yz_x_xz_yy, g_yz_x_xz_yz, g_yz_x_xz_zz, g_z_0_0_0_y_x_xz_xx, g_z_0_0_0_y_x_xz_xy, g_z_0_0_0_y_x_xz_xz, g_z_0_0_0_y_x_xz_yy, g_z_0_0_0_y_x_xz_yz, g_z_0_0_0_y_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_x_xz_xx[i] = 2.0 * g_yz_x_xz_xx[i] * a_exp;

        g_z_0_0_0_y_x_xz_xy[i] = 2.0 * g_yz_x_xz_xy[i] * a_exp;

        g_z_0_0_0_y_x_xz_xz[i] = 2.0 * g_yz_x_xz_xz[i] * a_exp;

        g_z_0_0_0_y_x_xz_yy[i] = 2.0 * g_yz_x_xz_yy[i] * a_exp;

        g_z_0_0_0_y_x_xz_yz[i] = 2.0 * g_yz_x_xz_yz[i] * a_exp;

        g_z_0_0_0_y_x_xz_zz[i] = 2.0 * g_yz_x_xz_zz[i] * a_exp;
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_yz_x_yy_xx, g_yz_x_yy_xy, g_yz_x_yy_xz, g_yz_x_yy_yy, g_yz_x_yy_yz, g_yz_x_yy_zz, g_z_0_0_0_y_x_yy_xx, g_z_0_0_0_y_x_yy_xy, g_z_0_0_0_y_x_yy_xz, g_z_0_0_0_y_x_yy_yy, g_z_0_0_0_y_x_yy_yz, g_z_0_0_0_y_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_x_yy_xx[i] = 2.0 * g_yz_x_yy_xx[i] * a_exp;

        g_z_0_0_0_y_x_yy_xy[i] = 2.0 * g_yz_x_yy_xy[i] * a_exp;

        g_z_0_0_0_y_x_yy_xz[i] = 2.0 * g_yz_x_yy_xz[i] * a_exp;

        g_z_0_0_0_y_x_yy_yy[i] = 2.0 * g_yz_x_yy_yy[i] * a_exp;

        g_z_0_0_0_y_x_yy_yz[i] = 2.0 * g_yz_x_yy_yz[i] * a_exp;

        g_z_0_0_0_y_x_yy_zz[i] = 2.0 * g_yz_x_yy_zz[i] * a_exp;
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_yz_x_yz_xx, g_yz_x_yz_xy, g_yz_x_yz_xz, g_yz_x_yz_yy, g_yz_x_yz_yz, g_yz_x_yz_zz, g_z_0_0_0_y_x_yz_xx, g_z_0_0_0_y_x_yz_xy, g_z_0_0_0_y_x_yz_xz, g_z_0_0_0_y_x_yz_yy, g_z_0_0_0_y_x_yz_yz, g_z_0_0_0_y_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_x_yz_xx[i] = 2.0 * g_yz_x_yz_xx[i] * a_exp;

        g_z_0_0_0_y_x_yz_xy[i] = 2.0 * g_yz_x_yz_xy[i] * a_exp;

        g_z_0_0_0_y_x_yz_xz[i] = 2.0 * g_yz_x_yz_xz[i] * a_exp;

        g_z_0_0_0_y_x_yz_yy[i] = 2.0 * g_yz_x_yz_yy[i] * a_exp;

        g_z_0_0_0_y_x_yz_yz[i] = 2.0 * g_yz_x_yz_yz[i] * a_exp;

        g_z_0_0_0_y_x_yz_zz[i] = 2.0 * g_yz_x_yz_zz[i] * a_exp;
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_yz_x_zz_xx, g_yz_x_zz_xy, g_yz_x_zz_xz, g_yz_x_zz_yy, g_yz_x_zz_yz, g_yz_x_zz_zz, g_z_0_0_0_y_x_zz_xx, g_z_0_0_0_y_x_zz_xy, g_z_0_0_0_y_x_zz_xz, g_z_0_0_0_y_x_zz_yy, g_z_0_0_0_y_x_zz_yz, g_z_0_0_0_y_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_x_zz_xx[i] = 2.0 * g_yz_x_zz_xx[i] * a_exp;

        g_z_0_0_0_y_x_zz_xy[i] = 2.0 * g_yz_x_zz_xy[i] * a_exp;

        g_z_0_0_0_y_x_zz_xz[i] = 2.0 * g_yz_x_zz_xz[i] * a_exp;

        g_z_0_0_0_y_x_zz_yy[i] = 2.0 * g_yz_x_zz_yy[i] * a_exp;

        g_z_0_0_0_y_x_zz_yz[i] = 2.0 * g_yz_x_zz_yz[i] * a_exp;

        g_z_0_0_0_y_x_zz_zz[i] = 2.0 * g_yz_x_zz_zz[i] * a_exp;
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_yz_y_xx_xx, g_yz_y_xx_xy, g_yz_y_xx_xz, g_yz_y_xx_yy, g_yz_y_xx_yz, g_yz_y_xx_zz, g_z_0_0_0_y_y_xx_xx, g_z_0_0_0_y_y_xx_xy, g_z_0_0_0_y_y_xx_xz, g_z_0_0_0_y_y_xx_yy, g_z_0_0_0_y_y_xx_yz, g_z_0_0_0_y_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_y_xx_xx[i] = 2.0 * g_yz_y_xx_xx[i] * a_exp;

        g_z_0_0_0_y_y_xx_xy[i] = 2.0 * g_yz_y_xx_xy[i] * a_exp;

        g_z_0_0_0_y_y_xx_xz[i] = 2.0 * g_yz_y_xx_xz[i] * a_exp;

        g_z_0_0_0_y_y_xx_yy[i] = 2.0 * g_yz_y_xx_yy[i] * a_exp;

        g_z_0_0_0_y_y_xx_yz[i] = 2.0 * g_yz_y_xx_yz[i] * a_exp;

        g_z_0_0_0_y_y_xx_zz[i] = 2.0 * g_yz_y_xx_zz[i] * a_exp;
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_yz_y_xy_xx, g_yz_y_xy_xy, g_yz_y_xy_xz, g_yz_y_xy_yy, g_yz_y_xy_yz, g_yz_y_xy_zz, g_z_0_0_0_y_y_xy_xx, g_z_0_0_0_y_y_xy_xy, g_z_0_0_0_y_y_xy_xz, g_z_0_0_0_y_y_xy_yy, g_z_0_0_0_y_y_xy_yz, g_z_0_0_0_y_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_y_xy_xx[i] = 2.0 * g_yz_y_xy_xx[i] * a_exp;

        g_z_0_0_0_y_y_xy_xy[i] = 2.0 * g_yz_y_xy_xy[i] * a_exp;

        g_z_0_0_0_y_y_xy_xz[i] = 2.0 * g_yz_y_xy_xz[i] * a_exp;

        g_z_0_0_0_y_y_xy_yy[i] = 2.0 * g_yz_y_xy_yy[i] * a_exp;

        g_z_0_0_0_y_y_xy_yz[i] = 2.0 * g_yz_y_xy_yz[i] * a_exp;

        g_z_0_0_0_y_y_xy_zz[i] = 2.0 * g_yz_y_xy_zz[i] * a_exp;
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_yz_y_xz_xx, g_yz_y_xz_xy, g_yz_y_xz_xz, g_yz_y_xz_yy, g_yz_y_xz_yz, g_yz_y_xz_zz, g_z_0_0_0_y_y_xz_xx, g_z_0_0_0_y_y_xz_xy, g_z_0_0_0_y_y_xz_xz, g_z_0_0_0_y_y_xz_yy, g_z_0_0_0_y_y_xz_yz, g_z_0_0_0_y_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_y_xz_xx[i] = 2.0 * g_yz_y_xz_xx[i] * a_exp;

        g_z_0_0_0_y_y_xz_xy[i] = 2.0 * g_yz_y_xz_xy[i] * a_exp;

        g_z_0_0_0_y_y_xz_xz[i] = 2.0 * g_yz_y_xz_xz[i] * a_exp;

        g_z_0_0_0_y_y_xz_yy[i] = 2.0 * g_yz_y_xz_yy[i] * a_exp;

        g_z_0_0_0_y_y_xz_yz[i] = 2.0 * g_yz_y_xz_yz[i] * a_exp;

        g_z_0_0_0_y_y_xz_zz[i] = 2.0 * g_yz_y_xz_zz[i] * a_exp;
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_yz_y_yy_xx, g_yz_y_yy_xy, g_yz_y_yy_xz, g_yz_y_yy_yy, g_yz_y_yy_yz, g_yz_y_yy_zz, g_z_0_0_0_y_y_yy_xx, g_z_0_0_0_y_y_yy_xy, g_z_0_0_0_y_y_yy_xz, g_z_0_0_0_y_y_yy_yy, g_z_0_0_0_y_y_yy_yz, g_z_0_0_0_y_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_y_yy_xx[i] = 2.0 * g_yz_y_yy_xx[i] * a_exp;

        g_z_0_0_0_y_y_yy_xy[i] = 2.0 * g_yz_y_yy_xy[i] * a_exp;

        g_z_0_0_0_y_y_yy_xz[i] = 2.0 * g_yz_y_yy_xz[i] * a_exp;

        g_z_0_0_0_y_y_yy_yy[i] = 2.0 * g_yz_y_yy_yy[i] * a_exp;

        g_z_0_0_0_y_y_yy_yz[i] = 2.0 * g_yz_y_yy_yz[i] * a_exp;

        g_z_0_0_0_y_y_yy_zz[i] = 2.0 * g_yz_y_yy_zz[i] * a_exp;
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_yz_y_yz_xx, g_yz_y_yz_xy, g_yz_y_yz_xz, g_yz_y_yz_yy, g_yz_y_yz_yz, g_yz_y_yz_zz, g_z_0_0_0_y_y_yz_xx, g_z_0_0_0_y_y_yz_xy, g_z_0_0_0_y_y_yz_xz, g_z_0_0_0_y_y_yz_yy, g_z_0_0_0_y_y_yz_yz, g_z_0_0_0_y_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_y_yz_xx[i] = 2.0 * g_yz_y_yz_xx[i] * a_exp;

        g_z_0_0_0_y_y_yz_xy[i] = 2.0 * g_yz_y_yz_xy[i] * a_exp;

        g_z_0_0_0_y_y_yz_xz[i] = 2.0 * g_yz_y_yz_xz[i] * a_exp;

        g_z_0_0_0_y_y_yz_yy[i] = 2.0 * g_yz_y_yz_yy[i] * a_exp;

        g_z_0_0_0_y_y_yz_yz[i] = 2.0 * g_yz_y_yz_yz[i] * a_exp;

        g_z_0_0_0_y_y_yz_zz[i] = 2.0 * g_yz_y_yz_zz[i] * a_exp;
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_yz_y_zz_xx, g_yz_y_zz_xy, g_yz_y_zz_xz, g_yz_y_zz_yy, g_yz_y_zz_yz, g_yz_y_zz_zz, g_z_0_0_0_y_y_zz_xx, g_z_0_0_0_y_y_zz_xy, g_z_0_0_0_y_y_zz_xz, g_z_0_0_0_y_y_zz_yy, g_z_0_0_0_y_y_zz_yz, g_z_0_0_0_y_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_y_zz_xx[i] = 2.0 * g_yz_y_zz_xx[i] * a_exp;

        g_z_0_0_0_y_y_zz_xy[i] = 2.0 * g_yz_y_zz_xy[i] * a_exp;

        g_z_0_0_0_y_y_zz_xz[i] = 2.0 * g_yz_y_zz_xz[i] * a_exp;

        g_z_0_0_0_y_y_zz_yy[i] = 2.0 * g_yz_y_zz_yy[i] * a_exp;

        g_z_0_0_0_y_y_zz_yz[i] = 2.0 * g_yz_y_zz_yz[i] * a_exp;

        g_z_0_0_0_y_y_zz_zz[i] = 2.0 * g_yz_y_zz_zz[i] * a_exp;
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_yz_z_xx_xx, g_yz_z_xx_xy, g_yz_z_xx_xz, g_yz_z_xx_yy, g_yz_z_xx_yz, g_yz_z_xx_zz, g_z_0_0_0_y_z_xx_xx, g_z_0_0_0_y_z_xx_xy, g_z_0_0_0_y_z_xx_xz, g_z_0_0_0_y_z_xx_yy, g_z_0_0_0_y_z_xx_yz, g_z_0_0_0_y_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_z_xx_xx[i] = 2.0 * g_yz_z_xx_xx[i] * a_exp;

        g_z_0_0_0_y_z_xx_xy[i] = 2.0 * g_yz_z_xx_xy[i] * a_exp;

        g_z_0_0_0_y_z_xx_xz[i] = 2.0 * g_yz_z_xx_xz[i] * a_exp;

        g_z_0_0_0_y_z_xx_yy[i] = 2.0 * g_yz_z_xx_yy[i] * a_exp;

        g_z_0_0_0_y_z_xx_yz[i] = 2.0 * g_yz_z_xx_yz[i] * a_exp;

        g_z_0_0_0_y_z_xx_zz[i] = 2.0 * g_yz_z_xx_zz[i] * a_exp;
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_yz_z_xy_xx, g_yz_z_xy_xy, g_yz_z_xy_xz, g_yz_z_xy_yy, g_yz_z_xy_yz, g_yz_z_xy_zz, g_z_0_0_0_y_z_xy_xx, g_z_0_0_0_y_z_xy_xy, g_z_0_0_0_y_z_xy_xz, g_z_0_0_0_y_z_xy_yy, g_z_0_0_0_y_z_xy_yz, g_z_0_0_0_y_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_z_xy_xx[i] = 2.0 * g_yz_z_xy_xx[i] * a_exp;

        g_z_0_0_0_y_z_xy_xy[i] = 2.0 * g_yz_z_xy_xy[i] * a_exp;

        g_z_0_0_0_y_z_xy_xz[i] = 2.0 * g_yz_z_xy_xz[i] * a_exp;

        g_z_0_0_0_y_z_xy_yy[i] = 2.0 * g_yz_z_xy_yy[i] * a_exp;

        g_z_0_0_0_y_z_xy_yz[i] = 2.0 * g_yz_z_xy_yz[i] * a_exp;

        g_z_0_0_0_y_z_xy_zz[i] = 2.0 * g_yz_z_xy_zz[i] * a_exp;
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_yz_z_xz_xx, g_yz_z_xz_xy, g_yz_z_xz_xz, g_yz_z_xz_yy, g_yz_z_xz_yz, g_yz_z_xz_zz, g_z_0_0_0_y_z_xz_xx, g_z_0_0_0_y_z_xz_xy, g_z_0_0_0_y_z_xz_xz, g_z_0_0_0_y_z_xz_yy, g_z_0_0_0_y_z_xz_yz, g_z_0_0_0_y_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_z_xz_xx[i] = 2.0 * g_yz_z_xz_xx[i] * a_exp;

        g_z_0_0_0_y_z_xz_xy[i] = 2.0 * g_yz_z_xz_xy[i] * a_exp;

        g_z_0_0_0_y_z_xz_xz[i] = 2.0 * g_yz_z_xz_xz[i] * a_exp;

        g_z_0_0_0_y_z_xz_yy[i] = 2.0 * g_yz_z_xz_yy[i] * a_exp;

        g_z_0_0_0_y_z_xz_yz[i] = 2.0 * g_yz_z_xz_yz[i] * a_exp;

        g_z_0_0_0_y_z_xz_zz[i] = 2.0 * g_yz_z_xz_zz[i] * a_exp;
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_yz_z_yy_xx, g_yz_z_yy_xy, g_yz_z_yy_xz, g_yz_z_yy_yy, g_yz_z_yy_yz, g_yz_z_yy_zz, g_z_0_0_0_y_z_yy_xx, g_z_0_0_0_y_z_yy_xy, g_z_0_0_0_y_z_yy_xz, g_z_0_0_0_y_z_yy_yy, g_z_0_0_0_y_z_yy_yz, g_z_0_0_0_y_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_z_yy_xx[i] = 2.0 * g_yz_z_yy_xx[i] * a_exp;

        g_z_0_0_0_y_z_yy_xy[i] = 2.0 * g_yz_z_yy_xy[i] * a_exp;

        g_z_0_0_0_y_z_yy_xz[i] = 2.0 * g_yz_z_yy_xz[i] * a_exp;

        g_z_0_0_0_y_z_yy_yy[i] = 2.0 * g_yz_z_yy_yy[i] * a_exp;

        g_z_0_0_0_y_z_yy_yz[i] = 2.0 * g_yz_z_yy_yz[i] * a_exp;

        g_z_0_0_0_y_z_yy_zz[i] = 2.0 * g_yz_z_yy_zz[i] * a_exp;
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_yz_z_yz_xx, g_yz_z_yz_xy, g_yz_z_yz_xz, g_yz_z_yz_yy, g_yz_z_yz_yz, g_yz_z_yz_zz, g_z_0_0_0_y_z_yz_xx, g_z_0_0_0_y_z_yz_xy, g_z_0_0_0_y_z_yz_xz, g_z_0_0_0_y_z_yz_yy, g_z_0_0_0_y_z_yz_yz, g_z_0_0_0_y_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_z_yz_xx[i] = 2.0 * g_yz_z_yz_xx[i] * a_exp;

        g_z_0_0_0_y_z_yz_xy[i] = 2.0 * g_yz_z_yz_xy[i] * a_exp;

        g_z_0_0_0_y_z_yz_xz[i] = 2.0 * g_yz_z_yz_xz[i] * a_exp;

        g_z_0_0_0_y_z_yz_yy[i] = 2.0 * g_yz_z_yz_yy[i] * a_exp;

        g_z_0_0_0_y_z_yz_yz[i] = 2.0 * g_yz_z_yz_yz[i] * a_exp;

        g_z_0_0_0_y_z_yz_zz[i] = 2.0 * g_yz_z_yz_zz[i] * a_exp;
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_yz_z_zz_xx, g_yz_z_zz_xy, g_yz_z_zz_xz, g_yz_z_zz_yy, g_yz_z_zz_yz, g_yz_z_zz_zz, g_z_0_0_0_y_z_zz_xx, g_z_0_0_0_y_z_zz_xy, g_z_0_0_0_y_z_zz_xz, g_z_0_0_0_y_z_zz_yy, g_z_0_0_0_y_z_zz_yz, g_z_0_0_0_y_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_z_zz_xx[i] = 2.0 * g_yz_z_zz_xx[i] * a_exp;

        g_z_0_0_0_y_z_zz_xy[i] = 2.0 * g_yz_z_zz_xy[i] * a_exp;

        g_z_0_0_0_y_z_zz_xz[i] = 2.0 * g_yz_z_zz_xz[i] * a_exp;

        g_z_0_0_0_y_z_zz_yy[i] = 2.0 * g_yz_z_zz_yy[i] * a_exp;

        g_z_0_0_0_y_z_zz_yz[i] = 2.0 * g_yz_z_zz_yz[i] * a_exp;

        g_z_0_0_0_y_z_zz_zz[i] = 2.0 * g_yz_z_zz_zz[i] * a_exp;
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_0_x_xx_xx, g_0_x_xx_xy, g_0_x_xx_xz, g_0_x_xx_yy, g_0_x_xx_yz, g_0_x_xx_zz, g_z_0_0_0_z_x_xx_xx, g_z_0_0_0_z_x_xx_xy, g_z_0_0_0_z_x_xx_xz, g_z_0_0_0_z_x_xx_yy, g_z_0_0_0_z_x_xx_yz, g_z_0_0_0_z_x_xx_zz, g_zz_x_xx_xx, g_zz_x_xx_xy, g_zz_x_xx_xz, g_zz_x_xx_yy, g_zz_x_xx_yz, g_zz_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_x_xx_xx[i] = -g_0_x_xx_xx[i] + 2.0 * g_zz_x_xx_xx[i] * a_exp;

        g_z_0_0_0_z_x_xx_xy[i] = -g_0_x_xx_xy[i] + 2.0 * g_zz_x_xx_xy[i] * a_exp;

        g_z_0_0_0_z_x_xx_xz[i] = -g_0_x_xx_xz[i] + 2.0 * g_zz_x_xx_xz[i] * a_exp;

        g_z_0_0_0_z_x_xx_yy[i] = -g_0_x_xx_yy[i] + 2.0 * g_zz_x_xx_yy[i] * a_exp;

        g_z_0_0_0_z_x_xx_yz[i] = -g_0_x_xx_yz[i] + 2.0 * g_zz_x_xx_yz[i] * a_exp;

        g_z_0_0_0_z_x_xx_zz[i] = -g_0_x_xx_zz[i] + 2.0 * g_zz_x_xx_zz[i] * a_exp;
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_0_x_xy_xx, g_0_x_xy_xy, g_0_x_xy_xz, g_0_x_xy_yy, g_0_x_xy_yz, g_0_x_xy_zz, g_z_0_0_0_z_x_xy_xx, g_z_0_0_0_z_x_xy_xy, g_z_0_0_0_z_x_xy_xz, g_z_0_0_0_z_x_xy_yy, g_z_0_0_0_z_x_xy_yz, g_z_0_0_0_z_x_xy_zz, g_zz_x_xy_xx, g_zz_x_xy_xy, g_zz_x_xy_xz, g_zz_x_xy_yy, g_zz_x_xy_yz, g_zz_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_x_xy_xx[i] = -g_0_x_xy_xx[i] + 2.0 * g_zz_x_xy_xx[i] * a_exp;

        g_z_0_0_0_z_x_xy_xy[i] = -g_0_x_xy_xy[i] + 2.0 * g_zz_x_xy_xy[i] * a_exp;

        g_z_0_0_0_z_x_xy_xz[i] = -g_0_x_xy_xz[i] + 2.0 * g_zz_x_xy_xz[i] * a_exp;

        g_z_0_0_0_z_x_xy_yy[i] = -g_0_x_xy_yy[i] + 2.0 * g_zz_x_xy_yy[i] * a_exp;

        g_z_0_0_0_z_x_xy_yz[i] = -g_0_x_xy_yz[i] + 2.0 * g_zz_x_xy_yz[i] * a_exp;

        g_z_0_0_0_z_x_xy_zz[i] = -g_0_x_xy_zz[i] + 2.0 * g_zz_x_xy_zz[i] * a_exp;
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_0_x_xz_xx, g_0_x_xz_xy, g_0_x_xz_xz, g_0_x_xz_yy, g_0_x_xz_yz, g_0_x_xz_zz, g_z_0_0_0_z_x_xz_xx, g_z_0_0_0_z_x_xz_xy, g_z_0_0_0_z_x_xz_xz, g_z_0_0_0_z_x_xz_yy, g_z_0_0_0_z_x_xz_yz, g_z_0_0_0_z_x_xz_zz, g_zz_x_xz_xx, g_zz_x_xz_xy, g_zz_x_xz_xz, g_zz_x_xz_yy, g_zz_x_xz_yz, g_zz_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_x_xz_xx[i] = -g_0_x_xz_xx[i] + 2.0 * g_zz_x_xz_xx[i] * a_exp;

        g_z_0_0_0_z_x_xz_xy[i] = -g_0_x_xz_xy[i] + 2.0 * g_zz_x_xz_xy[i] * a_exp;

        g_z_0_0_0_z_x_xz_xz[i] = -g_0_x_xz_xz[i] + 2.0 * g_zz_x_xz_xz[i] * a_exp;

        g_z_0_0_0_z_x_xz_yy[i] = -g_0_x_xz_yy[i] + 2.0 * g_zz_x_xz_yy[i] * a_exp;

        g_z_0_0_0_z_x_xz_yz[i] = -g_0_x_xz_yz[i] + 2.0 * g_zz_x_xz_yz[i] * a_exp;

        g_z_0_0_0_z_x_xz_zz[i] = -g_0_x_xz_zz[i] + 2.0 * g_zz_x_xz_zz[i] * a_exp;
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_0_x_yy_xx, g_0_x_yy_xy, g_0_x_yy_xz, g_0_x_yy_yy, g_0_x_yy_yz, g_0_x_yy_zz, g_z_0_0_0_z_x_yy_xx, g_z_0_0_0_z_x_yy_xy, g_z_0_0_0_z_x_yy_xz, g_z_0_0_0_z_x_yy_yy, g_z_0_0_0_z_x_yy_yz, g_z_0_0_0_z_x_yy_zz, g_zz_x_yy_xx, g_zz_x_yy_xy, g_zz_x_yy_xz, g_zz_x_yy_yy, g_zz_x_yy_yz, g_zz_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_x_yy_xx[i] = -g_0_x_yy_xx[i] + 2.0 * g_zz_x_yy_xx[i] * a_exp;

        g_z_0_0_0_z_x_yy_xy[i] = -g_0_x_yy_xy[i] + 2.0 * g_zz_x_yy_xy[i] * a_exp;

        g_z_0_0_0_z_x_yy_xz[i] = -g_0_x_yy_xz[i] + 2.0 * g_zz_x_yy_xz[i] * a_exp;

        g_z_0_0_0_z_x_yy_yy[i] = -g_0_x_yy_yy[i] + 2.0 * g_zz_x_yy_yy[i] * a_exp;

        g_z_0_0_0_z_x_yy_yz[i] = -g_0_x_yy_yz[i] + 2.0 * g_zz_x_yy_yz[i] * a_exp;

        g_z_0_0_0_z_x_yy_zz[i] = -g_0_x_yy_zz[i] + 2.0 * g_zz_x_yy_zz[i] * a_exp;
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_0_x_yz_xx, g_0_x_yz_xy, g_0_x_yz_xz, g_0_x_yz_yy, g_0_x_yz_yz, g_0_x_yz_zz, g_z_0_0_0_z_x_yz_xx, g_z_0_0_0_z_x_yz_xy, g_z_0_0_0_z_x_yz_xz, g_z_0_0_0_z_x_yz_yy, g_z_0_0_0_z_x_yz_yz, g_z_0_0_0_z_x_yz_zz, g_zz_x_yz_xx, g_zz_x_yz_xy, g_zz_x_yz_xz, g_zz_x_yz_yy, g_zz_x_yz_yz, g_zz_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_x_yz_xx[i] = -g_0_x_yz_xx[i] + 2.0 * g_zz_x_yz_xx[i] * a_exp;

        g_z_0_0_0_z_x_yz_xy[i] = -g_0_x_yz_xy[i] + 2.0 * g_zz_x_yz_xy[i] * a_exp;

        g_z_0_0_0_z_x_yz_xz[i] = -g_0_x_yz_xz[i] + 2.0 * g_zz_x_yz_xz[i] * a_exp;

        g_z_0_0_0_z_x_yz_yy[i] = -g_0_x_yz_yy[i] + 2.0 * g_zz_x_yz_yy[i] * a_exp;

        g_z_0_0_0_z_x_yz_yz[i] = -g_0_x_yz_yz[i] + 2.0 * g_zz_x_yz_yz[i] * a_exp;

        g_z_0_0_0_z_x_yz_zz[i] = -g_0_x_yz_zz[i] + 2.0 * g_zz_x_yz_zz[i] * a_exp;
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_0_x_zz_xx, g_0_x_zz_xy, g_0_x_zz_xz, g_0_x_zz_yy, g_0_x_zz_yz, g_0_x_zz_zz, g_z_0_0_0_z_x_zz_xx, g_z_0_0_0_z_x_zz_xy, g_z_0_0_0_z_x_zz_xz, g_z_0_0_0_z_x_zz_yy, g_z_0_0_0_z_x_zz_yz, g_z_0_0_0_z_x_zz_zz, g_zz_x_zz_xx, g_zz_x_zz_xy, g_zz_x_zz_xz, g_zz_x_zz_yy, g_zz_x_zz_yz, g_zz_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_x_zz_xx[i] = -g_0_x_zz_xx[i] + 2.0 * g_zz_x_zz_xx[i] * a_exp;

        g_z_0_0_0_z_x_zz_xy[i] = -g_0_x_zz_xy[i] + 2.0 * g_zz_x_zz_xy[i] * a_exp;

        g_z_0_0_0_z_x_zz_xz[i] = -g_0_x_zz_xz[i] + 2.0 * g_zz_x_zz_xz[i] * a_exp;

        g_z_0_0_0_z_x_zz_yy[i] = -g_0_x_zz_yy[i] + 2.0 * g_zz_x_zz_yy[i] * a_exp;

        g_z_0_0_0_z_x_zz_yz[i] = -g_0_x_zz_yz[i] + 2.0 * g_zz_x_zz_yz[i] * a_exp;

        g_z_0_0_0_z_x_zz_zz[i] = -g_0_x_zz_zz[i] + 2.0 * g_zz_x_zz_zz[i] * a_exp;
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_0_y_xx_xx, g_0_y_xx_xy, g_0_y_xx_xz, g_0_y_xx_yy, g_0_y_xx_yz, g_0_y_xx_zz, g_z_0_0_0_z_y_xx_xx, g_z_0_0_0_z_y_xx_xy, g_z_0_0_0_z_y_xx_xz, g_z_0_0_0_z_y_xx_yy, g_z_0_0_0_z_y_xx_yz, g_z_0_0_0_z_y_xx_zz, g_zz_y_xx_xx, g_zz_y_xx_xy, g_zz_y_xx_xz, g_zz_y_xx_yy, g_zz_y_xx_yz, g_zz_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_y_xx_xx[i] = -g_0_y_xx_xx[i] + 2.0 * g_zz_y_xx_xx[i] * a_exp;

        g_z_0_0_0_z_y_xx_xy[i] = -g_0_y_xx_xy[i] + 2.0 * g_zz_y_xx_xy[i] * a_exp;

        g_z_0_0_0_z_y_xx_xz[i] = -g_0_y_xx_xz[i] + 2.0 * g_zz_y_xx_xz[i] * a_exp;

        g_z_0_0_0_z_y_xx_yy[i] = -g_0_y_xx_yy[i] + 2.0 * g_zz_y_xx_yy[i] * a_exp;

        g_z_0_0_0_z_y_xx_yz[i] = -g_0_y_xx_yz[i] + 2.0 * g_zz_y_xx_yz[i] * a_exp;

        g_z_0_0_0_z_y_xx_zz[i] = -g_0_y_xx_zz[i] + 2.0 * g_zz_y_xx_zz[i] * a_exp;
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_0_y_xy_xx, g_0_y_xy_xy, g_0_y_xy_xz, g_0_y_xy_yy, g_0_y_xy_yz, g_0_y_xy_zz, g_z_0_0_0_z_y_xy_xx, g_z_0_0_0_z_y_xy_xy, g_z_0_0_0_z_y_xy_xz, g_z_0_0_0_z_y_xy_yy, g_z_0_0_0_z_y_xy_yz, g_z_0_0_0_z_y_xy_zz, g_zz_y_xy_xx, g_zz_y_xy_xy, g_zz_y_xy_xz, g_zz_y_xy_yy, g_zz_y_xy_yz, g_zz_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_y_xy_xx[i] = -g_0_y_xy_xx[i] + 2.0 * g_zz_y_xy_xx[i] * a_exp;

        g_z_0_0_0_z_y_xy_xy[i] = -g_0_y_xy_xy[i] + 2.0 * g_zz_y_xy_xy[i] * a_exp;

        g_z_0_0_0_z_y_xy_xz[i] = -g_0_y_xy_xz[i] + 2.0 * g_zz_y_xy_xz[i] * a_exp;

        g_z_0_0_0_z_y_xy_yy[i] = -g_0_y_xy_yy[i] + 2.0 * g_zz_y_xy_yy[i] * a_exp;

        g_z_0_0_0_z_y_xy_yz[i] = -g_0_y_xy_yz[i] + 2.0 * g_zz_y_xy_yz[i] * a_exp;

        g_z_0_0_0_z_y_xy_zz[i] = -g_0_y_xy_zz[i] + 2.0 * g_zz_y_xy_zz[i] * a_exp;
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_0_y_xz_xx, g_0_y_xz_xy, g_0_y_xz_xz, g_0_y_xz_yy, g_0_y_xz_yz, g_0_y_xz_zz, g_z_0_0_0_z_y_xz_xx, g_z_0_0_0_z_y_xz_xy, g_z_0_0_0_z_y_xz_xz, g_z_0_0_0_z_y_xz_yy, g_z_0_0_0_z_y_xz_yz, g_z_0_0_0_z_y_xz_zz, g_zz_y_xz_xx, g_zz_y_xz_xy, g_zz_y_xz_xz, g_zz_y_xz_yy, g_zz_y_xz_yz, g_zz_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_y_xz_xx[i] = -g_0_y_xz_xx[i] + 2.0 * g_zz_y_xz_xx[i] * a_exp;

        g_z_0_0_0_z_y_xz_xy[i] = -g_0_y_xz_xy[i] + 2.0 * g_zz_y_xz_xy[i] * a_exp;

        g_z_0_0_0_z_y_xz_xz[i] = -g_0_y_xz_xz[i] + 2.0 * g_zz_y_xz_xz[i] * a_exp;

        g_z_0_0_0_z_y_xz_yy[i] = -g_0_y_xz_yy[i] + 2.0 * g_zz_y_xz_yy[i] * a_exp;

        g_z_0_0_0_z_y_xz_yz[i] = -g_0_y_xz_yz[i] + 2.0 * g_zz_y_xz_yz[i] * a_exp;

        g_z_0_0_0_z_y_xz_zz[i] = -g_0_y_xz_zz[i] + 2.0 * g_zz_y_xz_zz[i] * a_exp;
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_0_y_yy_xx, g_0_y_yy_xy, g_0_y_yy_xz, g_0_y_yy_yy, g_0_y_yy_yz, g_0_y_yy_zz, g_z_0_0_0_z_y_yy_xx, g_z_0_0_0_z_y_yy_xy, g_z_0_0_0_z_y_yy_xz, g_z_0_0_0_z_y_yy_yy, g_z_0_0_0_z_y_yy_yz, g_z_0_0_0_z_y_yy_zz, g_zz_y_yy_xx, g_zz_y_yy_xy, g_zz_y_yy_xz, g_zz_y_yy_yy, g_zz_y_yy_yz, g_zz_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_y_yy_xx[i] = -g_0_y_yy_xx[i] + 2.0 * g_zz_y_yy_xx[i] * a_exp;

        g_z_0_0_0_z_y_yy_xy[i] = -g_0_y_yy_xy[i] + 2.0 * g_zz_y_yy_xy[i] * a_exp;

        g_z_0_0_0_z_y_yy_xz[i] = -g_0_y_yy_xz[i] + 2.0 * g_zz_y_yy_xz[i] * a_exp;

        g_z_0_0_0_z_y_yy_yy[i] = -g_0_y_yy_yy[i] + 2.0 * g_zz_y_yy_yy[i] * a_exp;

        g_z_0_0_0_z_y_yy_yz[i] = -g_0_y_yy_yz[i] + 2.0 * g_zz_y_yy_yz[i] * a_exp;

        g_z_0_0_0_z_y_yy_zz[i] = -g_0_y_yy_zz[i] + 2.0 * g_zz_y_yy_zz[i] * a_exp;
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_0_y_yz_xx, g_0_y_yz_xy, g_0_y_yz_xz, g_0_y_yz_yy, g_0_y_yz_yz, g_0_y_yz_zz, g_z_0_0_0_z_y_yz_xx, g_z_0_0_0_z_y_yz_xy, g_z_0_0_0_z_y_yz_xz, g_z_0_0_0_z_y_yz_yy, g_z_0_0_0_z_y_yz_yz, g_z_0_0_0_z_y_yz_zz, g_zz_y_yz_xx, g_zz_y_yz_xy, g_zz_y_yz_xz, g_zz_y_yz_yy, g_zz_y_yz_yz, g_zz_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_y_yz_xx[i] = -g_0_y_yz_xx[i] + 2.0 * g_zz_y_yz_xx[i] * a_exp;

        g_z_0_0_0_z_y_yz_xy[i] = -g_0_y_yz_xy[i] + 2.0 * g_zz_y_yz_xy[i] * a_exp;

        g_z_0_0_0_z_y_yz_xz[i] = -g_0_y_yz_xz[i] + 2.0 * g_zz_y_yz_xz[i] * a_exp;

        g_z_0_0_0_z_y_yz_yy[i] = -g_0_y_yz_yy[i] + 2.0 * g_zz_y_yz_yy[i] * a_exp;

        g_z_0_0_0_z_y_yz_yz[i] = -g_0_y_yz_yz[i] + 2.0 * g_zz_y_yz_yz[i] * a_exp;

        g_z_0_0_0_z_y_yz_zz[i] = -g_0_y_yz_zz[i] + 2.0 * g_zz_y_yz_zz[i] * a_exp;
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_0_y_zz_xx, g_0_y_zz_xy, g_0_y_zz_xz, g_0_y_zz_yy, g_0_y_zz_yz, g_0_y_zz_zz, g_z_0_0_0_z_y_zz_xx, g_z_0_0_0_z_y_zz_xy, g_z_0_0_0_z_y_zz_xz, g_z_0_0_0_z_y_zz_yy, g_z_0_0_0_z_y_zz_yz, g_z_0_0_0_z_y_zz_zz, g_zz_y_zz_xx, g_zz_y_zz_xy, g_zz_y_zz_xz, g_zz_y_zz_yy, g_zz_y_zz_yz, g_zz_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_y_zz_xx[i] = -g_0_y_zz_xx[i] + 2.0 * g_zz_y_zz_xx[i] * a_exp;

        g_z_0_0_0_z_y_zz_xy[i] = -g_0_y_zz_xy[i] + 2.0 * g_zz_y_zz_xy[i] * a_exp;

        g_z_0_0_0_z_y_zz_xz[i] = -g_0_y_zz_xz[i] + 2.0 * g_zz_y_zz_xz[i] * a_exp;

        g_z_0_0_0_z_y_zz_yy[i] = -g_0_y_zz_yy[i] + 2.0 * g_zz_y_zz_yy[i] * a_exp;

        g_z_0_0_0_z_y_zz_yz[i] = -g_0_y_zz_yz[i] + 2.0 * g_zz_y_zz_yz[i] * a_exp;

        g_z_0_0_0_z_y_zz_zz[i] = -g_0_y_zz_zz[i] + 2.0 * g_zz_y_zz_zz[i] * a_exp;
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_0_z_xx_xx, g_0_z_xx_xy, g_0_z_xx_xz, g_0_z_xx_yy, g_0_z_xx_yz, g_0_z_xx_zz, g_z_0_0_0_z_z_xx_xx, g_z_0_0_0_z_z_xx_xy, g_z_0_0_0_z_z_xx_xz, g_z_0_0_0_z_z_xx_yy, g_z_0_0_0_z_z_xx_yz, g_z_0_0_0_z_z_xx_zz, g_zz_z_xx_xx, g_zz_z_xx_xy, g_zz_z_xx_xz, g_zz_z_xx_yy, g_zz_z_xx_yz, g_zz_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_z_xx_xx[i] = -g_0_z_xx_xx[i] + 2.0 * g_zz_z_xx_xx[i] * a_exp;

        g_z_0_0_0_z_z_xx_xy[i] = -g_0_z_xx_xy[i] + 2.0 * g_zz_z_xx_xy[i] * a_exp;

        g_z_0_0_0_z_z_xx_xz[i] = -g_0_z_xx_xz[i] + 2.0 * g_zz_z_xx_xz[i] * a_exp;

        g_z_0_0_0_z_z_xx_yy[i] = -g_0_z_xx_yy[i] + 2.0 * g_zz_z_xx_yy[i] * a_exp;

        g_z_0_0_0_z_z_xx_yz[i] = -g_0_z_xx_yz[i] + 2.0 * g_zz_z_xx_yz[i] * a_exp;

        g_z_0_0_0_z_z_xx_zz[i] = -g_0_z_xx_zz[i] + 2.0 * g_zz_z_xx_zz[i] * a_exp;
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_0_z_xy_xx, g_0_z_xy_xy, g_0_z_xy_xz, g_0_z_xy_yy, g_0_z_xy_yz, g_0_z_xy_zz, g_z_0_0_0_z_z_xy_xx, g_z_0_0_0_z_z_xy_xy, g_z_0_0_0_z_z_xy_xz, g_z_0_0_0_z_z_xy_yy, g_z_0_0_0_z_z_xy_yz, g_z_0_0_0_z_z_xy_zz, g_zz_z_xy_xx, g_zz_z_xy_xy, g_zz_z_xy_xz, g_zz_z_xy_yy, g_zz_z_xy_yz, g_zz_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_z_xy_xx[i] = -g_0_z_xy_xx[i] + 2.0 * g_zz_z_xy_xx[i] * a_exp;

        g_z_0_0_0_z_z_xy_xy[i] = -g_0_z_xy_xy[i] + 2.0 * g_zz_z_xy_xy[i] * a_exp;

        g_z_0_0_0_z_z_xy_xz[i] = -g_0_z_xy_xz[i] + 2.0 * g_zz_z_xy_xz[i] * a_exp;

        g_z_0_0_0_z_z_xy_yy[i] = -g_0_z_xy_yy[i] + 2.0 * g_zz_z_xy_yy[i] * a_exp;

        g_z_0_0_0_z_z_xy_yz[i] = -g_0_z_xy_yz[i] + 2.0 * g_zz_z_xy_yz[i] * a_exp;

        g_z_0_0_0_z_z_xy_zz[i] = -g_0_z_xy_zz[i] + 2.0 * g_zz_z_xy_zz[i] * a_exp;
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_0_z_xz_xx, g_0_z_xz_xy, g_0_z_xz_xz, g_0_z_xz_yy, g_0_z_xz_yz, g_0_z_xz_zz, g_z_0_0_0_z_z_xz_xx, g_z_0_0_0_z_z_xz_xy, g_z_0_0_0_z_z_xz_xz, g_z_0_0_0_z_z_xz_yy, g_z_0_0_0_z_z_xz_yz, g_z_0_0_0_z_z_xz_zz, g_zz_z_xz_xx, g_zz_z_xz_xy, g_zz_z_xz_xz, g_zz_z_xz_yy, g_zz_z_xz_yz, g_zz_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_z_xz_xx[i] = -g_0_z_xz_xx[i] + 2.0 * g_zz_z_xz_xx[i] * a_exp;

        g_z_0_0_0_z_z_xz_xy[i] = -g_0_z_xz_xy[i] + 2.0 * g_zz_z_xz_xy[i] * a_exp;

        g_z_0_0_0_z_z_xz_xz[i] = -g_0_z_xz_xz[i] + 2.0 * g_zz_z_xz_xz[i] * a_exp;

        g_z_0_0_0_z_z_xz_yy[i] = -g_0_z_xz_yy[i] + 2.0 * g_zz_z_xz_yy[i] * a_exp;

        g_z_0_0_0_z_z_xz_yz[i] = -g_0_z_xz_yz[i] + 2.0 * g_zz_z_xz_yz[i] * a_exp;

        g_z_0_0_0_z_z_xz_zz[i] = -g_0_z_xz_zz[i] + 2.0 * g_zz_z_xz_zz[i] * a_exp;
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_0_z_yy_xx, g_0_z_yy_xy, g_0_z_yy_xz, g_0_z_yy_yy, g_0_z_yy_yz, g_0_z_yy_zz, g_z_0_0_0_z_z_yy_xx, g_z_0_0_0_z_z_yy_xy, g_z_0_0_0_z_z_yy_xz, g_z_0_0_0_z_z_yy_yy, g_z_0_0_0_z_z_yy_yz, g_z_0_0_0_z_z_yy_zz, g_zz_z_yy_xx, g_zz_z_yy_xy, g_zz_z_yy_xz, g_zz_z_yy_yy, g_zz_z_yy_yz, g_zz_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_z_yy_xx[i] = -g_0_z_yy_xx[i] + 2.0 * g_zz_z_yy_xx[i] * a_exp;

        g_z_0_0_0_z_z_yy_xy[i] = -g_0_z_yy_xy[i] + 2.0 * g_zz_z_yy_xy[i] * a_exp;

        g_z_0_0_0_z_z_yy_xz[i] = -g_0_z_yy_xz[i] + 2.0 * g_zz_z_yy_xz[i] * a_exp;

        g_z_0_0_0_z_z_yy_yy[i] = -g_0_z_yy_yy[i] + 2.0 * g_zz_z_yy_yy[i] * a_exp;

        g_z_0_0_0_z_z_yy_yz[i] = -g_0_z_yy_yz[i] + 2.0 * g_zz_z_yy_yz[i] * a_exp;

        g_z_0_0_0_z_z_yy_zz[i] = -g_0_z_yy_zz[i] + 2.0 * g_zz_z_yy_zz[i] * a_exp;
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_0_z_yz_xx, g_0_z_yz_xy, g_0_z_yz_xz, g_0_z_yz_yy, g_0_z_yz_yz, g_0_z_yz_zz, g_z_0_0_0_z_z_yz_xx, g_z_0_0_0_z_z_yz_xy, g_z_0_0_0_z_z_yz_xz, g_z_0_0_0_z_z_yz_yy, g_z_0_0_0_z_z_yz_yz, g_z_0_0_0_z_z_yz_zz, g_zz_z_yz_xx, g_zz_z_yz_xy, g_zz_z_yz_xz, g_zz_z_yz_yy, g_zz_z_yz_yz, g_zz_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_z_yz_xx[i] = -g_0_z_yz_xx[i] + 2.0 * g_zz_z_yz_xx[i] * a_exp;

        g_z_0_0_0_z_z_yz_xy[i] = -g_0_z_yz_xy[i] + 2.0 * g_zz_z_yz_xy[i] * a_exp;

        g_z_0_0_0_z_z_yz_xz[i] = -g_0_z_yz_xz[i] + 2.0 * g_zz_z_yz_xz[i] * a_exp;

        g_z_0_0_0_z_z_yz_yy[i] = -g_0_z_yz_yy[i] + 2.0 * g_zz_z_yz_yy[i] * a_exp;

        g_z_0_0_0_z_z_yz_yz[i] = -g_0_z_yz_yz[i] + 2.0 * g_zz_z_yz_yz[i] * a_exp;

        g_z_0_0_0_z_z_yz_zz[i] = -g_0_z_yz_zz[i] + 2.0 * g_zz_z_yz_zz[i] * a_exp;
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_0_z_zz_xx, g_0_z_zz_xy, g_0_z_zz_xz, g_0_z_zz_yy, g_0_z_zz_yz, g_0_z_zz_zz, g_z_0_0_0_z_z_zz_xx, g_z_0_0_0_z_z_zz_xy, g_z_0_0_0_z_z_zz_xz, g_z_0_0_0_z_z_zz_yy, g_z_0_0_0_z_z_zz_yz, g_z_0_0_0_z_z_zz_zz, g_zz_z_zz_xx, g_zz_z_zz_xy, g_zz_z_zz_xz, g_zz_z_zz_yy, g_zz_z_zz_yz, g_zz_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_z_zz_xx[i] = -g_0_z_zz_xx[i] + 2.0 * g_zz_z_zz_xx[i] * a_exp;

        g_z_0_0_0_z_z_zz_xy[i] = -g_0_z_zz_xy[i] + 2.0 * g_zz_z_zz_xy[i] * a_exp;

        g_z_0_0_0_z_z_zz_xz[i] = -g_0_z_zz_xz[i] + 2.0 * g_zz_z_zz_xz[i] * a_exp;

        g_z_0_0_0_z_z_zz_yy[i] = -g_0_z_zz_yy[i] + 2.0 * g_zz_z_zz_yy[i] * a_exp;

        g_z_0_0_0_z_z_zz_yz[i] = -g_0_z_zz_yz[i] + 2.0 * g_zz_z_zz_yz[i] * a_exp;

        g_z_0_0_0_z_z_zz_zz[i] = -g_0_z_zz_zz[i] + 2.0 * g_zz_z_zz_zz[i] * a_exp;
    }
}

} // t4c_geom namespace

