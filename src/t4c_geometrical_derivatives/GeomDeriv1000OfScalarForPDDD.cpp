#include "GeomDeriv1000OfScalarForPDDD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_pddd_0(CSimdArray<double>& buffer_1000_pddd,
                     const CSimdArray<double>& buffer_sddd,
                     const CSimdArray<double>& buffer_dddd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_pddd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sddd

    auto g_0_xx_xx_xx = buffer_sddd[0];

    auto g_0_xx_xx_xy = buffer_sddd[1];

    auto g_0_xx_xx_xz = buffer_sddd[2];

    auto g_0_xx_xx_yy = buffer_sddd[3];

    auto g_0_xx_xx_yz = buffer_sddd[4];

    auto g_0_xx_xx_zz = buffer_sddd[5];

    auto g_0_xx_xy_xx = buffer_sddd[6];

    auto g_0_xx_xy_xy = buffer_sddd[7];

    auto g_0_xx_xy_xz = buffer_sddd[8];

    auto g_0_xx_xy_yy = buffer_sddd[9];

    auto g_0_xx_xy_yz = buffer_sddd[10];

    auto g_0_xx_xy_zz = buffer_sddd[11];

    auto g_0_xx_xz_xx = buffer_sddd[12];

    auto g_0_xx_xz_xy = buffer_sddd[13];

    auto g_0_xx_xz_xz = buffer_sddd[14];

    auto g_0_xx_xz_yy = buffer_sddd[15];

    auto g_0_xx_xz_yz = buffer_sddd[16];

    auto g_0_xx_xz_zz = buffer_sddd[17];

    auto g_0_xx_yy_xx = buffer_sddd[18];

    auto g_0_xx_yy_xy = buffer_sddd[19];

    auto g_0_xx_yy_xz = buffer_sddd[20];

    auto g_0_xx_yy_yy = buffer_sddd[21];

    auto g_0_xx_yy_yz = buffer_sddd[22];

    auto g_0_xx_yy_zz = buffer_sddd[23];

    auto g_0_xx_yz_xx = buffer_sddd[24];

    auto g_0_xx_yz_xy = buffer_sddd[25];

    auto g_0_xx_yz_xz = buffer_sddd[26];

    auto g_0_xx_yz_yy = buffer_sddd[27];

    auto g_0_xx_yz_yz = buffer_sddd[28];

    auto g_0_xx_yz_zz = buffer_sddd[29];

    auto g_0_xx_zz_xx = buffer_sddd[30];

    auto g_0_xx_zz_xy = buffer_sddd[31];

    auto g_0_xx_zz_xz = buffer_sddd[32];

    auto g_0_xx_zz_yy = buffer_sddd[33];

    auto g_0_xx_zz_yz = buffer_sddd[34];

    auto g_0_xx_zz_zz = buffer_sddd[35];

    auto g_0_xy_xx_xx = buffer_sddd[36];

    auto g_0_xy_xx_xy = buffer_sddd[37];

    auto g_0_xy_xx_xz = buffer_sddd[38];

    auto g_0_xy_xx_yy = buffer_sddd[39];

    auto g_0_xy_xx_yz = buffer_sddd[40];

    auto g_0_xy_xx_zz = buffer_sddd[41];

    auto g_0_xy_xy_xx = buffer_sddd[42];

    auto g_0_xy_xy_xy = buffer_sddd[43];

    auto g_0_xy_xy_xz = buffer_sddd[44];

    auto g_0_xy_xy_yy = buffer_sddd[45];

    auto g_0_xy_xy_yz = buffer_sddd[46];

    auto g_0_xy_xy_zz = buffer_sddd[47];

    auto g_0_xy_xz_xx = buffer_sddd[48];

    auto g_0_xy_xz_xy = buffer_sddd[49];

    auto g_0_xy_xz_xz = buffer_sddd[50];

    auto g_0_xy_xz_yy = buffer_sddd[51];

    auto g_0_xy_xz_yz = buffer_sddd[52];

    auto g_0_xy_xz_zz = buffer_sddd[53];

    auto g_0_xy_yy_xx = buffer_sddd[54];

    auto g_0_xy_yy_xy = buffer_sddd[55];

    auto g_0_xy_yy_xz = buffer_sddd[56];

    auto g_0_xy_yy_yy = buffer_sddd[57];

    auto g_0_xy_yy_yz = buffer_sddd[58];

    auto g_0_xy_yy_zz = buffer_sddd[59];

    auto g_0_xy_yz_xx = buffer_sddd[60];

    auto g_0_xy_yz_xy = buffer_sddd[61];

    auto g_0_xy_yz_xz = buffer_sddd[62];

    auto g_0_xy_yz_yy = buffer_sddd[63];

    auto g_0_xy_yz_yz = buffer_sddd[64];

    auto g_0_xy_yz_zz = buffer_sddd[65];

    auto g_0_xy_zz_xx = buffer_sddd[66];

    auto g_0_xy_zz_xy = buffer_sddd[67];

    auto g_0_xy_zz_xz = buffer_sddd[68];

    auto g_0_xy_zz_yy = buffer_sddd[69];

    auto g_0_xy_zz_yz = buffer_sddd[70];

    auto g_0_xy_zz_zz = buffer_sddd[71];

    auto g_0_xz_xx_xx = buffer_sddd[72];

    auto g_0_xz_xx_xy = buffer_sddd[73];

    auto g_0_xz_xx_xz = buffer_sddd[74];

    auto g_0_xz_xx_yy = buffer_sddd[75];

    auto g_0_xz_xx_yz = buffer_sddd[76];

    auto g_0_xz_xx_zz = buffer_sddd[77];

    auto g_0_xz_xy_xx = buffer_sddd[78];

    auto g_0_xz_xy_xy = buffer_sddd[79];

    auto g_0_xz_xy_xz = buffer_sddd[80];

    auto g_0_xz_xy_yy = buffer_sddd[81];

    auto g_0_xz_xy_yz = buffer_sddd[82];

    auto g_0_xz_xy_zz = buffer_sddd[83];

    auto g_0_xz_xz_xx = buffer_sddd[84];

    auto g_0_xz_xz_xy = buffer_sddd[85];

    auto g_0_xz_xz_xz = buffer_sddd[86];

    auto g_0_xz_xz_yy = buffer_sddd[87];

    auto g_0_xz_xz_yz = buffer_sddd[88];

    auto g_0_xz_xz_zz = buffer_sddd[89];

    auto g_0_xz_yy_xx = buffer_sddd[90];

    auto g_0_xz_yy_xy = buffer_sddd[91];

    auto g_0_xz_yy_xz = buffer_sddd[92];

    auto g_0_xz_yy_yy = buffer_sddd[93];

    auto g_0_xz_yy_yz = buffer_sddd[94];

    auto g_0_xz_yy_zz = buffer_sddd[95];

    auto g_0_xz_yz_xx = buffer_sddd[96];

    auto g_0_xz_yz_xy = buffer_sddd[97];

    auto g_0_xz_yz_xz = buffer_sddd[98];

    auto g_0_xz_yz_yy = buffer_sddd[99];

    auto g_0_xz_yz_yz = buffer_sddd[100];

    auto g_0_xz_yz_zz = buffer_sddd[101];

    auto g_0_xz_zz_xx = buffer_sddd[102];

    auto g_0_xz_zz_xy = buffer_sddd[103];

    auto g_0_xz_zz_xz = buffer_sddd[104];

    auto g_0_xz_zz_yy = buffer_sddd[105];

    auto g_0_xz_zz_yz = buffer_sddd[106];

    auto g_0_xz_zz_zz = buffer_sddd[107];

    auto g_0_yy_xx_xx = buffer_sddd[108];

    auto g_0_yy_xx_xy = buffer_sddd[109];

    auto g_0_yy_xx_xz = buffer_sddd[110];

    auto g_0_yy_xx_yy = buffer_sddd[111];

    auto g_0_yy_xx_yz = buffer_sddd[112];

    auto g_0_yy_xx_zz = buffer_sddd[113];

    auto g_0_yy_xy_xx = buffer_sddd[114];

    auto g_0_yy_xy_xy = buffer_sddd[115];

    auto g_0_yy_xy_xz = buffer_sddd[116];

    auto g_0_yy_xy_yy = buffer_sddd[117];

    auto g_0_yy_xy_yz = buffer_sddd[118];

    auto g_0_yy_xy_zz = buffer_sddd[119];

    auto g_0_yy_xz_xx = buffer_sddd[120];

    auto g_0_yy_xz_xy = buffer_sddd[121];

    auto g_0_yy_xz_xz = buffer_sddd[122];

    auto g_0_yy_xz_yy = buffer_sddd[123];

    auto g_0_yy_xz_yz = buffer_sddd[124];

    auto g_0_yy_xz_zz = buffer_sddd[125];

    auto g_0_yy_yy_xx = buffer_sddd[126];

    auto g_0_yy_yy_xy = buffer_sddd[127];

    auto g_0_yy_yy_xz = buffer_sddd[128];

    auto g_0_yy_yy_yy = buffer_sddd[129];

    auto g_0_yy_yy_yz = buffer_sddd[130];

    auto g_0_yy_yy_zz = buffer_sddd[131];

    auto g_0_yy_yz_xx = buffer_sddd[132];

    auto g_0_yy_yz_xy = buffer_sddd[133];

    auto g_0_yy_yz_xz = buffer_sddd[134];

    auto g_0_yy_yz_yy = buffer_sddd[135];

    auto g_0_yy_yz_yz = buffer_sddd[136];

    auto g_0_yy_yz_zz = buffer_sddd[137];

    auto g_0_yy_zz_xx = buffer_sddd[138];

    auto g_0_yy_zz_xy = buffer_sddd[139];

    auto g_0_yy_zz_xz = buffer_sddd[140];

    auto g_0_yy_zz_yy = buffer_sddd[141];

    auto g_0_yy_zz_yz = buffer_sddd[142];

    auto g_0_yy_zz_zz = buffer_sddd[143];

    auto g_0_yz_xx_xx = buffer_sddd[144];

    auto g_0_yz_xx_xy = buffer_sddd[145];

    auto g_0_yz_xx_xz = buffer_sddd[146];

    auto g_0_yz_xx_yy = buffer_sddd[147];

    auto g_0_yz_xx_yz = buffer_sddd[148];

    auto g_0_yz_xx_zz = buffer_sddd[149];

    auto g_0_yz_xy_xx = buffer_sddd[150];

    auto g_0_yz_xy_xy = buffer_sddd[151];

    auto g_0_yz_xy_xz = buffer_sddd[152];

    auto g_0_yz_xy_yy = buffer_sddd[153];

    auto g_0_yz_xy_yz = buffer_sddd[154];

    auto g_0_yz_xy_zz = buffer_sddd[155];

    auto g_0_yz_xz_xx = buffer_sddd[156];

    auto g_0_yz_xz_xy = buffer_sddd[157];

    auto g_0_yz_xz_xz = buffer_sddd[158];

    auto g_0_yz_xz_yy = buffer_sddd[159];

    auto g_0_yz_xz_yz = buffer_sddd[160];

    auto g_0_yz_xz_zz = buffer_sddd[161];

    auto g_0_yz_yy_xx = buffer_sddd[162];

    auto g_0_yz_yy_xy = buffer_sddd[163];

    auto g_0_yz_yy_xz = buffer_sddd[164];

    auto g_0_yz_yy_yy = buffer_sddd[165];

    auto g_0_yz_yy_yz = buffer_sddd[166];

    auto g_0_yz_yy_zz = buffer_sddd[167];

    auto g_0_yz_yz_xx = buffer_sddd[168];

    auto g_0_yz_yz_xy = buffer_sddd[169];

    auto g_0_yz_yz_xz = buffer_sddd[170];

    auto g_0_yz_yz_yy = buffer_sddd[171];

    auto g_0_yz_yz_yz = buffer_sddd[172];

    auto g_0_yz_yz_zz = buffer_sddd[173];

    auto g_0_yz_zz_xx = buffer_sddd[174];

    auto g_0_yz_zz_xy = buffer_sddd[175];

    auto g_0_yz_zz_xz = buffer_sddd[176];

    auto g_0_yz_zz_yy = buffer_sddd[177];

    auto g_0_yz_zz_yz = buffer_sddd[178];

    auto g_0_yz_zz_zz = buffer_sddd[179];

    auto g_0_zz_xx_xx = buffer_sddd[180];

    auto g_0_zz_xx_xy = buffer_sddd[181];

    auto g_0_zz_xx_xz = buffer_sddd[182];

    auto g_0_zz_xx_yy = buffer_sddd[183];

    auto g_0_zz_xx_yz = buffer_sddd[184];

    auto g_0_zz_xx_zz = buffer_sddd[185];

    auto g_0_zz_xy_xx = buffer_sddd[186];

    auto g_0_zz_xy_xy = buffer_sddd[187];

    auto g_0_zz_xy_xz = buffer_sddd[188];

    auto g_0_zz_xy_yy = buffer_sddd[189];

    auto g_0_zz_xy_yz = buffer_sddd[190];

    auto g_0_zz_xy_zz = buffer_sddd[191];

    auto g_0_zz_xz_xx = buffer_sddd[192];

    auto g_0_zz_xz_xy = buffer_sddd[193];

    auto g_0_zz_xz_xz = buffer_sddd[194];

    auto g_0_zz_xz_yy = buffer_sddd[195];

    auto g_0_zz_xz_yz = buffer_sddd[196];

    auto g_0_zz_xz_zz = buffer_sddd[197];

    auto g_0_zz_yy_xx = buffer_sddd[198];

    auto g_0_zz_yy_xy = buffer_sddd[199];

    auto g_0_zz_yy_xz = buffer_sddd[200];

    auto g_0_zz_yy_yy = buffer_sddd[201];

    auto g_0_zz_yy_yz = buffer_sddd[202];

    auto g_0_zz_yy_zz = buffer_sddd[203];

    auto g_0_zz_yz_xx = buffer_sddd[204];

    auto g_0_zz_yz_xy = buffer_sddd[205];

    auto g_0_zz_yz_xz = buffer_sddd[206];

    auto g_0_zz_yz_yy = buffer_sddd[207];

    auto g_0_zz_yz_yz = buffer_sddd[208];

    auto g_0_zz_yz_zz = buffer_sddd[209];

    auto g_0_zz_zz_xx = buffer_sddd[210];

    auto g_0_zz_zz_xy = buffer_sddd[211];

    auto g_0_zz_zz_xz = buffer_sddd[212];

    auto g_0_zz_zz_yy = buffer_sddd[213];

    auto g_0_zz_zz_yz = buffer_sddd[214];

    auto g_0_zz_zz_zz = buffer_sddd[215];

    /// Set up components of auxilary buffer : buffer_dddd

    auto g_xx_xx_xx_xx = buffer_dddd[0];

    auto g_xx_xx_xx_xy = buffer_dddd[1];

    auto g_xx_xx_xx_xz = buffer_dddd[2];

    auto g_xx_xx_xx_yy = buffer_dddd[3];

    auto g_xx_xx_xx_yz = buffer_dddd[4];

    auto g_xx_xx_xx_zz = buffer_dddd[5];

    auto g_xx_xx_xy_xx = buffer_dddd[6];

    auto g_xx_xx_xy_xy = buffer_dddd[7];

    auto g_xx_xx_xy_xz = buffer_dddd[8];

    auto g_xx_xx_xy_yy = buffer_dddd[9];

    auto g_xx_xx_xy_yz = buffer_dddd[10];

    auto g_xx_xx_xy_zz = buffer_dddd[11];

    auto g_xx_xx_xz_xx = buffer_dddd[12];

    auto g_xx_xx_xz_xy = buffer_dddd[13];

    auto g_xx_xx_xz_xz = buffer_dddd[14];

    auto g_xx_xx_xz_yy = buffer_dddd[15];

    auto g_xx_xx_xz_yz = buffer_dddd[16];

    auto g_xx_xx_xz_zz = buffer_dddd[17];

    auto g_xx_xx_yy_xx = buffer_dddd[18];

    auto g_xx_xx_yy_xy = buffer_dddd[19];

    auto g_xx_xx_yy_xz = buffer_dddd[20];

    auto g_xx_xx_yy_yy = buffer_dddd[21];

    auto g_xx_xx_yy_yz = buffer_dddd[22];

    auto g_xx_xx_yy_zz = buffer_dddd[23];

    auto g_xx_xx_yz_xx = buffer_dddd[24];

    auto g_xx_xx_yz_xy = buffer_dddd[25];

    auto g_xx_xx_yz_xz = buffer_dddd[26];

    auto g_xx_xx_yz_yy = buffer_dddd[27];

    auto g_xx_xx_yz_yz = buffer_dddd[28];

    auto g_xx_xx_yz_zz = buffer_dddd[29];

    auto g_xx_xx_zz_xx = buffer_dddd[30];

    auto g_xx_xx_zz_xy = buffer_dddd[31];

    auto g_xx_xx_zz_xz = buffer_dddd[32];

    auto g_xx_xx_zz_yy = buffer_dddd[33];

    auto g_xx_xx_zz_yz = buffer_dddd[34];

    auto g_xx_xx_zz_zz = buffer_dddd[35];

    auto g_xx_xy_xx_xx = buffer_dddd[36];

    auto g_xx_xy_xx_xy = buffer_dddd[37];

    auto g_xx_xy_xx_xz = buffer_dddd[38];

    auto g_xx_xy_xx_yy = buffer_dddd[39];

    auto g_xx_xy_xx_yz = buffer_dddd[40];

    auto g_xx_xy_xx_zz = buffer_dddd[41];

    auto g_xx_xy_xy_xx = buffer_dddd[42];

    auto g_xx_xy_xy_xy = buffer_dddd[43];

    auto g_xx_xy_xy_xz = buffer_dddd[44];

    auto g_xx_xy_xy_yy = buffer_dddd[45];

    auto g_xx_xy_xy_yz = buffer_dddd[46];

    auto g_xx_xy_xy_zz = buffer_dddd[47];

    auto g_xx_xy_xz_xx = buffer_dddd[48];

    auto g_xx_xy_xz_xy = buffer_dddd[49];

    auto g_xx_xy_xz_xz = buffer_dddd[50];

    auto g_xx_xy_xz_yy = buffer_dddd[51];

    auto g_xx_xy_xz_yz = buffer_dddd[52];

    auto g_xx_xy_xz_zz = buffer_dddd[53];

    auto g_xx_xy_yy_xx = buffer_dddd[54];

    auto g_xx_xy_yy_xy = buffer_dddd[55];

    auto g_xx_xy_yy_xz = buffer_dddd[56];

    auto g_xx_xy_yy_yy = buffer_dddd[57];

    auto g_xx_xy_yy_yz = buffer_dddd[58];

    auto g_xx_xy_yy_zz = buffer_dddd[59];

    auto g_xx_xy_yz_xx = buffer_dddd[60];

    auto g_xx_xy_yz_xy = buffer_dddd[61];

    auto g_xx_xy_yz_xz = buffer_dddd[62];

    auto g_xx_xy_yz_yy = buffer_dddd[63];

    auto g_xx_xy_yz_yz = buffer_dddd[64];

    auto g_xx_xy_yz_zz = buffer_dddd[65];

    auto g_xx_xy_zz_xx = buffer_dddd[66];

    auto g_xx_xy_zz_xy = buffer_dddd[67];

    auto g_xx_xy_zz_xz = buffer_dddd[68];

    auto g_xx_xy_zz_yy = buffer_dddd[69];

    auto g_xx_xy_zz_yz = buffer_dddd[70];

    auto g_xx_xy_zz_zz = buffer_dddd[71];

    auto g_xx_xz_xx_xx = buffer_dddd[72];

    auto g_xx_xz_xx_xy = buffer_dddd[73];

    auto g_xx_xz_xx_xz = buffer_dddd[74];

    auto g_xx_xz_xx_yy = buffer_dddd[75];

    auto g_xx_xz_xx_yz = buffer_dddd[76];

    auto g_xx_xz_xx_zz = buffer_dddd[77];

    auto g_xx_xz_xy_xx = buffer_dddd[78];

    auto g_xx_xz_xy_xy = buffer_dddd[79];

    auto g_xx_xz_xy_xz = buffer_dddd[80];

    auto g_xx_xz_xy_yy = buffer_dddd[81];

    auto g_xx_xz_xy_yz = buffer_dddd[82];

    auto g_xx_xz_xy_zz = buffer_dddd[83];

    auto g_xx_xz_xz_xx = buffer_dddd[84];

    auto g_xx_xz_xz_xy = buffer_dddd[85];

    auto g_xx_xz_xz_xz = buffer_dddd[86];

    auto g_xx_xz_xz_yy = buffer_dddd[87];

    auto g_xx_xz_xz_yz = buffer_dddd[88];

    auto g_xx_xz_xz_zz = buffer_dddd[89];

    auto g_xx_xz_yy_xx = buffer_dddd[90];

    auto g_xx_xz_yy_xy = buffer_dddd[91];

    auto g_xx_xz_yy_xz = buffer_dddd[92];

    auto g_xx_xz_yy_yy = buffer_dddd[93];

    auto g_xx_xz_yy_yz = buffer_dddd[94];

    auto g_xx_xz_yy_zz = buffer_dddd[95];

    auto g_xx_xz_yz_xx = buffer_dddd[96];

    auto g_xx_xz_yz_xy = buffer_dddd[97];

    auto g_xx_xz_yz_xz = buffer_dddd[98];

    auto g_xx_xz_yz_yy = buffer_dddd[99];

    auto g_xx_xz_yz_yz = buffer_dddd[100];

    auto g_xx_xz_yz_zz = buffer_dddd[101];

    auto g_xx_xz_zz_xx = buffer_dddd[102];

    auto g_xx_xz_zz_xy = buffer_dddd[103];

    auto g_xx_xz_zz_xz = buffer_dddd[104];

    auto g_xx_xz_zz_yy = buffer_dddd[105];

    auto g_xx_xz_zz_yz = buffer_dddd[106];

    auto g_xx_xz_zz_zz = buffer_dddd[107];

    auto g_xx_yy_xx_xx = buffer_dddd[108];

    auto g_xx_yy_xx_xy = buffer_dddd[109];

    auto g_xx_yy_xx_xz = buffer_dddd[110];

    auto g_xx_yy_xx_yy = buffer_dddd[111];

    auto g_xx_yy_xx_yz = buffer_dddd[112];

    auto g_xx_yy_xx_zz = buffer_dddd[113];

    auto g_xx_yy_xy_xx = buffer_dddd[114];

    auto g_xx_yy_xy_xy = buffer_dddd[115];

    auto g_xx_yy_xy_xz = buffer_dddd[116];

    auto g_xx_yy_xy_yy = buffer_dddd[117];

    auto g_xx_yy_xy_yz = buffer_dddd[118];

    auto g_xx_yy_xy_zz = buffer_dddd[119];

    auto g_xx_yy_xz_xx = buffer_dddd[120];

    auto g_xx_yy_xz_xy = buffer_dddd[121];

    auto g_xx_yy_xz_xz = buffer_dddd[122];

    auto g_xx_yy_xz_yy = buffer_dddd[123];

    auto g_xx_yy_xz_yz = buffer_dddd[124];

    auto g_xx_yy_xz_zz = buffer_dddd[125];

    auto g_xx_yy_yy_xx = buffer_dddd[126];

    auto g_xx_yy_yy_xy = buffer_dddd[127];

    auto g_xx_yy_yy_xz = buffer_dddd[128];

    auto g_xx_yy_yy_yy = buffer_dddd[129];

    auto g_xx_yy_yy_yz = buffer_dddd[130];

    auto g_xx_yy_yy_zz = buffer_dddd[131];

    auto g_xx_yy_yz_xx = buffer_dddd[132];

    auto g_xx_yy_yz_xy = buffer_dddd[133];

    auto g_xx_yy_yz_xz = buffer_dddd[134];

    auto g_xx_yy_yz_yy = buffer_dddd[135];

    auto g_xx_yy_yz_yz = buffer_dddd[136];

    auto g_xx_yy_yz_zz = buffer_dddd[137];

    auto g_xx_yy_zz_xx = buffer_dddd[138];

    auto g_xx_yy_zz_xy = buffer_dddd[139];

    auto g_xx_yy_zz_xz = buffer_dddd[140];

    auto g_xx_yy_zz_yy = buffer_dddd[141];

    auto g_xx_yy_zz_yz = buffer_dddd[142];

    auto g_xx_yy_zz_zz = buffer_dddd[143];

    auto g_xx_yz_xx_xx = buffer_dddd[144];

    auto g_xx_yz_xx_xy = buffer_dddd[145];

    auto g_xx_yz_xx_xz = buffer_dddd[146];

    auto g_xx_yz_xx_yy = buffer_dddd[147];

    auto g_xx_yz_xx_yz = buffer_dddd[148];

    auto g_xx_yz_xx_zz = buffer_dddd[149];

    auto g_xx_yz_xy_xx = buffer_dddd[150];

    auto g_xx_yz_xy_xy = buffer_dddd[151];

    auto g_xx_yz_xy_xz = buffer_dddd[152];

    auto g_xx_yz_xy_yy = buffer_dddd[153];

    auto g_xx_yz_xy_yz = buffer_dddd[154];

    auto g_xx_yz_xy_zz = buffer_dddd[155];

    auto g_xx_yz_xz_xx = buffer_dddd[156];

    auto g_xx_yz_xz_xy = buffer_dddd[157];

    auto g_xx_yz_xz_xz = buffer_dddd[158];

    auto g_xx_yz_xz_yy = buffer_dddd[159];

    auto g_xx_yz_xz_yz = buffer_dddd[160];

    auto g_xx_yz_xz_zz = buffer_dddd[161];

    auto g_xx_yz_yy_xx = buffer_dddd[162];

    auto g_xx_yz_yy_xy = buffer_dddd[163];

    auto g_xx_yz_yy_xz = buffer_dddd[164];

    auto g_xx_yz_yy_yy = buffer_dddd[165];

    auto g_xx_yz_yy_yz = buffer_dddd[166];

    auto g_xx_yz_yy_zz = buffer_dddd[167];

    auto g_xx_yz_yz_xx = buffer_dddd[168];

    auto g_xx_yz_yz_xy = buffer_dddd[169];

    auto g_xx_yz_yz_xz = buffer_dddd[170];

    auto g_xx_yz_yz_yy = buffer_dddd[171];

    auto g_xx_yz_yz_yz = buffer_dddd[172];

    auto g_xx_yz_yz_zz = buffer_dddd[173];

    auto g_xx_yz_zz_xx = buffer_dddd[174];

    auto g_xx_yz_zz_xy = buffer_dddd[175];

    auto g_xx_yz_zz_xz = buffer_dddd[176];

    auto g_xx_yz_zz_yy = buffer_dddd[177];

    auto g_xx_yz_zz_yz = buffer_dddd[178];

    auto g_xx_yz_zz_zz = buffer_dddd[179];

    auto g_xx_zz_xx_xx = buffer_dddd[180];

    auto g_xx_zz_xx_xy = buffer_dddd[181];

    auto g_xx_zz_xx_xz = buffer_dddd[182];

    auto g_xx_zz_xx_yy = buffer_dddd[183];

    auto g_xx_zz_xx_yz = buffer_dddd[184];

    auto g_xx_zz_xx_zz = buffer_dddd[185];

    auto g_xx_zz_xy_xx = buffer_dddd[186];

    auto g_xx_zz_xy_xy = buffer_dddd[187];

    auto g_xx_zz_xy_xz = buffer_dddd[188];

    auto g_xx_zz_xy_yy = buffer_dddd[189];

    auto g_xx_zz_xy_yz = buffer_dddd[190];

    auto g_xx_zz_xy_zz = buffer_dddd[191];

    auto g_xx_zz_xz_xx = buffer_dddd[192];

    auto g_xx_zz_xz_xy = buffer_dddd[193];

    auto g_xx_zz_xz_xz = buffer_dddd[194];

    auto g_xx_zz_xz_yy = buffer_dddd[195];

    auto g_xx_zz_xz_yz = buffer_dddd[196];

    auto g_xx_zz_xz_zz = buffer_dddd[197];

    auto g_xx_zz_yy_xx = buffer_dddd[198];

    auto g_xx_zz_yy_xy = buffer_dddd[199];

    auto g_xx_zz_yy_xz = buffer_dddd[200];

    auto g_xx_zz_yy_yy = buffer_dddd[201];

    auto g_xx_zz_yy_yz = buffer_dddd[202];

    auto g_xx_zz_yy_zz = buffer_dddd[203];

    auto g_xx_zz_yz_xx = buffer_dddd[204];

    auto g_xx_zz_yz_xy = buffer_dddd[205];

    auto g_xx_zz_yz_xz = buffer_dddd[206];

    auto g_xx_zz_yz_yy = buffer_dddd[207];

    auto g_xx_zz_yz_yz = buffer_dddd[208];

    auto g_xx_zz_yz_zz = buffer_dddd[209];

    auto g_xx_zz_zz_xx = buffer_dddd[210];

    auto g_xx_zz_zz_xy = buffer_dddd[211];

    auto g_xx_zz_zz_xz = buffer_dddd[212];

    auto g_xx_zz_zz_yy = buffer_dddd[213];

    auto g_xx_zz_zz_yz = buffer_dddd[214];

    auto g_xx_zz_zz_zz = buffer_dddd[215];

    auto g_xy_xx_xx_xx = buffer_dddd[216];

    auto g_xy_xx_xx_xy = buffer_dddd[217];

    auto g_xy_xx_xx_xz = buffer_dddd[218];

    auto g_xy_xx_xx_yy = buffer_dddd[219];

    auto g_xy_xx_xx_yz = buffer_dddd[220];

    auto g_xy_xx_xx_zz = buffer_dddd[221];

    auto g_xy_xx_xy_xx = buffer_dddd[222];

    auto g_xy_xx_xy_xy = buffer_dddd[223];

    auto g_xy_xx_xy_xz = buffer_dddd[224];

    auto g_xy_xx_xy_yy = buffer_dddd[225];

    auto g_xy_xx_xy_yz = buffer_dddd[226];

    auto g_xy_xx_xy_zz = buffer_dddd[227];

    auto g_xy_xx_xz_xx = buffer_dddd[228];

    auto g_xy_xx_xz_xy = buffer_dddd[229];

    auto g_xy_xx_xz_xz = buffer_dddd[230];

    auto g_xy_xx_xz_yy = buffer_dddd[231];

    auto g_xy_xx_xz_yz = buffer_dddd[232];

    auto g_xy_xx_xz_zz = buffer_dddd[233];

    auto g_xy_xx_yy_xx = buffer_dddd[234];

    auto g_xy_xx_yy_xy = buffer_dddd[235];

    auto g_xy_xx_yy_xz = buffer_dddd[236];

    auto g_xy_xx_yy_yy = buffer_dddd[237];

    auto g_xy_xx_yy_yz = buffer_dddd[238];

    auto g_xy_xx_yy_zz = buffer_dddd[239];

    auto g_xy_xx_yz_xx = buffer_dddd[240];

    auto g_xy_xx_yz_xy = buffer_dddd[241];

    auto g_xy_xx_yz_xz = buffer_dddd[242];

    auto g_xy_xx_yz_yy = buffer_dddd[243];

    auto g_xy_xx_yz_yz = buffer_dddd[244];

    auto g_xy_xx_yz_zz = buffer_dddd[245];

    auto g_xy_xx_zz_xx = buffer_dddd[246];

    auto g_xy_xx_zz_xy = buffer_dddd[247];

    auto g_xy_xx_zz_xz = buffer_dddd[248];

    auto g_xy_xx_zz_yy = buffer_dddd[249];

    auto g_xy_xx_zz_yz = buffer_dddd[250];

    auto g_xy_xx_zz_zz = buffer_dddd[251];

    auto g_xy_xy_xx_xx = buffer_dddd[252];

    auto g_xy_xy_xx_xy = buffer_dddd[253];

    auto g_xy_xy_xx_xz = buffer_dddd[254];

    auto g_xy_xy_xx_yy = buffer_dddd[255];

    auto g_xy_xy_xx_yz = buffer_dddd[256];

    auto g_xy_xy_xx_zz = buffer_dddd[257];

    auto g_xy_xy_xy_xx = buffer_dddd[258];

    auto g_xy_xy_xy_xy = buffer_dddd[259];

    auto g_xy_xy_xy_xz = buffer_dddd[260];

    auto g_xy_xy_xy_yy = buffer_dddd[261];

    auto g_xy_xy_xy_yz = buffer_dddd[262];

    auto g_xy_xy_xy_zz = buffer_dddd[263];

    auto g_xy_xy_xz_xx = buffer_dddd[264];

    auto g_xy_xy_xz_xy = buffer_dddd[265];

    auto g_xy_xy_xz_xz = buffer_dddd[266];

    auto g_xy_xy_xz_yy = buffer_dddd[267];

    auto g_xy_xy_xz_yz = buffer_dddd[268];

    auto g_xy_xy_xz_zz = buffer_dddd[269];

    auto g_xy_xy_yy_xx = buffer_dddd[270];

    auto g_xy_xy_yy_xy = buffer_dddd[271];

    auto g_xy_xy_yy_xz = buffer_dddd[272];

    auto g_xy_xy_yy_yy = buffer_dddd[273];

    auto g_xy_xy_yy_yz = buffer_dddd[274];

    auto g_xy_xy_yy_zz = buffer_dddd[275];

    auto g_xy_xy_yz_xx = buffer_dddd[276];

    auto g_xy_xy_yz_xy = buffer_dddd[277];

    auto g_xy_xy_yz_xz = buffer_dddd[278];

    auto g_xy_xy_yz_yy = buffer_dddd[279];

    auto g_xy_xy_yz_yz = buffer_dddd[280];

    auto g_xy_xy_yz_zz = buffer_dddd[281];

    auto g_xy_xy_zz_xx = buffer_dddd[282];

    auto g_xy_xy_zz_xy = buffer_dddd[283];

    auto g_xy_xy_zz_xz = buffer_dddd[284];

    auto g_xy_xy_zz_yy = buffer_dddd[285];

    auto g_xy_xy_zz_yz = buffer_dddd[286];

    auto g_xy_xy_zz_zz = buffer_dddd[287];

    auto g_xy_xz_xx_xx = buffer_dddd[288];

    auto g_xy_xz_xx_xy = buffer_dddd[289];

    auto g_xy_xz_xx_xz = buffer_dddd[290];

    auto g_xy_xz_xx_yy = buffer_dddd[291];

    auto g_xy_xz_xx_yz = buffer_dddd[292];

    auto g_xy_xz_xx_zz = buffer_dddd[293];

    auto g_xy_xz_xy_xx = buffer_dddd[294];

    auto g_xy_xz_xy_xy = buffer_dddd[295];

    auto g_xy_xz_xy_xz = buffer_dddd[296];

    auto g_xy_xz_xy_yy = buffer_dddd[297];

    auto g_xy_xz_xy_yz = buffer_dddd[298];

    auto g_xy_xz_xy_zz = buffer_dddd[299];

    auto g_xy_xz_xz_xx = buffer_dddd[300];

    auto g_xy_xz_xz_xy = buffer_dddd[301];

    auto g_xy_xz_xz_xz = buffer_dddd[302];

    auto g_xy_xz_xz_yy = buffer_dddd[303];

    auto g_xy_xz_xz_yz = buffer_dddd[304];

    auto g_xy_xz_xz_zz = buffer_dddd[305];

    auto g_xy_xz_yy_xx = buffer_dddd[306];

    auto g_xy_xz_yy_xy = buffer_dddd[307];

    auto g_xy_xz_yy_xz = buffer_dddd[308];

    auto g_xy_xz_yy_yy = buffer_dddd[309];

    auto g_xy_xz_yy_yz = buffer_dddd[310];

    auto g_xy_xz_yy_zz = buffer_dddd[311];

    auto g_xy_xz_yz_xx = buffer_dddd[312];

    auto g_xy_xz_yz_xy = buffer_dddd[313];

    auto g_xy_xz_yz_xz = buffer_dddd[314];

    auto g_xy_xz_yz_yy = buffer_dddd[315];

    auto g_xy_xz_yz_yz = buffer_dddd[316];

    auto g_xy_xz_yz_zz = buffer_dddd[317];

    auto g_xy_xz_zz_xx = buffer_dddd[318];

    auto g_xy_xz_zz_xy = buffer_dddd[319];

    auto g_xy_xz_zz_xz = buffer_dddd[320];

    auto g_xy_xz_zz_yy = buffer_dddd[321];

    auto g_xy_xz_zz_yz = buffer_dddd[322];

    auto g_xy_xz_zz_zz = buffer_dddd[323];

    auto g_xy_yy_xx_xx = buffer_dddd[324];

    auto g_xy_yy_xx_xy = buffer_dddd[325];

    auto g_xy_yy_xx_xz = buffer_dddd[326];

    auto g_xy_yy_xx_yy = buffer_dddd[327];

    auto g_xy_yy_xx_yz = buffer_dddd[328];

    auto g_xy_yy_xx_zz = buffer_dddd[329];

    auto g_xy_yy_xy_xx = buffer_dddd[330];

    auto g_xy_yy_xy_xy = buffer_dddd[331];

    auto g_xy_yy_xy_xz = buffer_dddd[332];

    auto g_xy_yy_xy_yy = buffer_dddd[333];

    auto g_xy_yy_xy_yz = buffer_dddd[334];

    auto g_xy_yy_xy_zz = buffer_dddd[335];

    auto g_xy_yy_xz_xx = buffer_dddd[336];

    auto g_xy_yy_xz_xy = buffer_dddd[337];

    auto g_xy_yy_xz_xz = buffer_dddd[338];

    auto g_xy_yy_xz_yy = buffer_dddd[339];

    auto g_xy_yy_xz_yz = buffer_dddd[340];

    auto g_xy_yy_xz_zz = buffer_dddd[341];

    auto g_xy_yy_yy_xx = buffer_dddd[342];

    auto g_xy_yy_yy_xy = buffer_dddd[343];

    auto g_xy_yy_yy_xz = buffer_dddd[344];

    auto g_xy_yy_yy_yy = buffer_dddd[345];

    auto g_xy_yy_yy_yz = buffer_dddd[346];

    auto g_xy_yy_yy_zz = buffer_dddd[347];

    auto g_xy_yy_yz_xx = buffer_dddd[348];

    auto g_xy_yy_yz_xy = buffer_dddd[349];

    auto g_xy_yy_yz_xz = buffer_dddd[350];

    auto g_xy_yy_yz_yy = buffer_dddd[351];

    auto g_xy_yy_yz_yz = buffer_dddd[352];

    auto g_xy_yy_yz_zz = buffer_dddd[353];

    auto g_xy_yy_zz_xx = buffer_dddd[354];

    auto g_xy_yy_zz_xy = buffer_dddd[355];

    auto g_xy_yy_zz_xz = buffer_dddd[356];

    auto g_xy_yy_zz_yy = buffer_dddd[357];

    auto g_xy_yy_zz_yz = buffer_dddd[358];

    auto g_xy_yy_zz_zz = buffer_dddd[359];

    auto g_xy_yz_xx_xx = buffer_dddd[360];

    auto g_xy_yz_xx_xy = buffer_dddd[361];

    auto g_xy_yz_xx_xz = buffer_dddd[362];

    auto g_xy_yz_xx_yy = buffer_dddd[363];

    auto g_xy_yz_xx_yz = buffer_dddd[364];

    auto g_xy_yz_xx_zz = buffer_dddd[365];

    auto g_xy_yz_xy_xx = buffer_dddd[366];

    auto g_xy_yz_xy_xy = buffer_dddd[367];

    auto g_xy_yz_xy_xz = buffer_dddd[368];

    auto g_xy_yz_xy_yy = buffer_dddd[369];

    auto g_xy_yz_xy_yz = buffer_dddd[370];

    auto g_xy_yz_xy_zz = buffer_dddd[371];

    auto g_xy_yz_xz_xx = buffer_dddd[372];

    auto g_xy_yz_xz_xy = buffer_dddd[373];

    auto g_xy_yz_xz_xz = buffer_dddd[374];

    auto g_xy_yz_xz_yy = buffer_dddd[375];

    auto g_xy_yz_xz_yz = buffer_dddd[376];

    auto g_xy_yz_xz_zz = buffer_dddd[377];

    auto g_xy_yz_yy_xx = buffer_dddd[378];

    auto g_xy_yz_yy_xy = buffer_dddd[379];

    auto g_xy_yz_yy_xz = buffer_dddd[380];

    auto g_xy_yz_yy_yy = buffer_dddd[381];

    auto g_xy_yz_yy_yz = buffer_dddd[382];

    auto g_xy_yz_yy_zz = buffer_dddd[383];

    auto g_xy_yz_yz_xx = buffer_dddd[384];

    auto g_xy_yz_yz_xy = buffer_dddd[385];

    auto g_xy_yz_yz_xz = buffer_dddd[386];

    auto g_xy_yz_yz_yy = buffer_dddd[387];

    auto g_xy_yz_yz_yz = buffer_dddd[388];

    auto g_xy_yz_yz_zz = buffer_dddd[389];

    auto g_xy_yz_zz_xx = buffer_dddd[390];

    auto g_xy_yz_zz_xy = buffer_dddd[391];

    auto g_xy_yz_zz_xz = buffer_dddd[392];

    auto g_xy_yz_zz_yy = buffer_dddd[393];

    auto g_xy_yz_zz_yz = buffer_dddd[394];

    auto g_xy_yz_zz_zz = buffer_dddd[395];

    auto g_xy_zz_xx_xx = buffer_dddd[396];

    auto g_xy_zz_xx_xy = buffer_dddd[397];

    auto g_xy_zz_xx_xz = buffer_dddd[398];

    auto g_xy_zz_xx_yy = buffer_dddd[399];

    auto g_xy_zz_xx_yz = buffer_dddd[400];

    auto g_xy_zz_xx_zz = buffer_dddd[401];

    auto g_xy_zz_xy_xx = buffer_dddd[402];

    auto g_xy_zz_xy_xy = buffer_dddd[403];

    auto g_xy_zz_xy_xz = buffer_dddd[404];

    auto g_xy_zz_xy_yy = buffer_dddd[405];

    auto g_xy_zz_xy_yz = buffer_dddd[406];

    auto g_xy_zz_xy_zz = buffer_dddd[407];

    auto g_xy_zz_xz_xx = buffer_dddd[408];

    auto g_xy_zz_xz_xy = buffer_dddd[409];

    auto g_xy_zz_xz_xz = buffer_dddd[410];

    auto g_xy_zz_xz_yy = buffer_dddd[411];

    auto g_xy_zz_xz_yz = buffer_dddd[412];

    auto g_xy_zz_xz_zz = buffer_dddd[413];

    auto g_xy_zz_yy_xx = buffer_dddd[414];

    auto g_xy_zz_yy_xy = buffer_dddd[415];

    auto g_xy_zz_yy_xz = buffer_dddd[416];

    auto g_xy_zz_yy_yy = buffer_dddd[417];

    auto g_xy_zz_yy_yz = buffer_dddd[418];

    auto g_xy_zz_yy_zz = buffer_dddd[419];

    auto g_xy_zz_yz_xx = buffer_dddd[420];

    auto g_xy_zz_yz_xy = buffer_dddd[421];

    auto g_xy_zz_yz_xz = buffer_dddd[422];

    auto g_xy_zz_yz_yy = buffer_dddd[423];

    auto g_xy_zz_yz_yz = buffer_dddd[424];

    auto g_xy_zz_yz_zz = buffer_dddd[425];

    auto g_xy_zz_zz_xx = buffer_dddd[426];

    auto g_xy_zz_zz_xy = buffer_dddd[427];

    auto g_xy_zz_zz_xz = buffer_dddd[428];

    auto g_xy_zz_zz_yy = buffer_dddd[429];

    auto g_xy_zz_zz_yz = buffer_dddd[430];

    auto g_xy_zz_zz_zz = buffer_dddd[431];

    auto g_xz_xx_xx_xx = buffer_dddd[432];

    auto g_xz_xx_xx_xy = buffer_dddd[433];

    auto g_xz_xx_xx_xz = buffer_dddd[434];

    auto g_xz_xx_xx_yy = buffer_dddd[435];

    auto g_xz_xx_xx_yz = buffer_dddd[436];

    auto g_xz_xx_xx_zz = buffer_dddd[437];

    auto g_xz_xx_xy_xx = buffer_dddd[438];

    auto g_xz_xx_xy_xy = buffer_dddd[439];

    auto g_xz_xx_xy_xz = buffer_dddd[440];

    auto g_xz_xx_xy_yy = buffer_dddd[441];

    auto g_xz_xx_xy_yz = buffer_dddd[442];

    auto g_xz_xx_xy_zz = buffer_dddd[443];

    auto g_xz_xx_xz_xx = buffer_dddd[444];

    auto g_xz_xx_xz_xy = buffer_dddd[445];

    auto g_xz_xx_xz_xz = buffer_dddd[446];

    auto g_xz_xx_xz_yy = buffer_dddd[447];

    auto g_xz_xx_xz_yz = buffer_dddd[448];

    auto g_xz_xx_xz_zz = buffer_dddd[449];

    auto g_xz_xx_yy_xx = buffer_dddd[450];

    auto g_xz_xx_yy_xy = buffer_dddd[451];

    auto g_xz_xx_yy_xz = buffer_dddd[452];

    auto g_xz_xx_yy_yy = buffer_dddd[453];

    auto g_xz_xx_yy_yz = buffer_dddd[454];

    auto g_xz_xx_yy_zz = buffer_dddd[455];

    auto g_xz_xx_yz_xx = buffer_dddd[456];

    auto g_xz_xx_yz_xy = buffer_dddd[457];

    auto g_xz_xx_yz_xz = buffer_dddd[458];

    auto g_xz_xx_yz_yy = buffer_dddd[459];

    auto g_xz_xx_yz_yz = buffer_dddd[460];

    auto g_xz_xx_yz_zz = buffer_dddd[461];

    auto g_xz_xx_zz_xx = buffer_dddd[462];

    auto g_xz_xx_zz_xy = buffer_dddd[463];

    auto g_xz_xx_zz_xz = buffer_dddd[464];

    auto g_xz_xx_zz_yy = buffer_dddd[465];

    auto g_xz_xx_zz_yz = buffer_dddd[466];

    auto g_xz_xx_zz_zz = buffer_dddd[467];

    auto g_xz_xy_xx_xx = buffer_dddd[468];

    auto g_xz_xy_xx_xy = buffer_dddd[469];

    auto g_xz_xy_xx_xz = buffer_dddd[470];

    auto g_xz_xy_xx_yy = buffer_dddd[471];

    auto g_xz_xy_xx_yz = buffer_dddd[472];

    auto g_xz_xy_xx_zz = buffer_dddd[473];

    auto g_xz_xy_xy_xx = buffer_dddd[474];

    auto g_xz_xy_xy_xy = buffer_dddd[475];

    auto g_xz_xy_xy_xz = buffer_dddd[476];

    auto g_xz_xy_xy_yy = buffer_dddd[477];

    auto g_xz_xy_xy_yz = buffer_dddd[478];

    auto g_xz_xy_xy_zz = buffer_dddd[479];

    auto g_xz_xy_xz_xx = buffer_dddd[480];

    auto g_xz_xy_xz_xy = buffer_dddd[481];

    auto g_xz_xy_xz_xz = buffer_dddd[482];

    auto g_xz_xy_xz_yy = buffer_dddd[483];

    auto g_xz_xy_xz_yz = buffer_dddd[484];

    auto g_xz_xy_xz_zz = buffer_dddd[485];

    auto g_xz_xy_yy_xx = buffer_dddd[486];

    auto g_xz_xy_yy_xy = buffer_dddd[487];

    auto g_xz_xy_yy_xz = buffer_dddd[488];

    auto g_xz_xy_yy_yy = buffer_dddd[489];

    auto g_xz_xy_yy_yz = buffer_dddd[490];

    auto g_xz_xy_yy_zz = buffer_dddd[491];

    auto g_xz_xy_yz_xx = buffer_dddd[492];

    auto g_xz_xy_yz_xy = buffer_dddd[493];

    auto g_xz_xy_yz_xz = buffer_dddd[494];

    auto g_xz_xy_yz_yy = buffer_dddd[495];

    auto g_xz_xy_yz_yz = buffer_dddd[496];

    auto g_xz_xy_yz_zz = buffer_dddd[497];

    auto g_xz_xy_zz_xx = buffer_dddd[498];

    auto g_xz_xy_zz_xy = buffer_dddd[499];

    auto g_xz_xy_zz_xz = buffer_dddd[500];

    auto g_xz_xy_zz_yy = buffer_dddd[501];

    auto g_xz_xy_zz_yz = buffer_dddd[502];

    auto g_xz_xy_zz_zz = buffer_dddd[503];

    auto g_xz_xz_xx_xx = buffer_dddd[504];

    auto g_xz_xz_xx_xy = buffer_dddd[505];

    auto g_xz_xz_xx_xz = buffer_dddd[506];

    auto g_xz_xz_xx_yy = buffer_dddd[507];

    auto g_xz_xz_xx_yz = buffer_dddd[508];

    auto g_xz_xz_xx_zz = buffer_dddd[509];

    auto g_xz_xz_xy_xx = buffer_dddd[510];

    auto g_xz_xz_xy_xy = buffer_dddd[511];

    auto g_xz_xz_xy_xz = buffer_dddd[512];

    auto g_xz_xz_xy_yy = buffer_dddd[513];

    auto g_xz_xz_xy_yz = buffer_dddd[514];

    auto g_xz_xz_xy_zz = buffer_dddd[515];

    auto g_xz_xz_xz_xx = buffer_dddd[516];

    auto g_xz_xz_xz_xy = buffer_dddd[517];

    auto g_xz_xz_xz_xz = buffer_dddd[518];

    auto g_xz_xz_xz_yy = buffer_dddd[519];

    auto g_xz_xz_xz_yz = buffer_dddd[520];

    auto g_xz_xz_xz_zz = buffer_dddd[521];

    auto g_xz_xz_yy_xx = buffer_dddd[522];

    auto g_xz_xz_yy_xy = buffer_dddd[523];

    auto g_xz_xz_yy_xz = buffer_dddd[524];

    auto g_xz_xz_yy_yy = buffer_dddd[525];

    auto g_xz_xz_yy_yz = buffer_dddd[526];

    auto g_xz_xz_yy_zz = buffer_dddd[527];

    auto g_xz_xz_yz_xx = buffer_dddd[528];

    auto g_xz_xz_yz_xy = buffer_dddd[529];

    auto g_xz_xz_yz_xz = buffer_dddd[530];

    auto g_xz_xz_yz_yy = buffer_dddd[531];

    auto g_xz_xz_yz_yz = buffer_dddd[532];

    auto g_xz_xz_yz_zz = buffer_dddd[533];

    auto g_xz_xz_zz_xx = buffer_dddd[534];

    auto g_xz_xz_zz_xy = buffer_dddd[535];

    auto g_xz_xz_zz_xz = buffer_dddd[536];

    auto g_xz_xz_zz_yy = buffer_dddd[537];

    auto g_xz_xz_zz_yz = buffer_dddd[538];

    auto g_xz_xz_zz_zz = buffer_dddd[539];

    auto g_xz_yy_xx_xx = buffer_dddd[540];

    auto g_xz_yy_xx_xy = buffer_dddd[541];

    auto g_xz_yy_xx_xz = buffer_dddd[542];

    auto g_xz_yy_xx_yy = buffer_dddd[543];

    auto g_xz_yy_xx_yz = buffer_dddd[544];

    auto g_xz_yy_xx_zz = buffer_dddd[545];

    auto g_xz_yy_xy_xx = buffer_dddd[546];

    auto g_xz_yy_xy_xy = buffer_dddd[547];

    auto g_xz_yy_xy_xz = buffer_dddd[548];

    auto g_xz_yy_xy_yy = buffer_dddd[549];

    auto g_xz_yy_xy_yz = buffer_dddd[550];

    auto g_xz_yy_xy_zz = buffer_dddd[551];

    auto g_xz_yy_xz_xx = buffer_dddd[552];

    auto g_xz_yy_xz_xy = buffer_dddd[553];

    auto g_xz_yy_xz_xz = buffer_dddd[554];

    auto g_xz_yy_xz_yy = buffer_dddd[555];

    auto g_xz_yy_xz_yz = buffer_dddd[556];

    auto g_xz_yy_xz_zz = buffer_dddd[557];

    auto g_xz_yy_yy_xx = buffer_dddd[558];

    auto g_xz_yy_yy_xy = buffer_dddd[559];

    auto g_xz_yy_yy_xz = buffer_dddd[560];

    auto g_xz_yy_yy_yy = buffer_dddd[561];

    auto g_xz_yy_yy_yz = buffer_dddd[562];

    auto g_xz_yy_yy_zz = buffer_dddd[563];

    auto g_xz_yy_yz_xx = buffer_dddd[564];

    auto g_xz_yy_yz_xy = buffer_dddd[565];

    auto g_xz_yy_yz_xz = buffer_dddd[566];

    auto g_xz_yy_yz_yy = buffer_dddd[567];

    auto g_xz_yy_yz_yz = buffer_dddd[568];

    auto g_xz_yy_yz_zz = buffer_dddd[569];

    auto g_xz_yy_zz_xx = buffer_dddd[570];

    auto g_xz_yy_zz_xy = buffer_dddd[571];

    auto g_xz_yy_zz_xz = buffer_dddd[572];

    auto g_xz_yy_zz_yy = buffer_dddd[573];

    auto g_xz_yy_zz_yz = buffer_dddd[574];

    auto g_xz_yy_zz_zz = buffer_dddd[575];

    auto g_xz_yz_xx_xx = buffer_dddd[576];

    auto g_xz_yz_xx_xy = buffer_dddd[577];

    auto g_xz_yz_xx_xz = buffer_dddd[578];

    auto g_xz_yz_xx_yy = buffer_dddd[579];

    auto g_xz_yz_xx_yz = buffer_dddd[580];

    auto g_xz_yz_xx_zz = buffer_dddd[581];

    auto g_xz_yz_xy_xx = buffer_dddd[582];

    auto g_xz_yz_xy_xy = buffer_dddd[583];

    auto g_xz_yz_xy_xz = buffer_dddd[584];

    auto g_xz_yz_xy_yy = buffer_dddd[585];

    auto g_xz_yz_xy_yz = buffer_dddd[586];

    auto g_xz_yz_xy_zz = buffer_dddd[587];

    auto g_xz_yz_xz_xx = buffer_dddd[588];

    auto g_xz_yz_xz_xy = buffer_dddd[589];

    auto g_xz_yz_xz_xz = buffer_dddd[590];

    auto g_xz_yz_xz_yy = buffer_dddd[591];

    auto g_xz_yz_xz_yz = buffer_dddd[592];

    auto g_xz_yz_xz_zz = buffer_dddd[593];

    auto g_xz_yz_yy_xx = buffer_dddd[594];

    auto g_xz_yz_yy_xy = buffer_dddd[595];

    auto g_xz_yz_yy_xz = buffer_dddd[596];

    auto g_xz_yz_yy_yy = buffer_dddd[597];

    auto g_xz_yz_yy_yz = buffer_dddd[598];

    auto g_xz_yz_yy_zz = buffer_dddd[599];

    auto g_xz_yz_yz_xx = buffer_dddd[600];

    auto g_xz_yz_yz_xy = buffer_dddd[601];

    auto g_xz_yz_yz_xz = buffer_dddd[602];

    auto g_xz_yz_yz_yy = buffer_dddd[603];

    auto g_xz_yz_yz_yz = buffer_dddd[604];

    auto g_xz_yz_yz_zz = buffer_dddd[605];

    auto g_xz_yz_zz_xx = buffer_dddd[606];

    auto g_xz_yz_zz_xy = buffer_dddd[607];

    auto g_xz_yz_zz_xz = buffer_dddd[608];

    auto g_xz_yz_zz_yy = buffer_dddd[609];

    auto g_xz_yz_zz_yz = buffer_dddd[610];

    auto g_xz_yz_zz_zz = buffer_dddd[611];

    auto g_xz_zz_xx_xx = buffer_dddd[612];

    auto g_xz_zz_xx_xy = buffer_dddd[613];

    auto g_xz_zz_xx_xz = buffer_dddd[614];

    auto g_xz_zz_xx_yy = buffer_dddd[615];

    auto g_xz_zz_xx_yz = buffer_dddd[616];

    auto g_xz_zz_xx_zz = buffer_dddd[617];

    auto g_xz_zz_xy_xx = buffer_dddd[618];

    auto g_xz_zz_xy_xy = buffer_dddd[619];

    auto g_xz_zz_xy_xz = buffer_dddd[620];

    auto g_xz_zz_xy_yy = buffer_dddd[621];

    auto g_xz_zz_xy_yz = buffer_dddd[622];

    auto g_xz_zz_xy_zz = buffer_dddd[623];

    auto g_xz_zz_xz_xx = buffer_dddd[624];

    auto g_xz_zz_xz_xy = buffer_dddd[625];

    auto g_xz_zz_xz_xz = buffer_dddd[626];

    auto g_xz_zz_xz_yy = buffer_dddd[627];

    auto g_xz_zz_xz_yz = buffer_dddd[628];

    auto g_xz_zz_xz_zz = buffer_dddd[629];

    auto g_xz_zz_yy_xx = buffer_dddd[630];

    auto g_xz_zz_yy_xy = buffer_dddd[631];

    auto g_xz_zz_yy_xz = buffer_dddd[632];

    auto g_xz_zz_yy_yy = buffer_dddd[633];

    auto g_xz_zz_yy_yz = buffer_dddd[634];

    auto g_xz_zz_yy_zz = buffer_dddd[635];

    auto g_xz_zz_yz_xx = buffer_dddd[636];

    auto g_xz_zz_yz_xy = buffer_dddd[637];

    auto g_xz_zz_yz_xz = buffer_dddd[638];

    auto g_xz_zz_yz_yy = buffer_dddd[639];

    auto g_xz_zz_yz_yz = buffer_dddd[640];

    auto g_xz_zz_yz_zz = buffer_dddd[641];

    auto g_xz_zz_zz_xx = buffer_dddd[642];

    auto g_xz_zz_zz_xy = buffer_dddd[643];

    auto g_xz_zz_zz_xz = buffer_dddd[644];

    auto g_xz_zz_zz_yy = buffer_dddd[645];

    auto g_xz_zz_zz_yz = buffer_dddd[646];

    auto g_xz_zz_zz_zz = buffer_dddd[647];

    auto g_yy_xx_xx_xx = buffer_dddd[648];

    auto g_yy_xx_xx_xy = buffer_dddd[649];

    auto g_yy_xx_xx_xz = buffer_dddd[650];

    auto g_yy_xx_xx_yy = buffer_dddd[651];

    auto g_yy_xx_xx_yz = buffer_dddd[652];

    auto g_yy_xx_xx_zz = buffer_dddd[653];

    auto g_yy_xx_xy_xx = buffer_dddd[654];

    auto g_yy_xx_xy_xy = buffer_dddd[655];

    auto g_yy_xx_xy_xz = buffer_dddd[656];

    auto g_yy_xx_xy_yy = buffer_dddd[657];

    auto g_yy_xx_xy_yz = buffer_dddd[658];

    auto g_yy_xx_xy_zz = buffer_dddd[659];

    auto g_yy_xx_xz_xx = buffer_dddd[660];

    auto g_yy_xx_xz_xy = buffer_dddd[661];

    auto g_yy_xx_xz_xz = buffer_dddd[662];

    auto g_yy_xx_xz_yy = buffer_dddd[663];

    auto g_yy_xx_xz_yz = buffer_dddd[664];

    auto g_yy_xx_xz_zz = buffer_dddd[665];

    auto g_yy_xx_yy_xx = buffer_dddd[666];

    auto g_yy_xx_yy_xy = buffer_dddd[667];

    auto g_yy_xx_yy_xz = buffer_dddd[668];

    auto g_yy_xx_yy_yy = buffer_dddd[669];

    auto g_yy_xx_yy_yz = buffer_dddd[670];

    auto g_yy_xx_yy_zz = buffer_dddd[671];

    auto g_yy_xx_yz_xx = buffer_dddd[672];

    auto g_yy_xx_yz_xy = buffer_dddd[673];

    auto g_yy_xx_yz_xz = buffer_dddd[674];

    auto g_yy_xx_yz_yy = buffer_dddd[675];

    auto g_yy_xx_yz_yz = buffer_dddd[676];

    auto g_yy_xx_yz_zz = buffer_dddd[677];

    auto g_yy_xx_zz_xx = buffer_dddd[678];

    auto g_yy_xx_zz_xy = buffer_dddd[679];

    auto g_yy_xx_zz_xz = buffer_dddd[680];

    auto g_yy_xx_zz_yy = buffer_dddd[681];

    auto g_yy_xx_zz_yz = buffer_dddd[682];

    auto g_yy_xx_zz_zz = buffer_dddd[683];

    auto g_yy_xy_xx_xx = buffer_dddd[684];

    auto g_yy_xy_xx_xy = buffer_dddd[685];

    auto g_yy_xy_xx_xz = buffer_dddd[686];

    auto g_yy_xy_xx_yy = buffer_dddd[687];

    auto g_yy_xy_xx_yz = buffer_dddd[688];

    auto g_yy_xy_xx_zz = buffer_dddd[689];

    auto g_yy_xy_xy_xx = buffer_dddd[690];

    auto g_yy_xy_xy_xy = buffer_dddd[691];

    auto g_yy_xy_xy_xz = buffer_dddd[692];

    auto g_yy_xy_xy_yy = buffer_dddd[693];

    auto g_yy_xy_xy_yz = buffer_dddd[694];

    auto g_yy_xy_xy_zz = buffer_dddd[695];

    auto g_yy_xy_xz_xx = buffer_dddd[696];

    auto g_yy_xy_xz_xy = buffer_dddd[697];

    auto g_yy_xy_xz_xz = buffer_dddd[698];

    auto g_yy_xy_xz_yy = buffer_dddd[699];

    auto g_yy_xy_xz_yz = buffer_dddd[700];

    auto g_yy_xy_xz_zz = buffer_dddd[701];

    auto g_yy_xy_yy_xx = buffer_dddd[702];

    auto g_yy_xy_yy_xy = buffer_dddd[703];

    auto g_yy_xy_yy_xz = buffer_dddd[704];

    auto g_yy_xy_yy_yy = buffer_dddd[705];

    auto g_yy_xy_yy_yz = buffer_dddd[706];

    auto g_yy_xy_yy_zz = buffer_dddd[707];

    auto g_yy_xy_yz_xx = buffer_dddd[708];

    auto g_yy_xy_yz_xy = buffer_dddd[709];

    auto g_yy_xy_yz_xz = buffer_dddd[710];

    auto g_yy_xy_yz_yy = buffer_dddd[711];

    auto g_yy_xy_yz_yz = buffer_dddd[712];

    auto g_yy_xy_yz_zz = buffer_dddd[713];

    auto g_yy_xy_zz_xx = buffer_dddd[714];

    auto g_yy_xy_zz_xy = buffer_dddd[715];

    auto g_yy_xy_zz_xz = buffer_dddd[716];

    auto g_yy_xy_zz_yy = buffer_dddd[717];

    auto g_yy_xy_zz_yz = buffer_dddd[718];

    auto g_yy_xy_zz_zz = buffer_dddd[719];

    auto g_yy_xz_xx_xx = buffer_dddd[720];

    auto g_yy_xz_xx_xy = buffer_dddd[721];

    auto g_yy_xz_xx_xz = buffer_dddd[722];

    auto g_yy_xz_xx_yy = buffer_dddd[723];

    auto g_yy_xz_xx_yz = buffer_dddd[724];

    auto g_yy_xz_xx_zz = buffer_dddd[725];

    auto g_yy_xz_xy_xx = buffer_dddd[726];

    auto g_yy_xz_xy_xy = buffer_dddd[727];

    auto g_yy_xz_xy_xz = buffer_dddd[728];

    auto g_yy_xz_xy_yy = buffer_dddd[729];

    auto g_yy_xz_xy_yz = buffer_dddd[730];

    auto g_yy_xz_xy_zz = buffer_dddd[731];

    auto g_yy_xz_xz_xx = buffer_dddd[732];

    auto g_yy_xz_xz_xy = buffer_dddd[733];

    auto g_yy_xz_xz_xz = buffer_dddd[734];

    auto g_yy_xz_xz_yy = buffer_dddd[735];

    auto g_yy_xz_xz_yz = buffer_dddd[736];

    auto g_yy_xz_xz_zz = buffer_dddd[737];

    auto g_yy_xz_yy_xx = buffer_dddd[738];

    auto g_yy_xz_yy_xy = buffer_dddd[739];

    auto g_yy_xz_yy_xz = buffer_dddd[740];

    auto g_yy_xz_yy_yy = buffer_dddd[741];

    auto g_yy_xz_yy_yz = buffer_dddd[742];

    auto g_yy_xz_yy_zz = buffer_dddd[743];

    auto g_yy_xz_yz_xx = buffer_dddd[744];

    auto g_yy_xz_yz_xy = buffer_dddd[745];

    auto g_yy_xz_yz_xz = buffer_dddd[746];

    auto g_yy_xz_yz_yy = buffer_dddd[747];

    auto g_yy_xz_yz_yz = buffer_dddd[748];

    auto g_yy_xz_yz_zz = buffer_dddd[749];

    auto g_yy_xz_zz_xx = buffer_dddd[750];

    auto g_yy_xz_zz_xy = buffer_dddd[751];

    auto g_yy_xz_zz_xz = buffer_dddd[752];

    auto g_yy_xz_zz_yy = buffer_dddd[753];

    auto g_yy_xz_zz_yz = buffer_dddd[754];

    auto g_yy_xz_zz_zz = buffer_dddd[755];

    auto g_yy_yy_xx_xx = buffer_dddd[756];

    auto g_yy_yy_xx_xy = buffer_dddd[757];

    auto g_yy_yy_xx_xz = buffer_dddd[758];

    auto g_yy_yy_xx_yy = buffer_dddd[759];

    auto g_yy_yy_xx_yz = buffer_dddd[760];

    auto g_yy_yy_xx_zz = buffer_dddd[761];

    auto g_yy_yy_xy_xx = buffer_dddd[762];

    auto g_yy_yy_xy_xy = buffer_dddd[763];

    auto g_yy_yy_xy_xz = buffer_dddd[764];

    auto g_yy_yy_xy_yy = buffer_dddd[765];

    auto g_yy_yy_xy_yz = buffer_dddd[766];

    auto g_yy_yy_xy_zz = buffer_dddd[767];

    auto g_yy_yy_xz_xx = buffer_dddd[768];

    auto g_yy_yy_xz_xy = buffer_dddd[769];

    auto g_yy_yy_xz_xz = buffer_dddd[770];

    auto g_yy_yy_xz_yy = buffer_dddd[771];

    auto g_yy_yy_xz_yz = buffer_dddd[772];

    auto g_yy_yy_xz_zz = buffer_dddd[773];

    auto g_yy_yy_yy_xx = buffer_dddd[774];

    auto g_yy_yy_yy_xy = buffer_dddd[775];

    auto g_yy_yy_yy_xz = buffer_dddd[776];

    auto g_yy_yy_yy_yy = buffer_dddd[777];

    auto g_yy_yy_yy_yz = buffer_dddd[778];

    auto g_yy_yy_yy_zz = buffer_dddd[779];

    auto g_yy_yy_yz_xx = buffer_dddd[780];

    auto g_yy_yy_yz_xy = buffer_dddd[781];

    auto g_yy_yy_yz_xz = buffer_dddd[782];

    auto g_yy_yy_yz_yy = buffer_dddd[783];

    auto g_yy_yy_yz_yz = buffer_dddd[784];

    auto g_yy_yy_yz_zz = buffer_dddd[785];

    auto g_yy_yy_zz_xx = buffer_dddd[786];

    auto g_yy_yy_zz_xy = buffer_dddd[787];

    auto g_yy_yy_zz_xz = buffer_dddd[788];

    auto g_yy_yy_zz_yy = buffer_dddd[789];

    auto g_yy_yy_zz_yz = buffer_dddd[790];

    auto g_yy_yy_zz_zz = buffer_dddd[791];

    auto g_yy_yz_xx_xx = buffer_dddd[792];

    auto g_yy_yz_xx_xy = buffer_dddd[793];

    auto g_yy_yz_xx_xz = buffer_dddd[794];

    auto g_yy_yz_xx_yy = buffer_dddd[795];

    auto g_yy_yz_xx_yz = buffer_dddd[796];

    auto g_yy_yz_xx_zz = buffer_dddd[797];

    auto g_yy_yz_xy_xx = buffer_dddd[798];

    auto g_yy_yz_xy_xy = buffer_dddd[799];

    auto g_yy_yz_xy_xz = buffer_dddd[800];

    auto g_yy_yz_xy_yy = buffer_dddd[801];

    auto g_yy_yz_xy_yz = buffer_dddd[802];

    auto g_yy_yz_xy_zz = buffer_dddd[803];

    auto g_yy_yz_xz_xx = buffer_dddd[804];

    auto g_yy_yz_xz_xy = buffer_dddd[805];

    auto g_yy_yz_xz_xz = buffer_dddd[806];

    auto g_yy_yz_xz_yy = buffer_dddd[807];

    auto g_yy_yz_xz_yz = buffer_dddd[808];

    auto g_yy_yz_xz_zz = buffer_dddd[809];

    auto g_yy_yz_yy_xx = buffer_dddd[810];

    auto g_yy_yz_yy_xy = buffer_dddd[811];

    auto g_yy_yz_yy_xz = buffer_dddd[812];

    auto g_yy_yz_yy_yy = buffer_dddd[813];

    auto g_yy_yz_yy_yz = buffer_dddd[814];

    auto g_yy_yz_yy_zz = buffer_dddd[815];

    auto g_yy_yz_yz_xx = buffer_dddd[816];

    auto g_yy_yz_yz_xy = buffer_dddd[817];

    auto g_yy_yz_yz_xz = buffer_dddd[818];

    auto g_yy_yz_yz_yy = buffer_dddd[819];

    auto g_yy_yz_yz_yz = buffer_dddd[820];

    auto g_yy_yz_yz_zz = buffer_dddd[821];

    auto g_yy_yz_zz_xx = buffer_dddd[822];

    auto g_yy_yz_zz_xy = buffer_dddd[823];

    auto g_yy_yz_zz_xz = buffer_dddd[824];

    auto g_yy_yz_zz_yy = buffer_dddd[825];

    auto g_yy_yz_zz_yz = buffer_dddd[826];

    auto g_yy_yz_zz_zz = buffer_dddd[827];

    auto g_yy_zz_xx_xx = buffer_dddd[828];

    auto g_yy_zz_xx_xy = buffer_dddd[829];

    auto g_yy_zz_xx_xz = buffer_dddd[830];

    auto g_yy_zz_xx_yy = buffer_dddd[831];

    auto g_yy_zz_xx_yz = buffer_dddd[832];

    auto g_yy_zz_xx_zz = buffer_dddd[833];

    auto g_yy_zz_xy_xx = buffer_dddd[834];

    auto g_yy_zz_xy_xy = buffer_dddd[835];

    auto g_yy_zz_xy_xz = buffer_dddd[836];

    auto g_yy_zz_xy_yy = buffer_dddd[837];

    auto g_yy_zz_xy_yz = buffer_dddd[838];

    auto g_yy_zz_xy_zz = buffer_dddd[839];

    auto g_yy_zz_xz_xx = buffer_dddd[840];

    auto g_yy_zz_xz_xy = buffer_dddd[841];

    auto g_yy_zz_xz_xz = buffer_dddd[842];

    auto g_yy_zz_xz_yy = buffer_dddd[843];

    auto g_yy_zz_xz_yz = buffer_dddd[844];

    auto g_yy_zz_xz_zz = buffer_dddd[845];

    auto g_yy_zz_yy_xx = buffer_dddd[846];

    auto g_yy_zz_yy_xy = buffer_dddd[847];

    auto g_yy_zz_yy_xz = buffer_dddd[848];

    auto g_yy_zz_yy_yy = buffer_dddd[849];

    auto g_yy_zz_yy_yz = buffer_dddd[850];

    auto g_yy_zz_yy_zz = buffer_dddd[851];

    auto g_yy_zz_yz_xx = buffer_dddd[852];

    auto g_yy_zz_yz_xy = buffer_dddd[853];

    auto g_yy_zz_yz_xz = buffer_dddd[854];

    auto g_yy_zz_yz_yy = buffer_dddd[855];

    auto g_yy_zz_yz_yz = buffer_dddd[856];

    auto g_yy_zz_yz_zz = buffer_dddd[857];

    auto g_yy_zz_zz_xx = buffer_dddd[858];

    auto g_yy_zz_zz_xy = buffer_dddd[859];

    auto g_yy_zz_zz_xz = buffer_dddd[860];

    auto g_yy_zz_zz_yy = buffer_dddd[861];

    auto g_yy_zz_zz_yz = buffer_dddd[862];

    auto g_yy_zz_zz_zz = buffer_dddd[863];

    auto g_yz_xx_xx_xx = buffer_dddd[864];

    auto g_yz_xx_xx_xy = buffer_dddd[865];

    auto g_yz_xx_xx_xz = buffer_dddd[866];

    auto g_yz_xx_xx_yy = buffer_dddd[867];

    auto g_yz_xx_xx_yz = buffer_dddd[868];

    auto g_yz_xx_xx_zz = buffer_dddd[869];

    auto g_yz_xx_xy_xx = buffer_dddd[870];

    auto g_yz_xx_xy_xy = buffer_dddd[871];

    auto g_yz_xx_xy_xz = buffer_dddd[872];

    auto g_yz_xx_xy_yy = buffer_dddd[873];

    auto g_yz_xx_xy_yz = buffer_dddd[874];

    auto g_yz_xx_xy_zz = buffer_dddd[875];

    auto g_yz_xx_xz_xx = buffer_dddd[876];

    auto g_yz_xx_xz_xy = buffer_dddd[877];

    auto g_yz_xx_xz_xz = buffer_dddd[878];

    auto g_yz_xx_xz_yy = buffer_dddd[879];

    auto g_yz_xx_xz_yz = buffer_dddd[880];

    auto g_yz_xx_xz_zz = buffer_dddd[881];

    auto g_yz_xx_yy_xx = buffer_dddd[882];

    auto g_yz_xx_yy_xy = buffer_dddd[883];

    auto g_yz_xx_yy_xz = buffer_dddd[884];

    auto g_yz_xx_yy_yy = buffer_dddd[885];

    auto g_yz_xx_yy_yz = buffer_dddd[886];

    auto g_yz_xx_yy_zz = buffer_dddd[887];

    auto g_yz_xx_yz_xx = buffer_dddd[888];

    auto g_yz_xx_yz_xy = buffer_dddd[889];

    auto g_yz_xx_yz_xz = buffer_dddd[890];

    auto g_yz_xx_yz_yy = buffer_dddd[891];

    auto g_yz_xx_yz_yz = buffer_dddd[892];

    auto g_yz_xx_yz_zz = buffer_dddd[893];

    auto g_yz_xx_zz_xx = buffer_dddd[894];

    auto g_yz_xx_zz_xy = buffer_dddd[895];

    auto g_yz_xx_zz_xz = buffer_dddd[896];

    auto g_yz_xx_zz_yy = buffer_dddd[897];

    auto g_yz_xx_zz_yz = buffer_dddd[898];

    auto g_yz_xx_zz_zz = buffer_dddd[899];

    auto g_yz_xy_xx_xx = buffer_dddd[900];

    auto g_yz_xy_xx_xy = buffer_dddd[901];

    auto g_yz_xy_xx_xz = buffer_dddd[902];

    auto g_yz_xy_xx_yy = buffer_dddd[903];

    auto g_yz_xy_xx_yz = buffer_dddd[904];

    auto g_yz_xy_xx_zz = buffer_dddd[905];

    auto g_yz_xy_xy_xx = buffer_dddd[906];

    auto g_yz_xy_xy_xy = buffer_dddd[907];

    auto g_yz_xy_xy_xz = buffer_dddd[908];

    auto g_yz_xy_xy_yy = buffer_dddd[909];

    auto g_yz_xy_xy_yz = buffer_dddd[910];

    auto g_yz_xy_xy_zz = buffer_dddd[911];

    auto g_yz_xy_xz_xx = buffer_dddd[912];

    auto g_yz_xy_xz_xy = buffer_dddd[913];

    auto g_yz_xy_xz_xz = buffer_dddd[914];

    auto g_yz_xy_xz_yy = buffer_dddd[915];

    auto g_yz_xy_xz_yz = buffer_dddd[916];

    auto g_yz_xy_xz_zz = buffer_dddd[917];

    auto g_yz_xy_yy_xx = buffer_dddd[918];

    auto g_yz_xy_yy_xy = buffer_dddd[919];

    auto g_yz_xy_yy_xz = buffer_dddd[920];

    auto g_yz_xy_yy_yy = buffer_dddd[921];

    auto g_yz_xy_yy_yz = buffer_dddd[922];

    auto g_yz_xy_yy_zz = buffer_dddd[923];

    auto g_yz_xy_yz_xx = buffer_dddd[924];

    auto g_yz_xy_yz_xy = buffer_dddd[925];

    auto g_yz_xy_yz_xz = buffer_dddd[926];

    auto g_yz_xy_yz_yy = buffer_dddd[927];

    auto g_yz_xy_yz_yz = buffer_dddd[928];

    auto g_yz_xy_yz_zz = buffer_dddd[929];

    auto g_yz_xy_zz_xx = buffer_dddd[930];

    auto g_yz_xy_zz_xy = buffer_dddd[931];

    auto g_yz_xy_zz_xz = buffer_dddd[932];

    auto g_yz_xy_zz_yy = buffer_dddd[933];

    auto g_yz_xy_zz_yz = buffer_dddd[934];

    auto g_yz_xy_zz_zz = buffer_dddd[935];

    auto g_yz_xz_xx_xx = buffer_dddd[936];

    auto g_yz_xz_xx_xy = buffer_dddd[937];

    auto g_yz_xz_xx_xz = buffer_dddd[938];

    auto g_yz_xz_xx_yy = buffer_dddd[939];

    auto g_yz_xz_xx_yz = buffer_dddd[940];

    auto g_yz_xz_xx_zz = buffer_dddd[941];

    auto g_yz_xz_xy_xx = buffer_dddd[942];

    auto g_yz_xz_xy_xy = buffer_dddd[943];

    auto g_yz_xz_xy_xz = buffer_dddd[944];

    auto g_yz_xz_xy_yy = buffer_dddd[945];

    auto g_yz_xz_xy_yz = buffer_dddd[946];

    auto g_yz_xz_xy_zz = buffer_dddd[947];

    auto g_yz_xz_xz_xx = buffer_dddd[948];

    auto g_yz_xz_xz_xy = buffer_dddd[949];

    auto g_yz_xz_xz_xz = buffer_dddd[950];

    auto g_yz_xz_xz_yy = buffer_dddd[951];

    auto g_yz_xz_xz_yz = buffer_dddd[952];

    auto g_yz_xz_xz_zz = buffer_dddd[953];

    auto g_yz_xz_yy_xx = buffer_dddd[954];

    auto g_yz_xz_yy_xy = buffer_dddd[955];

    auto g_yz_xz_yy_xz = buffer_dddd[956];

    auto g_yz_xz_yy_yy = buffer_dddd[957];

    auto g_yz_xz_yy_yz = buffer_dddd[958];

    auto g_yz_xz_yy_zz = buffer_dddd[959];

    auto g_yz_xz_yz_xx = buffer_dddd[960];

    auto g_yz_xz_yz_xy = buffer_dddd[961];

    auto g_yz_xz_yz_xz = buffer_dddd[962];

    auto g_yz_xz_yz_yy = buffer_dddd[963];

    auto g_yz_xz_yz_yz = buffer_dddd[964];

    auto g_yz_xz_yz_zz = buffer_dddd[965];

    auto g_yz_xz_zz_xx = buffer_dddd[966];

    auto g_yz_xz_zz_xy = buffer_dddd[967];

    auto g_yz_xz_zz_xz = buffer_dddd[968];

    auto g_yz_xz_zz_yy = buffer_dddd[969];

    auto g_yz_xz_zz_yz = buffer_dddd[970];

    auto g_yz_xz_zz_zz = buffer_dddd[971];

    auto g_yz_yy_xx_xx = buffer_dddd[972];

    auto g_yz_yy_xx_xy = buffer_dddd[973];

    auto g_yz_yy_xx_xz = buffer_dddd[974];

    auto g_yz_yy_xx_yy = buffer_dddd[975];

    auto g_yz_yy_xx_yz = buffer_dddd[976];

    auto g_yz_yy_xx_zz = buffer_dddd[977];

    auto g_yz_yy_xy_xx = buffer_dddd[978];

    auto g_yz_yy_xy_xy = buffer_dddd[979];

    auto g_yz_yy_xy_xz = buffer_dddd[980];

    auto g_yz_yy_xy_yy = buffer_dddd[981];

    auto g_yz_yy_xy_yz = buffer_dddd[982];

    auto g_yz_yy_xy_zz = buffer_dddd[983];

    auto g_yz_yy_xz_xx = buffer_dddd[984];

    auto g_yz_yy_xz_xy = buffer_dddd[985];

    auto g_yz_yy_xz_xz = buffer_dddd[986];

    auto g_yz_yy_xz_yy = buffer_dddd[987];

    auto g_yz_yy_xz_yz = buffer_dddd[988];

    auto g_yz_yy_xz_zz = buffer_dddd[989];

    auto g_yz_yy_yy_xx = buffer_dddd[990];

    auto g_yz_yy_yy_xy = buffer_dddd[991];

    auto g_yz_yy_yy_xz = buffer_dddd[992];

    auto g_yz_yy_yy_yy = buffer_dddd[993];

    auto g_yz_yy_yy_yz = buffer_dddd[994];

    auto g_yz_yy_yy_zz = buffer_dddd[995];

    auto g_yz_yy_yz_xx = buffer_dddd[996];

    auto g_yz_yy_yz_xy = buffer_dddd[997];

    auto g_yz_yy_yz_xz = buffer_dddd[998];

    auto g_yz_yy_yz_yy = buffer_dddd[999];

    auto g_yz_yy_yz_yz = buffer_dddd[1000];

    auto g_yz_yy_yz_zz = buffer_dddd[1001];

    auto g_yz_yy_zz_xx = buffer_dddd[1002];

    auto g_yz_yy_zz_xy = buffer_dddd[1003];

    auto g_yz_yy_zz_xz = buffer_dddd[1004];

    auto g_yz_yy_zz_yy = buffer_dddd[1005];

    auto g_yz_yy_zz_yz = buffer_dddd[1006];

    auto g_yz_yy_zz_zz = buffer_dddd[1007];

    auto g_yz_yz_xx_xx = buffer_dddd[1008];

    auto g_yz_yz_xx_xy = buffer_dddd[1009];

    auto g_yz_yz_xx_xz = buffer_dddd[1010];

    auto g_yz_yz_xx_yy = buffer_dddd[1011];

    auto g_yz_yz_xx_yz = buffer_dddd[1012];

    auto g_yz_yz_xx_zz = buffer_dddd[1013];

    auto g_yz_yz_xy_xx = buffer_dddd[1014];

    auto g_yz_yz_xy_xy = buffer_dddd[1015];

    auto g_yz_yz_xy_xz = buffer_dddd[1016];

    auto g_yz_yz_xy_yy = buffer_dddd[1017];

    auto g_yz_yz_xy_yz = buffer_dddd[1018];

    auto g_yz_yz_xy_zz = buffer_dddd[1019];

    auto g_yz_yz_xz_xx = buffer_dddd[1020];

    auto g_yz_yz_xz_xy = buffer_dddd[1021];

    auto g_yz_yz_xz_xz = buffer_dddd[1022];

    auto g_yz_yz_xz_yy = buffer_dddd[1023];

    auto g_yz_yz_xz_yz = buffer_dddd[1024];

    auto g_yz_yz_xz_zz = buffer_dddd[1025];

    auto g_yz_yz_yy_xx = buffer_dddd[1026];

    auto g_yz_yz_yy_xy = buffer_dddd[1027];

    auto g_yz_yz_yy_xz = buffer_dddd[1028];

    auto g_yz_yz_yy_yy = buffer_dddd[1029];

    auto g_yz_yz_yy_yz = buffer_dddd[1030];

    auto g_yz_yz_yy_zz = buffer_dddd[1031];

    auto g_yz_yz_yz_xx = buffer_dddd[1032];

    auto g_yz_yz_yz_xy = buffer_dddd[1033];

    auto g_yz_yz_yz_xz = buffer_dddd[1034];

    auto g_yz_yz_yz_yy = buffer_dddd[1035];

    auto g_yz_yz_yz_yz = buffer_dddd[1036];

    auto g_yz_yz_yz_zz = buffer_dddd[1037];

    auto g_yz_yz_zz_xx = buffer_dddd[1038];

    auto g_yz_yz_zz_xy = buffer_dddd[1039];

    auto g_yz_yz_zz_xz = buffer_dddd[1040];

    auto g_yz_yz_zz_yy = buffer_dddd[1041];

    auto g_yz_yz_zz_yz = buffer_dddd[1042];

    auto g_yz_yz_zz_zz = buffer_dddd[1043];

    auto g_yz_zz_xx_xx = buffer_dddd[1044];

    auto g_yz_zz_xx_xy = buffer_dddd[1045];

    auto g_yz_zz_xx_xz = buffer_dddd[1046];

    auto g_yz_zz_xx_yy = buffer_dddd[1047];

    auto g_yz_zz_xx_yz = buffer_dddd[1048];

    auto g_yz_zz_xx_zz = buffer_dddd[1049];

    auto g_yz_zz_xy_xx = buffer_dddd[1050];

    auto g_yz_zz_xy_xy = buffer_dddd[1051];

    auto g_yz_zz_xy_xz = buffer_dddd[1052];

    auto g_yz_zz_xy_yy = buffer_dddd[1053];

    auto g_yz_zz_xy_yz = buffer_dddd[1054];

    auto g_yz_zz_xy_zz = buffer_dddd[1055];

    auto g_yz_zz_xz_xx = buffer_dddd[1056];

    auto g_yz_zz_xz_xy = buffer_dddd[1057];

    auto g_yz_zz_xz_xz = buffer_dddd[1058];

    auto g_yz_zz_xz_yy = buffer_dddd[1059];

    auto g_yz_zz_xz_yz = buffer_dddd[1060];

    auto g_yz_zz_xz_zz = buffer_dddd[1061];

    auto g_yz_zz_yy_xx = buffer_dddd[1062];

    auto g_yz_zz_yy_xy = buffer_dddd[1063];

    auto g_yz_zz_yy_xz = buffer_dddd[1064];

    auto g_yz_zz_yy_yy = buffer_dddd[1065];

    auto g_yz_zz_yy_yz = buffer_dddd[1066];

    auto g_yz_zz_yy_zz = buffer_dddd[1067];

    auto g_yz_zz_yz_xx = buffer_dddd[1068];

    auto g_yz_zz_yz_xy = buffer_dddd[1069];

    auto g_yz_zz_yz_xz = buffer_dddd[1070];

    auto g_yz_zz_yz_yy = buffer_dddd[1071];

    auto g_yz_zz_yz_yz = buffer_dddd[1072];

    auto g_yz_zz_yz_zz = buffer_dddd[1073];

    auto g_yz_zz_zz_xx = buffer_dddd[1074];

    auto g_yz_zz_zz_xy = buffer_dddd[1075];

    auto g_yz_zz_zz_xz = buffer_dddd[1076];

    auto g_yz_zz_zz_yy = buffer_dddd[1077];

    auto g_yz_zz_zz_yz = buffer_dddd[1078];

    auto g_yz_zz_zz_zz = buffer_dddd[1079];

    auto g_zz_xx_xx_xx = buffer_dddd[1080];

    auto g_zz_xx_xx_xy = buffer_dddd[1081];

    auto g_zz_xx_xx_xz = buffer_dddd[1082];

    auto g_zz_xx_xx_yy = buffer_dddd[1083];

    auto g_zz_xx_xx_yz = buffer_dddd[1084];

    auto g_zz_xx_xx_zz = buffer_dddd[1085];

    auto g_zz_xx_xy_xx = buffer_dddd[1086];

    auto g_zz_xx_xy_xy = buffer_dddd[1087];

    auto g_zz_xx_xy_xz = buffer_dddd[1088];

    auto g_zz_xx_xy_yy = buffer_dddd[1089];

    auto g_zz_xx_xy_yz = buffer_dddd[1090];

    auto g_zz_xx_xy_zz = buffer_dddd[1091];

    auto g_zz_xx_xz_xx = buffer_dddd[1092];

    auto g_zz_xx_xz_xy = buffer_dddd[1093];

    auto g_zz_xx_xz_xz = buffer_dddd[1094];

    auto g_zz_xx_xz_yy = buffer_dddd[1095];

    auto g_zz_xx_xz_yz = buffer_dddd[1096];

    auto g_zz_xx_xz_zz = buffer_dddd[1097];

    auto g_zz_xx_yy_xx = buffer_dddd[1098];

    auto g_zz_xx_yy_xy = buffer_dddd[1099];

    auto g_zz_xx_yy_xz = buffer_dddd[1100];

    auto g_zz_xx_yy_yy = buffer_dddd[1101];

    auto g_zz_xx_yy_yz = buffer_dddd[1102];

    auto g_zz_xx_yy_zz = buffer_dddd[1103];

    auto g_zz_xx_yz_xx = buffer_dddd[1104];

    auto g_zz_xx_yz_xy = buffer_dddd[1105];

    auto g_zz_xx_yz_xz = buffer_dddd[1106];

    auto g_zz_xx_yz_yy = buffer_dddd[1107];

    auto g_zz_xx_yz_yz = buffer_dddd[1108];

    auto g_zz_xx_yz_zz = buffer_dddd[1109];

    auto g_zz_xx_zz_xx = buffer_dddd[1110];

    auto g_zz_xx_zz_xy = buffer_dddd[1111];

    auto g_zz_xx_zz_xz = buffer_dddd[1112];

    auto g_zz_xx_zz_yy = buffer_dddd[1113];

    auto g_zz_xx_zz_yz = buffer_dddd[1114];

    auto g_zz_xx_zz_zz = buffer_dddd[1115];

    auto g_zz_xy_xx_xx = buffer_dddd[1116];

    auto g_zz_xy_xx_xy = buffer_dddd[1117];

    auto g_zz_xy_xx_xz = buffer_dddd[1118];

    auto g_zz_xy_xx_yy = buffer_dddd[1119];

    auto g_zz_xy_xx_yz = buffer_dddd[1120];

    auto g_zz_xy_xx_zz = buffer_dddd[1121];

    auto g_zz_xy_xy_xx = buffer_dddd[1122];

    auto g_zz_xy_xy_xy = buffer_dddd[1123];

    auto g_zz_xy_xy_xz = buffer_dddd[1124];

    auto g_zz_xy_xy_yy = buffer_dddd[1125];

    auto g_zz_xy_xy_yz = buffer_dddd[1126];

    auto g_zz_xy_xy_zz = buffer_dddd[1127];

    auto g_zz_xy_xz_xx = buffer_dddd[1128];

    auto g_zz_xy_xz_xy = buffer_dddd[1129];

    auto g_zz_xy_xz_xz = buffer_dddd[1130];

    auto g_zz_xy_xz_yy = buffer_dddd[1131];

    auto g_zz_xy_xz_yz = buffer_dddd[1132];

    auto g_zz_xy_xz_zz = buffer_dddd[1133];

    auto g_zz_xy_yy_xx = buffer_dddd[1134];

    auto g_zz_xy_yy_xy = buffer_dddd[1135];

    auto g_zz_xy_yy_xz = buffer_dddd[1136];

    auto g_zz_xy_yy_yy = buffer_dddd[1137];

    auto g_zz_xy_yy_yz = buffer_dddd[1138];

    auto g_zz_xy_yy_zz = buffer_dddd[1139];

    auto g_zz_xy_yz_xx = buffer_dddd[1140];

    auto g_zz_xy_yz_xy = buffer_dddd[1141];

    auto g_zz_xy_yz_xz = buffer_dddd[1142];

    auto g_zz_xy_yz_yy = buffer_dddd[1143];

    auto g_zz_xy_yz_yz = buffer_dddd[1144];

    auto g_zz_xy_yz_zz = buffer_dddd[1145];

    auto g_zz_xy_zz_xx = buffer_dddd[1146];

    auto g_zz_xy_zz_xy = buffer_dddd[1147];

    auto g_zz_xy_zz_xz = buffer_dddd[1148];

    auto g_zz_xy_zz_yy = buffer_dddd[1149];

    auto g_zz_xy_zz_yz = buffer_dddd[1150];

    auto g_zz_xy_zz_zz = buffer_dddd[1151];

    auto g_zz_xz_xx_xx = buffer_dddd[1152];

    auto g_zz_xz_xx_xy = buffer_dddd[1153];

    auto g_zz_xz_xx_xz = buffer_dddd[1154];

    auto g_zz_xz_xx_yy = buffer_dddd[1155];

    auto g_zz_xz_xx_yz = buffer_dddd[1156];

    auto g_zz_xz_xx_zz = buffer_dddd[1157];

    auto g_zz_xz_xy_xx = buffer_dddd[1158];

    auto g_zz_xz_xy_xy = buffer_dddd[1159];

    auto g_zz_xz_xy_xz = buffer_dddd[1160];

    auto g_zz_xz_xy_yy = buffer_dddd[1161];

    auto g_zz_xz_xy_yz = buffer_dddd[1162];

    auto g_zz_xz_xy_zz = buffer_dddd[1163];

    auto g_zz_xz_xz_xx = buffer_dddd[1164];

    auto g_zz_xz_xz_xy = buffer_dddd[1165];

    auto g_zz_xz_xz_xz = buffer_dddd[1166];

    auto g_zz_xz_xz_yy = buffer_dddd[1167];

    auto g_zz_xz_xz_yz = buffer_dddd[1168];

    auto g_zz_xz_xz_zz = buffer_dddd[1169];

    auto g_zz_xz_yy_xx = buffer_dddd[1170];

    auto g_zz_xz_yy_xy = buffer_dddd[1171];

    auto g_zz_xz_yy_xz = buffer_dddd[1172];

    auto g_zz_xz_yy_yy = buffer_dddd[1173];

    auto g_zz_xz_yy_yz = buffer_dddd[1174];

    auto g_zz_xz_yy_zz = buffer_dddd[1175];

    auto g_zz_xz_yz_xx = buffer_dddd[1176];

    auto g_zz_xz_yz_xy = buffer_dddd[1177];

    auto g_zz_xz_yz_xz = buffer_dddd[1178];

    auto g_zz_xz_yz_yy = buffer_dddd[1179];

    auto g_zz_xz_yz_yz = buffer_dddd[1180];

    auto g_zz_xz_yz_zz = buffer_dddd[1181];

    auto g_zz_xz_zz_xx = buffer_dddd[1182];

    auto g_zz_xz_zz_xy = buffer_dddd[1183];

    auto g_zz_xz_zz_xz = buffer_dddd[1184];

    auto g_zz_xz_zz_yy = buffer_dddd[1185];

    auto g_zz_xz_zz_yz = buffer_dddd[1186];

    auto g_zz_xz_zz_zz = buffer_dddd[1187];

    auto g_zz_yy_xx_xx = buffer_dddd[1188];

    auto g_zz_yy_xx_xy = buffer_dddd[1189];

    auto g_zz_yy_xx_xz = buffer_dddd[1190];

    auto g_zz_yy_xx_yy = buffer_dddd[1191];

    auto g_zz_yy_xx_yz = buffer_dddd[1192];

    auto g_zz_yy_xx_zz = buffer_dddd[1193];

    auto g_zz_yy_xy_xx = buffer_dddd[1194];

    auto g_zz_yy_xy_xy = buffer_dddd[1195];

    auto g_zz_yy_xy_xz = buffer_dddd[1196];

    auto g_zz_yy_xy_yy = buffer_dddd[1197];

    auto g_zz_yy_xy_yz = buffer_dddd[1198];

    auto g_zz_yy_xy_zz = buffer_dddd[1199];

    auto g_zz_yy_xz_xx = buffer_dddd[1200];

    auto g_zz_yy_xz_xy = buffer_dddd[1201];

    auto g_zz_yy_xz_xz = buffer_dddd[1202];

    auto g_zz_yy_xz_yy = buffer_dddd[1203];

    auto g_zz_yy_xz_yz = buffer_dddd[1204];

    auto g_zz_yy_xz_zz = buffer_dddd[1205];

    auto g_zz_yy_yy_xx = buffer_dddd[1206];

    auto g_zz_yy_yy_xy = buffer_dddd[1207];

    auto g_zz_yy_yy_xz = buffer_dddd[1208];

    auto g_zz_yy_yy_yy = buffer_dddd[1209];

    auto g_zz_yy_yy_yz = buffer_dddd[1210];

    auto g_zz_yy_yy_zz = buffer_dddd[1211];

    auto g_zz_yy_yz_xx = buffer_dddd[1212];

    auto g_zz_yy_yz_xy = buffer_dddd[1213];

    auto g_zz_yy_yz_xz = buffer_dddd[1214];

    auto g_zz_yy_yz_yy = buffer_dddd[1215];

    auto g_zz_yy_yz_yz = buffer_dddd[1216];

    auto g_zz_yy_yz_zz = buffer_dddd[1217];

    auto g_zz_yy_zz_xx = buffer_dddd[1218];

    auto g_zz_yy_zz_xy = buffer_dddd[1219];

    auto g_zz_yy_zz_xz = buffer_dddd[1220];

    auto g_zz_yy_zz_yy = buffer_dddd[1221];

    auto g_zz_yy_zz_yz = buffer_dddd[1222];

    auto g_zz_yy_zz_zz = buffer_dddd[1223];

    auto g_zz_yz_xx_xx = buffer_dddd[1224];

    auto g_zz_yz_xx_xy = buffer_dddd[1225];

    auto g_zz_yz_xx_xz = buffer_dddd[1226];

    auto g_zz_yz_xx_yy = buffer_dddd[1227];

    auto g_zz_yz_xx_yz = buffer_dddd[1228];

    auto g_zz_yz_xx_zz = buffer_dddd[1229];

    auto g_zz_yz_xy_xx = buffer_dddd[1230];

    auto g_zz_yz_xy_xy = buffer_dddd[1231];

    auto g_zz_yz_xy_xz = buffer_dddd[1232];

    auto g_zz_yz_xy_yy = buffer_dddd[1233];

    auto g_zz_yz_xy_yz = buffer_dddd[1234];

    auto g_zz_yz_xy_zz = buffer_dddd[1235];

    auto g_zz_yz_xz_xx = buffer_dddd[1236];

    auto g_zz_yz_xz_xy = buffer_dddd[1237];

    auto g_zz_yz_xz_xz = buffer_dddd[1238];

    auto g_zz_yz_xz_yy = buffer_dddd[1239];

    auto g_zz_yz_xz_yz = buffer_dddd[1240];

    auto g_zz_yz_xz_zz = buffer_dddd[1241];

    auto g_zz_yz_yy_xx = buffer_dddd[1242];

    auto g_zz_yz_yy_xy = buffer_dddd[1243];

    auto g_zz_yz_yy_xz = buffer_dddd[1244];

    auto g_zz_yz_yy_yy = buffer_dddd[1245];

    auto g_zz_yz_yy_yz = buffer_dddd[1246];

    auto g_zz_yz_yy_zz = buffer_dddd[1247];

    auto g_zz_yz_yz_xx = buffer_dddd[1248];

    auto g_zz_yz_yz_xy = buffer_dddd[1249];

    auto g_zz_yz_yz_xz = buffer_dddd[1250];

    auto g_zz_yz_yz_yy = buffer_dddd[1251];

    auto g_zz_yz_yz_yz = buffer_dddd[1252];

    auto g_zz_yz_yz_zz = buffer_dddd[1253];

    auto g_zz_yz_zz_xx = buffer_dddd[1254];

    auto g_zz_yz_zz_xy = buffer_dddd[1255];

    auto g_zz_yz_zz_xz = buffer_dddd[1256];

    auto g_zz_yz_zz_yy = buffer_dddd[1257];

    auto g_zz_yz_zz_yz = buffer_dddd[1258];

    auto g_zz_yz_zz_zz = buffer_dddd[1259];

    auto g_zz_zz_xx_xx = buffer_dddd[1260];

    auto g_zz_zz_xx_xy = buffer_dddd[1261];

    auto g_zz_zz_xx_xz = buffer_dddd[1262];

    auto g_zz_zz_xx_yy = buffer_dddd[1263];

    auto g_zz_zz_xx_yz = buffer_dddd[1264];

    auto g_zz_zz_xx_zz = buffer_dddd[1265];

    auto g_zz_zz_xy_xx = buffer_dddd[1266];

    auto g_zz_zz_xy_xy = buffer_dddd[1267];

    auto g_zz_zz_xy_xz = buffer_dddd[1268];

    auto g_zz_zz_xy_yy = buffer_dddd[1269];

    auto g_zz_zz_xy_yz = buffer_dddd[1270];

    auto g_zz_zz_xy_zz = buffer_dddd[1271];

    auto g_zz_zz_xz_xx = buffer_dddd[1272];

    auto g_zz_zz_xz_xy = buffer_dddd[1273];

    auto g_zz_zz_xz_xz = buffer_dddd[1274];

    auto g_zz_zz_xz_yy = buffer_dddd[1275];

    auto g_zz_zz_xz_yz = buffer_dddd[1276];

    auto g_zz_zz_xz_zz = buffer_dddd[1277];

    auto g_zz_zz_yy_xx = buffer_dddd[1278];

    auto g_zz_zz_yy_xy = buffer_dddd[1279];

    auto g_zz_zz_yy_xz = buffer_dddd[1280];

    auto g_zz_zz_yy_yy = buffer_dddd[1281];

    auto g_zz_zz_yy_yz = buffer_dddd[1282];

    auto g_zz_zz_yy_zz = buffer_dddd[1283];

    auto g_zz_zz_yz_xx = buffer_dddd[1284];

    auto g_zz_zz_yz_xy = buffer_dddd[1285];

    auto g_zz_zz_yz_xz = buffer_dddd[1286];

    auto g_zz_zz_yz_yy = buffer_dddd[1287];

    auto g_zz_zz_yz_yz = buffer_dddd[1288];

    auto g_zz_zz_yz_zz = buffer_dddd[1289];

    auto g_zz_zz_zz_xx = buffer_dddd[1290];

    auto g_zz_zz_zz_xy = buffer_dddd[1291];

    auto g_zz_zz_zz_xz = buffer_dddd[1292];

    auto g_zz_zz_zz_yy = buffer_dddd[1293];

    auto g_zz_zz_zz_yz = buffer_dddd[1294];

    auto g_zz_zz_zz_zz = buffer_dddd[1295];

    /// Set up components of integrals buffer : buffer_1000_pddd

    auto g_x_0_0_0_x_xx_xx_xx = buffer_1000_pddd[0];

    auto g_x_0_0_0_x_xx_xx_xy = buffer_1000_pddd[1];

    auto g_x_0_0_0_x_xx_xx_xz = buffer_1000_pddd[2];

    auto g_x_0_0_0_x_xx_xx_yy = buffer_1000_pddd[3];

    auto g_x_0_0_0_x_xx_xx_yz = buffer_1000_pddd[4];

    auto g_x_0_0_0_x_xx_xx_zz = buffer_1000_pddd[5];

    auto g_x_0_0_0_x_xx_xy_xx = buffer_1000_pddd[6];

    auto g_x_0_0_0_x_xx_xy_xy = buffer_1000_pddd[7];

    auto g_x_0_0_0_x_xx_xy_xz = buffer_1000_pddd[8];

    auto g_x_0_0_0_x_xx_xy_yy = buffer_1000_pddd[9];

    auto g_x_0_0_0_x_xx_xy_yz = buffer_1000_pddd[10];

    auto g_x_0_0_0_x_xx_xy_zz = buffer_1000_pddd[11];

    auto g_x_0_0_0_x_xx_xz_xx = buffer_1000_pddd[12];

    auto g_x_0_0_0_x_xx_xz_xy = buffer_1000_pddd[13];

    auto g_x_0_0_0_x_xx_xz_xz = buffer_1000_pddd[14];

    auto g_x_0_0_0_x_xx_xz_yy = buffer_1000_pddd[15];

    auto g_x_0_0_0_x_xx_xz_yz = buffer_1000_pddd[16];

    auto g_x_0_0_0_x_xx_xz_zz = buffer_1000_pddd[17];

    auto g_x_0_0_0_x_xx_yy_xx = buffer_1000_pddd[18];

    auto g_x_0_0_0_x_xx_yy_xy = buffer_1000_pddd[19];

    auto g_x_0_0_0_x_xx_yy_xz = buffer_1000_pddd[20];

    auto g_x_0_0_0_x_xx_yy_yy = buffer_1000_pddd[21];

    auto g_x_0_0_0_x_xx_yy_yz = buffer_1000_pddd[22];

    auto g_x_0_0_0_x_xx_yy_zz = buffer_1000_pddd[23];

    auto g_x_0_0_0_x_xx_yz_xx = buffer_1000_pddd[24];

    auto g_x_0_0_0_x_xx_yz_xy = buffer_1000_pddd[25];

    auto g_x_0_0_0_x_xx_yz_xz = buffer_1000_pddd[26];

    auto g_x_0_0_0_x_xx_yz_yy = buffer_1000_pddd[27];

    auto g_x_0_0_0_x_xx_yz_yz = buffer_1000_pddd[28];

    auto g_x_0_0_0_x_xx_yz_zz = buffer_1000_pddd[29];

    auto g_x_0_0_0_x_xx_zz_xx = buffer_1000_pddd[30];

    auto g_x_0_0_0_x_xx_zz_xy = buffer_1000_pddd[31];

    auto g_x_0_0_0_x_xx_zz_xz = buffer_1000_pddd[32];

    auto g_x_0_0_0_x_xx_zz_yy = buffer_1000_pddd[33];

    auto g_x_0_0_0_x_xx_zz_yz = buffer_1000_pddd[34];

    auto g_x_0_0_0_x_xx_zz_zz = buffer_1000_pddd[35];

    auto g_x_0_0_0_x_xy_xx_xx = buffer_1000_pddd[36];

    auto g_x_0_0_0_x_xy_xx_xy = buffer_1000_pddd[37];

    auto g_x_0_0_0_x_xy_xx_xz = buffer_1000_pddd[38];

    auto g_x_0_0_0_x_xy_xx_yy = buffer_1000_pddd[39];

    auto g_x_0_0_0_x_xy_xx_yz = buffer_1000_pddd[40];

    auto g_x_0_0_0_x_xy_xx_zz = buffer_1000_pddd[41];

    auto g_x_0_0_0_x_xy_xy_xx = buffer_1000_pddd[42];

    auto g_x_0_0_0_x_xy_xy_xy = buffer_1000_pddd[43];

    auto g_x_0_0_0_x_xy_xy_xz = buffer_1000_pddd[44];

    auto g_x_0_0_0_x_xy_xy_yy = buffer_1000_pddd[45];

    auto g_x_0_0_0_x_xy_xy_yz = buffer_1000_pddd[46];

    auto g_x_0_0_0_x_xy_xy_zz = buffer_1000_pddd[47];

    auto g_x_0_0_0_x_xy_xz_xx = buffer_1000_pddd[48];

    auto g_x_0_0_0_x_xy_xz_xy = buffer_1000_pddd[49];

    auto g_x_0_0_0_x_xy_xz_xz = buffer_1000_pddd[50];

    auto g_x_0_0_0_x_xy_xz_yy = buffer_1000_pddd[51];

    auto g_x_0_0_0_x_xy_xz_yz = buffer_1000_pddd[52];

    auto g_x_0_0_0_x_xy_xz_zz = buffer_1000_pddd[53];

    auto g_x_0_0_0_x_xy_yy_xx = buffer_1000_pddd[54];

    auto g_x_0_0_0_x_xy_yy_xy = buffer_1000_pddd[55];

    auto g_x_0_0_0_x_xy_yy_xz = buffer_1000_pddd[56];

    auto g_x_0_0_0_x_xy_yy_yy = buffer_1000_pddd[57];

    auto g_x_0_0_0_x_xy_yy_yz = buffer_1000_pddd[58];

    auto g_x_0_0_0_x_xy_yy_zz = buffer_1000_pddd[59];

    auto g_x_0_0_0_x_xy_yz_xx = buffer_1000_pddd[60];

    auto g_x_0_0_0_x_xy_yz_xy = buffer_1000_pddd[61];

    auto g_x_0_0_0_x_xy_yz_xz = buffer_1000_pddd[62];

    auto g_x_0_0_0_x_xy_yz_yy = buffer_1000_pddd[63];

    auto g_x_0_0_0_x_xy_yz_yz = buffer_1000_pddd[64];

    auto g_x_0_0_0_x_xy_yz_zz = buffer_1000_pddd[65];

    auto g_x_0_0_0_x_xy_zz_xx = buffer_1000_pddd[66];

    auto g_x_0_0_0_x_xy_zz_xy = buffer_1000_pddd[67];

    auto g_x_0_0_0_x_xy_zz_xz = buffer_1000_pddd[68];

    auto g_x_0_0_0_x_xy_zz_yy = buffer_1000_pddd[69];

    auto g_x_0_0_0_x_xy_zz_yz = buffer_1000_pddd[70];

    auto g_x_0_0_0_x_xy_zz_zz = buffer_1000_pddd[71];

    auto g_x_0_0_0_x_xz_xx_xx = buffer_1000_pddd[72];

    auto g_x_0_0_0_x_xz_xx_xy = buffer_1000_pddd[73];

    auto g_x_0_0_0_x_xz_xx_xz = buffer_1000_pddd[74];

    auto g_x_0_0_0_x_xz_xx_yy = buffer_1000_pddd[75];

    auto g_x_0_0_0_x_xz_xx_yz = buffer_1000_pddd[76];

    auto g_x_0_0_0_x_xz_xx_zz = buffer_1000_pddd[77];

    auto g_x_0_0_0_x_xz_xy_xx = buffer_1000_pddd[78];

    auto g_x_0_0_0_x_xz_xy_xy = buffer_1000_pddd[79];

    auto g_x_0_0_0_x_xz_xy_xz = buffer_1000_pddd[80];

    auto g_x_0_0_0_x_xz_xy_yy = buffer_1000_pddd[81];

    auto g_x_0_0_0_x_xz_xy_yz = buffer_1000_pddd[82];

    auto g_x_0_0_0_x_xz_xy_zz = buffer_1000_pddd[83];

    auto g_x_0_0_0_x_xz_xz_xx = buffer_1000_pddd[84];

    auto g_x_0_0_0_x_xz_xz_xy = buffer_1000_pddd[85];

    auto g_x_0_0_0_x_xz_xz_xz = buffer_1000_pddd[86];

    auto g_x_0_0_0_x_xz_xz_yy = buffer_1000_pddd[87];

    auto g_x_0_0_0_x_xz_xz_yz = buffer_1000_pddd[88];

    auto g_x_0_0_0_x_xz_xz_zz = buffer_1000_pddd[89];

    auto g_x_0_0_0_x_xz_yy_xx = buffer_1000_pddd[90];

    auto g_x_0_0_0_x_xz_yy_xy = buffer_1000_pddd[91];

    auto g_x_0_0_0_x_xz_yy_xz = buffer_1000_pddd[92];

    auto g_x_0_0_0_x_xz_yy_yy = buffer_1000_pddd[93];

    auto g_x_0_0_0_x_xz_yy_yz = buffer_1000_pddd[94];

    auto g_x_0_0_0_x_xz_yy_zz = buffer_1000_pddd[95];

    auto g_x_0_0_0_x_xz_yz_xx = buffer_1000_pddd[96];

    auto g_x_0_0_0_x_xz_yz_xy = buffer_1000_pddd[97];

    auto g_x_0_0_0_x_xz_yz_xz = buffer_1000_pddd[98];

    auto g_x_0_0_0_x_xz_yz_yy = buffer_1000_pddd[99];

    auto g_x_0_0_0_x_xz_yz_yz = buffer_1000_pddd[100];

    auto g_x_0_0_0_x_xz_yz_zz = buffer_1000_pddd[101];

    auto g_x_0_0_0_x_xz_zz_xx = buffer_1000_pddd[102];

    auto g_x_0_0_0_x_xz_zz_xy = buffer_1000_pddd[103];

    auto g_x_0_0_0_x_xz_zz_xz = buffer_1000_pddd[104];

    auto g_x_0_0_0_x_xz_zz_yy = buffer_1000_pddd[105];

    auto g_x_0_0_0_x_xz_zz_yz = buffer_1000_pddd[106];

    auto g_x_0_0_0_x_xz_zz_zz = buffer_1000_pddd[107];

    auto g_x_0_0_0_x_yy_xx_xx = buffer_1000_pddd[108];

    auto g_x_0_0_0_x_yy_xx_xy = buffer_1000_pddd[109];

    auto g_x_0_0_0_x_yy_xx_xz = buffer_1000_pddd[110];

    auto g_x_0_0_0_x_yy_xx_yy = buffer_1000_pddd[111];

    auto g_x_0_0_0_x_yy_xx_yz = buffer_1000_pddd[112];

    auto g_x_0_0_0_x_yy_xx_zz = buffer_1000_pddd[113];

    auto g_x_0_0_0_x_yy_xy_xx = buffer_1000_pddd[114];

    auto g_x_0_0_0_x_yy_xy_xy = buffer_1000_pddd[115];

    auto g_x_0_0_0_x_yy_xy_xz = buffer_1000_pddd[116];

    auto g_x_0_0_0_x_yy_xy_yy = buffer_1000_pddd[117];

    auto g_x_0_0_0_x_yy_xy_yz = buffer_1000_pddd[118];

    auto g_x_0_0_0_x_yy_xy_zz = buffer_1000_pddd[119];

    auto g_x_0_0_0_x_yy_xz_xx = buffer_1000_pddd[120];

    auto g_x_0_0_0_x_yy_xz_xy = buffer_1000_pddd[121];

    auto g_x_0_0_0_x_yy_xz_xz = buffer_1000_pddd[122];

    auto g_x_0_0_0_x_yy_xz_yy = buffer_1000_pddd[123];

    auto g_x_0_0_0_x_yy_xz_yz = buffer_1000_pddd[124];

    auto g_x_0_0_0_x_yy_xz_zz = buffer_1000_pddd[125];

    auto g_x_0_0_0_x_yy_yy_xx = buffer_1000_pddd[126];

    auto g_x_0_0_0_x_yy_yy_xy = buffer_1000_pddd[127];

    auto g_x_0_0_0_x_yy_yy_xz = buffer_1000_pddd[128];

    auto g_x_0_0_0_x_yy_yy_yy = buffer_1000_pddd[129];

    auto g_x_0_0_0_x_yy_yy_yz = buffer_1000_pddd[130];

    auto g_x_0_0_0_x_yy_yy_zz = buffer_1000_pddd[131];

    auto g_x_0_0_0_x_yy_yz_xx = buffer_1000_pddd[132];

    auto g_x_0_0_0_x_yy_yz_xy = buffer_1000_pddd[133];

    auto g_x_0_0_0_x_yy_yz_xz = buffer_1000_pddd[134];

    auto g_x_0_0_0_x_yy_yz_yy = buffer_1000_pddd[135];

    auto g_x_0_0_0_x_yy_yz_yz = buffer_1000_pddd[136];

    auto g_x_0_0_0_x_yy_yz_zz = buffer_1000_pddd[137];

    auto g_x_0_0_0_x_yy_zz_xx = buffer_1000_pddd[138];

    auto g_x_0_0_0_x_yy_zz_xy = buffer_1000_pddd[139];

    auto g_x_0_0_0_x_yy_zz_xz = buffer_1000_pddd[140];

    auto g_x_0_0_0_x_yy_zz_yy = buffer_1000_pddd[141];

    auto g_x_0_0_0_x_yy_zz_yz = buffer_1000_pddd[142];

    auto g_x_0_0_0_x_yy_zz_zz = buffer_1000_pddd[143];

    auto g_x_0_0_0_x_yz_xx_xx = buffer_1000_pddd[144];

    auto g_x_0_0_0_x_yz_xx_xy = buffer_1000_pddd[145];

    auto g_x_0_0_0_x_yz_xx_xz = buffer_1000_pddd[146];

    auto g_x_0_0_0_x_yz_xx_yy = buffer_1000_pddd[147];

    auto g_x_0_0_0_x_yz_xx_yz = buffer_1000_pddd[148];

    auto g_x_0_0_0_x_yz_xx_zz = buffer_1000_pddd[149];

    auto g_x_0_0_0_x_yz_xy_xx = buffer_1000_pddd[150];

    auto g_x_0_0_0_x_yz_xy_xy = buffer_1000_pddd[151];

    auto g_x_0_0_0_x_yz_xy_xz = buffer_1000_pddd[152];

    auto g_x_0_0_0_x_yz_xy_yy = buffer_1000_pddd[153];

    auto g_x_0_0_0_x_yz_xy_yz = buffer_1000_pddd[154];

    auto g_x_0_0_0_x_yz_xy_zz = buffer_1000_pddd[155];

    auto g_x_0_0_0_x_yz_xz_xx = buffer_1000_pddd[156];

    auto g_x_0_0_0_x_yz_xz_xy = buffer_1000_pddd[157];

    auto g_x_0_0_0_x_yz_xz_xz = buffer_1000_pddd[158];

    auto g_x_0_0_0_x_yz_xz_yy = buffer_1000_pddd[159];

    auto g_x_0_0_0_x_yz_xz_yz = buffer_1000_pddd[160];

    auto g_x_0_0_0_x_yz_xz_zz = buffer_1000_pddd[161];

    auto g_x_0_0_0_x_yz_yy_xx = buffer_1000_pddd[162];

    auto g_x_0_0_0_x_yz_yy_xy = buffer_1000_pddd[163];

    auto g_x_0_0_0_x_yz_yy_xz = buffer_1000_pddd[164];

    auto g_x_0_0_0_x_yz_yy_yy = buffer_1000_pddd[165];

    auto g_x_0_0_0_x_yz_yy_yz = buffer_1000_pddd[166];

    auto g_x_0_0_0_x_yz_yy_zz = buffer_1000_pddd[167];

    auto g_x_0_0_0_x_yz_yz_xx = buffer_1000_pddd[168];

    auto g_x_0_0_0_x_yz_yz_xy = buffer_1000_pddd[169];

    auto g_x_0_0_0_x_yz_yz_xz = buffer_1000_pddd[170];

    auto g_x_0_0_0_x_yz_yz_yy = buffer_1000_pddd[171];

    auto g_x_0_0_0_x_yz_yz_yz = buffer_1000_pddd[172];

    auto g_x_0_0_0_x_yz_yz_zz = buffer_1000_pddd[173];

    auto g_x_0_0_0_x_yz_zz_xx = buffer_1000_pddd[174];

    auto g_x_0_0_0_x_yz_zz_xy = buffer_1000_pddd[175];

    auto g_x_0_0_0_x_yz_zz_xz = buffer_1000_pddd[176];

    auto g_x_0_0_0_x_yz_zz_yy = buffer_1000_pddd[177];

    auto g_x_0_0_0_x_yz_zz_yz = buffer_1000_pddd[178];

    auto g_x_0_0_0_x_yz_zz_zz = buffer_1000_pddd[179];

    auto g_x_0_0_0_x_zz_xx_xx = buffer_1000_pddd[180];

    auto g_x_0_0_0_x_zz_xx_xy = buffer_1000_pddd[181];

    auto g_x_0_0_0_x_zz_xx_xz = buffer_1000_pddd[182];

    auto g_x_0_0_0_x_zz_xx_yy = buffer_1000_pddd[183];

    auto g_x_0_0_0_x_zz_xx_yz = buffer_1000_pddd[184];

    auto g_x_0_0_0_x_zz_xx_zz = buffer_1000_pddd[185];

    auto g_x_0_0_0_x_zz_xy_xx = buffer_1000_pddd[186];

    auto g_x_0_0_0_x_zz_xy_xy = buffer_1000_pddd[187];

    auto g_x_0_0_0_x_zz_xy_xz = buffer_1000_pddd[188];

    auto g_x_0_0_0_x_zz_xy_yy = buffer_1000_pddd[189];

    auto g_x_0_0_0_x_zz_xy_yz = buffer_1000_pddd[190];

    auto g_x_0_0_0_x_zz_xy_zz = buffer_1000_pddd[191];

    auto g_x_0_0_0_x_zz_xz_xx = buffer_1000_pddd[192];

    auto g_x_0_0_0_x_zz_xz_xy = buffer_1000_pddd[193];

    auto g_x_0_0_0_x_zz_xz_xz = buffer_1000_pddd[194];

    auto g_x_0_0_0_x_zz_xz_yy = buffer_1000_pddd[195];

    auto g_x_0_0_0_x_zz_xz_yz = buffer_1000_pddd[196];

    auto g_x_0_0_0_x_zz_xz_zz = buffer_1000_pddd[197];

    auto g_x_0_0_0_x_zz_yy_xx = buffer_1000_pddd[198];

    auto g_x_0_0_0_x_zz_yy_xy = buffer_1000_pddd[199];

    auto g_x_0_0_0_x_zz_yy_xz = buffer_1000_pddd[200];

    auto g_x_0_0_0_x_zz_yy_yy = buffer_1000_pddd[201];

    auto g_x_0_0_0_x_zz_yy_yz = buffer_1000_pddd[202];

    auto g_x_0_0_0_x_zz_yy_zz = buffer_1000_pddd[203];

    auto g_x_0_0_0_x_zz_yz_xx = buffer_1000_pddd[204];

    auto g_x_0_0_0_x_zz_yz_xy = buffer_1000_pddd[205];

    auto g_x_0_0_0_x_zz_yz_xz = buffer_1000_pddd[206];

    auto g_x_0_0_0_x_zz_yz_yy = buffer_1000_pddd[207];

    auto g_x_0_0_0_x_zz_yz_yz = buffer_1000_pddd[208];

    auto g_x_0_0_0_x_zz_yz_zz = buffer_1000_pddd[209];

    auto g_x_0_0_0_x_zz_zz_xx = buffer_1000_pddd[210];

    auto g_x_0_0_0_x_zz_zz_xy = buffer_1000_pddd[211];

    auto g_x_0_0_0_x_zz_zz_xz = buffer_1000_pddd[212];

    auto g_x_0_0_0_x_zz_zz_yy = buffer_1000_pddd[213];

    auto g_x_0_0_0_x_zz_zz_yz = buffer_1000_pddd[214];

    auto g_x_0_0_0_x_zz_zz_zz = buffer_1000_pddd[215];

    auto g_x_0_0_0_y_xx_xx_xx = buffer_1000_pddd[216];

    auto g_x_0_0_0_y_xx_xx_xy = buffer_1000_pddd[217];

    auto g_x_0_0_0_y_xx_xx_xz = buffer_1000_pddd[218];

    auto g_x_0_0_0_y_xx_xx_yy = buffer_1000_pddd[219];

    auto g_x_0_0_0_y_xx_xx_yz = buffer_1000_pddd[220];

    auto g_x_0_0_0_y_xx_xx_zz = buffer_1000_pddd[221];

    auto g_x_0_0_0_y_xx_xy_xx = buffer_1000_pddd[222];

    auto g_x_0_0_0_y_xx_xy_xy = buffer_1000_pddd[223];

    auto g_x_0_0_0_y_xx_xy_xz = buffer_1000_pddd[224];

    auto g_x_0_0_0_y_xx_xy_yy = buffer_1000_pddd[225];

    auto g_x_0_0_0_y_xx_xy_yz = buffer_1000_pddd[226];

    auto g_x_0_0_0_y_xx_xy_zz = buffer_1000_pddd[227];

    auto g_x_0_0_0_y_xx_xz_xx = buffer_1000_pddd[228];

    auto g_x_0_0_0_y_xx_xz_xy = buffer_1000_pddd[229];

    auto g_x_0_0_0_y_xx_xz_xz = buffer_1000_pddd[230];

    auto g_x_0_0_0_y_xx_xz_yy = buffer_1000_pddd[231];

    auto g_x_0_0_0_y_xx_xz_yz = buffer_1000_pddd[232];

    auto g_x_0_0_0_y_xx_xz_zz = buffer_1000_pddd[233];

    auto g_x_0_0_0_y_xx_yy_xx = buffer_1000_pddd[234];

    auto g_x_0_0_0_y_xx_yy_xy = buffer_1000_pddd[235];

    auto g_x_0_0_0_y_xx_yy_xz = buffer_1000_pddd[236];

    auto g_x_0_0_0_y_xx_yy_yy = buffer_1000_pddd[237];

    auto g_x_0_0_0_y_xx_yy_yz = buffer_1000_pddd[238];

    auto g_x_0_0_0_y_xx_yy_zz = buffer_1000_pddd[239];

    auto g_x_0_0_0_y_xx_yz_xx = buffer_1000_pddd[240];

    auto g_x_0_0_0_y_xx_yz_xy = buffer_1000_pddd[241];

    auto g_x_0_0_0_y_xx_yz_xz = buffer_1000_pddd[242];

    auto g_x_0_0_0_y_xx_yz_yy = buffer_1000_pddd[243];

    auto g_x_0_0_0_y_xx_yz_yz = buffer_1000_pddd[244];

    auto g_x_0_0_0_y_xx_yz_zz = buffer_1000_pddd[245];

    auto g_x_0_0_0_y_xx_zz_xx = buffer_1000_pddd[246];

    auto g_x_0_0_0_y_xx_zz_xy = buffer_1000_pddd[247];

    auto g_x_0_0_0_y_xx_zz_xz = buffer_1000_pddd[248];

    auto g_x_0_0_0_y_xx_zz_yy = buffer_1000_pddd[249];

    auto g_x_0_0_0_y_xx_zz_yz = buffer_1000_pddd[250];

    auto g_x_0_0_0_y_xx_zz_zz = buffer_1000_pddd[251];

    auto g_x_0_0_0_y_xy_xx_xx = buffer_1000_pddd[252];

    auto g_x_0_0_0_y_xy_xx_xy = buffer_1000_pddd[253];

    auto g_x_0_0_0_y_xy_xx_xz = buffer_1000_pddd[254];

    auto g_x_0_0_0_y_xy_xx_yy = buffer_1000_pddd[255];

    auto g_x_0_0_0_y_xy_xx_yz = buffer_1000_pddd[256];

    auto g_x_0_0_0_y_xy_xx_zz = buffer_1000_pddd[257];

    auto g_x_0_0_0_y_xy_xy_xx = buffer_1000_pddd[258];

    auto g_x_0_0_0_y_xy_xy_xy = buffer_1000_pddd[259];

    auto g_x_0_0_0_y_xy_xy_xz = buffer_1000_pddd[260];

    auto g_x_0_0_0_y_xy_xy_yy = buffer_1000_pddd[261];

    auto g_x_0_0_0_y_xy_xy_yz = buffer_1000_pddd[262];

    auto g_x_0_0_0_y_xy_xy_zz = buffer_1000_pddd[263];

    auto g_x_0_0_0_y_xy_xz_xx = buffer_1000_pddd[264];

    auto g_x_0_0_0_y_xy_xz_xy = buffer_1000_pddd[265];

    auto g_x_0_0_0_y_xy_xz_xz = buffer_1000_pddd[266];

    auto g_x_0_0_0_y_xy_xz_yy = buffer_1000_pddd[267];

    auto g_x_0_0_0_y_xy_xz_yz = buffer_1000_pddd[268];

    auto g_x_0_0_0_y_xy_xz_zz = buffer_1000_pddd[269];

    auto g_x_0_0_0_y_xy_yy_xx = buffer_1000_pddd[270];

    auto g_x_0_0_0_y_xy_yy_xy = buffer_1000_pddd[271];

    auto g_x_0_0_0_y_xy_yy_xz = buffer_1000_pddd[272];

    auto g_x_0_0_0_y_xy_yy_yy = buffer_1000_pddd[273];

    auto g_x_0_0_0_y_xy_yy_yz = buffer_1000_pddd[274];

    auto g_x_0_0_0_y_xy_yy_zz = buffer_1000_pddd[275];

    auto g_x_0_0_0_y_xy_yz_xx = buffer_1000_pddd[276];

    auto g_x_0_0_0_y_xy_yz_xy = buffer_1000_pddd[277];

    auto g_x_0_0_0_y_xy_yz_xz = buffer_1000_pddd[278];

    auto g_x_0_0_0_y_xy_yz_yy = buffer_1000_pddd[279];

    auto g_x_0_0_0_y_xy_yz_yz = buffer_1000_pddd[280];

    auto g_x_0_0_0_y_xy_yz_zz = buffer_1000_pddd[281];

    auto g_x_0_0_0_y_xy_zz_xx = buffer_1000_pddd[282];

    auto g_x_0_0_0_y_xy_zz_xy = buffer_1000_pddd[283];

    auto g_x_0_0_0_y_xy_zz_xz = buffer_1000_pddd[284];

    auto g_x_0_0_0_y_xy_zz_yy = buffer_1000_pddd[285];

    auto g_x_0_0_0_y_xy_zz_yz = buffer_1000_pddd[286];

    auto g_x_0_0_0_y_xy_zz_zz = buffer_1000_pddd[287];

    auto g_x_0_0_0_y_xz_xx_xx = buffer_1000_pddd[288];

    auto g_x_0_0_0_y_xz_xx_xy = buffer_1000_pddd[289];

    auto g_x_0_0_0_y_xz_xx_xz = buffer_1000_pddd[290];

    auto g_x_0_0_0_y_xz_xx_yy = buffer_1000_pddd[291];

    auto g_x_0_0_0_y_xz_xx_yz = buffer_1000_pddd[292];

    auto g_x_0_0_0_y_xz_xx_zz = buffer_1000_pddd[293];

    auto g_x_0_0_0_y_xz_xy_xx = buffer_1000_pddd[294];

    auto g_x_0_0_0_y_xz_xy_xy = buffer_1000_pddd[295];

    auto g_x_0_0_0_y_xz_xy_xz = buffer_1000_pddd[296];

    auto g_x_0_0_0_y_xz_xy_yy = buffer_1000_pddd[297];

    auto g_x_0_0_0_y_xz_xy_yz = buffer_1000_pddd[298];

    auto g_x_0_0_0_y_xz_xy_zz = buffer_1000_pddd[299];

    auto g_x_0_0_0_y_xz_xz_xx = buffer_1000_pddd[300];

    auto g_x_0_0_0_y_xz_xz_xy = buffer_1000_pddd[301];

    auto g_x_0_0_0_y_xz_xz_xz = buffer_1000_pddd[302];

    auto g_x_0_0_0_y_xz_xz_yy = buffer_1000_pddd[303];

    auto g_x_0_0_0_y_xz_xz_yz = buffer_1000_pddd[304];

    auto g_x_0_0_0_y_xz_xz_zz = buffer_1000_pddd[305];

    auto g_x_0_0_0_y_xz_yy_xx = buffer_1000_pddd[306];

    auto g_x_0_0_0_y_xz_yy_xy = buffer_1000_pddd[307];

    auto g_x_0_0_0_y_xz_yy_xz = buffer_1000_pddd[308];

    auto g_x_0_0_0_y_xz_yy_yy = buffer_1000_pddd[309];

    auto g_x_0_0_0_y_xz_yy_yz = buffer_1000_pddd[310];

    auto g_x_0_0_0_y_xz_yy_zz = buffer_1000_pddd[311];

    auto g_x_0_0_0_y_xz_yz_xx = buffer_1000_pddd[312];

    auto g_x_0_0_0_y_xz_yz_xy = buffer_1000_pddd[313];

    auto g_x_0_0_0_y_xz_yz_xz = buffer_1000_pddd[314];

    auto g_x_0_0_0_y_xz_yz_yy = buffer_1000_pddd[315];

    auto g_x_0_0_0_y_xz_yz_yz = buffer_1000_pddd[316];

    auto g_x_0_0_0_y_xz_yz_zz = buffer_1000_pddd[317];

    auto g_x_0_0_0_y_xz_zz_xx = buffer_1000_pddd[318];

    auto g_x_0_0_0_y_xz_zz_xy = buffer_1000_pddd[319];

    auto g_x_0_0_0_y_xz_zz_xz = buffer_1000_pddd[320];

    auto g_x_0_0_0_y_xz_zz_yy = buffer_1000_pddd[321];

    auto g_x_0_0_0_y_xz_zz_yz = buffer_1000_pddd[322];

    auto g_x_0_0_0_y_xz_zz_zz = buffer_1000_pddd[323];

    auto g_x_0_0_0_y_yy_xx_xx = buffer_1000_pddd[324];

    auto g_x_0_0_0_y_yy_xx_xy = buffer_1000_pddd[325];

    auto g_x_0_0_0_y_yy_xx_xz = buffer_1000_pddd[326];

    auto g_x_0_0_0_y_yy_xx_yy = buffer_1000_pddd[327];

    auto g_x_0_0_0_y_yy_xx_yz = buffer_1000_pddd[328];

    auto g_x_0_0_0_y_yy_xx_zz = buffer_1000_pddd[329];

    auto g_x_0_0_0_y_yy_xy_xx = buffer_1000_pddd[330];

    auto g_x_0_0_0_y_yy_xy_xy = buffer_1000_pddd[331];

    auto g_x_0_0_0_y_yy_xy_xz = buffer_1000_pddd[332];

    auto g_x_0_0_0_y_yy_xy_yy = buffer_1000_pddd[333];

    auto g_x_0_0_0_y_yy_xy_yz = buffer_1000_pddd[334];

    auto g_x_0_0_0_y_yy_xy_zz = buffer_1000_pddd[335];

    auto g_x_0_0_0_y_yy_xz_xx = buffer_1000_pddd[336];

    auto g_x_0_0_0_y_yy_xz_xy = buffer_1000_pddd[337];

    auto g_x_0_0_0_y_yy_xz_xz = buffer_1000_pddd[338];

    auto g_x_0_0_0_y_yy_xz_yy = buffer_1000_pddd[339];

    auto g_x_0_0_0_y_yy_xz_yz = buffer_1000_pddd[340];

    auto g_x_0_0_0_y_yy_xz_zz = buffer_1000_pddd[341];

    auto g_x_0_0_0_y_yy_yy_xx = buffer_1000_pddd[342];

    auto g_x_0_0_0_y_yy_yy_xy = buffer_1000_pddd[343];

    auto g_x_0_0_0_y_yy_yy_xz = buffer_1000_pddd[344];

    auto g_x_0_0_0_y_yy_yy_yy = buffer_1000_pddd[345];

    auto g_x_0_0_0_y_yy_yy_yz = buffer_1000_pddd[346];

    auto g_x_0_0_0_y_yy_yy_zz = buffer_1000_pddd[347];

    auto g_x_0_0_0_y_yy_yz_xx = buffer_1000_pddd[348];

    auto g_x_0_0_0_y_yy_yz_xy = buffer_1000_pddd[349];

    auto g_x_0_0_0_y_yy_yz_xz = buffer_1000_pddd[350];

    auto g_x_0_0_0_y_yy_yz_yy = buffer_1000_pddd[351];

    auto g_x_0_0_0_y_yy_yz_yz = buffer_1000_pddd[352];

    auto g_x_0_0_0_y_yy_yz_zz = buffer_1000_pddd[353];

    auto g_x_0_0_0_y_yy_zz_xx = buffer_1000_pddd[354];

    auto g_x_0_0_0_y_yy_zz_xy = buffer_1000_pddd[355];

    auto g_x_0_0_0_y_yy_zz_xz = buffer_1000_pddd[356];

    auto g_x_0_0_0_y_yy_zz_yy = buffer_1000_pddd[357];

    auto g_x_0_0_0_y_yy_zz_yz = buffer_1000_pddd[358];

    auto g_x_0_0_0_y_yy_zz_zz = buffer_1000_pddd[359];

    auto g_x_0_0_0_y_yz_xx_xx = buffer_1000_pddd[360];

    auto g_x_0_0_0_y_yz_xx_xy = buffer_1000_pddd[361];

    auto g_x_0_0_0_y_yz_xx_xz = buffer_1000_pddd[362];

    auto g_x_0_0_0_y_yz_xx_yy = buffer_1000_pddd[363];

    auto g_x_0_0_0_y_yz_xx_yz = buffer_1000_pddd[364];

    auto g_x_0_0_0_y_yz_xx_zz = buffer_1000_pddd[365];

    auto g_x_0_0_0_y_yz_xy_xx = buffer_1000_pddd[366];

    auto g_x_0_0_0_y_yz_xy_xy = buffer_1000_pddd[367];

    auto g_x_0_0_0_y_yz_xy_xz = buffer_1000_pddd[368];

    auto g_x_0_0_0_y_yz_xy_yy = buffer_1000_pddd[369];

    auto g_x_0_0_0_y_yz_xy_yz = buffer_1000_pddd[370];

    auto g_x_0_0_0_y_yz_xy_zz = buffer_1000_pddd[371];

    auto g_x_0_0_0_y_yz_xz_xx = buffer_1000_pddd[372];

    auto g_x_0_0_0_y_yz_xz_xy = buffer_1000_pddd[373];

    auto g_x_0_0_0_y_yz_xz_xz = buffer_1000_pddd[374];

    auto g_x_0_0_0_y_yz_xz_yy = buffer_1000_pddd[375];

    auto g_x_0_0_0_y_yz_xz_yz = buffer_1000_pddd[376];

    auto g_x_0_0_0_y_yz_xz_zz = buffer_1000_pddd[377];

    auto g_x_0_0_0_y_yz_yy_xx = buffer_1000_pddd[378];

    auto g_x_0_0_0_y_yz_yy_xy = buffer_1000_pddd[379];

    auto g_x_0_0_0_y_yz_yy_xz = buffer_1000_pddd[380];

    auto g_x_0_0_0_y_yz_yy_yy = buffer_1000_pddd[381];

    auto g_x_0_0_0_y_yz_yy_yz = buffer_1000_pddd[382];

    auto g_x_0_0_0_y_yz_yy_zz = buffer_1000_pddd[383];

    auto g_x_0_0_0_y_yz_yz_xx = buffer_1000_pddd[384];

    auto g_x_0_0_0_y_yz_yz_xy = buffer_1000_pddd[385];

    auto g_x_0_0_0_y_yz_yz_xz = buffer_1000_pddd[386];

    auto g_x_0_0_0_y_yz_yz_yy = buffer_1000_pddd[387];

    auto g_x_0_0_0_y_yz_yz_yz = buffer_1000_pddd[388];

    auto g_x_0_0_0_y_yz_yz_zz = buffer_1000_pddd[389];

    auto g_x_0_0_0_y_yz_zz_xx = buffer_1000_pddd[390];

    auto g_x_0_0_0_y_yz_zz_xy = buffer_1000_pddd[391];

    auto g_x_0_0_0_y_yz_zz_xz = buffer_1000_pddd[392];

    auto g_x_0_0_0_y_yz_zz_yy = buffer_1000_pddd[393];

    auto g_x_0_0_0_y_yz_zz_yz = buffer_1000_pddd[394];

    auto g_x_0_0_0_y_yz_zz_zz = buffer_1000_pddd[395];

    auto g_x_0_0_0_y_zz_xx_xx = buffer_1000_pddd[396];

    auto g_x_0_0_0_y_zz_xx_xy = buffer_1000_pddd[397];

    auto g_x_0_0_0_y_zz_xx_xz = buffer_1000_pddd[398];

    auto g_x_0_0_0_y_zz_xx_yy = buffer_1000_pddd[399];

    auto g_x_0_0_0_y_zz_xx_yz = buffer_1000_pddd[400];

    auto g_x_0_0_0_y_zz_xx_zz = buffer_1000_pddd[401];

    auto g_x_0_0_0_y_zz_xy_xx = buffer_1000_pddd[402];

    auto g_x_0_0_0_y_zz_xy_xy = buffer_1000_pddd[403];

    auto g_x_0_0_0_y_zz_xy_xz = buffer_1000_pddd[404];

    auto g_x_0_0_0_y_zz_xy_yy = buffer_1000_pddd[405];

    auto g_x_0_0_0_y_zz_xy_yz = buffer_1000_pddd[406];

    auto g_x_0_0_0_y_zz_xy_zz = buffer_1000_pddd[407];

    auto g_x_0_0_0_y_zz_xz_xx = buffer_1000_pddd[408];

    auto g_x_0_0_0_y_zz_xz_xy = buffer_1000_pddd[409];

    auto g_x_0_0_0_y_zz_xz_xz = buffer_1000_pddd[410];

    auto g_x_0_0_0_y_zz_xz_yy = buffer_1000_pddd[411];

    auto g_x_0_0_0_y_zz_xz_yz = buffer_1000_pddd[412];

    auto g_x_0_0_0_y_zz_xz_zz = buffer_1000_pddd[413];

    auto g_x_0_0_0_y_zz_yy_xx = buffer_1000_pddd[414];

    auto g_x_0_0_0_y_zz_yy_xy = buffer_1000_pddd[415];

    auto g_x_0_0_0_y_zz_yy_xz = buffer_1000_pddd[416];

    auto g_x_0_0_0_y_zz_yy_yy = buffer_1000_pddd[417];

    auto g_x_0_0_0_y_zz_yy_yz = buffer_1000_pddd[418];

    auto g_x_0_0_0_y_zz_yy_zz = buffer_1000_pddd[419];

    auto g_x_0_0_0_y_zz_yz_xx = buffer_1000_pddd[420];

    auto g_x_0_0_0_y_zz_yz_xy = buffer_1000_pddd[421];

    auto g_x_0_0_0_y_zz_yz_xz = buffer_1000_pddd[422];

    auto g_x_0_0_0_y_zz_yz_yy = buffer_1000_pddd[423];

    auto g_x_0_0_0_y_zz_yz_yz = buffer_1000_pddd[424];

    auto g_x_0_0_0_y_zz_yz_zz = buffer_1000_pddd[425];

    auto g_x_0_0_0_y_zz_zz_xx = buffer_1000_pddd[426];

    auto g_x_0_0_0_y_zz_zz_xy = buffer_1000_pddd[427];

    auto g_x_0_0_0_y_zz_zz_xz = buffer_1000_pddd[428];

    auto g_x_0_0_0_y_zz_zz_yy = buffer_1000_pddd[429];

    auto g_x_0_0_0_y_zz_zz_yz = buffer_1000_pddd[430];

    auto g_x_0_0_0_y_zz_zz_zz = buffer_1000_pddd[431];

    auto g_x_0_0_0_z_xx_xx_xx = buffer_1000_pddd[432];

    auto g_x_0_0_0_z_xx_xx_xy = buffer_1000_pddd[433];

    auto g_x_0_0_0_z_xx_xx_xz = buffer_1000_pddd[434];

    auto g_x_0_0_0_z_xx_xx_yy = buffer_1000_pddd[435];

    auto g_x_0_0_0_z_xx_xx_yz = buffer_1000_pddd[436];

    auto g_x_0_0_0_z_xx_xx_zz = buffer_1000_pddd[437];

    auto g_x_0_0_0_z_xx_xy_xx = buffer_1000_pddd[438];

    auto g_x_0_0_0_z_xx_xy_xy = buffer_1000_pddd[439];

    auto g_x_0_0_0_z_xx_xy_xz = buffer_1000_pddd[440];

    auto g_x_0_0_0_z_xx_xy_yy = buffer_1000_pddd[441];

    auto g_x_0_0_0_z_xx_xy_yz = buffer_1000_pddd[442];

    auto g_x_0_0_0_z_xx_xy_zz = buffer_1000_pddd[443];

    auto g_x_0_0_0_z_xx_xz_xx = buffer_1000_pddd[444];

    auto g_x_0_0_0_z_xx_xz_xy = buffer_1000_pddd[445];

    auto g_x_0_0_0_z_xx_xz_xz = buffer_1000_pddd[446];

    auto g_x_0_0_0_z_xx_xz_yy = buffer_1000_pddd[447];

    auto g_x_0_0_0_z_xx_xz_yz = buffer_1000_pddd[448];

    auto g_x_0_0_0_z_xx_xz_zz = buffer_1000_pddd[449];

    auto g_x_0_0_0_z_xx_yy_xx = buffer_1000_pddd[450];

    auto g_x_0_0_0_z_xx_yy_xy = buffer_1000_pddd[451];

    auto g_x_0_0_0_z_xx_yy_xz = buffer_1000_pddd[452];

    auto g_x_0_0_0_z_xx_yy_yy = buffer_1000_pddd[453];

    auto g_x_0_0_0_z_xx_yy_yz = buffer_1000_pddd[454];

    auto g_x_0_0_0_z_xx_yy_zz = buffer_1000_pddd[455];

    auto g_x_0_0_0_z_xx_yz_xx = buffer_1000_pddd[456];

    auto g_x_0_0_0_z_xx_yz_xy = buffer_1000_pddd[457];

    auto g_x_0_0_0_z_xx_yz_xz = buffer_1000_pddd[458];

    auto g_x_0_0_0_z_xx_yz_yy = buffer_1000_pddd[459];

    auto g_x_0_0_0_z_xx_yz_yz = buffer_1000_pddd[460];

    auto g_x_0_0_0_z_xx_yz_zz = buffer_1000_pddd[461];

    auto g_x_0_0_0_z_xx_zz_xx = buffer_1000_pddd[462];

    auto g_x_0_0_0_z_xx_zz_xy = buffer_1000_pddd[463];

    auto g_x_0_0_0_z_xx_zz_xz = buffer_1000_pddd[464];

    auto g_x_0_0_0_z_xx_zz_yy = buffer_1000_pddd[465];

    auto g_x_0_0_0_z_xx_zz_yz = buffer_1000_pddd[466];

    auto g_x_0_0_0_z_xx_zz_zz = buffer_1000_pddd[467];

    auto g_x_0_0_0_z_xy_xx_xx = buffer_1000_pddd[468];

    auto g_x_0_0_0_z_xy_xx_xy = buffer_1000_pddd[469];

    auto g_x_0_0_0_z_xy_xx_xz = buffer_1000_pddd[470];

    auto g_x_0_0_0_z_xy_xx_yy = buffer_1000_pddd[471];

    auto g_x_0_0_0_z_xy_xx_yz = buffer_1000_pddd[472];

    auto g_x_0_0_0_z_xy_xx_zz = buffer_1000_pddd[473];

    auto g_x_0_0_0_z_xy_xy_xx = buffer_1000_pddd[474];

    auto g_x_0_0_0_z_xy_xy_xy = buffer_1000_pddd[475];

    auto g_x_0_0_0_z_xy_xy_xz = buffer_1000_pddd[476];

    auto g_x_0_0_0_z_xy_xy_yy = buffer_1000_pddd[477];

    auto g_x_0_0_0_z_xy_xy_yz = buffer_1000_pddd[478];

    auto g_x_0_0_0_z_xy_xy_zz = buffer_1000_pddd[479];

    auto g_x_0_0_0_z_xy_xz_xx = buffer_1000_pddd[480];

    auto g_x_0_0_0_z_xy_xz_xy = buffer_1000_pddd[481];

    auto g_x_0_0_0_z_xy_xz_xz = buffer_1000_pddd[482];

    auto g_x_0_0_0_z_xy_xz_yy = buffer_1000_pddd[483];

    auto g_x_0_0_0_z_xy_xz_yz = buffer_1000_pddd[484];

    auto g_x_0_0_0_z_xy_xz_zz = buffer_1000_pddd[485];

    auto g_x_0_0_0_z_xy_yy_xx = buffer_1000_pddd[486];

    auto g_x_0_0_0_z_xy_yy_xy = buffer_1000_pddd[487];

    auto g_x_0_0_0_z_xy_yy_xz = buffer_1000_pddd[488];

    auto g_x_0_0_0_z_xy_yy_yy = buffer_1000_pddd[489];

    auto g_x_0_0_0_z_xy_yy_yz = buffer_1000_pddd[490];

    auto g_x_0_0_0_z_xy_yy_zz = buffer_1000_pddd[491];

    auto g_x_0_0_0_z_xy_yz_xx = buffer_1000_pddd[492];

    auto g_x_0_0_0_z_xy_yz_xy = buffer_1000_pddd[493];

    auto g_x_0_0_0_z_xy_yz_xz = buffer_1000_pddd[494];

    auto g_x_0_0_0_z_xy_yz_yy = buffer_1000_pddd[495];

    auto g_x_0_0_0_z_xy_yz_yz = buffer_1000_pddd[496];

    auto g_x_0_0_0_z_xy_yz_zz = buffer_1000_pddd[497];

    auto g_x_0_0_0_z_xy_zz_xx = buffer_1000_pddd[498];

    auto g_x_0_0_0_z_xy_zz_xy = buffer_1000_pddd[499];

    auto g_x_0_0_0_z_xy_zz_xz = buffer_1000_pddd[500];

    auto g_x_0_0_0_z_xy_zz_yy = buffer_1000_pddd[501];

    auto g_x_0_0_0_z_xy_zz_yz = buffer_1000_pddd[502];

    auto g_x_0_0_0_z_xy_zz_zz = buffer_1000_pddd[503];

    auto g_x_0_0_0_z_xz_xx_xx = buffer_1000_pddd[504];

    auto g_x_0_0_0_z_xz_xx_xy = buffer_1000_pddd[505];

    auto g_x_0_0_0_z_xz_xx_xz = buffer_1000_pddd[506];

    auto g_x_0_0_0_z_xz_xx_yy = buffer_1000_pddd[507];

    auto g_x_0_0_0_z_xz_xx_yz = buffer_1000_pddd[508];

    auto g_x_0_0_0_z_xz_xx_zz = buffer_1000_pddd[509];

    auto g_x_0_0_0_z_xz_xy_xx = buffer_1000_pddd[510];

    auto g_x_0_0_0_z_xz_xy_xy = buffer_1000_pddd[511];

    auto g_x_0_0_0_z_xz_xy_xz = buffer_1000_pddd[512];

    auto g_x_0_0_0_z_xz_xy_yy = buffer_1000_pddd[513];

    auto g_x_0_0_0_z_xz_xy_yz = buffer_1000_pddd[514];

    auto g_x_0_0_0_z_xz_xy_zz = buffer_1000_pddd[515];

    auto g_x_0_0_0_z_xz_xz_xx = buffer_1000_pddd[516];

    auto g_x_0_0_0_z_xz_xz_xy = buffer_1000_pddd[517];

    auto g_x_0_0_0_z_xz_xz_xz = buffer_1000_pddd[518];

    auto g_x_0_0_0_z_xz_xz_yy = buffer_1000_pddd[519];

    auto g_x_0_0_0_z_xz_xz_yz = buffer_1000_pddd[520];

    auto g_x_0_0_0_z_xz_xz_zz = buffer_1000_pddd[521];

    auto g_x_0_0_0_z_xz_yy_xx = buffer_1000_pddd[522];

    auto g_x_0_0_0_z_xz_yy_xy = buffer_1000_pddd[523];

    auto g_x_0_0_0_z_xz_yy_xz = buffer_1000_pddd[524];

    auto g_x_0_0_0_z_xz_yy_yy = buffer_1000_pddd[525];

    auto g_x_0_0_0_z_xz_yy_yz = buffer_1000_pddd[526];

    auto g_x_0_0_0_z_xz_yy_zz = buffer_1000_pddd[527];

    auto g_x_0_0_0_z_xz_yz_xx = buffer_1000_pddd[528];

    auto g_x_0_0_0_z_xz_yz_xy = buffer_1000_pddd[529];

    auto g_x_0_0_0_z_xz_yz_xz = buffer_1000_pddd[530];

    auto g_x_0_0_0_z_xz_yz_yy = buffer_1000_pddd[531];

    auto g_x_0_0_0_z_xz_yz_yz = buffer_1000_pddd[532];

    auto g_x_0_0_0_z_xz_yz_zz = buffer_1000_pddd[533];

    auto g_x_0_0_0_z_xz_zz_xx = buffer_1000_pddd[534];

    auto g_x_0_0_0_z_xz_zz_xy = buffer_1000_pddd[535];

    auto g_x_0_0_0_z_xz_zz_xz = buffer_1000_pddd[536];

    auto g_x_0_0_0_z_xz_zz_yy = buffer_1000_pddd[537];

    auto g_x_0_0_0_z_xz_zz_yz = buffer_1000_pddd[538];

    auto g_x_0_0_0_z_xz_zz_zz = buffer_1000_pddd[539];

    auto g_x_0_0_0_z_yy_xx_xx = buffer_1000_pddd[540];

    auto g_x_0_0_0_z_yy_xx_xy = buffer_1000_pddd[541];

    auto g_x_0_0_0_z_yy_xx_xz = buffer_1000_pddd[542];

    auto g_x_0_0_0_z_yy_xx_yy = buffer_1000_pddd[543];

    auto g_x_0_0_0_z_yy_xx_yz = buffer_1000_pddd[544];

    auto g_x_0_0_0_z_yy_xx_zz = buffer_1000_pddd[545];

    auto g_x_0_0_0_z_yy_xy_xx = buffer_1000_pddd[546];

    auto g_x_0_0_0_z_yy_xy_xy = buffer_1000_pddd[547];

    auto g_x_0_0_0_z_yy_xy_xz = buffer_1000_pddd[548];

    auto g_x_0_0_0_z_yy_xy_yy = buffer_1000_pddd[549];

    auto g_x_0_0_0_z_yy_xy_yz = buffer_1000_pddd[550];

    auto g_x_0_0_0_z_yy_xy_zz = buffer_1000_pddd[551];

    auto g_x_0_0_0_z_yy_xz_xx = buffer_1000_pddd[552];

    auto g_x_0_0_0_z_yy_xz_xy = buffer_1000_pddd[553];

    auto g_x_0_0_0_z_yy_xz_xz = buffer_1000_pddd[554];

    auto g_x_0_0_0_z_yy_xz_yy = buffer_1000_pddd[555];

    auto g_x_0_0_0_z_yy_xz_yz = buffer_1000_pddd[556];

    auto g_x_0_0_0_z_yy_xz_zz = buffer_1000_pddd[557];

    auto g_x_0_0_0_z_yy_yy_xx = buffer_1000_pddd[558];

    auto g_x_0_0_0_z_yy_yy_xy = buffer_1000_pddd[559];

    auto g_x_0_0_0_z_yy_yy_xz = buffer_1000_pddd[560];

    auto g_x_0_0_0_z_yy_yy_yy = buffer_1000_pddd[561];

    auto g_x_0_0_0_z_yy_yy_yz = buffer_1000_pddd[562];

    auto g_x_0_0_0_z_yy_yy_zz = buffer_1000_pddd[563];

    auto g_x_0_0_0_z_yy_yz_xx = buffer_1000_pddd[564];

    auto g_x_0_0_0_z_yy_yz_xy = buffer_1000_pddd[565];

    auto g_x_0_0_0_z_yy_yz_xz = buffer_1000_pddd[566];

    auto g_x_0_0_0_z_yy_yz_yy = buffer_1000_pddd[567];

    auto g_x_0_0_0_z_yy_yz_yz = buffer_1000_pddd[568];

    auto g_x_0_0_0_z_yy_yz_zz = buffer_1000_pddd[569];

    auto g_x_0_0_0_z_yy_zz_xx = buffer_1000_pddd[570];

    auto g_x_0_0_0_z_yy_zz_xy = buffer_1000_pddd[571];

    auto g_x_0_0_0_z_yy_zz_xz = buffer_1000_pddd[572];

    auto g_x_0_0_0_z_yy_zz_yy = buffer_1000_pddd[573];

    auto g_x_0_0_0_z_yy_zz_yz = buffer_1000_pddd[574];

    auto g_x_0_0_0_z_yy_zz_zz = buffer_1000_pddd[575];

    auto g_x_0_0_0_z_yz_xx_xx = buffer_1000_pddd[576];

    auto g_x_0_0_0_z_yz_xx_xy = buffer_1000_pddd[577];

    auto g_x_0_0_0_z_yz_xx_xz = buffer_1000_pddd[578];

    auto g_x_0_0_0_z_yz_xx_yy = buffer_1000_pddd[579];

    auto g_x_0_0_0_z_yz_xx_yz = buffer_1000_pddd[580];

    auto g_x_0_0_0_z_yz_xx_zz = buffer_1000_pddd[581];

    auto g_x_0_0_0_z_yz_xy_xx = buffer_1000_pddd[582];

    auto g_x_0_0_0_z_yz_xy_xy = buffer_1000_pddd[583];

    auto g_x_0_0_0_z_yz_xy_xz = buffer_1000_pddd[584];

    auto g_x_0_0_0_z_yz_xy_yy = buffer_1000_pddd[585];

    auto g_x_0_0_0_z_yz_xy_yz = buffer_1000_pddd[586];

    auto g_x_0_0_0_z_yz_xy_zz = buffer_1000_pddd[587];

    auto g_x_0_0_0_z_yz_xz_xx = buffer_1000_pddd[588];

    auto g_x_0_0_0_z_yz_xz_xy = buffer_1000_pddd[589];

    auto g_x_0_0_0_z_yz_xz_xz = buffer_1000_pddd[590];

    auto g_x_0_0_0_z_yz_xz_yy = buffer_1000_pddd[591];

    auto g_x_0_0_0_z_yz_xz_yz = buffer_1000_pddd[592];

    auto g_x_0_0_0_z_yz_xz_zz = buffer_1000_pddd[593];

    auto g_x_0_0_0_z_yz_yy_xx = buffer_1000_pddd[594];

    auto g_x_0_0_0_z_yz_yy_xy = buffer_1000_pddd[595];

    auto g_x_0_0_0_z_yz_yy_xz = buffer_1000_pddd[596];

    auto g_x_0_0_0_z_yz_yy_yy = buffer_1000_pddd[597];

    auto g_x_0_0_0_z_yz_yy_yz = buffer_1000_pddd[598];

    auto g_x_0_0_0_z_yz_yy_zz = buffer_1000_pddd[599];

    auto g_x_0_0_0_z_yz_yz_xx = buffer_1000_pddd[600];

    auto g_x_0_0_0_z_yz_yz_xy = buffer_1000_pddd[601];

    auto g_x_0_0_0_z_yz_yz_xz = buffer_1000_pddd[602];

    auto g_x_0_0_0_z_yz_yz_yy = buffer_1000_pddd[603];

    auto g_x_0_0_0_z_yz_yz_yz = buffer_1000_pddd[604];

    auto g_x_0_0_0_z_yz_yz_zz = buffer_1000_pddd[605];

    auto g_x_0_0_0_z_yz_zz_xx = buffer_1000_pddd[606];

    auto g_x_0_0_0_z_yz_zz_xy = buffer_1000_pddd[607];

    auto g_x_0_0_0_z_yz_zz_xz = buffer_1000_pddd[608];

    auto g_x_0_0_0_z_yz_zz_yy = buffer_1000_pddd[609];

    auto g_x_0_0_0_z_yz_zz_yz = buffer_1000_pddd[610];

    auto g_x_0_0_0_z_yz_zz_zz = buffer_1000_pddd[611];

    auto g_x_0_0_0_z_zz_xx_xx = buffer_1000_pddd[612];

    auto g_x_0_0_0_z_zz_xx_xy = buffer_1000_pddd[613];

    auto g_x_0_0_0_z_zz_xx_xz = buffer_1000_pddd[614];

    auto g_x_0_0_0_z_zz_xx_yy = buffer_1000_pddd[615];

    auto g_x_0_0_0_z_zz_xx_yz = buffer_1000_pddd[616];

    auto g_x_0_0_0_z_zz_xx_zz = buffer_1000_pddd[617];

    auto g_x_0_0_0_z_zz_xy_xx = buffer_1000_pddd[618];

    auto g_x_0_0_0_z_zz_xy_xy = buffer_1000_pddd[619];

    auto g_x_0_0_0_z_zz_xy_xz = buffer_1000_pddd[620];

    auto g_x_0_0_0_z_zz_xy_yy = buffer_1000_pddd[621];

    auto g_x_0_0_0_z_zz_xy_yz = buffer_1000_pddd[622];

    auto g_x_0_0_0_z_zz_xy_zz = buffer_1000_pddd[623];

    auto g_x_0_0_0_z_zz_xz_xx = buffer_1000_pddd[624];

    auto g_x_0_0_0_z_zz_xz_xy = buffer_1000_pddd[625];

    auto g_x_0_0_0_z_zz_xz_xz = buffer_1000_pddd[626];

    auto g_x_0_0_0_z_zz_xz_yy = buffer_1000_pddd[627];

    auto g_x_0_0_0_z_zz_xz_yz = buffer_1000_pddd[628];

    auto g_x_0_0_0_z_zz_xz_zz = buffer_1000_pddd[629];

    auto g_x_0_0_0_z_zz_yy_xx = buffer_1000_pddd[630];

    auto g_x_0_0_0_z_zz_yy_xy = buffer_1000_pddd[631];

    auto g_x_0_0_0_z_zz_yy_xz = buffer_1000_pddd[632];

    auto g_x_0_0_0_z_zz_yy_yy = buffer_1000_pddd[633];

    auto g_x_0_0_0_z_zz_yy_yz = buffer_1000_pddd[634];

    auto g_x_0_0_0_z_zz_yy_zz = buffer_1000_pddd[635];

    auto g_x_0_0_0_z_zz_yz_xx = buffer_1000_pddd[636];

    auto g_x_0_0_0_z_zz_yz_xy = buffer_1000_pddd[637];

    auto g_x_0_0_0_z_zz_yz_xz = buffer_1000_pddd[638];

    auto g_x_0_0_0_z_zz_yz_yy = buffer_1000_pddd[639];

    auto g_x_0_0_0_z_zz_yz_yz = buffer_1000_pddd[640];

    auto g_x_0_0_0_z_zz_yz_zz = buffer_1000_pddd[641];

    auto g_x_0_0_0_z_zz_zz_xx = buffer_1000_pddd[642];

    auto g_x_0_0_0_z_zz_zz_xy = buffer_1000_pddd[643];

    auto g_x_0_0_0_z_zz_zz_xz = buffer_1000_pddd[644];

    auto g_x_0_0_0_z_zz_zz_yy = buffer_1000_pddd[645];

    auto g_x_0_0_0_z_zz_zz_yz = buffer_1000_pddd[646];

    auto g_x_0_0_0_z_zz_zz_zz = buffer_1000_pddd[647];

    auto g_y_0_0_0_x_xx_xx_xx = buffer_1000_pddd[648];

    auto g_y_0_0_0_x_xx_xx_xy = buffer_1000_pddd[649];

    auto g_y_0_0_0_x_xx_xx_xz = buffer_1000_pddd[650];

    auto g_y_0_0_0_x_xx_xx_yy = buffer_1000_pddd[651];

    auto g_y_0_0_0_x_xx_xx_yz = buffer_1000_pddd[652];

    auto g_y_0_0_0_x_xx_xx_zz = buffer_1000_pddd[653];

    auto g_y_0_0_0_x_xx_xy_xx = buffer_1000_pddd[654];

    auto g_y_0_0_0_x_xx_xy_xy = buffer_1000_pddd[655];

    auto g_y_0_0_0_x_xx_xy_xz = buffer_1000_pddd[656];

    auto g_y_0_0_0_x_xx_xy_yy = buffer_1000_pddd[657];

    auto g_y_0_0_0_x_xx_xy_yz = buffer_1000_pddd[658];

    auto g_y_0_0_0_x_xx_xy_zz = buffer_1000_pddd[659];

    auto g_y_0_0_0_x_xx_xz_xx = buffer_1000_pddd[660];

    auto g_y_0_0_0_x_xx_xz_xy = buffer_1000_pddd[661];

    auto g_y_0_0_0_x_xx_xz_xz = buffer_1000_pddd[662];

    auto g_y_0_0_0_x_xx_xz_yy = buffer_1000_pddd[663];

    auto g_y_0_0_0_x_xx_xz_yz = buffer_1000_pddd[664];

    auto g_y_0_0_0_x_xx_xz_zz = buffer_1000_pddd[665];

    auto g_y_0_0_0_x_xx_yy_xx = buffer_1000_pddd[666];

    auto g_y_0_0_0_x_xx_yy_xy = buffer_1000_pddd[667];

    auto g_y_0_0_0_x_xx_yy_xz = buffer_1000_pddd[668];

    auto g_y_0_0_0_x_xx_yy_yy = buffer_1000_pddd[669];

    auto g_y_0_0_0_x_xx_yy_yz = buffer_1000_pddd[670];

    auto g_y_0_0_0_x_xx_yy_zz = buffer_1000_pddd[671];

    auto g_y_0_0_0_x_xx_yz_xx = buffer_1000_pddd[672];

    auto g_y_0_0_0_x_xx_yz_xy = buffer_1000_pddd[673];

    auto g_y_0_0_0_x_xx_yz_xz = buffer_1000_pddd[674];

    auto g_y_0_0_0_x_xx_yz_yy = buffer_1000_pddd[675];

    auto g_y_0_0_0_x_xx_yz_yz = buffer_1000_pddd[676];

    auto g_y_0_0_0_x_xx_yz_zz = buffer_1000_pddd[677];

    auto g_y_0_0_0_x_xx_zz_xx = buffer_1000_pddd[678];

    auto g_y_0_0_0_x_xx_zz_xy = buffer_1000_pddd[679];

    auto g_y_0_0_0_x_xx_zz_xz = buffer_1000_pddd[680];

    auto g_y_0_0_0_x_xx_zz_yy = buffer_1000_pddd[681];

    auto g_y_0_0_0_x_xx_zz_yz = buffer_1000_pddd[682];

    auto g_y_0_0_0_x_xx_zz_zz = buffer_1000_pddd[683];

    auto g_y_0_0_0_x_xy_xx_xx = buffer_1000_pddd[684];

    auto g_y_0_0_0_x_xy_xx_xy = buffer_1000_pddd[685];

    auto g_y_0_0_0_x_xy_xx_xz = buffer_1000_pddd[686];

    auto g_y_0_0_0_x_xy_xx_yy = buffer_1000_pddd[687];

    auto g_y_0_0_0_x_xy_xx_yz = buffer_1000_pddd[688];

    auto g_y_0_0_0_x_xy_xx_zz = buffer_1000_pddd[689];

    auto g_y_0_0_0_x_xy_xy_xx = buffer_1000_pddd[690];

    auto g_y_0_0_0_x_xy_xy_xy = buffer_1000_pddd[691];

    auto g_y_0_0_0_x_xy_xy_xz = buffer_1000_pddd[692];

    auto g_y_0_0_0_x_xy_xy_yy = buffer_1000_pddd[693];

    auto g_y_0_0_0_x_xy_xy_yz = buffer_1000_pddd[694];

    auto g_y_0_0_0_x_xy_xy_zz = buffer_1000_pddd[695];

    auto g_y_0_0_0_x_xy_xz_xx = buffer_1000_pddd[696];

    auto g_y_0_0_0_x_xy_xz_xy = buffer_1000_pddd[697];

    auto g_y_0_0_0_x_xy_xz_xz = buffer_1000_pddd[698];

    auto g_y_0_0_0_x_xy_xz_yy = buffer_1000_pddd[699];

    auto g_y_0_0_0_x_xy_xz_yz = buffer_1000_pddd[700];

    auto g_y_0_0_0_x_xy_xz_zz = buffer_1000_pddd[701];

    auto g_y_0_0_0_x_xy_yy_xx = buffer_1000_pddd[702];

    auto g_y_0_0_0_x_xy_yy_xy = buffer_1000_pddd[703];

    auto g_y_0_0_0_x_xy_yy_xz = buffer_1000_pddd[704];

    auto g_y_0_0_0_x_xy_yy_yy = buffer_1000_pddd[705];

    auto g_y_0_0_0_x_xy_yy_yz = buffer_1000_pddd[706];

    auto g_y_0_0_0_x_xy_yy_zz = buffer_1000_pddd[707];

    auto g_y_0_0_0_x_xy_yz_xx = buffer_1000_pddd[708];

    auto g_y_0_0_0_x_xy_yz_xy = buffer_1000_pddd[709];

    auto g_y_0_0_0_x_xy_yz_xz = buffer_1000_pddd[710];

    auto g_y_0_0_0_x_xy_yz_yy = buffer_1000_pddd[711];

    auto g_y_0_0_0_x_xy_yz_yz = buffer_1000_pddd[712];

    auto g_y_0_0_0_x_xy_yz_zz = buffer_1000_pddd[713];

    auto g_y_0_0_0_x_xy_zz_xx = buffer_1000_pddd[714];

    auto g_y_0_0_0_x_xy_zz_xy = buffer_1000_pddd[715];

    auto g_y_0_0_0_x_xy_zz_xz = buffer_1000_pddd[716];

    auto g_y_0_0_0_x_xy_zz_yy = buffer_1000_pddd[717];

    auto g_y_0_0_0_x_xy_zz_yz = buffer_1000_pddd[718];

    auto g_y_0_0_0_x_xy_zz_zz = buffer_1000_pddd[719];

    auto g_y_0_0_0_x_xz_xx_xx = buffer_1000_pddd[720];

    auto g_y_0_0_0_x_xz_xx_xy = buffer_1000_pddd[721];

    auto g_y_0_0_0_x_xz_xx_xz = buffer_1000_pddd[722];

    auto g_y_0_0_0_x_xz_xx_yy = buffer_1000_pddd[723];

    auto g_y_0_0_0_x_xz_xx_yz = buffer_1000_pddd[724];

    auto g_y_0_0_0_x_xz_xx_zz = buffer_1000_pddd[725];

    auto g_y_0_0_0_x_xz_xy_xx = buffer_1000_pddd[726];

    auto g_y_0_0_0_x_xz_xy_xy = buffer_1000_pddd[727];

    auto g_y_0_0_0_x_xz_xy_xz = buffer_1000_pddd[728];

    auto g_y_0_0_0_x_xz_xy_yy = buffer_1000_pddd[729];

    auto g_y_0_0_0_x_xz_xy_yz = buffer_1000_pddd[730];

    auto g_y_0_0_0_x_xz_xy_zz = buffer_1000_pddd[731];

    auto g_y_0_0_0_x_xz_xz_xx = buffer_1000_pddd[732];

    auto g_y_0_0_0_x_xz_xz_xy = buffer_1000_pddd[733];

    auto g_y_0_0_0_x_xz_xz_xz = buffer_1000_pddd[734];

    auto g_y_0_0_0_x_xz_xz_yy = buffer_1000_pddd[735];

    auto g_y_0_0_0_x_xz_xz_yz = buffer_1000_pddd[736];

    auto g_y_0_0_0_x_xz_xz_zz = buffer_1000_pddd[737];

    auto g_y_0_0_0_x_xz_yy_xx = buffer_1000_pddd[738];

    auto g_y_0_0_0_x_xz_yy_xy = buffer_1000_pddd[739];

    auto g_y_0_0_0_x_xz_yy_xz = buffer_1000_pddd[740];

    auto g_y_0_0_0_x_xz_yy_yy = buffer_1000_pddd[741];

    auto g_y_0_0_0_x_xz_yy_yz = buffer_1000_pddd[742];

    auto g_y_0_0_0_x_xz_yy_zz = buffer_1000_pddd[743];

    auto g_y_0_0_0_x_xz_yz_xx = buffer_1000_pddd[744];

    auto g_y_0_0_0_x_xz_yz_xy = buffer_1000_pddd[745];

    auto g_y_0_0_0_x_xz_yz_xz = buffer_1000_pddd[746];

    auto g_y_0_0_0_x_xz_yz_yy = buffer_1000_pddd[747];

    auto g_y_0_0_0_x_xz_yz_yz = buffer_1000_pddd[748];

    auto g_y_0_0_0_x_xz_yz_zz = buffer_1000_pddd[749];

    auto g_y_0_0_0_x_xz_zz_xx = buffer_1000_pddd[750];

    auto g_y_0_0_0_x_xz_zz_xy = buffer_1000_pddd[751];

    auto g_y_0_0_0_x_xz_zz_xz = buffer_1000_pddd[752];

    auto g_y_0_0_0_x_xz_zz_yy = buffer_1000_pddd[753];

    auto g_y_0_0_0_x_xz_zz_yz = buffer_1000_pddd[754];

    auto g_y_0_0_0_x_xz_zz_zz = buffer_1000_pddd[755];

    auto g_y_0_0_0_x_yy_xx_xx = buffer_1000_pddd[756];

    auto g_y_0_0_0_x_yy_xx_xy = buffer_1000_pddd[757];

    auto g_y_0_0_0_x_yy_xx_xz = buffer_1000_pddd[758];

    auto g_y_0_0_0_x_yy_xx_yy = buffer_1000_pddd[759];

    auto g_y_0_0_0_x_yy_xx_yz = buffer_1000_pddd[760];

    auto g_y_0_0_0_x_yy_xx_zz = buffer_1000_pddd[761];

    auto g_y_0_0_0_x_yy_xy_xx = buffer_1000_pddd[762];

    auto g_y_0_0_0_x_yy_xy_xy = buffer_1000_pddd[763];

    auto g_y_0_0_0_x_yy_xy_xz = buffer_1000_pddd[764];

    auto g_y_0_0_0_x_yy_xy_yy = buffer_1000_pddd[765];

    auto g_y_0_0_0_x_yy_xy_yz = buffer_1000_pddd[766];

    auto g_y_0_0_0_x_yy_xy_zz = buffer_1000_pddd[767];

    auto g_y_0_0_0_x_yy_xz_xx = buffer_1000_pddd[768];

    auto g_y_0_0_0_x_yy_xz_xy = buffer_1000_pddd[769];

    auto g_y_0_0_0_x_yy_xz_xz = buffer_1000_pddd[770];

    auto g_y_0_0_0_x_yy_xz_yy = buffer_1000_pddd[771];

    auto g_y_0_0_0_x_yy_xz_yz = buffer_1000_pddd[772];

    auto g_y_0_0_0_x_yy_xz_zz = buffer_1000_pddd[773];

    auto g_y_0_0_0_x_yy_yy_xx = buffer_1000_pddd[774];

    auto g_y_0_0_0_x_yy_yy_xy = buffer_1000_pddd[775];

    auto g_y_0_0_0_x_yy_yy_xz = buffer_1000_pddd[776];

    auto g_y_0_0_0_x_yy_yy_yy = buffer_1000_pddd[777];

    auto g_y_0_0_0_x_yy_yy_yz = buffer_1000_pddd[778];

    auto g_y_0_0_0_x_yy_yy_zz = buffer_1000_pddd[779];

    auto g_y_0_0_0_x_yy_yz_xx = buffer_1000_pddd[780];

    auto g_y_0_0_0_x_yy_yz_xy = buffer_1000_pddd[781];

    auto g_y_0_0_0_x_yy_yz_xz = buffer_1000_pddd[782];

    auto g_y_0_0_0_x_yy_yz_yy = buffer_1000_pddd[783];

    auto g_y_0_0_0_x_yy_yz_yz = buffer_1000_pddd[784];

    auto g_y_0_0_0_x_yy_yz_zz = buffer_1000_pddd[785];

    auto g_y_0_0_0_x_yy_zz_xx = buffer_1000_pddd[786];

    auto g_y_0_0_0_x_yy_zz_xy = buffer_1000_pddd[787];

    auto g_y_0_0_0_x_yy_zz_xz = buffer_1000_pddd[788];

    auto g_y_0_0_0_x_yy_zz_yy = buffer_1000_pddd[789];

    auto g_y_0_0_0_x_yy_zz_yz = buffer_1000_pddd[790];

    auto g_y_0_0_0_x_yy_zz_zz = buffer_1000_pddd[791];

    auto g_y_0_0_0_x_yz_xx_xx = buffer_1000_pddd[792];

    auto g_y_0_0_0_x_yz_xx_xy = buffer_1000_pddd[793];

    auto g_y_0_0_0_x_yz_xx_xz = buffer_1000_pddd[794];

    auto g_y_0_0_0_x_yz_xx_yy = buffer_1000_pddd[795];

    auto g_y_0_0_0_x_yz_xx_yz = buffer_1000_pddd[796];

    auto g_y_0_0_0_x_yz_xx_zz = buffer_1000_pddd[797];

    auto g_y_0_0_0_x_yz_xy_xx = buffer_1000_pddd[798];

    auto g_y_0_0_0_x_yz_xy_xy = buffer_1000_pddd[799];

    auto g_y_0_0_0_x_yz_xy_xz = buffer_1000_pddd[800];

    auto g_y_0_0_0_x_yz_xy_yy = buffer_1000_pddd[801];

    auto g_y_0_0_0_x_yz_xy_yz = buffer_1000_pddd[802];

    auto g_y_0_0_0_x_yz_xy_zz = buffer_1000_pddd[803];

    auto g_y_0_0_0_x_yz_xz_xx = buffer_1000_pddd[804];

    auto g_y_0_0_0_x_yz_xz_xy = buffer_1000_pddd[805];

    auto g_y_0_0_0_x_yz_xz_xz = buffer_1000_pddd[806];

    auto g_y_0_0_0_x_yz_xz_yy = buffer_1000_pddd[807];

    auto g_y_0_0_0_x_yz_xz_yz = buffer_1000_pddd[808];

    auto g_y_0_0_0_x_yz_xz_zz = buffer_1000_pddd[809];

    auto g_y_0_0_0_x_yz_yy_xx = buffer_1000_pddd[810];

    auto g_y_0_0_0_x_yz_yy_xy = buffer_1000_pddd[811];

    auto g_y_0_0_0_x_yz_yy_xz = buffer_1000_pddd[812];

    auto g_y_0_0_0_x_yz_yy_yy = buffer_1000_pddd[813];

    auto g_y_0_0_0_x_yz_yy_yz = buffer_1000_pddd[814];

    auto g_y_0_0_0_x_yz_yy_zz = buffer_1000_pddd[815];

    auto g_y_0_0_0_x_yz_yz_xx = buffer_1000_pddd[816];

    auto g_y_0_0_0_x_yz_yz_xy = buffer_1000_pddd[817];

    auto g_y_0_0_0_x_yz_yz_xz = buffer_1000_pddd[818];

    auto g_y_0_0_0_x_yz_yz_yy = buffer_1000_pddd[819];

    auto g_y_0_0_0_x_yz_yz_yz = buffer_1000_pddd[820];

    auto g_y_0_0_0_x_yz_yz_zz = buffer_1000_pddd[821];

    auto g_y_0_0_0_x_yz_zz_xx = buffer_1000_pddd[822];

    auto g_y_0_0_0_x_yz_zz_xy = buffer_1000_pddd[823];

    auto g_y_0_0_0_x_yz_zz_xz = buffer_1000_pddd[824];

    auto g_y_0_0_0_x_yz_zz_yy = buffer_1000_pddd[825];

    auto g_y_0_0_0_x_yz_zz_yz = buffer_1000_pddd[826];

    auto g_y_0_0_0_x_yz_zz_zz = buffer_1000_pddd[827];

    auto g_y_0_0_0_x_zz_xx_xx = buffer_1000_pddd[828];

    auto g_y_0_0_0_x_zz_xx_xy = buffer_1000_pddd[829];

    auto g_y_0_0_0_x_zz_xx_xz = buffer_1000_pddd[830];

    auto g_y_0_0_0_x_zz_xx_yy = buffer_1000_pddd[831];

    auto g_y_0_0_0_x_zz_xx_yz = buffer_1000_pddd[832];

    auto g_y_0_0_0_x_zz_xx_zz = buffer_1000_pddd[833];

    auto g_y_0_0_0_x_zz_xy_xx = buffer_1000_pddd[834];

    auto g_y_0_0_0_x_zz_xy_xy = buffer_1000_pddd[835];

    auto g_y_0_0_0_x_zz_xy_xz = buffer_1000_pddd[836];

    auto g_y_0_0_0_x_zz_xy_yy = buffer_1000_pddd[837];

    auto g_y_0_0_0_x_zz_xy_yz = buffer_1000_pddd[838];

    auto g_y_0_0_0_x_zz_xy_zz = buffer_1000_pddd[839];

    auto g_y_0_0_0_x_zz_xz_xx = buffer_1000_pddd[840];

    auto g_y_0_0_0_x_zz_xz_xy = buffer_1000_pddd[841];

    auto g_y_0_0_0_x_zz_xz_xz = buffer_1000_pddd[842];

    auto g_y_0_0_0_x_zz_xz_yy = buffer_1000_pddd[843];

    auto g_y_0_0_0_x_zz_xz_yz = buffer_1000_pddd[844];

    auto g_y_0_0_0_x_zz_xz_zz = buffer_1000_pddd[845];

    auto g_y_0_0_0_x_zz_yy_xx = buffer_1000_pddd[846];

    auto g_y_0_0_0_x_zz_yy_xy = buffer_1000_pddd[847];

    auto g_y_0_0_0_x_zz_yy_xz = buffer_1000_pddd[848];

    auto g_y_0_0_0_x_zz_yy_yy = buffer_1000_pddd[849];

    auto g_y_0_0_0_x_zz_yy_yz = buffer_1000_pddd[850];

    auto g_y_0_0_0_x_zz_yy_zz = buffer_1000_pddd[851];

    auto g_y_0_0_0_x_zz_yz_xx = buffer_1000_pddd[852];

    auto g_y_0_0_0_x_zz_yz_xy = buffer_1000_pddd[853];

    auto g_y_0_0_0_x_zz_yz_xz = buffer_1000_pddd[854];

    auto g_y_0_0_0_x_zz_yz_yy = buffer_1000_pddd[855];

    auto g_y_0_0_0_x_zz_yz_yz = buffer_1000_pddd[856];

    auto g_y_0_0_0_x_zz_yz_zz = buffer_1000_pddd[857];

    auto g_y_0_0_0_x_zz_zz_xx = buffer_1000_pddd[858];

    auto g_y_0_0_0_x_zz_zz_xy = buffer_1000_pddd[859];

    auto g_y_0_0_0_x_zz_zz_xz = buffer_1000_pddd[860];

    auto g_y_0_0_0_x_zz_zz_yy = buffer_1000_pddd[861];

    auto g_y_0_0_0_x_zz_zz_yz = buffer_1000_pddd[862];

    auto g_y_0_0_0_x_zz_zz_zz = buffer_1000_pddd[863];

    auto g_y_0_0_0_y_xx_xx_xx = buffer_1000_pddd[864];

    auto g_y_0_0_0_y_xx_xx_xy = buffer_1000_pddd[865];

    auto g_y_0_0_0_y_xx_xx_xz = buffer_1000_pddd[866];

    auto g_y_0_0_0_y_xx_xx_yy = buffer_1000_pddd[867];

    auto g_y_0_0_0_y_xx_xx_yz = buffer_1000_pddd[868];

    auto g_y_0_0_0_y_xx_xx_zz = buffer_1000_pddd[869];

    auto g_y_0_0_0_y_xx_xy_xx = buffer_1000_pddd[870];

    auto g_y_0_0_0_y_xx_xy_xy = buffer_1000_pddd[871];

    auto g_y_0_0_0_y_xx_xy_xz = buffer_1000_pddd[872];

    auto g_y_0_0_0_y_xx_xy_yy = buffer_1000_pddd[873];

    auto g_y_0_0_0_y_xx_xy_yz = buffer_1000_pddd[874];

    auto g_y_0_0_0_y_xx_xy_zz = buffer_1000_pddd[875];

    auto g_y_0_0_0_y_xx_xz_xx = buffer_1000_pddd[876];

    auto g_y_0_0_0_y_xx_xz_xy = buffer_1000_pddd[877];

    auto g_y_0_0_0_y_xx_xz_xz = buffer_1000_pddd[878];

    auto g_y_0_0_0_y_xx_xz_yy = buffer_1000_pddd[879];

    auto g_y_0_0_0_y_xx_xz_yz = buffer_1000_pddd[880];

    auto g_y_0_0_0_y_xx_xz_zz = buffer_1000_pddd[881];

    auto g_y_0_0_0_y_xx_yy_xx = buffer_1000_pddd[882];

    auto g_y_0_0_0_y_xx_yy_xy = buffer_1000_pddd[883];

    auto g_y_0_0_0_y_xx_yy_xz = buffer_1000_pddd[884];

    auto g_y_0_0_0_y_xx_yy_yy = buffer_1000_pddd[885];

    auto g_y_0_0_0_y_xx_yy_yz = buffer_1000_pddd[886];

    auto g_y_0_0_0_y_xx_yy_zz = buffer_1000_pddd[887];

    auto g_y_0_0_0_y_xx_yz_xx = buffer_1000_pddd[888];

    auto g_y_0_0_0_y_xx_yz_xy = buffer_1000_pddd[889];

    auto g_y_0_0_0_y_xx_yz_xz = buffer_1000_pddd[890];

    auto g_y_0_0_0_y_xx_yz_yy = buffer_1000_pddd[891];

    auto g_y_0_0_0_y_xx_yz_yz = buffer_1000_pddd[892];

    auto g_y_0_0_0_y_xx_yz_zz = buffer_1000_pddd[893];

    auto g_y_0_0_0_y_xx_zz_xx = buffer_1000_pddd[894];

    auto g_y_0_0_0_y_xx_zz_xy = buffer_1000_pddd[895];

    auto g_y_0_0_0_y_xx_zz_xz = buffer_1000_pddd[896];

    auto g_y_0_0_0_y_xx_zz_yy = buffer_1000_pddd[897];

    auto g_y_0_0_0_y_xx_zz_yz = buffer_1000_pddd[898];

    auto g_y_0_0_0_y_xx_zz_zz = buffer_1000_pddd[899];

    auto g_y_0_0_0_y_xy_xx_xx = buffer_1000_pddd[900];

    auto g_y_0_0_0_y_xy_xx_xy = buffer_1000_pddd[901];

    auto g_y_0_0_0_y_xy_xx_xz = buffer_1000_pddd[902];

    auto g_y_0_0_0_y_xy_xx_yy = buffer_1000_pddd[903];

    auto g_y_0_0_0_y_xy_xx_yz = buffer_1000_pddd[904];

    auto g_y_0_0_0_y_xy_xx_zz = buffer_1000_pddd[905];

    auto g_y_0_0_0_y_xy_xy_xx = buffer_1000_pddd[906];

    auto g_y_0_0_0_y_xy_xy_xy = buffer_1000_pddd[907];

    auto g_y_0_0_0_y_xy_xy_xz = buffer_1000_pddd[908];

    auto g_y_0_0_0_y_xy_xy_yy = buffer_1000_pddd[909];

    auto g_y_0_0_0_y_xy_xy_yz = buffer_1000_pddd[910];

    auto g_y_0_0_0_y_xy_xy_zz = buffer_1000_pddd[911];

    auto g_y_0_0_0_y_xy_xz_xx = buffer_1000_pddd[912];

    auto g_y_0_0_0_y_xy_xz_xy = buffer_1000_pddd[913];

    auto g_y_0_0_0_y_xy_xz_xz = buffer_1000_pddd[914];

    auto g_y_0_0_0_y_xy_xz_yy = buffer_1000_pddd[915];

    auto g_y_0_0_0_y_xy_xz_yz = buffer_1000_pddd[916];

    auto g_y_0_0_0_y_xy_xz_zz = buffer_1000_pddd[917];

    auto g_y_0_0_0_y_xy_yy_xx = buffer_1000_pddd[918];

    auto g_y_0_0_0_y_xy_yy_xy = buffer_1000_pddd[919];

    auto g_y_0_0_0_y_xy_yy_xz = buffer_1000_pddd[920];

    auto g_y_0_0_0_y_xy_yy_yy = buffer_1000_pddd[921];

    auto g_y_0_0_0_y_xy_yy_yz = buffer_1000_pddd[922];

    auto g_y_0_0_0_y_xy_yy_zz = buffer_1000_pddd[923];

    auto g_y_0_0_0_y_xy_yz_xx = buffer_1000_pddd[924];

    auto g_y_0_0_0_y_xy_yz_xy = buffer_1000_pddd[925];

    auto g_y_0_0_0_y_xy_yz_xz = buffer_1000_pddd[926];

    auto g_y_0_0_0_y_xy_yz_yy = buffer_1000_pddd[927];

    auto g_y_0_0_0_y_xy_yz_yz = buffer_1000_pddd[928];

    auto g_y_0_0_0_y_xy_yz_zz = buffer_1000_pddd[929];

    auto g_y_0_0_0_y_xy_zz_xx = buffer_1000_pddd[930];

    auto g_y_0_0_0_y_xy_zz_xy = buffer_1000_pddd[931];

    auto g_y_0_0_0_y_xy_zz_xz = buffer_1000_pddd[932];

    auto g_y_0_0_0_y_xy_zz_yy = buffer_1000_pddd[933];

    auto g_y_0_0_0_y_xy_zz_yz = buffer_1000_pddd[934];

    auto g_y_0_0_0_y_xy_zz_zz = buffer_1000_pddd[935];

    auto g_y_0_0_0_y_xz_xx_xx = buffer_1000_pddd[936];

    auto g_y_0_0_0_y_xz_xx_xy = buffer_1000_pddd[937];

    auto g_y_0_0_0_y_xz_xx_xz = buffer_1000_pddd[938];

    auto g_y_0_0_0_y_xz_xx_yy = buffer_1000_pddd[939];

    auto g_y_0_0_0_y_xz_xx_yz = buffer_1000_pddd[940];

    auto g_y_0_0_0_y_xz_xx_zz = buffer_1000_pddd[941];

    auto g_y_0_0_0_y_xz_xy_xx = buffer_1000_pddd[942];

    auto g_y_0_0_0_y_xz_xy_xy = buffer_1000_pddd[943];

    auto g_y_0_0_0_y_xz_xy_xz = buffer_1000_pddd[944];

    auto g_y_0_0_0_y_xz_xy_yy = buffer_1000_pddd[945];

    auto g_y_0_0_0_y_xz_xy_yz = buffer_1000_pddd[946];

    auto g_y_0_0_0_y_xz_xy_zz = buffer_1000_pddd[947];

    auto g_y_0_0_0_y_xz_xz_xx = buffer_1000_pddd[948];

    auto g_y_0_0_0_y_xz_xz_xy = buffer_1000_pddd[949];

    auto g_y_0_0_0_y_xz_xz_xz = buffer_1000_pddd[950];

    auto g_y_0_0_0_y_xz_xz_yy = buffer_1000_pddd[951];

    auto g_y_0_0_0_y_xz_xz_yz = buffer_1000_pddd[952];

    auto g_y_0_0_0_y_xz_xz_zz = buffer_1000_pddd[953];

    auto g_y_0_0_0_y_xz_yy_xx = buffer_1000_pddd[954];

    auto g_y_0_0_0_y_xz_yy_xy = buffer_1000_pddd[955];

    auto g_y_0_0_0_y_xz_yy_xz = buffer_1000_pddd[956];

    auto g_y_0_0_0_y_xz_yy_yy = buffer_1000_pddd[957];

    auto g_y_0_0_0_y_xz_yy_yz = buffer_1000_pddd[958];

    auto g_y_0_0_0_y_xz_yy_zz = buffer_1000_pddd[959];

    auto g_y_0_0_0_y_xz_yz_xx = buffer_1000_pddd[960];

    auto g_y_0_0_0_y_xz_yz_xy = buffer_1000_pddd[961];

    auto g_y_0_0_0_y_xz_yz_xz = buffer_1000_pddd[962];

    auto g_y_0_0_0_y_xz_yz_yy = buffer_1000_pddd[963];

    auto g_y_0_0_0_y_xz_yz_yz = buffer_1000_pddd[964];

    auto g_y_0_0_0_y_xz_yz_zz = buffer_1000_pddd[965];

    auto g_y_0_0_0_y_xz_zz_xx = buffer_1000_pddd[966];

    auto g_y_0_0_0_y_xz_zz_xy = buffer_1000_pddd[967];

    auto g_y_0_0_0_y_xz_zz_xz = buffer_1000_pddd[968];

    auto g_y_0_0_0_y_xz_zz_yy = buffer_1000_pddd[969];

    auto g_y_0_0_0_y_xz_zz_yz = buffer_1000_pddd[970];

    auto g_y_0_0_0_y_xz_zz_zz = buffer_1000_pddd[971];

    auto g_y_0_0_0_y_yy_xx_xx = buffer_1000_pddd[972];

    auto g_y_0_0_0_y_yy_xx_xy = buffer_1000_pddd[973];

    auto g_y_0_0_0_y_yy_xx_xz = buffer_1000_pddd[974];

    auto g_y_0_0_0_y_yy_xx_yy = buffer_1000_pddd[975];

    auto g_y_0_0_0_y_yy_xx_yz = buffer_1000_pddd[976];

    auto g_y_0_0_0_y_yy_xx_zz = buffer_1000_pddd[977];

    auto g_y_0_0_0_y_yy_xy_xx = buffer_1000_pddd[978];

    auto g_y_0_0_0_y_yy_xy_xy = buffer_1000_pddd[979];

    auto g_y_0_0_0_y_yy_xy_xz = buffer_1000_pddd[980];

    auto g_y_0_0_0_y_yy_xy_yy = buffer_1000_pddd[981];

    auto g_y_0_0_0_y_yy_xy_yz = buffer_1000_pddd[982];

    auto g_y_0_0_0_y_yy_xy_zz = buffer_1000_pddd[983];

    auto g_y_0_0_0_y_yy_xz_xx = buffer_1000_pddd[984];

    auto g_y_0_0_0_y_yy_xz_xy = buffer_1000_pddd[985];

    auto g_y_0_0_0_y_yy_xz_xz = buffer_1000_pddd[986];

    auto g_y_0_0_0_y_yy_xz_yy = buffer_1000_pddd[987];

    auto g_y_0_0_0_y_yy_xz_yz = buffer_1000_pddd[988];

    auto g_y_0_0_0_y_yy_xz_zz = buffer_1000_pddd[989];

    auto g_y_0_0_0_y_yy_yy_xx = buffer_1000_pddd[990];

    auto g_y_0_0_0_y_yy_yy_xy = buffer_1000_pddd[991];

    auto g_y_0_0_0_y_yy_yy_xz = buffer_1000_pddd[992];

    auto g_y_0_0_0_y_yy_yy_yy = buffer_1000_pddd[993];

    auto g_y_0_0_0_y_yy_yy_yz = buffer_1000_pddd[994];

    auto g_y_0_0_0_y_yy_yy_zz = buffer_1000_pddd[995];

    auto g_y_0_0_0_y_yy_yz_xx = buffer_1000_pddd[996];

    auto g_y_0_0_0_y_yy_yz_xy = buffer_1000_pddd[997];

    auto g_y_0_0_0_y_yy_yz_xz = buffer_1000_pddd[998];

    auto g_y_0_0_0_y_yy_yz_yy = buffer_1000_pddd[999];

    auto g_y_0_0_0_y_yy_yz_yz = buffer_1000_pddd[1000];

    auto g_y_0_0_0_y_yy_yz_zz = buffer_1000_pddd[1001];

    auto g_y_0_0_0_y_yy_zz_xx = buffer_1000_pddd[1002];

    auto g_y_0_0_0_y_yy_zz_xy = buffer_1000_pddd[1003];

    auto g_y_0_0_0_y_yy_zz_xz = buffer_1000_pddd[1004];

    auto g_y_0_0_0_y_yy_zz_yy = buffer_1000_pddd[1005];

    auto g_y_0_0_0_y_yy_zz_yz = buffer_1000_pddd[1006];

    auto g_y_0_0_0_y_yy_zz_zz = buffer_1000_pddd[1007];

    auto g_y_0_0_0_y_yz_xx_xx = buffer_1000_pddd[1008];

    auto g_y_0_0_0_y_yz_xx_xy = buffer_1000_pddd[1009];

    auto g_y_0_0_0_y_yz_xx_xz = buffer_1000_pddd[1010];

    auto g_y_0_0_0_y_yz_xx_yy = buffer_1000_pddd[1011];

    auto g_y_0_0_0_y_yz_xx_yz = buffer_1000_pddd[1012];

    auto g_y_0_0_0_y_yz_xx_zz = buffer_1000_pddd[1013];

    auto g_y_0_0_0_y_yz_xy_xx = buffer_1000_pddd[1014];

    auto g_y_0_0_0_y_yz_xy_xy = buffer_1000_pddd[1015];

    auto g_y_0_0_0_y_yz_xy_xz = buffer_1000_pddd[1016];

    auto g_y_0_0_0_y_yz_xy_yy = buffer_1000_pddd[1017];

    auto g_y_0_0_0_y_yz_xy_yz = buffer_1000_pddd[1018];

    auto g_y_0_0_0_y_yz_xy_zz = buffer_1000_pddd[1019];

    auto g_y_0_0_0_y_yz_xz_xx = buffer_1000_pddd[1020];

    auto g_y_0_0_0_y_yz_xz_xy = buffer_1000_pddd[1021];

    auto g_y_0_0_0_y_yz_xz_xz = buffer_1000_pddd[1022];

    auto g_y_0_0_0_y_yz_xz_yy = buffer_1000_pddd[1023];

    auto g_y_0_0_0_y_yz_xz_yz = buffer_1000_pddd[1024];

    auto g_y_0_0_0_y_yz_xz_zz = buffer_1000_pddd[1025];

    auto g_y_0_0_0_y_yz_yy_xx = buffer_1000_pddd[1026];

    auto g_y_0_0_0_y_yz_yy_xy = buffer_1000_pddd[1027];

    auto g_y_0_0_0_y_yz_yy_xz = buffer_1000_pddd[1028];

    auto g_y_0_0_0_y_yz_yy_yy = buffer_1000_pddd[1029];

    auto g_y_0_0_0_y_yz_yy_yz = buffer_1000_pddd[1030];

    auto g_y_0_0_0_y_yz_yy_zz = buffer_1000_pddd[1031];

    auto g_y_0_0_0_y_yz_yz_xx = buffer_1000_pddd[1032];

    auto g_y_0_0_0_y_yz_yz_xy = buffer_1000_pddd[1033];

    auto g_y_0_0_0_y_yz_yz_xz = buffer_1000_pddd[1034];

    auto g_y_0_0_0_y_yz_yz_yy = buffer_1000_pddd[1035];

    auto g_y_0_0_0_y_yz_yz_yz = buffer_1000_pddd[1036];

    auto g_y_0_0_0_y_yz_yz_zz = buffer_1000_pddd[1037];

    auto g_y_0_0_0_y_yz_zz_xx = buffer_1000_pddd[1038];

    auto g_y_0_0_0_y_yz_zz_xy = buffer_1000_pddd[1039];

    auto g_y_0_0_0_y_yz_zz_xz = buffer_1000_pddd[1040];

    auto g_y_0_0_0_y_yz_zz_yy = buffer_1000_pddd[1041];

    auto g_y_0_0_0_y_yz_zz_yz = buffer_1000_pddd[1042];

    auto g_y_0_0_0_y_yz_zz_zz = buffer_1000_pddd[1043];

    auto g_y_0_0_0_y_zz_xx_xx = buffer_1000_pddd[1044];

    auto g_y_0_0_0_y_zz_xx_xy = buffer_1000_pddd[1045];

    auto g_y_0_0_0_y_zz_xx_xz = buffer_1000_pddd[1046];

    auto g_y_0_0_0_y_zz_xx_yy = buffer_1000_pddd[1047];

    auto g_y_0_0_0_y_zz_xx_yz = buffer_1000_pddd[1048];

    auto g_y_0_0_0_y_zz_xx_zz = buffer_1000_pddd[1049];

    auto g_y_0_0_0_y_zz_xy_xx = buffer_1000_pddd[1050];

    auto g_y_0_0_0_y_zz_xy_xy = buffer_1000_pddd[1051];

    auto g_y_0_0_0_y_zz_xy_xz = buffer_1000_pddd[1052];

    auto g_y_0_0_0_y_zz_xy_yy = buffer_1000_pddd[1053];

    auto g_y_0_0_0_y_zz_xy_yz = buffer_1000_pddd[1054];

    auto g_y_0_0_0_y_zz_xy_zz = buffer_1000_pddd[1055];

    auto g_y_0_0_0_y_zz_xz_xx = buffer_1000_pddd[1056];

    auto g_y_0_0_0_y_zz_xz_xy = buffer_1000_pddd[1057];

    auto g_y_0_0_0_y_zz_xz_xz = buffer_1000_pddd[1058];

    auto g_y_0_0_0_y_zz_xz_yy = buffer_1000_pddd[1059];

    auto g_y_0_0_0_y_zz_xz_yz = buffer_1000_pddd[1060];

    auto g_y_0_0_0_y_zz_xz_zz = buffer_1000_pddd[1061];

    auto g_y_0_0_0_y_zz_yy_xx = buffer_1000_pddd[1062];

    auto g_y_0_0_0_y_zz_yy_xy = buffer_1000_pddd[1063];

    auto g_y_0_0_0_y_zz_yy_xz = buffer_1000_pddd[1064];

    auto g_y_0_0_0_y_zz_yy_yy = buffer_1000_pddd[1065];

    auto g_y_0_0_0_y_zz_yy_yz = buffer_1000_pddd[1066];

    auto g_y_0_0_0_y_zz_yy_zz = buffer_1000_pddd[1067];

    auto g_y_0_0_0_y_zz_yz_xx = buffer_1000_pddd[1068];

    auto g_y_0_0_0_y_zz_yz_xy = buffer_1000_pddd[1069];

    auto g_y_0_0_0_y_zz_yz_xz = buffer_1000_pddd[1070];

    auto g_y_0_0_0_y_zz_yz_yy = buffer_1000_pddd[1071];

    auto g_y_0_0_0_y_zz_yz_yz = buffer_1000_pddd[1072];

    auto g_y_0_0_0_y_zz_yz_zz = buffer_1000_pddd[1073];

    auto g_y_0_0_0_y_zz_zz_xx = buffer_1000_pddd[1074];

    auto g_y_0_0_0_y_zz_zz_xy = buffer_1000_pddd[1075];

    auto g_y_0_0_0_y_zz_zz_xz = buffer_1000_pddd[1076];

    auto g_y_0_0_0_y_zz_zz_yy = buffer_1000_pddd[1077];

    auto g_y_0_0_0_y_zz_zz_yz = buffer_1000_pddd[1078];

    auto g_y_0_0_0_y_zz_zz_zz = buffer_1000_pddd[1079];

    auto g_y_0_0_0_z_xx_xx_xx = buffer_1000_pddd[1080];

    auto g_y_0_0_0_z_xx_xx_xy = buffer_1000_pddd[1081];

    auto g_y_0_0_0_z_xx_xx_xz = buffer_1000_pddd[1082];

    auto g_y_0_0_0_z_xx_xx_yy = buffer_1000_pddd[1083];

    auto g_y_0_0_0_z_xx_xx_yz = buffer_1000_pddd[1084];

    auto g_y_0_0_0_z_xx_xx_zz = buffer_1000_pddd[1085];

    auto g_y_0_0_0_z_xx_xy_xx = buffer_1000_pddd[1086];

    auto g_y_0_0_0_z_xx_xy_xy = buffer_1000_pddd[1087];

    auto g_y_0_0_0_z_xx_xy_xz = buffer_1000_pddd[1088];

    auto g_y_0_0_0_z_xx_xy_yy = buffer_1000_pddd[1089];

    auto g_y_0_0_0_z_xx_xy_yz = buffer_1000_pddd[1090];

    auto g_y_0_0_0_z_xx_xy_zz = buffer_1000_pddd[1091];

    auto g_y_0_0_0_z_xx_xz_xx = buffer_1000_pddd[1092];

    auto g_y_0_0_0_z_xx_xz_xy = buffer_1000_pddd[1093];

    auto g_y_0_0_0_z_xx_xz_xz = buffer_1000_pddd[1094];

    auto g_y_0_0_0_z_xx_xz_yy = buffer_1000_pddd[1095];

    auto g_y_0_0_0_z_xx_xz_yz = buffer_1000_pddd[1096];

    auto g_y_0_0_0_z_xx_xz_zz = buffer_1000_pddd[1097];

    auto g_y_0_0_0_z_xx_yy_xx = buffer_1000_pddd[1098];

    auto g_y_0_0_0_z_xx_yy_xy = buffer_1000_pddd[1099];

    auto g_y_0_0_0_z_xx_yy_xz = buffer_1000_pddd[1100];

    auto g_y_0_0_0_z_xx_yy_yy = buffer_1000_pddd[1101];

    auto g_y_0_0_0_z_xx_yy_yz = buffer_1000_pddd[1102];

    auto g_y_0_0_0_z_xx_yy_zz = buffer_1000_pddd[1103];

    auto g_y_0_0_0_z_xx_yz_xx = buffer_1000_pddd[1104];

    auto g_y_0_0_0_z_xx_yz_xy = buffer_1000_pddd[1105];

    auto g_y_0_0_0_z_xx_yz_xz = buffer_1000_pddd[1106];

    auto g_y_0_0_0_z_xx_yz_yy = buffer_1000_pddd[1107];

    auto g_y_0_0_0_z_xx_yz_yz = buffer_1000_pddd[1108];

    auto g_y_0_0_0_z_xx_yz_zz = buffer_1000_pddd[1109];

    auto g_y_0_0_0_z_xx_zz_xx = buffer_1000_pddd[1110];

    auto g_y_0_0_0_z_xx_zz_xy = buffer_1000_pddd[1111];

    auto g_y_0_0_0_z_xx_zz_xz = buffer_1000_pddd[1112];

    auto g_y_0_0_0_z_xx_zz_yy = buffer_1000_pddd[1113];

    auto g_y_0_0_0_z_xx_zz_yz = buffer_1000_pddd[1114];

    auto g_y_0_0_0_z_xx_zz_zz = buffer_1000_pddd[1115];

    auto g_y_0_0_0_z_xy_xx_xx = buffer_1000_pddd[1116];

    auto g_y_0_0_0_z_xy_xx_xy = buffer_1000_pddd[1117];

    auto g_y_0_0_0_z_xy_xx_xz = buffer_1000_pddd[1118];

    auto g_y_0_0_0_z_xy_xx_yy = buffer_1000_pddd[1119];

    auto g_y_0_0_0_z_xy_xx_yz = buffer_1000_pddd[1120];

    auto g_y_0_0_0_z_xy_xx_zz = buffer_1000_pddd[1121];

    auto g_y_0_0_0_z_xy_xy_xx = buffer_1000_pddd[1122];

    auto g_y_0_0_0_z_xy_xy_xy = buffer_1000_pddd[1123];

    auto g_y_0_0_0_z_xy_xy_xz = buffer_1000_pddd[1124];

    auto g_y_0_0_0_z_xy_xy_yy = buffer_1000_pddd[1125];

    auto g_y_0_0_0_z_xy_xy_yz = buffer_1000_pddd[1126];

    auto g_y_0_0_0_z_xy_xy_zz = buffer_1000_pddd[1127];

    auto g_y_0_0_0_z_xy_xz_xx = buffer_1000_pddd[1128];

    auto g_y_0_0_0_z_xy_xz_xy = buffer_1000_pddd[1129];

    auto g_y_0_0_0_z_xy_xz_xz = buffer_1000_pddd[1130];

    auto g_y_0_0_0_z_xy_xz_yy = buffer_1000_pddd[1131];

    auto g_y_0_0_0_z_xy_xz_yz = buffer_1000_pddd[1132];

    auto g_y_0_0_0_z_xy_xz_zz = buffer_1000_pddd[1133];

    auto g_y_0_0_0_z_xy_yy_xx = buffer_1000_pddd[1134];

    auto g_y_0_0_0_z_xy_yy_xy = buffer_1000_pddd[1135];

    auto g_y_0_0_0_z_xy_yy_xz = buffer_1000_pddd[1136];

    auto g_y_0_0_0_z_xy_yy_yy = buffer_1000_pddd[1137];

    auto g_y_0_0_0_z_xy_yy_yz = buffer_1000_pddd[1138];

    auto g_y_0_0_0_z_xy_yy_zz = buffer_1000_pddd[1139];

    auto g_y_0_0_0_z_xy_yz_xx = buffer_1000_pddd[1140];

    auto g_y_0_0_0_z_xy_yz_xy = buffer_1000_pddd[1141];

    auto g_y_0_0_0_z_xy_yz_xz = buffer_1000_pddd[1142];

    auto g_y_0_0_0_z_xy_yz_yy = buffer_1000_pddd[1143];

    auto g_y_0_0_0_z_xy_yz_yz = buffer_1000_pddd[1144];

    auto g_y_0_0_0_z_xy_yz_zz = buffer_1000_pddd[1145];

    auto g_y_0_0_0_z_xy_zz_xx = buffer_1000_pddd[1146];

    auto g_y_0_0_0_z_xy_zz_xy = buffer_1000_pddd[1147];

    auto g_y_0_0_0_z_xy_zz_xz = buffer_1000_pddd[1148];

    auto g_y_0_0_0_z_xy_zz_yy = buffer_1000_pddd[1149];

    auto g_y_0_0_0_z_xy_zz_yz = buffer_1000_pddd[1150];

    auto g_y_0_0_0_z_xy_zz_zz = buffer_1000_pddd[1151];

    auto g_y_0_0_0_z_xz_xx_xx = buffer_1000_pddd[1152];

    auto g_y_0_0_0_z_xz_xx_xy = buffer_1000_pddd[1153];

    auto g_y_0_0_0_z_xz_xx_xz = buffer_1000_pddd[1154];

    auto g_y_0_0_0_z_xz_xx_yy = buffer_1000_pddd[1155];

    auto g_y_0_0_0_z_xz_xx_yz = buffer_1000_pddd[1156];

    auto g_y_0_0_0_z_xz_xx_zz = buffer_1000_pddd[1157];

    auto g_y_0_0_0_z_xz_xy_xx = buffer_1000_pddd[1158];

    auto g_y_0_0_0_z_xz_xy_xy = buffer_1000_pddd[1159];

    auto g_y_0_0_0_z_xz_xy_xz = buffer_1000_pddd[1160];

    auto g_y_0_0_0_z_xz_xy_yy = buffer_1000_pddd[1161];

    auto g_y_0_0_0_z_xz_xy_yz = buffer_1000_pddd[1162];

    auto g_y_0_0_0_z_xz_xy_zz = buffer_1000_pddd[1163];

    auto g_y_0_0_0_z_xz_xz_xx = buffer_1000_pddd[1164];

    auto g_y_0_0_0_z_xz_xz_xy = buffer_1000_pddd[1165];

    auto g_y_0_0_0_z_xz_xz_xz = buffer_1000_pddd[1166];

    auto g_y_0_0_0_z_xz_xz_yy = buffer_1000_pddd[1167];

    auto g_y_0_0_0_z_xz_xz_yz = buffer_1000_pddd[1168];

    auto g_y_0_0_0_z_xz_xz_zz = buffer_1000_pddd[1169];

    auto g_y_0_0_0_z_xz_yy_xx = buffer_1000_pddd[1170];

    auto g_y_0_0_0_z_xz_yy_xy = buffer_1000_pddd[1171];

    auto g_y_0_0_0_z_xz_yy_xz = buffer_1000_pddd[1172];

    auto g_y_0_0_0_z_xz_yy_yy = buffer_1000_pddd[1173];

    auto g_y_0_0_0_z_xz_yy_yz = buffer_1000_pddd[1174];

    auto g_y_0_0_0_z_xz_yy_zz = buffer_1000_pddd[1175];

    auto g_y_0_0_0_z_xz_yz_xx = buffer_1000_pddd[1176];

    auto g_y_0_0_0_z_xz_yz_xy = buffer_1000_pddd[1177];

    auto g_y_0_0_0_z_xz_yz_xz = buffer_1000_pddd[1178];

    auto g_y_0_0_0_z_xz_yz_yy = buffer_1000_pddd[1179];

    auto g_y_0_0_0_z_xz_yz_yz = buffer_1000_pddd[1180];

    auto g_y_0_0_0_z_xz_yz_zz = buffer_1000_pddd[1181];

    auto g_y_0_0_0_z_xz_zz_xx = buffer_1000_pddd[1182];

    auto g_y_0_0_0_z_xz_zz_xy = buffer_1000_pddd[1183];

    auto g_y_0_0_0_z_xz_zz_xz = buffer_1000_pddd[1184];

    auto g_y_0_0_0_z_xz_zz_yy = buffer_1000_pddd[1185];

    auto g_y_0_0_0_z_xz_zz_yz = buffer_1000_pddd[1186];

    auto g_y_0_0_0_z_xz_zz_zz = buffer_1000_pddd[1187];

    auto g_y_0_0_0_z_yy_xx_xx = buffer_1000_pddd[1188];

    auto g_y_0_0_0_z_yy_xx_xy = buffer_1000_pddd[1189];

    auto g_y_0_0_0_z_yy_xx_xz = buffer_1000_pddd[1190];

    auto g_y_0_0_0_z_yy_xx_yy = buffer_1000_pddd[1191];

    auto g_y_0_0_0_z_yy_xx_yz = buffer_1000_pddd[1192];

    auto g_y_0_0_0_z_yy_xx_zz = buffer_1000_pddd[1193];

    auto g_y_0_0_0_z_yy_xy_xx = buffer_1000_pddd[1194];

    auto g_y_0_0_0_z_yy_xy_xy = buffer_1000_pddd[1195];

    auto g_y_0_0_0_z_yy_xy_xz = buffer_1000_pddd[1196];

    auto g_y_0_0_0_z_yy_xy_yy = buffer_1000_pddd[1197];

    auto g_y_0_0_0_z_yy_xy_yz = buffer_1000_pddd[1198];

    auto g_y_0_0_0_z_yy_xy_zz = buffer_1000_pddd[1199];

    auto g_y_0_0_0_z_yy_xz_xx = buffer_1000_pddd[1200];

    auto g_y_0_0_0_z_yy_xz_xy = buffer_1000_pddd[1201];

    auto g_y_0_0_0_z_yy_xz_xz = buffer_1000_pddd[1202];

    auto g_y_0_0_0_z_yy_xz_yy = buffer_1000_pddd[1203];

    auto g_y_0_0_0_z_yy_xz_yz = buffer_1000_pddd[1204];

    auto g_y_0_0_0_z_yy_xz_zz = buffer_1000_pddd[1205];

    auto g_y_0_0_0_z_yy_yy_xx = buffer_1000_pddd[1206];

    auto g_y_0_0_0_z_yy_yy_xy = buffer_1000_pddd[1207];

    auto g_y_0_0_0_z_yy_yy_xz = buffer_1000_pddd[1208];

    auto g_y_0_0_0_z_yy_yy_yy = buffer_1000_pddd[1209];

    auto g_y_0_0_0_z_yy_yy_yz = buffer_1000_pddd[1210];

    auto g_y_0_0_0_z_yy_yy_zz = buffer_1000_pddd[1211];

    auto g_y_0_0_0_z_yy_yz_xx = buffer_1000_pddd[1212];

    auto g_y_0_0_0_z_yy_yz_xy = buffer_1000_pddd[1213];

    auto g_y_0_0_0_z_yy_yz_xz = buffer_1000_pddd[1214];

    auto g_y_0_0_0_z_yy_yz_yy = buffer_1000_pddd[1215];

    auto g_y_0_0_0_z_yy_yz_yz = buffer_1000_pddd[1216];

    auto g_y_0_0_0_z_yy_yz_zz = buffer_1000_pddd[1217];

    auto g_y_0_0_0_z_yy_zz_xx = buffer_1000_pddd[1218];

    auto g_y_0_0_0_z_yy_zz_xy = buffer_1000_pddd[1219];

    auto g_y_0_0_0_z_yy_zz_xz = buffer_1000_pddd[1220];

    auto g_y_0_0_0_z_yy_zz_yy = buffer_1000_pddd[1221];

    auto g_y_0_0_0_z_yy_zz_yz = buffer_1000_pddd[1222];

    auto g_y_0_0_0_z_yy_zz_zz = buffer_1000_pddd[1223];

    auto g_y_0_0_0_z_yz_xx_xx = buffer_1000_pddd[1224];

    auto g_y_0_0_0_z_yz_xx_xy = buffer_1000_pddd[1225];

    auto g_y_0_0_0_z_yz_xx_xz = buffer_1000_pddd[1226];

    auto g_y_0_0_0_z_yz_xx_yy = buffer_1000_pddd[1227];

    auto g_y_0_0_0_z_yz_xx_yz = buffer_1000_pddd[1228];

    auto g_y_0_0_0_z_yz_xx_zz = buffer_1000_pddd[1229];

    auto g_y_0_0_0_z_yz_xy_xx = buffer_1000_pddd[1230];

    auto g_y_0_0_0_z_yz_xy_xy = buffer_1000_pddd[1231];

    auto g_y_0_0_0_z_yz_xy_xz = buffer_1000_pddd[1232];

    auto g_y_0_0_0_z_yz_xy_yy = buffer_1000_pddd[1233];

    auto g_y_0_0_0_z_yz_xy_yz = buffer_1000_pddd[1234];

    auto g_y_0_0_0_z_yz_xy_zz = buffer_1000_pddd[1235];

    auto g_y_0_0_0_z_yz_xz_xx = buffer_1000_pddd[1236];

    auto g_y_0_0_0_z_yz_xz_xy = buffer_1000_pddd[1237];

    auto g_y_0_0_0_z_yz_xz_xz = buffer_1000_pddd[1238];

    auto g_y_0_0_0_z_yz_xz_yy = buffer_1000_pddd[1239];

    auto g_y_0_0_0_z_yz_xz_yz = buffer_1000_pddd[1240];

    auto g_y_0_0_0_z_yz_xz_zz = buffer_1000_pddd[1241];

    auto g_y_0_0_0_z_yz_yy_xx = buffer_1000_pddd[1242];

    auto g_y_0_0_0_z_yz_yy_xy = buffer_1000_pddd[1243];

    auto g_y_0_0_0_z_yz_yy_xz = buffer_1000_pddd[1244];

    auto g_y_0_0_0_z_yz_yy_yy = buffer_1000_pddd[1245];

    auto g_y_0_0_0_z_yz_yy_yz = buffer_1000_pddd[1246];

    auto g_y_0_0_0_z_yz_yy_zz = buffer_1000_pddd[1247];

    auto g_y_0_0_0_z_yz_yz_xx = buffer_1000_pddd[1248];

    auto g_y_0_0_0_z_yz_yz_xy = buffer_1000_pddd[1249];

    auto g_y_0_0_0_z_yz_yz_xz = buffer_1000_pddd[1250];

    auto g_y_0_0_0_z_yz_yz_yy = buffer_1000_pddd[1251];

    auto g_y_0_0_0_z_yz_yz_yz = buffer_1000_pddd[1252];

    auto g_y_0_0_0_z_yz_yz_zz = buffer_1000_pddd[1253];

    auto g_y_0_0_0_z_yz_zz_xx = buffer_1000_pddd[1254];

    auto g_y_0_0_0_z_yz_zz_xy = buffer_1000_pddd[1255];

    auto g_y_0_0_0_z_yz_zz_xz = buffer_1000_pddd[1256];

    auto g_y_0_0_0_z_yz_zz_yy = buffer_1000_pddd[1257];

    auto g_y_0_0_0_z_yz_zz_yz = buffer_1000_pddd[1258];

    auto g_y_0_0_0_z_yz_zz_zz = buffer_1000_pddd[1259];

    auto g_y_0_0_0_z_zz_xx_xx = buffer_1000_pddd[1260];

    auto g_y_0_0_0_z_zz_xx_xy = buffer_1000_pddd[1261];

    auto g_y_0_0_0_z_zz_xx_xz = buffer_1000_pddd[1262];

    auto g_y_0_0_0_z_zz_xx_yy = buffer_1000_pddd[1263];

    auto g_y_0_0_0_z_zz_xx_yz = buffer_1000_pddd[1264];

    auto g_y_0_0_0_z_zz_xx_zz = buffer_1000_pddd[1265];

    auto g_y_0_0_0_z_zz_xy_xx = buffer_1000_pddd[1266];

    auto g_y_0_0_0_z_zz_xy_xy = buffer_1000_pddd[1267];

    auto g_y_0_0_0_z_zz_xy_xz = buffer_1000_pddd[1268];

    auto g_y_0_0_0_z_zz_xy_yy = buffer_1000_pddd[1269];

    auto g_y_0_0_0_z_zz_xy_yz = buffer_1000_pddd[1270];

    auto g_y_0_0_0_z_zz_xy_zz = buffer_1000_pddd[1271];

    auto g_y_0_0_0_z_zz_xz_xx = buffer_1000_pddd[1272];

    auto g_y_0_0_0_z_zz_xz_xy = buffer_1000_pddd[1273];

    auto g_y_0_0_0_z_zz_xz_xz = buffer_1000_pddd[1274];

    auto g_y_0_0_0_z_zz_xz_yy = buffer_1000_pddd[1275];

    auto g_y_0_0_0_z_zz_xz_yz = buffer_1000_pddd[1276];

    auto g_y_0_0_0_z_zz_xz_zz = buffer_1000_pddd[1277];

    auto g_y_0_0_0_z_zz_yy_xx = buffer_1000_pddd[1278];

    auto g_y_0_0_0_z_zz_yy_xy = buffer_1000_pddd[1279];

    auto g_y_0_0_0_z_zz_yy_xz = buffer_1000_pddd[1280];

    auto g_y_0_0_0_z_zz_yy_yy = buffer_1000_pddd[1281];

    auto g_y_0_0_0_z_zz_yy_yz = buffer_1000_pddd[1282];

    auto g_y_0_0_0_z_zz_yy_zz = buffer_1000_pddd[1283];

    auto g_y_0_0_0_z_zz_yz_xx = buffer_1000_pddd[1284];

    auto g_y_0_0_0_z_zz_yz_xy = buffer_1000_pddd[1285];

    auto g_y_0_0_0_z_zz_yz_xz = buffer_1000_pddd[1286];

    auto g_y_0_0_0_z_zz_yz_yy = buffer_1000_pddd[1287];

    auto g_y_0_0_0_z_zz_yz_yz = buffer_1000_pddd[1288];

    auto g_y_0_0_0_z_zz_yz_zz = buffer_1000_pddd[1289];

    auto g_y_0_0_0_z_zz_zz_xx = buffer_1000_pddd[1290];

    auto g_y_0_0_0_z_zz_zz_xy = buffer_1000_pddd[1291];

    auto g_y_0_0_0_z_zz_zz_xz = buffer_1000_pddd[1292];

    auto g_y_0_0_0_z_zz_zz_yy = buffer_1000_pddd[1293];

    auto g_y_0_0_0_z_zz_zz_yz = buffer_1000_pddd[1294];

    auto g_y_0_0_0_z_zz_zz_zz = buffer_1000_pddd[1295];

    auto g_z_0_0_0_x_xx_xx_xx = buffer_1000_pddd[1296];

    auto g_z_0_0_0_x_xx_xx_xy = buffer_1000_pddd[1297];

    auto g_z_0_0_0_x_xx_xx_xz = buffer_1000_pddd[1298];

    auto g_z_0_0_0_x_xx_xx_yy = buffer_1000_pddd[1299];

    auto g_z_0_0_0_x_xx_xx_yz = buffer_1000_pddd[1300];

    auto g_z_0_0_0_x_xx_xx_zz = buffer_1000_pddd[1301];

    auto g_z_0_0_0_x_xx_xy_xx = buffer_1000_pddd[1302];

    auto g_z_0_0_0_x_xx_xy_xy = buffer_1000_pddd[1303];

    auto g_z_0_0_0_x_xx_xy_xz = buffer_1000_pddd[1304];

    auto g_z_0_0_0_x_xx_xy_yy = buffer_1000_pddd[1305];

    auto g_z_0_0_0_x_xx_xy_yz = buffer_1000_pddd[1306];

    auto g_z_0_0_0_x_xx_xy_zz = buffer_1000_pddd[1307];

    auto g_z_0_0_0_x_xx_xz_xx = buffer_1000_pddd[1308];

    auto g_z_0_0_0_x_xx_xz_xy = buffer_1000_pddd[1309];

    auto g_z_0_0_0_x_xx_xz_xz = buffer_1000_pddd[1310];

    auto g_z_0_0_0_x_xx_xz_yy = buffer_1000_pddd[1311];

    auto g_z_0_0_0_x_xx_xz_yz = buffer_1000_pddd[1312];

    auto g_z_0_0_0_x_xx_xz_zz = buffer_1000_pddd[1313];

    auto g_z_0_0_0_x_xx_yy_xx = buffer_1000_pddd[1314];

    auto g_z_0_0_0_x_xx_yy_xy = buffer_1000_pddd[1315];

    auto g_z_0_0_0_x_xx_yy_xz = buffer_1000_pddd[1316];

    auto g_z_0_0_0_x_xx_yy_yy = buffer_1000_pddd[1317];

    auto g_z_0_0_0_x_xx_yy_yz = buffer_1000_pddd[1318];

    auto g_z_0_0_0_x_xx_yy_zz = buffer_1000_pddd[1319];

    auto g_z_0_0_0_x_xx_yz_xx = buffer_1000_pddd[1320];

    auto g_z_0_0_0_x_xx_yz_xy = buffer_1000_pddd[1321];

    auto g_z_0_0_0_x_xx_yz_xz = buffer_1000_pddd[1322];

    auto g_z_0_0_0_x_xx_yz_yy = buffer_1000_pddd[1323];

    auto g_z_0_0_0_x_xx_yz_yz = buffer_1000_pddd[1324];

    auto g_z_0_0_0_x_xx_yz_zz = buffer_1000_pddd[1325];

    auto g_z_0_0_0_x_xx_zz_xx = buffer_1000_pddd[1326];

    auto g_z_0_0_0_x_xx_zz_xy = buffer_1000_pddd[1327];

    auto g_z_0_0_0_x_xx_zz_xz = buffer_1000_pddd[1328];

    auto g_z_0_0_0_x_xx_zz_yy = buffer_1000_pddd[1329];

    auto g_z_0_0_0_x_xx_zz_yz = buffer_1000_pddd[1330];

    auto g_z_0_0_0_x_xx_zz_zz = buffer_1000_pddd[1331];

    auto g_z_0_0_0_x_xy_xx_xx = buffer_1000_pddd[1332];

    auto g_z_0_0_0_x_xy_xx_xy = buffer_1000_pddd[1333];

    auto g_z_0_0_0_x_xy_xx_xz = buffer_1000_pddd[1334];

    auto g_z_0_0_0_x_xy_xx_yy = buffer_1000_pddd[1335];

    auto g_z_0_0_0_x_xy_xx_yz = buffer_1000_pddd[1336];

    auto g_z_0_0_0_x_xy_xx_zz = buffer_1000_pddd[1337];

    auto g_z_0_0_0_x_xy_xy_xx = buffer_1000_pddd[1338];

    auto g_z_0_0_0_x_xy_xy_xy = buffer_1000_pddd[1339];

    auto g_z_0_0_0_x_xy_xy_xz = buffer_1000_pddd[1340];

    auto g_z_0_0_0_x_xy_xy_yy = buffer_1000_pddd[1341];

    auto g_z_0_0_0_x_xy_xy_yz = buffer_1000_pddd[1342];

    auto g_z_0_0_0_x_xy_xy_zz = buffer_1000_pddd[1343];

    auto g_z_0_0_0_x_xy_xz_xx = buffer_1000_pddd[1344];

    auto g_z_0_0_0_x_xy_xz_xy = buffer_1000_pddd[1345];

    auto g_z_0_0_0_x_xy_xz_xz = buffer_1000_pddd[1346];

    auto g_z_0_0_0_x_xy_xz_yy = buffer_1000_pddd[1347];

    auto g_z_0_0_0_x_xy_xz_yz = buffer_1000_pddd[1348];

    auto g_z_0_0_0_x_xy_xz_zz = buffer_1000_pddd[1349];

    auto g_z_0_0_0_x_xy_yy_xx = buffer_1000_pddd[1350];

    auto g_z_0_0_0_x_xy_yy_xy = buffer_1000_pddd[1351];

    auto g_z_0_0_0_x_xy_yy_xz = buffer_1000_pddd[1352];

    auto g_z_0_0_0_x_xy_yy_yy = buffer_1000_pddd[1353];

    auto g_z_0_0_0_x_xy_yy_yz = buffer_1000_pddd[1354];

    auto g_z_0_0_0_x_xy_yy_zz = buffer_1000_pddd[1355];

    auto g_z_0_0_0_x_xy_yz_xx = buffer_1000_pddd[1356];

    auto g_z_0_0_0_x_xy_yz_xy = buffer_1000_pddd[1357];

    auto g_z_0_0_0_x_xy_yz_xz = buffer_1000_pddd[1358];

    auto g_z_0_0_0_x_xy_yz_yy = buffer_1000_pddd[1359];

    auto g_z_0_0_0_x_xy_yz_yz = buffer_1000_pddd[1360];

    auto g_z_0_0_0_x_xy_yz_zz = buffer_1000_pddd[1361];

    auto g_z_0_0_0_x_xy_zz_xx = buffer_1000_pddd[1362];

    auto g_z_0_0_0_x_xy_zz_xy = buffer_1000_pddd[1363];

    auto g_z_0_0_0_x_xy_zz_xz = buffer_1000_pddd[1364];

    auto g_z_0_0_0_x_xy_zz_yy = buffer_1000_pddd[1365];

    auto g_z_0_0_0_x_xy_zz_yz = buffer_1000_pddd[1366];

    auto g_z_0_0_0_x_xy_zz_zz = buffer_1000_pddd[1367];

    auto g_z_0_0_0_x_xz_xx_xx = buffer_1000_pddd[1368];

    auto g_z_0_0_0_x_xz_xx_xy = buffer_1000_pddd[1369];

    auto g_z_0_0_0_x_xz_xx_xz = buffer_1000_pddd[1370];

    auto g_z_0_0_0_x_xz_xx_yy = buffer_1000_pddd[1371];

    auto g_z_0_0_0_x_xz_xx_yz = buffer_1000_pddd[1372];

    auto g_z_0_0_0_x_xz_xx_zz = buffer_1000_pddd[1373];

    auto g_z_0_0_0_x_xz_xy_xx = buffer_1000_pddd[1374];

    auto g_z_0_0_0_x_xz_xy_xy = buffer_1000_pddd[1375];

    auto g_z_0_0_0_x_xz_xy_xz = buffer_1000_pddd[1376];

    auto g_z_0_0_0_x_xz_xy_yy = buffer_1000_pddd[1377];

    auto g_z_0_0_0_x_xz_xy_yz = buffer_1000_pddd[1378];

    auto g_z_0_0_0_x_xz_xy_zz = buffer_1000_pddd[1379];

    auto g_z_0_0_0_x_xz_xz_xx = buffer_1000_pddd[1380];

    auto g_z_0_0_0_x_xz_xz_xy = buffer_1000_pddd[1381];

    auto g_z_0_0_0_x_xz_xz_xz = buffer_1000_pddd[1382];

    auto g_z_0_0_0_x_xz_xz_yy = buffer_1000_pddd[1383];

    auto g_z_0_0_0_x_xz_xz_yz = buffer_1000_pddd[1384];

    auto g_z_0_0_0_x_xz_xz_zz = buffer_1000_pddd[1385];

    auto g_z_0_0_0_x_xz_yy_xx = buffer_1000_pddd[1386];

    auto g_z_0_0_0_x_xz_yy_xy = buffer_1000_pddd[1387];

    auto g_z_0_0_0_x_xz_yy_xz = buffer_1000_pddd[1388];

    auto g_z_0_0_0_x_xz_yy_yy = buffer_1000_pddd[1389];

    auto g_z_0_0_0_x_xz_yy_yz = buffer_1000_pddd[1390];

    auto g_z_0_0_0_x_xz_yy_zz = buffer_1000_pddd[1391];

    auto g_z_0_0_0_x_xz_yz_xx = buffer_1000_pddd[1392];

    auto g_z_0_0_0_x_xz_yz_xy = buffer_1000_pddd[1393];

    auto g_z_0_0_0_x_xz_yz_xz = buffer_1000_pddd[1394];

    auto g_z_0_0_0_x_xz_yz_yy = buffer_1000_pddd[1395];

    auto g_z_0_0_0_x_xz_yz_yz = buffer_1000_pddd[1396];

    auto g_z_0_0_0_x_xz_yz_zz = buffer_1000_pddd[1397];

    auto g_z_0_0_0_x_xz_zz_xx = buffer_1000_pddd[1398];

    auto g_z_0_0_0_x_xz_zz_xy = buffer_1000_pddd[1399];

    auto g_z_0_0_0_x_xz_zz_xz = buffer_1000_pddd[1400];

    auto g_z_0_0_0_x_xz_zz_yy = buffer_1000_pddd[1401];

    auto g_z_0_0_0_x_xz_zz_yz = buffer_1000_pddd[1402];

    auto g_z_0_0_0_x_xz_zz_zz = buffer_1000_pddd[1403];

    auto g_z_0_0_0_x_yy_xx_xx = buffer_1000_pddd[1404];

    auto g_z_0_0_0_x_yy_xx_xy = buffer_1000_pddd[1405];

    auto g_z_0_0_0_x_yy_xx_xz = buffer_1000_pddd[1406];

    auto g_z_0_0_0_x_yy_xx_yy = buffer_1000_pddd[1407];

    auto g_z_0_0_0_x_yy_xx_yz = buffer_1000_pddd[1408];

    auto g_z_0_0_0_x_yy_xx_zz = buffer_1000_pddd[1409];

    auto g_z_0_0_0_x_yy_xy_xx = buffer_1000_pddd[1410];

    auto g_z_0_0_0_x_yy_xy_xy = buffer_1000_pddd[1411];

    auto g_z_0_0_0_x_yy_xy_xz = buffer_1000_pddd[1412];

    auto g_z_0_0_0_x_yy_xy_yy = buffer_1000_pddd[1413];

    auto g_z_0_0_0_x_yy_xy_yz = buffer_1000_pddd[1414];

    auto g_z_0_0_0_x_yy_xy_zz = buffer_1000_pddd[1415];

    auto g_z_0_0_0_x_yy_xz_xx = buffer_1000_pddd[1416];

    auto g_z_0_0_0_x_yy_xz_xy = buffer_1000_pddd[1417];

    auto g_z_0_0_0_x_yy_xz_xz = buffer_1000_pddd[1418];

    auto g_z_0_0_0_x_yy_xz_yy = buffer_1000_pddd[1419];

    auto g_z_0_0_0_x_yy_xz_yz = buffer_1000_pddd[1420];

    auto g_z_0_0_0_x_yy_xz_zz = buffer_1000_pddd[1421];

    auto g_z_0_0_0_x_yy_yy_xx = buffer_1000_pddd[1422];

    auto g_z_0_0_0_x_yy_yy_xy = buffer_1000_pddd[1423];

    auto g_z_0_0_0_x_yy_yy_xz = buffer_1000_pddd[1424];

    auto g_z_0_0_0_x_yy_yy_yy = buffer_1000_pddd[1425];

    auto g_z_0_0_0_x_yy_yy_yz = buffer_1000_pddd[1426];

    auto g_z_0_0_0_x_yy_yy_zz = buffer_1000_pddd[1427];

    auto g_z_0_0_0_x_yy_yz_xx = buffer_1000_pddd[1428];

    auto g_z_0_0_0_x_yy_yz_xy = buffer_1000_pddd[1429];

    auto g_z_0_0_0_x_yy_yz_xz = buffer_1000_pddd[1430];

    auto g_z_0_0_0_x_yy_yz_yy = buffer_1000_pddd[1431];

    auto g_z_0_0_0_x_yy_yz_yz = buffer_1000_pddd[1432];

    auto g_z_0_0_0_x_yy_yz_zz = buffer_1000_pddd[1433];

    auto g_z_0_0_0_x_yy_zz_xx = buffer_1000_pddd[1434];

    auto g_z_0_0_0_x_yy_zz_xy = buffer_1000_pddd[1435];

    auto g_z_0_0_0_x_yy_zz_xz = buffer_1000_pddd[1436];

    auto g_z_0_0_0_x_yy_zz_yy = buffer_1000_pddd[1437];

    auto g_z_0_0_0_x_yy_zz_yz = buffer_1000_pddd[1438];

    auto g_z_0_0_0_x_yy_zz_zz = buffer_1000_pddd[1439];

    auto g_z_0_0_0_x_yz_xx_xx = buffer_1000_pddd[1440];

    auto g_z_0_0_0_x_yz_xx_xy = buffer_1000_pddd[1441];

    auto g_z_0_0_0_x_yz_xx_xz = buffer_1000_pddd[1442];

    auto g_z_0_0_0_x_yz_xx_yy = buffer_1000_pddd[1443];

    auto g_z_0_0_0_x_yz_xx_yz = buffer_1000_pddd[1444];

    auto g_z_0_0_0_x_yz_xx_zz = buffer_1000_pddd[1445];

    auto g_z_0_0_0_x_yz_xy_xx = buffer_1000_pddd[1446];

    auto g_z_0_0_0_x_yz_xy_xy = buffer_1000_pddd[1447];

    auto g_z_0_0_0_x_yz_xy_xz = buffer_1000_pddd[1448];

    auto g_z_0_0_0_x_yz_xy_yy = buffer_1000_pddd[1449];

    auto g_z_0_0_0_x_yz_xy_yz = buffer_1000_pddd[1450];

    auto g_z_0_0_0_x_yz_xy_zz = buffer_1000_pddd[1451];

    auto g_z_0_0_0_x_yz_xz_xx = buffer_1000_pddd[1452];

    auto g_z_0_0_0_x_yz_xz_xy = buffer_1000_pddd[1453];

    auto g_z_0_0_0_x_yz_xz_xz = buffer_1000_pddd[1454];

    auto g_z_0_0_0_x_yz_xz_yy = buffer_1000_pddd[1455];

    auto g_z_0_0_0_x_yz_xz_yz = buffer_1000_pddd[1456];

    auto g_z_0_0_0_x_yz_xz_zz = buffer_1000_pddd[1457];

    auto g_z_0_0_0_x_yz_yy_xx = buffer_1000_pddd[1458];

    auto g_z_0_0_0_x_yz_yy_xy = buffer_1000_pddd[1459];

    auto g_z_0_0_0_x_yz_yy_xz = buffer_1000_pddd[1460];

    auto g_z_0_0_0_x_yz_yy_yy = buffer_1000_pddd[1461];

    auto g_z_0_0_0_x_yz_yy_yz = buffer_1000_pddd[1462];

    auto g_z_0_0_0_x_yz_yy_zz = buffer_1000_pddd[1463];

    auto g_z_0_0_0_x_yz_yz_xx = buffer_1000_pddd[1464];

    auto g_z_0_0_0_x_yz_yz_xy = buffer_1000_pddd[1465];

    auto g_z_0_0_0_x_yz_yz_xz = buffer_1000_pddd[1466];

    auto g_z_0_0_0_x_yz_yz_yy = buffer_1000_pddd[1467];

    auto g_z_0_0_0_x_yz_yz_yz = buffer_1000_pddd[1468];

    auto g_z_0_0_0_x_yz_yz_zz = buffer_1000_pddd[1469];

    auto g_z_0_0_0_x_yz_zz_xx = buffer_1000_pddd[1470];

    auto g_z_0_0_0_x_yz_zz_xy = buffer_1000_pddd[1471];

    auto g_z_0_0_0_x_yz_zz_xz = buffer_1000_pddd[1472];

    auto g_z_0_0_0_x_yz_zz_yy = buffer_1000_pddd[1473];

    auto g_z_0_0_0_x_yz_zz_yz = buffer_1000_pddd[1474];

    auto g_z_0_0_0_x_yz_zz_zz = buffer_1000_pddd[1475];

    auto g_z_0_0_0_x_zz_xx_xx = buffer_1000_pddd[1476];

    auto g_z_0_0_0_x_zz_xx_xy = buffer_1000_pddd[1477];

    auto g_z_0_0_0_x_zz_xx_xz = buffer_1000_pddd[1478];

    auto g_z_0_0_0_x_zz_xx_yy = buffer_1000_pddd[1479];

    auto g_z_0_0_0_x_zz_xx_yz = buffer_1000_pddd[1480];

    auto g_z_0_0_0_x_zz_xx_zz = buffer_1000_pddd[1481];

    auto g_z_0_0_0_x_zz_xy_xx = buffer_1000_pddd[1482];

    auto g_z_0_0_0_x_zz_xy_xy = buffer_1000_pddd[1483];

    auto g_z_0_0_0_x_zz_xy_xz = buffer_1000_pddd[1484];

    auto g_z_0_0_0_x_zz_xy_yy = buffer_1000_pddd[1485];

    auto g_z_0_0_0_x_zz_xy_yz = buffer_1000_pddd[1486];

    auto g_z_0_0_0_x_zz_xy_zz = buffer_1000_pddd[1487];

    auto g_z_0_0_0_x_zz_xz_xx = buffer_1000_pddd[1488];

    auto g_z_0_0_0_x_zz_xz_xy = buffer_1000_pddd[1489];

    auto g_z_0_0_0_x_zz_xz_xz = buffer_1000_pddd[1490];

    auto g_z_0_0_0_x_zz_xz_yy = buffer_1000_pddd[1491];

    auto g_z_0_0_0_x_zz_xz_yz = buffer_1000_pddd[1492];

    auto g_z_0_0_0_x_zz_xz_zz = buffer_1000_pddd[1493];

    auto g_z_0_0_0_x_zz_yy_xx = buffer_1000_pddd[1494];

    auto g_z_0_0_0_x_zz_yy_xy = buffer_1000_pddd[1495];

    auto g_z_0_0_0_x_zz_yy_xz = buffer_1000_pddd[1496];

    auto g_z_0_0_0_x_zz_yy_yy = buffer_1000_pddd[1497];

    auto g_z_0_0_0_x_zz_yy_yz = buffer_1000_pddd[1498];

    auto g_z_0_0_0_x_zz_yy_zz = buffer_1000_pddd[1499];

    auto g_z_0_0_0_x_zz_yz_xx = buffer_1000_pddd[1500];

    auto g_z_0_0_0_x_zz_yz_xy = buffer_1000_pddd[1501];

    auto g_z_0_0_0_x_zz_yz_xz = buffer_1000_pddd[1502];

    auto g_z_0_0_0_x_zz_yz_yy = buffer_1000_pddd[1503];

    auto g_z_0_0_0_x_zz_yz_yz = buffer_1000_pddd[1504];

    auto g_z_0_0_0_x_zz_yz_zz = buffer_1000_pddd[1505];

    auto g_z_0_0_0_x_zz_zz_xx = buffer_1000_pddd[1506];

    auto g_z_0_0_0_x_zz_zz_xy = buffer_1000_pddd[1507];

    auto g_z_0_0_0_x_zz_zz_xz = buffer_1000_pddd[1508];

    auto g_z_0_0_0_x_zz_zz_yy = buffer_1000_pddd[1509];

    auto g_z_0_0_0_x_zz_zz_yz = buffer_1000_pddd[1510];

    auto g_z_0_0_0_x_zz_zz_zz = buffer_1000_pddd[1511];

    auto g_z_0_0_0_y_xx_xx_xx = buffer_1000_pddd[1512];

    auto g_z_0_0_0_y_xx_xx_xy = buffer_1000_pddd[1513];

    auto g_z_0_0_0_y_xx_xx_xz = buffer_1000_pddd[1514];

    auto g_z_0_0_0_y_xx_xx_yy = buffer_1000_pddd[1515];

    auto g_z_0_0_0_y_xx_xx_yz = buffer_1000_pddd[1516];

    auto g_z_0_0_0_y_xx_xx_zz = buffer_1000_pddd[1517];

    auto g_z_0_0_0_y_xx_xy_xx = buffer_1000_pddd[1518];

    auto g_z_0_0_0_y_xx_xy_xy = buffer_1000_pddd[1519];

    auto g_z_0_0_0_y_xx_xy_xz = buffer_1000_pddd[1520];

    auto g_z_0_0_0_y_xx_xy_yy = buffer_1000_pddd[1521];

    auto g_z_0_0_0_y_xx_xy_yz = buffer_1000_pddd[1522];

    auto g_z_0_0_0_y_xx_xy_zz = buffer_1000_pddd[1523];

    auto g_z_0_0_0_y_xx_xz_xx = buffer_1000_pddd[1524];

    auto g_z_0_0_0_y_xx_xz_xy = buffer_1000_pddd[1525];

    auto g_z_0_0_0_y_xx_xz_xz = buffer_1000_pddd[1526];

    auto g_z_0_0_0_y_xx_xz_yy = buffer_1000_pddd[1527];

    auto g_z_0_0_0_y_xx_xz_yz = buffer_1000_pddd[1528];

    auto g_z_0_0_0_y_xx_xz_zz = buffer_1000_pddd[1529];

    auto g_z_0_0_0_y_xx_yy_xx = buffer_1000_pddd[1530];

    auto g_z_0_0_0_y_xx_yy_xy = buffer_1000_pddd[1531];

    auto g_z_0_0_0_y_xx_yy_xz = buffer_1000_pddd[1532];

    auto g_z_0_0_0_y_xx_yy_yy = buffer_1000_pddd[1533];

    auto g_z_0_0_0_y_xx_yy_yz = buffer_1000_pddd[1534];

    auto g_z_0_0_0_y_xx_yy_zz = buffer_1000_pddd[1535];

    auto g_z_0_0_0_y_xx_yz_xx = buffer_1000_pddd[1536];

    auto g_z_0_0_0_y_xx_yz_xy = buffer_1000_pddd[1537];

    auto g_z_0_0_0_y_xx_yz_xz = buffer_1000_pddd[1538];

    auto g_z_0_0_0_y_xx_yz_yy = buffer_1000_pddd[1539];

    auto g_z_0_0_0_y_xx_yz_yz = buffer_1000_pddd[1540];

    auto g_z_0_0_0_y_xx_yz_zz = buffer_1000_pddd[1541];

    auto g_z_0_0_0_y_xx_zz_xx = buffer_1000_pddd[1542];

    auto g_z_0_0_0_y_xx_zz_xy = buffer_1000_pddd[1543];

    auto g_z_0_0_0_y_xx_zz_xz = buffer_1000_pddd[1544];

    auto g_z_0_0_0_y_xx_zz_yy = buffer_1000_pddd[1545];

    auto g_z_0_0_0_y_xx_zz_yz = buffer_1000_pddd[1546];

    auto g_z_0_0_0_y_xx_zz_zz = buffer_1000_pddd[1547];

    auto g_z_0_0_0_y_xy_xx_xx = buffer_1000_pddd[1548];

    auto g_z_0_0_0_y_xy_xx_xy = buffer_1000_pddd[1549];

    auto g_z_0_0_0_y_xy_xx_xz = buffer_1000_pddd[1550];

    auto g_z_0_0_0_y_xy_xx_yy = buffer_1000_pddd[1551];

    auto g_z_0_0_0_y_xy_xx_yz = buffer_1000_pddd[1552];

    auto g_z_0_0_0_y_xy_xx_zz = buffer_1000_pddd[1553];

    auto g_z_0_0_0_y_xy_xy_xx = buffer_1000_pddd[1554];

    auto g_z_0_0_0_y_xy_xy_xy = buffer_1000_pddd[1555];

    auto g_z_0_0_0_y_xy_xy_xz = buffer_1000_pddd[1556];

    auto g_z_0_0_0_y_xy_xy_yy = buffer_1000_pddd[1557];

    auto g_z_0_0_0_y_xy_xy_yz = buffer_1000_pddd[1558];

    auto g_z_0_0_0_y_xy_xy_zz = buffer_1000_pddd[1559];

    auto g_z_0_0_0_y_xy_xz_xx = buffer_1000_pddd[1560];

    auto g_z_0_0_0_y_xy_xz_xy = buffer_1000_pddd[1561];

    auto g_z_0_0_0_y_xy_xz_xz = buffer_1000_pddd[1562];

    auto g_z_0_0_0_y_xy_xz_yy = buffer_1000_pddd[1563];

    auto g_z_0_0_0_y_xy_xz_yz = buffer_1000_pddd[1564];

    auto g_z_0_0_0_y_xy_xz_zz = buffer_1000_pddd[1565];

    auto g_z_0_0_0_y_xy_yy_xx = buffer_1000_pddd[1566];

    auto g_z_0_0_0_y_xy_yy_xy = buffer_1000_pddd[1567];

    auto g_z_0_0_0_y_xy_yy_xz = buffer_1000_pddd[1568];

    auto g_z_0_0_0_y_xy_yy_yy = buffer_1000_pddd[1569];

    auto g_z_0_0_0_y_xy_yy_yz = buffer_1000_pddd[1570];

    auto g_z_0_0_0_y_xy_yy_zz = buffer_1000_pddd[1571];

    auto g_z_0_0_0_y_xy_yz_xx = buffer_1000_pddd[1572];

    auto g_z_0_0_0_y_xy_yz_xy = buffer_1000_pddd[1573];

    auto g_z_0_0_0_y_xy_yz_xz = buffer_1000_pddd[1574];

    auto g_z_0_0_0_y_xy_yz_yy = buffer_1000_pddd[1575];

    auto g_z_0_0_0_y_xy_yz_yz = buffer_1000_pddd[1576];

    auto g_z_0_0_0_y_xy_yz_zz = buffer_1000_pddd[1577];

    auto g_z_0_0_0_y_xy_zz_xx = buffer_1000_pddd[1578];

    auto g_z_0_0_0_y_xy_zz_xy = buffer_1000_pddd[1579];

    auto g_z_0_0_0_y_xy_zz_xz = buffer_1000_pddd[1580];

    auto g_z_0_0_0_y_xy_zz_yy = buffer_1000_pddd[1581];

    auto g_z_0_0_0_y_xy_zz_yz = buffer_1000_pddd[1582];

    auto g_z_0_0_0_y_xy_zz_zz = buffer_1000_pddd[1583];

    auto g_z_0_0_0_y_xz_xx_xx = buffer_1000_pddd[1584];

    auto g_z_0_0_0_y_xz_xx_xy = buffer_1000_pddd[1585];

    auto g_z_0_0_0_y_xz_xx_xz = buffer_1000_pddd[1586];

    auto g_z_0_0_0_y_xz_xx_yy = buffer_1000_pddd[1587];

    auto g_z_0_0_0_y_xz_xx_yz = buffer_1000_pddd[1588];

    auto g_z_0_0_0_y_xz_xx_zz = buffer_1000_pddd[1589];

    auto g_z_0_0_0_y_xz_xy_xx = buffer_1000_pddd[1590];

    auto g_z_0_0_0_y_xz_xy_xy = buffer_1000_pddd[1591];

    auto g_z_0_0_0_y_xz_xy_xz = buffer_1000_pddd[1592];

    auto g_z_0_0_0_y_xz_xy_yy = buffer_1000_pddd[1593];

    auto g_z_0_0_0_y_xz_xy_yz = buffer_1000_pddd[1594];

    auto g_z_0_0_0_y_xz_xy_zz = buffer_1000_pddd[1595];

    auto g_z_0_0_0_y_xz_xz_xx = buffer_1000_pddd[1596];

    auto g_z_0_0_0_y_xz_xz_xy = buffer_1000_pddd[1597];

    auto g_z_0_0_0_y_xz_xz_xz = buffer_1000_pddd[1598];

    auto g_z_0_0_0_y_xz_xz_yy = buffer_1000_pddd[1599];

    auto g_z_0_0_0_y_xz_xz_yz = buffer_1000_pddd[1600];

    auto g_z_0_0_0_y_xz_xz_zz = buffer_1000_pddd[1601];

    auto g_z_0_0_0_y_xz_yy_xx = buffer_1000_pddd[1602];

    auto g_z_0_0_0_y_xz_yy_xy = buffer_1000_pddd[1603];

    auto g_z_0_0_0_y_xz_yy_xz = buffer_1000_pddd[1604];

    auto g_z_0_0_0_y_xz_yy_yy = buffer_1000_pddd[1605];

    auto g_z_0_0_0_y_xz_yy_yz = buffer_1000_pddd[1606];

    auto g_z_0_0_0_y_xz_yy_zz = buffer_1000_pddd[1607];

    auto g_z_0_0_0_y_xz_yz_xx = buffer_1000_pddd[1608];

    auto g_z_0_0_0_y_xz_yz_xy = buffer_1000_pddd[1609];

    auto g_z_0_0_0_y_xz_yz_xz = buffer_1000_pddd[1610];

    auto g_z_0_0_0_y_xz_yz_yy = buffer_1000_pddd[1611];

    auto g_z_0_0_0_y_xz_yz_yz = buffer_1000_pddd[1612];

    auto g_z_0_0_0_y_xz_yz_zz = buffer_1000_pddd[1613];

    auto g_z_0_0_0_y_xz_zz_xx = buffer_1000_pddd[1614];

    auto g_z_0_0_0_y_xz_zz_xy = buffer_1000_pddd[1615];

    auto g_z_0_0_0_y_xz_zz_xz = buffer_1000_pddd[1616];

    auto g_z_0_0_0_y_xz_zz_yy = buffer_1000_pddd[1617];

    auto g_z_0_0_0_y_xz_zz_yz = buffer_1000_pddd[1618];

    auto g_z_0_0_0_y_xz_zz_zz = buffer_1000_pddd[1619];

    auto g_z_0_0_0_y_yy_xx_xx = buffer_1000_pddd[1620];

    auto g_z_0_0_0_y_yy_xx_xy = buffer_1000_pddd[1621];

    auto g_z_0_0_0_y_yy_xx_xz = buffer_1000_pddd[1622];

    auto g_z_0_0_0_y_yy_xx_yy = buffer_1000_pddd[1623];

    auto g_z_0_0_0_y_yy_xx_yz = buffer_1000_pddd[1624];

    auto g_z_0_0_0_y_yy_xx_zz = buffer_1000_pddd[1625];

    auto g_z_0_0_0_y_yy_xy_xx = buffer_1000_pddd[1626];

    auto g_z_0_0_0_y_yy_xy_xy = buffer_1000_pddd[1627];

    auto g_z_0_0_0_y_yy_xy_xz = buffer_1000_pddd[1628];

    auto g_z_0_0_0_y_yy_xy_yy = buffer_1000_pddd[1629];

    auto g_z_0_0_0_y_yy_xy_yz = buffer_1000_pddd[1630];

    auto g_z_0_0_0_y_yy_xy_zz = buffer_1000_pddd[1631];

    auto g_z_0_0_0_y_yy_xz_xx = buffer_1000_pddd[1632];

    auto g_z_0_0_0_y_yy_xz_xy = buffer_1000_pddd[1633];

    auto g_z_0_0_0_y_yy_xz_xz = buffer_1000_pddd[1634];

    auto g_z_0_0_0_y_yy_xz_yy = buffer_1000_pddd[1635];

    auto g_z_0_0_0_y_yy_xz_yz = buffer_1000_pddd[1636];

    auto g_z_0_0_0_y_yy_xz_zz = buffer_1000_pddd[1637];

    auto g_z_0_0_0_y_yy_yy_xx = buffer_1000_pddd[1638];

    auto g_z_0_0_0_y_yy_yy_xy = buffer_1000_pddd[1639];

    auto g_z_0_0_0_y_yy_yy_xz = buffer_1000_pddd[1640];

    auto g_z_0_0_0_y_yy_yy_yy = buffer_1000_pddd[1641];

    auto g_z_0_0_0_y_yy_yy_yz = buffer_1000_pddd[1642];

    auto g_z_0_0_0_y_yy_yy_zz = buffer_1000_pddd[1643];

    auto g_z_0_0_0_y_yy_yz_xx = buffer_1000_pddd[1644];

    auto g_z_0_0_0_y_yy_yz_xy = buffer_1000_pddd[1645];

    auto g_z_0_0_0_y_yy_yz_xz = buffer_1000_pddd[1646];

    auto g_z_0_0_0_y_yy_yz_yy = buffer_1000_pddd[1647];

    auto g_z_0_0_0_y_yy_yz_yz = buffer_1000_pddd[1648];

    auto g_z_0_0_0_y_yy_yz_zz = buffer_1000_pddd[1649];

    auto g_z_0_0_0_y_yy_zz_xx = buffer_1000_pddd[1650];

    auto g_z_0_0_0_y_yy_zz_xy = buffer_1000_pddd[1651];

    auto g_z_0_0_0_y_yy_zz_xz = buffer_1000_pddd[1652];

    auto g_z_0_0_0_y_yy_zz_yy = buffer_1000_pddd[1653];

    auto g_z_0_0_0_y_yy_zz_yz = buffer_1000_pddd[1654];

    auto g_z_0_0_0_y_yy_zz_zz = buffer_1000_pddd[1655];

    auto g_z_0_0_0_y_yz_xx_xx = buffer_1000_pddd[1656];

    auto g_z_0_0_0_y_yz_xx_xy = buffer_1000_pddd[1657];

    auto g_z_0_0_0_y_yz_xx_xz = buffer_1000_pddd[1658];

    auto g_z_0_0_0_y_yz_xx_yy = buffer_1000_pddd[1659];

    auto g_z_0_0_0_y_yz_xx_yz = buffer_1000_pddd[1660];

    auto g_z_0_0_0_y_yz_xx_zz = buffer_1000_pddd[1661];

    auto g_z_0_0_0_y_yz_xy_xx = buffer_1000_pddd[1662];

    auto g_z_0_0_0_y_yz_xy_xy = buffer_1000_pddd[1663];

    auto g_z_0_0_0_y_yz_xy_xz = buffer_1000_pddd[1664];

    auto g_z_0_0_0_y_yz_xy_yy = buffer_1000_pddd[1665];

    auto g_z_0_0_0_y_yz_xy_yz = buffer_1000_pddd[1666];

    auto g_z_0_0_0_y_yz_xy_zz = buffer_1000_pddd[1667];

    auto g_z_0_0_0_y_yz_xz_xx = buffer_1000_pddd[1668];

    auto g_z_0_0_0_y_yz_xz_xy = buffer_1000_pddd[1669];

    auto g_z_0_0_0_y_yz_xz_xz = buffer_1000_pddd[1670];

    auto g_z_0_0_0_y_yz_xz_yy = buffer_1000_pddd[1671];

    auto g_z_0_0_0_y_yz_xz_yz = buffer_1000_pddd[1672];

    auto g_z_0_0_0_y_yz_xz_zz = buffer_1000_pddd[1673];

    auto g_z_0_0_0_y_yz_yy_xx = buffer_1000_pddd[1674];

    auto g_z_0_0_0_y_yz_yy_xy = buffer_1000_pddd[1675];

    auto g_z_0_0_0_y_yz_yy_xz = buffer_1000_pddd[1676];

    auto g_z_0_0_0_y_yz_yy_yy = buffer_1000_pddd[1677];

    auto g_z_0_0_0_y_yz_yy_yz = buffer_1000_pddd[1678];

    auto g_z_0_0_0_y_yz_yy_zz = buffer_1000_pddd[1679];

    auto g_z_0_0_0_y_yz_yz_xx = buffer_1000_pddd[1680];

    auto g_z_0_0_0_y_yz_yz_xy = buffer_1000_pddd[1681];

    auto g_z_0_0_0_y_yz_yz_xz = buffer_1000_pddd[1682];

    auto g_z_0_0_0_y_yz_yz_yy = buffer_1000_pddd[1683];

    auto g_z_0_0_0_y_yz_yz_yz = buffer_1000_pddd[1684];

    auto g_z_0_0_0_y_yz_yz_zz = buffer_1000_pddd[1685];

    auto g_z_0_0_0_y_yz_zz_xx = buffer_1000_pddd[1686];

    auto g_z_0_0_0_y_yz_zz_xy = buffer_1000_pddd[1687];

    auto g_z_0_0_0_y_yz_zz_xz = buffer_1000_pddd[1688];

    auto g_z_0_0_0_y_yz_zz_yy = buffer_1000_pddd[1689];

    auto g_z_0_0_0_y_yz_zz_yz = buffer_1000_pddd[1690];

    auto g_z_0_0_0_y_yz_zz_zz = buffer_1000_pddd[1691];

    auto g_z_0_0_0_y_zz_xx_xx = buffer_1000_pddd[1692];

    auto g_z_0_0_0_y_zz_xx_xy = buffer_1000_pddd[1693];

    auto g_z_0_0_0_y_zz_xx_xz = buffer_1000_pddd[1694];

    auto g_z_0_0_0_y_zz_xx_yy = buffer_1000_pddd[1695];

    auto g_z_0_0_0_y_zz_xx_yz = buffer_1000_pddd[1696];

    auto g_z_0_0_0_y_zz_xx_zz = buffer_1000_pddd[1697];

    auto g_z_0_0_0_y_zz_xy_xx = buffer_1000_pddd[1698];

    auto g_z_0_0_0_y_zz_xy_xy = buffer_1000_pddd[1699];

    auto g_z_0_0_0_y_zz_xy_xz = buffer_1000_pddd[1700];

    auto g_z_0_0_0_y_zz_xy_yy = buffer_1000_pddd[1701];

    auto g_z_0_0_0_y_zz_xy_yz = buffer_1000_pddd[1702];

    auto g_z_0_0_0_y_zz_xy_zz = buffer_1000_pddd[1703];

    auto g_z_0_0_0_y_zz_xz_xx = buffer_1000_pddd[1704];

    auto g_z_0_0_0_y_zz_xz_xy = buffer_1000_pddd[1705];

    auto g_z_0_0_0_y_zz_xz_xz = buffer_1000_pddd[1706];

    auto g_z_0_0_0_y_zz_xz_yy = buffer_1000_pddd[1707];

    auto g_z_0_0_0_y_zz_xz_yz = buffer_1000_pddd[1708];

    auto g_z_0_0_0_y_zz_xz_zz = buffer_1000_pddd[1709];

    auto g_z_0_0_0_y_zz_yy_xx = buffer_1000_pddd[1710];

    auto g_z_0_0_0_y_zz_yy_xy = buffer_1000_pddd[1711];

    auto g_z_0_0_0_y_zz_yy_xz = buffer_1000_pddd[1712];

    auto g_z_0_0_0_y_zz_yy_yy = buffer_1000_pddd[1713];

    auto g_z_0_0_0_y_zz_yy_yz = buffer_1000_pddd[1714];

    auto g_z_0_0_0_y_zz_yy_zz = buffer_1000_pddd[1715];

    auto g_z_0_0_0_y_zz_yz_xx = buffer_1000_pddd[1716];

    auto g_z_0_0_0_y_zz_yz_xy = buffer_1000_pddd[1717];

    auto g_z_0_0_0_y_zz_yz_xz = buffer_1000_pddd[1718];

    auto g_z_0_0_0_y_zz_yz_yy = buffer_1000_pddd[1719];

    auto g_z_0_0_0_y_zz_yz_yz = buffer_1000_pddd[1720];

    auto g_z_0_0_0_y_zz_yz_zz = buffer_1000_pddd[1721];

    auto g_z_0_0_0_y_zz_zz_xx = buffer_1000_pddd[1722];

    auto g_z_0_0_0_y_zz_zz_xy = buffer_1000_pddd[1723];

    auto g_z_0_0_0_y_zz_zz_xz = buffer_1000_pddd[1724];

    auto g_z_0_0_0_y_zz_zz_yy = buffer_1000_pddd[1725];

    auto g_z_0_0_0_y_zz_zz_yz = buffer_1000_pddd[1726];

    auto g_z_0_0_0_y_zz_zz_zz = buffer_1000_pddd[1727];

    auto g_z_0_0_0_z_xx_xx_xx = buffer_1000_pddd[1728];

    auto g_z_0_0_0_z_xx_xx_xy = buffer_1000_pddd[1729];

    auto g_z_0_0_0_z_xx_xx_xz = buffer_1000_pddd[1730];

    auto g_z_0_0_0_z_xx_xx_yy = buffer_1000_pddd[1731];

    auto g_z_0_0_0_z_xx_xx_yz = buffer_1000_pddd[1732];

    auto g_z_0_0_0_z_xx_xx_zz = buffer_1000_pddd[1733];

    auto g_z_0_0_0_z_xx_xy_xx = buffer_1000_pddd[1734];

    auto g_z_0_0_0_z_xx_xy_xy = buffer_1000_pddd[1735];

    auto g_z_0_0_0_z_xx_xy_xz = buffer_1000_pddd[1736];

    auto g_z_0_0_0_z_xx_xy_yy = buffer_1000_pddd[1737];

    auto g_z_0_0_0_z_xx_xy_yz = buffer_1000_pddd[1738];

    auto g_z_0_0_0_z_xx_xy_zz = buffer_1000_pddd[1739];

    auto g_z_0_0_0_z_xx_xz_xx = buffer_1000_pddd[1740];

    auto g_z_0_0_0_z_xx_xz_xy = buffer_1000_pddd[1741];

    auto g_z_0_0_0_z_xx_xz_xz = buffer_1000_pddd[1742];

    auto g_z_0_0_0_z_xx_xz_yy = buffer_1000_pddd[1743];

    auto g_z_0_0_0_z_xx_xz_yz = buffer_1000_pddd[1744];

    auto g_z_0_0_0_z_xx_xz_zz = buffer_1000_pddd[1745];

    auto g_z_0_0_0_z_xx_yy_xx = buffer_1000_pddd[1746];

    auto g_z_0_0_0_z_xx_yy_xy = buffer_1000_pddd[1747];

    auto g_z_0_0_0_z_xx_yy_xz = buffer_1000_pddd[1748];

    auto g_z_0_0_0_z_xx_yy_yy = buffer_1000_pddd[1749];

    auto g_z_0_0_0_z_xx_yy_yz = buffer_1000_pddd[1750];

    auto g_z_0_0_0_z_xx_yy_zz = buffer_1000_pddd[1751];

    auto g_z_0_0_0_z_xx_yz_xx = buffer_1000_pddd[1752];

    auto g_z_0_0_0_z_xx_yz_xy = buffer_1000_pddd[1753];

    auto g_z_0_0_0_z_xx_yz_xz = buffer_1000_pddd[1754];

    auto g_z_0_0_0_z_xx_yz_yy = buffer_1000_pddd[1755];

    auto g_z_0_0_0_z_xx_yz_yz = buffer_1000_pddd[1756];

    auto g_z_0_0_0_z_xx_yz_zz = buffer_1000_pddd[1757];

    auto g_z_0_0_0_z_xx_zz_xx = buffer_1000_pddd[1758];

    auto g_z_0_0_0_z_xx_zz_xy = buffer_1000_pddd[1759];

    auto g_z_0_0_0_z_xx_zz_xz = buffer_1000_pddd[1760];

    auto g_z_0_0_0_z_xx_zz_yy = buffer_1000_pddd[1761];

    auto g_z_0_0_0_z_xx_zz_yz = buffer_1000_pddd[1762];

    auto g_z_0_0_0_z_xx_zz_zz = buffer_1000_pddd[1763];

    auto g_z_0_0_0_z_xy_xx_xx = buffer_1000_pddd[1764];

    auto g_z_0_0_0_z_xy_xx_xy = buffer_1000_pddd[1765];

    auto g_z_0_0_0_z_xy_xx_xz = buffer_1000_pddd[1766];

    auto g_z_0_0_0_z_xy_xx_yy = buffer_1000_pddd[1767];

    auto g_z_0_0_0_z_xy_xx_yz = buffer_1000_pddd[1768];

    auto g_z_0_0_0_z_xy_xx_zz = buffer_1000_pddd[1769];

    auto g_z_0_0_0_z_xy_xy_xx = buffer_1000_pddd[1770];

    auto g_z_0_0_0_z_xy_xy_xy = buffer_1000_pddd[1771];

    auto g_z_0_0_0_z_xy_xy_xz = buffer_1000_pddd[1772];

    auto g_z_0_0_0_z_xy_xy_yy = buffer_1000_pddd[1773];

    auto g_z_0_0_0_z_xy_xy_yz = buffer_1000_pddd[1774];

    auto g_z_0_0_0_z_xy_xy_zz = buffer_1000_pddd[1775];

    auto g_z_0_0_0_z_xy_xz_xx = buffer_1000_pddd[1776];

    auto g_z_0_0_0_z_xy_xz_xy = buffer_1000_pddd[1777];

    auto g_z_0_0_0_z_xy_xz_xz = buffer_1000_pddd[1778];

    auto g_z_0_0_0_z_xy_xz_yy = buffer_1000_pddd[1779];

    auto g_z_0_0_0_z_xy_xz_yz = buffer_1000_pddd[1780];

    auto g_z_0_0_0_z_xy_xz_zz = buffer_1000_pddd[1781];

    auto g_z_0_0_0_z_xy_yy_xx = buffer_1000_pddd[1782];

    auto g_z_0_0_0_z_xy_yy_xy = buffer_1000_pddd[1783];

    auto g_z_0_0_0_z_xy_yy_xz = buffer_1000_pddd[1784];

    auto g_z_0_0_0_z_xy_yy_yy = buffer_1000_pddd[1785];

    auto g_z_0_0_0_z_xy_yy_yz = buffer_1000_pddd[1786];

    auto g_z_0_0_0_z_xy_yy_zz = buffer_1000_pddd[1787];

    auto g_z_0_0_0_z_xy_yz_xx = buffer_1000_pddd[1788];

    auto g_z_0_0_0_z_xy_yz_xy = buffer_1000_pddd[1789];

    auto g_z_0_0_0_z_xy_yz_xz = buffer_1000_pddd[1790];

    auto g_z_0_0_0_z_xy_yz_yy = buffer_1000_pddd[1791];

    auto g_z_0_0_0_z_xy_yz_yz = buffer_1000_pddd[1792];

    auto g_z_0_0_0_z_xy_yz_zz = buffer_1000_pddd[1793];

    auto g_z_0_0_0_z_xy_zz_xx = buffer_1000_pddd[1794];

    auto g_z_0_0_0_z_xy_zz_xy = buffer_1000_pddd[1795];

    auto g_z_0_0_0_z_xy_zz_xz = buffer_1000_pddd[1796];

    auto g_z_0_0_0_z_xy_zz_yy = buffer_1000_pddd[1797];

    auto g_z_0_0_0_z_xy_zz_yz = buffer_1000_pddd[1798];

    auto g_z_0_0_0_z_xy_zz_zz = buffer_1000_pddd[1799];

    auto g_z_0_0_0_z_xz_xx_xx = buffer_1000_pddd[1800];

    auto g_z_0_0_0_z_xz_xx_xy = buffer_1000_pddd[1801];

    auto g_z_0_0_0_z_xz_xx_xz = buffer_1000_pddd[1802];

    auto g_z_0_0_0_z_xz_xx_yy = buffer_1000_pddd[1803];

    auto g_z_0_0_0_z_xz_xx_yz = buffer_1000_pddd[1804];

    auto g_z_0_0_0_z_xz_xx_zz = buffer_1000_pddd[1805];

    auto g_z_0_0_0_z_xz_xy_xx = buffer_1000_pddd[1806];

    auto g_z_0_0_0_z_xz_xy_xy = buffer_1000_pddd[1807];

    auto g_z_0_0_0_z_xz_xy_xz = buffer_1000_pddd[1808];

    auto g_z_0_0_0_z_xz_xy_yy = buffer_1000_pddd[1809];

    auto g_z_0_0_0_z_xz_xy_yz = buffer_1000_pddd[1810];

    auto g_z_0_0_0_z_xz_xy_zz = buffer_1000_pddd[1811];

    auto g_z_0_0_0_z_xz_xz_xx = buffer_1000_pddd[1812];

    auto g_z_0_0_0_z_xz_xz_xy = buffer_1000_pddd[1813];

    auto g_z_0_0_0_z_xz_xz_xz = buffer_1000_pddd[1814];

    auto g_z_0_0_0_z_xz_xz_yy = buffer_1000_pddd[1815];

    auto g_z_0_0_0_z_xz_xz_yz = buffer_1000_pddd[1816];

    auto g_z_0_0_0_z_xz_xz_zz = buffer_1000_pddd[1817];

    auto g_z_0_0_0_z_xz_yy_xx = buffer_1000_pddd[1818];

    auto g_z_0_0_0_z_xz_yy_xy = buffer_1000_pddd[1819];

    auto g_z_0_0_0_z_xz_yy_xz = buffer_1000_pddd[1820];

    auto g_z_0_0_0_z_xz_yy_yy = buffer_1000_pddd[1821];

    auto g_z_0_0_0_z_xz_yy_yz = buffer_1000_pddd[1822];

    auto g_z_0_0_0_z_xz_yy_zz = buffer_1000_pddd[1823];

    auto g_z_0_0_0_z_xz_yz_xx = buffer_1000_pddd[1824];

    auto g_z_0_0_0_z_xz_yz_xy = buffer_1000_pddd[1825];

    auto g_z_0_0_0_z_xz_yz_xz = buffer_1000_pddd[1826];

    auto g_z_0_0_0_z_xz_yz_yy = buffer_1000_pddd[1827];

    auto g_z_0_0_0_z_xz_yz_yz = buffer_1000_pddd[1828];

    auto g_z_0_0_0_z_xz_yz_zz = buffer_1000_pddd[1829];

    auto g_z_0_0_0_z_xz_zz_xx = buffer_1000_pddd[1830];

    auto g_z_0_0_0_z_xz_zz_xy = buffer_1000_pddd[1831];

    auto g_z_0_0_0_z_xz_zz_xz = buffer_1000_pddd[1832];

    auto g_z_0_0_0_z_xz_zz_yy = buffer_1000_pddd[1833];

    auto g_z_0_0_0_z_xz_zz_yz = buffer_1000_pddd[1834];

    auto g_z_0_0_0_z_xz_zz_zz = buffer_1000_pddd[1835];

    auto g_z_0_0_0_z_yy_xx_xx = buffer_1000_pddd[1836];

    auto g_z_0_0_0_z_yy_xx_xy = buffer_1000_pddd[1837];

    auto g_z_0_0_0_z_yy_xx_xz = buffer_1000_pddd[1838];

    auto g_z_0_0_0_z_yy_xx_yy = buffer_1000_pddd[1839];

    auto g_z_0_0_0_z_yy_xx_yz = buffer_1000_pddd[1840];

    auto g_z_0_0_0_z_yy_xx_zz = buffer_1000_pddd[1841];

    auto g_z_0_0_0_z_yy_xy_xx = buffer_1000_pddd[1842];

    auto g_z_0_0_0_z_yy_xy_xy = buffer_1000_pddd[1843];

    auto g_z_0_0_0_z_yy_xy_xz = buffer_1000_pddd[1844];

    auto g_z_0_0_0_z_yy_xy_yy = buffer_1000_pddd[1845];

    auto g_z_0_0_0_z_yy_xy_yz = buffer_1000_pddd[1846];

    auto g_z_0_0_0_z_yy_xy_zz = buffer_1000_pddd[1847];

    auto g_z_0_0_0_z_yy_xz_xx = buffer_1000_pddd[1848];

    auto g_z_0_0_0_z_yy_xz_xy = buffer_1000_pddd[1849];

    auto g_z_0_0_0_z_yy_xz_xz = buffer_1000_pddd[1850];

    auto g_z_0_0_0_z_yy_xz_yy = buffer_1000_pddd[1851];

    auto g_z_0_0_0_z_yy_xz_yz = buffer_1000_pddd[1852];

    auto g_z_0_0_0_z_yy_xz_zz = buffer_1000_pddd[1853];

    auto g_z_0_0_0_z_yy_yy_xx = buffer_1000_pddd[1854];

    auto g_z_0_0_0_z_yy_yy_xy = buffer_1000_pddd[1855];

    auto g_z_0_0_0_z_yy_yy_xz = buffer_1000_pddd[1856];

    auto g_z_0_0_0_z_yy_yy_yy = buffer_1000_pddd[1857];

    auto g_z_0_0_0_z_yy_yy_yz = buffer_1000_pddd[1858];

    auto g_z_0_0_0_z_yy_yy_zz = buffer_1000_pddd[1859];

    auto g_z_0_0_0_z_yy_yz_xx = buffer_1000_pddd[1860];

    auto g_z_0_0_0_z_yy_yz_xy = buffer_1000_pddd[1861];

    auto g_z_0_0_0_z_yy_yz_xz = buffer_1000_pddd[1862];

    auto g_z_0_0_0_z_yy_yz_yy = buffer_1000_pddd[1863];

    auto g_z_0_0_0_z_yy_yz_yz = buffer_1000_pddd[1864];

    auto g_z_0_0_0_z_yy_yz_zz = buffer_1000_pddd[1865];

    auto g_z_0_0_0_z_yy_zz_xx = buffer_1000_pddd[1866];

    auto g_z_0_0_0_z_yy_zz_xy = buffer_1000_pddd[1867];

    auto g_z_0_0_0_z_yy_zz_xz = buffer_1000_pddd[1868];

    auto g_z_0_0_0_z_yy_zz_yy = buffer_1000_pddd[1869];

    auto g_z_0_0_0_z_yy_zz_yz = buffer_1000_pddd[1870];

    auto g_z_0_0_0_z_yy_zz_zz = buffer_1000_pddd[1871];

    auto g_z_0_0_0_z_yz_xx_xx = buffer_1000_pddd[1872];

    auto g_z_0_0_0_z_yz_xx_xy = buffer_1000_pddd[1873];

    auto g_z_0_0_0_z_yz_xx_xz = buffer_1000_pddd[1874];

    auto g_z_0_0_0_z_yz_xx_yy = buffer_1000_pddd[1875];

    auto g_z_0_0_0_z_yz_xx_yz = buffer_1000_pddd[1876];

    auto g_z_0_0_0_z_yz_xx_zz = buffer_1000_pddd[1877];

    auto g_z_0_0_0_z_yz_xy_xx = buffer_1000_pddd[1878];

    auto g_z_0_0_0_z_yz_xy_xy = buffer_1000_pddd[1879];

    auto g_z_0_0_0_z_yz_xy_xz = buffer_1000_pddd[1880];

    auto g_z_0_0_0_z_yz_xy_yy = buffer_1000_pddd[1881];

    auto g_z_0_0_0_z_yz_xy_yz = buffer_1000_pddd[1882];

    auto g_z_0_0_0_z_yz_xy_zz = buffer_1000_pddd[1883];

    auto g_z_0_0_0_z_yz_xz_xx = buffer_1000_pddd[1884];

    auto g_z_0_0_0_z_yz_xz_xy = buffer_1000_pddd[1885];

    auto g_z_0_0_0_z_yz_xz_xz = buffer_1000_pddd[1886];

    auto g_z_0_0_0_z_yz_xz_yy = buffer_1000_pddd[1887];

    auto g_z_0_0_0_z_yz_xz_yz = buffer_1000_pddd[1888];

    auto g_z_0_0_0_z_yz_xz_zz = buffer_1000_pddd[1889];

    auto g_z_0_0_0_z_yz_yy_xx = buffer_1000_pddd[1890];

    auto g_z_0_0_0_z_yz_yy_xy = buffer_1000_pddd[1891];

    auto g_z_0_0_0_z_yz_yy_xz = buffer_1000_pddd[1892];

    auto g_z_0_0_0_z_yz_yy_yy = buffer_1000_pddd[1893];

    auto g_z_0_0_0_z_yz_yy_yz = buffer_1000_pddd[1894];

    auto g_z_0_0_0_z_yz_yy_zz = buffer_1000_pddd[1895];

    auto g_z_0_0_0_z_yz_yz_xx = buffer_1000_pddd[1896];

    auto g_z_0_0_0_z_yz_yz_xy = buffer_1000_pddd[1897];

    auto g_z_0_0_0_z_yz_yz_xz = buffer_1000_pddd[1898];

    auto g_z_0_0_0_z_yz_yz_yy = buffer_1000_pddd[1899];

    auto g_z_0_0_0_z_yz_yz_yz = buffer_1000_pddd[1900];

    auto g_z_0_0_0_z_yz_yz_zz = buffer_1000_pddd[1901];

    auto g_z_0_0_0_z_yz_zz_xx = buffer_1000_pddd[1902];

    auto g_z_0_0_0_z_yz_zz_xy = buffer_1000_pddd[1903];

    auto g_z_0_0_0_z_yz_zz_xz = buffer_1000_pddd[1904];

    auto g_z_0_0_0_z_yz_zz_yy = buffer_1000_pddd[1905];

    auto g_z_0_0_0_z_yz_zz_yz = buffer_1000_pddd[1906];

    auto g_z_0_0_0_z_yz_zz_zz = buffer_1000_pddd[1907];

    auto g_z_0_0_0_z_zz_xx_xx = buffer_1000_pddd[1908];

    auto g_z_0_0_0_z_zz_xx_xy = buffer_1000_pddd[1909];

    auto g_z_0_0_0_z_zz_xx_xz = buffer_1000_pddd[1910];

    auto g_z_0_0_0_z_zz_xx_yy = buffer_1000_pddd[1911];

    auto g_z_0_0_0_z_zz_xx_yz = buffer_1000_pddd[1912];

    auto g_z_0_0_0_z_zz_xx_zz = buffer_1000_pddd[1913];

    auto g_z_0_0_0_z_zz_xy_xx = buffer_1000_pddd[1914];

    auto g_z_0_0_0_z_zz_xy_xy = buffer_1000_pddd[1915];

    auto g_z_0_0_0_z_zz_xy_xz = buffer_1000_pddd[1916];

    auto g_z_0_0_0_z_zz_xy_yy = buffer_1000_pddd[1917];

    auto g_z_0_0_0_z_zz_xy_yz = buffer_1000_pddd[1918];

    auto g_z_0_0_0_z_zz_xy_zz = buffer_1000_pddd[1919];

    auto g_z_0_0_0_z_zz_xz_xx = buffer_1000_pddd[1920];

    auto g_z_0_0_0_z_zz_xz_xy = buffer_1000_pddd[1921];

    auto g_z_0_0_0_z_zz_xz_xz = buffer_1000_pddd[1922];

    auto g_z_0_0_0_z_zz_xz_yy = buffer_1000_pddd[1923];

    auto g_z_0_0_0_z_zz_xz_yz = buffer_1000_pddd[1924];

    auto g_z_0_0_0_z_zz_xz_zz = buffer_1000_pddd[1925];

    auto g_z_0_0_0_z_zz_yy_xx = buffer_1000_pddd[1926];

    auto g_z_0_0_0_z_zz_yy_xy = buffer_1000_pddd[1927];

    auto g_z_0_0_0_z_zz_yy_xz = buffer_1000_pddd[1928];

    auto g_z_0_0_0_z_zz_yy_yy = buffer_1000_pddd[1929];

    auto g_z_0_0_0_z_zz_yy_yz = buffer_1000_pddd[1930];

    auto g_z_0_0_0_z_zz_yy_zz = buffer_1000_pddd[1931];

    auto g_z_0_0_0_z_zz_yz_xx = buffer_1000_pddd[1932];

    auto g_z_0_0_0_z_zz_yz_xy = buffer_1000_pddd[1933];

    auto g_z_0_0_0_z_zz_yz_xz = buffer_1000_pddd[1934];

    auto g_z_0_0_0_z_zz_yz_yy = buffer_1000_pddd[1935];

    auto g_z_0_0_0_z_zz_yz_yz = buffer_1000_pddd[1936];

    auto g_z_0_0_0_z_zz_yz_zz = buffer_1000_pddd[1937];

    auto g_z_0_0_0_z_zz_zz_xx = buffer_1000_pddd[1938];

    auto g_z_0_0_0_z_zz_zz_xy = buffer_1000_pddd[1939];

    auto g_z_0_0_0_z_zz_zz_xz = buffer_1000_pddd[1940];

    auto g_z_0_0_0_z_zz_zz_yy = buffer_1000_pddd[1941];

    auto g_z_0_0_0_z_zz_zz_yz = buffer_1000_pddd[1942];

    auto g_z_0_0_0_z_zz_zz_zz = buffer_1000_pddd[1943];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_xx_xx_xx, g_0_xx_xx_xy, g_0_xx_xx_xz, g_0_xx_xx_yy, g_0_xx_xx_yz, g_0_xx_xx_zz, g_x_0_0_0_x_xx_xx_xx, g_x_0_0_0_x_xx_xx_xy, g_x_0_0_0_x_xx_xx_xz, g_x_0_0_0_x_xx_xx_yy, g_x_0_0_0_x_xx_xx_yz, g_x_0_0_0_x_xx_xx_zz, g_xx_xx_xx_xx, g_xx_xx_xx_xy, g_xx_xx_xx_xz, g_xx_xx_xx_yy, g_xx_xx_xx_yz, g_xx_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xx_xx_xx[i] = -g_0_xx_xx_xx[i] + 2.0 * g_xx_xx_xx_xx[i] * a_exp;

        g_x_0_0_0_x_xx_xx_xy[i] = -g_0_xx_xx_xy[i] + 2.0 * g_xx_xx_xx_xy[i] * a_exp;

        g_x_0_0_0_x_xx_xx_xz[i] = -g_0_xx_xx_xz[i] + 2.0 * g_xx_xx_xx_xz[i] * a_exp;

        g_x_0_0_0_x_xx_xx_yy[i] = -g_0_xx_xx_yy[i] + 2.0 * g_xx_xx_xx_yy[i] * a_exp;

        g_x_0_0_0_x_xx_xx_yz[i] = -g_0_xx_xx_yz[i] + 2.0 * g_xx_xx_xx_yz[i] * a_exp;

        g_x_0_0_0_x_xx_xx_zz[i] = -g_0_xx_xx_zz[i] + 2.0 * g_xx_xx_xx_zz[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_xx_xy_xx, g_0_xx_xy_xy, g_0_xx_xy_xz, g_0_xx_xy_yy, g_0_xx_xy_yz, g_0_xx_xy_zz, g_x_0_0_0_x_xx_xy_xx, g_x_0_0_0_x_xx_xy_xy, g_x_0_0_0_x_xx_xy_xz, g_x_0_0_0_x_xx_xy_yy, g_x_0_0_0_x_xx_xy_yz, g_x_0_0_0_x_xx_xy_zz, g_xx_xx_xy_xx, g_xx_xx_xy_xy, g_xx_xx_xy_xz, g_xx_xx_xy_yy, g_xx_xx_xy_yz, g_xx_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xx_xy_xx[i] = -g_0_xx_xy_xx[i] + 2.0 * g_xx_xx_xy_xx[i] * a_exp;

        g_x_0_0_0_x_xx_xy_xy[i] = -g_0_xx_xy_xy[i] + 2.0 * g_xx_xx_xy_xy[i] * a_exp;

        g_x_0_0_0_x_xx_xy_xz[i] = -g_0_xx_xy_xz[i] + 2.0 * g_xx_xx_xy_xz[i] * a_exp;

        g_x_0_0_0_x_xx_xy_yy[i] = -g_0_xx_xy_yy[i] + 2.0 * g_xx_xx_xy_yy[i] * a_exp;

        g_x_0_0_0_x_xx_xy_yz[i] = -g_0_xx_xy_yz[i] + 2.0 * g_xx_xx_xy_yz[i] * a_exp;

        g_x_0_0_0_x_xx_xy_zz[i] = -g_0_xx_xy_zz[i] + 2.0 * g_xx_xx_xy_zz[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_xx_xz_xx, g_0_xx_xz_xy, g_0_xx_xz_xz, g_0_xx_xz_yy, g_0_xx_xz_yz, g_0_xx_xz_zz, g_x_0_0_0_x_xx_xz_xx, g_x_0_0_0_x_xx_xz_xy, g_x_0_0_0_x_xx_xz_xz, g_x_0_0_0_x_xx_xz_yy, g_x_0_0_0_x_xx_xz_yz, g_x_0_0_0_x_xx_xz_zz, g_xx_xx_xz_xx, g_xx_xx_xz_xy, g_xx_xx_xz_xz, g_xx_xx_xz_yy, g_xx_xx_xz_yz, g_xx_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xx_xz_xx[i] = -g_0_xx_xz_xx[i] + 2.0 * g_xx_xx_xz_xx[i] * a_exp;

        g_x_0_0_0_x_xx_xz_xy[i] = -g_0_xx_xz_xy[i] + 2.0 * g_xx_xx_xz_xy[i] * a_exp;

        g_x_0_0_0_x_xx_xz_xz[i] = -g_0_xx_xz_xz[i] + 2.0 * g_xx_xx_xz_xz[i] * a_exp;

        g_x_0_0_0_x_xx_xz_yy[i] = -g_0_xx_xz_yy[i] + 2.0 * g_xx_xx_xz_yy[i] * a_exp;

        g_x_0_0_0_x_xx_xz_yz[i] = -g_0_xx_xz_yz[i] + 2.0 * g_xx_xx_xz_yz[i] * a_exp;

        g_x_0_0_0_x_xx_xz_zz[i] = -g_0_xx_xz_zz[i] + 2.0 * g_xx_xx_xz_zz[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_xx_yy_xx, g_0_xx_yy_xy, g_0_xx_yy_xz, g_0_xx_yy_yy, g_0_xx_yy_yz, g_0_xx_yy_zz, g_x_0_0_0_x_xx_yy_xx, g_x_0_0_0_x_xx_yy_xy, g_x_0_0_0_x_xx_yy_xz, g_x_0_0_0_x_xx_yy_yy, g_x_0_0_0_x_xx_yy_yz, g_x_0_0_0_x_xx_yy_zz, g_xx_xx_yy_xx, g_xx_xx_yy_xy, g_xx_xx_yy_xz, g_xx_xx_yy_yy, g_xx_xx_yy_yz, g_xx_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xx_yy_xx[i] = -g_0_xx_yy_xx[i] + 2.0 * g_xx_xx_yy_xx[i] * a_exp;

        g_x_0_0_0_x_xx_yy_xy[i] = -g_0_xx_yy_xy[i] + 2.0 * g_xx_xx_yy_xy[i] * a_exp;

        g_x_0_0_0_x_xx_yy_xz[i] = -g_0_xx_yy_xz[i] + 2.0 * g_xx_xx_yy_xz[i] * a_exp;

        g_x_0_0_0_x_xx_yy_yy[i] = -g_0_xx_yy_yy[i] + 2.0 * g_xx_xx_yy_yy[i] * a_exp;

        g_x_0_0_0_x_xx_yy_yz[i] = -g_0_xx_yy_yz[i] + 2.0 * g_xx_xx_yy_yz[i] * a_exp;

        g_x_0_0_0_x_xx_yy_zz[i] = -g_0_xx_yy_zz[i] + 2.0 * g_xx_xx_yy_zz[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_xx_yz_xx, g_0_xx_yz_xy, g_0_xx_yz_xz, g_0_xx_yz_yy, g_0_xx_yz_yz, g_0_xx_yz_zz, g_x_0_0_0_x_xx_yz_xx, g_x_0_0_0_x_xx_yz_xy, g_x_0_0_0_x_xx_yz_xz, g_x_0_0_0_x_xx_yz_yy, g_x_0_0_0_x_xx_yz_yz, g_x_0_0_0_x_xx_yz_zz, g_xx_xx_yz_xx, g_xx_xx_yz_xy, g_xx_xx_yz_xz, g_xx_xx_yz_yy, g_xx_xx_yz_yz, g_xx_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xx_yz_xx[i] = -g_0_xx_yz_xx[i] + 2.0 * g_xx_xx_yz_xx[i] * a_exp;

        g_x_0_0_0_x_xx_yz_xy[i] = -g_0_xx_yz_xy[i] + 2.0 * g_xx_xx_yz_xy[i] * a_exp;

        g_x_0_0_0_x_xx_yz_xz[i] = -g_0_xx_yz_xz[i] + 2.0 * g_xx_xx_yz_xz[i] * a_exp;

        g_x_0_0_0_x_xx_yz_yy[i] = -g_0_xx_yz_yy[i] + 2.0 * g_xx_xx_yz_yy[i] * a_exp;

        g_x_0_0_0_x_xx_yz_yz[i] = -g_0_xx_yz_yz[i] + 2.0 * g_xx_xx_yz_yz[i] * a_exp;

        g_x_0_0_0_x_xx_yz_zz[i] = -g_0_xx_yz_zz[i] + 2.0 * g_xx_xx_yz_zz[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_xx_zz_xx, g_0_xx_zz_xy, g_0_xx_zz_xz, g_0_xx_zz_yy, g_0_xx_zz_yz, g_0_xx_zz_zz, g_x_0_0_0_x_xx_zz_xx, g_x_0_0_0_x_xx_zz_xy, g_x_0_0_0_x_xx_zz_xz, g_x_0_0_0_x_xx_zz_yy, g_x_0_0_0_x_xx_zz_yz, g_x_0_0_0_x_xx_zz_zz, g_xx_xx_zz_xx, g_xx_xx_zz_xy, g_xx_xx_zz_xz, g_xx_xx_zz_yy, g_xx_xx_zz_yz, g_xx_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xx_zz_xx[i] = -g_0_xx_zz_xx[i] + 2.0 * g_xx_xx_zz_xx[i] * a_exp;

        g_x_0_0_0_x_xx_zz_xy[i] = -g_0_xx_zz_xy[i] + 2.0 * g_xx_xx_zz_xy[i] * a_exp;

        g_x_0_0_0_x_xx_zz_xz[i] = -g_0_xx_zz_xz[i] + 2.0 * g_xx_xx_zz_xz[i] * a_exp;

        g_x_0_0_0_x_xx_zz_yy[i] = -g_0_xx_zz_yy[i] + 2.0 * g_xx_xx_zz_yy[i] * a_exp;

        g_x_0_0_0_x_xx_zz_yz[i] = -g_0_xx_zz_yz[i] + 2.0 * g_xx_xx_zz_yz[i] * a_exp;

        g_x_0_0_0_x_xx_zz_zz[i] = -g_0_xx_zz_zz[i] + 2.0 * g_xx_xx_zz_zz[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_0_xy_xx_xx, g_0_xy_xx_xy, g_0_xy_xx_xz, g_0_xy_xx_yy, g_0_xy_xx_yz, g_0_xy_xx_zz, g_x_0_0_0_x_xy_xx_xx, g_x_0_0_0_x_xy_xx_xy, g_x_0_0_0_x_xy_xx_xz, g_x_0_0_0_x_xy_xx_yy, g_x_0_0_0_x_xy_xx_yz, g_x_0_0_0_x_xy_xx_zz, g_xx_xy_xx_xx, g_xx_xy_xx_xy, g_xx_xy_xx_xz, g_xx_xy_xx_yy, g_xx_xy_xx_yz, g_xx_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xy_xx_xx[i] = -g_0_xy_xx_xx[i] + 2.0 * g_xx_xy_xx_xx[i] * a_exp;

        g_x_0_0_0_x_xy_xx_xy[i] = -g_0_xy_xx_xy[i] + 2.0 * g_xx_xy_xx_xy[i] * a_exp;

        g_x_0_0_0_x_xy_xx_xz[i] = -g_0_xy_xx_xz[i] + 2.0 * g_xx_xy_xx_xz[i] * a_exp;

        g_x_0_0_0_x_xy_xx_yy[i] = -g_0_xy_xx_yy[i] + 2.0 * g_xx_xy_xx_yy[i] * a_exp;

        g_x_0_0_0_x_xy_xx_yz[i] = -g_0_xy_xx_yz[i] + 2.0 * g_xx_xy_xx_yz[i] * a_exp;

        g_x_0_0_0_x_xy_xx_zz[i] = -g_0_xy_xx_zz[i] + 2.0 * g_xx_xy_xx_zz[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_0_xy_xy_xx, g_0_xy_xy_xy, g_0_xy_xy_xz, g_0_xy_xy_yy, g_0_xy_xy_yz, g_0_xy_xy_zz, g_x_0_0_0_x_xy_xy_xx, g_x_0_0_0_x_xy_xy_xy, g_x_0_0_0_x_xy_xy_xz, g_x_0_0_0_x_xy_xy_yy, g_x_0_0_0_x_xy_xy_yz, g_x_0_0_0_x_xy_xy_zz, g_xx_xy_xy_xx, g_xx_xy_xy_xy, g_xx_xy_xy_xz, g_xx_xy_xy_yy, g_xx_xy_xy_yz, g_xx_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xy_xy_xx[i] = -g_0_xy_xy_xx[i] + 2.0 * g_xx_xy_xy_xx[i] * a_exp;

        g_x_0_0_0_x_xy_xy_xy[i] = -g_0_xy_xy_xy[i] + 2.0 * g_xx_xy_xy_xy[i] * a_exp;

        g_x_0_0_0_x_xy_xy_xz[i] = -g_0_xy_xy_xz[i] + 2.0 * g_xx_xy_xy_xz[i] * a_exp;

        g_x_0_0_0_x_xy_xy_yy[i] = -g_0_xy_xy_yy[i] + 2.0 * g_xx_xy_xy_yy[i] * a_exp;

        g_x_0_0_0_x_xy_xy_yz[i] = -g_0_xy_xy_yz[i] + 2.0 * g_xx_xy_xy_yz[i] * a_exp;

        g_x_0_0_0_x_xy_xy_zz[i] = -g_0_xy_xy_zz[i] + 2.0 * g_xx_xy_xy_zz[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_0_xy_xz_xx, g_0_xy_xz_xy, g_0_xy_xz_xz, g_0_xy_xz_yy, g_0_xy_xz_yz, g_0_xy_xz_zz, g_x_0_0_0_x_xy_xz_xx, g_x_0_0_0_x_xy_xz_xy, g_x_0_0_0_x_xy_xz_xz, g_x_0_0_0_x_xy_xz_yy, g_x_0_0_0_x_xy_xz_yz, g_x_0_0_0_x_xy_xz_zz, g_xx_xy_xz_xx, g_xx_xy_xz_xy, g_xx_xy_xz_xz, g_xx_xy_xz_yy, g_xx_xy_xz_yz, g_xx_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xy_xz_xx[i] = -g_0_xy_xz_xx[i] + 2.0 * g_xx_xy_xz_xx[i] * a_exp;

        g_x_0_0_0_x_xy_xz_xy[i] = -g_0_xy_xz_xy[i] + 2.0 * g_xx_xy_xz_xy[i] * a_exp;

        g_x_0_0_0_x_xy_xz_xz[i] = -g_0_xy_xz_xz[i] + 2.0 * g_xx_xy_xz_xz[i] * a_exp;

        g_x_0_0_0_x_xy_xz_yy[i] = -g_0_xy_xz_yy[i] + 2.0 * g_xx_xy_xz_yy[i] * a_exp;

        g_x_0_0_0_x_xy_xz_yz[i] = -g_0_xy_xz_yz[i] + 2.0 * g_xx_xy_xz_yz[i] * a_exp;

        g_x_0_0_0_x_xy_xz_zz[i] = -g_0_xy_xz_zz[i] + 2.0 * g_xx_xy_xz_zz[i] * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_0_xy_yy_xx, g_0_xy_yy_xy, g_0_xy_yy_xz, g_0_xy_yy_yy, g_0_xy_yy_yz, g_0_xy_yy_zz, g_x_0_0_0_x_xy_yy_xx, g_x_0_0_0_x_xy_yy_xy, g_x_0_0_0_x_xy_yy_xz, g_x_0_0_0_x_xy_yy_yy, g_x_0_0_0_x_xy_yy_yz, g_x_0_0_0_x_xy_yy_zz, g_xx_xy_yy_xx, g_xx_xy_yy_xy, g_xx_xy_yy_xz, g_xx_xy_yy_yy, g_xx_xy_yy_yz, g_xx_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xy_yy_xx[i] = -g_0_xy_yy_xx[i] + 2.0 * g_xx_xy_yy_xx[i] * a_exp;

        g_x_0_0_0_x_xy_yy_xy[i] = -g_0_xy_yy_xy[i] + 2.0 * g_xx_xy_yy_xy[i] * a_exp;

        g_x_0_0_0_x_xy_yy_xz[i] = -g_0_xy_yy_xz[i] + 2.0 * g_xx_xy_yy_xz[i] * a_exp;

        g_x_0_0_0_x_xy_yy_yy[i] = -g_0_xy_yy_yy[i] + 2.0 * g_xx_xy_yy_yy[i] * a_exp;

        g_x_0_0_0_x_xy_yy_yz[i] = -g_0_xy_yy_yz[i] + 2.0 * g_xx_xy_yy_yz[i] * a_exp;

        g_x_0_0_0_x_xy_yy_zz[i] = -g_0_xy_yy_zz[i] + 2.0 * g_xx_xy_yy_zz[i] * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_0_xy_yz_xx, g_0_xy_yz_xy, g_0_xy_yz_xz, g_0_xy_yz_yy, g_0_xy_yz_yz, g_0_xy_yz_zz, g_x_0_0_0_x_xy_yz_xx, g_x_0_0_0_x_xy_yz_xy, g_x_0_0_0_x_xy_yz_xz, g_x_0_0_0_x_xy_yz_yy, g_x_0_0_0_x_xy_yz_yz, g_x_0_0_0_x_xy_yz_zz, g_xx_xy_yz_xx, g_xx_xy_yz_xy, g_xx_xy_yz_xz, g_xx_xy_yz_yy, g_xx_xy_yz_yz, g_xx_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xy_yz_xx[i] = -g_0_xy_yz_xx[i] + 2.0 * g_xx_xy_yz_xx[i] * a_exp;

        g_x_0_0_0_x_xy_yz_xy[i] = -g_0_xy_yz_xy[i] + 2.0 * g_xx_xy_yz_xy[i] * a_exp;

        g_x_0_0_0_x_xy_yz_xz[i] = -g_0_xy_yz_xz[i] + 2.0 * g_xx_xy_yz_xz[i] * a_exp;

        g_x_0_0_0_x_xy_yz_yy[i] = -g_0_xy_yz_yy[i] + 2.0 * g_xx_xy_yz_yy[i] * a_exp;

        g_x_0_0_0_x_xy_yz_yz[i] = -g_0_xy_yz_yz[i] + 2.0 * g_xx_xy_yz_yz[i] * a_exp;

        g_x_0_0_0_x_xy_yz_zz[i] = -g_0_xy_yz_zz[i] + 2.0 * g_xx_xy_yz_zz[i] * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_0_xy_zz_xx, g_0_xy_zz_xy, g_0_xy_zz_xz, g_0_xy_zz_yy, g_0_xy_zz_yz, g_0_xy_zz_zz, g_x_0_0_0_x_xy_zz_xx, g_x_0_0_0_x_xy_zz_xy, g_x_0_0_0_x_xy_zz_xz, g_x_0_0_0_x_xy_zz_yy, g_x_0_0_0_x_xy_zz_yz, g_x_0_0_0_x_xy_zz_zz, g_xx_xy_zz_xx, g_xx_xy_zz_xy, g_xx_xy_zz_xz, g_xx_xy_zz_yy, g_xx_xy_zz_yz, g_xx_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xy_zz_xx[i] = -g_0_xy_zz_xx[i] + 2.0 * g_xx_xy_zz_xx[i] * a_exp;

        g_x_0_0_0_x_xy_zz_xy[i] = -g_0_xy_zz_xy[i] + 2.0 * g_xx_xy_zz_xy[i] * a_exp;

        g_x_0_0_0_x_xy_zz_xz[i] = -g_0_xy_zz_xz[i] + 2.0 * g_xx_xy_zz_xz[i] * a_exp;

        g_x_0_0_0_x_xy_zz_yy[i] = -g_0_xy_zz_yy[i] + 2.0 * g_xx_xy_zz_yy[i] * a_exp;

        g_x_0_0_0_x_xy_zz_yz[i] = -g_0_xy_zz_yz[i] + 2.0 * g_xx_xy_zz_yz[i] * a_exp;

        g_x_0_0_0_x_xy_zz_zz[i] = -g_0_xy_zz_zz[i] + 2.0 * g_xx_xy_zz_zz[i] * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_0_xz_xx_xx, g_0_xz_xx_xy, g_0_xz_xx_xz, g_0_xz_xx_yy, g_0_xz_xx_yz, g_0_xz_xx_zz, g_x_0_0_0_x_xz_xx_xx, g_x_0_0_0_x_xz_xx_xy, g_x_0_0_0_x_xz_xx_xz, g_x_0_0_0_x_xz_xx_yy, g_x_0_0_0_x_xz_xx_yz, g_x_0_0_0_x_xz_xx_zz, g_xx_xz_xx_xx, g_xx_xz_xx_xy, g_xx_xz_xx_xz, g_xx_xz_xx_yy, g_xx_xz_xx_yz, g_xx_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xz_xx_xx[i] = -g_0_xz_xx_xx[i] + 2.0 * g_xx_xz_xx_xx[i] * a_exp;

        g_x_0_0_0_x_xz_xx_xy[i] = -g_0_xz_xx_xy[i] + 2.0 * g_xx_xz_xx_xy[i] * a_exp;

        g_x_0_0_0_x_xz_xx_xz[i] = -g_0_xz_xx_xz[i] + 2.0 * g_xx_xz_xx_xz[i] * a_exp;

        g_x_0_0_0_x_xz_xx_yy[i] = -g_0_xz_xx_yy[i] + 2.0 * g_xx_xz_xx_yy[i] * a_exp;

        g_x_0_0_0_x_xz_xx_yz[i] = -g_0_xz_xx_yz[i] + 2.0 * g_xx_xz_xx_yz[i] * a_exp;

        g_x_0_0_0_x_xz_xx_zz[i] = -g_0_xz_xx_zz[i] + 2.0 * g_xx_xz_xx_zz[i] * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_0_xz_xy_xx, g_0_xz_xy_xy, g_0_xz_xy_xz, g_0_xz_xy_yy, g_0_xz_xy_yz, g_0_xz_xy_zz, g_x_0_0_0_x_xz_xy_xx, g_x_0_0_0_x_xz_xy_xy, g_x_0_0_0_x_xz_xy_xz, g_x_0_0_0_x_xz_xy_yy, g_x_0_0_0_x_xz_xy_yz, g_x_0_0_0_x_xz_xy_zz, g_xx_xz_xy_xx, g_xx_xz_xy_xy, g_xx_xz_xy_xz, g_xx_xz_xy_yy, g_xx_xz_xy_yz, g_xx_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xz_xy_xx[i] = -g_0_xz_xy_xx[i] + 2.0 * g_xx_xz_xy_xx[i] * a_exp;

        g_x_0_0_0_x_xz_xy_xy[i] = -g_0_xz_xy_xy[i] + 2.0 * g_xx_xz_xy_xy[i] * a_exp;

        g_x_0_0_0_x_xz_xy_xz[i] = -g_0_xz_xy_xz[i] + 2.0 * g_xx_xz_xy_xz[i] * a_exp;

        g_x_0_0_0_x_xz_xy_yy[i] = -g_0_xz_xy_yy[i] + 2.0 * g_xx_xz_xy_yy[i] * a_exp;

        g_x_0_0_0_x_xz_xy_yz[i] = -g_0_xz_xy_yz[i] + 2.0 * g_xx_xz_xy_yz[i] * a_exp;

        g_x_0_0_0_x_xz_xy_zz[i] = -g_0_xz_xy_zz[i] + 2.0 * g_xx_xz_xy_zz[i] * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_0_xz_xz_xx, g_0_xz_xz_xy, g_0_xz_xz_xz, g_0_xz_xz_yy, g_0_xz_xz_yz, g_0_xz_xz_zz, g_x_0_0_0_x_xz_xz_xx, g_x_0_0_0_x_xz_xz_xy, g_x_0_0_0_x_xz_xz_xz, g_x_0_0_0_x_xz_xz_yy, g_x_0_0_0_x_xz_xz_yz, g_x_0_0_0_x_xz_xz_zz, g_xx_xz_xz_xx, g_xx_xz_xz_xy, g_xx_xz_xz_xz, g_xx_xz_xz_yy, g_xx_xz_xz_yz, g_xx_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xz_xz_xx[i] = -g_0_xz_xz_xx[i] + 2.0 * g_xx_xz_xz_xx[i] * a_exp;

        g_x_0_0_0_x_xz_xz_xy[i] = -g_0_xz_xz_xy[i] + 2.0 * g_xx_xz_xz_xy[i] * a_exp;

        g_x_0_0_0_x_xz_xz_xz[i] = -g_0_xz_xz_xz[i] + 2.0 * g_xx_xz_xz_xz[i] * a_exp;

        g_x_0_0_0_x_xz_xz_yy[i] = -g_0_xz_xz_yy[i] + 2.0 * g_xx_xz_xz_yy[i] * a_exp;

        g_x_0_0_0_x_xz_xz_yz[i] = -g_0_xz_xz_yz[i] + 2.0 * g_xx_xz_xz_yz[i] * a_exp;

        g_x_0_0_0_x_xz_xz_zz[i] = -g_0_xz_xz_zz[i] + 2.0 * g_xx_xz_xz_zz[i] * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_0_xz_yy_xx, g_0_xz_yy_xy, g_0_xz_yy_xz, g_0_xz_yy_yy, g_0_xz_yy_yz, g_0_xz_yy_zz, g_x_0_0_0_x_xz_yy_xx, g_x_0_0_0_x_xz_yy_xy, g_x_0_0_0_x_xz_yy_xz, g_x_0_0_0_x_xz_yy_yy, g_x_0_0_0_x_xz_yy_yz, g_x_0_0_0_x_xz_yy_zz, g_xx_xz_yy_xx, g_xx_xz_yy_xy, g_xx_xz_yy_xz, g_xx_xz_yy_yy, g_xx_xz_yy_yz, g_xx_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xz_yy_xx[i] = -g_0_xz_yy_xx[i] + 2.0 * g_xx_xz_yy_xx[i] * a_exp;

        g_x_0_0_0_x_xz_yy_xy[i] = -g_0_xz_yy_xy[i] + 2.0 * g_xx_xz_yy_xy[i] * a_exp;

        g_x_0_0_0_x_xz_yy_xz[i] = -g_0_xz_yy_xz[i] + 2.0 * g_xx_xz_yy_xz[i] * a_exp;

        g_x_0_0_0_x_xz_yy_yy[i] = -g_0_xz_yy_yy[i] + 2.0 * g_xx_xz_yy_yy[i] * a_exp;

        g_x_0_0_0_x_xz_yy_yz[i] = -g_0_xz_yy_yz[i] + 2.0 * g_xx_xz_yy_yz[i] * a_exp;

        g_x_0_0_0_x_xz_yy_zz[i] = -g_0_xz_yy_zz[i] + 2.0 * g_xx_xz_yy_zz[i] * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_0_xz_yz_xx, g_0_xz_yz_xy, g_0_xz_yz_xz, g_0_xz_yz_yy, g_0_xz_yz_yz, g_0_xz_yz_zz, g_x_0_0_0_x_xz_yz_xx, g_x_0_0_0_x_xz_yz_xy, g_x_0_0_0_x_xz_yz_xz, g_x_0_0_0_x_xz_yz_yy, g_x_0_0_0_x_xz_yz_yz, g_x_0_0_0_x_xz_yz_zz, g_xx_xz_yz_xx, g_xx_xz_yz_xy, g_xx_xz_yz_xz, g_xx_xz_yz_yy, g_xx_xz_yz_yz, g_xx_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xz_yz_xx[i] = -g_0_xz_yz_xx[i] + 2.0 * g_xx_xz_yz_xx[i] * a_exp;

        g_x_0_0_0_x_xz_yz_xy[i] = -g_0_xz_yz_xy[i] + 2.0 * g_xx_xz_yz_xy[i] * a_exp;

        g_x_0_0_0_x_xz_yz_xz[i] = -g_0_xz_yz_xz[i] + 2.0 * g_xx_xz_yz_xz[i] * a_exp;

        g_x_0_0_0_x_xz_yz_yy[i] = -g_0_xz_yz_yy[i] + 2.0 * g_xx_xz_yz_yy[i] * a_exp;

        g_x_0_0_0_x_xz_yz_yz[i] = -g_0_xz_yz_yz[i] + 2.0 * g_xx_xz_yz_yz[i] * a_exp;

        g_x_0_0_0_x_xz_yz_zz[i] = -g_0_xz_yz_zz[i] + 2.0 * g_xx_xz_yz_zz[i] * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_0_xz_zz_xx, g_0_xz_zz_xy, g_0_xz_zz_xz, g_0_xz_zz_yy, g_0_xz_zz_yz, g_0_xz_zz_zz, g_x_0_0_0_x_xz_zz_xx, g_x_0_0_0_x_xz_zz_xy, g_x_0_0_0_x_xz_zz_xz, g_x_0_0_0_x_xz_zz_yy, g_x_0_0_0_x_xz_zz_yz, g_x_0_0_0_x_xz_zz_zz, g_xx_xz_zz_xx, g_xx_xz_zz_xy, g_xx_xz_zz_xz, g_xx_xz_zz_yy, g_xx_xz_zz_yz, g_xx_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xz_zz_xx[i] = -g_0_xz_zz_xx[i] + 2.0 * g_xx_xz_zz_xx[i] * a_exp;

        g_x_0_0_0_x_xz_zz_xy[i] = -g_0_xz_zz_xy[i] + 2.0 * g_xx_xz_zz_xy[i] * a_exp;

        g_x_0_0_0_x_xz_zz_xz[i] = -g_0_xz_zz_xz[i] + 2.0 * g_xx_xz_zz_xz[i] * a_exp;

        g_x_0_0_0_x_xz_zz_yy[i] = -g_0_xz_zz_yy[i] + 2.0 * g_xx_xz_zz_yy[i] * a_exp;

        g_x_0_0_0_x_xz_zz_yz[i] = -g_0_xz_zz_yz[i] + 2.0 * g_xx_xz_zz_yz[i] * a_exp;

        g_x_0_0_0_x_xz_zz_zz[i] = -g_0_xz_zz_zz[i] + 2.0 * g_xx_xz_zz_zz[i] * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_0_yy_xx_xx, g_0_yy_xx_xy, g_0_yy_xx_xz, g_0_yy_xx_yy, g_0_yy_xx_yz, g_0_yy_xx_zz, g_x_0_0_0_x_yy_xx_xx, g_x_0_0_0_x_yy_xx_xy, g_x_0_0_0_x_yy_xx_xz, g_x_0_0_0_x_yy_xx_yy, g_x_0_0_0_x_yy_xx_yz, g_x_0_0_0_x_yy_xx_zz, g_xx_yy_xx_xx, g_xx_yy_xx_xy, g_xx_yy_xx_xz, g_xx_yy_xx_yy, g_xx_yy_xx_yz, g_xx_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yy_xx_xx[i] = -g_0_yy_xx_xx[i] + 2.0 * g_xx_yy_xx_xx[i] * a_exp;

        g_x_0_0_0_x_yy_xx_xy[i] = -g_0_yy_xx_xy[i] + 2.0 * g_xx_yy_xx_xy[i] * a_exp;

        g_x_0_0_0_x_yy_xx_xz[i] = -g_0_yy_xx_xz[i] + 2.0 * g_xx_yy_xx_xz[i] * a_exp;

        g_x_0_0_0_x_yy_xx_yy[i] = -g_0_yy_xx_yy[i] + 2.0 * g_xx_yy_xx_yy[i] * a_exp;

        g_x_0_0_0_x_yy_xx_yz[i] = -g_0_yy_xx_yz[i] + 2.0 * g_xx_yy_xx_yz[i] * a_exp;

        g_x_0_0_0_x_yy_xx_zz[i] = -g_0_yy_xx_zz[i] + 2.0 * g_xx_yy_xx_zz[i] * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_0_yy_xy_xx, g_0_yy_xy_xy, g_0_yy_xy_xz, g_0_yy_xy_yy, g_0_yy_xy_yz, g_0_yy_xy_zz, g_x_0_0_0_x_yy_xy_xx, g_x_0_0_0_x_yy_xy_xy, g_x_0_0_0_x_yy_xy_xz, g_x_0_0_0_x_yy_xy_yy, g_x_0_0_0_x_yy_xy_yz, g_x_0_0_0_x_yy_xy_zz, g_xx_yy_xy_xx, g_xx_yy_xy_xy, g_xx_yy_xy_xz, g_xx_yy_xy_yy, g_xx_yy_xy_yz, g_xx_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yy_xy_xx[i] = -g_0_yy_xy_xx[i] + 2.0 * g_xx_yy_xy_xx[i] * a_exp;

        g_x_0_0_0_x_yy_xy_xy[i] = -g_0_yy_xy_xy[i] + 2.0 * g_xx_yy_xy_xy[i] * a_exp;

        g_x_0_0_0_x_yy_xy_xz[i] = -g_0_yy_xy_xz[i] + 2.0 * g_xx_yy_xy_xz[i] * a_exp;

        g_x_0_0_0_x_yy_xy_yy[i] = -g_0_yy_xy_yy[i] + 2.0 * g_xx_yy_xy_yy[i] * a_exp;

        g_x_0_0_0_x_yy_xy_yz[i] = -g_0_yy_xy_yz[i] + 2.0 * g_xx_yy_xy_yz[i] * a_exp;

        g_x_0_0_0_x_yy_xy_zz[i] = -g_0_yy_xy_zz[i] + 2.0 * g_xx_yy_xy_zz[i] * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_0_yy_xz_xx, g_0_yy_xz_xy, g_0_yy_xz_xz, g_0_yy_xz_yy, g_0_yy_xz_yz, g_0_yy_xz_zz, g_x_0_0_0_x_yy_xz_xx, g_x_0_0_0_x_yy_xz_xy, g_x_0_0_0_x_yy_xz_xz, g_x_0_0_0_x_yy_xz_yy, g_x_0_0_0_x_yy_xz_yz, g_x_0_0_0_x_yy_xz_zz, g_xx_yy_xz_xx, g_xx_yy_xz_xy, g_xx_yy_xz_xz, g_xx_yy_xz_yy, g_xx_yy_xz_yz, g_xx_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yy_xz_xx[i] = -g_0_yy_xz_xx[i] + 2.0 * g_xx_yy_xz_xx[i] * a_exp;

        g_x_0_0_0_x_yy_xz_xy[i] = -g_0_yy_xz_xy[i] + 2.0 * g_xx_yy_xz_xy[i] * a_exp;

        g_x_0_0_0_x_yy_xz_xz[i] = -g_0_yy_xz_xz[i] + 2.0 * g_xx_yy_xz_xz[i] * a_exp;

        g_x_0_0_0_x_yy_xz_yy[i] = -g_0_yy_xz_yy[i] + 2.0 * g_xx_yy_xz_yy[i] * a_exp;

        g_x_0_0_0_x_yy_xz_yz[i] = -g_0_yy_xz_yz[i] + 2.0 * g_xx_yy_xz_yz[i] * a_exp;

        g_x_0_0_0_x_yy_xz_zz[i] = -g_0_yy_xz_zz[i] + 2.0 * g_xx_yy_xz_zz[i] * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_0_yy_yy_xx, g_0_yy_yy_xy, g_0_yy_yy_xz, g_0_yy_yy_yy, g_0_yy_yy_yz, g_0_yy_yy_zz, g_x_0_0_0_x_yy_yy_xx, g_x_0_0_0_x_yy_yy_xy, g_x_0_0_0_x_yy_yy_xz, g_x_0_0_0_x_yy_yy_yy, g_x_0_0_0_x_yy_yy_yz, g_x_0_0_0_x_yy_yy_zz, g_xx_yy_yy_xx, g_xx_yy_yy_xy, g_xx_yy_yy_xz, g_xx_yy_yy_yy, g_xx_yy_yy_yz, g_xx_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yy_yy_xx[i] = -g_0_yy_yy_xx[i] + 2.0 * g_xx_yy_yy_xx[i] * a_exp;

        g_x_0_0_0_x_yy_yy_xy[i] = -g_0_yy_yy_xy[i] + 2.0 * g_xx_yy_yy_xy[i] * a_exp;

        g_x_0_0_0_x_yy_yy_xz[i] = -g_0_yy_yy_xz[i] + 2.0 * g_xx_yy_yy_xz[i] * a_exp;

        g_x_0_0_0_x_yy_yy_yy[i] = -g_0_yy_yy_yy[i] + 2.0 * g_xx_yy_yy_yy[i] * a_exp;

        g_x_0_0_0_x_yy_yy_yz[i] = -g_0_yy_yy_yz[i] + 2.0 * g_xx_yy_yy_yz[i] * a_exp;

        g_x_0_0_0_x_yy_yy_zz[i] = -g_0_yy_yy_zz[i] + 2.0 * g_xx_yy_yy_zz[i] * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_0_yy_yz_xx, g_0_yy_yz_xy, g_0_yy_yz_xz, g_0_yy_yz_yy, g_0_yy_yz_yz, g_0_yy_yz_zz, g_x_0_0_0_x_yy_yz_xx, g_x_0_0_0_x_yy_yz_xy, g_x_0_0_0_x_yy_yz_xz, g_x_0_0_0_x_yy_yz_yy, g_x_0_0_0_x_yy_yz_yz, g_x_0_0_0_x_yy_yz_zz, g_xx_yy_yz_xx, g_xx_yy_yz_xy, g_xx_yy_yz_xz, g_xx_yy_yz_yy, g_xx_yy_yz_yz, g_xx_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yy_yz_xx[i] = -g_0_yy_yz_xx[i] + 2.0 * g_xx_yy_yz_xx[i] * a_exp;

        g_x_0_0_0_x_yy_yz_xy[i] = -g_0_yy_yz_xy[i] + 2.0 * g_xx_yy_yz_xy[i] * a_exp;

        g_x_0_0_0_x_yy_yz_xz[i] = -g_0_yy_yz_xz[i] + 2.0 * g_xx_yy_yz_xz[i] * a_exp;

        g_x_0_0_0_x_yy_yz_yy[i] = -g_0_yy_yz_yy[i] + 2.0 * g_xx_yy_yz_yy[i] * a_exp;

        g_x_0_0_0_x_yy_yz_yz[i] = -g_0_yy_yz_yz[i] + 2.0 * g_xx_yy_yz_yz[i] * a_exp;

        g_x_0_0_0_x_yy_yz_zz[i] = -g_0_yy_yz_zz[i] + 2.0 * g_xx_yy_yz_zz[i] * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_0_yy_zz_xx, g_0_yy_zz_xy, g_0_yy_zz_xz, g_0_yy_zz_yy, g_0_yy_zz_yz, g_0_yy_zz_zz, g_x_0_0_0_x_yy_zz_xx, g_x_0_0_0_x_yy_zz_xy, g_x_0_0_0_x_yy_zz_xz, g_x_0_0_0_x_yy_zz_yy, g_x_0_0_0_x_yy_zz_yz, g_x_0_0_0_x_yy_zz_zz, g_xx_yy_zz_xx, g_xx_yy_zz_xy, g_xx_yy_zz_xz, g_xx_yy_zz_yy, g_xx_yy_zz_yz, g_xx_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yy_zz_xx[i] = -g_0_yy_zz_xx[i] + 2.0 * g_xx_yy_zz_xx[i] * a_exp;

        g_x_0_0_0_x_yy_zz_xy[i] = -g_0_yy_zz_xy[i] + 2.0 * g_xx_yy_zz_xy[i] * a_exp;

        g_x_0_0_0_x_yy_zz_xz[i] = -g_0_yy_zz_xz[i] + 2.0 * g_xx_yy_zz_xz[i] * a_exp;

        g_x_0_0_0_x_yy_zz_yy[i] = -g_0_yy_zz_yy[i] + 2.0 * g_xx_yy_zz_yy[i] * a_exp;

        g_x_0_0_0_x_yy_zz_yz[i] = -g_0_yy_zz_yz[i] + 2.0 * g_xx_yy_zz_yz[i] * a_exp;

        g_x_0_0_0_x_yy_zz_zz[i] = -g_0_yy_zz_zz[i] + 2.0 * g_xx_yy_zz_zz[i] * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_0_yz_xx_xx, g_0_yz_xx_xy, g_0_yz_xx_xz, g_0_yz_xx_yy, g_0_yz_xx_yz, g_0_yz_xx_zz, g_x_0_0_0_x_yz_xx_xx, g_x_0_0_0_x_yz_xx_xy, g_x_0_0_0_x_yz_xx_xz, g_x_0_0_0_x_yz_xx_yy, g_x_0_0_0_x_yz_xx_yz, g_x_0_0_0_x_yz_xx_zz, g_xx_yz_xx_xx, g_xx_yz_xx_xy, g_xx_yz_xx_xz, g_xx_yz_xx_yy, g_xx_yz_xx_yz, g_xx_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yz_xx_xx[i] = -g_0_yz_xx_xx[i] + 2.0 * g_xx_yz_xx_xx[i] * a_exp;

        g_x_0_0_0_x_yz_xx_xy[i] = -g_0_yz_xx_xy[i] + 2.0 * g_xx_yz_xx_xy[i] * a_exp;

        g_x_0_0_0_x_yz_xx_xz[i] = -g_0_yz_xx_xz[i] + 2.0 * g_xx_yz_xx_xz[i] * a_exp;

        g_x_0_0_0_x_yz_xx_yy[i] = -g_0_yz_xx_yy[i] + 2.0 * g_xx_yz_xx_yy[i] * a_exp;

        g_x_0_0_0_x_yz_xx_yz[i] = -g_0_yz_xx_yz[i] + 2.0 * g_xx_yz_xx_yz[i] * a_exp;

        g_x_0_0_0_x_yz_xx_zz[i] = -g_0_yz_xx_zz[i] + 2.0 * g_xx_yz_xx_zz[i] * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_0_yz_xy_xx, g_0_yz_xy_xy, g_0_yz_xy_xz, g_0_yz_xy_yy, g_0_yz_xy_yz, g_0_yz_xy_zz, g_x_0_0_0_x_yz_xy_xx, g_x_0_0_0_x_yz_xy_xy, g_x_0_0_0_x_yz_xy_xz, g_x_0_0_0_x_yz_xy_yy, g_x_0_0_0_x_yz_xy_yz, g_x_0_0_0_x_yz_xy_zz, g_xx_yz_xy_xx, g_xx_yz_xy_xy, g_xx_yz_xy_xz, g_xx_yz_xy_yy, g_xx_yz_xy_yz, g_xx_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yz_xy_xx[i] = -g_0_yz_xy_xx[i] + 2.0 * g_xx_yz_xy_xx[i] * a_exp;

        g_x_0_0_0_x_yz_xy_xy[i] = -g_0_yz_xy_xy[i] + 2.0 * g_xx_yz_xy_xy[i] * a_exp;

        g_x_0_0_0_x_yz_xy_xz[i] = -g_0_yz_xy_xz[i] + 2.0 * g_xx_yz_xy_xz[i] * a_exp;

        g_x_0_0_0_x_yz_xy_yy[i] = -g_0_yz_xy_yy[i] + 2.0 * g_xx_yz_xy_yy[i] * a_exp;

        g_x_0_0_0_x_yz_xy_yz[i] = -g_0_yz_xy_yz[i] + 2.0 * g_xx_yz_xy_yz[i] * a_exp;

        g_x_0_0_0_x_yz_xy_zz[i] = -g_0_yz_xy_zz[i] + 2.0 * g_xx_yz_xy_zz[i] * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_0_yz_xz_xx, g_0_yz_xz_xy, g_0_yz_xz_xz, g_0_yz_xz_yy, g_0_yz_xz_yz, g_0_yz_xz_zz, g_x_0_0_0_x_yz_xz_xx, g_x_0_0_0_x_yz_xz_xy, g_x_0_0_0_x_yz_xz_xz, g_x_0_0_0_x_yz_xz_yy, g_x_0_0_0_x_yz_xz_yz, g_x_0_0_0_x_yz_xz_zz, g_xx_yz_xz_xx, g_xx_yz_xz_xy, g_xx_yz_xz_xz, g_xx_yz_xz_yy, g_xx_yz_xz_yz, g_xx_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yz_xz_xx[i] = -g_0_yz_xz_xx[i] + 2.0 * g_xx_yz_xz_xx[i] * a_exp;

        g_x_0_0_0_x_yz_xz_xy[i] = -g_0_yz_xz_xy[i] + 2.0 * g_xx_yz_xz_xy[i] * a_exp;

        g_x_0_0_0_x_yz_xz_xz[i] = -g_0_yz_xz_xz[i] + 2.0 * g_xx_yz_xz_xz[i] * a_exp;

        g_x_0_0_0_x_yz_xz_yy[i] = -g_0_yz_xz_yy[i] + 2.0 * g_xx_yz_xz_yy[i] * a_exp;

        g_x_0_0_0_x_yz_xz_yz[i] = -g_0_yz_xz_yz[i] + 2.0 * g_xx_yz_xz_yz[i] * a_exp;

        g_x_0_0_0_x_yz_xz_zz[i] = -g_0_yz_xz_zz[i] + 2.0 * g_xx_yz_xz_zz[i] * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_0_yz_yy_xx, g_0_yz_yy_xy, g_0_yz_yy_xz, g_0_yz_yy_yy, g_0_yz_yy_yz, g_0_yz_yy_zz, g_x_0_0_0_x_yz_yy_xx, g_x_0_0_0_x_yz_yy_xy, g_x_0_0_0_x_yz_yy_xz, g_x_0_0_0_x_yz_yy_yy, g_x_0_0_0_x_yz_yy_yz, g_x_0_0_0_x_yz_yy_zz, g_xx_yz_yy_xx, g_xx_yz_yy_xy, g_xx_yz_yy_xz, g_xx_yz_yy_yy, g_xx_yz_yy_yz, g_xx_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yz_yy_xx[i] = -g_0_yz_yy_xx[i] + 2.0 * g_xx_yz_yy_xx[i] * a_exp;

        g_x_0_0_0_x_yz_yy_xy[i] = -g_0_yz_yy_xy[i] + 2.0 * g_xx_yz_yy_xy[i] * a_exp;

        g_x_0_0_0_x_yz_yy_xz[i] = -g_0_yz_yy_xz[i] + 2.0 * g_xx_yz_yy_xz[i] * a_exp;

        g_x_0_0_0_x_yz_yy_yy[i] = -g_0_yz_yy_yy[i] + 2.0 * g_xx_yz_yy_yy[i] * a_exp;

        g_x_0_0_0_x_yz_yy_yz[i] = -g_0_yz_yy_yz[i] + 2.0 * g_xx_yz_yy_yz[i] * a_exp;

        g_x_0_0_0_x_yz_yy_zz[i] = -g_0_yz_yy_zz[i] + 2.0 * g_xx_yz_yy_zz[i] * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_0_yz_yz_xx, g_0_yz_yz_xy, g_0_yz_yz_xz, g_0_yz_yz_yy, g_0_yz_yz_yz, g_0_yz_yz_zz, g_x_0_0_0_x_yz_yz_xx, g_x_0_0_0_x_yz_yz_xy, g_x_0_0_0_x_yz_yz_xz, g_x_0_0_0_x_yz_yz_yy, g_x_0_0_0_x_yz_yz_yz, g_x_0_0_0_x_yz_yz_zz, g_xx_yz_yz_xx, g_xx_yz_yz_xy, g_xx_yz_yz_xz, g_xx_yz_yz_yy, g_xx_yz_yz_yz, g_xx_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yz_yz_xx[i] = -g_0_yz_yz_xx[i] + 2.0 * g_xx_yz_yz_xx[i] * a_exp;

        g_x_0_0_0_x_yz_yz_xy[i] = -g_0_yz_yz_xy[i] + 2.0 * g_xx_yz_yz_xy[i] * a_exp;

        g_x_0_0_0_x_yz_yz_xz[i] = -g_0_yz_yz_xz[i] + 2.0 * g_xx_yz_yz_xz[i] * a_exp;

        g_x_0_0_0_x_yz_yz_yy[i] = -g_0_yz_yz_yy[i] + 2.0 * g_xx_yz_yz_yy[i] * a_exp;

        g_x_0_0_0_x_yz_yz_yz[i] = -g_0_yz_yz_yz[i] + 2.0 * g_xx_yz_yz_yz[i] * a_exp;

        g_x_0_0_0_x_yz_yz_zz[i] = -g_0_yz_yz_zz[i] + 2.0 * g_xx_yz_yz_zz[i] * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_0_yz_zz_xx, g_0_yz_zz_xy, g_0_yz_zz_xz, g_0_yz_zz_yy, g_0_yz_zz_yz, g_0_yz_zz_zz, g_x_0_0_0_x_yz_zz_xx, g_x_0_0_0_x_yz_zz_xy, g_x_0_0_0_x_yz_zz_xz, g_x_0_0_0_x_yz_zz_yy, g_x_0_0_0_x_yz_zz_yz, g_x_0_0_0_x_yz_zz_zz, g_xx_yz_zz_xx, g_xx_yz_zz_xy, g_xx_yz_zz_xz, g_xx_yz_zz_yy, g_xx_yz_zz_yz, g_xx_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yz_zz_xx[i] = -g_0_yz_zz_xx[i] + 2.0 * g_xx_yz_zz_xx[i] * a_exp;

        g_x_0_0_0_x_yz_zz_xy[i] = -g_0_yz_zz_xy[i] + 2.0 * g_xx_yz_zz_xy[i] * a_exp;

        g_x_0_0_0_x_yz_zz_xz[i] = -g_0_yz_zz_xz[i] + 2.0 * g_xx_yz_zz_xz[i] * a_exp;

        g_x_0_0_0_x_yz_zz_yy[i] = -g_0_yz_zz_yy[i] + 2.0 * g_xx_yz_zz_yy[i] * a_exp;

        g_x_0_0_0_x_yz_zz_yz[i] = -g_0_yz_zz_yz[i] + 2.0 * g_xx_yz_zz_yz[i] * a_exp;

        g_x_0_0_0_x_yz_zz_zz[i] = -g_0_yz_zz_zz[i] + 2.0 * g_xx_yz_zz_zz[i] * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_0_zz_xx_xx, g_0_zz_xx_xy, g_0_zz_xx_xz, g_0_zz_xx_yy, g_0_zz_xx_yz, g_0_zz_xx_zz, g_x_0_0_0_x_zz_xx_xx, g_x_0_0_0_x_zz_xx_xy, g_x_0_0_0_x_zz_xx_xz, g_x_0_0_0_x_zz_xx_yy, g_x_0_0_0_x_zz_xx_yz, g_x_0_0_0_x_zz_xx_zz, g_xx_zz_xx_xx, g_xx_zz_xx_xy, g_xx_zz_xx_xz, g_xx_zz_xx_yy, g_xx_zz_xx_yz, g_xx_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_zz_xx_xx[i] = -g_0_zz_xx_xx[i] + 2.0 * g_xx_zz_xx_xx[i] * a_exp;

        g_x_0_0_0_x_zz_xx_xy[i] = -g_0_zz_xx_xy[i] + 2.0 * g_xx_zz_xx_xy[i] * a_exp;

        g_x_0_0_0_x_zz_xx_xz[i] = -g_0_zz_xx_xz[i] + 2.0 * g_xx_zz_xx_xz[i] * a_exp;

        g_x_0_0_0_x_zz_xx_yy[i] = -g_0_zz_xx_yy[i] + 2.0 * g_xx_zz_xx_yy[i] * a_exp;

        g_x_0_0_0_x_zz_xx_yz[i] = -g_0_zz_xx_yz[i] + 2.0 * g_xx_zz_xx_yz[i] * a_exp;

        g_x_0_0_0_x_zz_xx_zz[i] = -g_0_zz_xx_zz[i] + 2.0 * g_xx_zz_xx_zz[i] * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_0_zz_xy_xx, g_0_zz_xy_xy, g_0_zz_xy_xz, g_0_zz_xy_yy, g_0_zz_xy_yz, g_0_zz_xy_zz, g_x_0_0_0_x_zz_xy_xx, g_x_0_0_0_x_zz_xy_xy, g_x_0_0_0_x_zz_xy_xz, g_x_0_0_0_x_zz_xy_yy, g_x_0_0_0_x_zz_xy_yz, g_x_0_0_0_x_zz_xy_zz, g_xx_zz_xy_xx, g_xx_zz_xy_xy, g_xx_zz_xy_xz, g_xx_zz_xy_yy, g_xx_zz_xy_yz, g_xx_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_zz_xy_xx[i] = -g_0_zz_xy_xx[i] + 2.0 * g_xx_zz_xy_xx[i] * a_exp;

        g_x_0_0_0_x_zz_xy_xy[i] = -g_0_zz_xy_xy[i] + 2.0 * g_xx_zz_xy_xy[i] * a_exp;

        g_x_0_0_0_x_zz_xy_xz[i] = -g_0_zz_xy_xz[i] + 2.0 * g_xx_zz_xy_xz[i] * a_exp;

        g_x_0_0_0_x_zz_xy_yy[i] = -g_0_zz_xy_yy[i] + 2.0 * g_xx_zz_xy_yy[i] * a_exp;

        g_x_0_0_0_x_zz_xy_yz[i] = -g_0_zz_xy_yz[i] + 2.0 * g_xx_zz_xy_yz[i] * a_exp;

        g_x_0_0_0_x_zz_xy_zz[i] = -g_0_zz_xy_zz[i] + 2.0 * g_xx_zz_xy_zz[i] * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_0_zz_xz_xx, g_0_zz_xz_xy, g_0_zz_xz_xz, g_0_zz_xz_yy, g_0_zz_xz_yz, g_0_zz_xz_zz, g_x_0_0_0_x_zz_xz_xx, g_x_0_0_0_x_zz_xz_xy, g_x_0_0_0_x_zz_xz_xz, g_x_0_0_0_x_zz_xz_yy, g_x_0_0_0_x_zz_xz_yz, g_x_0_0_0_x_zz_xz_zz, g_xx_zz_xz_xx, g_xx_zz_xz_xy, g_xx_zz_xz_xz, g_xx_zz_xz_yy, g_xx_zz_xz_yz, g_xx_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_zz_xz_xx[i] = -g_0_zz_xz_xx[i] + 2.0 * g_xx_zz_xz_xx[i] * a_exp;

        g_x_0_0_0_x_zz_xz_xy[i] = -g_0_zz_xz_xy[i] + 2.0 * g_xx_zz_xz_xy[i] * a_exp;

        g_x_0_0_0_x_zz_xz_xz[i] = -g_0_zz_xz_xz[i] + 2.0 * g_xx_zz_xz_xz[i] * a_exp;

        g_x_0_0_0_x_zz_xz_yy[i] = -g_0_zz_xz_yy[i] + 2.0 * g_xx_zz_xz_yy[i] * a_exp;

        g_x_0_0_0_x_zz_xz_yz[i] = -g_0_zz_xz_yz[i] + 2.0 * g_xx_zz_xz_yz[i] * a_exp;

        g_x_0_0_0_x_zz_xz_zz[i] = -g_0_zz_xz_zz[i] + 2.0 * g_xx_zz_xz_zz[i] * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_0_zz_yy_xx, g_0_zz_yy_xy, g_0_zz_yy_xz, g_0_zz_yy_yy, g_0_zz_yy_yz, g_0_zz_yy_zz, g_x_0_0_0_x_zz_yy_xx, g_x_0_0_0_x_zz_yy_xy, g_x_0_0_0_x_zz_yy_xz, g_x_0_0_0_x_zz_yy_yy, g_x_0_0_0_x_zz_yy_yz, g_x_0_0_0_x_zz_yy_zz, g_xx_zz_yy_xx, g_xx_zz_yy_xy, g_xx_zz_yy_xz, g_xx_zz_yy_yy, g_xx_zz_yy_yz, g_xx_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_zz_yy_xx[i] = -g_0_zz_yy_xx[i] + 2.0 * g_xx_zz_yy_xx[i] * a_exp;

        g_x_0_0_0_x_zz_yy_xy[i] = -g_0_zz_yy_xy[i] + 2.0 * g_xx_zz_yy_xy[i] * a_exp;

        g_x_0_0_0_x_zz_yy_xz[i] = -g_0_zz_yy_xz[i] + 2.0 * g_xx_zz_yy_xz[i] * a_exp;

        g_x_0_0_0_x_zz_yy_yy[i] = -g_0_zz_yy_yy[i] + 2.0 * g_xx_zz_yy_yy[i] * a_exp;

        g_x_0_0_0_x_zz_yy_yz[i] = -g_0_zz_yy_yz[i] + 2.0 * g_xx_zz_yy_yz[i] * a_exp;

        g_x_0_0_0_x_zz_yy_zz[i] = -g_0_zz_yy_zz[i] + 2.0 * g_xx_zz_yy_zz[i] * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_0_zz_yz_xx, g_0_zz_yz_xy, g_0_zz_yz_xz, g_0_zz_yz_yy, g_0_zz_yz_yz, g_0_zz_yz_zz, g_x_0_0_0_x_zz_yz_xx, g_x_0_0_0_x_zz_yz_xy, g_x_0_0_0_x_zz_yz_xz, g_x_0_0_0_x_zz_yz_yy, g_x_0_0_0_x_zz_yz_yz, g_x_0_0_0_x_zz_yz_zz, g_xx_zz_yz_xx, g_xx_zz_yz_xy, g_xx_zz_yz_xz, g_xx_zz_yz_yy, g_xx_zz_yz_yz, g_xx_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_zz_yz_xx[i] = -g_0_zz_yz_xx[i] + 2.0 * g_xx_zz_yz_xx[i] * a_exp;

        g_x_0_0_0_x_zz_yz_xy[i] = -g_0_zz_yz_xy[i] + 2.0 * g_xx_zz_yz_xy[i] * a_exp;

        g_x_0_0_0_x_zz_yz_xz[i] = -g_0_zz_yz_xz[i] + 2.0 * g_xx_zz_yz_xz[i] * a_exp;

        g_x_0_0_0_x_zz_yz_yy[i] = -g_0_zz_yz_yy[i] + 2.0 * g_xx_zz_yz_yy[i] * a_exp;

        g_x_0_0_0_x_zz_yz_yz[i] = -g_0_zz_yz_yz[i] + 2.0 * g_xx_zz_yz_yz[i] * a_exp;

        g_x_0_0_0_x_zz_yz_zz[i] = -g_0_zz_yz_zz[i] + 2.0 * g_xx_zz_yz_zz[i] * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_0_zz_zz_xx, g_0_zz_zz_xy, g_0_zz_zz_xz, g_0_zz_zz_yy, g_0_zz_zz_yz, g_0_zz_zz_zz, g_x_0_0_0_x_zz_zz_xx, g_x_0_0_0_x_zz_zz_xy, g_x_0_0_0_x_zz_zz_xz, g_x_0_0_0_x_zz_zz_yy, g_x_0_0_0_x_zz_zz_yz, g_x_0_0_0_x_zz_zz_zz, g_xx_zz_zz_xx, g_xx_zz_zz_xy, g_xx_zz_zz_xz, g_xx_zz_zz_yy, g_xx_zz_zz_yz, g_xx_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_zz_zz_xx[i] = -g_0_zz_zz_xx[i] + 2.0 * g_xx_zz_zz_xx[i] * a_exp;

        g_x_0_0_0_x_zz_zz_xy[i] = -g_0_zz_zz_xy[i] + 2.0 * g_xx_zz_zz_xy[i] * a_exp;

        g_x_0_0_0_x_zz_zz_xz[i] = -g_0_zz_zz_xz[i] + 2.0 * g_xx_zz_zz_xz[i] * a_exp;

        g_x_0_0_0_x_zz_zz_yy[i] = -g_0_zz_zz_yy[i] + 2.0 * g_xx_zz_zz_yy[i] * a_exp;

        g_x_0_0_0_x_zz_zz_yz[i] = -g_0_zz_zz_yz[i] + 2.0 * g_xx_zz_zz_yz[i] * a_exp;

        g_x_0_0_0_x_zz_zz_zz[i] = -g_0_zz_zz_zz[i] + 2.0 * g_xx_zz_zz_zz[i] * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_0_0_0_y_xx_xx_xx, g_x_0_0_0_y_xx_xx_xy, g_x_0_0_0_y_xx_xx_xz, g_x_0_0_0_y_xx_xx_yy, g_x_0_0_0_y_xx_xx_yz, g_x_0_0_0_y_xx_xx_zz, g_xy_xx_xx_xx, g_xy_xx_xx_xy, g_xy_xx_xx_xz, g_xy_xx_xx_yy, g_xy_xx_xx_yz, g_xy_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xx_xx_xx[i] = 2.0 * g_xy_xx_xx_xx[i] * a_exp;

        g_x_0_0_0_y_xx_xx_xy[i] = 2.0 * g_xy_xx_xx_xy[i] * a_exp;

        g_x_0_0_0_y_xx_xx_xz[i] = 2.0 * g_xy_xx_xx_xz[i] * a_exp;

        g_x_0_0_0_y_xx_xx_yy[i] = 2.0 * g_xy_xx_xx_yy[i] * a_exp;

        g_x_0_0_0_y_xx_xx_yz[i] = 2.0 * g_xy_xx_xx_yz[i] * a_exp;

        g_x_0_0_0_y_xx_xx_zz[i] = 2.0 * g_xy_xx_xx_zz[i] * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_0_0_0_y_xx_xy_xx, g_x_0_0_0_y_xx_xy_xy, g_x_0_0_0_y_xx_xy_xz, g_x_0_0_0_y_xx_xy_yy, g_x_0_0_0_y_xx_xy_yz, g_x_0_0_0_y_xx_xy_zz, g_xy_xx_xy_xx, g_xy_xx_xy_xy, g_xy_xx_xy_xz, g_xy_xx_xy_yy, g_xy_xx_xy_yz, g_xy_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xx_xy_xx[i] = 2.0 * g_xy_xx_xy_xx[i] * a_exp;

        g_x_0_0_0_y_xx_xy_xy[i] = 2.0 * g_xy_xx_xy_xy[i] * a_exp;

        g_x_0_0_0_y_xx_xy_xz[i] = 2.0 * g_xy_xx_xy_xz[i] * a_exp;

        g_x_0_0_0_y_xx_xy_yy[i] = 2.0 * g_xy_xx_xy_yy[i] * a_exp;

        g_x_0_0_0_y_xx_xy_yz[i] = 2.0 * g_xy_xx_xy_yz[i] * a_exp;

        g_x_0_0_0_y_xx_xy_zz[i] = 2.0 * g_xy_xx_xy_zz[i] * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_0_0_0_y_xx_xz_xx, g_x_0_0_0_y_xx_xz_xy, g_x_0_0_0_y_xx_xz_xz, g_x_0_0_0_y_xx_xz_yy, g_x_0_0_0_y_xx_xz_yz, g_x_0_0_0_y_xx_xz_zz, g_xy_xx_xz_xx, g_xy_xx_xz_xy, g_xy_xx_xz_xz, g_xy_xx_xz_yy, g_xy_xx_xz_yz, g_xy_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xx_xz_xx[i] = 2.0 * g_xy_xx_xz_xx[i] * a_exp;

        g_x_0_0_0_y_xx_xz_xy[i] = 2.0 * g_xy_xx_xz_xy[i] * a_exp;

        g_x_0_0_0_y_xx_xz_xz[i] = 2.0 * g_xy_xx_xz_xz[i] * a_exp;

        g_x_0_0_0_y_xx_xz_yy[i] = 2.0 * g_xy_xx_xz_yy[i] * a_exp;

        g_x_0_0_0_y_xx_xz_yz[i] = 2.0 * g_xy_xx_xz_yz[i] * a_exp;

        g_x_0_0_0_y_xx_xz_zz[i] = 2.0 * g_xy_xx_xz_zz[i] * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_0_0_0_y_xx_yy_xx, g_x_0_0_0_y_xx_yy_xy, g_x_0_0_0_y_xx_yy_xz, g_x_0_0_0_y_xx_yy_yy, g_x_0_0_0_y_xx_yy_yz, g_x_0_0_0_y_xx_yy_zz, g_xy_xx_yy_xx, g_xy_xx_yy_xy, g_xy_xx_yy_xz, g_xy_xx_yy_yy, g_xy_xx_yy_yz, g_xy_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xx_yy_xx[i] = 2.0 * g_xy_xx_yy_xx[i] * a_exp;

        g_x_0_0_0_y_xx_yy_xy[i] = 2.0 * g_xy_xx_yy_xy[i] * a_exp;

        g_x_0_0_0_y_xx_yy_xz[i] = 2.0 * g_xy_xx_yy_xz[i] * a_exp;

        g_x_0_0_0_y_xx_yy_yy[i] = 2.0 * g_xy_xx_yy_yy[i] * a_exp;

        g_x_0_0_0_y_xx_yy_yz[i] = 2.0 * g_xy_xx_yy_yz[i] * a_exp;

        g_x_0_0_0_y_xx_yy_zz[i] = 2.0 * g_xy_xx_yy_zz[i] * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_0_0_0_y_xx_yz_xx, g_x_0_0_0_y_xx_yz_xy, g_x_0_0_0_y_xx_yz_xz, g_x_0_0_0_y_xx_yz_yy, g_x_0_0_0_y_xx_yz_yz, g_x_0_0_0_y_xx_yz_zz, g_xy_xx_yz_xx, g_xy_xx_yz_xy, g_xy_xx_yz_xz, g_xy_xx_yz_yy, g_xy_xx_yz_yz, g_xy_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xx_yz_xx[i] = 2.0 * g_xy_xx_yz_xx[i] * a_exp;

        g_x_0_0_0_y_xx_yz_xy[i] = 2.0 * g_xy_xx_yz_xy[i] * a_exp;

        g_x_0_0_0_y_xx_yz_xz[i] = 2.0 * g_xy_xx_yz_xz[i] * a_exp;

        g_x_0_0_0_y_xx_yz_yy[i] = 2.0 * g_xy_xx_yz_yy[i] * a_exp;

        g_x_0_0_0_y_xx_yz_yz[i] = 2.0 * g_xy_xx_yz_yz[i] * a_exp;

        g_x_0_0_0_y_xx_yz_zz[i] = 2.0 * g_xy_xx_yz_zz[i] * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_0_0_0_y_xx_zz_xx, g_x_0_0_0_y_xx_zz_xy, g_x_0_0_0_y_xx_zz_xz, g_x_0_0_0_y_xx_zz_yy, g_x_0_0_0_y_xx_zz_yz, g_x_0_0_0_y_xx_zz_zz, g_xy_xx_zz_xx, g_xy_xx_zz_xy, g_xy_xx_zz_xz, g_xy_xx_zz_yy, g_xy_xx_zz_yz, g_xy_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xx_zz_xx[i] = 2.0 * g_xy_xx_zz_xx[i] * a_exp;

        g_x_0_0_0_y_xx_zz_xy[i] = 2.0 * g_xy_xx_zz_xy[i] * a_exp;

        g_x_0_0_0_y_xx_zz_xz[i] = 2.0 * g_xy_xx_zz_xz[i] * a_exp;

        g_x_0_0_0_y_xx_zz_yy[i] = 2.0 * g_xy_xx_zz_yy[i] * a_exp;

        g_x_0_0_0_y_xx_zz_yz[i] = 2.0 * g_xy_xx_zz_yz[i] * a_exp;

        g_x_0_0_0_y_xx_zz_zz[i] = 2.0 * g_xy_xx_zz_zz[i] * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_0_0_0_y_xy_xx_xx, g_x_0_0_0_y_xy_xx_xy, g_x_0_0_0_y_xy_xx_xz, g_x_0_0_0_y_xy_xx_yy, g_x_0_0_0_y_xy_xx_yz, g_x_0_0_0_y_xy_xx_zz, g_xy_xy_xx_xx, g_xy_xy_xx_xy, g_xy_xy_xx_xz, g_xy_xy_xx_yy, g_xy_xy_xx_yz, g_xy_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xy_xx_xx[i] = 2.0 * g_xy_xy_xx_xx[i] * a_exp;

        g_x_0_0_0_y_xy_xx_xy[i] = 2.0 * g_xy_xy_xx_xy[i] * a_exp;

        g_x_0_0_0_y_xy_xx_xz[i] = 2.0 * g_xy_xy_xx_xz[i] * a_exp;

        g_x_0_0_0_y_xy_xx_yy[i] = 2.0 * g_xy_xy_xx_yy[i] * a_exp;

        g_x_0_0_0_y_xy_xx_yz[i] = 2.0 * g_xy_xy_xx_yz[i] * a_exp;

        g_x_0_0_0_y_xy_xx_zz[i] = 2.0 * g_xy_xy_xx_zz[i] * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_0_0_0_y_xy_xy_xx, g_x_0_0_0_y_xy_xy_xy, g_x_0_0_0_y_xy_xy_xz, g_x_0_0_0_y_xy_xy_yy, g_x_0_0_0_y_xy_xy_yz, g_x_0_0_0_y_xy_xy_zz, g_xy_xy_xy_xx, g_xy_xy_xy_xy, g_xy_xy_xy_xz, g_xy_xy_xy_yy, g_xy_xy_xy_yz, g_xy_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xy_xy_xx[i] = 2.0 * g_xy_xy_xy_xx[i] * a_exp;

        g_x_0_0_0_y_xy_xy_xy[i] = 2.0 * g_xy_xy_xy_xy[i] * a_exp;

        g_x_0_0_0_y_xy_xy_xz[i] = 2.0 * g_xy_xy_xy_xz[i] * a_exp;

        g_x_0_0_0_y_xy_xy_yy[i] = 2.0 * g_xy_xy_xy_yy[i] * a_exp;

        g_x_0_0_0_y_xy_xy_yz[i] = 2.0 * g_xy_xy_xy_yz[i] * a_exp;

        g_x_0_0_0_y_xy_xy_zz[i] = 2.0 * g_xy_xy_xy_zz[i] * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_0_0_0_y_xy_xz_xx, g_x_0_0_0_y_xy_xz_xy, g_x_0_0_0_y_xy_xz_xz, g_x_0_0_0_y_xy_xz_yy, g_x_0_0_0_y_xy_xz_yz, g_x_0_0_0_y_xy_xz_zz, g_xy_xy_xz_xx, g_xy_xy_xz_xy, g_xy_xy_xz_xz, g_xy_xy_xz_yy, g_xy_xy_xz_yz, g_xy_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xy_xz_xx[i] = 2.0 * g_xy_xy_xz_xx[i] * a_exp;

        g_x_0_0_0_y_xy_xz_xy[i] = 2.0 * g_xy_xy_xz_xy[i] * a_exp;

        g_x_0_0_0_y_xy_xz_xz[i] = 2.0 * g_xy_xy_xz_xz[i] * a_exp;

        g_x_0_0_0_y_xy_xz_yy[i] = 2.0 * g_xy_xy_xz_yy[i] * a_exp;

        g_x_0_0_0_y_xy_xz_yz[i] = 2.0 * g_xy_xy_xz_yz[i] * a_exp;

        g_x_0_0_0_y_xy_xz_zz[i] = 2.0 * g_xy_xy_xz_zz[i] * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_0_0_0_y_xy_yy_xx, g_x_0_0_0_y_xy_yy_xy, g_x_0_0_0_y_xy_yy_xz, g_x_0_0_0_y_xy_yy_yy, g_x_0_0_0_y_xy_yy_yz, g_x_0_0_0_y_xy_yy_zz, g_xy_xy_yy_xx, g_xy_xy_yy_xy, g_xy_xy_yy_xz, g_xy_xy_yy_yy, g_xy_xy_yy_yz, g_xy_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xy_yy_xx[i] = 2.0 * g_xy_xy_yy_xx[i] * a_exp;

        g_x_0_0_0_y_xy_yy_xy[i] = 2.0 * g_xy_xy_yy_xy[i] * a_exp;

        g_x_0_0_0_y_xy_yy_xz[i] = 2.0 * g_xy_xy_yy_xz[i] * a_exp;

        g_x_0_0_0_y_xy_yy_yy[i] = 2.0 * g_xy_xy_yy_yy[i] * a_exp;

        g_x_0_0_0_y_xy_yy_yz[i] = 2.0 * g_xy_xy_yy_yz[i] * a_exp;

        g_x_0_0_0_y_xy_yy_zz[i] = 2.0 * g_xy_xy_yy_zz[i] * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_0_0_0_y_xy_yz_xx, g_x_0_0_0_y_xy_yz_xy, g_x_0_0_0_y_xy_yz_xz, g_x_0_0_0_y_xy_yz_yy, g_x_0_0_0_y_xy_yz_yz, g_x_0_0_0_y_xy_yz_zz, g_xy_xy_yz_xx, g_xy_xy_yz_xy, g_xy_xy_yz_xz, g_xy_xy_yz_yy, g_xy_xy_yz_yz, g_xy_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xy_yz_xx[i] = 2.0 * g_xy_xy_yz_xx[i] * a_exp;

        g_x_0_0_0_y_xy_yz_xy[i] = 2.0 * g_xy_xy_yz_xy[i] * a_exp;

        g_x_0_0_0_y_xy_yz_xz[i] = 2.0 * g_xy_xy_yz_xz[i] * a_exp;

        g_x_0_0_0_y_xy_yz_yy[i] = 2.0 * g_xy_xy_yz_yy[i] * a_exp;

        g_x_0_0_0_y_xy_yz_yz[i] = 2.0 * g_xy_xy_yz_yz[i] * a_exp;

        g_x_0_0_0_y_xy_yz_zz[i] = 2.0 * g_xy_xy_yz_zz[i] * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_0_0_0_y_xy_zz_xx, g_x_0_0_0_y_xy_zz_xy, g_x_0_0_0_y_xy_zz_xz, g_x_0_0_0_y_xy_zz_yy, g_x_0_0_0_y_xy_zz_yz, g_x_0_0_0_y_xy_zz_zz, g_xy_xy_zz_xx, g_xy_xy_zz_xy, g_xy_xy_zz_xz, g_xy_xy_zz_yy, g_xy_xy_zz_yz, g_xy_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xy_zz_xx[i] = 2.0 * g_xy_xy_zz_xx[i] * a_exp;

        g_x_0_0_0_y_xy_zz_xy[i] = 2.0 * g_xy_xy_zz_xy[i] * a_exp;

        g_x_0_0_0_y_xy_zz_xz[i] = 2.0 * g_xy_xy_zz_xz[i] * a_exp;

        g_x_0_0_0_y_xy_zz_yy[i] = 2.0 * g_xy_xy_zz_yy[i] * a_exp;

        g_x_0_0_0_y_xy_zz_yz[i] = 2.0 * g_xy_xy_zz_yz[i] * a_exp;

        g_x_0_0_0_y_xy_zz_zz[i] = 2.0 * g_xy_xy_zz_zz[i] * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_0_0_0_y_xz_xx_xx, g_x_0_0_0_y_xz_xx_xy, g_x_0_0_0_y_xz_xx_xz, g_x_0_0_0_y_xz_xx_yy, g_x_0_0_0_y_xz_xx_yz, g_x_0_0_0_y_xz_xx_zz, g_xy_xz_xx_xx, g_xy_xz_xx_xy, g_xy_xz_xx_xz, g_xy_xz_xx_yy, g_xy_xz_xx_yz, g_xy_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xz_xx_xx[i] = 2.0 * g_xy_xz_xx_xx[i] * a_exp;

        g_x_0_0_0_y_xz_xx_xy[i] = 2.0 * g_xy_xz_xx_xy[i] * a_exp;

        g_x_0_0_0_y_xz_xx_xz[i] = 2.0 * g_xy_xz_xx_xz[i] * a_exp;

        g_x_0_0_0_y_xz_xx_yy[i] = 2.0 * g_xy_xz_xx_yy[i] * a_exp;

        g_x_0_0_0_y_xz_xx_yz[i] = 2.0 * g_xy_xz_xx_yz[i] * a_exp;

        g_x_0_0_0_y_xz_xx_zz[i] = 2.0 * g_xy_xz_xx_zz[i] * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_0_0_0_y_xz_xy_xx, g_x_0_0_0_y_xz_xy_xy, g_x_0_0_0_y_xz_xy_xz, g_x_0_0_0_y_xz_xy_yy, g_x_0_0_0_y_xz_xy_yz, g_x_0_0_0_y_xz_xy_zz, g_xy_xz_xy_xx, g_xy_xz_xy_xy, g_xy_xz_xy_xz, g_xy_xz_xy_yy, g_xy_xz_xy_yz, g_xy_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xz_xy_xx[i] = 2.0 * g_xy_xz_xy_xx[i] * a_exp;

        g_x_0_0_0_y_xz_xy_xy[i] = 2.0 * g_xy_xz_xy_xy[i] * a_exp;

        g_x_0_0_0_y_xz_xy_xz[i] = 2.0 * g_xy_xz_xy_xz[i] * a_exp;

        g_x_0_0_0_y_xz_xy_yy[i] = 2.0 * g_xy_xz_xy_yy[i] * a_exp;

        g_x_0_0_0_y_xz_xy_yz[i] = 2.0 * g_xy_xz_xy_yz[i] * a_exp;

        g_x_0_0_0_y_xz_xy_zz[i] = 2.0 * g_xy_xz_xy_zz[i] * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_0_0_0_y_xz_xz_xx, g_x_0_0_0_y_xz_xz_xy, g_x_0_0_0_y_xz_xz_xz, g_x_0_0_0_y_xz_xz_yy, g_x_0_0_0_y_xz_xz_yz, g_x_0_0_0_y_xz_xz_zz, g_xy_xz_xz_xx, g_xy_xz_xz_xy, g_xy_xz_xz_xz, g_xy_xz_xz_yy, g_xy_xz_xz_yz, g_xy_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xz_xz_xx[i] = 2.0 * g_xy_xz_xz_xx[i] * a_exp;

        g_x_0_0_0_y_xz_xz_xy[i] = 2.0 * g_xy_xz_xz_xy[i] * a_exp;

        g_x_0_0_0_y_xz_xz_xz[i] = 2.0 * g_xy_xz_xz_xz[i] * a_exp;

        g_x_0_0_0_y_xz_xz_yy[i] = 2.0 * g_xy_xz_xz_yy[i] * a_exp;

        g_x_0_0_0_y_xz_xz_yz[i] = 2.0 * g_xy_xz_xz_yz[i] * a_exp;

        g_x_0_0_0_y_xz_xz_zz[i] = 2.0 * g_xy_xz_xz_zz[i] * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_0_0_0_y_xz_yy_xx, g_x_0_0_0_y_xz_yy_xy, g_x_0_0_0_y_xz_yy_xz, g_x_0_0_0_y_xz_yy_yy, g_x_0_0_0_y_xz_yy_yz, g_x_0_0_0_y_xz_yy_zz, g_xy_xz_yy_xx, g_xy_xz_yy_xy, g_xy_xz_yy_xz, g_xy_xz_yy_yy, g_xy_xz_yy_yz, g_xy_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xz_yy_xx[i] = 2.0 * g_xy_xz_yy_xx[i] * a_exp;

        g_x_0_0_0_y_xz_yy_xy[i] = 2.0 * g_xy_xz_yy_xy[i] * a_exp;

        g_x_0_0_0_y_xz_yy_xz[i] = 2.0 * g_xy_xz_yy_xz[i] * a_exp;

        g_x_0_0_0_y_xz_yy_yy[i] = 2.0 * g_xy_xz_yy_yy[i] * a_exp;

        g_x_0_0_0_y_xz_yy_yz[i] = 2.0 * g_xy_xz_yy_yz[i] * a_exp;

        g_x_0_0_0_y_xz_yy_zz[i] = 2.0 * g_xy_xz_yy_zz[i] * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_0_0_0_y_xz_yz_xx, g_x_0_0_0_y_xz_yz_xy, g_x_0_0_0_y_xz_yz_xz, g_x_0_0_0_y_xz_yz_yy, g_x_0_0_0_y_xz_yz_yz, g_x_0_0_0_y_xz_yz_zz, g_xy_xz_yz_xx, g_xy_xz_yz_xy, g_xy_xz_yz_xz, g_xy_xz_yz_yy, g_xy_xz_yz_yz, g_xy_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xz_yz_xx[i] = 2.0 * g_xy_xz_yz_xx[i] * a_exp;

        g_x_0_0_0_y_xz_yz_xy[i] = 2.0 * g_xy_xz_yz_xy[i] * a_exp;

        g_x_0_0_0_y_xz_yz_xz[i] = 2.0 * g_xy_xz_yz_xz[i] * a_exp;

        g_x_0_0_0_y_xz_yz_yy[i] = 2.0 * g_xy_xz_yz_yy[i] * a_exp;

        g_x_0_0_0_y_xz_yz_yz[i] = 2.0 * g_xy_xz_yz_yz[i] * a_exp;

        g_x_0_0_0_y_xz_yz_zz[i] = 2.0 * g_xy_xz_yz_zz[i] * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_0_0_0_y_xz_zz_xx, g_x_0_0_0_y_xz_zz_xy, g_x_0_0_0_y_xz_zz_xz, g_x_0_0_0_y_xz_zz_yy, g_x_0_0_0_y_xz_zz_yz, g_x_0_0_0_y_xz_zz_zz, g_xy_xz_zz_xx, g_xy_xz_zz_xy, g_xy_xz_zz_xz, g_xy_xz_zz_yy, g_xy_xz_zz_yz, g_xy_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xz_zz_xx[i] = 2.0 * g_xy_xz_zz_xx[i] * a_exp;

        g_x_0_0_0_y_xz_zz_xy[i] = 2.0 * g_xy_xz_zz_xy[i] * a_exp;

        g_x_0_0_0_y_xz_zz_xz[i] = 2.0 * g_xy_xz_zz_xz[i] * a_exp;

        g_x_0_0_0_y_xz_zz_yy[i] = 2.0 * g_xy_xz_zz_yy[i] * a_exp;

        g_x_0_0_0_y_xz_zz_yz[i] = 2.0 * g_xy_xz_zz_yz[i] * a_exp;

        g_x_0_0_0_y_xz_zz_zz[i] = 2.0 * g_xy_xz_zz_zz[i] * a_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_x_0_0_0_y_yy_xx_xx, g_x_0_0_0_y_yy_xx_xy, g_x_0_0_0_y_yy_xx_xz, g_x_0_0_0_y_yy_xx_yy, g_x_0_0_0_y_yy_xx_yz, g_x_0_0_0_y_yy_xx_zz, g_xy_yy_xx_xx, g_xy_yy_xx_xy, g_xy_yy_xx_xz, g_xy_yy_xx_yy, g_xy_yy_xx_yz, g_xy_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yy_xx_xx[i] = 2.0 * g_xy_yy_xx_xx[i] * a_exp;

        g_x_0_0_0_y_yy_xx_xy[i] = 2.0 * g_xy_yy_xx_xy[i] * a_exp;

        g_x_0_0_0_y_yy_xx_xz[i] = 2.0 * g_xy_yy_xx_xz[i] * a_exp;

        g_x_0_0_0_y_yy_xx_yy[i] = 2.0 * g_xy_yy_xx_yy[i] * a_exp;

        g_x_0_0_0_y_yy_xx_yz[i] = 2.0 * g_xy_yy_xx_yz[i] * a_exp;

        g_x_0_0_0_y_yy_xx_zz[i] = 2.0 * g_xy_yy_xx_zz[i] * a_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_x_0_0_0_y_yy_xy_xx, g_x_0_0_0_y_yy_xy_xy, g_x_0_0_0_y_yy_xy_xz, g_x_0_0_0_y_yy_xy_yy, g_x_0_0_0_y_yy_xy_yz, g_x_0_0_0_y_yy_xy_zz, g_xy_yy_xy_xx, g_xy_yy_xy_xy, g_xy_yy_xy_xz, g_xy_yy_xy_yy, g_xy_yy_xy_yz, g_xy_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yy_xy_xx[i] = 2.0 * g_xy_yy_xy_xx[i] * a_exp;

        g_x_0_0_0_y_yy_xy_xy[i] = 2.0 * g_xy_yy_xy_xy[i] * a_exp;

        g_x_0_0_0_y_yy_xy_xz[i] = 2.0 * g_xy_yy_xy_xz[i] * a_exp;

        g_x_0_0_0_y_yy_xy_yy[i] = 2.0 * g_xy_yy_xy_yy[i] * a_exp;

        g_x_0_0_0_y_yy_xy_yz[i] = 2.0 * g_xy_yy_xy_yz[i] * a_exp;

        g_x_0_0_0_y_yy_xy_zz[i] = 2.0 * g_xy_yy_xy_zz[i] * a_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_x_0_0_0_y_yy_xz_xx, g_x_0_0_0_y_yy_xz_xy, g_x_0_0_0_y_yy_xz_xz, g_x_0_0_0_y_yy_xz_yy, g_x_0_0_0_y_yy_xz_yz, g_x_0_0_0_y_yy_xz_zz, g_xy_yy_xz_xx, g_xy_yy_xz_xy, g_xy_yy_xz_xz, g_xy_yy_xz_yy, g_xy_yy_xz_yz, g_xy_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yy_xz_xx[i] = 2.0 * g_xy_yy_xz_xx[i] * a_exp;

        g_x_0_0_0_y_yy_xz_xy[i] = 2.0 * g_xy_yy_xz_xy[i] * a_exp;

        g_x_0_0_0_y_yy_xz_xz[i] = 2.0 * g_xy_yy_xz_xz[i] * a_exp;

        g_x_0_0_0_y_yy_xz_yy[i] = 2.0 * g_xy_yy_xz_yy[i] * a_exp;

        g_x_0_0_0_y_yy_xz_yz[i] = 2.0 * g_xy_yy_xz_yz[i] * a_exp;

        g_x_0_0_0_y_yy_xz_zz[i] = 2.0 * g_xy_yy_xz_zz[i] * a_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_x_0_0_0_y_yy_yy_xx, g_x_0_0_0_y_yy_yy_xy, g_x_0_0_0_y_yy_yy_xz, g_x_0_0_0_y_yy_yy_yy, g_x_0_0_0_y_yy_yy_yz, g_x_0_0_0_y_yy_yy_zz, g_xy_yy_yy_xx, g_xy_yy_yy_xy, g_xy_yy_yy_xz, g_xy_yy_yy_yy, g_xy_yy_yy_yz, g_xy_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yy_yy_xx[i] = 2.0 * g_xy_yy_yy_xx[i] * a_exp;

        g_x_0_0_0_y_yy_yy_xy[i] = 2.0 * g_xy_yy_yy_xy[i] * a_exp;

        g_x_0_0_0_y_yy_yy_xz[i] = 2.0 * g_xy_yy_yy_xz[i] * a_exp;

        g_x_0_0_0_y_yy_yy_yy[i] = 2.0 * g_xy_yy_yy_yy[i] * a_exp;

        g_x_0_0_0_y_yy_yy_yz[i] = 2.0 * g_xy_yy_yy_yz[i] * a_exp;

        g_x_0_0_0_y_yy_yy_zz[i] = 2.0 * g_xy_yy_yy_zz[i] * a_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_x_0_0_0_y_yy_yz_xx, g_x_0_0_0_y_yy_yz_xy, g_x_0_0_0_y_yy_yz_xz, g_x_0_0_0_y_yy_yz_yy, g_x_0_0_0_y_yy_yz_yz, g_x_0_0_0_y_yy_yz_zz, g_xy_yy_yz_xx, g_xy_yy_yz_xy, g_xy_yy_yz_xz, g_xy_yy_yz_yy, g_xy_yy_yz_yz, g_xy_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yy_yz_xx[i] = 2.0 * g_xy_yy_yz_xx[i] * a_exp;

        g_x_0_0_0_y_yy_yz_xy[i] = 2.0 * g_xy_yy_yz_xy[i] * a_exp;

        g_x_0_0_0_y_yy_yz_xz[i] = 2.0 * g_xy_yy_yz_xz[i] * a_exp;

        g_x_0_0_0_y_yy_yz_yy[i] = 2.0 * g_xy_yy_yz_yy[i] * a_exp;

        g_x_0_0_0_y_yy_yz_yz[i] = 2.0 * g_xy_yy_yz_yz[i] * a_exp;

        g_x_0_0_0_y_yy_yz_zz[i] = 2.0 * g_xy_yy_yz_zz[i] * a_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_x_0_0_0_y_yy_zz_xx, g_x_0_0_0_y_yy_zz_xy, g_x_0_0_0_y_yy_zz_xz, g_x_0_0_0_y_yy_zz_yy, g_x_0_0_0_y_yy_zz_yz, g_x_0_0_0_y_yy_zz_zz, g_xy_yy_zz_xx, g_xy_yy_zz_xy, g_xy_yy_zz_xz, g_xy_yy_zz_yy, g_xy_yy_zz_yz, g_xy_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yy_zz_xx[i] = 2.0 * g_xy_yy_zz_xx[i] * a_exp;

        g_x_0_0_0_y_yy_zz_xy[i] = 2.0 * g_xy_yy_zz_xy[i] * a_exp;

        g_x_0_0_0_y_yy_zz_xz[i] = 2.0 * g_xy_yy_zz_xz[i] * a_exp;

        g_x_0_0_0_y_yy_zz_yy[i] = 2.0 * g_xy_yy_zz_yy[i] * a_exp;

        g_x_0_0_0_y_yy_zz_yz[i] = 2.0 * g_xy_yy_zz_yz[i] * a_exp;

        g_x_0_0_0_y_yy_zz_zz[i] = 2.0 * g_xy_yy_zz_zz[i] * a_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_x_0_0_0_y_yz_xx_xx, g_x_0_0_0_y_yz_xx_xy, g_x_0_0_0_y_yz_xx_xz, g_x_0_0_0_y_yz_xx_yy, g_x_0_0_0_y_yz_xx_yz, g_x_0_0_0_y_yz_xx_zz, g_xy_yz_xx_xx, g_xy_yz_xx_xy, g_xy_yz_xx_xz, g_xy_yz_xx_yy, g_xy_yz_xx_yz, g_xy_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yz_xx_xx[i] = 2.0 * g_xy_yz_xx_xx[i] * a_exp;

        g_x_0_0_0_y_yz_xx_xy[i] = 2.0 * g_xy_yz_xx_xy[i] * a_exp;

        g_x_0_0_0_y_yz_xx_xz[i] = 2.0 * g_xy_yz_xx_xz[i] * a_exp;

        g_x_0_0_0_y_yz_xx_yy[i] = 2.0 * g_xy_yz_xx_yy[i] * a_exp;

        g_x_0_0_0_y_yz_xx_yz[i] = 2.0 * g_xy_yz_xx_yz[i] * a_exp;

        g_x_0_0_0_y_yz_xx_zz[i] = 2.0 * g_xy_yz_xx_zz[i] * a_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_x_0_0_0_y_yz_xy_xx, g_x_0_0_0_y_yz_xy_xy, g_x_0_0_0_y_yz_xy_xz, g_x_0_0_0_y_yz_xy_yy, g_x_0_0_0_y_yz_xy_yz, g_x_0_0_0_y_yz_xy_zz, g_xy_yz_xy_xx, g_xy_yz_xy_xy, g_xy_yz_xy_xz, g_xy_yz_xy_yy, g_xy_yz_xy_yz, g_xy_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yz_xy_xx[i] = 2.0 * g_xy_yz_xy_xx[i] * a_exp;

        g_x_0_0_0_y_yz_xy_xy[i] = 2.0 * g_xy_yz_xy_xy[i] * a_exp;

        g_x_0_0_0_y_yz_xy_xz[i] = 2.0 * g_xy_yz_xy_xz[i] * a_exp;

        g_x_0_0_0_y_yz_xy_yy[i] = 2.0 * g_xy_yz_xy_yy[i] * a_exp;

        g_x_0_0_0_y_yz_xy_yz[i] = 2.0 * g_xy_yz_xy_yz[i] * a_exp;

        g_x_0_0_0_y_yz_xy_zz[i] = 2.0 * g_xy_yz_xy_zz[i] * a_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_x_0_0_0_y_yz_xz_xx, g_x_0_0_0_y_yz_xz_xy, g_x_0_0_0_y_yz_xz_xz, g_x_0_0_0_y_yz_xz_yy, g_x_0_0_0_y_yz_xz_yz, g_x_0_0_0_y_yz_xz_zz, g_xy_yz_xz_xx, g_xy_yz_xz_xy, g_xy_yz_xz_xz, g_xy_yz_xz_yy, g_xy_yz_xz_yz, g_xy_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yz_xz_xx[i] = 2.0 * g_xy_yz_xz_xx[i] * a_exp;

        g_x_0_0_0_y_yz_xz_xy[i] = 2.0 * g_xy_yz_xz_xy[i] * a_exp;

        g_x_0_0_0_y_yz_xz_xz[i] = 2.0 * g_xy_yz_xz_xz[i] * a_exp;

        g_x_0_0_0_y_yz_xz_yy[i] = 2.0 * g_xy_yz_xz_yy[i] * a_exp;

        g_x_0_0_0_y_yz_xz_yz[i] = 2.0 * g_xy_yz_xz_yz[i] * a_exp;

        g_x_0_0_0_y_yz_xz_zz[i] = 2.0 * g_xy_yz_xz_zz[i] * a_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_x_0_0_0_y_yz_yy_xx, g_x_0_0_0_y_yz_yy_xy, g_x_0_0_0_y_yz_yy_xz, g_x_0_0_0_y_yz_yy_yy, g_x_0_0_0_y_yz_yy_yz, g_x_0_0_0_y_yz_yy_zz, g_xy_yz_yy_xx, g_xy_yz_yy_xy, g_xy_yz_yy_xz, g_xy_yz_yy_yy, g_xy_yz_yy_yz, g_xy_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yz_yy_xx[i] = 2.0 * g_xy_yz_yy_xx[i] * a_exp;

        g_x_0_0_0_y_yz_yy_xy[i] = 2.0 * g_xy_yz_yy_xy[i] * a_exp;

        g_x_0_0_0_y_yz_yy_xz[i] = 2.0 * g_xy_yz_yy_xz[i] * a_exp;

        g_x_0_0_0_y_yz_yy_yy[i] = 2.0 * g_xy_yz_yy_yy[i] * a_exp;

        g_x_0_0_0_y_yz_yy_yz[i] = 2.0 * g_xy_yz_yy_yz[i] * a_exp;

        g_x_0_0_0_y_yz_yy_zz[i] = 2.0 * g_xy_yz_yy_zz[i] * a_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_x_0_0_0_y_yz_yz_xx, g_x_0_0_0_y_yz_yz_xy, g_x_0_0_0_y_yz_yz_xz, g_x_0_0_0_y_yz_yz_yy, g_x_0_0_0_y_yz_yz_yz, g_x_0_0_0_y_yz_yz_zz, g_xy_yz_yz_xx, g_xy_yz_yz_xy, g_xy_yz_yz_xz, g_xy_yz_yz_yy, g_xy_yz_yz_yz, g_xy_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yz_yz_xx[i] = 2.0 * g_xy_yz_yz_xx[i] * a_exp;

        g_x_0_0_0_y_yz_yz_xy[i] = 2.0 * g_xy_yz_yz_xy[i] * a_exp;

        g_x_0_0_0_y_yz_yz_xz[i] = 2.0 * g_xy_yz_yz_xz[i] * a_exp;

        g_x_0_0_0_y_yz_yz_yy[i] = 2.0 * g_xy_yz_yz_yy[i] * a_exp;

        g_x_0_0_0_y_yz_yz_yz[i] = 2.0 * g_xy_yz_yz_yz[i] * a_exp;

        g_x_0_0_0_y_yz_yz_zz[i] = 2.0 * g_xy_yz_yz_zz[i] * a_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_x_0_0_0_y_yz_zz_xx, g_x_0_0_0_y_yz_zz_xy, g_x_0_0_0_y_yz_zz_xz, g_x_0_0_0_y_yz_zz_yy, g_x_0_0_0_y_yz_zz_yz, g_x_0_0_0_y_yz_zz_zz, g_xy_yz_zz_xx, g_xy_yz_zz_xy, g_xy_yz_zz_xz, g_xy_yz_zz_yy, g_xy_yz_zz_yz, g_xy_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yz_zz_xx[i] = 2.0 * g_xy_yz_zz_xx[i] * a_exp;

        g_x_0_0_0_y_yz_zz_xy[i] = 2.0 * g_xy_yz_zz_xy[i] * a_exp;

        g_x_0_0_0_y_yz_zz_xz[i] = 2.0 * g_xy_yz_zz_xz[i] * a_exp;

        g_x_0_0_0_y_yz_zz_yy[i] = 2.0 * g_xy_yz_zz_yy[i] * a_exp;

        g_x_0_0_0_y_yz_zz_yz[i] = 2.0 * g_xy_yz_zz_yz[i] * a_exp;

        g_x_0_0_0_y_yz_zz_zz[i] = 2.0 * g_xy_yz_zz_zz[i] * a_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_x_0_0_0_y_zz_xx_xx, g_x_0_0_0_y_zz_xx_xy, g_x_0_0_0_y_zz_xx_xz, g_x_0_0_0_y_zz_xx_yy, g_x_0_0_0_y_zz_xx_yz, g_x_0_0_0_y_zz_xx_zz, g_xy_zz_xx_xx, g_xy_zz_xx_xy, g_xy_zz_xx_xz, g_xy_zz_xx_yy, g_xy_zz_xx_yz, g_xy_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_zz_xx_xx[i] = 2.0 * g_xy_zz_xx_xx[i] * a_exp;

        g_x_0_0_0_y_zz_xx_xy[i] = 2.0 * g_xy_zz_xx_xy[i] * a_exp;

        g_x_0_0_0_y_zz_xx_xz[i] = 2.0 * g_xy_zz_xx_xz[i] * a_exp;

        g_x_0_0_0_y_zz_xx_yy[i] = 2.0 * g_xy_zz_xx_yy[i] * a_exp;

        g_x_0_0_0_y_zz_xx_yz[i] = 2.0 * g_xy_zz_xx_yz[i] * a_exp;

        g_x_0_0_0_y_zz_xx_zz[i] = 2.0 * g_xy_zz_xx_zz[i] * a_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_x_0_0_0_y_zz_xy_xx, g_x_0_0_0_y_zz_xy_xy, g_x_0_0_0_y_zz_xy_xz, g_x_0_0_0_y_zz_xy_yy, g_x_0_0_0_y_zz_xy_yz, g_x_0_0_0_y_zz_xy_zz, g_xy_zz_xy_xx, g_xy_zz_xy_xy, g_xy_zz_xy_xz, g_xy_zz_xy_yy, g_xy_zz_xy_yz, g_xy_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_zz_xy_xx[i] = 2.0 * g_xy_zz_xy_xx[i] * a_exp;

        g_x_0_0_0_y_zz_xy_xy[i] = 2.0 * g_xy_zz_xy_xy[i] * a_exp;

        g_x_0_0_0_y_zz_xy_xz[i] = 2.0 * g_xy_zz_xy_xz[i] * a_exp;

        g_x_0_0_0_y_zz_xy_yy[i] = 2.0 * g_xy_zz_xy_yy[i] * a_exp;

        g_x_0_0_0_y_zz_xy_yz[i] = 2.0 * g_xy_zz_xy_yz[i] * a_exp;

        g_x_0_0_0_y_zz_xy_zz[i] = 2.0 * g_xy_zz_xy_zz[i] * a_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_x_0_0_0_y_zz_xz_xx, g_x_0_0_0_y_zz_xz_xy, g_x_0_0_0_y_zz_xz_xz, g_x_0_0_0_y_zz_xz_yy, g_x_0_0_0_y_zz_xz_yz, g_x_0_0_0_y_zz_xz_zz, g_xy_zz_xz_xx, g_xy_zz_xz_xy, g_xy_zz_xz_xz, g_xy_zz_xz_yy, g_xy_zz_xz_yz, g_xy_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_zz_xz_xx[i] = 2.0 * g_xy_zz_xz_xx[i] * a_exp;

        g_x_0_0_0_y_zz_xz_xy[i] = 2.0 * g_xy_zz_xz_xy[i] * a_exp;

        g_x_0_0_0_y_zz_xz_xz[i] = 2.0 * g_xy_zz_xz_xz[i] * a_exp;

        g_x_0_0_0_y_zz_xz_yy[i] = 2.0 * g_xy_zz_xz_yy[i] * a_exp;

        g_x_0_0_0_y_zz_xz_yz[i] = 2.0 * g_xy_zz_xz_yz[i] * a_exp;

        g_x_0_0_0_y_zz_xz_zz[i] = 2.0 * g_xy_zz_xz_zz[i] * a_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_x_0_0_0_y_zz_yy_xx, g_x_0_0_0_y_zz_yy_xy, g_x_0_0_0_y_zz_yy_xz, g_x_0_0_0_y_zz_yy_yy, g_x_0_0_0_y_zz_yy_yz, g_x_0_0_0_y_zz_yy_zz, g_xy_zz_yy_xx, g_xy_zz_yy_xy, g_xy_zz_yy_xz, g_xy_zz_yy_yy, g_xy_zz_yy_yz, g_xy_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_zz_yy_xx[i] = 2.0 * g_xy_zz_yy_xx[i] * a_exp;

        g_x_0_0_0_y_zz_yy_xy[i] = 2.0 * g_xy_zz_yy_xy[i] * a_exp;

        g_x_0_0_0_y_zz_yy_xz[i] = 2.0 * g_xy_zz_yy_xz[i] * a_exp;

        g_x_0_0_0_y_zz_yy_yy[i] = 2.0 * g_xy_zz_yy_yy[i] * a_exp;

        g_x_0_0_0_y_zz_yy_yz[i] = 2.0 * g_xy_zz_yy_yz[i] * a_exp;

        g_x_0_0_0_y_zz_yy_zz[i] = 2.0 * g_xy_zz_yy_zz[i] * a_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_x_0_0_0_y_zz_yz_xx, g_x_0_0_0_y_zz_yz_xy, g_x_0_0_0_y_zz_yz_xz, g_x_0_0_0_y_zz_yz_yy, g_x_0_0_0_y_zz_yz_yz, g_x_0_0_0_y_zz_yz_zz, g_xy_zz_yz_xx, g_xy_zz_yz_xy, g_xy_zz_yz_xz, g_xy_zz_yz_yy, g_xy_zz_yz_yz, g_xy_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_zz_yz_xx[i] = 2.0 * g_xy_zz_yz_xx[i] * a_exp;

        g_x_0_0_0_y_zz_yz_xy[i] = 2.0 * g_xy_zz_yz_xy[i] * a_exp;

        g_x_0_0_0_y_zz_yz_xz[i] = 2.0 * g_xy_zz_yz_xz[i] * a_exp;

        g_x_0_0_0_y_zz_yz_yy[i] = 2.0 * g_xy_zz_yz_yy[i] * a_exp;

        g_x_0_0_0_y_zz_yz_yz[i] = 2.0 * g_xy_zz_yz_yz[i] * a_exp;

        g_x_0_0_0_y_zz_yz_zz[i] = 2.0 * g_xy_zz_yz_zz[i] * a_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_x_0_0_0_y_zz_zz_xx, g_x_0_0_0_y_zz_zz_xy, g_x_0_0_0_y_zz_zz_xz, g_x_0_0_0_y_zz_zz_yy, g_x_0_0_0_y_zz_zz_yz, g_x_0_0_0_y_zz_zz_zz, g_xy_zz_zz_xx, g_xy_zz_zz_xy, g_xy_zz_zz_xz, g_xy_zz_zz_yy, g_xy_zz_zz_yz, g_xy_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_zz_zz_xx[i] = 2.0 * g_xy_zz_zz_xx[i] * a_exp;

        g_x_0_0_0_y_zz_zz_xy[i] = 2.0 * g_xy_zz_zz_xy[i] * a_exp;

        g_x_0_0_0_y_zz_zz_xz[i] = 2.0 * g_xy_zz_zz_xz[i] * a_exp;

        g_x_0_0_0_y_zz_zz_yy[i] = 2.0 * g_xy_zz_zz_yy[i] * a_exp;

        g_x_0_0_0_y_zz_zz_yz[i] = 2.0 * g_xy_zz_zz_yz[i] * a_exp;

        g_x_0_0_0_y_zz_zz_zz[i] = 2.0 * g_xy_zz_zz_zz[i] * a_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_x_0_0_0_z_xx_xx_xx, g_x_0_0_0_z_xx_xx_xy, g_x_0_0_0_z_xx_xx_xz, g_x_0_0_0_z_xx_xx_yy, g_x_0_0_0_z_xx_xx_yz, g_x_0_0_0_z_xx_xx_zz, g_xz_xx_xx_xx, g_xz_xx_xx_xy, g_xz_xx_xx_xz, g_xz_xx_xx_yy, g_xz_xx_xx_yz, g_xz_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xx_xx_xx[i] = 2.0 * g_xz_xx_xx_xx[i] * a_exp;

        g_x_0_0_0_z_xx_xx_xy[i] = 2.0 * g_xz_xx_xx_xy[i] * a_exp;

        g_x_0_0_0_z_xx_xx_xz[i] = 2.0 * g_xz_xx_xx_xz[i] * a_exp;

        g_x_0_0_0_z_xx_xx_yy[i] = 2.0 * g_xz_xx_xx_yy[i] * a_exp;

        g_x_0_0_0_z_xx_xx_yz[i] = 2.0 * g_xz_xx_xx_yz[i] * a_exp;

        g_x_0_0_0_z_xx_xx_zz[i] = 2.0 * g_xz_xx_xx_zz[i] * a_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_x_0_0_0_z_xx_xy_xx, g_x_0_0_0_z_xx_xy_xy, g_x_0_0_0_z_xx_xy_xz, g_x_0_0_0_z_xx_xy_yy, g_x_0_0_0_z_xx_xy_yz, g_x_0_0_0_z_xx_xy_zz, g_xz_xx_xy_xx, g_xz_xx_xy_xy, g_xz_xx_xy_xz, g_xz_xx_xy_yy, g_xz_xx_xy_yz, g_xz_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xx_xy_xx[i] = 2.0 * g_xz_xx_xy_xx[i] * a_exp;

        g_x_0_0_0_z_xx_xy_xy[i] = 2.0 * g_xz_xx_xy_xy[i] * a_exp;

        g_x_0_0_0_z_xx_xy_xz[i] = 2.0 * g_xz_xx_xy_xz[i] * a_exp;

        g_x_0_0_0_z_xx_xy_yy[i] = 2.0 * g_xz_xx_xy_yy[i] * a_exp;

        g_x_0_0_0_z_xx_xy_yz[i] = 2.0 * g_xz_xx_xy_yz[i] * a_exp;

        g_x_0_0_0_z_xx_xy_zz[i] = 2.0 * g_xz_xx_xy_zz[i] * a_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_x_0_0_0_z_xx_xz_xx, g_x_0_0_0_z_xx_xz_xy, g_x_0_0_0_z_xx_xz_xz, g_x_0_0_0_z_xx_xz_yy, g_x_0_0_0_z_xx_xz_yz, g_x_0_0_0_z_xx_xz_zz, g_xz_xx_xz_xx, g_xz_xx_xz_xy, g_xz_xx_xz_xz, g_xz_xx_xz_yy, g_xz_xx_xz_yz, g_xz_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xx_xz_xx[i] = 2.0 * g_xz_xx_xz_xx[i] * a_exp;

        g_x_0_0_0_z_xx_xz_xy[i] = 2.0 * g_xz_xx_xz_xy[i] * a_exp;

        g_x_0_0_0_z_xx_xz_xz[i] = 2.0 * g_xz_xx_xz_xz[i] * a_exp;

        g_x_0_0_0_z_xx_xz_yy[i] = 2.0 * g_xz_xx_xz_yy[i] * a_exp;

        g_x_0_0_0_z_xx_xz_yz[i] = 2.0 * g_xz_xx_xz_yz[i] * a_exp;

        g_x_0_0_0_z_xx_xz_zz[i] = 2.0 * g_xz_xx_xz_zz[i] * a_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_x_0_0_0_z_xx_yy_xx, g_x_0_0_0_z_xx_yy_xy, g_x_0_0_0_z_xx_yy_xz, g_x_0_0_0_z_xx_yy_yy, g_x_0_0_0_z_xx_yy_yz, g_x_0_0_0_z_xx_yy_zz, g_xz_xx_yy_xx, g_xz_xx_yy_xy, g_xz_xx_yy_xz, g_xz_xx_yy_yy, g_xz_xx_yy_yz, g_xz_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xx_yy_xx[i] = 2.0 * g_xz_xx_yy_xx[i] * a_exp;

        g_x_0_0_0_z_xx_yy_xy[i] = 2.0 * g_xz_xx_yy_xy[i] * a_exp;

        g_x_0_0_0_z_xx_yy_xz[i] = 2.0 * g_xz_xx_yy_xz[i] * a_exp;

        g_x_0_0_0_z_xx_yy_yy[i] = 2.0 * g_xz_xx_yy_yy[i] * a_exp;

        g_x_0_0_0_z_xx_yy_yz[i] = 2.0 * g_xz_xx_yy_yz[i] * a_exp;

        g_x_0_0_0_z_xx_yy_zz[i] = 2.0 * g_xz_xx_yy_zz[i] * a_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_x_0_0_0_z_xx_yz_xx, g_x_0_0_0_z_xx_yz_xy, g_x_0_0_0_z_xx_yz_xz, g_x_0_0_0_z_xx_yz_yy, g_x_0_0_0_z_xx_yz_yz, g_x_0_0_0_z_xx_yz_zz, g_xz_xx_yz_xx, g_xz_xx_yz_xy, g_xz_xx_yz_xz, g_xz_xx_yz_yy, g_xz_xx_yz_yz, g_xz_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xx_yz_xx[i] = 2.0 * g_xz_xx_yz_xx[i] * a_exp;

        g_x_0_0_0_z_xx_yz_xy[i] = 2.0 * g_xz_xx_yz_xy[i] * a_exp;

        g_x_0_0_0_z_xx_yz_xz[i] = 2.0 * g_xz_xx_yz_xz[i] * a_exp;

        g_x_0_0_0_z_xx_yz_yy[i] = 2.0 * g_xz_xx_yz_yy[i] * a_exp;

        g_x_0_0_0_z_xx_yz_yz[i] = 2.0 * g_xz_xx_yz_yz[i] * a_exp;

        g_x_0_0_0_z_xx_yz_zz[i] = 2.0 * g_xz_xx_yz_zz[i] * a_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_x_0_0_0_z_xx_zz_xx, g_x_0_0_0_z_xx_zz_xy, g_x_0_0_0_z_xx_zz_xz, g_x_0_0_0_z_xx_zz_yy, g_x_0_0_0_z_xx_zz_yz, g_x_0_0_0_z_xx_zz_zz, g_xz_xx_zz_xx, g_xz_xx_zz_xy, g_xz_xx_zz_xz, g_xz_xx_zz_yy, g_xz_xx_zz_yz, g_xz_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xx_zz_xx[i] = 2.0 * g_xz_xx_zz_xx[i] * a_exp;

        g_x_0_0_0_z_xx_zz_xy[i] = 2.0 * g_xz_xx_zz_xy[i] * a_exp;

        g_x_0_0_0_z_xx_zz_xz[i] = 2.0 * g_xz_xx_zz_xz[i] * a_exp;

        g_x_0_0_0_z_xx_zz_yy[i] = 2.0 * g_xz_xx_zz_yy[i] * a_exp;

        g_x_0_0_0_z_xx_zz_yz[i] = 2.0 * g_xz_xx_zz_yz[i] * a_exp;

        g_x_0_0_0_z_xx_zz_zz[i] = 2.0 * g_xz_xx_zz_zz[i] * a_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_x_0_0_0_z_xy_xx_xx, g_x_0_0_0_z_xy_xx_xy, g_x_0_0_0_z_xy_xx_xz, g_x_0_0_0_z_xy_xx_yy, g_x_0_0_0_z_xy_xx_yz, g_x_0_0_0_z_xy_xx_zz, g_xz_xy_xx_xx, g_xz_xy_xx_xy, g_xz_xy_xx_xz, g_xz_xy_xx_yy, g_xz_xy_xx_yz, g_xz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xy_xx_xx[i] = 2.0 * g_xz_xy_xx_xx[i] * a_exp;

        g_x_0_0_0_z_xy_xx_xy[i] = 2.0 * g_xz_xy_xx_xy[i] * a_exp;

        g_x_0_0_0_z_xy_xx_xz[i] = 2.0 * g_xz_xy_xx_xz[i] * a_exp;

        g_x_0_0_0_z_xy_xx_yy[i] = 2.0 * g_xz_xy_xx_yy[i] * a_exp;

        g_x_0_0_0_z_xy_xx_yz[i] = 2.0 * g_xz_xy_xx_yz[i] * a_exp;

        g_x_0_0_0_z_xy_xx_zz[i] = 2.0 * g_xz_xy_xx_zz[i] * a_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_x_0_0_0_z_xy_xy_xx, g_x_0_0_0_z_xy_xy_xy, g_x_0_0_0_z_xy_xy_xz, g_x_0_0_0_z_xy_xy_yy, g_x_0_0_0_z_xy_xy_yz, g_x_0_0_0_z_xy_xy_zz, g_xz_xy_xy_xx, g_xz_xy_xy_xy, g_xz_xy_xy_xz, g_xz_xy_xy_yy, g_xz_xy_xy_yz, g_xz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xy_xy_xx[i] = 2.0 * g_xz_xy_xy_xx[i] * a_exp;

        g_x_0_0_0_z_xy_xy_xy[i] = 2.0 * g_xz_xy_xy_xy[i] * a_exp;

        g_x_0_0_0_z_xy_xy_xz[i] = 2.0 * g_xz_xy_xy_xz[i] * a_exp;

        g_x_0_0_0_z_xy_xy_yy[i] = 2.0 * g_xz_xy_xy_yy[i] * a_exp;

        g_x_0_0_0_z_xy_xy_yz[i] = 2.0 * g_xz_xy_xy_yz[i] * a_exp;

        g_x_0_0_0_z_xy_xy_zz[i] = 2.0 * g_xz_xy_xy_zz[i] * a_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_x_0_0_0_z_xy_xz_xx, g_x_0_0_0_z_xy_xz_xy, g_x_0_0_0_z_xy_xz_xz, g_x_0_0_0_z_xy_xz_yy, g_x_0_0_0_z_xy_xz_yz, g_x_0_0_0_z_xy_xz_zz, g_xz_xy_xz_xx, g_xz_xy_xz_xy, g_xz_xy_xz_xz, g_xz_xy_xz_yy, g_xz_xy_xz_yz, g_xz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xy_xz_xx[i] = 2.0 * g_xz_xy_xz_xx[i] * a_exp;

        g_x_0_0_0_z_xy_xz_xy[i] = 2.0 * g_xz_xy_xz_xy[i] * a_exp;

        g_x_0_0_0_z_xy_xz_xz[i] = 2.0 * g_xz_xy_xz_xz[i] * a_exp;

        g_x_0_0_0_z_xy_xz_yy[i] = 2.0 * g_xz_xy_xz_yy[i] * a_exp;

        g_x_0_0_0_z_xy_xz_yz[i] = 2.0 * g_xz_xy_xz_yz[i] * a_exp;

        g_x_0_0_0_z_xy_xz_zz[i] = 2.0 * g_xz_xy_xz_zz[i] * a_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_x_0_0_0_z_xy_yy_xx, g_x_0_0_0_z_xy_yy_xy, g_x_0_0_0_z_xy_yy_xz, g_x_0_0_0_z_xy_yy_yy, g_x_0_0_0_z_xy_yy_yz, g_x_0_0_0_z_xy_yy_zz, g_xz_xy_yy_xx, g_xz_xy_yy_xy, g_xz_xy_yy_xz, g_xz_xy_yy_yy, g_xz_xy_yy_yz, g_xz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xy_yy_xx[i] = 2.0 * g_xz_xy_yy_xx[i] * a_exp;

        g_x_0_0_0_z_xy_yy_xy[i] = 2.0 * g_xz_xy_yy_xy[i] * a_exp;

        g_x_0_0_0_z_xy_yy_xz[i] = 2.0 * g_xz_xy_yy_xz[i] * a_exp;

        g_x_0_0_0_z_xy_yy_yy[i] = 2.0 * g_xz_xy_yy_yy[i] * a_exp;

        g_x_0_0_0_z_xy_yy_yz[i] = 2.0 * g_xz_xy_yy_yz[i] * a_exp;

        g_x_0_0_0_z_xy_yy_zz[i] = 2.0 * g_xz_xy_yy_zz[i] * a_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_x_0_0_0_z_xy_yz_xx, g_x_0_0_0_z_xy_yz_xy, g_x_0_0_0_z_xy_yz_xz, g_x_0_0_0_z_xy_yz_yy, g_x_0_0_0_z_xy_yz_yz, g_x_0_0_0_z_xy_yz_zz, g_xz_xy_yz_xx, g_xz_xy_yz_xy, g_xz_xy_yz_xz, g_xz_xy_yz_yy, g_xz_xy_yz_yz, g_xz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xy_yz_xx[i] = 2.0 * g_xz_xy_yz_xx[i] * a_exp;

        g_x_0_0_0_z_xy_yz_xy[i] = 2.0 * g_xz_xy_yz_xy[i] * a_exp;

        g_x_0_0_0_z_xy_yz_xz[i] = 2.0 * g_xz_xy_yz_xz[i] * a_exp;

        g_x_0_0_0_z_xy_yz_yy[i] = 2.0 * g_xz_xy_yz_yy[i] * a_exp;

        g_x_0_0_0_z_xy_yz_yz[i] = 2.0 * g_xz_xy_yz_yz[i] * a_exp;

        g_x_0_0_0_z_xy_yz_zz[i] = 2.0 * g_xz_xy_yz_zz[i] * a_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_x_0_0_0_z_xy_zz_xx, g_x_0_0_0_z_xy_zz_xy, g_x_0_0_0_z_xy_zz_xz, g_x_0_0_0_z_xy_zz_yy, g_x_0_0_0_z_xy_zz_yz, g_x_0_0_0_z_xy_zz_zz, g_xz_xy_zz_xx, g_xz_xy_zz_xy, g_xz_xy_zz_xz, g_xz_xy_zz_yy, g_xz_xy_zz_yz, g_xz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xy_zz_xx[i] = 2.0 * g_xz_xy_zz_xx[i] * a_exp;

        g_x_0_0_0_z_xy_zz_xy[i] = 2.0 * g_xz_xy_zz_xy[i] * a_exp;

        g_x_0_0_0_z_xy_zz_xz[i] = 2.0 * g_xz_xy_zz_xz[i] * a_exp;

        g_x_0_0_0_z_xy_zz_yy[i] = 2.0 * g_xz_xy_zz_yy[i] * a_exp;

        g_x_0_0_0_z_xy_zz_yz[i] = 2.0 * g_xz_xy_zz_yz[i] * a_exp;

        g_x_0_0_0_z_xy_zz_zz[i] = 2.0 * g_xz_xy_zz_zz[i] * a_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_x_0_0_0_z_xz_xx_xx, g_x_0_0_0_z_xz_xx_xy, g_x_0_0_0_z_xz_xx_xz, g_x_0_0_0_z_xz_xx_yy, g_x_0_0_0_z_xz_xx_yz, g_x_0_0_0_z_xz_xx_zz, g_xz_xz_xx_xx, g_xz_xz_xx_xy, g_xz_xz_xx_xz, g_xz_xz_xx_yy, g_xz_xz_xx_yz, g_xz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xz_xx_xx[i] = 2.0 * g_xz_xz_xx_xx[i] * a_exp;

        g_x_0_0_0_z_xz_xx_xy[i] = 2.0 * g_xz_xz_xx_xy[i] * a_exp;

        g_x_0_0_0_z_xz_xx_xz[i] = 2.0 * g_xz_xz_xx_xz[i] * a_exp;

        g_x_0_0_0_z_xz_xx_yy[i] = 2.0 * g_xz_xz_xx_yy[i] * a_exp;

        g_x_0_0_0_z_xz_xx_yz[i] = 2.0 * g_xz_xz_xx_yz[i] * a_exp;

        g_x_0_0_0_z_xz_xx_zz[i] = 2.0 * g_xz_xz_xx_zz[i] * a_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_x_0_0_0_z_xz_xy_xx, g_x_0_0_0_z_xz_xy_xy, g_x_0_0_0_z_xz_xy_xz, g_x_0_0_0_z_xz_xy_yy, g_x_0_0_0_z_xz_xy_yz, g_x_0_0_0_z_xz_xy_zz, g_xz_xz_xy_xx, g_xz_xz_xy_xy, g_xz_xz_xy_xz, g_xz_xz_xy_yy, g_xz_xz_xy_yz, g_xz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xz_xy_xx[i] = 2.0 * g_xz_xz_xy_xx[i] * a_exp;

        g_x_0_0_0_z_xz_xy_xy[i] = 2.0 * g_xz_xz_xy_xy[i] * a_exp;

        g_x_0_0_0_z_xz_xy_xz[i] = 2.0 * g_xz_xz_xy_xz[i] * a_exp;

        g_x_0_0_0_z_xz_xy_yy[i] = 2.0 * g_xz_xz_xy_yy[i] * a_exp;

        g_x_0_0_0_z_xz_xy_yz[i] = 2.0 * g_xz_xz_xy_yz[i] * a_exp;

        g_x_0_0_0_z_xz_xy_zz[i] = 2.0 * g_xz_xz_xy_zz[i] * a_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_x_0_0_0_z_xz_xz_xx, g_x_0_0_0_z_xz_xz_xy, g_x_0_0_0_z_xz_xz_xz, g_x_0_0_0_z_xz_xz_yy, g_x_0_0_0_z_xz_xz_yz, g_x_0_0_0_z_xz_xz_zz, g_xz_xz_xz_xx, g_xz_xz_xz_xy, g_xz_xz_xz_xz, g_xz_xz_xz_yy, g_xz_xz_xz_yz, g_xz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xz_xz_xx[i] = 2.0 * g_xz_xz_xz_xx[i] * a_exp;

        g_x_0_0_0_z_xz_xz_xy[i] = 2.0 * g_xz_xz_xz_xy[i] * a_exp;

        g_x_0_0_0_z_xz_xz_xz[i] = 2.0 * g_xz_xz_xz_xz[i] * a_exp;

        g_x_0_0_0_z_xz_xz_yy[i] = 2.0 * g_xz_xz_xz_yy[i] * a_exp;

        g_x_0_0_0_z_xz_xz_yz[i] = 2.0 * g_xz_xz_xz_yz[i] * a_exp;

        g_x_0_0_0_z_xz_xz_zz[i] = 2.0 * g_xz_xz_xz_zz[i] * a_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_x_0_0_0_z_xz_yy_xx, g_x_0_0_0_z_xz_yy_xy, g_x_0_0_0_z_xz_yy_xz, g_x_0_0_0_z_xz_yy_yy, g_x_0_0_0_z_xz_yy_yz, g_x_0_0_0_z_xz_yy_zz, g_xz_xz_yy_xx, g_xz_xz_yy_xy, g_xz_xz_yy_xz, g_xz_xz_yy_yy, g_xz_xz_yy_yz, g_xz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xz_yy_xx[i] = 2.0 * g_xz_xz_yy_xx[i] * a_exp;

        g_x_0_0_0_z_xz_yy_xy[i] = 2.0 * g_xz_xz_yy_xy[i] * a_exp;

        g_x_0_0_0_z_xz_yy_xz[i] = 2.0 * g_xz_xz_yy_xz[i] * a_exp;

        g_x_0_0_0_z_xz_yy_yy[i] = 2.0 * g_xz_xz_yy_yy[i] * a_exp;

        g_x_0_0_0_z_xz_yy_yz[i] = 2.0 * g_xz_xz_yy_yz[i] * a_exp;

        g_x_0_0_0_z_xz_yy_zz[i] = 2.0 * g_xz_xz_yy_zz[i] * a_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_x_0_0_0_z_xz_yz_xx, g_x_0_0_0_z_xz_yz_xy, g_x_0_0_0_z_xz_yz_xz, g_x_0_0_0_z_xz_yz_yy, g_x_0_0_0_z_xz_yz_yz, g_x_0_0_0_z_xz_yz_zz, g_xz_xz_yz_xx, g_xz_xz_yz_xy, g_xz_xz_yz_xz, g_xz_xz_yz_yy, g_xz_xz_yz_yz, g_xz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xz_yz_xx[i] = 2.0 * g_xz_xz_yz_xx[i] * a_exp;

        g_x_0_0_0_z_xz_yz_xy[i] = 2.0 * g_xz_xz_yz_xy[i] * a_exp;

        g_x_0_0_0_z_xz_yz_xz[i] = 2.0 * g_xz_xz_yz_xz[i] * a_exp;

        g_x_0_0_0_z_xz_yz_yy[i] = 2.0 * g_xz_xz_yz_yy[i] * a_exp;

        g_x_0_0_0_z_xz_yz_yz[i] = 2.0 * g_xz_xz_yz_yz[i] * a_exp;

        g_x_0_0_0_z_xz_yz_zz[i] = 2.0 * g_xz_xz_yz_zz[i] * a_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_x_0_0_0_z_xz_zz_xx, g_x_0_0_0_z_xz_zz_xy, g_x_0_0_0_z_xz_zz_xz, g_x_0_0_0_z_xz_zz_yy, g_x_0_0_0_z_xz_zz_yz, g_x_0_0_0_z_xz_zz_zz, g_xz_xz_zz_xx, g_xz_xz_zz_xy, g_xz_xz_zz_xz, g_xz_xz_zz_yy, g_xz_xz_zz_yz, g_xz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xz_zz_xx[i] = 2.0 * g_xz_xz_zz_xx[i] * a_exp;

        g_x_0_0_0_z_xz_zz_xy[i] = 2.0 * g_xz_xz_zz_xy[i] * a_exp;

        g_x_0_0_0_z_xz_zz_xz[i] = 2.0 * g_xz_xz_zz_xz[i] * a_exp;

        g_x_0_0_0_z_xz_zz_yy[i] = 2.0 * g_xz_xz_zz_yy[i] * a_exp;

        g_x_0_0_0_z_xz_zz_yz[i] = 2.0 * g_xz_xz_zz_yz[i] * a_exp;

        g_x_0_0_0_z_xz_zz_zz[i] = 2.0 * g_xz_xz_zz_zz[i] * a_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_x_0_0_0_z_yy_xx_xx, g_x_0_0_0_z_yy_xx_xy, g_x_0_0_0_z_yy_xx_xz, g_x_0_0_0_z_yy_xx_yy, g_x_0_0_0_z_yy_xx_yz, g_x_0_0_0_z_yy_xx_zz, g_xz_yy_xx_xx, g_xz_yy_xx_xy, g_xz_yy_xx_xz, g_xz_yy_xx_yy, g_xz_yy_xx_yz, g_xz_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yy_xx_xx[i] = 2.0 * g_xz_yy_xx_xx[i] * a_exp;

        g_x_0_0_0_z_yy_xx_xy[i] = 2.0 * g_xz_yy_xx_xy[i] * a_exp;

        g_x_0_0_0_z_yy_xx_xz[i] = 2.0 * g_xz_yy_xx_xz[i] * a_exp;

        g_x_0_0_0_z_yy_xx_yy[i] = 2.0 * g_xz_yy_xx_yy[i] * a_exp;

        g_x_0_0_0_z_yy_xx_yz[i] = 2.0 * g_xz_yy_xx_yz[i] * a_exp;

        g_x_0_0_0_z_yy_xx_zz[i] = 2.0 * g_xz_yy_xx_zz[i] * a_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_x_0_0_0_z_yy_xy_xx, g_x_0_0_0_z_yy_xy_xy, g_x_0_0_0_z_yy_xy_xz, g_x_0_0_0_z_yy_xy_yy, g_x_0_0_0_z_yy_xy_yz, g_x_0_0_0_z_yy_xy_zz, g_xz_yy_xy_xx, g_xz_yy_xy_xy, g_xz_yy_xy_xz, g_xz_yy_xy_yy, g_xz_yy_xy_yz, g_xz_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yy_xy_xx[i] = 2.0 * g_xz_yy_xy_xx[i] * a_exp;

        g_x_0_0_0_z_yy_xy_xy[i] = 2.0 * g_xz_yy_xy_xy[i] * a_exp;

        g_x_0_0_0_z_yy_xy_xz[i] = 2.0 * g_xz_yy_xy_xz[i] * a_exp;

        g_x_0_0_0_z_yy_xy_yy[i] = 2.0 * g_xz_yy_xy_yy[i] * a_exp;

        g_x_0_0_0_z_yy_xy_yz[i] = 2.0 * g_xz_yy_xy_yz[i] * a_exp;

        g_x_0_0_0_z_yy_xy_zz[i] = 2.0 * g_xz_yy_xy_zz[i] * a_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_x_0_0_0_z_yy_xz_xx, g_x_0_0_0_z_yy_xz_xy, g_x_0_0_0_z_yy_xz_xz, g_x_0_0_0_z_yy_xz_yy, g_x_0_0_0_z_yy_xz_yz, g_x_0_0_0_z_yy_xz_zz, g_xz_yy_xz_xx, g_xz_yy_xz_xy, g_xz_yy_xz_xz, g_xz_yy_xz_yy, g_xz_yy_xz_yz, g_xz_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yy_xz_xx[i] = 2.0 * g_xz_yy_xz_xx[i] * a_exp;

        g_x_0_0_0_z_yy_xz_xy[i] = 2.0 * g_xz_yy_xz_xy[i] * a_exp;

        g_x_0_0_0_z_yy_xz_xz[i] = 2.0 * g_xz_yy_xz_xz[i] * a_exp;

        g_x_0_0_0_z_yy_xz_yy[i] = 2.0 * g_xz_yy_xz_yy[i] * a_exp;

        g_x_0_0_0_z_yy_xz_yz[i] = 2.0 * g_xz_yy_xz_yz[i] * a_exp;

        g_x_0_0_0_z_yy_xz_zz[i] = 2.0 * g_xz_yy_xz_zz[i] * a_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_x_0_0_0_z_yy_yy_xx, g_x_0_0_0_z_yy_yy_xy, g_x_0_0_0_z_yy_yy_xz, g_x_0_0_0_z_yy_yy_yy, g_x_0_0_0_z_yy_yy_yz, g_x_0_0_0_z_yy_yy_zz, g_xz_yy_yy_xx, g_xz_yy_yy_xy, g_xz_yy_yy_xz, g_xz_yy_yy_yy, g_xz_yy_yy_yz, g_xz_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yy_yy_xx[i] = 2.0 * g_xz_yy_yy_xx[i] * a_exp;

        g_x_0_0_0_z_yy_yy_xy[i] = 2.0 * g_xz_yy_yy_xy[i] * a_exp;

        g_x_0_0_0_z_yy_yy_xz[i] = 2.0 * g_xz_yy_yy_xz[i] * a_exp;

        g_x_0_0_0_z_yy_yy_yy[i] = 2.0 * g_xz_yy_yy_yy[i] * a_exp;

        g_x_0_0_0_z_yy_yy_yz[i] = 2.0 * g_xz_yy_yy_yz[i] * a_exp;

        g_x_0_0_0_z_yy_yy_zz[i] = 2.0 * g_xz_yy_yy_zz[i] * a_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_x_0_0_0_z_yy_yz_xx, g_x_0_0_0_z_yy_yz_xy, g_x_0_0_0_z_yy_yz_xz, g_x_0_0_0_z_yy_yz_yy, g_x_0_0_0_z_yy_yz_yz, g_x_0_0_0_z_yy_yz_zz, g_xz_yy_yz_xx, g_xz_yy_yz_xy, g_xz_yy_yz_xz, g_xz_yy_yz_yy, g_xz_yy_yz_yz, g_xz_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yy_yz_xx[i] = 2.0 * g_xz_yy_yz_xx[i] * a_exp;

        g_x_0_0_0_z_yy_yz_xy[i] = 2.0 * g_xz_yy_yz_xy[i] * a_exp;

        g_x_0_0_0_z_yy_yz_xz[i] = 2.0 * g_xz_yy_yz_xz[i] * a_exp;

        g_x_0_0_0_z_yy_yz_yy[i] = 2.0 * g_xz_yy_yz_yy[i] * a_exp;

        g_x_0_0_0_z_yy_yz_yz[i] = 2.0 * g_xz_yy_yz_yz[i] * a_exp;

        g_x_0_0_0_z_yy_yz_zz[i] = 2.0 * g_xz_yy_yz_zz[i] * a_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_x_0_0_0_z_yy_zz_xx, g_x_0_0_0_z_yy_zz_xy, g_x_0_0_0_z_yy_zz_xz, g_x_0_0_0_z_yy_zz_yy, g_x_0_0_0_z_yy_zz_yz, g_x_0_0_0_z_yy_zz_zz, g_xz_yy_zz_xx, g_xz_yy_zz_xy, g_xz_yy_zz_xz, g_xz_yy_zz_yy, g_xz_yy_zz_yz, g_xz_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yy_zz_xx[i] = 2.0 * g_xz_yy_zz_xx[i] * a_exp;

        g_x_0_0_0_z_yy_zz_xy[i] = 2.0 * g_xz_yy_zz_xy[i] * a_exp;

        g_x_0_0_0_z_yy_zz_xz[i] = 2.0 * g_xz_yy_zz_xz[i] * a_exp;

        g_x_0_0_0_z_yy_zz_yy[i] = 2.0 * g_xz_yy_zz_yy[i] * a_exp;

        g_x_0_0_0_z_yy_zz_yz[i] = 2.0 * g_xz_yy_zz_yz[i] * a_exp;

        g_x_0_0_0_z_yy_zz_zz[i] = 2.0 * g_xz_yy_zz_zz[i] * a_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_x_0_0_0_z_yz_xx_xx, g_x_0_0_0_z_yz_xx_xy, g_x_0_0_0_z_yz_xx_xz, g_x_0_0_0_z_yz_xx_yy, g_x_0_0_0_z_yz_xx_yz, g_x_0_0_0_z_yz_xx_zz, g_xz_yz_xx_xx, g_xz_yz_xx_xy, g_xz_yz_xx_xz, g_xz_yz_xx_yy, g_xz_yz_xx_yz, g_xz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yz_xx_xx[i] = 2.0 * g_xz_yz_xx_xx[i] * a_exp;

        g_x_0_0_0_z_yz_xx_xy[i] = 2.0 * g_xz_yz_xx_xy[i] * a_exp;

        g_x_0_0_0_z_yz_xx_xz[i] = 2.0 * g_xz_yz_xx_xz[i] * a_exp;

        g_x_0_0_0_z_yz_xx_yy[i] = 2.0 * g_xz_yz_xx_yy[i] * a_exp;

        g_x_0_0_0_z_yz_xx_yz[i] = 2.0 * g_xz_yz_xx_yz[i] * a_exp;

        g_x_0_0_0_z_yz_xx_zz[i] = 2.0 * g_xz_yz_xx_zz[i] * a_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_x_0_0_0_z_yz_xy_xx, g_x_0_0_0_z_yz_xy_xy, g_x_0_0_0_z_yz_xy_xz, g_x_0_0_0_z_yz_xy_yy, g_x_0_0_0_z_yz_xy_yz, g_x_0_0_0_z_yz_xy_zz, g_xz_yz_xy_xx, g_xz_yz_xy_xy, g_xz_yz_xy_xz, g_xz_yz_xy_yy, g_xz_yz_xy_yz, g_xz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yz_xy_xx[i] = 2.0 * g_xz_yz_xy_xx[i] * a_exp;

        g_x_0_0_0_z_yz_xy_xy[i] = 2.0 * g_xz_yz_xy_xy[i] * a_exp;

        g_x_0_0_0_z_yz_xy_xz[i] = 2.0 * g_xz_yz_xy_xz[i] * a_exp;

        g_x_0_0_0_z_yz_xy_yy[i] = 2.0 * g_xz_yz_xy_yy[i] * a_exp;

        g_x_0_0_0_z_yz_xy_yz[i] = 2.0 * g_xz_yz_xy_yz[i] * a_exp;

        g_x_0_0_0_z_yz_xy_zz[i] = 2.0 * g_xz_yz_xy_zz[i] * a_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_x_0_0_0_z_yz_xz_xx, g_x_0_0_0_z_yz_xz_xy, g_x_0_0_0_z_yz_xz_xz, g_x_0_0_0_z_yz_xz_yy, g_x_0_0_0_z_yz_xz_yz, g_x_0_0_0_z_yz_xz_zz, g_xz_yz_xz_xx, g_xz_yz_xz_xy, g_xz_yz_xz_xz, g_xz_yz_xz_yy, g_xz_yz_xz_yz, g_xz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yz_xz_xx[i] = 2.0 * g_xz_yz_xz_xx[i] * a_exp;

        g_x_0_0_0_z_yz_xz_xy[i] = 2.0 * g_xz_yz_xz_xy[i] * a_exp;

        g_x_0_0_0_z_yz_xz_xz[i] = 2.0 * g_xz_yz_xz_xz[i] * a_exp;

        g_x_0_0_0_z_yz_xz_yy[i] = 2.0 * g_xz_yz_xz_yy[i] * a_exp;

        g_x_0_0_0_z_yz_xz_yz[i] = 2.0 * g_xz_yz_xz_yz[i] * a_exp;

        g_x_0_0_0_z_yz_xz_zz[i] = 2.0 * g_xz_yz_xz_zz[i] * a_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_x_0_0_0_z_yz_yy_xx, g_x_0_0_0_z_yz_yy_xy, g_x_0_0_0_z_yz_yy_xz, g_x_0_0_0_z_yz_yy_yy, g_x_0_0_0_z_yz_yy_yz, g_x_0_0_0_z_yz_yy_zz, g_xz_yz_yy_xx, g_xz_yz_yy_xy, g_xz_yz_yy_xz, g_xz_yz_yy_yy, g_xz_yz_yy_yz, g_xz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yz_yy_xx[i] = 2.0 * g_xz_yz_yy_xx[i] * a_exp;

        g_x_0_0_0_z_yz_yy_xy[i] = 2.0 * g_xz_yz_yy_xy[i] * a_exp;

        g_x_0_0_0_z_yz_yy_xz[i] = 2.0 * g_xz_yz_yy_xz[i] * a_exp;

        g_x_0_0_0_z_yz_yy_yy[i] = 2.0 * g_xz_yz_yy_yy[i] * a_exp;

        g_x_0_0_0_z_yz_yy_yz[i] = 2.0 * g_xz_yz_yy_yz[i] * a_exp;

        g_x_0_0_0_z_yz_yy_zz[i] = 2.0 * g_xz_yz_yy_zz[i] * a_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_x_0_0_0_z_yz_yz_xx, g_x_0_0_0_z_yz_yz_xy, g_x_0_0_0_z_yz_yz_xz, g_x_0_0_0_z_yz_yz_yy, g_x_0_0_0_z_yz_yz_yz, g_x_0_0_0_z_yz_yz_zz, g_xz_yz_yz_xx, g_xz_yz_yz_xy, g_xz_yz_yz_xz, g_xz_yz_yz_yy, g_xz_yz_yz_yz, g_xz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yz_yz_xx[i] = 2.0 * g_xz_yz_yz_xx[i] * a_exp;

        g_x_0_0_0_z_yz_yz_xy[i] = 2.0 * g_xz_yz_yz_xy[i] * a_exp;

        g_x_0_0_0_z_yz_yz_xz[i] = 2.0 * g_xz_yz_yz_xz[i] * a_exp;

        g_x_0_0_0_z_yz_yz_yy[i] = 2.0 * g_xz_yz_yz_yy[i] * a_exp;

        g_x_0_0_0_z_yz_yz_yz[i] = 2.0 * g_xz_yz_yz_yz[i] * a_exp;

        g_x_0_0_0_z_yz_yz_zz[i] = 2.0 * g_xz_yz_yz_zz[i] * a_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_x_0_0_0_z_yz_zz_xx, g_x_0_0_0_z_yz_zz_xy, g_x_0_0_0_z_yz_zz_xz, g_x_0_0_0_z_yz_zz_yy, g_x_0_0_0_z_yz_zz_yz, g_x_0_0_0_z_yz_zz_zz, g_xz_yz_zz_xx, g_xz_yz_zz_xy, g_xz_yz_zz_xz, g_xz_yz_zz_yy, g_xz_yz_zz_yz, g_xz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yz_zz_xx[i] = 2.0 * g_xz_yz_zz_xx[i] * a_exp;

        g_x_0_0_0_z_yz_zz_xy[i] = 2.0 * g_xz_yz_zz_xy[i] * a_exp;

        g_x_0_0_0_z_yz_zz_xz[i] = 2.0 * g_xz_yz_zz_xz[i] * a_exp;

        g_x_0_0_0_z_yz_zz_yy[i] = 2.0 * g_xz_yz_zz_yy[i] * a_exp;

        g_x_0_0_0_z_yz_zz_yz[i] = 2.0 * g_xz_yz_zz_yz[i] * a_exp;

        g_x_0_0_0_z_yz_zz_zz[i] = 2.0 * g_xz_yz_zz_zz[i] * a_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_x_0_0_0_z_zz_xx_xx, g_x_0_0_0_z_zz_xx_xy, g_x_0_0_0_z_zz_xx_xz, g_x_0_0_0_z_zz_xx_yy, g_x_0_0_0_z_zz_xx_yz, g_x_0_0_0_z_zz_xx_zz, g_xz_zz_xx_xx, g_xz_zz_xx_xy, g_xz_zz_xx_xz, g_xz_zz_xx_yy, g_xz_zz_xx_yz, g_xz_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_zz_xx_xx[i] = 2.0 * g_xz_zz_xx_xx[i] * a_exp;

        g_x_0_0_0_z_zz_xx_xy[i] = 2.0 * g_xz_zz_xx_xy[i] * a_exp;

        g_x_0_0_0_z_zz_xx_xz[i] = 2.0 * g_xz_zz_xx_xz[i] * a_exp;

        g_x_0_0_0_z_zz_xx_yy[i] = 2.0 * g_xz_zz_xx_yy[i] * a_exp;

        g_x_0_0_0_z_zz_xx_yz[i] = 2.0 * g_xz_zz_xx_yz[i] * a_exp;

        g_x_0_0_0_z_zz_xx_zz[i] = 2.0 * g_xz_zz_xx_zz[i] * a_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_x_0_0_0_z_zz_xy_xx, g_x_0_0_0_z_zz_xy_xy, g_x_0_0_0_z_zz_xy_xz, g_x_0_0_0_z_zz_xy_yy, g_x_0_0_0_z_zz_xy_yz, g_x_0_0_0_z_zz_xy_zz, g_xz_zz_xy_xx, g_xz_zz_xy_xy, g_xz_zz_xy_xz, g_xz_zz_xy_yy, g_xz_zz_xy_yz, g_xz_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_zz_xy_xx[i] = 2.0 * g_xz_zz_xy_xx[i] * a_exp;

        g_x_0_0_0_z_zz_xy_xy[i] = 2.0 * g_xz_zz_xy_xy[i] * a_exp;

        g_x_0_0_0_z_zz_xy_xz[i] = 2.0 * g_xz_zz_xy_xz[i] * a_exp;

        g_x_0_0_0_z_zz_xy_yy[i] = 2.0 * g_xz_zz_xy_yy[i] * a_exp;

        g_x_0_0_0_z_zz_xy_yz[i] = 2.0 * g_xz_zz_xy_yz[i] * a_exp;

        g_x_0_0_0_z_zz_xy_zz[i] = 2.0 * g_xz_zz_xy_zz[i] * a_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_x_0_0_0_z_zz_xz_xx, g_x_0_0_0_z_zz_xz_xy, g_x_0_0_0_z_zz_xz_xz, g_x_0_0_0_z_zz_xz_yy, g_x_0_0_0_z_zz_xz_yz, g_x_0_0_0_z_zz_xz_zz, g_xz_zz_xz_xx, g_xz_zz_xz_xy, g_xz_zz_xz_xz, g_xz_zz_xz_yy, g_xz_zz_xz_yz, g_xz_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_zz_xz_xx[i] = 2.0 * g_xz_zz_xz_xx[i] * a_exp;

        g_x_0_0_0_z_zz_xz_xy[i] = 2.0 * g_xz_zz_xz_xy[i] * a_exp;

        g_x_0_0_0_z_zz_xz_xz[i] = 2.0 * g_xz_zz_xz_xz[i] * a_exp;

        g_x_0_0_0_z_zz_xz_yy[i] = 2.0 * g_xz_zz_xz_yy[i] * a_exp;

        g_x_0_0_0_z_zz_xz_yz[i] = 2.0 * g_xz_zz_xz_yz[i] * a_exp;

        g_x_0_0_0_z_zz_xz_zz[i] = 2.0 * g_xz_zz_xz_zz[i] * a_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_x_0_0_0_z_zz_yy_xx, g_x_0_0_0_z_zz_yy_xy, g_x_0_0_0_z_zz_yy_xz, g_x_0_0_0_z_zz_yy_yy, g_x_0_0_0_z_zz_yy_yz, g_x_0_0_0_z_zz_yy_zz, g_xz_zz_yy_xx, g_xz_zz_yy_xy, g_xz_zz_yy_xz, g_xz_zz_yy_yy, g_xz_zz_yy_yz, g_xz_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_zz_yy_xx[i] = 2.0 * g_xz_zz_yy_xx[i] * a_exp;

        g_x_0_0_0_z_zz_yy_xy[i] = 2.0 * g_xz_zz_yy_xy[i] * a_exp;

        g_x_0_0_0_z_zz_yy_xz[i] = 2.0 * g_xz_zz_yy_xz[i] * a_exp;

        g_x_0_0_0_z_zz_yy_yy[i] = 2.0 * g_xz_zz_yy_yy[i] * a_exp;

        g_x_0_0_0_z_zz_yy_yz[i] = 2.0 * g_xz_zz_yy_yz[i] * a_exp;

        g_x_0_0_0_z_zz_yy_zz[i] = 2.0 * g_xz_zz_yy_zz[i] * a_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_x_0_0_0_z_zz_yz_xx, g_x_0_0_0_z_zz_yz_xy, g_x_0_0_0_z_zz_yz_xz, g_x_0_0_0_z_zz_yz_yy, g_x_0_0_0_z_zz_yz_yz, g_x_0_0_0_z_zz_yz_zz, g_xz_zz_yz_xx, g_xz_zz_yz_xy, g_xz_zz_yz_xz, g_xz_zz_yz_yy, g_xz_zz_yz_yz, g_xz_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_zz_yz_xx[i] = 2.0 * g_xz_zz_yz_xx[i] * a_exp;

        g_x_0_0_0_z_zz_yz_xy[i] = 2.0 * g_xz_zz_yz_xy[i] * a_exp;

        g_x_0_0_0_z_zz_yz_xz[i] = 2.0 * g_xz_zz_yz_xz[i] * a_exp;

        g_x_0_0_0_z_zz_yz_yy[i] = 2.0 * g_xz_zz_yz_yy[i] * a_exp;

        g_x_0_0_0_z_zz_yz_yz[i] = 2.0 * g_xz_zz_yz_yz[i] * a_exp;

        g_x_0_0_0_z_zz_yz_zz[i] = 2.0 * g_xz_zz_yz_zz[i] * a_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_x_0_0_0_z_zz_zz_xx, g_x_0_0_0_z_zz_zz_xy, g_x_0_0_0_z_zz_zz_xz, g_x_0_0_0_z_zz_zz_yy, g_x_0_0_0_z_zz_zz_yz, g_x_0_0_0_z_zz_zz_zz, g_xz_zz_zz_xx, g_xz_zz_zz_xy, g_xz_zz_zz_xz, g_xz_zz_zz_yy, g_xz_zz_zz_yz, g_xz_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_zz_zz_xx[i] = 2.0 * g_xz_zz_zz_xx[i] * a_exp;

        g_x_0_0_0_z_zz_zz_xy[i] = 2.0 * g_xz_zz_zz_xy[i] * a_exp;

        g_x_0_0_0_z_zz_zz_xz[i] = 2.0 * g_xz_zz_zz_xz[i] * a_exp;

        g_x_0_0_0_z_zz_zz_yy[i] = 2.0 * g_xz_zz_zz_yy[i] * a_exp;

        g_x_0_0_0_z_zz_zz_yz[i] = 2.0 * g_xz_zz_zz_yz[i] * a_exp;

        g_x_0_0_0_z_zz_zz_zz[i] = 2.0 * g_xz_zz_zz_zz[i] * a_exp;
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_xy_xx_xx_xx, g_xy_xx_xx_xy, g_xy_xx_xx_xz, g_xy_xx_xx_yy, g_xy_xx_xx_yz, g_xy_xx_xx_zz, g_y_0_0_0_x_xx_xx_xx, g_y_0_0_0_x_xx_xx_xy, g_y_0_0_0_x_xx_xx_xz, g_y_0_0_0_x_xx_xx_yy, g_y_0_0_0_x_xx_xx_yz, g_y_0_0_0_x_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xx_xx_xx[i] = 2.0 * g_xy_xx_xx_xx[i] * a_exp;

        g_y_0_0_0_x_xx_xx_xy[i] = 2.0 * g_xy_xx_xx_xy[i] * a_exp;

        g_y_0_0_0_x_xx_xx_xz[i] = 2.0 * g_xy_xx_xx_xz[i] * a_exp;

        g_y_0_0_0_x_xx_xx_yy[i] = 2.0 * g_xy_xx_xx_yy[i] * a_exp;

        g_y_0_0_0_x_xx_xx_yz[i] = 2.0 * g_xy_xx_xx_yz[i] * a_exp;

        g_y_0_0_0_x_xx_xx_zz[i] = 2.0 * g_xy_xx_xx_zz[i] * a_exp;
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_xy_xx_xy_xx, g_xy_xx_xy_xy, g_xy_xx_xy_xz, g_xy_xx_xy_yy, g_xy_xx_xy_yz, g_xy_xx_xy_zz, g_y_0_0_0_x_xx_xy_xx, g_y_0_0_0_x_xx_xy_xy, g_y_0_0_0_x_xx_xy_xz, g_y_0_0_0_x_xx_xy_yy, g_y_0_0_0_x_xx_xy_yz, g_y_0_0_0_x_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xx_xy_xx[i] = 2.0 * g_xy_xx_xy_xx[i] * a_exp;

        g_y_0_0_0_x_xx_xy_xy[i] = 2.0 * g_xy_xx_xy_xy[i] * a_exp;

        g_y_0_0_0_x_xx_xy_xz[i] = 2.0 * g_xy_xx_xy_xz[i] * a_exp;

        g_y_0_0_0_x_xx_xy_yy[i] = 2.0 * g_xy_xx_xy_yy[i] * a_exp;

        g_y_0_0_0_x_xx_xy_yz[i] = 2.0 * g_xy_xx_xy_yz[i] * a_exp;

        g_y_0_0_0_x_xx_xy_zz[i] = 2.0 * g_xy_xx_xy_zz[i] * a_exp;
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_xy_xx_xz_xx, g_xy_xx_xz_xy, g_xy_xx_xz_xz, g_xy_xx_xz_yy, g_xy_xx_xz_yz, g_xy_xx_xz_zz, g_y_0_0_0_x_xx_xz_xx, g_y_0_0_0_x_xx_xz_xy, g_y_0_0_0_x_xx_xz_xz, g_y_0_0_0_x_xx_xz_yy, g_y_0_0_0_x_xx_xz_yz, g_y_0_0_0_x_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xx_xz_xx[i] = 2.0 * g_xy_xx_xz_xx[i] * a_exp;

        g_y_0_0_0_x_xx_xz_xy[i] = 2.0 * g_xy_xx_xz_xy[i] * a_exp;

        g_y_0_0_0_x_xx_xz_xz[i] = 2.0 * g_xy_xx_xz_xz[i] * a_exp;

        g_y_0_0_0_x_xx_xz_yy[i] = 2.0 * g_xy_xx_xz_yy[i] * a_exp;

        g_y_0_0_0_x_xx_xz_yz[i] = 2.0 * g_xy_xx_xz_yz[i] * a_exp;

        g_y_0_0_0_x_xx_xz_zz[i] = 2.0 * g_xy_xx_xz_zz[i] * a_exp;
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_xy_xx_yy_xx, g_xy_xx_yy_xy, g_xy_xx_yy_xz, g_xy_xx_yy_yy, g_xy_xx_yy_yz, g_xy_xx_yy_zz, g_y_0_0_0_x_xx_yy_xx, g_y_0_0_0_x_xx_yy_xy, g_y_0_0_0_x_xx_yy_xz, g_y_0_0_0_x_xx_yy_yy, g_y_0_0_0_x_xx_yy_yz, g_y_0_0_0_x_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xx_yy_xx[i] = 2.0 * g_xy_xx_yy_xx[i] * a_exp;

        g_y_0_0_0_x_xx_yy_xy[i] = 2.0 * g_xy_xx_yy_xy[i] * a_exp;

        g_y_0_0_0_x_xx_yy_xz[i] = 2.0 * g_xy_xx_yy_xz[i] * a_exp;

        g_y_0_0_0_x_xx_yy_yy[i] = 2.0 * g_xy_xx_yy_yy[i] * a_exp;

        g_y_0_0_0_x_xx_yy_yz[i] = 2.0 * g_xy_xx_yy_yz[i] * a_exp;

        g_y_0_0_0_x_xx_yy_zz[i] = 2.0 * g_xy_xx_yy_zz[i] * a_exp;
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_xy_xx_yz_xx, g_xy_xx_yz_xy, g_xy_xx_yz_xz, g_xy_xx_yz_yy, g_xy_xx_yz_yz, g_xy_xx_yz_zz, g_y_0_0_0_x_xx_yz_xx, g_y_0_0_0_x_xx_yz_xy, g_y_0_0_0_x_xx_yz_xz, g_y_0_0_0_x_xx_yz_yy, g_y_0_0_0_x_xx_yz_yz, g_y_0_0_0_x_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xx_yz_xx[i] = 2.0 * g_xy_xx_yz_xx[i] * a_exp;

        g_y_0_0_0_x_xx_yz_xy[i] = 2.0 * g_xy_xx_yz_xy[i] * a_exp;

        g_y_0_0_0_x_xx_yz_xz[i] = 2.0 * g_xy_xx_yz_xz[i] * a_exp;

        g_y_0_0_0_x_xx_yz_yy[i] = 2.0 * g_xy_xx_yz_yy[i] * a_exp;

        g_y_0_0_0_x_xx_yz_yz[i] = 2.0 * g_xy_xx_yz_yz[i] * a_exp;

        g_y_0_0_0_x_xx_yz_zz[i] = 2.0 * g_xy_xx_yz_zz[i] * a_exp;
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_xy_xx_zz_xx, g_xy_xx_zz_xy, g_xy_xx_zz_xz, g_xy_xx_zz_yy, g_xy_xx_zz_yz, g_xy_xx_zz_zz, g_y_0_0_0_x_xx_zz_xx, g_y_0_0_0_x_xx_zz_xy, g_y_0_0_0_x_xx_zz_xz, g_y_0_0_0_x_xx_zz_yy, g_y_0_0_0_x_xx_zz_yz, g_y_0_0_0_x_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xx_zz_xx[i] = 2.0 * g_xy_xx_zz_xx[i] * a_exp;

        g_y_0_0_0_x_xx_zz_xy[i] = 2.0 * g_xy_xx_zz_xy[i] * a_exp;

        g_y_0_0_0_x_xx_zz_xz[i] = 2.0 * g_xy_xx_zz_xz[i] * a_exp;

        g_y_0_0_0_x_xx_zz_yy[i] = 2.0 * g_xy_xx_zz_yy[i] * a_exp;

        g_y_0_0_0_x_xx_zz_yz[i] = 2.0 * g_xy_xx_zz_yz[i] * a_exp;

        g_y_0_0_0_x_xx_zz_zz[i] = 2.0 * g_xy_xx_zz_zz[i] * a_exp;
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_xy_xy_xx_xx, g_xy_xy_xx_xy, g_xy_xy_xx_xz, g_xy_xy_xx_yy, g_xy_xy_xx_yz, g_xy_xy_xx_zz, g_y_0_0_0_x_xy_xx_xx, g_y_0_0_0_x_xy_xx_xy, g_y_0_0_0_x_xy_xx_xz, g_y_0_0_0_x_xy_xx_yy, g_y_0_0_0_x_xy_xx_yz, g_y_0_0_0_x_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xy_xx_xx[i] = 2.0 * g_xy_xy_xx_xx[i] * a_exp;

        g_y_0_0_0_x_xy_xx_xy[i] = 2.0 * g_xy_xy_xx_xy[i] * a_exp;

        g_y_0_0_0_x_xy_xx_xz[i] = 2.0 * g_xy_xy_xx_xz[i] * a_exp;

        g_y_0_0_0_x_xy_xx_yy[i] = 2.0 * g_xy_xy_xx_yy[i] * a_exp;

        g_y_0_0_0_x_xy_xx_yz[i] = 2.0 * g_xy_xy_xx_yz[i] * a_exp;

        g_y_0_0_0_x_xy_xx_zz[i] = 2.0 * g_xy_xy_xx_zz[i] * a_exp;
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_xy_xy_xy_xx, g_xy_xy_xy_xy, g_xy_xy_xy_xz, g_xy_xy_xy_yy, g_xy_xy_xy_yz, g_xy_xy_xy_zz, g_y_0_0_0_x_xy_xy_xx, g_y_0_0_0_x_xy_xy_xy, g_y_0_0_0_x_xy_xy_xz, g_y_0_0_0_x_xy_xy_yy, g_y_0_0_0_x_xy_xy_yz, g_y_0_0_0_x_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xy_xy_xx[i] = 2.0 * g_xy_xy_xy_xx[i] * a_exp;

        g_y_0_0_0_x_xy_xy_xy[i] = 2.0 * g_xy_xy_xy_xy[i] * a_exp;

        g_y_0_0_0_x_xy_xy_xz[i] = 2.0 * g_xy_xy_xy_xz[i] * a_exp;

        g_y_0_0_0_x_xy_xy_yy[i] = 2.0 * g_xy_xy_xy_yy[i] * a_exp;

        g_y_0_0_0_x_xy_xy_yz[i] = 2.0 * g_xy_xy_xy_yz[i] * a_exp;

        g_y_0_0_0_x_xy_xy_zz[i] = 2.0 * g_xy_xy_xy_zz[i] * a_exp;
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_xy_xy_xz_xx, g_xy_xy_xz_xy, g_xy_xy_xz_xz, g_xy_xy_xz_yy, g_xy_xy_xz_yz, g_xy_xy_xz_zz, g_y_0_0_0_x_xy_xz_xx, g_y_0_0_0_x_xy_xz_xy, g_y_0_0_0_x_xy_xz_xz, g_y_0_0_0_x_xy_xz_yy, g_y_0_0_0_x_xy_xz_yz, g_y_0_0_0_x_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xy_xz_xx[i] = 2.0 * g_xy_xy_xz_xx[i] * a_exp;

        g_y_0_0_0_x_xy_xz_xy[i] = 2.0 * g_xy_xy_xz_xy[i] * a_exp;

        g_y_0_0_0_x_xy_xz_xz[i] = 2.0 * g_xy_xy_xz_xz[i] * a_exp;

        g_y_0_0_0_x_xy_xz_yy[i] = 2.0 * g_xy_xy_xz_yy[i] * a_exp;

        g_y_0_0_0_x_xy_xz_yz[i] = 2.0 * g_xy_xy_xz_yz[i] * a_exp;

        g_y_0_0_0_x_xy_xz_zz[i] = 2.0 * g_xy_xy_xz_zz[i] * a_exp;
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_xy_xy_yy_xx, g_xy_xy_yy_xy, g_xy_xy_yy_xz, g_xy_xy_yy_yy, g_xy_xy_yy_yz, g_xy_xy_yy_zz, g_y_0_0_0_x_xy_yy_xx, g_y_0_0_0_x_xy_yy_xy, g_y_0_0_0_x_xy_yy_xz, g_y_0_0_0_x_xy_yy_yy, g_y_0_0_0_x_xy_yy_yz, g_y_0_0_0_x_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xy_yy_xx[i] = 2.0 * g_xy_xy_yy_xx[i] * a_exp;

        g_y_0_0_0_x_xy_yy_xy[i] = 2.0 * g_xy_xy_yy_xy[i] * a_exp;

        g_y_0_0_0_x_xy_yy_xz[i] = 2.0 * g_xy_xy_yy_xz[i] * a_exp;

        g_y_0_0_0_x_xy_yy_yy[i] = 2.0 * g_xy_xy_yy_yy[i] * a_exp;

        g_y_0_0_0_x_xy_yy_yz[i] = 2.0 * g_xy_xy_yy_yz[i] * a_exp;

        g_y_0_0_0_x_xy_yy_zz[i] = 2.0 * g_xy_xy_yy_zz[i] * a_exp;
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_xy_xy_yz_xx, g_xy_xy_yz_xy, g_xy_xy_yz_xz, g_xy_xy_yz_yy, g_xy_xy_yz_yz, g_xy_xy_yz_zz, g_y_0_0_0_x_xy_yz_xx, g_y_0_0_0_x_xy_yz_xy, g_y_0_0_0_x_xy_yz_xz, g_y_0_0_0_x_xy_yz_yy, g_y_0_0_0_x_xy_yz_yz, g_y_0_0_0_x_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xy_yz_xx[i] = 2.0 * g_xy_xy_yz_xx[i] * a_exp;

        g_y_0_0_0_x_xy_yz_xy[i] = 2.0 * g_xy_xy_yz_xy[i] * a_exp;

        g_y_0_0_0_x_xy_yz_xz[i] = 2.0 * g_xy_xy_yz_xz[i] * a_exp;

        g_y_0_0_0_x_xy_yz_yy[i] = 2.0 * g_xy_xy_yz_yy[i] * a_exp;

        g_y_0_0_0_x_xy_yz_yz[i] = 2.0 * g_xy_xy_yz_yz[i] * a_exp;

        g_y_0_0_0_x_xy_yz_zz[i] = 2.0 * g_xy_xy_yz_zz[i] * a_exp;
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_xy_xy_zz_xx, g_xy_xy_zz_xy, g_xy_xy_zz_xz, g_xy_xy_zz_yy, g_xy_xy_zz_yz, g_xy_xy_zz_zz, g_y_0_0_0_x_xy_zz_xx, g_y_0_0_0_x_xy_zz_xy, g_y_0_0_0_x_xy_zz_xz, g_y_0_0_0_x_xy_zz_yy, g_y_0_0_0_x_xy_zz_yz, g_y_0_0_0_x_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xy_zz_xx[i] = 2.0 * g_xy_xy_zz_xx[i] * a_exp;

        g_y_0_0_0_x_xy_zz_xy[i] = 2.0 * g_xy_xy_zz_xy[i] * a_exp;

        g_y_0_0_0_x_xy_zz_xz[i] = 2.0 * g_xy_xy_zz_xz[i] * a_exp;

        g_y_0_0_0_x_xy_zz_yy[i] = 2.0 * g_xy_xy_zz_yy[i] * a_exp;

        g_y_0_0_0_x_xy_zz_yz[i] = 2.0 * g_xy_xy_zz_yz[i] * a_exp;

        g_y_0_0_0_x_xy_zz_zz[i] = 2.0 * g_xy_xy_zz_zz[i] * a_exp;
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_xy_xz_xx_xx, g_xy_xz_xx_xy, g_xy_xz_xx_xz, g_xy_xz_xx_yy, g_xy_xz_xx_yz, g_xy_xz_xx_zz, g_y_0_0_0_x_xz_xx_xx, g_y_0_0_0_x_xz_xx_xy, g_y_0_0_0_x_xz_xx_xz, g_y_0_0_0_x_xz_xx_yy, g_y_0_0_0_x_xz_xx_yz, g_y_0_0_0_x_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xz_xx_xx[i] = 2.0 * g_xy_xz_xx_xx[i] * a_exp;

        g_y_0_0_0_x_xz_xx_xy[i] = 2.0 * g_xy_xz_xx_xy[i] * a_exp;

        g_y_0_0_0_x_xz_xx_xz[i] = 2.0 * g_xy_xz_xx_xz[i] * a_exp;

        g_y_0_0_0_x_xz_xx_yy[i] = 2.0 * g_xy_xz_xx_yy[i] * a_exp;

        g_y_0_0_0_x_xz_xx_yz[i] = 2.0 * g_xy_xz_xx_yz[i] * a_exp;

        g_y_0_0_0_x_xz_xx_zz[i] = 2.0 * g_xy_xz_xx_zz[i] * a_exp;
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_xy_xz_xy_xx, g_xy_xz_xy_xy, g_xy_xz_xy_xz, g_xy_xz_xy_yy, g_xy_xz_xy_yz, g_xy_xz_xy_zz, g_y_0_0_0_x_xz_xy_xx, g_y_0_0_0_x_xz_xy_xy, g_y_0_0_0_x_xz_xy_xz, g_y_0_0_0_x_xz_xy_yy, g_y_0_0_0_x_xz_xy_yz, g_y_0_0_0_x_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xz_xy_xx[i] = 2.0 * g_xy_xz_xy_xx[i] * a_exp;

        g_y_0_0_0_x_xz_xy_xy[i] = 2.0 * g_xy_xz_xy_xy[i] * a_exp;

        g_y_0_0_0_x_xz_xy_xz[i] = 2.0 * g_xy_xz_xy_xz[i] * a_exp;

        g_y_0_0_0_x_xz_xy_yy[i] = 2.0 * g_xy_xz_xy_yy[i] * a_exp;

        g_y_0_0_0_x_xz_xy_yz[i] = 2.0 * g_xy_xz_xy_yz[i] * a_exp;

        g_y_0_0_0_x_xz_xy_zz[i] = 2.0 * g_xy_xz_xy_zz[i] * a_exp;
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_xy_xz_xz_xx, g_xy_xz_xz_xy, g_xy_xz_xz_xz, g_xy_xz_xz_yy, g_xy_xz_xz_yz, g_xy_xz_xz_zz, g_y_0_0_0_x_xz_xz_xx, g_y_0_0_0_x_xz_xz_xy, g_y_0_0_0_x_xz_xz_xz, g_y_0_0_0_x_xz_xz_yy, g_y_0_0_0_x_xz_xz_yz, g_y_0_0_0_x_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xz_xz_xx[i] = 2.0 * g_xy_xz_xz_xx[i] * a_exp;

        g_y_0_0_0_x_xz_xz_xy[i] = 2.0 * g_xy_xz_xz_xy[i] * a_exp;

        g_y_0_0_0_x_xz_xz_xz[i] = 2.0 * g_xy_xz_xz_xz[i] * a_exp;

        g_y_0_0_0_x_xz_xz_yy[i] = 2.0 * g_xy_xz_xz_yy[i] * a_exp;

        g_y_0_0_0_x_xz_xz_yz[i] = 2.0 * g_xy_xz_xz_yz[i] * a_exp;

        g_y_0_0_0_x_xz_xz_zz[i] = 2.0 * g_xy_xz_xz_zz[i] * a_exp;
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_xy_xz_yy_xx, g_xy_xz_yy_xy, g_xy_xz_yy_xz, g_xy_xz_yy_yy, g_xy_xz_yy_yz, g_xy_xz_yy_zz, g_y_0_0_0_x_xz_yy_xx, g_y_0_0_0_x_xz_yy_xy, g_y_0_0_0_x_xz_yy_xz, g_y_0_0_0_x_xz_yy_yy, g_y_0_0_0_x_xz_yy_yz, g_y_0_0_0_x_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xz_yy_xx[i] = 2.0 * g_xy_xz_yy_xx[i] * a_exp;

        g_y_0_0_0_x_xz_yy_xy[i] = 2.0 * g_xy_xz_yy_xy[i] * a_exp;

        g_y_0_0_0_x_xz_yy_xz[i] = 2.0 * g_xy_xz_yy_xz[i] * a_exp;

        g_y_0_0_0_x_xz_yy_yy[i] = 2.0 * g_xy_xz_yy_yy[i] * a_exp;

        g_y_0_0_0_x_xz_yy_yz[i] = 2.0 * g_xy_xz_yy_yz[i] * a_exp;

        g_y_0_0_0_x_xz_yy_zz[i] = 2.0 * g_xy_xz_yy_zz[i] * a_exp;
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_xy_xz_yz_xx, g_xy_xz_yz_xy, g_xy_xz_yz_xz, g_xy_xz_yz_yy, g_xy_xz_yz_yz, g_xy_xz_yz_zz, g_y_0_0_0_x_xz_yz_xx, g_y_0_0_0_x_xz_yz_xy, g_y_0_0_0_x_xz_yz_xz, g_y_0_0_0_x_xz_yz_yy, g_y_0_0_0_x_xz_yz_yz, g_y_0_0_0_x_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xz_yz_xx[i] = 2.0 * g_xy_xz_yz_xx[i] * a_exp;

        g_y_0_0_0_x_xz_yz_xy[i] = 2.0 * g_xy_xz_yz_xy[i] * a_exp;

        g_y_0_0_0_x_xz_yz_xz[i] = 2.0 * g_xy_xz_yz_xz[i] * a_exp;

        g_y_0_0_0_x_xz_yz_yy[i] = 2.0 * g_xy_xz_yz_yy[i] * a_exp;

        g_y_0_0_0_x_xz_yz_yz[i] = 2.0 * g_xy_xz_yz_yz[i] * a_exp;

        g_y_0_0_0_x_xz_yz_zz[i] = 2.0 * g_xy_xz_yz_zz[i] * a_exp;
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_xy_xz_zz_xx, g_xy_xz_zz_xy, g_xy_xz_zz_xz, g_xy_xz_zz_yy, g_xy_xz_zz_yz, g_xy_xz_zz_zz, g_y_0_0_0_x_xz_zz_xx, g_y_0_0_0_x_xz_zz_xy, g_y_0_0_0_x_xz_zz_xz, g_y_0_0_0_x_xz_zz_yy, g_y_0_0_0_x_xz_zz_yz, g_y_0_0_0_x_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xz_zz_xx[i] = 2.0 * g_xy_xz_zz_xx[i] * a_exp;

        g_y_0_0_0_x_xz_zz_xy[i] = 2.0 * g_xy_xz_zz_xy[i] * a_exp;

        g_y_0_0_0_x_xz_zz_xz[i] = 2.0 * g_xy_xz_zz_xz[i] * a_exp;

        g_y_0_0_0_x_xz_zz_yy[i] = 2.0 * g_xy_xz_zz_yy[i] * a_exp;

        g_y_0_0_0_x_xz_zz_yz[i] = 2.0 * g_xy_xz_zz_yz[i] * a_exp;

        g_y_0_0_0_x_xz_zz_zz[i] = 2.0 * g_xy_xz_zz_zz[i] * a_exp;
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_xy_yy_xx_xx, g_xy_yy_xx_xy, g_xy_yy_xx_xz, g_xy_yy_xx_yy, g_xy_yy_xx_yz, g_xy_yy_xx_zz, g_y_0_0_0_x_yy_xx_xx, g_y_0_0_0_x_yy_xx_xy, g_y_0_0_0_x_yy_xx_xz, g_y_0_0_0_x_yy_xx_yy, g_y_0_0_0_x_yy_xx_yz, g_y_0_0_0_x_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yy_xx_xx[i] = 2.0 * g_xy_yy_xx_xx[i] * a_exp;

        g_y_0_0_0_x_yy_xx_xy[i] = 2.0 * g_xy_yy_xx_xy[i] * a_exp;

        g_y_0_0_0_x_yy_xx_xz[i] = 2.0 * g_xy_yy_xx_xz[i] * a_exp;

        g_y_0_0_0_x_yy_xx_yy[i] = 2.0 * g_xy_yy_xx_yy[i] * a_exp;

        g_y_0_0_0_x_yy_xx_yz[i] = 2.0 * g_xy_yy_xx_yz[i] * a_exp;

        g_y_0_0_0_x_yy_xx_zz[i] = 2.0 * g_xy_yy_xx_zz[i] * a_exp;
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_xy_yy_xy_xx, g_xy_yy_xy_xy, g_xy_yy_xy_xz, g_xy_yy_xy_yy, g_xy_yy_xy_yz, g_xy_yy_xy_zz, g_y_0_0_0_x_yy_xy_xx, g_y_0_0_0_x_yy_xy_xy, g_y_0_0_0_x_yy_xy_xz, g_y_0_0_0_x_yy_xy_yy, g_y_0_0_0_x_yy_xy_yz, g_y_0_0_0_x_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yy_xy_xx[i] = 2.0 * g_xy_yy_xy_xx[i] * a_exp;

        g_y_0_0_0_x_yy_xy_xy[i] = 2.0 * g_xy_yy_xy_xy[i] * a_exp;

        g_y_0_0_0_x_yy_xy_xz[i] = 2.0 * g_xy_yy_xy_xz[i] * a_exp;

        g_y_0_0_0_x_yy_xy_yy[i] = 2.0 * g_xy_yy_xy_yy[i] * a_exp;

        g_y_0_0_0_x_yy_xy_yz[i] = 2.0 * g_xy_yy_xy_yz[i] * a_exp;

        g_y_0_0_0_x_yy_xy_zz[i] = 2.0 * g_xy_yy_xy_zz[i] * a_exp;
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_xy_yy_xz_xx, g_xy_yy_xz_xy, g_xy_yy_xz_xz, g_xy_yy_xz_yy, g_xy_yy_xz_yz, g_xy_yy_xz_zz, g_y_0_0_0_x_yy_xz_xx, g_y_0_0_0_x_yy_xz_xy, g_y_0_0_0_x_yy_xz_xz, g_y_0_0_0_x_yy_xz_yy, g_y_0_0_0_x_yy_xz_yz, g_y_0_0_0_x_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yy_xz_xx[i] = 2.0 * g_xy_yy_xz_xx[i] * a_exp;

        g_y_0_0_0_x_yy_xz_xy[i] = 2.0 * g_xy_yy_xz_xy[i] * a_exp;

        g_y_0_0_0_x_yy_xz_xz[i] = 2.0 * g_xy_yy_xz_xz[i] * a_exp;

        g_y_0_0_0_x_yy_xz_yy[i] = 2.0 * g_xy_yy_xz_yy[i] * a_exp;

        g_y_0_0_0_x_yy_xz_yz[i] = 2.0 * g_xy_yy_xz_yz[i] * a_exp;

        g_y_0_0_0_x_yy_xz_zz[i] = 2.0 * g_xy_yy_xz_zz[i] * a_exp;
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_xy_yy_yy_xx, g_xy_yy_yy_xy, g_xy_yy_yy_xz, g_xy_yy_yy_yy, g_xy_yy_yy_yz, g_xy_yy_yy_zz, g_y_0_0_0_x_yy_yy_xx, g_y_0_0_0_x_yy_yy_xy, g_y_0_0_0_x_yy_yy_xz, g_y_0_0_0_x_yy_yy_yy, g_y_0_0_0_x_yy_yy_yz, g_y_0_0_0_x_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yy_yy_xx[i] = 2.0 * g_xy_yy_yy_xx[i] * a_exp;

        g_y_0_0_0_x_yy_yy_xy[i] = 2.0 * g_xy_yy_yy_xy[i] * a_exp;

        g_y_0_0_0_x_yy_yy_xz[i] = 2.0 * g_xy_yy_yy_xz[i] * a_exp;

        g_y_0_0_0_x_yy_yy_yy[i] = 2.0 * g_xy_yy_yy_yy[i] * a_exp;

        g_y_0_0_0_x_yy_yy_yz[i] = 2.0 * g_xy_yy_yy_yz[i] * a_exp;

        g_y_0_0_0_x_yy_yy_zz[i] = 2.0 * g_xy_yy_yy_zz[i] * a_exp;
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_xy_yy_yz_xx, g_xy_yy_yz_xy, g_xy_yy_yz_xz, g_xy_yy_yz_yy, g_xy_yy_yz_yz, g_xy_yy_yz_zz, g_y_0_0_0_x_yy_yz_xx, g_y_0_0_0_x_yy_yz_xy, g_y_0_0_0_x_yy_yz_xz, g_y_0_0_0_x_yy_yz_yy, g_y_0_0_0_x_yy_yz_yz, g_y_0_0_0_x_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yy_yz_xx[i] = 2.0 * g_xy_yy_yz_xx[i] * a_exp;

        g_y_0_0_0_x_yy_yz_xy[i] = 2.0 * g_xy_yy_yz_xy[i] * a_exp;

        g_y_0_0_0_x_yy_yz_xz[i] = 2.0 * g_xy_yy_yz_xz[i] * a_exp;

        g_y_0_0_0_x_yy_yz_yy[i] = 2.0 * g_xy_yy_yz_yy[i] * a_exp;

        g_y_0_0_0_x_yy_yz_yz[i] = 2.0 * g_xy_yy_yz_yz[i] * a_exp;

        g_y_0_0_0_x_yy_yz_zz[i] = 2.0 * g_xy_yy_yz_zz[i] * a_exp;
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_xy_yy_zz_xx, g_xy_yy_zz_xy, g_xy_yy_zz_xz, g_xy_yy_zz_yy, g_xy_yy_zz_yz, g_xy_yy_zz_zz, g_y_0_0_0_x_yy_zz_xx, g_y_0_0_0_x_yy_zz_xy, g_y_0_0_0_x_yy_zz_xz, g_y_0_0_0_x_yy_zz_yy, g_y_0_0_0_x_yy_zz_yz, g_y_0_0_0_x_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yy_zz_xx[i] = 2.0 * g_xy_yy_zz_xx[i] * a_exp;

        g_y_0_0_0_x_yy_zz_xy[i] = 2.0 * g_xy_yy_zz_xy[i] * a_exp;

        g_y_0_0_0_x_yy_zz_xz[i] = 2.0 * g_xy_yy_zz_xz[i] * a_exp;

        g_y_0_0_0_x_yy_zz_yy[i] = 2.0 * g_xy_yy_zz_yy[i] * a_exp;

        g_y_0_0_0_x_yy_zz_yz[i] = 2.0 * g_xy_yy_zz_yz[i] * a_exp;

        g_y_0_0_0_x_yy_zz_zz[i] = 2.0 * g_xy_yy_zz_zz[i] * a_exp;
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_xy_yz_xx_xx, g_xy_yz_xx_xy, g_xy_yz_xx_xz, g_xy_yz_xx_yy, g_xy_yz_xx_yz, g_xy_yz_xx_zz, g_y_0_0_0_x_yz_xx_xx, g_y_0_0_0_x_yz_xx_xy, g_y_0_0_0_x_yz_xx_xz, g_y_0_0_0_x_yz_xx_yy, g_y_0_0_0_x_yz_xx_yz, g_y_0_0_0_x_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yz_xx_xx[i] = 2.0 * g_xy_yz_xx_xx[i] * a_exp;

        g_y_0_0_0_x_yz_xx_xy[i] = 2.0 * g_xy_yz_xx_xy[i] * a_exp;

        g_y_0_0_0_x_yz_xx_xz[i] = 2.0 * g_xy_yz_xx_xz[i] * a_exp;

        g_y_0_0_0_x_yz_xx_yy[i] = 2.0 * g_xy_yz_xx_yy[i] * a_exp;

        g_y_0_0_0_x_yz_xx_yz[i] = 2.0 * g_xy_yz_xx_yz[i] * a_exp;

        g_y_0_0_0_x_yz_xx_zz[i] = 2.0 * g_xy_yz_xx_zz[i] * a_exp;
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_xy_yz_xy_xx, g_xy_yz_xy_xy, g_xy_yz_xy_xz, g_xy_yz_xy_yy, g_xy_yz_xy_yz, g_xy_yz_xy_zz, g_y_0_0_0_x_yz_xy_xx, g_y_0_0_0_x_yz_xy_xy, g_y_0_0_0_x_yz_xy_xz, g_y_0_0_0_x_yz_xy_yy, g_y_0_0_0_x_yz_xy_yz, g_y_0_0_0_x_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yz_xy_xx[i] = 2.0 * g_xy_yz_xy_xx[i] * a_exp;

        g_y_0_0_0_x_yz_xy_xy[i] = 2.0 * g_xy_yz_xy_xy[i] * a_exp;

        g_y_0_0_0_x_yz_xy_xz[i] = 2.0 * g_xy_yz_xy_xz[i] * a_exp;

        g_y_0_0_0_x_yz_xy_yy[i] = 2.0 * g_xy_yz_xy_yy[i] * a_exp;

        g_y_0_0_0_x_yz_xy_yz[i] = 2.0 * g_xy_yz_xy_yz[i] * a_exp;

        g_y_0_0_0_x_yz_xy_zz[i] = 2.0 * g_xy_yz_xy_zz[i] * a_exp;
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_xy_yz_xz_xx, g_xy_yz_xz_xy, g_xy_yz_xz_xz, g_xy_yz_xz_yy, g_xy_yz_xz_yz, g_xy_yz_xz_zz, g_y_0_0_0_x_yz_xz_xx, g_y_0_0_0_x_yz_xz_xy, g_y_0_0_0_x_yz_xz_xz, g_y_0_0_0_x_yz_xz_yy, g_y_0_0_0_x_yz_xz_yz, g_y_0_0_0_x_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yz_xz_xx[i] = 2.0 * g_xy_yz_xz_xx[i] * a_exp;

        g_y_0_0_0_x_yz_xz_xy[i] = 2.0 * g_xy_yz_xz_xy[i] * a_exp;

        g_y_0_0_0_x_yz_xz_xz[i] = 2.0 * g_xy_yz_xz_xz[i] * a_exp;

        g_y_0_0_0_x_yz_xz_yy[i] = 2.0 * g_xy_yz_xz_yy[i] * a_exp;

        g_y_0_0_0_x_yz_xz_yz[i] = 2.0 * g_xy_yz_xz_yz[i] * a_exp;

        g_y_0_0_0_x_yz_xz_zz[i] = 2.0 * g_xy_yz_xz_zz[i] * a_exp;
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_xy_yz_yy_xx, g_xy_yz_yy_xy, g_xy_yz_yy_xz, g_xy_yz_yy_yy, g_xy_yz_yy_yz, g_xy_yz_yy_zz, g_y_0_0_0_x_yz_yy_xx, g_y_0_0_0_x_yz_yy_xy, g_y_0_0_0_x_yz_yy_xz, g_y_0_0_0_x_yz_yy_yy, g_y_0_0_0_x_yz_yy_yz, g_y_0_0_0_x_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yz_yy_xx[i] = 2.0 * g_xy_yz_yy_xx[i] * a_exp;

        g_y_0_0_0_x_yz_yy_xy[i] = 2.0 * g_xy_yz_yy_xy[i] * a_exp;

        g_y_0_0_0_x_yz_yy_xz[i] = 2.0 * g_xy_yz_yy_xz[i] * a_exp;

        g_y_0_0_0_x_yz_yy_yy[i] = 2.0 * g_xy_yz_yy_yy[i] * a_exp;

        g_y_0_0_0_x_yz_yy_yz[i] = 2.0 * g_xy_yz_yy_yz[i] * a_exp;

        g_y_0_0_0_x_yz_yy_zz[i] = 2.0 * g_xy_yz_yy_zz[i] * a_exp;
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_xy_yz_yz_xx, g_xy_yz_yz_xy, g_xy_yz_yz_xz, g_xy_yz_yz_yy, g_xy_yz_yz_yz, g_xy_yz_yz_zz, g_y_0_0_0_x_yz_yz_xx, g_y_0_0_0_x_yz_yz_xy, g_y_0_0_0_x_yz_yz_xz, g_y_0_0_0_x_yz_yz_yy, g_y_0_0_0_x_yz_yz_yz, g_y_0_0_0_x_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yz_yz_xx[i] = 2.0 * g_xy_yz_yz_xx[i] * a_exp;

        g_y_0_0_0_x_yz_yz_xy[i] = 2.0 * g_xy_yz_yz_xy[i] * a_exp;

        g_y_0_0_0_x_yz_yz_xz[i] = 2.0 * g_xy_yz_yz_xz[i] * a_exp;

        g_y_0_0_0_x_yz_yz_yy[i] = 2.0 * g_xy_yz_yz_yy[i] * a_exp;

        g_y_0_0_0_x_yz_yz_yz[i] = 2.0 * g_xy_yz_yz_yz[i] * a_exp;

        g_y_0_0_0_x_yz_yz_zz[i] = 2.0 * g_xy_yz_yz_zz[i] * a_exp;
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_xy_yz_zz_xx, g_xy_yz_zz_xy, g_xy_yz_zz_xz, g_xy_yz_zz_yy, g_xy_yz_zz_yz, g_xy_yz_zz_zz, g_y_0_0_0_x_yz_zz_xx, g_y_0_0_0_x_yz_zz_xy, g_y_0_0_0_x_yz_zz_xz, g_y_0_0_0_x_yz_zz_yy, g_y_0_0_0_x_yz_zz_yz, g_y_0_0_0_x_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yz_zz_xx[i] = 2.0 * g_xy_yz_zz_xx[i] * a_exp;

        g_y_0_0_0_x_yz_zz_xy[i] = 2.0 * g_xy_yz_zz_xy[i] * a_exp;

        g_y_0_0_0_x_yz_zz_xz[i] = 2.0 * g_xy_yz_zz_xz[i] * a_exp;

        g_y_0_0_0_x_yz_zz_yy[i] = 2.0 * g_xy_yz_zz_yy[i] * a_exp;

        g_y_0_0_0_x_yz_zz_yz[i] = 2.0 * g_xy_yz_zz_yz[i] * a_exp;

        g_y_0_0_0_x_yz_zz_zz[i] = 2.0 * g_xy_yz_zz_zz[i] * a_exp;
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_xy_zz_xx_xx, g_xy_zz_xx_xy, g_xy_zz_xx_xz, g_xy_zz_xx_yy, g_xy_zz_xx_yz, g_xy_zz_xx_zz, g_y_0_0_0_x_zz_xx_xx, g_y_0_0_0_x_zz_xx_xy, g_y_0_0_0_x_zz_xx_xz, g_y_0_0_0_x_zz_xx_yy, g_y_0_0_0_x_zz_xx_yz, g_y_0_0_0_x_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_zz_xx_xx[i] = 2.0 * g_xy_zz_xx_xx[i] * a_exp;

        g_y_0_0_0_x_zz_xx_xy[i] = 2.0 * g_xy_zz_xx_xy[i] * a_exp;

        g_y_0_0_0_x_zz_xx_xz[i] = 2.0 * g_xy_zz_xx_xz[i] * a_exp;

        g_y_0_0_0_x_zz_xx_yy[i] = 2.0 * g_xy_zz_xx_yy[i] * a_exp;

        g_y_0_0_0_x_zz_xx_yz[i] = 2.0 * g_xy_zz_xx_yz[i] * a_exp;

        g_y_0_0_0_x_zz_xx_zz[i] = 2.0 * g_xy_zz_xx_zz[i] * a_exp;
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_xy_zz_xy_xx, g_xy_zz_xy_xy, g_xy_zz_xy_xz, g_xy_zz_xy_yy, g_xy_zz_xy_yz, g_xy_zz_xy_zz, g_y_0_0_0_x_zz_xy_xx, g_y_0_0_0_x_zz_xy_xy, g_y_0_0_0_x_zz_xy_xz, g_y_0_0_0_x_zz_xy_yy, g_y_0_0_0_x_zz_xy_yz, g_y_0_0_0_x_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_zz_xy_xx[i] = 2.0 * g_xy_zz_xy_xx[i] * a_exp;

        g_y_0_0_0_x_zz_xy_xy[i] = 2.0 * g_xy_zz_xy_xy[i] * a_exp;

        g_y_0_0_0_x_zz_xy_xz[i] = 2.0 * g_xy_zz_xy_xz[i] * a_exp;

        g_y_0_0_0_x_zz_xy_yy[i] = 2.0 * g_xy_zz_xy_yy[i] * a_exp;

        g_y_0_0_0_x_zz_xy_yz[i] = 2.0 * g_xy_zz_xy_yz[i] * a_exp;

        g_y_0_0_0_x_zz_xy_zz[i] = 2.0 * g_xy_zz_xy_zz[i] * a_exp;
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_xy_zz_xz_xx, g_xy_zz_xz_xy, g_xy_zz_xz_xz, g_xy_zz_xz_yy, g_xy_zz_xz_yz, g_xy_zz_xz_zz, g_y_0_0_0_x_zz_xz_xx, g_y_0_0_0_x_zz_xz_xy, g_y_0_0_0_x_zz_xz_xz, g_y_0_0_0_x_zz_xz_yy, g_y_0_0_0_x_zz_xz_yz, g_y_0_0_0_x_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_zz_xz_xx[i] = 2.0 * g_xy_zz_xz_xx[i] * a_exp;

        g_y_0_0_0_x_zz_xz_xy[i] = 2.0 * g_xy_zz_xz_xy[i] * a_exp;

        g_y_0_0_0_x_zz_xz_xz[i] = 2.0 * g_xy_zz_xz_xz[i] * a_exp;

        g_y_0_0_0_x_zz_xz_yy[i] = 2.0 * g_xy_zz_xz_yy[i] * a_exp;

        g_y_0_0_0_x_zz_xz_yz[i] = 2.0 * g_xy_zz_xz_yz[i] * a_exp;

        g_y_0_0_0_x_zz_xz_zz[i] = 2.0 * g_xy_zz_xz_zz[i] * a_exp;
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_xy_zz_yy_xx, g_xy_zz_yy_xy, g_xy_zz_yy_xz, g_xy_zz_yy_yy, g_xy_zz_yy_yz, g_xy_zz_yy_zz, g_y_0_0_0_x_zz_yy_xx, g_y_0_0_0_x_zz_yy_xy, g_y_0_0_0_x_zz_yy_xz, g_y_0_0_0_x_zz_yy_yy, g_y_0_0_0_x_zz_yy_yz, g_y_0_0_0_x_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_zz_yy_xx[i] = 2.0 * g_xy_zz_yy_xx[i] * a_exp;

        g_y_0_0_0_x_zz_yy_xy[i] = 2.0 * g_xy_zz_yy_xy[i] * a_exp;

        g_y_0_0_0_x_zz_yy_xz[i] = 2.0 * g_xy_zz_yy_xz[i] * a_exp;

        g_y_0_0_0_x_zz_yy_yy[i] = 2.0 * g_xy_zz_yy_yy[i] * a_exp;

        g_y_0_0_0_x_zz_yy_yz[i] = 2.0 * g_xy_zz_yy_yz[i] * a_exp;

        g_y_0_0_0_x_zz_yy_zz[i] = 2.0 * g_xy_zz_yy_zz[i] * a_exp;
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_xy_zz_yz_xx, g_xy_zz_yz_xy, g_xy_zz_yz_xz, g_xy_zz_yz_yy, g_xy_zz_yz_yz, g_xy_zz_yz_zz, g_y_0_0_0_x_zz_yz_xx, g_y_0_0_0_x_zz_yz_xy, g_y_0_0_0_x_zz_yz_xz, g_y_0_0_0_x_zz_yz_yy, g_y_0_0_0_x_zz_yz_yz, g_y_0_0_0_x_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_zz_yz_xx[i] = 2.0 * g_xy_zz_yz_xx[i] * a_exp;

        g_y_0_0_0_x_zz_yz_xy[i] = 2.0 * g_xy_zz_yz_xy[i] * a_exp;

        g_y_0_0_0_x_zz_yz_xz[i] = 2.0 * g_xy_zz_yz_xz[i] * a_exp;

        g_y_0_0_0_x_zz_yz_yy[i] = 2.0 * g_xy_zz_yz_yy[i] * a_exp;

        g_y_0_0_0_x_zz_yz_yz[i] = 2.0 * g_xy_zz_yz_yz[i] * a_exp;

        g_y_0_0_0_x_zz_yz_zz[i] = 2.0 * g_xy_zz_yz_zz[i] * a_exp;
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_xy_zz_zz_xx, g_xy_zz_zz_xy, g_xy_zz_zz_xz, g_xy_zz_zz_yy, g_xy_zz_zz_yz, g_xy_zz_zz_zz, g_y_0_0_0_x_zz_zz_xx, g_y_0_0_0_x_zz_zz_xy, g_y_0_0_0_x_zz_zz_xz, g_y_0_0_0_x_zz_zz_yy, g_y_0_0_0_x_zz_zz_yz, g_y_0_0_0_x_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_zz_zz_xx[i] = 2.0 * g_xy_zz_zz_xx[i] * a_exp;

        g_y_0_0_0_x_zz_zz_xy[i] = 2.0 * g_xy_zz_zz_xy[i] * a_exp;

        g_y_0_0_0_x_zz_zz_xz[i] = 2.0 * g_xy_zz_zz_xz[i] * a_exp;

        g_y_0_0_0_x_zz_zz_yy[i] = 2.0 * g_xy_zz_zz_yy[i] * a_exp;

        g_y_0_0_0_x_zz_zz_yz[i] = 2.0 * g_xy_zz_zz_yz[i] * a_exp;

        g_y_0_0_0_x_zz_zz_zz[i] = 2.0 * g_xy_zz_zz_zz[i] * a_exp;
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_0_xx_xx_xx, g_0_xx_xx_xy, g_0_xx_xx_xz, g_0_xx_xx_yy, g_0_xx_xx_yz, g_0_xx_xx_zz, g_y_0_0_0_y_xx_xx_xx, g_y_0_0_0_y_xx_xx_xy, g_y_0_0_0_y_xx_xx_xz, g_y_0_0_0_y_xx_xx_yy, g_y_0_0_0_y_xx_xx_yz, g_y_0_0_0_y_xx_xx_zz, g_yy_xx_xx_xx, g_yy_xx_xx_xy, g_yy_xx_xx_xz, g_yy_xx_xx_yy, g_yy_xx_xx_yz, g_yy_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xx_xx_xx[i] = -g_0_xx_xx_xx[i] + 2.0 * g_yy_xx_xx_xx[i] * a_exp;

        g_y_0_0_0_y_xx_xx_xy[i] = -g_0_xx_xx_xy[i] + 2.0 * g_yy_xx_xx_xy[i] * a_exp;

        g_y_0_0_0_y_xx_xx_xz[i] = -g_0_xx_xx_xz[i] + 2.0 * g_yy_xx_xx_xz[i] * a_exp;

        g_y_0_0_0_y_xx_xx_yy[i] = -g_0_xx_xx_yy[i] + 2.0 * g_yy_xx_xx_yy[i] * a_exp;

        g_y_0_0_0_y_xx_xx_yz[i] = -g_0_xx_xx_yz[i] + 2.0 * g_yy_xx_xx_yz[i] * a_exp;

        g_y_0_0_0_y_xx_xx_zz[i] = -g_0_xx_xx_zz[i] + 2.0 * g_yy_xx_xx_zz[i] * a_exp;
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_0_xx_xy_xx, g_0_xx_xy_xy, g_0_xx_xy_xz, g_0_xx_xy_yy, g_0_xx_xy_yz, g_0_xx_xy_zz, g_y_0_0_0_y_xx_xy_xx, g_y_0_0_0_y_xx_xy_xy, g_y_0_0_0_y_xx_xy_xz, g_y_0_0_0_y_xx_xy_yy, g_y_0_0_0_y_xx_xy_yz, g_y_0_0_0_y_xx_xy_zz, g_yy_xx_xy_xx, g_yy_xx_xy_xy, g_yy_xx_xy_xz, g_yy_xx_xy_yy, g_yy_xx_xy_yz, g_yy_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xx_xy_xx[i] = -g_0_xx_xy_xx[i] + 2.0 * g_yy_xx_xy_xx[i] * a_exp;

        g_y_0_0_0_y_xx_xy_xy[i] = -g_0_xx_xy_xy[i] + 2.0 * g_yy_xx_xy_xy[i] * a_exp;

        g_y_0_0_0_y_xx_xy_xz[i] = -g_0_xx_xy_xz[i] + 2.0 * g_yy_xx_xy_xz[i] * a_exp;

        g_y_0_0_0_y_xx_xy_yy[i] = -g_0_xx_xy_yy[i] + 2.0 * g_yy_xx_xy_yy[i] * a_exp;

        g_y_0_0_0_y_xx_xy_yz[i] = -g_0_xx_xy_yz[i] + 2.0 * g_yy_xx_xy_yz[i] * a_exp;

        g_y_0_0_0_y_xx_xy_zz[i] = -g_0_xx_xy_zz[i] + 2.0 * g_yy_xx_xy_zz[i] * a_exp;
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_0_xx_xz_xx, g_0_xx_xz_xy, g_0_xx_xz_xz, g_0_xx_xz_yy, g_0_xx_xz_yz, g_0_xx_xz_zz, g_y_0_0_0_y_xx_xz_xx, g_y_0_0_0_y_xx_xz_xy, g_y_0_0_0_y_xx_xz_xz, g_y_0_0_0_y_xx_xz_yy, g_y_0_0_0_y_xx_xz_yz, g_y_0_0_0_y_xx_xz_zz, g_yy_xx_xz_xx, g_yy_xx_xz_xy, g_yy_xx_xz_xz, g_yy_xx_xz_yy, g_yy_xx_xz_yz, g_yy_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xx_xz_xx[i] = -g_0_xx_xz_xx[i] + 2.0 * g_yy_xx_xz_xx[i] * a_exp;

        g_y_0_0_0_y_xx_xz_xy[i] = -g_0_xx_xz_xy[i] + 2.0 * g_yy_xx_xz_xy[i] * a_exp;

        g_y_0_0_0_y_xx_xz_xz[i] = -g_0_xx_xz_xz[i] + 2.0 * g_yy_xx_xz_xz[i] * a_exp;

        g_y_0_0_0_y_xx_xz_yy[i] = -g_0_xx_xz_yy[i] + 2.0 * g_yy_xx_xz_yy[i] * a_exp;

        g_y_0_0_0_y_xx_xz_yz[i] = -g_0_xx_xz_yz[i] + 2.0 * g_yy_xx_xz_yz[i] * a_exp;

        g_y_0_0_0_y_xx_xz_zz[i] = -g_0_xx_xz_zz[i] + 2.0 * g_yy_xx_xz_zz[i] * a_exp;
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_0_xx_yy_xx, g_0_xx_yy_xy, g_0_xx_yy_xz, g_0_xx_yy_yy, g_0_xx_yy_yz, g_0_xx_yy_zz, g_y_0_0_0_y_xx_yy_xx, g_y_0_0_0_y_xx_yy_xy, g_y_0_0_0_y_xx_yy_xz, g_y_0_0_0_y_xx_yy_yy, g_y_0_0_0_y_xx_yy_yz, g_y_0_0_0_y_xx_yy_zz, g_yy_xx_yy_xx, g_yy_xx_yy_xy, g_yy_xx_yy_xz, g_yy_xx_yy_yy, g_yy_xx_yy_yz, g_yy_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xx_yy_xx[i] = -g_0_xx_yy_xx[i] + 2.0 * g_yy_xx_yy_xx[i] * a_exp;

        g_y_0_0_0_y_xx_yy_xy[i] = -g_0_xx_yy_xy[i] + 2.0 * g_yy_xx_yy_xy[i] * a_exp;

        g_y_0_0_0_y_xx_yy_xz[i] = -g_0_xx_yy_xz[i] + 2.0 * g_yy_xx_yy_xz[i] * a_exp;

        g_y_0_0_0_y_xx_yy_yy[i] = -g_0_xx_yy_yy[i] + 2.0 * g_yy_xx_yy_yy[i] * a_exp;

        g_y_0_0_0_y_xx_yy_yz[i] = -g_0_xx_yy_yz[i] + 2.0 * g_yy_xx_yy_yz[i] * a_exp;

        g_y_0_0_0_y_xx_yy_zz[i] = -g_0_xx_yy_zz[i] + 2.0 * g_yy_xx_yy_zz[i] * a_exp;
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_0_xx_yz_xx, g_0_xx_yz_xy, g_0_xx_yz_xz, g_0_xx_yz_yy, g_0_xx_yz_yz, g_0_xx_yz_zz, g_y_0_0_0_y_xx_yz_xx, g_y_0_0_0_y_xx_yz_xy, g_y_0_0_0_y_xx_yz_xz, g_y_0_0_0_y_xx_yz_yy, g_y_0_0_0_y_xx_yz_yz, g_y_0_0_0_y_xx_yz_zz, g_yy_xx_yz_xx, g_yy_xx_yz_xy, g_yy_xx_yz_xz, g_yy_xx_yz_yy, g_yy_xx_yz_yz, g_yy_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xx_yz_xx[i] = -g_0_xx_yz_xx[i] + 2.0 * g_yy_xx_yz_xx[i] * a_exp;

        g_y_0_0_0_y_xx_yz_xy[i] = -g_0_xx_yz_xy[i] + 2.0 * g_yy_xx_yz_xy[i] * a_exp;

        g_y_0_0_0_y_xx_yz_xz[i] = -g_0_xx_yz_xz[i] + 2.0 * g_yy_xx_yz_xz[i] * a_exp;

        g_y_0_0_0_y_xx_yz_yy[i] = -g_0_xx_yz_yy[i] + 2.0 * g_yy_xx_yz_yy[i] * a_exp;

        g_y_0_0_0_y_xx_yz_yz[i] = -g_0_xx_yz_yz[i] + 2.0 * g_yy_xx_yz_yz[i] * a_exp;

        g_y_0_0_0_y_xx_yz_zz[i] = -g_0_xx_yz_zz[i] + 2.0 * g_yy_xx_yz_zz[i] * a_exp;
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_0_xx_zz_xx, g_0_xx_zz_xy, g_0_xx_zz_xz, g_0_xx_zz_yy, g_0_xx_zz_yz, g_0_xx_zz_zz, g_y_0_0_0_y_xx_zz_xx, g_y_0_0_0_y_xx_zz_xy, g_y_0_0_0_y_xx_zz_xz, g_y_0_0_0_y_xx_zz_yy, g_y_0_0_0_y_xx_zz_yz, g_y_0_0_0_y_xx_zz_zz, g_yy_xx_zz_xx, g_yy_xx_zz_xy, g_yy_xx_zz_xz, g_yy_xx_zz_yy, g_yy_xx_zz_yz, g_yy_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xx_zz_xx[i] = -g_0_xx_zz_xx[i] + 2.0 * g_yy_xx_zz_xx[i] * a_exp;

        g_y_0_0_0_y_xx_zz_xy[i] = -g_0_xx_zz_xy[i] + 2.0 * g_yy_xx_zz_xy[i] * a_exp;

        g_y_0_0_0_y_xx_zz_xz[i] = -g_0_xx_zz_xz[i] + 2.0 * g_yy_xx_zz_xz[i] * a_exp;

        g_y_0_0_0_y_xx_zz_yy[i] = -g_0_xx_zz_yy[i] + 2.0 * g_yy_xx_zz_yy[i] * a_exp;

        g_y_0_0_0_y_xx_zz_yz[i] = -g_0_xx_zz_yz[i] + 2.0 * g_yy_xx_zz_yz[i] * a_exp;

        g_y_0_0_0_y_xx_zz_zz[i] = -g_0_xx_zz_zz[i] + 2.0 * g_yy_xx_zz_zz[i] * a_exp;
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_0_xy_xx_xx, g_0_xy_xx_xy, g_0_xy_xx_xz, g_0_xy_xx_yy, g_0_xy_xx_yz, g_0_xy_xx_zz, g_y_0_0_0_y_xy_xx_xx, g_y_0_0_0_y_xy_xx_xy, g_y_0_0_0_y_xy_xx_xz, g_y_0_0_0_y_xy_xx_yy, g_y_0_0_0_y_xy_xx_yz, g_y_0_0_0_y_xy_xx_zz, g_yy_xy_xx_xx, g_yy_xy_xx_xy, g_yy_xy_xx_xz, g_yy_xy_xx_yy, g_yy_xy_xx_yz, g_yy_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xy_xx_xx[i] = -g_0_xy_xx_xx[i] + 2.0 * g_yy_xy_xx_xx[i] * a_exp;

        g_y_0_0_0_y_xy_xx_xy[i] = -g_0_xy_xx_xy[i] + 2.0 * g_yy_xy_xx_xy[i] * a_exp;

        g_y_0_0_0_y_xy_xx_xz[i] = -g_0_xy_xx_xz[i] + 2.0 * g_yy_xy_xx_xz[i] * a_exp;

        g_y_0_0_0_y_xy_xx_yy[i] = -g_0_xy_xx_yy[i] + 2.0 * g_yy_xy_xx_yy[i] * a_exp;

        g_y_0_0_0_y_xy_xx_yz[i] = -g_0_xy_xx_yz[i] + 2.0 * g_yy_xy_xx_yz[i] * a_exp;

        g_y_0_0_0_y_xy_xx_zz[i] = -g_0_xy_xx_zz[i] + 2.0 * g_yy_xy_xx_zz[i] * a_exp;
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_0_xy_xy_xx, g_0_xy_xy_xy, g_0_xy_xy_xz, g_0_xy_xy_yy, g_0_xy_xy_yz, g_0_xy_xy_zz, g_y_0_0_0_y_xy_xy_xx, g_y_0_0_0_y_xy_xy_xy, g_y_0_0_0_y_xy_xy_xz, g_y_0_0_0_y_xy_xy_yy, g_y_0_0_0_y_xy_xy_yz, g_y_0_0_0_y_xy_xy_zz, g_yy_xy_xy_xx, g_yy_xy_xy_xy, g_yy_xy_xy_xz, g_yy_xy_xy_yy, g_yy_xy_xy_yz, g_yy_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xy_xy_xx[i] = -g_0_xy_xy_xx[i] + 2.0 * g_yy_xy_xy_xx[i] * a_exp;

        g_y_0_0_0_y_xy_xy_xy[i] = -g_0_xy_xy_xy[i] + 2.0 * g_yy_xy_xy_xy[i] * a_exp;

        g_y_0_0_0_y_xy_xy_xz[i] = -g_0_xy_xy_xz[i] + 2.0 * g_yy_xy_xy_xz[i] * a_exp;

        g_y_0_0_0_y_xy_xy_yy[i] = -g_0_xy_xy_yy[i] + 2.0 * g_yy_xy_xy_yy[i] * a_exp;

        g_y_0_0_0_y_xy_xy_yz[i] = -g_0_xy_xy_yz[i] + 2.0 * g_yy_xy_xy_yz[i] * a_exp;

        g_y_0_0_0_y_xy_xy_zz[i] = -g_0_xy_xy_zz[i] + 2.0 * g_yy_xy_xy_zz[i] * a_exp;
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_0_xy_xz_xx, g_0_xy_xz_xy, g_0_xy_xz_xz, g_0_xy_xz_yy, g_0_xy_xz_yz, g_0_xy_xz_zz, g_y_0_0_0_y_xy_xz_xx, g_y_0_0_0_y_xy_xz_xy, g_y_0_0_0_y_xy_xz_xz, g_y_0_0_0_y_xy_xz_yy, g_y_0_0_0_y_xy_xz_yz, g_y_0_0_0_y_xy_xz_zz, g_yy_xy_xz_xx, g_yy_xy_xz_xy, g_yy_xy_xz_xz, g_yy_xy_xz_yy, g_yy_xy_xz_yz, g_yy_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xy_xz_xx[i] = -g_0_xy_xz_xx[i] + 2.0 * g_yy_xy_xz_xx[i] * a_exp;

        g_y_0_0_0_y_xy_xz_xy[i] = -g_0_xy_xz_xy[i] + 2.0 * g_yy_xy_xz_xy[i] * a_exp;

        g_y_0_0_0_y_xy_xz_xz[i] = -g_0_xy_xz_xz[i] + 2.0 * g_yy_xy_xz_xz[i] * a_exp;

        g_y_0_0_0_y_xy_xz_yy[i] = -g_0_xy_xz_yy[i] + 2.0 * g_yy_xy_xz_yy[i] * a_exp;

        g_y_0_0_0_y_xy_xz_yz[i] = -g_0_xy_xz_yz[i] + 2.0 * g_yy_xy_xz_yz[i] * a_exp;

        g_y_0_0_0_y_xy_xz_zz[i] = -g_0_xy_xz_zz[i] + 2.0 * g_yy_xy_xz_zz[i] * a_exp;
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_0_xy_yy_xx, g_0_xy_yy_xy, g_0_xy_yy_xz, g_0_xy_yy_yy, g_0_xy_yy_yz, g_0_xy_yy_zz, g_y_0_0_0_y_xy_yy_xx, g_y_0_0_0_y_xy_yy_xy, g_y_0_0_0_y_xy_yy_xz, g_y_0_0_0_y_xy_yy_yy, g_y_0_0_0_y_xy_yy_yz, g_y_0_0_0_y_xy_yy_zz, g_yy_xy_yy_xx, g_yy_xy_yy_xy, g_yy_xy_yy_xz, g_yy_xy_yy_yy, g_yy_xy_yy_yz, g_yy_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xy_yy_xx[i] = -g_0_xy_yy_xx[i] + 2.0 * g_yy_xy_yy_xx[i] * a_exp;

        g_y_0_0_0_y_xy_yy_xy[i] = -g_0_xy_yy_xy[i] + 2.0 * g_yy_xy_yy_xy[i] * a_exp;

        g_y_0_0_0_y_xy_yy_xz[i] = -g_0_xy_yy_xz[i] + 2.0 * g_yy_xy_yy_xz[i] * a_exp;

        g_y_0_0_0_y_xy_yy_yy[i] = -g_0_xy_yy_yy[i] + 2.0 * g_yy_xy_yy_yy[i] * a_exp;

        g_y_0_0_0_y_xy_yy_yz[i] = -g_0_xy_yy_yz[i] + 2.0 * g_yy_xy_yy_yz[i] * a_exp;

        g_y_0_0_0_y_xy_yy_zz[i] = -g_0_xy_yy_zz[i] + 2.0 * g_yy_xy_yy_zz[i] * a_exp;
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_0_xy_yz_xx, g_0_xy_yz_xy, g_0_xy_yz_xz, g_0_xy_yz_yy, g_0_xy_yz_yz, g_0_xy_yz_zz, g_y_0_0_0_y_xy_yz_xx, g_y_0_0_0_y_xy_yz_xy, g_y_0_0_0_y_xy_yz_xz, g_y_0_0_0_y_xy_yz_yy, g_y_0_0_0_y_xy_yz_yz, g_y_0_0_0_y_xy_yz_zz, g_yy_xy_yz_xx, g_yy_xy_yz_xy, g_yy_xy_yz_xz, g_yy_xy_yz_yy, g_yy_xy_yz_yz, g_yy_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xy_yz_xx[i] = -g_0_xy_yz_xx[i] + 2.0 * g_yy_xy_yz_xx[i] * a_exp;

        g_y_0_0_0_y_xy_yz_xy[i] = -g_0_xy_yz_xy[i] + 2.0 * g_yy_xy_yz_xy[i] * a_exp;

        g_y_0_0_0_y_xy_yz_xz[i] = -g_0_xy_yz_xz[i] + 2.0 * g_yy_xy_yz_xz[i] * a_exp;

        g_y_0_0_0_y_xy_yz_yy[i] = -g_0_xy_yz_yy[i] + 2.0 * g_yy_xy_yz_yy[i] * a_exp;

        g_y_0_0_0_y_xy_yz_yz[i] = -g_0_xy_yz_yz[i] + 2.0 * g_yy_xy_yz_yz[i] * a_exp;

        g_y_0_0_0_y_xy_yz_zz[i] = -g_0_xy_yz_zz[i] + 2.0 * g_yy_xy_yz_zz[i] * a_exp;
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_0_xy_zz_xx, g_0_xy_zz_xy, g_0_xy_zz_xz, g_0_xy_zz_yy, g_0_xy_zz_yz, g_0_xy_zz_zz, g_y_0_0_0_y_xy_zz_xx, g_y_0_0_0_y_xy_zz_xy, g_y_0_0_0_y_xy_zz_xz, g_y_0_0_0_y_xy_zz_yy, g_y_0_0_0_y_xy_zz_yz, g_y_0_0_0_y_xy_zz_zz, g_yy_xy_zz_xx, g_yy_xy_zz_xy, g_yy_xy_zz_xz, g_yy_xy_zz_yy, g_yy_xy_zz_yz, g_yy_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xy_zz_xx[i] = -g_0_xy_zz_xx[i] + 2.0 * g_yy_xy_zz_xx[i] * a_exp;

        g_y_0_0_0_y_xy_zz_xy[i] = -g_0_xy_zz_xy[i] + 2.0 * g_yy_xy_zz_xy[i] * a_exp;

        g_y_0_0_0_y_xy_zz_xz[i] = -g_0_xy_zz_xz[i] + 2.0 * g_yy_xy_zz_xz[i] * a_exp;

        g_y_0_0_0_y_xy_zz_yy[i] = -g_0_xy_zz_yy[i] + 2.0 * g_yy_xy_zz_yy[i] * a_exp;

        g_y_0_0_0_y_xy_zz_yz[i] = -g_0_xy_zz_yz[i] + 2.0 * g_yy_xy_zz_yz[i] * a_exp;

        g_y_0_0_0_y_xy_zz_zz[i] = -g_0_xy_zz_zz[i] + 2.0 * g_yy_xy_zz_zz[i] * a_exp;
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_0_xz_xx_xx, g_0_xz_xx_xy, g_0_xz_xx_xz, g_0_xz_xx_yy, g_0_xz_xx_yz, g_0_xz_xx_zz, g_y_0_0_0_y_xz_xx_xx, g_y_0_0_0_y_xz_xx_xy, g_y_0_0_0_y_xz_xx_xz, g_y_0_0_0_y_xz_xx_yy, g_y_0_0_0_y_xz_xx_yz, g_y_0_0_0_y_xz_xx_zz, g_yy_xz_xx_xx, g_yy_xz_xx_xy, g_yy_xz_xx_xz, g_yy_xz_xx_yy, g_yy_xz_xx_yz, g_yy_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xz_xx_xx[i] = -g_0_xz_xx_xx[i] + 2.0 * g_yy_xz_xx_xx[i] * a_exp;

        g_y_0_0_0_y_xz_xx_xy[i] = -g_0_xz_xx_xy[i] + 2.0 * g_yy_xz_xx_xy[i] * a_exp;

        g_y_0_0_0_y_xz_xx_xz[i] = -g_0_xz_xx_xz[i] + 2.0 * g_yy_xz_xx_xz[i] * a_exp;

        g_y_0_0_0_y_xz_xx_yy[i] = -g_0_xz_xx_yy[i] + 2.0 * g_yy_xz_xx_yy[i] * a_exp;

        g_y_0_0_0_y_xz_xx_yz[i] = -g_0_xz_xx_yz[i] + 2.0 * g_yy_xz_xx_yz[i] * a_exp;

        g_y_0_0_0_y_xz_xx_zz[i] = -g_0_xz_xx_zz[i] + 2.0 * g_yy_xz_xx_zz[i] * a_exp;
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_0_xz_xy_xx, g_0_xz_xy_xy, g_0_xz_xy_xz, g_0_xz_xy_yy, g_0_xz_xy_yz, g_0_xz_xy_zz, g_y_0_0_0_y_xz_xy_xx, g_y_0_0_0_y_xz_xy_xy, g_y_0_0_0_y_xz_xy_xz, g_y_0_0_0_y_xz_xy_yy, g_y_0_0_0_y_xz_xy_yz, g_y_0_0_0_y_xz_xy_zz, g_yy_xz_xy_xx, g_yy_xz_xy_xy, g_yy_xz_xy_xz, g_yy_xz_xy_yy, g_yy_xz_xy_yz, g_yy_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xz_xy_xx[i] = -g_0_xz_xy_xx[i] + 2.0 * g_yy_xz_xy_xx[i] * a_exp;

        g_y_0_0_0_y_xz_xy_xy[i] = -g_0_xz_xy_xy[i] + 2.0 * g_yy_xz_xy_xy[i] * a_exp;

        g_y_0_0_0_y_xz_xy_xz[i] = -g_0_xz_xy_xz[i] + 2.0 * g_yy_xz_xy_xz[i] * a_exp;

        g_y_0_0_0_y_xz_xy_yy[i] = -g_0_xz_xy_yy[i] + 2.0 * g_yy_xz_xy_yy[i] * a_exp;

        g_y_0_0_0_y_xz_xy_yz[i] = -g_0_xz_xy_yz[i] + 2.0 * g_yy_xz_xy_yz[i] * a_exp;

        g_y_0_0_0_y_xz_xy_zz[i] = -g_0_xz_xy_zz[i] + 2.0 * g_yy_xz_xy_zz[i] * a_exp;
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_0_xz_xz_xx, g_0_xz_xz_xy, g_0_xz_xz_xz, g_0_xz_xz_yy, g_0_xz_xz_yz, g_0_xz_xz_zz, g_y_0_0_0_y_xz_xz_xx, g_y_0_0_0_y_xz_xz_xy, g_y_0_0_0_y_xz_xz_xz, g_y_0_0_0_y_xz_xz_yy, g_y_0_0_0_y_xz_xz_yz, g_y_0_0_0_y_xz_xz_zz, g_yy_xz_xz_xx, g_yy_xz_xz_xy, g_yy_xz_xz_xz, g_yy_xz_xz_yy, g_yy_xz_xz_yz, g_yy_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xz_xz_xx[i] = -g_0_xz_xz_xx[i] + 2.0 * g_yy_xz_xz_xx[i] * a_exp;

        g_y_0_0_0_y_xz_xz_xy[i] = -g_0_xz_xz_xy[i] + 2.0 * g_yy_xz_xz_xy[i] * a_exp;

        g_y_0_0_0_y_xz_xz_xz[i] = -g_0_xz_xz_xz[i] + 2.0 * g_yy_xz_xz_xz[i] * a_exp;

        g_y_0_0_0_y_xz_xz_yy[i] = -g_0_xz_xz_yy[i] + 2.0 * g_yy_xz_xz_yy[i] * a_exp;

        g_y_0_0_0_y_xz_xz_yz[i] = -g_0_xz_xz_yz[i] + 2.0 * g_yy_xz_xz_yz[i] * a_exp;

        g_y_0_0_0_y_xz_xz_zz[i] = -g_0_xz_xz_zz[i] + 2.0 * g_yy_xz_xz_zz[i] * a_exp;
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_0_xz_yy_xx, g_0_xz_yy_xy, g_0_xz_yy_xz, g_0_xz_yy_yy, g_0_xz_yy_yz, g_0_xz_yy_zz, g_y_0_0_0_y_xz_yy_xx, g_y_0_0_0_y_xz_yy_xy, g_y_0_0_0_y_xz_yy_xz, g_y_0_0_0_y_xz_yy_yy, g_y_0_0_0_y_xz_yy_yz, g_y_0_0_0_y_xz_yy_zz, g_yy_xz_yy_xx, g_yy_xz_yy_xy, g_yy_xz_yy_xz, g_yy_xz_yy_yy, g_yy_xz_yy_yz, g_yy_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xz_yy_xx[i] = -g_0_xz_yy_xx[i] + 2.0 * g_yy_xz_yy_xx[i] * a_exp;

        g_y_0_0_0_y_xz_yy_xy[i] = -g_0_xz_yy_xy[i] + 2.0 * g_yy_xz_yy_xy[i] * a_exp;

        g_y_0_0_0_y_xz_yy_xz[i] = -g_0_xz_yy_xz[i] + 2.0 * g_yy_xz_yy_xz[i] * a_exp;

        g_y_0_0_0_y_xz_yy_yy[i] = -g_0_xz_yy_yy[i] + 2.0 * g_yy_xz_yy_yy[i] * a_exp;

        g_y_0_0_0_y_xz_yy_yz[i] = -g_0_xz_yy_yz[i] + 2.0 * g_yy_xz_yy_yz[i] * a_exp;

        g_y_0_0_0_y_xz_yy_zz[i] = -g_0_xz_yy_zz[i] + 2.0 * g_yy_xz_yy_zz[i] * a_exp;
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_0_xz_yz_xx, g_0_xz_yz_xy, g_0_xz_yz_xz, g_0_xz_yz_yy, g_0_xz_yz_yz, g_0_xz_yz_zz, g_y_0_0_0_y_xz_yz_xx, g_y_0_0_0_y_xz_yz_xy, g_y_0_0_0_y_xz_yz_xz, g_y_0_0_0_y_xz_yz_yy, g_y_0_0_0_y_xz_yz_yz, g_y_0_0_0_y_xz_yz_zz, g_yy_xz_yz_xx, g_yy_xz_yz_xy, g_yy_xz_yz_xz, g_yy_xz_yz_yy, g_yy_xz_yz_yz, g_yy_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xz_yz_xx[i] = -g_0_xz_yz_xx[i] + 2.0 * g_yy_xz_yz_xx[i] * a_exp;

        g_y_0_0_0_y_xz_yz_xy[i] = -g_0_xz_yz_xy[i] + 2.0 * g_yy_xz_yz_xy[i] * a_exp;

        g_y_0_0_0_y_xz_yz_xz[i] = -g_0_xz_yz_xz[i] + 2.0 * g_yy_xz_yz_xz[i] * a_exp;

        g_y_0_0_0_y_xz_yz_yy[i] = -g_0_xz_yz_yy[i] + 2.0 * g_yy_xz_yz_yy[i] * a_exp;

        g_y_0_0_0_y_xz_yz_yz[i] = -g_0_xz_yz_yz[i] + 2.0 * g_yy_xz_yz_yz[i] * a_exp;

        g_y_0_0_0_y_xz_yz_zz[i] = -g_0_xz_yz_zz[i] + 2.0 * g_yy_xz_yz_zz[i] * a_exp;
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_0_xz_zz_xx, g_0_xz_zz_xy, g_0_xz_zz_xz, g_0_xz_zz_yy, g_0_xz_zz_yz, g_0_xz_zz_zz, g_y_0_0_0_y_xz_zz_xx, g_y_0_0_0_y_xz_zz_xy, g_y_0_0_0_y_xz_zz_xz, g_y_0_0_0_y_xz_zz_yy, g_y_0_0_0_y_xz_zz_yz, g_y_0_0_0_y_xz_zz_zz, g_yy_xz_zz_xx, g_yy_xz_zz_xy, g_yy_xz_zz_xz, g_yy_xz_zz_yy, g_yy_xz_zz_yz, g_yy_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xz_zz_xx[i] = -g_0_xz_zz_xx[i] + 2.0 * g_yy_xz_zz_xx[i] * a_exp;

        g_y_0_0_0_y_xz_zz_xy[i] = -g_0_xz_zz_xy[i] + 2.0 * g_yy_xz_zz_xy[i] * a_exp;

        g_y_0_0_0_y_xz_zz_xz[i] = -g_0_xz_zz_xz[i] + 2.0 * g_yy_xz_zz_xz[i] * a_exp;

        g_y_0_0_0_y_xz_zz_yy[i] = -g_0_xz_zz_yy[i] + 2.0 * g_yy_xz_zz_yy[i] * a_exp;

        g_y_0_0_0_y_xz_zz_yz[i] = -g_0_xz_zz_yz[i] + 2.0 * g_yy_xz_zz_yz[i] * a_exp;

        g_y_0_0_0_y_xz_zz_zz[i] = -g_0_xz_zz_zz[i] + 2.0 * g_yy_xz_zz_zz[i] * a_exp;
    }
    // integrals block (972-978)

    #pragma omp simd aligned(g_0_yy_xx_xx, g_0_yy_xx_xy, g_0_yy_xx_xz, g_0_yy_xx_yy, g_0_yy_xx_yz, g_0_yy_xx_zz, g_y_0_0_0_y_yy_xx_xx, g_y_0_0_0_y_yy_xx_xy, g_y_0_0_0_y_yy_xx_xz, g_y_0_0_0_y_yy_xx_yy, g_y_0_0_0_y_yy_xx_yz, g_y_0_0_0_y_yy_xx_zz, g_yy_yy_xx_xx, g_yy_yy_xx_xy, g_yy_yy_xx_xz, g_yy_yy_xx_yy, g_yy_yy_xx_yz, g_yy_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yy_xx_xx[i] = -g_0_yy_xx_xx[i] + 2.0 * g_yy_yy_xx_xx[i] * a_exp;

        g_y_0_0_0_y_yy_xx_xy[i] = -g_0_yy_xx_xy[i] + 2.0 * g_yy_yy_xx_xy[i] * a_exp;

        g_y_0_0_0_y_yy_xx_xz[i] = -g_0_yy_xx_xz[i] + 2.0 * g_yy_yy_xx_xz[i] * a_exp;

        g_y_0_0_0_y_yy_xx_yy[i] = -g_0_yy_xx_yy[i] + 2.0 * g_yy_yy_xx_yy[i] * a_exp;

        g_y_0_0_0_y_yy_xx_yz[i] = -g_0_yy_xx_yz[i] + 2.0 * g_yy_yy_xx_yz[i] * a_exp;

        g_y_0_0_0_y_yy_xx_zz[i] = -g_0_yy_xx_zz[i] + 2.0 * g_yy_yy_xx_zz[i] * a_exp;
    }
    // integrals block (978-984)

    #pragma omp simd aligned(g_0_yy_xy_xx, g_0_yy_xy_xy, g_0_yy_xy_xz, g_0_yy_xy_yy, g_0_yy_xy_yz, g_0_yy_xy_zz, g_y_0_0_0_y_yy_xy_xx, g_y_0_0_0_y_yy_xy_xy, g_y_0_0_0_y_yy_xy_xz, g_y_0_0_0_y_yy_xy_yy, g_y_0_0_0_y_yy_xy_yz, g_y_0_0_0_y_yy_xy_zz, g_yy_yy_xy_xx, g_yy_yy_xy_xy, g_yy_yy_xy_xz, g_yy_yy_xy_yy, g_yy_yy_xy_yz, g_yy_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yy_xy_xx[i] = -g_0_yy_xy_xx[i] + 2.0 * g_yy_yy_xy_xx[i] * a_exp;

        g_y_0_0_0_y_yy_xy_xy[i] = -g_0_yy_xy_xy[i] + 2.0 * g_yy_yy_xy_xy[i] * a_exp;

        g_y_0_0_0_y_yy_xy_xz[i] = -g_0_yy_xy_xz[i] + 2.0 * g_yy_yy_xy_xz[i] * a_exp;

        g_y_0_0_0_y_yy_xy_yy[i] = -g_0_yy_xy_yy[i] + 2.0 * g_yy_yy_xy_yy[i] * a_exp;

        g_y_0_0_0_y_yy_xy_yz[i] = -g_0_yy_xy_yz[i] + 2.0 * g_yy_yy_xy_yz[i] * a_exp;

        g_y_0_0_0_y_yy_xy_zz[i] = -g_0_yy_xy_zz[i] + 2.0 * g_yy_yy_xy_zz[i] * a_exp;
    }
    // integrals block (984-990)

    #pragma omp simd aligned(g_0_yy_xz_xx, g_0_yy_xz_xy, g_0_yy_xz_xz, g_0_yy_xz_yy, g_0_yy_xz_yz, g_0_yy_xz_zz, g_y_0_0_0_y_yy_xz_xx, g_y_0_0_0_y_yy_xz_xy, g_y_0_0_0_y_yy_xz_xz, g_y_0_0_0_y_yy_xz_yy, g_y_0_0_0_y_yy_xz_yz, g_y_0_0_0_y_yy_xz_zz, g_yy_yy_xz_xx, g_yy_yy_xz_xy, g_yy_yy_xz_xz, g_yy_yy_xz_yy, g_yy_yy_xz_yz, g_yy_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yy_xz_xx[i] = -g_0_yy_xz_xx[i] + 2.0 * g_yy_yy_xz_xx[i] * a_exp;

        g_y_0_0_0_y_yy_xz_xy[i] = -g_0_yy_xz_xy[i] + 2.0 * g_yy_yy_xz_xy[i] * a_exp;

        g_y_0_0_0_y_yy_xz_xz[i] = -g_0_yy_xz_xz[i] + 2.0 * g_yy_yy_xz_xz[i] * a_exp;

        g_y_0_0_0_y_yy_xz_yy[i] = -g_0_yy_xz_yy[i] + 2.0 * g_yy_yy_xz_yy[i] * a_exp;

        g_y_0_0_0_y_yy_xz_yz[i] = -g_0_yy_xz_yz[i] + 2.0 * g_yy_yy_xz_yz[i] * a_exp;

        g_y_0_0_0_y_yy_xz_zz[i] = -g_0_yy_xz_zz[i] + 2.0 * g_yy_yy_xz_zz[i] * a_exp;
    }
    // integrals block (990-996)

    #pragma omp simd aligned(g_0_yy_yy_xx, g_0_yy_yy_xy, g_0_yy_yy_xz, g_0_yy_yy_yy, g_0_yy_yy_yz, g_0_yy_yy_zz, g_y_0_0_0_y_yy_yy_xx, g_y_0_0_0_y_yy_yy_xy, g_y_0_0_0_y_yy_yy_xz, g_y_0_0_0_y_yy_yy_yy, g_y_0_0_0_y_yy_yy_yz, g_y_0_0_0_y_yy_yy_zz, g_yy_yy_yy_xx, g_yy_yy_yy_xy, g_yy_yy_yy_xz, g_yy_yy_yy_yy, g_yy_yy_yy_yz, g_yy_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yy_yy_xx[i] = -g_0_yy_yy_xx[i] + 2.0 * g_yy_yy_yy_xx[i] * a_exp;

        g_y_0_0_0_y_yy_yy_xy[i] = -g_0_yy_yy_xy[i] + 2.0 * g_yy_yy_yy_xy[i] * a_exp;

        g_y_0_0_0_y_yy_yy_xz[i] = -g_0_yy_yy_xz[i] + 2.0 * g_yy_yy_yy_xz[i] * a_exp;

        g_y_0_0_0_y_yy_yy_yy[i] = -g_0_yy_yy_yy[i] + 2.0 * g_yy_yy_yy_yy[i] * a_exp;

        g_y_0_0_0_y_yy_yy_yz[i] = -g_0_yy_yy_yz[i] + 2.0 * g_yy_yy_yy_yz[i] * a_exp;

        g_y_0_0_0_y_yy_yy_zz[i] = -g_0_yy_yy_zz[i] + 2.0 * g_yy_yy_yy_zz[i] * a_exp;
    }
    // integrals block (996-1002)

    #pragma omp simd aligned(g_0_yy_yz_xx, g_0_yy_yz_xy, g_0_yy_yz_xz, g_0_yy_yz_yy, g_0_yy_yz_yz, g_0_yy_yz_zz, g_y_0_0_0_y_yy_yz_xx, g_y_0_0_0_y_yy_yz_xy, g_y_0_0_0_y_yy_yz_xz, g_y_0_0_0_y_yy_yz_yy, g_y_0_0_0_y_yy_yz_yz, g_y_0_0_0_y_yy_yz_zz, g_yy_yy_yz_xx, g_yy_yy_yz_xy, g_yy_yy_yz_xz, g_yy_yy_yz_yy, g_yy_yy_yz_yz, g_yy_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yy_yz_xx[i] = -g_0_yy_yz_xx[i] + 2.0 * g_yy_yy_yz_xx[i] * a_exp;

        g_y_0_0_0_y_yy_yz_xy[i] = -g_0_yy_yz_xy[i] + 2.0 * g_yy_yy_yz_xy[i] * a_exp;

        g_y_0_0_0_y_yy_yz_xz[i] = -g_0_yy_yz_xz[i] + 2.0 * g_yy_yy_yz_xz[i] * a_exp;

        g_y_0_0_0_y_yy_yz_yy[i] = -g_0_yy_yz_yy[i] + 2.0 * g_yy_yy_yz_yy[i] * a_exp;

        g_y_0_0_0_y_yy_yz_yz[i] = -g_0_yy_yz_yz[i] + 2.0 * g_yy_yy_yz_yz[i] * a_exp;

        g_y_0_0_0_y_yy_yz_zz[i] = -g_0_yy_yz_zz[i] + 2.0 * g_yy_yy_yz_zz[i] * a_exp;
    }
    // integrals block (1002-1008)

    #pragma omp simd aligned(g_0_yy_zz_xx, g_0_yy_zz_xy, g_0_yy_zz_xz, g_0_yy_zz_yy, g_0_yy_zz_yz, g_0_yy_zz_zz, g_y_0_0_0_y_yy_zz_xx, g_y_0_0_0_y_yy_zz_xy, g_y_0_0_0_y_yy_zz_xz, g_y_0_0_0_y_yy_zz_yy, g_y_0_0_0_y_yy_zz_yz, g_y_0_0_0_y_yy_zz_zz, g_yy_yy_zz_xx, g_yy_yy_zz_xy, g_yy_yy_zz_xz, g_yy_yy_zz_yy, g_yy_yy_zz_yz, g_yy_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yy_zz_xx[i] = -g_0_yy_zz_xx[i] + 2.0 * g_yy_yy_zz_xx[i] * a_exp;

        g_y_0_0_0_y_yy_zz_xy[i] = -g_0_yy_zz_xy[i] + 2.0 * g_yy_yy_zz_xy[i] * a_exp;

        g_y_0_0_0_y_yy_zz_xz[i] = -g_0_yy_zz_xz[i] + 2.0 * g_yy_yy_zz_xz[i] * a_exp;

        g_y_0_0_0_y_yy_zz_yy[i] = -g_0_yy_zz_yy[i] + 2.0 * g_yy_yy_zz_yy[i] * a_exp;

        g_y_0_0_0_y_yy_zz_yz[i] = -g_0_yy_zz_yz[i] + 2.0 * g_yy_yy_zz_yz[i] * a_exp;

        g_y_0_0_0_y_yy_zz_zz[i] = -g_0_yy_zz_zz[i] + 2.0 * g_yy_yy_zz_zz[i] * a_exp;
    }
    // integrals block (1008-1014)

    #pragma omp simd aligned(g_0_yz_xx_xx, g_0_yz_xx_xy, g_0_yz_xx_xz, g_0_yz_xx_yy, g_0_yz_xx_yz, g_0_yz_xx_zz, g_y_0_0_0_y_yz_xx_xx, g_y_0_0_0_y_yz_xx_xy, g_y_0_0_0_y_yz_xx_xz, g_y_0_0_0_y_yz_xx_yy, g_y_0_0_0_y_yz_xx_yz, g_y_0_0_0_y_yz_xx_zz, g_yy_yz_xx_xx, g_yy_yz_xx_xy, g_yy_yz_xx_xz, g_yy_yz_xx_yy, g_yy_yz_xx_yz, g_yy_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yz_xx_xx[i] = -g_0_yz_xx_xx[i] + 2.0 * g_yy_yz_xx_xx[i] * a_exp;

        g_y_0_0_0_y_yz_xx_xy[i] = -g_0_yz_xx_xy[i] + 2.0 * g_yy_yz_xx_xy[i] * a_exp;

        g_y_0_0_0_y_yz_xx_xz[i] = -g_0_yz_xx_xz[i] + 2.0 * g_yy_yz_xx_xz[i] * a_exp;

        g_y_0_0_0_y_yz_xx_yy[i] = -g_0_yz_xx_yy[i] + 2.0 * g_yy_yz_xx_yy[i] * a_exp;

        g_y_0_0_0_y_yz_xx_yz[i] = -g_0_yz_xx_yz[i] + 2.0 * g_yy_yz_xx_yz[i] * a_exp;

        g_y_0_0_0_y_yz_xx_zz[i] = -g_0_yz_xx_zz[i] + 2.0 * g_yy_yz_xx_zz[i] * a_exp;
    }
    // integrals block (1014-1020)

    #pragma omp simd aligned(g_0_yz_xy_xx, g_0_yz_xy_xy, g_0_yz_xy_xz, g_0_yz_xy_yy, g_0_yz_xy_yz, g_0_yz_xy_zz, g_y_0_0_0_y_yz_xy_xx, g_y_0_0_0_y_yz_xy_xy, g_y_0_0_0_y_yz_xy_xz, g_y_0_0_0_y_yz_xy_yy, g_y_0_0_0_y_yz_xy_yz, g_y_0_0_0_y_yz_xy_zz, g_yy_yz_xy_xx, g_yy_yz_xy_xy, g_yy_yz_xy_xz, g_yy_yz_xy_yy, g_yy_yz_xy_yz, g_yy_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yz_xy_xx[i] = -g_0_yz_xy_xx[i] + 2.0 * g_yy_yz_xy_xx[i] * a_exp;

        g_y_0_0_0_y_yz_xy_xy[i] = -g_0_yz_xy_xy[i] + 2.0 * g_yy_yz_xy_xy[i] * a_exp;

        g_y_0_0_0_y_yz_xy_xz[i] = -g_0_yz_xy_xz[i] + 2.0 * g_yy_yz_xy_xz[i] * a_exp;

        g_y_0_0_0_y_yz_xy_yy[i] = -g_0_yz_xy_yy[i] + 2.0 * g_yy_yz_xy_yy[i] * a_exp;

        g_y_0_0_0_y_yz_xy_yz[i] = -g_0_yz_xy_yz[i] + 2.0 * g_yy_yz_xy_yz[i] * a_exp;

        g_y_0_0_0_y_yz_xy_zz[i] = -g_0_yz_xy_zz[i] + 2.0 * g_yy_yz_xy_zz[i] * a_exp;
    }
    // integrals block (1020-1026)

    #pragma omp simd aligned(g_0_yz_xz_xx, g_0_yz_xz_xy, g_0_yz_xz_xz, g_0_yz_xz_yy, g_0_yz_xz_yz, g_0_yz_xz_zz, g_y_0_0_0_y_yz_xz_xx, g_y_0_0_0_y_yz_xz_xy, g_y_0_0_0_y_yz_xz_xz, g_y_0_0_0_y_yz_xz_yy, g_y_0_0_0_y_yz_xz_yz, g_y_0_0_0_y_yz_xz_zz, g_yy_yz_xz_xx, g_yy_yz_xz_xy, g_yy_yz_xz_xz, g_yy_yz_xz_yy, g_yy_yz_xz_yz, g_yy_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yz_xz_xx[i] = -g_0_yz_xz_xx[i] + 2.0 * g_yy_yz_xz_xx[i] * a_exp;

        g_y_0_0_0_y_yz_xz_xy[i] = -g_0_yz_xz_xy[i] + 2.0 * g_yy_yz_xz_xy[i] * a_exp;

        g_y_0_0_0_y_yz_xz_xz[i] = -g_0_yz_xz_xz[i] + 2.0 * g_yy_yz_xz_xz[i] * a_exp;

        g_y_0_0_0_y_yz_xz_yy[i] = -g_0_yz_xz_yy[i] + 2.0 * g_yy_yz_xz_yy[i] * a_exp;

        g_y_0_0_0_y_yz_xz_yz[i] = -g_0_yz_xz_yz[i] + 2.0 * g_yy_yz_xz_yz[i] * a_exp;

        g_y_0_0_0_y_yz_xz_zz[i] = -g_0_yz_xz_zz[i] + 2.0 * g_yy_yz_xz_zz[i] * a_exp;
    }
    // integrals block (1026-1032)

    #pragma omp simd aligned(g_0_yz_yy_xx, g_0_yz_yy_xy, g_0_yz_yy_xz, g_0_yz_yy_yy, g_0_yz_yy_yz, g_0_yz_yy_zz, g_y_0_0_0_y_yz_yy_xx, g_y_0_0_0_y_yz_yy_xy, g_y_0_0_0_y_yz_yy_xz, g_y_0_0_0_y_yz_yy_yy, g_y_0_0_0_y_yz_yy_yz, g_y_0_0_0_y_yz_yy_zz, g_yy_yz_yy_xx, g_yy_yz_yy_xy, g_yy_yz_yy_xz, g_yy_yz_yy_yy, g_yy_yz_yy_yz, g_yy_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yz_yy_xx[i] = -g_0_yz_yy_xx[i] + 2.0 * g_yy_yz_yy_xx[i] * a_exp;

        g_y_0_0_0_y_yz_yy_xy[i] = -g_0_yz_yy_xy[i] + 2.0 * g_yy_yz_yy_xy[i] * a_exp;

        g_y_0_0_0_y_yz_yy_xz[i] = -g_0_yz_yy_xz[i] + 2.0 * g_yy_yz_yy_xz[i] * a_exp;

        g_y_0_0_0_y_yz_yy_yy[i] = -g_0_yz_yy_yy[i] + 2.0 * g_yy_yz_yy_yy[i] * a_exp;

        g_y_0_0_0_y_yz_yy_yz[i] = -g_0_yz_yy_yz[i] + 2.0 * g_yy_yz_yy_yz[i] * a_exp;

        g_y_0_0_0_y_yz_yy_zz[i] = -g_0_yz_yy_zz[i] + 2.0 * g_yy_yz_yy_zz[i] * a_exp;
    }
    // integrals block (1032-1038)

    #pragma omp simd aligned(g_0_yz_yz_xx, g_0_yz_yz_xy, g_0_yz_yz_xz, g_0_yz_yz_yy, g_0_yz_yz_yz, g_0_yz_yz_zz, g_y_0_0_0_y_yz_yz_xx, g_y_0_0_0_y_yz_yz_xy, g_y_0_0_0_y_yz_yz_xz, g_y_0_0_0_y_yz_yz_yy, g_y_0_0_0_y_yz_yz_yz, g_y_0_0_0_y_yz_yz_zz, g_yy_yz_yz_xx, g_yy_yz_yz_xy, g_yy_yz_yz_xz, g_yy_yz_yz_yy, g_yy_yz_yz_yz, g_yy_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yz_yz_xx[i] = -g_0_yz_yz_xx[i] + 2.0 * g_yy_yz_yz_xx[i] * a_exp;

        g_y_0_0_0_y_yz_yz_xy[i] = -g_0_yz_yz_xy[i] + 2.0 * g_yy_yz_yz_xy[i] * a_exp;

        g_y_0_0_0_y_yz_yz_xz[i] = -g_0_yz_yz_xz[i] + 2.0 * g_yy_yz_yz_xz[i] * a_exp;

        g_y_0_0_0_y_yz_yz_yy[i] = -g_0_yz_yz_yy[i] + 2.0 * g_yy_yz_yz_yy[i] * a_exp;

        g_y_0_0_0_y_yz_yz_yz[i] = -g_0_yz_yz_yz[i] + 2.0 * g_yy_yz_yz_yz[i] * a_exp;

        g_y_0_0_0_y_yz_yz_zz[i] = -g_0_yz_yz_zz[i] + 2.0 * g_yy_yz_yz_zz[i] * a_exp;
    }
    // integrals block (1038-1044)

    #pragma omp simd aligned(g_0_yz_zz_xx, g_0_yz_zz_xy, g_0_yz_zz_xz, g_0_yz_zz_yy, g_0_yz_zz_yz, g_0_yz_zz_zz, g_y_0_0_0_y_yz_zz_xx, g_y_0_0_0_y_yz_zz_xy, g_y_0_0_0_y_yz_zz_xz, g_y_0_0_0_y_yz_zz_yy, g_y_0_0_0_y_yz_zz_yz, g_y_0_0_0_y_yz_zz_zz, g_yy_yz_zz_xx, g_yy_yz_zz_xy, g_yy_yz_zz_xz, g_yy_yz_zz_yy, g_yy_yz_zz_yz, g_yy_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yz_zz_xx[i] = -g_0_yz_zz_xx[i] + 2.0 * g_yy_yz_zz_xx[i] * a_exp;

        g_y_0_0_0_y_yz_zz_xy[i] = -g_0_yz_zz_xy[i] + 2.0 * g_yy_yz_zz_xy[i] * a_exp;

        g_y_0_0_0_y_yz_zz_xz[i] = -g_0_yz_zz_xz[i] + 2.0 * g_yy_yz_zz_xz[i] * a_exp;

        g_y_0_0_0_y_yz_zz_yy[i] = -g_0_yz_zz_yy[i] + 2.0 * g_yy_yz_zz_yy[i] * a_exp;

        g_y_0_0_0_y_yz_zz_yz[i] = -g_0_yz_zz_yz[i] + 2.0 * g_yy_yz_zz_yz[i] * a_exp;

        g_y_0_0_0_y_yz_zz_zz[i] = -g_0_yz_zz_zz[i] + 2.0 * g_yy_yz_zz_zz[i] * a_exp;
    }
    // integrals block (1044-1050)

    #pragma omp simd aligned(g_0_zz_xx_xx, g_0_zz_xx_xy, g_0_zz_xx_xz, g_0_zz_xx_yy, g_0_zz_xx_yz, g_0_zz_xx_zz, g_y_0_0_0_y_zz_xx_xx, g_y_0_0_0_y_zz_xx_xy, g_y_0_0_0_y_zz_xx_xz, g_y_0_0_0_y_zz_xx_yy, g_y_0_0_0_y_zz_xx_yz, g_y_0_0_0_y_zz_xx_zz, g_yy_zz_xx_xx, g_yy_zz_xx_xy, g_yy_zz_xx_xz, g_yy_zz_xx_yy, g_yy_zz_xx_yz, g_yy_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_zz_xx_xx[i] = -g_0_zz_xx_xx[i] + 2.0 * g_yy_zz_xx_xx[i] * a_exp;

        g_y_0_0_0_y_zz_xx_xy[i] = -g_0_zz_xx_xy[i] + 2.0 * g_yy_zz_xx_xy[i] * a_exp;

        g_y_0_0_0_y_zz_xx_xz[i] = -g_0_zz_xx_xz[i] + 2.0 * g_yy_zz_xx_xz[i] * a_exp;

        g_y_0_0_0_y_zz_xx_yy[i] = -g_0_zz_xx_yy[i] + 2.0 * g_yy_zz_xx_yy[i] * a_exp;

        g_y_0_0_0_y_zz_xx_yz[i] = -g_0_zz_xx_yz[i] + 2.0 * g_yy_zz_xx_yz[i] * a_exp;

        g_y_0_0_0_y_zz_xx_zz[i] = -g_0_zz_xx_zz[i] + 2.0 * g_yy_zz_xx_zz[i] * a_exp;
    }
    // integrals block (1050-1056)

    #pragma omp simd aligned(g_0_zz_xy_xx, g_0_zz_xy_xy, g_0_zz_xy_xz, g_0_zz_xy_yy, g_0_zz_xy_yz, g_0_zz_xy_zz, g_y_0_0_0_y_zz_xy_xx, g_y_0_0_0_y_zz_xy_xy, g_y_0_0_0_y_zz_xy_xz, g_y_0_0_0_y_zz_xy_yy, g_y_0_0_0_y_zz_xy_yz, g_y_0_0_0_y_zz_xy_zz, g_yy_zz_xy_xx, g_yy_zz_xy_xy, g_yy_zz_xy_xz, g_yy_zz_xy_yy, g_yy_zz_xy_yz, g_yy_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_zz_xy_xx[i] = -g_0_zz_xy_xx[i] + 2.0 * g_yy_zz_xy_xx[i] * a_exp;

        g_y_0_0_0_y_zz_xy_xy[i] = -g_0_zz_xy_xy[i] + 2.0 * g_yy_zz_xy_xy[i] * a_exp;

        g_y_0_0_0_y_zz_xy_xz[i] = -g_0_zz_xy_xz[i] + 2.0 * g_yy_zz_xy_xz[i] * a_exp;

        g_y_0_0_0_y_zz_xy_yy[i] = -g_0_zz_xy_yy[i] + 2.0 * g_yy_zz_xy_yy[i] * a_exp;

        g_y_0_0_0_y_zz_xy_yz[i] = -g_0_zz_xy_yz[i] + 2.0 * g_yy_zz_xy_yz[i] * a_exp;

        g_y_0_0_0_y_zz_xy_zz[i] = -g_0_zz_xy_zz[i] + 2.0 * g_yy_zz_xy_zz[i] * a_exp;
    }
    // integrals block (1056-1062)

    #pragma omp simd aligned(g_0_zz_xz_xx, g_0_zz_xz_xy, g_0_zz_xz_xz, g_0_zz_xz_yy, g_0_zz_xz_yz, g_0_zz_xz_zz, g_y_0_0_0_y_zz_xz_xx, g_y_0_0_0_y_zz_xz_xy, g_y_0_0_0_y_zz_xz_xz, g_y_0_0_0_y_zz_xz_yy, g_y_0_0_0_y_zz_xz_yz, g_y_0_0_0_y_zz_xz_zz, g_yy_zz_xz_xx, g_yy_zz_xz_xy, g_yy_zz_xz_xz, g_yy_zz_xz_yy, g_yy_zz_xz_yz, g_yy_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_zz_xz_xx[i] = -g_0_zz_xz_xx[i] + 2.0 * g_yy_zz_xz_xx[i] * a_exp;

        g_y_0_0_0_y_zz_xz_xy[i] = -g_0_zz_xz_xy[i] + 2.0 * g_yy_zz_xz_xy[i] * a_exp;

        g_y_0_0_0_y_zz_xz_xz[i] = -g_0_zz_xz_xz[i] + 2.0 * g_yy_zz_xz_xz[i] * a_exp;

        g_y_0_0_0_y_zz_xz_yy[i] = -g_0_zz_xz_yy[i] + 2.0 * g_yy_zz_xz_yy[i] * a_exp;

        g_y_0_0_0_y_zz_xz_yz[i] = -g_0_zz_xz_yz[i] + 2.0 * g_yy_zz_xz_yz[i] * a_exp;

        g_y_0_0_0_y_zz_xz_zz[i] = -g_0_zz_xz_zz[i] + 2.0 * g_yy_zz_xz_zz[i] * a_exp;
    }
    // integrals block (1062-1068)

    #pragma omp simd aligned(g_0_zz_yy_xx, g_0_zz_yy_xy, g_0_zz_yy_xz, g_0_zz_yy_yy, g_0_zz_yy_yz, g_0_zz_yy_zz, g_y_0_0_0_y_zz_yy_xx, g_y_0_0_0_y_zz_yy_xy, g_y_0_0_0_y_zz_yy_xz, g_y_0_0_0_y_zz_yy_yy, g_y_0_0_0_y_zz_yy_yz, g_y_0_0_0_y_zz_yy_zz, g_yy_zz_yy_xx, g_yy_zz_yy_xy, g_yy_zz_yy_xz, g_yy_zz_yy_yy, g_yy_zz_yy_yz, g_yy_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_zz_yy_xx[i] = -g_0_zz_yy_xx[i] + 2.0 * g_yy_zz_yy_xx[i] * a_exp;

        g_y_0_0_0_y_zz_yy_xy[i] = -g_0_zz_yy_xy[i] + 2.0 * g_yy_zz_yy_xy[i] * a_exp;

        g_y_0_0_0_y_zz_yy_xz[i] = -g_0_zz_yy_xz[i] + 2.0 * g_yy_zz_yy_xz[i] * a_exp;

        g_y_0_0_0_y_zz_yy_yy[i] = -g_0_zz_yy_yy[i] + 2.0 * g_yy_zz_yy_yy[i] * a_exp;

        g_y_0_0_0_y_zz_yy_yz[i] = -g_0_zz_yy_yz[i] + 2.0 * g_yy_zz_yy_yz[i] * a_exp;

        g_y_0_0_0_y_zz_yy_zz[i] = -g_0_zz_yy_zz[i] + 2.0 * g_yy_zz_yy_zz[i] * a_exp;
    }
    // integrals block (1068-1074)

    #pragma omp simd aligned(g_0_zz_yz_xx, g_0_zz_yz_xy, g_0_zz_yz_xz, g_0_zz_yz_yy, g_0_zz_yz_yz, g_0_zz_yz_zz, g_y_0_0_0_y_zz_yz_xx, g_y_0_0_0_y_zz_yz_xy, g_y_0_0_0_y_zz_yz_xz, g_y_0_0_0_y_zz_yz_yy, g_y_0_0_0_y_zz_yz_yz, g_y_0_0_0_y_zz_yz_zz, g_yy_zz_yz_xx, g_yy_zz_yz_xy, g_yy_zz_yz_xz, g_yy_zz_yz_yy, g_yy_zz_yz_yz, g_yy_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_zz_yz_xx[i] = -g_0_zz_yz_xx[i] + 2.0 * g_yy_zz_yz_xx[i] * a_exp;

        g_y_0_0_0_y_zz_yz_xy[i] = -g_0_zz_yz_xy[i] + 2.0 * g_yy_zz_yz_xy[i] * a_exp;

        g_y_0_0_0_y_zz_yz_xz[i] = -g_0_zz_yz_xz[i] + 2.0 * g_yy_zz_yz_xz[i] * a_exp;

        g_y_0_0_0_y_zz_yz_yy[i] = -g_0_zz_yz_yy[i] + 2.0 * g_yy_zz_yz_yy[i] * a_exp;

        g_y_0_0_0_y_zz_yz_yz[i] = -g_0_zz_yz_yz[i] + 2.0 * g_yy_zz_yz_yz[i] * a_exp;

        g_y_0_0_0_y_zz_yz_zz[i] = -g_0_zz_yz_zz[i] + 2.0 * g_yy_zz_yz_zz[i] * a_exp;
    }
    // integrals block (1074-1080)

    #pragma omp simd aligned(g_0_zz_zz_xx, g_0_zz_zz_xy, g_0_zz_zz_xz, g_0_zz_zz_yy, g_0_zz_zz_yz, g_0_zz_zz_zz, g_y_0_0_0_y_zz_zz_xx, g_y_0_0_0_y_zz_zz_xy, g_y_0_0_0_y_zz_zz_xz, g_y_0_0_0_y_zz_zz_yy, g_y_0_0_0_y_zz_zz_yz, g_y_0_0_0_y_zz_zz_zz, g_yy_zz_zz_xx, g_yy_zz_zz_xy, g_yy_zz_zz_xz, g_yy_zz_zz_yy, g_yy_zz_zz_yz, g_yy_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_zz_zz_xx[i] = -g_0_zz_zz_xx[i] + 2.0 * g_yy_zz_zz_xx[i] * a_exp;

        g_y_0_0_0_y_zz_zz_xy[i] = -g_0_zz_zz_xy[i] + 2.0 * g_yy_zz_zz_xy[i] * a_exp;

        g_y_0_0_0_y_zz_zz_xz[i] = -g_0_zz_zz_xz[i] + 2.0 * g_yy_zz_zz_xz[i] * a_exp;

        g_y_0_0_0_y_zz_zz_yy[i] = -g_0_zz_zz_yy[i] + 2.0 * g_yy_zz_zz_yy[i] * a_exp;

        g_y_0_0_0_y_zz_zz_yz[i] = -g_0_zz_zz_yz[i] + 2.0 * g_yy_zz_zz_yz[i] * a_exp;

        g_y_0_0_0_y_zz_zz_zz[i] = -g_0_zz_zz_zz[i] + 2.0 * g_yy_zz_zz_zz[i] * a_exp;
    }
    // integrals block (1080-1086)

    #pragma omp simd aligned(g_y_0_0_0_z_xx_xx_xx, g_y_0_0_0_z_xx_xx_xy, g_y_0_0_0_z_xx_xx_xz, g_y_0_0_0_z_xx_xx_yy, g_y_0_0_0_z_xx_xx_yz, g_y_0_0_0_z_xx_xx_zz, g_yz_xx_xx_xx, g_yz_xx_xx_xy, g_yz_xx_xx_xz, g_yz_xx_xx_yy, g_yz_xx_xx_yz, g_yz_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xx_xx_xx[i] = 2.0 * g_yz_xx_xx_xx[i] * a_exp;

        g_y_0_0_0_z_xx_xx_xy[i] = 2.0 * g_yz_xx_xx_xy[i] * a_exp;

        g_y_0_0_0_z_xx_xx_xz[i] = 2.0 * g_yz_xx_xx_xz[i] * a_exp;

        g_y_0_0_0_z_xx_xx_yy[i] = 2.0 * g_yz_xx_xx_yy[i] * a_exp;

        g_y_0_0_0_z_xx_xx_yz[i] = 2.0 * g_yz_xx_xx_yz[i] * a_exp;

        g_y_0_0_0_z_xx_xx_zz[i] = 2.0 * g_yz_xx_xx_zz[i] * a_exp;
    }
    // integrals block (1086-1092)

    #pragma omp simd aligned(g_y_0_0_0_z_xx_xy_xx, g_y_0_0_0_z_xx_xy_xy, g_y_0_0_0_z_xx_xy_xz, g_y_0_0_0_z_xx_xy_yy, g_y_0_0_0_z_xx_xy_yz, g_y_0_0_0_z_xx_xy_zz, g_yz_xx_xy_xx, g_yz_xx_xy_xy, g_yz_xx_xy_xz, g_yz_xx_xy_yy, g_yz_xx_xy_yz, g_yz_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xx_xy_xx[i] = 2.0 * g_yz_xx_xy_xx[i] * a_exp;

        g_y_0_0_0_z_xx_xy_xy[i] = 2.0 * g_yz_xx_xy_xy[i] * a_exp;

        g_y_0_0_0_z_xx_xy_xz[i] = 2.0 * g_yz_xx_xy_xz[i] * a_exp;

        g_y_0_0_0_z_xx_xy_yy[i] = 2.0 * g_yz_xx_xy_yy[i] * a_exp;

        g_y_0_0_0_z_xx_xy_yz[i] = 2.0 * g_yz_xx_xy_yz[i] * a_exp;

        g_y_0_0_0_z_xx_xy_zz[i] = 2.0 * g_yz_xx_xy_zz[i] * a_exp;
    }
    // integrals block (1092-1098)

    #pragma omp simd aligned(g_y_0_0_0_z_xx_xz_xx, g_y_0_0_0_z_xx_xz_xy, g_y_0_0_0_z_xx_xz_xz, g_y_0_0_0_z_xx_xz_yy, g_y_0_0_0_z_xx_xz_yz, g_y_0_0_0_z_xx_xz_zz, g_yz_xx_xz_xx, g_yz_xx_xz_xy, g_yz_xx_xz_xz, g_yz_xx_xz_yy, g_yz_xx_xz_yz, g_yz_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xx_xz_xx[i] = 2.0 * g_yz_xx_xz_xx[i] * a_exp;

        g_y_0_0_0_z_xx_xz_xy[i] = 2.0 * g_yz_xx_xz_xy[i] * a_exp;

        g_y_0_0_0_z_xx_xz_xz[i] = 2.0 * g_yz_xx_xz_xz[i] * a_exp;

        g_y_0_0_0_z_xx_xz_yy[i] = 2.0 * g_yz_xx_xz_yy[i] * a_exp;

        g_y_0_0_0_z_xx_xz_yz[i] = 2.0 * g_yz_xx_xz_yz[i] * a_exp;

        g_y_0_0_0_z_xx_xz_zz[i] = 2.0 * g_yz_xx_xz_zz[i] * a_exp;
    }
    // integrals block (1098-1104)

    #pragma omp simd aligned(g_y_0_0_0_z_xx_yy_xx, g_y_0_0_0_z_xx_yy_xy, g_y_0_0_0_z_xx_yy_xz, g_y_0_0_0_z_xx_yy_yy, g_y_0_0_0_z_xx_yy_yz, g_y_0_0_0_z_xx_yy_zz, g_yz_xx_yy_xx, g_yz_xx_yy_xy, g_yz_xx_yy_xz, g_yz_xx_yy_yy, g_yz_xx_yy_yz, g_yz_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xx_yy_xx[i] = 2.0 * g_yz_xx_yy_xx[i] * a_exp;

        g_y_0_0_0_z_xx_yy_xy[i] = 2.0 * g_yz_xx_yy_xy[i] * a_exp;

        g_y_0_0_0_z_xx_yy_xz[i] = 2.0 * g_yz_xx_yy_xz[i] * a_exp;

        g_y_0_0_0_z_xx_yy_yy[i] = 2.0 * g_yz_xx_yy_yy[i] * a_exp;

        g_y_0_0_0_z_xx_yy_yz[i] = 2.0 * g_yz_xx_yy_yz[i] * a_exp;

        g_y_0_0_0_z_xx_yy_zz[i] = 2.0 * g_yz_xx_yy_zz[i] * a_exp;
    }
    // integrals block (1104-1110)

    #pragma omp simd aligned(g_y_0_0_0_z_xx_yz_xx, g_y_0_0_0_z_xx_yz_xy, g_y_0_0_0_z_xx_yz_xz, g_y_0_0_0_z_xx_yz_yy, g_y_0_0_0_z_xx_yz_yz, g_y_0_0_0_z_xx_yz_zz, g_yz_xx_yz_xx, g_yz_xx_yz_xy, g_yz_xx_yz_xz, g_yz_xx_yz_yy, g_yz_xx_yz_yz, g_yz_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xx_yz_xx[i] = 2.0 * g_yz_xx_yz_xx[i] * a_exp;

        g_y_0_0_0_z_xx_yz_xy[i] = 2.0 * g_yz_xx_yz_xy[i] * a_exp;

        g_y_0_0_0_z_xx_yz_xz[i] = 2.0 * g_yz_xx_yz_xz[i] * a_exp;

        g_y_0_0_0_z_xx_yz_yy[i] = 2.0 * g_yz_xx_yz_yy[i] * a_exp;

        g_y_0_0_0_z_xx_yz_yz[i] = 2.0 * g_yz_xx_yz_yz[i] * a_exp;

        g_y_0_0_0_z_xx_yz_zz[i] = 2.0 * g_yz_xx_yz_zz[i] * a_exp;
    }
    // integrals block (1110-1116)

    #pragma omp simd aligned(g_y_0_0_0_z_xx_zz_xx, g_y_0_0_0_z_xx_zz_xy, g_y_0_0_0_z_xx_zz_xz, g_y_0_0_0_z_xx_zz_yy, g_y_0_0_0_z_xx_zz_yz, g_y_0_0_0_z_xx_zz_zz, g_yz_xx_zz_xx, g_yz_xx_zz_xy, g_yz_xx_zz_xz, g_yz_xx_zz_yy, g_yz_xx_zz_yz, g_yz_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xx_zz_xx[i] = 2.0 * g_yz_xx_zz_xx[i] * a_exp;

        g_y_0_0_0_z_xx_zz_xy[i] = 2.0 * g_yz_xx_zz_xy[i] * a_exp;

        g_y_0_0_0_z_xx_zz_xz[i] = 2.0 * g_yz_xx_zz_xz[i] * a_exp;

        g_y_0_0_0_z_xx_zz_yy[i] = 2.0 * g_yz_xx_zz_yy[i] * a_exp;

        g_y_0_0_0_z_xx_zz_yz[i] = 2.0 * g_yz_xx_zz_yz[i] * a_exp;

        g_y_0_0_0_z_xx_zz_zz[i] = 2.0 * g_yz_xx_zz_zz[i] * a_exp;
    }
    // integrals block (1116-1122)

    #pragma omp simd aligned(g_y_0_0_0_z_xy_xx_xx, g_y_0_0_0_z_xy_xx_xy, g_y_0_0_0_z_xy_xx_xz, g_y_0_0_0_z_xy_xx_yy, g_y_0_0_0_z_xy_xx_yz, g_y_0_0_0_z_xy_xx_zz, g_yz_xy_xx_xx, g_yz_xy_xx_xy, g_yz_xy_xx_xz, g_yz_xy_xx_yy, g_yz_xy_xx_yz, g_yz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xy_xx_xx[i] = 2.0 * g_yz_xy_xx_xx[i] * a_exp;

        g_y_0_0_0_z_xy_xx_xy[i] = 2.0 * g_yz_xy_xx_xy[i] * a_exp;

        g_y_0_0_0_z_xy_xx_xz[i] = 2.0 * g_yz_xy_xx_xz[i] * a_exp;

        g_y_0_0_0_z_xy_xx_yy[i] = 2.0 * g_yz_xy_xx_yy[i] * a_exp;

        g_y_0_0_0_z_xy_xx_yz[i] = 2.0 * g_yz_xy_xx_yz[i] * a_exp;

        g_y_0_0_0_z_xy_xx_zz[i] = 2.0 * g_yz_xy_xx_zz[i] * a_exp;
    }
    // integrals block (1122-1128)

    #pragma omp simd aligned(g_y_0_0_0_z_xy_xy_xx, g_y_0_0_0_z_xy_xy_xy, g_y_0_0_0_z_xy_xy_xz, g_y_0_0_0_z_xy_xy_yy, g_y_0_0_0_z_xy_xy_yz, g_y_0_0_0_z_xy_xy_zz, g_yz_xy_xy_xx, g_yz_xy_xy_xy, g_yz_xy_xy_xz, g_yz_xy_xy_yy, g_yz_xy_xy_yz, g_yz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xy_xy_xx[i] = 2.0 * g_yz_xy_xy_xx[i] * a_exp;

        g_y_0_0_0_z_xy_xy_xy[i] = 2.0 * g_yz_xy_xy_xy[i] * a_exp;

        g_y_0_0_0_z_xy_xy_xz[i] = 2.0 * g_yz_xy_xy_xz[i] * a_exp;

        g_y_0_0_0_z_xy_xy_yy[i] = 2.0 * g_yz_xy_xy_yy[i] * a_exp;

        g_y_0_0_0_z_xy_xy_yz[i] = 2.0 * g_yz_xy_xy_yz[i] * a_exp;

        g_y_0_0_0_z_xy_xy_zz[i] = 2.0 * g_yz_xy_xy_zz[i] * a_exp;
    }
    // integrals block (1128-1134)

    #pragma omp simd aligned(g_y_0_0_0_z_xy_xz_xx, g_y_0_0_0_z_xy_xz_xy, g_y_0_0_0_z_xy_xz_xz, g_y_0_0_0_z_xy_xz_yy, g_y_0_0_0_z_xy_xz_yz, g_y_0_0_0_z_xy_xz_zz, g_yz_xy_xz_xx, g_yz_xy_xz_xy, g_yz_xy_xz_xz, g_yz_xy_xz_yy, g_yz_xy_xz_yz, g_yz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xy_xz_xx[i] = 2.0 * g_yz_xy_xz_xx[i] * a_exp;

        g_y_0_0_0_z_xy_xz_xy[i] = 2.0 * g_yz_xy_xz_xy[i] * a_exp;

        g_y_0_0_0_z_xy_xz_xz[i] = 2.0 * g_yz_xy_xz_xz[i] * a_exp;

        g_y_0_0_0_z_xy_xz_yy[i] = 2.0 * g_yz_xy_xz_yy[i] * a_exp;

        g_y_0_0_0_z_xy_xz_yz[i] = 2.0 * g_yz_xy_xz_yz[i] * a_exp;

        g_y_0_0_0_z_xy_xz_zz[i] = 2.0 * g_yz_xy_xz_zz[i] * a_exp;
    }
    // integrals block (1134-1140)

    #pragma omp simd aligned(g_y_0_0_0_z_xy_yy_xx, g_y_0_0_0_z_xy_yy_xy, g_y_0_0_0_z_xy_yy_xz, g_y_0_0_0_z_xy_yy_yy, g_y_0_0_0_z_xy_yy_yz, g_y_0_0_0_z_xy_yy_zz, g_yz_xy_yy_xx, g_yz_xy_yy_xy, g_yz_xy_yy_xz, g_yz_xy_yy_yy, g_yz_xy_yy_yz, g_yz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xy_yy_xx[i] = 2.0 * g_yz_xy_yy_xx[i] * a_exp;

        g_y_0_0_0_z_xy_yy_xy[i] = 2.0 * g_yz_xy_yy_xy[i] * a_exp;

        g_y_0_0_0_z_xy_yy_xz[i] = 2.0 * g_yz_xy_yy_xz[i] * a_exp;

        g_y_0_0_0_z_xy_yy_yy[i] = 2.0 * g_yz_xy_yy_yy[i] * a_exp;

        g_y_0_0_0_z_xy_yy_yz[i] = 2.0 * g_yz_xy_yy_yz[i] * a_exp;

        g_y_0_0_0_z_xy_yy_zz[i] = 2.0 * g_yz_xy_yy_zz[i] * a_exp;
    }
    // integrals block (1140-1146)

    #pragma omp simd aligned(g_y_0_0_0_z_xy_yz_xx, g_y_0_0_0_z_xy_yz_xy, g_y_0_0_0_z_xy_yz_xz, g_y_0_0_0_z_xy_yz_yy, g_y_0_0_0_z_xy_yz_yz, g_y_0_0_0_z_xy_yz_zz, g_yz_xy_yz_xx, g_yz_xy_yz_xy, g_yz_xy_yz_xz, g_yz_xy_yz_yy, g_yz_xy_yz_yz, g_yz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xy_yz_xx[i] = 2.0 * g_yz_xy_yz_xx[i] * a_exp;

        g_y_0_0_0_z_xy_yz_xy[i] = 2.0 * g_yz_xy_yz_xy[i] * a_exp;

        g_y_0_0_0_z_xy_yz_xz[i] = 2.0 * g_yz_xy_yz_xz[i] * a_exp;

        g_y_0_0_0_z_xy_yz_yy[i] = 2.0 * g_yz_xy_yz_yy[i] * a_exp;

        g_y_0_0_0_z_xy_yz_yz[i] = 2.0 * g_yz_xy_yz_yz[i] * a_exp;

        g_y_0_0_0_z_xy_yz_zz[i] = 2.0 * g_yz_xy_yz_zz[i] * a_exp;
    }
    // integrals block (1146-1152)

    #pragma omp simd aligned(g_y_0_0_0_z_xy_zz_xx, g_y_0_0_0_z_xy_zz_xy, g_y_0_0_0_z_xy_zz_xz, g_y_0_0_0_z_xy_zz_yy, g_y_0_0_0_z_xy_zz_yz, g_y_0_0_0_z_xy_zz_zz, g_yz_xy_zz_xx, g_yz_xy_zz_xy, g_yz_xy_zz_xz, g_yz_xy_zz_yy, g_yz_xy_zz_yz, g_yz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xy_zz_xx[i] = 2.0 * g_yz_xy_zz_xx[i] * a_exp;

        g_y_0_0_0_z_xy_zz_xy[i] = 2.0 * g_yz_xy_zz_xy[i] * a_exp;

        g_y_0_0_0_z_xy_zz_xz[i] = 2.0 * g_yz_xy_zz_xz[i] * a_exp;

        g_y_0_0_0_z_xy_zz_yy[i] = 2.0 * g_yz_xy_zz_yy[i] * a_exp;

        g_y_0_0_0_z_xy_zz_yz[i] = 2.0 * g_yz_xy_zz_yz[i] * a_exp;

        g_y_0_0_0_z_xy_zz_zz[i] = 2.0 * g_yz_xy_zz_zz[i] * a_exp;
    }
    // integrals block (1152-1158)

    #pragma omp simd aligned(g_y_0_0_0_z_xz_xx_xx, g_y_0_0_0_z_xz_xx_xy, g_y_0_0_0_z_xz_xx_xz, g_y_0_0_0_z_xz_xx_yy, g_y_0_0_0_z_xz_xx_yz, g_y_0_0_0_z_xz_xx_zz, g_yz_xz_xx_xx, g_yz_xz_xx_xy, g_yz_xz_xx_xz, g_yz_xz_xx_yy, g_yz_xz_xx_yz, g_yz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xz_xx_xx[i] = 2.0 * g_yz_xz_xx_xx[i] * a_exp;

        g_y_0_0_0_z_xz_xx_xy[i] = 2.0 * g_yz_xz_xx_xy[i] * a_exp;

        g_y_0_0_0_z_xz_xx_xz[i] = 2.0 * g_yz_xz_xx_xz[i] * a_exp;

        g_y_0_0_0_z_xz_xx_yy[i] = 2.0 * g_yz_xz_xx_yy[i] * a_exp;

        g_y_0_0_0_z_xz_xx_yz[i] = 2.0 * g_yz_xz_xx_yz[i] * a_exp;

        g_y_0_0_0_z_xz_xx_zz[i] = 2.0 * g_yz_xz_xx_zz[i] * a_exp;
    }
    // integrals block (1158-1164)

    #pragma omp simd aligned(g_y_0_0_0_z_xz_xy_xx, g_y_0_0_0_z_xz_xy_xy, g_y_0_0_0_z_xz_xy_xz, g_y_0_0_0_z_xz_xy_yy, g_y_0_0_0_z_xz_xy_yz, g_y_0_0_0_z_xz_xy_zz, g_yz_xz_xy_xx, g_yz_xz_xy_xy, g_yz_xz_xy_xz, g_yz_xz_xy_yy, g_yz_xz_xy_yz, g_yz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xz_xy_xx[i] = 2.0 * g_yz_xz_xy_xx[i] * a_exp;

        g_y_0_0_0_z_xz_xy_xy[i] = 2.0 * g_yz_xz_xy_xy[i] * a_exp;

        g_y_0_0_0_z_xz_xy_xz[i] = 2.0 * g_yz_xz_xy_xz[i] * a_exp;

        g_y_0_0_0_z_xz_xy_yy[i] = 2.0 * g_yz_xz_xy_yy[i] * a_exp;

        g_y_0_0_0_z_xz_xy_yz[i] = 2.0 * g_yz_xz_xy_yz[i] * a_exp;

        g_y_0_0_0_z_xz_xy_zz[i] = 2.0 * g_yz_xz_xy_zz[i] * a_exp;
    }
    // integrals block (1164-1170)

    #pragma omp simd aligned(g_y_0_0_0_z_xz_xz_xx, g_y_0_0_0_z_xz_xz_xy, g_y_0_0_0_z_xz_xz_xz, g_y_0_0_0_z_xz_xz_yy, g_y_0_0_0_z_xz_xz_yz, g_y_0_0_0_z_xz_xz_zz, g_yz_xz_xz_xx, g_yz_xz_xz_xy, g_yz_xz_xz_xz, g_yz_xz_xz_yy, g_yz_xz_xz_yz, g_yz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xz_xz_xx[i] = 2.0 * g_yz_xz_xz_xx[i] * a_exp;

        g_y_0_0_0_z_xz_xz_xy[i] = 2.0 * g_yz_xz_xz_xy[i] * a_exp;

        g_y_0_0_0_z_xz_xz_xz[i] = 2.0 * g_yz_xz_xz_xz[i] * a_exp;

        g_y_0_0_0_z_xz_xz_yy[i] = 2.0 * g_yz_xz_xz_yy[i] * a_exp;

        g_y_0_0_0_z_xz_xz_yz[i] = 2.0 * g_yz_xz_xz_yz[i] * a_exp;

        g_y_0_0_0_z_xz_xz_zz[i] = 2.0 * g_yz_xz_xz_zz[i] * a_exp;
    }
    // integrals block (1170-1176)

    #pragma omp simd aligned(g_y_0_0_0_z_xz_yy_xx, g_y_0_0_0_z_xz_yy_xy, g_y_0_0_0_z_xz_yy_xz, g_y_0_0_0_z_xz_yy_yy, g_y_0_0_0_z_xz_yy_yz, g_y_0_0_0_z_xz_yy_zz, g_yz_xz_yy_xx, g_yz_xz_yy_xy, g_yz_xz_yy_xz, g_yz_xz_yy_yy, g_yz_xz_yy_yz, g_yz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xz_yy_xx[i] = 2.0 * g_yz_xz_yy_xx[i] * a_exp;

        g_y_0_0_0_z_xz_yy_xy[i] = 2.0 * g_yz_xz_yy_xy[i] * a_exp;

        g_y_0_0_0_z_xz_yy_xz[i] = 2.0 * g_yz_xz_yy_xz[i] * a_exp;

        g_y_0_0_0_z_xz_yy_yy[i] = 2.0 * g_yz_xz_yy_yy[i] * a_exp;

        g_y_0_0_0_z_xz_yy_yz[i] = 2.0 * g_yz_xz_yy_yz[i] * a_exp;

        g_y_0_0_0_z_xz_yy_zz[i] = 2.0 * g_yz_xz_yy_zz[i] * a_exp;
    }
    // integrals block (1176-1182)

    #pragma omp simd aligned(g_y_0_0_0_z_xz_yz_xx, g_y_0_0_0_z_xz_yz_xy, g_y_0_0_0_z_xz_yz_xz, g_y_0_0_0_z_xz_yz_yy, g_y_0_0_0_z_xz_yz_yz, g_y_0_0_0_z_xz_yz_zz, g_yz_xz_yz_xx, g_yz_xz_yz_xy, g_yz_xz_yz_xz, g_yz_xz_yz_yy, g_yz_xz_yz_yz, g_yz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xz_yz_xx[i] = 2.0 * g_yz_xz_yz_xx[i] * a_exp;

        g_y_0_0_0_z_xz_yz_xy[i] = 2.0 * g_yz_xz_yz_xy[i] * a_exp;

        g_y_0_0_0_z_xz_yz_xz[i] = 2.0 * g_yz_xz_yz_xz[i] * a_exp;

        g_y_0_0_0_z_xz_yz_yy[i] = 2.0 * g_yz_xz_yz_yy[i] * a_exp;

        g_y_0_0_0_z_xz_yz_yz[i] = 2.0 * g_yz_xz_yz_yz[i] * a_exp;

        g_y_0_0_0_z_xz_yz_zz[i] = 2.0 * g_yz_xz_yz_zz[i] * a_exp;
    }
    // integrals block (1182-1188)

    #pragma omp simd aligned(g_y_0_0_0_z_xz_zz_xx, g_y_0_0_0_z_xz_zz_xy, g_y_0_0_0_z_xz_zz_xz, g_y_0_0_0_z_xz_zz_yy, g_y_0_0_0_z_xz_zz_yz, g_y_0_0_0_z_xz_zz_zz, g_yz_xz_zz_xx, g_yz_xz_zz_xy, g_yz_xz_zz_xz, g_yz_xz_zz_yy, g_yz_xz_zz_yz, g_yz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xz_zz_xx[i] = 2.0 * g_yz_xz_zz_xx[i] * a_exp;

        g_y_0_0_0_z_xz_zz_xy[i] = 2.0 * g_yz_xz_zz_xy[i] * a_exp;

        g_y_0_0_0_z_xz_zz_xz[i] = 2.0 * g_yz_xz_zz_xz[i] * a_exp;

        g_y_0_0_0_z_xz_zz_yy[i] = 2.0 * g_yz_xz_zz_yy[i] * a_exp;

        g_y_0_0_0_z_xz_zz_yz[i] = 2.0 * g_yz_xz_zz_yz[i] * a_exp;

        g_y_0_0_0_z_xz_zz_zz[i] = 2.0 * g_yz_xz_zz_zz[i] * a_exp;
    }
    // integrals block (1188-1194)

    #pragma omp simd aligned(g_y_0_0_0_z_yy_xx_xx, g_y_0_0_0_z_yy_xx_xy, g_y_0_0_0_z_yy_xx_xz, g_y_0_0_0_z_yy_xx_yy, g_y_0_0_0_z_yy_xx_yz, g_y_0_0_0_z_yy_xx_zz, g_yz_yy_xx_xx, g_yz_yy_xx_xy, g_yz_yy_xx_xz, g_yz_yy_xx_yy, g_yz_yy_xx_yz, g_yz_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yy_xx_xx[i] = 2.0 * g_yz_yy_xx_xx[i] * a_exp;

        g_y_0_0_0_z_yy_xx_xy[i] = 2.0 * g_yz_yy_xx_xy[i] * a_exp;

        g_y_0_0_0_z_yy_xx_xz[i] = 2.0 * g_yz_yy_xx_xz[i] * a_exp;

        g_y_0_0_0_z_yy_xx_yy[i] = 2.0 * g_yz_yy_xx_yy[i] * a_exp;

        g_y_0_0_0_z_yy_xx_yz[i] = 2.0 * g_yz_yy_xx_yz[i] * a_exp;

        g_y_0_0_0_z_yy_xx_zz[i] = 2.0 * g_yz_yy_xx_zz[i] * a_exp;
    }
    // integrals block (1194-1200)

    #pragma omp simd aligned(g_y_0_0_0_z_yy_xy_xx, g_y_0_0_0_z_yy_xy_xy, g_y_0_0_0_z_yy_xy_xz, g_y_0_0_0_z_yy_xy_yy, g_y_0_0_0_z_yy_xy_yz, g_y_0_0_0_z_yy_xy_zz, g_yz_yy_xy_xx, g_yz_yy_xy_xy, g_yz_yy_xy_xz, g_yz_yy_xy_yy, g_yz_yy_xy_yz, g_yz_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yy_xy_xx[i] = 2.0 * g_yz_yy_xy_xx[i] * a_exp;

        g_y_0_0_0_z_yy_xy_xy[i] = 2.0 * g_yz_yy_xy_xy[i] * a_exp;

        g_y_0_0_0_z_yy_xy_xz[i] = 2.0 * g_yz_yy_xy_xz[i] * a_exp;

        g_y_0_0_0_z_yy_xy_yy[i] = 2.0 * g_yz_yy_xy_yy[i] * a_exp;

        g_y_0_0_0_z_yy_xy_yz[i] = 2.0 * g_yz_yy_xy_yz[i] * a_exp;

        g_y_0_0_0_z_yy_xy_zz[i] = 2.0 * g_yz_yy_xy_zz[i] * a_exp;
    }
    // integrals block (1200-1206)

    #pragma omp simd aligned(g_y_0_0_0_z_yy_xz_xx, g_y_0_0_0_z_yy_xz_xy, g_y_0_0_0_z_yy_xz_xz, g_y_0_0_0_z_yy_xz_yy, g_y_0_0_0_z_yy_xz_yz, g_y_0_0_0_z_yy_xz_zz, g_yz_yy_xz_xx, g_yz_yy_xz_xy, g_yz_yy_xz_xz, g_yz_yy_xz_yy, g_yz_yy_xz_yz, g_yz_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yy_xz_xx[i] = 2.0 * g_yz_yy_xz_xx[i] * a_exp;

        g_y_0_0_0_z_yy_xz_xy[i] = 2.0 * g_yz_yy_xz_xy[i] * a_exp;

        g_y_0_0_0_z_yy_xz_xz[i] = 2.0 * g_yz_yy_xz_xz[i] * a_exp;

        g_y_0_0_0_z_yy_xz_yy[i] = 2.0 * g_yz_yy_xz_yy[i] * a_exp;

        g_y_0_0_0_z_yy_xz_yz[i] = 2.0 * g_yz_yy_xz_yz[i] * a_exp;

        g_y_0_0_0_z_yy_xz_zz[i] = 2.0 * g_yz_yy_xz_zz[i] * a_exp;
    }
    // integrals block (1206-1212)

    #pragma omp simd aligned(g_y_0_0_0_z_yy_yy_xx, g_y_0_0_0_z_yy_yy_xy, g_y_0_0_0_z_yy_yy_xz, g_y_0_0_0_z_yy_yy_yy, g_y_0_0_0_z_yy_yy_yz, g_y_0_0_0_z_yy_yy_zz, g_yz_yy_yy_xx, g_yz_yy_yy_xy, g_yz_yy_yy_xz, g_yz_yy_yy_yy, g_yz_yy_yy_yz, g_yz_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yy_yy_xx[i] = 2.0 * g_yz_yy_yy_xx[i] * a_exp;

        g_y_0_0_0_z_yy_yy_xy[i] = 2.0 * g_yz_yy_yy_xy[i] * a_exp;

        g_y_0_0_0_z_yy_yy_xz[i] = 2.0 * g_yz_yy_yy_xz[i] * a_exp;

        g_y_0_0_0_z_yy_yy_yy[i] = 2.0 * g_yz_yy_yy_yy[i] * a_exp;

        g_y_0_0_0_z_yy_yy_yz[i] = 2.0 * g_yz_yy_yy_yz[i] * a_exp;

        g_y_0_0_0_z_yy_yy_zz[i] = 2.0 * g_yz_yy_yy_zz[i] * a_exp;
    }
    // integrals block (1212-1218)

    #pragma omp simd aligned(g_y_0_0_0_z_yy_yz_xx, g_y_0_0_0_z_yy_yz_xy, g_y_0_0_0_z_yy_yz_xz, g_y_0_0_0_z_yy_yz_yy, g_y_0_0_0_z_yy_yz_yz, g_y_0_0_0_z_yy_yz_zz, g_yz_yy_yz_xx, g_yz_yy_yz_xy, g_yz_yy_yz_xz, g_yz_yy_yz_yy, g_yz_yy_yz_yz, g_yz_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yy_yz_xx[i] = 2.0 * g_yz_yy_yz_xx[i] * a_exp;

        g_y_0_0_0_z_yy_yz_xy[i] = 2.0 * g_yz_yy_yz_xy[i] * a_exp;

        g_y_0_0_0_z_yy_yz_xz[i] = 2.0 * g_yz_yy_yz_xz[i] * a_exp;

        g_y_0_0_0_z_yy_yz_yy[i] = 2.0 * g_yz_yy_yz_yy[i] * a_exp;

        g_y_0_0_0_z_yy_yz_yz[i] = 2.0 * g_yz_yy_yz_yz[i] * a_exp;

        g_y_0_0_0_z_yy_yz_zz[i] = 2.0 * g_yz_yy_yz_zz[i] * a_exp;
    }
    // integrals block (1218-1224)

    #pragma omp simd aligned(g_y_0_0_0_z_yy_zz_xx, g_y_0_0_0_z_yy_zz_xy, g_y_0_0_0_z_yy_zz_xz, g_y_0_0_0_z_yy_zz_yy, g_y_0_0_0_z_yy_zz_yz, g_y_0_0_0_z_yy_zz_zz, g_yz_yy_zz_xx, g_yz_yy_zz_xy, g_yz_yy_zz_xz, g_yz_yy_zz_yy, g_yz_yy_zz_yz, g_yz_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yy_zz_xx[i] = 2.0 * g_yz_yy_zz_xx[i] * a_exp;

        g_y_0_0_0_z_yy_zz_xy[i] = 2.0 * g_yz_yy_zz_xy[i] * a_exp;

        g_y_0_0_0_z_yy_zz_xz[i] = 2.0 * g_yz_yy_zz_xz[i] * a_exp;

        g_y_0_0_0_z_yy_zz_yy[i] = 2.0 * g_yz_yy_zz_yy[i] * a_exp;

        g_y_0_0_0_z_yy_zz_yz[i] = 2.0 * g_yz_yy_zz_yz[i] * a_exp;

        g_y_0_0_0_z_yy_zz_zz[i] = 2.0 * g_yz_yy_zz_zz[i] * a_exp;
    }
    // integrals block (1224-1230)

    #pragma omp simd aligned(g_y_0_0_0_z_yz_xx_xx, g_y_0_0_0_z_yz_xx_xy, g_y_0_0_0_z_yz_xx_xz, g_y_0_0_0_z_yz_xx_yy, g_y_0_0_0_z_yz_xx_yz, g_y_0_0_0_z_yz_xx_zz, g_yz_yz_xx_xx, g_yz_yz_xx_xy, g_yz_yz_xx_xz, g_yz_yz_xx_yy, g_yz_yz_xx_yz, g_yz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yz_xx_xx[i] = 2.0 * g_yz_yz_xx_xx[i] * a_exp;

        g_y_0_0_0_z_yz_xx_xy[i] = 2.0 * g_yz_yz_xx_xy[i] * a_exp;

        g_y_0_0_0_z_yz_xx_xz[i] = 2.0 * g_yz_yz_xx_xz[i] * a_exp;

        g_y_0_0_0_z_yz_xx_yy[i] = 2.0 * g_yz_yz_xx_yy[i] * a_exp;

        g_y_0_0_0_z_yz_xx_yz[i] = 2.0 * g_yz_yz_xx_yz[i] * a_exp;

        g_y_0_0_0_z_yz_xx_zz[i] = 2.0 * g_yz_yz_xx_zz[i] * a_exp;
    }
    // integrals block (1230-1236)

    #pragma omp simd aligned(g_y_0_0_0_z_yz_xy_xx, g_y_0_0_0_z_yz_xy_xy, g_y_0_0_0_z_yz_xy_xz, g_y_0_0_0_z_yz_xy_yy, g_y_0_0_0_z_yz_xy_yz, g_y_0_0_0_z_yz_xy_zz, g_yz_yz_xy_xx, g_yz_yz_xy_xy, g_yz_yz_xy_xz, g_yz_yz_xy_yy, g_yz_yz_xy_yz, g_yz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yz_xy_xx[i] = 2.0 * g_yz_yz_xy_xx[i] * a_exp;

        g_y_0_0_0_z_yz_xy_xy[i] = 2.0 * g_yz_yz_xy_xy[i] * a_exp;

        g_y_0_0_0_z_yz_xy_xz[i] = 2.0 * g_yz_yz_xy_xz[i] * a_exp;

        g_y_0_0_0_z_yz_xy_yy[i] = 2.0 * g_yz_yz_xy_yy[i] * a_exp;

        g_y_0_0_0_z_yz_xy_yz[i] = 2.0 * g_yz_yz_xy_yz[i] * a_exp;

        g_y_0_0_0_z_yz_xy_zz[i] = 2.0 * g_yz_yz_xy_zz[i] * a_exp;
    }
    // integrals block (1236-1242)

    #pragma omp simd aligned(g_y_0_0_0_z_yz_xz_xx, g_y_0_0_0_z_yz_xz_xy, g_y_0_0_0_z_yz_xz_xz, g_y_0_0_0_z_yz_xz_yy, g_y_0_0_0_z_yz_xz_yz, g_y_0_0_0_z_yz_xz_zz, g_yz_yz_xz_xx, g_yz_yz_xz_xy, g_yz_yz_xz_xz, g_yz_yz_xz_yy, g_yz_yz_xz_yz, g_yz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yz_xz_xx[i] = 2.0 * g_yz_yz_xz_xx[i] * a_exp;

        g_y_0_0_0_z_yz_xz_xy[i] = 2.0 * g_yz_yz_xz_xy[i] * a_exp;

        g_y_0_0_0_z_yz_xz_xz[i] = 2.0 * g_yz_yz_xz_xz[i] * a_exp;

        g_y_0_0_0_z_yz_xz_yy[i] = 2.0 * g_yz_yz_xz_yy[i] * a_exp;

        g_y_0_0_0_z_yz_xz_yz[i] = 2.0 * g_yz_yz_xz_yz[i] * a_exp;

        g_y_0_0_0_z_yz_xz_zz[i] = 2.0 * g_yz_yz_xz_zz[i] * a_exp;
    }
    // integrals block (1242-1248)

    #pragma omp simd aligned(g_y_0_0_0_z_yz_yy_xx, g_y_0_0_0_z_yz_yy_xy, g_y_0_0_0_z_yz_yy_xz, g_y_0_0_0_z_yz_yy_yy, g_y_0_0_0_z_yz_yy_yz, g_y_0_0_0_z_yz_yy_zz, g_yz_yz_yy_xx, g_yz_yz_yy_xy, g_yz_yz_yy_xz, g_yz_yz_yy_yy, g_yz_yz_yy_yz, g_yz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yz_yy_xx[i] = 2.0 * g_yz_yz_yy_xx[i] * a_exp;

        g_y_0_0_0_z_yz_yy_xy[i] = 2.0 * g_yz_yz_yy_xy[i] * a_exp;

        g_y_0_0_0_z_yz_yy_xz[i] = 2.0 * g_yz_yz_yy_xz[i] * a_exp;

        g_y_0_0_0_z_yz_yy_yy[i] = 2.0 * g_yz_yz_yy_yy[i] * a_exp;

        g_y_0_0_0_z_yz_yy_yz[i] = 2.0 * g_yz_yz_yy_yz[i] * a_exp;

        g_y_0_0_0_z_yz_yy_zz[i] = 2.0 * g_yz_yz_yy_zz[i] * a_exp;
    }
    // integrals block (1248-1254)

    #pragma omp simd aligned(g_y_0_0_0_z_yz_yz_xx, g_y_0_0_0_z_yz_yz_xy, g_y_0_0_0_z_yz_yz_xz, g_y_0_0_0_z_yz_yz_yy, g_y_0_0_0_z_yz_yz_yz, g_y_0_0_0_z_yz_yz_zz, g_yz_yz_yz_xx, g_yz_yz_yz_xy, g_yz_yz_yz_xz, g_yz_yz_yz_yy, g_yz_yz_yz_yz, g_yz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yz_yz_xx[i] = 2.0 * g_yz_yz_yz_xx[i] * a_exp;

        g_y_0_0_0_z_yz_yz_xy[i] = 2.0 * g_yz_yz_yz_xy[i] * a_exp;

        g_y_0_0_0_z_yz_yz_xz[i] = 2.0 * g_yz_yz_yz_xz[i] * a_exp;

        g_y_0_0_0_z_yz_yz_yy[i] = 2.0 * g_yz_yz_yz_yy[i] * a_exp;

        g_y_0_0_0_z_yz_yz_yz[i] = 2.0 * g_yz_yz_yz_yz[i] * a_exp;

        g_y_0_0_0_z_yz_yz_zz[i] = 2.0 * g_yz_yz_yz_zz[i] * a_exp;
    }
    // integrals block (1254-1260)

    #pragma omp simd aligned(g_y_0_0_0_z_yz_zz_xx, g_y_0_0_0_z_yz_zz_xy, g_y_0_0_0_z_yz_zz_xz, g_y_0_0_0_z_yz_zz_yy, g_y_0_0_0_z_yz_zz_yz, g_y_0_0_0_z_yz_zz_zz, g_yz_yz_zz_xx, g_yz_yz_zz_xy, g_yz_yz_zz_xz, g_yz_yz_zz_yy, g_yz_yz_zz_yz, g_yz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yz_zz_xx[i] = 2.0 * g_yz_yz_zz_xx[i] * a_exp;

        g_y_0_0_0_z_yz_zz_xy[i] = 2.0 * g_yz_yz_zz_xy[i] * a_exp;

        g_y_0_0_0_z_yz_zz_xz[i] = 2.0 * g_yz_yz_zz_xz[i] * a_exp;

        g_y_0_0_0_z_yz_zz_yy[i] = 2.0 * g_yz_yz_zz_yy[i] * a_exp;

        g_y_0_0_0_z_yz_zz_yz[i] = 2.0 * g_yz_yz_zz_yz[i] * a_exp;

        g_y_0_0_0_z_yz_zz_zz[i] = 2.0 * g_yz_yz_zz_zz[i] * a_exp;
    }
    // integrals block (1260-1266)

    #pragma omp simd aligned(g_y_0_0_0_z_zz_xx_xx, g_y_0_0_0_z_zz_xx_xy, g_y_0_0_0_z_zz_xx_xz, g_y_0_0_0_z_zz_xx_yy, g_y_0_0_0_z_zz_xx_yz, g_y_0_0_0_z_zz_xx_zz, g_yz_zz_xx_xx, g_yz_zz_xx_xy, g_yz_zz_xx_xz, g_yz_zz_xx_yy, g_yz_zz_xx_yz, g_yz_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_zz_xx_xx[i] = 2.0 * g_yz_zz_xx_xx[i] * a_exp;

        g_y_0_0_0_z_zz_xx_xy[i] = 2.0 * g_yz_zz_xx_xy[i] * a_exp;

        g_y_0_0_0_z_zz_xx_xz[i] = 2.0 * g_yz_zz_xx_xz[i] * a_exp;

        g_y_0_0_0_z_zz_xx_yy[i] = 2.0 * g_yz_zz_xx_yy[i] * a_exp;

        g_y_0_0_0_z_zz_xx_yz[i] = 2.0 * g_yz_zz_xx_yz[i] * a_exp;

        g_y_0_0_0_z_zz_xx_zz[i] = 2.0 * g_yz_zz_xx_zz[i] * a_exp;
    }
    // integrals block (1266-1272)

    #pragma omp simd aligned(g_y_0_0_0_z_zz_xy_xx, g_y_0_0_0_z_zz_xy_xy, g_y_0_0_0_z_zz_xy_xz, g_y_0_0_0_z_zz_xy_yy, g_y_0_0_0_z_zz_xy_yz, g_y_0_0_0_z_zz_xy_zz, g_yz_zz_xy_xx, g_yz_zz_xy_xy, g_yz_zz_xy_xz, g_yz_zz_xy_yy, g_yz_zz_xy_yz, g_yz_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_zz_xy_xx[i] = 2.0 * g_yz_zz_xy_xx[i] * a_exp;

        g_y_0_0_0_z_zz_xy_xy[i] = 2.0 * g_yz_zz_xy_xy[i] * a_exp;

        g_y_0_0_0_z_zz_xy_xz[i] = 2.0 * g_yz_zz_xy_xz[i] * a_exp;

        g_y_0_0_0_z_zz_xy_yy[i] = 2.0 * g_yz_zz_xy_yy[i] * a_exp;

        g_y_0_0_0_z_zz_xy_yz[i] = 2.0 * g_yz_zz_xy_yz[i] * a_exp;

        g_y_0_0_0_z_zz_xy_zz[i] = 2.0 * g_yz_zz_xy_zz[i] * a_exp;
    }
    // integrals block (1272-1278)

    #pragma omp simd aligned(g_y_0_0_0_z_zz_xz_xx, g_y_0_0_0_z_zz_xz_xy, g_y_0_0_0_z_zz_xz_xz, g_y_0_0_0_z_zz_xz_yy, g_y_0_0_0_z_zz_xz_yz, g_y_0_0_0_z_zz_xz_zz, g_yz_zz_xz_xx, g_yz_zz_xz_xy, g_yz_zz_xz_xz, g_yz_zz_xz_yy, g_yz_zz_xz_yz, g_yz_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_zz_xz_xx[i] = 2.0 * g_yz_zz_xz_xx[i] * a_exp;

        g_y_0_0_0_z_zz_xz_xy[i] = 2.0 * g_yz_zz_xz_xy[i] * a_exp;

        g_y_0_0_0_z_zz_xz_xz[i] = 2.0 * g_yz_zz_xz_xz[i] * a_exp;

        g_y_0_0_0_z_zz_xz_yy[i] = 2.0 * g_yz_zz_xz_yy[i] * a_exp;

        g_y_0_0_0_z_zz_xz_yz[i] = 2.0 * g_yz_zz_xz_yz[i] * a_exp;

        g_y_0_0_0_z_zz_xz_zz[i] = 2.0 * g_yz_zz_xz_zz[i] * a_exp;
    }
    // integrals block (1278-1284)

    #pragma omp simd aligned(g_y_0_0_0_z_zz_yy_xx, g_y_0_0_0_z_zz_yy_xy, g_y_0_0_0_z_zz_yy_xz, g_y_0_0_0_z_zz_yy_yy, g_y_0_0_0_z_zz_yy_yz, g_y_0_0_0_z_zz_yy_zz, g_yz_zz_yy_xx, g_yz_zz_yy_xy, g_yz_zz_yy_xz, g_yz_zz_yy_yy, g_yz_zz_yy_yz, g_yz_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_zz_yy_xx[i] = 2.0 * g_yz_zz_yy_xx[i] * a_exp;

        g_y_0_0_0_z_zz_yy_xy[i] = 2.0 * g_yz_zz_yy_xy[i] * a_exp;

        g_y_0_0_0_z_zz_yy_xz[i] = 2.0 * g_yz_zz_yy_xz[i] * a_exp;

        g_y_0_0_0_z_zz_yy_yy[i] = 2.0 * g_yz_zz_yy_yy[i] * a_exp;

        g_y_0_0_0_z_zz_yy_yz[i] = 2.0 * g_yz_zz_yy_yz[i] * a_exp;

        g_y_0_0_0_z_zz_yy_zz[i] = 2.0 * g_yz_zz_yy_zz[i] * a_exp;
    }
    // integrals block (1284-1290)

    #pragma omp simd aligned(g_y_0_0_0_z_zz_yz_xx, g_y_0_0_0_z_zz_yz_xy, g_y_0_0_0_z_zz_yz_xz, g_y_0_0_0_z_zz_yz_yy, g_y_0_0_0_z_zz_yz_yz, g_y_0_0_0_z_zz_yz_zz, g_yz_zz_yz_xx, g_yz_zz_yz_xy, g_yz_zz_yz_xz, g_yz_zz_yz_yy, g_yz_zz_yz_yz, g_yz_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_zz_yz_xx[i] = 2.0 * g_yz_zz_yz_xx[i] * a_exp;

        g_y_0_0_0_z_zz_yz_xy[i] = 2.0 * g_yz_zz_yz_xy[i] * a_exp;

        g_y_0_0_0_z_zz_yz_xz[i] = 2.0 * g_yz_zz_yz_xz[i] * a_exp;

        g_y_0_0_0_z_zz_yz_yy[i] = 2.0 * g_yz_zz_yz_yy[i] * a_exp;

        g_y_0_0_0_z_zz_yz_yz[i] = 2.0 * g_yz_zz_yz_yz[i] * a_exp;

        g_y_0_0_0_z_zz_yz_zz[i] = 2.0 * g_yz_zz_yz_zz[i] * a_exp;
    }
    // integrals block (1290-1296)

    #pragma omp simd aligned(g_y_0_0_0_z_zz_zz_xx, g_y_0_0_0_z_zz_zz_xy, g_y_0_0_0_z_zz_zz_xz, g_y_0_0_0_z_zz_zz_yy, g_y_0_0_0_z_zz_zz_yz, g_y_0_0_0_z_zz_zz_zz, g_yz_zz_zz_xx, g_yz_zz_zz_xy, g_yz_zz_zz_xz, g_yz_zz_zz_yy, g_yz_zz_zz_yz, g_yz_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_zz_zz_xx[i] = 2.0 * g_yz_zz_zz_xx[i] * a_exp;

        g_y_0_0_0_z_zz_zz_xy[i] = 2.0 * g_yz_zz_zz_xy[i] * a_exp;

        g_y_0_0_0_z_zz_zz_xz[i] = 2.0 * g_yz_zz_zz_xz[i] * a_exp;

        g_y_0_0_0_z_zz_zz_yy[i] = 2.0 * g_yz_zz_zz_yy[i] * a_exp;

        g_y_0_0_0_z_zz_zz_yz[i] = 2.0 * g_yz_zz_zz_yz[i] * a_exp;

        g_y_0_0_0_z_zz_zz_zz[i] = 2.0 * g_yz_zz_zz_zz[i] * a_exp;
    }
    // integrals block (1296-1302)

    #pragma omp simd aligned(g_xz_xx_xx_xx, g_xz_xx_xx_xy, g_xz_xx_xx_xz, g_xz_xx_xx_yy, g_xz_xx_xx_yz, g_xz_xx_xx_zz, g_z_0_0_0_x_xx_xx_xx, g_z_0_0_0_x_xx_xx_xy, g_z_0_0_0_x_xx_xx_xz, g_z_0_0_0_x_xx_xx_yy, g_z_0_0_0_x_xx_xx_yz, g_z_0_0_0_x_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xx_xx_xx[i] = 2.0 * g_xz_xx_xx_xx[i] * a_exp;

        g_z_0_0_0_x_xx_xx_xy[i] = 2.0 * g_xz_xx_xx_xy[i] * a_exp;

        g_z_0_0_0_x_xx_xx_xz[i] = 2.0 * g_xz_xx_xx_xz[i] * a_exp;

        g_z_0_0_0_x_xx_xx_yy[i] = 2.0 * g_xz_xx_xx_yy[i] * a_exp;

        g_z_0_0_0_x_xx_xx_yz[i] = 2.0 * g_xz_xx_xx_yz[i] * a_exp;

        g_z_0_0_0_x_xx_xx_zz[i] = 2.0 * g_xz_xx_xx_zz[i] * a_exp;
    }
    // integrals block (1302-1308)

    #pragma omp simd aligned(g_xz_xx_xy_xx, g_xz_xx_xy_xy, g_xz_xx_xy_xz, g_xz_xx_xy_yy, g_xz_xx_xy_yz, g_xz_xx_xy_zz, g_z_0_0_0_x_xx_xy_xx, g_z_0_0_0_x_xx_xy_xy, g_z_0_0_0_x_xx_xy_xz, g_z_0_0_0_x_xx_xy_yy, g_z_0_0_0_x_xx_xy_yz, g_z_0_0_0_x_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xx_xy_xx[i] = 2.0 * g_xz_xx_xy_xx[i] * a_exp;

        g_z_0_0_0_x_xx_xy_xy[i] = 2.0 * g_xz_xx_xy_xy[i] * a_exp;

        g_z_0_0_0_x_xx_xy_xz[i] = 2.0 * g_xz_xx_xy_xz[i] * a_exp;

        g_z_0_0_0_x_xx_xy_yy[i] = 2.0 * g_xz_xx_xy_yy[i] * a_exp;

        g_z_0_0_0_x_xx_xy_yz[i] = 2.0 * g_xz_xx_xy_yz[i] * a_exp;

        g_z_0_0_0_x_xx_xy_zz[i] = 2.0 * g_xz_xx_xy_zz[i] * a_exp;
    }
    // integrals block (1308-1314)

    #pragma omp simd aligned(g_xz_xx_xz_xx, g_xz_xx_xz_xy, g_xz_xx_xz_xz, g_xz_xx_xz_yy, g_xz_xx_xz_yz, g_xz_xx_xz_zz, g_z_0_0_0_x_xx_xz_xx, g_z_0_0_0_x_xx_xz_xy, g_z_0_0_0_x_xx_xz_xz, g_z_0_0_0_x_xx_xz_yy, g_z_0_0_0_x_xx_xz_yz, g_z_0_0_0_x_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xx_xz_xx[i] = 2.0 * g_xz_xx_xz_xx[i] * a_exp;

        g_z_0_0_0_x_xx_xz_xy[i] = 2.0 * g_xz_xx_xz_xy[i] * a_exp;

        g_z_0_0_0_x_xx_xz_xz[i] = 2.0 * g_xz_xx_xz_xz[i] * a_exp;

        g_z_0_0_0_x_xx_xz_yy[i] = 2.0 * g_xz_xx_xz_yy[i] * a_exp;

        g_z_0_0_0_x_xx_xz_yz[i] = 2.0 * g_xz_xx_xz_yz[i] * a_exp;

        g_z_0_0_0_x_xx_xz_zz[i] = 2.0 * g_xz_xx_xz_zz[i] * a_exp;
    }
    // integrals block (1314-1320)

    #pragma omp simd aligned(g_xz_xx_yy_xx, g_xz_xx_yy_xy, g_xz_xx_yy_xz, g_xz_xx_yy_yy, g_xz_xx_yy_yz, g_xz_xx_yy_zz, g_z_0_0_0_x_xx_yy_xx, g_z_0_0_0_x_xx_yy_xy, g_z_0_0_0_x_xx_yy_xz, g_z_0_0_0_x_xx_yy_yy, g_z_0_0_0_x_xx_yy_yz, g_z_0_0_0_x_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xx_yy_xx[i] = 2.0 * g_xz_xx_yy_xx[i] * a_exp;

        g_z_0_0_0_x_xx_yy_xy[i] = 2.0 * g_xz_xx_yy_xy[i] * a_exp;

        g_z_0_0_0_x_xx_yy_xz[i] = 2.0 * g_xz_xx_yy_xz[i] * a_exp;

        g_z_0_0_0_x_xx_yy_yy[i] = 2.0 * g_xz_xx_yy_yy[i] * a_exp;

        g_z_0_0_0_x_xx_yy_yz[i] = 2.0 * g_xz_xx_yy_yz[i] * a_exp;

        g_z_0_0_0_x_xx_yy_zz[i] = 2.0 * g_xz_xx_yy_zz[i] * a_exp;
    }
    // integrals block (1320-1326)

    #pragma omp simd aligned(g_xz_xx_yz_xx, g_xz_xx_yz_xy, g_xz_xx_yz_xz, g_xz_xx_yz_yy, g_xz_xx_yz_yz, g_xz_xx_yz_zz, g_z_0_0_0_x_xx_yz_xx, g_z_0_0_0_x_xx_yz_xy, g_z_0_0_0_x_xx_yz_xz, g_z_0_0_0_x_xx_yz_yy, g_z_0_0_0_x_xx_yz_yz, g_z_0_0_0_x_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xx_yz_xx[i] = 2.0 * g_xz_xx_yz_xx[i] * a_exp;

        g_z_0_0_0_x_xx_yz_xy[i] = 2.0 * g_xz_xx_yz_xy[i] * a_exp;

        g_z_0_0_0_x_xx_yz_xz[i] = 2.0 * g_xz_xx_yz_xz[i] * a_exp;

        g_z_0_0_0_x_xx_yz_yy[i] = 2.0 * g_xz_xx_yz_yy[i] * a_exp;

        g_z_0_0_0_x_xx_yz_yz[i] = 2.0 * g_xz_xx_yz_yz[i] * a_exp;

        g_z_0_0_0_x_xx_yz_zz[i] = 2.0 * g_xz_xx_yz_zz[i] * a_exp;
    }
    // integrals block (1326-1332)

    #pragma omp simd aligned(g_xz_xx_zz_xx, g_xz_xx_zz_xy, g_xz_xx_zz_xz, g_xz_xx_zz_yy, g_xz_xx_zz_yz, g_xz_xx_zz_zz, g_z_0_0_0_x_xx_zz_xx, g_z_0_0_0_x_xx_zz_xy, g_z_0_0_0_x_xx_zz_xz, g_z_0_0_0_x_xx_zz_yy, g_z_0_0_0_x_xx_zz_yz, g_z_0_0_0_x_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xx_zz_xx[i] = 2.0 * g_xz_xx_zz_xx[i] * a_exp;

        g_z_0_0_0_x_xx_zz_xy[i] = 2.0 * g_xz_xx_zz_xy[i] * a_exp;

        g_z_0_0_0_x_xx_zz_xz[i] = 2.0 * g_xz_xx_zz_xz[i] * a_exp;

        g_z_0_0_0_x_xx_zz_yy[i] = 2.0 * g_xz_xx_zz_yy[i] * a_exp;

        g_z_0_0_0_x_xx_zz_yz[i] = 2.0 * g_xz_xx_zz_yz[i] * a_exp;

        g_z_0_0_0_x_xx_zz_zz[i] = 2.0 * g_xz_xx_zz_zz[i] * a_exp;
    }
    // integrals block (1332-1338)

    #pragma omp simd aligned(g_xz_xy_xx_xx, g_xz_xy_xx_xy, g_xz_xy_xx_xz, g_xz_xy_xx_yy, g_xz_xy_xx_yz, g_xz_xy_xx_zz, g_z_0_0_0_x_xy_xx_xx, g_z_0_0_0_x_xy_xx_xy, g_z_0_0_0_x_xy_xx_xz, g_z_0_0_0_x_xy_xx_yy, g_z_0_0_0_x_xy_xx_yz, g_z_0_0_0_x_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xy_xx_xx[i] = 2.0 * g_xz_xy_xx_xx[i] * a_exp;

        g_z_0_0_0_x_xy_xx_xy[i] = 2.0 * g_xz_xy_xx_xy[i] * a_exp;

        g_z_0_0_0_x_xy_xx_xz[i] = 2.0 * g_xz_xy_xx_xz[i] * a_exp;

        g_z_0_0_0_x_xy_xx_yy[i] = 2.0 * g_xz_xy_xx_yy[i] * a_exp;

        g_z_0_0_0_x_xy_xx_yz[i] = 2.0 * g_xz_xy_xx_yz[i] * a_exp;

        g_z_0_0_0_x_xy_xx_zz[i] = 2.0 * g_xz_xy_xx_zz[i] * a_exp;
    }
    // integrals block (1338-1344)

    #pragma omp simd aligned(g_xz_xy_xy_xx, g_xz_xy_xy_xy, g_xz_xy_xy_xz, g_xz_xy_xy_yy, g_xz_xy_xy_yz, g_xz_xy_xy_zz, g_z_0_0_0_x_xy_xy_xx, g_z_0_0_0_x_xy_xy_xy, g_z_0_0_0_x_xy_xy_xz, g_z_0_0_0_x_xy_xy_yy, g_z_0_0_0_x_xy_xy_yz, g_z_0_0_0_x_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xy_xy_xx[i] = 2.0 * g_xz_xy_xy_xx[i] * a_exp;

        g_z_0_0_0_x_xy_xy_xy[i] = 2.0 * g_xz_xy_xy_xy[i] * a_exp;

        g_z_0_0_0_x_xy_xy_xz[i] = 2.0 * g_xz_xy_xy_xz[i] * a_exp;

        g_z_0_0_0_x_xy_xy_yy[i] = 2.0 * g_xz_xy_xy_yy[i] * a_exp;

        g_z_0_0_0_x_xy_xy_yz[i] = 2.0 * g_xz_xy_xy_yz[i] * a_exp;

        g_z_0_0_0_x_xy_xy_zz[i] = 2.0 * g_xz_xy_xy_zz[i] * a_exp;
    }
    // integrals block (1344-1350)

    #pragma omp simd aligned(g_xz_xy_xz_xx, g_xz_xy_xz_xy, g_xz_xy_xz_xz, g_xz_xy_xz_yy, g_xz_xy_xz_yz, g_xz_xy_xz_zz, g_z_0_0_0_x_xy_xz_xx, g_z_0_0_0_x_xy_xz_xy, g_z_0_0_0_x_xy_xz_xz, g_z_0_0_0_x_xy_xz_yy, g_z_0_0_0_x_xy_xz_yz, g_z_0_0_0_x_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xy_xz_xx[i] = 2.0 * g_xz_xy_xz_xx[i] * a_exp;

        g_z_0_0_0_x_xy_xz_xy[i] = 2.0 * g_xz_xy_xz_xy[i] * a_exp;

        g_z_0_0_0_x_xy_xz_xz[i] = 2.0 * g_xz_xy_xz_xz[i] * a_exp;

        g_z_0_0_0_x_xy_xz_yy[i] = 2.0 * g_xz_xy_xz_yy[i] * a_exp;

        g_z_0_0_0_x_xy_xz_yz[i] = 2.0 * g_xz_xy_xz_yz[i] * a_exp;

        g_z_0_0_0_x_xy_xz_zz[i] = 2.0 * g_xz_xy_xz_zz[i] * a_exp;
    }
    // integrals block (1350-1356)

    #pragma omp simd aligned(g_xz_xy_yy_xx, g_xz_xy_yy_xy, g_xz_xy_yy_xz, g_xz_xy_yy_yy, g_xz_xy_yy_yz, g_xz_xy_yy_zz, g_z_0_0_0_x_xy_yy_xx, g_z_0_0_0_x_xy_yy_xy, g_z_0_0_0_x_xy_yy_xz, g_z_0_0_0_x_xy_yy_yy, g_z_0_0_0_x_xy_yy_yz, g_z_0_0_0_x_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xy_yy_xx[i] = 2.0 * g_xz_xy_yy_xx[i] * a_exp;

        g_z_0_0_0_x_xy_yy_xy[i] = 2.0 * g_xz_xy_yy_xy[i] * a_exp;

        g_z_0_0_0_x_xy_yy_xz[i] = 2.0 * g_xz_xy_yy_xz[i] * a_exp;

        g_z_0_0_0_x_xy_yy_yy[i] = 2.0 * g_xz_xy_yy_yy[i] * a_exp;

        g_z_0_0_0_x_xy_yy_yz[i] = 2.0 * g_xz_xy_yy_yz[i] * a_exp;

        g_z_0_0_0_x_xy_yy_zz[i] = 2.0 * g_xz_xy_yy_zz[i] * a_exp;
    }
    // integrals block (1356-1362)

    #pragma omp simd aligned(g_xz_xy_yz_xx, g_xz_xy_yz_xy, g_xz_xy_yz_xz, g_xz_xy_yz_yy, g_xz_xy_yz_yz, g_xz_xy_yz_zz, g_z_0_0_0_x_xy_yz_xx, g_z_0_0_0_x_xy_yz_xy, g_z_0_0_0_x_xy_yz_xz, g_z_0_0_0_x_xy_yz_yy, g_z_0_0_0_x_xy_yz_yz, g_z_0_0_0_x_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xy_yz_xx[i] = 2.0 * g_xz_xy_yz_xx[i] * a_exp;

        g_z_0_0_0_x_xy_yz_xy[i] = 2.0 * g_xz_xy_yz_xy[i] * a_exp;

        g_z_0_0_0_x_xy_yz_xz[i] = 2.0 * g_xz_xy_yz_xz[i] * a_exp;

        g_z_0_0_0_x_xy_yz_yy[i] = 2.0 * g_xz_xy_yz_yy[i] * a_exp;

        g_z_0_0_0_x_xy_yz_yz[i] = 2.0 * g_xz_xy_yz_yz[i] * a_exp;

        g_z_0_0_0_x_xy_yz_zz[i] = 2.0 * g_xz_xy_yz_zz[i] * a_exp;
    }
    // integrals block (1362-1368)

    #pragma omp simd aligned(g_xz_xy_zz_xx, g_xz_xy_zz_xy, g_xz_xy_zz_xz, g_xz_xy_zz_yy, g_xz_xy_zz_yz, g_xz_xy_zz_zz, g_z_0_0_0_x_xy_zz_xx, g_z_0_0_0_x_xy_zz_xy, g_z_0_0_0_x_xy_zz_xz, g_z_0_0_0_x_xy_zz_yy, g_z_0_0_0_x_xy_zz_yz, g_z_0_0_0_x_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xy_zz_xx[i] = 2.0 * g_xz_xy_zz_xx[i] * a_exp;

        g_z_0_0_0_x_xy_zz_xy[i] = 2.0 * g_xz_xy_zz_xy[i] * a_exp;

        g_z_0_0_0_x_xy_zz_xz[i] = 2.0 * g_xz_xy_zz_xz[i] * a_exp;

        g_z_0_0_0_x_xy_zz_yy[i] = 2.0 * g_xz_xy_zz_yy[i] * a_exp;

        g_z_0_0_0_x_xy_zz_yz[i] = 2.0 * g_xz_xy_zz_yz[i] * a_exp;

        g_z_0_0_0_x_xy_zz_zz[i] = 2.0 * g_xz_xy_zz_zz[i] * a_exp;
    }
    // integrals block (1368-1374)

    #pragma omp simd aligned(g_xz_xz_xx_xx, g_xz_xz_xx_xy, g_xz_xz_xx_xz, g_xz_xz_xx_yy, g_xz_xz_xx_yz, g_xz_xz_xx_zz, g_z_0_0_0_x_xz_xx_xx, g_z_0_0_0_x_xz_xx_xy, g_z_0_0_0_x_xz_xx_xz, g_z_0_0_0_x_xz_xx_yy, g_z_0_0_0_x_xz_xx_yz, g_z_0_0_0_x_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xz_xx_xx[i] = 2.0 * g_xz_xz_xx_xx[i] * a_exp;

        g_z_0_0_0_x_xz_xx_xy[i] = 2.0 * g_xz_xz_xx_xy[i] * a_exp;

        g_z_0_0_0_x_xz_xx_xz[i] = 2.0 * g_xz_xz_xx_xz[i] * a_exp;

        g_z_0_0_0_x_xz_xx_yy[i] = 2.0 * g_xz_xz_xx_yy[i] * a_exp;

        g_z_0_0_0_x_xz_xx_yz[i] = 2.0 * g_xz_xz_xx_yz[i] * a_exp;

        g_z_0_0_0_x_xz_xx_zz[i] = 2.0 * g_xz_xz_xx_zz[i] * a_exp;
    }
    // integrals block (1374-1380)

    #pragma omp simd aligned(g_xz_xz_xy_xx, g_xz_xz_xy_xy, g_xz_xz_xy_xz, g_xz_xz_xy_yy, g_xz_xz_xy_yz, g_xz_xz_xy_zz, g_z_0_0_0_x_xz_xy_xx, g_z_0_0_0_x_xz_xy_xy, g_z_0_0_0_x_xz_xy_xz, g_z_0_0_0_x_xz_xy_yy, g_z_0_0_0_x_xz_xy_yz, g_z_0_0_0_x_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xz_xy_xx[i] = 2.0 * g_xz_xz_xy_xx[i] * a_exp;

        g_z_0_0_0_x_xz_xy_xy[i] = 2.0 * g_xz_xz_xy_xy[i] * a_exp;

        g_z_0_0_0_x_xz_xy_xz[i] = 2.0 * g_xz_xz_xy_xz[i] * a_exp;

        g_z_0_0_0_x_xz_xy_yy[i] = 2.0 * g_xz_xz_xy_yy[i] * a_exp;

        g_z_0_0_0_x_xz_xy_yz[i] = 2.0 * g_xz_xz_xy_yz[i] * a_exp;

        g_z_0_0_0_x_xz_xy_zz[i] = 2.0 * g_xz_xz_xy_zz[i] * a_exp;
    }
    // integrals block (1380-1386)

    #pragma omp simd aligned(g_xz_xz_xz_xx, g_xz_xz_xz_xy, g_xz_xz_xz_xz, g_xz_xz_xz_yy, g_xz_xz_xz_yz, g_xz_xz_xz_zz, g_z_0_0_0_x_xz_xz_xx, g_z_0_0_0_x_xz_xz_xy, g_z_0_0_0_x_xz_xz_xz, g_z_0_0_0_x_xz_xz_yy, g_z_0_0_0_x_xz_xz_yz, g_z_0_0_0_x_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xz_xz_xx[i] = 2.0 * g_xz_xz_xz_xx[i] * a_exp;

        g_z_0_0_0_x_xz_xz_xy[i] = 2.0 * g_xz_xz_xz_xy[i] * a_exp;

        g_z_0_0_0_x_xz_xz_xz[i] = 2.0 * g_xz_xz_xz_xz[i] * a_exp;

        g_z_0_0_0_x_xz_xz_yy[i] = 2.0 * g_xz_xz_xz_yy[i] * a_exp;

        g_z_0_0_0_x_xz_xz_yz[i] = 2.0 * g_xz_xz_xz_yz[i] * a_exp;

        g_z_0_0_0_x_xz_xz_zz[i] = 2.0 * g_xz_xz_xz_zz[i] * a_exp;
    }
    // integrals block (1386-1392)

    #pragma omp simd aligned(g_xz_xz_yy_xx, g_xz_xz_yy_xy, g_xz_xz_yy_xz, g_xz_xz_yy_yy, g_xz_xz_yy_yz, g_xz_xz_yy_zz, g_z_0_0_0_x_xz_yy_xx, g_z_0_0_0_x_xz_yy_xy, g_z_0_0_0_x_xz_yy_xz, g_z_0_0_0_x_xz_yy_yy, g_z_0_0_0_x_xz_yy_yz, g_z_0_0_0_x_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xz_yy_xx[i] = 2.0 * g_xz_xz_yy_xx[i] * a_exp;

        g_z_0_0_0_x_xz_yy_xy[i] = 2.0 * g_xz_xz_yy_xy[i] * a_exp;

        g_z_0_0_0_x_xz_yy_xz[i] = 2.0 * g_xz_xz_yy_xz[i] * a_exp;

        g_z_0_0_0_x_xz_yy_yy[i] = 2.0 * g_xz_xz_yy_yy[i] * a_exp;

        g_z_0_0_0_x_xz_yy_yz[i] = 2.0 * g_xz_xz_yy_yz[i] * a_exp;

        g_z_0_0_0_x_xz_yy_zz[i] = 2.0 * g_xz_xz_yy_zz[i] * a_exp;
    }
    // integrals block (1392-1398)

    #pragma omp simd aligned(g_xz_xz_yz_xx, g_xz_xz_yz_xy, g_xz_xz_yz_xz, g_xz_xz_yz_yy, g_xz_xz_yz_yz, g_xz_xz_yz_zz, g_z_0_0_0_x_xz_yz_xx, g_z_0_0_0_x_xz_yz_xy, g_z_0_0_0_x_xz_yz_xz, g_z_0_0_0_x_xz_yz_yy, g_z_0_0_0_x_xz_yz_yz, g_z_0_0_0_x_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xz_yz_xx[i] = 2.0 * g_xz_xz_yz_xx[i] * a_exp;

        g_z_0_0_0_x_xz_yz_xy[i] = 2.0 * g_xz_xz_yz_xy[i] * a_exp;

        g_z_0_0_0_x_xz_yz_xz[i] = 2.0 * g_xz_xz_yz_xz[i] * a_exp;

        g_z_0_0_0_x_xz_yz_yy[i] = 2.0 * g_xz_xz_yz_yy[i] * a_exp;

        g_z_0_0_0_x_xz_yz_yz[i] = 2.0 * g_xz_xz_yz_yz[i] * a_exp;

        g_z_0_0_0_x_xz_yz_zz[i] = 2.0 * g_xz_xz_yz_zz[i] * a_exp;
    }
    // integrals block (1398-1404)

    #pragma omp simd aligned(g_xz_xz_zz_xx, g_xz_xz_zz_xy, g_xz_xz_zz_xz, g_xz_xz_zz_yy, g_xz_xz_zz_yz, g_xz_xz_zz_zz, g_z_0_0_0_x_xz_zz_xx, g_z_0_0_0_x_xz_zz_xy, g_z_0_0_0_x_xz_zz_xz, g_z_0_0_0_x_xz_zz_yy, g_z_0_0_0_x_xz_zz_yz, g_z_0_0_0_x_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xz_zz_xx[i] = 2.0 * g_xz_xz_zz_xx[i] * a_exp;

        g_z_0_0_0_x_xz_zz_xy[i] = 2.0 * g_xz_xz_zz_xy[i] * a_exp;

        g_z_0_0_0_x_xz_zz_xz[i] = 2.0 * g_xz_xz_zz_xz[i] * a_exp;

        g_z_0_0_0_x_xz_zz_yy[i] = 2.0 * g_xz_xz_zz_yy[i] * a_exp;

        g_z_0_0_0_x_xz_zz_yz[i] = 2.0 * g_xz_xz_zz_yz[i] * a_exp;

        g_z_0_0_0_x_xz_zz_zz[i] = 2.0 * g_xz_xz_zz_zz[i] * a_exp;
    }
    // integrals block (1404-1410)

    #pragma omp simd aligned(g_xz_yy_xx_xx, g_xz_yy_xx_xy, g_xz_yy_xx_xz, g_xz_yy_xx_yy, g_xz_yy_xx_yz, g_xz_yy_xx_zz, g_z_0_0_0_x_yy_xx_xx, g_z_0_0_0_x_yy_xx_xy, g_z_0_0_0_x_yy_xx_xz, g_z_0_0_0_x_yy_xx_yy, g_z_0_0_0_x_yy_xx_yz, g_z_0_0_0_x_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yy_xx_xx[i] = 2.0 * g_xz_yy_xx_xx[i] * a_exp;

        g_z_0_0_0_x_yy_xx_xy[i] = 2.0 * g_xz_yy_xx_xy[i] * a_exp;

        g_z_0_0_0_x_yy_xx_xz[i] = 2.0 * g_xz_yy_xx_xz[i] * a_exp;

        g_z_0_0_0_x_yy_xx_yy[i] = 2.0 * g_xz_yy_xx_yy[i] * a_exp;

        g_z_0_0_0_x_yy_xx_yz[i] = 2.0 * g_xz_yy_xx_yz[i] * a_exp;

        g_z_0_0_0_x_yy_xx_zz[i] = 2.0 * g_xz_yy_xx_zz[i] * a_exp;
    }
    // integrals block (1410-1416)

    #pragma omp simd aligned(g_xz_yy_xy_xx, g_xz_yy_xy_xy, g_xz_yy_xy_xz, g_xz_yy_xy_yy, g_xz_yy_xy_yz, g_xz_yy_xy_zz, g_z_0_0_0_x_yy_xy_xx, g_z_0_0_0_x_yy_xy_xy, g_z_0_0_0_x_yy_xy_xz, g_z_0_0_0_x_yy_xy_yy, g_z_0_0_0_x_yy_xy_yz, g_z_0_0_0_x_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yy_xy_xx[i] = 2.0 * g_xz_yy_xy_xx[i] * a_exp;

        g_z_0_0_0_x_yy_xy_xy[i] = 2.0 * g_xz_yy_xy_xy[i] * a_exp;

        g_z_0_0_0_x_yy_xy_xz[i] = 2.0 * g_xz_yy_xy_xz[i] * a_exp;

        g_z_0_0_0_x_yy_xy_yy[i] = 2.0 * g_xz_yy_xy_yy[i] * a_exp;

        g_z_0_0_0_x_yy_xy_yz[i] = 2.0 * g_xz_yy_xy_yz[i] * a_exp;

        g_z_0_0_0_x_yy_xy_zz[i] = 2.0 * g_xz_yy_xy_zz[i] * a_exp;
    }
    // integrals block (1416-1422)

    #pragma omp simd aligned(g_xz_yy_xz_xx, g_xz_yy_xz_xy, g_xz_yy_xz_xz, g_xz_yy_xz_yy, g_xz_yy_xz_yz, g_xz_yy_xz_zz, g_z_0_0_0_x_yy_xz_xx, g_z_0_0_0_x_yy_xz_xy, g_z_0_0_0_x_yy_xz_xz, g_z_0_0_0_x_yy_xz_yy, g_z_0_0_0_x_yy_xz_yz, g_z_0_0_0_x_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yy_xz_xx[i] = 2.0 * g_xz_yy_xz_xx[i] * a_exp;

        g_z_0_0_0_x_yy_xz_xy[i] = 2.0 * g_xz_yy_xz_xy[i] * a_exp;

        g_z_0_0_0_x_yy_xz_xz[i] = 2.0 * g_xz_yy_xz_xz[i] * a_exp;

        g_z_0_0_0_x_yy_xz_yy[i] = 2.0 * g_xz_yy_xz_yy[i] * a_exp;

        g_z_0_0_0_x_yy_xz_yz[i] = 2.0 * g_xz_yy_xz_yz[i] * a_exp;

        g_z_0_0_0_x_yy_xz_zz[i] = 2.0 * g_xz_yy_xz_zz[i] * a_exp;
    }
    // integrals block (1422-1428)

    #pragma omp simd aligned(g_xz_yy_yy_xx, g_xz_yy_yy_xy, g_xz_yy_yy_xz, g_xz_yy_yy_yy, g_xz_yy_yy_yz, g_xz_yy_yy_zz, g_z_0_0_0_x_yy_yy_xx, g_z_0_0_0_x_yy_yy_xy, g_z_0_0_0_x_yy_yy_xz, g_z_0_0_0_x_yy_yy_yy, g_z_0_0_0_x_yy_yy_yz, g_z_0_0_0_x_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yy_yy_xx[i] = 2.0 * g_xz_yy_yy_xx[i] * a_exp;

        g_z_0_0_0_x_yy_yy_xy[i] = 2.0 * g_xz_yy_yy_xy[i] * a_exp;

        g_z_0_0_0_x_yy_yy_xz[i] = 2.0 * g_xz_yy_yy_xz[i] * a_exp;

        g_z_0_0_0_x_yy_yy_yy[i] = 2.0 * g_xz_yy_yy_yy[i] * a_exp;

        g_z_0_0_0_x_yy_yy_yz[i] = 2.0 * g_xz_yy_yy_yz[i] * a_exp;

        g_z_0_0_0_x_yy_yy_zz[i] = 2.0 * g_xz_yy_yy_zz[i] * a_exp;
    }
    // integrals block (1428-1434)

    #pragma omp simd aligned(g_xz_yy_yz_xx, g_xz_yy_yz_xy, g_xz_yy_yz_xz, g_xz_yy_yz_yy, g_xz_yy_yz_yz, g_xz_yy_yz_zz, g_z_0_0_0_x_yy_yz_xx, g_z_0_0_0_x_yy_yz_xy, g_z_0_0_0_x_yy_yz_xz, g_z_0_0_0_x_yy_yz_yy, g_z_0_0_0_x_yy_yz_yz, g_z_0_0_0_x_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yy_yz_xx[i] = 2.0 * g_xz_yy_yz_xx[i] * a_exp;

        g_z_0_0_0_x_yy_yz_xy[i] = 2.0 * g_xz_yy_yz_xy[i] * a_exp;

        g_z_0_0_0_x_yy_yz_xz[i] = 2.0 * g_xz_yy_yz_xz[i] * a_exp;

        g_z_0_0_0_x_yy_yz_yy[i] = 2.0 * g_xz_yy_yz_yy[i] * a_exp;

        g_z_0_0_0_x_yy_yz_yz[i] = 2.0 * g_xz_yy_yz_yz[i] * a_exp;

        g_z_0_0_0_x_yy_yz_zz[i] = 2.0 * g_xz_yy_yz_zz[i] * a_exp;
    }
    // integrals block (1434-1440)

    #pragma omp simd aligned(g_xz_yy_zz_xx, g_xz_yy_zz_xy, g_xz_yy_zz_xz, g_xz_yy_zz_yy, g_xz_yy_zz_yz, g_xz_yy_zz_zz, g_z_0_0_0_x_yy_zz_xx, g_z_0_0_0_x_yy_zz_xy, g_z_0_0_0_x_yy_zz_xz, g_z_0_0_0_x_yy_zz_yy, g_z_0_0_0_x_yy_zz_yz, g_z_0_0_0_x_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yy_zz_xx[i] = 2.0 * g_xz_yy_zz_xx[i] * a_exp;

        g_z_0_0_0_x_yy_zz_xy[i] = 2.0 * g_xz_yy_zz_xy[i] * a_exp;

        g_z_0_0_0_x_yy_zz_xz[i] = 2.0 * g_xz_yy_zz_xz[i] * a_exp;

        g_z_0_0_0_x_yy_zz_yy[i] = 2.0 * g_xz_yy_zz_yy[i] * a_exp;

        g_z_0_0_0_x_yy_zz_yz[i] = 2.0 * g_xz_yy_zz_yz[i] * a_exp;

        g_z_0_0_0_x_yy_zz_zz[i] = 2.0 * g_xz_yy_zz_zz[i] * a_exp;
    }
    // integrals block (1440-1446)

    #pragma omp simd aligned(g_xz_yz_xx_xx, g_xz_yz_xx_xy, g_xz_yz_xx_xz, g_xz_yz_xx_yy, g_xz_yz_xx_yz, g_xz_yz_xx_zz, g_z_0_0_0_x_yz_xx_xx, g_z_0_0_0_x_yz_xx_xy, g_z_0_0_0_x_yz_xx_xz, g_z_0_0_0_x_yz_xx_yy, g_z_0_0_0_x_yz_xx_yz, g_z_0_0_0_x_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yz_xx_xx[i] = 2.0 * g_xz_yz_xx_xx[i] * a_exp;

        g_z_0_0_0_x_yz_xx_xy[i] = 2.0 * g_xz_yz_xx_xy[i] * a_exp;

        g_z_0_0_0_x_yz_xx_xz[i] = 2.0 * g_xz_yz_xx_xz[i] * a_exp;

        g_z_0_0_0_x_yz_xx_yy[i] = 2.0 * g_xz_yz_xx_yy[i] * a_exp;

        g_z_0_0_0_x_yz_xx_yz[i] = 2.0 * g_xz_yz_xx_yz[i] * a_exp;

        g_z_0_0_0_x_yz_xx_zz[i] = 2.0 * g_xz_yz_xx_zz[i] * a_exp;
    }
    // integrals block (1446-1452)

    #pragma omp simd aligned(g_xz_yz_xy_xx, g_xz_yz_xy_xy, g_xz_yz_xy_xz, g_xz_yz_xy_yy, g_xz_yz_xy_yz, g_xz_yz_xy_zz, g_z_0_0_0_x_yz_xy_xx, g_z_0_0_0_x_yz_xy_xy, g_z_0_0_0_x_yz_xy_xz, g_z_0_0_0_x_yz_xy_yy, g_z_0_0_0_x_yz_xy_yz, g_z_0_0_0_x_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yz_xy_xx[i] = 2.0 * g_xz_yz_xy_xx[i] * a_exp;

        g_z_0_0_0_x_yz_xy_xy[i] = 2.0 * g_xz_yz_xy_xy[i] * a_exp;

        g_z_0_0_0_x_yz_xy_xz[i] = 2.0 * g_xz_yz_xy_xz[i] * a_exp;

        g_z_0_0_0_x_yz_xy_yy[i] = 2.0 * g_xz_yz_xy_yy[i] * a_exp;

        g_z_0_0_0_x_yz_xy_yz[i] = 2.0 * g_xz_yz_xy_yz[i] * a_exp;

        g_z_0_0_0_x_yz_xy_zz[i] = 2.0 * g_xz_yz_xy_zz[i] * a_exp;
    }
    // integrals block (1452-1458)

    #pragma omp simd aligned(g_xz_yz_xz_xx, g_xz_yz_xz_xy, g_xz_yz_xz_xz, g_xz_yz_xz_yy, g_xz_yz_xz_yz, g_xz_yz_xz_zz, g_z_0_0_0_x_yz_xz_xx, g_z_0_0_0_x_yz_xz_xy, g_z_0_0_0_x_yz_xz_xz, g_z_0_0_0_x_yz_xz_yy, g_z_0_0_0_x_yz_xz_yz, g_z_0_0_0_x_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yz_xz_xx[i] = 2.0 * g_xz_yz_xz_xx[i] * a_exp;

        g_z_0_0_0_x_yz_xz_xy[i] = 2.0 * g_xz_yz_xz_xy[i] * a_exp;

        g_z_0_0_0_x_yz_xz_xz[i] = 2.0 * g_xz_yz_xz_xz[i] * a_exp;

        g_z_0_0_0_x_yz_xz_yy[i] = 2.0 * g_xz_yz_xz_yy[i] * a_exp;

        g_z_0_0_0_x_yz_xz_yz[i] = 2.0 * g_xz_yz_xz_yz[i] * a_exp;

        g_z_0_0_0_x_yz_xz_zz[i] = 2.0 * g_xz_yz_xz_zz[i] * a_exp;
    }
    // integrals block (1458-1464)

    #pragma omp simd aligned(g_xz_yz_yy_xx, g_xz_yz_yy_xy, g_xz_yz_yy_xz, g_xz_yz_yy_yy, g_xz_yz_yy_yz, g_xz_yz_yy_zz, g_z_0_0_0_x_yz_yy_xx, g_z_0_0_0_x_yz_yy_xy, g_z_0_0_0_x_yz_yy_xz, g_z_0_0_0_x_yz_yy_yy, g_z_0_0_0_x_yz_yy_yz, g_z_0_0_0_x_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yz_yy_xx[i] = 2.0 * g_xz_yz_yy_xx[i] * a_exp;

        g_z_0_0_0_x_yz_yy_xy[i] = 2.0 * g_xz_yz_yy_xy[i] * a_exp;

        g_z_0_0_0_x_yz_yy_xz[i] = 2.0 * g_xz_yz_yy_xz[i] * a_exp;

        g_z_0_0_0_x_yz_yy_yy[i] = 2.0 * g_xz_yz_yy_yy[i] * a_exp;

        g_z_0_0_0_x_yz_yy_yz[i] = 2.0 * g_xz_yz_yy_yz[i] * a_exp;

        g_z_0_0_0_x_yz_yy_zz[i] = 2.0 * g_xz_yz_yy_zz[i] * a_exp;
    }
    // integrals block (1464-1470)

    #pragma omp simd aligned(g_xz_yz_yz_xx, g_xz_yz_yz_xy, g_xz_yz_yz_xz, g_xz_yz_yz_yy, g_xz_yz_yz_yz, g_xz_yz_yz_zz, g_z_0_0_0_x_yz_yz_xx, g_z_0_0_0_x_yz_yz_xy, g_z_0_0_0_x_yz_yz_xz, g_z_0_0_0_x_yz_yz_yy, g_z_0_0_0_x_yz_yz_yz, g_z_0_0_0_x_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yz_yz_xx[i] = 2.0 * g_xz_yz_yz_xx[i] * a_exp;

        g_z_0_0_0_x_yz_yz_xy[i] = 2.0 * g_xz_yz_yz_xy[i] * a_exp;

        g_z_0_0_0_x_yz_yz_xz[i] = 2.0 * g_xz_yz_yz_xz[i] * a_exp;

        g_z_0_0_0_x_yz_yz_yy[i] = 2.0 * g_xz_yz_yz_yy[i] * a_exp;

        g_z_0_0_0_x_yz_yz_yz[i] = 2.0 * g_xz_yz_yz_yz[i] * a_exp;

        g_z_0_0_0_x_yz_yz_zz[i] = 2.0 * g_xz_yz_yz_zz[i] * a_exp;
    }
    // integrals block (1470-1476)

    #pragma omp simd aligned(g_xz_yz_zz_xx, g_xz_yz_zz_xy, g_xz_yz_zz_xz, g_xz_yz_zz_yy, g_xz_yz_zz_yz, g_xz_yz_zz_zz, g_z_0_0_0_x_yz_zz_xx, g_z_0_0_0_x_yz_zz_xy, g_z_0_0_0_x_yz_zz_xz, g_z_0_0_0_x_yz_zz_yy, g_z_0_0_0_x_yz_zz_yz, g_z_0_0_0_x_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yz_zz_xx[i] = 2.0 * g_xz_yz_zz_xx[i] * a_exp;

        g_z_0_0_0_x_yz_zz_xy[i] = 2.0 * g_xz_yz_zz_xy[i] * a_exp;

        g_z_0_0_0_x_yz_zz_xz[i] = 2.0 * g_xz_yz_zz_xz[i] * a_exp;

        g_z_0_0_0_x_yz_zz_yy[i] = 2.0 * g_xz_yz_zz_yy[i] * a_exp;

        g_z_0_0_0_x_yz_zz_yz[i] = 2.0 * g_xz_yz_zz_yz[i] * a_exp;

        g_z_0_0_0_x_yz_zz_zz[i] = 2.0 * g_xz_yz_zz_zz[i] * a_exp;
    }
    // integrals block (1476-1482)

    #pragma omp simd aligned(g_xz_zz_xx_xx, g_xz_zz_xx_xy, g_xz_zz_xx_xz, g_xz_zz_xx_yy, g_xz_zz_xx_yz, g_xz_zz_xx_zz, g_z_0_0_0_x_zz_xx_xx, g_z_0_0_0_x_zz_xx_xy, g_z_0_0_0_x_zz_xx_xz, g_z_0_0_0_x_zz_xx_yy, g_z_0_0_0_x_zz_xx_yz, g_z_0_0_0_x_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_zz_xx_xx[i] = 2.0 * g_xz_zz_xx_xx[i] * a_exp;

        g_z_0_0_0_x_zz_xx_xy[i] = 2.0 * g_xz_zz_xx_xy[i] * a_exp;

        g_z_0_0_0_x_zz_xx_xz[i] = 2.0 * g_xz_zz_xx_xz[i] * a_exp;

        g_z_0_0_0_x_zz_xx_yy[i] = 2.0 * g_xz_zz_xx_yy[i] * a_exp;

        g_z_0_0_0_x_zz_xx_yz[i] = 2.0 * g_xz_zz_xx_yz[i] * a_exp;

        g_z_0_0_0_x_zz_xx_zz[i] = 2.0 * g_xz_zz_xx_zz[i] * a_exp;
    }
    // integrals block (1482-1488)

    #pragma omp simd aligned(g_xz_zz_xy_xx, g_xz_zz_xy_xy, g_xz_zz_xy_xz, g_xz_zz_xy_yy, g_xz_zz_xy_yz, g_xz_zz_xy_zz, g_z_0_0_0_x_zz_xy_xx, g_z_0_0_0_x_zz_xy_xy, g_z_0_0_0_x_zz_xy_xz, g_z_0_0_0_x_zz_xy_yy, g_z_0_0_0_x_zz_xy_yz, g_z_0_0_0_x_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_zz_xy_xx[i] = 2.0 * g_xz_zz_xy_xx[i] * a_exp;

        g_z_0_0_0_x_zz_xy_xy[i] = 2.0 * g_xz_zz_xy_xy[i] * a_exp;

        g_z_0_0_0_x_zz_xy_xz[i] = 2.0 * g_xz_zz_xy_xz[i] * a_exp;

        g_z_0_0_0_x_zz_xy_yy[i] = 2.0 * g_xz_zz_xy_yy[i] * a_exp;

        g_z_0_0_0_x_zz_xy_yz[i] = 2.0 * g_xz_zz_xy_yz[i] * a_exp;

        g_z_0_0_0_x_zz_xy_zz[i] = 2.0 * g_xz_zz_xy_zz[i] * a_exp;
    }
    // integrals block (1488-1494)

    #pragma omp simd aligned(g_xz_zz_xz_xx, g_xz_zz_xz_xy, g_xz_zz_xz_xz, g_xz_zz_xz_yy, g_xz_zz_xz_yz, g_xz_zz_xz_zz, g_z_0_0_0_x_zz_xz_xx, g_z_0_0_0_x_zz_xz_xy, g_z_0_0_0_x_zz_xz_xz, g_z_0_0_0_x_zz_xz_yy, g_z_0_0_0_x_zz_xz_yz, g_z_0_0_0_x_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_zz_xz_xx[i] = 2.0 * g_xz_zz_xz_xx[i] * a_exp;

        g_z_0_0_0_x_zz_xz_xy[i] = 2.0 * g_xz_zz_xz_xy[i] * a_exp;

        g_z_0_0_0_x_zz_xz_xz[i] = 2.0 * g_xz_zz_xz_xz[i] * a_exp;

        g_z_0_0_0_x_zz_xz_yy[i] = 2.0 * g_xz_zz_xz_yy[i] * a_exp;

        g_z_0_0_0_x_zz_xz_yz[i] = 2.0 * g_xz_zz_xz_yz[i] * a_exp;

        g_z_0_0_0_x_zz_xz_zz[i] = 2.0 * g_xz_zz_xz_zz[i] * a_exp;
    }
    // integrals block (1494-1500)

    #pragma omp simd aligned(g_xz_zz_yy_xx, g_xz_zz_yy_xy, g_xz_zz_yy_xz, g_xz_zz_yy_yy, g_xz_zz_yy_yz, g_xz_zz_yy_zz, g_z_0_0_0_x_zz_yy_xx, g_z_0_0_0_x_zz_yy_xy, g_z_0_0_0_x_zz_yy_xz, g_z_0_0_0_x_zz_yy_yy, g_z_0_0_0_x_zz_yy_yz, g_z_0_0_0_x_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_zz_yy_xx[i] = 2.0 * g_xz_zz_yy_xx[i] * a_exp;

        g_z_0_0_0_x_zz_yy_xy[i] = 2.0 * g_xz_zz_yy_xy[i] * a_exp;

        g_z_0_0_0_x_zz_yy_xz[i] = 2.0 * g_xz_zz_yy_xz[i] * a_exp;

        g_z_0_0_0_x_zz_yy_yy[i] = 2.0 * g_xz_zz_yy_yy[i] * a_exp;

        g_z_0_0_0_x_zz_yy_yz[i] = 2.0 * g_xz_zz_yy_yz[i] * a_exp;

        g_z_0_0_0_x_zz_yy_zz[i] = 2.0 * g_xz_zz_yy_zz[i] * a_exp;
    }
    // integrals block (1500-1506)

    #pragma omp simd aligned(g_xz_zz_yz_xx, g_xz_zz_yz_xy, g_xz_zz_yz_xz, g_xz_zz_yz_yy, g_xz_zz_yz_yz, g_xz_zz_yz_zz, g_z_0_0_0_x_zz_yz_xx, g_z_0_0_0_x_zz_yz_xy, g_z_0_0_0_x_zz_yz_xz, g_z_0_0_0_x_zz_yz_yy, g_z_0_0_0_x_zz_yz_yz, g_z_0_0_0_x_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_zz_yz_xx[i] = 2.0 * g_xz_zz_yz_xx[i] * a_exp;

        g_z_0_0_0_x_zz_yz_xy[i] = 2.0 * g_xz_zz_yz_xy[i] * a_exp;

        g_z_0_0_0_x_zz_yz_xz[i] = 2.0 * g_xz_zz_yz_xz[i] * a_exp;

        g_z_0_0_0_x_zz_yz_yy[i] = 2.0 * g_xz_zz_yz_yy[i] * a_exp;

        g_z_0_0_0_x_zz_yz_yz[i] = 2.0 * g_xz_zz_yz_yz[i] * a_exp;

        g_z_0_0_0_x_zz_yz_zz[i] = 2.0 * g_xz_zz_yz_zz[i] * a_exp;
    }
    // integrals block (1506-1512)

    #pragma omp simd aligned(g_xz_zz_zz_xx, g_xz_zz_zz_xy, g_xz_zz_zz_xz, g_xz_zz_zz_yy, g_xz_zz_zz_yz, g_xz_zz_zz_zz, g_z_0_0_0_x_zz_zz_xx, g_z_0_0_0_x_zz_zz_xy, g_z_0_0_0_x_zz_zz_xz, g_z_0_0_0_x_zz_zz_yy, g_z_0_0_0_x_zz_zz_yz, g_z_0_0_0_x_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_zz_zz_xx[i] = 2.0 * g_xz_zz_zz_xx[i] * a_exp;

        g_z_0_0_0_x_zz_zz_xy[i] = 2.0 * g_xz_zz_zz_xy[i] * a_exp;

        g_z_0_0_0_x_zz_zz_xz[i] = 2.0 * g_xz_zz_zz_xz[i] * a_exp;

        g_z_0_0_0_x_zz_zz_yy[i] = 2.0 * g_xz_zz_zz_yy[i] * a_exp;

        g_z_0_0_0_x_zz_zz_yz[i] = 2.0 * g_xz_zz_zz_yz[i] * a_exp;

        g_z_0_0_0_x_zz_zz_zz[i] = 2.0 * g_xz_zz_zz_zz[i] * a_exp;
    }
    // integrals block (1512-1518)

    #pragma omp simd aligned(g_yz_xx_xx_xx, g_yz_xx_xx_xy, g_yz_xx_xx_xz, g_yz_xx_xx_yy, g_yz_xx_xx_yz, g_yz_xx_xx_zz, g_z_0_0_0_y_xx_xx_xx, g_z_0_0_0_y_xx_xx_xy, g_z_0_0_0_y_xx_xx_xz, g_z_0_0_0_y_xx_xx_yy, g_z_0_0_0_y_xx_xx_yz, g_z_0_0_0_y_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xx_xx_xx[i] = 2.0 * g_yz_xx_xx_xx[i] * a_exp;

        g_z_0_0_0_y_xx_xx_xy[i] = 2.0 * g_yz_xx_xx_xy[i] * a_exp;

        g_z_0_0_0_y_xx_xx_xz[i] = 2.0 * g_yz_xx_xx_xz[i] * a_exp;

        g_z_0_0_0_y_xx_xx_yy[i] = 2.0 * g_yz_xx_xx_yy[i] * a_exp;

        g_z_0_0_0_y_xx_xx_yz[i] = 2.0 * g_yz_xx_xx_yz[i] * a_exp;

        g_z_0_0_0_y_xx_xx_zz[i] = 2.0 * g_yz_xx_xx_zz[i] * a_exp;
    }
    // integrals block (1518-1524)

    #pragma omp simd aligned(g_yz_xx_xy_xx, g_yz_xx_xy_xy, g_yz_xx_xy_xz, g_yz_xx_xy_yy, g_yz_xx_xy_yz, g_yz_xx_xy_zz, g_z_0_0_0_y_xx_xy_xx, g_z_0_0_0_y_xx_xy_xy, g_z_0_0_0_y_xx_xy_xz, g_z_0_0_0_y_xx_xy_yy, g_z_0_0_0_y_xx_xy_yz, g_z_0_0_0_y_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xx_xy_xx[i] = 2.0 * g_yz_xx_xy_xx[i] * a_exp;

        g_z_0_0_0_y_xx_xy_xy[i] = 2.0 * g_yz_xx_xy_xy[i] * a_exp;

        g_z_0_0_0_y_xx_xy_xz[i] = 2.0 * g_yz_xx_xy_xz[i] * a_exp;

        g_z_0_0_0_y_xx_xy_yy[i] = 2.0 * g_yz_xx_xy_yy[i] * a_exp;

        g_z_0_0_0_y_xx_xy_yz[i] = 2.0 * g_yz_xx_xy_yz[i] * a_exp;

        g_z_0_0_0_y_xx_xy_zz[i] = 2.0 * g_yz_xx_xy_zz[i] * a_exp;
    }
    // integrals block (1524-1530)

    #pragma omp simd aligned(g_yz_xx_xz_xx, g_yz_xx_xz_xy, g_yz_xx_xz_xz, g_yz_xx_xz_yy, g_yz_xx_xz_yz, g_yz_xx_xz_zz, g_z_0_0_0_y_xx_xz_xx, g_z_0_0_0_y_xx_xz_xy, g_z_0_0_0_y_xx_xz_xz, g_z_0_0_0_y_xx_xz_yy, g_z_0_0_0_y_xx_xz_yz, g_z_0_0_0_y_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xx_xz_xx[i] = 2.0 * g_yz_xx_xz_xx[i] * a_exp;

        g_z_0_0_0_y_xx_xz_xy[i] = 2.0 * g_yz_xx_xz_xy[i] * a_exp;

        g_z_0_0_0_y_xx_xz_xz[i] = 2.0 * g_yz_xx_xz_xz[i] * a_exp;

        g_z_0_0_0_y_xx_xz_yy[i] = 2.0 * g_yz_xx_xz_yy[i] * a_exp;

        g_z_0_0_0_y_xx_xz_yz[i] = 2.0 * g_yz_xx_xz_yz[i] * a_exp;

        g_z_0_0_0_y_xx_xz_zz[i] = 2.0 * g_yz_xx_xz_zz[i] * a_exp;
    }
    // integrals block (1530-1536)

    #pragma omp simd aligned(g_yz_xx_yy_xx, g_yz_xx_yy_xy, g_yz_xx_yy_xz, g_yz_xx_yy_yy, g_yz_xx_yy_yz, g_yz_xx_yy_zz, g_z_0_0_0_y_xx_yy_xx, g_z_0_0_0_y_xx_yy_xy, g_z_0_0_0_y_xx_yy_xz, g_z_0_0_0_y_xx_yy_yy, g_z_0_0_0_y_xx_yy_yz, g_z_0_0_0_y_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xx_yy_xx[i] = 2.0 * g_yz_xx_yy_xx[i] * a_exp;

        g_z_0_0_0_y_xx_yy_xy[i] = 2.0 * g_yz_xx_yy_xy[i] * a_exp;

        g_z_0_0_0_y_xx_yy_xz[i] = 2.0 * g_yz_xx_yy_xz[i] * a_exp;

        g_z_0_0_0_y_xx_yy_yy[i] = 2.0 * g_yz_xx_yy_yy[i] * a_exp;

        g_z_0_0_0_y_xx_yy_yz[i] = 2.0 * g_yz_xx_yy_yz[i] * a_exp;

        g_z_0_0_0_y_xx_yy_zz[i] = 2.0 * g_yz_xx_yy_zz[i] * a_exp;
    }
    // integrals block (1536-1542)

    #pragma omp simd aligned(g_yz_xx_yz_xx, g_yz_xx_yz_xy, g_yz_xx_yz_xz, g_yz_xx_yz_yy, g_yz_xx_yz_yz, g_yz_xx_yz_zz, g_z_0_0_0_y_xx_yz_xx, g_z_0_0_0_y_xx_yz_xy, g_z_0_0_0_y_xx_yz_xz, g_z_0_0_0_y_xx_yz_yy, g_z_0_0_0_y_xx_yz_yz, g_z_0_0_0_y_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xx_yz_xx[i] = 2.0 * g_yz_xx_yz_xx[i] * a_exp;

        g_z_0_0_0_y_xx_yz_xy[i] = 2.0 * g_yz_xx_yz_xy[i] * a_exp;

        g_z_0_0_0_y_xx_yz_xz[i] = 2.0 * g_yz_xx_yz_xz[i] * a_exp;

        g_z_0_0_0_y_xx_yz_yy[i] = 2.0 * g_yz_xx_yz_yy[i] * a_exp;

        g_z_0_0_0_y_xx_yz_yz[i] = 2.0 * g_yz_xx_yz_yz[i] * a_exp;

        g_z_0_0_0_y_xx_yz_zz[i] = 2.0 * g_yz_xx_yz_zz[i] * a_exp;
    }
    // integrals block (1542-1548)

    #pragma omp simd aligned(g_yz_xx_zz_xx, g_yz_xx_zz_xy, g_yz_xx_zz_xz, g_yz_xx_zz_yy, g_yz_xx_zz_yz, g_yz_xx_zz_zz, g_z_0_0_0_y_xx_zz_xx, g_z_0_0_0_y_xx_zz_xy, g_z_0_0_0_y_xx_zz_xz, g_z_0_0_0_y_xx_zz_yy, g_z_0_0_0_y_xx_zz_yz, g_z_0_0_0_y_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xx_zz_xx[i] = 2.0 * g_yz_xx_zz_xx[i] * a_exp;

        g_z_0_0_0_y_xx_zz_xy[i] = 2.0 * g_yz_xx_zz_xy[i] * a_exp;

        g_z_0_0_0_y_xx_zz_xz[i] = 2.0 * g_yz_xx_zz_xz[i] * a_exp;

        g_z_0_0_0_y_xx_zz_yy[i] = 2.0 * g_yz_xx_zz_yy[i] * a_exp;

        g_z_0_0_0_y_xx_zz_yz[i] = 2.0 * g_yz_xx_zz_yz[i] * a_exp;

        g_z_0_0_0_y_xx_zz_zz[i] = 2.0 * g_yz_xx_zz_zz[i] * a_exp;
    }
    // integrals block (1548-1554)

    #pragma omp simd aligned(g_yz_xy_xx_xx, g_yz_xy_xx_xy, g_yz_xy_xx_xz, g_yz_xy_xx_yy, g_yz_xy_xx_yz, g_yz_xy_xx_zz, g_z_0_0_0_y_xy_xx_xx, g_z_0_0_0_y_xy_xx_xy, g_z_0_0_0_y_xy_xx_xz, g_z_0_0_0_y_xy_xx_yy, g_z_0_0_0_y_xy_xx_yz, g_z_0_0_0_y_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xy_xx_xx[i] = 2.0 * g_yz_xy_xx_xx[i] * a_exp;

        g_z_0_0_0_y_xy_xx_xy[i] = 2.0 * g_yz_xy_xx_xy[i] * a_exp;

        g_z_0_0_0_y_xy_xx_xz[i] = 2.0 * g_yz_xy_xx_xz[i] * a_exp;

        g_z_0_0_0_y_xy_xx_yy[i] = 2.0 * g_yz_xy_xx_yy[i] * a_exp;

        g_z_0_0_0_y_xy_xx_yz[i] = 2.0 * g_yz_xy_xx_yz[i] * a_exp;

        g_z_0_0_0_y_xy_xx_zz[i] = 2.0 * g_yz_xy_xx_zz[i] * a_exp;
    }
    // integrals block (1554-1560)

    #pragma omp simd aligned(g_yz_xy_xy_xx, g_yz_xy_xy_xy, g_yz_xy_xy_xz, g_yz_xy_xy_yy, g_yz_xy_xy_yz, g_yz_xy_xy_zz, g_z_0_0_0_y_xy_xy_xx, g_z_0_0_0_y_xy_xy_xy, g_z_0_0_0_y_xy_xy_xz, g_z_0_0_0_y_xy_xy_yy, g_z_0_0_0_y_xy_xy_yz, g_z_0_0_0_y_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xy_xy_xx[i] = 2.0 * g_yz_xy_xy_xx[i] * a_exp;

        g_z_0_0_0_y_xy_xy_xy[i] = 2.0 * g_yz_xy_xy_xy[i] * a_exp;

        g_z_0_0_0_y_xy_xy_xz[i] = 2.0 * g_yz_xy_xy_xz[i] * a_exp;

        g_z_0_0_0_y_xy_xy_yy[i] = 2.0 * g_yz_xy_xy_yy[i] * a_exp;

        g_z_0_0_0_y_xy_xy_yz[i] = 2.0 * g_yz_xy_xy_yz[i] * a_exp;

        g_z_0_0_0_y_xy_xy_zz[i] = 2.0 * g_yz_xy_xy_zz[i] * a_exp;
    }
    // integrals block (1560-1566)

    #pragma omp simd aligned(g_yz_xy_xz_xx, g_yz_xy_xz_xy, g_yz_xy_xz_xz, g_yz_xy_xz_yy, g_yz_xy_xz_yz, g_yz_xy_xz_zz, g_z_0_0_0_y_xy_xz_xx, g_z_0_0_0_y_xy_xz_xy, g_z_0_0_0_y_xy_xz_xz, g_z_0_0_0_y_xy_xz_yy, g_z_0_0_0_y_xy_xz_yz, g_z_0_0_0_y_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xy_xz_xx[i] = 2.0 * g_yz_xy_xz_xx[i] * a_exp;

        g_z_0_0_0_y_xy_xz_xy[i] = 2.0 * g_yz_xy_xz_xy[i] * a_exp;

        g_z_0_0_0_y_xy_xz_xz[i] = 2.0 * g_yz_xy_xz_xz[i] * a_exp;

        g_z_0_0_0_y_xy_xz_yy[i] = 2.0 * g_yz_xy_xz_yy[i] * a_exp;

        g_z_0_0_0_y_xy_xz_yz[i] = 2.0 * g_yz_xy_xz_yz[i] * a_exp;

        g_z_0_0_0_y_xy_xz_zz[i] = 2.0 * g_yz_xy_xz_zz[i] * a_exp;
    }
    // integrals block (1566-1572)

    #pragma omp simd aligned(g_yz_xy_yy_xx, g_yz_xy_yy_xy, g_yz_xy_yy_xz, g_yz_xy_yy_yy, g_yz_xy_yy_yz, g_yz_xy_yy_zz, g_z_0_0_0_y_xy_yy_xx, g_z_0_0_0_y_xy_yy_xy, g_z_0_0_0_y_xy_yy_xz, g_z_0_0_0_y_xy_yy_yy, g_z_0_0_0_y_xy_yy_yz, g_z_0_0_0_y_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xy_yy_xx[i] = 2.0 * g_yz_xy_yy_xx[i] * a_exp;

        g_z_0_0_0_y_xy_yy_xy[i] = 2.0 * g_yz_xy_yy_xy[i] * a_exp;

        g_z_0_0_0_y_xy_yy_xz[i] = 2.0 * g_yz_xy_yy_xz[i] * a_exp;

        g_z_0_0_0_y_xy_yy_yy[i] = 2.0 * g_yz_xy_yy_yy[i] * a_exp;

        g_z_0_0_0_y_xy_yy_yz[i] = 2.0 * g_yz_xy_yy_yz[i] * a_exp;

        g_z_0_0_0_y_xy_yy_zz[i] = 2.0 * g_yz_xy_yy_zz[i] * a_exp;
    }
    // integrals block (1572-1578)

    #pragma omp simd aligned(g_yz_xy_yz_xx, g_yz_xy_yz_xy, g_yz_xy_yz_xz, g_yz_xy_yz_yy, g_yz_xy_yz_yz, g_yz_xy_yz_zz, g_z_0_0_0_y_xy_yz_xx, g_z_0_0_0_y_xy_yz_xy, g_z_0_0_0_y_xy_yz_xz, g_z_0_0_0_y_xy_yz_yy, g_z_0_0_0_y_xy_yz_yz, g_z_0_0_0_y_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xy_yz_xx[i] = 2.0 * g_yz_xy_yz_xx[i] * a_exp;

        g_z_0_0_0_y_xy_yz_xy[i] = 2.0 * g_yz_xy_yz_xy[i] * a_exp;

        g_z_0_0_0_y_xy_yz_xz[i] = 2.0 * g_yz_xy_yz_xz[i] * a_exp;

        g_z_0_0_0_y_xy_yz_yy[i] = 2.0 * g_yz_xy_yz_yy[i] * a_exp;

        g_z_0_0_0_y_xy_yz_yz[i] = 2.0 * g_yz_xy_yz_yz[i] * a_exp;

        g_z_0_0_0_y_xy_yz_zz[i] = 2.0 * g_yz_xy_yz_zz[i] * a_exp;
    }
    // integrals block (1578-1584)

    #pragma omp simd aligned(g_yz_xy_zz_xx, g_yz_xy_zz_xy, g_yz_xy_zz_xz, g_yz_xy_zz_yy, g_yz_xy_zz_yz, g_yz_xy_zz_zz, g_z_0_0_0_y_xy_zz_xx, g_z_0_0_0_y_xy_zz_xy, g_z_0_0_0_y_xy_zz_xz, g_z_0_0_0_y_xy_zz_yy, g_z_0_0_0_y_xy_zz_yz, g_z_0_0_0_y_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xy_zz_xx[i] = 2.0 * g_yz_xy_zz_xx[i] * a_exp;

        g_z_0_0_0_y_xy_zz_xy[i] = 2.0 * g_yz_xy_zz_xy[i] * a_exp;

        g_z_0_0_0_y_xy_zz_xz[i] = 2.0 * g_yz_xy_zz_xz[i] * a_exp;

        g_z_0_0_0_y_xy_zz_yy[i] = 2.0 * g_yz_xy_zz_yy[i] * a_exp;

        g_z_0_0_0_y_xy_zz_yz[i] = 2.0 * g_yz_xy_zz_yz[i] * a_exp;

        g_z_0_0_0_y_xy_zz_zz[i] = 2.0 * g_yz_xy_zz_zz[i] * a_exp;
    }
    // integrals block (1584-1590)

    #pragma omp simd aligned(g_yz_xz_xx_xx, g_yz_xz_xx_xy, g_yz_xz_xx_xz, g_yz_xz_xx_yy, g_yz_xz_xx_yz, g_yz_xz_xx_zz, g_z_0_0_0_y_xz_xx_xx, g_z_0_0_0_y_xz_xx_xy, g_z_0_0_0_y_xz_xx_xz, g_z_0_0_0_y_xz_xx_yy, g_z_0_0_0_y_xz_xx_yz, g_z_0_0_0_y_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xz_xx_xx[i] = 2.0 * g_yz_xz_xx_xx[i] * a_exp;

        g_z_0_0_0_y_xz_xx_xy[i] = 2.0 * g_yz_xz_xx_xy[i] * a_exp;

        g_z_0_0_0_y_xz_xx_xz[i] = 2.0 * g_yz_xz_xx_xz[i] * a_exp;

        g_z_0_0_0_y_xz_xx_yy[i] = 2.0 * g_yz_xz_xx_yy[i] * a_exp;

        g_z_0_0_0_y_xz_xx_yz[i] = 2.0 * g_yz_xz_xx_yz[i] * a_exp;

        g_z_0_0_0_y_xz_xx_zz[i] = 2.0 * g_yz_xz_xx_zz[i] * a_exp;
    }
    // integrals block (1590-1596)

    #pragma omp simd aligned(g_yz_xz_xy_xx, g_yz_xz_xy_xy, g_yz_xz_xy_xz, g_yz_xz_xy_yy, g_yz_xz_xy_yz, g_yz_xz_xy_zz, g_z_0_0_0_y_xz_xy_xx, g_z_0_0_0_y_xz_xy_xy, g_z_0_0_0_y_xz_xy_xz, g_z_0_0_0_y_xz_xy_yy, g_z_0_0_0_y_xz_xy_yz, g_z_0_0_0_y_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xz_xy_xx[i] = 2.0 * g_yz_xz_xy_xx[i] * a_exp;

        g_z_0_0_0_y_xz_xy_xy[i] = 2.0 * g_yz_xz_xy_xy[i] * a_exp;

        g_z_0_0_0_y_xz_xy_xz[i] = 2.0 * g_yz_xz_xy_xz[i] * a_exp;

        g_z_0_0_0_y_xz_xy_yy[i] = 2.0 * g_yz_xz_xy_yy[i] * a_exp;

        g_z_0_0_0_y_xz_xy_yz[i] = 2.0 * g_yz_xz_xy_yz[i] * a_exp;

        g_z_0_0_0_y_xz_xy_zz[i] = 2.0 * g_yz_xz_xy_zz[i] * a_exp;
    }
    // integrals block (1596-1602)

    #pragma omp simd aligned(g_yz_xz_xz_xx, g_yz_xz_xz_xy, g_yz_xz_xz_xz, g_yz_xz_xz_yy, g_yz_xz_xz_yz, g_yz_xz_xz_zz, g_z_0_0_0_y_xz_xz_xx, g_z_0_0_0_y_xz_xz_xy, g_z_0_0_0_y_xz_xz_xz, g_z_0_0_0_y_xz_xz_yy, g_z_0_0_0_y_xz_xz_yz, g_z_0_0_0_y_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xz_xz_xx[i] = 2.0 * g_yz_xz_xz_xx[i] * a_exp;

        g_z_0_0_0_y_xz_xz_xy[i] = 2.0 * g_yz_xz_xz_xy[i] * a_exp;

        g_z_0_0_0_y_xz_xz_xz[i] = 2.0 * g_yz_xz_xz_xz[i] * a_exp;

        g_z_0_0_0_y_xz_xz_yy[i] = 2.0 * g_yz_xz_xz_yy[i] * a_exp;

        g_z_0_0_0_y_xz_xz_yz[i] = 2.0 * g_yz_xz_xz_yz[i] * a_exp;

        g_z_0_0_0_y_xz_xz_zz[i] = 2.0 * g_yz_xz_xz_zz[i] * a_exp;
    }
    // integrals block (1602-1608)

    #pragma omp simd aligned(g_yz_xz_yy_xx, g_yz_xz_yy_xy, g_yz_xz_yy_xz, g_yz_xz_yy_yy, g_yz_xz_yy_yz, g_yz_xz_yy_zz, g_z_0_0_0_y_xz_yy_xx, g_z_0_0_0_y_xz_yy_xy, g_z_0_0_0_y_xz_yy_xz, g_z_0_0_0_y_xz_yy_yy, g_z_0_0_0_y_xz_yy_yz, g_z_0_0_0_y_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xz_yy_xx[i] = 2.0 * g_yz_xz_yy_xx[i] * a_exp;

        g_z_0_0_0_y_xz_yy_xy[i] = 2.0 * g_yz_xz_yy_xy[i] * a_exp;

        g_z_0_0_0_y_xz_yy_xz[i] = 2.0 * g_yz_xz_yy_xz[i] * a_exp;

        g_z_0_0_0_y_xz_yy_yy[i] = 2.0 * g_yz_xz_yy_yy[i] * a_exp;

        g_z_0_0_0_y_xz_yy_yz[i] = 2.0 * g_yz_xz_yy_yz[i] * a_exp;

        g_z_0_0_0_y_xz_yy_zz[i] = 2.0 * g_yz_xz_yy_zz[i] * a_exp;
    }
    // integrals block (1608-1614)

    #pragma omp simd aligned(g_yz_xz_yz_xx, g_yz_xz_yz_xy, g_yz_xz_yz_xz, g_yz_xz_yz_yy, g_yz_xz_yz_yz, g_yz_xz_yz_zz, g_z_0_0_0_y_xz_yz_xx, g_z_0_0_0_y_xz_yz_xy, g_z_0_0_0_y_xz_yz_xz, g_z_0_0_0_y_xz_yz_yy, g_z_0_0_0_y_xz_yz_yz, g_z_0_0_0_y_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xz_yz_xx[i] = 2.0 * g_yz_xz_yz_xx[i] * a_exp;

        g_z_0_0_0_y_xz_yz_xy[i] = 2.0 * g_yz_xz_yz_xy[i] * a_exp;

        g_z_0_0_0_y_xz_yz_xz[i] = 2.0 * g_yz_xz_yz_xz[i] * a_exp;

        g_z_0_0_0_y_xz_yz_yy[i] = 2.0 * g_yz_xz_yz_yy[i] * a_exp;

        g_z_0_0_0_y_xz_yz_yz[i] = 2.0 * g_yz_xz_yz_yz[i] * a_exp;

        g_z_0_0_0_y_xz_yz_zz[i] = 2.0 * g_yz_xz_yz_zz[i] * a_exp;
    }
    // integrals block (1614-1620)

    #pragma omp simd aligned(g_yz_xz_zz_xx, g_yz_xz_zz_xy, g_yz_xz_zz_xz, g_yz_xz_zz_yy, g_yz_xz_zz_yz, g_yz_xz_zz_zz, g_z_0_0_0_y_xz_zz_xx, g_z_0_0_0_y_xz_zz_xy, g_z_0_0_0_y_xz_zz_xz, g_z_0_0_0_y_xz_zz_yy, g_z_0_0_0_y_xz_zz_yz, g_z_0_0_0_y_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xz_zz_xx[i] = 2.0 * g_yz_xz_zz_xx[i] * a_exp;

        g_z_0_0_0_y_xz_zz_xy[i] = 2.0 * g_yz_xz_zz_xy[i] * a_exp;

        g_z_0_0_0_y_xz_zz_xz[i] = 2.0 * g_yz_xz_zz_xz[i] * a_exp;

        g_z_0_0_0_y_xz_zz_yy[i] = 2.0 * g_yz_xz_zz_yy[i] * a_exp;

        g_z_0_0_0_y_xz_zz_yz[i] = 2.0 * g_yz_xz_zz_yz[i] * a_exp;

        g_z_0_0_0_y_xz_zz_zz[i] = 2.0 * g_yz_xz_zz_zz[i] * a_exp;
    }
    // integrals block (1620-1626)

    #pragma omp simd aligned(g_yz_yy_xx_xx, g_yz_yy_xx_xy, g_yz_yy_xx_xz, g_yz_yy_xx_yy, g_yz_yy_xx_yz, g_yz_yy_xx_zz, g_z_0_0_0_y_yy_xx_xx, g_z_0_0_0_y_yy_xx_xy, g_z_0_0_0_y_yy_xx_xz, g_z_0_0_0_y_yy_xx_yy, g_z_0_0_0_y_yy_xx_yz, g_z_0_0_0_y_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yy_xx_xx[i] = 2.0 * g_yz_yy_xx_xx[i] * a_exp;

        g_z_0_0_0_y_yy_xx_xy[i] = 2.0 * g_yz_yy_xx_xy[i] * a_exp;

        g_z_0_0_0_y_yy_xx_xz[i] = 2.0 * g_yz_yy_xx_xz[i] * a_exp;

        g_z_0_0_0_y_yy_xx_yy[i] = 2.0 * g_yz_yy_xx_yy[i] * a_exp;

        g_z_0_0_0_y_yy_xx_yz[i] = 2.0 * g_yz_yy_xx_yz[i] * a_exp;

        g_z_0_0_0_y_yy_xx_zz[i] = 2.0 * g_yz_yy_xx_zz[i] * a_exp;
    }
    // integrals block (1626-1632)

    #pragma omp simd aligned(g_yz_yy_xy_xx, g_yz_yy_xy_xy, g_yz_yy_xy_xz, g_yz_yy_xy_yy, g_yz_yy_xy_yz, g_yz_yy_xy_zz, g_z_0_0_0_y_yy_xy_xx, g_z_0_0_0_y_yy_xy_xy, g_z_0_0_0_y_yy_xy_xz, g_z_0_0_0_y_yy_xy_yy, g_z_0_0_0_y_yy_xy_yz, g_z_0_0_0_y_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yy_xy_xx[i] = 2.0 * g_yz_yy_xy_xx[i] * a_exp;

        g_z_0_0_0_y_yy_xy_xy[i] = 2.0 * g_yz_yy_xy_xy[i] * a_exp;

        g_z_0_0_0_y_yy_xy_xz[i] = 2.0 * g_yz_yy_xy_xz[i] * a_exp;

        g_z_0_0_0_y_yy_xy_yy[i] = 2.0 * g_yz_yy_xy_yy[i] * a_exp;

        g_z_0_0_0_y_yy_xy_yz[i] = 2.0 * g_yz_yy_xy_yz[i] * a_exp;

        g_z_0_0_0_y_yy_xy_zz[i] = 2.0 * g_yz_yy_xy_zz[i] * a_exp;
    }
    // integrals block (1632-1638)

    #pragma omp simd aligned(g_yz_yy_xz_xx, g_yz_yy_xz_xy, g_yz_yy_xz_xz, g_yz_yy_xz_yy, g_yz_yy_xz_yz, g_yz_yy_xz_zz, g_z_0_0_0_y_yy_xz_xx, g_z_0_0_0_y_yy_xz_xy, g_z_0_0_0_y_yy_xz_xz, g_z_0_0_0_y_yy_xz_yy, g_z_0_0_0_y_yy_xz_yz, g_z_0_0_0_y_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yy_xz_xx[i] = 2.0 * g_yz_yy_xz_xx[i] * a_exp;

        g_z_0_0_0_y_yy_xz_xy[i] = 2.0 * g_yz_yy_xz_xy[i] * a_exp;

        g_z_0_0_0_y_yy_xz_xz[i] = 2.0 * g_yz_yy_xz_xz[i] * a_exp;

        g_z_0_0_0_y_yy_xz_yy[i] = 2.0 * g_yz_yy_xz_yy[i] * a_exp;

        g_z_0_0_0_y_yy_xz_yz[i] = 2.0 * g_yz_yy_xz_yz[i] * a_exp;

        g_z_0_0_0_y_yy_xz_zz[i] = 2.0 * g_yz_yy_xz_zz[i] * a_exp;
    }
    // integrals block (1638-1644)

    #pragma omp simd aligned(g_yz_yy_yy_xx, g_yz_yy_yy_xy, g_yz_yy_yy_xz, g_yz_yy_yy_yy, g_yz_yy_yy_yz, g_yz_yy_yy_zz, g_z_0_0_0_y_yy_yy_xx, g_z_0_0_0_y_yy_yy_xy, g_z_0_0_0_y_yy_yy_xz, g_z_0_0_0_y_yy_yy_yy, g_z_0_0_0_y_yy_yy_yz, g_z_0_0_0_y_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yy_yy_xx[i] = 2.0 * g_yz_yy_yy_xx[i] * a_exp;

        g_z_0_0_0_y_yy_yy_xy[i] = 2.0 * g_yz_yy_yy_xy[i] * a_exp;

        g_z_0_0_0_y_yy_yy_xz[i] = 2.0 * g_yz_yy_yy_xz[i] * a_exp;

        g_z_0_0_0_y_yy_yy_yy[i] = 2.0 * g_yz_yy_yy_yy[i] * a_exp;

        g_z_0_0_0_y_yy_yy_yz[i] = 2.0 * g_yz_yy_yy_yz[i] * a_exp;

        g_z_0_0_0_y_yy_yy_zz[i] = 2.0 * g_yz_yy_yy_zz[i] * a_exp;
    }
    // integrals block (1644-1650)

    #pragma omp simd aligned(g_yz_yy_yz_xx, g_yz_yy_yz_xy, g_yz_yy_yz_xz, g_yz_yy_yz_yy, g_yz_yy_yz_yz, g_yz_yy_yz_zz, g_z_0_0_0_y_yy_yz_xx, g_z_0_0_0_y_yy_yz_xy, g_z_0_0_0_y_yy_yz_xz, g_z_0_0_0_y_yy_yz_yy, g_z_0_0_0_y_yy_yz_yz, g_z_0_0_0_y_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yy_yz_xx[i] = 2.0 * g_yz_yy_yz_xx[i] * a_exp;

        g_z_0_0_0_y_yy_yz_xy[i] = 2.0 * g_yz_yy_yz_xy[i] * a_exp;

        g_z_0_0_0_y_yy_yz_xz[i] = 2.0 * g_yz_yy_yz_xz[i] * a_exp;

        g_z_0_0_0_y_yy_yz_yy[i] = 2.0 * g_yz_yy_yz_yy[i] * a_exp;

        g_z_0_0_0_y_yy_yz_yz[i] = 2.0 * g_yz_yy_yz_yz[i] * a_exp;

        g_z_0_0_0_y_yy_yz_zz[i] = 2.0 * g_yz_yy_yz_zz[i] * a_exp;
    }
    // integrals block (1650-1656)

    #pragma omp simd aligned(g_yz_yy_zz_xx, g_yz_yy_zz_xy, g_yz_yy_zz_xz, g_yz_yy_zz_yy, g_yz_yy_zz_yz, g_yz_yy_zz_zz, g_z_0_0_0_y_yy_zz_xx, g_z_0_0_0_y_yy_zz_xy, g_z_0_0_0_y_yy_zz_xz, g_z_0_0_0_y_yy_zz_yy, g_z_0_0_0_y_yy_zz_yz, g_z_0_0_0_y_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yy_zz_xx[i] = 2.0 * g_yz_yy_zz_xx[i] * a_exp;

        g_z_0_0_0_y_yy_zz_xy[i] = 2.0 * g_yz_yy_zz_xy[i] * a_exp;

        g_z_0_0_0_y_yy_zz_xz[i] = 2.0 * g_yz_yy_zz_xz[i] * a_exp;

        g_z_0_0_0_y_yy_zz_yy[i] = 2.0 * g_yz_yy_zz_yy[i] * a_exp;

        g_z_0_0_0_y_yy_zz_yz[i] = 2.0 * g_yz_yy_zz_yz[i] * a_exp;

        g_z_0_0_0_y_yy_zz_zz[i] = 2.0 * g_yz_yy_zz_zz[i] * a_exp;
    }
    // integrals block (1656-1662)

    #pragma omp simd aligned(g_yz_yz_xx_xx, g_yz_yz_xx_xy, g_yz_yz_xx_xz, g_yz_yz_xx_yy, g_yz_yz_xx_yz, g_yz_yz_xx_zz, g_z_0_0_0_y_yz_xx_xx, g_z_0_0_0_y_yz_xx_xy, g_z_0_0_0_y_yz_xx_xz, g_z_0_0_0_y_yz_xx_yy, g_z_0_0_0_y_yz_xx_yz, g_z_0_0_0_y_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yz_xx_xx[i] = 2.0 * g_yz_yz_xx_xx[i] * a_exp;

        g_z_0_0_0_y_yz_xx_xy[i] = 2.0 * g_yz_yz_xx_xy[i] * a_exp;

        g_z_0_0_0_y_yz_xx_xz[i] = 2.0 * g_yz_yz_xx_xz[i] * a_exp;

        g_z_0_0_0_y_yz_xx_yy[i] = 2.0 * g_yz_yz_xx_yy[i] * a_exp;

        g_z_0_0_0_y_yz_xx_yz[i] = 2.0 * g_yz_yz_xx_yz[i] * a_exp;

        g_z_0_0_0_y_yz_xx_zz[i] = 2.0 * g_yz_yz_xx_zz[i] * a_exp;
    }
    // integrals block (1662-1668)

    #pragma omp simd aligned(g_yz_yz_xy_xx, g_yz_yz_xy_xy, g_yz_yz_xy_xz, g_yz_yz_xy_yy, g_yz_yz_xy_yz, g_yz_yz_xy_zz, g_z_0_0_0_y_yz_xy_xx, g_z_0_0_0_y_yz_xy_xy, g_z_0_0_0_y_yz_xy_xz, g_z_0_0_0_y_yz_xy_yy, g_z_0_0_0_y_yz_xy_yz, g_z_0_0_0_y_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yz_xy_xx[i] = 2.0 * g_yz_yz_xy_xx[i] * a_exp;

        g_z_0_0_0_y_yz_xy_xy[i] = 2.0 * g_yz_yz_xy_xy[i] * a_exp;

        g_z_0_0_0_y_yz_xy_xz[i] = 2.0 * g_yz_yz_xy_xz[i] * a_exp;

        g_z_0_0_0_y_yz_xy_yy[i] = 2.0 * g_yz_yz_xy_yy[i] * a_exp;

        g_z_0_0_0_y_yz_xy_yz[i] = 2.0 * g_yz_yz_xy_yz[i] * a_exp;

        g_z_0_0_0_y_yz_xy_zz[i] = 2.0 * g_yz_yz_xy_zz[i] * a_exp;
    }
    // integrals block (1668-1674)

    #pragma omp simd aligned(g_yz_yz_xz_xx, g_yz_yz_xz_xy, g_yz_yz_xz_xz, g_yz_yz_xz_yy, g_yz_yz_xz_yz, g_yz_yz_xz_zz, g_z_0_0_0_y_yz_xz_xx, g_z_0_0_0_y_yz_xz_xy, g_z_0_0_0_y_yz_xz_xz, g_z_0_0_0_y_yz_xz_yy, g_z_0_0_0_y_yz_xz_yz, g_z_0_0_0_y_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yz_xz_xx[i] = 2.0 * g_yz_yz_xz_xx[i] * a_exp;

        g_z_0_0_0_y_yz_xz_xy[i] = 2.0 * g_yz_yz_xz_xy[i] * a_exp;

        g_z_0_0_0_y_yz_xz_xz[i] = 2.0 * g_yz_yz_xz_xz[i] * a_exp;

        g_z_0_0_0_y_yz_xz_yy[i] = 2.0 * g_yz_yz_xz_yy[i] * a_exp;

        g_z_0_0_0_y_yz_xz_yz[i] = 2.0 * g_yz_yz_xz_yz[i] * a_exp;

        g_z_0_0_0_y_yz_xz_zz[i] = 2.0 * g_yz_yz_xz_zz[i] * a_exp;
    }
    // integrals block (1674-1680)

    #pragma omp simd aligned(g_yz_yz_yy_xx, g_yz_yz_yy_xy, g_yz_yz_yy_xz, g_yz_yz_yy_yy, g_yz_yz_yy_yz, g_yz_yz_yy_zz, g_z_0_0_0_y_yz_yy_xx, g_z_0_0_0_y_yz_yy_xy, g_z_0_0_0_y_yz_yy_xz, g_z_0_0_0_y_yz_yy_yy, g_z_0_0_0_y_yz_yy_yz, g_z_0_0_0_y_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yz_yy_xx[i] = 2.0 * g_yz_yz_yy_xx[i] * a_exp;

        g_z_0_0_0_y_yz_yy_xy[i] = 2.0 * g_yz_yz_yy_xy[i] * a_exp;

        g_z_0_0_0_y_yz_yy_xz[i] = 2.0 * g_yz_yz_yy_xz[i] * a_exp;

        g_z_0_0_0_y_yz_yy_yy[i] = 2.0 * g_yz_yz_yy_yy[i] * a_exp;

        g_z_0_0_0_y_yz_yy_yz[i] = 2.0 * g_yz_yz_yy_yz[i] * a_exp;

        g_z_0_0_0_y_yz_yy_zz[i] = 2.0 * g_yz_yz_yy_zz[i] * a_exp;
    }
    // integrals block (1680-1686)

    #pragma omp simd aligned(g_yz_yz_yz_xx, g_yz_yz_yz_xy, g_yz_yz_yz_xz, g_yz_yz_yz_yy, g_yz_yz_yz_yz, g_yz_yz_yz_zz, g_z_0_0_0_y_yz_yz_xx, g_z_0_0_0_y_yz_yz_xy, g_z_0_0_0_y_yz_yz_xz, g_z_0_0_0_y_yz_yz_yy, g_z_0_0_0_y_yz_yz_yz, g_z_0_0_0_y_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yz_yz_xx[i] = 2.0 * g_yz_yz_yz_xx[i] * a_exp;

        g_z_0_0_0_y_yz_yz_xy[i] = 2.0 * g_yz_yz_yz_xy[i] * a_exp;

        g_z_0_0_0_y_yz_yz_xz[i] = 2.0 * g_yz_yz_yz_xz[i] * a_exp;

        g_z_0_0_0_y_yz_yz_yy[i] = 2.0 * g_yz_yz_yz_yy[i] * a_exp;

        g_z_0_0_0_y_yz_yz_yz[i] = 2.0 * g_yz_yz_yz_yz[i] * a_exp;

        g_z_0_0_0_y_yz_yz_zz[i] = 2.0 * g_yz_yz_yz_zz[i] * a_exp;
    }
    // integrals block (1686-1692)

    #pragma omp simd aligned(g_yz_yz_zz_xx, g_yz_yz_zz_xy, g_yz_yz_zz_xz, g_yz_yz_zz_yy, g_yz_yz_zz_yz, g_yz_yz_zz_zz, g_z_0_0_0_y_yz_zz_xx, g_z_0_0_0_y_yz_zz_xy, g_z_0_0_0_y_yz_zz_xz, g_z_0_0_0_y_yz_zz_yy, g_z_0_0_0_y_yz_zz_yz, g_z_0_0_0_y_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yz_zz_xx[i] = 2.0 * g_yz_yz_zz_xx[i] * a_exp;

        g_z_0_0_0_y_yz_zz_xy[i] = 2.0 * g_yz_yz_zz_xy[i] * a_exp;

        g_z_0_0_0_y_yz_zz_xz[i] = 2.0 * g_yz_yz_zz_xz[i] * a_exp;

        g_z_0_0_0_y_yz_zz_yy[i] = 2.0 * g_yz_yz_zz_yy[i] * a_exp;

        g_z_0_0_0_y_yz_zz_yz[i] = 2.0 * g_yz_yz_zz_yz[i] * a_exp;

        g_z_0_0_0_y_yz_zz_zz[i] = 2.0 * g_yz_yz_zz_zz[i] * a_exp;
    }
    // integrals block (1692-1698)

    #pragma omp simd aligned(g_yz_zz_xx_xx, g_yz_zz_xx_xy, g_yz_zz_xx_xz, g_yz_zz_xx_yy, g_yz_zz_xx_yz, g_yz_zz_xx_zz, g_z_0_0_0_y_zz_xx_xx, g_z_0_0_0_y_zz_xx_xy, g_z_0_0_0_y_zz_xx_xz, g_z_0_0_0_y_zz_xx_yy, g_z_0_0_0_y_zz_xx_yz, g_z_0_0_0_y_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_zz_xx_xx[i] = 2.0 * g_yz_zz_xx_xx[i] * a_exp;

        g_z_0_0_0_y_zz_xx_xy[i] = 2.0 * g_yz_zz_xx_xy[i] * a_exp;

        g_z_0_0_0_y_zz_xx_xz[i] = 2.0 * g_yz_zz_xx_xz[i] * a_exp;

        g_z_0_0_0_y_zz_xx_yy[i] = 2.0 * g_yz_zz_xx_yy[i] * a_exp;

        g_z_0_0_0_y_zz_xx_yz[i] = 2.0 * g_yz_zz_xx_yz[i] * a_exp;

        g_z_0_0_0_y_zz_xx_zz[i] = 2.0 * g_yz_zz_xx_zz[i] * a_exp;
    }
    // integrals block (1698-1704)

    #pragma omp simd aligned(g_yz_zz_xy_xx, g_yz_zz_xy_xy, g_yz_zz_xy_xz, g_yz_zz_xy_yy, g_yz_zz_xy_yz, g_yz_zz_xy_zz, g_z_0_0_0_y_zz_xy_xx, g_z_0_0_0_y_zz_xy_xy, g_z_0_0_0_y_zz_xy_xz, g_z_0_0_0_y_zz_xy_yy, g_z_0_0_0_y_zz_xy_yz, g_z_0_0_0_y_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_zz_xy_xx[i] = 2.0 * g_yz_zz_xy_xx[i] * a_exp;

        g_z_0_0_0_y_zz_xy_xy[i] = 2.0 * g_yz_zz_xy_xy[i] * a_exp;

        g_z_0_0_0_y_zz_xy_xz[i] = 2.0 * g_yz_zz_xy_xz[i] * a_exp;

        g_z_0_0_0_y_zz_xy_yy[i] = 2.0 * g_yz_zz_xy_yy[i] * a_exp;

        g_z_0_0_0_y_zz_xy_yz[i] = 2.0 * g_yz_zz_xy_yz[i] * a_exp;

        g_z_0_0_0_y_zz_xy_zz[i] = 2.0 * g_yz_zz_xy_zz[i] * a_exp;
    }
    // integrals block (1704-1710)

    #pragma omp simd aligned(g_yz_zz_xz_xx, g_yz_zz_xz_xy, g_yz_zz_xz_xz, g_yz_zz_xz_yy, g_yz_zz_xz_yz, g_yz_zz_xz_zz, g_z_0_0_0_y_zz_xz_xx, g_z_0_0_0_y_zz_xz_xy, g_z_0_0_0_y_zz_xz_xz, g_z_0_0_0_y_zz_xz_yy, g_z_0_0_0_y_zz_xz_yz, g_z_0_0_0_y_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_zz_xz_xx[i] = 2.0 * g_yz_zz_xz_xx[i] * a_exp;

        g_z_0_0_0_y_zz_xz_xy[i] = 2.0 * g_yz_zz_xz_xy[i] * a_exp;

        g_z_0_0_0_y_zz_xz_xz[i] = 2.0 * g_yz_zz_xz_xz[i] * a_exp;

        g_z_0_0_0_y_zz_xz_yy[i] = 2.0 * g_yz_zz_xz_yy[i] * a_exp;

        g_z_0_0_0_y_zz_xz_yz[i] = 2.0 * g_yz_zz_xz_yz[i] * a_exp;

        g_z_0_0_0_y_zz_xz_zz[i] = 2.0 * g_yz_zz_xz_zz[i] * a_exp;
    }
    // integrals block (1710-1716)

    #pragma omp simd aligned(g_yz_zz_yy_xx, g_yz_zz_yy_xy, g_yz_zz_yy_xz, g_yz_zz_yy_yy, g_yz_zz_yy_yz, g_yz_zz_yy_zz, g_z_0_0_0_y_zz_yy_xx, g_z_0_0_0_y_zz_yy_xy, g_z_0_0_0_y_zz_yy_xz, g_z_0_0_0_y_zz_yy_yy, g_z_0_0_0_y_zz_yy_yz, g_z_0_0_0_y_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_zz_yy_xx[i] = 2.0 * g_yz_zz_yy_xx[i] * a_exp;

        g_z_0_0_0_y_zz_yy_xy[i] = 2.0 * g_yz_zz_yy_xy[i] * a_exp;

        g_z_0_0_0_y_zz_yy_xz[i] = 2.0 * g_yz_zz_yy_xz[i] * a_exp;

        g_z_0_0_0_y_zz_yy_yy[i] = 2.0 * g_yz_zz_yy_yy[i] * a_exp;

        g_z_0_0_0_y_zz_yy_yz[i] = 2.0 * g_yz_zz_yy_yz[i] * a_exp;

        g_z_0_0_0_y_zz_yy_zz[i] = 2.0 * g_yz_zz_yy_zz[i] * a_exp;
    }
    // integrals block (1716-1722)

    #pragma omp simd aligned(g_yz_zz_yz_xx, g_yz_zz_yz_xy, g_yz_zz_yz_xz, g_yz_zz_yz_yy, g_yz_zz_yz_yz, g_yz_zz_yz_zz, g_z_0_0_0_y_zz_yz_xx, g_z_0_0_0_y_zz_yz_xy, g_z_0_0_0_y_zz_yz_xz, g_z_0_0_0_y_zz_yz_yy, g_z_0_0_0_y_zz_yz_yz, g_z_0_0_0_y_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_zz_yz_xx[i] = 2.0 * g_yz_zz_yz_xx[i] * a_exp;

        g_z_0_0_0_y_zz_yz_xy[i] = 2.0 * g_yz_zz_yz_xy[i] * a_exp;

        g_z_0_0_0_y_zz_yz_xz[i] = 2.0 * g_yz_zz_yz_xz[i] * a_exp;

        g_z_0_0_0_y_zz_yz_yy[i] = 2.0 * g_yz_zz_yz_yy[i] * a_exp;

        g_z_0_0_0_y_zz_yz_yz[i] = 2.0 * g_yz_zz_yz_yz[i] * a_exp;

        g_z_0_0_0_y_zz_yz_zz[i] = 2.0 * g_yz_zz_yz_zz[i] * a_exp;
    }
    // integrals block (1722-1728)

    #pragma omp simd aligned(g_yz_zz_zz_xx, g_yz_zz_zz_xy, g_yz_zz_zz_xz, g_yz_zz_zz_yy, g_yz_zz_zz_yz, g_yz_zz_zz_zz, g_z_0_0_0_y_zz_zz_xx, g_z_0_0_0_y_zz_zz_xy, g_z_0_0_0_y_zz_zz_xz, g_z_0_0_0_y_zz_zz_yy, g_z_0_0_0_y_zz_zz_yz, g_z_0_0_0_y_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_zz_zz_xx[i] = 2.0 * g_yz_zz_zz_xx[i] * a_exp;

        g_z_0_0_0_y_zz_zz_xy[i] = 2.0 * g_yz_zz_zz_xy[i] * a_exp;

        g_z_0_0_0_y_zz_zz_xz[i] = 2.0 * g_yz_zz_zz_xz[i] * a_exp;

        g_z_0_0_0_y_zz_zz_yy[i] = 2.0 * g_yz_zz_zz_yy[i] * a_exp;

        g_z_0_0_0_y_zz_zz_yz[i] = 2.0 * g_yz_zz_zz_yz[i] * a_exp;

        g_z_0_0_0_y_zz_zz_zz[i] = 2.0 * g_yz_zz_zz_zz[i] * a_exp;
    }
    // integrals block (1728-1734)

    #pragma omp simd aligned(g_0_xx_xx_xx, g_0_xx_xx_xy, g_0_xx_xx_xz, g_0_xx_xx_yy, g_0_xx_xx_yz, g_0_xx_xx_zz, g_z_0_0_0_z_xx_xx_xx, g_z_0_0_0_z_xx_xx_xy, g_z_0_0_0_z_xx_xx_xz, g_z_0_0_0_z_xx_xx_yy, g_z_0_0_0_z_xx_xx_yz, g_z_0_0_0_z_xx_xx_zz, g_zz_xx_xx_xx, g_zz_xx_xx_xy, g_zz_xx_xx_xz, g_zz_xx_xx_yy, g_zz_xx_xx_yz, g_zz_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xx_xx_xx[i] = -g_0_xx_xx_xx[i] + 2.0 * g_zz_xx_xx_xx[i] * a_exp;

        g_z_0_0_0_z_xx_xx_xy[i] = -g_0_xx_xx_xy[i] + 2.0 * g_zz_xx_xx_xy[i] * a_exp;

        g_z_0_0_0_z_xx_xx_xz[i] = -g_0_xx_xx_xz[i] + 2.0 * g_zz_xx_xx_xz[i] * a_exp;

        g_z_0_0_0_z_xx_xx_yy[i] = -g_0_xx_xx_yy[i] + 2.0 * g_zz_xx_xx_yy[i] * a_exp;

        g_z_0_0_0_z_xx_xx_yz[i] = -g_0_xx_xx_yz[i] + 2.0 * g_zz_xx_xx_yz[i] * a_exp;

        g_z_0_0_0_z_xx_xx_zz[i] = -g_0_xx_xx_zz[i] + 2.0 * g_zz_xx_xx_zz[i] * a_exp;
    }
    // integrals block (1734-1740)

    #pragma omp simd aligned(g_0_xx_xy_xx, g_0_xx_xy_xy, g_0_xx_xy_xz, g_0_xx_xy_yy, g_0_xx_xy_yz, g_0_xx_xy_zz, g_z_0_0_0_z_xx_xy_xx, g_z_0_0_0_z_xx_xy_xy, g_z_0_0_0_z_xx_xy_xz, g_z_0_0_0_z_xx_xy_yy, g_z_0_0_0_z_xx_xy_yz, g_z_0_0_0_z_xx_xy_zz, g_zz_xx_xy_xx, g_zz_xx_xy_xy, g_zz_xx_xy_xz, g_zz_xx_xy_yy, g_zz_xx_xy_yz, g_zz_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xx_xy_xx[i] = -g_0_xx_xy_xx[i] + 2.0 * g_zz_xx_xy_xx[i] * a_exp;

        g_z_0_0_0_z_xx_xy_xy[i] = -g_0_xx_xy_xy[i] + 2.0 * g_zz_xx_xy_xy[i] * a_exp;

        g_z_0_0_0_z_xx_xy_xz[i] = -g_0_xx_xy_xz[i] + 2.0 * g_zz_xx_xy_xz[i] * a_exp;

        g_z_0_0_0_z_xx_xy_yy[i] = -g_0_xx_xy_yy[i] + 2.0 * g_zz_xx_xy_yy[i] * a_exp;

        g_z_0_0_0_z_xx_xy_yz[i] = -g_0_xx_xy_yz[i] + 2.0 * g_zz_xx_xy_yz[i] * a_exp;

        g_z_0_0_0_z_xx_xy_zz[i] = -g_0_xx_xy_zz[i] + 2.0 * g_zz_xx_xy_zz[i] * a_exp;
    }
    // integrals block (1740-1746)

    #pragma omp simd aligned(g_0_xx_xz_xx, g_0_xx_xz_xy, g_0_xx_xz_xz, g_0_xx_xz_yy, g_0_xx_xz_yz, g_0_xx_xz_zz, g_z_0_0_0_z_xx_xz_xx, g_z_0_0_0_z_xx_xz_xy, g_z_0_0_0_z_xx_xz_xz, g_z_0_0_0_z_xx_xz_yy, g_z_0_0_0_z_xx_xz_yz, g_z_0_0_0_z_xx_xz_zz, g_zz_xx_xz_xx, g_zz_xx_xz_xy, g_zz_xx_xz_xz, g_zz_xx_xz_yy, g_zz_xx_xz_yz, g_zz_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xx_xz_xx[i] = -g_0_xx_xz_xx[i] + 2.0 * g_zz_xx_xz_xx[i] * a_exp;

        g_z_0_0_0_z_xx_xz_xy[i] = -g_0_xx_xz_xy[i] + 2.0 * g_zz_xx_xz_xy[i] * a_exp;

        g_z_0_0_0_z_xx_xz_xz[i] = -g_0_xx_xz_xz[i] + 2.0 * g_zz_xx_xz_xz[i] * a_exp;

        g_z_0_0_0_z_xx_xz_yy[i] = -g_0_xx_xz_yy[i] + 2.0 * g_zz_xx_xz_yy[i] * a_exp;

        g_z_0_0_0_z_xx_xz_yz[i] = -g_0_xx_xz_yz[i] + 2.0 * g_zz_xx_xz_yz[i] * a_exp;

        g_z_0_0_0_z_xx_xz_zz[i] = -g_0_xx_xz_zz[i] + 2.0 * g_zz_xx_xz_zz[i] * a_exp;
    }
    // integrals block (1746-1752)

    #pragma omp simd aligned(g_0_xx_yy_xx, g_0_xx_yy_xy, g_0_xx_yy_xz, g_0_xx_yy_yy, g_0_xx_yy_yz, g_0_xx_yy_zz, g_z_0_0_0_z_xx_yy_xx, g_z_0_0_0_z_xx_yy_xy, g_z_0_0_0_z_xx_yy_xz, g_z_0_0_0_z_xx_yy_yy, g_z_0_0_0_z_xx_yy_yz, g_z_0_0_0_z_xx_yy_zz, g_zz_xx_yy_xx, g_zz_xx_yy_xy, g_zz_xx_yy_xz, g_zz_xx_yy_yy, g_zz_xx_yy_yz, g_zz_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xx_yy_xx[i] = -g_0_xx_yy_xx[i] + 2.0 * g_zz_xx_yy_xx[i] * a_exp;

        g_z_0_0_0_z_xx_yy_xy[i] = -g_0_xx_yy_xy[i] + 2.0 * g_zz_xx_yy_xy[i] * a_exp;

        g_z_0_0_0_z_xx_yy_xz[i] = -g_0_xx_yy_xz[i] + 2.0 * g_zz_xx_yy_xz[i] * a_exp;

        g_z_0_0_0_z_xx_yy_yy[i] = -g_0_xx_yy_yy[i] + 2.0 * g_zz_xx_yy_yy[i] * a_exp;

        g_z_0_0_0_z_xx_yy_yz[i] = -g_0_xx_yy_yz[i] + 2.0 * g_zz_xx_yy_yz[i] * a_exp;

        g_z_0_0_0_z_xx_yy_zz[i] = -g_0_xx_yy_zz[i] + 2.0 * g_zz_xx_yy_zz[i] * a_exp;
    }
    // integrals block (1752-1758)

    #pragma omp simd aligned(g_0_xx_yz_xx, g_0_xx_yz_xy, g_0_xx_yz_xz, g_0_xx_yz_yy, g_0_xx_yz_yz, g_0_xx_yz_zz, g_z_0_0_0_z_xx_yz_xx, g_z_0_0_0_z_xx_yz_xy, g_z_0_0_0_z_xx_yz_xz, g_z_0_0_0_z_xx_yz_yy, g_z_0_0_0_z_xx_yz_yz, g_z_0_0_0_z_xx_yz_zz, g_zz_xx_yz_xx, g_zz_xx_yz_xy, g_zz_xx_yz_xz, g_zz_xx_yz_yy, g_zz_xx_yz_yz, g_zz_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xx_yz_xx[i] = -g_0_xx_yz_xx[i] + 2.0 * g_zz_xx_yz_xx[i] * a_exp;

        g_z_0_0_0_z_xx_yz_xy[i] = -g_0_xx_yz_xy[i] + 2.0 * g_zz_xx_yz_xy[i] * a_exp;

        g_z_0_0_0_z_xx_yz_xz[i] = -g_0_xx_yz_xz[i] + 2.0 * g_zz_xx_yz_xz[i] * a_exp;

        g_z_0_0_0_z_xx_yz_yy[i] = -g_0_xx_yz_yy[i] + 2.0 * g_zz_xx_yz_yy[i] * a_exp;

        g_z_0_0_0_z_xx_yz_yz[i] = -g_0_xx_yz_yz[i] + 2.0 * g_zz_xx_yz_yz[i] * a_exp;

        g_z_0_0_0_z_xx_yz_zz[i] = -g_0_xx_yz_zz[i] + 2.0 * g_zz_xx_yz_zz[i] * a_exp;
    }
    // integrals block (1758-1764)

    #pragma omp simd aligned(g_0_xx_zz_xx, g_0_xx_zz_xy, g_0_xx_zz_xz, g_0_xx_zz_yy, g_0_xx_zz_yz, g_0_xx_zz_zz, g_z_0_0_0_z_xx_zz_xx, g_z_0_0_0_z_xx_zz_xy, g_z_0_0_0_z_xx_zz_xz, g_z_0_0_0_z_xx_zz_yy, g_z_0_0_0_z_xx_zz_yz, g_z_0_0_0_z_xx_zz_zz, g_zz_xx_zz_xx, g_zz_xx_zz_xy, g_zz_xx_zz_xz, g_zz_xx_zz_yy, g_zz_xx_zz_yz, g_zz_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xx_zz_xx[i] = -g_0_xx_zz_xx[i] + 2.0 * g_zz_xx_zz_xx[i] * a_exp;

        g_z_0_0_0_z_xx_zz_xy[i] = -g_0_xx_zz_xy[i] + 2.0 * g_zz_xx_zz_xy[i] * a_exp;

        g_z_0_0_0_z_xx_zz_xz[i] = -g_0_xx_zz_xz[i] + 2.0 * g_zz_xx_zz_xz[i] * a_exp;

        g_z_0_0_0_z_xx_zz_yy[i] = -g_0_xx_zz_yy[i] + 2.0 * g_zz_xx_zz_yy[i] * a_exp;

        g_z_0_0_0_z_xx_zz_yz[i] = -g_0_xx_zz_yz[i] + 2.0 * g_zz_xx_zz_yz[i] * a_exp;

        g_z_0_0_0_z_xx_zz_zz[i] = -g_0_xx_zz_zz[i] + 2.0 * g_zz_xx_zz_zz[i] * a_exp;
    }
    // integrals block (1764-1770)

    #pragma omp simd aligned(g_0_xy_xx_xx, g_0_xy_xx_xy, g_0_xy_xx_xz, g_0_xy_xx_yy, g_0_xy_xx_yz, g_0_xy_xx_zz, g_z_0_0_0_z_xy_xx_xx, g_z_0_0_0_z_xy_xx_xy, g_z_0_0_0_z_xy_xx_xz, g_z_0_0_0_z_xy_xx_yy, g_z_0_0_0_z_xy_xx_yz, g_z_0_0_0_z_xy_xx_zz, g_zz_xy_xx_xx, g_zz_xy_xx_xy, g_zz_xy_xx_xz, g_zz_xy_xx_yy, g_zz_xy_xx_yz, g_zz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xy_xx_xx[i] = -g_0_xy_xx_xx[i] + 2.0 * g_zz_xy_xx_xx[i] * a_exp;

        g_z_0_0_0_z_xy_xx_xy[i] = -g_0_xy_xx_xy[i] + 2.0 * g_zz_xy_xx_xy[i] * a_exp;

        g_z_0_0_0_z_xy_xx_xz[i] = -g_0_xy_xx_xz[i] + 2.0 * g_zz_xy_xx_xz[i] * a_exp;

        g_z_0_0_0_z_xy_xx_yy[i] = -g_0_xy_xx_yy[i] + 2.0 * g_zz_xy_xx_yy[i] * a_exp;

        g_z_0_0_0_z_xy_xx_yz[i] = -g_0_xy_xx_yz[i] + 2.0 * g_zz_xy_xx_yz[i] * a_exp;

        g_z_0_0_0_z_xy_xx_zz[i] = -g_0_xy_xx_zz[i] + 2.0 * g_zz_xy_xx_zz[i] * a_exp;
    }
    // integrals block (1770-1776)

    #pragma omp simd aligned(g_0_xy_xy_xx, g_0_xy_xy_xy, g_0_xy_xy_xz, g_0_xy_xy_yy, g_0_xy_xy_yz, g_0_xy_xy_zz, g_z_0_0_0_z_xy_xy_xx, g_z_0_0_0_z_xy_xy_xy, g_z_0_0_0_z_xy_xy_xz, g_z_0_0_0_z_xy_xy_yy, g_z_0_0_0_z_xy_xy_yz, g_z_0_0_0_z_xy_xy_zz, g_zz_xy_xy_xx, g_zz_xy_xy_xy, g_zz_xy_xy_xz, g_zz_xy_xy_yy, g_zz_xy_xy_yz, g_zz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xy_xy_xx[i] = -g_0_xy_xy_xx[i] + 2.0 * g_zz_xy_xy_xx[i] * a_exp;

        g_z_0_0_0_z_xy_xy_xy[i] = -g_0_xy_xy_xy[i] + 2.0 * g_zz_xy_xy_xy[i] * a_exp;

        g_z_0_0_0_z_xy_xy_xz[i] = -g_0_xy_xy_xz[i] + 2.0 * g_zz_xy_xy_xz[i] * a_exp;

        g_z_0_0_0_z_xy_xy_yy[i] = -g_0_xy_xy_yy[i] + 2.0 * g_zz_xy_xy_yy[i] * a_exp;

        g_z_0_0_0_z_xy_xy_yz[i] = -g_0_xy_xy_yz[i] + 2.0 * g_zz_xy_xy_yz[i] * a_exp;

        g_z_0_0_0_z_xy_xy_zz[i] = -g_0_xy_xy_zz[i] + 2.0 * g_zz_xy_xy_zz[i] * a_exp;
    }
    // integrals block (1776-1782)

    #pragma omp simd aligned(g_0_xy_xz_xx, g_0_xy_xz_xy, g_0_xy_xz_xz, g_0_xy_xz_yy, g_0_xy_xz_yz, g_0_xy_xz_zz, g_z_0_0_0_z_xy_xz_xx, g_z_0_0_0_z_xy_xz_xy, g_z_0_0_0_z_xy_xz_xz, g_z_0_0_0_z_xy_xz_yy, g_z_0_0_0_z_xy_xz_yz, g_z_0_0_0_z_xy_xz_zz, g_zz_xy_xz_xx, g_zz_xy_xz_xy, g_zz_xy_xz_xz, g_zz_xy_xz_yy, g_zz_xy_xz_yz, g_zz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xy_xz_xx[i] = -g_0_xy_xz_xx[i] + 2.0 * g_zz_xy_xz_xx[i] * a_exp;

        g_z_0_0_0_z_xy_xz_xy[i] = -g_0_xy_xz_xy[i] + 2.0 * g_zz_xy_xz_xy[i] * a_exp;

        g_z_0_0_0_z_xy_xz_xz[i] = -g_0_xy_xz_xz[i] + 2.0 * g_zz_xy_xz_xz[i] * a_exp;

        g_z_0_0_0_z_xy_xz_yy[i] = -g_0_xy_xz_yy[i] + 2.0 * g_zz_xy_xz_yy[i] * a_exp;

        g_z_0_0_0_z_xy_xz_yz[i] = -g_0_xy_xz_yz[i] + 2.0 * g_zz_xy_xz_yz[i] * a_exp;

        g_z_0_0_0_z_xy_xz_zz[i] = -g_0_xy_xz_zz[i] + 2.0 * g_zz_xy_xz_zz[i] * a_exp;
    }
    // integrals block (1782-1788)

    #pragma omp simd aligned(g_0_xy_yy_xx, g_0_xy_yy_xy, g_0_xy_yy_xz, g_0_xy_yy_yy, g_0_xy_yy_yz, g_0_xy_yy_zz, g_z_0_0_0_z_xy_yy_xx, g_z_0_0_0_z_xy_yy_xy, g_z_0_0_0_z_xy_yy_xz, g_z_0_0_0_z_xy_yy_yy, g_z_0_0_0_z_xy_yy_yz, g_z_0_0_0_z_xy_yy_zz, g_zz_xy_yy_xx, g_zz_xy_yy_xy, g_zz_xy_yy_xz, g_zz_xy_yy_yy, g_zz_xy_yy_yz, g_zz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xy_yy_xx[i] = -g_0_xy_yy_xx[i] + 2.0 * g_zz_xy_yy_xx[i] * a_exp;

        g_z_0_0_0_z_xy_yy_xy[i] = -g_0_xy_yy_xy[i] + 2.0 * g_zz_xy_yy_xy[i] * a_exp;

        g_z_0_0_0_z_xy_yy_xz[i] = -g_0_xy_yy_xz[i] + 2.0 * g_zz_xy_yy_xz[i] * a_exp;

        g_z_0_0_0_z_xy_yy_yy[i] = -g_0_xy_yy_yy[i] + 2.0 * g_zz_xy_yy_yy[i] * a_exp;

        g_z_0_0_0_z_xy_yy_yz[i] = -g_0_xy_yy_yz[i] + 2.0 * g_zz_xy_yy_yz[i] * a_exp;

        g_z_0_0_0_z_xy_yy_zz[i] = -g_0_xy_yy_zz[i] + 2.0 * g_zz_xy_yy_zz[i] * a_exp;
    }
    // integrals block (1788-1794)

    #pragma omp simd aligned(g_0_xy_yz_xx, g_0_xy_yz_xy, g_0_xy_yz_xz, g_0_xy_yz_yy, g_0_xy_yz_yz, g_0_xy_yz_zz, g_z_0_0_0_z_xy_yz_xx, g_z_0_0_0_z_xy_yz_xy, g_z_0_0_0_z_xy_yz_xz, g_z_0_0_0_z_xy_yz_yy, g_z_0_0_0_z_xy_yz_yz, g_z_0_0_0_z_xy_yz_zz, g_zz_xy_yz_xx, g_zz_xy_yz_xy, g_zz_xy_yz_xz, g_zz_xy_yz_yy, g_zz_xy_yz_yz, g_zz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xy_yz_xx[i] = -g_0_xy_yz_xx[i] + 2.0 * g_zz_xy_yz_xx[i] * a_exp;

        g_z_0_0_0_z_xy_yz_xy[i] = -g_0_xy_yz_xy[i] + 2.0 * g_zz_xy_yz_xy[i] * a_exp;

        g_z_0_0_0_z_xy_yz_xz[i] = -g_0_xy_yz_xz[i] + 2.0 * g_zz_xy_yz_xz[i] * a_exp;

        g_z_0_0_0_z_xy_yz_yy[i] = -g_0_xy_yz_yy[i] + 2.0 * g_zz_xy_yz_yy[i] * a_exp;

        g_z_0_0_0_z_xy_yz_yz[i] = -g_0_xy_yz_yz[i] + 2.0 * g_zz_xy_yz_yz[i] * a_exp;

        g_z_0_0_0_z_xy_yz_zz[i] = -g_0_xy_yz_zz[i] + 2.0 * g_zz_xy_yz_zz[i] * a_exp;
    }
    // integrals block (1794-1800)

    #pragma omp simd aligned(g_0_xy_zz_xx, g_0_xy_zz_xy, g_0_xy_zz_xz, g_0_xy_zz_yy, g_0_xy_zz_yz, g_0_xy_zz_zz, g_z_0_0_0_z_xy_zz_xx, g_z_0_0_0_z_xy_zz_xy, g_z_0_0_0_z_xy_zz_xz, g_z_0_0_0_z_xy_zz_yy, g_z_0_0_0_z_xy_zz_yz, g_z_0_0_0_z_xy_zz_zz, g_zz_xy_zz_xx, g_zz_xy_zz_xy, g_zz_xy_zz_xz, g_zz_xy_zz_yy, g_zz_xy_zz_yz, g_zz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xy_zz_xx[i] = -g_0_xy_zz_xx[i] + 2.0 * g_zz_xy_zz_xx[i] * a_exp;

        g_z_0_0_0_z_xy_zz_xy[i] = -g_0_xy_zz_xy[i] + 2.0 * g_zz_xy_zz_xy[i] * a_exp;

        g_z_0_0_0_z_xy_zz_xz[i] = -g_0_xy_zz_xz[i] + 2.0 * g_zz_xy_zz_xz[i] * a_exp;

        g_z_0_0_0_z_xy_zz_yy[i] = -g_0_xy_zz_yy[i] + 2.0 * g_zz_xy_zz_yy[i] * a_exp;

        g_z_0_0_0_z_xy_zz_yz[i] = -g_0_xy_zz_yz[i] + 2.0 * g_zz_xy_zz_yz[i] * a_exp;

        g_z_0_0_0_z_xy_zz_zz[i] = -g_0_xy_zz_zz[i] + 2.0 * g_zz_xy_zz_zz[i] * a_exp;
    }
    // integrals block (1800-1806)

    #pragma omp simd aligned(g_0_xz_xx_xx, g_0_xz_xx_xy, g_0_xz_xx_xz, g_0_xz_xx_yy, g_0_xz_xx_yz, g_0_xz_xx_zz, g_z_0_0_0_z_xz_xx_xx, g_z_0_0_0_z_xz_xx_xy, g_z_0_0_0_z_xz_xx_xz, g_z_0_0_0_z_xz_xx_yy, g_z_0_0_0_z_xz_xx_yz, g_z_0_0_0_z_xz_xx_zz, g_zz_xz_xx_xx, g_zz_xz_xx_xy, g_zz_xz_xx_xz, g_zz_xz_xx_yy, g_zz_xz_xx_yz, g_zz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xz_xx_xx[i] = -g_0_xz_xx_xx[i] + 2.0 * g_zz_xz_xx_xx[i] * a_exp;

        g_z_0_0_0_z_xz_xx_xy[i] = -g_0_xz_xx_xy[i] + 2.0 * g_zz_xz_xx_xy[i] * a_exp;

        g_z_0_0_0_z_xz_xx_xz[i] = -g_0_xz_xx_xz[i] + 2.0 * g_zz_xz_xx_xz[i] * a_exp;

        g_z_0_0_0_z_xz_xx_yy[i] = -g_0_xz_xx_yy[i] + 2.0 * g_zz_xz_xx_yy[i] * a_exp;

        g_z_0_0_0_z_xz_xx_yz[i] = -g_0_xz_xx_yz[i] + 2.0 * g_zz_xz_xx_yz[i] * a_exp;

        g_z_0_0_0_z_xz_xx_zz[i] = -g_0_xz_xx_zz[i] + 2.0 * g_zz_xz_xx_zz[i] * a_exp;
    }
    // integrals block (1806-1812)

    #pragma omp simd aligned(g_0_xz_xy_xx, g_0_xz_xy_xy, g_0_xz_xy_xz, g_0_xz_xy_yy, g_0_xz_xy_yz, g_0_xz_xy_zz, g_z_0_0_0_z_xz_xy_xx, g_z_0_0_0_z_xz_xy_xy, g_z_0_0_0_z_xz_xy_xz, g_z_0_0_0_z_xz_xy_yy, g_z_0_0_0_z_xz_xy_yz, g_z_0_0_0_z_xz_xy_zz, g_zz_xz_xy_xx, g_zz_xz_xy_xy, g_zz_xz_xy_xz, g_zz_xz_xy_yy, g_zz_xz_xy_yz, g_zz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xz_xy_xx[i] = -g_0_xz_xy_xx[i] + 2.0 * g_zz_xz_xy_xx[i] * a_exp;

        g_z_0_0_0_z_xz_xy_xy[i] = -g_0_xz_xy_xy[i] + 2.0 * g_zz_xz_xy_xy[i] * a_exp;

        g_z_0_0_0_z_xz_xy_xz[i] = -g_0_xz_xy_xz[i] + 2.0 * g_zz_xz_xy_xz[i] * a_exp;

        g_z_0_0_0_z_xz_xy_yy[i] = -g_0_xz_xy_yy[i] + 2.0 * g_zz_xz_xy_yy[i] * a_exp;

        g_z_0_0_0_z_xz_xy_yz[i] = -g_0_xz_xy_yz[i] + 2.0 * g_zz_xz_xy_yz[i] * a_exp;

        g_z_0_0_0_z_xz_xy_zz[i] = -g_0_xz_xy_zz[i] + 2.0 * g_zz_xz_xy_zz[i] * a_exp;
    }
    // integrals block (1812-1818)

    #pragma omp simd aligned(g_0_xz_xz_xx, g_0_xz_xz_xy, g_0_xz_xz_xz, g_0_xz_xz_yy, g_0_xz_xz_yz, g_0_xz_xz_zz, g_z_0_0_0_z_xz_xz_xx, g_z_0_0_0_z_xz_xz_xy, g_z_0_0_0_z_xz_xz_xz, g_z_0_0_0_z_xz_xz_yy, g_z_0_0_0_z_xz_xz_yz, g_z_0_0_0_z_xz_xz_zz, g_zz_xz_xz_xx, g_zz_xz_xz_xy, g_zz_xz_xz_xz, g_zz_xz_xz_yy, g_zz_xz_xz_yz, g_zz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xz_xz_xx[i] = -g_0_xz_xz_xx[i] + 2.0 * g_zz_xz_xz_xx[i] * a_exp;

        g_z_0_0_0_z_xz_xz_xy[i] = -g_0_xz_xz_xy[i] + 2.0 * g_zz_xz_xz_xy[i] * a_exp;

        g_z_0_0_0_z_xz_xz_xz[i] = -g_0_xz_xz_xz[i] + 2.0 * g_zz_xz_xz_xz[i] * a_exp;

        g_z_0_0_0_z_xz_xz_yy[i] = -g_0_xz_xz_yy[i] + 2.0 * g_zz_xz_xz_yy[i] * a_exp;

        g_z_0_0_0_z_xz_xz_yz[i] = -g_0_xz_xz_yz[i] + 2.0 * g_zz_xz_xz_yz[i] * a_exp;

        g_z_0_0_0_z_xz_xz_zz[i] = -g_0_xz_xz_zz[i] + 2.0 * g_zz_xz_xz_zz[i] * a_exp;
    }
    // integrals block (1818-1824)

    #pragma omp simd aligned(g_0_xz_yy_xx, g_0_xz_yy_xy, g_0_xz_yy_xz, g_0_xz_yy_yy, g_0_xz_yy_yz, g_0_xz_yy_zz, g_z_0_0_0_z_xz_yy_xx, g_z_0_0_0_z_xz_yy_xy, g_z_0_0_0_z_xz_yy_xz, g_z_0_0_0_z_xz_yy_yy, g_z_0_0_0_z_xz_yy_yz, g_z_0_0_0_z_xz_yy_zz, g_zz_xz_yy_xx, g_zz_xz_yy_xy, g_zz_xz_yy_xz, g_zz_xz_yy_yy, g_zz_xz_yy_yz, g_zz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xz_yy_xx[i] = -g_0_xz_yy_xx[i] + 2.0 * g_zz_xz_yy_xx[i] * a_exp;

        g_z_0_0_0_z_xz_yy_xy[i] = -g_0_xz_yy_xy[i] + 2.0 * g_zz_xz_yy_xy[i] * a_exp;

        g_z_0_0_0_z_xz_yy_xz[i] = -g_0_xz_yy_xz[i] + 2.0 * g_zz_xz_yy_xz[i] * a_exp;

        g_z_0_0_0_z_xz_yy_yy[i] = -g_0_xz_yy_yy[i] + 2.0 * g_zz_xz_yy_yy[i] * a_exp;

        g_z_0_0_0_z_xz_yy_yz[i] = -g_0_xz_yy_yz[i] + 2.0 * g_zz_xz_yy_yz[i] * a_exp;

        g_z_0_0_0_z_xz_yy_zz[i] = -g_0_xz_yy_zz[i] + 2.0 * g_zz_xz_yy_zz[i] * a_exp;
    }
    // integrals block (1824-1830)

    #pragma omp simd aligned(g_0_xz_yz_xx, g_0_xz_yz_xy, g_0_xz_yz_xz, g_0_xz_yz_yy, g_0_xz_yz_yz, g_0_xz_yz_zz, g_z_0_0_0_z_xz_yz_xx, g_z_0_0_0_z_xz_yz_xy, g_z_0_0_0_z_xz_yz_xz, g_z_0_0_0_z_xz_yz_yy, g_z_0_0_0_z_xz_yz_yz, g_z_0_0_0_z_xz_yz_zz, g_zz_xz_yz_xx, g_zz_xz_yz_xy, g_zz_xz_yz_xz, g_zz_xz_yz_yy, g_zz_xz_yz_yz, g_zz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xz_yz_xx[i] = -g_0_xz_yz_xx[i] + 2.0 * g_zz_xz_yz_xx[i] * a_exp;

        g_z_0_0_0_z_xz_yz_xy[i] = -g_0_xz_yz_xy[i] + 2.0 * g_zz_xz_yz_xy[i] * a_exp;

        g_z_0_0_0_z_xz_yz_xz[i] = -g_0_xz_yz_xz[i] + 2.0 * g_zz_xz_yz_xz[i] * a_exp;

        g_z_0_0_0_z_xz_yz_yy[i] = -g_0_xz_yz_yy[i] + 2.0 * g_zz_xz_yz_yy[i] * a_exp;

        g_z_0_0_0_z_xz_yz_yz[i] = -g_0_xz_yz_yz[i] + 2.0 * g_zz_xz_yz_yz[i] * a_exp;

        g_z_0_0_0_z_xz_yz_zz[i] = -g_0_xz_yz_zz[i] + 2.0 * g_zz_xz_yz_zz[i] * a_exp;
    }
    // integrals block (1830-1836)

    #pragma omp simd aligned(g_0_xz_zz_xx, g_0_xz_zz_xy, g_0_xz_zz_xz, g_0_xz_zz_yy, g_0_xz_zz_yz, g_0_xz_zz_zz, g_z_0_0_0_z_xz_zz_xx, g_z_0_0_0_z_xz_zz_xy, g_z_0_0_0_z_xz_zz_xz, g_z_0_0_0_z_xz_zz_yy, g_z_0_0_0_z_xz_zz_yz, g_z_0_0_0_z_xz_zz_zz, g_zz_xz_zz_xx, g_zz_xz_zz_xy, g_zz_xz_zz_xz, g_zz_xz_zz_yy, g_zz_xz_zz_yz, g_zz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xz_zz_xx[i] = -g_0_xz_zz_xx[i] + 2.0 * g_zz_xz_zz_xx[i] * a_exp;

        g_z_0_0_0_z_xz_zz_xy[i] = -g_0_xz_zz_xy[i] + 2.0 * g_zz_xz_zz_xy[i] * a_exp;

        g_z_0_0_0_z_xz_zz_xz[i] = -g_0_xz_zz_xz[i] + 2.0 * g_zz_xz_zz_xz[i] * a_exp;

        g_z_0_0_0_z_xz_zz_yy[i] = -g_0_xz_zz_yy[i] + 2.0 * g_zz_xz_zz_yy[i] * a_exp;

        g_z_0_0_0_z_xz_zz_yz[i] = -g_0_xz_zz_yz[i] + 2.0 * g_zz_xz_zz_yz[i] * a_exp;

        g_z_0_0_0_z_xz_zz_zz[i] = -g_0_xz_zz_zz[i] + 2.0 * g_zz_xz_zz_zz[i] * a_exp;
    }
    // integrals block (1836-1842)

    #pragma omp simd aligned(g_0_yy_xx_xx, g_0_yy_xx_xy, g_0_yy_xx_xz, g_0_yy_xx_yy, g_0_yy_xx_yz, g_0_yy_xx_zz, g_z_0_0_0_z_yy_xx_xx, g_z_0_0_0_z_yy_xx_xy, g_z_0_0_0_z_yy_xx_xz, g_z_0_0_0_z_yy_xx_yy, g_z_0_0_0_z_yy_xx_yz, g_z_0_0_0_z_yy_xx_zz, g_zz_yy_xx_xx, g_zz_yy_xx_xy, g_zz_yy_xx_xz, g_zz_yy_xx_yy, g_zz_yy_xx_yz, g_zz_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yy_xx_xx[i] = -g_0_yy_xx_xx[i] + 2.0 * g_zz_yy_xx_xx[i] * a_exp;

        g_z_0_0_0_z_yy_xx_xy[i] = -g_0_yy_xx_xy[i] + 2.0 * g_zz_yy_xx_xy[i] * a_exp;

        g_z_0_0_0_z_yy_xx_xz[i] = -g_0_yy_xx_xz[i] + 2.0 * g_zz_yy_xx_xz[i] * a_exp;

        g_z_0_0_0_z_yy_xx_yy[i] = -g_0_yy_xx_yy[i] + 2.0 * g_zz_yy_xx_yy[i] * a_exp;

        g_z_0_0_0_z_yy_xx_yz[i] = -g_0_yy_xx_yz[i] + 2.0 * g_zz_yy_xx_yz[i] * a_exp;

        g_z_0_0_0_z_yy_xx_zz[i] = -g_0_yy_xx_zz[i] + 2.0 * g_zz_yy_xx_zz[i] * a_exp;
    }
    // integrals block (1842-1848)

    #pragma omp simd aligned(g_0_yy_xy_xx, g_0_yy_xy_xy, g_0_yy_xy_xz, g_0_yy_xy_yy, g_0_yy_xy_yz, g_0_yy_xy_zz, g_z_0_0_0_z_yy_xy_xx, g_z_0_0_0_z_yy_xy_xy, g_z_0_0_0_z_yy_xy_xz, g_z_0_0_0_z_yy_xy_yy, g_z_0_0_0_z_yy_xy_yz, g_z_0_0_0_z_yy_xy_zz, g_zz_yy_xy_xx, g_zz_yy_xy_xy, g_zz_yy_xy_xz, g_zz_yy_xy_yy, g_zz_yy_xy_yz, g_zz_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yy_xy_xx[i] = -g_0_yy_xy_xx[i] + 2.0 * g_zz_yy_xy_xx[i] * a_exp;

        g_z_0_0_0_z_yy_xy_xy[i] = -g_0_yy_xy_xy[i] + 2.0 * g_zz_yy_xy_xy[i] * a_exp;

        g_z_0_0_0_z_yy_xy_xz[i] = -g_0_yy_xy_xz[i] + 2.0 * g_zz_yy_xy_xz[i] * a_exp;

        g_z_0_0_0_z_yy_xy_yy[i] = -g_0_yy_xy_yy[i] + 2.0 * g_zz_yy_xy_yy[i] * a_exp;

        g_z_0_0_0_z_yy_xy_yz[i] = -g_0_yy_xy_yz[i] + 2.0 * g_zz_yy_xy_yz[i] * a_exp;

        g_z_0_0_0_z_yy_xy_zz[i] = -g_0_yy_xy_zz[i] + 2.0 * g_zz_yy_xy_zz[i] * a_exp;
    }
    // integrals block (1848-1854)

    #pragma omp simd aligned(g_0_yy_xz_xx, g_0_yy_xz_xy, g_0_yy_xz_xz, g_0_yy_xz_yy, g_0_yy_xz_yz, g_0_yy_xz_zz, g_z_0_0_0_z_yy_xz_xx, g_z_0_0_0_z_yy_xz_xy, g_z_0_0_0_z_yy_xz_xz, g_z_0_0_0_z_yy_xz_yy, g_z_0_0_0_z_yy_xz_yz, g_z_0_0_0_z_yy_xz_zz, g_zz_yy_xz_xx, g_zz_yy_xz_xy, g_zz_yy_xz_xz, g_zz_yy_xz_yy, g_zz_yy_xz_yz, g_zz_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yy_xz_xx[i] = -g_0_yy_xz_xx[i] + 2.0 * g_zz_yy_xz_xx[i] * a_exp;

        g_z_0_0_0_z_yy_xz_xy[i] = -g_0_yy_xz_xy[i] + 2.0 * g_zz_yy_xz_xy[i] * a_exp;

        g_z_0_0_0_z_yy_xz_xz[i] = -g_0_yy_xz_xz[i] + 2.0 * g_zz_yy_xz_xz[i] * a_exp;

        g_z_0_0_0_z_yy_xz_yy[i] = -g_0_yy_xz_yy[i] + 2.0 * g_zz_yy_xz_yy[i] * a_exp;

        g_z_0_0_0_z_yy_xz_yz[i] = -g_0_yy_xz_yz[i] + 2.0 * g_zz_yy_xz_yz[i] * a_exp;

        g_z_0_0_0_z_yy_xz_zz[i] = -g_0_yy_xz_zz[i] + 2.0 * g_zz_yy_xz_zz[i] * a_exp;
    }
    // integrals block (1854-1860)

    #pragma omp simd aligned(g_0_yy_yy_xx, g_0_yy_yy_xy, g_0_yy_yy_xz, g_0_yy_yy_yy, g_0_yy_yy_yz, g_0_yy_yy_zz, g_z_0_0_0_z_yy_yy_xx, g_z_0_0_0_z_yy_yy_xy, g_z_0_0_0_z_yy_yy_xz, g_z_0_0_0_z_yy_yy_yy, g_z_0_0_0_z_yy_yy_yz, g_z_0_0_0_z_yy_yy_zz, g_zz_yy_yy_xx, g_zz_yy_yy_xy, g_zz_yy_yy_xz, g_zz_yy_yy_yy, g_zz_yy_yy_yz, g_zz_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yy_yy_xx[i] = -g_0_yy_yy_xx[i] + 2.0 * g_zz_yy_yy_xx[i] * a_exp;

        g_z_0_0_0_z_yy_yy_xy[i] = -g_0_yy_yy_xy[i] + 2.0 * g_zz_yy_yy_xy[i] * a_exp;

        g_z_0_0_0_z_yy_yy_xz[i] = -g_0_yy_yy_xz[i] + 2.0 * g_zz_yy_yy_xz[i] * a_exp;

        g_z_0_0_0_z_yy_yy_yy[i] = -g_0_yy_yy_yy[i] + 2.0 * g_zz_yy_yy_yy[i] * a_exp;

        g_z_0_0_0_z_yy_yy_yz[i] = -g_0_yy_yy_yz[i] + 2.0 * g_zz_yy_yy_yz[i] * a_exp;

        g_z_0_0_0_z_yy_yy_zz[i] = -g_0_yy_yy_zz[i] + 2.0 * g_zz_yy_yy_zz[i] * a_exp;
    }
    // integrals block (1860-1866)

    #pragma omp simd aligned(g_0_yy_yz_xx, g_0_yy_yz_xy, g_0_yy_yz_xz, g_0_yy_yz_yy, g_0_yy_yz_yz, g_0_yy_yz_zz, g_z_0_0_0_z_yy_yz_xx, g_z_0_0_0_z_yy_yz_xy, g_z_0_0_0_z_yy_yz_xz, g_z_0_0_0_z_yy_yz_yy, g_z_0_0_0_z_yy_yz_yz, g_z_0_0_0_z_yy_yz_zz, g_zz_yy_yz_xx, g_zz_yy_yz_xy, g_zz_yy_yz_xz, g_zz_yy_yz_yy, g_zz_yy_yz_yz, g_zz_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yy_yz_xx[i] = -g_0_yy_yz_xx[i] + 2.0 * g_zz_yy_yz_xx[i] * a_exp;

        g_z_0_0_0_z_yy_yz_xy[i] = -g_0_yy_yz_xy[i] + 2.0 * g_zz_yy_yz_xy[i] * a_exp;

        g_z_0_0_0_z_yy_yz_xz[i] = -g_0_yy_yz_xz[i] + 2.0 * g_zz_yy_yz_xz[i] * a_exp;

        g_z_0_0_0_z_yy_yz_yy[i] = -g_0_yy_yz_yy[i] + 2.0 * g_zz_yy_yz_yy[i] * a_exp;

        g_z_0_0_0_z_yy_yz_yz[i] = -g_0_yy_yz_yz[i] + 2.0 * g_zz_yy_yz_yz[i] * a_exp;

        g_z_0_0_0_z_yy_yz_zz[i] = -g_0_yy_yz_zz[i] + 2.0 * g_zz_yy_yz_zz[i] * a_exp;
    }
    // integrals block (1866-1872)

    #pragma omp simd aligned(g_0_yy_zz_xx, g_0_yy_zz_xy, g_0_yy_zz_xz, g_0_yy_zz_yy, g_0_yy_zz_yz, g_0_yy_zz_zz, g_z_0_0_0_z_yy_zz_xx, g_z_0_0_0_z_yy_zz_xy, g_z_0_0_0_z_yy_zz_xz, g_z_0_0_0_z_yy_zz_yy, g_z_0_0_0_z_yy_zz_yz, g_z_0_0_0_z_yy_zz_zz, g_zz_yy_zz_xx, g_zz_yy_zz_xy, g_zz_yy_zz_xz, g_zz_yy_zz_yy, g_zz_yy_zz_yz, g_zz_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yy_zz_xx[i] = -g_0_yy_zz_xx[i] + 2.0 * g_zz_yy_zz_xx[i] * a_exp;

        g_z_0_0_0_z_yy_zz_xy[i] = -g_0_yy_zz_xy[i] + 2.0 * g_zz_yy_zz_xy[i] * a_exp;

        g_z_0_0_0_z_yy_zz_xz[i] = -g_0_yy_zz_xz[i] + 2.0 * g_zz_yy_zz_xz[i] * a_exp;

        g_z_0_0_0_z_yy_zz_yy[i] = -g_0_yy_zz_yy[i] + 2.0 * g_zz_yy_zz_yy[i] * a_exp;

        g_z_0_0_0_z_yy_zz_yz[i] = -g_0_yy_zz_yz[i] + 2.0 * g_zz_yy_zz_yz[i] * a_exp;

        g_z_0_0_0_z_yy_zz_zz[i] = -g_0_yy_zz_zz[i] + 2.0 * g_zz_yy_zz_zz[i] * a_exp;
    }
    // integrals block (1872-1878)

    #pragma omp simd aligned(g_0_yz_xx_xx, g_0_yz_xx_xy, g_0_yz_xx_xz, g_0_yz_xx_yy, g_0_yz_xx_yz, g_0_yz_xx_zz, g_z_0_0_0_z_yz_xx_xx, g_z_0_0_0_z_yz_xx_xy, g_z_0_0_0_z_yz_xx_xz, g_z_0_0_0_z_yz_xx_yy, g_z_0_0_0_z_yz_xx_yz, g_z_0_0_0_z_yz_xx_zz, g_zz_yz_xx_xx, g_zz_yz_xx_xy, g_zz_yz_xx_xz, g_zz_yz_xx_yy, g_zz_yz_xx_yz, g_zz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yz_xx_xx[i] = -g_0_yz_xx_xx[i] + 2.0 * g_zz_yz_xx_xx[i] * a_exp;

        g_z_0_0_0_z_yz_xx_xy[i] = -g_0_yz_xx_xy[i] + 2.0 * g_zz_yz_xx_xy[i] * a_exp;

        g_z_0_0_0_z_yz_xx_xz[i] = -g_0_yz_xx_xz[i] + 2.0 * g_zz_yz_xx_xz[i] * a_exp;

        g_z_0_0_0_z_yz_xx_yy[i] = -g_0_yz_xx_yy[i] + 2.0 * g_zz_yz_xx_yy[i] * a_exp;

        g_z_0_0_0_z_yz_xx_yz[i] = -g_0_yz_xx_yz[i] + 2.0 * g_zz_yz_xx_yz[i] * a_exp;

        g_z_0_0_0_z_yz_xx_zz[i] = -g_0_yz_xx_zz[i] + 2.0 * g_zz_yz_xx_zz[i] * a_exp;
    }
    // integrals block (1878-1884)

    #pragma omp simd aligned(g_0_yz_xy_xx, g_0_yz_xy_xy, g_0_yz_xy_xz, g_0_yz_xy_yy, g_0_yz_xy_yz, g_0_yz_xy_zz, g_z_0_0_0_z_yz_xy_xx, g_z_0_0_0_z_yz_xy_xy, g_z_0_0_0_z_yz_xy_xz, g_z_0_0_0_z_yz_xy_yy, g_z_0_0_0_z_yz_xy_yz, g_z_0_0_0_z_yz_xy_zz, g_zz_yz_xy_xx, g_zz_yz_xy_xy, g_zz_yz_xy_xz, g_zz_yz_xy_yy, g_zz_yz_xy_yz, g_zz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yz_xy_xx[i] = -g_0_yz_xy_xx[i] + 2.0 * g_zz_yz_xy_xx[i] * a_exp;

        g_z_0_0_0_z_yz_xy_xy[i] = -g_0_yz_xy_xy[i] + 2.0 * g_zz_yz_xy_xy[i] * a_exp;

        g_z_0_0_0_z_yz_xy_xz[i] = -g_0_yz_xy_xz[i] + 2.0 * g_zz_yz_xy_xz[i] * a_exp;

        g_z_0_0_0_z_yz_xy_yy[i] = -g_0_yz_xy_yy[i] + 2.0 * g_zz_yz_xy_yy[i] * a_exp;

        g_z_0_0_0_z_yz_xy_yz[i] = -g_0_yz_xy_yz[i] + 2.0 * g_zz_yz_xy_yz[i] * a_exp;

        g_z_0_0_0_z_yz_xy_zz[i] = -g_0_yz_xy_zz[i] + 2.0 * g_zz_yz_xy_zz[i] * a_exp;
    }
    // integrals block (1884-1890)

    #pragma omp simd aligned(g_0_yz_xz_xx, g_0_yz_xz_xy, g_0_yz_xz_xz, g_0_yz_xz_yy, g_0_yz_xz_yz, g_0_yz_xz_zz, g_z_0_0_0_z_yz_xz_xx, g_z_0_0_0_z_yz_xz_xy, g_z_0_0_0_z_yz_xz_xz, g_z_0_0_0_z_yz_xz_yy, g_z_0_0_0_z_yz_xz_yz, g_z_0_0_0_z_yz_xz_zz, g_zz_yz_xz_xx, g_zz_yz_xz_xy, g_zz_yz_xz_xz, g_zz_yz_xz_yy, g_zz_yz_xz_yz, g_zz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yz_xz_xx[i] = -g_0_yz_xz_xx[i] + 2.0 * g_zz_yz_xz_xx[i] * a_exp;

        g_z_0_0_0_z_yz_xz_xy[i] = -g_0_yz_xz_xy[i] + 2.0 * g_zz_yz_xz_xy[i] * a_exp;

        g_z_0_0_0_z_yz_xz_xz[i] = -g_0_yz_xz_xz[i] + 2.0 * g_zz_yz_xz_xz[i] * a_exp;

        g_z_0_0_0_z_yz_xz_yy[i] = -g_0_yz_xz_yy[i] + 2.0 * g_zz_yz_xz_yy[i] * a_exp;

        g_z_0_0_0_z_yz_xz_yz[i] = -g_0_yz_xz_yz[i] + 2.0 * g_zz_yz_xz_yz[i] * a_exp;

        g_z_0_0_0_z_yz_xz_zz[i] = -g_0_yz_xz_zz[i] + 2.0 * g_zz_yz_xz_zz[i] * a_exp;
    }
    // integrals block (1890-1896)

    #pragma omp simd aligned(g_0_yz_yy_xx, g_0_yz_yy_xy, g_0_yz_yy_xz, g_0_yz_yy_yy, g_0_yz_yy_yz, g_0_yz_yy_zz, g_z_0_0_0_z_yz_yy_xx, g_z_0_0_0_z_yz_yy_xy, g_z_0_0_0_z_yz_yy_xz, g_z_0_0_0_z_yz_yy_yy, g_z_0_0_0_z_yz_yy_yz, g_z_0_0_0_z_yz_yy_zz, g_zz_yz_yy_xx, g_zz_yz_yy_xy, g_zz_yz_yy_xz, g_zz_yz_yy_yy, g_zz_yz_yy_yz, g_zz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yz_yy_xx[i] = -g_0_yz_yy_xx[i] + 2.0 * g_zz_yz_yy_xx[i] * a_exp;

        g_z_0_0_0_z_yz_yy_xy[i] = -g_0_yz_yy_xy[i] + 2.0 * g_zz_yz_yy_xy[i] * a_exp;

        g_z_0_0_0_z_yz_yy_xz[i] = -g_0_yz_yy_xz[i] + 2.0 * g_zz_yz_yy_xz[i] * a_exp;

        g_z_0_0_0_z_yz_yy_yy[i] = -g_0_yz_yy_yy[i] + 2.0 * g_zz_yz_yy_yy[i] * a_exp;

        g_z_0_0_0_z_yz_yy_yz[i] = -g_0_yz_yy_yz[i] + 2.0 * g_zz_yz_yy_yz[i] * a_exp;

        g_z_0_0_0_z_yz_yy_zz[i] = -g_0_yz_yy_zz[i] + 2.0 * g_zz_yz_yy_zz[i] * a_exp;
    }
    // integrals block (1896-1902)

    #pragma omp simd aligned(g_0_yz_yz_xx, g_0_yz_yz_xy, g_0_yz_yz_xz, g_0_yz_yz_yy, g_0_yz_yz_yz, g_0_yz_yz_zz, g_z_0_0_0_z_yz_yz_xx, g_z_0_0_0_z_yz_yz_xy, g_z_0_0_0_z_yz_yz_xz, g_z_0_0_0_z_yz_yz_yy, g_z_0_0_0_z_yz_yz_yz, g_z_0_0_0_z_yz_yz_zz, g_zz_yz_yz_xx, g_zz_yz_yz_xy, g_zz_yz_yz_xz, g_zz_yz_yz_yy, g_zz_yz_yz_yz, g_zz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yz_yz_xx[i] = -g_0_yz_yz_xx[i] + 2.0 * g_zz_yz_yz_xx[i] * a_exp;

        g_z_0_0_0_z_yz_yz_xy[i] = -g_0_yz_yz_xy[i] + 2.0 * g_zz_yz_yz_xy[i] * a_exp;

        g_z_0_0_0_z_yz_yz_xz[i] = -g_0_yz_yz_xz[i] + 2.0 * g_zz_yz_yz_xz[i] * a_exp;

        g_z_0_0_0_z_yz_yz_yy[i] = -g_0_yz_yz_yy[i] + 2.0 * g_zz_yz_yz_yy[i] * a_exp;

        g_z_0_0_0_z_yz_yz_yz[i] = -g_0_yz_yz_yz[i] + 2.0 * g_zz_yz_yz_yz[i] * a_exp;

        g_z_0_0_0_z_yz_yz_zz[i] = -g_0_yz_yz_zz[i] + 2.0 * g_zz_yz_yz_zz[i] * a_exp;
    }
    // integrals block (1902-1908)

    #pragma omp simd aligned(g_0_yz_zz_xx, g_0_yz_zz_xy, g_0_yz_zz_xz, g_0_yz_zz_yy, g_0_yz_zz_yz, g_0_yz_zz_zz, g_z_0_0_0_z_yz_zz_xx, g_z_0_0_0_z_yz_zz_xy, g_z_0_0_0_z_yz_zz_xz, g_z_0_0_0_z_yz_zz_yy, g_z_0_0_0_z_yz_zz_yz, g_z_0_0_0_z_yz_zz_zz, g_zz_yz_zz_xx, g_zz_yz_zz_xy, g_zz_yz_zz_xz, g_zz_yz_zz_yy, g_zz_yz_zz_yz, g_zz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yz_zz_xx[i] = -g_0_yz_zz_xx[i] + 2.0 * g_zz_yz_zz_xx[i] * a_exp;

        g_z_0_0_0_z_yz_zz_xy[i] = -g_0_yz_zz_xy[i] + 2.0 * g_zz_yz_zz_xy[i] * a_exp;

        g_z_0_0_0_z_yz_zz_xz[i] = -g_0_yz_zz_xz[i] + 2.0 * g_zz_yz_zz_xz[i] * a_exp;

        g_z_0_0_0_z_yz_zz_yy[i] = -g_0_yz_zz_yy[i] + 2.0 * g_zz_yz_zz_yy[i] * a_exp;

        g_z_0_0_0_z_yz_zz_yz[i] = -g_0_yz_zz_yz[i] + 2.0 * g_zz_yz_zz_yz[i] * a_exp;

        g_z_0_0_0_z_yz_zz_zz[i] = -g_0_yz_zz_zz[i] + 2.0 * g_zz_yz_zz_zz[i] * a_exp;
    }
    // integrals block (1908-1914)

    #pragma omp simd aligned(g_0_zz_xx_xx, g_0_zz_xx_xy, g_0_zz_xx_xz, g_0_zz_xx_yy, g_0_zz_xx_yz, g_0_zz_xx_zz, g_z_0_0_0_z_zz_xx_xx, g_z_0_0_0_z_zz_xx_xy, g_z_0_0_0_z_zz_xx_xz, g_z_0_0_0_z_zz_xx_yy, g_z_0_0_0_z_zz_xx_yz, g_z_0_0_0_z_zz_xx_zz, g_zz_zz_xx_xx, g_zz_zz_xx_xy, g_zz_zz_xx_xz, g_zz_zz_xx_yy, g_zz_zz_xx_yz, g_zz_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_zz_xx_xx[i] = -g_0_zz_xx_xx[i] + 2.0 * g_zz_zz_xx_xx[i] * a_exp;

        g_z_0_0_0_z_zz_xx_xy[i] = -g_0_zz_xx_xy[i] + 2.0 * g_zz_zz_xx_xy[i] * a_exp;

        g_z_0_0_0_z_zz_xx_xz[i] = -g_0_zz_xx_xz[i] + 2.0 * g_zz_zz_xx_xz[i] * a_exp;

        g_z_0_0_0_z_zz_xx_yy[i] = -g_0_zz_xx_yy[i] + 2.0 * g_zz_zz_xx_yy[i] * a_exp;

        g_z_0_0_0_z_zz_xx_yz[i] = -g_0_zz_xx_yz[i] + 2.0 * g_zz_zz_xx_yz[i] * a_exp;

        g_z_0_0_0_z_zz_xx_zz[i] = -g_0_zz_xx_zz[i] + 2.0 * g_zz_zz_xx_zz[i] * a_exp;
    }
    // integrals block (1914-1920)

    #pragma omp simd aligned(g_0_zz_xy_xx, g_0_zz_xy_xy, g_0_zz_xy_xz, g_0_zz_xy_yy, g_0_zz_xy_yz, g_0_zz_xy_zz, g_z_0_0_0_z_zz_xy_xx, g_z_0_0_0_z_zz_xy_xy, g_z_0_0_0_z_zz_xy_xz, g_z_0_0_0_z_zz_xy_yy, g_z_0_0_0_z_zz_xy_yz, g_z_0_0_0_z_zz_xy_zz, g_zz_zz_xy_xx, g_zz_zz_xy_xy, g_zz_zz_xy_xz, g_zz_zz_xy_yy, g_zz_zz_xy_yz, g_zz_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_zz_xy_xx[i] = -g_0_zz_xy_xx[i] + 2.0 * g_zz_zz_xy_xx[i] * a_exp;

        g_z_0_0_0_z_zz_xy_xy[i] = -g_0_zz_xy_xy[i] + 2.0 * g_zz_zz_xy_xy[i] * a_exp;

        g_z_0_0_0_z_zz_xy_xz[i] = -g_0_zz_xy_xz[i] + 2.0 * g_zz_zz_xy_xz[i] * a_exp;

        g_z_0_0_0_z_zz_xy_yy[i] = -g_0_zz_xy_yy[i] + 2.0 * g_zz_zz_xy_yy[i] * a_exp;

        g_z_0_0_0_z_zz_xy_yz[i] = -g_0_zz_xy_yz[i] + 2.0 * g_zz_zz_xy_yz[i] * a_exp;

        g_z_0_0_0_z_zz_xy_zz[i] = -g_0_zz_xy_zz[i] + 2.0 * g_zz_zz_xy_zz[i] * a_exp;
    }
    // integrals block (1920-1926)

    #pragma omp simd aligned(g_0_zz_xz_xx, g_0_zz_xz_xy, g_0_zz_xz_xz, g_0_zz_xz_yy, g_0_zz_xz_yz, g_0_zz_xz_zz, g_z_0_0_0_z_zz_xz_xx, g_z_0_0_0_z_zz_xz_xy, g_z_0_0_0_z_zz_xz_xz, g_z_0_0_0_z_zz_xz_yy, g_z_0_0_0_z_zz_xz_yz, g_z_0_0_0_z_zz_xz_zz, g_zz_zz_xz_xx, g_zz_zz_xz_xy, g_zz_zz_xz_xz, g_zz_zz_xz_yy, g_zz_zz_xz_yz, g_zz_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_zz_xz_xx[i] = -g_0_zz_xz_xx[i] + 2.0 * g_zz_zz_xz_xx[i] * a_exp;

        g_z_0_0_0_z_zz_xz_xy[i] = -g_0_zz_xz_xy[i] + 2.0 * g_zz_zz_xz_xy[i] * a_exp;

        g_z_0_0_0_z_zz_xz_xz[i] = -g_0_zz_xz_xz[i] + 2.0 * g_zz_zz_xz_xz[i] * a_exp;

        g_z_0_0_0_z_zz_xz_yy[i] = -g_0_zz_xz_yy[i] + 2.0 * g_zz_zz_xz_yy[i] * a_exp;

        g_z_0_0_0_z_zz_xz_yz[i] = -g_0_zz_xz_yz[i] + 2.0 * g_zz_zz_xz_yz[i] * a_exp;

        g_z_0_0_0_z_zz_xz_zz[i] = -g_0_zz_xz_zz[i] + 2.0 * g_zz_zz_xz_zz[i] * a_exp;
    }
    // integrals block (1926-1932)

    #pragma omp simd aligned(g_0_zz_yy_xx, g_0_zz_yy_xy, g_0_zz_yy_xz, g_0_zz_yy_yy, g_0_zz_yy_yz, g_0_zz_yy_zz, g_z_0_0_0_z_zz_yy_xx, g_z_0_0_0_z_zz_yy_xy, g_z_0_0_0_z_zz_yy_xz, g_z_0_0_0_z_zz_yy_yy, g_z_0_0_0_z_zz_yy_yz, g_z_0_0_0_z_zz_yy_zz, g_zz_zz_yy_xx, g_zz_zz_yy_xy, g_zz_zz_yy_xz, g_zz_zz_yy_yy, g_zz_zz_yy_yz, g_zz_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_zz_yy_xx[i] = -g_0_zz_yy_xx[i] + 2.0 * g_zz_zz_yy_xx[i] * a_exp;

        g_z_0_0_0_z_zz_yy_xy[i] = -g_0_zz_yy_xy[i] + 2.0 * g_zz_zz_yy_xy[i] * a_exp;

        g_z_0_0_0_z_zz_yy_xz[i] = -g_0_zz_yy_xz[i] + 2.0 * g_zz_zz_yy_xz[i] * a_exp;

        g_z_0_0_0_z_zz_yy_yy[i] = -g_0_zz_yy_yy[i] + 2.0 * g_zz_zz_yy_yy[i] * a_exp;

        g_z_0_0_0_z_zz_yy_yz[i] = -g_0_zz_yy_yz[i] + 2.0 * g_zz_zz_yy_yz[i] * a_exp;

        g_z_0_0_0_z_zz_yy_zz[i] = -g_0_zz_yy_zz[i] + 2.0 * g_zz_zz_yy_zz[i] * a_exp;
    }
    // integrals block (1932-1938)

    #pragma omp simd aligned(g_0_zz_yz_xx, g_0_zz_yz_xy, g_0_zz_yz_xz, g_0_zz_yz_yy, g_0_zz_yz_yz, g_0_zz_yz_zz, g_z_0_0_0_z_zz_yz_xx, g_z_0_0_0_z_zz_yz_xy, g_z_0_0_0_z_zz_yz_xz, g_z_0_0_0_z_zz_yz_yy, g_z_0_0_0_z_zz_yz_yz, g_z_0_0_0_z_zz_yz_zz, g_zz_zz_yz_xx, g_zz_zz_yz_xy, g_zz_zz_yz_xz, g_zz_zz_yz_yy, g_zz_zz_yz_yz, g_zz_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_zz_yz_xx[i] = -g_0_zz_yz_xx[i] + 2.0 * g_zz_zz_yz_xx[i] * a_exp;

        g_z_0_0_0_z_zz_yz_xy[i] = -g_0_zz_yz_xy[i] + 2.0 * g_zz_zz_yz_xy[i] * a_exp;

        g_z_0_0_0_z_zz_yz_xz[i] = -g_0_zz_yz_xz[i] + 2.0 * g_zz_zz_yz_xz[i] * a_exp;

        g_z_0_0_0_z_zz_yz_yy[i] = -g_0_zz_yz_yy[i] + 2.0 * g_zz_zz_yz_yy[i] * a_exp;

        g_z_0_0_0_z_zz_yz_yz[i] = -g_0_zz_yz_yz[i] + 2.0 * g_zz_zz_yz_yz[i] * a_exp;

        g_z_0_0_0_z_zz_yz_zz[i] = -g_0_zz_yz_zz[i] + 2.0 * g_zz_zz_yz_zz[i] * a_exp;
    }
    // integrals block (1938-1944)

    #pragma omp simd aligned(g_0_zz_zz_xx, g_0_zz_zz_xy, g_0_zz_zz_xz, g_0_zz_zz_yy, g_0_zz_zz_yz, g_0_zz_zz_zz, g_z_0_0_0_z_zz_zz_xx, g_z_0_0_0_z_zz_zz_xy, g_z_0_0_0_z_zz_zz_xz, g_z_0_0_0_z_zz_zz_yy, g_z_0_0_0_z_zz_zz_yz, g_z_0_0_0_z_zz_zz_zz, g_zz_zz_zz_xx, g_zz_zz_zz_xy, g_zz_zz_zz_xz, g_zz_zz_zz_yy, g_zz_zz_zz_yz, g_zz_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_zz_zz_xx[i] = -g_0_zz_zz_xx[i] + 2.0 * g_zz_zz_zz_xx[i] * a_exp;

        g_z_0_0_0_z_zz_zz_xy[i] = -g_0_zz_zz_xy[i] + 2.0 * g_zz_zz_zz_xy[i] * a_exp;

        g_z_0_0_0_z_zz_zz_xz[i] = -g_0_zz_zz_xz[i] + 2.0 * g_zz_zz_zz_xz[i] * a_exp;

        g_z_0_0_0_z_zz_zz_yy[i] = -g_0_zz_zz_yy[i] + 2.0 * g_zz_zz_zz_yy[i] * a_exp;

        g_z_0_0_0_z_zz_zz_yz[i] = -g_0_zz_zz_yz[i] + 2.0 * g_zz_zz_zz_yz[i] * a_exp;

        g_z_0_0_0_z_zz_zz_zz[i] = -g_0_zz_zz_zz[i] + 2.0 * g_zz_zz_zz_zz[i] * a_exp;
    }
}

} // t4c_geom namespace

