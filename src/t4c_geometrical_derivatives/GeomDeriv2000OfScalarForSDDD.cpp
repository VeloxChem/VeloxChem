#include "GeomDeriv2000OfScalarForSDDD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_sddd_0(CSimdArray<double>& buffer_2000_sddd,
                     const CSimdArray<double>& buffer_sddd,
                     const CSimdArray<double>& buffer_dddd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_sddd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_2000_sddd

    auto g_xx_0_0_0_0_xx_xx_xx = buffer_2000_sddd[0];

    auto g_xx_0_0_0_0_xx_xx_xy = buffer_2000_sddd[1];

    auto g_xx_0_0_0_0_xx_xx_xz = buffer_2000_sddd[2];

    auto g_xx_0_0_0_0_xx_xx_yy = buffer_2000_sddd[3];

    auto g_xx_0_0_0_0_xx_xx_yz = buffer_2000_sddd[4];

    auto g_xx_0_0_0_0_xx_xx_zz = buffer_2000_sddd[5];

    auto g_xx_0_0_0_0_xx_xy_xx = buffer_2000_sddd[6];

    auto g_xx_0_0_0_0_xx_xy_xy = buffer_2000_sddd[7];

    auto g_xx_0_0_0_0_xx_xy_xz = buffer_2000_sddd[8];

    auto g_xx_0_0_0_0_xx_xy_yy = buffer_2000_sddd[9];

    auto g_xx_0_0_0_0_xx_xy_yz = buffer_2000_sddd[10];

    auto g_xx_0_0_0_0_xx_xy_zz = buffer_2000_sddd[11];

    auto g_xx_0_0_0_0_xx_xz_xx = buffer_2000_sddd[12];

    auto g_xx_0_0_0_0_xx_xz_xy = buffer_2000_sddd[13];

    auto g_xx_0_0_0_0_xx_xz_xz = buffer_2000_sddd[14];

    auto g_xx_0_0_0_0_xx_xz_yy = buffer_2000_sddd[15];

    auto g_xx_0_0_0_0_xx_xz_yz = buffer_2000_sddd[16];

    auto g_xx_0_0_0_0_xx_xz_zz = buffer_2000_sddd[17];

    auto g_xx_0_0_0_0_xx_yy_xx = buffer_2000_sddd[18];

    auto g_xx_0_0_0_0_xx_yy_xy = buffer_2000_sddd[19];

    auto g_xx_0_0_0_0_xx_yy_xz = buffer_2000_sddd[20];

    auto g_xx_0_0_0_0_xx_yy_yy = buffer_2000_sddd[21];

    auto g_xx_0_0_0_0_xx_yy_yz = buffer_2000_sddd[22];

    auto g_xx_0_0_0_0_xx_yy_zz = buffer_2000_sddd[23];

    auto g_xx_0_0_0_0_xx_yz_xx = buffer_2000_sddd[24];

    auto g_xx_0_0_0_0_xx_yz_xy = buffer_2000_sddd[25];

    auto g_xx_0_0_0_0_xx_yz_xz = buffer_2000_sddd[26];

    auto g_xx_0_0_0_0_xx_yz_yy = buffer_2000_sddd[27];

    auto g_xx_0_0_0_0_xx_yz_yz = buffer_2000_sddd[28];

    auto g_xx_0_0_0_0_xx_yz_zz = buffer_2000_sddd[29];

    auto g_xx_0_0_0_0_xx_zz_xx = buffer_2000_sddd[30];

    auto g_xx_0_0_0_0_xx_zz_xy = buffer_2000_sddd[31];

    auto g_xx_0_0_0_0_xx_zz_xz = buffer_2000_sddd[32];

    auto g_xx_0_0_0_0_xx_zz_yy = buffer_2000_sddd[33];

    auto g_xx_0_0_0_0_xx_zz_yz = buffer_2000_sddd[34];

    auto g_xx_0_0_0_0_xx_zz_zz = buffer_2000_sddd[35];

    auto g_xx_0_0_0_0_xy_xx_xx = buffer_2000_sddd[36];

    auto g_xx_0_0_0_0_xy_xx_xy = buffer_2000_sddd[37];

    auto g_xx_0_0_0_0_xy_xx_xz = buffer_2000_sddd[38];

    auto g_xx_0_0_0_0_xy_xx_yy = buffer_2000_sddd[39];

    auto g_xx_0_0_0_0_xy_xx_yz = buffer_2000_sddd[40];

    auto g_xx_0_0_0_0_xy_xx_zz = buffer_2000_sddd[41];

    auto g_xx_0_0_0_0_xy_xy_xx = buffer_2000_sddd[42];

    auto g_xx_0_0_0_0_xy_xy_xy = buffer_2000_sddd[43];

    auto g_xx_0_0_0_0_xy_xy_xz = buffer_2000_sddd[44];

    auto g_xx_0_0_0_0_xy_xy_yy = buffer_2000_sddd[45];

    auto g_xx_0_0_0_0_xy_xy_yz = buffer_2000_sddd[46];

    auto g_xx_0_0_0_0_xy_xy_zz = buffer_2000_sddd[47];

    auto g_xx_0_0_0_0_xy_xz_xx = buffer_2000_sddd[48];

    auto g_xx_0_0_0_0_xy_xz_xy = buffer_2000_sddd[49];

    auto g_xx_0_0_0_0_xy_xz_xz = buffer_2000_sddd[50];

    auto g_xx_0_0_0_0_xy_xz_yy = buffer_2000_sddd[51];

    auto g_xx_0_0_0_0_xy_xz_yz = buffer_2000_sddd[52];

    auto g_xx_0_0_0_0_xy_xz_zz = buffer_2000_sddd[53];

    auto g_xx_0_0_0_0_xy_yy_xx = buffer_2000_sddd[54];

    auto g_xx_0_0_0_0_xy_yy_xy = buffer_2000_sddd[55];

    auto g_xx_0_0_0_0_xy_yy_xz = buffer_2000_sddd[56];

    auto g_xx_0_0_0_0_xy_yy_yy = buffer_2000_sddd[57];

    auto g_xx_0_0_0_0_xy_yy_yz = buffer_2000_sddd[58];

    auto g_xx_0_0_0_0_xy_yy_zz = buffer_2000_sddd[59];

    auto g_xx_0_0_0_0_xy_yz_xx = buffer_2000_sddd[60];

    auto g_xx_0_0_0_0_xy_yz_xy = buffer_2000_sddd[61];

    auto g_xx_0_0_0_0_xy_yz_xz = buffer_2000_sddd[62];

    auto g_xx_0_0_0_0_xy_yz_yy = buffer_2000_sddd[63];

    auto g_xx_0_0_0_0_xy_yz_yz = buffer_2000_sddd[64];

    auto g_xx_0_0_0_0_xy_yz_zz = buffer_2000_sddd[65];

    auto g_xx_0_0_0_0_xy_zz_xx = buffer_2000_sddd[66];

    auto g_xx_0_0_0_0_xy_zz_xy = buffer_2000_sddd[67];

    auto g_xx_0_0_0_0_xy_zz_xz = buffer_2000_sddd[68];

    auto g_xx_0_0_0_0_xy_zz_yy = buffer_2000_sddd[69];

    auto g_xx_0_0_0_0_xy_zz_yz = buffer_2000_sddd[70];

    auto g_xx_0_0_0_0_xy_zz_zz = buffer_2000_sddd[71];

    auto g_xx_0_0_0_0_xz_xx_xx = buffer_2000_sddd[72];

    auto g_xx_0_0_0_0_xz_xx_xy = buffer_2000_sddd[73];

    auto g_xx_0_0_0_0_xz_xx_xz = buffer_2000_sddd[74];

    auto g_xx_0_0_0_0_xz_xx_yy = buffer_2000_sddd[75];

    auto g_xx_0_0_0_0_xz_xx_yz = buffer_2000_sddd[76];

    auto g_xx_0_0_0_0_xz_xx_zz = buffer_2000_sddd[77];

    auto g_xx_0_0_0_0_xz_xy_xx = buffer_2000_sddd[78];

    auto g_xx_0_0_0_0_xz_xy_xy = buffer_2000_sddd[79];

    auto g_xx_0_0_0_0_xz_xy_xz = buffer_2000_sddd[80];

    auto g_xx_0_0_0_0_xz_xy_yy = buffer_2000_sddd[81];

    auto g_xx_0_0_0_0_xz_xy_yz = buffer_2000_sddd[82];

    auto g_xx_0_0_0_0_xz_xy_zz = buffer_2000_sddd[83];

    auto g_xx_0_0_0_0_xz_xz_xx = buffer_2000_sddd[84];

    auto g_xx_0_0_0_0_xz_xz_xy = buffer_2000_sddd[85];

    auto g_xx_0_0_0_0_xz_xz_xz = buffer_2000_sddd[86];

    auto g_xx_0_0_0_0_xz_xz_yy = buffer_2000_sddd[87];

    auto g_xx_0_0_0_0_xz_xz_yz = buffer_2000_sddd[88];

    auto g_xx_0_0_0_0_xz_xz_zz = buffer_2000_sddd[89];

    auto g_xx_0_0_0_0_xz_yy_xx = buffer_2000_sddd[90];

    auto g_xx_0_0_0_0_xz_yy_xy = buffer_2000_sddd[91];

    auto g_xx_0_0_0_0_xz_yy_xz = buffer_2000_sddd[92];

    auto g_xx_0_0_0_0_xz_yy_yy = buffer_2000_sddd[93];

    auto g_xx_0_0_0_0_xz_yy_yz = buffer_2000_sddd[94];

    auto g_xx_0_0_0_0_xz_yy_zz = buffer_2000_sddd[95];

    auto g_xx_0_0_0_0_xz_yz_xx = buffer_2000_sddd[96];

    auto g_xx_0_0_0_0_xz_yz_xy = buffer_2000_sddd[97];

    auto g_xx_0_0_0_0_xz_yz_xz = buffer_2000_sddd[98];

    auto g_xx_0_0_0_0_xz_yz_yy = buffer_2000_sddd[99];

    auto g_xx_0_0_0_0_xz_yz_yz = buffer_2000_sddd[100];

    auto g_xx_0_0_0_0_xz_yz_zz = buffer_2000_sddd[101];

    auto g_xx_0_0_0_0_xz_zz_xx = buffer_2000_sddd[102];

    auto g_xx_0_0_0_0_xz_zz_xy = buffer_2000_sddd[103];

    auto g_xx_0_0_0_0_xz_zz_xz = buffer_2000_sddd[104];

    auto g_xx_0_0_0_0_xz_zz_yy = buffer_2000_sddd[105];

    auto g_xx_0_0_0_0_xz_zz_yz = buffer_2000_sddd[106];

    auto g_xx_0_0_0_0_xz_zz_zz = buffer_2000_sddd[107];

    auto g_xx_0_0_0_0_yy_xx_xx = buffer_2000_sddd[108];

    auto g_xx_0_0_0_0_yy_xx_xy = buffer_2000_sddd[109];

    auto g_xx_0_0_0_0_yy_xx_xz = buffer_2000_sddd[110];

    auto g_xx_0_0_0_0_yy_xx_yy = buffer_2000_sddd[111];

    auto g_xx_0_0_0_0_yy_xx_yz = buffer_2000_sddd[112];

    auto g_xx_0_0_0_0_yy_xx_zz = buffer_2000_sddd[113];

    auto g_xx_0_0_0_0_yy_xy_xx = buffer_2000_sddd[114];

    auto g_xx_0_0_0_0_yy_xy_xy = buffer_2000_sddd[115];

    auto g_xx_0_0_0_0_yy_xy_xz = buffer_2000_sddd[116];

    auto g_xx_0_0_0_0_yy_xy_yy = buffer_2000_sddd[117];

    auto g_xx_0_0_0_0_yy_xy_yz = buffer_2000_sddd[118];

    auto g_xx_0_0_0_0_yy_xy_zz = buffer_2000_sddd[119];

    auto g_xx_0_0_0_0_yy_xz_xx = buffer_2000_sddd[120];

    auto g_xx_0_0_0_0_yy_xz_xy = buffer_2000_sddd[121];

    auto g_xx_0_0_0_0_yy_xz_xz = buffer_2000_sddd[122];

    auto g_xx_0_0_0_0_yy_xz_yy = buffer_2000_sddd[123];

    auto g_xx_0_0_0_0_yy_xz_yz = buffer_2000_sddd[124];

    auto g_xx_0_0_0_0_yy_xz_zz = buffer_2000_sddd[125];

    auto g_xx_0_0_0_0_yy_yy_xx = buffer_2000_sddd[126];

    auto g_xx_0_0_0_0_yy_yy_xy = buffer_2000_sddd[127];

    auto g_xx_0_0_0_0_yy_yy_xz = buffer_2000_sddd[128];

    auto g_xx_0_0_0_0_yy_yy_yy = buffer_2000_sddd[129];

    auto g_xx_0_0_0_0_yy_yy_yz = buffer_2000_sddd[130];

    auto g_xx_0_0_0_0_yy_yy_zz = buffer_2000_sddd[131];

    auto g_xx_0_0_0_0_yy_yz_xx = buffer_2000_sddd[132];

    auto g_xx_0_0_0_0_yy_yz_xy = buffer_2000_sddd[133];

    auto g_xx_0_0_0_0_yy_yz_xz = buffer_2000_sddd[134];

    auto g_xx_0_0_0_0_yy_yz_yy = buffer_2000_sddd[135];

    auto g_xx_0_0_0_0_yy_yz_yz = buffer_2000_sddd[136];

    auto g_xx_0_0_0_0_yy_yz_zz = buffer_2000_sddd[137];

    auto g_xx_0_0_0_0_yy_zz_xx = buffer_2000_sddd[138];

    auto g_xx_0_0_0_0_yy_zz_xy = buffer_2000_sddd[139];

    auto g_xx_0_0_0_0_yy_zz_xz = buffer_2000_sddd[140];

    auto g_xx_0_0_0_0_yy_zz_yy = buffer_2000_sddd[141];

    auto g_xx_0_0_0_0_yy_zz_yz = buffer_2000_sddd[142];

    auto g_xx_0_0_0_0_yy_zz_zz = buffer_2000_sddd[143];

    auto g_xx_0_0_0_0_yz_xx_xx = buffer_2000_sddd[144];

    auto g_xx_0_0_0_0_yz_xx_xy = buffer_2000_sddd[145];

    auto g_xx_0_0_0_0_yz_xx_xz = buffer_2000_sddd[146];

    auto g_xx_0_0_0_0_yz_xx_yy = buffer_2000_sddd[147];

    auto g_xx_0_0_0_0_yz_xx_yz = buffer_2000_sddd[148];

    auto g_xx_0_0_0_0_yz_xx_zz = buffer_2000_sddd[149];

    auto g_xx_0_0_0_0_yz_xy_xx = buffer_2000_sddd[150];

    auto g_xx_0_0_0_0_yz_xy_xy = buffer_2000_sddd[151];

    auto g_xx_0_0_0_0_yz_xy_xz = buffer_2000_sddd[152];

    auto g_xx_0_0_0_0_yz_xy_yy = buffer_2000_sddd[153];

    auto g_xx_0_0_0_0_yz_xy_yz = buffer_2000_sddd[154];

    auto g_xx_0_0_0_0_yz_xy_zz = buffer_2000_sddd[155];

    auto g_xx_0_0_0_0_yz_xz_xx = buffer_2000_sddd[156];

    auto g_xx_0_0_0_0_yz_xz_xy = buffer_2000_sddd[157];

    auto g_xx_0_0_0_0_yz_xz_xz = buffer_2000_sddd[158];

    auto g_xx_0_0_0_0_yz_xz_yy = buffer_2000_sddd[159];

    auto g_xx_0_0_0_0_yz_xz_yz = buffer_2000_sddd[160];

    auto g_xx_0_0_0_0_yz_xz_zz = buffer_2000_sddd[161];

    auto g_xx_0_0_0_0_yz_yy_xx = buffer_2000_sddd[162];

    auto g_xx_0_0_0_0_yz_yy_xy = buffer_2000_sddd[163];

    auto g_xx_0_0_0_0_yz_yy_xz = buffer_2000_sddd[164];

    auto g_xx_0_0_0_0_yz_yy_yy = buffer_2000_sddd[165];

    auto g_xx_0_0_0_0_yz_yy_yz = buffer_2000_sddd[166];

    auto g_xx_0_0_0_0_yz_yy_zz = buffer_2000_sddd[167];

    auto g_xx_0_0_0_0_yz_yz_xx = buffer_2000_sddd[168];

    auto g_xx_0_0_0_0_yz_yz_xy = buffer_2000_sddd[169];

    auto g_xx_0_0_0_0_yz_yz_xz = buffer_2000_sddd[170];

    auto g_xx_0_0_0_0_yz_yz_yy = buffer_2000_sddd[171];

    auto g_xx_0_0_0_0_yz_yz_yz = buffer_2000_sddd[172];

    auto g_xx_0_0_0_0_yz_yz_zz = buffer_2000_sddd[173];

    auto g_xx_0_0_0_0_yz_zz_xx = buffer_2000_sddd[174];

    auto g_xx_0_0_0_0_yz_zz_xy = buffer_2000_sddd[175];

    auto g_xx_0_0_0_0_yz_zz_xz = buffer_2000_sddd[176];

    auto g_xx_0_0_0_0_yz_zz_yy = buffer_2000_sddd[177];

    auto g_xx_0_0_0_0_yz_zz_yz = buffer_2000_sddd[178];

    auto g_xx_0_0_0_0_yz_zz_zz = buffer_2000_sddd[179];

    auto g_xx_0_0_0_0_zz_xx_xx = buffer_2000_sddd[180];

    auto g_xx_0_0_0_0_zz_xx_xy = buffer_2000_sddd[181];

    auto g_xx_0_0_0_0_zz_xx_xz = buffer_2000_sddd[182];

    auto g_xx_0_0_0_0_zz_xx_yy = buffer_2000_sddd[183];

    auto g_xx_0_0_0_0_zz_xx_yz = buffer_2000_sddd[184];

    auto g_xx_0_0_0_0_zz_xx_zz = buffer_2000_sddd[185];

    auto g_xx_0_0_0_0_zz_xy_xx = buffer_2000_sddd[186];

    auto g_xx_0_0_0_0_zz_xy_xy = buffer_2000_sddd[187];

    auto g_xx_0_0_0_0_zz_xy_xz = buffer_2000_sddd[188];

    auto g_xx_0_0_0_0_zz_xy_yy = buffer_2000_sddd[189];

    auto g_xx_0_0_0_0_zz_xy_yz = buffer_2000_sddd[190];

    auto g_xx_0_0_0_0_zz_xy_zz = buffer_2000_sddd[191];

    auto g_xx_0_0_0_0_zz_xz_xx = buffer_2000_sddd[192];

    auto g_xx_0_0_0_0_zz_xz_xy = buffer_2000_sddd[193];

    auto g_xx_0_0_0_0_zz_xz_xz = buffer_2000_sddd[194];

    auto g_xx_0_0_0_0_zz_xz_yy = buffer_2000_sddd[195];

    auto g_xx_0_0_0_0_zz_xz_yz = buffer_2000_sddd[196];

    auto g_xx_0_0_0_0_zz_xz_zz = buffer_2000_sddd[197];

    auto g_xx_0_0_0_0_zz_yy_xx = buffer_2000_sddd[198];

    auto g_xx_0_0_0_0_zz_yy_xy = buffer_2000_sddd[199];

    auto g_xx_0_0_0_0_zz_yy_xz = buffer_2000_sddd[200];

    auto g_xx_0_0_0_0_zz_yy_yy = buffer_2000_sddd[201];

    auto g_xx_0_0_0_0_zz_yy_yz = buffer_2000_sddd[202];

    auto g_xx_0_0_0_0_zz_yy_zz = buffer_2000_sddd[203];

    auto g_xx_0_0_0_0_zz_yz_xx = buffer_2000_sddd[204];

    auto g_xx_0_0_0_0_zz_yz_xy = buffer_2000_sddd[205];

    auto g_xx_0_0_0_0_zz_yz_xz = buffer_2000_sddd[206];

    auto g_xx_0_0_0_0_zz_yz_yy = buffer_2000_sddd[207];

    auto g_xx_0_0_0_0_zz_yz_yz = buffer_2000_sddd[208];

    auto g_xx_0_0_0_0_zz_yz_zz = buffer_2000_sddd[209];

    auto g_xx_0_0_0_0_zz_zz_xx = buffer_2000_sddd[210];

    auto g_xx_0_0_0_0_zz_zz_xy = buffer_2000_sddd[211];

    auto g_xx_0_0_0_0_zz_zz_xz = buffer_2000_sddd[212];

    auto g_xx_0_0_0_0_zz_zz_yy = buffer_2000_sddd[213];

    auto g_xx_0_0_0_0_zz_zz_yz = buffer_2000_sddd[214];

    auto g_xx_0_0_0_0_zz_zz_zz = buffer_2000_sddd[215];

    auto g_xy_0_0_0_0_xx_xx_xx = buffer_2000_sddd[216];

    auto g_xy_0_0_0_0_xx_xx_xy = buffer_2000_sddd[217];

    auto g_xy_0_0_0_0_xx_xx_xz = buffer_2000_sddd[218];

    auto g_xy_0_0_0_0_xx_xx_yy = buffer_2000_sddd[219];

    auto g_xy_0_0_0_0_xx_xx_yz = buffer_2000_sddd[220];

    auto g_xy_0_0_0_0_xx_xx_zz = buffer_2000_sddd[221];

    auto g_xy_0_0_0_0_xx_xy_xx = buffer_2000_sddd[222];

    auto g_xy_0_0_0_0_xx_xy_xy = buffer_2000_sddd[223];

    auto g_xy_0_0_0_0_xx_xy_xz = buffer_2000_sddd[224];

    auto g_xy_0_0_0_0_xx_xy_yy = buffer_2000_sddd[225];

    auto g_xy_0_0_0_0_xx_xy_yz = buffer_2000_sddd[226];

    auto g_xy_0_0_0_0_xx_xy_zz = buffer_2000_sddd[227];

    auto g_xy_0_0_0_0_xx_xz_xx = buffer_2000_sddd[228];

    auto g_xy_0_0_0_0_xx_xz_xy = buffer_2000_sddd[229];

    auto g_xy_0_0_0_0_xx_xz_xz = buffer_2000_sddd[230];

    auto g_xy_0_0_0_0_xx_xz_yy = buffer_2000_sddd[231];

    auto g_xy_0_0_0_0_xx_xz_yz = buffer_2000_sddd[232];

    auto g_xy_0_0_0_0_xx_xz_zz = buffer_2000_sddd[233];

    auto g_xy_0_0_0_0_xx_yy_xx = buffer_2000_sddd[234];

    auto g_xy_0_0_0_0_xx_yy_xy = buffer_2000_sddd[235];

    auto g_xy_0_0_0_0_xx_yy_xz = buffer_2000_sddd[236];

    auto g_xy_0_0_0_0_xx_yy_yy = buffer_2000_sddd[237];

    auto g_xy_0_0_0_0_xx_yy_yz = buffer_2000_sddd[238];

    auto g_xy_0_0_0_0_xx_yy_zz = buffer_2000_sddd[239];

    auto g_xy_0_0_0_0_xx_yz_xx = buffer_2000_sddd[240];

    auto g_xy_0_0_0_0_xx_yz_xy = buffer_2000_sddd[241];

    auto g_xy_0_0_0_0_xx_yz_xz = buffer_2000_sddd[242];

    auto g_xy_0_0_0_0_xx_yz_yy = buffer_2000_sddd[243];

    auto g_xy_0_0_0_0_xx_yz_yz = buffer_2000_sddd[244];

    auto g_xy_0_0_0_0_xx_yz_zz = buffer_2000_sddd[245];

    auto g_xy_0_0_0_0_xx_zz_xx = buffer_2000_sddd[246];

    auto g_xy_0_0_0_0_xx_zz_xy = buffer_2000_sddd[247];

    auto g_xy_0_0_0_0_xx_zz_xz = buffer_2000_sddd[248];

    auto g_xy_0_0_0_0_xx_zz_yy = buffer_2000_sddd[249];

    auto g_xy_0_0_0_0_xx_zz_yz = buffer_2000_sddd[250];

    auto g_xy_0_0_0_0_xx_zz_zz = buffer_2000_sddd[251];

    auto g_xy_0_0_0_0_xy_xx_xx = buffer_2000_sddd[252];

    auto g_xy_0_0_0_0_xy_xx_xy = buffer_2000_sddd[253];

    auto g_xy_0_0_0_0_xy_xx_xz = buffer_2000_sddd[254];

    auto g_xy_0_0_0_0_xy_xx_yy = buffer_2000_sddd[255];

    auto g_xy_0_0_0_0_xy_xx_yz = buffer_2000_sddd[256];

    auto g_xy_0_0_0_0_xy_xx_zz = buffer_2000_sddd[257];

    auto g_xy_0_0_0_0_xy_xy_xx = buffer_2000_sddd[258];

    auto g_xy_0_0_0_0_xy_xy_xy = buffer_2000_sddd[259];

    auto g_xy_0_0_0_0_xy_xy_xz = buffer_2000_sddd[260];

    auto g_xy_0_0_0_0_xy_xy_yy = buffer_2000_sddd[261];

    auto g_xy_0_0_0_0_xy_xy_yz = buffer_2000_sddd[262];

    auto g_xy_0_0_0_0_xy_xy_zz = buffer_2000_sddd[263];

    auto g_xy_0_0_0_0_xy_xz_xx = buffer_2000_sddd[264];

    auto g_xy_0_0_0_0_xy_xz_xy = buffer_2000_sddd[265];

    auto g_xy_0_0_0_0_xy_xz_xz = buffer_2000_sddd[266];

    auto g_xy_0_0_0_0_xy_xz_yy = buffer_2000_sddd[267];

    auto g_xy_0_0_0_0_xy_xz_yz = buffer_2000_sddd[268];

    auto g_xy_0_0_0_0_xy_xz_zz = buffer_2000_sddd[269];

    auto g_xy_0_0_0_0_xy_yy_xx = buffer_2000_sddd[270];

    auto g_xy_0_0_0_0_xy_yy_xy = buffer_2000_sddd[271];

    auto g_xy_0_0_0_0_xy_yy_xz = buffer_2000_sddd[272];

    auto g_xy_0_0_0_0_xy_yy_yy = buffer_2000_sddd[273];

    auto g_xy_0_0_0_0_xy_yy_yz = buffer_2000_sddd[274];

    auto g_xy_0_0_0_0_xy_yy_zz = buffer_2000_sddd[275];

    auto g_xy_0_0_0_0_xy_yz_xx = buffer_2000_sddd[276];

    auto g_xy_0_0_0_0_xy_yz_xy = buffer_2000_sddd[277];

    auto g_xy_0_0_0_0_xy_yz_xz = buffer_2000_sddd[278];

    auto g_xy_0_0_0_0_xy_yz_yy = buffer_2000_sddd[279];

    auto g_xy_0_0_0_0_xy_yz_yz = buffer_2000_sddd[280];

    auto g_xy_0_0_0_0_xy_yz_zz = buffer_2000_sddd[281];

    auto g_xy_0_0_0_0_xy_zz_xx = buffer_2000_sddd[282];

    auto g_xy_0_0_0_0_xy_zz_xy = buffer_2000_sddd[283];

    auto g_xy_0_0_0_0_xy_zz_xz = buffer_2000_sddd[284];

    auto g_xy_0_0_0_0_xy_zz_yy = buffer_2000_sddd[285];

    auto g_xy_0_0_0_0_xy_zz_yz = buffer_2000_sddd[286];

    auto g_xy_0_0_0_0_xy_zz_zz = buffer_2000_sddd[287];

    auto g_xy_0_0_0_0_xz_xx_xx = buffer_2000_sddd[288];

    auto g_xy_0_0_0_0_xz_xx_xy = buffer_2000_sddd[289];

    auto g_xy_0_0_0_0_xz_xx_xz = buffer_2000_sddd[290];

    auto g_xy_0_0_0_0_xz_xx_yy = buffer_2000_sddd[291];

    auto g_xy_0_0_0_0_xz_xx_yz = buffer_2000_sddd[292];

    auto g_xy_0_0_0_0_xz_xx_zz = buffer_2000_sddd[293];

    auto g_xy_0_0_0_0_xz_xy_xx = buffer_2000_sddd[294];

    auto g_xy_0_0_0_0_xz_xy_xy = buffer_2000_sddd[295];

    auto g_xy_0_0_0_0_xz_xy_xz = buffer_2000_sddd[296];

    auto g_xy_0_0_0_0_xz_xy_yy = buffer_2000_sddd[297];

    auto g_xy_0_0_0_0_xz_xy_yz = buffer_2000_sddd[298];

    auto g_xy_0_0_0_0_xz_xy_zz = buffer_2000_sddd[299];

    auto g_xy_0_0_0_0_xz_xz_xx = buffer_2000_sddd[300];

    auto g_xy_0_0_0_0_xz_xz_xy = buffer_2000_sddd[301];

    auto g_xy_0_0_0_0_xz_xz_xz = buffer_2000_sddd[302];

    auto g_xy_0_0_0_0_xz_xz_yy = buffer_2000_sddd[303];

    auto g_xy_0_0_0_0_xz_xz_yz = buffer_2000_sddd[304];

    auto g_xy_0_0_0_0_xz_xz_zz = buffer_2000_sddd[305];

    auto g_xy_0_0_0_0_xz_yy_xx = buffer_2000_sddd[306];

    auto g_xy_0_0_0_0_xz_yy_xy = buffer_2000_sddd[307];

    auto g_xy_0_0_0_0_xz_yy_xz = buffer_2000_sddd[308];

    auto g_xy_0_0_0_0_xz_yy_yy = buffer_2000_sddd[309];

    auto g_xy_0_0_0_0_xz_yy_yz = buffer_2000_sddd[310];

    auto g_xy_0_0_0_0_xz_yy_zz = buffer_2000_sddd[311];

    auto g_xy_0_0_0_0_xz_yz_xx = buffer_2000_sddd[312];

    auto g_xy_0_0_0_0_xz_yz_xy = buffer_2000_sddd[313];

    auto g_xy_0_0_0_0_xz_yz_xz = buffer_2000_sddd[314];

    auto g_xy_0_0_0_0_xz_yz_yy = buffer_2000_sddd[315];

    auto g_xy_0_0_0_0_xz_yz_yz = buffer_2000_sddd[316];

    auto g_xy_0_0_0_0_xz_yz_zz = buffer_2000_sddd[317];

    auto g_xy_0_0_0_0_xz_zz_xx = buffer_2000_sddd[318];

    auto g_xy_0_0_0_0_xz_zz_xy = buffer_2000_sddd[319];

    auto g_xy_0_0_0_0_xz_zz_xz = buffer_2000_sddd[320];

    auto g_xy_0_0_0_0_xz_zz_yy = buffer_2000_sddd[321];

    auto g_xy_0_0_0_0_xz_zz_yz = buffer_2000_sddd[322];

    auto g_xy_0_0_0_0_xz_zz_zz = buffer_2000_sddd[323];

    auto g_xy_0_0_0_0_yy_xx_xx = buffer_2000_sddd[324];

    auto g_xy_0_0_0_0_yy_xx_xy = buffer_2000_sddd[325];

    auto g_xy_0_0_0_0_yy_xx_xz = buffer_2000_sddd[326];

    auto g_xy_0_0_0_0_yy_xx_yy = buffer_2000_sddd[327];

    auto g_xy_0_0_0_0_yy_xx_yz = buffer_2000_sddd[328];

    auto g_xy_0_0_0_0_yy_xx_zz = buffer_2000_sddd[329];

    auto g_xy_0_0_0_0_yy_xy_xx = buffer_2000_sddd[330];

    auto g_xy_0_0_0_0_yy_xy_xy = buffer_2000_sddd[331];

    auto g_xy_0_0_0_0_yy_xy_xz = buffer_2000_sddd[332];

    auto g_xy_0_0_0_0_yy_xy_yy = buffer_2000_sddd[333];

    auto g_xy_0_0_0_0_yy_xy_yz = buffer_2000_sddd[334];

    auto g_xy_0_0_0_0_yy_xy_zz = buffer_2000_sddd[335];

    auto g_xy_0_0_0_0_yy_xz_xx = buffer_2000_sddd[336];

    auto g_xy_0_0_0_0_yy_xz_xy = buffer_2000_sddd[337];

    auto g_xy_0_0_0_0_yy_xz_xz = buffer_2000_sddd[338];

    auto g_xy_0_0_0_0_yy_xz_yy = buffer_2000_sddd[339];

    auto g_xy_0_0_0_0_yy_xz_yz = buffer_2000_sddd[340];

    auto g_xy_0_0_0_0_yy_xz_zz = buffer_2000_sddd[341];

    auto g_xy_0_0_0_0_yy_yy_xx = buffer_2000_sddd[342];

    auto g_xy_0_0_0_0_yy_yy_xy = buffer_2000_sddd[343];

    auto g_xy_0_0_0_0_yy_yy_xz = buffer_2000_sddd[344];

    auto g_xy_0_0_0_0_yy_yy_yy = buffer_2000_sddd[345];

    auto g_xy_0_0_0_0_yy_yy_yz = buffer_2000_sddd[346];

    auto g_xy_0_0_0_0_yy_yy_zz = buffer_2000_sddd[347];

    auto g_xy_0_0_0_0_yy_yz_xx = buffer_2000_sddd[348];

    auto g_xy_0_0_0_0_yy_yz_xy = buffer_2000_sddd[349];

    auto g_xy_0_0_0_0_yy_yz_xz = buffer_2000_sddd[350];

    auto g_xy_0_0_0_0_yy_yz_yy = buffer_2000_sddd[351];

    auto g_xy_0_0_0_0_yy_yz_yz = buffer_2000_sddd[352];

    auto g_xy_0_0_0_0_yy_yz_zz = buffer_2000_sddd[353];

    auto g_xy_0_0_0_0_yy_zz_xx = buffer_2000_sddd[354];

    auto g_xy_0_0_0_0_yy_zz_xy = buffer_2000_sddd[355];

    auto g_xy_0_0_0_0_yy_zz_xz = buffer_2000_sddd[356];

    auto g_xy_0_0_0_0_yy_zz_yy = buffer_2000_sddd[357];

    auto g_xy_0_0_0_0_yy_zz_yz = buffer_2000_sddd[358];

    auto g_xy_0_0_0_0_yy_zz_zz = buffer_2000_sddd[359];

    auto g_xy_0_0_0_0_yz_xx_xx = buffer_2000_sddd[360];

    auto g_xy_0_0_0_0_yz_xx_xy = buffer_2000_sddd[361];

    auto g_xy_0_0_0_0_yz_xx_xz = buffer_2000_sddd[362];

    auto g_xy_0_0_0_0_yz_xx_yy = buffer_2000_sddd[363];

    auto g_xy_0_0_0_0_yz_xx_yz = buffer_2000_sddd[364];

    auto g_xy_0_0_0_0_yz_xx_zz = buffer_2000_sddd[365];

    auto g_xy_0_0_0_0_yz_xy_xx = buffer_2000_sddd[366];

    auto g_xy_0_0_0_0_yz_xy_xy = buffer_2000_sddd[367];

    auto g_xy_0_0_0_0_yz_xy_xz = buffer_2000_sddd[368];

    auto g_xy_0_0_0_0_yz_xy_yy = buffer_2000_sddd[369];

    auto g_xy_0_0_0_0_yz_xy_yz = buffer_2000_sddd[370];

    auto g_xy_0_0_0_0_yz_xy_zz = buffer_2000_sddd[371];

    auto g_xy_0_0_0_0_yz_xz_xx = buffer_2000_sddd[372];

    auto g_xy_0_0_0_0_yz_xz_xy = buffer_2000_sddd[373];

    auto g_xy_0_0_0_0_yz_xz_xz = buffer_2000_sddd[374];

    auto g_xy_0_0_0_0_yz_xz_yy = buffer_2000_sddd[375];

    auto g_xy_0_0_0_0_yz_xz_yz = buffer_2000_sddd[376];

    auto g_xy_0_0_0_0_yz_xz_zz = buffer_2000_sddd[377];

    auto g_xy_0_0_0_0_yz_yy_xx = buffer_2000_sddd[378];

    auto g_xy_0_0_0_0_yz_yy_xy = buffer_2000_sddd[379];

    auto g_xy_0_0_0_0_yz_yy_xz = buffer_2000_sddd[380];

    auto g_xy_0_0_0_0_yz_yy_yy = buffer_2000_sddd[381];

    auto g_xy_0_0_0_0_yz_yy_yz = buffer_2000_sddd[382];

    auto g_xy_0_0_0_0_yz_yy_zz = buffer_2000_sddd[383];

    auto g_xy_0_0_0_0_yz_yz_xx = buffer_2000_sddd[384];

    auto g_xy_0_0_0_0_yz_yz_xy = buffer_2000_sddd[385];

    auto g_xy_0_0_0_0_yz_yz_xz = buffer_2000_sddd[386];

    auto g_xy_0_0_0_0_yz_yz_yy = buffer_2000_sddd[387];

    auto g_xy_0_0_0_0_yz_yz_yz = buffer_2000_sddd[388];

    auto g_xy_0_0_0_0_yz_yz_zz = buffer_2000_sddd[389];

    auto g_xy_0_0_0_0_yz_zz_xx = buffer_2000_sddd[390];

    auto g_xy_0_0_0_0_yz_zz_xy = buffer_2000_sddd[391];

    auto g_xy_0_0_0_0_yz_zz_xz = buffer_2000_sddd[392];

    auto g_xy_0_0_0_0_yz_zz_yy = buffer_2000_sddd[393];

    auto g_xy_0_0_0_0_yz_zz_yz = buffer_2000_sddd[394];

    auto g_xy_0_0_0_0_yz_zz_zz = buffer_2000_sddd[395];

    auto g_xy_0_0_0_0_zz_xx_xx = buffer_2000_sddd[396];

    auto g_xy_0_0_0_0_zz_xx_xy = buffer_2000_sddd[397];

    auto g_xy_0_0_0_0_zz_xx_xz = buffer_2000_sddd[398];

    auto g_xy_0_0_0_0_zz_xx_yy = buffer_2000_sddd[399];

    auto g_xy_0_0_0_0_zz_xx_yz = buffer_2000_sddd[400];

    auto g_xy_0_0_0_0_zz_xx_zz = buffer_2000_sddd[401];

    auto g_xy_0_0_0_0_zz_xy_xx = buffer_2000_sddd[402];

    auto g_xy_0_0_0_0_zz_xy_xy = buffer_2000_sddd[403];

    auto g_xy_0_0_0_0_zz_xy_xz = buffer_2000_sddd[404];

    auto g_xy_0_0_0_0_zz_xy_yy = buffer_2000_sddd[405];

    auto g_xy_0_0_0_0_zz_xy_yz = buffer_2000_sddd[406];

    auto g_xy_0_0_0_0_zz_xy_zz = buffer_2000_sddd[407];

    auto g_xy_0_0_0_0_zz_xz_xx = buffer_2000_sddd[408];

    auto g_xy_0_0_0_0_zz_xz_xy = buffer_2000_sddd[409];

    auto g_xy_0_0_0_0_zz_xz_xz = buffer_2000_sddd[410];

    auto g_xy_0_0_0_0_zz_xz_yy = buffer_2000_sddd[411];

    auto g_xy_0_0_0_0_zz_xz_yz = buffer_2000_sddd[412];

    auto g_xy_0_0_0_0_zz_xz_zz = buffer_2000_sddd[413];

    auto g_xy_0_0_0_0_zz_yy_xx = buffer_2000_sddd[414];

    auto g_xy_0_0_0_0_zz_yy_xy = buffer_2000_sddd[415];

    auto g_xy_0_0_0_0_zz_yy_xz = buffer_2000_sddd[416];

    auto g_xy_0_0_0_0_zz_yy_yy = buffer_2000_sddd[417];

    auto g_xy_0_0_0_0_zz_yy_yz = buffer_2000_sddd[418];

    auto g_xy_0_0_0_0_zz_yy_zz = buffer_2000_sddd[419];

    auto g_xy_0_0_0_0_zz_yz_xx = buffer_2000_sddd[420];

    auto g_xy_0_0_0_0_zz_yz_xy = buffer_2000_sddd[421];

    auto g_xy_0_0_0_0_zz_yz_xz = buffer_2000_sddd[422];

    auto g_xy_0_0_0_0_zz_yz_yy = buffer_2000_sddd[423];

    auto g_xy_0_0_0_0_zz_yz_yz = buffer_2000_sddd[424];

    auto g_xy_0_0_0_0_zz_yz_zz = buffer_2000_sddd[425];

    auto g_xy_0_0_0_0_zz_zz_xx = buffer_2000_sddd[426];

    auto g_xy_0_0_0_0_zz_zz_xy = buffer_2000_sddd[427];

    auto g_xy_0_0_0_0_zz_zz_xz = buffer_2000_sddd[428];

    auto g_xy_0_0_0_0_zz_zz_yy = buffer_2000_sddd[429];

    auto g_xy_0_0_0_0_zz_zz_yz = buffer_2000_sddd[430];

    auto g_xy_0_0_0_0_zz_zz_zz = buffer_2000_sddd[431];

    auto g_xz_0_0_0_0_xx_xx_xx = buffer_2000_sddd[432];

    auto g_xz_0_0_0_0_xx_xx_xy = buffer_2000_sddd[433];

    auto g_xz_0_0_0_0_xx_xx_xz = buffer_2000_sddd[434];

    auto g_xz_0_0_0_0_xx_xx_yy = buffer_2000_sddd[435];

    auto g_xz_0_0_0_0_xx_xx_yz = buffer_2000_sddd[436];

    auto g_xz_0_0_0_0_xx_xx_zz = buffer_2000_sddd[437];

    auto g_xz_0_0_0_0_xx_xy_xx = buffer_2000_sddd[438];

    auto g_xz_0_0_0_0_xx_xy_xy = buffer_2000_sddd[439];

    auto g_xz_0_0_0_0_xx_xy_xz = buffer_2000_sddd[440];

    auto g_xz_0_0_0_0_xx_xy_yy = buffer_2000_sddd[441];

    auto g_xz_0_0_0_0_xx_xy_yz = buffer_2000_sddd[442];

    auto g_xz_0_0_0_0_xx_xy_zz = buffer_2000_sddd[443];

    auto g_xz_0_0_0_0_xx_xz_xx = buffer_2000_sddd[444];

    auto g_xz_0_0_0_0_xx_xz_xy = buffer_2000_sddd[445];

    auto g_xz_0_0_0_0_xx_xz_xz = buffer_2000_sddd[446];

    auto g_xz_0_0_0_0_xx_xz_yy = buffer_2000_sddd[447];

    auto g_xz_0_0_0_0_xx_xz_yz = buffer_2000_sddd[448];

    auto g_xz_0_0_0_0_xx_xz_zz = buffer_2000_sddd[449];

    auto g_xz_0_0_0_0_xx_yy_xx = buffer_2000_sddd[450];

    auto g_xz_0_0_0_0_xx_yy_xy = buffer_2000_sddd[451];

    auto g_xz_0_0_0_0_xx_yy_xz = buffer_2000_sddd[452];

    auto g_xz_0_0_0_0_xx_yy_yy = buffer_2000_sddd[453];

    auto g_xz_0_0_0_0_xx_yy_yz = buffer_2000_sddd[454];

    auto g_xz_0_0_0_0_xx_yy_zz = buffer_2000_sddd[455];

    auto g_xz_0_0_0_0_xx_yz_xx = buffer_2000_sddd[456];

    auto g_xz_0_0_0_0_xx_yz_xy = buffer_2000_sddd[457];

    auto g_xz_0_0_0_0_xx_yz_xz = buffer_2000_sddd[458];

    auto g_xz_0_0_0_0_xx_yz_yy = buffer_2000_sddd[459];

    auto g_xz_0_0_0_0_xx_yz_yz = buffer_2000_sddd[460];

    auto g_xz_0_0_0_0_xx_yz_zz = buffer_2000_sddd[461];

    auto g_xz_0_0_0_0_xx_zz_xx = buffer_2000_sddd[462];

    auto g_xz_0_0_0_0_xx_zz_xy = buffer_2000_sddd[463];

    auto g_xz_0_0_0_0_xx_zz_xz = buffer_2000_sddd[464];

    auto g_xz_0_0_0_0_xx_zz_yy = buffer_2000_sddd[465];

    auto g_xz_0_0_0_0_xx_zz_yz = buffer_2000_sddd[466];

    auto g_xz_0_0_0_0_xx_zz_zz = buffer_2000_sddd[467];

    auto g_xz_0_0_0_0_xy_xx_xx = buffer_2000_sddd[468];

    auto g_xz_0_0_0_0_xy_xx_xy = buffer_2000_sddd[469];

    auto g_xz_0_0_0_0_xy_xx_xz = buffer_2000_sddd[470];

    auto g_xz_0_0_0_0_xy_xx_yy = buffer_2000_sddd[471];

    auto g_xz_0_0_0_0_xy_xx_yz = buffer_2000_sddd[472];

    auto g_xz_0_0_0_0_xy_xx_zz = buffer_2000_sddd[473];

    auto g_xz_0_0_0_0_xy_xy_xx = buffer_2000_sddd[474];

    auto g_xz_0_0_0_0_xy_xy_xy = buffer_2000_sddd[475];

    auto g_xz_0_0_0_0_xy_xy_xz = buffer_2000_sddd[476];

    auto g_xz_0_0_0_0_xy_xy_yy = buffer_2000_sddd[477];

    auto g_xz_0_0_0_0_xy_xy_yz = buffer_2000_sddd[478];

    auto g_xz_0_0_0_0_xy_xy_zz = buffer_2000_sddd[479];

    auto g_xz_0_0_0_0_xy_xz_xx = buffer_2000_sddd[480];

    auto g_xz_0_0_0_0_xy_xz_xy = buffer_2000_sddd[481];

    auto g_xz_0_0_0_0_xy_xz_xz = buffer_2000_sddd[482];

    auto g_xz_0_0_0_0_xy_xz_yy = buffer_2000_sddd[483];

    auto g_xz_0_0_0_0_xy_xz_yz = buffer_2000_sddd[484];

    auto g_xz_0_0_0_0_xy_xz_zz = buffer_2000_sddd[485];

    auto g_xz_0_0_0_0_xy_yy_xx = buffer_2000_sddd[486];

    auto g_xz_0_0_0_0_xy_yy_xy = buffer_2000_sddd[487];

    auto g_xz_0_0_0_0_xy_yy_xz = buffer_2000_sddd[488];

    auto g_xz_0_0_0_0_xy_yy_yy = buffer_2000_sddd[489];

    auto g_xz_0_0_0_0_xy_yy_yz = buffer_2000_sddd[490];

    auto g_xz_0_0_0_0_xy_yy_zz = buffer_2000_sddd[491];

    auto g_xz_0_0_0_0_xy_yz_xx = buffer_2000_sddd[492];

    auto g_xz_0_0_0_0_xy_yz_xy = buffer_2000_sddd[493];

    auto g_xz_0_0_0_0_xy_yz_xz = buffer_2000_sddd[494];

    auto g_xz_0_0_0_0_xy_yz_yy = buffer_2000_sddd[495];

    auto g_xz_0_0_0_0_xy_yz_yz = buffer_2000_sddd[496];

    auto g_xz_0_0_0_0_xy_yz_zz = buffer_2000_sddd[497];

    auto g_xz_0_0_0_0_xy_zz_xx = buffer_2000_sddd[498];

    auto g_xz_0_0_0_0_xy_zz_xy = buffer_2000_sddd[499];

    auto g_xz_0_0_0_0_xy_zz_xz = buffer_2000_sddd[500];

    auto g_xz_0_0_0_0_xy_zz_yy = buffer_2000_sddd[501];

    auto g_xz_0_0_0_0_xy_zz_yz = buffer_2000_sddd[502];

    auto g_xz_0_0_0_0_xy_zz_zz = buffer_2000_sddd[503];

    auto g_xz_0_0_0_0_xz_xx_xx = buffer_2000_sddd[504];

    auto g_xz_0_0_0_0_xz_xx_xy = buffer_2000_sddd[505];

    auto g_xz_0_0_0_0_xz_xx_xz = buffer_2000_sddd[506];

    auto g_xz_0_0_0_0_xz_xx_yy = buffer_2000_sddd[507];

    auto g_xz_0_0_0_0_xz_xx_yz = buffer_2000_sddd[508];

    auto g_xz_0_0_0_0_xz_xx_zz = buffer_2000_sddd[509];

    auto g_xz_0_0_0_0_xz_xy_xx = buffer_2000_sddd[510];

    auto g_xz_0_0_0_0_xz_xy_xy = buffer_2000_sddd[511];

    auto g_xz_0_0_0_0_xz_xy_xz = buffer_2000_sddd[512];

    auto g_xz_0_0_0_0_xz_xy_yy = buffer_2000_sddd[513];

    auto g_xz_0_0_0_0_xz_xy_yz = buffer_2000_sddd[514];

    auto g_xz_0_0_0_0_xz_xy_zz = buffer_2000_sddd[515];

    auto g_xz_0_0_0_0_xz_xz_xx = buffer_2000_sddd[516];

    auto g_xz_0_0_0_0_xz_xz_xy = buffer_2000_sddd[517];

    auto g_xz_0_0_0_0_xz_xz_xz = buffer_2000_sddd[518];

    auto g_xz_0_0_0_0_xz_xz_yy = buffer_2000_sddd[519];

    auto g_xz_0_0_0_0_xz_xz_yz = buffer_2000_sddd[520];

    auto g_xz_0_0_0_0_xz_xz_zz = buffer_2000_sddd[521];

    auto g_xz_0_0_0_0_xz_yy_xx = buffer_2000_sddd[522];

    auto g_xz_0_0_0_0_xz_yy_xy = buffer_2000_sddd[523];

    auto g_xz_0_0_0_0_xz_yy_xz = buffer_2000_sddd[524];

    auto g_xz_0_0_0_0_xz_yy_yy = buffer_2000_sddd[525];

    auto g_xz_0_0_0_0_xz_yy_yz = buffer_2000_sddd[526];

    auto g_xz_0_0_0_0_xz_yy_zz = buffer_2000_sddd[527];

    auto g_xz_0_0_0_0_xz_yz_xx = buffer_2000_sddd[528];

    auto g_xz_0_0_0_0_xz_yz_xy = buffer_2000_sddd[529];

    auto g_xz_0_0_0_0_xz_yz_xz = buffer_2000_sddd[530];

    auto g_xz_0_0_0_0_xz_yz_yy = buffer_2000_sddd[531];

    auto g_xz_0_0_0_0_xz_yz_yz = buffer_2000_sddd[532];

    auto g_xz_0_0_0_0_xz_yz_zz = buffer_2000_sddd[533];

    auto g_xz_0_0_0_0_xz_zz_xx = buffer_2000_sddd[534];

    auto g_xz_0_0_0_0_xz_zz_xy = buffer_2000_sddd[535];

    auto g_xz_0_0_0_0_xz_zz_xz = buffer_2000_sddd[536];

    auto g_xz_0_0_0_0_xz_zz_yy = buffer_2000_sddd[537];

    auto g_xz_0_0_0_0_xz_zz_yz = buffer_2000_sddd[538];

    auto g_xz_0_0_0_0_xz_zz_zz = buffer_2000_sddd[539];

    auto g_xz_0_0_0_0_yy_xx_xx = buffer_2000_sddd[540];

    auto g_xz_0_0_0_0_yy_xx_xy = buffer_2000_sddd[541];

    auto g_xz_0_0_0_0_yy_xx_xz = buffer_2000_sddd[542];

    auto g_xz_0_0_0_0_yy_xx_yy = buffer_2000_sddd[543];

    auto g_xz_0_0_0_0_yy_xx_yz = buffer_2000_sddd[544];

    auto g_xz_0_0_0_0_yy_xx_zz = buffer_2000_sddd[545];

    auto g_xz_0_0_0_0_yy_xy_xx = buffer_2000_sddd[546];

    auto g_xz_0_0_0_0_yy_xy_xy = buffer_2000_sddd[547];

    auto g_xz_0_0_0_0_yy_xy_xz = buffer_2000_sddd[548];

    auto g_xz_0_0_0_0_yy_xy_yy = buffer_2000_sddd[549];

    auto g_xz_0_0_0_0_yy_xy_yz = buffer_2000_sddd[550];

    auto g_xz_0_0_0_0_yy_xy_zz = buffer_2000_sddd[551];

    auto g_xz_0_0_0_0_yy_xz_xx = buffer_2000_sddd[552];

    auto g_xz_0_0_0_0_yy_xz_xy = buffer_2000_sddd[553];

    auto g_xz_0_0_0_0_yy_xz_xz = buffer_2000_sddd[554];

    auto g_xz_0_0_0_0_yy_xz_yy = buffer_2000_sddd[555];

    auto g_xz_0_0_0_0_yy_xz_yz = buffer_2000_sddd[556];

    auto g_xz_0_0_0_0_yy_xz_zz = buffer_2000_sddd[557];

    auto g_xz_0_0_0_0_yy_yy_xx = buffer_2000_sddd[558];

    auto g_xz_0_0_0_0_yy_yy_xy = buffer_2000_sddd[559];

    auto g_xz_0_0_0_0_yy_yy_xz = buffer_2000_sddd[560];

    auto g_xz_0_0_0_0_yy_yy_yy = buffer_2000_sddd[561];

    auto g_xz_0_0_0_0_yy_yy_yz = buffer_2000_sddd[562];

    auto g_xz_0_0_0_0_yy_yy_zz = buffer_2000_sddd[563];

    auto g_xz_0_0_0_0_yy_yz_xx = buffer_2000_sddd[564];

    auto g_xz_0_0_0_0_yy_yz_xy = buffer_2000_sddd[565];

    auto g_xz_0_0_0_0_yy_yz_xz = buffer_2000_sddd[566];

    auto g_xz_0_0_0_0_yy_yz_yy = buffer_2000_sddd[567];

    auto g_xz_0_0_0_0_yy_yz_yz = buffer_2000_sddd[568];

    auto g_xz_0_0_0_0_yy_yz_zz = buffer_2000_sddd[569];

    auto g_xz_0_0_0_0_yy_zz_xx = buffer_2000_sddd[570];

    auto g_xz_0_0_0_0_yy_zz_xy = buffer_2000_sddd[571];

    auto g_xz_0_0_0_0_yy_zz_xz = buffer_2000_sddd[572];

    auto g_xz_0_0_0_0_yy_zz_yy = buffer_2000_sddd[573];

    auto g_xz_0_0_0_0_yy_zz_yz = buffer_2000_sddd[574];

    auto g_xz_0_0_0_0_yy_zz_zz = buffer_2000_sddd[575];

    auto g_xz_0_0_0_0_yz_xx_xx = buffer_2000_sddd[576];

    auto g_xz_0_0_0_0_yz_xx_xy = buffer_2000_sddd[577];

    auto g_xz_0_0_0_0_yz_xx_xz = buffer_2000_sddd[578];

    auto g_xz_0_0_0_0_yz_xx_yy = buffer_2000_sddd[579];

    auto g_xz_0_0_0_0_yz_xx_yz = buffer_2000_sddd[580];

    auto g_xz_0_0_0_0_yz_xx_zz = buffer_2000_sddd[581];

    auto g_xz_0_0_0_0_yz_xy_xx = buffer_2000_sddd[582];

    auto g_xz_0_0_0_0_yz_xy_xy = buffer_2000_sddd[583];

    auto g_xz_0_0_0_0_yz_xy_xz = buffer_2000_sddd[584];

    auto g_xz_0_0_0_0_yz_xy_yy = buffer_2000_sddd[585];

    auto g_xz_0_0_0_0_yz_xy_yz = buffer_2000_sddd[586];

    auto g_xz_0_0_0_0_yz_xy_zz = buffer_2000_sddd[587];

    auto g_xz_0_0_0_0_yz_xz_xx = buffer_2000_sddd[588];

    auto g_xz_0_0_0_0_yz_xz_xy = buffer_2000_sddd[589];

    auto g_xz_0_0_0_0_yz_xz_xz = buffer_2000_sddd[590];

    auto g_xz_0_0_0_0_yz_xz_yy = buffer_2000_sddd[591];

    auto g_xz_0_0_0_0_yz_xz_yz = buffer_2000_sddd[592];

    auto g_xz_0_0_0_0_yz_xz_zz = buffer_2000_sddd[593];

    auto g_xz_0_0_0_0_yz_yy_xx = buffer_2000_sddd[594];

    auto g_xz_0_0_0_0_yz_yy_xy = buffer_2000_sddd[595];

    auto g_xz_0_0_0_0_yz_yy_xz = buffer_2000_sddd[596];

    auto g_xz_0_0_0_0_yz_yy_yy = buffer_2000_sddd[597];

    auto g_xz_0_0_0_0_yz_yy_yz = buffer_2000_sddd[598];

    auto g_xz_0_0_0_0_yz_yy_zz = buffer_2000_sddd[599];

    auto g_xz_0_0_0_0_yz_yz_xx = buffer_2000_sddd[600];

    auto g_xz_0_0_0_0_yz_yz_xy = buffer_2000_sddd[601];

    auto g_xz_0_0_0_0_yz_yz_xz = buffer_2000_sddd[602];

    auto g_xz_0_0_0_0_yz_yz_yy = buffer_2000_sddd[603];

    auto g_xz_0_0_0_0_yz_yz_yz = buffer_2000_sddd[604];

    auto g_xz_0_0_0_0_yz_yz_zz = buffer_2000_sddd[605];

    auto g_xz_0_0_0_0_yz_zz_xx = buffer_2000_sddd[606];

    auto g_xz_0_0_0_0_yz_zz_xy = buffer_2000_sddd[607];

    auto g_xz_0_0_0_0_yz_zz_xz = buffer_2000_sddd[608];

    auto g_xz_0_0_0_0_yz_zz_yy = buffer_2000_sddd[609];

    auto g_xz_0_0_0_0_yz_zz_yz = buffer_2000_sddd[610];

    auto g_xz_0_0_0_0_yz_zz_zz = buffer_2000_sddd[611];

    auto g_xz_0_0_0_0_zz_xx_xx = buffer_2000_sddd[612];

    auto g_xz_0_0_0_0_zz_xx_xy = buffer_2000_sddd[613];

    auto g_xz_0_0_0_0_zz_xx_xz = buffer_2000_sddd[614];

    auto g_xz_0_0_0_0_zz_xx_yy = buffer_2000_sddd[615];

    auto g_xz_0_0_0_0_zz_xx_yz = buffer_2000_sddd[616];

    auto g_xz_0_0_0_0_zz_xx_zz = buffer_2000_sddd[617];

    auto g_xz_0_0_0_0_zz_xy_xx = buffer_2000_sddd[618];

    auto g_xz_0_0_0_0_zz_xy_xy = buffer_2000_sddd[619];

    auto g_xz_0_0_0_0_zz_xy_xz = buffer_2000_sddd[620];

    auto g_xz_0_0_0_0_zz_xy_yy = buffer_2000_sddd[621];

    auto g_xz_0_0_0_0_zz_xy_yz = buffer_2000_sddd[622];

    auto g_xz_0_0_0_0_zz_xy_zz = buffer_2000_sddd[623];

    auto g_xz_0_0_0_0_zz_xz_xx = buffer_2000_sddd[624];

    auto g_xz_0_0_0_0_zz_xz_xy = buffer_2000_sddd[625];

    auto g_xz_0_0_0_0_zz_xz_xz = buffer_2000_sddd[626];

    auto g_xz_0_0_0_0_zz_xz_yy = buffer_2000_sddd[627];

    auto g_xz_0_0_0_0_zz_xz_yz = buffer_2000_sddd[628];

    auto g_xz_0_0_0_0_zz_xz_zz = buffer_2000_sddd[629];

    auto g_xz_0_0_0_0_zz_yy_xx = buffer_2000_sddd[630];

    auto g_xz_0_0_0_0_zz_yy_xy = buffer_2000_sddd[631];

    auto g_xz_0_0_0_0_zz_yy_xz = buffer_2000_sddd[632];

    auto g_xz_0_0_0_0_zz_yy_yy = buffer_2000_sddd[633];

    auto g_xz_0_0_0_0_zz_yy_yz = buffer_2000_sddd[634];

    auto g_xz_0_0_0_0_zz_yy_zz = buffer_2000_sddd[635];

    auto g_xz_0_0_0_0_zz_yz_xx = buffer_2000_sddd[636];

    auto g_xz_0_0_0_0_zz_yz_xy = buffer_2000_sddd[637];

    auto g_xz_0_0_0_0_zz_yz_xz = buffer_2000_sddd[638];

    auto g_xz_0_0_0_0_zz_yz_yy = buffer_2000_sddd[639];

    auto g_xz_0_0_0_0_zz_yz_yz = buffer_2000_sddd[640];

    auto g_xz_0_0_0_0_zz_yz_zz = buffer_2000_sddd[641];

    auto g_xz_0_0_0_0_zz_zz_xx = buffer_2000_sddd[642];

    auto g_xz_0_0_0_0_zz_zz_xy = buffer_2000_sddd[643];

    auto g_xz_0_0_0_0_zz_zz_xz = buffer_2000_sddd[644];

    auto g_xz_0_0_0_0_zz_zz_yy = buffer_2000_sddd[645];

    auto g_xz_0_0_0_0_zz_zz_yz = buffer_2000_sddd[646];

    auto g_xz_0_0_0_0_zz_zz_zz = buffer_2000_sddd[647];

    auto g_yy_0_0_0_0_xx_xx_xx = buffer_2000_sddd[648];

    auto g_yy_0_0_0_0_xx_xx_xy = buffer_2000_sddd[649];

    auto g_yy_0_0_0_0_xx_xx_xz = buffer_2000_sddd[650];

    auto g_yy_0_0_0_0_xx_xx_yy = buffer_2000_sddd[651];

    auto g_yy_0_0_0_0_xx_xx_yz = buffer_2000_sddd[652];

    auto g_yy_0_0_0_0_xx_xx_zz = buffer_2000_sddd[653];

    auto g_yy_0_0_0_0_xx_xy_xx = buffer_2000_sddd[654];

    auto g_yy_0_0_0_0_xx_xy_xy = buffer_2000_sddd[655];

    auto g_yy_0_0_0_0_xx_xy_xz = buffer_2000_sddd[656];

    auto g_yy_0_0_0_0_xx_xy_yy = buffer_2000_sddd[657];

    auto g_yy_0_0_0_0_xx_xy_yz = buffer_2000_sddd[658];

    auto g_yy_0_0_0_0_xx_xy_zz = buffer_2000_sddd[659];

    auto g_yy_0_0_0_0_xx_xz_xx = buffer_2000_sddd[660];

    auto g_yy_0_0_0_0_xx_xz_xy = buffer_2000_sddd[661];

    auto g_yy_0_0_0_0_xx_xz_xz = buffer_2000_sddd[662];

    auto g_yy_0_0_0_0_xx_xz_yy = buffer_2000_sddd[663];

    auto g_yy_0_0_0_0_xx_xz_yz = buffer_2000_sddd[664];

    auto g_yy_0_0_0_0_xx_xz_zz = buffer_2000_sddd[665];

    auto g_yy_0_0_0_0_xx_yy_xx = buffer_2000_sddd[666];

    auto g_yy_0_0_0_0_xx_yy_xy = buffer_2000_sddd[667];

    auto g_yy_0_0_0_0_xx_yy_xz = buffer_2000_sddd[668];

    auto g_yy_0_0_0_0_xx_yy_yy = buffer_2000_sddd[669];

    auto g_yy_0_0_0_0_xx_yy_yz = buffer_2000_sddd[670];

    auto g_yy_0_0_0_0_xx_yy_zz = buffer_2000_sddd[671];

    auto g_yy_0_0_0_0_xx_yz_xx = buffer_2000_sddd[672];

    auto g_yy_0_0_0_0_xx_yz_xy = buffer_2000_sddd[673];

    auto g_yy_0_0_0_0_xx_yz_xz = buffer_2000_sddd[674];

    auto g_yy_0_0_0_0_xx_yz_yy = buffer_2000_sddd[675];

    auto g_yy_0_0_0_0_xx_yz_yz = buffer_2000_sddd[676];

    auto g_yy_0_0_0_0_xx_yz_zz = buffer_2000_sddd[677];

    auto g_yy_0_0_0_0_xx_zz_xx = buffer_2000_sddd[678];

    auto g_yy_0_0_0_0_xx_zz_xy = buffer_2000_sddd[679];

    auto g_yy_0_0_0_0_xx_zz_xz = buffer_2000_sddd[680];

    auto g_yy_0_0_0_0_xx_zz_yy = buffer_2000_sddd[681];

    auto g_yy_0_0_0_0_xx_zz_yz = buffer_2000_sddd[682];

    auto g_yy_0_0_0_0_xx_zz_zz = buffer_2000_sddd[683];

    auto g_yy_0_0_0_0_xy_xx_xx = buffer_2000_sddd[684];

    auto g_yy_0_0_0_0_xy_xx_xy = buffer_2000_sddd[685];

    auto g_yy_0_0_0_0_xy_xx_xz = buffer_2000_sddd[686];

    auto g_yy_0_0_0_0_xy_xx_yy = buffer_2000_sddd[687];

    auto g_yy_0_0_0_0_xy_xx_yz = buffer_2000_sddd[688];

    auto g_yy_0_0_0_0_xy_xx_zz = buffer_2000_sddd[689];

    auto g_yy_0_0_0_0_xy_xy_xx = buffer_2000_sddd[690];

    auto g_yy_0_0_0_0_xy_xy_xy = buffer_2000_sddd[691];

    auto g_yy_0_0_0_0_xy_xy_xz = buffer_2000_sddd[692];

    auto g_yy_0_0_0_0_xy_xy_yy = buffer_2000_sddd[693];

    auto g_yy_0_0_0_0_xy_xy_yz = buffer_2000_sddd[694];

    auto g_yy_0_0_0_0_xy_xy_zz = buffer_2000_sddd[695];

    auto g_yy_0_0_0_0_xy_xz_xx = buffer_2000_sddd[696];

    auto g_yy_0_0_0_0_xy_xz_xy = buffer_2000_sddd[697];

    auto g_yy_0_0_0_0_xy_xz_xz = buffer_2000_sddd[698];

    auto g_yy_0_0_0_0_xy_xz_yy = buffer_2000_sddd[699];

    auto g_yy_0_0_0_0_xy_xz_yz = buffer_2000_sddd[700];

    auto g_yy_0_0_0_0_xy_xz_zz = buffer_2000_sddd[701];

    auto g_yy_0_0_0_0_xy_yy_xx = buffer_2000_sddd[702];

    auto g_yy_0_0_0_0_xy_yy_xy = buffer_2000_sddd[703];

    auto g_yy_0_0_0_0_xy_yy_xz = buffer_2000_sddd[704];

    auto g_yy_0_0_0_0_xy_yy_yy = buffer_2000_sddd[705];

    auto g_yy_0_0_0_0_xy_yy_yz = buffer_2000_sddd[706];

    auto g_yy_0_0_0_0_xy_yy_zz = buffer_2000_sddd[707];

    auto g_yy_0_0_0_0_xy_yz_xx = buffer_2000_sddd[708];

    auto g_yy_0_0_0_0_xy_yz_xy = buffer_2000_sddd[709];

    auto g_yy_0_0_0_0_xy_yz_xz = buffer_2000_sddd[710];

    auto g_yy_0_0_0_0_xy_yz_yy = buffer_2000_sddd[711];

    auto g_yy_0_0_0_0_xy_yz_yz = buffer_2000_sddd[712];

    auto g_yy_0_0_0_0_xy_yz_zz = buffer_2000_sddd[713];

    auto g_yy_0_0_0_0_xy_zz_xx = buffer_2000_sddd[714];

    auto g_yy_0_0_0_0_xy_zz_xy = buffer_2000_sddd[715];

    auto g_yy_0_0_0_0_xy_zz_xz = buffer_2000_sddd[716];

    auto g_yy_0_0_0_0_xy_zz_yy = buffer_2000_sddd[717];

    auto g_yy_0_0_0_0_xy_zz_yz = buffer_2000_sddd[718];

    auto g_yy_0_0_0_0_xy_zz_zz = buffer_2000_sddd[719];

    auto g_yy_0_0_0_0_xz_xx_xx = buffer_2000_sddd[720];

    auto g_yy_0_0_0_0_xz_xx_xy = buffer_2000_sddd[721];

    auto g_yy_0_0_0_0_xz_xx_xz = buffer_2000_sddd[722];

    auto g_yy_0_0_0_0_xz_xx_yy = buffer_2000_sddd[723];

    auto g_yy_0_0_0_0_xz_xx_yz = buffer_2000_sddd[724];

    auto g_yy_0_0_0_0_xz_xx_zz = buffer_2000_sddd[725];

    auto g_yy_0_0_0_0_xz_xy_xx = buffer_2000_sddd[726];

    auto g_yy_0_0_0_0_xz_xy_xy = buffer_2000_sddd[727];

    auto g_yy_0_0_0_0_xz_xy_xz = buffer_2000_sddd[728];

    auto g_yy_0_0_0_0_xz_xy_yy = buffer_2000_sddd[729];

    auto g_yy_0_0_0_0_xz_xy_yz = buffer_2000_sddd[730];

    auto g_yy_0_0_0_0_xz_xy_zz = buffer_2000_sddd[731];

    auto g_yy_0_0_0_0_xz_xz_xx = buffer_2000_sddd[732];

    auto g_yy_0_0_0_0_xz_xz_xy = buffer_2000_sddd[733];

    auto g_yy_0_0_0_0_xz_xz_xz = buffer_2000_sddd[734];

    auto g_yy_0_0_0_0_xz_xz_yy = buffer_2000_sddd[735];

    auto g_yy_0_0_0_0_xz_xz_yz = buffer_2000_sddd[736];

    auto g_yy_0_0_0_0_xz_xz_zz = buffer_2000_sddd[737];

    auto g_yy_0_0_0_0_xz_yy_xx = buffer_2000_sddd[738];

    auto g_yy_0_0_0_0_xz_yy_xy = buffer_2000_sddd[739];

    auto g_yy_0_0_0_0_xz_yy_xz = buffer_2000_sddd[740];

    auto g_yy_0_0_0_0_xz_yy_yy = buffer_2000_sddd[741];

    auto g_yy_0_0_0_0_xz_yy_yz = buffer_2000_sddd[742];

    auto g_yy_0_0_0_0_xz_yy_zz = buffer_2000_sddd[743];

    auto g_yy_0_0_0_0_xz_yz_xx = buffer_2000_sddd[744];

    auto g_yy_0_0_0_0_xz_yz_xy = buffer_2000_sddd[745];

    auto g_yy_0_0_0_0_xz_yz_xz = buffer_2000_sddd[746];

    auto g_yy_0_0_0_0_xz_yz_yy = buffer_2000_sddd[747];

    auto g_yy_0_0_0_0_xz_yz_yz = buffer_2000_sddd[748];

    auto g_yy_0_0_0_0_xz_yz_zz = buffer_2000_sddd[749];

    auto g_yy_0_0_0_0_xz_zz_xx = buffer_2000_sddd[750];

    auto g_yy_0_0_0_0_xz_zz_xy = buffer_2000_sddd[751];

    auto g_yy_0_0_0_0_xz_zz_xz = buffer_2000_sddd[752];

    auto g_yy_0_0_0_0_xz_zz_yy = buffer_2000_sddd[753];

    auto g_yy_0_0_0_0_xz_zz_yz = buffer_2000_sddd[754];

    auto g_yy_0_0_0_0_xz_zz_zz = buffer_2000_sddd[755];

    auto g_yy_0_0_0_0_yy_xx_xx = buffer_2000_sddd[756];

    auto g_yy_0_0_0_0_yy_xx_xy = buffer_2000_sddd[757];

    auto g_yy_0_0_0_0_yy_xx_xz = buffer_2000_sddd[758];

    auto g_yy_0_0_0_0_yy_xx_yy = buffer_2000_sddd[759];

    auto g_yy_0_0_0_0_yy_xx_yz = buffer_2000_sddd[760];

    auto g_yy_0_0_0_0_yy_xx_zz = buffer_2000_sddd[761];

    auto g_yy_0_0_0_0_yy_xy_xx = buffer_2000_sddd[762];

    auto g_yy_0_0_0_0_yy_xy_xy = buffer_2000_sddd[763];

    auto g_yy_0_0_0_0_yy_xy_xz = buffer_2000_sddd[764];

    auto g_yy_0_0_0_0_yy_xy_yy = buffer_2000_sddd[765];

    auto g_yy_0_0_0_0_yy_xy_yz = buffer_2000_sddd[766];

    auto g_yy_0_0_0_0_yy_xy_zz = buffer_2000_sddd[767];

    auto g_yy_0_0_0_0_yy_xz_xx = buffer_2000_sddd[768];

    auto g_yy_0_0_0_0_yy_xz_xy = buffer_2000_sddd[769];

    auto g_yy_0_0_0_0_yy_xz_xz = buffer_2000_sddd[770];

    auto g_yy_0_0_0_0_yy_xz_yy = buffer_2000_sddd[771];

    auto g_yy_0_0_0_0_yy_xz_yz = buffer_2000_sddd[772];

    auto g_yy_0_0_0_0_yy_xz_zz = buffer_2000_sddd[773];

    auto g_yy_0_0_0_0_yy_yy_xx = buffer_2000_sddd[774];

    auto g_yy_0_0_0_0_yy_yy_xy = buffer_2000_sddd[775];

    auto g_yy_0_0_0_0_yy_yy_xz = buffer_2000_sddd[776];

    auto g_yy_0_0_0_0_yy_yy_yy = buffer_2000_sddd[777];

    auto g_yy_0_0_0_0_yy_yy_yz = buffer_2000_sddd[778];

    auto g_yy_0_0_0_0_yy_yy_zz = buffer_2000_sddd[779];

    auto g_yy_0_0_0_0_yy_yz_xx = buffer_2000_sddd[780];

    auto g_yy_0_0_0_0_yy_yz_xy = buffer_2000_sddd[781];

    auto g_yy_0_0_0_0_yy_yz_xz = buffer_2000_sddd[782];

    auto g_yy_0_0_0_0_yy_yz_yy = buffer_2000_sddd[783];

    auto g_yy_0_0_0_0_yy_yz_yz = buffer_2000_sddd[784];

    auto g_yy_0_0_0_0_yy_yz_zz = buffer_2000_sddd[785];

    auto g_yy_0_0_0_0_yy_zz_xx = buffer_2000_sddd[786];

    auto g_yy_0_0_0_0_yy_zz_xy = buffer_2000_sddd[787];

    auto g_yy_0_0_0_0_yy_zz_xz = buffer_2000_sddd[788];

    auto g_yy_0_0_0_0_yy_zz_yy = buffer_2000_sddd[789];

    auto g_yy_0_0_0_0_yy_zz_yz = buffer_2000_sddd[790];

    auto g_yy_0_0_0_0_yy_zz_zz = buffer_2000_sddd[791];

    auto g_yy_0_0_0_0_yz_xx_xx = buffer_2000_sddd[792];

    auto g_yy_0_0_0_0_yz_xx_xy = buffer_2000_sddd[793];

    auto g_yy_0_0_0_0_yz_xx_xz = buffer_2000_sddd[794];

    auto g_yy_0_0_0_0_yz_xx_yy = buffer_2000_sddd[795];

    auto g_yy_0_0_0_0_yz_xx_yz = buffer_2000_sddd[796];

    auto g_yy_0_0_0_0_yz_xx_zz = buffer_2000_sddd[797];

    auto g_yy_0_0_0_0_yz_xy_xx = buffer_2000_sddd[798];

    auto g_yy_0_0_0_0_yz_xy_xy = buffer_2000_sddd[799];

    auto g_yy_0_0_0_0_yz_xy_xz = buffer_2000_sddd[800];

    auto g_yy_0_0_0_0_yz_xy_yy = buffer_2000_sddd[801];

    auto g_yy_0_0_0_0_yz_xy_yz = buffer_2000_sddd[802];

    auto g_yy_0_0_0_0_yz_xy_zz = buffer_2000_sddd[803];

    auto g_yy_0_0_0_0_yz_xz_xx = buffer_2000_sddd[804];

    auto g_yy_0_0_0_0_yz_xz_xy = buffer_2000_sddd[805];

    auto g_yy_0_0_0_0_yz_xz_xz = buffer_2000_sddd[806];

    auto g_yy_0_0_0_0_yz_xz_yy = buffer_2000_sddd[807];

    auto g_yy_0_0_0_0_yz_xz_yz = buffer_2000_sddd[808];

    auto g_yy_0_0_0_0_yz_xz_zz = buffer_2000_sddd[809];

    auto g_yy_0_0_0_0_yz_yy_xx = buffer_2000_sddd[810];

    auto g_yy_0_0_0_0_yz_yy_xy = buffer_2000_sddd[811];

    auto g_yy_0_0_0_0_yz_yy_xz = buffer_2000_sddd[812];

    auto g_yy_0_0_0_0_yz_yy_yy = buffer_2000_sddd[813];

    auto g_yy_0_0_0_0_yz_yy_yz = buffer_2000_sddd[814];

    auto g_yy_0_0_0_0_yz_yy_zz = buffer_2000_sddd[815];

    auto g_yy_0_0_0_0_yz_yz_xx = buffer_2000_sddd[816];

    auto g_yy_0_0_0_0_yz_yz_xy = buffer_2000_sddd[817];

    auto g_yy_0_0_0_0_yz_yz_xz = buffer_2000_sddd[818];

    auto g_yy_0_0_0_0_yz_yz_yy = buffer_2000_sddd[819];

    auto g_yy_0_0_0_0_yz_yz_yz = buffer_2000_sddd[820];

    auto g_yy_0_0_0_0_yz_yz_zz = buffer_2000_sddd[821];

    auto g_yy_0_0_0_0_yz_zz_xx = buffer_2000_sddd[822];

    auto g_yy_0_0_0_0_yz_zz_xy = buffer_2000_sddd[823];

    auto g_yy_0_0_0_0_yz_zz_xz = buffer_2000_sddd[824];

    auto g_yy_0_0_0_0_yz_zz_yy = buffer_2000_sddd[825];

    auto g_yy_0_0_0_0_yz_zz_yz = buffer_2000_sddd[826];

    auto g_yy_0_0_0_0_yz_zz_zz = buffer_2000_sddd[827];

    auto g_yy_0_0_0_0_zz_xx_xx = buffer_2000_sddd[828];

    auto g_yy_0_0_0_0_zz_xx_xy = buffer_2000_sddd[829];

    auto g_yy_0_0_0_0_zz_xx_xz = buffer_2000_sddd[830];

    auto g_yy_0_0_0_0_zz_xx_yy = buffer_2000_sddd[831];

    auto g_yy_0_0_0_0_zz_xx_yz = buffer_2000_sddd[832];

    auto g_yy_0_0_0_0_zz_xx_zz = buffer_2000_sddd[833];

    auto g_yy_0_0_0_0_zz_xy_xx = buffer_2000_sddd[834];

    auto g_yy_0_0_0_0_zz_xy_xy = buffer_2000_sddd[835];

    auto g_yy_0_0_0_0_zz_xy_xz = buffer_2000_sddd[836];

    auto g_yy_0_0_0_0_zz_xy_yy = buffer_2000_sddd[837];

    auto g_yy_0_0_0_0_zz_xy_yz = buffer_2000_sddd[838];

    auto g_yy_0_0_0_0_zz_xy_zz = buffer_2000_sddd[839];

    auto g_yy_0_0_0_0_zz_xz_xx = buffer_2000_sddd[840];

    auto g_yy_0_0_0_0_zz_xz_xy = buffer_2000_sddd[841];

    auto g_yy_0_0_0_0_zz_xz_xz = buffer_2000_sddd[842];

    auto g_yy_0_0_0_0_zz_xz_yy = buffer_2000_sddd[843];

    auto g_yy_0_0_0_0_zz_xz_yz = buffer_2000_sddd[844];

    auto g_yy_0_0_0_0_zz_xz_zz = buffer_2000_sddd[845];

    auto g_yy_0_0_0_0_zz_yy_xx = buffer_2000_sddd[846];

    auto g_yy_0_0_0_0_zz_yy_xy = buffer_2000_sddd[847];

    auto g_yy_0_0_0_0_zz_yy_xz = buffer_2000_sddd[848];

    auto g_yy_0_0_0_0_zz_yy_yy = buffer_2000_sddd[849];

    auto g_yy_0_0_0_0_zz_yy_yz = buffer_2000_sddd[850];

    auto g_yy_0_0_0_0_zz_yy_zz = buffer_2000_sddd[851];

    auto g_yy_0_0_0_0_zz_yz_xx = buffer_2000_sddd[852];

    auto g_yy_0_0_0_0_zz_yz_xy = buffer_2000_sddd[853];

    auto g_yy_0_0_0_0_zz_yz_xz = buffer_2000_sddd[854];

    auto g_yy_0_0_0_0_zz_yz_yy = buffer_2000_sddd[855];

    auto g_yy_0_0_0_0_zz_yz_yz = buffer_2000_sddd[856];

    auto g_yy_0_0_0_0_zz_yz_zz = buffer_2000_sddd[857];

    auto g_yy_0_0_0_0_zz_zz_xx = buffer_2000_sddd[858];

    auto g_yy_0_0_0_0_zz_zz_xy = buffer_2000_sddd[859];

    auto g_yy_0_0_0_0_zz_zz_xz = buffer_2000_sddd[860];

    auto g_yy_0_0_0_0_zz_zz_yy = buffer_2000_sddd[861];

    auto g_yy_0_0_0_0_zz_zz_yz = buffer_2000_sddd[862];

    auto g_yy_0_0_0_0_zz_zz_zz = buffer_2000_sddd[863];

    auto g_yz_0_0_0_0_xx_xx_xx = buffer_2000_sddd[864];

    auto g_yz_0_0_0_0_xx_xx_xy = buffer_2000_sddd[865];

    auto g_yz_0_0_0_0_xx_xx_xz = buffer_2000_sddd[866];

    auto g_yz_0_0_0_0_xx_xx_yy = buffer_2000_sddd[867];

    auto g_yz_0_0_0_0_xx_xx_yz = buffer_2000_sddd[868];

    auto g_yz_0_0_0_0_xx_xx_zz = buffer_2000_sddd[869];

    auto g_yz_0_0_0_0_xx_xy_xx = buffer_2000_sddd[870];

    auto g_yz_0_0_0_0_xx_xy_xy = buffer_2000_sddd[871];

    auto g_yz_0_0_0_0_xx_xy_xz = buffer_2000_sddd[872];

    auto g_yz_0_0_0_0_xx_xy_yy = buffer_2000_sddd[873];

    auto g_yz_0_0_0_0_xx_xy_yz = buffer_2000_sddd[874];

    auto g_yz_0_0_0_0_xx_xy_zz = buffer_2000_sddd[875];

    auto g_yz_0_0_0_0_xx_xz_xx = buffer_2000_sddd[876];

    auto g_yz_0_0_0_0_xx_xz_xy = buffer_2000_sddd[877];

    auto g_yz_0_0_0_0_xx_xz_xz = buffer_2000_sddd[878];

    auto g_yz_0_0_0_0_xx_xz_yy = buffer_2000_sddd[879];

    auto g_yz_0_0_0_0_xx_xz_yz = buffer_2000_sddd[880];

    auto g_yz_0_0_0_0_xx_xz_zz = buffer_2000_sddd[881];

    auto g_yz_0_0_0_0_xx_yy_xx = buffer_2000_sddd[882];

    auto g_yz_0_0_0_0_xx_yy_xy = buffer_2000_sddd[883];

    auto g_yz_0_0_0_0_xx_yy_xz = buffer_2000_sddd[884];

    auto g_yz_0_0_0_0_xx_yy_yy = buffer_2000_sddd[885];

    auto g_yz_0_0_0_0_xx_yy_yz = buffer_2000_sddd[886];

    auto g_yz_0_0_0_0_xx_yy_zz = buffer_2000_sddd[887];

    auto g_yz_0_0_0_0_xx_yz_xx = buffer_2000_sddd[888];

    auto g_yz_0_0_0_0_xx_yz_xy = buffer_2000_sddd[889];

    auto g_yz_0_0_0_0_xx_yz_xz = buffer_2000_sddd[890];

    auto g_yz_0_0_0_0_xx_yz_yy = buffer_2000_sddd[891];

    auto g_yz_0_0_0_0_xx_yz_yz = buffer_2000_sddd[892];

    auto g_yz_0_0_0_0_xx_yz_zz = buffer_2000_sddd[893];

    auto g_yz_0_0_0_0_xx_zz_xx = buffer_2000_sddd[894];

    auto g_yz_0_0_0_0_xx_zz_xy = buffer_2000_sddd[895];

    auto g_yz_0_0_0_0_xx_zz_xz = buffer_2000_sddd[896];

    auto g_yz_0_0_0_0_xx_zz_yy = buffer_2000_sddd[897];

    auto g_yz_0_0_0_0_xx_zz_yz = buffer_2000_sddd[898];

    auto g_yz_0_0_0_0_xx_zz_zz = buffer_2000_sddd[899];

    auto g_yz_0_0_0_0_xy_xx_xx = buffer_2000_sddd[900];

    auto g_yz_0_0_0_0_xy_xx_xy = buffer_2000_sddd[901];

    auto g_yz_0_0_0_0_xy_xx_xz = buffer_2000_sddd[902];

    auto g_yz_0_0_0_0_xy_xx_yy = buffer_2000_sddd[903];

    auto g_yz_0_0_0_0_xy_xx_yz = buffer_2000_sddd[904];

    auto g_yz_0_0_0_0_xy_xx_zz = buffer_2000_sddd[905];

    auto g_yz_0_0_0_0_xy_xy_xx = buffer_2000_sddd[906];

    auto g_yz_0_0_0_0_xy_xy_xy = buffer_2000_sddd[907];

    auto g_yz_0_0_0_0_xy_xy_xz = buffer_2000_sddd[908];

    auto g_yz_0_0_0_0_xy_xy_yy = buffer_2000_sddd[909];

    auto g_yz_0_0_0_0_xy_xy_yz = buffer_2000_sddd[910];

    auto g_yz_0_0_0_0_xy_xy_zz = buffer_2000_sddd[911];

    auto g_yz_0_0_0_0_xy_xz_xx = buffer_2000_sddd[912];

    auto g_yz_0_0_0_0_xy_xz_xy = buffer_2000_sddd[913];

    auto g_yz_0_0_0_0_xy_xz_xz = buffer_2000_sddd[914];

    auto g_yz_0_0_0_0_xy_xz_yy = buffer_2000_sddd[915];

    auto g_yz_0_0_0_0_xy_xz_yz = buffer_2000_sddd[916];

    auto g_yz_0_0_0_0_xy_xz_zz = buffer_2000_sddd[917];

    auto g_yz_0_0_0_0_xy_yy_xx = buffer_2000_sddd[918];

    auto g_yz_0_0_0_0_xy_yy_xy = buffer_2000_sddd[919];

    auto g_yz_0_0_0_0_xy_yy_xz = buffer_2000_sddd[920];

    auto g_yz_0_0_0_0_xy_yy_yy = buffer_2000_sddd[921];

    auto g_yz_0_0_0_0_xy_yy_yz = buffer_2000_sddd[922];

    auto g_yz_0_0_0_0_xy_yy_zz = buffer_2000_sddd[923];

    auto g_yz_0_0_0_0_xy_yz_xx = buffer_2000_sddd[924];

    auto g_yz_0_0_0_0_xy_yz_xy = buffer_2000_sddd[925];

    auto g_yz_0_0_0_0_xy_yz_xz = buffer_2000_sddd[926];

    auto g_yz_0_0_0_0_xy_yz_yy = buffer_2000_sddd[927];

    auto g_yz_0_0_0_0_xy_yz_yz = buffer_2000_sddd[928];

    auto g_yz_0_0_0_0_xy_yz_zz = buffer_2000_sddd[929];

    auto g_yz_0_0_0_0_xy_zz_xx = buffer_2000_sddd[930];

    auto g_yz_0_0_0_0_xy_zz_xy = buffer_2000_sddd[931];

    auto g_yz_0_0_0_0_xy_zz_xz = buffer_2000_sddd[932];

    auto g_yz_0_0_0_0_xy_zz_yy = buffer_2000_sddd[933];

    auto g_yz_0_0_0_0_xy_zz_yz = buffer_2000_sddd[934];

    auto g_yz_0_0_0_0_xy_zz_zz = buffer_2000_sddd[935];

    auto g_yz_0_0_0_0_xz_xx_xx = buffer_2000_sddd[936];

    auto g_yz_0_0_0_0_xz_xx_xy = buffer_2000_sddd[937];

    auto g_yz_0_0_0_0_xz_xx_xz = buffer_2000_sddd[938];

    auto g_yz_0_0_0_0_xz_xx_yy = buffer_2000_sddd[939];

    auto g_yz_0_0_0_0_xz_xx_yz = buffer_2000_sddd[940];

    auto g_yz_0_0_0_0_xz_xx_zz = buffer_2000_sddd[941];

    auto g_yz_0_0_0_0_xz_xy_xx = buffer_2000_sddd[942];

    auto g_yz_0_0_0_0_xz_xy_xy = buffer_2000_sddd[943];

    auto g_yz_0_0_0_0_xz_xy_xz = buffer_2000_sddd[944];

    auto g_yz_0_0_0_0_xz_xy_yy = buffer_2000_sddd[945];

    auto g_yz_0_0_0_0_xz_xy_yz = buffer_2000_sddd[946];

    auto g_yz_0_0_0_0_xz_xy_zz = buffer_2000_sddd[947];

    auto g_yz_0_0_0_0_xz_xz_xx = buffer_2000_sddd[948];

    auto g_yz_0_0_0_0_xz_xz_xy = buffer_2000_sddd[949];

    auto g_yz_0_0_0_0_xz_xz_xz = buffer_2000_sddd[950];

    auto g_yz_0_0_0_0_xz_xz_yy = buffer_2000_sddd[951];

    auto g_yz_0_0_0_0_xz_xz_yz = buffer_2000_sddd[952];

    auto g_yz_0_0_0_0_xz_xz_zz = buffer_2000_sddd[953];

    auto g_yz_0_0_0_0_xz_yy_xx = buffer_2000_sddd[954];

    auto g_yz_0_0_0_0_xz_yy_xy = buffer_2000_sddd[955];

    auto g_yz_0_0_0_0_xz_yy_xz = buffer_2000_sddd[956];

    auto g_yz_0_0_0_0_xz_yy_yy = buffer_2000_sddd[957];

    auto g_yz_0_0_0_0_xz_yy_yz = buffer_2000_sddd[958];

    auto g_yz_0_0_0_0_xz_yy_zz = buffer_2000_sddd[959];

    auto g_yz_0_0_0_0_xz_yz_xx = buffer_2000_sddd[960];

    auto g_yz_0_0_0_0_xz_yz_xy = buffer_2000_sddd[961];

    auto g_yz_0_0_0_0_xz_yz_xz = buffer_2000_sddd[962];

    auto g_yz_0_0_0_0_xz_yz_yy = buffer_2000_sddd[963];

    auto g_yz_0_0_0_0_xz_yz_yz = buffer_2000_sddd[964];

    auto g_yz_0_0_0_0_xz_yz_zz = buffer_2000_sddd[965];

    auto g_yz_0_0_0_0_xz_zz_xx = buffer_2000_sddd[966];

    auto g_yz_0_0_0_0_xz_zz_xy = buffer_2000_sddd[967];

    auto g_yz_0_0_0_0_xz_zz_xz = buffer_2000_sddd[968];

    auto g_yz_0_0_0_0_xz_zz_yy = buffer_2000_sddd[969];

    auto g_yz_0_0_0_0_xz_zz_yz = buffer_2000_sddd[970];

    auto g_yz_0_0_0_0_xz_zz_zz = buffer_2000_sddd[971];

    auto g_yz_0_0_0_0_yy_xx_xx = buffer_2000_sddd[972];

    auto g_yz_0_0_0_0_yy_xx_xy = buffer_2000_sddd[973];

    auto g_yz_0_0_0_0_yy_xx_xz = buffer_2000_sddd[974];

    auto g_yz_0_0_0_0_yy_xx_yy = buffer_2000_sddd[975];

    auto g_yz_0_0_0_0_yy_xx_yz = buffer_2000_sddd[976];

    auto g_yz_0_0_0_0_yy_xx_zz = buffer_2000_sddd[977];

    auto g_yz_0_0_0_0_yy_xy_xx = buffer_2000_sddd[978];

    auto g_yz_0_0_0_0_yy_xy_xy = buffer_2000_sddd[979];

    auto g_yz_0_0_0_0_yy_xy_xz = buffer_2000_sddd[980];

    auto g_yz_0_0_0_0_yy_xy_yy = buffer_2000_sddd[981];

    auto g_yz_0_0_0_0_yy_xy_yz = buffer_2000_sddd[982];

    auto g_yz_0_0_0_0_yy_xy_zz = buffer_2000_sddd[983];

    auto g_yz_0_0_0_0_yy_xz_xx = buffer_2000_sddd[984];

    auto g_yz_0_0_0_0_yy_xz_xy = buffer_2000_sddd[985];

    auto g_yz_0_0_0_0_yy_xz_xz = buffer_2000_sddd[986];

    auto g_yz_0_0_0_0_yy_xz_yy = buffer_2000_sddd[987];

    auto g_yz_0_0_0_0_yy_xz_yz = buffer_2000_sddd[988];

    auto g_yz_0_0_0_0_yy_xz_zz = buffer_2000_sddd[989];

    auto g_yz_0_0_0_0_yy_yy_xx = buffer_2000_sddd[990];

    auto g_yz_0_0_0_0_yy_yy_xy = buffer_2000_sddd[991];

    auto g_yz_0_0_0_0_yy_yy_xz = buffer_2000_sddd[992];

    auto g_yz_0_0_0_0_yy_yy_yy = buffer_2000_sddd[993];

    auto g_yz_0_0_0_0_yy_yy_yz = buffer_2000_sddd[994];

    auto g_yz_0_0_0_0_yy_yy_zz = buffer_2000_sddd[995];

    auto g_yz_0_0_0_0_yy_yz_xx = buffer_2000_sddd[996];

    auto g_yz_0_0_0_0_yy_yz_xy = buffer_2000_sddd[997];

    auto g_yz_0_0_0_0_yy_yz_xz = buffer_2000_sddd[998];

    auto g_yz_0_0_0_0_yy_yz_yy = buffer_2000_sddd[999];

    auto g_yz_0_0_0_0_yy_yz_yz = buffer_2000_sddd[1000];

    auto g_yz_0_0_0_0_yy_yz_zz = buffer_2000_sddd[1001];

    auto g_yz_0_0_0_0_yy_zz_xx = buffer_2000_sddd[1002];

    auto g_yz_0_0_0_0_yy_zz_xy = buffer_2000_sddd[1003];

    auto g_yz_0_0_0_0_yy_zz_xz = buffer_2000_sddd[1004];

    auto g_yz_0_0_0_0_yy_zz_yy = buffer_2000_sddd[1005];

    auto g_yz_0_0_0_0_yy_zz_yz = buffer_2000_sddd[1006];

    auto g_yz_0_0_0_0_yy_zz_zz = buffer_2000_sddd[1007];

    auto g_yz_0_0_0_0_yz_xx_xx = buffer_2000_sddd[1008];

    auto g_yz_0_0_0_0_yz_xx_xy = buffer_2000_sddd[1009];

    auto g_yz_0_0_0_0_yz_xx_xz = buffer_2000_sddd[1010];

    auto g_yz_0_0_0_0_yz_xx_yy = buffer_2000_sddd[1011];

    auto g_yz_0_0_0_0_yz_xx_yz = buffer_2000_sddd[1012];

    auto g_yz_0_0_0_0_yz_xx_zz = buffer_2000_sddd[1013];

    auto g_yz_0_0_0_0_yz_xy_xx = buffer_2000_sddd[1014];

    auto g_yz_0_0_0_0_yz_xy_xy = buffer_2000_sddd[1015];

    auto g_yz_0_0_0_0_yz_xy_xz = buffer_2000_sddd[1016];

    auto g_yz_0_0_0_0_yz_xy_yy = buffer_2000_sddd[1017];

    auto g_yz_0_0_0_0_yz_xy_yz = buffer_2000_sddd[1018];

    auto g_yz_0_0_0_0_yz_xy_zz = buffer_2000_sddd[1019];

    auto g_yz_0_0_0_0_yz_xz_xx = buffer_2000_sddd[1020];

    auto g_yz_0_0_0_0_yz_xz_xy = buffer_2000_sddd[1021];

    auto g_yz_0_0_0_0_yz_xz_xz = buffer_2000_sddd[1022];

    auto g_yz_0_0_0_0_yz_xz_yy = buffer_2000_sddd[1023];

    auto g_yz_0_0_0_0_yz_xz_yz = buffer_2000_sddd[1024];

    auto g_yz_0_0_0_0_yz_xz_zz = buffer_2000_sddd[1025];

    auto g_yz_0_0_0_0_yz_yy_xx = buffer_2000_sddd[1026];

    auto g_yz_0_0_0_0_yz_yy_xy = buffer_2000_sddd[1027];

    auto g_yz_0_0_0_0_yz_yy_xz = buffer_2000_sddd[1028];

    auto g_yz_0_0_0_0_yz_yy_yy = buffer_2000_sddd[1029];

    auto g_yz_0_0_0_0_yz_yy_yz = buffer_2000_sddd[1030];

    auto g_yz_0_0_0_0_yz_yy_zz = buffer_2000_sddd[1031];

    auto g_yz_0_0_0_0_yz_yz_xx = buffer_2000_sddd[1032];

    auto g_yz_0_0_0_0_yz_yz_xy = buffer_2000_sddd[1033];

    auto g_yz_0_0_0_0_yz_yz_xz = buffer_2000_sddd[1034];

    auto g_yz_0_0_0_0_yz_yz_yy = buffer_2000_sddd[1035];

    auto g_yz_0_0_0_0_yz_yz_yz = buffer_2000_sddd[1036];

    auto g_yz_0_0_0_0_yz_yz_zz = buffer_2000_sddd[1037];

    auto g_yz_0_0_0_0_yz_zz_xx = buffer_2000_sddd[1038];

    auto g_yz_0_0_0_0_yz_zz_xy = buffer_2000_sddd[1039];

    auto g_yz_0_0_0_0_yz_zz_xz = buffer_2000_sddd[1040];

    auto g_yz_0_0_0_0_yz_zz_yy = buffer_2000_sddd[1041];

    auto g_yz_0_0_0_0_yz_zz_yz = buffer_2000_sddd[1042];

    auto g_yz_0_0_0_0_yz_zz_zz = buffer_2000_sddd[1043];

    auto g_yz_0_0_0_0_zz_xx_xx = buffer_2000_sddd[1044];

    auto g_yz_0_0_0_0_zz_xx_xy = buffer_2000_sddd[1045];

    auto g_yz_0_0_0_0_zz_xx_xz = buffer_2000_sddd[1046];

    auto g_yz_0_0_0_0_zz_xx_yy = buffer_2000_sddd[1047];

    auto g_yz_0_0_0_0_zz_xx_yz = buffer_2000_sddd[1048];

    auto g_yz_0_0_0_0_zz_xx_zz = buffer_2000_sddd[1049];

    auto g_yz_0_0_0_0_zz_xy_xx = buffer_2000_sddd[1050];

    auto g_yz_0_0_0_0_zz_xy_xy = buffer_2000_sddd[1051];

    auto g_yz_0_0_0_0_zz_xy_xz = buffer_2000_sddd[1052];

    auto g_yz_0_0_0_0_zz_xy_yy = buffer_2000_sddd[1053];

    auto g_yz_0_0_0_0_zz_xy_yz = buffer_2000_sddd[1054];

    auto g_yz_0_0_0_0_zz_xy_zz = buffer_2000_sddd[1055];

    auto g_yz_0_0_0_0_zz_xz_xx = buffer_2000_sddd[1056];

    auto g_yz_0_0_0_0_zz_xz_xy = buffer_2000_sddd[1057];

    auto g_yz_0_0_0_0_zz_xz_xz = buffer_2000_sddd[1058];

    auto g_yz_0_0_0_0_zz_xz_yy = buffer_2000_sddd[1059];

    auto g_yz_0_0_0_0_zz_xz_yz = buffer_2000_sddd[1060];

    auto g_yz_0_0_0_0_zz_xz_zz = buffer_2000_sddd[1061];

    auto g_yz_0_0_0_0_zz_yy_xx = buffer_2000_sddd[1062];

    auto g_yz_0_0_0_0_zz_yy_xy = buffer_2000_sddd[1063];

    auto g_yz_0_0_0_0_zz_yy_xz = buffer_2000_sddd[1064];

    auto g_yz_0_0_0_0_zz_yy_yy = buffer_2000_sddd[1065];

    auto g_yz_0_0_0_0_zz_yy_yz = buffer_2000_sddd[1066];

    auto g_yz_0_0_0_0_zz_yy_zz = buffer_2000_sddd[1067];

    auto g_yz_0_0_0_0_zz_yz_xx = buffer_2000_sddd[1068];

    auto g_yz_0_0_0_0_zz_yz_xy = buffer_2000_sddd[1069];

    auto g_yz_0_0_0_0_zz_yz_xz = buffer_2000_sddd[1070];

    auto g_yz_0_0_0_0_zz_yz_yy = buffer_2000_sddd[1071];

    auto g_yz_0_0_0_0_zz_yz_yz = buffer_2000_sddd[1072];

    auto g_yz_0_0_0_0_zz_yz_zz = buffer_2000_sddd[1073];

    auto g_yz_0_0_0_0_zz_zz_xx = buffer_2000_sddd[1074];

    auto g_yz_0_0_0_0_zz_zz_xy = buffer_2000_sddd[1075];

    auto g_yz_0_0_0_0_zz_zz_xz = buffer_2000_sddd[1076];

    auto g_yz_0_0_0_0_zz_zz_yy = buffer_2000_sddd[1077];

    auto g_yz_0_0_0_0_zz_zz_yz = buffer_2000_sddd[1078];

    auto g_yz_0_0_0_0_zz_zz_zz = buffer_2000_sddd[1079];

    auto g_zz_0_0_0_0_xx_xx_xx = buffer_2000_sddd[1080];

    auto g_zz_0_0_0_0_xx_xx_xy = buffer_2000_sddd[1081];

    auto g_zz_0_0_0_0_xx_xx_xz = buffer_2000_sddd[1082];

    auto g_zz_0_0_0_0_xx_xx_yy = buffer_2000_sddd[1083];

    auto g_zz_0_0_0_0_xx_xx_yz = buffer_2000_sddd[1084];

    auto g_zz_0_0_0_0_xx_xx_zz = buffer_2000_sddd[1085];

    auto g_zz_0_0_0_0_xx_xy_xx = buffer_2000_sddd[1086];

    auto g_zz_0_0_0_0_xx_xy_xy = buffer_2000_sddd[1087];

    auto g_zz_0_0_0_0_xx_xy_xz = buffer_2000_sddd[1088];

    auto g_zz_0_0_0_0_xx_xy_yy = buffer_2000_sddd[1089];

    auto g_zz_0_0_0_0_xx_xy_yz = buffer_2000_sddd[1090];

    auto g_zz_0_0_0_0_xx_xy_zz = buffer_2000_sddd[1091];

    auto g_zz_0_0_0_0_xx_xz_xx = buffer_2000_sddd[1092];

    auto g_zz_0_0_0_0_xx_xz_xy = buffer_2000_sddd[1093];

    auto g_zz_0_0_0_0_xx_xz_xz = buffer_2000_sddd[1094];

    auto g_zz_0_0_0_0_xx_xz_yy = buffer_2000_sddd[1095];

    auto g_zz_0_0_0_0_xx_xz_yz = buffer_2000_sddd[1096];

    auto g_zz_0_0_0_0_xx_xz_zz = buffer_2000_sddd[1097];

    auto g_zz_0_0_0_0_xx_yy_xx = buffer_2000_sddd[1098];

    auto g_zz_0_0_0_0_xx_yy_xy = buffer_2000_sddd[1099];

    auto g_zz_0_0_0_0_xx_yy_xz = buffer_2000_sddd[1100];

    auto g_zz_0_0_0_0_xx_yy_yy = buffer_2000_sddd[1101];

    auto g_zz_0_0_0_0_xx_yy_yz = buffer_2000_sddd[1102];

    auto g_zz_0_0_0_0_xx_yy_zz = buffer_2000_sddd[1103];

    auto g_zz_0_0_0_0_xx_yz_xx = buffer_2000_sddd[1104];

    auto g_zz_0_0_0_0_xx_yz_xy = buffer_2000_sddd[1105];

    auto g_zz_0_0_0_0_xx_yz_xz = buffer_2000_sddd[1106];

    auto g_zz_0_0_0_0_xx_yz_yy = buffer_2000_sddd[1107];

    auto g_zz_0_0_0_0_xx_yz_yz = buffer_2000_sddd[1108];

    auto g_zz_0_0_0_0_xx_yz_zz = buffer_2000_sddd[1109];

    auto g_zz_0_0_0_0_xx_zz_xx = buffer_2000_sddd[1110];

    auto g_zz_0_0_0_0_xx_zz_xy = buffer_2000_sddd[1111];

    auto g_zz_0_0_0_0_xx_zz_xz = buffer_2000_sddd[1112];

    auto g_zz_0_0_0_0_xx_zz_yy = buffer_2000_sddd[1113];

    auto g_zz_0_0_0_0_xx_zz_yz = buffer_2000_sddd[1114];

    auto g_zz_0_0_0_0_xx_zz_zz = buffer_2000_sddd[1115];

    auto g_zz_0_0_0_0_xy_xx_xx = buffer_2000_sddd[1116];

    auto g_zz_0_0_0_0_xy_xx_xy = buffer_2000_sddd[1117];

    auto g_zz_0_0_0_0_xy_xx_xz = buffer_2000_sddd[1118];

    auto g_zz_0_0_0_0_xy_xx_yy = buffer_2000_sddd[1119];

    auto g_zz_0_0_0_0_xy_xx_yz = buffer_2000_sddd[1120];

    auto g_zz_0_0_0_0_xy_xx_zz = buffer_2000_sddd[1121];

    auto g_zz_0_0_0_0_xy_xy_xx = buffer_2000_sddd[1122];

    auto g_zz_0_0_0_0_xy_xy_xy = buffer_2000_sddd[1123];

    auto g_zz_0_0_0_0_xy_xy_xz = buffer_2000_sddd[1124];

    auto g_zz_0_0_0_0_xy_xy_yy = buffer_2000_sddd[1125];

    auto g_zz_0_0_0_0_xy_xy_yz = buffer_2000_sddd[1126];

    auto g_zz_0_0_0_0_xy_xy_zz = buffer_2000_sddd[1127];

    auto g_zz_0_0_0_0_xy_xz_xx = buffer_2000_sddd[1128];

    auto g_zz_0_0_0_0_xy_xz_xy = buffer_2000_sddd[1129];

    auto g_zz_0_0_0_0_xy_xz_xz = buffer_2000_sddd[1130];

    auto g_zz_0_0_0_0_xy_xz_yy = buffer_2000_sddd[1131];

    auto g_zz_0_0_0_0_xy_xz_yz = buffer_2000_sddd[1132];

    auto g_zz_0_0_0_0_xy_xz_zz = buffer_2000_sddd[1133];

    auto g_zz_0_0_0_0_xy_yy_xx = buffer_2000_sddd[1134];

    auto g_zz_0_0_0_0_xy_yy_xy = buffer_2000_sddd[1135];

    auto g_zz_0_0_0_0_xy_yy_xz = buffer_2000_sddd[1136];

    auto g_zz_0_0_0_0_xy_yy_yy = buffer_2000_sddd[1137];

    auto g_zz_0_0_0_0_xy_yy_yz = buffer_2000_sddd[1138];

    auto g_zz_0_0_0_0_xy_yy_zz = buffer_2000_sddd[1139];

    auto g_zz_0_0_0_0_xy_yz_xx = buffer_2000_sddd[1140];

    auto g_zz_0_0_0_0_xy_yz_xy = buffer_2000_sddd[1141];

    auto g_zz_0_0_0_0_xy_yz_xz = buffer_2000_sddd[1142];

    auto g_zz_0_0_0_0_xy_yz_yy = buffer_2000_sddd[1143];

    auto g_zz_0_0_0_0_xy_yz_yz = buffer_2000_sddd[1144];

    auto g_zz_0_0_0_0_xy_yz_zz = buffer_2000_sddd[1145];

    auto g_zz_0_0_0_0_xy_zz_xx = buffer_2000_sddd[1146];

    auto g_zz_0_0_0_0_xy_zz_xy = buffer_2000_sddd[1147];

    auto g_zz_0_0_0_0_xy_zz_xz = buffer_2000_sddd[1148];

    auto g_zz_0_0_0_0_xy_zz_yy = buffer_2000_sddd[1149];

    auto g_zz_0_0_0_0_xy_zz_yz = buffer_2000_sddd[1150];

    auto g_zz_0_0_0_0_xy_zz_zz = buffer_2000_sddd[1151];

    auto g_zz_0_0_0_0_xz_xx_xx = buffer_2000_sddd[1152];

    auto g_zz_0_0_0_0_xz_xx_xy = buffer_2000_sddd[1153];

    auto g_zz_0_0_0_0_xz_xx_xz = buffer_2000_sddd[1154];

    auto g_zz_0_0_0_0_xz_xx_yy = buffer_2000_sddd[1155];

    auto g_zz_0_0_0_0_xz_xx_yz = buffer_2000_sddd[1156];

    auto g_zz_0_0_0_0_xz_xx_zz = buffer_2000_sddd[1157];

    auto g_zz_0_0_0_0_xz_xy_xx = buffer_2000_sddd[1158];

    auto g_zz_0_0_0_0_xz_xy_xy = buffer_2000_sddd[1159];

    auto g_zz_0_0_0_0_xz_xy_xz = buffer_2000_sddd[1160];

    auto g_zz_0_0_0_0_xz_xy_yy = buffer_2000_sddd[1161];

    auto g_zz_0_0_0_0_xz_xy_yz = buffer_2000_sddd[1162];

    auto g_zz_0_0_0_0_xz_xy_zz = buffer_2000_sddd[1163];

    auto g_zz_0_0_0_0_xz_xz_xx = buffer_2000_sddd[1164];

    auto g_zz_0_0_0_0_xz_xz_xy = buffer_2000_sddd[1165];

    auto g_zz_0_0_0_0_xz_xz_xz = buffer_2000_sddd[1166];

    auto g_zz_0_0_0_0_xz_xz_yy = buffer_2000_sddd[1167];

    auto g_zz_0_0_0_0_xz_xz_yz = buffer_2000_sddd[1168];

    auto g_zz_0_0_0_0_xz_xz_zz = buffer_2000_sddd[1169];

    auto g_zz_0_0_0_0_xz_yy_xx = buffer_2000_sddd[1170];

    auto g_zz_0_0_0_0_xz_yy_xy = buffer_2000_sddd[1171];

    auto g_zz_0_0_0_0_xz_yy_xz = buffer_2000_sddd[1172];

    auto g_zz_0_0_0_0_xz_yy_yy = buffer_2000_sddd[1173];

    auto g_zz_0_0_0_0_xz_yy_yz = buffer_2000_sddd[1174];

    auto g_zz_0_0_0_0_xz_yy_zz = buffer_2000_sddd[1175];

    auto g_zz_0_0_0_0_xz_yz_xx = buffer_2000_sddd[1176];

    auto g_zz_0_0_0_0_xz_yz_xy = buffer_2000_sddd[1177];

    auto g_zz_0_0_0_0_xz_yz_xz = buffer_2000_sddd[1178];

    auto g_zz_0_0_0_0_xz_yz_yy = buffer_2000_sddd[1179];

    auto g_zz_0_0_0_0_xz_yz_yz = buffer_2000_sddd[1180];

    auto g_zz_0_0_0_0_xz_yz_zz = buffer_2000_sddd[1181];

    auto g_zz_0_0_0_0_xz_zz_xx = buffer_2000_sddd[1182];

    auto g_zz_0_0_0_0_xz_zz_xy = buffer_2000_sddd[1183];

    auto g_zz_0_0_0_0_xz_zz_xz = buffer_2000_sddd[1184];

    auto g_zz_0_0_0_0_xz_zz_yy = buffer_2000_sddd[1185];

    auto g_zz_0_0_0_0_xz_zz_yz = buffer_2000_sddd[1186];

    auto g_zz_0_0_0_0_xz_zz_zz = buffer_2000_sddd[1187];

    auto g_zz_0_0_0_0_yy_xx_xx = buffer_2000_sddd[1188];

    auto g_zz_0_0_0_0_yy_xx_xy = buffer_2000_sddd[1189];

    auto g_zz_0_0_0_0_yy_xx_xz = buffer_2000_sddd[1190];

    auto g_zz_0_0_0_0_yy_xx_yy = buffer_2000_sddd[1191];

    auto g_zz_0_0_0_0_yy_xx_yz = buffer_2000_sddd[1192];

    auto g_zz_0_0_0_0_yy_xx_zz = buffer_2000_sddd[1193];

    auto g_zz_0_0_0_0_yy_xy_xx = buffer_2000_sddd[1194];

    auto g_zz_0_0_0_0_yy_xy_xy = buffer_2000_sddd[1195];

    auto g_zz_0_0_0_0_yy_xy_xz = buffer_2000_sddd[1196];

    auto g_zz_0_0_0_0_yy_xy_yy = buffer_2000_sddd[1197];

    auto g_zz_0_0_0_0_yy_xy_yz = buffer_2000_sddd[1198];

    auto g_zz_0_0_0_0_yy_xy_zz = buffer_2000_sddd[1199];

    auto g_zz_0_0_0_0_yy_xz_xx = buffer_2000_sddd[1200];

    auto g_zz_0_0_0_0_yy_xz_xy = buffer_2000_sddd[1201];

    auto g_zz_0_0_0_0_yy_xz_xz = buffer_2000_sddd[1202];

    auto g_zz_0_0_0_0_yy_xz_yy = buffer_2000_sddd[1203];

    auto g_zz_0_0_0_0_yy_xz_yz = buffer_2000_sddd[1204];

    auto g_zz_0_0_0_0_yy_xz_zz = buffer_2000_sddd[1205];

    auto g_zz_0_0_0_0_yy_yy_xx = buffer_2000_sddd[1206];

    auto g_zz_0_0_0_0_yy_yy_xy = buffer_2000_sddd[1207];

    auto g_zz_0_0_0_0_yy_yy_xz = buffer_2000_sddd[1208];

    auto g_zz_0_0_0_0_yy_yy_yy = buffer_2000_sddd[1209];

    auto g_zz_0_0_0_0_yy_yy_yz = buffer_2000_sddd[1210];

    auto g_zz_0_0_0_0_yy_yy_zz = buffer_2000_sddd[1211];

    auto g_zz_0_0_0_0_yy_yz_xx = buffer_2000_sddd[1212];

    auto g_zz_0_0_0_0_yy_yz_xy = buffer_2000_sddd[1213];

    auto g_zz_0_0_0_0_yy_yz_xz = buffer_2000_sddd[1214];

    auto g_zz_0_0_0_0_yy_yz_yy = buffer_2000_sddd[1215];

    auto g_zz_0_0_0_0_yy_yz_yz = buffer_2000_sddd[1216];

    auto g_zz_0_0_0_0_yy_yz_zz = buffer_2000_sddd[1217];

    auto g_zz_0_0_0_0_yy_zz_xx = buffer_2000_sddd[1218];

    auto g_zz_0_0_0_0_yy_zz_xy = buffer_2000_sddd[1219];

    auto g_zz_0_0_0_0_yy_zz_xz = buffer_2000_sddd[1220];

    auto g_zz_0_0_0_0_yy_zz_yy = buffer_2000_sddd[1221];

    auto g_zz_0_0_0_0_yy_zz_yz = buffer_2000_sddd[1222];

    auto g_zz_0_0_0_0_yy_zz_zz = buffer_2000_sddd[1223];

    auto g_zz_0_0_0_0_yz_xx_xx = buffer_2000_sddd[1224];

    auto g_zz_0_0_0_0_yz_xx_xy = buffer_2000_sddd[1225];

    auto g_zz_0_0_0_0_yz_xx_xz = buffer_2000_sddd[1226];

    auto g_zz_0_0_0_0_yz_xx_yy = buffer_2000_sddd[1227];

    auto g_zz_0_0_0_0_yz_xx_yz = buffer_2000_sddd[1228];

    auto g_zz_0_0_0_0_yz_xx_zz = buffer_2000_sddd[1229];

    auto g_zz_0_0_0_0_yz_xy_xx = buffer_2000_sddd[1230];

    auto g_zz_0_0_0_0_yz_xy_xy = buffer_2000_sddd[1231];

    auto g_zz_0_0_0_0_yz_xy_xz = buffer_2000_sddd[1232];

    auto g_zz_0_0_0_0_yz_xy_yy = buffer_2000_sddd[1233];

    auto g_zz_0_0_0_0_yz_xy_yz = buffer_2000_sddd[1234];

    auto g_zz_0_0_0_0_yz_xy_zz = buffer_2000_sddd[1235];

    auto g_zz_0_0_0_0_yz_xz_xx = buffer_2000_sddd[1236];

    auto g_zz_0_0_0_0_yz_xz_xy = buffer_2000_sddd[1237];

    auto g_zz_0_0_0_0_yz_xz_xz = buffer_2000_sddd[1238];

    auto g_zz_0_0_0_0_yz_xz_yy = buffer_2000_sddd[1239];

    auto g_zz_0_0_0_0_yz_xz_yz = buffer_2000_sddd[1240];

    auto g_zz_0_0_0_0_yz_xz_zz = buffer_2000_sddd[1241];

    auto g_zz_0_0_0_0_yz_yy_xx = buffer_2000_sddd[1242];

    auto g_zz_0_0_0_0_yz_yy_xy = buffer_2000_sddd[1243];

    auto g_zz_0_0_0_0_yz_yy_xz = buffer_2000_sddd[1244];

    auto g_zz_0_0_0_0_yz_yy_yy = buffer_2000_sddd[1245];

    auto g_zz_0_0_0_0_yz_yy_yz = buffer_2000_sddd[1246];

    auto g_zz_0_0_0_0_yz_yy_zz = buffer_2000_sddd[1247];

    auto g_zz_0_0_0_0_yz_yz_xx = buffer_2000_sddd[1248];

    auto g_zz_0_0_0_0_yz_yz_xy = buffer_2000_sddd[1249];

    auto g_zz_0_0_0_0_yz_yz_xz = buffer_2000_sddd[1250];

    auto g_zz_0_0_0_0_yz_yz_yy = buffer_2000_sddd[1251];

    auto g_zz_0_0_0_0_yz_yz_yz = buffer_2000_sddd[1252];

    auto g_zz_0_0_0_0_yz_yz_zz = buffer_2000_sddd[1253];

    auto g_zz_0_0_0_0_yz_zz_xx = buffer_2000_sddd[1254];

    auto g_zz_0_0_0_0_yz_zz_xy = buffer_2000_sddd[1255];

    auto g_zz_0_0_0_0_yz_zz_xz = buffer_2000_sddd[1256];

    auto g_zz_0_0_0_0_yz_zz_yy = buffer_2000_sddd[1257];

    auto g_zz_0_0_0_0_yz_zz_yz = buffer_2000_sddd[1258];

    auto g_zz_0_0_0_0_yz_zz_zz = buffer_2000_sddd[1259];

    auto g_zz_0_0_0_0_zz_xx_xx = buffer_2000_sddd[1260];

    auto g_zz_0_0_0_0_zz_xx_xy = buffer_2000_sddd[1261];

    auto g_zz_0_0_0_0_zz_xx_xz = buffer_2000_sddd[1262];

    auto g_zz_0_0_0_0_zz_xx_yy = buffer_2000_sddd[1263];

    auto g_zz_0_0_0_0_zz_xx_yz = buffer_2000_sddd[1264];

    auto g_zz_0_0_0_0_zz_xx_zz = buffer_2000_sddd[1265];

    auto g_zz_0_0_0_0_zz_xy_xx = buffer_2000_sddd[1266];

    auto g_zz_0_0_0_0_zz_xy_xy = buffer_2000_sddd[1267];

    auto g_zz_0_0_0_0_zz_xy_xz = buffer_2000_sddd[1268];

    auto g_zz_0_0_0_0_zz_xy_yy = buffer_2000_sddd[1269];

    auto g_zz_0_0_0_0_zz_xy_yz = buffer_2000_sddd[1270];

    auto g_zz_0_0_0_0_zz_xy_zz = buffer_2000_sddd[1271];

    auto g_zz_0_0_0_0_zz_xz_xx = buffer_2000_sddd[1272];

    auto g_zz_0_0_0_0_zz_xz_xy = buffer_2000_sddd[1273];

    auto g_zz_0_0_0_0_zz_xz_xz = buffer_2000_sddd[1274];

    auto g_zz_0_0_0_0_zz_xz_yy = buffer_2000_sddd[1275];

    auto g_zz_0_0_0_0_zz_xz_yz = buffer_2000_sddd[1276];

    auto g_zz_0_0_0_0_zz_xz_zz = buffer_2000_sddd[1277];

    auto g_zz_0_0_0_0_zz_yy_xx = buffer_2000_sddd[1278];

    auto g_zz_0_0_0_0_zz_yy_xy = buffer_2000_sddd[1279];

    auto g_zz_0_0_0_0_zz_yy_xz = buffer_2000_sddd[1280];

    auto g_zz_0_0_0_0_zz_yy_yy = buffer_2000_sddd[1281];

    auto g_zz_0_0_0_0_zz_yy_yz = buffer_2000_sddd[1282];

    auto g_zz_0_0_0_0_zz_yy_zz = buffer_2000_sddd[1283];

    auto g_zz_0_0_0_0_zz_yz_xx = buffer_2000_sddd[1284];

    auto g_zz_0_0_0_0_zz_yz_xy = buffer_2000_sddd[1285];

    auto g_zz_0_0_0_0_zz_yz_xz = buffer_2000_sddd[1286];

    auto g_zz_0_0_0_0_zz_yz_yy = buffer_2000_sddd[1287];

    auto g_zz_0_0_0_0_zz_yz_yz = buffer_2000_sddd[1288];

    auto g_zz_0_0_0_0_zz_yz_zz = buffer_2000_sddd[1289];

    auto g_zz_0_0_0_0_zz_zz_xx = buffer_2000_sddd[1290];

    auto g_zz_0_0_0_0_zz_zz_xy = buffer_2000_sddd[1291];

    auto g_zz_0_0_0_0_zz_zz_xz = buffer_2000_sddd[1292];

    auto g_zz_0_0_0_0_zz_zz_yy = buffer_2000_sddd[1293];

    auto g_zz_0_0_0_0_zz_zz_yz = buffer_2000_sddd[1294];

    auto g_zz_0_0_0_0_zz_zz_zz = buffer_2000_sddd[1295];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_xx_xx_xx, g_0_xx_xx_xy, g_0_xx_xx_xz, g_0_xx_xx_yy, g_0_xx_xx_yz, g_0_xx_xx_zz, g_xx_0_0_0_0_xx_xx_xx, g_xx_0_0_0_0_xx_xx_xy, g_xx_0_0_0_0_xx_xx_xz, g_xx_0_0_0_0_xx_xx_yy, g_xx_0_0_0_0_xx_xx_yz, g_xx_0_0_0_0_xx_xx_zz, g_xx_xx_xx_xx, g_xx_xx_xx_xy, g_xx_xx_xx_xz, g_xx_xx_xx_yy, g_xx_xx_xx_yz, g_xx_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xx_xx_xx[i] = -2.0 * g_0_xx_xx_xx[i] * a_exp + 4.0 * g_xx_xx_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_xx_xy[i] = -2.0 * g_0_xx_xx_xy[i] * a_exp + 4.0 * g_xx_xx_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_xx_xz[i] = -2.0 * g_0_xx_xx_xz[i] * a_exp + 4.0 * g_xx_xx_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_xx_yy[i] = -2.0 * g_0_xx_xx_yy[i] * a_exp + 4.0 * g_xx_xx_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_xx_yz[i] = -2.0 * g_0_xx_xx_yz[i] * a_exp + 4.0 * g_xx_xx_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_xx_zz[i] = -2.0 * g_0_xx_xx_zz[i] * a_exp + 4.0 * g_xx_xx_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_xx_xy_xx, g_0_xx_xy_xy, g_0_xx_xy_xz, g_0_xx_xy_yy, g_0_xx_xy_yz, g_0_xx_xy_zz, g_xx_0_0_0_0_xx_xy_xx, g_xx_0_0_0_0_xx_xy_xy, g_xx_0_0_0_0_xx_xy_xz, g_xx_0_0_0_0_xx_xy_yy, g_xx_0_0_0_0_xx_xy_yz, g_xx_0_0_0_0_xx_xy_zz, g_xx_xx_xy_xx, g_xx_xx_xy_xy, g_xx_xx_xy_xz, g_xx_xx_xy_yy, g_xx_xx_xy_yz, g_xx_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xx_xy_xx[i] = -2.0 * g_0_xx_xy_xx[i] * a_exp + 4.0 * g_xx_xx_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_xy_xy[i] = -2.0 * g_0_xx_xy_xy[i] * a_exp + 4.0 * g_xx_xx_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_xy_xz[i] = -2.0 * g_0_xx_xy_xz[i] * a_exp + 4.0 * g_xx_xx_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_xy_yy[i] = -2.0 * g_0_xx_xy_yy[i] * a_exp + 4.0 * g_xx_xx_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_xy_yz[i] = -2.0 * g_0_xx_xy_yz[i] * a_exp + 4.0 * g_xx_xx_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_xy_zz[i] = -2.0 * g_0_xx_xy_zz[i] * a_exp + 4.0 * g_xx_xx_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_xx_xz_xx, g_0_xx_xz_xy, g_0_xx_xz_xz, g_0_xx_xz_yy, g_0_xx_xz_yz, g_0_xx_xz_zz, g_xx_0_0_0_0_xx_xz_xx, g_xx_0_0_0_0_xx_xz_xy, g_xx_0_0_0_0_xx_xz_xz, g_xx_0_0_0_0_xx_xz_yy, g_xx_0_0_0_0_xx_xz_yz, g_xx_0_0_0_0_xx_xz_zz, g_xx_xx_xz_xx, g_xx_xx_xz_xy, g_xx_xx_xz_xz, g_xx_xx_xz_yy, g_xx_xx_xz_yz, g_xx_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xx_xz_xx[i] = -2.0 * g_0_xx_xz_xx[i] * a_exp + 4.0 * g_xx_xx_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_xz_xy[i] = -2.0 * g_0_xx_xz_xy[i] * a_exp + 4.0 * g_xx_xx_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_xz_xz[i] = -2.0 * g_0_xx_xz_xz[i] * a_exp + 4.0 * g_xx_xx_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_xz_yy[i] = -2.0 * g_0_xx_xz_yy[i] * a_exp + 4.0 * g_xx_xx_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_xz_yz[i] = -2.0 * g_0_xx_xz_yz[i] * a_exp + 4.0 * g_xx_xx_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_xz_zz[i] = -2.0 * g_0_xx_xz_zz[i] * a_exp + 4.0 * g_xx_xx_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_xx_yy_xx, g_0_xx_yy_xy, g_0_xx_yy_xz, g_0_xx_yy_yy, g_0_xx_yy_yz, g_0_xx_yy_zz, g_xx_0_0_0_0_xx_yy_xx, g_xx_0_0_0_0_xx_yy_xy, g_xx_0_0_0_0_xx_yy_xz, g_xx_0_0_0_0_xx_yy_yy, g_xx_0_0_0_0_xx_yy_yz, g_xx_0_0_0_0_xx_yy_zz, g_xx_xx_yy_xx, g_xx_xx_yy_xy, g_xx_xx_yy_xz, g_xx_xx_yy_yy, g_xx_xx_yy_yz, g_xx_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xx_yy_xx[i] = -2.0 * g_0_xx_yy_xx[i] * a_exp + 4.0 * g_xx_xx_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_yy_xy[i] = -2.0 * g_0_xx_yy_xy[i] * a_exp + 4.0 * g_xx_xx_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_yy_xz[i] = -2.0 * g_0_xx_yy_xz[i] * a_exp + 4.0 * g_xx_xx_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_yy_yy[i] = -2.0 * g_0_xx_yy_yy[i] * a_exp + 4.0 * g_xx_xx_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_yy_yz[i] = -2.0 * g_0_xx_yy_yz[i] * a_exp + 4.0 * g_xx_xx_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_yy_zz[i] = -2.0 * g_0_xx_yy_zz[i] * a_exp + 4.0 * g_xx_xx_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_xx_yz_xx, g_0_xx_yz_xy, g_0_xx_yz_xz, g_0_xx_yz_yy, g_0_xx_yz_yz, g_0_xx_yz_zz, g_xx_0_0_0_0_xx_yz_xx, g_xx_0_0_0_0_xx_yz_xy, g_xx_0_0_0_0_xx_yz_xz, g_xx_0_0_0_0_xx_yz_yy, g_xx_0_0_0_0_xx_yz_yz, g_xx_0_0_0_0_xx_yz_zz, g_xx_xx_yz_xx, g_xx_xx_yz_xy, g_xx_xx_yz_xz, g_xx_xx_yz_yy, g_xx_xx_yz_yz, g_xx_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xx_yz_xx[i] = -2.0 * g_0_xx_yz_xx[i] * a_exp + 4.0 * g_xx_xx_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_yz_xy[i] = -2.0 * g_0_xx_yz_xy[i] * a_exp + 4.0 * g_xx_xx_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_yz_xz[i] = -2.0 * g_0_xx_yz_xz[i] * a_exp + 4.0 * g_xx_xx_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_yz_yy[i] = -2.0 * g_0_xx_yz_yy[i] * a_exp + 4.0 * g_xx_xx_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_yz_yz[i] = -2.0 * g_0_xx_yz_yz[i] * a_exp + 4.0 * g_xx_xx_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_yz_zz[i] = -2.0 * g_0_xx_yz_zz[i] * a_exp + 4.0 * g_xx_xx_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_xx_zz_xx, g_0_xx_zz_xy, g_0_xx_zz_xz, g_0_xx_zz_yy, g_0_xx_zz_yz, g_0_xx_zz_zz, g_xx_0_0_0_0_xx_zz_xx, g_xx_0_0_0_0_xx_zz_xy, g_xx_0_0_0_0_xx_zz_xz, g_xx_0_0_0_0_xx_zz_yy, g_xx_0_0_0_0_xx_zz_yz, g_xx_0_0_0_0_xx_zz_zz, g_xx_xx_zz_xx, g_xx_xx_zz_xy, g_xx_xx_zz_xz, g_xx_xx_zz_yy, g_xx_xx_zz_yz, g_xx_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xx_zz_xx[i] = -2.0 * g_0_xx_zz_xx[i] * a_exp + 4.0 * g_xx_xx_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_zz_xy[i] = -2.0 * g_0_xx_zz_xy[i] * a_exp + 4.0 * g_xx_xx_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_zz_xz[i] = -2.0 * g_0_xx_zz_xz[i] * a_exp + 4.0 * g_xx_xx_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_zz_yy[i] = -2.0 * g_0_xx_zz_yy[i] * a_exp + 4.0 * g_xx_xx_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_zz_yz[i] = -2.0 * g_0_xx_zz_yz[i] * a_exp + 4.0 * g_xx_xx_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_zz_zz[i] = -2.0 * g_0_xx_zz_zz[i] * a_exp + 4.0 * g_xx_xx_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_0_xy_xx_xx, g_0_xy_xx_xy, g_0_xy_xx_xz, g_0_xy_xx_yy, g_0_xy_xx_yz, g_0_xy_xx_zz, g_xx_0_0_0_0_xy_xx_xx, g_xx_0_0_0_0_xy_xx_xy, g_xx_0_0_0_0_xy_xx_xz, g_xx_0_0_0_0_xy_xx_yy, g_xx_0_0_0_0_xy_xx_yz, g_xx_0_0_0_0_xy_xx_zz, g_xx_xy_xx_xx, g_xx_xy_xx_xy, g_xx_xy_xx_xz, g_xx_xy_xx_yy, g_xx_xy_xx_yz, g_xx_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xy_xx_xx[i] = -2.0 * g_0_xy_xx_xx[i] * a_exp + 4.0 * g_xx_xy_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_xx_xy[i] = -2.0 * g_0_xy_xx_xy[i] * a_exp + 4.0 * g_xx_xy_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_xx_xz[i] = -2.0 * g_0_xy_xx_xz[i] * a_exp + 4.0 * g_xx_xy_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_xx_yy[i] = -2.0 * g_0_xy_xx_yy[i] * a_exp + 4.0 * g_xx_xy_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_xx_yz[i] = -2.0 * g_0_xy_xx_yz[i] * a_exp + 4.0 * g_xx_xy_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_xx_zz[i] = -2.0 * g_0_xy_xx_zz[i] * a_exp + 4.0 * g_xx_xy_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_0_xy_xy_xx, g_0_xy_xy_xy, g_0_xy_xy_xz, g_0_xy_xy_yy, g_0_xy_xy_yz, g_0_xy_xy_zz, g_xx_0_0_0_0_xy_xy_xx, g_xx_0_0_0_0_xy_xy_xy, g_xx_0_0_0_0_xy_xy_xz, g_xx_0_0_0_0_xy_xy_yy, g_xx_0_0_0_0_xy_xy_yz, g_xx_0_0_0_0_xy_xy_zz, g_xx_xy_xy_xx, g_xx_xy_xy_xy, g_xx_xy_xy_xz, g_xx_xy_xy_yy, g_xx_xy_xy_yz, g_xx_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xy_xy_xx[i] = -2.0 * g_0_xy_xy_xx[i] * a_exp + 4.0 * g_xx_xy_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_xy_xy[i] = -2.0 * g_0_xy_xy_xy[i] * a_exp + 4.0 * g_xx_xy_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_xy_xz[i] = -2.0 * g_0_xy_xy_xz[i] * a_exp + 4.0 * g_xx_xy_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_xy_yy[i] = -2.0 * g_0_xy_xy_yy[i] * a_exp + 4.0 * g_xx_xy_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_xy_yz[i] = -2.0 * g_0_xy_xy_yz[i] * a_exp + 4.0 * g_xx_xy_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_xy_zz[i] = -2.0 * g_0_xy_xy_zz[i] * a_exp + 4.0 * g_xx_xy_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_0_xy_xz_xx, g_0_xy_xz_xy, g_0_xy_xz_xz, g_0_xy_xz_yy, g_0_xy_xz_yz, g_0_xy_xz_zz, g_xx_0_0_0_0_xy_xz_xx, g_xx_0_0_0_0_xy_xz_xy, g_xx_0_0_0_0_xy_xz_xz, g_xx_0_0_0_0_xy_xz_yy, g_xx_0_0_0_0_xy_xz_yz, g_xx_0_0_0_0_xy_xz_zz, g_xx_xy_xz_xx, g_xx_xy_xz_xy, g_xx_xy_xz_xz, g_xx_xy_xz_yy, g_xx_xy_xz_yz, g_xx_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xy_xz_xx[i] = -2.0 * g_0_xy_xz_xx[i] * a_exp + 4.0 * g_xx_xy_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_xz_xy[i] = -2.0 * g_0_xy_xz_xy[i] * a_exp + 4.0 * g_xx_xy_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_xz_xz[i] = -2.0 * g_0_xy_xz_xz[i] * a_exp + 4.0 * g_xx_xy_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_xz_yy[i] = -2.0 * g_0_xy_xz_yy[i] * a_exp + 4.0 * g_xx_xy_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_xz_yz[i] = -2.0 * g_0_xy_xz_yz[i] * a_exp + 4.0 * g_xx_xy_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_xz_zz[i] = -2.0 * g_0_xy_xz_zz[i] * a_exp + 4.0 * g_xx_xy_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_0_xy_yy_xx, g_0_xy_yy_xy, g_0_xy_yy_xz, g_0_xy_yy_yy, g_0_xy_yy_yz, g_0_xy_yy_zz, g_xx_0_0_0_0_xy_yy_xx, g_xx_0_0_0_0_xy_yy_xy, g_xx_0_0_0_0_xy_yy_xz, g_xx_0_0_0_0_xy_yy_yy, g_xx_0_0_0_0_xy_yy_yz, g_xx_0_0_0_0_xy_yy_zz, g_xx_xy_yy_xx, g_xx_xy_yy_xy, g_xx_xy_yy_xz, g_xx_xy_yy_yy, g_xx_xy_yy_yz, g_xx_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xy_yy_xx[i] = -2.0 * g_0_xy_yy_xx[i] * a_exp + 4.0 * g_xx_xy_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_yy_xy[i] = -2.0 * g_0_xy_yy_xy[i] * a_exp + 4.0 * g_xx_xy_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_yy_xz[i] = -2.0 * g_0_xy_yy_xz[i] * a_exp + 4.0 * g_xx_xy_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_yy_yy[i] = -2.0 * g_0_xy_yy_yy[i] * a_exp + 4.0 * g_xx_xy_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_yy_yz[i] = -2.0 * g_0_xy_yy_yz[i] * a_exp + 4.0 * g_xx_xy_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_yy_zz[i] = -2.0 * g_0_xy_yy_zz[i] * a_exp + 4.0 * g_xx_xy_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_0_xy_yz_xx, g_0_xy_yz_xy, g_0_xy_yz_xz, g_0_xy_yz_yy, g_0_xy_yz_yz, g_0_xy_yz_zz, g_xx_0_0_0_0_xy_yz_xx, g_xx_0_0_0_0_xy_yz_xy, g_xx_0_0_0_0_xy_yz_xz, g_xx_0_0_0_0_xy_yz_yy, g_xx_0_0_0_0_xy_yz_yz, g_xx_0_0_0_0_xy_yz_zz, g_xx_xy_yz_xx, g_xx_xy_yz_xy, g_xx_xy_yz_xz, g_xx_xy_yz_yy, g_xx_xy_yz_yz, g_xx_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xy_yz_xx[i] = -2.0 * g_0_xy_yz_xx[i] * a_exp + 4.0 * g_xx_xy_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_yz_xy[i] = -2.0 * g_0_xy_yz_xy[i] * a_exp + 4.0 * g_xx_xy_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_yz_xz[i] = -2.0 * g_0_xy_yz_xz[i] * a_exp + 4.0 * g_xx_xy_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_yz_yy[i] = -2.0 * g_0_xy_yz_yy[i] * a_exp + 4.0 * g_xx_xy_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_yz_yz[i] = -2.0 * g_0_xy_yz_yz[i] * a_exp + 4.0 * g_xx_xy_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_yz_zz[i] = -2.0 * g_0_xy_yz_zz[i] * a_exp + 4.0 * g_xx_xy_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_0_xy_zz_xx, g_0_xy_zz_xy, g_0_xy_zz_xz, g_0_xy_zz_yy, g_0_xy_zz_yz, g_0_xy_zz_zz, g_xx_0_0_0_0_xy_zz_xx, g_xx_0_0_0_0_xy_zz_xy, g_xx_0_0_0_0_xy_zz_xz, g_xx_0_0_0_0_xy_zz_yy, g_xx_0_0_0_0_xy_zz_yz, g_xx_0_0_0_0_xy_zz_zz, g_xx_xy_zz_xx, g_xx_xy_zz_xy, g_xx_xy_zz_xz, g_xx_xy_zz_yy, g_xx_xy_zz_yz, g_xx_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xy_zz_xx[i] = -2.0 * g_0_xy_zz_xx[i] * a_exp + 4.0 * g_xx_xy_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_zz_xy[i] = -2.0 * g_0_xy_zz_xy[i] * a_exp + 4.0 * g_xx_xy_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_zz_xz[i] = -2.0 * g_0_xy_zz_xz[i] * a_exp + 4.0 * g_xx_xy_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_zz_yy[i] = -2.0 * g_0_xy_zz_yy[i] * a_exp + 4.0 * g_xx_xy_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_zz_yz[i] = -2.0 * g_0_xy_zz_yz[i] * a_exp + 4.0 * g_xx_xy_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_zz_zz[i] = -2.0 * g_0_xy_zz_zz[i] * a_exp + 4.0 * g_xx_xy_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_0_xz_xx_xx, g_0_xz_xx_xy, g_0_xz_xx_xz, g_0_xz_xx_yy, g_0_xz_xx_yz, g_0_xz_xx_zz, g_xx_0_0_0_0_xz_xx_xx, g_xx_0_0_0_0_xz_xx_xy, g_xx_0_0_0_0_xz_xx_xz, g_xx_0_0_0_0_xz_xx_yy, g_xx_0_0_0_0_xz_xx_yz, g_xx_0_0_0_0_xz_xx_zz, g_xx_xz_xx_xx, g_xx_xz_xx_xy, g_xx_xz_xx_xz, g_xx_xz_xx_yy, g_xx_xz_xx_yz, g_xx_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xz_xx_xx[i] = -2.0 * g_0_xz_xx_xx[i] * a_exp + 4.0 * g_xx_xz_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_xx_xy[i] = -2.0 * g_0_xz_xx_xy[i] * a_exp + 4.0 * g_xx_xz_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_xx_xz[i] = -2.0 * g_0_xz_xx_xz[i] * a_exp + 4.0 * g_xx_xz_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_xx_yy[i] = -2.0 * g_0_xz_xx_yy[i] * a_exp + 4.0 * g_xx_xz_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_xx_yz[i] = -2.0 * g_0_xz_xx_yz[i] * a_exp + 4.0 * g_xx_xz_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_xx_zz[i] = -2.0 * g_0_xz_xx_zz[i] * a_exp + 4.0 * g_xx_xz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_0_xz_xy_xx, g_0_xz_xy_xy, g_0_xz_xy_xz, g_0_xz_xy_yy, g_0_xz_xy_yz, g_0_xz_xy_zz, g_xx_0_0_0_0_xz_xy_xx, g_xx_0_0_0_0_xz_xy_xy, g_xx_0_0_0_0_xz_xy_xz, g_xx_0_0_0_0_xz_xy_yy, g_xx_0_0_0_0_xz_xy_yz, g_xx_0_0_0_0_xz_xy_zz, g_xx_xz_xy_xx, g_xx_xz_xy_xy, g_xx_xz_xy_xz, g_xx_xz_xy_yy, g_xx_xz_xy_yz, g_xx_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xz_xy_xx[i] = -2.0 * g_0_xz_xy_xx[i] * a_exp + 4.0 * g_xx_xz_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_xy_xy[i] = -2.0 * g_0_xz_xy_xy[i] * a_exp + 4.0 * g_xx_xz_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_xy_xz[i] = -2.0 * g_0_xz_xy_xz[i] * a_exp + 4.0 * g_xx_xz_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_xy_yy[i] = -2.0 * g_0_xz_xy_yy[i] * a_exp + 4.0 * g_xx_xz_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_xy_yz[i] = -2.0 * g_0_xz_xy_yz[i] * a_exp + 4.0 * g_xx_xz_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_xy_zz[i] = -2.0 * g_0_xz_xy_zz[i] * a_exp + 4.0 * g_xx_xz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_0_xz_xz_xx, g_0_xz_xz_xy, g_0_xz_xz_xz, g_0_xz_xz_yy, g_0_xz_xz_yz, g_0_xz_xz_zz, g_xx_0_0_0_0_xz_xz_xx, g_xx_0_0_0_0_xz_xz_xy, g_xx_0_0_0_0_xz_xz_xz, g_xx_0_0_0_0_xz_xz_yy, g_xx_0_0_0_0_xz_xz_yz, g_xx_0_0_0_0_xz_xz_zz, g_xx_xz_xz_xx, g_xx_xz_xz_xy, g_xx_xz_xz_xz, g_xx_xz_xz_yy, g_xx_xz_xz_yz, g_xx_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xz_xz_xx[i] = -2.0 * g_0_xz_xz_xx[i] * a_exp + 4.0 * g_xx_xz_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_xz_xy[i] = -2.0 * g_0_xz_xz_xy[i] * a_exp + 4.0 * g_xx_xz_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_xz_xz[i] = -2.0 * g_0_xz_xz_xz[i] * a_exp + 4.0 * g_xx_xz_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_xz_yy[i] = -2.0 * g_0_xz_xz_yy[i] * a_exp + 4.0 * g_xx_xz_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_xz_yz[i] = -2.0 * g_0_xz_xz_yz[i] * a_exp + 4.0 * g_xx_xz_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_xz_zz[i] = -2.0 * g_0_xz_xz_zz[i] * a_exp + 4.0 * g_xx_xz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_0_xz_yy_xx, g_0_xz_yy_xy, g_0_xz_yy_xz, g_0_xz_yy_yy, g_0_xz_yy_yz, g_0_xz_yy_zz, g_xx_0_0_0_0_xz_yy_xx, g_xx_0_0_0_0_xz_yy_xy, g_xx_0_0_0_0_xz_yy_xz, g_xx_0_0_0_0_xz_yy_yy, g_xx_0_0_0_0_xz_yy_yz, g_xx_0_0_0_0_xz_yy_zz, g_xx_xz_yy_xx, g_xx_xz_yy_xy, g_xx_xz_yy_xz, g_xx_xz_yy_yy, g_xx_xz_yy_yz, g_xx_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xz_yy_xx[i] = -2.0 * g_0_xz_yy_xx[i] * a_exp + 4.0 * g_xx_xz_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_yy_xy[i] = -2.0 * g_0_xz_yy_xy[i] * a_exp + 4.0 * g_xx_xz_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_yy_xz[i] = -2.0 * g_0_xz_yy_xz[i] * a_exp + 4.0 * g_xx_xz_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_yy_yy[i] = -2.0 * g_0_xz_yy_yy[i] * a_exp + 4.0 * g_xx_xz_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_yy_yz[i] = -2.0 * g_0_xz_yy_yz[i] * a_exp + 4.0 * g_xx_xz_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_yy_zz[i] = -2.0 * g_0_xz_yy_zz[i] * a_exp + 4.0 * g_xx_xz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_0_xz_yz_xx, g_0_xz_yz_xy, g_0_xz_yz_xz, g_0_xz_yz_yy, g_0_xz_yz_yz, g_0_xz_yz_zz, g_xx_0_0_0_0_xz_yz_xx, g_xx_0_0_0_0_xz_yz_xy, g_xx_0_0_0_0_xz_yz_xz, g_xx_0_0_0_0_xz_yz_yy, g_xx_0_0_0_0_xz_yz_yz, g_xx_0_0_0_0_xz_yz_zz, g_xx_xz_yz_xx, g_xx_xz_yz_xy, g_xx_xz_yz_xz, g_xx_xz_yz_yy, g_xx_xz_yz_yz, g_xx_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xz_yz_xx[i] = -2.0 * g_0_xz_yz_xx[i] * a_exp + 4.0 * g_xx_xz_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_yz_xy[i] = -2.0 * g_0_xz_yz_xy[i] * a_exp + 4.0 * g_xx_xz_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_yz_xz[i] = -2.0 * g_0_xz_yz_xz[i] * a_exp + 4.0 * g_xx_xz_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_yz_yy[i] = -2.0 * g_0_xz_yz_yy[i] * a_exp + 4.0 * g_xx_xz_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_yz_yz[i] = -2.0 * g_0_xz_yz_yz[i] * a_exp + 4.0 * g_xx_xz_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_yz_zz[i] = -2.0 * g_0_xz_yz_zz[i] * a_exp + 4.0 * g_xx_xz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_0_xz_zz_xx, g_0_xz_zz_xy, g_0_xz_zz_xz, g_0_xz_zz_yy, g_0_xz_zz_yz, g_0_xz_zz_zz, g_xx_0_0_0_0_xz_zz_xx, g_xx_0_0_0_0_xz_zz_xy, g_xx_0_0_0_0_xz_zz_xz, g_xx_0_0_0_0_xz_zz_yy, g_xx_0_0_0_0_xz_zz_yz, g_xx_0_0_0_0_xz_zz_zz, g_xx_xz_zz_xx, g_xx_xz_zz_xy, g_xx_xz_zz_xz, g_xx_xz_zz_yy, g_xx_xz_zz_yz, g_xx_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xz_zz_xx[i] = -2.0 * g_0_xz_zz_xx[i] * a_exp + 4.0 * g_xx_xz_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_zz_xy[i] = -2.0 * g_0_xz_zz_xy[i] * a_exp + 4.0 * g_xx_xz_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_zz_xz[i] = -2.0 * g_0_xz_zz_xz[i] * a_exp + 4.0 * g_xx_xz_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_zz_yy[i] = -2.0 * g_0_xz_zz_yy[i] * a_exp + 4.0 * g_xx_xz_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_zz_yz[i] = -2.0 * g_0_xz_zz_yz[i] * a_exp + 4.0 * g_xx_xz_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_zz_zz[i] = -2.0 * g_0_xz_zz_zz[i] * a_exp + 4.0 * g_xx_xz_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_0_yy_xx_xx, g_0_yy_xx_xy, g_0_yy_xx_xz, g_0_yy_xx_yy, g_0_yy_xx_yz, g_0_yy_xx_zz, g_xx_0_0_0_0_yy_xx_xx, g_xx_0_0_0_0_yy_xx_xy, g_xx_0_0_0_0_yy_xx_xz, g_xx_0_0_0_0_yy_xx_yy, g_xx_0_0_0_0_yy_xx_yz, g_xx_0_0_0_0_yy_xx_zz, g_xx_yy_xx_xx, g_xx_yy_xx_xy, g_xx_yy_xx_xz, g_xx_yy_xx_yy, g_xx_yy_xx_yz, g_xx_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yy_xx_xx[i] = -2.0 * g_0_yy_xx_xx[i] * a_exp + 4.0 * g_xx_yy_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_xx_xy[i] = -2.0 * g_0_yy_xx_xy[i] * a_exp + 4.0 * g_xx_yy_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_xx_xz[i] = -2.0 * g_0_yy_xx_xz[i] * a_exp + 4.0 * g_xx_yy_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_xx_yy[i] = -2.0 * g_0_yy_xx_yy[i] * a_exp + 4.0 * g_xx_yy_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_xx_yz[i] = -2.0 * g_0_yy_xx_yz[i] * a_exp + 4.0 * g_xx_yy_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_xx_zz[i] = -2.0 * g_0_yy_xx_zz[i] * a_exp + 4.0 * g_xx_yy_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_0_yy_xy_xx, g_0_yy_xy_xy, g_0_yy_xy_xz, g_0_yy_xy_yy, g_0_yy_xy_yz, g_0_yy_xy_zz, g_xx_0_0_0_0_yy_xy_xx, g_xx_0_0_0_0_yy_xy_xy, g_xx_0_0_0_0_yy_xy_xz, g_xx_0_0_0_0_yy_xy_yy, g_xx_0_0_0_0_yy_xy_yz, g_xx_0_0_0_0_yy_xy_zz, g_xx_yy_xy_xx, g_xx_yy_xy_xy, g_xx_yy_xy_xz, g_xx_yy_xy_yy, g_xx_yy_xy_yz, g_xx_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yy_xy_xx[i] = -2.0 * g_0_yy_xy_xx[i] * a_exp + 4.0 * g_xx_yy_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_xy_xy[i] = -2.0 * g_0_yy_xy_xy[i] * a_exp + 4.0 * g_xx_yy_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_xy_xz[i] = -2.0 * g_0_yy_xy_xz[i] * a_exp + 4.0 * g_xx_yy_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_xy_yy[i] = -2.0 * g_0_yy_xy_yy[i] * a_exp + 4.0 * g_xx_yy_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_xy_yz[i] = -2.0 * g_0_yy_xy_yz[i] * a_exp + 4.0 * g_xx_yy_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_xy_zz[i] = -2.0 * g_0_yy_xy_zz[i] * a_exp + 4.0 * g_xx_yy_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_0_yy_xz_xx, g_0_yy_xz_xy, g_0_yy_xz_xz, g_0_yy_xz_yy, g_0_yy_xz_yz, g_0_yy_xz_zz, g_xx_0_0_0_0_yy_xz_xx, g_xx_0_0_0_0_yy_xz_xy, g_xx_0_0_0_0_yy_xz_xz, g_xx_0_0_0_0_yy_xz_yy, g_xx_0_0_0_0_yy_xz_yz, g_xx_0_0_0_0_yy_xz_zz, g_xx_yy_xz_xx, g_xx_yy_xz_xy, g_xx_yy_xz_xz, g_xx_yy_xz_yy, g_xx_yy_xz_yz, g_xx_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yy_xz_xx[i] = -2.0 * g_0_yy_xz_xx[i] * a_exp + 4.0 * g_xx_yy_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_xz_xy[i] = -2.0 * g_0_yy_xz_xy[i] * a_exp + 4.0 * g_xx_yy_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_xz_xz[i] = -2.0 * g_0_yy_xz_xz[i] * a_exp + 4.0 * g_xx_yy_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_xz_yy[i] = -2.0 * g_0_yy_xz_yy[i] * a_exp + 4.0 * g_xx_yy_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_xz_yz[i] = -2.0 * g_0_yy_xz_yz[i] * a_exp + 4.0 * g_xx_yy_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_xz_zz[i] = -2.0 * g_0_yy_xz_zz[i] * a_exp + 4.0 * g_xx_yy_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_0_yy_yy_xx, g_0_yy_yy_xy, g_0_yy_yy_xz, g_0_yy_yy_yy, g_0_yy_yy_yz, g_0_yy_yy_zz, g_xx_0_0_0_0_yy_yy_xx, g_xx_0_0_0_0_yy_yy_xy, g_xx_0_0_0_0_yy_yy_xz, g_xx_0_0_0_0_yy_yy_yy, g_xx_0_0_0_0_yy_yy_yz, g_xx_0_0_0_0_yy_yy_zz, g_xx_yy_yy_xx, g_xx_yy_yy_xy, g_xx_yy_yy_xz, g_xx_yy_yy_yy, g_xx_yy_yy_yz, g_xx_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yy_yy_xx[i] = -2.0 * g_0_yy_yy_xx[i] * a_exp + 4.0 * g_xx_yy_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_yy_xy[i] = -2.0 * g_0_yy_yy_xy[i] * a_exp + 4.0 * g_xx_yy_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_yy_xz[i] = -2.0 * g_0_yy_yy_xz[i] * a_exp + 4.0 * g_xx_yy_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_yy_yy[i] = -2.0 * g_0_yy_yy_yy[i] * a_exp + 4.0 * g_xx_yy_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_yy_yz[i] = -2.0 * g_0_yy_yy_yz[i] * a_exp + 4.0 * g_xx_yy_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_yy_zz[i] = -2.0 * g_0_yy_yy_zz[i] * a_exp + 4.0 * g_xx_yy_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_0_yy_yz_xx, g_0_yy_yz_xy, g_0_yy_yz_xz, g_0_yy_yz_yy, g_0_yy_yz_yz, g_0_yy_yz_zz, g_xx_0_0_0_0_yy_yz_xx, g_xx_0_0_0_0_yy_yz_xy, g_xx_0_0_0_0_yy_yz_xz, g_xx_0_0_0_0_yy_yz_yy, g_xx_0_0_0_0_yy_yz_yz, g_xx_0_0_0_0_yy_yz_zz, g_xx_yy_yz_xx, g_xx_yy_yz_xy, g_xx_yy_yz_xz, g_xx_yy_yz_yy, g_xx_yy_yz_yz, g_xx_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yy_yz_xx[i] = -2.0 * g_0_yy_yz_xx[i] * a_exp + 4.0 * g_xx_yy_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_yz_xy[i] = -2.0 * g_0_yy_yz_xy[i] * a_exp + 4.0 * g_xx_yy_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_yz_xz[i] = -2.0 * g_0_yy_yz_xz[i] * a_exp + 4.0 * g_xx_yy_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_yz_yy[i] = -2.0 * g_0_yy_yz_yy[i] * a_exp + 4.0 * g_xx_yy_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_yz_yz[i] = -2.0 * g_0_yy_yz_yz[i] * a_exp + 4.0 * g_xx_yy_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_yz_zz[i] = -2.0 * g_0_yy_yz_zz[i] * a_exp + 4.0 * g_xx_yy_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_0_yy_zz_xx, g_0_yy_zz_xy, g_0_yy_zz_xz, g_0_yy_zz_yy, g_0_yy_zz_yz, g_0_yy_zz_zz, g_xx_0_0_0_0_yy_zz_xx, g_xx_0_0_0_0_yy_zz_xy, g_xx_0_0_0_0_yy_zz_xz, g_xx_0_0_0_0_yy_zz_yy, g_xx_0_0_0_0_yy_zz_yz, g_xx_0_0_0_0_yy_zz_zz, g_xx_yy_zz_xx, g_xx_yy_zz_xy, g_xx_yy_zz_xz, g_xx_yy_zz_yy, g_xx_yy_zz_yz, g_xx_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yy_zz_xx[i] = -2.0 * g_0_yy_zz_xx[i] * a_exp + 4.0 * g_xx_yy_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_zz_xy[i] = -2.0 * g_0_yy_zz_xy[i] * a_exp + 4.0 * g_xx_yy_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_zz_xz[i] = -2.0 * g_0_yy_zz_xz[i] * a_exp + 4.0 * g_xx_yy_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_zz_yy[i] = -2.0 * g_0_yy_zz_yy[i] * a_exp + 4.0 * g_xx_yy_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_zz_yz[i] = -2.0 * g_0_yy_zz_yz[i] * a_exp + 4.0 * g_xx_yy_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_zz_zz[i] = -2.0 * g_0_yy_zz_zz[i] * a_exp + 4.0 * g_xx_yy_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_0_yz_xx_xx, g_0_yz_xx_xy, g_0_yz_xx_xz, g_0_yz_xx_yy, g_0_yz_xx_yz, g_0_yz_xx_zz, g_xx_0_0_0_0_yz_xx_xx, g_xx_0_0_0_0_yz_xx_xy, g_xx_0_0_0_0_yz_xx_xz, g_xx_0_0_0_0_yz_xx_yy, g_xx_0_0_0_0_yz_xx_yz, g_xx_0_0_0_0_yz_xx_zz, g_xx_yz_xx_xx, g_xx_yz_xx_xy, g_xx_yz_xx_xz, g_xx_yz_xx_yy, g_xx_yz_xx_yz, g_xx_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yz_xx_xx[i] = -2.0 * g_0_yz_xx_xx[i] * a_exp + 4.0 * g_xx_yz_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_xx_xy[i] = -2.0 * g_0_yz_xx_xy[i] * a_exp + 4.0 * g_xx_yz_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_xx_xz[i] = -2.0 * g_0_yz_xx_xz[i] * a_exp + 4.0 * g_xx_yz_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_xx_yy[i] = -2.0 * g_0_yz_xx_yy[i] * a_exp + 4.0 * g_xx_yz_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_xx_yz[i] = -2.0 * g_0_yz_xx_yz[i] * a_exp + 4.0 * g_xx_yz_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_xx_zz[i] = -2.0 * g_0_yz_xx_zz[i] * a_exp + 4.0 * g_xx_yz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_0_yz_xy_xx, g_0_yz_xy_xy, g_0_yz_xy_xz, g_0_yz_xy_yy, g_0_yz_xy_yz, g_0_yz_xy_zz, g_xx_0_0_0_0_yz_xy_xx, g_xx_0_0_0_0_yz_xy_xy, g_xx_0_0_0_0_yz_xy_xz, g_xx_0_0_0_0_yz_xy_yy, g_xx_0_0_0_0_yz_xy_yz, g_xx_0_0_0_0_yz_xy_zz, g_xx_yz_xy_xx, g_xx_yz_xy_xy, g_xx_yz_xy_xz, g_xx_yz_xy_yy, g_xx_yz_xy_yz, g_xx_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yz_xy_xx[i] = -2.0 * g_0_yz_xy_xx[i] * a_exp + 4.0 * g_xx_yz_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_xy_xy[i] = -2.0 * g_0_yz_xy_xy[i] * a_exp + 4.0 * g_xx_yz_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_xy_xz[i] = -2.0 * g_0_yz_xy_xz[i] * a_exp + 4.0 * g_xx_yz_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_xy_yy[i] = -2.0 * g_0_yz_xy_yy[i] * a_exp + 4.0 * g_xx_yz_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_xy_yz[i] = -2.0 * g_0_yz_xy_yz[i] * a_exp + 4.0 * g_xx_yz_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_xy_zz[i] = -2.0 * g_0_yz_xy_zz[i] * a_exp + 4.0 * g_xx_yz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_0_yz_xz_xx, g_0_yz_xz_xy, g_0_yz_xz_xz, g_0_yz_xz_yy, g_0_yz_xz_yz, g_0_yz_xz_zz, g_xx_0_0_0_0_yz_xz_xx, g_xx_0_0_0_0_yz_xz_xy, g_xx_0_0_0_0_yz_xz_xz, g_xx_0_0_0_0_yz_xz_yy, g_xx_0_0_0_0_yz_xz_yz, g_xx_0_0_0_0_yz_xz_zz, g_xx_yz_xz_xx, g_xx_yz_xz_xy, g_xx_yz_xz_xz, g_xx_yz_xz_yy, g_xx_yz_xz_yz, g_xx_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yz_xz_xx[i] = -2.0 * g_0_yz_xz_xx[i] * a_exp + 4.0 * g_xx_yz_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_xz_xy[i] = -2.0 * g_0_yz_xz_xy[i] * a_exp + 4.0 * g_xx_yz_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_xz_xz[i] = -2.0 * g_0_yz_xz_xz[i] * a_exp + 4.0 * g_xx_yz_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_xz_yy[i] = -2.0 * g_0_yz_xz_yy[i] * a_exp + 4.0 * g_xx_yz_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_xz_yz[i] = -2.0 * g_0_yz_xz_yz[i] * a_exp + 4.0 * g_xx_yz_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_xz_zz[i] = -2.0 * g_0_yz_xz_zz[i] * a_exp + 4.0 * g_xx_yz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_0_yz_yy_xx, g_0_yz_yy_xy, g_0_yz_yy_xz, g_0_yz_yy_yy, g_0_yz_yy_yz, g_0_yz_yy_zz, g_xx_0_0_0_0_yz_yy_xx, g_xx_0_0_0_0_yz_yy_xy, g_xx_0_0_0_0_yz_yy_xz, g_xx_0_0_0_0_yz_yy_yy, g_xx_0_0_0_0_yz_yy_yz, g_xx_0_0_0_0_yz_yy_zz, g_xx_yz_yy_xx, g_xx_yz_yy_xy, g_xx_yz_yy_xz, g_xx_yz_yy_yy, g_xx_yz_yy_yz, g_xx_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yz_yy_xx[i] = -2.0 * g_0_yz_yy_xx[i] * a_exp + 4.0 * g_xx_yz_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_yy_xy[i] = -2.0 * g_0_yz_yy_xy[i] * a_exp + 4.0 * g_xx_yz_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_yy_xz[i] = -2.0 * g_0_yz_yy_xz[i] * a_exp + 4.0 * g_xx_yz_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_yy_yy[i] = -2.0 * g_0_yz_yy_yy[i] * a_exp + 4.0 * g_xx_yz_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_yy_yz[i] = -2.0 * g_0_yz_yy_yz[i] * a_exp + 4.0 * g_xx_yz_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_yy_zz[i] = -2.0 * g_0_yz_yy_zz[i] * a_exp + 4.0 * g_xx_yz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_0_yz_yz_xx, g_0_yz_yz_xy, g_0_yz_yz_xz, g_0_yz_yz_yy, g_0_yz_yz_yz, g_0_yz_yz_zz, g_xx_0_0_0_0_yz_yz_xx, g_xx_0_0_0_0_yz_yz_xy, g_xx_0_0_0_0_yz_yz_xz, g_xx_0_0_0_0_yz_yz_yy, g_xx_0_0_0_0_yz_yz_yz, g_xx_0_0_0_0_yz_yz_zz, g_xx_yz_yz_xx, g_xx_yz_yz_xy, g_xx_yz_yz_xz, g_xx_yz_yz_yy, g_xx_yz_yz_yz, g_xx_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yz_yz_xx[i] = -2.0 * g_0_yz_yz_xx[i] * a_exp + 4.0 * g_xx_yz_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_yz_xy[i] = -2.0 * g_0_yz_yz_xy[i] * a_exp + 4.0 * g_xx_yz_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_yz_xz[i] = -2.0 * g_0_yz_yz_xz[i] * a_exp + 4.0 * g_xx_yz_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_yz_yy[i] = -2.0 * g_0_yz_yz_yy[i] * a_exp + 4.0 * g_xx_yz_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_yz_yz[i] = -2.0 * g_0_yz_yz_yz[i] * a_exp + 4.0 * g_xx_yz_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_yz_zz[i] = -2.0 * g_0_yz_yz_zz[i] * a_exp + 4.0 * g_xx_yz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_0_yz_zz_xx, g_0_yz_zz_xy, g_0_yz_zz_xz, g_0_yz_zz_yy, g_0_yz_zz_yz, g_0_yz_zz_zz, g_xx_0_0_0_0_yz_zz_xx, g_xx_0_0_0_0_yz_zz_xy, g_xx_0_0_0_0_yz_zz_xz, g_xx_0_0_0_0_yz_zz_yy, g_xx_0_0_0_0_yz_zz_yz, g_xx_0_0_0_0_yz_zz_zz, g_xx_yz_zz_xx, g_xx_yz_zz_xy, g_xx_yz_zz_xz, g_xx_yz_zz_yy, g_xx_yz_zz_yz, g_xx_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yz_zz_xx[i] = -2.0 * g_0_yz_zz_xx[i] * a_exp + 4.0 * g_xx_yz_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_zz_xy[i] = -2.0 * g_0_yz_zz_xy[i] * a_exp + 4.0 * g_xx_yz_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_zz_xz[i] = -2.0 * g_0_yz_zz_xz[i] * a_exp + 4.0 * g_xx_yz_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_zz_yy[i] = -2.0 * g_0_yz_zz_yy[i] * a_exp + 4.0 * g_xx_yz_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_zz_yz[i] = -2.0 * g_0_yz_zz_yz[i] * a_exp + 4.0 * g_xx_yz_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_zz_zz[i] = -2.0 * g_0_yz_zz_zz[i] * a_exp + 4.0 * g_xx_yz_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_0_zz_xx_xx, g_0_zz_xx_xy, g_0_zz_xx_xz, g_0_zz_xx_yy, g_0_zz_xx_yz, g_0_zz_xx_zz, g_xx_0_0_0_0_zz_xx_xx, g_xx_0_0_0_0_zz_xx_xy, g_xx_0_0_0_0_zz_xx_xz, g_xx_0_0_0_0_zz_xx_yy, g_xx_0_0_0_0_zz_xx_yz, g_xx_0_0_0_0_zz_xx_zz, g_xx_zz_xx_xx, g_xx_zz_xx_xy, g_xx_zz_xx_xz, g_xx_zz_xx_yy, g_xx_zz_xx_yz, g_xx_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_zz_xx_xx[i] = -2.0 * g_0_zz_xx_xx[i] * a_exp + 4.0 * g_xx_zz_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_xx_xy[i] = -2.0 * g_0_zz_xx_xy[i] * a_exp + 4.0 * g_xx_zz_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_xx_xz[i] = -2.0 * g_0_zz_xx_xz[i] * a_exp + 4.0 * g_xx_zz_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_xx_yy[i] = -2.0 * g_0_zz_xx_yy[i] * a_exp + 4.0 * g_xx_zz_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_xx_yz[i] = -2.0 * g_0_zz_xx_yz[i] * a_exp + 4.0 * g_xx_zz_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_xx_zz[i] = -2.0 * g_0_zz_xx_zz[i] * a_exp + 4.0 * g_xx_zz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_0_zz_xy_xx, g_0_zz_xy_xy, g_0_zz_xy_xz, g_0_zz_xy_yy, g_0_zz_xy_yz, g_0_zz_xy_zz, g_xx_0_0_0_0_zz_xy_xx, g_xx_0_0_0_0_zz_xy_xy, g_xx_0_0_0_0_zz_xy_xz, g_xx_0_0_0_0_zz_xy_yy, g_xx_0_0_0_0_zz_xy_yz, g_xx_0_0_0_0_zz_xy_zz, g_xx_zz_xy_xx, g_xx_zz_xy_xy, g_xx_zz_xy_xz, g_xx_zz_xy_yy, g_xx_zz_xy_yz, g_xx_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_zz_xy_xx[i] = -2.0 * g_0_zz_xy_xx[i] * a_exp + 4.0 * g_xx_zz_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_xy_xy[i] = -2.0 * g_0_zz_xy_xy[i] * a_exp + 4.0 * g_xx_zz_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_xy_xz[i] = -2.0 * g_0_zz_xy_xz[i] * a_exp + 4.0 * g_xx_zz_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_xy_yy[i] = -2.0 * g_0_zz_xy_yy[i] * a_exp + 4.0 * g_xx_zz_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_xy_yz[i] = -2.0 * g_0_zz_xy_yz[i] * a_exp + 4.0 * g_xx_zz_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_xy_zz[i] = -2.0 * g_0_zz_xy_zz[i] * a_exp + 4.0 * g_xx_zz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_0_zz_xz_xx, g_0_zz_xz_xy, g_0_zz_xz_xz, g_0_zz_xz_yy, g_0_zz_xz_yz, g_0_zz_xz_zz, g_xx_0_0_0_0_zz_xz_xx, g_xx_0_0_0_0_zz_xz_xy, g_xx_0_0_0_0_zz_xz_xz, g_xx_0_0_0_0_zz_xz_yy, g_xx_0_0_0_0_zz_xz_yz, g_xx_0_0_0_0_zz_xz_zz, g_xx_zz_xz_xx, g_xx_zz_xz_xy, g_xx_zz_xz_xz, g_xx_zz_xz_yy, g_xx_zz_xz_yz, g_xx_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_zz_xz_xx[i] = -2.0 * g_0_zz_xz_xx[i] * a_exp + 4.0 * g_xx_zz_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_xz_xy[i] = -2.0 * g_0_zz_xz_xy[i] * a_exp + 4.0 * g_xx_zz_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_xz_xz[i] = -2.0 * g_0_zz_xz_xz[i] * a_exp + 4.0 * g_xx_zz_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_xz_yy[i] = -2.0 * g_0_zz_xz_yy[i] * a_exp + 4.0 * g_xx_zz_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_xz_yz[i] = -2.0 * g_0_zz_xz_yz[i] * a_exp + 4.0 * g_xx_zz_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_xz_zz[i] = -2.0 * g_0_zz_xz_zz[i] * a_exp + 4.0 * g_xx_zz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_0_zz_yy_xx, g_0_zz_yy_xy, g_0_zz_yy_xz, g_0_zz_yy_yy, g_0_zz_yy_yz, g_0_zz_yy_zz, g_xx_0_0_0_0_zz_yy_xx, g_xx_0_0_0_0_zz_yy_xy, g_xx_0_0_0_0_zz_yy_xz, g_xx_0_0_0_0_zz_yy_yy, g_xx_0_0_0_0_zz_yy_yz, g_xx_0_0_0_0_zz_yy_zz, g_xx_zz_yy_xx, g_xx_zz_yy_xy, g_xx_zz_yy_xz, g_xx_zz_yy_yy, g_xx_zz_yy_yz, g_xx_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_zz_yy_xx[i] = -2.0 * g_0_zz_yy_xx[i] * a_exp + 4.0 * g_xx_zz_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_yy_xy[i] = -2.0 * g_0_zz_yy_xy[i] * a_exp + 4.0 * g_xx_zz_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_yy_xz[i] = -2.0 * g_0_zz_yy_xz[i] * a_exp + 4.0 * g_xx_zz_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_yy_yy[i] = -2.0 * g_0_zz_yy_yy[i] * a_exp + 4.0 * g_xx_zz_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_yy_yz[i] = -2.0 * g_0_zz_yy_yz[i] * a_exp + 4.0 * g_xx_zz_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_yy_zz[i] = -2.0 * g_0_zz_yy_zz[i] * a_exp + 4.0 * g_xx_zz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_0_zz_yz_xx, g_0_zz_yz_xy, g_0_zz_yz_xz, g_0_zz_yz_yy, g_0_zz_yz_yz, g_0_zz_yz_zz, g_xx_0_0_0_0_zz_yz_xx, g_xx_0_0_0_0_zz_yz_xy, g_xx_0_0_0_0_zz_yz_xz, g_xx_0_0_0_0_zz_yz_yy, g_xx_0_0_0_0_zz_yz_yz, g_xx_0_0_0_0_zz_yz_zz, g_xx_zz_yz_xx, g_xx_zz_yz_xy, g_xx_zz_yz_xz, g_xx_zz_yz_yy, g_xx_zz_yz_yz, g_xx_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_zz_yz_xx[i] = -2.0 * g_0_zz_yz_xx[i] * a_exp + 4.0 * g_xx_zz_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_yz_xy[i] = -2.0 * g_0_zz_yz_xy[i] * a_exp + 4.0 * g_xx_zz_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_yz_xz[i] = -2.0 * g_0_zz_yz_xz[i] * a_exp + 4.0 * g_xx_zz_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_yz_yy[i] = -2.0 * g_0_zz_yz_yy[i] * a_exp + 4.0 * g_xx_zz_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_yz_yz[i] = -2.0 * g_0_zz_yz_yz[i] * a_exp + 4.0 * g_xx_zz_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_yz_zz[i] = -2.0 * g_0_zz_yz_zz[i] * a_exp + 4.0 * g_xx_zz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_0_zz_zz_xx, g_0_zz_zz_xy, g_0_zz_zz_xz, g_0_zz_zz_yy, g_0_zz_zz_yz, g_0_zz_zz_zz, g_xx_0_0_0_0_zz_zz_xx, g_xx_0_0_0_0_zz_zz_xy, g_xx_0_0_0_0_zz_zz_xz, g_xx_0_0_0_0_zz_zz_yy, g_xx_0_0_0_0_zz_zz_yz, g_xx_0_0_0_0_zz_zz_zz, g_xx_zz_zz_xx, g_xx_zz_zz_xy, g_xx_zz_zz_xz, g_xx_zz_zz_yy, g_xx_zz_zz_yz, g_xx_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_zz_zz_xx[i] = -2.0 * g_0_zz_zz_xx[i] * a_exp + 4.0 * g_xx_zz_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_zz_xy[i] = -2.0 * g_0_zz_zz_xy[i] * a_exp + 4.0 * g_xx_zz_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_zz_xz[i] = -2.0 * g_0_zz_zz_xz[i] * a_exp + 4.0 * g_xx_zz_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_zz_yy[i] = -2.0 * g_0_zz_zz_yy[i] * a_exp + 4.0 * g_xx_zz_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_zz_yz[i] = -2.0 * g_0_zz_zz_yz[i] * a_exp + 4.0 * g_xx_zz_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_zz_zz[i] = -2.0 * g_0_zz_zz_zz[i] * a_exp + 4.0 * g_xx_zz_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_xy_0_0_0_0_xx_xx_xx, g_xy_0_0_0_0_xx_xx_xy, g_xy_0_0_0_0_xx_xx_xz, g_xy_0_0_0_0_xx_xx_yy, g_xy_0_0_0_0_xx_xx_yz, g_xy_0_0_0_0_xx_xx_zz, g_xy_xx_xx_xx, g_xy_xx_xx_xy, g_xy_xx_xx_xz, g_xy_xx_xx_yy, g_xy_xx_xx_yz, g_xy_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xx_xx_xx[i] = 4.0 * g_xy_xx_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_xx_xy[i] = 4.0 * g_xy_xx_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_xx_xz[i] = 4.0 * g_xy_xx_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_xx_yy[i] = 4.0 * g_xy_xx_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_xx_yz[i] = 4.0 * g_xy_xx_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_xx_zz[i] = 4.0 * g_xy_xx_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_xy_0_0_0_0_xx_xy_xx, g_xy_0_0_0_0_xx_xy_xy, g_xy_0_0_0_0_xx_xy_xz, g_xy_0_0_0_0_xx_xy_yy, g_xy_0_0_0_0_xx_xy_yz, g_xy_0_0_0_0_xx_xy_zz, g_xy_xx_xy_xx, g_xy_xx_xy_xy, g_xy_xx_xy_xz, g_xy_xx_xy_yy, g_xy_xx_xy_yz, g_xy_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xx_xy_xx[i] = 4.0 * g_xy_xx_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_xy_xy[i] = 4.0 * g_xy_xx_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_xy_xz[i] = 4.0 * g_xy_xx_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_xy_yy[i] = 4.0 * g_xy_xx_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_xy_yz[i] = 4.0 * g_xy_xx_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_xy_zz[i] = 4.0 * g_xy_xx_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_xy_0_0_0_0_xx_xz_xx, g_xy_0_0_0_0_xx_xz_xy, g_xy_0_0_0_0_xx_xz_xz, g_xy_0_0_0_0_xx_xz_yy, g_xy_0_0_0_0_xx_xz_yz, g_xy_0_0_0_0_xx_xz_zz, g_xy_xx_xz_xx, g_xy_xx_xz_xy, g_xy_xx_xz_xz, g_xy_xx_xz_yy, g_xy_xx_xz_yz, g_xy_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xx_xz_xx[i] = 4.0 * g_xy_xx_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_xz_xy[i] = 4.0 * g_xy_xx_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_xz_xz[i] = 4.0 * g_xy_xx_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_xz_yy[i] = 4.0 * g_xy_xx_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_xz_yz[i] = 4.0 * g_xy_xx_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_xz_zz[i] = 4.0 * g_xy_xx_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_xy_0_0_0_0_xx_yy_xx, g_xy_0_0_0_0_xx_yy_xy, g_xy_0_0_0_0_xx_yy_xz, g_xy_0_0_0_0_xx_yy_yy, g_xy_0_0_0_0_xx_yy_yz, g_xy_0_0_0_0_xx_yy_zz, g_xy_xx_yy_xx, g_xy_xx_yy_xy, g_xy_xx_yy_xz, g_xy_xx_yy_yy, g_xy_xx_yy_yz, g_xy_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xx_yy_xx[i] = 4.0 * g_xy_xx_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_yy_xy[i] = 4.0 * g_xy_xx_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_yy_xz[i] = 4.0 * g_xy_xx_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_yy_yy[i] = 4.0 * g_xy_xx_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_yy_yz[i] = 4.0 * g_xy_xx_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_yy_zz[i] = 4.0 * g_xy_xx_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_xy_0_0_0_0_xx_yz_xx, g_xy_0_0_0_0_xx_yz_xy, g_xy_0_0_0_0_xx_yz_xz, g_xy_0_0_0_0_xx_yz_yy, g_xy_0_0_0_0_xx_yz_yz, g_xy_0_0_0_0_xx_yz_zz, g_xy_xx_yz_xx, g_xy_xx_yz_xy, g_xy_xx_yz_xz, g_xy_xx_yz_yy, g_xy_xx_yz_yz, g_xy_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xx_yz_xx[i] = 4.0 * g_xy_xx_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_yz_xy[i] = 4.0 * g_xy_xx_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_yz_xz[i] = 4.0 * g_xy_xx_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_yz_yy[i] = 4.0 * g_xy_xx_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_yz_yz[i] = 4.0 * g_xy_xx_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_yz_zz[i] = 4.0 * g_xy_xx_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_xy_0_0_0_0_xx_zz_xx, g_xy_0_0_0_0_xx_zz_xy, g_xy_0_0_0_0_xx_zz_xz, g_xy_0_0_0_0_xx_zz_yy, g_xy_0_0_0_0_xx_zz_yz, g_xy_0_0_0_0_xx_zz_zz, g_xy_xx_zz_xx, g_xy_xx_zz_xy, g_xy_xx_zz_xz, g_xy_xx_zz_yy, g_xy_xx_zz_yz, g_xy_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xx_zz_xx[i] = 4.0 * g_xy_xx_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_zz_xy[i] = 4.0 * g_xy_xx_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_zz_xz[i] = 4.0 * g_xy_xx_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_zz_yy[i] = 4.0 * g_xy_xx_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_zz_yz[i] = 4.0 * g_xy_xx_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_zz_zz[i] = 4.0 * g_xy_xx_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_xy_0_0_0_0_xy_xx_xx, g_xy_0_0_0_0_xy_xx_xy, g_xy_0_0_0_0_xy_xx_xz, g_xy_0_0_0_0_xy_xx_yy, g_xy_0_0_0_0_xy_xx_yz, g_xy_0_0_0_0_xy_xx_zz, g_xy_xy_xx_xx, g_xy_xy_xx_xy, g_xy_xy_xx_xz, g_xy_xy_xx_yy, g_xy_xy_xx_yz, g_xy_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xy_xx_xx[i] = 4.0 * g_xy_xy_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_xx_xy[i] = 4.0 * g_xy_xy_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_xx_xz[i] = 4.0 * g_xy_xy_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_xx_yy[i] = 4.0 * g_xy_xy_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_xx_yz[i] = 4.0 * g_xy_xy_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_xx_zz[i] = 4.0 * g_xy_xy_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_xy_0_0_0_0_xy_xy_xx, g_xy_0_0_0_0_xy_xy_xy, g_xy_0_0_0_0_xy_xy_xz, g_xy_0_0_0_0_xy_xy_yy, g_xy_0_0_0_0_xy_xy_yz, g_xy_0_0_0_0_xy_xy_zz, g_xy_xy_xy_xx, g_xy_xy_xy_xy, g_xy_xy_xy_xz, g_xy_xy_xy_yy, g_xy_xy_xy_yz, g_xy_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xy_xy_xx[i] = 4.0 * g_xy_xy_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_xy_xy[i] = 4.0 * g_xy_xy_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_xy_xz[i] = 4.0 * g_xy_xy_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_xy_yy[i] = 4.0 * g_xy_xy_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_xy_yz[i] = 4.0 * g_xy_xy_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_xy_zz[i] = 4.0 * g_xy_xy_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_xy_0_0_0_0_xy_xz_xx, g_xy_0_0_0_0_xy_xz_xy, g_xy_0_0_0_0_xy_xz_xz, g_xy_0_0_0_0_xy_xz_yy, g_xy_0_0_0_0_xy_xz_yz, g_xy_0_0_0_0_xy_xz_zz, g_xy_xy_xz_xx, g_xy_xy_xz_xy, g_xy_xy_xz_xz, g_xy_xy_xz_yy, g_xy_xy_xz_yz, g_xy_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xy_xz_xx[i] = 4.0 * g_xy_xy_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_xz_xy[i] = 4.0 * g_xy_xy_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_xz_xz[i] = 4.0 * g_xy_xy_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_xz_yy[i] = 4.0 * g_xy_xy_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_xz_yz[i] = 4.0 * g_xy_xy_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_xz_zz[i] = 4.0 * g_xy_xy_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_xy_0_0_0_0_xy_yy_xx, g_xy_0_0_0_0_xy_yy_xy, g_xy_0_0_0_0_xy_yy_xz, g_xy_0_0_0_0_xy_yy_yy, g_xy_0_0_0_0_xy_yy_yz, g_xy_0_0_0_0_xy_yy_zz, g_xy_xy_yy_xx, g_xy_xy_yy_xy, g_xy_xy_yy_xz, g_xy_xy_yy_yy, g_xy_xy_yy_yz, g_xy_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xy_yy_xx[i] = 4.0 * g_xy_xy_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_yy_xy[i] = 4.0 * g_xy_xy_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_yy_xz[i] = 4.0 * g_xy_xy_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_yy_yy[i] = 4.0 * g_xy_xy_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_yy_yz[i] = 4.0 * g_xy_xy_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_yy_zz[i] = 4.0 * g_xy_xy_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_xy_0_0_0_0_xy_yz_xx, g_xy_0_0_0_0_xy_yz_xy, g_xy_0_0_0_0_xy_yz_xz, g_xy_0_0_0_0_xy_yz_yy, g_xy_0_0_0_0_xy_yz_yz, g_xy_0_0_0_0_xy_yz_zz, g_xy_xy_yz_xx, g_xy_xy_yz_xy, g_xy_xy_yz_xz, g_xy_xy_yz_yy, g_xy_xy_yz_yz, g_xy_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xy_yz_xx[i] = 4.0 * g_xy_xy_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_yz_xy[i] = 4.0 * g_xy_xy_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_yz_xz[i] = 4.0 * g_xy_xy_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_yz_yy[i] = 4.0 * g_xy_xy_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_yz_yz[i] = 4.0 * g_xy_xy_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_yz_zz[i] = 4.0 * g_xy_xy_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_xy_0_0_0_0_xy_zz_xx, g_xy_0_0_0_0_xy_zz_xy, g_xy_0_0_0_0_xy_zz_xz, g_xy_0_0_0_0_xy_zz_yy, g_xy_0_0_0_0_xy_zz_yz, g_xy_0_0_0_0_xy_zz_zz, g_xy_xy_zz_xx, g_xy_xy_zz_xy, g_xy_xy_zz_xz, g_xy_xy_zz_yy, g_xy_xy_zz_yz, g_xy_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xy_zz_xx[i] = 4.0 * g_xy_xy_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_zz_xy[i] = 4.0 * g_xy_xy_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_zz_xz[i] = 4.0 * g_xy_xy_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_zz_yy[i] = 4.0 * g_xy_xy_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_zz_yz[i] = 4.0 * g_xy_xy_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_zz_zz[i] = 4.0 * g_xy_xy_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_xy_0_0_0_0_xz_xx_xx, g_xy_0_0_0_0_xz_xx_xy, g_xy_0_0_0_0_xz_xx_xz, g_xy_0_0_0_0_xz_xx_yy, g_xy_0_0_0_0_xz_xx_yz, g_xy_0_0_0_0_xz_xx_zz, g_xy_xz_xx_xx, g_xy_xz_xx_xy, g_xy_xz_xx_xz, g_xy_xz_xx_yy, g_xy_xz_xx_yz, g_xy_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xz_xx_xx[i] = 4.0 * g_xy_xz_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_xx_xy[i] = 4.0 * g_xy_xz_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_xx_xz[i] = 4.0 * g_xy_xz_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_xx_yy[i] = 4.0 * g_xy_xz_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_xx_yz[i] = 4.0 * g_xy_xz_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_xx_zz[i] = 4.0 * g_xy_xz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_xy_0_0_0_0_xz_xy_xx, g_xy_0_0_0_0_xz_xy_xy, g_xy_0_0_0_0_xz_xy_xz, g_xy_0_0_0_0_xz_xy_yy, g_xy_0_0_0_0_xz_xy_yz, g_xy_0_0_0_0_xz_xy_zz, g_xy_xz_xy_xx, g_xy_xz_xy_xy, g_xy_xz_xy_xz, g_xy_xz_xy_yy, g_xy_xz_xy_yz, g_xy_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xz_xy_xx[i] = 4.0 * g_xy_xz_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_xy_xy[i] = 4.0 * g_xy_xz_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_xy_xz[i] = 4.0 * g_xy_xz_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_xy_yy[i] = 4.0 * g_xy_xz_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_xy_yz[i] = 4.0 * g_xy_xz_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_xy_zz[i] = 4.0 * g_xy_xz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_xy_0_0_0_0_xz_xz_xx, g_xy_0_0_0_0_xz_xz_xy, g_xy_0_0_0_0_xz_xz_xz, g_xy_0_0_0_0_xz_xz_yy, g_xy_0_0_0_0_xz_xz_yz, g_xy_0_0_0_0_xz_xz_zz, g_xy_xz_xz_xx, g_xy_xz_xz_xy, g_xy_xz_xz_xz, g_xy_xz_xz_yy, g_xy_xz_xz_yz, g_xy_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xz_xz_xx[i] = 4.0 * g_xy_xz_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_xz_xy[i] = 4.0 * g_xy_xz_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_xz_xz[i] = 4.0 * g_xy_xz_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_xz_yy[i] = 4.0 * g_xy_xz_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_xz_yz[i] = 4.0 * g_xy_xz_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_xz_zz[i] = 4.0 * g_xy_xz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_xy_0_0_0_0_xz_yy_xx, g_xy_0_0_0_0_xz_yy_xy, g_xy_0_0_0_0_xz_yy_xz, g_xy_0_0_0_0_xz_yy_yy, g_xy_0_0_0_0_xz_yy_yz, g_xy_0_0_0_0_xz_yy_zz, g_xy_xz_yy_xx, g_xy_xz_yy_xy, g_xy_xz_yy_xz, g_xy_xz_yy_yy, g_xy_xz_yy_yz, g_xy_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xz_yy_xx[i] = 4.0 * g_xy_xz_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_yy_xy[i] = 4.0 * g_xy_xz_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_yy_xz[i] = 4.0 * g_xy_xz_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_yy_yy[i] = 4.0 * g_xy_xz_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_yy_yz[i] = 4.0 * g_xy_xz_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_yy_zz[i] = 4.0 * g_xy_xz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_xy_0_0_0_0_xz_yz_xx, g_xy_0_0_0_0_xz_yz_xy, g_xy_0_0_0_0_xz_yz_xz, g_xy_0_0_0_0_xz_yz_yy, g_xy_0_0_0_0_xz_yz_yz, g_xy_0_0_0_0_xz_yz_zz, g_xy_xz_yz_xx, g_xy_xz_yz_xy, g_xy_xz_yz_xz, g_xy_xz_yz_yy, g_xy_xz_yz_yz, g_xy_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xz_yz_xx[i] = 4.0 * g_xy_xz_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_yz_xy[i] = 4.0 * g_xy_xz_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_yz_xz[i] = 4.0 * g_xy_xz_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_yz_yy[i] = 4.0 * g_xy_xz_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_yz_yz[i] = 4.0 * g_xy_xz_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_yz_zz[i] = 4.0 * g_xy_xz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_xy_0_0_0_0_xz_zz_xx, g_xy_0_0_0_0_xz_zz_xy, g_xy_0_0_0_0_xz_zz_xz, g_xy_0_0_0_0_xz_zz_yy, g_xy_0_0_0_0_xz_zz_yz, g_xy_0_0_0_0_xz_zz_zz, g_xy_xz_zz_xx, g_xy_xz_zz_xy, g_xy_xz_zz_xz, g_xy_xz_zz_yy, g_xy_xz_zz_yz, g_xy_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xz_zz_xx[i] = 4.0 * g_xy_xz_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_zz_xy[i] = 4.0 * g_xy_xz_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_zz_xz[i] = 4.0 * g_xy_xz_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_zz_yy[i] = 4.0 * g_xy_xz_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_zz_yz[i] = 4.0 * g_xy_xz_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_zz_zz[i] = 4.0 * g_xy_xz_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_xy_0_0_0_0_yy_xx_xx, g_xy_0_0_0_0_yy_xx_xy, g_xy_0_0_0_0_yy_xx_xz, g_xy_0_0_0_0_yy_xx_yy, g_xy_0_0_0_0_yy_xx_yz, g_xy_0_0_0_0_yy_xx_zz, g_xy_yy_xx_xx, g_xy_yy_xx_xy, g_xy_yy_xx_xz, g_xy_yy_xx_yy, g_xy_yy_xx_yz, g_xy_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yy_xx_xx[i] = 4.0 * g_xy_yy_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_xx_xy[i] = 4.0 * g_xy_yy_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_xx_xz[i] = 4.0 * g_xy_yy_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_xx_yy[i] = 4.0 * g_xy_yy_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_xx_yz[i] = 4.0 * g_xy_yy_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_xx_zz[i] = 4.0 * g_xy_yy_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_xy_0_0_0_0_yy_xy_xx, g_xy_0_0_0_0_yy_xy_xy, g_xy_0_0_0_0_yy_xy_xz, g_xy_0_0_0_0_yy_xy_yy, g_xy_0_0_0_0_yy_xy_yz, g_xy_0_0_0_0_yy_xy_zz, g_xy_yy_xy_xx, g_xy_yy_xy_xy, g_xy_yy_xy_xz, g_xy_yy_xy_yy, g_xy_yy_xy_yz, g_xy_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yy_xy_xx[i] = 4.0 * g_xy_yy_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_xy_xy[i] = 4.0 * g_xy_yy_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_xy_xz[i] = 4.0 * g_xy_yy_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_xy_yy[i] = 4.0 * g_xy_yy_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_xy_yz[i] = 4.0 * g_xy_yy_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_xy_zz[i] = 4.0 * g_xy_yy_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_xy_0_0_0_0_yy_xz_xx, g_xy_0_0_0_0_yy_xz_xy, g_xy_0_0_0_0_yy_xz_xz, g_xy_0_0_0_0_yy_xz_yy, g_xy_0_0_0_0_yy_xz_yz, g_xy_0_0_0_0_yy_xz_zz, g_xy_yy_xz_xx, g_xy_yy_xz_xy, g_xy_yy_xz_xz, g_xy_yy_xz_yy, g_xy_yy_xz_yz, g_xy_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yy_xz_xx[i] = 4.0 * g_xy_yy_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_xz_xy[i] = 4.0 * g_xy_yy_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_xz_xz[i] = 4.0 * g_xy_yy_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_xz_yy[i] = 4.0 * g_xy_yy_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_xz_yz[i] = 4.0 * g_xy_yy_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_xz_zz[i] = 4.0 * g_xy_yy_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_xy_0_0_0_0_yy_yy_xx, g_xy_0_0_0_0_yy_yy_xy, g_xy_0_0_0_0_yy_yy_xz, g_xy_0_0_0_0_yy_yy_yy, g_xy_0_0_0_0_yy_yy_yz, g_xy_0_0_0_0_yy_yy_zz, g_xy_yy_yy_xx, g_xy_yy_yy_xy, g_xy_yy_yy_xz, g_xy_yy_yy_yy, g_xy_yy_yy_yz, g_xy_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yy_yy_xx[i] = 4.0 * g_xy_yy_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_yy_xy[i] = 4.0 * g_xy_yy_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_yy_xz[i] = 4.0 * g_xy_yy_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_yy_yy[i] = 4.0 * g_xy_yy_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_yy_yz[i] = 4.0 * g_xy_yy_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_yy_zz[i] = 4.0 * g_xy_yy_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_xy_0_0_0_0_yy_yz_xx, g_xy_0_0_0_0_yy_yz_xy, g_xy_0_0_0_0_yy_yz_xz, g_xy_0_0_0_0_yy_yz_yy, g_xy_0_0_0_0_yy_yz_yz, g_xy_0_0_0_0_yy_yz_zz, g_xy_yy_yz_xx, g_xy_yy_yz_xy, g_xy_yy_yz_xz, g_xy_yy_yz_yy, g_xy_yy_yz_yz, g_xy_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yy_yz_xx[i] = 4.0 * g_xy_yy_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_yz_xy[i] = 4.0 * g_xy_yy_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_yz_xz[i] = 4.0 * g_xy_yy_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_yz_yy[i] = 4.0 * g_xy_yy_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_yz_yz[i] = 4.0 * g_xy_yy_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_yz_zz[i] = 4.0 * g_xy_yy_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_xy_0_0_0_0_yy_zz_xx, g_xy_0_0_0_0_yy_zz_xy, g_xy_0_0_0_0_yy_zz_xz, g_xy_0_0_0_0_yy_zz_yy, g_xy_0_0_0_0_yy_zz_yz, g_xy_0_0_0_0_yy_zz_zz, g_xy_yy_zz_xx, g_xy_yy_zz_xy, g_xy_yy_zz_xz, g_xy_yy_zz_yy, g_xy_yy_zz_yz, g_xy_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yy_zz_xx[i] = 4.0 * g_xy_yy_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_zz_xy[i] = 4.0 * g_xy_yy_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_zz_xz[i] = 4.0 * g_xy_yy_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_zz_yy[i] = 4.0 * g_xy_yy_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_zz_yz[i] = 4.0 * g_xy_yy_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_zz_zz[i] = 4.0 * g_xy_yy_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_xy_0_0_0_0_yz_xx_xx, g_xy_0_0_0_0_yz_xx_xy, g_xy_0_0_0_0_yz_xx_xz, g_xy_0_0_0_0_yz_xx_yy, g_xy_0_0_0_0_yz_xx_yz, g_xy_0_0_0_0_yz_xx_zz, g_xy_yz_xx_xx, g_xy_yz_xx_xy, g_xy_yz_xx_xz, g_xy_yz_xx_yy, g_xy_yz_xx_yz, g_xy_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yz_xx_xx[i] = 4.0 * g_xy_yz_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_xx_xy[i] = 4.0 * g_xy_yz_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_xx_xz[i] = 4.0 * g_xy_yz_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_xx_yy[i] = 4.0 * g_xy_yz_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_xx_yz[i] = 4.0 * g_xy_yz_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_xx_zz[i] = 4.0 * g_xy_yz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_xy_0_0_0_0_yz_xy_xx, g_xy_0_0_0_0_yz_xy_xy, g_xy_0_0_0_0_yz_xy_xz, g_xy_0_0_0_0_yz_xy_yy, g_xy_0_0_0_0_yz_xy_yz, g_xy_0_0_0_0_yz_xy_zz, g_xy_yz_xy_xx, g_xy_yz_xy_xy, g_xy_yz_xy_xz, g_xy_yz_xy_yy, g_xy_yz_xy_yz, g_xy_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yz_xy_xx[i] = 4.0 * g_xy_yz_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_xy_xy[i] = 4.0 * g_xy_yz_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_xy_xz[i] = 4.0 * g_xy_yz_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_xy_yy[i] = 4.0 * g_xy_yz_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_xy_yz[i] = 4.0 * g_xy_yz_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_xy_zz[i] = 4.0 * g_xy_yz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_xy_0_0_0_0_yz_xz_xx, g_xy_0_0_0_0_yz_xz_xy, g_xy_0_0_0_0_yz_xz_xz, g_xy_0_0_0_0_yz_xz_yy, g_xy_0_0_0_0_yz_xz_yz, g_xy_0_0_0_0_yz_xz_zz, g_xy_yz_xz_xx, g_xy_yz_xz_xy, g_xy_yz_xz_xz, g_xy_yz_xz_yy, g_xy_yz_xz_yz, g_xy_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yz_xz_xx[i] = 4.0 * g_xy_yz_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_xz_xy[i] = 4.0 * g_xy_yz_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_xz_xz[i] = 4.0 * g_xy_yz_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_xz_yy[i] = 4.0 * g_xy_yz_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_xz_yz[i] = 4.0 * g_xy_yz_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_xz_zz[i] = 4.0 * g_xy_yz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_xy_0_0_0_0_yz_yy_xx, g_xy_0_0_0_0_yz_yy_xy, g_xy_0_0_0_0_yz_yy_xz, g_xy_0_0_0_0_yz_yy_yy, g_xy_0_0_0_0_yz_yy_yz, g_xy_0_0_0_0_yz_yy_zz, g_xy_yz_yy_xx, g_xy_yz_yy_xy, g_xy_yz_yy_xz, g_xy_yz_yy_yy, g_xy_yz_yy_yz, g_xy_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yz_yy_xx[i] = 4.0 * g_xy_yz_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_yy_xy[i] = 4.0 * g_xy_yz_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_yy_xz[i] = 4.0 * g_xy_yz_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_yy_yy[i] = 4.0 * g_xy_yz_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_yy_yz[i] = 4.0 * g_xy_yz_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_yy_zz[i] = 4.0 * g_xy_yz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_xy_0_0_0_0_yz_yz_xx, g_xy_0_0_0_0_yz_yz_xy, g_xy_0_0_0_0_yz_yz_xz, g_xy_0_0_0_0_yz_yz_yy, g_xy_0_0_0_0_yz_yz_yz, g_xy_0_0_0_0_yz_yz_zz, g_xy_yz_yz_xx, g_xy_yz_yz_xy, g_xy_yz_yz_xz, g_xy_yz_yz_yy, g_xy_yz_yz_yz, g_xy_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yz_yz_xx[i] = 4.0 * g_xy_yz_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_yz_xy[i] = 4.0 * g_xy_yz_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_yz_xz[i] = 4.0 * g_xy_yz_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_yz_yy[i] = 4.0 * g_xy_yz_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_yz_yz[i] = 4.0 * g_xy_yz_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_yz_zz[i] = 4.0 * g_xy_yz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_xy_0_0_0_0_yz_zz_xx, g_xy_0_0_0_0_yz_zz_xy, g_xy_0_0_0_0_yz_zz_xz, g_xy_0_0_0_0_yz_zz_yy, g_xy_0_0_0_0_yz_zz_yz, g_xy_0_0_0_0_yz_zz_zz, g_xy_yz_zz_xx, g_xy_yz_zz_xy, g_xy_yz_zz_xz, g_xy_yz_zz_yy, g_xy_yz_zz_yz, g_xy_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yz_zz_xx[i] = 4.0 * g_xy_yz_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_zz_xy[i] = 4.0 * g_xy_yz_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_zz_xz[i] = 4.0 * g_xy_yz_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_zz_yy[i] = 4.0 * g_xy_yz_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_zz_yz[i] = 4.0 * g_xy_yz_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_zz_zz[i] = 4.0 * g_xy_yz_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_xy_0_0_0_0_zz_xx_xx, g_xy_0_0_0_0_zz_xx_xy, g_xy_0_0_0_0_zz_xx_xz, g_xy_0_0_0_0_zz_xx_yy, g_xy_0_0_0_0_zz_xx_yz, g_xy_0_0_0_0_zz_xx_zz, g_xy_zz_xx_xx, g_xy_zz_xx_xy, g_xy_zz_xx_xz, g_xy_zz_xx_yy, g_xy_zz_xx_yz, g_xy_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_zz_xx_xx[i] = 4.0 * g_xy_zz_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_xx_xy[i] = 4.0 * g_xy_zz_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_xx_xz[i] = 4.0 * g_xy_zz_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_xx_yy[i] = 4.0 * g_xy_zz_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_xx_yz[i] = 4.0 * g_xy_zz_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_xx_zz[i] = 4.0 * g_xy_zz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_xy_0_0_0_0_zz_xy_xx, g_xy_0_0_0_0_zz_xy_xy, g_xy_0_0_0_0_zz_xy_xz, g_xy_0_0_0_0_zz_xy_yy, g_xy_0_0_0_0_zz_xy_yz, g_xy_0_0_0_0_zz_xy_zz, g_xy_zz_xy_xx, g_xy_zz_xy_xy, g_xy_zz_xy_xz, g_xy_zz_xy_yy, g_xy_zz_xy_yz, g_xy_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_zz_xy_xx[i] = 4.0 * g_xy_zz_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_xy_xy[i] = 4.0 * g_xy_zz_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_xy_xz[i] = 4.0 * g_xy_zz_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_xy_yy[i] = 4.0 * g_xy_zz_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_xy_yz[i] = 4.0 * g_xy_zz_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_xy_zz[i] = 4.0 * g_xy_zz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_xy_0_0_0_0_zz_xz_xx, g_xy_0_0_0_0_zz_xz_xy, g_xy_0_0_0_0_zz_xz_xz, g_xy_0_0_0_0_zz_xz_yy, g_xy_0_0_0_0_zz_xz_yz, g_xy_0_0_0_0_zz_xz_zz, g_xy_zz_xz_xx, g_xy_zz_xz_xy, g_xy_zz_xz_xz, g_xy_zz_xz_yy, g_xy_zz_xz_yz, g_xy_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_zz_xz_xx[i] = 4.0 * g_xy_zz_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_xz_xy[i] = 4.0 * g_xy_zz_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_xz_xz[i] = 4.0 * g_xy_zz_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_xz_yy[i] = 4.0 * g_xy_zz_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_xz_yz[i] = 4.0 * g_xy_zz_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_xz_zz[i] = 4.0 * g_xy_zz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_xy_0_0_0_0_zz_yy_xx, g_xy_0_0_0_0_zz_yy_xy, g_xy_0_0_0_0_zz_yy_xz, g_xy_0_0_0_0_zz_yy_yy, g_xy_0_0_0_0_zz_yy_yz, g_xy_0_0_0_0_zz_yy_zz, g_xy_zz_yy_xx, g_xy_zz_yy_xy, g_xy_zz_yy_xz, g_xy_zz_yy_yy, g_xy_zz_yy_yz, g_xy_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_zz_yy_xx[i] = 4.0 * g_xy_zz_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_yy_xy[i] = 4.0 * g_xy_zz_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_yy_xz[i] = 4.0 * g_xy_zz_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_yy_yy[i] = 4.0 * g_xy_zz_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_yy_yz[i] = 4.0 * g_xy_zz_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_yy_zz[i] = 4.0 * g_xy_zz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_xy_0_0_0_0_zz_yz_xx, g_xy_0_0_0_0_zz_yz_xy, g_xy_0_0_0_0_zz_yz_xz, g_xy_0_0_0_0_zz_yz_yy, g_xy_0_0_0_0_zz_yz_yz, g_xy_0_0_0_0_zz_yz_zz, g_xy_zz_yz_xx, g_xy_zz_yz_xy, g_xy_zz_yz_xz, g_xy_zz_yz_yy, g_xy_zz_yz_yz, g_xy_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_zz_yz_xx[i] = 4.0 * g_xy_zz_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_yz_xy[i] = 4.0 * g_xy_zz_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_yz_xz[i] = 4.0 * g_xy_zz_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_yz_yy[i] = 4.0 * g_xy_zz_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_yz_yz[i] = 4.0 * g_xy_zz_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_yz_zz[i] = 4.0 * g_xy_zz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_xy_0_0_0_0_zz_zz_xx, g_xy_0_0_0_0_zz_zz_xy, g_xy_0_0_0_0_zz_zz_xz, g_xy_0_0_0_0_zz_zz_yy, g_xy_0_0_0_0_zz_zz_yz, g_xy_0_0_0_0_zz_zz_zz, g_xy_zz_zz_xx, g_xy_zz_zz_xy, g_xy_zz_zz_xz, g_xy_zz_zz_yy, g_xy_zz_zz_yz, g_xy_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_zz_zz_xx[i] = 4.0 * g_xy_zz_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_zz_xy[i] = 4.0 * g_xy_zz_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_zz_xz[i] = 4.0 * g_xy_zz_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_zz_yy[i] = 4.0 * g_xy_zz_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_zz_yz[i] = 4.0 * g_xy_zz_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_zz_zz[i] = 4.0 * g_xy_zz_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_xz_0_0_0_0_xx_xx_xx, g_xz_0_0_0_0_xx_xx_xy, g_xz_0_0_0_0_xx_xx_xz, g_xz_0_0_0_0_xx_xx_yy, g_xz_0_0_0_0_xx_xx_yz, g_xz_0_0_0_0_xx_xx_zz, g_xz_xx_xx_xx, g_xz_xx_xx_xy, g_xz_xx_xx_xz, g_xz_xx_xx_yy, g_xz_xx_xx_yz, g_xz_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xx_xx_xx[i] = 4.0 * g_xz_xx_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_xx_xy[i] = 4.0 * g_xz_xx_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_xx_xz[i] = 4.0 * g_xz_xx_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_xx_yy[i] = 4.0 * g_xz_xx_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_xx_yz[i] = 4.0 * g_xz_xx_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_xx_zz[i] = 4.0 * g_xz_xx_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_xz_0_0_0_0_xx_xy_xx, g_xz_0_0_0_0_xx_xy_xy, g_xz_0_0_0_0_xx_xy_xz, g_xz_0_0_0_0_xx_xy_yy, g_xz_0_0_0_0_xx_xy_yz, g_xz_0_0_0_0_xx_xy_zz, g_xz_xx_xy_xx, g_xz_xx_xy_xy, g_xz_xx_xy_xz, g_xz_xx_xy_yy, g_xz_xx_xy_yz, g_xz_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xx_xy_xx[i] = 4.0 * g_xz_xx_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_xy_xy[i] = 4.0 * g_xz_xx_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_xy_xz[i] = 4.0 * g_xz_xx_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_xy_yy[i] = 4.0 * g_xz_xx_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_xy_yz[i] = 4.0 * g_xz_xx_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_xy_zz[i] = 4.0 * g_xz_xx_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_xz_0_0_0_0_xx_xz_xx, g_xz_0_0_0_0_xx_xz_xy, g_xz_0_0_0_0_xx_xz_xz, g_xz_0_0_0_0_xx_xz_yy, g_xz_0_0_0_0_xx_xz_yz, g_xz_0_0_0_0_xx_xz_zz, g_xz_xx_xz_xx, g_xz_xx_xz_xy, g_xz_xx_xz_xz, g_xz_xx_xz_yy, g_xz_xx_xz_yz, g_xz_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xx_xz_xx[i] = 4.0 * g_xz_xx_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_xz_xy[i] = 4.0 * g_xz_xx_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_xz_xz[i] = 4.0 * g_xz_xx_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_xz_yy[i] = 4.0 * g_xz_xx_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_xz_yz[i] = 4.0 * g_xz_xx_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_xz_zz[i] = 4.0 * g_xz_xx_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_xz_0_0_0_0_xx_yy_xx, g_xz_0_0_0_0_xx_yy_xy, g_xz_0_0_0_0_xx_yy_xz, g_xz_0_0_0_0_xx_yy_yy, g_xz_0_0_0_0_xx_yy_yz, g_xz_0_0_0_0_xx_yy_zz, g_xz_xx_yy_xx, g_xz_xx_yy_xy, g_xz_xx_yy_xz, g_xz_xx_yy_yy, g_xz_xx_yy_yz, g_xz_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xx_yy_xx[i] = 4.0 * g_xz_xx_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_yy_xy[i] = 4.0 * g_xz_xx_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_yy_xz[i] = 4.0 * g_xz_xx_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_yy_yy[i] = 4.0 * g_xz_xx_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_yy_yz[i] = 4.0 * g_xz_xx_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_yy_zz[i] = 4.0 * g_xz_xx_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_xz_0_0_0_0_xx_yz_xx, g_xz_0_0_0_0_xx_yz_xy, g_xz_0_0_0_0_xx_yz_xz, g_xz_0_0_0_0_xx_yz_yy, g_xz_0_0_0_0_xx_yz_yz, g_xz_0_0_0_0_xx_yz_zz, g_xz_xx_yz_xx, g_xz_xx_yz_xy, g_xz_xx_yz_xz, g_xz_xx_yz_yy, g_xz_xx_yz_yz, g_xz_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xx_yz_xx[i] = 4.0 * g_xz_xx_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_yz_xy[i] = 4.0 * g_xz_xx_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_yz_xz[i] = 4.0 * g_xz_xx_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_yz_yy[i] = 4.0 * g_xz_xx_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_yz_yz[i] = 4.0 * g_xz_xx_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_yz_zz[i] = 4.0 * g_xz_xx_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_xz_0_0_0_0_xx_zz_xx, g_xz_0_0_0_0_xx_zz_xy, g_xz_0_0_0_0_xx_zz_xz, g_xz_0_0_0_0_xx_zz_yy, g_xz_0_0_0_0_xx_zz_yz, g_xz_0_0_0_0_xx_zz_zz, g_xz_xx_zz_xx, g_xz_xx_zz_xy, g_xz_xx_zz_xz, g_xz_xx_zz_yy, g_xz_xx_zz_yz, g_xz_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xx_zz_xx[i] = 4.0 * g_xz_xx_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_zz_xy[i] = 4.0 * g_xz_xx_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_zz_xz[i] = 4.0 * g_xz_xx_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_zz_yy[i] = 4.0 * g_xz_xx_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_zz_yz[i] = 4.0 * g_xz_xx_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_zz_zz[i] = 4.0 * g_xz_xx_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_xz_0_0_0_0_xy_xx_xx, g_xz_0_0_0_0_xy_xx_xy, g_xz_0_0_0_0_xy_xx_xz, g_xz_0_0_0_0_xy_xx_yy, g_xz_0_0_0_0_xy_xx_yz, g_xz_0_0_0_0_xy_xx_zz, g_xz_xy_xx_xx, g_xz_xy_xx_xy, g_xz_xy_xx_xz, g_xz_xy_xx_yy, g_xz_xy_xx_yz, g_xz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xy_xx_xx[i] = 4.0 * g_xz_xy_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_xx_xy[i] = 4.0 * g_xz_xy_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_xx_xz[i] = 4.0 * g_xz_xy_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_xx_yy[i] = 4.0 * g_xz_xy_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_xx_yz[i] = 4.0 * g_xz_xy_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_xx_zz[i] = 4.0 * g_xz_xy_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_xz_0_0_0_0_xy_xy_xx, g_xz_0_0_0_0_xy_xy_xy, g_xz_0_0_0_0_xy_xy_xz, g_xz_0_0_0_0_xy_xy_yy, g_xz_0_0_0_0_xy_xy_yz, g_xz_0_0_0_0_xy_xy_zz, g_xz_xy_xy_xx, g_xz_xy_xy_xy, g_xz_xy_xy_xz, g_xz_xy_xy_yy, g_xz_xy_xy_yz, g_xz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xy_xy_xx[i] = 4.0 * g_xz_xy_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_xy_xy[i] = 4.0 * g_xz_xy_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_xy_xz[i] = 4.0 * g_xz_xy_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_xy_yy[i] = 4.0 * g_xz_xy_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_xy_yz[i] = 4.0 * g_xz_xy_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_xy_zz[i] = 4.0 * g_xz_xy_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_xz_0_0_0_0_xy_xz_xx, g_xz_0_0_0_0_xy_xz_xy, g_xz_0_0_0_0_xy_xz_xz, g_xz_0_0_0_0_xy_xz_yy, g_xz_0_0_0_0_xy_xz_yz, g_xz_0_0_0_0_xy_xz_zz, g_xz_xy_xz_xx, g_xz_xy_xz_xy, g_xz_xy_xz_xz, g_xz_xy_xz_yy, g_xz_xy_xz_yz, g_xz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xy_xz_xx[i] = 4.0 * g_xz_xy_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_xz_xy[i] = 4.0 * g_xz_xy_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_xz_xz[i] = 4.0 * g_xz_xy_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_xz_yy[i] = 4.0 * g_xz_xy_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_xz_yz[i] = 4.0 * g_xz_xy_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_xz_zz[i] = 4.0 * g_xz_xy_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_xz_0_0_0_0_xy_yy_xx, g_xz_0_0_0_0_xy_yy_xy, g_xz_0_0_0_0_xy_yy_xz, g_xz_0_0_0_0_xy_yy_yy, g_xz_0_0_0_0_xy_yy_yz, g_xz_0_0_0_0_xy_yy_zz, g_xz_xy_yy_xx, g_xz_xy_yy_xy, g_xz_xy_yy_xz, g_xz_xy_yy_yy, g_xz_xy_yy_yz, g_xz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xy_yy_xx[i] = 4.0 * g_xz_xy_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_yy_xy[i] = 4.0 * g_xz_xy_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_yy_xz[i] = 4.0 * g_xz_xy_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_yy_yy[i] = 4.0 * g_xz_xy_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_yy_yz[i] = 4.0 * g_xz_xy_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_yy_zz[i] = 4.0 * g_xz_xy_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_xz_0_0_0_0_xy_yz_xx, g_xz_0_0_0_0_xy_yz_xy, g_xz_0_0_0_0_xy_yz_xz, g_xz_0_0_0_0_xy_yz_yy, g_xz_0_0_0_0_xy_yz_yz, g_xz_0_0_0_0_xy_yz_zz, g_xz_xy_yz_xx, g_xz_xy_yz_xy, g_xz_xy_yz_xz, g_xz_xy_yz_yy, g_xz_xy_yz_yz, g_xz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xy_yz_xx[i] = 4.0 * g_xz_xy_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_yz_xy[i] = 4.0 * g_xz_xy_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_yz_xz[i] = 4.0 * g_xz_xy_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_yz_yy[i] = 4.0 * g_xz_xy_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_yz_yz[i] = 4.0 * g_xz_xy_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_yz_zz[i] = 4.0 * g_xz_xy_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_xz_0_0_0_0_xy_zz_xx, g_xz_0_0_0_0_xy_zz_xy, g_xz_0_0_0_0_xy_zz_xz, g_xz_0_0_0_0_xy_zz_yy, g_xz_0_0_0_0_xy_zz_yz, g_xz_0_0_0_0_xy_zz_zz, g_xz_xy_zz_xx, g_xz_xy_zz_xy, g_xz_xy_zz_xz, g_xz_xy_zz_yy, g_xz_xy_zz_yz, g_xz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xy_zz_xx[i] = 4.0 * g_xz_xy_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_zz_xy[i] = 4.0 * g_xz_xy_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_zz_xz[i] = 4.0 * g_xz_xy_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_zz_yy[i] = 4.0 * g_xz_xy_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_zz_yz[i] = 4.0 * g_xz_xy_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_zz_zz[i] = 4.0 * g_xz_xy_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_xz_0_0_0_0_xz_xx_xx, g_xz_0_0_0_0_xz_xx_xy, g_xz_0_0_0_0_xz_xx_xz, g_xz_0_0_0_0_xz_xx_yy, g_xz_0_0_0_0_xz_xx_yz, g_xz_0_0_0_0_xz_xx_zz, g_xz_xz_xx_xx, g_xz_xz_xx_xy, g_xz_xz_xx_xz, g_xz_xz_xx_yy, g_xz_xz_xx_yz, g_xz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xz_xx_xx[i] = 4.0 * g_xz_xz_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_xx_xy[i] = 4.0 * g_xz_xz_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_xx_xz[i] = 4.0 * g_xz_xz_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_xx_yy[i] = 4.0 * g_xz_xz_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_xx_yz[i] = 4.0 * g_xz_xz_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_xx_zz[i] = 4.0 * g_xz_xz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_xz_0_0_0_0_xz_xy_xx, g_xz_0_0_0_0_xz_xy_xy, g_xz_0_0_0_0_xz_xy_xz, g_xz_0_0_0_0_xz_xy_yy, g_xz_0_0_0_0_xz_xy_yz, g_xz_0_0_0_0_xz_xy_zz, g_xz_xz_xy_xx, g_xz_xz_xy_xy, g_xz_xz_xy_xz, g_xz_xz_xy_yy, g_xz_xz_xy_yz, g_xz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xz_xy_xx[i] = 4.0 * g_xz_xz_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_xy_xy[i] = 4.0 * g_xz_xz_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_xy_xz[i] = 4.0 * g_xz_xz_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_xy_yy[i] = 4.0 * g_xz_xz_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_xy_yz[i] = 4.0 * g_xz_xz_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_xy_zz[i] = 4.0 * g_xz_xz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_xz_0_0_0_0_xz_xz_xx, g_xz_0_0_0_0_xz_xz_xy, g_xz_0_0_0_0_xz_xz_xz, g_xz_0_0_0_0_xz_xz_yy, g_xz_0_0_0_0_xz_xz_yz, g_xz_0_0_0_0_xz_xz_zz, g_xz_xz_xz_xx, g_xz_xz_xz_xy, g_xz_xz_xz_xz, g_xz_xz_xz_yy, g_xz_xz_xz_yz, g_xz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xz_xz_xx[i] = 4.0 * g_xz_xz_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_xz_xy[i] = 4.0 * g_xz_xz_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_xz_xz[i] = 4.0 * g_xz_xz_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_xz_yy[i] = 4.0 * g_xz_xz_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_xz_yz[i] = 4.0 * g_xz_xz_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_xz_zz[i] = 4.0 * g_xz_xz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_xz_0_0_0_0_xz_yy_xx, g_xz_0_0_0_0_xz_yy_xy, g_xz_0_0_0_0_xz_yy_xz, g_xz_0_0_0_0_xz_yy_yy, g_xz_0_0_0_0_xz_yy_yz, g_xz_0_0_0_0_xz_yy_zz, g_xz_xz_yy_xx, g_xz_xz_yy_xy, g_xz_xz_yy_xz, g_xz_xz_yy_yy, g_xz_xz_yy_yz, g_xz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xz_yy_xx[i] = 4.0 * g_xz_xz_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_yy_xy[i] = 4.0 * g_xz_xz_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_yy_xz[i] = 4.0 * g_xz_xz_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_yy_yy[i] = 4.0 * g_xz_xz_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_yy_yz[i] = 4.0 * g_xz_xz_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_yy_zz[i] = 4.0 * g_xz_xz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_xz_0_0_0_0_xz_yz_xx, g_xz_0_0_0_0_xz_yz_xy, g_xz_0_0_0_0_xz_yz_xz, g_xz_0_0_0_0_xz_yz_yy, g_xz_0_0_0_0_xz_yz_yz, g_xz_0_0_0_0_xz_yz_zz, g_xz_xz_yz_xx, g_xz_xz_yz_xy, g_xz_xz_yz_xz, g_xz_xz_yz_yy, g_xz_xz_yz_yz, g_xz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xz_yz_xx[i] = 4.0 * g_xz_xz_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_yz_xy[i] = 4.0 * g_xz_xz_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_yz_xz[i] = 4.0 * g_xz_xz_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_yz_yy[i] = 4.0 * g_xz_xz_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_yz_yz[i] = 4.0 * g_xz_xz_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_yz_zz[i] = 4.0 * g_xz_xz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_xz_0_0_0_0_xz_zz_xx, g_xz_0_0_0_0_xz_zz_xy, g_xz_0_0_0_0_xz_zz_xz, g_xz_0_0_0_0_xz_zz_yy, g_xz_0_0_0_0_xz_zz_yz, g_xz_0_0_0_0_xz_zz_zz, g_xz_xz_zz_xx, g_xz_xz_zz_xy, g_xz_xz_zz_xz, g_xz_xz_zz_yy, g_xz_xz_zz_yz, g_xz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xz_zz_xx[i] = 4.0 * g_xz_xz_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_zz_xy[i] = 4.0 * g_xz_xz_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_zz_xz[i] = 4.0 * g_xz_xz_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_zz_yy[i] = 4.0 * g_xz_xz_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_zz_yz[i] = 4.0 * g_xz_xz_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_zz_zz[i] = 4.0 * g_xz_xz_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_xz_0_0_0_0_yy_xx_xx, g_xz_0_0_0_0_yy_xx_xy, g_xz_0_0_0_0_yy_xx_xz, g_xz_0_0_0_0_yy_xx_yy, g_xz_0_0_0_0_yy_xx_yz, g_xz_0_0_0_0_yy_xx_zz, g_xz_yy_xx_xx, g_xz_yy_xx_xy, g_xz_yy_xx_xz, g_xz_yy_xx_yy, g_xz_yy_xx_yz, g_xz_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yy_xx_xx[i] = 4.0 * g_xz_yy_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_xx_xy[i] = 4.0 * g_xz_yy_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_xx_xz[i] = 4.0 * g_xz_yy_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_xx_yy[i] = 4.0 * g_xz_yy_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_xx_yz[i] = 4.0 * g_xz_yy_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_xx_zz[i] = 4.0 * g_xz_yy_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_xz_0_0_0_0_yy_xy_xx, g_xz_0_0_0_0_yy_xy_xy, g_xz_0_0_0_0_yy_xy_xz, g_xz_0_0_0_0_yy_xy_yy, g_xz_0_0_0_0_yy_xy_yz, g_xz_0_0_0_0_yy_xy_zz, g_xz_yy_xy_xx, g_xz_yy_xy_xy, g_xz_yy_xy_xz, g_xz_yy_xy_yy, g_xz_yy_xy_yz, g_xz_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yy_xy_xx[i] = 4.0 * g_xz_yy_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_xy_xy[i] = 4.0 * g_xz_yy_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_xy_xz[i] = 4.0 * g_xz_yy_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_xy_yy[i] = 4.0 * g_xz_yy_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_xy_yz[i] = 4.0 * g_xz_yy_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_xy_zz[i] = 4.0 * g_xz_yy_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_xz_0_0_0_0_yy_xz_xx, g_xz_0_0_0_0_yy_xz_xy, g_xz_0_0_0_0_yy_xz_xz, g_xz_0_0_0_0_yy_xz_yy, g_xz_0_0_0_0_yy_xz_yz, g_xz_0_0_0_0_yy_xz_zz, g_xz_yy_xz_xx, g_xz_yy_xz_xy, g_xz_yy_xz_xz, g_xz_yy_xz_yy, g_xz_yy_xz_yz, g_xz_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yy_xz_xx[i] = 4.0 * g_xz_yy_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_xz_xy[i] = 4.0 * g_xz_yy_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_xz_xz[i] = 4.0 * g_xz_yy_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_xz_yy[i] = 4.0 * g_xz_yy_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_xz_yz[i] = 4.0 * g_xz_yy_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_xz_zz[i] = 4.0 * g_xz_yy_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_xz_0_0_0_0_yy_yy_xx, g_xz_0_0_0_0_yy_yy_xy, g_xz_0_0_0_0_yy_yy_xz, g_xz_0_0_0_0_yy_yy_yy, g_xz_0_0_0_0_yy_yy_yz, g_xz_0_0_0_0_yy_yy_zz, g_xz_yy_yy_xx, g_xz_yy_yy_xy, g_xz_yy_yy_xz, g_xz_yy_yy_yy, g_xz_yy_yy_yz, g_xz_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yy_yy_xx[i] = 4.0 * g_xz_yy_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_yy_xy[i] = 4.0 * g_xz_yy_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_yy_xz[i] = 4.0 * g_xz_yy_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_yy_yy[i] = 4.0 * g_xz_yy_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_yy_yz[i] = 4.0 * g_xz_yy_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_yy_zz[i] = 4.0 * g_xz_yy_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_xz_0_0_0_0_yy_yz_xx, g_xz_0_0_0_0_yy_yz_xy, g_xz_0_0_0_0_yy_yz_xz, g_xz_0_0_0_0_yy_yz_yy, g_xz_0_0_0_0_yy_yz_yz, g_xz_0_0_0_0_yy_yz_zz, g_xz_yy_yz_xx, g_xz_yy_yz_xy, g_xz_yy_yz_xz, g_xz_yy_yz_yy, g_xz_yy_yz_yz, g_xz_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yy_yz_xx[i] = 4.0 * g_xz_yy_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_yz_xy[i] = 4.0 * g_xz_yy_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_yz_xz[i] = 4.0 * g_xz_yy_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_yz_yy[i] = 4.0 * g_xz_yy_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_yz_yz[i] = 4.0 * g_xz_yy_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_yz_zz[i] = 4.0 * g_xz_yy_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_xz_0_0_0_0_yy_zz_xx, g_xz_0_0_0_0_yy_zz_xy, g_xz_0_0_0_0_yy_zz_xz, g_xz_0_0_0_0_yy_zz_yy, g_xz_0_0_0_0_yy_zz_yz, g_xz_0_0_0_0_yy_zz_zz, g_xz_yy_zz_xx, g_xz_yy_zz_xy, g_xz_yy_zz_xz, g_xz_yy_zz_yy, g_xz_yy_zz_yz, g_xz_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yy_zz_xx[i] = 4.0 * g_xz_yy_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_zz_xy[i] = 4.0 * g_xz_yy_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_zz_xz[i] = 4.0 * g_xz_yy_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_zz_yy[i] = 4.0 * g_xz_yy_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_zz_yz[i] = 4.0 * g_xz_yy_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_zz_zz[i] = 4.0 * g_xz_yy_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_xz_0_0_0_0_yz_xx_xx, g_xz_0_0_0_0_yz_xx_xy, g_xz_0_0_0_0_yz_xx_xz, g_xz_0_0_0_0_yz_xx_yy, g_xz_0_0_0_0_yz_xx_yz, g_xz_0_0_0_0_yz_xx_zz, g_xz_yz_xx_xx, g_xz_yz_xx_xy, g_xz_yz_xx_xz, g_xz_yz_xx_yy, g_xz_yz_xx_yz, g_xz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yz_xx_xx[i] = 4.0 * g_xz_yz_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_xx_xy[i] = 4.0 * g_xz_yz_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_xx_xz[i] = 4.0 * g_xz_yz_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_xx_yy[i] = 4.0 * g_xz_yz_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_xx_yz[i] = 4.0 * g_xz_yz_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_xx_zz[i] = 4.0 * g_xz_yz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_xz_0_0_0_0_yz_xy_xx, g_xz_0_0_0_0_yz_xy_xy, g_xz_0_0_0_0_yz_xy_xz, g_xz_0_0_0_0_yz_xy_yy, g_xz_0_0_0_0_yz_xy_yz, g_xz_0_0_0_0_yz_xy_zz, g_xz_yz_xy_xx, g_xz_yz_xy_xy, g_xz_yz_xy_xz, g_xz_yz_xy_yy, g_xz_yz_xy_yz, g_xz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yz_xy_xx[i] = 4.0 * g_xz_yz_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_xy_xy[i] = 4.0 * g_xz_yz_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_xy_xz[i] = 4.0 * g_xz_yz_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_xy_yy[i] = 4.0 * g_xz_yz_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_xy_yz[i] = 4.0 * g_xz_yz_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_xy_zz[i] = 4.0 * g_xz_yz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_xz_0_0_0_0_yz_xz_xx, g_xz_0_0_0_0_yz_xz_xy, g_xz_0_0_0_0_yz_xz_xz, g_xz_0_0_0_0_yz_xz_yy, g_xz_0_0_0_0_yz_xz_yz, g_xz_0_0_0_0_yz_xz_zz, g_xz_yz_xz_xx, g_xz_yz_xz_xy, g_xz_yz_xz_xz, g_xz_yz_xz_yy, g_xz_yz_xz_yz, g_xz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yz_xz_xx[i] = 4.0 * g_xz_yz_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_xz_xy[i] = 4.0 * g_xz_yz_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_xz_xz[i] = 4.0 * g_xz_yz_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_xz_yy[i] = 4.0 * g_xz_yz_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_xz_yz[i] = 4.0 * g_xz_yz_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_xz_zz[i] = 4.0 * g_xz_yz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_xz_0_0_0_0_yz_yy_xx, g_xz_0_0_0_0_yz_yy_xy, g_xz_0_0_0_0_yz_yy_xz, g_xz_0_0_0_0_yz_yy_yy, g_xz_0_0_0_0_yz_yy_yz, g_xz_0_0_0_0_yz_yy_zz, g_xz_yz_yy_xx, g_xz_yz_yy_xy, g_xz_yz_yy_xz, g_xz_yz_yy_yy, g_xz_yz_yy_yz, g_xz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yz_yy_xx[i] = 4.0 * g_xz_yz_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_yy_xy[i] = 4.0 * g_xz_yz_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_yy_xz[i] = 4.0 * g_xz_yz_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_yy_yy[i] = 4.0 * g_xz_yz_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_yy_yz[i] = 4.0 * g_xz_yz_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_yy_zz[i] = 4.0 * g_xz_yz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_xz_0_0_0_0_yz_yz_xx, g_xz_0_0_0_0_yz_yz_xy, g_xz_0_0_0_0_yz_yz_xz, g_xz_0_0_0_0_yz_yz_yy, g_xz_0_0_0_0_yz_yz_yz, g_xz_0_0_0_0_yz_yz_zz, g_xz_yz_yz_xx, g_xz_yz_yz_xy, g_xz_yz_yz_xz, g_xz_yz_yz_yy, g_xz_yz_yz_yz, g_xz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yz_yz_xx[i] = 4.0 * g_xz_yz_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_yz_xy[i] = 4.0 * g_xz_yz_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_yz_xz[i] = 4.0 * g_xz_yz_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_yz_yy[i] = 4.0 * g_xz_yz_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_yz_yz[i] = 4.0 * g_xz_yz_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_yz_zz[i] = 4.0 * g_xz_yz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_xz_0_0_0_0_yz_zz_xx, g_xz_0_0_0_0_yz_zz_xy, g_xz_0_0_0_0_yz_zz_xz, g_xz_0_0_0_0_yz_zz_yy, g_xz_0_0_0_0_yz_zz_yz, g_xz_0_0_0_0_yz_zz_zz, g_xz_yz_zz_xx, g_xz_yz_zz_xy, g_xz_yz_zz_xz, g_xz_yz_zz_yy, g_xz_yz_zz_yz, g_xz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yz_zz_xx[i] = 4.0 * g_xz_yz_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_zz_xy[i] = 4.0 * g_xz_yz_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_zz_xz[i] = 4.0 * g_xz_yz_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_zz_yy[i] = 4.0 * g_xz_yz_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_zz_yz[i] = 4.0 * g_xz_yz_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_zz_zz[i] = 4.0 * g_xz_yz_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_xz_0_0_0_0_zz_xx_xx, g_xz_0_0_0_0_zz_xx_xy, g_xz_0_0_0_0_zz_xx_xz, g_xz_0_0_0_0_zz_xx_yy, g_xz_0_0_0_0_zz_xx_yz, g_xz_0_0_0_0_zz_xx_zz, g_xz_zz_xx_xx, g_xz_zz_xx_xy, g_xz_zz_xx_xz, g_xz_zz_xx_yy, g_xz_zz_xx_yz, g_xz_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_zz_xx_xx[i] = 4.0 * g_xz_zz_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_xx_xy[i] = 4.0 * g_xz_zz_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_xx_xz[i] = 4.0 * g_xz_zz_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_xx_yy[i] = 4.0 * g_xz_zz_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_xx_yz[i] = 4.0 * g_xz_zz_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_xx_zz[i] = 4.0 * g_xz_zz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_xz_0_0_0_0_zz_xy_xx, g_xz_0_0_0_0_zz_xy_xy, g_xz_0_0_0_0_zz_xy_xz, g_xz_0_0_0_0_zz_xy_yy, g_xz_0_0_0_0_zz_xy_yz, g_xz_0_0_0_0_zz_xy_zz, g_xz_zz_xy_xx, g_xz_zz_xy_xy, g_xz_zz_xy_xz, g_xz_zz_xy_yy, g_xz_zz_xy_yz, g_xz_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_zz_xy_xx[i] = 4.0 * g_xz_zz_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_xy_xy[i] = 4.0 * g_xz_zz_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_xy_xz[i] = 4.0 * g_xz_zz_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_xy_yy[i] = 4.0 * g_xz_zz_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_xy_yz[i] = 4.0 * g_xz_zz_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_xy_zz[i] = 4.0 * g_xz_zz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_xz_0_0_0_0_zz_xz_xx, g_xz_0_0_0_0_zz_xz_xy, g_xz_0_0_0_0_zz_xz_xz, g_xz_0_0_0_0_zz_xz_yy, g_xz_0_0_0_0_zz_xz_yz, g_xz_0_0_0_0_zz_xz_zz, g_xz_zz_xz_xx, g_xz_zz_xz_xy, g_xz_zz_xz_xz, g_xz_zz_xz_yy, g_xz_zz_xz_yz, g_xz_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_zz_xz_xx[i] = 4.0 * g_xz_zz_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_xz_xy[i] = 4.0 * g_xz_zz_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_xz_xz[i] = 4.0 * g_xz_zz_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_xz_yy[i] = 4.0 * g_xz_zz_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_xz_yz[i] = 4.0 * g_xz_zz_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_xz_zz[i] = 4.0 * g_xz_zz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_xz_0_0_0_0_zz_yy_xx, g_xz_0_0_0_0_zz_yy_xy, g_xz_0_0_0_0_zz_yy_xz, g_xz_0_0_0_0_zz_yy_yy, g_xz_0_0_0_0_zz_yy_yz, g_xz_0_0_0_0_zz_yy_zz, g_xz_zz_yy_xx, g_xz_zz_yy_xy, g_xz_zz_yy_xz, g_xz_zz_yy_yy, g_xz_zz_yy_yz, g_xz_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_zz_yy_xx[i] = 4.0 * g_xz_zz_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_yy_xy[i] = 4.0 * g_xz_zz_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_yy_xz[i] = 4.0 * g_xz_zz_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_yy_yy[i] = 4.0 * g_xz_zz_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_yy_yz[i] = 4.0 * g_xz_zz_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_yy_zz[i] = 4.0 * g_xz_zz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_xz_0_0_0_0_zz_yz_xx, g_xz_0_0_0_0_zz_yz_xy, g_xz_0_0_0_0_zz_yz_xz, g_xz_0_0_0_0_zz_yz_yy, g_xz_0_0_0_0_zz_yz_yz, g_xz_0_0_0_0_zz_yz_zz, g_xz_zz_yz_xx, g_xz_zz_yz_xy, g_xz_zz_yz_xz, g_xz_zz_yz_yy, g_xz_zz_yz_yz, g_xz_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_zz_yz_xx[i] = 4.0 * g_xz_zz_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_yz_xy[i] = 4.0 * g_xz_zz_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_yz_xz[i] = 4.0 * g_xz_zz_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_yz_yy[i] = 4.0 * g_xz_zz_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_yz_yz[i] = 4.0 * g_xz_zz_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_yz_zz[i] = 4.0 * g_xz_zz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_xz_0_0_0_0_zz_zz_xx, g_xz_0_0_0_0_zz_zz_xy, g_xz_0_0_0_0_zz_zz_xz, g_xz_0_0_0_0_zz_zz_yy, g_xz_0_0_0_0_zz_zz_yz, g_xz_0_0_0_0_zz_zz_zz, g_xz_zz_zz_xx, g_xz_zz_zz_xy, g_xz_zz_zz_xz, g_xz_zz_zz_yy, g_xz_zz_zz_yz, g_xz_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_zz_zz_xx[i] = 4.0 * g_xz_zz_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_zz_xy[i] = 4.0 * g_xz_zz_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_zz_xz[i] = 4.0 * g_xz_zz_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_zz_yy[i] = 4.0 * g_xz_zz_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_zz_yz[i] = 4.0 * g_xz_zz_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_zz_zz[i] = 4.0 * g_xz_zz_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_0_xx_xx_xx, g_0_xx_xx_xy, g_0_xx_xx_xz, g_0_xx_xx_yy, g_0_xx_xx_yz, g_0_xx_xx_zz, g_yy_0_0_0_0_xx_xx_xx, g_yy_0_0_0_0_xx_xx_xy, g_yy_0_0_0_0_xx_xx_xz, g_yy_0_0_0_0_xx_xx_yy, g_yy_0_0_0_0_xx_xx_yz, g_yy_0_0_0_0_xx_xx_zz, g_yy_xx_xx_xx, g_yy_xx_xx_xy, g_yy_xx_xx_xz, g_yy_xx_xx_yy, g_yy_xx_xx_yz, g_yy_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xx_xx_xx[i] = -2.0 * g_0_xx_xx_xx[i] * a_exp + 4.0 * g_yy_xx_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_xx_xy[i] = -2.0 * g_0_xx_xx_xy[i] * a_exp + 4.0 * g_yy_xx_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_xx_xz[i] = -2.0 * g_0_xx_xx_xz[i] * a_exp + 4.0 * g_yy_xx_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_xx_yy[i] = -2.0 * g_0_xx_xx_yy[i] * a_exp + 4.0 * g_yy_xx_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_xx_yz[i] = -2.0 * g_0_xx_xx_yz[i] * a_exp + 4.0 * g_yy_xx_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_xx_zz[i] = -2.0 * g_0_xx_xx_zz[i] * a_exp + 4.0 * g_yy_xx_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_0_xx_xy_xx, g_0_xx_xy_xy, g_0_xx_xy_xz, g_0_xx_xy_yy, g_0_xx_xy_yz, g_0_xx_xy_zz, g_yy_0_0_0_0_xx_xy_xx, g_yy_0_0_0_0_xx_xy_xy, g_yy_0_0_0_0_xx_xy_xz, g_yy_0_0_0_0_xx_xy_yy, g_yy_0_0_0_0_xx_xy_yz, g_yy_0_0_0_0_xx_xy_zz, g_yy_xx_xy_xx, g_yy_xx_xy_xy, g_yy_xx_xy_xz, g_yy_xx_xy_yy, g_yy_xx_xy_yz, g_yy_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xx_xy_xx[i] = -2.0 * g_0_xx_xy_xx[i] * a_exp + 4.0 * g_yy_xx_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_xy_xy[i] = -2.0 * g_0_xx_xy_xy[i] * a_exp + 4.0 * g_yy_xx_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_xy_xz[i] = -2.0 * g_0_xx_xy_xz[i] * a_exp + 4.0 * g_yy_xx_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_xy_yy[i] = -2.0 * g_0_xx_xy_yy[i] * a_exp + 4.0 * g_yy_xx_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_xy_yz[i] = -2.0 * g_0_xx_xy_yz[i] * a_exp + 4.0 * g_yy_xx_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_xy_zz[i] = -2.0 * g_0_xx_xy_zz[i] * a_exp + 4.0 * g_yy_xx_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_0_xx_xz_xx, g_0_xx_xz_xy, g_0_xx_xz_xz, g_0_xx_xz_yy, g_0_xx_xz_yz, g_0_xx_xz_zz, g_yy_0_0_0_0_xx_xz_xx, g_yy_0_0_0_0_xx_xz_xy, g_yy_0_0_0_0_xx_xz_xz, g_yy_0_0_0_0_xx_xz_yy, g_yy_0_0_0_0_xx_xz_yz, g_yy_0_0_0_0_xx_xz_zz, g_yy_xx_xz_xx, g_yy_xx_xz_xy, g_yy_xx_xz_xz, g_yy_xx_xz_yy, g_yy_xx_xz_yz, g_yy_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xx_xz_xx[i] = -2.0 * g_0_xx_xz_xx[i] * a_exp + 4.0 * g_yy_xx_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_xz_xy[i] = -2.0 * g_0_xx_xz_xy[i] * a_exp + 4.0 * g_yy_xx_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_xz_xz[i] = -2.0 * g_0_xx_xz_xz[i] * a_exp + 4.0 * g_yy_xx_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_xz_yy[i] = -2.0 * g_0_xx_xz_yy[i] * a_exp + 4.0 * g_yy_xx_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_xz_yz[i] = -2.0 * g_0_xx_xz_yz[i] * a_exp + 4.0 * g_yy_xx_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_xz_zz[i] = -2.0 * g_0_xx_xz_zz[i] * a_exp + 4.0 * g_yy_xx_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_0_xx_yy_xx, g_0_xx_yy_xy, g_0_xx_yy_xz, g_0_xx_yy_yy, g_0_xx_yy_yz, g_0_xx_yy_zz, g_yy_0_0_0_0_xx_yy_xx, g_yy_0_0_0_0_xx_yy_xy, g_yy_0_0_0_0_xx_yy_xz, g_yy_0_0_0_0_xx_yy_yy, g_yy_0_0_0_0_xx_yy_yz, g_yy_0_0_0_0_xx_yy_zz, g_yy_xx_yy_xx, g_yy_xx_yy_xy, g_yy_xx_yy_xz, g_yy_xx_yy_yy, g_yy_xx_yy_yz, g_yy_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xx_yy_xx[i] = -2.0 * g_0_xx_yy_xx[i] * a_exp + 4.0 * g_yy_xx_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_yy_xy[i] = -2.0 * g_0_xx_yy_xy[i] * a_exp + 4.0 * g_yy_xx_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_yy_xz[i] = -2.0 * g_0_xx_yy_xz[i] * a_exp + 4.0 * g_yy_xx_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_yy_yy[i] = -2.0 * g_0_xx_yy_yy[i] * a_exp + 4.0 * g_yy_xx_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_yy_yz[i] = -2.0 * g_0_xx_yy_yz[i] * a_exp + 4.0 * g_yy_xx_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_yy_zz[i] = -2.0 * g_0_xx_yy_zz[i] * a_exp + 4.0 * g_yy_xx_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_0_xx_yz_xx, g_0_xx_yz_xy, g_0_xx_yz_xz, g_0_xx_yz_yy, g_0_xx_yz_yz, g_0_xx_yz_zz, g_yy_0_0_0_0_xx_yz_xx, g_yy_0_0_0_0_xx_yz_xy, g_yy_0_0_0_0_xx_yz_xz, g_yy_0_0_0_0_xx_yz_yy, g_yy_0_0_0_0_xx_yz_yz, g_yy_0_0_0_0_xx_yz_zz, g_yy_xx_yz_xx, g_yy_xx_yz_xy, g_yy_xx_yz_xz, g_yy_xx_yz_yy, g_yy_xx_yz_yz, g_yy_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xx_yz_xx[i] = -2.0 * g_0_xx_yz_xx[i] * a_exp + 4.0 * g_yy_xx_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_yz_xy[i] = -2.0 * g_0_xx_yz_xy[i] * a_exp + 4.0 * g_yy_xx_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_yz_xz[i] = -2.0 * g_0_xx_yz_xz[i] * a_exp + 4.0 * g_yy_xx_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_yz_yy[i] = -2.0 * g_0_xx_yz_yy[i] * a_exp + 4.0 * g_yy_xx_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_yz_yz[i] = -2.0 * g_0_xx_yz_yz[i] * a_exp + 4.0 * g_yy_xx_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_yz_zz[i] = -2.0 * g_0_xx_yz_zz[i] * a_exp + 4.0 * g_yy_xx_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_0_xx_zz_xx, g_0_xx_zz_xy, g_0_xx_zz_xz, g_0_xx_zz_yy, g_0_xx_zz_yz, g_0_xx_zz_zz, g_yy_0_0_0_0_xx_zz_xx, g_yy_0_0_0_0_xx_zz_xy, g_yy_0_0_0_0_xx_zz_xz, g_yy_0_0_0_0_xx_zz_yy, g_yy_0_0_0_0_xx_zz_yz, g_yy_0_0_0_0_xx_zz_zz, g_yy_xx_zz_xx, g_yy_xx_zz_xy, g_yy_xx_zz_xz, g_yy_xx_zz_yy, g_yy_xx_zz_yz, g_yy_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xx_zz_xx[i] = -2.0 * g_0_xx_zz_xx[i] * a_exp + 4.0 * g_yy_xx_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_zz_xy[i] = -2.0 * g_0_xx_zz_xy[i] * a_exp + 4.0 * g_yy_xx_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_zz_xz[i] = -2.0 * g_0_xx_zz_xz[i] * a_exp + 4.0 * g_yy_xx_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_zz_yy[i] = -2.0 * g_0_xx_zz_yy[i] * a_exp + 4.0 * g_yy_xx_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_zz_yz[i] = -2.0 * g_0_xx_zz_yz[i] * a_exp + 4.0 * g_yy_xx_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_zz_zz[i] = -2.0 * g_0_xx_zz_zz[i] * a_exp + 4.0 * g_yy_xx_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_0_xy_xx_xx, g_0_xy_xx_xy, g_0_xy_xx_xz, g_0_xy_xx_yy, g_0_xy_xx_yz, g_0_xy_xx_zz, g_yy_0_0_0_0_xy_xx_xx, g_yy_0_0_0_0_xy_xx_xy, g_yy_0_0_0_0_xy_xx_xz, g_yy_0_0_0_0_xy_xx_yy, g_yy_0_0_0_0_xy_xx_yz, g_yy_0_0_0_0_xy_xx_zz, g_yy_xy_xx_xx, g_yy_xy_xx_xy, g_yy_xy_xx_xz, g_yy_xy_xx_yy, g_yy_xy_xx_yz, g_yy_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xy_xx_xx[i] = -2.0 * g_0_xy_xx_xx[i] * a_exp + 4.0 * g_yy_xy_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_xx_xy[i] = -2.0 * g_0_xy_xx_xy[i] * a_exp + 4.0 * g_yy_xy_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_xx_xz[i] = -2.0 * g_0_xy_xx_xz[i] * a_exp + 4.0 * g_yy_xy_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_xx_yy[i] = -2.0 * g_0_xy_xx_yy[i] * a_exp + 4.0 * g_yy_xy_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_xx_yz[i] = -2.0 * g_0_xy_xx_yz[i] * a_exp + 4.0 * g_yy_xy_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_xx_zz[i] = -2.0 * g_0_xy_xx_zz[i] * a_exp + 4.0 * g_yy_xy_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_0_xy_xy_xx, g_0_xy_xy_xy, g_0_xy_xy_xz, g_0_xy_xy_yy, g_0_xy_xy_yz, g_0_xy_xy_zz, g_yy_0_0_0_0_xy_xy_xx, g_yy_0_0_0_0_xy_xy_xy, g_yy_0_0_0_0_xy_xy_xz, g_yy_0_0_0_0_xy_xy_yy, g_yy_0_0_0_0_xy_xy_yz, g_yy_0_0_0_0_xy_xy_zz, g_yy_xy_xy_xx, g_yy_xy_xy_xy, g_yy_xy_xy_xz, g_yy_xy_xy_yy, g_yy_xy_xy_yz, g_yy_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xy_xy_xx[i] = -2.0 * g_0_xy_xy_xx[i] * a_exp + 4.0 * g_yy_xy_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_xy_xy[i] = -2.0 * g_0_xy_xy_xy[i] * a_exp + 4.0 * g_yy_xy_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_xy_xz[i] = -2.0 * g_0_xy_xy_xz[i] * a_exp + 4.0 * g_yy_xy_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_xy_yy[i] = -2.0 * g_0_xy_xy_yy[i] * a_exp + 4.0 * g_yy_xy_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_xy_yz[i] = -2.0 * g_0_xy_xy_yz[i] * a_exp + 4.0 * g_yy_xy_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_xy_zz[i] = -2.0 * g_0_xy_xy_zz[i] * a_exp + 4.0 * g_yy_xy_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_0_xy_xz_xx, g_0_xy_xz_xy, g_0_xy_xz_xz, g_0_xy_xz_yy, g_0_xy_xz_yz, g_0_xy_xz_zz, g_yy_0_0_0_0_xy_xz_xx, g_yy_0_0_0_0_xy_xz_xy, g_yy_0_0_0_0_xy_xz_xz, g_yy_0_0_0_0_xy_xz_yy, g_yy_0_0_0_0_xy_xz_yz, g_yy_0_0_0_0_xy_xz_zz, g_yy_xy_xz_xx, g_yy_xy_xz_xy, g_yy_xy_xz_xz, g_yy_xy_xz_yy, g_yy_xy_xz_yz, g_yy_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xy_xz_xx[i] = -2.0 * g_0_xy_xz_xx[i] * a_exp + 4.0 * g_yy_xy_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_xz_xy[i] = -2.0 * g_0_xy_xz_xy[i] * a_exp + 4.0 * g_yy_xy_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_xz_xz[i] = -2.0 * g_0_xy_xz_xz[i] * a_exp + 4.0 * g_yy_xy_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_xz_yy[i] = -2.0 * g_0_xy_xz_yy[i] * a_exp + 4.0 * g_yy_xy_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_xz_yz[i] = -2.0 * g_0_xy_xz_yz[i] * a_exp + 4.0 * g_yy_xy_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_xz_zz[i] = -2.0 * g_0_xy_xz_zz[i] * a_exp + 4.0 * g_yy_xy_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_0_xy_yy_xx, g_0_xy_yy_xy, g_0_xy_yy_xz, g_0_xy_yy_yy, g_0_xy_yy_yz, g_0_xy_yy_zz, g_yy_0_0_0_0_xy_yy_xx, g_yy_0_0_0_0_xy_yy_xy, g_yy_0_0_0_0_xy_yy_xz, g_yy_0_0_0_0_xy_yy_yy, g_yy_0_0_0_0_xy_yy_yz, g_yy_0_0_0_0_xy_yy_zz, g_yy_xy_yy_xx, g_yy_xy_yy_xy, g_yy_xy_yy_xz, g_yy_xy_yy_yy, g_yy_xy_yy_yz, g_yy_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xy_yy_xx[i] = -2.0 * g_0_xy_yy_xx[i] * a_exp + 4.0 * g_yy_xy_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_yy_xy[i] = -2.0 * g_0_xy_yy_xy[i] * a_exp + 4.0 * g_yy_xy_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_yy_xz[i] = -2.0 * g_0_xy_yy_xz[i] * a_exp + 4.0 * g_yy_xy_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_yy_yy[i] = -2.0 * g_0_xy_yy_yy[i] * a_exp + 4.0 * g_yy_xy_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_yy_yz[i] = -2.0 * g_0_xy_yy_yz[i] * a_exp + 4.0 * g_yy_xy_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_yy_zz[i] = -2.0 * g_0_xy_yy_zz[i] * a_exp + 4.0 * g_yy_xy_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_0_xy_yz_xx, g_0_xy_yz_xy, g_0_xy_yz_xz, g_0_xy_yz_yy, g_0_xy_yz_yz, g_0_xy_yz_zz, g_yy_0_0_0_0_xy_yz_xx, g_yy_0_0_0_0_xy_yz_xy, g_yy_0_0_0_0_xy_yz_xz, g_yy_0_0_0_0_xy_yz_yy, g_yy_0_0_0_0_xy_yz_yz, g_yy_0_0_0_0_xy_yz_zz, g_yy_xy_yz_xx, g_yy_xy_yz_xy, g_yy_xy_yz_xz, g_yy_xy_yz_yy, g_yy_xy_yz_yz, g_yy_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xy_yz_xx[i] = -2.0 * g_0_xy_yz_xx[i] * a_exp + 4.0 * g_yy_xy_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_yz_xy[i] = -2.0 * g_0_xy_yz_xy[i] * a_exp + 4.0 * g_yy_xy_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_yz_xz[i] = -2.0 * g_0_xy_yz_xz[i] * a_exp + 4.0 * g_yy_xy_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_yz_yy[i] = -2.0 * g_0_xy_yz_yy[i] * a_exp + 4.0 * g_yy_xy_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_yz_yz[i] = -2.0 * g_0_xy_yz_yz[i] * a_exp + 4.0 * g_yy_xy_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_yz_zz[i] = -2.0 * g_0_xy_yz_zz[i] * a_exp + 4.0 * g_yy_xy_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_0_xy_zz_xx, g_0_xy_zz_xy, g_0_xy_zz_xz, g_0_xy_zz_yy, g_0_xy_zz_yz, g_0_xy_zz_zz, g_yy_0_0_0_0_xy_zz_xx, g_yy_0_0_0_0_xy_zz_xy, g_yy_0_0_0_0_xy_zz_xz, g_yy_0_0_0_0_xy_zz_yy, g_yy_0_0_0_0_xy_zz_yz, g_yy_0_0_0_0_xy_zz_zz, g_yy_xy_zz_xx, g_yy_xy_zz_xy, g_yy_xy_zz_xz, g_yy_xy_zz_yy, g_yy_xy_zz_yz, g_yy_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xy_zz_xx[i] = -2.0 * g_0_xy_zz_xx[i] * a_exp + 4.0 * g_yy_xy_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_zz_xy[i] = -2.0 * g_0_xy_zz_xy[i] * a_exp + 4.0 * g_yy_xy_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_zz_xz[i] = -2.0 * g_0_xy_zz_xz[i] * a_exp + 4.0 * g_yy_xy_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_zz_yy[i] = -2.0 * g_0_xy_zz_yy[i] * a_exp + 4.0 * g_yy_xy_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_zz_yz[i] = -2.0 * g_0_xy_zz_yz[i] * a_exp + 4.0 * g_yy_xy_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_zz_zz[i] = -2.0 * g_0_xy_zz_zz[i] * a_exp + 4.0 * g_yy_xy_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_0_xz_xx_xx, g_0_xz_xx_xy, g_0_xz_xx_xz, g_0_xz_xx_yy, g_0_xz_xx_yz, g_0_xz_xx_zz, g_yy_0_0_0_0_xz_xx_xx, g_yy_0_0_0_0_xz_xx_xy, g_yy_0_0_0_0_xz_xx_xz, g_yy_0_0_0_0_xz_xx_yy, g_yy_0_0_0_0_xz_xx_yz, g_yy_0_0_0_0_xz_xx_zz, g_yy_xz_xx_xx, g_yy_xz_xx_xy, g_yy_xz_xx_xz, g_yy_xz_xx_yy, g_yy_xz_xx_yz, g_yy_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xz_xx_xx[i] = -2.0 * g_0_xz_xx_xx[i] * a_exp + 4.0 * g_yy_xz_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_xx_xy[i] = -2.0 * g_0_xz_xx_xy[i] * a_exp + 4.0 * g_yy_xz_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_xx_xz[i] = -2.0 * g_0_xz_xx_xz[i] * a_exp + 4.0 * g_yy_xz_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_xx_yy[i] = -2.0 * g_0_xz_xx_yy[i] * a_exp + 4.0 * g_yy_xz_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_xx_yz[i] = -2.0 * g_0_xz_xx_yz[i] * a_exp + 4.0 * g_yy_xz_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_xx_zz[i] = -2.0 * g_0_xz_xx_zz[i] * a_exp + 4.0 * g_yy_xz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_0_xz_xy_xx, g_0_xz_xy_xy, g_0_xz_xy_xz, g_0_xz_xy_yy, g_0_xz_xy_yz, g_0_xz_xy_zz, g_yy_0_0_0_0_xz_xy_xx, g_yy_0_0_0_0_xz_xy_xy, g_yy_0_0_0_0_xz_xy_xz, g_yy_0_0_0_0_xz_xy_yy, g_yy_0_0_0_0_xz_xy_yz, g_yy_0_0_0_0_xz_xy_zz, g_yy_xz_xy_xx, g_yy_xz_xy_xy, g_yy_xz_xy_xz, g_yy_xz_xy_yy, g_yy_xz_xy_yz, g_yy_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xz_xy_xx[i] = -2.0 * g_0_xz_xy_xx[i] * a_exp + 4.0 * g_yy_xz_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_xy_xy[i] = -2.0 * g_0_xz_xy_xy[i] * a_exp + 4.0 * g_yy_xz_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_xy_xz[i] = -2.0 * g_0_xz_xy_xz[i] * a_exp + 4.0 * g_yy_xz_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_xy_yy[i] = -2.0 * g_0_xz_xy_yy[i] * a_exp + 4.0 * g_yy_xz_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_xy_yz[i] = -2.0 * g_0_xz_xy_yz[i] * a_exp + 4.0 * g_yy_xz_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_xy_zz[i] = -2.0 * g_0_xz_xy_zz[i] * a_exp + 4.0 * g_yy_xz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_0_xz_xz_xx, g_0_xz_xz_xy, g_0_xz_xz_xz, g_0_xz_xz_yy, g_0_xz_xz_yz, g_0_xz_xz_zz, g_yy_0_0_0_0_xz_xz_xx, g_yy_0_0_0_0_xz_xz_xy, g_yy_0_0_0_0_xz_xz_xz, g_yy_0_0_0_0_xz_xz_yy, g_yy_0_0_0_0_xz_xz_yz, g_yy_0_0_0_0_xz_xz_zz, g_yy_xz_xz_xx, g_yy_xz_xz_xy, g_yy_xz_xz_xz, g_yy_xz_xz_yy, g_yy_xz_xz_yz, g_yy_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xz_xz_xx[i] = -2.0 * g_0_xz_xz_xx[i] * a_exp + 4.0 * g_yy_xz_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_xz_xy[i] = -2.0 * g_0_xz_xz_xy[i] * a_exp + 4.0 * g_yy_xz_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_xz_xz[i] = -2.0 * g_0_xz_xz_xz[i] * a_exp + 4.0 * g_yy_xz_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_xz_yy[i] = -2.0 * g_0_xz_xz_yy[i] * a_exp + 4.0 * g_yy_xz_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_xz_yz[i] = -2.0 * g_0_xz_xz_yz[i] * a_exp + 4.0 * g_yy_xz_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_xz_zz[i] = -2.0 * g_0_xz_xz_zz[i] * a_exp + 4.0 * g_yy_xz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_0_xz_yy_xx, g_0_xz_yy_xy, g_0_xz_yy_xz, g_0_xz_yy_yy, g_0_xz_yy_yz, g_0_xz_yy_zz, g_yy_0_0_0_0_xz_yy_xx, g_yy_0_0_0_0_xz_yy_xy, g_yy_0_0_0_0_xz_yy_xz, g_yy_0_0_0_0_xz_yy_yy, g_yy_0_0_0_0_xz_yy_yz, g_yy_0_0_0_0_xz_yy_zz, g_yy_xz_yy_xx, g_yy_xz_yy_xy, g_yy_xz_yy_xz, g_yy_xz_yy_yy, g_yy_xz_yy_yz, g_yy_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xz_yy_xx[i] = -2.0 * g_0_xz_yy_xx[i] * a_exp + 4.0 * g_yy_xz_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_yy_xy[i] = -2.0 * g_0_xz_yy_xy[i] * a_exp + 4.0 * g_yy_xz_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_yy_xz[i] = -2.0 * g_0_xz_yy_xz[i] * a_exp + 4.0 * g_yy_xz_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_yy_yy[i] = -2.0 * g_0_xz_yy_yy[i] * a_exp + 4.0 * g_yy_xz_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_yy_yz[i] = -2.0 * g_0_xz_yy_yz[i] * a_exp + 4.0 * g_yy_xz_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_yy_zz[i] = -2.0 * g_0_xz_yy_zz[i] * a_exp + 4.0 * g_yy_xz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_0_xz_yz_xx, g_0_xz_yz_xy, g_0_xz_yz_xz, g_0_xz_yz_yy, g_0_xz_yz_yz, g_0_xz_yz_zz, g_yy_0_0_0_0_xz_yz_xx, g_yy_0_0_0_0_xz_yz_xy, g_yy_0_0_0_0_xz_yz_xz, g_yy_0_0_0_0_xz_yz_yy, g_yy_0_0_0_0_xz_yz_yz, g_yy_0_0_0_0_xz_yz_zz, g_yy_xz_yz_xx, g_yy_xz_yz_xy, g_yy_xz_yz_xz, g_yy_xz_yz_yy, g_yy_xz_yz_yz, g_yy_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xz_yz_xx[i] = -2.0 * g_0_xz_yz_xx[i] * a_exp + 4.0 * g_yy_xz_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_yz_xy[i] = -2.0 * g_0_xz_yz_xy[i] * a_exp + 4.0 * g_yy_xz_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_yz_xz[i] = -2.0 * g_0_xz_yz_xz[i] * a_exp + 4.0 * g_yy_xz_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_yz_yy[i] = -2.0 * g_0_xz_yz_yy[i] * a_exp + 4.0 * g_yy_xz_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_yz_yz[i] = -2.0 * g_0_xz_yz_yz[i] * a_exp + 4.0 * g_yy_xz_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_yz_zz[i] = -2.0 * g_0_xz_yz_zz[i] * a_exp + 4.0 * g_yy_xz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_0_xz_zz_xx, g_0_xz_zz_xy, g_0_xz_zz_xz, g_0_xz_zz_yy, g_0_xz_zz_yz, g_0_xz_zz_zz, g_yy_0_0_0_0_xz_zz_xx, g_yy_0_0_0_0_xz_zz_xy, g_yy_0_0_0_0_xz_zz_xz, g_yy_0_0_0_0_xz_zz_yy, g_yy_0_0_0_0_xz_zz_yz, g_yy_0_0_0_0_xz_zz_zz, g_yy_xz_zz_xx, g_yy_xz_zz_xy, g_yy_xz_zz_xz, g_yy_xz_zz_yy, g_yy_xz_zz_yz, g_yy_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xz_zz_xx[i] = -2.0 * g_0_xz_zz_xx[i] * a_exp + 4.0 * g_yy_xz_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_zz_xy[i] = -2.0 * g_0_xz_zz_xy[i] * a_exp + 4.0 * g_yy_xz_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_zz_xz[i] = -2.0 * g_0_xz_zz_xz[i] * a_exp + 4.0 * g_yy_xz_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_zz_yy[i] = -2.0 * g_0_xz_zz_yy[i] * a_exp + 4.0 * g_yy_xz_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_zz_yz[i] = -2.0 * g_0_xz_zz_yz[i] * a_exp + 4.0 * g_yy_xz_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_zz_zz[i] = -2.0 * g_0_xz_zz_zz[i] * a_exp + 4.0 * g_yy_xz_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_0_yy_xx_xx, g_0_yy_xx_xy, g_0_yy_xx_xz, g_0_yy_xx_yy, g_0_yy_xx_yz, g_0_yy_xx_zz, g_yy_0_0_0_0_yy_xx_xx, g_yy_0_0_0_0_yy_xx_xy, g_yy_0_0_0_0_yy_xx_xz, g_yy_0_0_0_0_yy_xx_yy, g_yy_0_0_0_0_yy_xx_yz, g_yy_0_0_0_0_yy_xx_zz, g_yy_yy_xx_xx, g_yy_yy_xx_xy, g_yy_yy_xx_xz, g_yy_yy_xx_yy, g_yy_yy_xx_yz, g_yy_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yy_xx_xx[i] = -2.0 * g_0_yy_xx_xx[i] * a_exp + 4.0 * g_yy_yy_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_xx_xy[i] = -2.0 * g_0_yy_xx_xy[i] * a_exp + 4.0 * g_yy_yy_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_xx_xz[i] = -2.0 * g_0_yy_xx_xz[i] * a_exp + 4.0 * g_yy_yy_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_xx_yy[i] = -2.0 * g_0_yy_xx_yy[i] * a_exp + 4.0 * g_yy_yy_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_xx_yz[i] = -2.0 * g_0_yy_xx_yz[i] * a_exp + 4.0 * g_yy_yy_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_xx_zz[i] = -2.0 * g_0_yy_xx_zz[i] * a_exp + 4.0 * g_yy_yy_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_0_yy_xy_xx, g_0_yy_xy_xy, g_0_yy_xy_xz, g_0_yy_xy_yy, g_0_yy_xy_yz, g_0_yy_xy_zz, g_yy_0_0_0_0_yy_xy_xx, g_yy_0_0_0_0_yy_xy_xy, g_yy_0_0_0_0_yy_xy_xz, g_yy_0_0_0_0_yy_xy_yy, g_yy_0_0_0_0_yy_xy_yz, g_yy_0_0_0_0_yy_xy_zz, g_yy_yy_xy_xx, g_yy_yy_xy_xy, g_yy_yy_xy_xz, g_yy_yy_xy_yy, g_yy_yy_xy_yz, g_yy_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yy_xy_xx[i] = -2.0 * g_0_yy_xy_xx[i] * a_exp + 4.0 * g_yy_yy_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_xy_xy[i] = -2.0 * g_0_yy_xy_xy[i] * a_exp + 4.0 * g_yy_yy_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_xy_xz[i] = -2.0 * g_0_yy_xy_xz[i] * a_exp + 4.0 * g_yy_yy_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_xy_yy[i] = -2.0 * g_0_yy_xy_yy[i] * a_exp + 4.0 * g_yy_yy_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_xy_yz[i] = -2.0 * g_0_yy_xy_yz[i] * a_exp + 4.0 * g_yy_yy_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_xy_zz[i] = -2.0 * g_0_yy_xy_zz[i] * a_exp + 4.0 * g_yy_yy_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_0_yy_xz_xx, g_0_yy_xz_xy, g_0_yy_xz_xz, g_0_yy_xz_yy, g_0_yy_xz_yz, g_0_yy_xz_zz, g_yy_0_0_0_0_yy_xz_xx, g_yy_0_0_0_0_yy_xz_xy, g_yy_0_0_0_0_yy_xz_xz, g_yy_0_0_0_0_yy_xz_yy, g_yy_0_0_0_0_yy_xz_yz, g_yy_0_0_0_0_yy_xz_zz, g_yy_yy_xz_xx, g_yy_yy_xz_xy, g_yy_yy_xz_xz, g_yy_yy_xz_yy, g_yy_yy_xz_yz, g_yy_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yy_xz_xx[i] = -2.0 * g_0_yy_xz_xx[i] * a_exp + 4.0 * g_yy_yy_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_xz_xy[i] = -2.0 * g_0_yy_xz_xy[i] * a_exp + 4.0 * g_yy_yy_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_xz_xz[i] = -2.0 * g_0_yy_xz_xz[i] * a_exp + 4.0 * g_yy_yy_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_xz_yy[i] = -2.0 * g_0_yy_xz_yy[i] * a_exp + 4.0 * g_yy_yy_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_xz_yz[i] = -2.0 * g_0_yy_xz_yz[i] * a_exp + 4.0 * g_yy_yy_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_xz_zz[i] = -2.0 * g_0_yy_xz_zz[i] * a_exp + 4.0 * g_yy_yy_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_0_yy_yy_xx, g_0_yy_yy_xy, g_0_yy_yy_xz, g_0_yy_yy_yy, g_0_yy_yy_yz, g_0_yy_yy_zz, g_yy_0_0_0_0_yy_yy_xx, g_yy_0_0_0_0_yy_yy_xy, g_yy_0_0_0_0_yy_yy_xz, g_yy_0_0_0_0_yy_yy_yy, g_yy_0_0_0_0_yy_yy_yz, g_yy_0_0_0_0_yy_yy_zz, g_yy_yy_yy_xx, g_yy_yy_yy_xy, g_yy_yy_yy_xz, g_yy_yy_yy_yy, g_yy_yy_yy_yz, g_yy_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yy_yy_xx[i] = -2.0 * g_0_yy_yy_xx[i] * a_exp + 4.0 * g_yy_yy_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_yy_xy[i] = -2.0 * g_0_yy_yy_xy[i] * a_exp + 4.0 * g_yy_yy_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_yy_xz[i] = -2.0 * g_0_yy_yy_xz[i] * a_exp + 4.0 * g_yy_yy_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_yy_yy[i] = -2.0 * g_0_yy_yy_yy[i] * a_exp + 4.0 * g_yy_yy_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_yy_yz[i] = -2.0 * g_0_yy_yy_yz[i] * a_exp + 4.0 * g_yy_yy_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_yy_zz[i] = -2.0 * g_0_yy_yy_zz[i] * a_exp + 4.0 * g_yy_yy_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_0_yy_yz_xx, g_0_yy_yz_xy, g_0_yy_yz_xz, g_0_yy_yz_yy, g_0_yy_yz_yz, g_0_yy_yz_zz, g_yy_0_0_0_0_yy_yz_xx, g_yy_0_0_0_0_yy_yz_xy, g_yy_0_0_0_0_yy_yz_xz, g_yy_0_0_0_0_yy_yz_yy, g_yy_0_0_0_0_yy_yz_yz, g_yy_0_0_0_0_yy_yz_zz, g_yy_yy_yz_xx, g_yy_yy_yz_xy, g_yy_yy_yz_xz, g_yy_yy_yz_yy, g_yy_yy_yz_yz, g_yy_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yy_yz_xx[i] = -2.0 * g_0_yy_yz_xx[i] * a_exp + 4.0 * g_yy_yy_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_yz_xy[i] = -2.0 * g_0_yy_yz_xy[i] * a_exp + 4.0 * g_yy_yy_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_yz_xz[i] = -2.0 * g_0_yy_yz_xz[i] * a_exp + 4.0 * g_yy_yy_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_yz_yy[i] = -2.0 * g_0_yy_yz_yy[i] * a_exp + 4.0 * g_yy_yy_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_yz_yz[i] = -2.0 * g_0_yy_yz_yz[i] * a_exp + 4.0 * g_yy_yy_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_yz_zz[i] = -2.0 * g_0_yy_yz_zz[i] * a_exp + 4.0 * g_yy_yy_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_0_yy_zz_xx, g_0_yy_zz_xy, g_0_yy_zz_xz, g_0_yy_zz_yy, g_0_yy_zz_yz, g_0_yy_zz_zz, g_yy_0_0_0_0_yy_zz_xx, g_yy_0_0_0_0_yy_zz_xy, g_yy_0_0_0_0_yy_zz_xz, g_yy_0_0_0_0_yy_zz_yy, g_yy_0_0_0_0_yy_zz_yz, g_yy_0_0_0_0_yy_zz_zz, g_yy_yy_zz_xx, g_yy_yy_zz_xy, g_yy_yy_zz_xz, g_yy_yy_zz_yy, g_yy_yy_zz_yz, g_yy_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yy_zz_xx[i] = -2.0 * g_0_yy_zz_xx[i] * a_exp + 4.0 * g_yy_yy_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_zz_xy[i] = -2.0 * g_0_yy_zz_xy[i] * a_exp + 4.0 * g_yy_yy_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_zz_xz[i] = -2.0 * g_0_yy_zz_xz[i] * a_exp + 4.0 * g_yy_yy_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_zz_yy[i] = -2.0 * g_0_yy_zz_yy[i] * a_exp + 4.0 * g_yy_yy_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_zz_yz[i] = -2.0 * g_0_yy_zz_yz[i] * a_exp + 4.0 * g_yy_yy_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_zz_zz[i] = -2.0 * g_0_yy_zz_zz[i] * a_exp + 4.0 * g_yy_yy_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_0_yz_xx_xx, g_0_yz_xx_xy, g_0_yz_xx_xz, g_0_yz_xx_yy, g_0_yz_xx_yz, g_0_yz_xx_zz, g_yy_0_0_0_0_yz_xx_xx, g_yy_0_0_0_0_yz_xx_xy, g_yy_0_0_0_0_yz_xx_xz, g_yy_0_0_0_0_yz_xx_yy, g_yy_0_0_0_0_yz_xx_yz, g_yy_0_0_0_0_yz_xx_zz, g_yy_yz_xx_xx, g_yy_yz_xx_xy, g_yy_yz_xx_xz, g_yy_yz_xx_yy, g_yy_yz_xx_yz, g_yy_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yz_xx_xx[i] = -2.0 * g_0_yz_xx_xx[i] * a_exp + 4.0 * g_yy_yz_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_xx_xy[i] = -2.0 * g_0_yz_xx_xy[i] * a_exp + 4.0 * g_yy_yz_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_xx_xz[i] = -2.0 * g_0_yz_xx_xz[i] * a_exp + 4.0 * g_yy_yz_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_xx_yy[i] = -2.0 * g_0_yz_xx_yy[i] * a_exp + 4.0 * g_yy_yz_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_xx_yz[i] = -2.0 * g_0_yz_xx_yz[i] * a_exp + 4.0 * g_yy_yz_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_xx_zz[i] = -2.0 * g_0_yz_xx_zz[i] * a_exp + 4.0 * g_yy_yz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_0_yz_xy_xx, g_0_yz_xy_xy, g_0_yz_xy_xz, g_0_yz_xy_yy, g_0_yz_xy_yz, g_0_yz_xy_zz, g_yy_0_0_0_0_yz_xy_xx, g_yy_0_0_0_0_yz_xy_xy, g_yy_0_0_0_0_yz_xy_xz, g_yy_0_0_0_0_yz_xy_yy, g_yy_0_0_0_0_yz_xy_yz, g_yy_0_0_0_0_yz_xy_zz, g_yy_yz_xy_xx, g_yy_yz_xy_xy, g_yy_yz_xy_xz, g_yy_yz_xy_yy, g_yy_yz_xy_yz, g_yy_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yz_xy_xx[i] = -2.0 * g_0_yz_xy_xx[i] * a_exp + 4.0 * g_yy_yz_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_xy_xy[i] = -2.0 * g_0_yz_xy_xy[i] * a_exp + 4.0 * g_yy_yz_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_xy_xz[i] = -2.0 * g_0_yz_xy_xz[i] * a_exp + 4.0 * g_yy_yz_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_xy_yy[i] = -2.0 * g_0_yz_xy_yy[i] * a_exp + 4.0 * g_yy_yz_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_xy_yz[i] = -2.0 * g_0_yz_xy_yz[i] * a_exp + 4.0 * g_yy_yz_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_xy_zz[i] = -2.0 * g_0_yz_xy_zz[i] * a_exp + 4.0 * g_yy_yz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_0_yz_xz_xx, g_0_yz_xz_xy, g_0_yz_xz_xz, g_0_yz_xz_yy, g_0_yz_xz_yz, g_0_yz_xz_zz, g_yy_0_0_0_0_yz_xz_xx, g_yy_0_0_0_0_yz_xz_xy, g_yy_0_0_0_0_yz_xz_xz, g_yy_0_0_0_0_yz_xz_yy, g_yy_0_0_0_0_yz_xz_yz, g_yy_0_0_0_0_yz_xz_zz, g_yy_yz_xz_xx, g_yy_yz_xz_xy, g_yy_yz_xz_xz, g_yy_yz_xz_yy, g_yy_yz_xz_yz, g_yy_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yz_xz_xx[i] = -2.0 * g_0_yz_xz_xx[i] * a_exp + 4.0 * g_yy_yz_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_xz_xy[i] = -2.0 * g_0_yz_xz_xy[i] * a_exp + 4.0 * g_yy_yz_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_xz_xz[i] = -2.0 * g_0_yz_xz_xz[i] * a_exp + 4.0 * g_yy_yz_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_xz_yy[i] = -2.0 * g_0_yz_xz_yy[i] * a_exp + 4.0 * g_yy_yz_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_xz_yz[i] = -2.0 * g_0_yz_xz_yz[i] * a_exp + 4.0 * g_yy_yz_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_xz_zz[i] = -2.0 * g_0_yz_xz_zz[i] * a_exp + 4.0 * g_yy_yz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_0_yz_yy_xx, g_0_yz_yy_xy, g_0_yz_yy_xz, g_0_yz_yy_yy, g_0_yz_yy_yz, g_0_yz_yy_zz, g_yy_0_0_0_0_yz_yy_xx, g_yy_0_0_0_0_yz_yy_xy, g_yy_0_0_0_0_yz_yy_xz, g_yy_0_0_0_0_yz_yy_yy, g_yy_0_0_0_0_yz_yy_yz, g_yy_0_0_0_0_yz_yy_zz, g_yy_yz_yy_xx, g_yy_yz_yy_xy, g_yy_yz_yy_xz, g_yy_yz_yy_yy, g_yy_yz_yy_yz, g_yy_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yz_yy_xx[i] = -2.0 * g_0_yz_yy_xx[i] * a_exp + 4.0 * g_yy_yz_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_yy_xy[i] = -2.0 * g_0_yz_yy_xy[i] * a_exp + 4.0 * g_yy_yz_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_yy_xz[i] = -2.0 * g_0_yz_yy_xz[i] * a_exp + 4.0 * g_yy_yz_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_yy_yy[i] = -2.0 * g_0_yz_yy_yy[i] * a_exp + 4.0 * g_yy_yz_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_yy_yz[i] = -2.0 * g_0_yz_yy_yz[i] * a_exp + 4.0 * g_yy_yz_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_yy_zz[i] = -2.0 * g_0_yz_yy_zz[i] * a_exp + 4.0 * g_yy_yz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_0_yz_yz_xx, g_0_yz_yz_xy, g_0_yz_yz_xz, g_0_yz_yz_yy, g_0_yz_yz_yz, g_0_yz_yz_zz, g_yy_0_0_0_0_yz_yz_xx, g_yy_0_0_0_0_yz_yz_xy, g_yy_0_0_0_0_yz_yz_xz, g_yy_0_0_0_0_yz_yz_yy, g_yy_0_0_0_0_yz_yz_yz, g_yy_0_0_0_0_yz_yz_zz, g_yy_yz_yz_xx, g_yy_yz_yz_xy, g_yy_yz_yz_xz, g_yy_yz_yz_yy, g_yy_yz_yz_yz, g_yy_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yz_yz_xx[i] = -2.0 * g_0_yz_yz_xx[i] * a_exp + 4.0 * g_yy_yz_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_yz_xy[i] = -2.0 * g_0_yz_yz_xy[i] * a_exp + 4.0 * g_yy_yz_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_yz_xz[i] = -2.0 * g_0_yz_yz_xz[i] * a_exp + 4.0 * g_yy_yz_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_yz_yy[i] = -2.0 * g_0_yz_yz_yy[i] * a_exp + 4.0 * g_yy_yz_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_yz_yz[i] = -2.0 * g_0_yz_yz_yz[i] * a_exp + 4.0 * g_yy_yz_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_yz_zz[i] = -2.0 * g_0_yz_yz_zz[i] * a_exp + 4.0 * g_yy_yz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_0_yz_zz_xx, g_0_yz_zz_xy, g_0_yz_zz_xz, g_0_yz_zz_yy, g_0_yz_zz_yz, g_0_yz_zz_zz, g_yy_0_0_0_0_yz_zz_xx, g_yy_0_0_0_0_yz_zz_xy, g_yy_0_0_0_0_yz_zz_xz, g_yy_0_0_0_0_yz_zz_yy, g_yy_0_0_0_0_yz_zz_yz, g_yy_0_0_0_0_yz_zz_zz, g_yy_yz_zz_xx, g_yy_yz_zz_xy, g_yy_yz_zz_xz, g_yy_yz_zz_yy, g_yy_yz_zz_yz, g_yy_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yz_zz_xx[i] = -2.0 * g_0_yz_zz_xx[i] * a_exp + 4.0 * g_yy_yz_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_zz_xy[i] = -2.0 * g_0_yz_zz_xy[i] * a_exp + 4.0 * g_yy_yz_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_zz_xz[i] = -2.0 * g_0_yz_zz_xz[i] * a_exp + 4.0 * g_yy_yz_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_zz_yy[i] = -2.0 * g_0_yz_zz_yy[i] * a_exp + 4.0 * g_yy_yz_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_zz_yz[i] = -2.0 * g_0_yz_zz_yz[i] * a_exp + 4.0 * g_yy_yz_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_zz_zz[i] = -2.0 * g_0_yz_zz_zz[i] * a_exp + 4.0 * g_yy_yz_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_0_zz_xx_xx, g_0_zz_xx_xy, g_0_zz_xx_xz, g_0_zz_xx_yy, g_0_zz_xx_yz, g_0_zz_xx_zz, g_yy_0_0_0_0_zz_xx_xx, g_yy_0_0_0_0_zz_xx_xy, g_yy_0_0_0_0_zz_xx_xz, g_yy_0_0_0_0_zz_xx_yy, g_yy_0_0_0_0_zz_xx_yz, g_yy_0_0_0_0_zz_xx_zz, g_yy_zz_xx_xx, g_yy_zz_xx_xy, g_yy_zz_xx_xz, g_yy_zz_xx_yy, g_yy_zz_xx_yz, g_yy_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_zz_xx_xx[i] = -2.0 * g_0_zz_xx_xx[i] * a_exp + 4.0 * g_yy_zz_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_xx_xy[i] = -2.0 * g_0_zz_xx_xy[i] * a_exp + 4.0 * g_yy_zz_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_xx_xz[i] = -2.0 * g_0_zz_xx_xz[i] * a_exp + 4.0 * g_yy_zz_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_xx_yy[i] = -2.0 * g_0_zz_xx_yy[i] * a_exp + 4.0 * g_yy_zz_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_xx_yz[i] = -2.0 * g_0_zz_xx_yz[i] * a_exp + 4.0 * g_yy_zz_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_xx_zz[i] = -2.0 * g_0_zz_xx_zz[i] * a_exp + 4.0 * g_yy_zz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_0_zz_xy_xx, g_0_zz_xy_xy, g_0_zz_xy_xz, g_0_zz_xy_yy, g_0_zz_xy_yz, g_0_zz_xy_zz, g_yy_0_0_0_0_zz_xy_xx, g_yy_0_0_0_0_zz_xy_xy, g_yy_0_0_0_0_zz_xy_xz, g_yy_0_0_0_0_zz_xy_yy, g_yy_0_0_0_0_zz_xy_yz, g_yy_0_0_0_0_zz_xy_zz, g_yy_zz_xy_xx, g_yy_zz_xy_xy, g_yy_zz_xy_xz, g_yy_zz_xy_yy, g_yy_zz_xy_yz, g_yy_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_zz_xy_xx[i] = -2.0 * g_0_zz_xy_xx[i] * a_exp + 4.0 * g_yy_zz_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_xy_xy[i] = -2.0 * g_0_zz_xy_xy[i] * a_exp + 4.0 * g_yy_zz_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_xy_xz[i] = -2.0 * g_0_zz_xy_xz[i] * a_exp + 4.0 * g_yy_zz_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_xy_yy[i] = -2.0 * g_0_zz_xy_yy[i] * a_exp + 4.0 * g_yy_zz_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_xy_yz[i] = -2.0 * g_0_zz_xy_yz[i] * a_exp + 4.0 * g_yy_zz_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_xy_zz[i] = -2.0 * g_0_zz_xy_zz[i] * a_exp + 4.0 * g_yy_zz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_0_zz_xz_xx, g_0_zz_xz_xy, g_0_zz_xz_xz, g_0_zz_xz_yy, g_0_zz_xz_yz, g_0_zz_xz_zz, g_yy_0_0_0_0_zz_xz_xx, g_yy_0_0_0_0_zz_xz_xy, g_yy_0_0_0_0_zz_xz_xz, g_yy_0_0_0_0_zz_xz_yy, g_yy_0_0_0_0_zz_xz_yz, g_yy_0_0_0_0_zz_xz_zz, g_yy_zz_xz_xx, g_yy_zz_xz_xy, g_yy_zz_xz_xz, g_yy_zz_xz_yy, g_yy_zz_xz_yz, g_yy_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_zz_xz_xx[i] = -2.0 * g_0_zz_xz_xx[i] * a_exp + 4.0 * g_yy_zz_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_xz_xy[i] = -2.0 * g_0_zz_xz_xy[i] * a_exp + 4.0 * g_yy_zz_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_xz_xz[i] = -2.0 * g_0_zz_xz_xz[i] * a_exp + 4.0 * g_yy_zz_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_xz_yy[i] = -2.0 * g_0_zz_xz_yy[i] * a_exp + 4.0 * g_yy_zz_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_xz_yz[i] = -2.0 * g_0_zz_xz_yz[i] * a_exp + 4.0 * g_yy_zz_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_xz_zz[i] = -2.0 * g_0_zz_xz_zz[i] * a_exp + 4.0 * g_yy_zz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_0_zz_yy_xx, g_0_zz_yy_xy, g_0_zz_yy_xz, g_0_zz_yy_yy, g_0_zz_yy_yz, g_0_zz_yy_zz, g_yy_0_0_0_0_zz_yy_xx, g_yy_0_0_0_0_zz_yy_xy, g_yy_0_0_0_0_zz_yy_xz, g_yy_0_0_0_0_zz_yy_yy, g_yy_0_0_0_0_zz_yy_yz, g_yy_0_0_0_0_zz_yy_zz, g_yy_zz_yy_xx, g_yy_zz_yy_xy, g_yy_zz_yy_xz, g_yy_zz_yy_yy, g_yy_zz_yy_yz, g_yy_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_zz_yy_xx[i] = -2.0 * g_0_zz_yy_xx[i] * a_exp + 4.0 * g_yy_zz_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_yy_xy[i] = -2.0 * g_0_zz_yy_xy[i] * a_exp + 4.0 * g_yy_zz_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_yy_xz[i] = -2.0 * g_0_zz_yy_xz[i] * a_exp + 4.0 * g_yy_zz_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_yy_yy[i] = -2.0 * g_0_zz_yy_yy[i] * a_exp + 4.0 * g_yy_zz_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_yy_yz[i] = -2.0 * g_0_zz_yy_yz[i] * a_exp + 4.0 * g_yy_zz_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_yy_zz[i] = -2.0 * g_0_zz_yy_zz[i] * a_exp + 4.0 * g_yy_zz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_0_zz_yz_xx, g_0_zz_yz_xy, g_0_zz_yz_xz, g_0_zz_yz_yy, g_0_zz_yz_yz, g_0_zz_yz_zz, g_yy_0_0_0_0_zz_yz_xx, g_yy_0_0_0_0_zz_yz_xy, g_yy_0_0_0_0_zz_yz_xz, g_yy_0_0_0_0_zz_yz_yy, g_yy_0_0_0_0_zz_yz_yz, g_yy_0_0_0_0_zz_yz_zz, g_yy_zz_yz_xx, g_yy_zz_yz_xy, g_yy_zz_yz_xz, g_yy_zz_yz_yy, g_yy_zz_yz_yz, g_yy_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_zz_yz_xx[i] = -2.0 * g_0_zz_yz_xx[i] * a_exp + 4.0 * g_yy_zz_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_yz_xy[i] = -2.0 * g_0_zz_yz_xy[i] * a_exp + 4.0 * g_yy_zz_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_yz_xz[i] = -2.0 * g_0_zz_yz_xz[i] * a_exp + 4.0 * g_yy_zz_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_yz_yy[i] = -2.0 * g_0_zz_yz_yy[i] * a_exp + 4.0 * g_yy_zz_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_yz_yz[i] = -2.0 * g_0_zz_yz_yz[i] * a_exp + 4.0 * g_yy_zz_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_yz_zz[i] = -2.0 * g_0_zz_yz_zz[i] * a_exp + 4.0 * g_yy_zz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_0_zz_zz_xx, g_0_zz_zz_xy, g_0_zz_zz_xz, g_0_zz_zz_yy, g_0_zz_zz_yz, g_0_zz_zz_zz, g_yy_0_0_0_0_zz_zz_xx, g_yy_0_0_0_0_zz_zz_xy, g_yy_0_0_0_0_zz_zz_xz, g_yy_0_0_0_0_zz_zz_yy, g_yy_0_0_0_0_zz_zz_yz, g_yy_0_0_0_0_zz_zz_zz, g_yy_zz_zz_xx, g_yy_zz_zz_xy, g_yy_zz_zz_xz, g_yy_zz_zz_yy, g_yy_zz_zz_yz, g_yy_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_zz_zz_xx[i] = -2.0 * g_0_zz_zz_xx[i] * a_exp + 4.0 * g_yy_zz_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_zz_xy[i] = -2.0 * g_0_zz_zz_xy[i] * a_exp + 4.0 * g_yy_zz_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_zz_xz[i] = -2.0 * g_0_zz_zz_xz[i] * a_exp + 4.0 * g_yy_zz_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_zz_yy[i] = -2.0 * g_0_zz_zz_yy[i] * a_exp + 4.0 * g_yy_zz_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_zz_yz[i] = -2.0 * g_0_zz_zz_yz[i] * a_exp + 4.0 * g_yy_zz_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_zz_zz[i] = -2.0 * g_0_zz_zz_zz[i] * a_exp + 4.0 * g_yy_zz_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_yz_0_0_0_0_xx_xx_xx, g_yz_0_0_0_0_xx_xx_xy, g_yz_0_0_0_0_xx_xx_xz, g_yz_0_0_0_0_xx_xx_yy, g_yz_0_0_0_0_xx_xx_yz, g_yz_0_0_0_0_xx_xx_zz, g_yz_xx_xx_xx, g_yz_xx_xx_xy, g_yz_xx_xx_xz, g_yz_xx_xx_yy, g_yz_xx_xx_yz, g_yz_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xx_xx_xx[i] = 4.0 * g_yz_xx_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_xx_xy[i] = 4.0 * g_yz_xx_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_xx_xz[i] = 4.0 * g_yz_xx_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_xx_yy[i] = 4.0 * g_yz_xx_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_xx_yz[i] = 4.0 * g_yz_xx_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_xx_zz[i] = 4.0 * g_yz_xx_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_yz_0_0_0_0_xx_xy_xx, g_yz_0_0_0_0_xx_xy_xy, g_yz_0_0_0_0_xx_xy_xz, g_yz_0_0_0_0_xx_xy_yy, g_yz_0_0_0_0_xx_xy_yz, g_yz_0_0_0_0_xx_xy_zz, g_yz_xx_xy_xx, g_yz_xx_xy_xy, g_yz_xx_xy_xz, g_yz_xx_xy_yy, g_yz_xx_xy_yz, g_yz_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xx_xy_xx[i] = 4.0 * g_yz_xx_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_xy_xy[i] = 4.0 * g_yz_xx_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_xy_xz[i] = 4.0 * g_yz_xx_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_xy_yy[i] = 4.0 * g_yz_xx_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_xy_yz[i] = 4.0 * g_yz_xx_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_xy_zz[i] = 4.0 * g_yz_xx_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_yz_0_0_0_0_xx_xz_xx, g_yz_0_0_0_0_xx_xz_xy, g_yz_0_0_0_0_xx_xz_xz, g_yz_0_0_0_0_xx_xz_yy, g_yz_0_0_0_0_xx_xz_yz, g_yz_0_0_0_0_xx_xz_zz, g_yz_xx_xz_xx, g_yz_xx_xz_xy, g_yz_xx_xz_xz, g_yz_xx_xz_yy, g_yz_xx_xz_yz, g_yz_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xx_xz_xx[i] = 4.0 * g_yz_xx_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_xz_xy[i] = 4.0 * g_yz_xx_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_xz_xz[i] = 4.0 * g_yz_xx_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_xz_yy[i] = 4.0 * g_yz_xx_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_xz_yz[i] = 4.0 * g_yz_xx_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_xz_zz[i] = 4.0 * g_yz_xx_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_yz_0_0_0_0_xx_yy_xx, g_yz_0_0_0_0_xx_yy_xy, g_yz_0_0_0_0_xx_yy_xz, g_yz_0_0_0_0_xx_yy_yy, g_yz_0_0_0_0_xx_yy_yz, g_yz_0_0_0_0_xx_yy_zz, g_yz_xx_yy_xx, g_yz_xx_yy_xy, g_yz_xx_yy_xz, g_yz_xx_yy_yy, g_yz_xx_yy_yz, g_yz_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xx_yy_xx[i] = 4.0 * g_yz_xx_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_yy_xy[i] = 4.0 * g_yz_xx_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_yy_xz[i] = 4.0 * g_yz_xx_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_yy_yy[i] = 4.0 * g_yz_xx_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_yy_yz[i] = 4.0 * g_yz_xx_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_yy_zz[i] = 4.0 * g_yz_xx_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_yz_0_0_0_0_xx_yz_xx, g_yz_0_0_0_0_xx_yz_xy, g_yz_0_0_0_0_xx_yz_xz, g_yz_0_0_0_0_xx_yz_yy, g_yz_0_0_0_0_xx_yz_yz, g_yz_0_0_0_0_xx_yz_zz, g_yz_xx_yz_xx, g_yz_xx_yz_xy, g_yz_xx_yz_xz, g_yz_xx_yz_yy, g_yz_xx_yz_yz, g_yz_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xx_yz_xx[i] = 4.0 * g_yz_xx_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_yz_xy[i] = 4.0 * g_yz_xx_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_yz_xz[i] = 4.0 * g_yz_xx_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_yz_yy[i] = 4.0 * g_yz_xx_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_yz_yz[i] = 4.0 * g_yz_xx_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_yz_zz[i] = 4.0 * g_yz_xx_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_yz_0_0_0_0_xx_zz_xx, g_yz_0_0_0_0_xx_zz_xy, g_yz_0_0_0_0_xx_zz_xz, g_yz_0_0_0_0_xx_zz_yy, g_yz_0_0_0_0_xx_zz_yz, g_yz_0_0_0_0_xx_zz_zz, g_yz_xx_zz_xx, g_yz_xx_zz_xy, g_yz_xx_zz_xz, g_yz_xx_zz_yy, g_yz_xx_zz_yz, g_yz_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xx_zz_xx[i] = 4.0 * g_yz_xx_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_zz_xy[i] = 4.0 * g_yz_xx_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_zz_xz[i] = 4.0 * g_yz_xx_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_zz_yy[i] = 4.0 * g_yz_xx_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_zz_yz[i] = 4.0 * g_yz_xx_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_zz_zz[i] = 4.0 * g_yz_xx_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_yz_0_0_0_0_xy_xx_xx, g_yz_0_0_0_0_xy_xx_xy, g_yz_0_0_0_0_xy_xx_xz, g_yz_0_0_0_0_xy_xx_yy, g_yz_0_0_0_0_xy_xx_yz, g_yz_0_0_0_0_xy_xx_zz, g_yz_xy_xx_xx, g_yz_xy_xx_xy, g_yz_xy_xx_xz, g_yz_xy_xx_yy, g_yz_xy_xx_yz, g_yz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xy_xx_xx[i] = 4.0 * g_yz_xy_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_xx_xy[i] = 4.0 * g_yz_xy_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_xx_xz[i] = 4.0 * g_yz_xy_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_xx_yy[i] = 4.0 * g_yz_xy_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_xx_yz[i] = 4.0 * g_yz_xy_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_xx_zz[i] = 4.0 * g_yz_xy_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_yz_0_0_0_0_xy_xy_xx, g_yz_0_0_0_0_xy_xy_xy, g_yz_0_0_0_0_xy_xy_xz, g_yz_0_0_0_0_xy_xy_yy, g_yz_0_0_0_0_xy_xy_yz, g_yz_0_0_0_0_xy_xy_zz, g_yz_xy_xy_xx, g_yz_xy_xy_xy, g_yz_xy_xy_xz, g_yz_xy_xy_yy, g_yz_xy_xy_yz, g_yz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xy_xy_xx[i] = 4.0 * g_yz_xy_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_xy_xy[i] = 4.0 * g_yz_xy_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_xy_xz[i] = 4.0 * g_yz_xy_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_xy_yy[i] = 4.0 * g_yz_xy_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_xy_yz[i] = 4.0 * g_yz_xy_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_xy_zz[i] = 4.0 * g_yz_xy_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_yz_0_0_0_0_xy_xz_xx, g_yz_0_0_0_0_xy_xz_xy, g_yz_0_0_0_0_xy_xz_xz, g_yz_0_0_0_0_xy_xz_yy, g_yz_0_0_0_0_xy_xz_yz, g_yz_0_0_0_0_xy_xz_zz, g_yz_xy_xz_xx, g_yz_xy_xz_xy, g_yz_xy_xz_xz, g_yz_xy_xz_yy, g_yz_xy_xz_yz, g_yz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xy_xz_xx[i] = 4.0 * g_yz_xy_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_xz_xy[i] = 4.0 * g_yz_xy_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_xz_xz[i] = 4.0 * g_yz_xy_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_xz_yy[i] = 4.0 * g_yz_xy_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_xz_yz[i] = 4.0 * g_yz_xy_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_xz_zz[i] = 4.0 * g_yz_xy_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_yz_0_0_0_0_xy_yy_xx, g_yz_0_0_0_0_xy_yy_xy, g_yz_0_0_0_0_xy_yy_xz, g_yz_0_0_0_0_xy_yy_yy, g_yz_0_0_0_0_xy_yy_yz, g_yz_0_0_0_0_xy_yy_zz, g_yz_xy_yy_xx, g_yz_xy_yy_xy, g_yz_xy_yy_xz, g_yz_xy_yy_yy, g_yz_xy_yy_yz, g_yz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xy_yy_xx[i] = 4.0 * g_yz_xy_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_yy_xy[i] = 4.0 * g_yz_xy_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_yy_xz[i] = 4.0 * g_yz_xy_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_yy_yy[i] = 4.0 * g_yz_xy_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_yy_yz[i] = 4.0 * g_yz_xy_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_yy_zz[i] = 4.0 * g_yz_xy_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_yz_0_0_0_0_xy_yz_xx, g_yz_0_0_0_0_xy_yz_xy, g_yz_0_0_0_0_xy_yz_xz, g_yz_0_0_0_0_xy_yz_yy, g_yz_0_0_0_0_xy_yz_yz, g_yz_0_0_0_0_xy_yz_zz, g_yz_xy_yz_xx, g_yz_xy_yz_xy, g_yz_xy_yz_xz, g_yz_xy_yz_yy, g_yz_xy_yz_yz, g_yz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xy_yz_xx[i] = 4.0 * g_yz_xy_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_yz_xy[i] = 4.0 * g_yz_xy_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_yz_xz[i] = 4.0 * g_yz_xy_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_yz_yy[i] = 4.0 * g_yz_xy_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_yz_yz[i] = 4.0 * g_yz_xy_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_yz_zz[i] = 4.0 * g_yz_xy_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_yz_0_0_0_0_xy_zz_xx, g_yz_0_0_0_0_xy_zz_xy, g_yz_0_0_0_0_xy_zz_xz, g_yz_0_0_0_0_xy_zz_yy, g_yz_0_0_0_0_xy_zz_yz, g_yz_0_0_0_0_xy_zz_zz, g_yz_xy_zz_xx, g_yz_xy_zz_xy, g_yz_xy_zz_xz, g_yz_xy_zz_yy, g_yz_xy_zz_yz, g_yz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xy_zz_xx[i] = 4.0 * g_yz_xy_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_zz_xy[i] = 4.0 * g_yz_xy_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_zz_xz[i] = 4.0 * g_yz_xy_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_zz_yy[i] = 4.0 * g_yz_xy_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_zz_yz[i] = 4.0 * g_yz_xy_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_zz_zz[i] = 4.0 * g_yz_xy_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_yz_0_0_0_0_xz_xx_xx, g_yz_0_0_0_0_xz_xx_xy, g_yz_0_0_0_0_xz_xx_xz, g_yz_0_0_0_0_xz_xx_yy, g_yz_0_0_0_0_xz_xx_yz, g_yz_0_0_0_0_xz_xx_zz, g_yz_xz_xx_xx, g_yz_xz_xx_xy, g_yz_xz_xx_xz, g_yz_xz_xx_yy, g_yz_xz_xx_yz, g_yz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xz_xx_xx[i] = 4.0 * g_yz_xz_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_xx_xy[i] = 4.0 * g_yz_xz_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_xx_xz[i] = 4.0 * g_yz_xz_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_xx_yy[i] = 4.0 * g_yz_xz_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_xx_yz[i] = 4.0 * g_yz_xz_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_xx_zz[i] = 4.0 * g_yz_xz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_yz_0_0_0_0_xz_xy_xx, g_yz_0_0_0_0_xz_xy_xy, g_yz_0_0_0_0_xz_xy_xz, g_yz_0_0_0_0_xz_xy_yy, g_yz_0_0_0_0_xz_xy_yz, g_yz_0_0_0_0_xz_xy_zz, g_yz_xz_xy_xx, g_yz_xz_xy_xy, g_yz_xz_xy_xz, g_yz_xz_xy_yy, g_yz_xz_xy_yz, g_yz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xz_xy_xx[i] = 4.0 * g_yz_xz_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_xy_xy[i] = 4.0 * g_yz_xz_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_xy_xz[i] = 4.0 * g_yz_xz_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_xy_yy[i] = 4.0 * g_yz_xz_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_xy_yz[i] = 4.0 * g_yz_xz_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_xy_zz[i] = 4.0 * g_yz_xz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_yz_0_0_0_0_xz_xz_xx, g_yz_0_0_0_0_xz_xz_xy, g_yz_0_0_0_0_xz_xz_xz, g_yz_0_0_0_0_xz_xz_yy, g_yz_0_0_0_0_xz_xz_yz, g_yz_0_0_0_0_xz_xz_zz, g_yz_xz_xz_xx, g_yz_xz_xz_xy, g_yz_xz_xz_xz, g_yz_xz_xz_yy, g_yz_xz_xz_yz, g_yz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xz_xz_xx[i] = 4.0 * g_yz_xz_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_xz_xy[i] = 4.0 * g_yz_xz_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_xz_xz[i] = 4.0 * g_yz_xz_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_xz_yy[i] = 4.0 * g_yz_xz_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_xz_yz[i] = 4.0 * g_yz_xz_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_xz_zz[i] = 4.0 * g_yz_xz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_yz_0_0_0_0_xz_yy_xx, g_yz_0_0_0_0_xz_yy_xy, g_yz_0_0_0_0_xz_yy_xz, g_yz_0_0_0_0_xz_yy_yy, g_yz_0_0_0_0_xz_yy_yz, g_yz_0_0_0_0_xz_yy_zz, g_yz_xz_yy_xx, g_yz_xz_yy_xy, g_yz_xz_yy_xz, g_yz_xz_yy_yy, g_yz_xz_yy_yz, g_yz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xz_yy_xx[i] = 4.0 * g_yz_xz_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_yy_xy[i] = 4.0 * g_yz_xz_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_yy_xz[i] = 4.0 * g_yz_xz_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_yy_yy[i] = 4.0 * g_yz_xz_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_yy_yz[i] = 4.0 * g_yz_xz_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_yy_zz[i] = 4.0 * g_yz_xz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_yz_0_0_0_0_xz_yz_xx, g_yz_0_0_0_0_xz_yz_xy, g_yz_0_0_0_0_xz_yz_xz, g_yz_0_0_0_0_xz_yz_yy, g_yz_0_0_0_0_xz_yz_yz, g_yz_0_0_0_0_xz_yz_zz, g_yz_xz_yz_xx, g_yz_xz_yz_xy, g_yz_xz_yz_xz, g_yz_xz_yz_yy, g_yz_xz_yz_yz, g_yz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xz_yz_xx[i] = 4.0 * g_yz_xz_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_yz_xy[i] = 4.0 * g_yz_xz_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_yz_xz[i] = 4.0 * g_yz_xz_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_yz_yy[i] = 4.0 * g_yz_xz_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_yz_yz[i] = 4.0 * g_yz_xz_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_yz_zz[i] = 4.0 * g_yz_xz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_yz_0_0_0_0_xz_zz_xx, g_yz_0_0_0_0_xz_zz_xy, g_yz_0_0_0_0_xz_zz_xz, g_yz_0_0_0_0_xz_zz_yy, g_yz_0_0_0_0_xz_zz_yz, g_yz_0_0_0_0_xz_zz_zz, g_yz_xz_zz_xx, g_yz_xz_zz_xy, g_yz_xz_zz_xz, g_yz_xz_zz_yy, g_yz_xz_zz_yz, g_yz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xz_zz_xx[i] = 4.0 * g_yz_xz_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_zz_xy[i] = 4.0 * g_yz_xz_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_zz_xz[i] = 4.0 * g_yz_xz_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_zz_yy[i] = 4.0 * g_yz_xz_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_zz_yz[i] = 4.0 * g_yz_xz_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_zz_zz[i] = 4.0 * g_yz_xz_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (972-978)

    #pragma omp simd aligned(g_yz_0_0_0_0_yy_xx_xx, g_yz_0_0_0_0_yy_xx_xy, g_yz_0_0_0_0_yy_xx_xz, g_yz_0_0_0_0_yy_xx_yy, g_yz_0_0_0_0_yy_xx_yz, g_yz_0_0_0_0_yy_xx_zz, g_yz_yy_xx_xx, g_yz_yy_xx_xy, g_yz_yy_xx_xz, g_yz_yy_xx_yy, g_yz_yy_xx_yz, g_yz_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yy_xx_xx[i] = 4.0 * g_yz_yy_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_xx_xy[i] = 4.0 * g_yz_yy_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_xx_xz[i] = 4.0 * g_yz_yy_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_xx_yy[i] = 4.0 * g_yz_yy_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_xx_yz[i] = 4.0 * g_yz_yy_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_xx_zz[i] = 4.0 * g_yz_yy_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (978-984)

    #pragma omp simd aligned(g_yz_0_0_0_0_yy_xy_xx, g_yz_0_0_0_0_yy_xy_xy, g_yz_0_0_0_0_yy_xy_xz, g_yz_0_0_0_0_yy_xy_yy, g_yz_0_0_0_0_yy_xy_yz, g_yz_0_0_0_0_yy_xy_zz, g_yz_yy_xy_xx, g_yz_yy_xy_xy, g_yz_yy_xy_xz, g_yz_yy_xy_yy, g_yz_yy_xy_yz, g_yz_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yy_xy_xx[i] = 4.0 * g_yz_yy_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_xy_xy[i] = 4.0 * g_yz_yy_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_xy_xz[i] = 4.0 * g_yz_yy_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_xy_yy[i] = 4.0 * g_yz_yy_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_xy_yz[i] = 4.0 * g_yz_yy_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_xy_zz[i] = 4.0 * g_yz_yy_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (984-990)

    #pragma omp simd aligned(g_yz_0_0_0_0_yy_xz_xx, g_yz_0_0_0_0_yy_xz_xy, g_yz_0_0_0_0_yy_xz_xz, g_yz_0_0_0_0_yy_xz_yy, g_yz_0_0_0_0_yy_xz_yz, g_yz_0_0_0_0_yy_xz_zz, g_yz_yy_xz_xx, g_yz_yy_xz_xy, g_yz_yy_xz_xz, g_yz_yy_xz_yy, g_yz_yy_xz_yz, g_yz_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yy_xz_xx[i] = 4.0 * g_yz_yy_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_xz_xy[i] = 4.0 * g_yz_yy_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_xz_xz[i] = 4.0 * g_yz_yy_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_xz_yy[i] = 4.0 * g_yz_yy_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_xz_yz[i] = 4.0 * g_yz_yy_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_xz_zz[i] = 4.0 * g_yz_yy_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (990-996)

    #pragma omp simd aligned(g_yz_0_0_0_0_yy_yy_xx, g_yz_0_0_0_0_yy_yy_xy, g_yz_0_0_0_0_yy_yy_xz, g_yz_0_0_0_0_yy_yy_yy, g_yz_0_0_0_0_yy_yy_yz, g_yz_0_0_0_0_yy_yy_zz, g_yz_yy_yy_xx, g_yz_yy_yy_xy, g_yz_yy_yy_xz, g_yz_yy_yy_yy, g_yz_yy_yy_yz, g_yz_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yy_yy_xx[i] = 4.0 * g_yz_yy_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_yy_xy[i] = 4.0 * g_yz_yy_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_yy_xz[i] = 4.0 * g_yz_yy_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_yy_yy[i] = 4.0 * g_yz_yy_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_yy_yz[i] = 4.0 * g_yz_yy_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_yy_zz[i] = 4.0 * g_yz_yy_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (996-1002)

    #pragma omp simd aligned(g_yz_0_0_0_0_yy_yz_xx, g_yz_0_0_0_0_yy_yz_xy, g_yz_0_0_0_0_yy_yz_xz, g_yz_0_0_0_0_yy_yz_yy, g_yz_0_0_0_0_yy_yz_yz, g_yz_0_0_0_0_yy_yz_zz, g_yz_yy_yz_xx, g_yz_yy_yz_xy, g_yz_yy_yz_xz, g_yz_yy_yz_yy, g_yz_yy_yz_yz, g_yz_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yy_yz_xx[i] = 4.0 * g_yz_yy_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_yz_xy[i] = 4.0 * g_yz_yy_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_yz_xz[i] = 4.0 * g_yz_yy_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_yz_yy[i] = 4.0 * g_yz_yy_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_yz_yz[i] = 4.0 * g_yz_yy_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_yz_zz[i] = 4.0 * g_yz_yy_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1002-1008)

    #pragma omp simd aligned(g_yz_0_0_0_0_yy_zz_xx, g_yz_0_0_0_0_yy_zz_xy, g_yz_0_0_0_0_yy_zz_xz, g_yz_0_0_0_0_yy_zz_yy, g_yz_0_0_0_0_yy_zz_yz, g_yz_0_0_0_0_yy_zz_zz, g_yz_yy_zz_xx, g_yz_yy_zz_xy, g_yz_yy_zz_xz, g_yz_yy_zz_yy, g_yz_yy_zz_yz, g_yz_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yy_zz_xx[i] = 4.0 * g_yz_yy_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_zz_xy[i] = 4.0 * g_yz_yy_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_zz_xz[i] = 4.0 * g_yz_yy_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_zz_yy[i] = 4.0 * g_yz_yy_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_zz_yz[i] = 4.0 * g_yz_yy_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_zz_zz[i] = 4.0 * g_yz_yy_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1008-1014)

    #pragma omp simd aligned(g_yz_0_0_0_0_yz_xx_xx, g_yz_0_0_0_0_yz_xx_xy, g_yz_0_0_0_0_yz_xx_xz, g_yz_0_0_0_0_yz_xx_yy, g_yz_0_0_0_0_yz_xx_yz, g_yz_0_0_0_0_yz_xx_zz, g_yz_yz_xx_xx, g_yz_yz_xx_xy, g_yz_yz_xx_xz, g_yz_yz_xx_yy, g_yz_yz_xx_yz, g_yz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yz_xx_xx[i] = 4.0 * g_yz_yz_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_xx_xy[i] = 4.0 * g_yz_yz_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_xx_xz[i] = 4.0 * g_yz_yz_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_xx_yy[i] = 4.0 * g_yz_yz_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_xx_yz[i] = 4.0 * g_yz_yz_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_xx_zz[i] = 4.0 * g_yz_yz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1014-1020)

    #pragma omp simd aligned(g_yz_0_0_0_0_yz_xy_xx, g_yz_0_0_0_0_yz_xy_xy, g_yz_0_0_0_0_yz_xy_xz, g_yz_0_0_0_0_yz_xy_yy, g_yz_0_0_0_0_yz_xy_yz, g_yz_0_0_0_0_yz_xy_zz, g_yz_yz_xy_xx, g_yz_yz_xy_xy, g_yz_yz_xy_xz, g_yz_yz_xy_yy, g_yz_yz_xy_yz, g_yz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yz_xy_xx[i] = 4.0 * g_yz_yz_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_xy_xy[i] = 4.0 * g_yz_yz_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_xy_xz[i] = 4.0 * g_yz_yz_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_xy_yy[i] = 4.0 * g_yz_yz_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_xy_yz[i] = 4.0 * g_yz_yz_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_xy_zz[i] = 4.0 * g_yz_yz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1020-1026)

    #pragma omp simd aligned(g_yz_0_0_0_0_yz_xz_xx, g_yz_0_0_0_0_yz_xz_xy, g_yz_0_0_0_0_yz_xz_xz, g_yz_0_0_0_0_yz_xz_yy, g_yz_0_0_0_0_yz_xz_yz, g_yz_0_0_0_0_yz_xz_zz, g_yz_yz_xz_xx, g_yz_yz_xz_xy, g_yz_yz_xz_xz, g_yz_yz_xz_yy, g_yz_yz_xz_yz, g_yz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yz_xz_xx[i] = 4.0 * g_yz_yz_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_xz_xy[i] = 4.0 * g_yz_yz_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_xz_xz[i] = 4.0 * g_yz_yz_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_xz_yy[i] = 4.0 * g_yz_yz_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_xz_yz[i] = 4.0 * g_yz_yz_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_xz_zz[i] = 4.0 * g_yz_yz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1026-1032)

    #pragma omp simd aligned(g_yz_0_0_0_0_yz_yy_xx, g_yz_0_0_0_0_yz_yy_xy, g_yz_0_0_0_0_yz_yy_xz, g_yz_0_0_0_0_yz_yy_yy, g_yz_0_0_0_0_yz_yy_yz, g_yz_0_0_0_0_yz_yy_zz, g_yz_yz_yy_xx, g_yz_yz_yy_xy, g_yz_yz_yy_xz, g_yz_yz_yy_yy, g_yz_yz_yy_yz, g_yz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yz_yy_xx[i] = 4.0 * g_yz_yz_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_yy_xy[i] = 4.0 * g_yz_yz_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_yy_xz[i] = 4.0 * g_yz_yz_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_yy_yy[i] = 4.0 * g_yz_yz_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_yy_yz[i] = 4.0 * g_yz_yz_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_yy_zz[i] = 4.0 * g_yz_yz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1032-1038)

    #pragma omp simd aligned(g_yz_0_0_0_0_yz_yz_xx, g_yz_0_0_0_0_yz_yz_xy, g_yz_0_0_0_0_yz_yz_xz, g_yz_0_0_0_0_yz_yz_yy, g_yz_0_0_0_0_yz_yz_yz, g_yz_0_0_0_0_yz_yz_zz, g_yz_yz_yz_xx, g_yz_yz_yz_xy, g_yz_yz_yz_xz, g_yz_yz_yz_yy, g_yz_yz_yz_yz, g_yz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yz_yz_xx[i] = 4.0 * g_yz_yz_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_yz_xy[i] = 4.0 * g_yz_yz_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_yz_xz[i] = 4.0 * g_yz_yz_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_yz_yy[i] = 4.0 * g_yz_yz_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_yz_yz[i] = 4.0 * g_yz_yz_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_yz_zz[i] = 4.0 * g_yz_yz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1038-1044)

    #pragma omp simd aligned(g_yz_0_0_0_0_yz_zz_xx, g_yz_0_0_0_0_yz_zz_xy, g_yz_0_0_0_0_yz_zz_xz, g_yz_0_0_0_0_yz_zz_yy, g_yz_0_0_0_0_yz_zz_yz, g_yz_0_0_0_0_yz_zz_zz, g_yz_yz_zz_xx, g_yz_yz_zz_xy, g_yz_yz_zz_xz, g_yz_yz_zz_yy, g_yz_yz_zz_yz, g_yz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yz_zz_xx[i] = 4.0 * g_yz_yz_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_zz_xy[i] = 4.0 * g_yz_yz_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_zz_xz[i] = 4.0 * g_yz_yz_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_zz_yy[i] = 4.0 * g_yz_yz_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_zz_yz[i] = 4.0 * g_yz_yz_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_zz_zz[i] = 4.0 * g_yz_yz_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1044-1050)

    #pragma omp simd aligned(g_yz_0_0_0_0_zz_xx_xx, g_yz_0_0_0_0_zz_xx_xy, g_yz_0_0_0_0_zz_xx_xz, g_yz_0_0_0_0_zz_xx_yy, g_yz_0_0_0_0_zz_xx_yz, g_yz_0_0_0_0_zz_xx_zz, g_yz_zz_xx_xx, g_yz_zz_xx_xy, g_yz_zz_xx_xz, g_yz_zz_xx_yy, g_yz_zz_xx_yz, g_yz_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_zz_xx_xx[i] = 4.0 * g_yz_zz_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_xx_xy[i] = 4.0 * g_yz_zz_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_xx_xz[i] = 4.0 * g_yz_zz_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_xx_yy[i] = 4.0 * g_yz_zz_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_xx_yz[i] = 4.0 * g_yz_zz_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_xx_zz[i] = 4.0 * g_yz_zz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1050-1056)

    #pragma omp simd aligned(g_yz_0_0_0_0_zz_xy_xx, g_yz_0_0_0_0_zz_xy_xy, g_yz_0_0_0_0_zz_xy_xz, g_yz_0_0_0_0_zz_xy_yy, g_yz_0_0_0_0_zz_xy_yz, g_yz_0_0_0_0_zz_xy_zz, g_yz_zz_xy_xx, g_yz_zz_xy_xy, g_yz_zz_xy_xz, g_yz_zz_xy_yy, g_yz_zz_xy_yz, g_yz_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_zz_xy_xx[i] = 4.0 * g_yz_zz_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_xy_xy[i] = 4.0 * g_yz_zz_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_xy_xz[i] = 4.0 * g_yz_zz_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_xy_yy[i] = 4.0 * g_yz_zz_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_xy_yz[i] = 4.0 * g_yz_zz_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_xy_zz[i] = 4.0 * g_yz_zz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1056-1062)

    #pragma omp simd aligned(g_yz_0_0_0_0_zz_xz_xx, g_yz_0_0_0_0_zz_xz_xy, g_yz_0_0_0_0_zz_xz_xz, g_yz_0_0_0_0_zz_xz_yy, g_yz_0_0_0_0_zz_xz_yz, g_yz_0_0_0_0_zz_xz_zz, g_yz_zz_xz_xx, g_yz_zz_xz_xy, g_yz_zz_xz_xz, g_yz_zz_xz_yy, g_yz_zz_xz_yz, g_yz_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_zz_xz_xx[i] = 4.0 * g_yz_zz_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_xz_xy[i] = 4.0 * g_yz_zz_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_xz_xz[i] = 4.0 * g_yz_zz_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_xz_yy[i] = 4.0 * g_yz_zz_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_xz_yz[i] = 4.0 * g_yz_zz_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_xz_zz[i] = 4.0 * g_yz_zz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1062-1068)

    #pragma omp simd aligned(g_yz_0_0_0_0_zz_yy_xx, g_yz_0_0_0_0_zz_yy_xy, g_yz_0_0_0_0_zz_yy_xz, g_yz_0_0_0_0_zz_yy_yy, g_yz_0_0_0_0_zz_yy_yz, g_yz_0_0_0_0_zz_yy_zz, g_yz_zz_yy_xx, g_yz_zz_yy_xy, g_yz_zz_yy_xz, g_yz_zz_yy_yy, g_yz_zz_yy_yz, g_yz_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_zz_yy_xx[i] = 4.0 * g_yz_zz_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_yy_xy[i] = 4.0 * g_yz_zz_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_yy_xz[i] = 4.0 * g_yz_zz_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_yy_yy[i] = 4.0 * g_yz_zz_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_yy_yz[i] = 4.0 * g_yz_zz_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_yy_zz[i] = 4.0 * g_yz_zz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1068-1074)

    #pragma omp simd aligned(g_yz_0_0_0_0_zz_yz_xx, g_yz_0_0_0_0_zz_yz_xy, g_yz_0_0_0_0_zz_yz_xz, g_yz_0_0_0_0_zz_yz_yy, g_yz_0_0_0_0_zz_yz_yz, g_yz_0_0_0_0_zz_yz_zz, g_yz_zz_yz_xx, g_yz_zz_yz_xy, g_yz_zz_yz_xz, g_yz_zz_yz_yy, g_yz_zz_yz_yz, g_yz_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_zz_yz_xx[i] = 4.0 * g_yz_zz_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_yz_xy[i] = 4.0 * g_yz_zz_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_yz_xz[i] = 4.0 * g_yz_zz_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_yz_yy[i] = 4.0 * g_yz_zz_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_yz_yz[i] = 4.0 * g_yz_zz_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_yz_zz[i] = 4.0 * g_yz_zz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1074-1080)

    #pragma omp simd aligned(g_yz_0_0_0_0_zz_zz_xx, g_yz_0_0_0_0_zz_zz_xy, g_yz_0_0_0_0_zz_zz_xz, g_yz_0_0_0_0_zz_zz_yy, g_yz_0_0_0_0_zz_zz_yz, g_yz_0_0_0_0_zz_zz_zz, g_yz_zz_zz_xx, g_yz_zz_zz_xy, g_yz_zz_zz_xz, g_yz_zz_zz_yy, g_yz_zz_zz_yz, g_yz_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_zz_zz_xx[i] = 4.0 * g_yz_zz_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_zz_xy[i] = 4.0 * g_yz_zz_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_zz_xz[i] = 4.0 * g_yz_zz_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_zz_yy[i] = 4.0 * g_yz_zz_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_zz_yz[i] = 4.0 * g_yz_zz_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_zz_zz[i] = 4.0 * g_yz_zz_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1080-1086)

    #pragma omp simd aligned(g_0_xx_xx_xx, g_0_xx_xx_xy, g_0_xx_xx_xz, g_0_xx_xx_yy, g_0_xx_xx_yz, g_0_xx_xx_zz, g_zz_0_0_0_0_xx_xx_xx, g_zz_0_0_0_0_xx_xx_xy, g_zz_0_0_0_0_xx_xx_xz, g_zz_0_0_0_0_xx_xx_yy, g_zz_0_0_0_0_xx_xx_yz, g_zz_0_0_0_0_xx_xx_zz, g_zz_xx_xx_xx, g_zz_xx_xx_xy, g_zz_xx_xx_xz, g_zz_xx_xx_yy, g_zz_xx_xx_yz, g_zz_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xx_xx_xx[i] = -2.0 * g_0_xx_xx_xx[i] * a_exp + 4.0 * g_zz_xx_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_xx_xy[i] = -2.0 * g_0_xx_xx_xy[i] * a_exp + 4.0 * g_zz_xx_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_xx_xz[i] = -2.0 * g_0_xx_xx_xz[i] * a_exp + 4.0 * g_zz_xx_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_xx_yy[i] = -2.0 * g_0_xx_xx_yy[i] * a_exp + 4.0 * g_zz_xx_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_xx_yz[i] = -2.0 * g_0_xx_xx_yz[i] * a_exp + 4.0 * g_zz_xx_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_xx_zz[i] = -2.0 * g_0_xx_xx_zz[i] * a_exp + 4.0 * g_zz_xx_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1086-1092)

    #pragma omp simd aligned(g_0_xx_xy_xx, g_0_xx_xy_xy, g_0_xx_xy_xz, g_0_xx_xy_yy, g_0_xx_xy_yz, g_0_xx_xy_zz, g_zz_0_0_0_0_xx_xy_xx, g_zz_0_0_0_0_xx_xy_xy, g_zz_0_0_0_0_xx_xy_xz, g_zz_0_0_0_0_xx_xy_yy, g_zz_0_0_0_0_xx_xy_yz, g_zz_0_0_0_0_xx_xy_zz, g_zz_xx_xy_xx, g_zz_xx_xy_xy, g_zz_xx_xy_xz, g_zz_xx_xy_yy, g_zz_xx_xy_yz, g_zz_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xx_xy_xx[i] = -2.0 * g_0_xx_xy_xx[i] * a_exp + 4.0 * g_zz_xx_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_xy_xy[i] = -2.0 * g_0_xx_xy_xy[i] * a_exp + 4.0 * g_zz_xx_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_xy_xz[i] = -2.0 * g_0_xx_xy_xz[i] * a_exp + 4.0 * g_zz_xx_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_xy_yy[i] = -2.0 * g_0_xx_xy_yy[i] * a_exp + 4.0 * g_zz_xx_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_xy_yz[i] = -2.0 * g_0_xx_xy_yz[i] * a_exp + 4.0 * g_zz_xx_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_xy_zz[i] = -2.0 * g_0_xx_xy_zz[i] * a_exp + 4.0 * g_zz_xx_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1092-1098)

    #pragma omp simd aligned(g_0_xx_xz_xx, g_0_xx_xz_xy, g_0_xx_xz_xz, g_0_xx_xz_yy, g_0_xx_xz_yz, g_0_xx_xz_zz, g_zz_0_0_0_0_xx_xz_xx, g_zz_0_0_0_0_xx_xz_xy, g_zz_0_0_0_0_xx_xz_xz, g_zz_0_0_0_0_xx_xz_yy, g_zz_0_0_0_0_xx_xz_yz, g_zz_0_0_0_0_xx_xz_zz, g_zz_xx_xz_xx, g_zz_xx_xz_xy, g_zz_xx_xz_xz, g_zz_xx_xz_yy, g_zz_xx_xz_yz, g_zz_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xx_xz_xx[i] = -2.0 * g_0_xx_xz_xx[i] * a_exp + 4.0 * g_zz_xx_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_xz_xy[i] = -2.0 * g_0_xx_xz_xy[i] * a_exp + 4.0 * g_zz_xx_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_xz_xz[i] = -2.0 * g_0_xx_xz_xz[i] * a_exp + 4.0 * g_zz_xx_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_xz_yy[i] = -2.0 * g_0_xx_xz_yy[i] * a_exp + 4.0 * g_zz_xx_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_xz_yz[i] = -2.0 * g_0_xx_xz_yz[i] * a_exp + 4.0 * g_zz_xx_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_xz_zz[i] = -2.0 * g_0_xx_xz_zz[i] * a_exp + 4.0 * g_zz_xx_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1098-1104)

    #pragma omp simd aligned(g_0_xx_yy_xx, g_0_xx_yy_xy, g_0_xx_yy_xz, g_0_xx_yy_yy, g_0_xx_yy_yz, g_0_xx_yy_zz, g_zz_0_0_0_0_xx_yy_xx, g_zz_0_0_0_0_xx_yy_xy, g_zz_0_0_0_0_xx_yy_xz, g_zz_0_0_0_0_xx_yy_yy, g_zz_0_0_0_0_xx_yy_yz, g_zz_0_0_0_0_xx_yy_zz, g_zz_xx_yy_xx, g_zz_xx_yy_xy, g_zz_xx_yy_xz, g_zz_xx_yy_yy, g_zz_xx_yy_yz, g_zz_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xx_yy_xx[i] = -2.0 * g_0_xx_yy_xx[i] * a_exp + 4.0 * g_zz_xx_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_yy_xy[i] = -2.0 * g_0_xx_yy_xy[i] * a_exp + 4.0 * g_zz_xx_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_yy_xz[i] = -2.0 * g_0_xx_yy_xz[i] * a_exp + 4.0 * g_zz_xx_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_yy_yy[i] = -2.0 * g_0_xx_yy_yy[i] * a_exp + 4.0 * g_zz_xx_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_yy_yz[i] = -2.0 * g_0_xx_yy_yz[i] * a_exp + 4.0 * g_zz_xx_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_yy_zz[i] = -2.0 * g_0_xx_yy_zz[i] * a_exp + 4.0 * g_zz_xx_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1104-1110)

    #pragma omp simd aligned(g_0_xx_yz_xx, g_0_xx_yz_xy, g_0_xx_yz_xz, g_0_xx_yz_yy, g_0_xx_yz_yz, g_0_xx_yz_zz, g_zz_0_0_0_0_xx_yz_xx, g_zz_0_0_0_0_xx_yz_xy, g_zz_0_0_0_0_xx_yz_xz, g_zz_0_0_0_0_xx_yz_yy, g_zz_0_0_0_0_xx_yz_yz, g_zz_0_0_0_0_xx_yz_zz, g_zz_xx_yz_xx, g_zz_xx_yz_xy, g_zz_xx_yz_xz, g_zz_xx_yz_yy, g_zz_xx_yz_yz, g_zz_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xx_yz_xx[i] = -2.0 * g_0_xx_yz_xx[i] * a_exp + 4.0 * g_zz_xx_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_yz_xy[i] = -2.0 * g_0_xx_yz_xy[i] * a_exp + 4.0 * g_zz_xx_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_yz_xz[i] = -2.0 * g_0_xx_yz_xz[i] * a_exp + 4.0 * g_zz_xx_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_yz_yy[i] = -2.0 * g_0_xx_yz_yy[i] * a_exp + 4.0 * g_zz_xx_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_yz_yz[i] = -2.0 * g_0_xx_yz_yz[i] * a_exp + 4.0 * g_zz_xx_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_yz_zz[i] = -2.0 * g_0_xx_yz_zz[i] * a_exp + 4.0 * g_zz_xx_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1110-1116)

    #pragma omp simd aligned(g_0_xx_zz_xx, g_0_xx_zz_xy, g_0_xx_zz_xz, g_0_xx_zz_yy, g_0_xx_zz_yz, g_0_xx_zz_zz, g_zz_0_0_0_0_xx_zz_xx, g_zz_0_0_0_0_xx_zz_xy, g_zz_0_0_0_0_xx_zz_xz, g_zz_0_0_0_0_xx_zz_yy, g_zz_0_0_0_0_xx_zz_yz, g_zz_0_0_0_0_xx_zz_zz, g_zz_xx_zz_xx, g_zz_xx_zz_xy, g_zz_xx_zz_xz, g_zz_xx_zz_yy, g_zz_xx_zz_yz, g_zz_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xx_zz_xx[i] = -2.0 * g_0_xx_zz_xx[i] * a_exp + 4.0 * g_zz_xx_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_zz_xy[i] = -2.0 * g_0_xx_zz_xy[i] * a_exp + 4.0 * g_zz_xx_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_zz_xz[i] = -2.0 * g_0_xx_zz_xz[i] * a_exp + 4.0 * g_zz_xx_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_zz_yy[i] = -2.0 * g_0_xx_zz_yy[i] * a_exp + 4.0 * g_zz_xx_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_zz_yz[i] = -2.0 * g_0_xx_zz_yz[i] * a_exp + 4.0 * g_zz_xx_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_zz_zz[i] = -2.0 * g_0_xx_zz_zz[i] * a_exp + 4.0 * g_zz_xx_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1116-1122)

    #pragma omp simd aligned(g_0_xy_xx_xx, g_0_xy_xx_xy, g_0_xy_xx_xz, g_0_xy_xx_yy, g_0_xy_xx_yz, g_0_xy_xx_zz, g_zz_0_0_0_0_xy_xx_xx, g_zz_0_0_0_0_xy_xx_xy, g_zz_0_0_0_0_xy_xx_xz, g_zz_0_0_0_0_xy_xx_yy, g_zz_0_0_0_0_xy_xx_yz, g_zz_0_0_0_0_xy_xx_zz, g_zz_xy_xx_xx, g_zz_xy_xx_xy, g_zz_xy_xx_xz, g_zz_xy_xx_yy, g_zz_xy_xx_yz, g_zz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xy_xx_xx[i] = -2.0 * g_0_xy_xx_xx[i] * a_exp + 4.0 * g_zz_xy_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_xx_xy[i] = -2.0 * g_0_xy_xx_xy[i] * a_exp + 4.0 * g_zz_xy_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_xx_xz[i] = -2.0 * g_0_xy_xx_xz[i] * a_exp + 4.0 * g_zz_xy_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_xx_yy[i] = -2.0 * g_0_xy_xx_yy[i] * a_exp + 4.0 * g_zz_xy_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_xx_yz[i] = -2.0 * g_0_xy_xx_yz[i] * a_exp + 4.0 * g_zz_xy_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_xx_zz[i] = -2.0 * g_0_xy_xx_zz[i] * a_exp + 4.0 * g_zz_xy_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1122-1128)

    #pragma omp simd aligned(g_0_xy_xy_xx, g_0_xy_xy_xy, g_0_xy_xy_xz, g_0_xy_xy_yy, g_0_xy_xy_yz, g_0_xy_xy_zz, g_zz_0_0_0_0_xy_xy_xx, g_zz_0_0_0_0_xy_xy_xy, g_zz_0_0_0_0_xy_xy_xz, g_zz_0_0_0_0_xy_xy_yy, g_zz_0_0_0_0_xy_xy_yz, g_zz_0_0_0_0_xy_xy_zz, g_zz_xy_xy_xx, g_zz_xy_xy_xy, g_zz_xy_xy_xz, g_zz_xy_xy_yy, g_zz_xy_xy_yz, g_zz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xy_xy_xx[i] = -2.0 * g_0_xy_xy_xx[i] * a_exp + 4.0 * g_zz_xy_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_xy_xy[i] = -2.0 * g_0_xy_xy_xy[i] * a_exp + 4.0 * g_zz_xy_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_xy_xz[i] = -2.0 * g_0_xy_xy_xz[i] * a_exp + 4.0 * g_zz_xy_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_xy_yy[i] = -2.0 * g_0_xy_xy_yy[i] * a_exp + 4.0 * g_zz_xy_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_xy_yz[i] = -2.0 * g_0_xy_xy_yz[i] * a_exp + 4.0 * g_zz_xy_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_xy_zz[i] = -2.0 * g_0_xy_xy_zz[i] * a_exp + 4.0 * g_zz_xy_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1128-1134)

    #pragma omp simd aligned(g_0_xy_xz_xx, g_0_xy_xz_xy, g_0_xy_xz_xz, g_0_xy_xz_yy, g_0_xy_xz_yz, g_0_xy_xz_zz, g_zz_0_0_0_0_xy_xz_xx, g_zz_0_0_0_0_xy_xz_xy, g_zz_0_0_0_0_xy_xz_xz, g_zz_0_0_0_0_xy_xz_yy, g_zz_0_0_0_0_xy_xz_yz, g_zz_0_0_0_0_xy_xz_zz, g_zz_xy_xz_xx, g_zz_xy_xz_xy, g_zz_xy_xz_xz, g_zz_xy_xz_yy, g_zz_xy_xz_yz, g_zz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xy_xz_xx[i] = -2.0 * g_0_xy_xz_xx[i] * a_exp + 4.0 * g_zz_xy_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_xz_xy[i] = -2.0 * g_0_xy_xz_xy[i] * a_exp + 4.0 * g_zz_xy_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_xz_xz[i] = -2.0 * g_0_xy_xz_xz[i] * a_exp + 4.0 * g_zz_xy_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_xz_yy[i] = -2.0 * g_0_xy_xz_yy[i] * a_exp + 4.0 * g_zz_xy_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_xz_yz[i] = -2.0 * g_0_xy_xz_yz[i] * a_exp + 4.0 * g_zz_xy_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_xz_zz[i] = -2.0 * g_0_xy_xz_zz[i] * a_exp + 4.0 * g_zz_xy_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1134-1140)

    #pragma omp simd aligned(g_0_xy_yy_xx, g_0_xy_yy_xy, g_0_xy_yy_xz, g_0_xy_yy_yy, g_0_xy_yy_yz, g_0_xy_yy_zz, g_zz_0_0_0_0_xy_yy_xx, g_zz_0_0_0_0_xy_yy_xy, g_zz_0_0_0_0_xy_yy_xz, g_zz_0_0_0_0_xy_yy_yy, g_zz_0_0_0_0_xy_yy_yz, g_zz_0_0_0_0_xy_yy_zz, g_zz_xy_yy_xx, g_zz_xy_yy_xy, g_zz_xy_yy_xz, g_zz_xy_yy_yy, g_zz_xy_yy_yz, g_zz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xy_yy_xx[i] = -2.0 * g_0_xy_yy_xx[i] * a_exp + 4.0 * g_zz_xy_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_yy_xy[i] = -2.0 * g_0_xy_yy_xy[i] * a_exp + 4.0 * g_zz_xy_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_yy_xz[i] = -2.0 * g_0_xy_yy_xz[i] * a_exp + 4.0 * g_zz_xy_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_yy_yy[i] = -2.0 * g_0_xy_yy_yy[i] * a_exp + 4.0 * g_zz_xy_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_yy_yz[i] = -2.0 * g_0_xy_yy_yz[i] * a_exp + 4.0 * g_zz_xy_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_yy_zz[i] = -2.0 * g_0_xy_yy_zz[i] * a_exp + 4.0 * g_zz_xy_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1140-1146)

    #pragma omp simd aligned(g_0_xy_yz_xx, g_0_xy_yz_xy, g_0_xy_yz_xz, g_0_xy_yz_yy, g_0_xy_yz_yz, g_0_xy_yz_zz, g_zz_0_0_0_0_xy_yz_xx, g_zz_0_0_0_0_xy_yz_xy, g_zz_0_0_0_0_xy_yz_xz, g_zz_0_0_0_0_xy_yz_yy, g_zz_0_0_0_0_xy_yz_yz, g_zz_0_0_0_0_xy_yz_zz, g_zz_xy_yz_xx, g_zz_xy_yz_xy, g_zz_xy_yz_xz, g_zz_xy_yz_yy, g_zz_xy_yz_yz, g_zz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xy_yz_xx[i] = -2.0 * g_0_xy_yz_xx[i] * a_exp + 4.0 * g_zz_xy_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_yz_xy[i] = -2.0 * g_0_xy_yz_xy[i] * a_exp + 4.0 * g_zz_xy_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_yz_xz[i] = -2.0 * g_0_xy_yz_xz[i] * a_exp + 4.0 * g_zz_xy_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_yz_yy[i] = -2.0 * g_0_xy_yz_yy[i] * a_exp + 4.0 * g_zz_xy_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_yz_yz[i] = -2.0 * g_0_xy_yz_yz[i] * a_exp + 4.0 * g_zz_xy_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_yz_zz[i] = -2.0 * g_0_xy_yz_zz[i] * a_exp + 4.0 * g_zz_xy_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1146-1152)

    #pragma omp simd aligned(g_0_xy_zz_xx, g_0_xy_zz_xy, g_0_xy_zz_xz, g_0_xy_zz_yy, g_0_xy_zz_yz, g_0_xy_zz_zz, g_zz_0_0_0_0_xy_zz_xx, g_zz_0_0_0_0_xy_zz_xy, g_zz_0_0_0_0_xy_zz_xz, g_zz_0_0_0_0_xy_zz_yy, g_zz_0_0_0_0_xy_zz_yz, g_zz_0_0_0_0_xy_zz_zz, g_zz_xy_zz_xx, g_zz_xy_zz_xy, g_zz_xy_zz_xz, g_zz_xy_zz_yy, g_zz_xy_zz_yz, g_zz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xy_zz_xx[i] = -2.0 * g_0_xy_zz_xx[i] * a_exp + 4.0 * g_zz_xy_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_zz_xy[i] = -2.0 * g_0_xy_zz_xy[i] * a_exp + 4.0 * g_zz_xy_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_zz_xz[i] = -2.0 * g_0_xy_zz_xz[i] * a_exp + 4.0 * g_zz_xy_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_zz_yy[i] = -2.0 * g_0_xy_zz_yy[i] * a_exp + 4.0 * g_zz_xy_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_zz_yz[i] = -2.0 * g_0_xy_zz_yz[i] * a_exp + 4.0 * g_zz_xy_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_zz_zz[i] = -2.0 * g_0_xy_zz_zz[i] * a_exp + 4.0 * g_zz_xy_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1152-1158)

    #pragma omp simd aligned(g_0_xz_xx_xx, g_0_xz_xx_xy, g_0_xz_xx_xz, g_0_xz_xx_yy, g_0_xz_xx_yz, g_0_xz_xx_zz, g_zz_0_0_0_0_xz_xx_xx, g_zz_0_0_0_0_xz_xx_xy, g_zz_0_0_0_0_xz_xx_xz, g_zz_0_0_0_0_xz_xx_yy, g_zz_0_0_0_0_xz_xx_yz, g_zz_0_0_0_0_xz_xx_zz, g_zz_xz_xx_xx, g_zz_xz_xx_xy, g_zz_xz_xx_xz, g_zz_xz_xx_yy, g_zz_xz_xx_yz, g_zz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xz_xx_xx[i] = -2.0 * g_0_xz_xx_xx[i] * a_exp + 4.0 * g_zz_xz_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_xx_xy[i] = -2.0 * g_0_xz_xx_xy[i] * a_exp + 4.0 * g_zz_xz_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_xx_xz[i] = -2.0 * g_0_xz_xx_xz[i] * a_exp + 4.0 * g_zz_xz_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_xx_yy[i] = -2.0 * g_0_xz_xx_yy[i] * a_exp + 4.0 * g_zz_xz_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_xx_yz[i] = -2.0 * g_0_xz_xx_yz[i] * a_exp + 4.0 * g_zz_xz_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_xx_zz[i] = -2.0 * g_0_xz_xx_zz[i] * a_exp + 4.0 * g_zz_xz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1158-1164)

    #pragma omp simd aligned(g_0_xz_xy_xx, g_0_xz_xy_xy, g_0_xz_xy_xz, g_0_xz_xy_yy, g_0_xz_xy_yz, g_0_xz_xy_zz, g_zz_0_0_0_0_xz_xy_xx, g_zz_0_0_0_0_xz_xy_xy, g_zz_0_0_0_0_xz_xy_xz, g_zz_0_0_0_0_xz_xy_yy, g_zz_0_0_0_0_xz_xy_yz, g_zz_0_0_0_0_xz_xy_zz, g_zz_xz_xy_xx, g_zz_xz_xy_xy, g_zz_xz_xy_xz, g_zz_xz_xy_yy, g_zz_xz_xy_yz, g_zz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xz_xy_xx[i] = -2.0 * g_0_xz_xy_xx[i] * a_exp + 4.0 * g_zz_xz_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_xy_xy[i] = -2.0 * g_0_xz_xy_xy[i] * a_exp + 4.0 * g_zz_xz_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_xy_xz[i] = -2.0 * g_0_xz_xy_xz[i] * a_exp + 4.0 * g_zz_xz_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_xy_yy[i] = -2.0 * g_0_xz_xy_yy[i] * a_exp + 4.0 * g_zz_xz_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_xy_yz[i] = -2.0 * g_0_xz_xy_yz[i] * a_exp + 4.0 * g_zz_xz_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_xy_zz[i] = -2.0 * g_0_xz_xy_zz[i] * a_exp + 4.0 * g_zz_xz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1164-1170)

    #pragma omp simd aligned(g_0_xz_xz_xx, g_0_xz_xz_xy, g_0_xz_xz_xz, g_0_xz_xz_yy, g_0_xz_xz_yz, g_0_xz_xz_zz, g_zz_0_0_0_0_xz_xz_xx, g_zz_0_0_0_0_xz_xz_xy, g_zz_0_0_0_0_xz_xz_xz, g_zz_0_0_0_0_xz_xz_yy, g_zz_0_0_0_0_xz_xz_yz, g_zz_0_0_0_0_xz_xz_zz, g_zz_xz_xz_xx, g_zz_xz_xz_xy, g_zz_xz_xz_xz, g_zz_xz_xz_yy, g_zz_xz_xz_yz, g_zz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xz_xz_xx[i] = -2.0 * g_0_xz_xz_xx[i] * a_exp + 4.0 * g_zz_xz_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_xz_xy[i] = -2.0 * g_0_xz_xz_xy[i] * a_exp + 4.0 * g_zz_xz_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_xz_xz[i] = -2.0 * g_0_xz_xz_xz[i] * a_exp + 4.0 * g_zz_xz_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_xz_yy[i] = -2.0 * g_0_xz_xz_yy[i] * a_exp + 4.0 * g_zz_xz_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_xz_yz[i] = -2.0 * g_0_xz_xz_yz[i] * a_exp + 4.0 * g_zz_xz_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_xz_zz[i] = -2.0 * g_0_xz_xz_zz[i] * a_exp + 4.0 * g_zz_xz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1170-1176)

    #pragma omp simd aligned(g_0_xz_yy_xx, g_0_xz_yy_xy, g_0_xz_yy_xz, g_0_xz_yy_yy, g_0_xz_yy_yz, g_0_xz_yy_zz, g_zz_0_0_0_0_xz_yy_xx, g_zz_0_0_0_0_xz_yy_xy, g_zz_0_0_0_0_xz_yy_xz, g_zz_0_0_0_0_xz_yy_yy, g_zz_0_0_0_0_xz_yy_yz, g_zz_0_0_0_0_xz_yy_zz, g_zz_xz_yy_xx, g_zz_xz_yy_xy, g_zz_xz_yy_xz, g_zz_xz_yy_yy, g_zz_xz_yy_yz, g_zz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xz_yy_xx[i] = -2.0 * g_0_xz_yy_xx[i] * a_exp + 4.0 * g_zz_xz_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_yy_xy[i] = -2.0 * g_0_xz_yy_xy[i] * a_exp + 4.0 * g_zz_xz_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_yy_xz[i] = -2.0 * g_0_xz_yy_xz[i] * a_exp + 4.0 * g_zz_xz_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_yy_yy[i] = -2.0 * g_0_xz_yy_yy[i] * a_exp + 4.0 * g_zz_xz_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_yy_yz[i] = -2.0 * g_0_xz_yy_yz[i] * a_exp + 4.0 * g_zz_xz_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_yy_zz[i] = -2.0 * g_0_xz_yy_zz[i] * a_exp + 4.0 * g_zz_xz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1176-1182)

    #pragma omp simd aligned(g_0_xz_yz_xx, g_0_xz_yz_xy, g_0_xz_yz_xz, g_0_xz_yz_yy, g_0_xz_yz_yz, g_0_xz_yz_zz, g_zz_0_0_0_0_xz_yz_xx, g_zz_0_0_0_0_xz_yz_xy, g_zz_0_0_0_0_xz_yz_xz, g_zz_0_0_0_0_xz_yz_yy, g_zz_0_0_0_0_xz_yz_yz, g_zz_0_0_0_0_xz_yz_zz, g_zz_xz_yz_xx, g_zz_xz_yz_xy, g_zz_xz_yz_xz, g_zz_xz_yz_yy, g_zz_xz_yz_yz, g_zz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xz_yz_xx[i] = -2.0 * g_0_xz_yz_xx[i] * a_exp + 4.0 * g_zz_xz_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_yz_xy[i] = -2.0 * g_0_xz_yz_xy[i] * a_exp + 4.0 * g_zz_xz_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_yz_xz[i] = -2.0 * g_0_xz_yz_xz[i] * a_exp + 4.0 * g_zz_xz_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_yz_yy[i] = -2.0 * g_0_xz_yz_yy[i] * a_exp + 4.0 * g_zz_xz_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_yz_yz[i] = -2.0 * g_0_xz_yz_yz[i] * a_exp + 4.0 * g_zz_xz_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_yz_zz[i] = -2.0 * g_0_xz_yz_zz[i] * a_exp + 4.0 * g_zz_xz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1182-1188)

    #pragma omp simd aligned(g_0_xz_zz_xx, g_0_xz_zz_xy, g_0_xz_zz_xz, g_0_xz_zz_yy, g_0_xz_zz_yz, g_0_xz_zz_zz, g_zz_0_0_0_0_xz_zz_xx, g_zz_0_0_0_0_xz_zz_xy, g_zz_0_0_0_0_xz_zz_xz, g_zz_0_0_0_0_xz_zz_yy, g_zz_0_0_0_0_xz_zz_yz, g_zz_0_0_0_0_xz_zz_zz, g_zz_xz_zz_xx, g_zz_xz_zz_xy, g_zz_xz_zz_xz, g_zz_xz_zz_yy, g_zz_xz_zz_yz, g_zz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xz_zz_xx[i] = -2.0 * g_0_xz_zz_xx[i] * a_exp + 4.0 * g_zz_xz_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_zz_xy[i] = -2.0 * g_0_xz_zz_xy[i] * a_exp + 4.0 * g_zz_xz_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_zz_xz[i] = -2.0 * g_0_xz_zz_xz[i] * a_exp + 4.0 * g_zz_xz_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_zz_yy[i] = -2.0 * g_0_xz_zz_yy[i] * a_exp + 4.0 * g_zz_xz_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_zz_yz[i] = -2.0 * g_0_xz_zz_yz[i] * a_exp + 4.0 * g_zz_xz_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_zz_zz[i] = -2.0 * g_0_xz_zz_zz[i] * a_exp + 4.0 * g_zz_xz_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1188-1194)

    #pragma omp simd aligned(g_0_yy_xx_xx, g_0_yy_xx_xy, g_0_yy_xx_xz, g_0_yy_xx_yy, g_0_yy_xx_yz, g_0_yy_xx_zz, g_zz_0_0_0_0_yy_xx_xx, g_zz_0_0_0_0_yy_xx_xy, g_zz_0_0_0_0_yy_xx_xz, g_zz_0_0_0_0_yy_xx_yy, g_zz_0_0_0_0_yy_xx_yz, g_zz_0_0_0_0_yy_xx_zz, g_zz_yy_xx_xx, g_zz_yy_xx_xy, g_zz_yy_xx_xz, g_zz_yy_xx_yy, g_zz_yy_xx_yz, g_zz_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yy_xx_xx[i] = -2.0 * g_0_yy_xx_xx[i] * a_exp + 4.0 * g_zz_yy_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_xx_xy[i] = -2.0 * g_0_yy_xx_xy[i] * a_exp + 4.0 * g_zz_yy_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_xx_xz[i] = -2.0 * g_0_yy_xx_xz[i] * a_exp + 4.0 * g_zz_yy_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_xx_yy[i] = -2.0 * g_0_yy_xx_yy[i] * a_exp + 4.0 * g_zz_yy_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_xx_yz[i] = -2.0 * g_0_yy_xx_yz[i] * a_exp + 4.0 * g_zz_yy_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_xx_zz[i] = -2.0 * g_0_yy_xx_zz[i] * a_exp + 4.0 * g_zz_yy_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1194-1200)

    #pragma omp simd aligned(g_0_yy_xy_xx, g_0_yy_xy_xy, g_0_yy_xy_xz, g_0_yy_xy_yy, g_0_yy_xy_yz, g_0_yy_xy_zz, g_zz_0_0_0_0_yy_xy_xx, g_zz_0_0_0_0_yy_xy_xy, g_zz_0_0_0_0_yy_xy_xz, g_zz_0_0_0_0_yy_xy_yy, g_zz_0_0_0_0_yy_xy_yz, g_zz_0_0_0_0_yy_xy_zz, g_zz_yy_xy_xx, g_zz_yy_xy_xy, g_zz_yy_xy_xz, g_zz_yy_xy_yy, g_zz_yy_xy_yz, g_zz_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yy_xy_xx[i] = -2.0 * g_0_yy_xy_xx[i] * a_exp + 4.0 * g_zz_yy_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_xy_xy[i] = -2.0 * g_0_yy_xy_xy[i] * a_exp + 4.0 * g_zz_yy_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_xy_xz[i] = -2.0 * g_0_yy_xy_xz[i] * a_exp + 4.0 * g_zz_yy_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_xy_yy[i] = -2.0 * g_0_yy_xy_yy[i] * a_exp + 4.0 * g_zz_yy_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_xy_yz[i] = -2.0 * g_0_yy_xy_yz[i] * a_exp + 4.0 * g_zz_yy_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_xy_zz[i] = -2.0 * g_0_yy_xy_zz[i] * a_exp + 4.0 * g_zz_yy_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1200-1206)

    #pragma omp simd aligned(g_0_yy_xz_xx, g_0_yy_xz_xy, g_0_yy_xz_xz, g_0_yy_xz_yy, g_0_yy_xz_yz, g_0_yy_xz_zz, g_zz_0_0_0_0_yy_xz_xx, g_zz_0_0_0_0_yy_xz_xy, g_zz_0_0_0_0_yy_xz_xz, g_zz_0_0_0_0_yy_xz_yy, g_zz_0_0_0_0_yy_xz_yz, g_zz_0_0_0_0_yy_xz_zz, g_zz_yy_xz_xx, g_zz_yy_xz_xy, g_zz_yy_xz_xz, g_zz_yy_xz_yy, g_zz_yy_xz_yz, g_zz_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yy_xz_xx[i] = -2.0 * g_0_yy_xz_xx[i] * a_exp + 4.0 * g_zz_yy_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_xz_xy[i] = -2.0 * g_0_yy_xz_xy[i] * a_exp + 4.0 * g_zz_yy_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_xz_xz[i] = -2.0 * g_0_yy_xz_xz[i] * a_exp + 4.0 * g_zz_yy_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_xz_yy[i] = -2.0 * g_0_yy_xz_yy[i] * a_exp + 4.0 * g_zz_yy_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_xz_yz[i] = -2.0 * g_0_yy_xz_yz[i] * a_exp + 4.0 * g_zz_yy_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_xz_zz[i] = -2.0 * g_0_yy_xz_zz[i] * a_exp + 4.0 * g_zz_yy_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1206-1212)

    #pragma omp simd aligned(g_0_yy_yy_xx, g_0_yy_yy_xy, g_0_yy_yy_xz, g_0_yy_yy_yy, g_0_yy_yy_yz, g_0_yy_yy_zz, g_zz_0_0_0_0_yy_yy_xx, g_zz_0_0_0_0_yy_yy_xy, g_zz_0_0_0_0_yy_yy_xz, g_zz_0_0_0_0_yy_yy_yy, g_zz_0_0_0_0_yy_yy_yz, g_zz_0_0_0_0_yy_yy_zz, g_zz_yy_yy_xx, g_zz_yy_yy_xy, g_zz_yy_yy_xz, g_zz_yy_yy_yy, g_zz_yy_yy_yz, g_zz_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yy_yy_xx[i] = -2.0 * g_0_yy_yy_xx[i] * a_exp + 4.0 * g_zz_yy_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_yy_xy[i] = -2.0 * g_0_yy_yy_xy[i] * a_exp + 4.0 * g_zz_yy_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_yy_xz[i] = -2.0 * g_0_yy_yy_xz[i] * a_exp + 4.0 * g_zz_yy_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_yy_yy[i] = -2.0 * g_0_yy_yy_yy[i] * a_exp + 4.0 * g_zz_yy_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_yy_yz[i] = -2.0 * g_0_yy_yy_yz[i] * a_exp + 4.0 * g_zz_yy_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_yy_zz[i] = -2.0 * g_0_yy_yy_zz[i] * a_exp + 4.0 * g_zz_yy_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1212-1218)

    #pragma omp simd aligned(g_0_yy_yz_xx, g_0_yy_yz_xy, g_0_yy_yz_xz, g_0_yy_yz_yy, g_0_yy_yz_yz, g_0_yy_yz_zz, g_zz_0_0_0_0_yy_yz_xx, g_zz_0_0_0_0_yy_yz_xy, g_zz_0_0_0_0_yy_yz_xz, g_zz_0_0_0_0_yy_yz_yy, g_zz_0_0_0_0_yy_yz_yz, g_zz_0_0_0_0_yy_yz_zz, g_zz_yy_yz_xx, g_zz_yy_yz_xy, g_zz_yy_yz_xz, g_zz_yy_yz_yy, g_zz_yy_yz_yz, g_zz_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yy_yz_xx[i] = -2.0 * g_0_yy_yz_xx[i] * a_exp + 4.0 * g_zz_yy_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_yz_xy[i] = -2.0 * g_0_yy_yz_xy[i] * a_exp + 4.0 * g_zz_yy_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_yz_xz[i] = -2.0 * g_0_yy_yz_xz[i] * a_exp + 4.0 * g_zz_yy_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_yz_yy[i] = -2.0 * g_0_yy_yz_yy[i] * a_exp + 4.0 * g_zz_yy_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_yz_yz[i] = -2.0 * g_0_yy_yz_yz[i] * a_exp + 4.0 * g_zz_yy_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_yz_zz[i] = -2.0 * g_0_yy_yz_zz[i] * a_exp + 4.0 * g_zz_yy_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1218-1224)

    #pragma omp simd aligned(g_0_yy_zz_xx, g_0_yy_zz_xy, g_0_yy_zz_xz, g_0_yy_zz_yy, g_0_yy_zz_yz, g_0_yy_zz_zz, g_zz_0_0_0_0_yy_zz_xx, g_zz_0_0_0_0_yy_zz_xy, g_zz_0_0_0_0_yy_zz_xz, g_zz_0_0_0_0_yy_zz_yy, g_zz_0_0_0_0_yy_zz_yz, g_zz_0_0_0_0_yy_zz_zz, g_zz_yy_zz_xx, g_zz_yy_zz_xy, g_zz_yy_zz_xz, g_zz_yy_zz_yy, g_zz_yy_zz_yz, g_zz_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yy_zz_xx[i] = -2.0 * g_0_yy_zz_xx[i] * a_exp + 4.0 * g_zz_yy_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_zz_xy[i] = -2.0 * g_0_yy_zz_xy[i] * a_exp + 4.0 * g_zz_yy_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_zz_xz[i] = -2.0 * g_0_yy_zz_xz[i] * a_exp + 4.0 * g_zz_yy_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_zz_yy[i] = -2.0 * g_0_yy_zz_yy[i] * a_exp + 4.0 * g_zz_yy_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_zz_yz[i] = -2.0 * g_0_yy_zz_yz[i] * a_exp + 4.0 * g_zz_yy_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_zz_zz[i] = -2.0 * g_0_yy_zz_zz[i] * a_exp + 4.0 * g_zz_yy_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1224-1230)

    #pragma omp simd aligned(g_0_yz_xx_xx, g_0_yz_xx_xy, g_0_yz_xx_xz, g_0_yz_xx_yy, g_0_yz_xx_yz, g_0_yz_xx_zz, g_zz_0_0_0_0_yz_xx_xx, g_zz_0_0_0_0_yz_xx_xy, g_zz_0_0_0_0_yz_xx_xz, g_zz_0_0_0_0_yz_xx_yy, g_zz_0_0_0_0_yz_xx_yz, g_zz_0_0_0_0_yz_xx_zz, g_zz_yz_xx_xx, g_zz_yz_xx_xy, g_zz_yz_xx_xz, g_zz_yz_xx_yy, g_zz_yz_xx_yz, g_zz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yz_xx_xx[i] = -2.0 * g_0_yz_xx_xx[i] * a_exp + 4.0 * g_zz_yz_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_xx_xy[i] = -2.0 * g_0_yz_xx_xy[i] * a_exp + 4.0 * g_zz_yz_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_xx_xz[i] = -2.0 * g_0_yz_xx_xz[i] * a_exp + 4.0 * g_zz_yz_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_xx_yy[i] = -2.0 * g_0_yz_xx_yy[i] * a_exp + 4.0 * g_zz_yz_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_xx_yz[i] = -2.0 * g_0_yz_xx_yz[i] * a_exp + 4.0 * g_zz_yz_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_xx_zz[i] = -2.0 * g_0_yz_xx_zz[i] * a_exp + 4.0 * g_zz_yz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1230-1236)

    #pragma omp simd aligned(g_0_yz_xy_xx, g_0_yz_xy_xy, g_0_yz_xy_xz, g_0_yz_xy_yy, g_0_yz_xy_yz, g_0_yz_xy_zz, g_zz_0_0_0_0_yz_xy_xx, g_zz_0_0_0_0_yz_xy_xy, g_zz_0_0_0_0_yz_xy_xz, g_zz_0_0_0_0_yz_xy_yy, g_zz_0_0_0_0_yz_xy_yz, g_zz_0_0_0_0_yz_xy_zz, g_zz_yz_xy_xx, g_zz_yz_xy_xy, g_zz_yz_xy_xz, g_zz_yz_xy_yy, g_zz_yz_xy_yz, g_zz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yz_xy_xx[i] = -2.0 * g_0_yz_xy_xx[i] * a_exp + 4.0 * g_zz_yz_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_xy_xy[i] = -2.0 * g_0_yz_xy_xy[i] * a_exp + 4.0 * g_zz_yz_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_xy_xz[i] = -2.0 * g_0_yz_xy_xz[i] * a_exp + 4.0 * g_zz_yz_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_xy_yy[i] = -2.0 * g_0_yz_xy_yy[i] * a_exp + 4.0 * g_zz_yz_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_xy_yz[i] = -2.0 * g_0_yz_xy_yz[i] * a_exp + 4.0 * g_zz_yz_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_xy_zz[i] = -2.0 * g_0_yz_xy_zz[i] * a_exp + 4.0 * g_zz_yz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1236-1242)

    #pragma omp simd aligned(g_0_yz_xz_xx, g_0_yz_xz_xy, g_0_yz_xz_xz, g_0_yz_xz_yy, g_0_yz_xz_yz, g_0_yz_xz_zz, g_zz_0_0_0_0_yz_xz_xx, g_zz_0_0_0_0_yz_xz_xy, g_zz_0_0_0_0_yz_xz_xz, g_zz_0_0_0_0_yz_xz_yy, g_zz_0_0_0_0_yz_xz_yz, g_zz_0_0_0_0_yz_xz_zz, g_zz_yz_xz_xx, g_zz_yz_xz_xy, g_zz_yz_xz_xz, g_zz_yz_xz_yy, g_zz_yz_xz_yz, g_zz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yz_xz_xx[i] = -2.0 * g_0_yz_xz_xx[i] * a_exp + 4.0 * g_zz_yz_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_xz_xy[i] = -2.0 * g_0_yz_xz_xy[i] * a_exp + 4.0 * g_zz_yz_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_xz_xz[i] = -2.0 * g_0_yz_xz_xz[i] * a_exp + 4.0 * g_zz_yz_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_xz_yy[i] = -2.0 * g_0_yz_xz_yy[i] * a_exp + 4.0 * g_zz_yz_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_xz_yz[i] = -2.0 * g_0_yz_xz_yz[i] * a_exp + 4.0 * g_zz_yz_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_xz_zz[i] = -2.0 * g_0_yz_xz_zz[i] * a_exp + 4.0 * g_zz_yz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1242-1248)

    #pragma omp simd aligned(g_0_yz_yy_xx, g_0_yz_yy_xy, g_0_yz_yy_xz, g_0_yz_yy_yy, g_0_yz_yy_yz, g_0_yz_yy_zz, g_zz_0_0_0_0_yz_yy_xx, g_zz_0_0_0_0_yz_yy_xy, g_zz_0_0_0_0_yz_yy_xz, g_zz_0_0_0_0_yz_yy_yy, g_zz_0_0_0_0_yz_yy_yz, g_zz_0_0_0_0_yz_yy_zz, g_zz_yz_yy_xx, g_zz_yz_yy_xy, g_zz_yz_yy_xz, g_zz_yz_yy_yy, g_zz_yz_yy_yz, g_zz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yz_yy_xx[i] = -2.0 * g_0_yz_yy_xx[i] * a_exp + 4.0 * g_zz_yz_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_yy_xy[i] = -2.0 * g_0_yz_yy_xy[i] * a_exp + 4.0 * g_zz_yz_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_yy_xz[i] = -2.0 * g_0_yz_yy_xz[i] * a_exp + 4.0 * g_zz_yz_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_yy_yy[i] = -2.0 * g_0_yz_yy_yy[i] * a_exp + 4.0 * g_zz_yz_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_yy_yz[i] = -2.0 * g_0_yz_yy_yz[i] * a_exp + 4.0 * g_zz_yz_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_yy_zz[i] = -2.0 * g_0_yz_yy_zz[i] * a_exp + 4.0 * g_zz_yz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1248-1254)

    #pragma omp simd aligned(g_0_yz_yz_xx, g_0_yz_yz_xy, g_0_yz_yz_xz, g_0_yz_yz_yy, g_0_yz_yz_yz, g_0_yz_yz_zz, g_zz_0_0_0_0_yz_yz_xx, g_zz_0_0_0_0_yz_yz_xy, g_zz_0_0_0_0_yz_yz_xz, g_zz_0_0_0_0_yz_yz_yy, g_zz_0_0_0_0_yz_yz_yz, g_zz_0_0_0_0_yz_yz_zz, g_zz_yz_yz_xx, g_zz_yz_yz_xy, g_zz_yz_yz_xz, g_zz_yz_yz_yy, g_zz_yz_yz_yz, g_zz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yz_yz_xx[i] = -2.0 * g_0_yz_yz_xx[i] * a_exp + 4.0 * g_zz_yz_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_yz_xy[i] = -2.0 * g_0_yz_yz_xy[i] * a_exp + 4.0 * g_zz_yz_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_yz_xz[i] = -2.0 * g_0_yz_yz_xz[i] * a_exp + 4.0 * g_zz_yz_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_yz_yy[i] = -2.0 * g_0_yz_yz_yy[i] * a_exp + 4.0 * g_zz_yz_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_yz_yz[i] = -2.0 * g_0_yz_yz_yz[i] * a_exp + 4.0 * g_zz_yz_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_yz_zz[i] = -2.0 * g_0_yz_yz_zz[i] * a_exp + 4.0 * g_zz_yz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1254-1260)

    #pragma omp simd aligned(g_0_yz_zz_xx, g_0_yz_zz_xy, g_0_yz_zz_xz, g_0_yz_zz_yy, g_0_yz_zz_yz, g_0_yz_zz_zz, g_zz_0_0_0_0_yz_zz_xx, g_zz_0_0_0_0_yz_zz_xy, g_zz_0_0_0_0_yz_zz_xz, g_zz_0_0_0_0_yz_zz_yy, g_zz_0_0_0_0_yz_zz_yz, g_zz_0_0_0_0_yz_zz_zz, g_zz_yz_zz_xx, g_zz_yz_zz_xy, g_zz_yz_zz_xz, g_zz_yz_zz_yy, g_zz_yz_zz_yz, g_zz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yz_zz_xx[i] = -2.0 * g_0_yz_zz_xx[i] * a_exp + 4.0 * g_zz_yz_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_zz_xy[i] = -2.0 * g_0_yz_zz_xy[i] * a_exp + 4.0 * g_zz_yz_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_zz_xz[i] = -2.0 * g_0_yz_zz_xz[i] * a_exp + 4.0 * g_zz_yz_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_zz_yy[i] = -2.0 * g_0_yz_zz_yy[i] * a_exp + 4.0 * g_zz_yz_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_zz_yz[i] = -2.0 * g_0_yz_zz_yz[i] * a_exp + 4.0 * g_zz_yz_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_zz_zz[i] = -2.0 * g_0_yz_zz_zz[i] * a_exp + 4.0 * g_zz_yz_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1260-1266)

    #pragma omp simd aligned(g_0_zz_xx_xx, g_0_zz_xx_xy, g_0_zz_xx_xz, g_0_zz_xx_yy, g_0_zz_xx_yz, g_0_zz_xx_zz, g_zz_0_0_0_0_zz_xx_xx, g_zz_0_0_0_0_zz_xx_xy, g_zz_0_0_0_0_zz_xx_xz, g_zz_0_0_0_0_zz_xx_yy, g_zz_0_0_0_0_zz_xx_yz, g_zz_0_0_0_0_zz_xx_zz, g_zz_zz_xx_xx, g_zz_zz_xx_xy, g_zz_zz_xx_xz, g_zz_zz_xx_yy, g_zz_zz_xx_yz, g_zz_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_zz_xx_xx[i] = -2.0 * g_0_zz_xx_xx[i] * a_exp + 4.0 * g_zz_zz_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_xx_xy[i] = -2.0 * g_0_zz_xx_xy[i] * a_exp + 4.0 * g_zz_zz_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_xx_xz[i] = -2.0 * g_0_zz_xx_xz[i] * a_exp + 4.0 * g_zz_zz_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_xx_yy[i] = -2.0 * g_0_zz_xx_yy[i] * a_exp + 4.0 * g_zz_zz_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_xx_yz[i] = -2.0 * g_0_zz_xx_yz[i] * a_exp + 4.0 * g_zz_zz_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_xx_zz[i] = -2.0 * g_0_zz_xx_zz[i] * a_exp + 4.0 * g_zz_zz_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (1266-1272)

    #pragma omp simd aligned(g_0_zz_xy_xx, g_0_zz_xy_xy, g_0_zz_xy_xz, g_0_zz_xy_yy, g_0_zz_xy_yz, g_0_zz_xy_zz, g_zz_0_0_0_0_zz_xy_xx, g_zz_0_0_0_0_zz_xy_xy, g_zz_0_0_0_0_zz_xy_xz, g_zz_0_0_0_0_zz_xy_yy, g_zz_0_0_0_0_zz_xy_yz, g_zz_0_0_0_0_zz_xy_zz, g_zz_zz_xy_xx, g_zz_zz_xy_xy, g_zz_zz_xy_xz, g_zz_zz_xy_yy, g_zz_zz_xy_yz, g_zz_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_zz_xy_xx[i] = -2.0 * g_0_zz_xy_xx[i] * a_exp + 4.0 * g_zz_zz_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_xy_xy[i] = -2.0 * g_0_zz_xy_xy[i] * a_exp + 4.0 * g_zz_zz_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_xy_xz[i] = -2.0 * g_0_zz_xy_xz[i] * a_exp + 4.0 * g_zz_zz_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_xy_yy[i] = -2.0 * g_0_zz_xy_yy[i] * a_exp + 4.0 * g_zz_zz_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_xy_yz[i] = -2.0 * g_0_zz_xy_yz[i] * a_exp + 4.0 * g_zz_zz_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_xy_zz[i] = -2.0 * g_0_zz_xy_zz[i] * a_exp + 4.0 * g_zz_zz_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1272-1278)

    #pragma omp simd aligned(g_0_zz_xz_xx, g_0_zz_xz_xy, g_0_zz_xz_xz, g_0_zz_xz_yy, g_0_zz_xz_yz, g_0_zz_xz_zz, g_zz_0_0_0_0_zz_xz_xx, g_zz_0_0_0_0_zz_xz_xy, g_zz_0_0_0_0_zz_xz_xz, g_zz_0_0_0_0_zz_xz_yy, g_zz_0_0_0_0_zz_xz_yz, g_zz_0_0_0_0_zz_xz_zz, g_zz_zz_xz_xx, g_zz_zz_xz_xy, g_zz_zz_xz_xz, g_zz_zz_xz_yy, g_zz_zz_xz_yz, g_zz_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_zz_xz_xx[i] = -2.0 * g_0_zz_xz_xx[i] * a_exp + 4.0 * g_zz_zz_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_xz_xy[i] = -2.0 * g_0_zz_xz_xy[i] * a_exp + 4.0 * g_zz_zz_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_xz_xz[i] = -2.0 * g_0_zz_xz_xz[i] * a_exp + 4.0 * g_zz_zz_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_xz_yy[i] = -2.0 * g_0_zz_xz_yy[i] * a_exp + 4.0 * g_zz_zz_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_xz_yz[i] = -2.0 * g_0_zz_xz_yz[i] * a_exp + 4.0 * g_zz_zz_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_xz_zz[i] = -2.0 * g_0_zz_xz_zz[i] * a_exp + 4.0 * g_zz_zz_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1278-1284)

    #pragma omp simd aligned(g_0_zz_yy_xx, g_0_zz_yy_xy, g_0_zz_yy_xz, g_0_zz_yy_yy, g_0_zz_yy_yz, g_0_zz_yy_zz, g_zz_0_0_0_0_zz_yy_xx, g_zz_0_0_0_0_zz_yy_xy, g_zz_0_0_0_0_zz_yy_xz, g_zz_0_0_0_0_zz_yy_yy, g_zz_0_0_0_0_zz_yy_yz, g_zz_0_0_0_0_zz_yy_zz, g_zz_zz_yy_xx, g_zz_zz_yy_xy, g_zz_zz_yy_xz, g_zz_zz_yy_yy, g_zz_zz_yy_yz, g_zz_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_zz_yy_xx[i] = -2.0 * g_0_zz_yy_xx[i] * a_exp + 4.0 * g_zz_zz_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_yy_xy[i] = -2.0 * g_0_zz_yy_xy[i] * a_exp + 4.0 * g_zz_zz_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_yy_xz[i] = -2.0 * g_0_zz_yy_xz[i] * a_exp + 4.0 * g_zz_zz_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_yy_yy[i] = -2.0 * g_0_zz_yy_yy[i] * a_exp + 4.0 * g_zz_zz_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_yy_yz[i] = -2.0 * g_0_zz_yy_yz[i] * a_exp + 4.0 * g_zz_zz_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_yy_zz[i] = -2.0 * g_0_zz_yy_zz[i] * a_exp + 4.0 * g_zz_zz_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (1284-1290)

    #pragma omp simd aligned(g_0_zz_yz_xx, g_0_zz_yz_xy, g_0_zz_yz_xz, g_0_zz_yz_yy, g_0_zz_yz_yz, g_0_zz_yz_zz, g_zz_0_0_0_0_zz_yz_xx, g_zz_0_0_0_0_zz_yz_xy, g_zz_0_0_0_0_zz_yz_xz, g_zz_0_0_0_0_zz_yz_yy, g_zz_0_0_0_0_zz_yz_yz, g_zz_0_0_0_0_zz_yz_zz, g_zz_zz_yz_xx, g_zz_zz_yz_xy, g_zz_zz_yz_xz, g_zz_zz_yz_yy, g_zz_zz_yz_yz, g_zz_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_zz_yz_xx[i] = -2.0 * g_0_zz_yz_xx[i] * a_exp + 4.0 * g_zz_zz_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_yz_xy[i] = -2.0 * g_0_zz_yz_xy[i] * a_exp + 4.0 * g_zz_zz_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_yz_xz[i] = -2.0 * g_0_zz_yz_xz[i] * a_exp + 4.0 * g_zz_zz_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_yz_yy[i] = -2.0 * g_0_zz_yz_yy[i] * a_exp + 4.0 * g_zz_zz_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_yz_yz[i] = -2.0 * g_0_zz_yz_yz[i] * a_exp + 4.0 * g_zz_zz_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_yz_zz[i] = -2.0 * g_0_zz_yz_zz[i] * a_exp + 4.0 * g_zz_zz_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (1290-1296)

    #pragma omp simd aligned(g_0_zz_zz_xx, g_0_zz_zz_xy, g_0_zz_zz_xz, g_0_zz_zz_yy, g_0_zz_zz_yz, g_0_zz_zz_zz, g_zz_0_0_0_0_zz_zz_xx, g_zz_0_0_0_0_zz_zz_xy, g_zz_0_0_0_0_zz_zz_xz, g_zz_0_0_0_0_zz_zz_yy, g_zz_0_0_0_0_zz_zz_yz, g_zz_0_0_0_0_zz_zz_zz, g_zz_zz_zz_xx, g_zz_zz_zz_xy, g_zz_zz_zz_xz, g_zz_zz_zz_yy, g_zz_zz_zz_yz, g_zz_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_zz_zz_xx[i] = -2.0 * g_0_zz_zz_xx[i] * a_exp + 4.0 * g_zz_zz_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_zz_xy[i] = -2.0 * g_0_zz_zz_xy[i] * a_exp + 4.0 * g_zz_zz_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_zz_xz[i] = -2.0 * g_0_zz_zz_xz[i] * a_exp + 4.0 * g_zz_zz_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_zz_yy[i] = -2.0 * g_0_zz_zz_yy[i] * a_exp + 4.0 * g_zz_zz_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_zz_yz[i] = -2.0 * g_0_zz_zz_yz[i] * a_exp + 4.0 * g_zz_zz_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_zz_zz[i] = -2.0 * g_0_zz_zz_zz[i] * a_exp + 4.0 * g_zz_zz_zz_zz[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

