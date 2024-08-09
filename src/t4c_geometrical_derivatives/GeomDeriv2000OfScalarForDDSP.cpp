#include "GeomDeriv2000OfScalarForDDSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_ddsp_0(CSimdArray<double>& buffer_2000_ddsp,
                     const CSimdArray<double>& buffer_sdsp,
                     const CSimdArray<double>& buffer_ddsp,
                     const CSimdArray<double>& buffer_gdsp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_ddsp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sdsp

    auto g_0_xx_0_x = buffer_sdsp[0];

    auto g_0_xx_0_y = buffer_sdsp[1];

    auto g_0_xx_0_z = buffer_sdsp[2];

    auto g_0_xy_0_x = buffer_sdsp[3];

    auto g_0_xy_0_y = buffer_sdsp[4];

    auto g_0_xy_0_z = buffer_sdsp[5];

    auto g_0_xz_0_x = buffer_sdsp[6];

    auto g_0_xz_0_y = buffer_sdsp[7];

    auto g_0_xz_0_z = buffer_sdsp[8];

    auto g_0_yy_0_x = buffer_sdsp[9];

    auto g_0_yy_0_y = buffer_sdsp[10];

    auto g_0_yy_0_z = buffer_sdsp[11];

    auto g_0_yz_0_x = buffer_sdsp[12];

    auto g_0_yz_0_y = buffer_sdsp[13];

    auto g_0_yz_0_z = buffer_sdsp[14];

    auto g_0_zz_0_x = buffer_sdsp[15];

    auto g_0_zz_0_y = buffer_sdsp[16];

    auto g_0_zz_0_z = buffer_sdsp[17];

    /// Set up components of auxilary buffer : buffer_ddsp

    auto g_xx_xx_0_x = buffer_ddsp[0];

    auto g_xx_xx_0_y = buffer_ddsp[1];

    auto g_xx_xx_0_z = buffer_ddsp[2];

    auto g_xx_xy_0_x = buffer_ddsp[3];

    auto g_xx_xy_0_y = buffer_ddsp[4];

    auto g_xx_xy_0_z = buffer_ddsp[5];

    auto g_xx_xz_0_x = buffer_ddsp[6];

    auto g_xx_xz_0_y = buffer_ddsp[7];

    auto g_xx_xz_0_z = buffer_ddsp[8];

    auto g_xx_yy_0_x = buffer_ddsp[9];

    auto g_xx_yy_0_y = buffer_ddsp[10];

    auto g_xx_yy_0_z = buffer_ddsp[11];

    auto g_xx_yz_0_x = buffer_ddsp[12];

    auto g_xx_yz_0_y = buffer_ddsp[13];

    auto g_xx_yz_0_z = buffer_ddsp[14];

    auto g_xx_zz_0_x = buffer_ddsp[15];

    auto g_xx_zz_0_y = buffer_ddsp[16];

    auto g_xx_zz_0_z = buffer_ddsp[17];

    auto g_xy_xx_0_x = buffer_ddsp[18];

    auto g_xy_xx_0_y = buffer_ddsp[19];

    auto g_xy_xx_0_z = buffer_ddsp[20];

    auto g_xy_xy_0_x = buffer_ddsp[21];

    auto g_xy_xy_0_y = buffer_ddsp[22];

    auto g_xy_xy_0_z = buffer_ddsp[23];

    auto g_xy_xz_0_x = buffer_ddsp[24];

    auto g_xy_xz_0_y = buffer_ddsp[25];

    auto g_xy_xz_0_z = buffer_ddsp[26];

    auto g_xy_yy_0_x = buffer_ddsp[27];

    auto g_xy_yy_0_y = buffer_ddsp[28];

    auto g_xy_yy_0_z = buffer_ddsp[29];

    auto g_xy_yz_0_x = buffer_ddsp[30];

    auto g_xy_yz_0_y = buffer_ddsp[31];

    auto g_xy_yz_0_z = buffer_ddsp[32];

    auto g_xy_zz_0_x = buffer_ddsp[33];

    auto g_xy_zz_0_y = buffer_ddsp[34];

    auto g_xy_zz_0_z = buffer_ddsp[35];

    auto g_xz_xx_0_x = buffer_ddsp[36];

    auto g_xz_xx_0_y = buffer_ddsp[37];

    auto g_xz_xx_0_z = buffer_ddsp[38];

    auto g_xz_xy_0_x = buffer_ddsp[39];

    auto g_xz_xy_0_y = buffer_ddsp[40];

    auto g_xz_xy_0_z = buffer_ddsp[41];

    auto g_xz_xz_0_x = buffer_ddsp[42];

    auto g_xz_xz_0_y = buffer_ddsp[43];

    auto g_xz_xz_0_z = buffer_ddsp[44];

    auto g_xz_yy_0_x = buffer_ddsp[45];

    auto g_xz_yy_0_y = buffer_ddsp[46];

    auto g_xz_yy_0_z = buffer_ddsp[47];

    auto g_xz_yz_0_x = buffer_ddsp[48];

    auto g_xz_yz_0_y = buffer_ddsp[49];

    auto g_xz_yz_0_z = buffer_ddsp[50];

    auto g_xz_zz_0_x = buffer_ddsp[51];

    auto g_xz_zz_0_y = buffer_ddsp[52];

    auto g_xz_zz_0_z = buffer_ddsp[53];

    auto g_yy_xx_0_x = buffer_ddsp[54];

    auto g_yy_xx_0_y = buffer_ddsp[55];

    auto g_yy_xx_0_z = buffer_ddsp[56];

    auto g_yy_xy_0_x = buffer_ddsp[57];

    auto g_yy_xy_0_y = buffer_ddsp[58];

    auto g_yy_xy_0_z = buffer_ddsp[59];

    auto g_yy_xz_0_x = buffer_ddsp[60];

    auto g_yy_xz_0_y = buffer_ddsp[61];

    auto g_yy_xz_0_z = buffer_ddsp[62];

    auto g_yy_yy_0_x = buffer_ddsp[63];

    auto g_yy_yy_0_y = buffer_ddsp[64];

    auto g_yy_yy_0_z = buffer_ddsp[65];

    auto g_yy_yz_0_x = buffer_ddsp[66];

    auto g_yy_yz_0_y = buffer_ddsp[67];

    auto g_yy_yz_0_z = buffer_ddsp[68];

    auto g_yy_zz_0_x = buffer_ddsp[69];

    auto g_yy_zz_0_y = buffer_ddsp[70];

    auto g_yy_zz_0_z = buffer_ddsp[71];

    auto g_yz_xx_0_x = buffer_ddsp[72];

    auto g_yz_xx_0_y = buffer_ddsp[73];

    auto g_yz_xx_0_z = buffer_ddsp[74];

    auto g_yz_xy_0_x = buffer_ddsp[75];

    auto g_yz_xy_0_y = buffer_ddsp[76];

    auto g_yz_xy_0_z = buffer_ddsp[77];

    auto g_yz_xz_0_x = buffer_ddsp[78];

    auto g_yz_xz_0_y = buffer_ddsp[79];

    auto g_yz_xz_0_z = buffer_ddsp[80];

    auto g_yz_yy_0_x = buffer_ddsp[81];

    auto g_yz_yy_0_y = buffer_ddsp[82];

    auto g_yz_yy_0_z = buffer_ddsp[83];

    auto g_yz_yz_0_x = buffer_ddsp[84];

    auto g_yz_yz_0_y = buffer_ddsp[85];

    auto g_yz_yz_0_z = buffer_ddsp[86];

    auto g_yz_zz_0_x = buffer_ddsp[87];

    auto g_yz_zz_0_y = buffer_ddsp[88];

    auto g_yz_zz_0_z = buffer_ddsp[89];

    auto g_zz_xx_0_x = buffer_ddsp[90];

    auto g_zz_xx_0_y = buffer_ddsp[91];

    auto g_zz_xx_0_z = buffer_ddsp[92];

    auto g_zz_xy_0_x = buffer_ddsp[93];

    auto g_zz_xy_0_y = buffer_ddsp[94];

    auto g_zz_xy_0_z = buffer_ddsp[95];

    auto g_zz_xz_0_x = buffer_ddsp[96];

    auto g_zz_xz_0_y = buffer_ddsp[97];

    auto g_zz_xz_0_z = buffer_ddsp[98];

    auto g_zz_yy_0_x = buffer_ddsp[99];

    auto g_zz_yy_0_y = buffer_ddsp[100];

    auto g_zz_yy_0_z = buffer_ddsp[101];

    auto g_zz_yz_0_x = buffer_ddsp[102];

    auto g_zz_yz_0_y = buffer_ddsp[103];

    auto g_zz_yz_0_z = buffer_ddsp[104];

    auto g_zz_zz_0_x = buffer_ddsp[105];

    auto g_zz_zz_0_y = buffer_ddsp[106];

    auto g_zz_zz_0_z = buffer_ddsp[107];

    /// Set up components of auxilary buffer : buffer_gdsp

    auto g_xxxx_xx_0_x = buffer_gdsp[0];

    auto g_xxxx_xx_0_y = buffer_gdsp[1];

    auto g_xxxx_xx_0_z = buffer_gdsp[2];

    auto g_xxxx_xy_0_x = buffer_gdsp[3];

    auto g_xxxx_xy_0_y = buffer_gdsp[4];

    auto g_xxxx_xy_0_z = buffer_gdsp[5];

    auto g_xxxx_xz_0_x = buffer_gdsp[6];

    auto g_xxxx_xz_0_y = buffer_gdsp[7];

    auto g_xxxx_xz_0_z = buffer_gdsp[8];

    auto g_xxxx_yy_0_x = buffer_gdsp[9];

    auto g_xxxx_yy_0_y = buffer_gdsp[10];

    auto g_xxxx_yy_0_z = buffer_gdsp[11];

    auto g_xxxx_yz_0_x = buffer_gdsp[12];

    auto g_xxxx_yz_0_y = buffer_gdsp[13];

    auto g_xxxx_yz_0_z = buffer_gdsp[14];

    auto g_xxxx_zz_0_x = buffer_gdsp[15];

    auto g_xxxx_zz_0_y = buffer_gdsp[16];

    auto g_xxxx_zz_0_z = buffer_gdsp[17];

    auto g_xxxy_xx_0_x = buffer_gdsp[18];

    auto g_xxxy_xx_0_y = buffer_gdsp[19];

    auto g_xxxy_xx_0_z = buffer_gdsp[20];

    auto g_xxxy_xy_0_x = buffer_gdsp[21];

    auto g_xxxy_xy_0_y = buffer_gdsp[22];

    auto g_xxxy_xy_0_z = buffer_gdsp[23];

    auto g_xxxy_xz_0_x = buffer_gdsp[24];

    auto g_xxxy_xz_0_y = buffer_gdsp[25];

    auto g_xxxy_xz_0_z = buffer_gdsp[26];

    auto g_xxxy_yy_0_x = buffer_gdsp[27];

    auto g_xxxy_yy_0_y = buffer_gdsp[28];

    auto g_xxxy_yy_0_z = buffer_gdsp[29];

    auto g_xxxy_yz_0_x = buffer_gdsp[30];

    auto g_xxxy_yz_0_y = buffer_gdsp[31];

    auto g_xxxy_yz_0_z = buffer_gdsp[32];

    auto g_xxxy_zz_0_x = buffer_gdsp[33];

    auto g_xxxy_zz_0_y = buffer_gdsp[34];

    auto g_xxxy_zz_0_z = buffer_gdsp[35];

    auto g_xxxz_xx_0_x = buffer_gdsp[36];

    auto g_xxxz_xx_0_y = buffer_gdsp[37];

    auto g_xxxz_xx_0_z = buffer_gdsp[38];

    auto g_xxxz_xy_0_x = buffer_gdsp[39];

    auto g_xxxz_xy_0_y = buffer_gdsp[40];

    auto g_xxxz_xy_0_z = buffer_gdsp[41];

    auto g_xxxz_xz_0_x = buffer_gdsp[42];

    auto g_xxxz_xz_0_y = buffer_gdsp[43];

    auto g_xxxz_xz_0_z = buffer_gdsp[44];

    auto g_xxxz_yy_0_x = buffer_gdsp[45];

    auto g_xxxz_yy_0_y = buffer_gdsp[46];

    auto g_xxxz_yy_0_z = buffer_gdsp[47];

    auto g_xxxz_yz_0_x = buffer_gdsp[48];

    auto g_xxxz_yz_0_y = buffer_gdsp[49];

    auto g_xxxz_yz_0_z = buffer_gdsp[50];

    auto g_xxxz_zz_0_x = buffer_gdsp[51];

    auto g_xxxz_zz_0_y = buffer_gdsp[52];

    auto g_xxxz_zz_0_z = buffer_gdsp[53];

    auto g_xxyy_xx_0_x = buffer_gdsp[54];

    auto g_xxyy_xx_0_y = buffer_gdsp[55];

    auto g_xxyy_xx_0_z = buffer_gdsp[56];

    auto g_xxyy_xy_0_x = buffer_gdsp[57];

    auto g_xxyy_xy_0_y = buffer_gdsp[58];

    auto g_xxyy_xy_0_z = buffer_gdsp[59];

    auto g_xxyy_xz_0_x = buffer_gdsp[60];

    auto g_xxyy_xz_0_y = buffer_gdsp[61];

    auto g_xxyy_xz_0_z = buffer_gdsp[62];

    auto g_xxyy_yy_0_x = buffer_gdsp[63];

    auto g_xxyy_yy_0_y = buffer_gdsp[64];

    auto g_xxyy_yy_0_z = buffer_gdsp[65];

    auto g_xxyy_yz_0_x = buffer_gdsp[66];

    auto g_xxyy_yz_0_y = buffer_gdsp[67];

    auto g_xxyy_yz_0_z = buffer_gdsp[68];

    auto g_xxyy_zz_0_x = buffer_gdsp[69];

    auto g_xxyy_zz_0_y = buffer_gdsp[70];

    auto g_xxyy_zz_0_z = buffer_gdsp[71];

    auto g_xxyz_xx_0_x = buffer_gdsp[72];

    auto g_xxyz_xx_0_y = buffer_gdsp[73];

    auto g_xxyz_xx_0_z = buffer_gdsp[74];

    auto g_xxyz_xy_0_x = buffer_gdsp[75];

    auto g_xxyz_xy_0_y = buffer_gdsp[76];

    auto g_xxyz_xy_0_z = buffer_gdsp[77];

    auto g_xxyz_xz_0_x = buffer_gdsp[78];

    auto g_xxyz_xz_0_y = buffer_gdsp[79];

    auto g_xxyz_xz_0_z = buffer_gdsp[80];

    auto g_xxyz_yy_0_x = buffer_gdsp[81];

    auto g_xxyz_yy_0_y = buffer_gdsp[82];

    auto g_xxyz_yy_0_z = buffer_gdsp[83];

    auto g_xxyz_yz_0_x = buffer_gdsp[84];

    auto g_xxyz_yz_0_y = buffer_gdsp[85];

    auto g_xxyz_yz_0_z = buffer_gdsp[86];

    auto g_xxyz_zz_0_x = buffer_gdsp[87];

    auto g_xxyz_zz_0_y = buffer_gdsp[88];

    auto g_xxyz_zz_0_z = buffer_gdsp[89];

    auto g_xxzz_xx_0_x = buffer_gdsp[90];

    auto g_xxzz_xx_0_y = buffer_gdsp[91];

    auto g_xxzz_xx_0_z = buffer_gdsp[92];

    auto g_xxzz_xy_0_x = buffer_gdsp[93];

    auto g_xxzz_xy_0_y = buffer_gdsp[94];

    auto g_xxzz_xy_0_z = buffer_gdsp[95];

    auto g_xxzz_xz_0_x = buffer_gdsp[96];

    auto g_xxzz_xz_0_y = buffer_gdsp[97];

    auto g_xxzz_xz_0_z = buffer_gdsp[98];

    auto g_xxzz_yy_0_x = buffer_gdsp[99];

    auto g_xxzz_yy_0_y = buffer_gdsp[100];

    auto g_xxzz_yy_0_z = buffer_gdsp[101];

    auto g_xxzz_yz_0_x = buffer_gdsp[102];

    auto g_xxzz_yz_0_y = buffer_gdsp[103];

    auto g_xxzz_yz_0_z = buffer_gdsp[104];

    auto g_xxzz_zz_0_x = buffer_gdsp[105];

    auto g_xxzz_zz_0_y = buffer_gdsp[106];

    auto g_xxzz_zz_0_z = buffer_gdsp[107];

    auto g_xyyy_xx_0_x = buffer_gdsp[108];

    auto g_xyyy_xx_0_y = buffer_gdsp[109];

    auto g_xyyy_xx_0_z = buffer_gdsp[110];

    auto g_xyyy_xy_0_x = buffer_gdsp[111];

    auto g_xyyy_xy_0_y = buffer_gdsp[112];

    auto g_xyyy_xy_0_z = buffer_gdsp[113];

    auto g_xyyy_xz_0_x = buffer_gdsp[114];

    auto g_xyyy_xz_0_y = buffer_gdsp[115];

    auto g_xyyy_xz_0_z = buffer_gdsp[116];

    auto g_xyyy_yy_0_x = buffer_gdsp[117];

    auto g_xyyy_yy_0_y = buffer_gdsp[118];

    auto g_xyyy_yy_0_z = buffer_gdsp[119];

    auto g_xyyy_yz_0_x = buffer_gdsp[120];

    auto g_xyyy_yz_0_y = buffer_gdsp[121];

    auto g_xyyy_yz_0_z = buffer_gdsp[122];

    auto g_xyyy_zz_0_x = buffer_gdsp[123];

    auto g_xyyy_zz_0_y = buffer_gdsp[124];

    auto g_xyyy_zz_0_z = buffer_gdsp[125];

    auto g_xyyz_xx_0_x = buffer_gdsp[126];

    auto g_xyyz_xx_0_y = buffer_gdsp[127];

    auto g_xyyz_xx_0_z = buffer_gdsp[128];

    auto g_xyyz_xy_0_x = buffer_gdsp[129];

    auto g_xyyz_xy_0_y = buffer_gdsp[130];

    auto g_xyyz_xy_0_z = buffer_gdsp[131];

    auto g_xyyz_xz_0_x = buffer_gdsp[132];

    auto g_xyyz_xz_0_y = buffer_gdsp[133];

    auto g_xyyz_xz_0_z = buffer_gdsp[134];

    auto g_xyyz_yy_0_x = buffer_gdsp[135];

    auto g_xyyz_yy_0_y = buffer_gdsp[136];

    auto g_xyyz_yy_0_z = buffer_gdsp[137];

    auto g_xyyz_yz_0_x = buffer_gdsp[138];

    auto g_xyyz_yz_0_y = buffer_gdsp[139];

    auto g_xyyz_yz_0_z = buffer_gdsp[140];

    auto g_xyyz_zz_0_x = buffer_gdsp[141];

    auto g_xyyz_zz_0_y = buffer_gdsp[142];

    auto g_xyyz_zz_0_z = buffer_gdsp[143];

    auto g_xyzz_xx_0_x = buffer_gdsp[144];

    auto g_xyzz_xx_0_y = buffer_gdsp[145];

    auto g_xyzz_xx_0_z = buffer_gdsp[146];

    auto g_xyzz_xy_0_x = buffer_gdsp[147];

    auto g_xyzz_xy_0_y = buffer_gdsp[148];

    auto g_xyzz_xy_0_z = buffer_gdsp[149];

    auto g_xyzz_xz_0_x = buffer_gdsp[150];

    auto g_xyzz_xz_0_y = buffer_gdsp[151];

    auto g_xyzz_xz_0_z = buffer_gdsp[152];

    auto g_xyzz_yy_0_x = buffer_gdsp[153];

    auto g_xyzz_yy_0_y = buffer_gdsp[154];

    auto g_xyzz_yy_0_z = buffer_gdsp[155];

    auto g_xyzz_yz_0_x = buffer_gdsp[156];

    auto g_xyzz_yz_0_y = buffer_gdsp[157];

    auto g_xyzz_yz_0_z = buffer_gdsp[158];

    auto g_xyzz_zz_0_x = buffer_gdsp[159];

    auto g_xyzz_zz_0_y = buffer_gdsp[160];

    auto g_xyzz_zz_0_z = buffer_gdsp[161];

    auto g_xzzz_xx_0_x = buffer_gdsp[162];

    auto g_xzzz_xx_0_y = buffer_gdsp[163];

    auto g_xzzz_xx_0_z = buffer_gdsp[164];

    auto g_xzzz_xy_0_x = buffer_gdsp[165];

    auto g_xzzz_xy_0_y = buffer_gdsp[166];

    auto g_xzzz_xy_0_z = buffer_gdsp[167];

    auto g_xzzz_xz_0_x = buffer_gdsp[168];

    auto g_xzzz_xz_0_y = buffer_gdsp[169];

    auto g_xzzz_xz_0_z = buffer_gdsp[170];

    auto g_xzzz_yy_0_x = buffer_gdsp[171];

    auto g_xzzz_yy_0_y = buffer_gdsp[172];

    auto g_xzzz_yy_0_z = buffer_gdsp[173];

    auto g_xzzz_yz_0_x = buffer_gdsp[174];

    auto g_xzzz_yz_0_y = buffer_gdsp[175];

    auto g_xzzz_yz_0_z = buffer_gdsp[176];

    auto g_xzzz_zz_0_x = buffer_gdsp[177];

    auto g_xzzz_zz_0_y = buffer_gdsp[178];

    auto g_xzzz_zz_0_z = buffer_gdsp[179];

    auto g_yyyy_xx_0_x = buffer_gdsp[180];

    auto g_yyyy_xx_0_y = buffer_gdsp[181];

    auto g_yyyy_xx_0_z = buffer_gdsp[182];

    auto g_yyyy_xy_0_x = buffer_gdsp[183];

    auto g_yyyy_xy_0_y = buffer_gdsp[184];

    auto g_yyyy_xy_0_z = buffer_gdsp[185];

    auto g_yyyy_xz_0_x = buffer_gdsp[186];

    auto g_yyyy_xz_0_y = buffer_gdsp[187];

    auto g_yyyy_xz_0_z = buffer_gdsp[188];

    auto g_yyyy_yy_0_x = buffer_gdsp[189];

    auto g_yyyy_yy_0_y = buffer_gdsp[190];

    auto g_yyyy_yy_0_z = buffer_gdsp[191];

    auto g_yyyy_yz_0_x = buffer_gdsp[192];

    auto g_yyyy_yz_0_y = buffer_gdsp[193];

    auto g_yyyy_yz_0_z = buffer_gdsp[194];

    auto g_yyyy_zz_0_x = buffer_gdsp[195];

    auto g_yyyy_zz_0_y = buffer_gdsp[196];

    auto g_yyyy_zz_0_z = buffer_gdsp[197];

    auto g_yyyz_xx_0_x = buffer_gdsp[198];

    auto g_yyyz_xx_0_y = buffer_gdsp[199];

    auto g_yyyz_xx_0_z = buffer_gdsp[200];

    auto g_yyyz_xy_0_x = buffer_gdsp[201];

    auto g_yyyz_xy_0_y = buffer_gdsp[202];

    auto g_yyyz_xy_0_z = buffer_gdsp[203];

    auto g_yyyz_xz_0_x = buffer_gdsp[204];

    auto g_yyyz_xz_0_y = buffer_gdsp[205];

    auto g_yyyz_xz_0_z = buffer_gdsp[206];

    auto g_yyyz_yy_0_x = buffer_gdsp[207];

    auto g_yyyz_yy_0_y = buffer_gdsp[208];

    auto g_yyyz_yy_0_z = buffer_gdsp[209];

    auto g_yyyz_yz_0_x = buffer_gdsp[210];

    auto g_yyyz_yz_0_y = buffer_gdsp[211];

    auto g_yyyz_yz_0_z = buffer_gdsp[212];

    auto g_yyyz_zz_0_x = buffer_gdsp[213];

    auto g_yyyz_zz_0_y = buffer_gdsp[214];

    auto g_yyyz_zz_0_z = buffer_gdsp[215];

    auto g_yyzz_xx_0_x = buffer_gdsp[216];

    auto g_yyzz_xx_0_y = buffer_gdsp[217];

    auto g_yyzz_xx_0_z = buffer_gdsp[218];

    auto g_yyzz_xy_0_x = buffer_gdsp[219];

    auto g_yyzz_xy_0_y = buffer_gdsp[220];

    auto g_yyzz_xy_0_z = buffer_gdsp[221];

    auto g_yyzz_xz_0_x = buffer_gdsp[222];

    auto g_yyzz_xz_0_y = buffer_gdsp[223];

    auto g_yyzz_xz_0_z = buffer_gdsp[224];

    auto g_yyzz_yy_0_x = buffer_gdsp[225];

    auto g_yyzz_yy_0_y = buffer_gdsp[226];

    auto g_yyzz_yy_0_z = buffer_gdsp[227];

    auto g_yyzz_yz_0_x = buffer_gdsp[228];

    auto g_yyzz_yz_0_y = buffer_gdsp[229];

    auto g_yyzz_yz_0_z = buffer_gdsp[230];

    auto g_yyzz_zz_0_x = buffer_gdsp[231];

    auto g_yyzz_zz_0_y = buffer_gdsp[232];

    auto g_yyzz_zz_0_z = buffer_gdsp[233];

    auto g_yzzz_xx_0_x = buffer_gdsp[234];

    auto g_yzzz_xx_0_y = buffer_gdsp[235];

    auto g_yzzz_xx_0_z = buffer_gdsp[236];

    auto g_yzzz_xy_0_x = buffer_gdsp[237];

    auto g_yzzz_xy_0_y = buffer_gdsp[238];

    auto g_yzzz_xy_0_z = buffer_gdsp[239];

    auto g_yzzz_xz_0_x = buffer_gdsp[240];

    auto g_yzzz_xz_0_y = buffer_gdsp[241];

    auto g_yzzz_xz_0_z = buffer_gdsp[242];

    auto g_yzzz_yy_0_x = buffer_gdsp[243];

    auto g_yzzz_yy_0_y = buffer_gdsp[244];

    auto g_yzzz_yy_0_z = buffer_gdsp[245];

    auto g_yzzz_yz_0_x = buffer_gdsp[246];

    auto g_yzzz_yz_0_y = buffer_gdsp[247];

    auto g_yzzz_yz_0_z = buffer_gdsp[248];

    auto g_yzzz_zz_0_x = buffer_gdsp[249];

    auto g_yzzz_zz_0_y = buffer_gdsp[250];

    auto g_yzzz_zz_0_z = buffer_gdsp[251];

    auto g_zzzz_xx_0_x = buffer_gdsp[252];

    auto g_zzzz_xx_0_y = buffer_gdsp[253];

    auto g_zzzz_xx_0_z = buffer_gdsp[254];

    auto g_zzzz_xy_0_x = buffer_gdsp[255];

    auto g_zzzz_xy_0_y = buffer_gdsp[256];

    auto g_zzzz_xy_0_z = buffer_gdsp[257];

    auto g_zzzz_xz_0_x = buffer_gdsp[258];

    auto g_zzzz_xz_0_y = buffer_gdsp[259];

    auto g_zzzz_xz_0_z = buffer_gdsp[260];

    auto g_zzzz_yy_0_x = buffer_gdsp[261];

    auto g_zzzz_yy_0_y = buffer_gdsp[262];

    auto g_zzzz_yy_0_z = buffer_gdsp[263];

    auto g_zzzz_yz_0_x = buffer_gdsp[264];

    auto g_zzzz_yz_0_y = buffer_gdsp[265];

    auto g_zzzz_yz_0_z = buffer_gdsp[266];

    auto g_zzzz_zz_0_x = buffer_gdsp[267];

    auto g_zzzz_zz_0_y = buffer_gdsp[268];

    auto g_zzzz_zz_0_z = buffer_gdsp[269];

    /// Set up components of integrals buffer : buffer_2000_ddsp

    auto g_xx_0_0_0_xx_xx_0_x = buffer_2000_ddsp[0];

    auto g_xx_0_0_0_xx_xx_0_y = buffer_2000_ddsp[1];

    auto g_xx_0_0_0_xx_xx_0_z = buffer_2000_ddsp[2];

    auto g_xx_0_0_0_xx_xy_0_x = buffer_2000_ddsp[3];

    auto g_xx_0_0_0_xx_xy_0_y = buffer_2000_ddsp[4];

    auto g_xx_0_0_0_xx_xy_0_z = buffer_2000_ddsp[5];

    auto g_xx_0_0_0_xx_xz_0_x = buffer_2000_ddsp[6];

    auto g_xx_0_0_0_xx_xz_0_y = buffer_2000_ddsp[7];

    auto g_xx_0_0_0_xx_xz_0_z = buffer_2000_ddsp[8];

    auto g_xx_0_0_0_xx_yy_0_x = buffer_2000_ddsp[9];

    auto g_xx_0_0_0_xx_yy_0_y = buffer_2000_ddsp[10];

    auto g_xx_0_0_0_xx_yy_0_z = buffer_2000_ddsp[11];

    auto g_xx_0_0_0_xx_yz_0_x = buffer_2000_ddsp[12];

    auto g_xx_0_0_0_xx_yz_0_y = buffer_2000_ddsp[13];

    auto g_xx_0_0_0_xx_yz_0_z = buffer_2000_ddsp[14];

    auto g_xx_0_0_0_xx_zz_0_x = buffer_2000_ddsp[15];

    auto g_xx_0_0_0_xx_zz_0_y = buffer_2000_ddsp[16];

    auto g_xx_0_0_0_xx_zz_0_z = buffer_2000_ddsp[17];

    auto g_xx_0_0_0_xy_xx_0_x = buffer_2000_ddsp[18];

    auto g_xx_0_0_0_xy_xx_0_y = buffer_2000_ddsp[19];

    auto g_xx_0_0_0_xy_xx_0_z = buffer_2000_ddsp[20];

    auto g_xx_0_0_0_xy_xy_0_x = buffer_2000_ddsp[21];

    auto g_xx_0_0_0_xy_xy_0_y = buffer_2000_ddsp[22];

    auto g_xx_0_0_0_xy_xy_0_z = buffer_2000_ddsp[23];

    auto g_xx_0_0_0_xy_xz_0_x = buffer_2000_ddsp[24];

    auto g_xx_0_0_0_xy_xz_0_y = buffer_2000_ddsp[25];

    auto g_xx_0_0_0_xy_xz_0_z = buffer_2000_ddsp[26];

    auto g_xx_0_0_0_xy_yy_0_x = buffer_2000_ddsp[27];

    auto g_xx_0_0_0_xy_yy_0_y = buffer_2000_ddsp[28];

    auto g_xx_0_0_0_xy_yy_0_z = buffer_2000_ddsp[29];

    auto g_xx_0_0_0_xy_yz_0_x = buffer_2000_ddsp[30];

    auto g_xx_0_0_0_xy_yz_0_y = buffer_2000_ddsp[31];

    auto g_xx_0_0_0_xy_yz_0_z = buffer_2000_ddsp[32];

    auto g_xx_0_0_0_xy_zz_0_x = buffer_2000_ddsp[33];

    auto g_xx_0_0_0_xy_zz_0_y = buffer_2000_ddsp[34];

    auto g_xx_0_0_0_xy_zz_0_z = buffer_2000_ddsp[35];

    auto g_xx_0_0_0_xz_xx_0_x = buffer_2000_ddsp[36];

    auto g_xx_0_0_0_xz_xx_0_y = buffer_2000_ddsp[37];

    auto g_xx_0_0_0_xz_xx_0_z = buffer_2000_ddsp[38];

    auto g_xx_0_0_0_xz_xy_0_x = buffer_2000_ddsp[39];

    auto g_xx_0_0_0_xz_xy_0_y = buffer_2000_ddsp[40];

    auto g_xx_0_0_0_xz_xy_0_z = buffer_2000_ddsp[41];

    auto g_xx_0_0_0_xz_xz_0_x = buffer_2000_ddsp[42];

    auto g_xx_0_0_0_xz_xz_0_y = buffer_2000_ddsp[43];

    auto g_xx_0_0_0_xz_xz_0_z = buffer_2000_ddsp[44];

    auto g_xx_0_0_0_xz_yy_0_x = buffer_2000_ddsp[45];

    auto g_xx_0_0_0_xz_yy_0_y = buffer_2000_ddsp[46];

    auto g_xx_0_0_0_xz_yy_0_z = buffer_2000_ddsp[47];

    auto g_xx_0_0_0_xz_yz_0_x = buffer_2000_ddsp[48];

    auto g_xx_0_0_0_xz_yz_0_y = buffer_2000_ddsp[49];

    auto g_xx_0_0_0_xz_yz_0_z = buffer_2000_ddsp[50];

    auto g_xx_0_0_0_xz_zz_0_x = buffer_2000_ddsp[51];

    auto g_xx_0_0_0_xz_zz_0_y = buffer_2000_ddsp[52];

    auto g_xx_0_0_0_xz_zz_0_z = buffer_2000_ddsp[53];

    auto g_xx_0_0_0_yy_xx_0_x = buffer_2000_ddsp[54];

    auto g_xx_0_0_0_yy_xx_0_y = buffer_2000_ddsp[55];

    auto g_xx_0_0_0_yy_xx_0_z = buffer_2000_ddsp[56];

    auto g_xx_0_0_0_yy_xy_0_x = buffer_2000_ddsp[57];

    auto g_xx_0_0_0_yy_xy_0_y = buffer_2000_ddsp[58];

    auto g_xx_0_0_0_yy_xy_0_z = buffer_2000_ddsp[59];

    auto g_xx_0_0_0_yy_xz_0_x = buffer_2000_ddsp[60];

    auto g_xx_0_0_0_yy_xz_0_y = buffer_2000_ddsp[61];

    auto g_xx_0_0_0_yy_xz_0_z = buffer_2000_ddsp[62];

    auto g_xx_0_0_0_yy_yy_0_x = buffer_2000_ddsp[63];

    auto g_xx_0_0_0_yy_yy_0_y = buffer_2000_ddsp[64];

    auto g_xx_0_0_0_yy_yy_0_z = buffer_2000_ddsp[65];

    auto g_xx_0_0_0_yy_yz_0_x = buffer_2000_ddsp[66];

    auto g_xx_0_0_0_yy_yz_0_y = buffer_2000_ddsp[67];

    auto g_xx_0_0_0_yy_yz_0_z = buffer_2000_ddsp[68];

    auto g_xx_0_0_0_yy_zz_0_x = buffer_2000_ddsp[69];

    auto g_xx_0_0_0_yy_zz_0_y = buffer_2000_ddsp[70];

    auto g_xx_0_0_0_yy_zz_0_z = buffer_2000_ddsp[71];

    auto g_xx_0_0_0_yz_xx_0_x = buffer_2000_ddsp[72];

    auto g_xx_0_0_0_yz_xx_0_y = buffer_2000_ddsp[73];

    auto g_xx_0_0_0_yz_xx_0_z = buffer_2000_ddsp[74];

    auto g_xx_0_0_0_yz_xy_0_x = buffer_2000_ddsp[75];

    auto g_xx_0_0_0_yz_xy_0_y = buffer_2000_ddsp[76];

    auto g_xx_0_0_0_yz_xy_0_z = buffer_2000_ddsp[77];

    auto g_xx_0_0_0_yz_xz_0_x = buffer_2000_ddsp[78];

    auto g_xx_0_0_0_yz_xz_0_y = buffer_2000_ddsp[79];

    auto g_xx_0_0_0_yz_xz_0_z = buffer_2000_ddsp[80];

    auto g_xx_0_0_0_yz_yy_0_x = buffer_2000_ddsp[81];

    auto g_xx_0_0_0_yz_yy_0_y = buffer_2000_ddsp[82];

    auto g_xx_0_0_0_yz_yy_0_z = buffer_2000_ddsp[83];

    auto g_xx_0_0_0_yz_yz_0_x = buffer_2000_ddsp[84];

    auto g_xx_0_0_0_yz_yz_0_y = buffer_2000_ddsp[85];

    auto g_xx_0_0_0_yz_yz_0_z = buffer_2000_ddsp[86];

    auto g_xx_0_0_0_yz_zz_0_x = buffer_2000_ddsp[87];

    auto g_xx_0_0_0_yz_zz_0_y = buffer_2000_ddsp[88];

    auto g_xx_0_0_0_yz_zz_0_z = buffer_2000_ddsp[89];

    auto g_xx_0_0_0_zz_xx_0_x = buffer_2000_ddsp[90];

    auto g_xx_0_0_0_zz_xx_0_y = buffer_2000_ddsp[91];

    auto g_xx_0_0_0_zz_xx_0_z = buffer_2000_ddsp[92];

    auto g_xx_0_0_0_zz_xy_0_x = buffer_2000_ddsp[93];

    auto g_xx_0_0_0_zz_xy_0_y = buffer_2000_ddsp[94];

    auto g_xx_0_0_0_zz_xy_0_z = buffer_2000_ddsp[95];

    auto g_xx_0_0_0_zz_xz_0_x = buffer_2000_ddsp[96];

    auto g_xx_0_0_0_zz_xz_0_y = buffer_2000_ddsp[97];

    auto g_xx_0_0_0_zz_xz_0_z = buffer_2000_ddsp[98];

    auto g_xx_0_0_0_zz_yy_0_x = buffer_2000_ddsp[99];

    auto g_xx_0_0_0_zz_yy_0_y = buffer_2000_ddsp[100];

    auto g_xx_0_0_0_zz_yy_0_z = buffer_2000_ddsp[101];

    auto g_xx_0_0_0_zz_yz_0_x = buffer_2000_ddsp[102];

    auto g_xx_0_0_0_zz_yz_0_y = buffer_2000_ddsp[103];

    auto g_xx_0_0_0_zz_yz_0_z = buffer_2000_ddsp[104];

    auto g_xx_0_0_0_zz_zz_0_x = buffer_2000_ddsp[105];

    auto g_xx_0_0_0_zz_zz_0_y = buffer_2000_ddsp[106];

    auto g_xx_0_0_0_zz_zz_0_z = buffer_2000_ddsp[107];

    auto g_xy_0_0_0_xx_xx_0_x = buffer_2000_ddsp[108];

    auto g_xy_0_0_0_xx_xx_0_y = buffer_2000_ddsp[109];

    auto g_xy_0_0_0_xx_xx_0_z = buffer_2000_ddsp[110];

    auto g_xy_0_0_0_xx_xy_0_x = buffer_2000_ddsp[111];

    auto g_xy_0_0_0_xx_xy_0_y = buffer_2000_ddsp[112];

    auto g_xy_0_0_0_xx_xy_0_z = buffer_2000_ddsp[113];

    auto g_xy_0_0_0_xx_xz_0_x = buffer_2000_ddsp[114];

    auto g_xy_0_0_0_xx_xz_0_y = buffer_2000_ddsp[115];

    auto g_xy_0_0_0_xx_xz_0_z = buffer_2000_ddsp[116];

    auto g_xy_0_0_0_xx_yy_0_x = buffer_2000_ddsp[117];

    auto g_xy_0_0_0_xx_yy_0_y = buffer_2000_ddsp[118];

    auto g_xy_0_0_0_xx_yy_0_z = buffer_2000_ddsp[119];

    auto g_xy_0_0_0_xx_yz_0_x = buffer_2000_ddsp[120];

    auto g_xy_0_0_0_xx_yz_0_y = buffer_2000_ddsp[121];

    auto g_xy_0_0_0_xx_yz_0_z = buffer_2000_ddsp[122];

    auto g_xy_0_0_0_xx_zz_0_x = buffer_2000_ddsp[123];

    auto g_xy_0_0_0_xx_zz_0_y = buffer_2000_ddsp[124];

    auto g_xy_0_0_0_xx_zz_0_z = buffer_2000_ddsp[125];

    auto g_xy_0_0_0_xy_xx_0_x = buffer_2000_ddsp[126];

    auto g_xy_0_0_0_xy_xx_0_y = buffer_2000_ddsp[127];

    auto g_xy_0_0_0_xy_xx_0_z = buffer_2000_ddsp[128];

    auto g_xy_0_0_0_xy_xy_0_x = buffer_2000_ddsp[129];

    auto g_xy_0_0_0_xy_xy_0_y = buffer_2000_ddsp[130];

    auto g_xy_0_0_0_xy_xy_0_z = buffer_2000_ddsp[131];

    auto g_xy_0_0_0_xy_xz_0_x = buffer_2000_ddsp[132];

    auto g_xy_0_0_0_xy_xz_0_y = buffer_2000_ddsp[133];

    auto g_xy_0_0_0_xy_xz_0_z = buffer_2000_ddsp[134];

    auto g_xy_0_0_0_xy_yy_0_x = buffer_2000_ddsp[135];

    auto g_xy_0_0_0_xy_yy_0_y = buffer_2000_ddsp[136];

    auto g_xy_0_0_0_xy_yy_0_z = buffer_2000_ddsp[137];

    auto g_xy_0_0_0_xy_yz_0_x = buffer_2000_ddsp[138];

    auto g_xy_0_0_0_xy_yz_0_y = buffer_2000_ddsp[139];

    auto g_xy_0_0_0_xy_yz_0_z = buffer_2000_ddsp[140];

    auto g_xy_0_0_0_xy_zz_0_x = buffer_2000_ddsp[141];

    auto g_xy_0_0_0_xy_zz_0_y = buffer_2000_ddsp[142];

    auto g_xy_0_0_0_xy_zz_0_z = buffer_2000_ddsp[143];

    auto g_xy_0_0_0_xz_xx_0_x = buffer_2000_ddsp[144];

    auto g_xy_0_0_0_xz_xx_0_y = buffer_2000_ddsp[145];

    auto g_xy_0_0_0_xz_xx_0_z = buffer_2000_ddsp[146];

    auto g_xy_0_0_0_xz_xy_0_x = buffer_2000_ddsp[147];

    auto g_xy_0_0_0_xz_xy_0_y = buffer_2000_ddsp[148];

    auto g_xy_0_0_0_xz_xy_0_z = buffer_2000_ddsp[149];

    auto g_xy_0_0_0_xz_xz_0_x = buffer_2000_ddsp[150];

    auto g_xy_0_0_0_xz_xz_0_y = buffer_2000_ddsp[151];

    auto g_xy_0_0_0_xz_xz_0_z = buffer_2000_ddsp[152];

    auto g_xy_0_0_0_xz_yy_0_x = buffer_2000_ddsp[153];

    auto g_xy_0_0_0_xz_yy_0_y = buffer_2000_ddsp[154];

    auto g_xy_0_0_0_xz_yy_0_z = buffer_2000_ddsp[155];

    auto g_xy_0_0_0_xz_yz_0_x = buffer_2000_ddsp[156];

    auto g_xy_0_0_0_xz_yz_0_y = buffer_2000_ddsp[157];

    auto g_xy_0_0_0_xz_yz_0_z = buffer_2000_ddsp[158];

    auto g_xy_0_0_0_xz_zz_0_x = buffer_2000_ddsp[159];

    auto g_xy_0_0_0_xz_zz_0_y = buffer_2000_ddsp[160];

    auto g_xy_0_0_0_xz_zz_0_z = buffer_2000_ddsp[161];

    auto g_xy_0_0_0_yy_xx_0_x = buffer_2000_ddsp[162];

    auto g_xy_0_0_0_yy_xx_0_y = buffer_2000_ddsp[163];

    auto g_xy_0_0_0_yy_xx_0_z = buffer_2000_ddsp[164];

    auto g_xy_0_0_0_yy_xy_0_x = buffer_2000_ddsp[165];

    auto g_xy_0_0_0_yy_xy_0_y = buffer_2000_ddsp[166];

    auto g_xy_0_0_0_yy_xy_0_z = buffer_2000_ddsp[167];

    auto g_xy_0_0_0_yy_xz_0_x = buffer_2000_ddsp[168];

    auto g_xy_0_0_0_yy_xz_0_y = buffer_2000_ddsp[169];

    auto g_xy_0_0_0_yy_xz_0_z = buffer_2000_ddsp[170];

    auto g_xy_0_0_0_yy_yy_0_x = buffer_2000_ddsp[171];

    auto g_xy_0_0_0_yy_yy_0_y = buffer_2000_ddsp[172];

    auto g_xy_0_0_0_yy_yy_0_z = buffer_2000_ddsp[173];

    auto g_xy_0_0_0_yy_yz_0_x = buffer_2000_ddsp[174];

    auto g_xy_0_0_0_yy_yz_0_y = buffer_2000_ddsp[175];

    auto g_xy_0_0_0_yy_yz_0_z = buffer_2000_ddsp[176];

    auto g_xy_0_0_0_yy_zz_0_x = buffer_2000_ddsp[177];

    auto g_xy_0_0_0_yy_zz_0_y = buffer_2000_ddsp[178];

    auto g_xy_0_0_0_yy_zz_0_z = buffer_2000_ddsp[179];

    auto g_xy_0_0_0_yz_xx_0_x = buffer_2000_ddsp[180];

    auto g_xy_0_0_0_yz_xx_0_y = buffer_2000_ddsp[181];

    auto g_xy_0_0_0_yz_xx_0_z = buffer_2000_ddsp[182];

    auto g_xy_0_0_0_yz_xy_0_x = buffer_2000_ddsp[183];

    auto g_xy_0_0_0_yz_xy_0_y = buffer_2000_ddsp[184];

    auto g_xy_0_0_0_yz_xy_0_z = buffer_2000_ddsp[185];

    auto g_xy_0_0_0_yz_xz_0_x = buffer_2000_ddsp[186];

    auto g_xy_0_0_0_yz_xz_0_y = buffer_2000_ddsp[187];

    auto g_xy_0_0_0_yz_xz_0_z = buffer_2000_ddsp[188];

    auto g_xy_0_0_0_yz_yy_0_x = buffer_2000_ddsp[189];

    auto g_xy_0_0_0_yz_yy_0_y = buffer_2000_ddsp[190];

    auto g_xy_0_0_0_yz_yy_0_z = buffer_2000_ddsp[191];

    auto g_xy_0_0_0_yz_yz_0_x = buffer_2000_ddsp[192];

    auto g_xy_0_0_0_yz_yz_0_y = buffer_2000_ddsp[193];

    auto g_xy_0_0_0_yz_yz_0_z = buffer_2000_ddsp[194];

    auto g_xy_0_0_0_yz_zz_0_x = buffer_2000_ddsp[195];

    auto g_xy_0_0_0_yz_zz_0_y = buffer_2000_ddsp[196];

    auto g_xy_0_0_0_yz_zz_0_z = buffer_2000_ddsp[197];

    auto g_xy_0_0_0_zz_xx_0_x = buffer_2000_ddsp[198];

    auto g_xy_0_0_0_zz_xx_0_y = buffer_2000_ddsp[199];

    auto g_xy_0_0_0_zz_xx_0_z = buffer_2000_ddsp[200];

    auto g_xy_0_0_0_zz_xy_0_x = buffer_2000_ddsp[201];

    auto g_xy_0_0_0_zz_xy_0_y = buffer_2000_ddsp[202];

    auto g_xy_0_0_0_zz_xy_0_z = buffer_2000_ddsp[203];

    auto g_xy_0_0_0_zz_xz_0_x = buffer_2000_ddsp[204];

    auto g_xy_0_0_0_zz_xz_0_y = buffer_2000_ddsp[205];

    auto g_xy_0_0_0_zz_xz_0_z = buffer_2000_ddsp[206];

    auto g_xy_0_0_0_zz_yy_0_x = buffer_2000_ddsp[207];

    auto g_xy_0_0_0_zz_yy_0_y = buffer_2000_ddsp[208];

    auto g_xy_0_0_0_zz_yy_0_z = buffer_2000_ddsp[209];

    auto g_xy_0_0_0_zz_yz_0_x = buffer_2000_ddsp[210];

    auto g_xy_0_0_0_zz_yz_0_y = buffer_2000_ddsp[211];

    auto g_xy_0_0_0_zz_yz_0_z = buffer_2000_ddsp[212];

    auto g_xy_0_0_0_zz_zz_0_x = buffer_2000_ddsp[213];

    auto g_xy_0_0_0_zz_zz_0_y = buffer_2000_ddsp[214];

    auto g_xy_0_0_0_zz_zz_0_z = buffer_2000_ddsp[215];

    auto g_xz_0_0_0_xx_xx_0_x = buffer_2000_ddsp[216];

    auto g_xz_0_0_0_xx_xx_0_y = buffer_2000_ddsp[217];

    auto g_xz_0_0_0_xx_xx_0_z = buffer_2000_ddsp[218];

    auto g_xz_0_0_0_xx_xy_0_x = buffer_2000_ddsp[219];

    auto g_xz_0_0_0_xx_xy_0_y = buffer_2000_ddsp[220];

    auto g_xz_0_0_0_xx_xy_0_z = buffer_2000_ddsp[221];

    auto g_xz_0_0_0_xx_xz_0_x = buffer_2000_ddsp[222];

    auto g_xz_0_0_0_xx_xz_0_y = buffer_2000_ddsp[223];

    auto g_xz_0_0_0_xx_xz_0_z = buffer_2000_ddsp[224];

    auto g_xz_0_0_0_xx_yy_0_x = buffer_2000_ddsp[225];

    auto g_xz_0_0_0_xx_yy_0_y = buffer_2000_ddsp[226];

    auto g_xz_0_0_0_xx_yy_0_z = buffer_2000_ddsp[227];

    auto g_xz_0_0_0_xx_yz_0_x = buffer_2000_ddsp[228];

    auto g_xz_0_0_0_xx_yz_0_y = buffer_2000_ddsp[229];

    auto g_xz_0_0_0_xx_yz_0_z = buffer_2000_ddsp[230];

    auto g_xz_0_0_0_xx_zz_0_x = buffer_2000_ddsp[231];

    auto g_xz_0_0_0_xx_zz_0_y = buffer_2000_ddsp[232];

    auto g_xz_0_0_0_xx_zz_0_z = buffer_2000_ddsp[233];

    auto g_xz_0_0_0_xy_xx_0_x = buffer_2000_ddsp[234];

    auto g_xz_0_0_0_xy_xx_0_y = buffer_2000_ddsp[235];

    auto g_xz_0_0_0_xy_xx_0_z = buffer_2000_ddsp[236];

    auto g_xz_0_0_0_xy_xy_0_x = buffer_2000_ddsp[237];

    auto g_xz_0_0_0_xy_xy_0_y = buffer_2000_ddsp[238];

    auto g_xz_0_0_0_xy_xy_0_z = buffer_2000_ddsp[239];

    auto g_xz_0_0_0_xy_xz_0_x = buffer_2000_ddsp[240];

    auto g_xz_0_0_0_xy_xz_0_y = buffer_2000_ddsp[241];

    auto g_xz_0_0_0_xy_xz_0_z = buffer_2000_ddsp[242];

    auto g_xz_0_0_0_xy_yy_0_x = buffer_2000_ddsp[243];

    auto g_xz_0_0_0_xy_yy_0_y = buffer_2000_ddsp[244];

    auto g_xz_0_0_0_xy_yy_0_z = buffer_2000_ddsp[245];

    auto g_xz_0_0_0_xy_yz_0_x = buffer_2000_ddsp[246];

    auto g_xz_0_0_0_xy_yz_0_y = buffer_2000_ddsp[247];

    auto g_xz_0_0_0_xy_yz_0_z = buffer_2000_ddsp[248];

    auto g_xz_0_0_0_xy_zz_0_x = buffer_2000_ddsp[249];

    auto g_xz_0_0_0_xy_zz_0_y = buffer_2000_ddsp[250];

    auto g_xz_0_0_0_xy_zz_0_z = buffer_2000_ddsp[251];

    auto g_xz_0_0_0_xz_xx_0_x = buffer_2000_ddsp[252];

    auto g_xz_0_0_0_xz_xx_0_y = buffer_2000_ddsp[253];

    auto g_xz_0_0_0_xz_xx_0_z = buffer_2000_ddsp[254];

    auto g_xz_0_0_0_xz_xy_0_x = buffer_2000_ddsp[255];

    auto g_xz_0_0_0_xz_xy_0_y = buffer_2000_ddsp[256];

    auto g_xz_0_0_0_xz_xy_0_z = buffer_2000_ddsp[257];

    auto g_xz_0_0_0_xz_xz_0_x = buffer_2000_ddsp[258];

    auto g_xz_0_0_0_xz_xz_0_y = buffer_2000_ddsp[259];

    auto g_xz_0_0_0_xz_xz_0_z = buffer_2000_ddsp[260];

    auto g_xz_0_0_0_xz_yy_0_x = buffer_2000_ddsp[261];

    auto g_xz_0_0_0_xz_yy_0_y = buffer_2000_ddsp[262];

    auto g_xz_0_0_0_xz_yy_0_z = buffer_2000_ddsp[263];

    auto g_xz_0_0_0_xz_yz_0_x = buffer_2000_ddsp[264];

    auto g_xz_0_0_0_xz_yz_0_y = buffer_2000_ddsp[265];

    auto g_xz_0_0_0_xz_yz_0_z = buffer_2000_ddsp[266];

    auto g_xz_0_0_0_xz_zz_0_x = buffer_2000_ddsp[267];

    auto g_xz_0_0_0_xz_zz_0_y = buffer_2000_ddsp[268];

    auto g_xz_0_0_0_xz_zz_0_z = buffer_2000_ddsp[269];

    auto g_xz_0_0_0_yy_xx_0_x = buffer_2000_ddsp[270];

    auto g_xz_0_0_0_yy_xx_0_y = buffer_2000_ddsp[271];

    auto g_xz_0_0_0_yy_xx_0_z = buffer_2000_ddsp[272];

    auto g_xz_0_0_0_yy_xy_0_x = buffer_2000_ddsp[273];

    auto g_xz_0_0_0_yy_xy_0_y = buffer_2000_ddsp[274];

    auto g_xz_0_0_0_yy_xy_0_z = buffer_2000_ddsp[275];

    auto g_xz_0_0_0_yy_xz_0_x = buffer_2000_ddsp[276];

    auto g_xz_0_0_0_yy_xz_0_y = buffer_2000_ddsp[277];

    auto g_xz_0_0_0_yy_xz_0_z = buffer_2000_ddsp[278];

    auto g_xz_0_0_0_yy_yy_0_x = buffer_2000_ddsp[279];

    auto g_xz_0_0_0_yy_yy_0_y = buffer_2000_ddsp[280];

    auto g_xz_0_0_0_yy_yy_0_z = buffer_2000_ddsp[281];

    auto g_xz_0_0_0_yy_yz_0_x = buffer_2000_ddsp[282];

    auto g_xz_0_0_0_yy_yz_0_y = buffer_2000_ddsp[283];

    auto g_xz_0_0_0_yy_yz_0_z = buffer_2000_ddsp[284];

    auto g_xz_0_0_0_yy_zz_0_x = buffer_2000_ddsp[285];

    auto g_xz_0_0_0_yy_zz_0_y = buffer_2000_ddsp[286];

    auto g_xz_0_0_0_yy_zz_0_z = buffer_2000_ddsp[287];

    auto g_xz_0_0_0_yz_xx_0_x = buffer_2000_ddsp[288];

    auto g_xz_0_0_0_yz_xx_0_y = buffer_2000_ddsp[289];

    auto g_xz_0_0_0_yz_xx_0_z = buffer_2000_ddsp[290];

    auto g_xz_0_0_0_yz_xy_0_x = buffer_2000_ddsp[291];

    auto g_xz_0_0_0_yz_xy_0_y = buffer_2000_ddsp[292];

    auto g_xz_0_0_0_yz_xy_0_z = buffer_2000_ddsp[293];

    auto g_xz_0_0_0_yz_xz_0_x = buffer_2000_ddsp[294];

    auto g_xz_0_0_0_yz_xz_0_y = buffer_2000_ddsp[295];

    auto g_xz_0_0_0_yz_xz_0_z = buffer_2000_ddsp[296];

    auto g_xz_0_0_0_yz_yy_0_x = buffer_2000_ddsp[297];

    auto g_xz_0_0_0_yz_yy_0_y = buffer_2000_ddsp[298];

    auto g_xz_0_0_0_yz_yy_0_z = buffer_2000_ddsp[299];

    auto g_xz_0_0_0_yz_yz_0_x = buffer_2000_ddsp[300];

    auto g_xz_0_0_0_yz_yz_0_y = buffer_2000_ddsp[301];

    auto g_xz_0_0_0_yz_yz_0_z = buffer_2000_ddsp[302];

    auto g_xz_0_0_0_yz_zz_0_x = buffer_2000_ddsp[303];

    auto g_xz_0_0_0_yz_zz_0_y = buffer_2000_ddsp[304];

    auto g_xz_0_0_0_yz_zz_0_z = buffer_2000_ddsp[305];

    auto g_xz_0_0_0_zz_xx_0_x = buffer_2000_ddsp[306];

    auto g_xz_0_0_0_zz_xx_0_y = buffer_2000_ddsp[307];

    auto g_xz_0_0_0_zz_xx_0_z = buffer_2000_ddsp[308];

    auto g_xz_0_0_0_zz_xy_0_x = buffer_2000_ddsp[309];

    auto g_xz_0_0_0_zz_xy_0_y = buffer_2000_ddsp[310];

    auto g_xz_0_0_0_zz_xy_0_z = buffer_2000_ddsp[311];

    auto g_xz_0_0_0_zz_xz_0_x = buffer_2000_ddsp[312];

    auto g_xz_0_0_0_zz_xz_0_y = buffer_2000_ddsp[313];

    auto g_xz_0_0_0_zz_xz_0_z = buffer_2000_ddsp[314];

    auto g_xz_0_0_0_zz_yy_0_x = buffer_2000_ddsp[315];

    auto g_xz_0_0_0_zz_yy_0_y = buffer_2000_ddsp[316];

    auto g_xz_0_0_0_zz_yy_0_z = buffer_2000_ddsp[317];

    auto g_xz_0_0_0_zz_yz_0_x = buffer_2000_ddsp[318];

    auto g_xz_0_0_0_zz_yz_0_y = buffer_2000_ddsp[319];

    auto g_xz_0_0_0_zz_yz_0_z = buffer_2000_ddsp[320];

    auto g_xz_0_0_0_zz_zz_0_x = buffer_2000_ddsp[321];

    auto g_xz_0_0_0_zz_zz_0_y = buffer_2000_ddsp[322];

    auto g_xz_0_0_0_zz_zz_0_z = buffer_2000_ddsp[323];

    auto g_yy_0_0_0_xx_xx_0_x = buffer_2000_ddsp[324];

    auto g_yy_0_0_0_xx_xx_0_y = buffer_2000_ddsp[325];

    auto g_yy_0_0_0_xx_xx_0_z = buffer_2000_ddsp[326];

    auto g_yy_0_0_0_xx_xy_0_x = buffer_2000_ddsp[327];

    auto g_yy_0_0_0_xx_xy_0_y = buffer_2000_ddsp[328];

    auto g_yy_0_0_0_xx_xy_0_z = buffer_2000_ddsp[329];

    auto g_yy_0_0_0_xx_xz_0_x = buffer_2000_ddsp[330];

    auto g_yy_0_0_0_xx_xz_0_y = buffer_2000_ddsp[331];

    auto g_yy_0_0_0_xx_xz_0_z = buffer_2000_ddsp[332];

    auto g_yy_0_0_0_xx_yy_0_x = buffer_2000_ddsp[333];

    auto g_yy_0_0_0_xx_yy_0_y = buffer_2000_ddsp[334];

    auto g_yy_0_0_0_xx_yy_0_z = buffer_2000_ddsp[335];

    auto g_yy_0_0_0_xx_yz_0_x = buffer_2000_ddsp[336];

    auto g_yy_0_0_0_xx_yz_0_y = buffer_2000_ddsp[337];

    auto g_yy_0_0_0_xx_yz_0_z = buffer_2000_ddsp[338];

    auto g_yy_0_0_0_xx_zz_0_x = buffer_2000_ddsp[339];

    auto g_yy_0_0_0_xx_zz_0_y = buffer_2000_ddsp[340];

    auto g_yy_0_0_0_xx_zz_0_z = buffer_2000_ddsp[341];

    auto g_yy_0_0_0_xy_xx_0_x = buffer_2000_ddsp[342];

    auto g_yy_0_0_0_xy_xx_0_y = buffer_2000_ddsp[343];

    auto g_yy_0_0_0_xy_xx_0_z = buffer_2000_ddsp[344];

    auto g_yy_0_0_0_xy_xy_0_x = buffer_2000_ddsp[345];

    auto g_yy_0_0_0_xy_xy_0_y = buffer_2000_ddsp[346];

    auto g_yy_0_0_0_xy_xy_0_z = buffer_2000_ddsp[347];

    auto g_yy_0_0_0_xy_xz_0_x = buffer_2000_ddsp[348];

    auto g_yy_0_0_0_xy_xz_0_y = buffer_2000_ddsp[349];

    auto g_yy_0_0_0_xy_xz_0_z = buffer_2000_ddsp[350];

    auto g_yy_0_0_0_xy_yy_0_x = buffer_2000_ddsp[351];

    auto g_yy_0_0_0_xy_yy_0_y = buffer_2000_ddsp[352];

    auto g_yy_0_0_0_xy_yy_0_z = buffer_2000_ddsp[353];

    auto g_yy_0_0_0_xy_yz_0_x = buffer_2000_ddsp[354];

    auto g_yy_0_0_0_xy_yz_0_y = buffer_2000_ddsp[355];

    auto g_yy_0_0_0_xy_yz_0_z = buffer_2000_ddsp[356];

    auto g_yy_0_0_0_xy_zz_0_x = buffer_2000_ddsp[357];

    auto g_yy_0_0_0_xy_zz_0_y = buffer_2000_ddsp[358];

    auto g_yy_0_0_0_xy_zz_0_z = buffer_2000_ddsp[359];

    auto g_yy_0_0_0_xz_xx_0_x = buffer_2000_ddsp[360];

    auto g_yy_0_0_0_xz_xx_0_y = buffer_2000_ddsp[361];

    auto g_yy_0_0_0_xz_xx_0_z = buffer_2000_ddsp[362];

    auto g_yy_0_0_0_xz_xy_0_x = buffer_2000_ddsp[363];

    auto g_yy_0_0_0_xz_xy_0_y = buffer_2000_ddsp[364];

    auto g_yy_0_0_0_xz_xy_0_z = buffer_2000_ddsp[365];

    auto g_yy_0_0_0_xz_xz_0_x = buffer_2000_ddsp[366];

    auto g_yy_0_0_0_xz_xz_0_y = buffer_2000_ddsp[367];

    auto g_yy_0_0_0_xz_xz_0_z = buffer_2000_ddsp[368];

    auto g_yy_0_0_0_xz_yy_0_x = buffer_2000_ddsp[369];

    auto g_yy_0_0_0_xz_yy_0_y = buffer_2000_ddsp[370];

    auto g_yy_0_0_0_xz_yy_0_z = buffer_2000_ddsp[371];

    auto g_yy_0_0_0_xz_yz_0_x = buffer_2000_ddsp[372];

    auto g_yy_0_0_0_xz_yz_0_y = buffer_2000_ddsp[373];

    auto g_yy_0_0_0_xz_yz_0_z = buffer_2000_ddsp[374];

    auto g_yy_0_0_0_xz_zz_0_x = buffer_2000_ddsp[375];

    auto g_yy_0_0_0_xz_zz_0_y = buffer_2000_ddsp[376];

    auto g_yy_0_0_0_xz_zz_0_z = buffer_2000_ddsp[377];

    auto g_yy_0_0_0_yy_xx_0_x = buffer_2000_ddsp[378];

    auto g_yy_0_0_0_yy_xx_0_y = buffer_2000_ddsp[379];

    auto g_yy_0_0_0_yy_xx_0_z = buffer_2000_ddsp[380];

    auto g_yy_0_0_0_yy_xy_0_x = buffer_2000_ddsp[381];

    auto g_yy_0_0_0_yy_xy_0_y = buffer_2000_ddsp[382];

    auto g_yy_0_0_0_yy_xy_0_z = buffer_2000_ddsp[383];

    auto g_yy_0_0_0_yy_xz_0_x = buffer_2000_ddsp[384];

    auto g_yy_0_0_0_yy_xz_0_y = buffer_2000_ddsp[385];

    auto g_yy_0_0_0_yy_xz_0_z = buffer_2000_ddsp[386];

    auto g_yy_0_0_0_yy_yy_0_x = buffer_2000_ddsp[387];

    auto g_yy_0_0_0_yy_yy_0_y = buffer_2000_ddsp[388];

    auto g_yy_0_0_0_yy_yy_0_z = buffer_2000_ddsp[389];

    auto g_yy_0_0_0_yy_yz_0_x = buffer_2000_ddsp[390];

    auto g_yy_0_0_0_yy_yz_0_y = buffer_2000_ddsp[391];

    auto g_yy_0_0_0_yy_yz_0_z = buffer_2000_ddsp[392];

    auto g_yy_0_0_0_yy_zz_0_x = buffer_2000_ddsp[393];

    auto g_yy_0_0_0_yy_zz_0_y = buffer_2000_ddsp[394];

    auto g_yy_0_0_0_yy_zz_0_z = buffer_2000_ddsp[395];

    auto g_yy_0_0_0_yz_xx_0_x = buffer_2000_ddsp[396];

    auto g_yy_0_0_0_yz_xx_0_y = buffer_2000_ddsp[397];

    auto g_yy_0_0_0_yz_xx_0_z = buffer_2000_ddsp[398];

    auto g_yy_0_0_0_yz_xy_0_x = buffer_2000_ddsp[399];

    auto g_yy_0_0_0_yz_xy_0_y = buffer_2000_ddsp[400];

    auto g_yy_0_0_0_yz_xy_0_z = buffer_2000_ddsp[401];

    auto g_yy_0_0_0_yz_xz_0_x = buffer_2000_ddsp[402];

    auto g_yy_0_0_0_yz_xz_0_y = buffer_2000_ddsp[403];

    auto g_yy_0_0_0_yz_xz_0_z = buffer_2000_ddsp[404];

    auto g_yy_0_0_0_yz_yy_0_x = buffer_2000_ddsp[405];

    auto g_yy_0_0_0_yz_yy_0_y = buffer_2000_ddsp[406];

    auto g_yy_0_0_0_yz_yy_0_z = buffer_2000_ddsp[407];

    auto g_yy_0_0_0_yz_yz_0_x = buffer_2000_ddsp[408];

    auto g_yy_0_0_0_yz_yz_0_y = buffer_2000_ddsp[409];

    auto g_yy_0_0_0_yz_yz_0_z = buffer_2000_ddsp[410];

    auto g_yy_0_0_0_yz_zz_0_x = buffer_2000_ddsp[411];

    auto g_yy_0_0_0_yz_zz_0_y = buffer_2000_ddsp[412];

    auto g_yy_0_0_0_yz_zz_0_z = buffer_2000_ddsp[413];

    auto g_yy_0_0_0_zz_xx_0_x = buffer_2000_ddsp[414];

    auto g_yy_0_0_0_zz_xx_0_y = buffer_2000_ddsp[415];

    auto g_yy_0_0_0_zz_xx_0_z = buffer_2000_ddsp[416];

    auto g_yy_0_0_0_zz_xy_0_x = buffer_2000_ddsp[417];

    auto g_yy_0_0_0_zz_xy_0_y = buffer_2000_ddsp[418];

    auto g_yy_0_0_0_zz_xy_0_z = buffer_2000_ddsp[419];

    auto g_yy_0_0_0_zz_xz_0_x = buffer_2000_ddsp[420];

    auto g_yy_0_0_0_zz_xz_0_y = buffer_2000_ddsp[421];

    auto g_yy_0_0_0_zz_xz_0_z = buffer_2000_ddsp[422];

    auto g_yy_0_0_0_zz_yy_0_x = buffer_2000_ddsp[423];

    auto g_yy_0_0_0_zz_yy_0_y = buffer_2000_ddsp[424];

    auto g_yy_0_0_0_zz_yy_0_z = buffer_2000_ddsp[425];

    auto g_yy_0_0_0_zz_yz_0_x = buffer_2000_ddsp[426];

    auto g_yy_0_0_0_zz_yz_0_y = buffer_2000_ddsp[427];

    auto g_yy_0_0_0_zz_yz_0_z = buffer_2000_ddsp[428];

    auto g_yy_0_0_0_zz_zz_0_x = buffer_2000_ddsp[429];

    auto g_yy_0_0_0_zz_zz_0_y = buffer_2000_ddsp[430];

    auto g_yy_0_0_0_zz_zz_0_z = buffer_2000_ddsp[431];

    auto g_yz_0_0_0_xx_xx_0_x = buffer_2000_ddsp[432];

    auto g_yz_0_0_0_xx_xx_0_y = buffer_2000_ddsp[433];

    auto g_yz_0_0_0_xx_xx_0_z = buffer_2000_ddsp[434];

    auto g_yz_0_0_0_xx_xy_0_x = buffer_2000_ddsp[435];

    auto g_yz_0_0_0_xx_xy_0_y = buffer_2000_ddsp[436];

    auto g_yz_0_0_0_xx_xy_0_z = buffer_2000_ddsp[437];

    auto g_yz_0_0_0_xx_xz_0_x = buffer_2000_ddsp[438];

    auto g_yz_0_0_0_xx_xz_0_y = buffer_2000_ddsp[439];

    auto g_yz_0_0_0_xx_xz_0_z = buffer_2000_ddsp[440];

    auto g_yz_0_0_0_xx_yy_0_x = buffer_2000_ddsp[441];

    auto g_yz_0_0_0_xx_yy_0_y = buffer_2000_ddsp[442];

    auto g_yz_0_0_0_xx_yy_0_z = buffer_2000_ddsp[443];

    auto g_yz_0_0_0_xx_yz_0_x = buffer_2000_ddsp[444];

    auto g_yz_0_0_0_xx_yz_0_y = buffer_2000_ddsp[445];

    auto g_yz_0_0_0_xx_yz_0_z = buffer_2000_ddsp[446];

    auto g_yz_0_0_0_xx_zz_0_x = buffer_2000_ddsp[447];

    auto g_yz_0_0_0_xx_zz_0_y = buffer_2000_ddsp[448];

    auto g_yz_0_0_0_xx_zz_0_z = buffer_2000_ddsp[449];

    auto g_yz_0_0_0_xy_xx_0_x = buffer_2000_ddsp[450];

    auto g_yz_0_0_0_xy_xx_0_y = buffer_2000_ddsp[451];

    auto g_yz_0_0_0_xy_xx_0_z = buffer_2000_ddsp[452];

    auto g_yz_0_0_0_xy_xy_0_x = buffer_2000_ddsp[453];

    auto g_yz_0_0_0_xy_xy_0_y = buffer_2000_ddsp[454];

    auto g_yz_0_0_0_xy_xy_0_z = buffer_2000_ddsp[455];

    auto g_yz_0_0_0_xy_xz_0_x = buffer_2000_ddsp[456];

    auto g_yz_0_0_0_xy_xz_0_y = buffer_2000_ddsp[457];

    auto g_yz_0_0_0_xy_xz_0_z = buffer_2000_ddsp[458];

    auto g_yz_0_0_0_xy_yy_0_x = buffer_2000_ddsp[459];

    auto g_yz_0_0_0_xy_yy_0_y = buffer_2000_ddsp[460];

    auto g_yz_0_0_0_xy_yy_0_z = buffer_2000_ddsp[461];

    auto g_yz_0_0_0_xy_yz_0_x = buffer_2000_ddsp[462];

    auto g_yz_0_0_0_xy_yz_0_y = buffer_2000_ddsp[463];

    auto g_yz_0_0_0_xy_yz_0_z = buffer_2000_ddsp[464];

    auto g_yz_0_0_0_xy_zz_0_x = buffer_2000_ddsp[465];

    auto g_yz_0_0_0_xy_zz_0_y = buffer_2000_ddsp[466];

    auto g_yz_0_0_0_xy_zz_0_z = buffer_2000_ddsp[467];

    auto g_yz_0_0_0_xz_xx_0_x = buffer_2000_ddsp[468];

    auto g_yz_0_0_0_xz_xx_0_y = buffer_2000_ddsp[469];

    auto g_yz_0_0_0_xz_xx_0_z = buffer_2000_ddsp[470];

    auto g_yz_0_0_0_xz_xy_0_x = buffer_2000_ddsp[471];

    auto g_yz_0_0_0_xz_xy_0_y = buffer_2000_ddsp[472];

    auto g_yz_0_0_0_xz_xy_0_z = buffer_2000_ddsp[473];

    auto g_yz_0_0_0_xz_xz_0_x = buffer_2000_ddsp[474];

    auto g_yz_0_0_0_xz_xz_0_y = buffer_2000_ddsp[475];

    auto g_yz_0_0_0_xz_xz_0_z = buffer_2000_ddsp[476];

    auto g_yz_0_0_0_xz_yy_0_x = buffer_2000_ddsp[477];

    auto g_yz_0_0_0_xz_yy_0_y = buffer_2000_ddsp[478];

    auto g_yz_0_0_0_xz_yy_0_z = buffer_2000_ddsp[479];

    auto g_yz_0_0_0_xz_yz_0_x = buffer_2000_ddsp[480];

    auto g_yz_0_0_0_xz_yz_0_y = buffer_2000_ddsp[481];

    auto g_yz_0_0_0_xz_yz_0_z = buffer_2000_ddsp[482];

    auto g_yz_0_0_0_xz_zz_0_x = buffer_2000_ddsp[483];

    auto g_yz_0_0_0_xz_zz_0_y = buffer_2000_ddsp[484];

    auto g_yz_0_0_0_xz_zz_0_z = buffer_2000_ddsp[485];

    auto g_yz_0_0_0_yy_xx_0_x = buffer_2000_ddsp[486];

    auto g_yz_0_0_0_yy_xx_0_y = buffer_2000_ddsp[487];

    auto g_yz_0_0_0_yy_xx_0_z = buffer_2000_ddsp[488];

    auto g_yz_0_0_0_yy_xy_0_x = buffer_2000_ddsp[489];

    auto g_yz_0_0_0_yy_xy_0_y = buffer_2000_ddsp[490];

    auto g_yz_0_0_0_yy_xy_0_z = buffer_2000_ddsp[491];

    auto g_yz_0_0_0_yy_xz_0_x = buffer_2000_ddsp[492];

    auto g_yz_0_0_0_yy_xz_0_y = buffer_2000_ddsp[493];

    auto g_yz_0_0_0_yy_xz_0_z = buffer_2000_ddsp[494];

    auto g_yz_0_0_0_yy_yy_0_x = buffer_2000_ddsp[495];

    auto g_yz_0_0_0_yy_yy_0_y = buffer_2000_ddsp[496];

    auto g_yz_0_0_0_yy_yy_0_z = buffer_2000_ddsp[497];

    auto g_yz_0_0_0_yy_yz_0_x = buffer_2000_ddsp[498];

    auto g_yz_0_0_0_yy_yz_0_y = buffer_2000_ddsp[499];

    auto g_yz_0_0_0_yy_yz_0_z = buffer_2000_ddsp[500];

    auto g_yz_0_0_0_yy_zz_0_x = buffer_2000_ddsp[501];

    auto g_yz_0_0_0_yy_zz_0_y = buffer_2000_ddsp[502];

    auto g_yz_0_0_0_yy_zz_0_z = buffer_2000_ddsp[503];

    auto g_yz_0_0_0_yz_xx_0_x = buffer_2000_ddsp[504];

    auto g_yz_0_0_0_yz_xx_0_y = buffer_2000_ddsp[505];

    auto g_yz_0_0_0_yz_xx_0_z = buffer_2000_ddsp[506];

    auto g_yz_0_0_0_yz_xy_0_x = buffer_2000_ddsp[507];

    auto g_yz_0_0_0_yz_xy_0_y = buffer_2000_ddsp[508];

    auto g_yz_0_0_0_yz_xy_0_z = buffer_2000_ddsp[509];

    auto g_yz_0_0_0_yz_xz_0_x = buffer_2000_ddsp[510];

    auto g_yz_0_0_0_yz_xz_0_y = buffer_2000_ddsp[511];

    auto g_yz_0_0_0_yz_xz_0_z = buffer_2000_ddsp[512];

    auto g_yz_0_0_0_yz_yy_0_x = buffer_2000_ddsp[513];

    auto g_yz_0_0_0_yz_yy_0_y = buffer_2000_ddsp[514];

    auto g_yz_0_0_0_yz_yy_0_z = buffer_2000_ddsp[515];

    auto g_yz_0_0_0_yz_yz_0_x = buffer_2000_ddsp[516];

    auto g_yz_0_0_0_yz_yz_0_y = buffer_2000_ddsp[517];

    auto g_yz_0_0_0_yz_yz_0_z = buffer_2000_ddsp[518];

    auto g_yz_0_0_0_yz_zz_0_x = buffer_2000_ddsp[519];

    auto g_yz_0_0_0_yz_zz_0_y = buffer_2000_ddsp[520];

    auto g_yz_0_0_0_yz_zz_0_z = buffer_2000_ddsp[521];

    auto g_yz_0_0_0_zz_xx_0_x = buffer_2000_ddsp[522];

    auto g_yz_0_0_0_zz_xx_0_y = buffer_2000_ddsp[523];

    auto g_yz_0_0_0_zz_xx_0_z = buffer_2000_ddsp[524];

    auto g_yz_0_0_0_zz_xy_0_x = buffer_2000_ddsp[525];

    auto g_yz_0_0_0_zz_xy_0_y = buffer_2000_ddsp[526];

    auto g_yz_0_0_0_zz_xy_0_z = buffer_2000_ddsp[527];

    auto g_yz_0_0_0_zz_xz_0_x = buffer_2000_ddsp[528];

    auto g_yz_0_0_0_zz_xz_0_y = buffer_2000_ddsp[529];

    auto g_yz_0_0_0_zz_xz_0_z = buffer_2000_ddsp[530];

    auto g_yz_0_0_0_zz_yy_0_x = buffer_2000_ddsp[531];

    auto g_yz_0_0_0_zz_yy_0_y = buffer_2000_ddsp[532];

    auto g_yz_0_0_0_zz_yy_0_z = buffer_2000_ddsp[533];

    auto g_yz_0_0_0_zz_yz_0_x = buffer_2000_ddsp[534];

    auto g_yz_0_0_0_zz_yz_0_y = buffer_2000_ddsp[535];

    auto g_yz_0_0_0_zz_yz_0_z = buffer_2000_ddsp[536];

    auto g_yz_0_0_0_zz_zz_0_x = buffer_2000_ddsp[537];

    auto g_yz_0_0_0_zz_zz_0_y = buffer_2000_ddsp[538];

    auto g_yz_0_0_0_zz_zz_0_z = buffer_2000_ddsp[539];

    auto g_zz_0_0_0_xx_xx_0_x = buffer_2000_ddsp[540];

    auto g_zz_0_0_0_xx_xx_0_y = buffer_2000_ddsp[541];

    auto g_zz_0_0_0_xx_xx_0_z = buffer_2000_ddsp[542];

    auto g_zz_0_0_0_xx_xy_0_x = buffer_2000_ddsp[543];

    auto g_zz_0_0_0_xx_xy_0_y = buffer_2000_ddsp[544];

    auto g_zz_0_0_0_xx_xy_0_z = buffer_2000_ddsp[545];

    auto g_zz_0_0_0_xx_xz_0_x = buffer_2000_ddsp[546];

    auto g_zz_0_0_0_xx_xz_0_y = buffer_2000_ddsp[547];

    auto g_zz_0_0_0_xx_xz_0_z = buffer_2000_ddsp[548];

    auto g_zz_0_0_0_xx_yy_0_x = buffer_2000_ddsp[549];

    auto g_zz_0_0_0_xx_yy_0_y = buffer_2000_ddsp[550];

    auto g_zz_0_0_0_xx_yy_0_z = buffer_2000_ddsp[551];

    auto g_zz_0_0_0_xx_yz_0_x = buffer_2000_ddsp[552];

    auto g_zz_0_0_0_xx_yz_0_y = buffer_2000_ddsp[553];

    auto g_zz_0_0_0_xx_yz_0_z = buffer_2000_ddsp[554];

    auto g_zz_0_0_0_xx_zz_0_x = buffer_2000_ddsp[555];

    auto g_zz_0_0_0_xx_zz_0_y = buffer_2000_ddsp[556];

    auto g_zz_0_0_0_xx_zz_0_z = buffer_2000_ddsp[557];

    auto g_zz_0_0_0_xy_xx_0_x = buffer_2000_ddsp[558];

    auto g_zz_0_0_0_xy_xx_0_y = buffer_2000_ddsp[559];

    auto g_zz_0_0_0_xy_xx_0_z = buffer_2000_ddsp[560];

    auto g_zz_0_0_0_xy_xy_0_x = buffer_2000_ddsp[561];

    auto g_zz_0_0_0_xy_xy_0_y = buffer_2000_ddsp[562];

    auto g_zz_0_0_0_xy_xy_0_z = buffer_2000_ddsp[563];

    auto g_zz_0_0_0_xy_xz_0_x = buffer_2000_ddsp[564];

    auto g_zz_0_0_0_xy_xz_0_y = buffer_2000_ddsp[565];

    auto g_zz_0_0_0_xy_xz_0_z = buffer_2000_ddsp[566];

    auto g_zz_0_0_0_xy_yy_0_x = buffer_2000_ddsp[567];

    auto g_zz_0_0_0_xy_yy_0_y = buffer_2000_ddsp[568];

    auto g_zz_0_0_0_xy_yy_0_z = buffer_2000_ddsp[569];

    auto g_zz_0_0_0_xy_yz_0_x = buffer_2000_ddsp[570];

    auto g_zz_0_0_0_xy_yz_0_y = buffer_2000_ddsp[571];

    auto g_zz_0_0_0_xy_yz_0_z = buffer_2000_ddsp[572];

    auto g_zz_0_0_0_xy_zz_0_x = buffer_2000_ddsp[573];

    auto g_zz_0_0_0_xy_zz_0_y = buffer_2000_ddsp[574];

    auto g_zz_0_0_0_xy_zz_0_z = buffer_2000_ddsp[575];

    auto g_zz_0_0_0_xz_xx_0_x = buffer_2000_ddsp[576];

    auto g_zz_0_0_0_xz_xx_0_y = buffer_2000_ddsp[577];

    auto g_zz_0_0_0_xz_xx_0_z = buffer_2000_ddsp[578];

    auto g_zz_0_0_0_xz_xy_0_x = buffer_2000_ddsp[579];

    auto g_zz_0_0_0_xz_xy_0_y = buffer_2000_ddsp[580];

    auto g_zz_0_0_0_xz_xy_0_z = buffer_2000_ddsp[581];

    auto g_zz_0_0_0_xz_xz_0_x = buffer_2000_ddsp[582];

    auto g_zz_0_0_0_xz_xz_0_y = buffer_2000_ddsp[583];

    auto g_zz_0_0_0_xz_xz_0_z = buffer_2000_ddsp[584];

    auto g_zz_0_0_0_xz_yy_0_x = buffer_2000_ddsp[585];

    auto g_zz_0_0_0_xz_yy_0_y = buffer_2000_ddsp[586];

    auto g_zz_0_0_0_xz_yy_0_z = buffer_2000_ddsp[587];

    auto g_zz_0_0_0_xz_yz_0_x = buffer_2000_ddsp[588];

    auto g_zz_0_0_0_xz_yz_0_y = buffer_2000_ddsp[589];

    auto g_zz_0_0_0_xz_yz_0_z = buffer_2000_ddsp[590];

    auto g_zz_0_0_0_xz_zz_0_x = buffer_2000_ddsp[591];

    auto g_zz_0_0_0_xz_zz_0_y = buffer_2000_ddsp[592];

    auto g_zz_0_0_0_xz_zz_0_z = buffer_2000_ddsp[593];

    auto g_zz_0_0_0_yy_xx_0_x = buffer_2000_ddsp[594];

    auto g_zz_0_0_0_yy_xx_0_y = buffer_2000_ddsp[595];

    auto g_zz_0_0_0_yy_xx_0_z = buffer_2000_ddsp[596];

    auto g_zz_0_0_0_yy_xy_0_x = buffer_2000_ddsp[597];

    auto g_zz_0_0_0_yy_xy_0_y = buffer_2000_ddsp[598];

    auto g_zz_0_0_0_yy_xy_0_z = buffer_2000_ddsp[599];

    auto g_zz_0_0_0_yy_xz_0_x = buffer_2000_ddsp[600];

    auto g_zz_0_0_0_yy_xz_0_y = buffer_2000_ddsp[601];

    auto g_zz_0_0_0_yy_xz_0_z = buffer_2000_ddsp[602];

    auto g_zz_0_0_0_yy_yy_0_x = buffer_2000_ddsp[603];

    auto g_zz_0_0_0_yy_yy_0_y = buffer_2000_ddsp[604];

    auto g_zz_0_0_0_yy_yy_0_z = buffer_2000_ddsp[605];

    auto g_zz_0_0_0_yy_yz_0_x = buffer_2000_ddsp[606];

    auto g_zz_0_0_0_yy_yz_0_y = buffer_2000_ddsp[607];

    auto g_zz_0_0_0_yy_yz_0_z = buffer_2000_ddsp[608];

    auto g_zz_0_0_0_yy_zz_0_x = buffer_2000_ddsp[609];

    auto g_zz_0_0_0_yy_zz_0_y = buffer_2000_ddsp[610];

    auto g_zz_0_0_0_yy_zz_0_z = buffer_2000_ddsp[611];

    auto g_zz_0_0_0_yz_xx_0_x = buffer_2000_ddsp[612];

    auto g_zz_0_0_0_yz_xx_0_y = buffer_2000_ddsp[613];

    auto g_zz_0_0_0_yz_xx_0_z = buffer_2000_ddsp[614];

    auto g_zz_0_0_0_yz_xy_0_x = buffer_2000_ddsp[615];

    auto g_zz_0_0_0_yz_xy_0_y = buffer_2000_ddsp[616];

    auto g_zz_0_0_0_yz_xy_0_z = buffer_2000_ddsp[617];

    auto g_zz_0_0_0_yz_xz_0_x = buffer_2000_ddsp[618];

    auto g_zz_0_0_0_yz_xz_0_y = buffer_2000_ddsp[619];

    auto g_zz_0_0_0_yz_xz_0_z = buffer_2000_ddsp[620];

    auto g_zz_0_0_0_yz_yy_0_x = buffer_2000_ddsp[621];

    auto g_zz_0_0_0_yz_yy_0_y = buffer_2000_ddsp[622];

    auto g_zz_0_0_0_yz_yy_0_z = buffer_2000_ddsp[623];

    auto g_zz_0_0_0_yz_yz_0_x = buffer_2000_ddsp[624];

    auto g_zz_0_0_0_yz_yz_0_y = buffer_2000_ddsp[625];

    auto g_zz_0_0_0_yz_yz_0_z = buffer_2000_ddsp[626];

    auto g_zz_0_0_0_yz_zz_0_x = buffer_2000_ddsp[627];

    auto g_zz_0_0_0_yz_zz_0_y = buffer_2000_ddsp[628];

    auto g_zz_0_0_0_yz_zz_0_z = buffer_2000_ddsp[629];

    auto g_zz_0_0_0_zz_xx_0_x = buffer_2000_ddsp[630];

    auto g_zz_0_0_0_zz_xx_0_y = buffer_2000_ddsp[631];

    auto g_zz_0_0_0_zz_xx_0_z = buffer_2000_ddsp[632];

    auto g_zz_0_0_0_zz_xy_0_x = buffer_2000_ddsp[633];

    auto g_zz_0_0_0_zz_xy_0_y = buffer_2000_ddsp[634];

    auto g_zz_0_0_0_zz_xy_0_z = buffer_2000_ddsp[635];

    auto g_zz_0_0_0_zz_xz_0_x = buffer_2000_ddsp[636];

    auto g_zz_0_0_0_zz_xz_0_y = buffer_2000_ddsp[637];

    auto g_zz_0_0_0_zz_xz_0_z = buffer_2000_ddsp[638];

    auto g_zz_0_0_0_zz_yy_0_x = buffer_2000_ddsp[639];

    auto g_zz_0_0_0_zz_yy_0_y = buffer_2000_ddsp[640];

    auto g_zz_0_0_0_zz_yy_0_z = buffer_2000_ddsp[641];

    auto g_zz_0_0_0_zz_yz_0_x = buffer_2000_ddsp[642];

    auto g_zz_0_0_0_zz_yz_0_y = buffer_2000_ddsp[643];

    auto g_zz_0_0_0_zz_yz_0_z = buffer_2000_ddsp[644];

    auto g_zz_0_0_0_zz_zz_0_x = buffer_2000_ddsp[645];

    auto g_zz_0_0_0_zz_zz_0_y = buffer_2000_ddsp[646];

    auto g_zz_0_0_0_zz_zz_0_z = buffer_2000_ddsp[647];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_xx_0_0_0_xx_xx_0_x, g_xx_0_0_0_xx_xx_0_y, g_xx_0_0_0_xx_xx_0_z, g_xx_xx_0_x, g_xx_xx_0_y, g_xx_xx_0_z, g_xxxx_xx_0_x, g_xxxx_xx_0_y, g_xxxx_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_xx_0_x[i] = 2.0 * g_0_xx_0_x[i] - 10.0 * g_xx_xx_0_x[i] * a_exp + 4.0 * g_xxxx_xx_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xx_0_y[i] = 2.0 * g_0_xx_0_y[i] - 10.0 * g_xx_xx_0_y[i] * a_exp + 4.0 * g_xxxx_xx_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xx_0_z[i] = 2.0 * g_0_xx_0_z[i] - 10.0 * g_xx_xx_0_z[i] * a_exp + 4.0 * g_xxxx_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_xx_0_0_0_xx_xy_0_x, g_xx_0_0_0_xx_xy_0_y, g_xx_0_0_0_xx_xy_0_z, g_xx_xy_0_x, g_xx_xy_0_y, g_xx_xy_0_z, g_xxxx_xy_0_x, g_xxxx_xy_0_y, g_xxxx_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_xy_0_x[i] = 2.0 * g_0_xy_0_x[i] - 10.0 * g_xx_xy_0_x[i] * a_exp + 4.0 * g_xxxx_xy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xy_0_y[i] = 2.0 * g_0_xy_0_y[i] - 10.0 * g_xx_xy_0_y[i] * a_exp + 4.0 * g_xxxx_xy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xy_0_z[i] = 2.0 * g_0_xy_0_z[i] - 10.0 * g_xx_xy_0_z[i] * a_exp + 4.0 * g_xxxx_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_xx_0_0_0_xx_xz_0_x, g_xx_0_0_0_xx_xz_0_y, g_xx_0_0_0_xx_xz_0_z, g_xx_xz_0_x, g_xx_xz_0_y, g_xx_xz_0_z, g_xxxx_xz_0_x, g_xxxx_xz_0_y, g_xxxx_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_xz_0_x[i] = 2.0 * g_0_xz_0_x[i] - 10.0 * g_xx_xz_0_x[i] * a_exp + 4.0 * g_xxxx_xz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xz_0_y[i] = 2.0 * g_0_xz_0_y[i] - 10.0 * g_xx_xz_0_y[i] * a_exp + 4.0 * g_xxxx_xz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xz_0_z[i] = 2.0 * g_0_xz_0_z[i] - 10.0 * g_xx_xz_0_z[i] * a_exp + 4.0 * g_xxxx_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_xx_0_0_0_xx_yy_0_x, g_xx_0_0_0_xx_yy_0_y, g_xx_0_0_0_xx_yy_0_z, g_xx_yy_0_x, g_xx_yy_0_y, g_xx_yy_0_z, g_xxxx_yy_0_x, g_xxxx_yy_0_y, g_xxxx_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_yy_0_x[i] = 2.0 * g_0_yy_0_x[i] - 10.0 * g_xx_yy_0_x[i] * a_exp + 4.0 * g_xxxx_yy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yy_0_y[i] = 2.0 * g_0_yy_0_y[i] - 10.0 * g_xx_yy_0_y[i] * a_exp + 4.0 * g_xxxx_yy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yy_0_z[i] = 2.0 * g_0_yy_0_z[i] - 10.0 * g_xx_yy_0_z[i] * a_exp + 4.0 * g_xxxx_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_xx_0_0_0_xx_yz_0_x, g_xx_0_0_0_xx_yz_0_y, g_xx_0_0_0_xx_yz_0_z, g_xx_yz_0_x, g_xx_yz_0_y, g_xx_yz_0_z, g_xxxx_yz_0_x, g_xxxx_yz_0_y, g_xxxx_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_yz_0_x[i] = 2.0 * g_0_yz_0_x[i] - 10.0 * g_xx_yz_0_x[i] * a_exp + 4.0 * g_xxxx_yz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yz_0_y[i] = 2.0 * g_0_yz_0_y[i] - 10.0 * g_xx_yz_0_y[i] * a_exp + 4.0 * g_xxxx_yz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yz_0_z[i] = 2.0 * g_0_yz_0_z[i] - 10.0 * g_xx_yz_0_z[i] * a_exp + 4.0 * g_xxxx_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_xx_0_0_0_xx_zz_0_x, g_xx_0_0_0_xx_zz_0_y, g_xx_0_0_0_xx_zz_0_z, g_xx_zz_0_x, g_xx_zz_0_y, g_xx_zz_0_z, g_xxxx_zz_0_x, g_xxxx_zz_0_y, g_xxxx_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_zz_0_x[i] = 2.0 * g_0_zz_0_x[i] - 10.0 * g_xx_zz_0_x[i] * a_exp + 4.0 * g_xxxx_zz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_zz_0_y[i] = 2.0 * g_0_zz_0_y[i] - 10.0 * g_xx_zz_0_y[i] * a_exp + 4.0 * g_xxxx_zz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_zz_0_z[i] = 2.0 * g_0_zz_0_z[i] - 10.0 * g_xx_zz_0_z[i] * a_exp + 4.0 * g_xxxx_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_xx_0_0_0_xy_xx_0_x, g_xx_0_0_0_xy_xx_0_y, g_xx_0_0_0_xy_xx_0_z, g_xxxy_xx_0_x, g_xxxy_xx_0_y, g_xxxy_xx_0_z, g_xy_xx_0_x, g_xy_xx_0_y, g_xy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_xx_0_x[i] = -6.0 * g_xy_xx_0_x[i] * a_exp + 4.0 * g_xxxy_xx_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xx_0_y[i] = -6.0 * g_xy_xx_0_y[i] * a_exp + 4.0 * g_xxxy_xx_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xx_0_z[i] = -6.0 * g_xy_xx_0_z[i] * a_exp + 4.0 * g_xxxy_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_xx_0_0_0_xy_xy_0_x, g_xx_0_0_0_xy_xy_0_y, g_xx_0_0_0_xy_xy_0_z, g_xxxy_xy_0_x, g_xxxy_xy_0_y, g_xxxy_xy_0_z, g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_xy_0_x[i] = -6.0 * g_xy_xy_0_x[i] * a_exp + 4.0 * g_xxxy_xy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xy_0_y[i] = -6.0 * g_xy_xy_0_y[i] * a_exp + 4.0 * g_xxxy_xy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xy_0_z[i] = -6.0 * g_xy_xy_0_z[i] * a_exp + 4.0 * g_xxxy_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_xx_0_0_0_xy_xz_0_x, g_xx_0_0_0_xy_xz_0_y, g_xx_0_0_0_xy_xz_0_z, g_xxxy_xz_0_x, g_xxxy_xz_0_y, g_xxxy_xz_0_z, g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_xz_0_x[i] = -6.0 * g_xy_xz_0_x[i] * a_exp + 4.0 * g_xxxy_xz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xz_0_y[i] = -6.0 * g_xy_xz_0_y[i] * a_exp + 4.0 * g_xxxy_xz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xz_0_z[i] = -6.0 * g_xy_xz_0_z[i] * a_exp + 4.0 * g_xxxy_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_xx_0_0_0_xy_yy_0_x, g_xx_0_0_0_xy_yy_0_y, g_xx_0_0_0_xy_yy_0_z, g_xxxy_yy_0_x, g_xxxy_yy_0_y, g_xxxy_yy_0_z, g_xy_yy_0_x, g_xy_yy_0_y, g_xy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_yy_0_x[i] = -6.0 * g_xy_yy_0_x[i] * a_exp + 4.0 * g_xxxy_yy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yy_0_y[i] = -6.0 * g_xy_yy_0_y[i] * a_exp + 4.0 * g_xxxy_yy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yy_0_z[i] = -6.0 * g_xy_yy_0_z[i] * a_exp + 4.0 * g_xxxy_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_xx_0_0_0_xy_yz_0_x, g_xx_0_0_0_xy_yz_0_y, g_xx_0_0_0_xy_yz_0_z, g_xxxy_yz_0_x, g_xxxy_yz_0_y, g_xxxy_yz_0_z, g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_yz_0_x[i] = -6.0 * g_xy_yz_0_x[i] * a_exp + 4.0 * g_xxxy_yz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yz_0_y[i] = -6.0 * g_xy_yz_0_y[i] * a_exp + 4.0 * g_xxxy_yz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yz_0_z[i] = -6.0 * g_xy_yz_0_z[i] * a_exp + 4.0 * g_xxxy_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_xx_0_0_0_xy_zz_0_x, g_xx_0_0_0_xy_zz_0_y, g_xx_0_0_0_xy_zz_0_z, g_xxxy_zz_0_x, g_xxxy_zz_0_y, g_xxxy_zz_0_z, g_xy_zz_0_x, g_xy_zz_0_y, g_xy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_zz_0_x[i] = -6.0 * g_xy_zz_0_x[i] * a_exp + 4.0 * g_xxxy_zz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_zz_0_y[i] = -6.0 * g_xy_zz_0_y[i] * a_exp + 4.0 * g_xxxy_zz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_zz_0_z[i] = -6.0 * g_xy_zz_0_z[i] * a_exp + 4.0 * g_xxxy_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_xx_0_0_0_xz_xx_0_x, g_xx_0_0_0_xz_xx_0_y, g_xx_0_0_0_xz_xx_0_z, g_xxxz_xx_0_x, g_xxxz_xx_0_y, g_xxxz_xx_0_z, g_xz_xx_0_x, g_xz_xx_0_y, g_xz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_xx_0_x[i] = -6.0 * g_xz_xx_0_x[i] * a_exp + 4.0 * g_xxxz_xx_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xx_0_y[i] = -6.0 * g_xz_xx_0_y[i] * a_exp + 4.0 * g_xxxz_xx_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xx_0_z[i] = -6.0 * g_xz_xx_0_z[i] * a_exp + 4.0 * g_xxxz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_xx_0_0_0_xz_xy_0_x, g_xx_0_0_0_xz_xy_0_y, g_xx_0_0_0_xz_xy_0_z, g_xxxz_xy_0_x, g_xxxz_xy_0_y, g_xxxz_xy_0_z, g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_xy_0_x[i] = -6.0 * g_xz_xy_0_x[i] * a_exp + 4.0 * g_xxxz_xy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xy_0_y[i] = -6.0 * g_xz_xy_0_y[i] * a_exp + 4.0 * g_xxxz_xy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xy_0_z[i] = -6.0 * g_xz_xy_0_z[i] * a_exp + 4.0 * g_xxxz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_xx_0_0_0_xz_xz_0_x, g_xx_0_0_0_xz_xz_0_y, g_xx_0_0_0_xz_xz_0_z, g_xxxz_xz_0_x, g_xxxz_xz_0_y, g_xxxz_xz_0_z, g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_xz_0_x[i] = -6.0 * g_xz_xz_0_x[i] * a_exp + 4.0 * g_xxxz_xz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xz_0_y[i] = -6.0 * g_xz_xz_0_y[i] * a_exp + 4.0 * g_xxxz_xz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xz_0_z[i] = -6.0 * g_xz_xz_0_z[i] * a_exp + 4.0 * g_xxxz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_xx_0_0_0_xz_yy_0_x, g_xx_0_0_0_xz_yy_0_y, g_xx_0_0_0_xz_yy_0_z, g_xxxz_yy_0_x, g_xxxz_yy_0_y, g_xxxz_yy_0_z, g_xz_yy_0_x, g_xz_yy_0_y, g_xz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_yy_0_x[i] = -6.0 * g_xz_yy_0_x[i] * a_exp + 4.0 * g_xxxz_yy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yy_0_y[i] = -6.0 * g_xz_yy_0_y[i] * a_exp + 4.0 * g_xxxz_yy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yy_0_z[i] = -6.0 * g_xz_yy_0_z[i] * a_exp + 4.0 * g_xxxz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_xx_0_0_0_xz_yz_0_x, g_xx_0_0_0_xz_yz_0_y, g_xx_0_0_0_xz_yz_0_z, g_xxxz_yz_0_x, g_xxxz_yz_0_y, g_xxxz_yz_0_z, g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_yz_0_x[i] = -6.0 * g_xz_yz_0_x[i] * a_exp + 4.0 * g_xxxz_yz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yz_0_y[i] = -6.0 * g_xz_yz_0_y[i] * a_exp + 4.0 * g_xxxz_yz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yz_0_z[i] = -6.0 * g_xz_yz_0_z[i] * a_exp + 4.0 * g_xxxz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_xx_0_0_0_xz_zz_0_x, g_xx_0_0_0_xz_zz_0_y, g_xx_0_0_0_xz_zz_0_z, g_xxxz_zz_0_x, g_xxxz_zz_0_y, g_xxxz_zz_0_z, g_xz_zz_0_x, g_xz_zz_0_y, g_xz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_zz_0_x[i] = -6.0 * g_xz_zz_0_x[i] * a_exp + 4.0 * g_xxxz_zz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_zz_0_y[i] = -6.0 * g_xz_zz_0_y[i] * a_exp + 4.0 * g_xxxz_zz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_zz_0_z[i] = -6.0 * g_xz_zz_0_z[i] * a_exp + 4.0 * g_xxxz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_xx_0_0_0_yy_xx_0_x, g_xx_0_0_0_yy_xx_0_y, g_xx_0_0_0_yy_xx_0_z, g_xxyy_xx_0_x, g_xxyy_xx_0_y, g_xxyy_xx_0_z, g_yy_xx_0_x, g_yy_xx_0_y, g_yy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_xx_0_x[i] = -2.0 * g_yy_xx_0_x[i] * a_exp + 4.0 * g_xxyy_xx_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xx_0_y[i] = -2.0 * g_yy_xx_0_y[i] * a_exp + 4.0 * g_xxyy_xx_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xx_0_z[i] = -2.0 * g_yy_xx_0_z[i] * a_exp + 4.0 * g_xxyy_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_xx_0_0_0_yy_xy_0_x, g_xx_0_0_0_yy_xy_0_y, g_xx_0_0_0_yy_xy_0_z, g_xxyy_xy_0_x, g_xxyy_xy_0_y, g_xxyy_xy_0_z, g_yy_xy_0_x, g_yy_xy_0_y, g_yy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_xy_0_x[i] = -2.0 * g_yy_xy_0_x[i] * a_exp + 4.0 * g_xxyy_xy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xy_0_y[i] = -2.0 * g_yy_xy_0_y[i] * a_exp + 4.0 * g_xxyy_xy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xy_0_z[i] = -2.0 * g_yy_xy_0_z[i] * a_exp + 4.0 * g_xxyy_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_xx_0_0_0_yy_xz_0_x, g_xx_0_0_0_yy_xz_0_y, g_xx_0_0_0_yy_xz_0_z, g_xxyy_xz_0_x, g_xxyy_xz_0_y, g_xxyy_xz_0_z, g_yy_xz_0_x, g_yy_xz_0_y, g_yy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_xz_0_x[i] = -2.0 * g_yy_xz_0_x[i] * a_exp + 4.0 * g_xxyy_xz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xz_0_y[i] = -2.0 * g_yy_xz_0_y[i] * a_exp + 4.0 * g_xxyy_xz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xz_0_z[i] = -2.0 * g_yy_xz_0_z[i] * a_exp + 4.0 * g_xxyy_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_xx_0_0_0_yy_yy_0_x, g_xx_0_0_0_yy_yy_0_y, g_xx_0_0_0_yy_yy_0_z, g_xxyy_yy_0_x, g_xxyy_yy_0_y, g_xxyy_yy_0_z, g_yy_yy_0_x, g_yy_yy_0_y, g_yy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_yy_0_x[i] = -2.0 * g_yy_yy_0_x[i] * a_exp + 4.0 * g_xxyy_yy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yy_0_y[i] = -2.0 * g_yy_yy_0_y[i] * a_exp + 4.0 * g_xxyy_yy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yy_0_z[i] = -2.0 * g_yy_yy_0_z[i] * a_exp + 4.0 * g_xxyy_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_xx_0_0_0_yy_yz_0_x, g_xx_0_0_0_yy_yz_0_y, g_xx_0_0_0_yy_yz_0_z, g_xxyy_yz_0_x, g_xxyy_yz_0_y, g_xxyy_yz_0_z, g_yy_yz_0_x, g_yy_yz_0_y, g_yy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_yz_0_x[i] = -2.0 * g_yy_yz_0_x[i] * a_exp + 4.0 * g_xxyy_yz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yz_0_y[i] = -2.0 * g_yy_yz_0_y[i] * a_exp + 4.0 * g_xxyy_yz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yz_0_z[i] = -2.0 * g_yy_yz_0_z[i] * a_exp + 4.0 * g_xxyy_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_xx_0_0_0_yy_zz_0_x, g_xx_0_0_0_yy_zz_0_y, g_xx_0_0_0_yy_zz_0_z, g_xxyy_zz_0_x, g_xxyy_zz_0_y, g_xxyy_zz_0_z, g_yy_zz_0_x, g_yy_zz_0_y, g_yy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_zz_0_x[i] = -2.0 * g_yy_zz_0_x[i] * a_exp + 4.0 * g_xxyy_zz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_zz_0_y[i] = -2.0 * g_yy_zz_0_y[i] * a_exp + 4.0 * g_xxyy_zz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_zz_0_z[i] = -2.0 * g_yy_zz_0_z[i] * a_exp + 4.0 * g_xxyy_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_xx_0_0_0_yz_xx_0_x, g_xx_0_0_0_yz_xx_0_y, g_xx_0_0_0_yz_xx_0_z, g_xxyz_xx_0_x, g_xxyz_xx_0_y, g_xxyz_xx_0_z, g_yz_xx_0_x, g_yz_xx_0_y, g_yz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_xx_0_x[i] = -2.0 * g_yz_xx_0_x[i] * a_exp + 4.0 * g_xxyz_xx_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xx_0_y[i] = -2.0 * g_yz_xx_0_y[i] * a_exp + 4.0 * g_xxyz_xx_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xx_0_z[i] = -2.0 * g_yz_xx_0_z[i] * a_exp + 4.0 * g_xxyz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_xx_0_0_0_yz_xy_0_x, g_xx_0_0_0_yz_xy_0_y, g_xx_0_0_0_yz_xy_0_z, g_xxyz_xy_0_x, g_xxyz_xy_0_y, g_xxyz_xy_0_z, g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_xy_0_x[i] = -2.0 * g_yz_xy_0_x[i] * a_exp + 4.0 * g_xxyz_xy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xy_0_y[i] = -2.0 * g_yz_xy_0_y[i] * a_exp + 4.0 * g_xxyz_xy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xy_0_z[i] = -2.0 * g_yz_xy_0_z[i] * a_exp + 4.0 * g_xxyz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_xx_0_0_0_yz_xz_0_x, g_xx_0_0_0_yz_xz_0_y, g_xx_0_0_0_yz_xz_0_z, g_xxyz_xz_0_x, g_xxyz_xz_0_y, g_xxyz_xz_0_z, g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_xz_0_x[i] = -2.0 * g_yz_xz_0_x[i] * a_exp + 4.0 * g_xxyz_xz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xz_0_y[i] = -2.0 * g_yz_xz_0_y[i] * a_exp + 4.0 * g_xxyz_xz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xz_0_z[i] = -2.0 * g_yz_xz_0_z[i] * a_exp + 4.0 * g_xxyz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_xx_0_0_0_yz_yy_0_x, g_xx_0_0_0_yz_yy_0_y, g_xx_0_0_0_yz_yy_0_z, g_xxyz_yy_0_x, g_xxyz_yy_0_y, g_xxyz_yy_0_z, g_yz_yy_0_x, g_yz_yy_0_y, g_yz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_yy_0_x[i] = -2.0 * g_yz_yy_0_x[i] * a_exp + 4.0 * g_xxyz_yy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yy_0_y[i] = -2.0 * g_yz_yy_0_y[i] * a_exp + 4.0 * g_xxyz_yy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yy_0_z[i] = -2.0 * g_yz_yy_0_z[i] * a_exp + 4.0 * g_xxyz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_xx_0_0_0_yz_yz_0_x, g_xx_0_0_0_yz_yz_0_y, g_xx_0_0_0_yz_yz_0_z, g_xxyz_yz_0_x, g_xxyz_yz_0_y, g_xxyz_yz_0_z, g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_yz_0_x[i] = -2.0 * g_yz_yz_0_x[i] * a_exp + 4.0 * g_xxyz_yz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yz_0_y[i] = -2.0 * g_yz_yz_0_y[i] * a_exp + 4.0 * g_xxyz_yz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yz_0_z[i] = -2.0 * g_yz_yz_0_z[i] * a_exp + 4.0 * g_xxyz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_xx_0_0_0_yz_zz_0_x, g_xx_0_0_0_yz_zz_0_y, g_xx_0_0_0_yz_zz_0_z, g_xxyz_zz_0_x, g_xxyz_zz_0_y, g_xxyz_zz_0_z, g_yz_zz_0_x, g_yz_zz_0_y, g_yz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_zz_0_x[i] = -2.0 * g_yz_zz_0_x[i] * a_exp + 4.0 * g_xxyz_zz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_zz_0_y[i] = -2.0 * g_yz_zz_0_y[i] * a_exp + 4.0 * g_xxyz_zz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_zz_0_z[i] = -2.0 * g_yz_zz_0_z[i] * a_exp + 4.0 * g_xxyz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_xx_0_0_0_zz_xx_0_x, g_xx_0_0_0_zz_xx_0_y, g_xx_0_0_0_zz_xx_0_z, g_xxzz_xx_0_x, g_xxzz_xx_0_y, g_xxzz_xx_0_z, g_zz_xx_0_x, g_zz_xx_0_y, g_zz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_xx_0_x[i] = -2.0 * g_zz_xx_0_x[i] * a_exp + 4.0 * g_xxzz_xx_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xx_0_y[i] = -2.0 * g_zz_xx_0_y[i] * a_exp + 4.0 * g_xxzz_xx_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xx_0_z[i] = -2.0 * g_zz_xx_0_z[i] * a_exp + 4.0 * g_xxzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_xx_0_0_0_zz_xy_0_x, g_xx_0_0_0_zz_xy_0_y, g_xx_0_0_0_zz_xy_0_z, g_xxzz_xy_0_x, g_xxzz_xy_0_y, g_xxzz_xy_0_z, g_zz_xy_0_x, g_zz_xy_0_y, g_zz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_xy_0_x[i] = -2.0 * g_zz_xy_0_x[i] * a_exp + 4.0 * g_xxzz_xy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xy_0_y[i] = -2.0 * g_zz_xy_0_y[i] * a_exp + 4.0 * g_xxzz_xy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xy_0_z[i] = -2.0 * g_zz_xy_0_z[i] * a_exp + 4.0 * g_xxzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_xx_0_0_0_zz_xz_0_x, g_xx_0_0_0_zz_xz_0_y, g_xx_0_0_0_zz_xz_0_z, g_xxzz_xz_0_x, g_xxzz_xz_0_y, g_xxzz_xz_0_z, g_zz_xz_0_x, g_zz_xz_0_y, g_zz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_xz_0_x[i] = -2.0 * g_zz_xz_0_x[i] * a_exp + 4.0 * g_xxzz_xz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xz_0_y[i] = -2.0 * g_zz_xz_0_y[i] * a_exp + 4.0 * g_xxzz_xz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xz_0_z[i] = -2.0 * g_zz_xz_0_z[i] * a_exp + 4.0 * g_xxzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_xx_0_0_0_zz_yy_0_x, g_xx_0_0_0_zz_yy_0_y, g_xx_0_0_0_zz_yy_0_z, g_xxzz_yy_0_x, g_xxzz_yy_0_y, g_xxzz_yy_0_z, g_zz_yy_0_x, g_zz_yy_0_y, g_zz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_yy_0_x[i] = -2.0 * g_zz_yy_0_x[i] * a_exp + 4.0 * g_xxzz_yy_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yy_0_y[i] = -2.0 * g_zz_yy_0_y[i] * a_exp + 4.0 * g_xxzz_yy_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yy_0_z[i] = -2.0 * g_zz_yy_0_z[i] * a_exp + 4.0 * g_xxzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_xx_0_0_0_zz_yz_0_x, g_xx_0_0_0_zz_yz_0_y, g_xx_0_0_0_zz_yz_0_z, g_xxzz_yz_0_x, g_xxzz_yz_0_y, g_xxzz_yz_0_z, g_zz_yz_0_x, g_zz_yz_0_y, g_zz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_yz_0_x[i] = -2.0 * g_zz_yz_0_x[i] * a_exp + 4.0 * g_xxzz_yz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yz_0_y[i] = -2.0 * g_zz_yz_0_y[i] * a_exp + 4.0 * g_xxzz_yz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yz_0_z[i] = -2.0 * g_zz_yz_0_z[i] * a_exp + 4.0 * g_xxzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_xx_0_0_0_zz_zz_0_x, g_xx_0_0_0_zz_zz_0_y, g_xx_0_0_0_zz_zz_0_z, g_xxzz_zz_0_x, g_xxzz_zz_0_y, g_xxzz_zz_0_z, g_zz_zz_0_x, g_zz_zz_0_y, g_zz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_zz_0_x[i] = -2.0 * g_zz_zz_0_x[i] * a_exp + 4.0 * g_xxzz_zz_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_zz_0_y[i] = -2.0 * g_zz_zz_0_y[i] * a_exp + 4.0 * g_xxzz_zz_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_zz_0_z[i] = -2.0 * g_zz_zz_0_z[i] * a_exp + 4.0 * g_xxzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_xxxy_xx_0_x, g_xxxy_xx_0_y, g_xxxy_xx_0_z, g_xy_0_0_0_xx_xx_0_x, g_xy_0_0_0_xx_xx_0_y, g_xy_0_0_0_xx_xx_0_z, g_xy_xx_0_x, g_xy_xx_0_y, g_xy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_xx_0_x[i] = -4.0 * g_xy_xx_0_x[i] * a_exp + 4.0 * g_xxxy_xx_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xx_0_y[i] = -4.0 * g_xy_xx_0_y[i] * a_exp + 4.0 * g_xxxy_xx_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xx_0_z[i] = -4.0 * g_xy_xx_0_z[i] * a_exp + 4.0 * g_xxxy_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_xxxy_xy_0_x, g_xxxy_xy_0_y, g_xxxy_xy_0_z, g_xy_0_0_0_xx_xy_0_x, g_xy_0_0_0_xx_xy_0_y, g_xy_0_0_0_xx_xy_0_z, g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_xy_0_x[i] = -4.0 * g_xy_xy_0_x[i] * a_exp + 4.0 * g_xxxy_xy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xy_0_y[i] = -4.0 * g_xy_xy_0_y[i] * a_exp + 4.0 * g_xxxy_xy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xy_0_z[i] = -4.0 * g_xy_xy_0_z[i] * a_exp + 4.0 * g_xxxy_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_xxxy_xz_0_x, g_xxxy_xz_0_y, g_xxxy_xz_0_z, g_xy_0_0_0_xx_xz_0_x, g_xy_0_0_0_xx_xz_0_y, g_xy_0_0_0_xx_xz_0_z, g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_xz_0_x[i] = -4.0 * g_xy_xz_0_x[i] * a_exp + 4.0 * g_xxxy_xz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xz_0_y[i] = -4.0 * g_xy_xz_0_y[i] * a_exp + 4.0 * g_xxxy_xz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xz_0_z[i] = -4.0 * g_xy_xz_0_z[i] * a_exp + 4.0 * g_xxxy_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_xxxy_yy_0_x, g_xxxy_yy_0_y, g_xxxy_yy_0_z, g_xy_0_0_0_xx_yy_0_x, g_xy_0_0_0_xx_yy_0_y, g_xy_0_0_0_xx_yy_0_z, g_xy_yy_0_x, g_xy_yy_0_y, g_xy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_yy_0_x[i] = -4.0 * g_xy_yy_0_x[i] * a_exp + 4.0 * g_xxxy_yy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yy_0_y[i] = -4.0 * g_xy_yy_0_y[i] * a_exp + 4.0 * g_xxxy_yy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yy_0_z[i] = -4.0 * g_xy_yy_0_z[i] * a_exp + 4.0 * g_xxxy_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_xxxy_yz_0_x, g_xxxy_yz_0_y, g_xxxy_yz_0_z, g_xy_0_0_0_xx_yz_0_x, g_xy_0_0_0_xx_yz_0_y, g_xy_0_0_0_xx_yz_0_z, g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_yz_0_x[i] = -4.0 * g_xy_yz_0_x[i] * a_exp + 4.0 * g_xxxy_yz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yz_0_y[i] = -4.0 * g_xy_yz_0_y[i] * a_exp + 4.0 * g_xxxy_yz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yz_0_z[i] = -4.0 * g_xy_yz_0_z[i] * a_exp + 4.0 * g_xxxy_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_xxxy_zz_0_x, g_xxxy_zz_0_y, g_xxxy_zz_0_z, g_xy_0_0_0_xx_zz_0_x, g_xy_0_0_0_xx_zz_0_y, g_xy_0_0_0_xx_zz_0_z, g_xy_zz_0_x, g_xy_zz_0_y, g_xy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_zz_0_x[i] = -4.0 * g_xy_zz_0_x[i] * a_exp + 4.0 * g_xxxy_zz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_zz_0_y[i] = -4.0 * g_xy_zz_0_y[i] * a_exp + 4.0 * g_xxxy_zz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_zz_0_z[i] = -4.0 * g_xy_zz_0_z[i] * a_exp + 4.0 * g_xxxy_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_xx_xx_0_x, g_xx_xx_0_y, g_xx_xx_0_z, g_xxyy_xx_0_x, g_xxyy_xx_0_y, g_xxyy_xx_0_z, g_xy_0_0_0_xy_xx_0_x, g_xy_0_0_0_xy_xx_0_y, g_xy_0_0_0_xy_xx_0_z, g_yy_xx_0_x, g_yy_xx_0_y, g_yy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_xx_0_x[i] = g_0_xx_0_x[i] - 2.0 * g_yy_xx_0_x[i] * a_exp - 2.0 * g_xx_xx_0_x[i] * a_exp + 4.0 * g_xxyy_xx_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xx_0_y[i] = g_0_xx_0_y[i] - 2.0 * g_yy_xx_0_y[i] * a_exp - 2.0 * g_xx_xx_0_y[i] * a_exp + 4.0 * g_xxyy_xx_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xx_0_z[i] = g_0_xx_0_z[i] - 2.0 * g_yy_xx_0_z[i] * a_exp - 2.0 * g_xx_xx_0_z[i] * a_exp + 4.0 * g_xxyy_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_xx_xy_0_x, g_xx_xy_0_y, g_xx_xy_0_z, g_xxyy_xy_0_x, g_xxyy_xy_0_y, g_xxyy_xy_0_z, g_xy_0_0_0_xy_xy_0_x, g_xy_0_0_0_xy_xy_0_y, g_xy_0_0_0_xy_xy_0_z, g_yy_xy_0_x, g_yy_xy_0_y, g_yy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_xy_0_x[i] = g_0_xy_0_x[i] - 2.0 * g_yy_xy_0_x[i] * a_exp - 2.0 * g_xx_xy_0_x[i] * a_exp + 4.0 * g_xxyy_xy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xy_0_y[i] = g_0_xy_0_y[i] - 2.0 * g_yy_xy_0_y[i] * a_exp - 2.0 * g_xx_xy_0_y[i] * a_exp + 4.0 * g_xxyy_xy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xy_0_z[i] = g_0_xy_0_z[i] - 2.0 * g_yy_xy_0_z[i] * a_exp - 2.0 * g_xx_xy_0_z[i] * a_exp + 4.0 * g_xxyy_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_xx_xz_0_x, g_xx_xz_0_y, g_xx_xz_0_z, g_xxyy_xz_0_x, g_xxyy_xz_0_y, g_xxyy_xz_0_z, g_xy_0_0_0_xy_xz_0_x, g_xy_0_0_0_xy_xz_0_y, g_xy_0_0_0_xy_xz_0_z, g_yy_xz_0_x, g_yy_xz_0_y, g_yy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_xz_0_x[i] = g_0_xz_0_x[i] - 2.0 * g_yy_xz_0_x[i] * a_exp - 2.0 * g_xx_xz_0_x[i] * a_exp + 4.0 * g_xxyy_xz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xz_0_y[i] = g_0_xz_0_y[i] - 2.0 * g_yy_xz_0_y[i] * a_exp - 2.0 * g_xx_xz_0_y[i] * a_exp + 4.0 * g_xxyy_xz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xz_0_z[i] = g_0_xz_0_z[i] - 2.0 * g_yy_xz_0_z[i] * a_exp - 2.0 * g_xx_xz_0_z[i] * a_exp + 4.0 * g_xxyy_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_xx_yy_0_x, g_xx_yy_0_y, g_xx_yy_0_z, g_xxyy_yy_0_x, g_xxyy_yy_0_y, g_xxyy_yy_0_z, g_xy_0_0_0_xy_yy_0_x, g_xy_0_0_0_xy_yy_0_y, g_xy_0_0_0_xy_yy_0_z, g_yy_yy_0_x, g_yy_yy_0_y, g_yy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_yy_0_x[i] = g_0_yy_0_x[i] - 2.0 * g_yy_yy_0_x[i] * a_exp - 2.0 * g_xx_yy_0_x[i] * a_exp + 4.0 * g_xxyy_yy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yy_0_y[i] = g_0_yy_0_y[i] - 2.0 * g_yy_yy_0_y[i] * a_exp - 2.0 * g_xx_yy_0_y[i] * a_exp + 4.0 * g_xxyy_yy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yy_0_z[i] = g_0_yy_0_z[i] - 2.0 * g_yy_yy_0_z[i] * a_exp - 2.0 * g_xx_yy_0_z[i] * a_exp + 4.0 * g_xxyy_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_xx_yz_0_x, g_xx_yz_0_y, g_xx_yz_0_z, g_xxyy_yz_0_x, g_xxyy_yz_0_y, g_xxyy_yz_0_z, g_xy_0_0_0_xy_yz_0_x, g_xy_0_0_0_xy_yz_0_y, g_xy_0_0_0_xy_yz_0_z, g_yy_yz_0_x, g_yy_yz_0_y, g_yy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_yz_0_x[i] = g_0_yz_0_x[i] - 2.0 * g_yy_yz_0_x[i] * a_exp - 2.0 * g_xx_yz_0_x[i] * a_exp + 4.0 * g_xxyy_yz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yz_0_y[i] = g_0_yz_0_y[i] - 2.0 * g_yy_yz_0_y[i] * a_exp - 2.0 * g_xx_yz_0_y[i] * a_exp + 4.0 * g_xxyy_yz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yz_0_z[i] = g_0_yz_0_z[i] - 2.0 * g_yy_yz_0_z[i] * a_exp - 2.0 * g_xx_yz_0_z[i] * a_exp + 4.0 * g_xxyy_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_xx_zz_0_x, g_xx_zz_0_y, g_xx_zz_0_z, g_xxyy_zz_0_x, g_xxyy_zz_0_y, g_xxyy_zz_0_z, g_xy_0_0_0_xy_zz_0_x, g_xy_0_0_0_xy_zz_0_y, g_xy_0_0_0_xy_zz_0_z, g_yy_zz_0_x, g_yy_zz_0_y, g_yy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_zz_0_x[i] = g_0_zz_0_x[i] - 2.0 * g_yy_zz_0_x[i] * a_exp - 2.0 * g_xx_zz_0_x[i] * a_exp + 4.0 * g_xxyy_zz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_zz_0_y[i] = g_0_zz_0_y[i] - 2.0 * g_yy_zz_0_y[i] * a_exp - 2.0 * g_xx_zz_0_y[i] * a_exp + 4.0 * g_xxyy_zz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_zz_0_z[i] = g_0_zz_0_z[i] - 2.0 * g_yy_zz_0_z[i] * a_exp - 2.0 * g_xx_zz_0_z[i] * a_exp + 4.0 * g_xxyy_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_xxyz_xx_0_x, g_xxyz_xx_0_y, g_xxyz_xx_0_z, g_xy_0_0_0_xz_xx_0_x, g_xy_0_0_0_xz_xx_0_y, g_xy_0_0_0_xz_xx_0_z, g_yz_xx_0_x, g_yz_xx_0_y, g_yz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_xx_0_x[i] = -2.0 * g_yz_xx_0_x[i] * a_exp + 4.0 * g_xxyz_xx_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xx_0_y[i] = -2.0 * g_yz_xx_0_y[i] * a_exp + 4.0 * g_xxyz_xx_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xx_0_z[i] = -2.0 * g_yz_xx_0_z[i] * a_exp + 4.0 * g_xxyz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_xxyz_xy_0_x, g_xxyz_xy_0_y, g_xxyz_xy_0_z, g_xy_0_0_0_xz_xy_0_x, g_xy_0_0_0_xz_xy_0_y, g_xy_0_0_0_xz_xy_0_z, g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_xy_0_x[i] = -2.0 * g_yz_xy_0_x[i] * a_exp + 4.0 * g_xxyz_xy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xy_0_y[i] = -2.0 * g_yz_xy_0_y[i] * a_exp + 4.0 * g_xxyz_xy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xy_0_z[i] = -2.0 * g_yz_xy_0_z[i] * a_exp + 4.0 * g_xxyz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_xxyz_xz_0_x, g_xxyz_xz_0_y, g_xxyz_xz_0_z, g_xy_0_0_0_xz_xz_0_x, g_xy_0_0_0_xz_xz_0_y, g_xy_0_0_0_xz_xz_0_z, g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_xz_0_x[i] = -2.0 * g_yz_xz_0_x[i] * a_exp + 4.0 * g_xxyz_xz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xz_0_y[i] = -2.0 * g_yz_xz_0_y[i] * a_exp + 4.0 * g_xxyz_xz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xz_0_z[i] = -2.0 * g_yz_xz_0_z[i] * a_exp + 4.0 * g_xxyz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_xxyz_yy_0_x, g_xxyz_yy_0_y, g_xxyz_yy_0_z, g_xy_0_0_0_xz_yy_0_x, g_xy_0_0_0_xz_yy_0_y, g_xy_0_0_0_xz_yy_0_z, g_yz_yy_0_x, g_yz_yy_0_y, g_yz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_yy_0_x[i] = -2.0 * g_yz_yy_0_x[i] * a_exp + 4.0 * g_xxyz_yy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yy_0_y[i] = -2.0 * g_yz_yy_0_y[i] * a_exp + 4.0 * g_xxyz_yy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yy_0_z[i] = -2.0 * g_yz_yy_0_z[i] * a_exp + 4.0 * g_xxyz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_xxyz_yz_0_x, g_xxyz_yz_0_y, g_xxyz_yz_0_z, g_xy_0_0_0_xz_yz_0_x, g_xy_0_0_0_xz_yz_0_y, g_xy_0_0_0_xz_yz_0_z, g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_yz_0_x[i] = -2.0 * g_yz_yz_0_x[i] * a_exp + 4.0 * g_xxyz_yz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yz_0_y[i] = -2.0 * g_yz_yz_0_y[i] * a_exp + 4.0 * g_xxyz_yz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yz_0_z[i] = -2.0 * g_yz_yz_0_z[i] * a_exp + 4.0 * g_xxyz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_xxyz_zz_0_x, g_xxyz_zz_0_y, g_xxyz_zz_0_z, g_xy_0_0_0_xz_zz_0_x, g_xy_0_0_0_xz_zz_0_y, g_xy_0_0_0_xz_zz_0_z, g_yz_zz_0_x, g_yz_zz_0_y, g_yz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_zz_0_x[i] = -2.0 * g_yz_zz_0_x[i] * a_exp + 4.0 * g_xxyz_zz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_zz_0_y[i] = -2.0 * g_yz_zz_0_y[i] * a_exp + 4.0 * g_xxyz_zz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_zz_0_z[i] = -2.0 * g_yz_zz_0_z[i] * a_exp + 4.0 * g_xxyz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_xy_0_0_0_yy_xx_0_x, g_xy_0_0_0_yy_xx_0_y, g_xy_0_0_0_yy_xx_0_z, g_xy_xx_0_x, g_xy_xx_0_y, g_xy_xx_0_z, g_xyyy_xx_0_x, g_xyyy_xx_0_y, g_xyyy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_xx_0_x[i] = -4.0 * g_xy_xx_0_x[i] * a_exp + 4.0 * g_xyyy_xx_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xx_0_y[i] = -4.0 * g_xy_xx_0_y[i] * a_exp + 4.0 * g_xyyy_xx_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xx_0_z[i] = -4.0 * g_xy_xx_0_z[i] * a_exp + 4.0 * g_xyyy_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_xy_0_0_0_yy_xy_0_x, g_xy_0_0_0_yy_xy_0_y, g_xy_0_0_0_yy_xy_0_z, g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z, g_xyyy_xy_0_x, g_xyyy_xy_0_y, g_xyyy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_xy_0_x[i] = -4.0 * g_xy_xy_0_x[i] * a_exp + 4.0 * g_xyyy_xy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xy_0_y[i] = -4.0 * g_xy_xy_0_y[i] * a_exp + 4.0 * g_xyyy_xy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xy_0_z[i] = -4.0 * g_xy_xy_0_z[i] * a_exp + 4.0 * g_xyyy_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_xy_0_0_0_yy_xz_0_x, g_xy_0_0_0_yy_xz_0_y, g_xy_0_0_0_yy_xz_0_z, g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z, g_xyyy_xz_0_x, g_xyyy_xz_0_y, g_xyyy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_xz_0_x[i] = -4.0 * g_xy_xz_0_x[i] * a_exp + 4.0 * g_xyyy_xz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xz_0_y[i] = -4.0 * g_xy_xz_0_y[i] * a_exp + 4.0 * g_xyyy_xz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xz_0_z[i] = -4.0 * g_xy_xz_0_z[i] * a_exp + 4.0 * g_xyyy_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_xy_0_0_0_yy_yy_0_x, g_xy_0_0_0_yy_yy_0_y, g_xy_0_0_0_yy_yy_0_z, g_xy_yy_0_x, g_xy_yy_0_y, g_xy_yy_0_z, g_xyyy_yy_0_x, g_xyyy_yy_0_y, g_xyyy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_yy_0_x[i] = -4.0 * g_xy_yy_0_x[i] * a_exp + 4.0 * g_xyyy_yy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yy_0_y[i] = -4.0 * g_xy_yy_0_y[i] * a_exp + 4.0 * g_xyyy_yy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yy_0_z[i] = -4.0 * g_xy_yy_0_z[i] * a_exp + 4.0 * g_xyyy_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_xy_0_0_0_yy_yz_0_x, g_xy_0_0_0_yy_yz_0_y, g_xy_0_0_0_yy_yz_0_z, g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z, g_xyyy_yz_0_x, g_xyyy_yz_0_y, g_xyyy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_yz_0_x[i] = -4.0 * g_xy_yz_0_x[i] * a_exp + 4.0 * g_xyyy_yz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yz_0_y[i] = -4.0 * g_xy_yz_0_y[i] * a_exp + 4.0 * g_xyyy_yz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yz_0_z[i] = -4.0 * g_xy_yz_0_z[i] * a_exp + 4.0 * g_xyyy_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_xy_0_0_0_yy_zz_0_x, g_xy_0_0_0_yy_zz_0_y, g_xy_0_0_0_yy_zz_0_z, g_xy_zz_0_x, g_xy_zz_0_y, g_xy_zz_0_z, g_xyyy_zz_0_x, g_xyyy_zz_0_y, g_xyyy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_zz_0_x[i] = -4.0 * g_xy_zz_0_x[i] * a_exp + 4.0 * g_xyyy_zz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_zz_0_y[i] = -4.0 * g_xy_zz_0_y[i] * a_exp + 4.0 * g_xyyy_zz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_zz_0_z[i] = -4.0 * g_xy_zz_0_z[i] * a_exp + 4.0 * g_xyyy_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_xy_0_0_0_yz_xx_0_x, g_xy_0_0_0_yz_xx_0_y, g_xy_0_0_0_yz_xx_0_z, g_xyyz_xx_0_x, g_xyyz_xx_0_y, g_xyyz_xx_0_z, g_xz_xx_0_x, g_xz_xx_0_y, g_xz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_xx_0_x[i] = -2.0 * g_xz_xx_0_x[i] * a_exp + 4.0 * g_xyyz_xx_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xx_0_y[i] = -2.0 * g_xz_xx_0_y[i] * a_exp + 4.0 * g_xyyz_xx_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xx_0_z[i] = -2.0 * g_xz_xx_0_z[i] * a_exp + 4.0 * g_xyyz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_xy_0_0_0_yz_xy_0_x, g_xy_0_0_0_yz_xy_0_y, g_xy_0_0_0_yz_xy_0_z, g_xyyz_xy_0_x, g_xyyz_xy_0_y, g_xyyz_xy_0_z, g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_xy_0_x[i] = -2.0 * g_xz_xy_0_x[i] * a_exp + 4.0 * g_xyyz_xy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xy_0_y[i] = -2.0 * g_xz_xy_0_y[i] * a_exp + 4.0 * g_xyyz_xy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xy_0_z[i] = -2.0 * g_xz_xy_0_z[i] * a_exp + 4.0 * g_xyyz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_xy_0_0_0_yz_xz_0_x, g_xy_0_0_0_yz_xz_0_y, g_xy_0_0_0_yz_xz_0_z, g_xyyz_xz_0_x, g_xyyz_xz_0_y, g_xyyz_xz_0_z, g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_xz_0_x[i] = -2.0 * g_xz_xz_0_x[i] * a_exp + 4.0 * g_xyyz_xz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xz_0_y[i] = -2.0 * g_xz_xz_0_y[i] * a_exp + 4.0 * g_xyyz_xz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xz_0_z[i] = -2.0 * g_xz_xz_0_z[i] * a_exp + 4.0 * g_xyyz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_xy_0_0_0_yz_yy_0_x, g_xy_0_0_0_yz_yy_0_y, g_xy_0_0_0_yz_yy_0_z, g_xyyz_yy_0_x, g_xyyz_yy_0_y, g_xyyz_yy_0_z, g_xz_yy_0_x, g_xz_yy_0_y, g_xz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_yy_0_x[i] = -2.0 * g_xz_yy_0_x[i] * a_exp + 4.0 * g_xyyz_yy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yy_0_y[i] = -2.0 * g_xz_yy_0_y[i] * a_exp + 4.0 * g_xyyz_yy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yy_0_z[i] = -2.0 * g_xz_yy_0_z[i] * a_exp + 4.0 * g_xyyz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_xy_0_0_0_yz_yz_0_x, g_xy_0_0_0_yz_yz_0_y, g_xy_0_0_0_yz_yz_0_z, g_xyyz_yz_0_x, g_xyyz_yz_0_y, g_xyyz_yz_0_z, g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_yz_0_x[i] = -2.0 * g_xz_yz_0_x[i] * a_exp + 4.0 * g_xyyz_yz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yz_0_y[i] = -2.0 * g_xz_yz_0_y[i] * a_exp + 4.0 * g_xyyz_yz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yz_0_z[i] = -2.0 * g_xz_yz_0_z[i] * a_exp + 4.0 * g_xyyz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_xy_0_0_0_yz_zz_0_x, g_xy_0_0_0_yz_zz_0_y, g_xy_0_0_0_yz_zz_0_z, g_xyyz_zz_0_x, g_xyyz_zz_0_y, g_xyyz_zz_0_z, g_xz_zz_0_x, g_xz_zz_0_y, g_xz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_zz_0_x[i] = -2.0 * g_xz_zz_0_x[i] * a_exp + 4.0 * g_xyyz_zz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_zz_0_y[i] = -2.0 * g_xz_zz_0_y[i] * a_exp + 4.0 * g_xyyz_zz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_zz_0_z[i] = -2.0 * g_xz_zz_0_z[i] * a_exp + 4.0 * g_xyyz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_xy_0_0_0_zz_xx_0_x, g_xy_0_0_0_zz_xx_0_y, g_xy_0_0_0_zz_xx_0_z, g_xyzz_xx_0_x, g_xyzz_xx_0_y, g_xyzz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_xx_0_x[i] = 4.0 * g_xyzz_xx_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xx_0_y[i] = 4.0 * g_xyzz_xx_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xx_0_z[i] = 4.0 * g_xyzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_xy_0_0_0_zz_xy_0_x, g_xy_0_0_0_zz_xy_0_y, g_xy_0_0_0_zz_xy_0_z, g_xyzz_xy_0_x, g_xyzz_xy_0_y, g_xyzz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_xy_0_x[i] = 4.0 * g_xyzz_xy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xy_0_y[i] = 4.0 * g_xyzz_xy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xy_0_z[i] = 4.0 * g_xyzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_xy_0_0_0_zz_xz_0_x, g_xy_0_0_0_zz_xz_0_y, g_xy_0_0_0_zz_xz_0_z, g_xyzz_xz_0_x, g_xyzz_xz_0_y, g_xyzz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_xz_0_x[i] = 4.0 * g_xyzz_xz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xz_0_y[i] = 4.0 * g_xyzz_xz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xz_0_z[i] = 4.0 * g_xyzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_xy_0_0_0_zz_yy_0_x, g_xy_0_0_0_zz_yy_0_y, g_xy_0_0_0_zz_yy_0_z, g_xyzz_yy_0_x, g_xyzz_yy_0_y, g_xyzz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_yy_0_x[i] = 4.0 * g_xyzz_yy_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yy_0_y[i] = 4.0 * g_xyzz_yy_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yy_0_z[i] = 4.0 * g_xyzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_xy_0_0_0_zz_yz_0_x, g_xy_0_0_0_zz_yz_0_y, g_xy_0_0_0_zz_yz_0_z, g_xyzz_yz_0_x, g_xyzz_yz_0_y, g_xyzz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_yz_0_x[i] = 4.0 * g_xyzz_yz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yz_0_y[i] = 4.0 * g_xyzz_yz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yz_0_z[i] = 4.0 * g_xyzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_xy_0_0_0_zz_zz_0_x, g_xy_0_0_0_zz_zz_0_y, g_xy_0_0_0_zz_zz_0_z, g_xyzz_zz_0_x, g_xyzz_zz_0_y, g_xyzz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_zz_0_x[i] = 4.0 * g_xyzz_zz_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_zz_0_y[i] = 4.0 * g_xyzz_zz_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_zz_0_z[i] = 4.0 * g_xyzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_xxxz_xx_0_x, g_xxxz_xx_0_y, g_xxxz_xx_0_z, g_xz_0_0_0_xx_xx_0_x, g_xz_0_0_0_xx_xx_0_y, g_xz_0_0_0_xx_xx_0_z, g_xz_xx_0_x, g_xz_xx_0_y, g_xz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_xx_0_x[i] = -4.0 * g_xz_xx_0_x[i] * a_exp + 4.0 * g_xxxz_xx_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xx_0_y[i] = -4.0 * g_xz_xx_0_y[i] * a_exp + 4.0 * g_xxxz_xx_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xx_0_z[i] = -4.0 * g_xz_xx_0_z[i] * a_exp + 4.0 * g_xxxz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_xxxz_xy_0_x, g_xxxz_xy_0_y, g_xxxz_xy_0_z, g_xz_0_0_0_xx_xy_0_x, g_xz_0_0_0_xx_xy_0_y, g_xz_0_0_0_xx_xy_0_z, g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_xy_0_x[i] = -4.0 * g_xz_xy_0_x[i] * a_exp + 4.0 * g_xxxz_xy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xy_0_y[i] = -4.0 * g_xz_xy_0_y[i] * a_exp + 4.0 * g_xxxz_xy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xy_0_z[i] = -4.0 * g_xz_xy_0_z[i] * a_exp + 4.0 * g_xxxz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_xxxz_xz_0_x, g_xxxz_xz_0_y, g_xxxz_xz_0_z, g_xz_0_0_0_xx_xz_0_x, g_xz_0_0_0_xx_xz_0_y, g_xz_0_0_0_xx_xz_0_z, g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_xz_0_x[i] = -4.0 * g_xz_xz_0_x[i] * a_exp + 4.0 * g_xxxz_xz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xz_0_y[i] = -4.0 * g_xz_xz_0_y[i] * a_exp + 4.0 * g_xxxz_xz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xz_0_z[i] = -4.0 * g_xz_xz_0_z[i] * a_exp + 4.0 * g_xxxz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_xxxz_yy_0_x, g_xxxz_yy_0_y, g_xxxz_yy_0_z, g_xz_0_0_0_xx_yy_0_x, g_xz_0_0_0_xx_yy_0_y, g_xz_0_0_0_xx_yy_0_z, g_xz_yy_0_x, g_xz_yy_0_y, g_xz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_yy_0_x[i] = -4.0 * g_xz_yy_0_x[i] * a_exp + 4.0 * g_xxxz_yy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yy_0_y[i] = -4.0 * g_xz_yy_0_y[i] * a_exp + 4.0 * g_xxxz_yy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yy_0_z[i] = -4.0 * g_xz_yy_0_z[i] * a_exp + 4.0 * g_xxxz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_xxxz_yz_0_x, g_xxxz_yz_0_y, g_xxxz_yz_0_z, g_xz_0_0_0_xx_yz_0_x, g_xz_0_0_0_xx_yz_0_y, g_xz_0_0_0_xx_yz_0_z, g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_yz_0_x[i] = -4.0 * g_xz_yz_0_x[i] * a_exp + 4.0 * g_xxxz_yz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yz_0_y[i] = -4.0 * g_xz_yz_0_y[i] * a_exp + 4.0 * g_xxxz_yz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yz_0_z[i] = -4.0 * g_xz_yz_0_z[i] * a_exp + 4.0 * g_xxxz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_xxxz_zz_0_x, g_xxxz_zz_0_y, g_xxxz_zz_0_z, g_xz_0_0_0_xx_zz_0_x, g_xz_0_0_0_xx_zz_0_y, g_xz_0_0_0_xx_zz_0_z, g_xz_zz_0_x, g_xz_zz_0_y, g_xz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_zz_0_x[i] = -4.0 * g_xz_zz_0_x[i] * a_exp + 4.0 * g_xxxz_zz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_zz_0_y[i] = -4.0 * g_xz_zz_0_y[i] * a_exp + 4.0 * g_xxxz_zz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_zz_0_z[i] = -4.0 * g_xz_zz_0_z[i] * a_exp + 4.0 * g_xxxz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_xxyz_xx_0_x, g_xxyz_xx_0_y, g_xxyz_xx_0_z, g_xz_0_0_0_xy_xx_0_x, g_xz_0_0_0_xy_xx_0_y, g_xz_0_0_0_xy_xx_0_z, g_yz_xx_0_x, g_yz_xx_0_y, g_yz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_xx_0_x[i] = -2.0 * g_yz_xx_0_x[i] * a_exp + 4.0 * g_xxyz_xx_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xx_0_y[i] = -2.0 * g_yz_xx_0_y[i] * a_exp + 4.0 * g_xxyz_xx_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xx_0_z[i] = -2.0 * g_yz_xx_0_z[i] * a_exp + 4.0 * g_xxyz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_xxyz_xy_0_x, g_xxyz_xy_0_y, g_xxyz_xy_0_z, g_xz_0_0_0_xy_xy_0_x, g_xz_0_0_0_xy_xy_0_y, g_xz_0_0_0_xy_xy_0_z, g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_xy_0_x[i] = -2.0 * g_yz_xy_0_x[i] * a_exp + 4.0 * g_xxyz_xy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xy_0_y[i] = -2.0 * g_yz_xy_0_y[i] * a_exp + 4.0 * g_xxyz_xy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xy_0_z[i] = -2.0 * g_yz_xy_0_z[i] * a_exp + 4.0 * g_xxyz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_xxyz_xz_0_x, g_xxyz_xz_0_y, g_xxyz_xz_0_z, g_xz_0_0_0_xy_xz_0_x, g_xz_0_0_0_xy_xz_0_y, g_xz_0_0_0_xy_xz_0_z, g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_xz_0_x[i] = -2.0 * g_yz_xz_0_x[i] * a_exp + 4.0 * g_xxyz_xz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xz_0_y[i] = -2.0 * g_yz_xz_0_y[i] * a_exp + 4.0 * g_xxyz_xz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xz_0_z[i] = -2.0 * g_yz_xz_0_z[i] * a_exp + 4.0 * g_xxyz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (243-246)

    #pragma omp simd aligned(g_xxyz_yy_0_x, g_xxyz_yy_0_y, g_xxyz_yy_0_z, g_xz_0_0_0_xy_yy_0_x, g_xz_0_0_0_xy_yy_0_y, g_xz_0_0_0_xy_yy_0_z, g_yz_yy_0_x, g_yz_yy_0_y, g_yz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_yy_0_x[i] = -2.0 * g_yz_yy_0_x[i] * a_exp + 4.0 * g_xxyz_yy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yy_0_y[i] = -2.0 * g_yz_yy_0_y[i] * a_exp + 4.0 * g_xxyz_yy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yy_0_z[i] = -2.0 * g_yz_yy_0_z[i] * a_exp + 4.0 * g_xxyz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (246-249)

    #pragma omp simd aligned(g_xxyz_yz_0_x, g_xxyz_yz_0_y, g_xxyz_yz_0_z, g_xz_0_0_0_xy_yz_0_x, g_xz_0_0_0_xy_yz_0_y, g_xz_0_0_0_xy_yz_0_z, g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_yz_0_x[i] = -2.0 * g_yz_yz_0_x[i] * a_exp + 4.0 * g_xxyz_yz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yz_0_y[i] = -2.0 * g_yz_yz_0_y[i] * a_exp + 4.0 * g_xxyz_yz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yz_0_z[i] = -2.0 * g_yz_yz_0_z[i] * a_exp + 4.0 * g_xxyz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (249-252)

    #pragma omp simd aligned(g_xxyz_zz_0_x, g_xxyz_zz_0_y, g_xxyz_zz_0_z, g_xz_0_0_0_xy_zz_0_x, g_xz_0_0_0_xy_zz_0_y, g_xz_0_0_0_xy_zz_0_z, g_yz_zz_0_x, g_yz_zz_0_y, g_yz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_zz_0_x[i] = -2.0 * g_yz_zz_0_x[i] * a_exp + 4.0 * g_xxyz_zz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_zz_0_y[i] = -2.0 * g_yz_zz_0_y[i] * a_exp + 4.0 * g_xxyz_zz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_zz_0_z[i] = -2.0 * g_yz_zz_0_z[i] * a_exp + 4.0 * g_xxyz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (252-255)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_xx_xx_0_x, g_xx_xx_0_y, g_xx_xx_0_z, g_xxzz_xx_0_x, g_xxzz_xx_0_y, g_xxzz_xx_0_z, g_xz_0_0_0_xz_xx_0_x, g_xz_0_0_0_xz_xx_0_y, g_xz_0_0_0_xz_xx_0_z, g_zz_xx_0_x, g_zz_xx_0_y, g_zz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_xx_0_x[i] = g_0_xx_0_x[i] - 2.0 * g_zz_xx_0_x[i] * a_exp - 2.0 * g_xx_xx_0_x[i] * a_exp + 4.0 * g_xxzz_xx_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xx_0_y[i] = g_0_xx_0_y[i] - 2.0 * g_zz_xx_0_y[i] * a_exp - 2.0 * g_xx_xx_0_y[i] * a_exp + 4.0 * g_xxzz_xx_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xx_0_z[i] = g_0_xx_0_z[i] - 2.0 * g_zz_xx_0_z[i] * a_exp - 2.0 * g_xx_xx_0_z[i] * a_exp + 4.0 * g_xxzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (255-258)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_xx_xy_0_x, g_xx_xy_0_y, g_xx_xy_0_z, g_xxzz_xy_0_x, g_xxzz_xy_0_y, g_xxzz_xy_0_z, g_xz_0_0_0_xz_xy_0_x, g_xz_0_0_0_xz_xy_0_y, g_xz_0_0_0_xz_xy_0_z, g_zz_xy_0_x, g_zz_xy_0_y, g_zz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_xy_0_x[i] = g_0_xy_0_x[i] - 2.0 * g_zz_xy_0_x[i] * a_exp - 2.0 * g_xx_xy_0_x[i] * a_exp + 4.0 * g_xxzz_xy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xy_0_y[i] = g_0_xy_0_y[i] - 2.0 * g_zz_xy_0_y[i] * a_exp - 2.0 * g_xx_xy_0_y[i] * a_exp + 4.0 * g_xxzz_xy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xy_0_z[i] = g_0_xy_0_z[i] - 2.0 * g_zz_xy_0_z[i] * a_exp - 2.0 * g_xx_xy_0_z[i] * a_exp + 4.0 * g_xxzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (258-261)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_xx_xz_0_x, g_xx_xz_0_y, g_xx_xz_0_z, g_xxzz_xz_0_x, g_xxzz_xz_0_y, g_xxzz_xz_0_z, g_xz_0_0_0_xz_xz_0_x, g_xz_0_0_0_xz_xz_0_y, g_xz_0_0_0_xz_xz_0_z, g_zz_xz_0_x, g_zz_xz_0_y, g_zz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_xz_0_x[i] = g_0_xz_0_x[i] - 2.0 * g_zz_xz_0_x[i] * a_exp - 2.0 * g_xx_xz_0_x[i] * a_exp + 4.0 * g_xxzz_xz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xz_0_y[i] = g_0_xz_0_y[i] - 2.0 * g_zz_xz_0_y[i] * a_exp - 2.0 * g_xx_xz_0_y[i] * a_exp + 4.0 * g_xxzz_xz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xz_0_z[i] = g_0_xz_0_z[i] - 2.0 * g_zz_xz_0_z[i] * a_exp - 2.0 * g_xx_xz_0_z[i] * a_exp + 4.0 * g_xxzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (261-264)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_xx_yy_0_x, g_xx_yy_0_y, g_xx_yy_0_z, g_xxzz_yy_0_x, g_xxzz_yy_0_y, g_xxzz_yy_0_z, g_xz_0_0_0_xz_yy_0_x, g_xz_0_0_0_xz_yy_0_y, g_xz_0_0_0_xz_yy_0_z, g_zz_yy_0_x, g_zz_yy_0_y, g_zz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_yy_0_x[i] = g_0_yy_0_x[i] - 2.0 * g_zz_yy_0_x[i] * a_exp - 2.0 * g_xx_yy_0_x[i] * a_exp + 4.0 * g_xxzz_yy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yy_0_y[i] = g_0_yy_0_y[i] - 2.0 * g_zz_yy_0_y[i] * a_exp - 2.0 * g_xx_yy_0_y[i] * a_exp + 4.0 * g_xxzz_yy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yy_0_z[i] = g_0_yy_0_z[i] - 2.0 * g_zz_yy_0_z[i] * a_exp - 2.0 * g_xx_yy_0_z[i] * a_exp + 4.0 * g_xxzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (264-267)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_xx_yz_0_x, g_xx_yz_0_y, g_xx_yz_0_z, g_xxzz_yz_0_x, g_xxzz_yz_0_y, g_xxzz_yz_0_z, g_xz_0_0_0_xz_yz_0_x, g_xz_0_0_0_xz_yz_0_y, g_xz_0_0_0_xz_yz_0_z, g_zz_yz_0_x, g_zz_yz_0_y, g_zz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_yz_0_x[i] = g_0_yz_0_x[i] - 2.0 * g_zz_yz_0_x[i] * a_exp - 2.0 * g_xx_yz_0_x[i] * a_exp + 4.0 * g_xxzz_yz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yz_0_y[i] = g_0_yz_0_y[i] - 2.0 * g_zz_yz_0_y[i] * a_exp - 2.0 * g_xx_yz_0_y[i] * a_exp + 4.0 * g_xxzz_yz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yz_0_z[i] = g_0_yz_0_z[i] - 2.0 * g_zz_yz_0_z[i] * a_exp - 2.0 * g_xx_yz_0_z[i] * a_exp + 4.0 * g_xxzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (267-270)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_xx_zz_0_x, g_xx_zz_0_y, g_xx_zz_0_z, g_xxzz_zz_0_x, g_xxzz_zz_0_y, g_xxzz_zz_0_z, g_xz_0_0_0_xz_zz_0_x, g_xz_0_0_0_xz_zz_0_y, g_xz_0_0_0_xz_zz_0_z, g_zz_zz_0_x, g_zz_zz_0_y, g_zz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_zz_0_x[i] = g_0_zz_0_x[i] - 2.0 * g_zz_zz_0_x[i] * a_exp - 2.0 * g_xx_zz_0_x[i] * a_exp + 4.0 * g_xxzz_zz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_zz_0_y[i] = g_0_zz_0_y[i] - 2.0 * g_zz_zz_0_y[i] * a_exp - 2.0 * g_xx_zz_0_y[i] * a_exp + 4.0 * g_xxzz_zz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_zz_0_z[i] = g_0_zz_0_z[i] - 2.0 * g_zz_zz_0_z[i] * a_exp - 2.0 * g_xx_zz_0_z[i] * a_exp + 4.0 * g_xxzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (270-273)

    #pragma omp simd aligned(g_xyyz_xx_0_x, g_xyyz_xx_0_y, g_xyyz_xx_0_z, g_xz_0_0_0_yy_xx_0_x, g_xz_0_0_0_yy_xx_0_y, g_xz_0_0_0_yy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_xx_0_x[i] = 4.0 * g_xyyz_xx_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xx_0_y[i] = 4.0 * g_xyyz_xx_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xx_0_z[i] = 4.0 * g_xyyz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (273-276)

    #pragma omp simd aligned(g_xyyz_xy_0_x, g_xyyz_xy_0_y, g_xyyz_xy_0_z, g_xz_0_0_0_yy_xy_0_x, g_xz_0_0_0_yy_xy_0_y, g_xz_0_0_0_yy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_xy_0_x[i] = 4.0 * g_xyyz_xy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xy_0_y[i] = 4.0 * g_xyyz_xy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xy_0_z[i] = 4.0 * g_xyyz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (276-279)

    #pragma omp simd aligned(g_xyyz_xz_0_x, g_xyyz_xz_0_y, g_xyyz_xz_0_z, g_xz_0_0_0_yy_xz_0_x, g_xz_0_0_0_yy_xz_0_y, g_xz_0_0_0_yy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_xz_0_x[i] = 4.0 * g_xyyz_xz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xz_0_y[i] = 4.0 * g_xyyz_xz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xz_0_z[i] = 4.0 * g_xyyz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (279-282)

    #pragma omp simd aligned(g_xyyz_yy_0_x, g_xyyz_yy_0_y, g_xyyz_yy_0_z, g_xz_0_0_0_yy_yy_0_x, g_xz_0_0_0_yy_yy_0_y, g_xz_0_0_0_yy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_yy_0_x[i] = 4.0 * g_xyyz_yy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yy_0_y[i] = 4.0 * g_xyyz_yy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yy_0_z[i] = 4.0 * g_xyyz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (282-285)

    #pragma omp simd aligned(g_xyyz_yz_0_x, g_xyyz_yz_0_y, g_xyyz_yz_0_z, g_xz_0_0_0_yy_yz_0_x, g_xz_0_0_0_yy_yz_0_y, g_xz_0_0_0_yy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_yz_0_x[i] = 4.0 * g_xyyz_yz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yz_0_y[i] = 4.0 * g_xyyz_yz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yz_0_z[i] = 4.0 * g_xyyz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (285-288)

    #pragma omp simd aligned(g_xyyz_zz_0_x, g_xyyz_zz_0_y, g_xyyz_zz_0_z, g_xz_0_0_0_yy_zz_0_x, g_xz_0_0_0_yy_zz_0_y, g_xz_0_0_0_yy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_zz_0_x[i] = 4.0 * g_xyyz_zz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_zz_0_y[i] = 4.0 * g_xyyz_zz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_zz_0_z[i] = 4.0 * g_xyyz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (288-291)

    #pragma omp simd aligned(g_xy_xx_0_x, g_xy_xx_0_y, g_xy_xx_0_z, g_xyzz_xx_0_x, g_xyzz_xx_0_y, g_xyzz_xx_0_z, g_xz_0_0_0_yz_xx_0_x, g_xz_0_0_0_yz_xx_0_y, g_xz_0_0_0_yz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_xx_0_x[i] = -2.0 * g_xy_xx_0_x[i] * a_exp + 4.0 * g_xyzz_xx_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xx_0_y[i] = -2.0 * g_xy_xx_0_y[i] * a_exp + 4.0 * g_xyzz_xx_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xx_0_z[i] = -2.0 * g_xy_xx_0_z[i] * a_exp + 4.0 * g_xyzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (291-294)

    #pragma omp simd aligned(g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z, g_xyzz_xy_0_x, g_xyzz_xy_0_y, g_xyzz_xy_0_z, g_xz_0_0_0_yz_xy_0_x, g_xz_0_0_0_yz_xy_0_y, g_xz_0_0_0_yz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_xy_0_x[i] = -2.0 * g_xy_xy_0_x[i] * a_exp + 4.0 * g_xyzz_xy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xy_0_y[i] = -2.0 * g_xy_xy_0_y[i] * a_exp + 4.0 * g_xyzz_xy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xy_0_z[i] = -2.0 * g_xy_xy_0_z[i] * a_exp + 4.0 * g_xyzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (294-297)

    #pragma omp simd aligned(g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z, g_xyzz_xz_0_x, g_xyzz_xz_0_y, g_xyzz_xz_0_z, g_xz_0_0_0_yz_xz_0_x, g_xz_0_0_0_yz_xz_0_y, g_xz_0_0_0_yz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_xz_0_x[i] = -2.0 * g_xy_xz_0_x[i] * a_exp + 4.0 * g_xyzz_xz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xz_0_y[i] = -2.0 * g_xy_xz_0_y[i] * a_exp + 4.0 * g_xyzz_xz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xz_0_z[i] = -2.0 * g_xy_xz_0_z[i] * a_exp + 4.0 * g_xyzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (297-300)

    #pragma omp simd aligned(g_xy_yy_0_x, g_xy_yy_0_y, g_xy_yy_0_z, g_xyzz_yy_0_x, g_xyzz_yy_0_y, g_xyzz_yy_0_z, g_xz_0_0_0_yz_yy_0_x, g_xz_0_0_0_yz_yy_0_y, g_xz_0_0_0_yz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_yy_0_x[i] = -2.0 * g_xy_yy_0_x[i] * a_exp + 4.0 * g_xyzz_yy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yy_0_y[i] = -2.0 * g_xy_yy_0_y[i] * a_exp + 4.0 * g_xyzz_yy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yy_0_z[i] = -2.0 * g_xy_yy_0_z[i] * a_exp + 4.0 * g_xyzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (300-303)

    #pragma omp simd aligned(g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z, g_xyzz_yz_0_x, g_xyzz_yz_0_y, g_xyzz_yz_0_z, g_xz_0_0_0_yz_yz_0_x, g_xz_0_0_0_yz_yz_0_y, g_xz_0_0_0_yz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_yz_0_x[i] = -2.0 * g_xy_yz_0_x[i] * a_exp + 4.0 * g_xyzz_yz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yz_0_y[i] = -2.0 * g_xy_yz_0_y[i] * a_exp + 4.0 * g_xyzz_yz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yz_0_z[i] = -2.0 * g_xy_yz_0_z[i] * a_exp + 4.0 * g_xyzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (303-306)

    #pragma omp simd aligned(g_xy_zz_0_x, g_xy_zz_0_y, g_xy_zz_0_z, g_xyzz_zz_0_x, g_xyzz_zz_0_y, g_xyzz_zz_0_z, g_xz_0_0_0_yz_zz_0_x, g_xz_0_0_0_yz_zz_0_y, g_xz_0_0_0_yz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_zz_0_x[i] = -2.0 * g_xy_zz_0_x[i] * a_exp + 4.0 * g_xyzz_zz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_zz_0_y[i] = -2.0 * g_xy_zz_0_y[i] * a_exp + 4.0 * g_xyzz_zz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_zz_0_z[i] = -2.0 * g_xy_zz_0_z[i] * a_exp + 4.0 * g_xyzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (306-309)

    #pragma omp simd aligned(g_xz_0_0_0_zz_xx_0_x, g_xz_0_0_0_zz_xx_0_y, g_xz_0_0_0_zz_xx_0_z, g_xz_xx_0_x, g_xz_xx_0_y, g_xz_xx_0_z, g_xzzz_xx_0_x, g_xzzz_xx_0_y, g_xzzz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_xx_0_x[i] = -4.0 * g_xz_xx_0_x[i] * a_exp + 4.0 * g_xzzz_xx_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xx_0_y[i] = -4.0 * g_xz_xx_0_y[i] * a_exp + 4.0 * g_xzzz_xx_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xx_0_z[i] = -4.0 * g_xz_xx_0_z[i] * a_exp + 4.0 * g_xzzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (309-312)

    #pragma omp simd aligned(g_xz_0_0_0_zz_xy_0_x, g_xz_0_0_0_zz_xy_0_y, g_xz_0_0_0_zz_xy_0_z, g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z, g_xzzz_xy_0_x, g_xzzz_xy_0_y, g_xzzz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_xy_0_x[i] = -4.0 * g_xz_xy_0_x[i] * a_exp + 4.0 * g_xzzz_xy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xy_0_y[i] = -4.0 * g_xz_xy_0_y[i] * a_exp + 4.0 * g_xzzz_xy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xy_0_z[i] = -4.0 * g_xz_xy_0_z[i] * a_exp + 4.0 * g_xzzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (312-315)

    #pragma omp simd aligned(g_xz_0_0_0_zz_xz_0_x, g_xz_0_0_0_zz_xz_0_y, g_xz_0_0_0_zz_xz_0_z, g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z, g_xzzz_xz_0_x, g_xzzz_xz_0_y, g_xzzz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_xz_0_x[i] = -4.0 * g_xz_xz_0_x[i] * a_exp + 4.0 * g_xzzz_xz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xz_0_y[i] = -4.0 * g_xz_xz_0_y[i] * a_exp + 4.0 * g_xzzz_xz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xz_0_z[i] = -4.0 * g_xz_xz_0_z[i] * a_exp + 4.0 * g_xzzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (315-318)

    #pragma omp simd aligned(g_xz_0_0_0_zz_yy_0_x, g_xz_0_0_0_zz_yy_0_y, g_xz_0_0_0_zz_yy_0_z, g_xz_yy_0_x, g_xz_yy_0_y, g_xz_yy_0_z, g_xzzz_yy_0_x, g_xzzz_yy_0_y, g_xzzz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_yy_0_x[i] = -4.0 * g_xz_yy_0_x[i] * a_exp + 4.0 * g_xzzz_yy_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yy_0_y[i] = -4.0 * g_xz_yy_0_y[i] * a_exp + 4.0 * g_xzzz_yy_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yy_0_z[i] = -4.0 * g_xz_yy_0_z[i] * a_exp + 4.0 * g_xzzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (318-321)

    #pragma omp simd aligned(g_xz_0_0_0_zz_yz_0_x, g_xz_0_0_0_zz_yz_0_y, g_xz_0_0_0_zz_yz_0_z, g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z, g_xzzz_yz_0_x, g_xzzz_yz_0_y, g_xzzz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_yz_0_x[i] = -4.0 * g_xz_yz_0_x[i] * a_exp + 4.0 * g_xzzz_yz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yz_0_y[i] = -4.0 * g_xz_yz_0_y[i] * a_exp + 4.0 * g_xzzz_yz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yz_0_z[i] = -4.0 * g_xz_yz_0_z[i] * a_exp + 4.0 * g_xzzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (321-324)

    #pragma omp simd aligned(g_xz_0_0_0_zz_zz_0_x, g_xz_0_0_0_zz_zz_0_y, g_xz_0_0_0_zz_zz_0_z, g_xz_zz_0_x, g_xz_zz_0_y, g_xz_zz_0_z, g_xzzz_zz_0_x, g_xzzz_zz_0_y, g_xzzz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_zz_0_x[i] = -4.0 * g_xz_zz_0_x[i] * a_exp + 4.0 * g_xzzz_zz_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_zz_0_y[i] = -4.0 * g_xz_zz_0_y[i] * a_exp + 4.0 * g_xzzz_zz_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_zz_0_z[i] = -4.0 * g_xz_zz_0_z[i] * a_exp + 4.0 * g_xzzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (324-327)

    #pragma omp simd aligned(g_xx_xx_0_x, g_xx_xx_0_y, g_xx_xx_0_z, g_xxyy_xx_0_x, g_xxyy_xx_0_y, g_xxyy_xx_0_z, g_yy_0_0_0_xx_xx_0_x, g_yy_0_0_0_xx_xx_0_y, g_yy_0_0_0_xx_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_xx_0_x[i] = -2.0 * g_xx_xx_0_x[i] * a_exp + 4.0 * g_xxyy_xx_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xx_0_y[i] = -2.0 * g_xx_xx_0_y[i] * a_exp + 4.0 * g_xxyy_xx_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xx_0_z[i] = -2.0 * g_xx_xx_0_z[i] * a_exp + 4.0 * g_xxyy_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (327-330)

    #pragma omp simd aligned(g_xx_xy_0_x, g_xx_xy_0_y, g_xx_xy_0_z, g_xxyy_xy_0_x, g_xxyy_xy_0_y, g_xxyy_xy_0_z, g_yy_0_0_0_xx_xy_0_x, g_yy_0_0_0_xx_xy_0_y, g_yy_0_0_0_xx_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_xy_0_x[i] = -2.0 * g_xx_xy_0_x[i] * a_exp + 4.0 * g_xxyy_xy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xy_0_y[i] = -2.0 * g_xx_xy_0_y[i] * a_exp + 4.0 * g_xxyy_xy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xy_0_z[i] = -2.0 * g_xx_xy_0_z[i] * a_exp + 4.0 * g_xxyy_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (330-333)

    #pragma omp simd aligned(g_xx_xz_0_x, g_xx_xz_0_y, g_xx_xz_0_z, g_xxyy_xz_0_x, g_xxyy_xz_0_y, g_xxyy_xz_0_z, g_yy_0_0_0_xx_xz_0_x, g_yy_0_0_0_xx_xz_0_y, g_yy_0_0_0_xx_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_xz_0_x[i] = -2.0 * g_xx_xz_0_x[i] * a_exp + 4.0 * g_xxyy_xz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xz_0_y[i] = -2.0 * g_xx_xz_0_y[i] * a_exp + 4.0 * g_xxyy_xz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xz_0_z[i] = -2.0 * g_xx_xz_0_z[i] * a_exp + 4.0 * g_xxyy_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (333-336)

    #pragma omp simd aligned(g_xx_yy_0_x, g_xx_yy_0_y, g_xx_yy_0_z, g_xxyy_yy_0_x, g_xxyy_yy_0_y, g_xxyy_yy_0_z, g_yy_0_0_0_xx_yy_0_x, g_yy_0_0_0_xx_yy_0_y, g_yy_0_0_0_xx_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_yy_0_x[i] = -2.0 * g_xx_yy_0_x[i] * a_exp + 4.0 * g_xxyy_yy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yy_0_y[i] = -2.0 * g_xx_yy_0_y[i] * a_exp + 4.0 * g_xxyy_yy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yy_0_z[i] = -2.0 * g_xx_yy_0_z[i] * a_exp + 4.0 * g_xxyy_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (336-339)

    #pragma omp simd aligned(g_xx_yz_0_x, g_xx_yz_0_y, g_xx_yz_0_z, g_xxyy_yz_0_x, g_xxyy_yz_0_y, g_xxyy_yz_0_z, g_yy_0_0_0_xx_yz_0_x, g_yy_0_0_0_xx_yz_0_y, g_yy_0_0_0_xx_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_yz_0_x[i] = -2.0 * g_xx_yz_0_x[i] * a_exp + 4.0 * g_xxyy_yz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yz_0_y[i] = -2.0 * g_xx_yz_0_y[i] * a_exp + 4.0 * g_xxyy_yz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yz_0_z[i] = -2.0 * g_xx_yz_0_z[i] * a_exp + 4.0 * g_xxyy_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (339-342)

    #pragma omp simd aligned(g_xx_zz_0_x, g_xx_zz_0_y, g_xx_zz_0_z, g_xxyy_zz_0_x, g_xxyy_zz_0_y, g_xxyy_zz_0_z, g_yy_0_0_0_xx_zz_0_x, g_yy_0_0_0_xx_zz_0_y, g_yy_0_0_0_xx_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_zz_0_x[i] = -2.0 * g_xx_zz_0_x[i] * a_exp + 4.0 * g_xxyy_zz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_zz_0_y[i] = -2.0 * g_xx_zz_0_y[i] * a_exp + 4.0 * g_xxyy_zz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_zz_0_z[i] = -2.0 * g_xx_zz_0_z[i] * a_exp + 4.0 * g_xxyy_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (342-345)

    #pragma omp simd aligned(g_xy_xx_0_x, g_xy_xx_0_y, g_xy_xx_0_z, g_xyyy_xx_0_x, g_xyyy_xx_0_y, g_xyyy_xx_0_z, g_yy_0_0_0_xy_xx_0_x, g_yy_0_0_0_xy_xx_0_y, g_yy_0_0_0_xy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_xx_0_x[i] = -6.0 * g_xy_xx_0_x[i] * a_exp + 4.0 * g_xyyy_xx_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xx_0_y[i] = -6.0 * g_xy_xx_0_y[i] * a_exp + 4.0 * g_xyyy_xx_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xx_0_z[i] = -6.0 * g_xy_xx_0_z[i] * a_exp + 4.0 * g_xyyy_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (345-348)

    #pragma omp simd aligned(g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z, g_xyyy_xy_0_x, g_xyyy_xy_0_y, g_xyyy_xy_0_z, g_yy_0_0_0_xy_xy_0_x, g_yy_0_0_0_xy_xy_0_y, g_yy_0_0_0_xy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_xy_0_x[i] = -6.0 * g_xy_xy_0_x[i] * a_exp + 4.0 * g_xyyy_xy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xy_0_y[i] = -6.0 * g_xy_xy_0_y[i] * a_exp + 4.0 * g_xyyy_xy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xy_0_z[i] = -6.0 * g_xy_xy_0_z[i] * a_exp + 4.0 * g_xyyy_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (348-351)

    #pragma omp simd aligned(g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z, g_xyyy_xz_0_x, g_xyyy_xz_0_y, g_xyyy_xz_0_z, g_yy_0_0_0_xy_xz_0_x, g_yy_0_0_0_xy_xz_0_y, g_yy_0_0_0_xy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_xz_0_x[i] = -6.0 * g_xy_xz_0_x[i] * a_exp + 4.0 * g_xyyy_xz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xz_0_y[i] = -6.0 * g_xy_xz_0_y[i] * a_exp + 4.0 * g_xyyy_xz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xz_0_z[i] = -6.0 * g_xy_xz_0_z[i] * a_exp + 4.0 * g_xyyy_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (351-354)

    #pragma omp simd aligned(g_xy_yy_0_x, g_xy_yy_0_y, g_xy_yy_0_z, g_xyyy_yy_0_x, g_xyyy_yy_0_y, g_xyyy_yy_0_z, g_yy_0_0_0_xy_yy_0_x, g_yy_0_0_0_xy_yy_0_y, g_yy_0_0_0_xy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_yy_0_x[i] = -6.0 * g_xy_yy_0_x[i] * a_exp + 4.0 * g_xyyy_yy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yy_0_y[i] = -6.0 * g_xy_yy_0_y[i] * a_exp + 4.0 * g_xyyy_yy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yy_0_z[i] = -6.0 * g_xy_yy_0_z[i] * a_exp + 4.0 * g_xyyy_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (354-357)

    #pragma omp simd aligned(g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z, g_xyyy_yz_0_x, g_xyyy_yz_0_y, g_xyyy_yz_0_z, g_yy_0_0_0_xy_yz_0_x, g_yy_0_0_0_xy_yz_0_y, g_yy_0_0_0_xy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_yz_0_x[i] = -6.0 * g_xy_yz_0_x[i] * a_exp + 4.0 * g_xyyy_yz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yz_0_y[i] = -6.0 * g_xy_yz_0_y[i] * a_exp + 4.0 * g_xyyy_yz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yz_0_z[i] = -6.0 * g_xy_yz_0_z[i] * a_exp + 4.0 * g_xyyy_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (357-360)

    #pragma omp simd aligned(g_xy_zz_0_x, g_xy_zz_0_y, g_xy_zz_0_z, g_xyyy_zz_0_x, g_xyyy_zz_0_y, g_xyyy_zz_0_z, g_yy_0_0_0_xy_zz_0_x, g_yy_0_0_0_xy_zz_0_y, g_yy_0_0_0_xy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_zz_0_x[i] = -6.0 * g_xy_zz_0_x[i] * a_exp + 4.0 * g_xyyy_zz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_zz_0_y[i] = -6.0 * g_xy_zz_0_y[i] * a_exp + 4.0 * g_xyyy_zz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_zz_0_z[i] = -6.0 * g_xy_zz_0_z[i] * a_exp + 4.0 * g_xyyy_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (360-363)

    #pragma omp simd aligned(g_xyyz_xx_0_x, g_xyyz_xx_0_y, g_xyyz_xx_0_z, g_xz_xx_0_x, g_xz_xx_0_y, g_xz_xx_0_z, g_yy_0_0_0_xz_xx_0_x, g_yy_0_0_0_xz_xx_0_y, g_yy_0_0_0_xz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_xx_0_x[i] = -2.0 * g_xz_xx_0_x[i] * a_exp + 4.0 * g_xyyz_xx_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xx_0_y[i] = -2.0 * g_xz_xx_0_y[i] * a_exp + 4.0 * g_xyyz_xx_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xx_0_z[i] = -2.0 * g_xz_xx_0_z[i] * a_exp + 4.0 * g_xyyz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (363-366)

    #pragma omp simd aligned(g_xyyz_xy_0_x, g_xyyz_xy_0_y, g_xyyz_xy_0_z, g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z, g_yy_0_0_0_xz_xy_0_x, g_yy_0_0_0_xz_xy_0_y, g_yy_0_0_0_xz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_xy_0_x[i] = -2.0 * g_xz_xy_0_x[i] * a_exp + 4.0 * g_xyyz_xy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xy_0_y[i] = -2.0 * g_xz_xy_0_y[i] * a_exp + 4.0 * g_xyyz_xy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xy_0_z[i] = -2.0 * g_xz_xy_0_z[i] * a_exp + 4.0 * g_xyyz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (366-369)

    #pragma omp simd aligned(g_xyyz_xz_0_x, g_xyyz_xz_0_y, g_xyyz_xz_0_z, g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z, g_yy_0_0_0_xz_xz_0_x, g_yy_0_0_0_xz_xz_0_y, g_yy_0_0_0_xz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_xz_0_x[i] = -2.0 * g_xz_xz_0_x[i] * a_exp + 4.0 * g_xyyz_xz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xz_0_y[i] = -2.0 * g_xz_xz_0_y[i] * a_exp + 4.0 * g_xyyz_xz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xz_0_z[i] = -2.0 * g_xz_xz_0_z[i] * a_exp + 4.0 * g_xyyz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (369-372)

    #pragma omp simd aligned(g_xyyz_yy_0_x, g_xyyz_yy_0_y, g_xyyz_yy_0_z, g_xz_yy_0_x, g_xz_yy_0_y, g_xz_yy_0_z, g_yy_0_0_0_xz_yy_0_x, g_yy_0_0_0_xz_yy_0_y, g_yy_0_0_0_xz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_yy_0_x[i] = -2.0 * g_xz_yy_0_x[i] * a_exp + 4.0 * g_xyyz_yy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yy_0_y[i] = -2.0 * g_xz_yy_0_y[i] * a_exp + 4.0 * g_xyyz_yy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yy_0_z[i] = -2.0 * g_xz_yy_0_z[i] * a_exp + 4.0 * g_xyyz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (372-375)

    #pragma omp simd aligned(g_xyyz_yz_0_x, g_xyyz_yz_0_y, g_xyyz_yz_0_z, g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z, g_yy_0_0_0_xz_yz_0_x, g_yy_0_0_0_xz_yz_0_y, g_yy_0_0_0_xz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_yz_0_x[i] = -2.0 * g_xz_yz_0_x[i] * a_exp + 4.0 * g_xyyz_yz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yz_0_y[i] = -2.0 * g_xz_yz_0_y[i] * a_exp + 4.0 * g_xyyz_yz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yz_0_z[i] = -2.0 * g_xz_yz_0_z[i] * a_exp + 4.0 * g_xyyz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (375-378)

    #pragma omp simd aligned(g_xyyz_zz_0_x, g_xyyz_zz_0_y, g_xyyz_zz_0_z, g_xz_zz_0_x, g_xz_zz_0_y, g_xz_zz_0_z, g_yy_0_0_0_xz_zz_0_x, g_yy_0_0_0_xz_zz_0_y, g_yy_0_0_0_xz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_zz_0_x[i] = -2.0 * g_xz_zz_0_x[i] * a_exp + 4.0 * g_xyyz_zz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_zz_0_y[i] = -2.0 * g_xz_zz_0_y[i] * a_exp + 4.0 * g_xyyz_zz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_zz_0_z[i] = -2.0 * g_xz_zz_0_z[i] * a_exp + 4.0 * g_xyyz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (378-381)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_yy_0_0_0_yy_xx_0_x, g_yy_0_0_0_yy_xx_0_y, g_yy_0_0_0_yy_xx_0_z, g_yy_xx_0_x, g_yy_xx_0_y, g_yy_xx_0_z, g_yyyy_xx_0_x, g_yyyy_xx_0_y, g_yyyy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_xx_0_x[i] = 2.0 * g_0_xx_0_x[i] - 10.0 * g_yy_xx_0_x[i] * a_exp + 4.0 * g_yyyy_xx_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xx_0_y[i] = 2.0 * g_0_xx_0_y[i] - 10.0 * g_yy_xx_0_y[i] * a_exp + 4.0 * g_yyyy_xx_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xx_0_z[i] = 2.0 * g_0_xx_0_z[i] - 10.0 * g_yy_xx_0_z[i] * a_exp + 4.0 * g_yyyy_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (381-384)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_yy_0_0_0_yy_xy_0_x, g_yy_0_0_0_yy_xy_0_y, g_yy_0_0_0_yy_xy_0_z, g_yy_xy_0_x, g_yy_xy_0_y, g_yy_xy_0_z, g_yyyy_xy_0_x, g_yyyy_xy_0_y, g_yyyy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_xy_0_x[i] = 2.0 * g_0_xy_0_x[i] - 10.0 * g_yy_xy_0_x[i] * a_exp + 4.0 * g_yyyy_xy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xy_0_y[i] = 2.0 * g_0_xy_0_y[i] - 10.0 * g_yy_xy_0_y[i] * a_exp + 4.0 * g_yyyy_xy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xy_0_z[i] = 2.0 * g_0_xy_0_z[i] - 10.0 * g_yy_xy_0_z[i] * a_exp + 4.0 * g_yyyy_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (384-387)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_yy_0_0_0_yy_xz_0_x, g_yy_0_0_0_yy_xz_0_y, g_yy_0_0_0_yy_xz_0_z, g_yy_xz_0_x, g_yy_xz_0_y, g_yy_xz_0_z, g_yyyy_xz_0_x, g_yyyy_xz_0_y, g_yyyy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_xz_0_x[i] = 2.0 * g_0_xz_0_x[i] - 10.0 * g_yy_xz_0_x[i] * a_exp + 4.0 * g_yyyy_xz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xz_0_y[i] = 2.0 * g_0_xz_0_y[i] - 10.0 * g_yy_xz_0_y[i] * a_exp + 4.0 * g_yyyy_xz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xz_0_z[i] = 2.0 * g_0_xz_0_z[i] - 10.0 * g_yy_xz_0_z[i] * a_exp + 4.0 * g_yyyy_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (387-390)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_yy_0_0_0_yy_yy_0_x, g_yy_0_0_0_yy_yy_0_y, g_yy_0_0_0_yy_yy_0_z, g_yy_yy_0_x, g_yy_yy_0_y, g_yy_yy_0_z, g_yyyy_yy_0_x, g_yyyy_yy_0_y, g_yyyy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_yy_0_x[i] = 2.0 * g_0_yy_0_x[i] - 10.0 * g_yy_yy_0_x[i] * a_exp + 4.0 * g_yyyy_yy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yy_0_y[i] = 2.0 * g_0_yy_0_y[i] - 10.0 * g_yy_yy_0_y[i] * a_exp + 4.0 * g_yyyy_yy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yy_0_z[i] = 2.0 * g_0_yy_0_z[i] - 10.0 * g_yy_yy_0_z[i] * a_exp + 4.0 * g_yyyy_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (390-393)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_yy_0_0_0_yy_yz_0_x, g_yy_0_0_0_yy_yz_0_y, g_yy_0_0_0_yy_yz_0_z, g_yy_yz_0_x, g_yy_yz_0_y, g_yy_yz_0_z, g_yyyy_yz_0_x, g_yyyy_yz_0_y, g_yyyy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_yz_0_x[i] = 2.0 * g_0_yz_0_x[i] - 10.0 * g_yy_yz_0_x[i] * a_exp + 4.0 * g_yyyy_yz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yz_0_y[i] = 2.0 * g_0_yz_0_y[i] - 10.0 * g_yy_yz_0_y[i] * a_exp + 4.0 * g_yyyy_yz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yz_0_z[i] = 2.0 * g_0_yz_0_z[i] - 10.0 * g_yy_yz_0_z[i] * a_exp + 4.0 * g_yyyy_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (393-396)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_yy_0_0_0_yy_zz_0_x, g_yy_0_0_0_yy_zz_0_y, g_yy_0_0_0_yy_zz_0_z, g_yy_zz_0_x, g_yy_zz_0_y, g_yy_zz_0_z, g_yyyy_zz_0_x, g_yyyy_zz_0_y, g_yyyy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_zz_0_x[i] = 2.0 * g_0_zz_0_x[i] - 10.0 * g_yy_zz_0_x[i] * a_exp + 4.0 * g_yyyy_zz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_zz_0_y[i] = 2.0 * g_0_zz_0_y[i] - 10.0 * g_yy_zz_0_y[i] * a_exp + 4.0 * g_yyyy_zz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_zz_0_z[i] = 2.0 * g_0_zz_0_z[i] - 10.0 * g_yy_zz_0_z[i] * a_exp + 4.0 * g_yyyy_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (396-399)

    #pragma omp simd aligned(g_yy_0_0_0_yz_xx_0_x, g_yy_0_0_0_yz_xx_0_y, g_yy_0_0_0_yz_xx_0_z, g_yyyz_xx_0_x, g_yyyz_xx_0_y, g_yyyz_xx_0_z, g_yz_xx_0_x, g_yz_xx_0_y, g_yz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_xx_0_x[i] = -6.0 * g_yz_xx_0_x[i] * a_exp + 4.0 * g_yyyz_xx_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xx_0_y[i] = -6.0 * g_yz_xx_0_y[i] * a_exp + 4.0 * g_yyyz_xx_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xx_0_z[i] = -6.0 * g_yz_xx_0_z[i] * a_exp + 4.0 * g_yyyz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (399-402)

    #pragma omp simd aligned(g_yy_0_0_0_yz_xy_0_x, g_yy_0_0_0_yz_xy_0_y, g_yy_0_0_0_yz_xy_0_z, g_yyyz_xy_0_x, g_yyyz_xy_0_y, g_yyyz_xy_0_z, g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_xy_0_x[i] = -6.0 * g_yz_xy_0_x[i] * a_exp + 4.0 * g_yyyz_xy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xy_0_y[i] = -6.0 * g_yz_xy_0_y[i] * a_exp + 4.0 * g_yyyz_xy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xy_0_z[i] = -6.0 * g_yz_xy_0_z[i] * a_exp + 4.0 * g_yyyz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (402-405)

    #pragma omp simd aligned(g_yy_0_0_0_yz_xz_0_x, g_yy_0_0_0_yz_xz_0_y, g_yy_0_0_0_yz_xz_0_z, g_yyyz_xz_0_x, g_yyyz_xz_0_y, g_yyyz_xz_0_z, g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_xz_0_x[i] = -6.0 * g_yz_xz_0_x[i] * a_exp + 4.0 * g_yyyz_xz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xz_0_y[i] = -6.0 * g_yz_xz_0_y[i] * a_exp + 4.0 * g_yyyz_xz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xz_0_z[i] = -6.0 * g_yz_xz_0_z[i] * a_exp + 4.0 * g_yyyz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (405-408)

    #pragma omp simd aligned(g_yy_0_0_0_yz_yy_0_x, g_yy_0_0_0_yz_yy_0_y, g_yy_0_0_0_yz_yy_0_z, g_yyyz_yy_0_x, g_yyyz_yy_0_y, g_yyyz_yy_0_z, g_yz_yy_0_x, g_yz_yy_0_y, g_yz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_yy_0_x[i] = -6.0 * g_yz_yy_0_x[i] * a_exp + 4.0 * g_yyyz_yy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yy_0_y[i] = -6.0 * g_yz_yy_0_y[i] * a_exp + 4.0 * g_yyyz_yy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yy_0_z[i] = -6.0 * g_yz_yy_0_z[i] * a_exp + 4.0 * g_yyyz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (408-411)

    #pragma omp simd aligned(g_yy_0_0_0_yz_yz_0_x, g_yy_0_0_0_yz_yz_0_y, g_yy_0_0_0_yz_yz_0_z, g_yyyz_yz_0_x, g_yyyz_yz_0_y, g_yyyz_yz_0_z, g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_yz_0_x[i] = -6.0 * g_yz_yz_0_x[i] * a_exp + 4.0 * g_yyyz_yz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yz_0_y[i] = -6.0 * g_yz_yz_0_y[i] * a_exp + 4.0 * g_yyyz_yz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yz_0_z[i] = -6.0 * g_yz_yz_0_z[i] * a_exp + 4.0 * g_yyyz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (411-414)

    #pragma omp simd aligned(g_yy_0_0_0_yz_zz_0_x, g_yy_0_0_0_yz_zz_0_y, g_yy_0_0_0_yz_zz_0_z, g_yyyz_zz_0_x, g_yyyz_zz_0_y, g_yyyz_zz_0_z, g_yz_zz_0_x, g_yz_zz_0_y, g_yz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_zz_0_x[i] = -6.0 * g_yz_zz_0_x[i] * a_exp + 4.0 * g_yyyz_zz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_zz_0_y[i] = -6.0 * g_yz_zz_0_y[i] * a_exp + 4.0 * g_yyyz_zz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_zz_0_z[i] = -6.0 * g_yz_zz_0_z[i] * a_exp + 4.0 * g_yyyz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (414-417)

    #pragma omp simd aligned(g_yy_0_0_0_zz_xx_0_x, g_yy_0_0_0_zz_xx_0_y, g_yy_0_0_0_zz_xx_0_z, g_yyzz_xx_0_x, g_yyzz_xx_0_y, g_yyzz_xx_0_z, g_zz_xx_0_x, g_zz_xx_0_y, g_zz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_xx_0_x[i] = -2.0 * g_zz_xx_0_x[i] * a_exp + 4.0 * g_yyzz_xx_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xx_0_y[i] = -2.0 * g_zz_xx_0_y[i] * a_exp + 4.0 * g_yyzz_xx_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xx_0_z[i] = -2.0 * g_zz_xx_0_z[i] * a_exp + 4.0 * g_yyzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (417-420)

    #pragma omp simd aligned(g_yy_0_0_0_zz_xy_0_x, g_yy_0_0_0_zz_xy_0_y, g_yy_0_0_0_zz_xy_0_z, g_yyzz_xy_0_x, g_yyzz_xy_0_y, g_yyzz_xy_0_z, g_zz_xy_0_x, g_zz_xy_0_y, g_zz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_xy_0_x[i] = -2.0 * g_zz_xy_0_x[i] * a_exp + 4.0 * g_yyzz_xy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xy_0_y[i] = -2.0 * g_zz_xy_0_y[i] * a_exp + 4.0 * g_yyzz_xy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xy_0_z[i] = -2.0 * g_zz_xy_0_z[i] * a_exp + 4.0 * g_yyzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (420-423)

    #pragma omp simd aligned(g_yy_0_0_0_zz_xz_0_x, g_yy_0_0_0_zz_xz_0_y, g_yy_0_0_0_zz_xz_0_z, g_yyzz_xz_0_x, g_yyzz_xz_0_y, g_yyzz_xz_0_z, g_zz_xz_0_x, g_zz_xz_0_y, g_zz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_xz_0_x[i] = -2.0 * g_zz_xz_0_x[i] * a_exp + 4.0 * g_yyzz_xz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xz_0_y[i] = -2.0 * g_zz_xz_0_y[i] * a_exp + 4.0 * g_yyzz_xz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xz_0_z[i] = -2.0 * g_zz_xz_0_z[i] * a_exp + 4.0 * g_yyzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (423-426)

    #pragma omp simd aligned(g_yy_0_0_0_zz_yy_0_x, g_yy_0_0_0_zz_yy_0_y, g_yy_0_0_0_zz_yy_0_z, g_yyzz_yy_0_x, g_yyzz_yy_0_y, g_yyzz_yy_0_z, g_zz_yy_0_x, g_zz_yy_0_y, g_zz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_yy_0_x[i] = -2.0 * g_zz_yy_0_x[i] * a_exp + 4.0 * g_yyzz_yy_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yy_0_y[i] = -2.0 * g_zz_yy_0_y[i] * a_exp + 4.0 * g_yyzz_yy_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yy_0_z[i] = -2.0 * g_zz_yy_0_z[i] * a_exp + 4.0 * g_yyzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (426-429)

    #pragma omp simd aligned(g_yy_0_0_0_zz_yz_0_x, g_yy_0_0_0_zz_yz_0_y, g_yy_0_0_0_zz_yz_0_z, g_yyzz_yz_0_x, g_yyzz_yz_0_y, g_yyzz_yz_0_z, g_zz_yz_0_x, g_zz_yz_0_y, g_zz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_yz_0_x[i] = -2.0 * g_zz_yz_0_x[i] * a_exp + 4.0 * g_yyzz_yz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yz_0_y[i] = -2.0 * g_zz_yz_0_y[i] * a_exp + 4.0 * g_yyzz_yz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yz_0_z[i] = -2.0 * g_zz_yz_0_z[i] * a_exp + 4.0 * g_yyzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (429-432)

    #pragma omp simd aligned(g_yy_0_0_0_zz_zz_0_x, g_yy_0_0_0_zz_zz_0_y, g_yy_0_0_0_zz_zz_0_z, g_yyzz_zz_0_x, g_yyzz_zz_0_y, g_yyzz_zz_0_z, g_zz_zz_0_x, g_zz_zz_0_y, g_zz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_zz_0_x[i] = -2.0 * g_zz_zz_0_x[i] * a_exp + 4.0 * g_yyzz_zz_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_zz_0_y[i] = -2.0 * g_zz_zz_0_y[i] * a_exp + 4.0 * g_yyzz_zz_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_zz_0_z[i] = -2.0 * g_zz_zz_0_z[i] * a_exp + 4.0 * g_yyzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (432-435)

    #pragma omp simd aligned(g_xxyz_xx_0_x, g_xxyz_xx_0_y, g_xxyz_xx_0_z, g_yz_0_0_0_xx_xx_0_x, g_yz_0_0_0_xx_xx_0_y, g_yz_0_0_0_xx_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_xx_0_x[i] = 4.0 * g_xxyz_xx_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xx_0_y[i] = 4.0 * g_xxyz_xx_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xx_0_z[i] = 4.0 * g_xxyz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (435-438)

    #pragma omp simd aligned(g_xxyz_xy_0_x, g_xxyz_xy_0_y, g_xxyz_xy_0_z, g_yz_0_0_0_xx_xy_0_x, g_yz_0_0_0_xx_xy_0_y, g_yz_0_0_0_xx_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_xy_0_x[i] = 4.0 * g_xxyz_xy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xy_0_y[i] = 4.0 * g_xxyz_xy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xy_0_z[i] = 4.0 * g_xxyz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (438-441)

    #pragma omp simd aligned(g_xxyz_xz_0_x, g_xxyz_xz_0_y, g_xxyz_xz_0_z, g_yz_0_0_0_xx_xz_0_x, g_yz_0_0_0_xx_xz_0_y, g_yz_0_0_0_xx_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_xz_0_x[i] = 4.0 * g_xxyz_xz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xz_0_y[i] = 4.0 * g_xxyz_xz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xz_0_z[i] = 4.0 * g_xxyz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (441-444)

    #pragma omp simd aligned(g_xxyz_yy_0_x, g_xxyz_yy_0_y, g_xxyz_yy_0_z, g_yz_0_0_0_xx_yy_0_x, g_yz_0_0_0_xx_yy_0_y, g_yz_0_0_0_xx_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_yy_0_x[i] = 4.0 * g_xxyz_yy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yy_0_y[i] = 4.0 * g_xxyz_yy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yy_0_z[i] = 4.0 * g_xxyz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (444-447)

    #pragma omp simd aligned(g_xxyz_yz_0_x, g_xxyz_yz_0_y, g_xxyz_yz_0_z, g_yz_0_0_0_xx_yz_0_x, g_yz_0_0_0_xx_yz_0_y, g_yz_0_0_0_xx_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_yz_0_x[i] = 4.0 * g_xxyz_yz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yz_0_y[i] = 4.0 * g_xxyz_yz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yz_0_z[i] = 4.0 * g_xxyz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (447-450)

    #pragma omp simd aligned(g_xxyz_zz_0_x, g_xxyz_zz_0_y, g_xxyz_zz_0_z, g_yz_0_0_0_xx_zz_0_x, g_yz_0_0_0_xx_zz_0_y, g_yz_0_0_0_xx_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_zz_0_x[i] = 4.0 * g_xxyz_zz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_zz_0_y[i] = 4.0 * g_xxyz_zz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_zz_0_z[i] = 4.0 * g_xxyz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (450-453)

    #pragma omp simd aligned(g_xyyz_xx_0_x, g_xyyz_xx_0_y, g_xyyz_xx_0_z, g_xz_xx_0_x, g_xz_xx_0_y, g_xz_xx_0_z, g_yz_0_0_0_xy_xx_0_x, g_yz_0_0_0_xy_xx_0_y, g_yz_0_0_0_xy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_xx_0_x[i] = -2.0 * g_xz_xx_0_x[i] * a_exp + 4.0 * g_xyyz_xx_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xx_0_y[i] = -2.0 * g_xz_xx_0_y[i] * a_exp + 4.0 * g_xyyz_xx_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xx_0_z[i] = -2.0 * g_xz_xx_0_z[i] * a_exp + 4.0 * g_xyyz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (453-456)

    #pragma omp simd aligned(g_xyyz_xy_0_x, g_xyyz_xy_0_y, g_xyyz_xy_0_z, g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z, g_yz_0_0_0_xy_xy_0_x, g_yz_0_0_0_xy_xy_0_y, g_yz_0_0_0_xy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_xy_0_x[i] = -2.0 * g_xz_xy_0_x[i] * a_exp + 4.0 * g_xyyz_xy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xy_0_y[i] = -2.0 * g_xz_xy_0_y[i] * a_exp + 4.0 * g_xyyz_xy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xy_0_z[i] = -2.0 * g_xz_xy_0_z[i] * a_exp + 4.0 * g_xyyz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (456-459)

    #pragma omp simd aligned(g_xyyz_xz_0_x, g_xyyz_xz_0_y, g_xyyz_xz_0_z, g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z, g_yz_0_0_0_xy_xz_0_x, g_yz_0_0_0_xy_xz_0_y, g_yz_0_0_0_xy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_xz_0_x[i] = -2.0 * g_xz_xz_0_x[i] * a_exp + 4.0 * g_xyyz_xz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xz_0_y[i] = -2.0 * g_xz_xz_0_y[i] * a_exp + 4.0 * g_xyyz_xz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xz_0_z[i] = -2.0 * g_xz_xz_0_z[i] * a_exp + 4.0 * g_xyyz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (459-462)

    #pragma omp simd aligned(g_xyyz_yy_0_x, g_xyyz_yy_0_y, g_xyyz_yy_0_z, g_xz_yy_0_x, g_xz_yy_0_y, g_xz_yy_0_z, g_yz_0_0_0_xy_yy_0_x, g_yz_0_0_0_xy_yy_0_y, g_yz_0_0_0_xy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_yy_0_x[i] = -2.0 * g_xz_yy_0_x[i] * a_exp + 4.0 * g_xyyz_yy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yy_0_y[i] = -2.0 * g_xz_yy_0_y[i] * a_exp + 4.0 * g_xyyz_yy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yy_0_z[i] = -2.0 * g_xz_yy_0_z[i] * a_exp + 4.0 * g_xyyz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (462-465)

    #pragma omp simd aligned(g_xyyz_yz_0_x, g_xyyz_yz_0_y, g_xyyz_yz_0_z, g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z, g_yz_0_0_0_xy_yz_0_x, g_yz_0_0_0_xy_yz_0_y, g_yz_0_0_0_xy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_yz_0_x[i] = -2.0 * g_xz_yz_0_x[i] * a_exp + 4.0 * g_xyyz_yz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yz_0_y[i] = -2.0 * g_xz_yz_0_y[i] * a_exp + 4.0 * g_xyyz_yz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yz_0_z[i] = -2.0 * g_xz_yz_0_z[i] * a_exp + 4.0 * g_xyyz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (465-468)

    #pragma omp simd aligned(g_xyyz_zz_0_x, g_xyyz_zz_0_y, g_xyyz_zz_0_z, g_xz_zz_0_x, g_xz_zz_0_y, g_xz_zz_0_z, g_yz_0_0_0_xy_zz_0_x, g_yz_0_0_0_xy_zz_0_y, g_yz_0_0_0_xy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_zz_0_x[i] = -2.0 * g_xz_zz_0_x[i] * a_exp + 4.0 * g_xyyz_zz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_zz_0_y[i] = -2.0 * g_xz_zz_0_y[i] * a_exp + 4.0 * g_xyyz_zz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_zz_0_z[i] = -2.0 * g_xz_zz_0_z[i] * a_exp + 4.0 * g_xyyz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (468-471)

    #pragma omp simd aligned(g_xy_xx_0_x, g_xy_xx_0_y, g_xy_xx_0_z, g_xyzz_xx_0_x, g_xyzz_xx_0_y, g_xyzz_xx_0_z, g_yz_0_0_0_xz_xx_0_x, g_yz_0_0_0_xz_xx_0_y, g_yz_0_0_0_xz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_xx_0_x[i] = -2.0 * g_xy_xx_0_x[i] * a_exp + 4.0 * g_xyzz_xx_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xx_0_y[i] = -2.0 * g_xy_xx_0_y[i] * a_exp + 4.0 * g_xyzz_xx_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xx_0_z[i] = -2.0 * g_xy_xx_0_z[i] * a_exp + 4.0 * g_xyzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (471-474)

    #pragma omp simd aligned(g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z, g_xyzz_xy_0_x, g_xyzz_xy_0_y, g_xyzz_xy_0_z, g_yz_0_0_0_xz_xy_0_x, g_yz_0_0_0_xz_xy_0_y, g_yz_0_0_0_xz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_xy_0_x[i] = -2.0 * g_xy_xy_0_x[i] * a_exp + 4.0 * g_xyzz_xy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xy_0_y[i] = -2.0 * g_xy_xy_0_y[i] * a_exp + 4.0 * g_xyzz_xy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xy_0_z[i] = -2.0 * g_xy_xy_0_z[i] * a_exp + 4.0 * g_xyzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (474-477)

    #pragma omp simd aligned(g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z, g_xyzz_xz_0_x, g_xyzz_xz_0_y, g_xyzz_xz_0_z, g_yz_0_0_0_xz_xz_0_x, g_yz_0_0_0_xz_xz_0_y, g_yz_0_0_0_xz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_xz_0_x[i] = -2.0 * g_xy_xz_0_x[i] * a_exp + 4.0 * g_xyzz_xz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xz_0_y[i] = -2.0 * g_xy_xz_0_y[i] * a_exp + 4.0 * g_xyzz_xz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xz_0_z[i] = -2.0 * g_xy_xz_0_z[i] * a_exp + 4.0 * g_xyzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (477-480)

    #pragma omp simd aligned(g_xy_yy_0_x, g_xy_yy_0_y, g_xy_yy_0_z, g_xyzz_yy_0_x, g_xyzz_yy_0_y, g_xyzz_yy_0_z, g_yz_0_0_0_xz_yy_0_x, g_yz_0_0_0_xz_yy_0_y, g_yz_0_0_0_xz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_yy_0_x[i] = -2.0 * g_xy_yy_0_x[i] * a_exp + 4.0 * g_xyzz_yy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yy_0_y[i] = -2.0 * g_xy_yy_0_y[i] * a_exp + 4.0 * g_xyzz_yy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yy_0_z[i] = -2.0 * g_xy_yy_0_z[i] * a_exp + 4.0 * g_xyzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (480-483)

    #pragma omp simd aligned(g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z, g_xyzz_yz_0_x, g_xyzz_yz_0_y, g_xyzz_yz_0_z, g_yz_0_0_0_xz_yz_0_x, g_yz_0_0_0_xz_yz_0_y, g_yz_0_0_0_xz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_yz_0_x[i] = -2.0 * g_xy_yz_0_x[i] * a_exp + 4.0 * g_xyzz_yz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yz_0_y[i] = -2.0 * g_xy_yz_0_y[i] * a_exp + 4.0 * g_xyzz_yz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yz_0_z[i] = -2.0 * g_xy_yz_0_z[i] * a_exp + 4.0 * g_xyzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (483-486)

    #pragma omp simd aligned(g_xy_zz_0_x, g_xy_zz_0_y, g_xy_zz_0_z, g_xyzz_zz_0_x, g_xyzz_zz_0_y, g_xyzz_zz_0_z, g_yz_0_0_0_xz_zz_0_x, g_yz_0_0_0_xz_zz_0_y, g_yz_0_0_0_xz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_zz_0_x[i] = -2.0 * g_xy_zz_0_x[i] * a_exp + 4.0 * g_xyzz_zz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_zz_0_y[i] = -2.0 * g_xy_zz_0_y[i] * a_exp + 4.0 * g_xyzz_zz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_zz_0_z[i] = -2.0 * g_xy_zz_0_z[i] * a_exp + 4.0 * g_xyzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (486-489)

    #pragma omp simd aligned(g_yyyz_xx_0_x, g_yyyz_xx_0_y, g_yyyz_xx_0_z, g_yz_0_0_0_yy_xx_0_x, g_yz_0_0_0_yy_xx_0_y, g_yz_0_0_0_yy_xx_0_z, g_yz_xx_0_x, g_yz_xx_0_y, g_yz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_xx_0_x[i] = -4.0 * g_yz_xx_0_x[i] * a_exp + 4.0 * g_yyyz_xx_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xx_0_y[i] = -4.0 * g_yz_xx_0_y[i] * a_exp + 4.0 * g_yyyz_xx_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xx_0_z[i] = -4.0 * g_yz_xx_0_z[i] * a_exp + 4.0 * g_yyyz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (489-492)

    #pragma omp simd aligned(g_yyyz_xy_0_x, g_yyyz_xy_0_y, g_yyyz_xy_0_z, g_yz_0_0_0_yy_xy_0_x, g_yz_0_0_0_yy_xy_0_y, g_yz_0_0_0_yy_xy_0_z, g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_xy_0_x[i] = -4.0 * g_yz_xy_0_x[i] * a_exp + 4.0 * g_yyyz_xy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xy_0_y[i] = -4.0 * g_yz_xy_0_y[i] * a_exp + 4.0 * g_yyyz_xy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xy_0_z[i] = -4.0 * g_yz_xy_0_z[i] * a_exp + 4.0 * g_yyyz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (492-495)

    #pragma omp simd aligned(g_yyyz_xz_0_x, g_yyyz_xz_0_y, g_yyyz_xz_0_z, g_yz_0_0_0_yy_xz_0_x, g_yz_0_0_0_yy_xz_0_y, g_yz_0_0_0_yy_xz_0_z, g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_xz_0_x[i] = -4.0 * g_yz_xz_0_x[i] * a_exp + 4.0 * g_yyyz_xz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xz_0_y[i] = -4.0 * g_yz_xz_0_y[i] * a_exp + 4.0 * g_yyyz_xz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xz_0_z[i] = -4.0 * g_yz_xz_0_z[i] * a_exp + 4.0 * g_yyyz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (495-498)

    #pragma omp simd aligned(g_yyyz_yy_0_x, g_yyyz_yy_0_y, g_yyyz_yy_0_z, g_yz_0_0_0_yy_yy_0_x, g_yz_0_0_0_yy_yy_0_y, g_yz_0_0_0_yy_yy_0_z, g_yz_yy_0_x, g_yz_yy_0_y, g_yz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_yy_0_x[i] = -4.0 * g_yz_yy_0_x[i] * a_exp + 4.0 * g_yyyz_yy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yy_0_y[i] = -4.0 * g_yz_yy_0_y[i] * a_exp + 4.0 * g_yyyz_yy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yy_0_z[i] = -4.0 * g_yz_yy_0_z[i] * a_exp + 4.0 * g_yyyz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (498-501)

    #pragma omp simd aligned(g_yyyz_yz_0_x, g_yyyz_yz_0_y, g_yyyz_yz_0_z, g_yz_0_0_0_yy_yz_0_x, g_yz_0_0_0_yy_yz_0_y, g_yz_0_0_0_yy_yz_0_z, g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_yz_0_x[i] = -4.0 * g_yz_yz_0_x[i] * a_exp + 4.0 * g_yyyz_yz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yz_0_y[i] = -4.0 * g_yz_yz_0_y[i] * a_exp + 4.0 * g_yyyz_yz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yz_0_z[i] = -4.0 * g_yz_yz_0_z[i] * a_exp + 4.0 * g_yyyz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (501-504)

    #pragma omp simd aligned(g_yyyz_zz_0_x, g_yyyz_zz_0_y, g_yyyz_zz_0_z, g_yz_0_0_0_yy_zz_0_x, g_yz_0_0_0_yy_zz_0_y, g_yz_0_0_0_yy_zz_0_z, g_yz_zz_0_x, g_yz_zz_0_y, g_yz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_zz_0_x[i] = -4.0 * g_yz_zz_0_x[i] * a_exp + 4.0 * g_yyyz_zz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_zz_0_y[i] = -4.0 * g_yz_zz_0_y[i] * a_exp + 4.0 * g_yyyz_zz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_zz_0_z[i] = -4.0 * g_yz_zz_0_z[i] * a_exp + 4.0 * g_yyyz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (504-507)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_yy_xx_0_x, g_yy_xx_0_y, g_yy_xx_0_z, g_yyzz_xx_0_x, g_yyzz_xx_0_y, g_yyzz_xx_0_z, g_yz_0_0_0_yz_xx_0_x, g_yz_0_0_0_yz_xx_0_y, g_yz_0_0_0_yz_xx_0_z, g_zz_xx_0_x, g_zz_xx_0_y, g_zz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_xx_0_x[i] = g_0_xx_0_x[i] - 2.0 * g_zz_xx_0_x[i] * a_exp - 2.0 * g_yy_xx_0_x[i] * a_exp + 4.0 * g_yyzz_xx_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xx_0_y[i] = g_0_xx_0_y[i] - 2.0 * g_zz_xx_0_y[i] * a_exp - 2.0 * g_yy_xx_0_y[i] * a_exp + 4.0 * g_yyzz_xx_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xx_0_z[i] = g_0_xx_0_z[i] - 2.0 * g_zz_xx_0_z[i] * a_exp - 2.0 * g_yy_xx_0_z[i] * a_exp + 4.0 * g_yyzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (507-510)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_yy_xy_0_x, g_yy_xy_0_y, g_yy_xy_0_z, g_yyzz_xy_0_x, g_yyzz_xy_0_y, g_yyzz_xy_0_z, g_yz_0_0_0_yz_xy_0_x, g_yz_0_0_0_yz_xy_0_y, g_yz_0_0_0_yz_xy_0_z, g_zz_xy_0_x, g_zz_xy_0_y, g_zz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_xy_0_x[i] = g_0_xy_0_x[i] - 2.0 * g_zz_xy_0_x[i] * a_exp - 2.0 * g_yy_xy_0_x[i] * a_exp + 4.0 * g_yyzz_xy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xy_0_y[i] = g_0_xy_0_y[i] - 2.0 * g_zz_xy_0_y[i] * a_exp - 2.0 * g_yy_xy_0_y[i] * a_exp + 4.0 * g_yyzz_xy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xy_0_z[i] = g_0_xy_0_z[i] - 2.0 * g_zz_xy_0_z[i] * a_exp - 2.0 * g_yy_xy_0_z[i] * a_exp + 4.0 * g_yyzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (510-513)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_yy_xz_0_x, g_yy_xz_0_y, g_yy_xz_0_z, g_yyzz_xz_0_x, g_yyzz_xz_0_y, g_yyzz_xz_0_z, g_yz_0_0_0_yz_xz_0_x, g_yz_0_0_0_yz_xz_0_y, g_yz_0_0_0_yz_xz_0_z, g_zz_xz_0_x, g_zz_xz_0_y, g_zz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_xz_0_x[i] = g_0_xz_0_x[i] - 2.0 * g_zz_xz_0_x[i] * a_exp - 2.0 * g_yy_xz_0_x[i] * a_exp + 4.0 * g_yyzz_xz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xz_0_y[i] = g_0_xz_0_y[i] - 2.0 * g_zz_xz_0_y[i] * a_exp - 2.0 * g_yy_xz_0_y[i] * a_exp + 4.0 * g_yyzz_xz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xz_0_z[i] = g_0_xz_0_z[i] - 2.0 * g_zz_xz_0_z[i] * a_exp - 2.0 * g_yy_xz_0_z[i] * a_exp + 4.0 * g_yyzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (513-516)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_yy_yy_0_x, g_yy_yy_0_y, g_yy_yy_0_z, g_yyzz_yy_0_x, g_yyzz_yy_0_y, g_yyzz_yy_0_z, g_yz_0_0_0_yz_yy_0_x, g_yz_0_0_0_yz_yy_0_y, g_yz_0_0_0_yz_yy_0_z, g_zz_yy_0_x, g_zz_yy_0_y, g_zz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_yy_0_x[i] = g_0_yy_0_x[i] - 2.0 * g_zz_yy_0_x[i] * a_exp - 2.0 * g_yy_yy_0_x[i] * a_exp + 4.0 * g_yyzz_yy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yy_0_y[i] = g_0_yy_0_y[i] - 2.0 * g_zz_yy_0_y[i] * a_exp - 2.0 * g_yy_yy_0_y[i] * a_exp + 4.0 * g_yyzz_yy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yy_0_z[i] = g_0_yy_0_z[i] - 2.0 * g_zz_yy_0_z[i] * a_exp - 2.0 * g_yy_yy_0_z[i] * a_exp + 4.0 * g_yyzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (516-519)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_yy_yz_0_x, g_yy_yz_0_y, g_yy_yz_0_z, g_yyzz_yz_0_x, g_yyzz_yz_0_y, g_yyzz_yz_0_z, g_yz_0_0_0_yz_yz_0_x, g_yz_0_0_0_yz_yz_0_y, g_yz_0_0_0_yz_yz_0_z, g_zz_yz_0_x, g_zz_yz_0_y, g_zz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_yz_0_x[i] = g_0_yz_0_x[i] - 2.0 * g_zz_yz_0_x[i] * a_exp - 2.0 * g_yy_yz_0_x[i] * a_exp + 4.0 * g_yyzz_yz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yz_0_y[i] = g_0_yz_0_y[i] - 2.0 * g_zz_yz_0_y[i] * a_exp - 2.0 * g_yy_yz_0_y[i] * a_exp + 4.0 * g_yyzz_yz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yz_0_z[i] = g_0_yz_0_z[i] - 2.0 * g_zz_yz_0_z[i] * a_exp - 2.0 * g_yy_yz_0_z[i] * a_exp + 4.0 * g_yyzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (519-522)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_yy_zz_0_x, g_yy_zz_0_y, g_yy_zz_0_z, g_yyzz_zz_0_x, g_yyzz_zz_0_y, g_yyzz_zz_0_z, g_yz_0_0_0_yz_zz_0_x, g_yz_0_0_0_yz_zz_0_y, g_yz_0_0_0_yz_zz_0_z, g_zz_zz_0_x, g_zz_zz_0_y, g_zz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_zz_0_x[i] = g_0_zz_0_x[i] - 2.0 * g_zz_zz_0_x[i] * a_exp - 2.0 * g_yy_zz_0_x[i] * a_exp + 4.0 * g_yyzz_zz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_zz_0_y[i] = g_0_zz_0_y[i] - 2.0 * g_zz_zz_0_y[i] * a_exp - 2.0 * g_yy_zz_0_y[i] * a_exp + 4.0 * g_yyzz_zz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_zz_0_z[i] = g_0_zz_0_z[i] - 2.0 * g_zz_zz_0_z[i] * a_exp - 2.0 * g_yy_zz_0_z[i] * a_exp + 4.0 * g_yyzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (522-525)

    #pragma omp simd aligned(g_yz_0_0_0_zz_xx_0_x, g_yz_0_0_0_zz_xx_0_y, g_yz_0_0_0_zz_xx_0_z, g_yz_xx_0_x, g_yz_xx_0_y, g_yz_xx_0_z, g_yzzz_xx_0_x, g_yzzz_xx_0_y, g_yzzz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_xx_0_x[i] = -4.0 * g_yz_xx_0_x[i] * a_exp + 4.0 * g_yzzz_xx_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xx_0_y[i] = -4.0 * g_yz_xx_0_y[i] * a_exp + 4.0 * g_yzzz_xx_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xx_0_z[i] = -4.0 * g_yz_xx_0_z[i] * a_exp + 4.0 * g_yzzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (525-528)

    #pragma omp simd aligned(g_yz_0_0_0_zz_xy_0_x, g_yz_0_0_0_zz_xy_0_y, g_yz_0_0_0_zz_xy_0_z, g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z, g_yzzz_xy_0_x, g_yzzz_xy_0_y, g_yzzz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_xy_0_x[i] = -4.0 * g_yz_xy_0_x[i] * a_exp + 4.0 * g_yzzz_xy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xy_0_y[i] = -4.0 * g_yz_xy_0_y[i] * a_exp + 4.0 * g_yzzz_xy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xy_0_z[i] = -4.0 * g_yz_xy_0_z[i] * a_exp + 4.0 * g_yzzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (528-531)

    #pragma omp simd aligned(g_yz_0_0_0_zz_xz_0_x, g_yz_0_0_0_zz_xz_0_y, g_yz_0_0_0_zz_xz_0_z, g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z, g_yzzz_xz_0_x, g_yzzz_xz_0_y, g_yzzz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_xz_0_x[i] = -4.0 * g_yz_xz_0_x[i] * a_exp + 4.0 * g_yzzz_xz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xz_0_y[i] = -4.0 * g_yz_xz_0_y[i] * a_exp + 4.0 * g_yzzz_xz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xz_0_z[i] = -4.0 * g_yz_xz_0_z[i] * a_exp + 4.0 * g_yzzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (531-534)

    #pragma omp simd aligned(g_yz_0_0_0_zz_yy_0_x, g_yz_0_0_0_zz_yy_0_y, g_yz_0_0_0_zz_yy_0_z, g_yz_yy_0_x, g_yz_yy_0_y, g_yz_yy_0_z, g_yzzz_yy_0_x, g_yzzz_yy_0_y, g_yzzz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_yy_0_x[i] = -4.0 * g_yz_yy_0_x[i] * a_exp + 4.0 * g_yzzz_yy_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yy_0_y[i] = -4.0 * g_yz_yy_0_y[i] * a_exp + 4.0 * g_yzzz_yy_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yy_0_z[i] = -4.0 * g_yz_yy_0_z[i] * a_exp + 4.0 * g_yzzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (534-537)

    #pragma omp simd aligned(g_yz_0_0_0_zz_yz_0_x, g_yz_0_0_0_zz_yz_0_y, g_yz_0_0_0_zz_yz_0_z, g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z, g_yzzz_yz_0_x, g_yzzz_yz_0_y, g_yzzz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_yz_0_x[i] = -4.0 * g_yz_yz_0_x[i] * a_exp + 4.0 * g_yzzz_yz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yz_0_y[i] = -4.0 * g_yz_yz_0_y[i] * a_exp + 4.0 * g_yzzz_yz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yz_0_z[i] = -4.0 * g_yz_yz_0_z[i] * a_exp + 4.0 * g_yzzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (537-540)

    #pragma omp simd aligned(g_yz_0_0_0_zz_zz_0_x, g_yz_0_0_0_zz_zz_0_y, g_yz_0_0_0_zz_zz_0_z, g_yz_zz_0_x, g_yz_zz_0_y, g_yz_zz_0_z, g_yzzz_zz_0_x, g_yzzz_zz_0_y, g_yzzz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_zz_0_x[i] = -4.0 * g_yz_zz_0_x[i] * a_exp + 4.0 * g_yzzz_zz_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_zz_0_y[i] = -4.0 * g_yz_zz_0_y[i] * a_exp + 4.0 * g_yzzz_zz_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_zz_0_z[i] = -4.0 * g_yz_zz_0_z[i] * a_exp + 4.0 * g_yzzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (540-543)

    #pragma omp simd aligned(g_xx_xx_0_x, g_xx_xx_0_y, g_xx_xx_0_z, g_xxzz_xx_0_x, g_xxzz_xx_0_y, g_xxzz_xx_0_z, g_zz_0_0_0_xx_xx_0_x, g_zz_0_0_0_xx_xx_0_y, g_zz_0_0_0_xx_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_xx_0_x[i] = -2.0 * g_xx_xx_0_x[i] * a_exp + 4.0 * g_xxzz_xx_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xx_0_y[i] = -2.0 * g_xx_xx_0_y[i] * a_exp + 4.0 * g_xxzz_xx_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xx_0_z[i] = -2.0 * g_xx_xx_0_z[i] * a_exp + 4.0 * g_xxzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (543-546)

    #pragma omp simd aligned(g_xx_xy_0_x, g_xx_xy_0_y, g_xx_xy_0_z, g_xxzz_xy_0_x, g_xxzz_xy_0_y, g_xxzz_xy_0_z, g_zz_0_0_0_xx_xy_0_x, g_zz_0_0_0_xx_xy_0_y, g_zz_0_0_0_xx_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_xy_0_x[i] = -2.0 * g_xx_xy_0_x[i] * a_exp + 4.0 * g_xxzz_xy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xy_0_y[i] = -2.0 * g_xx_xy_0_y[i] * a_exp + 4.0 * g_xxzz_xy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xy_0_z[i] = -2.0 * g_xx_xy_0_z[i] * a_exp + 4.0 * g_xxzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (546-549)

    #pragma omp simd aligned(g_xx_xz_0_x, g_xx_xz_0_y, g_xx_xz_0_z, g_xxzz_xz_0_x, g_xxzz_xz_0_y, g_xxzz_xz_0_z, g_zz_0_0_0_xx_xz_0_x, g_zz_0_0_0_xx_xz_0_y, g_zz_0_0_0_xx_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_xz_0_x[i] = -2.0 * g_xx_xz_0_x[i] * a_exp + 4.0 * g_xxzz_xz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xz_0_y[i] = -2.0 * g_xx_xz_0_y[i] * a_exp + 4.0 * g_xxzz_xz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xz_0_z[i] = -2.0 * g_xx_xz_0_z[i] * a_exp + 4.0 * g_xxzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (549-552)

    #pragma omp simd aligned(g_xx_yy_0_x, g_xx_yy_0_y, g_xx_yy_0_z, g_xxzz_yy_0_x, g_xxzz_yy_0_y, g_xxzz_yy_0_z, g_zz_0_0_0_xx_yy_0_x, g_zz_0_0_0_xx_yy_0_y, g_zz_0_0_0_xx_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_yy_0_x[i] = -2.0 * g_xx_yy_0_x[i] * a_exp + 4.0 * g_xxzz_yy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yy_0_y[i] = -2.0 * g_xx_yy_0_y[i] * a_exp + 4.0 * g_xxzz_yy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yy_0_z[i] = -2.0 * g_xx_yy_0_z[i] * a_exp + 4.0 * g_xxzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (552-555)

    #pragma omp simd aligned(g_xx_yz_0_x, g_xx_yz_0_y, g_xx_yz_0_z, g_xxzz_yz_0_x, g_xxzz_yz_0_y, g_xxzz_yz_0_z, g_zz_0_0_0_xx_yz_0_x, g_zz_0_0_0_xx_yz_0_y, g_zz_0_0_0_xx_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_yz_0_x[i] = -2.0 * g_xx_yz_0_x[i] * a_exp + 4.0 * g_xxzz_yz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yz_0_y[i] = -2.0 * g_xx_yz_0_y[i] * a_exp + 4.0 * g_xxzz_yz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yz_0_z[i] = -2.0 * g_xx_yz_0_z[i] * a_exp + 4.0 * g_xxzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (555-558)

    #pragma omp simd aligned(g_xx_zz_0_x, g_xx_zz_0_y, g_xx_zz_0_z, g_xxzz_zz_0_x, g_xxzz_zz_0_y, g_xxzz_zz_0_z, g_zz_0_0_0_xx_zz_0_x, g_zz_0_0_0_xx_zz_0_y, g_zz_0_0_0_xx_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_zz_0_x[i] = -2.0 * g_xx_zz_0_x[i] * a_exp + 4.0 * g_xxzz_zz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_zz_0_y[i] = -2.0 * g_xx_zz_0_y[i] * a_exp + 4.0 * g_xxzz_zz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_zz_0_z[i] = -2.0 * g_xx_zz_0_z[i] * a_exp + 4.0 * g_xxzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (558-561)

    #pragma omp simd aligned(g_xy_xx_0_x, g_xy_xx_0_y, g_xy_xx_0_z, g_xyzz_xx_0_x, g_xyzz_xx_0_y, g_xyzz_xx_0_z, g_zz_0_0_0_xy_xx_0_x, g_zz_0_0_0_xy_xx_0_y, g_zz_0_0_0_xy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_xx_0_x[i] = -2.0 * g_xy_xx_0_x[i] * a_exp + 4.0 * g_xyzz_xx_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xx_0_y[i] = -2.0 * g_xy_xx_0_y[i] * a_exp + 4.0 * g_xyzz_xx_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xx_0_z[i] = -2.0 * g_xy_xx_0_z[i] * a_exp + 4.0 * g_xyzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (561-564)

    #pragma omp simd aligned(g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z, g_xyzz_xy_0_x, g_xyzz_xy_0_y, g_xyzz_xy_0_z, g_zz_0_0_0_xy_xy_0_x, g_zz_0_0_0_xy_xy_0_y, g_zz_0_0_0_xy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_xy_0_x[i] = -2.0 * g_xy_xy_0_x[i] * a_exp + 4.0 * g_xyzz_xy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xy_0_y[i] = -2.0 * g_xy_xy_0_y[i] * a_exp + 4.0 * g_xyzz_xy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xy_0_z[i] = -2.0 * g_xy_xy_0_z[i] * a_exp + 4.0 * g_xyzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (564-567)

    #pragma omp simd aligned(g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z, g_xyzz_xz_0_x, g_xyzz_xz_0_y, g_xyzz_xz_0_z, g_zz_0_0_0_xy_xz_0_x, g_zz_0_0_0_xy_xz_0_y, g_zz_0_0_0_xy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_xz_0_x[i] = -2.0 * g_xy_xz_0_x[i] * a_exp + 4.0 * g_xyzz_xz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xz_0_y[i] = -2.0 * g_xy_xz_0_y[i] * a_exp + 4.0 * g_xyzz_xz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xz_0_z[i] = -2.0 * g_xy_xz_0_z[i] * a_exp + 4.0 * g_xyzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (567-570)

    #pragma omp simd aligned(g_xy_yy_0_x, g_xy_yy_0_y, g_xy_yy_0_z, g_xyzz_yy_0_x, g_xyzz_yy_0_y, g_xyzz_yy_0_z, g_zz_0_0_0_xy_yy_0_x, g_zz_0_0_0_xy_yy_0_y, g_zz_0_0_0_xy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_yy_0_x[i] = -2.0 * g_xy_yy_0_x[i] * a_exp + 4.0 * g_xyzz_yy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yy_0_y[i] = -2.0 * g_xy_yy_0_y[i] * a_exp + 4.0 * g_xyzz_yy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yy_0_z[i] = -2.0 * g_xy_yy_0_z[i] * a_exp + 4.0 * g_xyzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (570-573)

    #pragma omp simd aligned(g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z, g_xyzz_yz_0_x, g_xyzz_yz_0_y, g_xyzz_yz_0_z, g_zz_0_0_0_xy_yz_0_x, g_zz_0_0_0_xy_yz_0_y, g_zz_0_0_0_xy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_yz_0_x[i] = -2.0 * g_xy_yz_0_x[i] * a_exp + 4.0 * g_xyzz_yz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yz_0_y[i] = -2.0 * g_xy_yz_0_y[i] * a_exp + 4.0 * g_xyzz_yz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yz_0_z[i] = -2.0 * g_xy_yz_0_z[i] * a_exp + 4.0 * g_xyzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (573-576)

    #pragma omp simd aligned(g_xy_zz_0_x, g_xy_zz_0_y, g_xy_zz_0_z, g_xyzz_zz_0_x, g_xyzz_zz_0_y, g_xyzz_zz_0_z, g_zz_0_0_0_xy_zz_0_x, g_zz_0_0_0_xy_zz_0_y, g_zz_0_0_0_xy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_zz_0_x[i] = -2.0 * g_xy_zz_0_x[i] * a_exp + 4.0 * g_xyzz_zz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_zz_0_y[i] = -2.0 * g_xy_zz_0_y[i] * a_exp + 4.0 * g_xyzz_zz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_zz_0_z[i] = -2.0 * g_xy_zz_0_z[i] * a_exp + 4.0 * g_xyzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (576-579)

    #pragma omp simd aligned(g_xz_xx_0_x, g_xz_xx_0_y, g_xz_xx_0_z, g_xzzz_xx_0_x, g_xzzz_xx_0_y, g_xzzz_xx_0_z, g_zz_0_0_0_xz_xx_0_x, g_zz_0_0_0_xz_xx_0_y, g_zz_0_0_0_xz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_xx_0_x[i] = -6.0 * g_xz_xx_0_x[i] * a_exp + 4.0 * g_xzzz_xx_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xx_0_y[i] = -6.0 * g_xz_xx_0_y[i] * a_exp + 4.0 * g_xzzz_xx_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xx_0_z[i] = -6.0 * g_xz_xx_0_z[i] * a_exp + 4.0 * g_xzzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (579-582)

    #pragma omp simd aligned(g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z, g_xzzz_xy_0_x, g_xzzz_xy_0_y, g_xzzz_xy_0_z, g_zz_0_0_0_xz_xy_0_x, g_zz_0_0_0_xz_xy_0_y, g_zz_0_0_0_xz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_xy_0_x[i] = -6.0 * g_xz_xy_0_x[i] * a_exp + 4.0 * g_xzzz_xy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xy_0_y[i] = -6.0 * g_xz_xy_0_y[i] * a_exp + 4.0 * g_xzzz_xy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xy_0_z[i] = -6.0 * g_xz_xy_0_z[i] * a_exp + 4.0 * g_xzzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (582-585)

    #pragma omp simd aligned(g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z, g_xzzz_xz_0_x, g_xzzz_xz_0_y, g_xzzz_xz_0_z, g_zz_0_0_0_xz_xz_0_x, g_zz_0_0_0_xz_xz_0_y, g_zz_0_0_0_xz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_xz_0_x[i] = -6.0 * g_xz_xz_0_x[i] * a_exp + 4.0 * g_xzzz_xz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xz_0_y[i] = -6.0 * g_xz_xz_0_y[i] * a_exp + 4.0 * g_xzzz_xz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xz_0_z[i] = -6.0 * g_xz_xz_0_z[i] * a_exp + 4.0 * g_xzzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (585-588)

    #pragma omp simd aligned(g_xz_yy_0_x, g_xz_yy_0_y, g_xz_yy_0_z, g_xzzz_yy_0_x, g_xzzz_yy_0_y, g_xzzz_yy_0_z, g_zz_0_0_0_xz_yy_0_x, g_zz_0_0_0_xz_yy_0_y, g_zz_0_0_0_xz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_yy_0_x[i] = -6.0 * g_xz_yy_0_x[i] * a_exp + 4.0 * g_xzzz_yy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yy_0_y[i] = -6.0 * g_xz_yy_0_y[i] * a_exp + 4.0 * g_xzzz_yy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yy_0_z[i] = -6.0 * g_xz_yy_0_z[i] * a_exp + 4.0 * g_xzzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (588-591)

    #pragma omp simd aligned(g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z, g_xzzz_yz_0_x, g_xzzz_yz_0_y, g_xzzz_yz_0_z, g_zz_0_0_0_xz_yz_0_x, g_zz_0_0_0_xz_yz_0_y, g_zz_0_0_0_xz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_yz_0_x[i] = -6.0 * g_xz_yz_0_x[i] * a_exp + 4.0 * g_xzzz_yz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yz_0_y[i] = -6.0 * g_xz_yz_0_y[i] * a_exp + 4.0 * g_xzzz_yz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yz_0_z[i] = -6.0 * g_xz_yz_0_z[i] * a_exp + 4.0 * g_xzzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (591-594)

    #pragma omp simd aligned(g_xz_zz_0_x, g_xz_zz_0_y, g_xz_zz_0_z, g_xzzz_zz_0_x, g_xzzz_zz_0_y, g_xzzz_zz_0_z, g_zz_0_0_0_xz_zz_0_x, g_zz_0_0_0_xz_zz_0_y, g_zz_0_0_0_xz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_zz_0_x[i] = -6.0 * g_xz_zz_0_x[i] * a_exp + 4.0 * g_xzzz_zz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_zz_0_y[i] = -6.0 * g_xz_zz_0_y[i] * a_exp + 4.0 * g_xzzz_zz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_zz_0_z[i] = -6.0 * g_xz_zz_0_z[i] * a_exp + 4.0 * g_xzzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (594-597)

    #pragma omp simd aligned(g_yy_xx_0_x, g_yy_xx_0_y, g_yy_xx_0_z, g_yyzz_xx_0_x, g_yyzz_xx_0_y, g_yyzz_xx_0_z, g_zz_0_0_0_yy_xx_0_x, g_zz_0_0_0_yy_xx_0_y, g_zz_0_0_0_yy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_xx_0_x[i] = -2.0 * g_yy_xx_0_x[i] * a_exp + 4.0 * g_yyzz_xx_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xx_0_y[i] = -2.0 * g_yy_xx_0_y[i] * a_exp + 4.0 * g_yyzz_xx_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xx_0_z[i] = -2.0 * g_yy_xx_0_z[i] * a_exp + 4.0 * g_yyzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (597-600)

    #pragma omp simd aligned(g_yy_xy_0_x, g_yy_xy_0_y, g_yy_xy_0_z, g_yyzz_xy_0_x, g_yyzz_xy_0_y, g_yyzz_xy_0_z, g_zz_0_0_0_yy_xy_0_x, g_zz_0_0_0_yy_xy_0_y, g_zz_0_0_0_yy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_xy_0_x[i] = -2.0 * g_yy_xy_0_x[i] * a_exp + 4.0 * g_yyzz_xy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xy_0_y[i] = -2.0 * g_yy_xy_0_y[i] * a_exp + 4.0 * g_yyzz_xy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xy_0_z[i] = -2.0 * g_yy_xy_0_z[i] * a_exp + 4.0 * g_yyzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (600-603)

    #pragma omp simd aligned(g_yy_xz_0_x, g_yy_xz_0_y, g_yy_xz_0_z, g_yyzz_xz_0_x, g_yyzz_xz_0_y, g_yyzz_xz_0_z, g_zz_0_0_0_yy_xz_0_x, g_zz_0_0_0_yy_xz_0_y, g_zz_0_0_0_yy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_xz_0_x[i] = -2.0 * g_yy_xz_0_x[i] * a_exp + 4.0 * g_yyzz_xz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xz_0_y[i] = -2.0 * g_yy_xz_0_y[i] * a_exp + 4.0 * g_yyzz_xz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xz_0_z[i] = -2.0 * g_yy_xz_0_z[i] * a_exp + 4.0 * g_yyzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (603-606)

    #pragma omp simd aligned(g_yy_yy_0_x, g_yy_yy_0_y, g_yy_yy_0_z, g_yyzz_yy_0_x, g_yyzz_yy_0_y, g_yyzz_yy_0_z, g_zz_0_0_0_yy_yy_0_x, g_zz_0_0_0_yy_yy_0_y, g_zz_0_0_0_yy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_yy_0_x[i] = -2.0 * g_yy_yy_0_x[i] * a_exp + 4.0 * g_yyzz_yy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yy_0_y[i] = -2.0 * g_yy_yy_0_y[i] * a_exp + 4.0 * g_yyzz_yy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yy_0_z[i] = -2.0 * g_yy_yy_0_z[i] * a_exp + 4.0 * g_yyzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (606-609)

    #pragma omp simd aligned(g_yy_yz_0_x, g_yy_yz_0_y, g_yy_yz_0_z, g_yyzz_yz_0_x, g_yyzz_yz_0_y, g_yyzz_yz_0_z, g_zz_0_0_0_yy_yz_0_x, g_zz_0_0_0_yy_yz_0_y, g_zz_0_0_0_yy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_yz_0_x[i] = -2.0 * g_yy_yz_0_x[i] * a_exp + 4.0 * g_yyzz_yz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yz_0_y[i] = -2.0 * g_yy_yz_0_y[i] * a_exp + 4.0 * g_yyzz_yz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yz_0_z[i] = -2.0 * g_yy_yz_0_z[i] * a_exp + 4.0 * g_yyzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (609-612)

    #pragma omp simd aligned(g_yy_zz_0_x, g_yy_zz_0_y, g_yy_zz_0_z, g_yyzz_zz_0_x, g_yyzz_zz_0_y, g_yyzz_zz_0_z, g_zz_0_0_0_yy_zz_0_x, g_zz_0_0_0_yy_zz_0_y, g_zz_0_0_0_yy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_zz_0_x[i] = -2.0 * g_yy_zz_0_x[i] * a_exp + 4.0 * g_yyzz_zz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_zz_0_y[i] = -2.0 * g_yy_zz_0_y[i] * a_exp + 4.0 * g_yyzz_zz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_zz_0_z[i] = -2.0 * g_yy_zz_0_z[i] * a_exp + 4.0 * g_yyzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (612-615)

    #pragma omp simd aligned(g_yz_xx_0_x, g_yz_xx_0_y, g_yz_xx_0_z, g_yzzz_xx_0_x, g_yzzz_xx_0_y, g_yzzz_xx_0_z, g_zz_0_0_0_yz_xx_0_x, g_zz_0_0_0_yz_xx_0_y, g_zz_0_0_0_yz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_xx_0_x[i] = -6.0 * g_yz_xx_0_x[i] * a_exp + 4.0 * g_yzzz_xx_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xx_0_y[i] = -6.0 * g_yz_xx_0_y[i] * a_exp + 4.0 * g_yzzz_xx_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xx_0_z[i] = -6.0 * g_yz_xx_0_z[i] * a_exp + 4.0 * g_yzzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (615-618)

    #pragma omp simd aligned(g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z, g_yzzz_xy_0_x, g_yzzz_xy_0_y, g_yzzz_xy_0_z, g_zz_0_0_0_yz_xy_0_x, g_zz_0_0_0_yz_xy_0_y, g_zz_0_0_0_yz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_xy_0_x[i] = -6.0 * g_yz_xy_0_x[i] * a_exp + 4.0 * g_yzzz_xy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xy_0_y[i] = -6.0 * g_yz_xy_0_y[i] * a_exp + 4.0 * g_yzzz_xy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xy_0_z[i] = -6.0 * g_yz_xy_0_z[i] * a_exp + 4.0 * g_yzzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (618-621)

    #pragma omp simd aligned(g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z, g_yzzz_xz_0_x, g_yzzz_xz_0_y, g_yzzz_xz_0_z, g_zz_0_0_0_yz_xz_0_x, g_zz_0_0_0_yz_xz_0_y, g_zz_0_0_0_yz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_xz_0_x[i] = -6.0 * g_yz_xz_0_x[i] * a_exp + 4.0 * g_yzzz_xz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xz_0_y[i] = -6.0 * g_yz_xz_0_y[i] * a_exp + 4.0 * g_yzzz_xz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xz_0_z[i] = -6.0 * g_yz_xz_0_z[i] * a_exp + 4.0 * g_yzzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (621-624)

    #pragma omp simd aligned(g_yz_yy_0_x, g_yz_yy_0_y, g_yz_yy_0_z, g_yzzz_yy_0_x, g_yzzz_yy_0_y, g_yzzz_yy_0_z, g_zz_0_0_0_yz_yy_0_x, g_zz_0_0_0_yz_yy_0_y, g_zz_0_0_0_yz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_yy_0_x[i] = -6.0 * g_yz_yy_0_x[i] * a_exp + 4.0 * g_yzzz_yy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yy_0_y[i] = -6.0 * g_yz_yy_0_y[i] * a_exp + 4.0 * g_yzzz_yy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yy_0_z[i] = -6.0 * g_yz_yy_0_z[i] * a_exp + 4.0 * g_yzzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (624-627)

    #pragma omp simd aligned(g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z, g_yzzz_yz_0_x, g_yzzz_yz_0_y, g_yzzz_yz_0_z, g_zz_0_0_0_yz_yz_0_x, g_zz_0_0_0_yz_yz_0_y, g_zz_0_0_0_yz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_yz_0_x[i] = -6.0 * g_yz_yz_0_x[i] * a_exp + 4.0 * g_yzzz_yz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yz_0_y[i] = -6.0 * g_yz_yz_0_y[i] * a_exp + 4.0 * g_yzzz_yz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yz_0_z[i] = -6.0 * g_yz_yz_0_z[i] * a_exp + 4.0 * g_yzzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (627-630)

    #pragma omp simd aligned(g_yz_zz_0_x, g_yz_zz_0_y, g_yz_zz_0_z, g_yzzz_zz_0_x, g_yzzz_zz_0_y, g_yzzz_zz_0_z, g_zz_0_0_0_yz_zz_0_x, g_zz_0_0_0_yz_zz_0_y, g_zz_0_0_0_yz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_zz_0_x[i] = -6.0 * g_yz_zz_0_x[i] * a_exp + 4.0 * g_yzzz_zz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_zz_0_y[i] = -6.0 * g_yz_zz_0_y[i] * a_exp + 4.0 * g_yzzz_zz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_zz_0_z[i] = -6.0 * g_yz_zz_0_z[i] * a_exp + 4.0 * g_yzzz_zz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (630-633)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_zz_0_0_0_zz_xx_0_x, g_zz_0_0_0_zz_xx_0_y, g_zz_0_0_0_zz_xx_0_z, g_zz_xx_0_x, g_zz_xx_0_y, g_zz_xx_0_z, g_zzzz_xx_0_x, g_zzzz_xx_0_y, g_zzzz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_xx_0_x[i] = 2.0 * g_0_xx_0_x[i] - 10.0 * g_zz_xx_0_x[i] * a_exp + 4.0 * g_zzzz_xx_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xx_0_y[i] = 2.0 * g_0_xx_0_y[i] - 10.0 * g_zz_xx_0_y[i] * a_exp + 4.0 * g_zzzz_xx_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xx_0_z[i] = 2.0 * g_0_xx_0_z[i] - 10.0 * g_zz_xx_0_z[i] * a_exp + 4.0 * g_zzzz_xx_0_z[i] * a_exp * a_exp;
    }
    // integrals block (633-636)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_zz_0_0_0_zz_xy_0_x, g_zz_0_0_0_zz_xy_0_y, g_zz_0_0_0_zz_xy_0_z, g_zz_xy_0_x, g_zz_xy_0_y, g_zz_xy_0_z, g_zzzz_xy_0_x, g_zzzz_xy_0_y, g_zzzz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_xy_0_x[i] = 2.0 * g_0_xy_0_x[i] - 10.0 * g_zz_xy_0_x[i] * a_exp + 4.0 * g_zzzz_xy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xy_0_y[i] = 2.0 * g_0_xy_0_y[i] - 10.0 * g_zz_xy_0_y[i] * a_exp + 4.0 * g_zzzz_xy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xy_0_z[i] = 2.0 * g_0_xy_0_z[i] - 10.0 * g_zz_xy_0_z[i] * a_exp + 4.0 * g_zzzz_xy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (636-639)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_zz_0_0_0_zz_xz_0_x, g_zz_0_0_0_zz_xz_0_y, g_zz_0_0_0_zz_xz_0_z, g_zz_xz_0_x, g_zz_xz_0_y, g_zz_xz_0_z, g_zzzz_xz_0_x, g_zzzz_xz_0_y, g_zzzz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_xz_0_x[i] = 2.0 * g_0_xz_0_x[i] - 10.0 * g_zz_xz_0_x[i] * a_exp + 4.0 * g_zzzz_xz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xz_0_y[i] = 2.0 * g_0_xz_0_y[i] - 10.0 * g_zz_xz_0_y[i] * a_exp + 4.0 * g_zzzz_xz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xz_0_z[i] = 2.0 * g_0_xz_0_z[i] - 10.0 * g_zz_xz_0_z[i] * a_exp + 4.0 * g_zzzz_xz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (639-642)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_zz_0_0_0_zz_yy_0_x, g_zz_0_0_0_zz_yy_0_y, g_zz_0_0_0_zz_yy_0_z, g_zz_yy_0_x, g_zz_yy_0_y, g_zz_yy_0_z, g_zzzz_yy_0_x, g_zzzz_yy_0_y, g_zzzz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_yy_0_x[i] = 2.0 * g_0_yy_0_x[i] - 10.0 * g_zz_yy_0_x[i] * a_exp + 4.0 * g_zzzz_yy_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yy_0_y[i] = 2.0 * g_0_yy_0_y[i] - 10.0 * g_zz_yy_0_y[i] * a_exp + 4.0 * g_zzzz_yy_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yy_0_z[i] = 2.0 * g_0_yy_0_z[i] - 10.0 * g_zz_yy_0_z[i] * a_exp + 4.0 * g_zzzz_yy_0_z[i] * a_exp * a_exp;
    }
    // integrals block (642-645)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_zz_0_0_0_zz_yz_0_x, g_zz_0_0_0_zz_yz_0_y, g_zz_0_0_0_zz_yz_0_z, g_zz_yz_0_x, g_zz_yz_0_y, g_zz_yz_0_z, g_zzzz_yz_0_x, g_zzzz_yz_0_y, g_zzzz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_yz_0_x[i] = 2.0 * g_0_yz_0_x[i] - 10.0 * g_zz_yz_0_x[i] * a_exp + 4.0 * g_zzzz_yz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yz_0_y[i] = 2.0 * g_0_yz_0_y[i] - 10.0 * g_zz_yz_0_y[i] * a_exp + 4.0 * g_zzzz_yz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yz_0_z[i] = 2.0 * g_0_yz_0_z[i] - 10.0 * g_zz_yz_0_z[i] * a_exp + 4.0 * g_zzzz_yz_0_z[i] * a_exp * a_exp;
    }
    // integrals block (645-648)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_zz_0_0_0_zz_zz_0_x, g_zz_0_0_0_zz_zz_0_y, g_zz_0_0_0_zz_zz_0_z, g_zz_zz_0_x, g_zz_zz_0_y, g_zz_zz_0_z, g_zzzz_zz_0_x, g_zzzz_zz_0_y, g_zzzz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_zz_0_x[i] = 2.0 * g_0_zz_0_x[i] - 10.0 * g_zz_zz_0_x[i] * a_exp + 4.0 * g_zzzz_zz_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_zz_0_y[i] = 2.0 * g_0_zz_0_y[i] - 10.0 * g_zz_zz_0_y[i] * a_exp + 4.0 * g_zzzz_zz_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_zz_0_z[i] = 2.0 * g_0_zz_0_z[i] - 10.0 * g_zz_zz_0_z[i] * a_exp + 4.0 * g_zzzz_zz_0_z[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

