#include "GeomDeriv1100OfScalarForPPDD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_ppdd_0(CSimdArray<double>& buffer_1100_ppdd,
                     const CSimdArray<double>& buffer_ssdd,
                     const CSimdArray<double>& buffer_sddd,
                     const CSimdArray<double>& buffer_dsdd,
                     const CSimdArray<double>& buffer_dddd,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_ppdd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_ssdd

    auto g_0_0_xx_xx = buffer_ssdd[0];

    auto g_0_0_xx_xy = buffer_ssdd[1];

    auto g_0_0_xx_xz = buffer_ssdd[2];

    auto g_0_0_xx_yy = buffer_ssdd[3];

    auto g_0_0_xx_yz = buffer_ssdd[4];

    auto g_0_0_xx_zz = buffer_ssdd[5];

    auto g_0_0_xy_xx = buffer_ssdd[6];

    auto g_0_0_xy_xy = buffer_ssdd[7];

    auto g_0_0_xy_xz = buffer_ssdd[8];

    auto g_0_0_xy_yy = buffer_ssdd[9];

    auto g_0_0_xy_yz = buffer_ssdd[10];

    auto g_0_0_xy_zz = buffer_ssdd[11];

    auto g_0_0_xz_xx = buffer_ssdd[12];

    auto g_0_0_xz_xy = buffer_ssdd[13];

    auto g_0_0_xz_xz = buffer_ssdd[14];

    auto g_0_0_xz_yy = buffer_ssdd[15];

    auto g_0_0_xz_yz = buffer_ssdd[16];

    auto g_0_0_xz_zz = buffer_ssdd[17];

    auto g_0_0_yy_xx = buffer_ssdd[18];

    auto g_0_0_yy_xy = buffer_ssdd[19];

    auto g_0_0_yy_xz = buffer_ssdd[20];

    auto g_0_0_yy_yy = buffer_ssdd[21];

    auto g_0_0_yy_yz = buffer_ssdd[22];

    auto g_0_0_yy_zz = buffer_ssdd[23];

    auto g_0_0_yz_xx = buffer_ssdd[24];

    auto g_0_0_yz_xy = buffer_ssdd[25];

    auto g_0_0_yz_xz = buffer_ssdd[26];

    auto g_0_0_yz_yy = buffer_ssdd[27];

    auto g_0_0_yz_yz = buffer_ssdd[28];

    auto g_0_0_yz_zz = buffer_ssdd[29];

    auto g_0_0_zz_xx = buffer_ssdd[30];

    auto g_0_0_zz_xy = buffer_ssdd[31];

    auto g_0_0_zz_xz = buffer_ssdd[32];

    auto g_0_0_zz_yy = buffer_ssdd[33];

    auto g_0_0_zz_yz = buffer_ssdd[34];

    auto g_0_0_zz_zz = buffer_ssdd[35];

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

    /// Set up components of auxilary buffer : buffer_dsdd

    auto g_xx_0_xx_xx = buffer_dsdd[0];

    auto g_xx_0_xx_xy = buffer_dsdd[1];

    auto g_xx_0_xx_xz = buffer_dsdd[2];

    auto g_xx_0_xx_yy = buffer_dsdd[3];

    auto g_xx_0_xx_yz = buffer_dsdd[4];

    auto g_xx_0_xx_zz = buffer_dsdd[5];

    auto g_xx_0_xy_xx = buffer_dsdd[6];

    auto g_xx_0_xy_xy = buffer_dsdd[7];

    auto g_xx_0_xy_xz = buffer_dsdd[8];

    auto g_xx_0_xy_yy = buffer_dsdd[9];

    auto g_xx_0_xy_yz = buffer_dsdd[10];

    auto g_xx_0_xy_zz = buffer_dsdd[11];

    auto g_xx_0_xz_xx = buffer_dsdd[12];

    auto g_xx_0_xz_xy = buffer_dsdd[13];

    auto g_xx_0_xz_xz = buffer_dsdd[14];

    auto g_xx_0_xz_yy = buffer_dsdd[15];

    auto g_xx_0_xz_yz = buffer_dsdd[16];

    auto g_xx_0_xz_zz = buffer_dsdd[17];

    auto g_xx_0_yy_xx = buffer_dsdd[18];

    auto g_xx_0_yy_xy = buffer_dsdd[19];

    auto g_xx_0_yy_xz = buffer_dsdd[20];

    auto g_xx_0_yy_yy = buffer_dsdd[21];

    auto g_xx_0_yy_yz = buffer_dsdd[22];

    auto g_xx_0_yy_zz = buffer_dsdd[23];

    auto g_xx_0_yz_xx = buffer_dsdd[24];

    auto g_xx_0_yz_xy = buffer_dsdd[25];

    auto g_xx_0_yz_xz = buffer_dsdd[26];

    auto g_xx_0_yz_yy = buffer_dsdd[27];

    auto g_xx_0_yz_yz = buffer_dsdd[28];

    auto g_xx_0_yz_zz = buffer_dsdd[29];

    auto g_xx_0_zz_xx = buffer_dsdd[30];

    auto g_xx_0_zz_xy = buffer_dsdd[31];

    auto g_xx_0_zz_xz = buffer_dsdd[32];

    auto g_xx_0_zz_yy = buffer_dsdd[33];

    auto g_xx_0_zz_yz = buffer_dsdd[34];

    auto g_xx_0_zz_zz = buffer_dsdd[35];

    auto g_xy_0_xx_xx = buffer_dsdd[36];

    auto g_xy_0_xx_xy = buffer_dsdd[37];

    auto g_xy_0_xx_xz = buffer_dsdd[38];

    auto g_xy_0_xx_yy = buffer_dsdd[39];

    auto g_xy_0_xx_yz = buffer_dsdd[40];

    auto g_xy_0_xx_zz = buffer_dsdd[41];

    auto g_xy_0_xy_xx = buffer_dsdd[42];

    auto g_xy_0_xy_xy = buffer_dsdd[43];

    auto g_xy_0_xy_xz = buffer_dsdd[44];

    auto g_xy_0_xy_yy = buffer_dsdd[45];

    auto g_xy_0_xy_yz = buffer_dsdd[46];

    auto g_xy_0_xy_zz = buffer_dsdd[47];

    auto g_xy_0_xz_xx = buffer_dsdd[48];

    auto g_xy_0_xz_xy = buffer_dsdd[49];

    auto g_xy_0_xz_xz = buffer_dsdd[50];

    auto g_xy_0_xz_yy = buffer_dsdd[51];

    auto g_xy_0_xz_yz = buffer_dsdd[52];

    auto g_xy_0_xz_zz = buffer_dsdd[53];

    auto g_xy_0_yy_xx = buffer_dsdd[54];

    auto g_xy_0_yy_xy = buffer_dsdd[55];

    auto g_xy_0_yy_xz = buffer_dsdd[56];

    auto g_xy_0_yy_yy = buffer_dsdd[57];

    auto g_xy_0_yy_yz = buffer_dsdd[58];

    auto g_xy_0_yy_zz = buffer_dsdd[59];

    auto g_xy_0_yz_xx = buffer_dsdd[60];

    auto g_xy_0_yz_xy = buffer_dsdd[61];

    auto g_xy_0_yz_xz = buffer_dsdd[62];

    auto g_xy_0_yz_yy = buffer_dsdd[63];

    auto g_xy_0_yz_yz = buffer_dsdd[64];

    auto g_xy_0_yz_zz = buffer_dsdd[65];

    auto g_xy_0_zz_xx = buffer_dsdd[66];

    auto g_xy_0_zz_xy = buffer_dsdd[67];

    auto g_xy_0_zz_xz = buffer_dsdd[68];

    auto g_xy_0_zz_yy = buffer_dsdd[69];

    auto g_xy_0_zz_yz = buffer_dsdd[70];

    auto g_xy_0_zz_zz = buffer_dsdd[71];

    auto g_xz_0_xx_xx = buffer_dsdd[72];

    auto g_xz_0_xx_xy = buffer_dsdd[73];

    auto g_xz_0_xx_xz = buffer_dsdd[74];

    auto g_xz_0_xx_yy = buffer_dsdd[75];

    auto g_xz_0_xx_yz = buffer_dsdd[76];

    auto g_xz_0_xx_zz = buffer_dsdd[77];

    auto g_xz_0_xy_xx = buffer_dsdd[78];

    auto g_xz_0_xy_xy = buffer_dsdd[79];

    auto g_xz_0_xy_xz = buffer_dsdd[80];

    auto g_xz_0_xy_yy = buffer_dsdd[81];

    auto g_xz_0_xy_yz = buffer_dsdd[82];

    auto g_xz_0_xy_zz = buffer_dsdd[83];

    auto g_xz_0_xz_xx = buffer_dsdd[84];

    auto g_xz_0_xz_xy = buffer_dsdd[85];

    auto g_xz_0_xz_xz = buffer_dsdd[86];

    auto g_xz_0_xz_yy = buffer_dsdd[87];

    auto g_xz_0_xz_yz = buffer_dsdd[88];

    auto g_xz_0_xz_zz = buffer_dsdd[89];

    auto g_xz_0_yy_xx = buffer_dsdd[90];

    auto g_xz_0_yy_xy = buffer_dsdd[91];

    auto g_xz_0_yy_xz = buffer_dsdd[92];

    auto g_xz_0_yy_yy = buffer_dsdd[93];

    auto g_xz_0_yy_yz = buffer_dsdd[94];

    auto g_xz_0_yy_zz = buffer_dsdd[95];

    auto g_xz_0_yz_xx = buffer_dsdd[96];

    auto g_xz_0_yz_xy = buffer_dsdd[97];

    auto g_xz_0_yz_xz = buffer_dsdd[98];

    auto g_xz_0_yz_yy = buffer_dsdd[99];

    auto g_xz_0_yz_yz = buffer_dsdd[100];

    auto g_xz_0_yz_zz = buffer_dsdd[101];

    auto g_xz_0_zz_xx = buffer_dsdd[102];

    auto g_xz_0_zz_xy = buffer_dsdd[103];

    auto g_xz_0_zz_xz = buffer_dsdd[104];

    auto g_xz_0_zz_yy = buffer_dsdd[105];

    auto g_xz_0_zz_yz = buffer_dsdd[106];

    auto g_xz_0_zz_zz = buffer_dsdd[107];

    auto g_yy_0_xx_xx = buffer_dsdd[108];

    auto g_yy_0_xx_xy = buffer_dsdd[109];

    auto g_yy_0_xx_xz = buffer_dsdd[110];

    auto g_yy_0_xx_yy = buffer_dsdd[111];

    auto g_yy_0_xx_yz = buffer_dsdd[112];

    auto g_yy_0_xx_zz = buffer_dsdd[113];

    auto g_yy_0_xy_xx = buffer_dsdd[114];

    auto g_yy_0_xy_xy = buffer_dsdd[115];

    auto g_yy_0_xy_xz = buffer_dsdd[116];

    auto g_yy_0_xy_yy = buffer_dsdd[117];

    auto g_yy_0_xy_yz = buffer_dsdd[118];

    auto g_yy_0_xy_zz = buffer_dsdd[119];

    auto g_yy_0_xz_xx = buffer_dsdd[120];

    auto g_yy_0_xz_xy = buffer_dsdd[121];

    auto g_yy_0_xz_xz = buffer_dsdd[122];

    auto g_yy_0_xz_yy = buffer_dsdd[123];

    auto g_yy_0_xz_yz = buffer_dsdd[124];

    auto g_yy_0_xz_zz = buffer_dsdd[125];

    auto g_yy_0_yy_xx = buffer_dsdd[126];

    auto g_yy_0_yy_xy = buffer_dsdd[127];

    auto g_yy_0_yy_xz = buffer_dsdd[128];

    auto g_yy_0_yy_yy = buffer_dsdd[129];

    auto g_yy_0_yy_yz = buffer_dsdd[130];

    auto g_yy_0_yy_zz = buffer_dsdd[131];

    auto g_yy_0_yz_xx = buffer_dsdd[132];

    auto g_yy_0_yz_xy = buffer_dsdd[133];

    auto g_yy_0_yz_xz = buffer_dsdd[134];

    auto g_yy_0_yz_yy = buffer_dsdd[135];

    auto g_yy_0_yz_yz = buffer_dsdd[136];

    auto g_yy_0_yz_zz = buffer_dsdd[137];

    auto g_yy_0_zz_xx = buffer_dsdd[138];

    auto g_yy_0_zz_xy = buffer_dsdd[139];

    auto g_yy_0_zz_xz = buffer_dsdd[140];

    auto g_yy_0_zz_yy = buffer_dsdd[141];

    auto g_yy_0_zz_yz = buffer_dsdd[142];

    auto g_yy_0_zz_zz = buffer_dsdd[143];

    auto g_yz_0_xx_xx = buffer_dsdd[144];

    auto g_yz_0_xx_xy = buffer_dsdd[145];

    auto g_yz_0_xx_xz = buffer_dsdd[146];

    auto g_yz_0_xx_yy = buffer_dsdd[147];

    auto g_yz_0_xx_yz = buffer_dsdd[148];

    auto g_yz_0_xx_zz = buffer_dsdd[149];

    auto g_yz_0_xy_xx = buffer_dsdd[150];

    auto g_yz_0_xy_xy = buffer_dsdd[151];

    auto g_yz_0_xy_xz = buffer_dsdd[152];

    auto g_yz_0_xy_yy = buffer_dsdd[153];

    auto g_yz_0_xy_yz = buffer_dsdd[154];

    auto g_yz_0_xy_zz = buffer_dsdd[155];

    auto g_yz_0_xz_xx = buffer_dsdd[156];

    auto g_yz_0_xz_xy = buffer_dsdd[157];

    auto g_yz_0_xz_xz = buffer_dsdd[158];

    auto g_yz_0_xz_yy = buffer_dsdd[159];

    auto g_yz_0_xz_yz = buffer_dsdd[160];

    auto g_yz_0_xz_zz = buffer_dsdd[161];

    auto g_yz_0_yy_xx = buffer_dsdd[162];

    auto g_yz_0_yy_xy = buffer_dsdd[163];

    auto g_yz_0_yy_xz = buffer_dsdd[164];

    auto g_yz_0_yy_yy = buffer_dsdd[165];

    auto g_yz_0_yy_yz = buffer_dsdd[166];

    auto g_yz_0_yy_zz = buffer_dsdd[167];

    auto g_yz_0_yz_xx = buffer_dsdd[168];

    auto g_yz_0_yz_xy = buffer_dsdd[169];

    auto g_yz_0_yz_xz = buffer_dsdd[170];

    auto g_yz_0_yz_yy = buffer_dsdd[171];

    auto g_yz_0_yz_yz = buffer_dsdd[172];

    auto g_yz_0_yz_zz = buffer_dsdd[173];

    auto g_yz_0_zz_xx = buffer_dsdd[174];

    auto g_yz_0_zz_xy = buffer_dsdd[175];

    auto g_yz_0_zz_xz = buffer_dsdd[176];

    auto g_yz_0_zz_yy = buffer_dsdd[177];

    auto g_yz_0_zz_yz = buffer_dsdd[178];

    auto g_yz_0_zz_zz = buffer_dsdd[179];

    auto g_zz_0_xx_xx = buffer_dsdd[180];

    auto g_zz_0_xx_xy = buffer_dsdd[181];

    auto g_zz_0_xx_xz = buffer_dsdd[182];

    auto g_zz_0_xx_yy = buffer_dsdd[183];

    auto g_zz_0_xx_yz = buffer_dsdd[184];

    auto g_zz_0_xx_zz = buffer_dsdd[185];

    auto g_zz_0_xy_xx = buffer_dsdd[186];

    auto g_zz_0_xy_xy = buffer_dsdd[187];

    auto g_zz_0_xy_xz = buffer_dsdd[188];

    auto g_zz_0_xy_yy = buffer_dsdd[189];

    auto g_zz_0_xy_yz = buffer_dsdd[190];

    auto g_zz_0_xy_zz = buffer_dsdd[191];

    auto g_zz_0_xz_xx = buffer_dsdd[192];

    auto g_zz_0_xz_xy = buffer_dsdd[193];

    auto g_zz_0_xz_xz = buffer_dsdd[194];

    auto g_zz_0_xz_yy = buffer_dsdd[195];

    auto g_zz_0_xz_yz = buffer_dsdd[196];

    auto g_zz_0_xz_zz = buffer_dsdd[197];

    auto g_zz_0_yy_xx = buffer_dsdd[198];

    auto g_zz_0_yy_xy = buffer_dsdd[199];

    auto g_zz_0_yy_xz = buffer_dsdd[200];

    auto g_zz_0_yy_yy = buffer_dsdd[201];

    auto g_zz_0_yy_yz = buffer_dsdd[202];

    auto g_zz_0_yy_zz = buffer_dsdd[203];

    auto g_zz_0_yz_xx = buffer_dsdd[204];

    auto g_zz_0_yz_xy = buffer_dsdd[205];

    auto g_zz_0_yz_xz = buffer_dsdd[206];

    auto g_zz_0_yz_yy = buffer_dsdd[207];

    auto g_zz_0_yz_yz = buffer_dsdd[208];

    auto g_zz_0_yz_zz = buffer_dsdd[209];

    auto g_zz_0_zz_xx = buffer_dsdd[210];

    auto g_zz_0_zz_xy = buffer_dsdd[211];

    auto g_zz_0_zz_xz = buffer_dsdd[212];

    auto g_zz_0_zz_yy = buffer_dsdd[213];

    auto g_zz_0_zz_yz = buffer_dsdd[214];

    auto g_zz_0_zz_zz = buffer_dsdd[215];

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

    /// Set up components of integrals buffer : buffer_1100_ppdd

    auto g_x_x_0_0_x_x_xx_xx = buffer_1100_ppdd[0];

    auto g_x_x_0_0_x_x_xx_xy = buffer_1100_ppdd[1];

    auto g_x_x_0_0_x_x_xx_xz = buffer_1100_ppdd[2];

    auto g_x_x_0_0_x_x_xx_yy = buffer_1100_ppdd[3];

    auto g_x_x_0_0_x_x_xx_yz = buffer_1100_ppdd[4];

    auto g_x_x_0_0_x_x_xx_zz = buffer_1100_ppdd[5];

    auto g_x_x_0_0_x_x_xy_xx = buffer_1100_ppdd[6];

    auto g_x_x_0_0_x_x_xy_xy = buffer_1100_ppdd[7];

    auto g_x_x_0_0_x_x_xy_xz = buffer_1100_ppdd[8];

    auto g_x_x_0_0_x_x_xy_yy = buffer_1100_ppdd[9];

    auto g_x_x_0_0_x_x_xy_yz = buffer_1100_ppdd[10];

    auto g_x_x_0_0_x_x_xy_zz = buffer_1100_ppdd[11];

    auto g_x_x_0_0_x_x_xz_xx = buffer_1100_ppdd[12];

    auto g_x_x_0_0_x_x_xz_xy = buffer_1100_ppdd[13];

    auto g_x_x_0_0_x_x_xz_xz = buffer_1100_ppdd[14];

    auto g_x_x_0_0_x_x_xz_yy = buffer_1100_ppdd[15];

    auto g_x_x_0_0_x_x_xz_yz = buffer_1100_ppdd[16];

    auto g_x_x_0_0_x_x_xz_zz = buffer_1100_ppdd[17];

    auto g_x_x_0_0_x_x_yy_xx = buffer_1100_ppdd[18];

    auto g_x_x_0_0_x_x_yy_xy = buffer_1100_ppdd[19];

    auto g_x_x_0_0_x_x_yy_xz = buffer_1100_ppdd[20];

    auto g_x_x_0_0_x_x_yy_yy = buffer_1100_ppdd[21];

    auto g_x_x_0_0_x_x_yy_yz = buffer_1100_ppdd[22];

    auto g_x_x_0_0_x_x_yy_zz = buffer_1100_ppdd[23];

    auto g_x_x_0_0_x_x_yz_xx = buffer_1100_ppdd[24];

    auto g_x_x_0_0_x_x_yz_xy = buffer_1100_ppdd[25];

    auto g_x_x_0_0_x_x_yz_xz = buffer_1100_ppdd[26];

    auto g_x_x_0_0_x_x_yz_yy = buffer_1100_ppdd[27];

    auto g_x_x_0_0_x_x_yz_yz = buffer_1100_ppdd[28];

    auto g_x_x_0_0_x_x_yz_zz = buffer_1100_ppdd[29];

    auto g_x_x_0_0_x_x_zz_xx = buffer_1100_ppdd[30];

    auto g_x_x_0_0_x_x_zz_xy = buffer_1100_ppdd[31];

    auto g_x_x_0_0_x_x_zz_xz = buffer_1100_ppdd[32];

    auto g_x_x_0_0_x_x_zz_yy = buffer_1100_ppdd[33];

    auto g_x_x_0_0_x_x_zz_yz = buffer_1100_ppdd[34];

    auto g_x_x_0_0_x_x_zz_zz = buffer_1100_ppdd[35];

    auto g_x_x_0_0_x_y_xx_xx = buffer_1100_ppdd[36];

    auto g_x_x_0_0_x_y_xx_xy = buffer_1100_ppdd[37];

    auto g_x_x_0_0_x_y_xx_xz = buffer_1100_ppdd[38];

    auto g_x_x_0_0_x_y_xx_yy = buffer_1100_ppdd[39];

    auto g_x_x_0_0_x_y_xx_yz = buffer_1100_ppdd[40];

    auto g_x_x_0_0_x_y_xx_zz = buffer_1100_ppdd[41];

    auto g_x_x_0_0_x_y_xy_xx = buffer_1100_ppdd[42];

    auto g_x_x_0_0_x_y_xy_xy = buffer_1100_ppdd[43];

    auto g_x_x_0_0_x_y_xy_xz = buffer_1100_ppdd[44];

    auto g_x_x_0_0_x_y_xy_yy = buffer_1100_ppdd[45];

    auto g_x_x_0_0_x_y_xy_yz = buffer_1100_ppdd[46];

    auto g_x_x_0_0_x_y_xy_zz = buffer_1100_ppdd[47];

    auto g_x_x_0_0_x_y_xz_xx = buffer_1100_ppdd[48];

    auto g_x_x_0_0_x_y_xz_xy = buffer_1100_ppdd[49];

    auto g_x_x_0_0_x_y_xz_xz = buffer_1100_ppdd[50];

    auto g_x_x_0_0_x_y_xz_yy = buffer_1100_ppdd[51];

    auto g_x_x_0_0_x_y_xz_yz = buffer_1100_ppdd[52];

    auto g_x_x_0_0_x_y_xz_zz = buffer_1100_ppdd[53];

    auto g_x_x_0_0_x_y_yy_xx = buffer_1100_ppdd[54];

    auto g_x_x_0_0_x_y_yy_xy = buffer_1100_ppdd[55];

    auto g_x_x_0_0_x_y_yy_xz = buffer_1100_ppdd[56];

    auto g_x_x_0_0_x_y_yy_yy = buffer_1100_ppdd[57];

    auto g_x_x_0_0_x_y_yy_yz = buffer_1100_ppdd[58];

    auto g_x_x_0_0_x_y_yy_zz = buffer_1100_ppdd[59];

    auto g_x_x_0_0_x_y_yz_xx = buffer_1100_ppdd[60];

    auto g_x_x_0_0_x_y_yz_xy = buffer_1100_ppdd[61];

    auto g_x_x_0_0_x_y_yz_xz = buffer_1100_ppdd[62];

    auto g_x_x_0_0_x_y_yz_yy = buffer_1100_ppdd[63];

    auto g_x_x_0_0_x_y_yz_yz = buffer_1100_ppdd[64];

    auto g_x_x_0_0_x_y_yz_zz = buffer_1100_ppdd[65];

    auto g_x_x_0_0_x_y_zz_xx = buffer_1100_ppdd[66];

    auto g_x_x_0_0_x_y_zz_xy = buffer_1100_ppdd[67];

    auto g_x_x_0_0_x_y_zz_xz = buffer_1100_ppdd[68];

    auto g_x_x_0_0_x_y_zz_yy = buffer_1100_ppdd[69];

    auto g_x_x_0_0_x_y_zz_yz = buffer_1100_ppdd[70];

    auto g_x_x_0_0_x_y_zz_zz = buffer_1100_ppdd[71];

    auto g_x_x_0_0_x_z_xx_xx = buffer_1100_ppdd[72];

    auto g_x_x_0_0_x_z_xx_xy = buffer_1100_ppdd[73];

    auto g_x_x_0_0_x_z_xx_xz = buffer_1100_ppdd[74];

    auto g_x_x_0_0_x_z_xx_yy = buffer_1100_ppdd[75];

    auto g_x_x_0_0_x_z_xx_yz = buffer_1100_ppdd[76];

    auto g_x_x_0_0_x_z_xx_zz = buffer_1100_ppdd[77];

    auto g_x_x_0_0_x_z_xy_xx = buffer_1100_ppdd[78];

    auto g_x_x_0_0_x_z_xy_xy = buffer_1100_ppdd[79];

    auto g_x_x_0_0_x_z_xy_xz = buffer_1100_ppdd[80];

    auto g_x_x_0_0_x_z_xy_yy = buffer_1100_ppdd[81];

    auto g_x_x_0_0_x_z_xy_yz = buffer_1100_ppdd[82];

    auto g_x_x_0_0_x_z_xy_zz = buffer_1100_ppdd[83];

    auto g_x_x_0_0_x_z_xz_xx = buffer_1100_ppdd[84];

    auto g_x_x_0_0_x_z_xz_xy = buffer_1100_ppdd[85];

    auto g_x_x_0_0_x_z_xz_xz = buffer_1100_ppdd[86];

    auto g_x_x_0_0_x_z_xz_yy = buffer_1100_ppdd[87];

    auto g_x_x_0_0_x_z_xz_yz = buffer_1100_ppdd[88];

    auto g_x_x_0_0_x_z_xz_zz = buffer_1100_ppdd[89];

    auto g_x_x_0_0_x_z_yy_xx = buffer_1100_ppdd[90];

    auto g_x_x_0_0_x_z_yy_xy = buffer_1100_ppdd[91];

    auto g_x_x_0_0_x_z_yy_xz = buffer_1100_ppdd[92];

    auto g_x_x_0_0_x_z_yy_yy = buffer_1100_ppdd[93];

    auto g_x_x_0_0_x_z_yy_yz = buffer_1100_ppdd[94];

    auto g_x_x_0_0_x_z_yy_zz = buffer_1100_ppdd[95];

    auto g_x_x_0_0_x_z_yz_xx = buffer_1100_ppdd[96];

    auto g_x_x_0_0_x_z_yz_xy = buffer_1100_ppdd[97];

    auto g_x_x_0_0_x_z_yz_xz = buffer_1100_ppdd[98];

    auto g_x_x_0_0_x_z_yz_yy = buffer_1100_ppdd[99];

    auto g_x_x_0_0_x_z_yz_yz = buffer_1100_ppdd[100];

    auto g_x_x_0_0_x_z_yz_zz = buffer_1100_ppdd[101];

    auto g_x_x_0_0_x_z_zz_xx = buffer_1100_ppdd[102];

    auto g_x_x_0_0_x_z_zz_xy = buffer_1100_ppdd[103];

    auto g_x_x_0_0_x_z_zz_xz = buffer_1100_ppdd[104];

    auto g_x_x_0_0_x_z_zz_yy = buffer_1100_ppdd[105];

    auto g_x_x_0_0_x_z_zz_yz = buffer_1100_ppdd[106];

    auto g_x_x_0_0_x_z_zz_zz = buffer_1100_ppdd[107];

    auto g_x_x_0_0_y_x_xx_xx = buffer_1100_ppdd[108];

    auto g_x_x_0_0_y_x_xx_xy = buffer_1100_ppdd[109];

    auto g_x_x_0_0_y_x_xx_xz = buffer_1100_ppdd[110];

    auto g_x_x_0_0_y_x_xx_yy = buffer_1100_ppdd[111];

    auto g_x_x_0_0_y_x_xx_yz = buffer_1100_ppdd[112];

    auto g_x_x_0_0_y_x_xx_zz = buffer_1100_ppdd[113];

    auto g_x_x_0_0_y_x_xy_xx = buffer_1100_ppdd[114];

    auto g_x_x_0_0_y_x_xy_xy = buffer_1100_ppdd[115];

    auto g_x_x_0_0_y_x_xy_xz = buffer_1100_ppdd[116];

    auto g_x_x_0_0_y_x_xy_yy = buffer_1100_ppdd[117];

    auto g_x_x_0_0_y_x_xy_yz = buffer_1100_ppdd[118];

    auto g_x_x_0_0_y_x_xy_zz = buffer_1100_ppdd[119];

    auto g_x_x_0_0_y_x_xz_xx = buffer_1100_ppdd[120];

    auto g_x_x_0_0_y_x_xz_xy = buffer_1100_ppdd[121];

    auto g_x_x_0_0_y_x_xz_xz = buffer_1100_ppdd[122];

    auto g_x_x_0_0_y_x_xz_yy = buffer_1100_ppdd[123];

    auto g_x_x_0_0_y_x_xz_yz = buffer_1100_ppdd[124];

    auto g_x_x_0_0_y_x_xz_zz = buffer_1100_ppdd[125];

    auto g_x_x_0_0_y_x_yy_xx = buffer_1100_ppdd[126];

    auto g_x_x_0_0_y_x_yy_xy = buffer_1100_ppdd[127];

    auto g_x_x_0_0_y_x_yy_xz = buffer_1100_ppdd[128];

    auto g_x_x_0_0_y_x_yy_yy = buffer_1100_ppdd[129];

    auto g_x_x_0_0_y_x_yy_yz = buffer_1100_ppdd[130];

    auto g_x_x_0_0_y_x_yy_zz = buffer_1100_ppdd[131];

    auto g_x_x_0_0_y_x_yz_xx = buffer_1100_ppdd[132];

    auto g_x_x_0_0_y_x_yz_xy = buffer_1100_ppdd[133];

    auto g_x_x_0_0_y_x_yz_xz = buffer_1100_ppdd[134];

    auto g_x_x_0_0_y_x_yz_yy = buffer_1100_ppdd[135];

    auto g_x_x_0_0_y_x_yz_yz = buffer_1100_ppdd[136];

    auto g_x_x_0_0_y_x_yz_zz = buffer_1100_ppdd[137];

    auto g_x_x_0_0_y_x_zz_xx = buffer_1100_ppdd[138];

    auto g_x_x_0_0_y_x_zz_xy = buffer_1100_ppdd[139];

    auto g_x_x_0_0_y_x_zz_xz = buffer_1100_ppdd[140];

    auto g_x_x_0_0_y_x_zz_yy = buffer_1100_ppdd[141];

    auto g_x_x_0_0_y_x_zz_yz = buffer_1100_ppdd[142];

    auto g_x_x_0_0_y_x_zz_zz = buffer_1100_ppdd[143];

    auto g_x_x_0_0_y_y_xx_xx = buffer_1100_ppdd[144];

    auto g_x_x_0_0_y_y_xx_xy = buffer_1100_ppdd[145];

    auto g_x_x_0_0_y_y_xx_xz = buffer_1100_ppdd[146];

    auto g_x_x_0_0_y_y_xx_yy = buffer_1100_ppdd[147];

    auto g_x_x_0_0_y_y_xx_yz = buffer_1100_ppdd[148];

    auto g_x_x_0_0_y_y_xx_zz = buffer_1100_ppdd[149];

    auto g_x_x_0_0_y_y_xy_xx = buffer_1100_ppdd[150];

    auto g_x_x_0_0_y_y_xy_xy = buffer_1100_ppdd[151];

    auto g_x_x_0_0_y_y_xy_xz = buffer_1100_ppdd[152];

    auto g_x_x_0_0_y_y_xy_yy = buffer_1100_ppdd[153];

    auto g_x_x_0_0_y_y_xy_yz = buffer_1100_ppdd[154];

    auto g_x_x_0_0_y_y_xy_zz = buffer_1100_ppdd[155];

    auto g_x_x_0_0_y_y_xz_xx = buffer_1100_ppdd[156];

    auto g_x_x_0_0_y_y_xz_xy = buffer_1100_ppdd[157];

    auto g_x_x_0_0_y_y_xz_xz = buffer_1100_ppdd[158];

    auto g_x_x_0_0_y_y_xz_yy = buffer_1100_ppdd[159];

    auto g_x_x_0_0_y_y_xz_yz = buffer_1100_ppdd[160];

    auto g_x_x_0_0_y_y_xz_zz = buffer_1100_ppdd[161];

    auto g_x_x_0_0_y_y_yy_xx = buffer_1100_ppdd[162];

    auto g_x_x_0_0_y_y_yy_xy = buffer_1100_ppdd[163];

    auto g_x_x_0_0_y_y_yy_xz = buffer_1100_ppdd[164];

    auto g_x_x_0_0_y_y_yy_yy = buffer_1100_ppdd[165];

    auto g_x_x_0_0_y_y_yy_yz = buffer_1100_ppdd[166];

    auto g_x_x_0_0_y_y_yy_zz = buffer_1100_ppdd[167];

    auto g_x_x_0_0_y_y_yz_xx = buffer_1100_ppdd[168];

    auto g_x_x_0_0_y_y_yz_xy = buffer_1100_ppdd[169];

    auto g_x_x_0_0_y_y_yz_xz = buffer_1100_ppdd[170];

    auto g_x_x_0_0_y_y_yz_yy = buffer_1100_ppdd[171];

    auto g_x_x_0_0_y_y_yz_yz = buffer_1100_ppdd[172];

    auto g_x_x_0_0_y_y_yz_zz = buffer_1100_ppdd[173];

    auto g_x_x_0_0_y_y_zz_xx = buffer_1100_ppdd[174];

    auto g_x_x_0_0_y_y_zz_xy = buffer_1100_ppdd[175];

    auto g_x_x_0_0_y_y_zz_xz = buffer_1100_ppdd[176];

    auto g_x_x_0_0_y_y_zz_yy = buffer_1100_ppdd[177];

    auto g_x_x_0_0_y_y_zz_yz = buffer_1100_ppdd[178];

    auto g_x_x_0_0_y_y_zz_zz = buffer_1100_ppdd[179];

    auto g_x_x_0_0_y_z_xx_xx = buffer_1100_ppdd[180];

    auto g_x_x_0_0_y_z_xx_xy = buffer_1100_ppdd[181];

    auto g_x_x_0_0_y_z_xx_xz = buffer_1100_ppdd[182];

    auto g_x_x_0_0_y_z_xx_yy = buffer_1100_ppdd[183];

    auto g_x_x_0_0_y_z_xx_yz = buffer_1100_ppdd[184];

    auto g_x_x_0_0_y_z_xx_zz = buffer_1100_ppdd[185];

    auto g_x_x_0_0_y_z_xy_xx = buffer_1100_ppdd[186];

    auto g_x_x_0_0_y_z_xy_xy = buffer_1100_ppdd[187];

    auto g_x_x_0_0_y_z_xy_xz = buffer_1100_ppdd[188];

    auto g_x_x_0_0_y_z_xy_yy = buffer_1100_ppdd[189];

    auto g_x_x_0_0_y_z_xy_yz = buffer_1100_ppdd[190];

    auto g_x_x_0_0_y_z_xy_zz = buffer_1100_ppdd[191];

    auto g_x_x_0_0_y_z_xz_xx = buffer_1100_ppdd[192];

    auto g_x_x_0_0_y_z_xz_xy = buffer_1100_ppdd[193];

    auto g_x_x_0_0_y_z_xz_xz = buffer_1100_ppdd[194];

    auto g_x_x_0_0_y_z_xz_yy = buffer_1100_ppdd[195];

    auto g_x_x_0_0_y_z_xz_yz = buffer_1100_ppdd[196];

    auto g_x_x_0_0_y_z_xz_zz = buffer_1100_ppdd[197];

    auto g_x_x_0_0_y_z_yy_xx = buffer_1100_ppdd[198];

    auto g_x_x_0_0_y_z_yy_xy = buffer_1100_ppdd[199];

    auto g_x_x_0_0_y_z_yy_xz = buffer_1100_ppdd[200];

    auto g_x_x_0_0_y_z_yy_yy = buffer_1100_ppdd[201];

    auto g_x_x_0_0_y_z_yy_yz = buffer_1100_ppdd[202];

    auto g_x_x_0_0_y_z_yy_zz = buffer_1100_ppdd[203];

    auto g_x_x_0_0_y_z_yz_xx = buffer_1100_ppdd[204];

    auto g_x_x_0_0_y_z_yz_xy = buffer_1100_ppdd[205];

    auto g_x_x_0_0_y_z_yz_xz = buffer_1100_ppdd[206];

    auto g_x_x_0_0_y_z_yz_yy = buffer_1100_ppdd[207];

    auto g_x_x_0_0_y_z_yz_yz = buffer_1100_ppdd[208];

    auto g_x_x_0_0_y_z_yz_zz = buffer_1100_ppdd[209];

    auto g_x_x_0_0_y_z_zz_xx = buffer_1100_ppdd[210];

    auto g_x_x_0_0_y_z_zz_xy = buffer_1100_ppdd[211];

    auto g_x_x_0_0_y_z_zz_xz = buffer_1100_ppdd[212];

    auto g_x_x_0_0_y_z_zz_yy = buffer_1100_ppdd[213];

    auto g_x_x_0_0_y_z_zz_yz = buffer_1100_ppdd[214];

    auto g_x_x_0_0_y_z_zz_zz = buffer_1100_ppdd[215];

    auto g_x_x_0_0_z_x_xx_xx = buffer_1100_ppdd[216];

    auto g_x_x_0_0_z_x_xx_xy = buffer_1100_ppdd[217];

    auto g_x_x_0_0_z_x_xx_xz = buffer_1100_ppdd[218];

    auto g_x_x_0_0_z_x_xx_yy = buffer_1100_ppdd[219];

    auto g_x_x_0_0_z_x_xx_yz = buffer_1100_ppdd[220];

    auto g_x_x_0_0_z_x_xx_zz = buffer_1100_ppdd[221];

    auto g_x_x_0_0_z_x_xy_xx = buffer_1100_ppdd[222];

    auto g_x_x_0_0_z_x_xy_xy = buffer_1100_ppdd[223];

    auto g_x_x_0_0_z_x_xy_xz = buffer_1100_ppdd[224];

    auto g_x_x_0_0_z_x_xy_yy = buffer_1100_ppdd[225];

    auto g_x_x_0_0_z_x_xy_yz = buffer_1100_ppdd[226];

    auto g_x_x_0_0_z_x_xy_zz = buffer_1100_ppdd[227];

    auto g_x_x_0_0_z_x_xz_xx = buffer_1100_ppdd[228];

    auto g_x_x_0_0_z_x_xz_xy = buffer_1100_ppdd[229];

    auto g_x_x_0_0_z_x_xz_xz = buffer_1100_ppdd[230];

    auto g_x_x_0_0_z_x_xz_yy = buffer_1100_ppdd[231];

    auto g_x_x_0_0_z_x_xz_yz = buffer_1100_ppdd[232];

    auto g_x_x_0_0_z_x_xz_zz = buffer_1100_ppdd[233];

    auto g_x_x_0_0_z_x_yy_xx = buffer_1100_ppdd[234];

    auto g_x_x_0_0_z_x_yy_xy = buffer_1100_ppdd[235];

    auto g_x_x_0_0_z_x_yy_xz = buffer_1100_ppdd[236];

    auto g_x_x_0_0_z_x_yy_yy = buffer_1100_ppdd[237];

    auto g_x_x_0_0_z_x_yy_yz = buffer_1100_ppdd[238];

    auto g_x_x_0_0_z_x_yy_zz = buffer_1100_ppdd[239];

    auto g_x_x_0_0_z_x_yz_xx = buffer_1100_ppdd[240];

    auto g_x_x_0_0_z_x_yz_xy = buffer_1100_ppdd[241];

    auto g_x_x_0_0_z_x_yz_xz = buffer_1100_ppdd[242];

    auto g_x_x_0_0_z_x_yz_yy = buffer_1100_ppdd[243];

    auto g_x_x_0_0_z_x_yz_yz = buffer_1100_ppdd[244];

    auto g_x_x_0_0_z_x_yz_zz = buffer_1100_ppdd[245];

    auto g_x_x_0_0_z_x_zz_xx = buffer_1100_ppdd[246];

    auto g_x_x_0_0_z_x_zz_xy = buffer_1100_ppdd[247];

    auto g_x_x_0_0_z_x_zz_xz = buffer_1100_ppdd[248];

    auto g_x_x_0_0_z_x_zz_yy = buffer_1100_ppdd[249];

    auto g_x_x_0_0_z_x_zz_yz = buffer_1100_ppdd[250];

    auto g_x_x_0_0_z_x_zz_zz = buffer_1100_ppdd[251];

    auto g_x_x_0_0_z_y_xx_xx = buffer_1100_ppdd[252];

    auto g_x_x_0_0_z_y_xx_xy = buffer_1100_ppdd[253];

    auto g_x_x_0_0_z_y_xx_xz = buffer_1100_ppdd[254];

    auto g_x_x_0_0_z_y_xx_yy = buffer_1100_ppdd[255];

    auto g_x_x_0_0_z_y_xx_yz = buffer_1100_ppdd[256];

    auto g_x_x_0_0_z_y_xx_zz = buffer_1100_ppdd[257];

    auto g_x_x_0_0_z_y_xy_xx = buffer_1100_ppdd[258];

    auto g_x_x_0_0_z_y_xy_xy = buffer_1100_ppdd[259];

    auto g_x_x_0_0_z_y_xy_xz = buffer_1100_ppdd[260];

    auto g_x_x_0_0_z_y_xy_yy = buffer_1100_ppdd[261];

    auto g_x_x_0_0_z_y_xy_yz = buffer_1100_ppdd[262];

    auto g_x_x_0_0_z_y_xy_zz = buffer_1100_ppdd[263];

    auto g_x_x_0_0_z_y_xz_xx = buffer_1100_ppdd[264];

    auto g_x_x_0_0_z_y_xz_xy = buffer_1100_ppdd[265];

    auto g_x_x_0_0_z_y_xz_xz = buffer_1100_ppdd[266];

    auto g_x_x_0_0_z_y_xz_yy = buffer_1100_ppdd[267];

    auto g_x_x_0_0_z_y_xz_yz = buffer_1100_ppdd[268];

    auto g_x_x_0_0_z_y_xz_zz = buffer_1100_ppdd[269];

    auto g_x_x_0_0_z_y_yy_xx = buffer_1100_ppdd[270];

    auto g_x_x_0_0_z_y_yy_xy = buffer_1100_ppdd[271];

    auto g_x_x_0_0_z_y_yy_xz = buffer_1100_ppdd[272];

    auto g_x_x_0_0_z_y_yy_yy = buffer_1100_ppdd[273];

    auto g_x_x_0_0_z_y_yy_yz = buffer_1100_ppdd[274];

    auto g_x_x_0_0_z_y_yy_zz = buffer_1100_ppdd[275];

    auto g_x_x_0_0_z_y_yz_xx = buffer_1100_ppdd[276];

    auto g_x_x_0_0_z_y_yz_xy = buffer_1100_ppdd[277];

    auto g_x_x_0_0_z_y_yz_xz = buffer_1100_ppdd[278];

    auto g_x_x_0_0_z_y_yz_yy = buffer_1100_ppdd[279];

    auto g_x_x_0_0_z_y_yz_yz = buffer_1100_ppdd[280];

    auto g_x_x_0_0_z_y_yz_zz = buffer_1100_ppdd[281];

    auto g_x_x_0_0_z_y_zz_xx = buffer_1100_ppdd[282];

    auto g_x_x_0_0_z_y_zz_xy = buffer_1100_ppdd[283];

    auto g_x_x_0_0_z_y_zz_xz = buffer_1100_ppdd[284];

    auto g_x_x_0_0_z_y_zz_yy = buffer_1100_ppdd[285];

    auto g_x_x_0_0_z_y_zz_yz = buffer_1100_ppdd[286];

    auto g_x_x_0_0_z_y_zz_zz = buffer_1100_ppdd[287];

    auto g_x_x_0_0_z_z_xx_xx = buffer_1100_ppdd[288];

    auto g_x_x_0_0_z_z_xx_xy = buffer_1100_ppdd[289];

    auto g_x_x_0_0_z_z_xx_xz = buffer_1100_ppdd[290];

    auto g_x_x_0_0_z_z_xx_yy = buffer_1100_ppdd[291];

    auto g_x_x_0_0_z_z_xx_yz = buffer_1100_ppdd[292];

    auto g_x_x_0_0_z_z_xx_zz = buffer_1100_ppdd[293];

    auto g_x_x_0_0_z_z_xy_xx = buffer_1100_ppdd[294];

    auto g_x_x_0_0_z_z_xy_xy = buffer_1100_ppdd[295];

    auto g_x_x_0_0_z_z_xy_xz = buffer_1100_ppdd[296];

    auto g_x_x_0_0_z_z_xy_yy = buffer_1100_ppdd[297];

    auto g_x_x_0_0_z_z_xy_yz = buffer_1100_ppdd[298];

    auto g_x_x_0_0_z_z_xy_zz = buffer_1100_ppdd[299];

    auto g_x_x_0_0_z_z_xz_xx = buffer_1100_ppdd[300];

    auto g_x_x_0_0_z_z_xz_xy = buffer_1100_ppdd[301];

    auto g_x_x_0_0_z_z_xz_xz = buffer_1100_ppdd[302];

    auto g_x_x_0_0_z_z_xz_yy = buffer_1100_ppdd[303];

    auto g_x_x_0_0_z_z_xz_yz = buffer_1100_ppdd[304];

    auto g_x_x_0_0_z_z_xz_zz = buffer_1100_ppdd[305];

    auto g_x_x_0_0_z_z_yy_xx = buffer_1100_ppdd[306];

    auto g_x_x_0_0_z_z_yy_xy = buffer_1100_ppdd[307];

    auto g_x_x_0_0_z_z_yy_xz = buffer_1100_ppdd[308];

    auto g_x_x_0_0_z_z_yy_yy = buffer_1100_ppdd[309];

    auto g_x_x_0_0_z_z_yy_yz = buffer_1100_ppdd[310];

    auto g_x_x_0_0_z_z_yy_zz = buffer_1100_ppdd[311];

    auto g_x_x_0_0_z_z_yz_xx = buffer_1100_ppdd[312];

    auto g_x_x_0_0_z_z_yz_xy = buffer_1100_ppdd[313];

    auto g_x_x_0_0_z_z_yz_xz = buffer_1100_ppdd[314];

    auto g_x_x_0_0_z_z_yz_yy = buffer_1100_ppdd[315];

    auto g_x_x_0_0_z_z_yz_yz = buffer_1100_ppdd[316];

    auto g_x_x_0_0_z_z_yz_zz = buffer_1100_ppdd[317];

    auto g_x_x_0_0_z_z_zz_xx = buffer_1100_ppdd[318];

    auto g_x_x_0_0_z_z_zz_xy = buffer_1100_ppdd[319];

    auto g_x_x_0_0_z_z_zz_xz = buffer_1100_ppdd[320];

    auto g_x_x_0_0_z_z_zz_yy = buffer_1100_ppdd[321];

    auto g_x_x_0_0_z_z_zz_yz = buffer_1100_ppdd[322];

    auto g_x_x_0_0_z_z_zz_zz = buffer_1100_ppdd[323];

    auto g_x_y_0_0_x_x_xx_xx = buffer_1100_ppdd[324];

    auto g_x_y_0_0_x_x_xx_xy = buffer_1100_ppdd[325];

    auto g_x_y_0_0_x_x_xx_xz = buffer_1100_ppdd[326];

    auto g_x_y_0_0_x_x_xx_yy = buffer_1100_ppdd[327];

    auto g_x_y_0_0_x_x_xx_yz = buffer_1100_ppdd[328];

    auto g_x_y_0_0_x_x_xx_zz = buffer_1100_ppdd[329];

    auto g_x_y_0_0_x_x_xy_xx = buffer_1100_ppdd[330];

    auto g_x_y_0_0_x_x_xy_xy = buffer_1100_ppdd[331];

    auto g_x_y_0_0_x_x_xy_xz = buffer_1100_ppdd[332];

    auto g_x_y_0_0_x_x_xy_yy = buffer_1100_ppdd[333];

    auto g_x_y_0_0_x_x_xy_yz = buffer_1100_ppdd[334];

    auto g_x_y_0_0_x_x_xy_zz = buffer_1100_ppdd[335];

    auto g_x_y_0_0_x_x_xz_xx = buffer_1100_ppdd[336];

    auto g_x_y_0_0_x_x_xz_xy = buffer_1100_ppdd[337];

    auto g_x_y_0_0_x_x_xz_xz = buffer_1100_ppdd[338];

    auto g_x_y_0_0_x_x_xz_yy = buffer_1100_ppdd[339];

    auto g_x_y_0_0_x_x_xz_yz = buffer_1100_ppdd[340];

    auto g_x_y_0_0_x_x_xz_zz = buffer_1100_ppdd[341];

    auto g_x_y_0_0_x_x_yy_xx = buffer_1100_ppdd[342];

    auto g_x_y_0_0_x_x_yy_xy = buffer_1100_ppdd[343];

    auto g_x_y_0_0_x_x_yy_xz = buffer_1100_ppdd[344];

    auto g_x_y_0_0_x_x_yy_yy = buffer_1100_ppdd[345];

    auto g_x_y_0_0_x_x_yy_yz = buffer_1100_ppdd[346];

    auto g_x_y_0_0_x_x_yy_zz = buffer_1100_ppdd[347];

    auto g_x_y_0_0_x_x_yz_xx = buffer_1100_ppdd[348];

    auto g_x_y_0_0_x_x_yz_xy = buffer_1100_ppdd[349];

    auto g_x_y_0_0_x_x_yz_xz = buffer_1100_ppdd[350];

    auto g_x_y_0_0_x_x_yz_yy = buffer_1100_ppdd[351];

    auto g_x_y_0_0_x_x_yz_yz = buffer_1100_ppdd[352];

    auto g_x_y_0_0_x_x_yz_zz = buffer_1100_ppdd[353];

    auto g_x_y_0_0_x_x_zz_xx = buffer_1100_ppdd[354];

    auto g_x_y_0_0_x_x_zz_xy = buffer_1100_ppdd[355];

    auto g_x_y_0_0_x_x_zz_xz = buffer_1100_ppdd[356];

    auto g_x_y_0_0_x_x_zz_yy = buffer_1100_ppdd[357];

    auto g_x_y_0_0_x_x_zz_yz = buffer_1100_ppdd[358];

    auto g_x_y_0_0_x_x_zz_zz = buffer_1100_ppdd[359];

    auto g_x_y_0_0_x_y_xx_xx = buffer_1100_ppdd[360];

    auto g_x_y_0_0_x_y_xx_xy = buffer_1100_ppdd[361];

    auto g_x_y_0_0_x_y_xx_xz = buffer_1100_ppdd[362];

    auto g_x_y_0_0_x_y_xx_yy = buffer_1100_ppdd[363];

    auto g_x_y_0_0_x_y_xx_yz = buffer_1100_ppdd[364];

    auto g_x_y_0_0_x_y_xx_zz = buffer_1100_ppdd[365];

    auto g_x_y_0_0_x_y_xy_xx = buffer_1100_ppdd[366];

    auto g_x_y_0_0_x_y_xy_xy = buffer_1100_ppdd[367];

    auto g_x_y_0_0_x_y_xy_xz = buffer_1100_ppdd[368];

    auto g_x_y_0_0_x_y_xy_yy = buffer_1100_ppdd[369];

    auto g_x_y_0_0_x_y_xy_yz = buffer_1100_ppdd[370];

    auto g_x_y_0_0_x_y_xy_zz = buffer_1100_ppdd[371];

    auto g_x_y_0_0_x_y_xz_xx = buffer_1100_ppdd[372];

    auto g_x_y_0_0_x_y_xz_xy = buffer_1100_ppdd[373];

    auto g_x_y_0_0_x_y_xz_xz = buffer_1100_ppdd[374];

    auto g_x_y_0_0_x_y_xz_yy = buffer_1100_ppdd[375];

    auto g_x_y_0_0_x_y_xz_yz = buffer_1100_ppdd[376];

    auto g_x_y_0_0_x_y_xz_zz = buffer_1100_ppdd[377];

    auto g_x_y_0_0_x_y_yy_xx = buffer_1100_ppdd[378];

    auto g_x_y_0_0_x_y_yy_xy = buffer_1100_ppdd[379];

    auto g_x_y_0_0_x_y_yy_xz = buffer_1100_ppdd[380];

    auto g_x_y_0_0_x_y_yy_yy = buffer_1100_ppdd[381];

    auto g_x_y_0_0_x_y_yy_yz = buffer_1100_ppdd[382];

    auto g_x_y_0_0_x_y_yy_zz = buffer_1100_ppdd[383];

    auto g_x_y_0_0_x_y_yz_xx = buffer_1100_ppdd[384];

    auto g_x_y_0_0_x_y_yz_xy = buffer_1100_ppdd[385];

    auto g_x_y_0_0_x_y_yz_xz = buffer_1100_ppdd[386];

    auto g_x_y_0_0_x_y_yz_yy = buffer_1100_ppdd[387];

    auto g_x_y_0_0_x_y_yz_yz = buffer_1100_ppdd[388];

    auto g_x_y_0_0_x_y_yz_zz = buffer_1100_ppdd[389];

    auto g_x_y_0_0_x_y_zz_xx = buffer_1100_ppdd[390];

    auto g_x_y_0_0_x_y_zz_xy = buffer_1100_ppdd[391];

    auto g_x_y_0_0_x_y_zz_xz = buffer_1100_ppdd[392];

    auto g_x_y_0_0_x_y_zz_yy = buffer_1100_ppdd[393];

    auto g_x_y_0_0_x_y_zz_yz = buffer_1100_ppdd[394];

    auto g_x_y_0_0_x_y_zz_zz = buffer_1100_ppdd[395];

    auto g_x_y_0_0_x_z_xx_xx = buffer_1100_ppdd[396];

    auto g_x_y_0_0_x_z_xx_xy = buffer_1100_ppdd[397];

    auto g_x_y_0_0_x_z_xx_xz = buffer_1100_ppdd[398];

    auto g_x_y_0_0_x_z_xx_yy = buffer_1100_ppdd[399];

    auto g_x_y_0_0_x_z_xx_yz = buffer_1100_ppdd[400];

    auto g_x_y_0_0_x_z_xx_zz = buffer_1100_ppdd[401];

    auto g_x_y_0_0_x_z_xy_xx = buffer_1100_ppdd[402];

    auto g_x_y_0_0_x_z_xy_xy = buffer_1100_ppdd[403];

    auto g_x_y_0_0_x_z_xy_xz = buffer_1100_ppdd[404];

    auto g_x_y_0_0_x_z_xy_yy = buffer_1100_ppdd[405];

    auto g_x_y_0_0_x_z_xy_yz = buffer_1100_ppdd[406];

    auto g_x_y_0_0_x_z_xy_zz = buffer_1100_ppdd[407];

    auto g_x_y_0_0_x_z_xz_xx = buffer_1100_ppdd[408];

    auto g_x_y_0_0_x_z_xz_xy = buffer_1100_ppdd[409];

    auto g_x_y_0_0_x_z_xz_xz = buffer_1100_ppdd[410];

    auto g_x_y_0_0_x_z_xz_yy = buffer_1100_ppdd[411];

    auto g_x_y_0_0_x_z_xz_yz = buffer_1100_ppdd[412];

    auto g_x_y_0_0_x_z_xz_zz = buffer_1100_ppdd[413];

    auto g_x_y_0_0_x_z_yy_xx = buffer_1100_ppdd[414];

    auto g_x_y_0_0_x_z_yy_xy = buffer_1100_ppdd[415];

    auto g_x_y_0_0_x_z_yy_xz = buffer_1100_ppdd[416];

    auto g_x_y_0_0_x_z_yy_yy = buffer_1100_ppdd[417];

    auto g_x_y_0_0_x_z_yy_yz = buffer_1100_ppdd[418];

    auto g_x_y_0_0_x_z_yy_zz = buffer_1100_ppdd[419];

    auto g_x_y_0_0_x_z_yz_xx = buffer_1100_ppdd[420];

    auto g_x_y_0_0_x_z_yz_xy = buffer_1100_ppdd[421];

    auto g_x_y_0_0_x_z_yz_xz = buffer_1100_ppdd[422];

    auto g_x_y_0_0_x_z_yz_yy = buffer_1100_ppdd[423];

    auto g_x_y_0_0_x_z_yz_yz = buffer_1100_ppdd[424];

    auto g_x_y_0_0_x_z_yz_zz = buffer_1100_ppdd[425];

    auto g_x_y_0_0_x_z_zz_xx = buffer_1100_ppdd[426];

    auto g_x_y_0_0_x_z_zz_xy = buffer_1100_ppdd[427];

    auto g_x_y_0_0_x_z_zz_xz = buffer_1100_ppdd[428];

    auto g_x_y_0_0_x_z_zz_yy = buffer_1100_ppdd[429];

    auto g_x_y_0_0_x_z_zz_yz = buffer_1100_ppdd[430];

    auto g_x_y_0_0_x_z_zz_zz = buffer_1100_ppdd[431];

    auto g_x_y_0_0_y_x_xx_xx = buffer_1100_ppdd[432];

    auto g_x_y_0_0_y_x_xx_xy = buffer_1100_ppdd[433];

    auto g_x_y_0_0_y_x_xx_xz = buffer_1100_ppdd[434];

    auto g_x_y_0_0_y_x_xx_yy = buffer_1100_ppdd[435];

    auto g_x_y_0_0_y_x_xx_yz = buffer_1100_ppdd[436];

    auto g_x_y_0_0_y_x_xx_zz = buffer_1100_ppdd[437];

    auto g_x_y_0_0_y_x_xy_xx = buffer_1100_ppdd[438];

    auto g_x_y_0_0_y_x_xy_xy = buffer_1100_ppdd[439];

    auto g_x_y_0_0_y_x_xy_xz = buffer_1100_ppdd[440];

    auto g_x_y_0_0_y_x_xy_yy = buffer_1100_ppdd[441];

    auto g_x_y_0_0_y_x_xy_yz = buffer_1100_ppdd[442];

    auto g_x_y_0_0_y_x_xy_zz = buffer_1100_ppdd[443];

    auto g_x_y_0_0_y_x_xz_xx = buffer_1100_ppdd[444];

    auto g_x_y_0_0_y_x_xz_xy = buffer_1100_ppdd[445];

    auto g_x_y_0_0_y_x_xz_xz = buffer_1100_ppdd[446];

    auto g_x_y_0_0_y_x_xz_yy = buffer_1100_ppdd[447];

    auto g_x_y_0_0_y_x_xz_yz = buffer_1100_ppdd[448];

    auto g_x_y_0_0_y_x_xz_zz = buffer_1100_ppdd[449];

    auto g_x_y_0_0_y_x_yy_xx = buffer_1100_ppdd[450];

    auto g_x_y_0_0_y_x_yy_xy = buffer_1100_ppdd[451];

    auto g_x_y_0_0_y_x_yy_xz = buffer_1100_ppdd[452];

    auto g_x_y_0_0_y_x_yy_yy = buffer_1100_ppdd[453];

    auto g_x_y_0_0_y_x_yy_yz = buffer_1100_ppdd[454];

    auto g_x_y_0_0_y_x_yy_zz = buffer_1100_ppdd[455];

    auto g_x_y_0_0_y_x_yz_xx = buffer_1100_ppdd[456];

    auto g_x_y_0_0_y_x_yz_xy = buffer_1100_ppdd[457];

    auto g_x_y_0_0_y_x_yz_xz = buffer_1100_ppdd[458];

    auto g_x_y_0_0_y_x_yz_yy = buffer_1100_ppdd[459];

    auto g_x_y_0_0_y_x_yz_yz = buffer_1100_ppdd[460];

    auto g_x_y_0_0_y_x_yz_zz = buffer_1100_ppdd[461];

    auto g_x_y_0_0_y_x_zz_xx = buffer_1100_ppdd[462];

    auto g_x_y_0_0_y_x_zz_xy = buffer_1100_ppdd[463];

    auto g_x_y_0_0_y_x_zz_xz = buffer_1100_ppdd[464];

    auto g_x_y_0_0_y_x_zz_yy = buffer_1100_ppdd[465];

    auto g_x_y_0_0_y_x_zz_yz = buffer_1100_ppdd[466];

    auto g_x_y_0_0_y_x_zz_zz = buffer_1100_ppdd[467];

    auto g_x_y_0_0_y_y_xx_xx = buffer_1100_ppdd[468];

    auto g_x_y_0_0_y_y_xx_xy = buffer_1100_ppdd[469];

    auto g_x_y_0_0_y_y_xx_xz = buffer_1100_ppdd[470];

    auto g_x_y_0_0_y_y_xx_yy = buffer_1100_ppdd[471];

    auto g_x_y_0_0_y_y_xx_yz = buffer_1100_ppdd[472];

    auto g_x_y_0_0_y_y_xx_zz = buffer_1100_ppdd[473];

    auto g_x_y_0_0_y_y_xy_xx = buffer_1100_ppdd[474];

    auto g_x_y_0_0_y_y_xy_xy = buffer_1100_ppdd[475];

    auto g_x_y_0_0_y_y_xy_xz = buffer_1100_ppdd[476];

    auto g_x_y_0_0_y_y_xy_yy = buffer_1100_ppdd[477];

    auto g_x_y_0_0_y_y_xy_yz = buffer_1100_ppdd[478];

    auto g_x_y_0_0_y_y_xy_zz = buffer_1100_ppdd[479];

    auto g_x_y_0_0_y_y_xz_xx = buffer_1100_ppdd[480];

    auto g_x_y_0_0_y_y_xz_xy = buffer_1100_ppdd[481];

    auto g_x_y_0_0_y_y_xz_xz = buffer_1100_ppdd[482];

    auto g_x_y_0_0_y_y_xz_yy = buffer_1100_ppdd[483];

    auto g_x_y_0_0_y_y_xz_yz = buffer_1100_ppdd[484];

    auto g_x_y_0_0_y_y_xz_zz = buffer_1100_ppdd[485];

    auto g_x_y_0_0_y_y_yy_xx = buffer_1100_ppdd[486];

    auto g_x_y_0_0_y_y_yy_xy = buffer_1100_ppdd[487];

    auto g_x_y_0_0_y_y_yy_xz = buffer_1100_ppdd[488];

    auto g_x_y_0_0_y_y_yy_yy = buffer_1100_ppdd[489];

    auto g_x_y_0_0_y_y_yy_yz = buffer_1100_ppdd[490];

    auto g_x_y_0_0_y_y_yy_zz = buffer_1100_ppdd[491];

    auto g_x_y_0_0_y_y_yz_xx = buffer_1100_ppdd[492];

    auto g_x_y_0_0_y_y_yz_xy = buffer_1100_ppdd[493];

    auto g_x_y_0_0_y_y_yz_xz = buffer_1100_ppdd[494];

    auto g_x_y_0_0_y_y_yz_yy = buffer_1100_ppdd[495];

    auto g_x_y_0_0_y_y_yz_yz = buffer_1100_ppdd[496];

    auto g_x_y_0_0_y_y_yz_zz = buffer_1100_ppdd[497];

    auto g_x_y_0_0_y_y_zz_xx = buffer_1100_ppdd[498];

    auto g_x_y_0_0_y_y_zz_xy = buffer_1100_ppdd[499];

    auto g_x_y_0_0_y_y_zz_xz = buffer_1100_ppdd[500];

    auto g_x_y_0_0_y_y_zz_yy = buffer_1100_ppdd[501];

    auto g_x_y_0_0_y_y_zz_yz = buffer_1100_ppdd[502];

    auto g_x_y_0_0_y_y_zz_zz = buffer_1100_ppdd[503];

    auto g_x_y_0_0_y_z_xx_xx = buffer_1100_ppdd[504];

    auto g_x_y_0_0_y_z_xx_xy = buffer_1100_ppdd[505];

    auto g_x_y_0_0_y_z_xx_xz = buffer_1100_ppdd[506];

    auto g_x_y_0_0_y_z_xx_yy = buffer_1100_ppdd[507];

    auto g_x_y_0_0_y_z_xx_yz = buffer_1100_ppdd[508];

    auto g_x_y_0_0_y_z_xx_zz = buffer_1100_ppdd[509];

    auto g_x_y_0_0_y_z_xy_xx = buffer_1100_ppdd[510];

    auto g_x_y_0_0_y_z_xy_xy = buffer_1100_ppdd[511];

    auto g_x_y_0_0_y_z_xy_xz = buffer_1100_ppdd[512];

    auto g_x_y_0_0_y_z_xy_yy = buffer_1100_ppdd[513];

    auto g_x_y_0_0_y_z_xy_yz = buffer_1100_ppdd[514];

    auto g_x_y_0_0_y_z_xy_zz = buffer_1100_ppdd[515];

    auto g_x_y_0_0_y_z_xz_xx = buffer_1100_ppdd[516];

    auto g_x_y_0_0_y_z_xz_xy = buffer_1100_ppdd[517];

    auto g_x_y_0_0_y_z_xz_xz = buffer_1100_ppdd[518];

    auto g_x_y_0_0_y_z_xz_yy = buffer_1100_ppdd[519];

    auto g_x_y_0_0_y_z_xz_yz = buffer_1100_ppdd[520];

    auto g_x_y_0_0_y_z_xz_zz = buffer_1100_ppdd[521];

    auto g_x_y_0_0_y_z_yy_xx = buffer_1100_ppdd[522];

    auto g_x_y_0_0_y_z_yy_xy = buffer_1100_ppdd[523];

    auto g_x_y_0_0_y_z_yy_xz = buffer_1100_ppdd[524];

    auto g_x_y_0_0_y_z_yy_yy = buffer_1100_ppdd[525];

    auto g_x_y_0_0_y_z_yy_yz = buffer_1100_ppdd[526];

    auto g_x_y_0_0_y_z_yy_zz = buffer_1100_ppdd[527];

    auto g_x_y_0_0_y_z_yz_xx = buffer_1100_ppdd[528];

    auto g_x_y_0_0_y_z_yz_xy = buffer_1100_ppdd[529];

    auto g_x_y_0_0_y_z_yz_xz = buffer_1100_ppdd[530];

    auto g_x_y_0_0_y_z_yz_yy = buffer_1100_ppdd[531];

    auto g_x_y_0_0_y_z_yz_yz = buffer_1100_ppdd[532];

    auto g_x_y_0_0_y_z_yz_zz = buffer_1100_ppdd[533];

    auto g_x_y_0_0_y_z_zz_xx = buffer_1100_ppdd[534];

    auto g_x_y_0_0_y_z_zz_xy = buffer_1100_ppdd[535];

    auto g_x_y_0_0_y_z_zz_xz = buffer_1100_ppdd[536];

    auto g_x_y_0_0_y_z_zz_yy = buffer_1100_ppdd[537];

    auto g_x_y_0_0_y_z_zz_yz = buffer_1100_ppdd[538];

    auto g_x_y_0_0_y_z_zz_zz = buffer_1100_ppdd[539];

    auto g_x_y_0_0_z_x_xx_xx = buffer_1100_ppdd[540];

    auto g_x_y_0_0_z_x_xx_xy = buffer_1100_ppdd[541];

    auto g_x_y_0_0_z_x_xx_xz = buffer_1100_ppdd[542];

    auto g_x_y_0_0_z_x_xx_yy = buffer_1100_ppdd[543];

    auto g_x_y_0_0_z_x_xx_yz = buffer_1100_ppdd[544];

    auto g_x_y_0_0_z_x_xx_zz = buffer_1100_ppdd[545];

    auto g_x_y_0_0_z_x_xy_xx = buffer_1100_ppdd[546];

    auto g_x_y_0_0_z_x_xy_xy = buffer_1100_ppdd[547];

    auto g_x_y_0_0_z_x_xy_xz = buffer_1100_ppdd[548];

    auto g_x_y_0_0_z_x_xy_yy = buffer_1100_ppdd[549];

    auto g_x_y_0_0_z_x_xy_yz = buffer_1100_ppdd[550];

    auto g_x_y_0_0_z_x_xy_zz = buffer_1100_ppdd[551];

    auto g_x_y_0_0_z_x_xz_xx = buffer_1100_ppdd[552];

    auto g_x_y_0_0_z_x_xz_xy = buffer_1100_ppdd[553];

    auto g_x_y_0_0_z_x_xz_xz = buffer_1100_ppdd[554];

    auto g_x_y_0_0_z_x_xz_yy = buffer_1100_ppdd[555];

    auto g_x_y_0_0_z_x_xz_yz = buffer_1100_ppdd[556];

    auto g_x_y_0_0_z_x_xz_zz = buffer_1100_ppdd[557];

    auto g_x_y_0_0_z_x_yy_xx = buffer_1100_ppdd[558];

    auto g_x_y_0_0_z_x_yy_xy = buffer_1100_ppdd[559];

    auto g_x_y_0_0_z_x_yy_xz = buffer_1100_ppdd[560];

    auto g_x_y_0_0_z_x_yy_yy = buffer_1100_ppdd[561];

    auto g_x_y_0_0_z_x_yy_yz = buffer_1100_ppdd[562];

    auto g_x_y_0_0_z_x_yy_zz = buffer_1100_ppdd[563];

    auto g_x_y_0_0_z_x_yz_xx = buffer_1100_ppdd[564];

    auto g_x_y_0_0_z_x_yz_xy = buffer_1100_ppdd[565];

    auto g_x_y_0_0_z_x_yz_xz = buffer_1100_ppdd[566];

    auto g_x_y_0_0_z_x_yz_yy = buffer_1100_ppdd[567];

    auto g_x_y_0_0_z_x_yz_yz = buffer_1100_ppdd[568];

    auto g_x_y_0_0_z_x_yz_zz = buffer_1100_ppdd[569];

    auto g_x_y_0_0_z_x_zz_xx = buffer_1100_ppdd[570];

    auto g_x_y_0_0_z_x_zz_xy = buffer_1100_ppdd[571];

    auto g_x_y_0_0_z_x_zz_xz = buffer_1100_ppdd[572];

    auto g_x_y_0_0_z_x_zz_yy = buffer_1100_ppdd[573];

    auto g_x_y_0_0_z_x_zz_yz = buffer_1100_ppdd[574];

    auto g_x_y_0_0_z_x_zz_zz = buffer_1100_ppdd[575];

    auto g_x_y_0_0_z_y_xx_xx = buffer_1100_ppdd[576];

    auto g_x_y_0_0_z_y_xx_xy = buffer_1100_ppdd[577];

    auto g_x_y_0_0_z_y_xx_xz = buffer_1100_ppdd[578];

    auto g_x_y_0_0_z_y_xx_yy = buffer_1100_ppdd[579];

    auto g_x_y_0_0_z_y_xx_yz = buffer_1100_ppdd[580];

    auto g_x_y_0_0_z_y_xx_zz = buffer_1100_ppdd[581];

    auto g_x_y_0_0_z_y_xy_xx = buffer_1100_ppdd[582];

    auto g_x_y_0_0_z_y_xy_xy = buffer_1100_ppdd[583];

    auto g_x_y_0_0_z_y_xy_xz = buffer_1100_ppdd[584];

    auto g_x_y_0_0_z_y_xy_yy = buffer_1100_ppdd[585];

    auto g_x_y_0_0_z_y_xy_yz = buffer_1100_ppdd[586];

    auto g_x_y_0_0_z_y_xy_zz = buffer_1100_ppdd[587];

    auto g_x_y_0_0_z_y_xz_xx = buffer_1100_ppdd[588];

    auto g_x_y_0_0_z_y_xz_xy = buffer_1100_ppdd[589];

    auto g_x_y_0_0_z_y_xz_xz = buffer_1100_ppdd[590];

    auto g_x_y_0_0_z_y_xz_yy = buffer_1100_ppdd[591];

    auto g_x_y_0_0_z_y_xz_yz = buffer_1100_ppdd[592];

    auto g_x_y_0_0_z_y_xz_zz = buffer_1100_ppdd[593];

    auto g_x_y_0_0_z_y_yy_xx = buffer_1100_ppdd[594];

    auto g_x_y_0_0_z_y_yy_xy = buffer_1100_ppdd[595];

    auto g_x_y_0_0_z_y_yy_xz = buffer_1100_ppdd[596];

    auto g_x_y_0_0_z_y_yy_yy = buffer_1100_ppdd[597];

    auto g_x_y_0_0_z_y_yy_yz = buffer_1100_ppdd[598];

    auto g_x_y_0_0_z_y_yy_zz = buffer_1100_ppdd[599];

    auto g_x_y_0_0_z_y_yz_xx = buffer_1100_ppdd[600];

    auto g_x_y_0_0_z_y_yz_xy = buffer_1100_ppdd[601];

    auto g_x_y_0_0_z_y_yz_xz = buffer_1100_ppdd[602];

    auto g_x_y_0_0_z_y_yz_yy = buffer_1100_ppdd[603];

    auto g_x_y_0_0_z_y_yz_yz = buffer_1100_ppdd[604];

    auto g_x_y_0_0_z_y_yz_zz = buffer_1100_ppdd[605];

    auto g_x_y_0_0_z_y_zz_xx = buffer_1100_ppdd[606];

    auto g_x_y_0_0_z_y_zz_xy = buffer_1100_ppdd[607];

    auto g_x_y_0_0_z_y_zz_xz = buffer_1100_ppdd[608];

    auto g_x_y_0_0_z_y_zz_yy = buffer_1100_ppdd[609];

    auto g_x_y_0_0_z_y_zz_yz = buffer_1100_ppdd[610];

    auto g_x_y_0_0_z_y_zz_zz = buffer_1100_ppdd[611];

    auto g_x_y_0_0_z_z_xx_xx = buffer_1100_ppdd[612];

    auto g_x_y_0_0_z_z_xx_xy = buffer_1100_ppdd[613];

    auto g_x_y_0_0_z_z_xx_xz = buffer_1100_ppdd[614];

    auto g_x_y_0_0_z_z_xx_yy = buffer_1100_ppdd[615];

    auto g_x_y_0_0_z_z_xx_yz = buffer_1100_ppdd[616];

    auto g_x_y_0_0_z_z_xx_zz = buffer_1100_ppdd[617];

    auto g_x_y_0_0_z_z_xy_xx = buffer_1100_ppdd[618];

    auto g_x_y_0_0_z_z_xy_xy = buffer_1100_ppdd[619];

    auto g_x_y_0_0_z_z_xy_xz = buffer_1100_ppdd[620];

    auto g_x_y_0_0_z_z_xy_yy = buffer_1100_ppdd[621];

    auto g_x_y_0_0_z_z_xy_yz = buffer_1100_ppdd[622];

    auto g_x_y_0_0_z_z_xy_zz = buffer_1100_ppdd[623];

    auto g_x_y_0_0_z_z_xz_xx = buffer_1100_ppdd[624];

    auto g_x_y_0_0_z_z_xz_xy = buffer_1100_ppdd[625];

    auto g_x_y_0_0_z_z_xz_xz = buffer_1100_ppdd[626];

    auto g_x_y_0_0_z_z_xz_yy = buffer_1100_ppdd[627];

    auto g_x_y_0_0_z_z_xz_yz = buffer_1100_ppdd[628];

    auto g_x_y_0_0_z_z_xz_zz = buffer_1100_ppdd[629];

    auto g_x_y_0_0_z_z_yy_xx = buffer_1100_ppdd[630];

    auto g_x_y_0_0_z_z_yy_xy = buffer_1100_ppdd[631];

    auto g_x_y_0_0_z_z_yy_xz = buffer_1100_ppdd[632];

    auto g_x_y_0_0_z_z_yy_yy = buffer_1100_ppdd[633];

    auto g_x_y_0_0_z_z_yy_yz = buffer_1100_ppdd[634];

    auto g_x_y_0_0_z_z_yy_zz = buffer_1100_ppdd[635];

    auto g_x_y_0_0_z_z_yz_xx = buffer_1100_ppdd[636];

    auto g_x_y_0_0_z_z_yz_xy = buffer_1100_ppdd[637];

    auto g_x_y_0_0_z_z_yz_xz = buffer_1100_ppdd[638];

    auto g_x_y_0_0_z_z_yz_yy = buffer_1100_ppdd[639];

    auto g_x_y_0_0_z_z_yz_yz = buffer_1100_ppdd[640];

    auto g_x_y_0_0_z_z_yz_zz = buffer_1100_ppdd[641];

    auto g_x_y_0_0_z_z_zz_xx = buffer_1100_ppdd[642];

    auto g_x_y_0_0_z_z_zz_xy = buffer_1100_ppdd[643];

    auto g_x_y_0_0_z_z_zz_xz = buffer_1100_ppdd[644];

    auto g_x_y_0_0_z_z_zz_yy = buffer_1100_ppdd[645];

    auto g_x_y_0_0_z_z_zz_yz = buffer_1100_ppdd[646];

    auto g_x_y_0_0_z_z_zz_zz = buffer_1100_ppdd[647];

    auto g_x_z_0_0_x_x_xx_xx = buffer_1100_ppdd[648];

    auto g_x_z_0_0_x_x_xx_xy = buffer_1100_ppdd[649];

    auto g_x_z_0_0_x_x_xx_xz = buffer_1100_ppdd[650];

    auto g_x_z_0_0_x_x_xx_yy = buffer_1100_ppdd[651];

    auto g_x_z_0_0_x_x_xx_yz = buffer_1100_ppdd[652];

    auto g_x_z_0_0_x_x_xx_zz = buffer_1100_ppdd[653];

    auto g_x_z_0_0_x_x_xy_xx = buffer_1100_ppdd[654];

    auto g_x_z_0_0_x_x_xy_xy = buffer_1100_ppdd[655];

    auto g_x_z_0_0_x_x_xy_xz = buffer_1100_ppdd[656];

    auto g_x_z_0_0_x_x_xy_yy = buffer_1100_ppdd[657];

    auto g_x_z_0_0_x_x_xy_yz = buffer_1100_ppdd[658];

    auto g_x_z_0_0_x_x_xy_zz = buffer_1100_ppdd[659];

    auto g_x_z_0_0_x_x_xz_xx = buffer_1100_ppdd[660];

    auto g_x_z_0_0_x_x_xz_xy = buffer_1100_ppdd[661];

    auto g_x_z_0_0_x_x_xz_xz = buffer_1100_ppdd[662];

    auto g_x_z_0_0_x_x_xz_yy = buffer_1100_ppdd[663];

    auto g_x_z_0_0_x_x_xz_yz = buffer_1100_ppdd[664];

    auto g_x_z_0_0_x_x_xz_zz = buffer_1100_ppdd[665];

    auto g_x_z_0_0_x_x_yy_xx = buffer_1100_ppdd[666];

    auto g_x_z_0_0_x_x_yy_xy = buffer_1100_ppdd[667];

    auto g_x_z_0_0_x_x_yy_xz = buffer_1100_ppdd[668];

    auto g_x_z_0_0_x_x_yy_yy = buffer_1100_ppdd[669];

    auto g_x_z_0_0_x_x_yy_yz = buffer_1100_ppdd[670];

    auto g_x_z_0_0_x_x_yy_zz = buffer_1100_ppdd[671];

    auto g_x_z_0_0_x_x_yz_xx = buffer_1100_ppdd[672];

    auto g_x_z_0_0_x_x_yz_xy = buffer_1100_ppdd[673];

    auto g_x_z_0_0_x_x_yz_xz = buffer_1100_ppdd[674];

    auto g_x_z_0_0_x_x_yz_yy = buffer_1100_ppdd[675];

    auto g_x_z_0_0_x_x_yz_yz = buffer_1100_ppdd[676];

    auto g_x_z_0_0_x_x_yz_zz = buffer_1100_ppdd[677];

    auto g_x_z_0_0_x_x_zz_xx = buffer_1100_ppdd[678];

    auto g_x_z_0_0_x_x_zz_xy = buffer_1100_ppdd[679];

    auto g_x_z_0_0_x_x_zz_xz = buffer_1100_ppdd[680];

    auto g_x_z_0_0_x_x_zz_yy = buffer_1100_ppdd[681];

    auto g_x_z_0_0_x_x_zz_yz = buffer_1100_ppdd[682];

    auto g_x_z_0_0_x_x_zz_zz = buffer_1100_ppdd[683];

    auto g_x_z_0_0_x_y_xx_xx = buffer_1100_ppdd[684];

    auto g_x_z_0_0_x_y_xx_xy = buffer_1100_ppdd[685];

    auto g_x_z_0_0_x_y_xx_xz = buffer_1100_ppdd[686];

    auto g_x_z_0_0_x_y_xx_yy = buffer_1100_ppdd[687];

    auto g_x_z_0_0_x_y_xx_yz = buffer_1100_ppdd[688];

    auto g_x_z_0_0_x_y_xx_zz = buffer_1100_ppdd[689];

    auto g_x_z_0_0_x_y_xy_xx = buffer_1100_ppdd[690];

    auto g_x_z_0_0_x_y_xy_xy = buffer_1100_ppdd[691];

    auto g_x_z_0_0_x_y_xy_xz = buffer_1100_ppdd[692];

    auto g_x_z_0_0_x_y_xy_yy = buffer_1100_ppdd[693];

    auto g_x_z_0_0_x_y_xy_yz = buffer_1100_ppdd[694];

    auto g_x_z_0_0_x_y_xy_zz = buffer_1100_ppdd[695];

    auto g_x_z_0_0_x_y_xz_xx = buffer_1100_ppdd[696];

    auto g_x_z_0_0_x_y_xz_xy = buffer_1100_ppdd[697];

    auto g_x_z_0_0_x_y_xz_xz = buffer_1100_ppdd[698];

    auto g_x_z_0_0_x_y_xz_yy = buffer_1100_ppdd[699];

    auto g_x_z_0_0_x_y_xz_yz = buffer_1100_ppdd[700];

    auto g_x_z_0_0_x_y_xz_zz = buffer_1100_ppdd[701];

    auto g_x_z_0_0_x_y_yy_xx = buffer_1100_ppdd[702];

    auto g_x_z_0_0_x_y_yy_xy = buffer_1100_ppdd[703];

    auto g_x_z_0_0_x_y_yy_xz = buffer_1100_ppdd[704];

    auto g_x_z_0_0_x_y_yy_yy = buffer_1100_ppdd[705];

    auto g_x_z_0_0_x_y_yy_yz = buffer_1100_ppdd[706];

    auto g_x_z_0_0_x_y_yy_zz = buffer_1100_ppdd[707];

    auto g_x_z_0_0_x_y_yz_xx = buffer_1100_ppdd[708];

    auto g_x_z_0_0_x_y_yz_xy = buffer_1100_ppdd[709];

    auto g_x_z_0_0_x_y_yz_xz = buffer_1100_ppdd[710];

    auto g_x_z_0_0_x_y_yz_yy = buffer_1100_ppdd[711];

    auto g_x_z_0_0_x_y_yz_yz = buffer_1100_ppdd[712];

    auto g_x_z_0_0_x_y_yz_zz = buffer_1100_ppdd[713];

    auto g_x_z_0_0_x_y_zz_xx = buffer_1100_ppdd[714];

    auto g_x_z_0_0_x_y_zz_xy = buffer_1100_ppdd[715];

    auto g_x_z_0_0_x_y_zz_xz = buffer_1100_ppdd[716];

    auto g_x_z_0_0_x_y_zz_yy = buffer_1100_ppdd[717];

    auto g_x_z_0_0_x_y_zz_yz = buffer_1100_ppdd[718];

    auto g_x_z_0_0_x_y_zz_zz = buffer_1100_ppdd[719];

    auto g_x_z_0_0_x_z_xx_xx = buffer_1100_ppdd[720];

    auto g_x_z_0_0_x_z_xx_xy = buffer_1100_ppdd[721];

    auto g_x_z_0_0_x_z_xx_xz = buffer_1100_ppdd[722];

    auto g_x_z_0_0_x_z_xx_yy = buffer_1100_ppdd[723];

    auto g_x_z_0_0_x_z_xx_yz = buffer_1100_ppdd[724];

    auto g_x_z_0_0_x_z_xx_zz = buffer_1100_ppdd[725];

    auto g_x_z_0_0_x_z_xy_xx = buffer_1100_ppdd[726];

    auto g_x_z_0_0_x_z_xy_xy = buffer_1100_ppdd[727];

    auto g_x_z_0_0_x_z_xy_xz = buffer_1100_ppdd[728];

    auto g_x_z_0_0_x_z_xy_yy = buffer_1100_ppdd[729];

    auto g_x_z_0_0_x_z_xy_yz = buffer_1100_ppdd[730];

    auto g_x_z_0_0_x_z_xy_zz = buffer_1100_ppdd[731];

    auto g_x_z_0_0_x_z_xz_xx = buffer_1100_ppdd[732];

    auto g_x_z_0_0_x_z_xz_xy = buffer_1100_ppdd[733];

    auto g_x_z_0_0_x_z_xz_xz = buffer_1100_ppdd[734];

    auto g_x_z_0_0_x_z_xz_yy = buffer_1100_ppdd[735];

    auto g_x_z_0_0_x_z_xz_yz = buffer_1100_ppdd[736];

    auto g_x_z_0_0_x_z_xz_zz = buffer_1100_ppdd[737];

    auto g_x_z_0_0_x_z_yy_xx = buffer_1100_ppdd[738];

    auto g_x_z_0_0_x_z_yy_xy = buffer_1100_ppdd[739];

    auto g_x_z_0_0_x_z_yy_xz = buffer_1100_ppdd[740];

    auto g_x_z_0_0_x_z_yy_yy = buffer_1100_ppdd[741];

    auto g_x_z_0_0_x_z_yy_yz = buffer_1100_ppdd[742];

    auto g_x_z_0_0_x_z_yy_zz = buffer_1100_ppdd[743];

    auto g_x_z_0_0_x_z_yz_xx = buffer_1100_ppdd[744];

    auto g_x_z_0_0_x_z_yz_xy = buffer_1100_ppdd[745];

    auto g_x_z_0_0_x_z_yz_xz = buffer_1100_ppdd[746];

    auto g_x_z_0_0_x_z_yz_yy = buffer_1100_ppdd[747];

    auto g_x_z_0_0_x_z_yz_yz = buffer_1100_ppdd[748];

    auto g_x_z_0_0_x_z_yz_zz = buffer_1100_ppdd[749];

    auto g_x_z_0_0_x_z_zz_xx = buffer_1100_ppdd[750];

    auto g_x_z_0_0_x_z_zz_xy = buffer_1100_ppdd[751];

    auto g_x_z_0_0_x_z_zz_xz = buffer_1100_ppdd[752];

    auto g_x_z_0_0_x_z_zz_yy = buffer_1100_ppdd[753];

    auto g_x_z_0_0_x_z_zz_yz = buffer_1100_ppdd[754];

    auto g_x_z_0_0_x_z_zz_zz = buffer_1100_ppdd[755];

    auto g_x_z_0_0_y_x_xx_xx = buffer_1100_ppdd[756];

    auto g_x_z_0_0_y_x_xx_xy = buffer_1100_ppdd[757];

    auto g_x_z_0_0_y_x_xx_xz = buffer_1100_ppdd[758];

    auto g_x_z_0_0_y_x_xx_yy = buffer_1100_ppdd[759];

    auto g_x_z_0_0_y_x_xx_yz = buffer_1100_ppdd[760];

    auto g_x_z_0_0_y_x_xx_zz = buffer_1100_ppdd[761];

    auto g_x_z_0_0_y_x_xy_xx = buffer_1100_ppdd[762];

    auto g_x_z_0_0_y_x_xy_xy = buffer_1100_ppdd[763];

    auto g_x_z_0_0_y_x_xy_xz = buffer_1100_ppdd[764];

    auto g_x_z_0_0_y_x_xy_yy = buffer_1100_ppdd[765];

    auto g_x_z_0_0_y_x_xy_yz = buffer_1100_ppdd[766];

    auto g_x_z_0_0_y_x_xy_zz = buffer_1100_ppdd[767];

    auto g_x_z_0_0_y_x_xz_xx = buffer_1100_ppdd[768];

    auto g_x_z_0_0_y_x_xz_xy = buffer_1100_ppdd[769];

    auto g_x_z_0_0_y_x_xz_xz = buffer_1100_ppdd[770];

    auto g_x_z_0_0_y_x_xz_yy = buffer_1100_ppdd[771];

    auto g_x_z_0_0_y_x_xz_yz = buffer_1100_ppdd[772];

    auto g_x_z_0_0_y_x_xz_zz = buffer_1100_ppdd[773];

    auto g_x_z_0_0_y_x_yy_xx = buffer_1100_ppdd[774];

    auto g_x_z_0_0_y_x_yy_xy = buffer_1100_ppdd[775];

    auto g_x_z_0_0_y_x_yy_xz = buffer_1100_ppdd[776];

    auto g_x_z_0_0_y_x_yy_yy = buffer_1100_ppdd[777];

    auto g_x_z_0_0_y_x_yy_yz = buffer_1100_ppdd[778];

    auto g_x_z_0_0_y_x_yy_zz = buffer_1100_ppdd[779];

    auto g_x_z_0_0_y_x_yz_xx = buffer_1100_ppdd[780];

    auto g_x_z_0_0_y_x_yz_xy = buffer_1100_ppdd[781];

    auto g_x_z_0_0_y_x_yz_xz = buffer_1100_ppdd[782];

    auto g_x_z_0_0_y_x_yz_yy = buffer_1100_ppdd[783];

    auto g_x_z_0_0_y_x_yz_yz = buffer_1100_ppdd[784];

    auto g_x_z_0_0_y_x_yz_zz = buffer_1100_ppdd[785];

    auto g_x_z_0_0_y_x_zz_xx = buffer_1100_ppdd[786];

    auto g_x_z_0_0_y_x_zz_xy = buffer_1100_ppdd[787];

    auto g_x_z_0_0_y_x_zz_xz = buffer_1100_ppdd[788];

    auto g_x_z_0_0_y_x_zz_yy = buffer_1100_ppdd[789];

    auto g_x_z_0_0_y_x_zz_yz = buffer_1100_ppdd[790];

    auto g_x_z_0_0_y_x_zz_zz = buffer_1100_ppdd[791];

    auto g_x_z_0_0_y_y_xx_xx = buffer_1100_ppdd[792];

    auto g_x_z_0_0_y_y_xx_xy = buffer_1100_ppdd[793];

    auto g_x_z_0_0_y_y_xx_xz = buffer_1100_ppdd[794];

    auto g_x_z_0_0_y_y_xx_yy = buffer_1100_ppdd[795];

    auto g_x_z_0_0_y_y_xx_yz = buffer_1100_ppdd[796];

    auto g_x_z_0_0_y_y_xx_zz = buffer_1100_ppdd[797];

    auto g_x_z_0_0_y_y_xy_xx = buffer_1100_ppdd[798];

    auto g_x_z_0_0_y_y_xy_xy = buffer_1100_ppdd[799];

    auto g_x_z_0_0_y_y_xy_xz = buffer_1100_ppdd[800];

    auto g_x_z_0_0_y_y_xy_yy = buffer_1100_ppdd[801];

    auto g_x_z_0_0_y_y_xy_yz = buffer_1100_ppdd[802];

    auto g_x_z_0_0_y_y_xy_zz = buffer_1100_ppdd[803];

    auto g_x_z_0_0_y_y_xz_xx = buffer_1100_ppdd[804];

    auto g_x_z_0_0_y_y_xz_xy = buffer_1100_ppdd[805];

    auto g_x_z_0_0_y_y_xz_xz = buffer_1100_ppdd[806];

    auto g_x_z_0_0_y_y_xz_yy = buffer_1100_ppdd[807];

    auto g_x_z_0_0_y_y_xz_yz = buffer_1100_ppdd[808];

    auto g_x_z_0_0_y_y_xz_zz = buffer_1100_ppdd[809];

    auto g_x_z_0_0_y_y_yy_xx = buffer_1100_ppdd[810];

    auto g_x_z_0_0_y_y_yy_xy = buffer_1100_ppdd[811];

    auto g_x_z_0_0_y_y_yy_xz = buffer_1100_ppdd[812];

    auto g_x_z_0_0_y_y_yy_yy = buffer_1100_ppdd[813];

    auto g_x_z_0_0_y_y_yy_yz = buffer_1100_ppdd[814];

    auto g_x_z_0_0_y_y_yy_zz = buffer_1100_ppdd[815];

    auto g_x_z_0_0_y_y_yz_xx = buffer_1100_ppdd[816];

    auto g_x_z_0_0_y_y_yz_xy = buffer_1100_ppdd[817];

    auto g_x_z_0_0_y_y_yz_xz = buffer_1100_ppdd[818];

    auto g_x_z_0_0_y_y_yz_yy = buffer_1100_ppdd[819];

    auto g_x_z_0_0_y_y_yz_yz = buffer_1100_ppdd[820];

    auto g_x_z_0_0_y_y_yz_zz = buffer_1100_ppdd[821];

    auto g_x_z_0_0_y_y_zz_xx = buffer_1100_ppdd[822];

    auto g_x_z_0_0_y_y_zz_xy = buffer_1100_ppdd[823];

    auto g_x_z_0_0_y_y_zz_xz = buffer_1100_ppdd[824];

    auto g_x_z_0_0_y_y_zz_yy = buffer_1100_ppdd[825];

    auto g_x_z_0_0_y_y_zz_yz = buffer_1100_ppdd[826];

    auto g_x_z_0_0_y_y_zz_zz = buffer_1100_ppdd[827];

    auto g_x_z_0_0_y_z_xx_xx = buffer_1100_ppdd[828];

    auto g_x_z_0_0_y_z_xx_xy = buffer_1100_ppdd[829];

    auto g_x_z_0_0_y_z_xx_xz = buffer_1100_ppdd[830];

    auto g_x_z_0_0_y_z_xx_yy = buffer_1100_ppdd[831];

    auto g_x_z_0_0_y_z_xx_yz = buffer_1100_ppdd[832];

    auto g_x_z_0_0_y_z_xx_zz = buffer_1100_ppdd[833];

    auto g_x_z_0_0_y_z_xy_xx = buffer_1100_ppdd[834];

    auto g_x_z_0_0_y_z_xy_xy = buffer_1100_ppdd[835];

    auto g_x_z_0_0_y_z_xy_xz = buffer_1100_ppdd[836];

    auto g_x_z_0_0_y_z_xy_yy = buffer_1100_ppdd[837];

    auto g_x_z_0_0_y_z_xy_yz = buffer_1100_ppdd[838];

    auto g_x_z_0_0_y_z_xy_zz = buffer_1100_ppdd[839];

    auto g_x_z_0_0_y_z_xz_xx = buffer_1100_ppdd[840];

    auto g_x_z_0_0_y_z_xz_xy = buffer_1100_ppdd[841];

    auto g_x_z_0_0_y_z_xz_xz = buffer_1100_ppdd[842];

    auto g_x_z_0_0_y_z_xz_yy = buffer_1100_ppdd[843];

    auto g_x_z_0_0_y_z_xz_yz = buffer_1100_ppdd[844];

    auto g_x_z_0_0_y_z_xz_zz = buffer_1100_ppdd[845];

    auto g_x_z_0_0_y_z_yy_xx = buffer_1100_ppdd[846];

    auto g_x_z_0_0_y_z_yy_xy = buffer_1100_ppdd[847];

    auto g_x_z_0_0_y_z_yy_xz = buffer_1100_ppdd[848];

    auto g_x_z_0_0_y_z_yy_yy = buffer_1100_ppdd[849];

    auto g_x_z_0_0_y_z_yy_yz = buffer_1100_ppdd[850];

    auto g_x_z_0_0_y_z_yy_zz = buffer_1100_ppdd[851];

    auto g_x_z_0_0_y_z_yz_xx = buffer_1100_ppdd[852];

    auto g_x_z_0_0_y_z_yz_xy = buffer_1100_ppdd[853];

    auto g_x_z_0_0_y_z_yz_xz = buffer_1100_ppdd[854];

    auto g_x_z_0_0_y_z_yz_yy = buffer_1100_ppdd[855];

    auto g_x_z_0_0_y_z_yz_yz = buffer_1100_ppdd[856];

    auto g_x_z_0_0_y_z_yz_zz = buffer_1100_ppdd[857];

    auto g_x_z_0_0_y_z_zz_xx = buffer_1100_ppdd[858];

    auto g_x_z_0_0_y_z_zz_xy = buffer_1100_ppdd[859];

    auto g_x_z_0_0_y_z_zz_xz = buffer_1100_ppdd[860];

    auto g_x_z_0_0_y_z_zz_yy = buffer_1100_ppdd[861];

    auto g_x_z_0_0_y_z_zz_yz = buffer_1100_ppdd[862];

    auto g_x_z_0_0_y_z_zz_zz = buffer_1100_ppdd[863];

    auto g_x_z_0_0_z_x_xx_xx = buffer_1100_ppdd[864];

    auto g_x_z_0_0_z_x_xx_xy = buffer_1100_ppdd[865];

    auto g_x_z_0_0_z_x_xx_xz = buffer_1100_ppdd[866];

    auto g_x_z_0_0_z_x_xx_yy = buffer_1100_ppdd[867];

    auto g_x_z_0_0_z_x_xx_yz = buffer_1100_ppdd[868];

    auto g_x_z_0_0_z_x_xx_zz = buffer_1100_ppdd[869];

    auto g_x_z_0_0_z_x_xy_xx = buffer_1100_ppdd[870];

    auto g_x_z_0_0_z_x_xy_xy = buffer_1100_ppdd[871];

    auto g_x_z_0_0_z_x_xy_xz = buffer_1100_ppdd[872];

    auto g_x_z_0_0_z_x_xy_yy = buffer_1100_ppdd[873];

    auto g_x_z_0_0_z_x_xy_yz = buffer_1100_ppdd[874];

    auto g_x_z_0_0_z_x_xy_zz = buffer_1100_ppdd[875];

    auto g_x_z_0_0_z_x_xz_xx = buffer_1100_ppdd[876];

    auto g_x_z_0_0_z_x_xz_xy = buffer_1100_ppdd[877];

    auto g_x_z_0_0_z_x_xz_xz = buffer_1100_ppdd[878];

    auto g_x_z_0_0_z_x_xz_yy = buffer_1100_ppdd[879];

    auto g_x_z_0_0_z_x_xz_yz = buffer_1100_ppdd[880];

    auto g_x_z_0_0_z_x_xz_zz = buffer_1100_ppdd[881];

    auto g_x_z_0_0_z_x_yy_xx = buffer_1100_ppdd[882];

    auto g_x_z_0_0_z_x_yy_xy = buffer_1100_ppdd[883];

    auto g_x_z_0_0_z_x_yy_xz = buffer_1100_ppdd[884];

    auto g_x_z_0_0_z_x_yy_yy = buffer_1100_ppdd[885];

    auto g_x_z_0_0_z_x_yy_yz = buffer_1100_ppdd[886];

    auto g_x_z_0_0_z_x_yy_zz = buffer_1100_ppdd[887];

    auto g_x_z_0_0_z_x_yz_xx = buffer_1100_ppdd[888];

    auto g_x_z_0_0_z_x_yz_xy = buffer_1100_ppdd[889];

    auto g_x_z_0_0_z_x_yz_xz = buffer_1100_ppdd[890];

    auto g_x_z_0_0_z_x_yz_yy = buffer_1100_ppdd[891];

    auto g_x_z_0_0_z_x_yz_yz = buffer_1100_ppdd[892];

    auto g_x_z_0_0_z_x_yz_zz = buffer_1100_ppdd[893];

    auto g_x_z_0_0_z_x_zz_xx = buffer_1100_ppdd[894];

    auto g_x_z_0_0_z_x_zz_xy = buffer_1100_ppdd[895];

    auto g_x_z_0_0_z_x_zz_xz = buffer_1100_ppdd[896];

    auto g_x_z_0_0_z_x_zz_yy = buffer_1100_ppdd[897];

    auto g_x_z_0_0_z_x_zz_yz = buffer_1100_ppdd[898];

    auto g_x_z_0_0_z_x_zz_zz = buffer_1100_ppdd[899];

    auto g_x_z_0_0_z_y_xx_xx = buffer_1100_ppdd[900];

    auto g_x_z_0_0_z_y_xx_xy = buffer_1100_ppdd[901];

    auto g_x_z_0_0_z_y_xx_xz = buffer_1100_ppdd[902];

    auto g_x_z_0_0_z_y_xx_yy = buffer_1100_ppdd[903];

    auto g_x_z_0_0_z_y_xx_yz = buffer_1100_ppdd[904];

    auto g_x_z_0_0_z_y_xx_zz = buffer_1100_ppdd[905];

    auto g_x_z_0_0_z_y_xy_xx = buffer_1100_ppdd[906];

    auto g_x_z_0_0_z_y_xy_xy = buffer_1100_ppdd[907];

    auto g_x_z_0_0_z_y_xy_xz = buffer_1100_ppdd[908];

    auto g_x_z_0_0_z_y_xy_yy = buffer_1100_ppdd[909];

    auto g_x_z_0_0_z_y_xy_yz = buffer_1100_ppdd[910];

    auto g_x_z_0_0_z_y_xy_zz = buffer_1100_ppdd[911];

    auto g_x_z_0_0_z_y_xz_xx = buffer_1100_ppdd[912];

    auto g_x_z_0_0_z_y_xz_xy = buffer_1100_ppdd[913];

    auto g_x_z_0_0_z_y_xz_xz = buffer_1100_ppdd[914];

    auto g_x_z_0_0_z_y_xz_yy = buffer_1100_ppdd[915];

    auto g_x_z_0_0_z_y_xz_yz = buffer_1100_ppdd[916];

    auto g_x_z_0_0_z_y_xz_zz = buffer_1100_ppdd[917];

    auto g_x_z_0_0_z_y_yy_xx = buffer_1100_ppdd[918];

    auto g_x_z_0_0_z_y_yy_xy = buffer_1100_ppdd[919];

    auto g_x_z_0_0_z_y_yy_xz = buffer_1100_ppdd[920];

    auto g_x_z_0_0_z_y_yy_yy = buffer_1100_ppdd[921];

    auto g_x_z_0_0_z_y_yy_yz = buffer_1100_ppdd[922];

    auto g_x_z_0_0_z_y_yy_zz = buffer_1100_ppdd[923];

    auto g_x_z_0_0_z_y_yz_xx = buffer_1100_ppdd[924];

    auto g_x_z_0_0_z_y_yz_xy = buffer_1100_ppdd[925];

    auto g_x_z_0_0_z_y_yz_xz = buffer_1100_ppdd[926];

    auto g_x_z_0_0_z_y_yz_yy = buffer_1100_ppdd[927];

    auto g_x_z_0_0_z_y_yz_yz = buffer_1100_ppdd[928];

    auto g_x_z_0_0_z_y_yz_zz = buffer_1100_ppdd[929];

    auto g_x_z_0_0_z_y_zz_xx = buffer_1100_ppdd[930];

    auto g_x_z_0_0_z_y_zz_xy = buffer_1100_ppdd[931];

    auto g_x_z_0_0_z_y_zz_xz = buffer_1100_ppdd[932];

    auto g_x_z_0_0_z_y_zz_yy = buffer_1100_ppdd[933];

    auto g_x_z_0_0_z_y_zz_yz = buffer_1100_ppdd[934];

    auto g_x_z_0_0_z_y_zz_zz = buffer_1100_ppdd[935];

    auto g_x_z_0_0_z_z_xx_xx = buffer_1100_ppdd[936];

    auto g_x_z_0_0_z_z_xx_xy = buffer_1100_ppdd[937];

    auto g_x_z_0_0_z_z_xx_xz = buffer_1100_ppdd[938];

    auto g_x_z_0_0_z_z_xx_yy = buffer_1100_ppdd[939];

    auto g_x_z_0_0_z_z_xx_yz = buffer_1100_ppdd[940];

    auto g_x_z_0_0_z_z_xx_zz = buffer_1100_ppdd[941];

    auto g_x_z_0_0_z_z_xy_xx = buffer_1100_ppdd[942];

    auto g_x_z_0_0_z_z_xy_xy = buffer_1100_ppdd[943];

    auto g_x_z_0_0_z_z_xy_xz = buffer_1100_ppdd[944];

    auto g_x_z_0_0_z_z_xy_yy = buffer_1100_ppdd[945];

    auto g_x_z_0_0_z_z_xy_yz = buffer_1100_ppdd[946];

    auto g_x_z_0_0_z_z_xy_zz = buffer_1100_ppdd[947];

    auto g_x_z_0_0_z_z_xz_xx = buffer_1100_ppdd[948];

    auto g_x_z_0_0_z_z_xz_xy = buffer_1100_ppdd[949];

    auto g_x_z_0_0_z_z_xz_xz = buffer_1100_ppdd[950];

    auto g_x_z_0_0_z_z_xz_yy = buffer_1100_ppdd[951];

    auto g_x_z_0_0_z_z_xz_yz = buffer_1100_ppdd[952];

    auto g_x_z_0_0_z_z_xz_zz = buffer_1100_ppdd[953];

    auto g_x_z_0_0_z_z_yy_xx = buffer_1100_ppdd[954];

    auto g_x_z_0_0_z_z_yy_xy = buffer_1100_ppdd[955];

    auto g_x_z_0_0_z_z_yy_xz = buffer_1100_ppdd[956];

    auto g_x_z_0_0_z_z_yy_yy = buffer_1100_ppdd[957];

    auto g_x_z_0_0_z_z_yy_yz = buffer_1100_ppdd[958];

    auto g_x_z_0_0_z_z_yy_zz = buffer_1100_ppdd[959];

    auto g_x_z_0_0_z_z_yz_xx = buffer_1100_ppdd[960];

    auto g_x_z_0_0_z_z_yz_xy = buffer_1100_ppdd[961];

    auto g_x_z_0_0_z_z_yz_xz = buffer_1100_ppdd[962];

    auto g_x_z_0_0_z_z_yz_yy = buffer_1100_ppdd[963];

    auto g_x_z_0_0_z_z_yz_yz = buffer_1100_ppdd[964];

    auto g_x_z_0_0_z_z_yz_zz = buffer_1100_ppdd[965];

    auto g_x_z_0_0_z_z_zz_xx = buffer_1100_ppdd[966];

    auto g_x_z_0_0_z_z_zz_xy = buffer_1100_ppdd[967];

    auto g_x_z_0_0_z_z_zz_xz = buffer_1100_ppdd[968];

    auto g_x_z_0_0_z_z_zz_yy = buffer_1100_ppdd[969];

    auto g_x_z_0_0_z_z_zz_yz = buffer_1100_ppdd[970];

    auto g_x_z_0_0_z_z_zz_zz = buffer_1100_ppdd[971];

    auto g_y_x_0_0_x_x_xx_xx = buffer_1100_ppdd[972];

    auto g_y_x_0_0_x_x_xx_xy = buffer_1100_ppdd[973];

    auto g_y_x_0_0_x_x_xx_xz = buffer_1100_ppdd[974];

    auto g_y_x_0_0_x_x_xx_yy = buffer_1100_ppdd[975];

    auto g_y_x_0_0_x_x_xx_yz = buffer_1100_ppdd[976];

    auto g_y_x_0_0_x_x_xx_zz = buffer_1100_ppdd[977];

    auto g_y_x_0_0_x_x_xy_xx = buffer_1100_ppdd[978];

    auto g_y_x_0_0_x_x_xy_xy = buffer_1100_ppdd[979];

    auto g_y_x_0_0_x_x_xy_xz = buffer_1100_ppdd[980];

    auto g_y_x_0_0_x_x_xy_yy = buffer_1100_ppdd[981];

    auto g_y_x_0_0_x_x_xy_yz = buffer_1100_ppdd[982];

    auto g_y_x_0_0_x_x_xy_zz = buffer_1100_ppdd[983];

    auto g_y_x_0_0_x_x_xz_xx = buffer_1100_ppdd[984];

    auto g_y_x_0_0_x_x_xz_xy = buffer_1100_ppdd[985];

    auto g_y_x_0_0_x_x_xz_xz = buffer_1100_ppdd[986];

    auto g_y_x_0_0_x_x_xz_yy = buffer_1100_ppdd[987];

    auto g_y_x_0_0_x_x_xz_yz = buffer_1100_ppdd[988];

    auto g_y_x_0_0_x_x_xz_zz = buffer_1100_ppdd[989];

    auto g_y_x_0_0_x_x_yy_xx = buffer_1100_ppdd[990];

    auto g_y_x_0_0_x_x_yy_xy = buffer_1100_ppdd[991];

    auto g_y_x_0_0_x_x_yy_xz = buffer_1100_ppdd[992];

    auto g_y_x_0_0_x_x_yy_yy = buffer_1100_ppdd[993];

    auto g_y_x_0_0_x_x_yy_yz = buffer_1100_ppdd[994];

    auto g_y_x_0_0_x_x_yy_zz = buffer_1100_ppdd[995];

    auto g_y_x_0_0_x_x_yz_xx = buffer_1100_ppdd[996];

    auto g_y_x_0_0_x_x_yz_xy = buffer_1100_ppdd[997];

    auto g_y_x_0_0_x_x_yz_xz = buffer_1100_ppdd[998];

    auto g_y_x_0_0_x_x_yz_yy = buffer_1100_ppdd[999];

    auto g_y_x_0_0_x_x_yz_yz = buffer_1100_ppdd[1000];

    auto g_y_x_0_0_x_x_yz_zz = buffer_1100_ppdd[1001];

    auto g_y_x_0_0_x_x_zz_xx = buffer_1100_ppdd[1002];

    auto g_y_x_0_0_x_x_zz_xy = buffer_1100_ppdd[1003];

    auto g_y_x_0_0_x_x_zz_xz = buffer_1100_ppdd[1004];

    auto g_y_x_0_0_x_x_zz_yy = buffer_1100_ppdd[1005];

    auto g_y_x_0_0_x_x_zz_yz = buffer_1100_ppdd[1006];

    auto g_y_x_0_0_x_x_zz_zz = buffer_1100_ppdd[1007];

    auto g_y_x_0_0_x_y_xx_xx = buffer_1100_ppdd[1008];

    auto g_y_x_0_0_x_y_xx_xy = buffer_1100_ppdd[1009];

    auto g_y_x_0_0_x_y_xx_xz = buffer_1100_ppdd[1010];

    auto g_y_x_0_0_x_y_xx_yy = buffer_1100_ppdd[1011];

    auto g_y_x_0_0_x_y_xx_yz = buffer_1100_ppdd[1012];

    auto g_y_x_0_0_x_y_xx_zz = buffer_1100_ppdd[1013];

    auto g_y_x_0_0_x_y_xy_xx = buffer_1100_ppdd[1014];

    auto g_y_x_0_0_x_y_xy_xy = buffer_1100_ppdd[1015];

    auto g_y_x_0_0_x_y_xy_xz = buffer_1100_ppdd[1016];

    auto g_y_x_0_0_x_y_xy_yy = buffer_1100_ppdd[1017];

    auto g_y_x_0_0_x_y_xy_yz = buffer_1100_ppdd[1018];

    auto g_y_x_0_0_x_y_xy_zz = buffer_1100_ppdd[1019];

    auto g_y_x_0_0_x_y_xz_xx = buffer_1100_ppdd[1020];

    auto g_y_x_0_0_x_y_xz_xy = buffer_1100_ppdd[1021];

    auto g_y_x_0_0_x_y_xz_xz = buffer_1100_ppdd[1022];

    auto g_y_x_0_0_x_y_xz_yy = buffer_1100_ppdd[1023];

    auto g_y_x_0_0_x_y_xz_yz = buffer_1100_ppdd[1024];

    auto g_y_x_0_0_x_y_xz_zz = buffer_1100_ppdd[1025];

    auto g_y_x_0_0_x_y_yy_xx = buffer_1100_ppdd[1026];

    auto g_y_x_0_0_x_y_yy_xy = buffer_1100_ppdd[1027];

    auto g_y_x_0_0_x_y_yy_xz = buffer_1100_ppdd[1028];

    auto g_y_x_0_0_x_y_yy_yy = buffer_1100_ppdd[1029];

    auto g_y_x_0_0_x_y_yy_yz = buffer_1100_ppdd[1030];

    auto g_y_x_0_0_x_y_yy_zz = buffer_1100_ppdd[1031];

    auto g_y_x_0_0_x_y_yz_xx = buffer_1100_ppdd[1032];

    auto g_y_x_0_0_x_y_yz_xy = buffer_1100_ppdd[1033];

    auto g_y_x_0_0_x_y_yz_xz = buffer_1100_ppdd[1034];

    auto g_y_x_0_0_x_y_yz_yy = buffer_1100_ppdd[1035];

    auto g_y_x_0_0_x_y_yz_yz = buffer_1100_ppdd[1036];

    auto g_y_x_0_0_x_y_yz_zz = buffer_1100_ppdd[1037];

    auto g_y_x_0_0_x_y_zz_xx = buffer_1100_ppdd[1038];

    auto g_y_x_0_0_x_y_zz_xy = buffer_1100_ppdd[1039];

    auto g_y_x_0_0_x_y_zz_xz = buffer_1100_ppdd[1040];

    auto g_y_x_0_0_x_y_zz_yy = buffer_1100_ppdd[1041];

    auto g_y_x_0_0_x_y_zz_yz = buffer_1100_ppdd[1042];

    auto g_y_x_0_0_x_y_zz_zz = buffer_1100_ppdd[1043];

    auto g_y_x_0_0_x_z_xx_xx = buffer_1100_ppdd[1044];

    auto g_y_x_0_0_x_z_xx_xy = buffer_1100_ppdd[1045];

    auto g_y_x_0_0_x_z_xx_xz = buffer_1100_ppdd[1046];

    auto g_y_x_0_0_x_z_xx_yy = buffer_1100_ppdd[1047];

    auto g_y_x_0_0_x_z_xx_yz = buffer_1100_ppdd[1048];

    auto g_y_x_0_0_x_z_xx_zz = buffer_1100_ppdd[1049];

    auto g_y_x_0_0_x_z_xy_xx = buffer_1100_ppdd[1050];

    auto g_y_x_0_0_x_z_xy_xy = buffer_1100_ppdd[1051];

    auto g_y_x_0_0_x_z_xy_xz = buffer_1100_ppdd[1052];

    auto g_y_x_0_0_x_z_xy_yy = buffer_1100_ppdd[1053];

    auto g_y_x_0_0_x_z_xy_yz = buffer_1100_ppdd[1054];

    auto g_y_x_0_0_x_z_xy_zz = buffer_1100_ppdd[1055];

    auto g_y_x_0_0_x_z_xz_xx = buffer_1100_ppdd[1056];

    auto g_y_x_0_0_x_z_xz_xy = buffer_1100_ppdd[1057];

    auto g_y_x_0_0_x_z_xz_xz = buffer_1100_ppdd[1058];

    auto g_y_x_0_0_x_z_xz_yy = buffer_1100_ppdd[1059];

    auto g_y_x_0_0_x_z_xz_yz = buffer_1100_ppdd[1060];

    auto g_y_x_0_0_x_z_xz_zz = buffer_1100_ppdd[1061];

    auto g_y_x_0_0_x_z_yy_xx = buffer_1100_ppdd[1062];

    auto g_y_x_0_0_x_z_yy_xy = buffer_1100_ppdd[1063];

    auto g_y_x_0_0_x_z_yy_xz = buffer_1100_ppdd[1064];

    auto g_y_x_0_0_x_z_yy_yy = buffer_1100_ppdd[1065];

    auto g_y_x_0_0_x_z_yy_yz = buffer_1100_ppdd[1066];

    auto g_y_x_0_0_x_z_yy_zz = buffer_1100_ppdd[1067];

    auto g_y_x_0_0_x_z_yz_xx = buffer_1100_ppdd[1068];

    auto g_y_x_0_0_x_z_yz_xy = buffer_1100_ppdd[1069];

    auto g_y_x_0_0_x_z_yz_xz = buffer_1100_ppdd[1070];

    auto g_y_x_0_0_x_z_yz_yy = buffer_1100_ppdd[1071];

    auto g_y_x_0_0_x_z_yz_yz = buffer_1100_ppdd[1072];

    auto g_y_x_0_0_x_z_yz_zz = buffer_1100_ppdd[1073];

    auto g_y_x_0_0_x_z_zz_xx = buffer_1100_ppdd[1074];

    auto g_y_x_0_0_x_z_zz_xy = buffer_1100_ppdd[1075];

    auto g_y_x_0_0_x_z_zz_xz = buffer_1100_ppdd[1076];

    auto g_y_x_0_0_x_z_zz_yy = buffer_1100_ppdd[1077];

    auto g_y_x_0_0_x_z_zz_yz = buffer_1100_ppdd[1078];

    auto g_y_x_0_0_x_z_zz_zz = buffer_1100_ppdd[1079];

    auto g_y_x_0_0_y_x_xx_xx = buffer_1100_ppdd[1080];

    auto g_y_x_0_0_y_x_xx_xy = buffer_1100_ppdd[1081];

    auto g_y_x_0_0_y_x_xx_xz = buffer_1100_ppdd[1082];

    auto g_y_x_0_0_y_x_xx_yy = buffer_1100_ppdd[1083];

    auto g_y_x_0_0_y_x_xx_yz = buffer_1100_ppdd[1084];

    auto g_y_x_0_0_y_x_xx_zz = buffer_1100_ppdd[1085];

    auto g_y_x_0_0_y_x_xy_xx = buffer_1100_ppdd[1086];

    auto g_y_x_0_0_y_x_xy_xy = buffer_1100_ppdd[1087];

    auto g_y_x_0_0_y_x_xy_xz = buffer_1100_ppdd[1088];

    auto g_y_x_0_0_y_x_xy_yy = buffer_1100_ppdd[1089];

    auto g_y_x_0_0_y_x_xy_yz = buffer_1100_ppdd[1090];

    auto g_y_x_0_0_y_x_xy_zz = buffer_1100_ppdd[1091];

    auto g_y_x_0_0_y_x_xz_xx = buffer_1100_ppdd[1092];

    auto g_y_x_0_0_y_x_xz_xy = buffer_1100_ppdd[1093];

    auto g_y_x_0_0_y_x_xz_xz = buffer_1100_ppdd[1094];

    auto g_y_x_0_0_y_x_xz_yy = buffer_1100_ppdd[1095];

    auto g_y_x_0_0_y_x_xz_yz = buffer_1100_ppdd[1096];

    auto g_y_x_0_0_y_x_xz_zz = buffer_1100_ppdd[1097];

    auto g_y_x_0_0_y_x_yy_xx = buffer_1100_ppdd[1098];

    auto g_y_x_0_0_y_x_yy_xy = buffer_1100_ppdd[1099];

    auto g_y_x_0_0_y_x_yy_xz = buffer_1100_ppdd[1100];

    auto g_y_x_0_0_y_x_yy_yy = buffer_1100_ppdd[1101];

    auto g_y_x_0_0_y_x_yy_yz = buffer_1100_ppdd[1102];

    auto g_y_x_0_0_y_x_yy_zz = buffer_1100_ppdd[1103];

    auto g_y_x_0_0_y_x_yz_xx = buffer_1100_ppdd[1104];

    auto g_y_x_0_0_y_x_yz_xy = buffer_1100_ppdd[1105];

    auto g_y_x_0_0_y_x_yz_xz = buffer_1100_ppdd[1106];

    auto g_y_x_0_0_y_x_yz_yy = buffer_1100_ppdd[1107];

    auto g_y_x_0_0_y_x_yz_yz = buffer_1100_ppdd[1108];

    auto g_y_x_0_0_y_x_yz_zz = buffer_1100_ppdd[1109];

    auto g_y_x_0_0_y_x_zz_xx = buffer_1100_ppdd[1110];

    auto g_y_x_0_0_y_x_zz_xy = buffer_1100_ppdd[1111];

    auto g_y_x_0_0_y_x_zz_xz = buffer_1100_ppdd[1112];

    auto g_y_x_0_0_y_x_zz_yy = buffer_1100_ppdd[1113];

    auto g_y_x_0_0_y_x_zz_yz = buffer_1100_ppdd[1114];

    auto g_y_x_0_0_y_x_zz_zz = buffer_1100_ppdd[1115];

    auto g_y_x_0_0_y_y_xx_xx = buffer_1100_ppdd[1116];

    auto g_y_x_0_0_y_y_xx_xy = buffer_1100_ppdd[1117];

    auto g_y_x_0_0_y_y_xx_xz = buffer_1100_ppdd[1118];

    auto g_y_x_0_0_y_y_xx_yy = buffer_1100_ppdd[1119];

    auto g_y_x_0_0_y_y_xx_yz = buffer_1100_ppdd[1120];

    auto g_y_x_0_0_y_y_xx_zz = buffer_1100_ppdd[1121];

    auto g_y_x_0_0_y_y_xy_xx = buffer_1100_ppdd[1122];

    auto g_y_x_0_0_y_y_xy_xy = buffer_1100_ppdd[1123];

    auto g_y_x_0_0_y_y_xy_xz = buffer_1100_ppdd[1124];

    auto g_y_x_0_0_y_y_xy_yy = buffer_1100_ppdd[1125];

    auto g_y_x_0_0_y_y_xy_yz = buffer_1100_ppdd[1126];

    auto g_y_x_0_0_y_y_xy_zz = buffer_1100_ppdd[1127];

    auto g_y_x_0_0_y_y_xz_xx = buffer_1100_ppdd[1128];

    auto g_y_x_0_0_y_y_xz_xy = buffer_1100_ppdd[1129];

    auto g_y_x_0_0_y_y_xz_xz = buffer_1100_ppdd[1130];

    auto g_y_x_0_0_y_y_xz_yy = buffer_1100_ppdd[1131];

    auto g_y_x_0_0_y_y_xz_yz = buffer_1100_ppdd[1132];

    auto g_y_x_0_0_y_y_xz_zz = buffer_1100_ppdd[1133];

    auto g_y_x_0_0_y_y_yy_xx = buffer_1100_ppdd[1134];

    auto g_y_x_0_0_y_y_yy_xy = buffer_1100_ppdd[1135];

    auto g_y_x_0_0_y_y_yy_xz = buffer_1100_ppdd[1136];

    auto g_y_x_0_0_y_y_yy_yy = buffer_1100_ppdd[1137];

    auto g_y_x_0_0_y_y_yy_yz = buffer_1100_ppdd[1138];

    auto g_y_x_0_0_y_y_yy_zz = buffer_1100_ppdd[1139];

    auto g_y_x_0_0_y_y_yz_xx = buffer_1100_ppdd[1140];

    auto g_y_x_0_0_y_y_yz_xy = buffer_1100_ppdd[1141];

    auto g_y_x_0_0_y_y_yz_xz = buffer_1100_ppdd[1142];

    auto g_y_x_0_0_y_y_yz_yy = buffer_1100_ppdd[1143];

    auto g_y_x_0_0_y_y_yz_yz = buffer_1100_ppdd[1144];

    auto g_y_x_0_0_y_y_yz_zz = buffer_1100_ppdd[1145];

    auto g_y_x_0_0_y_y_zz_xx = buffer_1100_ppdd[1146];

    auto g_y_x_0_0_y_y_zz_xy = buffer_1100_ppdd[1147];

    auto g_y_x_0_0_y_y_zz_xz = buffer_1100_ppdd[1148];

    auto g_y_x_0_0_y_y_zz_yy = buffer_1100_ppdd[1149];

    auto g_y_x_0_0_y_y_zz_yz = buffer_1100_ppdd[1150];

    auto g_y_x_0_0_y_y_zz_zz = buffer_1100_ppdd[1151];

    auto g_y_x_0_0_y_z_xx_xx = buffer_1100_ppdd[1152];

    auto g_y_x_0_0_y_z_xx_xy = buffer_1100_ppdd[1153];

    auto g_y_x_0_0_y_z_xx_xz = buffer_1100_ppdd[1154];

    auto g_y_x_0_0_y_z_xx_yy = buffer_1100_ppdd[1155];

    auto g_y_x_0_0_y_z_xx_yz = buffer_1100_ppdd[1156];

    auto g_y_x_0_0_y_z_xx_zz = buffer_1100_ppdd[1157];

    auto g_y_x_0_0_y_z_xy_xx = buffer_1100_ppdd[1158];

    auto g_y_x_0_0_y_z_xy_xy = buffer_1100_ppdd[1159];

    auto g_y_x_0_0_y_z_xy_xz = buffer_1100_ppdd[1160];

    auto g_y_x_0_0_y_z_xy_yy = buffer_1100_ppdd[1161];

    auto g_y_x_0_0_y_z_xy_yz = buffer_1100_ppdd[1162];

    auto g_y_x_0_0_y_z_xy_zz = buffer_1100_ppdd[1163];

    auto g_y_x_0_0_y_z_xz_xx = buffer_1100_ppdd[1164];

    auto g_y_x_0_0_y_z_xz_xy = buffer_1100_ppdd[1165];

    auto g_y_x_0_0_y_z_xz_xz = buffer_1100_ppdd[1166];

    auto g_y_x_0_0_y_z_xz_yy = buffer_1100_ppdd[1167];

    auto g_y_x_0_0_y_z_xz_yz = buffer_1100_ppdd[1168];

    auto g_y_x_0_0_y_z_xz_zz = buffer_1100_ppdd[1169];

    auto g_y_x_0_0_y_z_yy_xx = buffer_1100_ppdd[1170];

    auto g_y_x_0_0_y_z_yy_xy = buffer_1100_ppdd[1171];

    auto g_y_x_0_0_y_z_yy_xz = buffer_1100_ppdd[1172];

    auto g_y_x_0_0_y_z_yy_yy = buffer_1100_ppdd[1173];

    auto g_y_x_0_0_y_z_yy_yz = buffer_1100_ppdd[1174];

    auto g_y_x_0_0_y_z_yy_zz = buffer_1100_ppdd[1175];

    auto g_y_x_0_0_y_z_yz_xx = buffer_1100_ppdd[1176];

    auto g_y_x_0_0_y_z_yz_xy = buffer_1100_ppdd[1177];

    auto g_y_x_0_0_y_z_yz_xz = buffer_1100_ppdd[1178];

    auto g_y_x_0_0_y_z_yz_yy = buffer_1100_ppdd[1179];

    auto g_y_x_0_0_y_z_yz_yz = buffer_1100_ppdd[1180];

    auto g_y_x_0_0_y_z_yz_zz = buffer_1100_ppdd[1181];

    auto g_y_x_0_0_y_z_zz_xx = buffer_1100_ppdd[1182];

    auto g_y_x_0_0_y_z_zz_xy = buffer_1100_ppdd[1183];

    auto g_y_x_0_0_y_z_zz_xz = buffer_1100_ppdd[1184];

    auto g_y_x_0_0_y_z_zz_yy = buffer_1100_ppdd[1185];

    auto g_y_x_0_0_y_z_zz_yz = buffer_1100_ppdd[1186];

    auto g_y_x_0_0_y_z_zz_zz = buffer_1100_ppdd[1187];

    auto g_y_x_0_0_z_x_xx_xx = buffer_1100_ppdd[1188];

    auto g_y_x_0_0_z_x_xx_xy = buffer_1100_ppdd[1189];

    auto g_y_x_0_0_z_x_xx_xz = buffer_1100_ppdd[1190];

    auto g_y_x_0_0_z_x_xx_yy = buffer_1100_ppdd[1191];

    auto g_y_x_0_0_z_x_xx_yz = buffer_1100_ppdd[1192];

    auto g_y_x_0_0_z_x_xx_zz = buffer_1100_ppdd[1193];

    auto g_y_x_0_0_z_x_xy_xx = buffer_1100_ppdd[1194];

    auto g_y_x_0_0_z_x_xy_xy = buffer_1100_ppdd[1195];

    auto g_y_x_0_0_z_x_xy_xz = buffer_1100_ppdd[1196];

    auto g_y_x_0_0_z_x_xy_yy = buffer_1100_ppdd[1197];

    auto g_y_x_0_0_z_x_xy_yz = buffer_1100_ppdd[1198];

    auto g_y_x_0_0_z_x_xy_zz = buffer_1100_ppdd[1199];

    auto g_y_x_0_0_z_x_xz_xx = buffer_1100_ppdd[1200];

    auto g_y_x_0_0_z_x_xz_xy = buffer_1100_ppdd[1201];

    auto g_y_x_0_0_z_x_xz_xz = buffer_1100_ppdd[1202];

    auto g_y_x_0_0_z_x_xz_yy = buffer_1100_ppdd[1203];

    auto g_y_x_0_0_z_x_xz_yz = buffer_1100_ppdd[1204];

    auto g_y_x_0_0_z_x_xz_zz = buffer_1100_ppdd[1205];

    auto g_y_x_0_0_z_x_yy_xx = buffer_1100_ppdd[1206];

    auto g_y_x_0_0_z_x_yy_xy = buffer_1100_ppdd[1207];

    auto g_y_x_0_0_z_x_yy_xz = buffer_1100_ppdd[1208];

    auto g_y_x_0_0_z_x_yy_yy = buffer_1100_ppdd[1209];

    auto g_y_x_0_0_z_x_yy_yz = buffer_1100_ppdd[1210];

    auto g_y_x_0_0_z_x_yy_zz = buffer_1100_ppdd[1211];

    auto g_y_x_0_0_z_x_yz_xx = buffer_1100_ppdd[1212];

    auto g_y_x_0_0_z_x_yz_xy = buffer_1100_ppdd[1213];

    auto g_y_x_0_0_z_x_yz_xz = buffer_1100_ppdd[1214];

    auto g_y_x_0_0_z_x_yz_yy = buffer_1100_ppdd[1215];

    auto g_y_x_0_0_z_x_yz_yz = buffer_1100_ppdd[1216];

    auto g_y_x_0_0_z_x_yz_zz = buffer_1100_ppdd[1217];

    auto g_y_x_0_0_z_x_zz_xx = buffer_1100_ppdd[1218];

    auto g_y_x_0_0_z_x_zz_xy = buffer_1100_ppdd[1219];

    auto g_y_x_0_0_z_x_zz_xz = buffer_1100_ppdd[1220];

    auto g_y_x_0_0_z_x_zz_yy = buffer_1100_ppdd[1221];

    auto g_y_x_0_0_z_x_zz_yz = buffer_1100_ppdd[1222];

    auto g_y_x_0_0_z_x_zz_zz = buffer_1100_ppdd[1223];

    auto g_y_x_0_0_z_y_xx_xx = buffer_1100_ppdd[1224];

    auto g_y_x_0_0_z_y_xx_xy = buffer_1100_ppdd[1225];

    auto g_y_x_0_0_z_y_xx_xz = buffer_1100_ppdd[1226];

    auto g_y_x_0_0_z_y_xx_yy = buffer_1100_ppdd[1227];

    auto g_y_x_0_0_z_y_xx_yz = buffer_1100_ppdd[1228];

    auto g_y_x_0_0_z_y_xx_zz = buffer_1100_ppdd[1229];

    auto g_y_x_0_0_z_y_xy_xx = buffer_1100_ppdd[1230];

    auto g_y_x_0_0_z_y_xy_xy = buffer_1100_ppdd[1231];

    auto g_y_x_0_0_z_y_xy_xz = buffer_1100_ppdd[1232];

    auto g_y_x_0_0_z_y_xy_yy = buffer_1100_ppdd[1233];

    auto g_y_x_0_0_z_y_xy_yz = buffer_1100_ppdd[1234];

    auto g_y_x_0_0_z_y_xy_zz = buffer_1100_ppdd[1235];

    auto g_y_x_0_0_z_y_xz_xx = buffer_1100_ppdd[1236];

    auto g_y_x_0_0_z_y_xz_xy = buffer_1100_ppdd[1237];

    auto g_y_x_0_0_z_y_xz_xz = buffer_1100_ppdd[1238];

    auto g_y_x_0_0_z_y_xz_yy = buffer_1100_ppdd[1239];

    auto g_y_x_0_0_z_y_xz_yz = buffer_1100_ppdd[1240];

    auto g_y_x_0_0_z_y_xz_zz = buffer_1100_ppdd[1241];

    auto g_y_x_0_0_z_y_yy_xx = buffer_1100_ppdd[1242];

    auto g_y_x_0_0_z_y_yy_xy = buffer_1100_ppdd[1243];

    auto g_y_x_0_0_z_y_yy_xz = buffer_1100_ppdd[1244];

    auto g_y_x_0_0_z_y_yy_yy = buffer_1100_ppdd[1245];

    auto g_y_x_0_0_z_y_yy_yz = buffer_1100_ppdd[1246];

    auto g_y_x_0_0_z_y_yy_zz = buffer_1100_ppdd[1247];

    auto g_y_x_0_0_z_y_yz_xx = buffer_1100_ppdd[1248];

    auto g_y_x_0_0_z_y_yz_xy = buffer_1100_ppdd[1249];

    auto g_y_x_0_0_z_y_yz_xz = buffer_1100_ppdd[1250];

    auto g_y_x_0_0_z_y_yz_yy = buffer_1100_ppdd[1251];

    auto g_y_x_0_0_z_y_yz_yz = buffer_1100_ppdd[1252];

    auto g_y_x_0_0_z_y_yz_zz = buffer_1100_ppdd[1253];

    auto g_y_x_0_0_z_y_zz_xx = buffer_1100_ppdd[1254];

    auto g_y_x_0_0_z_y_zz_xy = buffer_1100_ppdd[1255];

    auto g_y_x_0_0_z_y_zz_xz = buffer_1100_ppdd[1256];

    auto g_y_x_0_0_z_y_zz_yy = buffer_1100_ppdd[1257];

    auto g_y_x_0_0_z_y_zz_yz = buffer_1100_ppdd[1258];

    auto g_y_x_0_0_z_y_zz_zz = buffer_1100_ppdd[1259];

    auto g_y_x_0_0_z_z_xx_xx = buffer_1100_ppdd[1260];

    auto g_y_x_0_0_z_z_xx_xy = buffer_1100_ppdd[1261];

    auto g_y_x_0_0_z_z_xx_xz = buffer_1100_ppdd[1262];

    auto g_y_x_0_0_z_z_xx_yy = buffer_1100_ppdd[1263];

    auto g_y_x_0_0_z_z_xx_yz = buffer_1100_ppdd[1264];

    auto g_y_x_0_0_z_z_xx_zz = buffer_1100_ppdd[1265];

    auto g_y_x_0_0_z_z_xy_xx = buffer_1100_ppdd[1266];

    auto g_y_x_0_0_z_z_xy_xy = buffer_1100_ppdd[1267];

    auto g_y_x_0_0_z_z_xy_xz = buffer_1100_ppdd[1268];

    auto g_y_x_0_0_z_z_xy_yy = buffer_1100_ppdd[1269];

    auto g_y_x_0_0_z_z_xy_yz = buffer_1100_ppdd[1270];

    auto g_y_x_0_0_z_z_xy_zz = buffer_1100_ppdd[1271];

    auto g_y_x_0_0_z_z_xz_xx = buffer_1100_ppdd[1272];

    auto g_y_x_0_0_z_z_xz_xy = buffer_1100_ppdd[1273];

    auto g_y_x_0_0_z_z_xz_xz = buffer_1100_ppdd[1274];

    auto g_y_x_0_0_z_z_xz_yy = buffer_1100_ppdd[1275];

    auto g_y_x_0_0_z_z_xz_yz = buffer_1100_ppdd[1276];

    auto g_y_x_0_0_z_z_xz_zz = buffer_1100_ppdd[1277];

    auto g_y_x_0_0_z_z_yy_xx = buffer_1100_ppdd[1278];

    auto g_y_x_0_0_z_z_yy_xy = buffer_1100_ppdd[1279];

    auto g_y_x_0_0_z_z_yy_xz = buffer_1100_ppdd[1280];

    auto g_y_x_0_0_z_z_yy_yy = buffer_1100_ppdd[1281];

    auto g_y_x_0_0_z_z_yy_yz = buffer_1100_ppdd[1282];

    auto g_y_x_0_0_z_z_yy_zz = buffer_1100_ppdd[1283];

    auto g_y_x_0_0_z_z_yz_xx = buffer_1100_ppdd[1284];

    auto g_y_x_0_0_z_z_yz_xy = buffer_1100_ppdd[1285];

    auto g_y_x_0_0_z_z_yz_xz = buffer_1100_ppdd[1286];

    auto g_y_x_0_0_z_z_yz_yy = buffer_1100_ppdd[1287];

    auto g_y_x_0_0_z_z_yz_yz = buffer_1100_ppdd[1288];

    auto g_y_x_0_0_z_z_yz_zz = buffer_1100_ppdd[1289];

    auto g_y_x_0_0_z_z_zz_xx = buffer_1100_ppdd[1290];

    auto g_y_x_0_0_z_z_zz_xy = buffer_1100_ppdd[1291];

    auto g_y_x_0_0_z_z_zz_xz = buffer_1100_ppdd[1292];

    auto g_y_x_0_0_z_z_zz_yy = buffer_1100_ppdd[1293];

    auto g_y_x_0_0_z_z_zz_yz = buffer_1100_ppdd[1294];

    auto g_y_x_0_0_z_z_zz_zz = buffer_1100_ppdd[1295];

    auto g_y_y_0_0_x_x_xx_xx = buffer_1100_ppdd[1296];

    auto g_y_y_0_0_x_x_xx_xy = buffer_1100_ppdd[1297];

    auto g_y_y_0_0_x_x_xx_xz = buffer_1100_ppdd[1298];

    auto g_y_y_0_0_x_x_xx_yy = buffer_1100_ppdd[1299];

    auto g_y_y_0_0_x_x_xx_yz = buffer_1100_ppdd[1300];

    auto g_y_y_0_0_x_x_xx_zz = buffer_1100_ppdd[1301];

    auto g_y_y_0_0_x_x_xy_xx = buffer_1100_ppdd[1302];

    auto g_y_y_0_0_x_x_xy_xy = buffer_1100_ppdd[1303];

    auto g_y_y_0_0_x_x_xy_xz = buffer_1100_ppdd[1304];

    auto g_y_y_0_0_x_x_xy_yy = buffer_1100_ppdd[1305];

    auto g_y_y_0_0_x_x_xy_yz = buffer_1100_ppdd[1306];

    auto g_y_y_0_0_x_x_xy_zz = buffer_1100_ppdd[1307];

    auto g_y_y_0_0_x_x_xz_xx = buffer_1100_ppdd[1308];

    auto g_y_y_0_0_x_x_xz_xy = buffer_1100_ppdd[1309];

    auto g_y_y_0_0_x_x_xz_xz = buffer_1100_ppdd[1310];

    auto g_y_y_0_0_x_x_xz_yy = buffer_1100_ppdd[1311];

    auto g_y_y_0_0_x_x_xz_yz = buffer_1100_ppdd[1312];

    auto g_y_y_0_0_x_x_xz_zz = buffer_1100_ppdd[1313];

    auto g_y_y_0_0_x_x_yy_xx = buffer_1100_ppdd[1314];

    auto g_y_y_0_0_x_x_yy_xy = buffer_1100_ppdd[1315];

    auto g_y_y_0_0_x_x_yy_xz = buffer_1100_ppdd[1316];

    auto g_y_y_0_0_x_x_yy_yy = buffer_1100_ppdd[1317];

    auto g_y_y_0_0_x_x_yy_yz = buffer_1100_ppdd[1318];

    auto g_y_y_0_0_x_x_yy_zz = buffer_1100_ppdd[1319];

    auto g_y_y_0_0_x_x_yz_xx = buffer_1100_ppdd[1320];

    auto g_y_y_0_0_x_x_yz_xy = buffer_1100_ppdd[1321];

    auto g_y_y_0_0_x_x_yz_xz = buffer_1100_ppdd[1322];

    auto g_y_y_0_0_x_x_yz_yy = buffer_1100_ppdd[1323];

    auto g_y_y_0_0_x_x_yz_yz = buffer_1100_ppdd[1324];

    auto g_y_y_0_0_x_x_yz_zz = buffer_1100_ppdd[1325];

    auto g_y_y_0_0_x_x_zz_xx = buffer_1100_ppdd[1326];

    auto g_y_y_0_0_x_x_zz_xy = buffer_1100_ppdd[1327];

    auto g_y_y_0_0_x_x_zz_xz = buffer_1100_ppdd[1328];

    auto g_y_y_0_0_x_x_zz_yy = buffer_1100_ppdd[1329];

    auto g_y_y_0_0_x_x_zz_yz = buffer_1100_ppdd[1330];

    auto g_y_y_0_0_x_x_zz_zz = buffer_1100_ppdd[1331];

    auto g_y_y_0_0_x_y_xx_xx = buffer_1100_ppdd[1332];

    auto g_y_y_0_0_x_y_xx_xy = buffer_1100_ppdd[1333];

    auto g_y_y_0_0_x_y_xx_xz = buffer_1100_ppdd[1334];

    auto g_y_y_0_0_x_y_xx_yy = buffer_1100_ppdd[1335];

    auto g_y_y_0_0_x_y_xx_yz = buffer_1100_ppdd[1336];

    auto g_y_y_0_0_x_y_xx_zz = buffer_1100_ppdd[1337];

    auto g_y_y_0_0_x_y_xy_xx = buffer_1100_ppdd[1338];

    auto g_y_y_0_0_x_y_xy_xy = buffer_1100_ppdd[1339];

    auto g_y_y_0_0_x_y_xy_xz = buffer_1100_ppdd[1340];

    auto g_y_y_0_0_x_y_xy_yy = buffer_1100_ppdd[1341];

    auto g_y_y_0_0_x_y_xy_yz = buffer_1100_ppdd[1342];

    auto g_y_y_0_0_x_y_xy_zz = buffer_1100_ppdd[1343];

    auto g_y_y_0_0_x_y_xz_xx = buffer_1100_ppdd[1344];

    auto g_y_y_0_0_x_y_xz_xy = buffer_1100_ppdd[1345];

    auto g_y_y_0_0_x_y_xz_xz = buffer_1100_ppdd[1346];

    auto g_y_y_0_0_x_y_xz_yy = buffer_1100_ppdd[1347];

    auto g_y_y_0_0_x_y_xz_yz = buffer_1100_ppdd[1348];

    auto g_y_y_0_0_x_y_xz_zz = buffer_1100_ppdd[1349];

    auto g_y_y_0_0_x_y_yy_xx = buffer_1100_ppdd[1350];

    auto g_y_y_0_0_x_y_yy_xy = buffer_1100_ppdd[1351];

    auto g_y_y_0_0_x_y_yy_xz = buffer_1100_ppdd[1352];

    auto g_y_y_0_0_x_y_yy_yy = buffer_1100_ppdd[1353];

    auto g_y_y_0_0_x_y_yy_yz = buffer_1100_ppdd[1354];

    auto g_y_y_0_0_x_y_yy_zz = buffer_1100_ppdd[1355];

    auto g_y_y_0_0_x_y_yz_xx = buffer_1100_ppdd[1356];

    auto g_y_y_0_0_x_y_yz_xy = buffer_1100_ppdd[1357];

    auto g_y_y_0_0_x_y_yz_xz = buffer_1100_ppdd[1358];

    auto g_y_y_0_0_x_y_yz_yy = buffer_1100_ppdd[1359];

    auto g_y_y_0_0_x_y_yz_yz = buffer_1100_ppdd[1360];

    auto g_y_y_0_0_x_y_yz_zz = buffer_1100_ppdd[1361];

    auto g_y_y_0_0_x_y_zz_xx = buffer_1100_ppdd[1362];

    auto g_y_y_0_0_x_y_zz_xy = buffer_1100_ppdd[1363];

    auto g_y_y_0_0_x_y_zz_xz = buffer_1100_ppdd[1364];

    auto g_y_y_0_0_x_y_zz_yy = buffer_1100_ppdd[1365];

    auto g_y_y_0_0_x_y_zz_yz = buffer_1100_ppdd[1366];

    auto g_y_y_0_0_x_y_zz_zz = buffer_1100_ppdd[1367];

    auto g_y_y_0_0_x_z_xx_xx = buffer_1100_ppdd[1368];

    auto g_y_y_0_0_x_z_xx_xy = buffer_1100_ppdd[1369];

    auto g_y_y_0_0_x_z_xx_xz = buffer_1100_ppdd[1370];

    auto g_y_y_0_0_x_z_xx_yy = buffer_1100_ppdd[1371];

    auto g_y_y_0_0_x_z_xx_yz = buffer_1100_ppdd[1372];

    auto g_y_y_0_0_x_z_xx_zz = buffer_1100_ppdd[1373];

    auto g_y_y_0_0_x_z_xy_xx = buffer_1100_ppdd[1374];

    auto g_y_y_0_0_x_z_xy_xy = buffer_1100_ppdd[1375];

    auto g_y_y_0_0_x_z_xy_xz = buffer_1100_ppdd[1376];

    auto g_y_y_0_0_x_z_xy_yy = buffer_1100_ppdd[1377];

    auto g_y_y_0_0_x_z_xy_yz = buffer_1100_ppdd[1378];

    auto g_y_y_0_0_x_z_xy_zz = buffer_1100_ppdd[1379];

    auto g_y_y_0_0_x_z_xz_xx = buffer_1100_ppdd[1380];

    auto g_y_y_0_0_x_z_xz_xy = buffer_1100_ppdd[1381];

    auto g_y_y_0_0_x_z_xz_xz = buffer_1100_ppdd[1382];

    auto g_y_y_0_0_x_z_xz_yy = buffer_1100_ppdd[1383];

    auto g_y_y_0_0_x_z_xz_yz = buffer_1100_ppdd[1384];

    auto g_y_y_0_0_x_z_xz_zz = buffer_1100_ppdd[1385];

    auto g_y_y_0_0_x_z_yy_xx = buffer_1100_ppdd[1386];

    auto g_y_y_0_0_x_z_yy_xy = buffer_1100_ppdd[1387];

    auto g_y_y_0_0_x_z_yy_xz = buffer_1100_ppdd[1388];

    auto g_y_y_0_0_x_z_yy_yy = buffer_1100_ppdd[1389];

    auto g_y_y_0_0_x_z_yy_yz = buffer_1100_ppdd[1390];

    auto g_y_y_0_0_x_z_yy_zz = buffer_1100_ppdd[1391];

    auto g_y_y_0_0_x_z_yz_xx = buffer_1100_ppdd[1392];

    auto g_y_y_0_0_x_z_yz_xy = buffer_1100_ppdd[1393];

    auto g_y_y_0_0_x_z_yz_xz = buffer_1100_ppdd[1394];

    auto g_y_y_0_0_x_z_yz_yy = buffer_1100_ppdd[1395];

    auto g_y_y_0_0_x_z_yz_yz = buffer_1100_ppdd[1396];

    auto g_y_y_0_0_x_z_yz_zz = buffer_1100_ppdd[1397];

    auto g_y_y_0_0_x_z_zz_xx = buffer_1100_ppdd[1398];

    auto g_y_y_0_0_x_z_zz_xy = buffer_1100_ppdd[1399];

    auto g_y_y_0_0_x_z_zz_xz = buffer_1100_ppdd[1400];

    auto g_y_y_0_0_x_z_zz_yy = buffer_1100_ppdd[1401];

    auto g_y_y_0_0_x_z_zz_yz = buffer_1100_ppdd[1402];

    auto g_y_y_0_0_x_z_zz_zz = buffer_1100_ppdd[1403];

    auto g_y_y_0_0_y_x_xx_xx = buffer_1100_ppdd[1404];

    auto g_y_y_0_0_y_x_xx_xy = buffer_1100_ppdd[1405];

    auto g_y_y_0_0_y_x_xx_xz = buffer_1100_ppdd[1406];

    auto g_y_y_0_0_y_x_xx_yy = buffer_1100_ppdd[1407];

    auto g_y_y_0_0_y_x_xx_yz = buffer_1100_ppdd[1408];

    auto g_y_y_0_0_y_x_xx_zz = buffer_1100_ppdd[1409];

    auto g_y_y_0_0_y_x_xy_xx = buffer_1100_ppdd[1410];

    auto g_y_y_0_0_y_x_xy_xy = buffer_1100_ppdd[1411];

    auto g_y_y_0_0_y_x_xy_xz = buffer_1100_ppdd[1412];

    auto g_y_y_0_0_y_x_xy_yy = buffer_1100_ppdd[1413];

    auto g_y_y_0_0_y_x_xy_yz = buffer_1100_ppdd[1414];

    auto g_y_y_0_0_y_x_xy_zz = buffer_1100_ppdd[1415];

    auto g_y_y_0_0_y_x_xz_xx = buffer_1100_ppdd[1416];

    auto g_y_y_0_0_y_x_xz_xy = buffer_1100_ppdd[1417];

    auto g_y_y_0_0_y_x_xz_xz = buffer_1100_ppdd[1418];

    auto g_y_y_0_0_y_x_xz_yy = buffer_1100_ppdd[1419];

    auto g_y_y_0_0_y_x_xz_yz = buffer_1100_ppdd[1420];

    auto g_y_y_0_0_y_x_xz_zz = buffer_1100_ppdd[1421];

    auto g_y_y_0_0_y_x_yy_xx = buffer_1100_ppdd[1422];

    auto g_y_y_0_0_y_x_yy_xy = buffer_1100_ppdd[1423];

    auto g_y_y_0_0_y_x_yy_xz = buffer_1100_ppdd[1424];

    auto g_y_y_0_0_y_x_yy_yy = buffer_1100_ppdd[1425];

    auto g_y_y_0_0_y_x_yy_yz = buffer_1100_ppdd[1426];

    auto g_y_y_0_0_y_x_yy_zz = buffer_1100_ppdd[1427];

    auto g_y_y_0_0_y_x_yz_xx = buffer_1100_ppdd[1428];

    auto g_y_y_0_0_y_x_yz_xy = buffer_1100_ppdd[1429];

    auto g_y_y_0_0_y_x_yz_xz = buffer_1100_ppdd[1430];

    auto g_y_y_0_0_y_x_yz_yy = buffer_1100_ppdd[1431];

    auto g_y_y_0_0_y_x_yz_yz = buffer_1100_ppdd[1432];

    auto g_y_y_0_0_y_x_yz_zz = buffer_1100_ppdd[1433];

    auto g_y_y_0_0_y_x_zz_xx = buffer_1100_ppdd[1434];

    auto g_y_y_0_0_y_x_zz_xy = buffer_1100_ppdd[1435];

    auto g_y_y_0_0_y_x_zz_xz = buffer_1100_ppdd[1436];

    auto g_y_y_0_0_y_x_zz_yy = buffer_1100_ppdd[1437];

    auto g_y_y_0_0_y_x_zz_yz = buffer_1100_ppdd[1438];

    auto g_y_y_0_0_y_x_zz_zz = buffer_1100_ppdd[1439];

    auto g_y_y_0_0_y_y_xx_xx = buffer_1100_ppdd[1440];

    auto g_y_y_0_0_y_y_xx_xy = buffer_1100_ppdd[1441];

    auto g_y_y_0_0_y_y_xx_xz = buffer_1100_ppdd[1442];

    auto g_y_y_0_0_y_y_xx_yy = buffer_1100_ppdd[1443];

    auto g_y_y_0_0_y_y_xx_yz = buffer_1100_ppdd[1444];

    auto g_y_y_0_0_y_y_xx_zz = buffer_1100_ppdd[1445];

    auto g_y_y_0_0_y_y_xy_xx = buffer_1100_ppdd[1446];

    auto g_y_y_0_0_y_y_xy_xy = buffer_1100_ppdd[1447];

    auto g_y_y_0_0_y_y_xy_xz = buffer_1100_ppdd[1448];

    auto g_y_y_0_0_y_y_xy_yy = buffer_1100_ppdd[1449];

    auto g_y_y_0_0_y_y_xy_yz = buffer_1100_ppdd[1450];

    auto g_y_y_0_0_y_y_xy_zz = buffer_1100_ppdd[1451];

    auto g_y_y_0_0_y_y_xz_xx = buffer_1100_ppdd[1452];

    auto g_y_y_0_0_y_y_xz_xy = buffer_1100_ppdd[1453];

    auto g_y_y_0_0_y_y_xz_xz = buffer_1100_ppdd[1454];

    auto g_y_y_0_0_y_y_xz_yy = buffer_1100_ppdd[1455];

    auto g_y_y_0_0_y_y_xz_yz = buffer_1100_ppdd[1456];

    auto g_y_y_0_0_y_y_xz_zz = buffer_1100_ppdd[1457];

    auto g_y_y_0_0_y_y_yy_xx = buffer_1100_ppdd[1458];

    auto g_y_y_0_0_y_y_yy_xy = buffer_1100_ppdd[1459];

    auto g_y_y_0_0_y_y_yy_xz = buffer_1100_ppdd[1460];

    auto g_y_y_0_0_y_y_yy_yy = buffer_1100_ppdd[1461];

    auto g_y_y_0_0_y_y_yy_yz = buffer_1100_ppdd[1462];

    auto g_y_y_0_0_y_y_yy_zz = buffer_1100_ppdd[1463];

    auto g_y_y_0_0_y_y_yz_xx = buffer_1100_ppdd[1464];

    auto g_y_y_0_0_y_y_yz_xy = buffer_1100_ppdd[1465];

    auto g_y_y_0_0_y_y_yz_xz = buffer_1100_ppdd[1466];

    auto g_y_y_0_0_y_y_yz_yy = buffer_1100_ppdd[1467];

    auto g_y_y_0_0_y_y_yz_yz = buffer_1100_ppdd[1468];

    auto g_y_y_0_0_y_y_yz_zz = buffer_1100_ppdd[1469];

    auto g_y_y_0_0_y_y_zz_xx = buffer_1100_ppdd[1470];

    auto g_y_y_0_0_y_y_zz_xy = buffer_1100_ppdd[1471];

    auto g_y_y_0_0_y_y_zz_xz = buffer_1100_ppdd[1472];

    auto g_y_y_0_0_y_y_zz_yy = buffer_1100_ppdd[1473];

    auto g_y_y_0_0_y_y_zz_yz = buffer_1100_ppdd[1474];

    auto g_y_y_0_0_y_y_zz_zz = buffer_1100_ppdd[1475];

    auto g_y_y_0_0_y_z_xx_xx = buffer_1100_ppdd[1476];

    auto g_y_y_0_0_y_z_xx_xy = buffer_1100_ppdd[1477];

    auto g_y_y_0_0_y_z_xx_xz = buffer_1100_ppdd[1478];

    auto g_y_y_0_0_y_z_xx_yy = buffer_1100_ppdd[1479];

    auto g_y_y_0_0_y_z_xx_yz = buffer_1100_ppdd[1480];

    auto g_y_y_0_0_y_z_xx_zz = buffer_1100_ppdd[1481];

    auto g_y_y_0_0_y_z_xy_xx = buffer_1100_ppdd[1482];

    auto g_y_y_0_0_y_z_xy_xy = buffer_1100_ppdd[1483];

    auto g_y_y_0_0_y_z_xy_xz = buffer_1100_ppdd[1484];

    auto g_y_y_0_0_y_z_xy_yy = buffer_1100_ppdd[1485];

    auto g_y_y_0_0_y_z_xy_yz = buffer_1100_ppdd[1486];

    auto g_y_y_0_0_y_z_xy_zz = buffer_1100_ppdd[1487];

    auto g_y_y_0_0_y_z_xz_xx = buffer_1100_ppdd[1488];

    auto g_y_y_0_0_y_z_xz_xy = buffer_1100_ppdd[1489];

    auto g_y_y_0_0_y_z_xz_xz = buffer_1100_ppdd[1490];

    auto g_y_y_0_0_y_z_xz_yy = buffer_1100_ppdd[1491];

    auto g_y_y_0_0_y_z_xz_yz = buffer_1100_ppdd[1492];

    auto g_y_y_0_0_y_z_xz_zz = buffer_1100_ppdd[1493];

    auto g_y_y_0_0_y_z_yy_xx = buffer_1100_ppdd[1494];

    auto g_y_y_0_0_y_z_yy_xy = buffer_1100_ppdd[1495];

    auto g_y_y_0_0_y_z_yy_xz = buffer_1100_ppdd[1496];

    auto g_y_y_0_0_y_z_yy_yy = buffer_1100_ppdd[1497];

    auto g_y_y_0_0_y_z_yy_yz = buffer_1100_ppdd[1498];

    auto g_y_y_0_0_y_z_yy_zz = buffer_1100_ppdd[1499];

    auto g_y_y_0_0_y_z_yz_xx = buffer_1100_ppdd[1500];

    auto g_y_y_0_0_y_z_yz_xy = buffer_1100_ppdd[1501];

    auto g_y_y_0_0_y_z_yz_xz = buffer_1100_ppdd[1502];

    auto g_y_y_0_0_y_z_yz_yy = buffer_1100_ppdd[1503];

    auto g_y_y_0_0_y_z_yz_yz = buffer_1100_ppdd[1504];

    auto g_y_y_0_0_y_z_yz_zz = buffer_1100_ppdd[1505];

    auto g_y_y_0_0_y_z_zz_xx = buffer_1100_ppdd[1506];

    auto g_y_y_0_0_y_z_zz_xy = buffer_1100_ppdd[1507];

    auto g_y_y_0_0_y_z_zz_xz = buffer_1100_ppdd[1508];

    auto g_y_y_0_0_y_z_zz_yy = buffer_1100_ppdd[1509];

    auto g_y_y_0_0_y_z_zz_yz = buffer_1100_ppdd[1510];

    auto g_y_y_0_0_y_z_zz_zz = buffer_1100_ppdd[1511];

    auto g_y_y_0_0_z_x_xx_xx = buffer_1100_ppdd[1512];

    auto g_y_y_0_0_z_x_xx_xy = buffer_1100_ppdd[1513];

    auto g_y_y_0_0_z_x_xx_xz = buffer_1100_ppdd[1514];

    auto g_y_y_0_0_z_x_xx_yy = buffer_1100_ppdd[1515];

    auto g_y_y_0_0_z_x_xx_yz = buffer_1100_ppdd[1516];

    auto g_y_y_0_0_z_x_xx_zz = buffer_1100_ppdd[1517];

    auto g_y_y_0_0_z_x_xy_xx = buffer_1100_ppdd[1518];

    auto g_y_y_0_0_z_x_xy_xy = buffer_1100_ppdd[1519];

    auto g_y_y_0_0_z_x_xy_xz = buffer_1100_ppdd[1520];

    auto g_y_y_0_0_z_x_xy_yy = buffer_1100_ppdd[1521];

    auto g_y_y_0_0_z_x_xy_yz = buffer_1100_ppdd[1522];

    auto g_y_y_0_0_z_x_xy_zz = buffer_1100_ppdd[1523];

    auto g_y_y_0_0_z_x_xz_xx = buffer_1100_ppdd[1524];

    auto g_y_y_0_0_z_x_xz_xy = buffer_1100_ppdd[1525];

    auto g_y_y_0_0_z_x_xz_xz = buffer_1100_ppdd[1526];

    auto g_y_y_0_0_z_x_xz_yy = buffer_1100_ppdd[1527];

    auto g_y_y_0_0_z_x_xz_yz = buffer_1100_ppdd[1528];

    auto g_y_y_0_0_z_x_xz_zz = buffer_1100_ppdd[1529];

    auto g_y_y_0_0_z_x_yy_xx = buffer_1100_ppdd[1530];

    auto g_y_y_0_0_z_x_yy_xy = buffer_1100_ppdd[1531];

    auto g_y_y_0_0_z_x_yy_xz = buffer_1100_ppdd[1532];

    auto g_y_y_0_0_z_x_yy_yy = buffer_1100_ppdd[1533];

    auto g_y_y_0_0_z_x_yy_yz = buffer_1100_ppdd[1534];

    auto g_y_y_0_0_z_x_yy_zz = buffer_1100_ppdd[1535];

    auto g_y_y_0_0_z_x_yz_xx = buffer_1100_ppdd[1536];

    auto g_y_y_0_0_z_x_yz_xy = buffer_1100_ppdd[1537];

    auto g_y_y_0_0_z_x_yz_xz = buffer_1100_ppdd[1538];

    auto g_y_y_0_0_z_x_yz_yy = buffer_1100_ppdd[1539];

    auto g_y_y_0_0_z_x_yz_yz = buffer_1100_ppdd[1540];

    auto g_y_y_0_0_z_x_yz_zz = buffer_1100_ppdd[1541];

    auto g_y_y_0_0_z_x_zz_xx = buffer_1100_ppdd[1542];

    auto g_y_y_0_0_z_x_zz_xy = buffer_1100_ppdd[1543];

    auto g_y_y_0_0_z_x_zz_xz = buffer_1100_ppdd[1544];

    auto g_y_y_0_0_z_x_zz_yy = buffer_1100_ppdd[1545];

    auto g_y_y_0_0_z_x_zz_yz = buffer_1100_ppdd[1546];

    auto g_y_y_0_0_z_x_zz_zz = buffer_1100_ppdd[1547];

    auto g_y_y_0_0_z_y_xx_xx = buffer_1100_ppdd[1548];

    auto g_y_y_0_0_z_y_xx_xy = buffer_1100_ppdd[1549];

    auto g_y_y_0_0_z_y_xx_xz = buffer_1100_ppdd[1550];

    auto g_y_y_0_0_z_y_xx_yy = buffer_1100_ppdd[1551];

    auto g_y_y_0_0_z_y_xx_yz = buffer_1100_ppdd[1552];

    auto g_y_y_0_0_z_y_xx_zz = buffer_1100_ppdd[1553];

    auto g_y_y_0_0_z_y_xy_xx = buffer_1100_ppdd[1554];

    auto g_y_y_0_0_z_y_xy_xy = buffer_1100_ppdd[1555];

    auto g_y_y_0_0_z_y_xy_xz = buffer_1100_ppdd[1556];

    auto g_y_y_0_0_z_y_xy_yy = buffer_1100_ppdd[1557];

    auto g_y_y_0_0_z_y_xy_yz = buffer_1100_ppdd[1558];

    auto g_y_y_0_0_z_y_xy_zz = buffer_1100_ppdd[1559];

    auto g_y_y_0_0_z_y_xz_xx = buffer_1100_ppdd[1560];

    auto g_y_y_0_0_z_y_xz_xy = buffer_1100_ppdd[1561];

    auto g_y_y_0_0_z_y_xz_xz = buffer_1100_ppdd[1562];

    auto g_y_y_0_0_z_y_xz_yy = buffer_1100_ppdd[1563];

    auto g_y_y_0_0_z_y_xz_yz = buffer_1100_ppdd[1564];

    auto g_y_y_0_0_z_y_xz_zz = buffer_1100_ppdd[1565];

    auto g_y_y_0_0_z_y_yy_xx = buffer_1100_ppdd[1566];

    auto g_y_y_0_0_z_y_yy_xy = buffer_1100_ppdd[1567];

    auto g_y_y_0_0_z_y_yy_xz = buffer_1100_ppdd[1568];

    auto g_y_y_0_0_z_y_yy_yy = buffer_1100_ppdd[1569];

    auto g_y_y_0_0_z_y_yy_yz = buffer_1100_ppdd[1570];

    auto g_y_y_0_0_z_y_yy_zz = buffer_1100_ppdd[1571];

    auto g_y_y_0_0_z_y_yz_xx = buffer_1100_ppdd[1572];

    auto g_y_y_0_0_z_y_yz_xy = buffer_1100_ppdd[1573];

    auto g_y_y_0_0_z_y_yz_xz = buffer_1100_ppdd[1574];

    auto g_y_y_0_0_z_y_yz_yy = buffer_1100_ppdd[1575];

    auto g_y_y_0_0_z_y_yz_yz = buffer_1100_ppdd[1576];

    auto g_y_y_0_0_z_y_yz_zz = buffer_1100_ppdd[1577];

    auto g_y_y_0_0_z_y_zz_xx = buffer_1100_ppdd[1578];

    auto g_y_y_0_0_z_y_zz_xy = buffer_1100_ppdd[1579];

    auto g_y_y_0_0_z_y_zz_xz = buffer_1100_ppdd[1580];

    auto g_y_y_0_0_z_y_zz_yy = buffer_1100_ppdd[1581];

    auto g_y_y_0_0_z_y_zz_yz = buffer_1100_ppdd[1582];

    auto g_y_y_0_0_z_y_zz_zz = buffer_1100_ppdd[1583];

    auto g_y_y_0_0_z_z_xx_xx = buffer_1100_ppdd[1584];

    auto g_y_y_0_0_z_z_xx_xy = buffer_1100_ppdd[1585];

    auto g_y_y_0_0_z_z_xx_xz = buffer_1100_ppdd[1586];

    auto g_y_y_0_0_z_z_xx_yy = buffer_1100_ppdd[1587];

    auto g_y_y_0_0_z_z_xx_yz = buffer_1100_ppdd[1588];

    auto g_y_y_0_0_z_z_xx_zz = buffer_1100_ppdd[1589];

    auto g_y_y_0_0_z_z_xy_xx = buffer_1100_ppdd[1590];

    auto g_y_y_0_0_z_z_xy_xy = buffer_1100_ppdd[1591];

    auto g_y_y_0_0_z_z_xy_xz = buffer_1100_ppdd[1592];

    auto g_y_y_0_0_z_z_xy_yy = buffer_1100_ppdd[1593];

    auto g_y_y_0_0_z_z_xy_yz = buffer_1100_ppdd[1594];

    auto g_y_y_0_0_z_z_xy_zz = buffer_1100_ppdd[1595];

    auto g_y_y_0_0_z_z_xz_xx = buffer_1100_ppdd[1596];

    auto g_y_y_0_0_z_z_xz_xy = buffer_1100_ppdd[1597];

    auto g_y_y_0_0_z_z_xz_xz = buffer_1100_ppdd[1598];

    auto g_y_y_0_0_z_z_xz_yy = buffer_1100_ppdd[1599];

    auto g_y_y_0_0_z_z_xz_yz = buffer_1100_ppdd[1600];

    auto g_y_y_0_0_z_z_xz_zz = buffer_1100_ppdd[1601];

    auto g_y_y_0_0_z_z_yy_xx = buffer_1100_ppdd[1602];

    auto g_y_y_0_0_z_z_yy_xy = buffer_1100_ppdd[1603];

    auto g_y_y_0_0_z_z_yy_xz = buffer_1100_ppdd[1604];

    auto g_y_y_0_0_z_z_yy_yy = buffer_1100_ppdd[1605];

    auto g_y_y_0_0_z_z_yy_yz = buffer_1100_ppdd[1606];

    auto g_y_y_0_0_z_z_yy_zz = buffer_1100_ppdd[1607];

    auto g_y_y_0_0_z_z_yz_xx = buffer_1100_ppdd[1608];

    auto g_y_y_0_0_z_z_yz_xy = buffer_1100_ppdd[1609];

    auto g_y_y_0_0_z_z_yz_xz = buffer_1100_ppdd[1610];

    auto g_y_y_0_0_z_z_yz_yy = buffer_1100_ppdd[1611];

    auto g_y_y_0_0_z_z_yz_yz = buffer_1100_ppdd[1612];

    auto g_y_y_0_0_z_z_yz_zz = buffer_1100_ppdd[1613];

    auto g_y_y_0_0_z_z_zz_xx = buffer_1100_ppdd[1614];

    auto g_y_y_0_0_z_z_zz_xy = buffer_1100_ppdd[1615];

    auto g_y_y_0_0_z_z_zz_xz = buffer_1100_ppdd[1616];

    auto g_y_y_0_0_z_z_zz_yy = buffer_1100_ppdd[1617];

    auto g_y_y_0_0_z_z_zz_yz = buffer_1100_ppdd[1618];

    auto g_y_y_0_0_z_z_zz_zz = buffer_1100_ppdd[1619];

    auto g_y_z_0_0_x_x_xx_xx = buffer_1100_ppdd[1620];

    auto g_y_z_0_0_x_x_xx_xy = buffer_1100_ppdd[1621];

    auto g_y_z_0_0_x_x_xx_xz = buffer_1100_ppdd[1622];

    auto g_y_z_0_0_x_x_xx_yy = buffer_1100_ppdd[1623];

    auto g_y_z_0_0_x_x_xx_yz = buffer_1100_ppdd[1624];

    auto g_y_z_0_0_x_x_xx_zz = buffer_1100_ppdd[1625];

    auto g_y_z_0_0_x_x_xy_xx = buffer_1100_ppdd[1626];

    auto g_y_z_0_0_x_x_xy_xy = buffer_1100_ppdd[1627];

    auto g_y_z_0_0_x_x_xy_xz = buffer_1100_ppdd[1628];

    auto g_y_z_0_0_x_x_xy_yy = buffer_1100_ppdd[1629];

    auto g_y_z_0_0_x_x_xy_yz = buffer_1100_ppdd[1630];

    auto g_y_z_0_0_x_x_xy_zz = buffer_1100_ppdd[1631];

    auto g_y_z_0_0_x_x_xz_xx = buffer_1100_ppdd[1632];

    auto g_y_z_0_0_x_x_xz_xy = buffer_1100_ppdd[1633];

    auto g_y_z_0_0_x_x_xz_xz = buffer_1100_ppdd[1634];

    auto g_y_z_0_0_x_x_xz_yy = buffer_1100_ppdd[1635];

    auto g_y_z_0_0_x_x_xz_yz = buffer_1100_ppdd[1636];

    auto g_y_z_0_0_x_x_xz_zz = buffer_1100_ppdd[1637];

    auto g_y_z_0_0_x_x_yy_xx = buffer_1100_ppdd[1638];

    auto g_y_z_0_0_x_x_yy_xy = buffer_1100_ppdd[1639];

    auto g_y_z_0_0_x_x_yy_xz = buffer_1100_ppdd[1640];

    auto g_y_z_0_0_x_x_yy_yy = buffer_1100_ppdd[1641];

    auto g_y_z_0_0_x_x_yy_yz = buffer_1100_ppdd[1642];

    auto g_y_z_0_0_x_x_yy_zz = buffer_1100_ppdd[1643];

    auto g_y_z_0_0_x_x_yz_xx = buffer_1100_ppdd[1644];

    auto g_y_z_0_0_x_x_yz_xy = buffer_1100_ppdd[1645];

    auto g_y_z_0_0_x_x_yz_xz = buffer_1100_ppdd[1646];

    auto g_y_z_0_0_x_x_yz_yy = buffer_1100_ppdd[1647];

    auto g_y_z_0_0_x_x_yz_yz = buffer_1100_ppdd[1648];

    auto g_y_z_0_0_x_x_yz_zz = buffer_1100_ppdd[1649];

    auto g_y_z_0_0_x_x_zz_xx = buffer_1100_ppdd[1650];

    auto g_y_z_0_0_x_x_zz_xy = buffer_1100_ppdd[1651];

    auto g_y_z_0_0_x_x_zz_xz = buffer_1100_ppdd[1652];

    auto g_y_z_0_0_x_x_zz_yy = buffer_1100_ppdd[1653];

    auto g_y_z_0_0_x_x_zz_yz = buffer_1100_ppdd[1654];

    auto g_y_z_0_0_x_x_zz_zz = buffer_1100_ppdd[1655];

    auto g_y_z_0_0_x_y_xx_xx = buffer_1100_ppdd[1656];

    auto g_y_z_0_0_x_y_xx_xy = buffer_1100_ppdd[1657];

    auto g_y_z_0_0_x_y_xx_xz = buffer_1100_ppdd[1658];

    auto g_y_z_0_0_x_y_xx_yy = buffer_1100_ppdd[1659];

    auto g_y_z_0_0_x_y_xx_yz = buffer_1100_ppdd[1660];

    auto g_y_z_0_0_x_y_xx_zz = buffer_1100_ppdd[1661];

    auto g_y_z_0_0_x_y_xy_xx = buffer_1100_ppdd[1662];

    auto g_y_z_0_0_x_y_xy_xy = buffer_1100_ppdd[1663];

    auto g_y_z_0_0_x_y_xy_xz = buffer_1100_ppdd[1664];

    auto g_y_z_0_0_x_y_xy_yy = buffer_1100_ppdd[1665];

    auto g_y_z_0_0_x_y_xy_yz = buffer_1100_ppdd[1666];

    auto g_y_z_0_0_x_y_xy_zz = buffer_1100_ppdd[1667];

    auto g_y_z_0_0_x_y_xz_xx = buffer_1100_ppdd[1668];

    auto g_y_z_0_0_x_y_xz_xy = buffer_1100_ppdd[1669];

    auto g_y_z_0_0_x_y_xz_xz = buffer_1100_ppdd[1670];

    auto g_y_z_0_0_x_y_xz_yy = buffer_1100_ppdd[1671];

    auto g_y_z_0_0_x_y_xz_yz = buffer_1100_ppdd[1672];

    auto g_y_z_0_0_x_y_xz_zz = buffer_1100_ppdd[1673];

    auto g_y_z_0_0_x_y_yy_xx = buffer_1100_ppdd[1674];

    auto g_y_z_0_0_x_y_yy_xy = buffer_1100_ppdd[1675];

    auto g_y_z_0_0_x_y_yy_xz = buffer_1100_ppdd[1676];

    auto g_y_z_0_0_x_y_yy_yy = buffer_1100_ppdd[1677];

    auto g_y_z_0_0_x_y_yy_yz = buffer_1100_ppdd[1678];

    auto g_y_z_0_0_x_y_yy_zz = buffer_1100_ppdd[1679];

    auto g_y_z_0_0_x_y_yz_xx = buffer_1100_ppdd[1680];

    auto g_y_z_0_0_x_y_yz_xy = buffer_1100_ppdd[1681];

    auto g_y_z_0_0_x_y_yz_xz = buffer_1100_ppdd[1682];

    auto g_y_z_0_0_x_y_yz_yy = buffer_1100_ppdd[1683];

    auto g_y_z_0_0_x_y_yz_yz = buffer_1100_ppdd[1684];

    auto g_y_z_0_0_x_y_yz_zz = buffer_1100_ppdd[1685];

    auto g_y_z_0_0_x_y_zz_xx = buffer_1100_ppdd[1686];

    auto g_y_z_0_0_x_y_zz_xy = buffer_1100_ppdd[1687];

    auto g_y_z_0_0_x_y_zz_xz = buffer_1100_ppdd[1688];

    auto g_y_z_0_0_x_y_zz_yy = buffer_1100_ppdd[1689];

    auto g_y_z_0_0_x_y_zz_yz = buffer_1100_ppdd[1690];

    auto g_y_z_0_0_x_y_zz_zz = buffer_1100_ppdd[1691];

    auto g_y_z_0_0_x_z_xx_xx = buffer_1100_ppdd[1692];

    auto g_y_z_0_0_x_z_xx_xy = buffer_1100_ppdd[1693];

    auto g_y_z_0_0_x_z_xx_xz = buffer_1100_ppdd[1694];

    auto g_y_z_0_0_x_z_xx_yy = buffer_1100_ppdd[1695];

    auto g_y_z_0_0_x_z_xx_yz = buffer_1100_ppdd[1696];

    auto g_y_z_0_0_x_z_xx_zz = buffer_1100_ppdd[1697];

    auto g_y_z_0_0_x_z_xy_xx = buffer_1100_ppdd[1698];

    auto g_y_z_0_0_x_z_xy_xy = buffer_1100_ppdd[1699];

    auto g_y_z_0_0_x_z_xy_xz = buffer_1100_ppdd[1700];

    auto g_y_z_0_0_x_z_xy_yy = buffer_1100_ppdd[1701];

    auto g_y_z_0_0_x_z_xy_yz = buffer_1100_ppdd[1702];

    auto g_y_z_0_0_x_z_xy_zz = buffer_1100_ppdd[1703];

    auto g_y_z_0_0_x_z_xz_xx = buffer_1100_ppdd[1704];

    auto g_y_z_0_0_x_z_xz_xy = buffer_1100_ppdd[1705];

    auto g_y_z_0_0_x_z_xz_xz = buffer_1100_ppdd[1706];

    auto g_y_z_0_0_x_z_xz_yy = buffer_1100_ppdd[1707];

    auto g_y_z_0_0_x_z_xz_yz = buffer_1100_ppdd[1708];

    auto g_y_z_0_0_x_z_xz_zz = buffer_1100_ppdd[1709];

    auto g_y_z_0_0_x_z_yy_xx = buffer_1100_ppdd[1710];

    auto g_y_z_0_0_x_z_yy_xy = buffer_1100_ppdd[1711];

    auto g_y_z_0_0_x_z_yy_xz = buffer_1100_ppdd[1712];

    auto g_y_z_0_0_x_z_yy_yy = buffer_1100_ppdd[1713];

    auto g_y_z_0_0_x_z_yy_yz = buffer_1100_ppdd[1714];

    auto g_y_z_0_0_x_z_yy_zz = buffer_1100_ppdd[1715];

    auto g_y_z_0_0_x_z_yz_xx = buffer_1100_ppdd[1716];

    auto g_y_z_0_0_x_z_yz_xy = buffer_1100_ppdd[1717];

    auto g_y_z_0_0_x_z_yz_xz = buffer_1100_ppdd[1718];

    auto g_y_z_0_0_x_z_yz_yy = buffer_1100_ppdd[1719];

    auto g_y_z_0_0_x_z_yz_yz = buffer_1100_ppdd[1720];

    auto g_y_z_0_0_x_z_yz_zz = buffer_1100_ppdd[1721];

    auto g_y_z_0_0_x_z_zz_xx = buffer_1100_ppdd[1722];

    auto g_y_z_0_0_x_z_zz_xy = buffer_1100_ppdd[1723];

    auto g_y_z_0_0_x_z_zz_xz = buffer_1100_ppdd[1724];

    auto g_y_z_0_0_x_z_zz_yy = buffer_1100_ppdd[1725];

    auto g_y_z_0_0_x_z_zz_yz = buffer_1100_ppdd[1726];

    auto g_y_z_0_0_x_z_zz_zz = buffer_1100_ppdd[1727];

    auto g_y_z_0_0_y_x_xx_xx = buffer_1100_ppdd[1728];

    auto g_y_z_0_0_y_x_xx_xy = buffer_1100_ppdd[1729];

    auto g_y_z_0_0_y_x_xx_xz = buffer_1100_ppdd[1730];

    auto g_y_z_0_0_y_x_xx_yy = buffer_1100_ppdd[1731];

    auto g_y_z_0_0_y_x_xx_yz = buffer_1100_ppdd[1732];

    auto g_y_z_0_0_y_x_xx_zz = buffer_1100_ppdd[1733];

    auto g_y_z_0_0_y_x_xy_xx = buffer_1100_ppdd[1734];

    auto g_y_z_0_0_y_x_xy_xy = buffer_1100_ppdd[1735];

    auto g_y_z_0_0_y_x_xy_xz = buffer_1100_ppdd[1736];

    auto g_y_z_0_0_y_x_xy_yy = buffer_1100_ppdd[1737];

    auto g_y_z_0_0_y_x_xy_yz = buffer_1100_ppdd[1738];

    auto g_y_z_0_0_y_x_xy_zz = buffer_1100_ppdd[1739];

    auto g_y_z_0_0_y_x_xz_xx = buffer_1100_ppdd[1740];

    auto g_y_z_0_0_y_x_xz_xy = buffer_1100_ppdd[1741];

    auto g_y_z_0_0_y_x_xz_xz = buffer_1100_ppdd[1742];

    auto g_y_z_0_0_y_x_xz_yy = buffer_1100_ppdd[1743];

    auto g_y_z_0_0_y_x_xz_yz = buffer_1100_ppdd[1744];

    auto g_y_z_0_0_y_x_xz_zz = buffer_1100_ppdd[1745];

    auto g_y_z_0_0_y_x_yy_xx = buffer_1100_ppdd[1746];

    auto g_y_z_0_0_y_x_yy_xy = buffer_1100_ppdd[1747];

    auto g_y_z_0_0_y_x_yy_xz = buffer_1100_ppdd[1748];

    auto g_y_z_0_0_y_x_yy_yy = buffer_1100_ppdd[1749];

    auto g_y_z_0_0_y_x_yy_yz = buffer_1100_ppdd[1750];

    auto g_y_z_0_0_y_x_yy_zz = buffer_1100_ppdd[1751];

    auto g_y_z_0_0_y_x_yz_xx = buffer_1100_ppdd[1752];

    auto g_y_z_0_0_y_x_yz_xy = buffer_1100_ppdd[1753];

    auto g_y_z_0_0_y_x_yz_xz = buffer_1100_ppdd[1754];

    auto g_y_z_0_0_y_x_yz_yy = buffer_1100_ppdd[1755];

    auto g_y_z_0_0_y_x_yz_yz = buffer_1100_ppdd[1756];

    auto g_y_z_0_0_y_x_yz_zz = buffer_1100_ppdd[1757];

    auto g_y_z_0_0_y_x_zz_xx = buffer_1100_ppdd[1758];

    auto g_y_z_0_0_y_x_zz_xy = buffer_1100_ppdd[1759];

    auto g_y_z_0_0_y_x_zz_xz = buffer_1100_ppdd[1760];

    auto g_y_z_0_0_y_x_zz_yy = buffer_1100_ppdd[1761];

    auto g_y_z_0_0_y_x_zz_yz = buffer_1100_ppdd[1762];

    auto g_y_z_0_0_y_x_zz_zz = buffer_1100_ppdd[1763];

    auto g_y_z_0_0_y_y_xx_xx = buffer_1100_ppdd[1764];

    auto g_y_z_0_0_y_y_xx_xy = buffer_1100_ppdd[1765];

    auto g_y_z_0_0_y_y_xx_xz = buffer_1100_ppdd[1766];

    auto g_y_z_0_0_y_y_xx_yy = buffer_1100_ppdd[1767];

    auto g_y_z_0_0_y_y_xx_yz = buffer_1100_ppdd[1768];

    auto g_y_z_0_0_y_y_xx_zz = buffer_1100_ppdd[1769];

    auto g_y_z_0_0_y_y_xy_xx = buffer_1100_ppdd[1770];

    auto g_y_z_0_0_y_y_xy_xy = buffer_1100_ppdd[1771];

    auto g_y_z_0_0_y_y_xy_xz = buffer_1100_ppdd[1772];

    auto g_y_z_0_0_y_y_xy_yy = buffer_1100_ppdd[1773];

    auto g_y_z_0_0_y_y_xy_yz = buffer_1100_ppdd[1774];

    auto g_y_z_0_0_y_y_xy_zz = buffer_1100_ppdd[1775];

    auto g_y_z_0_0_y_y_xz_xx = buffer_1100_ppdd[1776];

    auto g_y_z_0_0_y_y_xz_xy = buffer_1100_ppdd[1777];

    auto g_y_z_0_0_y_y_xz_xz = buffer_1100_ppdd[1778];

    auto g_y_z_0_0_y_y_xz_yy = buffer_1100_ppdd[1779];

    auto g_y_z_0_0_y_y_xz_yz = buffer_1100_ppdd[1780];

    auto g_y_z_0_0_y_y_xz_zz = buffer_1100_ppdd[1781];

    auto g_y_z_0_0_y_y_yy_xx = buffer_1100_ppdd[1782];

    auto g_y_z_0_0_y_y_yy_xy = buffer_1100_ppdd[1783];

    auto g_y_z_0_0_y_y_yy_xz = buffer_1100_ppdd[1784];

    auto g_y_z_0_0_y_y_yy_yy = buffer_1100_ppdd[1785];

    auto g_y_z_0_0_y_y_yy_yz = buffer_1100_ppdd[1786];

    auto g_y_z_0_0_y_y_yy_zz = buffer_1100_ppdd[1787];

    auto g_y_z_0_0_y_y_yz_xx = buffer_1100_ppdd[1788];

    auto g_y_z_0_0_y_y_yz_xy = buffer_1100_ppdd[1789];

    auto g_y_z_0_0_y_y_yz_xz = buffer_1100_ppdd[1790];

    auto g_y_z_0_0_y_y_yz_yy = buffer_1100_ppdd[1791];

    auto g_y_z_0_0_y_y_yz_yz = buffer_1100_ppdd[1792];

    auto g_y_z_0_0_y_y_yz_zz = buffer_1100_ppdd[1793];

    auto g_y_z_0_0_y_y_zz_xx = buffer_1100_ppdd[1794];

    auto g_y_z_0_0_y_y_zz_xy = buffer_1100_ppdd[1795];

    auto g_y_z_0_0_y_y_zz_xz = buffer_1100_ppdd[1796];

    auto g_y_z_0_0_y_y_zz_yy = buffer_1100_ppdd[1797];

    auto g_y_z_0_0_y_y_zz_yz = buffer_1100_ppdd[1798];

    auto g_y_z_0_0_y_y_zz_zz = buffer_1100_ppdd[1799];

    auto g_y_z_0_0_y_z_xx_xx = buffer_1100_ppdd[1800];

    auto g_y_z_0_0_y_z_xx_xy = buffer_1100_ppdd[1801];

    auto g_y_z_0_0_y_z_xx_xz = buffer_1100_ppdd[1802];

    auto g_y_z_0_0_y_z_xx_yy = buffer_1100_ppdd[1803];

    auto g_y_z_0_0_y_z_xx_yz = buffer_1100_ppdd[1804];

    auto g_y_z_0_0_y_z_xx_zz = buffer_1100_ppdd[1805];

    auto g_y_z_0_0_y_z_xy_xx = buffer_1100_ppdd[1806];

    auto g_y_z_0_0_y_z_xy_xy = buffer_1100_ppdd[1807];

    auto g_y_z_0_0_y_z_xy_xz = buffer_1100_ppdd[1808];

    auto g_y_z_0_0_y_z_xy_yy = buffer_1100_ppdd[1809];

    auto g_y_z_0_0_y_z_xy_yz = buffer_1100_ppdd[1810];

    auto g_y_z_0_0_y_z_xy_zz = buffer_1100_ppdd[1811];

    auto g_y_z_0_0_y_z_xz_xx = buffer_1100_ppdd[1812];

    auto g_y_z_0_0_y_z_xz_xy = buffer_1100_ppdd[1813];

    auto g_y_z_0_0_y_z_xz_xz = buffer_1100_ppdd[1814];

    auto g_y_z_0_0_y_z_xz_yy = buffer_1100_ppdd[1815];

    auto g_y_z_0_0_y_z_xz_yz = buffer_1100_ppdd[1816];

    auto g_y_z_0_0_y_z_xz_zz = buffer_1100_ppdd[1817];

    auto g_y_z_0_0_y_z_yy_xx = buffer_1100_ppdd[1818];

    auto g_y_z_0_0_y_z_yy_xy = buffer_1100_ppdd[1819];

    auto g_y_z_0_0_y_z_yy_xz = buffer_1100_ppdd[1820];

    auto g_y_z_0_0_y_z_yy_yy = buffer_1100_ppdd[1821];

    auto g_y_z_0_0_y_z_yy_yz = buffer_1100_ppdd[1822];

    auto g_y_z_0_0_y_z_yy_zz = buffer_1100_ppdd[1823];

    auto g_y_z_0_0_y_z_yz_xx = buffer_1100_ppdd[1824];

    auto g_y_z_0_0_y_z_yz_xy = buffer_1100_ppdd[1825];

    auto g_y_z_0_0_y_z_yz_xz = buffer_1100_ppdd[1826];

    auto g_y_z_0_0_y_z_yz_yy = buffer_1100_ppdd[1827];

    auto g_y_z_0_0_y_z_yz_yz = buffer_1100_ppdd[1828];

    auto g_y_z_0_0_y_z_yz_zz = buffer_1100_ppdd[1829];

    auto g_y_z_0_0_y_z_zz_xx = buffer_1100_ppdd[1830];

    auto g_y_z_0_0_y_z_zz_xy = buffer_1100_ppdd[1831];

    auto g_y_z_0_0_y_z_zz_xz = buffer_1100_ppdd[1832];

    auto g_y_z_0_0_y_z_zz_yy = buffer_1100_ppdd[1833];

    auto g_y_z_0_0_y_z_zz_yz = buffer_1100_ppdd[1834];

    auto g_y_z_0_0_y_z_zz_zz = buffer_1100_ppdd[1835];

    auto g_y_z_0_0_z_x_xx_xx = buffer_1100_ppdd[1836];

    auto g_y_z_0_0_z_x_xx_xy = buffer_1100_ppdd[1837];

    auto g_y_z_0_0_z_x_xx_xz = buffer_1100_ppdd[1838];

    auto g_y_z_0_0_z_x_xx_yy = buffer_1100_ppdd[1839];

    auto g_y_z_0_0_z_x_xx_yz = buffer_1100_ppdd[1840];

    auto g_y_z_0_0_z_x_xx_zz = buffer_1100_ppdd[1841];

    auto g_y_z_0_0_z_x_xy_xx = buffer_1100_ppdd[1842];

    auto g_y_z_0_0_z_x_xy_xy = buffer_1100_ppdd[1843];

    auto g_y_z_0_0_z_x_xy_xz = buffer_1100_ppdd[1844];

    auto g_y_z_0_0_z_x_xy_yy = buffer_1100_ppdd[1845];

    auto g_y_z_0_0_z_x_xy_yz = buffer_1100_ppdd[1846];

    auto g_y_z_0_0_z_x_xy_zz = buffer_1100_ppdd[1847];

    auto g_y_z_0_0_z_x_xz_xx = buffer_1100_ppdd[1848];

    auto g_y_z_0_0_z_x_xz_xy = buffer_1100_ppdd[1849];

    auto g_y_z_0_0_z_x_xz_xz = buffer_1100_ppdd[1850];

    auto g_y_z_0_0_z_x_xz_yy = buffer_1100_ppdd[1851];

    auto g_y_z_0_0_z_x_xz_yz = buffer_1100_ppdd[1852];

    auto g_y_z_0_0_z_x_xz_zz = buffer_1100_ppdd[1853];

    auto g_y_z_0_0_z_x_yy_xx = buffer_1100_ppdd[1854];

    auto g_y_z_0_0_z_x_yy_xy = buffer_1100_ppdd[1855];

    auto g_y_z_0_0_z_x_yy_xz = buffer_1100_ppdd[1856];

    auto g_y_z_0_0_z_x_yy_yy = buffer_1100_ppdd[1857];

    auto g_y_z_0_0_z_x_yy_yz = buffer_1100_ppdd[1858];

    auto g_y_z_0_0_z_x_yy_zz = buffer_1100_ppdd[1859];

    auto g_y_z_0_0_z_x_yz_xx = buffer_1100_ppdd[1860];

    auto g_y_z_0_0_z_x_yz_xy = buffer_1100_ppdd[1861];

    auto g_y_z_0_0_z_x_yz_xz = buffer_1100_ppdd[1862];

    auto g_y_z_0_0_z_x_yz_yy = buffer_1100_ppdd[1863];

    auto g_y_z_0_0_z_x_yz_yz = buffer_1100_ppdd[1864];

    auto g_y_z_0_0_z_x_yz_zz = buffer_1100_ppdd[1865];

    auto g_y_z_0_0_z_x_zz_xx = buffer_1100_ppdd[1866];

    auto g_y_z_0_0_z_x_zz_xy = buffer_1100_ppdd[1867];

    auto g_y_z_0_0_z_x_zz_xz = buffer_1100_ppdd[1868];

    auto g_y_z_0_0_z_x_zz_yy = buffer_1100_ppdd[1869];

    auto g_y_z_0_0_z_x_zz_yz = buffer_1100_ppdd[1870];

    auto g_y_z_0_0_z_x_zz_zz = buffer_1100_ppdd[1871];

    auto g_y_z_0_0_z_y_xx_xx = buffer_1100_ppdd[1872];

    auto g_y_z_0_0_z_y_xx_xy = buffer_1100_ppdd[1873];

    auto g_y_z_0_0_z_y_xx_xz = buffer_1100_ppdd[1874];

    auto g_y_z_0_0_z_y_xx_yy = buffer_1100_ppdd[1875];

    auto g_y_z_0_0_z_y_xx_yz = buffer_1100_ppdd[1876];

    auto g_y_z_0_0_z_y_xx_zz = buffer_1100_ppdd[1877];

    auto g_y_z_0_0_z_y_xy_xx = buffer_1100_ppdd[1878];

    auto g_y_z_0_0_z_y_xy_xy = buffer_1100_ppdd[1879];

    auto g_y_z_0_0_z_y_xy_xz = buffer_1100_ppdd[1880];

    auto g_y_z_0_0_z_y_xy_yy = buffer_1100_ppdd[1881];

    auto g_y_z_0_0_z_y_xy_yz = buffer_1100_ppdd[1882];

    auto g_y_z_0_0_z_y_xy_zz = buffer_1100_ppdd[1883];

    auto g_y_z_0_0_z_y_xz_xx = buffer_1100_ppdd[1884];

    auto g_y_z_0_0_z_y_xz_xy = buffer_1100_ppdd[1885];

    auto g_y_z_0_0_z_y_xz_xz = buffer_1100_ppdd[1886];

    auto g_y_z_0_0_z_y_xz_yy = buffer_1100_ppdd[1887];

    auto g_y_z_0_0_z_y_xz_yz = buffer_1100_ppdd[1888];

    auto g_y_z_0_0_z_y_xz_zz = buffer_1100_ppdd[1889];

    auto g_y_z_0_0_z_y_yy_xx = buffer_1100_ppdd[1890];

    auto g_y_z_0_0_z_y_yy_xy = buffer_1100_ppdd[1891];

    auto g_y_z_0_0_z_y_yy_xz = buffer_1100_ppdd[1892];

    auto g_y_z_0_0_z_y_yy_yy = buffer_1100_ppdd[1893];

    auto g_y_z_0_0_z_y_yy_yz = buffer_1100_ppdd[1894];

    auto g_y_z_0_0_z_y_yy_zz = buffer_1100_ppdd[1895];

    auto g_y_z_0_0_z_y_yz_xx = buffer_1100_ppdd[1896];

    auto g_y_z_0_0_z_y_yz_xy = buffer_1100_ppdd[1897];

    auto g_y_z_0_0_z_y_yz_xz = buffer_1100_ppdd[1898];

    auto g_y_z_0_0_z_y_yz_yy = buffer_1100_ppdd[1899];

    auto g_y_z_0_0_z_y_yz_yz = buffer_1100_ppdd[1900];

    auto g_y_z_0_0_z_y_yz_zz = buffer_1100_ppdd[1901];

    auto g_y_z_0_0_z_y_zz_xx = buffer_1100_ppdd[1902];

    auto g_y_z_0_0_z_y_zz_xy = buffer_1100_ppdd[1903];

    auto g_y_z_0_0_z_y_zz_xz = buffer_1100_ppdd[1904];

    auto g_y_z_0_0_z_y_zz_yy = buffer_1100_ppdd[1905];

    auto g_y_z_0_0_z_y_zz_yz = buffer_1100_ppdd[1906];

    auto g_y_z_0_0_z_y_zz_zz = buffer_1100_ppdd[1907];

    auto g_y_z_0_0_z_z_xx_xx = buffer_1100_ppdd[1908];

    auto g_y_z_0_0_z_z_xx_xy = buffer_1100_ppdd[1909];

    auto g_y_z_0_0_z_z_xx_xz = buffer_1100_ppdd[1910];

    auto g_y_z_0_0_z_z_xx_yy = buffer_1100_ppdd[1911];

    auto g_y_z_0_0_z_z_xx_yz = buffer_1100_ppdd[1912];

    auto g_y_z_0_0_z_z_xx_zz = buffer_1100_ppdd[1913];

    auto g_y_z_0_0_z_z_xy_xx = buffer_1100_ppdd[1914];

    auto g_y_z_0_0_z_z_xy_xy = buffer_1100_ppdd[1915];

    auto g_y_z_0_0_z_z_xy_xz = buffer_1100_ppdd[1916];

    auto g_y_z_0_0_z_z_xy_yy = buffer_1100_ppdd[1917];

    auto g_y_z_0_0_z_z_xy_yz = buffer_1100_ppdd[1918];

    auto g_y_z_0_0_z_z_xy_zz = buffer_1100_ppdd[1919];

    auto g_y_z_0_0_z_z_xz_xx = buffer_1100_ppdd[1920];

    auto g_y_z_0_0_z_z_xz_xy = buffer_1100_ppdd[1921];

    auto g_y_z_0_0_z_z_xz_xz = buffer_1100_ppdd[1922];

    auto g_y_z_0_0_z_z_xz_yy = buffer_1100_ppdd[1923];

    auto g_y_z_0_0_z_z_xz_yz = buffer_1100_ppdd[1924];

    auto g_y_z_0_0_z_z_xz_zz = buffer_1100_ppdd[1925];

    auto g_y_z_0_0_z_z_yy_xx = buffer_1100_ppdd[1926];

    auto g_y_z_0_0_z_z_yy_xy = buffer_1100_ppdd[1927];

    auto g_y_z_0_0_z_z_yy_xz = buffer_1100_ppdd[1928];

    auto g_y_z_0_0_z_z_yy_yy = buffer_1100_ppdd[1929];

    auto g_y_z_0_0_z_z_yy_yz = buffer_1100_ppdd[1930];

    auto g_y_z_0_0_z_z_yy_zz = buffer_1100_ppdd[1931];

    auto g_y_z_0_0_z_z_yz_xx = buffer_1100_ppdd[1932];

    auto g_y_z_0_0_z_z_yz_xy = buffer_1100_ppdd[1933];

    auto g_y_z_0_0_z_z_yz_xz = buffer_1100_ppdd[1934];

    auto g_y_z_0_0_z_z_yz_yy = buffer_1100_ppdd[1935];

    auto g_y_z_0_0_z_z_yz_yz = buffer_1100_ppdd[1936];

    auto g_y_z_0_0_z_z_yz_zz = buffer_1100_ppdd[1937];

    auto g_y_z_0_0_z_z_zz_xx = buffer_1100_ppdd[1938];

    auto g_y_z_0_0_z_z_zz_xy = buffer_1100_ppdd[1939];

    auto g_y_z_0_0_z_z_zz_xz = buffer_1100_ppdd[1940];

    auto g_y_z_0_0_z_z_zz_yy = buffer_1100_ppdd[1941];

    auto g_y_z_0_0_z_z_zz_yz = buffer_1100_ppdd[1942];

    auto g_y_z_0_0_z_z_zz_zz = buffer_1100_ppdd[1943];

    auto g_z_x_0_0_x_x_xx_xx = buffer_1100_ppdd[1944];

    auto g_z_x_0_0_x_x_xx_xy = buffer_1100_ppdd[1945];

    auto g_z_x_0_0_x_x_xx_xz = buffer_1100_ppdd[1946];

    auto g_z_x_0_0_x_x_xx_yy = buffer_1100_ppdd[1947];

    auto g_z_x_0_0_x_x_xx_yz = buffer_1100_ppdd[1948];

    auto g_z_x_0_0_x_x_xx_zz = buffer_1100_ppdd[1949];

    auto g_z_x_0_0_x_x_xy_xx = buffer_1100_ppdd[1950];

    auto g_z_x_0_0_x_x_xy_xy = buffer_1100_ppdd[1951];

    auto g_z_x_0_0_x_x_xy_xz = buffer_1100_ppdd[1952];

    auto g_z_x_0_0_x_x_xy_yy = buffer_1100_ppdd[1953];

    auto g_z_x_0_0_x_x_xy_yz = buffer_1100_ppdd[1954];

    auto g_z_x_0_0_x_x_xy_zz = buffer_1100_ppdd[1955];

    auto g_z_x_0_0_x_x_xz_xx = buffer_1100_ppdd[1956];

    auto g_z_x_0_0_x_x_xz_xy = buffer_1100_ppdd[1957];

    auto g_z_x_0_0_x_x_xz_xz = buffer_1100_ppdd[1958];

    auto g_z_x_0_0_x_x_xz_yy = buffer_1100_ppdd[1959];

    auto g_z_x_0_0_x_x_xz_yz = buffer_1100_ppdd[1960];

    auto g_z_x_0_0_x_x_xz_zz = buffer_1100_ppdd[1961];

    auto g_z_x_0_0_x_x_yy_xx = buffer_1100_ppdd[1962];

    auto g_z_x_0_0_x_x_yy_xy = buffer_1100_ppdd[1963];

    auto g_z_x_0_0_x_x_yy_xz = buffer_1100_ppdd[1964];

    auto g_z_x_0_0_x_x_yy_yy = buffer_1100_ppdd[1965];

    auto g_z_x_0_0_x_x_yy_yz = buffer_1100_ppdd[1966];

    auto g_z_x_0_0_x_x_yy_zz = buffer_1100_ppdd[1967];

    auto g_z_x_0_0_x_x_yz_xx = buffer_1100_ppdd[1968];

    auto g_z_x_0_0_x_x_yz_xy = buffer_1100_ppdd[1969];

    auto g_z_x_0_0_x_x_yz_xz = buffer_1100_ppdd[1970];

    auto g_z_x_0_0_x_x_yz_yy = buffer_1100_ppdd[1971];

    auto g_z_x_0_0_x_x_yz_yz = buffer_1100_ppdd[1972];

    auto g_z_x_0_0_x_x_yz_zz = buffer_1100_ppdd[1973];

    auto g_z_x_0_0_x_x_zz_xx = buffer_1100_ppdd[1974];

    auto g_z_x_0_0_x_x_zz_xy = buffer_1100_ppdd[1975];

    auto g_z_x_0_0_x_x_zz_xz = buffer_1100_ppdd[1976];

    auto g_z_x_0_0_x_x_zz_yy = buffer_1100_ppdd[1977];

    auto g_z_x_0_0_x_x_zz_yz = buffer_1100_ppdd[1978];

    auto g_z_x_0_0_x_x_zz_zz = buffer_1100_ppdd[1979];

    auto g_z_x_0_0_x_y_xx_xx = buffer_1100_ppdd[1980];

    auto g_z_x_0_0_x_y_xx_xy = buffer_1100_ppdd[1981];

    auto g_z_x_0_0_x_y_xx_xz = buffer_1100_ppdd[1982];

    auto g_z_x_0_0_x_y_xx_yy = buffer_1100_ppdd[1983];

    auto g_z_x_0_0_x_y_xx_yz = buffer_1100_ppdd[1984];

    auto g_z_x_0_0_x_y_xx_zz = buffer_1100_ppdd[1985];

    auto g_z_x_0_0_x_y_xy_xx = buffer_1100_ppdd[1986];

    auto g_z_x_0_0_x_y_xy_xy = buffer_1100_ppdd[1987];

    auto g_z_x_0_0_x_y_xy_xz = buffer_1100_ppdd[1988];

    auto g_z_x_0_0_x_y_xy_yy = buffer_1100_ppdd[1989];

    auto g_z_x_0_0_x_y_xy_yz = buffer_1100_ppdd[1990];

    auto g_z_x_0_0_x_y_xy_zz = buffer_1100_ppdd[1991];

    auto g_z_x_0_0_x_y_xz_xx = buffer_1100_ppdd[1992];

    auto g_z_x_0_0_x_y_xz_xy = buffer_1100_ppdd[1993];

    auto g_z_x_0_0_x_y_xz_xz = buffer_1100_ppdd[1994];

    auto g_z_x_0_0_x_y_xz_yy = buffer_1100_ppdd[1995];

    auto g_z_x_0_0_x_y_xz_yz = buffer_1100_ppdd[1996];

    auto g_z_x_0_0_x_y_xz_zz = buffer_1100_ppdd[1997];

    auto g_z_x_0_0_x_y_yy_xx = buffer_1100_ppdd[1998];

    auto g_z_x_0_0_x_y_yy_xy = buffer_1100_ppdd[1999];

    auto g_z_x_0_0_x_y_yy_xz = buffer_1100_ppdd[2000];

    auto g_z_x_0_0_x_y_yy_yy = buffer_1100_ppdd[2001];

    auto g_z_x_0_0_x_y_yy_yz = buffer_1100_ppdd[2002];

    auto g_z_x_0_0_x_y_yy_zz = buffer_1100_ppdd[2003];

    auto g_z_x_0_0_x_y_yz_xx = buffer_1100_ppdd[2004];

    auto g_z_x_0_0_x_y_yz_xy = buffer_1100_ppdd[2005];

    auto g_z_x_0_0_x_y_yz_xz = buffer_1100_ppdd[2006];

    auto g_z_x_0_0_x_y_yz_yy = buffer_1100_ppdd[2007];

    auto g_z_x_0_0_x_y_yz_yz = buffer_1100_ppdd[2008];

    auto g_z_x_0_0_x_y_yz_zz = buffer_1100_ppdd[2009];

    auto g_z_x_0_0_x_y_zz_xx = buffer_1100_ppdd[2010];

    auto g_z_x_0_0_x_y_zz_xy = buffer_1100_ppdd[2011];

    auto g_z_x_0_0_x_y_zz_xz = buffer_1100_ppdd[2012];

    auto g_z_x_0_0_x_y_zz_yy = buffer_1100_ppdd[2013];

    auto g_z_x_0_0_x_y_zz_yz = buffer_1100_ppdd[2014];

    auto g_z_x_0_0_x_y_zz_zz = buffer_1100_ppdd[2015];

    auto g_z_x_0_0_x_z_xx_xx = buffer_1100_ppdd[2016];

    auto g_z_x_0_0_x_z_xx_xy = buffer_1100_ppdd[2017];

    auto g_z_x_0_0_x_z_xx_xz = buffer_1100_ppdd[2018];

    auto g_z_x_0_0_x_z_xx_yy = buffer_1100_ppdd[2019];

    auto g_z_x_0_0_x_z_xx_yz = buffer_1100_ppdd[2020];

    auto g_z_x_0_0_x_z_xx_zz = buffer_1100_ppdd[2021];

    auto g_z_x_0_0_x_z_xy_xx = buffer_1100_ppdd[2022];

    auto g_z_x_0_0_x_z_xy_xy = buffer_1100_ppdd[2023];

    auto g_z_x_0_0_x_z_xy_xz = buffer_1100_ppdd[2024];

    auto g_z_x_0_0_x_z_xy_yy = buffer_1100_ppdd[2025];

    auto g_z_x_0_0_x_z_xy_yz = buffer_1100_ppdd[2026];

    auto g_z_x_0_0_x_z_xy_zz = buffer_1100_ppdd[2027];

    auto g_z_x_0_0_x_z_xz_xx = buffer_1100_ppdd[2028];

    auto g_z_x_0_0_x_z_xz_xy = buffer_1100_ppdd[2029];

    auto g_z_x_0_0_x_z_xz_xz = buffer_1100_ppdd[2030];

    auto g_z_x_0_0_x_z_xz_yy = buffer_1100_ppdd[2031];

    auto g_z_x_0_0_x_z_xz_yz = buffer_1100_ppdd[2032];

    auto g_z_x_0_0_x_z_xz_zz = buffer_1100_ppdd[2033];

    auto g_z_x_0_0_x_z_yy_xx = buffer_1100_ppdd[2034];

    auto g_z_x_0_0_x_z_yy_xy = buffer_1100_ppdd[2035];

    auto g_z_x_0_0_x_z_yy_xz = buffer_1100_ppdd[2036];

    auto g_z_x_0_0_x_z_yy_yy = buffer_1100_ppdd[2037];

    auto g_z_x_0_0_x_z_yy_yz = buffer_1100_ppdd[2038];

    auto g_z_x_0_0_x_z_yy_zz = buffer_1100_ppdd[2039];

    auto g_z_x_0_0_x_z_yz_xx = buffer_1100_ppdd[2040];

    auto g_z_x_0_0_x_z_yz_xy = buffer_1100_ppdd[2041];

    auto g_z_x_0_0_x_z_yz_xz = buffer_1100_ppdd[2042];

    auto g_z_x_0_0_x_z_yz_yy = buffer_1100_ppdd[2043];

    auto g_z_x_0_0_x_z_yz_yz = buffer_1100_ppdd[2044];

    auto g_z_x_0_0_x_z_yz_zz = buffer_1100_ppdd[2045];

    auto g_z_x_0_0_x_z_zz_xx = buffer_1100_ppdd[2046];

    auto g_z_x_0_0_x_z_zz_xy = buffer_1100_ppdd[2047];

    auto g_z_x_0_0_x_z_zz_xz = buffer_1100_ppdd[2048];

    auto g_z_x_0_0_x_z_zz_yy = buffer_1100_ppdd[2049];

    auto g_z_x_0_0_x_z_zz_yz = buffer_1100_ppdd[2050];

    auto g_z_x_0_0_x_z_zz_zz = buffer_1100_ppdd[2051];

    auto g_z_x_0_0_y_x_xx_xx = buffer_1100_ppdd[2052];

    auto g_z_x_0_0_y_x_xx_xy = buffer_1100_ppdd[2053];

    auto g_z_x_0_0_y_x_xx_xz = buffer_1100_ppdd[2054];

    auto g_z_x_0_0_y_x_xx_yy = buffer_1100_ppdd[2055];

    auto g_z_x_0_0_y_x_xx_yz = buffer_1100_ppdd[2056];

    auto g_z_x_0_0_y_x_xx_zz = buffer_1100_ppdd[2057];

    auto g_z_x_0_0_y_x_xy_xx = buffer_1100_ppdd[2058];

    auto g_z_x_0_0_y_x_xy_xy = buffer_1100_ppdd[2059];

    auto g_z_x_0_0_y_x_xy_xz = buffer_1100_ppdd[2060];

    auto g_z_x_0_0_y_x_xy_yy = buffer_1100_ppdd[2061];

    auto g_z_x_0_0_y_x_xy_yz = buffer_1100_ppdd[2062];

    auto g_z_x_0_0_y_x_xy_zz = buffer_1100_ppdd[2063];

    auto g_z_x_0_0_y_x_xz_xx = buffer_1100_ppdd[2064];

    auto g_z_x_0_0_y_x_xz_xy = buffer_1100_ppdd[2065];

    auto g_z_x_0_0_y_x_xz_xz = buffer_1100_ppdd[2066];

    auto g_z_x_0_0_y_x_xz_yy = buffer_1100_ppdd[2067];

    auto g_z_x_0_0_y_x_xz_yz = buffer_1100_ppdd[2068];

    auto g_z_x_0_0_y_x_xz_zz = buffer_1100_ppdd[2069];

    auto g_z_x_0_0_y_x_yy_xx = buffer_1100_ppdd[2070];

    auto g_z_x_0_0_y_x_yy_xy = buffer_1100_ppdd[2071];

    auto g_z_x_0_0_y_x_yy_xz = buffer_1100_ppdd[2072];

    auto g_z_x_0_0_y_x_yy_yy = buffer_1100_ppdd[2073];

    auto g_z_x_0_0_y_x_yy_yz = buffer_1100_ppdd[2074];

    auto g_z_x_0_0_y_x_yy_zz = buffer_1100_ppdd[2075];

    auto g_z_x_0_0_y_x_yz_xx = buffer_1100_ppdd[2076];

    auto g_z_x_0_0_y_x_yz_xy = buffer_1100_ppdd[2077];

    auto g_z_x_0_0_y_x_yz_xz = buffer_1100_ppdd[2078];

    auto g_z_x_0_0_y_x_yz_yy = buffer_1100_ppdd[2079];

    auto g_z_x_0_0_y_x_yz_yz = buffer_1100_ppdd[2080];

    auto g_z_x_0_0_y_x_yz_zz = buffer_1100_ppdd[2081];

    auto g_z_x_0_0_y_x_zz_xx = buffer_1100_ppdd[2082];

    auto g_z_x_0_0_y_x_zz_xy = buffer_1100_ppdd[2083];

    auto g_z_x_0_0_y_x_zz_xz = buffer_1100_ppdd[2084];

    auto g_z_x_0_0_y_x_zz_yy = buffer_1100_ppdd[2085];

    auto g_z_x_0_0_y_x_zz_yz = buffer_1100_ppdd[2086];

    auto g_z_x_0_0_y_x_zz_zz = buffer_1100_ppdd[2087];

    auto g_z_x_0_0_y_y_xx_xx = buffer_1100_ppdd[2088];

    auto g_z_x_0_0_y_y_xx_xy = buffer_1100_ppdd[2089];

    auto g_z_x_0_0_y_y_xx_xz = buffer_1100_ppdd[2090];

    auto g_z_x_0_0_y_y_xx_yy = buffer_1100_ppdd[2091];

    auto g_z_x_0_0_y_y_xx_yz = buffer_1100_ppdd[2092];

    auto g_z_x_0_0_y_y_xx_zz = buffer_1100_ppdd[2093];

    auto g_z_x_0_0_y_y_xy_xx = buffer_1100_ppdd[2094];

    auto g_z_x_0_0_y_y_xy_xy = buffer_1100_ppdd[2095];

    auto g_z_x_0_0_y_y_xy_xz = buffer_1100_ppdd[2096];

    auto g_z_x_0_0_y_y_xy_yy = buffer_1100_ppdd[2097];

    auto g_z_x_0_0_y_y_xy_yz = buffer_1100_ppdd[2098];

    auto g_z_x_0_0_y_y_xy_zz = buffer_1100_ppdd[2099];

    auto g_z_x_0_0_y_y_xz_xx = buffer_1100_ppdd[2100];

    auto g_z_x_0_0_y_y_xz_xy = buffer_1100_ppdd[2101];

    auto g_z_x_0_0_y_y_xz_xz = buffer_1100_ppdd[2102];

    auto g_z_x_0_0_y_y_xz_yy = buffer_1100_ppdd[2103];

    auto g_z_x_0_0_y_y_xz_yz = buffer_1100_ppdd[2104];

    auto g_z_x_0_0_y_y_xz_zz = buffer_1100_ppdd[2105];

    auto g_z_x_0_0_y_y_yy_xx = buffer_1100_ppdd[2106];

    auto g_z_x_0_0_y_y_yy_xy = buffer_1100_ppdd[2107];

    auto g_z_x_0_0_y_y_yy_xz = buffer_1100_ppdd[2108];

    auto g_z_x_0_0_y_y_yy_yy = buffer_1100_ppdd[2109];

    auto g_z_x_0_0_y_y_yy_yz = buffer_1100_ppdd[2110];

    auto g_z_x_0_0_y_y_yy_zz = buffer_1100_ppdd[2111];

    auto g_z_x_0_0_y_y_yz_xx = buffer_1100_ppdd[2112];

    auto g_z_x_0_0_y_y_yz_xy = buffer_1100_ppdd[2113];

    auto g_z_x_0_0_y_y_yz_xz = buffer_1100_ppdd[2114];

    auto g_z_x_0_0_y_y_yz_yy = buffer_1100_ppdd[2115];

    auto g_z_x_0_0_y_y_yz_yz = buffer_1100_ppdd[2116];

    auto g_z_x_0_0_y_y_yz_zz = buffer_1100_ppdd[2117];

    auto g_z_x_0_0_y_y_zz_xx = buffer_1100_ppdd[2118];

    auto g_z_x_0_0_y_y_zz_xy = buffer_1100_ppdd[2119];

    auto g_z_x_0_0_y_y_zz_xz = buffer_1100_ppdd[2120];

    auto g_z_x_0_0_y_y_zz_yy = buffer_1100_ppdd[2121];

    auto g_z_x_0_0_y_y_zz_yz = buffer_1100_ppdd[2122];

    auto g_z_x_0_0_y_y_zz_zz = buffer_1100_ppdd[2123];

    auto g_z_x_0_0_y_z_xx_xx = buffer_1100_ppdd[2124];

    auto g_z_x_0_0_y_z_xx_xy = buffer_1100_ppdd[2125];

    auto g_z_x_0_0_y_z_xx_xz = buffer_1100_ppdd[2126];

    auto g_z_x_0_0_y_z_xx_yy = buffer_1100_ppdd[2127];

    auto g_z_x_0_0_y_z_xx_yz = buffer_1100_ppdd[2128];

    auto g_z_x_0_0_y_z_xx_zz = buffer_1100_ppdd[2129];

    auto g_z_x_0_0_y_z_xy_xx = buffer_1100_ppdd[2130];

    auto g_z_x_0_0_y_z_xy_xy = buffer_1100_ppdd[2131];

    auto g_z_x_0_0_y_z_xy_xz = buffer_1100_ppdd[2132];

    auto g_z_x_0_0_y_z_xy_yy = buffer_1100_ppdd[2133];

    auto g_z_x_0_0_y_z_xy_yz = buffer_1100_ppdd[2134];

    auto g_z_x_0_0_y_z_xy_zz = buffer_1100_ppdd[2135];

    auto g_z_x_0_0_y_z_xz_xx = buffer_1100_ppdd[2136];

    auto g_z_x_0_0_y_z_xz_xy = buffer_1100_ppdd[2137];

    auto g_z_x_0_0_y_z_xz_xz = buffer_1100_ppdd[2138];

    auto g_z_x_0_0_y_z_xz_yy = buffer_1100_ppdd[2139];

    auto g_z_x_0_0_y_z_xz_yz = buffer_1100_ppdd[2140];

    auto g_z_x_0_0_y_z_xz_zz = buffer_1100_ppdd[2141];

    auto g_z_x_0_0_y_z_yy_xx = buffer_1100_ppdd[2142];

    auto g_z_x_0_0_y_z_yy_xy = buffer_1100_ppdd[2143];

    auto g_z_x_0_0_y_z_yy_xz = buffer_1100_ppdd[2144];

    auto g_z_x_0_0_y_z_yy_yy = buffer_1100_ppdd[2145];

    auto g_z_x_0_0_y_z_yy_yz = buffer_1100_ppdd[2146];

    auto g_z_x_0_0_y_z_yy_zz = buffer_1100_ppdd[2147];

    auto g_z_x_0_0_y_z_yz_xx = buffer_1100_ppdd[2148];

    auto g_z_x_0_0_y_z_yz_xy = buffer_1100_ppdd[2149];

    auto g_z_x_0_0_y_z_yz_xz = buffer_1100_ppdd[2150];

    auto g_z_x_0_0_y_z_yz_yy = buffer_1100_ppdd[2151];

    auto g_z_x_0_0_y_z_yz_yz = buffer_1100_ppdd[2152];

    auto g_z_x_0_0_y_z_yz_zz = buffer_1100_ppdd[2153];

    auto g_z_x_0_0_y_z_zz_xx = buffer_1100_ppdd[2154];

    auto g_z_x_0_0_y_z_zz_xy = buffer_1100_ppdd[2155];

    auto g_z_x_0_0_y_z_zz_xz = buffer_1100_ppdd[2156];

    auto g_z_x_0_0_y_z_zz_yy = buffer_1100_ppdd[2157];

    auto g_z_x_0_0_y_z_zz_yz = buffer_1100_ppdd[2158];

    auto g_z_x_0_0_y_z_zz_zz = buffer_1100_ppdd[2159];

    auto g_z_x_0_0_z_x_xx_xx = buffer_1100_ppdd[2160];

    auto g_z_x_0_0_z_x_xx_xy = buffer_1100_ppdd[2161];

    auto g_z_x_0_0_z_x_xx_xz = buffer_1100_ppdd[2162];

    auto g_z_x_0_0_z_x_xx_yy = buffer_1100_ppdd[2163];

    auto g_z_x_0_0_z_x_xx_yz = buffer_1100_ppdd[2164];

    auto g_z_x_0_0_z_x_xx_zz = buffer_1100_ppdd[2165];

    auto g_z_x_0_0_z_x_xy_xx = buffer_1100_ppdd[2166];

    auto g_z_x_0_0_z_x_xy_xy = buffer_1100_ppdd[2167];

    auto g_z_x_0_0_z_x_xy_xz = buffer_1100_ppdd[2168];

    auto g_z_x_0_0_z_x_xy_yy = buffer_1100_ppdd[2169];

    auto g_z_x_0_0_z_x_xy_yz = buffer_1100_ppdd[2170];

    auto g_z_x_0_0_z_x_xy_zz = buffer_1100_ppdd[2171];

    auto g_z_x_0_0_z_x_xz_xx = buffer_1100_ppdd[2172];

    auto g_z_x_0_0_z_x_xz_xy = buffer_1100_ppdd[2173];

    auto g_z_x_0_0_z_x_xz_xz = buffer_1100_ppdd[2174];

    auto g_z_x_0_0_z_x_xz_yy = buffer_1100_ppdd[2175];

    auto g_z_x_0_0_z_x_xz_yz = buffer_1100_ppdd[2176];

    auto g_z_x_0_0_z_x_xz_zz = buffer_1100_ppdd[2177];

    auto g_z_x_0_0_z_x_yy_xx = buffer_1100_ppdd[2178];

    auto g_z_x_0_0_z_x_yy_xy = buffer_1100_ppdd[2179];

    auto g_z_x_0_0_z_x_yy_xz = buffer_1100_ppdd[2180];

    auto g_z_x_0_0_z_x_yy_yy = buffer_1100_ppdd[2181];

    auto g_z_x_0_0_z_x_yy_yz = buffer_1100_ppdd[2182];

    auto g_z_x_0_0_z_x_yy_zz = buffer_1100_ppdd[2183];

    auto g_z_x_0_0_z_x_yz_xx = buffer_1100_ppdd[2184];

    auto g_z_x_0_0_z_x_yz_xy = buffer_1100_ppdd[2185];

    auto g_z_x_0_0_z_x_yz_xz = buffer_1100_ppdd[2186];

    auto g_z_x_0_0_z_x_yz_yy = buffer_1100_ppdd[2187];

    auto g_z_x_0_0_z_x_yz_yz = buffer_1100_ppdd[2188];

    auto g_z_x_0_0_z_x_yz_zz = buffer_1100_ppdd[2189];

    auto g_z_x_0_0_z_x_zz_xx = buffer_1100_ppdd[2190];

    auto g_z_x_0_0_z_x_zz_xy = buffer_1100_ppdd[2191];

    auto g_z_x_0_0_z_x_zz_xz = buffer_1100_ppdd[2192];

    auto g_z_x_0_0_z_x_zz_yy = buffer_1100_ppdd[2193];

    auto g_z_x_0_0_z_x_zz_yz = buffer_1100_ppdd[2194];

    auto g_z_x_0_0_z_x_zz_zz = buffer_1100_ppdd[2195];

    auto g_z_x_0_0_z_y_xx_xx = buffer_1100_ppdd[2196];

    auto g_z_x_0_0_z_y_xx_xy = buffer_1100_ppdd[2197];

    auto g_z_x_0_0_z_y_xx_xz = buffer_1100_ppdd[2198];

    auto g_z_x_0_0_z_y_xx_yy = buffer_1100_ppdd[2199];

    auto g_z_x_0_0_z_y_xx_yz = buffer_1100_ppdd[2200];

    auto g_z_x_0_0_z_y_xx_zz = buffer_1100_ppdd[2201];

    auto g_z_x_0_0_z_y_xy_xx = buffer_1100_ppdd[2202];

    auto g_z_x_0_0_z_y_xy_xy = buffer_1100_ppdd[2203];

    auto g_z_x_0_0_z_y_xy_xz = buffer_1100_ppdd[2204];

    auto g_z_x_0_0_z_y_xy_yy = buffer_1100_ppdd[2205];

    auto g_z_x_0_0_z_y_xy_yz = buffer_1100_ppdd[2206];

    auto g_z_x_0_0_z_y_xy_zz = buffer_1100_ppdd[2207];

    auto g_z_x_0_0_z_y_xz_xx = buffer_1100_ppdd[2208];

    auto g_z_x_0_0_z_y_xz_xy = buffer_1100_ppdd[2209];

    auto g_z_x_0_0_z_y_xz_xz = buffer_1100_ppdd[2210];

    auto g_z_x_0_0_z_y_xz_yy = buffer_1100_ppdd[2211];

    auto g_z_x_0_0_z_y_xz_yz = buffer_1100_ppdd[2212];

    auto g_z_x_0_0_z_y_xz_zz = buffer_1100_ppdd[2213];

    auto g_z_x_0_0_z_y_yy_xx = buffer_1100_ppdd[2214];

    auto g_z_x_0_0_z_y_yy_xy = buffer_1100_ppdd[2215];

    auto g_z_x_0_0_z_y_yy_xz = buffer_1100_ppdd[2216];

    auto g_z_x_0_0_z_y_yy_yy = buffer_1100_ppdd[2217];

    auto g_z_x_0_0_z_y_yy_yz = buffer_1100_ppdd[2218];

    auto g_z_x_0_0_z_y_yy_zz = buffer_1100_ppdd[2219];

    auto g_z_x_0_0_z_y_yz_xx = buffer_1100_ppdd[2220];

    auto g_z_x_0_0_z_y_yz_xy = buffer_1100_ppdd[2221];

    auto g_z_x_0_0_z_y_yz_xz = buffer_1100_ppdd[2222];

    auto g_z_x_0_0_z_y_yz_yy = buffer_1100_ppdd[2223];

    auto g_z_x_0_0_z_y_yz_yz = buffer_1100_ppdd[2224];

    auto g_z_x_0_0_z_y_yz_zz = buffer_1100_ppdd[2225];

    auto g_z_x_0_0_z_y_zz_xx = buffer_1100_ppdd[2226];

    auto g_z_x_0_0_z_y_zz_xy = buffer_1100_ppdd[2227];

    auto g_z_x_0_0_z_y_zz_xz = buffer_1100_ppdd[2228];

    auto g_z_x_0_0_z_y_zz_yy = buffer_1100_ppdd[2229];

    auto g_z_x_0_0_z_y_zz_yz = buffer_1100_ppdd[2230];

    auto g_z_x_0_0_z_y_zz_zz = buffer_1100_ppdd[2231];

    auto g_z_x_0_0_z_z_xx_xx = buffer_1100_ppdd[2232];

    auto g_z_x_0_0_z_z_xx_xy = buffer_1100_ppdd[2233];

    auto g_z_x_0_0_z_z_xx_xz = buffer_1100_ppdd[2234];

    auto g_z_x_0_0_z_z_xx_yy = buffer_1100_ppdd[2235];

    auto g_z_x_0_0_z_z_xx_yz = buffer_1100_ppdd[2236];

    auto g_z_x_0_0_z_z_xx_zz = buffer_1100_ppdd[2237];

    auto g_z_x_0_0_z_z_xy_xx = buffer_1100_ppdd[2238];

    auto g_z_x_0_0_z_z_xy_xy = buffer_1100_ppdd[2239];

    auto g_z_x_0_0_z_z_xy_xz = buffer_1100_ppdd[2240];

    auto g_z_x_0_0_z_z_xy_yy = buffer_1100_ppdd[2241];

    auto g_z_x_0_0_z_z_xy_yz = buffer_1100_ppdd[2242];

    auto g_z_x_0_0_z_z_xy_zz = buffer_1100_ppdd[2243];

    auto g_z_x_0_0_z_z_xz_xx = buffer_1100_ppdd[2244];

    auto g_z_x_0_0_z_z_xz_xy = buffer_1100_ppdd[2245];

    auto g_z_x_0_0_z_z_xz_xz = buffer_1100_ppdd[2246];

    auto g_z_x_0_0_z_z_xz_yy = buffer_1100_ppdd[2247];

    auto g_z_x_0_0_z_z_xz_yz = buffer_1100_ppdd[2248];

    auto g_z_x_0_0_z_z_xz_zz = buffer_1100_ppdd[2249];

    auto g_z_x_0_0_z_z_yy_xx = buffer_1100_ppdd[2250];

    auto g_z_x_0_0_z_z_yy_xy = buffer_1100_ppdd[2251];

    auto g_z_x_0_0_z_z_yy_xz = buffer_1100_ppdd[2252];

    auto g_z_x_0_0_z_z_yy_yy = buffer_1100_ppdd[2253];

    auto g_z_x_0_0_z_z_yy_yz = buffer_1100_ppdd[2254];

    auto g_z_x_0_0_z_z_yy_zz = buffer_1100_ppdd[2255];

    auto g_z_x_0_0_z_z_yz_xx = buffer_1100_ppdd[2256];

    auto g_z_x_0_0_z_z_yz_xy = buffer_1100_ppdd[2257];

    auto g_z_x_0_0_z_z_yz_xz = buffer_1100_ppdd[2258];

    auto g_z_x_0_0_z_z_yz_yy = buffer_1100_ppdd[2259];

    auto g_z_x_0_0_z_z_yz_yz = buffer_1100_ppdd[2260];

    auto g_z_x_0_0_z_z_yz_zz = buffer_1100_ppdd[2261];

    auto g_z_x_0_0_z_z_zz_xx = buffer_1100_ppdd[2262];

    auto g_z_x_0_0_z_z_zz_xy = buffer_1100_ppdd[2263];

    auto g_z_x_0_0_z_z_zz_xz = buffer_1100_ppdd[2264];

    auto g_z_x_0_0_z_z_zz_yy = buffer_1100_ppdd[2265];

    auto g_z_x_0_0_z_z_zz_yz = buffer_1100_ppdd[2266];

    auto g_z_x_0_0_z_z_zz_zz = buffer_1100_ppdd[2267];

    auto g_z_y_0_0_x_x_xx_xx = buffer_1100_ppdd[2268];

    auto g_z_y_0_0_x_x_xx_xy = buffer_1100_ppdd[2269];

    auto g_z_y_0_0_x_x_xx_xz = buffer_1100_ppdd[2270];

    auto g_z_y_0_0_x_x_xx_yy = buffer_1100_ppdd[2271];

    auto g_z_y_0_0_x_x_xx_yz = buffer_1100_ppdd[2272];

    auto g_z_y_0_0_x_x_xx_zz = buffer_1100_ppdd[2273];

    auto g_z_y_0_0_x_x_xy_xx = buffer_1100_ppdd[2274];

    auto g_z_y_0_0_x_x_xy_xy = buffer_1100_ppdd[2275];

    auto g_z_y_0_0_x_x_xy_xz = buffer_1100_ppdd[2276];

    auto g_z_y_0_0_x_x_xy_yy = buffer_1100_ppdd[2277];

    auto g_z_y_0_0_x_x_xy_yz = buffer_1100_ppdd[2278];

    auto g_z_y_0_0_x_x_xy_zz = buffer_1100_ppdd[2279];

    auto g_z_y_0_0_x_x_xz_xx = buffer_1100_ppdd[2280];

    auto g_z_y_0_0_x_x_xz_xy = buffer_1100_ppdd[2281];

    auto g_z_y_0_0_x_x_xz_xz = buffer_1100_ppdd[2282];

    auto g_z_y_0_0_x_x_xz_yy = buffer_1100_ppdd[2283];

    auto g_z_y_0_0_x_x_xz_yz = buffer_1100_ppdd[2284];

    auto g_z_y_0_0_x_x_xz_zz = buffer_1100_ppdd[2285];

    auto g_z_y_0_0_x_x_yy_xx = buffer_1100_ppdd[2286];

    auto g_z_y_0_0_x_x_yy_xy = buffer_1100_ppdd[2287];

    auto g_z_y_0_0_x_x_yy_xz = buffer_1100_ppdd[2288];

    auto g_z_y_0_0_x_x_yy_yy = buffer_1100_ppdd[2289];

    auto g_z_y_0_0_x_x_yy_yz = buffer_1100_ppdd[2290];

    auto g_z_y_0_0_x_x_yy_zz = buffer_1100_ppdd[2291];

    auto g_z_y_0_0_x_x_yz_xx = buffer_1100_ppdd[2292];

    auto g_z_y_0_0_x_x_yz_xy = buffer_1100_ppdd[2293];

    auto g_z_y_0_0_x_x_yz_xz = buffer_1100_ppdd[2294];

    auto g_z_y_0_0_x_x_yz_yy = buffer_1100_ppdd[2295];

    auto g_z_y_0_0_x_x_yz_yz = buffer_1100_ppdd[2296];

    auto g_z_y_0_0_x_x_yz_zz = buffer_1100_ppdd[2297];

    auto g_z_y_0_0_x_x_zz_xx = buffer_1100_ppdd[2298];

    auto g_z_y_0_0_x_x_zz_xy = buffer_1100_ppdd[2299];

    auto g_z_y_0_0_x_x_zz_xz = buffer_1100_ppdd[2300];

    auto g_z_y_0_0_x_x_zz_yy = buffer_1100_ppdd[2301];

    auto g_z_y_0_0_x_x_zz_yz = buffer_1100_ppdd[2302];

    auto g_z_y_0_0_x_x_zz_zz = buffer_1100_ppdd[2303];

    auto g_z_y_0_0_x_y_xx_xx = buffer_1100_ppdd[2304];

    auto g_z_y_0_0_x_y_xx_xy = buffer_1100_ppdd[2305];

    auto g_z_y_0_0_x_y_xx_xz = buffer_1100_ppdd[2306];

    auto g_z_y_0_0_x_y_xx_yy = buffer_1100_ppdd[2307];

    auto g_z_y_0_0_x_y_xx_yz = buffer_1100_ppdd[2308];

    auto g_z_y_0_0_x_y_xx_zz = buffer_1100_ppdd[2309];

    auto g_z_y_0_0_x_y_xy_xx = buffer_1100_ppdd[2310];

    auto g_z_y_0_0_x_y_xy_xy = buffer_1100_ppdd[2311];

    auto g_z_y_0_0_x_y_xy_xz = buffer_1100_ppdd[2312];

    auto g_z_y_0_0_x_y_xy_yy = buffer_1100_ppdd[2313];

    auto g_z_y_0_0_x_y_xy_yz = buffer_1100_ppdd[2314];

    auto g_z_y_0_0_x_y_xy_zz = buffer_1100_ppdd[2315];

    auto g_z_y_0_0_x_y_xz_xx = buffer_1100_ppdd[2316];

    auto g_z_y_0_0_x_y_xz_xy = buffer_1100_ppdd[2317];

    auto g_z_y_0_0_x_y_xz_xz = buffer_1100_ppdd[2318];

    auto g_z_y_0_0_x_y_xz_yy = buffer_1100_ppdd[2319];

    auto g_z_y_0_0_x_y_xz_yz = buffer_1100_ppdd[2320];

    auto g_z_y_0_0_x_y_xz_zz = buffer_1100_ppdd[2321];

    auto g_z_y_0_0_x_y_yy_xx = buffer_1100_ppdd[2322];

    auto g_z_y_0_0_x_y_yy_xy = buffer_1100_ppdd[2323];

    auto g_z_y_0_0_x_y_yy_xz = buffer_1100_ppdd[2324];

    auto g_z_y_0_0_x_y_yy_yy = buffer_1100_ppdd[2325];

    auto g_z_y_0_0_x_y_yy_yz = buffer_1100_ppdd[2326];

    auto g_z_y_0_0_x_y_yy_zz = buffer_1100_ppdd[2327];

    auto g_z_y_0_0_x_y_yz_xx = buffer_1100_ppdd[2328];

    auto g_z_y_0_0_x_y_yz_xy = buffer_1100_ppdd[2329];

    auto g_z_y_0_0_x_y_yz_xz = buffer_1100_ppdd[2330];

    auto g_z_y_0_0_x_y_yz_yy = buffer_1100_ppdd[2331];

    auto g_z_y_0_0_x_y_yz_yz = buffer_1100_ppdd[2332];

    auto g_z_y_0_0_x_y_yz_zz = buffer_1100_ppdd[2333];

    auto g_z_y_0_0_x_y_zz_xx = buffer_1100_ppdd[2334];

    auto g_z_y_0_0_x_y_zz_xy = buffer_1100_ppdd[2335];

    auto g_z_y_0_0_x_y_zz_xz = buffer_1100_ppdd[2336];

    auto g_z_y_0_0_x_y_zz_yy = buffer_1100_ppdd[2337];

    auto g_z_y_0_0_x_y_zz_yz = buffer_1100_ppdd[2338];

    auto g_z_y_0_0_x_y_zz_zz = buffer_1100_ppdd[2339];

    auto g_z_y_0_0_x_z_xx_xx = buffer_1100_ppdd[2340];

    auto g_z_y_0_0_x_z_xx_xy = buffer_1100_ppdd[2341];

    auto g_z_y_0_0_x_z_xx_xz = buffer_1100_ppdd[2342];

    auto g_z_y_0_0_x_z_xx_yy = buffer_1100_ppdd[2343];

    auto g_z_y_0_0_x_z_xx_yz = buffer_1100_ppdd[2344];

    auto g_z_y_0_0_x_z_xx_zz = buffer_1100_ppdd[2345];

    auto g_z_y_0_0_x_z_xy_xx = buffer_1100_ppdd[2346];

    auto g_z_y_0_0_x_z_xy_xy = buffer_1100_ppdd[2347];

    auto g_z_y_0_0_x_z_xy_xz = buffer_1100_ppdd[2348];

    auto g_z_y_0_0_x_z_xy_yy = buffer_1100_ppdd[2349];

    auto g_z_y_0_0_x_z_xy_yz = buffer_1100_ppdd[2350];

    auto g_z_y_0_0_x_z_xy_zz = buffer_1100_ppdd[2351];

    auto g_z_y_0_0_x_z_xz_xx = buffer_1100_ppdd[2352];

    auto g_z_y_0_0_x_z_xz_xy = buffer_1100_ppdd[2353];

    auto g_z_y_0_0_x_z_xz_xz = buffer_1100_ppdd[2354];

    auto g_z_y_0_0_x_z_xz_yy = buffer_1100_ppdd[2355];

    auto g_z_y_0_0_x_z_xz_yz = buffer_1100_ppdd[2356];

    auto g_z_y_0_0_x_z_xz_zz = buffer_1100_ppdd[2357];

    auto g_z_y_0_0_x_z_yy_xx = buffer_1100_ppdd[2358];

    auto g_z_y_0_0_x_z_yy_xy = buffer_1100_ppdd[2359];

    auto g_z_y_0_0_x_z_yy_xz = buffer_1100_ppdd[2360];

    auto g_z_y_0_0_x_z_yy_yy = buffer_1100_ppdd[2361];

    auto g_z_y_0_0_x_z_yy_yz = buffer_1100_ppdd[2362];

    auto g_z_y_0_0_x_z_yy_zz = buffer_1100_ppdd[2363];

    auto g_z_y_0_0_x_z_yz_xx = buffer_1100_ppdd[2364];

    auto g_z_y_0_0_x_z_yz_xy = buffer_1100_ppdd[2365];

    auto g_z_y_0_0_x_z_yz_xz = buffer_1100_ppdd[2366];

    auto g_z_y_0_0_x_z_yz_yy = buffer_1100_ppdd[2367];

    auto g_z_y_0_0_x_z_yz_yz = buffer_1100_ppdd[2368];

    auto g_z_y_0_0_x_z_yz_zz = buffer_1100_ppdd[2369];

    auto g_z_y_0_0_x_z_zz_xx = buffer_1100_ppdd[2370];

    auto g_z_y_0_0_x_z_zz_xy = buffer_1100_ppdd[2371];

    auto g_z_y_0_0_x_z_zz_xz = buffer_1100_ppdd[2372];

    auto g_z_y_0_0_x_z_zz_yy = buffer_1100_ppdd[2373];

    auto g_z_y_0_0_x_z_zz_yz = buffer_1100_ppdd[2374];

    auto g_z_y_0_0_x_z_zz_zz = buffer_1100_ppdd[2375];

    auto g_z_y_0_0_y_x_xx_xx = buffer_1100_ppdd[2376];

    auto g_z_y_0_0_y_x_xx_xy = buffer_1100_ppdd[2377];

    auto g_z_y_0_0_y_x_xx_xz = buffer_1100_ppdd[2378];

    auto g_z_y_0_0_y_x_xx_yy = buffer_1100_ppdd[2379];

    auto g_z_y_0_0_y_x_xx_yz = buffer_1100_ppdd[2380];

    auto g_z_y_0_0_y_x_xx_zz = buffer_1100_ppdd[2381];

    auto g_z_y_0_0_y_x_xy_xx = buffer_1100_ppdd[2382];

    auto g_z_y_0_0_y_x_xy_xy = buffer_1100_ppdd[2383];

    auto g_z_y_0_0_y_x_xy_xz = buffer_1100_ppdd[2384];

    auto g_z_y_0_0_y_x_xy_yy = buffer_1100_ppdd[2385];

    auto g_z_y_0_0_y_x_xy_yz = buffer_1100_ppdd[2386];

    auto g_z_y_0_0_y_x_xy_zz = buffer_1100_ppdd[2387];

    auto g_z_y_0_0_y_x_xz_xx = buffer_1100_ppdd[2388];

    auto g_z_y_0_0_y_x_xz_xy = buffer_1100_ppdd[2389];

    auto g_z_y_0_0_y_x_xz_xz = buffer_1100_ppdd[2390];

    auto g_z_y_0_0_y_x_xz_yy = buffer_1100_ppdd[2391];

    auto g_z_y_0_0_y_x_xz_yz = buffer_1100_ppdd[2392];

    auto g_z_y_0_0_y_x_xz_zz = buffer_1100_ppdd[2393];

    auto g_z_y_0_0_y_x_yy_xx = buffer_1100_ppdd[2394];

    auto g_z_y_0_0_y_x_yy_xy = buffer_1100_ppdd[2395];

    auto g_z_y_0_0_y_x_yy_xz = buffer_1100_ppdd[2396];

    auto g_z_y_0_0_y_x_yy_yy = buffer_1100_ppdd[2397];

    auto g_z_y_0_0_y_x_yy_yz = buffer_1100_ppdd[2398];

    auto g_z_y_0_0_y_x_yy_zz = buffer_1100_ppdd[2399];

    auto g_z_y_0_0_y_x_yz_xx = buffer_1100_ppdd[2400];

    auto g_z_y_0_0_y_x_yz_xy = buffer_1100_ppdd[2401];

    auto g_z_y_0_0_y_x_yz_xz = buffer_1100_ppdd[2402];

    auto g_z_y_0_0_y_x_yz_yy = buffer_1100_ppdd[2403];

    auto g_z_y_0_0_y_x_yz_yz = buffer_1100_ppdd[2404];

    auto g_z_y_0_0_y_x_yz_zz = buffer_1100_ppdd[2405];

    auto g_z_y_0_0_y_x_zz_xx = buffer_1100_ppdd[2406];

    auto g_z_y_0_0_y_x_zz_xy = buffer_1100_ppdd[2407];

    auto g_z_y_0_0_y_x_zz_xz = buffer_1100_ppdd[2408];

    auto g_z_y_0_0_y_x_zz_yy = buffer_1100_ppdd[2409];

    auto g_z_y_0_0_y_x_zz_yz = buffer_1100_ppdd[2410];

    auto g_z_y_0_0_y_x_zz_zz = buffer_1100_ppdd[2411];

    auto g_z_y_0_0_y_y_xx_xx = buffer_1100_ppdd[2412];

    auto g_z_y_0_0_y_y_xx_xy = buffer_1100_ppdd[2413];

    auto g_z_y_0_0_y_y_xx_xz = buffer_1100_ppdd[2414];

    auto g_z_y_0_0_y_y_xx_yy = buffer_1100_ppdd[2415];

    auto g_z_y_0_0_y_y_xx_yz = buffer_1100_ppdd[2416];

    auto g_z_y_0_0_y_y_xx_zz = buffer_1100_ppdd[2417];

    auto g_z_y_0_0_y_y_xy_xx = buffer_1100_ppdd[2418];

    auto g_z_y_0_0_y_y_xy_xy = buffer_1100_ppdd[2419];

    auto g_z_y_0_0_y_y_xy_xz = buffer_1100_ppdd[2420];

    auto g_z_y_0_0_y_y_xy_yy = buffer_1100_ppdd[2421];

    auto g_z_y_0_0_y_y_xy_yz = buffer_1100_ppdd[2422];

    auto g_z_y_0_0_y_y_xy_zz = buffer_1100_ppdd[2423];

    auto g_z_y_0_0_y_y_xz_xx = buffer_1100_ppdd[2424];

    auto g_z_y_0_0_y_y_xz_xy = buffer_1100_ppdd[2425];

    auto g_z_y_0_0_y_y_xz_xz = buffer_1100_ppdd[2426];

    auto g_z_y_0_0_y_y_xz_yy = buffer_1100_ppdd[2427];

    auto g_z_y_0_0_y_y_xz_yz = buffer_1100_ppdd[2428];

    auto g_z_y_0_0_y_y_xz_zz = buffer_1100_ppdd[2429];

    auto g_z_y_0_0_y_y_yy_xx = buffer_1100_ppdd[2430];

    auto g_z_y_0_0_y_y_yy_xy = buffer_1100_ppdd[2431];

    auto g_z_y_0_0_y_y_yy_xz = buffer_1100_ppdd[2432];

    auto g_z_y_0_0_y_y_yy_yy = buffer_1100_ppdd[2433];

    auto g_z_y_0_0_y_y_yy_yz = buffer_1100_ppdd[2434];

    auto g_z_y_0_0_y_y_yy_zz = buffer_1100_ppdd[2435];

    auto g_z_y_0_0_y_y_yz_xx = buffer_1100_ppdd[2436];

    auto g_z_y_0_0_y_y_yz_xy = buffer_1100_ppdd[2437];

    auto g_z_y_0_0_y_y_yz_xz = buffer_1100_ppdd[2438];

    auto g_z_y_0_0_y_y_yz_yy = buffer_1100_ppdd[2439];

    auto g_z_y_0_0_y_y_yz_yz = buffer_1100_ppdd[2440];

    auto g_z_y_0_0_y_y_yz_zz = buffer_1100_ppdd[2441];

    auto g_z_y_0_0_y_y_zz_xx = buffer_1100_ppdd[2442];

    auto g_z_y_0_0_y_y_zz_xy = buffer_1100_ppdd[2443];

    auto g_z_y_0_0_y_y_zz_xz = buffer_1100_ppdd[2444];

    auto g_z_y_0_0_y_y_zz_yy = buffer_1100_ppdd[2445];

    auto g_z_y_0_0_y_y_zz_yz = buffer_1100_ppdd[2446];

    auto g_z_y_0_0_y_y_zz_zz = buffer_1100_ppdd[2447];

    auto g_z_y_0_0_y_z_xx_xx = buffer_1100_ppdd[2448];

    auto g_z_y_0_0_y_z_xx_xy = buffer_1100_ppdd[2449];

    auto g_z_y_0_0_y_z_xx_xz = buffer_1100_ppdd[2450];

    auto g_z_y_0_0_y_z_xx_yy = buffer_1100_ppdd[2451];

    auto g_z_y_0_0_y_z_xx_yz = buffer_1100_ppdd[2452];

    auto g_z_y_0_0_y_z_xx_zz = buffer_1100_ppdd[2453];

    auto g_z_y_0_0_y_z_xy_xx = buffer_1100_ppdd[2454];

    auto g_z_y_0_0_y_z_xy_xy = buffer_1100_ppdd[2455];

    auto g_z_y_0_0_y_z_xy_xz = buffer_1100_ppdd[2456];

    auto g_z_y_0_0_y_z_xy_yy = buffer_1100_ppdd[2457];

    auto g_z_y_0_0_y_z_xy_yz = buffer_1100_ppdd[2458];

    auto g_z_y_0_0_y_z_xy_zz = buffer_1100_ppdd[2459];

    auto g_z_y_0_0_y_z_xz_xx = buffer_1100_ppdd[2460];

    auto g_z_y_0_0_y_z_xz_xy = buffer_1100_ppdd[2461];

    auto g_z_y_0_0_y_z_xz_xz = buffer_1100_ppdd[2462];

    auto g_z_y_0_0_y_z_xz_yy = buffer_1100_ppdd[2463];

    auto g_z_y_0_0_y_z_xz_yz = buffer_1100_ppdd[2464];

    auto g_z_y_0_0_y_z_xz_zz = buffer_1100_ppdd[2465];

    auto g_z_y_0_0_y_z_yy_xx = buffer_1100_ppdd[2466];

    auto g_z_y_0_0_y_z_yy_xy = buffer_1100_ppdd[2467];

    auto g_z_y_0_0_y_z_yy_xz = buffer_1100_ppdd[2468];

    auto g_z_y_0_0_y_z_yy_yy = buffer_1100_ppdd[2469];

    auto g_z_y_0_0_y_z_yy_yz = buffer_1100_ppdd[2470];

    auto g_z_y_0_0_y_z_yy_zz = buffer_1100_ppdd[2471];

    auto g_z_y_0_0_y_z_yz_xx = buffer_1100_ppdd[2472];

    auto g_z_y_0_0_y_z_yz_xy = buffer_1100_ppdd[2473];

    auto g_z_y_0_0_y_z_yz_xz = buffer_1100_ppdd[2474];

    auto g_z_y_0_0_y_z_yz_yy = buffer_1100_ppdd[2475];

    auto g_z_y_0_0_y_z_yz_yz = buffer_1100_ppdd[2476];

    auto g_z_y_0_0_y_z_yz_zz = buffer_1100_ppdd[2477];

    auto g_z_y_0_0_y_z_zz_xx = buffer_1100_ppdd[2478];

    auto g_z_y_0_0_y_z_zz_xy = buffer_1100_ppdd[2479];

    auto g_z_y_0_0_y_z_zz_xz = buffer_1100_ppdd[2480];

    auto g_z_y_0_0_y_z_zz_yy = buffer_1100_ppdd[2481];

    auto g_z_y_0_0_y_z_zz_yz = buffer_1100_ppdd[2482];

    auto g_z_y_0_0_y_z_zz_zz = buffer_1100_ppdd[2483];

    auto g_z_y_0_0_z_x_xx_xx = buffer_1100_ppdd[2484];

    auto g_z_y_0_0_z_x_xx_xy = buffer_1100_ppdd[2485];

    auto g_z_y_0_0_z_x_xx_xz = buffer_1100_ppdd[2486];

    auto g_z_y_0_0_z_x_xx_yy = buffer_1100_ppdd[2487];

    auto g_z_y_0_0_z_x_xx_yz = buffer_1100_ppdd[2488];

    auto g_z_y_0_0_z_x_xx_zz = buffer_1100_ppdd[2489];

    auto g_z_y_0_0_z_x_xy_xx = buffer_1100_ppdd[2490];

    auto g_z_y_0_0_z_x_xy_xy = buffer_1100_ppdd[2491];

    auto g_z_y_0_0_z_x_xy_xz = buffer_1100_ppdd[2492];

    auto g_z_y_0_0_z_x_xy_yy = buffer_1100_ppdd[2493];

    auto g_z_y_0_0_z_x_xy_yz = buffer_1100_ppdd[2494];

    auto g_z_y_0_0_z_x_xy_zz = buffer_1100_ppdd[2495];

    auto g_z_y_0_0_z_x_xz_xx = buffer_1100_ppdd[2496];

    auto g_z_y_0_0_z_x_xz_xy = buffer_1100_ppdd[2497];

    auto g_z_y_0_0_z_x_xz_xz = buffer_1100_ppdd[2498];

    auto g_z_y_0_0_z_x_xz_yy = buffer_1100_ppdd[2499];

    auto g_z_y_0_0_z_x_xz_yz = buffer_1100_ppdd[2500];

    auto g_z_y_0_0_z_x_xz_zz = buffer_1100_ppdd[2501];

    auto g_z_y_0_0_z_x_yy_xx = buffer_1100_ppdd[2502];

    auto g_z_y_0_0_z_x_yy_xy = buffer_1100_ppdd[2503];

    auto g_z_y_0_0_z_x_yy_xz = buffer_1100_ppdd[2504];

    auto g_z_y_0_0_z_x_yy_yy = buffer_1100_ppdd[2505];

    auto g_z_y_0_0_z_x_yy_yz = buffer_1100_ppdd[2506];

    auto g_z_y_0_0_z_x_yy_zz = buffer_1100_ppdd[2507];

    auto g_z_y_0_0_z_x_yz_xx = buffer_1100_ppdd[2508];

    auto g_z_y_0_0_z_x_yz_xy = buffer_1100_ppdd[2509];

    auto g_z_y_0_0_z_x_yz_xz = buffer_1100_ppdd[2510];

    auto g_z_y_0_0_z_x_yz_yy = buffer_1100_ppdd[2511];

    auto g_z_y_0_0_z_x_yz_yz = buffer_1100_ppdd[2512];

    auto g_z_y_0_0_z_x_yz_zz = buffer_1100_ppdd[2513];

    auto g_z_y_0_0_z_x_zz_xx = buffer_1100_ppdd[2514];

    auto g_z_y_0_0_z_x_zz_xy = buffer_1100_ppdd[2515];

    auto g_z_y_0_0_z_x_zz_xz = buffer_1100_ppdd[2516];

    auto g_z_y_0_0_z_x_zz_yy = buffer_1100_ppdd[2517];

    auto g_z_y_0_0_z_x_zz_yz = buffer_1100_ppdd[2518];

    auto g_z_y_0_0_z_x_zz_zz = buffer_1100_ppdd[2519];

    auto g_z_y_0_0_z_y_xx_xx = buffer_1100_ppdd[2520];

    auto g_z_y_0_0_z_y_xx_xy = buffer_1100_ppdd[2521];

    auto g_z_y_0_0_z_y_xx_xz = buffer_1100_ppdd[2522];

    auto g_z_y_0_0_z_y_xx_yy = buffer_1100_ppdd[2523];

    auto g_z_y_0_0_z_y_xx_yz = buffer_1100_ppdd[2524];

    auto g_z_y_0_0_z_y_xx_zz = buffer_1100_ppdd[2525];

    auto g_z_y_0_0_z_y_xy_xx = buffer_1100_ppdd[2526];

    auto g_z_y_0_0_z_y_xy_xy = buffer_1100_ppdd[2527];

    auto g_z_y_0_0_z_y_xy_xz = buffer_1100_ppdd[2528];

    auto g_z_y_0_0_z_y_xy_yy = buffer_1100_ppdd[2529];

    auto g_z_y_0_0_z_y_xy_yz = buffer_1100_ppdd[2530];

    auto g_z_y_0_0_z_y_xy_zz = buffer_1100_ppdd[2531];

    auto g_z_y_0_0_z_y_xz_xx = buffer_1100_ppdd[2532];

    auto g_z_y_0_0_z_y_xz_xy = buffer_1100_ppdd[2533];

    auto g_z_y_0_0_z_y_xz_xz = buffer_1100_ppdd[2534];

    auto g_z_y_0_0_z_y_xz_yy = buffer_1100_ppdd[2535];

    auto g_z_y_0_0_z_y_xz_yz = buffer_1100_ppdd[2536];

    auto g_z_y_0_0_z_y_xz_zz = buffer_1100_ppdd[2537];

    auto g_z_y_0_0_z_y_yy_xx = buffer_1100_ppdd[2538];

    auto g_z_y_0_0_z_y_yy_xy = buffer_1100_ppdd[2539];

    auto g_z_y_0_0_z_y_yy_xz = buffer_1100_ppdd[2540];

    auto g_z_y_0_0_z_y_yy_yy = buffer_1100_ppdd[2541];

    auto g_z_y_0_0_z_y_yy_yz = buffer_1100_ppdd[2542];

    auto g_z_y_0_0_z_y_yy_zz = buffer_1100_ppdd[2543];

    auto g_z_y_0_0_z_y_yz_xx = buffer_1100_ppdd[2544];

    auto g_z_y_0_0_z_y_yz_xy = buffer_1100_ppdd[2545];

    auto g_z_y_0_0_z_y_yz_xz = buffer_1100_ppdd[2546];

    auto g_z_y_0_0_z_y_yz_yy = buffer_1100_ppdd[2547];

    auto g_z_y_0_0_z_y_yz_yz = buffer_1100_ppdd[2548];

    auto g_z_y_0_0_z_y_yz_zz = buffer_1100_ppdd[2549];

    auto g_z_y_0_0_z_y_zz_xx = buffer_1100_ppdd[2550];

    auto g_z_y_0_0_z_y_zz_xy = buffer_1100_ppdd[2551];

    auto g_z_y_0_0_z_y_zz_xz = buffer_1100_ppdd[2552];

    auto g_z_y_0_0_z_y_zz_yy = buffer_1100_ppdd[2553];

    auto g_z_y_0_0_z_y_zz_yz = buffer_1100_ppdd[2554];

    auto g_z_y_0_0_z_y_zz_zz = buffer_1100_ppdd[2555];

    auto g_z_y_0_0_z_z_xx_xx = buffer_1100_ppdd[2556];

    auto g_z_y_0_0_z_z_xx_xy = buffer_1100_ppdd[2557];

    auto g_z_y_0_0_z_z_xx_xz = buffer_1100_ppdd[2558];

    auto g_z_y_0_0_z_z_xx_yy = buffer_1100_ppdd[2559];

    auto g_z_y_0_0_z_z_xx_yz = buffer_1100_ppdd[2560];

    auto g_z_y_0_0_z_z_xx_zz = buffer_1100_ppdd[2561];

    auto g_z_y_0_0_z_z_xy_xx = buffer_1100_ppdd[2562];

    auto g_z_y_0_0_z_z_xy_xy = buffer_1100_ppdd[2563];

    auto g_z_y_0_0_z_z_xy_xz = buffer_1100_ppdd[2564];

    auto g_z_y_0_0_z_z_xy_yy = buffer_1100_ppdd[2565];

    auto g_z_y_0_0_z_z_xy_yz = buffer_1100_ppdd[2566];

    auto g_z_y_0_0_z_z_xy_zz = buffer_1100_ppdd[2567];

    auto g_z_y_0_0_z_z_xz_xx = buffer_1100_ppdd[2568];

    auto g_z_y_0_0_z_z_xz_xy = buffer_1100_ppdd[2569];

    auto g_z_y_0_0_z_z_xz_xz = buffer_1100_ppdd[2570];

    auto g_z_y_0_0_z_z_xz_yy = buffer_1100_ppdd[2571];

    auto g_z_y_0_0_z_z_xz_yz = buffer_1100_ppdd[2572];

    auto g_z_y_0_0_z_z_xz_zz = buffer_1100_ppdd[2573];

    auto g_z_y_0_0_z_z_yy_xx = buffer_1100_ppdd[2574];

    auto g_z_y_0_0_z_z_yy_xy = buffer_1100_ppdd[2575];

    auto g_z_y_0_0_z_z_yy_xz = buffer_1100_ppdd[2576];

    auto g_z_y_0_0_z_z_yy_yy = buffer_1100_ppdd[2577];

    auto g_z_y_0_0_z_z_yy_yz = buffer_1100_ppdd[2578];

    auto g_z_y_0_0_z_z_yy_zz = buffer_1100_ppdd[2579];

    auto g_z_y_0_0_z_z_yz_xx = buffer_1100_ppdd[2580];

    auto g_z_y_0_0_z_z_yz_xy = buffer_1100_ppdd[2581];

    auto g_z_y_0_0_z_z_yz_xz = buffer_1100_ppdd[2582];

    auto g_z_y_0_0_z_z_yz_yy = buffer_1100_ppdd[2583];

    auto g_z_y_0_0_z_z_yz_yz = buffer_1100_ppdd[2584];

    auto g_z_y_0_0_z_z_yz_zz = buffer_1100_ppdd[2585];

    auto g_z_y_0_0_z_z_zz_xx = buffer_1100_ppdd[2586];

    auto g_z_y_0_0_z_z_zz_xy = buffer_1100_ppdd[2587];

    auto g_z_y_0_0_z_z_zz_xz = buffer_1100_ppdd[2588];

    auto g_z_y_0_0_z_z_zz_yy = buffer_1100_ppdd[2589];

    auto g_z_y_0_0_z_z_zz_yz = buffer_1100_ppdd[2590];

    auto g_z_y_0_0_z_z_zz_zz = buffer_1100_ppdd[2591];

    auto g_z_z_0_0_x_x_xx_xx = buffer_1100_ppdd[2592];

    auto g_z_z_0_0_x_x_xx_xy = buffer_1100_ppdd[2593];

    auto g_z_z_0_0_x_x_xx_xz = buffer_1100_ppdd[2594];

    auto g_z_z_0_0_x_x_xx_yy = buffer_1100_ppdd[2595];

    auto g_z_z_0_0_x_x_xx_yz = buffer_1100_ppdd[2596];

    auto g_z_z_0_0_x_x_xx_zz = buffer_1100_ppdd[2597];

    auto g_z_z_0_0_x_x_xy_xx = buffer_1100_ppdd[2598];

    auto g_z_z_0_0_x_x_xy_xy = buffer_1100_ppdd[2599];

    auto g_z_z_0_0_x_x_xy_xz = buffer_1100_ppdd[2600];

    auto g_z_z_0_0_x_x_xy_yy = buffer_1100_ppdd[2601];

    auto g_z_z_0_0_x_x_xy_yz = buffer_1100_ppdd[2602];

    auto g_z_z_0_0_x_x_xy_zz = buffer_1100_ppdd[2603];

    auto g_z_z_0_0_x_x_xz_xx = buffer_1100_ppdd[2604];

    auto g_z_z_0_0_x_x_xz_xy = buffer_1100_ppdd[2605];

    auto g_z_z_0_0_x_x_xz_xz = buffer_1100_ppdd[2606];

    auto g_z_z_0_0_x_x_xz_yy = buffer_1100_ppdd[2607];

    auto g_z_z_0_0_x_x_xz_yz = buffer_1100_ppdd[2608];

    auto g_z_z_0_0_x_x_xz_zz = buffer_1100_ppdd[2609];

    auto g_z_z_0_0_x_x_yy_xx = buffer_1100_ppdd[2610];

    auto g_z_z_0_0_x_x_yy_xy = buffer_1100_ppdd[2611];

    auto g_z_z_0_0_x_x_yy_xz = buffer_1100_ppdd[2612];

    auto g_z_z_0_0_x_x_yy_yy = buffer_1100_ppdd[2613];

    auto g_z_z_0_0_x_x_yy_yz = buffer_1100_ppdd[2614];

    auto g_z_z_0_0_x_x_yy_zz = buffer_1100_ppdd[2615];

    auto g_z_z_0_0_x_x_yz_xx = buffer_1100_ppdd[2616];

    auto g_z_z_0_0_x_x_yz_xy = buffer_1100_ppdd[2617];

    auto g_z_z_0_0_x_x_yz_xz = buffer_1100_ppdd[2618];

    auto g_z_z_0_0_x_x_yz_yy = buffer_1100_ppdd[2619];

    auto g_z_z_0_0_x_x_yz_yz = buffer_1100_ppdd[2620];

    auto g_z_z_0_0_x_x_yz_zz = buffer_1100_ppdd[2621];

    auto g_z_z_0_0_x_x_zz_xx = buffer_1100_ppdd[2622];

    auto g_z_z_0_0_x_x_zz_xy = buffer_1100_ppdd[2623];

    auto g_z_z_0_0_x_x_zz_xz = buffer_1100_ppdd[2624];

    auto g_z_z_0_0_x_x_zz_yy = buffer_1100_ppdd[2625];

    auto g_z_z_0_0_x_x_zz_yz = buffer_1100_ppdd[2626];

    auto g_z_z_0_0_x_x_zz_zz = buffer_1100_ppdd[2627];

    auto g_z_z_0_0_x_y_xx_xx = buffer_1100_ppdd[2628];

    auto g_z_z_0_0_x_y_xx_xy = buffer_1100_ppdd[2629];

    auto g_z_z_0_0_x_y_xx_xz = buffer_1100_ppdd[2630];

    auto g_z_z_0_0_x_y_xx_yy = buffer_1100_ppdd[2631];

    auto g_z_z_0_0_x_y_xx_yz = buffer_1100_ppdd[2632];

    auto g_z_z_0_0_x_y_xx_zz = buffer_1100_ppdd[2633];

    auto g_z_z_0_0_x_y_xy_xx = buffer_1100_ppdd[2634];

    auto g_z_z_0_0_x_y_xy_xy = buffer_1100_ppdd[2635];

    auto g_z_z_0_0_x_y_xy_xz = buffer_1100_ppdd[2636];

    auto g_z_z_0_0_x_y_xy_yy = buffer_1100_ppdd[2637];

    auto g_z_z_0_0_x_y_xy_yz = buffer_1100_ppdd[2638];

    auto g_z_z_0_0_x_y_xy_zz = buffer_1100_ppdd[2639];

    auto g_z_z_0_0_x_y_xz_xx = buffer_1100_ppdd[2640];

    auto g_z_z_0_0_x_y_xz_xy = buffer_1100_ppdd[2641];

    auto g_z_z_0_0_x_y_xz_xz = buffer_1100_ppdd[2642];

    auto g_z_z_0_0_x_y_xz_yy = buffer_1100_ppdd[2643];

    auto g_z_z_0_0_x_y_xz_yz = buffer_1100_ppdd[2644];

    auto g_z_z_0_0_x_y_xz_zz = buffer_1100_ppdd[2645];

    auto g_z_z_0_0_x_y_yy_xx = buffer_1100_ppdd[2646];

    auto g_z_z_0_0_x_y_yy_xy = buffer_1100_ppdd[2647];

    auto g_z_z_0_0_x_y_yy_xz = buffer_1100_ppdd[2648];

    auto g_z_z_0_0_x_y_yy_yy = buffer_1100_ppdd[2649];

    auto g_z_z_0_0_x_y_yy_yz = buffer_1100_ppdd[2650];

    auto g_z_z_0_0_x_y_yy_zz = buffer_1100_ppdd[2651];

    auto g_z_z_0_0_x_y_yz_xx = buffer_1100_ppdd[2652];

    auto g_z_z_0_0_x_y_yz_xy = buffer_1100_ppdd[2653];

    auto g_z_z_0_0_x_y_yz_xz = buffer_1100_ppdd[2654];

    auto g_z_z_0_0_x_y_yz_yy = buffer_1100_ppdd[2655];

    auto g_z_z_0_0_x_y_yz_yz = buffer_1100_ppdd[2656];

    auto g_z_z_0_0_x_y_yz_zz = buffer_1100_ppdd[2657];

    auto g_z_z_0_0_x_y_zz_xx = buffer_1100_ppdd[2658];

    auto g_z_z_0_0_x_y_zz_xy = buffer_1100_ppdd[2659];

    auto g_z_z_0_0_x_y_zz_xz = buffer_1100_ppdd[2660];

    auto g_z_z_0_0_x_y_zz_yy = buffer_1100_ppdd[2661];

    auto g_z_z_0_0_x_y_zz_yz = buffer_1100_ppdd[2662];

    auto g_z_z_0_0_x_y_zz_zz = buffer_1100_ppdd[2663];

    auto g_z_z_0_0_x_z_xx_xx = buffer_1100_ppdd[2664];

    auto g_z_z_0_0_x_z_xx_xy = buffer_1100_ppdd[2665];

    auto g_z_z_0_0_x_z_xx_xz = buffer_1100_ppdd[2666];

    auto g_z_z_0_0_x_z_xx_yy = buffer_1100_ppdd[2667];

    auto g_z_z_0_0_x_z_xx_yz = buffer_1100_ppdd[2668];

    auto g_z_z_0_0_x_z_xx_zz = buffer_1100_ppdd[2669];

    auto g_z_z_0_0_x_z_xy_xx = buffer_1100_ppdd[2670];

    auto g_z_z_0_0_x_z_xy_xy = buffer_1100_ppdd[2671];

    auto g_z_z_0_0_x_z_xy_xz = buffer_1100_ppdd[2672];

    auto g_z_z_0_0_x_z_xy_yy = buffer_1100_ppdd[2673];

    auto g_z_z_0_0_x_z_xy_yz = buffer_1100_ppdd[2674];

    auto g_z_z_0_0_x_z_xy_zz = buffer_1100_ppdd[2675];

    auto g_z_z_0_0_x_z_xz_xx = buffer_1100_ppdd[2676];

    auto g_z_z_0_0_x_z_xz_xy = buffer_1100_ppdd[2677];

    auto g_z_z_0_0_x_z_xz_xz = buffer_1100_ppdd[2678];

    auto g_z_z_0_0_x_z_xz_yy = buffer_1100_ppdd[2679];

    auto g_z_z_0_0_x_z_xz_yz = buffer_1100_ppdd[2680];

    auto g_z_z_0_0_x_z_xz_zz = buffer_1100_ppdd[2681];

    auto g_z_z_0_0_x_z_yy_xx = buffer_1100_ppdd[2682];

    auto g_z_z_0_0_x_z_yy_xy = buffer_1100_ppdd[2683];

    auto g_z_z_0_0_x_z_yy_xz = buffer_1100_ppdd[2684];

    auto g_z_z_0_0_x_z_yy_yy = buffer_1100_ppdd[2685];

    auto g_z_z_0_0_x_z_yy_yz = buffer_1100_ppdd[2686];

    auto g_z_z_0_0_x_z_yy_zz = buffer_1100_ppdd[2687];

    auto g_z_z_0_0_x_z_yz_xx = buffer_1100_ppdd[2688];

    auto g_z_z_0_0_x_z_yz_xy = buffer_1100_ppdd[2689];

    auto g_z_z_0_0_x_z_yz_xz = buffer_1100_ppdd[2690];

    auto g_z_z_0_0_x_z_yz_yy = buffer_1100_ppdd[2691];

    auto g_z_z_0_0_x_z_yz_yz = buffer_1100_ppdd[2692];

    auto g_z_z_0_0_x_z_yz_zz = buffer_1100_ppdd[2693];

    auto g_z_z_0_0_x_z_zz_xx = buffer_1100_ppdd[2694];

    auto g_z_z_0_0_x_z_zz_xy = buffer_1100_ppdd[2695];

    auto g_z_z_0_0_x_z_zz_xz = buffer_1100_ppdd[2696];

    auto g_z_z_0_0_x_z_zz_yy = buffer_1100_ppdd[2697];

    auto g_z_z_0_0_x_z_zz_yz = buffer_1100_ppdd[2698];

    auto g_z_z_0_0_x_z_zz_zz = buffer_1100_ppdd[2699];

    auto g_z_z_0_0_y_x_xx_xx = buffer_1100_ppdd[2700];

    auto g_z_z_0_0_y_x_xx_xy = buffer_1100_ppdd[2701];

    auto g_z_z_0_0_y_x_xx_xz = buffer_1100_ppdd[2702];

    auto g_z_z_0_0_y_x_xx_yy = buffer_1100_ppdd[2703];

    auto g_z_z_0_0_y_x_xx_yz = buffer_1100_ppdd[2704];

    auto g_z_z_0_0_y_x_xx_zz = buffer_1100_ppdd[2705];

    auto g_z_z_0_0_y_x_xy_xx = buffer_1100_ppdd[2706];

    auto g_z_z_0_0_y_x_xy_xy = buffer_1100_ppdd[2707];

    auto g_z_z_0_0_y_x_xy_xz = buffer_1100_ppdd[2708];

    auto g_z_z_0_0_y_x_xy_yy = buffer_1100_ppdd[2709];

    auto g_z_z_0_0_y_x_xy_yz = buffer_1100_ppdd[2710];

    auto g_z_z_0_0_y_x_xy_zz = buffer_1100_ppdd[2711];

    auto g_z_z_0_0_y_x_xz_xx = buffer_1100_ppdd[2712];

    auto g_z_z_0_0_y_x_xz_xy = buffer_1100_ppdd[2713];

    auto g_z_z_0_0_y_x_xz_xz = buffer_1100_ppdd[2714];

    auto g_z_z_0_0_y_x_xz_yy = buffer_1100_ppdd[2715];

    auto g_z_z_0_0_y_x_xz_yz = buffer_1100_ppdd[2716];

    auto g_z_z_0_0_y_x_xz_zz = buffer_1100_ppdd[2717];

    auto g_z_z_0_0_y_x_yy_xx = buffer_1100_ppdd[2718];

    auto g_z_z_0_0_y_x_yy_xy = buffer_1100_ppdd[2719];

    auto g_z_z_0_0_y_x_yy_xz = buffer_1100_ppdd[2720];

    auto g_z_z_0_0_y_x_yy_yy = buffer_1100_ppdd[2721];

    auto g_z_z_0_0_y_x_yy_yz = buffer_1100_ppdd[2722];

    auto g_z_z_0_0_y_x_yy_zz = buffer_1100_ppdd[2723];

    auto g_z_z_0_0_y_x_yz_xx = buffer_1100_ppdd[2724];

    auto g_z_z_0_0_y_x_yz_xy = buffer_1100_ppdd[2725];

    auto g_z_z_0_0_y_x_yz_xz = buffer_1100_ppdd[2726];

    auto g_z_z_0_0_y_x_yz_yy = buffer_1100_ppdd[2727];

    auto g_z_z_0_0_y_x_yz_yz = buffer_1100_ppdd[2728];

    auto g_z_z_0_0_y_x_yz_zz = buffer_1100_ppdd[2729];

    auto g_z_z_0_0_y_x_zz_xx = buffer_1100_ppdd[2730];

    auto g_z_z_0_0_y_x_zz_xy = buffer_1100_ppdd[2731];

    auto g_z_z_0_0_y_x_zz_xz = buffer_1100_ppdd[2732];

    auto g_z_z_0_0_y_x_zz_yy = buffer_1100_ppdd[2733];

    auto g_z_z_0_0_y_x_zz_yz = buffer_1100_ppdd[2734];

    auto g_z_z_0_0_y_x_zz_zz = buffer_1100_ppdd[2735];

    auto g_z_z_0_0_y_y_xx_xx = buffer_1100_ppdd[2736];

    auto g_z_z_0_0_y_y_xx_xy = buffer_1100_ppdd[2737];

    auto g_z_z_0_0_y_y_xx_xz = buffer_1100_ppdd[2738];

    auto g_z_z_0_0_y_y_xx_yy = buffer_1100_ppdd[2739];

    auto g_z_z_0_0_y_y_xx_yz = buffer_1100_ppdd[2740];

    auto g_z_z_0_0_y_y_xx_zz = buffer_1100_ppdd[2741];

    auto g_z_z_0_0_y_y_xy_xx = buffer_1100_ppdd[2742];

    auto g_z_z_0_0_y_y_xy_xy = buffer_1100_ppdd[2743];

    auto g_z_z_0_0_y_y_xy_xz = buffer_1100_ppdd[2744];

    auto g_z_z_0_0_y_y_xy_yy = buffer_1100_ppdd[2745];

    auto g_z_z_0_0_y_y_xy_yz = buffer_1100_ppdd[2746];

    auto g_z_z_0_0_y_y_xy_zz = buffer_1100_ppdd[2747];

    auto g_z_z_0_0_y_y_xz_xx = buffer_1100_ppdd[2748];

    auto g_z_z_0_0_y_y_xz_xy = buffer_1100_ppdd[2749];

    auto g_z_z_0_0_y_y_xz_xz = buffer_1100_ppdd[2750];

    auto g_z_z_0_0_y_y_xz_yy = buffer_1100_ppdd[2751];

    auto g_z_z_0_0_y_y_xz_yz = buffer_1100_ppdd[2752];

    auto g_z_z_0_0_y_y_xz_zz = buffer_1100_ppdd[2753];

    auto g_z_z_0_0_y_y_yy_xx = buffer_1100_ppdd[2754];

    auto g_z_z_0_0_y_y_yy_xy = buffer_1100_ppdd[2755];

    auto g_z_z_0_0_y_y_yy_xz = buffer_1100_ppdd[2756];

    auto g_z_z_0_0_y_y_yy_yy = buffer_1100_ppdd[2757];

    auto g_z_z_0_0_y_y_yy_yz = buffer_1100_ppdd[2758];

    auto g_z_z_0_0_y_y_yy_zz = buffer_1100_ppdd[2759];

    auto g_z_z_0_0_y_y_yz_xx = buffer_1100_ppdd[2760];

    auto g_z_z_0_0_y_y_yz_xy = buffer_1100_ppdd[2761];

    auto g_z_z_0_0_y_y_yz_xz = buffer_1100_ppdd[2762];

    auto g_z_z_0_0_y_y_yz_yy = buffer_1100_ppdd[2763];

    auto g_z_z_0_0_y_y_yz_yz = buffer_1100_ppdd[2764];

    auto g_z_z_0_0_y_y_yz_zz = buffer_1100_ppdd[2765];

    auto g_z_z_0_0_y_y_zz_xx = buffer_1100_ppdd[2766];

    auto g_z_z_0_0_y_y_zz_xy = buffer_1100_ppdd[2767];

    auto g_z_z_0_0_y_y_zz_xz = buffer_1100_ppdd[2768];

    auto g_z_z_0_0_y_y_zz_yy = buffer_1100_ppdd[2769];

    auto g_z_z_0_0_y_y_zz_yz = buffer_1100_ppdd[2770];

    auto g_z_z_0_0_y_y_zz_zz = buffer_1100_ppdd[2771];

    auto g_z_z_0_0_y_z_xx_xx = buffer_1100_ppdd[2772];

    auto g_z_z_0_0_y_z_xx_xy = buffer_1100_ppdd[2773];

    auto g_z_z_0_0_y_z_xx_xz = buffer_1100_ppdd[2774];

    auto g_z_z_0_0_y_z_xx_yy = buffer_1100_ppdd[2775];

    auto g_z_z_0_0_y_z_xx_yz = buffer_1100_ppdd[2776];

    auto g_z_z_0_0_y_z_xx_zz = buffer_1100_ppdd[2777];

    auto g_z_z_0_0_y_z_xy_xx = buffer_1100_ppdd[2778];

    auto g_z_z_0_0_y_z_xy_xy = buffer_1100_ppdd[2779];

    auto g_z_z_0_0_y_z_xy_xz = buffer_1100_ppdd[2780];

    auto g_z_z_0_0_y_z_xy_yy = buffer_1100_ppdd[2781];

    auto g_z_z_0_0_y_z_xy_yz = buffer_1100_ppdd[2782];

    auto g_z_z_0_0_y_z_xy_zz = buffer_1100_ppdd[2783];

    auto g_z_z_0_0_y_z_xz_xx = buffer_1100_ppdd[2784];

    auto g_z_z_0_0_y_z_xz_xy = buffer_1100_ppdd[2785];

    auto g_z_z_0_0_y_z_xz_xz = buffer_1100_ppdd[2786];

    auto g_z_z_0_0_y_z_xz_yy = buffer_1100_ppdd[2787];

    auto g_z_z_0_0_y_z_xz_yz = buffer_1100_ppdd[2788];

    auto g_z_z_0_0_y_z_xz_zz = buffer_1100_ppdd[2789];

    auto g_z_z_0_0_y_z_yy_xx = buffer_1100_ppdd[2790];

    auto g_z_z_0_0_y_z_yy_xy = buffer_1100_ppdd[2791];

    auto g_z_z_0_0_y_z_yy_xz = buffer_1100_ppdd[2792];

    auto g_z_z_0_0_y_z_yy_yy = buffer_1100_ppdd[2793];

    auto g_z_z_0_0_y_z_yy_yz = buffer_1100_ppdd[2794];

    auto g_z_z_0_0_y_z_yy_zz = buffer_1100_ppdd[2795];

    auto g_z_z_0_0_y_z_yz_xx = buffer_1100_ppdd[2796];

    auto g_z_z_0_0_y_z_yz_xy = buffer_1100_ppdd[2797];

    auto g_z_z_0_0_y_z_yz_xz = buffer_1100_ppdd[2798];

    auto g_z_z_0_0_y_z_yz_yy = buffer_1100_ppdd[2799];

    auto g_z_z_0_0_y_z_yz_yz = buffer_1100_ppdd[2800];

    auto g_z_z_0_0_y_z_yz_zz = buffer_1100_ppdd[2801];

    auto g_z_z_0_0_y_z_zz_xx = buffer_1100_ppdd[2802];

    auto g_z_z_0_0_y_z_zz_xy = buffer_1100_ppdd[2803];

    auto g_z_z_0_0_y_z_zz_xz = buffer_1100_ppdd[2804];

    auto g_z_z_0_0_y_z_zz_yy = buffer_1100_ppdd[2805];

    auto g_z_z_0_0_y_z_zz_yz = buffer_1100_ppdd[2806];

    auto g_z_z_0_0_y_z_zz_zz = buffer_1100_ppdd[2807];

    auto g_z_z_0_0_z_x_xx_xx = buffer_1100_ppdd[2808];

    auto g_z_z_0_0_z_x_xx_xy = buffer_1100_ppdd[2809];

    auto g_z_z_0_0_z_x_xx_xz = buffer_1100_ppdd[2810];

    auto g_z_z_0_0_z_x_xx_yy = buffer_1100_ppdd[2811];

    auto g_z_z_0_0_z_x_xx_yz = buffer_1100_ppdd[2812];

    auto g_z_z_0_0_z_x_xx_zz = buffer_1100_ppdd[2813];

    auto g_z_z_0_0_z_x_xy_xx = buffer_1100_ppdd[2814];

    auto g_z_z_0_0_z_x_xy_xy = buffer_1100_ppdd[2815];

    auto g_z_z_0_0_z_x_xy_xz = buffer_1100_ppdd[2816];

    auto g_z_z_0_0_z_x_xy_yy = buffer_1100_ppdd[2817];

    auto g_z_z_0_0_z_x_xy_yz = buffer_1100_ppdd[2818];

    auto g_z_z_0_0_z_x_xy_zz = buffer_1100_ppdd[2819];

    auto g_z_z_0_0_z_x_xz_xx = buffer_1100_ppdd[2820];

    auto g_z_z_0_0_z_x_xz_xy = buffer_1100_ppdd[2821];

    auto g_z_z_0_0_z_x_xz_xz = buffer_1100_ppdd[2822];

    auto g_z_z_0_0_z_x_xz_yy = buffer_1100_ppdd[2823];

    auto g_z_z_0_0_z_x_xz_yz = buffer_1100_ppdd[2824];

    auto g_z_z_0_0_z_x_xz_zz = buffer_1100_ppdd[2825];

    auto g_z_z_0_0_z_x_yy_xx = buffer_1100_ppdd[2826];

    auto g_z_z_0_0_z_x_yy_xy = buffer_1100_ppdd[2827];

    auto g_z_z_0_0_z_x_yy_xz = buffer_1100_ppdd[2828];

    auto g_z_z_0_0_z_x_yy_yy = buffer_1100_ppdd[2829];

    auto g_z_z_0_0_z_x_yy_yz = buffer_1100_ppdd[2830];

    auto g_z_z_0_0_z_x_yy_zz = buffer_1100_ppdd[2831];

    auto g_z_z_0_0_z_x_yz_xx = buffer_1100_ppdd[2832];

    auto g_z_z_0_0_z_x_yz_xy = buffer_1100_ppdd[2833];

    auto g_z_z_0_0_z_x_yz_xz = buffer_1100_ppdd[2834];

    auto g_z_z_0_0_z_x_yz_yy = buffer_1100_ppdd[2835];

    auto g_z_z_0_0_z_x_yz_yz = buffer_1100_ppdd[2836];

    auto g_z_z_0_0_z_x_yz_zz = buffer_1100_ppdd[2837];

    auto g_z_z_0_0_z_x_zz_xx = buffer_1100_ppdd[2838];

    auto g_z_z_0_0_z_x_zz_xy = buffer_1100_ppdd[2839];

    auto g_z_z_0_0_z_x_zz_xz = buffer_1100_ppdd[2840];

    auto g_z_z_0_0_z_x_zz_yy = buffer_1100_ppdd[2841];

    auto g_z_z_0_0_z_x_zz_yz = buffer_1100_ppdd[2842];

    auto g_z_z_0_0_z_x_zz_zz = buffer_1100_ppdd[2843];

    auto g_z_z_0_0_z_y_xx_xx = buffer_1100_ppdd[2844];

    auto g_z_z_0_0_z_y_xx_xy = buffer_1100_ppdd[2845];

    auto g_z_z_0_0_z_y_xx_xz = buffer_1100_ppdd[2846];

    auto g_z_z_0_0_z_y_xx_yy = buffer_1100_ppdd[2847];

    auto g_z_z_0_0_z_y_xx_yz = buffer_1100_ppdd[2848];

    auto g_z_z_0_0_z_y_xx_zz = buffer_1100_ppdd[2849];

    auto g_z_z_0_0_z_y_xy_xx = buffer_1100_ppdd[2850];

    auto g_z_z_0_0_z_y_xy_xy = buffer_1100_ppdd[2851];

    auto g_z_z_0_0_z_y_xy_xz = buffer_1100_ppdd[2852];

    auto g_z_z_0_0_z_y_xy_yy = buffer_1100_ppdd[2853];

    auto g_z_z_0_0_z_y_xy_yz = buffer_1100_ppdd[2854];

    auto g_z_z_0_0_z_y_xy_zz = buffer_1100_ppdd[2855];

    auto g_z_z_0_0_z_y_xz_xx = buffer_1100_ppdd[2856];

    auto g_z_z_0_0_z_y_xz_xy = buffer_1100_ppdd[2857];

    auto g_z_z_0_0_z_y_xz_xz = buffer_1100_ppdd[2858];

    auto g_z_z_0_0_z_y_xz_yy = buffer_1100_ppdd[2859];

    auto g_z_z_0_0_z_y_xz_yz = buffer_1100_ppdd[2860];

    auto g_z_z_0_0_z_y_xz_zz = buffer_1100_ppdd[2861];

    auto g_z_z_0_0_z_y_yy_xx = buffer_1100_ppdd[2862];

    auto g_z_z_0_0_z_y_yy_xy = buffer_1100_ppdd[2863];

    auto g_z_z_0_0_z_y_yy_xz = buffer_1100_ppdd[2864];

    auto g_z_z_0_0_z_y_yy_yy = buffer_1100_ppdd[2865];

    auto g_z_z_0_0_z_y_yy_yz = buffer_1100_ppdd[2866];

    auto g_z_z_0_0_z_y_yy_zz = buffer_1100_ppdd[2867];

    auto g_z_z_0_0_z_y_yz_xx = buffer_1100_ppdd[2868];

    auto g_z_z_0_0_z_y_yz_xy = buffer_1100_ppdd[2869];

    auto g_z_z_0_0_z_y_yz_xz = buffer_1100_ppdd[2870];

    auto g_z_z_0_0_z_y_yz_yy = buffer_1100_ppdd[2871];

    auto g_z_z_0_0_z_y_yz_yz = buffer_1100_ppdd[2872];

    auto g_z_z_0_0_z_y_yz_zz = buffer_1100_ppdd[2873];

    auto g_z_z_0_0_z_y_zz_xx = buffer_1100_ppdd[2874];

    auto g_z_z_0_0_z_y_zz_xy = buffer_1100_ppdd[2875];

    auto g_z_z_0_0_z_y_zz_xz = buffer_1100_ppdd[2876];

    auto g_z_z_0_0_z_y_zz_yy = buffer_1100_ppdd[2877];

    auto g_z_z_0_0_z_y_zz_yz = buffer_1100_ppdd[2878];

    auto g_z_z_0_0_z_y_zz_zz = buffer_1100_ppdd[2879];

    auto g_z_z_0_0_z_z_xx_xx = buffer_1100_ppdd[2880];

    auto g_z_z_0_0_z_z_xx_xy = buffer_1100_ppdd[2881];

    auto g_z_z_0_0_z_z_xx_xz = buffer_1100_ppdd[2882];

    auto g_z_z_0_0_z_z_xx_yy = buffer_1100_ppdd[2883];

    auto g_z_z_0_0_z_z_xx_yz = buffer_1100_ppdd[2884];

    auto g_z_z_0_0_z_z_xx_zz = buffer_1100_ppdd[2885];

    auto g_z_z_0_0_z_z_xy_xx = buffer_1100_ppdd[2886];

    auto g_z_z_0_0_z_z_xy_xy = buffer_1100_ppdd[2887];

    auto g_z_z_0_0_z_z_xy_xz = buffer_1100_ppdd[2888];

    auto g_z_z_0_0_z_z_xy_yy = buffer_1100_ppdd[2889];

    auto g_z_z_0_0_z_z_xy_yz = buffer_1100_ppdd[2890];

    auto g_z_z_0_0_z_z_xy_zz = buffer_1100_ppdd[2891];

    auto g_z_z_0_0_z_z_xz_xx = buffer_1100_ppdd[2892];

    auto g_z_z_0_0_z_z_xz_xy = buffer_1100_ppdd[2893];

    auto g_z_z_0_0_z_z_xz_xz = buffer_1100_ppdd[2894];

    auto g_z_z_0_0_z_z_xz_yy = buffer_1100_ppdd[2895];

    auto g_z_z_0_0_z_z_xz_yz = buffer_1100_ppdd[2896];

    auto g_z_z_0_0_z_z_xz_zz = buffer_1100_ppdd[2897];

    auto g_z_z_0_0_z_z_yy_xx = buffer_1100_ppdd[2898];

    auto g_z_z_0_0_z_z_yy_xy = buffer_1100_ppdd[2899];

    auto g_z_z_0_0_z_z_yy_xz = buffer_1100_ppdd[2900];

    auto g_z_z_0_0_z_z_yy_yy = buffer_1100_ppdd[2901];

    auto g_z_z_0_0_z_z_yy_yz = buffer_1100_ppdd[2902];

    auto g_z_z_0_0_z_z_yy_zz = buffer_1100_ppdd[2903];

    auto g_z_z_0_0_z_z_yz_xx = buffer_1100_ppdd[2904];

    auto g_z_z_0_0_z_z_yz_xy = buffer_1100_ppdd[2905];

    auto g_z_z_0_0_z_z_yz_xz = buffer_1100_ppdd[2906];

    auto g_z_z_0_0_z_z_yz_yy = buffer_1100_ppdd[2907];

    auto g_z_z_0_0_z_z_yz_yz = buffer_1100_ppdd[2908];

    auto g_z_z_0_0_z_z_yz_zz = buffer_1100_ppdd[2909];

    auto g_z_z_0_0_z_z_zz_xx = buffer_1100_ppdd[2910];

    auto g_z_z_0_0_z_z_zz_xy = buffer_1100_ppdd[2911];

    auto g_z_z_0_0_z_z_zz_xz = buffer_1100_ppdd[2912];

    auto g_z_z_0_0_z_z_zz_yy = buffer_1100_ppdd[2913];

    auto g_z_z_0_0_z_z_zz_yz = buffer_1100_ppdd[2914];

    auto g_z_z_0_0_z_z_zz_zz = buffer_1100_ppdd[2915];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_0_xx_xx, g_0_0_xx_xy, g_0_0_xx_xz, g_0_0_xx_yy, g_0_0_xx_yz, g_0_0_xx_zz, g_0_xx_xx_xx, g_0_xx_xx_xy, g_0_xx_xx_xz, g_0_xx_xx_yy, g_0_xx_xx_yz, g_0_xx_xx_zz, g_x_x_0_0_x_x_xx_xx, g_x_x_0_0_x_x_xx_xy, g_x_x_0_0_x_x_xx_xz, g_x_x_0_0_x_x_xx_yy, g_x_x_0_0_x_x_xx_yz, g_x_x_0_0_x_x_xx_zz, g_xx_0_xx_xx, g_xx_0_xx_xy, g_xx_0_xx_xz, g_xx_0_xx_yy, g_xx_0_xx_yz, g_xx_0_xx_zz, g_xx_xx_xx_xx, g_xx_xx_xx_xy, g_xx_xx_xx_xz, g_xx_xx_xx_yy, g_xx_xx_xx_yz, g_xx_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_x_xx_xx[i] = g_0_0_xx_xx[i] - 2.0 * g_0_xx_xx_xx[i] * b_exp - 2.0 * g_xx_0_xx_xx[i] * a_exp + 4.0 * g_xx_xx_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_xx_xy[i] = g_0_0_xx_xy[i] - 2.0 * g_0_xx_xx_xy[i] * b_exp - 2.0 * g_xx_0_xx_xy[i] * a_exp + 4.0 * g_xx_xx_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_xx_xz[i] = g_0_0_xx_xz[i] - 2.0 * g_0_xx_xx_xz[i] * b_exp - 2.0 * g_xx_0_xx_xz[i] * a_exp + 4.0 * g_xx_xx_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_xx_yy[i] = g_0_0_xx_yy[i] - 2.0 * g_0_xx_xx_yy[i] * b_exp - 2.0 * g_xx_0_xx_yy[i] * a_exp + 4.0 * g_xx_xx_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_xx_yz[i] = g_0_0_xx_yz[i] - 2.0 * g_0_xx_xx_yz[i] * b_exp - 2.0 * g_xx_0_xx_yz[i] * a_exp + 4.0 * g_xx_xx_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_xx_zz[i] = g_0_0_xx_zz[i] - 2.0 * g_0_xx_xx_zz[i] * b_exp - 2.0 * g_xx_0_xx_zz[i] * a_exp + 4.0 * g_xx_xx_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_0_xy_xx, g_0_0_xy_xy, g_0_0_xy_xz, g_0_0_xy_yy, g_0_0_xy_yz, g_0_0_xy_zz, g_0_xx_xy_xx, g_0_xx_xy_xy, g_0_xx_xy_xz, g_0_xx_xy_yy, g_0_xx_xy_yz, g_0_xx_xy_zz, g_x_x_0_0_x_x_xy_xx, g_x_x_0_0_x_x_xy_xy, g_x_x_0_0_x_x_xy_xz, g_x_x_0_0_x_x_xy_yy, g_x_x_0_0_x_x_xy_yz, g_x_x_0_0_x_x_xy_zz, g_xx_0_xy_xx, g_xx_0_xy_xy, g_xx_0_xy_xz, g_xx_0_xy_yy, g_xx_0_xy_yz, g_xx_0_xy_zz, g_xx_xx_xy_xx, g_xx_xx_xy_xy, g_xx_xx_xy_xz, g_xx_xx_xy_yy, g_xx_xx_xy_yz, g_xx_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_x_xy_xx[i] = g_0_0_xy_xx[i] - 2.0 * g_0_xx_xy_xx[i] * b_exp - 2.0 * g_xx_0_xy_xx[i] * a_exp + 4.0 * g_xx_xx_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_xy_xy[i] = g_0_0_xy_xy[i] - 2.0 * g_0_xx_xy_xy[i] * b_exp - 2.0 * g_xx_0_xy_xy[i] * a_exp + 4.0 * g_xx_xx_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_xy_xz[i] = g_0_0_xy_xz[i] - 2.0 * g_0_xx_xy_xz[i] * b_exp - 2.0 * g_xx_0_xy_xz[i] * a_exp + 4.0 * g_xx_xx_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_xy_yy[i] = g_0_0_xy_yy[i] - 2.0 * g_0_xx_xy_yy[i] * b_exp - 2.0 * g_xx_0_xy_yy[i] * a_exp + 4.0 * g_xx_xx_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_xy_yz[i] = g_0_0_xy_yz[i] - 2.0 * g_0_xx_xy_yz[i] * b_exp - 2.0 * g_xx_0_xy_yz[i] * a_exp + 4.0 * g_xx_xx_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_xy_zz[i] = g_0_0_xy_zz[i] - 2.0 * g_0_xx_xy_zz[i] * b_exp - 2.0 * g_xx_0_xy_zz[i] * a_exp + 4.0 * g_xx_xx_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_0_xz_xx, g_0_0_xz_xy, g_0_0_xz_xz, g_0_0_xz_yy, g_0_0_xz_yz, g_0_0_xz_zz, g_0_xx_xz_xx, g_0_xx_xz_xy, g_0_xx_xz_xz, g_0_xx_xz_yy, g_0_xx_xz_yz, g_0_xx_xz_zz, g_x_x_0_0_x_x_xz_xx, g_x_x_0_0_x_x_xz_xy, g_x_x_0_0_x_x_xz_xz, g_x_x_0_0_x_x_xz_yy, g_x_x_0_0_x_x_xz_yz, g_x_x_0_0_x_x_xz_zz, g_xx_0_xz_xx, g_xx_0_xz_xy, g_xx_0_xz_xz, g_xx_0_xz_yy, g_xx_0_xz_yz, g_xx_0_xz_zz, g_xx_xx_xz_xx, g_xx_xx_xz_xy, g_xx_xx_xz_xz, g_xx_xx_xz_yy, g_xx_xx_xz_yz, g_xx_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_x_xz_xx[i] = g_0_0_xz_xx[i] - 2.0 * g_0_xx_xz_xx[i] * b_exp - 2.0 * g_xx_0_xz_xx[i] * a_exp + 4.0 * g_xx_xx_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_xz_xy[i] = g_0_0_xz_xy[i] - 2.0 * g_0_xx_xz_xy[i] * b_exp - 2.0 * g_xx_0_xz_xy[i] * a_exp + 4.0 * g_xx_xx_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_xz_xz[i] = g_0_0_xz_xz[i] - 2.0 * g_0_xx_xz_xz[i] * b_exp - 2.0 * g_xx_0_xz_xz[i] * a_exp + 4.0 * g_xx_xx_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_xz_yy[i] = g_0_0_xz_yy[i] - 2.0 * g_0_xx_xz_yy[i] * b_exp - 2.0 * g_xx_0_xz_yy[i] * a_exp + 4.0 * g_xx_xx_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_xz_yz[i] = g_0_0_xz_yz[i] - 2.0 * g_0_xx_xz_yz[i] * b_exp - 2.0 * g_xx_0_xz_yz[i] * a_exp + 4.0 * g_xx_xx_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_xz_zz[i] = g_0_0_xz_zz[i] - 2.0 * g_0_xx_xz_zz[i] * b_exp - 2.0 * g_xx_0_xz_zz[i] * a_exp + 4.0 * g_xx_xx_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_0_yy_xx, g_0_0_yy_xy, g_0_0_yy_xz, g_0_0_yy_yy, g_0_0_yy_yz, g_0_0_yy_zz, g_0_xx_yy_xx, g_0_xx_yy_xy, g_0_xx_yy_xz, g_0_xx_yy_yy, g_0_xx_yy_yz, g_0_xx_yy_zz, g_x_x_0_0_x_x_yy_xx, g_x_x_0_0_x_x_yy_xy, g_x_x_0_0_x_x_yy_xz, g_x_x_0_0_x_x_yy_yy, g_x_x_0_0_x_x_yy_yz, g_x_x_0_0_x_x_yy_zz, g_xx_0_yy_xx, g_xx_0_yy_xy, g_xx_0_yy_xz, g_xx_0_yy_yy, g_xx_0_yy_yz, g_xx_0_yy_zz, g_xx_xx_yy_xx, g_xx_xx_yy_xy, g_xx_xx_yy_xz, g_xx_xx_yy_yy, g_xx_xx_yy_yz, g_xx_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_x_yy_xx[i] = g_0_0_yy_xx[i] - 2.0 * g_0_xx_yy_xx[i] * b_exp - 2.0 * g_xx_0_yy_xx[i] * a_exp + 4.0 * g_xx_xx_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_yy_xy[i] = g_0_0_yy_xy[i] - 2.0 * g_0_xx_yy_xy[i] * b_exp - 2.0 * g_xx_0_yy_xy[i] * a_exp + 4.0 * g_xx_xx_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_yy_xz[i] = g_0_0_yy_xz[i] - 2.0 * g_0_xx_yy_xz[i] * b_exp - 2.0 * g_xx_0_yy_xz[i] * a_exp + 4.0 * g_xx_xx_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_yy_yy[i] = g_0_0_yy_yy[i] - 2.0 * g_0_xx_yy_yy[i] * b_exp - 2.0 * g_xx_0_yy_yy[i] * a_exp + 4.0 * g_xx_xx_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_yy_yz[i] = g_0_0_yy_yz[i] - 2.0 * g_0_xx_yy_yz[i] * b_exp - 2.0 * g_xx_0_yy_yz[i] * a_exp + 4.0 * g_xx_xx_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_yy_zz[i] = g_0_0_yy_zz[i] - 2.0 * g_0_xx_yy_zz[i] * b_exp - 2.0 * g_xx_0_yy_zz[i] * a_exp + 4.0 * g_xx_xx_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_0_yz_xx, g_0_0_yz_xy, g_0_0_yz_xz, g_0_0_yz_yy, g_0_0_yz_yz, g_0_0_yz_zz, g_0_xx_yz_xx, g_0_xx_yz_xy, g_0_xx_yz_xz, g_0_xx_yz_yy, g_0_xx_yz_yz, g_0_xx_yz_zz, g_x_x_0_0_x_x_yz_xx, g_x_x_0_0_x_x_yz_xy, g_x_x_0_0_x_x_yz_xz, g_x_x_0_0_x_x_yz_yy, g_x_x_0_0_x_x_yz_yz, g_x_x_0_0_x_x_yz_zz, g_xx_0_yz_xx, g_xx_0_yz_xy, g_xx_0_yz_xz, g_xx_0_yz_yy, g_xx_0_yz_yz, g_xx_0_yz_zz, g_xx_xx_yz_xx, g_xx_xx_yz_xy, g_xx_xx_yz_xz, g_xx_xx_yz_yy, g_xx_xx_yz_yz, g_xx_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_x_yz_xx[i] = g_0_0_yz_xx[i] - 2.0 * g_0_xx_yz_xx[i] * b_exp - 2.0 * g_xx_0_yz_xx[i] * a_exp + 4.0 * g_xx_xx_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_yz_xy[i] = g_0_0_yz_xy[i] - 2.0 * g_0_xx_yz_xy[i] * b_exp - 2.0 * g_xx_0_yz_xy[i] * a_exp + 4.0 * g_xx_xx_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_yz_xz[i] = g_0_0_yz_xz[i] - 2.0 * g_0_xx_yz_xz[i] * b_exp - 2.0 * g_xx_0_yz_xz[i] * a_exp + 4.0 * g_xx_xx_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_yz_yy[i] = g_0_0_yz_yy[i] - 2.0 * g_0_xx_yz_yy[i] * b_exp - 2.0 * g_xx_0_yz_yy[i] * a_exp + 4.0 * g_xx_xx_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_yz_yz[i] = g_0_0_yz_yz[i] - 2.0 * g_0_xx_yz_yz[i] * b_exp - 2.0 * g_xx_0_yz_yz[i] * a_exp + 4.0 * g_xx_xx_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_yz_zz[i] = g_0_0_yz_zz[i] - 2.0 * g_0_xx_yz_zz[i] * b_exp - 2.0 * g_xx_0_yz_zz[i] * a_exp + 4.0 * g_xx_xx_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_0_zz_xx, g_0_0_zz_xy, g_0_0_zz_xz, g_0_0_zz_yy, g_0_0_zz_yz, g_0_0_zz_zz, g_0_xx_zz_xx, g_0_xx_zz_xy, g_0_xx_zz_xz, g_0_xx_zz_yy, g_0_xx_zz_yz, g_0_xx_zz_zz, g_x_x_0_0_x_x_zz_xx, g_x_x_0_0_x_x_zz_xy, g_x_x_0_0_x_x_zz_xz, g_x_x_0_0_x_x_zz_yy, g_x_x_0_0_x_x_zz_yz, g_x_x_0_0_x_x_zz_zz, g_xx_0_zz_xx, g_xx_0_zz_xy, g_xx_0_zz_xz, g_xx_0_zz_yy, g_xx_0_zz_yz, g_xx_0_zz_zz, g_xx_xx_zz_xx, g_xx_xx_zz_xy, g_xx_xx_zz_xz, g_xx_xx_zz_yy, g_xx_xx_zz_yz, g_xx_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_x_zz_xx[i] = g_0_0_zz_xx[i] - 2.0 * g_0_xx_zz_xx[i] * b_exp - 2.0 * g_xx_0_zz_xx[i] * a_exp + 4.0 * g_xx_xx_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_zz_xy[i] = g_0_0_zz_xy[i] - 2.0 * g_0_xx_zz_xy[i] * b_exp - 2.0 * g_xx_0_zz_xy[i] * a_exp + 4.0 * g_xx_xx_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_zz_xz[i] = g_0_0_zz_xz[i] - 2.0 * g_0_xx_zz_xz[i] * b_exp - 2.0 * g_xx_0_zz_xz[i] * a_exp + 4.0 * g_xx_xx_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_zz_yy[i] = g_0_0_zz_yy[i] - 2.0 * g_0_xx_zz_yy[i] * b_exp - 2.0 * g_xx_0_zz_yy[i] * a_exp + 4.0 * g_xx_xx_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_zz_yz[i] = g_0_0_zz_yz[i] - 2.0 * g_0_xx_zz_yz[i] * b_exp - 2.0 * g_xx_0_zz_yz[i] * a_exp + 4.0 * g_xx_xx_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_zz_zz[i] = g_0_0_zz_zz[i] - 2.0 * g_0_xx_zz_zz[i] * b_exp - 2.0 * g_xx_0_zz_zz[i] * a_exp + 4.0 * g_xx_xx_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_0_xy_xx_xx, g_0_xy_xx_xy, g_0_xy_xx_xz, g_0_xy_xx_yy, g_0_xy_xx_yz, g_0_xy_xx_zz, g_x_x_0_0_x_y_xx_xx, g_x_x_0_0_x_y_xx_xy, g_x_x_0_0_x_y_xx_xz, g_x_x_0_0_x_y_xx_yy, g_x_x_0_0_x_y_xx_yz, g_x_x_0_0_x_y_xx_zz, g_xx_xy_xx_xx, g_xx_xy_xx_xy, g_xx_xy_xx_xz, g_xx_xy_xx_yy, g_xx_xy_xx_yz, g_xx_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_y_xx_xx[i] = -2.0 * g_0_xy_xx_xx[i] * b_exp + 4.0 * g_xx_xy_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_xx_xy[i] = -2.0 * g_0_xy_xx_xy[i] * b_exp + 4.0 * g_xx_xy_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_xx_xz[i] = -2.0 * g_0_xy_xx_xz[i] * b_exp + 4.0 * g_xx_xy_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_xx_yy[i] = -2.0 * g_0_xy_xx_yy[i] * b_exp + 4.0 * g_xx_xy_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_xx_yz[i] = -2.0 * g_0_xy_xx_yz[i] * b_exp + 4.0 * g_xx_xy_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_xx_zz[i] = -2.0 * g_0_xy_xx_zz[i] * b_exp + 4.0 * g_xx_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_0_xy_xy_xx, g_0_xy_xy_xy, g_0_xy_xy_xz, g_0_xy_xy_yy, g_0_xy_xy_yz, g_0_xy_xy_zz, g_x_x_0_0_x_y_xy_xx, g_x_x_0_0_x_y_xy_xy, g_x_x_0_0_x_y_xy_xz, g_x_x_0_0_x_y_xy_yy, g_x_x_0_0_x_y_xy_yz, g_x_x_0_0_x_y_xy_zz, g_xx_xy_xy_xx, g_xx_xy_xy_xy, g_xx_xy_xy_xz, g_xx_xy_xy_yy, g_xx_xy_xy_yz, g_xx_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_y_xy_xx[i] = -2.0 * g_0_xy_xy_xx[i] * b_exp + 4.0 * g_xx_xy_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_xy_xy[i] = -2.0 * g_0_xy_xy_xy[i] * b_exp + 4.0 * g_xx_xy_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_xy_xz[i] = -2.0 * g_0_xy_xy_xz[i] * b_exp + 4.0 * g_xx_xy_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_xy_yy[i] = -2.0 * g_0_xy_xy_yy[i] * b_exp + 4.0 * g_xx_xy_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_xy_yz[i] = -2.0 * g_0_xy_xy_yz[i] * b_exp + 4.0 * g_xx_xy_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_xy_zz[i] = -2.0 * g_0_xy_xy_zz[i] * b_exp + 4.0 * g_xx_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_0_xy_xz_xx, g_0_xy_xz_xy, g_0_xy_xz_xz, g_0_xy_xz_yy, g_0_xy_xz_yz, g_0_xy_xz_zz, g_x_x_0_0_x_y_xz_xx, g_x_x_0_0_x_y_xz_xy, g_x_x_0_0_x_y_xz_xz, g_x_x_0_0_x_y_xz_yy, g_x_x_0_0_x_y_xz_yz, g_x_x_0_0_x_y_xz_zz, g_xx_xy_xz_xx, g_xx_xy_xz_xy, g_xx_xy_xz_xz, g_xx_xy_xz_yy, g_xx_xy_xz_yz, g_xx_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_y_xz_xx[i] = -2.0 * g_0_xy_xz_xx[i] * b_exp + 4.0 * g_xx_xy_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_xz_xy[i] = -2.0 * g_0_xy_xz_xy[i] * b_exp + 4.0 * g_xx_xy_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_xz_xz[i] = -2.0 * g_0_xy_xz_xz[i] * b_exp + 4.0 * g_xx_xy_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_xz_yy[i] = -2.0 * g_0_xy_xz_yy[i] * b_exp + 4.0 * g_xx_xy_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_xz_yz[i] = -2.0 * g_0_xy_xz_yz[i] * b_exp + 4.0 * g_xx_xy_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_xz_zz[i] = -2.0 * g_0_xy_xz_zz[i] * b_exp + 4.0 * g_xx_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_0_xy_yy_xx, g_0_xy_yy_xy, g_0_xy_yy_xz, g_0_xy_yy_yy, g_0_xy_yy_yz, g_0_xy_yy_zz, g_x_x_0_0_x_y_yy_xx, g_x_x_0_0_x_y_yy_xy, g_x_x_0_0_x_y_yy_xz, g_x_x_0_0_x_y_yy_yy, g_x_x_0_0_x_y_yy_yz, g_x_x_0_0_x_y_yy_zz, g_xx_xy_yy_xx, g_xx_xy_yy_xy, g_xx_xy_yy_xz, g_xx_xy_yy_yy, g_xx_xy_yy_yz, g_xx_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_y_yy_xx[i] = -2.0 * g_0_xy_yy_xx[i] * b_exp + 4.0 * g_xx_xy_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_yy_xy[i] = -2.0 * g_0_xy_yy_xy[i] * b_exp + 4.0 * g_xx_xy_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_yy_xz[i] = -2.0 * g_0_xy_yy_xz[i] * b_exp + 4.0 * g_xx_xy_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_yy_yy[i] = -2.0 * g_0_xy_yy_yy[i] * b_exp + 4.0 * g_xx_xy_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_yy_yz[i] = -2.0 * g_0_xy_yy_yz[i] * b_exp + 4.0 * g_xx_xy_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_yy_zz[i] = -2.0 * g_0_xy_yy_zz[i] * b_exp + 4.0 * g_xx_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_0_xy_yz_xx, g_0_xy_yz_xy, g_0_xy_yz_xz, g_0_xy_yz_yy, g_0_xy_yz_yz, g_0_xy_yz_zz, g_x_x_0_0_x_y_yz_xx, g_x_x_0_0_x_y_yz_xy, g_x_x_0_0_x_y_yz_xz, g_x_x_0_0_x_y_yz_yy, g_x_x_0_0_x_y_yz_yz, g_x_x_0_0_x_y_yz_zz, g_xx_xy_yz_xx, g_xx_xy_yz_xy, g_xx_xy_yz_xz, g_xx_xy_yz_yy, g_xx_xy_yz_yz, g_xx_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_y_yz_xx[i] = -2.0 * g_0_xy_yz_xx[i] * b_exp + 4.0 * g_xx_xy_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_yz_xy[i] = -2.0 * g_0_xy_yz_xy[i] * b_exp + 4.0 * g_xx_xy_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_yz_xz[i] = -2.0 * g_0_xy_yz_xz[i] * b_exp + 4.0 * g_xx_xy_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_yz_yy[i] = -2.0 * g_0_xy_yz_yy[i] * b_exp + 4.0 * g_xx_xy_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_yz_yz[i] = -2.0 * g_0_xy_yz_yz[i] * b_exp + 4.0 * g_xx_xy_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_yz_zz[i] = -2.0 * g_0_xy_yz_zz[i] * b_exp + 4.0 * g_xx_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_0_xy_zz_xx, g_0_xy_zz_xy, g_0_xy_zz_xz, g_0_xy_zz_yy, g_0_xy_zz_yz, g_0_xy_zz_zz, g_x_x_0_0_x_y_zz_xx, g_x_x_0_0_x_y_zz_xy, g_x_x_0_0_x_y_zz_xz, g_x_x_0_0_x_y_zz_yy, g_x_x_0_0_x_y_zz_yz, g_x_x_0_0_x_y_zz_zz, g_xx_xy_zz_xx, g_xx_xy_zz_xy, g_xx_xy_zz_xz, g_xx_xy_zz_yy, g_xx_xy_zz_yz, g_xx_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_y_zz_xx[i] = -2.0 * g_0_xy_zz_xx[i] * b_exp + 4.0 * g_xx_xy_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_zz_xy[i] = -2.0 * g_0_xy_zz_xy[i] * b_exp + 4.0 * g_xx_xy_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_zz_xz[i] = -2.0 * g_0_xy_zz_xz[i] * b_exp + 4.0 * g_xx_xy_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_zz_yy[i] = -2.0 * g_0_xy_zz_yy[i] * b_exp + 4.0 * g_xx_xy_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_zz_yz[i] = -2.0 * g_0_xy_zz_yz[i] * b_exp + 4.0 * g_xx_xy_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_zz_zz[i] = -2.0 * g_0_xy_zz_zz[i] * b_exp + 4.0 * g_xx_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_0_xz_xx_xx, g_0_xz_xx_xy, g_0_xz_xx_xz, g_0_xz_xx_yy, g_0_xz_xx_yz, g_0_xz_xx_zz, g_x_x_0_0_x_z_xx_xx, g_x_x_0_0_x_z_xx_xy, g_x_x_0_0_x_z_xx_xz, g_x_x_0_0_x_z_xx_yy, g_x_x_0_0_x_z_xx_yz, g_x_x_0_0_x_z_xx_zz, g_xx_xz_xx_xx, g_xx_xz_xx_xy, g_xx_xz_xx_xz, g_xx_xz_xx_yy, g_xx_xz_xx_yz, g_xx_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_z_xx_xx[i] = -2.0 * g_0_xz_xx_xx[i] * b_exp + 4.0 * g_xx_xz_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_xx_xy[i] = -2.0 * g_0_xz_xx_xy[i] * b_exp + 4.0 * g_xx_xz_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_xx_xz[i] = -2.0 * g_0_xz_xx_xz[i] * b_exp + 4.0 * g_xx_xz_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_xx_yy[i] = -2.0 * g_0_xz_xx_yy[i] * b_exp + 4.0 * g_xx_xz_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_xx_yz[i] = -2.0 * g_0_xz_xx_yz[i] * b_exp + 4.0 * g_xx_xz_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_xx_zz[i] = -2.0 * g_0_xz_xx_zz[i] * b_exp + 4.0 * g_xx_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_0_xz_xy_xx, g_0_xz_xy_xy, g_0_xz_xy_xz, g_0_xz_xy_yy, g_0_xz_xy_yz, g_0_xz_xy_zz, g_x_x_0_0_x_z_xy_xx, g_x_x_0_0_x_z_xy_xy, g_x_x_0_0_x_z_xy_xz, g_x_x_0_0_x_z_xy_yy, g_x_x_0_0_x_z_xy_yz, g_x_x_0_0_x_z_xy_zz, g_xx_xz_xy_xx, g_xx_xz_xy_xy, g_xx_xz_xy_xz, g_xx_xz_xy_yy, g_xx_xz_xy_yz, g_xx_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_z_xy_xx[i] = -2.0 * g_0_xz_xy_xx[i] * b_exp + 4.0 * g_xx_xz_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_xy_xy[i] = -2.0 * g_0_xz_xy_xy[i] * b_exp + 4.0 * g_xx_xz_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_xy_xz[i] = -2.0 * g_0_xz_xy_xz[i] * b_exp + 4.0 * g_xx_xz_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_xy_yy[i] = -2.0 * g_0_xz_xy_yy[i] * b_exp + 4.0 * g_xx_xz_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_xy_yz[i] = -2.0 * g_0_xz_xy_yz[i] * b_exp + 4.0 * g_xx_xz_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_xy_zz[i] = -2.0 * g_0_xz_xy_zz[i] * b_exp + 4.0 * g_xx_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_0_xz_xz_xx, g_0_xz_xz_xy, g_0_xz_xz_xz, g_0_xz_xz_yy, g_0_xz_xz_yz, g_0_xz_xz_zz, g_x_x_0_0_x_z_xz_xx, g_x_x_0_0_x_z_xz_xy, g_x_x_0_0_x_z_xz_xz, g_x_x_0_0_x_z_xz_yy, g_x_x_0_0_x_z_xz_yz, g_x_x_0_0_x_z_xz_zz, g_xx_xz_xz_xx, g_xx_xz_xz_xy, g_xx_xz_xz_xz, g_xx_xz_xz_yy, g_xx_xz_xz_yz, g_xx_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_z_xz_xx[i] = -2.0 * g_0_xz_xz_xx[i] * b_exp + 4.0 * g_xx_xz_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_xz_xy[i] = -2.0 * g_0_xz_xz_xy[i] * b_exp + 4.0 * g_xx_xz_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_xz_xz[i] = -2.0 * g_0_xz_xz_xz[i] * b_exp + 4.0 * g_xx_xz_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_xz_yy[i] = -2.0 * g_0_xz_xz_yy[i] * b_exp + 4.0 * g_xx_xz_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_xz_yz[i] = -2.0 * g_0_xz_xz_yz[i] * b_exp + 4.0 * g_xx_xz_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_xz_zz[i] = -2.0 * g_0_xz_xz_zz[i] * b_exp + 4.0 * g_xx_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_0_xz_yy_xx, g_0_xz_yy_xy, g_0_xz_yy_xz, g_0_xz_yy_yy, g_0_xz_yy_yz, g_0_xz_yy_zz, g_x_x_0_0_x_z_yy_xx, g_x_x_0_0_x_z_yy_xy, g_x_x_0_0_x_z_yy_xz, g_x_x_0_0_x_z_yy_yy, g_x_x_0_0_x_z_yy_yz, g_x_x_0_0_x_z_yy_zz, g_xx_xz_yy_xx, g_xx_xz_yy_xy, g_xx_xz_yy_xz, g_xx_xz_yy_yy, g_xx_xz_yy_yz, g_xx_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_z_yy_xx[i] = -2.0 * g_0_xz_yy_xx[i] * b_exp + 4.0 * g_xx_xz_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_yy_xy[i] = -2.0 * g_0_xz_yy_xy[i] * b_exp + 4.0 * g_xx_xz_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_yy_xz[i] = -2.0 * g_0_xz_yy_xz[i] * b_exp + 4.0 * g_xx_xz_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_yy_yy[i] = -2.0 * g_0_xz_yy_yy[i] * b_exp + 4.0 * g_xx_xz_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_yy_yz[i] = -2.0 * g_0_xz_yy_yz[i] * b_exp + 4.0 * g_xx_xz_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_yy_zz[i] = -2.0 * g_0_xz_yy_zz[i] * b_exp + 4.0 * g_xx_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_0_xz_yz_xx, g_0_xz_yz_xy, g_0_xz_yz_xz, g_0_xz_yz_yy, g_0_xz_yz_yz, g_0_xz_yz_zz, g_x_x_0_0_x_z_yz_xx, g_x_x_0_0_x_z_yz_xy, g_x_x_0_0_x_z_yz_xz, g_x_x_0_0_x_z_yz_yy, g_x_x_0_0_x_z_yz_yz, g_x_x_0_0_x_z_yz_zz, g_xx_xz_yz_xx, g_xx_xz_yz_xy, g_xx_xz_yz_xz, g_xx_xz_yz_yy, g_xx_xz_yz_yz, g_xx_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_z_yz_xx[i] = -2.0 * g_0_xz_yz_xx[i] * b_exp + 4.0 * g_xx_xz_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_yz_xy[i] = -2.0 * g_0_xz_yz_xy[i] * b_exp + 4.0 * g_xx_xz_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_yz_xz[i] = -2.0 * g_0_xz_yz_xz[i] * b_exp + 4.0 * g_xx_xz_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_yz_yy[i] = -2.0 * g_0_xz_yz_yy[i] * b_exp + 4.0 * g_xx_xz_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_yz_yz[i] = -2.0 * g_0_xz_yz_yz[i] * b_exp + 4.0 * g_xx_xz_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_yz_zz[i] = -2.0 * g_0_xz_yz_zz[i] * b_exp + 4.0 * g_xx_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_0_xz_zz_xx, g_0_xz_zz_xy, g_0_xz_zz_xz, g_0_xz_zz_yy, g_0_xz_zz_yz, g_0_xz_zz_zz, g_x_x_0_0_x_z_zz_xx, g_x_x_0_0_x_z_zz_xy, g_x_x_0_0_x_z_zz_xz, g_x_x_0_0_x_z_zz_yy, g_x_x_0_0_x_z_zz_yz, g_x_x_0_0_x_z_zz_zz, g_xx_xz_zz_xx, g_xx_xz_zz_xy, g_xx_xz_zz_xz, g_xx_xz_zz_yy, g_xx_xz_zz_yz, g_xx_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_z_zz_xx[i] = -2.0 * g_0_xz_zz_xx[i] * b_exp + 4.0 * g_xx_xz_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_zz_xy[i] = -2.0 * g_0_xz_zz_xy[i] * b_exp + 4.0 * g_xx_xz_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_zz_xz[i] = -2.0 * g_0_xz_zz_xz[i] * b_exp + 4.0 * g_xx_xz_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_zz_yy[i] = -2.0 * g_0_xz_zz_yy[i] * b_exp + 4.0 * g_xx_xz_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_zz_yz[i] = -2.0 * g_0_xz_zz_yz[i] * b_exp + 4.0 * g_xx_xz_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_zz_zz[i] = -2.0 * g_0_xz_zz_zz[i] * b_exp + 4.0 * g_xx_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_x_x_0_0_y_x_xx_xx, g_x_x_0_0_y_x_xx_xy, g_x_x_0_0_y_x_xx_xz, g_x_x_0_0_y_x_xx_yy, g_x_x_0_0_y_x_xx_yz, g_x_x_0_0_y_x_xx_zz, g_xy_0_xx_xx, g_xy_0_xx_xy, g_xy_0_xx_xz, g_xy_0_xx_yy, g_xy_0_xx_yz, g_xy_0_xx_zz, g_xy_xx_xx_xx, g_xy_xx_xx_xy, g_xy_xx_xx_xz, g_xy_xx_xx_yy, g_xy_xx_xx_yz, g_xy_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_x_xx_xx[i] = -2.0 * g_xy_0_xx_xx[i] * a_exp + 4.0 * g_xy_xx_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_xx_xy[i] = -2.0 * g_xy_0_xx_xy[i] * a_exp + 4.0 * g_xy_xx_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_xx_xz[i] = -2.0 * g_xy_0_xx_xz[i] * a_exp + 4.0 * g_xy_xx_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_xx_yy[i] = -2.0 * g_xy_0_xx_yy[i] * a_exp + 4.0 * g_xy_xx_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_xx_yz[i] = -2.0 * g_xy_0_xx_yz[i] * a_exp + 4.0 * g_xy_xx_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_xx_zz[i] = -2.0 * g_xy_0_xx_zz[i] * a_exp + 4.0 * g_xy_xx_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_x_x_0_0_y_x_xy_xx, g_x_x_0_0_y_x_xy_xy, g_x_x_0_0_y_x_xy_xz, g_x_x_0_0_y_x_xy_yy, g_x_x_0_0_y_x_xy_yz, g_x_x_0_0_y_x_xy_zz, g_xy_0_xy_xx, g_xy_0_xy_xy, g_xy_0_xy_xz, g_xy_0_xy_yy, g_xy_0_xy_yz, g_xy_0_xy_zz, g_xy_xx_xy_xx, g_xy_xx_xy_xy, g_xy_xx_xy_xz, g_xy_xx_xy_yy, g_xy_xx_xy_yz, g_xy_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_x_xy_xx[i] = -2.0 * g_xy_0_xy_xx[i] * a_exp + 4.0 * g_xy_xx_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_xy_xy[i] = -2.0 * g_xy_0_xy_xy[i] * a_exp + 4.0 * g_xy_xx_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_xy_xz[i] = -2.0 * g_xy_0_xy_xz[i] * a_exp + 4.0 * g_xy_xx_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_xy_yy[i] = -2.0 * g_xy_0_xy_yy[i] * a_exp + 4.0 * g_xy_xx_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_xy_yz[i] = -2.0 * g_xy_0_xy_yz[i] * a_exp + 4.0 * g_xy_xx_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_xy_zz[i] = -2.0 * g_xy_0_xy_zz[i] * a_exp + 4.0 * g_xy_xx_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_x_x_0_0_y_x_xz_xx, g_x_x_0_0_y_x_xz_xy, g_x_x_0_0_y_x_xz_xz, g_x_x_0_0_y_x_xz_yy, g_x_x_0_0_y_x_xz_yz, g_x_x_0_0_y_x_xz_zz, g_xy_0_xz_xx, g_xy_0_xz_xy, g_xy_0_xz_xz, g_xy_0_xz_yy, g_xy_0_xz_yz, g_xy_0_xz_zz, g_xy_xx_xz_xx, g_xy_xx_xz_xy, g_xy_xx_xz_xz, g_xy_xx_xz_yy, g_xy_xx_xz_yz, g_xy_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_x_xz_xx[i] = -2.0 * g_xy_0_xz_xx[i] * a_exp + 4.0 * g_xy_xx_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_xz_xy[i] = -2.0 * g_xy_0_xz_xy[i] * a_exp + 4.0 * g_xy_xx_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_xz_xz[i] = -2.0 * g_xy_0_xz_xz[i] * a_exp + 4.0 * g_xy_xx_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_xz_yy[i] = -2.0 * g_xy_0_xz_yy[i] * a_exp + 4.0 * g_xy_xx_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_xz_yz[i] = -2.0 * g_xy_0_xz_yz[i] * a_exp + 4.0 * g_xy_xx_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_xz_zz[i] = -2.0 * g_xy_0_xz_zz[i] * a_exp + 4.0 * g_xy_xx_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_x_0_0_y_x_yy_xx, g_x_x_0_0_y_x_yy_xy, g_x_x_0_0_y_x_yy_xz, g_x_x_0_0_y_x_yy_yy, g_x_x_0_0_y_x_yy_yz, g_x_x_0_0_y_x_yy_zz, g_xy_0_yy_xx, g_xy_0_yy_xy, g_xy_0_yy_xz, g_xy_0_yy_yy, g_xy_0_yy_yz, g_xy_0_yy_zz, g_xy_xx_yy_xx, g_xy_xx_yy_xy, g_xy_xx_yy_xz, g_xy_xx_yy_yy, g_xy_xx_yy_yz, g_xy_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_x_yy_xx[i] = -2.0 * g_xy_0_yy_xx[i] * a_exp + 4.0 * g_xy_xx_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_yy_xy[i] = -2.0 * g_xy_0_yy_xy[i] * a_exp + 4.0 * g_xy_xx_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_yy_xz[i] = -2.0 * g_xy_0_yy_xz[i] * a_exp + 4.0 * g_xy_xx_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_yy_yy[i] = -2.0 * g_xy_0_yy_yy[i] * a_exp + 4.0 * g_xy_xx_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_yy_yz[i] = -2.0 * g_xy_0_yy_yz[i] * a_exp + 4.0 * g_xy_xx_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_yy_zz[i] = -2.0 * g_xy_0_yy_zz[i] * a_exp + 4.0 * g_xy_xx_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_x_0_0_y_x_yz_xx, g_x_x_0_0_y_x_yz_xy, g_x_x_0_0_y_x_yz_xz, g_x_x_0_0_y_x_yz_yy, g_x_x_0_0_y_x_yz_yz, g_x_x_0_0_y_x_yz_zz, g_xy_0_yz_xx, g_xy_0_yz_xy, g_xy_0_yz_xz, g_xy_0_yz_yy, g_xy_0_yz_yz, g_xy_0_yz_zz, g_xy_xx_yz_xx, g_xy_xx_yz_xy, g_xy_xx_yz_xz, g_xy_xx_yz_yy, g_xy_xx_yz_yz, g_xy_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_x_yz_xx[i] = -2.0 * g_xy_0_yz_xx[i] * a_exp + 4.0 * g_xy_xx_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_yz_xy[i] = -2.0 * g_xy_0_yz_xy[i] * a_exp + 4.0 * g_xy_xx_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_yz_xz[i] = -2.0 * g_xy_0_yz_xz[i] * a_exp + 4.0 * g_xy_xx_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_yz_yy[i] = -2.0 * g_xy_0_yz_yy[i] * a_exp + 4.0 * g_xy_xx_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_yz_yz[i] = -2.0 * g_xy_0_yz_yz[i] * a_exp + 4.0 * g_xy_xx_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_yz_zz[i] = -2.0 * g_xy_0_yz_zz[i] * a_exp + 4.0 * g_xy_xx_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_x_0_0_y_x_zz_xx, g_x_x_0_0_y_x_zz_xy, g_x_x_0_0_y_x_zz_xz, g_x_x_0_0_y_x_zz_yy, g_x_x_0_0_y_x_zz_yz, g_x_x_0_0_y_x_zz_zz, g_xy_0_zz_xx, g_xy_0_zz_xy, g_xy_0_zz_xz, g_xy_0_zz_yy, g_xy_0_zz_yz, g_xy_0_zz_zz, g_xy_xx_zz_xx, g_xy_xx_zz_xy, g_xy_xx_zz_xz, g_xy_xx_zz_yy, g_xy_xx_zz_yz, g_xy_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_x_zz_xx[i] = -2.0 * g_xy_0_zz_xx[i] * a_exp + 4.0 * g_xy_xx_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_zz_xy[i] = -2.0 * g_xy_0_zz_xy[i] * a_exp + 4.0 * g_xy_xx_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_zz_xz[i] = -2.0 * g_xy_0_zz_xz[i] * a_exp + 4.0 * g_xy_xx_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_zz_yy[i] = -2.0 * g_xy_0_zz_yy[i] * a_exp + 4.0 * g_xy_xx_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_zz_yz[i] = -2.0 * g_xy_0_zz_yz[i] * a_exp + 4.0 * g_xy_xx_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_zz_zz[i] = -2.0 * g_xy_0_zz_zz[i] * a_exp + 4.0 * g_xy_xx_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_x_0_0_y_y_xx_xx, g_x_x_0_0_y_y_xx_xy, g_x_x_0_0_y_y_xx_xz, g_x_x_0_0_y_y_xx_yy, g_x_x_0_0_y_y_xx_yz, g_x_x_0_0_y_y_xx_zz, g_xy_xy_xx_xx, g_xy_xy_xx_xy, g_xy_xy_xx_xz, g_xy_xy_xx_yy, g_xy_xy_xx_yz, g_xy_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_y_xx_xx[i] = 4.0 * g_xy_xy_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_xx_xy[i] = 4.0 * g_xy_xy_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_xx_xz[i] = 4.0 * g_xy_xy_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_xx_yy[i] = 4.0 * g_xy_xy_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_xx_yz[i] = 4.0 * g_xy_xy_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_xx_zz[i] = 4.0 * g_xy_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_x_0_0_y_y_xy_xx, g_x_x_0_0_y_y_xy_xy, g_x_x_0_0_y_y_xy_xz, g_x_x_0_0_y_y_xy_yy, g_x_x_0_0_y_y_xy_yz, g_x_x_0_0_y_y_xy_zz, g_xy_xy_xy_xx, g_xy_xy_xy_xy, g_xy_xy_xy_xz, g_xy_xy_xy_yy, g_xy_xy_xy_yz, g_xy_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_y_xy_xx[i] = 4.0 * g_xy_xy_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_xy_xy[i] = 4.0 * g_xy_xy_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_xy_xz[i] = 4.0 * g_xy_xy_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_xy_yy[i] = 4.0 * g_xy_xy_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_xy_yz[i] = 4.0 * g_xy_xy_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_xy_zz[i] = 4.0 * g_xy_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_x_0_0_y_y_xz_xx, g_x_x_0_0_y_y_xz_xy, g_x_x_0_0_y_y_xz_xz, g_x_x_0_0_y_y_xz_yy, g_x_x_0_0_y_y_xz_yz, g_x_x_0_0_y_y_xz_zz, g_xy_xy_xz_xx, g_xy_xy_xz_xy, g_xy_xy_xz_xz, g_xy_xy_xz_yy, g_xy_xy_xz_yz, g_xy_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_y_xz_xx[i] = 4.0 * g_xy_xy_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_xz_xy[i] = 4.0 * g_xy_xy_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_xz_xz[i] = 4.0 * g_xy_xy_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_xz_yy[i] = 4.0 * g_xy_xy_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_xz_yz[i] = 4.0 * g_xy_xy_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_xz_zz[i] = 4.0 * g_xy_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_x_x_0_0_y_y_yy_xx, g_x_x_0_0_y_y_yy_xy, g_x_x_0_0_y_y_yy_xz, g_x_x_0_0_y_y_yy_yy, g_x_x_0_0_y_y_yy_yz, g_x_x_0_0_y_y_yy_zz, g_xy_xy_yy_xx, g_xy_xy_yy_xy, g_xy_xy_yy_xz, g_xy_xy_yy_yy, g_xy_xy_yy_yz, g_xy_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_y_yy_xx[i] = 4.0 * g_xy_xy_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_yy_xy[i] = 4.0 * g_xy_xy_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_yy_xz[i] = 4.0 * g_xy_xy_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_yy_yy[i] = 4.0 * g_xy_xy_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_yy_yz[i] = 4.0 * g_xy_xy_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_yy_zz[i] = 4.0 * g_xy_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_x_x_0_0_y_y_yz_xx, g_x_x_0_0_y_y_yz_xy, g_x_x_0_0_y_y_yz_xz, g_x_x_0_0_y_y_yz_yy, g_x_x_0_0_y_y_yz_yz, g_x_x_0_0_y_y_yz_zz, g_xy_xy_yz_xx, g_xy_xy_yz_xy, g_xy_xy_yz_xz, g_xy_xy_yz_yy, g_xy_xy_yz_yz, g_xy_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_y_yz_xx[i] = 4.0 * g_xy_xy_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_yz_xy[i] = 4.0 * g_xy_xy_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_yz_xz[i] = 4.0 * g_xy_xy_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_yz_yy[i] = 4.0 * g_xy_xy_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_yz_yz[i] = 4.0 * g_xy_xy_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_yz_zz[i] = 4.0 * g_xy_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_x_x_0_0_y_y_zz_xx, g_x_x_0_0_y_y_zz_xy, g_x_x_0_0_y_y_zz_xz, g_x_x_0_0_y_y_zz_yy, g_x_x_0_0_y_y_zz_yz, g_x_x_0_0_y_y_zz_zz, g_xy_xy_zz_xx, g_xy_xy_zz_xy, g_xy_xy_zz_xz, g_xy_xy_zz_yy, g_xy_xy_zz_yz, g_xy_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_y_zz_xx[i] = 4.0 * g_xy_xy_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_zz_xy[i] = 4.0 * g_xy_xy_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_zz_xz[i] = 4.0 * g_xy_xy_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_zz_yy[i] = 4.0 * g_xy_xy_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_zz_yz[i] = 4.0 * g_xy_xy_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_zz_zz[i] = 4.0 * g_xy_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_x_x_0_0_y_z_xx_xx, g_x_x_0_0_y_z_xx_xy, g_x_x_0_0_y_z_xx_xz, g_x_x_0_0_y_z_xx_yy, g_x_x_0_0_y_z_xx_yz, g_x_x_0_0_y_z_xx_zz, g_xy_xz_xx_xx, g_xy_xz_xx_xy, g_xy_xz_xx_xz, g_xy_xz_xx_yy, g_xy_xz_xx_yz, g_xy_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_z_xx_xx[i] = 4.0 * g_xy_xz_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_xx_xy[i] = 4.0 * g_xy_xz_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_xx_xz[i] = 4.0 * g_xy_xz_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_xx_yy[i] = 4.0 * g_xy_xz_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_xx_yz[i] = 4.0 * g_xy_xz_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_xx_zz[i] = 4.0 * g_xy_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_x_x_0_0_y_z_xy_xx, g_x_x_0_0_y_z_xy_xy, g_x_x_0_0_y_z_xy_xz, g_x_x_0_0_y_z_xy_yy, g_x_x_0_0_y_z_xy_yz, g_x_x_0_0_y_z_xy_zz, g_xy_xz_xy_xx, g_xy_xz_xy_xy, g_xy_xz_xy_xz, g_xy_xz_xy_yy, g_xy_xz_xy_yz, g_xy_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_z_xy_xx[i] = 4.0 * g_xy_xz_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_xy_xy[i] = 4.0 * g_xy_xz_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_xy_xz[i] = 4.0 * g_xy_xz_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_xy_yy[i] = 4.0 * g_xy_xz_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_xy_yz[i] = 4.0 * g_xy_xz_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_xy_zz[i] = 4.0 * g_xy_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_x_x_0_0_y_z_xz_xx, g_x_x_0_0_y_z_xz_xy, g_x_x_0_0_y_z_xz_xz, g_x_x_0_0_y_z_xz_yy, g_x_x_0_0_y_z_xz_yz, g_x_x_0_0_y_z_xz_zz, g_xy_xz_xz_xx, g_xy_xz_xz_xy, g_xy_xz_xz_xz, g_xy_xz_xz_yy, g_xy_xz_xz_yz, g_xy_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_z_xz_xx[i] = 4.0 * g_xy_xz_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_xz_xy[i] = 4.0 * g_xy_xz_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_xz_xz[i] = 4.0 * g_xy_xz_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_xz_yy[i] = 4.0 * g_xy_xz_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_xz_yz[i] = 4.0 * g_xy_xz_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_xz_zz[i] = 4.0 * g_xy_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_x_x_0_0_y_z_yy_xx, g_x_x_0_0_y_z_yy_xy, g_x_x_0_0_y_z_yy_xz, g_x_x_0_0_y_z_yy_yy, g_x_x_0_0_y_z_yy_yz, g_x_x_0_0_y_z_yy_zz, g_xy_xz_yy_xx, g_xy_xz_yy_xy, g_xy_xz_yy_xz, g_xy_xz_yy_yy, g_xy_xz_yy_yz, g_xy_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_z_yy_xx[i] = 4.0 * g_xy_xz_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_yy_xy[i] = 4.0 * g_xy_xz_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_yy_xz[i] = 4.0 * g_xy_xz_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_yy_yy[i] = 4.0 * g_xy_xz_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_yy_yz[i] = 4.0 * g_xy_xz_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_yy_zz[i] = 4.0 * g_xy_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_x_x_0_0_y_z_yz_xx, g_x_x_0_0_y_z_yz_xy, g_x_x_0_0_y_z_yz_xz, g_x_x_0_0_y_z_yz_yy, g_x_x_0_0_y_z_yz_yz, g_x_x_0_0_y_z_yz_zz, g_xy_xz_yz_xx, g_xy_xz_yz_xy, g_xy_xz_yz_xz, g_xy_xz_yz_yy, g_xy_xz_yz_yz, g_xy_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_z_yz_xx[i] = 4.0 * g_xy_xz_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_yz_xy[i] = 4.0 * g_xy_xz_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_yz_xz[i] = 4.0 * g_xy_xz_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_yz_yy[i] = 4.0 * g_xy_xz_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_yz_yz[i] = 4.0 * g_xy_xz_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_yz_zz[i] = 4.0 * g_xy_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_x_x_0_0_y_z_zz_xx, g_x_x_0_0_y_z_zz_xy, g_x_x_0_0_y_z_zz_xz, g_x_x_0_0_y_z_zz_yy, g_x_x_0_0_y_z_zz_yz, g_x_x_0_0_y_z_zz_zz, g_xy_xz_zz_xx, g_xy_xz_zz_xy, g_xy_xz_zz_xz, g_xy_xz_zz_yy, g_xy_xz_zz_yz, g_xy_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_z_zz_xx[i] = 4.0 * g_xy_xz_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_zz_xy[i] = 4.0 * g_xy_xz_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_zz_xz[i] = 4.0 * g_xy_xz_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_zz_yy[i] = 4.0 * g_xy_xz_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_zz_yz[i] = 4.0 * g_xy_xz_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_zz_zz[i] = 4.0 * g_xy_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_x_x_0_0_z_x_xx_xx, g_x_x_0_0_z_x_xx_xy, g_x_x_0_0_z_x_xx_xz, g_x_x_0_0_z_x_xx_yy, g_x_x_0_0_z_x_xx_yz, g_x_x_0_0_z_x_xx_zz, g_xz_0_xx_xx, g_xz_0_xx_xy, g_xz_0_xx_xz, g_xz_0_xx_yy, g_xz_0_xx_yz, g_xz_0_xx_zz, g_xz_xx_xx_xx, g_xz_xx_xx_xy, g_xz_xx_xx_xz, g_xz_xx_xx_yy, g_xz_xx_xx_yz, g_xz_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_x_xx_xx[i] = -2.0 * g_xz_0_xx_xx[i] * a_exp + 4.0 * g_xz_xx_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_xx_xy[i] = -2.0 * g_xz_0_xx_xy[i] * a_exp + 4.0 * g_xz_xx_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_xx_xz[i] = -2.0 * g_xz_0_xx_xz[i] * a_exp + 4.0 * g_xz_xx_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_xx_yy[i] = -2.0 * g_xz_0_xx_yy[i] * a_exp + 4.0 * g_xz_xx_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_xx_yz[i] = -2.0 * g_xz_0_xx_yz[i] * a_exp + 4.0 * g_xz_xx_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_xx_zz[i] = -2.0 * g_xz_0_xx_zz[i] * a_exp + 4.0 * g_xz_xx_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_x_x_0_0_z_x_xy_xx, g_x_x_0_0_z_x_xy_xy, g_x_x_0_0_z_x_xy_xz, g_x_x_0_0_z_x_xy_yy, g_x_x_0_0_z_x_xy_yz, g_x_x_0_0_z_x_xy_zz, g_xz_0_xy_xx, g_xz_0_xy_xy, g_xz_0_xy_xz, g_xz_0_xy_yy, g_xz_0_xy_yz, g_xz_0_xy_zz, g_xz_xx_xy_xx, g_xz_xx_xy_xy, g_xz_xx_xy_xz, g_xz_xx_xy_yy, g_xz_xx_xy_yz, g_xz_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_x_xy_xx[i] = -2.0 * g_xz_0_xy_xx[i] * a_exp + 4.0 * g_xz_xx_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_xy_xy[i] = -2.0 * g_xz_0_xy_xy[i] * a_exp + 4.0 * g_xz_xx_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_xy_xz[i] = -2.0 * g_xz_0_xy_xz[i] * a_exp + 4.0 * g_xz_xx_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_xy_yy[i] = -2.0 * g_xz_0_xy_yy[i] * a_exp + 4.0 * g_xz_xx_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_xy_yz[i] = -2.0 * g_xz_0_xy_yz[i] * a_exp + 4.0 * g_xz_xx_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_xy_zz[i] = -2.0 * g_xz_0_xy_zz[i] * a_exp + 4.0 * g_xz_xx_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_x_x_0_0_z_x_xz_xx, g_x_x_0_0_z_x_xz_xy, g_x_x_0_0_z_x_xz_xz, g_x_x_0_0_z_x_xz_yy, g_x_x_0_0_z_x_xz_yz, g_x_x_0_0_z_x_xz_zz, g_xz_0_xz_xx, g_xz_0_xz_xy, g_xz_0_xz_xz, g_xz_0_xz_yy, g_xz_0_xz_yz, g_xz_0_xz_zz, g_xz_xx_xz_xx, g_xz_xx_xz_xy, g_xz_xx_xz_xz, g_xz_xx_xz_yy, g_xz_xx_xz_yz, g_xz_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_x_xz_xx[i] = -2.0 * g_xz_0_xz_xx[i] * a_exp + 4.0 * g_xz_xx_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_xz_xy[i] = -2.0 * g_xz_0_xz_xy[i] * a_exp + 4.0 * g_xz_xx_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_xz_xz[i] = -2.0 * g_xz_0_xz_xz[i] * a_exp + 4.0 * g_xz_xx_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_xz_yy[i] = -2.0 * g_xz_0_xz_yy[i] * a_exp + 4.0 * g_xz_xx_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_xz_yz[i] = -2.0 * g_xz_0_xz_yz[i] * a_exp + 4.0 * g_xz_xx_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_xz_zz[i] = -2.0 * g_xz_0_xz_zz[i] * a_exp + 4.0 * g_xz_xx_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_x_x_0_0_z_x_yy_xx, g_x_x_0_0_z_x_yy_xy, g_x_x_0_0_z_x_yy_xz, g_x_x_0_0_z_x_yy_yy, g_x_x_0_0_z_x_yy_yz, g_x_x_0_0_z_x_yy_zz, g_xz_0_yy_xx, g_xz_0_yy_xy, g_xz_0_yy_xz, g_xz_0_yy_yy, g_xz_0_yy_yz, g_xz_0_yy_zz, g_xz_xx_yy_xx, g_xz_xx_yy_xy, g_xz_xx_yy_xz, g_xz_xx_yy_yy, g_xz_xx_yy_yz, g_xz_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_x_yy_xx[i] = -2.0 * g_xz_0_yy_xx[i] * a_exp + 4.0 * g_xz_xx_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_yy_xy[i] = -2.0 * g_xz_0_yy_xy[i] * a_exp + 4.0 * g_xz_xx_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_yy_xz[i] = -2.0 * g_xz_0_yy_xz[i] * a_exp + 4.0 * g_xz_xx_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_yy_yy[i] = -2.0 * g_xz_0_yy_yy[i] * a_exp + 4.0 * g_xz_xx_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_yy_yz[i] = -2.0 * g_xz_0_yy_yz[i] * a_exp + 4.0 * g_xz_xx_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_yy_zz[i] = -2.0 * g_xz_0_yy_zz[i] * a_exp + 4.0 * g_xz_xx_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_x_x_0_0_z_x_yz_xx, g_x_x_0_0_z_x_yz_xy, g_x_x_0_0_z_x_yz_xz, g_x_x_0_0_z_x_yz_yy, g_x_x_0_0_z_x_yz_yz, g_x_x_0_0_z_x_yz_zz, g_xz_0_yz_xx, g_xz_0_yz_xy, g_xz_0_yz_xz, g_xz_0_yz_yy, g_xz_0_yz_yz, g_xz_0_yz_zz, g_xz_xx_yz_xx, g_xz_xx_yz_xy, g_xz_xx_yz_xz, g_xz_xx_yz_yy, g_xz_xx_yz_yz, g_xz_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_x_yz_xx[i] = -2.0 * g_xz_0_yz_xx[i] * a_exp + 4.0 * g_xz_xx_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_yz_xy[i] = -2.0 * g_xz_0_yz_xy[i] * a_exp + 4.0 * g_xz_xx_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_yz_xz[i] = -2.0 * g_xz_0_yz_xz[i] * a_exp + 4.0 * g_xz_xx_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_yz_yy[i] = -2.0 * g_xz_0_yz_yy[i] * a_exp + 4.0 * g_xz_xx_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_yz_yz[i] = -2.0 * g_xz_0_yz_yz[i] * a_exp + 4.0 * g_xz_xx_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_yz_zz[i] = -2.0 * g_xz_0_yz_zz[i] * a_exp + 4.0 * g_xz_xx_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_x_x_0_0_z_x_zz_xx, g_x_x_0_0_z_x_zz_xy, g_x_x_0_0_z_x_zz_xz, g_x_x_0_0_z_x_zz_yy, g_x_x_0_0_z_x_zz_yz, g_x_x_0_0_z_x_zz_zz, g_xz_0_zz_xx, g_xz_0_zz_xy, g_xz_0_zz_xz, g_xz_0_zz_yy, g_xz_0_zz_yz, g_xz_0_zz_zz, g_xz_xx_zz_xx, g_xz_xx_zz_xy, g_xz_xx_zz_xz, g_xz_xx_zz_yy, g_xz_xx_zz_yz, g_xz_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_x_zz_xx[i] = -2.0 * g_xz_0_zz_xx[i] * a_exp + 4.0 * g_xz_xx_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_zz_xy[i] = -2.0 * g_xz_0_zz_xy[i] * a_exp + 4.0 * g_xz_xx_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_zz_xz[i] = -2.0 * g_xz_0_zz_xz[i] * a_exp + 4.0 * g_xz_xx_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_zz_yy[i] = -2.0 * g_xz_0_zz_yy[i] * a_exp + 4.0 * g_xz_xx_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_zz_yz[i] = -2.0 * g_xz_0_zz_yz[i] * a_exp + 4.0 * g_xz_xx_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_zz_zz[i] = -2.0 * g_xz_0_zz_zz[i] * a_exp + 4.0 * g_xz_xx_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_x_x_0_0_z_y_xx_xx, g_x_x_0_0_z_y_xx_xy, g_x_x_0_0_z_y_xx_xz, g_x_x_0_0_z_y_xx_yy, g_x_x_0_0_z_y_xx_yz, g_x_x_0_0_z_y_xx_zz, g_xz_xy_xx_xx, g_xz_xy_xx_xy, g_xz_xy_xx_xz, g_xz_xy_xx_yy, g_xz_xy_xx_yz, g_xz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_y_xx_xx[i] = 4.0 * g_xz_xy_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_xx_xy[i] = 4.0 * g_xz_xy_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_xx_xz[i] = 4.0 * g_xz_xy_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_xx_yy[i] = 4.0 * g_xz_xy_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_xx_yz[i] = 4.0 * g_xz_xy_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_xx_zz[i] = 4.0 * g_xz_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_x_x_0_0_z_y_xy_xx, g_x_x_0_0_z_y_xy_xy, g_x_x_0_0_z_y_xy_xz, g_x_x_0_0_z_y_xy_yy, g_x_x_0_0_z_y_xy_yz, g_x_x_0_0_z_y_xy_zz, g_xz_xy_xy_xx, g_xz_xy_xy_xy, g_xz_xy_xy_xz, g_xz_xy_xy_yy, g_xz_xy_xy_yz, g_xz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_y_xy_xx[i] = 4.0 * g_xz_xy_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_xy_xy[i] = 4.0 * g_xz_xy_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_xy_xz[i] = 4.0 * g_xz_xy_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_xy_yy[i] = 4.0 * g_xz_xy_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_xy_yz[i] = 4.0 * g_xz_xy_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_xy_zz[i] = 4.0 * g_xz_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_x_x_0_0_z_y_xz_xx, g_x_x_0_0_z_y_xz_xy, g_x_x_0_0_z_y_xz_xz, g_x_x_0_0_z_y_xz_yy, g_x_x_0_0_z_y_xz_yz, g_x_x_0_0_z_y_xz_zz, g_xz_xy_xz_xx, g_xz_xy_xz_xy, g_xz_xy_xz_xz, g_xz_xy_xz_yy, g_xz_xy_xz_yz, g_xz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_y_xz_xx[i] = 4.0 * g_xz_xy_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_xz_xy[i] = 4.0 * g_xz_xy_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_xz_xz[i] = 4.0 * g_xz_xy_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_xz_yy[i] = 4.0 * g_xz_xy_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_xz_yz[i] = 4.0 * g_xz_xy_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_xz_zz[i] = 4.0 * g_xz_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_x_x_0_0_z_y_yy_xx, g_x_x_0_0_z_y_yy_xy, g_x_x_0_0_z_y_yy_xz, g_x_x_0_0_z_y_yy_yy, g_x_x_0_0_z_y_yy_yz, g_x_x_0_0_z_y_yy_zz, g_xz_xy_yy_xx, g_xz_xy_yy_xy, g_xz_xy_yy_xz, g_xz_xy_yy_yy, g_xz_xy_yy_yz, g_xz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_y_yy_xx[i] = 4.0 * g_xz_xy_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_yy_xy[i] = 4.0 * g_xz_xy_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_yy_xz[i] = 4.0 * g_xz_xy_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_yy_yy[i] = 4.0 * g_xz_xy_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_yy_yz[i] = 4.0 * g_xz_xy_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_yy_zz[i] = 4.0 * g_xz_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_x_x_0_0_z_y_yz_xx, g_x_x_0_0_z_y_yz_xy, g_x_x_0_0_z_y_yz_xz, g_x_x_0_0_z_y_yz_yy, g_x_x_0_0_z_y_yz_yz, g_x_x_0_0_z_y_yz_zz, g_xz_xy_yz_xx, g_xz_xy_yz_xy, g_xz_xy_yz_xz, g_xz_xy_yz_yy, g_xz_xy_yz_yz, g_xz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_y_yz_xx[i] = 4.0 * g_xz_xy_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_yz_xy[i] = 4.0 * g_xz_xy_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_yz_xz[i] = 4.0 * g_xz_xy_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_yz_yy[i] = 4.0 * g_xz_xy_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_yz_yz[i] = 4.0 * g_xz_xy_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_yz_zz[i] = 4.0 * g_xz_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_x_x_0_0_z_y_zz_xx, g_x_x_0_0_z_y_zz_xy, g_x_x_0_0_z_y_zz_xz, g_x_x_0_0_z_y_zz_yy, g_x_x_0_0_z_y_zz_yz, g_x_x_0_0_z_y_zz_zz, g_xz_xy_zz_xx, g_xz_xy_zz_xy, g_xz_xy_zz_xz, g_xz_xy_zz_yy, g_xz_xy_zz_yz, g_xz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_y_zz_xx[i] = 4.0 * g_xz_xy_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_zz_xy[i] = 4.0 * g_xz_xy_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_zz_xz[i] = 4.0 * g_xz_xy_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_zz_yy[i] = 4.0 * g_xz_xy_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_zz_yz[i] = 4.0 * g_xz_xy_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_zz_zz[i] = 4.0 * g_xz_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_x_x_0_0_z_z_xx_xx, g_x_x_0_0_z_z_xx_xy, g_x_x_0_0_z_z_xx_xz, g_x_x_0_0_z_z_xx_yy, g_x_x_0_0_z_z_xx_yz, g_x_x_0_0_z_z_xx_zz, g_xz_xz_xx_xx, g_xz_xz_xx_xy, g_xz_xz_xx_xz, g_xz_xz_xx_yy, g_xz_xz_xx_yz, g_xz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_z_xx_xx[i] = 4.0 * g_xz_xz_xx_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_xx_xy[i] = 4.0 * g_xz_xz_xx_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_xx_xz[i] = 4.0 * g_xz_xz_xx_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_xx_yy[i] = 4.0 * g_xz_xz_xx_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_xx_yz[i] = 4.0 * g_xz_xz_xx_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_xx_zz[i] = 4.0 * g_xz_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_x_x_0_0_z_z_xy_xx, g_x_x_0_0_z_z_xy_xy, g_x_x_0_0_z_z_xy_xz, g_x_x_0_0_z_z_xy_yy, g_x_x_0_0_z_z_xy_yz, g_x_x_0_0_z_z_xy_zz, g_xz_xz_xy_xx, g_xz_xz_xy_xy, g_xz_xz_xy_xz, g_xz_xz_xy_yy, g_xz_xz_xy_yz, g_xz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_z_xy_xx[i] = 4.0 * g_xz_xz_xy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_xy_xy[i] = 4.0 * g_xz_xz_xy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_xy_xz[i] = 4.0 * g_xz_xz_xy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_xy_yy[i] = 4.0 * g_xz_xz_xy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_xy_yz[i] = 4.0 * g_xz_xz_xy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_xy_zz[i] = 4.0 * g_xz_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_x_x_0_0_z_z_xz_xx, g_x_x_0_0_z_z_xz_xy, g_x_x_0_0_z_z_xz_xz, g_x_x_0_0_z_z_xz_yy, g_x_x_0_0_z_z_xz_yz, g_x_x_0_0_z_z_xz_zz, g_xz_xz_xz_xx, g_xz_xz_xz_xy, g_xz_xz_xz_xz, g_xz_xz_xz_yy, g_xz_xz_xz_yz, g_xz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_z_xz_xx[i] = 4.0 * g_xz_xz_xz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_xz_xy[i] = 4.0 * g_xz_xz_xz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_xz_xz[i] = 4.0 * g_xz_xz_xz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_xz_yy[i] = 4.0 * g_xz_xz_xz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_xz_yz[i] = 4.0 * g_xz_xz_xz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_xz_zz[i] = 4.0 * g_xz_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_x_x_0_0_z_z_yy_xx, g_x_x_0_0_z_z_yy_xy, g_x_x_0_0_z_z_yy_xz, g_x_x_0_0_z_z_yy_yy, g_x_x_0_0_z_z_yy_yz, g_x_x_0_0_z_z_yy_zz, g_xz_xz_yy_xx, g_xz_xz_yy_xy, g_xz_xz_yy_xz, g_xz_xz_yy_yy, g_xz_xz_yy_yz, g_xz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_z_yy_xx[i] = 4.0 * g_xz_xz_yy_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_yy_xy[i] = 4.0 * g_xz_xz_yy_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_yy_xz[i] = 4.0 * g_xz_xz_yy_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_yy_yy[i] = 4.0 * g_xz_xz_yy_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_yy_yz[i] = 4.0 * g_xz_xz_yy_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_yy_zz[i] = 4.0 * g_xz_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_x_x_0_0_z_z_yz_xx, g_x_x_0_0_z_z_yz_xy, g_x_x_0_0_z_z_yz_xz, g_x_x_0_0_z_z_yz_yy, g_x_x_0_0_z_z_yz_yz, g_x_x_0_0_z_z_yz_zz, g_xz_xz_yz_xx, g_xz_xz_yz_xy, g_xz_xz_yz_xz, g_xz_xz_yz_yy, g_xz_xz_yz_yz, g_xz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_z_yz_xx[i] = 4.0 * g_xz_xz_yz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_yz_xy[i] = 4.0 * g_xz_xz_yz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_yz_xz[i] = 4.0 * g_xz_xz_yz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_yz_yy[i] = 4.0 * g_xz_xz_yz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_yz_yz[i] = 4.0 * g_xz_xz_yz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_yz_zz[i] = 4.0 * g_xz_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_x_x_0_0_z_z_zz_xx, g_x_x_0_0_z_z_zz_xy, g_x_x_0_0_z_z_zz_xz, g_x_x_0_0_z_z_zz_yy, g_x_x_0_0_z_z_zz_yz, g_x_x_0_0_z_z_zz_zz, g_xz_xz_zz_xx, g_xz_xz_zz_xy, g_xz_xz_zz_xz, g_xz_xz_zz_yy, g_xz_xz_zz_yz, g_xz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_z_zz_xx[i] = 4.0 * g_xz_xz_zz_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_zz_xy[i] = 4.0 * g_xz_xz_zz_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_zz_xz[i] = 4.0 * g_xz_xz_zz_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_zz_yy[i] = 4.0 * g_xz_xz_zz_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_zz_yz[i] = 4.0 * g_xz_xz_zz_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_zz_zz[i] = 4.0 * g_xz_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_0_xy_xx_xx, g_0_xy_xx_xy, g_0_xy_xx_xz, g_0_xy_xx_yy, g_0_xy_xx_yz, g_0_xy_xx_zz, g_x_y_0_0_x_x_xx_xx, g_x_y_0_0_x_x_xx_xy, g_x_y_0_0_x_x_xx_xz, g_x_y_0_0_x_x_xx_yy, g_x_y_0_0_x_x_xx_yz, g_x_y_0_0_x_x_xx_zz, g_xx_xy_xx_xx, g_xx_xy_xx_xy, g_xx_xy_xx_xz, g_xx_xy_xx_yy, g_xx_xy_xx_yz, g_xx_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_x_xx_xx[i] = -2.0 * g_0_xy_xx_xx[i] * b_exp + 4.0 * g_xx_xy_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_xx_xy[i] = -2.0 * g_0_xy_xx_xy[i] * b_exp + 4.0 * g_xx_xy_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_xx_xz[i] = -2.0 * g_0_xy_xx_xz[i] * b_exp + 4.0 * g_xx_xy_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_xx_yy[i] = -2.0 * g_0_xy_xx_yy[i] * b_exp + 4.0 * g_xx_xy_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_xx_yz[i] = -2.0 * g_0_xy_xx_yz[i] * b_exp + 4.0 * g_xx_xy_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_xx_zz[i] = -2.0 * g_0_xy_xx_zz[i] * b_exp + 4.0 * g_xx_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_0_xy_xy_xx, g_0_xy_xy_xy, g_0_xy_xy_xz, g_0_xy_xy_yy, g_0_xy_xy_yz, g_0_xy_xy_zz, g_x_y_0_0_x_x_xy_xx, g_x_y_0_0_x_x_xy_xy, g_x_y_0_0_x_x_xy_xz, g_x_y_0_0_x_x_xy_yy, g_x_y_0_0_x_x_xy_yz, g_x_y_0_0_x_x_xy_zz, g_xx_xy_xy_xx, g_xx_xy_xy_xy, g_xx_xy_xy_xz, g_xx_xy_xy_yy, g_xx_xy_xy_yz, g_xx_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_x_xy_xx[i] = -2.0 * g_0_xy_xy_xx[i] * b_exp + 4.0 * g_xx_xy_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_xy_xy[i] = -2.0 * g_0_xy_xy_xy[i] * b_exp + 4.0 * g_xx_xy_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_xy_xz[i] = -2.0 * g_0_xy_xy_xz[i] * b_exp + 4.0 * g_xx_xy_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_xy_yy[i] = -2.0 * g_0_xy_xy_yy[i] * b_exp + 4.0 * g_xx_xy_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_xy_yz[i] = -2.0 * g_0_xy_xy_yz[i] * b_exp + 4.0 * g_xx_xy_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_xy_zz[i] = -2.0 * g_0_xy_xy_zz[i] * b_exp + 4.0 * g_xx_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_0_xy_xz_xx, g_0_xy_xz_xy, g_0_xy_xz_xz, g_0_xy_xz_yy, g_0_xy_xz_yz, g_0_xy_xz_zz, g_x_y_0_0_x_x_xz_xx, g_x_y_0_0_x_x_xz_xy, g_x_y_0_0_x_x_xz_xz, g_x_y_0_0_x_x_xz_yy, g_x_y_0_0_x_x_xz_yz, g_x_y_0_0_x_x_xz_zz, g_xx_xy_xz_xx, g_xx_xy_xz_xy, g_xx_xy_xz_xz, g_xx_xy_xz_yy, g_xx_xy_xz_yz, g_xx_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_x_xz_xx[i] = -2.0 * g_0_xy_xz_xx[i] * b_exp + 4.0 * g_xx_xy_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_xz_xy[i] = -2.0 * g_0_xy_xz_xy[i] * b_exp + 4.0 * g_xx_xy_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_xz_xz[i] = -2.0 * g_0_xy_xz_xz[i] * b_exp + 4.0 * g_xx_xy_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_xz_yy[i] = -2.0 * g_0_xy_xz_yy[i] * b_exp + 4.0 * g_xx_xy_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_xz_yz[i] = -2.0 * g_0_xy_xz_yz[i] * b_exp + 4.0 * g_xx_xy_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_xz_zz[i] = -2.0 * g_0_xy_xz_zz[i] * b_exp + 4.0 * g_xx_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_0_xy_yy_xx, g_0_xy_yy_xy, g_0_xy_yy_xz, g_0_xy_yy_yy, g_0_xy_yy_yz, g_0_xy_yy_zz, g_x_y_0_0_x_x_yy_xx, g_x_y_0_0_x_x_yy_xy, g_x_y_0_0_x_x_yy_xz, g_x_y_0_0_x_x_yy_yy, g_x_y_0_0_x_x_yy_yz, g_x_y_0_0_x_x_yy_zz, g_xx_xy_yy_xx, g_xx_xy_yy_xy, g_xx_xy_yy_xz, g_xx_xy_yy_yy, g_xx_xy_yy_yz, g_xx_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_x_yy_xx[i] = -2.0 * g_0_xy_yy_xx[i] * b_exp + 4.0 * g_xx_xy_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_yy_xy[i] = -2.0 * g_0_xy_yy_xy[i] * b_exp + 4.0 * g_xx_xy_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_yy_xz[i] = -2.0 * g_0_xy_yy_xz[i] * b_exp + 4.0 * g_xx_xy_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_yy_yy[i] = -2.0 * g_0_xy_yy_yy[i] * b_exp + 4.0 * g_xx_xy_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_yy_yz[i] = -2.0 * g_0_xy_yy_yz[i] * b_exp + 4.0 * g_xx_xy_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_yy_zz[i] = -2.0 * g_0_xy_yy_zz[i] * b_exp + 4.0 * g_xx_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_0_xy_yz_xx, g_0_xy_yz_xy, g_0_xy_yz_xz, g_0_xy_yz_yy, g_0_xy_yz_yz, g_0_xy_yz_zz, g_x_y_0_0_x_x_yz_xx, g_x_y_0_0_x_x_yz_xy, g_x_y_0_0_x_x_yz_xz, g_x_y_0_0_x_x_yz_yy, g_x_y_0_0_x_x_yz_yz, g_x_y_0_0_x_x_yz_zz, g_xx_xy_yz_xx, g_xx_xy_yz_xy, g_xx_xy_yz_xz, g_xx_xy_yz_yy, g_xx_xy_yz_yz, g_xx_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_x_yz_xx[i] = -2.0 * g_0_xy_yz_xx[i] * b_exp + 4.0 * g_xx_xy_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_yz_xy[i] = -2.0 * g_0_xy_yz_xy[i] * b_exp + 4.0 * g_xx_xy_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_yz_xz[i] = -2.0 * g_0_xy_yz_xz[i] * b_exp + 4.0 * g_xx_xy_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_yz_yy[i] = -2.0 * g_0_xy_yz_yy[i] * b_exp + 4.0 * g_xx_xy_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_yz_yz[i] = -2.0 * g_0_xy_yz_yz[i] * b_exp + 4.0 * g_xx_xy_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_yz_zz[i] = -2.0 * g_0_xy_yz_zz[i] * b_exp + 4.0 * g_xx_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_0_xy_zz_xx, g_0_xy_zz_xy, g_0_xy_zz_xz, g_0_xy_zz_yy, g_0_xy_zz_yz, g_0_xy_zz_zz, g_x_y_0_0_x_x_zz_xx, g_x_y_0_0_x_x_zz_xy, g_x_y_0_0_x_x_zz_xz, g_x_y_0_0_x_x_zz_yy, g_x_y_0_0_x_x_zz_yz, g_x_y_0_0_x_x_zz_zz, g_xx_xy_zz_xx, g_xx_xy_zz_xy, g_xx_xy_zz_xz, g_xx_xy_zz_yy, g_xx_xy_zz_yz, g_xx_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_x_zz_xx[i] = -2.0 * g_0_xy_zz_xx[i] * b_exp + 4.0 * g_xx_xy_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_zz_xy[i] = -2.0 * g_0_xy_zz_xy[i] * b_exp + 4.0 * g_xx_xy_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_zz_xz[i] = -2.0 * g_0_xy_zz_xz[i] * b_exp + 4.0 * g_xx_xy_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_zz_yy[i] = -2.0 * g_0_xy_zz_yy[i] * b_exp + 4.0 * g_xx_xy_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_zz_yz[i] = -2.0 * g_0_xy_zz_yz[i] * b_exp + 4.0 * g_xx_xy_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_zz_zz[i] = -2.0 * g_0_xy_zz_zz[i] * b_exp + 4.0 * g_xx_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_0_0_xx_xx, g_0_0_xx_xy, g_0_0_xx_xz, g_0_0_xx_yy, g_0_0_xx_yz, g_0_0_xx_zz, g_0_yy_xx_xx, g_0_yy_xx_xy, g_0_yy_xx_xz, g_0_yy_xx_yy, g_0_yy_xx_yz, g_0_yy_xx_zz, g_x_y_0_0_x_y_xx_xx, g_x_y_0_0_x_y_xx_xy, g_x_y_0_0_x_y_xx_xz, g_x_y_0_0_x_y_xx_yy, g_x_y_0_0_x_y_xx_yz, g_x_y_0_0_x_y_xx_zz, g_xx_0_xx_xx, g_xx_0_xx_xy, g_xx_0_xx_xz, g_xx_0_xx_yy, g_xx_0_xx_yz, g_xx_0_xx_zz, g_xx_yy_xx_xx, g_xx_yy_xx_xy, g_xx_yy_xx_xz, g_xx_yy_xx_yy, g_xx_yy_xx_yz, g_xx_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_y_xx_xx[i] = g_0_0_xx_xx[i] - 2.0 * g_0_yy_xx_xx[i] * b_exp - 2.0 * g_xx_0_xx_xx[i] * a_exp + 4.0 * g_xx_yy_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_xx_xy[i] = g_0_0_xx_xy[i] - 2.0 * g_0_yy_xx_xy[i] * b_exp - 2.0 * g_xx_0_xx_xy[i] * a_exp + 4.0 * g_xx_yy_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_xx_xz[i] = g_0_0_xx_xz[i] - 2.0 * g_0_yy_xx_xz[i] * b_exp - 2.0 * g_xx_0_xx_xz[i] * a_exp + 4.0 * g_xx_yy_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_xx_yy[i] = g_0_0_xx_yy[i] - 2.0 * g_0_yy_xx_yy[i] * b_exp - 2.0 * g_xx_0_xx_yy[i] * a_exp + 4.0 * g_xx_yy_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_xx_yz[i] = g_0_0_xx_yz[i] - 2.0 * g_0_yy_xx_yz[i] * b_exp - 2.0 * g_xx_0_xx_yz[i] * a_exp + 4.0 * g_xx_yy_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_xx_zz[i] = g_0_0_xx_zz[i] - 2.0 * g_0_yy_xx_zz[i] * b_exp - 2.0 * g_xx_0_xx_zz[i] * a_exp + 4.0 * g_xx_yy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_0_0_xy_xx, g_0_0_xy_xy, g_0_0_xy_xz, g_0_0_xy_yy, g_0_0_xy_yz, g_0_0_xy_zz, g_0_yy_xy_xx, g_0_yy_xy_xy, g_0_yy_xy_xz, g_0_yy_xy_yy, g_0_yy_xy_yz, g_0_yy_xy_zz, g_x_y_0_0_x_y_xy_xx, g_x_y_0_0_x_y_xy_xy, g_x_y_0_0_x_y_xy_xz, g_x_y_0_0_x_y_xy_yy, g_x_y_0_0_x_y_xy_yz, g_x_y_0_0_x_y_xy_zz, g_xx_0_xy_xx, g_xx_0_xy_xy, g_xx_0_xy_xz, g_xx_0_xy_yy, g_xx_0_xy_yz, g_xx_0_xy_zz, g_xx_yy_xy_xx, g_xx_yy_xy_xy, g_xx_yy_xy_xz, g_xx_yy_xy_yy, g_xx_yy_xy_yz, g_xx_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_y_xy_xx[i] = g_0_0_xy_xx[i] - 2.0 * g_0_yy_xy_xx[i] * b_exp - 2.0 * g_xx_0_xy_xx[i] * a_exp + 4.0 * g_xx_yy_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_xy_xy[i] = g_0_0_xy_xy[i] - 2.0 * g_0_yy_xy_xy[i] * b_exp - 2.0 * g_xx_0_xy_xy[i] * a_exp + 4.0 * g_xx_yy_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_xy_xz[i] = g_0_0_xy_xz[i] - 2.0 * g_0_yy_xy_xz[i] * b_exp - 2.0 * g_xx_0_xy_xz[i] * a_exp + 4.0 * g_xx_yy_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_xy_yy[i] = g_0_0_xy_yy[i] - 2.0 * g_0_yy_xy_yy[i] * b_exp - 2.0 * g_xx_0_xy_yy[i] * a_exp + 4.0 * g_xx_yy_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_xy_yz[i] = g_0_0_xy_yz[i] - 2.0 * g_0_yy_xy_yz[i] * b_exp - 2.0 * g_xx_0_xy_yz[i] * a_exp + 4.0 * g_xx_yy_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_xy_zz[i] = g_0_0_xy_zz[i] - 2.0 * g_0_yy_xy_zz[i] * b_exp - 2.0 * g_xx_0_xy_zz[i] * a_exp + 4.0 * g_xx_yy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_0_0_xz_xx, g_0_0_xz_xy, g_0_0_xz_xz, g_0_0_xz_yy, g_0_0_xz_yz, g_0_0_xz_zz, g_0_yy_xz_xx, g_0_yy_xz_xy, g_0_yy_xz_xz, g_0_yy_xz_yy, g_0_yy_xz_yz, g_0_yy_xz_zz, g_x_y_0_0_x_y_xz_xx, g_x_y_0_0_x_y_xz_xy, g_x_y_0_0_x_y_xz_xz, g_x_y_0_0_x_y_xz_yy, g_x_y_0_0_x_y_xz_yz, g_x_y_0_0_x_y_xz_zz, g_xx_0_xz_xx, g_xx_0_xz_xy, g_xx_0_xz_xz, g_xx_0_xz_yy, g_xx_0_xz_yz, g_xx_0_xz_zz, g_xx_yy_xz_xx, g_xx_yy_xz_xy, g_xx_yy_xz_xz, g_xx_yy_xz_yy, g_xx_yy_xz_yz, g_xx_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_y_xz_xx[i] = g_0_0_xz_xx[i] - 2.0 * g_0_yy_xz_xx[i] * b_exp - 2.0 * g_xx_0_xz_xx[i] * a_exp + 4.0 * g_xx_yy_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_xz_xy[i] = g_0_0_xz_xy[i] - 2.0 * g_0_yy_xz_xy[i] * b_exp - 2.0 * g_xx_0_xz_xy[i] * a_exp + 4.0 * g_xx_yy_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_xz_xz[i] = g_0_0_xz_xz[i] - 2.0 * g_0_yy_xz_xz[i] * b_exp - 2.0 * g_xx_0_xz_xz[i] * a_exp + 4.0 * g_xx_yy_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_xz_yy[i] = g_0_0_xz_yy[i] - 2.0 * g_0_yy_xz_yy[i] * b_exp - 2.0 * g_xx_0_xz_yy[i] * a_exp + 4.0 * g_xx_yy_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_xz_yz[i] = g_0_0_xz_yz[i] - 2.0 * g_0_yy_xz_yz[i] * b_exp - 2.0 * g_xx_0_xz_yz[i] * a_exp + 4.0 * g_xx_yy_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_xz_zz[i] = g_0_0_xz_zz[i] - 2.0 * g_0_yy_xz_zz[i] * b_exp - 2.0 * g_xx_0_xz_zz[i] * a_exp + 4.0 * g_xx_yy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_0_0_yy_xx, g_0_0_yy_xy, g_0_0_yy_xz, g_0_0_yy_yy, g_0_0_yy_yz, g_0_0_yy_zz, g_0_yy_yy_xx, g_0_yy_yy_xy, g_0_yy_yy_xz, g_0_yy_yy_yy, g_0_yy_yy_yz, g_0_yy_yy_zz, g_x_y_0_0_x_y_yy_xx, g_x_y_0_0_x_y_yy_xy, g_x_y_0_0_x_y_yy_xz, g_x_y_0_0_x_y_yy_yy, g_x_y_0_0_x_y_yy_yz, g_x_y_0_0_x_y_yy_zz, g_xx_0_yy_xx, g_xx_0_yy_xy, g_xx_0_yy_xz, g_xx_0_yy_yy, g_xx_0_yy_yz, g_xx_0_yy_zz, g_xx_yy_yy_xx, g_xx_yy_yy_xy, g_xx_yy_yy_xz, g_xx_yy_yy_yy, g_xx_yy_yy_yz, g_xx_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_y_yy_xx[i] = g_0_0_yy_xx[i] - 2.0 * g_0_yy_yy_xx[i] * b_exp - 2.0 * g_xx_0_yy_xx[i] * a_exp + 4.0 * g_xx_yy_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_yy_xy[i] = g_0_0_yy_xy[i] - 2.0 * g_0_yy_yy_xy[i] * b_exp - 2.0 * g_xx_0_yy_xy[i] * a_exp + 4.0 * g_xx_yy_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_yy_xz[i] = g_0_0_yy_xz[i] - 2.0 * g_0_yy_yy_xz[i] * b_exp - 2.0 * g_xx_0_yy_xz[i] * a_exp + 4.0 * g_xx_yy_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_yy_yy[i] = g_0_0_yy_yy[i] - 2.0 * g_0_yy_yy_yy[i] * b_exp - 2.0 * g_xx_0_yy_yy[i] * a_exp + 4.0 * g_xx_yy_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_yy_yz[i] = g_0_0_yy_yz[i] - 2.0 * g_0_yy_yy_yz[i] * b_exp - 2.0 * g_xx_0_yy_yz[i] * a_exp + 4.0 * g_xx_yy_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_yy_zz[i] = g_0_0_yy_zz[i] - 2.0 * g_0_yy_yy_zz[i] * b_exp - 2.0 * g_xx_0_yy_zz[i] * a_exp + 4.0 * g_xx_yy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_0_0_yz_xx, g_0_0_yz_xy, g_0_0_yz_xz, g_0_0_yz_yy, g_0_0_yz_yz, g_0_0_yz_zz, g_0_yy_yz_xx, g_0_yy_yz_xy, g_0_yy_yz_xz, g_0_yy_yz_yy, g_0_yy_yz_yz, g_0_yy_yz_zz, g_x_y_0_0_x_y_yz_xx, g_x_y_0_0_x_y_yz_xy, g_x_y_0_0_x_y_yz_xz, g_x_y_0_0_x_y_yz_yy, g_x_y_0_0_x_y_yz_yz, g_x_y_0_0_x_y_yz_zz, g_xx_0_yz_xx, g_xx_0_yz_xy, g_xx_0_yz_xz, g_xx_0_yz_yy, g_xx_0_yz_yz, g_xx_0_yz_zz, g_xx_yy_yz_xx, g_xx_yy_yz_xy, g_xx_yy_yz_xz, g_xx_yy_yz_yy, g_xx_yy_yz_yz, g_xx_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_y_yz_xx[i] = g_0_0_yz_xx[i] - 2.0 * g_0_yy_yz_xx[i] * b_exp - 2.0 * g_xx_0_yz_xx[i] * a_exp + 4.0 * g_xx_yy_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_yz_xy[i] = g_0_0_yz_xy[i] - 2.0 * g_0_yy_yz_xy[i] * b_exp - 2.0 * g_xx_0_yz_xy[i] * a_exp + 4.0 * g_xx_yy_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_yz_xz[i] = g_0_0_yz_xz[i] - 2.0 * g_0_yy_yz_xz[i] * b_exp - 2.0 * g_xx_0_yz_xz[i] * a_exp + 4.0 * g_xx_yy_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_yz_yy[i] = g_0_0_yz_yy[i] - 2.0 * g_0_yy_yz_yy[i] * b_exp - 2.0 * g_xx_0_yz_yy[i] * a_exp + 4.0 * g_xx_yy_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_yz_yz[i] = g_0_0_yz_yz[i] - 2.0 * g_0_yy_yz_yz[i] * b_exp - 2.0 * g_xx_0_yz_yz[i] * a_exp + 4.0 * g_xx_yy_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_yz_zz[i] = g_0_0_yz_zz[i] - 2.0 * g_0_yy_yz_zz[i] * b_exp - 2.0 * g_xx_0_yz_zz[i] * a_exp + 4.0 * g_xx_yy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_0_0_zz_xx, g_0_0_zz_xy, g_0_0_zz_xz, g_0_0_zz_yy, g_0_0_zz_yz, g_0_0_zz_zz, g_0_yy_zz_xx, g_0_yy_zz_xy, g_0_yy_zz_xz, g_0_yy_zz_yy, g_0_yy_zz_yz, g_0_yy_zz_zz, g_x_y_0_0_x_y_zz_xx, g_x_y_0_0_x_y_zz_xy, g_x_y_0_0_x_y_zz_xz, g_x_y_0_0_x_y_zz_yy, g_x_y_0_0_x_y_zz_yz, g_x_y_0_0_x_y_zz_zz, g_xx_0_zz_xx, g_xx_0_zz_xy, g_xx_0_zz_xz, g_xx_0_zz_yy, g_xx_0_zz_yz, g_xx_0_zz_zz, g_xx_yy_zz_xx, g_xx_yy_zz_xy, g_xx_yy_zz_xz, g_xx_yy_zz_yy, g_xx_yy_zz_yz, g_xx_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_y_zz_xx[i] = g_0_0_zz_xx[i] - 2.0 * g_0_yy_zz_xx[i] * b_exp - 2.0 * g_xx_0_zz_xx[i] * a_exp + 4.0 * g_xx_yy_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_zz_xy[i] = g_0_0_zz_xy[i] - 2.0 * g_0_yy_zz_xy[i] * b_exp - 2.0 * g_xx_0_zz_xy[i] * a_exp + 4.0 * g_xx_yy_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_zz_xz[i] = g_0_0_zz_xz[i] - 2.0 * g_0_yy_zz_xz[i] * b_exp - 2.0 * g_xx_0_zz_xz[i] * a_exp + 4.0 * g_xx_yy_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_zz_yy[i] = g_0_0_zz_yy[i] - 2.0 * g_0_yy_zz_yy[i] * b_exp - 2.0 * g_xx_0_zz_yy[i] * a_exp + 4.0 * g_xx_yy_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_zz_yz[i] = g_0_0_zz_yz[i] - 2.0 * g_0_yy_zz_yz[i] * b_exp - 2.0 * g_xx_0_zz_yz[i] * a_exp + 4.0 * g_xx_yy_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_zz_zz[i] = g_0_0_zz_zz[i] - 2.0 * g_0_yy_zz_zz[i] * b_exp - 2.0 * g_xx_0_zz_zz[i] * a_exp + 4.0 * g_xx_yy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_0_yz_xx_xx, g_0_yz_xx_xy, g_0_yz_xx_xz, g_0_yz_xx_yy, g_0_yz_xx_yz, g_0_yz_xx_zz, g_x_y_0_0_x_z_xx_xx, g_x_y_0_0_x_z_xx_xy, g_x_y_0_0_x_z_xx_xz, g_x_y_0_0_x_z_xx_yy, g_x_y_0_0_x_z_xx_yz, g_x_y_0_0_x_z_xx_zz, g_xx_yz_xx_xx, g_xx_yz_xx_xy, g_xx_yz_xx_xz, g_xx_yz_xx_yy, g_xx_yz_xx_yz, g_xx_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_z_xx_xx[i] = -2.0 * g_0_yz_xx_xx[i] * b_exp + 4.0 * g_xx_yz_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_xx_xy[i] = -2.0 * g_0_yz_xx_xy[i] * b_exp + 4.0 * g_xx_yz_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_xx_xz[i] = -2.0 * g_0_yz_xx_xz[i] * b_exp + 4.0 * g_xx_yz_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_xx_yy[i] = -2.0 * g_0_yz_xx_yy[i] * b_exp + 4.0 * g_xx_yz_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_xx_yz[i] = -2.0 * g_0_yz_xx_yz[i] * b_exp + 4.0 * g_xx_yz_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_xx_zz[i] = -2.0 * g_0_yz_xx_zz[i] * b_exp + 4.0 * g_xx_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_0_yz_xy_xx, g_0_yz_xy_xy, g_0_yz_xy_xz, g_0_yz_xy_yy, g_0_yz_xy_yz, g_0_yz_xy_zz, g_x_y_0_0_x_z_xy_xx, g_x_y_0_0_x_z_xy_xy, g_x_y_0_0_x_z_xy_xz, g_x_y_0_0_x_z_xy_yy, g_x_y_0_0_x_z_xy_yz, g_x_y_0_0_x_z_xy_zz, g_xx_yz_xy_xx, g_xx_yz_xy_xy, g_xx_yz_xy_xz, g_xx_yz_xy_yy, g_xx_yz_xy_yz, g_xx_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_z_xy_xx[i] = -2.0 * g_0_yz_xy_xx[i] * b_exp + 4.0 * g_xx_yz_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_xy_xy[i] = -2.0 * g_0_yz_xy_xy[i] * b_exp + 4.0 * g_xx_yz_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_xy_xz[i] = -2.0 * g_0_yz_xy_xz[i] * b_exp + 4.0 * g_xx_yz_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_xy_yy[i] = -2.0 * g_0_yz_xy_yy[i] * b_exp + 4.0 * g_xx_yz_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_xy_yz[i] = -2.0 * g_0_yz_xy_yz[i] * b_exp + 4.0 * g_xx_yz_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_xy_zz[i] = -2.0 * g_0_yz_xy_zz[i] * b_exp + 4.0 * g_xx_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_0_yz_xz_xx, g_0_yz_xz_xy, g_0_yz_xz_xz, g_0_yz_xz_yy, g_0_yz_xz_yz, g_0_yz_xz_zz, g_x_y_0_0_x_z_xz_xx, g_x_y_0_0_x_z_xz_xy, g_x_y_0_0_x_z_xz_xz, g_x_y_0_0_x_z_xz_yy, g_x_y_0_0_x_z_xz_yz, g_x_y_0_0_x_z_xz_zz, g_xx_yz_xz_xx, g_xx_yz_xz_xy, g_xx_yz_xz_xz, g_xx_yz_xz_yy, g_xx_yz_xz_yz, g_xx_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_z_xz_xx[i] = -2.0 * g_0_yz_xz_xx[i] * b_exp + 4.0 * g_xx_yz_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_xz_xy[i] = -2.0 * g_0_yz_xz_xy[i] * b_exp + 4.0 * g_xx_yz_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_xz_xz[i] = -2.0 * g_0_yz_xz_xz[i] * b_exp + 4.0 * g_xx_yz_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_xz_yy[i] = -2.0 * g_0_yz_xz_yy[i] * b_exp + 4.0 * g_xx_yz_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_xz_yz[i] = -2.0 * g_0_yz_xz_yz[i] * b_exp + 4.0 * g_xx_yz_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_xz_zz[i] = -2.0 * g_0_yz_xz_zz[i] * b_exp + 4.0 * g_xx_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_0_yz_yy_xx, g_0_yz_yy_xy, g_0_yz_yy_xz, g_0_yz_yy_yy, g_0_yz_yy_yz, g_0_yz_yy_zz, g_x_y_0_0_x_z_yy_xx, g_x_y_0_0_x_z_yy_xy, g_x_y_0_0_x_z_yy_xz, g_x_y_0_0_x_z_yy_yy, g_x_y_0_0_x_z_yy_yz, g_x_y_0_0_x_z_yy_zz, g_xx_yz_yy_xx, g_xx_yz_yy_xy, g_xx_yz_yy_xz, g_xx_yz_yy_yy, g_xx_yz_yy_yz, g_xx_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_z_yy_xx[i] = -2.0 * g_0_yz_yy_xx[i] * b_exp + 4.0 * g_xx_yz_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_yy_xy[i] = -2.0 * g_0_yz_yy_xy[i] * b_exp + 4.0 * g_xx_yz_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_yy_xz[i] = -2.0 * g_0_yz_yy_xz[i] * b_exp + 4.0 * g_xx_yz_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_yy_yy[i] = -2.0 * g_0_yz_yy_yy[i] * b_exp + 4.0 * g_xx_yz_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_yy_yz[i] = -2.0 * g_0_yz_yy_yz[i] * b_exp + 4.0 * g_xx_yz_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_yy_zz[i] = -2.0 * g_0_yz_yy_zz[i] * b_exp + 4.0 * g_xx_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_0_yz_yz_xx, g_0_yz_yz_xy, g_0_yz_yz_xz, g_0_yz_yz_yy, g_0_yz_yz_yz, g_0_yz_yz_zz, g_x_y_0_0_x_z_yz_xx, g_x_y_0_0_x_z_yz_xy, g_x_y_0_0_x_z_yz_xz, g_x_y_0_0_x_z_yz_yy, g_x_y_0_0_x_z_yz_yz, g_x_y_0_0_x_z_yz_zz, g_xx_yz_yz_xx, g_xx_yz_yz_xy, g_xx_yz_yz_xz, g_xx_yz_yz_yy, g_xx_yz_yz_yz, g_xx_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_z_yz_xx[i] = -2.0 * g_0_yz_yz_xx[i] * b_exp + 4.0 * g_xx_yz_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_yz_xy[i] = -2.0 * g_0_yz_yz_xy[i] * b_exp + 4.0 * g_xx_yz_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_yz_xz[i] = -2.0 * g_0_yz_yz_xz[i] * b_exp + 4.0 * g_xx_yz_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_yz_yy[i] = -2.0 * g_0_yz_yz_yy[i] * b_exp + 4.0 * g_xx_yz_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_yz_yz[i] = -2.0 * g_0_yz_yz_yz[i] * b_exp + 4.0 * g_xx_yz_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_yz_zz[i] = -2.0 * g_0_yz_yz_zz[i] * b_exp + 4.0 * g_xx_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_0_yz_zz_xx, g_0_yz_zz_xy, g_0_yz_zz_xz, g_0_yz_zz_yy, g_0_yz_zz_yz, g_0_yz_zz_zz, g_x_y_0_0_x_z_zz_xx, g_x_y_0_0_x_z_zz_xy, g_x_y_0_0_x_z_zz_xz, g_x_y_0_0_x_z_zz_yy, g_x_y_0_0_x_z_zz_yz, g_x_y_0_0_x_z_zz_zz, g_xx_yz_zz_xx, g_xx_yz_zz_xy, g_xx_yz_zz_xz, g_xx_yz_zz_yy, g_xx_yz_zz_yz, g_xx_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_z_zz_xx[i] = -2.0 * g_0_yz_zz_xx[i] * b_exp + 4.0 * g_xx_yz_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_zz_xy[i] = -2.0 * g_0_yz_zz_xy[i] * b_exp + 4.0 * g_xx_yz_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_zz_xz[i] = -2.0 * g_0_yz_zz_xz[i] * b_exp + 4.0 * g_xx_yz_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_zz_yy[i] = -2.0 * g_0_yz_zz_yy[i] * b_exp + 4.0 * g_xx_yz_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_zz_yz[i] = -2.0 * g_0_yz_zz_yz[i] * b_exp + 4.0 * g_xx_yz_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_zz_zz[i] = -2.0 * g_0_yz_zz_zz[i] * b_exp + 4.0 * g_xx_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_x_y_0_0_y_x_xx_xx, g_x_y_0_0_y_x_xx_xy, g_x_y_0_0_y_x_xx_xz, g_x_y_0_0_y_x_xx_yy, g_x_y_0_0_y_x_xx_yz, g_x_y_0_0_y_x_xx_zz, g_xy_xy_xx_xx, g_xy_xy_xx_xy, g_xy_xy_xx_xz, g_xy_xy_xx_yy, g_xy_xy_xx_yz, g_xy_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_x_xx_xx[i] = 4.0 * g_xy_xy_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_xx_xy[i] = 4.0 * g_xy_xy_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_xx_xz[i] = 4.0 * g_xy_xy_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_xx_yy[i] = 4.0 * g_xy_xy_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_xx_yz[i] = 4.0 * g_xy_xy_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_xx_zz[i] = 4.0 * g_xy_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_x_y_0_0_y_x_xy_xx, g_x_y_0_0_y_x_xy_xy, g_x_y_0_0_y_x_xy_xz, g_x_y_0_0_y_x_xy_yy, g_x_y_0_0_y_x_xy_yz, g_x_y_0_0_y_x_xy_zz, g_xy_xy_xy_xx, g_xy_xy_xy_xy, g_xy_xy_xy_xz, g_xy_xy_xy_yy, g_xy_xy_xy_yz, g_xy_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_x_xy_xx[i] = 4.0 * g_xy_xy_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_xy_xy[i] = 4.0 * g_xy_xy_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_xy_xz[i] = 4.0 * g_xy_xy_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_xy_yy[i] = 4.0 * g_xy_xy_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_xy_yz[i] = 4.0 * g_xy_xy_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_xy_zz[i] = 4.0 * g_xy_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_x_y_0_0_y_x_xz_xx, g_x_y_0_0_y_x_xz_xy, g_x_y_0_0_y_x_xz_xz, g_x_y_0_0_y_x_xz_yy, g_x_y_0_0_y_x_xz_yz, g_x_y_0_0_y_x_xz_zz, g_xy_xy_xz_xx, g_xy_xy_xz_xy, g_xy_xy_xz_xz, g_xy_xy_xz_yy, g_xy_xy_xz_yz, g_xy_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_x_xz_xx[i] = 4.0 * g_xy_xy_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_xz_xy[i] = 4.0 * g_xy_xy_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_xz_xz[i] = 4.0 * g_xy_xy_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_xz_yy[i] = 4.0 * g_xy_xy_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_xz_yz[i] = 4.0 * g_xy_xy_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_xz_zz[i] = 4.0 * g_xy_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_x_y_0_0_y_x_yy_xx, g_x_y_0_0_y_x_yy_xy, g_x_y_0_0_y_x_yy_xz, g_x_y_0_0_y_x_yy_yy, g_x_y_0_0_y_x_yy_yz, g_x_y_0_0_y_x_yy_zz, g_xy_xy_yy_xx, g_xy_xy_yy_xy, g_xy_xy_yy_xz, g_xy_xy_yy_yy, g_xy_xy_yy_yz, g_xy_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_x_yy_xx[i] = 4.0 * g_xy_xy_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_yy_xy[i] = 4.0 * g_xy_xy_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_yy_xz[i] = 4.0 * g_xy_xy_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_yy_yy[i] = 4.0 * g_xy_xy_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_yy_yz[i] = 4.0 * g_xy_xy_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_yy_zz[i] = 4.0 * g_xy_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_x_y_0_0_y_x_yz_xx, g_x_y_0_0_y_x_yz_xy, g_x_y_0_0_y_x_yz_xz, g_x_y_0_0_y_x_yz_yy, g_x_y_0_0_y_x_yz_yz, g_x_y_0_0_y_x_yz_zz, g_xy_xy_yz_xx, g_xy_xy_yz_xy, g_xy_xy_yz_xz, g_xy_xy_yz_yy, g_xy_xy_yz_yz, g_xy_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_x_yz_xx[i] = 4.0 * g_xy_xy_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_yz_xy[i] = 4.0 * g_xy_xy_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_yz_xz[i] = 4.0 * g_xy_xy_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_yz_yy[i] = 4.0 * g_xy_xy_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_yz_yz[i] = 4.0 * g_xy_xy_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_yz_zz[i] = 4.0 * g_xy_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_x_y_0_0_y_x_zz_xx, g_x_y_0_0_y_x_zz_xy, g_x_y_0_0_y_x_zz_xz, g_x_y_0_0_y_x_zz_yy, g_x_y_0_0_y_x_zz_yz, g_x_y_0_0_y_x_zz_zz, g_xy_xy_zz_xx, g_xy_xy_zz_xy, g_xy_xy_zz_xz, g_xy_xy_zz_yy, g_xy_xy_zz_yz, g_xy_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_x_zz_xx[i] = 4.0 * g_xy_xy_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_zz_xy[i] = 4.0 * g_xy_xy_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_zz_xz[i] = 4.0 * g_xy_xy_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_zz_yy[i] = 4.0 * g_xy_xy_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_zz_yz[i] = 4.0 * g_xy_xy_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_zz_zz[i] = 4.0 * g_xy_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_x_y_0_0_y_y_xx_xx, g_x_y_0_0_y_y_xx_xy, g_x_y_0_0_y_y_xx_xz, g_x_y_0_0_y_y_xx_yy, g_x_y_0_0_y_y_xx_yz, g_x_y_0_0_y_y_xx_zz, g_xy_0_xx_xx, g_xy_0_xx_xy, g_xy_0_xx_xz, g_xy_0_xx_yy, g_xy_0_xx_yz, g_xy_0_xx_zz, g_xy_yy_xx_xx, g_xy_yy_xx_xy, g_xy_yy_xx_xz, g_xy_yy_xx_yy, g_xy_yy_xx_yz, g_xy_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_y_xx_xx[i] = -2.0 * g_xy_0_xx_xx[i] * a_exp + 4.0 * g_xy_yy_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_xx_xy[i] = -2.0 * g_xy_0_xx_xy[i] * a_exp + 4.0 * g_xy_yy_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_xx_xz[i] = -2.0 * g_xy_0_xx_xz[i] * a_exp + 4.0 * g_xy_yy_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_xx_yy[i] = -2.0 * g_xy_0_xx_yy[i] * a_exp + 4.0 * g_xy_yy_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_xx_yz[i] = -2.0 * g_xy_0_xx_yz[i] * a_exp + 4.0 * g_xy_yy_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_xx_zz[i] = -2.0 * g_xy_0_xx_zz[i] * a_exp + 4.0 * g_xy_yy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_x_y_0_0_y_y_xy_xx, g_x_y_0_0_y_y_xy_xy, g_x_y_0_0_y_y_xy_xz, g_x_y_0_0_y_y_xy_yy, g_x_y_0_0_y_y_xy_yz, g_x_y_0_0_y_y_xy_zz, g_xy_0_xy_xx, g_xy_0_xy_xy, g_xy_0_xy_xz, g_xy_0_xy_yy, g_xy_0_xy_yz, g_xy_0_xy_zz, g_xy_yy_xy_xx, g_xy_yy_xy_xy, g_xy_yy_xy_xz, g_xy_yy_xy_yy, g_xy_yy_xy_yz, g_xy_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_y_xy_xx[i] = -2.0 * g_xy_0_xy_xx[i] * a_exp + 4.0 * g_xy_yy_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_xy_xy[i] = -2.0 * g_xy_0_xy_xy[i] * a_exp + 4.0 * g_xy_yy_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_xy_xz[i] = -2.0 * g_xy_0_xy_xz[i] * a_exp + 4.0 * g_xy_yy_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_xy_yy[i] = -2.0 * g_xy_0_xy_yy[i] * a_exp + 4.0 * g_xy_yy_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_xy_yz[i] = -2.0 * g_xy_0_xy_yz[i] * a_exp + 4.0 * g_xy_yy_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_xy_zz[i] = -2.0 * g_xy_0_xy_zz[i] * a_exp + 4.0 * g_xy_yy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_x_y_0_0_y_y_xz_xx, g_x_y_0_0_y_y_xz_xy, g_x_y_0_0_y_y_xz_xz, g_x_y_0_0_y_y_xz_yy, g_x_y_0_0_y_y_xz_yz, g_x_y_0_0_y_y_xz_zz, g_xy_0_xz_xx, g_xy_0_xz_xy, g_xy_0_xz_xz, g_xy_0_xz_yy, g_xy_0_xz_yz, g_xy_0_xz_zz, g_xy_yy_xz_xx, g_xy_yy_xz_xy, g_xy_yy_xz_xz, g_xy_yy_xz_yy, g_xy_yy_xz_yz, g_xy_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_y_xz_xx[i] = -2.0 * g_xy_0_xz_xx[i] * a_exp + 4.0 * g_xy_yy_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_xz_xy[i] = -2.0 * g_xy_0_xz_xy[i] * a_exp + 4.0 * g_xy_yy_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_xz_xz[i] = -2.0 * g_xy_0_xz_xz[i] * a_exp + 4.0 * g_xy_yy_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_xz_yy[i] = -2.0 * g_xy_0_xz_yy[i] * a_exp + 4.0 * g_xy_yy_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_xz_yz[i] = -2.0 * g_xy_0_xz_yz[i] * a_exp + 4.0 * g_xy_yy_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_xz_zz[i] = -2.0 * g_xy_0_xz_zz[i] * a_exp + 4.0 * g_xy_yy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_x_y_0_0_y_y_yy_xx, g_x_y_0_0_y_y_yy_xy, g_x_y_0_0_y_y_yy_xz, g_x_y_0_0_y_y_yy_yy, g_x_y_0_0_y_y_yy_yz, g_x_y_0_0_y_y_yy_zz, g_xy_0_yy_xx, g_xy_0_yy_xy, g_xy_0_yy_xz, g_xy_0_yy_yy, g_xy_0_yy_yz, g_xy_0_yy_zz, g_xy_yy_yy_xx, g_xy_yy_yy_xy, g_xy_yy_yy_xz, g_xy_yy_yy_yy, g_xy_yy_yy_yz, g_xy_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_y_yy_xx[i] = -2.0 * g_xy_0_yy_xx[i] * a_exp + 4.0 * g_xy_yy_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_yy_xy[i] = -2.0 * g_xy_0_yy_xy[i] * a_exp + 4.0 * g_xy_yy_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_yy_xz[i] = -2.0 * g_xy_0_yy_xz[i] * a_exp + 4.0 * g_xy_yy_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_yy_yy[i] = -2.0 * g_xy_0_yy_yy[i] * a_exp + 4.0 * g_xy_yy_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_yy_yz[i] = -2.0 * g_xy_0_yy_yz[i] * a_exp + 4.0 * g_xy_yy_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_yy_zz[i] = -2.0 * g_xy_0_yy_zz[i] * a_exp + 4.0 * g_xy_yy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_x_y_0_0_y_y_yz_xx, g_x_y_0_0_y_y_yz_xy, g_x_y_0_0_y_y_yz_xz, g_x_y_0_0_y_y_yz_yy, g_x_y_0_0_y_y_yz_yz, g_x_y_0_0_y_y_yz_zz, g_xy_0_yz_xx, g_xy_0_yz_xy, g_xy_0_yz_xz, g_xy_0_yz_yy, g_xy_0_yz_yz, g_xy_0_yz_zz, g_xy_yy_yz_xx, g_xy_yy_yz_xy, g_xy_yy_yz_xz, g_xy_yy_yz_yy, g_xy_yy_yz_yz, g_xy_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_y_yz_xx[i] = -2.0 * g_xy_0_yz_xx[i] * a_exp + 4.0 * g_xy_yy_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_yz_xy[i] = -2.0 * g_xy_0_yz_xy[i] * a_exp + 4.0 * g_xy_yy_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_yz_xz[i] = -2.0 * g_xy_0_yz_xz[i] * a_exp + 4.0 * g_xy_yy_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_yz_yy[i] = -2.0 * g_xy_0_yz_yy[i] * a_exp + 4.0 * g_xy_yy_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_yz_yz[i] = -2.0 * g_xy_0_yz_yz[i] * a_exp + 4.0 * g_xy_yy_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_yz_zz[i] = -2.0 * g_xy_0_yz_zz[i] * a_exp + 4.0 * g_xy_yy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_x_y_0_0_y_y_zz_xx, g_x_y_0_0_y_y_zz_xy, g_x_y_0_0_y_y_zz_xz, g_x_y_0_0_y_y_zz_yy, g_x_y_0_0_y_y_zz_yz, g_x_y_0_0_y_y_zz_zz, g_xy_0_zz_xx, g_xy_0_zz_xy, g_xy_0_zz_xz, g_xy_0_zz_yy, g_xy_0_zz_yz, g_xy_0_zz_zz, g_xy_yy_zz_xx, g_xy_yy_zz_xy, g_xy_yy_zz_xz, g_xy_yy_zz_yy, g_xy_yy_zz_yz, g_xy_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_y_zz_xx[i] = -2.0 * g_xy_0_zz_xx[i] * a_exp + 4.0 * g_xy_yy_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_zz_xy[i] = -2.0 * g_xy_0_zz_xy[i] * a_exp + 4.0 * g_xy_yy_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_zz_xz[i] = -2.0 * g_xy_0_zz_xz[i] * a_exp + 4.0 * g_xy_yy_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_zz_yy[i] = -2.0 * g_xy_0_zz_yy[i] * a_exp + 4.0 * g_xy_yy_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_zz_yz[i] = -2.0 * g_xy_0_zz_yz[i] * a_exp + 4.0 * g_xy_yy_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_zz_zz[i] = -2.0 * g_xy_0_zz_zz[i] * a_exp + 4.0 * g_xy_yy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_x_y_0_0_y_z_xx_xx, g_x_y_0_0_y_z_xx_xy, g_x_y_0_0_y_z_xx_xz, g_x_y_0_0_y_z_xx_yy, g_x_y_0_0_y_z_xx_yz, g_x_y_0_0_y_z_xx_zz, g_xy_yz_xx_xx, g_xy_yz_xx_xy, g_xy_yz_xx_xz, g_xy_yz_xx_yy, g_xy_yz_xx_yz, g_xy_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_z_xx_xx[i] = 4.0 * g_xy_yz_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_xx_xy[i] = 4.0 * g_xy_yz_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_xx_xz[i] = 4.0 * g_xy_yz_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_xx_yy[i] = 4.0 * g_xy_yz_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_xx_yz[i] = 4.0 * g_xy_yz_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_xx_zz[i] = 4.0 * g_xy_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_x_y_0_0_y_z_xy_xx, g_x_y_0_0_y_z_xy_xy, g_x_y_0_0_y_z_xy_xz, g_x_y_0_0_y_z_xy_yy, g_x_y_0_0_y_z_xy_yz, g_x_y_0_0_y_z_xy_zz, g_xy_yz_xy_xx, g_xy_yz_xy_xy, g_xy_yz_xy_xz, g_xy_yz_xy_yy, g_xy_yz_xy_yz, g_xy_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_z_xy_xx[i] = 4.0 * g_xy_yz_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_xy_xy[i] = 4.0 * g_xy_yz_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_xy_xz[i] = 4.0 * g_xy_yz_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_xy_yy[i] = 4.0 * g_xy_yz_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_xy_yz[i] = 4.0 * g_xy_yz_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_xy_zz[i] = 4.0 * g_xy_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_x_y_0_0_y_z_xz_xx, g_x_y_0_0_y_z_xz_xy, g_x_y_0_0_y_z_xz_xz, g_x_y_0_0_y_z_xz_yy, g_x_y_0_0_y_z_xz_yz, g_x_y_0_0_y_z_xz_zz, g_xy_yz_xz_xx, g_xy_yz_xz_xy, g_xy_yz_xz_xz, g_xy_yz_xz_yy, g_xy_yz_xz_yz, g_xy_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_z_xz_xx[i] = 4.0 * g_xy_yz_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_xz_xy[i] = 4.0 * g_xy_yz_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_xz_xz[i] = 4.0 * g_xy_yz_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_xz_yy[i] = 4.0 * g_xy_yz_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_xz_yz[i] = 4.0 * g_xy_yz_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_xz_zz[i] = 4.0 * g_xy_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_x_y_0_0_y_z_yy_xx, g_x_y_0_0_y_z_yy_xy, g_x_y_0_0_y_z_yy_xz, g_x_y_0_0_y_z_yy_yy, g_x_y_0_0_y_z_yy_yz, g_x_y_0_0_y_z_yy_zz, g_xy_yz_yy_xx, g_xy_yz_yy_xy, g_xy_yz_yy_xz, g_xy_yz_yy_yy, g_xy_yz_yy_yz, g_xy_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_z_yy_xx[i] = 4.0 * g_xy_yz_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_yy_xy[i] = 4.0 * g_xy_yz_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_yy_xz[i] = 4.0 * g_xy_yz_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_yy_yy[i] = 4.0 * g_xy_yz_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_yy_yz[i] = 4.0 * g_xy_yz_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_yy_zz[i] = 4.0 * g_xy_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_x_y_0_0_y_z_yz_xx, g_x_y_0_0_y_z_yz_xy, g_x_y_0_0_y_z_yz_xz, g_x_y_0_0_y_z_yz_yy, g_x_y_0_0_y_z_yz_yz, g_x_y_0_0_y_z_yz_zz, g_xy_yz_yz_xx, g_xy_yz_yz_xy, g_xy_yz_yz_xz, g_xy_yz_yz_yy, g_xy_yz_yz_yz, g_xy_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_z_yz_xx[i] = 4.0 * g_xy_yz_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_yz_xy[i] = 4.0 * g_xy_yz_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_yz_xz[i] = 4.0 * g_xy_yz_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_yz_yy[i] = 4.0 * g_xy_yz_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_yz_yz[i] = 4.0 * g_xy_yz_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_yz_zz[i] = 4.0 * g_xy_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_x_y_0_0_y_z_zz_xx, g_x_y_0_0_y_z_zz_xy, g_x_y_0_0_y_z_zz_xz, g_x_y_0_0_y_z_zz_yy, g_x_y_0_0_y_z_zz_yz, g_x_y_0_0_y_z_zz_zz, g_xy_yz_zz_xx, g_xy_yz_zz_xy, g_xy_yz_zz_xz, g_xy_yz_zz_yy, g_xy_yz_zz_yz, g_xy_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_z_zz_xx[i] = 4.0 * g_xy_yz_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_zz_xy[i] = 4.0 * g_xy_yz_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_zz_xz[i] = 4.0 * g_xy_yz_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_zz_yy[i] = 4.0 * g_xy_yz_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_zz_yz[i] = 4.0 * g_xy_yz_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_zz_zz[i] = 4.0 * g_xy_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_x_y_0_0_z_x_xx_xx, g_x_y_0_0_z_x_xx_xy, g_x_y_0_0_z_x_xx_xz, g_x_y_0_0_z_x_xx_yy, g_x_y_0_0_z_x_xx_yz, g_x_y_0_0_z_x_xx_zz, g_xz_xy_xx_xx, g_xz_xy_xx_xy, g_xz_xy_xx_xz, g_xz_xy_xx_yy, g_xz_xy_xx_yz, g_xz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_x_xx_xx[i] = 4.0 * g_xz_xy_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_xx_xy[i] = 4.0 * g_xz_xy_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_xx_xz[i] = 4.0 * g_xz_xy_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_xx_yy[i] = 4.0 * g_xz_xy_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_xx_yz[i] = 4.0 * g_xz_xy_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_xx_zz[i] = 4.0 * g_xz_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_x_y_0_0_z_x_xy_xx, g_x_y_0_0_z_x_xy_xy, g_x_y_0_0_z_x_xy_xz, g_x_y_0_0_z_x_xy_yy, g_x_y_0_0_z_x_xy_yz, g_x_y_0_0_z_x_xy_zz, g_xz_xy_xy_xx, g_xz_xy_xy_xy, g_xz_xy_xy_xz, g_xz_xy_xy_yy, g_xz_xy_xy_yz, g_xz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_x_xy_xx[i] = 4.0 * g_xz_xy_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_xy_xy[i] = 4.0 * g_xz_xy_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_xy_xz[i] = 4.0 * g_xz_xy_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_xy_yy[i] = 4.0 * g_xz_xy_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_xy_yz[i] = 4.0 * g_xz_xy_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_xy_zz[i] = 4.0 * g_xz_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_x_y_0_0_z_x_xz_xx, g_x_y_0_0_z_x_xz_xy, g_x_y_0_0_z_x_xz_xz, g_x_y_0_0_z_x_xz_yy, g_x_y_0_0_z_x_xz_yz, g_x_y_0_0_z_x_xz_zz, g_xz_xy_xz_xx, g_xz_xy_xz_xy, g_xz_xy_xz_xz, g_xz_xy_xz_yy, g_xz_xy_xz_yz, g_xz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_x_xz_xx[i] = 4.0 * g_xz_xy_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_xz_xy[i] = 4.0 * g_xz_xy_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_xz_xz[i] = 4.0 * g_xz_xy_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_xz_yy[i] = 4.0 * g_xz_xy_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_xz_yz[i] = 4.0 * g_xz_xy_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_xz_zz[i] = 4.0 * g_xz_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_x_y_0_0_z_x_yy_xx, g_x_y_0_0_z_x_yy_xy, g_x_y_0_0_z_x_yy_xz, g_x_y_0_0_z_x_yy_yy, g_x_y_0_0_z_x_yy_yz, g_x_y_0_0_z_x_yy_zz, g_xz_xy_yy_xx, g_xz_xy_yy_xy, g_xz_xy_yy_xz, g_xz_xy_yy_yy, g_xz_xy_yy_yz, g_xz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_x_yy_xx[i] = 4.0 * g_xz_xy_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_yy_xy[i] = 4.0 * g_xz_xy_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_yy_xz[i] = 4.0 * g_xz_xy_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_yy_yy[i] = 4.0 * g_xz_xy_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_yy_yz[i] = 4.0 * g_xz_xy_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_yy_zz[i] = 4.0 * g_xz_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_x_y_0_0_z_x_yz_xx, g_x_y_0_0_z_x_yz_xy, g_x_y_0_0_z_x_yz_xz, g_x_y_0_0_z_x_yz_yy, g_x_y_0_0_z_x_yz_yz, g_x_y_0_0_z_x_yz_zz, g_xz_xy_yz_xx, g_xz_xy_yz_xy, g_xz_xy_yz_xz, g_xz_xy_yz_yy, g_xz_xy_yz_yz, g_xz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_x_yz_xx[i] = 4.0 * g_xz_xy_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_yz_xy[i] = 4.0 * g_xz_xy_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_yz_xz[i] = 4.0 * g_xz_xy_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_yz_yy[i] = 4.0 * g_xz_xy_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_yz_yz[i] = 4.0 * g_xz_xy_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_yz_zz[i] = 4.0 * g_xz_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_x_y_0_0_z_x_zz_xx, g_x_y_0_0_z_x_zz_xy, g_x_y_0_0_z_x_zz_xz, g_x_y_0_0_z_x_zz_yy, g_x_y_0_0_z_x_zz_yz, g_x_y_0_0_z_x_zz_zz, g_xz_xy_zz_xx, g_xz_xy_zz_xy, g_xz_xy_zz_xz, g_xz_xy_zz_yy, g_xz_xy_zz_yz, g_xz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_x_zz_xx[i] = 4.0 * g_xz_xy_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_zz_xy[i] = 4.0 * g_xz_xy_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_zz_xz[i] = 4.0 * g_xz_xy_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_zz_yy[i] = 4.0 * g_xz_xy_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_zz_yz[i] = 4.0 * g_xz_xy_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_zz_zz[i] = 4.0 * g_xz_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_x_y_0_0_z_y_xx_xx, g_x_y_0_0_z_y_xx_xy, g_x_y_0_0_z_y_xx_xz, g_x_y_0_0_z_y_xx_yy, g_x_y_0_0_z_y_xx_yz, g_x_y_0_0_z_y_xx_zz, g_xz_0_xx_xx, g_xz_0_xx_xy, g_xz_0_xx_xz, g_xz_0_xx_yy, g_xz_0_xx_yz, g_xz_0_xx_zz, g_xz_yy_xx_xx, g_xz_yy_xx_xy, g_xz_yy_xx_xz, g_xz_yy_xx_yy, g_xz_yy_xx_yz, g_xz_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_y_xx_xx[i] = -2.0 * g_xz_0_xx_xx[i] * a_exp + 4.0 * g_xz_yy_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_xx_xy[i] = -2.0 * g_xz_0_xx_xy[i] * a_exp + 4.0 * g_xz_yy_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_xx_xz[i] = -2.0 * g_xz_0_xx_xz[i] * a_exp + 4.0 * g_xz_yy_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_xx_yy[i] = -2.0 * g_xz_0_xx_yy[i] * a_exp + 4.0 * g_xz_yy_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_xx_yz[i] = -2.0 * g_xz_0_xx_yz[i] * a_exp + 4.0 * g_xz_yy_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_xx_zz[i] = -2.0 * g_xz_0_xx_zz[i] * a_exp + 4.0 * g_xz_yy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_x_y_0_0_z_y_xy_xx, g_x_y_0_0_z_y_xy_xy, g_x_y_0_0_z_y_xy_xz, g_x_y_0_0_z_y_xy_yy, g_x_y_0_0_z_y_xy_yz, g_x_y_0_0_z_y_xy_zz, g_xz_0_xy_xx, g_xz_0_xy_xy, g_xz_0_xy_xz, g_xz_0_xy_yy, g_xz_0_xy_yz, g_xz_0_xy_zz, g_xz_yy_xy_xx, g_xz_yy_xy_xy, g_xz_yy_xy_xz, g_xz_yy_xy_yy, g_xz_yy_xy_yz, g_xz_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_y_xy_xx[i] = -2.0 * g_xz_0_xy_xx[i] * a_exp + 4.0 * g_xz_yy_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_xy_xy[i] = -2.0 * g_xz_0_xy_xy[i] * a_exp + 4.0 * g_xz_yy_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_xy_xz[i] = -2.0 * g_xz_0_xy_xz[i] * a_exp + 4.0 * g_xz_yy_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_xy_yy[i] = -2.0 * g_xz_0_xy_yy[i] * a_exp + 4.0 * g_xz_yy_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_xy_yz[i] = -2.0 * g_xz_0_xy_yz[i] * a_exp + 4.0 * g_xz_yy_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_xy_zz[i] = -2.0 * g_xz_0_xy_zz[i] * a_exp + 4.0 * g_xz_yy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_x_y_0_0_z_y_xz_xx, g_x_y_0_0_z_y_xz_xy, g_x_y_0_0_z_y_xz_xz, g_x_y_0_0_z_y_xz_yy, g_x_y_0_0_z_y_xz_yz, g_x_y_0_0_z_y_xz_zz, g_xz_0_xz_xx, g_xz_0_xz_xy, g_xz_0_xz_xz, g_xz_0_xz_yy, g_xz_0_xz_yz, g_xz_0_xz_zz, g_xz_yy_xz_xx, g_xz_yy_xz_xy, g_xz_yy_xz_xz, g_xz_yy_xz_yy, g_xz_yy_xz_yz, g_xz_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_y_xz_xx[i] = -2.0 * g_xz_0_xz_xx[i] * a_exp + 4.0 * g_xz_yy_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_xz_xy[i] = -2.0 * g_xz_0_xz_xy[i] * a_exp + 4.0 * g_xz_yy_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_xz_xz[i] = -2.0 * g_xz_0_xz_xz[i] * a_exp + 4.0 * g_xz_yy_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_xz_yy[i] = -2.0 * g_xz_0_xz_yy[i] * a_exp + 4.0 * g_xz_yy_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_xz_yz[i] = -2.0 * g_xz_0_xz_yz[i] * a_exp + 4.0 * g_xz_yy_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_xz_zz[i] = -2.0 * g_xz_0_xz_zz[i] * a_exp + 4.0 * g_xz_yy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_x_y_0_0_z_y_yy_xx, g_x_y_0_0_z_y_yy_xy, g_x_y_0_0_z_y_yy_xz, g_x_y_0_0_z_y_yy_yy, g_x_y_0_0_z_y_yy_yz, g_x_y_0_0_z_y_yy_zz, g_xz_0_yy_xx, g_xz_0_yy_xy, g_xz_0_yy_xz, g_xz_0_yy_yy, g_xz_0_yy_yz, g_xz_0_yy_zz, g_xz_yy_yy_xx, g_xz_yy_yy_xy, g_xz_yy_yy_xz, g_xz_yy_yy_yy, g_xz_yy_yy_yz, g_xz_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_y_yy_xx[i] = -2.0 * g_xz_0_yy_xx[i] * a_exp + 4.0 * g_xz_yy_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_yy_xy[i] = -2.0 * g_xz_0_yy_xy[i] * a_exp + 4.0 * g_xz_yy_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_yy_xz[i] = -2.0 * g_xz_0_yy_xz[i] * a_exp + 4.0 * g_xz_yy_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_yy_yy[i] = -2.0 * g_xz_0_yy_yy[i] * a_exp + 4.0 * g_xz_yy_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_yy_yz[i] = -2.0 * g_xz_0_yy_yz[i] * a_exp + 4.0 * g_xz_yy_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_yy_zz[i] = -2.0 * g_xz_0_yy_zz[i] * a_exp + 4.0 * g_xz_yy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_x_y_0_0_z_y_yz_xx, g_x_y_0_0_z_y_yz_xy, g_x_y_0_0_z_y_yz_xz, g_x_y_0_0_z_y_yz_yy, g_x_y_0_0_z_y_yz_yz, g_x_y_0_0_z_y_yz_zz, g_xz_0_yz_xx, g_xz_0_yz_xy, g_xz_0_yz_xz, g_xz_0_yz_yy, g_xz_0_yz_yz, g_xz_0_yz_zz, g_xz_yy_yz_xx, g_xz_yy_yz_xy, g_xz_yy_yz_xz, g_xz_yy_yz_yy, g_xz_yy_yz_yz, g_xz_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_y_yz_xx[i] = -2.0 * g_xz_0_yz_xx[i] * a_exp + 4.0 * g_xz_yy_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_yz_xy[i] = -2.0 * g_xz_0_yz_xy[i] * a_exp + 4.0 * g_xz_yy_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_yz_xz[i] = -2.0 * g_xz_0_yz_xz[i] * a_exp + 4.0 * g_xz_yy_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_yz_yy[i] = -2.0 * g_xz_0_yz_yy[i] * a_exp + 4.0 * g_xz_yy_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_yz_yz[i] = -2.0 * g_xz_0_yz_yz[i] * a_exp + 4.0 * g_xz_yy_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_yz_zz[i] = -2.0 * g_xz_0_yz_zz[i] * a_exp + 4.0 * g_xz_yy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_x_y_0_0_z_y_zz_xx, g_x_y_0_0_z_y_zz_xy, g_x_y_0_0_z_y_zz_xz, g_x_y_0_0_z_y_zz_yy, g_x_y_0_0_z_y_zz_yz, g_x_y_0_0_z_y_zz_zz, g_xz_0_zz_xx, g_xz_0_zz_xy, g_xz_0_zz_xz, g_xz_0_zz_yy, g_xz_0_zz_yz, g_xz_0_zz_zz, g_xz_yy_zz_xx, g_xz_yy_zz_xy, g_xz_yy_zz_xz, g_xz_yy_zz_yy, g_xz_yy_zz_yz, g_xz_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_y_zz_xx[i] = -2.0 * g_xz_0_zz_xx[i] * a_exp + 4.0 * g_xz_yy_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_zz_xy[i] = -2.0 * g_xz_0_zz_xy[i] * a_exp + 4.0 * g_xz_yy_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_zz_xz[i] = -2.0 * g_xz_0_zz_xz[i] * a_exp + 4.0 * g_xz_yy_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_zz_yy[i] = -2.0 * g_xz_0_zz_yy[i] * a_exp + 4.0 * g_xz_yy_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_zz_yz[i] = -2.0 * g_xz_0_zz_yz[i] * a_exp + 4.0 * g_xz_yy_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_zz_zz[i] = -2.0 * g_xz_0_zz_zz[i] * a_exp + 4.0 * g_xz_yy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_x_y_0_0_z_z_xx_xx, g_x_y_0_0_z_z_xx_xy, g_x_y_0_0_z_z_xx_xz, g_x_y_0_0_z_z_xx_yy, g_x_y_0_0_z_z_xx_yz, g_x_y_0_0_z_z_xx_zz, g_xz_yz_xx_xx, g_xz_yz_xx_xy, g_xz_yz_xx_xz, g_xz_yz_xx_yy, g_xz_yz_xx_yz, g_xz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_z_xx_xx[i] = 4.0 * g_xz_yz_xx_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_xx_xy[i] = 4.0 * g_xz_yz_xx_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_xx_xz[i] = 4.0 * g_xz_yz_xx_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_xx_yy[i] = 4.0 * g_xz_yz_xx_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_xx_yz[i] = 4.0 * g_xz_yz_xx_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_xx_zz[i] = 4.0 * g_xz_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_x_y_0_0_z_z_xy_xx, g_x_y_0_0_z_z_xy_xy, g_x_y_0_0_z_z_xy_xz, g_x_y_0_0_z_z_xy_yy, g_x_y_0_0_z_z_xy_yz, g_x_y_0_0_z_z_xy_zz, g_xz_yz_xy_xx, g_xz_yz_xy_xy, g_xz_yz_xy_xz, g_xz_yz_xy_yy, g_xz_yz_xy_yz, g_xz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_z_xy_xx[i] = 4.0 * g_xz_yz_xy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_xy_xy[i] = 4.0 * g_xz_yz_xy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_xy_xz[i] = 4.0 * g_xz_yz_xy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_xy_yy[i] = 4.0 * g_xz_yz_xy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_xy_yz[i] = 4.0 * g_xz_yz_xy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_xy_zz[i] = 4.0 * g_xz_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_x_y_0_0_z_z_xz_xx, g_x_y_0_0_z_z_xz_xy, g_x_y_0_0_z_z_xz_xz, g_x_y_0_0_z_z_xz_yy, g_x_y_0_0_z_z_xz_yz, g_x_y_0_0_z_z_xz_zz, g_xz_yz_xz_xx, g_xz_yz_xz_xy, g_xz_yz_xz_xz, g_xz_yz_xz_yy, g_xz_yz_xz_yz, g_xz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_z_xz_xx[i] = 4.0 * g_xz_yz_xz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_xz_xy[i] = 4.0 * g_xz_yz_xz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_xz_xz[i] = 4.0 * g_xz_yz_xz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_xz_yy[i] = 4.0 * g_xz_yz_xz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_xz_yz[i] = 4.0 * g_xz_yz_xz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_xz_zz[i] = 4.0 * g_xz_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_x_y_0_0_z_z_yy_xx, g_x_y_0_0_z_z_yy_xy, g_x_y_0_0_z_z_yy_xz, g_x_y_0_0_z_z_yy_yy, g_x_y_0_0_z_z_yy_yz, g_x_y_0_0_z_z_yy_zz, g_xz_yz_yy_xx, g_xz_yz_yy_xy, g_xz_yz_yy_xz, g_xz_yz_yy_yy, g_xz_yz_yy_yz, g_xz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_z_yy_xx[i] = 4.0 * g_xz_yz_yy_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_yy_xy[i] = 4.0 * g_xz_yz_yy_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_yy_xz[i] = 4.0 * g_xz_yz_yy_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_yy_yy[i] = 4.0 * g_xz_yz_yy_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_yy_yz[i] = 4.0 * g_xz_yz_yy_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_yy_zz[i] = 4.0 * g_xz_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_x_y_0_0_z_z_yz_xx, g_x_y_0_0_z_z_yz_xy, g_x_y_0_0_z_z_yz_xz, g_x_y_0_0_z_z_yz_yy, g_x_y_0_0_z_z_yz_yz, g_x_y_0_0_z_z_yz_zz, g_xz_yz_yz_xx, g_xz_yz_yz_xy, g_xz_yz_yz_xz, g_xz_yz_yz_yy, g_xz_yz_yz_yz, g_xz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_z_yz_xx[i] = 4.0 * g_xz_yz_yz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_yz_xy[i] = 4.0 * g_xz_yz_yz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_yz_xz[i] = 4.0 * g_xz_yz_yz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_yz_yy[i] = 4.0 * g_xz_yz_yz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_yz_yz[i] = 4.0 * g_xz_yz_yz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_yz_zz[i] = 4.0 * g_xz_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_x_y_0_0_z_z_zz_xx, g_x_y_0_0_z_z_zz_xy, g_x_y_0_0_z_z_zz_xz, g_x_y_0_0_z_z_zz_yy, g_x_y_0_0_z_z_zz_yz, g_x_y_0_0_z_z_zz_zz, g_xz_yz_zz_xx, g_xz_yz_zz_xy, g_xz_yz_zz_xz, g_xz_yz_zz_yy, g_xz_yz_zz_yz, g_xz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_z_zz_xx[i] = 4.0 * g_xz_yz_zz_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_zz_xy[i] = 4.0 * g_xz_yz_zz_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_zz_xz[i] = 4.0 * g_xz_yz_zz_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_zz_yy[i] = 4.0 * g_xz_yz_zz_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_zz_yz[i] = 4.0 * g_xz_yz_zz_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_zz_zz[i] = 4.0 * g_xz_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_0_xz_xx_xx, g_0_xz_xx_xy, g_0_xz_xx_xz, g_0_xz_xx_yy, g_0_xz_xx_yz, g_0_xz_xx_zz, g_x_z_0_0_x_x_xx_xx, g_x_z_0_0_x_x_xx_xy, g_x_z_0_0_x_x_xx_xz, g_x_z_0_0_x_x_xx_yy, g_x_z_0_0_x_x_xx_yz, g_x_z_0_0_x_x_xx_zz, g_xx_xz_xx_xx, g_xx_xz_xx_xy, g_xx_xz_xx_xz, g_xx_xz_xx_yy, g_xx_xz_xx_yz, g_xx_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_x_xx_xx[i] = -2.0 * g_0_xz_xx_xx[i] * b_exp + 4.0 * g_xx_xz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_xx_xy[i] = -2.0 * g_0_xz_xx_xy[i] * b_exp + 4.0 * g_xx_xz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_xx_xz[i] = -2.0 * g_0_xz_xx_xz[i] * b_exp + 4.0 * g_xx_xz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_xx_yy[i] = -2.0 * g_0_xz_xx_yy[i] * b_exp + 4.0 * g_xx_xz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_xx_yz[i] = -2.0 * g_0_xz_xx_yz[i] * b_exp + 4.0 * g_xx_xz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_xx_zz[i] = -2.0 * g_0_xz_xx_zz[i] * b_exp + 4.0 * g_xx_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_0_xz_xy_xx, g_0_xz_xy_xy, g_0_xz_xy_xz, g_0_xz_xy_yy, g_0_xz_xy_yz, g_0_xz_xy_zz, g_x_z_0_0_x_x_xy_xx, g_x_z_0_0_x_x_xy_xy, g_x_z_0_0_x_x_xy_xz, g_x_z_0_0_x_x_xy_yy, g_x_z_0_0_x_x_xy_yz, g_x_z_0_0_x_x_xy_zz, g_xx_xz_xy_xx, g_xx_xz_xy_xy, g_xx_xz_xy_xz, g_xx_xz_xy_yy, g_xx_xz_xy_yz, g_xx_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_x_xy_xx[i] = -2.0 * g_0_xz_xy_xx[i] * b_exp + 4.0 * g_xx_xz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_xy_xy[i] = -2.0 * g_0_xz_xy_xy[i] * b_exp + 4.0 * g_xx_xz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_xy_xz[i] = -2.0 * g_0_xz_xy_xz[i] * b_exp + 4.0 * g_xx_xz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_xy_yy[i] = -2.0 * g_0_xz_xy_yy[i] * b_exp + 4.0 * g_xx_xz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_xy_yz[i] = -2.0 * g_0_xz_xy_yz[i] * b_exp + 4.0 * g_xx_xz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_xy_zz[i] = -2.0 * g_0_xz_xy_zz[i] * b_exp + 4.0 * g_xx_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_0_xz_xz_xx, g_0_xz_xz_xy, g_0_xz_xz_xz, g_0_xz_xz_yy, g_0_xz_xz_yz, g_0_xz_xz_zz, g_x_z_0_0_x_x_xz_xx, g_x_z_0_0_x_x_xz_xy, g_x_z_0_0_x_x_xz_xz, g_x_z_0_0_x_x_xz_yy, g_x_z_0_0_x_x_xz_yz, g_x_z_0_0_x_x_xz_zz, g_xx_xz_xz_xx, g_xx_xz_xz_xy, g_xx_xz_xz_xz, g_xx_xz_xz_yy, g_xx_xz_xz_yz, g_xx_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_x_xz_xx[i] = -2.0 * g_0_xz_xz_xx[i] * b_exp + 4.0 * g_xx_xz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_xz_xy[i] = -2.0 * g_0_xz_xz_xy[i] * b_exp + 4.0 * g_xx_xz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_xz_xz[i] = -2.0 * g_0_xz_xz_xz[i] * b_exp + 4.0 * g_xx_xz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_xz_yy[i] = -2.0 * g_0_xz_xz_yy[i] * b_exp + 4.0 * g_xx_xz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_xz_yz[i] = -2.0 * g_0_xz_xz_yz[i] * b_exp + 4.0 * g_xx_xz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_xz_zz[i] = -2.0 * g_0_xz_xz_zz[i] * b_exp + 4.0 * g_xx_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_0_xz_yy_xx, g_0_xz_yy_xy, g_0_xz_yy_xz, g_0_xz_yy_yy, g_0_xz_yy_yz, g_0_xz_yy_zz, g_x_z_0_0_x_x_yy_xx, g_x_z_0_0_x_x_yy_xy, g_x_z_0_0_x_x_yy_xz, g_x_z_0_0_x_x_yy_yy, g_x_z_0_0_x_x_yy_yz, g_x_z_0_0_x_x_yy_zz, g_xx_xz_yy_xx, g_xx_xz_yy_xy, g_xx_xz_yy_xz, g_xx_xz_yy_yy, g_xx_xz_yy_yz, g_xx_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_x_yy_xx[i] = -2.0 * g_0_xz_yy_xx[i] * b_exp + 4.0 * g_xx_xz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_yy_xy[i] = -2.0 * g_0_xz_yy_xy[i] * b_exp + 4.0 * g_xx_xz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_yy_xz[i] = -2.0 * g_0_xz_yy_xz[i] * b_exp + 4.0 * g_xx_xz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_yy_yy[i] = -2.0 * g_0_xz_yy_yy[i] * b_exp + 4.0 * g_xx_xz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_yy_yz[i] = -2.0 * g_0_xz_yy_yz[i] * b_exp + 4.0 * g_xx_xz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_yy_zz[i] = -2.0 * g_0_xz_yy_zz[i] * b_exp + 4.0 * g_xx_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_0_xz_yz_xx, g_0_xz_yz_xy, g_0_xz_yz_xz, g_0_xz_yz_yy, g_0_xz_yz_yz, g_0_xz_yz_zz, g_x_z_0_0_x_x_yz_xx, g_x_z_0_0_x_x_yz_xy, g_x_z_0_0_x_x_yz_xz, g_x_z_0_0_x_x_yz_yy, g_x_z_0_0_x_x_yz_yz, g_x_z_0_0_x_x_yz_zz, g_xx_xz_yz_xx, g_xx_xz_yz_xy, g_xx_xz_yz_xz, g_xx_xz_yz_yy, g_xx_xz_yz_yz, g_xx_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_x_yz_xx[i] = -2.0 * g_0_xz_yz_xx[i] * b_exp + 4.0 * g_xx_xz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_yz_xy[i] = -2.0 * g_0_xz_yz_xy[i] * b_exp + 4.0 * g_xx_xz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_yz_xz[i] = -2.0 * g_0_xz_yz_xz[i] * b_exp + 4.0 * g_xx_xz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_yz_yy[i] = -2.0 * g_0_xz_yz_yy[i] * b_exp + 4.0 * g_xx_xz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_yz_yz[i] = -2.0 * g_0_xz_yz_yz[i] * b_exp + 4.0 * g_xx_xz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_yz_zz[i] = -2.0 * g_0_xz_yz_zz[i] * b_exp + 4.0 * g_xx_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_0_xz_zz_xx, g_0_xz_zz_xy, g_0_xz_zz_xz, g_0_xz_zz_yy, g_0_xz_zz_yz, g_0_xz_zz_zz, g_x_z_0_0_x_x_zz_xx, g_x_z_0_0_x_x_zz_xy, g_x_z_0_0_x_x_zz_xz, g_x_z_0_0_x_x_zz_yy, g_x_z_0_0_x_x_zz_yz, g_x_z_0_0_x_x_zz_zz, g_xx_xz_zz_xx, g_xx_xz_zz_xy, g_xx_xz_zz_xz, g_xx_xz_zz_yy, g_xx_xz_zz_yz, g_xx_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_x_zz_xx[i] = -2.0 * g_0_xz_zz_xx[i] * b_exp + 4.0 * g_xx_xz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_zz_xy[i] = -2.0 * g_0_xz_zz_xy[i] * b_exp + 4.0 * g_xx_xz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_zz_xz[i] = -2.0 * g_0_xz_zz_xz[i] * b_exp + 4.0 * g_xx_xz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_zz_yy[i] = -2.0 * g_0_xz_zz_yy[i] * b_exp + 4.0 * g_xx_xz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_zz_yz[i] = -2.0 * g_0_xz_zz_yz[i] * b_exp + 4.0 * g_xx_xz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_zz_zz[i] = -2.0 * g_0_xz_zz_zz[i] * b_exp + 4.0 * g_xx_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_0_yz_xx_xx, g_0_yz_xx_xy, g_0_yz_xx_xz, g_0_yz_xx_yy, g_0_yz_xx_yz, g_0_yz_xx_zz, g_x_z_0_0_x_y_xx_xx, g_x_z_0_0_x_y_xx_xy, g_x_z_0_0_x_y_xx_xz, g_x_z_0_0_x_y_xx_yy, g_x_z_0_0_x_y_xx_yz, g_x_z_0_0_x_y_xx_zz, g_xx_yz_xx_xx, g_xx_yz_xx_xy, g_xx_yz_xx_xz, g_xx_yz_xx_yy, g_xx_yz_xx_yz, g_xx_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_y_xx_xx[i] = -2.0 * g_0_yz_xx_xx[i] * b_exp + 4.0 * g_xx_yz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_xx_xy[i] = -2.0 * g_0_yz_xx_xy[i] * b_exp + 4.0 * g_xx_yz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_xx_xz[i] = -2.0 * g_0_yz_xx_xz[i] * b_exp + 4.0 * g_xx_yz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_xx_yy[i] = -2.0 * g_0_yz_xx_yy[i] * b_exp + 4.0 * g_xx_yz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_xx_yz[i] = -2.0 * g_0_yz_xx_yz[i] * b_exp + 4.0 * g_xx_yz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_xx_zz[i] = -2.0 * g_0_yz_xx_zz[i] * b_exp + 4.0 * g_xx_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_0_yz_xy_xx, g_0_yz_xy_xy, g_0_yz_xy_xz, g_0_yz_xy_yy, g_0_yz_xy_yz, g_0_yz_xy_zz, g_x_z_0_0_x_y_xy_xx, g_x_z_0_0_x_y_xy_xy, g_x_z_0_0_x_y_xy_xz, g_x_z_0_0_x_y_xy_yy, g_x_z_0_0_x_y_xy_yz, g_x_z_0_0_x_y_xy_zz, g_xx_yz_xy_xx, g_xx_yz_xy_xy, g_xx_yz_xy_xz, g_xx_yz_xy_yy, g_xx_yz_xy_yz, g_xx_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_y_xy_xx[i] = -2.0 * g_0_yz_xy_xx[i] * b_exp + 4.0 * g_xx_yz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_xy_xy[i] = -2.0 * g_0_yz_xy_xy[i] * b_exp + 4.0 * g_xx_yz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_xy_xz[i] = -2.0 * g_0_yz_xy_xz[i] * b_exp + 4.0 * g_xx_yz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_xy_yy[i] = -2.0 * g_0_yz_xy_yy[i] * b_exp + 4.0 * g_xx_yz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_xy_yz[i] = -2.0 * g_0_yz_xy_yz[i] * b_exp + 4.0 * g_xx_yz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_xy_zz[i] = -2.0 * g_0_yz_xy_zz[i] * b_exp + 4.0 * g_xx_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_0_yz_xz_xx, g_0_yz_xz_xy, g_0_yz_xz_xz, g_0_yz_xz_yy, g_0_yz_xz_yz, g_0_yz_xz_zz, g_x_z_0_0_x_y_xz_xx, g_x_z_0_0_x_y_xz_xy, g_x_z_0_0_x_y_xz_xz, g_x_z_0_0_x_y_xz_yy, g_x_z_0_0_x_y_xz_yz, g_x_z_0_0_x_y_xz_zz, g_xx_yz_xz_xx, g_xx_yz_xz_xy, g_xx_yz_xz_xz, g_xx_yz_xz_yy, g_xx_yz_xz_yz, g_xx_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_y_xz_xx[i] = -2.0 * g_0_yz_xz_xx[i] * b_exp + 4.0 * g_xx_yz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_xz_xy[i] = -2.0 * g_0_yz_xz_xy[i] * b_exp + 4.0 * g_xx_yz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_xz_xz[i] = -2.0 * g_0_yz_xz_xz[i] * b_exp + 4.0 * g_xx_yz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_xz_yy[i] = -2.0 * g_0_yz_xz_yy[i] * b_exp + 4.0 * g_xx_yz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_xz_yz[i] = -2.0 * g_0_yz_xz_yz[i] * b_exp + 4.0 * g_xx_yz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_xz_zz[i] = -2.0 * g_0_yz_xz_zz[i] * b_exp + 4.0 * g_xx_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_0_yz_yy_xx, g_0_yz_yy_xy, g_0_yz_yy_xz, g_0_yz_yy_yy, g_0_yz_yy_yz, g_0_yz_yy_zz, g_x_z_0_0_x_y_yy_xx, g_x_z_0_0_x_y_yy_xy, g_x_z_0_0_x_y_yy_xz, g_x_z_0_0_x_y_yy_yy, g_x_z_0_0_x_y_yy_yz, g_x_z_0_0_x_y_yy_zz, g_xx_yz_yy_xx, g_xx_yz_yy_xy, g_xx_yz_yy_xz, g_xx_yz_yy_yy, g_xx_yz_yy_yz, g_xx_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_y_yy_xx[i] = -2.0 * g_0_yz_yy_xx[i] * b_exp + 4.0 * g_xx_yz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_yy_xy[i] = -2.0 * g_0_yz_yy_xy[i] * b_exp + 4.0 * g_xx_yz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_yy_xz[i] = -2.0 * g_0_yz_yy_xz[i] * b_exp + 4.0 * g_xx_yz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_yy_yy[i] = -2.0 * g_0_yz_yy_yy[i] * b_exp + 4.0 * g_xx_yz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_yy_yz[i] = -2.0 * g_0_yz_yy_yz[i] * b_exp + 4.0 * g_xx_yz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_yy_zz[i] = -2.0 * g_0_yz_yy_zz[i] * b_exp + 4.0 * g_xx_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_0_yz_yz_xx, g_0_yz_yz_xy, g_0_yz_yz_xz, g_0_yz_yz_yy, g_0_yz_yz_yz, g_0_yz_yz_zz, g_x_z_0_0_x_y_yz_xx, g_x_z_0_0_x_y_yz_xy, g_x_z_0_0_x_y_yz_xz, g_x_z_0_0_x_y_yz_yy, g_x_z_0_0_x_y_yz_yz, g_x_z_0_0_x_y_yz_zz, g_xx_yz_yz_xx, g_xx_yz_yz_xy, g_xx_yz_yz_xz, g_xx_yz_yz_yy, g_xx_yz_yz_yz, g_xx_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_y_yz_xx[i] = -2.0 * g_0_yz_yz_xx[i] * b_exp + 4.0 * g_xx_yz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_yz_xy[i] = -2.0 * g_0_yz_yz_xy[i] * b_exp + 4.0 * g_xx_yz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_yz_xz[i] = -2.0 * g_0_yz_yz_xz[i] * b_exp + 4.0 * g_xx_yz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_yz_yy[i] = -2.0 * g_0_yz_yz_yy[i] * b_exp + 4.0 * g_xx_yz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_yz_yz[i] = -2.0 * g_0_yz_yz_yz[i] * b_exp + 4.0 * g_xx_yz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_yz_zz[i] = -2.0 * g_0_yz_yz_zz[i] * b_exp + 4.0 * g_xx_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_0_yz_zz_xx, g_0_yz_zz_xy, g_0_yz_zz_xz, g_0_yz_zz_yy, g_0_yz_zz_yz, g_0_yz_zz_zz, g_x_z_0_0_x_y_zz_xx, g_x_z_0_0_x_y_zz_xy, g_x_z_0_0_x_y_zz_xz, g_x_z_0_0_x_y_zz_yy, g_x_z_0_0_x_y_zz_yz, g_x_z_0_0_x_y_zz_zz, g_xx_yz_zz_xx, g_xx_yz_zz_xy, g_xx_yz_zz_xz, g_xx_yz_zz_yy, g_xx_yz_zz_yz, g_xx_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_y_zz_xx[i] = -2.0 * g_0_yz_zz_xx[i] * b_exp + 4.0 * g_xx_yz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_zz_xy[i] = -2.0 * g_0_yz_zz_xy[i] * b_exp + 4.0 * g_xx_yz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_zz_xz[i] = -2.0 * g_0_yz_zz_xz[i] * b_exp + 4.0 * g_xx_yz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_zz_yy[i] = -2.0 * g_0_yz_zz_yy[i] * b_exp + 4.0 * g_xx_yz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_zz_yz[i] = -2.0 * g_0_yz_zz_yz[i] * b_exp + 4.0 * g_xx_yz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_zz_zz[i] = -2.0 * g_0_yz_zz_zz[i] * b_exp + 4.0 * g_xx_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_0_0_xx_xx, g_0_0_xx_xy, g_0_0_xx_xz, g_0_0_xx_yy, g_0_0_xx_yz, g_0_0_xx_zz, g_0_zz_xx_xx, g_0_zz_xx_xy, g_0_zz_xx_xz, g_0_zz_xx_yy, g_0_zz_xx_yz, g_0_zz_xx_zz, g_x_z_0_0_x_z_xx_xx, g_x_z_0_0_x_z_xx_xy, g_x_z_0_0_x_z_xx_xz, g_x_z_0_0_x_z_xx_yy, g_x_z_0_0_x_z_xx_yz, g_x_z_0_0_x_z_xx_zz, g_xx_0_xx_xx, g_xx_0_xx_xy, g_xx_0_xx_xz, g_xx_0_xx_yy, g_xx_0_xx_yz, g_xx_0_xx_zz, g_xx_zz_xx_xx, g_xx_zz_xx_xy, g_xx_zz_xx_xz, g_xx_zz_xx_yy, g_xx_zz_xx_yz, g_xx_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_z_xx_xx[i] = g_0_0_xx_xx[i] - 2.0 * g_0_zz_xx_xx[i] * b_exp - 2.0 * g_xx_0_xx_xx[i] * a_exp + 4.0 * g_xx_zz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_xx_xy[i] = g_0_0_xx_xy[i] - 2.0 * g_0_zz_xx_xy[i] * b_exp - 2.0 * g_xx_0_xx_xy[i] * a_exp + 4.0 * g_xx_zz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_xx_xz[i] = g_0_0_xx_xz[i] - 2.0 * g_0_zz_xx_xz[i] * b_exp - 2.0 * g_xx_0_xx_xz[i] * a_exp + 4.0 * g_xx_zz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_xx_yy[i] = g_0_0_xx_yy[i] - 2.0 * g_0_zz_xx_yy[i] * b_exp - 2.0 * g_xx_0_xx_yy[i] * a_exp + 4.0 * g_xx_zz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_xx_yz[i] = g_0_0_xx_yz[i] - 2.0 * g_0_zz_xx_yz[i] * b_exp - 2.0 * g_xx_0_xx_yz[i] * a_exp + 4.0 * g_xx_zz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_xx_zz[i] = g_0_0_xx_zz[i] - 2.0 * g_0_zz_xx_zz[i] * b_exp - 2.0 * g_xx_0_xx_zz[i] * a_exp + 4.0 * g_xx_zz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_0_0_xy_xx, g_0_0_xy_xy, g_0_0_xy_xz, g_0_0_xy_yy, g_0_0_xy_yz, g_0_0_xy_zz, g_0_zz_xy_xx, g_0_zz_xy_xy, g_0_zz_xy_xz, g_0_zz_xy_yy, g_0_zz_xy_yz, g_0_zz_xy_zz, g_x_z_0_0_x_z_xy_xx, g_x_z_0_0_x_z_xy_xy, g_x_z_0_0_x_z_xy_xz, g_x_z_0_0_x_z_xy_yy, g_x_z_0_0_x_z_xy_yz, g_x_z_0_0_x_z_xy_zz, g_xx_0_xy_xx, g_xx_0_xy_xy, g_xx_0_xy_xz, g_xx_0_xy_yy, g_xx_0_xy_yz, g_xx_0_xy_zz, g_xx_zz_xy_xx, g_xx_zz_xy_xy, g_xx_zz_xy_xz, g_xx_zz_xy_yy, g_xx_zz_xy_yz, g_xx_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_z_xy_xx[i] = g_0_0_xy_xx[i] - 2.0 * g_0_zz_xy_xx[i] * b_exp - 2.0 * g_xx_0_xy_xx[i] * a_exp + 4.0 * g_xx_zz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_xy_xy[i] = g_0_0_xy_xy[i] - 2.0 * g_0_zz_xy_xy[i] * b_exp - 2.0 * g_xx_0_xy_xy[i] * a_exp + 4.0 * g_xx_zz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_xy_xz[i] = g_0_0_xy_xz[i] - 2.0 * g_0_zz_xy_xz[i] * b_exp - 2.0 * g_xx_0_xy_xz[i] * a_exp + 4.0 * g_xx_zz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_xy_yy[i] = g_0_0_xy_yy[i] - 2.0 * g_0_zz_xy_yy[i] * b_exp - 2.0 * g_xx_0_xy_yy[i] * a_exp + 4.0 * g_xx_zz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_xy_yz[i] = g_0_0_xy_yz[i] - 2.0 * g_0_zz_xy_yz[i] * b_exp - 2.0 * g_xx_0_xy_yz[i] * a_exp + 4.0 * g_xx_zz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_xy_zz[i] = g_0_0_xy_zz[i] - 2.0 * g_0_zz_xy_zz[i] * b_exp - 2.0 * g_xx_0_xy_zz[i] * a_exp + 4.0 * g_xx_zz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_0_0_xz_xx, g_0_0_xz_xy, g_0_0_xz_xz, g_0_0_xz_yy, g_0_0_xz_yz, g_0_0_xz_zz, g_0_zz_xz_xx, g_0_zz_xz_xy, g_0_zz_xz_xz, g_0_zz_xz_yy, g_0_zz_xz_yz, g_0_zz_xz_zz, g_x_z_0_0_x_z_xz_xx, g_x_z_0_0_x_z_xz_xy, g_x_z_0_0_x_z_xz_xz, g_x_z_0_0_x_z_xz_yy, g_x_z_0_0_x_z_xz_yz, g_x_z_0_0_x_z_xz_zz, g_xx_0_xz_xx, g_xx_0_xz_xy, g_xx_0_xz_xz, g_xx_0_xz_yy, g_xx_0_xz_yz, g_xx_0_xz_zz, g_xx_zz_xz_xx, g_xx_zz_xz_xy, g_xx_zz_xz_xz, g_xx_zz_xz_yy, g_xx_zz_xz_yz, g_xx_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_z_xz_xx[i] = g_0_0_xz_xx[i] - 2.0 * g_0_zz_xz_xx[i] * b_exp - 2.0 * g_xx_0_xz_xx[i] * a_exp + 4.0 * g_xx_zz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_xz_xy[i] = g_0_0_xz_xy[i] - 2.0 * g_0_zz_xz_xy[i] * b_exp - 2.0 * g_xx_0_xz_xy[i] * a_exp + 4.0 * g_xx_zz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_xz_xz[i] = g_0_0_xz_xz[i] - 2.0 * g_0_zz_xz_xz[i] * b_exp - 2.0 * g_xx_0_xz_xz[i] * a_exp + 4.0 * g_xx_zz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_xz_yy[i] = g_0_0_xz_yy[i] - 2.0 * g_0_zz_xz_yy[i] * b_exp - 2.0 * g_xx_0_xz_yy[i] * a_exp + 4.0 * g_xx_zz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_xz_yz[i] = g_0_0_xz_yz[i] - 2.0 * g_0_zz_xz_yz[i] * b_exp - 2.0 * g_xx_0_xz_yz[i] * a_exp + 4.0 * g_xx_zz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_xz_zz[i] = g_0_0_xz_zz[i] - 2.0 * g_0_zz_xz_zz[i] * b_exp - 2.0 * g_xx_0_xz_zz[i] * a_exp + 4.0 * g_xx_zz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_0_0_yy_xx, g_0_0_yy_xy, g_0_0_yy_xz, g_0_0_yy_yy, g_0_0_yy_yz, g_0_0_yy_zz, g_0_zz_yy_xx, g_0_zz_yy_xy, g_0_zz_yy_xz, g_0_zz_yy_yy, g_0_zz_yy_yz, g_0_zz_yy_zz, g_x_z_0_0_x_z_yy_xx, g_x_z_0_0_x_z_yy_xy, g_x_z_0_0_x_z_yy_xz, g_x_z_0_0_x_z_yy_yy, g_x_z_0_0_x_z_yy_yz, g_x_z_0_0_x_z_yy_zz, g_xx_0_yy_xx, g_xx_0_yy_xy, g_xx_0_yy_xz, g_xx_0_yy_yy, g_xx_0_yy_yz, g_xx_0_yy_zz, g_xx_zz_yy_xx, g_xx_zz_yy_xy, g_xx_zz_yy_xz, g_xx_zz_yy_yy, g_xx_zz_yy_yz, g_xx_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_z_yy_xx[i] = g_0_0_yy_xx[i] - 2.0 * g_0_zz_yy_xx[i] * b_exp - 2.0 * g_xx_0_yy_xx[i] * a_exp + 4.0 * g_xx_zz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_yy_xy[i] = g_0_0_yy_xy[i] - 2.0 * g_0_zz_yy_xy[i] * b_exp - 2.0 * g_xx_0_yy_xy[i] * a_exp + 4.0 * g_xx_zz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_yy_xz[i] = g_0_0_yy_xz[i] - 2.0 * g_0_zz_yy_xz[i] * b_exp - 2.0 * g_xx_0_yy_xz[i] * a_exp + 4.0 * g_xx_zz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_yy_yy[i] = g_0_0_yy_yy[i] - 2.0 * g_0_zz_yy_yy[i] * b_exp - 2.0 * g_xx_0_yy_yy[i] * a_exp + 4.0 * g_xx_zz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_yy_yz[i] = g_0_0_yy_yz[i] - 2.0 * g_0_zz_yy_yz[i] * b_exp - 2.0 * g_xx_0_yy_yz[i] * a_exp + 4.0 * g_xx_zz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_yy_zz[i] = g_0_0_yy_zz[i] - 2.0 * g_0_zz_yy_zz[i] * b_exp - 2.0 * g_xx_0_yy_zz[i] * a_exp + 4.0 * g_xx_zz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_0_0_yz_xx, g_0_0_yz_xy, g_0_0_yz_xz, g_0_0_yz_yy, g_0_0_yz_yz, g_0_0_yz_zz, g_0_zz_yz_xx, g_0_zz_yz_xy, g_0_zz_yz_xz, g_0_zz_yz_yy, g_0_zz_yz_yz, g_0_zz_yz_zz, g_x_z_0_0_x_z_yz_xx, g_x_z_0_0_x_z_yz_xy, g_x_z_0_0_x_z_yz_xz, g_x_z_0_0_x_z_yz_yy, g_x_z_0_0_x_z_yz_yz, g_x_z_0_0_x_z_yz_zz, g_xx_0_yz_xx, g_xx_0_yz_xy, g_xx_0_yz_xz, g_xx_0_yz_yy, g_xx_0_yz_yz, g_xx_0_yz_zz, g_xx_zz_yz_xx, g_xx_zz_yz_xy, g_xx_zz_yz_xz, g_xx_zz_yz_yy, g_xx_zz_yz_yz, g_xx_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_z_yz_xx[i] = g_0_0_yz_xx[i] - 2.0 * g_0_zz_yz_xx[i] * b_exp - 2.0 * g_xx_0_yz_xx[i] * a_exp + 4.0 * g_xx_zz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_yz_xy[i] = g_0_0_yz_xy[i] - 2.0 * g_0_zz_yz_xy[i] * b_exp - 2.0 * g_xx_0_yz_xy[i] * a_exp + 4.0 * g_xx_zz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_yz_xz[i] = g_0_0_yz_xz[i] - 2.0 * g_0_zz_yz_xz[i] * b_exp - 2.0 * g_xx_0_yz_xz[i] * a_exp + 4.0 * g_xx_zz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_yz_yy[i] = g_0_0_yz_yy[i] - 2.0 * g_0_zz_yz_yy[i] * b_exp - 2.0 * g_xx_0_yz_yy[i] * a_exp + 4.0 * g_xx_zz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_yz_yz[i] = g_0_0_yz_yz[i] - 2.0 * g_0_zz_yz_yz[i] * b_exp - 2.0 * g_xx_0_yz_yz[i] * a_exp + 4.0 * g_xx_zz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_yz_zz[i] = g_0_0_yz_zz[i] - 2.0 * g_0_zz_yz_zz[i] * b_exp - 2.0 * g_xx_0_yz_zz[i] * a_exp + 4.0 * g_xx_zz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_0_0_zz_xx, g_0_0_zz_xy, g_0_0_zz_xz, g_0_0_zz_yy, g_0_0_zz_yz, g_0_0_zz_zz, g_0_zz_zz_xx, g_0_zz_zz_xy, g_0_zz_zz_xz, g_0_zz_zz_yy, g_0_zz_zz_yz, g_0_zz_zz_zz, g_x_z_0_0_x_z_zz_xx, g_x_z_0_0_x_z_zz_xy, g_x_z_0_0_x_z_zz_xz, g_x_z_0_0_x_z_zz_yy, g_x_z_0_0_x_z_zz_yz, g_x_z_0_0_x_z_zz_zz, g_xx_0_zz_xx, g_xx_0_zz_xy, g_xx_0_zz_xz, g_xx_0_zz_yy, g_xx_0_zz_yz, g_xx_0_zz_zz, g_xx_zz_zz_xx, g_xx_zz_zz_xy, g_xx_zz_zz_xz, g_xx_zz_zz_yy, g_xx_zz_zz_yz, g_xx_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_z_zz_xx[i] = g_0_0_zz_xx[i] - 2.0 * g_0_zz_zz_xx[i] * b_exp - 2.0 * g_xx_0_zz_xx[i] * a_exp + 4.0 * g_xx_zz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_zz_xy[i] = g_0_0_zz_xy[i] - 2.0 * g_0_zz_zz_xy[i] * b_exp - 2.0 * g_xx_0_zz_xy[i] * a_exp + 4.0 * g_xx_zz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_zz_xz[i] = g_0_0_zz_xz[i] - 2.0 * g_0_zz_zz_xz[i] * b_exp - 2.0 * g_xx_0_zz_xz[i] * a_exp + 4.0 * g_xx_zz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_zz_yy[i] = g_0_0_zz_yy[i] - 2.0 * g_0_zz_zz_yy[i] * b_exp - 2.0 * g_xx_0_zz_yy[i] * a_exp + 4.0 * g_xx_zz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_zz_yz[i] = g_0_0_zz_yz[i] - 2.0 * g_0_zz_zz_yz[i] * b_exp - 2.0 * g_xx_0_zz_yz[i] * a_exp + 4.0 * g_xx_zz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_zz_zz[i] = g_0_0_zz_zz[i] - 2.0 * g_0_zz_zz_zz[i] * b_exp - 2.0 * g_xx_0_zz_zz[i] * a_exp + 4.0 * g_xx_zz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_x_z_0_0_y_x_xx_xx, g_x_z_0_0_y_x_xx_xy, g_x_z_0_0_y_x_xx_xz, g_x_z_0_0_y_x_xx_yy, g_x_z_0_0_y_x_xx_yz, g_x_z_0_0_y_x_xx_zz, g_xy_xz_xx_xx, g_xy_xz_xx_xy, g_xy_xz_xx_xz, g_xy_xz_xx_yy, g_xy_xz_xx_yz, g_xy_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_x_xx_xx[i] = 4.0 * g_xy_xz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_xx_xy[i] = 4.0 * g_xy_xz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_xx_xz[i] = 4.0 * g_xy_xz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_xx_yy[i] = 4.0 * g_xy_xz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_xx_yz[i] = 4.0 * g_xy_xz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_xx_zz[i] = 4.0 * g_xy_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_x_z_0_0_y_x_xy_xx, g_x_z_0_0_y_x_xy_xy, g_x_z_0_0_y_x_xy_xz, g_x_z_0_0_y_x_xy_yy, g_x_z_0_0_y_x_xy_yz, g_x_z_0_0_y_x_xy_zz, g_xy_xz_xy_xx, g_xy_xz_xy_xy, g_xy_xz_xy_xz, g_xy_xz_xy_yy, g_xy_xz_xy_yz, g_xy_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_x_xy_xx[i] = 4.0 * g_xy_xz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_xy_xy[i] = 4.0 * g_xy_xz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_xy_xz[i] = 4.0 * g_xy_xz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_xy_yy[i] = 4.0 * g_xy_xz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_xy_yz[i] = 4.0 * g_xy_xz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_xy_zz[i] = 4.0 * g_xy_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_x_z_0_0_y_x_xz_xx, g_x_z_0_0_y_x_xz_xy, g_x_z_0_0_y_x_xz_xz, g_x_z_0_0_y_x_xz_yy, g_x_z_0_0_y_x_xz_yz, g_x_z_0_0_y_x_xz_zz, g_xy_xz_xz_xx, g_xy_xz_xz_xy, g_xy_xz_xz_xz, g_xy_xz_xz_yy, g_xy_xz_xz_yz, g_xy_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_x_xz_xx[i] = 4.0 * g_xy_xz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_xz_xy[i] = 4.0 * g_xy_xz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_xz_xz[i] = 4.0 * g_xy_xz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_xz_yy[i] = 4.0 * g_xy_xz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_xz_yz[i] = 4.0 * g_xy_xz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_xz_zz[i] = 4.0 * g_xy_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_x_z_0_0_y_x_yy_xx, g_x_z_0_0_y_x_yy_xy, g_x_z_0_0_y_x_yy_xz, g_x_z_0_0_y_x_yy_yy, g_x_z_0_0_y_x_yy_yz, g_x_z_0_0_y_x_yy_zz, g_xy_xz_yy_xx, g_xy_xz_yy_xy, g_xy_xz_yy_xz, g_xy_xz_yy_yy, g_xy_xz_yy_yz, g_xy_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_x_yy_xx[i] = 4.0 * g_xy_xz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_yy_xy[i] = 4.0 * g_xy_xz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_yy_xz[i] = 4.0 * g_xy_xz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_yy_yy[i] = 4.0 * g_xy_xz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_yy_yz[i] = 4.0 * g_xy_xz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_yy_zz[i] = 4.0 * g_xy_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_x_z_0_0_y_x_yz_xx, g_x_z_0_0_y_x_yz_xy, g_x_z_0_0_y_x_yz_xz, g_x_z_0_0_y_x_yz_yy, g_x_z_0_0_y_x_yz_yz, g_x_z_0_0_y_x_yz_zz, g_xy_xz_yz_xx, g_xy_xz_yz_xy, g_xy_xz_yz_xz, g_xy_xz_yz_yy, g_xy_xz_yz_yz, g_xy_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_x_yz_xx[i] = 4.0 * g_xy_xz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_yz_xy[i] = 4.0 * g_xy_xz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_yz_xz[i] = 4.0 * g_xy_xz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_yz_yy[i] = 4.0 * g_xy_xz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_yz_yz[i] = 4.0 * g_xy_xz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_yz_zz[i] = 4.0 * g_xy_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_x_z_0_0_y_x_zz_xx, g_x_z_0_0_y_x_zz_xy, g_x_z_0_0_y_x_zz_xz, g_x_z_0_0_y_x_zz_yy, g_x_z_0_0_y_x_zz_yz, g_x_z_0_0_y_x_zz_zz, g_xy_xz_zz_xx, g_xy_xz_zz_xy, g_xy_xz_zz_xz, g_xy_xz_zz_yy, g_xy_xz_zz_yz, g_xy_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_x_zz_xx[i] = 4.0 * g_xy_xz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_zz_xy[i] = 4.0 * g_xy_xz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_zz_xz[i] = 4.0 * g_xy_xz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_zz_yy[i] = 4.0 * g_xy_xz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_zz_yz[i] = 4.0 * g_xy_xz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_zz_zz[i] = 4.0 * g_xy_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_x_z_0_0_y_y_xx_xx, g_x_z_0_0_y_y_xx_xy, g_x_z_0_0_y_y_xx_xz, g_x_z_0_0_y_y_xx_yy, g_x_z_0_0_y_y_xx_yz, g_x_z_0_0_y_y_xx_zz, g_xy_yz_xx_xx, g_xy_yz_xx_xy, g_xy_yz_xx_xz, g_xy_yz_xx_yy, g_xy_yz_xx_yz, g_xy_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_y_xx_xx[i] = 4.0 * g_xy_yz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_xx_xy[i] = 4.0 * g_xy_yz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_xx_xz[i] = 4.0 * g_xy_yz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_xx_yy[i] = 4.0 * g_xy_yz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_xx_yz[i] = 4.0 * g_xy_yz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_xx_zz[i] = 4.0 * g_xy_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_x_z_0_0_y_y_xy_xx, g_x_z_0_0_y_y_xy_xy, g_x_z_0_0_y_y_xy_xz, g_x_z_0_0_y_y_xy_yy, g_x_z_0_0_y_y_xy_yz, g_x_z_0_0_y_y_xy_zz, g_xy_yz_xy_xx, g_xy_yz_xy_xy, g_xy_yz_xy_xz, g_xy_yz_xy_yy, g_xy_yz_xy_yz, g_xy_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_y_xy_xx[i] = 4.0 * g_xy_yz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_xy_xy[i] = 4.0 * g_xy_yz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_xy_xz[i] = 4.0 * g_xy_yz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_xy_yy[i] = 4.0 * g_xy_yz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_xy_yz[i] = 4.0 * g_xy_yz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_xy_zz[i] = 4.0 * g_xy_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_x_z_0_0_y_y_xz_xx, g_x_z_0_0_y_y_xz_xy, g_x_z_0_0_y_y_xz_xz, g_x_z_0_0_y_y_xz_yy, g_x_z_0_0_y_y_xz_yz, g_x_z_0_0_y_y_xz_zz, g_xy_yz_xz_xx, g_xy_yz_xz_xy, g_xy_yz_xz_xz, g_xy_yz_xz_yy, g_xy_yz_xz_yz, g_xy_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_y_xz_xx[i] = 4.0 * g_xy_yz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_xz_xy[i] = 4.0 * g_xy_yz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_xz_xz[i] = 4.0 * g_xy_yz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_xz_yy[i] = 4.0 * g_xy_yz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_xz_yz[i] = 4.0 * g_xy_yz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_xz_zz[i] = 4.0 * g_xy_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_x_z_0_0_y_y_yy_xx, g_x_z_0_0_y_y_yy_xy, g_x_z_0_0_y_y_yy_xz, g_x_z_0_0_y_y_yy_yy, g_x_z_0_0_y_y_yy_yz, g_x_z_0_0_y_y_yy_zz, g_xy_yz_yy_xx, g_xy_yz_yy_xy, g_xy_yz_yy_xz, g_xy_yz_yy_yy, g_xy_yz_yy_yz, g_xy_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_y_yy_xx[i] = 4.0 * g_xy_yz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_yy_xy[i] = 4.0 * g_xy_yz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_yy_xz[i] = 4.0 * g_xy_yz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_yy_yy[i] = 4.0 * g_xy_yz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_yy_yz[i] = 4.0 * g_xy_yz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_yy_zz[i] = 4.0 * g_xy_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_x_z_0_0_y_y_yz_xx, g_x_z_0_0_y_y_yz_xy, g_x_z_0_0_y_y_yz_xz, g_x_z_0_0_y_y_yz_yy, g_x_z_0_0_y_y_yz_yz, g_x_z_0_0_y_y_yz_zz, g_xy_yz_yz_xx, g_xy_yz_yz_xy, g_xy_yz_yz_xz, g_xy_yz_yz_yy, g_xy_yz_yz_yz, g_xy_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_y_yz_xx[i] = 4.0 * g_xy_yz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_yz_xy[i] = 4.0 * g_xy_yz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_yz_xz[i] = 4.0 * g_xy_yz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_yz_yy[i] = 4.0 * g_xy_yz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_yz_yz[i] = 4.0 * g_xy_yz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_yz_zz[i] = 4.0 * g_xy_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_x_z_0_0_y_y_zz_xx, g_x_z_0_0_y_y_zz_xy, g_x_z_0_0_y_y_zz_xz, g_x_z_0_0_y_y_zz_yy, g_x_z_0_0_y_y_zz_yz, g_x_z_0_0_y_y_zz_zz, g_xy_yz_zz_xx, g_xy_yz_zz_xy, g_xy_yz_zz_xz, g_xy_yz_zz_yy, g_xy_yz_zz_yz, g_xy_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_y_zz_xx[i] = 4.0 * g_xy_yz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_zz_xy[i] = 4.0 * g_xy_yz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_zz_xz[i] = 4.0 * g_xy_yz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_zz_yy[i] = 4.0 * g_xy_yz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_zz_yz[i] = 4.0 * g_xy_yz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_zz_zz[i] = 4.0 * g_xy_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_x_z_0_0_y_z_xx_xx, g_x_z_0_0_y_z_xx_xy, g_x_z_0_0_y_z_xx_xz, g_x_z_0_0_y_z_xx_yy, g_x_z_0_0_y_z_xx_yz, g_x_z_0_0_y_z_xx_zz, g_xy_0_xx_xx, g_xy_0_xx_xy, g_xy_0_xx_xz, g_xy_0_xx_yy, g_xy_0_xx_yz, g_xy_0_xx_zz, g_xy_zz_xx_xx, g_xy_zz_xx_xy, g_xy_zz_xx_xz, g_xy_zz_xx_yy, g_xy_zz_xx_yz, g_xy_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_z_xx_xx[i] = -2.0 * g_xy_0_xx_xx[i] * a_exp + 4.0 * g_xy_zz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_xx_xy[i] = -2.0 * g_xy_0_xx_xy[i] * a_exp + 4.0 * g_xy_zz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_xx_xz[i] = -2.0 * g_xy_0_xx_xz[i] * a_exp + 4.0 * g_xy_zz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_xx_yy[i] = -2.0 * g_xy_0_xx_yy[i] * a_exp + 4.0 * g_xy_zz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_xx_yz[i] = -2.0 * g_xy_0_xx_yz[i] * a_exp + 4.0 * g_xy_zz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_xx_zz[i] = -2.0 * g_xy_0_xx_zz[i] * a_exp + 4.0 * g_xy_zz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_x_z_0_0_y_z_xy_xx, g_x_z_0_0_y_z_xy_xy, g_x_z_0_0_y_z_xy_xz, g_x_z_0_0_y_z_xy_yy, g_x_z_0_0_y_z_xy_yz, g_x_z_0_0_y_z_xy_zz, g_xy_0_xy_xx, g_xy_0_xy_xy, g_xy_0_xy_xz, g_xy_0_xy_yy, g_xy_0_xy_yz, g_xy_0_xy_zz, g_xy_zz_xy_xx, g_xy_zz_xy_xy, g_xy_zz_xy_xz, g_xy_zz_xy_yy, g_xy_zz_xy_yz, g_xy_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_z_xy_xx[i] = -2.0 * g_xy_0_xy_xx[i] * a_exp + 4.0 * g_xy_zz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_xy_xy[i] = -2.0 * g_xy_0_xy_xy[i] * a_exp + 4.0 * g_xy_zz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_xy_xz[i] = -2.0 * g_xy_0_xy_xz[i] * a_exp + 4.0 * g_xy_zz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_xy_yy[i] = -2.0 * g_xy_0_xy_yy[i] * a_exp + 4.0 * g_xy_zz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_xy_yz[i] = -2.0 * g_xy_0_xy_yz[i] * a_exp + 4.0 * g_xy_zz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_xy_zz[i] = -2.0 * g_xy_0_xy_zz[i] * a_exp + 4.0 * g_xy_zz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_x_z_0_0_y_z_xz_xx, g_x_z_0_0_y_z_xz_xy, g_x_z_0_0_y_z_xz_xz, g_x_z_0_0_y_z_xz_yy, g_x_z_0_0_y_z_xz_yz, g_x_z_0_0_y_z_xz_zz, g_xy_0_xz_xx, g_xy_0_xz_xy, g_xy_0_xz_xz, g_xy_0_xz_yy, g_xy_0_xz_yz, g_xy_0_xz_zz, g_xy_zz_xz_xx, g_xy_zz_xz_xy, g_xy_zz_xz_xz, g_xy_zz_xz_yy, g_xy_zz_xz_yz, g_xy_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_z_xz_xx[i] = -2.0 * g_xy_0_xz_xx[i] * a_exp + 4.0 * g_xy_zz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_xz_xy[i] = -2.0 * g_xy_0_xz_xy[i] * a_exp + 4.0 * g_xy_zz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_xz_xz[i] = -2.0 * g_xy_0_xz_xz[i] * a_exp + 4.0 * g_xy_zz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_xz_yy[i] = -2.0 * g_xy_0_xz_yy[i] * a_exp + 4.0 * g_xy_zz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_xz_yz[i] = -2.0 * g_xy_0_xz_yz[i] * a_exp + 4.0 * g_xy_zz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_xz_zz[i] = -2.0 * g_xy_0_xz_zz[i] * a_exp + 4.0 * g_xy_zz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_x_z_0_0_y_z_yy_xx, g_x_z_0_0_y_z_yy_xy, g_x_z_0_0_y_z_yy_xz, g_x_z_0_0_y_z_yy_yy, g_x_z_0_0_y_z_yy_yz, g_x_z_0_0_y_z_yy_zz, g_xy_0_yy_xx, g_xy_0_yy_xy, g_xy_0_yy_xz, g_xy_0_yy_yy, g_xy_0_yy_yz, g_xy_0_yy_zz, g_xy_zz_yy_xx, g_xy_zz_yy_xy, g_xy_zz_yy_xz, g_xy_zz_yy_yy, g_xy_zz_yy_yz, g_xy_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_z_yy_xx[i] = -2.0 * g_xy_0_yy_xx[i] * a_exp + 4.0 * g_xy_zz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_yy_xy[i] = -2.0 * g_xy_0_yy_xy[i] * a_exp + 4.0 * g_xy_zz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_yy_xz[i] = -2.0 * g_xy_0_yy_xz[i] * a_exp + 4.0 * g_xy_zz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_yy_yy[i] = -2.0 * g_xy_0_yy_yy[i] * a_exp + 4.0 * g_xy_zz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_yy_yz[i] = -2.0 * g_xy_0_yy_yz[i] * a_exp + 4.0 * g_xy_zz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_yy_zz[i] = -2.0 * g_xy_0_yy_zz[i] * a_exp + 4.0 * g_xy_zz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_x_z_0_0_y_z_yz_xx, g_x_z_0_0_y_z_yz_xy, g_x_z_0_0_y_z_yz_xz, g_x_z_0_0_y_z_yz_yy, g_x_z_0_0_y_z_yz_yz, g_x_z_0_0_y_z_yz_zz, g_xy_0_yz_xx, g_xy_0_yz_xy, g_xy_0_yz_xz, g_xy_0_yz_yy, g_xy_0_yz_yz, g_xy_0_yz_zz, g_xy_zz_yz_xx, g_xy_zz_yz_xy, g_xy_zz_yz_xz, g_xy_zz_yz_yy, g_xy_zz_yz_yz, g_xy_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_z_yz_xx[i] = -2.0 * g_xy_0_yz_xx[i] * a_exp + 4.0 * g_xy_zz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_yz_xy[i] = -2.0 * g_xy_0_yz_xy[i] * a_exp + 4.0 * g_xy_zz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_yz_xz[i] = -2.0 * g_xy_0_yz_xz[i] * a_exp + 4.0 * g_xy_zz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_yz_yy[i] = -2.0 * g_xy_0_yz_yy[i] * a_exp + 4.0 * g_xy_zz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_yz_yz[i] = -2.0 * g_xy_0_yz_yz[i] * a_exp + 4.0 * g_xy_zz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_yz_zz[i] = -2.0 * g_xy_0_yz_zz[i] * a_exp + 4.0 * g_xy_zz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_x_z_0_0_y_z_zz_xx, g_x_z_0_0_y_z_zz_xy, g_x_z_0_0_y_z_zz_xz, g_x_z_0_0_y_z_zz_yy, g_x_z_0_0_y_z_zz_yz, g_x_z_0_0_y_z_zz_zz, g_xy_0_zz_xx, g_xy_0_zz_xy, g_xy_0_zz_xz, g_xy_0_zz_yy, g_xy_0_zz_yz, g_xy_0_zz_zz, g_xy_zz_zz_xx, g_xy_zz_zz_xy, g_xy_zz_zz_xz, g_xy_zz_zz_yy, g_xy_zz_zz_yz, g_xy_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_z_zz_xx[i] = -2.0 * g_xy_0_zz_xx[i] * a_exp + 4.0 * g_xy_zz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_zz_xy[i] = -2.0 * g_xy_0_zz_xy[i] * a_exp + 4.0 * g_xy_zz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_zz_xz[i] = -2.0 * g_xy_0_zz_xz[i] * a_exp + 4.0 * g_xy_zz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_zz_yy[i] = -2.0 * g_xy_0_zz_yy[i] * a_exp + 4.0 * g_xy_zz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_zz_yz[i] = -2.0 * g_xy_0_zz_yz[i] * a_exp + 4.0 * g_xy_zz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_zz_zz[i] = -2.0 * g_xy_0_zz_zz[i] * a_exp + 4.0 * g_xy_zz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_x_z_0_0_z_x_xx_xx, g_x_z_0_0_z_x_xx_xy, g_x_z_0_0_z_x_xx_xz, g_x_z_0_0_z_x_xx_yy, g_x_z_0_0_z_x_xx_yz, g_x_z_0_0_z_x_xx_zz, g_xz_xz_xx_xx, g_xz_xz_xx_xy, g_xz_xz_xx_xz, g_xz_xz_xx_yy, g_xz_xz_xx_yz, g_xz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_x_xx_xx[i] = 4.0 * g_xz_xz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_xx_xy[i] = 4.0 * g_xz_xz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_xx_xz[i] = 4.0 * g_xz_xz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_xx_yy[i] = 4.0 * g_xz_xz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_xx_yz[i] = 4.0 * g_xz_xz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_xx_zz[i] = 4.0 * g_xz_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_x_z_0_0_z_x_xy_xx, g_x_z_0_0_z_x_xy_xy, g_x_z_0_0_z_x_xy_xz, g_x_z_0_0_z_x_xy_yy, g_x_z_0_0_z_x_xy_yz, g_x_z_0_0_z_x_xy_zz, g_xz_xz_xy_xx, g_xz_xz_xy_xy, g_xz_xz_xy_xz, g_xz_xz_xy_yy, g_xz_xz_xy_yz, g_xz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_x_xy_xx[i] = 4.0 * g_xz_xz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_xy_xy[i] = 4.0 * g_xz_xz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_xy_xz[i] = 4.0 * g_xz_xz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_xy_yy[i] = 4.0 * g_xz_xz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_xy_yz[i] = 4.0 * g_xz_xz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_xy_zz[i] = 4.0 * g_xz_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_x_z_0_0_z_x_xz_xx, g_x_z_0_0_z_x_xz_xy, g_x_z_0_0_z_x_xz_xz, g_x_z_0_0_z_x_xz_yy, g_x_z_0_0_z_x_xz_yz, g_x_z_0_0_z_x_xz_zz, g_xz_xz_xz_xx, g_xz_xz_xz_xy, g_xz_xz_xz_xz, g_xz_xz_xz_yy, g_xz_xz_xz_yz, g_xz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_x_xz_xx[i] = 4.0 * g_xz_xz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_xz_xy[i] = 4.0 * g_xz_xz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_xz_xz[i] = 4.0 * g_xz_xz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_xz_yy[i] = 4.0 * g_xz_xz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_xz_yz[i] = 4.0 * g_xz_xz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_xz_zz[i] = 4.0 * g_xz_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_x_z_0_0_z_x_yy_xx, g_x_z_0_0_z_x_yy_xy, g_x_z_0_0_z_x_yy_xz, g_x_z_0_0_z_x_yy_yy, g_x_z_0_0_z_x_yy_yz, g_x_z_0_0_z_x_yy_zz, g_xz_xz_yy_xx, g_xz_xz_yy_xy, g_xz_xz_yy_xz, g_xz_xz_yy_yy, g_xz_xz_yy_yz, g_xz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_x_yy_xx[i] = 4.0 * g_xz_xz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_yy_xy[i] = 4.0 * g_xz_xz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_yy_xz[i] = 4.0 * g_xz_xz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_yy_yy[i] = 4.0 * g_xz_xz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_yy_yz[i] = 4.0 * g_xz_xz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_yy_zz[i] = 4.0 * g_xz_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_x_z_0_0_z_x_yz_xx, g_x_z_0_0_z_x_yz_xy, g_x_z_0_0_z_x_yz_xz, g_x_z_0_0_z_x_yz_yy, g_x_z_0_0_z_x_yz_yz, g_x_z_0_0_z_x_yz_zz, g_xz_xz_yz_xx, g_xz_xz_yz_xy, g_xz_xz_yz_xz, g_xz_xz_yz_yy, g_xz_xz_yz_yz, g_xz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_x_yz_xx[i] = 4.0 * g_xz_xz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_yz_xy[i] = 4.0 * g_xz_xz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_yz_xz[i] = 4.0 * g_xz_xz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_yz_yy[i] = 4.0 * g_xz_xz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_yz_yz[i] = 4.0 * g_xz_xz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_yz_zz[i] = 4.0 * g_xz_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_x_z_0_0_z_x_zz_xx, g_x_z_0_0_z_x_zz_xy, g_x_z_0_0_z_x_zz_xz, g_x_z_0_0_z_x_zz_yy, g_x_z_0_0_z_x_zz_yz, g_x_z_0_0_z_x_zz_zz, g_xz_xz_zz_xx, g_xz_xz_zz_xy, g_xz_xz_zz_xz, g_xz_xz_zz_yy, g_xz_xz_zz_yz, g_xz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_x_zz_xx[i] = 4.0 * g_xz_xz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_zz_xy[i] = 4.0 * g_xz_xz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_zz_xz[i] = 4.0 * g_xz_xz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_zz_yy[i] = 4.0 * g_xz_xz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_zz_yz[i] = 4.0 * g_xz_xz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_zz_zz[i] = 4.0 * g_xz_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_x_z_0_0_z_y_xx_xx, g_x_z_0_0_z_y_xx_xy, g_x_z_0_0_z_y_xx_xz, g_x_z_0_0_z_y_xx_yy, g_x_z_0_0_z_y_xx_yz, g_x_z_0_0_z_y_xx_zz, g_xz_yz_xx_xx, g_xz_yz_xx_xy, g_xz_yz_xx_xz, g_xz_yz_xx_yy, g_xz_yz_xx_yz, g_xz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_y_xx_xx[i] = 4.0 * g_xz_yz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_xx_xy[i] = 4.0 * g_xz_yz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_xx_xz[i] = 4.0 * g_xz_yz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_xx_yy[i] = 4.0 * g_xz_yz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_xx_yz[i] = 4.0 * g_xz_yz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_xx_zz[i] = 4.0 * g_xz_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_x_z_0_0_z_y_xy_xx, g_x_z_0_0_z_y_xy_xy, g_x_z_0_0_z_y_xy_xz, g_x_z_0_0_z_y_xy_yy, g_x_z_0_0_z_y_xy_yz, g_x_z_0_0_z_y_xy_zz, g_xz_yz_xy_xx, g_xz_yz_xy_xy, g_xz_yz_xy_xz, g_xz_yz_xy_yy, g_xz_yz_xy_yz, g_xz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_y_xy_xx[i] = 4.0 * g_xz_yz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_xy_xy[i] = 4.0 * g_xz_yz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_xy_xz[i] = 4.0 * g_xz_yz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_xy_yy[i] = 4.0 * g_xz_yz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_xy_yz[i] = 4.0 * g_xz_yz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_xy_zz[i] = 4.0 * g_xz_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_x_z_0_0_z_y_xz_xx, g_x_z_0_0_z_y_xz_xy, g_x_z_0_0_z_y_xz_xz, g_x_z_0_0_z_y_xz_yy, g_x_z_0_0_z_y_xz_yz, g_x_z_0_0_z_y_xz_zz, g_xz_yz_xz_xx, g_xz_yz_xz_xy, g_xz_yz_xz_xz, g_xz_yz_xz_yy, g_xz_yz_xz_yz, g_xz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_y_xz_xx[i] = 4.0 * g_xz_yz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_xz_xy[i] = 4.0 * g_xz_yz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_xz_xz[i] = 4.0 * g_xz_yz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_xz_yy[i] = 4.0 * g_xz_yz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_xz_yz[i] = 4.0 * g_xz_yz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_xz_zz[i] = 4.0 * g_xz_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_x_z_0_0_z_y_yy_xx, g_x_z_0_0_z_y_yy_xy, g_x_z_0_0_z_y_yy_xz, g_x_z_0_0_z_y_yy_yy, g_x_z_0_0_z_y_yy_yz, g_x_z_0_0_z_y_yy_zz, g_xz_yz_yy_xx, g_xz_yz_yy_xy, g_xz_yz_yy_xz, g_xz_yz_yy_yy, g_xz_yz_yy_yz, g_xz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_y_yy_xx[i] = 4.0 * g_xz_yz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_yy_xy[i] = 4.0 * g_xz_yz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_yy_xz[i] = 4.0 * g_xz_yz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_yy_yy[i] = 4.0 * g_xz_yz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_yy_yz[i] = 4.0 * g_xz_yz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_yy_zz[i] = 4.0 * g_xz_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_x_z_0_0_z_y_yz_xx, g_x_z_0_0_z_y_yz_xy, g_x_z_0_0_z_y_yz_xz, g_x_z_0_0_z_y_yz_yy, g_x_z_0_0_z_y_yz_yz, g_x_z_0_0_z_y_yz_zz, g_xz_yz_yz_xx, g_xz_yz_yz_xy, g_xz_yz_yz_xz, g_xz_yz_yz_yy, g_xz_yz_yz_yz, g_xz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_y_yz_xx[i] = 4.0 * g_xz_yz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_yz_xy[i] = 4.0 * g_xz_yz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_yz_xz[i] = 4.0 * g_xz_yz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_yz_yy[i] = 4.0 * g_xz_yz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_yz_yz[i] = 4.0 * g_xz_yz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_yz_zz[i] = 4.0 * g_xz_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_x_z_0_0_z_y_zz_xx, g_x_z_0_0_z_y_zz_xy, g_x_z_0_0_z_y_zz_xz, g_x_z_0_0_z_y_zz_yy, g_x_z_0_0_z_y_zz_yz, g_x_z_0_0_z_y_zz_zz, g_xz_yz_zz_xx, g_xz_yz_zz_xy, g_xz_yz_zz_xz, g_xz_yz_zz_yy, g_xz_yz_zz_yz, g_xz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_y_zz_xx[i] = 4.0 * g_xz_yz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_zz_xy[i] = 4.0 * g_xz_yz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_zz_xz[i] = 4.0 * g_xz_yz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_zz_yy[i] = 4.0 * g_xz_yz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_zz_yz[i] = 4.0 * g_xz_yz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_zz_zz[i] = 4.0 * g_xz_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_x_z_0_0_z_z_xx_xx, g_x_z_0_0_z_z_xx_xy, g_x_z_0_0_z_z_xx_xz, g_x_z_0_0_z_z_xx_yy, g_x_z_0_0_z_z_xx_yz, g_x_z_0_0_z_z_xx_zz, g_xz_0_xx_xx, g_xz_0_xx_xy, g_xz_0_xx_xz, g_xz_0_xx_yy, g_xz_0_xx_yz, g_xz_0_xx_zz, g_xz_zz_xx_xx, g_xz_zz_xx_xy, g_xz_zz_xx_xz, g_xz_zz_xx_yy, g_xz_zz_xx_yz, g_xz_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_z_xx_xx[i] = -2.0 * g_xz_0_xx_xx[i] * a_exp + 4.0 * g_xz_zz_xx_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_xx_xy[i] = -2.0 * g_xz_0_xx_xy[i] * a_exp + 4.0 * g_xz_zz_xx_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_xx_xz[i] = -2.0 * g_xz_0_xx_xz[i] * a_exp + 4.0 * g_xz_zz_xx_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_xx_yy[i] = -2.0 * g_xz_0_xx_yy[i] * a_exp + 4.0 * g_xz_zz_xx_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_xx_yz[i] = -2.0 * g_xz_0_xx_yz[i] * a_exp + 4.0 * g_xz_zz_xx_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_xx_zz[i] = -2.0 * g_xz_0_xx_zz[i] * a_exp + 4.0 * g_xz_zz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_x_z_0_0_z_z_xy_xx, g_x_z_0_0_z_z_xy_xy, g_x_z_0_0_z_z_xy_xz, g_x_z_0_0_z_z_xy_yy, g_x_z_0_0_z_z_xy_yz, g_x_z_0_0_z_z_xy_zz, g_xz_0_xy_xx, g_xz_0_xy_xy, g_xz_0_xy_xz, g_xz_0_xy_yy, g_xz_0_xy_yz, g_xz_0_xy_zz, g_xz_zz_xy_xx, g_xz_zz_xy_xy, g_xz_zz_xy_xz, g_xz_zz_xy_yy, g_xz_zz_xy_yz, g_xz_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_z_xy_xx[i] = -2.0 * g_xz_0_xy_xx[i] * a_exp + 4.0 * g_xz_zz_xy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_xy_xy[i] = -2.0 * g_xz_0_xy_xy[i] * a_exp + 4.0 * g_xz_zz_xy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_xy_xz[i] = -2.0 * g_xz_0_xy_xz[i] * a_exp + 4.0 * g_xz_zz_xy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_xy_yy[i] = -2.0 * g_xz_0_xy_yy[i] * a_exp + 4.0 * g_xz_zz_xy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_xy_yz[i] = -2.0 * g_xz_0_xy_yz[i] * a_exp + 4.0 * g_xz_zz_xy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_xy_zz[i] = -2.0 * g_xz_0_xy_zz[i] * a_exp + 4.0 * g_xz_zz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_x_z_0_0_z_z_xz_xx, g_x_z_0_0_z_z_xz_xy, g_x_z_0_0_z_z_xz_xz, g_x_z_0_0_z_z_xz_yy, g_x_z_0_0_z_z_xz_yz, g_x_z_0_0_z_z_xz_zz, g_xz_0_xz_xx, g_xz_0_xz_xy, g_xz_0_xz_xz, g_xz_0_xz_yy, g_xz_0_xz_yz, g_xz_0_xz_zz, g_xz_zz_xz_xx, g_xz_zz_xz_xy, g_xz_zz_xz_xz, g_xz_zz_xz_yy, g_xz_zz_xz_yz, g_xz_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_z_xz_xx[i] = -2.0 * g_xz_0_xz_xx[i] * a_exp + 4.0 * g_xz_zz_xz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_xz_xy[i] = -2.0 * g_xz_0_xz_xy[i] * a_exp + 4.0 * g_xz_zz_xz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_xz_xz[i] = -2.0 * g_xz_0_xz_xz[i] * a_exp + 4.0 * g_xz_zz_xz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_xz_yy[i] = -2.0 * g_xz_0_xz_yy[i] * a_exp + 4.0 * g_xz_zz_xz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_xz_yz[i] = -2.0 * g_xz_0_xz_yz[i] * a_exp + 4.0 * g_xz_zz_xz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_xz_zz[i] = -2.0 * g_xz_0_xz_zz[i] * a_exp + 4.0 * g_xz_zz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_x_z_0_0_z_z_yy_xx, g_x_z_0_0_z_z_yy_xy, g_x_z_0_0_z_z_yy_xz, g_x_z_0_0_z_z_yy_yy, g_x_z_0_0_z_z_yy_yz, g_x_z_0_0_z_z_yy_zz, g_xz_0_yy_xx, g_xz_0_yy_xy, g_xz_0_yy_xz, g_xz_0_yy_yy, g_xz_0_yy_yz, g_xz_0_yy_zz, g_xz_zz_yy_xx, g_xz_zz_yy_xy, g_xz_zz_yy_xz, g_xz_zz_yy_yy, g_xz_zz_yy_yz, g_xz_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_z_yy_xx[i] = -2.0 * g_xz_0_yy_xx[i] * a_exp + 4.0 * g_xz_zz_yy_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_yy_xy[i] = -2.0 * g_xz_0_yy_xy[i] * a_exp + 4.0 * g_xz_zz_yy_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_yy_xz[i] = -2.0 * g_xz_0_yy_xz[i] * a_exp + 4.0 * g_xz_zz_yy_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_yy_yy[i] = -2.0 * g_xz_0_yy_yy[i] * a_exp + 4.0 * g_xz_zz_yy_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_yy_yz[i] = -2.0 * g_xz_0_yy_yz[i] * a_exp + 4.0 * g_xz_zz_yy_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_yy_zz[i] = -2.0 * g_xz_0_yy_zz[i] * a_exp + 4.0 * g_xz_zz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_x_z_0_0_z_z_yz_xx, g_x_z_0_0_z_z_yz_xy, g_x_z_0_0_z_z_yz_xz, g_x_z_0_0_z_z_yz_yy, g_x_z_0_0_z_z_yz_yz, g_x_z_0_0_z_z_yz_zz, g_xz_0_yz_xx, g_xz_0_yz_xy, g_xz_0_yz_xz, g_xz_0_yz_yy, g_xz_0_yz_yz, g_xz_0_yz_zz, g_xz_zz_yz_xx, g_xz_zz_yz_xy, g_xz_zz_yz_xz, g_xz_zz_yz_yy, g_xz_zz_yz_yz, g_xz_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_z_yz_xx[i] = -2.0 * g_xz_0_yz_xx[i] * a_exp + 4.0 * g_xz_zz_yz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_yz_xy[i] = -2.0 * g_xz_0_yz_xy[i] * a_exp + 4.0 * g_xz_zz_yz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_yz_xz[i] = -2.0 * g_xz_0_yz_xz[i] * a_exp + 4.0 * g_xz_zz_yz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_yz_yy[i] = -2.0 * g_xz_0_yz_yy[i] * a_exp + 4.0 * g_xz_zz_yz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_yz_yz[i] = -2.0 * g_xz_0_yz_yz[i] * a_exp + 4.0 * g_xz_zz_yz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_yz_zz[i] = -2.0 * g_xz_0_yz_zz[i] * a_exp + 4.0 * g_xz_zz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_x_z_0_0_z_z_zz_xx, g_x_z_0_0_z_z_zz_xy, g_x_z_0_0_z_z_zz_xz, g_x_z_0_0_z_z_zz_yy, g_x_z_0_0_z_z_zz_yz, g_x_z_0_0_z_z_zz_zz, g_xz_0_zz_xx, g_xz_0_zz_xy, g_xz_0_zz_xz, g_xz_0_zz_yy, g_xz_0_zz_yz, g_xz_0_zz_zz, g_xz_zz_zz_xx, g_xz_zz_zz_xy, g_xz_zz_zz_xz, g_xz_zz_zz_yy, g_xz_zz_zz_yz, g_xz_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_z_zz_xx[i] = -2.0 * g_xz_0_zz_xx[i] * a_exp + 4.0 * g_xz_zz_zz_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_zz_xy[i] = -2.0 * g_xz_0_zz_xy[i] * a_exp + 4.0 * g_xz_zz_zz_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_zz_xz[i] = -2.0 * g_xz_0_zz_xz[i] * a_exp + 4.0 * g_xz_zz_zz_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_zz_yy[i] = -2.0 * g_xz_0_zz_yy[i] * a_exp + 4.0 * g_xz_zz_zz_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_zz_yz[i] = -2.0 * g_xz_0_zz_yz[i] * a_exp + 4.0 * g_xz_zz_zz_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_zz_zz[i] = -2.0 * g_xz_0_zz_zz[i] * a_exp + 4.0 * g_xz_zz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (972-978)

    #pragma omp simd aligned(g_xy_0_xx_xx, g_xy_0_xx_xy, g_xy_0_xx_xz, g_xy_0_xx_yy, g_xy_0_xx_yz, g_xy_0_xx_zz, g_xy_xx_xx_xx, g_xy_xx_xx_xy, g_xy_xx_xx_xz, g_xy_xx_xx_yy, g_xy_xx_xx_yz, g_xy_xx_xx_zz, g_y_x_0_0_x_x_xx_xx, g_y_x_0_0_x_x_xx_xy, g_y_x_0_0_x_x_xx_xz, g_y_x_0_0_x_x_xx_yy, g_y_x_0_0_x_x_xx_yz, g_y_x_0_0_x_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_x_xx_xx[i] = -2.0 * g_xy_0_xx_xx[i] * a_exp + 4.0 * g_xy_xx_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_xx_xy[i] = -2.0 * g_xy_0_xx_xy[i] * a_exp + 4.0 * g_xy_xx_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_xx_xz[i] = -2.0 * g_xy_0_xx_xz[i] * a_exp + 4.0 * g_xy_xx_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_xx_yy[i] = -2.0 * g_xy_0_xx_yy[i] * a_exp + 4.0 * g_xy_xx_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_xx_yz[i] = -2.0 * g_xy_0_xx_yz[i] * a_exp + 4.0 * g_xy_xx_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_xx_zz[i] = -2.0 * g_xy_0_xx_zz[i] * a_exp + 4.0 * g_xy_xx_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (978-984)

    #pragma omp simd aligned(g_xy_0_xy_xx, g_xy_0_xy_xy, g_xy_0_xy_xz, g_xy_0_xy_yy, g_xy_0_xy_yz, g_xy_0_xy_zz, g_xy_xx_xy_xx, g_xy_xx_xy_xy, g_xy_xx_xy_xz, g_xy_xx_xy_yy, g_xy_xx_xy_yz, g_xy_xx_xy_zz, g_y_x_0_0_x_x_xy_xx, g_y_x_0_0_x_x_xy_xy, g_y_x_0_0_x_x_xy_xz, g_y_x_0_0_x_x_xy_yy, g_y_x_0_0_x_x_xy_yz, g_y_x_0_0_x_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_x_xy_xx[i] = -2.0 * g_xy_0_xy_xx[i] * a_exp + 4.0 * g_xy_xx_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_xy_xy[i] = -2.0 * g_xy_0_xy_xy[i] * a_exp + 4.0 * g_xy_xx_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_xy_xz[i] = -2.0 * g_xy_0_xy_xz[i] * a_exp + 4.0 * g_xy_xx_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_xy_yy[i] = -2.0 * g_xy_0_xy_yy[i] * a_exp + 4.0 * g_xy_xx_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_xy_yz[i] = -2.0 * g_xy_0_xy_yz[i] * a_exp + 4.0 * g_xy_xx_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_xy_zz[i] = -2.0 * g_xy_0_xy_zz[i] * a_exp + 4.0 * g_xy_xx_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (984-990)

    #pragma omp simd aligned(g_xy_0_xz_xx, g_xy_0_xz_xy, g_xy_0_xz_xz, g_xy_0_xz_yy, g_xy_0_xz_yz, g_xy_0_xz_zz, g_xy_xx_xz_xx, g_xy_xx_xz_xy, g_xy_xx_xz_xz, g_xy_xx_xz_yy, g_xy_xx_xz_yz, g_xy_xx_xz_zz, g_y_x_0_0_x_x_xz_xx, g_y_x_0_0_x_x_xz_xy, g_y_x_0_0_x_x_xz_xz, g_y_x_0_0_x_x_xz_yy, g_y_x_0_0_x_x_xz_yz, g_y_x_0_0_x_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_x_xz_xx[i] = -2.0 * g_xy_0_xz_xx[i] * a_exp + 4.0 * g_xy_xx_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_xz_xy[i] = -2.0 * g_xy_0_xz_xy[i] * a_exp + 4.0 * g_xy_xx_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_xz_xz[i] = -2.0 * g_xy_0_xz_xz[i] * a_exp + 4.0 * g_xy_xx_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_xz_yy[i] = -2.0 * g_xy_0_xz_yy[i] * a_exp + 4.0 * g_xy_xx_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_xz_yz[i] = -2.0 * g_xy_0_xz_yz[i] * a_exp + 4.0 * g_xy_xx_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_xz_zz[i] = -2.0 * g_xy_0_xz_zz[i] * a_exp + 4.0 * g_xy_xx_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (990-996)

    #pragma omp simd aligned(g_xy_0_yy_xx, g_xy_0_yy_xy, g_xy_0_yy_xz, g_xy_0_yy_yy, g_xy_0_yy_yz, g_xy_0_yy_zz, g_xy_xx_yy_xx, g_xy_xx_yy_xy, g_xy_xx_yy_xz, g_xy_xx_yy_yy, g_xy_xx_yy_yz, g_xy_xx_yy_zz, g_y_x_0_0_x_x_yy_xx, g_y_x_0_0_x_x_yy_xy, g_y_x_0_0_x_x_yy_xz, g_y_x_0_0_x_x_yy_yy, g_y_x_0_0_x_x_yy_yz, g_y_x_0_0_x_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_x_yy_xx[i] = -2.0 * g_xy_0_yy_xx[i] * a_exp + 4.0 * g_xy_xx_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_yy_xy[i] = -2.0 * g_xy_0_yy_xy[i] * a_exp + 4.0 * g_xy_xx_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_yy_xz[i] = -2.0 * g_xy_0_yy_xz[i] * a_exp + 4.0 * g_xy_xx_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_yy_yy[i] = -2.0 * g_xy_0_yy_yy[i] * a_exp + 4.0 * g_xy_xx_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_yy_yz[i] = -2.0 * g_xy_0_yy_yz[i] * a_exp + 4.0 * g_xy_xx_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_yy_zz[i] = -2.0 * g_xy_0_yy_zz[i] * a_exp + 4.0 * g_xy_xx_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (996-1002)

    #pragma omp simd aligned(g_xy_0_yz_xx, g_xy_0_yz_xy, g_xy_0_yz_xz, g_xy_0_yz_yy, g_xy_0_yz_yz, g_xy_0_yz_zz, g_xy_xx_yz_xx, g_xy_xx_yz_xy, g_xy_xx_yz_xz, g_xy_xx_yz_yy, g_xy_xx_yz_yz, g_xy_xx_yz_zz, g_y_x_0_0_x_x_yz_xx, g_y_x_0_0_x_x_yz_xy, g_y_x_0_0_x_x_yz_xz, g_y_x_0_0_x_x_yz_yy, g_y_x_0_0_x_x_yz_yz, g_y_x_0_0_x_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_x_yz_xx[i] = -2.0 * g_xy_0_yz_xx[i] * a_exp + 4.0 * g_xy_xx_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_yz_xy[i] = -2.0 * g_xy_0_yz_xy[i] * a_exp + 4.0 * g_xy_xx_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_yz_xz[i] = -2.0 * g_xy_0_yz_xz[i] * a_exp + 4.0 * g_xy_xx_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_yz_yy[i] = -2.0 * g_xy_0_yz_yy[i] * a_exp + 4.0 * g_xy_xx_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_yz_yz[i] = -2.0 * g_xy_0_yz_yz[i] * a_exp + 4.0 * g_xy_xx_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_yz_zz[i] = -2.0 * g_xy_0_yz_zz[i] * a_exp + 4.0 * g_xy_xx_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1002-1008)

    #pragma omp simd aligned(g_xy_0_zz_xx, g_xy_0_zz_xy, g_xy_0_zz_xz, g_xy_0_zz_yy, g_xy_0_zz_yz, g_xy_0_zz_zz, g_xy_xx_zz_xx, g_xy_xx_zz_xy, g_xy_xx_zz_xz, g_xy_xx_zz_yy, g_xy_xx_zz_yz, g_xy_xx_zz_zz, g_y_x_0_0_x_x_zz_xx, g_y_x_0_0_x_x_zz_xy, g_y_x_0_0_x_x_zz_xz, g_y_x_0_0_x_x_zz_yy, g_y_x_0_0_x_x_zz_yz, g_y_x_0_0_x_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_x_zz_xx[i] = -2.0 * g_xy_0_zz_xx[i] * a_exp + 4.0 * g_xy_xx_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_zz_xy[i] = -2.0 * g_xy_0_zz_xy[i] * a_exp + 4.0 * g_xy_xx_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_zz_xz[i] = -2.0 * g_xy_0_zz_xz[i] * a_exp + 4.0 * g_xy_xx_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_zz_yy[i] = -2.0 * g_xy_0_zz_yy[i] * a_exp + 4.0 * g_xy_xx_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_zz_yz[i] = -2.0 * g_xy_0_zz_yz[i] * a_exp + 4.0 * g_xy_xx_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_zz_zz[i] = -2.0 * g_xy_0_zz_zz[i] * a_exp + 4.0 * g_xy_xx_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1008-1014)

    #pragma omp simd aligned(g_xy_xy_xx_xx, g_xy_xy_xx_xy, g_xy_xy_xx_xz, g_xy_xy_xx_yy, g_xy_xy_xx_yz, g_xy_xy_xx_zz, g_y_x_0_0_x_y_xx_xx, g_y_x_0_0_x_y_xx_xy, g_y_x_0_0_x_y_xx_xz, g_y_x_0_0_x_y_xx_yy, g_y_x_0_0_x_y_xx_yz, g_y_x_0_0_x_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_y_xx_xx[i] = 4.0 * g_xy_xy_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_xx_xy[i] = 4.0 * g_xy_xy_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_xx_xz[i] = 4.0 * g_xy_xy_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_xx_yy[i] = 4.0 * g_xy_xy_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_xx_yz[i] = 4.0 * g_xy_xy_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_xx_zz[i] = 4.0 * g_xy_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1014-1020)

    #pragma omp simd aligned(g_xy_xy_xy_xx, g_xy_xy_xy_xy, g_xy_xy_xy_xz, g_xy_xy_xy_yy, g_xy_xy_xy_yz, g_xy_xy_xy_zz, g_y_x_0_0_x_y_xy_xx, g_y_x_0_0_x_y_xy_xy, g_y_x_0_0_x_y_xy_xz, g_y_x_0_0_x_y_xy_yy, g_y_x_0_0_x_y_xy_yz, g_y_x_0_0_x_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_y_xy_xx[i] = 4.0 * g_xy_xy_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_xy_xy[i] = 4.0 * g_xy_xy_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_xy_xz[i] = 4.0 * g_xy_xy_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_xy_yy[i] = 4.0 * g_xy_xy_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_xy_yz[i] = 4.0 * g_xy_xy_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_xy_zz[i] = 4.0 * g_xy_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1020-1026)

    #pragma omp simd aligned(g_xy_xy_xz_xx, g_xy_xy_xz_xy, g_xy_xy_xz_xz, g_xy_xy_xz_yy, g_xy_xy_xz_yz, g_xy_xy_xz_zz, g_y_x_0_0_x_y_xz_xx, g_y_x_0_0_x_y_xz_xy, g_y_x_0_0_x_y_xz_xz, g_y_x_0_0_x_y_xz_yy, g_y_x_0_0_x_y_xz_yz, g_y_x_0_0_x_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_y_xz_xx[i] = 4.0 * g_xy_xy_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_xz_xy[i] = 4.0 * g_xy_xy_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_xz_xz[i] = 4.0 * g_xy_xy_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_xz_yy[i] = 4.0 * g_xy_xy_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_xz_yz[i] = 4.0 * g_xy_xy_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_xz_zz[i] = 4.0 * g_xy_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1026-1032)

    #pragma omp simd aligned(g_xy_xy_yy_xx, g_xy_xy_yy_xy, g_xy_xy_yy_xz, g_xy_xy_yy_yy, g_xy_xy_yy_yz, g_xy_xy_yy_zz, g_y_x_0_0_x_y_yy_xx, g_y_x_0_0_x_y_yy_xy, g_y_x_0_0_x_y_yy_xz, g_y_x_0_0_x_y_yy_yy, g_y_x_0_0_x_y_yy_yz, g_y_x_0_0_x_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_y_yy_xx[i] = 4.0 * g_xy_xy_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_yy_xy[i] = 4.0 * g_xy_xy_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_yy_xz[i] = 4.0 * g_xy_xy_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_yy_yy[i] = 4.0 * g_xy_xy_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_yy_yz[i] = 4.0 * g_xy_xy_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_yy_zz[i] = 4.0 * g_xy_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1032-1038)

    #pragma omp simd aligned(g_xy_xy_yz_xx, g_xy_xy_yz_xy, g_xy_xy_yz_xz, g_xy_xy_yz_yy, g_xy_xy_yz_yz, g_xy_xy_yz_zz, g_y_x_0_0_x_y_yz_xx, g_y_x_0_0_x_y_yz_xy, g_y_x_0_0_x_y_yz_xz, g_y_x_0_0_x_y_yz_yy, g_y_x_0_0_x_y_yz_yz, g_y_x_0_0_x_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_y_yz_xx[i] = 4.0 * g_xy_xy_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_yz_xy[i] = 4.0 * g_xy_xy_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_yz_xz[i] = 4.0 * g_xy_xy_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_yz_yy[i] = 4.0 * g_xy_xy_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_yz_yz[i] = 4.0 * g_xy_xy_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_yz_zz[i] = 4.0 * g_xy_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1038-1044)

    #pragma omp simd aligned(g_xy_xy_zz_xx, g_xy_xy_zz_xy, g_xy_xy_zz_xz, g_xy_xy_zz_yy, g_xy_xy_zz_yz, g_xy_xy_zz_zz, g_y_x_0_0_x_y_zz_xx, g_y_x_0_0_x_y_zz_xy, g_y_x_0_0_x_y_zz_xz, g_y_x_0_0_x_y_zz_yy, g_y_x_0_0_x_y_zz_yz, g_y_x_0_0_x_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_y_zz_xx[i] = 4.0 * g_xy_xy_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_zz_xy[i] = 4.0 * g_xy_xy_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_zz_xz[i] = 4.0 * g_xy_xy_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_zz_yy[i] = 4.0 * g_xy_xy_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_zz_yz[i] = 4.0 * g_xy_xy_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_zz_zz[i] = 4.0 * g_xy_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1044-1050)

    #pragma omp simd aligned(g_xy_xz_xx_xx, g_xy_xz_xx_xy, g_xy_xz_xx_xz, g_xy_xz_xx_yy, g_xy_xz_xx_yz, g_xy_xz_xx_zz, g_y_x_0_0_x_z_xx_xx, g_y_x_0_0_x_z_xx_xy, g_y_x_0_0_x_z_xx_xz, g_y_x_0_0_x_z_xx_yy, g_y_x_0_0_x_z_xx_yz, g_y_x_0_0_x_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_z_xx_xx[i] = 4.0 * g_xy_xz_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_xx_xy[i] = 4.0 * g_xy_xz_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_xx_xz[i] = 4.0 * g_xy_xz_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_xx_yy[i] = 4.0 * g_xy_xz_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_xx_yz[i] = 4.0 * g_xy_xz_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_xx_zz[i] = 4.0 * g_xy_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1050-1056)

    #pragma omp simd aligned(g_xy_xz_xy_xx, g_xy_xz_xy_xy, g_xy_xz_xy_xz, g_xy_xz_xy_yy, g_xy_xz_xy_yz, g_xy_xz_xy_zz, g_y_x_0_0_x_z_xy_xx, g_y_x_0_0_x_z_xy_xy, g_y_x_0_0_x_z_xy_xz, g_y_x_0_0_x_z_xy_yy, g_y_x_0_0_x_z_xy_yz, g_y_x_0_0_x_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_z_xy_xx[i] = 4.0 * g_xy_xz_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_xy_xy[i] = 4.0 * g_xy_xz_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_xy_xz[i] = 4.0 * g_xy_xz_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_xy_yy[i] = 4.0 * g_xy_xz_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_xy_yz[i] = 4.0 * g_xy_xz_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_xy_zz[i] = 4.0 * g_xy_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1056-1062)

    #pragma omp simd aligned(g_xy_xz_xz_xx, g_xy_xz_xz_xy, g_xy_xz_xz_xz, g_xy_xz_xz_yy, g_xy_xz_xz_yz, g_xy_xz_xz_zz, g_y_x_0_0_x_z_xz_xx, g_y_x_0_0_x_z_xz_xy, g_y_x_0_0_x_z_xz_xz, g_y_x_0_0_x_z_xz_yy, g_y_x_0_0_x_z_xz_yz, g_y_x_0_0_x_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_z_xz_xx[i] = 4.0 * g_xy_xz_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_xz_xy[i] = 4.0 * g_xy_xz_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_xz_xz[i] = 4.0 * g_xy_xz_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_xz_yy[i] = 4.0 * g_xy_xz_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_xz_yz[i] = 4.0 * g_xy_xz_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_xz_zz[i] = 4.0 * g_xy_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1062-1068)

    #pragma omp simd aligned(g_xy_xz_yy_xx, g_xy_xz_yy_xy, g_xy_xz_yy_xz, g_xy_xz_yy_yy, g_xy_xz_yy_yz, g_xy_xz_yy_zz, g_y_x_0_0_x_z_yy_xx, g_y_x_0_0_x_z_yy_xy, g_y_x_0_0_x_z_yy_xz, g_y_x_0_0_x_z_yy_yy, g_y_x_0_0_x_z_yy_yz, g_y_x_0_0_x_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_z_yy_xx[i] = 4.0 * g_xy_xz_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_yy_xy[i] = 4.0 * g_xy_xz_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_yy_xz[i] = 4.0 * g_xy_xz_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_yy_yy[i] = 4.0 * g_xy_xz_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_yy_yz[i] = 4.0 * g_xy_xz_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_yy_zz[i] = 4.0 * g_xy_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1068-1074)

    #pragma omp simd aligned(g_xy_xz_yz_xx, g_xy_xz_yz_xy, g_xy_xz_yz_xz, g_xy_xz_yz_yy, g_xy_xz_yz_yz, g_xy_xz_yz_zz, g_y_x_0_0_x_z_yz_xx, g_y_x_0_0_x_z_yz_xy, g_y_x_0_0_x_z_yz_xz, g_y_x_0_0_x_z_yz_yy, g_y_x_0_0_x_z_yz_yz, g_y_x_0_0_x_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_z_yz_xx[i] = 4.0 * g_xy_xz_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_yz_xy[i] = 4.0 * g_xy_xz_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_yz_xz[i] = 4.0 * g_xy_xz_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_yz_yy[i] = 4.0 * g_xy_xz_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_yz_yz[i] = 4.0 * g_xy_xz_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_yz_zz[i] = 4.0 * g_xy_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1074-1080)

    #pragma omp simd aligned(g_xy_xz_zz_xx, g_xy_xz_zz_xy, g_xy_xz_zz_xz, g_xy_xz_zz_yy, g_xy_xz_zz_yz, g_xy_xz_zz_zz, g_y_x_0_0_x_z_zz_xx, g_y_x_0_0_x_z_zz_xy, g_y_x_0_0_x_z_zz_xz, g_y_x_0_0_x_z_zz_yy, g_y_x_0_0_x_z_zz_yz, g_y_x_0_0_x_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_z_zz_xx[i] = 4.0 * g_xy_xz_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_zz_xy[i] = 4.0 * g_xy_xz_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_zz_xz[i] = 4.0 * g_xy_xz_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_zz_yy[i] = 4.0 * g_xy_xz_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_zz_yz[i] = 4.0 * g_xy_xz_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_zz_zz[i] = 4.0 * g_xy_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1080-1086)

    #pragma omp simd aligned(g_0_0_xx_xx, g_0_0_xx_xy, g_0_0_xx_xz, g_0_0_xx_yy, g_0_0_xx_yz, g_0_0_xx_zz, g_0_xx_xx_xx, g_0_xx_xx_xy, g_0_xx_xx_xz, g_0_xx_xx_yy, g_0_xx_xx_yz, g_0_xx_xx_zz, g_y_x_0_0_y_x_xx_xx, g_y_x_0_0_y_x_xx_xy, g_y_x_0_0_y_x_xx_xz, g_y_x_0_0_y_x_xx_yy, g_y_x_0_0_y_x_xx_yz, g_y_x_0_0_y_x_xx_zz, g_yy_0_xx_xx, g_yy_0_xx_xy, g_yy_0_xx_xz, g_yy_0_xx_yy, g_yy_0_xx_yz, g_yy_0_xx_zz, g_yy_xx_xx_xx, g_yy_xx_xx_xy, g_yy_xx_xx_xz, g_yy_xx_xx_yy, g_yy_xx_xx_yz, g_yy_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_x_xx_xx[i] = g_0_0_xx_xx[i] - 2.0 * g_0_xx_xx_xx[i] * b_exp - 2.0 * g_yy_0_xx_xx[i] * a_exp + 4.0 * g_yy_xx_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_xx_xy[i] = g_0_0_xx_xy[i] - 2.0 * g_0_xx_xx_xy[i] * b_exp - 2.0 * g_yy_0_xx_xy[i] * a_exp + 4.0 * g_yy_xx_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_xx_xz[i] = g_0_0_xx_xz[i] - 2.0 * g_0_xx_xx_xz[i] * b_exp - 2.0 * g_yy_0_xx_xz[i] * a_exp + 4.0 * g_yy_xx_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_xx_yy[i] = g_0_0_xx_yy[i] - 2.0 * g_0_xx_xx_yy[i] * b_exp - 2.0 * g_yy_0_xx_yy[i] * a_exp + 4.0 * g_yy_xx_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_xx_yz[i] = g_0_0_xx_yz[i] - 2.0 * g_0_xx_xx_yz[i] * b_exp - 2.0 * g_yy_0_xx_yz[i] * a_exp + 4.0 * g_yy_xx_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_xx_zz[i] = g_0_0_xx_zz[i] - 2.0 * g_0_xx_xx_zz[i] * b_exp - 2.0 * g_yy_0_xx_zz[i] * a_exp + 4.0 * g_yy_xx_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1086-1092)

    #pragma omp simd aligned(g_0_0_xy_xx, g_0_0_xy_xy, g_0_0_xy_xz, g_0_0_xy_yy, g_0_0_xy_yz, g_0_0_xy_zz, g_0_xx_xy_xx, g_0_xx_xy_xy, g_0_xx_xy_xz, g_0_xx_xy_yy, g_0_xx_xy_yz, g_0_xx_xy_zz, g_y_x_0_0_y_x_xy_xx, g_y_x_0_0_y_x_xy_xy, g_y_x_0_0_y_x_xy_xz, g_y_x_0_0_y_x_xy_yy, g_y_x_0_0_y_x_xy_yz, g_y_x_0_0_y_x_xy_zz, g_yy_0_xy_xx, g_yy_0_xy_xy, g_yy_0_xy_xz, g_yy_0_xy_yy, g_yy_0_xy_yz, g_yy_0_xy_zz, g_yy_xx_xy_xx, g_yy_xx_xy_xy, g_yy_xx_xy_xz, g_yy_xx_xy_yy, g_yy_xx_xy_yz, g_yy_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_x_xy_xx[i] = g_0_0_xy_xx[i] - 2.0 * g_0_xx_xy_xx[i] * b_exp - 2.0 * g_yy_0_xy_xx[i] * a_exp + 4.0 * g_yy_xx_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_xy_xy[i] = g_0_0_xy_xy[i] - 2.0 * g_0_xx_xy_xy[i] * b_exp - 2.0 * g_yy_0_xy_xy[i] * a_exp + 4.0 * g_yy_xx_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_xy_xz[i] = g_0_0_xy_xz[i] - 2.0 * g_0_xx_xy_xz[i] * b_exp - 2.0 * g_yy_0_xy_xz[i] * a_exp + 4.0 * g_yy_xx_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_xy_yy[i] = g_0_0_xy_yy[i] - 2.0 * g_0_xx_xy_yy[i] * b_exp - 2.0 * g_yy_0_xy_yy[i] * a_exp + 4.0 * g_yy_xx_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_xy_yz[i] = g_0_0_xy_yz[i] - 2.0 * g_0_xx_xy_yz[i] * b_exp - 2.0 * g_yy_0_xy_yz[i] * a_exp + 4.0 * g_yy_xx_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_xy_zz[i] = g_0_0_xy_zz[i] - 2.0 * g_0_xx_xy_zz[i] * b_exp - 2.0 * g_yy_0_xy_zz[i] * a_exp + 4.0 * g_yy_xx_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1092-1098)

    #pragma omp simd aligned(g_0_0_xz_xx, g_0_0_xz_xy, g_0_0_xz_xz, g_0_0_xz_yy, g_0_0_xz_yz, g_0_0_xz_zz, g_0_xx_xz_xx, g_0_xx_xz_xy, g_0_xx_xz_xz, g_0_xx_xz_yy, g_0_xx_xz_yz, g_0_xx_xz_zz, g_y_x_0_0_y_x_xz_xx, g_y_x_0_0_y_x_xz_xy, g_y_x_0_0_y_x_xz_xz, g_y_x_0_0_y_x_xz_yy, g_y_x_0_0_y_x_xz_yz, g_y_x_0_0_y_x_xz_zz, g_yy_0_xz_xx, g_yy_0_xz_xy, g_yy_0_xz_xz, g_yy_0_xz_yy, g_yy_0_xz_yz, g_yy_0_xz_zz, g_yy_xx_xz_xx, g_yy_xx_xz_xy, g_yy_xx_xz_xz, g_yy_xx_xz_yy, g_yy_xx_xz_yz, g_yy_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_x_xz_xx[i] = g_0_0_xz_xx[i] - 2.0 * g_0_xx_xz_xx[i] * b_exp - 2.0 * g_yy_0_xz_xx[i] * a_exp + 4.0 * g_yy_xx_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_xz_xy[i] = g_0_0_xz_xy[i] - 2.0 * g_0_xx_xz_xy[i] * b_exp - 2.0 * g_yy_0_xz_xy[i] * a_exp + 4.0 * g_yy_xx_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_xz_xz[i] = g_0_0_xz_xz[i] - 2.0 * g_0_xx_xz_xz[i] * b_exp - 2.0 * g_yy_0_xz_xz[i] * a_exp + 4.0 * g_yy_xx_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_xz_yy[i] = g_0_0_xz_yy[i] - 2.0 * g_0_xx_xz_yy[i] * b_exp - 2.0 * g_yy_0_xz_yy[i] * a_exp + 4.0 * g_yy_xx_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_xz_yz[i] = g_0_0_xz_yz[i] - 2.0 * g_0_xx_xz_yz[i] * b_exp - 2.0 * g_yy_0_xz_yz[i] * a_exp + 4.0 * g_yy_xx_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_xz_zz[i] = g_0_0_xz_zz[i] - 2.0 * g_0_xx_xz_zz[i] * b_exp - 2.0 * g_yy_0_xz_zz[i] * a_exp + 4.0 * g_yy_xx_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1098-1104)

    #pragma omp simd aligned(g_0_0_yy_xx, g_0_0_yy_xy, g_0_0_yy_xz, g_0_0_yy_yy, g_0_0_yy_yz, g_0_0_yy_zz, g_0_xx_yy_xx, g_0_xx_yy_xy, g_0_xx_yy_xz, g_0_xx_yy_yy, g_0_xx_yy_yz, g_0_xx_yy_zz, g_y_x_0_0_y_x_yy_xx, g_y_x_0_0_y_x_yy_xy, g_y_x_0_0_y_x_yy_xz, g_y_x_0_0_y_x_yy_yy, g_y_x_0_0_y_x_yy_yz, g_y_x_0_0_y_x_yy_zz, g_yy_0_yy_xx, g_yy_0_yy_xy, g_yy_0_yy_xz, g_yy_0_yy_yy, g_yy_0_yy_yz, g_yy_0_yy_zz, g_yy_xx_yy_xx, g_yy_xx_yy_xy, g_yy_xx_yy_xz, g_yy_xx_yy_yy, g_yy_xx_yy_yz, g_yy_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_x_yy_xx[i] = g_0_0_yy_xx[i] - 2.0 * g_0_xx_yy_xx[i] * b_exp - 2.0 * g_yy_0_yy_xx[i] * a_exp + 4.0 * g_yy_xx_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_yy_xy[i] = g_0_0_yy_xy[i] - 2.0 * g_0_xx_yy_xy[i] * b_exp - 2.0 * g_yy_0_yy_xy[i] * a_exp + 4.0 * g_yy_xx_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_yy_xz[i] = g_0_0_yy_xz[i] - 2.0 * g_0_xx_yy_xz[i] * b_exp - 2.0 * g_yy_0_yy_xz[i] * a_exp + 4.0 * g_yy_xx_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_yy_yy[i] = g_0_0_yy_yy[i] - 2.0 * g_0_xx_yy_yy[i] * b_exp - 2.0 * g_yy_0_yy_yy[i] * a_exp + 4.0 * g_yy_xx_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_yy_yz[i] = g_0_0_yy_yz[i] - 2.0 * g_0_xx_yy_yz[i] * b_exp - 2.0 * g_yy_0_yy_yz[i] * a_exp + 4.0 * g_yy_xx_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_yy_zz[i] = g_0_0_yy_zz[i] - 2.0 * g_0_xx_yy_zz[i] * b_exp - 2.0 * g_yy_0_yy_zz[i] * a_exp + 4.0 * g_yy_xx_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1104-1110)

    #pragma omp simd aligned(g_0_0_yz_xx, g_0_0_yz_xy, g_0_0_yz_xz, g_0_0_yz_yy, g_0_0_yz_yz, g_0_0_yz_zz, g_0_xx_yz_xx, g_0_xx_yz_xy, g_0_xx_yz_xz, g_0_xx_yz_yy, g_0_xx_yz_yz, g_0_xx_yz_zz, g_y_x_0_0_y_x_yz_xx, g_y_x_0_0_y_x_yz_xy, g_y_x_0_0_y_x_yz_xz, g_y_x_0_0_y_x_yz_yy, g_y_x_0_0_y_x_yz_yz, g_y_x_0_0_y_x_yz_zz, g_yy_0_yz_xx, g_yy_0_yz_xy, g_yy_0_yz_xz, g_yy_0_yz_yy, g_yy_0_yz_yz, g_yy_0_yz_zz, g_yy_xx_yz_xx, g_yy_xx_yz_xy, g_yy_xx_yz_xz, g_yy_xx_yz_yy, g_yy_xx_yz_yz, g_yy_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_x_yz_xx[i] = g_0_0_yz_xx[i] - 2.0 * g_0_xx_yz_xx[i] * b_exp - 2.0 * g_yy_0_yz_xx[i] * a_exp + 4.0 * g_yy_xx_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_yz_xy[i] = g_0_0_yz_xy[i] - 2.0 * g_0_xx_yz_xy[i] * b_exp - 2.0 * g_yy_0_yz_xy[i] * a_exp + 4.0 * g_yy_xx_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_yz_xz[i] = g_0_0_yz_xz[i] - 2.0 * g_0_xx_yz_xz[i] * b_exp - 2.0 * g_yy_0_yz_xz[i] * a_exp + 4.0 * g_yy_xx_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_yz_yy[i] = g_0_0_yz_yy[i] - 2.0 * g_0_xx_yz_yy[i] * b_exp - 2.0 * g_yy_0_yz_yy[i] * a_exp + 4.0 * g_yy_xx_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_yz_yz[i] = g_0_0_yz_yz[i] - 2.0 * g_0_xx_yz_yz[i] * b_exp - 2.0 * g_yy_0_yz_yz[i] * a_exp + 4.0 * g_yy_xx_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_yz_zz[i] = g_0_0_yz_zz[i] - 2.0 * g_0_xx_yz_zz[i] * b_exp - 2.0 * g_yy_0_yz_zz[i] * a_exp + 4.0 * g_yy_xx_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1110-1116)

    #pragma omp simd aligned(g_0_0_zz_xx, g_0_0_zz_xy, g_0_0_zz_xz, g_0_0_zz_yy, g_0_0_zz_yz, g_0_0_zz_zz, g_0_xx_zz_xx, g_0_xx_zz_xy, g_0_xx_zz_xz, g_0_xx_zz_yy, g_0_xx_zz_yz, g_0_xx_zz_zz, g_y_x_0_0_y_x_zz_xx, g_y_x_0_0_y_x_zz_xy, g_y_x_0_0_y_x_zz_xz, g_y_x_0_0_y_x_zz_yy, g_y_x_0_0_y_x_zz_yz, g_y_x_0_0_y_x_zz_zz, g_yy_0_zz_xx, g_yy_0_zz_xy, g_yy_0_zz_xz, g_yy_0_zz_yy, g_yy_0_zz_yz, g_yy_0_zz_zz, g_yy_xx_zz_xx, g_yy_xx_zz_xy, g_yy_xx_zz_xz, g_yy_xx_zz_yy, g_yy_xx_zz_yz, g_yy_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_x_zz_xx[i] = g_0_0_zz_xx[i] - 2.0 * g_0_xx_zz_xx[i] * b_exp - 2.0 * g_yy_0_zz_xx[i] * a_exp + 4.0 * g_yy_xx_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_zz_xy[i] = g_0_0_zz_xy[i] - 2.0 * g_0_xx_zz_xy[i] * b_exp - 2.0 * g_yy_0_zz_xy[i] * a_exp + 4.0 * g_yy_xx_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_zz_xz[i] = g_0_0_zz_xz[i] - 2.0 * g_0_xx_zz_xz[i] * b_exp - 2.0 * g_yy_0_zz_xz[i] * a_exp + 4.0 * g_yy_xx_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_zz_yy[i] = g_0_0_zz_yy[i] - 2.0 * g_0_xx_zz_yy[i] * b_exp - 2.0 * g_yy_0_zz_yy[i] * a_exp + 4.0 * g_yy_xx_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_zz_yz[i] = g_0_0_zz_yz[i] - 2.0 * g_0_xx_zz_yz[i] * b_exp - 2.0 * g_yy_0_zz_yz[i] * a_exp + 4.0 * g_yy_xx_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_zz_zz[i] = g_0_0_zz_zz[i] - 2.0 * g_0_xx_zz_zz[i] * b_exp - 2.0 * g_yy_0_zz_zz[i] * a_exp + 4.0 * g_yy_xx_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1116-1122)

    #pragma omp simd aligned(g_0_xy_xx_xx, g_0_xy_xx_xy, g_0_xy_xx_xz, g_0_xy_xx_yy, g_0_xy_xx_yz, g_0_xy_xx_zz, g_y_x_0_0_y_y_xx_xx, g_y_x_0_0_y_y_xx_xy, g_y_x_0_0_y_y_xx_xz, g_y_x_0_0_y_y_xx_yy, g_y_x_0_0_y_y_xx_yz, g_y_x_0_0_y_y_xx_zz, g_yy_xy_xx_xx, g_yy_xy_xx_xy, g_yy_xy_xx_xz, g_yy_xy_xx_yy, g_yy_xy_xx_yz, g_yy_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_y_xx_xx[i] = -2.0 * g_0_xy_xx_xx[i] * b_exp + 4.0 * g_yy_xy_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_xx_xy[i] = -2.0 * g_0_xy_xx_xy[i] * b_exp + 4.0 * g_yy_xy_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_xx_xz[i] = -2.0 * g_0_xy_xx_xz[i] * b_exp + 4.0 * g_yy_xy_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_xx_yy[i] = -2.0 * g_0_xy_xx_yy[i] * b_exp + 4.0 * g_yy_xy_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_xx_yz[i] = -2.0 * g_0_xy_xx_yz[i] * b_exp + 4.0 * g_yy_xy_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_xx_zz[i] = -2.0 * g_0_xy_xx_zz[i] * b_exp + 4.0 * g_yy_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1122-1128)

    #pragma omp simd aligned(g_0_xy_xy_xx, g_0_xy_xy_xy, g_0_xy_xy_xz, g_0_xy_xy_yy, g_0_xy_xy_yz, g_0_xy_xy_zz, g_y_x_0_0_y_y_xy_xx, g_y_x_0_0_y_y_xy_xy, g_y_x_0_0_y_y_xy_xz, g_y_x_0_0_y_y_xy_yy, g_y_x_0_0_y_y_xy_yz, g_y_x_0_0_y_y_xy_zz, g_yy_xy_xy_xx, g_yy_xy_xy_xy, g_yy_xy_xy_xz, g_yy_xy_xy_yy, g_yy_xy_xy_yz, g_yy_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_y_xy_xx[i] = -2.0 * g_0_xy_xy_xx[i] * b_exp + 4.0 * g_yy_xy_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_xy_xy[i] = -2.0 * g_0_xy_xy_xy[i] * b_exp + 4.0 * g_yy_xy_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_xy_xz[i] = -2.0 * g_0_xy_xy_xz[i] * b_exp + 4.0 * g_yy_xy_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_xy_yy[i] = -2.0 * g_0_xy_xy_yy[i] * b_exp + 4.0 * g_yy_xy_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_xy_yz[i] = -2.0 * g_0_xy_xy_yz[i] * b_exp + 4.0 * g_yy_xy_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_xy_zz[i] = -2.0 * g_0_xy_xy_zz[i] * b_exp + 4.0 * g_yy_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1128-1134)

    #pragma omp simd aligned(g_0_xy_xz_xx, g_0_xy_xz_xy, g_0_xy_xz_xz, g_0_xy_xz_yy, g_0_xy_xz_yz, g_0_xy_xz_zz, g_y_x_0_0_y_y_xz_xx, g_y_x_0_0_y_y_xz_xy, g_y_x_0_0_y_y_xz_xz, g_y_x_0_0_y_y_xz_yy, g_y_x_0_0_y_y_xz_yz, g_y_x_0_0_y_y_xz_zz, g_yy_xy_xz_xx, g_yy_xy_xz_xy, g_yy_xy_xz_xz, g_yy_xy_xz_yy, g_yy_xy_xz_yz, g_yy_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_y_xz_xx[i] = -2.0 * g_0_xy_xz_xx[i] * b_exp + 4.0 * g_yy_xy_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_xz_xy[i] = -2.0 * g_0_xy_xz_xy[i] * b_exp + 4.0 * g_yy_xy_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_xz_xz[i] = -2.0 * g_0_xy_xz_xz[i] * b_exp + 4.0 * g_yy_xy_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_xz_yy[i] = -2.0 * g_0_xy_xz_yy[i] * b_exp + 4.0 * g_yy_xy_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_xz_yz[i] = -2.0 * g_0_xy_xz_yz[i] * b_exp + 4.0 * g_yy_xy_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_xz_zz[i] = -2.0 * g_0_xy_xz_zz[i] * b_exp + 4.0 * g_yy_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1134-1140)

    #pragma omp simd aligned(g_0_xy_yy_xx, g_0_xy_yy_xy, g_0_xy_yy_xz, g_0_xy_yy_yy, g_0_xy_yy_yz, g_0_xy_yy_zz, g_y_x_0_0_y_y_yy_xx, g_y_x_0_0_y_y_yy_xy, g_y_x_0_0_y_y_yy_xz, g_y_x_0_0_y_y_yy_yy, g_y_x_0_0_y_y_yy_yz, g_y_x_0_0_y_y_yy_zz, g_yy_xy_yy_xx, g_yy_xy_yy_xy, g_yy_xy_yy_xz, g_yy_xy_yy_yy, g_yy_xy_yy_yz, g_yy_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_y_yy_xx[i] = -2.0 * g_0_xy_yy_xx[i] * b_exp + 4.0 * g_yy_xy_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_yy_xy[i] = -2.0 * g_0_xy_yy_xy[i] * b_exp + 4.0 * g_yy_xy_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_yy_xz[i] = -2.0 * g_0_xy_yy_xz[i] * b_exp + 4.0 * g_yy_xy_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_yy_yy[i] = -2.0 * g_0_xy_yy_yy[i] * b_exp + 4.0 * g_yy_xy_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_yy_yz[i] = -2.0 * g_0_xy_yy_yz[i] * b_exp + 4.0 * g_yy_xy_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_yy_zz[i] = -2.0 * g_0_xy_yy_zz[i] * b_exp + 4.0 * g_yy_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1140-1146)

    #pragma omp simd aligned(g_0_xy_yz_xx, g_0_xy_yz_xy, g_0_xy_yz_xz, g_0_xy_yz_yy, g_0_xy_yz_yz, g_0_xy_yz_zz, g_y_x_0_0_y_y_yz_xx, g_y_x_0_0_y_y_yz_xy, g_y_x_0_0_y_y_yz_xz, g_y_x_0_0_y_y_yz_yy, g_y_x_0_0_y_y_yz_yz, g_y_x_0_0_y_y_yz_zz, g_yy_xy_yz_xx, g_yy_xy_yz_xy, g_yy_xy_yz_xz, g_yy_xy_yz_yy, g_yy_xy_yz_yz, g_yy_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_y_yz_xx[i] = -2.0 * g_0_xy_yz_xx[i] * b_exp + 4.0 * g_yy_xy_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_yz_xy[i] = -2.0 * g_0_xy_yz_xy[i] * b_exp + 4.0 * g_yy_xy_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_yz_xz[i] = -2.0 * g_0_xy_yz_xz[i] * b_exp + 4.0 * g_yy_xy_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_yz_yy[i] = -2.0 * g_0_xy_yz_yy[i] * b_exp + 4.0 * g_yy_xy_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_yz_yz[i] = -2.0 * g_0_xy_yz_yz[i] * b_exp + 4.0 * g_yy_xy_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_yz_zz[i] = -2.0 * g_0_xy_yz_zz[i] * b_exp + 4.0 * g_yy_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1146-1152)

    #pragma omp simd aligned(g_0_xy_zz_xx, g_0_xy_zz_xy, g_0_xy_zz_xz, g_0_xy_zz_yy, g_0_xy_zz_yz, g_0_xy_zz_zz, g_y_x_0_0_y_y_zz_xx, g_y_x_0_0_y_y_zz_xy, g_y_x_0_0_y_y_zz_xz, g_y_x_0_0_y_y_zz_yy, g_y_x_0_0_y_y_zz_yz, g_y_x_0_0_y_y_zz_zz, g_yy_xy_zz_xx, g_yy_xy_zz_xy, g_yy_xy_zz_xz, g_yy_xy_zz_yy, g_yy_xy_zz_yz, g_yy_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_y_zz_xx[i] = -2.0 * g_0_xy_zz_xx[i] * b_exp + 4.0 * g_yy_xy_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_zz_xy[i] = -2.0 * g_0_xy_zz_xy[i] * b_exp + 4.0 * g_yy_xy_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_zz_xz[i] = -2.0 * g_0_xy_zz_xz[i] * b_exp + 4.0 * g_yy_xy_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_zz_yy[i] = -2.0 * g_0_xy_zz_yy[i] * b_exp + 4.0 * g_yy_xy_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_zz_yz[i] = -2.0 * g_0_xy_zz_yz[i] * b_exp + 4.0 * g_yy_xy_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_zz_zz[i] = -2.0 * g_0_xy_zz_zz[i] * b_exp + 4.0 * g_yy_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1152-1158)

    #pragma omp simd aligned(g_0_xz_xx_xx, g_0_xz_xx_xy, g_0_xz_xx_xz, g_0_xz_xx_yy, g_0_xz_xx_yz, g_0_xz_xx_zz, g_y_x_0_0_y_z_xx_xx, g_y_x_0_0_y_z_xx_xy, g_y_x_0_0_y_z_xx_xz, g_y_x_0_0_y_z_xx_yy, g_y_x_0_0_y_z_xx_yz, g_y_x_0_0_y_z_xx_zz, g_yy_xz_xx_xx, g_yy_xz_xx_xy, g_yy_xz_xx_xz, g_yy_xz_xx_yy, g_yy_xz_xx_yz, g_yy_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_z_xx_xx[i] = -2.0 * g_0_xz_xx_xx[i] * b_exp + 4.0 * g_yy_xz_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_xx_xy[i] = -2.0 * g_0_xz_xx_xy[i] * b_exp + 4.0 * g_yy_xz_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_xx_xz[i] = -2.0 * g_0_xz_xx_xz[i] * b_exp + 4.0 * g_yy_xz_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_xx_yy[i] = -2.0 * g_0_xz_xx_yy[i] * b_exp + 4.0 * g_yy_xz_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_xx_yz[i] = -2.0 * g_0_xz_xx_yz[i] * b_exp + 4.0 * g_yy_xz_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_xx_zz[i] = -2.0 * g_0_xz_xx_zz[i] * b_exp + 4.0 * g_yy_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1158-1164)

    #pragma omp simd aligned(g_0_xz_xy_xx, g_0_xz_xy_xy, g_0_xz_xy_xz, g_0_xz_xy_yy, g_0_xz_xy_yz, g_0_xz_xy_zz, g_y_x_0_0_y_z_xy_xx, g_y_x_0_0_y_z_xy_xy, g_y_x_0_0_y_z_xy_xz, g_y_x_0_0_y_z_xy_yy, g_y_x_0_0_y_z_xy_yz, g_y_x_0_0_y_z_xy_zz, g_yy_xz_xy_xx, g_yy_xz_xy_xy, g_yy_xz_xy_xz, g_yy_xz_xy_yy, g_yy_xz_xy_yz, g_yy_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_z_xy_xx[i] = -2.0 * g_0_xz_xy_xx[i] * b_exp + 4.0 * g_yy_xz_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_xy_xy[i] = -2.0 * g_0_xz_xy_xy[i] * b_exp + 4.0 * g_yy_xz_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_xy_xz[i] = -2.0 * g_0_xz_xy_xz[i] * b_exp + 4.0 * g_yy_xz_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_xy_yy[i] = -2.0 * g_0_xz_xy_yy[i] * b_exp + 4.0 * g_yy_xz_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_xy_yz[i] = -2.0 * g_0_xz_xy_yz[i] * b_exp + 4.0 * g_yy_xz_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_xy_zz[i] = -2.0 * g_0_xz_xy_zz[i] * b_exp + 4.0 * g_yy_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1164-1170)

    #pragma omp simd aligned(g_0_xz_xz_xx, g_0_xz_xz_xy, g_0_xz_xz_xz, g_0_xz_xz_yy, g_0_xz_xz_yz, g_0_xz_xz_zz, g_y_x_0_0_y_z_xz_xx, g_y_x_0_0_y_z_xz_xy, g_y_x_0_0_y_z_xz_xz, g_y_x_0_0_y_z_xz_yy, g_y_x_0_0_y_z_xz_yz, g_y_x_0_0_y_z_xz_zz, g_yy_xz_xz_xx, g_yy_xz_xz_xy, g_yy_xz_xz_xz, g_yy_xz_xz_yy, g_yy_xz_xz_yz, g_yy_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_z_xz_xx[i] = -2.0 * g_0_xz_xz_xx[i] * b_exp + 4.0 * g_yy_xz_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_xz_xy[i] = -2.0 * g_0_xz_xz_xy[i] * b_exp + 4.0 * g_yy_xz_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_xz_xz[i] = -2.0 * g_0_xz_xz_xz[i] * b_exp + 4.0 * g_yy_xz_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_xz_yy[i] = -2.0 * g_0_xz_xz_yy[i] * b_exp + 4.0 * g_yy_xz_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_xz_yz[i] = -2.0 * g_0_xz_xz_yz[i] * b_exp + 4.0 * g_yy_xz_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_xz_zz[i] = -2.0 * g_0_xz_xz_zz[i] * b_exp + 4.0 * g_yy_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1170-1176)

    #pragma omp simd aligned(g_0_xz_yy_xx, g_0_xz_yy_xy, g_0_xz_yy_xz, g_0_xz_yy_yy, g_0_xz_yy_yz, g_0_xz_yy_zz, g_y_x_0_0_y_z_yy_xx, g_y_x_0_0_y_z_yy_xy, g_y_x_0_0_y_z_yy_xz, g_y_x_0_0_y_z_yy_yy, g_y_x_0_0_y_z_yy_yz, g_y_x_0_0_y_z_yy_zz, g_yy_xz_yy_xx, g_yy_xz_yy_xy, g_yy_xz_yy_xz, g_yy_xz_yy_yy, g_yy_xz_yy_yz, g_yy_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_z_yy_xx[i] = -2.0 * g_0_xz_yy_xx[i] * b_exp + 4.0 * g_yy_xz_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_yy_xy[i] = -2.0 * g_0_xz_yy_xy[i] * b_exp + 4.0 * g_yy_xz_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_yy_xz[i] = -2.0 * g_0_xz_yy_xz[i] * b_exp + 4.0 * g_yy_xz_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_yy_yy[i] = -2.0 * g_0_xz_yy_yy[i] * b_exp + 4.0 * g_yy_xz_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_yy_yz[i] = -2.0 * g_0_xz_yy_yz[i] * b_exp + 4.0 * g_yy_xz_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_yy_zz[i] = -2.0 * g_0_xz_yy_zz[i] * b_exp + 4.0 * g_yy_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1176-1182)

    #pragma omp simd aligned(g_0_xz_yz_xx, g_0_xz_yz_xy, g_0_xz_yz_xz, g_0_xz_yz_yy, g_0_xz_yz_yz, g_0_xz_yz_zz, g_y_x_0_0_y_z_yz_xx, g_y_x_0_0_y_z_yz_xy, g_y_x_0_0_y_z_yz_xz, g_y_x_0_0_y_z_yz_yy, g_y_x_0_0_y_z_yz_yz, g_y_x_0_0_y_z_yz_zz, g_yy_xz_yz_xx, g_yy_xz_yz_xy, g_yy_xz_yz_xz, g_yy_xz_yz_yy, g_yy_xz_yz_yz, g_yy_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_z_yz_xx[i] = -2.0 * g_0_xz_yz_xx[i] * b_exp + 4.0 * g_yy_xz_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_yz_xy[i] = -2.0 * g_0_xz_yz_xy[i] * b_exp + 4.0 * g_yy_xz_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_yz_xz[i] = -2.0 * g_0_xz_yz_xz[i] * b_exp + 4.0 * g_yy_xz_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_yz_yy[i] = -2.0 * g_0_xz_yz_yy[i] * b_exp + 4.0 * g_yy_xz_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_yz_yz[i] = -2.0 * g_0_xz_yz_yz[i] * b_exp + 4.0 * g_yy_xz_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_yz_zz[i] = -2.0 * g_0_xz_yz_zz[i] * b_exp + 4.0 * g_yy_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1182-1188)

    #pragma omp simd aligned(g_0_xz_zz_xx, g_0_xz_zz_xy, g_0_xz_zz_xz, g_0_xz_zz_yy, g_0_xz_zz_yz, g_0_xz_zz_zz, g_y_x_0_0_y_z_zz_xx, g_y_x_0_0_y_z_zz_xy, g_y_x_0_0_y_z_zz_xz, g_y_x_0_0_y_z_zz_yy, g_y_x_0_0_y_z_zz_yz, g_y_x_0_0_y_z_zz_zz, g_yy_xz_zz_xx, g_yy_xz_zz_xy, g_yy_xz_zz_xz, g_yy_xz_zz_yy, g_yy_xz_zz_yz, g_yy_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_z_zz_xx[i] = -2.0 * g_0_xz_zz_xx[i] * b_exp + 4.0 * g_yy_xz_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_zz_xy[i] = -2.0 * g_0_xz_zz_xy[i] * b_exp + 4.0 * g_yy_xz_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_zz_xz[i] = -2.0 * g_0_xz_zz_xz[i] * b_exp + 4.0 * g_yy_xz_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_zz_yy[i] = -2.0 * g_0_xz_zz_yy[i] * b_exp + 4.0 * g_yy_xz_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_zz_yz[i] = -2.0 * g_0_xz_zz_yz[i] * b_exp + 4.0 * g_yy_xz_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_zz_zz[i] = -2.0 * g_0_xz_zz_zz[i] * b_exp + 4.0 * g_yy_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1188-1194)

    #pragma omp simd aligned(g_y_x_0_0_z_x_xx_xx, g_y_x_0_0_z_x_xx_xy, g_y_x_0_0_z_x_xx_xz, g_y_x_0_0_z_x_xx_yy, g_y_x_0_0_z_x_xx_yz, g_y_x_0_0_z_x_xx_zz, g_yz_0_xx_xx, g_yz_0_xx_xy, g_yz_0_xx_xz, g_yz_0_xx_yy, g_yz_0_xx_yz, g_yz_0_xx_zz, g_yz_xx_xx_xx, g_yz_xx_xx_xy, g_yz_xx_xx_xz, g_yz_xx_xx_yy, g_yz_xx_xx_yz, g_yz_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_x_xx_xx[i] = -2.0 * g_yz_0_xx_xx[i] * a_exp + 4.0 * g_yz_xx_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_xx_xy[i] = -2.0 * g_yz_0_xx_xy[i] * a_exp + 4.0 * g_yz_xx_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_xx_xz[i] = -2.0 * g_yz_0_xx_xz[i] * a_exp + 4.0 * g_yz_xx_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_xx_yy[i] = -2.0 * g_yz_0_xx_yy[i] * a_exp + 4.0 * g_yz_xx_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_xx_yz[i] = -2.0 * g_yz_0_xx_yz[i] * a_exp + 4.0 * g_yz_xx_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_xx_zz[i] = -2.0 * g_yz_0_xx_zz[i] * a_exp + 4.0 * g_yz_xx_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1194-1200)

    #pragma omp simd aligned(g_y_x_0_0_z_x_xy_xx, g_y_x_0_0_z_x_xy_xy, g_y_x_0_0_z_x_xy_xz, g_y_x_0_0_z_x_xy_yy, g_y_x_0_0_z_x_xy_yz, g_y_x_0_0_z_x_xy_zz, g_yz_0_xy_xx, g_yz_0_xy_xy, g_yz_0_xy_xz, g_yz_0_xy_yy, g_yz_0_xy_yz, g_yz_0_xy_zz, g_yz_xx_xy_xx, g_yz_xx_xy_xy, g_yz_xx_xy_xz, g_yz_xx_xy_yy, g_yz_xx_xy_yz, g_yz_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_x_xy_xx[i] = -2.0 * g_yz_0_xy_xx[i] * a_exp + 4.0 * g_yz_xx_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_xy_xy[i] = -2.0 * g_yz_0_xy_xy[i] * a_exp + 4.0 * g_yz_xx_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_xy_xz[i] = -2.0 * g_yz_0_xy_xz[i] * a_exp + 4.0 * g_yz_xx_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_xy_yy[i] = -2.0 * g_yz_0_xy_yy[i] * a_exp + 4.0 * g_yz_xx_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_xy_yz[i] = -2.0 * g_yz_0_xy_yz[i] * a_exp + 4.0 * g_yz_xx_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_xy_zz[i] = -2.0 * g_yz_0_xy_zz[i] * a_exp + 4.0 * g_yz_xx_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1200-1206)

    #pragma omp simd aligned(g_y_x_0_0_z_x_xz_xx, g_y_x_0_0_z_x_xz_xy, g_y_x_0_0_z_x_xz_xz, g_y_x_0_0_z_x_xz_yy, g_y_x_0_0_z_x_xz_yz, g_y_x_0_0_z_x_xz_zz, g_yz_0_xz_xx, g_yz_0_xz_xy, g_yz_0_xz_xz, g_yz_0_xz_yy, g_yz_0_xz_yz, g_yz_0_xz_zz, g_yz_xx_xz_xx, g_yz_xx_xz_xy, g_yz_xx_xz_xz, g_yz_xx_xz_yy, g_yz_xx_xz_yz, g_yz_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_x_xz_xx[i] = -2.0 * g_yz_0_xz_xx[i] * a_exp + 4.0 * g_yz_xx_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_xz_xy[i] = -2.0 * g_yz_0_xz_xy[i] * a_exp + 4.0 * g_yz_xx_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_xz_xz[i] = -2.0 * g_yz_0_xz_xz[i] * a_exp + 4.0 * g_yz_xx_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_xz_yy[i] = -2.0 * g_yz_0_xz_yy[i] * a_exp + 4.0 * g_yz_xx_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_xz_yz[i] = -2.0 * g_yz_0_xz_yz[i] * a_exp + 4.0 * g_yz_xx_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_xz_zz[i] = -2.0 * g_yz_0_xz_zz[i] * a_exp + 4.0 * g_yz_xx_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1206-1212)

    #pragma omp simd aligned(g_y_x_0_0_z_x_yy_xx, g_y_x_0_0_z_x_yy_xy, g_y_x_0_0_z_x_yy_xz, g_y_x_0_0_z_x_yy_yy, g_y_x_0_0_z_x_yy_yz, g_y_x_0_0_z_x_yy_zz, g_yz_0_yy_xx, g_yz_0_yy_xy, g_yz_0_yy_xz, g_yz_0_yy_yy, g_yz_0_yy_yz, g_yz_0_yy_zz, g_yz_xx_yy_xx, g_yz_xx_yy_xy, g_yz_xx_yy_xz, g_yz_xx_yy_yy, g_yz_xx_yy_yz, g_yz_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_x_yy_xx[i] = -2.0 * g_yz_0_yy_xx[i] * a_exp + 4.0 * g_yz_xx_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_yy_xy[i] = -2.0 * g_yz_0_yy_xy[i] * a_exp + 4.0 * g_yz_xx_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_yy_xz[i] = -2.0 * g_yz_0_yy_xz[i] * a_exp + 4.0 * g_yz_xx_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_yy_yy[i] = -2.0 * g_yz_0_yy_yy[i] * a_exp + 4.0 * g_yz_xx_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_yy_yz[i] = -2.0 * g_yz_0_yy_yz[i] * a_exp + 4.0 * g_yz_xx_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_yy_zz[i] = -2.0 * g_yz_0_yy_zz[i] * a_exp + 4.0 * g_yz_xx_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1212-1218)

    #pragma omp simd aligned(g_y_x_0_0_z_x_yz_xx, g_y_x_0_0_z_x_yz_xy, g_y_x_0_0_z_x_yz_xz, g_y_x_0_0_z_x_yz_yy, g_y_x_0_0_z_x_yz_yz, g_y_x_0_0_z_x_yz_zz, g_yz_0_yz_xx, g_yz_0_yz_xy, g_yz_0_yz_xz, g_yz_0_yz_yy, g_yz_0_yz_yz, g_yz_0_yz_zz, g_yz_xx_yz_xx, g_yz_xx_yz_xy, g_yz_xx_yz_xz, g_yz_xx_yz_yy, g_yz_xx_yz_yz, g_yz_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_x_yz_xx[i] = -2.0 * g_yz_0_yz_xx[i] * a_exp + 4.0 * g_yz_xx_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_yz_xy[i] = -2.0 * g_yz_0_yz_xy[i] * a_exp + 4.0 * g_yz_xx_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_yz_xz[i] = -2.0 * g_yz_0_yz_xz[i] * a_exp + 4.0 * g_yz_xx_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_yz_yy[i] = -2.0 * g_yz_0_yz_yy[i] * a_exp + 4.0 * g_yz_xx_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_yz_yz[i] = -2.0 * g_yz_0_yz_yz[i] * a_exp + 4.0 * g_yz_xx_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_yz_zz[i] = -2.0 * g_yz_0_yz_zz[i] * a_exp + 4.0 * g_yz_xx_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1218-1224)

    #pragma omp simd aligned(g_y_x_0_0_z_x_zz_xx, g_y_x_0_0_z_x_zz_xy, g_y_x_0_0_z_x_zz_xz, g_y_x_0_0_z_x_zz_yy, g_y_x_0_0_z_x_zz_yz, g_y_x_0_0_z_x_zz_zz, g_yz_0_zz_xx, g_yz_0_zz_xy, g_yz_0_zz_xz, g_yz_0_zz_yy, g_yz_0_zz_yz, g_yz_0_zz_zz, g_yz_xx_zz_xx, g_yz_xx_zz_xy, g_yz_xx_zz_xz, g_yz_xx_zz_yy, g_yz_xx_zz_yz, g_yz_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_x_zz_xx[i] = -2.0 * g_yz_0_zz_xx[i] * a_exp + 4.0 * g_yz_xx_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_zz_xy[i] = -2.0 * g_yz_0_zz_xy[i] * a_exp + 4.0 * g_yz_xx_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_zz_xz[i] = -2.0 * g_yz_0_zz_xz[i] * a_exp + 4.0 * g_yz_xx_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_zz_yy[i] = -2.0 * g_yz_0_zz_yy[i] * a_exp + 4.0 * g_yz_xx_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_zz_yz[i] = -2.0 * g_yz_0_zz_yz[i] * a_exp + 4.0 * g_yz_xx_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_zz_zz[i] = -2.0 * g_yz_0_zz_zz[i] * a_exp + 4.0 * g_yz_xx_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1224-1230)

    #pragma omp simd aligned(g_y_x_0_0_z_y_xx_xx, g_y_x_0_0_z_y_xx_xy, g_y_x_0_0_z_y_xx_xz, g_y_x_0_0_z_y_xx_yy, g_y_x_0_0_z_y_xx_yz, g_y_x_0_0_z_y_xx_zz, g_yz_xy_xx_xx, g_yz_xy_xx_xy, g_yz_xy_xx_xz, g_yz_xy_xx_yy, g_yz_xy_xx_yz, g_yz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_y_xx_xx[i] = 4.0 * g_yz_xy_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_xx_xy[i] = 4.0 * g_yz_xy_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_xx_xz[i] = 4.0 * g_yz_xy_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_xx_yy[i] = 4.0 * g_yz_xy_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_xx_yz[i] = 4.0 * g_yz_xy_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_xx_zz[i] = 4.0 * g_yz_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1230-1236)

    #pragma omp simd aligned(g_y_x_0_0_z_y_xy_xx, g_y_x_0_0_z_y_xy_xy, g_y_x_0_0_z_y_xy_xz, g_y_x_0_0_z_y_xy_yy, g_y_x_0_0_z_y_xy_yz, g_y_x_0_0_z_y_xy_zz, g_yz_xy_xy_xx, g_yz_xy_xy_xy, g_yz_xy_xy_xz, g_yz_xy_xy_yy, g_yz_xy_xy_yz, g_yz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_y_xy_xx[i] = 4.0 * g_yz_xy_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_xy_xy[i] = 4.0 * g_yz_xy_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_xy_xz[i] = 4.0 * g_yz_xy_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_xy_yy[i] = 4.0 * g_yz_xy_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_xy_yz[i] = 4.0 * g_yz_xy_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_xy_zz[i] = 4.0 * g_yz_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1236-1242)

    #pragma omp simd aligned(g_y_x_0_0_z_y_xz_xx, g_y_x_0_0_z_y_xz_xy, g_y_x_0_0_z_y_xz_xz, g_y_x_0_0_z_y_xz_yy, g_y_x_0_0_z_y_xz_yz, g_y_x_0_0_z_y_xz_zz, g_yz_xy_xz_xx, g_yz_xy_xz_xy, g_yz_xy_xz_xz, g_yz_xy_xz_yy, g_yz_xy_xz_yz, g_yz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_y_xz_xx[i] = 4.0 * g_yz_xy_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_xz_xy[i] = 4.0 * g_yz_xy_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_xz_xz[i] = 4.0 * g_yz_xy_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_xz_yy[i] = 4.0 * g_yz_xy_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_xz_yz[i] = 4.0 * g_yz_xy_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_xz_zz[i] = 4.0 * g_yz_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1242-1248)

    #pragma omp simd aligned(g_y_x_0_0_z_y_yy_xx, g_y_x_0_0_z_y_yy_xy, g_y_x_0_0_z_y_yy_xz, g_y_x_0_0_z_y_yy_yy, g_y_x_0_0_z_y_yy_yz, g_y_x_0_0_z_y_yy_zz, g_yz_xy_yy_xx, g_yz_xy_yy_xy, g_yz_xy_yy_xz, g_yz_xy_yy_yy, g_yz_xy_yy_yz, g_yz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_y_yy_xx[i] = 4.0 * g_yz_xy_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_yy_xy[i] = 4.0 * g_yz_xy_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_yy_xz[i] = 4.0 * g_yz_xy_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_yy_yy[i] = 4.0 * g_yz_xy_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_yy_yz[i] = 4.0 * g_yz_xy_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_yy_zz[i] = 4.0 * g_yz_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1248-1254)

    #pragma omp simd aligned(g_y_x_0_0_z_y_yz_xx, g_y_x_0_0_z_y_yz_xy, g_y_x_0_0_z_y_yz_xz, g_y_x_0_0_z_y_yz_yy, g_y_x_0_0_z_y_yz_yz, g_y_x_0_0_z_y_yz_zz, g_yz_xy_yz_xx, g_yz_xy_yz_xy, g_yz_xy_yz_xz, g_yz_xy_yz_yy, g_yz_xy_yz_yz, g_yz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_y_yz_xx[i] = 4.0 * g_yz_xy_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_yz_xy[i] = 4.0 * g_yz_xy_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_yz_xz[i] = 4.0 * g_yz_xy_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_yz_yy[i] = 4.0 * g_yz_xy_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_yz_yz[i] = 4.0 * g_yz_xy_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_yz_zz[i] = 4.0 * g_yz_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1254-1260)

    #pragma omp simd aligned(g_y_x_0_0_z_y_zz_xx, g_y_x_0_0_z_y_zz_xy, g_y_x_0_0_z_y_zz_xz, g_y_x_0_0_z_y_zz_yy, g_y_x_0_0_z_y_zz_yz, g_y_x_0_0_z_y_zz_zz, g_yz_xy_zz_xx, g_yz_xy_zz_xy, g_yz_xy_zz_xz, g_yz_xy_zz_yy, g_yz_xy_zz_yz, g_yz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_y_zz_xx[i] = 4.0 * g_yz_xy_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_zz_xy[i] = 4.0 * g_yz_xy_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_zz_xz[i] = 4.0 * g_yz_xy_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_zz_yy[i] = 4.0 * g_yz_xy_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_zz_yz[i] = 4.0 * g_yz_xy_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_zz_zz[i] = 4.0 * g_yz_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1260-1266)

    #pragma omp simd aligned(g_y_x_0_0_z_z_xx_xx, g_y_x_0_0_z_z_xx_xy, g_y_x_0_0_z_z_xx_xz, g_y_x_0_0_z_z_xx_yy, g_y_x_0_0_z_z_xx_yz, g_y_x_0_0_z_z_xx_zz, g_yz_xz_xx_xx, g_yz_xz_xx_xy, g_yz_xz_xx_xz, g_yz_xz_xx_yy, g_yz_xz_xx_yz, g_yz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_z_xx_xx[i] = 4.0 * g_yz_xz_xx_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_xx_xy[i] = 4.0 * g_yz_xz_xx_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_xx_xz[i] = 4.0 * g_yz_xz_xx_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_xx_yy[i] = 4.0 * g_yz_xz_xx_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_xx_yz[i] = 4.0 * g_yz_xz_xx_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_xx_zz[i] = 4.0 * g_yz_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1266-1272)

    #pragma omp simd aligned(g_y_x_0_0_z_z_xy_xx, g_y_x_0_0_z_z_xy_xy, g_y_x_0_0_z_z_xy_xz, g_y_x_0_0_z_z_xy_yy, g_y_x_0_0_z_z_xy_yz, g_y_x_0_0_z_z_xy_zz, g_yz_xz_xy_xx, g_yz_xz_xy_xy, g_yz_xz_xy_xz, g_yz_xz_xy_yy, g_yz_xz_xy_yz, g_yz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_z_xy_xx[i] = 4.0 * g_yz_xz_xy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_xy_xy[i] = 4.0 * g_yz_xz_xy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_xy_xz[i] = 4.0 * g_yz_xz_xy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_xy_yy[i] = 4.0 * g_yz_xz_xy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_xy_yz[i] = 4.0 * g_yz_xz_xy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_xy_zz[i] = 4.0 * g_yz_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1272-1278)

    #pragma omp simd aligned(g_y_x_0_0_z_z_xz_xx, g_y_x_0_0_z_z_xz_xy, g_y_x_0_0_z_z_xz_xz, g_y_x_0_0_z_z_xz_yy, g_y_x_0_0_z_z_xz_yz, g_y_x_0_0_z_z_xz_zz, g_yz_xz_xz_xx, g_yz_xz_xz_xy, g_yz_xz_xz_xz, g_yz_xz_xz_yy, g_yz_xz_xz_yz, g_yz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_z_xz_xx[i] = 4.0 * g_yz_xz_xz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_xz_xy[i] = 4.0 * g_yz_xz_xz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_xz_xz[i] = 4.0 * g_yz_xz_xz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_xz_yy[i] = 4.0 * g_yz_xz_xz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_xz_yz[i] = 4.0 * g_yz_xz_xz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_xz_zz[i] = 4.0 * g_yz_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1278-1284)

    #pragma omp simd aligned(g_y_x_0_0_z_z_yy_xx, g_y_x_0_0_z_z_yy_xy, g_y_x_0_0_z_z_yy_xz, g_y_x_0_0_z_z_yy_yy, g_y_x_0_0_z_z_yy_yz, g_y_x_0_0_z_z_yy_zz, g_yz_xz_yy_xx, g_yz_xz_yy_xy, g_yz_xz_yy_xz, g_yz_xz_yy_yy, g_yz_xz_yy_yz, g_yz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_z_yy_xx[i] = 4.0 * g_yz_xz_yy_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_yy_xy[i] = 4.0 * g_yz_xz_yy_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_yy_xz[i] = 4.0 * g_yz_xz_yy_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_yy_yy[i] = 4.0 * g_yz_xz_yy_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_yy_yz[i] = 4.0 * g_yz_xz_yy_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_yy_zz[i] = 4.0 * g_yz_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1284-1290)

    #pragma omp simd aligned(g_y_x_0_0_z_z_yz_xx, g_y_x_0_0_z_z_yz_xy, g_y_x_0_0_z_z_yz_xz, g_y_x_0_0_z_z_yz_yy, g_y_x_0_0_z_z_yz_yz, g_y_x_0_0_z_z_yz_zz, g_yz_xz_yz_xx, g_yz_xz_yz_xy, g_yz_xz_yz_xz, g_yz_xz_yz_yy, g_yz_xz_yz_yz, g_yz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_z_yz_xx[i] = 4.0 * g_yz_xz_yz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_yz_xy[i] = 4.0 * g_yz_xz_yz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_yz_xz[i] = 4.0 * g_yz_xz_yz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_yz_yy[i] = 4.0 * g_yz_xz_yz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_yz_yz[i] = 4.0 * g_yz_xz_yz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_yz_zz[i] = 4.0 * g_yz_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1290-1296)

    #pragma omp simd aligned(g_y_x_0_0_z_z_zz_xx, g_y_x_0_0_z_z_zz_xy, g_y_x_0_0_z_z_zz_xz, g_y_x_0_0_z_z_zz_yy, g_y_x_0_0_z_z_zz_yz, g_y_x_0_0_z_z_zz_zz, g_yz_xz_zz_xx, g_yz_xz_zz_xy, g_yz_xz_zz_xz, g_yz_xz_zz_yy, g_yz_xz_zz_yz, g_yz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_z_zz_xx[i] = 4.0 * g_yz_xz_zz_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_zz_xy[i] = 4.0 * g_yz_xz_zz_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_zz_xz[i] = 4.0 * g_yz_xz_zz_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_zz_yy[i] = 4.0 * g_yz_xz_zz_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_zz_yz[i] = 4.0 * g_yz_xz_zz_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_zz_zz[i] = 4.0 * g_yz_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1296-1302)

    #pragma omp simd aligned(g_xy_xy_xx_xx, g_xy_xy_xx_xy, g_xy_xy_xx_xz, g_xy_xy_xx_yy, g_xy_xy_xx_yz, g_xy_xy_xx_zz, g_y_y_0_0_x_x_xx_xx, g_y_y_0_0_x_x_xx_xy, g_y_y_0_0_x_x_xx_xz, g_y_y_0_0_x_x_xx_yy, g_y_y_0_0_x_x_xx_yz, g_y_y_0_0_x_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_x_xx_xx[i] = 4.0 * g_xy_xy_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_xx_xy[i] = 4.0 * g_xy_xy_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_xx_xz[i] = 4.0 * g_xy_xy_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_xx_yy[i] = 4.0 * g_xy_xy_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_xx_yz[i] = 4.0 * g_xy_xy_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_xx_zz[i] = 4.0 * g_xy_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1302-1308)

    #pragma omp simd aligned(g_xy_xy_xy_xx, g_xy_xy_xy_xy, g_xy_xy_xy_xz, g_xy_xy_xy_yy, g_xy_xy_xy_yz, g_xy_xy_xy_zz, g_y_y_0_0_x_x_xy_xx, g_y_y_0_0_x_x_xy_xy, g_y_y_0_0_x_x_xy_xz, g_y_y_0_0_x_x_xy_yy, g_y_y_0_0_x_x_xy_yz, g_y_y_0_0_x_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_x_xy_xx[i] = 4.0 * g_xy_xy_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_xy_xy[i] = 4.0 * g_xy_xy_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_xy_xz[i] = 4.0 * g_xy_xy_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_xy_yy[i] = 4.0 * g_xy_xy_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_xy_yz[i] = 4.0 * g_xy_xy_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_xy_zz[i] = 4.0 * g_xy_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1308-1314)

    #pragma omp simd aligned(g_xy_xy_xz_xx, g_xy_xy_xz_xy, g_xy_xy_xz_xz, g_xy_xy_xz_yy, g_xy_xy_xz_yz, g_xy_xy_xz_zz, g_y_y_0_0_x_x_xz_xx, g_y_y_0_0_x_x_xz_xy, g_y_y_0_0_x_x_xz_xz, g_y_y_0_0_x_x_xz_yy, g_y_y_0_0_x_x_xz_yz, g_y_y_0_0_x_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_x_xz_xx[i] = 4.0 * g_xy_xy_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_xz_xy[i] = 4.0 * g_xy_xy_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_xz_xz[i] = 4.0 * g_xy_xy_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_xz_yy[i] = 4.0 * g_xy_xy_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_xz_yz[i] = 4.0 * g_xy_xy_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_xz_zz[i] = 4.0 * g_xy_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1314-1320)

    #pragma omp simd aligned(g_xy_xy_yy_xx, g_xy_xy_yy_xy, g_xy_xy_yy_xz, g_xy_xy_yy_yy, g_xy_xy_yy_yz, g_xy_xy_yy_zz, g_y_y_0_0_x_x_yy_xx, g_y_y_0_0_x_x_yy_xy, g_y_y_0_0_x_x_yy_xz, g_y_y_0_0_x_x_yy_yy, g_y_y_0_0_x_x_yy_yz, g_y_y_0_0_x_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_x_yy_xx[i] = 4.0 * g_xy_xy_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_yy_xy[i] = 4.0 * g_xy_xy_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_yy_xz[i] = 4.0 * g_xy_xy_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_yy_yy[i] = 4.0 * g_xy_xy_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_yy_yz[i] = 4.0 * g_xy_xy_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_yy_zz[i] = 4.0 * g_xy_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1320-1326)

    #pragma omp simd aligned(g_xy_xy_yz_xx, g_xy_xy_yz_xy, g_xy_xy_yz_xz, g_xy_xy_yz_yy, g_xy_xy_yz_yz, g_xy_xy_yz_zz, g_y_y_0_0_x_x_yz_xx, g_y_y_0_0_x_x_yz_xy, g_y_y_0_0_x_x_yz_xz, g_y_y_0_0_x_x_yz_yy, g_y_y_0_0_x_x_yz_yz, g_y_y_0_0_x_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_x_yz_xx[i] = 4.0 * g_xy_xy_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_yz_xy[i] = 4.0 * g_xy_xy_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_yz_xz[i] = 4.0 * g_xy_xy_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_yz_yy[i] = 4.0 * g_xy_xy_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_yz_yz[i] = 4.0 * g_xy_xy_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_yz_zz[i] = 4.0 * g_xy_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1326-1332)

    #pragma omp simd aligned(g_xy_xy_zz_xx, g_xy_xy_zz_xy, g_xy_xy_zz_xz, g_xy_xy_zz_yy, g_xy_xy_zz_yz, g_xy_xy_zz_zz, g_y_y_0_0_x_x_zz_xx, g_y_y_0_0_x_x_zz_xy, g_y_y_0_0_x_x_zz_xz, g_y_y_0_0_x_x_zz_yy, g_y_y_0_0_x_x_zz_yz, g_y_y_0_0_x_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_x_zz_xx[i] = 4.0 * g_xy_xy_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_zz_xy[i] = 4.0 * g_xy_xy_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_zz_xz[i] = 4.0 * g_xy_xy_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_zz_yy[i] = 4.0 * g_xy_xy_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_zz_yz[i] = 4.0 * g_xy_xy_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_zz_zz[i] = 4.0 * g_xy_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1332-1338)

    #pragma omp simd aligned(g_xy_0_xx_xx, g_xy_0_xx_xy, g_xy_0_xx_xz, g_xy_0_xx_yy, g_xy_0_xx_yz, g_xy_0_xx_zz, g_xy_yy_xx_xx, g_xy_yy_xx_xy, g_xy_yy_xx_xz, g_xy_yy_xx_yy, g_xy_yy_xx_yz, g_xy_yy_xx_zz, g_y_y_0_0_x_y_xx_xx, g_y_y_0_0_x_y_xx_xy, g_y_y_0_0_x_y_xx_xz, g_y_y_0_0_x_y_xx_yy, g_y_y_0_0_x_y_xx_yz, g_y_y_0_0_x_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_y_xx_xx[i] = -2.0 * g_xy_0_xx_xx[i] * a_exp + 4.0 * g_xy_yy_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_xx_xy[i] = -2.0 * g_xy_0_xx_xy[i] * a_exp + 4.0 * g_xy_yy_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_xx_xz[i] = -2.0 * g_xy_0_xx_xz[i] * a_exp + 4.0 * g_xy_yy_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_xx_yy[i] = -2.0 * g_xy_0_xx_yy[i] * a_exp + 4.0 * g_xy_yy_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_xx_yz[i] = -2.0 * g_xy_0_xx_yz[i] * a_exp + 4.0 * g_xy_yy_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_xx_zz[i] = -2.0 * g_xy_0_xx_zz[i] * a_exp + 4.0 * g_xy_yy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1338-1344)

    #pragma omp simd aligned(g_xy_0_xy_xx, g_xy_0_xy_xy, g_xy_0_xy_xz, g_xy_0_xy_yy, g_xy_0_xy_yz, g_xy_0_xy_zz, g_xy_yy_xy_xx, g_xy_yy_xy_xy, g_xy_yy_xy_xz, g_xy_yy_xy_yy, g_xy_yy_xy_yz, g_xy_yy_xy_zz, g_y_y_0_0_x_y_xy_xx, g_y_y_0_0_x_y_xy_xy, g_y_y_0_0_x_y_xy_xz, g_y_y_0_0_x_y_xy_yy, g_y_y_0_0_x_y_xy_yz, g_y_y_0_0_x_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_y_xy_xx[i] = -2.0 * g_xy_0_xy_xx[i] * a_exp + 4.0 * g_xy_yy_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_xy_xy[i] = -2.0 * g_xy_0_xy_xy[i] * a_exp + 4.0 * g_xy_yy_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_xy_xz[i] = -2.0 * g_xy_0_xy_xz[i] * a_exp + 4.0 * g_xy_yy_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_xy_yy[i] = -2.0 * g_xy_0_xy_yy[i] * a_exp + 4.0 * g_xy_yy_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_xy_yz[i] = -2.0 * g_xy_0_xy_yz[i] * a_exp + 4.0 * g_xy_yy_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_xy_zz[i] = -2.0 * g_xy_0_xy_zz[i] * a_exp + 4.0 * g_xy_yy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1344-1350)

    #pragma omp simd aligned(g_xy_0_xz_xx, g_xy_0_xz_xy, g_xy_0_xz_xz, g_xy_0_xz_yy, g_xy_0_xz_yz, g_xy_0_xz_zz, g_xy_yy_xz_xx, g_xy_yy_xz_xy, g_xy_yy_xz_xz, g_xy_yy_xz_yy, g_xy_yy_xz_yz, g_xy_yy_xz_zz, g_y_y_0_0_x_y_xz_xx, g_y_y_0_0_x_y_xz_xy, g_y_y_0_0_x_y_xz_xz, g_y_y_0_0_x_y_xz_yy, g_y_y_0_0_x_y_xz_yz, g_y_y_0_0_x_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_y_xz_xx[i] = -2.0 * g_xy_0_xz_xx[i] * a_exp + 4.0 * g_xy_yy_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_xz_xy[i] = -2.0 * g_xy_0_xz_xy[i] * a_exp + 4.0 * g_xy_yy_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_xz_xz[i] = -2.0 * g_xy_0_xz_xz[i] * a_exp + 4.0 * g_xy_yy_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_xz_yy[i] = -2.0 * g_xy_0_xz_yy[i] * a_exp + 4.0 * g_xy_yy_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_xz_yz[i] = -2.0 * g_xy_0_xz_yz[i] * a_exp + 4.0 * g_xy_yy_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_xz_zz[i] = -2.0 * g_xy_0_xz_zz[i] * a_exp + 4.0 * g_xy_yy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1350-1356)

    #pragma omp simd aligned(g_xy_0_yy_xx, g_xy_0_yy_xy, g_xy_0_yy_xz, g_xy_0_yy_yy, g_xy_0_yy_yz, g_xy_0_yy_zz, g_xy_yy_yy_xx, g_xy_yy_yy_xy, g_xy_yy_yy_xz, g_xy_yy_yy_yy, g_xy_yy_yy_yz, g_xy_yy_yy_zz, g_y_y_0_0_x_y_yy_xx, g_y_y_0_0_x_y_yy_xy, g_y_y_0_0_x_y_yy_xz, g_y_y_0_0_x_y_yy_yy, g_y_y_0_0_x_y_yy_yz, g_y_y_0_0_x_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_y_yy_xx[i] = -2.0 * g_xy_0_yy_xx[i] * a_exp + 4.0 * g_xy_yy_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_yy_xy[i] = -2.0 * g_xy_0_yy_xy[i] * a_exp + 4.0 * g_xy_yy_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_yy_xz[i] = -2.0 * g_xy_0_yy_xz[i] * a_exp + 4.0 * g_xy_yy_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_yy_yy[i] = -2.0 * g_xy_0_yy_yy[i] * a_exp + 4.0 * g_xy_yy_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_yy_yz[i] = -2.0 * g_xy_0_yy_yz[i] * a_exp + 4.0 * g_xy_yy_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_yy_zz[i] = -2.0 * g_xy_0_yy_zz[i] * a_exp + 4.0 * g_xy_yy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1356-1362)

    #pragma omp simd aligned(g_xy_0_yz_xx, g_xy_0_yz_xy, g_xy_0_yz_xz, g_xy_0_yz_yy, g_xy_0_yz_yz, g_xy_0_yz_zz, g_xy_yy_yz_xx, g_xy_yy_yz_xy, g_xy_yy_yz_xz, g_xy_yy_yz_yy, g_xy_yy_yz_yz, g_xy_yy_yz_zz, g_y_y_0_0_x_y_yz_xx, g_y_y_0_0_x_y_yz_xy, g_y_y_0_0_x_y_yz_xz, g_y_y_0_0_x_y_yz_yy, g_y_y_0_0_x_y_yz_yz, g_y_y_0_0_x_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_y_yz_xx[i] = -2.0 * g_xy_0_yz_xx[i] * a_exp + 4.0 * g_xy_yy_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_yz_xy[i] = -2.0 * g_xy_0_yz_xy[i] * a_exp + 4.0 * g_xy_yy_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_yz_xz[i] = -2.0 * g_xy_0_yz_xz[i] * a_exp + 4.0 * g_xy_yy_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_yz_yy[i] = -2.0 * g_xy_0_yz_yy[i] * a_exp + 4.0 * g_xy_yy_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_yz_yz[i] = -2.0 * g_xy_0_yz_yz[i] * a_exp + 4.0 * g_xy_yy_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_yz_zz[i] = -2.0 * g_xy_0_yz_zz[i] * a_exp + 4.0 * g_xy_yy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1362-1368)

    #pragma omp simd aligned(g_xy_0_zz_xx, g_xy_0_zz_xy, g_xy_0_zz_xz, g_xy_0_zz_yy, g_xy_0_zz_yz, g_xy_0_zz_zz, g_xy_yy_zz_xx, g_xy_yy_zz_xy, g_xy_yy_zz_xz, g_xy_yy_zz_yy, g_xy_yy_zz_yz, g_xy_yy_zz_zz, g_y_y_0_0_x_y_zz_xx, g_y_y_0_0_x_y_zz_xy, g_y_y_0_0_x_y_zz_xz, g_y_y_0_0_x_y_zz_yy, g_y_y_0_0_x_y_zz_yz, g_y_y_0_0_x_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_y_zz_xx[i] = -2.0 * g_xy_0_zz_xx[i] * a_exp + 4.0 * g_xy_yy_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_zz_xy[i] = -2.0 * g_xy_0_zz_xy[i] * a_exp + 4.0 * g_xy_yy_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_zz_xz[i] = -2.0 * g_xy_0_zz_xz[i] * a_exp + 4.0 * g_xy_yy_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_zz_yy[i] = -2.0 * g_xy_0_zz_yy[i] * a_exp + 4.0 * g_xy_yy_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_zz_yz[i] = -2.0 * g_xy_0_zz_yz[i] * a_exp + 4.0 * g_xy_yy_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_zz_zz[i] = -2.0 * g_xy_0_zz_zz[i] * a_exp + 4.0 * g_xy_yy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1368-1374)

    #pragma omp simd aligned(g_xy_yz_xx_xx, g_xy_yz_xx_xy, g_xy_yz_xx_xz, g_xy_yz_xx_yy, g_xy_yz_xx_yz, g_xy_yz_xx_zz, g_y_y_0_0_x_z_xx_xx, g_y_y_0_0_x_z_xx_xy, g_y_y_0_0_x_z_xx_xz, g_y_y_0_0_x_z_xx_yy, g_y_y_0_0_x_z_xx_yz, g_y_y_0_0_x_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_z_xx_xx[i] = 4.0 * g_xy_yz_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_xx_xy[i] = 4.0 * g_xy_yz_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_xx_xz[i] = 4.0 * g_xy_yz_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_xx_yy[i] = 4.0 * g_xy_yz_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_xx_yz[i] = 4.0 * g_xy_yz_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_xx_zz[i] = 4.0 * g_xy_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1374-1380)

    #pragma omp simd aligned(g_xy_yz_xy_xx, g_xy_yz_xy_xy, g_xy_yz_xy_xz, g_xy_yz_xy_yy, g_xy_yz_xy_yz, g_xy_yz_xy_zz, g_y_y_0_0_x_z_xy_xx, g_y_y_0_0_x_z_xy_xy, g_y_y_0_0_x_z_xy_xz, g_y_y_0_0_x_z_xy_yy, g_y_y_0_0_x_z_xy_yz, g_y_y_0_0_x_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_z_xy_xx[i] = 4.0 * g_xy_yz_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_xy_xy[i] = 4.0 * g_xy_yz_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_xy_xz[i] = 4.0 * g_xy_yz_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_xy_yy[i] = 4.0 * g_xy_yz_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_xy_yz[i] = 4.0 * g_xy_yz_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_xy_zz[i] = 4.0 * g_xy_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1380-1386)

    #pragma omp simd aligned(g_xy_yz_xz_xx, g_xy_yz_xz_xy, g_xy_yz_xz_xz, g_xy_yz_xz_yy, g_xy_yz_xz_yz, g_xy_yz_xz_zz, g_y_y_0_0_x_z_xz_xx, g_y_y_0_0_x_z_xz_xy, g_y_y_0_0_x_z_xz_xz, g_y_y_0_0_x_z_xz_yy, g_y_y_0_0_x_z_xz_yz, g_y_y_0_0_x_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_z_xz_xx[i] = 4.0 * g_xy_yz_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_xz_xy[i] = 4.0 * g_xy_yz_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_xz_xz[i] = 4.0 * g_xy_yz_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_xz_yy[i] = 4.0 * g_xy_yz_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_xz_yz[i] = 4.0 * g_xy_yz_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_xz_zz[i] = 4.0 * g_xy_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1386-1392)

    #pragma omp simd aligned(g_xy_yz_yy_xx, g_xy_yz_yy_xy, g_xy_yz_yy_xz, g_xy_yz_yy_yy, g_xy_yz_yy_yz, g_xy_yz_yy_zz, g_y_y_0_0_x_z_yy_xx, g_y_y_0_0_x_z_yy_xy, g_y_y_0_0_x_z_yy_xz, g_y_y_0_0_x_z_yy_yy, g_y_y_0_0_x_z_yy_yz, g_y_y_0_0_x_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_z_yy_xx[i] = 4.0 * g_xy_yz_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_yy_xy[i] = 4.0 * g_xy_yz_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_yy_xz[i] = 4.0 * g_xy_yz_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_yy_yy[i] = 4.0 * g_xy_yz_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_yy_yz[i] = 4.0 * g_xy_yz_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_yy_zz[i] = 4.0 * g_xy_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1392-1398)

    #pragma omp simd aligned(g_xy_yz_yz_xx, g_xy_yz_yz_xy, g_xy_yz_yz_xz, g_xy_yz_yz_yy, g_xy_yz_yz_yz, g_xy_yz_yz_zz, g_y_y_0_0_x_z_yz_xx, g_y_y_0_0_x_z_yz_xy, g_y_y_0_0_x_z_yz_xz, g_y_y_0_0_x_z_yz_yy, g_y_y_0_0_x_z_yz_yz, g_y_y_0_0_x_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_z_yz_xx[i] = 4.0 * g_xy_yz_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_yz_xy[i] = 4.0 * g_xy_yz_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_yz_xz[i] = 4.0 * g_xy_yz_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_yz_yy[i] = 4.0 * g_xy_yz_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_yz_yz[i] = 4.0 * g_xy_yz_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_yz_zz[i] = 4.0 * g_xy_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1398-1404)

    #pragma omp simd aligned(g_xy_yz_zz_xx, g_xy_yz_zz_xy, g_xy_yz_zz_xz, g_xy_yz_zz_yy, g_xy_yz_zz_yz, g_xy_yz_zz_zz, g_y_y_0_0_x_z_zz_xx, g_y_y_0_0_x_z_zz_xy, g_y_y_0_0_x_z_zz_xz, g_y_y_0_0_x_z_zz_yy, g_y_y_0_0_x_z_zz_yz, g_y_y_0_0_x_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_z_zz_xx[i] = 4.0 * g_xy_yz_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_zz_xy[i] = 4.0 * g_xy_yz_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_zz_xz[i] = 4.0 * g_xy_yz_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_zz_yy[i] = 4.0 * g_xy_yz_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_zz_yz[i] = 4.0 * g_xy_yz_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_zz_zz[i] = 4.0 * g_xy_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1404-1410)

    #pragma omp simd aligned(g_0_xy_xx_xx, g_0_xy_xx_xy, g_0_xy_xx_xz, g_0_xy_xx_yy, g_0_xy_xx_yz, g_0_xy_xx_zz, g_y_y_0_0_y_x_xx_xx, g_y_y_0_0_y_x_xx_xy, g_y_y_0_0_y_x_xx_xz, g_y_y_0_0_y_x_xx_yy, g_y_y_0_0_y_x_xx_yz, g_y_y_0_0_y_x_xx_zz, g_yy_xy_xx_xx, g_yy_xy_xx_xy, g_yy_xy_xx_xz, g_yy_xy_xx_yy, g_yy_xy_xx_yz, g_yy_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_x_xx_xx[i] = -2.0 * g_0_xy_xx_xx[i] * b_exp + 4.0 * g_yy_xy_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_xx_xy[i] = -2.0 * g_0_xy_xx_xy[i] * b_exp + 4.0 * g_yy_xy_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_xx_xz[i] = -2.0 * g_0_xy_xx_xz[i] * b_exp + 4.0 * g_yy_xy_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_xx_yy[i] = -2.0 * g_0_xy_xx_yy[i] * b_exp + 4.0 * g_yy_xy_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_xx_yz[i] = -2.0 * g_0_xy_xx_yz[i] * b_exp + 4.0 * g_yy_xy_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_xx_zz[i] = -2.0 * g_0_xy_xx_zz[i] * b_exp + 4.0 * g_yy_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1410-1416)

    #pragma omp simd aligned(g_0_xy_xy_xx, g_0_xy_xy_xy, g_0_xy_xy_xz, g_0_xy_xy_yy, g_0_xy_xy_yz, g_0_xy_xy_zz, g_y_y_0_0_y_x_xy_xx, g_y_y_0_0_y_x_xy_xy, g_y_y_0_0_y_x_xy_xz, g_y_y_0_0_y_x_xy_yy, g_y_y_0_0_y_x_xy_yz, g_y_y_0_0_y_x_xy_zz, g_yy_xy_xy_xx, g_yy_xy_xy_xy, g_yy_xy_xy_xz, g_yy_xy_xy_yy, g_yy_xy_xy_yz, g_yy_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_x_xy_xx[i] = -2.0 * g_0_xy_xy_xx[i] * b_exp + 4.0 * g_yy_xy_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_xy_xy[i] = -2.0 * g_0_xy_xy_xy[i] * b_exp + 4.0 * g_yy_xy_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_xy_xz[i] = -2.0 * g_0_xy_xy_xz[i] * b_exp + 4.0 * g_yy_xy_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_xy_yy[i] = -2.0 * g_0_xy_xy_yy[i] * b_exp + 4.0 * g_yy_xy_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_xy_yz[i] = -2.0 * g_0_xy_xy_yz[i] * b_exp + 4.0 * g_yy_xy_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_xy_zz[i] = -2.0 * g_0_xy_xy_zz[i] * b_exp + 4.0 * g_yy_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1416-1422)

    #pragma omp simd aligned(g_0_xy_xz_xx, g_0_xy_xz_xy, g_0_xy_xz_xz, g_0_xy_xz_yy, g_0_xy_xz_yz, g_0_xy_xz_zz, g_y_y_0_0_y_x_xz_xx, g_y_y_0_0_y_x_xz_xy, g_y_y_0_0_y_x_xz_xz, g_y_y_0_0_y_x_xz_yy, g_y_y_0_0_y_x_xz_yz, g_y_y_0_0_y_x_xz_zz, g_yy_xy_xz_xx, g_yy_xy_xz_xy, g_yy_xy_xz_xz, g_yy_xy_xz_yy, g_yy_xy_xz_yz, g_yy_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_x_xz_xx[i] = -2.0 * g_0_xy_xz_xx[i] * b_exp + 4.0 * g_yy_xy_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_xz_xy[i] = -2.0 * g_0_xy_xz_xy[i] * b_exp + 4.0 * g_yy_xy_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_xz_xz[i] = -2.0 * g_0_xy_xz_xz[i] * b_exp + 4.0 * g_yy_xy_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_xz_yy[i] = -2.0 * g_0_xy_xz_yy[i] * b_exp + 4.0 * g_yy_xy_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_xz_yz[i] = -2.0 * g_0_xy_xz_yz[i] * b_exp + 4.0 * g_yy_xy_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_xz_zz[i] = -2.0 * g_0_xy_xz_zz[i] * b_exp + 4.0 * g_yy_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1422-1428)

    #pragma omp simd aligned(g_0_xy_yy_xx, g_0_xy_yy_xy, g_0_xy_yy_xz, g_0_xy_yy_yy, g_0_xy_yy_yz, g_0_xy_yy_zz, g_y_y_0_0_y_x_yy_xx, g_y_y_0_0_y_x_yy_xy, g_y_y_0_0_y_x_yy_xz, g_y_y_0_0_y_x_yy_yy, g_y_y_0_0_y_x_yy_yz, g_y_y_0_0_y_x_yy_zz, g_yy_xy_yy_xx, g_yy_xy_yy_xy, g_yy_xy_yy_xz, g_yy_xy_yy_yy, g_yy_xy_yy_yz, g_yy_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_x_yy_xx[i] = -2.0 * g_0_xy_yy_xx[i] * b_exp + 4.0 * g_yy_xy_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_yy_xy[i] = -2.0 * g_0_xy_yy_xy[i] * b_exp + 4.0 * g_yy_xy_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_yy_xz[i] = -2.0 * g_0_xy_yy_xz[i] * b_exp + 4.0 * g_yy_xy_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_yy_yy[i] = -2.0 * g_0_xy_yy_yy[i] * b_exp + 4.0 * g_yy_xy_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_yy_yz[i] = -2.0 * g_0_xy_yy_yz[i] * b_exp + 4.0 * g_yy_xy_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_yy_zz[i] = -2.0 * g_0_xy_yy_zz[i] * b_exp + 4.0 * g_yy_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1428-1434)

    #pragma omp simd aligned(g_0_xy_yz_xx, g_0_xy_yz_xy, g_0_xy_yz_xz, g_0_xy_yz_yy, g_0_xy_yz_yz, g_0_xy_yz_zz, g_y_y_0_0_y_x_yz_xx, g_y_y_0_0_y_x_yz_xy, g_y_y_0_0_y_x_yz_xz, g_y_y_0_0_y_x_yz_yy, g_y_y_0_0_y_x_yz_yz, g_y_y_0_0_y_x_yz_zz, g_yy_xy_yz_xx, g_yy_xy_yz_xy, g_yy_xy_yz_xz, g_yy_xy_yz_yy, g_yy_xy_yz_yz, g_yy_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_x_yz_xx[i] = -2.0 * g_0_xy_yz_xx[i] * b_exp + 4.0 * g_yy_xy_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_yz_xy[i] = -2.0 * g_0_xy_yz_xy[i] * b_exp + 4.0 * g_yy_xy_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_yz_xz[i] = -2.0 * g_0_xy_yz_xz[i] * b_exp + 4.0 * g_yy_xy_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_yz_yy[i] = -2.0 * g_0_xy_yz_yy[i] * b_exp + 4.0 * g_yy_xy_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_yz_yz[i] = -2.0 * g_0_xy_yz_yz[i] * b_exp + 4.0 * g_yy_xy_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_yz_zz[i] = -2.0 * g_0_xy_yz_zz[i] * b_exp + 4.0 * g_yy_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1434-1440)

    #pragma omp simd aligned(g_0_xy_zz_xx, g_0_xy_zz_xy, g_0_xy_zz_xz, g_0_xy_zz_yy, g_0_xy_zz_yz, g_0_xy_zz_zz, g_y_y_0_0_y_x_zz_xx, g_y_y_0_0_y_x_zz_xy, g_y_y_0_0_y_x_zz_xz, g_y_y_0_0_y_x_zz_yy, g_y_y_0_0_y_x_zz_yz, g_y_y_0_0_y_x_zz_zz, g_yy_xy_zz_xx, g_yy_xy_zz_xy, g_yy_xy_zz_xz, g_yy_xy_zz_yy, g_yy_xy_zz_yz, g_yy_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_x_zz_xx[i] = -2.0 * g_0_xy_zz_xx[i] * b_exp + 4.0 * g_yy_xy_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_zz_xy[i] = -2.0 * g_0_xy_zz_xy[i] * b_exp + 4.0 * g_yy_xy_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_zz_xz[i] = -2.0 * g_0_xy_zz_xz[i] * b_exp + 4.0 * g_yy_xy_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_zz_yy[i] = -2.0 * g_0_xy_zz_yy[i] * b_exp + 4.0 * g_yy_xy_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_zz_yz[i] = -2.0 * g_0_xy_zz_yz[i] * b_exp + 4.0 * g_yy_xy_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_zz_zz[i] = -2.0 * g_0_xy_zz_zz[i] * b_exp + 4.0 * g_yy_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1440-1446)

    #pragma omp simd aligned(g_0_0_xx_xx, g_0_0_xx_xy, g_0_0_xx_xz, g_0_0_xx_yy, g_0_0_xx_yz, g_0_0_xx_zz, g_0_yy_xx_xx, g_0_yy_xx_xy, g_0_yy_xx_xz, g_0_yy_xx_yy, g_0_yy_xx_yz, g_0_yy_xx_zz, g_y_y_0_0_y_y_xx_xx, g_y_y_0_0_y_y_xx_xy, g_y_y_0_0_y_y_xx_xz, g_y_y_0_0_y_y_xx_yy, g_y_y_0_0_y_y_xx_yz, g_y_y_0_0_y_y_xx_zz, g_yy_0_xx_xx, g_yy_0_xx_xy, g_yy_0_xx_xz, g_yy_0_xx_yy, g_yy_0_xx_yz, g_yy_0_xx_zz, g_yy_yy_xx_xx, g_yy_yy_xx_xy, g_yy_yy_xx_xz, g_yy_yy_xx_yy, g_yy_yy_xx_yz, g_yy_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_y_xx_xx[i] = g_0_0_xx_xx[i] - 2.0 * g_0_yy_xx_xx[i] * b_exp - 2.0 * g_yy_0_xx_xx[i] * a_exp + 4.0 * g_yy_yy_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_xx_xy[i] = g_0_0_xx_xy[i] - 2.0 * g_0_yy_xx_xy[i] * b_exp - 2.0 * g_yy_0_xx_xy[i] * a_exp + 4.0 * g_yy_yy_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_xx_xz[i] = g_0_0_xx_xz[i] - 2.0 * g_0_yy_xx_xz[i] * b_exp - 2.0 * g_yy_0_xx_xz[i] * a_exp + 4.0 * g_yy_yy_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_xx_yy[i] = g_0_0_xx_yy[i] - 2.0 * g_0_yy_xx_yy[i] * b_exp - 2.0 * g_yy_0_xx_yy[i] * a_exp + 4.0 * g_yy_yy_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_xx_yz[i] = g_0_0_xx_yz[i] - 2.0 * g_0_yy_xx_yz[i] * b_exp - 2.0 * g_yy_0_xx_yz[i] * a_exp + 4.0 * g_yy_yy_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_xx_zz[i] = g_0_0_xx_zz[i] - 2.0 * g_0_yy_xx_zz[i] * b_exp - 2.0 * g_yy_0_xx_zz[i] * a_exp + 4.0 * g_yy_yy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1446-1452)

    #pragma omp simd aligned(g_0_0_xy_xx, g_0_0_xy_xy, g_0_0_xy_xz, g_0_0_xy_yy, g_0_0_xy_yz, g_0_0_xy_zz, g_0_yy_xy_xx, g_0_yy_xy_xy, g_0_yy_xy_xz, g_0_yy_xy_yy, g_0_yy_xy_yz, g_0_yy_xy_zz, g_y_y_0_0_y_y_xy_xx, g_y_y_0_0_y_y_xy_xy, g_y_y_0_0_y_y_xy_xz, g_y_y_0_0_y_y_xy_yy, g_y_y_0_0_y_y_xy_yz, g_y_y_0_0_y_y_xy_zz, g_yy_0_xy_xx, g_yy_0_xy_xy, g_yy_0_xy_xz, g_yy_0_xy_yy, g_yy_0_xy_yz, g_yy_0_xy_zz, g_yy_yy_xy_xx, g_yy_yy_xy_xy, g_yy_yy_xy_xz, g_yy_yy_xy_yy, g_yy_yy_xy_yz, g_yy_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_y_xy_xx[i] = g_0_0_xy_xx[i] - 2.0 * g_0_yy_xy_xx[i] * b_exp - 2.0 * g_yy_0_xy_xx[i] * a_exp + 4.0 * g_yy_yy_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_xy_xy[i] = g_0_0_xy_xy[i] - 2.0 * g_0_yy_xy_xy[i] * b_exp - 2.0 * g_yy_0_xy_xy[i] * a_exp + 4.0 * g_yy_yy_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_xy_xz[i] = g_0_0_xy_xz[i] - 2.0 * g_0_yy_xy_xz[i] * b_exp - 2.0 * g_yy_0_xy_xz[i] * a_exp + 4.0 * g_yy_yy_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_xy_yy[i] = g_0_0_xy_yy[i] - 2.0 * g_0_yy_xy_yy[i] * b_exp - 2.0 * g_yy_0_xy_yy[i] * a_exp + 4.0 * g_yy_yy_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_xy_yz[i] = g_0_0_xy_yz[i] - 2.0 * g_0_yy_xy_yz[i] * b_exp - 2.0 * g_yy_0_xy_yz[i] * a_exp + 4.0 * g_yy_yy_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_xy_zz[i] = g_0_0_xy_zz[i] - 2.0 * g_0_yy_xy_zz[i] * b_exp - 2.0 * g_yy_0_xy_zz[i] * a_exp + 4.0 * g_yy_yy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1452-1458)

    #pragma omp simd aligned(g_0_0_xz_xx, g_0_0_xz_xy, g_0_0_xz_xz, g_0_0_xz_yy, g_0_0_xz_yz, g_0_0_xz_zz, g_0_yy_xz_xx, g_0_yy_xz_xy, g_0_yy_xz_xz, g_0_yy_xz_yy, g_0_yy_xz_yz, g_0_yy_xz_zz, g_y_y_0_0_y_y_xz_xx, g_y_y_0_0_y_y_xz_xy, g_y_y_0_0_y_y_xz_xz, g_y_y_0_0_y_y_xz_yy, g_y_y_0_0_y_y_xz_yz, g_y_y_0_0_y_y_xz_zz, g_yy_0_xz_xx, g_yy_0_xz_xy, g_yy_0_xz_xz, g_yy_0_xz_yy, g_yy_0_xz_yz, g_yy_0_xz_zz, g_yy_yy_xz_xx, g_yy_yy_xz_xy, g_yy_yy_xz_xz, g_yy_yy_xz_yy, g_yy_yy_xz_yz, g_yy_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_y_xz_xx[i] = g_0_0_xz_xx[i] - 2.0 * g_0_yy_xz_xx[i] * b_exp - 2.0 * g_yy_0_xz_xx[i] * a_exp + 4.0 * g_yy_yy_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_xz_xy[i] = g_0_0_xz_xy[i] - 2.0 * g_0_yy_xz_xy[i] * b_exp - 2.0 * g_yy_0_xz_xy[i] * a_exp + 4.0 * g_yy_yy_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_xz_xz[i] = g_0_0_xz_xz[i] - 2.0 * g_0_yy_xz_xz[i] * b_exp - 2.0 * g_yy_0_xz_xz[i] * a_exp + 4.0 * g_yy_yy_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_xz_yy[i] = g_0_0_xz_yy[i] - 2.0 * g_0_yy_xz_yy[i] * b_exp - 2.0 * g_yy_0_xz_yy[i] * a_exp + 4.0 * g_yy_yy_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_xz_yz[i] = g_0_0_xz_yz[i] - 2.0 * g_0_yy_xz_yz[i] * b_exp - 2.0 * g_yy_0_xz_yz[i] * a_exp + 4.0 * g_yy_yy_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_xz_zz[i] = g_0_0_xz_zz[i] - 2.0 * g_0_yy_xz_zz[i] * b_exp - 2.0 * g_yy_0_xz_zz[i] * a_exp + 4.0 * g_yy_yy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1458-1464)

    #pragma omp simd aligned(g_0_0_yy_xx, g_0_0_yy_xy, g_0_0_yy_xz, g_0_0_yy_yy, g_0_0_yy_yz, g_0_0_yy_zz, g_0_yy_yy_xx, g_0_yy_yy_xy, g_0_yy_yy_xz, g_0_yy_yy_yy, g_0_yy_yy_yz, g_0_yy_yy_zz, g_y_y_0_0_y_y_yy_xx, g_y_y_0_0_y_y_yy_xy, g_y_y_0_0_y_y_yy_xz, g_y_y_0_0_y_y_yy_yy, g_y_y_0_0_y_y_yy_yz, g_y_y_0_0_y_y_yy_zz, g_yy_0_yy_xx, g_yy_0_yy_xy, g_yy_0_yy_xz, g_yy_0_yy_yy, g_yy_0_yy_yz, g_yy_0_yy_zz, g_yy_yy_yy_xx, g_yy_yy_yy_xy, g_yy_yy_yy_xz, g_yy_yy_yy_yy, g_yy_yy_yy_yz, g_yy_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_y_yy_xx[i] = g_0_0_yy_xx[i] - 2.0 * g_0_yy_yy_xx[i] * b_exp - 2.0 * g_yy_0_yy_xx[i] * a_exp + 4.0 * g_yy_yy_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_yy_xy[i] = g_0_0_yy_xy[i] - 2.0 * g_0_yy_yy_xy[i] * b_exp - 2.0 * g_yy_0_yy_xy[i] * a_exp + 4.0 * g_yy_yy_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_yy_xz[i] = g_0_0_yy_xz[i] - 2.0 * g_0_yy_yy_xz[i] * b_exp - 2.0 * g_yy_0_yy_xz[i] * a_exp + 4.0 * g_yy_yy_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_yy_yy[i] = g_0_0_yy_yy[i] - 2.0 * g_0_yy_yy_yy[i] * b_exp - 2.0 * g_yy_0_yy_yy[i] * a_exp + 4.0 * g_yy_yy_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_yy_yz[i] = g_0_0_yy_yz[i] - 2.0 * g_0_yy_yy_yz[i] * b_exp - 2.0 * g_yy_0_yy_yz[i] * a_exp + 4.0 * g_yy_yy_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_yy_zz[i] = g_0_0_yy_zz[i] - 2.0 * g_0_yy_yy_zz[i] * b_exp - 2.0 * g_yy_0_yy_zz[i] * a_exp + 4.0 * g_yy_yy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1464-1470)

    #pragma omp simd aligned(g_0_0_yz_xx, g_0_0_yz_xy, g_0_0_yz_xz, g_0_0_yz_yy, g_0_0_yz_yz, g_0_0_yz_zz, g_0_yy_yz_xx, g_0_yy_yz_xy, g_0_yy_yz_xz, g_0_yy_yz_yy, g_0_yy_yz_yz, g_0_yy_yz_zz, g_y_y_0_0_y_y_yz_xx, g_y_y_0_0_y_y_yz_xy, g_y_y_0_0_y_y_yz_xz, g_y_y_0_0_y_y_yz_yy, g_y_y_0_0_y_y_yz_yz, g_y_y_0_0_y_y_yz_zz, g_yy_0_yz_xx, g_yy_0_yz_xy, g_yy_0_yz_xz, g_yy_0_yz_yy, g_yy_0_yz_yz, g_yy_0_yz_zz, g_yy_yy_yz_xx, g_yy_yy_yz_xy, g_yy_yy_yz_xz, g_yy_yy_yz_yy, g_yy_yy_yz_yz, g_yy_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_y_yz_xx[i] = g_0_0_yz_xx[i] - 2.0 * g_0_yy_yz_xx[i] * b_exp - 2.0 * g_yy_0_yz_xx[i] * a_exp + 4.0 * g_yy_yy_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_yz_xy[i] = g_0_0_yz_xy[i] - 2.0 * g_0_yy_yz_xy[i] * b_exp - 2.0 * g_yy_0_yz_xy[i] * a_exp + 4.0 * g_yy_yy_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_yz_xz[i] = g_0_0_yz_xz[i] - 2.0 * g_0_yy_yz_xz[i] * b_exp - 2.0 * g_yy_0_yz_xz[i] * a_exp + 4.0 * g_yy_yy_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_yz_yy[i] = g_0_0_yz_yy[i] - 2.0 * g_0_yy_yz_yy[i] * b_exp - 2.0 * g_yy_0_yz_yy[i] * a_exp + 4.0 * g_yy_yy_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_yz_yz[i] = g_0_0_yz_yz[i] - 2.0 * g_0_yy_yz_yz[i] * b_exp - 2.0 * g_yy_0_yz_yz[i] * a_exp + 4.0 * g_yy_yy_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_yz_zz[i] = g_0_0_yz_zz[i] - 2.0 * g_0_yy_yz_zz[i] * b_exp - 2.0 * g_yy_0_yz_zz[i] * a_exp + 4.0 * g_yy_yy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1470-1476)

    #pragma omp simd aligned(g_0_0_zz_xx, g_0_0_zz_xy, g_0_0_zz_xz, g_0_0_zz_yy, g_0_0_zz_yz, g_0_0_zz_zz, g_0_yy_zz_xx, g_0_yy_zz_xy, g_0_yy_zz_xz, g_0_yy_zz_yy, g_0_yy_zz_yz, g_0_yy_zz_zz, g_y_y_0_0_y_y_zz_xx, g_y_y_0_0_y_y_zz_xy, g_y_y_0_0_y_y_zz_xz, g_y_y_0_0_y_y_zz_yy, g_y_y_0_0_y_y_zz_yz, g_y_y_0_0_y_y_zz_zz, g_yy_0_zz_xx, g_yy_0_zz_xy, g_yy_0_zz_xz, g_yy_0_zz_yy, g_yy_0_zz_yz, g_yy_0_zz_zz, g_yy_yy_zz_xx, g_yy_yy_zz_xy, g_yy_yy_zz_xz, g_yy_yy_zz_yy, g_yy_yy_zz_yz, g_yy_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_y_zz_xx[i] = g_0_0_zz_xx[i] - 2.0 * g_0_yy_zz_xx[i] * b_exp - 2.0 * g_yy_0_zz_xx[i] * a_exp + 4.0 * g_yy_yy_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_zz_xy[i] = g_0_0_zz_xy[i] - 2.0 * g_0_yy_zz_xy[i] * b_exp - 2.0 * g_yy_0_zz_xy[i] * a_exp + 4.0 * g_yy_yy_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_zz_xz[i] = g_0_0_zz_xz[i] - 2.0 * g_0_yy_zz_xz[i] * b_exp - 2.0 * g_yy_0_zz_xz[i] * a_exp + 4.0 * g_yy_yy_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_zz_yy[i] = g_0_0_zz_yy[i] - 2.0 * g_0_yy_zz_yy[i] * b_exp - 2.0 * g_yy_0_zz_yy[i] * a_exp + 4.0 * g_yy_yy_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_zz_yz[i] = g_0_0_zz_yz[i] - 2.0 * g_0_yy_zz_yz[i] * b_exp - 2.0 * g_yy_0_zz_yz[i] * a_exp + 4.0 * g_yy_yy_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_zz_zz[i] = g_0_0_zz_zz[i] - 2.0 * g_0_yy_zz_zz[i] * b_exp - 2.0 * g_yy_0_zz_zz[i] * a_exp + 4.0 * g_yy_yy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1476-1482)

    #pragma omp simd aligned(g_0_yz_xx_xx, g_0_yz_xx_xy, g_0_yz_xx_xz, g_0_yz_xx_yy, g_0_yz_xx_yz, g_0_yz_xx_zz, g_y_y_0_0_y_z_xx_xx, g_y_y_0_0_y_z_xx_xy, g_y_y_0_0_y_z_xx_xz, g_y_y_0_0_y_z_xx_yy, g_y_y_0_0_y_z_xx_yz, g_y_y_0_0_y_z_xx_zz, g_yy_yz_xx_xx, g_yy_yz_xx_xy, g_yy_yz_xx_xz, g_yy_yz_xx_yy, g_yy_yz_xx_yz, g_yy_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_z_xx_xx[i] = -2.0 * g_0_yz_xx_xx[i] * b_exp + 4.0 * g_yy_yz_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_xx_xy[i] = -2.0 * g_0_yz_xx_xy[i] * b_exp + 4.0 * g_yy_yz_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_xx_xz[i] = -2.0 * g_0_yz_xx_xz[i] * b_exp + 4.0 * g_yy_yz_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_xx_yy[i] = -2.0 * g_0_yz_xx_yy[i] * b_exp + 4.0 * g_yy_yz_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_xx_yz[i] = -2.0 * g_0_yz_xx_yz[i] * b_exp + 4.0 * g_yy_yz_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_xx_zz[i] = -2.0 * g_0_yz_xx_zz[i] * b_exp + 4.0 * g_yy_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1482-1488)

    #pragma omp simd aligned(g_0_yz_xy_xx, g_0_yz_xy_xy, g_0_yz_xy_xz, g_0_yz_xy_yy, g_0_yz_xy_yz, g_0_yz_xy_zz, g_y_y_0_0_y_z_xy_xx, g_y_y_0_0_y_z_xy_xy, g_y_y_0_0_y_z_xy_xz, g_y_y_0_0_y_z_xy_yy, g_y_y_0_0_y_z_xy_yz, g_y_y_0_0_y_z_xy_zz, g_yy_yz_xy_xx, g_yy_yz_xy_xy, g_yy_yz_xy_xz, g_yy_yz_xy_yy, g_yy_yz_xy_yz, g_yy_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_z_xy_xx[i] = -2.0 * g_0_yz_xy_xx[i] * b_exp + 4.0 * g_yy_yz_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_xy_xy[i] = -2.0 * g_0_yz_xy_xy[i] * b_exp + 4.0 * g_yy_yz_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_xy_xz[i] = -2.0 * g_0_yz_xy_xz[i] * b_exp + 4.0 * g_yy_yz_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_xy_yy[i] = -2.0 * g_0_yz_xy_yy[i] * b_exp + 4.0 * g_yy_yz_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_xy_yz[i] = -2.0 * g_0_yz_xy_yz[i] * b_exp + 4.0 * g_yy_yz_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_xy_zz[i] = -2.0 * g_0_yz_xy_zz[i] * b_exp + 4.0 * g_yy_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1488-1494)

    #pragma omp simd aligned(g_0_yz_xz_xx, g_0_yz_xz_xy, g_0_yz_xz_xz, g_0_yz_xz_yy, g_0_yz_xz_yz, g_0_yz_xz_zz, g_y_y_0_0_y_z_xz_xx, g_y_y_0_0_y_z_xz_xy, g_y_y_0_0_y_z_xz_xz, g_y_y_0_0_y_z_xz_yy, g_y_y_0_0_y_z_xz_yz, g_y_y_0_0_y_z_xz_zz, g_yy_yz_xz_xx, g_yy_yz_xz_xy, g_yy_yz_xz_xz, g_yy_yz_xz_yy, g_yy_yz_xz_yz, g_yy_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_z_xz_xx[i] = -2.0 * g_0_yz_xz_xx[i] * b_exp + 4.0 * g_yy_yz_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_xz_xy[i] = -2.0 * g_0_yz_xz_xy[i] * b_exp + 4.0 * g_yy_yz_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_xz_xz[i] = -2.0 * g_0_yz_xz_xz[i] * b_exp + 4.0 * g_yy_yz_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_xz_yy[i] = -2.0 * g_0_yz_xz_yy[i] * b_exp + 4.0 * g_yy_yz_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_xz_yz[i] = -2.0 * g_0_yz_xz_yz[i] * b_exp + 4.0 * g_yy_yz_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_xz_zz[i] = -2.0 * g_0_yz_xz_zz[i] * b_exp + 4.0 * g_yy_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1494-1500)

    #pragma omp simd aligned(g_0_yz_yy_xx, g_0_yz_yy_xy, g_0_yz_yy_xz, g_0_yz_yy_yy, g_0_yz_yy_yz, g_0_yz_yy_zz, g_y_y_0_0_y_z_yy_xx, g_y_y_0_0_y_z_yy_xy, g_y_y_0_0_y_z_yy_xz, g_y_y_0_0_y_z_yy_yy, g_y_y_0_0_y_z_yy_yz, g_y_y_0_0_y_z_yy_zz, g_yy_yz_yy_xx, g_yy_yz_yy_xy, g_yy_yz_yy_xz, g_yy_yz_yy_yy, g_yy_yz_yy_yz, g_yy_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_z_yy_xx[i] = -2.0 * g_0_yz_yy_xx[i] * b_exp + 4.0 * g_yy_yz_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_yy_xy[i] = -2.0 * g_0_yz_yy_xy[i] * b_exp + 4.0 * g_yy_yz_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_yy_xz[i] = -2.0 * g_0_yz_yy_xz[i] * b_exp + 4.0 * g_yy_yz_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_yy_yy[i] = -2.0 * g_0_yz_yy_yy[i] * b_exp + 4.0 * g_yy_yz_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_yy_yz[i] = -2.0 * g_0_yz_yy_yz[i] * b_exp + 4.0 * g_yy_yz_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_yy_zz[i] = -2.0 * g_0_yz_yy_zz[i] * b_exp + 4.0 * g_yy_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1500-1506)

    #pragma omp simd aligned(g_0_yz_yz_xx, g_0_yz_yz_xy, g_0_yz_yz_xz, g_0_yz_yz_yy, g_0_yz_yz_yz, g_0_yz_yz_zz, g_y_y_0_0_y_z_yz_xx, g_y_y_0_0_y_z_yz_xy, g_y_y_0_0_y_z_yz_xz, g_y_y_0_0_y_z_yz_yy, g_y_y_0_0_y_z_yz_yz, g_y_y_0_0_y_z_yz_zz, g_yy_yz_yz_xx, g_yy_yz_yz_xy, g_yy_yz_yz_xz, g_yy_yz_yz_yy, g_yy_yz_yz_yz, g_yy_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_z_yz_xx[i] = -2.0 * g_0_yz_yz_xx[i] * b_exp + 4.0 * g_yy_yz_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_yz_xy[i] = -2.0 * g_0_yz_yz_xy[i] * b_exp + 4.0 * g_yy_yz_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_yz_xz[i] = -2.0 * g_0_yz_yz_xz[i] * b_exp + 4.0 * g_yy_yz_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_yz_yy[i] = -2.0 * g_0_yz_yz_yy[i] * b_exp + 4.0 * g_yy_yz_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_yz_yz[i] = -2.0 * g_0_yz_yz_yz[i] * b_exp + 4.0 * g_yy_yz_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_yz_zz[i] = -2.0 * g_0_yz_yz_zz[i] * b_exp + 4.0 * g_yy_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1506-1512)

    #pragma omp simd aligned(g_0_yz_zz_xx, g_0_yz_zz_xy, g_0_yz_zz_xz, g_0_yz_zz_yy, g_0_yz_zz_yz, g_0_yz_zz_zz, g_y_y_0_0_y_z_zz_xx, g_y_y_0_0_y_z_zz_xy, g_y_y_0_0_y_z_zz_xz, g_y_y_0_0_y_z_zz_yy, g_y_y_0_0_y_z_zz_yz, g_y_y_0_0_y_z_zz_zz, g_yy_yz_zz_xx, g_yy_yz_zz_xy, g_yy_yz_zz_xz, g_yy_yz_zz_yy, g_yy_yz_zz_yz, g_yy_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_z_zz_xx[i] = -2.0 * g_0_yz_zz_xx[i] * b_exp + 4.0 * g_yy_yz_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_zz_xy[i] = -2.0 * g_0_yz_zz_xy[i] * b_exp + 4.0 * g_yy_yz_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_zz_xz[i] = -2.0 * g_0_yz_zz_xz[i] * b_exp + 4.0 * g_yy_yz_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_zz_yy[i] = -2.0 * g_0_yz_zz_yy[i] * b_exp + 4.0 * g_yy_yz_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_zz_yz[i] = -2.0 * g_0_yz_zz_yz[i] * b_exp + 4.0 * g_yy_yz_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_zz_zz[i] = -2.0 * g_0_yz_zz_zz[i] * b_exp + 4.0 * g_yy_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1512-1518)

    #pragma omp simd aligned(g_y_y_0_0_z_x_xx_xx, g_y_y_0_0_z_x_xx_xy, g_y_y_0_0_z_x_xx_xz, g_y_y_0_0_z_x_xx_yy, g_y_y_0_0_z_x_xx_yz, g_y_y_0_0_z_x_xx_zz, g_yz_xy_xx_xx, g_yz_xy_xx_xy, g_yz_xy_xx_xz, g_yz_xy_xx_yy, g_yz_xy_xx_yz, g_yz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_x_xx_xx[i] = 4.0 * g_yz_xy_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_xx_xy[i] = 4.0 * g_yz_xy_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_xx_xz[i] = 4.0 * g_yz_xy_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_xx_yy[i] = 4.0 * g_yz_xy_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_xx_yz[i] = 4.0 * g_yz_xy_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_xx_zz[i] = 4.0 * g_yz_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1518-1524)

    #pragma omp simd aligned(g_y_y_0_0_z_x_xy_xx, g_y_y_0_0_z_x_xy_xy, g_y_y_0_0_z_x_xy_xz, g_y_y_0_0_z_x_xy_yy, g_y_y_0_0_z_x_xy_yz, g_y_y_0_0_z_x_xy_zz, g_yz_xy_xy_xx, g_yz_xy_xy_xy, g_yz_xy_xy_xz, g_yz_xy_xy_yy, g_yz_xy_xy_yz, g_yz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_x_xy_xx[i] = 4.0 * g_yz_xy_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_xy_xy[i] = 4.0 * g_yz_xy_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_xy_xz[i] = 4.0 * g_yz_xy_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_xy_yy[i] = 4.0 * g_yz_xy_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_xy_yz[i] = 4.0 * g_yz_xy_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_xy_zz[i] = 4.0 * g_yz_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1524-1530)

    #pragma omp simd aligned(g_y_y_0_0_z_x_xz_xx, g_y_y_0_0_z_x_xz_xy, g_y_y_0_0_z_x_xz_xz, g_y_y_0_0_z_x_xz_yy, g_y_y_0_0_z_x_xz_yz, g_y_y_0_0_z_x_xz_zz, g_yz_xy_xz_xx, g_yz_xy_xz_xy, g_yz_xy_xz_xz, g_yz_xy_xz_yy, g_yz_xy_xz_yz, g_yz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_x_xz_xx[i] = 4.0 * g_yz_xy_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_xz_xy[i] = 4.0 * g_yz_xy_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_xz_xz[i] = 4.0 * g_yz_xy_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_xz_yy[i] = 4.0 * g_yz_xy_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_xz_yz[i] = 4.0 * g_yz_xy_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_xz_zz[i] = 4.0 * g_yz_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1530-1536)

    #pragma omp simd aligned(g_y_y_0_0_z_x_yy_xx, g_y_y_0_0_z_x_yy_xy, g_y_y_0_0_z_x_yy_xz, g_y_y_0_0_z_x_yy_yy, g_y_y_0_0_z_x_yy_yz, g_y_y_0_0_z_x_yy_zz, g_yz_xy_yy_xx, g_yz_xy_yy_xy, g_yz_xy_yy_xz, g_yz_xy_yy_yy, g_yz_xy_yy_yz, g_yz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_x_yy_xx[i] = 4.0 * g_yz_xy_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_yy_xy[i] = 4.0 * g_yz_xy_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_yy_xz[i] = 4.0 * g_yz_xy_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_yy_yy[i] = 4.0 * g_yz_xy_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_yy_yz[i] = 4.0 * g_yz_xy_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_yy_zz[i] = 4.0 * g_yz_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1536-1542)

    #pragma omp simd aligned(g_y_y_0_0_z_x_yz_xx, g_y_y_0_0_z_x_yz_xy, g_y_y_0_0_z_x_yz_xz, g_y_y_0_0_z_x_yz_yy, g_y_y_0_0_z_x_yz_yz, g_y_y_0_0_z_x_yz_zz, g_yz_xy_yz_xx, g_yz_xy_yz_xy, g_yz_xy_yz_xz, g_yz_xy_yz_yy, g_yz_xy_yz_yz, g_yz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_x_yz_xx[i] = 4.0 * g_yz_xy_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_yz_xy[i] = 4.0 * g_yz_xy_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_yz_xz[i] = 4.0 * g_yz_xy_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_yz_yy[i] = 4.0 * g_yz_xy_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_yz_yz[i] = 4.0 * g_yz_xy_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_yz_zz[i] = 4.0 * g_yz_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1542-1548)

    #pragma omp simd aligned(g_y_y_0_0_z_x_zz_xx, g_y_y_0_0_z_x_zz_xy, g_y_y_0_0_z_x_zz_xz, g_y_y_0_0_z_x_zz_yy, g_y_y_0_0_z_x_zz_yz, g_y_y_0_0_z_x_zz_zz, g_yz_xy_zz_xx, g_yz_xy_zz_xy, g_yz_xy_zz_xz, g_yz_xy_zz_yy, g_yz_xy_zz_yz, g_yz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_x_zz_xx[i] = 4.0 * g_yz_xy_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_zz_xy[i] = 4.0 * g_yz_xy_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_zz_xz[i] = 4.0 * g_yz_xy_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_zz_yy[i] = 4.0 * g_yz_xy_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_zz_yz[i] = 4.0 * g_yz_xy_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_zz_zz[i] = 4.0 * g_yz_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1548-1554)

    #pragma omp simd aligned(g_y_y_0_0_z_y_xx_xx, g_y_y_0_0_z_y_xx_xy, g_y_y_0_0_z_y_xx_xz, g_y_y_0_0_z_y_xx_yy, g_y_y_0_0_z_y_xx_yz, g_y_y_0_0_z_y_xx_zz, g_yz_0_xx_xx, g_yz_0_xx_xy, g_yz_0_xx_xz, g_yz_0_xx_yy, g_yz_0_xx_yz, g_yz_0_xx_zz, g_yz_yy_xx_xx, g_yz_yy_xx_xy, g_yz_yy_xx_xz, g_yz_yy_xx_yy, g_yz_yy_xx_yz, g_yz_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_y_xx_xx[i] = -2.0 * g_yz_0_xx_xx[i] * a_exp + 4.0 * g_yz_yy_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_xx_xy[i] = -2.0 * g_yz_0_xx_xy[i] * a_exp + 4.0 * g_yz_yy_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_xx_xz[i] = -2.0 * g_yz_0_xx_xz[i] * a_exp + 4.0 * g_yz_yy_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_xx_yy[i] = -2.0 * g_yz_0_xx_yy[i] * a_exp + 4.0 * g_yz_yy_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_xx_yz[i] = -2.0 * g_yz_0_xx_yz[i] * a_exp + 4.0 * g_yz_yy_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_xx_zz[i] = -2.0 * g_yz_0_xx_zz[i] * a_exp + 4.0 * g_yz_yy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1554-1560)

    #pragma omp simd aligned(g_y_y_0_0_z_y_xy_xx, g_y_y_0_0_z_y_xy_xy, g_y_y_0_0_z_y_xy_xz, g_y_y_0_0_z_y_xy_yy, g_y_y_0_0_z_y_xy_yz, g_y_y_0_0_z_y_xy_zz, g_yz_0_xy_xx, g_yz_0_xy_xy, g_yz_0_xy_xz, g_yz_0_xy_yy, g_yz_0_xy_yz, g_yz_0_xy_zz, g_yz_yy_xy_xx, g_yz_yy_xy_xy, g_yz_yy_xy_xz, g_yz_yy_xy_yy, g_yz_yy_xy_yz, g_yz_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_y_xy_xx[i] = -2.0 * g_yz_0_xy_xx[i] * a_exp + 4.0 * g_yz_yy_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_xy_xy[i] = -2.0 * g_yz_0_xy_xy[i] * a_exp + 4.0 * g_yz_yy_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_xy_xz[i] = -2.0 * g_yz_0_xy_xz[i] * a_exp + 4.0 * g_yz_yy_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_xy_yy[i] = -2.0 * g_yz_0_xy_yy[i] * a_exp + 4.0 * g_yz_yy_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_xy_yz[i] = -2.0 * g_yz_0_xy_yz[i] * a_exp + 4.0 * g_yz_yy_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_xy_zz[i] = -2.0 * g_yz_0_xy_zz[i] * a_exp + 4.0 * g_yz_yy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1560-1566)

    #pragma omp simd aligned(g_y_y_0_0_z_y_xz_xx, g_y_y_0_0_z_y_xz_xy, g_y_y_0_0_z_y_xz_xz, g_y_y_0_0_z_y_xz_yy, g_y_y_0_0_z_y_xz_yz, g_y_y_0_0_z_y_xz_zz, g_yz_0_xz_xx, g_yz_0_xz_xy, g_yz_0_xz_xz, g_yz_0_xz_yy, g_yz_0_xz_yz, g_yz_0_xz_zz, g_yz_yy_xz_xx, g_yz_yy_xz_xy, g_yz_yy_xz_xz, g_yz_yy_xz_yy, g_yz_yy_xz_yz, g_yz_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_y_xz_xx[i] = -2.0 * g_yz_0_xz_xx[i] * a_exp + 4.0 * g_yz_yy_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_xz_xy[i] = -2.0 * g_yz_0_xz_xy[i] * a_exp + 4.0 * g_yz_yy_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_xz_xz[i] = -2.0 * g_yz_0_xz_xz[i] * a_exp + 4.0 * g_yz_yy_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_xz_yy[i] = -2.0 * g_yz_0_xz_yy[i] * a_exp + 4.0 * g_yz_yy_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_xz_yz[i] = -2.0 * g_yz_0_xz_yz[i] * a_exp + 4.0 * g_yz_yy_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_xz_zz[i] = -2.0 * g_yz_0_xz_zz[i] * a_exp + 4.0 * g_yz_yy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1566-1572)

    #pragma omp simd aligned(g_y_y_0_0_z_y_yy_xx, g_y_y_0_0_z_y_yy_xy, g_y_y_0_0_z_y_yy_xz, g_y_y_0_0_z_y_yy_yy, g_y_y_0_0_z_y_yy_yz, g_y_y_0_0_z_y_yy_zz, g_yz_0_yy_xx, g_yz_0_yy_xy, g_yz_0_yy_xz, g_yz_0_yy_yy, g_yz_0_yy_yz, g_yz_0_yy_zz, g_yz_yy_yy_xx, g_yz_yy_yy_xy, g_yz_yy_yy_xz, g_yz_yy_yy_yy, g_yz_yy_yy_yz, g_yz_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_y_yy_xx[i] = -2.0 * g_yz_0_yy_xx[i] * a_exp + 4.0 * g_yz_yy_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_yy_xy[i] = -2.0 * g_yz_0_yy_xy[i] * a_exp + 4.0 * g_yz_yy_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_yy_xz[i] = -2.0 * g_yz_0_yy_xz[i] * a_exp + 4.0 * g_yz_yy_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_yy_yy[i] = -2.0 * g_yz_0_yy_yy[i] * a_exp + 4.0 * g_yz_yy_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_yy_yz[i] = -2.0 * g_yz_0_yy_yz[i] * a_exp + 4.0 * g_yz_yy_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_yy_zz[i] = -2.0 * g_yz_0_yy_zz[i] * a_exp + 4.0 * g_yz_yy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1572-1578)

    #pragma omp simd aligned(g_y_y_0_0_z_y_yz_xx, g_y_y_0_0_z_y_yz_xy, g_y_y_0_0_z_y_yz_xz, g_y_y_0_0_z_y_yz_yy, g_y_y_0_0_z_y_yz_yz, g_y_y_0_0_z_y_yz_zz, g_yz_0_yz_xx, g_yz_0_yz_xy, g_yz_0_yz_xz, g_yz_0_yz_yy, g_yz_0_yz_yz, g_yz_0_yz_zz, g_yz_yy_yz_xx, g_yz_yy_yz_xy, g_yz_yy_yz_xz, g_yz_yy_yz_yy, g_yz_yy_yz_yz, g_yz_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_y_yz_xx[i] = -2.0 * g_yz_0_yz_xx[i] * a_exp + 4.0 * g_yz_yy_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_yz_xy[i] = -2.0 * g_yz_0_yz_xy[i] * a_exp + 4.0 * g_yz_yy_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_yz_xz[i] = -2.0 * g_yz_0_yz_xz[i] * a_exp + 4.0 * g_yz_yy_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_yz_yy[i] = -2.0 * g_yz_0_yz_yy[i] * a_exp + 4.0 * g_yz_yy_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_yz_yz[i] = -2.0 * g_yz_0_yz_yz[i] * a_exp + 4.0 * g_yz_yy_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_yz_zz[i] = -2.0 * g_yz_0_yz_zz[i] * a_exp + 4.0 * g_yz_yy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1578-1584)

    #pragma omp simd aligned(g_y_y_0_0_z_y_zz_xx, g_y_y_0_0_z_y_zz_xy, g_y_y_0_0_z_y_zz_xz, g_y_y_0_0_z_y_zz_yy, g_y_y_0_0_z_y_zz_yz, g_y_y_0_0_z_y_zz_zz, g_yz_0_zz_xx, g_yz_0_zz_xy, g_yz_0_zz_xz, g_yz_0_zz_yy, g_yz_0_zz_yz, g_yz_0_zz_zz, g_yz_yy_zz_xx, g_yz_yy_zz_xy, g_yz_yy_zz_xz, g_yz_yy_zz_yy, g_yz_yy_zz_yz, g_yz_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_y_zz_xx[i] = -2.0 * g_yz_0_zz_xx[i] * a_exp + 4.0 * g_yz_yy_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_zz_xy[i] = -2.0 * g_yz_0_zz_xy[i] * a_exp + 4.0 * g_yz_yy_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_zz_xz[i] = -2.0 * g_yz_0_zz_xz[i] * a_exp + 4.0 * g_yz_yy_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_zz_yy[i] = -2.0 * g_yz_0_zz_yy[i] * a_exp + 4.0 * g_yz_yy_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_zz_yz[i] = -2.0 * g_yz_0_zz_yz[i] * a_exp + 4.0 * g_yz_yy_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_zz_zz[i] = -2.0 * g_yz_0_zz_zz[i] * a_exp + 4.0 * g_yz_yy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1584-1590)

    #pragma omp simd aligned(g_y_y_0_0_z_z_xx_xx, g_y_y_0_0_z_z_xx_xy, g_y_y_0_0_z_z_xx_xz, g_y_y_0_0_z_z_xx_yy, g_y_y_0_0_z_z_xx_yz, g_y_y_0_0_z_z_xx_zz, g_yz_yz_xx_xx, g_yz_yz_xx_xy, g_yz_yz_xx_xz, g_yz_yz_xx_yy, g_yz_yz_xx_yz, g_yz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_z_xx_xx[i] = 4.0 * g_yz_yz_xx_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_xx_xy[i] = 4.0 * g_yz_yz_xx_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_xx_xz[i] = 4.0 * g_yz_yz_xx_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_xx_yy[i] = 4.0 * g_yz_yz_xx_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_xx_yz[i] = 4.0 * g_yz_yz_xx_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_xx_zz[i] = 4.0 * g_yz_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1590-1596)

    #pragma omp simd aligned(g_y_y_0_0_z_z_xy_xx, g_y_y_0_0_z_z_xy_xy, g_y_y_0_0_z_z_xy_xz, g_y_y_0_0_z_z_xy_yy, g_y_y_0_0_z_z_xy_yz, g_y_y_0_0_z_z_xy_zz, g_yz_yz_xy_xx, g_yz_yz_xy_xy, g_yz_yz_xy_xz, g_yz_yz_xy_yy, g_yz_yz_xy_yz, g_yz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_z_xy_xx[i] = 4.0 * g_yz_yz_xy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_xy_xy[i] = 4.0 * g_yz_yz_xy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_xy_xz[i] = 4.0 * g_yz_yz_xy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_xy_yy[i] = 4.0 * g_yz_yz_xy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_xy_yz[i] = 4.0 * g_yz_yz_xy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_xy_zz[i] = 4.0 * g_yz_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1596-1602)

    #pragma omp simd aligned(g_y_y_0_0_z_z_xz_xx, g_y_y_0_0_z_z_xz_xy, g_y_y_0_0_z_z_xz_xz, g_y_y_0_0_z_z_xz_yy, g_y_y_0_0_z_z_xz_yz, g_y_y_0_0_z_z_xz_zz, g_yz_yz_xz_xx, g_yz_yz_xz_xy, g_yz_yz_xz_xz, g_yz_yz_xz_yy, g_yz_yz_xz_yz, g_yz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_z_xz_xx[i] = 4.0 * g_yz_yz_xz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_xz_xy[i] = 4.0 * g_yz_yz_xz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_xz_xz[i] = 4.0 * g_yz_yz_xz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_xz_yy[i] = 4.0 * g_yz_yz_xz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_xz_yz[i] = 4.0 * g_yz_yz_xz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_xz_zz[i] = 4.0 * g_yz_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1602-1608)

    #pragma omp simd aligned(g_y_y_0_0_z_z_yy_xx, g_y_y_0_0_z_z_yy_xy, g_y_y_0_0_z_z_yy_xz, g_y_y_0_0_z_z_yy_yy, g_y_y_0_0_z_z_yy_yz, g_y_y_0_0_z_z_yy_zz, g_yz_yz_yy_xx, g_yz_yz_yy_xy, g_yz_yz_yy_xz, g_yz_yz_yy_yy, g_yz_yz_yy_yz, g_yz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_z_yy_xx[i] = 4.0 * g_yz_yz_yy_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_yy_xy[i] = 4.0 * g_yz_yz_yy_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_yy_xz[i] = 4.0 * g_yz_yz_yy_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_yy_yy[i] = 4.0 * g_yz_yz_yy_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_yy_yz[i] = 4.0 * g_yz_yz_yy_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_yy_zz[i] = 4.0 * g_yz_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1608-1614)

    #pragma omp simd aligned(g_y_y_0_0_z_z_yz_xx, g_y_y_0_0_z_z_yz_xy, g_y_y_0_0_z_z_yz_xz, g_y_y_0_0_z_z_yz_yy, g_y_y_0_0_z_z_yz_yz, g_y_y_0_0_z_z_yz_zz, g_yz_yz_yz_xx, g_yz_yz_yz_xy, g_yz_yz_yz_xz, g_yz_yz_yz_yy, g_yz_yz_yz_yz, g_yz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_z_yz_xx[i] = 4.0 * g_yz_yz_yz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_yz_xy[i] = 4.0 * g_yz_yz_yz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_yz_xz[i] = 4.0 * g_yz_yz_yz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_yz_yy[i] = 4.0 * g_yz_yz_yz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_yz_yz[i] = 4.0 * g_yz_yz_yz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_yz_zz[i] = 4.0 * g_yz_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1614-1620)

    #pragma omp simd aligned(g_y_y_0_0_z_z_zz_xx, g_y_y_0_0_z_z_zz_xy, g_y_y_0_0_z_z_zz_xz, g_y_y_0_0_z_z_zz_yy, g_y_y_0_0_z_z_zz_yz, g_y_y_0_0_z_z_zz_zz, g_yz_yz_zz_xx, g_yz_yz_zz_xy, g_yz_yz_zz_xz, g_yz_yz_zz_yy, g_yz_yz_zz_yz, g_yz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_z_zz_xx[i] = 4.0 * g_yz_yz_zz_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_zz_xy[i] = 4.0 * g_yz_yz_zz_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_zz_xz[i] = 4.0 * g_yz_yz_zz_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_zz_yy[i] = 4.0 * g_yz_yz_zz_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_zz_yz[i] = 4.0 * g_yz_yz_zz_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_zz_zz[i] = 4.0 * g_yz_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1620-1626)

    #pragma omp simd aligned(g_xy_xz_xx_xx, g_xy_xz_xx_xy, g_xy_xz_xx_xz, g_xy_xz_xx_yy, g_xy_xz_xx_yz, g_xy_xz_xx_zz, g_y_z_0_0_x_x_xx_xx, g_y_z_0_0_x_x_xx_xy, g_y_z_0_0_x_x_xx_xz, g_y_z_0_0_x_x_xx_yy, g_y_z_0_0_x_x_xx_yz, g_y_z_0_0_x_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_x_xx_xx[i] = 4.0 * g_xy_xz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_xx_xy[i] = 4.0 * g_xy_xz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_xx_xz[i] = 4.0 * g_xy_xz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_xx_yy[i] = 4.0 * g_xy_xz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_xx_yz[i] = 4.0 * g_xy_xz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_xx_zz[i] = 4.0 * g_xy_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1626-1632)

    #pragma omp simd aligned(g_xy_xz_xy_xx, g_xy_xz_xy_xy, g_xy_xz_xy_xz, g_xy_xz_xy_yy, g_xy_xz_xy_yz, g_xy_xz_xy_zz, g_y_z_0_0_x_x_xy_xx, g_y_z_0_0_x_x_xy_xy, g_y_z_0_0_x_x_xy_xz, g_y_z_0_0_x_x_xy_yy, g_y_z_0_0_x_x_xy_yz, g_y_z_0_0_x_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_x_xy_xx[i] = 4.0 * g_xy_xz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_xy_xy[i] = 4.0 * g_xy_xz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_xy_xz[i] = 4.0 * g_xy_xz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_xy_yy[i] = 4.0 * g_xy_xz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_xy_yz[i] = 4.0 * g_xy_xz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_xy_zz[i] = 4.0 * g_xy_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1632-1638)

    #pragma omp simd aligned(g_xy_xz_xz_xx, g_xy_xz_xz_xy, g_xy_xz_xz_xz, g_xy_xz_xz_yy, g_xy_xz_xz_yz, g_xy_xz_xz_zz, g_y_z_0_0_x_x_xz_xx, g_y_z_0_0_x_x_xz_xy, g_y_z_0_0_x_x_xz_xz, g_y_z_0_0_x_x_xz_yy, g_y_z_0_0_x_x_xz_yz, g_y_z_0_0_x_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_x_xz_xx[i] = 4.0 * g_xy_xz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_xz_xy[i] = 4.0 * g_xy_xz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_xz_xz[i] = 4.0 * g_xy_xz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_xz_yy[i] = 4.0 * g_xy_xz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_xz_yz[i] = 4.0 * g_xy_xz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_xz_zz[i] = 4.0 * g_xy_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1638-1644)

    #pragma omp simd aligned(g_xy_xz_yy_xx, g_xy_xz_yy_xy, g_xy_xz_yy_xz, g_xy_xz_yy_yy, g_xy_xz_yy_yz, g_xy_xz_yy_zz, g_y_z_0_0_x_x_yy_xx, g_y_z_0_0_x_x_yy_xy, g_y_z_0_0_x_x_yy_xz, g_y_z_0_0_x_x_yy_yy, g_y_z_0_0_x_x_yy_yz, g_y_z_0_0_x_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_x_yy_xx[i] = 4.0 * g_xy_xz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_yy_xy[i] = 4.0 * g_xy_xz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_yy_xz[i] = 4.0 * g_xy_xz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_yy_yy[i] = 4.0 * g_xy_xz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_yy_yz[i] = 4.0 * g_xy_xz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_yy_zz[i] = 4.0 * g_xy_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1644-1650)

    #pragma omp simd aligned(g_xy_xz_yz_xx, g_xy_xz_yz_xy, g_xy_xz_yz_xz, g_xy_xz_yz_yy, g_xy_xz_yz_yz, g_xy_xz_yz_zz, g_y_z_0_0_x_x_yz_xx, g_y_z_0_0_x_x_yz_xy, g_y_z_0_0_x_x_yz_xz, g_y_z_0_0_x_x_yz_yy, g_y_z_0_0_x_x_yz_yz, g_y_z_0_0_x_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_x_yz_xx[i] = 4.0 * g_xy_xz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_yz_xy[i] = 4.0 * g_xy_xz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_yz_xz[i] = 4.0 * g_xy_xz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_yz_yy[i] = 4.0 * g_xy_xz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_yz_yz[i] = 4.0 * g_xy_xz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_yz_zz[i] = 4.0 * g_xy_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1650-1656)

    #pragma omp simd aligned(g_xy_xz_zz_xx, g_xy_xz_zz_xy, g_xy_xz_zz_xz, g_xy_xz_zz_yy, g_xy_xz_zz_yz, g_xy_xz_zz_zz, g_y_z_0_0_x_x_zz_xx, g_y_z_0_0_x_x_zz_xy, g_y_z_0_0_x_x_zz_xz, g_y_z_0_0_x_x_zz_yy, g_y_z_0_0_x_x_zz_yz, g_y_z_0_0_x_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_x_zz_xx[i] = 4.0 * g_xy_xz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_zz_xy[i] = 4.0 * g_xy_xz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_zz_xz[i] = 4.0 * g_xy_xz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_zz_yy[i] = 4.0 * g_xy_xz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_zz_yz[i] = 4.0 * g_xy_xz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_zz_zz[i] = 4.0 * g_xy_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1656-1662)

    #pragma omp simd aligned(g_xy_yz_xx_xx, g_xy_yz_xx_xy, g_xy_yz_xx_xz, g_xy_yz_xx_yy, g_xy_yz_xx_yz, g_xy_yz_xx_zz, g_y_z_0_0_x_y_xx_xx, g_y_z_0_0_x_y_xx_xy, g_y_z_0_0_x_y_xx_xz, g_y_z_0_0_x_y_xx_yy, g_y_z_0_0_x_y_xx_yz, g_y_z_0_0_x_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_y_xx_xx[i] = 4.0 * g_xy_yz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_xx_xy[i] = 4.0 * g_xy_yz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_xx_xz[i] = 4.0 * g_xy_yz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_xx_yy[i] = 4.0 * g_xy_yz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_xx_yz[i] = 4.0 * g_xy_yz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_xx_zz[i] = 4.0 * g_xy_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1662-1668)

    #pragma omp simd aligned(g_xy_yz_xy_xx, g_xy_yz_xy_xy, g_xy_yz_xy_xz, g_xy_yz_xy_yy, g_xy_yz_xy_yz, g_xy_yz_xy_zz, g_y_z_0_0_x_y_xy_xx, g_y_z_0_0_x_y_xy_xy, g_y_z_0_0_x_y_xy_xz, g_y_z_0_0_x_y_xy_yy, g_y_z_0_0_x_y_xy_yz, g_y_z_0_0_x_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_y_xy_xx[i] = 4.0 * g_xy_yz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_xy_xy[i] = 4.0 * g_xy_yz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_xy_xz[i] = 4.0 * g_xy_yz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_xy_yy[i] = 4.0 * g_xy_yz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_xy_yz[i] = 4.0 * g_xy_yz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_xy_zz[i] = 4.0 * g_xy_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1668-1674)

    #pragma omp simd aligned(g_xy_yz_xz_xx, g_xy_yz_xz_xy, g_xy_yz_xz_xz, g_xy_yz_xz_yy, g_xy_yz_xz_yz, g_xy_yz_xz_zz, g_y_z_0_0_x_y_xz_xx, g_y_z_0_0_x_y_xz_xy, g_y_z_0_0_x_y_xz_xz, g_y_z_0_0_x_y_xz_yy, g_y_z_0_0_x_y_xz_yz, g_y_z_0_0_x_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_y_xz_xx[i] = 4.0 * g_xy_yz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_xz_xy[i] = 4.0 * g_xy_yz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_xz_xz[i] = 4.0 * g_xy_yz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_xz_yy[i] = 4.0 * g_xy_yz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_xz_yz[i] = 4.0 * g_xy_yz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_xz_zz[i] = 4.0 * g_xy_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1674-1680)

    #pragma omp simd aligned(g_xy_yz_yy_xx, g_xy_yz_yy_xy, g_xy_yz_yy_xz, g_xy_yz_yy_yy, g_xy_yz_yy_yz, g_xy_yz_yy_zz, g_y_z_0_0_x_y_yy_xx, g_y_z_0_0_x_y_yy_xy, g_y_z_0_0_x_y_yy_xz, g_y_z_0_0_x_y_yy_yy, g_y_z_0_0_x_y_yy_yz, g_y_z_0_0_x_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_y_yy_xx[i] = 4.0 * g_xy_yz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_yy_xy[i] = 4.0 * g_xy_yz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_yy_xz[i] = 4.0 * g_xy_yz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_yy_yy[i] = 4.0 * g_xy_yz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_yy_yz[i] = 4.0 * g_xy_yz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_yy_zz[i] = 4.0 * g_xy_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1680-1686)

    #pragma omp simd aligned(g_xy_yz_yz_xx, g_xy_yz_yz_xy, g_xy_yz_yz_xz, g_xy_yz_yz_yy, g_xy_yz_yz_yz, g_xy_yz_yz_zz, g_y_z_0_0_x_y_yz_xx, g_y_z_0_0_x_y_yz_xy, g_y_z_0_0_x_y_yz_xz, g_y_z_0_0_x_y_yz_yy, g_y_z_0_0_x_y_yz_yz, g_y_z_0_0_x_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_y_yz_xx[i] = 4.0 * g_xy_yz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_yz_xy[i] = 4.0 * g_xy_yz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_yz_xz[i] = 4.0 * g_xy_yz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_yz_yy[i] = 4.0 * g_xy_yz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_yz_yz[i] = 4.0 * g_xy_yz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_yz_zz[i] = 4.0 * g_xy_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1686-1692)

    #pragma omp simd aligned(g_xy_yz_zz_xx, g_xy_yz_zz_xy, g_xy_yz_zz_xz, g_xy_yz_zz_yy, g_xy_yz_zz_yz, g_xy_yz_zz_zz, g_y_z_0_0_x_y_zz_xx, g_y_z_0_0_x_y_zz_xy, g_y_z_0_0_x_y_zz_xz, g_y_z_0_0_x_y_zz_yy, g_y_z_0_0_x_y_zz_yz, g_y_z_0_0_x_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_y_zz_xx[i] = 4.0 * g_xy_yz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_zz_xy[i] = 4.0 * g_xy_yz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_zz_xz[i] = 4.0 * g_xy_yz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_zz_yy[i] = 4.0 * g_xy_yz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_zz_yz[i] = 4.0 * g_xy_yz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_zz_zz[i] = 4.0 * g_xy_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1692-1698)

    #pragma omp simd aligned(g_xy_0_xx_xx, g_xy_0_xx_xy, g_xy_0_xx_xz, g_xy_0_xx_yy, g_xy_0_xx_yz, g_xy_0_xx_zz, g_xy_zz_xx_xx, g_xy_zz_xx_xy, g_xy_zz_xx_xz, g_xy_zz_xx_yy, g_xy_zz_xx_yz, g_xy_zz_xx_zz, g_y_z_0_0_x_z_xx_xx, g_y_z_0_0_x_z_xx_xy, g_y_z_0_0_x_z_xx_xz, g_y_z_0_0_x_z_xx_yy, g_y_z_0_0_x_z_xx_yz, g_y_z_0_0_x_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_z_xx_xx[i] = -2.0 * g_xy_0_xx_xx[i] * a_exp + 4.0 * g_xy_zz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_xx_xy[i] = -2.0 * g_xy_0_xx_xy[i] * a_exp + 4.0 * g_xy_zz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_xx_xz[i] = -2.0 * g_xy_0_xx_xz[i] * a_exp + 4.0 * g_xy_zz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_xx_yy[i] = -2.0 * g_xy_0_xx_yy[i] * a_exp + 4.0 * g_xy_zz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_xx_yz[i] = -2.0 * g_xy_0_xx_yz[i] * a_exp + 4.0 * g_xy_zz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_xx_zz[i] = -2.0 * g_xy_0_xx_zz[i] * a_exp + 4.0 * g_xy_zz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1698-1704)

    #pragma omp simd aligned(g_xy_0_xy_xx, g_xy_0_xy_xy, g_xy_0_xy_xz, g_xy_0_xy_yy, g_xy_0_xy_yz, g_xy_0_xy_zz, g_xy_zz_xy_xx, g_xy_zz_xy_xy, g_xy_zz_xy_xz, g_xy_zz_xy_yy, g_xy_zz_xy_yz, g_xy_zz_xy_zz, g_y_z_0_0_x_z_xy_xx, g_y_z_0_0_x_z_xy_xy, g_y_z_0_0_x_z_xy_xz, g_y_z_0_0_x_z_xy_yy, g_y_z_0_0_x_z_xy_yz, g_y_z_0_0_x_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_z_xy_xx[i] = -2.0 * g_xy_0_xy_xx[i] * a_exp + 4.0 * g_xy_zz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_xy_xy[i] = -2.0 * g_xy_0_xy_xy[i] * a_exp + 4.0 * g_xy_zz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_xy_xz[i] = -2.0 * g_xy_0_xy_xz[i] * a_exp + 4.0 * g_xy_zz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_xy_yy[i] = -2.0 * g_xy_0_xy_yy[i] * a_exp + 4.0 * g_xy_zz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_xy_yz[i] = -2.0 * g_xy_0_xy_yz[i] * a_exp + 4.0 * g_xy_zz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_xy_zz[i] = -2.0 * g_xy_0_xy_zz[i] * a_exp + 4.0 * g_xy_zz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1704-1710)

    #pragma omp simd aligned(g_xy_0_xz_xx, g_xy_0_xz_xy, g_xy_0_xz_xz, g_xy_0_xz_yy, g_xy_0_xz_yz, g_xy_0_xz_zz, g_xy_zz_xz_xx, g_xy_zz_xz_xy, g_xy_zz_xz_xz, g_xy_zz_xz_yy, g_xy_zz_xz_yz, g_xy_zz_xz_zz, g_y_z_0_0_x_z_xz_xx, g_y_z_0_0_x_z_xz_xy, g_y_z_0_0_x_z_xz_xz, g_y_z_0_0_x_z_xz_yy, g_y_z_0_0_x_z_xz_yz, g_y_z_0_0_x_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_z_xz_xx[i] = -2.0 * g_xy_0_xz_xx[i] * a_exp + 4.0 * g_xy_zz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_xz_xy[i] = -2.0 * g_xy_0_xz_xy[i] * a_exp + 4.0 * g_xy_zz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_xz_xz[i] = -2.0 * g_xy_0_xz_xz[i] * a_exp + 4.0 * g_xy_zz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_xz_yy[i] = -2.0 * g_xy_0_xz_yy[i] * a_exp + 4.0 * g_xy_zz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_xz_yz[i] = -2.0 * g_xy_0_xz_yz[i] * a_exp + 4.0 * g_xy_zz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_xz_zz[i] = -2.0 * g_xy_0_xz_zz[i] * a_exp + 4.0 * g_xy_zz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1710-1716)

    #pragma omp simd aligned(g_xy_0_yy_xx, g_xy_0_yy_xy, g_xy_0_yy_xz, g_xy_0_yy_yy, g_xy_0_yy_yz, g_xy_0_yy_zz, g_xy_zz_yy_xx, g_xy_zz_yy_xy, g_xy_zz_yy_xz, g_xy_zz_yy_yy, g_xy_zz_yy_yz, g_xy_zz_yy_zz, g_y_z_0_0_x_z_yy_xx, g_y_z_0_0_x_z_yy_xy, g_y_z_0_0_x_z_yy_xz, g_y_z_0_0_x_z_yy_yy, g_y_z_0_0_x_z_yy_yz, g_y_z_0_0_x_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_z_yy_xx[i] = -2.0 * g_xy_0_yy_xx[i] * a_exp + 4.0 * g_xy_zz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_yy_xy[i] = -2.0 * g_xy_0_yy_xy[i] * a_exp + 4.0 * g_xy_zz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_yy_xz[i] = -2.0 * g_xy_0_yy_xz[i] * a_exp + 4.0 * g_xy_zz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_yy_yy[i] = -2.0 * g_xy_0_yy_yy[i] * a_exp + 4.0 * g_xy_zz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_yy_yz[i] = -2.0 * g_xy_0_yy_yz[i] * a_exp + 4.0 * g_xy_zz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_yy_zz[i] = -2.0 * g_xy_0_yy_zz[i] * a_exp + 4.0 * g_xy_zz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1716-1722)

    #pragma omp simd aligned(g_xy_0_yz_xx, g_xy_0_yz_xy, g_xy_0_yz_xz, g_xy_0_yz_yy, g_xy_0_yz_yz, g_xy_0_yz_zz, g_xy_zz_yz_xx, g_xy_zz_yz_xy, g_xy_zz_yz_xz, g_xy_zz_yz_yy, g_xy_zz_yz_yz, g_xy_zz_yz_zz, g_y_z_0_0_x_z_yz_xx, g_y_z_0_0_x_z_yz_xy, g_y_z_0_0_x_z_yz_xz, g_y_z_0_0_x_z_yz_yy, g_y_z_0_0_x_z_yz_yz, g_y_z_0_0_x_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_z_yz_xx[i] = -2.0 * g_xy_0_yz_xx[i] * a_exp + 4.0 * g_xy_zz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_yz_xy[i] = -2.0 * g_xy_0_yz_xy[i] * a_exp + 4.0 * g_xy_zz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_yz_xz[i] = -2.0 * g_xy_0_yz_xz[i] * a_exp + 4.0 * g_xy_zz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_yz_yy[i] = -2.0 * g_xy_0_yz_yy[i] * a_exp + 4.0 * g_xy_zz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_yz_yz[i] = -2.0 * g_xy_0_yz_yz[i] * a_exp + 4.0 * g_xy_zz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_yz_zz[i] = -2.0 * g_xy_0_yz_zz[i] * a_exp + 4.0 * g_xy_zz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1722-1728)

    #pragma omp simd aligned(g_xy_0_zz_xx, g_xy_0_zz_xy, g_xy_0_zz_xz, g_xy_0_zz_yy, g_xy_0_zz_yz, g_xy_0_zz_zz, g_xy_zz_zz_xx, g_xy_zz_zz_xy, g_xy_zz_zz_xz, g_xy_zz_zz_yy, g_xy_zz_zz_yz, g_xy_zz_zz_zz, g_y_z_0_0_x_z_zz_xx, g_y_z_0_0_x_z_zz_xy, g_y_z_0_0_x_z_zz_xz, g_y_z_0_0_x_z_zz_yy, g_y_z_0_0_x_z_zz_yz, g_y_z_0_0_x_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_z_zz_xx[i] = -2.0 * g_xy_0_zz_xx[i] * a_exp + 4.0 * g_xy_zz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_zz_xy[i] = -2.0 * g_xy_0_zz_xy[i] * a_exp + 4.0 * g_xy_zz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_zz_xz[i] = -2.0 * g_xy_0_zz_xz[i] * a_exp + 4.0 * g_xy_zz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_zz_yy[i] = -2.0 * g_xy_0_zz_yy[i] * a_exp + 4.0 * g_xy_zz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_zz_yz[i] = -2.0 * g_xy_0_zz_yz[i] * a_exp + 4.0 * g_xy_zz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_zz_zz[i] = -2.0 * g_xy_0_zz_zz[i] * a_exp + 4.0 * g_xy_zz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1728-1734)

    #pragma omp simd aligned(g_0_xz_xx_xx, g_0_xz_xx_xy, g_0_xz_xx_xz, g_0_xz_xx_yy, g_0_xz_xx_yz, g_0_xz_xx_zz, g_y_z_0_0_y_x_xx_xx, g_y_z_0_0_y_x_xx_xy, g_y_z_0_0_y_x_xx_xz, g_y_z_0_0_y_x_xx_yy, g_y_z_0_0_y_x_xx_yz, g_y_z_0_0_y_x_xx_zz, g_yy_xz_xx_xx, g_yy_xz_xx_xy, g_yy_xz_xx_xz, g_yy_xz_xx_yy, g_yy_xz_xx_yz, g_yy_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_x_xx_xx[i] = -2.0 * g_0_xz_xx_xx[i] * b_exp + 4.0 * g_yy_xz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_xx_xy[i] = -2.0 * g_0_xz_xx_xy[i] * b_exp + 4.0 * g_yy_xz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_xx_xz[i] = -2.0 * g_0_xz_xx_xz[i] * b_exp + 4.0 * g_yy_xz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_xx_yy[i] = -2.0 * g_0_xz_xx_yy[i] * b_exp + 4.0 * g_yy_xz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_xx_yz[i] = -2.0 * g_0_xz_xx_yz[i] * b_exp + 4.0 * g_yy_xz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_xx_zz[i] = -2.0 * g_0_xz_xx_zz[i] * b_exp + 4.0 * g_yy_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1734-1740)

    #pragma omp simd aligned(g_0_xz_xy_xx, g_0_xz_xy_xy, g_0_xz_xy_xz, g_0_xz_xy_yy, g_0_xz_xy_yz, g_0_xz_xy_zz, g_y_z_0_0_y_x_xy_xx, g_y_z_0_0_y_x_xy_xy, g_y_z_0_0_y_x_xy_xz, g_y_z_0_0_y_x_xy_yy, g_y_z_0_0_y_x_xy_yz, g_y_z_0_0_y_x_xy_zz, g_yy_xz_xy_xx, g_yy_xz_xy_xy, g_yy_xz_xy_xz, g_yy_xz_xy_yy, g_yy_xz_xy_yz, g_yy_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_x_xy_xx[i] = -2.0 * g_0_xz_xy_xx[i] * b_exp + 4.0 * g_yy_xz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_xy_xy[i] = -2.0 * g_0_xz_xy_xy[i] * b_exp + 4.0 * g_yy_xz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_xy_xz[i] = -2.0 * g_0_xz_xy_xz[i] * b_exp + 4.0 * g_yy_xz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_xy_yy[i] = -2.0 * g_0_xz_xy_yy[i] * b_exp + 4.0 * g_yy_xz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_xy_yz[i] = -2.0 * g_0_xz_xy_yz[i] * b_exp + 4.0 * g_yy_xz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_xy_zz[i] = -2.0 * g_0_xz_xy_zz[i] * b_exp + 4.0 * g_yy_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1740-1746)

    #pragma omp simd aligned(g_0_xz_xz_xx, g_0_xz_xz_xy, g_0_xz_xz_xz, g_0_xz_xz_yy, g_0_xz_xz_yz, g_0_xz_xz_zz, g_y_z_0_0_y_x_xz_xx, g_y_z_0_0_y_x_xz_xy, g_y_z_0_0_y_x_xz_xz, g_y_z_0_0_y_x_xz_yy, g_y_z_0_0_y_x_xz_yz, g_y_z_0_0_y_x_xz_zz, g_yy_xz_xz_xx, g_yy_xz_xz_xy, g_yy_xz_xz_xz, g_yy_xz_xz_yy, g_yy_xz_xz_yz, g_yy_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_x_xz_xx[i] = -2.0 * g_0_xz_xz_xx[i] * b_exp + 4.0 * g_yy_xz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_xz_xy[i] = -2.0 * g_0_xz_xz_xy[i] * b_exp + 4.0 * g_yy_xz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_xz_xz[i] = -2.0 * g_0_xz_xz_xz[i] * b_exp + 4.0 * g_yy_xz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_xz_yy[i] = -2.0 * g_0_xz_xz_yy[i] * b_exp + 4.0 * g_yy_xz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_xz_yz[i] = -2.0 * g_0_xz_xz_yz[i] * b_exp + 4.0 * g_yy_xz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_xz_zz[i] = -2.0 * g_0_xz_xz_zz[i] * b_exp + 4.0 * g_yy_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1746-1752)

    #pragma omp simd aligned(g_0_xz_yy_xx, g_0_xz_yy_xy, g_0_xz_yy_xz, g_0_xz_yy_yy, g_0_xz_yy_yz, g_0_xz_yy_zz, g_y_z_0_0_y_x_yy_xx, g_y_z_0_0_y_x_yy_xy, g_y_z_0_0_y_x_yy_xz, g_y_z_0_0_y_x_yy_yy, g_y_z_0_0_y_x_yy_yz, g_y_z_0_0_y_x_yy_zz, g_yy_xz_yy_xx, g_yy_xz_yy_xy, g_yy_xz_yy_xz, g_yy_xz_yy_yy, g_yy_xz_yy_yz, g_yy_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_x_yy_xx[i] = -2.0 * g_0_xz_yy_xx[i] * b_exp + 4.0 * g_yy_xz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_yy_xy[i] = -2.0 * g_0_xz_yy_xy[i] * b_exp + 4.0 * g_yy_xz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_yy_xz[i] = -2.0 * g_0_xz_yy_xz[i] * b_exp + 4.0 * g_yy_xz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_yy_yy[i] = -2.0 * g_0_xz_yy_yy[i] * b_exp + 4.0 * g_yy_xz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_yy_yz[i] = -2.0 * g_0_xz_yy_yz[i] * b_exp + 4.0 * g_yy_xz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_yy_zz[i] = -2.0 * g_0_xz_yy_zz[i] * b_exp + 4.0 * g_yy_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1752-1758)

    #pragma omp simd aligned(g_0_xz_yz_xx, g_0_xz_yz_xy, g_0_xz_yz_xz, g_0_xz_yz_yy, g_0_xz_yz_yz, g_0_xz_yz_zz, g_y_z_0_0_y_x_yz_xx, g_y_z_0_0_y_x_yz_xy, g_y_z_0_0_y_x_yz_xz, g_y_z_0_0_y_x_yz_yy, g_y_z_0_0_y_x_yz_yz, g_y_z_0_0_y_x_yz_zz, g_yy_xz_yz_xx, g_yy_xz_yz_xy, g_yy_xz_yz_xz, g_yy_xz_yz_yy, g_yy_xz_yz_yz, g_yy_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_x_yz_xx[i] = -2.0 * g_0_xz_yz_xx[i] * b_exp + 4.0 * g_yy_xz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_yz_xy[i] = -2.0 * g_0_xz_yz_xy[i] * b_exp + 4.0 * g_yy_xz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_yz_xz[i] = -2.0 * g_0_xz_yz_xz[i] * b_exp + 4.0 * g_yy_xz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_yz_yy[i] = -2.0 * g_0_xz_yz_yy[i] * b_exp + 4.0 * g_yy_xz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_yz_yz[i] = -2.0 * g_0_xz_yz_yz[i] * b_exp + 4.0 * g_yy_xz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_yz_zz[i] = -2.0 * g_0_xz_yz_zz[i] * b_exp + 4.0 * g_yy_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1758-1764)

    #pragma omp simd aligned(g_0_xz_zz_xx, g_0_xz_zz_xy, g_0_xz_zz_xz, g_0_xz_zz_yy, g_0_xz_zz_yz, g_0_xz_zz_zz, g_y_z_0_0_y_x_zz_xx, g_y_z_0_0_y_x_zz_xy, g_y_z_0_0_y_x_zz_xz, g_y_z_0_0_y_x_zz_yy, g_y_z_0_0_y_x_zz_yz, g_y_z_0_0_y_x_zz_zz, g_yy_xz_zz_xx, g_yy_xz_zz_xy, g_yy_xz_zz_xz, g_yy_xz_zz_yy, g_yy_xz_zz_yz, g_yy_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_x_zz_xx[i] = -2.0 * g_0_xz_zz_xx[i] * b_exp + 4.0 * g_yy_xz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_zz_xy[i] = -2.0 * g_0_xz_zz_xy[i] * b_exp + 4.0 * g_yy_xz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_zz_xz[i] = -2.0 * g_0_xz_zz_xz[i] * b_exp + 4.0 * g_yy_xz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_zz_yy[i] = -2.0 * g_0_xz_zz_yy[i] * b_exp + 4.0 * g_yy_xz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_zz_yz[i] = -2.0 * g_0_xz_zz_yz[i] * b_exp + 4.0 * g_yy_xz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_zz_zz[i] = -2.0 * g_0_xz_zz_zz[i] * b_exp + 4.0 * g_yy_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1764-1770)

    #pragma omp simd aligned(g_0_yz_xx_xx, g_0_yz_xx_xy, g_0_yz_xx_xz, g_0_yz_xx_yy, g_0_yz_xx_yz, g_0_yz_xx_zz, g_y_z_0_0_y_y_xx_xx, g_y_z_0_0_y_y_xx_xy, g_y_z_0_0_y_y_xx_xz, g_y_z_0_0_y_y_xx_yy, g_y_z_0_0_y_y_xx_yz, g_y_z_0_0_y_y_xx_zz, g_yy_yz_xx_xx, g_yy_yz_xx_xy, g_yy_yz_xx_xz, g_yy_yz_xx_yy, g_yy_yz_xx_yz, g_yy_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_y_xx_xx[i] = -2.0 * g_0_yz_xx_xx[i] * b_exp + 4.0 * g_yy_yz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_xx_xy[i] = -2.0 * g_0_yz_xx_xy[i] * b_exp + 4.0 * g_yy_yz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_xx_xz[i] = -2.0 * g_0_yz_xx_xz[i] * b_exp + 4.0 * g_yy_yz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_xx_yy[i] = -2.0 * g_0_yz_xx_yy[i] * b_exp + 4.0 * g_yy_yz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_xx_yz[i] = -2.0 * g_0_yz_xx_yz[i] * b_exp + 4.0 * g_yy_yz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_xx_zz[i] = -2.0 * g_0_yz_xx_zz[i] * b_exp + 4.0 * g_yy_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1770-1776)

    #pragma omp simd aligned(g_0_yz_xy_xx, g_0_yz_xy_xy, g_0_yz_xy_xz, g_0_yz_xy_yy, g_0_yz_xy_yz, g_0_yz_xy_zz, g_y_z_0_0_y_y_xy_xx, g_y_z_0_0_y_y_xy_xy, g_y_z_0_0_y_y_xy_xz, g_y_z_0_0_y_y_xy_yy, g_y_z_0_0_y_y_xy_yz, g_y_z_0_0_y_y_xy_zz, g_yy_yz_xy_xx, g_yy_yz_xy_xy, g_yy_yz_xy_xz, g_yy_yz_xy_yy, g_yy_yz_xy_yz, g_yy_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_y_xy_xx[i] = -2.0 * g_0_yz_xy_xx[i] * b_exp + 4.0 * g_yy_yz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_xy_xy[i] = -2.0 * g_0_yz_xy_xy[i] * b_exp + 4.0 * g_yy_yz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_xy_xz[i] = -2.0 * g_0_yz_xy_xz[i] * b_exp + 4.0 * g_yy_yz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_xy_yy[i] = -2.0 * g_0_yz_xy_yy[i] * b_exp + 4.0 * g_yy_yz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_xy_yz[i] = -2.0 * g_0_yz_xy_yz[i] * b_exp + 4.0 * g_yy_yz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_xy_zz[i] = -2.0 * g_0_yz_xy_zz[i] * b_exp + 4.0 * g_yy_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1776-1782)

    #pragma omp simd aligned(g_0_yz_xz_xx, g_0_yz_xz_xy, g_0_yz_xz_xz, g_0_yz_xz_yy, g_0_yz_xz_yz, g_0_yz_xz_zz, g_y_z_0_0_y_y_xz_xx, g_y_z_0_0_y_y_xz_xy, g_y_z_0_0_y_y_xz_xz, g_y_z_0_0_y_y_xz_yy, g_y_z_0_0_y_y_xz_yz, g_y_z_0_0_y_y_xz_zz, g_yy_yz_xz_xx, g_yy_yz_xz_xy, g_yy_yz_xz_xz, g_yy_yz_xz_yy, g_yy_yz_xz_yz, g_yy_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_y_xz_xx[i] = -2.0 * g_0_yz_xz_xx[i] * b_exp + 4.0 * g_yy_yz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_xz_xy[i] = -2.0 * g_0_yz_xz_xy[i] * b_exp + 4.0 * g_yy_yz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_xz_xz[i] = -2.0 * g_0_yz_xz_xz[i] * b_exp + 4.0 * g_yy_yz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_xz_yy[i] = -2.0 * g_0_yz_xz_yy[i] * b_exp + 4.0 * g_yy_yz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_xz_yz[i] = -2.0 * g_0_yz_xz_yz[i] * b_exp + 4.0 * g_yy_yz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_xz_zz[i] = -2.0 * g_0_yz_xz_zz[i] * b_exp + 4.0 * g_yy_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1782-1788)

    #pragma omp simd aligned(g_0_yz_yy_xx, g_0_yz_yy_xy, g_0_yz_yy_xz, g_0_yz_yy_yy, g_0_yz_yy_yz, g_0_yz_yy_zz, g_y_z_0_0_y_y_yy_xx, g_y_z_0_0_y_y_yy_xy, g_y_z_0_0_y_y_yy_xz, g_y_z_0_0_y_y_yy_yy, g_y_z_0_0_y_y_yy_yz, g_y_z_0_0_y_y_yy_zz, g_yy_yz_yy_xx, g_yy_yz_yy_xy, g_yy_yz_yy_xz, g_yy_yz_yy_yy, g_yy_yz_yy_yz, g_yy_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_y_yy_xx[i] = -2.0 * g_0_yz_yy_xx[i] * b_exp + 4.0 * g_yy_yz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_yy_xy[i] = -2.0 * g_0_yz_yy_xy[i] * b_exp + 4.0 * g_yy_yz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_yy_xz[i] = -2.0 * g_0_yz_yy_xz[i] * b_exp + 4.0 * g_yy_yz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_yy_yy[i] = -2.0 * g_0_yz_yy_yy[i] * b_exp + 4.0 * g_yy_yz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_yy_yz[i] = -2.0 * g_0_yz_yy_yz[i] * b_exp + 4.0 * g_yy_yz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_yy_zz[i] = -2.0 * g_0_yz_yy_zz[i] * b_exp + 4.0 * g_yy_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1788-1794)

    #pragma omp simd aligned(g_0_yz_yz_xx, g_0_yz_yz_xy, g_0_yz_yz_xz, g_0_yz_yz_yy, g_0_yz_yz_yz, g_0_yz_yz_zz, g_y_z_0_0_y_y_yz_xx, g_y_z_0_0_y_y_yz_xy, g_y_z_0_0_y_y_yz_xz, g_y_z_0_0_y_y_yz_yy, g_y_z_0_0_y_y_yz_yz, g_y_z_0_0_y_y_yz_zz, g_yy_yz_yz_xx, g_yy_yz_yz_xy, g_yy_yz_yz_xz, g_yy_yz_yz_yy, g_yy_yz_yz_yz, g_yy_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_y_yz_xx[i] = -2.0 * g_0_yz_yz_xx[i] * b_exp + 4.0 * g_yy_yz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_yz_xy[i] = -2.0 * g_0_yz_yz_xy[i] * b_exp + 4.0 * g_yy_yz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_yz_xz[i] = -2.0 * g_0_yz_yz_xz[i] * b_exp + 4.0 * g_yy_yz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_yz_yy[i] = -2.0 * g_0_yz_yz_yy[i] * b_exp + 4.0 * g_yy_yz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_yz_yz[i] = -2.0 * g_0_yz_yz_yz[i] * b_exp + 4.0 * g_yy_yz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_yz_zz[i] = -2.0 * g_0_yz_yz_zz[i] * b_exp + 4.0 * g_yy_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1794-1800)

    #pragma omp simd aligned(g_0_yz_zz_xx, g_0_yz_zz_xy, g_0_yz_zz_xz, g_0_yz_zz_yy, g_0_yz_zz_yz, g_0_yz_zz_zz, g_y_z_0_0_y_y_zz_xx, g_y_z_0_0_y_y_zz_xy, g_y_z_0_0_y_y_zz_xz, g_y_z_0_0_y_y_zz_yy, g_y_z_0_0_y_y_zz_yz, g_y_z_0_0_y_y_zz_zz, g_yy_yz_zz_xx, g_yy_yz_zz_xy, g_yy_yz_zz_xz, g_yy_yz_zz_yy, g_yy_yz_zz_yz, g_yy_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_y_zz_xx[i] = -2.0 * g_0_yz_zz_xx[i] * b_exp + 4.0 * g_yy_yz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_zz_xy[i] = -2.0 * g_0_yz_zz_xy[i] * b_exp + 4.0 * g_yy_yz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_zz_xz[i] = -2.0 * g_0_yz_zz_xz[i] * b_exp + 4.0 * g_yy_yz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_zz_yy[i] = -2.0 * g_0_yz_zz_yy[i] * b_exp + 4.0 * g_yy_yz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_zz_yz[i] = -2.0 * g_0_yz_zz_yz[i] * b_exp + 4.0 * g_yy_yz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_zz_zz[i] = -2.0 * g_0_yz_zz_zz[i] * b_exp + 4.0 * g_yy_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1800-1806)

    #pragma omp simd aligned(g_0_0_xx_xx, g_0_0_xx_xy, g_0_0_xx_xz, g_0_0_xx_yy, g_0_0_xx_yz, g_0_0_xx_zz, g_0_zz_xx_xx, g_0_zz_xx_xy, g_0_zz_xx_xz, g_0_zz_xx_yy, g_0_zz_xx_yz, g_0_zz_xx_zz, g_y_z_0_0_y_z_xx_xx, g_y_z_0_0_y_z_xx_xy, g_y_z_0_0_y_z_xx_xz, g_y_z_0_0_y_z_xx_yy, g_y_z_0_0_y_z_xx_yz, g_y_z_0_0_y_z_xx_zz, g_yy_0_xx_xx, g_yy_0_xx_xy, g_yy_0_xx_xz, g_yy_0_xx_yy, g_yy_0_xx_yz, g_yy_0_xx_zz, g_yy_zz_xx_xx, g_yy_zz_xx_xy, g_yy_zz_xx_xz, g_yy_zz_xx_yy, g_yy_zz_xx_yz, g_yy_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_z_xx_xx[i] = g_0_0_xx_xx[i] - 2.0 * g_0_zz_xx_xx[i] * b_exp - 2.0 * g_yy_0_xx_xx[i] * a_exp + 4.0 * g_yy_zz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_xx_xy[i] = g_0_0_xx_xy[i] - 2.0 * g_0_zz_xx_xy[i] * b_exp - 2.0 * g_yy_0_xx_xy[i] * a_exp + 4.0 * g_yy_zz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_xx_xz[i] = g_0_0_xx_xz[i] - 2.0 * g_0_zz_xx_xz[i] * b_exp - 2.0 * g_yy_0_xx_xz[i] * a_exp + 4.0 * g_yy_zz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_xx_yy[i] = g_0_0_xx_yy[i] - 2.0 * g_0_zz_xx_yy[i] * b_exp - 2.0 * g_yy_0_xx_yy[i] * a_exp + 4.0 * g_yy_zz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_xx_yz[i] = g_0_0_xx_yz[i] - 2.0 * g_0_zz_xx_yz[i] * b_exp - 2.0 * g_yy_0_xx_yz[i] * a_exp + 4.0 * g_yy_zz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_xx_zz[i] = g_0_0_xx_zz[i] - 2.0 * g_0_zz_xx_zz[i] * b_exp - 2.0 * g_yy_0_xx_zz[i] * a_exp + 4.0 * g_yy_zz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1806-1812)

    #pragma omp simd aligned(g_0_0_xy_xx, g_0_0_xy_xy, g_0_0_xy_xz, g_0_0_xy_yy, g_0_0_xy_yz, g_0_0_xy_zz, g_0_zz_xy_xx, g_0_zz_xy_xy, g_0_zz_xy_xz, g_0_zz_xy_yy, g_0_zz_xy_yz, g_0_zz_xy_zz, g_y_z_0_0_y_z_xy_xx, g_y_z_0_0_y_z_xy_xy, g_y_z_0_0_y_z_xy_xz, g_y_z_0_0_y_z_xy_yy, g_y_z_0_0_y_z_xy_yz, g_y_z_0_0_y_z_xy_zz, g_yy_0_xy_xx, g_yy_0_xy_xy, g_yy_0_xy_xz, g_yy_0_xy_yy, g_yy_0_xy_yz, g_yy_0_xy_zz, g_yy_zz_xy_xx, g_yy_zz_xy_xy, g_yy_zz_xy_xz, g_yy_zz_xy_yy, g_yy_zz_xy_yz, g_yy_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_z_xy_xx[i] = g_0_0_xy_xx[i] - 2.0 * g_0_zz_xy_xx[i] * b_exp - 2.0 * g_yy_0_xy_xx[i] * a_exp + 4.0 * g_yy_zz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_xy_xy[i] = g_0_0_xy_xy[i] - 2.0 * g_0_zz_xy_xy[i] * b_exp - 2.0 * g_yy_0_xy_xy[i] * a_exp + 4.0 * g_yy_zz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_xy_xz[i] = g_0_0_xy_xz[i] - 2.0 * g_0_zz_xy_xz[i] * b_exp - 2.0 * g_yy_0_xy_xz[i] * a_exp + 4.0 * g_yy_zz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_xy_yy[i] = g_0_0_xy_yy[i] - 2.0 * g_0_zz_xy_yy[i] * b_exp - 2.0 * g_yy_0_xy_yy[i] * a_exp + 4.0 * g_yy_zz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_xy_yz[i] = g_0_0_xy_yz[i] - 2.0 * g_0_zz_xy_yz[i] * b_exp - 2.0 * g_yy_0_xy_yz[i] * a_exp + 4.0 * g_yy_zz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_xy_zz[i] = g_0_0_xy_zz[i] - 2.0 * g_0_zz_xy_zz[i] * b_exp - 2.0 * g_yy_0_xy_zz[i] * a_exp + 4.0 * g_yy_zz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1812-1818)

    #pragma omp simd aligned(g_0_0_xz_xx, g_0_0_xz_xy, g_0_0_xz_xz, g_0_0_xz_yy, g_0_0_xz_yz, g_0_0_xz_zz, g_0_zz_xz_xx, g_0_zz_xz_xy, g_0_zz_xz_xz, g_0_zz_xz_yy, g_0_zz_xz_yz, g_0_zz_xz_zz, g_y_z_0_0_y_z_xz_xx, g_y_z_0_0_y_z_xz_xy, g_y_z_0_0_y_z_xz_xz, g_y_z_0_0_y_z_xz_yy, g_y_z_0_0_y_z_xz_yz, g_y_z_0_0_y_z_xz_zz, g_yy_0_xz_xx, g_yy_0_xz_xy, g_yy_0_xz_xz, g_yy_0_xz_yy, g_yy_0_xz_yz, g_yy_0_xz_zz, g_yy_zz_xz_xx, g_yy_zz_xz_xy, g_yy_zz_xz_xz, g_yy_zz_xz_yy, g_yy_zz_xz_yz, g_yy_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_z_xz_xx[i] = g_0_0_xz_xx[i] - 2.0 * g_0_zz_xz_xx[i] * b_exp - 2.0 * g_yy_0_xz_xx[i] * a_exp + 4.0 * g_yy_zz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_xz_xy[i] = g_0_0_xz_xy[i] - 2.0 * g_0_zz_xz_xy[i] * b_exp - 2.0 * g_yy_0_xz_xy[i] * a_exp + 4.0 * g_yy_zz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_xz_xz[i] = g_0_0_xz_xz[i] - 2.0 * g_0_zz_xz_xz[i] * b_exp - 2.0 * g_yy_0_xz_xz[i] * a_exp + 4.0 * g_yy_zz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_xz_yy[i] = g_0_0_xz_yy[i] - 2.0 * g_0_zz_xz_yy[i] * b_exp - 2.0 * g_yy_0_xz_yy[i] * a_exp + 4.0 * g_yy_zz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_xz_yz[i] = g_0_0_xz_yz[i] - 2.0 * g_0_zz_xz_yz[i] * b_exp - 2.0 * g_yy_0_xz_yz[i] * a_exp + 4.0 * g_yy_zz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_xz_zz[i] = g_0_0_xz_zz[i] - 2.0 * g_0_zz_xz_zz[i] * b_exp - 2.0 * g_yy_0_xz_zz[i] * a_exp + 4.0 * g_yy_zz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1818-1824)

    #pragma omp simd aligned(g_0_0_yy_xx, g_0_0_yy_xy, g_0_0_yy_xz, g_0_0_yy_yy, g_0_0_yy_yz, g_0_0_yy_zz, g_0_zz_yy_xx, g_0_zz_yy_xy, g_0_zz_yy_xz, g_0_zz_yy_yy, g_0_zz_yy_yz, g_0_zz_yy_zz, g_y_z_0_0_y_z_yy_xx, g_y_z_0_0_y_z_yy_xy, g_y_z_0_0_y_z_yy_xz, g_y_z_0_0_y_z_yy_yy, g_y_z_0_0_y_z_yy_yz, g_y_z_0_0_y_z_yy_zz, g_yy_0_yy_xx, g_yy_0_yy_xy, g_yy_0_yy_xz, g_yy_0_yy_yy, g_yy_0_yy_yz, g_yy_0_yy_zz, g_yy_zz_yy_xx, g_yy_zz_yy_xy, g_yy_zz_yy_xz, g_yy_zz_yy_yy, g_yy_zz_yy_yz, g_yy_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_z_yy_xx[i] = g_0_0_yy_xx[i] - 2.0 * g_0_zz_yy_xx[i] * b_exp - 2.0 * g_yy_0_yy_xx[i] * a_exp + 4.0 * g_yy_zz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_yy_xy[i] = g_0_0_yy_xy[i] - 2.0 * g_0_zz_yy_xy[i] * b_exp - 2.0 * g_yy_0_yy_xy[i] * a_exp + 4.0 * g_yy_zz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_yy_xz[i] = g_0_0_yy_xz[i] - 2.0 * g_0_zz_yy_xz[i] * b_exp - 2.0 * g_yy_0_yy_xz[i] * a_exp + 4.0 * g_yy_zz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_yy_yy[i] = g_0_0_yy_yy[i] - 2.0 * g_0_zz_yy_yy[i] * b_exp - 2.0 * g_yy_0_yy_yy[i] * a_exp + 4.0 * g_yy_zz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_yy_yz[i] = g_0_0_yy_yz[i] - 2.0 * g_0_zz_yy_yz[i] * b_exp - 2.0 * g_yy_0_yy_yz[i] * a_exp + 4.0 * g_yy_zz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_yy_zz[i] = g_0_0_yy_zz[i] - 2.0 * g_0_zz_yy_zz[i] * b_exp - 2.0 * g_yy_0_yy_zz[i] * a_exp + 4.0 * g_yy_zz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1824-1830)

    #pragma omp simd aligned(g_0_0_yz_xx, g_0_0_yz_xy, g_0_0_yz_xz, g_0_0_yz_yy, g_0_0_yz_yz, g_0_0_yz_zz, g_0_zz_yz_xx, g_0_zz_yz_xy, g_0_zz_yz_xz, g_0_zz_yz_yy, g_0_zz_yz_yz, g_0_zz_yz_zz, g_y_z_0_0_y_z_yz_xx, g_y_z_0_0_y_z_yz_xy, g_y_z_0_0_y_z_yz_xz, g_y_z_0_0_y_z_yz_yy, g_y_z_0_0_y_z_yz_yz, g_y_z_0_0_y_z_yz_zz, g_yy_0_yz_xx, g_yy_0_yz_xy, g_yy_0_yz_xz, g_yy_0_yz_yy, g_yy_0_yz_yz, g_yy_0_yz_zz, g_yy_zz_yz_xx, g_yy_zz_yz_xy, g_yy_zz_yz_xz, g_yy_zz_yz_yy, g_yy_zz_yz_yz, g_yy_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_z_yz_xx[i] = g_0_0_yz_xx[i] - 2.0 * g_0_zz_yz_xx[i] * b_exp - 2.0 * g_yy_0_yz_xx[i] * a_exp + 4.0 * g_yy_zz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_yz_xy[i] = g_0_0_yz_xy[i] - 2.0 * g_0_zz_yz_xy[i] * b_exp - 2.0 * g_yy_0_yz_xy[i] * a_exp + 4.0 * g_yy_zz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_yz_xz[i] = g_0_0_yz_xz[i] - 2.0 * g_0_zz_yz_xz[i] * b_exp - 2.0 * g_yy_0_yz_xz[i] * a_exp + 4.0 * g_yy_zz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_yz_yy[i] = g_0_0_yz_yy[i] - 2.0 * g_0_zz_yz_yy[i] * b_exp - 2.0 * g_yy_0_yz_yy[i] * a_exp + 4.0 * g_yy_zz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_yz_yz[i] = g_0_0_yz_yz[i] - 2.0 * g_0_zz_yz_yz[i] * b_exp - 2.0 * g_yy_0_yz_yz[i] * a_exp + 4.0 * g_yy_zz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_yz_zz[i] = g_0_0_yz_zz[i] - 2.0 * g_0_zz_yz_zz[i] * b_exp - 2.0 * g_yy_0_yz_zz[i] * a_exp + 4.0 * g_yy_zz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1830-1836)

    #pragma omp simd aligned(g_0_0_zz_xx, g_0_0_zz_xy, g_0_0_zz_xz, g_0_0_zz_yy, g_0_0_zz_yz, g_0_0_zz_zz, g_0_zz_zz_xx, g_0_zz_zz_xy, g_0_zz_zz_xz, g_0_zz_zz_yy, g_0_zz_zz_yz, g_0_zz_zz_zz, g_y_z_0_0_y_z_zz_xx, g_y_z_0_0_y_z_zz_xy, g_y_z_0_0_y_z_zz_xz, g_y_z_0_0_y_z_zz_yy, g_y_z_0_0_y_z_zz_yz, g_y_z_0_0_y_z_zz_zz, g_yy_0_zz_xx, g_yy_0_zz_xy, g_yy_0_zz_xz, g_yy_0_zz_yy, g_yy_0_zz_yz, g_yy_0_zz_zz, g_yy_zz_zz_xx, g_yy_zz_zz_xy, g_yy_zz_zz_xz, g_yy_zz_zz_yy, g_yy_zz_zz_yz, g_yy_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_z_zz_xx[i] = g_0_0_zz_xx[i] - 2.0 * g_0_zz_zz_xx[i] * b_exp - 2.0 * g_yy_0_zz_xx[i] * a_exp + 4.0 * g_yy_zz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_zz_xy[i] = g_0_0_zz_xy[i] - 2.0 * g_0_zz_zz_xy[i] * b_exp - 2.0 * g_yy_0_zz_xy[i] * a_exp + 4.0 * g_yy_zz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_zz_xz[i] = g_0_0_zz_xz[i] - 2.0 * g_0_zz_zz_xz[i] * b_exp - 2.0 * g_yy_0_zz_xz[i] * a_exp + 4.0 * g_yy_zz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_zz_yy[i] = g_0_0_zz_yy[i] - 2.0 * g_0_zz_zz_yy[i] * b_exp - 2.0 * g_yy_0_zz_yy[i] * a_exp + 4.0 * g_yy_zz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_zz_yz[i] = g_0_0_zz_yz[i] - 2.0 * g_0_zz_zz_yz[i] * b_exp - 2.0 * g_yy_0_zz_yz[i] * a_exp + 4.0 * g_yy_zz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_zz_zz[i] = g_0_0_zz_zz[i] - 2.0 * g_0_zz_zz_zz[i] * b_exp - 2.0 * g_yy_0_zz_zz[i] * a_exp + 4.0 * g_yy_zz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1836-1842)

    #pragma omp simd aligned(g_y_z_0_0_z_x_xx_xx, g_y_z_0_0_z_x_xx_xy, g_y_z_0_0_z_x_xx_xz, g_y_z_0_0_z_x_xx_yy, g_y_z_0_0_z_x_xx_yz, g_y_z_0_0_z_x_xx_zz, g_yz_xz_xx_xx, g_yz_xz_xx_xy, g_yz_xz_xx_xz, g_yz_xz_xx_yy, g_yz_xz_xx_yz, g_yz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_x_xx_xx[i] = 4.0 * g_yz_xz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_xx_xy[i] = 4.0 * g_yz_xz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_xx_xz[i] = 4.0 * g_yz_xz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_xx_yy[i] = 4.0 * g_yz_xz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_xx_yz[i] = 4.0 * g_yz_xz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_xx_zz[i] = 4.0 * g_yz_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1842-1848)

    #pragma omp simd aligned(g_y_z_0_0_z_x_xy_xx, g_y_z_0_0_z_x_xy_xy, g_y_z_0_0_z_x_xy_xz, g_y_z_0_0_z_x_xy_yy, g_y_z_0_0_z_x_xy_yz, g_y_z_0_0_z_x_xy_zz, g_yz_xz_xy_xx, g_yz_xz_xy_xy, g_yz_xz_xy_xz, g_yz_xz_xy_yy, g_yz_xz_xy_yz, g_yz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_x_xy_xx[i] = 4.0 * g_yz_xz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_xy_xy[i] = 4.0 * g_yz_xz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_xy_xz[i] = 4.0 * g_yz_xz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_xy_yy[i] = 4.0 * g_yz_xz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_xy_yz[i] = 4.0 * g_yz_xz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_xy_zz[i] = 4.0 * g_yz_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1848-1854)

    #pragma omp simd aligned(g_y_z_0_0_z_x_xz_xx, g_y_z_0_0_z_x_xz_xy, g_y_z_0_0_z_x_xz_xz, g_y_z_0_0_z_x_xz_yy, g_y_z_0_0_z_x_xz_yz, g_y_z_0_0_z_x_xz_zz, g_yz_xz_xz_xx, g_yz_xz_xz_xy, g_yz_xz_xz_xz, g_yz_xz_xz_yy, g_yz_xz_xz_yz, g_yz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_x_xz_xx[i] = 4.0 * g_yz_xz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_xz_xy[i] = 4.0 * g_yz_xz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_xz_xz[i] = 4.0 * g_yz_xz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_xz_yy[i] = 4.0 * g_yz_xz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_xz_yz[i] = 4.0 * g_yz_xz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_xz_zz[i] = 4.0 * g_yz_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1854-1860)

    #pragma omp simd aligned(g_y_z_0_0_z_x_yy_xx, g_y_z_0_0_z_x_yy_xy, g_y_z_0_0_z_x_yy_xz, g_y_z_0_0_z_x_yy_yy, g_y_z_0_0_z_x_yy_yz, g_y_z_0_0_z_x_yy_zz, g_yz_xz_yy_xx, g_yz_xz_yy_xy, g_yz_xz_yy_xz, g_yz_xz_yy_yy, g_yz_xz_yy_yz, g_yz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_x_yy_xx[i] = 4.0 * g_yz_xz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_yy_xy[i] = 4.0 * g_yz_xz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_yy_xz[i] = 4.0 * g_yz_xz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_yy_yy[i] = 4.0 * g_yz_xz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_yy_yz[i] = 4.0 * g_yz_xz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_yy_zz[i] = 4.0 * g_yz_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1860-1866)

    #pragma omp simd aligned(g_y_z_0_0_z_x_yz_xx, g_y_z_0_0_z_x_yz_xy, g_y_z_0_0_z_x_yz_xz, g_y_z_0_0_z_x_yz_yy, g_y_z_0_0_z_x_yz_yz, g_y_z_0_0_z_x_yz_zz, g_yz_xz_yz_xx, g_yz_xz_yz_xy, g_yz_xz_yz_xz, g_yz_xz_yz_yy, g_yz_xz_yz_yz, g_yz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_x_yz_xx[i] = 4.0 * g_yz_xz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_yz_xy[i] = 4.0 * g_yz_xz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_yz_xz[i] = 4.0 * g_yz_xz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_yz_yy[i] = 4.0 * g_yz_xz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_yz_yz[i] = 4.0 * g_yz_xz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_yz_zz[i] = 4.0 * g_yz_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1866-1872)

    #pragma omp simd aligned(g_y_z_0_0_z_x_zz_xx, g_y_z_0_0_z_x_zz_xy, g_y_z_0_0_z_x_zz_xz, g_y_z_0_0_z_x_zz_yy, g_y_z_0_0_z_x_zz_yz, g_y_z_0_0_z_x_zz_zz, g_yz_xz_zz_xx, g_yz_xz_zz_xy, g_yz_xz_zz_xz, g_yz_xz_zz_yy, g_yz_xz_zz_yz, g_yz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_x_zz_xx[i] = 4.0 * g_yz_xz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_zz_xy[i] = 4.0 * g_yz_xz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_zz_xz[i] = 4.0 * g_yz_xz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_zz_yy[i] = 4.0 * g_yz_xz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_zz_yz[i] = 4.0 * g_yz_xz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_zz_zz[i] = 4.0 * g_yz_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1872-1878)

    #pragma omp simd aligned(g_y_z_0_0_z_y_xx_xx, g_y_z_0_0_z_y_xx_xy, g_y_z_0_0_z_y_xx_xz, g_y_z_0_0_z_y_xx_yy, g_y_z_0_0_z_y_xx_yz, g_y_z_0_0_z_y_xx_zz, g_yz_yz_xx_xx, g_yz_yz_xx_xy, g_yz_yz_xx_xz, g_yz_yz_xx_yy, g_yz_yz_xx_yz, g_yz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_y_xx_xx[i] = 4.0 * g_yz_yz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_xx_xy[i] = 4.0 * g_yz_yz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_xx_xz[i] = 4.0 * g_yz_yz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_xx_yy[i] = 4.0 * g_yz_yz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_xx_yz[i] = 4.0 * g_yz_yz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_xx_zz[i] = 4.0 * g_yz_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1878-1884)

    #pragma omp simd aligned(g_y_z_0_0_z_y_xy_xx, g_y_z_0_0_z_y_xy_xy, g_y_z_0_0_z_y_xy_xz, g_y_z_0_0_z_y_xy_yy, g_y_z_0_0_z_y_xy_yz, g_y_z_0_0_z_y_xy_zz, g_yz_yz_xy_xx, g_yz_yz_xy_xy, g_yz_yz_xy_xz, g_yz_yz_xy_yy, g_yz_yz_xy_yz, g_yz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_y_xy_xx[i] = 4.0 * g_yz_yz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_xy_xy[i] = 4.0 * g_yz_yz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_xy_xz[i] = 4.0 * g_yz_yz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_xy_yy[i] = 4.0 * g_yz_yz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_xy_yz[i] = 4.0 * g_yz_yz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_xy_zz[i] = 4.0 * g_yz_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1884-1890)

    #pragma omp simd aligned(g_y_z_0_0_z_y_xz_xx, g_y_z_0_0_z_y_xz_xy, g_y_z_0_0_z_y_xz_xz, g_y_z_0_0_z_y_xz_yy, g_y_z_0_0_z_y_xz_yz, g_y_z_0_0_z_y_xz_zz, g_yz_yz_xz_xx, g_yz_yz_xz_xy, g_yz_yz_xz_xz, g_yz_yz_xz_yy, g_yz_yz_xz_yz, g_yz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_y_xz_xx[i] = 4.0 * g_yz_yz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_xz_xy[i] = 4.0 * g_yz_yz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_xz_xz[i] = 4.0 * g_yz_yz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_xz_yy[i] = 4.0 * g_yz_yz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_xz_yz[i] = 4.0 * g_yz_yz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_xz_zz[i] = 4.0 * g_yz_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1890-1896)

    #pragma omp simd aligned(g_y_z_0_0_z_y_yy_xx, g_y_z_0_0_z_y_yy_xy, g_y_z_0_0_z_y_yy_xz, g_y_z_0_0_z_y_yy_yy, g_y_z_0_0_z_y_yy_yz, g_y_z_0_0_z_y_yy_zz, g_yz_yz_yy_xx, g_yz_yz_yy_xy, g_yz_yz_yy_xz, g_yz_yz_yy_yy, g_yz_yz_yy_yz, g_yz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_y_yy_xx[i] = 4.0 * g_yz_yz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_yy_xy[i] = 4.0 * g_yz_yz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_yy_xz[i] = 4.0 * g_yz_yz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_yy_yy[i] = 4.0 * g_yz_yz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_yy_yz[i] = 4.0 * g_yz_yz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_yy_zz[i] = 4.0 * g_yz_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1896-1902)

    #pragma omp simd aligned(g_y_z_0_0_z_y_yz_xx, g_y_z_0_0_z_y_yz_xy, g_y_z_0_0_z_y_yz_xz, g_y_z_0_0_z_y_yz_yy, g_y_z_0_0_z_y_yz_yz, g_y_z_0_0_z_y_yz_zz, g_yz_yz_yz_xx, g_yz_yz_yz_xy, g_yz_yz_yz_xz, g_yz_yz_yz_yy, g_yz_yz_yz_yz, g_yz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_y_yz_xx[i] = 4.0 * g_yz_yz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_yz_xy[i] = 4.0 * g_yz_yz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_yz_xz[i] = 4.0 * g_yz_yz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_yz_yy[i] = 4.0 * g_yz_yz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_yz_yz[i] = 4.0 * g_yz_yz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_yz_zz[i] = 4.0 * g_yz_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1902-1908)

    #pragma omp simd aligned(g_y_z_0_0_z_y_zz_xx, g_y_z_0_0_z_y_zz_xy, g_y_z_0_0_z_y_zz_xz, g_y_z_0_0_z_y_zz_yy, g_y_z_0_0_z_y_zz_yz, g_y_z_0_0_z_y_zz_zz, g_yz_yz_zz_xx, g_yz_yz_zz_xy, g_yz_yz_zz_xz, g_yz_yz_zz_yy, g_yz_yz_zz_yz, g_yz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_y_zz_xx[i] = 4.0 * g_yz_yz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_zz_xy[i] = 4.0 * g_yz_yz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_zz_xz[i] = 4.0 * g_yz_yz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_zz_yy[i] = 4.0 * g_yz_yz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_zz_yz[i] = 4.0 * g_yz_yz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_zz_zz[i] = 4.0 * g_yz_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1908-1914)

    #pragma omp simd aligned(g_y_z_0_0_z_z_xx_xx, g_y_z_0_0_z_z_xx_xy, g_y_z_0_0_z_z_xx_xz, g_y_z_0_0_z_z_xx_yy, g_y_z_0_0_z_z_xx_yz, g_y_z_0_0_z_z_xx_zz, g_yz_0_xx_xx, g_yz_0_xx_xy, g_yz_0_xx_xz, g_yz_0_xx_yy, g_yz_0_xx_yz, g_yz_0_xx_zz, g_yz_zz_xx_xx, g_yz_zz_xx_xy, g_yz_zz_xx_xz, g_yz_zz_xx_yy, g_yz_zz_xx_yz, g_yz_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_z_xx_xx[i] = -2.0 * g_yz_0_xx_xx[i] * a_exp + 4.0 * g_yz_zz_xx_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_xx_xy[i] = -2.0 * g_yz_0_xx_xy[i] * a_exp + 4.0 * g_yz_zz_xx_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_xx_xz[i] = -2.0 * g_yz_0_xx_xz[i] * a_exp + 4.0 * g_yz_zz_xx_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_xx_yy[i] = -2.0 * g_yz_0_xx_yy[i] * a_exp + 4.0 * g_yz_zz_xx_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_xx_yz[i] = -2.0 * g_yz_0_xx_yz[i] * a_exp + 4.0 * g_yz_zz_xx_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_xx_zz[i] = -2.0 * g_yz_0_xx_zz[i] * a_exp + 4.0 * g_yz_zz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1914-1920)

    #pragma omp simd aligned(g_y_z_0_0_z_z_xy_xx, g_y_z_0_0_z_z_xy_xy, g_y_z_0_0_z_z_xy_xz, g_y_z_0_0_z_z_xy_yy, g_y_z_0_0_z_z_xy_yz, g_y_z_0_0_z_z_xy_zz, g_yz_0_xy_xx, g_yz_0_xy_xy, g_yz_0_xy_xz, g_yz_0_xy_yy, g_yz_0_xy_yz, g_yz_0_xy_zz, g_yz_zz_xy_xx, g_yz_zz_xy_xy, g_yz_zz_xy_xz, g_yz_zz_xy_yy, g_yz_zz_xy_yz, g_yz_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_z_xy_xx[i] = -2.0 * g_yz_0_xy_xx[i] * a_exp + 4.0 * g_yz_zz_xy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_xy_xy[i] = -2.0 * g_yz_0_xy_xy[i] * a_exp + 4.0 * g_yz_zz_xy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_xy_xz[i] = -2.0 * g_yz_0_xy_xz[i] * a_exp + 4.0 * g_yz_zz_xy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_xy_yy[i] = -2.0 * g_yz_0_xy_yy[i] * a_exp + 4.0 * g_yz_zz_xy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_xy_yz[i] = -2.0 * g_yz_0_xy_yz[i] * a_exp + 4.0 * g_yz_zz_xy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_xy_zz[i] = -2.0 * g_yz_0_xy_zz[i] * a_exp + 4.0 * g_yz_zz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1920-1926)

    #pragma omp simd aligned(g_y_z_0_0_z_z_xz_xx, g_y_z_0_0_z_z_xz_xy, g_y_z_0_0_z_z_xz_xz, g_y_z_0_0_z_z_xz_yy, g_y_z_0_0_z_z_xz_yz, g_y_z_0_0_z_z_xz_zz, g_yz_0_xz_xx, g_yz_0_xz_xy, g_yz_0_xz_xz, g_yz_0_xz_yy, g_yz_0_xz_yz, g_yz_0_xz_zz, g_yz_zz_xz_xx, g_yz_zz_xz_xy, g_yz_zz_xz_xz, g_yz_zz_xz_yy, g_yz_zz_xz_yz, g_yz_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_z_xz_xx[i] = -2.0 * g_yz_0_xz_xx[i] * a_exp + 4.0 * g_yz_zz_xz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_xz_xy[i] = -2.0 * g_yz_0_xz_xy[i] * a_exp + 4.0 * g_yz_zz_xz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_xz_xz[i] = -2.0 * g_yz_0_xz_xz[i] * a_exp + 4.0 * g_yz_zz_xz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_xz_yy[i] = -2.0 * g_yz_0_xz_yy[i] * a_exp + 4.0 * g_yz_zz_xz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_xz_yz[i] = -2.0 * g_yz_0_xz_yz[i] * a_exp + 4.0 * g_yz_zz_xz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_xz_zz[i] = -2.0 * g_yz_0_xz_zz[i] * a_exp + 4.0 * g_yz_zz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1926-1932)

    #pragma omp simd aligned(g_y_z_0_0_z_z_yy_xx, g_y_z_0_0_z_z_yy_xy, g_y_z_0_0_z_z_yy_xz, g_y_z_0_0_z_z_yy_yy, g_y_z_0_0_z_z_yy_yz, g_y_z_0_0_z_z_yy_zz, g_yz_0_yy_xx, g_yz_0_yy_xy, g_yz_0_yy_xz, g_yz_0_yy_yy, g_yz_0_yy_yz, g_yz_0_yy_zz, g_yz_zz_yy_xx, g_yz_zz_yy_xy, g_yz_zz_yy_xz, g_yz_zz_yy_yy, g_yz_zz_yy_yz, g_yz_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_z_yy_xx[i] = -2.0 * g_yz_0_yy_xx[i] * a_exp + 4.0 * g_yz_zz_yy_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_yy_xy[i] = -2.0 * g_yz_0_yy_xy[i] * a_exp + 4.0 * g_yz_zz_yy_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_yy_xz[i] = -2.0 * g_yz_0_yy_xz[i] * a_exp + 4.0 * g_yz_zz_yy_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_yy_yy[i] = -2.0 * g_yz_0_yy_yy[i] * a_exp + 4.0 * g_yz_zz_yy_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_yy_yz[i] = -2.0 * g_yz_0_yy_yz[i] * a_exp + 4.0 * g_yz_zz_yy_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_yy_zz[i] = -2.0 * g_yz_0_yy_zz[i] * a_exp + 4.0 * g_yz_zz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1932-1938)

    #pragma omp simd aligned(g_y_z_0_0_z_z_yz_xx, g_y_z_0_0_z_z_yz_xy, g_y_z_0_0_z_z_yz_xz, g_y_z_0_0_z_z_yz_yy, g_y_z_0_0_z_z_yz_yz, g_y_z_0_0_z_z_yz_zz, g_yz_0_yz_xx, g_yz_0_yz_xy, g_yz_0_yz_xz, g_yz_0_yz_yy, g_yz_0_yz_yz, g_yz_0_yz_zz, g_yz_zz_yz_xx, g_yz_zz_yz_xy, g_yz_zz_yz_xz, g_yz_zz_yz_yy, g_yz_zz_yz_yz, g_yz_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_z_yz_xx[i] = -2.0 * g_yz_0_yz_xx[i] * a_exp + 4.0 * g_yz_zz_yz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_yz_xy[i] = -2.0 * g_yz_0_yz_xy[i] * a_exp + 4.0 * g_yz_zz_yz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_yz_xz[i] = -2.0 * g_yz_0_yz_xz[i] * a_exp + 4.0 * g_yz_zz_yz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_yz_yy[i] = -2.0 * g_yz_0_yz_yy[i] * a_exp + 4.0 * g_yz_zz_yz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_yz_yz[i] = -2.0 * g_yz_0_yz_yz[i] * a_exp + 4.0 * g_yz_zz_yz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_yz_zz[i] = -2.0 * g_yz_0_yz_zz[i] * a_exp + 4.0 * g_yz_zz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1938-1944)

    #pragma omp simd aligned(g_y_z_0_0_z_z_zz_xx, g_y_z_0_0_z_z_zz_xy, g_y_z_0_0_z_z_zz_xz, g_y_z_0_0_z_z_zz_yy, g_y_z_0_0_z_z_zz_yz, g_y_z_0_0_z_z_zz_zz, g_yz_0_zz_xx, g_yz_0_zz_xy, g_yz_0_zz_xz, g_yz_0_zz_yy, g_yz_0_zz_yz, g_yz_0_zz_zz, g_yz_zz_zz_xx, g_yz_zz_zz_xy, g_yz_zz_zz_xz, g_yz_zz_zz_yy, g_yz_zz_zz_yz, g_yz_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_z_zz_xx[i] = -2.0 * g_yz_0_zz_xx[i] * a_exp + 4.0 * g_yz_zz_zz_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_zz_xy[i] = -2.0 * g_yz_0_zz_xy[i] * a_exp + 4.0 * g_yz_zz_zz_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_zz_xz[i] = -2.0 * g_yz_0_zz_xz[i] * a_exp + 4.0 * g_yz_zz_zz_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_zz_yy[i] = -2.0 * g_yz_0_zz_yy[i] * a_exp + 4.0 * g_yz_zz_zz_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_zz_yz[i] = -2.0 * g_yz_0_zz_yz[i] * a_exp + 4.0 * g_yz_zz_zz_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_zz_zz[i] = -2.0 * g_yz_0_zz_zz[i] * a_exp + 4.0 * g_yz_zz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1944-1950)

    #pragma omp simd aligned(g_xz_0_xx_xx, g_xz_0_xx_xy, g_xz_0_xx_xz, g_xz_0_xx_yy, g_xz_0_xx_yz, g_xz_0_xx_zz, g_xz_xx_xx_xx, g_xz_xx_xx_xy, g_xz_xx_xx_xz, g_xz_xx_xx_yy, g_xz_xx_xx_yz, g_xz_xx_xx_zz, g_z_x_0_0_x_x_xx_xx, g_z_x_0_0_x_x_xx_xy, g_z_x_0_0_x_x_xx_xz, g_z_x_0_0_x_x_xx_yy, g_z_x_0_0_x_x_xx_yz, g_z_x_0_0_x_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_x_xx_xx[i] = -2.0 * g_xz_0_xx_xx[i] * a_exp + 4.0 * g_xz_xx_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_xx_xy[i] = -2.0 * g_xz_0_xx_xy[i] * a_exp + 4.0 * g_xz_xx_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_xx_xz[i] = -2.0 * g_xz_0_xx_xz[i] * a_exp + 4.0 * g_xz_xx_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_xx_yy[i] = -2.0 * g_xz_0_xx_yy[i] * a_exp + 4.0 * g_xz_xx_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_xx_yz[i] = -2.0 * g_xz_0_xx_yz[i] * a_exp + 4.0 * g_xz_xx_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_xx_zz[i] = -2.0 * g_xz_0_xx_zz[i] * a_exp + 4.0 * g_xz_xx_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1950-1956)

    #pragma omp simd aligned(g_xz_0_xy_xx, g_xz_0_xy_xy, g_xz_0_xy_xz, g_xz_0_xy_yy, g_xz_0_xy_yz, g_xz_0_xy_zz, g_xz_xx_xy_xx, g_xz_xx_xy_xy, g_xz_xx_xy_xz, g_xz_xx_xy_yy, g_xz_xx_xy_yz, g_xz_xx_xy_zz, g_z_x_0_0_x_x_xy_xx, g_z_x_0_0_x_x_xy_xy, g_z_x_0_0_x_x_xy_xz, g_z_x_0_0_x_x_xy_yy, g_z_x_0_0_x_x_xy_yz, g_z_x_0_0_x_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_x_xy_xx[i] = -2.0 * g_xz_0_xy_xx[i] * a_exp + 4.0 * g_xz_xx_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_xy_xy[i] = -2.0 * g_xz_0_xy_xy[i] * a_exp + 4.0 * g_xz_xx_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_xy_xz[i] = -2.0 * g_xz_0_xy_xz[i] * a_exp + 4.0 * g_xz_xx_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_xy_yy[i] = -2.0 * g_xz_0_xy_yy[i] * a_exp + 4.0 * g_xz_xx_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_xy_yz[i] = -2.0 * g_xz_0_xy_yz[i] * a_exp + 4.0 * g_xz_xx_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_xy_zz[i] = -2.0 * g_xz_0_xy_zz[i] * a_exp + 4.0 * g_xz_xx_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1956-1962)

    #pragma omp simd aligned(g_xz_0_xz_xx, g_xz_0_xz_xy, g_xz_0_xz_xz, g_xz_0_xz_yy, g_xz_0_xz_yz, g_xz_0_xz_zz, g_xz_xx_xz_xx, g_xz_xx_xz_xy, g_xz_xx_xz_xz, g_xz_xx_xz_yy, g_xz_xx_xz_yz, g_xz_xx_xz_zz, g_z_x_0_0_x_x_xz_xx, g_z_x_0_0_x_x_xz_xy, g_z_x_0_0_x_x_xz_xz, g_z_x_0_0_x_x_xz_yy, g_z_x_0_0_x_x_xz_yz, g_z_x_0_0_x_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_x_xz_xx[i] = -2.0 * g_xz_0_xz_xx[i] * a_exp + 4.0 * g_xz_xx_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_xz_xy[i] = -2.0 * g_xz_0_xz_xy[i] * a_exp + 4.0 * g_xz_xx_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_xz_xz[i] = -2.0 * g_xz_0_xz_xz[i] * a_exp + 4.0 * g_xz_xx_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_xz_yy[i] = -2.0 * g_xz_0_xz_yy[i] * a_exp + 4.0 * g_xz_xx_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_xz_yz[i] = -2.0 * g_xz_0_xz_yz[i] * a_exp + 4.0 * g_xz_xx_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_xz_zz[i] = -2.0 * g_xz_0_xz_zz[i] * a_exp + 4.0 * g_xz_xx_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1962-1968)

    #pragma omp simd aligned(g_xz_0_yy_xx, g_xz_0_yy_xy, g_xz_0_yy_xz, g_xz_0_yy_yy, g_xz_0_yy_yz, g_xz_0_yy_zz, g_xz_xx_yy_xx, g_xz_xx_yy_xy, g_xz_xx_yy_xz, g_xz_xx_yy_yy, g_xz_xx_yy_yz, g_xz_xx_yy_zz, g_z_x_0_0_x_x_yy_xx, g_z_x_0_0_x_x_yy_xy, g_z_x_0_0_x_x_yy_xz, g_z_x_0_0_x_x_yy_yy, g_z_x_0_0_x_x_yy_yz, g_z_x_0_0_x_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_x_yy_xx[i] = -2.0 * g_xz_0_yy_xx[i] * a_exp + 4.0 * g_xz_xx_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_yy_xy[i] = -2.0 * g_xz_0_yy_xy[i] * a_exp + 4.0 * g_xz_xx_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_yy_xz[i] = -2.0 * g_xz_0_yy_xz[i] * a_exp + 4.0 * g_xz_xx_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_yy_yy[i] = -2.0 * g_xz_0_yy_yy[i] * a_exp + 4.0 * g_xz_xx_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_yy_yz[i] = -2.0 * g_xz_0_yy_yz[i] * a_exp + 4.0 * g_xz_xx_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_yy_zz[i] = -2.0 * g_xz_0_yy_zz[i] * a_exp + 4.0 * g_xz_xx_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1968-1974)

    #pragma omp simd aligned(g_xz_0_yz_xx, g_xz_0_yz_xy, g_xz_0_yz_xz, g_xz_0_yz_yy, g_xz_0_yz_yz, g_xz_0_yz_zz, g_xz_xx_yz_xx, g_xz_xx_yz_xy, g_xz_xx_yz_xz, g_xz_xx_yz_yy, g_xz_xx_yz_yz, g_xz_xx_yz_zz, g_z_x_0_0_x_x_yz_xx, g_z_x_0_0_x_x_yz_xy, g_z_x_0_0_x_x_yz_xz, g_z_x_0_0_x_x_yz_yy, g_z_x_0_0_x_x_yz_yz, g_z_x_0_0_x_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_x_yz_xx[i] = -2.0 * g_xz_0_yz_xx[i] * a_exp + 4.0 * g_xz_xx_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_yz_xy[i] = -2.0 * g_xz_0_yz_xy[i] * a_exp + 4.0 * g_xz_xx_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_yz_xz[i] = -2.0 * g_xz_0_yz_xz[i] * a_exp + 4.0 * g_xz_xx_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_yz_yy[i] = -2.0 * g_xz_0_yz_yy[i] * a_exp + 4.0 * g_xz_xx_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_yz_yz[i] = -2.0 * g_xz_0_yz_yz[i] * a_exp + 4.0 * g_xz_xx_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_yz_zz[i] = -2.0 * g_xz_0_yz_zz[i] * a_exp + 4.0 * g_xz_xx_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1974-1980)

    #pragma omp simd aligned(g_xz_0_zz_xx, g_xz_0_zz_xy, g_xz_0_zz_xz, g_xz_0_zz_yy, g_xz_0_zz_yz, g_xz_0_zz_zz, g_xz_xx_zz_xx, g_xz_xx_zz_xy, g_xz_xx_zz_xz, g_xz_xx_zz_yy, g_xz_xx_zz_yz, g_xz_xx_zz_zz, g_z_x_0_0_x_x_zz_xx, g_z_x_0_0_x_x_zz_xy, g_z_x_0_0_x_x_zz_xz, g_z_x_0_0_x_x_zz_yy, g_z_x_0_0_x_x_zz_yz, g_z_x_0_0_x_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_x_zz_xx[i] = -2.0 * g_xz_0_zz_xx[i] * a_exp + 4.0 * g_xz_xx_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_zz_xy[i] = -2.0 * g_xz_0_zz_xy[i] * a_exp + 4.0 * g_xz_xx_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_zz_xz[i] = -2.0 * g_xz_0_zz_xz[i] * a_exp + 4.0 * g_xz_xx_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_zz_yy[i] = -2.0 * g_xz_0_zz_yy[i] * a_exp + 4.0 * g_xz_xx_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_zz_yz[i] = -2.0 * g_xz_0_zz_yz[i] * a_exp + 4.0 * g_xz_xx_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_zz_zz[i] = -2.0 * g_xz_0_zz_zz[i] * a_exp + 4.0 * g_xz_xx_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1980-1986)

    #pragma omp simd aligned(g_xz_xy_xx_xx, g_xz_xy_xx_xy, g_xz_xy_xx_xz, g_xz_xy_xx_yy, g_xz_xy_xx_yz, g_xz_xy_xx_zz, g_z_x_0_0_x_y_xx_xx, g_z_x_0_0_x_y_xx_xy, g_z_x_0_0_x_y_xx_xz, g_z_x_0_0_x_y_xx_yy, g_z_x_0_0_x_y_xx_yz, g_z_x_0_0_x_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_y_xx_xx[i] = 4.0 * g_xz_xy_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_xx_xy[i] = 4.0 * g_xz_xy_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_xx_xz[i] = 4.0 * g_xz_xy_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_xx_yy[i] = 4.0 * g_xz_xy_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_xx_yz[i] = 4.0 * g_xz_xy_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_xx_zz[i] = 4.0 * g_xz_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (1986-1992)

    #pragma omp simd aligned(g_xz_xy_xy_xx, g_xz_xy_xy_xy, g_xz_xy_xy_xz, g_xz_xy_xy_yy, g_xz_xy_xy_yz, g_xz_xy_xy_zz, g_z_x_0_0_x_y_xy_xx, g_z_x_0_0_x_y_xy_xy, g_z_x_0_0_x_y_xy_xz, g_z_x_0_0_x_y_xy_yy, g_z_x_0_0_x_y_xy_yz, g_z_x_0_0_x_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_y_xy_xx[i] = 4.0 * g_xz_xy_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_xy_xy[i] = 4.0 * g_xz_xy_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_xy_xz[i] = 4.0 * g_xz_xy_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_xy_yy[i] = 4.0 * g_xz_xy_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_xy_yz[i] = 4.0 * g_xz_xy_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_xy_zz[i] = 4.0 * g_xz_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (1992-1998)

    #pragma omp simd aligned(g_xz_xy_xz_xx, g_xz_xy_xz_xy, g_xz_xy_xz_xz, g_xz_xy_xz_yy, g_xz_xy_xz_yz, g_xz_xy_xz_zz, g_z_x_0_0_x_y_xz_xx, g_z_x_0_0_x_y_xz_xy, g_z_x_0_0_x_y_xz_xz, g_z_x_0_0_x_y_xz_yy, g_z_x_0_0_x_y_xz_yz, g_z_x_0_0_x_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_y_xz_xx[i] = 4.0 * g_xz_xy_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_xz_xy[i] = 4.0 * g_xz_xy_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_xz_xz[i] = 4.0 * g_xz_xy_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_xz_yy[i] = 4.0 * g_xz_xy_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_xz_yz[i] = 4.0 * g_xz_xy_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_xz_zz[i] = 4.0 * g_xz_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (1998-2004)

    #pragma omp simd aligned(g_xz_xy_yy_xx, g_xz_xy_yy_xy, g_xz_xy_yy_xz, g_xz_xy_yy_yy, g_xz_xy_yy_yz, g_xz_xy_yy_zz, g_z_x_0_0_x_y_yy_xx, g_z_x_0_0_x_y_yy_xy, g_z_x_0_0_x_y_yy_xz, g_z_x_0_0_x_y_yy_yy, g_z_x_0_0_x_y_yy_yz, g_z_x_0_0_x_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_y_yy_xx[i] = 4.0 * g_xz_xy_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_yy_xy[i] = 4.0 * g_xz_xy_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_yy_xz[i] = 4.0 * g_xz_xy_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_yy_yy[i] = 4.0 * g_xz_xy_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_yy_yz[i] = 4.0 * g_xz_xy_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_yy_zz[i] = 4.0 * g_xz_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2004-2010)

    #pragma omp simd aligned(g_xz_xy_yz_xx, g_xz_xy_yz_xy, g_xz_xy_yz_xz, g_xz_xy_yz_yy, g_xz_xy_yz_yz, g_xz_xy_yz_zz, g_z_x_0_0_x_y_yz_xx, g_z_x_0_0_x_y_yz_xy, g_z_x_0_0_x_y_yz_xz, g_z_x_0_0_x_y_yz_yy, g_z_x_0_0_x_y_yz_yz, g_z_x_0_0_x_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_y_yz_xx[i] = 4.0 * g_xz_xy_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_yz_xy[i] = 4.0 * g_xz_xy_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_yz_xz[i] = 4.0 * g_xz_xy_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_yz_yy[i] = 4.0 * g_xz_xy_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_yz_yz[i] = 4.0 * g_xz_xy_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_yz_zz[i] = 4.0 * g_xz_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2010-2016)

    #pragma omp simd aligned(g_xz_xy_zz_xx, g_xz_xy_zz_xy, g_xz_xy_zz_xz, g_xz_xy_zz_yy, g_xz_xy_zz_yz, g_xz_xy_zz_zz, g_z_x_0_0_x_y_zz_xx, g_z_x_0_0_x_y_zz_xy, g_z_x_0_0_x_y_zz_xz, g_z_x_0_0_x_y_zz_yy, g_z_x_0_0_x_y_zz_yz, g_z_x_0_0_x_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_y_zz_xx[i] = 4.0 * g_xz_xy_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_zz_xy[i] = 4.0 * g_xz_xy_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_zz_xz[i] = 4.0 * g_xz_xy_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_zz_yy[i] = 4.0 * g_xz_xy_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_zz_yz[i] = 4.0 * g_xz_xy_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_zz_zz[i] = 4.0 * g_xz_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2016-2022)

    #pragma omp simd aligned(g_xz_xz_xx_xx, g_xz_xz_xx_xy, g_xz_xz_xx_xz, g_xz_xz_xx_yy, g_xz_xz_xx_yz, g_xz_xz_xx_zz, g_z_x_0_0_x_z_xx_xx, g_z_x_0_0_x_z_xx_xy, g_z_x_0_0_x_z_xx_xz, g_z_x_0_0_x_z_xx_yy, g_z_x_0_0_x_z_xx_yz, g_z_x_0_0_x_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_z_xx_xx[i] = 4.0 * g_xz_xz_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_xx_xy[i] = 4.0 * g_xz_xz_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_xx_xz[i] = 4.0 * g_xz_xz_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_xx_yy[i] = 4.0 * g_xz_xz_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_xx_yz[i] = 4.0 * g_xz_xz_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_xx_zz[i] = 4.0 * g_xz_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2022-2028)

    #pragma omp simd aligned(g_xz_xz_xy_xx, g_xz_xz_xy_xy, g_xz_xz_xy_xz, g_xz_xz_xy_yy, g_xz_xz_xy_yz, g_xz_xz_xy_zz, g_z_x_0_0_x_z_xy_xx, g_z_x_0_0_x_z_xy_xy, g_z_x_0_0_x_z_xy_xz, g_z_x_0_0_x_z_xy_yy, g_z_x_0_0_x_z_xy_yz, g_z_x_0_0_x_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_z_xy_xx[i] = 4.0 * g_xz_xz_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_xy_xy[i] = 4.0 * g_xz_xz_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_xy_xz[i] = 4.0 * g_xz_xz_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_xy_yy[i] = 4.0 * g_xz_xz_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_xy_yz[i] = 4.0 * g_xz_xz_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_xy_zz[i] = 4.0 * g_xz_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2028-2034)

    #pragma omp simd aligned(g_xz_xz_xz_xx, g_xz_xz_xz_xy, g_xz_xz_xz_xz, g_xz_xz_xz_yy, g_xz_xz_xz_yz, g_xz_xz_xz_zz, g_z_x_0_0_x_z_xz_xx, g_z_x_0_0_x_z_xz_xy, g_z_x_0_0_x_z_xz_xz, g_z_x_0_0_x_z_xz_yy, g_z_x_0_0_x_z_xz_yz, g_z_x_0_0_x_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_z_xz_xx[i] = 4.0 * g_xz_xz_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_xz_xy[i] = 4.0 * g_xz_xz_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_xz_xz[i] = 4.0 * g_xz_xz_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_xz_yy[i] = 4.0 * g_xz_xz_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_xz_yz[i] = 4.0 * g_xz_xz_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_xz_zz[i] = 4.0 * g_xz_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2034-2040)

    #pragma omp simd aligned(g_xz_xz_yy_xx, g_xz_xz_yy_xy, g_xz_xz_yy_xz, g_xz_xz_yy_yy, g_xz_xz_yy_yz, g_xz_xz_yy_zz, g_z_x_0_0_x_z_yy_xx, g_z_x_0_0_x_z_yy_xy, g_z_x_0_0_x_z_yy_xz, g_z_x_0_0_x_z_yy_yy, g_z_x_0_0_x_z_yy_yz, g_z_x_0_0_x_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_z_yy_xx[i] = 4.0 * g_xz_xz_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_yy_xy[i] = 4.0 * g_xz_xz_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_yy_xz[i] = 4.0 * g_xz_xz_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_yy_yy[i] = 4.0 * g_xz_xz_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_yy_yz[i] = 4.0 * g_xz_xz_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_yy_zz[i] = 4.0 * g_xz_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2040-2046)

    #pragma omp simd aligned(g_xz_xz_yz_xx, g_xz_xz_yz_xy, g_xz_xz_yz_xz, g_xz_xz_yz_yy, g_xz_xz_yz_yz, g_xz_xz_yz_zz, g_z_x_0_0_x_z_yz_xx, g_z_x_0_0_x_z_yz_xy, g_z_x_0_0_x_z_yz_xz, g_z_x_0_0_x_z_yz_yy, g_z_x_0_0_x_z_yz_yz, g_z_x_0_0_x_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_z_yz_xx[i] = 4.0 * g_xz_xz_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_yz_xy[i] = 4.0 * g_xz_xz_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_yz_xz[i] = 4.0 * g_xz_xz_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_yz_yy[i] = 4.0 * g_xz_xz_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_yz_yz[i] = 4.0 * g_xz_xz_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_yz_zz[i] = 4.0 * g_xz_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2046-2052)

    #pragma omp simd aligned(g_xz_xz_zz_xx, g_xz_xz_zz_xy, g_xz_xz_zz_xz, g_xz_xz_zz_yy, g_xz_xz_zz_yz, g_xz_xz_zz_zz, g_z_x_0_0_x_z_zz_xx, g_z_x_0_0_x_z_zz_xy, g_z_x_0_0_x_z_zz_xz, g_z_x_0_0_x_z_zz_yy, g_z_x_0_0_x_z_zz_yz, g_z_x_0_0_x_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_z_zz_xx[i] = 4.0 * g_xz_xz_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_zz_xy[i] = 4.0 * g_xz_xz_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_zz_xz[i] = 4.0 * g_xz_xz_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_zz_yy[i] = 4.0 * g_xz_xz_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_zz_yz[i] = 4.0 * g_xz_xz_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_zz_zz[i] = 4.0 * g_xz_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2052-2058)

    #pragma omp simd aligned(g_yz_0_xx_xx, g_yz_0_xx_xy, g_yz_0_xx_xz, g_yz_0_xx_yy, g_yz_0_xx_yz, g_yz_0_xx_zz, g_yz_xx_xx_xx, g_yz_xx_xx_xy, g_yz_xx_xx_xz, g_yz_xx_xx_yy, g_yz_xx_xx_yz, g_yz_xx_xx_zz, g_z_x_0_0_y_x_xx_xx, g_z_x_0_0_y_x_xx_xy, g_z_x_0_0_y_x_xx_xz, g_z_x_0_0_y_x_xx_yy, g_z_x_0_0_y_x_xx_yz, g_z_x_0_0_y_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_x_xx_xx[i] = -2.0 * g_yz_0_xx_xx[i] * a_exp + 4.0 * g_yz_xx_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_xx_xy[i] = -2.0 * g_yz_0_xx_xy[i] * a_exp + 4.0 * g_yz_xx_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_xx_xz[i] = -2.0 * g_yz_0_xx_xz[i] * a_exp + 4.0 * g_yz_xx_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_xx_yy[i] = -2.0 * g_yz_0_xx_yy[i] * a_exp + 4.0 * g_yz_xx_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_xx_yz[i] = -2.0 * g_yz_0_xx_yz[i] * a_exp + 4.0 * g_yz_xx_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_xx_zz[i] = -2.0 * g_yz_0_xx_zz[i] * a_exp + 4.0 * g_yz_xx_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2058-2064)

    #pragma omp simd aligned(g_yz_0_xy_xx, g_yz_0_xy_xy, g_yz_0_xy_xz, g_yz_0_xy_yy, g_yz_0_xy_yz, g_yz_0_xy_zz, g_yz_xx_xy_xx, g_yz_xx_xy_xy, g_yz_xx_xy_xz, g_yz_xx_xy_yy, g_yz_xx_xy_yz, g_yz_xx_xy_zz, g_z_x_0_0_y_x_xy_xx, g_z_x_0_0_y_x_xy_xy, g_z_x_0_0_y_x_xy_xz, g_z_x_0_0_y_x_xy_yy, g_z_x_0_0_y_x_xy_yz, g_z_x_0_0_y_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_x_xy_xx[i] = -2.0 * g_yz_0_xy_xx[i] * a_exp + 4.0 * g_yz_xx_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_xy_xy[i] = -2.0 * g_yz_0_xy_xy[i] * a_exp + 4.0 * g_yz_xx_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_xy_xz[i] = -2.0 * g_yz_0_xy_xz[i] * a_exp + 4.0 * g_yz_xx_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_xy_yy[i] = -2.0 * g_yz_0_xy_yy[i] * a_exp + 4.0 * g_yz_xx_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_xy_yz[i] = -2.0 * g_yz_0_xy_yz[i] * a_exp + 4.0 * g_yz_xx_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_xy_zz[i] = -2.0 * g_yz_0_xy_zz[i] * a_exp + 4.0 * g_yz_xx_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2064-2070)

    #pragma omp simd aligned(g_yz_0_xz_xx, g_yz_0_xz_xy, g_yz_0_xz_xz, g_yz_0_xz_yy, g_yz_0_xz_yz, g_yz_0_xz_zz, g_yz_xx_xz_xx, g_yz_xx_xz_xy, g_yz_xx_xz_xz, g_yz_xx_xz_yy, g_yz_xx_xz_yz, g_yz_xx_xz_zz, g_z_x_0_0_y_x_xz_xx, g_z_x_0_0_y_x_xz_xy, g_z_x_0_0_y_x_xz_xz, g_z_x_0_0_y_x_xz_yy, g_z_x_0_0_y_x_xz_yz, g_z_x_0_0_y_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_x_xz_xx[i] = -2.0 * g_yz_0_xz_xx[i] * a_exp + 4.0 * g_yz_xx_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_xz_xy[i] = -2.0 * g_yz_0_xz_xy[i] * a_exp + 4.0 * g_yz_xx_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_xz_xz[i] = -2.0 * g_yz_0_xz_xz[i] * a_exp + 4.0 * g_yz_xx_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_xz_yy[i] = -2.0 * g_yz_0_xz_yy[i] * a_exp + 4.0 * g_yz_xx_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_xz_yz[i] = -2.0 * g_yz_0_xz_yz[i] * a_exp + 4.0 * g_yz_xx_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_xz_zz[i] = -2.0 * g_yz_0_xz_zz[i] * a_exp + 4.0 * g_yz_xx_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2070-2076)

    #pragma omp simd aligned(g_yz_0_yy_xx, g_yz_0_yy_xy, g_yz_0_yy_xz, g_yz_0_yy_yy, g_yz_0_yy_yz, g_yz_0_yy_zz, g_yz_xx_yy_xx, g_yz_xx_yy_xy, g_yz_xx_yy_xz, g_yz_xx_yy_yy, g_yz_xx_yy_yz, g_yz_xx_yy_zz, g_z_x_0_0_y_x_yy_xx, g_z_x_0_0_y_x_yy_xy, g_z_x_0_0_y_x_yy_xz, g_z_x_0_0_y_x_yy_yy, g_z_x_0_0_y_x_yy_yz, g_z_x_0_0_y_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_x_yy_xx[i] = -2.0 * g_yz_0_yy_xx[i] * a_exp + 4.0 * g_yz_xx_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_yy_xy[i] = -2.0 * g_yz_0_yy_xy[i] * a_exp + 4.0 * g_yz_xx_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_yy_xz[i] = -2.0 * g_yz_0_yy_xz[i] * a_exp + 4.0 * g_yz_xx_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_yy_yy[i] = -2.0 * g_yz_0_yy_yy[i] * a_exp + 4.0 * g_yz_xx_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_yy_yz[i] = -2.0 * g_yz_0_yy_yz[i] * a_exp + 4.0 * g_yz_xx_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_yy_zz[i] = -2.0 * g_yz_0_yy_zz[i] * a_exp + 4.0 * g_yz_xx_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2076-2082)

    #pragma omp simd aligned(g_yz_0_yz_xx, g_yz_0_yz_xy, g_yz_0_yz_xz, g_yz_0_yz_yy, g_yz_0_yz_yz, g_yz_0_yz_zz, g_yz_xx_yz_xx, g_yz_xx_yz_xy, g_yz_xx_yz_xz, g_yz_xx_yz_yy, g_yz_xx_yz_yz, g_yz_xx_yz_zz, g_z_x_0_0_y_x_yz_xx, g_z_x_0_0_y_x_yz_xy, g_z_x_0_0_y_x_yz_xz, g_z_x_0_0_y_x_yz_yy, g_z_x_0_0_y_x_yz_yz, g_z_x_0_0_y_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_x_yz_xx[i] = -2.0 * g_yz_0_yz_xx[i] * a_exp + 4.0 * g_yz_xx_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_yz_xy[i] = -2.0 * g_yz_0_yz_xy[i] * a_exp + 4.0 * g_yz_xx_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_yz_xz[i] = -2.0 * g_yz_0_yz_xz[i] * a_exp + 4.0 * g_yz_xx_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_yz_yy[i] = -2.0 * g_yz_0_yz_yy[i] * a_exp + 4.0 * g_yz_xx_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_yz_yz[i] = -2.0 * g_yz_0_yz_yz[i] * a_exp + 4.0 * g_yz_xx_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_yz_zz[i] = -2.0 * g_yz_0_yz_zz[i] * a_exp + 4.0 * g_yz_xx_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2082-2088)

    #pragma omp simd aligned(g_yz_0_zz_xx, g_yz_0_zz_xy, g_yz_0_zz_xz, g_yz_0_zz_yy, g_yz_0_zz_yz, g_yz_0_zz_zz, g_yz_xx_zz_xx, g_yz_xx_zz_xy, g_yz_xx_zz_xz, g_yz_xx_zz_yy, g_yz_xx_zz_yz, g_yz_xx_zz_zz, g_z_x_0_0_y_x_zz_xx, g_z_x_0_0_y_x_zz_xy, g_z_x_0_0_y_x_zz_xz, g_z_x_0_0_y_x_zz_yy, g_z_x_0_0_y_x_zz_yz, g_z_x_0_0_y_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_x_zz_xx[i] = -2.0 * g_yz_0_zz_xx[i] * a_exp + 4.0 * g_yz_xx_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_zz_xy[i] = -2.0 * g_yz_0_zz_xy[i] * a_exp + 4.0 * g_yz_xx_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_zz_xz[i] = -2.0 * g_yz_0_zz_xz[i] * a_exp + 4.0 * g_yz_xx_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_zz_yy[i] = -2.0 * g_yz_0_zz_yy[i] * a_exp + 4.0 * g_yz_xx_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_zz_yz[i] = -2.0 * g_yz_0_zz_yz[i] * a_exp + 4.0 * g_yz_xx_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_zz_zz[i] = -2.0 * g_yz_0_zz_zz[i] * a_exp + 4.0 * g_yz_xx_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2088-2094)

    #pragma omp simd aligned(g_yz_xy_xx_xx, g_yz_xy_xx_xy, g_yz_xy_xx_xz, g_yz_xy_xx_yy, g_yz_xy_xx_yz, g_yz_xy_xx_zz, g_z_x_0_0_y_y_xx_xx, g_z_x_0_0_y_y_xx_xy, g_z_x_0_0_y_y_xx_xz, g_z_x_0_0_y_y_xx_yy, g_z_x_0_0_y_y_xx_yz, g_z_x_0_0_y_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_y_xx_xx[i] = 4.0 * g_yz_xy_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_xx_xy[i] = 4.0 * g_yz_xy_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_xx_xz[i] = 4.0 * g_yz_xy_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_xx_yy[i] = 4.0 * g_yz_xy_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_xx_yz[i] = 4.0 * g_yz_xy_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_xx_zz[i] = 4.0 * g_yz_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2094-2100)

    #pragma omp simd aligned(g_yz_xy_xy_xx, g_yz_xy_xy_xy, g_yz_xy_xy_xz, g_yz_xy_xy_yy, g_yz_xy_xy_yz, g_yz_xy_xy_zz, g_z_x_0_0_y_y_xy_xx, g_z_x_0_0_y_y_xy_xy, g_z_x_0_0_y_y_xy_xz, g_z_x_0_0_y_y_xy_yy, g_z_x_0_0_y_y_xy_yz, g_z_x_0_0_y_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_y_xy_xx[i] = 4.0 * g_yz_xy_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_xy_xy[i] = 4.0 * g_yz_xy_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_xy_xz[i] = 4.0 * g_yz_xy_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_xy_yy[i] = 4.0 * g_yz_xy_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_xy_yz[i] = 4.0 * g_yz_xy_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_xy_zz[i] = 4.0 * g_yz_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2100-2106)

    #pragma omp simd aligned(g_yz_xy_xz_xx, g_yz_xy_xz_xy, g_yz_xy_xz_xz, g_yz_xy_xz_yy, g_yz_xy_xz_yz, g_yz_xy_xz_zz, g_z_x_0_0_y_y_xz_xx, g_z_x_0_0_y_y_xz_xy, g_z_x_0_0_y_y_xz_xz, g_z_x_0_0_y_y_xz_yy, g_z_x_0_0_y_y_xz_yz, g_z_x_0_0_y_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_y_xz_xx[i] = 4.0 * g_yz_xy_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_xz_xy[i] = 4.0 * g_yz_xy_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_xz_xz[i] = 4.0 * g_yz_xy_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_xz_yy[i] = 4.0 * g_yz_xy_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_xz_yz[i] = 4.0 * g_yz_xy_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_xz_zz[i] = 4.0 * g_yz_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2106-2112)

    #pragma omp simd aligned(g_yz_xy_yy_xx, g_yz_xy_yy_xy, g_yz_xy_yy_xz, g_yz_xy_yy_yy, g_yz_xy_yy_yz, g_yz_xy_yy_zz, g_z_x_0_0_y_y_yy_xx, g_z_x_0_0_y_y_yy_xy, g_z_x_0_0_y_y_yy_xz, g_z_x_0_0_y_y_yy_yy, g_z_x_0_0_y_y_yy_yz, g_z_x_0_0_y_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_y_yy_xx[i] = 4.0 * g_yz_xy_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_yy_xy[i] = 4.0 * g_yz_xy_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_yy_xz[i] = 4.0 * g_yz_xy_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_yy_yy[i] = 4.0 * g_yz_xy_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_yy_yz[i] = 4.0 * g_yz_xy_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_yy_zz[i] = 4.0 * g_yz_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2112-2118)

    #pragma omp simd aligned(g_yz_xy_yz_xx, g_yz_xy_yz_xy, g_yz_xy_yz_xz, g_yz_xy_yz_yy, g_yz_xy_yz_yz, g_yz_xy_yz_zz, g_z_x_0_0_y_y_yz_xx, g_z_x_0_0_y_y_yz_xy, g_z_x_0_0_y_y_yz_xz, g_z_x_0_0_y_y_yz_yy, g_z_x_0_0_y_y_yz_yz, g_z_x_0_0_y_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_y_yz_xx[i] = 4.0 * g_yz_xy_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_yz_xy[i] = 4.0 * g_yz_xy_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_yz_xz[i] = 4.0 * g_yz_xy_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_yz_yy[i] = 4.0 * g_yz_xy_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_yz_yz[i] = 4.0 * g_yz_xy_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_yz_zz[i] = 4.0 * g_yz_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2118-2124)

    #pragma omp simd aligned(g_yz_xy_zz_xx, g_yz_xy_zz_xy, g_yz_xy_zz_xz, g_yz_xy_zz_yy, g_yz_xy_zz_yz, g_yz_xy_zz_zz, g_z_x_0_0_y_y_zz_xx, g_z_x_0_0_y_y_zz_xy, g_z_x_0_0_y_y_zz_xz, g_z_x_0_0_y_y_zz_yy, g_z_x_0_0_y_y_zz_yz, g_z_x_0_0_y_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_y_zz_xx[i] = 4.0 * g_yz_xy_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_zz_xy[i] = 4.0 * g_yz_xy_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_zz_xz[i] = 4.0 * g_yz_xy_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_zz_yy[i] = 4.0 * g_yz_xy_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_zz_yz[i] = 4.0 * g_yz_xy_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_zz_zz[i] = 4.0 * g_yz_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2124-2130)

    #pragma omp simd aligned(g_yz_xz_xx_xx, g_yz_xz_xx_xy, g_yz_xz_xx_xz, g_yz_xz_xx_yy, g_yz_xz_xx_yz, g_yz_xz_xx_zz, g_z_x_0_0_y_z_xx_xx, g_z_x_0_0_y_z_xx_xy, g_z_x_0_0_y_z_xx_xz, g_z_x_0_0_y_z_xx_yy, g_z_x_0_0_y_z_xx_yz, g_z_x_0_0_y_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_z_xx_xx[i] = 4.0 * g_yz_xz_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_xx_xy[i] = 4.0 * g_yz_xz_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_xx_xz[i] = 4.0 * g_yz_xz_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_xx_yy[i] = 4.0 * g_yz_xz_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_xx_yz[i] = 4.0 * g_yz_xz_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_xx_zz[i] = 4.0 * g_yz_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2130-2136)

    #pragma omp simd aligned(g_yz_xz_xy_xx, g_yz_xz_xy_xy, g_yz_xz_xy_xz, g_yz_xz_xy_yy, g_yz_xz_xy_yz, g_yz_xz_xy_zz, g_z_x_0_0_y_z_xy_xx, g_z_x_0_0_y_z_xy_xy, g_z_x_0_0_y_z_xy_xz, g_z_x_0_0_y_z_xy_yy, g_z_x_0_0_y_z_xy_yz, g_z_x_0_0_y_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_z_xy_xx[i] = 4.0 * g_yz_xz_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_xy_xy[i] = 4.0 * g_yz_xz_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_xy_xz[i] = 4.0 * g_yz_xz_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_xy_yy[i] = 4.0 * g_yz_xz_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_xy_yz[i] = 4.0 * g_yz_xz_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_xy_zz[i] = 4.0 * g_yz_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2136-2142)

    #pragma omp simd aligned(g_yz_xz_xz_xx, g_yz_xz_xz_xy, g_yz_xz_xz_xz, g_yz_xz_xz_yy, g_yz_xz_xz_yz, g_yz_xz_xz_zz, g_z_x_0_0_y_z_xz_xx, g_z_x_0_0_y_z_xz_xy, g_z_x_0_0_y_z_xz_xz, g_z_x_0_0_y_z_xz_yy, g_z_x_0_0_y_z_xz_yz, g_z_x_0_0_y_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_z_xz_xx[i] = 4.0 * g_yz_xz_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_xz_xy[i] = 4.0 * g_yz_xz_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_xz_xz[i] = 4.0 * g_yz_xz_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_xz_yy[i] = 4.0 * g_yz_xz_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_xz_yz[i] = 4.0 * g_yz_xz_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_xz_zz[i] = 4.0 * g_yz_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2142-2148)

    #pragma omp simd aligned(g_yz_xz_yy_xx, g_yz_xz_yy_xy, g_yz_xz_yy_xz, g_yz_xz_yy_yy, g_yz_xz_yy_yz, g_yz_xz_yy_zz, g_z_x_0_0_y_z_yy_xx, g_z_x_0_0_y_z_yy_xy, g_z_x_0_0_y_z_yy_xz, g_z_x_0_0_y_z_yy_yy, g_z_x_0_0_y_z_yy_yz, g_z_x_0_0_y_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_z_yy_xx[i] = 4.0 * g_yz_xz_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_yy_xy[i] = 4.0 * g_yz_xz_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_yy_xz[i] = 4.0 * g_yz_xz_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_yy_yy[i] = 4.0 * g_yz_xz_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_yy_yz[i] = 4.0 * g_yz_xz_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_yy_zz[i] = 4.0 * g_yz_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2148-2154)

    #pragma omp simd aligned(g_yz_xz_yz_xx, g_yz_xz_yz_xy, g_yz_xz_yz_xz, g_yz_xz_yz_yy, g_yz_xz_yz_yz, g_yz_xz_yz_zz, g_z_x_0_0_y_z_yz_xx, g_z_x_0_0_y_z_yz_xy, g_z_x_0_0_y_z_yz_xz, g_z_x_0_0_y_z_yz_yy, g_z_x_0_0_y_z_yz_yz, g_z_x_0_0_y_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_z_yz_xx[i] = 4.0 * g_yz_xz_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_yz_xy[i] = 4.0 * g_yz_xz_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_yz_xz[i] = 4.0 * g_yz_xz_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_yz_yy[i] = 4.0 * g_yz_xz_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_yz_yz[i] = 4.0 * g_yz_xz_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_yz_zz[i] = 4.0 * g_yz_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2154-2160)

    #pragma omp simd aligned(g_yz_xz_zz_xx, g_yz_xz_zz_xy, g_yz_xz_zz_xz, g_yz_xz_zz_yy, g_yz_xz_zz_yz, g_yz_xz_zz_zz, g_z_x_0_0_y_z_zz_xx, g_z_x_0_0_y_z_zz_xy, g_z_x_0_0_y_z_zz_xz, g_z_x_0_0_y_z_zz_yy, g_z_x_0_0_y_z_zz_yz, g_z_x_0_0_y_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_z_zz_xx[i] = 4.0 * g_yz_xz_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_zz_xy[i] = 4.0 * g_yz_xz_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_zz_xz[i] = 4.0 * g_yz_xz_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_zz_yy[i] = 4.0 * g_yz_xz_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_zz_yz[i] = 4.0 * g_yz_xz_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_zz_zz[i] = 4.0 * g_yz_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2160-2166)

    #pragma omp simd aligned(g_0_0_xx_xx, g_0_0_xx_xy, g_0_0_xx_xz, g_0_0_xx_yy, g_0_0_xx_yz, g_0_0_xx_zz, g_0_xx_xx_xx, g_0_xx_xx_xy, g_0_xx_xx_xz, g_0_xx_xx_yy, g_0_xx_xx_yz, g_0_xx_xx_zz, g_z_x_0_0_z_x_xx_xx, g_z_x_0_0_z_x_xx_xy, g_z_x_0_0_z_x_xx_xz, g_z_x_0_0_z_x_xx_yy, g_z_x_0_0_z_x_xx_yz, g_z_x_0_0_z_x_xx_zz, g_zz_0_xx_xx, g_zz_0_xx_xy, g_zz_0_xx_xz, g_zz_0_xx_yy, g_zz_0_xx_yz, g_zz_0_xx_zz, g_zz_xx_xx_xx, g_zz_xx_xx_xy, g_zz_xx_xx_xz, g_zz_xx_xx_yy, g_zz_xx_xx_yz, g_zz_xx_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_x_xx_xx[i] = g_0_0_xx_xx[i] - 2.0 * g_0_xx_xx_xx[i] * b_exp - 2.0 * g_zz_0_xx_xx[i] * a_exp + 4.0 * g_zz_xx_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_xx_xy[i] = g_0_0_xx_xy[i] - 2.0 * g_0_xx_xx_xy[i] * b_exp - 2.0 * g_zz_0_xx_xy[i] * a_exp + 4.0 * g_zz_xx_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_xx_xz[i] = g_0_0_xx_xz[i] - 2.0 * g_0_xx_xx_xz[i] * b_exp - 2.0 * g_zz_0_xx_xz[i] * a_exp + 4.0 * g_zz_xx_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_xx_yy[i] = g_0_0_xx_yy[i] - 2.0 * g_0_xx_xx_yy[i] * b_exp - 2.0 * g_zz_0_xx_yy[i] * a_exp + 4.0 * g_zz_xx_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_xx_yz[i] = g_0_0_xx_yz[i] - 2.0 * g_0_xx_xx_yz[i] * b_exp - 2.0 * g_zz_0_xx_yz[i] * a_exp + 4.0 * g_zz_xx_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_xx_zz[i] = g_0_0_xx_zz[i] - 2.0 * g_0_xx_xx_zz[i] * b_exp - 2.0 * g_zz_0_xx_zz[i] * a_exp + 4.0 * g_zz_xx_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2166-2172)

    #pragma omp simd aligned(g_0_0_xy_xx, g_0_0_xy_xy, g_0_0_xy_xz, g_0_0_xy_yy, g_0_0_xy_yz, g_0_0_xy_zz, g_0_xx_xy_xx, g_0_xx_xy_xy, g_0_xx_xy_xz, g_0_xx_xy_yy, g_0_xx_xy_yz, g_0_xx_xy_zz, g_z_x_0_0_z_x_xy_xx, g_z_x_0_0_z_x_xy_xy, g_z_x_0_0_z_x_xy_xz, g_z_x_0_0_z_x_xy_yy, g_z_x_0_0_z_x_xy_yz, g_z_x_0_0_z_x_xy_zz, g_zz_0_xy_xx, g_zz_0_xy_xy, g_zz_0_xy_xz, g_zz_0_xy_yy, g_zz_0_xy_yz, g_zz_0_xy_zz, g_zz_xx_xy_xx, g_zz_xx_xy_xy, g_zz_xx_xy_xz, g_zz_xx_xy_yy, g_zz_xx_xy_yz, g_zz_xx_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_x_xy_xx[i] = g_0_0_xy_xx[i] - 2.0 * g_0_xx_xy_xx[i] * b_exp - 2.0 * g_zz_0_xy_xx[i] * a_exp + 4.0 * g_zz_xx_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_xy_xy[i] = g_0_0_xy_xy[i] - 2.0 * g_0_xx_xy_xy[i] * b_exp - 2.0 * g_zz_0_xy_xy[i] * a_exp + 4.0 * g_zz_xx_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_xy_xz[i] = g_0_0_xy_xz[i] - 2.0 * g_0_xx_xy_xz[i] * b_exp - 2.0 * g_zz_0_xy_xz[i] * a_exp + 4.0 * g_zz_xx_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_xy_yy[i] = g_0_0_xy_yy[i] - 2.0 * g_0_xx_xy_yy[i] * b_exp - 2.0 * g_zz_0_xy_yy[i] * a_exp + 4.0 * g_zz_xx_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_xy_yz[i] = g_0_0_xy_yz[i] - 2.0 * g_0_xx_xy_yz[i] * b_exp - 2.0 * g_zz_0_xy_yz[i] * a_exp + 4.0 * g_zz_xx_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_xy_zz[i] = g_0_0_xy_zz[i] - 2.0 * g_0_xx_xy_zz[i] * b_exp - 2.0 * g_zz_0_xy_zz[i] * a_exp + 4.0 * g_zz_xx_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2172-2178)

    #pragma omp simd aligned(g_0_0_xz_xx, g_0_0_xz_xy, g_0_0_xz_xz, g_0_0_xz_yy, g_0_0_xz_yz, g_0_0_xz_zz, g_0_xx_xz_xx, g_0_xx_xz_xy, g_0_xx_xz_xz, g_0_xx_xz_yy, g_0_xx_xz_yz, g_0_xx_xz_zz, g_z_x_0_0_z_x_xz_xx, g_z_x_0_0_z_x_xz_xy, g_z_x_0_0_z_x_xz_xz, g_z_x_0_0_z_x_xz_yy, g_z_x_0_0_z_x_xz_yz, g_z_x_0_0_z_x_xz_zz, g_zz_0_xz_xx, g_zz_0_xz_xy, g_zz_0_xz_xz, g_zz_0_xz_yy, g_zz_0_xz_yz, g_zz_0_xz_zz, g_zz_xx_xz_xx, g_zz_xx_xz_xy, g_zz_xx_xz_xz, g_zz_xx_xz_yy, g_zz_xx_xz_yz, g_zz_xx_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_x_xz_xx[i] = g_0_0_xz_xx[i] - 2.0 * g_0_xx_xz_xx[i] * b_exp - 2.0 * g_zz_0_xz_xx[i] * a_exp + 4.0 * g_zz_xx_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_xz_xy[i] = g_0_0_xz_xy[i] - 2.0 * g_0_xx_xz_xy[i] * b_exp - 2.0 * g_zz_0_xz_xy[i] * a_exp + 4.0 * g_zz_xx_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_xz_xz[i] = g_0_0_xz_xz[i] - 2.0 * g_0_xx_xz_xz[i] * b_exp - 2.0 * g_zz_0_xz_xz[i] * a_exp + 4.0 * g_zz_xx_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_xz_yy[i] = g_0_0_xz_yy[i] - 2.0 * g_0_xx_xz_yy[i] * b_exp - 2.0 * g_zz_0_xz_yy[i] * a_exp + 4.0 * g_zz_xx_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_xz_yz[i] = g_0_0_xz_yz[i] - 2.0 * g_0_xx_xz_yz[i] * b_exp - 2.0 * g_zz_0_xz_yz[i] * a_exp + 4.0 * g_zz_xx_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_xz_zz[i] = g_0_0_xz_zz[i] - 2.0 * g_0_xx_xz_zz[i] * b_exp - 2.0 * g_zz_0_xz_zz[i] * a_exp + 4.0 * g_zz_xx_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2178-2184)

    #pragma omp simd aligned(g_0_0_yy_xx, g_0_0_yy_xy, g_0_0_yy_xz, g_0_0_yy_yy, g_0_0_yy_yz, g_0_0_yy_zz, g_0_xx_yy_xx, g_0_xx_yy_xy, g_0_xx_yy_xz, g_0_xx_yy_yy, g_0_xx_yy_yz, g_0_xx_yy_zz, g_z_x_0_0_z_x_yy_xx, g_z_x_0_0_z_x_yy_xy, g_z_x_0_0_z_x_yy_xz, g_z_x_0_0_z_x_yy_yy, g_z_x_0_0_z_x_yy_yz, g_z_x_0_0_z_x_yy_zz, g_zz_0_yy_xx, g_zz_0_yy_xy, g_zz_0_yy_xz, g_zz_0_yy_yy, g_zz_0_yy_yz, g_zz_0_yy_zz, g_zz_xx_yy_xx, g_zz_xx_yy_xy, g_zz_xx_yy_xz, g_zz_xx_yy_yy, g_zz_xx_yy_yz, g_zz_xx_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_x_yy_xx[i] = g_0_0_yy_xx[i] - 2.0 * g_0_xx_yy_xx[i] * b_exp - 2.0 * g_zz_0_yy_xx[i] * a_exp + 4.0 * g_zz_xx_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_yy_xy[i] = g_0_0_yy_xy[i] - 2.0 * g_0_xx_yy_xy[i] * b_exp - 2.0 * g_zz_0_yy_xy[i] * a_exp + 4.0 * g_zz_xx_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_yy_xz[i] = g_0_0_yy_xz[i] - 2.0 * g_0_xx_yy_xz[i] * b_exp - 2.0 * g_zz_0_yy_xz[i] * a_exp + 4.0 * g_zz_xx_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_yy_yy[i] = g_0_0_yy_yy[i] - 2.0 * g_0_xx_yy_yy[i] * b_exp - 2.0 * g_zz_0_yy_yy[i] * a_exp + 4.0 * g_zz_xx_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_yy_yz[i] = g_0_0_yy_yz[i] - 2.0 * g_0_xx_yy_yz[i] * b_exp - 2.0 * g_zz_0_yy_yz[i] * a_exp + 4.0 * g_zz_xx_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_yy_zz[i] = g_0_0_yy_zz[i] - 2.0 * g_0_xx_yy_zz[i] * b_exp - 2.0 * g_zz_0_yy_zz[i] * a_exp + 4.0 * g_zz_xx_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2184-2190)

    #pragma omp simd aligned(g_0_0_yz_xx, g_0_0_yz_xy, g_0_0_yz_xz, g_0_0_yz_yy, g_0_0_yz_yz, g_0_0_yz_zz, g_0_xx_yz_xx, g_0_xx_yz_xy, g_0_xx_yz_xz, g_0_xx_yz_yy, g_0_xx_yz_yz, g_0_xx_yz_zz, g_z_x_0_0_z_x_yz_xx, g_z_x_0_0_z_x_yz_xy, g_z_x_0_0_z_x_yz_xz, g_z_x_0_0_z_x_yz_yy, g_z_x_0_0_z_x_yz_yz, g_z_x_0_0_z_x_yz_zz, g_zz_0_yz_xx, g_zz_0_yz_xy, g_zz_0_yz_xz, g_zz_0_yz_yy, g_zz_0_yz_yz, g_zz_0_yz_zz, g_zz_xx_yz_xx, g_zz_xx_yz_xy, g_zz_xx_yz_xz, g_zz_xx_yz_yy, g_zz_xx_yz_yz, g_zz_xx_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_x_yz_xx[i] = g_0_0_yz_xx[i] - 2.0 * g_0_xx_yz_xx[i] * b_exp - 2.0 * g_zz_0_yz_xx[i] * a_exp + 4.0 * g_zz_xx_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_yz_xy[i] = g_0_0_yz_xy[i] - 2.0 * g_0_xx_yz_xy[i] * b_exp - 2.0 * g_zz_0_yz_xy[i] * a_exp + 4.0 * g_zz_xx_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_yz_xz[i] = g_0_0_yz_xz[i] - 2.0 * g_0_xx_yz_xz[i] * b_exp - 2.0 * g_zz_0_yz_xz[i] * a_exp + 4.0 * g_zz_xx_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_yz_yy[i] = g_0_0_yz_yy[i] - 2.0 * g_0_xx_yz_yy[i] * b_exp - 2.0 * g_zz_0_yz_yy[i] * a_exp + 4.0 * g_zz_xx_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_yz_yz[i] = g_0_0_yz_yz[i] - 2.0 * g_0_xx_yz_yz[i] * b_exp - 2.0 * g_zz_0_yz_yz[i] * a_exp + 4.0 * g_zz_xx_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_yz_zz[i] = g_0_0_yz_zz[i] - 2.0 * g_0_xx_yz_zz[i] * b_exp - 2.0 * g_zz_0_yz_zz[i] * a_exp + 4.0 * g_zz_xx_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2190-2196)

    #pragma omp simd aligned(g_0_0_zz_xx, g_0_0_zz_xy, g_0_0_zz_xz, g_0_0_zz_yy, g_0_0_zz_yz, g_0_0_zz_zz, g_0_xx_zz_xx, g_0_xx_zz_xy, g_0_xx_zz_xz, g_0_xx_zz_yy, g_0_xx_zz_yz, g_0_xx_zz_zz, g_z_x_0_0_z_x_zz_xx, g_z_x_0_0_z_x_zz_xy, g_z_x_0_0_z_x_zz_xz, g_z_x_0_0_z_x_zz_yy, g_z_x_0_0_z_x_zz_yz, g_z_x_0_0_z_x_zz_zz, g_zz_0_zz_xx, g_zz_0_zz_xy, g_zz_0_zz_xz, g_zz_0_zz_yy, g_zz_0_zz_yz, g_zz_0_zz_zz, g_zz_xx_zz_xx, g_zz_xx_zz_xy, g_zz_xx_zz_xz, g_zz_xx_zz_yy, g_zz_xx_zz_yz, g_zz_xx_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_x_zz_xx[i] = g_0_0_zz_xx[i] - 2.0 * g_0_xx_zz_xx[i] * b_exp - 2.0 * g_zz_0_zz_xx[i] * a_exp + 4.0 * g_zz_xx_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_zz_xy[i] = g_0_0_zz_xy[i] - 2.0 * g_0_xx_zz_xy[i] * b_exp - 2.0 * g_zz_0_zz_xy[i] * a_exp + 4.0 * g_zz_xx_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_zz_xz[i] = g_0_0_zz_xz[i] - 2.0 * g_0_xx_zz_xz[i] * b_exp - 2.0 * g_zz_0_zz_xz[i] * a_exp + 4.0 * g_zz_xx_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_zz_yy[i] = g_0_0_zz_yy[i] - 2.0 * g_0_xx_zz_yy[i] * b_exp - 2.0 * g_zz_0_zz_yy[i] * a_exp + 4.0 * g_zz_xx_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_zz_yz[i] = g_0_0_zz_yz[i] - 2.0 * g_0_xx_zz_yz[i] * b_exp - 2.0 * g_zz_0_zz_yz[i] * a_exp + 4.0 * g_zz_xx_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_zz_zz[i] = g_0_0_zz_zz[i] - 2.0 * g_0_xx_zz_zz[i] * b_exp - 2.0 * g_zz_0_zz_zz[i] * a_exp + 4.0 * g_zz_xx_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2196-2202)

    #pragma omp simd aligned(g_0_xy_xx_xx, g_0_xy_xx_xy, g_0_xy_xx_xz, g_0_xy_xx_yy, g_0_xy_xx_yz, g_0_xy_xx_zz, g_z_x_0_0_z_y_xx_xx, g_z_x_0_0_z_y_xx_xy, g_z_x_0_0_z_y_xx_xz, g_z_x_0_0_z_y_xx_yy, g_z_x_0_0_z_y_xx_yz, g_z_x_0_0_z_y_xx_zz, g_zz_xy_xx_xx, g_zz_xy_xx_xy, g_zz_xy_xx_xz, g_zz_xy_xx_yy, g_zz_xy_xx_yz, g_zz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_y_xx_xx[i] = -2.0 * g_0_xy_xx_xx[i] * b_exp + 4.0 * g_zz_xy_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_xx_xy[i] = -2.0 * g_0_xy_xx_xy[i] * b_exp + 4.0 * g_zz_xy_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_xx_xz[i] = -2.0 * g_0_xy_xx_xz[i] * b_exp + 4.0 * g_zz_xy_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_xx_yy[i] = -2.0 * g_0_xy_xx_yy[i] * b_exp + 4.0 * g_zz_xy_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_xx_yz[i] = -2.0 * g_0_xy_xx_yz[i] * b_exp + 4.0 * g_zz_xy_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_xx_zz[i] = -2.0 * g_0_xy_xx_zz[i] * b_exp + 4.0 * g_zz_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2202-2208)

    #pragma omp simd aligned(g_0_xy_xy_xx, g_0_xy_xy_xy, g_0_xy_xy_xz, g_0_xy_xy_yy, g_0_xy_xy_yz, g_0_xy_xy_zz, g_z_x_0_0_z_y_xy_xx, g_z_x_0_0_z_y_xy_xy, g_z_x_0_0_z_y_xy_xz, g_z_x_0_0_z_y_xy_yy, g_z_x_0_0_z_y_xy_yz, g_z_x_0_0_z_y_xy_zz, g_zz_xy_xy_xx, g_zz_xy_xy_xy, g_zz_xy_xy_xz, g_zz_xy_xy_yy, g_zz_xy_xy_yz, g_zz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_y_xy_xx[i] = -2.0 * g_0_xy_xy_xx[i] * b_exp + 4.0 * g_zz_xy_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_xy_xy[i] = -2.0 * g_0_xy_xy_xy[i] * b_exp + 4.0 * g_zz_xy_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_xy_xz[i] = -2.0 * g_0_xy_xy_xz[i] * b_exp + 4.0 * g_zz_xy_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_xy_yy[i] = -2.0 * g_0_xy_xy_yy[i] * b_exp + 4.0 * g_zz_xy_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_xy_yz[i] = -2.0 * g_0_xy_xy_yz[i] * b_exp + 4.0 * g_zz_xy_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_xy_zz[i] = -2.0 * g_0_xy_xy_zz[i] * b_exp + 4.0 * g_zz_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2208-2214)

    #pragma omp simd aligned(g_0_xy_xz_xx, g_0_xy_xz_xy, g_0_xy_xz_xz, g_0_xy_xz_yy, g_0_xy_xz_yz, g_0_xy_xz_zz, g_z_x_0_0_z_y_xz_xx, g_z_x_0_0_z_y_xz_xy, g_z_x_0_0_z_y_xz_xz, g_z_x_0_0_z_y_xz_yy, g_z_x_0_0_z_y_xz_yz, g_z_x_0_0_z_y_xz_zz, g_zz_xy_xz_xx, g_zz_xy_xz_xy, g_zz_xy_xz_xz, g_zz_xy_xz_yy, g_zz_xy_xz_yz, g_zz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_y_xz_xx[i] = -2.0 * g_0_xy_xz_xx[i] * b_exp + 4.0 * g_zz_xy_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_xz_xy[i] = -2.0 * g_0_xy_xz_xy[i] * b_exp + 4.0 * g_zz_xy_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_xz_xz[i] = -2.0 * g_0_xy_xz_xz[i] * b_exp + 4.0 * g_zz_xy_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_xz_yy[i] = -2.0 * g_0_xy_xz_yy[i] * b_exp + 4.0 * g_zz_xy_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_xz_yz[i] = -2.0 * g_0_xy_xz_yz[i] * b_exp + 4.0 * g_zz_xy_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_xz_zz[i] = -2.0 * g_0_xy_xz_zz[i] * b_exp + 4.0 * g_zz_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2214-2220)

    #pragma omp simd aligned(g_0_xy_yy_xx, g_0_xy_yy_xy, g_0_xy_yy_xz, g_0_xy_yy_yy, g_0_xy_yy_yz, g_0_xy_yy_zz, g_z_x_0_0_z_y_yy_xx, g_z_x_0_0_z_y_yy_xy, g_z_x_0_0_z_y_yy_xz, g_z_x_0_0_z_y_yy_yy, g_z_x_0_0_z_y_yy_yz, g_z_x_0_0_z_y_yy_zz, g_zz_xy_yy_xx, g_zz_xy_yy_xy, g_zz_xy_yy_xz, g_zz_xy_yy_yy, g_zz_xy_yy_yz, g_zz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_y_yy_xx[i] = -2.0 * g_0_xy_yy_xx[i] * b_exp + 4.0 * g_zz_xy_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_yy_xy[i] = -2.0 * g_0_xy_yy_xy[i] * b_exp + 4.0 * g_zz_xy_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_yy_xz[i] = -2.0 * g_0_xy_yy_xz[i] * b_exp + 4.0 * g_zz_xy_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_yy_yy[i] = -2.0 * g_0_xy_yy_yy[i] * b_exp + 4.0 * g_zz_xy_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_yy_yz[i] = -2.0 * g_0_xy_yy_yz[i] * b_exp + 4.0 * g_zz_xy_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_yy_zz[i] = -2.0 * g_0_xy_yy_zz[i] * b_exp + 4.0 * g_zz_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2220-2226)

    #pragma omp simd aligned(g_0_xy_yz_xx, g_0_xy_yz_xy, g_0_xy_yz_xz, g_0_xy_yz_yy, g_0_xy_yz_yz, g_0_xy_yz_zz, g_z_x_0_0_z_y_yz_xx, g_z_x_0_0_z_y_yz_xy, g_z_x_0_0_z_y_yz_xz, g_z_x_0_0_z_y_yz_yy, g_z_x_0_0_z_y_yz_yz, g_z_x_0_0_z_y_yz_zz, g_zz_xy_yz_xx, g_zz_xy_yz_xy, g_zz_xy_yz_xz, g_zz_xy_yz_yy, g_zz_xy_yz_yz, g_zz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_y_yz_xx[i] = -2.0 * g_0_xy_yz_xx[i] * b_exp + 4.0 * g_zz_xy_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_yz_xy[i] = -2.0 * g_0_xy_yz_xy[i] * b_exp + 4.0 * g_zz_xy_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_yz_xz[i] = -2.0 * g_0_xy_yz_xz[i] * b_exp + 4.0 * g_zz_xy_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_yz_yy[i] = -2.0 * g_0_xy_yz_yy[i] * b_exp + 4.0 * g_zz_xy_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_yz_yz[i] = -2.0 * g_0_xy_yz_yz[i] * b_exp + 4.0 * g_zz_xy_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_yz_zz[i] = -2.0 * g_0_xy_yz_zz[i] * b_exp + 4.0 * g_zz_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2226-2232)

    #pragma omp simd aligned(g_0_xy_zz_xx, g_0_xy_zz_xy, g_0_xy_zz_xz, g_0_xy_zz_yy, g_0_xy_zz_yz, g_0_xy_zz_zz, g_z_x_0_0_z_y_zz_xx, g_z_x_0_0_z_y_zz_xy, g_z_x_0_0_z_y_zz_xz, g_z_x_0_0_z_y_zz_yy, g_z_x_0_0_z_y_zz_yz, g_z_x_0_0_z_y_zz_zz, g_zz_xy_zz_xx, g_zz_xy_zz_xy, g_zz_xy_zz_xz, g_zz_xy_zz_yy, g_zz_xy_zz_yz, g_zz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_y_zz_xx[i] = -2.0 * g_0_xy_zz_xx[i] * b_exp + 4.0 * g_zz_xy_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_zz_xy[i] = -2.0 * g_0_xy_zz_xy[i] * b_exp + 4.0 * g_zz_xy_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_zz_xz[i] = -2.0 * g_0_xy_zz_xz[i] * b_exp + 4.0 * g_zz_xy_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_zz_yy[i] = -2.0 * g_0_xy_zz_yy[i] * b_exp + 4.0 * g_zz_xy_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_zz_yz[i] = -2.0 * g_0_xy_zz_yz[i] * b_exp + 4.0 * g_zz_xy_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_zz_zz[i] = -2.0 * g_0_xy_zz_zz[i] * b_exp + 4.0 * g_zz_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2232-2238)

    #pragma omp simd aligned(g_0_xz_xx_xx, g_0_xz_xx_xy, g_0_xz_xx_xz, g_0_xz_xx_yy, g_0_xz_xx_yz, g_0_xz_xx_zz, g_z_x_0_0_z_z_xx_xx, g_z_x_0_0_z_z_xx_xy, g_z_x_0_0_z_z_xx_xz, g_z_x_0_0_z_z_xx_yy, g_z_x_0_0_z_z_xx_yz, g_z_x_0_0_z_z_xx_zz, g_zz_xz_xx_xx, g_zz_xz_xx_xy, g_zz_xz_xx_xz, g_zz_xz_xx_yy, g_zz_xz_xx_yz, g_zz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_z_xx_xx[i] = -2.0 * g_0_xz_xx_xx[i] * b_exp + 4.0 * g_zz_xz_xx_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_xx_xy[i] = -2.0 * g_0_xz_xx_xy[i] * b_exp + 4.0 * g_zz_xz_xx_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_xx_xz[i] = -2.0 * g_0_xz_xx_xz[i] * b_exp + 4.0 * g_zz_xz_xx_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_xx_yy[i] = -2.0 * g_0_xz_xx_yy[i] * b_exp + 4.0 * g_zz_xz_xx_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_xx_yz[i] = -2.0 * g_0_xz_xx_yz[i] * b_exp + 4.0 * g_zz_xz_xx_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_xx_zz[i] = -2.0 * g_0_xz_xx_zz[i] * b_exp + 4.0 * g_zz_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2238-2244)

    #pragma omp simd aligned(g_0_xz_xy_xx, g_0_xz_xy_xy, g_0_xz_xy_xz, g_0_xz_xy_yy, g_0_xz_xy_yz, g_0_xz_xy_zz, g_z_x_0_0_z_z_xy_xx, g_z_x_0_0_z_z_xy_xy, g_z_x_0_0_z_z_xy_xz, g_z_x_0_0_z_z_xy_yy, g_z_x_0_0_z_z_xy_yz, g_z_x_0_0_z_z_xy_zz, g_zz_xz_xy_xx, g_zz_xz_xy_xy, g_zz_xz_xy_xz, g_zz_xz_xy_yy, g_zz_xz_xy_yz, g_zz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_z_xy_xx[i] = -2.0 * g_0_xz_xy_xx[i] * b_exp + 4.0 * g_zz_xz_xy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_xy_xy[i] = -2.0 * g_0_xz_xy_xy[i] * b_exp + 4.0 * g_zz_xz_xy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_xy_xz[i] = -2.0 * g_0_xz_xy_xz[i] * b_exp + 4.0 * g_zz_xz_xy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_xy_yy[i] = -2.0 * g_0_xz_xy_yy[i] * b_exp + 4.0 * g_zz_xz_xy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_xy_yz[i] = -2.0 * g_0_xz_xy_yz[i] * b_exp + 4.0 * g_zz_xz_xy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_xy_zz[i] = -2.0 * g_0_xz_xy_zz[i] * b_exp + 4.0 * g_zz_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2244-2250)

    #pragma omp simd aligned(g_0_xz_xz_xx, g_0_xz_xz_xy, g_0_xz_xz_xz, g_0_xz_xz_yy, g_0_xz_xz_yz, g_0_xz_xz_zz, g_z_x_0_0_z_z_xz_xx, g_z_x_0_0_z_z_xz_xy, g_z_x_0_0_z_z_xz_xz, g_z_x_0_0_z_z_xz_yy, g_z_x_0_0_z_z_xz_yz, g_z_x_0_0_z_z_xz_zz, g_zz_xz_xz_xx, g_zz_xz_xz_xy, g_zz_xz_xz_xz, g_zz_xz_xz_yy, g_zz_xz_xz_yz, g_zz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_z_xz_xx[i] = -2.0 * g_0_xz_xz_xx[i] * b_exp + 4.0 * g_zz_xz_xz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_xz_xy[i] = -2.0 * g_0_xz_xz_xy[i] * b_exp + 4.0 * g_zz_xz_xz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_xz_xz[i] = -2.0 * g_0_xz_xz_xz[i] * b_exp + 4.0 * g_zz_xz_xz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_xz_yy[i] = -2.0 * g_0_xz_xz_yy[i] * b_exp + 4.0 * g_zz_xz_xz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_xz_yz[i] = -2.0 * g_0_xz_xz_yz[i] * b_exp + 4.0 * g_zz_xz_xz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_xz_zz[i] = -2.0 * g_0_xz_xz_zz[i] * b_exp + 4.0 * g_zz_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2250-2256)

    #pragma omp simd aligned(g_0_xz_yy_xx, g_0_xz_yy_xy, g_0_xz_yy_xz, g_0_xz_yy_yy, g_0_xz_yy_yz, g_0_xz_yy_zz, g_z_x_0_0_z_z_yy_xx, g_z_x_0_0_z_z_yy_xy, g_z_x_0_0_z_z_yy_xz, g_z_x_0_0_z_z_yy_yy, g_z_x_0_0_z_z_yy_yz, g_z_x_0_0_z_z_yy_zz, g_zz_xz_yy_xx, g_zz_xz_yy_xy, g_zz_xz_yy_xz, g_zz_xz_yy_yy, g_zz_xz_yy_yz, g_zz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_z_yy_xx[i] = -2.0 * g_0_xz_yy_xx[i] * b_exp + 4.0 * g_zz_xz_yy_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_yy_xy[i] = -2.0 * g_0_xz_yy_xy[i] * b_exp + 4.0 * g_zz_xz_yy_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_yy_xz[i] = -2.0 * g_0_xz_yy_xz[i] * b_exp + 4.0 * g_zz_xz_yy_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_yy_yy[i] = -2.0 * g_0_xz_yy_yy[i] * b_exp + 4.0 * g_zz_xz_yy_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_yy_yz[i] = -2.0 * g_0_xz_yy_yz[i] * b_exp + 4.0 * g_zz_xz_yy_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_yy_zz[i] = -2.0 * g_0_xz_yy_zz[i] * b_exp + 4.0 * g_zz_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2256-2262)

    #pragma omp simd aligned(g_0_xz_yz_xx, g_0_xz_yz_xy, g_0_xz_yz_xz, g_0_xz_yz_yy, g_0_xz_yz_yz, g_0_xz_yz_zz, g_z_x_0_0_z_z_yz_xx, g_z_x_0_0_z_z_yz_xy, g_z_x_0_0_z_z_yz_xz, g_z_x_0_0_z_z_yz_yy, g_z_x_0_0_z_z_yz_yz, g_z_x_0_0_z_z_yz_zz, g_zz_xz_yz_xx, g_zz_xz_yz_xy, g_zz_xz_yz_xz, g_zz_xz_yz_yy, g_zz_xz_yz_yz, g_zz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_z_yz_xx[i] = -2.0 * g_0_xz_yz_xx[i] * b_exp + 4.0 * g_zz_xz_yz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_yz_xy[i] = -2.0 * g_0_xz_yz_xy[i] * b_exp + 4.0 * g_zz_xz_yz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_yz_xz[i] = -2.0 * g_0_xz_yz_xz[i] * b_exp + 4.0 * g_zz_xz_yz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_yz_yy[i] = -2.0 * g_0_xz_yz_yy[i] * b_exp + 4.0 * g_zz_xz_yz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_yz_yz[i] = -2.0 * g_0_xz_yz_yz[i] * b_exp + 4.0 * g_zz_xz_yz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_yz_zz[i] = -2.0 * g_0_xz_yz_zz[i] * b_exp + 4.0 * g_zz_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2262-2268)

    #pragma omp simd aligned(g_0_xz_zz_xx, g_0_xz_zz_xy, g_0_xz_zz_xz, g_0_xz_zz_yy, g_0_xz_zz_yz, g_0_xz_zz_zz, g_z_x_0_0_z_z_zz_xx, g_z_x_0_0_z_z_zz_xy, g_z_x_0_0_z_z_zz_xz, g_z_x_0_0_z_z_zz_yy, g_z_x_0_0_z_z_zz_yz, g_z_x_0_0_z_z_zz_zz, g_zz_xz_zz_xx, g_zz_xz_zz_xy, g_zz_xz_zz_xz, g_zz_xz_zz_yy, g_zz_xz_zz_yz, g_zz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_z_zz_xx[i] = -2.0 * g_0_xz_zz_xx[i] * b_exp + 4.0 * g_zz_xz_zz_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_zz_xy[i] = -2.0 * g_0_xz_zz_xy[i] * b_exp + 4.0 * g_zz_xz_zz_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_zz_xz[i] = -2.0 * g_0_xz_zz_xz[i] * b_exp + 4.0 * g_zz_xz_zz_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_zz_yy[i] = -2.0 * g_0_xz_zz_yy[i] * b_exp + 4.0 * g_zz_xz_zz_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_zz_yz[i] = -2.0 * g_0_xz_zz_yz[i] * b_exp + 4.0 * g_zz_xz_zz_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_zz_zz[i] = -2.0 * g_0_xz_zz_zz[i] * b_exp + 4.0 * g_zz_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2268-2274)

    #pragma omp simd aligned(g_xz_xy_xx_xx, g_xz_xy_xx_xy, g_xz_xy_xx_xz, g_xz_xy_xx_yy, g_xz_xy_xx_yz, g_xz_xy_xx_zz, g_z_y_0_0_x_x_xx_xx, g_z_y_0_0_x_x_xx_xy, g_z_y_0_0_x_x_xx_xz, g_z_y_0_0_x_x_xx_yy, g_z_y_0_0_x_x_xx_yz, g_z_y_0_0_x_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_x_xx_xx[i] = 4.0 * g_xz_xy_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_xx_xy[i] = 4.0 * g_xz_xy_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_xx_xz[i] = 4.0 * g_xz_xy_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_xx_yy[i] = 4.0 * g_xz_xy_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_xx_yz[i] = 4.0 * g_xz_xy_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_xx_zz[i] = 4.0 * g_xz_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2274-2280)

    #pragma omp simd aligned(g_xz_xy_xy_xx, g_xz_xy_xy_xy, g_xz_xy_xy_xz, g_xz_xy_xy_yy, g_xz_xy_xy_yz, g_xz_xy_xy_zz, g_z_y_0_0_x_x_xy_xx, g_z_y_0_0_x_x_xy_xy, g_z_y_0_0_x_x_xy_xz, g_z_y_0_0_x_x_xy_yy, g_z_y_0_0_x_x_xy_yz, g_z_y_0_0_x_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_x_xy_xx[i] = 4.0 * g_xz_xy_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_xy_xy[i] = 4.0 * g_xz_xy_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_xy_xz[i] = 4.0 * g_xz_xy_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_xy_yy[i] = 4.0 * g_xz_xy_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_xy_yz[i] = 4.0 * g_xz_xy_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_xy_zz[i] = 4.0 * g_xz_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2280-2286)

    #pragma omp simd aligned(g_xz_xy_xz_xx, g_xz_xy_xz_xy, g_xz_xy_xz_xz, g_xz_xy_xz_yy, g_xz_xy_xz_yz, g_xz_xy_xz_zz, g_z_y_0_0_x_x_xz_xx, g_z_y_0_0_x_x_xz_xy, g_z_y_0_0_x_x_xz_xz, g_z_y_0_0_x_x_xz_yy, g_z_y_0_0_x_x_xz_yz, g_z_y_0_0_x_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_x_xz_xx[i] = 4.0 * g_xz_xy_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_xz_xy[i] = 4.0 * g_xz_xy_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_xz_xz[i] = 4.0 * g_xz_xy_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_xz_yy[i] = 4.0 * g_xz_xy_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_xz_yz[i] = 4.0 * g_xz_xy_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_xz_zz[i] = 4.0 * g_xz_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2286-2292)

    #pragma omp simd aligned(g_xz_xy_yy_xx, g_xz_xy_yy_xy, g_xz_xy_yy_xz, g_xz_xy_yy_yy, g_xz_xy_yy_yz, g_xz_xy_yy_zz, g_z_y_0_0_x_x_yy_xx, g_z_y_0_0_x_x_yy_xy, g_z_y_0_0_x_x_yy_xz, g_z_y_0_0_x_x_yy_yy, g_z_y_0_0_x_x_yy_yz, g_z_y_0_0_x_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_x_yy_xx[i] = 4.0 * g_xz_xy_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_yy_xy[i] = 4.0 * g_xz_xy_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_yy_xz[i] = 4.0 * g_xz_xy_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_yy_yy[i] = 4.0 * g_xz_xy_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_yy_yz[i] = 4.0 * g_xz_xy_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_yy_zz[i] = 4.0 * g_xz_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2292-2298)

    #pragma omp simd aligned(g_xz_xy_yz_xx, g_xz_xy_yz_xy, g_xz_xy_yz_xz, g_xz_xy_yz_yy, g_xz_xy_yz_yz, g_xz_xy_yz_zz, g_z_y_0_0_x_x_yz_xx, g_z_y_0_0_x_x_yz_xy, g_z_y_0_0_x_x_yz_xz, g_z_y_0_0_x_x_yz_yy, g_z_y_0_0_x_x_yz_yz, g_z_y_0_0_x_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_x_yz_xx[i] = 4.0 * g_xz_xy_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_yz_xy[i] = 4.0 * g_xz_xy_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_yz_xz[i] = 4.0 * g_xz_xy_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_yz_yy[i] = 4.0 * g_xz_xy_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_yz_yz[i] = 4.0 * g_xz_xy_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_yz_zz[i] = 4.0 * g_xz_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2298-2304)

    #pragma omp simd aligned(g_xz_xy_zz_xx, g_xz_xy_zz_xy, g_xz_xy_zz_xz, g_xz_xy_zz_yy, g_xz_xy_zz_yz, g_xz_xy_zz_zz, g_z_y_0_0_x_x_zz_xx, g_z_y_0_0_x_x_zz_xy, g_z_y_0_0_x_x_zz_xz, g_z_y_0_0_x_x_zz_yy, g_z_y_0_0_x_x_zz_yz, g_z_y_0_0_x_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_x_zz_xx[i] = 4.0 * g_xz_xy_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_zz_xy[i] = 4.0 * g_xz_xy_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_zz_xz[i] = 4.0 * g_xz_xy_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_zz_yy[i] = 4.0 * g_xz_xy_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_zz_yz[i] = 4.0 * g_xz_xy_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_zz_zz[i] = 4.0 * g_xz_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2304-2310)

    #pragma omp simd aligned(g_xz_0_xx_xx, g_xz_0_xx_xy, g_xz_0_xx_xz, g_xz_0_xx_yy, g_xz_0_xx_yz, g_xz_0_xx_zz, g_xz_yy_xx_xx, g_xz_yy_xx_xy, g_xz_yy_xx_xz, g_xz_yy_xx_yy, g_xz_yy_xx_yz, g_xz_yy_xx_zz, g_z_y_0_0_x_y_xx_xx, g_z_y_0_0_x_y_xx_xy, g_z_y_0_0_x_y_xx_xz, g_z_y_0_0_x_y_xx_yy, g_z_y_0_0_x_y_xx_yz, g_z_y_0_0_x_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_y_xx_xx[i] = -2.0 * g_xz_0_xx_xx[i] * a_exp + 4.0 * g_xz_yy_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_xx_xy[i] = -2.0 * g_xz_0_xx_xy[i] * a_exp + 4.0 * g_xz_yy_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_xx_xz[i] = -2.0 * g_xz_0_xx_xz[i] * a_exp + 4.0 * g_xz_yy_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_xx_yy[i] = -2.0 * g_xz_0_xx_yy[i] * a_exp + 4.0 * g_xz_yy_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_xx_yz[i] = -2.0 * g_xz_0_xx_yz[i] * a_exp + 4.0 * g_xz_yy_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_xx_zz[i] = -2.0 * g_xz_0_xx_zz[i] * a_exp + 4.0 * g_xz_yy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2310-2316)

    #pragma omp simd aligned(g_xz_0_xy_xx, g_xz_0_xy_xy, g_xz_0_xy_xz, g_xz_0_xy_yy, g_xz_0_xy_yz, g_xz_0_xy_zz, g_xz_yy_xy_xx, g_xz_yy_xy_xy, g_xz_yy_xy_xz, g_xz_yy_xy_yy, g_xz_yy_xy_yz, g_xz_yy_xy_zz, g_z_y_0_0_x_y_xy_xx, g_z_y_0_0_x_y_xy_xy, g_z_y_0_0_x_y_xy_xz, g_z_y_0_0_x_y_xy_yy, g_z_y_0_0_x_y_xy_yz, g_z_y_0_0_x_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_y_xy_xx[i] = -2.0 * g_xz_0_xy_xx[i] * a_exp + 4.0 * g_xz_yy_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_xy_xy[i] = -2.0 * g_xz_0_xy_xy[i] * a_exp + 4.0 * g_xz_yy_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_xy_xz[i] = -2.0 * g_xz_0_xy_xz[i] * a_exp + 4.0 * g_xz_yy_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_xy_yy[i] = -2.0 * g_xz_0_xy_yy[i] * a_exp + 4.0 * g_xz_yy_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_xy_yz[i] = -2.0 * g_xz_0_xy_yz[i] * a_exp + 4.0 * g_xz_yy_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_xy_zz[i] = -2.0 * g_xz_0_xy_zz[i] * a_exp + 4.0 * g_xz_yy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2316-2322)

    #pragma omp simd aligned(g_xz_0_xz_xx, g_xz_0_xz_xy, g_xz_0_xz_xz, g_xz_0_xz_yy, g_xz_0_xz_yz, g_xz_0_xz_zz, g_xz_yy_xz_xx, g_xz_yy_xz_xy, g_xz_yy_xz_xz, g_xz_yy_xz_yy, g_xz_yy_xz_yz, g_xz_yy_xz_zz, g_z_y_0_0_x_y_xz_xx, g_z_y_0_0_x_y_xz_xy, g_z_y_0_0_x_y_xz_xz, g_z_y_0_0_x_y_xz_yy, g_z_y_0_0_x_y_xz_yz, g_z_y_0_0_x_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_y_xz_xx[i] = -2.0 * g_xz_0_xz_xx[i] * a_exp + 4.0 * g_xz_yy_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_xz_xy[i] = -2.0 * g_xz_0_xz_xy[i] * a_exp + 4.0 * g_xz_yy_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_xz_xz[i] = -2.0 * g_xz_0_xz_xz[i] * a_exp + 4.0 * g_xz_yy_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_xz_yy[i] = -2.0 * g_xz_0_xz_yy[i] * a_exp + 4.0 * g_xz_yy_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_xz_yz[i] = -2.0 * g_xz_0_xz_yz[i] * a_exp + 4.0 * g_xz_yy_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_xz_zz[i] = -2.0 * g_xz_0_xz_zz[i] * a_exp + 4.0 * g_xz_yy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2322-2328)

    #pragma omp simd aligned(g_xz_0_yy_xx, g_xz_0_yy_xy, g_xz_0_yy_xz, g_xz_0_yy_yy, g_xz_0_yy_yz, g_xz_0_yy_zz, g_xz_yy_yy_xx, g_xz_yy_yy_xy, g_xz_yy_yy_xz, g_xz_yy_yy_yy, g_xz_yy_yy_yz, g_xz_yy_yy_zz, g_z_y_0_0_x_y_yy_xx, g_z_y_0_0_x_y_yy_xy, g_z_y_0_0_x_y_yy_xz, g_z_y_0_0_x_y_yy_yy, g_z_y_0_0_x_y_yy_yz, g_z_y_0_0_x_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_y_yy_xx[i] = -2.0 * g_xz_0_yy_xx[i] * a_exp + 4.0 * g_xz_yy_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_yy_xy[i] = -2.0 * g_xz_0_yy_xy[i] * a_exp + 4.0 * g_xz_yy_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_yy_xz[i] = -2.0 * g_xz_0_yy_xz[i] * a_exp + 4.0 * g_xz_yy_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_yy_yy[i] = -2.0 * g_xz_0_yy_yy[i] * a_exp + 4.0 * g_xz_yy_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_yy_yz[i] = -2.0 * g_xz_0_yy_yz[i] * a_exp + 4.0 * g_xz_yy_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_yy_zz[i] = -2.0 * g_xz_0_yy_zz[i] * a_exp + 4.0 * g_xz_yy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2328-2334)

    #pragma omp simd aligned(g_xz_0_yz_xx, g_xz_0_yz_xy, g_xz_0_yz_xz, g_xz_0_yz_yy, g_xz_0_yz_yz, g_xz_0_yz_zz, g_xz_yy_yz_xx, g_xz_yy_yz_xy, g_xz_yy_yz_xz, g_xz_yy_yz_yy, g_xz_yy_yz_yz, g_xz_yy_yz_zz, g_z_y_0_0_x_y_yz_xx, g_z_y_0_0_x_y_yz_xy, g_z_y_0_0_x_y_yz_xz, g_z_y_0_0_x_y_yz_yy, g_z_y_0_0_x_y_yz_yz, g_z_y_0_0_x_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_y_yz_xx[i] = -2.0 * g_xz_0_yz_xx[i] * a_exp + 4.0 * g_xz_yy_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_yz_xy[i] = -2.0 * g_xz_0_yz_xy[i] * a_exp + 4.0 * g_xz_yy_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_yz_xz[i] = -2.0 * g_xz_0_yz_xz[i] * a_exp + 4.0 * g_xz_yy_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_yz_yy[i] = -2.0 * g_xz_0_yz_yy[i] * a_exp + 4.0 * g_xz_yy_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_yz_yz[i] = -2.0 * g_xz_0_yz_yz[i] * a_exp + 4.0 * g_xz_yy_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_yz_zz[i] = -2.0 * g_xz_0_yz_zz[i] * a_exp + 4.0 * g_xz_yy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2334-2340)

    #pragma omp simd aligned(g_xz_0_zz_xx, g_xz_0_zz_xy, g_xz_0_zz_xz, g_xz_0_zz_yy, g_xz_0_zz_yz, g_xz_0_zz_zz, g_xz_yy_zz_xx, g_xz_yy_zz_xy, g_xz_yy_zz_xz, g_xz_yy_zz_yy, g_xz_yy_zz_yz, g_xz_yy_zz_zz, g_z_y_0_0_x_y_zz_xx, g_z_y_0_0_x_y_zz_xy, g_z_y_0_0_x_y_zz_xz, g_z_y_0_0_x_y_zz_yy, g_z_y_0_0_x_y_zz_yz, g_z_y_0_0_x_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_y_zz_xx[i] = -2.0 * g_xz_0_zz_xx[i] * a_exp + 4.0 * g_xz_yy_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_zz_xy[i] = -2.0 * g_xz_0_zz_xy[i] * a_exp + 4.0 * g_xz_yy_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_zz_xz[i] = -2.0 * g_xz_0_zz_xz[i] * a_exp + 4.0 * g_xz_yy_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_zz_yy[i] = -2.0 * g_xz_0_zz_yy[i] * a_exp + 4.0 * g_xz_yy_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_zz_yz[i] = -2.0 * g_xz_0_zz_yz[i] * a_exp + 4.0 * g_xz_yy_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_zz_zz[i] = -2.0 * g_xz_0_zz_zz[i] * a_exp + 4.0 * g_xz_yy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2340-2346)

    #pragma omp simd aligned(g_xz_yz_xx_xx, g_xz_yz_xx_xy, g_xz_yz_xx_xz, g_xz_yz_xx_yy, g_xz_yz_xx_yz, g_xz_yz_xx_zz, g_z_y_0_0_x_z_xx_xx, g_z_y_0_0_x_z_xx_xy, g_z_y_0_0_x_z_xx_xz, g_z_y_0_0_x_z_xx_yy, g_z_y_0_0_x_z_xx_yz, g_z_y_0_0_x_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_z_xx_xx[i] = 4.0 * g_xz_yz_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_xx_xy[i] = 4.0 * g_xz_yz_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_xx_xz[i] = 4.0 * g_xz_yz_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_xx_yy[i] = 4.0 * g_xz_yz_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_xx_yz[i] = 4.0 * g_xz_yz_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_xx_zz[i] = 4.0 * g_xz_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2346-2352)

    #pragma omp simd aligned(g_xz_yz_xy_xx, g_xz_yz_xy_xy, g_xz_yz_xy_xz, g_xz_yz_xy_yy, g_xz_yz_xy_yz, g_xz_yz_xy_zz, g_z_y_0_0_x_z_xy_xx, g_z_y_0_0_x_z_xy_xy, g_z_y_0_0_x_z_xy_xz, g_z_y_0_0_x_z_xy_yy, g_z_y_0_0_x_z_xy_yz, g_z_y_0_0_x_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_z_xy_xx[i] = 4.0 * g_xz_yz_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_xy_xy[i] = 4.0 * g_xz_yz_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_xy_xz[i] = 4.0 * g_xz_yz_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_xy_yy[i] = 4.0 * g_xz_yz_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_xy_yz[i] = 4.0 * g_xz_yz_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_xy_zz[i] = 4.0 * g_xz_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2352-2358)

    #pragma omp simd aligned(g_xz_yz_xz_xx, g_xz_yz_xz_xy, g_xz_yz_xz_xz, g_xz_yz_xz_yy, g_xz_yz_xz_yz, g_xz_yz_xz_zz, g_z_y_0_0_x_z_xz_xx, g_z_y_0_0_x_z_xz_xy, g_z_y_0_0_x_z_xz_xz, g_z_y_0_0_x_z_xz_yy, g_z_y_0_0_x_z_xz_yz, g_z_y_0_0_x_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_z_xz_xx[i] = 4.0 * g_xz_yz_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_xz_xy[i] = 4.0 * g_xz_yz_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_xz_xz[i] = 4.0 * g_xz_yz_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_xz_yy[i] = 4.0 * g_xz_yz_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_xz_yz[i] = 4.0 * g_xz_yz_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_xz_zz[i] = 4.0 * g_xz_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2358-2364)

    #pragma omp simd aligned(g_xz_yz_yy_xx, g_xz_yz_yy_xy, g_xz_yz_yy_xz, g_xz_yz_yy_yy, g_xz_yz_yy_yz, g_xz_yz_yy_zz, g_z_y_0_0_x_z_yy_xx, g_z_y_0_0_x_z_yy_xy, g_z_y_0_0_x_z_yy_xz, g_z_y_0_0_x_z_yy_yy, g_z_y_0_0_x_z_yy_yz, g_z_y_0_0_x_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_z_yy_xx[i] = 4.0 * g_xz_yz_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_yy_xy[i] = 4.0 * g_xz_yz_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_yy_xz[i] = 4.0 * g_xz_yz_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_yy_yy[i] = 4.0 * g_xz_yz_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_yy_yz[i] = 4.0 * g_xz_yz_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_yy_zz[i] = 4.0 * g_xz_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2364-2370)

    #pragma omp simd aligned(g_xz_yz_yz_xx, g_xz_yz_yz_xy, g_xz_yz_yz_xz, g_xz_yz_yz_yy, g_xz_yz_yz_yz, g_xz_yz_yz_zz, g_z_y_0_0_x_z_yz_xx, g_z_y_0_0_x_z_yz_xy, g_z_y_0_0_x_z_yz_xz, g_z_y_0_0_x_z_yz_yy, g_z_y_0_0_x_z_yz_yz, g_z_y_0_0_x_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_z_yz_xx[i] = 4.0 * g_xz_yz_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_yz_xy[i] = 4.0 * g_xz_yz_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_yz_xz[i] = 4.0 * g_xz_yz_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_yz_yy[i] = 4.0 * g_xz_yz_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_yz_yz[i] = 4.0 * g_xz_yz_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_yz_zz[i] = 4.0 * g_xz_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2370-2376)

    #pragma omp simd aligned(g_xz_yz_zz_xx, g_xz_yz_zz_xy, g_xz_yz_zz_xz, g_xz_yz_zz_yy, g_xz_yz_zz_yz, g_xz_yz_zz_zz, g_z_y_0_0_x_z_zz_xx, g_z_y_0_0_x_z_zz_xy, g_z_y_0_0_x_z_zz_xz, g_z_y_0_0_x_z_zz_yy, g_z_y_0_0_x_z_zz_yz, g_z_y_0_0_x_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_z_zz_xx[i] = 4.0 * g_xz_yz_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_zz_xy[i] = 4.0 * g_xz_yz_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_zz_xz[i] = 4.0 * g_xz_yz_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_zz_yy[i] = 4.0 * g_xz_yz_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_zz_yz[i] = 4.0 * g_xz_yz_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_zz_zz[i] = 4.0 * g_xz_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2376-2382)

    #pragma omp simd aligned(g_yz_xy_xx_xx, g_yz_xy_xx_xy, g_yz_xy_xx_xz, g_yz_xy_xx_yy, g_yz_xy_xx_yz, g_yz_xy_xx_zz, g_z_y_0_0_y_x_xx_xx, g_z_y_0_0_y_x_xx_xy, g_z_y_0_0_y_x_xx_xz, g_z_y_0_0_y_x_xx_yy, g_z_y_0_0_y_x_xx_yz, g_z_y_0_0_y_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_x_xx_xx[i] = 4.0 * g_yz_xy_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_xx_xy[i] = 4.0 * g_yz_xy_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_xx_xz[i] = 4.0 * g_yz_xy_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_xx_yy[i] = 4.0 * g_yz_xy_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_xx_yz[i] = 4.0 * g_yz_xy_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_xx_zz[i] = 4.0 * g_yz_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2382-2388)

    #pragma omp simd aligned(g_yz_xy_xy_xx, g_yz_xy_xy_xy, g_yz_xy_xy_xz, g_yz_xy_xy_yy, g_yz_xy_xy_yz, g_yz_xy_xy_zz, g_z_y_0_0_y_x_xy_xx, g_z_y_0_0_y_x_xy_xy, g_z_y_0_0_y_x_xy_xz, g_z_y_0_0_y_x_xy_yy, g_z_y_0_0_y_x_xy_yz, g_z_y_0_0_y_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_x_xy_xx[i] = 4.0 * g_yz_xy_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_xy_xy[i] = 4.0 * g_yz_xy_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_xy_xz[i] = 4.0 * g_yz_xy_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_xy_yy[i] = 4.0 * g_yz_xy_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_xy_yz[i] = 4.0 * g_yz_xy_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_xy_zz[i] = 4.0 * g_yz_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2388-2394)

    #pragma omp simd aligned(g_yz_xy_xz_xx, g_yz_xy_xz_xy, g_yz_xy_xz_xz, g_yz_xy_xz_yy, g_yz_xy_xz_yz, g_yz_xy_xz_zz, g_z_y_0_0_y_x_xz_xx, g_z_y_0_0_y_x_xz_xy, g_z_y_0_0_y_x_xz_xz, g_z_y_0_0_y_x_xz_yy, g_z_y_0_0_y_x_xz_yz, g_z_y_0_0_y_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_x_xz_xx[i] = 4.0 * g_yz_xy_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_xz_xy[i] = 4.0 * g_yz_xy_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_xz_xz[i] = 4.0 * g_yz_xy_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_xz_yy[i] = 4.0 * g_yz_xy_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_xz_yz[i] = 4.0 * g_yz_xy_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_xz_zz[i] = 4.0 * g_yz_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2394-2400)

    #pragma omp simd aligned(g_yz_xy_yy_xx, g_yz_xy_yy_xy, g_yz_xy_yy_xz, g_yz_xy_yy_yy, g_yz_xy_yy_yz, g_yz_xy_yy_zz, g_z_y_0_0_y_x_yy_xx, g_z_y_0_0_y_x_yy_xy, g_z_y_0_0_y_x_yy_xz, g_z_y_0_0_y_x_yy_yy, g_z_y_0_0_y_x_yy_yz, g_z_y_0_0_y_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_x_yy_xx[i] = 4.0 * g_yz_xy_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_yy_xy[i] = 4.0 * g_yz_xy_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_yy_xz[i] = 4.0 * g_yz_xy_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_yy_yy[i] = 4.0 * g_yz_xy_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_yy_yz[i] = 4.0 * g_yz_xy_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_yy_zz[i] = 4.0 * g_yz_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2400-2406)

    #pragma omp simd aligned(g_yz_xy_yz_xx, g_yz_xy_yz_xy, g_yz_xy_yz_xz, g_yz_xy_yz_yy, g_yz_xy_yz_yz, g_yz_xy_yz_zz, g_z_y_0_0_y_x_yz_xx, g_z_y_0_0_y_x_yz_xy, g_z_y_0_0_y_x_yz_xz, g_z_y_0_0_y_x_yz_yy, g_z_y_0_0_y_x_yz_yz, g_z_y_0_0_y_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_x_yz_xx[i] = 4.0 * g_yz_xy_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_yz_xy[i] = 4.0 * g_yz_xy_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_yz_xz[i] = 4.0 * g_yz_xy_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_yz_yy[i] = 4.0 * g_yz_xy_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_yz_yz[i] = 4.0 * g_yz_xy_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_yz_zz[i] = 4.0 * g_yz_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2406-2412)

    #pragma omp simd aligned(g_yz_xy_zz_xx, g_yz_xy_zz_xy, g_yz_xy_zz_xz, g_yz_xy_zz_yy, g_yz_xy_zz_yz, g_yz_xy_zz_zz, g_z_y_0_0_y_x_zz_xx, g_z_y_0_0_y_x_zz_xy, g_z_y_0_0_y_x_zz_xz, g_z_y_0_0_y_x_zz_yy, g_z_y_0_0_y_x_zz_yz, g_z_y_0_0_y_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_x_zz_xx[i] = 4.0 * g_yz_xy_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_zz_xy[i] = 4.0 * g_yz_xy_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_zz_xz[i] = 4.0 * g_yz_xy_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_zz_yy[i] = 4.0 * g_yz_xy_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_zz_yz[i] = 4.0 * g_yz_xy_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_zz_zz[i] = 4.0 * g_yz_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2412-2418)

    #pragma omp simd aligned(g_yz_0_xx_xx, g_yz_0_xx_xy, g_yz_0_xx_xz, g_yz_0_xx_yy, g_yz_0_xx_yz, g_yz_0_xx_zz, g_yz_yy_xx_xx, g_yz_yy_xx_xy, g_yz_yy_xx_xz, g_yz_yy_xx_yy, g_yz_yy_xx_yz, g_yz_yy_xx_zz, g_z_y_0_0_y_y_xx_xx, g_z_y_0_0_y_y_xx_xy, g_z_y_0_0_y_y_xx_xz, g_z_y_0_0_y_y_xx_yy, g_z_y_0_0_y_y_xx_yz, g_z_y_0_0_y_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_y_xx_xx[i] = -2.0 * g_yz_0_xx_xx[i] * a_exp + 4.0 * g_yz_yy_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_xx_xy[i] = -2.0 * g_yz_0_xx_xy[i] * a_exp + 4.0 * g_yz_yy_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_xx_xz[i] = -2.0 * g_yz_0_xx_xz[i] * a_exp + 4.0 * g_yz_yy_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_xx_yy[i] = -2.0 * g_yz_0_xx_yy[i] * a_exp + 4.0 * g_yz_yy_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_xx_yz[i] = -2.0 * g_yz_0_xx_yz[i] * a_exp + 4.0 * g_yz_yy_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_xx_zz[i] = -2.0 * g_yz_0_xx_zz[i] * a_exp + 4.0 * g_yz_yy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2418-2424)

    #pragma omp simd aligned(g_yz_0_xy_xx, g_yz_0_xy_xy, g_yz_0_xy_xz, g_yz_0_xy_yy, g_yz_0_xy_yz, g_yz_0_xy_zz, g_yz_yy_xy_xx, g_yz_yy_xy_xy, g_yz_yy_xy_xz, g_yz_yy_xy_yy, g_yz_yy_xy_yz, g_yz_yy_xy_zz, g_z_y_0_0_y_y_xy_xx, g_z_y_0_0_y_y_xy_xy, g_z_y_0_0_y_y_xy_xz, g_z_y_0_0_y_y_xy_yy, g_z_y_0_0_y_y_xy_yz, g_z_y_0_0_y_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_y_xy_xx[i] = -2.0 * g_yz_0_xy_xx[i] * a_exp + 4.0 * g_yz_yy_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_xy_xy[i] = -2.0 * g_yz_0_xy_xy[i] * a_exp + 4.0 * g_yz_yy_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_xy_xz[i] = -2.0 * g_yz_0_xy_xz[i] * a_exp + 4.0 * g_yz_yy_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_xy_yy[i] = -2.0 * g_yz_0_xy_yy[i] * a_exp + 4.0 * g_yz_yy_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_xy_yz[i] = -2.0 * g_yz_0_xy_yz[i] * a_exp + 4.0 * g_yz_yy_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_xy_zz[i] = -2.0 * g_yz_0_xy_zz[i] * a_exp + 4.0 * g_yz_yy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2424-2430)

    #pragma omp simd aligned(g_yz_0_xz_xx, g_yz_0_xz_xy, g_yz_0_xz_xz, g_yz_0_xz_yy, g_yz_0_xz_yz, g_yz_0_xz_zz, g_yz_yy_xz_xx, g_yz_yy_xz_xy, g_yz_yy_xz_xz, g_yz_yy_xz_yy, g_yz_yy_xz_yz, g_yz_yy_xz_zz, g_z_y_0_0_y_y_xz_xx, g_z_y_0_0_y_y_xz_xy, g_z_y_0_0_y_y_xz_xz, g_z_y_0_0_y_y_xz_yy, g_z_y_0_0_y_y_xz_yz, g_z_y_0_0_y_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_y_xz_xx[i] = -2.0 * g_yz_0_xz_xx[i] * a_exp + 4.0 * g_yz_yy_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_xz_xy[i] = -2.0 * g_yz_0_xz_xy[i] * a_exp + 4.0 * g_yz_yy_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_xz_xz[i] = -2.0 * g_yz_0_xz_xz[i] * a_exp + 4.0 * g_yz_yy_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_xz_yy[i] = -2.0 * g_yz_0_xz_yy[i] * a_exp + 4.0 * g_yz_yy_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_xz_yz[i] = -2.0 * g_yz_0_xz_yz[i] * a_exp + 4.0 * g_yz_yy_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_xz_zz[i] = -2.0 * g_yz_0_xz_zz[i] * a_exp + 4.0 * g_yz_yy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2430-2436)

    #pragma omp simd aligned(g_yz_0_yy_xx, g_yz_0_yy_xy, g_yz_0_yy_xz, g_yz_0_yy_yy, g_yz_0_yy_yz, g_yz_0_yy_zz, g_yz_yy_yy_xx, g_yz_yy_yy_xy, g_yz_yy_yy_xz, g_yz_yy_yy_yy, g_yz_yy_yy_yz, g_yz_yy_yy_zz, g_z_y_0_0_y_y_yy_xx, g_z_y_0_0_y_y_yy_xy, g_z_y_0_0_y_y_yy_xz, g_z_y_0_0_y_y_yy_yy, g_z_y_0_0_y_y_yy_yz, g_z_y_0_0_y_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_y_yy_xx[i] = -2.0 * g_yz_0_yy_xx[i] * a_exp + 4.0 * g_yz_yy_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_yy_xy[i] = -2.0 * g_yz_0_yy_xy[i] * a_exp + 4.0 * g_yz_yy_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_yy_xz[i] = -2.0 * g_yz_0_yy_xz[i] * a_exp + 4.0 * g_yz_yy_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_yy_yy[i] = -2.0 * g_yz_0_yy_yy[i] * a_exp + 4.0 * g_yz_yy_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_yy_yz[i] = -2.0 * g_yz_0_yy_yz[i] * a_exp + 4.0 * g_yz_yy_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_yy_zz[i] = -2.0 * g_yz_0_yy_zz[i] * a_exp + 4.0 * g_yz_yy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2436-2442)

    #pragma omp simd aligned(g_yz_0_yz_xx, g_yz_0_yz_xy, g_yz_0_yz_xz, g_yz_0_yz_yy, g_yz_0_yz_yz, g_yz_0_yz_zz, g_yz_yy_yz_xx, g_yz_yy_yz_xy, g_yz_yy_yz_xz, g_yz_yy_yz_yy, g_yz_yy_yz_yz, g_yz_yy_yz_zz, g_z_y_0_0_y_y_yz_xx, g_z_y_0_0_y_y_yz_xy, g_z_y_0_0_y_y_yz_xz, g_z_y_0_0_y_y_yz_yy, g_z_y_0_0_y_y_yz_yz, g_z_y_0_0_y_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_y_yz_xx[i] = -2.0 * g_yz_0_yz_xx[i] * a_exp + 4.0 * g_yz_yy_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_yz_xy[i] = -2.0 * g_yz_0_yz_xy[i] * a_exp + 4.0 * g_yz_yy_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_yz_xz[i] = -2.0 * g_yz_0_yz_xz[i] * a_exp + 4.0 * g_yz_yy_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_yz_yy[i] = -2.0 * g_yz_0_yz_yy[i] * a_exp + 4.0 * g_yz_yy_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_yz_yz[i] = -2.0 * g_yz_0_yz_yz[i] * a_exp + 4.0 * g_yz_yy_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_yz_zz[i] = -2.0 * g_yz_0_yz_zz[i] * a_exp + 4.0 * g_yz_yy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2442-2448)

    #pragma omp simd aligned(g_yz_0_zz_xx, g_yz_0_zz_xy, g_yz_0_zz_xz, g_yz_0_zz_yy, g_yz_0_zz_yz, g_yz_0_zz_zz, g_yz_yy_zz_xx, g_yz_yy_zz_xy, g_yz_yy_zz_xz, g_yz_yy_zz_yy, g_yz_yy_zz_yz, g_yz_yy_zz_zz, g_z_y_0_0_y_y_zz_xx, g_z_y_0_0_y_y_zz_xy, g_z_y_0_0_y_y_zz_xz, g_z_y_0_0_y_y_zz_yy, g_z_y_0_0_y_y_zz_yz, g_z_y_0_0_y_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_y_zz_xx[i] = -2.0 * g_yz_0_zz_xx[i] * a_exp + 4.0 * g_yz_yy_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_zz_xy[i] = -2.0 * g_yz_0_zz_xy[i] * a_exp + 4.0 * g_yz_yy_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_zz_xz[i] = -2.0 * g_yz_0_zz_xz[i] * a_exp + 4.0 * g_yz_yy_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_zz_yy[i] = -2.0 * g_yz_0_zz_yy[i] * a_exp + 4.0 * g_yz_yy_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_zz_yz[i] = -2.0 * g_yz_0_zz_yz[i] * a_exp + 4.0 * g_yz_yy_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_zz_zz[i] = -2.0 * g_yz_0_zz_zz[i] * a_exp + 4.0 * g_yz_yy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2448-2454)

    #pragma omp simd aligned(g_yz_yz_xx_xx, g_yz_yz_xx_xy, g_yz_yz_xx_xz, g_yz_yz_xx_yy, g_yz_yz_xx_yz, g_yz_yz_xx_zz, g_z_y_0_0_y_z_xx_xx, g_z_y_0_0_y_z_xx_xy, g_z_y_0_0_y_z_xx_xz, g_z_y_0_0_y_z_xx_yy, g_z_y_0_0_y_z_xx_yz, g_z_y_0_0_y_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_z_xx_xx[i] = 4.0 * g_yz_yz_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_xx_xy[i] = 4.0 * g_yz_yz_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_xx_xz[i] = 4.0 * g_yz_yz_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_xx_yy[i] = 4.0 * g_yz_yz_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_xx_yz[i] = 4.0 * g_yz_yz_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_xx_zz[i] = 4.0 * g_yz_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2454-2460)

    #pragma omp simd aligned(g_yz_yz_xy_xx, g_yz_yz_xy_xy, g_yz_yz_xy_xz, g_yz_yz_xy_yy, g_yz_yz_xy_yz, g_yz_yz_xy_zz, g_z_y_0_0_y_z_xy_xx, g_z_y_0_0_y_z_xy_xy, g_z_y_0_0_y_z_xy_xz, g_z_y_0_0_y_z_xy_yy, g_z_y_0_0_y_z_xy_yz, g_z_y_0_0_y_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_z_xy_xx[i] = 4.0 * g_yz_yz_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_xy_xy[i] = 4.0 * g_yz_yz_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_xy_xz[i] = 4.0 * g_yz_yz_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_xy_yy[i] = 4.0 * g_yz_yz_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_xy_yz[i] = 4.0 * g_yz_yz_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_xy_zz[i] = 4.0 * g_yz_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2460-2466)

    #pragma omp simd aligned(g_yz_yz_xz_xx, g_yz_yz_xz_xy, g_yz_yz_xz_xz, g_yz_yz_xz_yy, g_yz_yz_xz_yz, g_yz_yz_xz_zz, g_z_y_0_0_y_z_xz_xx, g_z_y_0_0_y_z_xz_xy, g_z_y_0_0_y_z_xz_xz, g_z_y_0_0_y_z_xz_yy, g_z_y_0_0_y_z_xz_yz, g_z_y_0_0_y_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_z_xz_xx[i] = 4.0 * g_yz_yz_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_xz_xy[i] = 4.0 * g_yz_yz_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_xz_xz[i] = 4.0 * g_yz_yz_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_xz_yy[i] = 4.0 * g_yz_yz_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_xz_yz[i] = 4.0 * g_yz_yz_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_xz_zz[i] = 4.0 * g_yz_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2466-2472)

    #pragma omp simd aligned(g_yz_yz_yy_xx, g_yz_yz_yy_xy, g_yz_yz_yy_xz, g_yz_yz_yy_yy, g_yz_yz_yy_yz, g_yz_yz_yy_zz, g_z_y_0_0_y_z_yy_xx, g_z_y_0_0_y_z_yy_xy, g_z_y_0_0_y_z_yy_xz, g_z_y_0_0_y_z_yy_yy, g_z_y_0_0_y_z_yy_yz, g_z_y_0_0_y_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_z_yy_xx[i] = 4.0 * g_yz_yz_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_yy_xy[i] = 4.0 * g_yz_yz_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_yy_xz[i] = 4.0 * g_yz_yz_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_yy_yy[i] = 4.0 * g_yz_yz_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_yy_yz[i] = 4.0 * g_yz_yz_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_yy_zz[i] = 4.0 * g_yz_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2472-2478)

    #pragma omp simd aligned(g_yz_yz_yz_xx, g_yz_yz_yz_xy, g_yz_yz_yz_xz, g_yz_yz_yz_yy, g_yz_yz_yz_yz, g_yz_yz_yz_zz, g_z_y_0_0_y_z_yz_xx, g_z_y_0_0_y_z_yz_xy, g_z_y_0_0_y_z_yz_xz, g_z_y_0_0_y_z_yz_yy, g_z_y_0_0_y_z_yz_yz, g_z_y_0_0_y_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_z_yz_xx[i] = 4.0 * g_yz_yz_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_yz_xy[i] = 4.0 * g_yz_yz_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_yz_xz[i] = 4.0 * g_yz_yz_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_yz_yy[i] = 4.0 * g_yz_yz_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_yz_yz[i] = 4.0 * g_yz_yz_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_yz_zz[i] = 4.0 * g_yz_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2478-2484)

    #pragma omp simd aligned(g_yz_yz_zz_xx, g_yz_yz_zz_xy, g_yz_yz_zz_xz, g_yz_yz_zz_yy, g_yz_yz_zz_yz, g_yz_yz_zz_zz, g_z_y_0_0_y_z_zz_xx, g_z_y_0_0_y_z_zz_xy, g_z_y_0_0_y_z_zz_xz, g_z_y_0_0_y_z_zz_yy, g_z_y_0_0_y_z_zz_yz, g_z_y_0_0_y_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_z_zz_xx[i] = 4.0 * g_yz_yz_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_zz_xy[i] = 4.0 * g_yz_yz_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_zz_xz[i] = 4.0 * g_yz_yz_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_zz_yy[i] = 4.0 * g_yz_yz_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_zz_yz[i] = 4.0 * g_yz_yz_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_zz_zz[i] = 4.0 * g_yz_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2484-2490)

    #pragma omp simd aligned(g_0_xy_xx_xx, g_0_xy_xx_xy, g_0_xy_xx_xz, g_0_xy_xx_yy, g_0_xy_xx_yz, g_0_xy_xx_zz, g_z_y_0_0_z_x_xx_xx, g_z_y_0_0_z_x_xx_xy, g_z_y_0_0_z_x_xx_xz, g_z_y_0_0_z_x_xx_yy, g_z_y_0_0_z_x_xx_yz, g_z_y_0_0_z_x_xx_zz, g_zz_xy_xx_xx, g_zz_xy_xx_xy, g_zz_xy_xx_xz, g_zz_xy_xx_yy, g_zz_xy_xx_yz, g_zz_xy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_x_xx_xx[i] = -2.0 * g_0_xy_xx_xx[i] * b_exp + 4.0 * g_zz_xy_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_xx_xy[i] = -2.0 * g_0_xy_xx_xy[i] * b_exp + 4.0 * g_zz_xy_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_xx_xz[i] = -2.0 * g_0_xy_xx_xz[i] * b_exp + 4.0 * g_zz_xy_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_xx_yy[i] = -2.0 * g_0_xy_xx_yy[i] * b_exp + 4.0 * g_zz_xy_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_xx_yz[i] = -2.0 * g_0_xy_xx_yz[i] * b_exp + 4.0 * g_zz_xy_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_xx_zz[i] = -2.0 * g_0_xy_xx_zz[i] * b_exp + 4.0 * g_zz_xy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2490-2496)

    #pragma omp simd aligned(g_0_xy_xy_xx, g_0_xy_xy_xy, g_0_xy_xy_xz, g_0_xy_xy_yy, g_0_xy_xy_yz, g_0_xy_xy_zz, g_z_y_0_0_z_x_xy_xx, g_z_y_0_0_z_x_xy_xy, g_z_y_0_0_z_x_xy_xz, g_z_y_0_0_z_x_xy_yy, g_z_y_0_0_z_x_xy_yz, g_z_y_0_0_z_x_xy_zz, g_zz_xy_xy_xx, g_zz_xy_xy_xy, g_zz_xy_xy_xz, g_zz_xy_xy_yy, g_zz_xy_xy_yz, g_zz_xy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_x_xy_xx[i] = -2.0 * g_0_xy_xy_xx[i] * b_exp + 4.0 * g_zz_xy_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_xy_xy[i] = -2.0 * g_0_xy_xy_xy[i] * b_exp + 4.0 * g_zz_xy_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_xy_xz[i] = -2.0 * g_0_xy_xy_xz[i] * b_exp + 4.0 * g_zz_xy_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_xy_yy[i] = -2.0 * g_0_xy_xy_yy[i] * b_exp + 4.0 * g_zz_xy_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_xy_yz[i] = -2.0 * g_0_xy_xy_yz[i] * b_exp + 4.0 * g_zz_xy_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_xy_zz[i] = -2.0 * g_0_xy_xy_zz[i] * b_exp + 4.0 * g_zz_xy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2496-2502)

    #pragma omp simd aligned(g_0_xy_xz_xx, g_0_xy_xz_xy, g_0_xy_xz_xz, g_0_xy_xz_yy, g_0_xy_xz_yz, g_0_xy_xz_zz, g_z_y_0_0_z_x_xz_xx, g_z_y_0_0_z_x_xz_xy, g_z_y_0_0_z_x_xz_xz, g_z_y_0_0_z_x_xz_yy, g_z_y_0_0_z_x_xz_yz, g_z_y_0_0_z_x_xz_zz, g_zz_xy_xz_xx, g_zz_xy_xz_xy, g_zz_xy_xz_xz, g_zz_xy_xz_yy, g_zz_xy_xz_yz, g_zz_xy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_x_xz_xx[i] = -2.0 * g_0_xy_xz_xx[i] * b_exp + 4.0 * g_zz_xy_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_xz_xy[i] = -2.0 * g_0_xy_xz_xy[i] * b_exp + 4.0 * g_zz_xy_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_xz_xz[i] = -2.0 * g_0_xy_xz_xz[i] * b_exp + 4.0 * g_zz_xy_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_xz_yy[i] = -2.0 * g_0_xy_xz_yy[i] * b_exp + 4.0 * g_zz_xy_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_xz_yz[i] = -2.0 * g_0_xy_xz_yz[i] * b_exp + 4.0 * g_zz_xy_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_xz_zz[i] = -2.0 * g_0_xy_xz_zz[i] * b_exp + 4.0 * g_zz_xy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2502-2508)

    #pragma omp simd aligned(g_0_xy_yy_xx, g_0_xy_yy_xy, g_0_xy_yy_xz, g_0_xy_yy_yy, g_0_xy_yy_yz, g_0_xy_yy_zz, g_z_y_0_0_z_x_yy_xx, g_z_y_0_0_z_x_yy_xy, g_z_y_0_0_z_x_yy_xz, g_z_y_0_0_z_x_yy_yy, g_z_y_0_0_z_x_yy_yz, g_z_y_0_0_z_x_yy_zz, g_zz_xy_yy_xx, g_zz_xy_yy_xy, g_zz_xy_yy_xz, g_zz_xy_yy_yy, g_zz_xy_yy_yz, g_zz_xy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_x_yy_xx[i] = -2.0 * g_0_xy_yy_xx[i] * b_exp + 4.0 * g_zz_xy_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_yy_xy[i] = -2.0 * g_0_xy_yy_xy[i] * b_exp + 4.0 * g_zz_xy_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_yy_xz[i] = -2.0 * g_0_xy_yy_xz[i] * b_exp + 4.0 * g_zz_xy_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_yy_yy[i] = -2.0 * g_0_xy_yy_yy[i] * b_exp + 4.0 * g_zz_xy_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_yy_yz[i] = -2.0 * g_0_xy_yy_yz[i] * b_exp + 4.0 * g_zz_xy_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_yy_zz[i] = -2.0 * g_0_xy_yy_zz[i] * b_exp + 4.0 * g_zz_xy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2508-2514)

    #pragma omp simd aligned(g_0_xy_yz_xx, g_0_xy_yz_xy, g_0_xy_yz_xz, g_0_xy_yz_yy, g_0_xy_yz_yz, g_0_xy_yz_zz, g_z_y_0_0_z_x_yz_xx, g_z_y_0_0_z_x_yz_xy, g_z_y_0_0_z_x_yz_xz, g_z_y_0_0_z_x_yz_yy, g_z_y_0_0_z_x_yz_yz, g_z_y_0_0_z_x_yz_zz, g_zz_xy_yz_xx, g_zz_xy_yz_xy, g_zz_xy_yz_xz, g_zz_xy_yz_yy, g_zz_xy_yz_yz, g_zz_xy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_x_yz_xx[i] = -2.0 * g_0_xy_yz_xx[i] * b_exp + 4.0 * g_zz_xy_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_yz_xy[i] = -2.0 * g_0_xy_yz_xy[i] * b_exp + 4.0 * g_zz_xy_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_yz_xz[i] = -2.0 * g_0_xy_yz_xz[i] * b_exp + 4.0 * g_zz_xy_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_yz_yy[i] = -2.0 * g_0_xy_yz_yy[i] * b_exp + 4.0 * g_zz_xy_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_yz_yz[i] = -2.0 * g_0_xy_yz_yz[i] * b_exp + 4.0 * g_zz_xy_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_yz_zz[i] = -2.0 * g_0_xy_yz_zz[i] * b_exp + 4.0 * g_zz_xy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2514-2520)

    #pragma omp simd aligned(g_0_xy_zz_xx, g_0_xy_zz_xy, g_0_xy_zz_xz, g_0_xy_zz_yy, g_0_xy_zz_yz, g_0_xy_zz_zz, g_z_y_0_0_z_x_zz_xx, g_z_y_0_0_z_x_zz_xy, g_z_y_0_0_z_x_zz_xz, g_z_y_0_0_z_x_zz_yy, g_z_y_0_0_z_x_zz_yz, g_z_y_0_0_z_x_zz_zz, g_zz_xy_zz_xx, g_zz_xy_zz_xy, g_zz_xy_zz_xz, g_zz_xy_zz_yy, g_zz_xy_zz_yz, g_zz_xy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_x_zz_xx[i] = -2.0 * g_0_xy_zz_xx[i] * b_exp + 4.0 * g_zz_xy_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_zz_xy[i] = -2.0 * g_0_xy_zz_xy[i] * b_exp + 4.0 * g_zz_xy_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_zz_xz[i] = -2.0 * g_0_xy_zz_xz[i] * b_exp + 4.0 * g_zz_xy_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_zz_yy[i] = -2.0 * g_0_xy_zz_yy[i] * b_exp + 4.0 * g_zz_xy_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_zz_yz[i] = -2.0 * g_0_xy_zz_yz[i] * b_exp + 4.0 * g_zz_xy_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_zz_zz[i] = -2.0 * g_0_xy_zz_zz[i] * b_exp + 4.0 * g_zz_xy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2520-2526)

    #pragma omp simd aligned(g_0_0_xx_xx, g_0_0_xx_xy, g_0_0_xx_xz, g_0_0_xx_yy, g_0_0_xx_yz, g_0_0_xx_zz, g_0_yy_xx_xx, g_0_yy_xx_xy, g_0_yy_xx_xz, g_0_yy_xx_yy, g_0_yy_xx_yz, g_0_yy_xx_zz, g_z_y_0_0_z_y_xx_xx, g_z_y_0_0_z_y_xx_xy, g_z_y_0_0_z_y_xx_xz, g_z_y_0_0_z_y_xx_yy, g_z_y_0_0_z_y_xx_yz, g_z_y_0_0_z_y_xx_zz, g_zz_0_xx_xx, g_zz_0_xx_xy, g_zz_0_xx_xz, g_zz_0_xx_yy, g_zz_0_xx_yz, g_zz_0_xx_zz, g_zz_yy_xx_xx, g_zz_yy_xx_xy, g_zz_yy_xx_xz, g_zz_yy_xx_yy, g_zz_yy_xx_yz, g_zz_yy_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_y_xx_xx[i] = g_0_0_xx_xx[i] - 2.0 * g_0_yy_xx_xx[i] * b_exp - 2.0 * g_zz_0_xx_xx[i] * a_exp + 4.0 * g_zz_yy_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_xx_xy[i] = g_0_0_xx_xy[i] - 2.0 * g_0_yy_xx_xy[i] * b_exp - 2.0 * g_zz_0_xx_xy[i] * a_exp + 4.0 * g_zz_yy_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_xx_xz[i] = g_0_0_xx_xz[i] - 2.0 * g_0_yy_xx_xz[i] * b_exp - 2.0 * g_zz_0_xx_xz[i] * a_exp + 4.0 * g_zz_yy_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_xx_yy[i] = g_0_0_xx_yy[i] - 2.0 * g_0_yy_xx_yy[i] * b_exp - 2.0 * g_zz_0_xx_yy[i] * a_exp + 4.0 * g_zz_yy_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_xx_yz[i] = g_0_0_xx_yz[i] - 2.0 * g_0_yy_xx_yz[i] * b_exp - 2.0 * g_zz_0_xx_yz[i] * a_exp + 4.0 * g_zz_yy_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_xx_zz[i] = g_0_0_xx_zz[i] - 2.0 * g_0_yy_xx_zz[i] * b_exp - 2.0 * g_zz_0_xx_zz[i] * a_exp + 4.0 * g_zz_yy_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2526-2532)

    #pragma omp simd aligned(g_0_0_xy_xx, g_0_0_xy_xy, g_0_0_xy_xz, g_0_0_xy_yy, g_0_0_xy_yz, g_0_0_xy_zz, g_0_yy_xy_xx, g_0_yy_xy_xy, g_0_yy_xy_xz, g_0_yy_xy_yy, g_0_yy_xy_yz, g_0_yy_xy_zz, g_z_y_0_0_z_y_xy_xx, g_z_y_0_0_z_y_xy_xy, g_z_y_0_0_z_y_xy_xz, g_z_y_0_0_z_y_xy_yy, g_z_y_0_0_z_y_xy_yz, g_z_y_0_0_z_y_xy_zz, g_zz_0_xy_xx, g_zz_0_xy_xy, g_zz_0_xy_xz, g_zz_0_xy_yy, g_zz_0_xy_yz, g_zz_0_xy_zz, g_zz_yy_xy_xx, g_zz_yy_xy_xy, g_zz_yy_xy_xz, g_zz_yy_xy_yy, g_zz_yy_xy_yz, g_zz_yy_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_y_xy_xx[i] = g_0_0_xy_xx[i] - 2.0 * g_0_yy_xy_xx[i] * b_exp - 2.0 * g_zz_0_xy_xx[i] * a_exp + 4.0 * g_zz_yy_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_xy_xy[i] = g_0_0_xy_xy[i] - 2.0 * g_0_yy_xy_xy[i] * b_exp - 2.0 * g_zz_0_xy_xy[i] * a_exp + 4.0 * g_zz_yy_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_xy_xz[i] = g_0_0_xy_xz[i] - 2.0 * g_0_yy_xy_xz[i] * b_exp - 2.0 * g_zz_0_xy_xz[i] * a_exp + 4.0 * g_zz_yy_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_xy_yy[i] = g_0_0_xy_yy[i] - 2.0 * g_0_yy_xy_yy[i] * b_exp - 2.0 * g_zz_0_xy_yy[i] * a_exp + 4.0 * g_zz_yy_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_xy_yz[i] = g_0_0_xy_yz[i] - 2.0 * g_0_yy_xy_yz[i] * b_exp - 2.0 * g_zz_0_xy_yz[i] * a_exp + 4.0 * g_zz_yy_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_xy_zz[i] = g_0_0_xy_zz[i] - 2.0 * g_0_yy_xy_zz[i] * b_exp - 2.0 * g_zz_0_xy_zz[i] * a_exp + 4.0 * g_zz_yy_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2532-2538)

    #pragma omp simd aligned(g_0_0_xz_xx, g_0_0_xz_xy, g_0_0_xz_xz, g_0_0_xz_yy, g_0_0_xz_yz, g_0_0_xz_zz, g_0_yy_xz_xx, g_0_yy_xz_xy, g_0_yy_xz_xz, g_0_yy_xz_yy, g_0_yy_xz_yz, g_0_yy_xz_zz, g_z_y_0_0_z_y_xz_xx, g_z_y_0_0_z_y_xz_xy, g_z_y_0_0_z_y_xz_xz, g_z_y_0_0_z_y_xz_yy, g_z_y_0_0_z_y_xz_yz, g_z_y_0_0_z_y_xz_zz, g_zz_0_xz_xx, g_zz_0_xz_xy, g_zz_0_xz_xz, g_zz_0_xz_yy, g_zz_0_xz_yz, g_zz_0_xz_zz, g_zz_yy_xz_xx, g_zz_yy_xz_xy, g_zz_yy_xz_xz, g_zz_yy_xz_yy, g_zz_yy_xz_yz, g_zz_yy_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_y_xz_xx[i] = g_0_0_xz_xx[i] - 2.0 * g_0_yy_xz_xx[i] * b_exp - 2.0 * g_zz_0_xz_xx[i] * a_exp + 4.0 * g_zz_yy_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_xz_xy[i] = g_0_0_xz_xy[i] - 2.0 * g_0_yy_xz_xy[i] * b_exp - 2.0 * g_zz_0_xz_xy[i] * a_exp + 4.0 * g_zz_yy_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_xz_xz[i] = g_0_0_xz_xz[i] - 2.0 * g_0_yy_xz_xz[i] * b_exp - 2.0 * g_zz_0_xz_xz[i] * a_exp + 4.0 * g_zz_yy_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_xz_yy[i] = g_0_0_xz_yy[i] - 2.0 * g_0_yy_xz_yy[i] * b_exp - 2.0 * g_zz_0_xz_yy[i] * a_exp + 4.0 * g_zz_yy_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_xz_yz[i] = g_0_0_xz_yz[i] - 2.0 * g_0_yy_xz_yz[i] * b_exp - 2.0 * g_zz_0_xz_yz[i] * a_exp + 4.0 * g_zz_yy_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_xz_zz[i] = g_0_0_xz_zz[i] - 2.0 * g_0_yy_xz_zz[i] * b_exp - 2.0 * g_zz_0_xz_zz[i] * a_exp + 4.0 * g_zz_yy_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2538-2544)

    #pragma omp simd aligned(g_0_0_yy_xx, g_0_0_yy_xy, g_0_0_yy_xz, g_0_0_yy_yy, g_0_0_yy_yz, g_0_0_yy_zz, g_0_yy_yy_xx, g_0_yy_yy_xy, g_0_yy_yy_xz, g_0_yy_yy_yy, g_0_yy_yy_yz, g_0_yy_yy_zz, g_z_y_0_0_z_y_yy_xx, g_z_y_0_0_z_y_yy_xy, g_z_y_0_0_z_y_yy_xz, g_z_y_0_0_z_y_yy_yy, g_z_y_0_0_z_y_yy_yz, g_z_y_0_0_z_y_yy_zz, g_zz_0_yy_xx, g_zz_0_yy_xy, g_zz_0_yy_xz, g_zz_0_yy_yy, g_zz_0_yy_yz, g_zz_0_yy_zz, g_zz_yy_yy_xx, g_zz_yy_yy_xy, g_zz_yy_yy_xz, g_zz_yy_yy_yy, g_zz_yy_yy_yz, g_zz_yy_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_y_yy_xx[i] = g_0_0_yy_xx[i] - 2.0 * g_0_yy_yy_xx[i] * b_exp - 2.0 * g_zz_0_yy_xx[i] * a_exp + 4.0 * g_zz_yy_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_yy_xy[i] = g_0_0_yy_xy[i] - 2.0 * g_0_yy_yy_xy[i] * b_exp - 2.0 * g_zz_0_yy_xy[i] * a_exp + 4.0 * g_zz_yy_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_yy_xz[i] = g_0_0_yy_xz[i] - 2.0 * g_0_yy_yy_xz[i] * b_exp - 2.0 * g_zz_0_yy_xz[i] * a_exp + 4.0 * g_zz_yy_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_yy_yy[i] = g_0_0_yy_yy[i] - 2.0 * g_0_yy_yy_yy[i] * b_exp - 2.0 * g_zz_0_yy_yy[i] * a_exp + 4.0 * g_zz_yy_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_yy_yz[i] = g_0_0_yy_yz[i] - 2.0 * g_0_yy_yy_yz[i] * b_exp - 2.0 * g_zz_0_yy_yz[i] * a_exp + 4.0 * g_zz_yy_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_yy_zz[i] = g_0_0_yy_zz[i] - 2.0 * g_0_yy_yy_zz[i] * b_exp - 2.0 * g_zz_0_yy_zz[i] * a_exp + 4.0 * g_zz_yy_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2544-2550)

    #pragma omp simd aligned(g_0_0_yz_xx, g_0_0_yz_xy, g_0_0_yz_xz, g_0_0_yz_yy, g_0_0_yz_yz, g_0_0_yz_zz, g_0_yy_yz_xx, g_0_yy_yz_xy, g_0_yy_yz_xz, g_0_yy_yz_yy, g_0_yy_yz_yz, g_0_yy_yz_zz, g_z_y_0_0_z_y_yz_xx, g_z_y_0_0_z_y_yz_xy, g_z_y_0_0_z_y_yz_xz, g_z_y_0_0_z_y_yz_yy, g_z_y_0_0_z_y_yz_yz, g_z_y_0_0_z_y_yz_zz, g_zz_0_yz_xx, g_zz_0_yz_xy, g_zz_0_yz_xz, g_zz_0_yz_yy, g_zz_0_yz_yz, g_zz_0_yz_zz, g_zz_yy_yz_xx, g_zz_yy_yz_xy, g_zz_yy_yz_xz, g_zz_yy_yz_yy, g_zz_yy_yz_yz, g_zz_yy_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_y_yz_xx[i] = g_0_0_yz_xx[i] - 2.0 * g_0_yy_yz_xx[i] * b_exp - 2.0 * g_zz_0_yz_xx[i] * a_exp + 4.0 * g_zz_yy_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_yz_xy[i] = g_0_0_yz_xy[i] - 2.0 * g_0_yy_yz_xy[i] * b_exp - 2.0 * g_zz_0_yz_xy[i] * a_exp + 4.0 * g_zz_yy_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_yz_xz[i] = g_0_0_yz_xz[i] - 2.0 * g_0_yy_yz_xz[i] * b_exp - 2.0 * g_zz_0_yz_xz[i] * a_exp + 4.0 * g_zz_yy_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_yz_yy[i] = g_0_0_yz_yy[i] - 2.0 * g_0_yy_yz_yy[i] * b_exp - 2.0 * g_zz_0_yz_yy[i] * a_exp + 4.0 * g_zz_yy_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_yz_yz[i] = g_0_0_yz_yz[i] - 2.0 * g_0_yy_yz_yz[i] * b_exp - 2.0 * g_zz_0_yz_yz[i] * a_exp + 4.0 * g_zz_yy_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_yz_zz[i] = g_0_0_yz_zz[i] - 2.0 * g_0_yy_yz_zz[i] * b_exp - 2.0 * g_zz_0_yz_zz[i] * a_exp + 4.0 * g_zz_yy_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2550-2556)

    #pragma omp simd aligned(g_0_0_zz_xx, g_0_0_zz_xy, g_0_0_zz_xz, g_0_0_zz_yy, g_0_0_zz_yz, g_0_0_zz_zz, g_0_yy_zz_xx, g_0_yy_zz_xy, g_0_yy_zz_xz, g_0_yy_zz_yy, g_0_yy_zz_yz, g_0_yy_zz_zz, g_z_y_0_0_z_y_zz_xx, g_z_y_0_0_z_y_zz_xy, g_z_y_0_0_z_y_zz_xz, g_z_y_0_0_z_y_zz_yy, g_z_y_0_0_z_y_zz_yz, g_z_y_0_0_z_y_zz_zz, g_zz_0_zz_xx, g_zz_0_zz_xy, g_zz_0_zz_xz, g_zz_0_zz_yy, g_zz_0_zz_yz, g_zz_0_zz_zz, g_zz_yy_zz_xx, g_zz_yy_zz_xy, g_zz_yy_zz_xz, g_zz_yy_zz_yy, g_zz_yy_zz_yz, g_zz_yy_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_y_zz_xx[i] = g_0_0_zz_xx[i] - 2.0 * g_0_yy_zz_xx[i] * b_exp - 2.0 * g_zz_0_zz_xx[i] * a_exp + 4.0 * g_zz_yy_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_zz_xy[i] = g_0_0_zz_xy[i] - 2.0 * g_0_yy_zz_xy[i] * b_exp - 2.0 * g_zz_0_zz_xy[i] * a_exp + 4.0 * g_zz_yy_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_zz_xz[i] = g_0_0_zz_xz[i] - 2.0 * g_0_yy_zz_xz[i] * b_exp - 2.0 * g_zz_0_zz_xz[i] * a_exp + 4.0 * g_zz_yy_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_zz_yy[i] = g_0_0_zz_yy[i] - 2.0 * g_0_yy_zz_yy[i] * b_exp - 2.0 * g_zz_0_zz_yy[i] * a_exp + 4.0 * g_zz_yy_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_zz_yz[i] = g_0_0_zz_yz[i] - 2.0 * g_0_yy_zz_yz[i] * b_exp - 2.0 * g_zz_0_zz_yz[i] * a_exp + 4.0 * g_zz_yy_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_zz_zz[i] = g_0_0_zz_zz[i] - 2.0 * g_0_yy_zz_zz[i] * b_exp - 2.0 * g_zz_0_zz_zz[i] * a_exp + 4.0 * g_zz_yy_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2556-2562)

    #pragma omp simd aligned(g_0_yz_xx_xx, g_0_yz_xx_xy, g_0_yz_xx_xz, g_0_yz_xx_yy, g_0_yz_xx_yz, g_0_yz_xx_zz, g_z_y_0_0_z_z_xx_xx, g_z_y_0_0_z_z_xx_xy, g_z_y_0_0_z_z_xx_xz, g_z_y_0_0_z_z_xx_yy, g_z_y_0_0_z_z_xx_yz, g_z_y_0_0_z_z_xx_zz, g_zz_yz_xx_xx, g_zz_yz_xx_xy, g_zz_yz_xx_xz, g_zz_yz_xx_yy, g_zz_yz_xx_yz, g_zz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_z_xx_xx[i] = -2.0 * g_0_yz_xx_xx[i] * b_exp + 4.0 * g_zz_yz_xx_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_xx_xy[i] = -2.0 * g_0_yz_xx_xy[i] * b_exp + 4.0 * g_zz_yz_xx_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_xx_xz[i] = -2.0 * g_0_yz_xx_xz[i] * b_exp + 4.0 * g_zz_yz_xx_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_xx_yy[i] = -2.0 * g_0_yz_xx_yy[i] * b_exp + 4.0 * g_zz_yz_xx_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_xx_yz[i] = -2.0 * g_0_yz_xx_yz[i] * b_exp + 4.0 * g_zz_yz_xx_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_xx_zz[i] = -2.0 * g_0_yz_xx_zz[i] * b_exp + 4.0 * g_zz_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2562-2568)

    #pragma omp simd aligned(g_0_yz_xy_xx, g_0_yz_xy_xy, g_0_yz_xy_xz, g_0_yz_xy_yy, g_0_yz_xy_yz, g_0_yz_xy_zz, g_z_y_0_0_z_z_xy_xx, g_z_y_0_0_z_z_xy_xy, g_z_y_0_0_z_z_xy_xz, g_z_y_0_0_z_z_xy_yy, g_z_y_0_0_z_z_xy_yz, g_z_y_0_0_z_z_xy_zz, g_zz_yz_xy_xx, g_zz_yz_xy_xy, g_zz_yz_xy_xz, g_zz_yz_xy_yy, g_zz_yz_xy_yz, g_zz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_z_xy_xx[i] = -2.0 * g_0_yz_xy_xx[i] * b_exp + 4.0 * g_zz_yz_xy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_xy_xy[i] = -2.0 * g_0_yz_xy_xy[i] * b_exp + 4.0 * g_zz_yz_xy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_xy_xz[i] = -2.0 * g_0_yz_xy_xz[i] * b_exp + 4.0 * g_zz_yz_xy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_xy_yy[i] = -2.0 * g_0_yz_xy_yy[i] * b_exp + 4.0 * g_zz_yz_xy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_xy_yz[i] = -2.0 * g_0_yz_xy_yz[i] * b_exp + 4.0 * g_zz_yz_xy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_xy_zz[i] = -2.0 * g_0_yz_xy_zz[i] * b_exp + 4.0 * g_zz_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2568-2574)

    #pragma omp simd aligned(g_0_yz_xz_xx, g_0_yz_xz_xy, g_0_yz_xz_xz, g_0_yz_xz_yy, g_0_yz_xz_yz, g_0_yz_xz_zz, g_z_y_0_0_z_z_xz_xx, g_z_y_0_0_z_z_xz_xy, g_z_y_0_0_z_z_xz_xz, g_z_y_0_0_z_z_xz_yy, g_z_y_0_0_z_z_xz_yz, g_z_y_0_0_z_z_xz_zz, g_zz_yz_xz_xx, g_zz_yz_xz_xy, g_zz_yz_xz_xz, g_zz_yz_xz_yy, g_zz_yz_xz_yz, g_zz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_z_xz_xx[i] = -2.0 * g_0_yz_xz_xx[i] * b_exp + 4.0 * g_zz_yz_xz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_xz_xy[i] = -2.0 * g_0_yz_xz_xy[i] * b_exp + 4.0 * g_zz_yz_xz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_xz_xz[i] = -2.0 * g_0_yz_xz_xz[i] * b_exp + 4.0 * g_zz_yz_xz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_xz_yy[i] = -2.0 * g_0_yz_xz_yy[i] * b_exp + 4.0 * g_zz_yz_xz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_xz_yz[i] = -2.0 * g_0_yz_xz_yz[i] * b_exp + 4.0 * g_zz_yz_xz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_xz_zz[i] = -2.0 * g_0_yz_xz_zz[i] * b_exp + 4.0 * g_zz_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2574-2580)

    #pragma omp simd aligned(g_0_yz_yy_xx, g_0_yz_yy_xy, g_0_yz_yy_xz, g_0_yz_yy_yy, g_0_yz_yy_yz, g_0_yz_yy_zz, g_z_y_0_0_z_z_yy_xx, g_z_y_0_0_z_z_yy_xy, g_z_y_0_0_z_z_yy_xz, g_z_y_0_0_z_z_yy_yy, g_z_y_0_0_z_z_yy_yz, g_z_y_0_0_z_z_yy_zz, g_zz_yz_yy_xx, g_zz_yz_yy_xy, g_zz_yz_yy_xz, g_zz_yz_yy_yy, g_zz_yz_yy_yz, g_zz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_z_yy_xx[i] = -2.0 * g_0_yz_yy_xx[i] * b_exp + 4.0 * g_zz_yz_yy_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_yy_xy[i] = -2.0 * g_0_yz_yy_xy[i] * b_exp + 4.0 * g_zz_yz_yy_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_yy_xz[i] = -2.0 * g_0_yz_yy_xz[i] * b_exp + 4.0 * g_zz_yz_yy_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_yy_yy[i] = -2.0 * g_0_yz_yy_yy[i] * b_exp + 4.0 * g_zz_yz_yy_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_yy_yz[i] = -2.0 * g_0_yz_yy_yz[i] * b_exp + 4.0 * g_zz_yz_yy_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_yy_zz[i] = -2.0 * g_0_yz_yy_zz[i] * b_exp + 4.0 * g_zz_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2580-2586)

    #pragma omp simd aligned(g_0_yz_yz_xx, g_0_yz_yz_xy, g_0_yz_yz_xz, g_0_yz_yz_yy, g_0_yz_yz_yz, g_0_yz_yz_zz, g_z_y_0_0_z_z_yz_xx, g_z_y_0_0_z_z_yz_xy, g_z_y_0_0_z_z_yz_xz, g_z_y_0_0_z_z_yz_yy, g_z_y_0_0_z_z_yz_yz, g_z_y_0_0_z_z_yz_zz, g_zz_yz_yz_xx, g_zz_yz_yz_xy, g_zz_yz_yz_xz, g_zz_yz_yz_yy, g_zz_yz_yz_yz, g_zz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_z_yz_xx[i] = -2.0 * g_0_yz_yz_xx[i] * b_exp + 4.0 * g_zz_yz_yz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_yz_xy[i] = -2.0 * g_0_yz_yz_xy[i] * b_exp + 4.0 * g_zz_yz_yz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_yz_xz[i] = -2.0 * g_0_yz_yz_xz[i] * b_exp + 4.0 * g_zz_yz_yz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_yz_yy[i] = -2.0 * g_0_yz_yz_yy[i] * b_exp + 4.0 * g_zz_yz_yz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_yz_yz[i] = -2.0 * g_0_yz_yz_yz[i] * b_exp + 4.0 * g_zz_yz_yz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_yz_zz[i] = -2.0 * g_0_yz_yz_zz[i] * b_exp + 4.0 * g_zz_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2586-2592)

    #pragma omp simd aligned(g_0_yz_zz_xx, g_0_yz_zz_xy, g_0_yz_zz_xz, g_0_yz_zz_yy, g_0_yz_zz_yz, g_0_yz_zz_zz, g_z_y_0_0_z_z_zz_xx, g_z_y_0_0_z_z_zz_xy, g_z_y_0_0_z_z_zz_xz, g_z_y_0_0_z_z_zz_yy, g_z_y_0_0_z_z_zz_yz, g_z_y_0_0_z_z_zz_zz, g_zz_yz_zz_xx, g_zz_yz_zz_xy, g_zz_yz_zz_xz, g_zz_yz_zz_yy, g_zz_yz_zz_yz, g_zz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_z_zz_xx[i] = -2.0 * g_0_yz_zz_xx[i] * b_exp + 4.0 * g_zz_yz_zz_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_zz_xy[i] = -2.0 * g_0_yz_zz_xy[i] * b_exp + 4.0 * g_zz_yz_zz_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_zz_xz[i] = -2.0 * g_0_yz_zz_xz[i] * b_exp + 4.0 * g_zz_yz_zz_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_zz_yy[i] = -2.0 * g_0_yz_zz_yy[i] * b_exp + 4.0 * g_zz_yz_zz_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_zz_yz[i] = -2.0 * g_0_yz_zz_yz[i] * b_exp + 4.0 * g_zz_yz_zz_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_zz_zz[i] = -2.0 * g_0_yz_zz_zz[i] * b_exp + 4.0 * g_zz_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2592-2598)

    #pragma omp simd aligned(g_xz_xz_xx_xx, g_xz_xz_xx_xy, g_xz_xz_xx_xz, g_xz_xz_xx_yy, g_xz_xz_xx_yz, g_xz_xz_xx_zz, g_z_z_0_0_x_x_xx_xx, g_z_z_0_0_x_x_xx_xy, g_z_z_0_0_x_x_xx_xz, g_z_z_0_0_x_x_xx_yy, g_z_z_0_0_x_x_xx_yz, g_z_z_0_0_x_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_x_xx_xx[i] = 4.0 * g_xz_xz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_xx_xy[i] = 4.0 * g_xz_xz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_xx_xz[i] = 4.0 * g_xz_xz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_xx_yy[i] = 4.0 * g_xz_xz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_xx_yz[i] = 4.0 * g_xz_xz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_xx_zz[i] = 4.0 * g_xz_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2598-2604)

    #pragma omp simd aligned(g_xz_xz_xy_xx, g_xz_xz_xy_xy, g_xz_xz_xy_xz, g_xz_xz_xy_yy, g_xz_xz_xy_yz, g_xz_xz_xy_zz, g_z_z_0_0_x_x_xy_xx, g_z_z_0_0_x_x_xy_xy, g_z_z_0_0_x_x_xy_xz, g_z_z_0_0_x_x_xy_yy, g_z_z_0_0_x_x_xy_yz, g_z_z_0_0_x_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_x_xy_xx[i] = 4.0 * g_xz_xz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_xy_xy[i] = 4.0 * g_xz_xz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_xy_xz[i] = 4.0 * g_xz_xz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_xy_yy[i] = 4.0 * g_xz_xz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_xy_yz[i] = 4.0 * g_xz_xz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_xy_zz[i] = 4.0 * g_xz_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2604-2610)

    #pragma omp simd aligned(g_xz_xz_xz_xx, g_xz_xz_xz_xy, g_xz_xz_xz_xz, g_xz_xz_xz_yy, g_xz_xz_xz_yz, g_xz_xz_xz_zz, g_z_z_0_0_x_x_xz_xx, g_z_z_0_0_x_x_xz_xy, g_z_z_0_0_x_x_xz_xz, g_z_z_0_0_x_x_xz_yy, g_z_z_0_0_x_x_xz_yz, g_z_z_0_0_x_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_x_xz_xx[i] = 4.0 * g_xz_xz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_xz_xy[i] = 4.0 * g_xz_xz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_xz_xz[i] = 4.0 * g_xz_xz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_xz_yy[i] = 4.0 * g_xz_xz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_xz_yz[i] = 4.0 * g_xz_xz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_xz_zz[i] = 4.0 * g_xz_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2610-2616)

    #pragma omp simd aligned(g_xz_xz_yy_xx, g_xz_xz_yy_xy, g_xz_xz_yy_xz, g_xz_xz_yy_yy, g_xz_xz_yy_yz, g_xz_xz_yy_zz, g_z_z_0_0_x_x_yy_xx, g_z_z_0_0_x_x_yy_xy, g_z_z_0_0_x_x_yy_xz, g_z_z_0_0_x_x_yy_yy, g_z_z_0_0_x_x_yy_yz, g_z_z_0_0_x_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_x_yy_xx[i] = 4.0 * g_xz_xz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_yy_xy[i] = 4.0 * g_xz_xz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_yy_xz[i] = 4.0 * g_xz_xz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_yy_yy[i] = 4.0 * g_xz_xz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_yy_yz[i] = 4.0 * g_xz_xz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_yy_zz[i] = 4.0 * g_xz_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2616-2622)

    #pragma omp simd aligned(g_xz_xz_yz_xx, g_xz_xz_yz_xy, g_xz_xz_yz_xz, g_xz_xz_yz_yy, g_xz_xz_yz_yz, g_xz_xz_yz_zz, g_z_z_0_0_x_x_yz_xx, g_z_z_0_0_x_x_yz_xy, g_z_z_0_0_x_x_yz_xz, g_z_z_0_0_x_x_yz_yy, g_z_z_0_0_x_x_yz_yz, g_z_z_0_0_x_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_x_yz_xx[i] = 4.0 * g_xz_xz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_yz_xy[i] = 4.0 * g_xz_xz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_yz_xz[i] = 4.0 * g_xz_xz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_yz_yy[i] = 4.0 * g_xz_xz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_yz_yz[i] = 4.0 * g_xz_xz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_yz_zz[i] = 4.0 * g_xz_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2622-2628)

    #pragma omp simd aligned(g_xz_xz_zz_xx, g_xz_xz_zz_xy, g_xz_xz_zz_xz, g_xz_xz_zz_yy, g_xz_xz_zz_yz, g_xz_xz_zz_zz, g_z_z_0_0_x_x_zz_xx, g_z_z_0_0_x_x_zz_xy, g_z_z_0_0_x_x_zz_xz, g_z_z_0_0_x_x_zz_yy, g_z_z_0_0_x_x_zz_yz, g_z_z_0_0_x_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_x_zz_xx[i] = 4.0 * g_xz_xz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_zz_xy[i] = 4.0 * g_xz_xz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_zz_xz[i] = 4.0 * g_xz_xz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_zz_yy[i] = 4.0 * g_xz_xz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_zz_yz[i] = 4.0 * g_xz_xz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_zz_zz[i] = 4.0 * g_xz_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2628-2634)

    #pragma omp simd aligned(g_xz_yz_xx_xx, g_xz_yz_xx_xy, g_xz_yz_xx_xz, g_xz_yz_xx_yy, g_xz_yz_xx_yz, g_xz_yz_xx_zz, g_z_z_0_0_x_y_xx_xx, g_z_z_0_0_x_y_xx_xy, g_z_z_0_0_x_y_xx_xz, g_z_z_0_0_x_y_xx_yy, g_z_z_0_0_x_y_xx_yz, g_z_z_0_0_x_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_y_xx_xx[i] = 4.0 * g_xz_yz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_xx_xy[i] = 4.0 * g_xz_yz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_xx_xz[i] = 4.0 * g_xz_yz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_xx_yy[i] = 4.0 * g_xz_yz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_xx_yz[i] = 4.0 * g_xz_yz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_xx_zz[i] = 4.0 * g_xz_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2634-2640)

    #pragma omp simd aligned(g_xz_yz_xy_xx, g_xz_yz_xy_xy, g_xz_yz_xy_xz, g_xz_yz_xy_yy, g_xz_yz_xy_yz, g_xz_yz_xy_zz, g_z_z_0_0_x_y_xy_xx, g_z_z_0_0_x_y_xy_xy, g_z_z_0_0_x_y_xy_xz, g_z_z_0_0_x_y_xy_yy, g_z_z_0_0_x_y_xy_yz, g_z_z_0_0_x_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_y_xy_xx[i] = 4.0 * g_xz_yz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_xy_xy[i] = 4.0 * g_xz_yz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_xy_xz[i] = 4.0 * g_xz_yz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_xy_yy[i] = 4.0 * g_xz_yz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_xy_yz[i] = 4.0 * g_xz_yz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_xy_zz[i] = 4.0 * g_xz_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2640-2646)

    #pragma omp simd aligned(g_xz_yz_xz_xx, g_xz_yz_xz_xy, g_xz_yz_xz_xz, g_xz_yz_xz_yy, g_xz_yz_xz_yz, g_xz_yz_xz_zz, g_z_z_0_0_x_y_xz_xx, g_z_z_0_0_x_y_xz_xy, g_z_z_0_0_x_y_xz_xz, g_z_z_0_0_x_y_xz_yy, g_z_z_0_0_x_y_xz_yz, g_z_z_0_0_x_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_y_xz_xx[i] = 4.0 * g_xz_yz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_xz_xy[i] = 4.0 * g_xz_yz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_xz_xz[i] = 4.0 * g_xz_yz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_xz_yy[i] = 4.0 * g_xz_yz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_xz_yz[i] = 4.0 * g_xz_yz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_xz_zz[i] = 4.0 * g_xz_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2646-2652)

    #pragma omp simd aligned(g_xz_yz_yy_xx, g_xz_yz_yy_xy, g_xz_yz_yy_xz, g_xz_yz_yy_yy, g_xz_yz_yy_yz, g_xz_yz_yy_zz, g_z_z_0_0_x_y_yy_xx, g_z_z_0_0_x_y_yy_xy, g_z_z_0_0_x_y_yy_xz, g_z_z_0_0_x_y_yy_yy, g_z_z_0_0_x_y_yy_yz, g_z_z_0_0_x_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_y_yy_xx[i] = 4.0 * g_xz_yz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_yy_xy[i] = 4.0 * g_xz_yz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_yy_xz[i] = 4.0 * g_xz_yz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_yy_yy[i] = 4.0 * g_xz_yz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_yy_yz[i] = 4.0 * g_xz_yz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_yy_zz[i] = 4.0 * g_xz_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2652-2658)

    #pragma omp simd aligned(g_xz_yz_yz_xx, g_xz_yz_yz_xy, g_xz_yz_yz_xz, g_xz_yz_yz_yy, g_xz_yz_yz_yz, g_xz_yz_yz_zz, g_z_z_0_0_x_y_yz_xx, g_z_z_0_0_x_y_yz_xy, g_z_z_0_0_x_y_yz_xz, g_z_z_0_0_x_y_yz_yy, g_z_z_0_0_x_y_yz_yz, g_z_z_0_0_x_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_y_yz_xx[i] = 4.0 * g_xz_yz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_yz_xy[i] = 4.0 * g_xz_yz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_yz_xz[i] = 4.0 * g_xz_yz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_yz_yy[i] = 4.0 * g_xz_yz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_yz_yz[i] = 4.0 * g_xz_yz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_yz_zz[i] = 4.0 * g_xz_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2658-2664)

    #pragma omp simd aligned(g_xz_yz_zz_xx, g_xz_yz_zz_xy, g_xz_yz_zz_xz, g_xz_yz_zz_yy, g_xz_yz_zz_yz, g_xz_yz_zz_zz, g_z_z_0_0_x_y_zz_xx, g_z_z_0_0_x_y_zz_xy, g_z_z_0_0_x_y_zz_xz, g_z_z_0_0_x_y_zz_yy, g_z_z_0_0_x_y_zz_yz, g_z_z_0_0_x_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_y_zz_xx[i] = 4.0 * g_xz_yz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_zz_xy[i] = 4.0 * g_xz_yz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_zz_xz[i] = 4.0 * g_xz_yz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_zz_yy[i] = 4.0 * g_xz_yz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_zz_yz[i] = 4.0 * g_xz_yz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_zz_zz[i] = 4.0 * g_xz_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2664-2670)

    #pragma omp simd aligned(g_xz_0_xx_xx, g_xz_0_xx_xy, g_xz_0_xx_xz, g_xz_0_xx_yy, g_xz_0_xx_yz, g_xz_0_xx_zz, g_xz_zz_xx_xx, g_xz_zz_xx_xy, g_xz_zz_xx_xz, g_xz_zz_xx_yy, g_xz_zz_xx_yz, g_xz_zz_xx_zz, g_z_z_0_0_x_z_xx_xx, g_z_z_0_0_x_z_xx_xy, g_z_z_0_0_x_z_xx_xz, g_z_z_0_0_x_z_xx_yy, g_z_z_0_0_x_z_xx_yz, g_z_z_0_0_x_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_z_xx_xx[i] = -2.0 * g_xz_0_xx_xx[i] * a_exp + 4.0 * g_xz_zz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_xx_xy[i] = -2.0 * g_xz_0_xx_xy[i] * a_exp + 4.0 * g_xz_zz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_xx_xz[i] = -2.0 * g_xz_0_xx_xz[i] * a_exp + 4.0 * g_xz_zz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_xx_yy[i] = -2.0 * g_xz_0_xx_yy[i] * a_exp + 4.0 * g_xz_zz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_xx_yz[i] = -2.0 * g_xz_0_xx_yz[i] * a_exp + 4.0 * g_xz_zz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_xx_zz[i] = -2.0 * g_xz_0_xx_zz[i] * a_exp + 4.0 * g_xz_zz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2670-2676)

    #pragma omp simd aligned(g_xz_0_xy_xx, g_xz_0_xy_xy, g_xz_0_xy_xz, g_xz_0_xy_yy, g_xz_0_xy_yz, g_xz_0_xy_zz, g_xz_zz_xy_xx, g_xz_zz_xy_xy, g_xz_zz_xy_xz, g_xz_zz_xy_yy, g_xz_zz_xy_yz, g_xz_zz_xy_zz, g_z_z_0_0_x_z_xy_xx, g_z_z_0_0_x_z_xy_xy, g_z_z_0_0_x_z_xy_xz, g_z_z_0_0_x_z_xy_yy, g_z_z_0_0_x_z_xy_yz, g_z_z_0_0_x_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_z_xy_xx[i] = -2.0 * g_xz_0_xy_xx[i] * a_exp + 4.0 * g_xz_zz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_xy_xy[i] = -2.0 * g_xz_0_xy_xy[i] * a_exp + 4.0 * g_xz_zz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_xy_xz[i] = -2.0 * g_xz_0_xy_xz[i] * a_exp + 4.0 * g_xz_zz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_xy_yy[i] = -2.0 * g_xz_0_xy_yy[i] * a_exp + 4.0 * g_xz_zz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_xy_yz[i] = -2.0 * g_xz_0_xy_yz[i] * a_exp + 4.0 * g_xz_zz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_xy_zz[i] = -2.0 * g_xz_0_xy_zz[i] * a_exp + 4.0 * g_xz_zz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2676-2682)

    #pragma omp simd aligned(g_xz_0_xz_xx, g_xz_0_xz_xy, g_xz_0_xz_xz, g_xz_0_xz_yy, g_xz_0_xz_yz, g_xz_0_xz_zz, g_xz_zz_xz_xx, g_xz_zz_xz_xy, g_xz_zz_xz_xz, g_xz_zz_xz_yy, g_xz_zz_xz_yz, g_xz_zz_xz_zz, g_z_z_0_0_x_z_xz_xx, g_z_z_0_0_x_z_xz_xy, g_z_z_0_0_x_z_xz_xz, g_z_z_0_0_x_z_xz_yy, g_z_z_0_0_x_z_xz_yz, g_z_z_0_0_x_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_z_xz_xx[i] = -2.0 * g_xz_0_xz_xx[i] * a_exp + 4.0 * g_xz_zz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_xz_xy[i] = -2.0 * g_xz_0_xz_xy[i] * a_exp + 4.0 * g_xz_zz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_xz_xz[i] = -2.0 * g_xz_0_xz_xz[i] * a_exp + 4.0 * g_xz_zz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_xz_yy[i] = -2.0 * g_xz_0_xz_yy[i] * a_exp + 4.0 * g_xz_zz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_xz_yz[i] = -2.0 * g_xz_0_xz_yz[i] * a_exp + 4.0 * g_xz_zz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_xz_zz[i] = -2.0 * g_xz_0_xz_zz[i] * a_exp + 4.0 * g_xz_zz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2682-2688)

    #pragma omp simd aligned(g_xz_0_yy_xx, g_xz_0_yy_xy, g_xz_0_yy_xz, g_xz_0_yy_yy, g_xz_0_yy_yz, g_xz_0_yy_zz, g_xz_zz_yy_xx, g_xz_zz_yy_xy, g_xz_zz_yy_xz, g_xz_zz_yy_yy, g_xz_zz_yy_yz, g_xz_zz_yy_zz, g_z_z_0_0_x_z_yy_xx, g_z_z_0_0_x_z_yy_xy, g_z_z_0_0_x_z_yy_xz, g_z_z_0_0_x_z_yy_yy, g_z_z_0_0_x_z_yy_yz, g_z_z_0_0_x_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_z_yy_xx[i] = -2.0 * g_xz_0_yy_xx[i] * a_exp + 4.0 * g_xz_zz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_yy_xy[i] = -2.0 * g_xz_0_yy_xy[i] * a_exp + 4.0 * g_xz_zz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_yy_xz[i] = -2.0 * g_xz_0_yy_xz[i] * a_exp + 4.0 * g_xz_zz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_yy_yy[i] = -2.0 * g_xz_0_yy_yy[i] * a_exp + 4.0 * g_xz_zz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_yy_yz[i] = -2.0 * g_xz_0_yy_yz[i] * a_exp + 4.0 * g_xz_zz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_yy_zz[i] = -2.0 * g_xz_0_yy_zz[i] * a_exp + 4.0 * g_xz_zz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2688-2694)

    #pragma omp simd aligned(g_xz_0_yz_xx, g_xz_0_yz_xy, g_xz_0_yz_xz, g_xz_0_yz_yy, g_xz_0_yz_yz, g_xz_0_yz_zz, g_xz_zz_yz_xx, g_xz_zz_yz_xy, g_xz_zz_yz_xz, g_xz_zz_yz_yy, g_xz_zz_yz_yz, g_xz_zz_yz_zz, g_z_z_0_0_x_z_yz_xx, g_z_z_0_0_x_z_yz_xy, g_z_z_0_0_x_z_yz_xz, g_z_z_0_0_x_z_yz_yy, g_z_z_0_0_x_z_yz_yz, g_z_z_0_0_x_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_z_yz_xx[i] = -2.0 * g_xz_0_yz_xx[i] * a_exp + 4.0 * g_xz_zz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_yz_xy[i] = -2.0 * g_xz_0_yz_xy[i] * a_exp + 4.0 * g_xz_zz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_yz_xz[i] = -2.0 * g_xz_0_yz_xz[i] * a_exp + 4.0 * g_xz_zz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_yz_yy[i] = -2.0 * g_xz_0_yz_yy[i] * a_exp + 4.0 * g_xz_zz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_yz_yz[i] = -2.0 * g_xz_0_yz_yz[i] * a_exp + 4.0 * g_xz_zz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_yz_zz[i] = -2.0 * g_xz_0_yz_zz[i] * a_exp + 4.0 * g_xz_zz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2694-2700)

    #pragma omp simd aligned(g_xz_0_zz_xx, g_xz_0_zz_xy, g_xz_0_zz_xz, g_xz_0_zz_yy, g_xz_0_zz_yz, g_xz_0_zz_zz, g_xz_zz_zz_xx, g_xz_zz_zz_xy, g_xz_zz_zz_xz, g_xz_zz_zz_yy, g_xz_zz_zz_yz, g_xz_zz_zz_zz, g_z_z_0_0_x_z_zz_xx, g_z_z_0_0_x_z_zz_xy, g_z_z_0_0_x_z_zz_xz, g_z_z_0_0_x_z_zz_yy, g_z_z_0_0_x_z_zz_yz, g_z_z_0_0_x_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_z_zz_xx[i] = -2.0 * g_xz_0_zz_xx[i] * a_exp + 4.0 * g_xz_zz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_zz_xy[i] = -2.0 * g_xz_0_zz_xy[i] * a_exp + 4.0 * g_xz_zz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_zz_xz[i] = -2.0 * g_xz_0_zz_xz[i] * a_exp + 4.0 * g_xz_zz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_zz_yy[i] = -2.0 * g_xz_0_zz_yy[i] * a_exp + 4.0 * g_xz_zz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_zz_yz[i] = -2.0 * g_xz_0_zz_yz[i] * a_exp + 4.0 * g_xz_zz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_zz_zz[i] = -2.0 * g_xz_0_zz_zz[i] * a_exp + 4.0 * g_xz_zz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2700-2706)

    #pragma omp simd aligned(g_yz_xz_xx_xx, g_yz_xz_xx_xy, g_yz_xz_xx_xz, g_yz_xz_xx_yy, g_yz_xz_xx_yz, g_yz_xz_xx_zz, g_z_z_0_0_y_x_xx_xx, g_z_z_0_0_y_x_xx_xy, g_z_z_0_0_y_x_xx_xz, g_z_z_0_0_y_x_xx_yy, g_z_z_0_0_y_x_xx_yz, g_z_z_0_0_y_x_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_x_xx_xx[i] = 4.0 * g_yz_xz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_xx_xy[i] = 4.0 * g_yz_xz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_xx_xz[i] = 4.0 * g_yz_xz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_xx_yy[i] = 4.0 * g_yz_xz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_xx_yz[i] = 4.0 * g_yz_xz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_xx_zz[i] = 4.0 * g_yz_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2706-2712)

    #pragma omp simd aligned(g_yz_xz_xy_xx, g_yz_xz_xy_xy, g_yz_xz_xy_xz, g_yz_xz_xy_yy, g_yz_xz_xy_yz, g_yz_xz_xy_zz, g_z_z_0_0_y_x_xy_xx, g_z_z_0_0_y_x_xy_xy, g_z_z_0_0_y_x_xy_xz, g_z_z_0_0_y_x_xy_yy, g_z_z_0_0_y_x_xy_yz, g_z_z_0_0_y_x_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_x_xy_xx[i] = 4.0 * g_yz_xz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_xy_xy[i] = 4.0 * g_yz_xz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_xy_xz[i] = 4.0 * g_yz_xz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_xy_yy[i] = 4.0 * g_yz_xz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_xy_yz[i] = 4.0 * g_yz_xz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_xy_zz[i] = 4.0 * g_yz_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2712-2718)

    #pragma omp simd aligned(g_yz_xz_xz_xx, g_yz_xz_xz_xy, g_yz_xz_xz_xz, g_yz_xz_xz_yy, g_yz_xz_xz_yz, g_yz_xz_xz_zz, g_z_z_0_0_y_x_xz_xx, g_z_z_0_0_y_x_xz_xy, g_z_z_0_0_y_x_xz_xz, g_z_z_0_0_y_x_xz_yy, g_z_z_0_0_y_x_xz_yz, g_z_z_0_0_y_x_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_x_xz_xx[i] = 4.0 * g_yz_xz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_xz_xy[i] = 4.0 * g_yz_xz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_xz_xz[i] = 4.0 * g_yz_xz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_xz_yy[i] = 4.0 * g_yz_xz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_xz_yz[i] = 4.0 * g_yz_xz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_xz_zz[i] = 4.0 * g_yz_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2718-2724)

    #pragma omp simd aligned(g_yz_xz_yy_xx, g_yz_xz_yy_xy, g_yz_xz_yy_xz, g_yz_xz_yy_yy, g_yz_xz_yy_yz, g_yz_xz_yy_zz, g_z_z_0_0_y_x_yy_xx, g_z_z_0_0_y_x_yy_xy, g_z_z_0_0_y_x_yy_xz, g_z_z_0_0_y_x_yy_yy, g_z_z_0_0_y_x_yy_yz, g_z_z_0_0_y_x_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_x_yy_xx[i] = 4.0 * g_yz_xz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_yy_xy[i] = 4.0 * g_yz_xz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_yy_xz[i] = 4.0 * g_yz_xz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_yy_yy[i] = 4.0 * g_yz_xz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_yy_yz[i] = 4.0 * g_yz_xz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_yy_zz[i] = 4.0 * g_yz_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2724-2730)

    #pragma omp simd aligned(g_yz_xz_yz_xx, g_yz_xz_yz_xy, g_yz_xz_yz_xz, g_yz_xz_yz_yy, g_yz_xz_yz_yz, g_yz_xz_yz_zz, g_z_z_0_0_y_x_yz_xx, g_z_z_0_0_y_x_yz_xy, g_z_z_0_0_y_x_yz_xz, g_z_z_0_0_y_x_yz_yy, g_z_z_0_0_y_x_yz_yz, g_z_z_0_0_y_x_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_x_yz_xx[i] = 4.0 * g_yz_xz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_yz_xy[i] = 4.0 * g_yz_xz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_yz_xz[i] = 4.0 * g_yz_xz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_yz_yy[i] = 4.0 * g_yz_xz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_yz_yz[i] = 4.0 * g_yz_xz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_yz_zz[i] = 4.0 * g_yz_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2730-2736)

    #pragma omp simd aligned(g_yz_xz_zz_xx, g_yz_xz_zz_xy, g_yz_xz_zz_xz, g_yz_xz_zz_yy, g_yz_xz_zz_yz, g_yz_xz_zz_zz, g_z_z_0_0_y_x_zz_xx, g_z_z_0_0_y_x_zz_xy, g_z_z_0_0_y_x_zz_xz, g_z_z_0_0_y_x_zz_yy, g_z_z_0_0_y_x_zz_yz, g_z_z_0_0_y_x_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_x_zz_xx[i] = 4.0 * g_yz_xz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_zz_xy[i] = 4.0 * g_yz_xz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_zz_xz[i] = 4.0 * g_yz_xz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_zz_yy[i] = 4.0 * g_yz_xz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_zz_yz[i] = 4.0 * g_yz_xz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_zz_zz[i] = 4.0 * g_yz_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2736-2742)

    #pragma omp simd aligned(g_yz_yz_xx_xx, g_yz_yz_xx_xy, g_yz_yz_xx_xz, g_yz_yz_xx_yy, g_yz_yz_xx_yz, g_yz_yz_xx_zz, g_z_z_0_0_y_y_xx_xx, g_z_z_0_0_y_y_xx_xy, g_z_z_0_0_y_y_xx_xz, g_z_z_0_0_y_y_xx_yy, g_z_z_0_0_y_y_xx_yz, g_z_z_0_0_y_y_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_y_xx_xx[i] = 4.0 * g_yz_yz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_xx_xy[i] = 4.0 * g_yz_yz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_xx_xz[i] = 4.0 * g_yz_yz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_xx_yy[i] = 4.0 * g_yz_yz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_xx_yz[i] = 4.0 * g_yz_yz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_xx_zz[i] = 4.0 * g_yz_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2742-2748)

    #pragma omp simd aligned(g_yz_yz_xy_xx, g_yz_yz_xy_xy, g_yz_yz_xy_xz, g_yz_yz_xy_yy, g_yz_yz_xy_yz, g_yz_yz_xy_zz, g_z_z_0_0_y_y_xy_xx, g_z_z_0_0_y_y_xy_xy, g_z_z_0_0_y_y_xy_xz, g_z_z_0_0_y_y_xy_yy, g_z_z_0_0_y_y_xy_yz, g_z_z_0_0_y_y_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_y_xy_xx[i] = 4.0 * g_yz_yz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_xy_xy[i] = 4.0 * g_yz_yz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_xy_xz[i] = 4.0 * g_yz_yz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_xy_yy[i] = 4.0 * g_yz_yz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_xy_yz[i] = 4.0 * g_yz_yz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_xy_zz[i] = 4.0 * g_yz_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2748-2754)

    #pragma omp simd aligned(g_yz_yz_xz_xx, g_yz_yz_xz_xy, g_yz_yz_xz_xz, g_yz_yz_xz_yy, g_yz_yz_xz_yz, g_yz_yz_xz_zz, g_z_z_0_0_y_y_xz_xx, g_z_z_0_0_y_y_xz_xy, g_z_z_0_0_y_y_xz_xz, g_z_z_0_0_y_y_xz_yy, g_z_z_0_0_y_y_xz_yz, g_z_z_0_0_y_y_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_y_xz_xx[i] = 4.0 * g_yz_yz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_xz_xy[i] = 4.0 * g_yz_yz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_xz_xz[i] = 4.0 * g_yz_yz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_xz_yy[i] = 4.0 * g_yz_yz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_xz_yz[i] = 4.0 * g_yz_yz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_xz_zz[i] = 4.0 * g_yz_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2754-2760)

    #pragma omp simd aligned(g_yz_yz_yy_xx, g_yz_yz_yy_xy, g_yz_yz_yy_xz, g_yz_yz_yy_yy, g_yz_yz_yy_yz, g_yz_yz_yy_zz, g_z_z_0_0_y_y_yy_xx, g_z_z_0_0_y_y_yy_xy, g_z_z_0_0_y_y_yy_xz, g_z_z_0_0_y_y_yy_yy, g_z_z_0_0_y_y_yy_yz, g_z_z_0_0_y_y_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_y_yy_xx[i] = 4.0 * g_yz_yz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_yy_xy[i] = 4.0 * g_yz_yz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_yy_xz[i] = 4.0 * g_yz_yz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_yy_yy[i] = 4.0 * g_yz_yz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_yy_yz[i] = 4.0 * g_yz_yz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_yy_zz[i] = 4.0 * g_yz_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2760-2766)

    #pragma omp simd aligned(g_yz_yz_yz_xx, g_yz_yz_yz_xy, g_yz_yz_yz_xz, g_yz_yz_yz_yy, g_yz_yz_yz_yz, g_yz_yz_yz_zz, g_z_z_0_0_y_y_yz_xx, g_z_z_0_0_y_y_yz_xy, g_z_z_0_0_y_y_yz_xz, g_z_z_0_0_y_y_yz_yy, g_z_z_0_0_y_y_yz_yz, g_z_z_0_0_y_y_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_y_yz_xx[i] = 4.0 * g_yz_yz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_yz_xy[i] = 4.0 * g_yz_yz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_yz_xz[i] = 4.0 * g_yz_yz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_yz_yy[i] = 4.0 * g_yz_yz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_yz_yz[i] = 4.0 * g_yz_yz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_yz_zz[i] = 4.0 * g_yz_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2766-2772)

    #pragma omp simd aligned(g_yz_yz_zz_xx, g_yz_yz_zz_xy, g_yz_yz_zz_xz, g_yz_yz_zz_yy, g_yz_yz_zz_yz, g_yz_yz_zz_zz, g_z_z_0_0_y_y_zz_xx, g_z_z_0_0_y_y_zz_xy, g_z_z_0_0_y_y_zz_xz, g_z_z_0_0_y_y_zz_yy, g_z_z_0_0_y_y_zz_yz, g_z_z_0_0_y_y_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_y_zz_xx[i] = 4.0 * g_yz_yz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_zz_xy[i] = 4.0 * g_yz_yz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_zz_xz[i] = 4.0 * g_yz_yz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_zz_yy[i] = 4.0 * g_yz_yz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_zz_yz[i] = 4.0 * g_yz_yz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_zz_zz[i] = 4.0 * g_yz_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2772-2778)

    #pragma omp simd aligned(g_yz_0_xx_xx, g_yz_0_xx_xy, g_yz_0_xx_xz, g_yz_0_xx_yy, g_yz_0_xx_yz, g_yz_0_xx_zz, g_yz_zz_xx_xx, g_yz_zz_xx_xy, g_yz_zz_xx_xz, g_yz_zz_xx_yy, g_yz_zz_xx_yz, g_yz_zz_xx_zz, g_z_z_0_0_y_z_xx_xx, g_z_z_0_0_y_z_xx_xy, g_z_z_0_0_y_z_xx_xz, g_z_z_0_0_y_z_xx_yy, g_z_z_0_0_y_z_xx_yz, g_z_z_0_0_y_z_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_z_xx_xx[i] = -2.0 * g_yz_0_xx_xx[i] * a_exp + 4.0 * g_yz_zz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_xx_xy[i] = -2.0 * g_yz_0_xx_xy[i] * a_exp + 4.0 * g_yz_zz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_xx_xz[i] = -2.0 * g_yz_0_xx_xz[i] * a_exp + 4.0 * g_yz_zz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_xx_yy[i] = -2.0 * g_yz_0_xx_yy[i] * a_exp + 4.0 * g_yz_zz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_xx_yz[i] = -2.0 * g_yz_0_xx_yz[i] * a_exp + 4.0 * g_yz_zz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_xx_zz[i] = -2.0 * g_yz_0_xx_zz[i] * a_exp + 4.0 * g_yz_zz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2778-2784)

    #pragma omp simd aligned(g_yz_0_xy_xx, g_yz_0_xy_xy, g_yz_0_xy_xz, g_yz_0_xy_yy, g_yz_0_xy_yz, g_yz_0_xy_zz, g_yz_zz_xy_xx, g_yz_zz_xy_xy, g_yz_zz_xy_xz, g_yz_zz_xy_yy, g_yz_zz_xy_yz, g_yz_zz_xy_zz, g_z_z_0_0_y_z_xy_xx, g_z_z_0_0_y_z_xy_xy, g_z_z_0_0_y_z_xy_xz, g_z_z_0_0_y_z_xy_yy, g_z_z_0_0_y_z_xy_yz, g_z_z_0_0_y_z_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_z_xy_xx[i] = -2.0 * g_yz_0_xy_xx[i] * a_exp + 4.0 * g_yz_zz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_xy_xy[i] = -2.0 * g_yz_0_xy_xy[i] * a_exp + 4.0 * g_yz_zz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_xy_xz[i] = -2.0 * g_yz_0_xy_xz[i] * a_exp + 4.0 * g_yz_zz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_xy_yy[i] = -2.0 * g_yz_0_xy_yy[i] * a_exp + 4.0 * g_yz_zz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_xy_yz[i] = -2.0 * g_yz_0_xy_yz[i] * a_exp + 4.0 * g_yz_zz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_xy_zz[i] = -2.0 * g_yz_0_xy_zz[i] * a_exp + 4.0 * g_yz_zz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2784-2790)

    #pragma omp simd aligned(g_yz_0_xz_xx, g_yz_0_xz_xy, g_yz_0_xz_xz, g_yz_0_xz_yy, g_yz_0_xz_yz, g_yz_0_xz_zz, g_yz_zz_xz_xx, g_yz_zz_xz_xy, g_yz_zz_xz_xz, g_yz_zz_xz_yy, g_yz_zz_xz_yz, g_yz_zz_xz_zz, g_z_z_0_0_y_z_xz_xx, g_z_z_0_0_y_z_xz_xy, g_z_z_0_0_y_z_xz_xz, g_z_z_0_0_y_z_xz_yy, g_z_z_0_0_y_z_xz_yz, g_z_z_0_0_y_z_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_z_xz_xx[i] = -2.0 * g_yz_0_xz_xx[i] * a_exp + 4.0 * g_yz_zz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_xz_xy[i] = -2.0 * g_yz_0_xz_xy[i] * a_exp + 4.0 * g_yz_zz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_xz_xz[i] = -2.0 * g_yz_0_xz_xz[i] * a_exp + 4.0 * g_yz_zz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_xz_yy[i] = -2.0 * g_yz_0_xz_yy[i] * a_exp + 4.0 * g_yz_zz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_xz_yz[i] = -2.0 * g_yz_0_xz_yz[i] * a_exp + 4.0 * g_yz_zz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_xz_zz[i] = -2.0 * g_yz_0_xz_zz[i] * a_exp + 4.0 * g_yz_zz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2790-2796)

    #pragma omp simd aligned(g_yz_0_yy_xx, g_yz_0_yy_xy, g_yz_0_yy_xz, g_yz_0_yy_yy, g_yz_0_yy_yz, g_yz_0_yy_zz, g_yz_zz_yy_xx, g_yz_zz_yy_xy, g_yz_zz_yy_xz, g_yz_zz_yy_yy, g_yz_zz_yy_yz, g_yz_zz_yy_zz, g_z_z_0_0_y_z_yy_xx, g_z_z_0_0_y_z_yy_xy, g_z_z_0_0_y_z_yy_xz, g_z_z_0_0_y_z_yy_yy, g_z_z_0_0_y_z_yy_yz, g_z_z_0_0_y_z_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_z_yy_xx[i] = -2.0 * g_yz_0_yy_xx[i] * a_exp + 4.0 * g_yz_zz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_yy_xy[i] = -2.0 * g_yz_0_yy_xy[i] * a_exp + 4.0 * g_yz_zz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_yy_xz[i] = -2.0 * g_yz_0_yy_xz[i] * a_exp + 4.0 * g_yz_zz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_yy_yy[i] = -2.0 * g_yz_0_yy_yy[i] * a_exp + 4.0 * g_yz_zz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_yy_yz[i] = -2.0 * g_yz_0_yy_yz[i] * a_exp + 4.0 * g_yz_zz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_yy_zz[i] = -2.0 * g_yz_0_yy_zz[i] * a_exp + 4.0 * g_yz_zz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2796-2802)

    #pragma omp simd aligned(g_yz_0_yz_xx, g_yz_0_yz_xy, g_yz_0_yz_xz, g_yz_0_yz_yy, g_yz_0_yz_yz, g_yz_0_yz_zz, g_yz_zz_yz_xx, g_yz_zz_yz_xy, g_yz_zz_yz_xz, g_yz_zz_yz_yy, g_yz_zz_yz_yz, g_yz_zz_yz_zz, g_z_z_0_0_y_z_yz_xx, g_z_z_0_0_y_z_yz_xy, g_z_z_0_0_y_z_yz_xz, g_z_z_0_0_y_z_yz_yy, g_z_z_0_0_y_z_yz_yz, g_z_z_0_0_y_z_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_z_yz_xx[i] = -2.0 * g_yz_0_yz_xx[i] * a_exp + 4.0 * g_yz_zz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_yz_xy[i] = -2.0 * g_yz_0_yz_xy[i] * a_exp + 4.0 * g_yz_zz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_yz_xz[i] = -2.0 * g_yz_0_yz_xz[i] * a_exp + 4.0 * g_yz_zz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_yz_yy[i] = -2.0 * g_yz_0_yz_yy[i] * a_exp + 4.0 * g_yz_zz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_yz_yz[i] = -2.0 * g_yz_0_yz_yz[i] * a_exp + 4.0 * g_yz_zz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_yz_zz[i] = -2.0 * g_yz_0_yz_zz[i] * a_exp + 4.0 * g_yz_zz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2802-2808)

    #pragma omp simd aligned(g_yz_0_zz_xx, g_yz_0_zz_xy, g_yz_0_zz_xz, g_yz_0_zz_yy, g_yz_0_zz_yz, g_yz_0_zz_zz, g_yz_zz_zz_xx, g_yz_zz_zz_xy, g_yz_zz_zz_xz, g_yz_zz_zz_yy, g_yz_zz_zz_yz, g_yz_zz_zz_zz, g_z_z_0_0_y_z_zz_xx, g_z_z_0_0_y_z_zz_xy, g_z_z_0_0_y_z_zz_xz, g_z_z_0_0_y_z_zz_yy, g_z_z_0_0_y_z_zz_yz, g_z_z_0_0_y_z_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_z_zz_xx[i] = -2.0 * g_yz_0_zz_xx[i] * a_exp + 4.0 * g_yz_zz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_zz_xy[i] = -2.0 * g_yz_0_zz_xy[i] * a_exp + 4.0 * g_yz_zz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_zz_xz[i] = -2.0 * g_yz_0_zz_xz[i] * a_exp + 4.0 * g_yz_zz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_zz_yy[i] = -2.0 * g_yz_0_zz_yy[i] * a_exp + 4.0 * g_yz_zz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_zz_yz[i] = -2.0 * g_yz_0_zz_yz[i] * a_exp + 4.0 * g_yz_zz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_zz_zz[i] = -2.0 * g_yz_0_zz_zz[i] * a_exp + 4.0 * g_yz_zz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2808-2814)

    #pragma omp simd aligned(g_0_xz_xx_xx, g_0_xz_xx_xy, g_0_xz_xx_xz, g_0_xz_xx_yy, g_0_xz_xx_yz, g_0_xz_xx_zz, g_z_z_0_0_z_x_xx_xx, g_z_z_0_0_z_x_xx_xy, g_z_z_0_0_z_x_xx_xz, g_z_z_0_0_z_x_xx_yy, g_z_z_0_0_z_x_xx_yz, g_z_z_0_0_z_x_xx_zz, g_zz_xz_xx_xx, g_zz_xz_xx_xy, g_zz_xz_xx_xz, g_zz_xz_xx_yy, g_zz_xz_xx_yz, g_zz_xz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_x_xx_xx[i] = -2.0 * g_0_xz_xx_xx[i] * b_exp + 4.0 * g_zz_xz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_xx_xy[i] = -2.0 * g_0_xz_xx_xy[i] * b_exp + 4.0 * g_zz_xz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_xx_xz[i] = -2.0 * g_0_xz_xx_xz[i] * b_exp + 4.0 * g_zz_xz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_xx_yy[i] = -2.0 * g_0_xz_xx_yy[i] * b_exp + 4.0 * g_zz_xz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_xx_yz[i] = -2.0 * g_0_xz_xx_yz[i] * b_exp + 4.0 * g_zz_xz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_xx_zz[i] = -2.0 * g_0_xz_xx_zz[i] * b_exp + 4.0 * g_zz_xz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2814-2820)

    #pragma omp simd aligned(g_0_xz_xy_xx, g_0_xz_xy_xy, g_0_xz_xy_xz, g_0_xz_xy_yy, g_0_xz_xy_yz, g_0_xz_xy_zz, g_z_z_0_0_z_x_xy_xx, g_z_z_0_0_z_x_xy_xy, g_z_z_0_0_z_x_xy_xz, g_z_z_0_0_z_x_xy_yy, g_z_z_0_0_z_x_xy_yz, g_z_z_0_0_z_x_xy_zz, g_zz_xz_xy_xx, g_zz_xz_xy_xy, g_zz_xz_xy_xz, g_zz_xz_xy_yy, g_zz_xz_xy_yz, g_zz_xz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_x_xy_xx[i] = -2.0 * g_0_xz_xy_xx[i] * b_exp + 4.0 * g_zz_xz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_xy_xy[i] = -2.0 * g_0_xz_xy_xy[i] * b_exp + 4.0 * g_zz_xz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_xy_xz[i] = -2.0 * g_0_xz_xy_xz[i] * b_exp + 4.0 * g_zz_xz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_xy_yy[i] = -2.0 * g_0_xz_xy_yy[i] * b_exp + 4.0 * g_zz_xz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_xy_yz[i] = -2.0 * g_0_xz_xy_yz[i] * b_exp + 4.0 * g_zz_xz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_xy_zz[i] = -2.0 * g_0_xz_xy_zz[i] * b_exp + 4.0 * g_zz_xz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2820-2826)

    #pragma omp simd aligned(g_0_xz_xz_xx, g_0_xz_xz_xy, g_0_xz_xz_xz, g_0_xz_xz_yy, g_0_xz_xz_yz, g_0_xz_xz_zz, g_z_z_0_0_z_x_xz_xx, g_z_z_0_0_z_x_xz_xy, g_z_z_0_0_z_x_xz_xz, g_z_z_0_0_z_x_xz_yy, g_z_z_0_0_z_x_xz_yz, g_z_z_0_0_z_x_xz_zz, g_zz_xz_xz_xx, g_zz_xz_xz_xy, g_zz_xz_xz_xz, g_zz_xz_xz_yy, g_zz_xz_xz_yz, g_zz_xz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_x_xz_xx[i] = -2.0 * g_0_xz_xz_xx[i] * b_exp + 4.0 * g_zz_xz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_xz_xy[i] = -2.0 * g_0_xz_xz_xy[i] * b_exp + 4.0 * g_zz_xz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_xz_xz[i] = -2.0 * g_0_xz_xz_xz[i] * b_exp + 4.0 * g_zz_xz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_xz_yy[i] = -2.0 * g_0_xz_xz_yy[i] * b_exp + 4.0 * g_zz_xz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_xz_yz[i] = -2.0 * g_0_xz_xz_yz[i] * b_exp + 4.0 * g_zz_xz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_xz_zz[i] = -2.0 * g_0_xz_xz_zz[i] * b_exp + 4.0 * g_zz_xz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2826-2832)

    #pragma omp simd aligned(g_0_xz_yy_xx, g_0_xz_yy_xy, g_0_xz_yy_xz, g_0_xz_yy_yy, g_0_xz_yy_yz, g_0_xz_yy_zz, g_z_z_0_0_z_x_yy_xx, g_z_z_0_0_z_x_yy_xy, g_z_z_0_0_z_x_yy_xz, g_z_z_0_0_z_x_yy_yy, g_z_z_0_0_z_x_yy_yz, g_z_z_0_0_z_x_yy_zz, g_zz_xz_yy_xx, g_zz_xz_yy_xy, g_zz_xz_yy_xz, g_zz_xz_yy_yy, g_zz_xz_yy_yz, g_zz_xz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_x_yy_xx[i] = -2.0 * g_0_xz_yy_xx[i] * b_exp + 4.0 * g_zz_xz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_yy_xy[i] = -2.0 * g_0_xz_yy_xy[i] * b_exp + 4.0 * g_zz_xz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_yy_xz[i] = -2.0 * g_0_xz_yy_xz[i] * b_exp + 4.0 * g_zz_xz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_yy_yy[i] = -2.0 * g_0_xz_yy_yy[i] * b_exp + 4.0 * g_zz_xz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_yy_yz[i] = -2.0 * g_0_xz_yy_yz[i] * b_exp + 4.0 * g_zz_xz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_yy_zz[i] = -2.0 * g_0_xz_yy_zz[i] * b_exp + 4.0 * g_zz_xz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2832-2838)

    #pragma omp simd aligned(g_0_xz_yz_xx, g_0_xz_yz_xy, g_0_xz_yz_xz, g_0_xz_yz_yy, g_0_xz_yz_yz, g_0_xz_yz_zz, g_z_z_0_0_z_x_yz_xx, g_z_z_0_0_z_x_yz_xy, g_z_z_0_0_z_x_yz_xz, g_z_z_0_0_z_x_yz_yy, g_z_z_0_0_z_x_yz_yz, g_z_z_0_0_z_x_yz_zz, g_zz_xz_yz_xx, g_zz_xz_yz_xy, g_zz_xz_yz_xz, g_zz_xz_yz_yy, g_zz_xz_yz_yz, g_zz_xz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_x_yz_xx[i] = -2.0 * g_0_xz_yz_xx[i] * b_exp + 4.0 * g_zz_xz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_yz_xy[i] = -2.0 * g_0_xz_yz_xy[i] * b_exp + 4.0 * g_zz_xz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_yz_xz[i] = -2.0 * g_0_xz_yz_xz[i] * b_exp + 4.0 * g_zz_xz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_yz_yy[i] = -2.0 * g_0_xz_yz_yy[i] * b_exp + 4.0 * g_zz_xz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_yz_yz[i] = -2.0 * g_0_xz_yz_yz[i] * b_exp + 4.0 * g_zz_xz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_yz_zz[i] = -2.0 * g_0_xz_yz_zz[i] * b_exp + 4.0 * g_zz_xz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2838-2844)

    #pragma omp simd aligned(g_0_xz_zz_xx, g_0_xz_zz_xy, g_0_xz_zz_xz, g_0_xz_zz_yy, g_0_xz_zz_yz, g_0_xz_zz_zz, g_z_z_0_0_z_x_zz_xx, g_z_z_0_0_z_x_zz_xy, g_z_z_0_0_z_x_zz_xz, g_z_z_0_0_z_x_zz_yy, g_z_z_0_0_z_x_zz_yz, g_z_z_0_0_z_x_zz_zz, g_zz_xz_zz_xx, g_zz_xz_zz_xy, g_zz_xz_zz_xz, g_zz_xz_zz_yy, g_zz_xz_zz_yz, g_zz_xz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_x_zz_xx[i] = -2.0 * g_0_xz_zz_xx[i] * b_exp + 4.0 * g_zz_xz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_zz_xy[i] = -2.0 * g_0_xz_zz_xy[i] * b_exp + 4.0 * g_zz_xz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_zz_xz[i] = -2.0 * g_0_xz_zz_xz[i] * b_exp + 4.0 * g_zz_xz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_zz_yy[i] = -2.0 * g_0_xz_zz_yy[i] * b_exp + 4.0 * g_zz_xz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_zz_yz[i] = -2.0 * g_0_xz_zz_yz[i] * b_exp + 4.0 * g_zz_xz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_zz_zz[i] = -2.0 * g_0_xz_zz_zz[i] * b_exp + 4.0 * g_zz_xz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2844-2850)

    #pragma omp simd aligned(g_0_yz_xx_xx, g_0_yz_xx_xy, g_0_yz_xx_xz, g_0_yz_xx_yy, g_0_yz_xx_yz, g_0_yz_xx_zz, g_z_z_0_0_z_y_xx_xx, g_z_z_0_0_z_y_xx_xy, g_z_z_0_0_z_y_xx_xz, g_z_z_0_0_z_y_xx_yy, g_z_z_0_0_z_y_xx_yz, g_z_z_0_0_z_y_xx_zz, g_zz_yz_xx_xx, g_zz_yz_xx_xy, g_zz_yz_xx_xz, g_zz_yz_xx_yy, g_zz_yz_xx_yz, g_zz_yz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_y_xx_xx[i] = -2.0 * g_0_yz_xx_xx[i] * b_exp + 4.0 * g_zz_yz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_xx_xy[i] = -2.0 * g_0_yz_xx_xy[i] * b_exp + 4.0 * g_zz_yz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_xx_xz[i] = -2.0 * g_0_yz_xx_xz[i] * b_exp + 4.0 * g_zz_yz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_xx_yy[i] = -2.0 * g_0_yz_xx_yy[i] * b_exp + 4.0 * g_zz_yz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_xx_yz[i] = -2.0 * g_0_yz_xx_yz[i] * b_exp + 4.0 * g_zz_yz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_xx_zz[i] = -2.0 * g_0_yz_xx_zz[i] * b_exp + 4.0 * g_zz_yz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2850-2856)

    #pragma omp simd aligned(g_0_yz_xy_xx, g_0_yz_xy_xy, g_0_yz_xy_xz, g_0_yz_xy_yy, g_0_yz_xy_yz, g_0_yz_xy_zz, g_z_z_0_0_z_y_xy_xx, g_z_z_0_0_z_y_xy_xy, g_z_z_0_0_z_y_xy_xz, g_z_z_0_0_z_y_xy_yy, g_z_z_0_0_z_y_xy_yz, g_z_z_0_0_z_y_xy_zz, g_zz_yz_xy_xx, g_zz_yz_xy_xy, g_zz_yz_xy_xz, g_zz_yz_xy_yy, g_zz_yz_xy_yz, g_zz_yz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_y_xy_xx[i] = -2.0 * g_0_yz_xy_xx[i] * b_exp + 4.0 * g_zz_yz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_xy_xy[i] = -2.0 * g_0_yz_xy_xy[i] * b_exp + 4.0 * g_zz_yz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_xy_xz[i] = -2.0 * g_0_yz_xy_xz[i] * b_exp + 4.0 * g_zz_yz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_xy_yy[i] = -2.0 * g_0_yz_xy_yy[i] * b_exp + 4.0 * g_zz_yz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_xy_yz[i] = -2.0 * g_0_yz_xy_yz[i] * b_exp + 4.0 * g_zz_yz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_xy_zz[i] = -2.0 * g_0_yz_xy_zz[i] * b_exp + 4.0 * g_zz_yz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2856-2862)

    #pragma omp simd aligned(g_0_yz_xz_xx, g_0_yz_xz_xy, g_0_yz_xz_xz, g_0_yz_xz_yy, g_0_yz_xz_yz, g_0_yz_xz_zz, g_z_z_0_0_z_y_xz_xx, g_z_z_0_0_z_y_xz_xy, g_z_z_0_0_z_y_xz_xz, g_z_z_0_0_z_y_xz_yy, g_z_z_0_0_z_y_xz_yz, g_z_z_0_0_z_y_xz_zz, g_zz_yz_xz_xx, g_zz_yz_xz_xy, g_zz_yz_xz_xz, g_zz_yz_xz_yy, g_zz_yz_xz_yz, g_zz_yz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_y_xz_xx[i] = -2.0 * g_0_yz_xz_xx[i] * b_exp + 4.0 * g_zz_yz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_xz_xy[i] = -2.0 * g_0_yz_xz_xy[i] * b_exp + 4.0 * g_zz_yz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_xz_xz[i] = -2.0 * g_0_yz_xz_xz[i] * b_exp + 4.0 * g_zz_yz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_xz_yy[i] = -2.0 * g_0_yz_xz_yy[i] * b_exp + 4.0 * g_zz_yz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_xz_yz[i] = -2.0 * g_0_yz_xz_yz[i] * b_exp + 4.0 * g_zz_yz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_xz_zz[i] = -2.0 * g_0_yz_xz_zz[i] * b_exp + 4.0 * g_zz_yz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2862-2868)

    #pragma omp simd aligned(g_0_yz_yy_xx, g_0_yz_yy_xy, g_0_yz_yy_xz, g_0_yz_yy_yy, g_0_yz_yy_yz, g_0_yz_yy_zz, g_z_z_0_0_z_y_yy_xx, g_z_z_0_0_z_y_yy_xy, g_z_z_0_0_z_y_yy_xz, g_z_z_0_0_z_y_yy_yy, g_z_z_0_0_z_y_yy_yz, g_z_z_0_0_z_y_yy_zz, g_zz_yz_yy_xx, g_zz_yz_yy_xy, g_zz_yz_yy_xz, g_zz_yz_yy_yy, g_zz_yz_yy_yz, g_zz_yz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_y_yy_xx[i] = -2.0 * g_0_yz_yy_xx[i] * b_exp + 4.0 * g_zz_yz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_yy_xy[i] = -2.0 * g_0_yz_yy_xy[i] * b_exp + 4.0 * g_zz_yz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_yy_xz[i] = -2.0 * g_0_yz_yy_xz[i] * b_exp + 4.0 * g_zz_yz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_yy_yy[i] = -2.0 * g_0_yz_yy_yy[i] * b_exp + 4.0 * g_zz_yz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_yy_yz[i] = -2.0 * g_0_yz_yy_yz[i] * b_exp + 4.0 * g_zz_yz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_yy_zz[i] = -2.0 * g_0_yz_yy_zz[i] * b_exp + 4.0 * g_zz_yz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2868-2874)

    #pragma omp simd aligned(g_0_yz_yz_xx, g_0_yz_yz_xy, g_0_yz_yz_xz, g_0_yz_yz_yy, g_0_yz_yz_yz, g_0_yz_yz_zz, g_z_z_0_0_z_y_yz_xx, g_z_z_0_0_z_y_yz_xy, g_z_z_0_0_z_y_yz_xz, g_z_z_0_0_z_y_yz_yy, g_z_z_0_0_z_y_yz_yz, g_z_z_0_0_z_y_yz_zz, g_zz_yz_yz_xx, g_zz_yz_yz_xy, g_zz_yz_yz_xz, g_zz_yz_yz_yy, g_zz_yz_yz_yz, g_zz_yz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_y_yz_xx[i] = -2.0 * g_0_yz_yz_xx[i] * b_exp + 4.0 * g_zz_yz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_yz_xy[i] = -2.0 * g_0_yz_yz_xy[i] * b_exp + 4.0 * g_zz_yz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_yz_xz[i] = -2.0 * g_0_yz_yz_xz[i] * b_exp + 4.0 * g_zz_yz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_yz_yy[i] = -2.0 * g_0_yz_yz_yy[i] * b_exp + 4.0 * g_zz_yz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_yz_yz[i] = -2.0 * g_0_yz_yz_yz[i] * b_exp + 4.0 * g_zz_yz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_yz_zz[i] = -2.0 * g_0_yz_yz_zz[i] * b_exp + 4.0 * g_zz_yz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2874-2880)

    #pragma omp simd aligned(g_0_yz_zz_xx, g_0_yz_zz_xy, g_0_yz_zz_xz, g_0_yz_zz_yy, g_0_yz_zz_yz, g_0_yz_zz_zz, g_z_z_0_0_z_y_zz_xx, g_z_z_0_0_z_y_zz_xy, g_z_z_0_0_z_y_zz_xz, g_z_z_0_0_z_y_zz_yy, g_z_z_0_0_z_y_zz_yz, g_z_z_0_0_z_y_zz_zz, g_zz_yz_zz_xx, g_zz_yz_zz_xy, g_zz_yz_zz_xz, g_zz_yz_zz_yy, g_zz_yz_zz_yz, g_zz_yz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_y_zz_xx[i] = -2.0 * g_0_yz_zz_xx[i] * b_exp + 4.0 * g_zz_yz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_zz_xy[i] = -2.0 * g_0_yz_zz_xy[i] * b_exp + 4.0 * g_zz_yz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_zz_xz[i] = -2.0 * g_0_yz_zz_xz[i] * b_exp + 4.0 * g_zz_yz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_zz_yy[i] = -2.0 * g_0_yz_zz_yy[i] * b_exp + 4.0 * g_zz_yz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_zz_yz[i] = -2.0 * g_0_yz_zz_yz[i] * b_exp + 4.0 * g_zz_yz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_zz_zz[i] = -2.0 * g_0_yz_zz_zz[i] * b_exp + 4.0 * g_zz_yz_zz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2880-2886)

    #pragma omp simd aligned(g_0_0_xx_xx, g_0_0_xx_xy, g_0_0_xx_xz, g_0_0_xx_yy, g_0_0_xx_yz, g_0_0_xx_zz, g_0_zz_xx_xx, g_0_zz_xx_xy, g_0_zz_xx_xz, g_0_zz_xx_yy, g_0_zz_xx_yz, g_0_zz_xx_zz, g_z_z_0_0_z_z_xx_xx, g_z_z_0_0_z_z_xx_xy, g_z_z_0_0_z_z_xx_xz, g_z_z_0_0_z_z_xx_yy, g_z_z_0_0_z_z_xx_yz, g_z_z_0_0_z_z_xx_zz, g_zz_0_xx_xx, g_zz_0_xx_xy, g_zz_0_xx_xz, g_zz_0_xx_yy, g_zz_0_xx_yz, g_zz_0_xx_zz, g_zz_zz_xx_xx, g_zz_zz_xx_xy, g_zz_zz_xx_xz, g_zz_zz_xx_yy, g_zz_zz_xx_yz, g_zz_zz_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_z_xx_xx[i] = g_0_0_xx_xx[i] - 2.0 * g_0_zz_xx_xx[i] * b_exp - 2.0 * g_zz_0_xx_xx[i] * a_exp + 4.0 * g_zz_zz_xx_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_xx_xy[i] = g_0_0_xx_xy[i] - 2.0 * g_0_zz_xx_xy[i] * b_exp - 2.0 * g_zz_0_xx_xy[i] * a_exp + 4.0 * g_zz_zz_xx_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_xx_xz[i] = g_0_0_xx_xz[i] - 2.0 * g_0_zz_xx_xz[i] * b_exp - 2.0 * g_zz_0_xx_xz[i] * a_exp + 4.0 * g_zz_zz_xx_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_xx_yy[i] = g_0_0_xx_yy[i] - 2.0 * g_0_zz_xx_yy[i] * b_exp - 2.0 * g_zz_0_xx_yy[i] * a_exp + 4.0 * g_zz_zz_xx_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_xx_yz[i] = g_0_0_xx_yz[i] - 2.0 * g_0_zz_xx_yz[i] * b_exp - 2.0 * g_zz_0_xx_yz[i] * a_exp + 4.0 * g_zz_zz_xx_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_xx_zz[i] = g_0_0_xx_zz[i] - 2.0 * g_0_zz_xx_zz[i] * b_exp - 2.0 * g_zz_0_xx_zz[i] * a_exp + 4.0 * g_zz_zz_xx_zz[i] * a_exp * b_exp;
    }
    // integrals block (2886-2892)

    #pragma omp simd aligned(g_0_0_xy_xx, g_0_0_xy_xy, g_0_0_xy_xz, g_0_0_xy_yy, g_0_0_xy_yz, g_0_0_xy_zz, g_0_zz_xy_xx, g_0_zz_xy_xy, g_0_zz_xy_xz, g_0_zz_xy_yy, g_0_zz_xy_yz, g_0_zz_xy_zz, g_z_z_0_0_z_z_xy_xx, g_z_z_0_0_z_z_xy_xy, g_z_z_0_0_z_z_xy_xz, g_z_z_0_0_z_z_xy_yy, g_z_z_0_0_z_z_xy_yz, g_z_z_0_0_z_z_xy_zz, g_zz_0_xy_xx, g_zz_0_xy_xy, g_zz_0_xy_xz, g_zz_0_xy_yy, g_zz_0_xy_yz, g_zz_0_xy_zz, g_zz_zz_xy_xx, g_zz_zz_xy_xy, g_zz_zz_xy_xz, g_zz_zz_xy_yy, g_zz_zz_xy_yz, g_zz_zz_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_z_xy_xx[i] = g_0_0_xy_xx[i] - 2.0 * g_0_zz_xy_xx[i] * b_exp - 2.0 * g_zz_0_xy_xx[i] * a_exp + 4.0 * g_zz_zz_xy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_xy_xy[i] = g_0_0_xy_xy[i] - 2.0 * g_0_zz_xy_xy[i] * b_exp - 2.0 * g_zz_0_xy_xy[i] * a_exp + 4.0 * g_zz_zz_xy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_xy_xz[i] = g_0_0_xy_xz[i] - 2.0 * g_0_zz_xy_xz[i] * b_exp - 2.0 * g_zz_0_xy_xz[i] * a_exp + 4.0 * g_zz_zz_xy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_xy_yy[i] = g_0_0_xy_yy[i] - 2.0 * g_0_zz_xy_yy[i] * b_exp - 2.0 * g_zz_0_xy_yy[i] * a_exp + 4.0 * g_zz_zz_xy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_xy_yz[i] = g_0_0_xy_yz[i] - 2.0 * g_0_zz_xy_yz[i] * b_exp - 2.0 * g_zz_0_xy_yz[i] * a_exp + 4.0 * g_zz_zz_xy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_xy_zz[i] = g_0_0_xy_zz[i] - 2.0 * g_0_zz_xy_zz[i] * b_exp - 2.0 * g_zz_0_xy_zz[i] * a_exp + 4.0 * g_zz_zz_xy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2892-2898)

    #pragma omp simd aligned(g_0_0_xz_xx, g_0_0_xz_xy, g_0_0_xz_xz, g_0_0_xz_yy, g_0_0_xz_yz, g_0_0_xz_zz, g_0_zz_xz_xx, g_0_zz_xz_xy, g_0_zz_xz_xz, g_0_zz_xz_yy, g_0_zz_xz_yz, g_0_zz_xz_zz, g_z_z_0_0_z_z_xz_xx, g_z_z_0_0_z_z_xz_xy, g_z_z_0_0_z_z_xz_xz, g_z_z_0_0_z_z_xz_yy, g_z_z_0_0_z_z_xz_yz, g_z_z_0_0_z_z_xz_zz, g_zz_0_xz_xx, g_zz_0_xz_xy, g_zz_0_xz_xz, g_zz_0_xz_yy, g_zz_0_xz_yz, g_zz_0_xz_zz, g_zz_zz_xz_xx, g_zz_zz_xz_xy, g_zz_zz_xz_xz, g_zz_zz_xz_yy, g_zz_zz_xz_yz, g_zz_zz_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_z_xz_xx[i] = g_0_0_xz_xx[i] - 2.0 * g_0_zz_xz_xx[i] * b_exp - 2.0 * g_zz_0_xz_xx[i] * a_exp + 4.0 * g_zz_zz_xz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_xz_xy[i] = g_0_0_xz_xy[i] - 2.0 * g_0_zz_xz_xy[i] * b_exp - 2.0 * g_zz_0_xz_xy[i] * a_exp + 4.0 * g_zz_zz_xz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_xz_xz[i] = g_0_0_xz_xz[i] - 2.0 * g_0_zz_xz_xz[i] * b_exp - 2.0 * g_zz_0_xz_xz[i] * a_exp + 4.0 * g_zz_zz_xz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_xz_yy[i] = g_0_0_xz_yy[i] - 2.0 * g_0_zz_xz_yy[i] * b_exp - 2.0 * g_zz_0_xz_yy[i] * a_exp + 4.0 * g_zz_zz_xz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_xz_yz[i] = g_0_0_xz_yz[i] - 2.0 * g_0_zz_xz_yz[i] * b_exp - 2.0 * g_zz_0_xz_yz[i] * a_exp + 4.0 * g_zz_zz_xz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_xz_zz[i] = g_0_0_xz_zz[i] - 2.0 * g_0_zz_xz_zz[i] * b_exp - 2.0 * g_zz_0_xz_zz[i] * a_exp + 4.0 * g_zz_zz_xz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2898-2904)

    #pragma omp simd aligned(g_0_0_yy_xx, g_0_0_yy_xy, g_0_0_yy_xz, g_0_0_yy_yy, g_0_0_yy_yz, g_0_0_yy_zz, g_0_zz_yy_xx, g_0_zz_yy_xy, g_0_zz_yy_xz, g_0_zz_yy_yy, g_0_zz_yy_yz, g_0_zz_yy_zz, g_z_z_0_0_z_z_yy_xx, g_z_z_0_0_z_z_yy_xy, g_z_z_0_0_z_z_yy_xz, g_z_z_0_0_z_z_yy_yy, g_z_z_0_0_z_z_yy_yz, g_z_z_0_0_z_z_yy_zz, g_zz_0_yy_xx, g_zz_0_yy_xy, g_zz_0_yy_xz, g_zz_0_yy_yy, g_zz_0_yy_yz, g_zz_0_yy_zz, g_zz_zz_yy_xx, g_zz_zz_yy_xy, g_zz_zz_yy_xz, g_zz_zz_yy_yy, g_zz_zz_yy_yz, g_zz_zz_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_z_yy_xx[i] = g_0_0_yy_xx[i] - 2.0 * g_0_zz_yy_xx[i] * b_exp - 2.0 * g_zz_0_yy_xx[i] * a_exp + 4.0 * g_zz_zz_yy_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_yy_xy[i] = g_0_0_yy_xy[i] - 2.0 * g_0_zz_yy_xy[i] * b_exp - 2.0 * g_zz_0_yy_xy[i] * a_exp + 4.0 * g_zz_zz_yy_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_yy_xz[i] = g_0_0_yy_xz[i] - 2.0 * g_0_zz_yy_xz[i] * b_exp - 2.0 * g_zz_0_yy_xz[i] * a_exp + 4.0 * g_zz_zz_yy_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_yy_yy[i] = g_0_0_yy_yy[i] - 2.0 * g_0_zz_yy_yy[i] * b_exp - 2.0 * g_zz_0_yy_yy[i] * a_exp + 4.0 * g_zz_zz_yy_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_yy_yz[i] = g_0_0_yy_yz[i] - 2.0 * g_0_zz_yy_yz[i] * b_exp - 2.0 * g_zz_0_yy_yz[i] * a_exp + 4.0 * g_zz_zz_yy_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_yy_zz[i] = g_0_0_yy_zz[i] - 2.0 * g_0_zz_yy_zz[i] * b_exp - 2.0 * g_zz_0_yy_zz[i] * a_exp + 4.0 * g_zz_zz_yy_zz[i] * a_exp * b_exp;
    }
    // integrals block (2904-2910)

    #pragma omp simd aligned(g_0_0_yz_xx, g_0_0_yz_xy, g_0_0_yz_xz, g_0_0_yz_yy, g_0_0_yz_yz, g_0_0_yz_zz, g_0_zz_yz_xx, g_0_zz_yz_xy, g_0_zz_yz_xz, g_0_zz_yz_yy, g_0_zz_yz_yz, g_0_zz_yz_zz, g_z_z_0_0_z_z_yz_xx, g_z_z_0_0_z_z_yz_xy, g_z_z_0_0_z_z_yz_xz, g_z_z_0_0_z_z_yz_yy, g_z_z_0_0_z_z_yz_yz, g_z_z_0_0_z_z_yz_zz, g_zz_0_yz_xx, g_zz_0_yz_xy, g_zz_0_yz_xz, g_zz_0_yz_yy, g_zz_0_yz_yz, g_zz_0_yz_zz, g_zz_zz_yz_xx, g_zz_zz_yz_xy, g_zz_zz_yz_xz, g_zz_zz_yz_yy, g_zz_zz_yz_yz, g_zz_zz_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_z_yz_xx[i] = g_0_0_yz_xx[i] - 2.0 * g_0_zz_yz_xx[i] * b_exp - 2.0 * g_zz_0_yz_xx[i] * a_exp + 4.0 * g_zz_zz_yz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_yz_xy[i] = g_0_0_yz_xy[i] - 2.0 * g_0_zz_yz_xy[i] * b_exp - 2.0 * g_zz_0_yz_xy[i] * a_exp + 4.0 * g_zz_zz_yz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_yz_xz[i] = g_0_0_yz_xz[i] - 2.0 * g_0_zz_yz_xz[i] * b_exp - 2.0 * g_zz_0_yz_xz[i] * a_exp + 4.0 * g_zz_zz_yz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_yz_yy[i] = g_0_0_yz_yy[i] - 2.0 * g_0_zz_yz_yy[i] * b_exp - 2.0 * g_zz_0_yz_yy[i] * a_exp + 4.0 * g_zz_zz_yz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_yz_yz[i] = g_0_0_yz_yz[i] - 2.0 * g_0_zz_yz_yz[i] * b_exp - 2.0 * g_zz_0_yz_yz[i] * a_exp + 4.0 * g_zz_zz_yz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_yz_zz[i] = g_0_0_yz_zz[i] - 2.0 * g_0_zz_yz_zz[i] * b_exp - 2.0 * g_zz_0_yz_zz[i] * a_exp + 4.0 * g_zz_zz_yz_zz[i] * a_exp * b_exp;
    }
    // integrals block (2910-2916)

    #pragma omp simd aligned(g_0_0_zz_xx, g_0_0_zz_xy, g_0_0_zz_xz, g_0_0_zz_yy, g_0_0_zz_yz, g_0_0_zz_zz, g_0_zz_zz_xx, g_0_zz_zz_xy, g_0_zz_zz_xz, g_0_zz_zz_yy, g_0_zz_zz_yz, g_0_zz_zz_zz, g_z_z_0_0_z_z_zz_xx, g_z_z_0_0_z_z_zz_xy, g_z_z_0_0_z_z_zz_xz, g_z_z_0_0_z_z_zz_yy, g_z_z_0_0_z_z_zz_yz, g_z_z_0_0_z_z_zz_zz, g_zz_0_zz_xx, g_zz_0_zz_xy, g_zz_0_zz_xz, g_zz_0_zz_yy, g_zz_0_zz_yz, g_zz_0_zz_zz, g_zz_zz_zz_xx, g_zz_zz_zz_xy, g_zz_zz_zz_xz, g_zz_zz_zz_yy, g_zz_zz_zz_yz, g_zz_zz_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_z_zz_xx[i] = g_0_0_zz_xx[i] - 2.0 * g_0_zz_zz_xx[i] * b_exp - 2.0 * g_zz_0_zz_xx[i] * a_exp + 4.0 * g_zz_zz_zz_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_zz_xy[i] = g_0_0_zz_xy[i] - 2.0 * g_0_zz_zz_xy[i] * b_exp - 2.0 * g_zz_0_zz_xy[i] * a_exp + 4.0 * g_zz_zz_zz_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_zz_xz[i] = g_0_0_zz_xz[i] - 2.0 * g_0_zz_zz_xz[i] * b_exp - 2.0 * g_zz_0_zz_xz[i] * a_exp + 4.0 * g_zz_zz_zz_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_zz_yy[i] = g_0_0_zz_yy[i] - 2.0 * g_0_zz_zz_yy[i] * b_exp - 2.0 * g_zz_0_zz_yy[i] * a_exp + 4.0 * g_zz_zz_zz_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_zz_yz[i] = g_0_0_zz_yz[i] - 2.0 * g_0_zz_zz_yz[i] * b_exp - 2.0 * g_zz_0_zz_yz[i] * a_exp + 4.0 * g_zz_zz_zz_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_zz_zz[i] = g_0_0_zz_zz[i] - 2.0 * g_0_zz_zz_zz[i] * b_exp - 2.0 * g_zz_0_zz_zz[i] * a_exp + 4.0 * g_zz_zz_zz_zz[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

