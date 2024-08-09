#include "GeomDeriv2000OfScalarForSSDD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_ssdd_0(CSimdArray<double>& buffer_2000_ssdd,
                     const CSimdArray<double>& buffer_ssdd,
                     const CSimdArray<double>& buffer_dsdd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_ssdd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_2000_ssdd

    auto g_xx_0_0_0_0_0_xx_xx = buffer_2000_ssdd[0];

    auto g_xx_0_0_0_0_0_xx_xy = buffer_2000_ssdd[1];

    auto g_xx_0_0_0_0_0_xx_xz = buffer_2000_ssdd[2];

    auto g_xx_0_0_0_0_0_xx_yy = buffer_2000_ssdd[3];

    auto g_xx_0_0_0_0_0_xx_yz = buffer_2000_ssdd[4];

    auto g_xx_0_0_0_0_0_xx_zz = buffer_2000_ssdd[5];

    auto g_xx_0_0_0_0_0_xy_xx = buffer_2000_ssdd[6];

    auto g_xx_0_0_0_0_0_xy_xy = buffer_2000_ssdd[7];

    auto g_xx_0_0_0_0_0_xy_xz = buffer_2000_ssdd[8];

    auto g_xx_0_0_0_0_0_xy_yy = buffer_2000_ssdd[9];

    auto g_xx_0_0_0_0_0_xy_yz = buffer_2000_ssdd[10];

    auto g_xx_0_0_0_0_0_xy_zz = buffer_2000_ssdd[11];

    auto g_xx_0_0_0_0_0_xz_xx = buffer_2000_ssdd[12];

    auto g_xx_0_0_0_0_0_xz_xy = buffer_2000_ssdd[13];

    auto g_xx_0_0_0_0_0_xz_xz = buffer_2000_ssdd[14];

    auto g_xx_0_0_0_0_0_xz_yy = buffer_2000_ssdd[15];

    auto g_xx_0_0_0_0_0_xz_yz = buffer_2000_ssdd[16];

    auto g_xx_0_0_0_0_0_xz_zz = buffer_2000_ssdd[17];

    auto g_xx_0_0_0_0_0_yy_xx = buffer_2000_ssdd[18];

    auto g_xx_0_0_0_0_0_yy_xy = buffer_2000_ssdd[19];

    auto g_xx_0_0_0_0_0_yy_xz = buffer_2000_ssdd[20];

    auto g_xx_0_0_0_0_0_yy_yy = buffer_2000_ssdd[21];

    auto g_xx_0_0_0_0_0_yy_yz = buffer_2000_ssdd[22];

    auto g_xx_0_0_0_0_0_yy_zz = buffer_2000_ssdd[23];

    auto g_xx_0_0_0_0_0_yz_xx = buffer_2000_ssdd[24];

    auto g_xx_0_0_0_0_0_yz_xy = buffer_2000_ssdd[25];

    auto g_xx_0_0_0_0_0_yz_xz = buffer_2000_ssdd[26];

    auto g_xx_0_0_0_0_0_yz_yy = buffer_2000_ssdd[27];

    auto g_xx_0_0_0_0_0_yz_yz = buffer_2000_ssdd[28];

    auto g_xx_0_0_0_0_0_yz_zz = buffer_2000_ssdd[29];

    auto g_xx_0_0_0_0_0_zz_xx = buffer_2000_ssdd[30];

    auto g_xx_0_0_0_0_0_zz_xy = buffer_2000_ssdd[31];

    auto g_xx_0_0_0_0_0_zz_xz = buffer_2000_ssdd[32];

    auto g_xx_0_0_0_0_0_zz_yy = buffer_2000_ssdd[33];

    auto g_xx_0_0_0_0_0_zz_yz = buffer_2000_ssdd[34];

    auto g_xx_0_0_0_0_0_zz_zz = buffer_2000_ssdd[35];

    auto g_xy_0_0_0_0_0_xx_xx = buffer_2000_ssdd[36];

    auto g_xy_0_0_0_0_0_xx_xy = buffer_2000_ssdd[37];

    auto g_xy_0_0_0_0_0_xx_xz = buffer_2000_ssdd[38];

    auto g_xy_0_0_0_0_0_xx_yy = buffer_2000_ssdd[39];

    auto g_xy_0_0_0_0_0_xx_yz = buffer_2000_ssdd[40];

    auto g_xy_0_0_0_0_0_xx_zz = buffer_2000_ssdd[41];

    auto g_xy_0_0_0_0_0_xy_xx = buffer_2000_ssdd[42];

    auto g_xy_0_0_0_0_0_xy_xy = buffer_2000_ssdd[43];

    auto g_xy_0_0_0_0_0_xy_xz = buffer_2000_ssdd[44];

    auto g_xy_0_0_0_0_0_xy_yy = buffer_2000_ssdd[45];

    auto g_xy_0_0_0_0_0_xy_yz = buffer_2000_ssdd[46];

    auto g_xy_0_0_0_0_0_xy_zz = buffer_2000_ssdd[47];

    auto g_xy_0_0_0_0_0_xz_xx = buffer_2000_ssdd[48];

    auto g_xy_0_0_0_0_0_xz_xy = buffer_2000_ssdd[49];

    auto g_xy_0_0_0_0_0_xz_xz = buffer_2000_ssdd[50];

    auto g_xy_0_0_0_0_0_xz_yy = buffer_2000_ssdd[51];

    auto g_xy_0_0_0_0_0_xz_yz = buffer_2000_ssdd[52];

    auto g_xy_0_0_0_0_0_xz_zz = buffer_2000_ssdd[53];

    auto g_xy_0_0_0_0_0_yy_xx = buffer_2000_ssdd[54];

    auto g_xy_0_0_0_0_0_yy_xy = buffer_2000_ssdd[55];

    auto g_xy_0_0_0_0_0_yy_xz = buffer_2000_ssdd[56];

    auto g_xy_0_0_0_0_0_yy_yy = buffer_2000_ssdd[57];

    auto g_xy_0_0_0_0_0_yy_yz = buffer_2000_ssdd[58];

    auto g_xy_0_0_0_0_0_yy_zz = buffer_2000_ssdd[59];

    auto g_xy_0_0_0_0_0_yz_xx = buffer_2000_ssdd[60];

    auto g_xy_0_0_0_0_0_yz_xy = buffer_2000_ssdd[61];

    auto g_xy_0_0_0_0_0_yz_xz = buffer_2000_ssdd[62];

    auto g_xy_0_0_0_0_0_yz_yy = buffer_2000_ssdd[63];

    auto g_xy_0_0_0_0_0_yz_yz = buffer_2000_ssdd[64];

    auto g_xy_0_0_0_0_0_yz_zz = buffer_2000_ssdd[65];

    auto g_xy_0_0_0_0_0_zz_xx = buffer_2000_ssdd[66];

    auto g_xy_0_0_0_0_0_zz_xy = buffer_2000_ssdd[67];

    auto g_xy_0_0_0_0_0_zz_xz = buffer_2000_ssdd[68];

    auto g_xy_0_0_0_0_0_zz_yy = buffer_2000_ssdd[69];

    auto g_xy_0_0_0_0_0_zz_yz = buffer_2000_ssdd[70];

    auto g_xy_0_0_0_0_0_zz_zz = buffer_2000_ssdd[71];

    auto g_xz_0_0_0_0_0_xx_xx = buffer_2000_ssdd[72];

    auto g_xz_0_0_0_0_0_xx_xy = buffer_2000_ssdd[73];

    auto g_xz_0_0_0_0_0_xx_xz = buffer_2000_ssdd[74];

    auto g_xz_0_0_0_0_0_xx_yy = buffer_2000_ssdd[75];

    auto g_xz_0_0_0_0_0_xx_yz = buffer_2000_ssdd[76];

    auto g_xz_0_0_0_0_0_xx_zz = buffer_2000_ssdd[77];

    auto g_xz_0_0_0_0_0_xy_xx = buffer_2000_ssdd[78];

    auto g_xz_0_0_0_0_0_xy_xy = buffer_2000_ssdd[79];

    auto g_xz_0_0_0_0_0_xy_xz = buffer_2000_ssdd[80];

    auto g_xz_0_0_0_0_0_xy_yy = buffer_2000_ssdd[81];

    auto g_xz_0_0_0_0_0_xy_yz = buffer_2000_ssdd[82];

    auto g_xz_0_0_0_0_0_xy_zz = buffer_2000_ssdd[83];

    auto g_xz_0_0_0_0_0_xz_xx = buffer_2000_ssdd[84];

    auto g_xz_0_0_0_0_0_xz_xy = buffer_2000_ssdd[85];

    auto g_xz_0_0_0_0_0_xz_xz = buffer_2000_ssdd[86];

    auto g_xz_0_0_0_0_0_xz_yy = buffer_2000_ssdd[87];

    auto g_xz_0_0_0_0_0_xz_yz = buffer_2000_ssdd[88];

    auto g_xz_0_0_0_0_0_xz_zz = buffer_2000_ssdd[89];

    auto g_xz_0_0_0_0_0_yy_xx = buffer_2000_ssdd[90];

    auto g_xz_0_0_0_0_0_yy_xy = buffer_2000_ssdd[91];

    auto g_xz_0_0_0_0_0_yy_xz = buffer_2000_ssdd[92];

    auto g_xz_0_0_0_0_0_yy_yy = buffer_2000_ssdd[93];

    auto g_xz_0_0_0_0_0_yy_yz = buffer_2000_ssdd[94];

    auto g_xz_0_0_0_0_0_yy_zz = buffer_2000_ssdd[95];

    auto g_xz_0_0_0_0_0_yz_xx = buffer_2000_ssdd[96];

    auto g_xz_0_0_0_0_0_yz_xy = buffer_2000_ssdd[97];

    auto g_xz_0_0_0_0_0_yz_xz = buffer_2000_ssdd[98];

    auto g_xz_0_0_0_0_0_yz_yy = buffer_2000_ssdd[99];

    auto g_xz_0_0_0_0_0_yz_yz = buffer_2000_ssdd[100];

    auto g_xz_0_0_0_0_0_yz_zz = buffer_2000_ssdd[101];

    auto g_xz_0_0_0_0_0_zz_xx = buffer_2000_ssdd[102];

    auto g_xz_0_0_0_0_0_zz_xy = buffer_2000_ssdd[103];

    auto g_xz_0_0_0_0_0_zz_xz = buffer_2000_ssdd[104];

    auto g_xz_0_0_0_0_0_zz_yy = buffer_2000_ssdd[105];

    auto g_xz_0_0_0_0_0_zz_yz = buffer_2000_ssdd[106];

    auto g_xz_0_0_0_0_0_zz_zz = buffer_2000_ssdd[107];

    auto g_yy_0_0_0_0_0_xx_xx = buffer_2000_ssdd[108];

    auto g_yy_0_0_0_0_0_xx_xy = buffer_2000_ssdd[109];

    auto g_yy_0_0_0_0_0_xx_xz = buffer_2000_ssdd[110];

    auto g_yy_0_0_0_0_0_xx_yy = buffer_2000_ssdd[111];

    auto g_yy_0_0_0_0_0_xx_yz = buffer_2000_ssdd[112];

    auto g_yy_0_0_0_0_0_xx_zz = buffer_2000_ssdd[113];

    auto g_yy_0_0_0_0_0_xy_xx = buffer_2000_ssdd[114];

    auto g_yy_0_0_0_0_0_xy_xy = buffer_2000_ssdd[115];

    auto g_yy_0_0_0_0_0_xy_xz = buffer_2000_ssdd[116];

    auto g_yy_0_0_0_0_0_xy_yy = buffer_2000_ssdd[117];

    auto g_yy_0_0_0_0_0_xy_yz = buffer_2000_ssdd[118];

    auto g_yy_0_0_0_0_0_xy_zz = buffer_2000_ssdd[119];

    auto g_yy_0_0_0_0_0_xz_xx = buffer_2000_ssdd[120];

    auto g_yy_0_0_0_0_0_xz_xy = buffer_2000_ssdd[121];

    auto g_yy_0_0_0_0_0_xz_xz = buffer_2000_ssdd[122];

    auto g_yy_0_0_0_0_0_xz_yy = buffer_2000_ssdd[123];

    auto g_yy_0_0_0_0_0_xz_yz = buffer_2000_ssdd[124];

    auto g_yy_0_0_0_0_0_xz_zz = buffer_2000_ssdd[125];

    auto g_yy_0_0_0_0_0_yy_xx = buffer_2000_ssdd[126];

    auto g_yy_0_0_0_0_0_yy_xy = buffer_2000_ssdd[127];

    auto g_yy_0_0_0_0_0_yy_xz = buffer_2000_ssdd[128];

    auto g_yy_0_0_0_0_0_yy_yy = buffer_2000_ssdd[129];

    auto g_yy_0_0_0_0_0_yy_yz = buffer_2000_ssdd[130];

    auto g_yy_0_0_0_0_0_yy_zz = buffer_2000_ssdd[131];

    auto g_yy_0_0_0_0_0_yz_xx = buffer_2000_ssdd[132];

    auto g_yy_0_0_0_0_0_yz_xy = buffer_2000_ssdd[133];

    auto g_yy_0_0_0_0_0_yz_xz = buffer_2000_ssdd[134];

    auto g_yy_0_0_0_0_0_yz_yy = buffer_2000_ssdd[135];

    auto g_yy_0_0_0_0_0_yz_yz = buffer_2000_ssdd[136];

    auto g_yy_0_0_0_0_0_yz_zz = buffer_2000_ssdd[137];

    auto g_yy_0_0_0_0_0_zz_xx = buffer_2000_ssdd[138];

    auto g_yy_0_0_0_0_0_zz_xy = buffer_2000_ssdd[139];

    auto g_yy_0_0_0_0_0_zz_xz = buffer_2000_ssdd[140];

    auto g_yy_0_0_0_0_0_zz_yy = buffer_2000_ssdd[141];

    auto g_yy_0_0_0_0_0_zz_yz = buffer_2000_ssdd[142];

    auto g_yy_0_0_0_0_0_zz_zz = buffer_2000_ssdd[143];

    auto g_yz_0_0_0_0_0_xx_xx = buffer_2000_ssdd[144];

    auto g_yz_0_0_0_0_0_xx_xy = buffer_2000_ssdd[145];

    auto g_yz_0_0_0_0_0_xx_xz = buffer_2000_ssdd[146];

    auto g_yz_0_0_0_0_0_xx_yy = buffer_2000_ssdd[147];

    auto g_yz_0_0_0_0_0_xx_yz = buffer_2000_ssdd[148];

    auto g_yz_0_0_0_0_0_xx_zz = buffer_2000_ssdd[149];

    auto g_yz_0_0_0_0_0_xy_xx = buffer_2000_ssdd[150];

    auto g_yz_0_0_0_0_0_xy_xy = buffer_2000_ssdd[151];

    auto g_yz_0_0_0_0_0_xy_xz = buffer_2000_ssdd[152];

    auto g_yz_0_0_0_0_0_xy_yy = buffer_2000_ssdd[153];

    auto g_yz_0_0_0_0_0_xy_yz = buffer_2000_ssdd[154];

    auto g_yz_0_0_0_0_0_xy_zz = buffer_2000_ssdd[155];

    auto g_yz_0_0_0_0_0_xz_xx = buffer_2000_ssdd[156];

    auto g_yz_0_0_0_0_0_xz_xy = buffer_2000_ssdd[157];

    auto g_yz_0_0_0_0_0_xz_xz = buffer_2000_ssdd[158];

    auto g_yz_0_0_0_0_0_xz_yy = buffer_2000_ssdd[159];

    auto g_yz_0_0_0_0_0_xz_yz = buffer_2000_ssdd[160];

    auto g_yz_0_0_0_0_0_xz_zz = buffer_2000_ssdd[161];

    auto g_yz_0_0_0_0_0_yy_xx = buffer_2000_ssdd[162];

    auto g_yz_0_0_0_0_0_yy_xy = buffer_2000_ssdd[163];

    auto g_yz_0_0_0_0_0_yy_xz = buffer_2000_ssdd[164];

    auto g_yz_0_0_0_0_0_yy_yy = buffer_2000_ssdd[165];

    auto g_yz_0_0_0_0_0_yy_yz = buffer_2000_ssdd[166];

    auto g_yz_0_0_0_0_0_yy_zz = buffer_2000_ssdd[167];

    auto g_yz_0_0_0_0_0_yz_xx = buffer_2000_ssdd[168];

    auto g_yz_0_0_0_0_0_yz_xy = buffer_2000_ssdd[169];

    auto g_yz_0_0_0_0_0_yz_xz = buffer_2000_ssdd[170];

    auto g_yz_0_0_0_0_0_yz_yy = buffer_2000_ssdd[171];

    auto g_yz_0_0_0_0_0_yz_yz = buffer_2000_ssdd[172];

    auto g_yz_0_0_0_0_0_yz_zz = buffer_2000_ssdd[173];

    auto g_yz_0_0_0_0_0_zz_xx = buffer_2000_ssdd[174];

    auto g_yz_0_0_0_0_0_zz_xy = buffer_2000_ssdd[175];

    auto g_yz_0_0_0_0_0_zz_xz = buffer_2000_ssdd[176];

    auto g_yz_0_0_0_0_0_zz_yy = buffer_2000_ssdd[177];

    auto g_yz_0_0_0_0_0_zz_yz = buffer_2000_ssdd[178];

    auto g_yz_0_0_0_0_0_zz_zz = buffer_2000_ssdd[179];

    auto g_zz_0_0_0_0_0_xx_xx = buffer_2000_ssdd[180];

    auto g_zz_0_0_0_0_0_xx_xy = buffer_2000_ssdd[181];

    auto g_zz_0_0_0_0_0_xx_xz = buffer_2000_ssdd[182];

    auto g_zz_0_0_0_0_0_xx_yy = buffer_2000_ssdd[183];

    auto g_zz_0_0_0_0_0_xx_yz = buffer_2000_ssdd[184];

    auto g_zz_0_0_0_0_0_xx_zz = buffer_2000_ssdd[185];

    auto g_zz_0_0_0_0_0_xy_xx = buffer_2000_ssdd[186];

    auto g_zz_0_0_0_0_0_xy_xy = buffer_2000_ssdd[187];

    auto g_zz_0_0_0_0_0_xy_xz = buffer_2000_ssdd[188];

    auto g_zz_0_0_0_0_0_xy_yy = buffer_2000_ssdd[189];

    auto g_zz_0_0_0_0_0_xy_yz = buffer_2000_ssdd[190];

    auto g_zz_0_0_0_0_0_xy_zz = buffer_2000_ssdd[191];

    auto g_zz_0_0_0_0_0_xz_xx = buffer_2000_ssdd[192];

    auto g_zz_0_0_0_0_0_xz_xy = buffer_2000_ssdd[193];

    auto g_zz_0_0_0_0_0_xz_xz = buffer_2000_ssdd[194];

    auto g_zz_0_0_0_0_0_xz_yy = buffer_2000_ssdd[195];

    auto g_zz_0_0_0_0_0_xz_yz = buffer_2000_ssdd[196];

    auto g_zz_0_0_0_0_0_xz_zz = buffer_2000_ssdd[197];

    auto g_zz_0_0_0_0_0_yy_xx = buffer_2000_ssdd[198];

    auto g_zz_0_0_0_0_0_yy_xy = buffer_2000_ssdd[199];

    auto g_zz_0_0_0_0_0_yy_xz = buffer_2000_ssdd[200];

    auto g_zz_0_0_0_0_0_yy_yy = buffer_2000_ssdd[201];

    auto g_zz_0_0_0_0_0_yy_yz = buffer_2000_ssdd[202];

    auto g_zz_0_0_0_0_0_yy_zz = buffer_2000_ssdd[203];

    auto g_zz_0_0_0_0_0_yz_xx = buffer_2000_ssdd[204];

    auto g_zz_0_0_0_0_0_yz_xy = buffer_2000_ssdd[205];

    auto g_zz_0_0_0_0_0_yz_xz = buffer_2000_ssdd[206];

    auto g_zz_0_0_0_0_0_yz_yy = buffer_2000_ssdd[207];

    auto g_zz_0_0_0_0_0_yz_yz = buffer_2000_ssdd[208];

    auto g_zz_0_0_0_0_0_yz_zz = buffer_2000_ssdd[209];

    auto g_zz_0_0_0_0_0_zz_xx = buffer_2000_ssdd[210];

    auto g_zz_0_0_0_0_0_zz_xy = buffer_2000_ssdd[211];

    auto g_zz_0_0_0_0_0_zz_xz = buffer_2000_ssdd[212];

    auto g_zz_0_0_0_0_0_zz_yy = buffer_2000_ssdd[213];

    auto g_zz_0_0_0_0_0_zz_yz = buffer_2000_ssdd[214];

    auto g_zz_0_0_0_0_0_zz_zz = buffer_2000_ssdd[215];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_0_xx_xx, g_0_0_xx_xy, g_0_0_xx_xz, g_0_0_xx_yy, g_0_0_xx_yz, g_0_0_xx_zz, g_xx_0_0_0_0_0_xx_xx, g_xx_0_0_0_0_0_xx_xy, g_xx_0_0_0_0_0_xx_xz, g_xx_0_0_0_0_0_xx_yy, g_xx_0_0_0_0_0_xx_yz, g_xx_0_0_0_0_0_xx_zz, g_xx_0_xx_xx, g_xx_0_xx_xy, g_xx_0_xx_xz, g_xx_0_xx_yy, g_xx_0_xx_yz, g_xx_0_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_0_xx_xx[i] = -2.0 * g_0_0_xx_xx[i] * a_exp + 4.0 * g_xx_0_xx_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_xx_xy[i] = -2.0 * g_0_0_xx_xy[i] * a_exp + 4.0 * g_xx_0_xx_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_xx_xz[i] = -2.0 * g_0_0_xx_xz[i] * a_exp + 4.0 * g_xx_0_xx_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_xx_yy[i] = -2.0 * g_0_0_xx_yy[i] * a_exp + 4.0 * g_xx_0_xx_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_xx_yz[i] = -2.0 * g_0_0_xx_yz[i] * a_exp + 4.0 * g_xx_0_xx_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_xx_zz[i] = -2.0 * g_0_0_xx_zz[i] * a_exp + 4.0 * g_xx_0_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_0_xy_xx, g_0_0_xy_xy, g_0_0_xy_xz, g_0_0_xy_yy, g_0_0_xy_yz, g_0_0_xy_zz, g_xx_0_0_0_0_0_xy_xx, g_xx_0_0_0_0_0_xy_xy, g_xx_0_0_0_0_0_xy_xz, g_xx_0_0_0_0_0_xy_yy, g_xx_0_0_0_0_0_xy_yz, g_xx_0_0_0_0_0_xy_zz, g_xx_0_xy_xx, g_xx_0_xy_xy, g_xx_0_xy_xz, g_xx_0_xy_yy, g_xx_0_xy_yz, g_xx_0_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_0_xy_xx[i] = -2.0 * g_0_0_xy_xx[i] * a_exp + 4.0 * g_xx_0_xy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_xy_xy[i] = -2.0 * g_0_0_xy_xy[i] * a_exp + 4.0 * g_xx_0_xy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_xy_xz[i] = -2.0 * g_0_0_xy_xz[i] * a_exp + 4.0 * g_xx_0_xy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_xy_yy[i] = -2.0 * g_0_0_xy_yy[i] * a_exp + 4.0 * g_xx_0_xy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_xy_yz[i] = -2.0 * g_0_0_xy_yz[i] * a_exp + 4.0 * g_xx_0_xy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_xy_zz[i] = -2.0 * g_0_0_xy_zz[i] * a_exp + 4.0 * g_xx_0_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_0_xz_xx, g_0_0_xz_xy, g_0_0_xz_xz, g_0_0_xz_yy, g_0_0_xz_yz, g_0_0_xz_zz, g_xx_0_0_0_0_0_xz_xx, g_xx_0_0_0_0_0_xz_xy, g_xx_0_0_0_0_0_xz_xz, g_xx_0_0_0_0_0_xz_yy, g_xx_0_0_0_0_0_xz_yz, g_xx_0_0_0_0_0_xz_zz, g_xx_0_xz_xx, g_xx_0_xz_xy, g_xx_0_xz_xz, g_xx_0_xz_yy, g_xx_0_xz_yz, g_xx_0_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_0_xz_xx[i] = -2.0 * g_0_0_xz_xx[i] * a_exp + 4.0 * g_xx_0_xz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_xz_xy[i] = -2.0 * g_0_0_xz_xy[i] * a_exp + 4.0 * g_xx_0_xz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_xz_xz[i] = -2.0 * g_0_0_xz_xz[i] * a_exp + 4.0 * g_xx_0_xz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_xz_yy[i] = -2.0 * g_0_0_xz_yy[i] * a_exp + 4.0 * g_xx_0_xz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_xz_yz[i] = -2.0 * g_0_0_xz_yz[i] * a_exp + 4.0 * g_xx_0_xz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_xz_zz[i] = -2.0 * g_0_0_xz_zz[i] * a_exp + 4.0 * g_xx_0_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_0_yy_xx, g_0_0_yy_xy, g_0_0_yy_xz, g_0_0_yy_yy, g_0_0_yy_yz, g_0_0_yy_zz, g_xx_0_0_0_0_0_yy_xx, g_xx_0_0_0_0_0_yy_xy, g_xx_0_0_0_0_0_yy_xz, g_xx_0_0_0_0_0_yy_yy, g_xx_0_0_0_0_0_yy_yz, g_xx_0_0_0_0_0_yy_zz, g_xx_0_yy_xx, g_xx_0_yy_xy, g_xx_0_yy_xz, g_xx_0_yy_yy, g_xx_0_yy_yz, g_xx_0_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_0_yy_xx[i] = -2.0 * g_0_0_yy_xx[i] * a_exp + 4.0 * g_xx_0_yy_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_yy_xy[i] = -2.0 * g_0_0_yy_xy[i] * a_exp + 4.0 * g_xx_0_yy_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_yy_xz[i] = -2.0 * g_0_0_yy_xz[i] * a_exp + 4.0 * g_xx_0_yy_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_yy_yy[i] = -2.0 * g_0_0_yy_yy[i] * a_exp + 4.0 * g_xx_0_yy_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_yy_yz[i] = -2.0 * g_0_0_yy_yz[i] * a_exp + 4.0 * g_xx_0_yy_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_yy_zz[i] = -2.0 * g_0_0_yy_zz[i] * a_exp + 4.0 * g_xx_0_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_0_yz_xx, g_0_0_yz_xy, g_0_0_yz_xz, g_0_0_yz_yy, g_0_0_yz_yz, g_0_0_yz_zz, g_xx_0_0_0_0_0_yz_xx, g_xx_0_0_0_0_0_yz_xy, g_xx_0_0_0_0_0_yz_xz, g_xx_0_0_0_0_0_yz_yy, g_xx_0_0_0_0_0_yz_yz, g_xx_0_0_0_0_0_yz_zz, g_xx_0_yz_xx, g_xx_0_yz_xy, g_xx_0_yz_xz, g_xx_0_yz_yy, g_xx_0_yz_yz, g_xx_0_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_0_yz_xx[i] = -2.0 * g_0_0_yz_xx[i] * a_exp + 4.0 * g_xx_0_yz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_yz_xy[i] = -2.0 * g_0_0_yz_xy[i] * a_exp + 4.0 * g_xx_0_yz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_yz_xz[i] = -2.0 * g_0_0_yz_xz[i] * a_exp + 4.0 * g_xx_0_yz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_yz_yy[i] = -2.0 * g_0_0_yz_yy[i] * a_exp + 4.0 * g_xx_0_yz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_yz_yz[i] = -2.0 * g_0_0_yz_yz[i] * a_exp + 4.0 * g_xx_0_yz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_yz_zz[i] = -2.0 * g_0_0_yz_zz[i] * a_exp + 4.0 * g_xx_0_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_0_zz_xx, g_0_0_zz_xy, g_0_0_zz_xz, g_0_0_zz_yy, g_0_0_zz_yz, g_0_0_zz_zz, g_xx_0_0_0_0_0_zz_xx, g_xx_0_0_0_0_0_zz_xy, g_xx_0_0_0_0_0_zz_xz, g_xx_0_0_0_0_0_zz_yy, g_xx_0_0_0_0_0_zz_yz, g_xx_0_0_0_0_0_zz_zz, g_xx_0_zz_xx, g_xx_0_zz_xy, g_xx_0_zz_xz, g_xx_0_zz_yy, g_xx_0_zz_yz, g_xx_0_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_0_zz_xx[i] = -2.0 * g_0_0_zz_xx[i] * a_exp + 4.0 * g_xx_0_zz_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_zz_xy[i] = -2.0 * g_0_0_zz_xy[i] * a_exp + 4.0 * g_xx_0_zz_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_zz_xz[i] = -2.0 * g_0_0_zz_xz[i] * a_exp + 4.0 * g_xx_0_zz_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_zz_yy[i] = -2.0 * g_0_0_zz_yy[i] * a_exp + 4.0 * g_xx_0_zz_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_zz_yz[i] = -2.0 * g_0_0_zz_yz[i] * a_exp + 4.0 * g_xx_0_zz_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_zz_zz[i] = -2.0 * g_0_0_zz_zz[i] * a_exp + 4.0 * g_xx_0_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_xy_0_0_0_0_0_xx_xx, g_xy_0_0_0_0_0_xx_xy, g_xy_0_0_0_0_0_xx_xz, g_xy_0_0_0_0_0_xx_yy, g_xy_0_0_0_0_0_xx_yz, g_xy_0_0_0_0_0_xx_zz, g_xy_0_xx_xx, g_xy_0_xx_xy, g_xy_0_xx_xz, g_xy_0_xx_yy, g_xy_0_xx_yz, g_xy_0_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_0_xx_xx[i] = 4.0 * g_xy_0_xx_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_xx_xy[i] = 4.0 * g_xy_0_xx_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_xx_xz[i] = 4.0 * g_xy_0_xx_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_xx_yy[i] = 4.0 * g_xy_0_xx_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_xx_yz[i] = 4.0 * g_xy_0_xx_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_xx_zz[i] = 4.0 * g_xy_0_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_xy_0_0_0_0_0_xy_xx, g_xy_0_0_0_0_0_xy_xy, g_xy_0_0_0_0_0_xy_xz, g_xy_0_0_0_0_0_xy_yy, g_xy_0_0_0_0_0_xy_yz, g_xy_0_0_0_0_0_xy_zz, g_xy_0_xy_xx, g_xy_0_xy_xy, g_xy_0_xy_xz, g_xy_0_xy_yy, g_xy_0_xy_yz, g_xy_0_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_0_xy_xx[i] = 4.0 * g_xy_0_xy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_xy_xy[i] = 4.0 * g_xy_0_xy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_xy_xz[i] = 4.0 * g_xy_0_xy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_xy_yy[i] = 4.0 * g_xy_0_xy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_xy_yz[i] = 4.0 * g_xy_0_xy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_xy_zz[i] = 4.0 * g_xy_0_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_xy_0_0_0_0_0_xz_xx, g_xy_0_0_0_0_0_xz_xy, g_xy_0_0_0_0_0_xz_xz, g_xy_0_0_0_0_0_xz_yy, g_xy_0_0_0_0_0_xz_yz, g_xy_0_0_0_0_0_xz_zz, g_xy_0_xz_xx, g_xy_0_xz_xy, g_xy_0_xz_xz, g_xy_0_xz_yy, g_xy_0_xz_yz, g_xy_0_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_0_xz_xx[i] = 4.0 * g_xy_0_xz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_xz_xy[i] = 4.0 * g_xy_0_xz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_xz_xz[i] = 4.0 * g_xy_0_xz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_xz_yy[i] = 4.0 * g_xy_0_xz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_xz_yz[i] = 4.0 * g_xy_0_xz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_xz_zz[i] = 4.0 * g_xy_0_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_xy_0_0_0_0_0_yy_xx, g_xy_0_0_0_0_0_yy_xy, g_xy_0_0_0_0_0_yy_xz, g_xy_0_0_0_0_0_yy_yy, g_xy_0_0_0_0_0_yy_yz, g_xy_0_0_0_0_0_yy_zz, g_xy_0_yy_xx, g_xy_0_yy_xy, g_xy_0_yy_xz, g_xy_0_yy_yy, g_xy_0_yy_yz, g_xy_0_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_0_yy_xx[i] = 4.0 * g_xy_0_yy_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_yy_xy[i] = 4.0 * g_xy_0_yy_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_yy_xz[i] = 4.0 * g_xy_0_yy_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_yy_yy[i] = 4.0 * g_xy_0_yy_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_yy_yz[i] = 4.0 * g_xy_0_yy_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_yy_zz[i] = 4.0 * g_xy_0_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_xy_0_0_0_0_0_yz_xx, g_xy_0_0_0_0_0_yz_xy, g_xy_0_0_0_0_0_yz_xz, g_xy_0_0_0_0_0_yz_yy, g_xy_0_0_0_0_0_yz_yz, g_xy_0_0_0_0_0_yz_zz, g_xy_0_yz_xx, g_xy_0_yz_xy, g_xy_0_yz_xz, g_xy_0_yz_yy, g_xy_0_yz_yz, g_xy_0_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_0_yz_xx[i] = 4.0 * g_xy_0_yz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_yz_xy[i] = 4.0 * g_xy_0_yz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_yz_xz[i] = 4.0 * g_xy_0_yz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_yz_yy[i] = 4.0 * g_xy_0_yz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_yz_yz[i] = 4.0 * g_xy_0_yz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_yz_zz[i] = 4.0 * g_xy_0_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_xy_0_0_0_0_0_zz_xx, g_xy_0_0_0_0_0_zz_xy, g_xy_0_0_0_0_0_zz_xz, g_xy_0_0_0_0_0_zz_yy, g_xy_0_0_0_0_0_zz_yz, g_xy_0_0_0_0_0_zz_zz, g_xy_0_zz_xx, g_xy_0_zz_xy, g_xy_0_zz_xz, g_xy_0_zz_yy, g_xy_0_zz_yz, g_xy_0_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_0_zz_xx[i] = 4.0 * g_xy_0_zz_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_zz_xy[i] = 4.0 * g_xy_0_zz_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_zz_xz[i] = 4.0 * g_xy_0_zz_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_zz_yy[i] = 4.0 * g_xy_0_zz_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_zz_yz[i] = 4.0 * g_xy_0_zz_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_zz_zz[i] = 4.0 * g_xy_0_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_xz_0_0_0_0_0_xx_xx, g_xz_0_0_0_0_0_xx_xy, g_xz_0_0_0_0_0_xx_xz, g_xz_0_0_0_0_0_xx_yy, g_xz_0_0_0_0_0_xx_yz, g_xz_0_0_0_0_0_xx_zz, g_xz_0_xx_xx, g_xz_0_xx_xy, g_xz_0_xx_xz, g_xz_0_xx_yy, g_xz_0_xx_yz, g_xz_0_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_0_xx_xx[i] = 4.0 * g_xz_0_xx_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_xx_xy[i] = 4.0 * g_xz_0_xx_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_xx_xz[i] = 4.0 * g_xz_0_xx_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_xx_yy[i] = 4.0 * g_xz_0_xx_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_xx_yz[i] = 4.0 * g_xz_0_xx_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_xx_zz[i] = 4.0 * g_xz_0_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_xz_0_0_0_0_0_xy_xx, g_xz_0_0_0_0_0_xy_xy, g_xz_0_0_0_0_0_xy_xz, g_xz_0_0_0_0_0_xy_yy, g_xz_0_0_0_0_0_xy_yz, g_xz_0_0_0_0_0_xy_zz, g_xz_0_xy_xx, g_xz_0_xy_xy, g_xz_0_xy_xz, g_xz_0_xy_yy, g_xz_0_xy_yz, g_xz_0_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_0_xy_xx[i] = 4.0 * g_xz_0_xy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_xy_xy[i] = 4.0 * g_xz_0_xy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_xy_xz[i] = 4.0 * g_xz_0_xy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_xy_yy[i] = 4.0 * g_xz_0_xy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_xy_yz[i] = 4.0 * g_xz_0_xy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_xy_zz[i] = 4.0 * g_xz_0_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_xz_0_0_0_0_0_xz_xx, g_xz_0_0_0_0_0_xz_xy, g_xz_0_0_0_0_0_xz_xz, g_xz_0_0_0_0_0_xz_yy, g_xz_0_0_0_0_0_xz_yz, g_xz_0_0_0_0_0_xz_zz, g_xz_0_xz_xx, g_xz_0_xz_xy, g_xz_0_xz_xz, g_xz_0_xz_yy, g_xz_0_xz_yz, g_xz_0_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_0_xz_xx[i] = 4.0 * g_xz_0_xz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_xz_xy[i] = 4.0 * g_xz_0_xz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_xz_xz[i] = 4.0 * g_xz_0_xz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_xz_yy[i] = 4.0 * g_xz_0_xz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_xz_yz[i] = 4.0 * g_xz_0_xz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_xz_zz[i] = 4.0 * g_xz_0_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_xz_0_0_0_0_0_yy_xx, g_xz_0_0_0_0_0_yy_xy, g_xz_0_0_0_0_0_yy_xz, g_xz_0_0_0_0_0_yy_yy, g_xz_0_0_0_0_0_yy_yz, g_xz_0_0_0_0_0_yy_zz, g_xz_0_yy_xx, g_xz_0_yy_xy, g_xz_0_yy_xz, g_xz_0_yy_yy, g_xz_0_yy_yz, g_xz_0_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_0_yy_xx[i] = 4.0 * g_xz_0_yy_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_yy_xy[i] = 4.0 * g_xz_0_yy_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_yy_xz[i] = 4.0 * g_xz_0_yy_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_yy_yy[i] = 4.0 * g_xz_0_yy_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_yy_yz[i] = 4.0 * g_xz_0_yy_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_yy_zz[i] = 4.0 * g_xz_0_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_xz_0_0_0_0_0_yz_xx, g_xz_0_0_0_0_0_yz_xy, g_xz_0_0_0_0_0_yz_xz, g_xz_0_0_0_0_0_yz_yy, g_xz_0_0_0_0_0_yz_yz, g_xz_0_0_0_0_0_yz_zz, g_xz_0_yz_xx, g_xz_0_yz_xy, g_xz_0_yz_xz, g_xz_0_yz_yy, g_xz_0_yz_yz, g_xz_0_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_0_yz_xx[i] = 4.0 * g_xz_0_yz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_yz_xy[i] = 4.0 * g_xz_0_yz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_yz_xz[i] = 4.0 * g_xz_0_yz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_yz_yy[i] = 4.0 * g_xz_0_yz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_yz_yz[i] = 4.0 * g_xz_0_yz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_yz_zz[i] = 4.0 * g_xz_0_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_xz_0_0_0_0_0_zz_xx, g_xz_0_0_0_0_0_zz_xy, g_xz_0_0_0_0_0_zz_xz, g_xz_0_0_0_0_0_zz_yy, g_xz_0_0_0_0_0_zz_yz, g_xz_0_0_0_0_0_zz_zz, g_xz_0_zz_xx, g_xz_0_zz_xy, g_xz_0_zz_xz, g_xz_0_zz_yy, g_xz_0_zz_yz, g_xz_0_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_0_zz_xx[i] = 4.0 * g_xz_0_zz_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_zz_xy[i] = 4.0 * g_xz_0_zz_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_zz_xz[i] = 4.0 * g_xz_0_zz_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_zz_yy[i] = 4.0 * g_xz_0_zz_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_zz_yz[i] = 4.0 * g_xz_0_zz_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_zz_zz[i] = 4.0 * g_xz_0_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_0_0_xx_xx, g_0_0_xx_xy, g_0_0_xx_xz, g_0_0_xx_yy, g_0_0_xx_yz, g_0_0_xx_zz, g_yy_0_0_0_0_0_xx_xx, g_yy_0_0_0_0_0_xx_xy, g_yy_0_0_0_0_0_xx_xz, g_yy_0_0_0_0_0_xx_yy, g_yy_0_0_0_0_0_xx_yz, g_yy_0_0_0_0_0_xx_zz, g_yy_0_xx_xx, g_yy_0_xx_xy, g_yy_0_xx_xz, g_yy_0_xx_yy, g_yy_0_xx_yz, g_yy_0_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_0_xx_xx[i] = -2.0 * g_0_0_xx_xx[i] * a_exp + 4.0 * g_yy_0_xx_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_xx_xy[i] = -2.0 * g_0_0_xx_xy[i] * a_exp + 4.0 * g_yy_0_xx_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_xx_xz[i] = -2.0 * g_0_0_xx_xz[i] * a_exp + 4.0 * g_yy_0_xx_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_xx_yy[i] = -2.0 * g_0_0_xx_yy[i] * a_exp + 4.0 * g_yy_0_xx_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_xx_yz[i] = -2.0 * g_0_0_xx_yz[i] * a_exp + 4.0 * g_yy_0_xx_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_xx_zz[i] = -2.0 * g_0_0_xx_zz[i] * a_exp + 4.0 * g_yy_0_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_0_0_xy_xx, g_0_0_xy_xy, g_0_0_xy_xz, g_0_0_xy_yy, g_0_0_xy_yz, g_0_0_xy_zz, g_yy_0_0_0_0_0_xy_xx, g_yy_0_0_0_0_0_xy_xy, g_yy_0_0_0_0_0_xy_xz, g_yy_0_0_0_0_0_xy_yy, g_yy_0_0_0_0_0_xy_yz, g_yy_0_0_0_0_0_xy_zz, g_yy_0_xy_xx, g_yy_0_xy_xy, g_yy_0_xy_xz, g_yy_0_xy_yy, g_yy_0_xy_yz, g_yy_0_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_0_xy_xx[i] = -2.0 * g_0_0_xy_xx[i] * a_exp + 4.0 * g_yy_0_xy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_xy_xy[i] = -2.0 * g_0_0_xy_xy[i] * a_exp + 4.0 * g_yy_0_xy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_xy_xz[i] = -2.0 * g_0_0_xy_xz[i] * a_exp + 4.0 * g_yy_0_xy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_xy_yy[i] = -2.0 * g_0_0_xy_yy[i] * a_exp + 4.0 * g_yy_0_xy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_xy_yz[i] = -2.0 * g_0_0_xy_yz[i] * a_exp + 4.0 * g_yy_0_xy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_xy_zz[i] = -2.0 * g_0_0_xy_zz[i] * a_exp + 4.0 * g_yy_0_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_0_0_xz_xx, g_0_0_xz_xy, g_0_0_xz_xz, g_0_0_xz_yy, g_0_0_xz_yz, g_0_0_xz_zz, g_yy_0_0_0_0_0_xz_xx, g_yy_0_0_0_0_0_xz_xy, g_yy_0_0_0_0_0_xz_xz, g_yy_0_0_0_0_0_xz_yy, g_yy_0_0_0_0_0_xz_yz, g_yy_0_0_0_0_0_xz_zz, g_yy_0_xz_xx, g_yy_0_xz_xy, g_yy_0_xz_xz, g_yy_0_xz_yy, g_yy_0_xz_yz, g_yy_0_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_0_xz_xx[i] = -2.0 * g_0_0_xz_xx[i] * a_exp + 4.0 * g_yy_0_xz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_xz_xy[i] = -2.0 * g_0_0_xz_xy[i] * a_exp + 4.0 * g_yy_0_xz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_xz_xz[i] = -2.0 * g_0_0_xz_xz[i] * a_exp + 4.0 * g_yy_0_xz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_xz_yy[i] = -2.0 * g_0_0_xz_yy[i] * a_exp + 4.0 * g_yy_0_xz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_xz_yz[i] = -2.0 * g_0_0_xz_yz[i] * a_exp + 4.0 * g_yy_0_xz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_xz_zz[i] = -2.0 * g_0_0_xz_zz[i] * a_exp + 4.0 * g_yy_0_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_0_0_yy_xx, g_0_0_yy_xy, g_0_0_yy_xz, g_0_0_yy_yy, g_0_0_yy_yz, g_0_0_yy_zz, g_yy_0_0_0_0_0_yy_xx, g_yy_0_0_0_0_0_yy_xy, g_yy_0_0_0_0_0_yy_xz, g_yy_0_0_0_0_0_yy_yy, g_yy_0_0_0_0_0_yy_yz, g_yy_0_0_0_0_0_yy_zz, g_yy_0_yy_xx, g_yy_0_yy_xy, g_yy_0_yy_xz, g_yy_0_yy_yy, g_yy_0_yy_yz, g_yy_0_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_0_yy_xx[i] = -2.0 * g_0_0_yy_xx[i] * a_exp + 4.0 * g_yy_0_yy_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_yy_xy[i] = -2.0 * g_0_0_yy_xy[i] * a_exp + 4.0 * g_yy_0_yy_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_yy_xz[i] = -2.0 * g_0_0_yy_xz[i] * a_exp + 4.0 * g_yy_0_yy_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_yy_yy[i] = -2.0 * g_0_0_yy_yy[i] * a_exp + 4.0 * g_yy_0_yy_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_yy_yz[i] = -2.0 * g_0_0_yy_yz[i] * a_exp + 4.0 * g_yy_0_yy_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_yy_zz[i] = -2.0 * g_0_0_yy_zz[i] * a_exp + 4.0 * g_yy_0_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_0_0_yz_xx, g_0_0_yz_xy, g_0_0_yz_xz, g_0_0_yz_yy, g_0_0_yz_yz, g_0_0_yz_zz, g_yy_0_0_0_0_0_yz_xx, g_yy_0_0_0_0_0_yz_xy, g_yy_0_0_0_0_0_yz_xz, g_yy_0_0_0_0_0_yz_yy, g_yy_0_0_0_0_0_yz_yz, g_yy_0_0_0_0_0_yz_zz, g_yy_0_yz_xx, g_yy_0_yz_xy, g_yy_0_yz_xz, g_yy_0_yz_yy, g_yy_0_yz_yz, g_yy_0_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_0_yz_xx[i] = -2.0 * g_0_0_yz_xx[i] * a_exp + 4.0 * g_yy_0_yz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_yz_xy[i] = -2.0 * g_0_0_yz_xy[i] * a_exp + 4.0 * g_yy_0_yz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_yz_xz[i] = -2.0 * g_0_0_yz_xz[i] * a_exp + 4.0 * g_yy_0_yz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_yz_yy[i] = -2.0 * g_0_0_yz_yy[i] * a_exp + 4.0 * g_yy_0_yz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_yz_yz[i] = -2.0 * g_0_0_yz_yz[i] * a_exp + 4.0 * g_yy_0_yz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_yz_zz[i] = -2.0 * g_0_0_yz_zz[i] * a_exp + 4.0 * g_yy_0_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_0_0_zz_xx, g_0_0_zz_xy, g_0_0_zz_xz, g_0_0_zz_yy, g_0_0_zz_yz, g_0_0_zz_zz, g_yy_0_0_0_0_0_zz_xx, g_yy_0_0_0_0_0_zz_xy, g_yy_0_0_0_0_0_zz_xz, g_yy_0_0_0_0_0_zz_yy, g_yy_0_0_0_0_0_zz_yz, g_yy_0_0_0_0_0_zz_zz, g_yy_0_zz_xx, g_yy_0_zz_xy, g_yy_0_zz_xz, g_yy_0_zz_yy, g_yy_0_zz_yz, g_yy_0_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_0_zz_xx[i] = -2.0 * g_0_0_zz_xx[i] * a_exp + 4.0 * g_yy_0_zz_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_zz_xy[i] = -2.0 * g_0_0_zz_xy[i] * a_exp + 4.0 * g_yy_0_zz_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_zz_xz[i] = -2.0 * g_0_0_zz_xz[i] * a_exp + 4.0 * g_yy_0_zz_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_zz_yy[i] = -2.0 * g_0_0_zz_yy[i] * a_exp + 4.0 * g_yy_0_zz_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_zz_yz[i] = -2.0 * g_0_0_zz_yz[i] * a_exp + 4.0 * g_yy_0_zz_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_zz_zz[i] = -2.0 * g_0_0_zz_zz[i] * a_exp + 4.0 * g_yy_0_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_yz_0_0_0_0_0_xx_xx, g_yz_0_0_0_0_0_xx_xy, g_yz_0_0_0_0_0_xx_xz, g_yz_0_0_0_0_0_xx_yy, g_yz_0_0_0_0_0_xx_yz, g_yz_0_0_0_0_0_xx_zz, g_yz_0_xx_xx, g_yz_0_xx_xy, g_yz_0_xx_xz, g_yz_0_xx_yy, g_yz_0_xx_yz, g_yz_0_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_0_xx_xx[i] = 4.0 * g_yz_0_xx_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_xx_xy[i] = 4.0 * g_yz_0_xx_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_xx_xz[i] = 4.0 * g_yz_0_xx_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_xx_yy[i] = 4.0 * g_yz_0_xx_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_xx_yz[i] = 4.0 * g_yz_0_xx_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_xx_zz[i] = 4.0 * g_yz_0_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_yz_0_0_0_0_0_xy_xx, g_yz_0_0_0_0_0_xy_xy, g_yz_0_0_0_0_0_xy_xz, g_yz_0_0_0_0_0_xy_yy, g_yz_0_0_0_0_0_xy_yz, g_yz_0_0_0_0_0_xy_zz, g_yz_0_xy_xx, g_yz_0_xy_xy, g_yz_0_xy_xz, g_yz_0_xy_yy, g_yz_0_xy_yz, g_yz_0_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_0_xy_xx[i] = 4.0 * g_yz_0_xy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_xy_xy[i] = 4.0 * g_yz_0_xy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_xy_xz[i] = 4.0 * g_yz_0_xy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_xy_yy[i] = 4.0 * g_yz_0_xy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_xy_yz[i] = 4.0 * g_yz_0_xy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_xy_zz[i] = 4.0 * g_yz_0_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_yz_0_0_0_0_0_xz_xx, g_yz_0_0_0_0_0_xz_xy, g_yz_0_0_0_0_0_xz_xz, g_yz_0_0_0_0_0_xz_yy, g_yz_0_0_0_0_0_xz_yz, g_yz_0_0_0_0_0_xz_zz, g_yz_0_xz_xx, g_yz_0_xz_xy, g_yz_0_xz_xz, g_yz_0_xz_yy, g_yz_0_xz_yz, g_yz_0_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_0_xz_xx[i] = 4.0 * g_yz_0_xz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_xz_xy[i] = 4.0 * g_yz_0_xz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_xz_xz[i] = 4.0 * g_yz_0_xz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_xz_yy[i] = 4.0 * g_yz_0_xz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_xz_yz[i] = 4.0 * g_yz_0_xz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_xz_zz[i] = 4.0 * g_yz_0_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_yz_0_0_0_0_0_yy_xx, g_yz_0_0_0_0_0_yy_xy, g_yz_0_0_0_0_0_yy_xz, g_yz_0_0_0_0_0_yy_yy, g_yz_0_0_0_0_0_yy_yz, g_yz_0_0_0_0_0_yy_zz, g_yz_0_yy_xx, g_yz_0_yy_xy, g_yz_0_yy_xz, g_yz_0_yy_yy, g_yz_0_yy_yz, g_yz_0_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_0_yy_xx[i] = 4.0 * g_yz_0_yy_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_yy_xy[i] = 4.0 * g_yz_0_yy_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_yy_xz[i] = 4.0 * g_yz_0_yy_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_yy_yy[i] = 4.0 * g_yz_0_yy_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_yy_yz[i] = 4.0 * g_yz_0_yy_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_yy_zz[i] = 4.0 * g_yz_0_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_yz_0_0_0_0_0_yz_xx, g_yz_0_0_0_0_0_yz_xy, g_yz_0_0_0_0_0_yz_xz, g_yz_0_0_0_0_0_yz_yy, g_yz_0_0_0_0_0_yz_yz, g_yz_0_0_0_0_0_yz_zz, g_yz_0_yz_xx, g_yz_0_yz_xy, g_yz_0_yz_xz, g_yz_0_yz_yy, g_yz_0_yz_yz, g_yz_0_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_0_yz_xx[i] = 4.0 * g_yz_0_yz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_yz_xy[i] = 4.0 * g_yz_0_yz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_yz_xz[i] = 4.0 * g_yz_0_yz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_yz_yy[i] = 4.0 * g_yz_0_yz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_yz_yz[i] = 4.0 * g_yz_0_yz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_yz_zz[i] = 4.0 * g_yz_0_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_yz_0_0_0_0_0_zz_xx, g_yz_0_0_0_0_0_zz_xy, g_yz_0_0_0_0_0_zz_xz, g_yz_0_0_0_0_0_zz_yy, g_yz_0_0_0_0_0_zz_yz, g_yz_0_0_0_0_0_zz_zz, g_yz_0_zz_xx, g_yz_0_zz_xy, g_yz_0_zz_xz, g_yz_0_zz_yy, g_yz_0_zz_yz, g_yz_0_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_0_zz_xx[i] = 4.0 * g_yz_0_zz_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_zz_xy[i] = 4.0 * g_yz_0_zz_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_zz_xz[i] = 4.0 * g_yz_0_zz_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_zz_yy[i] = 4.0 * g_yz_0_zz_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_zz_yz[i] = 4.0 * g_yz_0_zz_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_zz_zz[i] = 4.0 * g_yz_0_zz_zz[i] * a_exp * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_0_0_xx_xx, g_0_0_xx_xy, g_0_0_xx_xz, g_0_0_xx_yy, g_0_0_xx_yz, g_0_0_xx_zz, g_zz_0_0_0_0_0_xx_xx, g_zz_0_0_0_0_0_xx_xy, g_zz_0_0_0_0_0_xx_xz, g_zz_0_0_0_0_0_xx_yy, g_zz_0_0_0_0_0_xx_yz, g_zz_0_0_0_0_0_xx_zz, g_zz_0_xx_xx, g_zz_0_xx_xy, g_zz_0_xx_xz, g_zz_0_xx_yy, g_zz_0_xx_yz, g_zz_0_xx_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_0_xx_xx[i] = -2.0 * g_0_0_xx_xx[i] * a_exp + 4.0 * g_zz_0_xx_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_xx_xy[i] = -2.0 * g_0_0_xx_xy[i] * a_exp + 4.0 * g_zz_0_xx_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_xx_xz[i] = -2.0 * g_0_0_xx_xz[i] * a_exp + 4.0 * g_zz_0_xx_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_xx_yy[i] = -2.0 * g_0_0_xx_yy[i] * a_exp + 4.0 * g_zz_0_xx_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_xx_yz[i] = -2.0 * g_0_0_xx_yz[i] * a_exp + 4.0 * g_zz_0_xx_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_xx_zz[i] = -2.0 * g_0_0_xx_zz[i] * a_exp + 4.0 * g_zz_0_xx_zz[i] * a_exp * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_0_0_xy_xx, g_0_0_xy_xy, g_0_0_xy_xz, g_0_0_xy_yy, g_0_0_xy_yz, g_0_0_xy_zz, g_zz_0_0_0_0_0_xy_xx, g_zz_0_0_0_0_0_xy_xy, g_zz_0_0_0_0_0_xy_xz, g_zz_0_0_0_0_0_xy_yy, g_zz_0_0_0_0_0_xy_yz, g_zz_0_0_0_0_0_xy_zz, g_zz_0_xy_xx, g_zz_0_xy_xy, g_zz_0_xy_xz, g_zz_0_xy_yy, g_zz_0_xy_yz, g_zz_0_xy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_0_xy_xx[i] = -2.0 * g_0_0_xy_xx[i] * a_exp + 4.0 * g_zz_0_xy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_xy_xy[i] = -2.0 * g_0_0_xy_xy[i] * a_exp + 4.0 * g_zz_0_xy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_xy_xz[i] = -2.0 * g_0_0_xy_xz[i] * a_exp + 4.0 * g_zz_0_xy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_xy_yy[i] = -2.0 * g_0_0_xy_yy[i] * a_exp + 4.0 * g_zz_0_xy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_xy_yz[i] = -2.0 * g_0_0_xy_yz[i] * a_exp + 4.0 * g_zz_0_xy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_xy_zz[i] = -2.0 * g_0_0_xy_zz[i] * a_exp + 4.0 * g_zz_0_xy_zz[i] * a_exp * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_0_0_xz_xx, g_0_0_xz_xy, g_0_0_xz_xz, g_0_0_xz_yy, g_0_0_xz_yz, g_0_0_xz_zz, g_zz_0_0_0_0_0_xz_xx, g_zz_0_0_0_0_0_xz_xy, g_zz_0_0_0_0_0_xz_xz, g_zz_0_0_0_0_0_xz_yy, g_zz_0_0_0_0_0_xz_yz, g_zz_0_0_0_0_0_xz_zz, g_zz_0_xz_xx, g_zz_0_xz_xy, g_zz_0_xz_xz, g_zz_0_xz_yy, g_zz_0_xz_yz, g_zz_0_xz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_0_xz_xx[i] = -2.0 * g_0_0_xz_xx[i] * a_exp + 4.0 * g_zz_0_xz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_xz_xy[i] = -2.0 * g_0_0_xz_xy[i] * a_exp + 4.0 * g_zz_0_xz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_xz_xz[i] = -2.0 * g_0_0_xz_xz[i] * a_exp + 4.0 * g_zz_0_xz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_xz_yy[i] = -2.0 * g_0_0_xz_yy[i] * a_exp + 4.0 * g_zz_0_xz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_xz_yz[i] = -2.0 * g_0_0_xz_yz[i] * a_exp + 4.0 * g_zz_0_xz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_xz_zz[i] = -2.0 * g_0_0_xz_zz[i] * a_exp + 4.0 * g_zz_0_xz_zz[i] * a_exp * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_0_0_yy_xx, g_0_0_yy_xy, g_0_0_yy_xz, g_0_0_yy_yy, g_0_0_yy_yz, g_0_0_yy_zz, g_zz_0_0_0_0_0_yy_xx, g_zz_0_0_0_0_0_yy_xy, g_zz_0_0_0_0_0_yy_xz, g_zz_0_0_0_0_0_yy_yy, g_zz_0_0_0_0_0_yy_yz, g_zz_0_0_0_0_0_yy_zz, g_zz_0_yy_xx, g_zz_0_yy_xy, g_zz_0_yy_xz, g_zz_0_yy_yy, g_zz_0_yy_yz, g_zz_0_yy_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_0_yy_xx[i] = -2.0 * g_0_0_yy_xx[i] * a_exp + 4.0 * g_zz_0_yy_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_yy_xy[i] = -2.0 * g_0_0_yy_xy[i] * a_exp + 4.0 * g_zz_0_yy_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_yy_xz[i] = -2.0 * g_0_0_yy_xz[i] * a_exp + 4.0 * g_zz_0_yy_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_yy_yy[i] = -2.0 * g_0_0_yy_yy[i] * a_exp + 4.0 * g_zz_0_yy_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_yy_yz[i] = -2.0 * g_0_0_yy_yz[i] * a_exp + 4.0 * g_zz_0_yy_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_yy_zz[i] = -2.0 * g_0_0_yy_zz[i] * a_exp + 4.0 * g_zz_0_yy_zz[i] * a_exp * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_0_0_yz_xx, g_0_0_yz_xy, g_0_0_yz_xz, g_0_0_yz_yy, g_0_0_yz_yz, g_0_0_yz_zz, g_zz_0_0_0_0_0_yz_xx, g_zz_0_0_0_0_0_yz_xy, g_zz_0_0_0_0_0_yz_xz, g_zz_0_0_0_0_0_yz_yy, g_zz_0_0_0_0_0_yz_yz, g_zz_0_0_0_0_0_yz_zz, g_zz_0_yz_xx, g_zz_0_yz_xy, g_zz_0_yz_xz, g_zz_0_yz_yy, g_zz_0_yz_yz, g_zz_0_yz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_0_yz_xx[i] = -2.0 * g_0_0_yz_xx[i] * a_exp + 4.0 * g_zz_0_yz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_yz_xy[i] = -2.0 * g_0_0_yz_xy[i] * a_exp + 4.0 * g_zz_0_yz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_yz_xz[i] = -2.0 * g_0_0_yz_xz[i] * a_exp + 4.0 * g_zz_0_yz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_yz_yy[i] = -2.0 * g_0_0_yz_yy[i] * a_exp + 4.0 * g_zz_0_yz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_yz_yz[i] = -2.0 * g_0_0_yz_yz[i] * a_exp + 4.0 * g_zz_0_yz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_yz_zz[i] = -2.0 * g_0_0_yz_zz[i] * a_exp + 4.0 * g_zz_0_yz_zz[i] * a_exp * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_0_0_zz_xx, g_0_0_zz_xy, g_0_0_zz_xz, g_0_0_zz_yy, g_0_0_zz_yz, g_0_0_zz_zz, g_zz_0_0_0_0_0_zz_xx, g_zz_0_0_0_0_0_zz_xy, g_zz_0_0_0_0_0_zz_xz, g_zz_0_0_0_0_0_zz_yy, g_zz_0_0_0_0_0_zz_yz, g_zz_0_0_0_0_0_zz_zz, g_zz_0_zz_xx, g_zz_0_zz_xy, g_zz_0_zz_xz, g_zz_0_zz_yy, g_zz_0_zz_yz, g_zz_0_zz_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_0_zz_xx[i] = -2.0 * g_0_0_zz_xx[i] * a_exp + 4.0 * g_zz_0_zz_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_zz_xy[i] = -2.0 * g_0_0_zz_xy[i] * a_exp + 4.0 * g_zz_0_zz_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_zz_xz[i] = -2.0 * g_0_0_zz_xz[i] * a_exp + 4.0 * g_zz_0_zz_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_zz_yy[i] = -2.0 * g_0_0_zz_yy[i] * a_exp + 4.0 * g_zz_0_zz_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_zz_yz[i] = -2.0 * g_0_0_zz_yz[i] * a_exp + 4.0 * g_zz_0_zz_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_zz_zz[i] = -2.0 * g_0_0_zz_zz[i] * a_exp + 4.0 * g_zz_0_zz_zz[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

