#include "GeomDeriv2000OfScalarForDDSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_ddsd_0(CSimdArray<double>& buffer_2000_ddsd,
                     const CSimdArray<double>& buffer_sdsd,
                     const CSimdArray<double>& buffer_ddsd,
                     const CSimdArray<double>& buffer_gdsd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_ddsd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sdsd

    auto g_0_xx_0_xx = buffer_sdsd[0];

    auto g_0_xx_0_xy = buffer_sdsd[1];

    auto g_0_xx_0_xz = buffer_sdsd[2];

    auto g_0_xx_0_yy = buffer_sdsd[3];

    auto g_0_xx_0_yz = buffer_sdsd[4];

    auto g_0_xx_0_zz = buffer_sdsd[5];

    auto g_0_xy_0_xx = buffer_sdsd[6];

    auto g_0_xy_0_xy = buffer_sdsd[7];

    auto g_0_xy_0_xz = buffer_sdsd[8];

    auto g_0_xy_0_yy = buffer_sdsd[9];

    auto g_0_xy_0_yz = buffer_sdsd[10];

    auto g_0_xy_0_zz = buffer_sdsd[11];

    auto g_0_xz_0_xx = buffer_sdsd[12];

    auto g_0_xz_0_xy = buffer_sdsd[13];

    auto g_0_xz_0_xz = buffer_sdsd[14];

    auto g_0_xz_0_yy = buffer_sdsd[15];

    auto g_0_xz_0_yz = buffer_sdsd[16];

    auto g_0_xz_0_zz = buffer_sdsd[17];

    auto g_0_yy_0_xx = buffer_sdsd[18];

    auto g_0_yy_0_xy = buffer_sdsd[19];

    auto g_0_yy_0_xz = buffer_sdsd[20];

    auto g_0_yy_0_yy = buffer_sdsd[21];

    auto g_0_yy_0_yz = buffer_sdsd[22];

    auto g_0_yy_0_zz = buffer_sdsd[23];

    auto g_0_yz_0_xx = buffer_sdsd[24];

    auto g_0_yz_0_xy = buffer_sdsd[25];

    auto g_0_yz_0_xz = buffer_sdsd[26];

    auto g_0_yz_0_yy = buffer_sdsd[27];

    auto g_0_yz_0_yz = buffer_sdsd[28];

    auto g_0_yz_0_zz = buffer_sdsd[29];

    auto g_0_zz_0_xx = buffer_sdsd[30];

    auto g_0_zz_0_xy = buffer_sdsd[31];

    auto g_0_zz_0_xz = buffer_sdsd[32];

    auto g_0_zz_0_yy = buffer_sdsd[33];

    auto g_0_zz_0_yz = buffer_sdsd[34];

    auto g_0_zz_0_zz = buffer_sdsd[35];

    /// Set up components of auxilary buffer : buffer_ddsd

    auto g_xx_xx_0_xx = buffer_ddsd[0];

    auto g_xx_xx_0_xy = buffer_ddsd[1];

    auto g_xx_xx_0_xz = buffer_ddsd[2];

    auto g_xx_xx_0_yy = buffer_ddsd[3];

    auto g_xx_xx_0_yz = buffer_ddsd[4];

    auto g_xx_xx_0_zz = buffer_ddsd[5];

    auto g_xx_xy_0_xx = buffer_ddsd[6];

    auto g_xx_xy_0_xy = buffer_ddsd[7];

    auto g_xx_xy_0_xz = buffer_ddsd[8];

    auto g_xx_xy_0_yy = buffer_ddsd[9];

    auto g_xx_xy_0_yz = buffer_ddsd[10];

    auto g_xx_xy_0_zz = buffer_ddsd[11];

    auto g_xx_xz_0_xx = buffer_ddsd[12];

    auto g_xx_xz_0_xy = buffer_ddsd[13];

    auto g_xx_xz_0_xz = buffer_ddsd[14];

    auto g_xx_xz_0_yy = buffer_ddsd[15];

    auto g_xx_xz_0_yz = buffer_ddsd[16];

    auto g_xx_xz_0_zz = buffer_ddsd[17];

    auto g_xx_yy_0_xx = buffer_ddsd[18];

    auto g_xx_yy_0_xy = buffer_ddsd[19];

    auto g_xx_yy_0_xz = buffer_ddsd[20];

    auto g_xx_yy_0_yy = buffer_ddsd[21];

    auto g_xx_yy_0_yz = buffer_ddsd[22];

    auto g_xx_yy_0_zz = buffer_ddsd[23];

    auto g_xx_yz_0_xx = buffer_ddsd[24];

    auto g_xx_yz_0_xy = buffer_ddsd[25];

    auto g_xx_yz_0_xz = buffer_ddsd[26];

    auto g_xx_yz_0_yy = buffer_ddsd[27];

    auto g_xx_yz_0_yz = buffer_ddsd[28];

    auto g_xx_yz_0_zz = buffer_ddsd[29];

    auto g_xx_zz_0_xx = buffer_ddsd[30];

    auto g_xx_zz_0_xy = buffer_ddsd[31];

    auto g_xx_zz_0_xz = buffer_ddsd[32];

    auto g_xx_zz_0_yy = buffer_ddsd[33];

    auto g_xx_zz_0_yz = buffer_ddsd[34];

    auto g_xx_zz_0_zz = buffer_ddsd[35];

    auto g_xy_xx_0_xx = buffer_ddsd[36];

    auto g_xy_xx_0_xy = buffer_ddsd[37];

    auto g_xy_xx_0_xz = buffer_ddsd[38];

    auto g_xy_xx_0_yy = buffer_ddsd[39];

    auto g_xy_xx_0_yz = buffer_ddsd[40];

    auto g_xy_xx_0_zz = buffer_ddsd[41];

    auto g_xy_xy_0_xx = buffer_ddsd[42];

    auto g_xy_xy_0_xy = buffer_ddsd[43];

    auto g_xy_xy_0_xz = buffer_ddsd[44];

    auto g_xy_xy_0_yy = buffer_ddsd[45];

    auto g_xy_xy_0_yz = buffer_ddsd[46];

    auto g_xy_xy_0_zz = buffer_ddsd[47];

    auto g_xy_xz_0_xx = buffer_ddsd[48];

    auto g_xy_xz_0_xy = buffer_ddsd[49];

    auto g_xy_xz_0_xz = buffer_ddsd[50];

    auto g_xy_xz_0_yy = buffer_ddsd[51];

    auto g_xy_xz_0_yz = buffer_ddsd[52];

    auto g_xy_xz_0_zz = buffer_ddsd[53];

    auto g_xy_yy_0_xx = buffer_ddsd[54];

    auto g_xy_yy_0_xy = buffer_ddsd[55];

    auto g_xy_yy_0_xz = buffer_ddsd[56];

    auto g_xy_yy_0_yy = buffer_ddsd[57];

    auto g_xy_yy_0_yz = buffer_ddsd[58];

    auto g_xy_yy_0_zz = buffer_ddsd[59];

    auto g_xy_yz_0_xx = buffer_ddsd[60];

    auto g_xy_yz_0_xy = buffer_ddsd[61];

    auto g_xy_yz_0_xz = buffer_ddsd[62];

    auto g_xy_yz_0_yy = buffer_ddsd[63];

    auto g_xy_yz_0_yz = buffer_ddsd[64];

    auto g_xy_yz_0_zz = buffer_ddsd[65];

    auto g_xy_zz_0_xx = buffer_ddsd[66];

    auto g_xy_zz_0_xy = buffer_ddsd[67];

    auto g_xy_zz_0_xz = buffer_ddsd[68];

    auto g_xy_zz_0_yy = buffer_ddsd[69];

    auto g_xy_zz_0_yz = buffer_ddsd[70];

    auto g_xy_zz_0_zz = buffer_ddsd[71];

    auto g_xz_xx_0_xx = buffer_ddsd[72];

    auto g_xz_xx_0_xy = buffer_ddsd[73];

    auto g_xz_xx_0_xz = buffer_ddsd[74];

    auto g_xz_xx_0_yy = buffer_ddsd[75];

    auto g_xz_xx_0_yz = buffer_ddsd[76];

    auto g_xz_xx_0_zz = buffer_ddsd[77];

    auto g_xz_xy_0_xx = buffer_ddsd[78];

    auto g_xz_xy_0_xy = buffer_ddsd[79];

    auto g_xz_xy_0_xz = buffer_ddsd[80];

    auto g_xz_xy_0_yy = buffer_ddsd[81];

    auto g_xz_xy_0_yz = buffer_ddsd[82];

    auto g_xz_xy_0_zz = buffer_ddsd[83];

    auto g_xz_xz_0_xx = buffer_ddsd[84];

    auto g_xz_xz_0_xy = buffer_ddsd[85];

    auto g_xz_xz_0_xz = buffer_ddsd[86];

    auto g_xz_xz_0_yy = buffer_ddsd[87];

    auto g_xz_xz_0_yz = buffer_ddsd[88];

    auto g_xz_xz_0_zz = buffer_ddsd[89];

    auto g_xz_yy_0_xx = buffer_ddsd[90];

    auto g_xz_yy_0_xy = buffer_ddsd[91];

    auto g_xz_yy_0_xz = buffer_ddsd[92];

    auto g_xz_yy_0_yy = buffer_ddsd[93];

    auto g_xz_yy_0_yz = buffer_ddsd[94];

    auto g_xz_yy_0_zz = buffer_ddsd[95];

    auto g_xz_yz_0_xx = buffer_ddsd[96];

    auto g_xz_yz_0_xy = buffer_ddsd[97];

    auto g_xz_yz_0_xz = buffer_ddsd[98];

    auto g_xz_yz_0_yy = buffer_ddsd[99];

    auto g_xz_yz_0_yz = buffer_ddsd[100];

    auto g_xz_yz_0_zz = buffer_ddsd[101];

    auto g_xz_zz_0_xx = buffer_ddsd[102];

    auto g_xz_zz_0_xy = buffer_ddsd[103];

    auto g_xz_zz_0_xz = buffer_ddsd[104];

    auto g_xz_zz_0_yy = buffer_ddsd[105];

    auto g_xz_zz_0_yz = buffer_ddsd[106];

    auto g_xz_zz_0_zz = buffer_ddsd[107];

    auto g_yy_xx_0_xx = buffer_ddsd[108];

    auto g_yy_xx_0_xy = buffer_ddsd[109];

    auto g_yy_xx_0_xz = buffer_ddsd[110];

    auto g_yy_xx_0_yy = buffer_ddsd[111];

    auto g_yy_xx_0_yz = buffer_ddsd[112];

    auto g_yy_xx_0_zz = buffer_ddsd[113];

    auto g_yy_xy_0_xx = buffer_ddsd[114];

    auto g_yy_xy_0_xy = buffer_ddsd[115];

    auto g_yy_xy_0_xz = buffer_ddsd[116];

    auto g_yy_xy_0_yy = buffer_ddsd[117];

    auto g_yy_xy_0_yz = buffer_ddsd[118];

    auto g_yy_xy_0_zz = buffer_ddsd[119];

    auto g_yy_xz_0_xx = buffer_ddsd[120];

    auto g_yy_xz_0_xy = buffer_ddsd[121];

    auto g_yy_xz_0_xz = buffer_ddsd[122];

    auto g_yy_xz_0_yy = buffer_ddsd[123];

    auto g_yy_xz_0_yz = buffer_ddsd[124];

    auto g_yy_xz_0_zz = buffer_ddsd[125];

    auto g_yy_yy_0_xx = buffer_ddsd[126];

    auto g_yy_yy_0_xy = buffer_ddsd[127];

    auto g_yy_yy_0_xz = buffer_ddsd[128];

    auto g_yy_yy_0_yy = buffer_ddsd[129];

    auto g_yy_yy_0_yz = buffer_ddsd[130];

    auto g_yy_yy_0_zz = buffer_ddsd[131];

    auto g_yy_yz_0_xx = buffer_ddsd[132];

    auto g_yy_yz_0_xy = buffer_ddsd[133];

    auto g_yy_yz_0_xz = buffer_ddsd[134];

    auto g_yy_yz_0_yy = buffer_ddsd[135];

    auto g_yy_yz_0_yz = buffer_ddsd[136];

    auto g_yy_yz_0_zz = buffer_ddsd[137];

    auto g_yy_zz_0_xx = buffer_ddsd[138];

    auto g_yy_zz_0_xy = buffer_ddsd[139];

    auto g_yy_zz_0_xz = buffer_ddsd[140];

    auto g_yy_zz_0_yy = buffer_ddsd[141];

    auto g_yy_zz_0_yz = buffer_ddsd[142];

    auto g_yy_zz_0_zz = buffer_ddsd[143];

    auto g_yz_xx_0_xx = buffer_ddsd[144];

    auto g_yz_xx_0_xy = buffer_ddsd[145];

    auto g_yz_xx_0_xz = buffer_ddsd[146];

    auto g_yz_xx_0_yy = buffer_ddsd[147];

    auto g_yz_xx_0_yz = buffer_ddsd[148];

    auto g_yz_xx_0_zz = buffer_ddsd[149];

    auto g_yz_xy_0_xx = buffer_ddsd[150];

    auto g_yz_xy_0_xy = buffer_ddsd[151];

    auto g_yz_xy_0_xz = buffer_ddsd[152];

    auto g_yz_xy_0_yy = buffer_ddsd[153];

    auto g_yz_xy_0_yz = buffer_ddsd[154];

    auto g_yz_xy_0_zz = buffer_ddsd[155];

    auto g_yz_xz_0_xx = buffer_ddsd[156];

    auto g_yz_xz_0_xy = buffer_ddsd[157];

    auto g_yz_xz_0_xz = buffer_ddsd[158];

    auto g_yz_xz_0_yy = buffer_ddsd[159];

    auto g_yz_xz_0_yz = buffer_ddsd[160];

    auto g_yz_xz_0_zz = buffer_ddsd[161];

    auto g_yz_yy_0_xx = buffer_ddsd[162];

    auto g_yz_yy_0_xy = buffer_ddsd[163];

    auto g_yz_yy_0_xz = buffer_ddsd[164];

    auto g_yz_yy_0_yy = buffer_ddsd[165];

    auto g_yz_yy_0_yz = buffer_ddsd[166];

    auto g_yz_yy_0_zz = buffer_ddsd[167];

    auto g_yz_yz_0_xx = buffer_ddsd[168];

    auto g_yz_yz_0_xy = buffer_ddsd[169];

    auto g_yz_yz_0_xz = buffer_ddsd[170];

    auto g_yz_yz_0_yy = buffer_ddsd[171];

    auto g_yz_yz_0_yz = buffer_ddsd[172];

    auto g_yz_yz_0_zz = buffer_ddsd[173];

    auto g_yz_zz_0_xx = buffer_ddsd[174];

    auto g_yz_zz_0_xy = buffer_ddsd[175];

    auto g_yz_zz_0_xz = buffer_ddsd[176];

    auto g_yz_zz_0_yy = buffer_ddsd[177];

    auto g_yz_zz_0_yz = buffer_ddsd[178];

    auto g_yz_zz_0_zz = buffer_ddsd[179];

    auto g_zz_xx_0_xx = buffer_ddsd[180];

    auto g_zz_xx_0_xy = buffer_ddsd[181];

    auto g_zz_xx_0_xz = buffer_ddsd[182];

    auto g_zz_xx_0_yy = buffer_ddsd[183];

    auto g_zz_xx_0_yz = buffer_ddsd[184];

    auto g_zz_xx_0_zz = buffer_ddsd[185];

    auto g_zz_xy_0_xx = buffer_ddsd[186];

    auto g_zz_xy_0_xy = buffer_ddsd[187];

    auto g_zz_xy_0_xz = buffer_ddsd[188];

    auto g_zz_xy_0_yy = buffer_ddsd[189];

    auto g_zz_xy_0_yz = buffer_ddsd[190];

    auto g_zz_xy_0_zz = buffer_ddsd[191];

    auto g_zz_xz_0_xx = buffer_ddsd[192];

    auto g_zz_xz_0_xy = buffer_ddsd[193];

    auto g_zz_xz_0_xz = buffer_ddsd[194];

    auto g_zz_xz_0_yy = buffer_ddsd[195];

    auto g_zz_xz_0_yz = buffer_ddsd[196];

    auto g_zz_xz_0_zz = buffer_ddsd[197];

    auto g_zz_yy_0_xx = buffer_ddsd[198];

    auto g_zz_yy_0_xy = buffer_ddsd[199];

    auto g_zz_yy_0_xz = buffer_ddsd[200];

    auto g_zz_yy_0_yy = buffer_ddsd[201];

    auto g_zz_yy_0_yz = buffer_ddsd[202];

    auto g_zz_yy_0_zz = buffer_ddsd[203];

    auto g_zz_yz_0_xx = buffer_ddsd[204];

    auto g_zz_yz_0_xy = buffer_ddsd[205];

    auto g_zz_yz_0_xz = buffer_ddsd[206];

    auto g_zz_yz_0_yy = buffer_ddsd[207];

    auto g_zz_yz_0_yz = buffer_ddsd[208];

    auto g_zz_yz_0_zz = buffer_ddsd[209];

    auto g_zz_zz_0_xx = buffer_ddsd[210];

    auto g_zz_zz_0_xy = buffer_ddsd[211];

    auto g_zz_zz_0_xz = buffer_ddsd[212];

    auto g_zz_zz_0_yy = buffer_ddsd[213];

    auto g_zz_zz_0_yz = buffer_ddsd[214];

    auto g_zz_zz_0_zz = buffer_ddsd[215];

    /// Set up components of auxilary buffer : buffer_gdsd

    auto g_xxxx_xx_0_xx = buffer_gdsd[0];

    auto g_xxxx_xx_0_xy = buffer_gdsd[1];

    auto g_xxxx_xx_0_xz = buffer_gdsd[2];

    auto g_xxxx_xx_0_yy = buffer_gdsd[3];

    auto g_xxxx_xx_0_yz = buffer_gdsd[4];

    auto g_xxxx_xx_0_zz = buffer_gdsd[5];

    auto g_xxxx_xy_0_xx = buffer_gdsd[6];

    auto g_xxxx_xy_0_xy = buffer_gdsd[7];

    auto g_xxxx_xy_0_xz = buffer_gdsd[8];

    auto g_xxxx_xy_0_yy = buffer_gdsd[9];

    auto g_xxxx_xy_0_yz = buffer_gdsd[10];

    auto g_xxxx_xy_0_zz = buffer_gdsd[11];

    auto g_xxxx_xz_0_xx = buffer_gdsd[12];

    auto g_xxxx_xz_0_xy = buffer_gdsd[13];

    auto g_xxxx_xz_0_xz = buffer_gdsd[14];

    auto g_xxxx_xz_0_yy = buffer_gdsd[15];

    auto g_xxxx_xz_0_yz = buffer_gdsd[16];

    auto g_xxxx_xz_0_zz = buffer_gdsd[17];

    auto g_xxxx_yy_0_xx = buffer_gdsd[18];

    auto g_xxxx_yy_0_xy = buffer_gdsd[19];

    auto g_xxxx_yy_0_xz = buffer_gdsd[20];

    auto g_xxxx_yy_0_yy = buffer_gdsd[21];

    auto g_xxxx_yy_0_yz = buffer_gdsd[22];

    auto g_xxxx_yy_0_zz = buffer_gdsd[23];

    auto g_xxxx_yz_0_xx = buffer_gdsd[24];

    auto g_xxxx_yz_0_xy = buffer_gdsd[25];

    auto g_xxxx_yz_0_xz = buffer_gdsd[26];

    auto g_xxxx_yz_0_yy = buffer_gdsd[27];

    auto g_xxxx_yz_0_yz = buffer_gdsd[28];

    auto g_xxxx_yz_0_zz = buffer_gdsd[29];

    auto g_xxxx_zz_0_xx = buffer_gdsd[30];

    auto g_xxxx_zz_0_xy = buffer_gdsd[31];

    auto g_xxxx_zz_0_xz = buffer_gdsd[32];

    auto g_xxxx_zz_0_yy = buffer_gdsd[33];

    auto g_xxxx_zz_0_yz = buffer_gdsd[34];

    auto g_xxxx_zz_0_zz = buffer_gdsd[35];

    auto g_xxxy_xx_0_xx = buffer_gdsd[36];

    auto g_xxxy_xx_0_xy = buffer_gdsd[37];

    auto g_xxxy_xx_0_xz = buffer_gdsd[38];

    auto g_xxxy_xx_0_yy = buffer_gdsd[39];

    auto g_xxxy_xx_0_yz = buffer_gdsd[40];

    auto g_xxxy_xx_0_zz = buffer_gdsd[41];

    auto g_xxxy_xy_0_xx = buffer_gdsd[42];

    auto g_xxxy_xy_0_xy = buffer_gdsd[43];

    auto g_xxxy_xy_0_xz = buffer_gdsd[44];

    auto g_xxxy_xy_0_yy = buffer_gdsd[45];

    auto g_xxxy_xy_0_yz = buffer_gdsd[46];

    auto g_xxxy_xy_0_zz = buffer_gdsd[47];

    auto g_xxxy_xz_0_xx = buffer_gdsd[48];

    auto g_xxxy_xz_0_xy = buffer_gdsd[49];

    auto g_xxxy_xz_0_xz = buffer_gdsd[50];

    auto g_xxxy_xz_0_yy = buffer_gdsd[51];

    auto g_xxxy_xz_0_yz = buffer_gdsd[52];

    auto g_xxxy_xz_0_zz = buffer_gdsd[53];

    auto g_xxxy_yy_0_xx = buffer_gdsd[54];

    auto g_xxxy_yy_0_xy = buffer_gdsd[55];

    auto g_xxxy_yy_0_xz = buffer_gdsd[56];

    auto g_xxxy_yy_0_yy = buffer_gdsd[57];

    auto g_xxxy_yy_0_yz = buffer_gdsd[58];

    auto g_xxxy_yy_0_zz = buffer_gdsd[59];

    auto g_xxxy_yz_0_xx = buffer_gdsd[60];

    auto g_xxxy_yz_0_xy = buffer_gdsd[61];

    auto g_xxxy_yz_0_xz = buffer_gdsd[62];

    auto g_xxxy_yz_0_yy = buffer_gdsd[63];

    auto g_xxxy_yz_0_yz = buffer_gdsd[64];

    auto g_xxxy_yz_0_zz = buffer_gdsd[65];

    auto g_xxxy_zz_0_xx = buffer_gdsd[66];

    auto g_xxxy_zz_0_xy = buffer_gdsd[67];

    auto g_xxxy_zz_0_xz = buffer_gdsd[68];

    auto g_xxxy_zz_0_yy = buffer_gdsd[69];

    auto g_xxxy_zz_0_yz = buffer_gdsd[70];

    auto g_xxxy_zz_0_zz = buffer_gdsd[71];

    auto g_xxxz_xx_0_xx = buffer_gdsd[72];

    auto g_xxxz_xx_0_xy = buffer_gdsd[73];

    auto g_xxxz_xx_0_xz = buffer_gdsd[74];

    auto g_xxxz_xx_0_yy = buffer_gdsd[75];

    auto g_xxxz_xx_0_yz = buffer_gdsd[76];

    auto g_xxxz_xx_0_zz = buffer_gdsd[77];

    auto g_xxxz_xy_0_xx = buffer_gdsd[78];

    auto g_xxxz_xy_0_xy = buffer_gdsd[79];

    auto g_xxxz_xy_0_xz = buffer_gdsd[80];

    auto g_xxxz_xy_0_yy = buffer_gdsd[81];

    auto g_xxxz_xy_0_yz = buffer_gdsd[82];

    auto g_xxxz_xy_0_zz = buffer_gdsd[83];

    auto g_xxxz_xz_0_xx = buffer_gdsd[84];

    auto g_xxxz_xz_0_xy = buffer_gdsd[85];

    auto g_xxxz_xz_0_xz = buffer_gdsd[86];

    auto g_xxxz_xz_0_yy = buffer_gdsd[87];

    auto g_xxxz_xz_0_yz = buffer_gdsd[88];

    auto g_xxxz_xz_0_zz = buffer_gdsd[89];

    auto g_xxxz_yy_0_xx = buffer_gdsd[90];

    auto g_xxxz_yy_0_xy = buffer_gdsd[91];

    auto g_xxxz_yy_0_xz = buffer_gdsd[92];

    auto g_xxxz_yy_0_yy = buffer_gdsd[93];

    auto g_xxxz_yy_0_yz = buffer_gdsd[94];

    auto g_xxxz_yy_0_zz = buffer_gdsd[95];

    auto g_xxxz_yz_0_xx = buffer_gdsd[96];

    auto g_xxxz_yz_0_xy = buffer_gdsd[97];

    auto g_xxxz_yz_0_xz = buffer_gdsd[98];

    auto g_xxxz_yz_0_yy = buffer_gdsd[99];

    auto g_xxxz_yz_0_yz = buffer_gdsd[100];

    auto g_xxxz_yz_0_zz = buffer_gdsd[101];

    auto g_xxxz_zz_0_xx = buffer_gdsd[102];

    auto g_xxxz_zz_0_xy = buffer_gdsd[103];

    auto g_xxxz_zz_0_xz = buffer_gdsd[104];

    auto g_xxxz_zz_0_yy = buffer_gdsd[105];

    auto g_xxxz_zz_0_yz = buffer_gdsd[106];

    auto g_xxxz_zz_0_zz = buffer_gdsd[107];

    auto g_xxyy_xx_0_xx = buffer_gdsd[108];

    auto g_xxyy_xx_0_xy = buffer_gdsd[109];

    auto g_xxyy_xx_0_xz = buffer_gdsd[110];

    auto g_xxyy_xx_0_yy = buffer_gdsd[111];

    auto g_xxyy_xx_0_yz = buffer_gdsd[112];

    auto g_xxyy_xx_0_zz = buffer_gdsd[113];

    auto g_xxyy_xy_0_xx = buffer_gdsd[114];

    auto g_xxyy_xy_0_xy = buffer_gdsd[115];

    auto g_xxyy_xy_0_xz = buffer_gdsd[116];

    auto g_xxyy_xy_0_yy = buffer_gdsd[117];

    auto g_xxyy_xy_0_yz = buffer_gdsd[118];

    auto g_xxyy_xy_0_zz = buffer_gdsd[119];

    auto g_xxyy_xz_0_xx = buffer_gdsd[120];

    auto g_xxyy_xz_0_xy = buffer_gdsd[121];

    auto g_xxyy_xz_0_xz = buffer_gdsd[122];

    auto g_xxyy_xz_0_yy = buffer_gdsd[123];

    auto g_xxyy_xz_0_yz = buffer_gdsd[124];

    auto g_xxyy_xz_0_zz = buffer_gdsd[125];

    auto g_xxyy_yy_0_xx = buffer_gdsd[126];

    auto g_xxyy_yy_0_xy = buffer_gdsd[127];

    auto g_xxyy_yy_0_xz = buffer_gdsd[128];

    auto g_xxyy_yy_0_yy = buffer_gdsd[129];

    auto g_xxyy_yy_0_yz = buffer_gdsd[130];

    auto g_xxyy_yy_0_zz = buffer_gdsd[131];

    auto g_xxyy_yz_0_xx = buffer_gdsd[132];

    auto g_xxyy_yz_0_xy = buffer_gdsd[133];

    auto g_xxyy_yz_0_xz = buffer_gdsd[134];

    auto g_xxyy_yz_0_yy = buffer_gdsd[135];

    auto g_xxyy_yz_0_yz = buffer_gdsd[136];

    auto g_xxyy_yz_0_zz = buffer_gdsd[137];

    auto g_xxyy_zz_0_xx = buffer_gdsd[138];

    auto g_xxyy_zz_0_xy = buffer_gdsd[139];

    auto g_xxyy_zz_0_xz = buffer_gdsd[140];

    auto g_xxyy_zz_0_yy = buffer_gdsd[141];

    auto g_xxyy_zz_0_yz = buffer_gdsd[142];

    auto g_xxyy_zz_0_zz = buffer_gdsd[143];

    auto g_xxyz_xx_0_xx = buffer_gdsd[144];

    auto g_xxyz_xx_0_xy = buffer_gdsd[145];

    auto g_xxyz_xx_0_xz = buffer_gdsd[146];

    auto g_xxyz_xx_0_yy = buffer_gdsd[147];

    auto g_xxyz_xx_0_yz = buffer_gdsd[148];

    auto g_xxyz_xx_0_zz = buffer_gdsd[149];

    auto g_xxyz_xy_0_xx = buffer_gdsd[150];

    auto g_xxyz_xy_0_xy = buffer_gdsd[151];

    auto g_xxyz_xy_0_xz = buffer_gdsd[152];

    auto g_xxyz_xy_0_yy = buffer_gdsd[153];

    auto g_xxyz_xy_0_yz = buffer_gdsd[154];

    auto g_xxyz_xy_0_zz = buffer_gdsd[155];

    auto g_xxyz_xz_0_xx = buffer_gdsd[156];

    auto g_xxyz_xz_0_xy = buffer_gdsd[157];

    auto g_xxyz_xz_0_xz = buffer_gdsd[158];

    auto g_xxyz_xz_0_yy = buffer_gdsd[159];

    auto g_xxyz_xz_0_yz = buffer_gdsd[160];

    auto g_xxyz_xz_0_zz = buffer_gdsd[161];

    auto g_xxyz_yy_0_xx = buffer_gdsd[162];

    auto g_xxyz_yy_0_xy = buffer_gdsd[163];

    auto g_xxyz_yy_0_xz = buffer_gdsd[164];

    auto g_xxyz_yy_0_yy = buffer_gdsd[165];

    auto g_xxyz_yy_0_yz = buffer_gdsd[166];

    auto g_xxyz_yy_0_zz = buffer_gdsd[167];

    auto g_xxyz_yz_0_xx = buffer_gdsd[168];

    auto g_xxyz_yz_0_xy = buffer_gdsd[169];

    auto g_xxyz_yz_0_xz = buffer_gdsd[170];

    auto g_xxyz_yz_0_yy = buffer_gdsd[171];

    auto g_xxyz_yz_0_yz = buffer_gdsd[172];

    auto g_xxyz_yz_0_zz = buffer_gdsd[173];

    auto g_xxyz_zz_0_xx = buffer_gdsd[174];

    auto g_xxyz_zz_0_xy = buffer_gdsd[175];

    auto g_xxyz_zz_0_xz = buffer_gdsd[176];

    auto g_xxyz_zz_0_yy = buffer_gdsd[177];

    auto g_xxyz_zz_0_yz = buffer_gdsd[178];

    auto g_xxyz_zz_0_zz = buffer_gdsd[179];

    auto g_xxzz_xx_0_xx = buffer_gdsd[180];

    auto g_xxzz_xx_0_xy = buffer_gdsd[181];

    auto g_xxzz_xx_0_xz = buffer_gdsd[182];

    auto g_xxzz_xx_0_yy = buffer_gdsd[183];

    auto g_xxzz_xx_0_yz = buffer_gdsd[184];

    auto g_xxzz_xx_0_zz = buffer_gdsd[185];

    auto g_xxzz_xy_0_xx = buffer_gdsd[186];

    auto g_xxzz_xy_0_xy = buffer_gdsd[187];

    auto g_xxzz_xy_0_xz = buffer_gdsd[188];

    auto g_xxzz_xy_0_yy = buffer_gdsd[189];

    auto g_xxzz_xy_0_yz = buffer_gdsd[190];

    auto g_xxzz_xy_0_zz = buffer_gdsd[191];

    auto g_xxzz_xz_0_xx = buffer_gdsd[192];

    auto g_xxzz_xz_0_xy = buffer_gdsd[193];

    auto g_xxzz_xz_0_xz = buffer_gdsd[194];

    auto g_xxzz_xz_0_yy = buffer_gdsd[195];

    auto g_xxzz_xz_0_yz = buffer_gdsd[196];

    auto g_xxzz_xz_0_zz = buffer_gdsd[197];

    auto g_xxzz_yy_0_xx = buffer_gdsd[198];

    auto g_xxzz_yy_0_xy = buffer_gdsd[199];

    auto g_xxzz_yy_0_xz = buffer_gdsd[200];

    auto g_xxzz_yy_0_yy = buffer_gdsd[201];

    auto g_xxzz_yy_0_yz = buffer_gdsd[202];

    auto g_xxzz_yy_0_zz = buffer_gdsd[203];

    auto g_xxzz_yz_0_xx = buffer_gdsd[204];

    auto g_xxzz_yz_0_xy = buffer_gdsd[205];

    auto g_xxzz_yz_0_xz = buffer_gdsd[206];

    auto g_xxzz_yz_0_yy = buffer_gdsd[207];

    auto g_xxzz_yz_0_yz = buffer_gdsd[208];

    auto g_xxzz_yz_0_zz = buffer_gdsd[209];

    auto g_xxzz_zz_0_xx = buffer_gdsd[210];

    auto g_xxzz_zz_0_xy = buffer_gdsd[211];

    auto g_xxzz_zz_0_xz = buffer_gdsd[212];

    auto g_xxzz_zz_0_yy = buffer_gdsd[213];

    auto g_xxzz_zz_0_yz = buffer_gdsd[214];

    auto g_xxzz_zz_0_zz = buffer_gdsd[215];

    auto g_xyyy_xx_0_xx = buffer_gdsd[216];

    auto g_xyyy_xx_0_xy = buffer_gdsd[217];

    auto g_xyyy_xx_0_xz = buffer_gdsd[218];

    auto g_xyyy_xx_0_yy = buffer_gdsd[219];

    auto g_xyyy_xx_0_yz = buffer_gdsd[220];

    auto g_xyyy_xx_0_zz = buffer_gdsd[221];

    auto g_xyyy_xy_0_xx = buffer_gdsd[222];

    auto g_xyyy_xy_0_xy = buffer_gdsd[223];

    auto g_xyyy_xy_0_xz = buffer_gdsd[224];

    auto g_xyyy_xy_0_yy = buffer_gdsd[225];

    auto g_xyyy_xy_0_yz = buffer_gdsd[226];

    auto g_xyyy_xy_0_zz = buffer_gdsd[227];

    auto g_xyyy_xz_0_xx = buffer_gdsd[228];

    auto g_xyyy_xz_0_xy = buffer_gdsd[229];

    auto g_xyyy_xz_0_xz = buffer_gdsd[230];

    auto g_xyyy_xz_0_yy = buffer_gdsd[231];

    auto g_xyyy_xz_0_yz = buffer_gdsd[232];

    auto g_xyyy_xz_0_zz = buffer_gdsd[233];

    auto g_xyyy_yy_0_xx = buffer_gdsd[234];

    auto g_xyyy_yy_0_xy = buffer_gdsd[235];

    auto g_xyyy_yy_0_xz = buffer_gdsd[236];

    auto g_xyyy_yy_0_yy = buffer_gdsd[237];

    auto g_xyyy_yy_0_yz = buffer_gdsd[238];

    auto g_xyyy_yy_0_zz = buffer_gdsd[239];

    auto g_xyyy_yz_0_xx = buffer_gdsd[240];

    auto g_xyyy_yz_0_xy = buffer_gdsd[241];

    auto g_xyyy_yz_0_xz = buffer_gdsd[242];

    auto g_xyyy_yz_0_yy = buffer_gdsd[243];

    auto g_xyyy_yz_0_yz = buffer_gdsd[244];

    auto g_xyyy_yz_0_zz = buffer_gdsd[245];

    auto g_xyyy_zz_0_xx = buffer_gdsd[246];

    auto g_xyyy_zz_0_xy = buffer_gdsd[247];

    auto g_xyyy_zz_0_xz = buffer_gdsd[248];

    auto g_xyyy_zz_0_yy = buffer_gdsd[249];

    auto g_xyyy_zz_0_yz = buffer_gdsd[250];

    auto g_xyyy_zz_0_zz = buffer_gdsd[251];

    auto g_xyyz_xx_0_xx = buffer_gdsd[252];

    auto g_xyyz_xx_0_xy = buffer_gdsd[253];

    auto g_xyyz_xx_0_xz = buffer_gdsd[254];

    auto g_xyyz_xx_0_yy = buffer_gdsd[255];

    auto g_xyyz_xx_0_yz = buffer_gdsd[256];

    auto g_xyyz_xx_0_zz = buffer_gdsd[257];

    auto g_xyyz_xy_0_xx = buffer_gdsd[258];

    auto g_xyyz_xy_0_xy = buffer_gdsd[259];

    auto g_xyyz_xy_0_xz = buffer_gdsd[260];

    auto g_xyyz_xy_0_yy = buffer_gdsd[261];

    auto g_xyyz_xy_0_yz = buffer_gdsd[262];

    auto g_xyyz_xy_0_zz = buffer_gdsd[263];

    auto g_xyyz_xz_0_xx = buffer_gdsd[264];

    auto g_xyyz_xz_0_xy = buffer_gdsd[265];

    auto g_xyyz_xz_0_xz = buffer_gdsd[266];

    auto g_xyyz_xz_0_yy = buffer_gdsd[267];

    auto g_xyyz_xz_0_yz = buffer_gdsd[268];

    auto g_xyyz_xz_0_zz = buffer_gdsd[269];

    auto g_xyyz_yy_0_xx = buffer_gdsd[270];

    auto g_xyyz_yy_0_xy = buffer_gdsd[271];

    auto g_xyyz_yy_0_xz = buffer_gdsd[272];

    auto g_xyyz_yy_0_yy = buffer_gdsd[273];

    auto g_xyyz_yy_0_yz = buffer_gdsd[274];

    auto g_xyyz_yy_0_zz = buffer_gdsd[275];

    auto g_xyyz_yz_0_xx = buffer_gdsd[276];

    auto g_xyyz_yz_0_xy = buffer_gdsd[277];

    auto g_xyyz_yz_0_xz = buffer_gdsd[278];

    auto g_xyyz_yz_0_yy = buffer_gdsd[279];

    auto g_xyyz_yz_0_yz = buffer_gdsd[280];

    auto g_xyyz_yz_0_zz = buffer_gdsd[281];

    auto g_xyyz_zz_0_xx = buffer_gdsd[282];

    auto g_xyyz_zz_0_xy = buffer_gdsd[283];

    auto g_xyyz_zz_0_xz = buffer_gdsd[284];

    auto g_xyyz_zz_0_yy = buffer_gdsd[285];

    auto g_xyyz_zz_0_yz = buffer_gdsd[286];

    auto g_xyyz_zz_0_zz = buffer_gdsd[287];

    auto g_xyzz_xx_0_xx = buffer_gdsd[288];

    auto g_xyzz_xx_0_xy = buffer_gdsd[289];

    auto g_xyzz_xx_0_xz = buffer_gdsd[290];

    auto g_xyzz_xx_0_yy = buffer_gdsd[291];

    auto g_xyzz_xx_0_yz = buffer_gdsd[292];

    auto g_xyzz_xx_0_zz = buffer_gdsd[293];

    auto g_xyzz_xy_0_xx = buffer_gdsd[294];

    auto g_xyzz_xy_0_xy = buffer_gdsd[295];

    auto g_xyzz_xy_0_xz = buffer_gdsd[296];

    auto g_xyzz_xy_0_yy = buffer_gdsd[297];

    auto g_xyzz_xy_0_yz = buffer_gdsd[298];

    auto g_xyzz_xy_0_zz = buffer_gdsd[299];

    auto g_xyzz_xz_0_xx = buffer_gdsd[300];

    auto g_xyzz_xz_0_xy = buffer_gdsd[301];

    auto g_xyzz_xz_0_xz = buffer_gdsd[302];

    auto g_xyzz_xz_0_yy = buffer_gdsd[303];

    auto g_xyzz_xz_0_yz = buffer_gdsd[304];

    auto g_xyzz_xz_0_zz = buffer_gdsd[305];

    auto g_xyzz_yy_0_xx = buffer_gdsd[306];

    auto g_xyzz_yy_0_xy = buffer_gdsd[307];

    auto g_xyzz_yy_0_xz = buffer_gdsd[308];

    auto g_xyzz_yy_0_yy = buffer_gdsd[309];

    auto g_xyzz_yy_0_yz = buffer_gdsd[310];

    auto g_xyzz_yy_0_zz = buffer_gdsd[311];

    auto g_xyzz_yz_0_xx = buffer_gdsd[312];

    auto g_xyzz_yz_0_xy = buffer_gdsd[313];

    auto g_xyzz_yz_0_xz = buffer_gdsd[314];

    auto g_xyzz_yz_0_yy = buffer_gdsd[315];

    auto g_xyzz_yz_0_yz = buffer_gdsd[316];

    auto g_xyzz_yz_0_zz = buffer_gdsd[317];

    auto g_xyzz_zz_0_xx = buffer_gdsd[318];

    auto g_xyzz_zz_0_xy = buffer_gdsd[319];

    auto g_xyzz_zz_0_xz = buffer_gdsd[320];

    auto g_xyzz_zz_0_yy = buffer_gdsd[321];

    auto g_xyzz_zz_0_yz = buffer_gdsd[322];

    auto g_xyzz_zz_0_zz = buffer_gdsd[323];

    auto g_xzzz_xx_0_xx = buffer_gdsd[324];

    auto g_xzzz_xx_0_xy = buffer_gdsd[325];

    auto g_xzzz_xx_0_xz = buffer_gdsd[326];

    auto g_xzzz_xx_0_yy = buffer_gdsd[327];

    auto g_xzzz_xx_0_yz = buffer_gdsd[328];

    auto g_xzzz_xx_0_zz = buffer_gdsd[329];

    auto g_xzzz_xy_0_xx = buffer_gdsd[330];

    auto g_xzzz_xy_0_xy = buffer_gdsd[331];

    auto g_xzzz_xy_0_xz = buffer_gdsd[332];

    auto g_xzzz_xy_0_yy = buffer_gdsd[333];

    auto g_xzzz_xy_0_yz = buffer_gdsd[334];

    auto g_xzzz_xy_0_zz = buffer_gdsd[335];

    auto g_xzzz_xz_0_xx = buffer_gdsd[336];

    auto g_xzzz_xz_0_xy = buffer_gdsd[337];

    auto g_xzzz_xz_0_xz = buffer_gdsd[338];

    auto g_xzzz_xz_0_yy = buffer_gdsd[339];

    auto g_xzzz_xz_0_yz = buffer_gdsd[340];

    auto g_xzzz_xz_0_zz = buffer_gdsd[341];

    auto g_xzzz_yy_0_xx = buffer_gdsd[342];

    auto g_xzzz_yy_0_xy = buffer_gdsd[343];

    auto g_xzzz_yy_0_xz = buffer_gdsd[344];

    auto g_xzzz_yy_0_yy = buffer_gdsd[345];

    auto g_xzzz_yy_0_yz = buffer_gdsd[346];

    auto g_xzzz_yy_0_zz = buffer_gdsd[347];

    auto g_xzzz_yz_0_xx = buffer_gdsd[348];

    auto g_xzzz_yz_0_xy = buffer_gdsd[349];

    auto g_xzzz_yz_0_xz = buffer_gdsd[350];

    auto g_xzzz_yz_0_yy = buffer_gdsd[351];

    auto g_xzzz_yz_0_yz = buffer_gdsd[352];

    auto g_xzzz_yz_0_zz = buffer_gdsd[353];

    auto g_xzzz_zz_0_xx = buffer_gdsd[354];

    auto g_xzzz_zz_0_xy = buffer_gdsd[355];

    auto g_xzzz_zz_0_xz = buffer_gdsd[356];

    auto g_xzzz_zz_0_yy = buffer_gdsd[357];

    auto g_xzzz_zz_0_yz = buffer_gdsd[358];

    auto g_xzzz_zz_0_zz = buffer_gdsd[359];

    auto g_yyyy_xx_0_xx = buffer_gdsd[360];

    auto g_yyyy_xx_0_xy = buffer_gdsd[361];

    auto g_yyyy_xx_0_xz = buffer_gdsd[362];

    auto g_yyyy_xx_0_yy = buffer_gdsd[363];

    auto g_yyyy_xx_0_yz = buffer_gdsd[364];

    auto g_yyyy_xx_0_zz = buffer_gdsd[365];

    auto g_yyyy_xy_0_xx = buffer_gdsd[366];

    auto g_yyyy_xy_0_xy = buffer_gdsd[367];

    auto g_yyyy_xy_0_xz = buffer_gdsd[368];

    auto g_yyyy_xy_0_yy = buffer_gdsd[369];

    auto g_yyyy_xy_0_yz = buffer_gdsd[370];

    auto g_yyyy_xy_0_zz = buffer_gdsd[371];

    auto g_yyyy_xz_0_xx = buffer_gdsd[372];

    auto g_yyyy_xz_0_xy = buffer_gdsd[373];

    auto g_yyyy_xz_0_xz = buffer_gdsd[374];

    auto g_yyyy_xz_0_yy = buffer_gdsd[375];

    auto g_yyyy_xz_0_yz = buffer_gdsd[376];

    auto g_yyyy_xz_0_zz = buffer_gdsd[377];

    auto g_yyyy_yy_0_xx = buffer_gdsd[378];

    auto g_yyyy_yy_0_xy = buffer_gdsd[379];

    auto g_yyyy_yy_0_xz = buffer_gdsd[380];

    auto g_yyyy_yy_0_yy = buffer_gdsd[381];

    auto g_yyyy_yy_0_yz = buffer_gdsd[382];

    auto g_yyyy_yy_0_zz = buffer_gdsd[383];

    auto g_yyyy_yz_0_xx = buffer_gdsd[384];

    auto g_yyyy_yz_0_xy = buffer_gdsd[385];

    auto g_yyyy_yz_0_xz = buffer_gdsd[386];

    auto g_yyyy_yz_0_yy = buffer_gdsd[387];

    auto g_yyyy_yz_0_yz = buffer_gdsd[388];

    auto g_yyyy_yz_0_zz = buffer_gdsd[389];

    auto g_yyyy_zz_0_xx = buffer_gdsd[390];

    auto g_yyyy_zz_0_xy = buffer_gdsd[391];

    auto g_yyyy_zz_0_xz = buffer_gdsd[392];

    auto g_yyyy_zz_0_yy = buffer_gdsd[393];

    auto g_yyyy_zz_0_yz = buffer_gdsd[394];

    auto g_yyyy_zz_0_zz = buffer_gdsd[395];

    auto g_yyyz_xx_0_xx = buffer_gdsd[396];

    auto g_yyyz_xx_0_xy = buffer_gdsd[397];

    auto g_yyyz_xx_0_xz = buffer_gdsd[398];

    auto g_yyyz_xx_0_yy = buffer_gdsd[399];

    auto g_yyyz_xx_0_yz = buffer_gdsd[400];

    auto g_yyyz_xx_0_zz = buffer_gdsd[401];

    auto g_yyyz_xy_0_xx = buffer_gdsd[402];

    auto g_yyyz_xy_0_xy = buffer_gdsd[403];

    auto g_yyyz_xy_0_xz = buffer_gdsd[404];

    auto g_yyyz_xy_0_yy = buffer_gdsd[405];

    auto g_yyyz_xy_0_yz = buffer_gdsd[406];

    auto g_yyyz_xy_0_zz = buffer_gdsd[407];

    auto g_yyyz_xz_0_xx = buffer_gdsd[408];

    auto g_yyyz_xz_0_xy = buffer_gdsd[409];

    auto g_yyyz_xz_0_xz = buffer_gdsd[410];

    auto g_yyyz_xz_0_yy = buffer_gdsd[411];

    auto g_yyyz_xz_0_yz = buffer_gdsd[412];

    auto g_yyyz_xz_0_zz = buffer_gdsd[413];

    auto g_yyyz_yy_0_xx = buffer_gdsd[414];

    auto g_yyyz_yy_0_xy = buffer_gdsd[415];

    auto g_yyyz_yy_0_xz = buffer_gdsd[416];

    auto g_yyyz_yy_0_yy = buffer_gdsd[417];

    auto g_yyyz_yy_0_yz = buffer_gdsd[418];

    auto g_yyyz_yy_0_zz = buffer_gdsd[419];

    auto g_yyyz_yz_0_xx = buffer_gdsd[420];

    auto g_yyyz_yz_0_xy = buffer_gdsd[421];

    auto g_yyyz_yz_0_xz = buffer_gdsd[422];

    auto g_yyyz_yz_0_yy = buffer_gdsd[423];

    auto g_yyyz_yz_0_yz = buffer_gdsd[424];

    auto g_yyyz_yz_0_zz = buffer_gdsd[425];

    auto g_yyyz_zz_0_xx = buffer_gdsd[426];

    auto g_yyyz_zz_0_xy = buffer_gdsd[427];

    auto g_yyyz_zz_0_xz = buffer_gdsd[428];

    auto g_yyyz_zz_0_yy = buffer_gdsd[429];

    auto g_yyyz_zz_0_yz = buffer_gdsd[430];

    auto g_yyyz_zz_0_zz = buffer_gdsd[431];

    auto g_yyzz_xx_0_xx = buffer_gdsd[432];

    auto g_yyzz_xx_0_xy = buffer_gdsd[433];

    auto g_yyzz_xx_0_xz = buffer_gdsd[434];

    auto g_yyzz_xx_0_yy = buffer_gdsd[435];

    auto g_yyzz_xx_0_yz = buffer_gdsd[436];

    auto g_yyzz_xx_0_zz = buffer_gdsd[437];

    auto g_yyzz_xy_0_xx = buffer_gdsd[438];

    auto g_yyzz_xy_0_xy = buffer_gdsd[439];

    auto g_yyzz_xy_0_xz = buffer_gdsd[440];

    auto g_yyzz_xy_0_yy = buffer_gdsd[441];

    auto g_yyzz_xy_0_yz = buffer_gdsd[442];

    auto g_yyzz_xy_0_zz = buffer_gdsd[443];

    auto g_yyzz_xz_0_xx = buffer_gdsd[444];

    auto g_yyzz_xz_0_xy = buffer_gdsd[445];

    auto g_yyzz_xz_0_xz = buffer_gdsd[446];

    auto g_yyzz_xz_0_yy = buffer_gdsd[447];

    auto g_yyzz_xz_0_yz = buffer_gdsd[448];

    auto g_yyzz_xz_0_zz = buffer_gdsd[449];

    auto g_yyzz_yy_0_xx = buffer_gdsd[450];

    auto g_yyzz_yy_0_xy = buffer_gdsd[451];

    auto g_yyzz_yy_0_xz = buffer_gdsd[452];

    auto g_yyzz_yy_0_yy = buffer_gdsd[453];

    auto g_yyzz_yy_0_yz = buffer_gdsd[454];

    auto g_yyzz_yy_0_zz = buffer_gdsd[455];

    auto g_yyzz_yz_0_xx = buffer_gdsd[456];

    auto g_yyzz_yz_0_xy = buffer_gdsd[457];

    auto g_yyzz_yz_0_xz = buffer_gdsd[458];

    auto g_yyzz_yz_0_yy = buffer_gdsd[459];

    auto g_yyzz_yz_0_yz = buffer_gdsd[460];

    auto g_yyzz_yz_0_zz = buffer_gdsd[461];

    auto g_yyzz_zz_0_xx = buffer_gdsd[462];

    auto g_yyzz_zz_0_xy = buffer_gdsd[463];

    auto g_yyzz_zz_0_xz = buffer_gdsd[464];

    auto g_yyzz_zz_0_yy = buffer_gdsd[465];

    auto g_yyzz_zz_0_yz = buffer_gdsd[466];

    auto g_yyzz_zz_0_zz = buffer_gdsd[467];

    auto g_yzzz_xx_0_xx = buffer_gdsd[468];

    auto g_yzzz_xx_0_xy = buffer_gdsd[469];

    auto g_yzzz_xx_0_xz = buffer_gdsd[470];

    auto g_yzzz_xx_0_yy = buffer_gdsd[471];

    auto g_yzzz_xx_0_yz = buffer_gdsd[472];

    auto g_yzzz_xx_0_zz = buffer_gdsd[473];

    auto g_yzzz_xy_0_xx = buffer_gdsd[474];

    auto g_yzzz_xy_0_xy = buffer_gdsd[475];

    auto g_yzzz_xy_0_xz = buffer_gdsd[476];

    auto g_yzzz_xy_0_yy = buffer_gdsd[477];

    auto g_yzzz_xy_0_yz = buffer_gdsd[478];

    auto g_yzzz_xy_0_zz = buffer_gdsd[479];

    auto g_yzzz_xz_0_xx = buffer_gdsd[480];

    auto g_yzzz_xz_0_xy = buffer_gdsd[481];

    auto g_yzzz_xz_0_xz = buffer_gdsd[482];

    auto g_yzzz_xz_0_yy = buffer_gdsd[483];

    auto g_yzzz_xz_0_yz = buffer_gdsd[484];

    auto g_yzzz_xz_0_zz = buffer_gdsd[485];

    auto g_yzzz_yy_0_xx = buffer_gdsd[486];

    auto g_yzzz_yy_0_xy = buffer_gdsd[487];

    auto g_yzzz_yy_0_xz = buffer_gdsd[488];

    auto g_yzzz_yy_0_yy = buffer_gdsd[489];

    auto g_yzzz_yy_0_yz = buffer_gdsd[490];

    auto g_yzzz_yy_0_zz = buffer_gdsd[491];

    auto g_yzzz_yz_0_xx = buffer_gdsd[492];

    auto g_yzzz_yz_0_xy = buffer_gdsd[493];

    auto g_yzzz_yz_0_xz = buffer_gdsd[494];

    auto g_yzzz_yz_0_yy = buffer_gdsd[495];

    auto g_yzzz_yz_0_yz = buffer_gdsd[496];

    auto g_yzzz_yz_0_zz = buffer_gdsd[497];

    auto g_yzzz_zz_0_xx = buffer_gdsd[498];

    auto g_yzzz_zz_0_xy = buffer_gdsd[499];

    auto g_yzzz_zz_0_xz = buffer_gdsd[500];

    auto g_yzzz_zz_0_yy = buffer_gdsd[501];

    auto g_yzzz_zz_0_yz = buffer_gdsd[502];

    auto g_yzzz_zz_0_zz = buffer_gdsd[503];

    auto g_zzzz_xx_0_xx = buffer_gdsd[504];

    auto g_zzzz_xx_0_xy = buffer_gdsd[505];

    auto g_zzzz_xx_0_xz = buffer_gdsd[506];

    auto g_zzzz_xx_0_yy = buffer_gdsd[507];

    auto g_zzzz_xx_0_yz = buffer_gdsd[508];

    auto g_zzzz_xx_0_zz = buffer_gdsd[509];

    auto g_zzzz_xy_0_xx = buffer_gdsd[510];

    auto g_zzzz_xy_0_xy = buffer_gdsd[511];

    auto g_zzzz_xy_0_xz = buffer_gdsd[512];

    auto g_zzzz_xy_0_yy = buffer_gdsd[513];

    auto g_zzzz_xy_0_yz = buffer_gdsd[514];

    auto g_zzzz_xy_0_zz = buffer_gdsd[515];

    auto g_zzzz_xz_0_xx = buffer_gdsd[516];

    auto g_zzzz_xz_0_xy = buffer_gdsd[517];

    auto g_zzzz_xz_0_xz = buffer_gdsd[518];

    auto g_zzzz_xz_0_yy = buffer_gdsd[519];

    auto g_zzzz_xz_0_yz = buffer_gdsd[520];

    auto g_zzzz_xz_0_zz = buffer_gdsd[521];

    auto g_zzzz_yy_0_xx = buffer_gdsd[522];

    auto g_zzzz_yy_0_xy = buffer_gdsd[523];

    auto g_zzzz_yy_0_xz = buffer_gdsd[524];

    auto g_zzzz_yy_0_yy = buffer_gdsd[525];

    auto g_zzzz_yy_0_yz = buffer_gdsd[526];

    auto g_zzzz_yy_0_zz = buffer_gdsd[527];

    auto g_zzzz_yz_0_xx = buffer_gdsd[528];

    auto g_zzzz_yz_0_xy = buffer_gdsd[529];

    auto g_zzzz_yz_0_xz = buffer_gdsd[530];

    auto g_zzzz_yz_0_yy = buffer_gdsd[531];

    auto g_zzzz_yz_0_yz = buffer_gdsd[532];

    auto g_zzzz_yz_0_zz = buffer_gdsd[533];

    auto g_zzzz_zz_0_xx = buffer_gdsd[534];

    auto g_zzzz_zz_0_xy = buffer_gdsd[535];

    auto g_zzzz_zz_0_xz = buffer_gdsd[536];

    auto g_zzzz_zz_0_yy = buffer_gdsd[537];

    auto g_zzzz_zz_0_yz = buffer_gdsd[538];

    auto g_zzzz_zz_0_zz = buffer_gdsd[539];

    /// Set up components of integrals buffer : buffer_2000_ddsd

    auto g_xx_0_0_0_xx_xx_0_xx = buffer_2000_ddsd[0];

    auto g_xx_0_0_0_xx_xx_0_xy = buffer_2000_ddsd[1];

    auto g_xx_0_0_0_xx_xx_0_xz = buffer_2000_ddsd[2];

    auto g_xx_0_0_0_xx_xx_0_yy = buffer_2000_ddsd[3];

    auto g_xx_0_0_0_xx_xx_0_yz = buffer_2000_ddsd[4];

    auto g_xx_0_0_0_xx_xx_0_zz = buffer_2000_ddsd[5];

    auto g_xx_0_0_0_xx_xy_0_xx = buffer_2000_ddsd[6];

    auto g_xx_0_0_0_xx_xy_0_xy = buffer_2000_ddsd[7];

    auto g_xx_0_0_0_xx_xy_0_xz = buffer_2000_ddsd[8];

    auto g_xx_0_0_0_xx_xy_0_yy = buffer_2000_ddsd[9];

    auto g_xx_0_0_0_xx_xy_0_yz = buffer_2000_ddsd[10];

    auto g_xx_0_0_0_xx_xy_0_zz = buffer_2000_ddsd[11];

    auto g_xx_0_0_0_xx_xz_0_xx = buffer_2000_ddsd[12];

    auto g_xx_0_0_0_xx_xz_0_xy = buffer_2000_ddsd[13];

    auto g_xx_0_0_0_xx_xz_0_xz = buffer_2000_ddsd[14];

    auto g_xx_0_0_0_xx_xz_0_yy = buffer_2000_ddsd[15];

    auto g_xx_0_0_0_xx_xz_0_yz = buffer_2000_ddsd[16];

    auto g_xx_0_0_0_xx_xz_0_zz = buffer_2000_ddsd[17];

    auto g_xx_0_0_0_xx_yy_0_xx = buffer_2000_ddsd[18];

    auto g_xx_0_0_0_xx_yy_0_xy = buffer_2000_ddsd[19];

    auto g_xx_0_0_0_xx_yy_0_xz = buffer_2000_ddsd[20];

    auto g_xx_0_0_0_xx_yy_0_yy = buffer_2000_ddsd[21];

    auto g_xx_0_0_0_xx_yy_0_yz = buffer_2000_ddsd[22];

    auto g_xx_0_0_0_xx_yy_0_zz = buffer_2000_ddsd[23];

    auto g_xx_0_0_0_xx_yz_0_xx = buffer_2000_ddsd[24];

    auto g_xx_0_0_0_xx_yz_0_xy = buffer_2000_ddsd[25];

    auto g_xx_0_0_0_xx_yz_0_xz = buffer_2000_ddsd[26];

    auto g_xx_0_0_0_xx_yz_0_yy = buffer_2000_ddsd[27];

    auto g_xx_0_0_0_xx_yz_0_yz = buffer_2000_ddsd[28];

    auto g_xx_0_0_0_xx_yz_0_zz = buffer_2000_ddsd[29];

    auto g_xx_0_0_0_xx_zz_0_xx = buffer_2000_ddsd[30];

    auto g_xx_0_0_0_xx_zz_0_xy = buffer_2000_ddsd[31];

    auto g_xx_0_0_0_xx_zz_0_xz = buffer_2000_ddsd[32];

    auto g_xx_0_0_0_xx_zz_0_yy = buffer_2000_ddsd[33];

    auto g_xx_0_0_0_xx_zz_0_yz = buffer_2000_ddsd[34];

    auto g_xx_0_0_0_xx_zz_0_zz = buffer_2000_ddsd[35];

    auto g_xx_0_0_0_xy_xx_0_xx = buffer_2000_ddsd[36];

    auto g_xx_0_0_0_xy_xx_0_xy = buffer_2000_ddsd[37];

    auto g_xx_0_0_0_xy_xx_0_xz = buffer_2000_ddsd[38];

    auto g_xx_0_0_0_xy_xx_0_yy = buffer_2000_ddsd[39];

    auto g_xx_0_0_0_xy_xx_0_yz = buffer_2000_ddsd[40];

    auto g_xx_0_0_0_xy_xx_0_zz = buffer_2000_ddsd[41];

    auto g_xx_0_0_0_xy_xy_0_xx = buffer_2000_ddsd[42];

    auto g_xx_0_0_0_xy_xy_0_xy = buffer_2000_ddsd[43];

    auto g_xx_0_0_0_xy_xy_0_xz = buffer_2000_ddsd[44];

    auto g_xx_0_0_0_xy_xy_0_yy = buffer_2000_ddsd[45];

    auto g_xx_0_0_0_xy_xy_0_yz = buffer_2000_ddsd[46];

    auto g_xx_0_0_0_xy_xy_0_zz = buffer_2000_ddsd[47];

    auto g_xx_0_0_0_xy_xz_0_xx = buffer_2000_ddsd[48];

    auto g_xx_0_0_0_xy_xz_0_xy = buffer_2000_ddsd[49];

    auto g_xx_0_0_0_xy_xz_0_xz = buffer_2000_ddsd[50];

    auto g_xx_0_0_0_xy_xz_0_yy = buffer_2000_ddsd[51];

    auto g_xx_0_0_0_xy_xz_0_yz = buffer_2000_ddsd[52];

    auto g_xx_0_0_0_xy_xz_0_zz = buffer_2000_ddsd[53];

    auto g_xx_0_0_0_xy_yy_0_xx = buffer_2000_ddsd[54];

    auto g_xx_0_0_0_xy_yy_0_xy = buffer_2000_ddsd[55];

    auto g_xx_0_0_0_xy_yy_0_xz = buffer_2000_ddsd[56];

    auto g_xx_0_0_0_xy_yy_0_yy = buffer_2000_ddsd[57];

    auto g_xx_0_0_0_xy_yy_0_yz = buffer_2000_ddsd[58];

    auto g_xx_0_0_0_xy_yy_0_zz = buffer_2000_ddsd[59];

    auto g_xx_0_0_0_xy_yz_0_xx = buffer_2000_ddsd[60];

    auto g_xx_0_0_0_xy_yz_0_xy = buffer_2000_ddsd[61];

    auto g_xx_0_0_0_xy_yz_0_xz = buffer_2000_ddsd[62];

    auto g_xx_0_0_0_xy_yz_0_yy = buffer_2000_ddsd[63];

    auto g_xx_0_0_0_xy_yz_0_yz = buffer_2000_ddsd[64];

    auto g_xx_0_0_0_xy_yz_0_zz = buffer_2000_ddsd[65];

    auto g_xx_0_0_0_xy_zz_0_xx = buffer_2000_ddsd[66];

    auto g_xx_0_0_0_xy_zz_0_xy = buffer_2000_ddsd[67];

    auto g_xx_0_0_0_xy_zz_0_xz = buffer_2000_ddsd[68];

    auto g_xx_0_0_0_xy_zz_0_yy = buffer_2000_ddsd[69];

    auto g_xx_0_0_0_xy_zz_0_yz = buffer_2000_ddsd[70];

    auto g_xx_0_0_0_xy_zz_0_zz = buffer_2000_ddsd[71];

    auto g_xx_0_0_0_xz_xx_0_xx = buffer_2000_ddsd[72];

    auto g_xx_0_0_0_xz_xx_0_xy = buffer_2000_ddsd[73];

    auto g_xx_0_0_0_xz_xx_0_xz = buffer_2000_ddsd[74];

    auto g_xx_0_0_0_xz_xx_0_yy = buffer_2000_ddsd[75];

    auto g_xx_0_0_0_xz_xx_0_yz = buffer_2000_ddsd[76];

    auto g_xx_0_0_0_xz_xx_0_zz = buffer_2000_ddsd[77];

    auto g_xx_0_0_0_xz_xy_0_xx = buffer_2000_ddsd[78];

    auto g_xx_0_0_0_xz_xy_0_xy = buffer_2000_ddsd[79];

    auto g_xx_0_0_0_xz_xy_0_xz = buffer_2000_ddsd[80];

    auto g_xx_0_0_0_xz_xy_0_yy = buffer_2000_ddsd[81];

    auto g_xx_0_0_0_xz_xy_0_yz = buffer_2000_ddsd[82];

    auto g_xx_0_0_0_xz_xy_0_zz = buffer_2000_ddsd[83];

    auto g_xx_0_0_0_xz_xz_0_xx = buffer_2000_ddsd[84];

    auto g_xx_0_0_0_xz_xz_0_xy = buffer_2000_ddsd[85];

    auto g_xx_0_0_0_xz_xz_0_xz = buffer_2000_ddsd[86];

    auto g_xx_0_0_0_xz_xz_0_yy = buffer_2000_ddsd[87];

    auto g_xx_0_0_0_xz_xz_0_yz = buffer_2000_ddsd[88];

    auto g_xx_0_0_0_xz_xz_0_zz = buffer_2000_ddsd[89];

    auto g_xx_0_0_0_xz_yy_0_xx = buffer_2000_ddsd[90];

    auto g_xx_0_0_0_xz_yy_0_xy = buffer_2000_ddsd[91];

    auto g_xx_0_0_0_xz_yy_0_xz = buffer_2000_ddsd[92];

    auto g_xx_0_0_0_xz_yy_0_yy = buffer_2000_ddsd[93];

    auto g_xx_0_0_0_xz_yy_0_yz = buffer_2000_ddsd[94];

    auto g_xx_0_0_0_xz_yy_0_zz = buffer_2000_ddsd[95];

    auto g_xx_0_0_0_xz_yz_0_xx = buffer_2000_ddsd[96];

    auto g_xx_0_0_0_xz_yz_0_xy = buffer_2000_ddsd[97];

    auto g_xx_0_0_0_xz_yz_0_xz = buffer_2000_ddsd[98];

    auto g_xx_0_0_0_xz_yz_0_yy = buffer_2000_ddsd[99];

    auto g_xx_0_0_0_xz_yz_0_yz = buffer_2000_ddsd[100];

    auto g_xx_0_0_0_xz_yz_0_zz = buffer_2000_ddsd[101];

    auto g_xx_0_0_0_xz_zz_0_xx = buffer_2000_ddsd[102];

    auto g_xx_0_0_0_xz_zz_0_xy = buffer_2000_ddsd[103];

    auto g_xx_0_0_0_xz_zz_0_xz = buffer_2000_ddsd[104];

    auto g_xx_0_0_0_xz_zz_0_yy = buffer_2000_ddsd[105];

    auto g_xx_0_0_0_xz_zz_0_yz = buffer_2000_ddsd[106];

    auto g_xx_0_0_0_xz_zz_0_zz = buffer_2000_ddsd[107];

    auto g_xx_0_0_0_yy_xx_0_xx = buffer_2000_ddsd[108];

    auto g_xx_0_0_0_yy_xx_0_xy = buffer_2000_ddsd[109];

    auto g_xx_0_0_0_yy_xx_0_xz = buffer_2000_ddsd[110];

    auto g_xx_0_0_0_yy_xx_0_yy = buffer_2000_ddsd[111];

    auto g_xx_0_0_0_yy_xx_0_yz = buffer_2000_ddsd[112];

    auto g_xx_0_0_0_yy_xx_0_zz = buffer_2000_ddsd[113];

    auto g_xx_0_0_0_yy_xy_0_xx = buffer_2000_ddsd[114];

    auto g_xx_0_0_0_yy_xy_0_xy = buffer_2000_ddsd[115];

    auto g_xx_0_0_0_yy_xy_0_xz = buffer_2000_ddsd[116];

    auto g_xx_0_0_0_yy_xy_0_yy = buffer_2000_ddsd[117];

    auto g_xx_0_0_0_yy_xy_0_yz = buffer_2000_ddsd[118];

    auto g_xx_0_0_0_yy_xy_0_zz = buffer_2000_ddsd[119];

    auto g_xx_0_0_0_yy_xz_0_xx = buffer_2000_ddsd[120];

    auto g_xx_0_0_0_yy_xz_0_xy = buffer_2000_ddsd[121];

    auto g_xx_0_0_0_yy_xz_0_xz = buffer_2000_ddsd[122];

    auto g_xx_0_0_0_yy_xz_0_yy = buffer_2000_ddsd[123];

    auto g_xx_0_0_0_yy_xz_0_yz = buffer_2000_ddsd[124];

    auto g_xx_0_0_0_yy_xz_0_zz = buffer_2000_ddsd[125];

    auto g_xx_0_0_0_yy_yy_0_xx = buffer_2000_ddsd[126];

    auto g_xx_0_0_0_yy_yy_0_xy = buffer_2000_ddsd[127];

    auto g_xx_0_0_0_yy_yy_0_xz = buffer_2000_ddsd[128];

    auto g_xx_0_0_0_yy_yy_0_yy = buffer_2000_ddsd[129];

    auto g_xx_0_0_0_yy_yy_0_yz = buffer_2000_ddsd[130];

    auto g_xx_0_0_0_yy_yy_0_zz = buffer_2000_ddsd[131];

    auto g_xx_0_0_0_yy_yz_0_xx = buffer_2000_ddsd[132];

    auto g_xx_0_0_0_yy_yz_0_xy = buffer_2000_ddsd[133];

    auto g_xx_0_0_0_yy_yz_0_xz = buffer_2000_ddsd[134];

    auto g_xx_0_0_0_yy_yz_0_yy = buffer_2000_ddsd[135];

    auto g_xx_0_0_0_yy_yz_0_yz = buffer_2000_ddsd[136];

    auto g_xx_0_0_0_yy_yz_0_zz = buffer_2000_ddsd[137];

    auto g_xx_0_0_0_yy_zz_0_xx = buffer_2000_ddsd[138];

    auto g_xx_0_0_0_yy_zz_0_xy = buffer_2000_ddsd[139];

    auto g_xx_0_0_0_yy_zz_0_xz = buffer_2000_ddsd[140];

    auto g_xx_0_0_0_yy_zz_0_yy = buffer_2000_ddsd[141];

    auto g_xx_0_0_0_yy_zz_0_yz = buffer_2000_ddsd[142];

    auto g_xx_0_0_0_yy_zz_0_zz = buffer_2000_ddsd[143];

    auto g_xx_0_0_0_yz_xx_0_xx = buffer_2000_ddsd[144];

    auto g_xx_0_0_0_yz_xx_0_xy = buffer_2000_ddsd[145];

    auto g_xx_0_0_0_yz_xx_0_xz = buffer_2000_ddsd[146];

    auto g_xx_0_0_0_yz_xx_0_yy = buffer_2000_ddsd[147];

    auto g_xx_0_0_0_yz_xx_0_yz = buffer_2000_ddsd[148];

    auto g_xx_0_0_0_yz_xx_0_zz = buffer_2000_ddsd[149];

    auto g_xx_0_0_0_yz_xy_0_xx = buffer_2000_ddsd[150];

    auto g_xx_0_0_0_yz_xy_0_xy = buffer_2000_ddsd[151];

    auto g_xx_0_0_0_yz_xy_0_xz = buffer_2000_ddsd[152];

    auto g_xx_0_0_0_yz_xy_0_yy = buffer_2000_ddsd[153];

    auto g_xx_0_0_0_yz_xy_0_yz = buffer_2000_ddsd[154];

    auto g_xx_0_0_0_yz_xy_0_zz = buffer_2000_ddsd[155];

    auto g_xx_0_0_0_yz_xz_0_xx = buffer_2000_ddsd[156];

    auto g_xx_0_0_0_yz_xz_0_xy = buffer_2000_ddsd[157];

    auto g_xx_0_0_0_yz_xz_0_xz = buffer_2000_ddsd[158];

    auto g_xx_0_0_0_yz_xz_0_yy = buffer_2000_ddsd[159];

    auto g_xx_0_0_0_yz_xz_0_yz = buffer_2000_ddsd[160];

    auto g_xx_0_0_0_yz_xz_0_zz = buffer_2000_ddsd[161];

    auto g_xx_0_0_0_yz_yy_0_xx = buffer_2000_ddsd[162];

    auto g_xx_0_0_0_yz_yy_0_xy = buffer_2000_ddsd[163];

    auto g_xx_0_0_0_yz_yy_0_xz = buffer_2000_ddsd[164];

    auto g_xx_0_0_0_yz_yy_0_yy = buffer_2000_ddsd[165];

    auto g_xx_0_0_0_yz_yy_0_yz = buffer_2000_ddsd[166];

    auto g_xx_0_0_0_yz_yy_0_zz = buffer_2000_ddsd[167];

    auto g_xx_0_0_0_yz_yz_0_xx = buffer_2000_ddsd[168];

    auto g_xx_0_0_0_yz_yz_0_xy = buffer_2000_ddsd[169];

    auto g_xx_0_0_0_yz_yz_0_xz = buffer_2000_ddsd[170];

    auto g_xx_0_0_0_yz_yz_0_yy = buffer_2000_ddsd[171];

    auto g_xx_0_0_0_yz_yz_0_yz = buffer_2000_ddsd[172];

    auto g_xx_0_0_0_yz_yz_0_zz = buffer_2000_ddsd[173];

    auto g_xx_0_0_0_yz_zz_0_xx = buffer_2000_ddsd[174];

    auto g_xx_0_0_0_yz_zz_0_xy = buffer_2000_ddsd[175];

    auto g_xx_0_0_0_yz_zz_0_xz = buffer_2000_ddsd[176];

    auto g_xx_0_0_0_yz_zz_0_yy = buffer_2000_ddsd[177];

    auto g_xx_0_0_0_yz_zz_0_yz = buffer_2000_ddsd[178];

    auto g_xx_0_0_0_yz_zz_0_zz = buffer_2000_ddsd[179];

    auto g_xx_0_0_0_zz_xx_0_xx = buffer_2000_ddsd[180];

    auto g_xx_0_0_0_zz_xx_0_xy = buffer_2000_ddsd[181];

    auto g_xx_0_0_0_zz_xx_0_xz = buffer_2000_ddsd[182];

    auto g_xx_0_0_0_zz_xx_0_yy = buffer_2000_ddsd[183];

    auto g_xx_0_0_0_zz_xx_0_yz = buffer_2000_ddsd[184];

    auto g_xx_0_0_0_zz_xx_0_zz = buffer_2000_ddsd[185];

    auto g_xx_0_0_0_zz_xy_0_xx = buffer_2000_ddsd[186];

    auto g_xx_0_0_0_zz_xy_0_xy = buffer_2000_ddsd[187];

    auto g_xx_0_0_0_zz_xy_0_xz = buffer_2000_ddsd[188];

    auto g_xx_0_0_0_zz_xy_0_yy = buffer_2000_ddsd[189];

    auto g_xx_0_0_0_zz_xy_0_yz = buffer_2000_ddsd[190];

    auto g_xx_0_0_0_zz_xy_0_zz = buffer_2000_ddsd[191];

    auto g_xx_0_0_0_zz_xz_0_xx = buffer_2000_ddsd[192];

    auto g_xx_0_0_0_zz_xz_0_xy = buffer_2000_ddsd[193];

    auto g_xx_0_0_0_zz_xz_0_xz = buffer_2000_ddsd[194];

    auto g_xx_0_0_0_zz_xz_0_yy = buffer_2000_ddsd[195];

    auto g_xx_0_0_0_zz_xz_0_yz = buffer_2000_ddsd[196];

    auto g_xx_0_0_0_zz_xz_0_zz = buffer_2000_ddsd[197];

    auto g_xx_0_0_0_zz_yy_0_xx = buffer_2000_ddsd[198];

    auto g_xx_0_0_0_zz_yy_0_xy = buffer_2000_ddsd[199];

    auto g_xx_0_0_0_zz_yy_0_xz = buffer_2000_ddsd[200];

    auto g_xx_0_0_0_zz_yy_0_yy = buffer_2000_ddsd[201];

    auto g_xx_0_0_0_zz_yy_0_yz = buffer_2000_ddsd[202];

    auto g_xx_0_0_0_zz_yy_0_zz = buffer_2000_ddsd[203];

    auto g_xx_0_0_0_zz_yz_0_xx = buffer_2000_ddsd[204];

    auto g_xx_0_0_0_zz_yz_0_xy = buffer_2000_ddsd[205];

    auto g_xx_0_0_0_zz_yz_0_xz = buffer_2000_ddsd[206];

    auto g_xx_0_0_0_zz_yz_0_yy = buffer_2000_ddsd[207];

    auto g_xx_0_0_0_zz_yz_0_yz = buffer_2000_ddsd[208];

    auto g_xx_0_0_0_zz_yz_0_zz = buffer_2000_ddsd[209];

    auto g_xx_0_0_0_zz_zz_0_xx = buffer_2000_ddsd[210];

    auto g_xx_0_0_0_zz_zz_0_xy = buffer_2000_ddsd[211];

    auto g_xx_0_0_0_zz_zz_0_xz = buffer_2000_ddsd[212];

    auto g_xx_0_0_0_zz_zz_0_yy = buffer_2000_ddsd[213];

    auto g_xx_0_0_0_zz_zz_0_yz = buffer_2000_ddsd[214];

    auto g_xx_0_0_0_zz_zz_0_zz = buffer_2000_ddsd[215];

    auto g_xy_0_0_0_xx_xx_0_xx = buffer_2000_ddsd[216];

    auto g_xy_0_0_0_xx_xx_0_xy = buffer_2000_ddsd[217];

    auto g_xy_0_0_0_xx_xx_0_xz = buffer_2000_ddsd[218];

    auto g_xy_0_0_0_xx_xx_0_yy = buffer_2000_ddsd[219];

    auto g_xy_0_0_0_xx_xx_0_yz = buffer_2000_ddsd[220];

    auto g_xy_0_0_0_xx_xx_0_zz = buffer_2000_ddsd[221];

    auto g_xy_0_0_0_xx_xy_0_xx = buffer_2000_ddsd[222];

    auto g_xy_0_0_0_xx_xy_0_xy = buffer_2000_ddsd[223];

    auto g_xy_0_0_0_xx_xy_0_xz = buffer_2000_ddsd[224];

    auto g_xy_0_0_0_xx_xy_0_yy = buffer_2000_ddsd[225];

    auto g_xy_0_0_0_xx_xy_0_yz = buffer_2000_ddsd[226];

    auto g_xy_0_0_0_xx_xy_0_zz = buffer_2000_ddsd[227];

    auto g_xy_0_0_0_xx_xz_0_xx = buffer_2000_ddsd[228];

    auto g_xy_0_0_0_xx_xz_0_xy = buffer_2000_ddsd[229];

    auto g_xy_0_0_0_xx_xz_0_xz = buffer_2000_ddsd[230];

    auto g_xy_0_0_0_xx_xz_0_yy = buffer_2000_ddsd[231];

    auto g_xy_0_0_0_xx_xz_0_yz = buffer_2000_ddsd[232];

    auto g_xy_0_0_0_xx_xz_0_zz = buffer_2000_ddsd[233];

    auto g_xy_0_0_0_xx_yy_0_xx = buffer_2000_ddsd[234];

    auto g_xy_0_0_0_xx_yy_0_xy = buffer_2000_ddsd[235];

    auto g_xy_0_0_0_xx_yy_0_xz = buffer_2000_ddsd[236];

    auto g_xy_0_0_0_xx_yy_0_yy = buffer_2000_ddsd[237];

    auto g_xy_0_0_0_xx_yy_0_yz = buffer_2000_ddsd[238];

    auto g_xy_0_0_0_xx_yy_0_zz = buffer_2000_ddsd[239];

    auto g_xy_0_0_0_xx_yz_0_xx = buffer_2000_ddsd[240];

    auto g_xy_0_0_0_xx_yz_0_xy = buffer_2000_ddsd[241];

    auto g_xy_0_0_0_xx_yz_0_xz = buffer_2000_ddsd[242];

    auto g_xy_0_0_0_xx_yz_0_yy = buffer_2000_ddsd[243];

    auto g_xy_0_0_0_xx_yz_0_yz = buffer_2000_ddsd[244];

    auto g_xy_0_0_0_xx_yz_0_zz = buffer_2000_ddsd[245];

    auto g_xy_0_0_0_xx_zz_0_xx = buffer_2000_ddsd[246];

    auto g_xy_0_0_0_xx_zz_0_xy = buffer_2000_ddsd[247];

    auto g_xy_0_0_0_xx_zz_0_xz = buffer_2000_ddsd[248];

    auto g_xy_0_0_0_xx_zz_0_yy = buffer_2000_ddsd[249];

    auto g_xy_0_0_0_xx_zz_0_yz = buffer_2000_ddsd[250];

    auto g_xy_0_0_0_xx_zz_0_zz = buffer_2000_ddsd[251];

    auto g_xy_0_0_0_xy_xx_0_xx = buffer_2000_ddsd[252];

    auto g_xy_0_0_0_xy_xx_0_xy = buffer_2000_ddsd[253];

    auto g_xy_0_0_0_xy_xx_0_xz = buffer_2000_ddsd[254];

    auto g_xy_0_0_0_xy_xx_0_yy = buffer_2000_ddsd[255];

    auto g_xy_0_0_0_xy_xx_0_yz = buffer_2000_ddsd[256];

    auto g_xy_0_0_0_xy_xx_0_zz = buffer_2000_ddsd[257];

    auto g_xy_0_0_0_xy_xy_0_xx = buffer_2000_ddsd[258];

    auto g_xy_0_0_0_xy_xy_0_xy = buffer_2000_ddsd[259];

    auto g_xy_0_0_0_xy_xy_0_xz = buffer_2000_ddsd[260];

    auto g_xy_0_0_0_xy_xy_0_yy = buffer_2000_ddsd[261];

    auto g_xy_0_0_0_xy_xy_0_yz = buffer_2000_ddsd[262];

    auto g_xy_0_0_0_xy_xy_0_zz = buffer_2000_ddsd[263];

    auto g_xy_0_0_0_xy_xz_0_xx = buffer_2000_ddsd[264];

    auto g_xy_0_0_0_xy_xz_0_xy = buffer_2000_ddsd[265];

    auto g_xy_0_0_0_xy_xz_0_xz = buffer_2000_ddsd[266];

    auto g_xy_0_0_0_xy_xz_0_yy = buffer_2000_ddsd[267];

    auto g_xy_0_0_0_xy_xz_0_yz = buffer_2000_ddsd[268];

    auto g_xy_0_0_0_xy_xz_0_zz = buffer_2000_ddsd[269];

    auto g_xy_0_0_0_xy_yy_0_xx = buffer_2000_ddsd[270];

    auto g_xy_0_0_0_xy_yy_0_xy = buffer_2000_ddsd[271];

    auto g_xy_0_0_0_xy_yy_0_xz = buffer_2000_ddsd[272];

    auto g_xy_0_0_0_xy_yy_0_yy = buffer_2000_ddsd[273];

    auto g_xy_0_0_0_xy_yy_0_yz = buffer_2000_ddsd[274];

    auto g_xy_0_0_0_xy_yy_0_zz = buffer_2000_ddsd[275];

    auto g_xy_0_0_0_xy_yz_0_xx = buffer_2000_ddsd[276];

    auto g_xy_0_0_0_xy_yz_0_xy = buffer_2000_ddsd[277];

    auto g_xy_0_0_0_xy_yz_0_xz = buffer_2000_ddsd[278];

    auto g_xy_0_0_0_xy_yz_0_yy = buffer_2000_ddsd[279];

    auto g_xy_0_0_0_xy_yz_0_yz = buffer_2000_ddsd[280];

    auto g_xy_0_0_0_xy_yz_0_zz = buffer_2000_ddsd[281];

    auto g_xy_0_0_0_xy_zz_0_xx = buffer_2000_ddsd[282];

    auto g_xy_0_0_0_xy_zz_0_xy = buffer_2000_ddsd[283];

    auto g_xy_0_0_0_xy_zz_0_xz = buffer_2000_ddsd[284];

    auto g_xy_0_0_0_xy_zz_0_yy = buffer_2000_ddsd[285];

    auto g_xy_0_0_0_xy_zz_0_yz = buffer_2000_ddsd[286];

    auto g_xy_0_0_0_xy_zz_0_zz = buffer_2000_ddsd[287];

    auto g_xy_0_0_0_xz_xx_0_xx = buffer_2000_ddsd[288];

    auto g_xy_0_0_0_xz_xx_0_xy = buffer_2000_ddsd[289];

    auto g_xy_0_0_0_xz_xx_0_xz = buffer_2000_ddsd[290];

    auto g_xy_0_0_0_xz_xx_0_yy = buffer_2000_ddsd[291];

    auto g_xy_0_0_0_xz_xx_0_yz = buffer_2000_ddsd[292];

    auto g_xy_0_0_0_xz_xx_0_zz = buffer_2000_ddsd[293];

    auto g_xy_0_0_0_xz_xy_0_xx = buffer_2000_ddsd[294];

    auto g_xy_0_0_0_xz_xy_0_xy = buffer_2000_ddsd[295];

    auto g_xy_0_0_0_xz_xy_0_xz = buffer_2000_ddsd[296];

    auto g_xy_0_0_0_xz_xy_0_yy = buffer_2000_ddsd[297];

    auto g_xy_0_0_0_xz_xy_0_yz = buffer_2000_ddsd[298];

    auto g_xy_0_0_0_xz_xy_0_zz = buffer_2000_ddsd[299];

    auto g_xy_0_0_0_xz_xz_0_xx = buffer_2000_ddsd[300];

    auto g_xy_0_0_0_xz_xz_0_xy = buffer_2000_ddsd[301];

    auto g_xy_0_0_0_xz_xz_0_xz = buffer_2000_ddsd[302];

    auto g_xy_0_0_0_xz_xz_0_yy = buffer_2000_ddsd[303];

    auto g_xy_0_0_0_xz_xz_0_yz = buffer_2000_ddsd[304];

    auto g_xy_0_0_0_xz_xz_0_zz = buffer_2000_ddsd[305];

    auto g_xy_0_0_0_xz_yy_0_xx = buffer_2000_ddsd[306];

    auto g_xy_0_0_0_xz_yy_0_xy = buffer_2000_ddsd[307];

    auto g_xy_0_0_0_xz_yy_0_xz = buffer_2000_ddsd[308];

    auto g_xy_0_0_0_xz_yy_0_yy = buffer_2000_ddsd[309];

    auto g_xy_0_0_0_xz_yy_0_yz = buffer_2000_ddsd[310];

    auto g_xy_0_0_0_xz_yy_0_zz = buffer_2000_ddsd[311];

    auto g_xy_0_0_0_xz_yz_0_xx = buffer_2000_ddsd[312];

    auto g_xy_0_0_0_xz_yz_0_xy = buffer_2000_ddsd[313];

    auto g_xy_0_0_0_xz_yz_0_xz = buffer_2000_ddsd[314];

    auto g_xy_0_0_0_xz_yz_0_yy = buffer_2000_ddsd[315];

    auto g_xy_0_0_0_xz_yz_0_yz = buffer_2000_ddsd[316];

    auto g_xy_0_0_0_xz_yz_0_zz = buffer_2000_ddsd[317];

    auto g_xy_0_0_0_xz_zz_0_xx = buffer_2000_ddsd[318];

    auto g_xy_0_0_0_xz_zz_0_xy = buffer_2000_ddsd[319];

    auto g_xy_0_0_0_xz_zz_0_xz = buffer_2000_ddsd[320];

    auto g_xy_0_0_0_xz_zz_0_yy = buffer_2000_ddsd[321];

    auto g_xy_0_0_0_xz_zz_0_yz = buffer_2000_ddsd[322];

    auto g_xy_0_0_0_xz_zz_0_zz = buffer_2000_ddsd[323];

    auto g_xy_0_0_0_yy_xx_0_xx = buffer_2000_ddsd[324];

    auto g_xy_0_0_0_yy_xx_0_xy = buffer_2000_ddsd[325];

    auto g_xy_0_0_0_yy_xx_0_xz = buffer_2000_ddsd[326];

    auto g_xy_0_0_0_yy_xx_0_yy = buffer_2000_ddsd[327];

    auto g_xy_0_0_0_yy_xx_0_yz = buffer_2000_ddsd[328];

    auto g_xy_0_0_0_yy_xx_0_zz = buffer_2000_ddsd[329];

    auto g_xy_0_0_0_yy_xy_0_xx = buffer_2000_ddsd[330];

    auto g_xy_0_0_0_yy_xy_0_xy = buffer_2000_ddsd[331];

    auto g_xy_0_0_0_yy_xy_0_xz = buffer_2000_ddsd[332];

    auto g_xy_0_0_0_yy_xy_0_yy = buffer_2000_ddsd[333];

    auto g_xy_0_0_0_yy_xy_0_yz = buffer_2000_ddsd[334];

    auto g_xy_0_0_0_yy_xy_0_zz = buffer_2000_ddsd[335];

    auto g_xy_0_0_0_yy_xz_0_xx = buffer_2000_ddsd[336];

    auto g_xy_0_0_0_yy_xz_0_xy = buffer_2000_ddsd[337];

    auto g_xy_0_0_0_yy_xz_0_xz = buffer_2000_ddsd[338];

    auto g_xy_0_0_0_yy_xz_0_yy = buffer_2000_ddsd[339];

    auto g_xy_0_0_0_yy_xz_0_yz = buffer_2000_ddsd[340];

    auto g_xy_0_0_0_yy_xz_0_zz = buffer_2000_ddsd[341];

    auto g_xy_0_0_0_yy_yy_0_xx = buffer_2000_ddsd[342];

    auto g_xy_0_0_0_yy_yy_0_xy = buffer_2000_ddsd[343];

    auto g_xy_0_0_0_yy_yy_0_xz = buffer_2000_ddsd[344];

    auto g_xy_0_0_0_yy_yy_0_yy = buffer_2000_ddsd[345];

    auto g_xy_0_0_0_yy_yy_0_yz = buffer_2000_ddsd[346];

    auto g_xy_0_0_0_yy_yy_0_zz = buffer_2000_ddsd[347];

    auto g_xy_0_0_0_yy_yz_0_xx = buffer_2000_ddsd[348];

    auto g_xy_0_0_0_yy_yz_0_xy = buffer_2000_ddsd[349];

    auto g_xy_0_0_0_yy_yz_0_xz = buffer_2000_ddsd[350];

    auto g_xy_0_0_0_yy_yz_0_yy = buffer_2000_ddsd[351];

    auto g_xy_0_0_0_yy_yz_0_yz = buffer_2000_ddsd[352];

    auto g_xy_0_0_0_yy_yz_0_zz = buffer_2000_ddsd[353];

    auto g_xy_0_0_0_yy_zz_0_xx = buffer_2000_ddsd[354];

    auto g_xy_0_0_0_yy_zz_0_xy = buffer_2000_ddsd[355];

    auto g_xy_0_0_0_yy_zz_0_xz = buffer_2000_ddsd[356];

    auto g_xy_0_0_0_yy_zz_0_yy = buffer_2000_ddsd[357];

    auto g_xy_0_0_0_yy_zz_0_yz = buffer_2000_ddsd[358];

    auto g_xy_0_0_0_yy_zz_0_zz = buffer_2000_ddsd[359];

    auto g_xy_0_0_0_yz_xx_0_xx = buffer_2000_ddsd[360];

    auto g_xy_0_0_0_yz_xx_0_xy = buffer_2000_ddsd[361];

    auto g_xy_0_0_0_yz_xx_0_xz = buffer_2000_ddsd[362];

    auto g_xy_0_0_0_yz_xx_0_yy = buffer_2000_ddsd[363];

    auto g_xy_0_0_0_yz_xx_0_yz = buffer_2000_ddsd[364];

    auto g_xy_0_0_0_yz_xx_0_zz = buffer_2000_ddsd[365];

    auto g_xy_0_0_0_yz_xy_0_xx = buffer_2000_ddsd[366];

    auto g_xy_0_0_0_yz_xy_0_xy = buffer_2000_ddsd[367];

    auto g_xy_0_0_0_yz_xy_0_xz = buffer_2000_ddsd[368];

    auto g_xy_0_0_0_yz_xy_0_yy = buffer_2000_ddsd[369];

    auto g_xy_0_0_0_yz_xy_0_yz = buffer_2000_ddsd[370];

    auto g_xy_0_0_0_yz_xy_0_zz = buffer_2000_ddsd[371];

    auto g_xy_0_0_0_yz_xz_0_xx = buffer_2000_ddsd[372];

    auto g_xy_0_0_0_yz_xz_0_xy = buffer_2000_ddsd[373];

    auto g_xy_0_0_0_yz_xz_0_xz = buffer_2000_ddsd[374];

    auto g_xy_0_0_0_yz_xz_0_yy = buffer_2000_ddsd[375];

    auto g_xy_0_0_0_yz_xz_0_yz = buffer_2000_ddsd[376];

    auto g_xy_0_0_0_yz_xz_0_zz = buffer_2000_ddsd[377];

    auto g_xy_0_0_0_yz_yy_0_xx = buffer_2000_ddsd[378];

    auto g_xy_0_0_0_yz_yy_0_xy = buffer_2000_ddsd[379];

    auto g_xy_0_0_0_yz_yy_0_xz = buffer_2000_ddsd[380];

    auto g_xy_0_0_0_yz_yy_0_yy = buffer_2000_ddsd[381];

    auto g_xy_0_0_0_yz_yy_0_yz = buffer_2000_ddsd[382];

    auto g_xy_0_0_0_yz_yy_0_zz = buffer_2000_ddsd[383];

    auto g_xy_0_0_0_yz_yz_0_xx = buffer_2000_ddsd[384];

    auto g_xy_0_0_0_yz_yz_0_xy = buffer_2000_ddsd[385];

    auto g_xy_0_0_0_yz_yz_0_xz = buffer_2000_ddsd[386];

    auto g_xy_0_0_0_yz_yz_0_yy = buffer_2000_ddsd[387];

    auto g_xy_0_0_0_yz_yz_0_yz = buffer_2000_ddsd[388];

    auto g_xy_0_0_0_yz_yz_0_zz = buffer_2000_ddsd[389];

    auto g_xy_0_0_0_yz_zz_0_xx = buffer_2000_ddsd[390];

    auto g_xy_0_0_0_yz_zz_0_xy = buffer_2000_ddsd[391];

    auto g_xy_0_0_0_yz_zz_0_xz = buffer_2000_ddsd[392];

    auto g_xy_0_0_0_yz_zz_0_yy = buffer_2000_ddsd[393];

    auto g_xy_0_0_0_yz_zz_0_yz = buffer_2000_ddsd[394];

    auto g_xy_0_0_0_yz_zz_0_zz = buffer_2000_ddsd[395];

    auto g_xy_0_0_0_zz_xx_0_xx = buffer_2000_ddsd[396];

    auto g_xy_0_0_0_zz_xx_0_xy = buffer_2000_ddsd[397];

    auto g_xy_0_0_0_zz_xx_0_xz = buffer_2000_ddsd[398];

    auto g_xy_0_0_0_zz_xx_0_yy = buffer_2000_ddsd[399];

    auto g_xy_0_0_0_zz_xx_0_yz = buffer_2000_ddsd[400];

    auto g_xy_0_0_0_zz_xx_0_zz = buffer_2000_ddsd[401];

    auto g_xy_0_0_0_zz_xy_0_xx = buffer_2000_ddsd[402];

    auto g_xy_0_0_0_zz_xy_0_xy = buffer_2000_ddsd[403];

    auto g_xy_0_0_0_zz_xy_0_xz = buffer_2000_ddsd[404];

    auto g_xy_0_0_0_zz_xy_0_yy = buffer_2000_ddsd[405];

    auto g_xy_0_0_0_zz_xy_0_yz = buffer_2000_ddsd[406];

    auto g_xy_0_0_0_zz_xy_0_zz = buffer_2000_ddsd[407];

    auto g_xy_0_0_0_zz_xz_0_xx = buffer_2000_ddsd[408];

    auto g_xy_0_0_0_zz_xz_0_xy = buffer_2000_ddsd[409];

    auto g_xy_0_0_0_zz_xz_0_xz = buffer_2000_ddsd[410];

    auto g_xy_0_0_0_zz_xz_0_yy = buffer_2000_ddsd[411];

    auto g_xy_0_0_0_zz_xz_0_yz = buffer_2000_ddsd[412];

    auto g_xy_0_0_0_zz_xz_0_zz = buffer_2000_ddsd[413];

    auto g_xy_0_0_0_zz_yy_0_xx = buffer_2000_ddsd[414];

    auto g_xy_0_0_0_zz_yy_0_xy = buffer_2000_ddsd[415];

    auto g_xy_0_0_0_zz_yy_0_xz = buffer_2000_ddsd[416];

    auto g_xy_0_0_0_zz_yy_0_yy = buffer_2000_ddsd[417];

    auto g_xy_0_0_0_zz_yy_0_yz = buffer_2000_ddsd[418];

    auto g_xy_0_0_0_zz_yy_0_zz = buffer_2000_ddsd[419];

    auto g_xy_0_0_0_zz_yz_0_xx = buffer_2000_ddsd[420];

    auto g_xy_0_0_0_zz_yz_0_xy = buffer_2000_ddsd[421];

    auto g_xy_0_0_0_zz_yz_0_xz = buffer_2000_ddsd[422];

    auto g_xy_0_0_0_zz_yz_0_yy = buffer_2000_ddsd[423];

    auto g_xy_0_0_0_zz_yz_0_yz = buffer_2000_ddsd[424];

    auto g_xy_0_0_0_zz_yz_0_zz = buffer_2000_ddsd[425];

    auto g_xy_0_0_0_zz_zz_0_xx = buffer_2000_ddsd[426];

    auto g_xy_0_0_0_zz_zz_0_xy = buffer_2000_ddsd[427];

    auto g_xy_0_0_0_zz_zz_0_xz = buffer_2000_ddsd[428];

    auto g_xy_0_0_0_zz_zz_0_yy = buffer_2000_ddsd[429];

    auto g_xy_0_0_0_zz_zz_0_yz = buffer_2000_ddsd[430];

    auto g_xy_0_0_0_zz_zz_0_zz = buffer_2000_ddsd[431];

    auto g_xz_0_0_0_xx_xx_0_xx = buffer_2000_ddsd[432];

    auto g_xz_0_0_0_xx_xx_0_xy = buffer_2000_ddsd[433];

    auto g_xz_0_0_0_xx_xx_0_xz = buffer_2000_ddsd[434];

    auto g_xz_0_0_0_xx_xx_0_yy = buffer_2000_ddsd[435];

    auto g_xz_0_0_0_xx_xx_0_yz = buffer_2000_ddsd[436];

    auto g_xz_0_0_0_xx_xx_0_zz = buffer_2000_ddsd[437];

    auto g_xz_0_0_0_xx_xy_0_xx = buffer_2000_ddsd[438];

    auto g_xz_0_0_0_xx_xy_0_xy = buffer_2000_ddsd[439];

    auto g_xz_0_0_0_xx_xy_0_xz = buffer_2000_ddsd[440];

    auto g_xz_0_0_0_xx_xy_0_yy = buffer_2000_ddsd[441];

    auto g_xz_0_0_0_xx_xy_0_yz = buffer_2000_ddsd[442];

    auto g_xz_0_0_0_xx_xy_0_zz = buffer_2000_ddsd[443];

    auto g_xz_0_0_0_xx_xz_0_xx = buffer_2000_ddsd[444];

    auto g_xz_0_0_0_xx_xz_0_xy = buffer_2000_ddsd[445];

    auto g_xz_0_0_0_xx_xz_0_xz = buffer_2000_ddsd[446];

    auto g_xz_0_0_0_xx_xz_0_yy = buffer_2000_ddsd[447];

    auto g_xz_0_0_0_xx_xz_0_yz = buffer_2000_ddsd[448];

    auto g_xz_0_0_0_xx_xz_0_zz = buffer_2000_ddsd[449];

    auto g_xz_0_0_0_xx_yy_0_xx = buffer_2000_ddsd[450];

    auto g_xz_0_0_0_xx_yy_0_xy = buffer_2000_ddsd[451];

    auto g_xz_0_0_0_xx_yy_0_xz = buffer_2000_ddsd[452];

    auto g_xz_0_0_0_xx_yy_0_yy = buffer_2000_ddsd[453];

    auto g_xz_0_0_0_xx_yy_0_yz = buffer_2000_ddsd[454];

    auto g_xz_0_0_0_xx_yy_0_zz = buffer_2000_ddsd[455];

    auto g_xz_0_0_0_xx_yz_0_xx = buffer_2000_ddsd[456];

    auto g_xz_0_0_0_xx_yz_0_xy = buffer_2000_ddsd[457];

    auto g_xz_0_0_0_xx_yz_0_xz = buffer_2000_ddsd[458];

    auto g_xz_0_0_0_xx_yz_0_yy = buffer_2000_ddsd[459];

    auto g_xz_0_0_0_xx_yz_0_yz = buffer_2000_ddsd[460];

    auto g_xz_0_0_0_xx_yz_0_zz = buffer_2000_ddsd[461];

    auto g_xz_0_0_0_xx_zz_0_xx = buffer_2000_ddsd[462];

    auto g_xz_0_0_0_xx_zz_0_xy = buffer_2000_ddsd[463];

    auto g_xz_0_0_0_xx_zz_0_xz = buffer_2000_ddsd[464];

    auto g_xz_0_0_0_xx_zz_0_yy = buffer_2000_ddsd[465];

    auto g_xz_0_0_0_xx_zz_0_yz = buffer_2000_ddsd[466];

    auto g_xz_0_0_0_xx_zz_0_zz = buffer_2000_ddsd[467];

    auto g_xz_0_0_0_xy_xx_0_xx = buffer_2000_ddsd[468];

    auto g_xz_0_0_0_xy_xx_0_xy = buffer_2000_ddsd[469];

    auto g_xz_0_0_0_xy_xx_0_xz = buffer_2000_ddsd[470];

    auto g_xz_0_0_0_xy_xx_0_yy = buffer_2000_ddsd[471];

    auto g_xz_0_0_0_xy_xx_0_yz = buffer_2000_ddsd[472];

    auto g_xz_0_0_0_xy_xx_0_zz = buffer_2000_ddsd[473];

    auto g_xz_0_0_0_xy_xy_0_xx = buffer_2000_ddsd[474];

    auto g_xz_0_0_0_xy_xy_0_xy = buffer_2000_ddsd[475];

    auto g_xz_0_0_0_xy_xy_0_xz = buffer_2000_ddsd[476];

    auto g_xz_0_0_0_xy_xy_0_yy = buffer_2000_ddsd[477];

    auto g_xz_0_0_0_xy_xy_0_yz = buffer_2000_ddsd[478];

    auto g_xz_0_0_0_xy_xy_0_zz = buffer_2000_ddsd[479];

    auto g_xz_0_0_0_xy_xz_0_xx = buffer_2000_ddsd[480];

    auto g_xz_0_0_0_xy_xz_0_xy = buffer_2000_ddsd[481];

    auto g_xz_0_0_0_xy_xz_0_xz = buffer_2000_ddsd[482];

    auto g_xz_0_0_0_xy_xz_0_yy = buffer_2000_ddsd[483];

    auto g_xz_0_0_0_xy_xz_0_yz = buffer_2000_ddsd[484];

    auto g_xz_0_0_0_xy_xz_0_zz = buffer_2000_ddsd[485];

    auto g_xz_0_0_0_xy_yy_0_xx = buffer_2000_ddsd[486];

    auto g_xz_0_0_0_xy_yy_0_xy = buffer_2000_ddsd[487];

    auto g_xz_0_0_0_xy_yy_0_xz = buffer_2000_ddsd[488];

    auto g_xz_0_0_0_xy_yy_0_yy = buffer_2000_ddsd[489];

    auto g_xz_0_0_0_xy_yy_0_yz = buffer_2000_ddsd[490];

    auto g_xz_0_0_0_xy_yy_0_zz = buffer_2000_ddsd[491];

    auto g_xz_0_0_0_xy_yz_0_xx = buffer_2000_ddsd[492];

    auto g_xz_0_0_0_xy_yz_0_xy = buffer_2000_ddsd[493];

    auto g_xz_0_0_0_xy_yz_0_xz = buffer_2000_ddsd[494];

    auto g_xz_0_0_0_xy_yz_0_yy = buffer_2000_ddsd[495];

    auto g_xz_0_0_0_xy_yz_0_yz = buffer_2000_ddsd[496];

    auto g_xz_0_0_0_xy_yz_0_zz = buffer_2000_ddsd[497];

    auto g_xz_0_0_0_xy_zz_0_xx = buffer_2000_ddsd[498];

    auto g_xz_0_0_0_xy_zz_0_xy = buffer_2000_ddsd[499];

    auto g_xz_0_0_0_xy_zz_0_xz = buffer_2000_ddsd[500];

    auto g_xz_0_0_0_xy_zz_0_yy = buffer_2000_ddsd[501];

    auto g_xz_0_0_0_xy_zz_0_yz = buffer_2000_ddsd[502];

    auto g_xz_0_0_0_xy_zz_0_zz = buffer_2000_ddsd[503];

    auto g_xz_0_0_0_xz_xx_0_xx = buffer_2000_ddsd[504];

    auto g_xz_0_0_0_xz_xx_0_xy = buffer_2000_ddsd[505];

    auto g_xz_0_0_0_xz_xx_0_xz = buffer_2000_ddsd[506];

    auto g_xz_0_0_0_xz_xx_0_yy = buffer_2000_ddsd[507];

    auto g_xz_0_0_0_xz_xx_0_yz = buffer_2000_ddsd[508];

    auto g_xz_0_0_0_xz_xx_0_zz = buffer_2000_ddsd[509];

    auto g_xz_0_0_0_xz_xy_0_xx = buffer_2000_ddsd[510];

    auto g_xz_0_0_0_xz_xy_0_xy = buffer_2000_ddsd[511];

    auto g_xz_0_0_0_xz_xy_0_xz = buffer_2000_ddsd[512];

    auto g_xz_0_0_0_xz_xy_0_yy = buffer_2000_ddsd[513];

    auto g_xz_0_0_0_xz_xy_0_yz = buffer_2000_ddsd[514];

    auto g_xz_0_0_0_xz_xy_0_zz = buffer_2000_ddsd[515];

    auto g_xz_0_0_0_xz_xz_0_xx = buffer_2000_ddsd[516];

    auto g_xz_0_0_0_xz_xz_0_xy = buffer_2000_ddsd[517];

    auto g_xz_0_0_0_xz_xz_0_xz = buffer_2000_ddsd[518];

    auto g_xz_0_0_0_xz_xz_0_yy = buffer_2000_ddsd[519];

    auto g_xz_0_0_0_xz_xz_0_yz = buffer_2000_ddsd[520];

    auto g_xz_0_0_0_xz_xz_0_zz = buffer_2000_ddsd[521];

    auto g_xz_0_0_0_xz_yy_0_xx = buffer_2000_ddsd[522];

    auto g_xz_0_0_0_xz_yy_0_xy = buffer_2000_ddsd[523];

    auto g_xz_0_0_0_xz_yy_0_xz = buffer_2000_ddsd[524];

    auto g_xz_0_0_0_xz_yy_0_yy = buffer_2000_ddsd[525];

    auto g_xz_0_0_0_xz_yy_0_yz = buffer_2000_ddsd[526];

    auto g_xz_0_0_0_xz_yy_0_zz = buffer_2000_ddsd[527];

    auto g_xz_0_0_0_xz_yz_0_xx = buffer_2000_ddsd[528];

    auto g_xz_0_0_0_xz_yz_0_xy = buffer_2000_ddsd[529];

    auto g_xz_0_0_0_xz_yz_0_xz = buffer_2000_ddsd[530];

    auto g_xz_0_0_0_xz_yz_0_yy = buffer_2000_ddsd[531];

    auto g_xz_0_0_0_xz_yz_0_yz = buffer_2000_ddsd[532];

    auto g_xz_0_0_0_xz_yz_0_zz = buffer_2000_ddsd[533];

    auto g_xz_0_0_0_xz_zz_0_xx = buffer_2000_ddsd[534];

    auto g_xz_0_0_0_xz_zz_0_xy = buffer_2000_ddsd[535];

    auto g_xz_0_0_0_xz_zz_0_xz = buffer_2000_ddsd[536];

    auto g_xz_0_0_0_xz_zz_0_yy = buffer_2000_ddsd[537];

    auto g_xz_0_0_0_xz_zz_0_yz = buffer_2000_ddsd[538];

    auto g_xz_0_0_0_xz_zz_0_zz = buffer_2000_ddsd[539];

    auto g_xz_0_0_0_yy_xx_0_xx = buffer_2000_ddsd[540];

    auto g_xz_0_0_0_yy_xx_0_xy = buffer_2000_ddsd[541];

    auto g_xz_0_0_0_yy_xx_0_xz = buffer_2000_ddsd[542];

    auto g_xz_0_0_0_yy_xx_0_yy = buffer_2000_ddsd[543];

    auto g_xz_0_0_0_yy_xx_0_yz = buffer_2000_ddsd[544];

    auto g_xz_0_0_0_yy_xx_0_zz = buffer_2000_ddsd[545];

    auto g_xz_0_0_0_yy_xy_0_xx = buffer_2000_ddsd[546];

    auto g_xz_0_0_0_yy_xy_0_xy = buffer_2000_ddsd[547];

    auto g_xz_0_0_0_yy_xy_0_xz = buffer_2000_ddsd[548];

    auto g_xz_0_0_0_yy_xy_0_yy = buffer_2000_ddsd[549];

    auto g_xz_0_0_0_yy_xy_0_yz = buffer_2000_ddsd[550];

    auto g_xz_0_0_0_yy_xy_0_zz = buffer_2000_ddsd[551];

    auto g_xz_0_0_0_yy_xz_0_xx = buffer_2000_ddsd[552];

    auto g_xz_0_0_0_yy_xz_0_xy = buffer_2000_ddsd[553];

    auto g_xz_0_0_0_yy_xz_0_xz = buffer_2000_ddsd[554];

    auto g_xz_0_0_0_yy_xz_0_yy = buffer_2000_ddsd[555];

    auto g_xz_0_0_0_yy_xz_0_yz = buffer_2000_ddsd[556];

    auto g_xz_0_0_0_yy_xz_0_zz = buffer_2000_ddsd[557];

    auto g_xz_0_0_0_yy_yy_0_xx = buffer_2000_ddsd[558];

    auto g_xz_0_0_0_yy_yy_0_xy = buffer_2000_ddsd[559];

    auto g_xz_0_0_0_yy_yy_0_xz = buffer_2000_ddsd[560];

    auto g_xz_0_0_0_yy_yy_0_yy = buffer_2000_ddsd[561];

    auto g_xz_0_0_0_yy_yy_0_yz = buffer_2000_ddsd[562];

    auto g_xz_0_0_0_yy_yy_0_zz = buffer_2000_ddsd[563];

    auto g_xz_0_0_0_yy_yz_0_xx = buffer_2000_ddsd[564];

    auto g_xz_0_0_0_yy_yz_0_xy = buffer_2000_ddsd[565];

    auto g_xz_0_0_0_yy_yz_0_xz = buffer_2000_ddsd[566];

    auto g_xz_0_0_0_yy_yz_0_yy = buffer_2000_ddsd[567];

    auto g_xz_0_0_0_yy_yz_0_yz = buffer_2000_ddsd[568];

    auto g_xz_0_0_0_yy_yz_0_zz = buffer_2000_ddsd[569];

    auto g_xz_0_0_0_yy_zz_0_xx = buffer_2000_ddsd[570];

    auto g_xz_0_0_0_yy_zz_0_xy = buffer_2000_ddsd[571];

    auto g_xz_0_0_0_yy_zz_0_xz = buffer_2000_ddsd[572];

    auto g_xz_0_0_0_yy_zz_0_yy = buffer_2000_ddsd[573];

    auto g_xz_0_0_0_yy_zz_0_yz = buffer_2000_ddsd[574];

    auto g_xz_0_0_0_yy_zz_0_zz = buffer_2000_ddsd[575];

    auto g_xz_0_0_0_yz_xx_0_xx = buffer_2000_ddsd[576];

    auto g_xz_0_0_0_yz_xx_0_xy = buffer_2000_ddsd[577];

    auto g_xz_0_0_0_yz_xx_0_xz = buffer_2000_ddsd[578];

    auto g_xz_0_0_0_yz_xx_0_yy = buffer_2000_ddsd[579];

    auto g_xz_0_0_0_yz_xx_0_yz = buffer_2000_ddsd[580];

    auto g_xz_0_0_0_yz_xx_0_zz = buffer_2000_ddsd[581];

    auto g_xz_0_0_0_yz_xy_0_xx = buffer_2000_ddsd[582];

    auto g_xz_0_0_0_yz_xy_0_xy = buffer_2000_ddsd[583];

    auto g_xz_0_0_0_yz_xy_0_xz = buffer_2000_ddsd[584];

    auto g_xz_0_0_0_yz_xy_0_yy = buffer_2000_ddsd[585];

    auto g_xz_0_0_0_yz_xy_0_yz = buffer_2000_ddsd[586];

    auto g_xz_0_0_0_yz_xy_0_zz = buffer_2000_ddsd[587];

    auto g_xz_0_0_0_yz_xz_0_xx = buffer_2000_ddsd[588];

    auto g_xz_0_0_0_yz_xz_0_xy = buffer_2000_ddsd[589];

    auto g_xz_0_0_0_yz_xz_0_xz = buffer_2000_ddsd[590];

    auto g_xz_0_0_0_yz_xz_0_yy = buffer_2000_ddsd[591];

    auto g_xz_0_0_0_yz_xz_0_yz = buffer_2000_ddsd[592];

    auto g_xz_0_0_0_yz_xz_0_zz = buffer_2000_ddsd[593];

    auto g_xz_0_0_0_yz_yy_0_xx = buffer_2000_ddsd[594];

    auto g_xz_0_0_0_yz_yy_0_xy = buffer_2000_ddsd[595];

    auto g_xz_0_0_0_yz_yy_0_xz = buffer_2000_ddsd[596];

    auto g_xz_0_0_0_yz_yy_0_yy = buffer_2000_ddsd[597];

    auto g_xz_0_0_0_yz_yy_0_yz = buffer_2000_ddsd[598];

    auto g_xz_0_0_0_yz_yy_0_zz = buffer_2000_ddsd[599];

    auto g_xz_0_0_0_yz_yz_0_xx = buffer_2000_ddsd[600];

    auto g_xz_0_0_0_yz_yz_0_xy = buffer_2000_ddsd[601];

    auto g_xz_0_0_0_yz_yz_0_xz = buffer_2000_ddsd[602];

    auto g_xz_0_0_0_yz_yz_0_yy = buffer_2000_ddsd[603];

    auto g_xz_0_0_0_yz_yz_0_yz = buffer_2000_ddsd[604];

    auto g_xz_0_0_0_yz_yz_0_zz = buffer_2000_ddsd[605];

    auto g_xz_0_0_0_yz_zz_0_xx = buffer_2000_ddsd[606];

    auto g_xz_0_0_0_yz_zz_0_xy = buffer_2000_ddsd[607];

    auto g_xz_0_0_0_yz_zz_0_xz = buffer_2000_ddsd[608];

    auto g_xz_0_0_0_yz_zz_0_yy = buffer_2000_ddsd[609];

    auto g_xz_0_0_0_yz_zz_0_yz = buffer_2000_ddsd[610];

    auto g_xz_0_0_0_yz_zz_0_zz = buffer_2000_ddsd[611];

    auto g_xz_0_0_0_zz_xx_0_xx = buffer_2000_ddsd[612];

    auto g_xz_0_0_0_zz_xx_0_xy = buffer_2000_ddsd[613];

    auto g_xz_0_0_0_zz_xx_0_xz = buffer_2000_ddsd[614];

    auto g_xz_0_0_0_zz_xx_0_yy = buffer_2000_ddsd[615];

    auto g_xz_0_0_0_zz_xx_0_yz = buffer_2000_ddsd[616];

    auto g_xz_0_0_0_zz_xx_0_zz = buffer_2000_ddsd[617];

    auto g_xz_0_0_0_zz_xy_0_xx = buffer_2000_ddsd[618];

    auto g_xz_0_0_0_zz_xy_0_xy = buffer_2000_ddsd[619];

    auto g_xz_0_0_0_zz_xy_0_xz = buffer_2000_ddsd[620];

    auto g_xz_0_0_0_zz_xy_0_yy = buffer_2000_ddsd[621];

    auto g_xz_0_0_0_zz_xy_0_yz = buffer_2000_ddsd[622];

    auto g_xz_0_0_0_zz_xy_0_zz = buffer_2000_ddsd[623];

    auto g_xz_0_0_0_zz_xz_0_xx = buffer_2000_ddsd[624];

    auto g_xz_0_0_0_zz_xz_0_xy = buffer_2000_ddsd[625];

    auto g_xz_0_0_0_zz_xz_0_xz = buffer_2000_ddsd[626];

    auto g_xz_0_0_0_zz_xz_0_yy = buffer_2000_ddsd[627];

    auto g_xz_0_0_0_zz_xz_0_yz = buffer_2000_ddsd[628];

    auto g_xz_0_0_0_zz_xz_0_zz = buffer_2000_ddsd[629];

    auto g_xz_0_0_0_zz_yy_0_xx = buffer_2000_ddsd[630];

    auto g_xz_0_0_0_zz_yy_0_xy = buffer_2000_ddsd[631];

    auto g_xz_0_0_0_zz_yy_0_xz = buffer_2000_ddsd[632];

    auto g_xz_0_0_0_zz_yy_0_yy = buffer_2000_ddsd[633];

    auto g_xz_0_0_0_zz_yy_0_yz = buffer_2000_ddsd[634];

    auto g_xz_0_0_0_zz_yy_0_zz = buffer_2000_ddsd[635];

    auto g_xz_0_0_0_zz_yz_0_xx = buffer_2000_ddsd[636];

    auto g_xz_0_0_0_zz_yz_0_xy = buffer_2000_ddsd[637];

    auto g_xz_0_0_0_zz_yz_0_xz = buffer_2000_ddsd[638];

    auto g_xz_0_0_0_zz_yz_0_yy = buffer_2000_ddsd[639];

    auto g_xz_0_0_0_zz_yz_0_yz = buffer_2000_ddsd[640];

    auto g_xz_0_0_0_zz_yz_0_zz = buffer_2000_ddsd[641];

    auto g_xz_0_0_0_zz_zz_0_xx = buffer_2000_ddsd[642];

    auto g_xz_0_0_0_zz_zz_0_xy = buffer_2000_ddsd[643];

    auto g_xz_0_0_0_zz_zz_0_xz = buffer_2000_ddsd[644];

    auto g_xz_0_0_0_zz_zz_0_yy = buffer_2000_ddsd[645];

    auto g_xz_0_0_0_zz_zz_0_yz = buffer_2000_ddsd[646];

    auto g_xz_0_0_0_zz_zz_0_zz = buffer_2000_ddsd[647];

    auto g_yy_0_0_0_xx_xx_0_xx = buffer_2000_ddsd[648];

    auto g_yy_0_0_0_xx_xx_0_xy = buffer_2000_ddsd[649];

    auto g_yy_0_0_0_xx_xx_0_xz = buffer_2000_ddsd[650];

    auto g_yy_0_0_0_xx_xx_0_yy = buffer_2000_ddsd[651];

    auto g_yy_0_0_0_xx_xx_0_yz = buffer_2000_ddsd[652];

    auto g_yy_0_0_0_xx_xx_0_zz = buffer_2000_ddsd[653];

    auto g_yy_0_0_0_xx_xy_0_xx = buffer_2000_ddsd[654];

    auto g_yy_0_0_0_xx_xy_0_xy = buffer_2000_ddsd[655];

    auto g_yy_0_0_0_xx_xy_0_xz = buffer_2000_ddsd[656];

    auto g_yy_0_0_0_xx_xy_0_yy = buffer_2000_ddsd[657];

    auto g_yy_0_0_0_xx_xy_0_yz = buffer_2000_ddsd[658];

    auto g_yy_0_0_0_xx_xy_0_zz = buffer_2000_ddsd[659];

    auto g_yy_0_0_0_xx_xz_0_xx = buffer_2000_ddsd[660];

    auto g_yy_0_0_0_xx_xz_0_xy = buffer_2000_ddsd[661];

    auto g_yy_0_0_0_xx_xz_0_xz = buffer_2000_ddsd[662];

    auto g_yy_0_0_0_xx_xz_0_yy = buffer_2000_ddsd[663];

    auto g_yy_0_0_0_xx_xz_0_yz = buffer_2000_ddsd[664];

    auto g_yy_0_0_0_xx_xz_0_zz = buffer_2000_ddsd[665];

    auto g_yy_0_0_0_xx_yy_0_xx = buffer_2000_ddsd[666];

    auto g_yy_0_0_0_xx_yy_0_xy = buffer_2000_ddsd[667];

    auto g_yy_0_0_0_xx_yy_0_xz = buffer_2000_ddsd[668];

    auto g_yy_0_0_0_xx_yy_0_yy = buffer_2000_ddsd[669];

    auto g_yy_0_0_0_xx_yy_0_yz = buffer_2000_ddsd[670];

    auto g_yy_0_0_0_xx_yy_0_zz = buffer_2000_ddsd[671];

    auto g_yy_0_0_0_xx_yz_0_xx = buffer_2000_ddsd[672];

    auto g_yy_0_0_0_xx_yz_0_xy = buffer_2000_ddsd[673];

    auto g_yy_0_0_0_xx_yz_0_xz = buffer_2000_ddsd[674];

    auto g_yy_0_0_0_xx_yz_0_yy = buffer_2000_ddsd[675];

    auto g_yy_0_0_0_xx_yz_0_yz = buffer_2000_ddsd[676];

    auto g_yy_0_0_0_xx_yz_0_zz = buffer_2000_ddsd[677];

    auto g_yy_0_0_0_xx_zz_0_xx = buffer_2000_ddsd[678];

    auto g_yy_0_0_0_xx_zz_0_xy = buffer_2000_ddsd[679];

    auto g_yy_0_0_0_xx_zz_0_xz = buffer_2000_ddsd[680];

    auto g_yy_0_0_0_xx_zz_0_yy = buffer_2000_ddsd[681];

    auto g_yy_0_0_0_xx_zz_0_yz = buffer_2000_ddsd[682];

    auto g_yy_0_0_0_xx_zz_0_zz = buffer_2000_ddsd[683];

    auto g_yy_0_0_0_xy_xx_0_xx = buffer_2000_ddsd[684];

    auto g_yy_0_0_0_xy_xx_0_xy = buffer_2000_ddsd[685];

    auto g_yy_0_0_0_xy_xx_0_xz = buffer_2000_ddsd[686];

    auto g_yy_0_0_0_xy_xx_0_yy = buffer_2000_ddsd[687];

    auto g_yy_0_0_0_xy_xx_0_yz = buffer_2000_ddsd[688];

    auto g_yy_0_0_0_xy_xx_0_zz = buffer_2000_ddsd[689];

    auto g_yy_0_0_0_xy_xy_0_xx = buffer_2000_ddsd[690];

    auto g_yy_0_0_0_xy_xy_0_xy = buffer_2000_ddsd[691];

    auto g_yy_0_0_0_xy_xy_0_xz = buffer_2000_ddsd[692];

    auto g_yy_0_0_0_xy_xy_0_yy = buffer_2000_ddsd[693];

    auto g_yy_0_0_0_xy_xy_0_yz = buffer_2000_ddsd[694];

    auto g_yy_0_0_0_xy_xy_0_zz = buffer_2000_ddsd[695];

    auto g_yy_0_0_0_xy_xz_0_xx = buffer_2000_ddsd[696];

    auto g_yy_0_0_0_xy_xz_0_xy = buffer_2000_ddsd[697];

    auto g_yy_0_0_0_xy_xz_0_xz = buffer_2000_ddsd[698];

    auto g_yy_0_0_0_xy_xz_0_yy = buffer_2000_ddsd[699];

    auto g_yy_0_0_0_xy_xz_0_yz = buffer_2000_ddsd[700];

    auto g_yy_0_0_0_xy_xz_0_zz = buffer_2000_ddsd[701];

    auto g_yy_0_0_0_xy_yy_0_xx = buffer_2000_ddsd[702];

    auto g_yy_0_0_0_xy_yy_0_xy = buffer_2000_ddsd[703];

    auto g_yy_0_0_0_xy_yy_0_xz = buffer_2000_ddsd[704];

    auto g_yy_0_0_0_xy_yy_0_yy = buffer_2000_ddsd[705];

    auto g_yy_0_0_0_xy_yy_0_yz = buffer_2000_ddsd[706];

    auto g_yy_0_0_0_xy_yy_0_zz = buffer_2000_ddsd[707];

    auto g_yy_0_0_0_xy_yz_0_xx = buffer_2000_ddsd[708];

    auto g_yy_0_0_0_xy_yz_0_xy = buffer_2000_ddsd[709];

    auto g_yy_0_0_0_xy_yz_0_xz = buffer_2000_ddsd[710];

    auto g_yy_0_0_0_xy_yz_0_yy = buffer_2000_ddsd[711];

    auto g_yy_0_0_0_xy_yz_0_yz = buffer_2000_ddsd[712];

    auto g_yy_0_0_0_xy_yz_0_zz = buffer_2000_ddsd[713];

    auto g_yy_0_0_0_xy_zz_0_xx = buffer_2000_ddsd[714];

    auto g_yy_0_0_0_xy_zz_0_xy = buffer_2000_ddsd[715];

    auto g_yy_0_0_0_xy_zz_0_xz = buffer_2000_ddsd[716];

    auto g_yy_0_0_0_xy_zz_0_yy = buffer_2000_ddsd[717];

    auto g_yy_0_0_0_xy_zz_0_yz = buffer_2000_ddsd[718];

    auto g_yy_0_0_0_xy_zz_0_zz = buffer_2000_ddsd[719];

    auto g_yy_0_0_0_xz_xx_0_xx = buffer_2000_ddsd[720];

    auto g_yy_0_0_0_xz_xx_0_xy = buffer_2000_ddsd[721];

    auto g_yy_0_0_0_xz_xx_0_xz = buffer_2000_ddsd[722];

    auto g_yy_0_0_0_xz_xx_0_yy = buffer_2000_ddsd[723];

    auto g_yy_0_0_0_xz_xx_0_yz = buffer_2000_ddsd[724];

    auto g_yy_0_0_0_xz_xx_0_zz = buffer_2000_ddsd[725];

    auto g_yy_0_0_0_xz_xy_0_xx = buffer_2000_ddsd[726];

    auto g_yy_0_0_0_xz_xy_0_xy = buffer_2000_ddsd[727];

    auto g_yy_0_0_0_xz_xy_0_xz = buffer_2000_ddsd[728];

    auto g_yy_0_0_0_xz_xy_0_yy = buffer_2000_ddsd[729];

    auto g_yy_0_0_0_xz_xy_0_yz = buffer_2000_ddsd[730];

    auto g_yy_0_0_0_xz_xy_0_zz = buffer_2000_ddsd[731];

    auto g_yy_0_0_0_xz_xz_0_xx = buffer_2000_ddsd[732];

    auto g_yy_0_0_0_xz_xz_0_xy = buffer_2000_ddsd[733];

    auto g_yy_0_0_0_xz_xz_0_xz = buffer_2000_ddsd[734];

    auto g_yy_0_0_0_xz_xz_0_yy = buffer_2000_ddsd[735];

    auto g_yy_0_0_0_xz_xz_0_yz = buffer_2000_ddsd[736];

    auto g_yy_0_0_0_xz_xz_0_zz = buffer_2000_ddsd[737];

    auto g_yy_0_0_0_xz_yy_0_xx = buffer_2000_ddsd[738];

    auto g_yy_0_0_0_xz_yy_0_xy = buffer_2000_ddsd[739];

    auto g_yy_0_0_0_xz_yy_0_xz = buffer_2000_ddsd[740];

    auto g_yy_0_0_0_xz_yy_0_yy = buffer_2000_ddsd[741];

    auto g_yy_0_0_0_xz_yy_0_yz = buffer_2000_ddsd[742];

    auto g_yy_0_0_0_xz_yy_0_zz = buffer_2000_ddsd[743];

    auto g_yy_0_0_0_xz_yz_0_xx = buffer_2000_ddsd[744];

    auto g_yy_0_0_0_xz_yz_0_xy = buffer_2000_ddsd[745];

    auto g_yy_0_0_0_xz_yz_0_xz = buffer_2000_ddsd[746];

    auto g_yy_0_0_0_xz_yz_0_yy = buffer_2000_ddsd[747];

    auto g_yy_0_0_0_xz_yz_0_yz = buffer_2000_ddsd[748];

    auto g_yy_0_0_0_xz_yz_0_zz = buffer_2000_ddsd[749];

    auto g_yy_0_0_0_xz_zz_0_xx = buffer_2000_ddsd[750];

    auto g_yy_0_0_0_xz_zz_0_xy = buffer_2000_ddsd[751];

    auto g_yy_0_0_0_xz_zz_0_xz = buffer_2000_ddsd[752];

    auto g_yy_0_0_0_xz_zz_0_yy = buffer_2000_ddsd[753];

    auto g_yy_0_0_0_xz_zz_0_yz = buffer_2000_ddsd[754];

    auto g_yy_0_0_0_xz_zz_0_zz = buffer_2000_ddsd[755];

    auto g_yy_0_0_0_yy_xx_0_xx = buffer_2000_ddsd[756];

    auto g_yy_0_0_0_yy_xx_0_xy = buffer_2000_ddsd[757];

    auto g_yy_0_0_0_yy_xx_0_xz = buffer_2000_ddsd[758];

    auto g_yy_0_0_0_yy_xx_0_yy = buffer_2000_ddsd[759];

    auto g_yy_0_0_0_yy_xx_0_yz = buffer_2000_ddsd[760];

    auto g_yy_0_0_0_yy_xx_0_zz = buffer_2000_ddsd[761];

    auto g_yy_0_0_0_yy_xy_0_xx = buffer_2000_ddsd[762];

    auto g_yy_0_0_0_yy_xy_0_xy = buffer_2000_ddsd[763];

    auto g_yy_0_0_0_yy_xy_0_xz = buffer_2000_ddsd[764];

    auto g_yy_0_0_0_yy_xy_0_yy = buffer_2000_ddsd[765];

    auto g_yy_0_0_0_yy_xy_0_yz = buffer_2000_ddsd[766];

    auto g_yy_0_0_0_yy_xy_0_zz = buffer_2000_ddsd[767];

    auto g_yy_0_0_0_yy_xz_0_xx = buffer_2000_ddsd[768];

    auto g_yy_0_0_0_yy_xz_0_xy = buffer_2000_ddsd[769];

    auto g_yy_0_0_0_yy_xz_0_xz = buffer_2000_ddsd[770];

    auto g_yy_0_0_0_yy_xz_0_yy = buffer_2000_ddsd[771];

    auto g_yy_0_0_0_yy_xz_0_yz = buffer_2000_ddsd[772];

    auto g_yy_0_0_0_yy_xz_0_zz = buffer_2000_ddsd[773];

    auto g_yy_0_0_0_yy_yy_0_xx = buffer_2000_ddsd[774];

    auto g_yy_0_0_0_yy_yy_0_xy = buffer_2000_ddsd[775];

    auto g_yy_0_0_0_yy_yy_0_xz = buffer_2000_ddsd[776];

    auto g_yy_0_0_0_yy_yy_0_yy = buffer_2000_ddsd[777];

    auto g_yy_0_0_0_yy_yy_0_yz = buffer_2000_ddsd[778];

    auto g_yy_0_0_0_yy_yy_0_zz = buffer_2000_ddsd[779];

    auto g_yy_0_0_0_yy_yz_0_xx = buffer_2000_ddsd[780];

    auto g_yy_0_0_0_yy_yz_0_xy = buffer_2000_ddsd[781];

    auto g_yy_0_0_0_yy_yz_0_xz = buffer_2000_ddsd[782];

    auto g_yy_0_0_0_yy_yz_0_yy = buffer_2000_ddsd[783];

    auto g_yy_0_0_0_yy_yz_0_yz = buffer_2000_ddsd[784];

    auto g_yy_0_0_0_yy_yz_0_zz = buffer_2000_ddsd[785];

    auto g_yy_0_0_0_yy_zz_0_xx = buffer_2000_ddsd[786];

    auto g_yy_0_0_0_yy_zz_0_xy = buffer_2000_ddsd[787];

    auto g_yy_0_0_0_yy_zz_0_xz = buffer_2000_ddsd[788];

    auto g_yy_0_0_0_yy_zz_0_yy = buffer_2000_ddsd[789];

    auto g_yy_0_0_0_yy_zz_0_yz = buffer_2000_ddsd[790];

    auto g_yy_0_0_0_yy_zz_0_zz = buffer_2000_ddsd[791];

    auto g_yy_0_0_0_yz_xx_0_xx = buffer_2000_ddsd[792];

    auto g_yy_0_0_0_yz_xx_0_xy = buffer_2000_ddsd[793];

    auto g_yy_0_0_0_yz_xx_0_xz = buffer_2000_ddsd[794];

    auto g_yy_0_0_0_yz_xx_0_yy = buffer_2000_ddsd[795];

    auto g_yy_0_0_0_yz_xx_0_yz = buffer_2000_ddsd[796];

    auto g_yy_0_0_0_yz_xx_0_zz = buffer_2000_ddsd[797];

    auto g_yy_0_0_0_yz_xy_0_xx = buffer_2000_ddsd[798];

    auto g_yy_0_0_0_yz_xy_0_xy = buffer_2000_ddsd[799];

    auto g_yy_0_0_0_yz_xy_0_xz = buffer_2000_ddsd[800];

    auto g_yy_0_0_0_yz_xy_0_yy = buffer_2000_ddsd[801];

    auto g_yy_0_0_0_yz_xy_0_yz = buffer_2000_ddsd[802];

    auto g_yy_0_0_0_yz_xy_0_zz = buffer_2000_ddsd[803];

    auto g_yy_0_0_0_yz_xz_0_xx = buffer_2000_ddsd[804];

    auto g_yy_0_0_0_yz_xz_0_xy = buffer_2000_ddsd[805];

    auto g_yy_0_0_0_yz_xz_0_xz = buffer_2000_ddsd[806];

    auto g_yy_0_0_0_yz_xz_0_yy = buffer_2000_ddsd[807];

    auto g_yy_0_0_0_yz_xz_0_yz = buffer_2000_ddsd[808];

    auto g_yy_0_0_0_yz_xz_0_zz = buffer_2000_ddsd[809];

    auto g_yy_0_0_0_yz_yy_0_xx = buffer_2000_ddsd[810];

    auto g_yy_0_0_0_yz_yy_0_xy = buffer_2000_ddsd[811];

    auto g_yy_0_0_0_yz_yy_0_xz = buffer_2000_ddsd[812];

    auto g_yy_0_0_0_yz_yy_0_yy = buffer_2000_ddsd[813];

    auto g_yy_0_0_0_yz_yy_0_yz = buffer_2000_ddsd[814];

    auto g_yy_0_0_0_yz_yy_0_zz = buffer_2000_ddsd[815];

    auto g_yy_0_0_0_yz_yz_0_xx = buffer_2000_ddsd[816];

    auto g_yy_0_0_0_yz_yz_0_xy = buffer_2000_ddsd[817];

    auto g_yy_0_0_0_yz_yz_0_xz = buffer_2000_ddsd[818];

    auto g_yy_0_0_0_yz_yz_0_yy = buffer_2000_ddsd[819];

    auto g_yy_0_0_0_yz_yz_0_yz = buffer_2000_ddsd[820];

    auto g_yy_0_0_0_yz_yz_0_zz = buffer_2000_ddsd[821];

    auto g_yy_0_0_0_yz_zz_0_xx = buffer_2000_ddsd[822];

    auto g_yy_0_0_0_yz_zz_0_xy = buffer_2000_ddsd[823];

    auto g_yy_0_0_0_yz_zz_0_xz = buffer_2000_ddsd[824];

    auto g_yy_0_0_0_yz_zz_0_yy = buffer_2000_ddsd[825];

    auto g_yy_0_0_0_yz_zz_0_yz = buffer_2000_ddsd[826];

    auto g_yy_0_0_0_yz_zz_0_zz = buffer_2000_ddsd[827];

    auto g_yy_0_0_0_zz_xx_0_xx = buffer_2000_ddsd[828];

    auto g_yy_0_0_0_zz_xx_0_xy = buffer_2000_ddsd[829];

    auto g_yy_0_0_0_zz_xx_0_xz = buffer_2000_ddsd[830];

    auto g_yy_0_0_0_zz_xx_0_yy = buffer_2000_ddsd[831];

    auto g_yy_0_0_0_zz_xx_0_yz = buffer_2000_ddsd[832];

    auto g_yy_0_0_0_zz_xx_0_zz = buffer_2000_ddsd[833];

    auto g_yy_0_0_0_zz_xy_0_xx = buffer_2000_ddsd[834];

    auto g_yy_0_0_0_zz_xy_0_xy = buffer_2000_ddsd[835];

    auto g_yy_0_0_0_zz_xy_0_xz = buffer_2000_ddsd[836];

    auto g_yy_0_0_0_zz_xy_0_yy = buffer_2000_ddsd[837];

    auto g_yy_0_0_0_zz_xy_0_yz = buffer_2000_ddsd[838];

    auto g_yy_0_0_0_zz_xy_0_zz = buffer_2000_ddsd[839];

    auto g_yy_0_0_0_zz_xz_0_xx = buffer_2000_ddsd[840];

    auto g_yy_0_0_0_zz_xz_0_xy = buffer_2000_ddsd[841];

    auto g_yy_0_0_0_zz_xz_0_xz = buffer_2000_ddsd[842];

    auto g_yy_0_0_0_zz_xz_0_yy = buffer_2000_ddsd[843];

    auto g_yy_0_0_0_zz_xz_0_yz = buffer_2000_ddsd[844];

    auto g_yy_0_0_0_zz_xz_0_zz = buffer_2000_ddsd[845];

    auto g_yy_0_0_0_zz_yy_0_xx = buffer_2000_ddsd[846];

    auto g_yy_0_0_0_zz_yy_0_xy = buffer_2000_ddsd[847];

    auto g_yy_0_0_0_zz_yy_0_xz = buffer_2000_ddsd[848];

    auto g_yy_0_0_0_zz_yy_0_yy = buffer_2000_ddsd[849];

    auto g_yy_0_0_0_zz_yy_0_yz = buffer_2000_ddsd[850];

    auto g_yy_0_0_0_zz_yy_0_zz = buffer_2000_ddsd[851];

    auto g_yy_0_0_0_zz_yz_0_xx = buffer_2000_ddsd[852];

    auto g_yy_0_0_0_zz_yz_0_xy = buffer_2000_ddsd[853];

    auto g_yy_0_0_0_zz_yz_0_xz = buffer_2000_ddsd[854];

    auto g_yy_0_0_0_zz_yz_0_yy = buffer_2000_ddsd[855];

    auto g_yy_0_0_0_zz_yz_0_yz = buffer_2000_ddsd[856];

    auto g_yy_0_0_0_zz_yz_0_zz = buffer_2000_ddsd[857];

    auto g_yy_0_0_0_zz_zz_0_xx = buffer_2000_ddsd[858];

    auto g_yy_0_0_0_zz_zz_0_xy = buffer_2000_ddsd[859];

    auto g_yy_0_0_0_zz_zz_0_xz = buffer_2000_ddsd[860];

    auto g_yy_0_0_0_zz_zz_0_yy = buffer_2000_ddsd[861];

    auto g_yy_0_0_0_zz_zz_0_yz = buffer_2000_ddsd[862];

    auto g_yy_0_0_0_zz_zz_0_zz = buffer_2000_ddsd[863];

    auto g_yz_0_0_0_xx_xx_0_xx = buffer_2000_ddsd[864];

    auto g_yz_0_0_0_xx_xx_0_xy = buffer_2000_ddsd[865];

    auto g_yz_0_0_0_xx_xx_0_xz = buffer_2000_ddsd[866];

    auto g_yz_0_0_0_xx_xx_0_yy = buffer_2000_ddsd[867];

    auto g_yz_0_0_0_xx_xx_0_yz = buffer_2000_ddsd[868];

    auto g_yz_0_0_0_xx_xx_0_zz = buffer_2000_ddsd[869];

    auto g_yz_0_0_0_xx_xy_0_xx = buffer_2000_ddsd[870];

    auto g_yz_0_0_0_xx_xy_0_xy = buffer_2000_ddsd[871];

    auto g_yz_0_0_0_xx_xy_0_xz = buffer_2000_ddsd[872];

    auto g_yz_0_0_0_xx_xy_0_yy = buffer_2000_ddsd[873];

    auto g_yz_0_0_0_xx_xy_0_yz = buffer_2000_ddsd[874];

    auto g_yz_0_0_0_xx_xy_0_zz = buffer_2000_ddsd[875];

    auto g_yz_0_0_0_xx_xz_0_xx = buffer_2000_ddsd[876];

    auto g_yz_0_0_0_xx_xz_0_xy = buffer_2000_ddsd[877];

    auto g_yz_0_0_0_xx_xz_0_xz = buffer_2000_ddsd[878];

    auto g_yz_0_0_0_xx_xz_0_yy = buffer_2000_ddsd[879];

    auto g_yz_0_0_0_xx_xz_0_yz = buffer_2000_ddsd[880];

    auto g_yz_0_0_0_xx_xz_0_zz = buffer_2000_ddsd[881];

    auto g_yz_0_0_0_xx_yy_0_xx = buffer_2000_ddsd[882];

    auto g_yz_0_0_0_xx_yy_0_xy = buffer_2000_ddsd[883];

    auto g_yz_0_0_0_xx_yy_0_xz = buffer_2000_ddsd[884];

    auto g_yz_0_0_0_xx_yy_0_yy = buffer_2000_ddsd[885];

    auto g_yz_0_0_0_xx_yy_0_yz = buffer_2000_ddsd[886];

    auto g_yz_0_0_0_xx_yy_0_zz = buffer_2000_ddsd[887];

    auto g_yz_0_0_0_xx_yz_0_xx = buffer_2000_ddsd[888];

    auto g_yz_0_0_0_xx_yz_0_xy = buffer_2000_ddsd[889];

    auto g_yz_0_0_0_xx_yz_0_xz = buffer_2000_ddsd[890];

    auto g_yz_0_0_0_xx_yz_0_yy = buffer_2000_ddsd[891];

    auto g_yz_0_0_0_xx_yz_0_yz = buffer_2000_ddsd[892];

    auto g_yz_0_0_0_xx_yz_0_zz = buffer_2000_ddsd[893];

    auto g_yz_0_0_0_xx_zz_0_xx = buffer_2000_ddsd[894];

    auto g_yz_0_0_0_xx_zz_0_xy = buffer_2000_ddsd[895];

    auto g_yz_0_0_0_xx_zz_0_xz = buffer_2000_ddsd[896];

    auto g_yz_0_0_0_xx_zz_0_yy = buffer_2000_ddsd[897];

    auto g_yz_0_0_0_xx_zz_0_yz = buffer_2000_ddsd[898];

    auto g_yz_0_0_0_xx_zz_0_zz = buffer_2000_ddsd[899];

    auto g_yz_0_0_0_xy_xx_0_xx = buffer_2000_ddsd[900];

    auto g_yz_0_0_0_xy_xx_0_xy = buffer_2000_ddsd[901];

    auto g_yz_0_0_0_xy_xx_0_xz = buffer_2000_ddsd[902];

    auto g_yz_0_0_0_xy_xx_0_yy = buffer_2000_ddsd[903];

    auto g_yz_0_0_0_xy_xx_0_yz = buffer_2000_ddsd[904];

    auto g_yz_0_0_0_xy_xx_0_zz = buffer_2000_ddsd[905];

    auto g_yz_0_0_0_xy_xy_0_xx = buffer_2000_ddsd[906];

    auto g_yz_0_0_0_xy_xy_0_xy = buffer_2000_ddsd[907];

    auto g_yz_0_0_0_xy_xy_0_xz = buffer_2000_ddsd[908];

    auto g_yz_0_0_0_xy_xy_0_yy = buffer_2000_ddsd[909];

    auto g_yz_0_0_0_xy_xy_0_yz = buffer_2000_ddsd[910];

    auto g_yz_0_0_0_xy_xy_0_zz = buffer_2000_ddsd[911];

    auto g_yz_0_0_0_xy_xz_0_xx = buffer_2000_ddsd[912];

    auto g_yz_0_0_0_xy_xz_0_xy = buffer_2000_ddsd[913];

    auto g_yz_0_0_0_xy_xz_0_xz = buffer_2000_ddsd[914];

    auto g_yz_0_0_0_xy_xz_0_yy = buffer_2000_ddsd[915];

    auto g_yz_0_0_0_xy_xz_0_yz = buffer_2000_ddsd[916];

    auto g_yz_0_0_0_xy_xz_0_zz = buffer_2000_ddsd[917];

    auto g_yz_0_0_0_xy_yy_0_xx = buffer_2000_ddsd[918];

    auto g_yz_0_0_0_xy_yy_0_xy = buffer_2000_ddsd[919];

    auto g_yz_0_0_0_xy_yy_0_xz = buffer_2000_ddsd[920];

    auto g_yz_0_0_0_xy_yy_0_yy = buffer_2000_ddsd[921];

    auto g_yz_0_0_0_xy_yy_0_yz = buffer_2000_ddsd[922];

    auto g_yz_0_0_0_xy_yy_0_zz = buffer_2000_ddsd[923];

    auto g_yz_0_0_0_xy_yz_0_xx = buffer_2000_ddsd[924];

    auto g_yz_0_0_0_xy_yz_0_xy = buffer_2000_ddsd[925];

    auto g_yz_0_0_0_xy_yz_0_xz = buffer_2000_ddsd[926];

    auto g_yz_0_0_0_xy_yz_0_yy = buffer_2000_ddsd[927];

    auto g_yz_0_0_0_xy_yz_0_yz = buffer_2000_ddsd[928];

    auto g_yz_0_0_0_xy_yz_0_zz = buffer_2000_ddsd[929];

    auto g_yz_0_0_0_xy_zz_0_xx = buffer_2000_ddsd[930];

    auto g_yz_0_0_0_xy_zz_0_xy = buffer_2000_ddsd[931];

    auto g_yz_0_0_0_xy_zz_0_xz = buffer_2000_ddsd[932];

    auto g_yz_0_0_0_xy_zz_0_yy = buffer_2000_ddsd[933];

    auto g_yz_0_0_0_xy_zz_0_yz = buffer_2000_ddsd[934];

    auto g_yz_0_0_0_xy_zz_0_zz = buffer_2000_ddsd[935];

    auto g_yz_0_0_0_xz_xx_0_xx = buffer_2000_ddsd[936];

    auto g_yz_0_0_0_xz_xx_0_xy = buffer_2000_ddsd[937];

    auto g_yz_0_0_0_xz_xx_0_xz = buffer_2000_ddsd[938];

    auto g_yz_0_0_0_xz_xx_0_yy = buffer_2000_ddsd[939];

    auto g_yz_0_0_0_xz_xx_0_yz = buffer_2000_ddsd[940];

    auto g_yz_0_0_0_xz_xx_0_zz = buffer_2000_ddsd[941];

    auto g_yz_0_0_0_xz_xy_0_xx = buffer_2000_ddsd[942];

    auto g_yz_0_0_0_xz_xy_0_xy = buffer_2000_ddsd[943];

    auto g_yz_0_0_0_xz_xy_0_xz = buffer_2000_ddsd[944];

    auto g_yz_0_0_0_xz_xy_0_yy = buffer_2000_ddsd[945];

    auto g_yz_0_0_0_xz_xy_0_yz = buffer_2000_ddsd[946];

    auto g_yz_0_0_0_xz_xy_0_zz = buffer_2000_ddsd[947];

    auto g_yz_0_0_0_xz_xz_0_xx = buffer_2000_ddsd[948];

    auto g_yz_0_0_0_xz_xz_0_xy = buffer_2000_ddsd[949];

    auto g_yz_0_0_0_xz_xz_0_xz = buffer_2000_ddsd[950];

    auto g_yz_0_0_0_xz_xz_0_yy = buffer_2000_ddsd[951];

    auto g_yz_0_0_0_xz_xz_0_yz = buffer_2000_ddsd[952];

    auto g_yz_0_0_0_xz_xz_0_zz = buffer_2000_ddsd[953];

    auto g_yz_0_0_0_xz_yy_0_xx = buffer_2000_ddsd[954];

    auto g_yz_0_0_0_xz_yy_0_xy = buffer_2000_ddsd[955];

    auto g_yz_0_0_0_xz_yy_0_xz = buffer_2000_ddsd[956];

    auto g_yz_0_0_0_xz_yy_0_yy = buffer_2000_ddsd[957];

    auto g_yz_0_0_0_xz_yy_0_yz = buffer_2000_ddsd[958];

    auto g_yz_0_0_0_xz_yy_0_zz = buffer_2000_ddsd[959];

    auto g_yz_0_0_0_xz_yz_0_xx = buffer_2000_ddsd[960];

    auto g_yz_0_0_0_xz_yz_0_xy = buffer_2000_ddsd[961];

    auto g_yz_0_0_0_xz_yz_0_xz = buffer_2000_ddsd[962];

    auto g_yz_0_0_0_xz_yz_0_yy = buffer_2000_ddsd[963];

    auto g_yz_0_0_0_xz_yz_0_yz = buffer_2000_ddsd[964];

    auto g_yz_0_0_0_xz_yz_0_zz = buffer_2000_ddsd[965];

    auto g_yz_0_0_0_xz_zz_0_xx = buffer_2000_ddsd[966];

    auto g_yz_0_0_0_xz_zz_0_xy = buffer_2000_ddsd[967];

    auto g_yz_0_0_0_xz_zz_0_xz = buffer_2000_ddsd[968];

    auto g_yz_0_0_0_xz_zz_0_yy = buffer_2000_ddsd[969];

    auto g_yz_0_0_0_xz_zz_0_yz = buffer_2000_ddsd[970];

    auto g_yz_0_0_0_xz_zz_0_zz = buffer_2000_ddsd[971];

    auto g_yz_0_0_0_yy_xx_0_xx = buffer_2000_ddsd[972];

    auto g_yz_0_0_0_yy_xx_0_xy = buffer_2000_ddsd[973];

    auto g_yz_0_0_0_yy_xx_0_xz = buffer_2000_ddsd[974];

    auto g_yz_0_0_0_yy_xx_0_yy = buffer_2000_ddsd[975];

    auto g_yz_0_0_0_yy_xx_0_yz = buffer_2000_ddsd[976];

    auto g_yz_0_0_0_yy_xx_0_zz = buffer_2000_ddsd[977];

    auto g_yz_0_0_0_yy_xy_0_xx = buffer_2000_ddsd[978];

    auto g_yz_0_0_0_yy_xy_0_xy = buffer_2000_ddsd[979];

    auto g_yz_0_0_0_yy_xy_0_xz = buffer_2000_ddsd[980];

    auto g_yz_0_0_0_yy_xy_0_yy = buffer_2000_ddsd[981];

    auto g_yz_0_0_0_yy_xy_0_yz = buffer_2000_ddsd[982];

    auto g_yz_0_0_0_yy_xy_0_zz = buffer_2000_ddsd[983];

    auto g_yz_0_0_0_yy_xz_0_xx = buffer_2000_ddsd[984];

    auto g_yz_0_0_0_yy_xz_0_xy = buffer_2000_ddsd[985];

    auto g_yz_0_0_0_yy_xz_0_xz = buffer_2000_ddsd[986];

    auto g_yz_0_0_0_yy_xz_0_yy = buffer_2000_ddsd[987];

    auto g_yz_0_0_0_yy_xz_0_yz = buffer_2000_ddsd[988];

    auto g_yz_0_0_0_yy_xz_0_zz = buffer_2000_ddsd[989];

    auto g_yz_0_0_0_yy_yy_0_xx = buffer_2000_ddsd[990];

    auto g_yz_0_0_0_yy_yy_0_xy = buffer_2000_ddsd[991];

    auto g_yz_0_0_0_yy_yy_0_xz = buffer_2000_ddsd[992];

    auto g_yz_0_0_0_yy_yy_0_yy = buffer_2000_ddsd[993];

    auto g_yz_0_0_0_yy_yy_0_yz = buffer_2000_ddsd[994];

    auto g_yz_0_0_0_yy_yy_0_zz = buffer_2000_ddsd[995];

    auto g_yz_0_0_0_yy_yz_0_xx = buffer_2000_ddsd[996];

    auto g_yz_0_0_0_yy_yz_0_xy = buffer_2000_ddsd[997];

    auto g_yz_0_0_0_yy_yz_0_xz = buffer_2000_ddsd[998];

    auto g_yz_0_0_0_yy_yz_0_yy = buffer_2000_ddsd[999];

    auto g_yz_0_0_0_yy_yz_0_yz = buffer_2000_ddsd[1000];

    auto g_yz_0_0_0_yy_yz_0_zz = buffer_2000_ddsd[1001];

    auto g_yz_0_0_0_yy_zz_0_xx = buffer_2000_ddsd[1002];

    auto g_yz_0_0_0_yy_zz_0_xy = buffer_2000_ddsd[1003];

    auto g_yz_0_0_0_yy_zz_0_xz = buffer_2000_ddsd[1004];

    auto g_yz_0_0_0_yy_zz_0_yy = buffer_2000_ddsd[1005];

    auto g_yz_0_0_0_yy_zz_0_yz = buffer_2000_ddsd[1006];

    auto g_yz_0_0_0_yy_zz_0_zz = buffer_2000_ddsd[1007];

    auto g_yz_0_0_0_yz_xx_0_xx = buffer_2000_ddsd[1008];

    auto g_yz_0_0_0_yz_xx_0_xy = buffer_2000_ddsd[1009];

    auto g_yz_0_0_0_yz_xx_0_xz = buffer_2000_ddsd[1010];

    auto g_yz_0_0_0_yz_xx_0_yy = buffer_2000_ddsd[1011];

    auto g_yz_0_0_0_yz_xx_0_yz = buffer_2000_ddsd[1012];

    auto g_yz_0_0_0_yz_xx_0_zz = buffer_2000_ddsd[1013];

    auto g_yz_0_0_0_yz_xy_0_xx = buffer_2000_ddsd[1014];

    auto g_yz_0_0_0_yz_xy_0_xy = buffer_2000_ddsd[1015];

    auto g_yz_0_0_0_yz_xy_0_xz = buffer_2000_ddsd[1016];

    auto g_yz_0_0_0_yz_xy_0_yy = buffer_2000_ddsd[1017];

    auto g_yz_0_0_0_yz_xy_0_yz = buffer_2000_ddsd[1018];

    auto g_yz_0_0_0_yz_xy_0_zz = buffer_2000_ddsd[1019];

    auto g_yz_0_0_0_yz_xz_0_xx = buffer_2000_ddsd[1020];

    auto g_yz_0_0_0_yz_xz_0_xy = buffer_2000_ddsd[1021];

    auto g_yz_0_0_0_yz_xz_0_xz = buffer_2000_ddsd[1022];

    auto g_yz_0_0_0_yz_xz_0_yy = buffer_2000_ddsd[1023];

    auto g_yz_0_0_0_yz_xz_0_yz = buffer_2000_ddsd[1024];

    auto g_yz_0_0_0_yz_xz_0_zz = buffer_2000_ddsd[1025];

    auto g_yz_0_0_0_yz_yy_0_xx = buffer_2000_ddsd[1026];

    auto g_yz_0_0_0_yz_yy_0_xy = buffer_2000_ddsd[1027];

    auto g_yz_0_0_0_yz_yy_0_xz = buffer_2000_ddsd[1028];

    auto g_yz_0_0_0_yz_yy_0_yy = buffer_2000_ddsd[1029];

    auto g_yz_0_0_0_yz_yy_0_yz = buffer_2000_ddsd[1030];

    auto g_yz_0_0_0_yz_yy_0_zz = buffer_2000_ddsd[1031];

    auto g_yz_0_0_0_yz_yz_0_xx = buffer_2000_ddsd[1032];

    auto g_yz_0_0_0_yz_yz_0_xy = buffer_2000_ddsd[1033];

    auto g_yz_0_0_0_yz_yz_0_xz = buffer_2000_ddsd[1034];

    auto g_yz_0_0_0_yz_yz_0_yy = buffer_2000_ddsd[1035];

    auto g_yz_0_0_0_yz_yz_0_yz = buffer_2000_ddsd[1036];

    auto g_yz_0_0_0_yz_yz_0_zz = buffer_2000_ddsd[1037];

    auto g_yz_0_0_0_yz_zz_0_xx = buffer_2000_ddsd[1038];

    auto g_yz_0_0_0_yz_zz_0_xy = buffer_2000_ddsd[1039];

    auto g_yz_0_0_0_yz_zz_0_xz = buffer_2000_ddsd[1040];

    auto g_yz_0_0_0_yz_zz_0_yy = buffer_2000_ddsd[1041];

    auto g_yz_0_0_0_yz_zz_0_yz = buffer_2000_ddsd[1042];

    auto g_yz_0_0_0_yz_zz_0_zz = buffer_2000_ddsd[1043];

    auto g_yz_0_0_0_zz_xx_0_xx = buffer_2000_ddsd[1044];

    auto g_yz_0_0_0_zz_xx_0_xy = buffer_2000_ddsd[1045];

    auto g_yz_0_0_0_zz_xx_0_xz = buffer_2000_ddsd[1046];

    auto g_yz_0_0_0_zz_xx_0_yy = buffer_2000_ddsd[1047];

    auto g_yz_0_0_0_zz_xx_0_yz = buffer_2000_ddsd[1048];

    auto g_yz_0_0_0_zz_xx_0_zz = buffer_2000_ddsd[1049];

    auto g_yz_0_0_0_zz_xy_0_xx = buffer_2000_ddsd[1050];

    auto g_yz_0_0_0_zz_xy_0_xy = buffer_2000_ddsd[1051];

    auto g_yz_0_0_0_zz_xy_0_xz = buffer_2000_ddsd[1052];

    auto g_yz_0_0_0_zz_xy_0_yy = buffer_2000_ddsd[1053];

    auto g_yz_0_0_0_zz_xy_0_yz = buffer_2000_ddsd[1054];

    auto g_yz_0_0_0_zz_xy_0_zz = buffer_2000_ddsd[1055];

    auto g_yz_0_0_0_zz_xz_0_xx = buffer_2000_ddsd[1056];

    auto g_yz_0_0_0_zz_xz_0_xy = buffer_2000_ddsd[1057];

    auto g_yz_0_0_0_zz_xz_0_xz = buffer_2000_ddsd[1058];

    auto g_yz_0_0_0_zz_xz_0_yy = buffer_2000_ddsd[1059];

    auto g_yz_0_0_0_zz_xz_0_yz = buffer_2000_ddsd[1060];

    auto g_yz_0_0_0_zz_xz_0_zz = buffer_2000_ddsd[1061];

    auto g_yz_0_0_0_zz_yy_0_xx = buffer_2000_ddsd[1062];

    auto g_yz_0_0_0_zz_yy_0_xy = buffer_2000_ddsd[1063];

    auto g_yz_0_0_0_zz_yy_0_xz = buffer_2000_ddsd[1064];

    auto g_yz_0_0_0_zz_yy_0_yy = buffer_2000_ddsd[1065];

    auto g_yz_0_0_0_zz_yy_0_yz = buffer_2000_ddsd[1066];

    auto g_yz_0_0_0_zz_yy_0_zz = buffer_2000_ddsd[1067];

    auto g_yz_0_0_0_zz_yz_0_xx = buffer_2000_ddsd[1068];

    auto g_yz_0_0_0_zz_yz_0_xy = buffer_2000_ddsd[1069];

    auto g_yz_0_0_0_zz_yz_0_xz = buffer_2000_ddsd[1070];

    auto g_yz_0_0_0_zz_yz_0_yy = buffer_2000_ddsd[1071];

    auto g_yz_0_0_0_zz_yz_0_yz = buffer_2000_ddsd[1072];

    auto g_yz_0_0_0_zz_yz_0_zz = buffer_2000_ddsd[1073];

    auto g_yz_0_0_0_zz_zz_0_xx = buffer_2000_ddsd[1074];

    auto g_yz_0_0_0_zz_zz_0_xy = buffer_2000_ddsd[1075];

    auto g_yz_0_0_0_zz_zz_0_xz = buffer_2000_ddsd[1076];

    auto g_yz_0_0_0_zz_zz_0_yy = buffer_2000_ddsd[1077];

    auto g_yz_0_0_0_zz_zz_0_yz = buffer_2000_ddsd[1078];

    auto g_yz_0_0_0_zz_zz_0_zz = buffer_2000_ddsd[1079];

    auto g_zz_0_0_0_xx_xx_0_xx = buffer_2000_ddsd[1080];

    auto g_zz_0_0_0_xx_xx_0_xy = buffer_2000_ddsd[1081];

    auto g_zz_0_0_0_xx_xx_0_xz = buffer_2000_ddsd[1082];

    auto g_zz_0_0_0_xx_xx_0_yy = buffer_2000_ddsd[1083];

    auto g_zz_0_0_0_xx_xx_0_yz = buffer_2000_ddsd[1084];

    auto g_zz_0_0_0_xx_xx_0_zz = buffer_2000_ddsd[1085];

    auto g_zz_0_0_0_xx_xy_0_xx = buffer_2000_ddsd[1086];

    auto g_zz_0_0_0_xx_xy_0_xy = buffer_2000_ddsd[1087];

    auto g_zz_0_0_0_xx_xy_0_xz = buffer_2000_ddsd[1088];

    auto g_zz_0_0_0_xx_xy_0_yy = buffer_2000_ddsd[1089];

    auto g_zz_0_0_0_xx_xy_0_yz = buffer_2000_ddsd[1090];

    auto g_zz_0_0_0_xx_xy_0_zz = buffer_2000_ddsd[1091];

    auto g_zz_0_0_0_xx_xz_0_xx = buffer_2000_ddsd[1092];

    auto g_zz_0_0_0_xx_xz_0_xy = buffer_2000_ddsd[1093];

    auto g_zz_0_0_0_xx_xz_0_xz = buffer_2000_ddsd[1094];

    auto g_zz_0_0_0_xx_xz_0_yy = buffer_2000_ddsd[1095];

    auto g_zz_0_0_0_xx_xz_0_yz = buffer_2000_ddsd[1096];

    auto g_zz_0_0_0_xx_xz_0_zz = buffer_2000_ddsd[1097];

    auto g_zz_0_0_0_xx_yy_0_xx = buffer_2000_ddsd[1098];

    auto g_zz_0_0_0_xx_yy_0_xy = buffer_2000_ddsd[1099];

    auto g_zz_0_0_0_xx_yy_0_xz = buffer_2000_ddsd[1100];

    auto g_zz_0_0_0_xx_yy_0_yy = buffer_2000_ddsd[1101];

    auto g_zz_0_0_0_xx_yy_0_yz = buffer_2000_ddsd[1102];

    auto g_zz_0_0_0_xx_yy_0_zz = buffer_2000_ddsd[1103];

    auto g_zz_0_0_0_xx_yz_0_xx = buffer_2000_ddsd[1104];

    auto g_zz_0_0_0_xx_yz_0_xy = buffer_2000_ddsd[1105];

    auto g_zz_0_0_0_xx_yz_0_xz = buffer_2000_ddsd[1106];

    auto g_zz_0_0_0_xx_yz_0_yy = buffer_2000_ddsd[1107];

    auto g_zz_0_0_0_xx_yz_0_yz = buffer_2000_ddsd[1108];

    auto g_zz_0_0_0_xx_yz_0_zz = buffer_2000_ddsd[1109];

    auto g_zz_0_0_0_xx_zz_0_xx = buffer_2000_ddsd[1110];

    auto g_zz_0_0_0_xx_zz_0_xy = buffer_2000_ddsd[1111];

    auto g_zz_0_0_0_xx_zz_0_xz = buffer_2000_ddsd[1112];

    auto g_zz_0_0_0_xx_zz_0_yy = buffer_2000_ddsd[1113];

    auto g_zz_0_0_0_xx_zz_0_yz = buffer_2000_ddsd[1114];

    auto g_zz_0_0_0_xx_zz_0_zz = buffer_2000_ddsd[1115];

    auto g_zz_0_0_0_xy_xx_0_xx = buffer_2000_ddsd[1116];

    auto g_zz_0_0_0_xy_xx_0_xy = buffer_2000_ddsd[1117];

    auto g_zz_0_0_0_xy_xx_0_xz = buffer_2000_ddsd[1118];

    auto g_zz_0_0_0_xy_xx_0_yy = buffer_2000_ddsd[1119];

    auto g_zz_0_0_0_xy_xx_0_yz = buffer_2000_ddsd[1120];

    auto g_zz_0_0_0_xy_xx_0_zz = buffer_2000_ddsd[1121];

    auto g_zz_0_0_0_xy_xy_0_xx = buffer_2000_ddsd[1122];

    auto g_zz_0_0_0_xy_xy_0_xy = buffer_2000_ddsd[1123];

    auto g_zz_0_0_0_xy_xy_0_xz = buffer_2000_ddsd[1124];

    auto g_zz_0_0_0_xy_xy_0_yy = buffer_2000_ddsd[1125];

    auto g_zz_0_0_0_xy_xy_0_yz = buffer_2000_ddsd[1126];

    auto g_zz_0_0_0_xy_xy_0_zz = buffer_2000_ddsd[1127];

    auto g_zz_0_0_0_xy_xz_0_xx = buffer_2000_ddsd[1128];

    auto g_zz_0_0_0_xy_xz_0_xy = buffer_2000_ddsd[1129];

    auto g_zz_0_0_0_xy_xz_0_xz = buffer_2000_ddsd[1130];

    auto g_zz_0_0_0_xy_xz_0_yy = buffer_2000_ddsd[1131];

    auto g_zz_0_0_0_xy_xz_0_yz = buffer_2000_ddsd[1132];

    auto g_zz_0_0_0_xy_xz_0_zz = buffer_2000_ddsd[1133];

    auto g_zz_0_0_0_xy_yy_0_xx = buffer_2000_ddsd[1134];

    auto g_zz_0_0_0_xy_yy_0_xy = buffer_2000_ddsd[1135];

    auto g_zz_0_0_0_xy_yy_0_xz = buffer_2000_ddsd[1136];

    auto g_zz_0_0_0_xy_yy_0_yy = buffer_2000_ddsd[1137];

    auto g_zz_0_0_0_xy_yy_0_yz = buffer_2000_ddsd[1138];

    auto g_zz_0_0_0_xy_yy_0_zz = buffer_2000_ddsd[1139];

    auto g_zz_0_0_0_xy_yz_0_xx = buffer_2000_ddsd[1140];

    auto g_zz_0_0_0_xy_yz_0_xy = buffer_2000_ddsd[1141];

    auto g_zz_0_0_0_xy_yz_0_xz = buffer_2000_ddsd[1142];

    auto g_zz_0_0_0_xy_yz_0_yy = buffer_2000_ddsd[1143];

    auto g_zz_0_0_0_xy_yz_0_yz = buffer_2000_ddsd[1144];

    auto g_zz_0_0_0_xy_yz_0_zz = buffer_2000_ddsd[1145];

    auto g_zz_0_0_0_xy_zz_0_xx = buffer_2000_ddsd[1146];

    auto g_zz_0_0_0_xy_zz_0_xy = buffer_2000_ddsd[1147];

    auto g_zz_0_0_0_xy_zz_0_xz = buffer_2000_ddsd[1148];

    auto g_zz_0_0_0_xy_zz_0_yy = buffer_2000_ddsd[1149];

    auto g_zz_0_0_0_xy_zz_0_yz = buffer_2000_ddsd[1150];

    auto g_zz_0_0_0_xy_zz_0_zz = buffer_2000_ddsd[1151];

    auto g_zz_0_0_0_xz_xx_0_xx = buffer_2000_ddsd[1152];

    auto g_zz_0_0_0_xz_xx_0_xy = buffer_2000_ddsd[1153];

    auto g_zz_0_0_0_xz_xx_0_xz = buffer_2000_ddsd[1154];

    auto g_zz_0_0_0_xz_xx_0_yy = buffer_2000_ddsd[1155];

    auto g_zz_0_0_0_xz_xx_0_yz = buffer_2000_ddsd[1156];

    auto g_zz_0_0_0_xz_xx_0_zz = buffer_2000_ddsd[1157];

    auto g_zz_0_0_0_xz_xy_0_xx = buffer_2000_ddsd[1158];

    auto g_zz_0_0_0_xz_xy_0_xy = buffer_2000_ddsd[1159];

    auto g_zz_0_0_0_xz_xy_0_xz = buffer_2000_ddsd[1160];

    auto g_zz_0_0_0_xz_xy_0_yy = buffer_2000_ddsd[1161];

    auto g_zz_0_0_0_xz_xy_0_yz = buffer_2000_ddsd[1162];

    auto g_zz_0_0_0_xz_xy_0_zz = buffer_2000_ddsd[1163];

    auto g_zz_0_0_0_xz_xz_0_xx = buffer_2000_ddsd[1164];

    auto g_zz_0_0_0_xz_xz_0_xy = buffer_2000_ddsd[1165];

    auto g_zz_0_0_0_xz_xz_0_xz = buffer_2000_ddsd[1166];

    auto g_zz_0_0_0_xz_xz_0_yy = buffer_2000_ddsd[1167];

    auto g_zz_0_0_0_xz_xz_0_yz = buffer_2000_ddsd[1168];

    auto g_zz_0_0_0_xz_xz_0_zz = buffer_2000_ddsd[1169];

    auto g_zz_0_0_0_xz_yy_0_xx = buffer_2000_ddsd[1170];

    auto g_zz_0_0_0_xz_yy_0_xy = buffer_2000_ddsd[1171];

    auto g_zz_0_0_0_xz_yy_0_xz = buffer_2000_ddsd[1172];

    auto g_zz_0_0_0_xz_yy_0_yy = buffer_2000_ddsd[1173];

    auto g_zz_0_0_0_xz_yy_0_yz = buffer_2000_ddsd[1174];

    auto g_zz_0_0_0_xz_yy_0_zz = buffer_2000_ddsd[1175];

    auto g_zz_0_0_0_xz_yz_0_xx = buffer_2000_ddsd[1176];

    auto g_zz_0_0_0_xz_yz_0_xy = buffer_2000_ddsd[1177];

    auto g_zz_0_0_0_xz_yz_0_xz = buffer_2000_ddsd[1178];

    auto g_zz_0_0_0_xz_yz_0_yy = buffer_2000_ddsd[1179];

    auto g_zz_0_0_0_xz_yz_0_yz = buffer_2000_ddsd[1180];

    auto g_zz_0_0_0_xz_yz_0_zz = buffer_2000_ddsd[1181];

    auto g_zz_0_0_0_xz_zz_0_xx = buffer_2000_ddsd[1182];

    auto g_zz_0_0_0_xz_zz_0_xy = buffer_2000_ddsd[1183];

    auto g_zz_0_0_0_xz_zz_0_xz = buffer_2000_ddsd[1184];

    auto g_zz_0_0_0_xz_zz_0_yy = buffer_2000_ddsd[1185];

    auto g_zz_0_0_0_xz_zz_0_yz = buffer_2000_ddsd[1186];

    auto g_zz_0_0_0_xz_zz_0_zz = buffer_2000_ddsd[1187];

    auto g_zz_0_0_0_yy_xx_0_xx = buffer_2000_ddsd[1188];

    auto g_zz_0_0_0_yy_xx_0_xy = buffer_2000_ddsd[1189];

    auto g_zz_0_0_0_yy_xx_0_xz = buffer_2000_ddsd[1190];

    auto g_zz_0_0_0_yy_xx_0_yy = buffer_2000_ddsd[1191];

    auto g_zz_0_0_0_yy_xx_0_yz = buffer_2000_ddsd[1192];

    auto g_zz_0_0_0_yy_xx_0_zz = buffer_2000_ddsd[1193];

    auto g_zz_0_0_0_yy_xy_0_xx = buffer_2000_ddsd[1194];

    auto g_zz_0_0_0_yy_xy_0_xy = buffer_2000_ddsd[1195];

    auto g_zz_0_0_0_yy_xy_0_xz = buffer_2000_ddsd[1196];

    auto g_zz_0_0_0_yy_xy_0_yy = buffer_2000_ddsd[1197];

    auto g_zz_0_0_0_yy_xy_0_yz = buffer_2000_ddsd[1198];

    auto g_zz_0_0_0_yy_xy_0_zz = buffer_2000_ddsd[1199];

    auto g_zz_0_0_0_yy_xz_0_xx = buffer_2000_ddsd[1200];

    auto g_zz_0_0_0_yy_xz_0_xy = buffer_2000_ddsd[1201];

    auto g_zz_0_0_0_yy_xz_0_xz = buffer_2000_ddsd[1202];

    auto g_zz_0_0_0_yy_xz_0_yy = buffer_2000_ddsd[1203];

    auto g_zz_0_0_0_yy_xz_0_yz = buffer_2000_ddsd[1204];

    auto g_zz_0_0_0_yy_xz_0_zz = buffer_2000_ddsd[1205];

    auto g_zz_0_0_0_yy_yy_0_xx = buffer_2000_ddsd[1206];

    auto g_zz_0_0_0_yy_yy_0_xy = buffer_2000_ddsd[1207];

    auto g_zz_0_0_0_yy_yy_0_xz = buffer_2000_ddsd[1208];

    auto g_zz_0_0_0_yy_yy_0_yy = buffer_2000_ddsd[1209];

    auto g_zz_0_0_0_yy_yy_0_yz = buffer_2000_ddsd[1210];

    auto g_zz_0_0_0_yy_yy_0_zz = buffer_2000_ddsd[1211];

    auto g_zz_0_0_0_yy_yz_0_xx = buffer_2000_ddsd[1212];

    auto g_zz_0_0_0_yy_yz_0_xy = buffer_2000_ddsd[1213];

    auto g_zz_0_0_0_yy_yz_0_xz = buffer_2000_ddsd[1214];

    auto g_zz_0_0_0_yy_yz_0_yy = buffer_2000_ddsd[1215];

    auto g_zz_0_0_0_yy_yz_0_yz = buffer_2000_ddsd[1216];

    auto g_zz_0_0_0_yy_yz_0_zz = buffer_2000_ddsd[1217];

    auto g_zz_0_0_0_yy_zz_0_xx = buffer_2000_ddsd[1218];

    auto g_zz_0_0_0_yy_zz_0_xy = buffer_2000_ddsd[1219];

    auto g_zz_0_0_0_yy_zz_0_xz = buffer_2000_ddsd[1220];

    auto g_zz_0_0_0_yy_zz_0_yy = buffer_2000_ddsd[1221];

    auto g_zz_0_0_0_yy_zz_0_yz = buffer_2000_ddsd[1222];

    auto g_zz_0_0_0_yy_zz_0_zz = buffer_2000_ddsd[1223];

    auto g_zz_0_0_0_yz_xx_0_xx = buffer_2000_ddsd[1224];

    auto g_zz_0_0_0_yz_xx_0_xy = buffer_2000_ddsd[1225];

    auto g_zz_0_0_0_yz_xx_0_xz = buffer_2000_ddsd[1226];

    auto g_zz_0_0_0_yz_xx_0_yy = buffer_2000_ddsd[1227];

    auto g_zz_0_0_0_yz_xx_0_yz = buffer_2000_ddsd[1228];

    auto g_zz_0_0_0_yz_xx_0_zz = buffer_2000_ddsd[1229];

    auto g_zz_0_0_0_yz_xy_0_xx = buffer_2000_ddsd[1230];

    auto g_zz_0_0_0_yz_xy_0_xy = buffer_2000_ddsd[1231];

    auto g_zz_0_0_0_yz_xy_0_xz = buffer_2000_ddsd[1232];

    auto g_zz_0_0_0_yz_xy_0_yy = buffer_2000_ddsd[1233];

    auto g_zz_0_0_0_yz_xy_0_yz = buffer_2000_ddsd[1234];

    auto g_zz_0_0_0_yz_xy_0_zz = buffer_2000_ddsd[1235];

    auto g_zz_0_0_0_yz_xz_0_xx = buffer_2000_ddsd[1236];

    auto g_zz_0_0_0_yz_xz_0_xy = buffer_2000_ddsd[1237];

    auto g_zz_0_0_0_yz_xz_0_xz = buffer_2000_ddsd[1238];

    auto g_zz_0_0_0_yz_xz_0_yy = buffer_2000_ddsd[1239];

    auto g_zz_0_0_0_yz_xz_0_yz = buffer_2000_ddsd[1240];

    auto g_zz_0_0_0_yz_xz_0_zz = buffer_2000_ddsd[1241];

    auto g_zz_0_0_0_yz_yy_0_xx = buffer_2000_ddsd[1242];

    auto g_zz_0_0_0_yz_yy_0_xy = buffer_2000_ddsd[1243];

    auto g_zz_0_0_0_yz_yy_0_xz = buffer_2000_ddsd[1244];

    auto g_zz_0_0_0_yz_yy_0_yy = buffer_2000_ddsd[1245];

    auto g_zz_0_0_0_yz_yy_0_yz = buffer_2000_ddsd[1246];

    auto g_zz_0_0_0_yz_yy_0_zz = buffer_2000_ddsd[1247];

    auto g_zz_0_0_0_yz_yz_0_xx = buffer_2000_ddsd[1248];

    auto g_zz_0_0_0_yz_yz_0_xy = buffer_2000_ddsd[1249];

    auto g_zz_0_0_0_yz_yz_0_xz = buffer_2000_ddsd[1250];

    auto g_zz_0_0_0_yz_yz_0_yy = buffer_2000_ddsd[1251];

    auto g_zz_0_0_0_yz_yz_0_yz = buffer_2000_ddsd[1252];

    auto g_zz_0_0_0_yz_yz_0_zz = buffer_2000_ddsd[1253];

    auto g_zz_0_0_0_yz_zz_0_xx = buffer_2000_ddsd[1254];

    auto g_zz_0_0_0_yz_zz_0_xy = buffer_2000_ddsd[1255];

    auto g_zz_0_0_0_yz_zz_0_xz = buffer_2000_ddsd[1256];

    auto g_zz_0_0_0_yz_zz_0_yy = buffer_2000_ddsd[1257];

    auto g_zz_0_0_0_yz_zz_0_yz = buffer_2000_ddsd[1258];

    auto g_zz_0_0_0_yz_zz_0_zz = buffer_2000_ddsd[1259];

    auto g_zz_0_0_0_zz_xx_0_xx = buffer_2000_ddsd[1260];

    auto g_zz_0_0_0_zz_xx_0_xy = buffer_2000_ddsd[1261];

    auto g_zz_0_0_0_zz_xx_0_xz = buffer_2000_ddsd[1262];

    auto g_zz_0_0_0_zz_xx_0_yy = buffer_2000_ddsd[1263];

    auto g_zz_0_0_0_zz_xx_0_yz = buffer_2000_ddsd[1264];

    auto g_zz_0_0_0_zz_xx_0_zz = buffer_2000_ddsd[1265];

    auto g_zz_0_0_0_zz_xy_0_xx = buffer_2000_ddsd[1266];

    auto g_zz_0_0_0_zz_xy_0_xy = buffer_2000_ddsd[1267];

    auto g_zz_0_0_0_zz_xy_0_xz = buffer_2000_ddsd[1268];

    auto g_zz_0_0_0_zz_xy_0_yy = buffer_2000_ddsd[1269];

    auto g_zz_0_0_0_zz_xy_0_yz = buffer_2000_ddsd[1270];

    auto g_zz_0_0_0_zz_xy_0_zz = buffer_2000_ddsd[1271];

    auto g_zz_0_0_0_zz_xz_0_xx = buffer_2000_ddsd[1272];

    auto g_zz_0_0_0_zz_xz_0_xy = buffer_2000_ddsd[1273];

    auto g_zz_0_0_0_zz_xz_0_xz = buffer_2000_ddsd[1274];

    auto g_zz_0_0_0_zz_xz_0_yy = buffer_2000_ddsd[1275];

    auto g_zz_0_0_0_zz_xz_0_yz = buffer_2000_ddsd[1276];

    auto g_zz_0_0_0_zz_xz_0_zz = buffer_2000_ddsd[1277];

    auto g_zz_0_0_0_zz_yy_0_xx = buffer_2000_ddsd[1278];

    auto g_zz_0_0_0_zz_yy_0_xy = buffer_2000_ddsd[1279];

    auto g_zz_0_0_0_zz_yy_0_xz = buffer_2000_ddsd[1280];

    auto g_zz_0_0_0_zz_yy_0_yy = buffer_2000_ddsd[1281];

    auto g_zz_0_0_0_zz_yy_0_yz = buffer_2000_ddsd[1282];

    auto g_zz_0_0_0_zz_yy_0_zz = buffer_2000_ddsd[1283];

    auto g_zz_0_0_0_zz_yz_0_xx = buffer_2000_ddsd[1284];

    auto g_zz_0_0_0_zz_yz_0_xy = buffer_2000_ddsd[1285];

    auto g_zz_0_0_0_zz_yz_0_xz = buffer_2000_ddsd[1286];

    auto g_zz_0_0_0_zz_yz_0_yy = buffer_2000_ddsd[1287];

    auto g_zz_0_0_0_zz_yz_0_yz = buffer_2000_ddsd[1288];

    auto g_zz_0_0_0_zz_yz_0_zz = buffer_2000_ddsd[1289];

    auto g_zz_0_0_0_zz_zz_0_xx = buffer_2000_ddsd[1290];

    auto g_zz_0_0_0_zz_zz_0_xy = buffer_2000_ddsd[1291];

    auto g_zz_0_0_0_zz_zz_0_xz = buffer_2000_ddsd[1292];

    auto g_zz_0_0_0_zz_zz_0_yy = buffer_2000_ddsd[1293];

    auto g_zz_0_0_0_zz_zz_0_yz = buffer_2000_ddsd[1294];

    auto g_zz_0_0_0_zz_zz_0_zz = buffer_2000_ddsd[1295];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_xx_0_xx, g_0_xx_0_xy, g_0_xx_0_xz, g_0_xx_0_yy, g_0_xx_0_yz, g_0_xx_0_zz, g_xx_0_0_0_xx_xx_0_xx, g_xx_0_0_0_xx_xx_0_xy, g_xx_0_0_0_xx_xx_0_xz, g_xx_0_0_0_xx_xx_0_yy, g_xx_0_0_0_xx_xx_0_yz, g_xx_0_0_0_xx_xx_0_zz, g_xx_xx_0_xx, g_xx_xx_0_xy, g_xx_xx_0_xz, g_xx_xx_0_yy, g_xx_xx_0_yz, g_xx_xx_0_zz, g_xxxx_xx_0_xx, g_xxxx_xx_0_xy, g_xxxx_xx_0_xz, g_xxxx_xx_0_yy, g_xxxx_xx_0_yz, g_xxxx_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_xx_0_xx[i] = 2.0 * g_0_xx_0_xx[i] - 10.0 * g_xx_xx_0_xx[i] * a_exp + 4.0 * g_xxxx_xx_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xx_0_xy[i] = 2.0 * g_0_xx_0_xy[i] - 10.0 * g_xx_xx_0_xy[i] * a_exp + 4.0 * g_xxxx_xx_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xx_0_xz[i] = 2.0 * g_0_xx_0_xz[i] - 10.0 * g_xx_xx_0_xz[i] * a_exp + 4.0 * g_xxxx_xx_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xx_0_yy[i] = 2.0 * g_0_xx_0_yy[i] - 10.0 * g_xx_xx_0_yy[i] * a_exp + 4.0 * g_xxxx_xx_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xx_0_yz[i] = 2.0 * g_0_xx_0_yz[i] - 10.0 * g_xx_xx_0_yz[i] * a_exp + 4.0 * g_xxxx_xx_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xx_0_zz[i] = 2.0 * g_0_xx_0_zz[i] - 10.0 * g_xx_xx_0_zz[i] * a_exp + 4.0 * g_xxxx_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_xx_0_0_0_xx_xy_0_xx, g_xx_0_0_0_xx_xy_0_xy, g_xx_0_0_0_xx_xy_0_xz, g_xx_0_0_0_xx_xy_0_yy, g_xx_0_0_0_xx_xy_0_yz, g_xx_0_0_0_xx_xy_0_zz, g_xx_xy_0_xx, g_xx_xy_0_xy, g_xx_xy_0_xz, g_xx_xy_0_yy, g_xx_xy_0_yz, g_xx_xy_0_zz, g_xxxx_xy_0_xx, g_xxxx_xy_0_xy, g_xxxx_xy_0_xz, g_xxxx_xy_0_yy, g_xxxx_xy_0_yz, g_xxxx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_xy_0_xx[i] = 2.0 * g_0_xy_0_xx[i] - 10.0 * g_xx_xy_0_xx[i] * a_exp + 4.0 * g_xxxx_xy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xy_0_xy[i] = 2.0 * g_0_xy_0_xy[i] - 10.0 * g_xx_xy_0_xy[i] * a_exp + 4.0 * g_xxxx_xy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xy_0_xz[i] = 2.0 * g_0_xy_0_xz[i] - 10.0 * g_xx_xy_0_xz[i] * a_exp + 4.0 * g_xxxx_xy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xy_0_yy[i] = 2.0 * g_0_xy_0_yy[i] - 10.0 * g_xx_xy_0_yy[i] * a_exp + 4.0 * g_xxxx_xy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xy_0_yz[i] = 2.0 * g_0_xy_0_yz[i] - 10.0 * g_xx_xy_0_yz[i] * a_exp + 4.0 * g_xxxx_xy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xy_0_zz[i] = 2.0 * g_0_xy_0_zz[i] - 10.0 * g_xx_xy_0_zz[i] * a_exp + 4.0 * g_xxxx_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_xx_0_0_0_xx_xz_0_xx, g_xx_0_0_0_xx_xz_0_xy, g_xx_0_0_0_xx_xz_0_xz, g_xx_0_0_0_xx_xz_0_yy, g_xx_0_0_0_xx_xz_0_yz, g_xx_0_0_0_xx_xz_0_zz, g_xx_xz_0_xx, g_xx_xz_0_xy, g_xx_xz_0_xz, g_xx_xz_0_yy, g_xx_xz_0_yz, g_xx_xz_0_zz, g_xxxx_xz_0_xx, g_xxxx_xz_0_xy, g_xxxx_xz_0_xz, g_xxxx_xz_0_yy, g_xxxx_xz_0_yz, g_xxxx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_xz_0_xx[i] = 2.0 * g_0_xz_0_xx[i] - 10.0 * g_xx_xz_0_xx[i] * a_exp + 4.0 * g_xxxx_xz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xz_0_xy[i] = 2.0 * g_0_xz_0_xy[i] - 10.0 * g_xx_xz_0_xy[i] * a_exp + 4.0 * g_xxxx_xz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xz_0_xz[i] = 2.0 * g_0_xz_0_xz[i] - 10.0 * g_xx_xz_0_xz[i] * a_exp + 4.0 * g_xxxx_xz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xz_0_yy[i] = 2.0 * g_0_xz_0_yy[i] - 10.0 * g_xx_xz_0_yy[i] * a_exp + 4.0 * g_xxxx_xz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xz_0_yz[i] = 2.0 * g_0_xz_0_yz[i] - 10.0 * g_xx_xz_0_yz[i] * a_exp + 4.0 * g_xxxx_xz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_xz_0_zz[i] = 2.0 * g_0_xz_0_zz[i] - 10.0 * g_xx_xz_0_zz[i] * a_exp + 4.0 * g_xxxx_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_yy_0_xx, g_0_yy_0_xy, g_0_yy_0_xz, g_0_yy_0_yy, g_0_yy_0_yz, g_0_yy_0_zz, g_xx_0_0_0_xx_yy_0_xx, g_xx_0_0_0_xx_yy_0_xy, g_xx_0_0_0_xx_yy_0_xz, g_xx_0_0_0_xx_yy_0_yy, g_xx_0_0_0_xx_yy_0_yz, g_xx_0_0_0_xx_yy_0_zz, g_xx_yy_0_xx, g_xx_yy_0_xy, g_xx_yy_0_xz, g_xx_yy_0_yy, g_xx_yy_0_yz, g_xx_yy_0_zz, g_xxxx_yy_0_xx, g_xxxx_yy_0_xy, g_xxxx_yy_0_xz, g_xxxx_yy_0_yy, g_xxxx_yy_0_yz, g_xxxx_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_yy_0_xx[i] = 2.0 * g_0_yy_0_xx[i] - 10.0 * g_xx_yy_0_xx[i] * a_exp + 4.0 * g_xxxx_yy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yy_0_xy[i] = 2.0 * g_0_yy_0_xy[i] - 10.0 * g_xx_yy_0_xy[i] * a_exp + 4.0 * g_xxxx_yy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yy_0_xz[i] = 2.0 * g_0_yy_0_xz[i] - 10.0 * g_xx_yy_0_xz[i] * a_exp + 4.0 * g_xxxx_yy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yy_0_yy[i] = 2.0 * g_0_yy_0_yy[i] - 10.0 * g_xx_yy_0_yy[i] * a_exp + 4.0 * g_xxxx_yy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yy_0_yz[i] = 2.0 * g_0_yy_0_yz[i] - 10.0 * g_xx_yy_0_yz[i] * a_exp + 4.0 * g_xxxx_yy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yy_0_zz[i] = 2.0 * g_0_yy_0_zz[i] - 10.0 * g_xx_yy_0_zz[i] * a_exp + 4.0 * g_xxxx_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_xx_0_0_0_xx_yz_0_xx, g_xx_0_0_0_xx_yz_0_xy, g_xx_0_0_0_xx_yz_0_xz, g_xx_0_0_0_xx_yz_0_yy, g_xx_0_0_0_xx_yz_0_yz, g_xx_0_0_0_xx_yz_0_zz, g_xx_yz_0_xx, g_xx_yz_0_xy, g_xx_yz_0_xz, g_xx_yz_0_yy, g_xx_yz_0_yz, g_xx_yz_0_zz, g_xxxx_yz_0_xx, g_xxxx_yz_0_xy, g_xxxx_yz_0_xz, g_xxxx_yz_0_yy, g_xxxx_yz_0_yz, g_xxxx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_yz_0_xx[i] = 2.0 * g_0_yz_0_xx[i] - 10.0 * g_xx_yz_0_xx[i] * a_exp + 4.0 * g_xxxx_yz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yz_0_xy[i] = 2.0 * g_0_yz_0_xy[i] - 10.0 * g_xx_yz_0_xy[i] * a_exp + 4.0 * g_xxxx_yz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yz_0_xz[i] = 2.0 * g_0_yz_0_xz[i] - 10.0 * g_xx_yz_0_xz[i] * a_exp + 4.0 * g_xxxx_yz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yz_0_yy[i] = 2.0 * g_0_yz_0_yy[i] - 10.0 * g_xx_yz_0_yy[i] * a_exp + 4.0 * g_xxxx_yz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yz_0_yz[i] = 2.0 * g_0_yz_0_yz[i] - 10.0 * g_xx_yz_0_yz[i] * a_exp + 4.0 * g_xxxx_yz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_yz_0_zz[i] = 2.0 * g_0_yz_0_zz[i] - 10.0 * g_xx_yz_0_zz[i] * a_exp + 4.0 * g_xxxx_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_zz_0_xx, g_0_zz_0_xy, g_0_zz_0_xz, g_0_zz_0_yy, g_0_zz_0_yz, g_0_zz_0_zz, g_xx_0_0_0_xx_zz_0_xx, g_xx_0_0_0_xx_zz_0_xy, g_xx_0_0_0_xx_zz_0_xz, g_xx_0_0_0_xx_zz_0_yy, g_xx_0_0_0_xx_zz_0_yz, g_xx_0_0_0_xx_zz_0_zz, g_xx_zz_0_xx, g_xx_zz_0_xy, g_xx_zz_0_xz, g_xx_zz_0_yy, g_xx_zz_0_yz, g_xx_zz_0_zz, g_xxxx_zz_0_xx, g_xxxx_zz_0_xy, g_xxxx_zz_0_xz, g_xxxx_zz_0_yy, g_xxxx_zz_0_yz, g_xxxx_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xx_zz_0_xx[i] = 2.0 * g_0_zz_0_xx[i] - 10.0 * g_xx_zz_0_xx[i] * a_exp + 4.0 * g_xxxx_zz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_zz_0_xy[i] = 2.0 * g_0_zz_0_xy[i] - 10.0 * g_xx_zz_0_xy[i] * a_exp + 4.0 * g_xxxx_zz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_zz_0_xz[i] = 2.0 * g_0_zz_0_xz[i] - 10.0 * g_xx_zz_0_xz[i] * a_exp + 4.0 * g_xxxx_zz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_zz_0_yy[i] = 2.0 * g_0_zz_0_yy[i] - 10.0 * g_xx_zz_0_yy[i] * a_exp + 4.0 * g_xxxx_zz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_zz_0_yz[i] = 2.0 * g_0_zz_0_yz[i] - 10.0 * g_xx_zz_0_yz[i] * a_exp + 4.0 * g_xxxx_zz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xx_zz_0_zz[i] = 2.0 * g_0_zz_0_zz[i] - 10.0 * g_xx_zz_0_zz[i] * a_exp + 4.0 * g_xxxx_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_xx_0_0_0_xy_xx_0_xx, g_xx_0_0_0_xy_xx_0_xy, g_xx_0_0_0_xy_xx_0_xz, g_xx_0_0_0_xy_xx_0_yy, g_xx_0_0_0_xy_xx_0_yz, g_xx_0_0_0_xy_xx_0_zz, g_xxxy_xx_0_xx, g_xxxy_xx_0_xy, g_xxxy_xx_0_xz, g_xxxy_xx_0_yy, g_xxxy_xx_0_yz, g_xxxy_xx_0_zz, g_xy_xx_0_xx, g_xy_xx_0_xy, g_xy_xx_0_xz, g_xy_xx_0_yy, g_xy_xx_0_yz, g_xy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_xx_0_xx[i] = -6.0 * g_xy_xx_0_xx[i] * a_exp + 4.0 * g_xxxy_xx_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xx_0_xy[i] = -6.0 * g_xy_xx_0_xy[i] * a_exp + 4.0 * g_xxxy_xx_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xx_0_xz[i] = -6.0 * g_xy_xx_0_xz[i] * a_exp + 4.0 * g_xxxy_xx_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xx_0_yy[i] = -6.0 * g_xy_xx_0_yy[i] * a_exp + 4.0 * g_xxxy_xx_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xx_0_yz[i] = -6.0 * g_xy_xx_0_yz[i] * a_exp + 4.0 * g_xxxy_xx_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xx_0_zz[i] = -6.0 * g_xy_xx_0_zz[i] * a_exp + 4.0 * g_xxxy_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_xx_0_0_0_xy_xy_0_xx, g_xx_0_0_0_xy_xy_0_xy, g_xx_0_0_0_xy_xy_0_xz, g_xx_0_0_0_xy_xy_0_yy, g_xx_0_0_0_xy_xy_0_yz, g_xx_0_0_0_xy_xy_0_zz, g_xxxy_xy_0_xx, g_xxxy_xy_0_xy, g_xxxy_xy_0_xz, g_xxxy_xy_0_yy, g_xxxy_xy_0_yz, g_xxxy_xy_0_zz, g_xy_xy_0_xx, g_xy_xy_0_xy, g_xy_xy_0_xz, g_xy_xy_0_yy, g_xy_xy_0_yz, g_xy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_xy_0_xx[i] = -6.0 * g_xy_xy_0_xx[i] * a_exp + 4.0 * g_xxxy_xy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xy_0_xy[i] = -6.0 * g_xy_xy_0_xy[i] * a_exp + 4.0 * g_xxxy_xy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xy_0_xz[i] = -6.0 * g_xy_xy_0_xz[i] * a_exp + 4.0 * g_xxxy_xy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xy_0_yy[i] = -6.0 * g_xy_xy_0_yy[i] * a_exp + 4.0 * g_xxxy_xy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xy_0_yz[i] = -6.0 * g_xy_xy_0_yz[i] * a_exp + 4.0 * g_xxxy_xy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xy_0_zz[i] = -6.0 * g_xy_xy_0_zz[i] * a_exp + 4.0 * g_xxxy_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_xx_0_0_0_xy_xz_0_xx, g_xx_0_0_0_xy_xz_0_xy, g_xx_0_0_0_xy_xz_0_xz, g_xx_0_0_0_xy_xz_0_yy, g_xx_0_0_0_xy_xz_0_yz, g_xx_0_0_0_xy_xz_0_zz, g_xxxy_xz_0_xx, g_xxxy_xz_0_xy, g_xxxy_xz_0_xz, g_xxxy_xz_0_yy, g_xxxy_xz_0_yz, g_xxxy_xz_0_zz, g_xy_xz_0_xx, g_xy_xz_0_xy, g_xy_xz_0_xz, g_xy_xz_0_yy, g_xy_xz_0_yz, g_xy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_xz_0_xx[i] = -6.0 * g_xy_xz_0_xx[i] * a_exp + 4.0 * g_xxxy_xz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xz_0_xy[i] = -6.0 * g_xy_xz_0_xy[i] * a_exp + 4.0 * g_xxxy_xz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xz_0_xz[i] = -6.0 * g_xy_xz_0_xz[i] * a_exp + 4.0 * g_xxxy_xz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xz_0_yy[i] = -6.0 * g_xy_xz_0_yy[i] * a_exp + 4.0 * g_xxxy_xz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xz_0_yz[i] = -6.0 * g_xy_xz_0_yz[i] * a_exp + 4.0 * g_xxxy_xz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_xz_0_zz[i] = -6.0 * g_xy_xz_0_zz[i] * a_exp + 4.0 * g_xxxy_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_xx_0_0_0_xy_yy_0_xx, g_xx_0_0_0_xy_yy_0_xy, g_xx_0_0_0_xy_yy_0_xz, g_xx_0_0_0_xy_yy_0_yy, g_xx_0_0_0_xy_yy_0_yz, g_xx_0_0_0_xy_yy_0_zz, g_xxxy_yy_0_xx, g_xxxy_yy_0_xy, g_xxxy_yy_0_xz, g_xxxy_yy_0_yy, g_xxxy_yy_0_yz, g_xxxy_yy_0_zz, g_xy_yy_0_xx, g_xy_yy_0_xy, g_xy_yy_0_xz, g_xy_yy_0_yy, g_xy_yy_0_yz, g_xy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_yy_0_xx[i] = -6.0 * g_xy_yy_0_xx[i] * a_exp + 4.0 * g_xxxy_yy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yy_0_xy[i] = -6.0 * g_xy_yy_0_xy[i] * a_exp + 4.0 * g_xxxy_yy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yy_0_xz[i] = -6.0 * g_xy_yy_0_xz[i] * a_exp + 4.0 * g_xxxy_yy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yy_0_yy[i] = -6.0 * g_xy_yy_0_yy[i] * a_exp + 4.0 * g_xxxy_yy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yy_0_yz[i] = -6.0 * g_xy_yy_0_yz[i] * a_exp + 4.0 * g_xxxy_yy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yy_0_zz[i] = -6.0 * g_xy_yy_0_zz[i] * a_exp + 4.0 * g_xxxy_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_xx_0_0_0_xy_yz_0_xx, g_xx_0_0_0_xy_yz_0_xy, g_xx_0_0_0_xy_yz_0_xz, g_xx_0_0_0_xy_yz_0_yy, g_xx_0_0_0_xy_yz_0_yz, g_xx_0_0_0_xy_yz_0_zz, g_xxxy_yz_0_xx, g_xxxy_yz_0_xy, g_xxxy_yz_0_xz, g_xxxy_yz_0_yy, g_xxxy_yz_0_yz, g_xxxy_yz_0_zz, g_xy_yz_0_xx, g_xy_yz_0_xy, g_xy_yz_0_xz, g_xy_yz_0_yy, g_xy_yz_0_yz, g_xy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_yz_0_xx[i] = -6.0 * g_xy_yz_0_xx[i] * a_exp + 4.0 * g_xxxy_yz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yz_0_xy[i] = -6.0 * g_xy_yz_0_xy[i] * a_exp + 4.0 * g_xxxy_yz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yz_0_xz[i] = -6.0 * g_xy_yz_0_xz[i] * a_exp + 4.0 * g_xxxy_yz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yz_0_yy[i] = -6.0 * g_xy_yz_0_yy[i] * a_exp + 4.0 * g_xxxy_yz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yz_0_yz[i] = -6.0 * g_xy_yz_0_yz[i] * a_exp + 4.0 * g_xxxy_yz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_yz_0_zz[i] = -6.0 * g_xy_yz_0_zz[i] * a_exp + 4.0 * g_xxxy_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_xx_0_0_0_xy_zz_0_xx, g_xx_0_0_0_xy_zz_0_xy, g_xx_0_0_0_xy_zz_0_xz, g_xx_0_0_0_xy_zz_0_yy, g_xx_0_0_0_xy_zz_0_yz, g_xx_0_0_0_xy_zz_0_zz, g_xxxy_zz_0_xx, g_xxxy_zz_0_xy, g_xxxy_zz_0_xz, g_xxxy_zz_0_yy, g_xxxy_zz_0_yz, g_xxxy_zz_0_zz, g_xy_zz_0_xx, g_xy_zz_0_xy, g_xy_zz_0_xz, g_xy_zz_0_yy, g_xy_zz_0_yz, g_xy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xy_zz_0_xx[i] = -6.0 * g_xy_zz_0_xx[i] * a_exp + 4.0 * g_xxxy_zz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_zz_0_xy[i] = -6.0 * g_xy_zz_0_xy[i] * a_exp + 4.0 * g_xxxy_zz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_zz_0_xz[i] = -6.0 * g_xy_zz_0_xz[i] * a_exp + 4.0 * g_xxxy_zz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_zz_0_yy[i] = -6.0 * g_xy_zz_0_yy[i] * a_exp + 4.0 * g_xxxy_zz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_zz_0_yz[i] = -6.0 * g_xy_zz_0_yz[i] * a_exp + 4.0 * g_xxxy_zz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xy_zz_0_zz[i] = -6.0 * g_xy_zz_0_zz[i] * a_exp + 4.0 * g_xxxy_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_xx_0_0_0_xz_xx_0_xx, g_xx_0_0_0_xz_xx_0_xy, g_xx_0_0_0_xz_xx_0_xz, g_xx_0_0_0_xz_xx_0_yy, g_xx_0_0_0_xz_xx_0_yz, g_xx_0_0_0_xz_xx_0_zz, g_xxxz_xx_0_xx, g_xxxz_xx_0_xy, g_xxxz_xx_0_xz, g_xxxz_xx_0_yy, g_xxxz_xx_0_yz, g_xxxz_xx_0_zz, g_xz_xx_0_xx, g_xz_xx_0_xy, g_xz_xx_0_xz, g_xz_xx_0_yy, g_xz_xx_0_yz, g_xz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_xx_0_xx[i] = -6.0 * g_xz_xx_0_xx[i] * a_exp + 4.0 * g_xxxz_xx_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xx_0_xy[i] = -6.0 * g_xz_xx_0_xy[i] * a_exp + 4.0 * g_xxxz_xx_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xx_0_xz[i] = -6.0 * g_xz_xx_0_xz[i] * a_exp + 4.0 * g_xxxz_xx_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xx_0_yy[i] = -6.0 * g_xz_xx_0_yy[i] * a_exp + 4.0 * g_xxxz_xx_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xx_0_yz[i] = -6.0 * g_xz_xx_0_yz[i] * a_exp + 4.0 * g_xxxz_xx_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xx_0_zz[i] = -6.0 * g_xz_xx_0_zz[i] * a_exp + 4.0 * g_xxxz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_xx_0_0_0_xz_xy_0_xx, g_xx_0_0_0_xz_xy_0_xy, g_xx_0_0_0_xz_xy_0_xz, g_xx_0_0_0_xz_xy_0_yy, g_xx_0_0_0_xz_xy_0_yz, g_xx_0_0_0_xz_xy_0_zz, g_xxxz_xy_0_xx, g_xxxz_xy_0_xy, g_xxxz_xy_0_xz, g_xxxz_xy_0_yy, g_xxxz_xy_0_yz, g_xxxz_xy_0_zz, g_xz_xy_0_xx, g_xz_xy_0_xy, g_xz_xy_0_xz, g_xz_xy_0_yy, g_xz_xy_0_yz, g_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_xy_0_xx[i] = -6.0 * g_xz_xy_0_xx[i] * a_exp + 4.0 * g_xxxz_xy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xy_0_xy[i] = -6.0 * g_xz_xy_0_xy[i] * a_exp + 4.0 * g_xxxz_xy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xy_0_xz[i] = -6.0 * g_xz_xy_0_xz[i] * a_exp + 4.0 * g_xxxz_xy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xy_0_yy[i] = -6.0 * g_xz_xy_0_yy[i] * a_exp + 4.0 * g_xxxz_xy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xy_0_yz[i] = -6.0 * g_xz_xy_0_yz[i] * a_exp + 4.0 * g_xxxz_xy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xy_0_zz[i] = -6.0 * g_xz_xy_0_zz[i] * a_exp + 4.0 * g_xxxz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_xx_0_0_0_xz_xz_0_xx, g_xx_0_0_0_xz_xz_0_xy, g_xx_0_0_0_xz_xz_0_xz, g_xx_0_0_0_xz_xz_0_yy, g_xx_0_0_0_xz_xz_0_yz, g_xx_0_0_0_xz_xz_0_zz, g_xxxz_xz_0_xx, g_xxxz_xz_0_xy, g_xxxz_xz_0_xz, g_xxxz_xz_0_yy, g_xxxz_xz_0_yz, g_xxxz_xz_0_zz, g_xz_xz_0_xx, g_xz_xz_0_xy, g_xz_xz_0_xz, g_xz_xz_0_yy, g_xz_xz_0_yz, g_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_xz_0_xx[i] = -6.0 * g_xz_xz_0_xx[i] * a_exp + 4.0 * g_xxxz_xz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xz_0_xy[i] = -6.0 * g_xz_xz_0_xy[i] * a_exp + 4.0 * g_xxxz_xz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xz_0_xz[i] = -6.0 * g_xz_xz_0_xz[i] * a_exp + 4.0 * g_xxxz_xz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xz_0_yy[i] = -6.0 * g_xz_xz_0_yy[i] * a_exp + 4.0 * g_xxxz_xz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xz_0_yz[i] = -6.0 * g_xz_xz_0_yz[i] * a_exp + 4.0 * g_xxxz_xz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_xz_0_zz[i] = -6.0 * g_xz_xz_0_zz[i] * a_exp + 4.0 * g_xxxz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_xx_0_0_0_xz_yy_0_xx, g_xx_0_0_0_xz_yy_0_xy, g_xx_0_0_0_xz_yy_0_xz, g_xx_0_0_0_xz_yy_0_yy, g_xx_0_0_0_xz_yy_0_yz, g_xx_0_0_0_xz_yy_0_zz, g_xxxz_yy_0_xx, g_xxxz_yy_0_xy, g_xxxz_yy_0_xz, g_xxxz_yy_0_yy, g_xxxz_yy_0_yz, g_xxxz_yy_0_zz, g_xz_yy_0_xx, g_xz_yy_0_xy, g_xz_yy_0_xz, g_xz_yy_0_yy, g_xz_yy_0_yz, g_xz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_yy_0_xx[i] = -6.0 * g_xz_yy_0_xx[i] * a_exp + 4.0 * g_xxxz_yy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yy_0_xy[i] = -6.0 * g_xz_yy_0_xy[i] * a_exp + 4.0 * g_xxxz_yy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yy_0_xz[i] = -6.0 * g_xz_yy_0_xz[i] * a_exp + 4.0 * g_xxxz_yy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yy_0_yy[i] = -6.0 * g_xz_yy_0_yy[i] * a_exp + 4.0 * g_xxxz_yy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yy_0_yz[i] = -6.0 * g_xz_yy_0_yz[i] * a_exp + 4.0 * g_xxxz_yy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yy_0_zz[i] = -6.0 * g_xz_yy_0_zz[i] * a_exp + 4.0 * g_xxxz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_xx_0_0_0_xz_yz_0_xx, g_xx_0_0_0_xz_yz_0_xy, g_xx_0_0_0_xz_yz_0_xz, g_xx_0_0_0_xz_yz_0_yy, g_xx_0_0_0_xz_yz_0_yz, g_xx_0_0_0_xz_yz_0_zz, g_xxxz_yz_0_xx, g_xxxz_yz_0_xy, g_xxxz_yz_0_xz, g_xxxz_yz_0_yy, g_xxxz_yz_0_yz, g_xxxz_yz_0_zz, g_xz_yz_0_xx, g_xz_yz_0_xy, g_xz_yz_0_xz, g_xz_yz_0_yy, g_xz_yz_0_yz, g_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_yz_0_xx[i] = -6.0 * g_xz_yz_0_xx[i] * a_exp + 4.0 * g_xxxz_yz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yz_0_xy[i] = -6.0 * g_xz_yz_0_xy[i] * a_exp + 4.0 * g_xxxz_yz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yz_0_xz[i] = -6.0 * g_xz_yz_0_xz[i] * a_exp + 4.0 * g_xxxz_yz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yz_0_yy[i] = -6.0 * g_xz_yz_0_yy[i] * a_exp + 4.0 * g_xxxz_yz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yz_0_yz[i] = -6.0 * g_xz_yz_0_yz[i] * a_exp + 4.0 * g_xxxz_yz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_yz_0_zz[i] = -6.0 * g_xz_yz_0_zz[i] * a_exp + 4.0 * g_xxxz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_xx_0_0_0_xz_zz_0_xx, g_xx_0_0_0_xz_zz_0_xy, g_xx_0_0_0_xz_zz_0_xz, g_xx_0_0_0_xz_zz_0_yy, g_xx_0_0_0_xz_zz_0_yz, g_xx_0_0_0_xz_zz_0_zz, g_xxxz_zz_0_xx, g_xxxz_zz_0_xy, g_xxxz_zz_0_xz, g_xxxz_zz_0_yy, g_xxxz_zz_0_yz, g_xxxz_zz_0_zz, g_xz_zz_0_xx, g_xz_zz_0_xy, g_xz_zz_0_xz, g_xz_zz_0_yy, g_xz_zz_0_yz, g_xz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_xz_zz_0_xx[i] = -6.0 * g_xz_zz_0_xx[i] * a_exp + 4.0 * g_xxxz_zz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_zz_0_xy[i] = -6.0 * g_xz_zz_0_xy[i] * a_exp + 4.0 * g_xxxz_zz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_zz_0_xz[i] = -6.0 * g_xz_zz_0_xz[i] * a_exp + 4.0 * g_xxxz_zz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_zz_0_yy[i] = -6.0 * g_xz_zz_0_yy[i] * a_exp + 4.0 * g_xxxz_zz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_zz_0_yz[i] = -6.0 * g_xz_zz_0_yz[i] * a_exp + 4.0 * g_xxxz_zz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_xz_zz_0_zz[i] = -6.0 * g_xz_zz_0_zz[i] * a_exp + 4.0 * g_xxxz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_xx_0_0_0_yy_xx_0_xx, g_xx_0_0_0_yy_xx_0_xy, g_xx_0_0_0_yy_xx_0_xz, g_xx_0_0_0_yy_xx_0_yy, g_xx_0_0_0_yy_xx_0_yz, g_xx_0_0_0_yy_xx_0_zz, g_xxyy_xx_0_xx, g_xxyy_xx_0_xy, g_xxyy_xx_0_xz, g_xxyy_xx_0_yy, g_xxyy_xx_0_yz, g_xxyy_xx_0_zz, g_yy_xx_0_xx, g_yy_xx_0_xy, g_yy_xx_0_xz, g_yy_xx_0_yy, g_yy_xx_0_yz, g_yy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_xx_0_xx[i] = -2.0 * g_yy_xx_0_xx[i] * a_exp + 4.0 * g_xxyy_xx_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xx_0_xy[i] = -2.0 * g_yy_xx_0_xy[i] * a_exp + 4.0 * g_xxyy_xx_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xx_0_xz[i] = -2.0 * g_yy_xx_0_xz[i] * a_exp + 4.0 * g_xxyy_xx_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xx_0_yy[i] = -2.0 * g_yy_xx_0_yy[i] * a_exp + 4.0 * g_xxyy_xx_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xx_0_yz[i] = -2.0 * g_yy_xx_0_yz[i] * a_exp + 4.0 * g_xxyy_xx_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xx_0_zz[i] = -2.0 * g_yy_xx_0_zz[i] * a_exp + 4.0 * g_xxyy_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_xx_0_0_0_yy_xy_0_xx, g_xx_0_0_0_yy_xy_0_xy, g_xx_0_0_0_yy_xy_0_xz, g_xx_0_0_0_yy_xy_0_yy, g_xx_0_0_0_yy_xy_0_yz, g_xx_0_0_0_yy_xy_0_zz, g_xxyy_xy_0_xx, g_xxyy_xy_0_xy, g_xxyy_xy_0_xz, g_xxyy_xy_0_yy, g_xxyy_xy_0_yz, g_xxyy_xy_0_zz, g_yy_xy_0_xx, g_yy_xy_0_xy, g_yy_xy_0_xz, g_yy_xy_0_yy, g_yy_xy_0_yz, g_yy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_xy_0_xx[i] = -2.0 * g_yy_xy_0_xx[i] * a_exp + 4.0 * g_xxyy_xy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xy_0_xy[i] = -2.0 * g_yy_xy_0_xy[i] * a_exp + 4.0 * g_xxyy_xy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xy_0_xz[i] = -2.0 * g_yy_xy_0_xz[i] * a_exp + 4.0 * g_xxyy_xy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xy_0_yy[i] = -2.0 * g_yy_xy_0_yy[i] * a_exp + 4.0 * g_xxyy_xy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xy_0_yz[i] = -2.0 * g_yy_xy_0_yz[i] * a_exp + 4.0 * g_xxyy_xy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xy_0_zz[i] = -2.0 * g_yy_xy_0_zz[i] * a_exp + 4.0 * g_xxyy_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_xx_0_0_0_yy_xz_0_xx, g_xx_0_0_0_yy_xz_0_xy, g_xx_0_0_0_yy_xz_0_xz, g_xx_0_0_0_yy_xz_0_yy, g_xx_0_0_0_yy_xz_0_yz, g_xx_0_0_0_yy_xz_0_zz, g_xxyy_xz_0_xx, g_xxyy_xz_0_xy, g_xxyy_xz_0_xz, g_xxyy_xz_0_yy, g_xxyy_xz_0_yz, g_xxyy_xz_0_zz, g_yy_xz_0_xx, g_yy_xz_0_xy, g_yy_xz_0_xz, g_yy_xz_0_yy, g_yy_xz_0_yz, g_yy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_xz_0_xx[i] = -2.0 * g_yy_xz_0_xx[i] * a_exp + 4.0 * g_xxyy_xz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xz_0_xy[i] = -2.0 * g_yy_xz_0_xy[i] * a_exp + 4.0 * g_xxyy_xz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xz_0_xz[i] = -2.0 * g_yy_xz_0_xz[i] * a_exp + 4.0 * g_xxyy_xz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xz_0_yy[i] = -2.0 * g_yy_xz_0_yy[i] * a_exp + 4.0 * g_xxyy_xz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xz_0_yz[i] = -2.0 * g_yy_xz_0_yz[i] * a_exp + 4.0 * g_xxyy_xz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_xz_0_zz[i] = -2.0 * g_yy_xz_0_zz[i] * a_exp + 4.0 * g_xxyy_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_xx_0_0_0_yy_yy_0_xx, g_xx_0_0_0_yy_yy_0_xy, g_xx_0_0_0_yy_yy_0_xz, g_xx_0_0_0_yy_yy_0_yy, g_xx_0_0_0_yy_yy_0_yz, g_xx_0_0_0_yy_yy_0_zz, g_xxyy_yy_0_xx, g_xxyy_yy_0_xy, g_xxyy_yy_0_xz, g_xxyy_yy_0_yy, g_xxyy_yy_0_yz, g_xxyy_yy_0_zz, g_yy_yy_0_xx, g_yy_yy_0_xy, g_yy_yy_0_xz, g_yy_yy_0_yy, g_yy_yy_0_yz, g_yy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_yy_0_xx[i] = -2.0 * g_yy_yy_0_xx[i] * a_exp + 4.0 * g_xxyy_yy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yy_0_xy[i] = -2.0 * g_yy_yy_0_xy[i] * a_exp + 4.0 * g_xxyy_yy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yy_0_xz[i] = -2.0 * g_yy_yy_0_xz[i] * a_exp + 4.0 * g_xxyy_yy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yy_0_yy[i] = -2.0 * g_yy_yy_0_yy[i] * a_exp + 4.0 * g_xxyy_yy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yy_0_yz[i] = -2.0 * g_yy_yy_0_yz[i] * a_exp + 4.0 * g_xxyy_yy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yy_0_zz[i] = -2.0 * g_yy_yy_0_zz[i] * a_exp + 4.0 * g_xxyy_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_xx_0_0_0_yy_yz_0_xx, g_xx_0_0_0_yy_yz_0_xy, g_xx_0_0_0_yy_yz_0_xz, g_xx_0_0_0_yy_yz_0_yy, g_xx_0_0_0_yy_yz_0_yz, g_xx_0_0_0_yy_yz_0_zz, g_xxyy_yz_0_xx, g_xxyy_yz_0_xy, g_xxyy_yz_0_xz, g_xxyy_yz_0_yy, g_xxyy_yz_0_yz, g_xxyy_yz_0_zz, g_yy_yz_0_xx, g_yy_yz_0_xy, g_yy_yz_0_xz, g_yy_yz_0_yy, g_yy_yz_0_yz, g_yy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_yz_0_xx[i] = -2.0 * g_yy_yz_0_xx[i] * a_exp + 4.0 * g_xxyy_yz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yz_0_xy[i] = -2.0 * g_yy_yz_0_xy[i] * a_exp + 4.0 * g_xxyy_yz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yz_0_xz[i] = -2.0 * g_yy_yz_0_xz[i] * a_exp + 4.0 * g_xxyy_yz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yz_0_yy[i] = -2.0 * g_yy_yz_0_yy[i] * a_exp + 4.0 * g_xxyy_yz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yz_0_yz[i] = -2.0 * g_yy_yz_0_yz[i] * a_exp + 4.0 * g_xxyy_yz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_yz_0_zz[i] = -2.0 * g_yy_yz_0_zz[i] * a_exp + 4.0 * g_xxyy_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_xx_0_0_0_yy_zz_0_xx, g_xx_0_0_0_yy_zz_0_xy, g_xx_0_0_0_yy_zz_0_xz, g_xx_0_0_0_yy_zz_0_yy, g_xx_0_0_0_yy_zz_0_yz, g_xx_0_0_0_yy_zz_0_zz, g_xxyy_zz_0_xx, g_xxyy_zz_0_xy, g_xxyy_zz_0_xz, g_xxyy_zz_0_yy, g_xxyy_zz_0_yz, g_xxyy_zz_0_zz, g_yy_zz_0_xx, g_yy_zz_0_xy, g_yy_zz_0_xz, g_yy_zz_0_yy, g_yy_zz_0_yz, g_yy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yy_zz_0_xx[i] = -2.0 * g_yy_zz_0_xx[i] * a_exp + 4.0 * g_xxyy_zz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_zz_0_xy[i] = -2.0 * g_yy_zz_0_xy[i] * a_exp + 4.0 * g_xxyy_zz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_zz_0_xz[i] = -2.0 * g_yy_zz_0_xz[i] * a_exp + 4.0 * g_xxyy_zz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_zz_0_yy[i] = -2.0 * g_yy_zz_0_yy[i] * a_exp + 4.0 * g_xxyy_zz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_zz_0_yz[i] = -2.0 * g_yy_zz_0_yz[i] * a_exp + 4.0 * g_xxyy_zz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_yy_zz_0_zz[i] = -2.0 * g_yy_zz_0_zz[i] * a_exp + 4.0 * g_xxyy_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_xx_0_0_0_yz_xx_0_xx, g_xx_0_0_0_yz_xx_0_xy, g_xx_0_0_0_yz_xx_0_xz, g_xx_0_0_0_yz_xx_0_yy, g_xx_0_0_0_yz_xx_0_yz, g_xx_0_0_0_yz_xx_0_zz, g_xxyz_xx_0_xx, g_xxyz_xx_0_xy, g_xxyz_xx_0_xz, g_xxyz_xx_0_yy, g_xxyz_xx_0_yz, g_xxyz_xx_0_zz, g_yz_xx_0_xx, g_yz_xx_0_xy, g_yz_xx_0_xz, g_yz_xx_0_yy, g_yz_xx_0_yz, g_yz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_xx_0_xx[i] = -2.0 * g_yz_xx_0_xx[i] * a_exp + 4.0 * g_xxyz_xx_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xx_0_xy[i] = -2.0 * g_yz_xx_0_xy[i] * a_exp + 4.0 * g_xxyz_xx_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xx_0_xz[i] = -2.0 * g_yz_xx_0_xz[i] * a_exp + 4.0 * g_xxyz_xx_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xx_0_yy[i] = -2.0 * g_yz_xx_0_yy[i] * a_exp + 4.0 * g_xxyz_xx_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xx_0_yz[i] = -2.0 * g_yz_xx_0_yz[i] * a_exp + 4.0 * g_xxyz_xx_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xx_0_zz[i] = -2.0 * g_yz_xx_0_zz[i] * a_exp + 4.0 * g_xxyz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_xx_0_0_0_yz_xy_0_xx, g_xx_0_0_0_yz_xy_0_xy, g_xx_0_0_0_yz_xy_0_xz, g_xx_0_0_0_yz_xy_0_yy, g_xx_0_0_0_yz_xy_0_yz, g_xx_0_0_0_yz_xy_0_zz, g_xxyz_xy_0_xx, g_xxyz_xy_0_xy, g_xxyz_xy_0_xz, g_xxyz_xy_0_yy, g_xxyz_xy_0_yz, g_xxyz_xy_0_zz, g_yz_xy_0_xx, g_yz_xy_0_xy, g_yz_xy_0_xz, g_yz_xy_0_yy, g_yz_xy_0_yz, g_yz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_xy_0_xx[i] = -2.0 * g_yz_xy_0_xx[i] * a_exp + 4.0 * g_xxyz_xy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xy_0_xy[i] = -2.0 * g_yz_xy_0_xy[i] * a_exp + 4.0 * g_xxyz_xy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xy_0_xz[i] = -2.0 * g_yz_xy_0_xz[i] * a_exp + 4.0 * g_xxyz_xy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xy_0_yy[i] = -2.0 * g_yz_xy_0_yy[i] * a_exp + 4.0 * g_xxyz_xy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xy_0_yz[i] = -2.0 * g_yz_xy_0_yz[i] * a_exp + 4.0 * g_xxyz_xy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xy_0_zz[i] = -2.0 * g_yz_xy_0_zz[i] * a_exp + 4.0 * g_xxyz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_xx_0_0_0_yz_xz_0_xx, g_xx_0_0_0_yz_xz_0_xy, g_xx_0_0_0_yz_xz_0_xz, g_xx_0_0_0_yz_xz_0_yy, g_xx_0_0_0_yz_xz_0_yz, g_xx_0_0_0_yz_xz_0_zz, g_xxyz_xz_0_xx, g_xxyz_xz_0_xy, g_xxyz_xz_0_xz, g_xxyz_xz_0_yy, g_xxyz_xz_0_yz, g_xxyz_xz_0_zz, g_yz_xz_0_xx, g_yz_xz_0_xy, g_yz_xz_0_xz, g_yz_xz_0_yy, g_yz_xz_0_yz, g_yz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_xz_0_xx[i] = -2.0 * g_yz_xz_0_xx[i] * a_exp + 4.0 * g_xxyz_xz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xz_0_xy[i] = -2.0 * g_yz_xz_0_xy[i] * a_exp + 4.0 * g_xxyz_xz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xz_0_xz[i] = -2.0 * g_yz_xz_0_xz[i] * a_exp + 4.0 * g_xxyz_xz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xz_0_yy[i] = -2.0 * g_yz_xz_0_yy[i] * a_exp + 4.0 * g_xxyz_xz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xz_0_yz[i] = -2.0 * g_yz_xz_0_yz[i] * a_exp + 4.0 * g_xxyz_xz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_xz_0_zz[i] = -2.0 * g_yz_xz_0_zz[i] * a_exp + 4.0 * g_xxyz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_xx_0_0_0_yz_yy_0_xx, g_xx_0_0_0_yz_yy_0_xy, g_xx_0_0_0_yz_yy_0_xz, g_xx_0_0_0_yz_yy_0_yy, g_xx_0_0_0_yz_yy_0_yz, g_xx_0_0_0_yz_yy_0_zz, g_xxyz_yy_0_xx, g_xxyz_yy_0_xy, g_xxyz_yy_0_xz, g_xxyz_yy_0_yy, g_xxyz_yy_0_yz, g_xxyz_yy_0_zz, g_yz_yy_0_xx, g_yz_yy_0_xy, g_yz_yy_0_xz, g_yz_yy_0_yy, g_yz_yy_0_yz, g_yz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_yy_0_xx[i] = -2.0 * g_yz_yy_0_xx[i] * a_exp + 4.0 * g_xxyz_yy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yy_0_xy[i] = -2.0 * g_yz_yy_0_xy[i] * a_exp + 4.0 * g_xxyz_yy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yy_0_xz[i] = -2.0 * g_yz_yy_0_xz[i] * a_exp + 4.0 * g_xxyz_yy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yy_0_yy[i] = -2.0 * g_yz_yy_0_yy[i] * a_exp + 4.0 * g_xxyz_yy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yy_0_yz[i] = -2.0 * g_yz_yy_0_yz[i] * a_exp + 4.0 * g_xxyz_yy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yy_0_zz[i] = -2.0 * g_yz_yy_0_zz[i] * a_exp + 4.0 * g_xxyz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_xx_0_0_0_yz_yz_0_xx, g_xx_0_0_0_yz_yz_0_xy, g_xx_0_0_0_yz_yz_0_xz, g_xx_0_0_0_yz_yz_0_yy, g_xx_0_0_0_yz_yz_0_yz, g_xx_0_0_0_yz_yz_0_zz, g_xxyz_yz_0_xx, g_xxyz_yz_0_xy, g_xxyz_yz_0_xz, g_xxyz_yz_0_yy, g_xxyz_yz_0_yz, g_xxyz_yz_0_zz, g_yz_yz_0_xx, g_yz_yz_0_xy, g_yz_yz_0_xz, g_yz_yz_0_yy, g_yz_yz_0_yz, g_yz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_yz_0_xx[i] = -2.0 * g_yz_yz_0_xx[i] * a_exp + 4.0 * g_xxyz_yz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yz_0_xy[i] = -2.0 * g_yz_yz_0_xy[i] * a_exp + 4.0 * g_xxyz_yz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yz_0_xz[i] = -2.0 * g_yz_yz_0_xz[i] * a_exp + 4.0 * g_xxyz_yz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yz_0_yy[i] = -2.0 * g_yz_yz_0_yy[i] * a_exp + 4.0 * g_xxyz_yz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yz_0_yz[i] = -2.0 * g_yz_yz_0_yz[i] * a_exp + 4.0 * g_xxyz_yz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_yz_0_zz[i] = -2.0 * g_yz_yz_0_zz[i] * a_exp + 4.0 * g_xxyz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_xx_0_0_0_yz_zz_0_xx, g_xx_0_0_0_yz_zz_0_xy, g_xx_0_0_0_yz_zz_0_xz, g_xx_0_0_0_yz_zz_0_yy, g_xx_0_0_0_yz_zz_0_yz, g_xx_0_0_0_yz_zz_0_zz, g_xxyz_zz_0_xx, g_xxyz_zz_0_xy, g_xxyz_zz_0_xz, g_xxyz_zz_0_yy, g_xxyz_zz_0_yz, g_xxyz_zz_0_zz, g_yz_zz_0_xx, g_yz_zz_0_xy, g_yz_zz_0_xz, g_yz_zz_0_yy, g_yz_zz_0_yz, g_yz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_yz_zz_0_xx[i] = -2.0 * g_yz_zz_0_xx[i] * a_exp + 4.0 * g_xxyz_zz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_zz_0_xy[i] = -2.0 * g_yz_zz_0_xy[i] * a_exp + 4.0 * g_xxyz_zz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_zz_0_xz[i] = -2.0 * g_yz_zz_0_xz[i] * a_exp + 4.0 * g_xxyz_zz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_zz_0_yy[i] = -2.0 * g_yz_zz_0_yy[i] * a_exp + 4.0 * g_xxyz_zz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_zz_0_yz[i] = -2.0 * g_yz_zz_0_yz[i] * a_exp + 4.0 * g_xxyz_zz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_yz_zz_0_zz[i] = -2.0 * g_yz_zz_0_zz[i] * a_exp + 4.0 * g_xxyz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_xx_0_0_0_zz_xx_0_xx, g_xx_0_0_0_zz_xx_0_xy, g_xx_0_0_0_zz_xx_0_xz, g_xx_0_0_0_zz_xx_0_yy, g_xx_0_0_0_zz_xx_0_yz, g_xx_0_0_0_zz_xx_0_zz, g_xxzz_xx_0_xx, g_xxzz_xx_0_xy, g_xxzz_xx_0_xz, g_xxzz_xx_0_yy, g_xxzz_xx_0_yz, g_xxzz_xx_0_zz, g_zz_xx_0_xx, g_zz_xx_0_xy, g_zz_xx_0_xz, g_zz_xx_0_yy, g_zz_xx_0_yz, g_zz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_xx_0_xx[i] = -2.0 * g_zz_xx_0_xx[i] * a_exp + 4.0 * g_xxzz_xx_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xx_0_xy[i] = -2.0 * g_zz_xx_0_xy[i] * a_exp + 4.0 * g_xxzz_xx_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xx_0_xz[i] = -2.0 * g_zz_xx_0_xz[i] * a_exp + 4.0 * g_xxzz_xx_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xx_0_yy[i] = -2.0 * g_zz_xx_0_yy[i] * a_exp + 4.0 * g_xxzz_xx_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xx_0_yz[i] = -2.0 * g_zz_xx_0_yz[i] * a_exp + 4.0 * g_xxzz_xx_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xx_0_zz[i] = -2.0 * g_zz_xx_0_zz[i] * a_exp + 4.0 * g_xxzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_xx_0_0_0_zz_xy_0_xx, g_xx_0_0_0_zz_xy_0_xy, g_xx_0_0_0_zz_xy_0_xz, g_xx_0_0_0_zz_xy_0_yy, g_xx_0_0_0_zz_xy_0_yz, g_xx_0_0_0_zz_xy_0_zz, g_xxzz_xy_0_xx, g_xxzz_xy_0_xy, g_xxzz_xy_0_xz, g_xxzz_xy_0_yy, g_xxzz_xy_0_yz, g_xxzz_xy_0_zz, g_zz_xy_0_xx, g_zz_xy_0_xy, g_zz_xy_0_xz, g_zz_xy_0_yy, g_zz_xy_0_yz, g_zz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_xy_0_xx[i] = -2.0 * g_zz_xy_0_xx[i] * a_exp + 4.0 * g_xxzz_xy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xy_0_xy[i] = -2.0 * g_zz_xy_0_xy[i] * a_exp + 4.0 * g_xxzz_xy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xy_0_xz[i] = -2.0 * g_zz_xy_0_xz[i] * a_exp + 4.0 * g_xxzz_xy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xy_0_yy[i] = -2.0 * g_zz_xy_0_yy[i] * a_exp + 4.0 * g_xxzz_xy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xy_0_yz[i] = -2.0 * g_zz_xy_0_yz[i] * a_exp + 4.0 * g_xxzz_xy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xy_0_zz[i] = -2.0 * g_zz_xy_0_zz[i] * a_exp + 4.0 * g_xxzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_xx_0_0_0_zz_xz_0_xx, g_xx_0_0_0_zz_xz_0_xy, g_xx_0_0_0_zz_xz_0_xz, g_xx_0_0_0_zz_xz_0_yy, g_xx_0_0_0_zz_xz_0_yz, g_xx_0_0_0_zz_xz_0_zz, g_xxzz_xz_0_xx, g_xxzz_xz_0_xy, g_xxzz_xz_0_xz, g_xxzz_xz_0_yy, g_xxzz_xz_0_yz, g_xxzz_xz_0_zz, g_zz_xz_0_xx, g_zz_xz_0_xy, g_zz_xz_0_xz, g_zz_xz_0_yy, g_zz_xz_0_yz, g_zz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_xz_0_xx[i] = -2.0 * g_zz_xz_0_xx[i] * a_exp + 4.0 * g_xxzz_xz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xz_0_xy[i] = -2.0 * g_zz_xz_0_xy[i] * a_exp + 4.0 * g_xxzz_xz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xz_0_xz[i] = -2.0 * g_zz_xz_0_xz[i] * a_exp + 4.0 * g_xxzz_xz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xz_0_yy[i] = -2.0 * g_zz_xz_0_yy[i] * a_exp + 4.0 * g_xxzz_xz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xz_0_yz[i] = -2.0 * g_zz_xz_0_yz[i] * a_exp + 4.0 * g_xxzz_xz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_xz_0_zz[i] = -2.0 * g_zz_xz_0_zz[i] * a_exp + 4.0 * g_xxzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_xx_0_0_0_zz_yy_0_xx, g_xx_0_0_0_zz_yy_0_xy, g_xx_0_0_0_zz_yy_0_xz, g_xx_0_0_0_zz_yy_0_yy, g_xx_0_0_0_zz_yy_0_yz, g_xx_0_0_0_zz_yy_0_zz, g_xxzz_yy_0_xx, g_xxzz_yy_0_xy, g_xxzz_yy_0_xz, g_xxzz_yy_0_yy, g_xxzz_yy_0_yz, g_xxzz_yy_0_zz, g_zz_yy_0_xx, g_zz_yy_0_xy, g_zz_yy_0_xz, g_zz_yy_0_yy, g_zz_yy_0_yz, g_zz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_yy_0_xx[i] = -2.0 * g_zz_yy_0_xx[i] * a_exp + 4.0 * g_xxzz_yy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yy_0_xy[i] = -2.0 * g_zz_yy_0_xy[i] * a_exp + 4.0 * g_xxzz_yy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yy_0_xz[i] = -2.0 * g_zz_yy_0_xz[i] * a_exp + 4.0 * g_xxzz_yy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yy_0_yy[i] = -2.0 * g_zz_yy_0_yy[i] * a_exp + 4.0 * g_xxzz_yy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yy_0_yz[i] = -2.0 * g_zz_yy_0_yz[i] * a_exp + 4.0 * g_xxzz_yy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yy_0_zz[i] = -2.0 * g_zz_yy_0_zz[i] * a_exp + 4.0 * g_xxzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_xx_0_0_0_zz_yz_0_xx, g_xx_0_0_0_zz_yz_0_xy, g_xx_0_0_0_zz_yz_0_xz, g_xx_0_0_0_zz_yz_0_yy, g_xx_0_0_0_zz_yz_0_yz, g_xx_0_0_0_zz_yz_0_zz, g_xxzz_yz_0_xx, g_xxzz_yz_0_xy, g_xxzz_yz_0_xz, g_xxzz_yz_0_yy, g_xxzz_yz_0_yz, g_xxzz_yz_0_zz, g_zz_yz_0_xx, g_zz_yz_0_xy, g_zz_yz_0_xz, g_zz_yz_0_yy, g_zz_yz_0_yz, g_zz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_yz_0_xx[i] = -2.0 * g_zz_yz_0_xx[i] * a_exp + 4.0 * g_xxzz_yz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yz_0_xy[i] = -2.0 * g_zz_yz_0_xy[i] * a_exp + 4.0 * g_xxzz_yz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yz_0_xz[i] = -2.0 * g_zz_yz_0_xz[i] * a_exp + 4.0 * g_xxzz_yz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yz_0_yy[i] = -2.0 * g_zz_yz_0_yy[i] * a_exp + 4.0 * g_xxzz_yz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yz_0_yz[i] = -2.0 * g_zz_yz_0_yz[i] * a_exp + 4.0 * g_xxzz_yz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_yz_0_zz[i] = -2.0 * g_zz_yz_0_zz[i] * a_exp + 4.0 * g_xxzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_xx_0_0_0_zz_zz_0_xx, g_xx_0_0_0_zz_zz_0_xy, g_xx_0_0_0_zz_zz_0_xz, g_xx_0_0_0_zz_zz_0_yy, g_xx_0_0_0_zz_zz_0_yz, g_xx_0_0_0_zz_zz_0_zz, g_xxzz_zz_0_xx, g_xxzz_zz_0_xy, g_xxzz_zz_0_xz, g_xxzz_zz_0_yy, g_xxzz_zz_0_yz, g_xxzz_zz_0_zz, g_zz_zz_0_xx, g_zz_zz_0_xy, g_zz_zz_0_xz, g_zz_zz_0_yy, g_zz_zz_0_yz, g_zz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_zz_zz_0_xx[i] = -2.0 * g_zz_zz_0_xx[i] * a_exp + 4.0 * g_xxzz_zz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_zz_0_xy[i] = -2.0 * g_zz_zz_0_xy[i] * a_exp + 4.0 * g_xxzz_zz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_zz_0_xz[i] = -2.0 * g_zz_zz_0_xz[i] * a_exp + 4.0 * g_xxzz_zz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_zz_0_yy[i] = -2.0 * g_zz_zz_0_yy[i] * a_exp + 4.0 * g_xxzz_zz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_zz_0_yz[i] = -2.0 * g_zz_zz_0_yz[i] * a_exp + 4.0 * g_xxzz_zz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_zz_zz_0_zz[i] = -2.0 * g_zz_zz_0_zz[i] * a_exp + 4.0 * g_xxzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_xxxy_xx_0_xx, g_xxxy_xx_0_xy, g_xxxy_xx_0_xz, g_xxxy_xx_0_yy, g_xxxy_xx_0_yz, g_xxxy_xx_0_zz, g_xy_0_0_0_xx_xx_0_xx, g_xy_0_0_0_xx_xx_0_xy, g_xy_0_0_0_xx_xx_0_xz, g_xy_0_0_0_xx_xx_0_yy, g_xy_0_0_0_xx_xx_0_yz, g_xy_0_0_0_xx_xx_0_zz, g_xy_xx_0_xx, g_xy_xx_0_xy, g_xy_xx_0_xz, g_xy_xx_0_yy, g_xy_xx_0_yz, g_xy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_xx_0_xx[i] = -4.0 * g_xy_xx_0_xx[i] * a_exp + 4.0 * g_xxxy_xx_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xx_0_xy[i] = -4.0 * g_xy_xx_0_xy[i] * a_exp + 4.0 * g_xxxy_xx_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xx_0_xz[i] = -4.0 * g_xy_xx_0_xz[i] * a_exp + 4.0 * g_xxxy_xx_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xx_0_yy[i] = -4.0 * g_xy_xx_0_yy[i] * a_exp + 4.0 * g_xxxy_xx_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xx_0_yz[i] = -4.0 * g_xy_xx_0_yz[i] * a_exp + 4.0 * g_xxxy_xx_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xx_0_zz[i] = -4.0 * g_xy_xx_0_zz[i] * a_exp + 4.0 * g_xxxy_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_xxxy_xy_0_xx, g_xxxy_xy_0_xy, g_xxxy_xy_0_xz, g_xxxy_xy_0_yy, g_xxxy_xy_0_yz, g_xxxy_xy_0_zz, g_xy_0_0_0_xx_xy_0_xx, g_xy_0_0_0_xx_xy_0_xy, g_xy_0_0_0_xx_xy_0_xz, g_xy_0_0_0_xx_xy_0_yy, g_xy_0_0_0_xx_xy_0_yz, g_xy_0_0_0_xx_xy_0_zz, g_xy_xy_0_xx, g_xy_xy_0_xy, g_xy_xy_0_xz, g_xy_xy_0_yy, g_xy_xy_0_yz, g_xy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_xy_0_xx[i] = -4.0 * g_xy_xy_0_xx[i] * a_exp + 4.0 * g_xxxy_xy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xy_0_xy[i] = -4.0 * g_xy_xy_0_xy[i] * a_exp + 4.0 * g_xxxy_xy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xy_0_xz[i] = -4.0 * g_xy_xy_0_xz[i] * a_exp + 4.0 * g_xxxy_xy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xy_0_yy[i] = -4.0 * g_xy_xy_0_yy[i] * a_exp + 4.0 * g_xxxy_xy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xy_0_yz[i] = -4.0 * g_xy_xy_0_yz[i] * a_exp + 4.0 * g_xxxy_xy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xy_0_zz[i] = -4.0 * g_xy_xy_0_zz[i] * a_exp + 4.0 * g_xxxy_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_xxxy_xz_0_xx, g_xxxy_xz_0_xy, g_xxxy_xz_0_xz, g_xxxy_xz_0_yy, g_xxxy_xz_0_yz, g_xxxy_xz_0_zz, g_xy_0_0_0_xx_xz_0_xx, g_xy_0_0_0_xx_xz_0_xy, g_xy_0_0_0_xx_xz_0_xz, g_xy_0_0_0_xx_xz_0_yy, g_xy_0_0_0_xx_xz_0_yz, g_xy_0_0_0_xx_xz_0_zz, g_xy_xz_0_xx, g_xy_xz_0_xy, g_xy_xz_0_xz, g_xy_xz_0_yy, g_xy_xz_0_yz, g_xy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_xz_0_xx[i] = -4.0 * g_xy_xz_0_xx[i] * a_exp + 4.0 * g_xxxy_xz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xz_0_xy[i] = -4.0 * g_xy_xz_0_xy[i] * a_exp + 4.0 * g_xxxy_xz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xz_0_xz[i] = -4.0 * g_xy_xz_0_xz[i] * a_exp + 4.0 * g_xxxy_xz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xz_0_yy[i] = -4.0 * g_xy_xz_0_yy[i] * a_exp + 4.0 * g_xxxy_xz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xz_0_yz[i] = -4.0 * g_xy_xz_0_yz[i] * a_exp + 4.0 * g_xxxy_xz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_xz_0_zz[i] = -4.0 * g_xy_xz_0_zz[i] * a_exp + 4.0 * g_xxxy_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_xxxy_yy_0_xx, g_xxxy_yy_0_xy, g_xxxy_yy_0_xz, g_xxxy_yy_0_yy, g_xxxy_yy_0_yz, g_xxxy_yy_0_zz, g_xy_0_0_0_xx_yy_0_xx, g_xy_0_0_0_xx_yy_0_xy, g_xy_0_0_0_xx_yy_0_xz, g_xy_0_0_0_xx_yy_0_yy, g_xy_0_0_0_xx_yy_0_yz, g_xy_0_0_0_xx_yy_0_zz, g_xy_yy_0_xx, g_xy_yy_0_xy, g_xy_yy_0_xz, g_xy_yy_0_yy, g_xy_yy_0_yz, g_xy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_yy_0_xx[i] = -4.0 * g_xy_yy_0_xx[i] * a_exp + 4.0 * g_xxxy_yy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yy_0_xy[i] = -4.0 * g_xy_yy_0_xy[i] * a_exp + 4.0 * g_xxxy_yy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yy_0_xz[i] = -4.0 * g_xy_yy_0_xz[i] * a_exp + 4.0 * g_xxxy_yy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yy_0_yy[i] = -4.0 * g_xy_yy_0_yy[i] * a_exp + 4.0 * g_xxxy_yy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yy_0_yz[i] = -4.0 * g_xy_yy_0_yz[i] * a_exp + 4.0 * g_xxxy_yy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yy_0_zz[i] = -4.0 * g_xy_yy_0_zz[i] * a_exp + 4.0 * g_xxxy_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_xxxy_yz_0_xx, g_xxxy_yz_0_xy, g_xxxy_yz_0_xz, g_xxxy_yz_0_yy, g_xxxy_yz_0_yz, g_xxxy_yz_0_zz, g_xy_0_0_0_xx_yz_0_xx, g_xy_0_0_0_xx_yz_0_xy, g_xy_0_0_0_xx_yz_0_xz, g_xy_0_0_0_xx_yz_0_yy, g_xy_0_0_0_xx_yz_0_yz, g_xy_0_0_0_xx_yz_0_zz, g_xy_yz_0_xx, g_xy_yz_0_xy, g_xy_yz_0_xz, g_xy_yz_0_yy, g_xy_yz_0_yz, g_xy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_yz_0_xx[i] = -4.0 * g_xy_yz_0_xx[i] * a_exp + 4.0 * g_xxxy_yz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yz_0_xy[i] = -4.0 * g_xy_yz_0_xy[i] * a_exp + 4.0 * g_xxxy_yz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yz_0_xz[i] = -4.0 * g_xy_yz_0_xz[i] * a_exp + 4.0 * g_xxxy_yz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yz_0_yy[i] = -4.0 * g_xy_yz_0_yy[i] * a_exp + 4.0 * g_xxxy_yz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yz_0_yz[i] = -4.0 * g_xy_yz_0_yz[i] * a_exp + 4.0 * g_xxxy_yz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_yz_0_zz[i] = -4.0 * g_xy_yz_0_zz[i] * a_exp + 4.0 * g_xxxy_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_xxxy_zz_0_xx, g_xxxy_zz_0_xy, g_xxxy_zz_0_xz, g_xxxy_zz_0_yy, g_xxxy_zz_0_yz, g_xxxy_zz_0_zz, g_xy_0_0_0_xx_zz_0_xx, g_xy_0_0_0_xx_zz_0_xy, g_xy_0_0_0_xx_zz_0_xz, g_xy_0_0_0_xx_zz_0_yy, g_xy_0_0_0_xx_zz_0_yz, g_xy_0_0_0_xx_zz_0_zz, g_xy_zz_0_xx, g_xy_zz_0_xy, g_xy_zz_0_xz, g_xy_zz_0_yy, g_xy_zz_0_yz, g_xy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xx_zz_0_xx[i] = -4.0 * g_xy_zz_0_xx[i] * a_exp + 4.0 * g_xxxy_zz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_zz_0_xy[i] = -4.0 * g_xy_zz_0_xy[i] * a_exp + 4.0 * g_xxxy_zz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_zz_0_xz[i] = -4.0 * g_xy_zz_0_xz[i] * a_exp + 4.0 * g_xxxy_zz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_zz_0_yy[i] = -4.0 * g_xy_zz_0_yy[i] * a_exp + 4.0 * g_xxxy_zz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_zz_0_yz[i] = -4.0 * g_xy_zz_0_yz[i] * a_exp + 4.0 * g_xxxy_zz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xx_zz_0_zz[i] = -4.0 * g_xy_zz_0_zz[i] * a_exp + 4.0 * g_xxxy_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_0_xx_0_xx, g_0_xx_0_xy, g_0_xx_0_xz, g_0_xx_0_yy, g_0_xx_0_yz, g_0_xx_0_zz, g_xx_xx_0_xx, g_xx_xx_0_xy, g_xx_xx_0_xz, g_xx_xx_0_yy, g_xx_xx_0_yz, g_xx_xx_0_zz, g_xxyy_xx_0_xx, g_xxyy_xx_0_xy, g_xxyy_xx_0_xz, g_xxyy_xx_0_yy, g_xxyy_xx_0_yz, g_xxyy_xx_0_zz, g_xy_0_0_0_xy_xx_0_xx, g_xy_0_0_0_xy_xx_0_xy, g_xy_0_0_0_xy_xx_0_xz, g_xy_0_0_0_xy_xx_0_yy, g_xy_0_0_0_xy_xx_0_yz, g_xy_0_0_0_xy_xx_0_zz, g_yy_xx_0_xx, g_yy_xx_0_xy, g_yy_xx_0_xz, g_yy_xx_0_yy, g_yy_xx_0_yz, g_yy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_xx_0_xx[i] = g_0_xx_0_xx[i] - 2.0 * g_yy_xx_0_xx[i] * a_exp - 2.0 * g_xx_xx_0_xx[i] * a_exp + 4.0 * g_xxyy_xx_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xx_0_xy[i] = g_0_xx_0_xy[i] - 2.0 * g_yy_xx_0_xy[i] * a_exp - 2.0 * g_xx_xx_0_xy[i] * a_exp + 4.0 * g_xxyy_xx_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xx_0_xz[i] = g_0_xx_0_xz[i] - 2.0 * g_yy_xx_0_xz[i] * a_exp - 2.0 * g_xx_xx_0_xz[i] * a_exp + 4.0 * g_xxyy_xx_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xx_0_yy[i] = g_0_xx_0_yy[i] - 2.0 * g_yy_xx_0_yy[i] * a_exp - 2.0 * g_xx_xx_0_yy[i] * a_exp + 4.0 * g_xxyy_xx_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xx_0_yz[i] = g_0_xx_0_yz[i] - 2.0 * g_yy_xx_0_yz[i] * a_exp - 2.0 * g_xx_xx_0_yz[i] * a_exp + 4.0 * g_xxyy_xx_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xx_0_zz[i] = g_0_xx_0_zz[i] - 2.0 * g_yy_xx_0_zz[i] * a_exp - 2.0 * g_xx_xx_0_zz[i] * a_exp + 4.0 * g_xxyy_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_xx_xy_0_xx, g_xx_xy_0_xy, g_xx_xy_0_xz, g_xx_xy_0_yy, g_xx_xy_0_yz, g_xx_xy_0_zz, g_xxyy_xy_0_xx, g_xxyy_xy_0_xy, g_xxyy_xy_0_xz, g_xxyy_xy_0_yy, g_xxyy_xy_0_yz, g_xxyy_xy_0_zz, g_xy_0_0_0_xy_xy_0_xx, g_xy_0_0_0_xy_xy_0_xy, g_xy_0_0_0_xy_xy_0_xz, g_xy_0_0_0_xy_xy_0_yy, g_xy_0_0_0_xy_xy_0_yz, g_xy_0_0_0_xy_xy_0_zz, g_yy_xy_0_xx, g_yy_xy_0_xy, g_yy_xy_0_xz, g_yy_xy_0_yy, g_yy_xy_0_yz, g_yy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_xy_0_xx[i] = g_0_xy_0_xx[i] - 2.0 * g_yy_xy_0_xx[i] * a_exp - 2.0 * g_xx_xy_0_xx[i] * a_exp + 4.0 * g_xxyy_xy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xy_0_xy[i] = g_0_xy_0_xy[i] - 2.0 * g_yy_xy_0_xy[i] * a_exp - 2.0 * g_xx_xy_0_xy[i] * a_exp + 4.0 * g_xxyy_xy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xy_0_xz[i] = g_0_xy_0_xz[i] - 2.0 * g_yy_xy_0_xz[i] * a_exp - 2.0 * g_xx_xy_0_xz[i] * a_exp + 4.0 * g_xxyy_xy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xy_0_yy[i] = g_0_xy_0_yy[i] - 2.0 * g_yy_xy_0_yy[i] * a_exp - 2.0 * g_xx_xy_0_yy[i] * a_exp + 4.0 * g_xxyy_xy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xy_0_yz[i] = g_0_xy_0_yz[i] - 2.0 * g_yy_xy_0_yz[i] * a_exp - 2.0 * g_xx_xy_0_yz[i] * a_exp + 4.0 * g_xxyy_xy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xy_0_zz[i] = g_0_xy_0_zz[i] - 2.0 * g_yy_xy_0_zz[i] * a_exp - 2.0 * g_xx_xy_0_zz[i] * a_exp + 4.0 * g_xxyy_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_xx_xz_0_xx, g_xx_xz_0_xy, g_xx_xz_0_xz, g_xx_xz_0_yy, g_xx_xz_0_yz, g_xx_xz_0_zz, g_xxyy_xz_0_xx, g_xxyy_xz_0_xy, g_xxyy_xz_0_xz, g_xxyy_xz_0_yy, g_xxyy_xz_0_yz, g_xxyy_xz_0_zz, g_xy_0_0_0_xy_xz_0_xx, g_xy_0_0_0_xy_xz_0_xy, g_xy_0_0_0_xy_xz_0_xz, g_xy_0_0_0_xy_xz_0_yy, g_xy_0_0_0_xy_xz_0_yz, g_xy_0_0_0_xy_xz_0_zz, g_yy_xz_0_xx, g_yy_xz_0_xy, g_yy_xz_0_xz, g_yy_xz_0_yy, g_yy_xz_0_yz, g_yy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_xz_0_xx[i] = g_0_xz_0_xx[i] - 2.0 * g_yy_xz_0_xx[i] * a_exp - 2.0 * g_xx_xz_0_xx[i] * a_exp + 4.0 * g_xxyy_xz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xz_0_xy[i] = g_0_xz_0_xy[i] - 2.0 * g_yy_xz_0_xy[i] * a_exp - 2.0 * g_xx_xz_0_xy[i] * a_exp + 4.0 * g_xxyy_xz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xz_0_xz[i] = g_0_xz_0_xz[i] - 2.0 * g_yy_xz_0_xz[i] * a_exp - 2.0 * g_xx_xz_0_xz[i] * a_exp + 4.0 * g_xxyy_xz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xz_0_yy[i] = g_0_xz_0_yy[i] - 2.0 * g_yy_xz_0_yy[i] * a_exp - 2.0 * g_xx_xz_0_yy[i] * a_exp + 4.0 * g_xxyy_xz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xz_0_yz[i] = g_0_xz_0_yz[i] - 2.0 * g_yy_xz_0_yz[i] * a_exp - 2.0 * g_xx_xz_0_yz[i] * a_exp + 4.0 * g_xxyy_xz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_xz_0_zz[i] = g_0_xz_0_zz[i] - 2.0 * g_yy_xz_0_zz[i] * a_exp - 2.0 * g_xx_xz_0_zz[i] * a_exp + 4.0 * g_xxyy_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_0_yy_0_xx, g_0_yy_0_xy, g_0_yy_0_xz, g_0_yy_0_yy, g_0_yy_0_yz, g_0_yy_0_zz, g_xx_yy_0_xx, g_xx_yy_0_xy, g_xx_yy_0_xz, g_xx_yy_0_yy, g_xx_yy_0_yz, g_xx_yy_0_zz, g_xxyy_yy_0_xx, g_xxyy_yy_0_xy, g_xxyy_yy_0_xz, g_xxyy_yy_0_yy, g_xxyy_yy_0_yz, g_xxyy_yy_0_zz, g_xy_0_0_0_xy_yy_0_xx, g_xy_0_0_0_xy_yy_0_xy, g_xy_0_0_0_xy_yy_0_xz, g_xy_0_0_0_xy_yy_0_yy, g_xy_0_0_0_xy_yy_0_yz, g_xy_0_0_0_xy_yy_0_zz, g_yy_yy_0_xx, g_yy_yy_0_xy, g_yy_yy_0_xz, g_yy_yy_0_yy, g_yy_yy_0_yz, g_yy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_yy_0_xx[i] = g_0_yy_0_xx[i] - 2.0 * g_yy_yy_0_xx[i] * a_exp - 2.0 * g_xx_yy_0_xx[i] * a_exp + 4.0 * g_xxyy_yy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yy_0_xy[i] = g_0_yy_0_xy[i] - 2.0 * g_yy_yy_0_xy[i] * a_exp - 2.0 * g_xx_yy_0_xy[i] * a_exp + 4.0 * g_xxyy_yy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yy_0_xz[i] = g_0_yy_0_xz[i] - 2.0 * g_yy_yy_0_xz[i] * a_exp - 2.0 * g_xx_yy_0_xz[i] * a_exp + 4.0 * g_xxyy_yy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yy_0_yy[i] = g_0_yy_0_yy[i] - 2.0 * g_yy_yy_0_yy[i] * a_exp - 2.0 * g_xx_yy_0_yy[i] * a_exp + 4.0 * g_xxyy_yy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yy_0_yz[i] = g_0_yy_0_yz[i] - 2.0 * g_yy_yy_0_yz[i] * a_exp - 2.0 * g_xx_yy_0_yz[i] * a_exp + 4.0 * g_xxyy_yy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yy_0_zz[i] = g_0_yy_0_zz[i] - 2.0 * g_yy_yy_0_zz[i] * a_exp - 2.0 * g_xx_yy_0_zz[i] * a_exp + 4.0 * g_xxyy_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_xx_yz_0_xx, g_xx_yz_0_xy, g_xx_yz_0_xz, g_xx_yz_0_yy, g_xx_yz_0_yz, g_xx_yz_0_zz, g_xxyy_yz_0_xx, g_xxyy_yz_0_xy, g_xxyy_yz_0_xz, g_xxyy_yz_0_yy, g_xxyy_yz_0_yz, g_xxyy_yz_0_zz, g_xy_0_0_0_xy_yz_0_xx, g_xy_0_0_0_xy_yz_0_xy, g_xy_0_0_0_xy_yz_0_xz, g_xy_0_0_0_xy_yz_0_yy, g_xy_0_0_0_xy_yz_0_yz, g_xy_0_0_0_xy_yz_0_zz, g_yy_yz_0_xx, g_yy_yz_0_xy, g_yy_yz_0_xz, g_yy_yz_0_yy, g_yy_yz_0_yz, g_yy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_yz_0_xx[i] = g_0_yz_0_xx[i] - 2.0 * g_yy_yz_0_xx[i] * a_exp - 2.0 * g_xx_yz_0_xx[i] * a_exp + 4.0 * g_xxyy_yz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yz_0_xy[i] = g_0_yz_0_xy[i] - 2.0 * g_yy_yz_0_xy[i] * a_exp - 2.0 * g_xx_yz_0_xy[i] * a_exp + 4.0 * g_xxyy_yz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yz_0_xz[i] = g_0_yz_0_xz[i] - 2.0 * g_yy_yz_0_xz[i] * a_exp - 2.0 * g_xx_yz_0_xz[i] * a_exp + 4.0 * g_xxyy_yz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yz_0_yy[i] = g_0_yz_0_yy[i] - 2.0 * g_yy_yz_0_yy[i] * a_exp - 2.0 * g_xx_yz_0_yy[i] * a_exp + 4.0 * g_xxyy_yz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yz_0_yz[i] = g_0_yz_0_yz[i] - 2.0 * g_yy_yz_0_yz[i] * a_exp - 2.0 * g_xx_yz_0_yz[i] * a_exp + 4.0 * g_xxyy_yz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_yz_0_zz[i] = g_0_yz_0_zz[i] - 2.0 * g_yy_yz_0_zz[i] * a_exp - 2.0 * g_xx_yz_0_zz[i] * a_exp + 4.0 * g_xxyy_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_0_zz_0_xx, g_0_zz_0_xy, g_0_zz_0_xz, g_0_zz_0_yy, g_0_zz_0_yz, g_0_zz_0_zz, g_xx_zz_0_xx, g_xx_zz_0_xy, g_xx_zz_0_xz, g_xx_zz_0_yy, g_xx_zz_0_yz, g_xx_zz_0_zz, g_xxyy_zz_0_xx, g_xxyy_zz_0_xy, g_xxyy_zz_0_xz, g_xxyy_zz_0_yy, g_xxyy_zz_0_yz, g_xxyy_zz_0_zz, g_xy_0_0_0_xy_zz_0_xx, g_xy_0_0_0_xy_zz_0_xy, g_xy_0_0_0_xy_zz_0_xz, g_xy_0_0_0_xy_zz_0_yy, g_xy_0_0_0_xy_zz_0_yz, g_xy_0_0_0_xy_zz_0_zz, g_yy_zz_0_xx, g_yy_zz_0_xy, g_yy_zz_0_xz, g_yy_zz_0_yy, g_yy_zz_0_yz, g_yy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xy_zz_0_xx[i] = g_0_zz_0_xx[i] - 2.0 * g_yy_zz_0_xx[i] * a_exp - 2.0 * g_xx_zz_0_xx[i] * a_exp + 4.0 * g_xxyy_zz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_zz_0_xy[i] = g_0_zz_0_xy[i] - 2.0 * g_yy_zz_0_xy[i] * a_exp - 2.0 * g_xx_zz_0_xy[i] * a_exp + 4.0 * g_xxyy_zz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_zz_0_xz[i] = g_0_zz_0_xz[i] - 2.0 * g_yy_zz_0_xz[i] * a_exp - 2.0 * g_xx_zz_0_xz[i] * a_exp + 4.0 * g_xxyy_zz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_zz_0_yy[i] = g_0_zz_0_yy[i] - 2.0 * g_yy_zz_0_yy[i] * a_exp - 2.0 * g_xx_zz_0_yy[i] * a_exp + 4.0 * g_xxyy_zz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_zz_0_yz[i] = g_0_zz_0_yz[i] - 2.0 * g_yy_zz_0_yz[i] * a_exp - 2.0 * g_xx_zz_0_yz[i] * a_exp + 4.0 * g_xxyy_zz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xy_zz_0_zz[i] = g_0_zz_0_zz[i] - 2.0 * g_yy_zz_0_zz[i] * a_exp - 2.0 * g_xx_zz_0_zz[i] * a_exp + 4.0 * g_xxyy_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_xxyz_xx_0_xx, g_xxyz_xx_0_xy, g_xxyz_xx_0_xz, g_xxyz_xx_0_yy, g_xxyz_xx_0_yz, g_xxyz_xx_0_zz, g_xy_0_0_0_xz_xx_0_xx, g_xy_0_0_0_xz_xx_0_xy, g_xy_0_0_0_xz_xx_0_xz, g_xy_0_0_0_xz_xx_0_yy, g_xy_0_0_0_xz_xx_0_yz, g_xy_0_0_0_xz_xx_0_zz, g_yz_xx_0_xx, g_yz_xx_0_xy, g_yz_xx_0_xz, g_yz_xx_0_yy, g_yz_xx_0_yz, g_yz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_xx_0_xx[i] = -2.0 * g_yz_xx_0_xx[i] * a_exp + 4.0 * g_xxyz_xx_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xx_0_xy[i] = -2.0 * g_yz_xx_0_xy[i] * a_exp + 4.0 * g_xxyz_xx_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xx_0_xz[i] = -2.0 * g_yz_xx_0_xz[i] * a_exp + 4.0 * g_xxyz_xx_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xx_0_yy[i] = -2.0 * g_yz_xx_0_yy[i] * a_exp + 4.0 * g_xxyz_xx_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xx_0_yz[i] = -2.0 * g_yz_xx_0_yz[i] * a_exp + 4.0 * g_xxyz_xx_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xx_0_zz[i] = -2.0 * g_yz_xx_0_zz[i] * a_exp + 4.0 * g_xxyz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_xxyz_xy_0_xx, g_xxyz_xy_0_xy, g_xxyz_xy_0_xz, g_xxyz_xy_0_yy, g_xxyz_xy_0_yz, g_xxyz_xy_0_zz, g_xy_0_0_0_xz_xy_0_xx, g_xy_0_0_0_xz_xy_0_xy, g_xy_0_0_0_xz_xy_0_xz, g_xy_0_0_0_xz_xy_0_yy, g_xy_0_0_0_xz_xy_0_yz, g_xy_0_0_0_xz_xy_0_zz, g_yz_xy_0_xx, g_yz_xy_0_xy, g_yz_xy_0_xz, g_yz_xy_0_yy, g_yz_xy_0_yz, g_yz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_xy_0_xx[i] = -2.0 * g_yz_xy_0_xx[i] * a_exp + 4.0 * g_xxyz_xy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xy_0_xy[i] = -2.0 * g_yz_xy_0_xy[i] * a_exp + 4.0 * g_xxyz_xy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xy_0_xz[i] = -2.0 * g_yz_xy_0_xz[i] * a_exp + 4.0 * g_xxyz_xy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xy_0_yy[i] = -2.0 * g_yz_xy_0_yy[i] * a_exp + 4.0 * g_xxyz_xy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xy_0_yz[i] = -2.0 * g_yz_xy_0_yz[i] * a_exp + 4.0 * g_xxyz_xy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xy_0_zz[i] = -2.0 * g_yz_xy_0_zz[i] * a_exp + 4.0 * g_xxyz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_xxyz_xz_0_xx, g_xxyz_xz_0_xy, g_xxyz_xz_0_xz, g_xxyz_xz_0_yy, g_xxyz_xz_0_yz, g_xxyz_xz_0_zz, g_xy_0_0_0_xz_xz_0_xx, g_xy_0_0_0_xz_xz_0_xy, g_xy_0_0_0_xz_xz_0_xz, g_xy_0_0_0_xz_xz_0_yy, g_xy_0_0_0_xz_xz_0_yz, g_xy_0_0_0_xz_xz_0_zz, g_yz_xz_0_xx, g_yz_xz_0_xy, g_yz_xz_0_xz, g_yz_xz_0_yy, g_yz_xz_0_yz, g_yz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_xz_0_xx[i] = -2.0 * g_yz_xz_0_xx[i] * a_exp + 4.0 * g_xxyz_xz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xz_0_xy[i] = -2.0 * g_yz_xz_0_xy[i] * a_exp + 4.0 * g_xxyz_xz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xz_0_xz[i] = -2.0 * g_yz_xz_0_xz[i] * a_exp + 4.0 * g_xxyz_xz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xz_0_yy[i] = -2.0 * g_yz_xz_0_yy[i] * a_exp + 4.0 * g_xxyz_xz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xz_0_yz[i] = -2.0 * g_yz_xz_0_yz[i] * a_exp + 4.0 * g_xxyz_xz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_xz_0_zz[i] = -2.0 * g_yz_xz_0_zz[i] * a_exp + 4.0 * g_xxyz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_xxyz_yy_0_xx, g_xxyz_yy_0_xy, g_xxyz_yy_0_xz, g_xxyz_yy_0_yy, g_xxyz_yy_0_yz, g_xxyz_yy_0_zz, g_xy_0_0_0_xz_yy_0_xx, g_xy_0_0_0_xz_yy_0_xy, g_xy_0_0_0_xz_yy_0_xz, g_xy_0_0_0_xz_yy_0_yy, g_xy_0_0_0_xz_yy_0_yz, g_xy_0_0_0_xz_yy_0_zz, g_yz_yy_0_xx, g_yz_yy_0_xy, g_yz_yy_0_xz, g_yz_yy_0_yy, g_yz_yy_0_yz, g_yz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_yy_0_xx[i] = -2.0 * g_yz_yy_0_xx[i] * a_exp + 4.0 * g_xxyz_yy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yy_0_xy[i] = -2.0 * g_yz_yy_0_xy[i] * a_exp + 4.0 * g_xxyz_yy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yy_0_xz[i] = -2.0 * g_yz_yy_0_xz[i] * a_exp + 4.0 * g_xxyz_yy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yy_0_yy[i] = -2.0 * g_yz_yy_0_yy[i] * a_exp + 4.0 * g_xxyz_yy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yy_0_yz[i] = -2.0 * g_yz_yy_0_yz[i] * a_exp + 4.0 * g_xxyz_yy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yy_0_zz[i] = -2.0 * g_yz_yy_0_zz[i] * a_exp + 4.0 * g_xxyz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_xxyz_yz_0_xx, g_xxyz_yz_0_xy, g_xxyz_yz_0_xz, g_xxyz_yz_0_yy, g_xxyz_yz_0_yz, g_xxyz_yz_0_zz, g_xy_0_0_0_xz_yz_0_xx, g_xy_0_0_0_xz_yz_0_xy, g_xy_0_0_0_xz_yz_0_xz, g_xy_0_0_0_xz_yz_0_yy, g_xy_0_0_0_xz_yz_0_yz, g_xy_0_0_0_xz_yz_0_zz, g_yz_yz_0_xx, g_yz_yz_0_xy, g_yz_yz_0_xz, g_yz_yz_0_yy, g_yz_yz_0_yz, g_yz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_yz_0_xx[i] = -2.0 * g_yz_yz_0_xx[i] * a_exp + 4.0 * g_xxyz_yz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yz_0_xy[i] = -2.0 * g_yz_yz_0_xy[i] * a_exp + 4.0 * g_xxyz_yz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yz_0_xz[i] = -2.0 * g_yz_yz_0_xz[i] * a_exp + 4.0 * g_xxyz_yz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yz_0_yy[i] = -2.0 * g_yz_yz_0_yy[i] * a_exp + 4.0 * g_xxyz_yz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yz_0_yz[i] = -2.0 * g_yz_yz_0_yz[i] * a_exp + 4.0 * g_xxyz_yz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_yz_0_zz[i] = -2.0 * g_yz_yz_0_zz[i] * a_exp + 4.0 * g_xxyz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_xxyz_zz_0_xx, g_xxyz_zz_0_xy, g_xxyz_zz_0_xz, g_xxyz_zz_0_yy, g_xxyz_zz_0_yz, g_xxyz_zz_0_zz, g_xy_0_0_0_xz_zz_0_xx, g_xy_0_0_0_xz_zz_0_xy, g_xy_0_0_0_xz_zz_0_xz, g_xy_0_0_0_xz_zz_0_yy, g_xy_0_0_0_xz_zz_0_yz, g_xy_0_0_0_xz_zz_0_zz, g_yz_zz_0_xx, g_yz_zz_0_xy, g_yz_zz_0_xz, g_yz_zz_0_yy, g_yz_zz_0_yz, g_yz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_xz_zz_0_xx[i] = -2.0 * g_yz_zz_0_xx[i] * a_exp + 4.0 * g_xxyz_zz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_zz_0_xy[i] = -2.0 * g_yz_zz_0_xy[i] * a_exp + 4.0 * g_xxyz_zz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_zz_0_xz[i] = -2.0 * g_yz_zz_0_xz[i] * a_exp + 4.0 * g_xxyz_zz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_zz_0_yy[i] = -2.0 * g_yz_zz_0_yy[i] * a_exp + 4.0 * g_xxyz_zz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_zz_0_yz[i] = -2.0 * g_yz_zz_0_yz[i] * a_exp + 4.0 * g_xxyz_zz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_xz_zz_0_zz[i] = -2.0 * g_yz_zz_0_zz[i] * a_exp + 4.0 * g_xxyz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_xy_0_0_0_yy_xx_0_xx, g_xy_0_0_0_yy_xx_0_xy, g_xy_0_0_0_yy_xx_0_xz, g_xy_0_0_0_yy_xx_0_yy, g_xy_0_0_0_yy_xx_0_yz, g_xy_0_0_0_yy_xx_0_zz, g_xy_xx_0_xx, g_xy_xx_0_xy, g_xy_xx_0_xz, g_xy_xx_0_yy, g_xy_xx_0_yz, g_xy_xx_0_zz, g_xyyy_xx_0_xx, g_xyyy_xx_0_xy, g_xyyy_xx_0_xz, g_xyyy_xx_0_yy, g_xyyy_xx_0_yz, g_xyyy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_xx_0_xx[i] = -4.0 * g_xy_xx_0_xx[i] * a_exp + 4.0 * g_xyyy_xx_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xx_0_xy[i] = -4.0 * g_xy_xx_0_xy[i] * a_exp + 4.0 * g_xyyy_xx_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xx_0_xz[i] = -4.0 * g_xy_xx_0_xz[i] * a_exp + 4.0 * g_xyyy_xx_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xx_0_yy[i] = -4.0 * g_xy_xx_0_yy[i] * a_exp + 4.0 * g_xyyy_xx_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xx_0_yz[i] = -4.0 * g_xy_xx_0_yz[i] * a_exp + 4.0 * g_xyyy_xx_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xx_0_zz[i] = -4.0 * g_xy_xx_0_zz[i] * a_exp + 4.0 * g_xyyy_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_xy_0_0_0_yy_xy_0_xx, g_xy_0_0_0_yy_xy_0_xy, g_xy_0_0_0_yy_xy_0_xz, g_xy_0_0_0_yy_xy_0_yy, g_xy_0_0_0_yy_xy_0_yz, g_xy_0_0_0_yy_xy_0_zz, g_xy_xy_0_xx, g_xy_xy_0_xy, g_xy_xy_0_xz, g_xy_xy_0_yy, g_xy_xy_0_yz, g_xy_xy_0_zz, g_xyyy_xy_0_xx, g_xyyy_xy_0_xy, g_xyyy_xy_0_xz, g_xyyy_xy_0_yy, g_xyyy_xy_0_yz, g_xyyy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_xy_0_xx[i] = -4.0 * g_xy_xy_0_xx[i] * a_exp + 4.0 * g_xyyy_xy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xy_0_xy[i] = -4.0 * g_xy_xy_0_xy[i] * a_exp + 4.0 * g_xyyy_xy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xy_0_xz[i] = -4.0 * g_xy_xy_0_xz[i] * a_exp + 4.0 * g_xyyy_xy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xy_0_yy[i] = -4.0 * g_xy_xy_0_yy[i] * a_exp + 4.0 * g_xyyy_xy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xy_0_yz[i] = -4.0 * g_xy_xy_0_yz[i] * a_exp + 4.0 * g_xyyy_xy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xy_0_zz[i] = -4.0 * g_xy_xy_0_zz[i] * a_exp + 4.0 * g_xyyy_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_xy_0_0_0_yy_xz_0_xx, g_xy_0_0_0_yy_xz_0_xy, g_xy_0_0_0_yy_xz_0_xz, g_xy_0_0_0_yy_xz_0_yy, g_xy_0_0_0_yy_xz_0_yz, g_xy_0_0_0_yy_xz_0_zz, g_xy_xz_0_xx, g_xy_xz_0_xy, g_xy_xz_0_xz, g_xy_xz_0_yy, g_xy_xz_0_yz, g_xy_xz_0_zz, g_xyyy_xz_0_xx, g_xyyy_xz_0_xy, g_xyyy_xz_0_xz, g_xyyy_xz_0_yy, g_xyyy_xz_0_yz, g_xyyy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_xz_0_xx[i] = -4.0 * g_xy_xz_0_xx[i] * a_exp + 4.0 * g_xyyy_xz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xz_0_xy[i] = -4.0 * g_xy_xz_0_xy[i] * a_exp + 4.0 * g_xyyy_xz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xz_0_xz[i] = -4.0 * g_xy_xz_0_xz[i] * a_exp + 4.0 * g_xyyy_xz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xz_0_yy[i] = -4.0 * g_xy_xz_0_yy[i] * a_exp + 4.0 * g_xyyy_xz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xz_0_yz[i] = -4.0 * g_xy_xz_0_yz[i] * a_exp + 4.0 * g_xyyy_xz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_xz_0_zz[i] = -4.0 * g_xy_xz_0_zz[i] * a_exp + 4.0 * g_xyyy_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_xy_0_0_0_yy_yy_0_xx, g_xy_0_0_0_yy_yy_0_xy, g_xy_0_0_0_yy_yy_0_xz, g_xy_0_0_0_yy_yy_0_yy, g_xy_0_0_0_yy_yy_0_yz, g_xy_0_0_0_yy_yy_0_zz, g_xy_yy_0_xx, g_xy_yy_0_xy, g_xy_yy_0_xz, g_xy_yy_0_yy, g_xy_yy_0_yz, g_xy_yy_0_zz, g_xyyy_yy_0_xx, g_xyyy_yy_0_xy, g_xyyy_yy_0_xz, g_xyyy_yy_0_yy, g_xyyy_yy_0_yz, g_xyyy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_yy_0_xx[i] = -4.0 * g_xy_yy_0_xx[i] * a_exp + 4.0 * g_xyyy_yy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yy_0_xy[i] = -4.0 * g_xy_yy_0_xy[i] * a_exp + 4.0 * g_xyyy_yy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yy_0_xz[i] = -4.0 * g_xy_yy_0_xz[i] * a_exp + 4.0 * g_xyyy_yy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yy_0_yy[i] = -4.0 * g_xy_yy_0_yy[i] * a_exp + 4.0 * g_xyyy_yy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yy_0_yz[i] = -4.0 * g_xy_yy_0_yz[i] * a_exp + 4.0 * g_xyyy_yy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yy_0_zz[i] = -4.0 * g_xy_yy_0_zz[i] * a_exp + 4.0 * g_xyyy_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_xy_0_0_0_yy_yz_0_xx, g_xy_0_0_0_yy_yz_0_xy, g_xy_0_0_0_yy_yz_0_xz, g_xy_0_0_0_yy_yz_0_yy, g_xy_0_0_0_yy_yz_0_yz, g_xy_0_0_0_yy_yz_0_zz, g_xy_yz_0_xx, g_xy_yz_0_xy, g_xy_yz_0_xz, g_xy_yz_0_yy, g_xy_yz_0_yz, g_xy_yz_0_zz, g_xyyy_yz_0_xx, g_xyyy_yz_0_xy, g_xyyy_yz_0_xz, g_xyyy_yz_0_yy, g_xyyy_yz_0_yz, g_xyyy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_yz_0_xx[i] = -4.0 * g_xy_yz_0_xx[i] * a_exp + 4.0 * g_xyyy_yz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yz_0_xy[i] = -4.0 * g_xy_yz_0_xy[i] * a_exp + 4.0 * g_xyyy_yz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yz_0_xz[i] = -4.0 * g_xy_yz_0_xz[i] * a_exp + 4.0 * g_xyyy_yz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yz_0_yy[i] = -4.0 * g_xy_yz_0_yy[i] * a_exp + 4.0 * g_xyyy_yz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yz_0_yz[i] = -4.0 * g_xy_yz_0_yz[i] * a_exp + 4.0 * g_xyyy_yz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_yz_0_zz[i] = -4.0 * g_xy_yz_0_zz[i] * a_exp + 4.0 * g_xyyy_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_xy_0_0_0_yy_zz_0_xx, g_xy_0_0_0_yy_zz_0_xy, g_xy_0_0_0_yy_zz_0_xz, g_xy_0_0_0_yy_zz_0_yy, g_xy_0_0_0_yy_zz_0_yz, g_xy_0_0_0_yy_zz_0_zz, g_xy_zz_0_xx, g_xy_zz_0_xy, g_xy_zz_0_xz, g_xy_zz_0_yy, g_xy_zz_0_yz, g_xy_zz_0_zz, g_xyyy_zz_0_xx, g_xyyy_zz_0_xy, g_xyyy_zz_0_xz, g_xyyy_zz_0_yy, g_xyyy_zz_0_yz, g_xyyy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yy_zz_0_xx[i] = -4.0 * g_xy_zz_0_xx[i] * a_exp + 4.0 * g_xyyy_zz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_zz_0_xy[i] = -4.0 * g_xy_zz_0_xy[i] * a_exp + 4.0 * g_xyyy_zz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_zz_0_xz[i] = -4.0 * g_xy_zz_0_xz[i] * a_exp + 4.0 * g_xyyy_zz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_zz_0_yy[i] = -4.0 * g_xy_zz_0_yy[i] * a_exp + 4.0 * g_xyyy_zz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_zz_0_yz[i] = -4.0 * g_xy_zz_0_yz[i] * a_exp + 4.0 * g_xyyy_zz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_yy_zz_0_zz[i] = -4.0 * g_xy_zz_0_zz[i] * a_exp + 4.0 * g_xyyy_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_xy_0_0_0_yz_xx_0_xx, g_xy_0_0_0_yz_xx_0_xy, g_xy_0_0_0_yz_xx_0_xz, g_xy_0_0_0_yz_xx_0_yy, g_xy_0_0_0_yz_xx_0_yz, g_xy_0_0_0_yz_xx_0_zz, g_xyyz_xx_0_xx, g_xyyz_xx_0_xy, g_xyyz_xx_0_xz, g_xyyz_xx_0_yy, g_xyyz_xx_0_yz, g_xyyz_xx_0_zz, g_xz_xx_0_xx, g_xz_xx_0_xy, g_xz_xx_0_xz, g_xz_xx_0_yy, g_xz_xx_0_yz, g_xz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_xx_0_xx[i] = -2.0 * g_xz_xx_0_xx[i] * a_exp + 4.0 * g_xyyz_xx_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xx_0_xy[i] = -2.0 * g_xz_xx_0_xy[i] * a_exp + 4.0 * g_xyyz_xx_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xx_0_xz[i] = -2.0 * g_xz_xx_0_xz[i] * a_exp + 4.0 * g_xyyz_xx_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xx_0_yy[i] = -2.0 * g_xz_xx_0_yy[i] * a_exp + 4.0 * g_xyyz_xx_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xx_0_yz[i] = -2.0 * g_xz_xx_0_yz[i] * a_exp + 4.0 * g_xyyz_xx_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xx_0_zz[i] = -2.0 * g_xz_xx_0_zz[i] * a_exp + 4.0 * g_xyyz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_xy_0_0_0_yz_xy_0_xx, g_xy_0_0_0_yz_xy_0_xy, g_xy_0_0_0_yz_xy_0_xz, g_xy_0_0_0_yz_xy_0_yy, g_xy_0_0_0_yz_xy_0_yz, g_xy_0_0_0_yz_xy_0_zz, g_xyyz_xy_0_xx, g_xyyz_xy_0_xy, g_xyyz_xy_0_xz, g_xyyz_xy_0_yy, g_xyyz_xy_0_yz, g_xyyz_xy_0_zz, g_xz_xy_0_xx, g_xz_xy_0_xy, g_xz_xy_0_xz, g_xz_xy_0_yy, g_xz_xy_0_yz, g_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_xy_0_xx[i] = -2.0 * g_xz_xy_0_xx[i] * a_exp + 4.0 * g_xyyz_xy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xy_0_xy[i] = -2.0 * g_xz_xy_0_xy[i] * a_exp + 4.0 * g_xyyz_xy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xy_0_xz[i] = -2.0 * g_xz_xy_0_xz[i] * a_exp + 4.0 * g_xyyz_xy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xy_0_yy[i] = -2.0 * g_xz_xy_0_yy[i] * a_exp + 4.0 * g_xyyz_xy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xy_0_yz[i] = -2.0 * g_xz_xy_0_yz[i] * a_exp + 4.0 * g_xyyz_xy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xy_0_zz[i] = -2.0 * g_xz_xy_0_zz[i] * a_exp + 4.0 * g_xyyz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_xy_0_0_0_yz_xz_0_xx, g_xy_0_0_0_yz_xz_0_xy, g_xy_0_0_0_yz_xz_0_xz, g_xy_0_0_0_yz_xz_0_yy, g_xy_0_0_0_yz_xz_0_yz, g_xy_0_0_0_yz_xz_0_zz, g_xyyz_xz_0_xx, g_xyyz_xz_0_xy, g_xyyz_xz_0_xz, g_xyyz_xz_0_yy, g_xyyz_xz_0_yz, g_xyyz_xz_0_zz, g_xz_xz_0_xx, g_xz_xz_0_xy, g_xz_xz_0_xz, g_xz_xz_0_yy, g_xz_xz_0_yz, g_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_xz_0_xx[i] = -2.0 * g_xz_xz_0_xx[i] * a_exp + 4.0 * g_xyyz_xz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xz_0_xy[i] = -2.0 * g_xz_xz_0_xy[i] * a_exp + 4.0 * g_xyyz_xz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xz_0_xz[i] = -2.0 * g_xz_xz_0_xz[i] * a_exp + 4.0 * g_xyyz_xz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xz_0_yy[i] = -2.0 * g_xz_xz_0_yy[i] * a_exp + 4.0 * g_xyyz_xz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xz_0_yz[i] = -2.0 * g_xz_xz_0_yz[i] * a_exp + 4.0 * g_xyyz_xz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_xz_0_zz[i] = -2.0 * g_xz_xz_0_zz[i] * a_exp + 4.0 * g_xyyz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_xy_0_0_0_yz_yy_0_xx, g_xy_0_0_0_yz_yy_0_xy, g_xy_0_0_0_yz_yy_0_xz, g_xy_0_0_0_yz_yy_0_yy, g_xy_0_0_0_yz_yy_0_yz, g_xy_0_0_0_yz_yy_0_zz, g_xyyz_yy_0_xx, g_xyyz_yy_0_xy, g_xyyz_yy_0_xz, g_xyyz_yy_0_yy, g_xyyz_yy_0_yz, g_xyyz_yy_0_zz, g_xz_yy_0_xx, g_xz_yy_0_xy, g_xz_yy_0_xz, g_xz_yy_0_yy, g_xz_yy_0_yz, g_xz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_yy_0_xx[i] = -2.0 * g_xz_yy_0_xx[i] * a_exp + 4.0 * g_xyyz_yy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yy_0_xy[i] = -2.0 * g_xz_yy_0_xy[i] * a_exp + 4.0 * g_xyyz_yy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yy_0_xz[i] = -2.0 * g_xz_yy_0_xz[i] * a_exp + 4.0 * g_xyyz_yy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yy_0_yy[i] = -2.0 * g_xz_yy_0_yy[i] * a_exp + 4.0 * g_xyyz_yy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yy_0_yz[i] = -2.0 * g_xz_yy_0_yz[i] * a_exp + 4.0 * g_xyyz_yy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yy_0_zz[i] = -2.0 * g_xz_yy_0_zz[i] * a_exp + 4.0 * g_xyyz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_xy_0_0_0_yz_yz_0_xx, g_xy_0_0_0_yz_yz_0_xy, g_xy_0_0_0_yz_yz_0_xz, g_xy_0_0_0_yz_yz_0_yy, g_xy_0_0_0_yz_yz_0_yz, g_xy_0_0_0_yz_yz_0_zz, g_xyyz_yz_0_xx, g_xyyz_yz_0_xy, g_xyyz_yz_0_xz, g_xyyz_yz_0_yy, g_xyyz_yz_0_yz, g_xyyz_yz_0_zz, g_xz_yz_0_xx, g_xz_yz_0_xy, g_xz_yz_0_xz, g_xz_yz_0_yy, g_xz_yz_0_yz, g_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_yz_0_xx[i] = -2.0 * g_xz_yz_0_xx[i] * a_exp + 4.0 * g_xyyz_yz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yz_0_xy[i] = -2.0 * g_xz_yz_0_xy[i] * a_exp + 4.0 * g_xyyz_yz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yz_0_xz[i] = -2.0 * g_xz_yz_0_xz[i] * a_exp + 4.0 * g_xyyz_yz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yz_0_yy[i] = -2.0 * g_xz_yz_0_yy[i] * a_exp + 4.0 * g_xyyz_yz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yz_0_yz[i] = -2.0 * g_xz_yz_0_yz[i] * a_exp + 4.0 * g_xyyz_yz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_yz_0_zz[i] = -2.0 * g_xz_yz_0_zz[i] * a_exp + 4.0 * g_xyyz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_xy_0_0_0_yz_zz_0_xx, g_xy_0_0_0_yz_zz_0_xy, g_xy_0_0_0_yz_zz_0_xz, g_xy_0_0_0_yz_zz_0_yy, g_xy_0_0_0_yz_zz_0_yz, g_xy_0_0_0_yz_zz_0_zz, g_xyyz_zz_0_xx, g_xyyz_zz_0_xy, g_xyyz_zz_0_xz, g_xyyz_zz_0_yy, g_xyyz_zz_0_yz, g_xyyz_zz_0_zz, g_xz_zz_0_xx, g_xz_zz_0_xy, g_xz_zz_0_xz, g_xz_zz_0_yy, g_xz_zz_0_yz, g_xz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_yz_zz_0_xx[i] = -2.0 * g_xz_zz_0_xx[i] * a_exp + 4.0 * g_xyyz_zz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_zz_0_xy[i] = -2.0 * g_xz_zz_0_xy[i] * a_exp + 4.0 * g_xyyz_zz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_zz_0_xz[i] = -2.0 * g_xz_zz_0_xz[i] * a_exp + 4.0 * g_xyyz_zz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_zz_0_yy[i] = -2.0 * g_xz_zz_0_yy[i] * a_exp + 4.0 * g_xyyz_zz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_zz_0_yz[i] = -2.0 * g_xz_zz_0_yz[i] * a_exp + 4.0 * g_xyyz_zz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_yz_zz_0_zz[i] = -2.0 * g_xz_zz_0_zz[i] * a_exp + 4.0 * g_xyyz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_xy_0_0_0_zz_xx_0_xx, g_xy_0_0_0_zz_xx_0_xy, g_xy_0_0_0_zz_xx_0_xz, g_xy_0_0_0_zz_xx_0_yy, g_xy_0_0_0_zz_xx_0_yz, g_xy_0_0_0_zz_xx_0_zz, g_xyzz_xx_0_xx, g_xyzz_xx_0_xy, g_xyzz_xx_0_xz, g_xyzz_xx_0_yy, g_xyzz_xx_0_yz, g_xyzz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_xx_0_xx[i] = 4.0 * g_xyzz_xx_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xx_0_xy[i] = 4.0 * g_xyzz_xx_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xx_0_xz[i] = 4.0 * g_xyzz_xx_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xx_0_yy[i] = 4.0 * g_xyzz_xx_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xx_0_yz[i] = 4.0 * g_xyzz_xx_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xx_0_zz[i] = 4.0 * g_xyzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_xy_0_0_0_zz_xy_0_xx, g_xy_0_0_0_zz_xy_0_xy, g_xy_0_0_0_zz_xy_0_xz, g_xy_0_0_0_zz_xy_0_yy, g_xy_0_0_0_zz_xy_0_yz, g_xy_0_0_0_zz_xy_0_zz, g_xyzz_xy_0_xx, g_xyzz_xy_0_xy, g_xyzz_xy_0_xz, g_xyzz_xy_0_yy, g_xyzz_xy_0_yz, g_xyzz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_xy_0_xx[i] = 4.0 * g_xyzz_xy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xy_0_xy[i] = 4.0 * g_xyzz_xy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xy_0_xz[i] = 4.0 * g_xyzz_xy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xy_0_yy[i] = 4.0 * g_xyzz_xy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xy_0_yz[i] = 4.0 * g_xyzz_xy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xy_0_zz[i] = 4.0 * g_xyzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_xy_0_0_0_zz_xz_0_xx, g_xy_0_0_0_zz_xz_0_xy, g_xy_0_0_0_zz_xz_0_xz, g_xy_0_0_0_zz_xz_0_yy, g_xy_0_0_0_zz_xz_0_yz, g_xy_0_0_0_zz_xz_0_zz, g_xyzz_xz_0_xx, g_xyzz_xz_0_xy, g_xyzz_xz_0_xz, g_xyzz_xz_0_yy, g_xyzz_xz_0_yz, g_xyzz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_xz_0_xx[i] = 4.0 * g_xyzz_xz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xz_0_xy[i] = 4.0 * g_xyzz_xz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xz_0_xz[i] = 4.0 * g_xyzz_xz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xz_0_yy[i] = 4.0 * g_xyzz_xz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xz_0_yz[i] = 4.0 * g_xyzz_xz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_xz_0_zz[i] = 4.0 * g_xyzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_xy_0_0_0_zz_yy_0_xx, g_xy_0_0_0_zz_yy_0_xy, g_xy_0_0_0_zz_yy_0_xz, g_xy_0_0_0_zz_yy_0_yy, g_xy_0_0_0_zz_yy_0_yz, g_xy_0_0_0_zz_yy_0_zz, g_xyzz_yy_0_xx, g_xyzz_yy_0_xy, g_xyzz_yy_0_xz, g_xyzz_yy_0_yy, g_xyzz_yy_0_yz, g_xyzz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_yy_0_xx[i] = 4.0 * g_xyzz_yy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yy_0_xy[i] = 4.0 * g_xyzz_yy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yy_0_xz[i] = 4.0 * g_xyzz_yy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yy_0_yy[i] = 4.0 * g_xyzz_yy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yy_0_yz[i] = 4.0 * g_xyzz_yy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yy_0_zz[i] = 4.0 * g_xyzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_xy_0_0_0_zz_yz_0_xx, g_xy_0_0_0_zz_yz_0_xy, g_xy_0_0_0_zz_yz_0_xz, g_xy_0_0_0_zz_yz_0_yy, g_xy_0_0_0_zz_yz_0_yz, g_xy_0_0_0_zz_yz_0_zz, g_xyzz_yz_0_xx, g_xyzz_yz_0_xy, g_xyzz_yz_0_xz, g_xyzz_yz_0_yy, g_xyzz_yz_0_yz, g_xyzz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_yz_0_xx[i] = 4.0 * g_xyzz_yz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yz_0_xy[i] = 4.0 * g_xyzz_yz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yz_0_xz[i] = 4.0 * g_xyzz_yz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yz_0_yy[i] = 4.0 * g_xyzz_yz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yz_0_yz[i] = 4.0 * g_xyzz_yz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_yz_0_zz[i] = 4.0 * g_xyzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_xy_0_0_0_zz_zz_0_xx, g_xy_0_0_0_zz_zz_0_xy, g_xy_0_0_0_zz_zz_0_xz, g_xy_0_0_0_zz_zz_0_yy, g_xy_0_0_0_zz_zz_0_yz, g_xy_0_0_0_zz_zz_0_zz, g_xyzz_zz_0_xx, g_xyzz_zz_0_xy, g_xyzz_zz_0_xz, g_xyzz_zz_0_yy, g_xyzz_zz_0_yz, g_xyzz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_zz_zz_0_xx[i] = 4.0 * g_xyzz_zz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_zz_0_xy[i] = 4.0 * g_xyzz_zz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_zz_0_xz[i] = 4.0 * g_xyzz_zz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_zz_0_yy[i] = 4.0 * g_xyzz_zz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_zz_0_yz[i] = 4.0 * g_xyzz_zz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_zz_zz_0_zz[i] = 4.0 * g_xyzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_xxxz_xx_0_xx, g_xxxz_xx_0_xy, g_xxxz_xx_0_xz, g_xxxz_xx_0_yy, g_xxxz_xx_0_yz, g_xxxz_xx_0_zz, g_xz_0_0_0_xx_xx_0_xx, g_xz_0_0_0_xx_xx_0_xy, g_xz_0_0_0_xx_xx_0_xz, g_xz_0_0_0_xx_xx_0_yy, g_xz_0_0_0_xx_xx_0_yz, g_xz_0_0_0_xx_xx_0_zz, g_xz_xx_0_xx, g_xz_xx_0_xy, g_xz_xx_0_xz, g_xz_xx_0_yy, g_xz_xx_0_yz, g_xz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_xx_0_xx[i] = -4.0 * g_xz_xx_0_xx[i] * a_exp + 4.0 * g_xxxz_xx_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xx_0_xy[i] = -4.0 * g_xz_xx_0_xy[i] * a_exp + 4.0 * g_xxxz_xx_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xx_0_xz[i] = -4.0 * g_xz_xx_0_xz[i] * a_exp + 4.0 * g_xxxz_xx_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xx_0_yy[i] = -4.0 * g_xz_xx_0_yy[i] * a_exp + 4.0 * g_xxxz_xx_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xx_0_yz[i] = -4.0 * g_xz_xx_0_yz[i] * a_exp + 4.0 * g_xxxz_xx_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xx_0_zz[i] = -4.0 * g_xz_xx_0_zz[i] * a_exp + 4.0 * g_xxxz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_xxxz_xy_0_xx, g_xxxz_xy_0_xy, g_xxxz_xy_0_xz, g_xxxz_xy_0_yy, g_xxxz_xy_0_yz, g_xxxz_xy_0_zz, g_xz_0_0_0_xx_xy_0_xx, g_xz_0_0_0_xx_xy_0_xy, g_xz_0_0_0_xx_xy_0_xz, g_xz_0_0_0_xx_xy_0_yy, g_xz_0_0_0_xx_xy_0_yz, g_xz_0_0_0_xx_xy_0_zz, g_xz_xy_0_xx, g_xz_xy_0_xy, g_xz_xy_0_xz, g_xz_xy_0_yy, g_xz_xy_0_yz, g_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_xy_0_xx[i] = -4.0 * g_xz_xy_0_xx[i] * a_exp + 4.0 * g_xxxz_xy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xy_0_xy[i] = -4.0 * g_xz_xy_0_xy[i] * a_exp + 4.0 * g_xxxz_xy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xy_0_xz[i] = -4.0 * g_xz_xy_0_xz[i] * a_exp + 4.0 * g_xxxz_xy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xy_0_yy[i] = -4.0 * g_xz_xy_0_yy[i] * a_exp + 4.0 * g_xxxz_xy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xy_0_yz[i] = -4.0 * g_xz_xy_0_yz[i] * a_exp + 4.0 * g_xxxz_xy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xy_0_zz[i] = -4.0 * g_xz_xy_0_zz[i] * a_exp + 4.0 * g_xxxz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_xxxz_xz_0_xx, g_xxxz_xz_0_xy, g_xxxz_xz_0_xz, g_xxxz_xz_0_yy, g_xxxz_xz_0_yz, g_xxxz_xz_0_zz, g_xz_0_0_0_xx_xz_0_xx, g_xz_0_0_0_xx_xz_0_xy, g_xz_0_0_0_xx_xz_0_xz, g_xz_0_0_0_xx_xz_0_yy, g_xz_0_0_0_xx_xz_0_yz, g_xz_0_0_0_xx_xz_0_zz, g_xz_xz_0_xx, g_xz_xz_0_xy, g_xz_xz_0_xz, g_xz_xz_0_yy, g_xz_xz_0_yz, g_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_xz_0_xx[i] = -4.0 * g_xz_xz_0_xx[i] * a_exp + 4.0 * g_xxxz_xz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xz_0_xy[i] = -4.0 * g_xz_xz_0_xy[i] * a_exp + 4.0 * g_xxxz_xz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xz_0_xz[i] = -4.0 * g_xz_xz_0_xz[i] * a_exp + 4.0 * g_xxxz_xz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xz_0_yy[i] = -4.0 * g_xz_xz_0_yy[i] * a_exp + 4.0 * g_xxxz_xz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xz_0_yz[i] = -4.0 * g_xz_xz_0_yz[i] * a_exp + 4.0 * g_xxxz_xz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_xz_0_zz[i] = -4.0 * g_xz_xz_0_zz[i] * a_exp + 4.0 * g_xxxz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_xxxz_yy_0_xx, g_xxxz_yy_0_xy, g_xxxz_yy_0_xz, g_xxxz_yy_0_yy, g_xxxz_yy_0_yz, g_xxxz_yy_0_zz, g_xz_0_0_0_xx_yy_0_xx, g_xz_0_0_0_xx_yy_0_xy, g_xz_0_0_0_xx_yy_0_xz, g_xz_0_0_0_xx_yy_0_yy, g_xz_0_0_0_xx_yy_0_yz, g_xz_0_0_0_xx_yy_0_zz, g_xz_yy_0_xx, g_xz_yy_0_xy, g_xz_yy_0_xz, g_xz_yy_0_yy, g_xz_yy_0_yz, g_xz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_yy_0_xx[i] = -4.0 * g_xz_yy_0_xx[i] * a_exp + 4.0 * g_xxxz_yy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yy_0_xy[i] = -4.0 * g_xz_yy_0_xy[i] * a_exp + 4.0 * g_xxxz_yy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yy_0_xz[i] = -4.0 * g_xz_yy_0_xz[i] * a_exp + 4.0 * g_xxxz_yy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yy_0_yy[i] = -4.0 * g_xz_yy_0_yy[i] * a_exp + 4.0 * g_xxxz_yy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yy_0_yz[i] = -4.0 * g_xz_yy_0_yz[i] * a_exp + 4.0 * g_xxxz_yy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yy_0_zz[i] = -4.0 * g_xz_yy_0_zz[i] * a_exp + 4.0 * g_xxxz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_xxxz_yz_0_xx, g_xxxz_yz_0_xy, g_xxxz_yz_0_xz, g_xxxz_yz_0_yy, g_xxxz_yz_0_yz, g_xxxz_yz_0_zz, g_xz_0_0_0_xx_yz_0_xx, g_xz_0_0_0_xx_yz_0_xy, g_xz_0_0_0_xx_yz_0_xz, g_xz_0_0_0_xx_yz_0_yy, g_xz_0_0_0_xx_yz_0_yz, g_xz_0_0_0_xx_yz_0_zz, g_xz_yz_0_xx, g_xz_yz_0_xy, g_xz_yz_0_xz, g_xz_yz_0_yy, g_xz_yz_0_yz, g_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_yz_0_xx[i] = -4.0 * g_xz_yz_0_xx[i] * a_exp + 4.0 * g_xxxz_yz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yz_0_xy[i] = -4.0 * g_xz_yz_0_xy[i] * a_exp + 4.0 * g_xxxz_yz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yz_0_xz[i] = -4.0 * g_xz_yz_0_xz[i] * a_exp + 4.0 * g_xxxz_yz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yz_0_yy[i] = -4.0 * g_xz_yz_0_yy[i] * a_exp + 4.0 * g_xxxz_yz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yz_0_yz[i] = -4.0 * g_xz_yz_0_yz[i] * a_exp + 4.0 * g_xxxz_yz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_yz_0_zz[i] = -4.0 * g_xz_yz_0_zz[i] * a_exp + 4.0 * g_xxxz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_xxxz_zz_0_xx, g_xxxz_zz_0_xy, g_xxxz_zz_0_xz, g_xxxz_zz_0_yy, g_xxxz_zz_0_yz, g_xxxz_zz_0_zz, g_xz_0_0_0_xx_zz_0_xx, g_xz_0_0_0_xx_zz_0_xy, g_xz_0_0_0_xx_zz_0_xz, g_xz_0_0_0_xx_zz_0_yy, g_xz_0_0_0_xx_zz_0_yz, g_xz_0_0_0_xx_zz_0_zz, g_xz_zz_0_xx, g_xz_zz_0_xy, g_xz_zz_0_xz, g_xz_zz_0_yy, g_xz_zz_0_yz, g_xz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xx_zz_0_xx[i] = -4.0 * g_xz_zz_0_xx[i] * a_exp + 4.0 * g_xxxz_zz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_zz_0_xy[i] = -4.0 * g_xz_zz_0_xy[i] * a_exp + 4.0 * g_xxxz_zz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_zz_0_xz[i] = -4.0 * g_xz_zz_0_xz[i] * a_exp + 4.0 * g_xxxz_zz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_zz_0_yy[i] = -4.0 * g_xz_zz_0_yy[i] * a_exp + 4.0 * g_xxxz_zz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_zz_0_yz[i] = -4.0 * g_xz_zz_0_yz[i] * a_exp + 4.0 * g_xxxz_zz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xx_zz_0_zz[i] = -4.0 * g_xz_zz_0_zz[i] * a_exp + 4.0 * g_xxxz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_xxyz_xx_0_xx, g_xxyz_xx_0_xy, g_xxyz_xx_0_xz, g_xxyz_xx_0_yy, g_xxyz_xx_0_yz, g_xxyz_xx_0_zz, g_xz_0_0_0_xy_xx_0_xx, g_xz_0_0_0_xy_xx_0_xy, g_xz_0_0_0_xy_xx_0_xz, g_xz_0_0_0_xy_xx_0_yy, g_xz_0_0_0_xy_xx_0_yz, g_xz_0_0_0_xy_xx_0_zz, g_yz_xx_0_xx, g_yz_xx_0_xy, g_yz_xx_0_xz, g_yz_xx_0_yy, g_yz_xx_0_yz, g_yz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_xx_0_xx[i] = -2.0 * g_yz_xx_0_xx[i] * a_exp + 4.0 * g_xxyz_xx_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xx_0_xy[i] = -2.0 * g_yz_xx_0_xy[i] * a_exp + 4.0 * g_xxyz_xx_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xx_0_xz[i] = -2.0 * g_yz_xx_0_xz[i] * a_exp + 4.0 * g_xxyz_xx_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xx_0_yy[i] = -2.0 * g_yz_xx_0_yy[i] * a_exp + 4.0 * g_xxyz_xx_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xx_0_yz[i] = -2.0 * g_yz_xx_0_yz[i] * a_exp + 4.0 * g_xxyz_xx_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xx_0_zz[i] = -2.0 * g_yz_xx_0_zz[i] * a_exp + 4.0 * g_xxyz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_xxyz_xy_0_xx, g_xxyz_xy_0_xy, g_xxyz_xy_0_xz, g_xxyz_xy_0_yy, g_xxyz_xy_0_yz, g_xxyz_xy_0_zz, g_xz_0_0_0_xy_xy_0_xx, g_xz_0_0_0_xy_xy_0_xy, g_xz_0_0_0_xy_xy_0_xz, g_xz_0_0_0_xy_xy_0_yy, g_xz_0_0_0_xy_xy_0_yz, g_xz_0_0_0_xy_xy_0_zz, g_yz_xy_0_xx, g_yz_xy_0_xy, g_yz_xy_0_xz, g_yz_xy_0_yy, g_yz_xy_0_yz, g_yz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_xy_0_xx[i] = -2.0 * g_yz_xy_0_xx[i] * a_exp + 4.0 * g_xxyz_xy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xy_0_xy[i] = -2.0 * g_yz_xy_0_xy[i] * a_exp + 4.0 * g_xxyz_xy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xy_0_xz[i] = -2.0 * g_yz_xy_0_xz[i] * a_exp + 4.0 * g_xxyz_xy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xy_0_yy[i] = -2.0 * g_yz_xy_0_yy[i] * a_exp + 4.0 * g_xxyz_xy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xy_0_yz[i] = -2.0 * g_yz_xy_0_yz[i] * a_exp + 4.0 * g_xxyz_xy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xy_0_zz[i] = -2.0 * g_yz_xy_0_zz[i] * a_exp + 4.0 * g_xxyz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_xxyz_xz_0_xx, g_xxyz_xz_0_xy, g_xxyz_xz_0_xz, g_xxyz_xz_0_yy, g_xxyz_xz_0_yz, g_xxyz_xz_0_zz, g_xz_0_0_0_xy_xz_0_xx, g_xz_0_0_0_xy_xz_0_xy, g_xz_0_0_0_xy_xz_0_xz, g_xz_0_0_0_xy_xz_0_yy, g_xz_0_0_0_xy_xz_0_yz, g_xz_0_0_0_xy_xz_0_zz, g_yz_xz_0_xx, g_yz_xz_0_xy, g_yz_xz_0_xz, g_yz_xz_0_yy, g_yz_xz_0_yz, g_yz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_xz_0_xx[i] = -2.0 * g_yz_xz_0_xx[i] * a_exp + 4.0 * g_xxyz_xz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xz_0_xy[i] = -2.0 * g_yz_xz_0_xy[i] * a_exp + 4.0 * g_xxyz_xz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xz_0_xz[i] = -2.0 * g_yz_xz_0_xz[i] * a_exp + 4.0 * g_xxyz_xz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xz_0_yy[i] = -2.0 * g_yz_xz_0_yy[i] * a_exp + 4.0 * g_xxyz_xz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xz_0_yz[i] = -2.0 * g_yz_xz_0_yz[i] * a_exp + 4.0 * g_xxyz_xz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_xz_0_zz[i] = -2.0 * g_yz_xz_0_zz[i] * a_exp + 4.0 * g_xxyz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (486-492)

    #pragma omp simd aligned(g_xxyz_yy_0_xx, g_xxyz_yy_0_xy, g_xxyz_yy_0_xz, g_xxyz_yy_0_yy, g_xxyz_yy_0_yz, g_xxyz_yy_0_zz, g_xz_0_0_0_xy_yy_0_xx, g_xz_0_0_0_xy_yy_0_xy, g_xz_0_0_0_xy_yy_0_xz, g_xz_0_0_0_xy_yy_0_yy, g_xz_0_0_0_xy_yy_0_yz, g_xz_0_0_0_xy_yy_0_zz, g_yz_yy_0_xx, g_yz_yy_0_xy, g_yz_yy_0_xz, g_yz_yy_0_yy, g_yz_yy_0_yz, g_yz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_yy_0_xx[i] = -2.0 * g_yz_yy_0_xx[i] * a_exp + 4.0 * g_xxyz_yy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yy_0_xy[i] = -2.0 * g_yz_yy_0_xy[i] * a_exp + 4.0 * g_xxyz_yy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yy_0_xz[i] = -2.0 * g_yz_yy_0_xz[i] * a_exp + 4.0 * g_xxyz_yy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yy_0_yy[i] = -2.0 * g_yz_yy_0_yy[i] * a_exp + 4.0 * g_xxyz_yy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yy_0_yz[i] = -2.0 * g_yz_yy_0_yz[i] * a_exp + 4.0 * g_xxyz_yy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yy_0_zz[i] = -2.0 * g_yz_yy_0_zz[i] * a_exp + 4.0 * g_xxyz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (492-498)

    #pragma omp simd aligned(g_xxyz_yz_0_xx, g_xxyz_yz_0_xy, g_xxyz_yz_0_xz, g_xxyz_yz_0_yy, g_xxyz_yz_0_yz, g_xxyz_yz_0_zz, g_xz_0_0_0_xy_yz_0_xx, g_xz_0_0_0_xy_yz_0_xy, g_xz_0_0_0_xy_yz_0_xz, g_xz_0_0_0_xy_yz_0_yy, g_xz_0_0_0_xy_yz_0_yz, g_xz_0_0_0_xy_yz_0_zz, g_yz_yz_0_xx, g_yz_yz_0_xy, g_yz_yz_0_xz, g_yz_yz_0_yy, g_yz_yz_0_yz, g_yz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_yz_0_xx[i] = -2.0 * g_yz_yz_0_xx[i] * a_exp + 4.0 * g_xxyz_yz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yz_0_xy[i] = -2.0 * g_yz_yz_0_xy[i] * a_exp + 4.0 * g_xxyz_yz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yz_0_xz[i] = -2.0 * g_yz_yz_0_xz[i] * a_exp + 4.0 * g_xxyz_yz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yz_0_yy[i] = -2.0 * g_yz_yz_0_yy[i] * a_exp + 4.0 * g_xxyz_yz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yz_0_yz[i] = -2.0 * g_yz_yz_0_yz[i] * a_exp + 4.0 * g_xxyz_yz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_yz_0_zz[i] = -2.0 * g_yz_yz_0_zz[i] * a_exp + 4.0 * g_xxyz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (498-504)

    #pragma omp simd aligned(g_xxyz_zz_0_xx, g_xxyz_zz_0_xy, g_xxyz_zz_0_xz, g_xxyz_zz_0_yy, g_xxyz_zz_0_yz, g_xxyz_zz_0_zz, g_xz_0_0_0_xy_zz_0_xx, g_xz_0_0_0_xy_zz_0_xy, g_xz_0_0_0_xy_zz_0_xz, g_xz_0_0_0_xy_zz_0_yy, g_xz_0_0_0_xy_zz_0_yz, g_xz_0_0_0_xy_zz_0_zz, g_yz_zz_0_xx, g_yz_zz_0_xy, g_yz_zz_0_xz, g_yz_zz_0_yy, g_yz_zz_0_yz, g_yz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xy_zz_0_xx[i] = -2.0 * g_yz_zz_0_xx[i] * a_exp + 4.0 * g_xxyz_zz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_zz_0_xy[i] = -2.0 * g_yz_zz_0_xy[i] * a_exp + 4.0 * g_xxyz_zz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_zz_0_xz[i] = -2.0 * g_yz_zz_0_xz[i] * a_exp + 4.0 * g_xxyz_zz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_zz_0_yy[i] = -2.0 * g_yz_zz_0_yy[i] * a_exp + 4.0 * g_xxyz_zz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_zz_0_yz[i] = -2.0 * g_yz_zz_0_yz[i] * a_exp + 4.0 * g_xxyz_zz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xy_zz_0_zz[i] = -2.0 * g_yz_zz_0_zz[i] * a_exp + 4.0 * g_xxyz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (504-510)

    #pragma omp simd aligned(g_0_xx_0_xx, g_0_xx_0_xy, g_0_xx_0_xz, g_0_xx_0_yy, g_0_xx_0_yz, g_0_xx_0_zz, g_xx_xx_0_xx, g_xx_xx_0_xy, g_xx_xx_0_xz, g_xx_xx_0_yy, g_xx_xx_0_yz, g_xx_xx_0_zz, g_xxzz_xx_0_xx, g_xxzz_xx_0_xy, g_xxzz_xx_0_xz, g_xxzz_xx_0_yy, g_xxzz_xx_0_yz, g_xxzz_xx_0_zz, g_xz_0_0_0_xz_xx_0_xx, g_xz_0_0_0_xz_xx_0_xy, g_xz_0_0_0_xz_xx_0_xz, g_xz_0_0_0_xz_xx_0_yy, g_xz_0_0_0_xz_xx_0_yz, g_xz_0_0_0_xz_xx_0_zz, g_zz_xx_0_xx, g_zz_xx_0_xy, g_zz_xx_0_xz, g_zz_xx_0_yy, g_zz_xx_0_yz, g_zz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_xx_0_xx[i] = g_0_xx_0_xx[i] - 2.0 * g_zz_xx_0_xx[i] * a_exp - 2.0 * g_xx_xx_0_xx[i] * a_exp + 4.0 * g_xxzz_xx_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xx_0_xy[i] = g_0_xx_0_xy[i] - 2.0 * g_zz_xx_0_xy[i] * a_exp - 2.0 * g_xx_xx_0_xy[i] * a_exp + 4.0 * g_xxzz_xx_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xx_0_xz[i] = g_0_xx_0_xz[i] - 2.0 * g_zz_xx_0_xz[i] * a_exp - 2.0 * g_xx_xx_0_xz[i] * a_exp + 4.0 * g_xxzz_xx_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xx_0_yy[i] = g_0_xx_0_yy[i] - 2.0 * g_zz_xx_0_yy[i] * a_exp - 2.0 * g_xx_xx_0_yy[i] * a_exp + 4.0 * g_xxzz_xx_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xx_0_yz[i] = g_0_xx_0_yz[i] - 2.0 * g_zz_xx_0_yz[i] * a_exp - 2.0 * g_xx_xx_0_yz[i] * a_exp + 4.0 * g_xxzz_xx_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xx_0_zz[i] = g_0_xx_0_zz[i] - 2.0 * g_zz_xx_0_zz[i] * a_exp - 2.0 * g_xx_xx_0_zz[i] * a_exp + 4.0 * g_xxzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (510-516)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_xx_xy_0_xx, g_xx_xy_0_xy, g_xx_xy_0_xz, g_xx_xy_0_yy, g_xx_xy_0_yz, g_xx_xy_0_zz, g_xxzz_xy_0_xx, g_xxzz_xy_0_xy, g_xxzz_xy_0_xz, g_xxzz_xy_0_yy, g_xxzz_xy_0_yz, g_xxzz_xy_0_zz, g_xz_0_0_0_xz_xy_0_xx, g_xz_0_0_0_xz_xy_0_xy, g_xz_0_0_0_xz_xy_0_xz, g_xz_0_0_0_xz_xy_0_yy, g_xz_0_0_0_xz_xy_0_yz, g_xz_0_0_0_xz_xy_0_zz, g_zz_xy_0_xx, g_zz_xy_0_xy, g_zz_xy_0_xz, g_zz_xy_0_yy, g_zz_xy_0_yz, g_zz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_xy_0_xx[i] = g_0_xy_0_xx[i] - 2.0 * g_zz_xy_0_xx[i] * a_exp - 2.0 * g_xx_xy_0_xx[i] * a_exp + 4.0 * g_xxzz_xy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xy_0_xy[i] = g_0_xy_0_xy[i] - 2.0 * g_zz_xy_0_xy[i] * a_exp - 2.0 * g_xx_xy_0_xy[i] * a_exp + 4.0 * g_xxzz_xy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xy_0_xz[i] = g_0_xy_0_xz[i] - 2.0 * g_zz_xy_0_xz[i] * a_exp - 2.0 * g_xx_xy_0_xz[i] * a_exp + 4.0 * g_xxzz_xy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xy_0_yy[i] = g_0_xy_0_yy[i] - 2.0 * g_zz_xy_0_yy[i] * a_exp - 2.0 * g_xx_xy_0_yy[i] * a_exp + 4.0 * g_xxzz_xy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xy_0_yz[i] = g_0_xy_0_yz[i] - 2.0 * g_zz_xy_0_yz[i] * a_exp - 2.0 * g_xx_xy_0_yz[i] * a_exp + 4.0 * g_xxzz_xy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xy_0_zz[i] = g_0_xy_0_zz[i] - 2.0 * g_zz_xy_0_zz[i] * a_exp - 2.0 * g_xx_xy_0_zz[i] * a_exp + 4.0 * g_xxzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (516-522)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_xx_xz_0_xx, g_xx_xz_0_xy, g_xx_xz_0_xz, g_xx_xz_0_yy, g_xx_xz_0_yz, g_xx_xz_0_zz, g_xxzz_xz_0_xx, g_xxzz_xz_0_xy, g_xxzz_xz_0_xz, g_xxzz_xz_0_yy, g_xxzz_xz_0_yz, g_xxzz_xz_0_zz, g_xz_0_0_0_xz_xz_0_xx, g_xz_0_0_0_xz_xz_0_xy, g_xz_0_0_0_xz_xz_0_xz, g_xz_0_0_0_xz_xz_0_yy, g_xz_0_0_0_xz_xz_0_yz, g_xz_0_0_0_xz_xz_0_zz, g_zz_xz_0_xx, g_zz_xz_0_xy, g_zz_xz_0_xz, g_zz_xz_0_yy, g_zz_xz_0_yz, g_zz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_xz_0_xx[i] = g_0_xz_0_xx[i] - 2.0 * g_zz_xz_0_xx[i] * a_exp - 2.0 * g_xx_xz_0_xx[i] * a_exp + 4.0 * g_xxzz_xz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xz_0_xy[i] = g_0_xz_0_xy[i] - 2.0 * g_zz_xz_0_xy[i] * a_exp - 2.0 * g_xx_xz_0_xy[i] * a_exp + 4.0 * g_xxzz_xz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xz_0_xz[i] = g_0_xz_0_xz[i] - 2.0 * g_zz_xz_0_xz[i] * a_exp - 2.0 * g_xx_xz_0_xz[i] * a_exp + 4.0 * g_xxzz_xz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xz_0_yy[i] = g_0_xz_0_yy[i] - 2.0 * g_zz_xz_0_yy[i] * a_exp - 2.0 * g_xx_xz_0_yy[i] * a_exp + 4.0 * g_xxzz_xz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xz_0_yz[i] = g_0_xz_0_yz[i] - 2.0 * g_zz_xz_0_yz[i] * a_exp - 2.0 * g_xx_xz_0_yz[i] * a_exp + 4.0 * g_xxzz_xz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_xz_0_zz[i] = g_0_xz_0_zz[i] - 2.0 * g_zz_xz_0_zz[i] * a_exp - 2.0 * g_xx_xz_0_zz[i] * a_exp + 4.0 * g_xxzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (522-528)

    #pragma omp simd aligned(g_0_yy_0_xx, g_0_yy_0_xy, g_0_yy_0_xz, g_0_yy_0_yy, g_0_yy_0_yz, g_0_yy_0_zz, g_xx_yy_0_xx, g_xx_yy_0_xy, g_xx_yy_0_xz, g_xx_yy_0_yy, g_xx_yy_0_yz, g_xx_yy_0_zz, g_xxzz_yy_0_xx, g_xxzz_yy_0_xy, g_xxzz_yy_0_xz, g_xxzz_yy_0_yy, g_xxzz_yy_0_yz, g_xxzz_yy_0_zz, g_xz_0_0_0_xz_yy_0_xx, g_xz_0_0_0_xz_yy_0_xy, g_xz_0_0_0_xz_yy_0_xz, g_xz_0_0_0_xz_yy_0_yy, g_xz_0_0_0_xz_yy_0_yz, g_xz_0_0_0_xz_yy_0_zz, g_zz_yy_0_xx, g_zz_yy_0_xy, g_zz_yy_0_xz, g_zz_yy_0_yy, g_zz_yy_0_yz, g_zz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_yy_0_xx[i] = g_0_yy_0_xx[i] - 2.0 * g_zz_yy_0_xx[i] * a_exp - 2.0 * g_xx_yy_0_xx[i] * a_exp + 4.0 * g_xxzz_yy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yy_0_xy[i] = g_0_yy_0_xy[i] - 2.0 * g_zz_yy_0_xy[i] * a_exp - 2.0 * g_xx_yy_0_xy[i] * a_exp + 4.0 * g_xxzz_yy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yy_0_xz[i] = g_0_yy_0_xz[i] - 2.0 * g_zz_yy_0_xz[i] * a_exp - 2.0 * g_xx_yy_0_xz[i] * a_exp + 4.0 * g_xxzz_yy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yy_0_yy[i] = g_0_yy_0_yy[i] - 2.0 * g_zz_yy_0_yy[i] * a_exp - 2.0 * g_xx_yy_0_yy[i] * a_exp + 4.0 * g_xxzz_yy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yy_0_yz[i] = g_0_yy_0_yz[i] - 2.0 * g_zz_yy_0_yz[i] * a_exp - 2.0 * g_xx_yy_0_yz[i] * a_exp + 4.0 * g_xxzz_yy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yy_0_zz[i] = g_0_yy_0_zz[i] - 2.0 * g_zz_yy_0_zz[i] * a_exp - 2.0 * g_xx_yy_0_zz[i] * a_exp + 4.0 * g_xxzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (528-534)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_xx_yz_0_xx, g_xx_yz_0_xy, g_xx_yz_0_xz, g_xx_yz_0_yy, g_xx_yz_0_yz, g_xx_yz_0_zz, g_xxzz_yz_0_xx, g_xxzz_yz_0_xy, g_xxzz_yz_0_xz, g_xxzz_yz_0_yy, g_xxzz_yz_0_yz, g_xxzz_yz_0_zz, g_xz_0_0_0_xz_yz_0_xx, g_xz_0_0_0_xz_yz_0_xy, g_xz_0_0_0_xz_yz_0_xz, g_xz_0_0_0_xz_yz_0_yy, g_xz_0_0_0_xz_yz_0_yz, g_xz_0_0_0_xz_yz_0_zz, g_zz_yz_0_xx, g_zz_yz_0_xy, g_zz_yz_0_xz, g_zz_yz_0_yy, g_zz_yz_0_yz, g_zz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_yz_0_xx[i] = g_0_yz_0_xx[i] - 2.0 * g_zz_yz_0_xx[i] * a_exp - 2.0 * g_xx_yz_0_xx[i] * a_exp + 4.0 * g_xxzz_yz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yz_0_xy[i] = g_0_yz_0_xy[i] - 2.0 * g_zz_yz_0_xy[i] * a_exp - 2.0 * g_xx_yz_0_xy[i] * a_exp + 4.0 * g_xxzz_yz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yz_0_xz[i] = g_0_yz_0_xz[i] - 2.0 * g_zz_yz_0_xz[i] * a_exp - 2.0 * g_xx_yz_0_xz[i] * a_exp + 4.0 * g_xxzz_yz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yz_0_yy[i] = g_0_yz_0_yy[i] - 2.0 * g_zz_yz_0_yy[i] * a_exp - 2.0 * g_xx_yz_0_yy[i] * a_exp + 4.0 * g_xxzz_yz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yz_0_yz[i] = g_0_yz_0_yz[i] - 2.0 * g_zz_yz_0_yz[i] * a_exp - 2.0 * g_xx_yz_0_yz[i] * a_exp + 4.0 * g_xxzz_yz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_yz_0_zz[i] = g_0_yz_0_zz[i] - 2.0 * g_zz_yz_0_zz[i] * a_exp - 2.0 * g_xx_yz_0_zz[i] * a_exp + 4.0 * g_xxzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (534-540)

    #pragma omp simd aligned(g_0_zz_0_xx, g_0_zz_0_xy, g_0_zz_0_xz, g_0_zz_0_yy, g_0_zz_0_yz, g_0_zz_0_zz, g_xx_zz_0_xx, g_xx_zz_0_xy, g_xx_zz_0_xz, g_xx_zz_0_yy, g_xx_zz_0_yz, g_xx_zz_0_zz, g_xxzz_zz_0_xx, g_xxzz_zz_0_xy, g_xxzz_zz_0_xz, g_xxzz_zz_0_yy, g_xxzz_zz_0_yz, g_xxzz_zz_0_zz, g_xz_0_0_0_xz_zz_0_xx, g_xz_0_0_0_xz_zz_0_xy, g_xz_0_0_0_xz_zz_0_xz, g_xz_0_0_0_xz_zz_0_yy, g_xz_0_0_0_xz_zz_0_yz, g_xz_0_0_0_xz_zz_0_zz, g_zz_zz_0_xx, g_zz_zz_0_xy, g_zz_zz_0_xz, g_zz_zz_0_yy, g_zz_zz_0_yz, g_zz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_xz_zz_0_xx[i] = g_0_zz_0_xx[i] - 2.0 * g_zz_zz_0_xx[i] * a_exp - 2.0 * g_xx_zz_0_xx[i] * a_exp + 4.0 * g_xxzz_zz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_zz_0_xy[i] = g_0_zz_0_xy[i] - 2.0 * g_zz_zz_0_xy[i] * a_exp - 2.0 * g_xx_zz_0_xy[i] * a_exp + 4.0 * g_xxzz_zz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_zz_0_xz[i] = g_0_zz_0_xz[i] - 2.0 * g_zz_zz_0_xz[i] * a_exp - 2.0 * g_xx_zz_0_xz[i] * a_exp + 4.0 * g_xxzz_zz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_zz_0_yy[i] = g_0_zz_0_yy[i] - 2.0 * g_zz_zz_0_yy[i] * a_exp - 2.0 * g_xx_zz_0_yy[i] * a_exp + 4.0 * g_xxzz_zz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_zz_0_yz[i] = g_0_zz_0_yz[i] - 2.0 * g_zz_zz_0_yz[i] * a_exp - 2.0 * g_xx_zz_0_yz[i] * a_exp + 4.0 * g_xxzz_zz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_xz_zz_0_zz[i] = g_0_zz_0_zz[i] - 2.0 * g_zz_zz_0_zz[i] * a_exp - 2.0 * g_xx_zz_0_zz[i] * a_exp + 4.0 * g_xxzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (540-546)

    #pragma omp simd aligned(g_xyyz_xx_0_xx, g_xyyz_xx_0_xy, g_xyyz_xx_0_xz, g_xyyz_xx_0_yy, g_xyyz_xx_0_yz, g_xyyz_xx_0_zz, g_xz_0_0_0_yy_xx_0_xx, g_xz_0_0_0_yy_xx_0_xy, g_xz_0_0_0_yy_xx_0_xz, g_xz_0_0_0_yy_xx_0_yy, g_xz_0_0_0_yy_xx_0_yz, g_xz_0_0_0_yy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_xx_0_xx[i] = 4.0 * g_xyyz_xx_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xx_0_xy[i] = 4.0 * g_xyyz_xx_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xx_0_xz[i] = 4.0 * g_xyyz_xx_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xx_0_yy[i] = 4.0 * g_xyyz_xx_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xx_0_yz[i] = 4.0 * g_xyyz_xx_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xx_0_zz[i] = 4.0 * g_xyyz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (546-552)

    #pragma omp simd aligned(g_xyyz_xy_0_xx, g_xyyz_xy_0_xy, g_xyyz_xy_0_xz, g_xyyz_xy_0_yy, g_xyyz_xy_0_yz, g_xyyz_xy_0_zz, g_xz_0_0_0_yy_xy_0_xx, g_xz_0_0_0_yy_xy_0_xy, g_xz_0_0_0_yy_xy_0_xz, g_xz_0_0_0_yy_xy_0_yy, g_xz_0_0_0_yy_xy_0_yz, g_xz_0_0_0_yy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_xy_0_xx[i] = 4.0 * g_xyyz_xy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xy_0_xy[i] = 4.0 * g_xyyz_xy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xy_0_xz[i] = 4.0 * g_xyyz_xy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xy_0_yy[i] = 4.0 * g_xyyz_xy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xy_0_yz[i] = 4.0 * g_xyyz_xy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xy_0_zz[i] = 4.0 * g_xyyz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (552-558)

    #pragma omp simd aligned(g_xyyz_xz_0_xx, g_xyyz_xz_0_xy, g_xyyz_xz_0_xz, g_xyyz_xz_0_yy, g_xyyz_xz_0_yz, g_xyyz_xz_0_zz, g_xz_0_0_0_yy_xz_0_xx, g_xz_0_0_0_yy_xz_0_xy, g_xz_0_0_0_yy_xz_0_xz, g_xz_0_0_0_yy_xz_0_yy, g_xz_0_0_0_yy_xz_0_yz, g_xz_0_0_0_yy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_xz_0_xx[i] = 4.0 * g_xyyz_xz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xz_0_xy[i] = 4.0 * g_xyyz_xz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xz_0_xz[i] = 4.0 * g_xyyz_xz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xz_0_yy[i] = 4.0 * g_xyyz_xz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xz_0_yz[i] = 4.0 * g_xyyz_xz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_xz_0_zz[i] = 4.0 * g_xyyz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (558-564)

    #pragma omp simd aligned(g_xyyz_yy_0_xx, g_xyyz_yy_0_xy, g_xyyz_yy_0_xz, g_xyyz_yy_0_yy, g_xyyz_yy_0_yz, g_xyyz_yy_0_zz, g_xz_0_0_0_yy_yy_0_xx, g_xz_0_0_0_yy_yy_0_xy, g_xz_0_0_0_yy_yy_0_xz, g_xz_0_0_0_yy_yy_0_yy, g_xz_0_0_0_yy_yy_0_yz, g_xz_0_0_0_yy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_yy_0_xx[i] = 4.0 * g_xyyz_yy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yy_0_xy[i] = 4.0 * g_xyyz_yy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yy_0_xz[i] = 4.0 * g_xyyz_yy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yy_0_yy[i] = 4.0 * g_xyyz_yy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yy_0_yz[i] = 4.0 * g_xyyz_yy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yy_0_zz[i] = 4.0 * g_xyyz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (564-570)

    #pragma omp simd aligned(g_xyyz_yz_0_xx, g_xyyz_yz_0_xy, g_xyyz_yz_0_xz, g_xyyz_yz_0_yy, g_xyyz_yz_0_yz, g_xyyz_yz_0_zz, g_xz_0_0_0_yy_yz_0_xx, g_xz_0_0_0_yy_yz_0_xy, g_xz_0_0_0_yy_yz_0_xz, g_xz_0_0_0_yy_yz_0_yy, g_xz_0_0_0_yy_yz_0_yz, g_xz_0_0_0_yy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_yz_0_xx[i] = 4.0 * g_xyyz_yz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yz_0_xy[i] = 4.0 * g_xyyz_yz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yz_0_xz[i] = 4.0 * g_xyyz_yz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yz_0_yy[i] = 4.0 * g_xyyz_yz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yz_0_yz[i] = 4.0 * g_xyyz_yz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_yz_0_zz[i] = 4.0 * g_xyyz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (570-576)

    #pragma omp simd aligned(g_xyyz_zz_0_xx, g_xyyz_zz_0_xy, g_xyyz_zz_0_xz, g_xyyz_zz_0_yy, g_xyyz_zz_0_yz, g_xyyz_zz_0_zz, g_xz_0_0_0_yy_zz_0_xx, g_xz_0_0_0_yy_zz_0_xy, g_xz_0_0_0_yy_zz_0_xz, g_xz_0_0_0_yy_zz_0_yy, g_xz_0_0_0_yy_zz_0_yz, g_xz_0_0_0_yy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yy_zz_0_xx[i] = 4.0 * g_xyyz_zz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_zz_0_xy[i] = 4.0 * g_xyyz_zz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_zz_0_xz[i] = 4.0 * g_xyyz_zz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_zz_0_yy[i] = 4.0 * g_xyyz_zz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_zz_0_yz[i] = 4.0 * g_xyyz_zz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_yy_zz_0_zz[i] = 4.0 * g_xyyz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (576-582)

    #pragma omp simd aligned(g_xy_xx_0_xx, g_xy_xx_0_xy, g_xy_xx_0_xz, g_xy_xx_0_yy, g_xy_xx_0_yz, g_xy_xx_0_zz, g_xyzz_xx_0_xx, g_xyzz_xx_0_xy, g_xyzz_xx_0_xz, g_xyzz_xx_0_yy, g_xyzz_xx_0_yz, g_xyzz_xx_0_zz, g_xz_0_0_0_yz_xx_0_xx, g_xz_0_0_0_yz_xx_0_xy, g_xz_0_0_0_yz_xx_0_xz, g_xz_0_0_0_yz_xx_0_yy, g_xz_0_0_0_yz_xx_0_yz, g_xz_0_0_0_yz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_xx_0_xx[i] = -2.0 * g_xy_xx_0_xx[i] * a_exp + 4.0 * g_xyzz_xx_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xx_0_xy[i] = -2.0 * g_xy_xx_0_xy[i] * a_exp + 4.0 * g_xyzz_xx_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xx_0_xz[i] = -2.0 * g_xy_xx_0_xz[i] * a_exp + 4.0 * g_xyzz_xx_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xx_0_yy[i] = -2.0 * g_xy_xx_0_yy[i] * a_exp + 4.0 * g_xyzz_xx_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xx_0_yz[i] = -2.0 * g_xy_xx_0_yz[i] * a_exp + 4.0 * g_xyzz_xx_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xx_0_zz[i] = -2.0 * g_xy_xx_0_zz[i] * a_exp + 4.0 * g_xyzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (582-588)

    #pragma omp simd aligned(g_xy_xy_0_xx, g_xy_xy_0_xy, g_xy_xy_0_xz, g_xy_xy_0_yy, g_xy_xy_0_yz, g_xy_xy_0_zz, g_xyzz_xy_0_xx, g_xyzz_xy_0_xy, g_xyzz_xy_0_xz, g_xyzz_xy_0_yy, g_xyzz_xy_0_yz, g_xyzz_xy_0_zz, g_xz_0_0_0_yz_xy_0_xx, g_xz_0_0_0_yz_xy_0_xy, g_xz_0_0_0_yz_xy_0_xz, g_xz_0_0_0_yz_xy_0_yy, g_xz_0_0_0_yz_xy_0_yz, g_xz_0_0_0_yz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_xy_0_xx[i] = -2.0 * g_xy_xy_0_xx[i] * a_exp + 4.0 * g_xyzz_xy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xy_0_xy[i] = -2.0 * g_xy_xy_0_xy[i] * a_exp + 4.0 * g_xyzz_xy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xy_0_xz[i] = -2.0 * g_xy_xy_0_xz[i] * a_exp + 4.0 * g_xyzz_xy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xy_0_yy[i] = -2.0 * g_xy_xy_0_yy[i] * a_exp + 4.0 * g_xyzz_xy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xy_0_yz[i] = -2.0 * g_xy_xy_0_yz[i] * a_exp + 4.0 * g_xyzz_xy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xy_0_zz[i] = -2.0 * g_xy_xy_0_zz[i] * a_exp + 4.0 * g_xyzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (588-594)

    #pragma omp simd aligned(g_xy_xz_0_xx, g_xy_xz_0_xy, g_xy_xz_0_xz, g_xy_xz_0_yy, g_xy_xz_0_yz, g_xy_xz_0_zz, g_xyzz_xz_0_xx, g_xyzz_xz_0_xy, g_xyzz_xz_0_xz, g_xyzz_xz_0_yy, g_xyzz_xz_0_yz, g_xyzz_xz_0_zz, g_xz_0_0_0_yz_xz_0_xx, g_xz_0_0_0_yz_xz_0_xy, g_xz_0_0_0_yz_xz_0_xz, g_xz_0_0_0_yz_xz_0_yy, g_xz_0_0_0_yz_xz_0_yz, g_xz_0_0_0_yz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_xz_0_xx[i] = -2.0 * g_xy_xz_0_xx[i] * a_exp + 4.0 * g_xyzz_xz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xz_0_xy[i] = -2.0 * g_xy_xz_0_xy[i] * a_exp + 4.0 * g_xyzz_xz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xz_0_xz[i] = -2.0 * g_xy_xz_0_xz[i] * a_exp + 4.0 * g_xyzz_xz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xz_0_yy[i] = -2.0 * g_xy_xz_0_yy[i] * a_exp + 4.0 * g_xyzz_xz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xz_0_yz[i] = -2.0 * g_xy_xz_0_yz[i] * a_exp + 4.0 * g_xyzz_xz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_xz_0_zz[i] = -2.0 * g_xy_xz_0_zz[i] * a_exp + 4.0 * g_xyzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (594-600)

    #pragma omp simd aligned(g_xy_yy_0_xx, g_xy_yy_0_xy, g_xy_yy_0_xz, g_xy_yy_0_yy, g_xy_yy_0_yz, g_xy_yy_0_zz, g_xyzz_yy_0_xx, g_xyzz_yy_0_xy, g_xyzz_yy_0_xz, g_xyzz_yy_0_yy, g_xyzz_yy_0_yz, g_xyzz_yy_0_zz, g_xz_0_0_0_yz_yy_0_xx, g_xz_0_0_0_yz_yy_0_xy, g_xz_0_0_0_yz_yy_0_xz, g_xz_0_0_0_yz_yy_0_yy, g_xz_0_0_0_yz_yy_0_yz, g_xz_0_0_0_yz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_yy_0_xx[i] = -2.0 * g_xy_yy_0_xx[i] * a_exp + 4.0 * g_xyzz_yy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yy_0_xy[i] = -2.0 * g_xy_yy_0_xy[i] * a_exp + 4.0 * g_xyzz_yy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yy_0_xz[i] = -2.0 * g_xy_yy_0_xz[i] * a_exp + 4.0 * g_xyzz_yy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yy_0_yy[i] = -2.0 * g_xy_yy_0_yy[i] * a_exp + 4.0 * g_xyzz_yy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yy_0_yz[i] = -2.0 * g_xy_yy_0_yz[i] * a_exp + 4.0 * g_xyzz_yy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yy_0_zz[i] = -2.0 * g_xy_yy_0_zz[i] * a_exp + 4.0 * g_xyzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (600-606)

    #pragma omp simd aligned(g_xy_yz_0_xx, g_xy_yz_0_xy, g_xy_yz_0_xz, g_xy_yz_0_yy, g_xy_yz_0_yz, g_xy_yz_0_zz, g_xyzz_yz_0_xx, g_xyzz_yz_0_xy, g_xyzz_yz_0_xz, g_xyzz_yz_0_yy, g_xyzz_yz_0_yz, g_xyzz_yz_0_zz, g_xz_0_0_0_yz_yz_0_xx, g_xz_0_0_0_yz_yz_0_xy, g_xz_0_0_0_yz_yz_0_xz, g_xz_0_0_0_yz_yz_0_yy, g_xz_0_0_0_yz_yz_0_yz, g_xz_0_0_0_yz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_yz_0_xx[i] = -2.0 * g_xy_yz_0_xx[i] * a_exp + 4.0 * g_xyzz_yz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yz_0_xy[i] = -2.0 * g_xy_yz_0_xy[i] * a_exp + 4.0 * g_xyzz_yz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yz_0_xz[i] = -2.0 * g_xy_yz_0_xz[i] * a_exp + 4.0 * g_xyzz_yz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yz_0_yy[i] = -2.0 * g_xy_yz_0_yy[i] * a_exp + 4.0 * g_xyzz_yz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yz_0_yz[i] = -2.0 * g_xy_yz_0_yz[i] * a_exp + 4.0 * g_xyzz_yz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_yz_0_zz[i] = -2.0 * g_xy_yz_0_zz[i] * a_exp + 4.0 * g_xyzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (606-612)

    #pragma omp simd aligned(g_xy_zz_0_xx, g_xy_zz_0_xy, g_xy_zz_0_xz, g_xy_zz_0_yy, g_xy_zz_0_yz, g_xy_zz_0_zz, g_xyzz_zz_0_xx, g_xyzz_zz_0_xy, g_xyzz_zz_0_xz, g_xyzz_zz_0_yy, g_xyzz_zz_0_yz, g_xyzz_zz_0_zz, g_xz_0_0_0_yz_zz_0_xx, g_xz_0_0_0_yz_zz_0_xy, g_xz_0_0_0_yz_zz_0_xz, g_xz_0_0_0_yz_zz_0_yy, g_xz_0_0_0_yz_zz_0_yz, g_xz_0_0_0_yz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_yz_zz_0_xx[i] = -2.0 * g_xy_zz_0_xx[i] * a_exp + 4.0 * g_xyzz_zz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_zz_0_xy[i] = -2.0 * g_xy_zz_0_xy[i] * a_exp + 4.0 * g_xyzz_zz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_zz_0_xz[i] = -2.0 * g_xy_zz_0_xz[i] * a_exp + 4.0 * g_xyzz_zz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_zz_0_yy[i] = -2.0 * g_xy_zz_0_yy[i] * a_exp + 4.0 * g_xyzz_zz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_zz_0_yz[i] = -2.0 * g_xy_zz_0_yz[i] * a_exp + 4.0 * g_xyzz_zz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_yz_zz_0_zz[i] = -2.0 * g_xy_zz_0_zz[i] * a_exp + 4.0 * g_xyzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (612-618)

    #pragma omp simd aligned(g_xz_0_0_0_zz_xx_0_xx, g_xz_0_0_0_zz_xx_0_xy, g_xz_0_0_0_zz_xx_0_xz, g_xz_0_0_0_zz_xx_0_yy, g_xz_0_0_0_zz_xx_0_yz, g_xz_0_0_0_zz_xx_0_zz, g_xz_xx_0_xx, g_xz_xx_0_xy, g_xz_xx_0_xz, g_xz_xx_0_yy, g_xz_xx_0_yz, g_xz_xx_0_zz, g_xzzz_xx_0_xx, g_xzzz_xx_0_xy, g_xzzz_xx_0_xz, g_xzzz_xx_0_yy, g_xzzz_xx_0_yz, g_xzzz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_xx_0_xx[i] = -4.0 * g_xz_xx_0_xx[i] * a_exp + 4.0 * g_xzzz_xx_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xx_0_xy[i] = -4.0 * g_xz_xx_0_xy[i] * a_exp + 4.0 * g_xzzz_xx_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xx_0_xz[i] = -4.0 * g_xz_xx_0_xz[i] * a_exp + 4.0 * g_xzzz_xx_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xx_0_yy[i] = -4.0 * g_xz_xx_0_yy[i] * a_exp + 4.0 * g_xzzz_xx_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xx_0_yz[i] = -4.0 * g_xz_xx_0_yz[i] * a_exp + 4.0 * g_xzzz_xx_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xx_0_zz[i] = -4.0 * g_xz_xx_0_zz[i] * a_exp + 4.0 * g_xzzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (618-624)

    #pragma omp simd aligned(g_xz_0_0_0_zz_xy_0_xx, g_xz_0_0_0_zz_xy_0_xy, g_xz_0_0_0_zz_xy_0_xz, g_xz_0_0_0_zz_xy_0_yy, g_xz_0_0_0_zz_xy_0_yz, g_xz_0_0_0_zz_xy_0_zz, g_xz_xy_0_xx, g_xz_xy_0_xy, g_xz_xy_0_xz, g_xz_xy_0_yy, g_xz_xy_0_yz, g_xz_xy_0_zz, g_xzzz_xy_0_xx, g_xzzz_xy_0_xy, g_xzzz_xy_0_xz, g_xzzz_xy_0_yy, g_xzzz_xy_0_yz, g_xzzz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_xy_0_xx[i] = -4.0 * g_xz_xy_0_xx[i] * a_exp + 4.0 * g_xzzz_xy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xy_0_xy[i] = -4.0 * g_xz_xy_0_xy[i] * a_exp + 4.0 * g_xzzz_xy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xy_0_xz[i] = -4.0 * g_xz_xy_0_xz[i] * a_exp + 4.0 * g_xzzz_xy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xy_0_yy[i] = -4.0 * g_xz_xy_0_yy[i] * a_exp + 4.0 * g_xzzz_xy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xy_0_yz[i] = -4.0 * g_xz_xy_0_yz[i] * a_exp + 4.0 * g_xzzz_xy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xy_0_zz[i] = -4.0 * g_xz_xy_0_zz[i] * a_exp + 4.0 * g_xzzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (624-630)

    #pragma omp simd aligned(g_xz_0_0_0_zz_xz_0_xx, g_xz_0_0_0_zz_xz_0_xy, g_xz_0_0_0_zz_xz_0_xz, g_xz_0_0_0_zz_xz_0_yy, g_xz_0_0_0_zz_xz_0_yz, g_xz_0_0_0_zz_xz_0_zz, g_xz_xz_0_xx, g_xz_xz_0_xy, g_xz_xz_0_xz, g_xz_xz_0_yy, g_xz_xz_0_yz, g_xz_xz_0_zz, g_xzzz_xz_0_xx, g_xzzz_xz_0_xy, g_xzzz_xz_0_xz, g_xzzz_xz_0_yy, g_xzzz_xz_0_yz, g_xzzz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_xz_0_xx[i] = -4.0 * g_xz_xz_0_xx[i] * a_exp + 4.0 * g_xzzz_xz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xz_0_xy[i] = -4.0 * g_xz_xz_0_xy[i] * a_exp + 4.0 * g_xzzz_xz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xz_0_xz[i] = -4.0 * g_xz_xz_0_xz[i] * a_exp + 4.0 * g_xzzz_xz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xz_0_yy[i] = -4.0 * g_xz_xz_0_yy[i] * a_exp + 4.0 * g_xzzz_xz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xz_0_yz[i] = -4.0 * g_xz_xz_0_yz[i] * a_exp + 4.0 * g_xzzz_xz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_xz_0_zz[i] = -4.0 * g_xz_xz_0_zz[i] * a_exp + 4.0 * g_xzzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (630-636)

    #pragma omp simd aligned(g_xz_0_0_0_zz_yy_0_xx, g_xz_0_0_0_zz_yy_0_xy, g_xz_0_0_0_zz_yy_0_xz, g_xz_0_0_0_zz_yy_0_yy, g_xz_0_0_0_zz_yy_0_yz, g_xz_0_0_0_zz_yy_0_zz, g_xz_yy_0_xx, g_xz_yy_0_xy, g_xz_yy_0_xz, g_xz_yy_0_yy, g_xz_yy_0_yz, g_xz_yy_0_zz, g_xzzz_yy_0_xx, g_xzzz_yy_0_xy, g_xzzz_yy_0_xz, g_xzzz_yy_0_yy, g_xzzz_yy_0_yz, g_xzzz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_yy_0_xx[i] = -4.0 * g_xz_yy_0_xx[i] * a_exp + 4.0 * g_xzzz_yy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yy_0_xy[i] = -4.0 * g_xz_yy_0_xy[i] * a_exp + 4.0 * g_xzzz_yy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yy_0_xz[i] = -4.0 * g_xz_yy_0_xz[i] * a_exp + 4.0 * g_xzzz_yy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yy_0_yy[i] = -4.0 * g_xz_yy_0_yy[i] * a_exp + 4.0 * g_xzzz_yy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yy_0_yz[i] = -4.0 * g_xz_yy_0_yz[i] * a_exp + 4.0 * g_xzzz_yy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yy_0_zz[i] = -4.0 * g_xz_yy_0_zz[i] * a_exp + 4.0 * g_xzzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (636-642)

    #pragma omp simd aligned(g_xz_0_0_0_zz_yz_0_xx, g_xz_0_0_0_zz_yz_0_xy, g_xz_0_0_0_zz_yz_0_xz, g_xz_0_0_0_zz_yz_0_yy, g_xz_0_0_0_zz_yz_0_yz, g_xz_0_0_0_zz_yz_0_zz, g_xz_yz_0_xx, g_xz_yz_0_xy, g_xz_yz_0_xz, g_xz_yz_0_yy, g_xz_yz_0_yz, g_xz_yz_0_zz, g_xzzz_yz_0_xx, g_xzzz_yz_0_xy, g_xzzz_yz_0_xz, g_xzzz_yz_0_yy, g_xzzz_yz_0_yz, g_xzzz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_yz_0_xx[i] = -4.0 * g_xz_yz_0_xx[i] * a_exp + 4.0 * g_xzzz_yz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yz_0_xy[i] = -4.0 * g_xz_yz_0_xy[i] * a_exp + 4.0 * g_xzzz_yz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yz_0_xz[i] = -4.0 * g_xz_yz_0_xz[i] * a_exp + 4.0 * g_xzzz_yz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yz_0_yy[i] = -4.0 * g_xz_yz_0_yy[i] * a_exp + 4.0 * g_xzzz_yz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yz_0_yz[i] = -4.0 * g_xz_yz_0_yz[i] * a_exp + 4.0 * g_xzzz_yz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_yz_0_zz[i] = -4.0 * g_xz_yz_0_zz[i] * a_exp + 4.0 * g_xzzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (642-648)

    #pragma omp simd aligned(g_xz_0_0_0_zz_zz_0_xx, g_xz_0_0_0_zz_zz_0_xy, g_xz_0_0_0_zz_zz_0_xz, g_xz_0_0_0_zz_zz_0_yy, g_xz_0_0_0_zz_zz_0_yz, g_xz_0_0_0_zz_zz_0_zz, g_xz_zz_0_xx, g_xz_zz_0_xy, g_xz_zz_0_xz, g_xz_zz_0_yy, g_xz_zz_0_yz, g_xz_zz_0_zz, g_xzzz_zz_0_xx, g_xzzz_zz_0_xy, g_xzzz_zz_0_xz, g_xzzz_zz_0_yy, g_xzzz_zz_0_yz, g_xzzz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_zz_zz_0_xx[i] = -4.0 * g_xz_zz_0_xx[i] * a_exp + 4.0 * g_xzzz_zz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_zz_0_xy[i] = -4.0 * g_xz_zz_0_xy[i] * a_exp + 4.0 * g_xzzz_zz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_zz_0_xz[i] = -4.0 * g_xz_zz_0_xz[i] * a_exp + 4.0 * g_xzzz_zz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_zz_0_yy[i] = -4.0 * g_xz_zz_0_yy[i] * a_exp + 4.0 * g_xzzz_zz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_zz_0_yz[i] = -4.0 * g_xz_zz_0_yz[i] * a_exp + 4.0 * g_xzzz_zz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_zz_zz_0_zz[i] = -4.0 * g_xz_zz_0_zz[i] * a_exp + 4.0 * g_xzzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (648-654)

    #pragma omp simd aligned(g_xx_xx_0_xx, g_xx_xx_0_xy, g_xx_xx_0_xz, g_xx_xx_0_yy, g_xx_xx_0_yz, g_xx_xx_0_zz, g_xxyy_xx_0_xx, g_xxyy_xx_0_xy, g_xxyy_xx_0_xz, g_xxyy_xx_0_yy, g_xxyy_xx_0_yz, g_xxyy_xx_0_zz, g_yy_0_0_0_xx_xx_0_xx, g_yy_0_0_0_xx_xx_0_xy, g_yy_0_0_0_xx_xx_0_xz, g_yy_0_0_0_xx_xx_0_yy, g_yy_0_0_0_xx_xx_0_yz, g_yy_0_0_0_xx_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_xx_0_xx[i] = -2.0 * g_xx_xx_0_xx[i] * a_exp + 4.0 * g_xxyy_xx_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xx_0_xy[i] = -2.0 * g_xx_xx_0_xy[i] * a_exp + 4.0 * g_xxyy_xx_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xx_0_xz[i] = -2.0 * g_xx_xx_0_xz[i] * a_exp + 4.0 * g_xxyy_xx_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xx_0_yy[i] = -2.0 * g_xx_xx_0_yy[i] * a_exp + 4.0 * g_xxyy_xx_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xx_0_yz[i] = -2.0 * g_xx_xx_0_yz[i] * a_exp + 4.0 * g_xxyy_xx_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xx_0_zz[i] = -2.0 * g_xx_xx_0_zz[i] * a_exp + 4.0 * g_xxyy_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (654-660)

    #pragma omp simd aligned(g_xx_xy_0_xx, g_xx_xy_0_xy, g_xx_xy_0_xz, g_xx_xy_0_yy, g_xx_xy_0_yz, g_xx_xy_0_zz, g_xxyy_xy_0_xx, g_xxyy_xy_0_xy, g_xxyy_xy_0_xz, g_xxyy_xy_0_yy, g_xxyy_xy_0_yz, g_xxyy_xy_0_zz, g_yy_0_0_0_xx_xy_0_xx, g_yy_0_0_0_xx_xy_0_xy, g_yy_0_0_0_xx_xy_0_xz, g_yy_0_0_0_xx_xy_0_yy, g_yy_0_0_0_xx_xy_0_yz, g_yy_0_0_0_xx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_xy_0_xx[i] = -2.0 * g_xx_xy_0_xx[i] * a_exp + 4.0 * g_xxyy_xy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xy_0_xy[i] = -2.0 * g_xx_xy_0_xy[i] * a_exp + 4.0 * g_xxyy_xy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xy_0_xz[i] = -2.0 * g_xx_xy_0_xz[i] * a_exp + 4.0 * g_xxyy_xy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xy_0_yy[i] = -2.0 * g_xx_xy_0_yy[i] * a_exp + 4.0 * g_xxyy_xy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xy_0_yz[i] = -2.0 * g_xx_xy_0_yz[i] * a_exp + 4.0 * g_xxyy_xy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xy_0_zz[i] = -2.0 * g_xx_xy_0_zz[i] * a_exp + 4.0 * g_xxyy_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (660-666)

    #pragma omp simd aligned(g_xx_xz_0_xx, g_xx_xz_0_xy, g_xx_xz_0_xz, g_xx_xz_0_yy, g_xx_xz_0_yz, g_xx_xz_0_zz, g_xxyy_xz_0_xx, g_xxyy_xz_0_xy, g_xxyy_xz_0_xz, g_xxyy_xz_0_yy, g_xxyy_xz_0_yz, g_xxyy_xz_0_zz, g_yy_0_0_0_xx_xz_0_xx, g_yy_0_0_0_xx_xz_0_xy, g_yy_0_0_0_xx_xz_0_xz, g_yy_0_0_0_xx_xz_0_yy, g_yy_0_0_0_xx_xz_0_yz, g_yy_0_0_0_xx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_xz_0_xx[i] = -2.0 * g_xx_xz_0_xx[i] * a_exp + 4.0 * g_xxyy_xz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xz_0_xy[i] = -2.0 * g_xx_xz_0_xy[i] * a_exp + 4.0 * g_xxyy_xz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xz_0_xz[i] = -2.0 * g_xx_xz_0_xz[i] * a_exp + 4.0 * g_xxyy_xz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xz_0_yy[i] = -2.0 * g_xx_xz_0_yy[i] * a_exp + 4.0 * g_xxyy_xz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xz_0_yz[i] = -2.0 * g_xx_xz_0_yz[i] * a_exp + 4.0 * g_xxyy_xz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_xz_0_zz[i] = -2.0 * g_xx_xz_0_zz[i] * a_exp + 4.0 * g_xxyy_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (666-672)

    #pragma omp simd aligned(g_xx_yy_0_xx, g_xx_yy_0_xy, g_xx_yy_0_xz, g_xx_yy_0_yy, g_xx_yy_0_yz, g_xx_yy_0_zz, g_xxyy_yy_0_xx, g_xxyy_yy_0_xy, g_xxyy_yy_0_xz, g_xxyy_yy_0_yy, g_xxyy_yy_0_yz, g_xxyy_yy_0_zz, g_yy_0_0_0_xx_yy_0_xx, g_yy_0_0_0_xx_yy_0_xy, g_yy_0_0_0_xx_yy_0_xz, g_yy_0_0_0_xx_yy_0_yy, g_yy_0_0_0_xx_yy_0_yz, g_yy_0_0_0_xx_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_yy_0_xx[i] = -2.0 * g_xx_yy_0_xx[i] * a_exp + 4.0 * g_xxyy_yy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yy_0_xy[i] = -2.0 * g_xx_yy_0_xy[i] * a_exp + 4.0 * g_xxyy_yy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yy_0_xz[i] = -2.0 * g_xx_yy_0_xz[i] * a_exp + 4.0 * g_xxyy_yy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yy_0_yy[i] = -2.0 * g_xx_yy_0_yy[i] * a_exp + 4.0 * g_xxyy_yy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yy_0_yz[i] = -2.0 * g_xx_yy_0_yz[i] * a_exp + 4.0 * g_xxyy_yy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yy_0_zz[i] = -2.0 * g_xx_yy_0_zz[i] * a_exp + 4.0 * g_xxyy_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (672-678)

    #pragma omp simd aligned(g_xx_yz_0_xx, g_xx_yz_0_xy, g_xx_yz_0_xz, g_xx_yz_0_yy, g_xx_yz_0_yz, g_xx_yz_0_zz, g_xxyy_yz_0_xx, g_xxyy_yz_0_xy, g_xxyy_yz_0_xz, g_xxyy_yz_0_yy, g_xxyy_yz_0_yz, g_xxyy_yz_0_zz, g_yy_0_0_0_xx_yz_0_xx, g_yy_0_0_0_xx_yz_0_xy, g_yy_0_0_0_xx_yz_0_xz, g_yy_0_0_0_xx_yz_0_yy, g_yy_0_0_0_xx_yz_0_yz, g_yy_0_0_0_xx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_yz_0_xx[i] = -2.0 * g_xx_yz_0_xx[i] * a_exp + 4.0 * g_xxyy_yz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yz_0_xy[i] = -2.0 * g_xx_yz_0_xy[i] * a_exp + 4.0 * g_xxyy_yz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yz_0_xz[i] = -2.0 * g_xx_yz_0_xz[i] * a_exp + 4.0 * g_xxyy_yz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yz_0_yy[i] = -2.0 * g_xx_yz_0_yy[i] * a_exp + 4.0 * g_xxyy_yz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yz_0_yz[i] = -2.0 * g_xx_yz_0_yz[i] * a_exp + 4.0 * g_xxyy_yz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_yz_0_zz[i] = -2.0 * g_xx_yz_0_zz[i] * a_exp + 4.0 * g_xxyy_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (678-684)

    #pragma omp simd aligned(g_xx_zz_0_xx, g_xx_zz_0_xy, g_xx_zz_0_xz, g_xx_zz_0_yy, g_xx_zz_0_yz, g_xx_zz_0_zz, g_xxyy_zz_0_xx, g_xxyy_zz_0_xy, g_xxyy_zz_0_xz, g_xxyy_zz_0_yy, g_xxyy_zz_0_yz, g_xxyy_zz_0_zz, g_yy_0_0_0_xx_zz_0_xx, g_yy_0_0_0_xx_zz_0_xy, g_yy_0_0_0_xx_zz_0_xz, g_yy_0_0_0_xx_zz_0_yy, g_yy_0_0_0_xx_zz_0_yz, g_yy_0_0_0_xx_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xx_zz_0_xx[i] = -2.0 * g_xx_zz_0_xx[i] * a_exp + 4.0 * g_xxyy_zz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_zz_0_xy[i] = -2.0 * g_xx_zz_0_xy[i] * a_exp + 4.0 * g_xxyy_zz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_zz_0_xz[i] = -2.0 * g_xx_zz_0_xz[i] * a_exp + 4.0 * g_xxyy_zz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_zz_0_yy[i] = -2.0 * g_xx_zz_0_yy[i] * a_exp + 4.0 * g_xxyy_zz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_zz_0_yz[i] = -2.0 * g_xx_zz_0_yz[i] * a_exp + 4.0 * g_xxyy_zz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xx_zz_0_zz[i] = -2.0 * g_xx_zz_0_zz[i] * a_exp + 4.0 * g_xxyy_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (684-690)

    #pragma omp simd aligned(g_xy_xx_0_xx, g_xy_xx_0_xy, g_xy_xx_0_xz, g_xy_xx_0_yy, g_xy_xx_0_yz, g_xy_xx_0_zz, g_xyyy_xx_0_xx, g_xyyy_xx_0_xy, g_xyyy_xx_0_xz, g_xyyy_xx_0_yy, g_xyyy_xx_0_yz, g_xyyy_xx_0_zz, g_yy_0_0_0_xy_xx_0_xx, g_yy_0_0_0_xy_xx_0_xy, g_yy_0_0_0_xy_xx_0_xz, g_yy_0_0_0_xy_xx_0_yy, g_yy_0_0_0_xy_xx_0_yz, g_yy_0_0_0_xy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_xx_0_xx[i] = -6.0 * g_xy_xx_0_xx[i] * a_exp + 4.0 * g_xyyy_xx_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xx_0_xy[i] = -6.0 * g_xy_xx_0_xy[i] * a_exp + 4.0 * g_xyyy_xx_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xx_0_xz[i] = -6.0 * g_xy_xx_0_xz[i] * a_exp + 4.0 * g_xyyy_xx_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xx_0_yy[i] = -6.0 * g_xy_xx_0_yy[i] * a_exp + 4.0 * g_xyyy_xx_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xx_0_yz[i] = -6.0 * g_xy_xx_0_yz[i] * a_exp + 4.0 * g_xyyy_xx_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xx_0_zz[i] = -6.0 * g_xy_xx_0_zz[i] * a_exp + 4.0 * g_xyyy_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (690-696)

    #pragma omp simd aligned(g_xy_xy_0_xx, g_xy_xy_0_xy, g_xy_xy_0_xz, g_xy_xy_0_yy, g_xy_xy_0_yz, g_xy_xy_0_zz, g_xyyy_xy_0_xx, g_xyyy_xy_0_xy, g_xyyy_xy_0_xz, g_xyyy_xy_0_yy, g_xyyy_xy_0_yz, g_xyyy_xy_0_zz, g_yy_0_0_0_xy_xy_0_xx, g_yy_0_0_0_xy_xy_0_xy, g_yy_0_0_0_xy_xy_0_xz, g_yy_0_0_0_xy_xy_0_yy, g_yy_0_0_0_xy_xy_0_yz, g_yy_0_0_0_xy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_xy_0_xx[i] = -6.0 * g_xy_xy_0_xx[i] * a_exp + 4.0 * g_xyyy_xy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xy_0_xy[i] = -6.0 * g_xy_xy_0_xy[i] * a_exp + 4.0 * g_xyyy_xy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xy_0_xz[i] = -6.0 * g_xy_xy_0_xz[i] * a_exp + 4.0 * g_xyyy_xy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xy_0_yy[i] = -6.0 * g_xy_xy_0_yy[i] * a_exp + 4.0 * g_xyyy_xy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xy_0_yz[i] = -6.0 * g_xy_xy_0_yz[i] * a_exp + 4.0 * g_xyyy_xy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xy_0_zz[i] = -6.0 * g_xy_xy_0_zz[i] * a_exp + 4.0 * g_xyyy_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (696-702)

    #pragma omp simd aligned(g_xy_xz_0_xx, g_xy_xz_0_xy, g_xy_xz_0_xz, g_xy_xz_0_yy, g_xy_xz_0_yz, g_xy_xz_0_zz, g_xyyy_xz_0_xx, g_xyyy_xz_0_xy, g_xyyy_xz_0_xz, g_xyyy_xz_0_yy, g_xyyy_xz_0_yz, g_xyyy_xz_0_zz, g_yy_0_0_0_xy_xz_0_xx, g_yy_0_0_0_xy_xz_0_xy, g_yy_0_0_0_xy_xz_0_xz, g_yy_0_0_0_xy_xz_0_yy, g_yy_0_0_0_xy_xz_0_yz, g_yy_0_0_0_xy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_xz_0_xx[i] = -6.0 * g_xy_xz_0_xx[i] * a_exp + 4.0 * g_xyyy_xz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xz_0_xy[i] = -6.0 * g_xy_xz_0_xy[i] * a_exp + 4.0 * g_xyyy_xz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xz_0_xz[i] = -6.0 * g_xy_xz_0_xz[i] * a_exp + 4.0 * g_xyyy_xz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xz_0_yy[i] = -6.0 * g_xy_xz_0_yy[i] * a_exp + 4.0 * g_xyyy_xz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xz_0_yz[i] = -6.0 * g_xy_xz_0_yz[i] * a_exp + 4.0 * g_xyyy_xz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_xz_0_zz[i] = -6.0 * g_xy_xz_0_zz[i] * a_exp + 4.0 * g_xyyy_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (702-708)

    #pragma omp simd aligned(g_xy_yy_0_xx, g_xy_yy_0_xy, g_xy_yy_0_xz, g_xy_yy_0_yy, g_xy_yy_0_yz, g_xy_yy_0_zz, g_xyyy_yy_0_xx, g_xyyy_yy_0_xy, g_xyyy_yy_0_xz, g_xyyy_yy_0_yy, g_xyyy_yy_0_yz, g_xyyy_yy_0_zz, g_yy_0_0_0_xy_yy_0_xx, g_yy_0_0_0_xy_yy_0_xy, g_yy_0_0_0_xy_yy_0_xz, g_yy_0_0_0_xy_yy_0_yy, g_yy_0_0_0_xy_yy_0_yz, g_yy_0_0_0_xy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_yy_0_xx[i] = -6.0 * g_xy_yy_0_xx[i] * a_exp + 4.0 * g_xyyy_yy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yy_0_xy[i] = -6.0 * g_xy_yy_0_xy[i] * a_exp + 4.0 * g_xyyy_yy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yy_0_xz[i] = -6.0 * g_xy_yy_0_xz[i] * a_exp + 4.0 * g_xyyy_yy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yy_0_yy[i] = -6.0 * g_xy_yy_0_yy[i] * a_exp + 4.0 * g_xyyy_yy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yy_0_yz[i] = -6.0 * g_xy_yy_0_yz[i] * a_exp + 4.0 * g_xyyy_yy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yy_0_zz[i] = -6.0 * g_xy_yy_0_zz[i] * a_exp + 4.0 * g_xyyy_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (708-714)

    #pragma omp simd aligned(g_xy_yz_0_xx, g_xy_yz_0_xy, g_xy_yz_0_xz, g_xy_yz_0_yy, g_xy_yz_0_yz, g_xy_yz_0_zz, g_xyyy_yz_0_xx, g_xyyy_yz_0_xy, g_xyyy_yz_0_xz, g_xyyy_yz_0_yy, g_xyyy_yz_0_yz, g_xyyy_yz_0_zz, g_yy_0_0_0_xy_yz_0_xx, g_yy_0_0_0_xy_yz_0_xy, g_yy_0_0_0_xy_yz_0_xz, g_yy_0_0_0_xy_yz_0_yy, g_yy_0_0_0_xy_yz_0_yz, g_yy_0_0_0_xy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_yz_0_xx[i] = -6.0 * g_xy_yz_0_xx[i] * a_exp + 4.0 * g_xyyy_yz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yz_0_xy[i] = -6.0 * g_xy_yz_0_xy[i] * a_exp + 4.0 * g_xyyy_yz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yz_0_xz[i] = -6.0 * g_xy_yz_0_xz[i] * a_exp + 4.0 * g_xyyy_yz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yz_0_yy[i] = -6.0 * g_xy_yz_0_yy[i] * a_exp + 4.0 * g_xyyy_yz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yz_0_yz[i] = -6.0 * g_xy_yz_0_yz[i] * a_exp + 4.0 * g_xyyy_yz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_yz_0_zz[i] = -6.0 * g_xy_yz_0_zz[i] * a_exp + 4.0 * g_xyyy_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (714-720)

    #pragma omp simd aligned(g_xy_zz_0_xx, g_xy_zz_0_xy, g_xy_zz_0_xz, g_xy_zz_0_yy, g_xy_zz_0_yz, g_xy_zz_0_zz, g_xyyy_zz_0_xx, g_xyyy_zz_0_xy, g_xyyy_zz_0_xz, g_xyyy_zz_0_yy, g_xyyy_zz_0_yz, g_xyyy_zz_0_zz, g_yy_0_0_0_xy_zz_0_xx, g_yy_0_0_0_xy_zz_0_xy, g_yy_0_0_0_xy_zz_0_xz, g_yy_0_0_0_xy_zz_0_yy, g_yy_0_0_0_xy_zz_0_yz, g_yy_0_0_0_xy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xy_zz_0_xx[i] = -6.0 * g_xy_zz_0_xx[i] * a_exp + 4.0 * g_xyyy_zz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_zz_0_xy[i] = -6.0 * g_xy_zz_0_xy[i] * a_exp + 4.0 * g_xyyy_zz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_zz_0_xz[i] = -6.0 * g_xy_zz_0_xz[i] * a_exp + 4.0 * g_xyyy_zz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_zz_0_yy[i] = -6.0 * g_xy_zz_0_yy[i] * a_exp + 4.0 * g_xyyy_zz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_zz_0_yz[i] = -6.0 * g_xy_zz_0_yz[i] * a_exp + 4.0 * g_xyyy_zz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xy_zz_0_zz[i] = -6.0 * g_xy_zz_0_zz[i] * a_exp + 4.0 * g_xyyy_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (720-726)

    #pragma omp simd aligned(g_xyyz_xx_0_xx, g_xyyz_xx_0_xy, g_xyyz_xx_0_xz, g_xyyz_xx_0_yy, g_xyyz_xx_0_yz, g_xyyz_xx_0_zz, g_xz_xx_0_xx, g_xz_xx_0_xy, g_xz_xx_0_xz, g_xz_xx_0_yy, g_xz_xx_0_yz, g_xz_xx_0_zz, g_yy_0_0_0_xz_xx_0_xx, g_yy_0_0_0_xz_xx_0_xy, g_yy_0_0_0_xz_xx_0_xz, g_yy_0_0_0_xz_xx_0_yy, g_yy_0_0_0_xz_xx_0_yz, g_yy_0_0_0_xz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_xx_0_xx[i] = -2.0 * g_xz_xx_0_xx[i] * a_exp + 4.0 * g_xyyz_xx_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xx_0_xy[i] = -2.0 * g_xz_xx_0_xy[i] * a_exp + 4.0 * g_xyyz_xx_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xx_0_xz[i] = -2.0 * g_xz_xx_0_xz[i] * a_exp + 4.0 * g_xyyz_xx_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xx_0_yy[i] = -2.0 * g_xz_xx_0_yy[i] * a_exp + 4.0 * g_xyyz_xx_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xx_0_yz[i] = -2.0 * g_xz_xx_0_yz[i] * a_exp + 4.0 * g_xyyz_xx_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xx_0_zz[i] = -2.0 * g_xz_xx_0_zz[i] * a_exp + 4.0 * g_xyyz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (726-732)

    #pragma omp simd aligned(g_xyyz_xy_0_xx, g_xyyz_xy_0_xy, g_xyyz_xy_0_xz, g_xyyz_xy_0_yy, g_xyyz_xy_0_yz, g_xyyz_xy_0_zz, g_xz_xy_0_xx, g_xz_xy_0_xy, g_xz_xy_0_xz, g_xz_xy_0_yy, g_xz_xy_0_yz, g_xz_xy_0_zz, g_yy_0_0_0_xz_xy_0_xx, g_yy_0_0_0_xz_xy_0_xy, g_yy_0_0_0_xz_xy_0_xz, g_yy_0_0_0_xz_xy_0_yy, g_yy_0_0_0_xz_xy_0_yz, g_yy_0_0_0_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_xy_0_xx[i] = -2.0 * g_xz_xy_0_xx[i] * a_exp + 4.0 * g_xyyz_xy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xy_0_xy[i] = -2.0 * g_xz_xy_0_xy[i] * a_exp + 4.0 * g_xyyz_xy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xy_0_xz[i] = -2.0 * g_xz_xy_0_xz[i] * a_exp + 4.0 * g_xyyz_xy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xy_0_yy[i] = -2.0 * g_xz_xy_0_yy[i] * a_exp + 4.0 * g_xyyz_xy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xy_0_yz[i] = -2.0 * g_xz_xy_0_yz[i] * a_exp + 4.0 * g_xyyz_xy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xy_0_zz[i] = -2.0 * g_xz_xy_0_zz[i] * a_exp + 4.0 * g_xyyz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (732-738)

    #pragma omp simd aligned(g_xyyz_xz_0_xx, g_xyyz_xz_0_xy, g_xyyz_xz_0_xz, g_xyyz_xz_0_yy, g_xyyz_xz_0_yz, g_xyyz_xz_0_zz, g_xz_xz_0_xx, g_xz_xz_0_xy, g_xz_xz_0_xz, g_xz_xz_0_yy, g_xz_xz_0_yz, g_xz_xz_0_zz, g_yy_0_0_0_xz_xz_0_xx, g_yy_0_0_0_xz_xz_0_xy, g_yy_0_0_0_xz_xz_0_xz, g_yy_0_0_0_xz_xz_0_yy, g_yy_0_0_0_xz_xz_0_yz, g_yy_0_0_0_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_xz_0_xx[i] = -2.0 * g_xz_xz_0_xx[i] * a_exp + 4.0 * g_xyyz_xz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xz_0_xy[i] = -2.0 * g_xz_xz_0_xy[i] * a_exp + 4.0 * g_xyyz_xz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xz_0_xz[i] = -2.0 * g_xz_xz_0_xz[i] * a_exp + 4.0 * g_xyyz_xz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xz_0_yy[i] = -2.0 * g_xz_xz_0_yy[i] * a_exp + 4.0 * g_xyyz_xz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xz_0_yz[i] = -2.0 * g_xz_xz_0_yz[i] * a_exp + 4.0 * g_xyyz_xz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_xz_0_zz[i] = -2.0 * g_xz_xz_0_zz[i] * a_exp + 4.0 * g_xyyz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (738-744)

    #pragma omp simd aligned(g_xyyz_yy_0_xx, g_xyyz_yy_0_xy, g_xyyz_yy_0_xz, g_xyyz_yy_0_yy, g_xyyz_yy_0_yz, g_xyyz_yy_0_zz, g_xz_yy_0_xx, g_xz_yy_0_xy, g_xz_yy_0_xz, g_xz_yy_0_yy, g_xz_yy_0_yz, g_xz_yy_0_zz, g_yy_0_0_0_xz_yy_0_xx, g_yy_0_0_0_xz_yy_0_xy, g_yy_0_0_0_xz_yy_0_xz, g_yy_0_0_0_xz_yy_0_yy, g_yy_0_0_0_xz_yy_0_yz, g_yy_0_0_0_xz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_yy_0_xx[i] = -2.0 * g_xz_yy_0_xx[i] * a_exp + 4.0 * g_xyyz_yy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yy_0_xy[i] = -2.0 * g_xz_yy_0_xy[i] * a_exp + 4.0 * g_xyyz_yy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yy_0_xz[i] = -2.0 * g_xz_yy_0_xz[i] * a_exp + 4.0 * g_xyyz_yy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yy_0_yy[i] = -2.0 * g_xz_yy_0_yy[i] * a_exp + 4.0 * g_xyyz_yy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yy_0_yz[i] = -2.0 * g_xz_yy_0_yz[i] * a_exp + 4.0 * g_xyyz_yy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yy_0_zz[i] = -2.0 * g_xz_yy_0_zz[i] * a_exp + 4.0 * g_xyyz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (744-750)

    #pragma omp simd aligned(g_xyyz_yz_0_xx, g_xyyz_yz_0_xy, g_xyyz_yz_0_xz, g_xyyz_yz_0_yy, g_xyyz_yz_0_yz, g_xyyz_yz_0_zz, g_xz_yz_0_xx, g_xz_yz_0_xy, g_xz_yz_0_xz, g_xz_yz_0_yy, g_xz_yz_0_yz, g_xz_yz_0_zz, g_yy_0_0_0_xz_yz_0_xx, g_yy_0_0_0_xz_yz_0_xy, g_yy_0_0_0_xz_yz_0_xz, g_yy_0_0_0_xz_yz_0_yy, g_yy_0_0_0_xz_yz_0_yz, g_yy_0_0_0_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_yz_0_xx[i] = -2.0 * g_xz_yz_0_xx[i] * a_exp + 4.0 * g_xyyz_yz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yz_0_xy[i] = -2.0 * g_xz_yz_0_xy[i] * a_exp + 4.0 * g_xyyz_yz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yz_0_xz[i] = -2.0 * g_xz_yz_0_xz[i] * a_exp + 4.0 * g_xyyz_yz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yz_0_yy[i] = -2.0 * g_xz_yz_0_yy[i] * a_exp + 4.0 * g_xyyz_yz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yz_0_yz[i] = -2.0 * g_xz_yz_0_yz[i] * a_exp + 4.0 * g_xyyz_yz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_yz_0_zz[i] = -2.0 * g_xz_yz_0_zz[i] * a_exp + 4.0 * g_xyyz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (750-756)

    #pragma omp simd aligned(g_xyyz_zz_0_xx, g_xyyz_zz_0_xy, g_xyyz_zz_0_xz, g_xyyz_zz_0_yy, g_xyyz_zz_0_yz, g_xyyz_zz_0_zz, g_xz_zz_0_xx, g_xz_zz_0_xy, g_xz_zz_0_xz, g_xz_zz_0_yy, g_xz_zz_0_yz, g_xz_zz_0_zz, g_yy_0_0_0_xz_zz_0_xx, g_yy_0_0_0_xz_zz_0_xy, g_yy_0_0_0_xz_zz_0_xz, g_yy_0_0_0_xz_zz_0_yy, g_yy_0_0_0_xz_zz_0_yz, g_yy_0_0_0_xz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_xz_zz_0_xx[i] = -2.0 * g_xz_zz_0_xx[i] * a_exp + 4.0 * g_xyyz_zz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_zz_0_xy[i] = -2.0 * g_xz_zz_0_xy[i] * a_exp + 4.0 * g_xyyz_zz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_zz_0_xz[i] = -2.0 * g_xz_zz_0_xz[i] * a_exp + 4.0 * g_xyyz_zz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_zz_0_yy[i] = -2.0 * g_xz_zz_0_yy[i] * a_exp + 4.0 * g_xyyz_zz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_zz_0_yz[i] = -2.0 * g_xz_zz_0_yz[i] * a_exp + 4.0 * g_xyyz_zz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_xz_zz_0_zz[i] = -2.0 * g_xz_zz_0_zz[i] * a_exp + 4.0 * g_xyyz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (756-762)

    #pragma omp simd aligned(g_0_xx_0_xx, g_0_xx_0_xy, g_0_xx_0_xz, g_0_xx_0_yy, g_0_xx_0_yz, g_0_xx_0_zz, g_yy_0_0_0_yy_xx_0_xx, g_yy_0_0_0_yy_xx_0_xy, g_yy_0_0_0_yy_xx_0_xz, g_yy_0_0_0_yy_xx_0_yy, g_yy_0_0_0_yy_xx_0_yz, g_yy_0_0_0_yy_xx_0_zz, g_yy_xx_0_xx, g_yy_xx_0_xy, g_yy_xx_0_xz, g_yy_xx_0_yy, g_yy_xx_0_yz, g_yy_xx_0_zz, g_yyyy_xx_0_xx, g_yyyy_xx_0_xy, g_yyyy_xx_0_xz, g_yyyy_xx_0_yy, g_yyyy_xx_0_yz, g_yyyy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_xx_0_xx[i] = 2.0 * g_0_xx_0_xx[i] - 10.0 * g_yy_xx_0_xx[i] * a_exp + 4.0 * g_yyyy_xx_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xx_0_xy[i] = 2.0 * g_0_xx_0_xy[i] - 10.0 * g_yy_xx_0_xy[i] * a_exp + 4.0 * g_yyyy_xx_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xx_0_xz[i] = 2.0 * g_0_xx_0_xz[i] - 10.0 * g_yy_xx_0_xz[i] * a_exp + 4.0 * g_yyyy_xx_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xx_0_yy[i] = 2.0 * g_0_xx_0_yy[i] - 10.0 * g_yy_xx_0_yy[i] * a_exp + 4.0 * g_yyyy_xx_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xx_0_yz[i] = 2.0 * g_0_xx_0_yz[i] - 10.0 * g_yy_xx_0_yz[i] * a_exp + 4.0 * g_yyyy_xx_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xx_0_zz[i] = 2.0 * g_0_xx_0_zz[i] - 10.0 * g_yy_xx_0_zz[i] * a_exp + 4.0 * g_yyyy_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (762-768)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_yy_0_0_0_yy_xy_0_xx, g_yy_0_0_0_yy_xy_0_xy, g_yy_0_0_0_yy_xy_0_xz, g_yy_0_0_0_yy_xy_0_yy, g_yy_0_0_0_yy_xy_0_yz, g_yy_0_0_0_yy_xy_0_zz, g_yy_xy_0_xx, g_yy_xy_0_xy, g_yy_xy_0_xz, g_yy_xy_0_yy, g_yy_xy_0_yz, g_yy_xy_0_zz, g_yyyy_xy_0_xx, g_yyyy_xy_0_xy, g_yyyy_xy_0_xz, g_yyyy_xy_0_yy, g_yyyy_xy_0_yz, g_yyyy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_xy_0_xx[i] = 2.0 * g_0_xy_0_xx[i] - 10.0 * g_yy_xy_0_xx[i] * a_exp + 4.0 * g_yyyy_xy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xy_0_xy[i] = 2.0 * g_0_xy_0_xy[i] - 10.0 * g_yy_xy_0_xy[i] * a_exp + 4.0 * g_yyyy_xy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xy_0_xz[i] = 2.0 * g_0_xy_0_xz[i] - 10.0 * g_yy_xy_0_xz[i] * a_exp + 4.0 * g_yyyy_xy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xy_0_yy[i] = 2.0 * g_0_xy_0_yy[i] - 10.0 * g_yy_xy_0_yy[i] * a_exp + 4.0 * g_yyyy_xy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xy_0_yz[i] = 2.0 * g_0_xy_0_yz[i] - 10.0 * g_yy_xy_0_yz[i] * a_exp + 4.0 * g_yyyy_xy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xy_0_zz[i] = 2.0 * g_0_xy_0_zz[i] - 10.0 * g_yy_xy_0_zz[i] * a_exp + 4.0 * g_yyyy_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (768-774)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_yy_0_0_0_yy_xz_0_xx, g_yy_0_0_0_yy_xz_0_xy, g_yy_0_0_0_yy_xz_0_xz, g_yy_0_0_0_yy_xz_0_yy, g_yy_0_0_0_yy_xz_0_yz, g_yy_0_0_0_yy_xz_0_zz, g_yy_xz_0_xx, g_yy_xz_0_xy, g_yy_xz_0_xz, g_yy_xz_0_yy, g_yy_xz_0_yz, g_yy_xz_0_zz, g_yyyy_xz_0_xx, g_yyyy_xz_0_xy, g_yyyy_xz_0_xz, g_yyyy_xz_0_yy, g_yyyy_xz_0_yz, g_yyyy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_xz_0_xx[i] = 2.0 * g_0_xz_0_xx[i] - 10.0 * g_yy_xz_0_xx[i] * a_exp + 4.0 * g_yyyy_xz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xz_0_xy[i] = 2.0 * g_0_xz_0_xy[i] - 10.0 * g_yy_xz_0_xy[i] * a_exp + 4.0 * g_yyyy_xz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xz_0_xz[i] = 2.0 * g_0_xz_0_xz[i] - 10.0 * g_yy_xz_0_xz[i] * a_exp + 4.0 * g_yyyy_xz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xz_0_yy[i] = 2.0 * g_0_xz_0_yy[i] - 10.0 * g_yy_xz_0_yy[i] * a_exp + 4.0 * g_yyyy_xz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xz_0_yz[i] = 2.0 * g_0_xz_0_yz[i] - 10.0 * g_yy_xz_0_yz[i] * a_exp + 4.0 * g_yyyy_xz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_xz_0_zz[i] = 2.0 * g_0_xz_0_zz[i] - 10.0 * g_yy_xz_0_zz[i] * a_exp + 4.0 * g_yyyy_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (774-780)

    #pragma omp simd aligned(g_0_yy_0_xx, g_0_yy_0_xy, g_0_yy_0_xz, g_0_yy_0_yy, g_0_yy_0_yz, g_0_yy_0_zz, g_yy_0_0_0_yy_yy_0_xx, g_yy_0_0_0_yy_yy_0_xy, g_yy_0_0_0_yy_yy_0_xz, g_yy_0_0_0_yy_yy_0_yy, g_yy_0_0_0_yy_yy_0_yz, g_yy_0_0_0_yy_yy_0_zz, g_yy_yy_0_xx, g_yy_yy_0_xy, g_yy_yy_0_xz, g_yy_yy_0_yy, g_yy_yy_0_yz, g_yy_yy_0_zz, g_yyyy_yy_0_xx, g_yyyy_yy_0_xy, g_yyyy_yy_0_xz, g_yyyy_yy_0_yy, g_yyyy_yy_0_yz, g_yyyy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_yy_0_xx[i] = 2.0 * g_0_yy_0_xx[i] - 10.0 * g_yy_yy_0_xx[i] * a_exp + 4.0 * g_yyyy_yy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yy_0_xy[i] = 2.0 * g_0_yy_0_xy[i] - 10.0 * g_yy_yy_0_xy[i] * a_exp + 4.0 * g_yyyy_yy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yy_0_xz[i] = 2.0 * g_0_yy_0_xz[i] - 10.0 * g_yy_yy_0_xz[i] * a_exp + 4.0 * g_yyyy_yy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yy_0_yy[i] = 2.0 * g_0_yy_0_yy[i] - 10.0 * g_yy_yy_0_yy[i] * a_exp + 4.0 * g_yyyy_yy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yy_0_yz[i] = 2.0 * g_0_yy_0_yz[i] - 10.0 * g_yy_yy_0_yz[i] * a_exp + 4.0 * g_yyyy_yy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yy_0_zz[i] = 2.0 * g_0_yy_0_zz[i] - 10.0 * g_yy_yy_0_zz[i] * a_exp + 4.0 * g_yyyy_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (780-786)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_yy_0_0_0_yy_yz_0_xx, g_yy_0_0_0_yy_yz_0_xy, g_yy_0_0_0_yy_yz_0_xz, g_yy_0_0_0_yy_yz_0_yy, g_yy_0_0_0_yy_yz_0_yz, g_yy_0_0_0_yy_yz_0_zz, g_yy_yz_0_xx, g_yy_yz_0_xy, g_yy_yz_0_xz, g_yy_yz_0_yy, g_yy_yz_0_yz, g_yy_yz_0_zz, g_yyyy_yz_0_xx, g_yyyy_yz_0_xy, g_yyyy_yz_0_xz, g_yyyy_yz_0_yy, g_yyyy_yz_0_yz, g_yyyy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_yz_0_xx[i] = 2.0 * g_0_yz_0_xx[i] - 10.0 * g_yy_yz_0_xx[i] * a_exp + 4.0 * g_yyyy_yz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yz_0_xy[i] = 2.0 * g_0_yz_0_xy[i] - 10.0 * g_yy_yz_0_xy[i] * a_exp + 4.0 * g_yyyy_yz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yz_0_xz[i] = 2.0 * g_0_yz_0_xz[i] - 10.0 * g_yy_yz_0_xz[i] * a_exp + 4.0 * g_yyyy_yz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yz_0_yy[i] = 2.0 * g_0_yz_0_yy[i] - 10.0 * g_yy_yz_0_yy[i] * a_exp + 4.0 * g_yyyy_yz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yz_0_yz[i] = 2.0 * g_0_yz_0_yz[i] - 10.0 * g_yy_yz_0_yz[i] * a_exp + 4.0 * g_yyyy_yz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_yz_0_zz[i] = 2.0 * g_0_yz_0_zz[i] - 10.0 * g_yy_yz_0_zz[i] * a_exp + 4.0 * g_yyyy_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (786-792)

    #pragma omp simd aligned(g_0_zz_0_xx, g_0_zz_0_xy, g_0_zz_0_xz, g_0_zz_0_yy, g_0_zz_0_yz, g_0_zz_0_zz, g_yy_0_0_0_yy_zz_0_xx, g_yy_0_0_0_yy_zz_0_xy, g_yy_0_0_0_yy_zz_0_xz, g_yy_0_0_0_yy_zz_0_yy, g_yy_0_0_0_yy_zz_0_yz, g_yy_0_0_0_yy_zz_0_zz, g_yy_zz_0_xx, g_yy_zz_0_xy, g_yy_zz_0_xz, g_yy_zz_0_yy, g_yy_zz_0_yz, g_yy_zz_0_zz, g_yyyy_zz_0_xx, g_yyyy_zz_0_xy, g_yyyy_zz_0_xz, g_yyyy_zz_0_yy, g_yyyy_zz_0_yz, g_yyyy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yy_zz_0_xx[i] = 2.0 * g_0_zz_0_xx[i] - 10.0 * g_yy_zz_0_xx[i] * a_exp + 4.0 * g_yyyy_zz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_zz_0_xy[i] = 2.0 * g_0_zz_0_xy[i] - 10.0 * g_yy_zz_0_xy[i] * a_exp + 4.0 * g_yyyy_zz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_zz_0_xz[i] = 2.0 * g_0_zz_0_xz[i] - 10.0 * g_yy_zz_0_xz[i] * a_exp + 4.0 * g_yyyy_zz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_zz_0_yy[i] = 2.0 * g_0_zz_0_yy[i] - 10.0 * g_yy_zz_0_yy[i] * a_exp + 4.0 * g_yyyy_zz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_zz_0_yz[i] = 2.0 * g_0_zz_0_yz[i] - 10.0 * g_yy_zz_0_yz[i] * a_exp + 4.0 * g_yyyy_zz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_yy_zz_0_zz[i] = 2.0 * g_0_zz_0_zz[i] - 10.0 * g_yy_zz_0_zz[i] * a_exp + 4.0 * g_yyyy_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (792-798)

    #pragma omp simd aligned(g_yy_0_0_0_yz_xx_0_xx, g_yy_0_0_0_yz_xx_0_xy, g_yy_0_0_0_yz_xx_0_xz, g_yy_0_0_0_yz_xx_0_yy, g_yy_0_0_0_yz_xx_0_yz, g_yy_0_0_0_yz_xx_0_zz, g_yyyz_xx_0_xx, g_yyyz_xx_0_xy, g_yyyz_xx_0_xz, g_yyyz_xx_0_yy, g_yyyz_xx_0_yz, g_yyyz_xx_0_zz, g_yz_xx_0_xx, g_yz_xx_0_xy, g_yz_xx_0_xz, g_yz_xx_0_yy, g_yz_xx_0_yz, g_yz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_xx_0_xx[i] = -6.0 * g_yz_xx_0_xx[i] * a_exp + 4.0 * g_yyyz_xx_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xx_0_xy[i] = -6.0 * g_yz_xx_0_xy[i] * a_exp + 4.0 * g_yyyz_xx_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xx_0_xz[i] = -6.0 * g_yz_xx_0_xz[i] * a_exp + 4.0 * g_yyyz_xx_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xx_0_yy[i] = -6.0 * g_yz_xx_0_yy[i] * a_exp + 4.0 * g_yyyz_xx_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xx_0_yz[i] = -6.0 * g_yz_xx_0_yz[i] * a_exp + 4.0 * g_yyyz_xx_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xx_0_zz[i] = -6.0 * g_yz_xx_0_zz[i] * a_exp + 4.0 * g_yyyz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (798-804)

    #pragma omp simd aligned(g_yy_0_0_0_yz_xy_0_xx, g_yy_0_0_0_yz_xy_0_xy, g_yy_0_0_0_yz_xy_0_xz, g_yy_0_0_0_yz_xy_0_yy, g_yy_0_0_0_yz_xy_0_yz, g_yy_0_0_0_yz_xy_0_zz, g_yyyz_xy_0_xx, g_yyyz_xy_0_xy, g_yyyz_xy_0_xz, g_yyyz_xy_0_yy, g_yyyz_xy_0_yz, g_yyyz_xy_0_zz, g_yz_xy_0_xx, g_yz_xy_0_xy, g_yz_xy_0_xz, g_yz_xy_0_yy, g_yz_xy_0_yz, g_yz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_xy_0_xx[i] = -6.0 * g_yz_xy_0_xx[i] * a_exp + 4.0 * g_yyyz_xy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xy_0_xy[i] = -6.0 * g_yz_xy_0_xy[i] * a_exp + 4.0 * g_yyyz_xy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xy_0_xz[i] = -6.0 * g_yz_xy_0_xz[i] * a_exp + 4.0 * g_yyyz_xy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xy_0_yy[i] = -6.0 * g_yz_xy_0_yy[i] * a_exp + 4.0 * g_yyyz_xy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xy_0_yz[i] = -6.0 * g_yz_xy_0_yz[i] * a_exp + 4.0 * g_yyyz_xy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xy_0_zz[i] = -6.0 * g_yz_xy_0_zz[i] * a_exp + 4.0 * g_yyyz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (804-810)

    #pragma omp simd aligned(g_yy_0_0_0_yz_xz_0_xx, g_yy_0_0_0_yz_xz_0_xy, g_yy_0_0_0_yz_xz_0_xz, g_yy_0_0_0_yz_xz_0_yy, g_yy_0_0_0_yz_xz_0_yz, g_yy_0_0_0_yz_xz_0_zz, g_yyyz_xz_0_xx, g_yyyz_xz_0_xy, g_yyyz_xz_0_xz, g_yyyz_xz_0_yy, g_yyyz_xz_0_yz, g_yyyz_xz_0_zz, g_yz_xz_0_xx, g_yz_xz_0_xy, g_yz_xz_0_xz, g_yz_xz_0_yy, g_yz_xz_0_yz, g_yz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_xz_0_xx[i] = -6.0 * g_yz_xz_0_xx[i] * a_exp + 4.0 * g_yyyz_xz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xz_0_xy[i] = -6.0 * g_yz_xz_0_xy[i] * a_exp + 4.0 * g_yyyz_xz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xz_0_xz[i] = -6.0 * g_yz_xz_0_xz[i] * a_exp + 4.0 * g_yyyz_xz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xz_0_yy[i] = -6.0 * g_yz_xz_0_yy[i] * a_exp + 4.0 * g_yyyz_xz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xz_0_yz[i] = -6.0 * g_yz_xz_0_yz[i] * a_exp + 4.0 * g_yyyz_xz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_xz_0_zz[i] = -6.0 * g_yz_xz_0_zz[i] * a_exp + 4.0 * g_yyyz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (810-816)

    #pragma omp simd aligned(g_yy_0_0_0_yz_yy_0_xx, g_yy_0_0_0_yz_yy_0_xy, g_yy_0_0_0_yz_yy_0_xz, g_yy_0_0_0_yz_yy_0_yy, g_yy_0_0_0_yz_yy_0_yz, g_yy_0_0_0_yz_yy_0_zz, g_yyyz_yy_0_xx, g_yyyz_yy_0_xy, g_yyyz_yy_0_xz, g_yyyz_yy_0_yy, g_yyyz_yy_0_yz, g_yyyz_yy_0_zz, g_yz_yy_0_xx, g_yz_yy_0_xy, g_yz_yy_0_xz, g_yz_yy_0_yy, g_yz_yy_0_yz, g_yz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_yy_0_xx[i] = -6.0 * g_yz_yy_0_xx[i] * a_exp + 4.0 * g_yyyz_yy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yy_0_xy[i] = -6.0 * g_yz_yy_0_xy[i] * a_exp + 4.0 * g_yyyz_yy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yy_0_xz[i] = -6.0 * g_yz_yy_0_xz[i] * a_exp + 4.0 * g_yyyz_yy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yy_0_yy[i] = -6.0 * g_yz_yy_0_yy[i] * a_exp + 4.0 * g_yyyz_yy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yy_0_yz[i] = -6.0 * g_yz_yy_0_yz[i] * a_exp + 4.0 * g_yyyz_yy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yy_0_zz[i] = -6.0 * g_yz_yy_0_zz[i] * a_exp + 4.0 * g_yyyz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (816-822)

    #pragma omp simd aligned(g_yy_0_0_0_yz_yz_0_xx, g_yy_0_0_0_yz_yz_0_xy, g_yy_0_0_0_yz_yz_0_xz, g_yy_0_0_0_yz_yz_0_yy, g_yy_0_0_0_yz_yz_0_yz, g_yy_0_0_0_yz_yz_0_zz, g_yyyz_yz_0_xx, g_yyyz_yz_0_xy, g_yyyz_yz_0_xz, g_yyyz_yz_0_yy, g_yyyz_yz_0_yz, g_yyyz_yz_0_zz, g_yz_yz_0_xx, g_yz_yz_0_xy, g_yz_yz_0_xz, g_yz_yz_0_yy, g_yz_yz_0_yz, g_yz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_yz_0_xx[i] = -6.0 * g_yz_yz_0_xx[i] * a_exp + 4.0 * g_yyyz_yz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yz_0_xy[i] = -6.0 * g_yz_yz_0_xy[i] * a_exp + 4.0 * g_yyyz_yz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yz_0_xz[i] = -6.0 * g_yz_yz_0_xz[i] * a_exp + 4.0 * g_yyyz_yz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yz_0_yy[i] = -6.0 * g_yz_yz_0_yy[i] * a_exp + 4.0 * g_yyyz_yz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yz_0_yz[i] = -6.0 * g_yz_yz_0_yz[i] * a_exp + 4.0 * g_yyyz_yz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_yz_0_zz[i] = -6.0 * g_yz_yz_0_zz[i] * a_exp + 4.0 * g_yyyz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (822-828)

    #pragma omp simd aligned(g_yy_0_0_0_yz_zz_0_xx, g_yy_0_0_0_yz_zz_0_xy, g_yy_0_0_0_yz_zz_0_xz, g_yy_0_0_0_yz_zz_0_yy, g_yy_0_0_0_yz_zz_0_yz, g_yy_0_0_0_yz_zz_0_zz, g_yyyz_zz_0_xx, g_yyyz_zz_0_xy, g_yyyz_zz_0_xz, g_yyyz_zz_0_yy, g_yyyz_zz_0_yz, g_yyyz_zz_0_zz, g_yz_zz_0_xx, g_yz_zz_0_xy, g_yz_zz_0_xz, g_yz_zz_0_yy, g_yz_zz_0_yz, g_yz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_yz_zz_0_xx[i] = -6.0 * g_yz_zz_0_xx[i] * a_exp + 4.0 * g_yyyz_zz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_zz_0_xy[i] = -6.0 * g_yz_zz_0_xy[i] * a_exp + 4.0 * g_yyyz_zz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_zz_0_xz[i] = -6.0 * g_yz_zz_0_xz[i] * a_exp + 4.0 * g_yyyz_zz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_zz_0_yy[i] = -6.0 * g_yz_zz_0_yy[i] * a_exp + 4.0 * g_yyyz_zz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_zz_0_yz[i] = -6.0 * g_yz_zz_0_yz[i] * a_exp + 4.0 * g_yyyz_zz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_yz_zz_0_zz[i] = -6.0 * g_yz_zz_0_zz[i] * a_exp + 4.0 * g_yyyz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (828-834)

    #pragma omp simd aligned(g_yy_0_0_0_zz_xx_0_xx, g_yy_0_0_0_zz_xx_0_xy, g_yy_0_0_0_zz_xx_0_xz, g_yy_0_0_0_zz_xx_0_yy, g_yy_0_0_0_zz_xx_0_yz, g_yy_0_0_0_zz_xx_0_zz, g_yyzz_xx_0_xx, g_yyzz_xx_0_xy, g_yyzz_xx_0_xz, g_yyzz_xx_0_yy, g_yyzz_xx_0_yz, g_yyzz_xx_0_zz, g_zz_xx_0_xx, g_zz_xx_0_xy, g_zz_xx_0_xz, g_zz_xx_0_yy, g_zz_xx_0_yz, g_zz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_xx_0_xx[i] = -2.0 * g_zz_xx_0_xx[i] * a_exp + 4.0 * g_yyzz_xx_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xx_0_xy[i] = -2.0 * g_zz_xx_0_xy[i] * a_exp + 4.0 * g_yyzz_xx_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xx_0_xz[i] = -2.0 * g_zz_xx_0_xz[i] * a_exp + 4.0 * g_yyzz_xx_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xx_0_yy[i] = -2.0 * g_zz_xx_0_yy[i] * a_exp + 4.0 * g_yyzz_xx_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xx_0_yz[i] = -2.0 * g_zz_xx_0_yz[i] * a_exp + 4.0 * g_yyzz_xx_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xx_0_zz[i] = -2.0 * g_zz_xx_0_zz[i] * a_exp + 4.0 * g_yyzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (834-840)

    #pragma omp simd aligned(g_yy_0_0_0_zz_xy_0_xx, g_yy_0_0_0_zz_xy_0_xy, g_yy_0_0_0_zz_xy_0_xz, g_yy_0_0_0_zz_xy_0_yy, g_yy_0_0_0_zz_xy_0_yz, g_yy_0_0_0_zz_xy_0_zz, g_yyzz_xy_0_xx, g_yyzz_xy_0_xy, g_yyzz_xy_0_xz, g_yyzz_xy_0_yy, g_yyzz_xy_0_yz, g_yyzz_xy_0_zz, g_zz_xy_0_xx, g_zz_xy_0_xy, g_zz_xy_0_xz, g_zz_xy_0_yy, g_zz_xy_0_yz, g_zz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_xy_0_xx[i] = -2.0 * g_zz_xy_0_xx[i] * a_exp + 4.0 * g_yyzz_xy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xy_0_xy[i] = -2.0 * g_zz_xy_0_xy[i] * a_exp + 4.0 * g_yyzz_xy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xy_0_xz[i] = -2.0 * g_zz_xy_0_xz[i] * a_exp + 4.0 * g_yyzz_xy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xy_0_yy[i] = -2.0 * g_zz_xy_0_yy[i] * a_exp + 4.0 * g_yyzz_xy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xy_0_yz[i] = -2.0 * g_zz_xy_0_yz[i] * a_exp + 4.0 * g_yyzz_xy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xy_0_zz[i] = -2.0 * g_zz_xy_0_zz[i] * a_exp + 4.0 * g_yyzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (840-846)

    #pragma omp simd aligned(g_yy_0_0_0_zz_xz_0_xx, g_yy_0_0_0_zz_xz_0_xy, g_yy_0_0_0_zz_xz_0_xz, g_yy_0_0_0_zz_xz_0_yy, g_yy_0_0_0_zz_xz_0_yz, g_yy_0_0_0_zz_xz_0_zz, g_yyzz_xz_0_xx, g_yyzz_xz_0_xy, g_yyzz_xz_0_xz, g_yyzz_xz_0_yy, g_yyzz_xz_0_yz, g_yyzz_xz_0_zz, g_zz_xz_0_xx, g_zz_xz_0_xy, g_zz_xz_0_xz, g_zz_xz_0_yy, g_zz_xz_0_yz, g_zz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_xz_0_xx[i] = -2.0 * g_zz_xz_0_xx[i] * a_exp + 4.0 * g_yyzz_xz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xz_0_xy[i] = -2.0 * g_zz_xz_0_xy[i] * a_exp + 4.0 * g_yyzz_xz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xz_0_xz[i] = -2.0 * g_zz_xz_0_xz[i] * a_exp + 4.0 * g_yyzz_xz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xz_0_yy[i] = -2.0 * g_zz_xz_0_yy[i] * a_exp + 4.0 * g_yyzz_xz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xz_0_yz[i] = -2.0 * g_zz_xz_0_yz[i] * a_exp + 4.0 * g_yyzz_xz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_xz_0_zz[i] = -2.0 * g_zz_xz_0_zz[i] * a_exp + 4.0 * g_yyzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (846-852)

    #pragma omp simd aligned(g_yy_0_0_0_zz_yy_0_xx, g_yy_0_0_0_zz_yy_0_xy, g_yy_0_0_0_zz_yy_0_xz, g_yy_0_0_0_zz_yy_0_yy, g_yy_0_0_0_zz_yy_0_yz, g_yy_0_0_0_zz_yy_0_zz, g_yyzz_yy_0_xx, g_yyzz_yy_0_xy, g_yyzz_yy_0_xz, g_yyzz_yy_0_yy, g_yyzz_yy_0_yz, g_yyzz_yy_0_zz, g_zz_yy_0_xx, g_zz_yy_0_xy, g_zz_yy_0_xz, g_zz_yy_0_yy, g_zz_yy_0_yz, g_zz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_yy_0_xx[i] = -2.0 * g_zz_yy_0_xx[i] * a_exp + 4.0 * g_yyzz_yy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yy_0_xy[i] = -2.0 * g_zz_yy_0_xy[i] * a_exp + 4.0 * g_yyzz_yy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yy_0_xz[i] = -2.0 * g_zz_yy_0_xz[i] * a_exp + 4.0 * g_yyzz_yy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yy_0_yy[i] = -2.0 * g_zz_yy_0_yy[i] * a_exp + 4.0 * g_yyzz_yy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yy_0_yz[i] = -2.0 * g_zz_yy_0_yz[i] * a_exp + 4.0 * g_yyzz_yy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yy_0_zz[i] = -2.0 * g_zz_yy_0_zz[i] * a_exp + 4.0 * g_yyzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (852-858)

    #pragma omp simd aligned(g_yy_0_0_0_zz_yz_0_xx, g_yy_0_0_0_zz_yz_0_xy, g_yy_0_0_0_zz_yz_0_xz, g_yy_0_0_0_zz_yz_0_yy, g_yy_0_0_0_zz_yz_0_yz, g_yy_0_0_0_zz_yz_0_zz, g_yyzz_yz_0_xx, g_yyzz_yz_0_xy, g_yyzz_yz_0_xz, g_yyzz_yz_0_yy, g_yyzz_yz_0_yz, g_yyzz_yz_0_zz, g_zz_yz_0_xx, g_zz_yz_0_xy, g_zz_yz_0_xz, g_zz_yz_0_yy, g_zz_yz_0_yz, g_zz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_yz_0_xx[i] = -2.0 * g_zz_yz_0_xx[i] * a_exp + 4.0 * g_yyzz_yz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yz_0_xy[i] = -2.0 * g_zz_yz_0_xy[i] * a_exp + 4.0 * g_yyzz_yz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yz_0_xz[i] = -2.0 * g_zz_yz_0_xz[i] * a_exp + 4.0 * g_yyzz_yz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yz_0_yy[i] = -2.0 * g_zz_yz_0_yy[i] * a_exp + 4.0 * g_yyzz_yz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yz_0_yz[i] = -2.0 * g_zz_yz_0_yz[i] * a_exp + 4.0 * g_yyzz_yz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_yz_0_zz[i] = -2.0 * g_zz_yz_0_zz[i] * a_exp + 4.0 * g_yyzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (858-864)

    #pragma omp simd aligned(g_yy_0_0_0_zz_zz_0_xx, g_yy_0_0_0_zz_zz_0_xy, g_yy_0_0_0_zz_zz_0_xz, g_yy_0_0_0_zz_zz_0_yy, g_yy_0_0_0_zz_zz_0_yz, g_yy_0_0_0_zz_zz_0_zz, g_yyzz_zz_0_xx, g_yyzz_zz_0_xy, g_yyzz_zz_0_xz, g_yyzz_zz_0_yy, g_yyzz_zz_0_yz, g_yyzz_zz_0_zz, g_zz_zz_0_xx, g_zz_zz_0_xy, g_zz_zz_0_xz, g_zz_zz_0_yy, g_zz_zz_0_yz, g_zz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_zz_zz_0_xx[i] = -2.0 * g_zz_zz_0_xx[i] * a_exp + 4.0 * g_yyzz_zz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_zz_0_xy[i] = -2.0 * g_zz_zz_0_xy[i] * a_exp + 4.0 * g_yyzz_zz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_zz_0_xz[i] = -2.0 * g_zz_zz_0_xz[i] * a_exp + 4.0 * g_yyzz_zz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_zz_0_yy[i] = -2.0 * g_zz_zz_0_yy[i] * a_exp + 4.0 * g_yyzz_zz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_zz_0_yz[i] = -2.0 * g_zz_zz_0_yz[i] * a_exp + 4.0 * g_yyzz_zz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_zz_zz_0_zz[i] = -2.0 * g_zz_zz_0_zz[i] * a_exp + 4.0 * g_yyzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (864-870)

    #pragma omp simd aligned(g_xxyz_xx_0_xx, g_xxyz_xx_0_xy, g_xxyz_xx_0_xz, g_xxyz_xx_0_yy, g_xxyz_xx_0_yz, g_xxyz_xx_0_zz, g_yz_0_0_0_xx_xx_0_xx, g_yz_0_0_0_xx_xx_0_xy, g_yz_0_0_0_xx_xx_0_xz, g_yz_0_0_0_xx_xx_0_yy, g_yz_0_0_0_xx_xx_0_yz, g_yz_0_0_0_xx_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_xx_0_xx[i] = 4.0 * g_xxyz_xx_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xx_0_xy[i] = 4.0 * g_xxyz_xx_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xx_0_xz[i] = 4.0 * g_xxyz_xx_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xx_0_yy[i] = 4.0 * g_xxyz_xx_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xx_0_yz[i] = 4.0 * g_xxyz_xx_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xx_0_zz[i] = 4.0 * g_xxyz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (870-876)

    #pragma omp simd aligned(g_xxyz_xy_0_xx, g_xxyz_xy_0_xy, g_xxyz_xy_0_xz, g_xxyz_xy_0_yy, g_xxyz_xy_0_yz, g_xxyz_xy_0_zz, g_yz_0_0_0_xx_xy_0_xx, g_yz_0_0_0_xx_xy_0_xy, g_yz_0_0_0_xx_xy_0_xz, g_yz_0_0_0_xx_xy_0_yy, g_yz_0_0_0_xx_xy_0_yz, g_yz_0_0_0_xx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_xy_0_xx[i] = 4.0 * g_xxyz_xy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xy_0_xy[i] = 4.0 * g_xxyz_xy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xy_0_xz[i] = 4.0 * g_xxyz_xy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xy_0_yy[i] = 4.0 * g_xxyz_xy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xy_0_yz[i] = 4.0 * g_xxyz_xy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xy_0_zz[i] = 4.0 * g_xxyz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (876-882)

    #pragma omp simd aligned(g_xxyz_xz_0_xx, g_xxyz_xz_0_xy, g_xxyz_xz_0_xz, g_xxyz_xz_0_yy, g_xxyz_xz_0_yz, g_xxyz_xz_0_zz, g_yz_0_0_0_xx_xz_0_xx, g_yz_0_0_0_xx_xz_0_xy, g_yz_0_0_0_xx_xz_0_xz, g_yz_0_0_0_xx_xz_0_yy, g_yz_0_0_0_xx_xz_0_yz, g_yz_0_0_0_xx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_xz_0_xx[i] = 4.0 * g_xxyz_xz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xz_0_xy[i] = 4.0 * g_xxyz_xz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xz_0_xz[i] = 4.0 * g_xxyz_xz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xz_0_yy[i] = 4.0 * g_xxyz_xz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xz_0_yz[i] = 4.0 * g_xxyz_xz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_xz_0_zz[i] = 4.0 * g_xxyz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (882-888)

    #pragma omp simd aligned(g_xxyz_yy_0_xx, g_xxyz_yy_0_xy, g_xxyz_yy_0_xz, g_xxyz_yy_0_yy, g_xxyz_yy_0_yz, g_xxyz_yy_0_zz, g_yz_0_0_0_xx_yy_0_xx, g_yz_0_0_0_xx_yy_0_xy, g_yz_0_0_0_xx_yy_0_xz, g_yz_0_0_0_xx_yy_0_yy, g_yz_0_0_0_xx_yy_0_yz, g_yz_0_0_0_xx_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_yy_0_xx[i] = 4.0 * g_xxyz_yy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yy_0_xy[i] = 4.0 * g_xxyz_yy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yy_0_xz[i] = 4.0 * g_xxyz_yy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yy_0_yy[i] = 4.0 * g_xxyz_yy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yy_0_yz[i] = 4.0 * g_xxyz_yy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yy_0_zz[i] = 4.0 * g_xxyz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (888-894)

    #pragma omp simd aligned(g_xxyz_yz_0_xx, g_xxyz_yz_0_xy, g_xxyz_yz_0_xz, g_xxyz_yz_0_yy, g_xxyz_yz_0_yz, g_xxyz_yz_0_zz, g_yz_0_0_0_xx_yz_0_xx, g_yz_0_0_0_xx_yz_0_xy, g_yz_0_0_0_xx_yz_0_xz, g_yz_0_0_0_xx_yz_0_yy, g_yz_0_0_0_xx_yz_0_yz, g_yz_0_0_0_xx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_yz_0_xx[i] = 4.0 * g_xxyz_yz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yz_0_xy[i] = 4.0 * g_xxyz_yz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yz_0_xz[i] = 4.0 * g_xxyz_yz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yz_0_yy[i] = 4.0 * g_xxyz_yz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yz_0_yz[i] = 4.0 * g_xxyz_yz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_yz_0_zz[i] = 4.0 * g_xxyz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (894-900)

    #pragma omp simd aligned(g_xxyz_zz_0_xx, g_xxyz_zz_0_xy, g_xxyz_zz_0_xz, g_xxyz_zz_0_yy, g_xxyz_zz_0_yz, g_xxyz_zz_0_zz, g_yz_0_0_0_xx_zz_0_xx, g_yz_0_0_0_xx_zz_0_xy, g_yz_0_0_0_xx_zz_0_xz, g_yz_0_0_0_xx_zz_0_yy, g_yz_0_0_0_xx_zz_0_yz, g_yz_0_0_0_xx_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xx_zz_0_xx[i] = 4.0 * g_xxyz_zz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_zz_0_xy[i] = 4.0 * g_xxyz_zz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_zz_0_xz[i] = 4.0 * g_xxyz_zz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_zz_0_yy[i] = 4.0 * g_xxyz_zz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_zz_0_yz[i] = 4.0 * g_xxyz_zz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xx_zz_0_zz[i] = 4.0 * g_xxyz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (900-906)

    #pragma omp simd aligned(g_xyyz_xx_0_xx, g_xyyz_xx_0_xy, g_xyyz_xx_0_xz, g_xyyz_xx_0_yy, g_xyyz_xx_0_yz, g_xyyz_xx_0_zz, g_xz_xx_0_xx, g_xz_xx_0_xy, g_xz_xx_0_xz, g_xz_xx_0_yy, g_xz_xx_0_yz, g_xz_xx_0_zz, g_yz_0_0_0_xy_xx_0_xx, g_yz_0_0_0_xy_xx_0_xy, g_yz_0_0_0_xy_xx_0_xz, g_yz_0_0_0_xy_xx_0_yy, g_yz_0_0_0_xy_xx_0_yz, g_yz_0_0_0_xy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_xx_0_xx[i] = -2.0 * g_xz_xx_0_xx[i] * a_exp + 4.0 * g_xyyz_xx_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xx_0_xy[i] = -2.0 * g_xz_xx_0_xy[i] * a_exp + 4.0 * g_xyyz_xx_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xx_0_xz[i] = -2.0 * g_xz_xx_0_xz[i] * a_exp + 4.0 * g_xyyz_xx_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xx_0_yy[i] = -2.0 * g_xz_xx_0_yy[i] * a_exp + 4.0 * g_xyyz_xx_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xx_0_yz[i] = -2.0 * g_xz_xx_0_yz[i] * a_exp + 4.0 * g_xyyz_xx_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xx_0_zz[i] = -2.0 * g_xz_xx_0_zz[i] * a_exp + 4.0 * g_xyyz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (906-912)

    #pragma omp simd aligned(g_xyyz_xy_0_xx, g_xyyz_xy_0_xy, g_xyyz_xy_0_xz, g_xyyz_xy_0_yy, g_xyyz_xy_0_yz, g_xyyz_xy_0_zz, g_xz_xy_0_xx, g_xz_xy_0_xy, g_xz_xy_0_xz, g_xz_xy_0_yy, g_xz_xy_0_yz, g_xz_xy_0_zz, g_yz_0_0_0_xy_xy_0_xx, g_yz_0_0_0_xy_xy_0_xy, g_yz_0_0_0_xy_xy_0_xz, g_yz_0_0_0_xy_xy_0_yy, g_yz_0_0_0_xy_xy_0_yz, g_yz_0_0_0_xy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_xy_0_xx[i] = -2.0 * g_xz_xy_0_xx[i] * a_exp + 4.0 * g_xyyz_xy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xy_0_xy[i] = -2.0 * g_xz_xy_0_xy[i] * a_exp + 4.0 * g_xyyz_xy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xy_0_xz[i] = -2.0 * g_xz_xy_0_xz[i] * a_exp + 4.0 * g_xyyz_xy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xy_0_yy[i] = -2.0 * g_xz_xy_0_yy[i] * a_exp + 4.0 * g_xyyz_xy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xy_0_yz[i] = -2.0 * g_xz_xy_0_yz[i] * a_exp + 4.0 * g_xyyz_xy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xy_0_zz[i] = -2.0 * g_xz_xy_0_zz[i] * a_exp + 4.0 * g_xyyz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (912-918)

    #pragma omp simd aligned(g_xyyz_xz_0_xx, g_xyyz_xz_0_xy, g_xyyz_xz_0_xz, g_xyyz_xz_0_yy, g_xyyz_xz_0_yz, g_xyyz_xz_0_zz, g_xz_xz_0_xx, g_xz_xz_0_xy, g_xz_xz_0_xz, g_xz_xz_0_yy, g_xz_xz_0_yz, g_xz_xz_0_zz, g_yz_0_0_0_xy_xz_0_xx, g_yz_0_0_0_xy_xz_0_xy, g_yz_0_0_0_xy_xz_0_xz, g_yz_0_0_0_xy_xz_0_yy, g_yz_0_0_0_xy_xz_0_yz, g_yz_0_0_0_xy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_xz_0_xx[i] = -2.0 * g_xz_xz_0_xx[i] * a_exp + 4.0 * g_xyyz_xz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xz_0_xy[i] = -2.0 * g_xz_xz_0_xy[i] * a_exp + 4.0 * g_xyyz_xz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xz_0_xz[i] = -2.0 * g_xz_xz_0_xz[i] * a_exp + 4.0 * g_xyyz_xz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xz_0_yy[i] = -2.0 * g_xz_xz_0_yy[i] * a_exp + 4.0 * g_xyyz_xz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xz_0_yz[i] = -2.0 * g_xz_xz_0_yz[i] * a_exp + 4.0 * g_xyyz_xz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_xz_0_zz[i] = -2.0 * g_xz_xz_0_zz[i] * a_exp + 4.0 * g_xyyz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (918-924)

    #pragma omp simd aligned(g_xyyz_yy_0_xx, g_xyyz_yy_0_xy, g_xyyz_yy_0_xz, g_xyyz_yy_0_yy, g_xyyz_yy_0_yz, g_xyyz_yy_0_zz, g_xz_yy_0_xx, g_xz_yy_0_xy, g_xz_yy_0_xz, g_xz_yy_0_yy, g_xz_yy_0_yz, g_xz_yy_0_zz, g_yz_0_0_0_xy_yy_0_xx, g_yz_0_0_0_xy_yy_0_xy, g_yz_0_0_0_xy_yy_0_xz, g_yz_0_0_0_xy_yy_0_yy, g_yz_0_0_0_xy_yy_0_yz, g_yz_0_0_0_xy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_yy_0_xx[i] = -2.0 * g_xz_yy_0_xx[i] * a_exp + 4.0 * g_xyyz_yy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yy_0_xy[i] = -2.0 * g_xz_yy_0_xy[i] * a_exp + 4.0 * g_xyyz_yy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yy_0_xz[i] = -2.0 * g_xz_yy_0_xz[i] * a_exp + 4.0 * g_xyyz_yy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yy_0_yy[i] = -2.0 * g_xz_yy_0_yy[i] * a_exp + 4.0 * g_xyyz_yy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yy_0_yz[i] = -2.0 * g_xz_yy_0_yz[i] * a_exp + 4.0 * g_xyyz_yy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yy_0_zz[i] = -2.0 * g_xz_yy_0_zz[i] * a_exp + 4.0 * g_xyyz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (924-930)

    #pragma omp simd aligned(g_xyyz_yz_0_xx, g_xyyz_yz_0_xy, g_xyyz_yz_0_xz, g_xyyz_yz_0_yy, g_xyyz_yz_0_yz, g_xyyz_yz_0_zz, g_xz_yz_0_xx, g_xz_yz_0_xy, g_xz_yz_0_xz, g_xz_yz_0_yy, g_xz_yz_0_yz, g_xz_yz_0_zz, g_yz_0_0_0_xy_yz_0_xx, g_yz_0_0_0_xy_yz_0_xy, g_yz_0_0_0_xy_yz_0_xz, g_yz_0_0_0_xy_yz_0_yy, g_yz_0_0_0_xy_yz_0_yz, g_yz_0_0_0_xy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_yz_0_xx[i] = -2.0 * g_xz_yz_0_xx[i] * a_exp + 4.0 * g_xyyz_yz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yz_0_xy[i] = -2.0 * g_xz_yz_0_xy[i] * a_exp + 4.0 * g_xyyz_yz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yz_0_xz[i] = -2.0 * g_xz_yz_0_xz[i] * a_exp + 4.0 * g_xyyz_yz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yz_0_yy[i] = -2.0 * g_xz_yz_0_yy[i] * a_exp + 4.0 * g_xyyz_yz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yz_0_yz[i] = -2.0 * g_xz_yz_0_yz[i] * a_exp + 4.0 * g_xyyz_yz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_yz_0_zz[i] = -2.0 * g_xz_yz_0_zz[i] * a_exp + 4.0 * g_xyyz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (930-936)

    #pragma omp simd aligned(g_xyyz_zz_0_xx, g_xyyz_zz_0_xy, g_xyyz_zz_0_xz, g_xyyz_zz_0_yy, g_xyyz_zz_0_yz, g_xyyz_zz_0_zz, g_xz_zz_0_xx, g_xz_zz_0_xy, g_xz_zz_0_xz, g_xz_zz_0_yy, g_xz_zz_0_yz, g_xz_zz_0_zz, g_yz_0_0_0_xy_zz_0_xx, g_yz_0_0_0_xy_zz_0_xy, g_yz_0_0_0_xy_zz_0_xz, g_yz_0_0_0_xy_zz_0_yy, g_yz_0_0_0_xy_zz_0_yz, g_yz_0_0_0_xy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xy_zz_0_xx[i] = -2.0 * g_xz_zz_0_xx[i] * a_exp + 4.0 * g_xyyz_zz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_zz_0_xy[i] = -2.0 * g_xz_zz_0_xy[i] * a_exp + 4.0 * g_xyyz_zz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_zz_0_xz[i] = -2.0 * g_xz_zz_0_xz[i] * a_exp + 4.0 * g_xyyz_zz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_zz_0_yy[i] = -2.0 * g_xz_zz_0_yy[i] * a_exp + 4.0 * g_xyyz_zz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_zz_0_yz[i] = -2.0 * g_xz_zz_0_yz[i] * a_exp + 4.0 * g_xyyz_zz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xy_zz_0_zz[i] = -2.0 * g_xz_zz_0_zz[i] * a_exp + 4.0 * g_xyyz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (936-942)

    #pragma omp simd aligned(g_xy_xx_0_xx, g_xy_xx_0_xy, g_xy_xx_0_xz, g_xy_xx_0_yy, g_xy_xx_0_yz, g_xy_xx_0_zz, g_xyzz_xx_0_xx, g_xyzz_xx_0_xy, g_xyzz_xx_0_xz, g_xyzz_xx_0_yy, g_xyzz_xx_0_yz, g_xyzz_xx_0_zz, g_yz_0_0_0_xz_xx_0_xx, g_yz_0_0_0_xz_xx_0_xy, g_yz_0_0_0_xz_xx_0_xz, g_yz_0_0_0_xz_xx_0_yy, g_yz_0_0_0_xz_xx_0_yz, g_yz_0_0_0_xz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_xx_0_xx[i] = -2.0 * g_xy_xx_0_xx[i] * a_exp + 4.0 * g_xyzz_xx_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xx_0_xy[i] = -2.0 * g_xy_xx_0_xy[i] * a_exp + 4.0 * g_xyzz_xx_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xx_0_xz[i] = -2.0 * g_xy_xx_0_xz[i] * a_exp + 4.0 * g_xyzz_xx_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xx_0_yy[i] = -2.0 * g_xy_xx_0_yy[i] * a_exp + 4.0 * g_xyzz_xx_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xx_0_yz[i] = -2.0 * g_xy_xx_0_yz[i] * a_exp + 4.0 * g_xyzz_xx_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xx_0_zz[i] = -2.0 * g_xy_xx_0_zz[i] * a_exp + 4.0 * g_xyzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (942-948)

    #pragma omp simd aligned(g_xy_xy_0_xx, g_xy_xy_0_xy, g_xy_xy_0_xz, g_xy_xy_0_yy, g_xy_xy_0_yz, g_xy_xy_0_zz, g_xyzz_xy_0_xx, g_xyzz_xy_0_xy, g_xyzz_xy_0_xz, g_xyzz_xy_0_yy, g_xyzz_xy_0_yz, g_xyzz_xy_0_zz, g_yz_0_0_0_xz_xy_0_xx, g_yz_0_0_0_xz_xy_0_xy, g_yz_0_0_0_xz_xy_0_xz, g_yz_0_0_0_xz_xy_0_yy, g_yz_0_0_0_xz_xy_0_yz, g_yz_0_0_0_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_xy_0_xx[i] = -2.0 * g_xy_xy_0_xx[i] * a_exp + 4.0 * g_xyzz_xy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xy_0_xy[i] = -2.0 * g_xy_xy_0_xy[i] * a_exp + 4.0 * g_xyzz_xy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xy_0_xz[i] = -2.0 * g_xy_xy_0_xz[i] * a_exp + 4.0 * g_xyzz_xy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xy_0_yy[i] = -2.0 * g_xy_xy_0_yy[i] * a_exp + 4.0 * g_xyzz_xy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xy_0_yz[i] = -2.0 * g_xy_xy_0_yz[i] * a_exp + 4.0 * g_xyzz_xy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xy_0_zz[i] = -2.0 * g_xy_xy_0_zz[i] * a_exp + 4.0 * g_xyzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (948-954)

    #pragma omp simd aligned(g_xy_xz_0_xx, g_xy_xz_0_xy, g_xy_xz_0_xz, g_xy_xz_0_yy, g_xy_xz_0_yz, g_xy_xz_0_zz, g_xyzz_xz_0_xx, g_xyzz_xz_0_xy, g_xyzz_xz_0_xz, g_xyzz_xz_0_yy, g_xyzz_xz_0_yz, g_xyzz_xz_0_zz, g_yz_0_0_0_xz_xz_0_xx, g_yz_0_0_0_xz_xz_0_xy, g_yz_0_0_0_xz_xz_0_xz, g_yz_0_0_0_xz_xz_0_yy, g_yz_0_0_0_xz_xz_0_yz, g_yz_0_0_0_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_xz_0_xx[i] = -2.0 * g_xy_xz_0_xx[i] * a_exp + 4.0 * g_xyzz_xz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xz_0_xy[i] = -2.0 * g_xy_xz_0_xy[i] * a_exp + 4.0 * g_xyzz_xz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xz_0_xz[i] = -2.0 * g_xy_xz_0_xz[i] * a_exp + 4.0 * g_xyzz_xz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xz_0_yy[i] = -2.0 * g_xy_xz_0_yy[i] * a_exp + 4.0 * g_xyzz_xz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xz_0_yz[i] = -2.0 * g_xy_xz_0_yz[i] * a_exp + 4.0 * g_xyzz_xz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_xz_0_zz[i] = -2.0 * g_xy_xz_0_zz[i] * a_exp + 4.0 * g_xyzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (954-960)

    #pragma omp simd aligned(g_xy_yy_0_xx, g_xy_yy_0_xy, g_xy_yy_0_xz, g_xy_yy_0_yy, g_xy_yy_0_yz, g_xy_yy_0_zz, g_xyzz_yy_0_xx, g_xyzz_yy_0_xy, g_xyzz_yy_0_xz, g_xyzz_yy_0_yy, g_xyzz_yy_0_yz, g_xyzz_yy_0_zz, g_yz_0_0_0_xz_yy_0_xx, g_yz_0_0_0_xz_yy_0_xy, g_yz_0_0_0_xz_yy_0_xz, g_yz_0_0_0_xz_yy_0_yy, g_yz_0_0_0_xz_yy_0_yz, g_yz_0_0_0_xz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_yy_0_xx[i] = -2.0 * g_xy_yy_0_xx[i] * a_exp + 4.0 * g_xyzz_yy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yy_0_xy[i] = -2.0 * g_xy_yy_0_xy[i] * a_exp + 4.0 * g_xyzz_yy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yy_0_xz[i] = -2.0 * g_xy_yy_0_xz[i] * a_exp + 4.0 * g_xyzz_yy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yy_0_yy[i] = -2.0 * g_xy_yy_0_yy[i] * a_exp + 4.0 * g_xyzz_yy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yy_0_yz[i] = -2.0 * g_xy_yy_0_yz[i] * a_exp + 4.0 * g_xyzz_yy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yy_0_zz[i] = -2.0 * g_xy_yy_0_zz[i] * a_exp + 4.0 * g_xyzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (960-966)

    #pragma omp simd aligned(g_xy_yz_0_xx, g_xy_yz_0_xy, g_xy_yz_0_xz, g_xy_yz_0_yy, g_xy_yz_0_yz, g_xy_yz_0_zz, g_xyzz_yz_0_xx, g_xyzz_yz_0_xy, g_xyzz_yz_0_xz, g_xyzz_yz_0_yy, g_xyzz_yz_0_yz, g_xyzz_yz_0_zz, g_yz_0_0_0_xz_yz_0_xx, g_yz_0_0_0_xz_yz_0_xy, g_yz_0_0_0_xz_yz_0_xz, g_yz_0_0_0_xz_yz_0_yy, g_yz_0_0_0_xz_yz_0_yz, g_yz_0_0_0_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_yz_0_xx[i] = -2.0 * g_xy_yz_0_xx[i] * a_exp + 4.0 * g_xyzz_yz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yz_0_xy[i] = -2.0 * g_xy_yz_0_xy[i] * a_exp + 4.0 * g_xyzz_yz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yz_0_xz[i] = -2.0 * g_xy_yz_0_xz[i] * a_exp + 4.0 * g_xyzz_yz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yz_0_yy[i] = -2.0 * g_xy_yz_0_yy[i] * a_exp + 4.0 * g_xyzz_yz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yz_0_yz[i] = -2.0 * g_xy_yz_0_yz[i] * a_exp + 4.0 * g_xyzz_yz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_yz_0_zz[i] = -2.0 * g_xy_yz_0_zz[i] * a_exp + 4.0 * g_xyzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (966-972)

    #pragma omp simd aligned(g_xy_zz_0_xx, g_xy_zz_0_xy, g_xy_zz_0_xz, g_xy_zz_0_yy, g_xy_zz_0_yz, g_xy_zz_0_zz, g_xyzz_zz_0_xx, g_xyzz_zz_0_xy, g_xyzz_zz_0_xz, g_xyzz_zz_0_yy, g_xyzz_zz_0_yz, g_xyzz_zz_0_zz, g_yz_0_0_0_xz_zz_0_xx, g_yz_0_0_0_xz_zz_0_xy, g_yz_0_0_0_xz_zz_0_xz, g_yz_0_0_0_xz_zz_0_yy, g_yz_0_0_0_xz_zz_0_yz, g_yz_0_0_0_xz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_xz_zz_0_xx[i] = -2.0 * g_xy_zz_0_xx[i] * a_exp + 4.0 * g_xyzz_zz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_zz_0_xy[i] = -2.0 * g_xy_zz_0_xy[i] * a_exp + 4.0 * g_xyzz_zz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_zz_0_xz[i] = -2.0 * g_xy_zz_0_xz[i] * a_exp + 4.0 * g_xyzz_zz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_zz_0_yy[i] = -2.0 * g_xy_zz_0_yy[i] * a_exp + 4.0 * g_xyzz_zz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_zz_0_yz[i] = -2.0 * g_xy_zz_0_yz[i] * a_exp + 4.0 * g_xyzz_zz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_xz_zz_0_zz[i] = -2.0 * g_xy_zz_0_zz[i] * a_exp + 4.0 * g_xyzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (972-978)

    #pragma omp simd aligned(g_yyyz_xx_0_xx, g_yyyz_xx_0_xy, g_yyyz_xx_0_xz, g_yyyz_xx_0_yy, g_yyyz_xx_0_yz, g_yyyz_xx_0_zz, g_yz_0_0_0_yy_xx_0_xx, g_yz_0_0_0_yy_xx_0_xy, g_yz_0_0_0_yy_xx_0_xz, g_yz_0_0_0_yy_xx_0_yy, g_yz_0_0_0_yy_xx_0_yz, g_yz_0_0_0_yy_xx_0_zz, g_yz_xx_0_xx, g_yz_xx_0_xy, g_yz_xx_0_xz, g_yz_xx_0_yy, g_yz_xx_0_yz, g_yz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_xx_0_xx[i] = -4.0 * g_yz_xx_0_xx[i] * a_exp + 4.0 * g_yyyz_xx_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xx_0_xy[i] = -4.0 * g_yz_xx_0_xy[i] * a_exp + 4.0 * g_yyyz_xx_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xx_0_xz[i] = -4.0 * g_yz_xx_0_xz[i] * a_exp + 4.0 * g_yyyz_xx_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xx_0_yy[i] = -4.0 * g_yz_xx_0_yy[i] * a_exp + 4.0 * g_yyyz_xx_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xx_0_yz[i] = -4.0 * g_yz_xx_0_yz[i] * a_exp + 4.0 * g_yyyz_xx_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xx_0_zz[i] = -4.0 * g_yz_xx_0_zz[i] * a_exp + 4.0 * g_yyyz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (978-984)

    #pragma omp simd aligned(g_yyyz_xy_0_xx, g_yyyz_xy_0_xy, g_yyyz_xy_0_xz, g_yyyz_xy_0_yy, g_yyyz_xy_0_yz, g_yyyz_xy_0_zz, g_yz_0_0_0_yy_xy_0_xx, g_yz_0_0_0_yy_xy_0_xy, g_yz_0_0_0_yy_xy_0_xz, g_yz_0_0_0_yy_xy_0_yy, g_yz_0_0_0_yy_xy_0_yz, g_yz_0_0_0_yy_xy_0_zz, g_yz_xy_0_xx, g_yz_xy_0_xy, g_yz_xy_0_xz, g_yz_xy_0_yy, g_yz_xy_0_yz, g_yz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_xy_0_xx[i] = -4.0 * g_yz_xy_0_xx[i] * a_exp + 4.0 * g_yyyz_xy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xy_0_xy[i] = -4.0 * g_yz_xy_0_xy[i] * a_exp + 4.0 * g_yyyz_xy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xy_0_xz[i] = -4.0 * g_yz_xy_0_xz[i] * a_exp + 4.0 * g_yyyz_xy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xy_0_yy[i] = -4.0 * g_yz_xy_0_yy[i] * a_exp + 4.0 * g_yyyz_xy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xy_0_yz[i] = -4.0 * g_yz_xy_0_yz[i] * a_exp + 4.0 * g_yyyz_xy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xy_0_zz[i] = -4.0 * g_yz_xy_0_zz[i] * a_exp + 4.0 * g_yyyz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (984-990)

    #pragma omp simd aligned(g_yyyz_xz_0_xx, g_yyyz_xz_0_xy, g_yyyz_xz_0_xz, g_yyyz_xz_0_yy, g_yyyz_xz_0_yz, g_yyyz_xz_0_zz, g_yz_0_0_0_yy_xz_0_xx, g_yz_0_0_0_yy_xz_0_xy, g_yz_0_0_0_yy_xz_0_xz, g_yz_0_0_0_yy_xz_0_yy, g_yz_0_0_0_yy_xz_0_yz, g_yz_0_0_0_yy_xz_0_zz, g_yz_xz_0_xx, g_yz_xz_0_xy, g_yz_xz_0_xz, g_yz_xz_0_yy, g_yz_xz_0_yz, g_yz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_xz_0_xx[i] = -4.0 * g_yz_xz_0_xx[i] * a_exp + 4.0 * g_yyyz_xz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xz_0_xy[i] = -4.0 * g_yz_xz_0_xy[i] * a_exp + 4.0 * g_yyyz_xz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xz_0_xz[i] = -4.0 * g_yz_xz_0_xz[i] * a_exp + 4.0 * g_yyyz_xz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xz_0_yy[i] = -4.0 * g_yz_xz_0_yy[i] * a_exp + 4.0 * g_yyyz_xz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xz_0_yz[i] = -4.0 * g_yz_xz_0_yz[i] * a_exp + 4.0 * g_yyyz_xz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_xz_0_zz[i] = -4.0 * g_yz_xz_0_zz[i] * a_exp + 4.0 * g_yyyz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (990-996)

    #pragma omp simd aligned(g_yyyz_yy_0_xx, g_yyyz_yy_0_xy, g_yyyz_yy_0_xz, g_yyyz_yy_0_yy, g_yyyz_yy_0_yz, g_yyyz_yy_0_zz, g_yz_0_0_0_yy_yy_0_xx, g_yz_0_0_0_yy_yy_0_xy, g_yz_0_0_0_yy_yy_0_xz, g_yz_0_0_0_yy_yy_0_yy, g_yz_0_0_0_yy_yy_0_yz, g_yz_0_0_0_yy_yy_0_zz, g_yz_yy_0_xx, g_yz_yy_0_xy, g_yz_yy_0_xz, g_yz_yy_0_yy, g_yz_yy_0_yz, g_yz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_yy_0_xx[i] = -4.0 * g_yz_yy_0_xx[i] * a_exp + 4.0 * g_yyyz_yy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yy_0_xy[i] = -4.0 * g_yz_yy_0_xy[i] * a_exp + 4.0 * g_yyyz_yy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yy_0_xz[i] = -4.0 * g_yz_yy_0_xz[i] * a_exp + 4.0 * g_yyyz_yy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yy_0_yy[i] = -4.0 * g_yz_yy_0_yy[i] * a_exp + 4.0 * g_yyyz_yy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yy_0_yz[i] = -4.0 * g_yz_yy_0_yz[i] * a_exp + 4.0 * g_yyyz_yy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yy_0_zz[i] = -4.0 * g_yz_yy_0_zz[i] * a_exp + 4.0 * g_yyyz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (996-1002)

    #pragma omp simd aligned(g_yyyz_yz_0_xx, g_yyyz_yz_0_xy, g_yyyz_yz_0_xz, g_yyyz_yz_0_yy, g_yyyz_yz_0_yz, g_yyyz_yz_0_zz, g_yz_0_0_0_yy_yz_0_xx, g_yz_0_0_0_yy_yz_0_xy, g_yz_0_0_0_yy_yz_0_xz, g_yz_0_0_0_yy_yz_0_yy, g_yz_0_0_0_yy_yz_0_yz, g_yz_0_0_0_yy_yz_0_zz, g_yz_yz_0_xx, g_yz_yz_0_xy, g_yz_yz_0_xz, g_yz_yz_0_yy, g_yz_yz_0_yz, g_yz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_yz_0_xx[i] = -4.0 * g_yz_yz_0_xx[i] * a_exp + 4.0 * g_yyyz_yz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yz_0_xy[i] = -4.0 * g_yz_yz_0_xy[i] * a_exp + 4.0 * g_yyyz_yz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yz_0_xz[i] = -4.0 * g_yz_yz_0_xz[i] * a_exp + 4.0 * g_yyyz_yz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yz_0_yy[i] = -4.0 * g_yz_yz_0_yy[i] * a_exp + 4.0 * g_yyyz_yz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yz_0_yz[i] = -4.0 * g_yz_yz_0_yz[i] * a_exp + 4.0 * g_yyyz_yz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_yz_0_zz[i] = -4.0 * g_yz_yz_0_zz[i] * a_exp + 4.0 * g_yyyz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1002-1008)

    #pragma omp simd aligned(g_yyyz_zz_0_xx, g_yyyz_zz_0_xy, g_yyyz_zz_0_xz, g_yyyz_zz_0_yy, g_yyyz_zz_0_yz, g_yyyz_zz_0_zz, g_yz_0_0_0_yy_zz_0_xx, g_yz_0_0_0_yy_zz_0_xy, g_yz_0_0_0_yy_zz_0_xz, g_yz_0_0_0_yy_zz_0_yy, g_yz_0_0_0_yy_zz_0_yz, g_yz_0_0_0_yy_zz_0_zz, g_yz_zz_0_xx, g_yz_zz_0_xy, g_yz_zz_0_xz, g_yz_zz_0_yy, g_yz_zz_0_yz, g_yz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yy_zz_0_xx[i] = -4.0 * g_yz_zz_0_xx[i] * a_exp + 4.0 * g_yyyz_zz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_zz_0_xy[i] = -4.0 * g_yz_zz_0_xy[i] * a_exp + 4.0 * g_yyyz_zz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_zz_0_xz[i] = -4.0 * g_yz_zz_0_xz[i] * a_exp + 4.0 * g_yyyz_zz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_zz_0_yy[i] = -4.0 * g_yz_zz_0_yy[i] * a_exp + 4.0 * g_yyyz_zz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_zz_0_yz[i] = -4.0 * g_yz_zz_0_yz[i] * a_exp + 4.0 * g_yyyz_zz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_yy_zz_0_zz[i] = -4.0 * g_yz_zz_0_zz[i] * a_exp + 4.0 * g_yyyz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1008-1014)

    #pragma omp simd aligned(g_0_xx_0_xx, g_0_xx_0_xy, g_0_xx_0_xz, g_0_xx_0_yy, g_0_xx_0_yz, g_0_xx_0_zz, g_yy_xx_0_xx, g_yy_xx_0_xy, g_yy_xx_0_xz, g_yy_xx_0_yy, g_yy_xx_0_yz, g_yy_xx_0_zz, g_yyzz_xx_0_xx, g_yyzz_xx_0_xy, g_yyzz_xx_0_xz, g_yyzz_xx_0_yy, g_yyzz_xx_0_yz, g_yyzz_xx_0_zz, g_yz_0_0_0_yz_xx_0_xx, g_yz_0_0_0_yz_xx_0_xy, g_yz_0_0_0_yz_xx_0_xz, g_yz_0_0_0_yz_xx_0_yy, g_yz_0_0_0_yz_xx_0_yz, g_yz_0_0_0_yz_xx_0_zz, g_zz_xx_0_xx, g_zz_xx_0_xy, g_zz_xx_0_xz, g_zz_xx_0_yy, g_zz_xx_0_yz, g_zz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_xx_0_xx[i] = g_0_xx_0_xx[i] - 2.0 * g_zz_xx_0_xx[i] * a_exp - 2.0 * g_yy_xx_0_xx[i] * a_exp + 4.0 * g_yyzz_xx_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xx_0_xy[i] = g_0_xx_0_xy[i] - 2.0 * g_zz_xx_0_xy[i] * a_exp - 2.0 * g_yy_xx_0_xy[i] * a_exp + 4.0 * g_yyzz_xx_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xx_0_xz[i] = g_0_xx_0_xz[i] - 2.0 * g_zz_xx_0_xz[i] * a_exp - 2.0 * g_yy_xx_0_xz[i] * a_exp + 4.0 * g_yyzz_xx_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xx_0_yy[i] = g_0_xx_0_yy[i] - 2.0 * g_zz_xx_0_yy[i] * a_exp - 2.0 * g_yy_xx_0_yy[i] * a_exp + 4.0 * g_yyzz_xx_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xx_0_yz[i] = g_0_xx_0_yz[i] - 2.0 * g_zz_xx_0_yz[i] * a_exp - 2.0 * g_yy_xx_0_yz[i] * a_exp + 4.0 * g_yyzz_xx_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xx_0_zz[i] = g_0_xx_0_zz[i] - 2.0 * g_zz_xx_0_zz[i] * a_exp - 2.0 * g_yy_xx_0_zz[i] * a_exp + 4.0 * g_yyzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1014-1020)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_yy_xy_0_xx, g_yy_xy_0_xy, g_yy_xy_0_xz, g_yy_xy_0_yy, g_yy_xy_0_yz, g_yy_xy_0_zz, g_yyzz_xy_0_xx, g_yyzz_xy_0_xy, g_yyzz_xy_0_xz, g_yyzz_xy_0_yy, g_yyzz_xy_0_yz, g_yyzz_xy_0_zz, g_yz_0_0_0_yz_xy_0_xx, g_yz_0_0_0_yz_xy_0_xy, g_yz_0_0_0_yz_xy_0_xz, g_yz_0_0_0_yz_xy_0_yy, g_yz_0_0_0_yz_xy_0_yz, g_yz_0_0_0_yz_xy_0_zz, g_zz_xy_0_xx, g_zz_xy_0_xy, g_zz_xy_0_xz, g_zz_xy_0_yy, g_zz_xy_0_yz, g_zz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_xy_0_xx[i] = g_0_xy_0_xx[i] - 2.0 * g_zz_xy_0_xx[i] * a_exp - 2.0 * g_yy_xy_0_xx[i] * a_exp + 4.0 * g_yyzz_xy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xy_0_xy[i] = g_0_xy_0_xy[i] - 2.0 * g_zz_xy_0_xy[i] * a_exp - 2.0 * g_yy_xy_0_xy[i] * a_exp + 4.0 * g_yyzz_xy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xy_0_xz[i] = g_0_xy_0_xz[i] - 2.0 * g_zz_xy_0_xz[i] * a_exp - 2.0 * g_yy_xy_0_xz[i] * a_exp + 4.0 * g_yyzz_xy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xy_0_yy[i] = g_0_xy_0_yy[i] - 2.0 * g_zz_xy_0_yy[i] * a_exp - 2.0 * g_yy_xy_0_yy[i] * a_exp + 4.0 * g_yyzz_xy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xy_0_yz[i] = g_0_xy_0_yz[i] - 2.0 * g_zz_xy_0_yz[i] * a_exp - 2.0 * g_yy_xy_0_yz[i] * a_exp + 4.0 * g_yyzz_xy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xy_0_zz[i] = g_0_xy_0_zz[i] - 2.0 * g_zz_xy_0_zz[i] * a_exp - 2.0 * g_yy_xy_0_zz[i] * a_exp + 4.0 * g_yyzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1020-1026)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_yy_xz_0_xx, g_yy_xz_0_xy, g_yy_xz_0_xz, g_yy_xz_0_yy, g_yy_xz_0_yz, g_yy_xz_0_zz, g_yyzz_xz_0_xx, g_yyzz_xz_0_xy, g_yyzz_xz_0_xz, g_yyzz_xz_0_yy, g_yyzz_xz_0_yz, g_yyzz_xz_0_zz, g_yz_0_0_0_yz_xz_0_xx, g_yz_0_0_0_yz_xz_0_xy, g_yz_0_0_0_yz_xz_0_xz, g_yz_0_0_0_yz_xz_0_yy, g_yz_0_0_0_yz_xz_0_yz, g_yz_0_0_0_yz_xz_0_zz, g_zz_xz_0_xx, g_zz_xz_0_xy, g_zz_xz_0_xz, g_zz_xz_0_yy, g_zz_xz_0_yz, g_zz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_xz_0_xx[i] = g_0_xz_0_xx[i] - 2.0 * g_zz_xz_0_xx[i] * a_exp - 2.0 * g_yy_xz_0_xx[i] * a_exp + 4.0 * g_yyzz_xz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xz_0_xy[i] = g_0_xz_0_xy[i] - 2.0 * g_zz_xz_0_xy[i] * a_exp - 2.0 * g_yy_xz_0_xy[i] * a_exp + 4.0 * g_yyzz_xz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xz_0_xz[i] = g_0_xz_0_xz[i] - 2.0 * g_zz_xz_0_xz[i] * a_exp - 2.0 * g_yy_xz_0_xz[i] * a_exp + 4.0 * g_yyzz_xz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xz_0_yy[i] = g_0_xz_0_yy[i] - 2.0 * g_zz_xz_0_yy[i] * a_exp - 2.0 * g_yy_xz_0_yy[i] * a_exp + 4.0 * g_yyzz_xz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xz_0_yz[i] = g_0_xz_0_yz[i] - 2.0 * g_zz_xz_0_yz[i] * a_exp - 2.0 * g_yy_xz_0_yz[i] * a_exp + 4.0 * g_yyzz_xz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_xz_0_zz[i] = g_0_xz_0_zz[i] - 2.0 * g_zz_xz_0_zz[i] * a_exp - 2.0 * g_yy_xz_0_zz[i] * a_exp + 4.0 * g_yyzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1026-1032)

    #pragma omp simd aligned(g_0_yy_0_xx, g_0_yy_0_xy, g_0_yy_0_xz, g_0_yy_0_yy, g_0_yy_0_yz, g_0_yy_0_zz, g_yy_yy_0_xx, g_yy_yy_0_xy, g_yy_yy_0_xz, g_yy_yy_0_yy, g_yy_yy_0_yz, g_yy_yy_0_zz, g_yyzz_yy_0_xx, g_yyzz_yy_0_xy, g_yyzz_yy_0_xz, g_yyzz_yy_0_yy, g_yyzz_yy_0_yz, g_yyzz_yy_0_zz, g_yz_0_0_0_yz_yy_0_xx, g_yz_0_0_0_yz_yy_0_xy, g_yz_0_0_0_yz_yy_0_xz, g_yz_0_0_0_yz_yy_0_yy, g_yz_0_0_0_yz_yy_0_yz, g_yz_0_0_0_yz_yy_0_zz, g_zz_yy_0_xx, g_zz_yy_0_xy, g_zz_yy_0_xz, g_zz_yy_0_yy, g_zz_yy_0_yz, g_zz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_yy_0_xx[i] = g_0_yy_0_xx[i] - 2.0 * g_zz_yy_0_xx[i] * a_exp - 2.0 * g_yy_yy_0_xx[i] * a_exp + 4.0 * g_yyzz_yy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yy_0_xy[i] = g_0_yy_0_xy[i] - 2.0 * g_zz_yy_0_xy[i] * a_exp - 2.0 * g_yy_yy_0_xy[i] * a_exp + 4.0 * g_yyzz_yy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yy_0_xz[i] = g_0_yy_0_xz[i] - 2.0 * g_zz_yy_0_xz[i] * a_exp - 2.0 * g_yy_yy_0_xz[i] * a_exp + 4.0 * g_yyzz_yy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yy_0_yy[i] = g_0_yy_0_yy[i] - 2.0 * g_zz_yy_0_yy[i] * a_exp - 2.0 * g_yy_yy_0_yy[i] * a_exp + 4.0 * g_yyzz_yy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yy_0_yz[i] = g_0_yy_0_yz[i] - 2.0 * g_zz_yy_0_yz[i] * a_exp - 2.0 * g_yy_yy_0_yz[i] * a_exp + 4.0 * g_yyzz_yy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yy_0_zz[i] = g_0_yy_0_zz[i] - 2.0 * g_zz_yy_0_zz[i] * a_exp - 2.0 * g_yy_yy_0_zz[i] * a_exp + 4.0 * g_yyzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1032-1038)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_yy_yz_0_xx, g_yy_yz_0_xy, g_yy_yz_0_xz, g_yy_yz_0_yy, g_yy_yz_0_yz, g_yy_yz_0_zz, g_yyzz_yz_0_xx, g_yyzz_yz_0_xy, g_yyzz_yz_0_xz, g_yyzz_yz_0_yy, g_yyzz_yz_0_yz, g_yyzz_yz_0_zz, g_yz_0_0_0_yz_yz_0_xx, g_yz_0_0_0_yz_yz_0_xy, g_yz_0_0_0_yz_yz_0_xz, g_yz_0_0_0_yz_yz_0_yy, g_yz_0_0_0_yz_yz_0_yz, g_yz_0_0_0_yz_yz_0_zz, g_zz_yz_0_xx, g_zz_yz_0_xy, g_zz_yz_0_xz, g_zz_yz_0_yy, g_zz_yz_0_yz, g_zz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_yz_0_xx[i] = g_0_yz_0_xx[i] - 2.0 * g_zz_yz_0_xx[i] * a_exp - 2.0 * g_yy_yz_0_xx[i] * a_exp + 4.0 * g_yyzz_yz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yz_0_xy[i] = g_0_yz_0_xy[i] - 2.0 * g_zz_yz_0_xy[i] * a_exp - 2.0 * g_yy_yz_0_xy[i] * a_exp + 4.0 * g_yyzz_yz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yz_0_xz[i] = g_0_yz_0_xz[i] - 2.0 * g_zz_yz_0_xz[i] * a_exp - 2.0 * g_yy_yz_0_xz[i] * a_exp + 4.0 * g_yyzz_yz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yz_0_yy[i] = g_0_yz_0_yy[i] - 2.0 * g_zz_yz_0_yy[i] * a_exp - 2.0 * g_yy_yz_0_yy[i] * a_exp + 4.0 * g_yyzz_yz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yz_0_yz[i] = g_0_yz_0_yz[i] - 2.0 * g_zz_yz_0_yz[i] * a_exp - 2.0 * g_yy_yz_0_yz[i] * a_exp + 4.0 * g_yyzz_yz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_yz_0_zz[i] = g_0_yz_0_zz[i] - 2.0 * g_zz_yz_0_zz[i] * a_exp - 2.0 * g_yy_yz_0_zz[i] * a_exp + 4.0 * g_yyzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1038-1044)

    #pragma omp simd aligned(g_0_zz_0_xx, g_0_zz_0_xy, g_0_zz_0_xz, g_0_zz_0_yy, g_0_zz_0_yz, g_0_zz_0_zz, g_yy_zz_0_xx, g_yy_zz_0_xy, g_yy_zz_0_xz, g_yy_zz_0_yy, g_yy_zz_0_yz, g_yy_zz_0_zz, g_yyzz_zz_0_xx, g_yyzz_zz_0_xy, g_yyzz_zz_0_xz, g_yyzz_zz_0_yy, g_yyzz_zz_0_yz, g_yyzz_zz_0_zz, g_yz_0_0_0_yz_zz_0_xx, g_yz_0_0_0_yz_zz_0_xy, g_yz_0_0_0_yz_zz_0_xz, g_yz_0_0_0_yz_zz_0_yy, g_yz_0_0_0_yz_zz_0_yz, g_yz_0_0_0_yz_zz_0_zz, g_zz_zz_0_xx, g_zz_zz_0_xy, g_zz_zz_0_xz, g_zz_zz_0_yy, g_zz_zz_0_yz, g_zz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_yz_zz_0_xx[i] = g_0_zz_0_xx[i] - 2.0 * g_zz_zz_0_xx[i] * a_exp - 2.0 * g_yy_zz_0_xx[i] * a_exp + 4.0 * g_yyzz_zz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_zz_0_xy[i] = g_0_zz_0_xy[i] - 2.0 * g_zz_zz_0_xy[i] * a_exp - 2.0 * g_yy_zz_0_xy[i] * a_exp + 4.0 * g_yyzz_zz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_zz_0_xz[i] = g_0_zz_0_xz[i] - 2.0 * g_zz_zz_0_xz[i] * a_exp - 2.0 * g_yy_zz_0_xz[i] * a_exp + 4.0 * g_yyzz_zz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_zz_0_yy[i] = g_0_zz_0_yy[i] - 2.0 * g_zz_zz_0_yy[i] * a_exp - 2.0 * g_yy_zz_0_yy[i] * a_exp + 4.0 * g_yyzz_zz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_zz_0_yz[i] = g_0_zz_0_yz[i] - 2.0 * g_zz_zz_0_yz[i] * a_exp - 2.0 * g_yy_zz_0_yz[i] * a_exp + 4.0 * g_yyzz_zz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_yz_zz_0_zz[i] = g_0_zz_0_zz[i] - 2.0 * g_zz_zz_0_zz[i] * a_exp - 2.0 * g_yy_zz_0_zz[i] * a_exp + 4.0 * g_yyzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1044-1050)

    #pragma omp simd aligned(g_yz_0_0_0_zz_xx_0_xx, g_yz_0_0_0_zz_xx_0_xy, g_yz_0_0_0_zz_xx_0_xz, g_yz_0_0_0_zz_xx_0_yy, g_yz_0_0_0_zz_xx_0_yz, g_yz_0_0_0_zz_xx_0_zz, g_yz_xx_0_xx, g_yz_xx_0_xy, g_yz_xx_0_xz, g_yz_xx_0_yy, g_yz_xx_0_yz, g_yz_xx_0_zz, g_yzzz_xx_0_xx, g_yzzz_xx_0_xy, g_yzzz_xx_0_xz, g_yzzz_xx_0_yy, g_yzzz_xx_0_yz, g_yzzz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_xx_0_xx[i] = -4.0 * g_yz_xx_0_xx[i] * a_exp + 4.0 * g_yzzz_xx_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xx_0_xy[i] = -4.0 * g_yz_xx_0_xy[i] * a_exp + 4.0 * g_yzzz_xx_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xx_0_xz[i] = -4.0 * g_yz_xx_0_xz[i] * a_exp + 4.0 * g_yzzz_xx_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xx_0_yy[i] = -4.0 * g_yz_xx_0_yy[i] * a_exp + 4.0 * g_yzzz_xx_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xx_0_yz[i] = -4.0 * g_yz_xx_0_yz[i] * a_exp + 4.0 * g_yzzz_xx_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xx_0_zz[i] = -4.0 * g_yz_xx_0_zz[i] * a_exp + 4.0 * g_yzzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1050-1056)

    #pragma omp simd aligned(g_yz_0_0_0_zz_xy_0_xx, g_yz_0_0_0_zz_xy_0_xy, g_yz_0_0_0_zz_xy_0_xz, g_yz_0_0_0_zz_xy_0_yy, g_yz_0_0_0_zz_xy_0_yz, g_yz_0_0_0_zz_xy_0_zz, g_yz_xy_0_xx, g_yz_xy_0_xy, g_yz_xy_0_xz, g_yz_xy_0_yy, g_yz_xy_0_yz, g_yz_xy_0_zz, g_yzzz_xy_0_xx, g_yzzz_xy_0_xy, g_yzzz_xy_0_xz, g_yzzz_xy_0_yy, g_yzzz_xy_0_yz, g_yzzz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_xy_0_xx[i] = -4.0 * g_yz_xy_0_xx[i] * a_exp + 4.0 * g_yzzz_xy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xy_0_xy[i] = -4.0 * g_yz_xy_0_xy[i] * a_exp + 4.0 * g_yzzz_xy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xy_0_xz[i] = -4.0 * g_yz_xy_0_xz[i] * a_exp + 4.0 * g_yzzz_xy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xy_0_yy[i] = -4.0 * g_yz_xy_0_yy[i] * a_exp + 4.0 * g_yzzz_xy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xy_0_yz[i] = -4.0 * g_yz_xy_0_yz[i] * a_exp + 4.0 * g_yzzz_xy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xy_0_zz[i] = -4.0 * g_yz_xy_0_zz[i] * a_exp + 4.0 * g_yzzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1056-1062)

    #pragma omp simd aligned(g_yz_0_0_0_zz_xz_0_xx, g_yz_0_0_0_zz_xz_0_xy, g_yz_0_0_0_zz_xz_0_xz, g_yz_0_0_0_zz_xz_0_yy, g_yz_0_0_0_zz_xz_0_yz, g_yz_0_0_0_zz_xz_0_zz, g_yz_xz_0_xx, g_yz_xz_0_xy, g_yz_xz_0_xz, g_yz_xz_0_yy, g_yz_xz_0_yz, g_yz_xz_0_zz, g_yzzz_xz_0_xx, g_yzzz_xz_0_xy, g_yzzz_xz_0_xz, g_yzzz_xz_0_yy, g_yzzz_xz_0_yz, g_yzzz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_xz_0_xx[i] = -4.0 * g_yz_xz_0_xx[i] * a_exp + 4.0 * g_yzzz_xz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xz_0_xy[i] = -4.0 * g_yz_xz_0_xy[i] * a_exp + 4.0 * g_yzzz_xz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xz_0_xz[i] = -4.0 * g_yz_xz_0_xz[i] * a_exp + 4.0 * g_yzzz_xz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xz_0_yy[i] = -4.0 * g_yz_xz_0_yy[i] * a_exp + 4.0 * g_yzzz_xz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xz_0_yz[i] = -4.0 * g_yz_xz_0_yz[i] * a_exp + 4.0 * g_yzzz_xz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_xz_0_zz[i] = -4.0 * g_yz_xz_0_zz[i] * a_exp + 4.0 * g_yzzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1062-1068)

    #pragma omp simd aligned(g_yz_0_0_0_zz_yy_0_xx, g_yz_0_0_0_zz_yy_0_xy, g_yz_0_0_0_zz_yy_0_xz, g_yz_0_0_0_zz_yy_0_yy, g_yz_0_0_0_zz_yy_0_yz, g_yz_0_0_0_zz_yy_0_zz, g_yz_yy_0_xx, g_yz_yy_0_xy, g_yz_yy_0_xz, g_yz_yy_0_yy, g_yz_yy_0_yz, g_yz_yy_0_zz, g_yzzz_yy_0_xx, g_yzzz_yy_0_xy, g_yzzz_yy_0_xz, g_yzzz_yy_0_yy, g_yzzz_yy_0_yz, g_yzzz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_yy_0_xx[i] = -4.0 * g_yz_yy_0_xx[i] * a_exp + 4.0 * g_yzzz_yy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yy_0_xy[i] = -4.0 * g_yz_yy_0_xy[i] * a_exp + 4.0 * g_yzzz_yy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yy_0_xz[i] = -4.0 * g_yz_yy_0_xz[i] * a_exp + 4.0 * g_yzzz_yy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yy_0_yy[i] = -4.0 * g_yz_yy_0_yy[i] * a_exp + 4.0 * g_yzzz_yy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yy_0_yz[i] = -4.0 * g_yz_yy_0_yz[i] * a_exp + 4.0 * g_yzzz_yy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yy_0_zz[i] = -4.0 * g_yz_yy_0_zz[i] * a_exp + 4.0 * g_yzzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1068-1074)

    #pragma omp simd aligned(g_yz_0_0_0_zz_yz_0_xx, g_yz_0_0_0_zz_yz_0_xy, g_yz_0_0_0_zz_yz_0_xz, g_yz_0_0_0_zz_yz_0_yy, g_yz_0_0_0_zz_yz_0_yz, g_yz_0_0_0_zz_yz_0_zz, g_yz_yz_0_xx, g_yz_yz_0_xy, g_yz_yz_0_xz, g_yz_yz_0_yy, g_yz_yz_0_yz, g_yz_yz_0_zz, g_yzzz_yz_0_xx, g_yzzz_yz_0_xy, g_yzzz_yz_0_xz, g_yzzz_yz_0_yy, g_yzzz_yz_0_yz, g_yzzz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_yz_0_xx[i] = -4.0 * g_yz_yz_0_xx[i] * a_exp + 4.0 * g_yzzz_yz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yz_0_xy[i] = -4.0 * g_yz_yz_0_xy[i] * a_exp + 4.0 * g_yzzz_yz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yz_0_xz[i] = -4.0 * g_yz_yz_0_xz[i] * a_exp + 4.0 * g_yzzz_yz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yz_0_yy[i] = -4.0 * g_yz_yz_0_yy[i] * a_exp + 4.0 * g_yzzz_yz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yz_0_yz[i] = -4.0 * g_yz_yz_0_yz[i] * a_exp + 4.0 * g_yzzz_yz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_yz_0_zz[i] = -4.0 * g_yz_yz_0_zz[i] * a_exp + 4.0 * g_yzzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1074-1080)

    #pragma omp simd aligned(g_yz_0_0_0_zz_zz_0_xx, g_yz_0_0_0_zz_zz_0_xy, g_yz_0_0_0_zz_zz_0_xz, g_yz_0_0_0_zz_zz_0_yy, g_yz_0_0_0_zz_zz_0_yz, g_yz_0_0_0_zz_zz_0_zz, g_yz_zz_0_xx, g_yz_zz_0_xy, g_yz_zz_0_xz, g_yz_zz_0_yy, g_yz_zz_0_yz, g_yz_zz_0_zz, g_yzzz_zz_0_xx, g_yzzz_zz_0_xy, g_yzzz_zz_0_xz, g_yzzz_zz_0_yy, g_yzzz_zz_0_yz, g_yzzz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_zz_zz_0_xx[i] = -4.0 * g_yz_zz_0_xx[i] * a_exp + 4.0 * g_yzzz_zz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_zz_0_xy[i] = -4.0 * g_yz_zz_0_xy[i] * a_exp + 4.0 * g_yzzz_zz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_zz_0_xz[i] = -4.0 * g_yz_zz_0_xz[i] * a_exp + 4.0 * g_yzzz_zz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_zz_0_yy[i] = -4.0 * g_yz_zz_0_yy[i] * a_exp + 4.0 * g_yzzz_zz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_zz_0_yz[i] = -4.0 * g_yz_zz_0_yz[i] * a_exp + 4.0 * g_yzzz_zz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_zz_zz_0_zz[i] = -4.0 * g_yz_zz_0_zz[i] * a_exp + 4.0 * g_yzzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1080-1086)

    #pragma omp simd aligned(g_xx_xx_0_xx, g_xx_xx_0_xy, g_xx_xx_0_xz, g_xx_xx_0_yy, g_xx_xx_0_yz, g_xx_xx_0_zz, g_xxzz_xx_0_xx, g_xxzz_xx_0_xy, g_xxzz_xx_0_xz, g_xxzz_xx_0_yy, g_xxzz_xx_0_yz, g_xxzz_xx_0_zz, g_zz_0_0_0_xx_xx_0_xx, g_zz_0_0_0_xx_xx_0_xy, g_zz_0_0_0_xx_xx_0_xz, g_zz_0_0_0_xx_xx_0_yy, g_zz_0_0_0_xx_xx_0_yz, g_zz_0_0_0_xx_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_xx_0_xx[i] = -2.0 * g_xx_xx_0_xx[i] * a_exp + 4.0 * g_xxzz_xx_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xx_0_xy[i] = -2.0 * g_xx_xx_0_xy[i] * a_exp + 4.0 * g_xxzz_xx_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xx_0_xz[i] = -2.0 * g_xx_xx_0_xz[i] * a_exp + 4.0 * g_xxzz_xx_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xx_0_yy[i] = -2.0 * g_xx_xx_0_yy[i] * a_exp + 4.0 * g_xxzz_xx_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xx_0_yz[i] = -2.0 * g_xx_xx_0_yz[i] * a_exp + 4.0 * g_xxzz_xx_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xx_0_zz[i] = -2.0 * g_xx_xx_0_zz[i] * a_exp + 4.0 * g_xxzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1086-1092)

    #pragma omp simd aligned(g_xx_xy_0_xx, g_xx_xy_0_xy, g_xx_xy_0_xz, g_xx_xy_0_yy, g_xx_xy_0_yz, g_xx_xy_0_zz, g_xxzz_xy_0_xx, g_xxzz_xy_0_xy, g_xxzz_xy_0_xz, g_xxzz_xy_0_yy, g_xxzz_xy_0_yz, g_xxzz_xy_0_zz, g_zz_0_0_0_xx_xy_0_xx, g_zz_0_0_0_xx_xy_0_xy, g_zz_0_0_0_xx_xy_0_xz, g_zz_0_0_0_xx_xy_0_yy, g_zz_0_0_0_xx_xy_0_yz, g_zz_0_0_0_xx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_xy_0_xx[i] = -2.0 * g_xx_xy_0_xx[i] * a_exp + 4.0 * g_xxzz_xy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xy_0_xy[i] = -2.0 * g_xx_xy_0_xy[i] * a_exp + 4.0 * g_xxzz_xy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xy_0_xz[i] = -2.0 * g_xx_xy_0_xz[i] * a_exp + 4.0 * g_xxzz_xy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xy_0_yy[i] = -2.0 * g_xx_xy_0_yy[i] * a_exp + 4.0 * g_xxzz_xy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xy_0_yz[i] = -2.0 * g_xx_xy_0_yz[i] * a_exp + 4.0 * g_xxzz_xy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xy_0_zz[i] = -2.0 * g_xx_xy_0_zz[i] * a_exp + 4.0 * g_xxzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1092-1098)

    #pragma omp simd aligned(g_xx_xz_0_xx, g_xx_xz_0_xy, g_xx_xz_0_xz, g_xx_xz_0_yy, g_xx_xz_0_yz, g_xx_xz_0_zz, g_xxzz_xz_0_xx, g_xxzz_xz_0_xy, g_xxzz_xz_0_xz, g_xxzz_xz_0_yy, g_xxzz_xz_0_yz, g_xxzz_xz_0_zz, g_zz_0_0_0_xx_xz_0_xx, g_zz_0_0_0_xx_xz_0_xy, g_zz_0_0_0_xx_xz_0_xz, g_zz_0_0_0_xx_xz_0_yy, g_zz_0_0_0_xx_xz_0_yz, g_zz_0_0_0_xx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_xz_0_xx[i] = -2.0 * g_xx_xz_0_xx[i] * a_exp + 4.0 * g_xxzz_xz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xz_0_xy[i] = -2.0 * g_xx_xz_0_xy[i] * a_exp + 4.0 * g_xxzz_xz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xz_0_xz[i] = -2.0 * g_xx_xz_0_xz[i] * a_exp + 4.0 * g_xxzz_xz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xz_0_yy[i] = -2.0 * g_xx_xz_0_yy[i] * a_exp + 4.0 * g_xxzz_xz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xz_0_yz[i] = -2.0 * g_xx_xz_0_yz[i] * a_exp + 4.0 * g_xxzz_xz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_xz_0_zz[i] = -2.0 * g_xx_xz_0_zz[i] * a_exp + 4.0 * g_xxzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1098-1104)

    #pragma omp simd aligned(g_xx_yy_0_xx, g_xx_yy_0_xy, g_xx_yy_0_xz, g_xx_yy_0_yy, g_xx_yy_0_yz, g_xx_yy_0_zz, g_xxzz_yy_0_xx, g_xxzz_yy_0_xy, g_xxzz_yy_0_xz, g_xxzz_yy_0_yy, g_xxzz_yy_0_yz, g_xxzz_yy_0_zz, g_zz_0_0_0_xx_yy_0_xx, g_zz_0_0_0_xx_yy_0_xy, g_zz_0_0_0_xx_yy_0_xz, g_zz_0_0_0_xx_yy_0_yy, g_zz_0_0_0_xx_yy_0_yz, g_zz_0_0_0_xx_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_yy_0_xx[i] = -2.0 * g_xx_yy_0_xx[i] * a_exp + 4.0 * g_xxzz_yy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yy_0_xy[i] = -2.0 * g_xx_yy_0_xy[i] * a_exp + 4.0 * g_xxzz_yy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yy_0_xz[i] = -2.0 * g_xx_yy_0_xz[i] * a_exp + 4.0 * g_xxzz_yy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yy_0_yy[i] = -2.0 * g_xx_yy_0_yy[i] * a_exp + 4.0 * g_xxzz_yy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yy_0_yz[i] = -2.0 * g_xx_yy_0_yz[i] * a_exp + 4.0 * g_xxzz_yy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yy_0_zz[i] = -2.0 * g_xx_yy_0_zz[i] * a_exp + 4.0 * g_xxzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1104-1110)

    #pragma omp simd aligned(g_xx_yz_0_xx, g_xx_yz_0_xy, g_xx_yz_0_xz, g_xx_yz_0_yy, g_xx_yz_0_yz, g_xx_yz_0_zz, g_xxzz_yz_0_xx, g_xxzz_yz_0_xy, g_xxzz_yz_0_xz, g_xxzz_yz_0_yy, g_xxzz_yz_0_yz, g_xxzz_yz_0_zz, g_zz_0_0_0_xx_yz_0_xx, g_zz_0_0_0_xx_yz_0_xy, g_zz_0_0_0_xx_yz_0_xz, g_zz_0_0_0_xx_yz_0_yy, g_zz_0_0_0_xx_yz_0_yz, g_zz_0_0_0_xx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_yz_0_xx[i] = -2.0 * g_xx_yz_0_xx[i] * a_exp + 4.0 * g_xxzz_yz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yz_0_xy[i] = -2.0 * g_xx_yz_0_xy[i] * a_exp + 4.0 * g_xxzz_yz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yz_0_xz[i] = -2.0 * g_xx_yz_0_xz[i] * a_exp + 4.0 * g_xxzz_yz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yz_0_yy[i] = -2.0 * g_xx_yz_0_yy[i] * a_exp + 4.0 * g_xxzz_yz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yz_0_yz[i] = -2.0 * g_xx_yz_0_yz[i] * a_exp + 4.0 * g_xxzz_yz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_yz_0_zz[i] = -2.0 * g_xx_yz_0_zz[i] * a_exp + 4.0 * g_xxzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1110-1116)

    #pragma omp simd aligned(g_xx_zz_0_xx, g_xx_zz_0_xy, g_xx_zz_0_xz, g_xx_zz_0_yy, g_xx_zz_0_yz, g_xx_zz_0_zz, g_xxzz_zz_0_xx, g_xxzz_zz_0_xy, g_xxzz_zz_0_xz, g_xxzz_zz_0_yy, g_xxzz_zz_0_yz, g_xxzz_zz_0_zz, g_zz_0_0_0_xx_zz_0_xx, g_zz_0_0_0_xx_zz_0_xy, g_zz_0_0_0_xx_zz_0_xz, g_zz_0_0_0_xx_zz_0_yy, g_zz_0_0_0_xx_zz_0_yz, g_zz_0_0_0_xx_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xx_zz_0_xx[i] = -2.0 * g_xx_zz_0_xx[i] * a_exp + 4.0 * g_xxzz_zz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_zz_0_xy[i] = -2.0 * g_xx_zz_0_xy[i] * a_exp + 4.0 * g_xxzz_zz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_zz_0_xz[i] = -2.0 * g_xx_zz_0_xz[i] * a_exp + 4.0 * g_xxzz_zz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_zz_0_yy[i] = -2.0 * g_xx_zz_0_yy[i] * a_exp + 4.0 * g_xxzz_zz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_zz_0_yz[i] = -2.0 * g_xx_zz_0_yz[i] * a_exp + 4.0 * g_xxzz_zz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xx_zz_0_zz[i] = -2.0 * g_xx_zz_0_zz[i] * a_exp + 4.0 * g_xxzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1116-1122)

    #pragma omp simd aligned(g_xy_xx_0_xx, g_xy_xx_0_xy, g_xy_xx_0_xz, g_xy_xx_0_yy, g_xy_xx_0_yz, g_xy_xx_0_zz, g_xyzz_xx_0_xx, g_xyzz_xx_0_xy, g_xyzz_xx_0_xz, g_xyzz_xx_0_yy, g_xyzz_xx_0_yz, g_xyzz_xx_0_zz, g_zz_0_0_0_xy_xx_0_xx, g_zz_0_0_0_xy_xx_0_xy, g_zz_0_0_0_xy_xx_0_xz, g_zz_0_0_0_xy_xx_0_yy, g_zz_0_0_0_xy_xx_0_yz, g_zz_0_0_0_xy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_xx_0_xx[i] = -2.0 * g_xy_xx_0_xx[i] * a_exp + 4.0 * g_xyzz_xx_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xx_0_xy[i] = -2.0 * g_xy_xx_0_xy[i] * a_exp + 4.0 * g_xyzz_xx_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xx_0_xz[i] = -2.0 * g_xy_xx_0_xz[i] * a_exp + 4.0 * g_xyzz_xx_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xx_0_yy[i] = -2.0 * g_xy_xx_0_yy[i] * a_exp + 4.0 * g_xyzz_xx_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xx_0_yz[i] = -2.0 * g_xy_xx_0_yz[i] * a_exp + 4.0 * g_xyzz_xx_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xx_0_zz[i] = -2.0 * g_xy_xx_0_zz[i] * a_exp + 4.0 * g_xyzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1122-1128)

    #pragma omp simd aligned(g_xy_xy_0_xx, g_xy_xy_0_xy, g_xy_xy_0_xz, g_xy_xy_0_yy, g_xy_xy_0_yz, g_xy_xy_0_zz, g_xyzz_xy_0_xx, g_xyzz_xy_0_xy, g_xyzz_xy_0_xz, g_xyzz_xy_0_yy, g_xyzz_xy_0_yz, g_xyzz_xy_0_zz, g_zz_0_0_0_xy_xy_0_xx, g_zz_0_0_0_xy_xy_0_xy, g_zz_0_0_0_xy_xy_0_xz, g_zz_0_0_0_xy_xy_0_yy, g_zz_0_0_0_xy_xy_0_yz, g_zz_0_0_0_xy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_xy_0_xx[i] = -2.0 * g_xy_xy_0_xx[i] * a_exp + 4.0 * g_xyzz_xy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xy_0_xy[i] = -2.0 * g_xy_xy_0_xy[i] * a_exp + 4.0 * g_xyzz_xy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xy_0_xz[i] = -2.0 * g_xy_xy_0_xz[i] * a_exp + 4.0 * g_xyzz_xy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xy_0_yy[i] = -2.0 * g_xy_xy_0_yy[i] * a_exp + 4.0 * g_xyzz_xy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xy_0_yz[i] = -2.0 * g_xy_xy_0_yz[i] * a_exp + 4.0 * g_xyzz_xy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xy_0_zz[i] = -2.0 * g_xy_xy_0_zz[i] * a_exp + 4.0 * g_xyzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1128-1134)

    #pragma omp simd aligned(g_xy_xz_0_xx, g_xy_xz_0_xy, g_xy_xz_0_xz, g_xy_xz_0_yy, g_xy_xz_0_yz, g_xy_xz_0_zz, g_xyzz_xz_0_xx, g_xyzz_xz_0_xy, g_xyzz_xz_0_xz, g_xyzz_xz_0_yy, g_xyzz_xz_0_yz, g_xyzz_xz_0_zz, g_zz_0_0_0_xy_xz_0_xx, g_zz_0_0_0_xy_xz_0_xy, g_zz_0_0_0_xy_xz_0_xz, g_zz_0_0_0_xy_xz_0_yy, g_zz_0_0_0_xy_xz_0_yz, g_zz_0_0_0_xy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_xz_0_xx[i] = -2.0 * g_xy_xz_0_xx[i] * a_exp + 4.0 * g_xyzz_xz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xz_0_xy[i] = -2.0 * g_xy_xz_0_xy[i] * a_exp + 4.0 * g_xyzz_xz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xz_0_xz[i] = -2.0 * g_xy_xz_0_xz[i] * a_exp + 4.0 * g_xyzz_xz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xz_0_yy[i] = -2.0 * g_xy_xz_0_yy[i] * a_exp + 4.0 * g_xyzz_xz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xz_0_yz[i] = -2.0 * g_xy_xz_0_yz[i] * a_exp + 4.0 * g_xyzz_xz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_xz_0_zz[i] = -2.0 * g_xy_xz_0_zz[i] * a_exp + 4.0 * g_xyzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1134-1140)

    #pragma omp simd aligned(g_xy_yy_0_xx, g_xy_yy_0_xy, g_xy_yy_0_xz, g_xy_yy_0_yy, g_xy_yy_0_yz, g_xy_yy_0_zz, g_xyzz_yy_0_xx, g_xyzz_yy_0_xy, g_xyzz_yy_0_xz, g_xyzz_yy_0_yy, g_xyzz_yy_0_yz, g_xyzz_yy_0_zz, g_zz_0_0_0_xy_yy_0_xx, g_zz_0_0_0_xy_yy_0_xy, g_zz_0_0_0_xy_yy_0_xz, g_zz_0_0_0_xy_yy_0_yy, g_zz_0_0_0_xy_yy_0_yz, g_zz_0_0_0_xy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_yy_0_xx[i] = -2.0 * g_xy_yy_0_xx[i] * a_exp + 4.0 * g_xyzz_yy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yy_0_xy[i] = -2.0 * g_xy_yy_0_xy[i] * a_exp + 4.0 * g_xyzz_yy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yy_0_xz[i] = -2.0 * g_xy_yy_0_xz[i] * a_exp + 4.0 * g_xyzz_yy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yy_0_yy[i] = -2.0 * g_xy_yy_0_yy[i] * a_exp + 4.0 * g_xyzz_yy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yy_0_yz[i] = -2.0 * g_xy_yy_0_yz[i] * a_exp + 4.0 * g_xyzz_yy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yy_0_zz[i] = -2.0 * g_xy_yy_0_zz[i] * a_exp + 4.0 * g_xyzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1140-1146)

    #pragma omp simd aligned(g_xy_yz_0_xx, g_xy_yz_0_xy, g_xy_yz_0_xz, g_xy_yz_0_yy, g_xy_yz_0_yz, g_xy_yz_0_zz, g_xyzz_yz_0_xx, g_xyzz_yz_0_xy, g_xyzz_yz_0_xz, g_xyzz_yz_0_yy, g_xyzz_yz_0_yz, g_xyzz_yz_0_zz, g_zz_0_0_0_xy_yz_0_xx, g_zz_0_0_0_xy_yz_0_xy, g_zz_0_0_0_xy_yz_0_xz, g_zz_0_0_0_xy_yz_0_yy, g_zz_0_0_0_xy_yz_0_yz, g_zz_0_0_0_xy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_yz_0_xx[i] = -2.0 * g_xy_yz_0_xx[i] * a_exp + 4.0 * g_xyzz_yz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yz_0_xy[i] = -2.0 * g_xy_yz_0_xy[i] * a_exp + 4.0 * g_xyzz_yz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yz_0_xz[i] = -2.0 * g_xy_yz_0_xz[i] * a_exp + 4.0 * g_xyzz_yz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yz_0_yy[i] = -2.0 * g_xy_yz_0_yy[i] * a_exp + 4.0 * g_xyzz_yz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yz_0_yz[i] = -2.0 * g_xy_yz_0_yz[i] * a_exp + 4.0 * g_xyzz_yz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_yz_0_zz[i] = -2.0 * g_xy_yz_0_zz[i] * a_exp + 4.0 * g_xyzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1146-1152)

    #pragma omp simd aligned(g_xy_zz_0_xx, g_xy_zz_0_xy, g_xy_zz_0_xz, g_xy_zz_0_yy, g_xy_zz_0_yz, g_xy_zz_0_zz, g_xyzz_zz_0_xx, g_xyzz_zz_0_xy, g_xyzz_zz_0_xz, g_xyzz_zz_0_yy, g_xyzz_zz_0_yz, g_xyzz_zz_0_zz, g_zz_0_0_0_xy_zz_0_xx, g_zz_0_0_0_xy_zz_0_xy, g_zz_0_0_0_xy_zz_0_xz, g_zz_0_0_0_xy_zz_0_yy, g_zz_0_0_0_xy_zz_0_yz, g_zz_0_0_0_xy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xy_zz_0_xx[i] = -2.0 * g_xy_zz_0_xx[i] * a_exp + 4.0 * g_xyzz_zz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_zz_0_xy[i] = -2.0 * g_xy_zz_0_xy[i] * a_exp + 4.0 * g_xyzz_zz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_zz_0_xz[i] = -2.0 * g_xy_zz_0_xz[i] * a_exp + 4.0 * g_xyzz_zz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_zz_0_yy[i] = -2.0 * g_xy_zz_0_yy[i] * a_exp + 4.0 * g_xyzz_zz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_zz_0_yz[i] = -2.0 * g_xy_zz_0_yz[i] * a_exp + 4.0 * g_xyzz_zz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xy_zz_0_zz[i] = -2.0 * g_xy_zz_0_zz[i] * a_exp + 4.0 * g_xyzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1152-1158)

    #pragma omp simd aligned(g_xz_xx_0_xx, g_xz_xx_0_xy, g_xz_xx_0_xz, g_xz_xx_0_yy, g_xz_xx_0_yz, g_xz_xx_0_zz, g_xzzz_xx_0_xx, g_xzzz_xx_0_xy, g_xzzz_xx_0_xz, g_xzzz_xx_0_yy, g_xzzz_xx_0_yz, g_xzzz_xx_0_zz, g_zz_0_0_0_xz_xx_0_xx, g_zz_0_0_0_xz_xx_0_xy, g_zz_0_0_0_xz_xx_0_xz, g_zz_0_0_0_xz_xx_0_yy, g_zz_0_0_0_xz_xx_0_yz, g_zz_0_0_0_xz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_xx_0_xx[i] = -6.0 * g_xz_xx_0_xx[i] * a_exp + 4.0 * g_xzzz_xx_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xx_0_xy[i] = -6.0 * g_xz_xx_0_xy[i] * a_exp + 4.0 * g_xzzz_xx_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xx_0_xz[i] = -6.0 * g_xz_xx_0_xz[i] * a_exp + 4.0 * g_xzzz_xx_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xx_0_yy[i] = -6.0 * g_xz_xx_0_yy[i] * a_exp + 4.0 * g_xzzz_xx_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xx_0_yz[i] = -6.0 * g_xz_xx_0_yz[i] * a_exp + 4.0 * g_xzzz_xx_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xx_0_zz[i] = -6.0 * g_xz_xx_0_zz[i] * a_exp + 4.0 * g_xzzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1158-1164)

    #pragma omp simd aligned(g_xz_xy_0_xx, g_xz_xy_0_xy, g_xz_xy_0_xz, g_xz_xy_0_yy, g_xz_xy_0_yz, g_xz_xy_0_zz, g_xzzz_xy_0_xx, g_xzzz_xy_0_xy, g_xzzz_xy_0_xz, g_xzzz_xy_0_yy, g_xzzz_xy_0_yz, g_xzzz_xy_0_zz, g_zz_0_0_0_xz_xy_0_xx, g_zz_0_0_0_xz_xy_0_xy, g_zz_0_0_0_xz_xy_0_xz, g_zz_0_0_0_xz_xy_0_yy, g_zz_0_0_0_xz_xy_0_yz, g_zz_0_0_0_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_xy_0_xx[i] = -6.0 * g_xz_xy_0_xx[i] * a_exp + 4.0 * g_xzzz_xy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xy_0_xy[i] = -6.0 * g_xz_xy_0_xy[i] * a_exp + 4.0 * g_xzzz_xy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xy_0_xz[i] = -6.0 * g_xz_xy_0_xz[i] * a_exp + 4.0 * g_xzzz_xy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xy_0_yy[i] = -6.0 * g_xz_xy_0_yy[i] * a_exp + 4.0 * g_xzzz_xy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xy_0_yz[i] = -6.0 * g_xz_xy_0_yz[i] * a_exp + 4.0 * g_xzzz_xy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xy_0_zz[i] = -6.0 * g_xz_xy_0_zz[i] * a_exp + 4.0 * g_xzzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1164-1170)

    #pragma omp simd aligned(g_xz_xz_0_xx, g_xz_xz_0_xy, g_xz_xz_0_xz, g_xz_xz_0_yy, g_xz_xz_0_yz, g_xz_xz_0_zz, g_xzzz_xz_0_xx, g_xzzz_xz_0_xy, g_xzzz_xz_0_xz, g_xzzz_xz_0_yy, g_xzzz_xz_0_yz, g_xzzz_xz_0_zz, g_zz_0_0_0_xz_xz_0_xx, g_zz_0_0_0_xz_xz_0_xy, g_zz_0_0_0_xz_xz_0_xz, g_zz_0_0_0_xz_xz_0_yy, g_zz_0_0_0_xz_xz_0_yz, g_zz_0_0_0_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_xz_0_xx[i] = -6.0 * g_xz_xz_0_xx[i] * a_exp + 4.0 * g_xzzz_xz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xz_0_xy[i] = -6.0 * g_xz_xz_0_xy[i] * a_exp + 4.0 * g_xzzz_xz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xz_0_xz[i] = -6.0 * g_xz_xz_0_xz[i] * a_exp + 4.0 * g_xzzz_xz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xz_0_yy[i] = -6.0 * g_xz_xz_0_yy[i] * a_exp + 4.0 * g_xzzz_xz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xz_0_yz[i] = -6.0 * g_xz_xz_0_yz[i] * a_exp + 4.0 * g_xzzz_xz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_xz_0_zz[i] = -6.0 * g_xz_xz_0_zz[i] * a_exp + 4.0 * g_xzzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1170-1176)

    #pragma omp simd aligned(g_xz_yy_0_xx, g_xz_yy_0_xy, g_xz_yy_0_xz, g_xz_yy_0_yy, g_xz_yy_0_yz, g_xz_yy_0_zz, g_xzzz_yy_0_xx, g_xzzz_yy_0_xy, g_xzzz_yy_0_xz, g_xzzz_yy_0_yy, g_xzzz_yy_0_yz, g_xzzz_yy_0_zz, g_zz_0_0_0_xz_yy_0_xx, g_zz_0_0_0_xz_yy_0_xy, g_zz_0_0_0_xz_yy_0_xz, g_zz_0_0_0_xz_yy_0_yy, g_zz_0_0_0_xz_yy_0_yz, g_zz_0_0_0_xz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_yy_0_xx[i] = -6.0 * g_xz_yy_0_xx[i] * a_exp + 4.0 * g_xzzz_yy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yy_0_xy[i] = -6.0 * g_xz_yy_0_xy[i] * a_exp + 4.0 * g_xzzz_yy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yy_0_xz[i] = -6.0 * g_xz_yy_0_xz[i] * a_exp + 4.0 * g_xzzz_yy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yy_0_yy[i] = -6.0 * g_xz_yy_0_yy[i] * a_exp + 4.0 * g_xzzz_yy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yy_0_yz[i] = -6.0 * g_xz_yy_0_yz[i] * a_exp + 4.0 * g_xzzz_yy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yy_0_zz[i] = -6.0 * g_xz_yy_0_zz[i] * a_exp + 4.0 * g_xzzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1176-1182)

    #pragma omp simd aligned(g_xz_yz_0_xx, g_xz_yz_0_xy, g_xz_yz_0_xz, g_xz_yz_0_yy, g_xz_yz_0_yz, g_xz_yz_0_zz, g_xzzz_yz_0_xx, g_xzzz_yz_0_xy, g_xzzz_yz_0_xz, g_xzzz_yz_0_yy, g_xzzz_yz_0_yz, g_xzzz_yz_0_zz, g_zz_0_0_0_xz_yz_0_xx, g_zz_0_0_0_xz_yz_0_xy, g_zz_0_0_0_xz_yz_0_xz, g_zz_0_0_0_xz_yz_0_yy, g_zz_0_0_0_xz_yz_0_yz, g_zz_0_0_0_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_yz_0_xx[i] = -6.0 * g_xz_yz_0_xx[i] * a_exp + 4.0 * g_xzzz_yz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yz_0_xy[i] = -6.0 * g_xz_yz_0_xy[i] * a_exp + 4.0 * g_xzzz_yz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yz_0_xz[i] = -6.0 * g_xz_yz_0_xz[i] * a_exp + 4.0 * g_xzzz_yz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yz_0_yy[i] = -6.0 * g_xz_yz_0_yy[i] * a_exp + 4.0 * g_xzzz_yz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yz_0_yz[i] = -6.0 * g_xz_yz_0_yz[i] * a_exp + 4.0 * g_xzzz_yz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_yz_0_zz[i] = -6.0 * g_xz_yz_0_zz[i] * a_exp + 4.0 * g_xzzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1182-1188)

    #pragma omp simd aligned(g_xz_zz_0_xx, g_xz_zz_0_xy, g_xz_zz_0_xz, g_xz_zz_0_yy, g_xz_zz_0_yz, g_xz_zz_0_zz, g_xzzz_zz_0_xx, g_xzzz_zz_0_xy, g_xzzz_zz_0_xz, g_xzzz_zz_0_yy, g_xzzz_zz_0_yz, g_xzzz_zz_0_zz, g_zz_0_0_0_xz_zz_0_xx, g_zz_0_0_0_xz_zz_0_xy, g_zz_0_0_0_xz_zz_0_xz, g_zz_0_0_0_xz_zz_0_yy, g_zz_0_0_0_xz_zz_0_yz, g_zz_0_0_0_xz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_xz_zz_0_xx[i] = -6.0 * g_xz_zz_0_xx[i] * a_exp + 4.0 * g_xzzz_zz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_zz_0_xy[i] = -6.0 * g_xz_zz_0_xy[i] * a_exp + 4.0 * g_xzzz_zz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_zz_0_xz[i] = -6.0 * g_xz_zz_0_xz[i] * a_exp + 4.0 * g_xzzz_zz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_zz_0_yy[i] = -6.0 * g_xz_zz_0_yy[i] * a_exp + 4.0 * g_xzzz_zz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_zz_0_yz[i] = -6.0 * g_xz_zz_0_yz[i] * a_exp + 4.0 * g_xzzz_zz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_xz_zz_0_zz[i] = -6.0 * g_xz_zz_0_zz[i] * a_exp + 4.0 * g_xzzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1188-1194)

    #pragma omp simd aligned(g_yy_xx_0_xx, g_yy_xx_0_xy, g_yy_xx_0_xz, g_yy_xx_0_yy, g_yy_xx_0_yz, g_yy_xx_0_zz, g_yyzz_xx_0_xx, g_yyzz_xx_0_xy, g_yyzz_xx_0_xz, g_yyzz_xx_0_yy, g_yyzz_xx_0_yz, g_yyzz_xx_0_zz, g_zz_0_0_0_yy_xx_0_xx, g_zz_0_0_0_yy_xx_0_xy, g_zz_0_0_0_yy_xx_0_xz, g_zz_0_0_0_yy_xx_0_yy, g_zz_0_0_0_yy_xx_0_yz, g_zz_0_0_0_yy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_xx_0_xx[i] = -2.0 * g_yy_xx_0_xx[i] * a_exp + 4.0 * g_yyzz_xx_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xx_0_xy[i] = -2.0 * g_yy_xx_0_xy[i] * a_exp + 4.0 * g_yyzz_xx_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xx_0_xz[i] = -2.0 * g_yy_xx_0_xz[i] * a_exp + 4.0 * g_yyzz_xx_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xx_0_yy[i] = -2.0 * g_yy_xx_0_yy[i] * a_exp + 4.0 * g_yyzz_xx_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xx_0_yz[i] = -2.0 * g_yy_xx_0_yz[i] * a_exp + 4.0 * g_yyzz_xx_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xx_0_zz[i] = -2.0 * g_yy_xx_0_zz[i] * a_exp + 4.0 * g_yyzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1194-1200)

    #pragma omp simd aligned(g_yy_xy_0_xx, g_yy_xy_0_xy, g_yy_xy_0_xz, g_yy_xy_0_yy, g_yy_xy_0_yz, g_yy_xy_0_zz, g_yyzz_xy_0_xx, g_yyzz_xy_0_xy, g_yyzz_xy_0_xz, g_yyzz_xy_0_yy, g_yyzz_xy_0_yz, g_yyzz_xy_0_zz, g_zz_0_0_0_yy_xy_0_xx, g_zz_0_0_0_yy_xy_0_xy, g_zz_0_0_0_yy_xy_0_xz, g_zz_0_0_0_yy_xy_0_yy, g_zz_0_0_0_yy_xy_0_yz, g_zz_0_0_0_yy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_xy_0_xx[i] = -2.0 * g_yy_xy_0_xx[i] * a_exp + 4.0 * g_yyzz_xy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xy_0_xy[i] = -2.0 * g_yy_xy_0_xy[i] * a_exp + 4.0 * g_yyzz_xy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xy_0_xz[i] = -2.0 * g_yy_xy_0_xz[i] * a_exp + 4.0 * g_yyzz_xy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xy_0_yy[i] = -2.0 * g_yy_xy_0_yy[i] * a_exp + 4.0 * g_yyzz_xy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xy_0_yz[i] = -2.0 * g_yy_xy_0_yz[i] * a_exp + 4.0 * g_yyzz_xy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xy_0_zz[i] = -2.0 * g_yy_xy_0_zz[i] * a_exp + 4.0 * g_yyzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1200-1206)

    #pragma omp simd aligned(g_yy_xz_0_xx, g_yy_xz_0_xy, g_yy_xz_0_xz, g_yy_xz_0_yy, g_yy_xz_0_yz, g_yy_xz_0_zz, g_yyzz_xz_0_xx, g_yyzz_xz_0_xy, g_yyzz_xz_0_xz, g_yyzz_xz_0_yy, g_yyzz_xz_0_yz, g_yyzz_xz_0_zz, g_zz_0_0_0_yy_xz_0_xx, g_zz_0_0_0_yy_xz_0_xy, g_zz_0_0_0_yy_xz_0_xz, g_zz_0_0_0_yy_xz_0_yy, g_zz_0_0_0_yy_xz_0_yz, g_zz_0_0_0_yy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_xz_0_xx[i] = -2.0 * g_yy_xz_0_xx[i] * a_exp + 4.0 * g_yyzz_xz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xz_0_xy[i] = -2.0 * g_yy_xz_0_xy[i] * a_exp + 4.0 * g_yyzz_xz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xz_0_xz[i] = -2.0 * g_yy_xz_0_xz[i] * a_exp + 4.0 * g_yyzz_xz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xz_0_yy[i] = -2.0 * g_yy_xz_0_yy[i] * a_exp + 4.0 * g_yyzz_xz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xz_0_yz[i] = -2.0 * g_yy_xz_0_yz[i] * a_exp + 4.0 * g_yyzz_xz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_xz_0_zz[i] = -2.0 * g_yy_xz_0_zz[i] * a_exp + 4.0 * g_yyzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1206-1212)

    #pragma omp simd aligned(g_yy_yy_0_xx, g_yy_yy_0_xy, g_yy_yy_0_xz, g_yy_yy_0_yy, g_yy_yy_0_yz, g_yy_yy_0_zz, g_yyzz_yy_0_xx, g_yyzz_yy_0_xy, g_yyzz_yy_0_xz, g_yyzz_yy_0_yy, g_yyzz_yy_0_yz, g_yyzz_yy_0_zz, g_zz_0_0_0_yy_yy_0_xx, g_zz_0_0_0_yy_yy_0_xy, g_zz_0_0_0_yy_yy_0_xz, g_zz_0_0_0_yy_yy_0_yy, g_zz_0_0_0_yy_yy_0_yz, g_zz_0_0_0_yy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_yy_0_xx[i] = -2.0 * g_yy_yy_0_xx[i] * a_exp + 4.0 * g_yyzz_yy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yy_0_xy[i] = -2.0 * g_yy_yy_0_xy[i] * a_exp + 4.0 * g_yyzz_yy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yy_0_xz[i] = -2.0 * g_yy_yy_0_xz[i] * a_exp + 4.0 * g_yyzz_yy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yy_0_yy[i] = -2.0 * g_yy_yy_0_yy[i] * a_exp + 4.0 * g_yyzz_yy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yy_0_yz[i] = -2.0 * g_yy_yy_0_yz[i] * a_exp + 4.0 * g_yyzz_yy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yy_0_zz[i] = -2.0 * g_yy_yy_0_zz[i] * a_exp + 4.0 * g_yyzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1212-1218)

    #pragma omp simd aligned(g_yy_yz_0_xx, g_yy_yz_0_xy, g_yy_yz_0_xz, g_yy_yz_0_yy, g_yy_yz_0_yz, g_yy_yz_0_zz, g_yyzz_yz_0_xx, g_yyzz_yz_0_xy, g_yyzz_yz_0_xz, g_yyzz_yz_0_yy, g_yyzz_yz_0_yz, g_yyzz_yz_0_zz, g_zz_0_0_0_yy_yz_0_xx, g_zz_0_0_0_yy_yz_0_xy, g_zz_0_0_0_yy_yz_0_xz, g_zz_0_0_0_yy_yz_0_yy, g_zz_0_0_0_yy_yz_0_yz, g_zz_0_0_0_yy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_yz_0_xx[i] = -2.0 * g_yy_yz_0_xx[i] * a_exp + 4.0 * g_yyzz_yz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yz_0_xy[i] = -2.0 * g_yy_yz_0_xy[i] * a_exp + 4.0 * g_yyzz_yz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yz_0_xz[i] = -2.0 * g_yy_yz_0_xz[i] * a_exp + 4.0 * g_yyzz_yz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yz_0_yy[i] = -2.0 * g_yy_yz_0_yy[i] * a_exp + 4.0 * g_yyzz_yz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yz_0_yz[i] = -2.0 * g_yy_yz_0_yz[i] * a_exp + 4.0 * g_yyzz_yz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_yz_0_zz[i] = -2.0 * g_yy_yz_0_zz[i] * a_exp + 4.0 * g_yyzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1218-1224)

    #pragma omp simd aligned(g_yy_zz_0_xx, g_yy_zz_0_xy, g_yy_zz_0_xz, g_yy_zz_0_yy, g_yy_zz_0_yz, g_yy_zz_0_zz, g_yyzz_zz_0_xx, g_yyzz_zz_0_xy, g_yyzz_zz_0_xz, g_yyzz_zz_0_yy, g_yyzz_zz_0_yz, g_yyzz_zz_0_zz, g_zz_0_0_0_yy_zz_0_xx, g_zz_0_0_0_yy_zz_0_xy, g_zz_0_0_0_yy_zz_0_xz, g_zz_0_0_0_yy_zz_0_yy, g_zz_0_0_0_yy_zz_0_yz, g_zz_0_0_0_yy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yy_zz_0_xx[i] = -2.0 * g_yy_zz_0_xx[i] * a_exp + 4.0 * g_yyzz_zz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_zz_0_xy[i] = -2.0 * g_yy_zz_0_xy[i] * a_exp + 4.0 * g_yyzz_zz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_zz_0_xz[i] = -2.0 * g_yy_zz_0_xz[i] * a_exp + 4.0 * g_yyzz_zz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_zz_0_yy[i] = -2.0 * g_yy_zz_0_yy[i] * a_exp + 4.0 * g_yyzz_zz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_zz_0_yz[i] = -2.0 * g_yy_zz_0_yz[i] * a_exp + 4.0 * g_yyzz_zz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_yy_zz_0_zz[i] = -2.0 * g_yy_zz_0_zz[i] * a_exp + 4.0 * g_yyzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1224-1230)

    #pragma omp simd aligned(g_yz_xx_0_xx, g_yz_xx_0_xy, g_yz_xx_0_xz, g_yz_xx_0_yy, g_yz_xx_0_yz, g_yz_xx_0_zz, g_yzzz_xx_0_xx, g_yzzz_xx_0_xy, g_yzzz_xx_0_xz, g_yzzz_xx_0_yy, g_yzzz_xx_0_yz, g_yzzz_xx_0_zz, g_zz_0_0_0_yz_xx_0_xx, g_zz_0_0_0_yz_xx_0_xy, g_zz_0_0_0_yz_xx_0_xz, g_zz_0_0_0_yz_xx_0_yy, g_zz_0_0_0_yz_xx_0_yz, g_zz_0_0_0_yz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_xx_0_xx[i] = -6.0 * g_yz_xx_0_xx[i] * a_exp + 4.0 * g_yzzz_xx_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xx_0_xy[i] = -6.0 * g_yz_xx_0_xy[i] * a_exp + 4.0 * g_yzzz_xx_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xx_0_xz[i] = -6.0 * g_yz_xx_0_xz[i] * a_exp + 4.0 * g_yzzz_xx_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xx_0_yy[i] = -6.0 * g_yz_xx_0_yy[i] * a_exp + 4.0 * g_yzzz_xx_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xx_0_yz[i] = -6.0 * g_yz_xx_0_yz[i] * a_exp + 4.0 * g_yzzz_xx_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xx_0_zz[i] = -6.0 * g_yz_xx_0_zz[i] * a_exp + 4.0 * g_yzzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1230-1236)

    #pragma omp simd aligned(g_yz_xy_0_xx, g_yz_xy_0_xy, g_yz_xy_0_xz, g_yz_xy_0_yy, g_yz_xy_0_yz, g_yz_xy_0_zz, g_yzzz_xy_0_xx, g_yzzz_xy_0_xy, g_yzzz_xy_0_xz, g_yzzz_xy_0_yy, g_yzzz_xy_0_yz, g_yzzz_xy_0_zz, g_zz_0_0_0_yz_xy_0_xx, g_zz_0_0_0_yz_xy_0_xy, g_zz_0_0_0_yz_xy_0_xz, g_zz_0_0_0_yz_xy_0_yy, g_zz_0_0_0_yz_xy_0_yz, g_zz_0_0_0_yz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_xy_0_xx[i] = -6.0 * g_yz_xy_0_xx[i] * a_exp + 4.0 * g_yzzz_xy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xy_0_xy[i] = -6.0 * g_yz_xy_0_xy[i] * a_exp + 4.0 * g_yzzz_xy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xy_0_xz[i] = -6.0 * g_yz_xy_0_xz[i] * a_exp + 4.0 * g_yzzz_xy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xy_0_yy[i] = -6.0 * g_yz_xy_0_yy[i] * a_exp + 4.0 * g_yzzz_xy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xy_0_yz[i] = -6.0 * g_yz_xy_0_yz[i] * a_exp + 4.0 * g_yzzz_xy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xy_0_zz[i] = -6.0 * g_yz_xy_0_zz[i] * a_exp + 4.0 * g_yzzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1236-1242)

    #pragma omp simd aligned(g_yz_xz_0_xx, g_yz_xz_0_xy, g_yz_xz_0_xz, g_yz_xz_0_yy, g_yz_xz_0_yz, g_yz_xz_0_zz, g_yzzz_xz_0_xx, g_yzzz_xz_0_xy, g_yzzz_xz_0_xz, g_yzzz_xz_0_yy, g_yzzz_xz_0_yz, g_yzzz_xz_0_zz, g_zz_0_0_0_yz_xz_0_xx, g_zz_0_0_0_yz_xz_0_xy, g_zz_0_0_0_yz_xz_0_xz, g_zz_0_0_0_yz_xz_0_yy, g_zz_0_0_0_yz_xz_0_yz, g_zz_0_0_0_yz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_xz_0_xx[i] = -6.0 * g_yz_xz_0_xx[i] * a_exp + 4.0 * g_yzzz_xz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xz_0_xy[i] = -6.0 * g_yz_xz_0_xy[i] * a_exp + 4.0 * g_yzzz_xz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xz_0_xz[i] = -6.0 * g_yz_xz_0_xz[i] * a_exp + 4.0 * g_yzzz_xz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xz_0_yy[i] = -6.0 * g_yz_xz_0_yy[i] * a_exp + 4.0 * g_yzzz_xz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xz_0_yz[i] = -6.0 * g_yz_xz_0_yz[i] * a_exp + 4.0 * g_yzzz_xz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_xz_0_zz[i] = -6.0 * g_yz_xz_0_zz[i] * a_exp + 4.0 * g_yzzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1242-1248)

    #pragma omp simd aligned(g_yz_yy_0_xx, g_yz_yy_0_xy, g_yz_yy_0_xz, g_yz_yy_0_yy, g_yz_yy_0_yz, g_yz_yy_0_zz, g_yzzz_yy_0_xx, g_yzzz_yy_0_xy, g_yzzz_yy_0_xz, g_yzzz_yy_0_yy, g_yzzz_yy_0_yz, g_yzzz_yy_0_zz, g_zz_0_0_0_yz_yy_0_xx, g_zz_0_0_0_yz_yy_0_xy, g_zz_0_0_0_yz_yy_0_xz, g_zz_0_0_0_yz_yy_0_yy, g_zz_0_0_0_yz_yy_0_yz, g_zz_0_0_0_yz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_yy_0_xx[i] = -6.0 * g_yz_yy_0_xx[i] * a_exp + 4.0 * g_yzzz_yy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yy_0_xy[i] = -6.0 * g_yz_yy_0_xy[i] * a_exp + 4.0 * g_yzzz_yy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yy_0_xz[i] = -6.0 * g_yz_yy_0_xz[i] * a_exp + 4.0 * g_yzzz_yy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yy_0_yy[i] = -6.0 * g_yz_yy_0_yy[i] * a_exp + 4.0 * g_yzzz_yy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yy_0_yz[i] = -6.0 * g_yz_yy_0_yz[i] * a_exp + 4.0 * g_yzzz_yy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yy_0_zz[i] = -6.0 * g_yz_yy_0_zz[i] * a_exp + 4.0 * g_yzzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1248-1254)

    #pragma omp simd aligned(g_yz_yz_0_xx, g_yz_yz_0_xy, g_yz_yz_0_xz, g_yz_yz_0_yy, g_yz_yz_0_yz, g_yz_yz_0_zz, g_yzzz_yz_0_xx, g_yzzz_yz_0_xy, g_yzzz_yz_0_xz, g_yzzz_yz_0_yy, g_yzzz_yz_0_yz, g_yzzz_yz_0_zz, g_zz_0_0_0_yz_yz_0_xx, g_zz_0_0_0_yz_yz_0_xy, g_zz_0_0_0_yz_yz_0_xz, g_zz_0_0_0_yz_yz_0_yy, g_zz_0_0_0_yz_yz_0_yz, g_zz_0_0_0_yz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_yz_0_xx[i] = -6.0 * g_yz_yz_0_xx[i] * a_exp + 4.0 * g_yzzz_yz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yz_0_xy[i] = -6.0 * g_yz_yz_0_xy[i] * a_exp + 4.0 * g_yzzz_yz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yz_0_xz[i] = -6.0 * g_yz_yz_0_xz[i] * a_exp + 4.0 * g_yzzz_yz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yz_0_yy[i] = -6.0 * g_yz_yz_0_yy[i] * a_exp + 4.0 * g_yzzz_yz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yz_0_yz[i] = -6.0 * g_yz_yz_0_yz[i] * a_exp + 4.0 * g_yzzz_yz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_yz_0_zz[i] = -6.0 * g_yz_yz_0_zz[i] * a_exp + 4.0 * g_yzzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1254-1260)

    #pragma omp simd aligned(g_yz_zz_0_xx, g_yz_zz_0_xy, g_yz_zz_0_xz, g_yz_zz_0_yy, g_yz_zz_0_yz, g_yz_zz_0_zz, g_yzzz_zz_0_xx, g_yzzz_zz_0_xy, g_yzzz_zz_0_xz, g_yzzz_zz_0_yy, g_yzzz_zz_0_yz, g_yzzz_zz_0_zz, g_zz_0_0_0_yz_zz_0_xx, g_zz_0_0_0_yz_zz_0_xy, g_zz_0_0_0_yz_zz_0_xz, g_zz_0_0_0_yz_zz_0_yy, g_zz_0_0_0_yz_zz_0_yz, g_zz_0_0_0_yz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_yz_zz_0_xx[i] = -6.0 * g_yz_zz_0_xx[i] * a_exp + 4.0 * g_yzzz_zz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_zz_0_xy[i] = -6.0 * g_yz_zz_0_xy[i] * a_exp + 4.0 * g_yzzz_zz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_zz_0_xz[i] = -6.0 * g_yz_zz_0_xz[i] * a_exp + 4.0 * g_yzzz_zz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_zz_0_yy[i] = -6.0 * g_yz_zz_0_yy[i] * a_exp + 4.0 * g_yzzz_zz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_zz_0_yz[i] = -6.0 * g_yz_zz_0_yz[i] * a_exp + 4.0 * g_yzzz_zz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_yz_zz_0_zz[i] = -6.0 * g_yz_zz_0_zz[i] * a_exp + 4.0 * g_yzzz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1260-1266)

    #pragma omp simd aligned(g_0_xx_0_xx, g_0_xx_0_xy, g_0_xx_0_xz, g_0_xx_0_yy, g_0_xx_0_yz, g_0_xx_0_zz, g_zz_0_0_0_zz_xx_0_xx, g_zz_0_0_0_zz_xx_0_xy, g_zz_0_0_0_zz_xx_0_xz, g_zz_0_0_0_zz_xx_0_yy, g_zz_0_0_0_zz_xx_0_yz, g_zz_0_0_0_zz_xx_0_zz, g_zz_xx_0_xx, g_zz_xx_0_xy, g_zz_xx_0_xz, g_zz_xx_0_yy, g_zz_xx_0_yz, g_zz_xx_0_zz, g_zzzz_xx_0_xx, g_zzzz_xx_0_xy, g_zzzz_xx_0_xz, g_zzzz_xx_0_yy, g_zzzz_xx_0_yz, g_zzzz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_xx_0_xx[i] = 2.0 * g_0_xx_0_xx[i] - 10.0 * g_zz_xx_0_xx[i] * a_exp + 4.0 * g_zzzz_xx_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xx_0_xy[i] = 2.0 * g_0_xx_0_xy[i] - 10.0 * g_zz_xx_0_xy[i] * a_exp + 4.0 * g_zzzz_xx_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xx_0_xz[i] = 2.0 * g_0_xx_0_xz[i] - 10.0 * g_zz_xx_0_xz[i] * a_exp + 4.0 * g_zzzz_xx_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xx_0_yy[i] = 2.0 * g_0_xx_0_yy[i] - 10.0 * g_zz_xx_0_yy[i] * a_exp + 4.0 * g_zzzz_xx_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xx_0_yz[i] = 2.0 * g_0_xx_0_yz[i] - 10.0 * g_zz_xx_0_yz[i] * a_exp + 4.0 * g_zzzz_xx_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xx_0_zz[i] = 2.0 * g_0_xx_0_zz[i] - 10.0 * g_zz_xx_0_zz[i] * a_exp + 4.0 * g_zzzz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1266-1272)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_zz_0_0_0_zz_xy_0_xx, g_zz_0_0_0_zz_xy_0_xy, g_zz_0_0_0_zz_xy_0_xz, g_zz_0_0_0_zz_xy_0_yy, g_zz_0_0_0_zz_xy_0_yz, g_zz_0_0_0_zz_xy_0_zz, g_zz_xy_0_xx, g_zz_xy_0_xy, g_zz_xy_0_xz, g_zz_xy_0_yy, g_zz_xy_0_yz, g_zz_xy_0_zz, g_zzzz_xy_0_xx, g_zzzz_xy_0_xy, g_zzzz_xy_0_xz, g_zzzz_xy_0_yy, g_zzzz_xy_0_yz, g_zzzz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_xy_0_xx[i] = 2.0 * g_0_xy_0_xx[i] - 10.0 * g_zz_xy_0_xx[i] * a_exp + 4.0 * g_zzzz_xy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xy_0_xy[i] = 2.0 * g_0_xy_0_xy[i] - 10.0 * g_zz_xy_0_xy[i] * a_exp + 4.0 * g_zzzz_xy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xy_0_xz[i] = 2.0 * g_0_xy_0_xz[i] - 10.0 * g_zz_xy_0_xz[i] * a_exp + 4.0 * g_zzzz_xy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xy_0_yy[i] = 2.0 * g_0_xy_0_yy[i] - 10.0 * g_zz_xy_0_yy[i] * a_exp + 4.0 * g_zzzz_xy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xy_0_yz[i] = 2.0 * g_0_xy_0_yz[i] - 10.0 * g_zz_xy_0_yz[i] * a_exp + 4.0 * g_zzzz_xy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xy_0_zz[i] = 2.0 * g_0_xy_0_zz[i] - 10.0 * g_zz_xy_0_zz[i] * a_exp + 4.0 * g_zzzz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1272-1278)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_zz_0_0_0_zz_xz_0_xx, g_zz_0_0_0_zz_xz_0_xy, g_zz_0_0_0_zz_xz_0_xz, g_zz_0_0_0_zz_xz_0_yy, g_zz_0_0_0_zz_xz_0_yz, g_zz_0_0_0_zz_xz_0_zz, g_zz_xz_0_xx, g_zz_xz_0_xy, g_zz_xz_0_xz, g_zz_xz_0_yy, g_zz_xz_0_yz, g_zz_xz_0_zz, g_zzzz_xz_0_xx, g_zzzz_xz_0_xy, g_zzzz_xz_0_xz, g_zzzz_xz_0_yy, g_zzzz_xz_0_yz, g_zzzz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_xz_0_xx[i] = 2.0 * g_0_xz_0_xx[i] - 10.0 * g_zz_xz_0_xx[i] * a_exp + 4.0 * g_zzzz_xz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xz_0_xy[i] = 2.0 * g_0_xz_0_xy[i] - 10.0 * g_zz_xz_0_xy[i] * a_exp + 4.0 * g_zzzz_xz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xz_0_xz[i] = 2.0 * g_0_xz_0_xz[i] - 10.0 * g_zz_xz_0_xz[i] * a_exp + 4.0 * g_zzzz_xz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xz_0_yy[i] = 2.0 * g_0_xz_0_yy[i] - 10.0 * g_zz_xz_0_yy[i] * a_exp + 4.0 * g_zzzz_xz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xz_0_yz[i] = 2.0 * g_0_xz_0_yz[i] - 10.0 * g_zz_xz_0_yz[i] * a_exp + 4.0 * g_zzzz_xz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_xz_0_zz[i] = 2.0 * g_0_xz_0_zz[i] - 10.0 * g_zz_xz_0_zz[i] * a_exp + 4.0 * g_zzzz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1278-1284)

    #pragma omp simd aligned(g_0_yy_0_xx, g_0_yy_0_xy, g_0_yy_0_xz, g_0_yy_0_yy, g_0_yy_0_yz, g_0_yy_0_zz, g_zz_0_0_0_zz_yy_0_xx, g_zz_0_0_0_zz_yy_0_xy, g_zz_0_0_0_zz_yy_0_xz, g_zz_0_0_0_zz_yy_0_yy, g_zz_0_0_0_zz_yy_0_yz, g_zz_0_0_0_zz_yy_0_zz, g_zz_yy_0_xx, g_zz_yy_0_xy, g_zz_yy_0_xz, g_zz_yy_0_yy, g_zz_yy_0_yz, g_zz_yy_0_zz, g_zzzz_yy_0_xx, g_zzzz_yy_0_xy, g_zzzz_yy_0_xz, g_zzzz_yy_0_yy, g_zzzz_yy_0_yz, g_zzzz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_yy_0_xx[i] = 2.0 * g_0_yy_0_xx[i] - 10.0 * g_zz_yy_0_xx[i] * a_exp + 4.0 * g_zzzz_yy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yy_0_xy[i] = 2.0 * g_0_yy_0_xy[i] - 10.0 * g_zz_yy_0_xy[i] * a_exp + 4.0 * g_zzzz_yy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yy_0_xz[i] = 2.0 * g_0_yy_0_xz[i] - 10.0 * g_zz_yy_0_xz[i] * a_exp + 4.0 * g_zzzz_yy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yy_0_yy[i] = 2.0 * g_0_yy_0_yy[i] - 10.0 * g_zz_yy_0_yy[i] * a_exp + 4.0 * g_zzzz_yy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yy_0_yz[i] = 2.0 * g_0_yy_0_yz[i] - 10.0 * g_zz_yy_0_yz[i] * a_exp + 4.0 * g_zzzz_yy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yy_0_zz[i] = 2.0 * g_0_yy_0_zz[i] - 10.0 * g_zz_yy_0_zz[i] * a_exp + 4.0 * g_zzzz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1284-1290)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_zz_0_0_0_zz_yz_0_xx, g_zz_0_0_0_zz_yz_0_xy, g_zz_0_0_0_zz_yz_0_xz, g_zz_0_0_0_zz_yz_0_yy, g_zz_0_0_0_zz_yz_0_yz, g_zz_0_0_0_zz_yz_0_zz, g_zz_yz_0_xx, g_zz_yz_0_xy, g_zz_yz_0_xz, g_zz_yz_0_yy, g_zz_yz_0_yz, g_zz_yz_0_zz, g_zzzz_yz_0_xx, g_zzzz_yz_0_xy, g_zzzz_yz_0_xz, g_zzzz_yz_0_yy, g_zzzz_yz_0_yz, g_zzzz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_yz_0_xx[i] = 2.0 * g_0_yz_0_xx[i] - 10.0 * g_zz_yz_0_xx[i] * a_exp + 4.0 * g_zzzz_yz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yz_0_xy[i] = 2.0 * g_0_yz_0_xy[i] - 10.0 * g_zz_yz_0_xy[i] * a_exp + 4.0 * g_zzzz_yz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yz_0_xz[i] = 2.0 * g_0_yz_0_xz[i] - 10.0 * g_zz_yz_0_xz[i] * a_exp + 4.0 * g_zzzz_yz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yz_0_yy[i] = 2.0 * g_0_yz_0_yy[i] - 10.0 * g_zz_yz_0_yy[i] * a_exp + 4.0 * g_zzzz_yz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yz_0_yz[i] = 2.0 * g_0_yz_0_yz[i] - 10.0 * g_zz_yz_0_yz[i] * a_exp + 4.0 * g_zzzz_yz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_yz_0_zz[i] = 2.0 * g_0_yz_0_zz[i] - 10.0 * g_zz_yz_0_zz[i] * a_exp + 4.0 * g_zzzz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (1290-1296)

    #pragma omp simd aligned(g_0_zz_0_xx, g_0_zz_0_xy, g_0_zz_0_xz, g_0_zz_0_yy, g_0_zz_0_yz, g_0_zz_0_zz, g_zz_0_0_0_zz_zz_0_xx, g_zz_0_0_0_zz_zz_0_xy, g_zz_0_0_0_zz_zz_0_xz, g_zz_0_0_0_zz_zz_0_yy, g_zz_0_0_0_zz_zz_0_yz, g_zz_0_0_0_zz_zz_0_zz, g_zz_zz_0_xx, g_zz_zz_0_xy, g_zz_zz_0_xz, g_zz_zz_0_yy, g_zz_zz_0_yz, g_zz_zz_0_zz, g_zzzz_zz_0_xx, g_zzzz_zz_0_xy, g_zzzz_zz_0_xz, g_zzzz_zz_0_yy, g_zzzz_zz_0_yz, g_zzzz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_zz_zz_0_xx[i] = 2.0 * g_0_zz_0_xx[i] - 10.0 * g_zz_zz_0_xx[i] * a_exp + 4.0 * g_zzzz_zz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_zz_0_xy[i] = 2.0 * g_0_zz_0_xy[i] - 10.0 * g_zz_zz_0_xy[i] * a_exp + 4.0 * g_zzzz_zz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_zz_0_xz[i] = 2.0 * g_0_zz_0_xz[i] - 10.0 * g_zz_zz_0_xz[i] * a_exp + 4.0 * g_zzzz_zz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_zz_0_yy[i] = 2.0 * g_0_zz_0_yy[i] - 10.0 * g_zz_zz_0_yy[i] * a_exp + 4.0 * g_zzzz_zz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_zz_0_yz[i] = 2.0 * g_0_zz_0_yz[i] - 10.0 * g_zz_zz_0_yz[i] * a_exp + 4.0 * g_zzzz_zz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_zz_zz_0_zz[i] = 2.0 * g_0_zz_0_zz[i] - 10.0 * g_zz_zz_0_zz[i] * a_exp + 4.0 * g_zzzz_zz_0_zz[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

