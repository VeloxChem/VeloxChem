#include "GeomDeriv2000OfScalarForSDSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_sdsd_0(CSimdArray<double>& buffer_2000_sdsd,
                     const CSimdArray<double>& buffer_sdsd,
                     const CSimdArray<double>& buffer_ddsd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_sdsd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_2000_sdsd

    auto g_xx_0_0_0_0_xx_0_xx = buffer_2000_sdsd[0];

    auto g_xx_0_0_0_0_xx_0_xy = buffer_2000_sdsd[1];

    auto g_xx_0_0_0_0_xx_0_xz = buffer_2000_sdsd[2];

    auto g_xx_0_0_0_0_xx_0_yy = buffer_2000_sdsd[3];

    auto g_xx_0_0_0_0_xx_0_yz = buffer_2000_sdsd[4];

    auto g_xx_0_0_0_0_xx_0_zz = buffer_2000_sdsd[5];

    auto g_xx_0_0_0_0_xy_0_xx = buffer_2000_sdsd[6];

    auto g_xx_0_0_0_0_xy_0_xy = buffer_2000_sdsd[7];

    auto g_xx_0_0_0_0_xy_0_xz = buffer_2000_sdsd[8];

    auto g_xx_0_0_0_0_xy_0_yy = buffer_2000_sdsd[9];

    auto g_xx_0_0_0_0_xy_0_yz = buffer_2000_sdsd[10];

    auto g_xx_0_0_0_0_xy_0_zz = buffer_2000_sdsd[11];

    auto g_xx_0_0_0_0_xz_0_xx = buffer_2000_sdsd[12];

    auto g_xx_0_0_0_0_xz_0_xy = buffer_2000_sdsd[13];

    auto g_xx_0_0_0_0_xz_0_xz = buffer_2000_sdsd[14];

    auto g_xx_0_0_0_0_xz_0_yy = buffer_2000_sdsd[15];

    auto g_xx_0_0_0_0_xz_0_yz = buffer_2000_sdsd[16];

    auto g_xx_0_0_0_0_xz_0_zz = buffer_2000_sdsd[17];

    auto g_xx_0_0_0_0_yy_0_xx = buffer_2000_sdsd[18];

    auto g_xx_0_0_0_0_yy_0_xy = buffer_2000_sdsd[19];

    auto g_xx_0_0_0_0_yy_0_xz = buffer_2000_sdsd[20];

    auto g_xx_0_0_0_0_yy_0_yy = buffer_2000_sdsd[21];

    auto g_xx_0_0_0_0_yy_0_yz = buffer_2000_sdsd[22];

    auto g_xx_0_0_0_0_yy_0_zz = buffer_2000_sdsd[23];

    auto g_xx_0_0_0_0_yz_0_xx = buffer_2000_sdsd[24];

    auto g_xx_0_0_0_0_yz_0_xy = buffer_2000_sdsd[25];

    auto g_xx_0_0_0_0_yz_0_xz = buffer_2000_sdsd[26];

    auto g_xx_0_0_0_0_yz_0_yy = buffer_2000_sdsd[27];

    auto g_xx_0_0_0_0_yz_0_yz = buffer_2000_sdsd[28];

    auto g_xx_0_0_0_0_yz_0_zz = buffer_2000_sdsd[29];

    auto g_xx_0_0_0_0_zz_0_xx = buffer_2000_sdsd[30];

    auto g_xx_0_0_0_0_zz_0_xy = buffer_2000_sdsd[31];

    auto g_xx_0_0_0_0_zz_0_xz = buffer_2000_sdsd[32];

    auto g_xx_0_0_0_0_zz_0_yy = buffer_2000_sdsd[33];

    auto g_xx_0_0_0_0_zz_0_yz = buffer_2000_sdsd[34];

    auto g_xx_0_0_0_0_zz_0_zz = buffer_2000_sdsd[35];

    auto g_xy_0_0_0_0_xx_0_xx = buffer_2000_sdsd[36];

    auto g_xy_0_0_0_0_xx_0_xy = buffer_2000_sdsd[37];

    auto g_xy_0_0_0_0_xx_0_xz = buffer_2000_sdsd[38];

    auto g_xy_0_0_0_0_xx_0_yy = buffer_2000_sdsd[39];

    auto g_xy_0_0_0_0_xx_0_yz = buffer_2000_sdsd[40];

    auto g_xy_0_0_0_0_xx_0_zz = buffer_2000_sdsd[41];

    auto g_xy_0_0_0_0_xy_0_xx = buffer_2000_sdsd[42];

    auto g_xy_0_0_0_0_xy_0_xy = buffer_2000_sdsd[43];

    auto g_xy_0_0_0_0_xy_0_xz = buffer_2000_sdsd[44];

    auto g_xy_0_0_0_0_xy_0_yy = buffer_2000_sdsd[45];

    auto g_xy_0_0_0_0_xy_0_yz = buffer_2000_sdsd[46];

    auto g_xy_0_0_0_0_xy_0_zz = buffer_2000_sdsd[47];

    auto g_xy_0_0_0_0_xz_0_xx = buffer_2000_sdsd[48];

    auto g_xy_0_0_0_0_xz_0_xy = buffer_2000_sdsd[49];

    auto g_xy_0_0_0_0_xz_0_xz = buffer_2000_sdsd[50];

    auto g_xy_0_0_0_0_xz_0_yy = buffer_2000_sdsd[51];

    auto g_xy_0_0_0_0_xz_0_yz = buffer_2000_sdsd[52];

    auto g_xy_0_0_0_0_xz_0_zz = buffer_2000_sdsd[53];

    auto g_xy_0_0_0_0_yy_0_xx = buffer_2000_sdsd[54];

    auto g_xy_0_0_0_0_yy_0_xy = buffer_2000_sdsd[55];

    auto g_xy_0_0_0_0_yy_0_xz = buffer_2000_sdsd[56];

    auto g_xy_0_0_0_0_yy_0_yy = buffer_2000_sdsd[57];

    auto g_xy_0_0_0_0_yy_0_yz = buffer_2000_sdsd[58];

    auto g_xy_0_0_0_0_yy_0_zz = buffer_2000_sdsd[59];

    auto g_xy_0_0_0_0_yz_0_xx = buffer_2000_sdsd[60];

    auto g_xy_0_0_0_0_yz_0_xy = buffer_2000_sdsd[61];

    auto g_xy_0_0_0_0_yz_0_xz = buffer_2000_sdsd[62];

    auto g_xy_0_0_0_0_yz_0_yy = buffer_2000_sdsd[63];

    auto g_xy_0_0_0_0_yz_0_yz = buffer_2000_sdsd[64];

    auto g_xy_0_0_0_0_yz_0_zz = buffer_2000_sdsd[65];

    auto g_xy_0_0_0_0_zz_0_xx = buffer_2000_sdsd[66];

    auto g_xy_0_0_0_0_zz_0_xy = buffer_2000_sdsd[67];

    auto g_xy_0_0_0_0_zz_0_xz = buffer_2000_sdsd[68];

    auto g_xy_0_0_0_0_zz_0_yy = buffer_2000_sdsd[69];

    auto g_xy_0_0_0_0_zz_0_yz = buffer_2000_sdsd[70];

    auto g_xy_0_0_0_0_zz_0_zz = buffer_2000_sdsd[71];

    auto g_xz_0_0_0_0_xx_0_xx = buffer_2000_sdsd[72];

    auto g_xz_0_0_0_0_xx_0_xy = buffer_2000_sdsd[73];

    auto g_xz_0_0_0_0_xx_0_xz = buffer_2000_sdsd[74];

    auto g_xz_0_0_0_0_xx_0_yy = buffer_2000_sdsd[75];

    auto g_xz_0_0_0_0_xx_0_yz = buffer_2000_sdsd[76];

    auto g_xz_0_0_0_0_xx_0_zz = buffer_2000_sdsd[77];

    auto g_xz_0_0_0_0_xy_0_xx = buffer_2000_sdsd[78];

    auto g_xz_0_0_0_0_xy_0_xy = buffer_2000_sdsd[79];

    auto g_xz_0_0_0_0_xy_0_xz = buffer_2000_sdsd[80];

    auto g_xz_0_0_0_0_xy_0_yy = buffer_2000_sdsd[81];

    auto g_xz_0_0_0_0_xy_0_yz = buffer_2000_sdsd[82];

    auto g_xz_0_0_0_0_xy_0_zz = buffer_2000_sdsd[83];

    auto g_xz_0_0_0_0_xz_0_xx = buffer_2000_sdsd[84];

    auto g_xz_0_0_0_0_xz_0_xy = buffer_2000_sdsd[85];

    auto g_xz_0_0_0_0_xz_0_xz = buffer_2000_sdsd[86];

    auto g_xz_0_0_0_0_xz_0_yy = buffer_2000_sdsd[87];

    auto g_xz_0_0_0_0_xz_0_yz = buffer_2000_sdsd[88];

    auto g_xz_0_0_0_0_xz_0_zz = buffer_2000_sdsd[89];

    auto g_xz_0_0_0_0_yy_0_xx = buffer_2000_sdsd[90];

    auto g_xz_0_0_0_0_yy_0_xy = buffer_2000_sdsd[91];

    auto g_xz_0_0_0_0_yy_0_xz = buffer_2000_sdsd[92];

    auto g_xz_0_0_0_0_yy_0_yy = buffer_2000_sdsd[93];

    auto g_xz_0_0_0_0_yy_0_yz = buffer_2000_sdsd[94];

    auto g_xz_0_0_0_0_yy_0_zz = buffer_2000_sdsd[95];

    auto g_xz_0_0_0_0_yz_0_xx = buffer_2000_sdsd[96];

    auto g_xz_0_0_0_0_yz_0_xy = buffer_2000_sdsd[97];

    auto g_xz_0_0_0_0_yz_0_xz = buffer_2000_sdsd[98];

    auto g_xz_0_0_0_0_yz_0_yy = buffer_2000_sdsd[99];

    auto g_xz_0_0_0_0_yz_0_yz = buffer_2000_sdsd[100];

    auto g_xz_0_0_0_0_yz_0_zz = buffer_2000_sdsd[101];

    auto g_xz_0_0_0_0_zz_0_xx = buffer_2000_sdsd[102];

    auto g_xz_0_0_0_0_zz_0_xy = buffer_2000_sdsd[103];

    auto g_xz_0_0_0_0_zz_0_xz = buffer_2000_sdsd[104];

    auto g_xz_0_0_0_0_zz_0_yy = buffer_2000_sdsd[105];

    auto g_xz_0_0_0_0_zz_0_yz = buffer_2000_sdsd[106];

    auto g_xz_0_0_0_0_zz_0_zz = buffer_2000_sdsd[107];

    auto g_yy_0_0_0_0_xx_0_xx = buffer_2000_sdsd[108];

    auto g_yy_0_0_0_0_xx_0_xy = buffer_2000_sdsd[109];

    auto g_yy_0_0_0_0_xx_0_xz = buffer_2000_sdsd[110];

    auto g_yy_0_0_0_0_xx_0_yy = buffer_2000_sdsd[111];

    auto g_yy_0_0_0_0_xx_0_yz = buffer_2000_sdsd[112];

    auto g_yy_0_0_0_0_xx_0_zz = buffer_2000_sdsd[113];

    auto g_yy_0_0_0_0_xy_0_xx = buffer_2000_sdsd[114];

    auto g_yy_0_0_0_0_xy_0_xy = buffer_2000_sdsd[115];

    auto g_yy_0_0_0_0_xy_0_xz = buffer_2000_sdsd[116];

    auto g_yy_0_0_0_0_xy_0_yy = buffer_2000_sdsd[117];

    auto g_yy_0_0_0_0_xy_0_yz = buffer_2000_sdsd[118];

    auto g_yy_0_0_0_0_xy_0_zz = buffer_2000_sdsd[119];

    auto g_yy_0_0_0_0_xz_0_xx = buffer_2000_sdsd[120];

    auto g_yy_0_0_0_0_xz_0_xy = buffer_2000_sdsd[121];

    auto g_yy_0_0_0_0_xz_0_xz = buffer_2000_sdsd[122];

    auto g_yy_0_0_0_0_xz_0_yy = buffer_2000_sdsd[123];

    auto g_yy_0_0_0_0_xz_0_yz = buffer_2000_sdsd[124];

    auto g_yy_0_0_0_0_xz_0_zz = buffer_2000_sdsd[125];

    auto g_yy_0_0_0_0_yy_0_xx = buffer_2000_sdsd[126];

    auto g_yy_0_0_0_0_yy_0_xy = buffer_2000_sdsd[127];

    auto g_yy_0_0_0_0_yy_0_xz = buffer_2000_sdsd[128];

    auto g_yy_0_0_0_0_yy_0_yy = buffer_2000_sdsd[129];

    auto g_yy_0_0_0_0_yy_0_yz = buffer_2000_sdsd[130];

    auto g_yy_0_0_0_0_yy_0_zz = buffer_2000_sdsd[131];

    auto g_yy_0_0_0_0_yz_0_xx = buffer_2000_sdsd[132];

    auto g_yy_0_0_0_0_yz_0_xy = buffer_2000_sdsd[133];

    auto g_yy_0_0_0_0_yz_0_xz = buffer_2000_sdsd[134];

    auto g_yy_0_0_0_0_yz_0_yy = buffer_2000_sdsd[135];

    auto g_yy_0_0_0_0_yz_0_yz = buffer_2000_sdsd[136];

    auto g_yy_0_0_0_0_yz_0_zz = buffer_2000_sdsd[137];

    auto g_yy_0_0_0_0_zz_0_xx = buffer_2000_sdsd[138];

    auto g_yy_0_0_0_0_zz_0_xy = buffer_2000_sdsd[139];

    auto g_yy_0_0_0_0_zz_0_xz = buffer_2000_sdsd[140];

    auto g_yy_0_0_0_0_zz_0_yy = buffer_2000_sdsd[141];

    auto g_yy_0_0_0_0_zz_0_yz = buffer_2000_sdsd[142];

    auto g_yy_0_0_0_0_zz_0_zz = buffer_2000_sdsd[143];

    auto g_yz_0_0_0_0_xx_0_xx = buffer_2000_sdsd[144];

    auto g_yz_0_0_0_0_xx_0_xy = buffer_2000_sdsd[145];

    auto g_yz_0_0_0_0_xx_0_xz = buffer_2000_sdsd[146];

    auto g_yz_0_0_0_0_xx_0_yy = buffer_2000_sdsd[147];

    auto g_yz_0_0_0_0_xx_0_yz = buffer_2000_sdsd[148];

    auto g_yz_0_0_0_0_xx_0_zz = buffer_2000_sdsd[149];

    auto g_yz_0_0_0_0_xy_0_xx = buffer_2000_sdsd[150];

    auto g_yz_0_0_0_0_xy_0_xy = buffer_2000_sdsd[151];

    auto g_yz_0_0_0_0_xy_0_xz = buffer_2000_sdsd[152];

    auto g_yz_0_0_0_0_xy_0_yy = buffer_2000_sdsd[153];

    auto g_yz_0_0_0_0_xy_0_yz = buffer_2000_sdsd[154];

    auto g_yz_0_0_0_0_xy_0_zz = buffer_2000_sdsd[155];

    auto g_yz_0_0_0_0_xz_0_xx = buffer_2000_sdsd[156];

    auto g_yz_0_0_0_0_xz_0_xy = buffer_2000_sdsd[157];

    auto g_yz_0_0_0_0_xz_0_xz = buffer_2000_sdsd[158];

    auto g_yz_0_0_0_0_xz_0_yy = buffer_2000_sdsd[159];

    auto g_yz_0_0_0_0_xz_0_yz = buffer_2000_sdsd[160];

    auto g_yz_0_0_0_0_xz_0_zz = buffer_2000_sdsd[161];

    auto g_yz_0_0_0_0_yy_0_xx = buffer_2000_sdsd[162];

    auto g_yz_0_0_0_0_yy_0_xy = buffer_2000_sdsd[163];

    auto g_yz_0_0_0_0_yy_0_xz = buffer_2000_sdsd[164];

    auto g_yz_0_0_0_0_yy_0_yy = buffer_2000_sdsd[165];

    auto g_yz_0_0_0_0_yy_0_yz = buffer_2000_sdsd[166];

    auto g_yz_0_0_0_0_yy_0_zz = buffer_2000_sdsd[167];

    auto g_yz_0_0_0_0_yz_0_xx = buffer_2000_sdsd[168];

    auto g_yz_0_0_0_0_yz_0_xy = buffer_2000_sdsd[169];

    auto g_yz_0_0_0_0_yz_0_xz = buffer_2000_sdsd[170];

    auto g_yz_0_0_0_0_yz_0_yy = buffer_2000_sdsd[171];

    auto g_yz_0_0_0_0_yz_0_yz = buffer_2000_sdsd[172];

    auto g_yz_0_0_0_0_yz_0_zz = buffer_2000_sdsd[173];

    auto g_yz_0_0_0_0_zz_0_xx = buffer_2000_sdsd[174];

    auto g_yz_0_0_0_0_zz_0_xy = buffer_2000_sdsd[175];

    auto g_yz_0_0_0_0_zz_0_xz = buffer_2000_sdsd[176];

    auto g_yz_0_0_0_0_zz_0_yy = buffer_2000_sdsd[177];

    auto g_yz_0_0_0_0_zz_0_yz = buffer_2000_sdsd[178];

    auto g_yz_0_0_0_0_zz_0_zz = buffer_2000_sdsd[179];

    auto g_zz_0_0_0_0_xx_0_xx = buffer_2000_sdsd[180];

    auto g_zz_0_0_0_0_xx_0_xy = buffer_2000_sdsd[181];

    auto g_zz_0_0_0_0_xx_0_xz = buffer_2000_sdsd[182];

    auto g_zz_0_0_0_0_xx_0_yy = buffer_2000_sdsd[183];

    auto g_zz_0_0_0_0_xx_0_yz = buffer_2000_sdsd[184];

    auto g_zz_0_0_0_0_xx_0_zz = buffer_2000_sdsd[185];

    auto g_zz_0_0_0_0_xy_0_xx = buffer_2000_sdsd[186];

    auto g_zz_0_0_0_0_xy_0_xy = buffer_2000_sdsd[187];

    auto g_zz_0_0_0_0_xy_0_xz = buffer_2000_sdsd[188];

    auto g_zz_0_0_0_0_xy_0_yy = buffer_2000_sdsd[189];

    auto g_zz_0_0_0_0_xy_0_yz = buffer_2000_sdsd[190];

    auto g_zz_0_0_0_0_xy_0_zz = buffer_2000_sdsd[191];

    auto g_zz_0_0_0_0_xz_0_xx = buffer_2000_sdsd[192];

    auto g_zz_0_0_0_0_xz_0_xy = buffer_2000_sdsd[193];

    auto g_zz_0_0_0_0_xz_0_xz = buffer_2000_sdsd[194];

    auto g_zz_0_0_0_0_xz_0_yy = buffer_2000_sdsd[195];

    auto g_zz_0_0_0_0_xz_0_yz = buffer_2000_sdsd[196];

    auto g_zz_0_0_0_0_xz_0_zz = buffer_2000_sdsd[197];

    auto g_zz_0_0_0_0_yy_0_xx = buffer_2000_sdsd[198];

    auto g_zz_0_0_0_0_yy_0_xy = buffer_2000_sdsd[199];

    auto g_zz_0_0_0_0_yy_0_xz = buffer_2000_sdsd[200];

    auto g_zz_0_0_0_0_yy_0_yy = buffer_2000_sdsd[201];

    auto g_zz_0_0_0_0_yy_0_yz = buffer_2000_sdsd[202];

    auto g_zz_0_0_0_0_yy_0_zz = buffer_2000_sdsd[203];

    auto g_zz_0_0_0_0_yz_0_xx = buffer_2000_sdsd[204];

    auto g_zz_0_0_0_0_yz_0_xy = buffer_2000_sdsd[205];

    auto g_zz_0_0_0_0_yz_0_xz = buffer_2000_sdsd[206];

    auto g_zz_0_0_0_0_yz_0_yy = buffer_2000_sdsd[207];

    auto g_zz_0_0_0_0_yz_0_yz = buffer_2000_sdsd[208];

    auto g_zz_0_0_0_0_yz_0_zz = buffer_2000_sdsd[209];

    auto g_zz_0_0_0_0_zz_0_xx = buffer_2000_sdsd[210];

    auto g_zz_0_0_0_0_zz_0_xy = buffer_2000_sdsd[211];

    auto g_zz_0_0_0_0_zz_0_xz = buffer_2000_sdsd[212];

    auto g_zz_0_0_0_0_zz_0_yy = buffer_2000_sdsd[213];

    auto g_zz_0_0_0_0_zz_0_yz = buffer_2000_sdsd[214];

    auto g_zz_0_0_0_0_zz_0_zz = buffer_2000_sdsd[215];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_xx_0_xx, g_0_xx_0_xy, g_0_xx_0_xz, g_0_xx_0_yy, g_0_xx_0_yz, g_0_xx_0_zz, g_xx_0_0_0_0_xx_0_xx, g_xx_0_0_0_0_xx_0_xy, g_xx_0_0_0_0_xx_0_xz, g_xx_0_0_0_0_xx_0_yy, g_xx_0_0_0_0_xx_0_yz, g_xx_0_0_0_0_xx_0_zz, g_xx_xx_0_xx, g_xx_xx_0_xy, g_xx_xx_0_xz, g_xx_xx_0_yy, g_xx_xx_0_yz, g_xx_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xx_0_xx[i] = -2.0 * g_0_xx_0_xx[i] * a_exp + 4.0 * g_xx_xx_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_0_xy[i] = -2.0 * g_0_xx_0_xy[i] * a_exp + 4.0 * g_xx_xx_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_0_xz[i] = -2.0 * g_0_xx_0_xz[i] * a_exp + 4.0 * g_xx_xx_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_0_yy[i] = -2.0 * g_0_xx_0_yy[i] * a_exp + 4.0 * g_xx_xx_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_0_yz[i] = -2.0 * g_0_xx_0_yz[i] * a_exp + 4.0 * g_xx_xx_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xx_0_zz[i] = -2.0 * g_0_xx_0_zz[i] * a_exp + 4.0 * g_xx_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_xx_0_0_0_0_xy_0_xx, g_xx_0_0_0_0_xy_0_xy, g_xx_0_0_0_0_xy_0_xz, g_xx_0_0_0_0_xy_0_yy, g_xx_0_0_0_0_xy_0_yz, g_xx_0_0_0_0_xy_0_zz, g_xx_xy_0_xx, g_xx_xy_0_xy, g_xx_xy_0_xz, g_xx_xy_0_yy, g_xx_xy_0_yz, g_xx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xy_0_xx[i] = -2.0 * g_0_xy_0_xx[i] * a_exp + 4.0 * g_xx_xy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_0_xy[i] = -2.0 * g_0_xy_0_xy[i] * a_exp + 4.0 * g_xx_xy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_0_xz[i] = -2.0 * g_0_xy_0_xz[i] * a_exp + 4.0 * g_xx_xy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_0_yy[i] = -2.0 * g_0_xy_0_yy[i] * a_exp + 4.0 * g_xx_xy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_0_yz[i] = -2.0 * g_0_xy_0_yz[i] * a_exp + 4.0 * g_xx_xy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xy_0_zz[i] = -2.0 * g_0_xy_0_zz[i] * a_exp + 4.0 * g_xx_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_xx_0_0_0_0_xz_0_xx, g_xx_0_0_0_0_xz_0_xy, g_xx_0_0_0_0_xz_0_xz, g_xx_0_0_0_0_xz_0_yy, g_xx_0_0_0_0_xz_0_yz, g_xx_0_0_0_0_xz_0_zz, g_xx_xz_0_xx, g_xx_xz_0_xy, g_xx_xz_0_xz, g_xx_xz_0_yy, g_xx_xz_0_yz, g_xx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_xz_0_xx[i] = -2.0 * g_0_xz_0_xx[i] * a_exp + 4.0 * g_xx_xz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_0_xy[i] = -2.0 * g_0_xz_0_xy[i] * a_exp + 4.0 * g_xx_xz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_0_xz[i] = -2.0 * g_0_xz_0_xz[i] * a_exp + 4.0 * g_xx_xz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_0_yy[i] = -2.0 * g_0_xz_0_yy[i] * a_exp + 4.0 * g_xx_xz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_0_yz[i] = -2.0 * g_0_xz_0_yz[i] * a_exp + 4.0 * g_xx_xz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_xz_0_zz[i] = -2.0 * g_0_xz_0_zz[i] * a_exp + 4.0 * g_xx_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_yy_0_xx, g_0_yy_0_xy, g_0_yy_0_xz, g_0_yy_0_yy, g_0_yy_0_yz, g_0_yy_0_zz, g_xx_0_0_0_0_yy_0_xx, g_xx_0_0_0_0_yy_0_xy, g_xx_0_0_0_0_yy_0_xz, g_xx_0_0_0_0_yy_0_yy, g_xx_0_0_0_0_yy_0_yz, g_xx_0_0_0_0_yy_0_zz, g_xx_yy_0_xx, g_xx_yy_0_xy, g_xx_yy_0_xz, g_xx_yy_0_yy, g_xx_yy_0_yz, g_xx_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yy_0_xx[i] = -2.0 * g_0_yy_0_xx[i] * a_exp + 4.0 * g_xx_yy_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_0_xy[i] = -2.0 * g_0_yy_0_xy[i] * a_exp + 4.0 * g_xx_yy_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_0_xz[i] = -2.0 * g_0_yy_0_xz[i] * a_exp + 4.0 * g_xx_yy_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_0_yy[i] = -2.0 * g_0_yy_0_yy[i] * a_exp + 4.0 * g_xx_yy_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_0_yz[i] = -2.0 * g_0_yy_0_yz[i] * a_exp + 4.0 * g_xx_yy_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yy_0_zz[i] = -2.0 * g_0_yy_0_zz[i] * a_exp + 4.0 * g_xx_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_xx_0_0_0_0_yz_0_xx, g_xx_0_0_0_0_yz_0_xy, g_xx_0_0_0_0_yz_0_xz, g_xx_0_0_0_0_yz_0_yy, g_xx_0_0_0_0_yz_0_yz, g_xx_0_0_0_0_yz_0_zz, g_xx_yz_0_xx, g_xx_yz_0_xy, g_xx_yz_0_xz, g_xx_yz_0_yy, g_xx_yz_0_yz, g_xx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_yz_0_xx[i] = -2.0 * g_0_yz_0_xx[i] * a_exp + 4.0 * g_xx_yz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_0_xy[i] = -2.0 * g_0_yz_0_xy[i] * a_exp + 4.0 * g_xx_yz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_0_xz[i] = -2.0 * g_0_yz_0_xz[i] * a_exp + 4.0 * g_xx_yz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_0_yy[i] = -2.0 * g_0_yz_0_yy[i] * a_exp + 4.0 * g_xx_yz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_0_yz[i] = -2.0 * g_0_yz_0_yz[i] * a_exp + 4.0 * g_xx_yz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_yz_0_zz[i] = -2.0 * g_0_yz_0_zz[i] * a_exp + 4.0 * g_xx_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_zz_0_xx, g_0_zz_0_xy, g_0_zz_0_xz, g_0_zz_0_yy, g_0_zz_0_yz, g_0_zz_0_zz, g_xx_0_0_0_0_zz_0_xx, g_xx_0_0_0_0_zz_0_xy, g_xx_0_0_0_0_zz_0_xz, g_xx_0_0_0_0_zz_0_yy, g_xx_0_0_0_0_zz_0_yz, g_xx_0_0_0_0_zz_0_zz, g_xx_zz_0_xx, g_xx_zz_0_xy, g_xx_zz_0_xz, g_xx_zz_0_yy, g_xx_zz_0_yz, g_xx_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_zz_0_xx[i] = -2.0 * g_0_zz_0_xx[i] * a_exp + 4.0 * g_xx_zz_0_xx[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_0_xy[i] = -2.0 * g_0_zz_0_xy[i] * a_exp + 4.0 * g_xx_zz_0_xy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_0_xz[i] = -2.0 * g_0_zz_0_xz[i] * a_exp + 4.0 * g_xx_zz_0_xz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_0_yy[i] = -2.0 * g_0_zz_0_yy[i] * a_exp + 4.0 * g_xx_zz_0_yy[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_0_yz[i] = -2.0 * g_0_zz_0_yz[i] * a_exp + 4.0 * g_xx_zz_0_yz[i] * a_exp * a_exp;

        g_xx_0_0_0_0_zz_0_zz[i] = -2.0 * g_0_zz_0_zz[i] * a_exp + 4.0 * g_xx_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_xy_0_0_0_0_xx_0_xx, g_xy_0_0_0_0_xx_0_xy, g_xy_0_0_0_0_xx_0_xz, g_xy_0_0_0_0_xx_0_yy, g_xy_0_0_0_0_xx_0_yz, g_xy_0_0_0_0_xx_0_zz, g_xy_xx_0_xx, g_xy_xx_0_xy, g_xy_xx_0_xz, g_xy_xx_0_yy, g_xy_xx_0_yz, g_xy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xx_0_xx[i] = 4.0 * g_xy_xx_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_0_xy[i] = 4.0 * g_xy_xx_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_0_xz[i] = 4.0 * g_xy_xx_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_0_yy[i] = 4.0 * g_xy_xx_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_0_yz[i] = 4.0 * g_xy_xx_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xx_0_zz[i] = 4.0 * g_xy_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_xy_0_0_0_0_xy_0_xx, g_xy_0_0_0_0_xy_0_xy, g_xy_0_0_0_0_xy_0_xz, g_xy_0_0_0_0_xy_0_yy, g_xy_0_0_0_0_xy_0_yz, g_xy_0_0_0_0_xy_0_zz, g_xy_xy_0_xx, g_xy_xy_0_xy, g_xy_xy_0_xz, g_xy_xy_0_yy, g_xy_xy_0_yz, g_xy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xy_0_xx[i] = 4.0 * g_xy_xy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_0_xy[i] = 4.0 * g_xy_xy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_0_xz[i] = 4.0 * g_xy_xy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_0_yy[i] = 4.0 * g_xy_xy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_0_yz[i] = 4.0 * g_xy_xy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xy_0_zz[i] = 4.0 * g_xy_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_xy_0_0_0_0_xz_0_xx, g_xy_0_0_0_0_xz_0_xy, g_xy_0_0_0_0_xz_0_xz, g_xy_0_0_0_0_xz_0_yy, g_xy_0_0_0_0_xz_0_yz, g_xy_0_0_0_0_xz_0_zz, g_xy_xz_0_xx, g_xy_xz_0_xy, g_xy_xz_0_xz, g_xy_xz_0_yy, g_xy_xz_0_yz, g_xy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_xz_0_xx[i] = 4.0 * g_xy_xz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_0_xy[i] = 4.0 * g_xy_xz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_0_xz[i] = 4.0 * g_xy_xz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_0_yy[i] = 4.0 * g_xy_xz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_0_yz[i] = 4.0 * g_xy_xz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_xz_0_zz[i] = 4.0 * g_xy_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_xy_0_0_0_0_yy_0_xx, g_xy_0_0_0_0_yy_0_xy, g_xy_0_0_0_0_yy_0_xz, g_xy_0_0_0_0_yy_0_yy, g_xy_0_0_0_0_yy_0_yz, g_xy_0_0_0_0_yy_0_zz, g_xy_yy_0_xx, g_xy_yy_0_xy, g_xy_yy_0_xz, g_xy_yy_0_yy, g_xy_yy_0_yz, g_xy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yy_0_xx[i] = 4.0 * g_xy_yy_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_0_xy[i] = 4.0 * g_xy_yy_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_0_xz[i] = 4.0 * g_xy_yy_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_0_yy[i] = 4.0 * g_xy_yy_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_0_yz[i] = 4.0 * g_xy_yy_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yy_0_zz[i] = 4.0 * g_xy_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_xy_0_0_0_0_yz_0_xx, g_xy_0_0_0_0_yz_0_xy, g_xy_0_0_0_0_yz_0_xz, g_xy_0_0_0_0_yz_0_yy, g_xy_0_0_0_0_yz_0_yz, g_xy_0_0_0_0_yz_0_zz, g_xy_yz_0_xx, g_xy_yz_0_xy, g_xy_yz_0_xz, g_xy_yz_0_yy, g_xy_yz_0_yz, g_xy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_yz_0_xx[i] = 4.0 * g_xy_yz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_0_xy[i] = 4.0 * g_xy_yz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_0_xz[i] = 4.0 * g_xy_yz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_0_yy[i] = 4.0 * g_xy_yz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_0_yz[i] = 4.0 * g_xy_yz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_yz_0_zz[i] = 4.0 * g_xy_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_xy_0_0_0_0_zz_0_xx, g_xy_0_0_0_0_zz_0_xy, g_xy_0_0_0_0_zz_0_xz, g_xy_0_0_0_0_zz_0_yy, g_xy_0_0_0_0_zz_0_yz, g_xy_0_0_0_0_zz_0_zz, g_xy_zz_0_xx, g_xy_zz_0_xy, g_xy_zz_0_xz, g_xy_zz_0_yy, g_xy_zz_0_yz, g_xy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_zz_0_xx[i] = 4.0 * g_xy_zz_0_xx[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_0_xy[i] = 4.0 * g_xy_zz_0_xy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_0_xz[i] = 4.0 * g_xy_zz_0_xz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_0_yy[i] = 4.0 * g_xy_zz_0_yy[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_0_yz[i] = 4.0 * g_xy_zz_0_yz[i] * a_exp * a_exp;

        g_xy_0_0_0_0_zz_0_zz[i] = 4.0 * g_xy_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_xz_0_0_0_0_xx_0_xx, g_xz_0_0_0_0_xx_0_xy, g_xz_0_0_0_0_xx_0_xz, g_xz_0_0_0_0_xx_0_yy, g_xz_0_0_0_0_xx_0_yz, g_xz_0_0_0_0_xx_0_zz, g_xz_xx_0_xx, g_xz_xx_0_xy, g_xz_xx_0_xz, g_xz_xx_0_yy, g_xz_xx_0_yz, g_xz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xx_0_xx[i] = 4.0 * g_xz_xx_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_0_xy[i] = 4.0 * g_xz_xx_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_0_xz[i] = 4.0 * g_xz_xx_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_0_yy[i] = 4.0 * g_xz_xx_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_0_yz[i] = 4.0 * g_xz_xx_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xx_0_zz[i] = 4.0 * g_xz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_xz_0_0_0_0_xy_0_xx, g_xz_0_0_0_0_xy_0_xy, g_xz_0_0_0_0_xy_0_xz, g_xz_0_0_0_0_xy_0_yy, g_xz_0_0_0_0_xy_0_yz, g_xz_0_0_0_0_xy_0_zz, g_xz_xy_0_xx, g_xz_xy_0_xy, g_xz_xy_0_xz, g_xz_xy_0_yy, g_xz_xy_0_yz, g_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xy_0_xx[i] = 4.0 * g_xz_xy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_0_xy[i] = 4.0 * g_xz_xy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_0_xz[i] = 4.0 * g_xz_xy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_0_yy[i] = 4.0 * g_xz_xy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_0_yz[i] = 4.0 * g_xz_xy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xy_0_zz[i] = 4.0 * g_xz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_xz_0_0_0_0_xz_0_xx, g_xz_0_0_0_0_xz_0_xy, g_xz_0_0_0_0_xz_0_xz, g_xz_0_0_0_0_xz_0_yy, g_xz_0_0_0_0_xz_0_yz, g_xz_0_0_0_0_xz_0_zz, g_xz_xz_0_xx, g_xz_xz_0_xy, g_xz_xz_0_xz, g_xz_xz_0_yy, g_xz_xz_0_yz, g_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_xz_0_xx[i] = 4.0 * g_xz_xz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_0_xy[i] = 4.0 * g_xz_xz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_0_xz[i] = 4.0 * g_xz_xz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_0_yy[i] = 4.0 * g_xz_xz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_0_yz[i] = 4.0 * g_xz_xz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_xz_0_zz[i] = 4.0 * g_xz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_xz_0_0_0_0_yy_0_xx, g_xz_0_0_0_0_yy_0_xy, g_xz_0_0_0_0_yy_0_xz, g_xz_0_0_0_0_yy_0_yy, g_xz_0_0_0_0_yy_0_yz, g_xz_0_0_0_0_yy_0_zz, g_xz_yy_0_xx, g_xz_yy_0_xy, g_xz_yy_0_xz, g_xz_yy_0_yy, g_xz_yy_0_yz, g_xz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yy_0_xx[i] = 4.0 * g_xz_yy_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_0_xy[i] = 4.0 * g_xz_yy_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_0_xz[i] = 4.0 * g_xz_yy_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_0_yy[i] = 4.0 * g_xz_yy_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_0_yz[i] = 4.0 * g_xz_yy_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yy_0_zz[i] = 4.0 * g_xz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_xz_0_0_0_0_yz_0_xx, g_xz_0_0_0_0_yz_0_xy, g_xz_0_0_0_0_yz_0_xz, g_xz_0_0_0_0_yz_0_yy, g_xz_0_0_0_0_yz_0_yz, g_xz_0_0_0_0_yz_0_zz, g_xz_yz_0_xx, g_xz_yz_0_xy, g_xz_yz_0_xz, g_xz_yz_0_yy, g_xz_yz_0_yz, g_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_yz_0_xx[i] = 4.0 * g_xz_yz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_0_xy[i] = 4.0 * g_xz_yz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_0_xz[i] = 4.0 * g_xz_yz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_0_yy[i] = 4.0 * g_xz_yz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_0_yz[i] = 4.0 * g_xz_yz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_yz_0_zz[i] = 4.0 * g_xz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_xz_0_0_0_0_zz_0_xx, g_xz_0_0_0_0_zz_0_xy, g_xz_0_0_0_0_zz_0_xz, g_xz_0_0_0_0_zz_0_yy, g_xz_0_0_0_0_zz_0_yz, g_xz_0_0_0_0_zz_0_zz, g_xz_zz_0_xx, g_xz_zz_0_xy, g_xz_zz_0_xz, g_xz_zz_0_yy, g_xz_zz_0_yz, g_xz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_zz_0_xx[i] = 4.0 * g_xz_zz_0_xx[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_0_xy[i] = 4.0 * g_xz_zz_0_xy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_0_xz[i] = 4.0 * g_xz_zz_0_xz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_0_yy[i] = 4.0 * g_xz_zz_0_yy[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_0_yz[i] = 4.0 * g_xz_zz_0_yz[i] * a_exp * a_exp;

        g_xz_0_0_0_0_zz_0_zz[i] = 4.0 * g_xz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_0_xx_0_xx, g_0_xx_0_xy, g_0_xx_0_xz, g_0_xx_0_yy, g_0_xx_0_yz, g_0_xx_0_zz, g_yy_0_0_0_0_xx_0_xx, g_yy_0_0_0_0_xx_0_xy, g_yy_0_0_0_0_xx_0_xz, g_yy_0_0_0_0_xx_0_yy, g_yy_0_0_0_0_xx_0_yz, g_yy_0_0_0_0_xx_0_zz, g_yy_xx_0_xx, g_yy_xx_0_xy, g_yy_xx_0_xz, g_yy_xx_0_yy, g_yy_xx_0_yz, g_yy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xx_0_xx[i] = -2.0 * g_0_xx_0_xx[i] * a_exp + 4.0 * g_yy_xx_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_0_xy[i] = -2.0 * g_0_xx_0_xy[i] * a_exp + 4.0 * g_yy_xx_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_0_xz[i] = -2.0 * g_0_xx_0_xz[i] * a_exp + 4.0 * g_yy_xx_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_0_yy[i] = -2.0 * g_0_xx_0_yy[i] * a_exp + 4.0 * g_yy_xx_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_0_yz[i] = -2.0 * g_0_xx_0_yz[i] * a_exp + 4.0 * g_yy_xx_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xx_0_zz[i] = -2.0 * g_0_xx_0_zz[i] * a_exp + 4.0 * g_yy_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_yy_0_0_0_0_xy_0_xx, g_yy_0_0_0_0_xy_0_xy, g_yy_0_0_0_0_xy_0_xz, g_yy_0_0_0_0_xy_0_yy, g_yy_0_0_0_0_xy_0_yz, g_yy_0_0_0_0_xy_0_zz, g_yy_xy_0_xx, g_yy_xy_0_xy, g_yy_xy_0_xz, g_yy_xy_0_yy, g_yy_xy_0_yz, g_yy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xy_0_xx[i] = -2.0 * g_0_xy_0_xx[i] * a_exp + 4.0 * g_yy_xy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_0_xy[i] = -2.0 * g_0_xy_0_xy[i] * a_exp + 4.0 * g_yy_xy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_0_xz[i] = -2.0 * g_0_xy_0_xz[i] * a_exp + 4.0 * g_yy_xy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_0_yy[i] = -2.0 * g_0_xy_0_yy[i] * a_exp + 4.0 * g_yy_xy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_0_yz[i] = -2.0 * g_0_xy_0_yz[i] * a_exp + 4.0 * g_yy_xy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xy_0_zz[i] = -2.0 * g_0_xy_0_zz[i] * a_exp + 4.0 * g_yy_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_yy_0_0_0_0_xz_0_xx, g_yy_0_0_0_0_xz_0_xy, g_yy_0_0_0_0_xz_0_xz, g_yy_0_0_0_0_xz_0_yy, g_yy_0_0_0_0_xz_0_yz, g_yy_0_0_0_0_xz_0_zz, g_yy_xz_0_xx, g_yy_xz_0_xy, g_yy_xz_0_xz, g_yy_xz_0_yy, g_yy_xz_0_yz, g_yy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_xz_0_xx[i] = -2.0 * g_0_xz_0_xx[i] * a_exp + 4.0 * g_yy_xz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_0_xy[i] = -2.0 * g_0_xz_0_xy[i] * a_exp + 4.0 * g_yy_xz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_0_xz[i] = -2.0 * g_0_xz_0_xz[i] * a_exp + 4.0 * g_yy_xz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_0_yy[i] = -2.0 * g_0_xz_0_yy[i] * a_exp + 4.0 * g_yy_xz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_0_yz[i] = -2.0 * g_0_xz_0_yz[i] * a_exp + 4.0 * g_yy_xz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_xz_0_zz[i] = -2.0 * g_0_xz_0_zz[i] * a_exp + 4.0 * g_yy_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_0_yy_0_xx, g_0_yy_0_xy, g_0_yy_0_xz, g_0_yy_0_yy, g_0_yy_0_yz, g_0_yy_0_zz, g_yy_0_0_0_0_yy_0_xx, g_yy_0_0_0_0_yy_0_xy, g_yy_0_0_0_0_yy_0_xz, g_yy_0_0_0_0_yy_0_yy, g_yy_0_0_0_0_yy_0_yz, g_yy_0_0_0_0_yy_0_zz, g_yy_yy_0_xx, g_yy_yy_0_xy, g_yy_yy_0_xz, g_yy_yy_0_yy, g_yy_yy_0_yz, g_yy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yy_0_xx[i] = -2.0 * g_0_yy_0_xx[i] * a_exp + 4.0 * g_yy_yy_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_0_xy[i] = -2.0 * g_0_yy_0_xy[i] * a_exp + 4.0 * g_yy_yy_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_0_xz[i] = -2.0 * g_0_yy_0_xz[i] * a_exp + 4.0 * g_yy_yy_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_0_yy[i] = -2.0 * g_0_yy_0_yy[i] * a_exp + 4.0 * g_yy_yy_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_0_yz[i] = -2.0 * g_0_yy_0_yz[i] * a_exp + 4.0 * g_yy_yy_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yy_0_zz[i] = -2.0 * g_0_yy_0_zz[i] * a_exp + 4.0 * g_yy_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_yy_0_0_0_0_yz_0_xx, g_yy_0_0_0_0_yz_0_xy, g_yy_0_0_0_0_yz_0_xz, g_yy_0_0_0_0_yz_0_yy, g_yy_0_0_0_0_yz_0_yz, g_yy_0_0_0_0_yz_0_zz, g_yy_yz_0_xx, g_yy_yz_0_xy, g_yy_yz_0_xz, g_yy_yz_0_yy, g_yy_yz_0_yz, g_yy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_yz_0_xx[i] = -2.0 * g_0_yz_0_xx[i] * a_exp + 4.0 * g_yy_yz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_0_xy[i] = -2.0 * g_0_yz_0_xy[i] * a_exp + 4.0 * g_yy_yz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_0_xz[i] = -2.0 * g_0_yz_0_xz[i] * a_exp + 4.0 * g_yy_yz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_0_yy[i] = -2.0 * g_0_yz_0_yy[i] * a_exp + 4.0 * g_yy_yz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_0_yz[i] = -2.0 * g_0_yz_0_yz[i] * a_exp + 4.0 * g_yy_yz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_yz_0_zz[i] = -2.0 * g_0_yz_0_zz[i] * a_exp + 4.0 * g_yy_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_0_zz_0_xx, g_0_zz_0_xy, g_0_zz_0_xz, g_0_zz_0_yy, g_0_zz_0_yz, g_0_zz_0_zz, g_yy_0_0_0_0_zz_0_xx, g_yy_0_0_0_0_zz_0_xy, g_yy_0_0_0_0_zz_0_xz, g_yy_0_0_0_0_zz_0_yy, g_yy_0_0_0_0_zz_0_yz, g_yy_0_0_0_0_zz_0_zz, g_yy_zz_0_xx, g_yy_zz_0_xy, g_yy_zz_0_xz, g_yy_zz_0_yy, g_yy_zz_0_yz, g_yy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_zz_0_xx[i] = -2.0 * g_0_zz_0_xx[i] * a_exp + 4.0 * g_yy_zz_0_xx[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_0_xy[i] = -2.0 * g_0_zz_0_xy[i] * a_exp + 4.0 * g_yy_zz_0_xy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_0_xz[i] = -2.0 * g_0_zz_0_xz[i] * a_exp + 4.0 * g_yy_zz_0_xz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_0_yy[i] = -2.0 * g_0_zz_0_yy[i] * a_exp + 4.0 * g_yy_zz_0_yy[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_0_yz[i] = -2.0 * g_0_zz_0_yz[i] * a_exp + 4.0 * g_yy_zz_0_yz[i] * a_exp * a_exp;

        g_yy_0_0_0_0_zz_0_zz[i] = -2.0 * g_0_zz_0_zz[i] * a_exp + 4.0 * g_yy_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_yz_0_0_0_0_xx_0_xx, g_yz_0_0_0_0_xx_0_xy, g_yz_0_0_0_0_xx_0_xz, g_yz_0_0_0_0_xx_0_yy, g_yz_0_0_0_0_xx_0_yz, g_yz_0_0_0_0_xx_0_zz, g_yz_xx_0_xx, g_yz_xx_0_xy, g_yz_xx_0_xz, g_yz_xx_0_yy, g_yz_xx_0_yz, g_yz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xx_0_xx[i] = 4.0 * g_yz_xx_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_0_xy[i] = 4.0 * g_yz_xx_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_0_xz[i] = 4.0 * g_yz_xx_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_0_yy[i] = 4.0 * g_yz_xx_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_0_yz[i] = 4.0 * g_yz_xx_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xx_0_zz[i] = 4.0 * g_yz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_yz_0_0_0_0_xy_0_xx, g_yz_0_0_0_0_xy_0_xy, g_yz_0_0_0_0_xy_0_xz, g_yz_0_0_0_0_xy_0_yy, g_yz_0_0_0_0_xy_0_yz, g_yz_0_0_0_0_xy_0_zz, g_yz_xy_0_xx, g_yz_xy_0_xy, g_yz_xy_0_xz, g_yz_xy_0_yy, g_yz_xy_0_yz, g_yz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xy_0_xx[i] = 4.0 * g_yz_xy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_0_xy[i] = 4.0 * g_yz_xy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_0_xz[i] = 4.0 * g_yz_xy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_0_yy[i] = 4.0 * g_yz_xy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_0_yz[i] = 4.0 * g_yz_xy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xy_0_zz[i] = 4.0 * g_yz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_yz_0_0_0_0_xz_0_xx, g_yz_0_0_0_0_xz_0_xy, g_yz_0_0_0_0_xz_0_xz, g_yz_0_0_0_0_xz_0_yy, g_yz_0_0_0_0_xz_0_yz, g_yz_0_0_0_0_xz_0_zz, g_yz_xz_0_xx, g_yz_xz_0_xy, g_yz_xz_0_xz, g_yz_xz_0_yy, g_yz_xz_0_yz, g_yz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_xz_0_xx[i] = 4.0 * g_yz_xz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_0_xy[i] = 4.0 * g_yz_xz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_0_xz[i] = 4.0 * g_yz_xz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_0_yy[i] = 4.0 * g_yz_xz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_0_yz[i] = 4.0 * g_yz_xz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_xz_0_zz[i] = 4.0 * g_yz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_yz_0_0_0_0_yy_0_xx, g_yz_0_0_0_0_yy_0_xy, g_yz_0_0_0_0_yy_0_xz, g_yz_0_0_0_0_yy_0_yy, g_yz_0_0_0_0_yy_0_yz, g_yz_0_0_0_0_yy_0_zz, g_yz_yy_0_xx, g_yz_yy_0_xy, g_yz_yy_0_xz, g_yz_yy_0_yy, g_yz_yy_0_yz, g_yz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yy_0_xx[i] = 4.0 * g_yz_yy_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_0_xy[i] = 4.0 * g_yz_yy_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_0_xz[i] = 4.0 * g_yz_yy_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_0_yy[i] = 4.0 * g_yz_yy_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_0_yz[i] = 4.0 * g_yz_yy_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yy_0_zz[i] = 4.0 * g_yz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_yz_0_0_0_0_yz_0_xx, g_yz_0_0_0_0_yz_0_xy, g_yz_0_0_0_0_yz_0_xz, g_yz_0_0_0_0_yz_0_yy, g_yz_0_0_0_0_yz_0_yz, g_yz_0_0_0_0_yz_0_zz, g_yz_yz_0_xx, g_yz_yz_0_xy, g_yz_yz_0_xz, g_yz_yz_0_yy, g_yz_yz_0_yz, g_yz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_yz_0_xx[i] = 4.0 * g_yz_yz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_0_xy[i] = 4.0 * g_yz_yz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_0_xz[i] = 4.0 * g_yz_yz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_0_yy[i] = 4.0 * g_yz_yz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_0_yz[i] = 4.0 * g_yz_yz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_yz_0_zz[i] = 4.0 * g_yz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_yz_0_0_0_0_zz_0_xx, g_yz_0_0_0_0_zz_0_xy, g_yz_0_0_0_0_zz_0_xz, g_yz_0_0_0_0_zz_0_yy, g_yz_0_0_0_0_zz_0_yz, g_yz_0_0_0_0_zz_0_zz, g_yz_zz_0_xx, g_yz_zz_0_xy, g_yz_zz_0_xz, g_yz_zz_0_yy, g_yz_zz_0_yz, g_yz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_zz_0_xx[i] = 4.0 * g_yz_zz_0_xx[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_0_xy[i] = 4.0 * g_yz_zz_0_xy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_0_xz[i] = 4.0 * g_yz_zz_0_xz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_0_yy[i] = 4.0 * g_yz_zz_0_yy[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_0_yz[i] = 4.0 * g_yz_zz_0_yz[i] * a_exp * a_exp;

        g_yz_0_0_0_0_zz_0_zz[i] = 4.0 * g_yz_zz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_0_xx_0_xx, g_0_xx_0_xy, g_0_xx_0_xz, g_0_xx_0_yy, g_0_xx_0_yz, g_0_xx_0_zz, g_zz_0_0_0_0_xx_0_xx, g_zz_0_0_0_0_xx_0_xy, g_zz_0_0_0_0_xx_0_xz, g_zz_0_0_0_0_xx_0_yy, g_zz_0_0_0_0_xx_0_yz, g_zz_0_0_0_0_xx_0_zz, g_zz_xx_0_xx, g_zz_xx_0_xy, g_zz_xx_0_xz, g_zz_xx_0_yy, g_zz_xx_0_yz, g_zz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xx_0_xx[i] = -2.0 * g_0_xx_0_xx[i] * a_exp + 4.0 * g_zz_xx_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_0_xy[i] = -2.0 * g_0_xx_0_xy[i] * a_exp + 4.0 * g_zz_xx_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_0_xz[i] = -2.0 * g_0_xx_0_xz[i] * a_exp + 4.0 * g_zz_xx_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_0_yy[i] = -2.0 * g_0_xx_0_yy[i] * a_exp + 4.0 * g_zz_xx_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_0_yz[i] = -2.0 * g_0_xx_0_yz[i] * a_exp + 4.0 * g_zz_xx_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xx_0_zz[i] = -2.0 * g_0_xx_0_zz[i] * a_exp + 4.0 * g_zz_xx_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_zz_0_0_0_0_xy_0_xx, g_zz_0_0_0_0_xy_0_xy, g_zz_0_0_0_0_xy_0_xz, g_zz_0_0_0_0_xy_0_yy, g_zz_0_0_0_0_xy_0_yz, g_zz_0_0_0_0_xy_0_zz, g_zz_xy_0_xx, g_zz_xy_0_xy, g_zz_xy_0_xz, g_zz_xy_0_yy, g_zz_xy_0_yz, g_zz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xy_0_xx[i] = -2.0 * g_0_xy_0_xx[i] * a_exp + 4.0 * g_zz_xy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_0_xy[i] = -2.0 * g_0_xy_0_xy[i] * a_exp + 4.0 * g_zz_xy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_0_xz[i] = -2.0 * g_0_xy_0_xz[i] * a_exp + 4.0 * g_zz_xy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_0_yy[i] = -2.0 * g_0_xy_0_yy[i] * a_exp + 4.0 * g_zz_xy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_0_yz[i] = -2.0 * g_0_xy_0_yz[i] * a_exp + 4.0 * g_zz_xy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xy_0_zz[i] = -2.0 * g_0_xy_0_zz[i] * a_exp + 4.0 * g_zz_xy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_zz_0_0_0_0_xz_0_xx, g_zz_0_0_0_0_xz_0_xy, g_zz_0_0_0_0_xz_0_xz, g_zz_0_0_0_0_xz_0_yy, g_zz_0_0_0_0_xz_0_yz, g_zz_0_0_0_0_xz_0_zz, g_zz_xz_0_xx, g_zz_xz_0_xy, g_zz_xz_0_xz, g_zz_xz_0_yy, g_zz_xz_0_yz, g_zz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_xz_0_xx[i] = -2.0 * g_0_xz_0_xx[i] * a_exp + 4.0 * g_zz_xz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_0_xy[i] = -2.0 * g_0_xz_0_xy[i] * a_exp + 4.0 * g_zz_xz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_0_xz[i] = -2.0 * g_0_xz_0_xz[i] * a_exp + 4.0 * g_zz_xz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_0_yy[i] = -2.0 * g_0_xz_0_yy[i] * a_exp + 4.0 * g_zz_xz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_0_yz[i] = -2.0 * g_0_xz_0_yz[i] * a_exp + 4.0 * g_zz_xz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_xz_0_zz[i] = -2.0 * g_0_xz_0_zz[i] * a_exp + 4.0 * g_zz_xz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_0_yy_0_xx, g_0_yy_0_xy, g_0_yy_0_xz, g_0_yy_0_yy, g_0_yy_0_yz, g_0_yy_0_zz, g_zz_0_0_0_0_yy_0_xx, g_zz_0_0_0_0_yy_0_xy, g_zz_0_0_0_0_yy_0_xz, g_zz_0_0_0_0_yy_0_yy, g_zz_0_0_0_0_yy_0_yz, g_zz_0_0_0_0_yy_0_zz, g_zz_yy_0_xx, g_zz_yy_0_xy, g_zz_yy_0_xz, g_zz_yy_0_yy, g_zz_yy_0_yz, g_zz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yy_0_xx[i] = -2.0 * g_0_yy_0_xx[i] * a_exp + 4.0 * g_zz_yy_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_0_xy[i] = -2.0 * g_0_yy_0_xy[i] * a_exp + 4.0 * g_zz_yy_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_0_xz[i] = -2.0 * g_0_yy_0_xz[i] * a_exp + 4.0 * g_zz_yy_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_0_yy[i] = -2.0 * g_0_yy_0_yy[i] * a_exp + 4.0 * g_zz_yy_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_0_yz[i] = -2.0 * g_0_yy_0_yz[i] * a_exp + 4.0 * g_zz_yy_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yy_0_zz[i] = -2.0 * g_0_yy_0_zz[i] * a_exp + 4.0 * g_zz_yy_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_zz_0_0_0_0_yz_0_xx, g_zz_0_0_0_0_yz_0_xy, g_zz_0_0_0_0_yz_0_xz, g_zz_0_0_0_0_yz_0_yy, g_zz_0_0_0_0_yz_0_yz, g_zz_0_0_0_0_yz_0_zz, g_zz_yz_0_xx, g_zz_yz_0_xy, g_zz_yz_0_xz, g_zz_yz_0_yy, g_zz_yz_0_yz, g_zz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_yz_0_xx[i] = -2.0 * g_0_yz_0_xx[i] * a_exp + 4.0 * g_zz_yz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_0_xy[i] = -2.0 * g_0_yz_0_xy[i] * a_exp + 4.0 * g_zz_yz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_0_xz[i] = -2.0 * g_0_yz_0_xz[i] * a_exp + 4.0 * g_zz_yz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_0_yy[i] = -2.0 * g_0_yz_0_yy[i] * a_exp + 4.0 * g_zz_yz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_0_yz[i] = -2.0 * g_0_yz_0_yz[i] * a_exp + 4.0 * g_zz_yz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_yz_0_zz[i] = -2.0 * g_0_yz_0_zz[i] * a_exp + 4.0 * g_zz_yz_0_zz[i] * a_exp * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_0_zz_0_xx, g_0_zz_0_xy, g_0_zz_0_xz, g_0_zz_0_yy, g_0_zz_0_yz, g_0_zz_0_zz, g_zz_0_0_0_0_zz_0_xx, g_zz_0_0_0_0_zz_0_xy, g_zz_0_0_0_0_zz_0_xz, g_zz_0_0_0_0_zz_0_yy, g_zz_0_0_0_0_zz_0_yz, g_zz_0_0_0_0_zz_0_zz, g_zz_zz_0_xx, g_zz_zz_0_xy, g_zz_zz_0_xz, g_zz_zz_0_yy, g_zz_zz_0_yz, g_zz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_zz_0_xx[i] = -2.0 * g_0_zz_0_xx[i] * a_exp + 4.0 * g_zz_zz_0_xx[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_0_xy[i] = -2.0 * g_0_zz_0_xy[i] * a_exp + 4.0 * g_zz_zz_0_xy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_0_xz[i] = -2.0 * g_0_zz_0_xz[i] * a_exp + 4.0 * g_zz_zz_0_xz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_0_yy[i] = -2.0 * g_0_zz_0_yy[i] * a_exp + 4.0 * g_zz_zz_0_yy[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_0_yz[i] = -2.0 * g_0_zz_0_yz[i] * a_exp + 4.0 * g_zz_zz_0_yz[i] * a_exp * a_exp;

        g_zz_0_0_0_0_zz_0_zz[i] = -2.0 * g_0_zz_0_zz[i] * a_exp + 4.0 * g_zz_zz_0_zz[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

