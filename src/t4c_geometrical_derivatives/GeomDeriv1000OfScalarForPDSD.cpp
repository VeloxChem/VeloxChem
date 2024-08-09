#include "GeomDeriv1000OfScalarForPDSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_pdsd_0(CSimdArray<double>& buffer_1000_pdsd,
                     const CSimdArray<double>& buffer_sdsd,
                     const CSimdArray<double>& buffer_ddsd,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_pdsd.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_pdsd

    auto g_x_0_0_0_x_xx_0_xx = buffer_1000_pdsd[0];

    auto g_x_0_0_0_x_xx_0_xy = buffer_1000_pdsd[1];

    auto g_x_0_0_0_x_xx_0_xz = buffer_1000_pdsd[2];

    auto g_x_0_0_0_x_xx_0_yy = buffer_1000_pdsd[3];

    auto g_x_0_0_0_x_xx_0_yz = buffer_1000_pdsd[4];

    auto g_x_0_0_0_x_xx_0_zz = buffer_1000_pdsd[5];

    auto g_x_0_0_0_x_xy_0_xx = buffer_1000_pdsd[6];

    auto g_x_0_0_0_x_xy_0_xy = buffer_1000_pdsd[7];

    auto g_x_0_0_0_x_xy_0_xz = buffer_1000_pdsd[8];

    auto g_x_0_0_0_x_xy_0_yy = buffer_1000_pdsd[9];

    auto g_x_0_0_0_x_xy_0_yz = buffer_1000_pdsd[10];

    auto g_x_0_0_0_x_xy_0_zz = buffer_1000_pdsd[11];

    auto g_x_0_0_0_x_xz_0_xx = buffer_1000_pdsd[12];

    auto g_x_0_0_0_x_xz_0_xy = buffer_1000_pdsd[13];

    auto g_x_0_0_0_x_xz_0_xz = buffer_1000_pdsd[14];

    auto g_x_0_0_0_x_xz_0_yy = buffer_1000_pdsd[15];

    auto g_x_0_0_0_x_xz_0_yz = buffer_1000_pdsd[16];

    auto g_x_0_0_0_x_xz_0_zz = buffer_1000_pdsd[17];

    auto g_x_0_0_0_x_yy_0_xx = buffer_1000_pdsd[18];

    auto g_x_0_0_0_x_yy_0_xy = buffer_1000_pdsd[19];

    auto g_x_0_0_0_x_yy_0_xz = buffer_1000_pdsd[20];

    auto g_x_0_0_0_x_yy_0_yy = buffer_1000_pdsd[21];

    auto g_x_0_0_0_x_yy_0_yz = buffer_1000_pdsd[22];

    auto g_x_0_0_0_x_yy_0_zz = buffer_1000_pdsd[23];

    auto g_x_0_0_0_x_yz_0_xx = buffer_1000_pdsd[24];

    auto g_x_0_0_0_x_yz_0_xy = buffer_1000_pdsd[25];

    auto g_x_0_0_0_x_yz_0_xz = buffer_1000_pdsd[26];

    auto g_x_0_0_0_x_yz_0_yy = buffer_1000_pdsd[27];

    auto g_x_0_0_0_x_yz_0_yz = buffer_1000_pdsd[28];

    auto g_x_0_0_0_x_yz_0_zz = buffer_1000_pdsd[29];

    auto g_x_0_0_0_x_zz_0_xx = buffer_1000_pdsd[30];

    auto g_x_0_0_0_x_zz_0_xy = buffer_1000_pdsd[31];

    auto g_x_0_0_0_x_zz_0_xz = buffer_1000_pdsd[32];

    auto g_x_0_0_0_x_zz_0_yy = buffer_1000_pdsd[33];

    auto g_x_0_0_0_x_zz_0_yz = buffer_1000_pdsd[34];

    auto g_x_0_0_0_x_zz_0_zz = buffer_1000_pdsd[35];

    auto g_x_0_0_0_y_xx_0_xx = buffer_1000_pdsd[36];

    auto g_x_0_0_0_y_xx_0_xy = buffer_1000_pdsd[37];

    auto g_x_0_0_0_y_xx_0_xz = buffer_1000_pdsd[38];

    auto g_x_0_0_0_y_xx_0_yy = buffer_1000_pdsd[39];

    auto g_x_0_0_0_y_xx_0_yz = buffer_1000_pdsd[40];

    auto g_x_0_0_0_y_xx_0_zz = buffer_1000_pdsd[41];

    auto g_x_0_0_0_y_xy_0_xx = buffer_1000_pdsd[42];

    auto g_x_0_0_0_y_xy_0_xy = buffer_1000_pdsd[43];

    auto g_x_0_0_0_y_xy_0_xz = buffer_1000_pdsd[44];

    auto g_x_0_0_0_y_xy_0_yy = buffer_1000_pdsd[45];

    auto g_x_0_0_0_y_xy_0_yz = buffer_1000_pdsd[46];

    auto g_x_0_0_0_y_xy_0_zz = buffer_1000_pdsd[47];

    auto g_x_0_0_0_y_xz_0_xx = buffer_1000_pdsd[48];

    auto g_x_0_0_0_y_xz_0_xy = buffer_1000_pdsd[49];

    auto g_x_0_0_0_y_xz_0_xz = buffer_1000_pdsd[50];

    auto g_x_0_0_0_y_xz_0_yy = buffer_1000_pdsd[51];

    auto g_x_0_0_0_y_xz_0_yz = buffer_1000_pdsd[52];

    auto g_x_0_0_0_y_xz_0_zz = buffer_1000_pdsd[53];

    auto g_x_0_0_0_y_yy_0_xx = buffer_1000_pdsd[54];

    auto g_x_0_0_0_y_yy_0_xy = buffer_1000_pdsd[55];

    auto g_x_0_0_0_y_yy_0_xz = buffer_1000_pdsd[56];

    auto g_x_0_0_0_y_yy_0_yy = buffer_1000_pdsd[57];

    auto g_x_0_0_0_y_yy_0_yz = buffer_1000_pdsd[58];

    auto g_x_0_0_0_y_yy_0_zz = buffer_1000_pdsd[59];

    auto g_x_0_0_0_y_yz_0_xx = buffer_1000_pdsd[60];

    auto g_x_0_0_0_y_yz_0_xy = buffer_1000_pdsd[61];

    auto g_x_0_0_0_y_yz_0_xz = buffer_1000_pdsd[62];

    auto g_x_0_0_0_y_yz_0_yy = buffer_1000_pdsd[63];

    auto g_x_0_0_0_y_yz_0_yz = buffer_1000_pdsd[64];

    auto g_x_0_0_0_y_yz_0_zz = buffer_1000_pdsd[65];

    auto g_x_0_0_0_y_zz_0_xx = buffer_1000_pdsd[66];

    auto g_x_0_0_0_y_zz_0_xy = buffer_1000_pdsd[67];

    auto g_x_0_0_0_y_zz_0_xz = buffer_1000_pdsd[68];

    auto g_x_0_0_0_y_zz_0_yy = buffer_1000_pdsd[69];

    auto g_x_0_0_0_y_zz_0_yz = buffer_1000_pdsd[70];

    auto g_x_0_0_0_y_zz_0_zz = buffer_1000_pdsd[71];

    auto g_x_0_0_0_z_xx_0_xx = buffer_1000_pdsd[72];

    auto g_x_0_0_0_z_xx_0_xy = buffer_1000_pdsd[73];

    auto g_x_0_0_0_z_xx_0_xz = buffer_1000_pdsd[74];

    auto g_x_0_0_0_z_xx_0_yy = buffer_1000_pdsd[75];

    auto g_x_0_0_0_z_xx_0_yz = buffer_1000_pdsd[76];

    auto g_x_0_0_0_z_xx_0_zz = buffer_1000_pdsd[77];

    auto g_x_0_0_0_z_xy_0_xx = buffer_1000_pdsd[78];

    auto g_x_0_0_0_z_xy_0_xy = buffer_1000_pdsd[79];

    auto g_x_0_0_0_z_xy_0_xz = buffer_1000_pdsd[80];

    auto g_x_0_0_0_z_xy_0_yy = buffer_1000_pdsd[81];

    auto g_x_0_0_0_z_xy_0_yz = buffer_1000_pdsd[82];

    auto g_x_0_0_0_z_xy_0_zz = buffer_1000_pdsd[83];

    auto g_x_0_0_0_z_xz_0_xx = buffer_1000_pdsd[84];

    auto g_x_0_0_0_z_xz_0_xy = buffer_1000_pdsd[85];

    auto g_x_0_0_0_z_xz_0_xz = buffer_1000_pdsd[86];

    auto g_x_0_0_0_z_xz_0_yy = buffer_1000_pdsd[87];

    auto g_x_0_0_0_z_xz_0_yz = buffer_1000_pdsd[88];

    auto g_x_0_0_0_z_xz_0_zz = buffer_1000_pdsd[89];

    auto g_x_0_0_0_z_yy_0_xx = buffer_1000_pdsd[90];

    auto g_x_0_0_0_z_yy_0_xy = buffer_1000_pdsd[91];

    auto g_x_0_0_0_z_yy_0_xz = buffer_1000_pdsd[92];

    auto g_x_0_0_0_z_yy_0_yy = buffer_1000_pdsd[93];

    auto g_x_0_0_0_z_yy_0_yz = buffer_1000_pdsd[94];

    auto g_x_0_0_0_z_yy_0_zz = buffer_1000_pdsd[95];

    auto g_x_0_0_0_z_yz_0_xx = buffer_1000_pdsd[96];

    auto g_x_0_0_0_z_yz_0_xy = buffer_1000_pdsd[97];

    auto g_x_0_0_0_z_yz_0_xz = buffer_1000_pdsd[98];

    auto g_x_0_0_0_z_yz_0_yy = buffer_1000_pdsd[99];

    auto g_x_0_0_0_z_yz_0_yz = buffer_1000_pdsd[100];

    auto g_x_0_0_0_z_yz_0_zz = buffer_1000_pdsd[101];

    auto g_x_0_0_0_z_zz_0_xx = buffer_1000_pdsd[102];

    auto g_x_0_0_0_z_zz_0_xy = buffer_1000_pdsd[103];

    auto g_x_0_0_0_z_zz_0_xz = buffer_1000_pdsd[104];

    auto g_x_0_0_0_z_zz_0_yy = buffer_1000_pdsd[105];

    auto g_x_0_0_0_z_zz_0_yz = buffer_1000_pdsd[106];

    auto g_x_0_0_0_z_zz_0_zz = buffer_1000_pdsd[107];

    auto g_y_0_0_0_x_xx_0_xx = buffer_1000_pdsd[108];

    auto g_y_0_0_0_x_xx_0_xy = buffer_1000_pdsd[109];

    auto g_y_0_0_0_x_xx_0_xz = buffer_1000_pdsd[110];

    auto g_y_0_0_0_x_xx_0_yy = buffer_1000_pdsd[111];

    auto g_y_0_0_0_x_xx_0_yz = buffer_1000_pdsd[112];

    auto g_y_0_0_0_x_xx_0_zz = buffer_1000_pdsd[113];

    auto g_y_0_0_0_x_xy_0_xx = buffer_1000_pdsd[114];

    auto g_y_0_0_0_x_xy_0_xy = buffer_1000_pdsd[115];

    auto g_y_0_0_0_x_xy_0_xz = buffer_1000_pdsd[116];

    auto g_y_0_0_0_x_xy_0_yy = buffer_1000_pdsd[117];

    auto g_y_0_0_0_x_xy_0_yz = buffer_1000_pdsd[118];

    auto g_y_0_0_0_x_xy_0_zz = buffer_1000_pdsd[119];

    auto g_y_0_0_0_x_xz_0_xx = buffer_1000_pdsd[120];

    auto g_y_0_0_0_x_xz_0_xy = buffer_1000_pdsd[121];

    auto g_y_0_0_0_x_xz_0_xz = buffer_1000_pdsd[122];

    auto g_y_0_0_0_x_xz_0_yy = buffer_1000_pdsd[123];

    auto g_y_0_0_0_x_xz_0_yz = buffer_1000_pdsd[124];

    auto g_y_0_0_0_x_xz_0_zz = buffer_1000_pdsd[125];

    auto g_y_0_0_0_x_yy_0_xx = buffer_1000_pdsd[126];

    auto g_y_0_0_0_x_yy_0_xy = buffer_1000_pdsd[127];

    auto g_y_0_0_0_x_yy_0_xz = buffer_1000_pdsd[128];

    auto g_y_0_0_0_x_yy_0_yy = buffer_1000_pdsd[129];

    auto g_y_0_0_0_x_yy_0_yz = buffer_1000_pdsd[130];

    auto g_y_0_0_0_x_yy_0_zz = buffer_1000_pdsd[131];

    auto g_y_0_0_0_x_yz_0_xx = buffer_1000_pdsd[132];

    auto g_y_0_0_0_x_yz_0_xy = buffer_1000_pdsd[133];

    auto g_y_0_0_0_x_yz_0_xz = buffer_1000_pdsd[134];

    auto g_y_0_0_0_x_yz_0_yy = buffer_1000_pdsd[135];

    auto g_y_0_0_0_x_yz_0_yz = buffer_1000_pdsd[136];

    auto g_y_0_0_0_x_yz_0_zz = buffer_1000_pdsd[137];

    auto g_y_0_0_0_x_zz_0_xx = buffer_1000_pdsd[138];

    auto g_y_0_0_0_x_zz_0_xy = buffer_1000_pdsd[139];

    auto g_y_0_0_0_x_zz_0_xz = buffer_1000_pdsd[140];

    auto g_y_0_0_0_x_zz_0_yy = buffer_1000_pdsd[141];

    auto g_y_0_0_0_x_zz_0_yz = buffer_1000_pdsd[142];

    auto g_y_0_0_0_x_zz_0_zz = buffer_1000_pdsd[143];

    auto g_y_0_0_0_y_xx_0_xx = buffer_1000_pdsd[144];

    auto g_y_0_0_0_y_xx_0_xy = buffer_1000_pdsd[145];

    auto g_y_0_0_0_y_xx_0_xz = buffer_1000_pdsd[146];

    auto g_y_0_0_0_y_xx_0_yy = buffer_1000_pdsd[147];

    auto g_y_0_0_0_y_xx_0_yz = buffer_1000_pdsd[148];

    auto g_y_0_0_0_y_xx_0_zz = buffer_1000_pdsd[149];

    auto g_y_0_0_0_y_xy_0_xx = buffer_1000_pdsd[150];

    auto g_y_0_0_0_y_xy_0_xy = buffer_1000_pdsd[151];

    auto g_y_0_0_0_y_xy_0_xz = buffer_1000_pdsd[152];

    auto g_y_0_0_0_y_xy_0_yy = buffer_1000_pdsd[153];

    auto g_y_0_0_0_y_xy_0_yz = buffer_1000_pdsd[154];

    auto g_y_0_0_0_y_xy_0_zz = buffer_1000_pdsd[155];

    auto g_y_0_0_0_y_xz_0_xx = buffer_1000_pdsd[156];

    auto g_y_0_0_0_y_xz_0_xy = buffer_1000_pdsd[157];

    auto g_y_0_0_0_y_xz_0_xz = buffer_1000_pdsd[158];

    auto g_y_0_0_0_y_xz_0_yy = buffer_1000_pdsd[159];

    auto g_y_0_0_0_y_xz_0_yz = buffer_1000_pdsd[160];

    auto g_y_0_0_0_y_xz_0_zz = buffer_1000_pdsd[161];

    auto g_y_0_0_0_y_yy_0_xx = buffer_1000_pdsd[162];

    auto g_y_0_0_0_y_yy_0_xy = buffer_1000_pdsd[163];

    auto g_y_0_0_0_y_yy_0_xz = buffer_1000_pdsd[164];

    auto g_y_0_0_0_y_yy_0_yy = buffer_1000_pdsd[165];

    auto g_y_0_0_0_y_yy_0_yz = buffer_1000_pdsd[166];

    auto g_y_0_0_0_y_yy_0_zz = buffer_1000_pdsd[167];

    auto g_y_0_0_0_y_yz_0_xx = buffer_1000_pdsd[168];

    auto g_y_0_0_0_y_yz_0_xy = buffer_1000_pdsd[169];

    auto g_y_0_0_0_y_yz_0_xz = buffer_1000_pdsd[170];

    auto g_y_0_0_0_y_yz_0_yy = buffer_1000_pdsd[171];

    auto g_y_0_0_0_y_yz_0_yz = buffer_1000_pdsd[172];

    auto g_y_0_0_0_y_yz_0_zz = buffer_1000_pdsd[173];

    auto g_y_0_0_0_y_zz_0_xx = buffer_1000_pdsd[174];

    auto g_y_0_0_0_y_zz_0_xy = buffer_1000_pdsd[175];

    auto g_y_0_0_0_y_zz_0_xz = buffer_1000_pdsd[176];

    auto g_y_0_0_0_y_zz_0_yy = buffer_1000_pdsd[177];

    auto g_y_0_0_0_y_zz_0_yz = buffer_1000_pdsd[178];

    auto g_y_0_0_0_y_zz_0_zz = buffer_1000_pdsd[179];

    auto g_y_0_0_0_z_xx_0_xx = buffer_1000_pdsd[180];

    auto g_y_0_0_0_z_xx_0_xy = buffer_1000_pdsd[181];

    auto g_y_0_0_0_z_xx_0_xz = buffer_1000_pdsd[182];

    auto g_y_0_0_0_z_xx_0_yy = buffer_1000_pdsd[183];

    auto g_y_0_0_0_z_xx_0_yz = buffer_1000_pdsd[184];

    auto g_y_0_0_0_z_xx_0_zz = buffer_1000_pdsd[185];

    auto g_y_0_0_0_z_xy_0_xx = buffer_1000_pdsd[186];

    auto g_y_0_0_0_z_xy_0_xy = buffer_1000_pdsd[187];

    auto g_y_0_0_0_z_xy_0_xz = buffer_1000_pdsd[188];

    auto g_y_0_0_0_z_xy_0_yy = buffer_1000_pdsd[189];

    auto g_y_0_0_0_z_xy_0_yz = buffer_1000_pdsd[190];

    auto g_y_0_0_0_z_xy_0_zz = buffer_1000_pdsd[191];

    auto g_y_0_0_0_z_xz_0_xx = buffer_1000_pdsd[192];

    auto g_y_0_0_0_z_xz_0_xy = buffer_1000_pdsd[193];

    auto g_y_0_0_0_z_xz_0_xz = buffer_1000_pdsd[194];

    auto g_y_0_0_0_z_xz_0_yy = buffer_1000_pdsd[195];

    auto g_y_0_0_0_z_xz_0_yz = buffer_1000_pdsd[196];

    auto g_y_0_0_0_z_xz_0_zz = buffer_1000_pdsd[197];

    auto g_y_0_0_0_z_yy_0_xx = buffer_1000_pdsd[198];

    auto g_y_0_0_0_z_yy_0_xy = buffer_1000_pdsd[199];

    auto g_y_0_0_0_z_yy_0_xz = buffer_1000_pdsd[200];

    auto g_y_0_0_0_z_yy_0_yy = buffer_1000_pdsd[201];

    auto g_y_0_0_0_z_yy_0_yz = buffer_1000_pdsd[202];

    auto g_y_0_0_0_z_yy_0_zz = buffer_1000_pdsd[203];

    auto g_y_0_0_0_z_yz_0_xx = buffer_1000_pdsd[204];

    auto g_y_0_0_0_z_yz_0_xy = buffer_1000_pdsd[205];

    auto g_y_0_0_0_z_yz_0_xz = buffer_1000_pdsd[206];

    auto g_y_0_0_0_z_yz_0_yy = buffer_1000_pdsd[207];

    auto g_y_0_0_0_z_yz_0_yz = buffer_1000_pdsd[208];

    auto g_y_0_0_0_z_yz_0_zz = buffer_1000_pdsd[209];

    auto g_y_0_0_0_z_zz_0_xx = buffer_1000_pdsd[210];

    auto g_y_0_0_0_z_zz_0_xy = buffer_1000_pdsd[211];

    auto g_y_0_0_0_z_zz_0_xz = buffer_1000_pdsd[212];

    auto g_y_0_0_0_z_zz_0_yy = buffer_1000_pdsd[213];

    auto g_y_0_0_0_z_zz_0_yz = buffer_1000_pdsd[214];

    auto g_y_0_0_0_z_zz_0_zz = buffer_1000_pdsd[215];

    auto g_z_0_0_0_x_xx_0_xx = buffer_1000_pdsd[216];

    auto g_z_0_0_0_x_xx_0_xy = buffer_1000_pdsd[217];

    auto g_z_0_0_0_x_xx_0_xz = buffer_1000_pdsd[218];

    auto g_z_0_0_0_x_xx_0_yy = buffer_1000_pdsd[219];

    auto g_z_0_0_0_x_xx_0_yz = buffer_1000_pdsd[220];

    auto g_z_0_0_0_x_xx_0_zz = buffer_1000_pdsd[221];

    auto g_z_0_0_0_x_xy_0_xx = buffer_1000_pdsd[222];

    auto g_z_0_0_0_x_xy_0_xy = buffer_1000_pdsd[223];

    auto g_z_0_0_0_x_xy_0_xz = buffer_1000_pdsd[224];

    auto g_z_0_0_0_x_xy_0_yy = buffer_1000_pdsd[225];

    auto g_z_0_0_0_x_xy_0_yz = buffer_1000_pdsd[226];

    auto g_z_0_0_0_x_xy_0_zz = buffer_1000_pdsd[227];

    auto g_z_0_0_0_x_xz_0_xx = buffer_1000_pdsd[228];

    auto g_z_0_0_0_x_xz_0_xy = buffer_1000_pdsd[229];

    auto g_z_0_0_0_x_xz_0_xz = buffer_1000_pdsd[230];

    auto g_z_0_0_0_x_xz_0_yy = buffer_1000_pdsd[231];

    auto g_z_0_0_0_x_xz_0_yz = buffer_1000_pdsd[232];

    auto g_z_0_0_0_x_xz_0_zz = buffer_1000_pdsd[233];

    auto g_z_0_0_0_x_yy_0_xx = buffer_1000_pdsd[234];

    auto g_z_0_0_0_x_yy_0_xy = buffer_1000_pdsd[235];

    auto g_z_0_0_0_x_yy_0_xz = buffer_1000_pdsd[236];

    auto g_z_0_0_0_x_yy_0_yy = buffer_1000_pdsd[237];

    auto g_z_0_0_0_x_yy_0_yz = buffer_1000_pdsd[238];

    auto g_z_0_0_0_x_yy_0_zz = buffer_1000_pdsd[239];

    auto g_z_0_0_0_x_yz_0_xx = buffer_1000_pdsd[240];

    auto g_z_0_0_0_x_yz_0_xy = buffer_1000_pdsd[241];

    auto g_z_0_0_0_x_yz_0_xz = buffer_1000_pdsd[242];

    auto g_z_0_0_0_x_yz_0_yy = buffer_1000_pdsd[243];

    auto g_z_0_0_0_x_yz_0_yz = buffer_1000_pdsd[244];

    auto g_z_0_0_0_x_yz_0_zz = buffer_1000_pdsd[245];

    auto g_z_0_0_0_x_zz_0_xx = buffer_1000_pdsd[246];

    auto g_z_0_0_0_x_zz_0_xy = buffer_1000_pdsd[247];

    auto g_z_0_0_0_x_zz_0_xz = buffer_1000_pdsd[248];

    auto g_z_0_0_0_x_zz_0_yy = buffer_1000_pdsd[249];

    auto g_z_0_0_0_x_zz_0_yz = buffer_1000_pdsd[250];

    auto g_z_0_0_0_x_zz_0_zz = buffer_1000_pdsd[251];

    auto g_z_0_0_0_y_xx_0_xx = buffer_1000_pdsd[252];

    auto g_z_0_0_0_y_xx_0_xy = buffer_1000_pdsd[253];

    auto g_z_0_0_0_y_xx_0_xz = buffer_1000_pdsd[254];

    auto g_z_0_0_0_y_xx_0_yy = buffer_1000_pdsd[255];

    auto g_z_0_0_0_y_xx_0_yz = buffer_1000_pdsd[256];

    auto g_z_0_0_0_y_xx_0_zz = buffer_1000_pdsd[257];

    auto g_z_0_0_0_y_xy_0_xx = buffer_1000_pdsd[258];

    auto g_z_0_0_0_y_xy_0_xy = buffer_1000_pdsd[259];

    auto g_z_0_0_0_y_xy_0_xz = buffer_1000_pdsd[260];

    auto g_z_0_0_0_y_xy_0_yy = buffer_1000_pdsd[261];

    auto g_z_0_0_0_y_xy_0_yz = buffer_1000_pdsd[262];

    auto g_z_0_0_0_y_xy_0_zz = buffer_1000_pdsd[263];

    auto g_z_0_0_0_y_xz_0_xx = buffer_1000_pdsd[264];

    auto g_z_0_0_0_y_xz_0_xy = buffer_1000_pdsd[265];

    auto g_z_0_0_0_y_xz_0_xz = buffer_1000_pdsd[266];

    auto g_z_0_0_0_y_xz_0_yy = buffer_1000_pdsd[267];

    auto g_z_0_0_0_y_xz_0_yz = buffer_1000_pdsd[268];

    auto g_z_0_0_0_y_xz_0_zz = buffer_1000_pdsd[269];

    auto g_z_0_0_0_y_yy_0_xx = buffer_1000_pdsd[270];

    auto g_z_0_0_0_y_yy_0_xy = buffer_1000_pdsd[271];

    auto g_z_0_0_0_y_yy_0_xz = buffer_1000_pdsd[272];

    auto g_z_0_0_0_y_yy_0_yy = buffer_1000_pdsd[273];

    auto g_z_0_0_0_y_yy_0_yz = buffer_1000_pdsd[274];

    auto g_z_0_0_0_y_yy_0_zz = buffer_1000_pdsd[275];

    auto g_z_0_0_0_y_yz_0_xx = buffer_1000_pdsd[276];

    auto g_z_0_0_0_y_yz_0_xy = buffer_1000_pdsd[277];

    auto g_z_0_0_0_y_yz_0_xz = buffer_1000_pdsd[278];

    auto g_z_0_0_0_y_yz_0_yy = buffer_1000_pdsd[279];

    auto g_z_0_0_0_y_yz_0_yz = buffer_1000_pdsd[280];

    auto g_z_0_0_0_y_yz_0_zz = buffer_1000_pdsd[281];

    auto g_z_0_0_0_y_zz_0_xx = buffer_1000_pdsd[282];

    auto g_z_0_0_0_y_zz_0_xy = buffer_1000_pdsd[283];

    auto g_z_0_0_0_y_zz_0_xz = buffer_1000_pdsd[284];

    auto g_z_0_0_0_y_zz_0_yy = buffer_1000_pdsd[285];

    auto g_z_0_0_0_y_zz_0_yz = buffer_1000_pdsd[286];

    auto g_z_0_0_0_y_zz_0_zz = buffer_1000_pdsd[287];

    auto g_z_0_0_0_z_xx_0_xx = buffer_1000_pdsd[288];

    auto g_z_0_0_0_z_xx_0_xy = buffer_1000_pdsd[289];

    auto g_z_0_0_0_z_xx_0_xz = buffer_1000_pdsd[290];

    auto g_z_0_0_0_z_xx_0_yy = buffer_1000_pdsd[291];

    auto g_z_0_0_0_z_xx_0_yz = buffer_1000_pdsd[292];

    auto g_z_0_0_0_z_xx_0_zz = buffer_1000_pdsd[293];

    auto g_z_0_0_0_z_xy_0_xx = buffer_1000_pdsd[294];

    auto g_z_0_0_0_z_xy_0_xy = buffer_1000_pdsd[295];

    auto g_z_0_0_0_z_xy_0_xz = buffer_1000_pdsd[296];

    auto g_z_0_0_0_z_xy_0_yy = buffer_1000_pdsd[297];

    auto g_z_0_0_0_z_xy_0_yz = buffer_1000_pdsd[298];

    auto g_z_0_0_0_z_xy_0_zz = buffer_1000_pdsd[299];

    auto g_z_0_0_0_z_xz_0_xx = buffer_1000_pdsd[300];

    auto g_z_0_0_0_z_xz_0_xy = buffer_1000_pdsd[301];

    auto g_z_0_0_0_z_xz_0_xz = buffer_1000_pdsd[302];

    auto g_z_0_0_0_z_xz_0_yy = buffer_1000_pdsd[303];

    auto g_z_0_0_0_z_xz_0_yz = buffer_1000_pdsd[304];

    auto g_z_0_0_0_z_xz_0_zz = buffer_1000_pdsd[305];

    auto g_z_0_0_0_z_yy_0_xx = buffer_1000_pdsd[306];

    auto g_z_0_0_0_z_yy_0_xy = buffer_1000_pdsd[307];

    auto g_z_0_0_0_z_yy_0_xz = buffer_1000_pdsd[308];

    auto g_z_0_0_0_z_yy_0_yy = buffer_1000_pdsd[309];

    auto g_z_0_0_0_z_yy_0_yz = buffer_1000_pdsd[310];

    auto g_z_0_0_0_z_yy_0_zz = buffer_1000_pdsd[311];

    auto g_z_0_0_0_z_yz_0_xx = buffer_1000_pdsd[312];

    auto g_z_0_0_0_z_yz_0_xy = buffer_1000_pdsd[313];

    auto g_z_0_0_0_z_yz_0_xz = buffer_1000_pdsd[314];

    auto g_z_0_0_0_z_yz_0_yy = buffer_1000_pdsd[315];

    auto g_z_0_0_0_z_yz_0_yz = buffer_1000_pdsd[316];

    auto g_z_0_0_0_z_yz_0_zz = buffer_1000_pdsd[317];

    auto g_z_0_0_0_z_zz_0_xx = buffer_1000_pdsd[318];

    auto g_z_0_0_0_z_zz_0_xy = buffer_1000_pdsd[319];

    auto g_z_0_0_0_z_zz_0_xz = buffer_1000_pdsd[320];

    auto g_z_0_0_0_z_zz_0_yy = buffer_1000_pdsd[321];

    auto g_z_0_0_0_z_zz_0_yz = buffer_1000_pdsd[322];

    auto g_z_0_0_0_z_zz_0_zz = buffer_1000_pdsd[323];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_xx_0_xx, g_0_xx_0_xy, g_0_xx_0_xz, g_0_xx_0_yy, g_0_xx_0_yz, g_0_xx_0_zz, g_x_0_0_0_x_xx_0_xx, g_x_0_0_0_x_xx_0_xy, g_x_0_0_0_x_xx_0_xz, g_x_0_0_0_x_xx_0_yy, g_x_0_0_0_x_xx_0_yz, g_x_0_0_0_x_xx_0_zz, g_xx_xx_0_xx, g_xx_xx_0_xy, g_xx_xx_0_xz, g_xx_xx_0_yy, g_xx_xx_0_yz, g_xx_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xx_0_xx[i] = -g_0_xx_0_xx[i] + 2.0 * g_xx_xx_0_xx[i] * a_exp;

        g_x_0_0_0_x_xx_0_xy[i] = -g_0_xx_0_xy[i] + 2.0 * g_xx_xx_0_xy[i] * a_exp;

        g_x_0_0_0_x_xx_0_xz[i] = -g_0_xx_0_xz[i] + 2.0 * g_xx_xx_0_xz[i] * a_exp;

        g_x_0_0_0_x_xx_0_yy[i] = -g_0_xx_0_yy[i] + 2.0 * g_xx_xx_0_yy[i] * a_exp;

        g_x_0_0_0_x_xx_0_yz[i] = -g_0_xx_0_yz[i] + 2.0 * g_xx_xx_0_yz[i] * a_exp;

        g_x_0_0_0_x_xx_0_zz[i] = -g_0_xx_0_zz[i] + 2.0 * g_xx_xx_0_zz[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_x_0_0_0_x_xy_0_xx, g_x_0_0_0_x_xy_0_xy, g_x_0_0_0_x_xy_0_xz, g_x_0_0_0_x_xy_0_yy, g_x_0_0_0_x_xy_0_yz, g_x_0_0_0_x_xy_0_zz, g_xx_xy_0_xx, g_xx_xy_0_xy, g_xx_xy_0_xz, g_xx_xy_0_yy, g_xx_xy_0_yz, g_xx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xy_0_xx[i] = -g_0_xy_0_xx[i] + 2.0 * g_xx_xy_0_xx[i] * a_exp;

        g_x_0_0_0_x_xy_0_xy[i] = -g_0_xy_0_xy[i] + 2.0 * g_xx_xy_0_xy[i] * a_exp;

        g_x_0_0_0_x_xy_0_xz[i] = -g_0_xy_0_xz[i] + 2.0 * g_xx_xy_0_xz[i] * a_exp;

        g_x_0_0_0_x_xy_0_yy[i] = -g_0_xy_0_yy[i] + 2.0 * g_xx_xy_0_yy[i] * a_exp;

        g_x_0_0_0_x_xy_0_yz[i] = -g_0_xy_0_yz[i] + 2.0 * g_xx_xy_0_yz[i] * a_exp;

        g_x_0_0_0_x_xy_0_zz[i] = -g_0_xy_0_zz[i] + 2.0 * g_xx_xy_0_zz[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_x_0_0_0_x_xz_0_xx, g_x_0_0_0_x_xz_0_xy, g_x_0_0_0_x_xz_0_xz, g_x_0_0_0_x_xz_0_yy, g_x_0_0_0_x_xz_0_yz, g_x_0_0_0_x_xz_0_zz, g_xx_xz_0_xx, g_xx_xz_0_xy, g_xx_xz_0_xz, g_xx_xz_0_yy, g_xx_xz_0_yz, g_xx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xz_0_xx[i] = -g_0_xz_0_xx[i] + 2.0 * g_xx_xz_0_xx[i] * a_exp;

        g_x_0_0_0_x_xz_0_xy[i] = -g_0_xz_0_xy[i] + 2.0 * g_xx_xz_0_xy[i] * a_exp;

        g_x_0_0_0_x_xz_0_xz[i] = -g_0_xz_0_xz[i] + 2.0 * g_xx_xz_0_xz[i] * a_exp;

        g_x_0_0_0_x_xz_0_yy[i] = -g_0_xz_0_yy[i] + 2.0 * g_xx_xz_0_yy[i] * a_exp;

        g_x_0_0_0_x_xz_0_yz[i] = -g_0_xz_0_yz[i] + 2.0 * g_xx_xz_0_yz[i] * a_exp;

        g_x_0_0_0_x_xz_0_zz[i] = -g_0_xz_0_zz[i] + 2.0 * g_xx_xz_0_zz[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_0_yy_0_xx, g_0_yy_0_xy, g_0_yy_0_xz, g_0_yy_0_yy, g_0_yy_0_yz, g_0_yy_0_zz, g_x_0_0_0_x_yy_0_xx, g_x_0_0_0_x_yy_0_xy, g_x_0_0_0_x_yy_0_xz, g_x_0_0_0_x_yy_0_yy, g_x_0_0_0_x_yy_0_yz, g_x_0_0_0_x_yy_0_zz, g_xx_yy_0_xx, g_xx_yy_0_xy, g_xx_yy_0_xz, g_xx_yy_0_yy, g_xx_yy_0_yz, g_xx_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yy_0_xx[i] = -g_0_yy_0_xx[i] + 2.0 * g_xx_yy_0_xx[i] * a_exp;

        g_x_0_0_0_x_yy_0_xy[i] = -g_0_yy_0_xy[i] + 2.0 * g_xx_yy_0_xy[i] * a_exp;

        g_x_0_0_0_x_yy_0_xz[i] = -g_0_yy_0_xz[i] + 2.0 * g_xx_yy_0_xz[i] * a_exp;

        g_x_0_0_0_x_yy_0_yy[i] = -g_0_yy_0_yy[i] + 2.0 * g_xx_yy_0_yy[i] * a_exp;

        g_x_0_0_0_x_yy_0_yz[i] = -g_0_yy_0_yz[i] + 2.0 * g_xx_yy_0_yz[i] * a_exp;

        g_x_0_0_0_x_yy_0_zz[i] = -g_0_yy_0_zz[i] + 2.0 * g_xx_yy_0_zz[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_x_0_0_0_x_yz_0_xx, g_x_0_0_0_x_yz_0_xy, g_x_0_0_0_x_yz_0_xz, g_x_0_0_0_x_yz_0_yy, g_x_0_0_0_x_yz_0_yz, g_x_0_0_0_x_yz_0_zz, g_xx_yz_0_xx, g_xx_yz_0_xy, g_xx_yz_0_xz, g_xx_yz_0_yy, g_xx_yz_0_yz, g_xx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yz_0_xx[i] = -g_0_yz_0_xx[i] + 2.0 * g_xx_yz_0_xx[i] * a_exp;

        g_x_0_0_0_x_yz_0_xy[i] = -g_0_yz_0_xy[i] + 2.0 * g_xx_yz_0_xy[i] * a_exp;

        g_x_0_0_0_x_yz_0_xz[i] = -g_0_yz_0_xz[i] + 2.0 * g_xx_yz_0_xz[i] * a_exp;

        g_x_0_0_0_x_yz_0_yy[i] = -g_0_yz_0_yy[i] + 2.0 * g_xx_yz_0_yy[i] * a_exp;

        g_x_0_0_0_x_yz_0_yz[i] = -g_0_yz_0_yz[i] + 2.0 * g_xx_yz_0_yz[i] * a_exp;

        g_x_0_0_0_x_yz_0_zz[i] = -g_0_yz_0_zz[i] + 2.0 * g_xx_yz_0_zz[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_0_zz_0_xx, g_0_zz_0_xy, g_0_zz_0_xz, g_0_zz_0_yy, g_0_zz_0_yz, g_0_zz_0_zz, g_x_0_0_0_x_zz_0_xx, g_x_0_0_0_x_zz_0_xy, g_x_0_0_0_x_zz_0_xz, g_x_0_0_0_x_zz_0_yy, g_x_0_0_0_x_zz_0_yz, g_x_0_0_0_x_zz_0_zz, g_xx_zz_0_xx, g_xx_zz_0_xy, g_xx_zz_0_xz, g_xx_zz_0_yy, g_xx_zz_0_yz, g_xx_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_zz_0_xx[i] = -g_0_zz_0_xx[i] + 2.0 * g_xx_zz_0_xx[i] * a_exp;

        g_x_0_0_0_x_zz_0_xy[i] = -g_0_zz_0_xy[i] + 2.0 * g_xx_zz_0_xy[i] * a_exp;

        g_x_0_0_0_x_zz_0_xz[i] = -g_0_zz_0_xz[i] + 2.0 * g_xx_zz_0_xz[i] * a_exp;

        g_x_0_0_0_x_zz_0_yy[i] = -g_0_zz_0_yy[i] + 2.0 * g_xx_zz_0_yy[i] * a_exp;

        g_x_0_0_0_x_zz_0_yz[i] = -g_0_zz_0_yz[i] + 2.0 * g_xx_zz_0_yz[i] * a_exp;

        g_x_0_0_0_x_zz_0_zz[i] = -g_0_zz_0_zz[i] + 2.0 * g_xx_zz_0_zz[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_0_0_0_y_xx_0_xx, g_x_0_0_0_y_xx_0_xy, g_x_0_0_0_y_xx_0_xz, g_x_0_0_0_y_xx_0_yy, g_x_0_0_0_y_xx_0_yz, g_x_0_0_0_y_xx_0_zz, g_xy_xx_0_xx, g_xy_xx_0_xy, g_xy_xx_0_xz, g_xy_xx_0_yy, g_xy_xx_0_yz, g_xy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xx_0_xx[i] = 2.0 * g_xy_xx_0_xx[i] * a_exp;

        g_x_0_0_0_y_xx_0_xy[i] = 2.0 * g_xy_xx_0_xy[i] * a_exp;

        g_x_0_0_0_y_xx_0_xz[i] = 2.0 * g_xy_xx_0_xz[i] * a_exp;

        g_x_0_0_0_y_xx_0_yy[i] = 2.0 * g_xy_xx_0_yy[i] * a_exp;

        g_x_0_0_0_y_xx_0_yz[i] = 2.0 * g_xy_xx_0_yz[i] * a_exp;

        g_x_0_0_0_y_xx_0_zz[i] = 2.0 * g_xy_xx_0_zz[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_0_0_0_y_xy_0_xx, g_x_0_0_0_y_xy_0_xy, g_x_0_0_0_y_xy_0_xz, g_x_0_0_0_y_xy_0_yy, g_x_0_0_0_y_xy_0_yz, g_x_0_0_0_y_xy_0_zz, g_xy_xy_0_xx, g_xy_xy_0_xy, g_xy_xy_0_xz, g_xy_xy_0_yy, g_xy_xy_0_yz, g_xy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xy_0_xx[i] = 2.0 * g_xy_xy_0_xx[i] * a_exp;

        g_x_0_0_0_y_xy_0_xy[i] = 2.0 * g_xy_xy_0_xy[i] * a_exp;

        g_x_0_0_0_y_xy_0_xz[i] = 2.0 * g_xy_xy_0_xz[i] * a_exp;

        g_x_0_0_0_y_xy_0_yy[i] = 2.0 * g_xy_xy_0_yy[i] * a_exp;

        g_x_0_0_0_y_xy_0_yz[i] = 2.0 * g_xy_xy_0_yz[i] * a_exp;

        g_x_0_0_0_y_xy_0_zz[i] = 2.0 * g_xy_xy_0_zz[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_0_0_y_xz_0_xx, g_x_0_0_0_y_xz_0_xy, g_x_0_0_0_y_xz_0_xz, g_x_0_0_0_y_xz_0_yy, g_x_0_0_0_y_xz_0_yz, g_x_0_0_0_y_xz_0_zz, g_xy_xz_0_xx, g_xy_xz_0_xy, g_xy_xz_0_xz, g_xy_xz_0_yy, g_xy_xz_0_yz, g_xy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xz_0_xx[i] = 2.0 * g_xy_xz_0_xx[i] * a_exp;

        g_x_0_0_0_y_xz_0_xy[i] = 2.0 * g_xy_xz_0_xy[i] * a_exp;

        g_x_0_0_0_y_xz_0_xz[i] = 2.0 * g_xy_xz_0_xz[i] * a_exp;

        g_x_0_0_0_y_xz_0_yy[i] = 2.0 * g_xy_xz_0_yy[i] * a_exp;

        g_x_0_0_0_y_xz_0_yz[i] = 2.0 * g_xy_xz_0_yz[i] * a_exp;

        g_x_0_0_0_y_xz_0_zz[i] = 2.0 * g_xy_xz_0_zz[i] * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_0_0_0_y_yy_0_xx, g_x_0_0_0_y_yy_0_xy, g_x_0_0_0_y_yy_0_xz, g_x_0_0_0_y_yy_0_yy, g_x_0_0_0_y_yy_0_yz, g_x_0_0_0_y_yy_0_zz, g_xy_yy_0_xx, g_xy_yy_0_xy, g_xy_yy_0_xz, g_xy_yy_0_yy, g_xy_yy_0_yz, g_xy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yy_0_xx[i] = 2.0 * g_xy_yy_0_xx[i] * a_exp;

        g_x_0_0_0_y_yy_0_xy[i] = 2.0 * g_xy_yy_0_xy[i] * a_exp;

        g_x_0_0_0_y_yy_0_xz[i] = 2.0 * g_xy_yy_0_xz[i] * a_exp;

        g_x_0_0_0_y_yy_0_yy[i] = 2.0 * g_xy_yy_0_yy[i] * a_exp;

        g_x_0_0_0_y_yy_0_yz[i] = 2.0 * g_xy_yy_0_yz[i] * a_exp;

        g_x_0_0_0_y_yy_0_zz[i] = 2.0 * g_xy_yy_0_zz[i] * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_x_0_0_0_y_yz_0_xx, g_x_0_0_0_y_yz_0_xy, g_x_0_0_0_y_yz_0_xz, g_x_0_0_0_y_yz_0_yy, g_x_0_0_0_y_yz_0_yz, g_x_0_0_0_y_yz_0_zz, g_xy_yz_0_xx, g_xy_yz_0_xy, g_xy_yz_0_xz, g_xy_yz_0_yy, g_xy_yz_0_yz, g_xy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yz_0_xx[i] = 2.0 * g_xy_yz_0_xx[i] * a_exp;

        g_x_0_0_0_y_yz_0_xy[i] = 2.0 * g_xy_yz_0_xy[i] * a_exp;

        g_x_0_0_0_y_yz_0_xz[i] = 2.0 * g_xy_yz_0_xz[i] * a_exp;

        g_x_0_0_0_y_yz_0_yy[i] = 2.0 * g_xy_yz_0_yy[i] * a_exp;

        g_x_0_0_0_y_yz_0_yz[i] = 2.0 * g_xy_yz_0_yz[i] * a_exp;

        g_x_0_0_0_y_yz_0_zz[i] = 2.0 * g_xy_yz_0_zz[i] * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_x_0_0_0_y_zz_0_xx, g_x_0_0_0_y_zz_0_xy, g_x_0_0_0_y_zz_0_xz, g_x_0_0_0_y_zz_0_yy, g_x_0_0_0_y_zz_0_yz, g_x_0_0_0_y_zz_0_zz, g_xy_zz_0_xx, g_xy_zz_0_xy, g_xy_zz_0_xz, g_xy_zz_0_yy, g_xy_zz_0_yz, g_xy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_zz_0_xx[i] = 2.0 * g_xy_zz_0_xx[i] * a_exp;

        g_x_0_0_0_y_zz_0_xy[i] = 2.0 * g_xy_zz_0_xy[i] * a_exp;

        g_x_0_0_0_y_zz_0_xz[i] = 2.0 * g_xy_zz_0_xz[i] * a_exp;

        g_x_0_0_0_y_zz_0_yy[i] = 2.0 * g_xy_zz_0_yy[i] * a_exp;

        g_x_0_0_0_y_zz_0_yz[i] = 2.0 * g_xy_zz_0_yz[i] * a_exp;

        g_x_0_0_0_y_zz_0_zz[i] = 2.0 * g_xy_zz_0_zz[i] * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_0_0_0_z_xx_0_xx, g_x_0_0_0_z_xx_0_xy, g_x_0_0_0_z_xx_0_xz, g_x_0_0_0_z_xx_0_yy, g_x_0_0_0_z_xx_0_yz, g_x_0_0_0_z_xx_0_zz, g_xz_xx_0_xx, g_xz_xx_0_xy, g_xz_xx_0_xz, g_xz_xx_0_yy, g_xz_xx_0_yz, g_xz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xx_0_xx[i] = 2.0 * g_xz_xx_0_xx[i] * a_exp;

        g_x_0_0_0_z_xx_0_xy[i] = 2.0 * g_xz_xx_0_xy[i] * a_exp;

        g_x_0_0_0_z_xx_0_xz[i] = 2.0 * g_xz_xx_0_xz[i] * a_exp;

        g_x_0_0_0_z_xx_0_yy[i] = 2.0 * g_xz_xx_0_yy[i] * a_exp;

        g_x_0_0_0_z_xx_0_yz[i] = 2.0 * g_xz_xx_0_yz[i] * a_exp;

        g_x_0_0_0_z_xx_0_zz[i] = 2.0 * g_xz_xx_0_zz[i] * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_0_0_0_z_xy_0_xx, g_x_0_0_0_z_xy_0_xy, g_x_0_0_0_z_xy_0_xz, g_x_0_0_0_z_xy_0_yy, g_x_0_0_0_z_xy_0_yz, g_x_0_0_0_z_xy_0_zz, g_xz_xy_0_xx, g_xz_xy_0_xy, g_xz_xy_0_xz, g_xz_xy_0_yy, g_xz_xy_0_yz, g_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xy_0_xx[i] = 2.0 * g_xz_xy_0_xx[i] * a_exp;

        g_x_0_0_0_z_xy_0_xy[i] = 2.0 * g_xz_xy_0_xy[i] * a_exp;

        g_x_0_0_0_z_xy_0_xz[i] = 2.0 * g_xz_xy_0_xz[i] * a_exp;

        g_x_0_0_0_z_xy_0_yy[i] = 2.0 * g_xz_xy_0_yy[i] * a_exp;

        g_x_0_0_0_z_xy_0_yz[i] = 2.0 * g_xz_xy_0_yz[i] * a_exp;

        g_x_0_0_0_z_xy_0_zz[i] = 2.0 * g_xz_xy_0_zz[i] * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_0_0_0_z_xz_0_xx, g_x_0_0_0_z_xz_0_xy, g_x_0_0_0_z_xz_0_xz, g_x_0_0_0_z_xz_0_yy, g_x_0_0_0_z_xz_0_yz, g_x_0_0_0_z_xz_0_zz, g_xz_xz_0_xx, g_xz_xz_0_xy, g_xz_xz_0_xz, g_xz_xz_0_yy, g_xz_xz_0_yz, g_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xz_0_xx[i] = 2.0 * g_xz_xz_0_xx[i] * a_exp;

        g_x_0_0_0_z_xz_0_xy[i] = 2.0 * g_xz_xz_0_xy[i] * a_exp;

        g_x_0_0_0_z_xz_0_xz[i] = 2.0 * g_xz_xz_0_xz[i] * a_exp;

        g_x_0_0_0_z_xz_0_yy[i] = 2.0 * g_xz_xz_0_yy[i] * a_exp;

        g_x_0_0_0_z_xz_0_yz[i] = 2.0 * g_xz_xz_0_yz[i] * a_exp;

        g_x_0_0_0_z_xz_0_zz[i] = 2.0 * g_xz_xz_0_zz[i] * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_0_0_0_z_yy_0_xx, g_x_0_0_0_z_yy_0_xy, g_x_0_0_0_z_yy_0_xz, g_x_0_0_0_z_yy_0_yy, g_x_0_0_0_z_yy_0_yz, g_x_0_0_0_z_yy_0_zz, g_xz_yy_0_xx, g_xz_yy_0_xy, g_xz_yy_0_xz, g_xz_yy_0_yy, g_xz_yy_0_yz, g_xz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yy_0_xx[i] = 2.0 * g_xz_yy_0_xx[i] * a_exp;

        g_x_0_0_0_z_yy_0_xy[i] = 2.0 * g_xz_yy_0_xy[i] * a_exp;

        g_x_0_0_0_z_yy_0_xz[i] = 2.0 * g_xz_yy_0_xz[i] * a_exp;

        g_x_0_0_0_z_yy_0_yy[i] = 2.0 * g_xz_yy_0_yy[i] * a_exp;

        g_x_0_0_0_z_yy_0_yz[i] = 2.0 * g_xz_yy_0_yz[i] * a_exp;

        g_x_0_0_0_z_yy_0_zz[i] = 2.0 * g_xz_yy_0_zz[i] * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_0_0_0_z_yz_0_xx, g_x_0_0_0_z_yz_0_xy, g_x_0_0_0_z_yz_0_xz, g_x_0_0_0_z_yz_0_yy, g_x_0_0_0_z_yz_0_yz, g_x_0_0_0_z_yz_0_zz, g_xz_yz_0_xx, g_xz_yz_0_xy, g_xz_yz_0_xz, g_xz_yz_0_yy, g_xz_yz_0_yz, g_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yz_0_xx[i] = 2.0 * g_xz_yz_0_xx[i] * a_exp;

        g_x_0_0_0_z_yz_0_xy[i] = 2.0 * g_xz_yz_0_xy[i] * a_exp;

        g_x_0_0_0_z_yz_0_xz[i] = 2.0 * g_xz_yz_0_xz[i] * a_exp;

        g_x_0_0_0_z_yz_0_yy[i] = 2.0 * g_xz_yz_0_yy[i] * a_exp;

        g_x_0_0_0_z_yz_0_yz[i] = 2.0 * g_xz_yz_0_yz[i] * a_exp;

        g_x_0_0_0_z_yz_0_zz[i] = 2.0 * g_xz_yz_0_zz[i] * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_0_0_0_z_zz_0_xx, g_x_0_0_0_z_zz_0_xy, g_x_0_0_0_z_zz_0_xz, g_x_0_0_0_z_zz_0_yy, g_x_0_0_0_z_zz_0_yz, g_x_0_0_0_z_zz_0_zz, g_xz_zz_0_xx, g_xz_zz_0_xy, g_xz_zz_0_xz, g_xz_zz_0_yy, g_xz_zz_0_yz, g_xz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_zz_0_xx[i] = 2.0 * g_xz_zz_0_xx[i] * a_exp;

        g_x_0_0_0_z_zz_0_xy[i] = 2.0 * g_xz_zz_0_xy[i] * a_exp;

        g_x_0_0_0_z_zz_0_xz[i] = 2.0 * g_xz_zz_0_xz[i] * a_exp;

        g_x_0_0_0_z_zz_0_yy[i] = 2.0 * g_xz_zz_0_yy[i] * a_exp;

        g_x_0_0_0_z_zz_0_yz[i] = 2.0 * g_xz_zz_0_yz[i] * a_exp;

        g_x_0_0_0_z_zz_0_zz[i] = 2.0 * g_xz_zz_0_zz[i] * a_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_xy_xx_0_xx, g_xy_xx_0_xy, g_xy_xx_0_xz, g_xy_xx_0_yy, g_xy_xx_0_yz, g_xy_xx_0_zz, g_y_0_0_0_x_xx_0_xx, g_y_0_0_0_x_xx_0_xy, g_y_0_0_0_x_xx_0_xz, g_y_0_0_0_x_xx_0_yy, g_y_0_0_0_x_xx_0_yz, g_y_0_0_0_x_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xx_0_xx[i] = 2.0 * g_xy_xx_0_xx[i] * a_exp;

        g_y_0_0_0_x_xx_0_xy[i] = 2.0 * g_xy_xx_0_xy[i] * a_exp;

        g_y_0_0_0_x_xx_0_xz[i] = 2.0 * g_xy_xx_0_xz[i] * a_exp;

        g_y_0_0_0_x_xx_0_yy[i] = 2.0 * g_xy_xx_0_yy[i] * a_exp;

        g_y_0_0_0_x_xx_0_yz[i] = 2.0 * g_xy_xx_0_yz[i] * a_exp;

        g_y_0_0_0_x_xx_0_zz[i] = 2.0 * g_xy_xx_0_zz[i] * a_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_xy_xy_0_xx, g_xy_xy_0_xy, g_xy_xy_0_xz, g_xy_xy_0_yy, g_xy_xy_0_yz, g_xy_xy_0_zz, g_y_0_0_0_x_xy_0_xx, g_y_0_0_0_x_xy_0_xy, g_y_0_0_0_x_xy_0_xz, g_y_0_0_0_x_xy_0_yy, g_y_0_0_0_x_xy_0_yz, g_y_0_0_0_x_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xy_0_xx[i] = 2.0 * g_xy_xy_0_xx[i] * a_exp;

        g_y_0_0_0_x_xy_0_xy[i] = 2.0 * g_xy_xy_0_xy[i] * a_exp;

        g_y_0_0_0_x_xy_0_xz[i] = 2.0 * g_xy_xy_0_xz[i] * a_exp;

        g_y_0_0_0_x_xy_0_yy[i] = 2.0 * g_xy_xy_0_yy[i] * a_exp;

        g_y_0_0_0_x_xy_0_yz[i] = 2.0 * g_xy_xy_0_yz[i] * a_exp;

        g_y_0_0_0_x_xy_0_zz[i] = 2.0 * g_xy_xy_0_zz[i] * a_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_xy_xz_0_xx, g_xy_xz_0_xy, g_xy_xz_0_xz, g_xy_xz_0_yy, g_xy_xz_0_yz, g_xy_xz_0_zz, g_y_0_0_0_x_xz_0_xx, g_y_0_0_0_x_xz_0_xy, g_y_0_0_0_x_xz_0_xz, g_y_0_0_0_x_xz_0_yy, g_y_0_0_0_x_xz_0_yz, g_y_0_0_0_x_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xz_0_xx[i] = 2.0 * g_xy_xz_0_xx[i] * a_exp;

        g_y_0_0_0_x_xz_0_xy[i] = 2.0 * g_xy_xz_0_xy[i] * a_exp;

        g_y_0_0_0_x_xz_0_xz[i] = 2.0 * g_xy_xz_0_xz[i] * a_exp;

        g_y_0_0_0_x_xz_0_yy[i] = 2.0 * g_xy_xz_0_yy[i] * a_exp;

        g_y_0_0_0_x_xz_0_yz[i] = 2.0 * g_xy_xz_0_yz[i] * a_exp;

        g_y_0_0_0_x_xz_0_zz[i] = 2.0 * g_xy_xz_0_zz[i] * a_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_xy_yy_0_xx, g_xy_yy_0_xy, g_xy_yy_0_xz, g_xy_yy_0_yy, g_xy_yy_0_yz, g_xy_yy_0_zz, g_y_0_0_0_x_yy_0_xx, g_y_0_0_0_x_yy_0_xy, g_y_0_0_0_x_yy_0_xz, g_y_0_0_0_x_yy_0_yy, g_y_0_0_0_x_yy_0_yz, g_y_0_0_0_x_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yy_0_xx[i] = 2.0 * g_xy_yy_0_xx[i] * a_exp;

        g_y_0_0_0_x_yy_0_xy[i] = 2.0 * g_xy_yy_0_xy[i] * a_exp;

        g_y_0_0_0_x_yy_0_xz[i] = 2.0 * g_xy_yy_0_xz[i] * a_exp;

        g_y_0_0_0_x_yy_0_yy[i] = 2.0 * g_xy_yy_0_yy[i] * a_exp;

        g_y_0_0_0_x_yy_0_yz[i] = 2.0 * g_xy_yy_0_yz[i] * a_exp;

        g_y_0_0_0_x_yy_0_zz[i] = 2.0 * g_xy_yy_0_zz[i] * a_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_xy_yz_0_xx, g_xy_yz_0_xy, g_xy_yz_0_xz, g_xy_yz_0_yy, g_xy_yz_0_yz, g_xy_yz_0_zz, g_y_0_0_0_x_yz_0_xx, g_y_0_0_0_x_yz_0_xy, g_y_0_0_0_x_yz_0_xz, g_y_0_0_0_x_yz_0_yy, g_y_0_0_0_x_yz_0_yz, g_y_0_0_0_x_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yz_0_xx[i] = 2.0 * g_xy_yz_0_xx[i] * a_exp;

        g_y_0_0_0_x_yz_0_xy[i] = 2.0 * g_xy_yz_0_xy[i] * a_exp;

        g_y_0_0_0_x_yz_0_xz[i] = 2.0 * g_xy_yz_0_xz[i] * a_exp;

        g_y_0_0_0_x_yz_0_yy[i] = 2.0 * g_xy_yz_0_yy[i] * a_exp;

        g_y_0_0_0_x_yz_0_yz[i] = 2.0 * g_xy_yz_0_yz[i] * a_exp;

        g_y_0_0_0_x_yz_0_zz[i] = 2.0 * g_xy_yz_0_zz[i] * a_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_xy_zz_0_xx, g_xy_zz_0_xy, g_xy_zz_0_xz, g_xy_zz_0_yy, g_xy_zz_0_yz, g_xy_zz_0_zz, g_y_0_0_0_x_zz_0_xx, g_y_0_0_0_x_zz_0_xy, g_y_0_0_0_x_zz_0_xz, g_y_0_0_0_x_zz_0_yy, g_y_0_0_0_x_zz_0_yz, g_y_0_0_0_x_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_zz_0_xx[i] = 2.0 * g_xy_zz_0_xx[i] * a_exp;

        g_y_0_0_0_x_zz_0_xy[i] = 2.0 * g_xy_zz_0_xy[i] * a_exp;

        g_y_0_0_0_x_zz_0_xz[i] = 2.0 * g_xy_zz_0_xz[i] * a_exp;

        g_y_0_0_0_x_zz_0_yy[i] = 2.0 * g_xy_zz_0_yy[i] * a_exp;

        g_y_0_0_0_x_zz_0_yz[i] = 2.0 * g_xy_zz_0_yz[i] * a_exp;

        g_y_0_0_0_x_zz_0_zz[i] = 2.0 * g_xy_zz_0_zz[i] * a_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_0_xx_0_xx, g_0_xx_0_xy, g_0_xx_0_xz, g_0_xx_0_yy, g_0_xx_0_yz, g_0_xx_0_zz, g_y_0_0_0_y_xx_0_xx, g_y_0_0_0_y_xx_0_xy, g_y_0_0_0_y_xx_0_xz, g_y_0_0_0_y_xx_0_yy, g_y_0_0_0_y_xx_0_yz, g_y_0_0_0_y_xx_0_zz, g_yy_xx_0_xx, g_yy_xx_0_xy, g_yy_xx_0_xz, g_yy_xx_0_yy, g_yy_xx_0_yz, g_yy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xx_0_xx[i] = -g_0_xx_0_xx[i] + 2.0 * g_yy_xx_0_xx[i] * a_exp;

        g_y_0_0_0_y_xx_0_xy[i] = -g_0_xx_0_xy[i] + 2.0 * g_yy_xx_0_xy[i] * a_exp;

        g_y_0_0_0_y_xx_0_xz[i] = -g_0_xx_0_xz[i] + 2.0 * g_yy_xx_0_xz[i] * a_exp;

        g_y_0_0_0_y_xx_0_yy[i] = -g_0_xx_0_yy[i] + 2.0 * g_yy_xx_0_yy[i] * a_exp;

        g_y_0_0_0_y_xx_0_yz[i] = -g_0_xx_0_yz[i] + 2.0 * g_yy_xx_0_yz[i] * a_exp;

        g_y_0_0_0_y_xx_0_zz[i] = -g_0_xx_0_zz[i] + 2.0 * g_yy_xx_0_zz[i] * a_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_y_0_0_0_y_xy_0_xx, g_y_0_0_0_y_xy_0_xy, g_y_0_0_0_y_xy_0_xz, g_y_0_0_0_y_xy_0_yy, g_y_0_0_0_y_xy_0_yz, g_y_0_0_0_y_xy_0_zz, g_yy_xy_0_xx, g_yy_xy_0_xy, g_yy_xy_0_xz, g_yy_xy_0_yy, g_yy_xy_0_yz, g_yy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xy_0_xx[i] = -g_0_xy_0_xx[i] + 2.0 * g_yy_xy_0_xx[i] * a_exp;

        g_y_0_0_0_y_xy_0_xy[i] = -g_0_xy_0_xy[i] + 2.0 * g_yy_xy_0_xy[i] * a_exp;

        g_y_0_0_0_y_xy_0_xz[i] = -g_0_xy_0_xz[i] + 2.0 * g_yy_xy_0_xz[i] * a_exp;

        g_y_0_0_0_y_xy_0_yy[i] = -g_0_xy_0_yy[i] + 2.0 * g_yy_xy_0_yy[i] * a_exp;

        g_y_0_0_0_y_xy_0_yz[i] = -g_0_xy_0_yz[i] + 2.0 * g_yy_xy_0_yz[i] * a_exp;

        g_y_0_0_0_y_xy_0_zz[i] = -g_0_xy_0_zz[i] + 2.0 * g_yy_xy_0_zz[i] * a_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_y_0_0_0_y_xz_0_xx, g_y_0_0_0_y_xz_0_xy, g_y_0_0_0_y_xz_0_xz, g_y_0_0_0_y_xz_0_yy, g_y_0_0_0_y_xz_0_yz, g_y_0_0_0_y_xz_0_zz, g_yy_xz_0_xx, g_yy_xz_0_xy, g_yy_xz_0_xz, g_yy_xz_0_yy, g_yy_xz_0_yz, g_yy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xz_0_xx[i] = -g_0_xz_0_xx[i] + 2.0 * g_yy_xz_0_xx[i] * a_exp;

        g_y_0_0_0_y_xz_0_xy[i] = -g_0_xz_0_xy[i] + 2.0 * g_yy_xz_0_xy[i] * a_exp;

        g_y_0_0_0_y_xz_0_xz[i] = -g_0_xz_0_xz[i] + 2.0 * g_yy_xz_0_xz[i] * a_exp;

        g_y_0_0_0_y_xz_0_yy[i] = -g_0_xz_0_yy[i] + 2.0 * g_yy_xz_0_yy[i] * a_exp;

        g_y_0_0_0_y_xz_0_yz[i] = -g_0_xz_0_yz[i] + 2.0 * g_yy_xz_0_yz[i] * a_exp;

        g_y_0_0_0_y_xz_0_zz[i] = -g_0_xz_0_zz[i] + 2.0 * g_yy_xz_0_zz[i] * a_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_0_yy_0_xx, g_0_yy_0_xy, g_0_yy_0_xz, g_0_yy_0_yy, g_0_yy_0_yz, g_0_yy_0_zz, g_y_0_0_0_y_yy_0_xx, g_y_0_0_0_y_yy_0_xy, g_y_0_0_0_y_yy_0_xz, g_y_0_0_0_y_yy_0_yy, g_y_0_0_0_y_yy_0_yz, g_y_0_0_0_y_yy_0_zz, g_yy_yy_0_xx, g_yy_yy_0_xy, g_yy_yy_0_xz, g_yy_yy_0_yy, g_yy_yy_0_yz, g_yy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yy_0_xx[i] = -g_0_yy_0_xx[i] + 2.0 * g_yy_yy_0_xx[i] * a_exp;

        g_y_0_0_0_y_yy_0_xy[i] = -g_0_yy_0_xy[i] + 2.0 * g_yy_yy_0_xy[i] * a_exp;

        g_y_0_0_0_y_yy_0_xz[i] = -g_0_yy_0_xz[i] + 2.0 * g_yy_yy_0_xz[i] * a_exp;

        g_y_0_0_0_y_yy_0_yy[i] = -g_0_yy_0_yy[i] + 2.0 * g_yy_yy_0_yy[i] * a_exp;

        g_y_0_0_0_y_yy_0_yz[i] = -g_0_yy_0_yz[i] + 2.0 * g_yy_yy_0_yz[i] * a_exp;

        g_y_0_0_0_y_yy_0_zz[i] = -g_0_yy_0_zz[i] + 2.0 * g_yy_yy_0_zz[i] * a_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_y_0_0_0_y_yz_0_xx, g_y_0_0_0_y_yz_0_xy, g_y_0_0_0_y_yz_0_xz, g_y_0_0_0_y_yz_0_yy, g_y_0_0_0_y_yz_0_yz, g_y_0_0_0_y_yz_0_zz, g_yy_yz_0_xx, g_yy_yz_0_xy, g_yy_yz_0_xz, g_yy_yz_0_yy, g_yy_yz_0_yz, g_yy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yz_0_xx[i] = -g_0_yz_0_xx[i] + 2.0 * g_yy_yz_0_xx[i] * a_exp;

        g_y_0_0_0_y_yz_0_xy[i] = -g_0_yz_0_xy[i] + 2.0 * g_yy_yz_0_xy[i] * a_exp;

        g_y_0_0_0_y_yz_0_xz[i] = -g_0_yz_0_xz[i] + 2.0 * g_yy_yz_0_xz[i] * a_exp;

        g_y_0_0_0_y_yz_0_yy[i] = -g_0_yz_0_yy[i] + 2.0 * g_yy_yz_0_yy[i] * a_exp;

        g_y_0_0_0_y_yz_0_yz[i] = -g_0_yz_0_yz[i] + 2.0 * g_yy_yz_0_yz[i] * a_exp;

        g_y_0_0_0_y_yz_0_zz[i] = -g_0_yz_0_zz[i] + 2.0 * g_yy_yz_0_zz[i] * a_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_0_zz_0_xx, g_0_zz_0_xy, g_0_zz_0_xz, g_0_zz_0_yy, g_0_zz_0_yz, g_0_zz_0_zz, g_y_0_0_0_y_zz_0_xx, g_y_0_0_0_y_zz_0_xy, g_y_0_0_0_y_zz_0_xz, g_y_0_0_0_y_zz_0_yy, g_y_0_0_0_y_zz_0_yz, g_y_0_0_0_y_zz_0_zz, g_yy_zz_0_xx, g_yy_zz_0_xy, g_yy_zz_0_xz, g_yy_zz_0_yy, g_yy_zz_0_yz, g_yy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_zz_0_xx[i] = -g_0_zz_0_xx[i] + 2.0 * g_yy_zz_0_xx[i] * a_exp;

        g_y_0_0_0_y_zz_0_xy[i] = -g_0_zz_0_xy[i] + 2.0 * g_yy_zz_0_xy[i] * a_exp;

        g_y_0_0_0_y_zz_0_xz[i] = -g_0_zz_0_xz[i] + 2.0 * g_yy_zz_0_xz[i] * a_exp;

        g_y_0_0_0_y_zz_0_yy[i] = -g_0_zz_0_yy[i] + 2.0 * g_yy_zz_0_yy[i] * a_exp;

        g_y_0_0_0_y_zz_0_yz[i] = -g_0_zz_0_yz[i] + 2.0 * g_yy_zz_0_yz[i] * a_exp;

        g_y_0_0_0_y_zz_0_zz[i] = -g_0_zz_0_zz[i] + 2.0 * g_yy_zz_0_zz[i] * a_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_y_0_0_0_z_xx_0_xx, g_y_0_0_0_z_xx_0_xy, g_y_0_0_0_z_xx_0_xz, g_y_0_0_0_z_xx_0_yy, g_y_0_0_0_z_xx_0_yz, g_y_0_0_0_z_xx_0_zz, g_yz_xx_0_xx, g_yz_xx_0_xy, g_yz_xx_0_xz, g_yz_xx_0_yy, g_yz_xx_0_yz, g_yz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xx_0_xx[i] = 2.0 * g_yz_xx_0_xx[i] * a_exp;

        g_y_0_0_0_z_xx_0_xy[i] = 2.0 * g_yz_xx_0_xy[i] * a_exp;

        g_y_0_0_0_z_xx_0_xz[i] = 2.0 * g_yz_xx_0_xz[i] * a_exp;

        g_y_0_0_0_z_xx_0_yy[i] = 2.0 * g_yz_xx_0_yy[i] * a_exp;

        g_y_0_0_0_z_xx_0_yz[i] = 2.0 * g_yz_xx_0_yz[i] * a_exp;

        g_y_0_0_0_z_xx_0_zz[i] = 2.0 * g_yz_xx_0_zz[i] * a_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_y_0_0_0_z_xy_0_xx, g_y_0_0_0_z_xy_0_xy, g_y_0_0_0_z_xy_0_xz, g_y_0_0_0_z_xy_0_yy, g_y_0_0_0_z_xy_0_yz, g_y_0_0_0_z_xy_0_zz, g_yz_xy_0_xx, g_yz_xy_0_xy, g_yz_xy_0_xz, g_yz_xy_0_yy, g_yz_xy_0_yz, g_yz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xy_0_xx[i] = 2.0 * g_yz_xy_0_xx[i] * a_exp;

        g_y_0_0_0_z_xy_0_xy[i] = 2.0 * g_yz_xy_0_xy[i] * a_exp;

        g_y_0_0_0_z_xy_0_xz[i] = 2.0 * g_yz_xy_0_xz[i] * a_exp;

        g_y_0_0_0_z_xy_0_yy[i] = 2.0 * g_yz_xy_0_yy[i] * a_exp;

        g_y_0_0_0_z_xy_0_yz[i] = 2.0 * g_yz_xy_0_yz[i] * a_exp;

        g_y_0_0_0_z_xy_0_zz[i] = 2.0 * g_yz_xy_0_zz[i] * a_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_y_0_0_0_z_xz_0_xx, g_y_0_0_0_z_xz_0_xy, g_y_0_0_0_z_xz_0_xz, g_y_0_0_0_z_xz_0_yy, g_y_0_0_0_z_xz_0_yz, g_y_0_0_0_z_xz_0_zz, g_yz_xz_0_xx, g_yz_xz_0_xy, g_yz_xz_0_xz, g_yz_xz_0_yy, g_yz_xz_0_yz, g_yz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xz_0_xx[i] = 2.0 * g_yz_xz_0_xx[i] * a_exp;

        g_y_0_0_0_z_xz_0_xy[i] = 2.0 * g_yz_xz_0_xy[i] * a_exp;

        g_y_0_0_0_z_xz_0_xz[i] = 2.0 * g_yz_xz_0_xz[i] * a_exp;

        g_y_0_0_0_z_xz_0_yy[i] = 2.0 * g_yz_xz_0_yy[i] * a_exp;

        g_y_0_0_0_z_xz_0_yz[i] = 2.0 * g_yz_xz_0_yz[i] * a_exp;

        g_y_0_0_0_z_xz_0_zz[i] = 2.0 * g_yz_xz_0_zz[i] * a_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_y_0_0_0_z_yy_0_xx, g_y_0_0_0_z_yy_0_xy, g_y_0_0_0_z_yy_0_xz, g_y_0_0_0_z_yy_0_yy, g_y_0_0_0_z_yy_0_yz, g_y_0_0_0_z_yy_0_zz, g_yz_yy_0_xx, g_yz_yy_0_xy, g_yz_yy_0_xz, g_yz_yy_0_yy, g_yz_yy_0_yz, g_yz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yy_0_xx[i] = 2.0 * g_yz_yy_0_xx[i] * a_exp;

        g_y_0_0_0_z_yy_0_xy[i] = 2.0 * g_yz_yy_0_xy[i] * a_exp;

        g_y_0_0_0_z_yy_0_xz[i] = 2.0 * g_yz_yy_0_xz[i] * a_exp;

        g_y_0_0_0_z_yy_0_yy[i] = 2.0 * g_yz_yy_0_yy[i] * a_exp;

        g_y_0_0_0_z_yy_0_yz[i] = 2.0 * g_yz_yy_0_yz[i] * a_exp;

        g_y_0_0_0_z_yy_0_zz[i] = 2.0 * g_yz_yy_0_zz[i] * a_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_y_0_0_0_z_yz_0_xx, g_y_0_0_0_z_yz_0_xy, g_y_0_0_0_z_yz_0_xz, g_y_0_0_0_z_yz_0_yy, g_y_0_0_0_z_yz_0_yz, g_y_0_0_0_z_yz_0_zz, g_yz_yz_0_xx, g_yz_yz_0_xy, g_yz_yz_0_xz, g_yz_yz_0_yy, g_yz_yz_0_yz, g_yz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yz_0_xx[i] = 2.0 * g_yz_yz_0_xx[i] * a_exp;

        g_y_0_0_0_z_yz_0_xy[i] = 2.0 * g_yz_yz_0_xy[i] * a_exp;

        g_y_0_0_0_z_yz_0_xz[i] = 2.0 * g_yz_yz_0_xz[i] * a_exp;

        g_y_0_0_0_z_yz_0_yy[i] = 2.0 * g_yz_yz_0_yy[i] * a_exp;

        g_y_0_0_0_z_yz_0_yz[i] = 2.0 * g_yz_yz_0_yz[i] * a_exp;

        g_y_0_0_0_z_yz_0_zz[i] = 2.0 * g_yz_yz_0_zz[i] * a_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_y_0_0_0_z_zz_0_xx, g_y_0_0_0_z_zz_0_xy, g_y_0_0_0_z_zz_0_xz, g_y_0_0_0_z_zz_0_yy, g_y_0_0_0_z_zz_0_yz, g_y_0_0_0_z_zz_0_zz, g_yz_zz_0_xx, g_yz_zz_0_xy, g_yz_zz_0_xz, g_yz_zz_0_yy, g_yz_zz_0_yz, g_yz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_zz_0_xx[i] = 2.0 * g_yz_zz_0_xx[i] * a_exp;

        g_y_0_0_0_z_zz_0_xy[i] = 2.0 * g_yz_zz_0_xy[i] * a_exp;

        g_y_0_0_0_z_zz_0_xz[i] = 2.0 * g_yz_zz_0_xz[i] * a_exp;

        g_y_0_0_0_z_zz_0_yy[i] = 2.0 * g_yz_zz_0_yy[i] * a_exp;

        g_y_0_0_0_z_zz_0_yz[i] = 2.0 * g_yz_zz_0_yz[i] * a_exp;

        g_y_0_0_0_z_zz_0_zz[i] = 2.0 * g_yz_zz_0_zz[i] * a_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_xz_xx_0_xx, g_xz_xx_0_xy, g_xz_xx_0_xz, g_xz_xx_0_yy, g_xz_xx_0_yz, g_xz_xx_0_zz, g_z_0_0_0_x_xx_0_xx, g_z_0_0_0_x_xx_0_xy, g_z_0_0_0_x_xx_0_xz, g_z_0_0_0_x_xx_0_yy, g_z_0_0_0_x_xx_0_yz, g_z_0_0_0_x_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xx_0_xx[i] = 2.0 * g_xz_xx_0_xx[i] * a_exp;

        g_z_0_0_0_x_xx_0_xy[i] = 2.0 * g_xz_xx_0_xy[i] * a_exp;

        g_z_0_0_0_x_xx_0_xz[i] = 2.0 * g_xz_xx_0_xz[i] * a_exp;

        g_z_0_0_0_x_xx_0_yy[i] = 2.0 * g_xz_xx_0_yy[i] * a_exp;

        g_z_0_0_0_x_xx_0_yz[i] = 2.0 * g_xz_xx_0_yz[i] * a_exp;

        g_z_0_0_0_x_xx_0_zz[i] = 2.0 * g_xz_xx_0_zz[i] * a_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_xz_xy_0_xx, g_xz_xy_0_xy, g_xz_xy_0_xz, g_xz_xy_0_yy, g_xz_xy_0_yz, g_xz_xy_0_zz, g_z_0_0_0_x_xy_0_xx, g_z_0_0_0_x_xy_0_xy, g_z_0_0_0_x_xy_0_xz, g_z_0_0_0_x_xy_0_yy, g_z_0_0_0_x_xy_0_yz, g_z_0_0_0_x_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xy_0_xx[i] = 2.0 * g_xz_xy_0_xx[i] * a_exp;

        g_z_0_0_0_x_xy_0_xy[i] = 2.0 * g_xz_xy_0_xy[i] * a_exp;

        g_z_0_0_0_x_xy_0_xz[i] = 2.0 * g_xz_xy_0_xz[i] * a_exp;

        g_z_0_0_0_x_xy_0_yy[i] = 2.0 * g_xz_xy_0_yy[i] * a_exp;

        g_z_0_0_0_x_xy_0_yz[i] = 2.0 * g_xz_xy_0_yz[i] * a_exp;

        g_z_0_0_0_x_xy_0_zz[i] = 2.0 * g_xz_xy_0_zz[i] * a_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_xz_xz_0_xx, g_xz_xz_0_xy, g_xz_xz_0_xz, g_xz_xz_0_yy, g_xz_xz_0_yz, g_xz_xz_0_zz, g_z_0_0_0_x_xz_0_xx, g_z_0_0_0_x_xz_0_xy, g_z_0_0_0_x_xz_0_xz, g_z_0_0_0_x_xz_0_yy, g_z_0_0_0_x_xz_0_yz, g_z_0_0_0_x_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xz_0_xx[i] = 2.0 * g_xz_xz_0_xx[i] * a_exp;

        g_z_0_0_0_x_xz_0_xy[i] = 2.0 * g_xz_xz_0_xy[i] * a_exp;

        g_z_0_0_0_x_xz_0_xz[i] = 2.0 * g_xz_xz_0_xz[i] * a_exp;

        g_z_0_0_0_x_xz_0_yy[i] = 2.0 * g_xz_xz_0_yy[i] * a_exp;

        g_z_0_0_0_x_xz_0_yz[i] = 2.0 * g_xz_xz_0_yz[i] * a_exp;

        g_z_0_0_0_x_xz_0_zz[i] = 2.0 * g_xz_xz_0_zz[i] * a_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_xz_yy_0_xx, g_xz_yy_0_xy, g_xz_yy_0_xz, g_xz_yy_0_yy, g_xz_yy_0_yz, g_xz_yy_0_zz, g_z_0_0_0_x_yy_0_xx, g_z_0_0_0_x_yy_0_xy, g_z_0_0_0_x_yy_0_xz, g_z_0_0_0_x_yy_0_yy, g_z_0_0_0_x_yy_0_yz, g_z_0_0_0_x_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yy_0_xx[i] = 2.0 * g_xz_yy_0_xx[i] * a_exp;

        g_z_0_0_0_x_yy_0_xy[i] = 2.0 * g_xz_yy_0_xy[i] * a_exp;

        g_z_0_0_0_x_yy_0_xz[i] = 2.0 * g_xz_yy_0_xz[i] * a_exp;

        g_z_0_0_0_x_yy_0_yy[i] = 2.0 * g_xz_yy_0_yy[i] * a_exp;

        g_z_0_0_0_x_yy_0_yz[i] = 2.0 * g_xz_yy_0_yz[i] * a_exp;

        g_z_0_0_0_x_yy_0_zz[i] = 2.0 * g_xz_yy_0_zz[i] * a_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_xz_yz_0_xx, g_xz_yz_0_xy, g_xz_yz_0_xz, g_xz_yz_0_yy, g_xz_yz_0_yz, g_xz_yz_0_zz, g_z_0_0_0_x_yz_0_xx, g_z_0_0_0_x_yz_0_xy, g_z_0_0_0_x_yz_0_xz, g_z_0_0_0_x_yz_0_yy, g_z_0_0_0_x_yz_0_yz, g_z_0_0_0_x_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yz_0_xx[i] = 2.0 * g_xz_yz_0_xx[i] * a_exp;

        g_z_0_0_0_x_yz_0_xy[i] = 2.0 * g_xz_yz_0_xy[i] * a_exp;

        g_z_0_0_0_x_yz_0_xz[i] = 2.0 * g_xz_yz_0_xz[i] * a_exp;

        g_z_0_0_0_x_yz_0_yy[i] = 2.0 * g_xz_yz_0_yy[i] * a_exp;

        g_z_0_0_0_x_yz_0_yz[i] = 2.0 * g_xz_yz_0_yz[i] * a_exp;

        g_z_0_0_0_x_yz_0_zz[i] = 2.0 * g_xz_yz_0_zz[i] * a_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_xz_zz_0_xx, g_xz_zz_0_xy, g_xz_zz_0_xz, g_xz_zz_0_yy, g_xz_zz_0_yz, g_xz_zz_0_zz, g_z_0_0_0_x_zz_0_xx, g_z_0_0_0_x_zz_0_xy, g_z_0_0_0_x_zz_0_xz, g_z_0_0_0_x_zz_0_yy, g_z_0_0_0_x_zz_0_yz, g_z_0_0_0_x_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_zz_0_xx[i] = 2.0 * g_xz_zz_0_xx[i] * a_exp;

        g_z_0_0_0_x_zz_0_xy[i] = 2.0 * g_xz_zz_0_xy[i] * a_exp;

        g_z_0_0_0_x_zz_0_xz[i] = 2.0 * g_xz_zz_0_xz[i] * a_exp;

        g_z_0_0_0_x_zz_0_yy[i] = 2.0 * g_xz_zz_0_yy[i] * a_exp;

        g_z_0_0_0_x_zz_0_yz[i] = 2.0 * g_xz_zz_0_yz[i] * a_exp;

        g_z_0_0_0_x_zz_0_zz[i] = 2.0 * g_xz_zz_0_zz[i] * a_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_yz_xx_0_xx, g_yz_xx_0_xy, g_yz_xx_0_xz, g_yz_xx_0_yy, g_yz_xx_0_yz, g_yz_xx_0_zz, g_z_0_0_0_y_xx_0_xx, g_z_0_0_0_y_xx_0_xy, g_z_0_0_0_y_xx_0_xz, g_z_0_0_0_y_xx_0_yy, g_z_0_0_0_y_xx_0_yz, g_z_0_0_0_y_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xx_0_xx[i] = 2.0 * g_yz_xx_0_xx[i] * a_exp;

        g_z_0_0_0_y_xx_0_xy[i] = 2.0 * g_yz_xx_0_xy[i] * a_exp;

        g_z_0_0_0_y_xx_0_xz[i] = 2.0 * g_yz_xx_0_xz[i] * a_exp;

        g_z_0_0_0_y_xx_0_yy[i] = 2.0 * g_yz_xx_0_yy[i] * a_exp;

        g_z_0_0_0_y_xx_0_yz[i] = 2.0 * g_yz_xx_0_yz[i] * a_exp;

        g_z_0_0_0_y_xx_0_zz[i] = 2.0 * g_yz_xx_0_zz[i] * a_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_yz_xy_0_xx, g_yz_xy_0_xy, g_yz_xy_0_xz, g_yz_xy_0_yy, g_yz_xy_0_yz, g_yz_xy_0_zz, g_z_0_0_0_y_xy_0_xx, g_z_0_0_0_y_xy_0_xy, g_z_0_0_0_y_xy_0_xz, g_z_0_0_0_y_xy_0_yy, g_z_0_0_0_y_xy_0_yz, g_z_0_0_0_y_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xy_0_xx[i] = 2.0 * g_yz_xy_0_xx[i] * a_exp;

        g_z_0_0_0_y_xy_0_xy[i] = 2.0 * g_yz_xy_0_xy[i] * a_exp;

        g_z_0_0_0_y_xy_0_xz[i] = 2.0 * g_yz_xy_0_xz[i] * a_exp;

        g_z_0_0_0_y_xy_0_yy[i] = 2.0 * g_yz_xy_0_yy[i] * a_exp;

        g_z_0_0_0_y_xy_0_yz[i] = 2.0 * g_yz_xy_0_yz[i] * a_exp;

        g_z_0_0_0_y_xy_0_zz[i] = 2.0 * g_yz_xy_0_zz[i] * a_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_yz_xz_0_xx, g_yz_xz_0_xy, g_yz_xz_0_xz, g_yz_xz_0_yy, g_yz_xz_0_yz, g_yz_xz_0_zz, g_z_0_0_0_y_xz_0_xx, g_z_0_0_0_y_xz_0_xy, g_z_0_0_0_y_xz_0_xz, g_z_0_0_0_y_xz_0_yy, g_z_0_0_0_y_xz_0_yz, g_z_0_0_0_y_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xz_0_xx[i] = 2.0 * g_yz_xz_0_xx[i] * a_exp;

        g_z_0_0_0_y_xz_0_xy[i] = 2.0 * g_yz_xz_0_xy[i] * a_exp;

        g_z_0_0_0_y_xz_0_xz[i] = 2.0 * g_yz_xz_0_xz[i] * a_exp;

        g_z_0_0_0_y_xz_0_yy[i] = 2.0 * g_yz_xz_0_yy[i] * a_exp;

        g_z_0_0_0_y_xz_0_yz[i] = 2.0 * g_yz_xz_0_yz[i] * a_exp;

        g_z_0_0_0_y_xz_0_zz[i] = 2.0 * g_yz_xz_0_zz[i] * a_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_yz_yy_0_xx, g_yz_yy_0_xy, g_yz_yy_0_xz, g_yz_yy_0_yy, g_yz_yy_0_yz, g_yz_yy_0_zz, g_z_0_0_0_y_yy_0_xx, g_z_0_0_0_y_yy_0_xy, g_z_0_0_0_y_yy_0_xz, g_z_0_0_0_y_yy_0_yy, g_z_0_0_0_y_yy_0_yz, g_z_0_0_0_y_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yy_0_xx[i] = 2.0 * g_yz_yy_0_xx[i] * a_exp;

        g_z_0_0_0_y_yy_0_xy[i] = 2.0 * g_yz_yy_0_xy[i] * a_exp;

        g_z_0_0_0_y_yy_0_xz[i] = 2.0 * g_yz_yy_0_xz[i] * a_exp;

        g_z_0_0_0_y_yy_0_yy[i] = 2.0 * g_yz_yy_0_yy[i] * a_exp;

        g_z_0_0_0_y_yy_0_yz[i] = 2.0 * g_yz_yy_0_yz[i] * a_exp;

        g_z_0_0_0_y_yy_0_zz[i] = 2.0 * g_yz_yy_0_zz[i] * a_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_yz_yz_0_xx, g_yz_yz_0_xy, g_yz_yz_0_xz, g_yz_yz_0_yy, g_yz_yz_0_yz, g_yz_yz_0_zz, g_z_0_0_0_y_yz_0_xx, g_z_0_0_0_y_yz_0_xy, g_z_0_0_0_y_yz_0_xz, g_z_0_0_0_y_yz_0_yy, g_z_0_0_0_y_yz_0_yz, g_z_0_0_0_y_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yz_0_xx[i] = 2.0 * g_yz_yz_0_xx[i] * a_exp;

        g_z_0_0_0_y_yz_0_xy[i] = 2.0 * g_yz_yz_0_xy[i] * a_exp;

        g_z_0_0_0_y_yz_0_xz[i] = 2.0 * g_yz_yz_0_xz[i] * a_exp;

        g_z_0_0_0_y_yz_0_yy[i] = 2.0 * g_yz_yz_0_yy[i] * a_exp;

        g_z_0_0_0_y_yz_0_yz[i] = 2.0 * g_yz_yz_0_yz[i] * a_exp;

        g_z_0_0_0_y_yz_0_zz[i] = 2.0 * g_yz_yz_0_zz[i] * a_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_yz_zz_0_xx, g_yz_zz_0_xy, g_yz_zz_0_xz, g_yz_zz_0_yy, g_yz_zz_0_yz, g_yz_zz_0_zz, g_z_0_0_0_y_zz_0_xx, g_z_0_0_0_y_zz_0_xy, g_z_0_0_0_y_zz_0_xz, g_z_0_0_0_y_zz_0_yy, g_z_0_0_0_y_zz_0_yz, g_z_0_0_0_y_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_zz_0_xx[i] = 2.0 * g_yz_zz_0_xx[i] * a_exp;

        g_z_0_0_0_y_zz_0_xy[i] = 2.0 * g_yz_zz_0_xy[i] * a_exp;

        g_z_0_0_0_y_zz_0_xz[i] = 2.0 * g_yz_zz_0_xz[i] * a_exp;

        g_z_0_0_0_y_zz_0_yy[i] = 2.0 * g_yz_zz_0_yy[i] * a_exp;

        g_z_0_0_0_y_zz_0_yz[i] = 2.0 * g_yz_zz_0_yz[i] * a_exp;

        g_z_0_0_0_y_zz_0_zz[i] = 2.0 * g_yz_zz_0_zz[i] * a_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_0_xx_0_xx, g_0_xx_0_xy, g_0_xx_0_xz, g_0_xx_0_yy, g_0_xx_0_yz, g_0_xx_0_zz, g_z_0_0_0_z_xx_0_xx, g_z_0_0_0_z_xx_0_xy, g_z_0_0_0_z_xx_0_xz, g_z_0_0_0_z_xx_0_yy, g_z_0_0_0_z_xx_0_yz, g_z_0_0_0_z_xx_0_zz, g_zz_xx_0_xx, g_zz_xx_0_xy, g_zz_xx_0_xz, g_zz_xx_0_yy, g_zz_xx_0_yz, g_zz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xx_0_xx[i] = -g_0_xx_0_xx[i] + 2.0 * g_zz_xx_0_xx[i] * a_exp;

        g_z_0_0_0_z_xx_0_xy[i] = -g_0_xx_0_xy[i] + 2.0 * g_zz_xx_0_xy[i] * a_exp;

        g_z_0_0_0_z_xx_0_xz[i] = -g_0_xx_0_xz[i] + 2.0 * g_zz_xx_0_xz[i] * a_exp;

        g_z_0_0_0_z_xx_0_yy[i] = -g_0_xx_0_yy[i] + 2.0 * g_zz_xx_0_yy[i] * a_exp;

        g_z_0_0_0_z_xx_0_yz[i] = -g_0_xx_0_yz[i] + 2.0 * g_zz_xx_0_yz[i] * a_exp;

        g_z_0_0_0_z_xx_0_zz[i] = -g_0_xx_0_zz[i] + 2.0 * g_zz_xx_0_zz[i] * a_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_z_0_0_0_z_xy_0_xx, g_z_0_0_0_z_xy_0_xy, g_z_0_0_0_z_xy_0_xz, g_z_0_0_0_z_xy_0_yy, g_z_0_0_0_z_xy_0_yz, g_z_0_0_0_z_xy_0_zz, g_zz_xy_0_xx, g_zz_xy_0_xy, g_zz_xy_0_xz, g_zz_xy_0_yy, g_zz_xy_0_yz, g_zz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xy_0_xx[i] = -g_0_xy_0_xx[i] + 2.0 * g_zz_xy_0_xx[i] * a_exp;

        g_z_0_0_0_z_xy_0_xy[i] = -g_0_xy_0_xy[i] + 2.0 * g_zz_xy_0_xy[i] * a_exp;

        g_z_0_0_0_z_xy_0_xz[i] = -g_0_xy_0_xz[i] + 2.0 * g_zz_xy_0_xz[i] * a_exp;

        g_z_0_0_0_z_xy_0_yy[i] = -g_0_xy_0_yy[i] + 2.0 * g_zz_xy_0_yy[i] * a_exp;

        g_z_0_0_0_z_xy_0_yz[i] = -g_0_xy_0_yz[i] + 2.0 * g_zz_xy_0_yz[i] * a_exp;

        g_z_0_0_0_z_xy_0_zz[i] = -g_0_xy_0_zz[i] + 2.0 * g_zz_xy_0_zz[i] * a_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_z_0_0_0_z_xz_0_xx, g_z_0_0_0_z_xz_0_xy, g_z_0_0_0_z_xz_0_xz, g_z_0_0_0_z_xz_0_yy, g_z_0_0_0_z_xz_0_yz, g_z_0_0_0_z_xz_0_zz, g_zz_xz_0_xx, g_zz_xz_0_xy, g_zz_xz_0_xz, g_zz_xz_0_yy, g_zz_xz_0_yz, g_zz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xz_0_xx[i] = -g_0_xz_0_xx[i] + 2.0 * g_zz_xz_0_xx[i] * a_exp;

        g_z_0_0_0_z_xz_0_xy[i] = -g_0_xz_0_xy[i] + 2.0 * g_zz_xz_0_xy[i] * a_exp;

        g_z_0_0_0_z_xz_0_xz[i] = -g_0_xz_0_xz[i] + 2.0 * g_zz_xz_0_xz[i] * a_exp;

        g_z_0_0_0_z_xz_0_yy[i] = -g_0_xz_0_yy[i] + 2.0 * g_zz_xz_0_yy[i] * a_exp;

        g_z_0_0_0_z_xz_0_yz[i] = -g_0_xz_0_yz[i] + 2.0 * g_zz_xz_0_yz[i] * a_exp;

        g_z_0_0_0_z_xz_0_zz[i] = -g_0_xz_0_zz[i] + 2.0 * g_zz_xz_0_zz[i] * a_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_0_yy_0_xx, g_0_yy_0_xy, g_0_yy_0_xz, g_0_yy_0_yy, g_0_yy_0_yz, g_0_yy_0_zz, g_z_0_0_0_z_yy_0_xx, g_z_0_0_0_z_yy_0_xy, g_z_0_0_0_z_yy_0_xz, g_z_0_0_0_z_yy_0_yy, g_z_0_0_0_z_yy_0_yz, g_z_0_0_0_z_yy_0_zz, g_zz_yy_0_xx, g_zz_yy_0_xy, g_zz_yy_0_xz, g_zz_yy_0_yy, g_zz_yy_0_yz, g_zz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yy_0_xx[i] = -g_0_yy_0_xx[i] + 2.0 * g_zz_yy_0_xx[i] * a_exp;

        g_z_0_0_0_z_yy_0_xy[i] = -g_0_yy_0_xy[i] + 2.0 * g_zz_yy_0_xy[i] * a_exp;

        g_z_0_0_0_z_yy_0_xz[i] = -g_0_yy_0_xz[i] + 2.0 * g_zz_yy_0_xz[i] * a_exp;

        g_z_0_0_0_z_yy_0_yy[i] = -g_0_yy_0_yy[i] + 2.0 * g_zz_yy_0_yy[i] * a_exp;

        g_z_0_0_0_z_yy_0_yz[i] = -g_0_yy_0_yz[i] + 2.0 * g_zz_yy_0_yz[i] * a_exp;

        g_z_0_0_0_z_yy_0_zz[i] = -g_0_yy_0_zz[i] + 2.0 * g_zz_yy_0_zz[i] * a_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_z_0_0_0_z_yz_0_xx, g_z_0_0_0_z_yz_0_xy, g_z_0_0_0_z_yz_0_xz, g_z_0_0_0_z_yz_0_yy, g_z_0_0_0_z_yz_0_yz, g_z_0_0_0_z_yz_0_zz, g_zz_yz_0_xx, g_zz_yz_0_xy, g_zz_yz_0_xz, g_zz_yz_0_yy, g_zz_yz_0_yz, g_zz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yz_0_xx[i] = -g_0_yz_0_xx[i] + 2.0 * g_zz_yz_0_xx[i] * a_exp;

        g_z_0_0_0_z_yz_0_xy[i] = -g_0_yz_0_xy[i] + 2.0 * g_zz_yz_0_xy[i] * a_exp;

        g_z_0_0_0_z_yz_0_xz[i] = -g_0_yz_0_xz[i] + 2.0 * g_zz_yz_0_xz[i] * a_exp;

        g_z_0_0_0_z_yz_0_yy[i] = -g_0_yz_0_yy[i] + 2.0 * g_zz_yz_0_yy[i] * a_exp;

        g_z_0_0_0_z_yz_0_yz[i] = -g_0_yz_0_yz[i] + 2.0 * g_zz_yz_0_yz[i] * a_exp;

        g_z_0_0_0_z_yz_0_zz[i] = -g_0_yz_0_zz[i] + 2.0 * g_zz_yz_0_zz[i] * a_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_0_zz_0_xx, g_0_zz_0_xy, g_0_zz_0_xz, g_0_zz_0_yy, g_0_zz_0_yz, g_0_zz_0_zz, g_z_0_0_0_z_zz_0_xx, g_z_0_0_0_z_zz_0_xy, g_z_0_0_0_z_zz_0_xz, g_z_0_0_0_z_zz_0_yy, g_z_0_0_0_z_zz_0_yz, g_z_0_0_0_z_zz_0_zz, g_zz_zz_0_xx, g_zz_zz_0_xy, g_zz_zz_0_xz, g_zz_zz_0_yy, g_zz_zz_0_yz, g_zz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_zz_0_xx[i] = -g_0_zz_0_xx[i] + 2.0 * g_zz_zz_0_xx[i] * a_exp;

        g_z_0_0_0_z_zz_0_xy[i] = -g_0_zz_0_xy[i] + 2.0 * g_zz_zz_0_xy[i] * a_exp;

        g_z_0_0_0_z_zz_0_xz[i] = -g_0_zz_0_xz[i] + 2.0 * g_zz_zz_0_xz[i] * a_exp;

        g_z_0_0_0_z_zz_0_yy[i] = -g_0_zz_0_yy[i] + 2.0 * g_zz_zz_0_yy[i] * a_exp;

        g_z_0_0_0_z_zz_0_yz[i] = -g_0_zz_0_yz[i] + 2.0 * g_zz_zz_0_yz[i] * a_exp;

        g_z_0_0_0_z_zz_0_zz[i] = -g_0_zz_0_zz[i] + 2.0 * g_zz_zz_0_zz[i] * a_exp;
    }
}

} // t4c_geom namespace

