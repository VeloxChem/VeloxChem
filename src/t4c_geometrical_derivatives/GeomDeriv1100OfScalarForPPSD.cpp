#include "GeomDeriv1100OfScalarForPPSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_ppsd_0(CSimdArray<double>& buffer_1100_ppsd,
                     const CSimdArray<double>& buffer_sssd,
                     const CSimdArray<double>& buffer_sdsd,
                     const CSimdArray<double>& buffer_dssd,
                     const CSimdArray<double>& buffer_ddsd,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_ppsd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sssd

    auto g_0_0_0_xx = buffer_sssd[0];

    auto g_0_0_0_xy = buffer_sssd[1];

    auto g_0_0_0_xz = buffer_sssd[2];

    auto g_0_0_0_yy = buffer_sssd[3];

    auto g_0_0_0_yz = buffer_sssd[4];

    auto g_0_0_0_zz = buffer_sssd[5];

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

    /// Set up components of auxilary buffer : buffer_dssd

    auto g_xx_0_0_xx = buffer_dssd[0];

    auto g_xx_0_0_xy = buffer_dssd[1];

    auto g_xx_0_0_xz = buffer_dssd[2];

    auto g_xx_0_0_yy = buffer_dssd[3];

    auto g_xx_0_0_yz = buffer_dssd[4];

    auto g_xx_0_0_zz = buffer_dssd[5];

    auto g_xy_0_0_xx = buffer_dssd[6];

    auto g_xy_0_0_xy = buffer_dssd[7];

    auto g_xy_0_0_xz = buffer_dssd[8];

    auto g_xy_0_0_yy = buffer_dssd[9];

    auto g_xy_0_0_yz = buffer_dssd[10];

    auto g_xy_0_0_zz = buffer_dssd[11];

    auto g_xz_0_0_xx = buffer_dssd[12];

    auto g_xz_0_0_xy = buffer_dssd[13];

    auto g_xz_0_0_xz = buffer_dssd[14];

    auto g_xz_0_0_yy = buffer_dssd[15];

    auto g_xz_0_0_yz = buffer_dssd[16];

    auto g_xz_0_0_zz = buffer_dssd[17];

    auto g_yy_0_0_xx = buffer_dssd[18];

    auto g_yy_0_0_xy = buffer_dssd[19];

    auto g_yy_0_0_xz = buffer_dssd[20];

    auto g_yy_0_0_yy = buffer_dssd[21];

    auto g_yy_0_0_yz = buffer_dssd[22];

    auto g_yy_0_0_zz = buffer_dssd[23];

    auto g_yz_0_0_xx = buffer_dssd[24];

    auto g_yz_0_0_xy = buffer_dssd[25];

    auto g_yz_0_0_xz = buffer_dssd[26];

    auto g_yz_0_0_yy = buffer_dssd[27];

    auto g_yz_0_0_yz = buffer_dssd[28];

    auto g_yz_0_0_zz = buffer_dssd[29];

    auto g_zz_0_0_xx = buffer_dssd[30];

    auto g_zz_0_0_xy = buffer_dssd[31];

    auto g_zz_0_0_xz = buffer_dssd[32];

    auto g_zz_0_0_yy = buffer_dssd[33];

    auto g_zz_0_0_yz = buffer_dssd[34];

    auto g_zz_0_0_zz = buffer_dssd[35];

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

    /// Set up components of integrals buffer : buffer_1100_ppsd

    auto g_x_x_0_0_x_x_0_xx = buffer_1100_ppsd[0];

    auto g_x_x_0_0_x_x_0_xy = buffer_1100_ppsd[1];

    auto g_x_x_0_0_x_x_0_xz = buffer_1100_ppsd[2];

    auto g_x_x_0_0_x_x_0_yy = buffer_1100_ppsd[3];

    auto g_x_x_0_0_x_x_0_yz = buffer_1100_ppsd[4];

    auto g_x_x_0_0_x_x_0_zz = buffer_1100_ppsd[5];

    auto g_x_x_0_0_x_y_0_xx = buffer_1100_ppsd[6];

    auto g_x_x_0_0_x_y_0_xy = buffer_1100_ppsd[7];

    auto g_x_x_0_0_x_y_0_xz = buffer_1100_ppsd[8];

    auto g_x_x_0_0_x_y_0_yy = buffer_1100_ppsd[9];

    auto g_x_x_0_0_x_y_0_yz = buffer_1100_ppsd[10];

    auto g_x_x_0_0_x_y_0_zz = buffer_1100_ppsd[11];

    auto g_x_x_0_0_x_z_0_xx = buffer_1100_ppsd[12];

    auto g_x_x_0_0_x_z_0_xy = buffer_1100_ppsd[13];

    auto g_x_x_0_0_x_z_0_xz = buffer_1100_ppsd[14];

    auto g_x_x_0_0_x_z_0_yy = buffer_1100_ppsd[15];

    auto g_x_x_0_0_x_z_0_yz = buffer_1100_ppsd[16];

    auto g_x_x_0_0_x_z_0_zz = buffer_1100_ppsd[17];

    auto g_x_x_0_0_y_x_0_xx = buffer_1100_ppsd[18];

    auto g_x_x_0_0_y_x_0_xy = buffer_1100_ppsd[19];

    auto g_x_x_0_0_y_x_0_xz = buffer_1100_ppsd[20];

    auto g_x_x_0_0_y_x_0_yy = buffer_1100_ppsd[21];

    auto g_x_x_0_0_y_x_0_yz = buffer_1100_ppsd[22];

    auto g_x_x_0_0_y_x_0_zz = buffer_1100_ppsd[23];

    auto g_x_x_0_0_y_y_0_xx = buffer_1100_ppsd[24];

    auto g_x_x_0_0_y_y_0_xy = buffer_1100_ppsd[25];

    auto g_x_x_0_0_y_y_0_xz = buffer_1100_ppsd[26];

    auto g_x_x_0_0_y_y_0_yy = buffer_1100_ppsd[27];

    auto g_x_x_0_0_y_y_0_yz = buffer_1100_ppsd[28];

    auto g_x_x_0_0_y_y_0_zz = buffer_1100_ppsd[29];

    auto g_x_x_0_0_y_z_0_xx = buffer_1100_ppsd[30];

    auto g_x_x_0_0_y_z_0_xy = buffer_1100_ppsd[31];

    auto g_x_x_0_0_y_z_0_xz = buffer_1100_ppsd[32];

    auto g_x_x_0_0_y_z_0_yy = buffer_1100_ppsd[33];

    auto g_x_x_0_0_y_z_0_yz = buffer_1100_ppsd[34];

    auto g_x_x_0_0_y_z_0_zz = buffer_1100_ppsd[35];

    auto g_x_x_0_0_z_x_0_xx = buffer_1100_ppsd[36];

    auto g_x_x_0_0_z_x_0_xy = buffer_1100_ppsd[37];

    auto g_x_x_0_0_z_x_0_xz = buffer_1100_ppsd[38];

    auto g_x_x_0_0_z_x_0_yy = buffer_1100_ppsd[39];

    auto g_x_x_0_0_z_x_0_yz = buffer_1100_ppsd[40];

    auto g_x_x_0_0_z_x_0_zz = buffer_1100_ppsd[41];

    auto g_x_x_0_0_z_y_0_xx = buffer_1100_ppsd[42];

    auto g_x_x_0_0_z_y_0_xy = buffer_1100_ppsd[43];

    auto g_x_x_0_0_z_y_0_xz = buffer_1100_ppsd[44];

    auto g_x_x_0_0_z_y_0_yy = buffer_1100_ppsd[45];

    auto g_x_x_0_0_z_y_0_yz = buffer_1100_ppsd[46];

    auto g_x_x_0_0_z_y_0_zz = buffer_1100_ppsd[47];

    auto g_x_x_0_0_z_z_0_xx = buffer_1100_ppsd[48];

    auto g_x_x_0_0_z_z_0_xy = buffer_1100_ppsd[49];

    auto g_x_x_0_0_z_z_0_xz = buffer_1100_ppsd[50];

    auto g_x_x_0_0_z_z_0_yy = buffer_1100_ppsd[51];

    auto g_x_x_0_0_z_z_0_yz = buffer_1100_ppsd[52];

    auto g_x_x_0_0_z_z_0_zz = buffer_1100_ppsd[53];

    auto g_x_y_0_0_x_x_0_xx = buffer_1100_ppsd[54];

    auto g_x_y_0_0_x_x_0_xy = buffer_1100_ppsd[55];

    auto g_x_y_0_0_x_x_0_xz = buffer_1100_ppsd[56];

    auto g_x_y_0_0_x_x_0_yy = buffer_1100_ppsd[57];

    auto g_x_y_0_0_x_x_0_yz = buffer_1100_ppsd[58];

    auto g_x_y_0_0_x_x_0_zz = buffer_1100_ppsd[59];

    auto g_x_y_0_0_x_y_0_xx = buffer_1100_ppsd[60];

    auto g_x_y_0_0_x_y_0_xy = buffer_1100_ppsd[61];

    auto g_x_y_0_0_x_y_0_xz = buffer_1100_ppsd[62];

    auto g_x_y_0_0_x_y_0_yy = buffer_1100_ppsd[63];

    auto g_x_y_0_0_x_y_0_yz = buffer_1100_ppsd[64];

    auto g_x_y_0_0_x_y_0_zz = buffer_1100_ppsd[65];

    auto g_x_y_0_0_x_z_0_xx = buffer_1100_ppsd[66];

    auto g_x_y_0_0_x_z_0_xy = buffer_1100_ppsd[67];

    auto g_x_y_0_0_x_z_0_xz = buffer_1100_ppsd[68];

    auto g_x_y_0_0_x_z_0_yy = buffer_1100_ppsd[69];

    auto g_x_y_0_0_x_z_0_yz = buffer_1100_ppsd[70];

    auto g_x_y_0_0_x_z_0_zz = buffer_1100_ppsd[71];

    auto g_x_y_0_0_y_x_0_xx = buffer_1100_ppsd[72];

    auto g_x_y_0_0_y_x_0_xy = buffer_1100_ppsd[73];

    auto g_x_y_0_0_y_x_0_xz = buffer_1100_ppsd[74];

    auto g_x_y_0_0_y_x_0_yy = buffer_1100_ppsd[75];

    auto g_x_y_0_0_y_x_0_yz = buffer_1100_ppsd[76];

    auto g_x_y_0_0_y_x_0_zz = buffer_1100_ppsd[77];

    auto g_x_y_0_0_y_y_0_xx = buffer_1100_ppsd[78];

    auto g_x_y_0_0_y_y_0_xy = buffer_1100_ppsd[79];

    auto g_x_y_0_0_y_y_0_xz = buffer_1100_ppsd[80];

    auto g_x_y_0_0_y_y_0_yy = buffer_1100_ppsd[81];

    auto g_x_y_0_0_y_y_0_yz = buffer_1100_ppsd[82];

    auto g_x_y_0_0_y_y_0_zz = buffer_1100_ppsd[83];

    auto g_x_y_0_0_y_z_0_xx = buffer_1100_ppsd[84];

    auto g_x_y_0_0_y_z_0_xy = buffer_1100_ppsd[85];

    auto g_x_y_0_0_y_z_0_xz = buffer_1100_ppsd[86];

    auto g_x_y_0_0_y_z_0_yy = buffer_1100_ppsd[87];

    auto g_x_y_0_0_y_z_0_yz = buffer_1100_ppsd[88];

    auto g_x_y_0_0_y_z_0_zz = buffer_1100_ppsd[89];

    auto g_x_y_0_0_z_x_0_xx = buffer_1100_ppsd[90];

    auto g_x_y_0_0_z_x_0_xy = buffer_1100_ppsd[91];

    auto g_x_y_0_0_z_x_0_xz = buffer_1100_ppsd[92];

    auto g_x_y_0_0_z_x_0_yy = buffer_1100_ppsd[93];

    auto g_x_y_0_0_z_x_0_yz = buffer_1100_ppsd[94];

    auto g_x_y_0_0_z_x_0_zz = buffer_1100_ppsd[95];

    auto g_x_y_0_0_z_y_0_xx = buffer_1100_ppsd[96];

    auto g_x_y_0_0_z_y_0_xy = buffer_1100_ppsd[97];

    auto g_x_y_0_0_z_y_0_xz = buffer_1100_ppsd[98];

    auto g_x_y_0_0_z_y_0_yy = buffer_1100_ppsd[99];

    auto g_x_y_0_0_z_y_0_yz = buffer_1100_ppsd[100];

    auto g_x_y_0_0_z_y_0_zz = buffer_1100_ppsd[101];

    auto g_x_y_0_0_z_z_0_xx = buffer_1100_ppsd[102];

    auto g_x_y_0_0_z_z_0_xy = buffer_1100_ppsd[103];

    auto g_x_y_0_0_z_z_0_xz = buffer_1100_ppsd[104];

    auto g_x_y_0_0_z_z_0_yy = buffer_1100_ppsd[105];

    auto g_x_y_0_0_z_z_0_yz = buffer_1100_ppsd[106];

    auto g_x_y_0_0_z_z_0_zz = buffer_1100_ppsd[107];

    auto g_x_z_0_0_x_x_0_xx = buffer_1100_ppsd[108];

    auto g_x_z_0_0_x_x_0_xy = buffer_1100_ppsd[109];

    auto g_x_z_0_0_x_x_0_xz = buffer_1100_ppsd[110];

    auto g_x_z_0_0_x_x_0_yy = buffer_1100_ppsd[111];

    auto g_x_z_0_0_x_x_0_yz = buffer_1100_ppsd[112];

    auto g_x_z_0_0_x_x_0_zz = buffer_1100_ppsd[113];

    auto g_x_z_0_0_x_y_0_xx = buffer_1100_ppsd[114];

    auto g_x_z_0_0_x_y_0_xy = buffer_1100_ppsd[115];

    auto g_x_z_0_0_x_y_0_xz = buffer_1100_ppsd[116];

    auto g_x_z_0_0_x_y_0_yy = buffer_1100_ppsd[117];

    auto g_x_z_0_0_x_y_0_yz = buffer_1100_ppsd[118];

    auto g_x_z_0_0_x_y_0_zz = buffer_1100_ppsd[119];

    auto g_x_z_0_0_x_z_0_xx = buffer_1100_ppsd[120];

    auto g_x_z_0_0_x_z_0_xy = buffer_1100_ppsd[121];

    auto g_x_z_0_0_x_z_0_xz = buffer_1100_ppsd[122];

    auto g_x_z_0_0_x_z_0_yy = buffer_1100_ppsd[123];

    auto g_x_z_0_0_x_z_0_yz = buffer_1100_ppsd[124];

    auto g_x_z_0_0_x_z_0_zz = buffer_1100_ppsd[125];

    auto g_x_z_0_0_y_x_0_xx = buffer_1100_ppsd[126];

    auto g_x_z_0_0_y_x_0_xy = buffer_1100_ppsd[127];

    auto g_x_z_0_0_y_x_0_xz = buffer_1100_ppsd[128];

    auto g_x_z_0_0_y_x_0_yy = buffer_1100_ppsd[129];

    auto g_x_z_0_0_y_x_0_yz = buffer_1100_ppsd[130];

    auto g_x_z_0_0_y_x_0_zz = buffer_1100_ppsd[131];

    auto g_x_z_0_0_y_y_0_xx = buffer_1100_ppsd[132];

    auto g_x_z_0_0_y_y_0_xy = buffer_1100_ppsd[133];

    auto g_x_z_0_0_y_y_0_xz = buffer_1100_ppsd[134];

    auto g_x_z_0_0_y_y_0_yy = buffer_1100_ppsd[135];

    auto g_x_z_0_0_y_y_0_yz = buffer_1100_ppsd[136];

    auto g_x_z_0_0_y_y_0_zz = buffer_1100_ppsd[137];

    auto g_x_z_0_0_y_z_0_xx = buffer_1100_ppsd[138];

    auto g_x_z_0_0_y_z_0_xy = buffer_1100_ppsd[139];

    auto g_x_z_0_0_y_z_0_xz = buffer_1100_ppsd[140];

    auto g_x_z_0_0_y_z_0_yy = buffer_1100_ppsd[141];

    auto g_x_z_0_0_y_z_0_yz = buffer_1100_ppsd[142];

    auto g_x_z_0_0_y_z_0_zz = buffer_1100_ppsd[143];

    auto g_x_z_0_0_z_x_0_xx = buffer_1100_ppsd[144];

    auto g_x_z_0_0_z_x_0_xy = buffer_1100_ppsd[145];

    auto g_x_z_0_0_z_x_0_xz = buffer_1100_ppsd[146];

    auto g_x_z_0_0_z_x_0_yy = buffer_1100_ppsd[147];

    auto g_x_z_0_0_z_x_0_yz = buffer_1100_ppsd[148];

    auto g_x_z_0_0_z_x_0_zz = buffer_1100_ppsd[149];

    auto g_x_z_0_0_z_y_0_xx = buffer_1100_ppsd[150];

    auto g_x_z_0_0_z_y_0_xy = buffer_1100_ppsd[151];

    auto g_x_z_0_0_z_y_0_xz = buffer_1100_ppsd[152];

    auto g_x_z_0_0_z_y_0_yy = buffer_1100_ppsd[153];

    auto g_x_z_0_0_z_y_0_yz = buffer_1100_ppsd[154];

    auto g_x_z_0_0_z_y_0_zz = buffer_1100_ppsd[155];

    auto g_x_z_0_0_z_z_0_xx = buffer_1100_ppsd[156];

    auto g_x_z_0_0_z_z_0_xy = buffer_1100_ppsd[157];

    auto g_x_z_0_0_z_z_0_xz = buffer_1100_ppsd[158];

    auto g_x_z_0_0_z_z_0_yy = buffer_1100_ppsd[159];

    auto g_x_z_0_0_z_z_0_yz = buffer_1100_ppsd[160];

    auto g_x_z_0_0_z_z_0_zz = buffer_1100_ppsd[161];

    auto g_y_x_0_0_x_x_0_xx = buffer_1100_ppsd[162];

    auto g_y_x_0_0_x_x_0_xy = buffer_1100_ppsd[163];

    auto g_y_x_0_0_x_x_0_xz = buffer_1100_ppsd[164];

    auto g_y_x_0_0_x_x_0_yy = buffer_1100_ppsd[165];

    auto g_y_x_0_0_x_x_0_yz = buffer_1100_ppsd[166];

    auto g_y_x_0_0_x_x_0_zz = buffer_1100_ppsd[167];

    auto g_y_x_0_0_x_y_0_xx = buffer_1100_ppsd[168];

    auto g_y_x_0_0_x_y_0_xy = buffer_1100_ppsd[169];

    auto g_y_x_0_0_x_y_0_xz = buffer_1100_ppsd[170];

    auto g_y_x_0_0_x_y_0_yy = buffer_1100_ppsd[171];

    auto g_y_x_0_0_x_y_0_yz = buffer_1100_ppsd[172];

    auto g_y_x_0_0_x_y_0_zz = buffer_1100_ppsd[173];

    auto g_y_x_0_0_x_z_0_xx = buffer_1100_ppsd[174];

    auto g_y_x_0_0_x_z_0_xy = buffer_1100_ppsd[175];

    auto g_y_x_0_0_x_z_0_xz = buffer_1100_ppsd[176];

    auto g_y_x_0_0_x_z_0_yy = buffer_1100_ppsd[177];

    auto g_y_x_0_0_x_z_0_yz = buffer_1100_ppsd[178];

    auto g_y_x_0_0_x_z_0_zz = buffer_1100_ppsd[179];

    auto g_y_x_0_0_y_x_0_xx = buffer_1100_ppsd[180];

    auto g_y_x_0_0_y_x_0_xy = buffer_1100_ppsd[181];

    auto g_y_x_0_0_y_x_0_xz = buffer_1100_ppsd[182];

    auto g_y_x_0_0_y_x_0_yy = buffer_1100_ppsd[183];

    auto g_y_x_0_0_y_x_0_yz = buffer_1100_ppsd[184];

    auto g_y_x_0_0_y_x_0_zz = buffer_1100_ppsd[185];

    auto g_y_x_0_0_y_y_0_xx = buffer_1100_ppsd[186];

    auto g_y_x_0_0_y_y_0_xy = buffer_1100_ppsd[187];

    auto g_y_x_0_0_y_y_0_xz = buffer_1100_ppsd[188];

    auto g_y_x_0_0_y_y_0_yy = buffer_1100_ppsd[189];

    auto g_y_x_0_0_y_y_0_yz = buffer_1100_ppsd[190];

    auto g_y_x_0_0_y_y_0_zz = buffer_1100_ppsd[191];

    auto g_y_x_0_0_y_z_0_xx = buffer_1100_ppsd[192];

    auto g_y_x_0_0_y_z_0_xy = buffer_1100_ppsd[193];

    auto g_y_x_0_0_y_z_0_xz = buffer_1100_ppsd[194];

    auto g_y_x_0_0_y_z_0_yy = buffer_1100_ppsd[195];

    auto g_y_x_0_0_y_z_0_yz = buffer_1100_ppsd[196];

    auto g_y_x_0_0_y_z_0_zz = buffer_1100_ppsd[197];

    auto g_y_x_0_0_z_x_0_xx = buffer_1100_ppsd[198];

    auto g_y_x_0_0_z_x_0_xy = buffer_1100_ppsd[199];

    auto g_y_x_0_0_z_x_0_xz = buffer_1100_ppsd[200];

    auto g_y_x_0_0_z_x_0_yy = buffer_1100_ppsd[201];

    auto g_y_x_0_0_z_x_0_yz = buffer_1100_ppsd[202];

    auto g_y_x_0_0_z_x_0_zz = buffer_1100_ppsd[203];

    auto g_y_x_0_0_z_y_0_xx = buffer_1100_ppsd[204];

    auto g_y_x_0_0_z_y_0_xy = buffer_1100_ppsd[205];

    auto g_y_x_0_0_z_y_0_xz = buffer_1100_ppsd[206];

    auto g_y_x_0_0_z_y_0_yy = buffer_1100_ppsd[207];

    auto g_y_x_0_0_z_y_0_yz = buffer_1100_ppsd[208];

    auto g_y_x_0_0_z_y_0_zz = buffer_1100_ppsd[209];

    auto g_y_x_0_0_z_z_0_xx = buffer_1100_ppsd[210];

    auto g_y_x_0_0_z_z_0_xy = buffer_1100_ppsd[211];

    auto g_y_x_0_0_z_z_0_xz = buffer_1100_ppsd[212];

    auto g_y_x_0_0_z_z_0_yy = buffer_1100_ppsd[213];

    auto g_y_x_0_0_z_z_0_yz = buffer_1100_ppsd[214];

    auto g_y_x_0_0_z_z_0_zz = buffer_1100_ppsd[215];

    auto g_y_y_0_0_x_x_0_xx = buffer_1100_ppsd[216];

    auto g_y_y_0_0_x_x_0_xy = buffer_1100_ppsd[217];

    auto g_y_y_0_0_x_x_0_xz = buffer_1100_ppsd[218];

    auto g_y_y_0_0_x_x_0_yy = buffer_1100_ppsd[219];

    auto g_y_y_0_0_x_x_0_yz = buffer_1100_ppsd[220];

    auto g_y_y_0_0_x_x_0_zz = buffer_1100_ppsd[221];

    auto g_y_y_0_0_x_y_0_xx = buffer_1100_ppsd[222];

    auto g_y_y_0_0_x_y_0_xy = buffer_1100_ppsd[223];

    auto g_y_y_0_0_x_y_0_xz = buffer_1100_ppsd[224];

    auto g_y_y_0_0_x_y_0_yy = buffer_1100_ppsd[225];

    auto g_y_y_0_0_x_y_0_yz = buffer_1100_ppsd[226];

    auto g_y_y_0_0_x_y_0_zz = buffer_1100_ppsd[227];

    auto g_y_y_0_0_x_z_0_xx = buffer_1100_ppsd[228];

    auto g_y_y_0_0_x_z_0_xy = buffer_1100_ppsd[229];

    auto g_y_y_0_0_x_z_0_xz = buffer_1100_ppsd[230];

    auto g_y_y_0_0_x_z_0_yy = buffer_1100_ppsd[231];

    auto g_y_y_0_0_x_z_0_yz = buffer_1100_ppsd[232];

    auto g_y_y_0_0_x_z_0_zz = buffer_1100_ppsd[233];

    auto g_y_y_0_0_y_x_0_xx = buffer_1100_ppsd[234];

    auto g_y_y_0_0_y_x_0_xy = buffer_1100_ppsd[235];

    auto g_y_y_0_0_y_x_0_xz = buffer_1100_ppsd[236];

    auto g_y_y_0_0_y_x_0_yy = buffer_1100_ppsd[237];

    auto g_y_y_0_0_y_x_0_yz = buffer_1100_ppsd[238];

    auto g_y_y_0_0_y_x_0_zz = buffer_1100_ppsd[239];

    auto g_y_y_0_0_y_y_0_xx = buffer_1100_ppsd[240];

    auto g_y_y_0_0_y_y_0_xy = buffer_1100_ppsd[241];

    auto g_y_y_0_0_y_y_0_xz = buffer_1100_ppsd[242];

    auto g_y_y_0_0_y_y_0_yy = buffer_1100_ppsd[243];

    auto g_y_y_0_0_y_y_0_yz = buffer_1100_ppsd[244];

    auto g_y_y_0_0_y_y_0_zz = buffer_1100_ppsd[245];

    auto g_y_y_0_0_y_z_0_xx = buffer_1100_ppsd[246];

    auto g_y_y_0_0_y_z_0_xy = buffer_1100_ppsd[247];

    auto g_y_y_0_0_y_z_0_xz = buffer_1100_ppsd[248];

    auto g_y_y_0_0_y_z_0_yy = buffer_1100_ppsd[249];

    auto g_y_y_0_0_y_z_0_yz = buffer_1100_ppsd[250];

    auto g_y_y_0_0_y_z_0_zz = buffer_1100_ppsd[251];

    auto g_y_y_0_0_z_x_0_xx = buffer_1100_ppsd[252];

    auto g_y_y_0_0_z_x_0_xy = buffer_1100_ppsd[253];

    auto g_y_y_0_0_z_x_0_xz = buffer_1100_ppsd[254];

    auto g_y_y_0_0_z_x_0_yy = buffer_1100_ppsd[255];

    auto g_y_y_0_0_z_x_0_yz = buffer_1100_ppsd[256];

    auto g_y_y_0_0_z_x_0_zz = buffer_1100_ppsd[257];

    auto g_y_y_0_0_z_y_0_xx = buffer_1100_ppsd[258];

    auto g_y_y_0_0_z_y_0_xy = buffer_1100_ppsd[259];

    auto g_y_y_0_0_z_y_0_xz = buffer_1100_ppsd[260];

    auto g_y_y_0_0_z_y_0_yy = buffer_1100_ppsd[261];

    auto g_y_y_0_0_z_y_0_yz = buffer_1100_ppsd[262];

    auto g_y_y_0_0_z_y_0_zz = buffer_1100_ppsd[263];

    auto g_y_y_0_0_z_z_0_xx = buffer_1100_ppsd[264];

    auto g_y_y_0_0_z_z_0_xy = buffer_1100_ppsd[265];

    auto g_y_y_0_0_z_z_0_xz = buffer_1100_ppsd[266];

    auto g_y_y_0_0_z_z_0_yy = buffer_1100_ppsd[267];

    auto g_y_y_0_0_z_z_0_yz = buffer_1100_ppsd[268];

    auto g_y_y_0_0_z_z_0_zz = buffer_1100_ppsd[269];

    auto g_y_z_0_0_x_x_0_xx = buffer_1100_ppsd[270];

    auto g_y_z_0_0_x_x_0_xy = buffer_1100_ppsd[271];

    auto g_y_z_0_0_x_x_0_xz = buffer_1100_ppsd[272];

    auto g_y_z_0_0_x_x_0_yy = buffer_1100_ppsd[273];

    auto g_y_z_0_0_x_x_0_yz = buffer_1100_ppsd[274];

    auto g_y_z_0_0_x_x_0_zz = buffer_1100_ppsd[275];

    auto g_y_z_0_0_x_y_0_xx = buffer_1100_ppsd[276];

    auto g_y_z_0_0_x_y_0_xy = buffer_1100_ppsd[277];

    auto g_y_z_0_0_x_y_0_xz = buffer_1100_ppsd[278];

    auto g_y_z_0_0_x_y_0_yy = buffer_1100_ppsd[279];

    auto g_y_z_0_0_x_y_0_yz = buffer_1100_ppsd[280];

    auto g_y_z_0_0_x_y_0_zz = buffer_1100_ppsd[281];

    auto g_y_z_0_0_x_z_0_xx = buffer_1100_ppsd[282];

    auto g_y_z_0_0_x_z_0_xy = buffer_1100_ppsd[283];

    auto g_y_z_0_0_x_z_0_xz = buffer_1100_ppsd[284];

    auto g_y_z_0_0_x_z_0_yy = buffer_1100_ppsd[285];

    auto g_y_z_0_0_x_z_0_yz = buffer_1100_ppsd[286];

    auto g_y_z_0_0_x_z_0_zz = buffer_1100_ppsd[287];

    auto g_y_z_0_0_y_x_0_xx = buffer_1100_ppsd[288];

    auto g_y_z_0_0_y_x_0_xy = buffer_1100_ppsd[289];

    auto g_y_z_0_0_y_x_0_xz = buffer_1100_ppsd[290];

    auto g_y_z_0_0_y_x_0_yy = buffer_1100_ppsd[291];

    auto g_y_z_0_0_y_x_0_yz = buffer_1100_ppsd[292];

    auto g_y_z_0_0_y_x_0_zz = buffer_1100_ppsd[293];

    auto g_y_z_0_0_y_y_0_xx = buffer_1100_ppsd[294];

    auto g_y_z_0_0_y_y_0_xy = buffer_1100_ppsd[295];

    auto g_y_z_0_0_y_y_0_xz = buffer_1100_ppsd[296];

    auto g_y_z_0_0_y_y_0_yy = buffer_1100_ppsd[297];

    auto g_y_z_0_0_y_y_0_yz = buffer_1100_ppsd[298];

    auto g_y_z_0_0_y_y_0_zz = buffer_1100_ppsd[299];

    auto g_y_z_0_0_y_z_0_xx = buffer_1100_ppsd[300];

    auto g_y_z_0_0_y_z_0_xy = buffer_1100_ppsd[301];

    auto g_y_z_0_0_y_z_0_xz = buffer_1100_ppsd[302];

    auto g_y_z_0_0_y_z_0_yy = buffer_1100_ppsd[303];

    auto g_y_z_0_0_y_z_0_yz = buffer_1100_ppsd[304];

    auto g_y_z_0_0_y_z_0_zz = buffer_1100_ppsd[305];

    auto g_y_z_0_0_z_x_0_xx = buffer_1100_ppsd[306];

    auto g_y_z_0_0_z_x_0_xy = buffer_1100_ppsd[307];

    auto g_y_z_0_0_z_x_0_xz = buffer_1100_ppsd[308];

    auto g_y_z_0_0_z_x_0_yy = buffer_1100_ppsd[309];

    auto g_y_z_0_0_z_x_0_yz = buffer_1100_ppsd[310];

    auto g_y_z_0_0_z_x_0_zz = buffer_1100_ppsd[311];

    auto g_y_z_0_0_z_y_0_xx = buffer_1100_ppsd[312];

    auto g_y_z_0_0_z_y_0_xy = buffer_1100_ppsd[313];

    auto g_y_z_0_0_z_y_0_xz = buffer_1100_ppsd[314];

    auto g_y_z_0_0_z_y_0_yy = buffer_1100_ppsd[315];

    auto g_y_z_0_0_z_y_0_yz = buffer_1100_ppsd[316];

    auto g_y_z_0_0_z_y_0_zz = buffer_1100_ppsd[317];

    auto g_y_z_0_0_z_z_0_xx = buffer_1100_ppsd[318];

    auto g_y_z_0_0_z_z_0_xy = buffer_1100_ppsd[319];

    auto g_y_z_0_0_z_z_0_xz = buffer_1100_ppsd[320];

    auto g_y_z_0_0_z_z_0_yy = buffer_1100_ppsd[321];

    auto g_y_z_0_0_z_z_0_yz = buffer_1100_ppsd[322];

    auto g_y_z_0_0_z_z_0_zz = buffer_1100_ppsd[323];

    auto g_z_x_0_0_x_x_0_xx = buffer_1100_ppsd[324];

    auto g_z_x_0_0_x_x_0_xy = buffer_1100_ppsd[325];

    auto g_z_x_0_0_x_x_0_xz = buffer_1100_ppsd[326];

    auto g_z_x_0_0_x_x_0_yy = buffer_1100_ppsd[327];

    auto g_z_x_0_0_x_x_0_yz = buffer_1100_ppsd[328];

    auto g_z_x_0_0_x_x_0_zz = buffer_1100_ppsd[329];

    auto g_z_x_0_0_x_y_0_xx = buffer_1100_ppsd[330];

    auto g_z_x_0_0_x_y_0_xy = buffer_1100_ppsd[331];

    auto g_z_x_0_0_x_y_0_xz = buffer_1100_ppsd[332];

    auto g_z_x_0_0_x_y_0_yy = buffer_1100_ppsd[333];

    auto g_z_x_0_0_x_y_0_yz = buffer_1100_ppsd[334];

    auto g_z_x_0_0_x_y_0_zz = buffer_1100_ppsd[335];

    auto g_z_x_0_0_x_z_0_xx = buffer_1100_ppsd[336];

    auto g_z_x_0_0_x_z_0_xy = buffer_1100_ppsd[337];

    auto g_z_x_0_0_x_z_0_xz = buffer_1100_ppsd[338];

    auto g_z_x_0_0_x_z_0_yy = buffer_1100_ppsd[339];

    auto g_z_x_0_0_x_z_0_yz = buffer_1100_ppsd[340];

    auto g_z_x_0_0_x_z_0_zz = buffer_1100_ppsd[341];

    auto g_z_x_0_0_y_x_0_xx = buffer_1100_ppsd[342];

    auto g_z_x_0_0_y_x_0_xy = buffer_1100_ppsd[343];

    auto g_z_x_0_0_y_x_0_xz = buffer_1100_ppsd[344];

    auto g_z_x_0_0_y_x_0_yy = buffer_1100_ppsd[345];

    auto g_z_x_0_0_y_x_0_yz = buffer_1100_ppsd[346];

    auto g_z_x_0_0_y_x_0_zz = buffer_1100_ppsd[347];

    auto g_z_x_0_0_y_y_0_xx = buffer_1100_ppsd[348];

    auto g_z_x_0_0_y_y_0_xy = buffer_1100_ppsd[349];

    auto g_z_x_0_0_y_y_0_xz = buffer_1100_ppsd[350];

    auto g_z_x_0_0_y_y_0_yy = buffer_1100_ppsd[351];

    auto g_z_x_0_0_y_y_0_yz = buffer_1100_ppsd[352];

    auto g_z_x_0_0_y_y_0_zz = buffer_1100_ppsd[353];

    auto g_z_x_0_0_y_z_0_xx = buffer_1100_ppsd[354];

    auto g_z_x_0_0_y_z_0_xy = buffer_1100_ppsd[355];

    auto g_z_x_0_0_y_z_0_xz = buffer_1100_ppsd[356];

    auto g_z_x_0_0_y_z_0_yy = buffer_1100_ppsd[357];

    auto g_z_x_0_0_y_z_0_yz = buffer_1100_ppsd[358];

    auto g_z_x_0_0_y_z_0_zz = buffer_1100_ppsd[359];

    auto g_z_x_0_0_z_x_0_xx = buffer_1100_ppsd[360];

    auto g_z_x_0_0_z_x_0_xy = buffer_1100_ppsd[361];

    auto g_z_x_0_0_z_x_0_xz = buffer_1100_ppsd[362];

    auto g_z_x_0_0_z_x_0_yy = buffer_1100_ppsd[363];

    auto g_z_x_0_0_z_x_0_yz = buffer_1100_ppsd[364];

    auto g_z_x_0_0_z_x_0_zz = buffer_1100_ppsd[365];

    auto g_z_x_0_0_z_y_0_xx = buffer_1100_ppsd[366];

    auto g_z_x_0_0_z_y_0_xy = buffer_1100_ppsd[367];

    auto g_z_x_0_0_z_y_0_xz = buffer_1100_ppsd[368];

    auto g_z_x_0_0_z_y_0_yy = buffer_1100_ppsd[369];

    auto g_z_x_0_0_z_y_0_yz = buffer_1100_ppsd[370];

    auto g_z_x_0_0_z_y_0_zz = buffer_1100_ppsd[371];

    auto g_z_x_0_0_z_z_0_xx = buffer_1100_ppsd[372];

    auto g_z_x_0_0_z_z_0_xy = buffer_1100_ppsd[373];

    auto g_z_x_0_0_z_z_0_xz = buffer_1100_ppsd[374];

    auto g_z_x_0_0_z_z_0_yy = buffer_1100_ppsd[375];

    auto g_z_x_0_0_z_z_0_yz = buffer_1100_ppsd[376];

    auto g_z_x_0_0_z_z_0_zz = buffer_1100_ppsd[377];

    auto g_z_y_0_0_x_x_0_xx = buffer_1100_ppsd[378];

    auto g_z_y_0_0_x_x_0_xy = buffer_1100_ppsd[379];

    auto g_z_y_0_0_x_x_0_xz = buffer_1100_ppsd[380];

    auto g_z_y_0_0_x_x_0_yy = buffer_1100_ppsd[381];

    auto g_z_y_0_0_x_x_0_yz = buffer_1100_ppsd[382];

    auto g_z_y_0_0_x_x_0_zz = buffer_1100_ppsd[383];

    auto g_z_y_0_0_x_y_0_xx = buffer_1100_ppsd[384];

    auto g_z_y_0_0_x_y_0_xy = buffer_1100_ppsd[385];

    auto g_z_y_0_0_x_y_0_xz = buffer_1100_ppsd[386];

    auto g_z_y_0_0_x_y_0_yy = buffer_1100_ppsd[387];

    auto g_z_y_0_0_x_y_0_yz = buffer_1100_ppsd[388];

    auto g_z_y_0_0_x_y_0_zz = buffer_1100_ppsd[389];

    auto g_z_y_0_0_x_z_0_xx = buffer_1100_ppsd[390];

    auto g_z_y_0_0_x_z_0_xy = buffer_1100_ppsd[391];

    auto g_z_y_0_0_x_z_0_xz = buffer_1100_ppsd[392];

    auto g_z_y_0_0_x_z_0_yy = buffer_1100_ppsd[393];

    auto g_z_y_0_0_x_z_0_yz = buffer_1100_ppsd[394];

    auto g_z_y_0_0_x_z_0_zz = buffer_1100_ppsd[395];

    auto g_z_y_0_0_y_x_0_xx = buffer_1100_ppsd[396];

    auto g_z_y_0_0_y_x_0_xy = buffer_1100_ppsd[397];

    auto g_z_y_0_0_y_x_0_xz = buffer_1100_ppsd[398];

    auto g_z_y_0_0_y_x_0_yy = buffer_1100_ppsd[399];

    auto g_z_y_0_0_y_x_0_yz = buffer_1100_ppsd[400];

    auto g_z_y_0_0_y_x_0_zz = buffer_1100_ppsd[401];

    auto g_z_y_0_0_y_y_0_xx = buffer_1100_ppsd[402];

    auto g_z_y_0_0_y_y_0_xy = buffer_1100_ppsd[403];

    auto g_z_y_0_0_y_y_0_xz = buffer_1100_ppsd[404];

    auto g_z_y_0_0_y_y_0_yy = buffer_1100_ppsd[405];

    auto g_z_y_0_0_y_y_0_yz = buffer_1100_ppsd[406];

    auto g_z_y_0_0_y_y_0_zz = buffer_1100_ppsd[407];

    auto g_z_y_0_0_y_z_0_xx = buffer_1100_ppsd[408];

    auto g_z_y_0_0_y_z_0_xy = buffer_1100_ppsd[409];

    auto g_z_y_0_0_y_z_0_xz = buffer_1100_ppsd[410];

    auto g_z_y_0_0_y_z_0_yy = buffer_1100_ppsd[411];

    auto g_z_y_0_0_y_z_0_yz = buffer_1100_ppsd[412];

    auto g_z_y_0_0_y_z_0_zz = buffer_1100_ppsd[413];

    auto g_z_y_0_0_z_x_0_xx = buffer_1100_ppsd[414];

    auto g_z_y_0_0_z_x_0_xy = buffer_1100_ppsd[415];

    auto g_z_y_0_0_z_x_0_xz = buffer_1100_ppsd[416];

    auto g_z_y_0_0_z_x_0_yy = buffer_1100_ppsd[417];

    auto g_z_y_0_0_z_x_0_yz = buffer_1100_ppsd[418];

    auto g_z_y_0_0_z_x_0_zz = buffer_1100_ppsd[419];

    auto g_z_y_0_0_z_y_0_xx = buffer_1100_ppsd[420];

    auto g_z_y_0_0_z_y_0_xy = buffer_1100_ppsd[421];

    auto g_z_y_0_0_z_y_0_xz = buffer_1100_ppsd[422];

    auto g_z_y_0_0_z_y_0_yy = buffer_1100_ppsd[423];

    auto g_z_y_0_0_z_y_0_yz = buffer_1100_ppsd[424];

    auto g_z_y_0_0_z_y_0_zz = buffer_1100_ppsd[425];

    auto g_z_y_0_0_z_z_0_xx = buffer_1100_ppsd[426];

    auto g_z_y_0_0_z_z_0_xy = buffer_1100_ppsd[427];

    auto g_z_y_0_0_z_z_0_xz = buffer_1100_ppsd[428];

    auto g_z_y_0_0_z_z_0_yy = buffer_1100_ppsd[429];

    auto g_z_y_0_0_z_z_0_yz = buffer_1100_ppsd[430];

    auto g_z_y_0_0_z_z_0_zz = buffer_1100_ppsd[431];

    auto g_z_z_0_0_x_x_0_xx = buffer_1100_ppsd[432];

    auto g_z_z_0_0_x_x_0_xy = buffer_1100_ppsd[433];

    auto g_z_z_0_0_x_x_0_xz = buffer_1100_ppsd[434];

    auto g_z_z_0_0_x_x_0_yy = buffer_1100_ppsd[435];

    auto g_z_z_0_0_x_x_0_yz = buffer_1100_ppsd[436];

    auto g_z_z_0_0_x_x_0_zz = buffer_1100_ppsd[437];

    auto g_z_z_0_0_x_y_0_xx = buffer_1100_ppsd[438];

    auto g_z_z_0_0_x_y_0_xy = buffer_1100_ppsd[439];

    auto g_z_z_0_0_x_y_0_xz = buffer_1100_ppsd[440];

    auto g_z_z_0_0_x_y_0_yy = buffer_1100_ppsd[441];

    auto g_z_z_0_0_x_y_0_yz = buffer_1100_ppsd[442];

    auto g_z_z_0_0_x_y_0_zz = buffer_1100_ppsd[443];

    auto g_z_z_0_0_x_z_0_xx = buffer_1100_ppsd[444];

    auto g_z_z_0_0_x_z_0_xy = buffer_1100_ppsd[445];

    auto g_z_z_0_0_x_z_0_xz = buffer_1100_ppsd[446];

    auto g_z_z_0_0_x_z_0_yy = buffer_1100_ppsd[447];

    auto g_z_z_0_0_x_z_0_yz = buffer_1100_ppsd[448];

    auto g_z_z_0_0_x_z_0_zz = buffer_1100_ppsd[449];

    auto g_z_z_0_0_y_x_0_xx = buffer_1100_ppsd[450];

    auto g_z_z_0_0_y_x_0_xy = buffer_1100_ppsd[451];

    auto g_z_z_0_0_y_x_0_xz = buffer_1100_ppsd[452];

    auto g_z_z_0_0_y_x_0_yy = buffer_1100_ppsd[453];

    auto g_z_z_0_0_y_x_0_yz = buffer_1100_ppsd[454];

    auto g_z_z_0_0_y_x_0_zz = buffer_1100_ppsd[455];

    auto g_z_z_0_0_y_y_0_xx = buffer_1100_ppsd[456];

    auto g_z_z_0_0_y_y_0_xy = buffer_1100_ppsd[457];

    auto g_z_z_0_0_y_y_0_xz = buffer_1100_ppsd[458];

    auto g_z_z_0_0_y_y_0_yy = buffer_1100_ppsd[459];

    auto g_z_z_0_0_y_y_0_yz = buffer_1100_ppsd[460];

    auto g_z_z_0_0_y_y_0_zz = buffer_1100_ppsd[461];

    auto g_z_z_0_0_y_z_0_xx = buffer_1100_ppsd[462];

    auto g_z_z_0_0_y_z_0_xy = buffer_1100_ppsd[463];

    auto g_z_z_0_0_y_z_0_xz = buffer_1100_ppsd[464];

    auto g_z_z_0_0_y_z_0_yy = buffer_1100_ppsd[465];

    auto g_z_z_0_0_y_z_0_yz = buffer_1100_ppsd[466];

    auto g_z_z_0_0_y_z_0_zz = buffer_1100_ppsd[467];

    auto g_z_z_0_0_z_x_0_xx = buffer_1100_ppsd[468];

    auto g_z_z_0_0_z_x_0_xy = buffer_1100_ppsd[469];

    auto g_z_z_0_0_z_x_0_xz = buffer_1100_ppsd[470];

    auto g_z_z_0_0_z_x_0_yy = buffer_1100_ppsd[471];

    auto g_z_z_0_0_z_x_0_yz = buffer_1100_ppsd[472];

    auto g_z_z_0_0_z_x_0_zz = buffer_1100_ppsd[473];

    auto g_z_z_0_0_z_y_0_xx = buffer_1100_ppsd[474];

    auto g_z_z_0_0_z_y_0_xy = buffer_1100_ppsd[475];

    auto g_z_z_0_0_z_y_0_xz = buffer_1100_ppsd[476];

    auto g_z_z_0_0_z_y_0_yy = buffer_1100_ppsd[477];

    auto g_z_z_0_0_z_y_0_yz = buffer_1100_ppsd[478];

    auto g_z_z_0_0_z_y_0_zz = buffer_1100_ppsd[479];

    auto g_z_z_0_0_z_z_0_xx = buffer_1100_ppsd[480];

    auto g_z_z_0_0_z_z_0_xy = buffer_1100_ppsd[481];

    auto g_z_z_0_0_z_z_0_xz = buffer_1100_ppsd[482];

    auto g_z_z_0_0_z_z_0_yy = buffer_1100_ppsd[483];

    auto g_z_z_0_0_z_z_0_yz = buffer_1100_ppsd[484];

    auto g_z_z_0_0_z_z_0_zz = buffer_1100_ppsd[485];

    // integrals block (0-6)

    #pragma omp simd aligned(g_0_0_0_xx, g_0_0_0_xy, g_0_0_0_xz, g_0_0_0_yy, g_0_0_0_yz, g_0_0_0_zz, g_0_xx_0_xx, g_0_xx_0_xy, g_0_xx_0_xz, g_0_xx_0_yy, g_0_xx_0_yz, g_0_xx_0_zz, g_x_x_0_0_x_x_0_xx, g_x_x_0_0_x_x_0_xy, g_x_x_0_0_x_x_0_xz, g_x_x_0_0_x_x_0_yy, g_x_x_0_0_x_x_0_yz, g_x_x_0_0_x_x_0_zz, g_xx_0_0_xx, g_xx_0_0_xy, g_xx_0_0_xz, g_xx_0_0_yy, g_xx_0_0_yz, g_xx_0_0_zz, g_xx_xx_0_xx, g_xx_xx_0_xy, g_xx_xx_0_xz, g_xx_xx_0_yy, g_xx_xx_0_yz, g_xx_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_x_0_xx[i] = g_0_0_0_xx[i] - 2.0 * g_0_xx_0_xx[i] * b_exp - 2.0 * g_xx_0_0_xx[i] * a_exp + 4.0 * g_xx_xx_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_0_xy[i] = g_0_0_0_xy[i] - 2.0 * g_0_xx_0_xy[i] * b_exp - 2.0 * g_xx_0_0_xy[i] * a_exp + 4.0 * g_xx_xx_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_0_xz[i] = g_0_0_0_xz[i] - 2.0 * g_0_xx_0_xz[i] * b_exp - 2.0 * g_xx_0_0_xz[i] * a_exp + 4.0 * g_xx_xx_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_0_yy[i] = g_0_0_0_yy[i] - 2.0 * g_0_xx_0_yy[i] * b_exp - 2.0 * g_xx_0_0_yy[i] * a_exp + 4.0 * g_xx_xx_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_0_yz[i] = g_0_0_0_yz[i] - 2.0 * g_0_xx_0_yz[i] * b_exp - 2.0 * g_xx_0_0_yz[i] * a_exp + 4.0 * g_xx_xx_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_0_zz[i] = g_0_0_0_zz[i] - 2.0 * g_0_xx_0_zz[i] * b_exp - 2.0 * g_xx_0_0_zz[i] * a_exp + 4.0 * g_xx_xx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_x_x_0_0_x_y_0_xx, g_x_x_0_0_x_y_0_xy, g_x_x_0_0_x_y_0_xz, g_x_x_0_0_x_y_0_yy, g_x_x_0_0_x_y_0_yz, g_x_x_0_0_x_y_0_zz, g_xx_xy_0_xx, g_xx_xy_0_xy, g_xx_xy_0_xz, g_xx_xy_0_yy, g_xx_xy_0_yz, g_xx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_y_0_xx[i] = -2.0 * g_0_xy_0_xx[i] * b_exp + 4.0 * g_xx_xy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_0_xy[i] = -2.0 * g_0_xy_0_xy[i] * b_exp + 4.0 * g_xx_xy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_0_xz[i] = -2.0 * g_0_xy_0_xz[i] * b_exp + 4.0 * g_xx_xy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_0_yy[i] = -2.0 * g_0_xy_0_yy[i] * b_exp + 4.0 * g_xx_xy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_0_yz[i] = -2.0 * g_0_xy_0_yz[i] * b_exp + 4.0 * g_xx_xy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_0_zz[i] = -2.0 * g_0_xy_0_zz[i] * b_exp + 4.0 * g_xx_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_x_x_0_0_x_z_0_xx, g_x_x_0_0_x_z_0_xy, g_x_x_0_0_x_z_0_xz, g_x_x_0_0_x_z_0_yy, g_x_x_0_0_x_z_0_yz, g_x_x_0_0_x_z_0_zz, g_xx_xz_0_xx, g_xx_xz_0_xy, g_xx_xz_0_xz, g_xx_xz_0_yy, g_xx_xz_0_yz, g_xx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_z_0_xx[i] = -2.0 * g_0_xz_0_xx[i] * b_exp + 4.0 * g_xx_xz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_0_xy[i] = -2.0 * g_0_xz_0_xy[i] * b_exp + 4.0 * g_xx_xz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_0_xz[i] = -2.0 * g_0_xz_0_xz[i] * b_exp + 4.0 * g_xx_xz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_0_yy[i] = -2.0 * g_0_xz_0_yy[i] * b_exp + 4.0 * g_xx_xz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_0_yz[i] = -2.0 * g_0_xz_0_yz[i] * b_exp + 4.0 * g_xx_xz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_0_zz[i] = -2.0 * g_0_xz_0_zz[i] * b_exp + 4.0 * g_xx_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_x_0_0_y_x_0_xx, g_x_x_0_0_y_x_0_xy, g_x_x_0_0_y_x_0_xz, g_x_x_0_0_y_x_0_yy, g_x_x_0_0_y_x_0_yz, g_x_x_0_0_y_x_0_zz, g_xy_0_0_xx, g_xy_0_0_xy, g_xy_0_0_xz, g_xy_0_0_yy, g_xy_0_0_yz, g_xy_0_0_zz, g_xy_xx_0_xx, g_xy_xx_0_xy, g_xy_xx_0_xz, g_xy_xx_0_yy, g_xy_xx_0_yz, g_xy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_x_0_xx[i] = -2.0 * g_xy_0_0_xx[i] * a_exp + 4.0 * g_xy_xx_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_0_xy[i] = -2.0 * g_xy_0_0_xy[i] * a_exp + 4.0 * g_xy_xx_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_0_xz[i] = -2.0 * g_xy_0_0_xz[i] * a_exp + 4.0 * g_xy_xx_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_0_yy[i] = -2.0 * g_xy_0_0_yy[i] * a_exp + 4.0 * g_xy_xx_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_0_yz[i] = -2.0 * g_xy_0_0_yz[i] * a_exp + 4.0 * g_xy_xx_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_0_zz[i] = -2.0 * g_xy_0_0_zz[i] * a_exp + 4.0 * g_xy_xx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_x_0_0_y_y_0_xx, g_x_x_0_0_y_y_0_xy, g_x_x_0_0_y_y_0_xz, g_x_x_0_0_y_y_0_yy, g_x_x_0_0_y_y_0_yz, g_x_x_0_0_y_y_0_zz, g_xy_xy_0_xx, g_xy_xy_0_xy, g_xy_xy_0_xz, g_xy_xy_0_yy, g_xy_xy_0_yz, g_xy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_y_0_xx[i] = 4.0 * g_xy_xy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_0_xy[i] = 4.0 * g_xy_xy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_0_xz[i] = 4.0 * g_xy_xy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_0_yy[i] = 4.0 * g_xy_xy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_0_yz[i] = 4.0 * g_xy_xy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_0_zz[i] = 4.0 * g_xy_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_x_0_0_y_z_0_xx, g_x_x_0_0_y_z_0_xy, g_x_x_0_0_y_z_0_xz, g_x_x_0_0_y_z_0_yy, g_x_x_0_0_y_z_0_yz, g_x_x_0_0_y_z_0_zz, g_xy_xz_0_xx, g_xy_xz_0_xy, g_xy_xz_0_xz, g_xy_xz_0_yy, g_xy_xz_0_yz, g_xy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_z_0_xx[i] = 4.0 * g_xy_xz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_0_xy[i] = 4.0 * g_xy_xz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_0_xz[i] = 4.0 * g_xy_xz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_0_yy[i] = 4.0 * g_xy_xz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_0_yz[i] = 4.0 * g_xy_xz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_0_zz[i] = 4.0 * g_xy_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_x_0_0_z_x_0_xx, g_x_x_0_0_z_x_0_xy, g_x_x_0_0_z_x_0_xz, g_x_x_0_0_z_x_0_yy, g_x_x_0_0_z_x_0_yz, g_x_x_0_0_z_x_0_zz, g_xz_0_0_xx, g_xz_0_0_xy, g_xz_0_0_xz, g_xz_0_0_yy, g_xz_0_0_yz, g_xz_0_0_zz, g_xz_xx_0_xx, g_xz_xx_0_xy, g_xz_xx_0_xz, g_xz_xx_0_yy, g_xz_xx_0_yz, g_xz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_x_0_xx[i] = -2.0 * g_xz_0_0_xx[i] * a_exp + 4.0 * g_xz_xx_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_0_xy[i] = -2.0 * g_xz_0_0_xy[i] * a_exp + 4.0 * g_xz_xx_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_0_xz[i] = -2.0 * g_xz_0_0_xz[i] * a_exp + 4.0 * g_xz_xx_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_0_yy[i] = -2.0 * g_xz_0_0_yy[i] * a_exp + 4.0 * g_xz_xx_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_0_yz[i] = -2.0 * g_xz_0_0_yz[i] * a_exp + 4.0 * g_xz_xx_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_0_zz[i] = -2.0 * g_xz_0_0_zz[i] * a_exp + 4.0 * g_xz_xx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_x_0_0_z_y_0_xx, g_x_x_0_0_z_y_0_xy, g_x_x_0_0_z_y_0_xz, g_x_x_0_0_z_y_0_yy, g_x_x_0_0_z_y_0_yz, g_x_x_0_0_z_y_0_zz, g_xz_xy_0_xx, g_xz_xy_0_xy, g_xz_xy_0_xz, g_xz_xy_0_yy, g_xz_xy_0_yz, g_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_y_0_xx[i] = 4.0 * g_xz_xy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_0_xy[i] = 4.0 * g_xz_xy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_0_xz[i] = 4.0 * g_xz_xy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_0_yy[i] = 4.0 * g_xz_xy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_0_yz[i] = 4.0 * g_xz_xy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_0_zz[i] = 4.0 * g_xz_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_x_0_0_z_z_0_xx, g_x_x_0_0_z_z_0_xy, g_x_x_0_0_z_z_0_xz, g_x_x_0_0_z_z_0_yy, g_x_x_0_0_z_z_0_yz, g_x_x_0_0_z_z_0_zz, g_xz_xz_0_xx, g_xz_xz_0_xy, g_xz_xz_0_xz, g_xz_xz_0_yy, g_xz_xz_0_yz, g_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_z_0_xx[i] = 4.0 * g_xz_xz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_0_xy[i] = 4.0 * g_xz_xz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_0_xz[i] = 4.0 * g_xz_xz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_0_yy[i] = 4.0 * g_xz_xz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_0_yz[i] = 4.0 * g_xz_xz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_0_zz[i] = 4.0 * g_xz_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_x_y_0_0_x_x_0_xx, g_x_y_0_0_x_x_0_xy, g_x_y_0_0_x_x_0_xz, g_x_y_0_0_x_x_0_yy, g_x_y_0_0_x_x_0_yz, g_x_y_0_0_x_x_0_zz, g_xx_xy_0_xx, g_xx_xy_0_xy, g_xx_xy_0_xz, g_xx_xy_0_yy, g_xx_xy_0_yz, g_xx_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_x_0_xx[i] = -2.0 * g_0_xy_0_xx[i] * b_exp + 4.0 * g_xx_xy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_0_xy[i] = -2.0 * g_0_xy_0_xy[i] * b_exp + 4.0 * g_xx_xy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_0_xz[i] = -2.0 * g_0_xy_0_xz[i] * b_exp + 4.0 * g_xx_xy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_0_yy[i] = -2.0 * g_0_xy_0_yy[i] * b_exp + 4.0 * g_xx_xy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_0_yz[i] = -2.0 * g_0_xy_0_yz[i] * b_exp + 4.0 * g_xx_xy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_0_zz[i] = -2.0 * g_0_xy_0_zz[i] * b_exp + 4.0 * g_xx_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_0_0_0_xx, g_0_0_0_xy, g_0_0_0_xz, g_0_0_0_yy, g_0_0_0_yz, g_0_0_0_zz, g_0_yy_0_xx, g_0_yy_0_xy, g_0_yy_0_xz, g_0_yy_0_yy, g_0_yy_0_yz, g_0_yy_0_zz, g_x_y_0_0_x_y_0_xx, g_x_y_0_0_x_y_0_xy, g_x_y_0_0_x_y_0_xz, g_x_y_0_0_x_y_0_yy, g_x_y_0_0_x_y_0_yz, g_x_y_0_0_x_y_0_zz, g_xx_0_0_xx, g_xx_0_0_xy, g_xx_0_0_xz, g_xx_0_0_yy, g_xx_0_0_yz, g_xx_0_0_zz, g_xx_yy_0_xx, g_xx_yy_0_xy, g_xx_yy_0_xz, g_xx_yy_0_yy, g_xx_yy_0_yz, g_xx_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_y_0_xx[i] = g_0_0_0_xx[i] - 2.0 * g_0_yy_0_xx[i] * b_exp - 2.0 * g_xx_0_0_xx[i] * a_exp + 4.0 * g_xx_yy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_0_xy[i] = g_0_0_0_xy[i] - 2.0 * g_0_yy_0_xy[i] * b_exp - 2.0 * g_xx_0_0_xy[i] * a_exp + 4.0 * g_xx_yy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_0_xz[i] = g_0_0_0_xz[i] - 2.0 * g_0_yy_0_xz[i] * b_exp - 2.0 * g_xx_0_0_xz[i] * a_exp + 4.0 * g_xx_yy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_0_yy[i] = g_0_0_0_yy[i] - 2.0 * g_0_yy_0_yy[i] * b_exp - 2.0 * g_xx_0_0_yy[i] * a_exp + 4.0 * g_xx_yy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_0_yz[i] = g_0_0_0_yz[i] - 2.0 * g_0_yy_0_yz[i] * b_exp - 2.0 * g_xx_0_0_yz[i] * a_exp + 4.0 * g_xx_yy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_0_zz[i] = g_0_0_0_zz[i] - 2.0 * g_0_yy_0_zz[i] * b_exp - 2.0 * g_xx_0_0_zz[i] * a_exp + 4.0 * g_xx_yy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_x_y_0_0_x_z_0_xx, g_x_y_0_0_x_z_0_xy, g_x_y_0_0_x_z_0_xz, g_x_y_0_0_x_z_0_yy, g_x_y_0_0_x_z_0_yz, g_x_y_0_0_x_z_0_zz, g_xx_yz_0_xx, g_xx_yz_0_xy, g_xx_yz_0_xz, g_xx_yz_0_yy, g_xx_yz_0_yz, g_xx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_z_0_xx[i] = -2.0 * g_0_yz_0_xx[i] * b_exp + 4.0 * g_xx_yz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_0_xy[i] = -2.0 * g_0_yz_0_xy[i] * b_exp + 4.0 * g_xx_yz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_0_xz[i] = -2.0 * g_0_yz_0_xz[i] * b_exp + 4.0 * g_xx_yz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_0_yy[i] = -2.0 * g_0_yz_0_yy[i] * b_exp + 4.0 * g_xx_yz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_0_yz[i] = -2.0 * g_0_yz_0_yz[i] * b_exp + 4.0 * g_xx_yz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_0_zz[i] = -2.0 * g_0_yz_0_zz[i] * b_exp + 4.0 * g_xx_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_x_y_0_0_y_x_0_xx, g_x_y_0_0_y_x_0_xy, g_x_y_0_0_y_x_0_xz, g_x_y_0_0_y_x_0_yy, g_x_y_0_0_y_x_0_yz, g_x_y_0_0_y_x_0_zz, g_xy_xy_0_xx, g_xy_xy_0_xy, g_xy_xy_0_xz, g_xy_xy_0_yy, g_xy_xy_0_yz, g_xy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_x_0_xx[i] = 4.0 * g_xy_xy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_0_xy[i] = 4.0 * g_xy_xy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_0_xz[i] = 4.0 * g_xy_xy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_0_yy[i] = 4.0 * g_xy_xy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_0_yz[i] = 4.0 * g_xy_xy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_0_zz[i] = 4.0 * g_xy_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_x_y_0_0_y_y_0_xx, g_x_y_0_0_y_y_0_xy, g_x_y_0_0_y_y_0_xz, g_x_y_0_0_y_y_0_yy, g_x_y_0_0_y_y_0_yz, g_x_y_0_0_y_y_0_zz, g_xy_0_0_xx, g_xy_0_0_xy, g_xy_0_0_xz, g_xy_0_0_yy, g_xy_0_0_yz, g_xy_0_0_zz, g_xy_yy_0_xx, g_xy_yy_0_xy, g_xy_yy_0_xz, g_xy_yy_0_yy, g_xy_yy_0_yz, g_xy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_y_0_xx[i] = -2.0 * g_xy_0_0_xx[i] * a_exp + 4.0 * g_xy_yy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_0_xy[i] = -2.0 * g_xy_0_0_xy[i] * a_exp + 4.0 * g_xy_yy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_0_xz[i] = -2.0 * g_xy_0_0_xz[i] * a_exp + 4.0 * g_xy_yy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_0_yy[i] = -2.0 * g_xy_0_0_yy[i] * a_exp + 4.0 * g_xy_yy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_0_yz[i] = -2.0 * g_xy_0_0_yz[i] * a_exp + 4.0 * g_xy_yy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_0_zz[i] = -2.0 * g_xy_0_0_zz[i] * a_exp + 4.0 * g_xy_yy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_y_0_0_y_z_0_xx, g_x_y_0_0_y_z_0_xy, g_x_y_0_0_y_z_0_xz, g_x_y_0_0_y_z_0_yy, g_x_y_0_0_y_z_0_yz, g_x_y_0_0_y_z_0_zz, g_xy_yz_0_xx, g_xy_yz_0_xy, g_xy_yz_0_xz, g_xy_yz_0_yy, g_xy_yz_0_yz, g_xy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_z_0_xx[i] = 4.0 * g_xy_yz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_0_xy[i] = 4.0 * g_xy_yz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_0_xz[i] = 4.0 * g_xy_yz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_0_yy[i] = 4.0 * g_xy_yz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_0_yz[i] = 4.0 * g_xy_yz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_0_zz[i] = 4.0 * g_xy_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_y_0_0_z_x_0_xx, g_x_y_0_0_z_x_0_xy, g_x_y_0_0_z_x_0_xz, g_x_y_0_0_z_x_0_yy, g_x_y_0_0_z_x_0_yz, g_x_y_0_0_z_x_0_zz, g_xz_xy_0_xx, g_xz_xy_0_xy, g_xz_xy_0_xz, g_xz_xy_0_yy, g_xz_xy_0_yz, g_xz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_x_0_xx[i] = 4.0 * g_xz_xy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_0_xy[i] = 4.0 * g_xz_xy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_0_xz[i] = 4.0 * g_xz_xy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_0_yy[i] = 4.0 * g_xz_xy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_0_yz[i] = 4.0 * g_xz_xy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_0_zz[i] = 4.0 * g_xz_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_x_y_0_0_z_y_0_xx, g_x_y_0_0_z_y_0_xy, g_x_y_0_0_z_y_0_xz, g_x_y_0_0_z_y_0_yy, g_x_y_0_0_z_y_0_yz, g_x_y_0_0_z_y_0_zz, g_xz_0_0_xx, g_xz_0_0_xy, g_xz_0_0_xz, g_xz_0_0_yy, g_xz_0_0_yz, g_xz_0_0_zz, g_xz_yy_0_xx, g_xz_yy_0_xy, g_xz_yy_0_xz, g_xz_yy_0_yy, g_xz_yy_0_yz, g_xz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_y_0_xx[i] = -2.0 * g_xz_0_0_xx[i] * a_exp + 4.0 * g_xz_yy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_0_xy[i] = -2.0 * g_xz_0_0_xy[i] * a_exp + 4.0 * g_xz_yy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_0_xz[i] = -2.0 * g_xz_0_0_xz[i] * a_exp + 4.0 * g_xz_yy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_0_yy[i] = -2.0 * g_xz_0_0_yy[i] * a_exp + 4.0 * g_xz_yy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_0_yz[i] = -2.0 * g_xz_0_0_yz[i] * a_exp + 4.0 * g_xz_yy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_0_zz[i] = -2.0 * g_xz_0_0_zz[i] * a_exp + 4.0 * g_xz_yy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_x_y_0_0_z_z_0_xx, g_x_y_0_0_z_z_0_xy, g_x_y_0_0_z_z_0_xz, g_x_y_0_0_z_z_0_yy, g_x_y_0_0_z_z_0_yz, g_x_y_0_0_z_z_0_zz, g_xz_yz_0_xx, g_xz_yz_0_xy, g_xz_yz_0_xz, g_xz_yz_0_yy, g_xz_yz_0_yz, g_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_z_0_xx[i] = 4.0 * g_xz_yz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_0_xy[i] = 4.0 * g_xz_yz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_0_xz[i] = 4.0 * g_xz_yz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_0_yy[i] = 4.0 * g_xz_yz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_0_yz[i] = 4.0 * g_xz_yz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_0_zz[i] = 4.0 * g_xz_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_x_z_0_0_x_x_0_xx, g_x_z_0_0_x_x_0_xy, g_x_z_0_0_x_x_0_xz, g_x_z_0_0_x_x_0_yy, g_x_z_0_0_x_x_0_yz, g_x_z_0_0_x_x_0_zz, g_xx_xz_0_xx, g_xx_xz_0_xy, g_xx_xz_0_xz, g_xx_xz_0_yy, g_xx_xz_0_yz, g_xx_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_x_0_xx[i] = -2.0 * g_0_xz_0_xx[i] * b_exp + 4.0 * g_xx_xz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_0_xy[i] = -2.0 * g_0_xz_0_xy[i] * b_exp + 4.0 * g_xx_xz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_0_xz[i] = -2.0 * g_0_xz_0_xz[i] * b_exp + 4.0 * g_xx_xz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_0_yy[i] = -2.0 * g_0_xz_0_yy[i] * b_exp + 4.0 * g_xx_xz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_0_yz[i] = -2.0 * g_0_xz_0_yz[i] * b_exp + 4.0 * g_xx_xz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_0_zz[i] = -2.0 * g_0_xz_0_zz[i] * b_exp + 4.0 * g_xx_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_x_z_0_0_x_y_0_xx, g_x_z_0_0_x_y_0_xy, g_x_z_0_0_x_y_0_xz, g_x_z_0_0_x_y_0_yy, g_x_z_0_0_x_y_0_yz, g_x_z_0_0_x_y_0_zz, g_xx_yz_0_xx, g_xx_yz_0_xy, g_xx_yz_0_xz, g_xx_yz_0_yy, g_xx_yz_0_yz, g_xx_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_y_0_xx[i] = -2.0 * g_0_yz_0_xx[i] * b_exp + 4.0 * g_xx_yz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_0_xy[i] = -2.0 * g_0_yz_0_xy[i] * b_exp + 4.0 * g_xx_yz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_0_xz[i] = -2.0 * g_0_yz_0_xz[i] * b_exp + 4.0 * g_xx_yz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_0_yy[i] = -2.0 * g_0_yz_0_yy[i] * b_exp + 4.0 * g_xx_yz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_0_yz[i] = -2.0 * g_0_yz_0_yz[i] * b_exp + 4.0 * g_xx_yz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_0_zz[i] = -2.0 * g_0_yz_0_zz[i] * b_exp + 4.0 * g_xx_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_0_0_0_xx, g_0_0_0_xy, g_0_0_0_xz, g_0_0_0_yy, g_0_0_0_yz, g_0_0_0_zz, g_0_zz_0_xx, g_0_zz_0_xy, g_0_zz_0_xz, g_0_zz_0_yy, g_0_zz_0_yz, g_0_zz_0_zz, g_x_z_0_0_x_z_0_xx, g_x_z_0_0_x_z_0_xy, g_x_z_0_0_x_z_0_xz, g_x_z_0_0_x_z_0_yy, g_x_z_0_0_x_z_0_yz, g_x_z_0_0_x_z_0_zz, g_xx_0_0_xx, g_xx_0_0_xy, g_xx_0_0_xz, g_xx_0_0_yy, g_xx_0_0_yz, g_xx_0_0_zz, g_xx_zz_0_xx, g_xx_zz_0_xy, g_xx_zz_0_xz, g_xx_zz_0_yy, g_xx_zz_0_yz, g_xx_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_z_0_xx[i] = g_0_0_0_xx[i] - 2.0 * g_0_zz_0_xx[i] * b_exp - 2.0 * g_xx_0_0_xx[i] * a_exp + 4.0 * g_xx_zz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_0_xy[i] = g_0_0_0_xy[i] - 2.0 * g_0_zz_0_xy[i] * b_exp - 2.0 * g_xx_0_0_xy[i] * a_exp + 4.0 * g_xx_zz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_0_xz[i] = g_0_0_0_xz[i] - 2.0 * g_0_zz_0_xz[i] * b_exp - 2.0 * g_xx_0_0_xz[i] * a_exp + 4.0 * g_xx_zz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_0_yy[i] = g_0_0_0_yy[i] - 2.0 * g_0_zz_0_yy[i] * b_exp - 2.0 * g_xx_0_0_yy[i] * a_exp + 4.0 * g_xx_zz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_0_yz[i] = g_0_0_0_yz[i] - 2.0 * g_0_zz_0_yz[i] * b_exp - 2.0 * g_xx_0_0_yz[i] * a_exp + 4.0 * g_xx_zz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_0_zz[i] = g_0_0_0_zz[i] - 2.0 * g_0_zz_0_zz[i] * b_exp - 2.0 * g_xx_0_0_zz[i] * a_exp + 4.0 * g_xx_zz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_x_z_0_0_y_x_0_xx, g_x_z_0_0_y_x_0_xy, g_x_z_0_0_y_x_0_xz, g_x_z_0_0_y_x_0_yy, g_x_z_0_0_y_x_0_yz, g_x_z_0_0_y_x_0_zz, g_xy_xz_0_xx, g_xy_xz_0_xy, g_xy_xz_0_xz, g_xy_xz_0_yy, g_xy_xz_0_yz, g_xy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_x_0_xx[i] = 4.0 * g_xy_xz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_0_xy[i] = 4.0 * g_xy_xz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_0_xz[i] = 4.0 * g_xy_xz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_0_yy[i] = 4.0 * g_xy_xz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_0_yz[i] = 4.0 * g_xy_xz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_0_zz[i] = 4.0 * g_xy_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_x_z_0_0_y_y_0_xx, g_x_z_0_0_y_y_0_xy, g_x_z_0_0_y_y_0_xz, g_x_z_0_0_y_y_0_yy, g_x_z_0_0_y_y_0_yz, g_x_z_0_0_y_y_0_zz, g_xy_yz_0_xx, g_xy_yz_0_xy, g_xy_yz_0_xz, g_xy_yz_0_yy, g_xy_yz_0_yz, g_xy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_y_0_xx[i] = 4.0 * g_xy_yz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_0_xy[i] = 4.0 * g_xy_yz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_0_xz[i] = 4.0 * g_xy_yz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_0_yy[i] = 4.0 * g_xy_yz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_0_yz[i] = 4.0 * g_xy_yz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_0_zz[i] = 4.0 * g_xy_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_x_z_0_0_y_z_0_xx, g_x_z_0_0_y_z_0_xy, g_x_z_0_0_y_z_0_xz, g_x_z_0_0_y_z_0_yy, g_x_z_0_0_y_z_0_yz, g_x_z_0_0_y_z_0_zz, g_xy_0_0_xx, g_xy_0_0_xy, g_xy_0_0_xz, g_xy_0_0_yy, g_xy_0_0_yz, g_xy_0_0_zz, g_xy_zz_0_xx, g_xy_zz_0_xy, g_xy_zz_0_xz, g_xy_zz_0_yy, g_xy_zz_0_yz, g_xy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_z_0_xx[i] = -2.0 * g_xy_0_0_xx[i] * a_exp + 4.0 * g_xy_zz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_0_xy[i] = -2.0 * g_xy_0_0_xy[i] * a_exp + 4.0 * g_xy_zz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_0_xz[i] = -2.0 * g_xy_0_0_xz[i] * a_exp + 4.0 * g_xy_zz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_0_yy[i] = -2.0 * g_xy_0_0_yy[i] * a_exp + 4.0 * g_xy_zz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_0_yz[i] = -2.0 * g_xy_0_0_yz[i] * a_exp + 4.0 * g_xy_zz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_0_zz[i] = -2.0 * g_xy_0_0_zz[i] * a_exp + 4.0 * g_xy_zz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_x_z_0_0_z_x_0_xx, g_x_z_0_0_z_x_0_xy, g_x_z_0_0_z_x_0_xz, g_x_z_0_0_z_x_0_yy, g_x_z_0_0_z_x_0_yz, g_x_z_0_0_z_x_0_zz, g_xz_xz_0_xx, g_xz_xz_0_xy, g_xz_xz_0_xz, g_xz_xz_0_yy, g_xz_xz_0_yz, g_xz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_x_0_xx[i] = 4.0 * g_xz_xz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_0_xy[i] = 4.0 * g_xz_xz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_0_xz[i] = 4.0 * g_xz_xz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_0_yy[i] = 4.0 * g_xz_xz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_0_yz[i] = 4.0 * g_xz_xz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_0_zz[i] = 4.0 * g_xz_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_x_z_0_0_z_y_0_xx, g_x_z_0_0_z_y_0_xy, g_x_z_0_0_z_y_0_xz, g_x_z_0_0_z_y_0_yy, g_x_z_0_0_z_y_0_yz, g_x_z_0_0_z_y_0_zz, g_xz_yz_0_xx, g_xz_yz_0_xy, g_xz_yz_0_xz, g_xz_yz_0_yy, g_xz_yz_0_yz, g_xz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_y_0_xx[i] = 4.0 * g_xz_yz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_0_xy[i] = 4.0 * g_xz_yz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_0_xz[i] = 4.0 * g_xz_yz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_0_yy[i] = 4.0 * g_xz_yz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_0_yz[i] = 4.0 * g_xz_yz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_0_zz[i] = 4.0 * g_xz_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_x_z_0_0_z_z_0_xx, g_x_z_0_0_z_z_0_xy, g_x_z_0_0_z_z_0_xz, g_x_z_0_0_z_z_0_yy, g_x_z_0_0_z_z_0_yz, g_x_z_0_0_z_z_0_zz, g_xz_0_0_xx, g_xz_0_0_xy, g_xz_0_0_xz, g_xz_0_0_yy, g_xz_0_0_yz, g_xz_0_0_zz, g_xz_zz_0_xx, g_xz_zz_0_xy, g_xz_zz_0_xz, g_xz_zz_0_yy, g_xz_zz_0_yz, g_xz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_z_0_xx[i] = -2.0 * g_xz_0_0_xx[i] * a_exp + 4.0 * g_xz_zz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_0_xy[i] = -2.0 * g_xz_0_0_xy[i] * a_exp + 4.0 * g_xz_zz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_0_xz[i] = -2.0 * g_xz_0_0_xz[i] * a_exp + 4.0 * g_xz_zz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_0_yy[i] = -2.0 * g_xz_0_0_yy[i] * a_exp + 4.0 * g_xz_zz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_0_yz[i] = -2.0 * g_xz_0_0_yz[i] * a_exp + 4.0 * g_xz_zz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_0_zz[i] = -2.0 * g_xz_0_0_zz[i] * a_exp + 4.0 * g_xz_zz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (162-168)

    #pragma omp simd aligned(g_xy_0_0_xx, g_xy_0_0_xy, g_xy_0_0_xz, g_xy_0_0_yy, g_xy_0_0_yz, g_xy_0_0_zz, g_xy_xx_0_xx, g_xy_xx_0_xy, g_xy_xx_0_xz, g_xy_xx_0_yy, g_xy_xx_0_yz, g_xy_xx_0_zz, g_y_x_0_0_x_x_0_xx, g_y_x_0_0_x_x_0_xy, g_y_x_0_0_x_x_0_xz, g_y_x_0_0_x_x_0_yy, g_y_x_0_0_x_x_0_yz, g_y_x_0_0_x_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_x_0_xx[i] = -2.0 * g_xy_0_0_xx[i] * a_exp + 4.0 * g_xy_xx_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_0_xy[i] = -2.0 * g_xy_0_0_xy[i] * a_exp + 4.0 * g_xy_xx_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_0_xz[i] = -2.0 * g_xy_0_0_xz[i] * a_exp + 4.0 * g_xy_xx_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_0_yy[i] = -2.0 * g_xy_0_0_yy[i] * a_exp + 4.0 * g_xy_xx_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_0_yz[i] = -2.0 * g_xy_0_0_yz[i] * a_exp + 4.0 * g_xy_xx_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_0_zz[i] = -2.0 * g_xy_0_0_zz[i] * a_exp + 4.0 * g_xy_xx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (168-174)

    #pragma omp simd aligned(g_xy_xy_0_xx, g_xy_xy_0_xy, g_xy_xy_0_xz, g_xy_xy_0_yy, g_xy_xy_0_yz, g_xy_xy_0_zz, g_y_x_0_0_x_y_0_xx, g_y_x_0_0_x_y_0_xy, g_y_x_0_0_x_y_0_xz, g_y_x_0_0_x_y_0_yy, g_y_x_0_0_x_y_0_yz, g_y_x_0_0_x_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_y_0_xx[i] = 4.0 * g_xy_xy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_0_xy[i] = 4.0 * g_xy_xy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_0_xz[i] = 4.0 * g_xy_xy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_0_yy[i] = 4.0 * g_xy_xy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_0_yz[i] = 4.0 * g_xy_xy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_0_zz[i] = 4.0 * g_xy_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (174-180)

    #pragma omp simd aligned(g_xy_xz_0_xx, g_xy_xz_0_xy, g_xy_xz_0_xz, g_xy_xz_0_yy, g_xy_xz_0_yz, g_xy_xz_0_zz, g_y_x_0_0_x_z_0_xx, g_y_x_0_0_x_z_0_xy, g_y_x_0_0_x_z_0_xz, g_y_x_0_0_x_z_0_yy, g_y_x_0_0_x_z_0_yz, g_y_x_0_0_x_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_z_0_xx[i] = 4.0 * g_xy_xz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_0_xy[i] = 4.0 * g_xy_xz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_0_xz[i] = 4.0 * g_xy_xz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_0_yy[i] = 4.0 * g_xy_xz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_0_yz[i] = 4.0 * g_xy_xz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_0_zz[i] = 4.0 * g_xy_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (180-186)

    #pragma omp simd aligned(g_0_0_0_xx, g_0_0_0_xy, g_0_0_0_xz, g_0_0_0_yy, g_0_0_0_yz, g_0_0_0_zz, g_0_xx_0_xx, g_0_xx_0_xy, g_0_xx_0_xz, g_0_xx_0_yy, g_0_xx_0_yz, g_0_xx_0_zz, g_y_x_0_0_y_x_0_xx, g_y_x_0_0_y_x_0_xy, g_y_x_0_0_y_x_0_xz, g_y_x_0_0_y_x_0_yy, g_y_x_0_0_y_x_0_yz, g_y_x_0_0_y_x_0_zz, g_yy_0_0_xx, g_yy_0_0_xy, g_yy_0_0_xz, g_yy_0_0_yy, g_yy_0_0_yz, g_yy_0_0_zz, g_yy_xx_0_xx, g_yy_xx_0_xy, g_yy_xx_0_xz, g_yy_xx_0_yy, g_yy_xx_0_yz, g_yy_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_x_0_xx[i] = g_0_0_0_xx[i] - 2.0 * g_0_xx_0_xx[i] * b_exp - 2.0 * g_yy_0_0_xx[i] * a_exp + 4.0 * g_yy_xx_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_0_xy[i] = g_0_0_0_xy[i] - 2.0 * g_0_xx_0_xy[i] * b_exp - 2.0 * g_yy_0_0_xy[i] * a_exp + 4.0 * g_yy_xx_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_0_xz[i] = g_0_0_0_xz[i] - 2.0 * g_0_xx_0_xz[i] * b_exp - 2.0 * g_yy_0_0_xz[i] * a_exp + 4.0 * g_yy_xx_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_0_yy[i] = g_0_0_0_yy[i] - 2.0 * g_0_xx_0_yy[i] * b_exp - 2.0 * g_yy_0_0_yy[i] * a_exp + 4.0 * g_yy_xx_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_0_yz[i] = g_0_0_0_yz[i] - 2.0 * g_0_xx_0_yz[i] * b_exp - 2.0 * g_yy_0_0_yz[i] * a_exp + 4.0 * g_yy_xx_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_0_zz[i] = g_0_0_0_zz[i] - 2.0 * g_0_xx_0_zz[i] * b_exp - 2.0 * g_yy_0_0_zz[i] * a_exp + 4.0 * g_yy_xx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (186-192)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_y_x_0_0_y_y_0_xx, g_y_x_0_0_y_y_0_xy, g_y_x_0_0_y_y_0_xz, g_y_x_0_0_y_y_0_yy, g_y_x_0_0_y_y_0_yz, g_y_x_0_0_y_y_0_zz, g_yy_xy_0_xx, g_yy_xy_0_xy, g_yy_xy_0_xz, g_yy_xy_0_yy, g_yy_xy_0_yz, g_yy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_y_0_xx[i] = -2.0 * g_0_xy_0_xx[i] * b_exp + 4.0 * g_yy_xy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_0_xy[i] = -2.0 * g_0_xy_0_xy[i] * b_exp + 4.0 * g_yy_xy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_0_xz[i] = -2.0 * g_0_xy_0_xz[i] * b_exp + 4.0 * g_yy_xy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_0_yy[i] = -2.0 * g_0_xy_0_yy[i] * b_exp + 4.0 * g_yy_xy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_0_yz[i] = -2.0 * g_0_xy_0_yz[i] * b_exp + 4.0 * g_yy_xy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_0_zz[i] = -2.0 * g_0_xy_0_zz[i] * b_exp + 4.0 * g_yy_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (192-198)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_y_x_0_0_y_z_0_xx, g_y_x_0_0_y_z_0_xy, g_y_x_0_0_y_z_0_xz, g_y_x_0_0_y_z_0_yy, g_y_x_0_0_y_z_0_yz, g_y_x_0_0_y_z_0_zz, g_yy_xz_0_xx, g_yy_xz_0_xy, g_yy_xz_0_xz, g_yy_xz_0_yy, g_yy_xz_0_yz, g_yy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_z_0_xx[i] = -2.0 * g_0_xz_0_xx[i] * b_exp + 4.0 * g_yy_xz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_0_xy[i] = -2.0 * g_0_xz_0_xy[i] * b_exp + 4.0 * g_yy_xz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_0_xz[i] = -2.0 * g_0_xz_0_xz[i] * b_exp + 4.0 * g_yy_xz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_0_yy[i] = -2.0 * g_0_xz_0_yy[i] * b_exp + 4.0 * g_yy_xz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_0_yz[i] = -2.0 * g_0_xz_0_yz[i] * b_exp + 4.0 * g_yy_xz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_0_zz[i] = -2.0 * g_0_xz_0_zz[i] * b_exp + 4.0 * g_yy_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (198-204)

    #pragma omp simd aligned(g_y_x_0_0_z_x_0_xx, g_y_x_0_0_z_x_0_xy, g_y_x_0_0_z_x_0_xz, g_y_x_0_0_z_x_0_yy, g_y_x_0_0_z_x_0_yz, g_y_x_0_0_z_x_0_zz, g_yz_0_0_xx, g_yz_0_0_xy, g_yz_0_0_xz, g_yz_0_0_yy, g_yz_0_0_yz, g_yz_0_0_zz, g_yz_xx_0_xx, g_yz_xx_0_xy, g_yz_xx_0_xz, g_yz_xx_0_yy, g_yz_xx_0_yz, g_yz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_x_0_xx[i] = -2.0 * g_yz_0_0_xx[i] * a_exp + 4.0 * g_yz_xx_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_0_xy[i] = -2.0 * g_yz_0_0_xy[i] * a_exp + 4.0 * g_yz_xx_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_0_xz[i] = -2.0 * g_yz_0_0_xz[i] * a_exp + 4.0 * g_yz_xx_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_0_yy[i] = -2.0 * g_yz_0_0_yy[i] * a_exp + 4.0 * g_yz_xx_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_0_yz[i] = -2.0 * g_yz_0_0_yz[i] * a_exp + 4.0 * g_yz_xx_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_0_zz[i] = -2.0 * g_yz_0_0_zz[i] * a_exp + 4.0 * g_yz_xx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (204-210)

    #pragma omp simd aligned(g_y_x_0_0_z_y_0_xx, g_y_x_0_0_z_y_0_xy, g_y_x_0_0_z_y_0_xz, g_y_x_0_0_z_y_0_yy, g_y_x_0_0_z_y_0_yz, g_y_x_0_0_z_y_0_zz, g_yz_xy_0_xx, g_yz_xy_0_xy, g_yz_xy_0_xz, g_yz_xy_0_yy, g_yz_xy_0_yz, g_yz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_y_0_xx[i] = 4.0 * g_yz_xy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_0_xy[i] = 4.0 * g_yz_xy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_0_xz[i] = 4.0 * g_yz_xy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_0_yy[i] = 4.0 * g_yz_xy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_0_yz[i] = 4.0 * g_yz_xy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_0_zz[i] = 4.0 * g_yz_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (210-216)

    #pragma omp simd aligned(g_y_x_0_0_z_z_0_xx, g_y_x_0_0_z_z_0_xy, g_y_x_0_0_z_z_0_xz, g_y_x_0_0_z_z_0_yy, g_y_x_0_0_z_z_0_yz, g_y_x_0_0_z_z_0_zz, g_yz_xz_0_xx, g_yz_xz_0_xy, g_yz_xz_0_xz, g_yz_xz_0_yy, g_yz_xz_0_yz, g_yz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_z_0_xx[i] = 4.0 * g_yz_xz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_0_xy[i] = 4.0 * g_yz_xz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_0_xz[i] = 4.0 * g_yz_xz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_0_yy[i] = 4.0 * g_yz_xz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_0_yz[i] = 4.0 * g_yz_xz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_0_zz[i] = 4.0 * g_yz_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (216-222)

    #pragma omp simd aligned(g_xy_xy_0_xx, g_xy_xy_0_xy, g_xy_xy_0_xz, g_xy_xy_0_yy, g_xy_xy_0_yz, g_xy_xy_0_zz, g_y_y_0_0_x_x_0_xx, g_y_y_0_0_x_x_0_xy, g_y_y_0_0_x_x_0_xz, g_y_y_0_0_x_x_0_yy, g_y_y_0_0_x_x_0_yz, g_y_y_0_0_x_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_x_0_xx[i] = 4.0 * g_xy_xy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_0_xy[i] = 4.0 * g_xy_xy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_0_xz[i] = 4.0 * g_xy_xy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_0_yy[i] = 4.0 * g_xy_xy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_0_yz[i] = 4.0 * g_xy_xy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_0_zz[i] = 4.0 * g_xy_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (222-228)

    #pragma omp simd aligned(g_xy_0_0_xx, g_xy_0_0_xy, g_xy_0_0_xz, g_xy_0_0_yy, g_xy_0_0_yz, g_xy_0_0_zz, g_xy_yy_0_xx, g_xy_yy_0_xy, g_xy_yy_0_xz, g_xy_yy_0_yy, g_xy_yy_0_yz, g_xy_yy_0_zz, g_y_y_0_0_x_y_0_xx, g_y_y_0_0_x_y_0_xy, g_y_y_0_0_x_y_0_xz, g_y_y_0_0_x_y_0_yy, g_y_y_0_0_x_y_0_yz, g_y_y_0_0_x_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_y_0_xx[i] = -2.0 * g_xy_0_0_xx[i] * a_exp + 4.0 * g_xy_yy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_0_xy[i] = -2.0 * g_xy_0_0_xy[i] * a_exp + 4.0 * g_xy_yy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_0_xz[i] = -2.0 * g_xy_0_0_xz[i] * a_exp + 4.0 * g_xy_yy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_0_yy[i] = -2.0 * g_xy_0_0_yy[i] * a_exp + 4.0 * g_xy_yy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_0_yz[i] = -2.0 * g_xy_0_0_yz[i] * a_exp + 4.0 * g_xy_yy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_0_zz[i] = -2.0 * g_xy_0_0_zz[i] * a_exp + 4.0 * g_xy_yy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (228-234)

    #pragma omp simd aligned(g_xy_yz_0_xx, g_xy_yz_0_xy, g_xy_yz_0_xz, g_xy_yz_0_yy, g_xy_yz_0_yz, g_xy_yz_0_zz, g_y_y_0_0_x_z_0_xx, g_y_y_0_0_x_z_0_xy, g_y_y_0_0_x_z_0_xz, g_y_y_0_0_x_z_0_yy, g_y_y_0_0_x_z_0_yz, g_y_y_0_0_x_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_z_0_xx[i] = 4.0 * g_xy_yz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_0_xy[i] = 4.0 * g_xy_yz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_0_xz[i] = 4.0 * g_xy_yz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_0_yy[i] = 4.0 * g_xy_yz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_0_yz[i] = 4.0 * g_xy_yz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_0_zz[i] = 4.0 * g_xy_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (234-240)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_y_y_0_0_y_x_0_xx, g_y_y_0_0_y_x_0_xy, g_y_y_0_0_y_x_0_xz, g_y_y_0_0_y_x_0_yy, g_y_y_0_0_y_x_0_yz, g_y_y_0_0_y_x_0_zz, g_yy_xy_0_xx, g_yy_xy_0_xy, g_yy_xy_0_xz, g_yy_xy_0_yy, g_yy_xy_0_yz, g_yy_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_x_0_xx[i] = -2.0 * g_0_xy_0_xx[i] * b_exp + 4.0 * g_yy_xy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_0_xy[i] = -2.0 * g_0_xy_0_xy[i] * b_exp + 4.0 * g_yy_xy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_0_xz[i] = -2.0 * g_0_xy_0_xz[i] * b_exp + 4.0 * g_yy_xy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_0_yy[i] = -2.0 * g_0_xy_0_yy[i] * b_exp + 4.0 * g_yy_xy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_0_yz[i] = -2.0 * g_0_xy_0_yz[i] * b_exp + 4.0 * g_yy_xy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_0_zz[i] = -2.0 * g_0_xy_0_zz[i] * b_exp + 4.0 * g_yy_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (240-246)

    #pragma omp simd aligned(g_0_0_0_xx, g_0_0_0_xy, g_0_0_0_xz, g_0_0_0_yy, g_0_0_0_yz, g_0_0_0_zz, g_0_yy_0_xx, g_0_yy_0_xy, g_0_yy_0_xz, g_0_yy_0_yy, g_0_yy_0_yz, g_0_yy_0_zz, g_y_y_0_0_y_y_0_xx, g_y_y_0_0_y_y_0_xy, g_y_y_0_0_y_y_0_xz, g_y_y_0_0_y_y_0_yy, g_y_y_0_0_y_y_0_yz, g_y_y_0_0_y_y_0_zz, g_yy_0_0_xx, g_yy_0_0_xy, g_yy_0_0_xz, g_yy_0_0_yy, g_yy_0_0_yz, g_yy_0_0_zz, g_yy_yy_0_xx, g_yy_yy_0_xy, g_yy_yy_0_xz, g_yy_yy_0_yy, g_yy_yy_0_yz, g_yy_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_y_0_xx[i] = g_0_0_0_xx[i] - 2.0 * g_0_yy_0_xx[i] * b_exp - 2.0 * g_yy_0_0_xx[i] * a_exp + 4.0 * g_yy_yy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_0_xy[i] = g_0_0_0_xy[i] - 2.0 * g_0_yy_0_xy[i] * b_exp - 2.0 * g_yy_0_0_xy[i] * a_exp + 4.0 * g_yy_yy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_0_xz[i] = g_0_0_0_xz[i] - 2.0 * g_0_yy_0_xz[i] * b_exp - 2.0 * g_yy_0_0_xz[i] * a_exp + 4.0 * g_yy_yy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_0_yy[i] = g_0_0_0_yy[i] - 2.0 * g_0_yy_0_yy[i] * b_exp - 2.0 * g_yy_0_0_yy[i] * a_exp + 4.0 * g_yy_yy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_0_yz[i] = g_0_0_0_yz[i] - 2.0 * g_0_yy_0_yz[i] * b_exp - 2.0 * g_yy_0_0_yz[i] * a_exp + 4.0 * g_yy_yy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_0_zz[i] = g_0_0_0_zz[i] - 2.0 * g_0_yy_0_zz[i] * b_exp - 2.0 * g_yy_0_0_zz[i] * a_exp + 4.0 * g_yy_yy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (246-252)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_y_y_0_0_y_z_0_xx, g_y_y_0_0_y_z_0_xy, g_y_y_0_0_y_z_0_xz, g_y_y_0_0_y_z_0_yy, g_y_y_0_0_y_z_0_yz, g_y_y_0_0_y_z_0_zz, g_yy_yz_0_xx, g_yy_yz_0_xy, g_yy_yz_0_xz, g_yy_yz_0_yy, g_yy_yz_0_yz, g_yy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_z_0_xx[i] = -2.0 * g_0_yz_0_xx[i] * b_exp + 4.0 * g_yy_yz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_0_xy[i] = -2.0 * g_0_yz_0_xy[i] * b_exp + 4.0 * g_yy_yz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_0_xz[i] = -2.0 * g_0_yz_0_xz[i] * b_exp + 4.0 * g_yy_yz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_0_yy[i] = -2.0 * g_0_yz_0_yy[i] * b_exp + 4.0 * g_yy_yz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_0_yz[i] = -2.0 * g_0_yz_0_yz[i] * b_exp + 4.0 * g_yy_yz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_0_zz[i] = -2.0 * g_0_yz_0_zz[i] * b_exp + 4.0 * g_yy_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (252-258)

    #pragma omp simd aligned(g_y_y_0_0_z_x_0_xx, g_y_y_0_0_z_x_0_xy, g_y_y_0_0_z_x_0_xz, g_y_y_0_0_z_x_0_yy, g_y_y_0_0_z_x_0_yz, g_y_y_0_0_z_x_0_zz, g_yz_xy_0_xx, g_yz_xy_0_xy, g_yz_xy_0_xz, g_yz_xy_0_yy, g_yz_xy_0_yz, g_yz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_x_0_xx[i] = 4.0 * g_yz_xy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_0_xy[i] = 4.0 * g_yz_xy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_0_xz[i] = 4.0 * g_yz_xy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_0_yy[i] = 4.0 * g_yz_xy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_0_yz[i] = 4.0 * g_yz_xy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_0_zz[i] = 4.0 * g_yz_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (258-264)

    #pragma omp simd aligned(g_y_y_0_0_z_y_0_xx, g_y_y_0_0_z_y_0_xy, g_y_y_0_0_z_y_0_xz, g_y_y_0_0_z_y_0_yy, g_y_y_0_0_z_y_0_yz, g_y_y_0_0_z_y_0_zz, g_yz_0_0_xx, g_yz_0_0_xy, g_yz_0_0_xz, g_yz_0_0_yy, g_yz_0_0_yz, g_yz_0_0_zz, g_yz_yy_0_xx, g_yz_yy_0_xy, g_yz_yy_0_xz, g_yz_yy_0_yy, g_yz_yy_0_yz, g_yz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_y_0_xx[i] = -2.0 * g_yz_0_0_xx[i] * a_exp + 4.0 * g_yz_yy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_0_xy[i] = -2.0 * g_yz_0_0_xy[i] * a_exp + 4.0 * g_yz_yy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_0_xz[i] = -2.0 * g_yz_0_0_xz[i] * a_exp + 4.0 * g_yz_yy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_0_yy[i] = -2.0 * g_yz_0_0_yy[i] * a_exp + 4.0 * g_yz_yy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_0_yz[i] = -2.0 * g_yz_0_0_yz[i] * a_exp + 4.0 * g_yz_yy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_0_zz[i] = -2.0 * g_yz_0_0_zz[i] * a_exp + 4.0 * g_yz_yy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (264-270)

    #pragma omp simd aligned(g_y_y_0_0_z_z_0_xx, g_y_y_0_0_z_z_0_xy, g_y_y_0_0_z_z_0_xz, g_y_y_0_0_z_z_0_yy, g_y_y_0_0_z_z_0_yz, g_y_y_0_0_z_z_0_zz, g_yz_yz_0_xx, g_yz_yz_0_xy, g_yz_yz_0_xz, g_yz_yz_0_yy, g_yz_yz_0_yz, g_yz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_z_0_xx[i] = 4.0 * g_yz_yz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_0_xy[i] = 4.0 * g_yz_yz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_0_xz[i] = 4.0 * g_yz_yz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_0_yy[i] = 4.0 * g_yz_yz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_0_yz[i] = 4.0 * g_yz_yz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_0_zz[i] = 4.0 * g_yz_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (270-276)

    #pragma omp simd aligned(g_xy_xz_0_xx, g_xy_xz_0_xy, g_xy_xz_0_xz, g_xy_xz_0_yy, g_xy_xz_0_yz, g_xy_xz_0_zz, g_y_z_0_0_x_x_0_xx, g_y_z_0_0_x_x_0_xy, g_y_z_0_0_x_x_0_xz, g_y_z_0_0_x_x_0_yy, g_y_z_0_0_x_x_0_yz, g_y_z_0_0_x_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_x_0_xx[i] = 4.0 * g_xy_xz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_0_xy[i] = 4.0 * g_xy_xz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_0_xz[i] = 4.0 * g_xy_xz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_0_yy[i] = 4.0 * g_xy_xz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_0_yz[i] = 4.0 * g_xy_xz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_0_zz[i] = 4.0 * g_xy_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (276-282)

    #pragma omp simd aligned(g_xy_yz_0_xx, g_xy_yz_0_xy, g_xy_yz_0_xz, g_xy_yz_0_yy, g_xy_yz_0_yz, g_xy_yz_0_zz, g_y_z_0_0_x_y_0_xx, g_y_z_0_0_x_y_0_xy, g_y_z_0_0_x_y_0_xz, g_y_z_0_0_x_y_0_yy, g_y_z_0_0_x_y_0_yz, g_y_z_0_0_x_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_y_0_xx[i] = 4.0 * g_xy_yz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_0_xy[i] = 4.0 * g_xy_yz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_0_xz[i] = 4.0 * g_xy_yz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_0_yy[i] = 4.0 * g_xy_yz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_0_yz[i] = 4.0 * g_xy_yz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_0_zz[i] = 4.0 * g_xy_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (282-288)

    #pragma omp simd aligned(g_xy_0_0_xx, g_xy_0_0_xy, g_xy_0_0_xz, g_xy_0_0_yy, g_xy_0_0_yz, g_xy_0_0_zz, g_xy_zz_0_xx, g_xy_zz_0_xy, g_xy_zz_0_xz, g_xy_zz_0_yy, g_xy_zz_0_yz, g_xy_zz_0_zz, g_y_z_0_0_x_z_0_xx, g_y_z_0_0_x_z_0_xy, g_y_z_0_0_x_z_0_xz, g_y_z_0_0_x_z_0_yy, g_y_z_0_0_x_z_0_yz, g_y_z_0_0_x_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_z_0_xx[i] = -2.0 * g_xy_0_0_xx[i] * a_exp + 4.0 * g_xy_zz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_0_xy[i] = -2.0 * g_xy_0_0_xy[i] * a_exp + 4.0 * g_xy_zz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_0_xz[i] = -2.0 * g_xy_0_0_xz[i] * a_exp + 4.0 * g_xy_zz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_0_yy[i] = -2.0 * g_xy_0_0_yy[i] * a_exp + 4.0 * g_xy_zz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_0_yz[i] = -2.0 * g_xy_0_0_yz[i] * a_exp + 4.0 * g_xy_zz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_0_zz[i] = -2.0 * g_xy_0_0_zz[i] * a_exp + 4.0 * g_xy_zz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (288-294)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_y_z_0_0_y_x_0_xx, g_y_z_0_0_y_x_0_xy, g_y_z_0_0_y_x_0_xz, g_y_z_0_0_y_x_0_yy, g_y_z_0_0_y_x_0_yz, g_y_z_0_0_y_x_0_zz, g_yy_xz_0_xx, g_yy_xz_0_xy, g_yy_xz_0_xz, g_yy_xz_0_yy, g_yy_xz_0_yz, g_yy_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_x_0_xx[i] = -2.0 * g_0_xz_0_xx[i] * b_exp + 4.0 * g_yy_xz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_0_xy[i] = -2.0 * g_0_xz_0_xy[i] * b_exp + 4.0 * g_yy_xz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_0_xz[i] = -2.0 * g_0_xz_0_xz[i] * b_exp + 4.0 * g_yy_xz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_0_yy[i] = -2.0 * g_0_xz_0_yy[i] * b_exp + 4.0 * g_yy_xz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_0_yz[i] = -2.0 * g_0_xz_0_yz[i] * b_exp + 4.0 * g_yy_xz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_0_zz[i] = -2.0 * g_0_xz_0_zz[i] * b_exp + 4.0 * g_yy_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (294-300)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_y_z_0_0_y_y_0_xx, g_y_z_0_0_y_y_0_xy, g_y_z_0_0_y_y_0_xz, g_y_z_0_0_y_y_0_yy, g_y_z_0_0_y_y_0_yz, g_y_z_0_0_y_y_0_zz, g_yy_yz_0_xx, g_yy_yz_0_xy, g_yy_yz_0_xz, g_yy_yz_0_yy, g_yy_yz_0_yz, g_yy_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_y_0_xx[i] = -2.0 * g_0_yz_0_xx[i] * b_exp + 4.0 * g_yy_yz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_0_xy[i] = -2.0 * g_0_yz_0_xy[i] * b_exp + 4.0 * g_yy_yz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_0_xz[i] = -2.0 * g_0_yz_0_xz[i] * b_exp + 4.0 * g_yy_yz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_0_yy[i] = -2.0 * g_0_yz_0_yy[i] * b_exp + 4.0 * g_yy_yz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_0_yz[i] = -2.0 * g_0_yz_0_yz[i] * b_exp + 4.0 * g_yy_yz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_0_zz[i] = -2.0 * g_0_yz_0_zz[i] * b_exp + 4.0 * g_yy_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (300-306)

    #pragma omp simd aligned(g_0_0_0_xx, g_0_0_0_xy, g_0_0_0_xz, g_0_0_0_yy, g_0_0_0_yz, g_0_0_0_zz, g_0_zz_0_xx, g_0_zz_0_xy, g_0_zz_0_xz, g_0_zz_0_yy, g_0_zz_0_yz, g_0_zz_0_zz, g_y_z_0_0_y_z_0_xx, g_y_z_0_0_y_z_0_xy, g_y_z_0_0_y_z_0_xz, g_y_z_0_0_y_z_0_yy, g_y_z_0_0_y_z_0_yz, g_y_z_0_0_y_z_0_zz, g_yy_0_0_xx, g_yy_0_0_xy, g_yy_0_0_xz, g_yy_0_0_yy, g_yy_0_0_yz, g_yy_0_0_zz, g_yy_zz_0_xx, g_yy_zz_0_xy, g_yy_zz_0_xz, g_yy_zz_0_yy, g_yy_zz_0_yz, g_yy_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_z_0_xx[i] = g_0_0_0_xx[i] - 2.0 * g_0_zz_0_xx[i] * b_exp - 2.0 * g_yy_0_0_xx[i] * a_exp + 4.0 * g_yy_zz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_0_xy[i] = g_0_0_0_xy[i] - 2.0 * g_0_zz_0_xy[i] * b_exp - 2.0 * g_yy_0_0_xy[i] * a_exp + 4.0 * g_yy_zz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_0_xz[i] = g_0_0_0_xz[i] - 2.0 * g_0_zz_0_xz[i] * b_exp - 2.0 * g_yy_0_0_xz[i] * a_exp + 4.0 * g_yy_zz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_0_yy[i] = g_0_0_0_yy[i] - 2.0 * g_0_zz_0_yy[i] * b_exp - 2.0 * g_yy_0_0_yy[i] * a_exp + 4.0 * g_yy_zz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_0_yz[i] = g_0_0_0_yz[i] - 2.0 * g_0_zz_0_yz[i] * b_exp - 2.0 * g_yy_0_0_yz[i] * a_exp + 4.0 * g_yy_zz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_0_zz[i] = g_0_0_0_zz[i] - 2.0 * g_0_zz_0_zz[i] * b_exp - 2.0 * g_yy_0_0_zz[i] * a_exp + 4.0 * g_yy_zz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (306-312)

    #pragma omp simd aligned(g_y_z_0_0_z_x_0_xx, g_y_z_0_0_z_x_0_xy, g_y_z_0_0_z_x_0_xz, g_y_z_0_0_z_x_0_yy, g_y_z_0_0_z_x_0_yz, g_y_z_0_0_z_x_0_zz, g_yz_xz_0_xx, g_yz_xz_0_xy, g_yz_xz_0_xz, g_yz_xz_0_yy, g_yz_xz_0_yz, g_yz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_x_0_xx[i] = 4.0 * g_yz_xz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_0_xy[i] = 4.0 * g_yz_xz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_0_xz[i] = 4.0 * g_yz_xz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_0_yy[i] = 4.0 * g_yz_xz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_0_yz[i] = 4.0 * g_yz_xz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_0_zz[i] = 4.0 * g_yz_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (312-318)

    #pragma omp simd aligned(g_y_z_0_0_z_y_0_xx, g_y_z_0_0_z_y_0_xy, g_y_z_0_0_z_y_0_xz, g_y_z_0_0_z_y_0_yy, g_y_z_0_0_z_y_0_yz, g_y_z_0_0_z_y_0_zz, g_yz_yz_0_xx, g_yz_yz_0_xy, g_yz_yz_0_xz, g_yz_yz_0_yy, g_yz_yz_0_yz, g_yz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_y_0_xx[i] = 4.0 * g_yz_yz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_0_xy[i] = 4.0 * g_yz_yz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_0_xz[i] = 4.0 * g_yz_yz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_0_yy[i] = 4.0 * g_yz_yz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_0_yz[i] = 4.0 * g_yz_yz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_0_zz[i] = 4.0 * g_yz_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (318-324)

    #pragma omp simd aligned(g_y_z_0_0_z_z_0_xx, g_y_z_0_0_z_z_0_xy, g_y_z_0_0_z_z_0_xz, g_y_z_0_0_z_z_0_yy, g_y_z_0_0_z_z_0_yz, g_y_z_0_0_z_z_0_zz, g_yz_0_0_xx, g_yz_0_0_xy, g_yz_0_0_xz, g_yz_0_0_yy, g_yz_0_0_yz, g_yz_0_0_zz, g_yz_zz_0_xx, g_yz_zz_0_xy, g_yz_zz_0_xz, g_yz_zz_0_yy, g_yz_zz_0_yz, g_yz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_z_0_xx[i] = -2.0 * g_yz_0_0_xx[i] * a_exp + 4.0 * g_yz_zz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_0_xy[i] = -2.0 * g_yz_0_0_xy[i] * a_exp + 4.0 * g_yz_zz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_0_xz[i] = -2.0 * g_yz_0_0_xz[i] * a_exp + 4.0 * g_yz_zz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_0_yy[i] = -2.0 * g_yz_0_0_yy[i] * a_exp + 4.0 * g_yz_zz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_0_yz[i] = -2.0 * g_yz_0_0_yz[i] * a_exp + 4.0 * g_yz_zz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_0_zz[i] = -2.0 * g_yz_0_0_zz[i] * a_exp + 4.0 * g_yz_zz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (324-330)

    #pragma omp simd aligned(g_xz_0_0_xx, g_xz_0_0_xy, g_xz_0_0_xz, g_xz_0_0_yy, g_xz_0_0_yz, g_xz_0_0_zz, g_xz_xx_0_xx, g_xz_xx_0_xy, g_xz_xx_0_xz, g_xz_xx_0_yy, g_xz_xx_0_yz, g_xz_xx_0_zz, g_z_x_0_0_x_x_0_xx, g_z_x_0_0_x_x_0_xy, g_z_x_0_0_x_x_0_xz, g_z_x_0_0_x_x_0_yy, g_z_x_0_0_x_x_0_yz, g_z_x_0_0_x_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_x_0_xx[i] = -2.0 * g_xz_0_0_xx[i] * a_exp + 4.0 * g_xz_xx_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_0_xy[i] = -2.0 * g_xz_0_0_xy[i] * a_exp + 4.0 * g_xz_xx_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_0_xz[i] = -2.0 * g_xz_0_0_xz[i] * a_exp + 4.0 * g_xz_xx_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_0_yy[i] = -2.0 * g_xz_0_0_yy[i] * a_exp + 4.0 * g_xz_xx_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_0_yz[i] = -2.0 * g_xz_0_0_yz[i] * a_exp + 4.0 * g_xz_xx_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_0_zz[i] = -2.0 * g_xz_0_0_zz[i] * a_exp + 4.0 * g_xz_xx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (330-336)

    #pragma omp simd aligned(g_xz_xy_0_xx, g_xz_xy_0_xy, g_xz_xy_0_xz, g_xz_xy_0_yy, g_xz_xy_0_yz, g_xz_xy_0_zz, g_z_x_0_0_x_y_0_xx, g_z_x_0_0_x_y_0_xy, g_z_x_0_0_x_y_0_xz, g_z_x_0_0_x_y_0_yy, g_z_x_0_0_x_y_0_yz, g_z_x_0_0_x_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_y_0_xx[i] = 4.0 * g_xz_xy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_0_xy[i] = 4.0 * g_xz_xy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_0_xz[i] = 4.0 * g_xz_xy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_0_yy[i] = 4.0 * g_xz_xy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_0_yz[i] = 4.0 * g_xz_xy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_0_zz[i] = 4.0 * g_xz_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (336-342)

    #pragma omp simd aligned(g_xz_xz_0_xx, g_xz_xz_0_xy, g_xz_xz_0_xz, g_xz_xz_0_yy, g_xz_xz_0_yz, g_xz_xz_0_zz, g_z_x_0_0_x_z_0_xx, g_z_x_0_0_x_z_0_xy, g_z_x_0_0_x_z_0_xz, g_z_x_0_0_x_z_0_yy, g_z_x_0_0_x_z_0_yz, g_z_x_0_0_x_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_z_0_xx[i] = 4.0 * g_xz_xz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_0_xy[i] = 4.0 * g_xz_xz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_0_xz[i] = 4.0 * g_xz_xz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_0_yy[i] = 4.0 * g_xz_xz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_0_yz[i] = 4.0 * g_xz_xz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_0_zz[i] = 4.0 * g_xz_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (342-348)

    #pragma omp simd aligned(g_yz_0_0_xx, g_yz_0_0_xy, g_yz_0_0_xz, g_yz_0_0_yy, g_yz_0_0_yz, g_yz_0_0_zz, g_yz_xx_0_xx, g_yz_xx_0_xy, g_yz_xx_0_xz, g_yz_xx_0_yy, g_yz_xx_0_yz, g_yz_xx_0_zz, g_z_x_0_0_y_x_0_xx, g_z_x_0_0_y_x_0_xy, g_z_x_0_0_y_x_0_xz, g_z_x_0_0_y_x_0_yy, g_z_x_0_0_y_x_0_yz, g_z_x_0_0_y_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_x_0_xx[i] = -2.0 * g_yz_0_0_xx[i] * a_exp + 4.0 * g_yz_xx_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_0_xy[i] = -2.0 * g_yz_0_0_xy[i] * a_exp + 4.0 * g_yz_xx_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_0_xz[i] = -2.0 * g_yz_0_0_xz[i] * a_exp + 4.0 * g_yz_xx_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_0_yy[i] = -2.0 * g_yz_0_0_yy[i] * a_exp + 4.0 * g_yz_xx_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_0_yz[i] = -2.0 * g_yz_0_0_yz[i] * a_exp + 4.0 * g_yz_xx_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_0_zz[i] = -2.0 * g_yz_0_0_zz[i] * a_exp + 4.0 * g_yz_xx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (348-354)

    #pragma omp simd aligned(g_yz_xy_0_xx, g_yz_xy_0_xy, g_yz_xy_0_xz, g_yz_xy_0_yy, g_yz_xy_0_yz, g_yz_xy_0_zz, g_z_x_0_0_y_y_0_xx, g_z_x_0_0_y_y_0_xy, g_z_x_0_0_y_y_0_xz, g_z_x_0_0_y_y_0_yy, g_z_x_0_0_y_y_0_yz, g_z_x_0_0_y_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_y_0_xx[i] = 4.0 * g_yz_xy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_0_xy[i] = 4.0 * g_yz_xy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_0_xz[i] = 4.0 * g_yz_xy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_0_yy[i] = 4.0 * g_yz_xy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_0_yz[i] = 4.0 * g_yz_xy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_0_zz[i] = 4.0 * g_yz_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (354-360)

    #pragma omp simd aligned(g_yz_xz_0_xx, g_yz_xz_0_xy, g_yz_xz_0_xz, g_yz_xz_0_yy, g_yz_xz_0_yz, g_yz_xz_0_zz, g_z_x_0_0_y_z_0_xx, g_z_x_0_0_y_z_0_xy, g_z_x_0_0_y_z_0_xz, g_z_x_0_0_y_z_0_yy, g_z_x_0_0_y_z_0_yz, g_z_x_0_0_y_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_z_0_xx[i] = 4.0 * g_yz_xz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_0_xy[i] = 4.0 * g_yz_xz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_0_xz[i] = 4.0 * g_yz_xz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_0_yy[i] = 4.0 * g_yz_xz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_0_yz[i] = 4.0 * g_yz_xz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_0_zz[i] = 4.0 * g_yz_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (360-366)

    #pragma omp simd aligned(g_0_0_0_xx, g_0_0_0_xy, g_0_0_0_xz, g_0_0_0_yy, g_0_0_0_yz, g_0_0_0_zz, g_0_xx_0_xx, g_0_xx_0_xy, g_0_xx_0_xz, g_0_xx_0_yy, g_0_xx_0_yz, g_0_xx_0_zz, g_z_x_0_0_z_x_0_xx, g_z_x_0_0_z_x_0_xy, g_z_x_0_0_z_x_0_xz, g_z_x_0_0_z_x_0_yy, g_z_x_0_0_z_x_0_yz, g_z_x_0_0_z_x_0_zz, g_zz_0_0_xx, g_zz_0_0_xy, g_zz_0_0_xz, g_zz_0_0_yy, g_zz_0_0_yz, g_zz_0_0_zz, g_zz_xx_0_xx, g_zz_xx_0_xy, g_zz_xx_0_xz, g_zz_xx_0_yy, g_zz_xx_0_yz, g_zz_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_x_0_xx[i] = g_0_0_0_xx[i] - 2.0 * g_0_xx_0_xx[i] * b_exp - 2.0 * g_zz_0_0_xx[i] * a_exp + 4.0 * g_zz_xx_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_0_xy[i] = g_0_0_0_xy[i] - 2.0 * g_0_xx_0_xy[i] * b_exp - 2.0 * g_zz_0_0_xy[i] * a_exp + 4.0 * g_zz_xx_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_0_xz[i] = g_0_0_0_xz[i] - 2.0 * g_0_xx_0_xz[i] * b_exp - 2.0 * g_zz_0_0_xz[i] * a_exp + 4.0 * g_zz_xx_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_0_yy[i] = g_0_0_0_yy[i] - 2.0 * g_0_xx_0_yy[i] * b_exp - 2.0 * g_zz_0_0_yy[i] * a_exp + 4.0 * g_zz_xx_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_0_yz[i] = g_0_0_0_yz[i] - 2.0 * g_0_xx_0_yz[i] * b_exp - 2.0 * g_zz_0_0_yz[i] * a_exp + 4.0 * g_zz_xx_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_0_zz[i] = g_0_0_0_zz[i] - 2.0 * g_0_xx_0_zz[i] * b_exp - 2.0 * g_zz_0_0_zz[i] * a_exp + 4.0 * g_zz_xx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (366-372)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_z_x_0_0_z_y_0_xx, g_z_x_0_0_z_y_0_xy, g_z_x_0_0_z_y_0_xz, g_z_x_0_0_z_y_0_yy, g_z_x_0_0_z_y_0_yz, g_z_x_0_0_z_y_0_zz, g_zz_xy_0_xx, g_zz_xy_0_xy, g_zz_xy_0_xz, g_zz_xy_0_yy, g_zz_xy_0_yz, g_zz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_y_0_xx[i] = -2.0 * g_0_xy_0_xx[i] * b_exp + 4.0 * g_zz_xy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_0_xy[i] = -2.0 * g_0_xy_0_xy[i] * b_exp + 4.0 * g_zz_xy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_0_xz[i] = -2.0 * g_0_xy_0_xz[i] * b_exp + 4.0 * g_zz_xy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_0_yy[i] = -2.0 * g_0_xy_0_yy[i] * b_exp + 4.0 * g_zz_xy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_0_yz[i] = -2.0 * g_0_xy_0_yz[i] * b_exp + 4.0 * g_zz_xy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_0_zz[i] = -2.0 * g_0_xy_0_zz[i] * b_exp + 4.0 * g_zz_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (372-378)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_z_x_0_0_z_z_0_xx, g_z_x_0_0_z_z_0_xy, g_z_x_0_0_z_z_0_xz, g_z_x_0_0_z_z_0_yy, g_z_x_0_0_z_z_0_yz, g_z_x_0_0_z_z_0_zz, g_zz_xz_0_xx, g_zz_xz_0_xy, g_zz_xz_0_xz, g_zz_xz_0_yy, g_zz_xz_0_yz, g_zz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_z_0_xx[i] = -2.0 * g_0_xz_0_xx[i] * b_exp + 4.0 * g_zz_xz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_0_xy[i] = -2.0 * g_0_xz_0_xy[i] * b_exp + 4.0 * g_zz_xz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_0_xz[i] = -2.0 * g_0_xz_0_xz[i] * b_exp + 4.0 * g_zz_xz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_0_yy[i] = -2.0 * g_0_xz_0_yy[i] * b_exp + 4.0 * g_zz_xz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_0_yz[i] = -2.0 * g_0_xz_0_yz[i] * b_exp + 4.0 * g_zz_xz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_0_zz[i] = -2.0 * g_0_xz_0_zz[i] * b_exp + 4.0 * g_zz_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (378-384)

    #pragma omp simd aligned(g_xz_xy_0_xx, g_xz_xy_0_xy, g_xz_xy_0_xz, g_xz_xy_0_yy, g_xz_xy_0_yz, g_xz_xy_0_zz, g_z_y_0_0_x_x_0_xx, g_z_y_0_0_x_x_0_xy, g_z_y_0_0_x_x_0_xz, g_z_y_0_0_x_x_0_yy, g_z_y_0_0_x_x_0_yz, g_z_y_0_0_x_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_x_0_xx[i] = 4.0 * g_xz_xy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_0_xy[i] = 4.0 * g_xz_xy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_0_xz[i] = 4.0 * g_xz_xy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_0_yy[i] = 4.0 * g_xz_xy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_0_yz[i] = 4.0 * g_xz_xy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_0_zz[i] = 4.0 * g_xz_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (384-390)

    #pragma omp simd aligned(g_xz_0_0_xx, g_xz_0_0_xy, g_xz_0_0_xz, g_xz_0_0_yy, g_xz_0_0_yz, g_xz_0_0_zz, g_xz_yy_0_xx, g_xz_yy_0_xy, g_xz_yy_0_xz, g_xz_yy_0_yy, g_xz_yy_0_yz, g_xz_yy_0_zz, g_z_y_0_0_x_y_0_xx, g_z_y_0_0_x_y_0_xy, g_z_y_0_0_x_y_0_xz, g_z_y_0_0_x_y_0_yy, g_z_y_0_0_x_y_0_yz, g_z_y_0_0_x_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_y_0_xx[i] = -2.0 * g_xz_0_0_xx[i] * a_exp + 4.0 * g_xz_yy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_0_xy[i] = -2.0 * g_xz_0_0_xy[i] * a_exp + 4.0 * g_xz_yy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_0_xz[i] = -2.0 * g_xz_0_0_xz[i] * a_exp + 4.0 * g_xz_yy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_0_yy[i] = -2.0 * g_xz_0_0_yy[i] * a_exp + 4.0 * g_xz_yy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_0_yz[i] = -2.0 * g_xz_0_0_yz[i] * a_exp + 4.0 * g_xz_yy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_0_zz[i] = -2.0 * g_xz_0_0_zz[i] * a_exp + 4.0 * g_xz_yy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (390-396)

    #pragma omp simd aligned(g_xz_yz_0_xx, g_xz_yz_0_xy, g_xz_yz_0_xz, g_xz_yz_0_yy, g_xz_yz_0_yz, g_xz_yz_0_zz, g_z_y_0_0_x_z_0_xx, g_z_y_0_0_x_z_0_xy, g_z_y_0_0_x_z_0_xz, g_z_y_0_0_x_z_0_yy, g_z_y_0_0_x_z_0_yz, g_z_y_0_0_x_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_z_0_xx[i] = 4.0 * g_xz_yz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_0_xy[i] = 4.0 * g_xz_yz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_0_xz[i] = 4.0 * g_xz_yz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_0_yy[i] = 4.0 * g_xz_yz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_0_yz[i] = 4.0 * g_xz_yz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_0_zz[i] = 4.0 * g_xz_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (396-402)

    #pragma omp simd aligned(g_yz_xy_0_xx, g_yz_xy_0_xy, g_yz_xy_0_xz, g_yz_xy_0_yy, g_yz_xy_0_yz, g_yz_xy_0_zz, g_z_y_0_0_y_x_0_xx, g_z_y_0_0_y_x_0_xy, g_z_y_0_0_y_x_0_xz, g_z_y_0_0_y_x_0_yy, g_z_y_0_0_y_x_0_yz, g_z_y_0_0_y_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_x_0_xx[i] = 4.0 * g_yz_xy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_0_xy[i] = 4.0 * g_yz_xy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_0_xz[i] = 4.0 * g_yz_xy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_0_yy[i] = 4.0 * g_yz_xy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_0_yz[i] = 4.0 * g_yz_xy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_0_zz[i] = 4.0 * g_yz_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (402-408)

    #pragma omp simd aligned(g_yz_0_0_xx, g_yz_0_0_xy, g_yz_0_0_xz, g_yz_0_0_yy, g_yz_0_0_yz, g_yz_0_0_zz, g_yz_yy_0_xx, g_yz_yy_0_xy, g_yz_yy_0_xz, g_yz_yy_0_yy, g_yz_yy_0_yz, g_yz_yy_0_zz, g_z_y_0_0_y_y_0_xx, g_z_y_0_0_y_y_0_xy, g_z_y_0_0_y_y_0_xz, g_z_y_0_0_y_y_0_yy, g_z_y_0_0_y_y_0_yz, g_z_y_0_0_y_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_y_0_xx[i] = -2.0 * g_yz_0_0_xx[i] * a_exp + 4.0 * g_yz_yy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_0_xy[i] = -2.0 * g_yz_0_0_xy[i] * a_exp + 4.0 * g_yz_yy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_0_xz[i] = -2.0 * g_yz_0_0_xz[i] * a_exp + 4.0 * g_yz_yy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_0_yy[i] = -2.0 * g_yz_0_0_yy[i] * a_exp + 4.0 * g_yz_yy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_0_yz[i] = -2.0 * g_yz_0_0_yz[i] * a_exp + 4.0 * g_yz_yy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_0_zz[i] = -2.0 * g_yz_0_0_zz[i] * a_exp + 4.0 * g_yz_yy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (408-414)

    #pragma omp simd aligned(g_yz_yz_0_xx, g_yz_yz_0_xy, g_yz_yz_0_xz, g_yz_yz_0_yy, g_yz_yz_0_yz, g_yz_yz_0_zz, g_z_y_0_0_y_z_0_xx, g_z_y_0_0_y_z_0_xy, g_z_y_0_0_y_z_0_xz, g_z_y_0_0_y_z_0_yy, g_z_y_0_0_y_z_0_yz, g_z_y_0_0_y_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_z_0_xx[i] = 4.0 * g_yz_yz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_0_xy[i] = 4.0 * g_yz_yz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_0_xz[i] = 4.0 * g_yz_yz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_0_yy[i] = 4.0 * g_yz_yz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_0_yz[i] = 4.0 * g_yz_yz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_0_zz[i] = 4.0 * g_yz_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (414-420)

    #pragma omp simd aligned(g_0_xy_0_xx, g_0_xy_0_xy, g_0_xy_0_xz, g_0_xy_0_yy, g_0_xy_0_yz, g_0_xy_0_zz, g_z_y_0_0_z_x_0_xx, g_z_y_0_0_z_x_0_xy, g_z_y_0_0_z_x_0_xz, g_z_y_0_0_z_x_0_yy, g_z_y_0_0_z_x_0_yz, g_z_y_0_0_z_x_0_zz, g_zz_xy_0_xx, g_zz_xy_0_xy, g_zz_xy_0_xz, g_zz_xy_0_yy, g_zz_xy_0_yz, g_zz_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_x_0_xx[i] = -2.0 * g_0_xy_0_xx[i] * b_exp + 4.0 * g_zz_xy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_0_xy[i] = -2.0 * g_0_xy_0_xy[i] * b_exp + 4.0 * g_zz_xy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_0_xz[i] = -2.0 * g_0_xy_0_xz[i] * b_exp + 4.0 * g_zz_xy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_0_yy[i] = -2.0 * g_0_xy_0_yy[i] * b_exp + 4.0 * g_zz_xy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_0_yz[i] = -2.0 * g_0_xy_0_yz[i] * b_exp + 4.0 * g_zz_xy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_0_zz[i] = -2.0 * g_0_xy_0_zz[i] * b_exp + 4.0 * g_zz_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (420-426)

    #pragma omp simd aligned(g_0_0_0_xx, g_0_0_0_xy, g_0_0_0_xz, g_0_0_0_yy, g_0_0_0_yz, g_0_0_0_zz, g_0_yy_0_xx, g_0_yy_0_xy, g_0_yy_0_xz, g_0_yy_0_yy, g_0_yy_0_yz, g_0_yy_0_zz, g_z_y_0_0_z_y_0_xx, g_z_y_0_0_z_y_0_xy, g_z_y_0_0_z_y_0_xz, g_z_y_0_0_z_y_0_yy, g_z_y_0_0_z_y_0_yz, g_z_y_0_0_z_y_0_zz, g_zz_0_0_xx, g_zz_0_0_xy, g_zz_0_0_xz, g_zz_0_0_yy, g_zz_0_0_yz, g_zz_0_0_zz, g_zz_yy_0_xx, g_zz_yy_0_xy, g_zz_yy_0_xz, g_zz_yy_0_yy, g_zz_yy_0_yz, g_zz_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_y_0_xx[i] = g_0_0_0_xx[i] - 2.0 * g_0_yy_0_xx[i] * b_exp - 2.0 * g_zz_0_0_xx[i] * a_exp + 4.0 * g_zz_yy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_0_xy[i] = g_0_0_0_xy[i] - 2.0 * g_0_yy_0_xy[i] * b_exp - 2.0 * g_zz_0_0_xy[i] * a_exp + 4.0 * g_zz_yy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_0_xz[i] = g_0_0_0_xz[i] - 2.0 * g_0_yy_0_xz[i] * b_exp - 2.0 * g_zz_0_0_xz[i] * a_exp + 4.0 * g_zz_yy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_0_yy[i] = g_0_0_0_yy[i] - 2.0 * g_0_yy_0_yy[i] * b_exp - 2.0 * g_zz_0_0_yy[i] * a_exp + 4.0 * g_zz_yy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_0_yz[i] = g_0_0_0_yz[i] - 2.0 * g_0_yy_0_yz[i] * b_exp - 2.0 * g_zz_0_0_yz[i] * a_exp + 4.0 * g_zz_yy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_0_zz[i] = g_0_0_0_zz[i] - 2.0 * g_0_yy_0_zz[i] * b_exp - 2.0 * g_zz_0_0_zz[i] * a_exp + 4.0 * g_zz_yy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (426-432)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_z_y_0_0_z_z_0_xx, g_z_y_0_0_z_z_0_xy, g_z_y_0_0_z_z_0_xz, g_z_y_0_0_z_z_0_yy, g_z_y_0_0_z_z_0_yz, g_z_y_0_0_z_z_0_zz, g_zz_yz_0_xx, g_zz_yz_0_xy, g_zz_yz_0_xz, g_zz_yz_0_yy, g_zz_yz_0_yz, g_zz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_z_0_xx[i] = -2.0 * g_0_yz_0_xx[i] * b_exp + 4.0 * g_zz_yz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_0_xy[i] = -2.0 * g_0_yz_0_xy[i] * b_exp + 4.0 * g_zz_yz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_0_xz[i] = -2.0 * g_0_yz_0_xz[i] * b_exp + 4.0 * g_zz_yz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_0_yy[i] = -2.0 * g_0_yz_0_yy[i] * b_exp + 4.0 * g_zz_yz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_0_yz[i] = -2.0 * g_0_yz_0_yz[i] * b_exp + 4.0 * g_zz_yz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_0_zz[i] = -2.0 * g_0_yz_0_zz[i] * b_exp + 4.0 * g_zz_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (432-438)

    #pragma omp simd aligned(g_xz_xz_0_xx, g_xz_xz_0_xy, g_xz_xz_0_xz, g_xz_xz_0_yy, g_xz_xz_0_yz, g_xz_xz_0_zz, g_z_z_0_0_x_x_0_xx, g_z_z_0_0_x_x_0_xy, g_z_z_0_0_x_x_0_xz, g_z_z_0_0_x_x_0_yy, g_z_z_0_0_x_x_0_yz, g_z_z_0_0_x_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_x_0_xx[i] = 4.0 * g_xz_xz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_0_xy[i] = 4.0 * g_xz_xz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_0_xz[i] = 4.0 * g_xz_xz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_0_yy[i] = 4.0 * g_xz_xz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_0_yz[i] = 4.0 * g_xz_xz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_0_zz[i] = 4.0 * g_xz_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (438-444)

    #pragma omp simd aligned(g_xz_yz_0_xx, g_xz_yz_0_xy, g_xz_yz_0_xz, g_xz_yz_0_yy, g_xz_yz_0_yz, g_xz_yz_0_zz, g_z_z_0_0_x_y_0_xx, g_z_z_0_0_x_y_0_xy, g_z_z_0_0_x_y_0_xz, g_z_z_0_0_x_y_0_yy, g_z_z_0_0_x_y_0_yz, g_z_z_0_0_x_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_y_0_xx[i] = 4.0 * g_xz_yz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_0_xy[i] = 4.0 * g_xz_yz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_0_xz[i] = 4.0 * g_xz_yz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_0_yy[i] = 4.0 * g_xz_yz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_0_yz[i] = 4.0 * g_xz_yz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_0_zz[i] = 4.0 * g_xz_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (444-450)

    #pragma omp simd aligned(g_xz_0_0_xx, g_xz_0_0_xy, g_xz_0_0_xz, g_xz_0_0_yy, g_xz_0_0_yz, g_xz_0_0_zz, g_xz_zz_0_xx, g_xz_zz_0_xy, g_xz_zz_0_xz, g_xz_zz_0_yy, g_xz_zz_0_yz, g_xz_zz_0_zz, g_z_z_0_0_x_z_0_xx, g_z_z_0_0_x_z_0_xy, g_z_z_0_0_x_z_0_xz, g_z_z_0_0_x_z_0_yy, g_z_z_0_0_x_z_0_yz, g_z_z_0_0_x_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_z_0_xx[i] = -2.0 * g_xz_0_0_xx[i] * a_exp + 4.0 * g_xz_zz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_0_xy[i] = -2.0 * g_xz_0_0_xy[i] * a_exp + 4.0 * g_xz_zz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_0_xz[i] = -2.0 * g_xz_0_0_xz[i] * a_exp + 4.0 * g_xz_zz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_0_yy[i] = -2.0 * g_xz_0_0_yy[i] * a_exp + 4.0 * g_xz_zz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_0_yz[i] = -2.0 * g_xz_0_0_yz[i] * a_exp + 4.0 * g_xz_zz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_0_zz[i] = -2.0 * g_xz_0_0_zz[i] * a_exp + 4.0 * g_xz_zz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (450-456)

    #pragma omp simd aligned(g_yz_xz_0_xx, g_yz_xz_0_xy, g_yz_xz_0_xz, g_yz_xz_0_yy, g_yz_xz_0_yz, g_yz_xz_0_zz, g_z_z_0_0_y_x_0_xx, g_z_z_0_0_y_x_0_xy, g_z_z_0_0_y_x_0_xz, g_z_z_0_0_y_x_0_yy, g_z_z_0_0_y_x_0_yz, g_z_z_0_0_y_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_x_0_xx[i] = 4.0 * g_yz_xz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_0_xy[i] = 4.0 * g_yz_xz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_0_xz[i] = 4.0 * g_yz_xz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_0_yy[i] = 4.0 * g_yz_xz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_0_yz[i] = 4.0 * g_yz_xz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_0_zz[i] = 4.0 * g_yz_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (456-462)

    #pragma omp simd aligned(g_yz_yz_0_xx, g_yz_yz_0_xy, g_yz_yz_0_xz, g_yz_yz_0_yy, g_yz_yz_0_yz, g_yz_yz_0_zz, g_z_z_0_0_y_y_0_xx, g_z_z_0_0_y_y_0_xy, g_z_z_0_0_y_y_0_xz, g_z_z_0_0_y_y_0_yy, g_z_z_0_0_y_y_0_yz, g_z_z_0_0_y_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_y_0_xx[i] = 4.0 * g_yz_yz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_0_xy[i] = 4.0 * g_yz_yz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_0_xz[i] = 4.0 * g_yz_yz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_0_yy[i] = 4.0 * g_yz_yz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_0_yz[i] = 4.0 * g_yz_yz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_0_zz[i] = 4.0 * g_yz_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (462-468)

    #pragma omp simd aligned(g_yz_0_0_xx, g_yz_0_0_xy, g_yz_0_0_xz, g_yz_0_0_yy, g_yz_0_0_yz, g_yz_0_0_zz, g_yz_zz_0_xx, g_yz_zz_0_xy, g_yz_zz_0_xz, g_yz_zz_0_yy, g_yz_zz_0_yz, g_yz_zz_0_zz, g_z_z_0_0_y_z_0_xx, g_z_z_0_0_y_z_0_xy, g_z_z_0_0_y_z_0_xz, g_z_z_0_0_y_z_0_yy, g_z_z_0_0_y_z_0_yz, g_z_z_0_0_y_z_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_z_0_xx[i] = -2.0 * g_yz_0_0_xx[i] * a_exp + 4.0 * g_yz_zz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_0_xy[i] = -2.0 * g_yz_0_0_xy[i] * a_exp + 4.0 * g_yz_zz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_0_xz[i] = -2.0 * g_yz_0_0_xz[i] * a_exp + 4.0 * g_yz_zz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_0_yy[i] = -2.0 * g_yz_0_0_yy[i] * a_exp + 4.0 * g_yz_zz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_0_yz[i] = -2.0 * g_yz_0_0_yz[i] * a_exp + 4.0 * g_yz_zz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_0_zz[i] = -2.0 * g_yz_0_0_zz[i] * a_exp + 4.0 * g_yz_zz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (468-474)

    #pragma omp simd aligned(g_0_xz_0_xx, g_0_xz_0_xy, g_0_xz_0_xz, g_0_xz_0_yy, g_0_xz_0_yz, g_0_xz_0_zz, g_z_z_0_0_z_x_0_xx, g_z_z_0_0_z_x_0_xy, g_z_z_0_0_z_x_0_xz, g_z_z_0_0_z_x_0_yy, g_z_z_0_0_z_x_0_yz, g_z_z_0_0_z_x_0_zz, g_zz_xz_0_xx, g_zz_xz_0_xy, g_zz_xz_0_xz, g_zz_xz_0_yy, g_zz_xz_0_yz, g_zz_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_x_0_xx[i] = -2.0 * g_0_xz_0_xx[i] * b_exp + 4.0 * g_zz_xz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_0_xy[i] = -2.0 * g_0_xz_0_xy[i] * b_exp + 4.0 * g_zz_xz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_0_xz[i] = -2.0 * g_0_xz_0_xz[i] * b_exp + 4.0 * g_zz_xz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_0_yy[i] = -2.0 * g_0_xz_0_yy[i] * b_exp + 4.0 * g_zz_xz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_0_yz[i] = -2.0 * g_0_xz_0_yz[i] * b_exp + 4.0 * g_zz_xz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_0_zz[i] = -2.0 * g_0_xz_0_zz[i] * b_exp + 4.0 * g_zz_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (474-480)

    #pragma omp simd aligned(g_0_yz_0_xx, g_0_yz_0_xy, g_0_yz_0_xz, g_0_yz_0_yy, g_0_yz_0_yz, g_0_yz_0_zz, g_z_z_0_0_z_y_0_xx, g_z_z_0_0_z_y_0_xy, g_z_z_0_0_z_y_0_xz, g_z_z_0_0_z_y_0_yy, g_z_z_0_0_z_y_0_yz, g_z_z_0_0_z_y_0_zz, g_zz_yz_0_xx, g_zz_yz_0_xy, g_zz_yz_0_xz, g_zz_yz_0_yy, g_zz_yz_0_yz, g_zz_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_y_0_xx[i] = -2.0 * g_0_yz_0_xx[i] * b_exp + 4.0 * g_zz_yz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_0_xy[i] = -2.0 * g_0_yz_0_xy[i] * b_exp + 4.0 * g_zz_yz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_0_xz[i] = -2.0 * g_0_yz_0_xz[i] * b_exp + 4.0 * g_zz_yz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_0_yy[i] = -2.0 * g_0_yz_0_yy[i] * b_exp + 4.0 * g_zz_yz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_0_yz[i] = -2.0 * g_0_yz_0_yz[i] * b_exp + 4.0 * g_zz_yz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_0_zz[i] = -2.0 * g_0_yz_0_zz[i] * b_exp + 4.0 * g_zz_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (480-486)

    #pragma omp simd aligned(g_0_0_0_xx, g_0_0_0_xy, g_0_0_0_xz, g_0_0_0_yy, g_0_0_0_yz, g_0_0_0_zz, g_0_zz_0_xx, g_0_zz_0_xy, g_0_zz_0_xz, g_0_zz_0_yy, g_0_zz_0_yz, g_0_zz_0_zz, g_z_z_0_0_z_z_0_xx, g_z_z_0_0_z_z_0_xy, g_z_z_0_0_z_z_0_xz, g_z_z_0_0_z_z_0_yy, g_z_z_0_0_z_z_0_yz, g_z_z_0_0_z_z_0_zz, g_zz_0_0_xx, g_zz_0_0_xy, g_zz_0_0_xz, g_zz_0_0_yy, g_zz_0_0_yz, g_zz_0_0_zz, g_zz_zz_0_xx, g_zz_zz_0_xy, g_zz_zz_0_xz, g_zz_zz_0_yy, g_zz_zz_0_yz, g_zz_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_z_0_xx[i] = g_0_0_0_xx[i] - 2.0 * g_0_zz_0_xx[i] * b_exp - 2.0 * g_zz_0_0_xx[i] * a_exp + 4.0 * g_zz_zz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_0_xy[i] = g_0_0_0_xy[i] - 2.0 * g_0_zz_0_xy[i] * b_exp - 2.0 * g_zz_0_0_xy[i] * a_exp + 4.0 * g_zz_zz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_0_xz[i] = g_0_0_0_xz[i] - 2.0 * g_0_zz_0_xz[i] * b_exp - 2.0 * g_zz_0_0_xz[i] * a_exp + 4.0 * g_zz_zz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_0_yy[i] = g_0_0_0_yy[i] - 2.0 * g_0_zz_0_yy[i] * b_exp - 2.0 * g_zz_0_0_yy[i] * a_exp + 4.0 * g_zz_zz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_0_yz[i] = g_0_0_0_yz[i] - 2.0 * g_0_zz_0_yz[i] * b_exp - 2.0 * g_zz_0_0_yz[i] * a_exp + 4.0 * g_zz_zz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_0_zz[i] = g_0_0_0_zz[i] - 2.0 * g_0_zz_0_zz[i] * b_exp - 2.0 * g_zz_0_0_zz[i] * a_exp + 4.0 * g_zz_zz_0_zz[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

