#include "GeomDeriv1100OfScalarForSPPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_sppp_0(CSimdArray<double>& buffer_1100_sppp,
                     const CSimdArray<double>& buffer_pspp,
                     const CSimdArray<double>& buffer_pdpp,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_sppp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pspp

    auto g_x_0_x_x = buffer_pspp[0];

    auto g_x_0_x_y = buffer_pspp[1];

    auto g_x_0_x_z = buffer_pspp[2];

    auto g_x_0_y_x = buffer_pspp[3];

    auto g_x_0_y_y = buffer_pspp[4];

    auto g_x_0_y_z = buffer_pspp[5];

    auto g_x_0_z_x = buffer_pspp[6];

    auto g_x_0_z_y = buffer_pspp[7];

    auto g_x_0_z_z = buffer_pspp[8];

    auto g_y_0_x_x = buffer_pspp[9];

    auto g_y_0_x_y = buffer_pspp[10];

    auto g_y_0_x_z = buffer_pspp[11];

    auto g_y_0_y_x = buffer_pspp[12];

    auto g_y_0_y_y = buffer_pspp[13];

    auto g_y_0_y_z = buffer_pspp[14];

    auto g_y_0_z_x = buffer_pspp[15];

    auto g_y_0_z_y = buffer_pspp[16];

    auto g_y_0_z_z = buffer_pspp[17];

    auto g_z_0_x_x = buffer_pspp[18];

    auto g_z_0_x_y = buffer_pspp[19];

    auto g_z_0_x_z = buffer_pspp[20];

    auto g_z_0_y_x = buffer_pspp[21];

    auto g_z_0_y_y = buffer_pspp[22];

    auto g_z_0_y_z = buffer_pspp[23];

    auto g_z_0_z_x = buffer_pspp[24];

    auto g_z_0_z_y = buffer_pspp[25];

    auto g_z_0_z_z = buffer_pspp[26];

    /// Set up components of auxilary buffer : buffer_pdpp

    auto g_x_xx_x_x = buffer_pdpp[0];

    auto g_x_xx_x_y = buffer_pdpp[1];

    auto g_x_xx_x_z = buffer_pdpp[2];

    auto g_x_xx_y_x = buffer_pdpp[3];

    auto g_x_xx_y_y = buffer_pdpp[4];

    auto g_x_xx_y_z = buffer_pdpp[5];

    auto g_x_xx_z_x = buffer_pdpp[6];

    auto g_x_xx_z_y = buffer_pdpp[7];

    auto g_x_xx_z_z = buffer_pdpp[8];

    auto g_x_xy_x_x = buffer_pdpp[9];

    auto g_x_xy_x_y = buffer_pdpp[10];

    auto g_x_xy_x_z = buffer_pdpp[11];

    auto g_x_xy_y_x = buffer_pdpp[12];

    auto g_x_xy_y_y = buffer_pdpp[13];

    auto g_x_xy_y_z = buffer_pdpp[14];

    auto g_x_xy_z_x = buffer_pdpp[15];

    auto g_x_xy_z_y = buffer_pdpp[16];

    auto g_x_xy_z_z = buffer_pdpp[17];

    auto g_x_xz_x_x = buffer_pdpp[18];

    auto g_x_xz_x_y = buffer_pdpp[19];

    auto g_x_xz_x_z = buffer_pdpp[20];

    auto g_x_xz_y_x = buffer_pdpp[21];

    auto g_x_xz_y_y = buffer_pdpp[22];

    auto g_x_xz_y_z = buffer_pdpp[23];

    auto g_x_xz_z_x = buffer_pdpp[24];

    auto g_x_xz_z_y = buffer_pdpp[25];

    auto g_x_xz_z_z = buffer_pdpp[26];

    auto g_x_yy_x_x = buffer_pdpp[27];

    auto g_x_yy_x_y = buffer_pdpp[28];

    auto g_x_yy_x_z = buffer_pdpp[29];

    auto g_x_yy_y_x = buffer_pdpp[30];

    auto g_x_yy_y_y = buffer_pdpp[31];

    auto g_x_yy_y_z = buffer_pdpp[32];

    auto g_x_yy_z_x = buffer_pdpp[33];

    auto g_x_yy_z_y = buffer_pdpp[34];

    auto g_x_yy_z_z = buffer_pdpp[35];

    auto g_x_yz_x_x = buffer_pdpp[36];

    auto g_x_yz_x_y = buffer_pdpp[37];

    auto g_x_yz_x_z = buffer_pdpp[38];

    auto g_x_yz_y_x = buffer_pdpp[39];

    auto g_x_yz_y_y = buffer_pdpp[40];

    auto g_x_yz_y_z = buffer_pdpp[41];

    auto g_x_yz_z_x = buffer_pdpp[42];

    auto g_x_yz_z_y = buffer_pdpp[43];

    auto g_x_yz_z_z = buffer_pdpp[44];

    auto g_x_zz_x_x = buffer_pdpp[45];

    auto g_x_zz_x_y = buffer_pdpp[46];

    auto g_x_zz_x_z = buffer_pdpp[47];

    auto g_x_zz_y_x = buffer_pdpp[48];

    auto g_x_zz_y_y = buffer_pdpp[49];

    auto g_x_zz_y_z = buffer_pdpp[50];

    auto g_x_zz_z_x = buffer_pdpp[51];

    auto g_x_zz_z_y = buffer_pdpp[52];

    auto g_x_zz_z_z = buffer_pdpp[53];

    auto g_y_xx_x_x = buffer_pdpp[54];

    auto g_y_xx_x_y = buffer_pdpp[55];

    auto g_y_xx_x_z = buffer_pdpp[56];

    auto g_y_xx_y_x = buffer_pdpp[57];

    auto g_y_xx_y_y = buffer_pdpp[58];

    auto g_y_xx_y_z = buffer_pdpp[59];

    auto g_y_xx_z_x = buffer_pdpp[60];

    auto g_y_xx_z_y = buffer_pdpp[61];

    auto g_y_xx_z_z = buffer_pdpp[62];

    auto g_y_xy_x_x = buffer_pdpp[63];

    auto g_y_xy_x_y = buffer_pdpp[64];

    auto g_y_xy_x_z = buffer_pdpp[65];

    auto g_y_xy_y_x = buffer_pdpp[66];

    auto g_y_xy_y_y = buffer_pdpp[67];

    auto g_y_xy_y_z = buffer_pdpp[68];

    auto g_y_xy_z_x = buffer_pdpp[69];

    auto g_y_xy_z_y = buffer_pdpp[70];

    auto g_y_xy_z_z = buffer_pdpp[71];

    auto g_y_xz_x_x = buffer_pdpp[72];

    auto g_y_xz_x_y = buffer_pdpp[73];

    auto g_y_xz_x_z = buffer_pdpp[74];

    auto g_y_xz_y_x = buffer_pdpp[75];

    auto g_y_xz_y_y = buffer_pdpp[76];

    auto g_y_xz_y_z = buffer_pdpp[77];

    auto g_y_xz_z_x = buffer_pdpp[78];

    auto g_y_xz_z_y = buffer_pdpp[79];

    auto g_y_xz_z_z = buffer_pdpp[80];

    auto g_y_yy_x_x = buffer_pdpp[81];

    auto g_y_yy_x_y = buffer_pdpp[82];

    auto g_y_yy_x_z = buffer_pdpp[83];

    auto g_y_yy_y_x = buffer_pdpp[84];

    auto g_y_yy_y_y = buffer_pdpp[85];

    auto g_y_yy_y_z = buffer_pdpp[86];

    auto g_y_yy_z_x = buffer_pdpp[87];

    auto g_y_yy_z_y = buffer_pdpp[88];

    auto g_y_yy_z_z = buffer_pdpp[89];

    auto g_y_yz_x_x = buffer_pdpp[90];

    auto g_y_yz_x_y = buffer_pdpp[91];

    auto g_y_yz_x_z = buffer_pdpp[92];

    auto g_y_yz_y_x = buffer_pdpp[93];

    auto g_y_yz_y_y = buffer_pdpp[94];

    auto g_y_yz_y_z = buffer_pdpp[95];

    auto g_y_yz_z_x = buffer_pdpp[96];

    auto g_y_yz_z_y = buffer_pdpp[97];

    auto g_y_yz_z_z = buffer_pdpp[98];

    auto g_y_zz_x_x = buffer_pdpp[99];

    auto g_y_zz_x_y = buffer_pdpp[100];

    auto g_y_zz_x_z = buffer_pdpp[101];

    auto g_y_zz_y_x = buffer_pdpp[102];

    auto g_y_zz_y_y = buffer_pdpp[103];

    auto g_y_zz_y_z = buffer_pdpp[104];

    auto g_y_zz_z_x = buffer_pdpp[105];

    auto g_y_zz_z_y = buffer_pdpp[106];

    auto g_y_zz_z_z = buffer_pdpp[107];

    auto g_z_xx_x_x = buffer_pdpp[108];

    auto g_z_xx_x_y = buffer_pdpp[109];

    auto g_z_xx_x_z = buffer_pdpp[110];

    auto g_z_xx_y_x = buffer_pdpp[111];

    auto g_z_xx_y_y = buffer_pdpp[112];

    auto g_z_xx_y_z = buffer_pdpp[113];

    auto g_z_xx_z_x = buffer_pdpp[114];

    auto g_z_xx_z_y = buffer_pdpp[115];

    auto g_z_xx_z_z = buffer_pdpp[116];

    auto g_z_xy_x_x = buffer_pdpp[117];

    auto g_z_xy_x_y = buffer_pdpp[118];

    auto g_z_xy_x_z = buffer_pdpp[119];

    auto g_z_xy_y_x = buffer_pdpp[120];

    auto g_z_xy_y_y = buffer_pdpp[121];

    auto g_z_xy_y_z = buffer_pdpp[122];

    auto g_z_xy_z_x = buffer_pdpp[123];

    auto g_z_xy_z_y = buffer_pdpp[124];

    auto g_z_xy_z_z = buffer_pdpp[125];

    auto g_z_xz_x_x = buffer_pdpp[126];

    auto g_z_xz_x_y = buffer_pdpp[127];

    auto g_z_xz_x_z = buffer_pdpp[128];

    auto g_z_xz_y_x = buffer_pdpp[129];

    auto g_z_xz_y_y = buffer_pdpp[130];

    auto g_z_xz_y_z = buffer_pdpp[131];

    auto g_z_xz_z_x = buffer_pdpp[132];

    auto g_z_xz_z_y = buffer_pdpp[133];

    auto g_z_xz_z_z = buffer_pdpp[134];

    auto g_z_yy_x_x = buffer_pdpp[135];

    auto g_z_yy_x_y = buffer_pdpp[136];

    auto g_z_yy_x_z = buffer_pdpp[137];

    auto g_z_yy_y_x = buffer_pdpp[138];

    auto g_z_yy_y_y = buffer_pdpp[139];

    auto g_z_yy_y_z = buffer_pdpp[140];

    auto g_z_yy_z_x = buffer_pdpp[141];

    auto g_z_yy_z_y = buffer_pdpp[142];

    auto g_z_yy_z_z = buffer_pdpp[143];

    auto g_z_yz_x_x = buffer_pdpp[144];

    auto g_z_yz_x_y = buffer_pdpp[145];

    auto g_z_yz_x_z = buffer_pdpp[146];

    auto g_z_yz_y_x = buffer_pdpp[147];

    auto g_z_yz_y_y = buffer_pdpp[148];

    auto g_z_yz_y_z = buffer_pdpp[149];

    auto g_z_yz_z_x = buffer_pdpp[150];

    auto g_z_yz_z_y = buffer_pdpp[151];

    auto g_z_yz_z_z = buffer_pdpp[152];

    auto g_z_zz_x_x = buffer_pdpp[153];

    auto g_z_zz_x_y = buffer_pdpp[154];

    auto g_z_zz_x_z = buffer_pdpp[155];

    auto g_z_zz_y_x = buffer_pdpp[156];

    auto g_z_zz_y_y = buffer_pdpp[157];

    auto g_z_zz_y_z = buffer_pdpp[158];

    auto g_z_zz_z_x = buffer_pdpp[159];

    auto g_z_zz_z_y = buffer_pdpp[160];

    auto g_z_zz_z_z = buffer_pdpp[161];

    /// Set up components of integrals buffer : buffer_1100_sppp

    auto g_x_x_0_0_0_x_x_x = buffer_1100_sppp[0];

    auto g_x_x_0_0_0_x_x_y = buffer_1100_sppp[1];

    auto g_x_x_0_0_0_x_x_z = buffer_1100_sppp[2];

    auto g_x_x_0_0_0_x_y_x = buffer_1100_sppp[3];

    auto g_x_x_0_0_0_x_y_y = buffer_1100_sppp[4];

    auto g_x_x_0_0_0_x_y_z = buffer_1100_sppp[5];

    auto g_x_x_0_0_0_x_z_x = buffer_1100_sppp[6];

    auto g_x_x_0_0_0_x_z_y = buffer_1100_sppp[7];

    auto g_x_x_0_0_0_x_z_z = buffer_1100_sppp[8];

    auto g_x_x_0_0_0_y_x_x = buffer_1100_sppp[9];

    auto g_x_x_0_0_0_y_x_y = buffer_1100_sppp[10];

    auto g_x_x_0_0_0_y_x_z = buffer_1100_sppp[11];

    auto g_x_x_0_0_0_y_y_x = buffer_1100_sppp[12];

    auto g_x_x_0_0_0_y_y_y = buffer_1100_sppp[13];

    auto g_x_x_0_0_0_y_y_z = buffer_1100_sppp[14];

    auto g_x_x_0_0_0_y_z_x = buffer_1100_sppp[15];

    auto g_x_x_0_0_0_y_z_y = buffer_1100_sppp[16];

    auto g_x_x_0_0_0_y_z_z = buffer_1100_sppp[17];

    auto g_x_x_0_0_0_z_x_x = buffer_1100_sppp[18];

    auto g_x_x_0_0_0_z_x_y = buffer_1100_sppp[19];

    auto g_x_x_0_0_0_z_x_z = buffer_1100_sppp[20];

    auto g_x_x_0_0_0_z_y_x = buffer_1100_sppp[21];

    auto g_x_x_0_0_0_z_y_y = buffer_1100_sppp[22];

    auto g_x_x_0_0_0_z_y_z = buffer_1100_sppp[23];

    auto g_x_x_0_0_0_z_z_x = buffer_1100_sppp[24];

    auto g_x_x_0_0_0_z_z_y = buffer_1100_sppp[25];

    auto g_x_x_0_0_0_z_z_z = buffer_1100_sppp[26];

    auto g_x_y_0_0_0_x_x_x = buffer_1100_sppp[27];

    auto g_x_y_0_0_0_x_x_y = buffer_1100_sppp[28];

    auto g_x_y_0_0_0_x_x_z = buffer_1100_sppp[29];

    auto g_x_y_0_0_0_x_y_x = buffer_1100_sppp[30];

    auto g_x_y_0_0_0_x_y_y = buffer_1100_sppp[31];

    auto g_x_y_0_0_0_x_y_z = buffer_1100_sppp[32];

    auto g_x_y_0_0_0_x_z_x = buffer_1100_sppp[33];

    auto g_x_y_0_0_0_x_z_y = buffer_1100_sppp[34];

    auto g_x_y_0_0_0_x_z_z = buffer_1100_sppp[35];

    auto g_x_y_0_0_0_y_x_x = buffer_1100_sppp[36];

    auto g_x_y_0_0_0_y_x_y = buffer_1100_sppp[37];

    auto g_x_y_0_0_0_y_x_z = buffer_1100_sppp[38];

    auto g_x_y_0_0_0_y_y_x = buffer_1100_sppp[39];

    auto g_x_y_0_0_0_y_y_y = buffer_1100_sppp[40];

    auto g_x_y_0_0_0_y_y_z = buffer_1100_sppp[41];

    auto g_x_y_0_0_0_y_z_x = buffer_1100_sppp[42];

    auto g_x_y_0_0_0_y_z_y = buffer_1100_sppp[43];

    auto g_x_y_0_0_0_y_z_z = buffer_1100_sppp[44];

    auto g_x_y_0_0_0_z_x_x = buffer_1100_sppp[45];

    auto g_x_y_0_0_0_z_x_y = buffer_1100_sppp[46];

    auto g_x_y_0_0_0_z_x_z = buffer_1100_sppp[47];

    auto g_x_y_0_0_0_z_y_x = buffer_1100_sppp[48];

    auto g_x_y_0_0_0_z_y_y = buffer_1100_sppp[49];

    auto g_x_y_0_0_0_z_y_z = buffer_1100_sppp[50];

    auto g_x_y_0_0_0_z_z_x = buffer_1100_sppp[51];

    auto g_x_y_0_0_0_z_z_y = buffer_1100_sppp[52];

    auto g_x_y_0_0_0_z_z_z = buffer_1100_sppp[53];

    auto g_x_z_0_0_0_x_x_x = buffer_1100_sppp[54];

    auto g_x_z_0_0_0_x_x_y = buffer_1100_sppp[55];

    auto g_x_z_0_0_0_x_x_z = buffer_1100_sppp[56];

    auto g_x_z_0_0_0_x_y_x = buffer_1100_sppp[57];

    auto g_x_z_0_0_0_x_y_y = buffer_1100_sppp[58];

    auto g_x_z_0_0_0_x_y_z = buffer_1100_sppp[59];

    auto g_x_z_0_0_0_x_z_x = buffer_1100_sppp[60];

    auto g_x_z_0_0_0_x_z_y = buffer_1100_sppp[61];

    auto g_x_z_0_0_0_x_z_z = buffer_1100_sppp[62];

    auto g_x_z_0_0_0_y_x_x = buffer_1100_sppp[63];

    auto g_x_z_0_0_0_y_x_y = buffer_1100_sppp[64];

    auto g_x_z_0_0_0_y_x_z = buffer_1100_sppp[65];

    auto g_x_z_0_0_0_y_y_x = buffer_1100_sppp[66];

    auto g_x_z_0_0_0_y_y_y = buffer_1100_sppp[67];

    auto g_x_z_0_0_0_y_y_z = buffer_1100_sppp[68];

    auto g_x_z_0_0_0_y_z_x = buffer_1100_sppp[69];

    auto g_x_z_0_0_0_y_z_y = buffer_1100_sppp[70];

    auto g_x_z_0_0_0_y_z_z = buffer_1100_sppp[71];

    auto g_x_z_0_0_0_z_x_x = buffer_1100_sppp[72];

    auto g_x_z_0_0_0_z_x_y = buffer_1100_sppp[73];

    auto g_x_z_0_0_0_z_x_z = buffer_1100_sppp[74];

    auto g_x_z_0_0_0_z_y_x = buffer_1100_sppp[75];

    auto g_x_z_0_0_0_z_y_y = buffer_1100_sppp[76];

    auto g_x_z_0_0_0_z_y_z = buffer_1100_sppp[77];

    auto g_x_z_0_0_0_z_z_x = buffer_1100_sppp[78];

    auto g_x_z_0_0_0_z_z_y = buffer_1100_sppp[79];

    auto g_x_z_0_0_0_z_z_z = buffer_1100_sppp[80];

    auto g_y_x_0_0_0_x_x_x = buffer_1100_sppp[81];

    auto g_y_x_0_0_0_x_x_y = buffer_1100_sppp[82];

    auto g_y_x_0_0_0_x_x_z = buffer_1100_sppp[83];

    auto g_y_x_0_0_0_x_y_x = buffer_1100_sppp[84];

    auto g_y_x_0_0_0_x_y_y = buffer_1100_sppp[85];

    auto g_y_x_0_0_0_x_y_z = buffer_1100_sppp[86];

    auto g_y_x_0_0_0_x_z_x = buffer_1100_sppp[87];

    auto g_y_x_0_0_0_x_z_y = buffer_1100_sppp[88];

    auto g_y_x_0_0_0_x_z_z = buffer_1100_sppp[89];

    auto g_y_x_0_0_0_y_x_x = buffer_1100_sppp[90];

    auto g_y_x_0_0_0_y_x_y = buffer_1100_sppp[91];

    auto g_y_x_0_0_0_y_x_z = buffer_1100_sppp[92];

    auto g_y_x_0_0_0_y_y_x = buffer_1100_sppp[93];

    auto g_y_x_0_0_0_y_y_y = buffer_1100_sppp[94];

    auto g_y_x_0_0_0_y_y_z = buffer_1100_sppp[95];

    auto g_y_x_0_0_0_y_z_x = buffer_1100_sppp[96];

    auto g_y_x_0_0_0_y_z_y = buffer_1100_sppp[97];

    auto g_y_x_0_0_0_y_z_z = buffer_1100_sppp[98];

    auto g_y_x_0_0_0_z_x_x = buffer_1100_sppp[99];

    auto g_y_x_0_0_0_z_x_y = buffer_1100_sppp[100];

    auto g_y_x_0_0_0_z_x_z = buffer_1100_sppp[101];

    auto g_y_x_0_0_0_z_y_x = buffer_1100_sppp[102];

    auto g_y_x_0_0_0_z_y_y = buffer_1100_sppp[103];

    auto g_y_x_0_0_0_z_y_z = buffer_1100_sppp[104];

    auto g_y_x_0_0_0_z_z_x = buffer_1100_sppp[105];

    auto g_y_x_0_0_0_z_z_y = buffer_1100_sppp[106];

    auto g_y_x_0_0_0_z_z_z = buffer_1100_sppp[107];

    auto g_y_y_0_0_0_x_x_x = buffer_1100_sppp[108];

    auto g_y_y_0_0_0_x_x_y = buffer_1100_sppp[109];

    auto g_y_y_0_0_0_x_x_z = buffer_1100_sppp[110];

    auto g_y_y_0_0_0_x_y_x = buffer_1100_sppp[111];

    auto g_y_y_0_0_0_x_y_y = buffer_1100_sppp[112];

    auto g_y_y_0_0_0_x_y_z = buffer_1100_sppp[113];

    auto g_y_y_0_0_0_x_z_x = buffer_1100_sppp[114];

    auto g_y_y_0_0_0_x_z_y = buffer_1100_sppp[115];

    auto g_y_y_0_0_0_x_z_z = buffer_1100_sppp[116];

    auto g_y_y_0_0_0_y_x_x = buffer_1100_sppp[117];

    auto g_y_y_0_0_0_y_x_y = buffer_1100_sppp[118];

    auto g_y_y_0_0_0_y_x_z = buffer_1100_sppp[119];

    auto g_y_y_0_0_0_y_y_x = buffer_1100_sppp[120];

    auto g_y_y_0_0_0_y_y_y = buffer_1100_sppp[121];

    auto g_y_y_0_0_0_y_y_z = buffer_1100_sppp[122];

    auto g_y_y_0_0_0_y_z_x = buffer_1100_sppp[123];

    auto g_y_y_0_0_0_y_z_y = buffer_1100_sppp[124];

    auto g_y_y_0_0_0_y_z_z = buffer_1100_sppp[125];

    auto g_y_y_0_0_0_z_x_x = buffer_1100_sppp[126];

    auto g_y_y_0_0_0_z_x_y = buffer_1100_sppp[127];

    auto g_y_y_0_0_0_z_x_z = buffer_1100_sppp[128];

    auto g_y_y_0_0_0_z_y_x = buffer_1100_sppp[129];

    auto g_y_y_0_0_0_z_y_y = buffer_1100_sppp[130];

    auto g_y_y_0_0_0_z_y_z = buffer_1100_sppp[131];

    auto g_y_y_0_0_0_z_z_x = buffer_1100_sppp[132];

    auto g_y_y_0_0_0_z_z_y = buffer_1100_sppp[133];

    auto g_y_y_0_0_0_z_z_z = buffer_1100_sppp[134];

    auto g_y_z_0_0_0_x_x_x = buffer_1100_sppp[135];

    auto g_y_z_0_0_0_x_x_y = buffer_1100_sppp[136];

    auto g_y_z_0_0_0_x_x_z = buffer_1100_sppp[137];

    auto g_y_z_0_0_0_x_y_x = buffer_1100_sppp[138];

    auto g_y_z_0_0_0_x_y_y = buffer_1100_sppp[139];

    auto g_y_z_0_0_0_x_y_z = buffer_1100_sppp[140];

    auto g_y_z_0_0_0_x_z_x = buffer_1100_sppp[141];

    auto g_y_z_0_0_0_x_z_y = buffer_1100_sppp[142];

    auto g_y_z_0_0_0_x_z_z = buffer_1100_sppp[143];

    auto g_y_z_0_0_0_y_x_x = buffer_1100_sppp[144];

    auto g_y_z_0_0_0_y_x_y = buffer_1100_sppp[145];

    auto g_y_z_0_0_0_y_x_z = buffer_1100_sppp[146];

    auto g_y_z_0_0_0_y_y_x = buffer_1100_sppp[147];

    auto g_y_z_0_0_0_y_y_y = buffer_1100_sppp[148];

    auto g_y_z_0_0_0_y_y_z = buffer_1100_sppp[149];

    auto g_y_z_0_0_0_y_z_x = buffer_1100_sppp[150];

    auto g_y_z_0_0_0_y_z_y = buffer_1100_sppp[151];

    auto g_y_z_0_0_0_y_z_z = buffer_1100_sppp[152];

    auto g_y_z_0_0_0_z_x_x = buffer_1100_sppp[153];

    auto g_y_z_0_0_0_z_x_y = buffer_1100_sppp[154];

    auto g_y_z_0_0_0_z_x_z = buffer_1100_sppp[155];

    auto g_y_z_0_0_0_z_y_x = buffer_1100_sppp[156];

    auto g_y_z_0_0_0_z_y_y = buffer_1100_sppp[157];

    auto g_y_z_0_0_0_z_y_z = buffer_1100_sppp[158];

    auto g_y_z_0_0_0_z_z_x = buffer_1100_sppp[159];

    auto g_y_z_0_0_0_z_z_y = buffer_1100_sppp[160];

    auto g_y_z_0_0_0_z_z_z = buffer_1100_sppp[161];

    auto g_z_x_0_0_0_x_x_x = buffer_1100_sppp[162];

    auto g_z_x_0_0_0_x_x_y = buffer_1100_sppp[163];

    auto g_z_x_0_0_0_x_x_z = buffer_1100_sppp[164];

    auto g_z_x_0_0_0_x_y_x = buffer_1100_sppp[165];

    auto g_z_x_0_0_0_x_y_y = buffer_1100_sppp[166];

    auto g_z_x_0_0_0_x_y_z = buffer_1100_sppp[167];

    auto g_z_x_0_0_0_x_z_x = buffer_1100_sppp[168];

    auto g_z_x_0_0_0_x_z_y = buffer_1100_sppp[169];

    auto g_z_x_0_0_0_x_z_z = buffer_1100_sppp[170];

    auto g_z_x_0_0_0_y_x_x = buffer_1100_sppp[171];

    auto g_z_x_0_0_0_y_x_y = buffer_1100_sppp[172];

    auto g_z_x_0_0_0_y_x_z = buffer_1100_sppp[173];

    auto g_z_x_0_0_0_y_y_x = buffer_1100_sppp[174];

    auto g_z_x_0_0_0_y_y_y = buffer_1100_sppp[175];

    auto g_z_x_0_0_0_y_y_z = buffer_1100_sppp[176];

    auto g_z_x_0_0_0_y_z_x = buffer_1100_sppp[177];

    auto g_z_x_0_0_0_y_z_y = buffer_1100_sppp[178];

    auto g_z_x_0_0_0_y_z_z = buffer_1100_sppp[179];

    auto g_z_x_0_0_0_z_x_x = buffer_1100_sppp[180];

    auto g_z_x_0_0_0_z_x_y = buffer_1100_sppp[181];

    auto g_z_x_0_0_0_z_x_z = buffer_1100_sppp[182];

    auto g_z_x_0_0_0_z_y_x = buffer_1100_sppp[183];

    auto g_z_x_0_0_0_z_y_y = buffer_1100_sppp[184];

    auto g_z_x_0_0_0_z_y_z = buffer_1100_sppp[185];

    auto g_z_x_0_0_0_z_z_x = buffer_1100_sppp[186];

    auto g_z_x_0_0_0_z_z_y = buffer_1100_sppp[187];

    auto g_z_x_0_0_0_z_z_z = buffer_1100_sppp[188];

    auto g_z_y_0_0_0_x_x_x = buffer_1100_sppp[189];

    auto g_z_y_0_0_0_x_x_y = buffer_1100_sppp[190];

    auto g_z_y_0_0_0_x_x_z = buffer_1100_sppp[191];

    auto g_z_y_0_0_0_x_y_x = buffer_1100_sppp[192];

    auto g_z_y_0_0_0_x_y_y = buffer_1100_sppp[193];

    auto g_z_y_0_0_0_x_y_z = buffer_1100_sppp[194];

    auto g_z_y_0_0_0_x_z_x = buffer_1100_sppp[195];

    auto g_z_y_0_0_0_x_z_y = buffer_1100_sppp[196];

    auto g_z_y_0_0_0_x_z_z = buffer_1100_sppp[197];

    auto g_z_y_0_0_0_y_x_x = buffer_1100_sppp[198];

    auto g_z_y_0_0_0_y_x_y = buffer_1100_sppp[199];

    auto g_z_y_0_0_0_y_x_z = buffer_1100_sppp[200];

    auto g_z_y_0_0_0_y_y_x = buffer_1100_sppp[201];

    auto g_z_y_0_0_0_y_y_y = buffer_1100_sppp[202];

    auto g_z_y_0_0_0_y_y_z = buffer_1100_sppp[203];

    auto g_z_y_0_0_0_y_z_x = buffer_1100_sppp[204];

    auto g_z_y_0_0_0_y_z_y = buffer_1100_sppp[205];

    auto g_z_y_0_0_0_y_z_z = buffer_1100_sppp[206];

    auto g_z_y_0_0_0_z_x_x = buffer_1100_sppp[207];

    auto g_z_y_0_0_0_z_x_y = buffer_1100_sppp[208];

    auto g_z_y_0_0_0_z_x_z = buffer_1100_sppp[209];

    auto g_z_y_0_0_0_z_y_x = buffer_1100_sppp[210];

    auto g_z_y_0_0_0_z_y_y = buffer_1100_sppp[211];

    auto g_z_y_0_0_0_z_y_z = buffer_1100_sppp[212];

    auto g_z_y_0_0_0_z_z_x = buffer_1100_sppp[213];

    auto g_z_y_0_0_0_z_z_y = buffer_1100_sppp[214];

    auto g_z_y_0_0_0_z_z_z = buffer_1100_sppp[215];

    auto g_z_z_0_0_0_x_x_x = buffer_1100_sppp[216];

    auto g_z_z_0_0_0_x_x_y = buffer_1100_sppp[217];

    auto g_z_z_0_0_0_x_x_z = buffer_1100_sppp[218];

    auto g_z_z_0_0_0_x_y_x = buffer_1100_sppp[219];

    auto g_z_z_0_0_0_x_y_y = buffer_1100_sppp[220];

    auto g_z_z_0_0_0_x_y_z = buffer_1100_sppp[221];

    auto g_z_z_0_0_0_x_z_x = buffer_1100_sppp[222];

    auto g_z_z_0_0_0_x_z_y = buffer_1100_sppp[223];

    auto g_z_z_0_0_0_x_z_z = buffer_1100_sppp[224];

    auto g_z_z_0_0_0_y_x_x = buffer_1100_sppp[225];

    auto g_z_z_0_0_0_y_x_y = buffer_1100_sppp[226];

    auto g_z_z_0_0_0_y_x_z = buffer_1100_sppp[227];

    auto g_z_z_0_0_0_y_y_x = buffer_1100_sppp[228];

    auto g_z_z_0_0_0_y_y_y = buffer_1100_sppp[229];

    auto g_z_z_0_0_0_y_y_z = buffer_1100_sppp[230];

    auto g_z_z_0_0_0_y_z_x = buffer_1100_sppp[231];

    auto g_z_z_0_0_0_y_z_y = buffer_1100_sppp[232];

    auto g_z_z_0_0_0_y_z_z = buffer_1100_sppp[233];

    auto g_z_z_0_0_0_z_x_x = buffer_1100_sppp[234];

    auto g_z_z_0_0_0_z_x_y = buffer_1100_sppp[235];

    auto g_z_z_0_0_0_z_x_z = buffer_1100_sppp[236];

    auto g_z_z_0_0_0_z_y_x = buffer_1100_sppp[237];

    auto g_z_z_0_0_0_z_y_y = buffer_1100_sppp[238];

    auto g_z_z_0_0_0_z_y_z = buffer_1100_sppp[239];

    auto g_z_z_0_0_0_z_z_x = buffer_1100_sppp[240];

    auto g_z_z_0_0_0_z_z_y = buffer_1100_sppp[241];

    auto g_z_z_0_0_0_z_z_z = buffer_1100_sppp[242];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_x_x, g_x_0_x_y, g_x_0_x_z, g_x_x_0_0_0_x_x_x, g_x_x_0_0_0_x_x_y, g_x_x_0_0_0_x_x_z, g_x_xx_x_x, g_x_xx_x_y, g_x_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_x_x_x[i] = -2.0 * g_x_0_x_x[i] * a_exp + 4.0 * g_x_xx_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_x_y[i] = -2.0 * g_x_0_x_y[i] * a_exp + 4.0 * g_x_xx_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_x_z[i] = -2.0 * g_x_0_x_z[i] * a_exp + 4.0 * g_x_xx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_y_x, g_x_0_y_y, g_x_0_y_z, g_x_x_0_0_0_x_y_x, g_x_x_0_0_0_x_y_y, g_x_x_0_0_0_x_y_z, g_x_xx_y_x, g_x_xx_y_y, g_x_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_x_y_x[i] = -2.0 * g_x_0_y_x[i] * a_exp + 4.0 * g_x_xx_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_y_y[i] = -2.0 * g_x_0_y_y[i] * a_exp + 4.0 * g_x_xx_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_y_z[i] = -2.0 * g_x_0_y_z[i] * a_exp + 4.0 * g_x_xx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_z_x, g_x_0_z_y, g_x_0_z_z, g_x_x_0_0_0_x_z_x, g_x_x_0_0_0_x_z_y, g_x_x_0_0_0_x_z_z, g_x_xx_z_x, g_x_xx_z_y, g_x_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_x_z_x[i] = -2.0 * g_x_0_z_x[i] * a_exp + 4.0 * g_x_xx_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_z_y[i] = -2.0 * g_x_0_z_y[i] * a_exp + 4.0 * g_x_xx_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_z_z[i] = -2.0 * g_x_0_z_z[i] * a_exp + 4.0 * g_x_xx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_x_0_0_0_y_x_x, g_x_x_0_0_0_y_x_y, g_x_x_0_0_0_y_x_z, g_x_xy_x_x, g_x_xy_x_y, g_x_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_y_x_x[i] = 4.0 * g_x_xy_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_x_y[i] = 4.0 * g_x_xy_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_x_z[i] = 4.0 * g_x_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_x_0_0_0_y_y_x, g_x_x_0_0_0_y_y_y, g_x_x_0_0_0_y_y_z, g_x_xy_y_x, g_x_xy_y_y, g_x_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_y_y_x[i] = 4.0 * g_x_xy_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_y_y[i] = 4.0 * g_x_xy_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_y_z[i] = 4.0 * g_x_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_x_0_0_0_y_z_x, g_x_x_0_0_0_y_z_y, g_x_x_0_0_0_y_z_z, g_x_xy_z_x, g_x_xy_z_y, g_x_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_y_z_x[i] = 4.0 * g_x_xy_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_z_y[i] = 4.0 * g_x_xy_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_z_z[i] = 4.0 * g_x_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_x_0_0_0_z_x_x, g_x_x_0_0_0_z_x_y, g_x_x_0_0_0_z_x_z, g_x_xz_x_x, g_x_xz_x_y, g_x_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_z_x_x[i] = 4.0 * g_x_xz_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_x_y[i] = 4.0 * g_x_xz_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_x_z[i] = 4.0 * g_x_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_x_0_0_0_z_y_x, g_x_x_0_0_0_z_y_y, g_x_x_0_0_0_z_y_z, g_x_xz_y_x, g_x_xz_y_y, g_x_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_z_y_x[i] = 4.0 * g_x_xz_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_y_y[i] = 4.0 * g_x_xz_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_y_z[i] = 4.0 * g_x_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_x_0_0_0_z_z_x, g_x_x_0_0_0_z_z_y, g_x_x_0_0_0_z_z_z, g_x_xz_z_x, g_x_xz_z_y, g_x_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_z_z_x[i] = 4.0 * g_x_xz_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_z_y[i] = 4.0 * g_x_xz_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_z_z[i] = 4.0 * g_x_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_xy_x_x, g_x_xy_x_y, g_x_xy_x_z, g_x_y_0_0_0_x_x_x, g_x_y_0_0_0_x_x_y, g_x_y_0_0_0_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_x_x_x[i] = 4.0 * g_x_xy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_x_y[i] = 4.0 * g_x_xy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_x_z[i] = 4.0 * g_x_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_x_xy_y_x, g_x_xy_y_y, g_x_xy_y_z, g_x_y_0_0_0_x_y_x, g_x_y_0_0_0_x_y_y, g_x_y_0_0_0_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_x_y_x[i] = 4.0 * g_x_xy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_y_y[i] = 4.0 * g_x_xy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_y_z[i] = 4.0 * g_x_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_x_xy_z_x, g_x_xy_z_y, g_x_xy_z_z, g_x_y_0_0_0_x_z_x, g_x_y_0_0_0_x_z_y, g_x_y_0_0_0_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_x_z_x[i] = 4.0 * g_x_xy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_z_y[i] = 4.0 * g_x_xy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_z_z[i] = 4.0 * g_x_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_0_x_x, g_x_0_x_y, g_x_0_x_z, g_x_y_0_0_0_y_x_x, g_x_y_0_0_0_y_x_y, g_x_y_0_0_0_y_x_z, g_x_yy_x_x, g_x_yy_x_y, g_x_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_y_x_x[i] = -2.0 * g_x_0_x_x[i] * a_exp + 4.0 * g_x_yy_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_x_y[i] = -2.0 * g_x_0_x_y[i] * a_exp + 4.0 * g_x_yy_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_x_z[i] = -2.0 * g_x_0_x_z[i] * a_exp + 4.0 * g_x_yy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_0_y_x, g_x_0_y_y, g_x_0_y_z, g_x_y_0_0_0_y_y_x, g_x_y_0_0_0_y_y_y, g_x_y_0_0_0_y_y_z, g_x_yy_y_x, g_x_yy_y_y, g_x_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_y_y_x[i] = -2.0 * g_x_0_y_x[i] * a_exp + 4.0 * g_x_yy_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_y_y[i] = -2.0 * g_x_0_y_y[i] * a_exp + 4.0 * g_x_yy_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_y_z[i] = -2.0 * g_x_0_y_z[i] * a_exp + 4.0 * g_x_yy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_0_z_x, g_x_0_z_y, g_x_0_z_z, g_x_y_0_0_0_y_z_x, g_x_y_0_0_0_y_z_y, g_x_y_0_0_0_y_z_z, g_x_yy_z_x, g_x_yy_z_y, g_x_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_y_z_x[i] = -2.0 * g_x_0_z_x[i] * a_exp + 4.0 * g_x_yy_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_z_y[i] = -2.0 * g_x_0_z_y[i] * a_exp + 4.0 * g_x_yy_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_z_z[i] = -2.0 * g_x_0_z_z[i] * a_exp + 4.0 * g_x_yy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_y_0_0_0_z_x_x, g_x_y_0_0_0_z_x_y, g_x_y_0_0_0_z_x_z, g_x_yz_x_x, g_x_yz_x_y, g_x_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_z_x_x[i] = 4.0 * g_x_yz_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_x_y[i] = 4.0 * g_x_yz_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_x_z[i] = 4.0 * g_x_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_y_0_0_0_z_y_x, g_x_y_0_0_0_z_y_y, g_x_y_0_0_0_z_y_z, g_x_yz_y_x, g_x_yz_y_y, g_x_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_z_y_x[i] = 4.0 * g_x_yz_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_y_y[i] = 4.0 * g_x_yz_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_y_z[i] = 4.0 * g_x_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_y_0_0_0_z_z_x, g_x_y_0_0_0_z_z_y, g_x_y_0_0_0_z_z_z, g_x_yz_z_x, g_x_yz_z_y, g_x_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_z_z_x[i] = 4.0 * g_x_yz_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_z_y[i] = 4.0 * g_x_yz_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_z_z[i] = 4.0 * g_x_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_x_xz_x_x, g_x_xz_x_y, g_x_xz_x_z, g_x_z_0_0_0_x_x_x, g_x_z_0_0_0_x_x_y, g_x_z_0_0_0_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_x_x_x[i] = 4.0 * g_x_xz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_x_y[i] = 4.0 * g_x_xz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_x_z[i] = 4.0 * g_x_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_x_xz_y_x, g_x_xz_y_y, g_x_xz_y_z, g_x_z_0_0_0_x_y_x, g_x_z_0_0_0_x_y_y, g_x_z_0_0_0_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_x_y_x[i] = 4.0 * g_x_xz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_y_y[i] = 4.0 * g_x_xz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_y_z[i] = 4.0 * g_x_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_x_xz_z_x, g_x_xz_z_y, g_x_xz_z_z, g_x_z_0_0_0_x_z_x, g_x_z_0_0_0_x_z_y, g_x_z_0_0_0_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_x_z_x[i] = 4.0 * g_x_xz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_z_y[i] = 4.0 * g_x_xz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_z_z[i] = 4.0 * g_x_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_x_yz_x_x, g_x_yz_x_y, g_x_yz_x_z, g_x_z_0_0_0_y_x_x, g_x_z_0_0_0_y_x_y, g_x_z_0_0_0_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_y_x_x[i] = 4.0 * g_x_yz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_x_y[i] = 4.0 * g_x_yz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_x_z[i] = 4.0 * g_x_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_x_yz_y_x, g_x_yz_y_y, g_x_yz_y_z, g_x_z_0_0_0_y_y_x, g_x_z_0_0_0_y_y_y, g_x_z_0_0_0_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_y_y_x[i] = 4.0 * g_x_yz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_y_y[i] = 4.0 * g_x_yz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_y_z[i] = 4.0 * g_x_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_x_yz_z_x, g_x_yz_z_y, g_x_yz_z_z, g_x_z_0_0_0_y_z_x, g_x_z_0_0_0_y_z_y, g_x_z_0_0_0_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_y_z_x[i] = 4.0 * g_x_yz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_z_y[i] = 4.0 * g_x_yz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_z_z[i] = 4.0 * g_x_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_0_x_x, g_x_0_x_y, g_x_0_x_z, g_x_z_0_0_0_z_x_x, g_x_z_0_0_0_z_x_y, g_x_z_0_0_0_z_x_z, g_x_zz_x_x, g_x_zz_x_y, g_x_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_z_x_x[i] = -2.0 * g_x_0_x_x[i] * a_exp + 4.0 * g_x_zz_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_x_y[i] = -2.0 * g_x_0_x_y[i] * a_exp + 4.0 * g_x_zz_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_x_z[i] = -2.0 * g_x_0_x_z[i] * a_exp + 4.0 * g_x_zz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_0_y_x, g_x_0_y_y, g_x_0_y_z, g_x_z_0_0_0_z_y_x, g_x_z_0_0_0_z_y_y, g_x_z_0_0_0_z_y_z, g_x_zz_y_x, g_x_zz_y_y, g_x_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_z_y_x[i] = -2.0 * g_x_0_y_x[i] * a_exp + 4.0 * g_x_zz_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_y_y[i] = -2.0 * g_x_0_y_y[i] * a_exp + 4.0 * g_x_zz_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_y_z[i] = -2.0 * g_x_0_y_z[i] * a_exp + 4.0 * g_x_zz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_0_z_x, g_x_0_z_y, g_x_0_z_z, g_x_z_0_0_0_z_z_x, g_x_z_0_0_0_z_z_y, g_x_z_0_0_0_z_z_z, g_x_zz_z_x, g_x_zz_z_y, g_x_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_z_z_x[i] = -2.0 * g_x_0_z_x[i] * a_exp + 4.0 * g_x_zz_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_z_y[i] = -2.0 * g_x_0_z_y[i] * a_exp + 4.0 * g_x_zz_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_z_z[i] = -2.0 * g_x_0_z_z[i] * a_exp + 4.0 * g_x_zz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_y_0_x_x, g_y_0_x_y, g_y_0_x_z, g_y_x_0_0_0_x_x_x, g_y_x_0_0_0_x_x_y, g_y_x_0_0_0_x_x_z, g_y_xx_x_x, g_y_xx_x_y, g_y_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_x_x_x[i] = -2.0 * g_y_0_x_x[i] * a_exp + 4.0 * g_y_xx_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_x_y[i] = -2.0 * g_y_0_x_y[i] * a_exp + 4.0 * g_y_xx_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_x_z[i] = -2.0 * g_y_0_x_z[i] * a_exp + 4.0 * g_y_xx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_y_0_y_x, g_y_0_y_y, g_y_0_y_z, g_y_x_0_0_0_x_y_x, g_y_x_0_0_0_x_y_y, g_y_x_0_0_0_x_y_z, g_y_xx_y_x, g_y_xx_y_y, g_y_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_x_y_x[i] = -2.0 * g_y_0_y_x[i] * a_exp + 4.0 * g_y_xx_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_y_y[i] = -2.0 * g_y_0_y_y[i] * a_exp + 4.0 * g_y_xx_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_y_z[i] = -2.0 * g_y_0_y_z[i] * a_exp + 4.0 * g_y_xx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_y_0_z_x, g_y_0_z_y, g_y_0_z_z, g_y_x_0_0_0_x_z_x, g_y_x_0_0_0_x_z_y, g_y_x_0_0_0_x_z_z, g_y_xx_z_x, g_y_xx_z_y, g_y_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_x_z_x[i] = -2.0 * g_y_0_z_x[i] * a_exp + 4.0 * g_y_xx_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_z_y[i] = -2.0 * g_y_0_z_y[i] * a_exp + 4.0 * g_y_xx_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_z_z[i] = -2.0 * g_y_0_z_z[i] * a_exp + 4.0 * g_y_xx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_y_x_0_0_0_y_x_x, g_y_x_0_0_0_y_x_y, g_y_x_0_0_0_y_x_z, g_y_xy_x_x, g_y_xy_x_y, g_y_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_y_x_x[i] = 4.0 * g_y_xy_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_x_y[i] = 4.0 * g_y_xy_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_x_z[i] = 4.0 * g_y_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_y_x_0_0_0_y_y_x, g_y_x_0_0_0_y_y_y, g_y_x_0_0_0_y_y_z, g_y_xy_y_x, g_y_xy_y_y, g_y_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_y_y_x[i] = 4.0 * g_y_xy_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_y_y[i] = 4.0 * g_y_xy_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_y_z[i] = 4.0 * g_y_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_y_x_0_0_0_y_z_x, g_y_x_0_0_0_y_z_y, g_y_x_0_0_0_y_z_z, g_y_xy_z_x, g_y_xy_z_y, g_y_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_y_z_x[i] = 4.0 * g_y_xy_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_z_y[i] = 4.0 * g_y_xy_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_z_z[i] = 4.0 * g_y_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_y_x_0_0_0_z_x_x, g_y_x_0_0_0_z_x_y, g_y_x_0_0_0_z_x_z, g_y_xz_x_x, g_y_xz_x_y, g_y_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_z_x_x[i] = 4.0 * g_y_xz_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_x_y[i] = 4.0 * g_y_xz_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_x_z[i] = 4.0 * g_y_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_y_x_0_0_0_z_y_x, g_y_x_0_0_0_z_y_y, g_y_x_0_0_0_z_y_z, g_y_xz_y_x, g_y_xz_y_y, g_y_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_z_y_x[i] = 4.0 * g_y_xz_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_y_y[i] = 4.0 * g_y_xz_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_y_z[i] = 4.0 * g_y_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_y_x_0_0_0_z_z_x, g_y_x_0_0_0_z_z_y, g_y_x_0_0_0_z_z_z, g_y_xz_z_x, g_y_xz_z_y, g_y_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_z_z_x[i] = 4.0 * g_y_xz_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_z_y[i] = 4.0 * g_y_xz_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_z_z[i] = 4.0 * g_y_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_y_xy_x_x, g_y_xy_x_y, g_y_xy_x_z, g_y_y_0_0_0_x_x_x, g_y_y_0_0_0_x_x_y, g_y_y_0_0_0_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_x_x_x[i] = 4.0 * g_y_xy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_x_y[i] = 4.0 * g_y_xy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_x_z[i] = 4.0 * g_y_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_y_xy_y_x, g_y_xy_y_y, g_y_xy_y_z, g_y_y_0_0_0_x_y_x, g_y_y_0_0_0_x_y_y, g_y_y_0_0_0_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_x_y_x[i] = 4.0 * g_y_xy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_y_y[i] = 4.0 * g_y_xy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_y_z[i] = 4.0 * g_y_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_y_xy_z_x, g_y_xy_z_y, g_y_xy_z_z, g_y_y_0_0_0_x_z_x, g_y_y_0_0_0_x_z_y, g_y_y_0_0_0_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_x_z_x[i] = 4.0 * g_y_xy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_z_y[i] = 4.0 * g_y_xy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_z_z[i] = 4.0 * g_y_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_y_0_x_x, g_y_0_x_y, g_y_0_x_z, g_y_y_0_0_0_y_x_x, g_y_y_0_0_0_y_x_y, g_y_y_0_0_0_y_x_z, g_y_yy_x_x, g_y_yy_x_y, g_y_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_y_x_x[i] = -2.0 * g_y_0_x_x[i] * a_exp + 4.0 * g_y_yy_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_x_y[i] = -2.0 * g_y_0_x_y[i] * a_exp + 4.0 * g_y_yy_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_x_z[i] = -2.0 * g_y_0_x_z[i] * a_exp + 4.0 * g_y_yy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_y_0_y_x, g_y_0_y_y, g_y_0_y_z, g_y_y_0_0_0_y_y_x, g_y_y_0_0_0_y_y_y, g_y_y_0_0_0_y_y_z, g_y_yy_y_x, g_y_yy_y_y, g_y_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_y_y_x[i] = -2.0 * g_y_0_y_x[i] * a_exp + 4.0 * g_y_yy_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_y_y[i] = -2.0 * g_y_0_y_y[i] * a_exp + 4.0 * g_y_yy_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_y_z[i] = -2.0 * g_y_0_y_z[i] * a_exp + 4.0 * g_y_yy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_y_0_z_x, g_y_0_z_y, g_y_0_z_z, g_y_y_0_0_0_y_z_x, g_y_y_0_0_0_y_z_y, g_y_y_0_0_0_y_z_z, g_y_yy_z_x, g_y_yy_z_y, g_y_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_y_z_x[i] = -2.0 * g_y_0_z_x[i] * a_exp + 4.0 * g_y_yy_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_z_y[i] = -2.0 * g_y_0_z_y[i] * a_exp + 4.0 * g_y_yy_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_z_z[i] = -2.0 * g_y_0_z_z[i] * a_exp + 4.0 * g_y_yy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_y_y_0_0_0_z_x_x, g_y_y_0_0_0_z_x_y, g_y_y_0_0_0_z_x_z, g_y_yz_x_x, g_y_yz_x_y, g_y_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_z_x_x[i] = 4.0 * g_y_yz_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_x_y[i] = 4.0 * g_y_yz_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_x_z[i] = 4.0 * g_y_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_y_y_0_0_0_z_y_x, g_y_y_0_0_0_z_y_y, g_y_y_0_0_0_z_y_z, g_y_yz_y_x, g_y_yz_y_y, g_y_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_z_y_x[i] = 4.0 * g_y_yz_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_y_y[i] = 4.0 * g_y_yz_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_y_z[i] = 4.0 * g_y_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_y_y_0_0_0_z_z_x, g_y_y_0_0_0_z_z_y, g_y_y_0_0_0_z_z_z, g_y_yz_z_x, g_y_yz_z_y, g_y_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_z_z_x[i] = 4.0 * g_y_yz_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_z_y[i] = 4.0 * g_y_yz_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_z_z[i] = 4.0 * g_y_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_y_xz_x_x, g_y_xz_x_y, g_y_xz_x_z, g_y_z_0_0_0_x_x_x, g_y_z_0_0_0_x_x_y, g_y_z_0_0_0_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_x_x_x[i] = 4.0 * g_y_xz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_x_y[i] = 4.0 * g_y_xz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_x_z[i] = 4.0 * g_y_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_y_xz_y_x, g_y_xz_y_y, g_y_xz_y_z, g_y_z_0_0_0_x_y_x, g_y_z_0_0_0_x_y_y, g_y_z_0_0_0_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_x_y_x[i] = 4.0 * g_y_xz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_y_y[i] = 4.0 * g_y_xz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_y_z[i] = 4.0 * g_y_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_y_xz_z_x, g_y_xz_z_y, g_y_xz_z_z, g_y_z_0_0_0_x_z_x, g_y_z_0_0_0_x_z_y, g_y_z_0_0_0_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_x_z_x[i] = 4.0 * g_y_xz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_z_y[i] = 4.0 * g_y_xz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_z_z[i] = 4.0 * g_y_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_y_yz_x_x, g_y_yz_x_y, g_y_yz_x_z, g_y_z_0_0_0_y_x_x, g_y_z_0_0_0_y_x_y, g_y_z_0_0_0_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_y_x_x[i] = 4.0 * g_y_yz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_x_y[i] = 4.0 * g_y_yz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_x_z[i] = 4.0 * g_y_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_y_yz_y_x, g_y_yz_y_y, g_y_yz_y_z, g_y_z_0_0_0_y_y_x, g_y_z_0_0_0_y_y_y, g_y_z_0_0_0_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_y_y_x[i] = 4.0 * g_y_yz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_y_y[i] = 4.0 * g_y_yz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_y_z[i] = 4.0 * g_y_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_y_yz_z_x, g_y_yz_z_y, g_y_yz_z_z, g_y_z_0_0_0_y_z_x, g_y_z_0_0_0_y_z_y, g_y_z_0_0_0_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_y_z_x[i] = 4.0 * g_y_yz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_z_y[i] = 4.0 * g_y_yz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_z_z[i] = 4.0 * g_y_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_y_0_x_x, g_y_0_x_y, g_y_0_x_z, g_y_z_0_0_0_z_x_x, g_y_z_0_0_0_z_x_y, g_y_z_0_0_0_z_x_z, g_y_zz_x_x, g_y_zz_x_y, g_y_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_z_x_x[i] = -2.0 * g_y_0_x_x[i] * a_exp + 4.0 * g_y_zz_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_x_y[i] = -2.0 * g_y_0_x_y[i] * a_exp + 4.0 * g_y_zz_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_x_z[i] = -2.0 * g_y_0_x_z[i] * a_exp + 4.0 * g_y_zz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_y_0_y_x, g_y_0_y_y, g_y_0_y_z, g_y_z_0_0_0_z_y_x, g_y_z_0_0_0_z_y_y, g_y_z_0_0_0_z_y_z, g_y_zz_y_x, g_y_zz_y_y, g_y_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_z_y_x[i] = -2.0 * g_y_0_y_x[i] * a_exp + 4.0 * g_y_zz_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_y_y[i] = -2.0 * g_y_0_y_y[i] * a_exp + 4.0 * g_y_zz_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_y_z[i] = -2.0 * g_y_0_y_z[i] * a_exp + 4.0 * g_y_zz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_y_0_z_x, g_y_0_z_y, g_y_0_z_z, g_y_z_0_0_0_z_z_x, g_y_z_0_0_0_z_z_y, g_y_z_0_0_0_z_z_z, g_y_zz_z_x, g_y_zz_z_y, g_y_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_z_z_x[i] = -2.0 * g_y_0_z_x[i] * a_exp + 4.0 * g_y_zz_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_z_y[i] = -2.0 * g_y_0_z_y[i] * a_exp + 4.0 * g_y_zz_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_z_z[i] = -2.0 * g_y_0_z_z[i] * a_exp + 4.0 * g_y_zz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_z_0_x_x, g_z_0_x_y, g_z_0_x_z, g_z_x_0_0_0_x_x_x, g_z_x_0_0_0_x_x_y, g_z_x_0_0_0_x_x_z, g_z_xx_x_x, g_z_xx_x_y, g_z_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_x_x_x[i] = -2.0 * g_z_0_x_x[i] * a_exp + 4.0 * g_z_xx_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_x_y[i] = -2.0 * g_z_0_x_y[i] * a_exp + 4.0 * g_z_xx_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_x_z[i] = -2.0 * g_z_0_x_z[i] * a_exp + 4.0 * g_z_xx_x_z[i] * a_exp * b_exp;
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_z_0_y_x, g_z_0_y_y, g_z_0_y_z, g_z_x_0_0_0_x_y_x, g_z_x_0_0_0_x_y_y, g_z_x_0_0_0_x_y_z, g_z_xx_y_x, g_z_xx_y_y, g_z_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_x_y_x[i] = -2.0 * g_z_0_y_x[i] * a_exp + 4.0 * g_z_xx_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_y_y[i] = -2.0 * g_z_0_y_y[i] * a_exp + 4.0 * g_z_xx_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_y_z[i] = -2.0 * g_z_0_y_z[i] * a_exp + 4.0 * g_z_xx_y_z[i] * a_exp * b_exp;
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_z_0_z_x, g_z_0_z_y, g_z_0_z_z, g_z_x_0_0_0_x_z_x, g_z_x_0_0_0_x_z_y, g_z_x_0_0_0_x_z_z, g_z_xx_z_x, g_z_xx_z_y, g_z_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_x_z_x[i] = -2.0 * g_z_0_z_x[i] * a_exp + 4.0 * g_z_xx_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_z_y[i] = -2.0 * g_z_0_z_y[i] * a_exp + 4.0 * g_z_xx_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_z_z[i] = -2.0 * g_z_0_z_z[i] * a_exp + 4.0 * g_z_xx_z_z[i] * a_exp * b_exp;
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_z_x_0_0_0_y_x_x, g_z_x_0_0_0_y_x_y, g_z_x_0_0_0_y_x_z, g_z_xy_x_x, g_z_xy_x_y, g_z_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_y_x_x[i] = 4.0 * g_z_xy_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_x_y[i] = 4.0 * g_z_xy_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_x_z[i] = 4.0 * g_z_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_z_x_0_0_0_y_y_x, g_z_x_0_0_0_y_y_y, g_z_x_0_0_0_y_y_z, g_z_xy_y_x, g_z_xy_y_y, g_z_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_y_y_x[i] = 4.0 * g_z_xy_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_y_y[i] = 4.0 * g_z_xy_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_y_z[i] = 4.0 * g_z_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_z_x_0_0_0_y_z_x, g_z_x_0_0_0_y_z_y, g_z_x_0_0_0_y_z_z, g_z_xy_z_x, g_z_xy_z_y, g_z_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_y_z_x[i] = 4.0 * g_z_xy_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_z_y[i] = 4.0 * g_z_xy_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_z_z[i] = 4.0 * g_z_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_z_x_0_0_0_z_x_x, g_z_x_0_0_0_z_x_y, g_z_x_0_0_0_z_x_z, g_z_xz_x_x, g_z_xz_x_y, g_z_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_z_x_x[i] = 4.0 * g_z_xz_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_x_y[i] = 4.0 * g_z_xz_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_x_z[i] = 4.0 * g_z_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_z_x_0_0_0_z_y_x, g_z_x_0_0_0_z_y_y, g_z_x_0_0_0_z_y_z, g_z_xz_y_x, g_z_xz_y_y, g_z_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_z_y_x[i] = 4.0 * g_z_xz_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_y_y[i] = 4.0 * g_z_xz_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_y_z[i] = 4.0 * g_z_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_z_x_0_0_0_z_z_x, g_z_x_0_0_0_z_z_y, g_z_x_0_0_0_z_z_z, g_z_xz_z_x, g_z_xz_z_y, g_z_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_z_z_x[i] = 4.0 * g_z_xz_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_z_y[i] = 4.0 * g_z_xz_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_z_z[i] = 4.0 * g_z_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_z_xy_x_x, g_z_xy_x_y, g_z_xy_x_z, g_z_y_0_0_0_x_x_x, g_z_y_0_0_0_x_x_y, g_z_y_0_0_0_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_x_x_x[i] = 4.0 * g_z_xy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_x_y[i] = 4.0 * g_z_xy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_x_z[i] = 4.0 * g_z_xy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_z_xy_y_x, g_z_xy_y_y, g_z_xy_y_z, g_z_y_0_0_0_x_y_x, g_z_y_0_0_0_x_y_y, g_z_y_0_0_0_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_x_y_x[i] = 4.0 * g_z_xy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_y_y[i] = 4.0 * g_z_xy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_y_z[i] = 4.0 * g_z_xy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_z_xy_z_x, g_z_xy_z_y, g_z_xy_z_z, g_z_y_0_0_0_x_z_x, g_z_y_0_0_0_x_z_y, g_z_y_0_0_0_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_x_z_x[i] = 4.0 * g_z_xy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_z_y[i] = 4.0 * g_z_xy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_z_z[i] = 4.0 * g_z_xy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_z_0_x_x, g_z_0_x_y, g_z_0_x_z, g_z_y_0_0_0_y_x_x, g_z_y_0_0_0_y_x_y, g_z_y_0_0_0_y_x_z, g_z_yy_x_x, g_z_yy_x_y, g_z_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_y_x_x[i] = -2.0 * g_z_0_x_x[i] * a_exp + 4.0 * g_z_yy_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_x_y[i] = -2.0 * g_z_0_x_y[i] * a_exp + 4.0 * g_z_yy_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_x_z[i] = -2.0 * g_z_0_x_z[i] * a_exp + 4.0 * g_z_yy_x_z[i] * a_exp * b_exp;
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_z_0_y_x, g_z_0_y_y, g_z_0_y_z, g_z_y_0_0_0_y_y_x, g_z_y_0_0_0_y_y_y, g_z_y_0_0_0_y_y_z, g_z_yy_y_x, g_z_yy_y_y, g_z_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_y_y_x[i] = -2.0 * g_z_0_y_x[i] * a_exp + 4.0 * g_z_yy_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_y_y[i] = -2.0 * g_z_0_y_y[i] * a_exp + 4.0 * g_z_yy_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_y_z[i] = -2.0 * g_z_0_y_z[i] * a_exp + 4.0 * g_z_yy_y_z[i] * a_exp * b_exp;
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_z_0_z_x, g_z_0_z_y, g_z_0_z_z, g_z_y_0_0_0_y_z_x, g_z_y_0_0_0_y_z_y, g_z_y_0_0_0_y_z_z, g_z_yy_z_x, g_z_yy_z_y, g_z_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_y_z_x[i] = -2.0 * g_z_0_z_x[i] * a_exp + 4.0 * g_z_yy_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_z_y[i] = -2.0 * g_z_0_z_y[i] * a_exp + 4.0 * g_z_yy_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_z_z[i] = -2.0 * g_z_0_z_z[i] * a_exp + 4.0 * g_z_yy_z_z[i] * a_exp * b_exp;
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_z_y_0_0_0_z_x_x, g_z_y_0_0_0_z_x_y, g_z_y_0_0_0_z_x_z, g_z_yz_x_x, g_z_yz_x_y, g_z_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_z_x_x[i] = 4.0 * g_z_yz_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_x_y[i] = 4.0 * g_z_yz_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_x_z[i] = 4.0 * g_z_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_z_y_0_0_0_z_y_x, g_z_y_0_0_0_z_y_y, g_z_y_0_0_0_z_y_z, g_z_yz_y_x, g_z_yz_y_y, g_z_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_z_y_x[i] = 4.0 * g_z_yz_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_y_y[i] = 4.0 * g_z_yz_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_y_z[i] = 4.0 * g_z_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_z_y_0_0_0_z_z_x, g_z_y_0_0_0_z_z_y, g_z_y_0_0_0_z_z_z, g_z_yz_z_x, g_z_yz_z_y, g_z_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_z_z_x[i] = 4.0 * g_z_yz_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_z_y[i] = 4.0 * g_z_yz_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_z_z[i] = 4.0 * g_z_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_z_xz_x_x, g_z_xz_x_y, g_z_xz_x_z, g_z_z_0_0_0_x_x_x, g_z_z_0_0_0_x_x_y, g_z_z_0_0_0_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_x_x_x[i] = 4.0 * g_z_xz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_x_y[i] = 4.0 * g_z_xz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_x_z[i] = 4.0 * g_z_xz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_z_xz_y_x, g_z_xz_y_y, g_z_xz_y_z, g_z_z_0_0_0_x_y_x, g_z_z_0_0_0_x_y_y, g_z_z_0_0_0_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_x_y_x[i] = 4.0 * g_z_xz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_y_y[i] = 4.0 * g_z_xz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_y_z[i] = 4.0 * g_z_xz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_z_xz_z_x, g_z_xz_z_y, g_z_xz_z_z, g_z_z_0_0_0_x_z_x, g_z_z_0_0_0_x_z_y, g_z_z_0_0_0_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_x_z_x[i] = 4.0 * g_z_xz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_z_y[i] = 4.0 * g_z_xz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_z_z[i] = 4.0 * g_z_xz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_z_yz_x_x, g_z_yz_x_y, g_z_yz_x_z, g_z_z_0_0_0_y_x_x, g_z_z_0_0_0_y_x_y, g_z_z_0_0_0_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_y_x_x[i] = 4.0 * g_z_yz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_x_y[i] = 4.0 * g_z_yz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_x_z[i] = 4.0 * g_z_yz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_z_yz_y_x, g_z_yz_y_y, g_z_yz_y_z, g_z_z_0_0_0_y_y_x, g_z_z_0_0_0_y_y_y, g_z_z_0_0_0_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_y_y_x[i] = 4.0 * g_z_yz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_y_y[i] = 4.0 * g_z_yz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_y_z[i] = 4.0 * g_z_yz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_z_yz_z_x, g_z_yz_z_y, g_z_yz_z_z, g_z_z_0_0_0_y_z_x, g_z_z_0_0_0_y_z_y, g_z_z_0_0_0_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_y_z_x[i] = 4.0 * g_z_yz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_z_y[i] = 4.0 * g_z_yz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_z_z[i] = 4.0 * g_z_yz_z_z[i] * a_exp * b_exp;
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_z_0_x_x, g_z_0_x_y, g_z_0_x_z, g_z_z_0_0_0_z_x_x, g_z_z_0_0_0_z_x_y, g_z_z_0_0_0_z_x_z, g_z_zz_x_x, g_z_zz_x_y, g_z_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_z_x_x[i] = -2.0 * g_z_0_x_x[i] * a_exp + 4.0 * g_z_zz_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_x_y[i] = -2.0 * g_z_0_x_y[i] * a_exp + 4.0 * g_z_zz_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_x_z[i] = -2.0 * g_z_0_x_z[i] * a_exp + 4.0 * g_z_zz_x_z[i] * a_exp * b_exp;
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_z_0_y_x, g_z_0_y_y, g_z_0_y_z, g_z_z_0_0_0_z_y_x, g_z_z_0_0_0_z_y_y, g_z_z_0_0_0_z_y_z, g_z_zz_y_x, g_z_zz_y_y, g_z_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_z_y_x[i] = -2.0 * g_z_0_y_x[i] * a_exp + 4.0 * g_z_zz_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_y_y[i] = -2.0 * g_z_0_y_y[i] * a_exp + 4.0 * g_z_zz_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_y_z[i] = -2.0 * g_z_0_y_z[i] * a_exp + 4.0 * g_z_zz_y_z[i] * a_exp * b_exp;
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_z_0_z_x, g_z_0_z_y, g_z_0_z_z, g_z_z_0_0_0_z_z_x, g_z_z_0_0_0_z_z_y, g_z_z_0_0_0_z_z_z, g_z_zz_z_x, g_z_zz_z_y, g_z_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_z_z_x[i] = -2.0 * g_z_0_z_x[i] * a_exp + 4.0 * g_z_zz_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_z_y[i] = -2.0 * g_z_0_z_y[i] * a_exp + 4.0 * g_z_zz_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_z_z[i] = -2.0 * g_z_0_z_z[i] * a_exp + 4.0 * g_z_zz_z_z[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

