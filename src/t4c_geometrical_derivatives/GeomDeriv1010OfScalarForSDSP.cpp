#include "GeomDeriv1010OfScalarForSDSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_sdsp_0(CSimdArray<double>& buffer_1010_sdsp,
                     const CSimdArray<double>& buffer_pdpp,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_sdsp.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1010_sdsp

    auto g_x_0_x_0_0_xx_0_x = buffer_1010_sdsp[0];

    auto g_x_0_x_0_0_xx_0_y = buffer_1010_sdsp[1];

    auto g_x_0_x_0_0_xx_0_z = buffer_1010_sdsp[2];

    auto g_x_0_x_0_0_xy_0_x = buffer_1010_sdsp[3];

    auto g_x_0_x_0_0_xy_0_y = buffer_1010_sdsp[4];

    auto g_x_0_x_0_0_xy_0_z = buffer_1010_sdsp[5];

    auto g_x_0_x_0_0_xz_0_x = buffer_1010_sdsp[6];

    auto g_x_0_x_0_0_xz_0_y = buffer_1010_sdsp[7];

    auto g_x_0_x_0_0_xz_0_z = buffer_1010_sdsp[8];

    auto g_x_0_x_0_0_yy_0_x = buffer_1010_sdsp[9];

    auto g_x_0_x_0_0_yy_0_y = buffer_1010_sdsp[10];

    auto g_x_0_x_0_0_yy_0_z = buffer_1010_sdsp[11];

    auto g_x_0_x_0_0_yz_0_x = buffer_1010_sdsp[12];

    auto g_x_0_x_0_0_yz_0_y = buffer_1010_sdsp[13];

    auto g_x_0_x_0_0_yz_0_z = buffer_1010_sdsp[14];

    auto g_x_0_x_0_0_zz_0_x = buffer_1010_sdsp[15];

    auto g_x_0_x_0_0_zz_0_y = buffer_1010_sdsp[16];

    auto g_x_0_x_0_0_zz_0_z = buffer_1010_sdsp[17];

    auto g_x_0_y_0_0_xx_0_x = buffer_1010_sdsp[18];

    auto g_x_0_y_0_0_xx_0_y = buffer_1010_sdsp[19];

    auto g_x_0_y_0_0_xx_0_z = buffer_1010_sdsp[20];

    auto g_x_0_y_0_0_xy_0_x = buffer_1010_sdsp[21];

    auto g_x_0_y_0_0_xy_0_y = buffer_1010_sdsp[22];

    auto g_x_0_y_0_0_xy_0_z = buffer_1010_sdsp[23];

    auto g_x_0_y_0_0_xz_0_x = buffer_1010_sdsp[24];

    auto g_x_0_y_0_0_xz_0_y = buffer_1010_sdsp[25];

    auto g_x_0_y_0_0_xz_0_z = buffer_1010_sdsp[26];

    auto g_x_0_y_0_0_yy_0_x = buffer_1010_sdsp[27];

    auto g_x_0_y_0_0_yy_0_y = buffer_1010_sdsp[28];

    auto g_x_0_y_0_0_yy_0_z = buffer_1010_sdsp[29];

    auto g_x_0_y_0_0_yz_0_x = buffer_1010_sdsp[30];

    auto g_x_0_y_0_0_yz_0_y = buffer_1010_sdsp[31];

    auto g_x_0_y_0_0_yz_0_z = buffer_1010_sdsp[32];

    auto g_x_0_y_0_0_zz_0_x = buffer_1010_sdsp[33];

    auto g_x_0_y_0_0_zz_0_y = buffer_1010_sdsp[34];

    auto g_x_0_y_0_0_zz_0_z = buffer_1010_sdsp[35];

    auto g_x_0_z_0_0_xx_0_x = buffer_1010_sdsp[36];

    auto g_x_0_z_0_0_xx_0_y = buffer_1010_sdsp[37];

    auto g_x_0_z_0_0_xx_0_z = buffer_1010_sdsp[38];

    auto g_x_0_z_0_0_xy_0_x = buffer_1010_sdsp[39];

    auto g_x_0_z_0_0_xy_0_y = buffer_1010_sdsp[40];

    auto g_x_0_z_0_0_xy_0_z = buffer_1010_sdsp[41];

    auto g_x_0_z_0_0_xz_0_x = buffer_1010_sdsp[42];

    auto g_x_0_z_0_0_xz_0_y = buffer_1010_sdsp[43];

    auto g_x_0_z_0_0_xz_0_z = buffer_1010_sdsp[44];

    auto g_x_0_z_0_0_yy_0_x = buffer_1010_sdsp[45];

    auto g_x_0_z_0_0_yy_0_y = buffer_1010_sdsp[46];

    auto g_x_0_z_0_0_yy_0_z = buffer_1010_sdsp[47];

    auto g_x_0_z_0_0_yz_0_x = buffer_1010_sdsp[48];

    auto g_x_0_z_0_0_yz_0_y = buffer_1010_sdsp[49];

    auto g_x_0_z_0_0_yz_0_z = buffer_1010_sdsp[50];

    auto g_x_0_z_0_0_zz_0_x = buffer_1010_sdsp[51];

    auto g_x_0_z_0_0_zz_0_y = buffer_1010_sdsp[52];

    auto g_x_0_z_0_0_zz_0_z = buffer_1010_sdsp[53];

    auto g_y_0_x_0_0_xx_0_x = buffer_1010_sdsp[54];

    auto g_y_0_x_0_0_xx_0_y = buffer_1010_sdsp[55];

    auto g_y_0_x_0_0_xx_0_z = buffer_1010_sdsp[56];

    auto g_y_0_x_0_0_xy_0_x = buffer_1010_sdsp[57];

    auto g_y_0_x_0_0_xy_0_y = buffer_1010_sdsp[58];

    auto g_y_0_x_0_0_xy_0_z = buffer_1010_sdsp[59];

    auto g_y_0_x_0_0_xz_0_x = buffer_1010_sdsp[60];

    auto g_y_0_x_0_0_xz_0_y = buffer_1010_sdsp[61];

    auto g_y_0_x_0_0_xz_0_z = buffer_1010_sdsp[62];

    auto g_y_0_x_0_0_yy_0_x = buffer_1010_sdsp[63];

    auto g_y_0_x_0_0_yy_0_y = buffer_1010_sdsp[64];

    auto g_y_0_x_0_0_yy_0_z = buffer_1010_sdsp[65];

    auto g_y_0_x_0_0_yz_0_x = buffer_1010_sdsp[66];

    auto g_y_0_x_0_0_yz_0_y = buffer_1010_sdsp[67];

    auto g_y_0_x_0_0_yz_0_z = buffer_1010_sdsp[68];

    auto g_y_0_x_0_0_zz_0_x = buffer_1010_sdsp[69];

    auto g_y_0_x_0_0_zz_0_y = buffer_1010_sdsp[70];

    auto g_y_0_x_0_0_zz_0_z = buffer_1010_sdsp[71];

    auto g_y_0_y_0_0_xx_0_x = buffer_1010_sdsp[72];

    auto g_y_0_y_0_0_xx_0_y = buffer_1010_sdsp[73];

    auto g_y_0_y_0_0_xx_0_z = buffer_1010_sdsp[74];

    auto g_y_0_y_0_0_xy_0_x = buffer_1010_sdsp[75];

    auto g_y_0_y_0_0_xy_0_y = buffer_1010_sdsp[76];

    auto g_y_0_y_0_0_xy_0_z = buffer_1010_sdsp[77];

    auto g_y_0_y_0_0_xz_0_x = buffer_1010_sdsp[78];

    auto g_y_0_y_0_0_xz_0_y = buffer_1010_sdsp[79];

    auto g_y_0_y_0_0_xz_0_z = buffer_1010_sdsp[80];

    auto g_y_0_y_0_0_yy_0_x = buffer_1010_sdsp[81];

    auto g_y_0_y_0_0_yy_0_y = buffer_1010_sdsp[82];

    auto g_y_0_y_0_0_yy_0_z = buffer_1010_sdsp[83];

    auto g_y_0_y_0_0_yz_0_x = buffer_1010_sdsp[84];

    auto g_y_0_y_0_0_yz_0_y = buffer_1010_sdsp[85];

    auto g_y_0_y_0_0_yz_0_z = buffer_1010_sdsp[86];

    auto g_y_0_y_0_0_zz_0_x = buffer_1010_sdsp[87];

    auto g_y_0_y_0_0_zz_0_y = buffer_1010_sdsp[88];

    auto g_y_0_y_0_0_zz_0_z = buffer_1010_sdsp[89];

    auto g_y_0_z_0_0_xx_0_x = buffer_1010_sdsp[90];

    auto g_y_0_z_0_0_xx_0_y = buffer_1010_sdsp[91];

    auto g_y_0_z_0_0_xx_0_z = buffer_1010_sdsp[92];

    auto g_y_0_z_0_0_xy_0_x = buffer_1010_sdsp[93];

    auto g_y_0_z_0_0_xy_0_y = buffer_1010_sdsp[94];

    auto g_y_0_z_0_0_xy_0_z = buffer_1010_sdsp[95];

    auto g_y_0_z_0_0_xz_0_x = buffer_1010_sdsp[96];

    auto g_y_0_z_0_0_xz_0_y = buffer_1010_sdsp[97];

    auto g_y_0_z_0_0_xz_0_z = buffer_1010_sdsp[98];

    auto g_y_0_z_0_0_yy_0_x = buffer_1010_sdsp[99];

    auto g_y_0_z_0_0_yy_0_y = buffer_1010_sdsp[100];

    auto g_y_0_z_0_0_yy_0_z = buffer_1010_sdsp[101];

    auto g_y_0_z_0_0_yz_0_x = buffer_1010_sdsp[102];

    auto g_y_0_z_0_0_yz_0_y = buffer_1010_sdsp[103];

    auto g_y_0_z_0_0_yz_0_z = buffer_1010_sdsp[104];

    auto g_y_0_z_0_0_zz_0_x = buffer_1010_sdsp[105];

    auto g_y_0_z_0_0_zz_0_y = buffer_1010_sdsp[106];

    auto g_y_0_z_0_0_zz_0_z = buffer_1010_sdsp[107];

    auto g_z_0_x_0_0_xx_0_x = buffer_1010_sdsp[108];

    auto g_z_0_x_0_0_xx_0_y = buffer_1010_sdsp[109];

    auto g_z_0_x_0_0_xx_0_z = buffer_1010_sdsp[110];

    auto g_z_0_x_0_0_xy_0_x = buffer_1010_sdsp[111];

    auto g_z_0_x_0_0_xy_0_y = buffer_1010_sdsp[112];

    auto g_z_0_x_0_0_xy_0_z = buffer_1010_sdsp[113];

    auto g_z_0_x_0_0_xz_0_x = buffer_1010_sdsp[114];

    auto g_z_0_x_0_0_xz_0_y = buffer_1010_sdsp[115];

    auto g_z_0_x_0_0_xz_0_z = buffer_1010_sdsp[116];

    auto g_z_0_x_0_0_yy_0_x = buffer_1010_sdsp[117];

    auto g_z_0_x_0_0_yy_0_y = buffer_1010_sdsp[118];

    auto g_z_0_x_0_0_yy_0_z = buffer_1010_sdsp[119];

    auto g_z_0_x_0_0_yz_0_x = buffer_1010_sdsp[120];

    auto g_z_0_x_0_0_yz_0_y = buffer_1010_sdsp[121];

    auto g_z_0_x_0_0_yz_0_z = buffer_1010_sdsp[122];

    auto g_z_0_x_0_0_zz_0_x = buffer_1010_sdsp[123];

    auto g_z_0_x_0_0_zz_0_y = buffer_1010_sdsp[124];

    auto g_z_0_x_0_0_zz_0_z = buffer_1010_sdsp[125];

    auto g_z_0_y_0_0_xx_0_x = buffer_1010_sdsp[126];

    auto g_z_0_y_0_0_xx_0_y = buffer_1010_sdsp[127];

    auto g_z_0_y_0_0_xx_0_z = buffer_1010_sdsp[128];

    auto g_z_0_y_0_0_xy_0_x = buffer_1010_sdsp[129];

    auto g_z_0_y_0_0_xy_0_y = buffer_1010_sdsp[130];

    auto g_z_0_y_0_0_xy_0_z = buffer_1010_sdsp[131];

    auto g_z_0_y_0_0_xz_0_x = buffer_1010_sdsp[132];

    auto g_z_0_y_0_0_xz_0_y = buffer_1010_sdsp[133];

    auto g_z_0_y_0_0_xz_0_z = buffer_1010_sdsp[134];

    auto g_z_0_y_0_0_yy_0_x = buffer_1010_sdsp[135];

    auto g_z_0_y_0_0_yy_0_y = buffer_1010_sdsp[136];

    auto g_z_0_y_0_0_yy_0_z = buffer_1010_sdsp[137];

    auto g_z_0_y_0_0_yz_0_x = buffer_1010_sdsp[138];

    auto g_z_0_y_0_0_yz_0_y = buffer_1010_sdsp[139];

    auto g_z_0_y_0_0_yz_0_z = buffer_1010_sdsp[140];

    auto g_z_0_y_0_0_zz_0_x = buffer_1010_sdsp[141];

    auto g_z_0_y_0_0_zz_0_y = buffer_1010_sdsp[142];

    auto g_z_0_y_0_0_zz_0_z = buffer_1010_sdsp[143];

    auto g_z_0_z_0_0_xx_0_x = buffer_1010_sdsp[144];

    auto g_z_0_z_0_0_xx_0_y = buffer_1010_sdsp[145];

    auto g_z_0_z_0_0_xx_0_z = buffer_1010_sdsp[146];

    auto g_z_0_z_0_0_xy_0_x = buffer_1010_sdsp[147];

    auto g_z_0_z_0_0_xy_0_y = buffer_1010_sdsp[148];

    auto g_z_0_z_0_0_xy_0_z = buffer_1010_sdsp[149];

    auto g_z_0_z_0_0_xz_0_x = buffer_1010_sdsp[150];

    auto g_z_0_z_0_0_xz_0_y = buffer_1010_sdsp[151];

    auto g_z_0_z_0_0_xz_0_z = buffer_1010_sdsp[152];

    auto g_z_0_z_0_0_yy_0_x = buffer_1010_sdsp[153];

    auto g_z_0_z_0_0_yy_0_y = buffer_1010_sdsp[154];

    auto g_z_0_z_0_0_yy_0_z = buffer_1010_sdsp[155];

    auto g_z_0_z_0_0_yz_0_x = buffer_1010_sdsp[156];

    auto g_z_0_z_0_0_yz_0_y = buffer_1010_sdsp[157];

    auto g_z_0_z_0_0_yz_0_z = buffer_1010_sdsp[158];

    auto g_z_0_z_0_0_zz_0_x = buffer_1010_sdsp[159];

    auto g_z_0_z_0_0_zz_0_y = buffer_1010_sdsp[160];

    auto g_z_0_z_0_0_zz_0_z = buffer_1010_sdsp[161];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_x_0_0_xx_0_x, g_x_0_x_0_0_xx_0_y, g_x_0_x_0_0_xx_0_z, g_x_xx_x_x, g_x_xx_x_y, g_x_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xx_0_x[i] = 4.0 * g_x_xx_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_0_y[i] = 4.0 * g_x_xx_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xx_0_z[i] = 4.0 * g_x_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_x_0_0_xy_0_x, g_x_0_x_0_0_xy_0_y, g_x_0_x_0_0_xy_0_z, g_x_xy_x_x, g_x_xy_x_y, g_x_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xy_0_x[i] = 4.0 * g_x_xy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_0_y[i] = 4.0 * g_x_xy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xy_0_z[i] = 4.0 * g_x_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_x_0_0_xz_0_x, g_x_0_x_0_0_xz_0_y, g_x_0_x_0_0_xz_0_z, g_x_xz_x_x, g_x_xz_x_y, g_x_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_xz_0_x[i] = 4.0 * g_x_xz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_0_y[i] = 4.0 * g_x_xz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_xz_0_z[i] = 4.0 * g_x_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_0_x_0_0_yy_0_x, g_x_0_x_0_0_yy_0_y, g_x_0_x_0_0_yy_0_z, g_x_yy_x_x, g_x_yy_x_y, g_x_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yy_0_x[i] = 4.0 * g_x_yy_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_0_y[i] = 4.0 * g_x_yy_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yy_0_z[i] = 4.0 * g_x_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_0_x_0_0_yz_0_x, g_x_0_x_0_0_yz_0_y, g_x_0_x_0_0_yz_0_z, g_x_yz_x_x, g_x_yz_x_y, g_x_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_yz_0_x[i] = 4.0 * g_x_yz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_0_y[i] = 4.0 * g_x_yz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_yz_0_z[i] = 4.0 * g_x_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_0_x_0_0_zz_0_x, g_x_0_x_0_0_zz_0_y, g_x_0_x_0_0_zz_0_z, g_x_zz_x_x, g_x_zz_x_y, g_x_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_zz_0_x[i] = 4.0 * g_x_zz_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_0_y[i] = 4.0 * g_x_zz_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_zz_0_z[i] = 4.0 * g_x_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_0_y_0_0_xx_0_x, g_x_0_y_0_0_xx_0_y, g_x_0_y_0_0_xx_0_z, g_x_xx_y_x, g_x_xx_y_y, g_x_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xx_0_x[i] = 4.0 * g_x_xx_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_0_y[i] = 4.0 * g_x_xx_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xx_0_z[i] = 4.0 * g_x_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_0_y_0_0_xy_0_x, g_x_0_y_0_0_xy_0_y, g_x_0_y_0_0_xy_0_z, g_x_xy_y_x, g_x_xy_y_y, g_x_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xy_0_x[i] = 4.0 * g_x_xy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_0_y[i] = 4.0 * g_x_xy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xy_0_z[i] = 4.0 * g_x_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_0_y_0_0_xz_0_x, g_x_0_y_0_0_xz_0_y, g_x_0_y_0_0_xz_0_z, g_x_xz_y_x, g_x_xz_y_y, g_x_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_xz_0_x[i] = 4.0 * g_x_xz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_0_y[i] = 4.0 * g_x_xz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_xz_0_z[i] = 4.0 * g_x_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_0_y_0_0_yy_0_x, g_x_0_y_0_0_yy_0_y, g_x_0_y_0_0_yy_0_z, g_x_yy_y_x, g_x_yy_y_y, g_x_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yy_0_x[i] = 4.0 * g_x_yy_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_0_y[i] = 4.0 * g_x_yy_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yy_0_z[i] = 4.0 * g_x_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_x_0_y_0_0_yz_0_x, g_x_0_y_0_0_yz_0_y, g_x_0_y_0_0_yz_0_z, g_x_yz_y_x, g_x_yz_y_y, g_x_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_yz_0_x[i] = 4.0 * g_x_yz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_0_y[i] = 4.0 * g_x_yz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_yz_0_z[i] = 4.0 * g_x_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_x_0_y_0_0_zz_0_x, g_x_0_y_0_0_zz_0_y, g_x_0_y_0_0_zz_0_z, g_x_zz_y_x, g_x_zz_y_y, g_x_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_zz_0_x[i] = 4.0 * g_x_zz_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_0_y[i] = 4.0 * g_x_zz_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_zz_0_z[i] = 4.0 * g_x_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_0_z_0_0_xx_0_x, g_x_0_z_0_0_xx_0_y, g_x_0_z_0_0_xx_0_z, g_x_xx_z_x, g_x_xx_z_y, g_x_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xx_0_x[i] = 4.0 * g_x_xx_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_0_y[i] = 4.0 * g_x_xx_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xx_0_z[i] = 4.0 * g_x_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_0_z_0_0_xy_0_x, g_x_0_z_0_0_xy_0_y, g_x_0_z_0_0_xy_0_z, g_x_xy_z_x, g_x_xy_z_y, g_x_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xy_0_x[i] = 4.0 * g_x_xy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_0_y[i] = 4.0 * g_x_xy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xy_0_z[i] = 4.0 * g_x_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_0_z_0_0_xz_0_x, g_x_0_z_0_0_xz_0_y, g_x_0_z_0_0_xz_0_z, g_x_xz_z_x, g_x_xz_z_y, g_x_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_xz_0_x[i] = 4.0 * g_x_xz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_0_y[i] = 4.0 * g_x_xz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_xz_0_z[i] = 4.0 * g_x_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_0_z_0_0_yy_0_x, g_x_0_z_0_0_yy_0_y, g_x_0_z_0_0_yy_0_z, g_x_yy_z_x, g_x_yy_z_y, g_x_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yy_0_x[i] = 4.0 * g_x_yy_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_0_y[i] = 4.0 * g_x_yy_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yy_0_z[i] = 4.0 * g_x_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_0_z_0_0_yz_0_x, g_x_0_z_0_0_yz_0_y, g_x_0_z_0_0_yz_0_z, g_x_yz_z_x, g_x_yz_z_y, g_x_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_yz_0_x[i] = 4.0 * g_x_yz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_0_y[i] = 4.0 * g_x_yz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_yz_0_z[i] = 4.0 * g_x_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_0_z_0_0_zz_0_x, g_x_0_z_0_0_zz_0_y, g_x_0_z_0_0_zz_0_z, g_x_zz_z_x, g_x_zz_z_y, g_x_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_zz_0_x[i] = 4.0 * g_x_zz_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_0_y[i] = 4.0 * g_x_zz_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_zz_0_z[i] = 4.0 * g_x_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_y_0_x_0_0_xx_0_x, g_y_0_x_0_0_xx_0_y, g_y_0_x_0_0_xx_0_z, g_y_xx_x_x, g_y_xx_x_y, g_y_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xx_0_x[i] = 4.0 * g_y_xx_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_0_y[i] = 4.0 * g_y_xx_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xx_0_z[i] = 4.0 * g_y_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_y_0_x_0_0_xy_0_x, g_y_0_x_0_0_xy_0_y, g_y_0_x_0_0_xy_0_z, g_y_xy_x_x, g_y_xy_x_y, g_y_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xy_0_x[i] = 4.0 * g_y_xy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_0_y[i] = 4.0 * g_y_xy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xy_0_z[i] = 4.0 * g_y_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_y_0_x_0_0_xz_0_x, g_y_0_x_0_0_xz_0_y, g_y_0_x_0_0_xz_0_z, g_y_xz_x_x, g_y_xz_x_y, g_y_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_xz_0_x[i] = 4.0 * g_y_xz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_0_y[i] = 4.0 * g_y_xz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_xz_0_z[i] = 4.0 * g_y_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_y_0_x_0_0_yy_0_x, g_y_0_x_0_0_yy_0_y, g_y_0_x_0_0_yy_0_z, g_y_yy_x_x, g_y_yy_x_y, g_y_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yy_0_x[i] = 4.0 * g_y_yy_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_0_y[i] = 4.0 * g_y_yy_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yy_0_z[i] = 4.0 * g_y_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_y_0_x_0_0_yz_0_x, g_y_0_x_0_0_yz_0_y, g_y_0_x_0_0_yz_0_z, g_y_yz_x_x, g_y_yz_x_y, g_y_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_yz_0_x[i] = 4.0 * g_y_yz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_0_y[i] = 4.0 * g_y_yz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_yz_0_z[i] = 4.0 * g_y_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_y_0_x_0_0_zz_0_x, g_y_0_x_0_0_zz_0_y, g_y_0_x_0_0_zz_0_z, g_y_zz_x_x, g_y_zz_x_y, g_y_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_zz_0_x[i] = 4.0 * g_y_zz_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_0_y[i] = 4.0 * g_y_zz_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_zz_0_z[i] = 4.0 * g_y_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_y_0_y_0_0_xx_0_x, g_y_0_y_0_0_xx_0_y, g_y_0_y_0_0_xx_0_z, g_y_xx_y_x, g_y_xx_y_y, g_y_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xx_0_x[i] = 4.0 * g_y_xx_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_0_y[i] = 4.0 * g_y_xx_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xx_0_z[i] = 4.0 * g_y_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_y_0_y_0_0_xy_0_x, g_y_0_y_0_0_xy_0_y, g_y_0_y_0_0_xy_0_z, g_y_xy_y_x, g_y_xy_y_y, g_y_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xy_0_x[i] = 4.0 * g_y_xy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_0_y[i] = 4.0 * g_y_xy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xy_0_z[i] = 4.0 * g_y_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_y_0_y_0_0_xz_0_x, g_y_0_y_0_0_xz_0_y, g_y_0_y_0_0_xz_0_z, g_y_xz_y_x, g_y_xz_y_y, g_y_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_xz_0_x[i] = 4.0 * g_y_xz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_0_y[i] = 4.0 * g_y_xz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_xz_0_z[i] = 4.0 * g_y_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_y_0_y_0_0_yy_0_x, g_y_0_y_0_0_yy_0_y, g_y_0_y_0_0_yy_0_z, g_y_yy_y_x, g_y_yy_y_y, g_y_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yy_0_x[i] = 4.0 * g_y_yy_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_0_y[i] = 4.0 * g_y_yy_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yy_0_z[i] = 4.0 * g_y_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_y_0_y_0_0_yz_0_x, g_y_0_y_0_0_yz_0_y, g_y_0_y_0_0_yz_0_z, g_y_yz_y_x, g_y_yz_y_y, g_y_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_yz_0_x[i] = 4.0 * g_y_yz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_0_y[i] = 4.0 * g_y_yz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_yz_0_z[i] = 4.0 * g_y_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_y_0_y_0_0_zz_0_x, g_y_0_y_0_0_zz_0_y, g_y_0_y_0_0_zz_0_z, g_y_zz_y_x, g_y_zz_y_y, g_y_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_zz_0_x[i] = 4.0 * g_y_zz_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_0_y[i] = 4.0 * g_y_zz_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_zz_0_z[i] = 4.0 * g_y_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_y_0_z_0_0_xx_0_x, g_y_0_z_0_0_xx_0_y, g_y_0_z_0_0_xx_0_z, g_y_xx_z_x, g_y_xx_z_y, g_y_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xx_0_x[i] = 4.0 * g_y_xx_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_0_y[i] = 4.0 * g_y_xx_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xx_0_z[i] = 4.0 * g_y_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_y_0_z_0_0_xy_0_x, g_y_0_z_0_0_xy_0_y, g_y_0_z_0_0_xy_0_z, g_y_xy_z_x, g_y_xy_z_y, g_y_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xy_0_x[i] = 4.0 * g_y_xy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_0_y[i] = 4.0 * g_y_xy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xy_0_z[i] = 4.0 * g_y_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_y_0_z_0_0_xz_0_x, g_y_0_z_0_0_xz_0_y, g_y_0_z_0_0_xz_0_z, g_y_xz_z_x, g_y_xz_z_y, g_y_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_xz_0_x[i] = 4.0 * g_y_xz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_0_y[i] = 4.0 * g_y_xz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_xz_0_z[i] = 4.0 * g_y_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_y_0_z_0_0_yy_0_x, g_y_0_z_0_0_yy_0_y, g_y_0_z_0_0_yy_0_z, g_y_yy_z_x, g_y_yy_z_y, g_y_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yy_0_x[i] = 4.0 * g_y_yy_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_0_y[i] = 4.0 * g_y_yy_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yy_0_z[i] = 4.0 * g_y_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_y_0_z_0_0_yz_0_x, g_y_0_z_0_0_yz_0_y, g_y_0_z_0_0_yz_0_z, g_y_yz_z_x, g_y_yz_z_y, g_y_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_yz_0_x[i] = 4.0 * g_y_yz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_0_y[i] = 4.0 * g_y_yz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_yz_0_z[i] = 4.0 * g_y_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_y_0_z_0_0_zz_0_x, g_y_0_z_0_0_zz_0_y, g_y_0_z_0_0_zz_0_z, g_y_zz_z_x, g_y_zz_z_y, g_y_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_zz_0_x[i] = 4.0 * g_y_zz_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_0_y[i] = 4.0 * g_y_zz_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_zz_0_z[i] = 4.0 * g_y_zz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_z_0_x_0_0_xx_0_x, g_z_0_x_0_0_xx_0_y, g_z_0_x_0_0_xx_0_z, g_z_xx_x_x, g_z_xx_x_y, g_z_xx_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xx_0_x[i] = 4.0 * g_z_xx_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_0_y[i] = 4.0 * g_z_xx_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xx_0_z[i] = 4.0 * g_z_xx_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_z_0_x_0_0_xy_0_x, g_z_0_x_0_0_xy_0_y, g_z_0_x_0_0_xy_0_z, g_z_xy_x_x, g_z_xy_x_y, g_z_xy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xy_0_x[i] = 4.0 * g_z_xy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_0_y[i] = 4.0 * g_z_xy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xy_0_z[i] = 4.0 * g_z_xy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_z_0_x_0_0_xz_0_x, g_z_0_x_0_0_xz_0_y, g_z_0_x_0_0_xz_0_z, g_z_xz_x_x, g_z_xz_x_y, g_z_xz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_xz_0_x[i] = 4.0 * g_z_xz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_0_y[i] = 4.0 * g_z_xz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_xz_0_z[i] = 4.0 * g_z_xz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_z_0_x_0_0_yy_0_x, g_z_0_x_0_0_yy_0_y, g_z_0_x_0_0_yy_0_z, g_z_yy_x_x, g_z_yy_x_y, g_z_yy_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yy_0_x[i] = 4.0 * g_z_yy_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_0_y[i] = 4.0 * g_z_yy_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yy_0_z[i] = 4.0 * g_z_yy_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_z_0_x_0_0_yz_0_x, g_z_0_x_0_0_yz_0_y, g_z_0_x_0_0_yz_0_z, g_z_yz_x_x, g_z_yz_x_y, g_z_yz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_yz_0_x[i] = 4.0 * g_z_yz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_0_y[i] = 4.0 * g_z_yz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_yz_0_z[i] = 4.0 * g_z_yz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_z_0_x_0_0_zz_0_x, g_z_0_x_0_0_zz_0_y, g_z_0_x_0_0_zz_0_z, g_z_zz_x_x, g_z_zz_x_y, g_z_zz_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_zz_0_x[i] = 4.0 * g_z_zz_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_0_y[i] = 4.0 * g_z_zz_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_zz_0_z[i] = 4.0 * g_z_zz_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_z_0_y_0_0_xx_0_x, g_z_0_y_0_0_xx_0_y, g_z_0_y_0_0_xx_0_z, g_z_xx_y_x, g_z_xx_y_y, g_z_xx_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xx_0_x[i] = 4.0 * g_z_xx_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_0_y[i] = 4.0 * g_z_xx_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xx_0_z[i] = 4.0 * g_z_xx_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_z_0_y_0_0_xy_0_x, g_z_0_y_0_0_xy_0_y, g_z_0_y_0_0_xy_0_z, g_z_xy_y_x, g_z_xy_y_y, g_z_xy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xy_0_x[i] = 4.0 * g_z_xy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_0_y[i] = 4.0 * g_z_xy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xy_0_z[i] = 4.0 * g_z_xy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_z_0_y_0_0_xz_0_x, g_z_0_y_0_0_xz_0_y, g_z_0_y_0_0_xz_0_z, g_z_xz_y_x, g_z_xz_y_y, g_z_xz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_xz_0_x[i] = 4.0 * g_z_xz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_0_y[i] = 4.0 * g_z_xz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_xz_0_z[i] = 4.0 * g_z_xz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_z_0_y_0_0_yy_0_x, g_z_0_y_0_0_yy_0_y, g_z_0_y_0_0_yy_0_z, g_z_yy_y_x, g_z_yy_y_y, g_z_yy_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yy_0_x[i] = 4.0 * g_z_yy_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_0_y[i] = 4.0 * g_z_yy_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yy_0_z[i] = 4.0 * g_z_yy_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_z_0_y_0_0_yz_0_x, g_z_0_y_0_0_yz_0_y, g_z_0_y_0_0_yz_0_z, g_z_yz_y_x, g_z_yz_y_y, g_z_yz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_yz_0_x[i] = 4.0 * g_z_yz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_0_y[i] = 4.0 * g_z_yz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_yz_0_z[i] = 4.0 * g_z_yz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_z_0_y_0_0_zz_0_x, g_z_0_y_0_0_zz_0_y, g_z_0_y_0_0_zz_0_z, g_z_zz_y_x, g_z_zz_y_y, g_z_zz_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_zz_0_x[i] = 4.0 * g_z_zz_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_0_y[i] = 4.0 * g_z_zz_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_zz_0_z[i] = 4.0 * g_z_zz_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_z_0_z_0_0_xx_0_x, g_z_0_z_0_0_xx_0_y, g_z_0_z_0_0_xx_0_z, g_z_xx_z_x, g_z_xx_z_y, g_z_xx_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xx_0_x[i] = 4.0 * g_z_xx_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_0_y[i] = 4.0 * g_z_xx_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xx_0_z[i] = 4.0 * g_z_xx_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_z_0_z_0_0_xy_0_x, g_z_0_z_0_0_xy_0_y, g_z_0_z_0_0_xy_0_z, g_z_xy_z_x, g_z_xy_z_y, g_z_xy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xy_0_x[i] = 4.0 * g_z_xy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_0_y[i] = 4.0 * g_z_xy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xy_0_z[i] = 4.0 * g_z_xy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_z_0_z_0_0_xz_0_x, g_z_0_z_0_0_xz_0_y, g_z_0_z_0_0_xz_0_z, g_z_xz_z_x, g_z_xz_z_y, g_z_xz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_xz_0_x[i] = 4.0 * g_z_xz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_0_y[i] = 4.0 * g_z_xz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_xz_0_z[i] = 4.0 * g_z_xz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_z_0_z_0_0_yy_0_x, g_z_0_z_0_0_yy_0_y, g_z_0_z_0_0_yy_0_z, g_z_yy_z_x, g_z_yy_z_y, g_z_yy_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yy_0_x[i] = 4.0 * g_z_yy_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_0_y[i] = 4.0 * g_z_yy_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yy_0_z[i] = 4.0 * g_z_yy_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_z_0_z_0_0_yz_0_x, g_z_0_z_0_0_yz_0_y, g_z_0_z_0_0_yz_0_z, g_z_yz_z_x, g_z_yz_z_y, g_z_yz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_yz_0_x[i] = 4.0 * g_z_yz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_0_y[i] = 4.0 * g_z_yz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_yz_0_z[i] = 4.0 * g_z_yz_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_z_0_z_0_0_zz_0_x, g_z_0_z_0_0_zz_0_y, g_z_0_z_0_0_zz_0_z, g_z_zz_z_x, g_z_zz_z_y, g_z_zz_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_zz_0_x[i] = 4.0 * g_z_zz_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_0_y[i] = 4.0 * g_z_zz_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_zz_0_z[i] = 4.0 * g_z_zz_z_z[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace

