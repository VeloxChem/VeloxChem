#include "GeomDeriv1000OfScalarForSDPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_sdpp_0(CSimdArray<double>& buffer_1000_sdpp,
                     const CSimdArray<double>& buffer_pdpp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_sdpp.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_sdpp

    auto g_x_0_0_0_0_xx_x_x = buffer_1000_sdpp[0];

    auto g_x_0_0_0_0_xx_x_y = buffer_1000_sdpp[1];

    auto g_x_0_0_0_0_xx_x_z = buffer_1000_sdpp[2];

    auto g_x_0_0_0_0_xx_y_x = buffer_1000_sdpp[3];

    auto g_x_0_0_0_0_xx_y_y = buffer_1000_sdpp[4];

    auto g_x_0_0_0_0_xx_y_z = buffer_1000_sdpp[5];

    auto g_x_0_0_0_0_xx_z_x = buffer_1000_sdpp[6];

    auto g_x_0_0_0_0_xx_z_y = buffer_1000_sdpp[7];

    auto g_x_0_0_0_0_xx_z_z = buffer_1000_sdpp[8];

    auto g_x_0_0_0_0_xy_x_x = buffer_1000_sdpp[9];

    auto g_x_0_0_0_0_xy_x_y = buffer_1000_sdpp[10];

    auto g_x_0_0_0_0_xy_x_z = buffer_1000_sdpp[11];

    auto g_x_0_0_0_0_xy_y_x = buffer_1000_sdpp[12];

    auto g_x_0_0_0_0_xy_y_y = buffer_1000_sdpp[13];

    auto g_x_0_0_0_0_xy_y_z = buffer_1000_sdpp[14];

    auto g_x_0_0_0_0_xy_z_x = buffer_1000_sdpp[15];

    auto g_x_0_0_0_0_xy_z_y = buffer_1000_sdpp[16];

    auto g_x_0_0_0_0_xy_z_z = buffer_1000_sdpp[17];

    auto g_x_0_0_0_0_xz_x_x = buffer_1000_sdpp[18];

    auto g_x_0_0_0_0_xz_x_y = buffer_1000_sdpp[19];

    auto g_x_0_0_0_0_xz_x_z = buffer_1000_sdpp[20];

    auto g_x_0_0_0_0_xz_y_x = buffer_1000_sdpp[21];

    auto g_x_0_0_0_0_xz_y_y = buffer_1000_sdpp[22];

    auto g_x_0_0_0_0_xz_y_z = buffer_1000_sdpp[23];

    auto g_x_0_0_0_0_xz_z_x = buffer_1000_sdpp[24];

    auto g_x_0_0_0_0_xz_z_y = buffer_1000_sdpp[25];

    auto g_x_0_0_0_0_xz_z_z = buffer_1000_sdpp[26];

    auto g_x_0_0_0_0_yy_x_x = buffer_1000_sdpp[27];

    auto g_x_0_0_0_0_yy_x_y = buffer_1000_sdpp[28];

    auto g_x_0_0_0_0_yy_x_z = buffer_1000_sdpp[29];

    auto g_x_0_0_0_0_yy_y_x = buffer_1000_sdpp[30];

    auto g_x_0_0_0_0_yy_y_y = buffer_1000_sdpp[31];

    auto g_x_0_0_0_0_yy_y_z = buffer_1000_sdpp[32];

    auto g_x_0_0_0_0_yy_z_x = buffer_1000_sdpp[33];

    auto g_x_0_0_0_0_yy_z_y = buffer_1000_sdpp[34];

    auto g_x_0_0_0_0_yy_z_z = buffer_1000_sdpp[35];

    auto g_x_0_0_0_0_yz_x_x = buffer_1000_sdpp[36];

    auto g_x_0_0_0_0_yz_x_y = buffer_1000_sdpp[37];

    auto g_x_0_0_0_0_yz_x_z = buffer_1000_sdpp[38];

    auto g_x_0_0_0_0_yz_y_x = buffer_1000_sdpp[39];

    auto g_x_0_0_0_0_yz_y_y = buffer_1000_sdpp[40];

    auto g_x_0_0_0_0_yz_y_z = buffer_1000_sdpp[41];

    auto g_x_0_0_0_0_yz_z_x = buffer_1000_sdpp[42];

    auto g_x_0_0_0_0_yz_z_y = buffer_1000_sdpp[43];

    auto g_x_0_0_0_0_yz_z_z = buffer_1000_sdpp[44];

    auto g_x_0_0_0_0_zz_x_x = buffer_1000_sdpp[45];

    auto g_x_0_0_0_0_zz_x_y = buffer_1000_sdpp[46];

    auto g_x_0_0_0_0_zz_x_z = buffer_1000_sdpp[47];

    auto g_x_0_0_0_0_zz_y_x = buffer_1000_sdpp[48];

    auto g_x_0_0_0_0_zz_y_y = buffer_1000_sdpp[49];

    auto g_x_0_0_0_0_zz_y_z = buffer_1000_sdpp[50];

    auto g_x_0_0_0_0_zz_z_x = buffer_1000_sdpp[51];

    auto g_x_0_0_0_0_zz_z_y = buffer_1000_sdpp[52];

    auto g_x_0_0_0_0_zz_z_z = buffer_1000_sdpp[53];

    auto g_y_0_0_0_0_xx_x_x = buffer_1000_sdpp[54];

    auto g_y_0_0_0_0_xx_x_y = buffer_1000_sdpp[55];

    auto g_y_0_0_0_0_xx_x_z = buffer_1000_sdpp[56];

    auto g_y_0_0_0_0_xx_y_x = buffer_1000_sdpp[57];

    auto g_y_0_0_0_0_xx_y_y = buffer_1000_sdpp[58];

    auto g_y_0_0_0_0_xx_y_z = buffer_1000_sdpp[59];

    auto g_y_0_0_0_0_xx_z_x = buffer_1000_sdpp[60];

    auto g_y_0_0_0_0_xx_z_y = buffer_1000_sdpp[61];

    auto g_y_0_0_0_0_xx_z_z = buffer_1000_sdpp[62];

    auto g_y_0_0_0_0_xy_x_x = buffer_1000_sdpp[63];

    auto g_y_0_0_0_0_xy_x_y = buffer_1000_sdpp[64];

    auto g_y_0_0_0_0_xy_x_z = buffer_1000_sdpp[65];

    auto g_y_0_0_0_0_xy_y_x = buffer_1000_sdpp[66];

    auto g_y_0_0_0_0_xy_y_y = buffer_1000_sdpp[67];

    auto g_y_0_0_0_0_xy_y_z = buffer_1000_sdpp[68];

    auto g_y_0_0_0_0_xy_z_x = buffer_1000_sdpp[69];

    auto g_y_0_0_0_0_xy_z_y = buffer_1000_sdpp[70];

    auto g_y_0_0_0_0_xy_z_z = buffer_1000_sdpp[71];

    auto g_y_0_0_0_0_xz_x_x = buffer_1000_sdpp[72];

    auto g_y_0_0_0_0_xz_x_y = buffer_1000_sdpp[73];

    auto g_y_0_0_0_0_xz_x_z = buffer_1000_sdpp[74];

    auto g_y_0_0_0_0_xz_y_x = buffer_1000_sdpp[75];

    auto g_y_0_0_0_0_xz_y_y = buffer_1000_sdpp[76];

    auto g_y_0_0_0_0_xz_y_z = buffer_1000_sdpp[77];

    auto g_y_0_0_0_0_xz_z_x = buffer_1000_sdpp[78];

    auto g_y_0_0_0_0_xz_z_y = buffer_1000_sdpp[79];

    auto g_y_0_0_0_0_xz_z_z = buffer_1000_sdpp[80];

    auto g_y_0_0_0_0_yy_x_x = buffer_1000_sdpp[81];

    auto g_y_0_0_0_0_yy_x_y = buffer_1000_sdpp[82];

    auto g_y_0_0_0_0_yy_x_z = buffer_1000_sdpp[83];

    auto g_y_0_0_0_0_yy_y_x = buffer_1000_sdpp[84];

    auto g_y_0_0_0_0_yy_y_y = buffer_1000_sdpp[85];

    auto g_y_0_0_0_0_yy_y_z = buffer_1000_sdpp[86];

    auto g_y_0_0_0_0_yy_z_x = buffer_1000_sdpp[87];

    auto g_y_0_0_0_0_yy_z_y = buffer_1000_sdpp[88];

    auto g_y_0_0_0_0_yy_z_z = buffer_1000_sdpp[89];

    auto g_y_0_0_0_0_yz_x_x = buffer_1000_sdpp[90];

    auto g_y_0_0_0_0_yz_x_y = buffer_1000_sdpp[91];

    auto g_y_0_0_0_0_yz_x_z = buffer_1000_sdpp[92];

    auto g_y_0_0_0_0_yz_y_x = buffer_1000_sdpp[93];

    auto g_y_0_0_0_0_yz_y_y = buffer_1000_sdpp[94];

    auto g_y_0_0_0_0_yz_y_z = buffer_1000_sdpp[95];

    auto g_y_0_0_0_0_yz_z_x = buffer_1000_sdpp[96];

    auto g_y_0_0_0_0_yz_z_y = buffer_1000_sdpp[97];

    auto g_y_0_0_0_0_yz_z_z = buffer_1000_sdpp[98];

    auto g_y_0_0_0_0_zz_x_x = buffer_1000_sdpp[99];

    auto g_y_0_0_0_0_zz_x_y = buffer_1000_sdpp[100];

    auto g_y_0_0_0_0_zz_x_z = buffer_1000_sdpp[101];

    auto g_y_0_0_0_0_zz_y_x = buffer_1000_sdpp[102];

    auto g_y_0_0_0_0_zz_y_y = buffer_1000_sdpp[103];

    auto g_y_0_0_0_0_zz_y_z = buffer_1000_sdpp[104];

    auto g_y_0_0_0_0_zz_z_x = buffer_1000_sdpp[105];

    auto g_y_0_0_0_0_zz_z_y = buffer_1000_sdpp[106];

    auto g_y_0_0_0_0_zz_z_z = buffer_1000_sdpp[107];

    auto g_z_0_0_0_0_xx_x_x = buffer_1000_sdpp[108];

    auto g_z_0_0_0_0_xx_x_y = buffer_1000_sdpp[109];

    auto g_z_0_0_0_0_xx_x_z = buffer_1000_sdpp[110];

    auto g_z_0_0_0_0_xx_y_x = buffer_1000_sdpp[111];

    auto g_z_0_0_0_0_xx_y_y = buffer_1000_sdpp[112];

    auto g_z_0_0_0_0_xx_y_z = buffer_1000_sdpp[113];

    auto g_z_0_0_0_0_xx_z_x = buffer_1000_sdpp[114];

    auto g_z_0_0_0_0_xx_z_y = buffer_1000_sdpp[115];

    auto g_z_0_0_0_0_xx_z_z = buffer_1000_sdpp[116];

    auto g_z_0_0_0_0_xy_x_x = buffer_1000_sdpp[117];

    auto g_z_0_0_0_0_xy_x_y = buffer_1000_sdpp[118];

    auto g_z_0_0_0_0_xy_x_z = buffer_1000_sdpp[119];

    auto g_z_0_0_0_0_xy_y_x = buffer_1000_sdpp[120];

    auto g_z_0_0_0_0_xy_y_y = buffer_1000_sdpp[121];

    auto g_z_0_0_0_0_xy_y_z = buffer_1000_sdpp[122];

    auto g_z_0_0_0_0_xy_z_x = buffer_1000_sdpp[123];

    auto g_z_0_0_0_0_xy_z_y = buffer_1000_sdpp[124];

    auto g_z_0_0_0_0_xy_z_z = buffer_1000_sdpp[125];

    auto g_z_0_0_0_0_xz_x_x = buffer_1000_sdpp[126];

    auto g_z_0_0_0_0_xz_x_y = buffer_1000_sdpp[127];

    auto g_z_0_0_0_0_xz_x_z = buffer_1000_sdpp[128];

    auto g_z_0_0_0_0_xz_y_x = buffer_1000_sdpp[129];

    auto g_z_0_0_0_0_xz_y_y = buffer_1000_sdpp[130];

    auto g_z_0_0_0_0_xz_y_z = buffer_1000_sdpp[131];

    auto g_z_0_0_0_0_xz_z_x = buffer_1000_sdpp[132];

    auto g_z_0_0_0_0_xz_z_y = buffer_1000_sdpp[133];

    auto g_z_0_0_0_0_xz_z_z = buffer_1000_sdpp[134];

    auto g_z_0_0_0_0_yy_x_x = buffer_1000_sdpp[135];

    auto g_z_0_0_0_0_yy_x_y = buffer_1000_sdpp[136];

    auto g_z_0_0_0_0_yy_x_z = buffer_1000_sdpp[137];

    auto g_z_0_0_0_0_yy_y_x = buffer_1000_sdpp[138];

    auto g_z_0_0_0_0_yy_y_y = buffer_1000_sdpp[139];

    auto g_z_0_0_0_0_yy_y_z = buffer_1000_sdpp[140];

    auto g_z_0_0_0_0_yy_z_x = buffer_1000_sdpp[141];

    auto g_z_0_0_0_0_yy_z_y = buffer_1000_sdpp[142];

    auto g_z_0_0_0_0_yy_z_z = buffer_1000_sdpp[143];

    auto g_z_0_0_0_0_yz_x_x = buffer_1000_sdpp[144];

    auto g_z_0_0_0_0_yz_x_y = buffer_1000_sdpp[145];

    auto g_z_0_0_0_0_yz_x_z = buffer_1000_sdpp[146];

    auto g_z_0_0_0_0_yz_y_x = buffer_1000_sdpp[147];

    auto g_z_0_0_0_0_yz_y_y = buffer_1000_sdpp[148];

    auto g_z_0_0_0_0_yz_y_z = buffer_1000_sdpp[149];

    auto g_z_0_0_0_0_yz_z_x = buffer_1000_sdpp[150];

    auto g_z_0_0_0_0_yz_z_y = buffer_1000_sdpp[151];

    auto g_z_0_0_0_0_yz_z_z = buffer_1000_sdpp[152];

    auto g_z_0_0_0_0_zz_x_x = buffer_1000_sdpp[153];

    auto g_z_0_0_0_0_zz_x_y = buffer_1000_sdpp[154];

    auto g_z_0_0_0_0_zz_x_z = buffer_1000_sdpp[155];

    auto g_z_0_0_0_0_zz_y_x = buffer_1000_sdpp[156];

    auto g_z_0_0_0_0_zz_y_y = buffer_1000_sdpp[157];

    auto g_z_0_0_0_0_zz_y_z = buffer_1000_sdpp[158];

    auto g_z_0_0_0_0_zz_z_x = buffer_1000_sdpp[159];

    auto g_z_0_0_0_0_zz_z_y = buffer_1000_sdpp[160];

    auto g_z_0_0_0_0_zz_z_z = buffer_1000_sdpp[161];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_0_0_0_xx_x_x, g_x_0_0_0_0_xx_x_y, g_x_0_0_0_0_xx_x_z, g_x_xx_x_x, g_x_xx_x_y, g_x_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xx_x_x[i] = 2.0 * g_x_xx_x_x[i] * a_exp;

        g_x_0_0_0_0_xx_x_y[i] = 2.0 * g_x_xx_x_y[i] * a_exp;

        g_x_0_0_0_0_xx_x_z[i] = 2.0 * g_x_xx_x_z[i] * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_0_0_0_xx_y_x, g_x_0_0_0_0_xx_y_y, g_x_0_0_0_0_xx_y_z, g_x_xx_y_x, g_x_xx_y_y, g_x_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xx_y_x[i] = 2.0 * g_x_xx_y_x[i] * a_exp;

        g_x_0_0_0_0_xx_y_y[i] = 2.0 * g_x_xx_y_y[i] * a_exp;

        g_x_0_0_0_0_xx_y_z[i] = 2.0 * g_x_xx_y_z[i] * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_0_0_0_xx_z_x, g_x_0_0_0_0_xx_z_y, g_x_0_0_0_0_xx_z_z, g_x_xx_z_x, g_x_xx_z_y, g_x_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xx_z_x[i] = 2.0 * g_x_xx_z_x[i] * a_exp;

        g_x_0_0_0_0_xx_z_y[i] = 2.0 * g_x_xx_z_y[i] * a_exp;

        g_x_0_0_0_0_xx_z_z[i] = 2.0 * g_x_xx_z_z[i] * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_0_0_0_0_xy_x_x, g_x_0_0_0_0_xy_x_y, g_x_0_0_0_0_xy_x_z, g_x_xy_x_x, g_x_xy_x_y, g_x_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xy_x_x[i] = 2.0 * g_x_xy_x_x[i] * a_exp;

        g_x_0_0_0_0_xy_x_y[i] = 2.0 * g_x_xy_x_y[i] * a_exp;

        g_x_0_0_0_0_xy_x_z[i] = 2.0 * g_x_xy_x_z[i] * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_0_0_0_0_xy_y_x, g_x_0_0_0_0_xy_y_y, g_x_0_0_0_0_xy_y_z, g_x_xy_y_x, g_x_xy_y_y, g_x_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xy_y_x[i] = 2.0 * g_x_xy_y_x[i] * a_exp;

        g_x_0_0_0_0_xy_y_y[i] = 2.0 * g_x_xy_y_y[i] * a_exp;

        g_x_0_0_0_0_xy_y_z[i] = 2.0 * g_x_xy_y_z[i] * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_0_0_0_0_xy_z_x, g_x_0_0_0_0_xy_z_y, g_x_0_0_0_0_xy_z_z, g_x_xy_z_x, g_x_xy_z_y, g_x_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xy_z_x[i] = 2.0 * g_x_xy_z_x[i] * a_exp;

        g_x_0_0_0_0_xy_z_y[i] = 2.0 * g_x_xy_z_y[i] * a_exp;

        g_x_0_0_0_0_xy_z_z[i] = 2.0 * g_x_xy_z_z[i] * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_0_0_0_0_xz_x_x, g_x_0_0_0_0_xz_x_y, g_x_0_0_0_0_xz_x_z, g_x_xz_x_x, g_x_xz_x_y, g_x_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xz_x_x[i] = 2.0 * g_x_xz_x_x[i] * a_exp;

        g_x_0_0_0_0_xz_x_y[i] = 2.0 * g_x_xz_x_y[i] * a_exp;

        g_x_0_0_0_0_xz_x_z[i] = 2.0 * g_x_xz_x_z[i] * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_0_0_0_0_xz_y_x, g_x_0_0_0_0_xz_y_y, g_x_0_0_0_0_xz_y_z, g_x_xz_y_x, g_x_xz_y_y, g_x_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xz_y_x[i] = 2.0 * g_x_xz_y_x[i] * a_exp;

        g_x_0_0_0_0_xz_y_y[i] = 2.0 * g_x_xz_y_y[i] * a_exp;

        g_x_0_0_0_0_xz_y_z[i] = 2.0 * g_x_xz_y_z[i] * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_0_0_0_0_xz_z_x, g_x_0_0_0_0_xz_z_y, g_x_0_0_0_0_xz_z_z, g_x_xz_z_x, g_x_xz_z_y, g_x_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xz_z_x[i] = 2.0 * g_x_xz_z_x[i] * a_exp;

        g_x_0_0_0_0_xz_z_y[i] = 2.0 * g_x_xz_z_y[i] * a_exp;

        g_x_0_0_0_0_xz_z_z[i] = 2.0 * g_x_xz_z_z[i] * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_0_0_0_0_yy_x_x, g_x_0_0_0_0_yy_x_y, g_x_0_0_0_0_yy_x_z, g_x_yy_x_x, g_x_yy_x_y, g_x_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yy_x_x[i] = 2.0 * g_x_yy_x_x[i] * a_exp;

        g_x_0_0_0_0_yy_x_y[i] = 2.0 * g_x_yy_x_y[i] * a_exp;

        g_x_0_0_0_0_yy_x_z[i] = 2.0 * g_x_yy_x_z[i] * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_x_0_0_0_0_yy_y_x, g_x_0_0_0_0_yy_y_y, g_x_0_0_0_0_yy_y_z, g_x_yy_y_x, g_x_yy_y_y, g_x_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yy_y_x[i] = 2.0 * g_x_yy_y_x[i] * a_exp;

        g_x_0_0_0_0_yy_y_y[i] = 2.0 * g_x_yy_y_y[i] * a_exp;

        g_x_0_0_0_0_yy_y_z[i] = 2.0 * g_x_yy_y_z[i] * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_x_0_0_0_0_yy_z_x, g_x_0_0_0_0_yy_z_y, g_x_0_0_0_0_yy_z_z, g_x_yy_z_x, g_x_yy_z_y, g_x_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yy_z_x[i] = 2.0 * g_x_yy_z_x[i] * a_exp;

        g_x_0_0_0_0_yy_z_y[i] = 2.0 * g_x_yy_z_y[i] * a_exp;

        g_x_0_0_0_0_yy_z_z[i] = 2.0 * g_x_yy_z_z[i] * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_0_0_0_0_yz_x_x, g_x_0_0_0_0_yz_x_y, g_x_0_0_0_0_yz_x_z, g_x_yz_x_x, g_x_yz_x_y, g_x_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yz_x_x[i] = 2.0 * g_x_yz_x_x[i] * a_exp;

        g_x_0_0_0_0_yz_x_y[i] = 2.0 * g_x_yz_x_y[i] * a_exp;

        g_x_0_0_0_0_yz_x_z[i] = 2.0 * g_x_yz_x_z[i] * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_0_0_0_0_yz_y_x, g_x_0_0_0_0_yz_y_y, g_x_0_0_0_0_yz_y_z, g_x_yz_y_x, g_x_yz_y_y, g_x_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yz_y_x[i] = 2.0 * g_x_yz_y_x[i] * a_exp;

        g_x_0_0_0_0_yz_y_y[i] = 2.0 * g_x_yz_y_y[i] * a_exp;

        g_x_0_0_0_0_yz_y_z[i] = 2.0 * g_x_yz_y_z[i] * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_0_0_0_0_yz_z_x, g_x_0_0_0_0_yz_z_y, g_x_0_0_0_0_yz_z_z, g_x_yz_z_x, g_x_yz_z_y, g_x_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yz_z_x[i] = 2.0 * g_x_yz_z_x[i] * a_exp;

        g_x_0_0_0_0_yz_z_y[i] = 2.0 * g_x_yz_z_y[i] * a_exp;

        g_x_0_0_0_0_yz_z_z[i] = 2.0 * g_x_yz_z_z[i] * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_0_0_0_0_zz_x_x, g_x_0_0_0_0_zz_x_y, g_x_0_0_0_0_zz_x_z, g_x_zz_x_x, g_x_zz_x_y, g_x_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_zz_x_x[i] = 2.0 * g_x_zz_x_x[i] * a_exp;

        g_x_0_0_0_0_zz_x_y[i] = 2.0 * g_x_zz_x_y[i] * a_exp;

        g_x_0_0_0_0_zz_x_z[i] = 2.0 * g_x_zz_x_z[i] * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_0_0_0_0_zz_y_x, g_x_0_0_0_0_zz_y_y, g_x_0_0_0_0_zz_y_z, g_x_zz_y_x, g_x_zz_y_y, g_x_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_zz_y_x[i] = 2.0 * g_x_zz_y_x[i] * a_exp;

        g_x_0_0_0_0_zz_y_y[i] = 2.0 * g_x_zz_y_y[i] * a_exp;

        g_x_0_0_0_0_zz_y_z[i] = 2.0 * g_x_zz_y_z[i] * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_0_0_0_0_zz_z_x, g_x_0_0_0_0_zz_z_y, g_x_0_0_0_0_zz_z_z, g_x_zz_z_x, g_x_zz_z_y, g_x_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_zz_z_x[i] = 2.0 * g_x_zz_z_x[i] * a_exp;

        g_x_0_0_0_0_zz_z_y[i] = 2.0 * g_x_zz_z_y[i] * a_exp;

        g_x_0_0_0_0_zz_z_z[i] = 2.0 * g_x_zz_z_z[i] * a_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_y_0_0_0_0_xx_x_x, g_y_0_0_0_0_xx_x_y, g_y_0_0_0_0_xx_x_z, g_y_xx_x_x, g_y_xx_x_y, g_y_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xx_x_x[i] = 2.0 * g_y_xx_x_x[i] * a_exp;

        g_y_0_0_0_0_xx_x_y[i] = 2.0 * g_y_xx_x_y[i] * a_exp;

        g_y_0_0_0_0_xx_x_z[i] = 2.0 * g_y_xx_x_z[i] * a_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_y_0_0_0_0_xx_y_x, g_y_0_0_0_0_xx_y_y, g_y_0_0_0_0_xx_y_z, g_y_xx_y_x, g_y_xx_y_y, g_y_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xx_y_x[i] = 2.0 * g_y_xx_y_x[i] * a_exp;

        g_y_0_0_0_0_xx_y_y[i] = 2.0 * g_y_xx_y_y[i] * a_exp;

        g_y_0_0_0_0_xx_y_z[i] = 2.0 * g_y_xx_y_z[i] * a_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_y_0_0_0_0_xx_z_x, g_y_0_0_0_0_xx_z_y, g_y_0_0_0_0_xx_z_z, g_y_xx_z_x, g_y_xx_z_y, g_y_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xx_z_x[i] = 2.0 * g_y_xx_z_x[i] * a_exp;

        g_y_0_0_0_0_xx_z_y[i] = 2.0 * g_y_xx_z_y[i] * a_exp;

        g_y_0_0_0_0_xx_z_z[i] = 2.0 * g_y_xx_z_z[i] * a_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_y_0_0_0_0_xy_x_x, g_y_0_0_0_0_xy_x_y, g_y_0_0_0_0_xy_x_z, g_y_xy_x_x, g_y_xy_x_y, g_y_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xy_x_x[i] = 2.0 * g_y_xy_x_x[i] * a_exp;

        g_y_0_0_0_0_xy_x_y[i] = 2.0 * g_y_xy_x_y[i] * a_exp;

        g_y_0_0_0_0_xy_x_z[i] = 2.0 * g_y_xy_x_z[i] * a_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_y_0_0_0_0_xy_y_x, g_y_0_0_0_0_xy_y_y, g_y_0_0_0_0_xy_y_z, g_y_xy_y_x, g_y_xy_y_y, g_y_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xy_y_x[i] = 2.0 * g_y_xy_y_x[i] * a_exp;

        g_y_0_0_0_0_xy_y_y[i] = 2.0 * g_y_xy_y_y[i] * a_exp;

        g_y_0_0_0_0_xy_y_z[i] = 2.0 * g_y_xy_y_z[i] * a_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_y_0_0_0_0_xy_z_x, g_y_0_0_0_0_xy_z_y, g_y_0_0_0_0_xy_z_z, g_y_xy_z_x, g_y_xy_z_y, g_y_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xy_z_x[i] = 2.0 * g_y_xy_z_x[i] * a_exp;

        g_y_0_0_0_0_xy_z_y[i] = 2.0 * g_y_xy_z_y[i] * a_exp;

        g_y_0_0_0_0_xy_z_z[i] = 2.0 * g_y_xy_z_z[i] * a_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_y_0_0_0_0_xz_x_x, g_y_0_0_0_0_xz_x_y, g_y_0_0_0_0_xz_x_z, g_y_xz_x_x, g_y_xz_x_y, g_y_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xz_x_x[i] = 2.0 * g_y_xz_x_x[i] * a_exp;

        g_y_0_0_0_0_xz_x_y[i] = 2.0 * g_y_xz_x_y[i] * a_exp;

        g_y_0_0_0_0_xz_x_z[i] = 2.0 * g_y_xz_x_z[i] * a_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_y_0_0_0_0_xz_y_x, g_y_0_0_0_0_xz_y_y, g_y_0_0_0_0_xz_y_z, g_y_xz_y_x, g_y_xz_y_y, g_y_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xz_y_x[i] = 2.0 * g_y_xz_y_x[i] * a_exp;

        g_y_0_0_0_0_xz_y_y[i] = 2.0 * g_y_xz_y_y[i] * a_exp;

        g_y_0_0_0_0_xz_y_z[i] = 2.0 * g_y_xz_y_z[i] * a_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_y_0_0_0_0_xz_z_x, g_y_0_0_0_0_xz_z_y, g_y_0_0_0_0_xz_z_z, g_y_xz_z_x, g_y_xz_z_y, g_y_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xz_z_x[i] = 2.0 * g_y_xz_z_x[i] * a_exp;

        g_y_0_0_0_0_xz_z_y[i] = 2.0 * g_y_xz_z_y[i] * a_exp;

        g_y_0_0_0_0_xz_z_z[i] = 2.0 * g_y_xz_z_z[i] * a_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_y_0_0_0_0_yy_x_x, g_y_0_0_0_0_yy_x_y, g_y_0_0_0_0_yy_x_z, g_y_yy_x_x, g_y_yy_x_y, g_y_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yy_x_x[i] = 2.0 * g_y_yy_x_x[i] * a_exp;

        g_y_0_0_0_0_yy_x_y[i] = 2.0 * g_y_yy_x_y[i] * a_exp;

        g_y_0_0_0_0_yy_x_z[i] = 2.0 * g_y_yy_x_z[i] * a_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_y_0_0_0_0_yy_y_x, g_y_0_0_0_0_yy_y_y, g_y_0_0_0_0_yy_y_z, g_y_yy_y_x, g_y_yy_y_y, g_y_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yy_y_x[i] = 2.0 * g_y_yy_y_x[i] * a_exp;

        g_y_0_0_0_0_yy_y_y[i] = 2.0 * g_y_yy_y_y[i] * a_exp;

        g_y_0_0_0_0_yy_y_z[i] = 2.0 * g_y_yy_y_z[i] * a_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_y_0_0_0_0_yy_z_x, g_y_0_0_0_0_yy_z_y, g_y_0_0_0_0_yy_z_z, g_y_yy_z_x, g_y_yy_z_y, g_y_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yy_z_x[i] = 2.0 * g_y_yy_z_x[i] * a_exp;

        g_y_0_0_0_0_yy_z_y[i] = 2.0 * g_y_yy_z_y[i] * a_exp;

        g_y_0_0_0_0_yy_z_z[i] = 2.0 * g_y_yy_z_z[i] * a_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_y_0_0_0_0_yz_x_x, g_y_0_0_0_0_yz_x_y, g_y_0_0_0_0_yz_x_z, g_y_yz_x_x, g_y_yz_x_y, g_y_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yz_x_x[i] = 2.0 * g_y_yz_x_x[i] * a_exp;

        g_y_0_0_0_0_yz_x_y[i] = 2.0 * g_y_yz_x_y[i] * a_exp;

        g_y_0_0_0_0_yz_x_z[i] = 2.0 * g_y_yz_x_z[i] * a_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_y_0_0_0_0_yz_y_x, g_y_0_0_0_0_yz_y_y, g_y_0_0_0_0_yz_y_z, g_y_yz_y_x, g_y_yz_y_y, g_y_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yz_y_x[i] = 2.0 * g_y_yz_y_x[i] * a_exp;

        g_y_0_0_0_0_yz_y_y[i] = 2.0 * g_y_yz_y_y[i] * a_exp;

        g_y_0_0_0_0_yz_y_z[i] = 2.0 * g_y_yz_y_z[i] * a_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_y_0_0_0_0_yz_z_x, g_y_0_0_0_0_yz_z_y, g_y_0_0_0_0_yz_z_z, g_y_yz_z_x, g_y_yz_z_y, g_y_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yz_z_x[i] = 2.0 * g_y_yz_z_x[i] * a_exp;

        g_y_0_0_0_0_yz_z_y[i] = 2.0 * g_y_yz_z_y[i] * a_exp;

        g_y_0_0_0_0_yz_z_z[i] = 2.0 * g_y_yz_z_z[i] * a_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_y_0_0_0_0_zz_x_x, g_y_0_0_0_0_zz_x_y, g_y_0_0_0_0_zz_x_z, g_y_zz_x_x, g_y_zz_x_y, g_y_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_zz_x_x[i] = 2.0 * g_y_zz_x_x[i] * a_exp;

        g_y_0_0_0_0_zz_x_y[i] = 2.0 * g_y_zz_x_y[i] * a_exp;

        g_y_0_0_0_0_zz_x_z[i] = 2.0 * g_y_zz_x_z[i] * a_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_y_0_0_0_0_zz_y_x, g_y_0_0_0_0_zz_y_y, g_y_0_0_0_0_zz_y_z, g_y_zz_y_x, g_y_zz_y_y, g_y_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_zz_y_x[i] = 2.0 * g_y_zz_y_x[i] * a_exp;

        g_y_0_0_0_0_zz_y_y[i] = 2.0 * g_y_zz_y_y[i] * a_exp;

        g_y_0_0_0_0_zz_y_z[i] = 2.0 * g_y_zz_y_z[i] * a_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_y_0_0_0_0_zz_z_x, g_y_0_0_0_0_zz_z_y, g_y_0_0_0_0_zz_z_z, g_y_zz_z_x, g_y_zz_z_y, g_y_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_zz_z_x[i] = 2.0 * g_y_zz_z_x[i] * a_exp;

        g_y_0_0_0_0_zz_z_y[i] = 2.0 * g_y_zz_z_y[i] * a_exp;

        g_y_0_0_0_0_zz_z_z[i] = 2.0 * g_y_zz_z_z[i] * a_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_z_0_0_0_0_xx_x_x, g_z_0_0_0_0_xx_x_y, g_z_0_0_0_0_xx_x_z, g_z_xx_x_x, g_z_xx_x_y, g_z_xx_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xx_x_x[i] = 2.0 * g_z_xx_x_x[i] * a_exp;

        g_z_0_0_0_0_xx_x_y[i] = 2.0 * g_z_xx_x_y[i] * a_exp;

        g_z_0_0_0_0_xx_x_z[i] = 2.0 * g_z_xx_x_z[i] * a_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_z_0_0_0_0_xx_y_x, g_z_0_0_0_0_xx_y_y, g_z_0_0_0_0_xx_y_z, g_z_xx_y_x, g_z_xx_y_y, g_z_xx_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xx_y_x[i] = 2.0 * g_z_xx_y_x[i] * a_exp;

        g_z_0_0_0_0_xx_y_y[i] = 2.0 * g_z_xx_y_y[i] * a_exp;

        g_z_0_0_0_0_xx_y_z[i] = 2.0 * g_z_xx_y_z[i] * a_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_z_0_0_0_0_xx_z_x, g_z_0_0_0_0_xx_z_y, g_z_0_0_0_0_xx_z_z, g_z_xx_z_x, g_z_xx_z_y, g_z_xx_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xx_z_x[i] = 2.0 * g_z_xx_z_x[i] * a_exp;

        g_z_0_0_0_0_xx_z_y[i] = 2.0 * g_z_xx_z_y[i] * a_exp;

        g_z_0_0_0_0_xx_z_z[i] = 2.0 * g_z_xx_z_z[i] * a_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_z_0_0_0_0_xy_x_x, g_z_0_0_0_0_xy_x_y, g_z_0_0_0_0_xy_x_z, g_z_xy_x_x, g_z_xy_x_y, g_z_xy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xy_x_x[i] = 2.0 * g_z_xy_x_x[i] * a_exp;

        g_z_0_0_0_0_xy_x_y[i] = 2.0 * g_z_xy_x_y[i] * a_exp;

        g_z_0_0_0_0_xy_x_z[i] = 2.0 * g_z_xy_x_z[i] * a_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_z_0_0_0_0_xy_y_x, g_z_0_0_0_0_xy_y_y, g_z_0_0_0_0_xy_y_z, g_z_xy_y_x, g_z_xy_y_y, g_z_xy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xy_y_x[i] = 2.0 * g_z_xy_y_x[i] * a_exp;

        g_z_0_0_0_0_xy_y_y[i] = 2.0 * g_z_xy_y_y[i] * a_exp;

        g_z_0_0_0_0_xy_y_z[i] = 2.0 * g_z_xy_y_z[i] * a_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_z_0_0_0_0_xy_z_x, g_z_0_0_0_0_xy_z_y, g_z_0_0_0_0_xy_z_z, g_z_xy_z_x, g_z_xy_z_y, g_z_xy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xy_z_x[i] = 2.0 * g_z_xy_z_x[i] * a_exp;

        g_z_0_0_0_0_xy_z_y[i] = 2.0 * g_z_xy_z_y[i] * a_exp;

        g_z_0_0_0_0_xy_z_z[i] = 2.0 * g_z_xy_z_z[i] * a_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_z_0_0_0_0_xz_x_x, g_z_0_0_0_0_xz_x_y, g_z_0_0_0_0_xz_x_z, g_z_xz_x_x, g_z_xz_x_y, g_z_xz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xz_x_x[i] = 2.0 * g_z_xz_x_x[i] * a_exp;

        g_z_0_0_0_0_xz_x_y[i] = 2.0 * g_z_xz_x_y[i] * a_exp;

        g_z_0_0_0_0_xz_x_z[i] = 2.0 * g_z_xz_x_z[i] * a_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_z_0_0_0_0_xz_y_x, g_z_0_0_0_0_xz_y_y, g_z_0_0_0_0_xz_y_z, g_z_xz_y_x, g_z_xz_y_y, g_z_xz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xz_y_x[i] = 2.0 * g_z_xz_y_x[i] * a_exp;

        g_z_0_0_0_0_xz_y_y[i] = 2.0 * g_z_xz_y_y[i] * a_exp;

        g_z_0_0_0_0_xz_y_z[i] = 2.0 * g_z_xz_y_z[i] * a_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_z_0_0_0_0_xz_z_x, g_z_0_0_0_0_xz_z_y, g_z_0_0_0_0_xz_z_z, g_z_xz_z_x, g_z_xz_z_y, g_z_xz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xz_z_x[i] = 2.0 * g_z_xz_z_x[i] * a_exp;

        g_z_0_0_0_0_xz_z_y[i] = 2.0 * g_z_xz_z_y[i] * a_exp;

        g_z_0_0_0_0_xz_z_z[i] = 2.0 * g_z_xz_z_z[i] * a_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_z_0_0_0_0_yy_x_x, g_z_0_0_0_0_yy_x_y, g_z_0_0_0_0_yy_x_z, g_z_yy_x_x, g_z_yy_x_y, g_z_yy_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yy_x_x[i] = 2.0 * g_z_yy_x_x[i] * a_exp;

        g_z_0_0_0_0_yy_x_y[i] = 2.0 * g_z_yy_x_y[i] * a_exp;

        g_z_0_0_0_0_yy_x_z[i] = 2.0 * g_z_yy_x_z[i] * a_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_z_0_0_0_0_yy_y_x, g_z_0_0_0_0_yy_y_y, g_z_0_0_0_0_yy_y_z, g_z_yy_y_x, g_z_yy_y_y, g_z_yy_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yy_y_x[i] = 2.0 * g_z_yy_y_x[i] * a_exp;

        g_z_0_0_0_0_yy_y_y[i] = 2.0 * g_z_yy_y_y[i] * a_exp;

        g_z_0_0_0_0_yy_y_z[i] = 2.0 * g_z_yy_y_z[i] * a_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_z_0_0_0_0_yy_z_x, g_z_0_0_0_0_yy_z_y, g_z_0_0_0_0_yy_z_z, g_z_yy_z_x, g_z_yy_z_y, g_z_yy_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yy_z_x[i] = 2.0 * g_z_yy_z_x[i] * a_exp;

        g_z_0_0_0_0_yy_z_y[i] = 2.0 * g_z_yy_z_y[i] * a_exp;

        g_z_0_0_0_0_yy_z_z[i] = 2.0 * g_z_yy_z_z[i] * a_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_z_0_0_0_0_yz_x_x, g_z_0_0_0_0_yz_x_y, g_z_0_0_0_0_yz_x_z, g_z_yz_x_x, g_z_yz_x_y, g_z_yz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yz_x_x[i] = 2.0 * g_z_yz_x_x[i] * a_exp;

        g_z_0_0_0_0_yz_x_y[i] = 2.0 * g_z_yz_x_y[i] * a_exp;

        g_z_0_0_0_0_yz_x_z[i] = 2.0 * g_z_yz_x_z[i] * a_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_z_0_0_0_0_yz_y_x, g_z_0_0_0_0_yz_y_y, g_z_0_0_0_0_yz_y_z, g_z_yz_y_x, g_z_yz_y_y, g_z_yz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yz_y_x[i] = 2.0 * g_z_yz_y_x[i] * a_exp;

        g_z_0_0_0_0_yz_y_y[i] = 2.0 * g_z_yz_y_y[i] * a_exp;

        g_z_0_0_0_0_yz_y_z[i] = 2.0 * g_z_yz_y_z[i] * a_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_z_0_0_0_0_yz_z_x, g_z_0_0_0_0_yz_z_y, g_z_0_0_0_0_yz_z_z, g_z_yz_z_x, g_z_yz_z_y, g_z_yz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yz_z_x[i] = 2.0 * g_z_yz_z_x[i] * a_exp;

        g_z_0_0_0_0_yz_z_y[i] = 2.0 * g_z_yz_z_y[i] * a_exp;

        g_z_0_0_0_0_yz_z_z[i] = 2.0 * g_z_yz_z_z[i] * a_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_z_0_0_0_0_zz_x_x, g_z_0_0_0_0_zz_x_y, g_z_0_0_0_0_zz_x_z, g_z_zz_x_x, g_z_zz_x_y, g_z_zz_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_zz_x_x[i] = 2.0 * g_z_zz_x_x[i] * a_exp;

        g_z_0_0_0_0_zz_x_y[i] = 2.0 * g_z_zz_x_y[i] * a_exp;

        g_z_0_0_0_0_zz_x_z[i] = 2.0 * g_z_zz_x_z[i] * a_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_z_0_0_0_0_zz_y_x, g_z_0_0_0_0_zz_y_y, g_z_0_0_0_0_zz_y_z, g_z_zz_y_x, g_z_zz_y_y, g_z_zz_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_zz_y_x[i] = 2.0 * g_z_zz_y_x[i] * a_exp;

        g_z_0_0_0_0_zz_y_y[i] = 2.0 * g_z_zz_y_y[i] * a_exp;

        g_z_0_0_0_0_zz_y_z[i] = 2.0 * g_z_zz_y_z[i] * a_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_z_0_0_0_0_zz_z_x, g_z_0_0_0_0_zz_z_y, g_z_0_0_0_0_zz_z_z, g_z_zz_z_x, g_z_zz_z_y, g_z_zz_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_zz_z_x[i] = 2.0 * g_z_zz_z_x[i] * a_exp;

        g_z_0_0_0_0_zz_z_y[i] = 2.0 * g_z_zz_z_y[i] * a_exp;

        g_z_0_0_0_0_zz_z_z[i] = 2.0 * g_z_zz_z_z[i] * a_exp;
    }
}

} // t4c_geom namespace

