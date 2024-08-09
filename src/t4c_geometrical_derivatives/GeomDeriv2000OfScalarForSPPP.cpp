#include "GeomDeriv2000OfScalarForSPPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_sppp_0(CSimdArray<double>& buffer_2000_sppp,
                     const CSimdArray<double>& buffer_sppp,
                     const CSimdArray<double>& buffer_dppp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_sppp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sppp

    auto g_0_x_x_x = buffer_sppp[0];

    auto g_0_x_x_y = buffer_sppp[1];

    auto g_0_x_x_z = buffer_sppp[2];

    auto g_0_x_y_x = buffer_sppp[3];

    auto g_0_x_y_y = buffer_sppp[4];

    auto g_0_x_y_z = buffer_sppp[5];

    auto g_0_x_z_x = buffer_sppp[6];

    auto g_0_x_z_y = buffer_sppp[7];

    auto g_0_x_z_z = buffer_sppp[8];

    auto g_0_y_x_x = buffer_sppp[9];

    auto g_0_y_x_y = buffer_sppp[10];

    auto g_0_y_x_z = buffer_sppp[11];

    auto g_0_y_y_x = buffer_sppp[12];

    auto g_0_y_y_y = buffer_sppp[13];

    auto g_0_y_y_z = buffer_sppp[14];

    auto g_0_y_z_x = buffer_sppp[15];

    auto g_0_y_z_y = buffer_sppp[16];

    auto g_0_y_z_z = buffer_sppp[17];

    auto g_0_z_x_x = buffer_sppp[18];

    auto g_0_z_x_y = buffer_sppp[19];

    auto g_0_z_x_z = buffer_sppp[20];

    auto g_0_z_y_x = buffer_sppp[21];

    auto g_0_z_y_y = buffer_sppp[22];

    auto g_0_z_y_z = buffer_sppp[23];

    auto g_0_z_z_x = buffer_sppp[24];

    auto g_0_z_z_y = buffer_sppp[25];

    auto g_0_z_z_z = buffer_sppp[26];

    /// Set up components of auxilary buffer : buffer_dppp

    auto g_xx_x_x_x = buffer_dppp[0];

    auto g_xx_x_x_y = buffer_dppp[1];

    auto g_xx_x_x_z = buffer_dppp[2];

    auto g_xx_x_y_x = buffer_dppp[3];

    auto g_xx_x_y_y = buffer_dppp[4];

    auto g_xx_x_y_z = buffer_dppp[5];

    auto g_xx_x_z_x = buffer_dppp[6];

    auto g_xx_x_z_y = buffer_dppp[7];

    auto g_xx_x_z_z = buffer_dppp[8];

    auto g_xx_y_x_x = buffer_dppp[9];

    auto g_xx_y_x_y = buffer_dppp[10];

    auto g_xx_y_x_z = buffer_dppp[11];

    auto g_xx_y_y_x = buffer_dppp[12];

    auto g_xx_y_y_y = buffer_dppp[13];

    auto g_xx_y_y_z = buffer_dppp[14];

    auto g_xx_y_z_x = buffer_dppp[15];

    auto g_xx_y_z_y = buffer_dppp[16];

    auto g_xx_y_z_z = buffer_dppp[17];

    auto g_xx_z_x_x = buffer_dppp[18];

    auto g_xx_z_x_y = buffer_dppp[19];

    auto g_xx_z_x_z = buffer_dppp[20];

    auto g_xx_z_y_x = buffer_dppp[21];

    auto g_xx_z_y_y = buffer_dppp[22];

    auto g_xx_z_y_z = buffer_dppp[23];

    auto g_xx_z_z_x = buffer_dppp[24];

    auto g_xx_z_z_y = buffer_dppp[25];

    auto g_xx_z_z_z = buffer_dppp[26];

    auto g_xy_x_x_x = buffer_dppp[27];

    auto g_xy_x_x_y = buffer_dppp[28];

    auto g_xy_x_x_z = buffer_dppp[29];

    auto g_xy_x_y_x = buffer_dppp[30];

    auto g_xy_x_y_y = buffer_dppp[31];

    auto g_xy_x_y_z = buffer_dppp[32];

    auto g_xy_x_z_x = buffer_dppp[33];

    auto g_xy_x_z_y = buffer_dppp[34];

    auto g_xy_x_z_z = buffer_dppp[35];

    auto g_xy_y_x_x = buffer_dppp[36];

    auto g_xy_y_x_y = buffer_dppp[37];

    auto g_xy_y_x_z = buffer_dppp[38];

    auto g_xy_y_y_x = buffer_dppp[39];

    auto g_xy_y_y_y = buffer_dppp[40];

    auto g_xy_y_y_z = buffer_dppp[41];

    auto g_xy_y_z_x = buffer_dppp[42];

    auto g_xy_y_z_y = buffer_dppp[43];

    auto g_xy_y_z_z = buffer_dppp[44];

    auto g_xy_z_x_x = buffer_dppp[45];

    auto g_xy_z_x_y = buffer_dppp[46];

    auto g_xy_z_x_z = buffer_dppp[47];

    auto g_xy_z_y_x = buffer_dppp[48];

    auto g_xy_z_y_y = buffer_dppp[49];

    auto g_xy_z_y_z = buffer_dppp[50];

    auto g_xy_z_z_x = buffer_dppp[51];

    auto g_xy_z_z_y = buffer_dppp[52];

    auto g_xy_z_z_z = buffer_dppp[53];

    auto g_xz_x_x_x = buffer_dppp[54];

    auto g_xz_x_x_y = buffer_dppp[55];

    auto g_xz_x_x_z = buffer_dppp[56];

    auto g_xz_x_y_x = buffer_dppp[57];

    auto g_xz_x_y_y = buffer_dppp[58];

    auto g_xz_x_y_z = buffer_dppp[59];

    auto g_xz_x_z_x = buffer_dppp[60];

    auto g_xz_x_z_y = buffer_dppp[61];

    auto g_xz_x_z_z = buffer_dppp[62];

    auto g_xz_y_x_x = buffer_dppp[63];

    auto g_xz_y_x_y = buffer_dppp[64];

    auto g_xz_y_x_z = buffer_dppp[65];

    auto g_xz_y_y_x = buffer_dppp[66];

    auto g_xz_y_y_y = buffer_dppp[67];

    auto g_xz_y_y_z = buffer_dppp[68];

    auto g_xz_y_z_x = buffer_dppp[69];

    auto g_xz_y_z_y = buffer_dppp[70];

    auto g_xz_y_z_z = buffer_dppp[71];

    auto g_xz_z_x_x = buffer_dppp[72];

    auto g_xz_z_x_y = buffer_dppp[73];

    auto g_xz_z_x_z = buffer_dppp[74];

    auto g_xz_z_y_x = buffer_dppp[75];

    auto g_xz_z_y_y = buffer_dppp[76];

    auto g_xz_z_y_z = buffer_dppp[77];

    auto g_xz_z_z_x = buffer_dppp[78];

    auto g_xz_z_z_y = buffer_dppp[79];

    auto g_xz_z_z_z = buffer_dppp[80];

    auto g_yy_x_x_x = buffer_dppp[81];

    auto g_yy_x_x_y = buffer_dppp[82];

    auto g_yy_x_x_z = buffer_dppp[83];

    auto g_yy_x_y_x = buffer_dppp[84];

    auto g_yy_x_y_y = buffer_dppp[85];

    auto g_yy_x_y_z = buffer_dppp[86];

    auto g_yy_x_z_x = buffer_dppp[87];

    auto g_yy_x_z_y = buffer_dppp[88];

    auto g_yy_x_z_z = buffer_dppp[89];

    auto g_yy_y_x_x = buffer_dppp[90];

    auto g_yy_y_x_y = buffer_dppp[91];

    auto g_yy_y_x_z = buffer_dppp[92];

    auto g_yy_y_y_x = buffer_dppp[93];

    auto g_yy_y_y_y = buffer_dppp[94];

    auto g_yy_y_y_z = buffer_dppp[95];

    auto g_yy_y_z_x = buffer_dppp[96];

    auto g_yy_y_z_y = buffer_dppp[97];

    auto g_yy_y_z_z = buffer_dppp[98];

    auto g_yy_z_x_x = buffer_dppp[99];

    auto g_yy_z_x_y = buffer_dppp[100];

    auto g_yy_z_x_z = buffer_dppp[101];

    auto g_yy_z_y_x = buffer_dppp[102];

    auto g_yy_z_y_y = buffer_dppp[103];

    auto g_yy_z_y_z = buffer_dppp[104];

    auto g_yy_z_z_x = buffer_dppp[105];

    auto g_yy_z_z_y = buffer_dppp[106];

    auto g_yy_z_z_z = buffer_dppp[107];

    auto g_yz_x_x_x = buffer_dppp[108];

    auto g_yz_x_x_y = buffer_dppp[109];

    auto g_yz_x_x_z = buffer_dppp[110];

    auto g_yz_x_y_x = buffer_dppp[111];

    auto g_yz_x_y_y = buffer_dppp[112];

    auto g_yz_x_y_z = buffer_dppp[113];

    auto g_yz_x_z_x = buffer_dppp[114];

    auto g_yz_x_z_y = buffer_dppp[115];

    auto g_yz_x_z_z = buffer_dppp[116];

    auto g_yz_y_x_x = buffer_dppp[117];

    auto g_yz_y_x_y = buffer_dppp[118];

    auto g_yz_y_x_z = buffer_dppp[119];

    auto g_yz_y_y_x = buffer_dppp[120];

    auto g_yz_y_y_y = buffer_dppp[121];

    auto g_yz_y_y_z = buffer_dppp[122];

    auto g_yz_y_z_x = buffer_dppp[123];

    auto g_yz_y_z_y = buffer_dppp[124];

    auto g_yz_y_z_z = buffer_dppp[125];

    auto g_yz_z_x_x = buffer_dppp[126];

    auto g_yz_z_x_y = buffer_dppp[127];

    auto g_yz_z_x_z = buffer_dppp[128];

    auto g_yz_z_y_x = buffer_dppp[129];

    auto g_yz_z_y_y = buffer_dppp[130];

    auto g_yz_z_y_z = buffer_dppp[131];

    auto g_yz_z_z_x = buffer_dppp[132];

    auto g_yz_z_z_y = buffer_dppp[133];

    auto g_yz_z_z_z = buffer_dppp[134];

    auto g_zz_x_x_x = buffer_dppp[135];

    auto g_zz_x_x_y = buffer_dppp[136];

    auto g_zz_x_x_z = buffer_dppp[137];

    auto g_zz_x_y_x = buffer_dppp[138];

    auto g_zz_x_y_y = buffer_dppp[139];

    auto g_zz_x_y_z = buffer_dppp[140];

    auto g_zz_x_z_x = buffer_dppp[141];

    auto g_zz_x_z_y = buffer_dppp[142];

    auto g_zz_x_z_z = buffer_dppp[143];

    auto g_zz_y_x_x = buffer_dppp[144];

    auto g_zz_y_x_y = buffer_dppp[145];

    auto g_zz_y_x_z = buffer_dppp[146];

    auto g_zz_y_y_x = buffer_dppp[147];

    auto g_zz_y_y_y = buffer_dppp[148];

    auto g_zz_y_y_z = buffer_dppp[149];

    auto g_zz_y_z_x = buffer_dppp[150];

    auto g_zz_y_z_y = buffer_dppp[151];

    auto g_zz_y_z_z = buffer_dppp[152];

    auto g_zz_z_x_x = buffer_dppp[153];

    auto g_zz_z_x_y = buffer_dppp[154];

    auto g_zz_z_x_z = buffer_dppp[155];

    auto g_zz_z_y_x = buffer_dppp[156];

    auto g_zz_z_y_y = buffer_dppp[157];

    auto g_zz_z_y_z = buffer_dppp[158];

    auto g_zz_z_z_x = buffer_dppp[159];

    auto g_zz_z_z_y = buffer_dppp[160];

    auto g_zz_z_z_z = buffer_dppp[161];

    /// Set up components of integrals buffer : buffer_2000_sppp

    auto g_xx_0_0_0_0_x_x_x = buffer_2000_sppp[0];

    auto g_xx_0_0_0_0_x_x_y = buffer_2000_sppp[1];

    auto g_xx_0_0_0_0_x_x_z = buffer_2000_sppp[2];

    auto g_xx_0_0_0_0_x_y_x = buffer_2000_sppp[3];

    auto g_xx_0_0_0_0_x_y_y = buffer_2000_sppp[4];

    auto g_xx_0_0_0_0_x_y_z = buffer_2000_sppp[5];

    auto g_xx_0_0_0_0_x_z_x = buffer_2000_sppp[6];

    auto g_xx_0_0_0_0_x_z_y = buffer_2000_sppp[7];

    auto g_xx_0_0_0_0_x_z_z = buffer_2000_sppp[8];

    auto g_xx_0_0_0_0_y_x_x = buffer_2000_sppp[9];

    auto g_xx_0_0_0_0_y_x_y = buffer_2000_sppp[10];

    auto g_xx_0_0_0_0_y_x_z = buffer_2000_sppp[11];

    auto g_xx_0_0_0_0_y_y_x = buffer_2000_sppp[12];

    auto g_xx_0_0_0_0_y_y_y = buffer_2000_sppp[13];

    auto g_xx_0_0_0_0_y_y_z = buffer_2000_sppp[14];

    auto g_xx_0_0_0_0_y_z_x = buffer_2000_sppp[15];

    auto g_xx_0_0_0_0_y_z_y = buffer_2000_sppp[16];

    auto g_xx_0_0_0_0_y_z_z = buffer_2000_sppp[17];

    auto g_xx_0_0_0_0_z_x_x = buffer_2000_sppp[18];

    auto g_xx_0_0_0_0_z_x_y = buffer_2000_sppp[19];

    auto g_xx_0_0_0_0_z_x_z = buffer_2000_sppp[20];

    auto g_xx_0_0_0_0_z_y_x = buffer_2000_sppp[21];

    auto g_xx_0_0_0_0_z_y_y = buffer_2000_sppp[22];

    auto g_xx_0_0_0_0_z_y_z = buffer_2000_sppp[23];

    auto g_xx_0_0_0_0_z_z_x = buffer_2000_sppp[24];

    auto g_xx_0_0_0_0_z_z_y = buffer_2000_sppp[25];

    auto g_xx_0_0_0_0_z_z_z = buffer_2000_sppp[26];

    auto g_xy_0_0_0_0_x_x_x = buffer_2000_sppp[27];

    auto g_xy_0_0_0_0_x_x_y = buffer_2000_sppp[28];

    auto g_xy_0_0_0_0_x_x_z = buffer_2000_sppp[29];

    auto g_xy_0_0_0_0_x_y_x = buffer_2000_sppp[30];

    auto g_xy_0_0_0_0_x_y_y = buffer_2000_sppp[31];

    auto g_xy_0_0_0_0_x_y_z = buffer_2000_sppp[32];

    auto g_xy_0_0_0_0_x_z_x = buffer_2000_sppp[33];

    auto g_xy_0_0_0_0_x_z_y = buffer_2000_sppp[34];

    auto g_xy_0_0_0_0_x_z_z = buffer_2000_sppp[35];

    auto g_xy_0_0_0_0_y_x_x = buffer_2000_sppp[36];

    auto g_xy_0_0_0_0_y_x_y = buffer_2000_sppp[37];

    auto g_xy_0_0_0_0_y_x_z = buffer_2000_sppp[38];

    auto g_xy_0_0_0_0_y_y_x = buffer_2000_sppp[39];

    auto g_xy_0_0_0_0_y_y_y = buffer_2000_sppp[40];

    auto g_xy_0_0_0_0_y_y_z = buffer_2000_sppp[41];

    auto g_xy_0_0_0_0_y_z_x = buffer_2000_sppp[42];

    auto g_xy_0_0_0_0_y_z_y = buffer_2000_sppp[43];

    auto g_xy_0_0_0_0_y_z_z = buffer_2000_sppp[44];

    auto g_xy_0_0_0_0_z_x_x = buffer_2000_sppp[45];

    auto g_xy_0_0_0_0_z_x_y = buffer_2000_sppp[46];

    auto g_xy_0_0_0_0_z_x_z = buffer_2000_sppp[47];

    auto g_xy_0_0_0_0_z_y_x = buffer_2000_sppp[48];

    auto g_xy_0_0_0_0_z_y_y = buffer_2000_sppp[49];

    auto g_xy_0_0_0_0_z_y_z = buffer_2000_sppp[50];

    auto g_xy_0_0_0_0_z_z_x = buffer_2000_sppp[51];

    auto g_xy_0_0_0_0_z_z_y = buffer_2000_sppp[52];

    auto g_xy_0_0_0_0_z_z_z = buffer_2000_sppp[53];

    auto g_xz_0_0_0_0_x_x_x = buffer_2000_sppp[54];

    auto g_xz_0_0_0_0_x_x_y = buffer_2000_sppp[55];

    auto g_xz_0_0_0_0_x_x_z = buffer_2000_sppp[56];

    auto g_xz_0_0_0_0_x_y_x = buffer_2000_sppp[57];

    auto g_xz_0_0_0_0_x_y_y = buffer_2000_sppp[58];

    auto g_xz_0_0_0_0_x_y_z = buffer_2000_sppp[59];

    auto g_xz_0_0_0_0_x_z_x = buffer_2000_sppp[60];

    auto g_xz_0_0_0_0_x_z_y = buffer_2000_sppp[61];

    auto g_xz_0_0_0_0_x_z_z = buffer_2000_sppp[62];

    auto g_xz_0_0_0_0_y_x_x = buffer_2000_sppp[63];

    auto g_xz_0_0_0_0_y_x_y = buffer_2000_sppp[64];

    auto g_xz_0_0_0_0_y_x_z = buffer_2000_sppp[65];

    auto g_xz_0_0_0_0_y_y_x = buffer_2000_sppp[66];

    auto g_xz_0_0_0_0_y_y_y = buffer_2000_sppp[67];

    auto g_xz_0_0_0_0_y_y_z = buffer_2000_sppp[68];

    auto g_xz_0_0_0_0_y_z_x = buffer_2000_sppp[69];

    auto g_xz_0_0_0_0_y_z_y = buffer_2000_sppp[70];

    auto g_xz_0_0_0_0_y_z_z = buffer_2000_sppp[71];

    auto g_xz_0_0_0_0_z_x_x = buffer_2000_sppp[72];

    auto g_xz_0_0_0_0_z_x_y = buffer_2000_sppp[73];

    auto g_xz_0_0_0_0_z_x_z = buffer_2000_sppp[74];

    auto g_xz_0_0_0_0_z_y_x = buffer_2000_sppp[75];

    auto g_xz_0_0_0_0_z_y_y = buffer_2000_sppp[76];

    auto g_xz_0_0_0_0_z_y_z = buffer_2000_sppp[77];

    auto g_xz_0_0_0_0_z_z_x = buffer_2000_sppp[78];

    auto g_xz_0_0_0_0_z_z_y = buffer_2000_sppp[79];

    auto g_xz_0_0_0_0_z_z_z = buffer_2000_sppp[80];

    auto g_yy_0_0_0_0_x_x_x = buffer_2000_sppp[81];

    auto g_yy_0_0_0_0_x_x_y = buffer_2000_sppp[82];

    auto g_yy_0_0_0_0_x_x_z = buffer_2000_sppp[83];

    auto g_yy_0_0_0_0_x_y_x = buffer_2000_sppp[84];

    auto g_yy_0_0_0_0_x_y_y = buffer_2000_sppp[85];

    auto g_yy_0_0_0_0_x_y_z = buffer_2000_sppp[86];

    auto g_yy_0_0_0_0_x_z_x = buffer_2000_sppp[87];

    auto g_yy_0_0_0_0_x_z_y = buffer_2000_sppp[88];

    auto g_yy_0_0_0_0_x_z_z = buffer_2000_sppp[89];

    auto g_yy_0_0_0_0_y_x_x = buffer_2000_sppp[90];

    auto g_yy_0_0_0_0_y_x_y = buffer_2000_sppp[91];

    auto g_yy_0_0_0_0_y_x_z = buffer_2000_sppp[92];

    auto g_yy_0_0_0_0_y_y_x = buffer_2000_sppp[93];

    auto g_yy_0_0_0_0_y_y_y = buffer_2000_sppp[94];

    auto g_yy_0_0_0_0_y_y_z = buffer_2000_sppp[95];

    auto g_yy_0_0_0_0_y_z_x = buffer_2000_sppp[96];

    auto g_yy_0_0_0_0_y_z_y = buffer_2000_sppp[97];

    auto g_yy_0_0_0_0_y_z_z = buffer_2000_sppp[98];

    auto g_yy_0_0_0_0_z_x_x = buffer_2000_sppp[99];

    auto g_yy_0_0_0_0_z_x_y = buffer_2000_sppp[100];

    auto g_yy_0_0_0_0_z_x_z = buffer_2000_sppp[101];

    auto g_yy_0_0_0_0_z_y_x = buffer_2000_sppp[102];

    auto g_yy_0_0_0_0_z_y_y = buffer_2000_sppp[103];

    auto g_yy_0_0_0_0_z_y_z = buffer_2000_sppp[104];

    auto g_yy_0_0_0_0_z_z_x = buffer_2000_sppp[105];

    auto g_yy_0_0_0_0_z_z_y = buffer_2000_sppp[106];

    auto g_yy_0_0_0_0_z_z_z = buffer_2000_sppp[107];

    auto g_yz_0_0_0_0_x_x_x = buffer_2000_sppp[108];

    auto g_yz_0_0_0_0_x_x_y = buffer_2000_sppp[109];

    auto g_yz_0_0_0_0_x_x_z = buffer_2000_sppp[110];

    auto g_yz_0_0_0_0_x_y_x = buffer_2000_sppp[111];

    auto g_yz_0_0_0_0_x_y_y = buffer_2000_sppp[112];

    auto g_yz_0_0_0_0_x_y_z = buffer_2000_sppp[113];

    auto g_yz_0_0_0_0_x_z_x = buffer_2000_sppp[114];

    auto g_yz_0_0_0_0_x_z_y = buffer_2000_sppp[115];

    auto g_yz_0_0_0_0_x_z_z = buffer_2000_sppp[116];

    auto g_yz_0_0_0_0_y_x_x = buffer_2000_sppp[117];

    auto g_yz_0_0_0_0_y_x_y = buffer_2000_sppp[118];

    auto g_yz_0_0_0_0_y_x_z = buffer_2000_sppp[119];

    auto g_yz_0_0_0_0_y_y_x = buffer_2000_sppp[120];

    auto g_yz_0_0_0_0_y_y_y = buffer_2000_sppp[121];

    auto g_yz_0_0_0_0_y_y_z = buffer_2000_sppp[122];

    auto g_yz_0_0_0_0_y_z_x = buffer_2000_sppp[123];

    auto g_yz_0_0_0_0_y_z_y = buffer_2000_sppp[124];

    auto g_yz_0_0_0_0_y_z_z = buffer_2000_sppp[125];

    auto g_yz_0_0_0_0_z_x_x = buffer_2000_sppp[126];

    auto g_yz_0_0_0_0_z_x_y = buffer_2000_sppp[127];

    auto g_yz_0_0_0_0_z_x_z = buffer_2000_sppp[128];

    auto g_yz_0_0_0_0_z_y_x = buffer_2000_sppp[129];

    auto g_yz_0_0_0_0_z_y_y = buffer_2000_sppp[130];

    auto g_yz_0_0_0_0_z_y_z = buffer_2000_sppp[131];

    auto g_yz_0_0_0_0_z_z_x = buffer_2000_sppp[132];

    auto g_yz_0_0_0_0_z_z_y = buffer_2000_sppp[133];

    auto g_yz_0_0_0_0_z_z_z = buffer_2000_sppp[134];

    auto g_zz_0_0_0_0_x_x_x = buffer_2000_sppp[135];

    auto g_zz_0_0_0_0_x_x_y = buffer_2000_sppp[136];

    auto g_zz_0_0_0_0_x_x_z = buffer_2000_sppp[137];

    auto g_zz_0_0_0_0_x_y_x = buffer_2000_sppp[138];

    auto g_zz_0_0_0_0_x_y_y = buffer_2000_sppp[139];

    auto g_zz_0_0_0_0_x_y_z = buffer_2000_sppp[140];

    auto g_zz_0_0_0_0_x_z_x = buffer_2000_sppp[141];

    auto g_zz_0_0_0_0_x_z_y = buffer_2000_sppp[142];

    auto g_zz_0_0_0_0_x_z_z = buffer_2000_sppp[143];

    auto g_zz_0_0_0_0_y_x_x = buffer_2000_sppp[144];

    auto g_zz_0_0_0_0_y_x_y = buffer_2000_sppp[145];

    auto g_zz_0_0_0_0_y_x_z = buffer_2000_sppp[146];

    auto g_zz_0_0_0_0_y_y_x = buffer_2000_sppp[147];

    auto g_zz_0_0_0_0_y_y_y = buffer_2000_sppp[148];

    auto g_zz_0_0_0_0_y_y_z = buffer_2000_sppp[149];

    auto g_zz_0_0_0_0_y_z_x = buffer_2000_sppp[150];

    auto g_zz_0_0_0_0_y_z_y = buffer_2000_sppp[151];

    auto g_zz_0_0_0_0_y_z_z = buffer_2000_sppp[152];

    auto g_zz_0_0_0_0_z_x_x = buffer_2000_sppp[153];

    auto g_zz_0_0_0_0_z_x_y = buffer_2000_sppp[154];

    auto g_zz_0_0_0_0_z_x_z = buffer_2000_sppp[155];

    auto g_zz_0_0_0_0_z_y_x = buffer_2000_sppp[156];

    auto g_zz_0_0_0_0_z_y_y = buffer_2000_sppp[157];

    auto g_zz_0_0_0_0_z_y_z = buffer_2000_sppp[158];

    auto g_zz_0_0_0_0_z_z_x = buffer_2000_sppp[159];

    auto g_zz_0_0_0_0_z_z_y = buffer_2000_sppp[160];

    auto g_zz_0_0_0_0_z_z_z = buffer_2000_sppp[161];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_y, g_0_x_x_z, g_xx_0_0_0_0_x_x_x, g_xx_0_0_0_0_x_x_y, g_xx_0_0_0_0_x_x_z, g_xx_x_x_x, g_xx_x_x_y, g_xx_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_x_x_x[i] = -2.0 * g_0_x_x_x[i] * a_exp + 4.0 * g_xx_x_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_x_y[i] = -2.0 * g_0_x_x_y[i] * a_exp + 4.0 * g_xx_x_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_x_z[i] = -2.0 * g_0_x_x_z[i] * a_exp + 4.0 * g_xx_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_0_x_y_x, g_0_x_y_y, g_0_x_y_z, g_xx_0_0_0_0_x_y_x, g_xx_0_0_0_0_x_y_y, g_xx_0_0_0_0_x_y_z, g_xx_x_y_x, g_xx_x_y_y, g_xx_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_x_y_x[i] = -2.0 * g_0_x_y_x[i] * a_exp + 4.0 * g_xx_x_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_y_y[i] = -2.0 * g_0_x_y_y[i] * a_exp + 4.0 * g_xx_x_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_y_z[i] = -2.0 * g_0_x_y_z[i] * a_exp + 4.0 * g_xx_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_0_x_z_x, g_0_x_z_y, g_0_x_z_z, g_xx_0_0_0_0_x_z_x, g_xx_0_0_0_0_x_z_y, g_xx_0_0_0_0_x_z_z, g_xx_x_z_x, g_xx_x_z_y, g_xx_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_x_z_x[i] = -2.0 * g_0_x_z_x[i] * a_exp + 4.0 * g_xx_x_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_z_y[i] = -2.0 * g_0_x_z_y[i] * a_exp + 4.0 * g_xx_x_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_z_z[i] = -2.0 * g_0_x_z_z[i] * a_exp + 4.0 * g_xx_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_0_y_x_x, g_0_y_x_y, g_0_y_x_z, g_xx_0_0_0_0_y_x_x, g_xx_0_0_0_0_y_x_y, g_xx_0_0_0_0_y_x_z, g_xx_y_x_x, g_xx_y_x_y, g_xx_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_y_x_x[i] = -2.0 * g_0_y_x_x[i] * a_exp + 4.0 * g_xx_y_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_x_y[i] = -2.0 * g_0_y_x_y[i] * a_exp + 4.0 * g_xx_y_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_x_z[i] = -2.0 * g_0_y_x_z[i] * a_exp + 4.0 * g_xx_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_0_y_y_x, g_0_y_y_y, g_0_y_y_z, g_xx_0_0_0_0_y_y_x, g_xx_0_0_0_0_y_y_y, g_xx_0_0_0_0_y_y_z, g_xx_y_y_x, g_xx_y_y_y, g_xx_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_y_y_x[i] = -2.0 * g_0_y_y_x[i] * a_exp + 4.0 * g_xx_y_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_y_y[i] = -2.0 * g_0_y_y_y[i] * a_exp + 4.0 * g_xx_y_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_y_z[i] = -2.0 * g_0_y_y_z[i] * a_exp + 4.0 * g_xx_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_0_y_z_x, g_0_y_z_y, g_0_y_z_z, g_xx_0_0_0_0_y_z_x, g_xx_0_0_0_0_y_z_y, g_xx_0_0_0_0_y_z_z, g_xx_y_z_x, g_xx_y_z_y, g_xx_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_y_z_x[i] = -2.0 * g_0_y_z_x[i] * a_exp + 4.0 * g_xx_y_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_z_y[i] = -2.0 * g_0_y_z_y[i] * a_exp + 4.0 * g_xx_y_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_z_z[i] = -2.0 * g_0_y_z_z[i] * a_exp + 4.0 * g_xx_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_0_z_x_x, g_0_z_x_y, g_0_z_x_z, g_xx_0_0_0_0_z_x_x, g_xx_0_0_0_0_z_x_y, g_xx_0_0_0_0_z_x_z, g_xx_z_x_x, g_xx_z_x_y, g_xx_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_z_x_x[i] = -2.0 * g_0_z_x_x[i] * a_exp + 4.0 * g_xx_z_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_x_y[i] = -2.0 * g_0_z_x_y[i] * a_exp + 4.0 * g_xx_z_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_x_z[i] = -2.0 * g_0_z_x_z[i] * a_exp + 4.0 * g_xx_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_0_z_y_x, g_0_z_y_y, g_0_z_y_z, g_xx_0_0_0_0_z_y_x, g_xx_0_0_0_0_z_y_y, g_xx_0_0_0_0_z_y_z, g_xx_z_y_x, g_xx_z_y_y, g_xx_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_z_y_x[i] = -2.0 * g_0_z_y_x[i] * a_exp + 4.0 * g_xx_z_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_y_y[i] = -2.0 * g_0_z_y_y[i] * a_exp + 4.0 * g_xx_z_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_y_z[i] = -2.0 * g_0_z_y_z[i] * a_exp + 4.0 * g_xx_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_xx_0_0_0_0_z_z_x, g_xx_0_0_0_0_z_z_y, g_xx_0_0_0_0_z_z_z, g_xx_z_z_x, g_xx_z_z_y, g_xx_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_z_z_x[i] = -2.0 * g_0_z_z_x[i] * a_exp + 4.0 * g_xx_z_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_z_y[i] = -2.0 * g_0_z_z_y[i] * a_exp + 4.0 * g_xx_z_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_z_z[i] = -2.0 * g_0_z_z_z[i] * a_exp + 4.0 * g_xx_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_xy_0_0_0_0_x_x_x, g_xy_0_0_0_0_x_x_y, g_xy_0_0_0_0_x_x_z, g_xy_x_x_x, g_xy_x_x_y, g_xy_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_x_x_x[i] = 4.0 * g_xy_x_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_x_y[i] = 4.0 * g_xy_x_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_x_z[i] = 4.0 * g_xy_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_xy_0_0_0_0_x_y_x, g_xy_0_0_0_0_x_y_y, g_xy_0_0_0_0_x_y_z, g_xy_x_y_x, g_xy_x_y_y, g_xy_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_x_y_x[i] = 4.0 * g_xy_x_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_y_y[i] = 4.0 * g_xy_x_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_y_z[i] = 4.0 * g_xy_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_xy_0_0_0_0_x_z_x, g_xy_0_0_0_0_x_z_y, g_xy_0_0_0_0_x_z_z, g_xy_x_z_x, g_xy_x_z_y, g_xy_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_x_z_x[i] = 4.0 * g_xy_x_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_z_y[i] = 4.0 * g_xy_x_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_z_z[i] = 4.0 * g_xy_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_xy_0_0_0_0_y_x_x, g_xy_0_0_0_0_y_x_y, g_xy_0_0_0_0_y_x_z, g_xy_y_x_x, g_xy_y_x_y, g_xy_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_y_x_x[i] = 4.0 * g_xy_y_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_x_y[i] = 4.0 * g_xy_y_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_x_z[i] = 4.0 * g_xy_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_xy_0_0_0_0_y_y_x, g_xy_0_0_0_0_y_y_y, g_xy_0_0_0_0_y_y_z, g_xy_y_y_x, g_xy_y_y_y, g_xy_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_y_y_x[i] = 4.0 * g_xy_y_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_y_y[i] = 4.0 * g_xy_y_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_y_z[i] = 4.0 * g_xy_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_xy_0_0_0_0_y_z_x, g_xy_0_0_0_0_y_z_y, g_xy_0_0_0_0_y_z_z, g_xy_y_z_x, g_xy_y_z_y, g_xy_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_y_z_x[i] = 4.0 * g_xy_y_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_z_y[i] = 4.0 * g_xy_y_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_z_z[i] = 4.0 * g_xy_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_xy_0_0_0_0_z_x_x, g_xy_0_0_0_0_z_x_y, g_xy_0_0_0_0_z_x_z, g_xy_z_x_x, g_xy_z_x_y, g_xy_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_z_x_x[i] = 4.0 * g_xy_z_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_x_y[i] = 4.0 * g_xy_z_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_x_z[i] = 4.0 * g_xy_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_xy_0_0_0_0_z_y_x, g_xy_0_0_0_0_z_y_y, g_xy_0_0_0_0_z_y_z, g_xy_z_y_x, g_xy_z_y_y, g_xy_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_z_y_x[i] = 4.0 * g_xy_z_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_y_y[i] = 4.0 * g_xy_z_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_y_z[i] = 4.0 * g_xy_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_xy_0_0_0_0_z_z_x, g_xy_0_0_0_0_z_z_y, g_xy_0_0_0_0_z_z_z, g_xy_z_z_x, g_xy_z_z_y, g_xy_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_z_z_x[i] = 4.0 * g_xy_z_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_z_y[i] = 4.0 * g_xy_z_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_z_z[i] = 4.0 * g_xy_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_xz_0_0_0_0_x_x_x, g_xz_0_0_0_0_x_x_y, g_xz_0_0_0_0_x_x_z, g_xz_x_x_x, g_xz_x_x_y, g_xz_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_x_x_x[i] = 4.0 * g_xz_x_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_x_y[i] = 4.0 * g_xz_x_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_x_z[i] = 4.0 * g_xz_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_xz_0_0_0_0_x_y_x, g_xz_0_0_0_0_x_y_y, g_xz_0_0_0_0_x_y_z, g_xz_x_y_x, g_xz_x_y_y, g_xz_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_x_y_x[i] = 4.0 * g_xz_x_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_y_y[i] = 4.0 * g_xz_x_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_y_z[i] = 4.0 * g_xz_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_xz_0_0_0_0_x_z_x, g_xz_0_0_0_0_x_z_y, g_xz_0_0_0_0_x_z_z, g_xz_x_z_x, g_xz_x_z_y, g_xz_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_x_z_x[i] = 4.0 * g_xz_x_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_z_y[i] = 4.0 * g_xz_x_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_z_z[i] = 4.0 * g_xz_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_xz_0_0_0_0_y_x_x, g_xz_0_0_0_0_y_x_y, g_xz_0_0_0_0_y_x_z, g_xz_y_x_x, g_xz_y_x_y, g_xz_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_y_x_x[i] = 4.0 * g_xz_y_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_x_y[i] = 4.0 * g_xz_y_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_x_z[i] = 4.0 * g_xz_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_xz_0_0_0_0_y_y_x, g_xz_0_0_0_0_y_y_y, g_xz_0_0_0_0_y_y_z, g_xz_y_y_x, g_xz_y_y_y, g_xz_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_y_y_x[i] = 4.0 * g_xz_y_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_y_y[i] = 4.0 * g_xz_y_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_y_z[i] = 4.0 * g_xz_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_xz_0_0_0_0_y_z_x, g_xz_0_0_0_0_y_z_y, g_xz_0_0_0_0_y_z_z, g_xz_y_z_x, g_xz_y_z_y, g_xz_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_y_z_x[i] = 4.0 * g_xz_y_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_z_y[i] = 4.0 * g_xz_y_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_z_z[i] = 4.0 * g_xz_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_xz_0_0_0_0_z_x_x, g_xz_0_0_0_0_z_x_y, g_xz_0_0_0_0_z_x_z, g_xz_z_x_x, g_xz_z_x_y, g_xz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_z_x_x[i] = 4.0 * g_xz_z_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_x_y[i] = 4.0 * g_xz_z_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_x_z[i] = 4.0 * g_xz_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_xz_0_0_0_0_z_y_x, g_xz_0_0_0_0_z_y_y, g_xz_0_0_0_0_z_y_z, g_xz_z_y_x, g_xz_z_y_y, g_xz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_z_y_x[i] = 4.0 * g_xz_z_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_y_y[i] = 4.0 * g_xz_z_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_y_z[i] = 4.0 * g_xz_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_xz_0_0_0_0_z_z_x, g_xz_0_0_0_0_z_z_y, g_xz_0_0_0_0_z_z_z, g_xz_z_z_x, g_xz_z_z_y, g_xz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_z_z_x[i] = 4.0 * g_xz_z_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_z_y[i] = 4.0 * g_xz_z_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_z_z[i] = 4.0 * g_xz_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_y, g_0_x_x_z, g_yy_0_0_0_0_x_x_x, g_yy_0_0_0_0_x_x_y, g_yy_0_0_0_0_x_x_z, g_yy_x_x_x, g_yy_x_x_y, g_yy_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_x_x_x[i] = -2.0 * g_0_x_x_x[i] * a_exp + 4.0 * g_yy_x_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_x_y[i] = -2.0 * g_0_x_x_y[i] * a_exp + 4.0 * g_yy_x_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_x_z[i] = -2.0 * g_0_x_x_z[i] * a_exp + 4.0 * g_yy_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_0_x_y_x, g_0_x_y_y, g_0_x_y_z, g_yy_0_0_0_0_x_y_x, g_yy_0_0_0_0_x_y_y, g_yy_0_0_0_0_x_y_z, g_yy_x_y_x, g_yy_x_y_y, g_yy_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_x_y_x[i] = -2.0 * g_0_x_y_x[i] * a_exp + 4.0 * g_yy_x_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_y_y[i] = -2.0 * g_0_x_y_y[i] * a_exp + 4.0 * g_yy_x_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_y_z[i] = -2.0 * g_0_x_y_z[i] * a_exp + 4.0 * g_yy_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_0_x_z_x, g_0_x_z_y, g_0_x_z_z, g_yy_0_0_0_0_x_z_x, g_yy_0_0_0_0_x_z_y, g_yy_0_0_0_0_x_z_z, g_yy_x_z_x, g_yy_x_z_y, g_yy_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_x_z_x[i] = -2.0 * g_0_x_z_x[i] * a_exp + 4.0 * g_yy_x_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_z_y[i] = -2.0 * g_0_x_z_y[i] * a_exp + 4.0 * g_yy_x_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_z_z[i] = -2.0 * g_0_x_z_z[i] * a_exp + 4.0 * g_yy_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_0_y_x_x, g_0_y_x_y, g_0_y_x_z, g_yy_0_0_0_0_y_x_x, g_yy_0_0_0_0_y_x_y, g_yy_0_0_0_0_y_x_z, g_yy_y_x_x, g_yy_y_x_y, g_yy_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_y_x_x[i] = -2.0 * g_0_y_x_x[i] * a_exp + 4.0 * g_yy_y_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_x_y[i] = -2.0 * g_0_y_x_y[i] * a_exp + 4.0 * g_yy_y_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_x_z[i] = -2.0 * g_0_y_x_z[i] * a_exp + 4.0 * g_yy_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_0_y_y_x, g_0_y_y_y, g_0_y_y_z, g_yy_0_0_0_0_y_y_x, g_yy_0_0_0_0_y_y_y, g_yy_0_0_0_0_y_y_z, g_yy_y_y_x, g_yy_y_y_y, g_yy_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_y_y_x[i] = -2.0 * g_0_y_y_x[i] * a_exp + 4.0 * g_yy_y_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_y_y[i] = -2.0 * g_0_y_y_y[i] * a_exp + 4.0 * g_yy_y_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_y_z[i] = -2.0 * g_0_y_y_z[i] * a_exp + 4.0 * g_yy_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_0_y_z_x, g_0_y_z_y, g_0_y_z_z, g_yy_0_0_0_0_y_z_x, g_yy_0_0_0_0_y_z_y, g_yy_0_0_0_0_y_z_z, g_yy_y_z_x, g_yy_y_z_y, g_yy_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_y_z_x[i] = -2.0 * g_0_y_z_x[i] * a_exp + 4.0 * g_yy_y_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_z_y[i] = -2.0 * g_0_y_z_y[i] * a_exp + 4.0 * g_yy_y_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_z_z[i] = -2.0 * g_0_y_z_z[i] * a_exp + 4.0 * g_yy_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_0_z_x_x, g_0_z_x_y, g_0_z_x_z, g_yy_0_0_0_0_z_x_x, g_yy_0_0_0_0_z_x_y, g_yy_0_0_0_0_z_x_z, g_yy_z_x_x, g_yy_z_x_y, g_yy_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_z_x_x[i] = -2.0 * g_0_z_x_x[i] * a_exp + 4.0 * g_yy_z_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_x_y[i] = -2.0 * g_0_z_x_y[i] * a_exp + 4.0 * g_yy_z_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_x_z[i] = -2.0 * g_0_z_x_z[i] * a_exp + 4.0 * g_yy_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_0_z_y_x, g_0_z_y_y, g_0_z_y_z, g_yy_0_0_0_0_z_y_x, g_yy_0_0_0_0_z_y_y, g_yy_0_0_0_0_z_y_z, g_yy_z_y_x, g_yy_z_y_y, g_yy_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_z_y_x[i] = -2.0 * g_0_z_y_x[i] * a_exp + 4.0 * g_yy_z_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_y_y[i] = -2.0 * g_0_z_y_y[i] * a_exp + 4.0 * g_yy_z_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_y_z[i] = -2.0 * g_0_z_y_z[i] * a_exp + 4.0 * g_yy_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_yy_0_0_0_0_z_z_x, g_yy_0_0_0_0_z_z_y, g_yy_0_0_0_0_z_z_z, g_yy_z_z_x, g_yy_z_z_y, g_yy_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_z_z_x[i] = -2.0 * g_0_z_z_x[i] * a_exp + 4.0 * g_yy_z_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_z_y[i] = -2.0 * g_0_z_z_y[i] * a_exp + 4.0 * g_yy_z_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_z_z[i] = -2.0 * g_0_z_z_z[i] * a_exp + 4.0 * g_yy_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_yz_0_0_0_0_x_x_x, g_yz_0_0_0_0_x_x_y, g_yz_0_0_0_0_x_x_z, g_yz_x_x_x, g_yz_x_x_y, g_yz_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_x_x_x[i] = 4.0 * g_yz_x_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_x_y[i] = 4.0 * g_yz_x_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_x_z[i] = 4.0 * g_yz_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_yz_0_0_0_0_x_y_x, g_yz_0_0_0_0_x_y_y, g_yz_0_0_0_0_x_y_z, g_yz_x_y_x, g_yz_x_y_y, g_yz_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_x_y_x[i] = 4.0 * g_yz_x_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_y_y[i] = 4.0 * g_yz_x_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_y_z[i] = 4.0 * g_yz_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_yz_0_0_0_0_x_z_x, g_yz_0_0_0_0_x_z_y, g_yz_0_0_0_0_x_z_z, g_yz_x_z_x, g_yz_x_z_y, g_yz_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_x_z_x[i] = 4.0 * g_yz_x_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_z_y[i] = 4.0 * g_yz_x_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_z_z[i] = 4.0 * g_yz_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_yz_0_0_0_0_y_x_x, g_yz_0_0_0_0_y_x_y, g_yz_0_0_0_0_y_x_z, g_yz_y_x_x, g_yz_y_x_y, g_yz_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_y_x_x[i] = 4.0 * g_yz_y_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_x_y[i] = 4.0 * g_yz_y_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_x_z[i] = 4.0 * g_yz_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_yz_0_0_0_0_y_y_x, g_yz_0_0_0_0_y_y_y, g_yz_0_0_0_0_y_y_z, g_yz_y_y_x, g_yz_y_y_y, g_yz_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_y_y_x[i] = 4.0 * g_yz_y_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_y_y[i] = 4.0 * g_yz_y_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_y_z[i] = 4.0 * g_yz_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_yz_0_0_0_0_y_z_x, g_yz_0_0_0_0_y_z_y, g_yz_0_0_0_0_y_z_z, g_yz_y_z_x, g_yz_y_z_y, g_yz_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_y_z_x[i] = 4.0 * g_yz_y_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_z_y[i] = 4.0 * g_yz_y_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_z_z[i] = 4.0 * g_yz_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_yz_0_0_0_0_z_x_x, g_yz_0_0_0_0_z_x_y, g_yz_0_0_0_0_z_x_z, g_yz_z_x_x, g_yz_z_x_y, g_yz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_z_x_x[i] = 4.0 * g_yz_z_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_x_y[i] = 4.0 * g_yz_z_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_x_z[i] = 4.0 * g_yz_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_yz_0_0_0_0_z_y_x, g_yz_0_0_0_0_z_y_y, g_yz_0_0_0_0_z_y_z, g_yz_z_y_x, g_yz_z_y_y, g_yz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_z_y_x[i] = 4.0 * g_yz_z_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_y_y[i] = 4.0 * g_yz_z_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_y_z[i] = 4.0 * g_yz_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_yz_0_0_0_0_z_z_x, g_yz_0_0_0_0_z_z_y, g_yz_0_0_0_0_z_z_z, g_yz_z_z_x, g_yz_z_z_y, g_yz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_z_z_x[i] = 4.0 * g_yz_z_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_z_y[i] = 4.0 * g_yz_z_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_z_z[i] = 4.0 * g_yz_z_z_z[i] * a_exp * a_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_y, g_0_x_x_z, g_zz_0_0_0_0_x_x_x, g_zz_0_0_0_0_x_x_y, g_zz_0_0_0_0_x_x_z, g_zz_x_x_x, g_zz_x_x_y, g_zz_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_x_x_x[i] = -2.0 * g_0_x_x_x[i] * a_exp + 4.0 * g_zz_x_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_x_y[i] = -2.0 * g_0_x_x_y[i] * a_exp + 4.0 * g_zz_x_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_x_z[i] = -2.0 * g_0_x_x_z[i] * a_exp + 4.0 * g_zz_x_x_z[i] * a_exp * a_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_0_x_y_x, g_0_x_y_y, g_0_x_y_z, g_zz_0_0_0_0_x_y_x, g_zz_0_0_0_0_x_y_y, g_zz_0_0_0_0_x_y_z, g_zz_x_y_x, g_zz_x_y_y, g_zz_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_x_y_x[i] = -2.0 * g_0_x_y_x[i] * a_exp + 4.0 * g_zz_x_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_y_y[i] = -2.0 * g_0_x_y_y[i] * a_exp + 4.0 * g_zz_x_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_y_z[i] = -2.0 * g_0_x_y_z[i] * a_exp + 4.0 * g_zz_x_y_z[i] * a_exp * a_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_0_x_z_x, g_0_x_z_y, g_0_x_z_z, g_zz_0_0_0_0_x_z_x, g_zz_0_0_0_0_x_z_y, g_zz_0_0_0_0_x_z_z, g_zz_x_z_x, g_zz_x_z_y, g_zz_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_x_z_x[i] = -2.0 * g_0_x_z_x[i] * a_exp + 4.0 * g_zz_x_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_z_y[i] = -2.0 * g_0_x_z_y[i] * a_exp + 4.0 * g_zz_x_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_z_z[i] = -2.0 * g_0_x_z_z[i] * a_exp + 4.0 * g_zz_x_z_z[i] * a_exp * a_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_0_y_x_x, g_0_y_x_y, g_0_y_x_z, g_zz_0_0_0_0_y_x_x, g_zz_0_0_0_0_y_x_y, g_zz_0_0_0_0_y_x_z, g_zz_y_x_x, g_zz_y_x_y, g_zz_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_y_x_x[i] = -2.0 * g_0_y_x_x[i] * a_exp + 4.0 * g_zz_y_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_x_y[i] = -2.0 * g_0_y_x_y[i] * a_exp + 4.0 * g_zz_y_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_x_z[i] = -2.0 * g_0_y_x_z[i] * a_exp + 4.0 * g_zz_y_x_z[i] * a_exp * a_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_0_y_y_x, g_0_y_y_y, g_0_y_y_z, g_zz_0_0_0_0_y_y_x, g_zz_0_0_0_0_y_y_y, g_zz_0_0_0_0_y_y_z, g_zz_y_y_x, g_zz_y_y_y, g_zz_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_y_y_x[i] = -2.0 * g_0_y_y_x[i] * a_exp + 4.0 * g_zz_y_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_y_y[i] = -2.0 * g_0_y_y_y[i] * a_exp + 4.0 * g_zz_y_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_y_z[i] = -2.0 * g_0_y_y_z[i] * a_exp + 4.0 * g_zz_y_y_z[i] * a_exp * a_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_0_y_z_x, g_0_y_z_y, g_0_y_z_z, g_zz_0_0_0_0_y_z_x, g_zz_0_0_0_0_y_z_y, g_zz_0_0_0_0_y_z_z, g_zz_y_z_x, g_zz_y_z_y, g_zz_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_y_z_x[i] = -2.0 * g_0_y_z_x[i] * a_exp + 4.0 * g_zz_y_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_z_y[i] = -2.0 * g_0_y_z_y[i] * a_exp + 4.0 * g_zz_y_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_z_z[i] = -2.0 * g_0_y_z_z[i] * a_exp + 4.0 * g_zz_y_z_z[i] * a_exp * a_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_0_z_x_x, g_0_z_x_y, g_0_z_x_z, g_zz_0_0_0_0_z_x_x, g_zz_0_0_0_0_z_x_y, g_zz_0_0_0_0_z_x_z, g_zz_z_x_x, g_zz_z_x_y, g_zz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_z_x_x[i] = -2.0 * g_0_z_x_x[i] * a_exp + 4.0 * g_zz_z_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_x_y[i] = -2.0 * g_0_z_x_y[i] * a_exp + 4.0 * g_zz_z_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_x_z[i] = -2.0 * g_0_z_x_z[i] * a_exp + 4.0 * g_zz_z_x_z[i] * a_exp * a_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_0_z_y_x, g_0_z_y_y, g_0_z_y_z, g_zz_0_0_0_0_z_y_x, g_zz_0_0_0_0_z_y_y, g_zz_0_0_0_0_z_y_z, g_zz_z_y_x, g_zz_z_y_y, g_zz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_z_y_x[i] = -2.0 * g_0_z_y_x[i] * a_exp + 4.0 * g_zz_z_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_y_y[i] = -2.0 * g_0_z_y_y[i] * a_exp + 4.0 * g_zz_z_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_y_z[i] = -2.0 * g_0_z_y_z[i] * a_exp + 4.0 * g_zz_z_y_z[i] * a_exp * a_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_zz_0_0_0_0_z_z_x, g_zz_0_0_0_0_z_z_y, g_zz_0_0_0_0_z_z_z, g_zz_z_z_x, g_zz_z_z_y, g_zz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_z_z_x[i] = -2.0 * g_0_z_z_x[i] * a_exp + 4.0 * g_zz_z_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_z_y[i] = -2.0 * g_0_z_z_y[i] * a_exp + 4.0 * g_zz_z_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_z_z[i] = -2.0 * g_0_z_z_z[i] * a_exp + 4.0 * g_zz_z_z_z[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

