#include "GeomDeriv1000OfScalarForPPPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_pppp_0(CSimdArray<double>& buffer_1000_pppp,
                     const CSimdArray<double>& buffer_sppp,
                     const CSimdArray<double>& buffer_dppp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_pppp.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_pppp

    auto g_x_0_0_0_x_x_x_x = buffer_1000_pppp[0];

    auto g_x_0_0_0_x_x_x_y = buffer_1000_pppp[1];

    auto g_x_0_0_0_x_x_x_z = buffer_1000_pppp[2];

    auto g_x_0_0_0_x_x_y_x = buffer_1000_pppp[3];

    auto g_x_0_0_0_x_x_y_y = buffer_1000_pppp[4];

    auto g_x_0_0_0_x_x_y_z = buffer_1000_pppp[5];

    auto g_x_0_0_0_x_x_z_x = buffer_1000_pppp[6];

    auto g_x_0_0_0_x_x_z_y = buffer_1000_pppp[7];

    auto g_x_0_0_0_x_x_z_z = buffer_1000_pppp[8];

    auto g_x_0_0_0_x_y_x_x = buffer_1000_pppp[9];

    auto g_x_0_0_0_x_y_x_y = buffer_1000_pppp[10];

    auto g_x_0_0_0_x_y_x_z = buffer_1000_pppp[11];

    auto g_x_0_0_0_x_y_y_x = buffer_1000_pppp[12];

    auto g_x_0_0_0_x_y_y_y = buffer_1000_pppp[13];

    auto g_x_0_0_0_x_y_y_z = buffer_1000_pppp[14];

    auto g_x_0_0_0_x_y_z_x = buffer_1000_pppp[15];

    auto g_x_0_0_0_x_y_z_y = buffer_1000_pppp[16];

    auto g_x_0_0_0_x_y_z_z = buffer_1000_pppp[17];

    auto g_x_0_0_0_x_z_x_x = buffer_1000_pppp[18];

    auto g_x_0_0_0_x_z_x_y = buffer_1000_pppp[19];

    auto g_x_0_0_0_x_z_x_z = buffer_1000_pppp[20];

    auto g_x_0_0_0_x_z_y_x = buffer_1000_pppp[21];

    auto g_x_0_0_0_x_z_y_y = buffer_1000_pppp[22];

    auto g_x_0_0_0_x_z_y_z = buffer_1000_pppp[23];

    auto g_x_0_0_0_x_z_z_x = buffer_1000_pppp[24];

    auto g_x_0_0_0_x_z_z_y = buffer_1000_pppp[25];

    auto g_x_0_0_0_x_z_z_z = buffer_1000_pppp[26];

    auto g_x_0_0_0_y_x_x_x = buffer_1000_pppp[27];

    auto g_x_0_0_0_y_x_x_y = buffer_1000_pppp[28];

    auto g_x_0_0_0_y_x_x_z = buffer_1000_pppp[29];

    auto g_x_0_0_0_y_x_y_x = buffer_1000_pppp[30];

    auto g_x_0_0_0_y_x_y_y = buffer_1000_pppp[31];

    auto g_x_0_0_0_y_x_y_z = buffer_1000_pppp[32];

    auto g_x_0_0_0_y_x_z_x = buffer_1000_pppp[33];

    auto g_x_0_0_0_y_x_z_y = buffer_1000_pppp[34];

    auto g_x_0_0_0_y_x_z_z = buffer_1000_pppp[35];

    auto g_x_0_0_0_y_y_x_x = buffer_1000_pppp[36];

    auto g_x_0_0_0_y_y_x_y = buffer_1000_pppp[37];

    auto g_x_0_0_0_y_y_x_z = buffer_1000_pppp[38];

    auto g_x_0_0_0_y_y_y_x = buffer_1000_pppp[39];

    auto g_x_0_0_0_y_y_y_y = buffer_1000_pppp[40];

    auto g_x_0_0_0_y_y_y_z = buffer_1000_pppp[41];

    auto g_x_0_0_0_y_y_z_x = buffer_1000_pppp[42];

    auto g_x_0_0_0_y_y_z_y = buffer_1000_pppp[43];

    auto g_x_0_0_0_y_y_z_z = buffer_1000_pppp[44];

    auto g_x_0_0_0_y_z_x_x = buffer_1000_pppp[45];

    auto g_x_0_0_0_y_z_x_y = buffer_1000_pppp[46];

    auto g_x_0_0_0_y_z_x_z = buffer_1000_pppp[47];

    auto g_x_0_0_0_y_z_y_x = buffer_1000_pppp[48];

    auto g_x_0_0_0_y_z_y_y = buffer_1000_pppp[49];

    auto g_x_0_0_0_y_z_y_z = buffer_1000_pppp[50];

    auto g_x_0_0_0_y_z_z_x = buffer_1000_pppp[51];

    auto g_x_0_0_0_y_z_z_y = buffer_1000_pppp[52];

    auto g_x_0_0_0_y_z_z_z = buffer_1000_pppp[53];

    auto g_x_0_0_0_z_x_x_x = buffer_1000_pppp[54];

    auto g_x_0_0_0_z_x_x_y = buffer_1000_pppp[55];

    auto g_x_0_0_0_z_x_x_z = buffer_1000_pppp[56];

    auto g_x_0_0_0_z_x_y_x = buffer_1000_pppp[57];

    auto g_x_0_0_0_z_x_y_y = buffer_1000_pppp[58];

    auto g_x_0_0_0_z_x_y_z = buffer_1000_pppp[59];

    auto g_x_0_0_0_z_x_z_x = buffer_1000_pppp[60];

    auto g_x_0_0_0_z_x_z_y = buffer_1000_pppp[61];

    auto g_x_0_0_0_z_x_z_z = buffer_1000_pppp[62];

    auto g_x_0_0_0_z_y_x_x = buffer_1000_pppp[63];

    auto g_x_0_0_0_z_y_x_y = buffer_1000_pppp[64];

    auto g_x_0_0_0_z_y_x_z = buffer_1000_pppp[65];

    auto g_x_0_0_0_z_y_y_x = buffer_1000_pppp[66];

    auto g_x_0_0_0_z_y_y_y = buffer_1000_pppp[67];

    auto g_x_0_0_0_z_y_y_z = buffer_1000_pppp[68];

    auto g_x_0_0_0_z_y_z_x = buffer_1000_pppp[69];

    auto g_x_0_0_0_z_y_z_y = buffer_1000_pppp[70];

    auto g_x_0_0_0_z_y_z_z = buffer_1000_pppp[71];

    auto g_x_0_0_0_z_z_x_x = buffer_1000_pppp[72];

    auto g_x_0_0_0_z_z_x_y = buffer_1000_pppp[73];

    auto g_x_0_0_0_z_z_x_z = buffer_1000_pppp[74];

    auto g_x_0_0_0_z_z_y_x = buffer_1000_pppp[75];

    auto g_x_0_0_0_z_z_y_y = buffer_1000_pppp[76];

    auto g_x_0_0_0_z_z_y_z = buffer_1000_pppp[77];

    auto g_x_0_0_0_z_z_z_x = buffer_1000_pppp[78];

    auto g_x_0_0_0_z_z_z_y = buffer_1000_pppp[79];

    auto g_x_0_0_0_z_z_z_z = buffer_1000_pppp[80];

    auto g_y_0_0_0_x_x_x_x = buffer_1000_pppp[81];

    auto g_y_0_0_0_x_x_x_y = buffer_1000_pppp[82];

    auto g_y_0_0_0_x_x_x_z = buffer_1000_pppp[83];

    auto g_y_0_0_0_x_x_y_x = buffer_1000_pppp[84];

    auto g_y_0_0_0_x_x_y_y = buffer_1000_pppp[85];

    auto g_y_0_0_0_x_x_y_z = buffer_1000_pppp[86];

    auto g_y_0_0_0_x_x_z_x = buffer_1000_pppp[87];

    auto g_y_0_0_0_x_x_z_y = buffer_1000_pppp[88];

    auto g_y_0_0_0_x_x_z_z = buffer_1000_pppp[89];

    auto g_y_0_0_0_x_y_x_x = buffer_1000_pppp[90];

    auto g_y_0_0_0_x_y_x_y = buffer_1000_pppp[91];

    auto g_y_0_0_0_x_y_x_z = buffer_1000_pppp[92];

    auto g_y_0_0_0_x_y_y_x = buffer_1000_pppp[93];

    auto g_y_0_0_0_x_y_y_y = buffer_1000_pppp[94];

    auto g_y_0_0_0_x_y_y_z = buffer_1000_pppp[95];

    auto g_y_0_0_0_x_y_z_x = buffer_1000_pppp[96];

    auto g_y_0_0_0_x_y_z_y = buffer_1000_pppp[97];

    auto g_y_0_0_0_x_y_z_z = buffer_1000_pppp[98];

    auto g_y_0_0_0_x_z_x_x = buffer_1000_pppp[99];

    auto g_y_0_0_0_x_z_x_y = buffer_1000_pppp[100];

    auto g_y_0_0_0_x_z_x_z = buffer_1000_pppp[101];

    auto g_y_0_0_0_x_z_y_x = buffer_1000_pppp[102];

    auto g_y_0_0_0_x_z_y_y = buffer_1000_pppp[103];

    auto g_y_0_0_0_x_z_y_z = buffer_1000_pppp[104];

    auto g_y_0_0_0_x_z_z_x = buffer_1000_pppp[105];

    auto g_y_0_0_0_x_z_z_y = buffer_1000_pppp[106];

    auto g_y_0_0_0_x_z_z_z = buffer_1000_pppp[107];

    auto g_y_0_0_0_y_x_x_x = buffer_1000_pppp[108];

    auto g_y_0_0_0_y_x_x_y = buffer_1000_pppp[109];

    auto g_y_0_0_0_y_x_x_z = buffer_1000_pppp[110];

    auto g_y_0_0_0_y_x_y_x = buffer_1000_pppp[111];

    auto g_y_0_0_0_y_x_y_y = buffer_1000_pppp[112];

    auto g_y_0_0_0_y_x_y_z = buffer_1000_pppp[113];

    auto g_y_0_0_0_y_x_z_x = buffer_1000_pppp[114];

    auto g_y_0_0_0_y_x_z_y = buffer_1000_pppp[115];

    auto g_y_0_0_0_y_x_z_z = buffer_1000_pppp[116];

    auto g_y_0_0_0_y_y_x_x = buffer_1000_pppp[117];

    auto g_y_0_0_0_y_y_x_y = buffer_1000_pppp[118];

    auto g_y_0_0_0_y_y_x_z = buffer_1000_pppp[119];

    auto g_y_0_0_0_y_y_y_x = buffer_1000_pppp[120];

    auto g_y_0_0_0_y_y_y_y = buffer_1000_pppp[121];

    auto g_y_0_0_0_y_y_y_z = buffer_1000_pppp[122];

    auto g_y_0_0_0_y_y_z_x = buffer_1000_pppp[123];

    auto g_y_0_0_0_y_y_z_y = buffer_1000_pppp[124];

    auto g_y_0_0_0_y_y_z_z = buffer_1000_pppp[125];

    auto g_y_0_0_0_y_z_x_x = buffer_1000_pppp[126];

    auto g_y_0_0_0_y_z_x_y = buffer_1000_pppp[127];

    auto g_y_0_0_0_y_z_x_z = buffer_1000_pppp[128];

    auto g_y_0_0_0_y_z_y_x = buffer_1000_pppp[129];

    auto g_y_0_0_0_y_z_y_y = buffer_1000_pppp[130];

    auto g_y_0_0_0_y_z_y_z = buffer_1000_pppp[131];

    auto g_y_0_0_0_y_z_z_x = buffer_1000_pppp[132];

    auto g_y_0_0_0_y_z_z_y = buffer_1000_pppp[133];

    auto g_y_0_0_0_y_z_z_z = buffer_1000_pppp[134];

    auto g_y_0_0_0_z_x_x_x = buffer_1000_pppp[135];

    auto g_y_0_0_0_z_x_x_y = buffer_1000_pppp[136];

    auto g_y_0_0_0_z_x_x_z = buffer_1000_pppp[137];

    auto g_y_0_0_0_z_x_y_x = buffer_1000_pppp[138];

    auto g_y_0_0_0_z_x_y_y = buffer_1000_pppp[139];

    auto g_y_0_0_0_z_x_y_z = buffer_1000_pppp[140];

    auto g_y_0_0_0_z_x_z_x = buffer_1000_pppp[141];

    auto g_y_0_0_0_z_x_z_y = buffer_1000_pppp[142];

    auto g_y_0_0_0_z_x_z_z = buffer_1000_pppp[143];

    auto g_y_0_0_0_z_y_x_x = buffer_1000_pppp[144];

    auto g_y_0_0_0_z_y_x_y = buffer_1000_pppp[145];

    auto g_y_0_0_0_z_y_x_z = buffer_1000_pppp[146];

    auto g_y_0_0_0_z_y_y_x = buffer_1000_pppp[147];

    auto g_y_0_0_0_z_y_y_y = buffer_1000_pppp[148];

    auto g_y_0_0_0_z_y_y_z = buffer_1000_pppp[149];

    auto g_y_0_0_0_z_y_z_x = buffer_1000_pppp[150];

    auto g_y_0_0_0_z_y_z_y = buffer_1000_pppp[151];

    auto g_y_0_0_0_z_y_z_z = buffer_1000_pppp[152];

    auto g_y_0_0_0_z_z_x_x = buffer_1000_pppp[153];

    auto g_y_0_0_0_z_z_x_y = buffer_1000_pppp[154];

    auto g_y_0_0_0_z_z_x_z = buffer_1000_pppp[155];

    auto g_y_0_0_0_z_z_y_x = buffer_1000_pppp[156];

    auto g_y_0_0_0_z_z_y_y = buffer_1000_pppp[157];

    auto g_y_0_0_0_z_z_y_z = buffer_1000_pppp[158];

    auto g_y_0_0_0_z_z_z_x = buffer_1000_pppp[159];

    auto g_y_0_0_0_z_z_z_y = buffer_1000_pppp[160];

    auto g_y_0_0_0_z_z_z_z = buffer_1000_pppp[161];

    auto g_z_0_0_0_x_x_x_x = buffer_1000_pppp[162];

    auto g_z_0_0_0_x_x_x_y = buffer_1000_pppp[163];

    auto g_z_0_0_0_x_x_x_z = buffer_1000_pppp[164];

    auto g_z_0_0_0_x_x_y_x = buffer_1000_pppp[165];

    auto g_z_0_0_0_x_x_y_y = buffer_1000_pppp[166];

    auto g_z_0_0_0_x_x_y_z = buffer_1000_pppp[167];

    auto g_z_0_0_0_x_x_z_x = buffer_1000_pppp[168];

    auto g_z_0_0_0_x_x_z_y = buffer_1000_pppp[169];

    auto g_z_0_0_0_x_x_z_z = buffer_1000_pppp[170];

    auto g_z_0_0_0_x_y_x_x = buffer_1000_pppp[171];

    auto g_z_0_0_0_x_y_x_y = buffer_1000_pppp[172];

    auto g_z_0_0_0_x_y_x_z = buffer_1000_pppp[173];

    auto g_z_0_0_0_x_y_y_x = buffer_1000_pppp[174];

    auto g_z_0_0_0_x_y_y_y = buffer_1000_pppp[175];

    auto g_z_0_0_0_x_y_y_z = buffer_1000_pppp[176];

    auto g_z_0_0_0_x_y_z_x = buffer_1000_pppp[177];

    auto g_z_0_0_0_x_y_z_y = buffer_1000_pppp[178];

    auto g_z_0_0_0_x_y_z_z = buffer_1000_pppp[179];

    auto g_z_0_0_0_x_z_x_x = buffer_1000_pppp[180];

    auto g_z_0_0_0_x_z_x_y = buffer_1000_pppp[181];

    auto g_z_0_0_0_x_z_x_z = buffer_1000_pppp[182];

    auto g_z_0_0_0_x_z_y_x = buffer_1000_pppp[183];

    auto g_z_0_0_0_x_z_y_y = buffer_1000_pppp[184];

    auto g_z_0_0_0_x_z_y_z = buffer_1000_pppp[185];

    auto g_z_0_0_0_x_z_z_x = buffer_1000_pppp[186];

    auto g_z_0_0_0_x_z_z_y = buffer_1000_pppp[187];

    auto g_z_0_0_0_x_z_z_z = buffer_1000_pppp[188];

    auto g_z_0_0_0_y_x_x_x = buffer_1000_pppp[189];

    auto g_z_0_0_0_y_x_x_y = buffer_1000_pppp[190];

    auto g_z_0_0_0_y_x_x_z = buffer_1000_pppp[191];

    auto g_z_0_0_0_y_x_y_x = buffer_1000_pppp[192];

    auto g_z_0_0_0_y_x_y_y = buffer_1000_pppp[193];

    auto g_z_0_0_0_y_x_y_z = buffer_1000_pppp[194];

    auto g_z_0_0_0_y_x_z_x = buffer_1000_pppp[195];

    auto g_z_0_0_0_y_x_z_y = buffer_1000_pppp[196];

    auto g_z_0_0_0_y_x_z_z = buffer_1000_pppp[197];

    auto g_z_0_0_0_y_y_x_x = buffer_1000_pppp[198];

    auto g_z_0_0_0_y_y_x_y = buffer_1000_pppp[199];

    auto g_z_0_0_0_y_y_x_z = buffer_1000_pppp[200];

    auto g_z_0_0_0_y_y_y_x = buffer_1000_pppp[201];

    auto g_z_0_0_0_y_y_y_y = buffer_1000_pppp[202];

    auto g_z_0_0_0_y_y_y_z = buffer_1000_pppp[203];

    auto g_z_0_0_0_y_y_z_x = buffer_1000_pppp[204];

    auto g_z_0_0_0_y_y_z_y = buffer_1000_pppp[205];

    auto g_z_0_0_0_y_y_z_z = buffer_1000_pppp[206];

    auto g_z_0_0_0_y_z_x_x = buffer_1000_pppp[207];

    auto g_z_0_0_0_y_z_x_y = buffer_1000_pppp[208];

    auto g_z_0_0_0_y_z_x_z = buffer_1000_pppp[209];

    auto g_z_0_0_0_y_z_y_x = buffer_1000_pppp[210];

    auto g_z_0_0_0_y_z_y_y = buffer_1000_pppp[211];

    auto g_z_0_0_0_y_z_y_z = buffer_1000_pppp[212];

    auto g_z_0_0_0_y_z_z_x = buffer_1000_pppp[213];

    auto g_z_0_0_0_y_z_z_y = buffer_1000_pppp[214];

    auto g_z_0_0_0_y_z_z_z = buffer_1000_pppp[215];

    auto g_z_0_0_0_z_x_x_x = buffer_1000_pppp[216];

    auto g_z_0_0_0_z_x_x_y = buffer_1000_pppp[217];

    auto g_z_0_0_0_z_x_x_z = buffer_1000_pppp[218];

    auto g_z_0_0_0_z_x_y_x = buffer_1000_pppp[219];

    auto g_z_0_0_0_z_x_y_y = buffer_1000_pppp[220];

    auto g_z_0_0_0_z_x_y_z = buffer_1000_pppp[221];

    auto g_z_0_0_0_z_x_z_x = buffer_1000_pppp[222];

    auto g_z_0_0_0_z_x_z_y = buffer_1000_pppp[223];

    auto g_z_0_0_0_z_x_z_z = buffer_1000_pppp[224];

    auto g_z_0_0_0_z_y_x_x = buffer_1000_pppp[225];

    auto g_z_0_0_0_z_y_x_y = buffer_1000_pppp[226];

    auto g_z_0_0_0_z_y_x_z = buffer_1000_pppp[227];

    auto g_z_0_0_0_z_y_y_x = buffer_1000_pppp[228];

    auto g_z_0_0_0_z_y_y_y = buffer_1000_pppp[229];

    auto g_z_0_0_0_z_y_y_z = buffer_1000_pppp[230];

    auto g_z_0_0_0_z_y_z_x = buffer_1000_pppp[231];

    auto g_z_0_0_0_z_y_z_y = buffer_1000_pppp[232];

    auto g_z_0_0_0_z_y_z_z = buffer_1000_pppp[233];

    auto g_z_0_0_0_z_z_x_x = buffer_1000_pppp[234];

    auto g_z_0_0_0_z_z_x_y = buffer_1000_pppp[235];

    auto g_z_0_0_0_z_z_x_z = buffer_1000_pppp[236];

    auto g_z_0_0_0_z_z_y_x = buffer_1000_pppp[237];

    auto g_z_0_0_0_z_z_y_y = buffer_1000_pppp[238];

    auto g_z_0_0_0_z_z_y_z = buffer_1000_pppp[239];

    auto g_z_0_0_0_z_z_z_x = buffer_1000_pppp[240];

    auto g_z_0_0_0_z_z_z_y = buffer_1000_pppp[241];

    auto g_z_0_0_0_z_z_z_z = buffer_1000_pppp[242];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_y, g_0_x_x_z, g_x_0_0_0_x_x_x_x, g_x_0_0_0_x_x_x_y, g_x_0_0_0_x_x_x_z, g_xx_x_x_x, g_xx_x_x_y, g_xx_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_x_x_x[i] = -g_0_x_x_x[i] + 2.0 * g_xx_x_x_x[i] * a_exp;

        g_x_0_0_0_x_x_x_y[i] = -g_0_x_x_y[i] + 2.0 * g_xx_x_x_y[i] * a_exp;

        g_x_0_0_0_x_x_x_z[i] = -g_0_x_x_z[i] + 2.0 * g_xx_x_x_z[i] * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_0_x_y_x, g_0_x_y_y, g_0_x_y_z, g_x_0_0_0_x_x_y_x, g_x_0_0_0_x_x_y_y, g_x_0_0_0_x_x_y_z, g_xx_x_y_x, g_xx_x_y_y, g_xx_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_x_y_x[i] = -g_0_x_y_x[i] + 2.0 * g_xx_x_y_x[i] * a_exp;

        g_x_0_0_0_x_x_y_y[i] = -g_0_x_y_y[i] + 2.0 * g_xx_x_y_y[i] * a_exp;

        g_x_0_0_0_x_x_y_z[i] = -g_0_x_y_z[i] + 2.0 * g_xx_x_y_z[i] * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_0_x_z_x, g_0_x_z_y, g_0_x_z_z, g_x_0_0_0_x_x_z_x, g_x_0_0_0_x_x_z_y, g_x_0_0_0_x_x_z_z, g_xx_x_z_x, g_xx_x_z_y, g_xx_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_x_z_x[i] = -g_0_x_z_x[i] + 2.0 * g_xx_x_z_x[i] * a_exp;

        g_x_0_0_0_x_x_z_y[i] = -g_0_x_z_y[i] + 2.0 * g_xx_x_z_y[i] * a_exp;

        g_x_0_0_0_x_x_z_z[i] = -g_0_x_z_z[i] + 2.0 * g_xx_x_z_z[i] * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_0_y_x_x, g_0_y_x_y, g_0_y_x_z, g_x_0_0_0_x_y_x_x, g_x_0_0_0_x_y_x_y, g_x_0_0_0_x_y_x_z, g_xx_y_x_x, g_xx_y_x_y, g_xx_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_y_x_x[i] = -g_0_y_x_x[i] + 2.0 * g_xx_y_x_x[i] * a_exp;

        g_x_0_0_0_x_y_x_y[i] = -g_0_y_x_y[i] + 2.0 * g_xx_y_x_y[i] * a_exp;

        g_x_0_0_0_x_y_x_z[i] = -g_0_y_x_z[i] + 2.0 * g_xx_y_x_z[i] * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_0_y_y_x, g_0_y_y_y, g_0_y_y_z, g_x_0_0_0_x_y_y_x, g_x_0_0_0_x_y_y_y, g_x_0_0_0_x_y_y_z, g_xx_y_y_x, g_xx_y_y_y, g_xx_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_y_y_x[i] = -g_0_y_y_x[i] + 2.0 * g_xx_y_y_x[i] * a_exp;

        g_x_0_0_0_x_y_y_y[i] = -g_0_y_y_y[i] + 2.0 * g_xx_y_y_y[i] * a_exp;

        g_x_0_0_0_x_y_y_z[i] = -g_0_y_y_z[i] + 2.0 * g_xx_y_y_z[i] * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_0_y_z_x, g_0_y_z_y, g_0_y_z_z, g_x_0_0_0_x_y_z_x, g_x_0_0_0_x_y_z_y, g_x_0_0_0_x_y_z_z, g_xx_y_z_x, g_xx_y_z_y, g_xx_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_y_z_x[i] = -g_0_y_z_x[i] + 2.0 * g_xx_y_z_x[i] * a_exp;

        g_x_0_0_0_x_y_z_y[i] = -g_0_y_z_y[i] + 2.0 * g_xx_y_z_y[i] * a_exp;

        g_x_0_0_0_x_y_z_z[i] = -g_0_y_z_z[i] + 2.0 * g_xx_y_z_z[i] * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_0_z_x_x, g_0_z_x_y, g_0_z_x_z, g_x_0_0_0_x_z_x_x, g_x_0_0_0_x_z_x_y, g_x_0_0_0_x_z_x_z, g_xx_z_x_x, g_xx_z_x_y, g_xx_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_z_x_x[i] = -g_0_z_x_x[i] + 2.0 * g_xx_z_x_x[i] * a_exp;

        g_x_0_0_0_x_z_x_y[i] = -g_0_z_x_y[i] + 2.0 * g_xx_z_x_y[i] * a_exp;

        g_x_0_0_0_x_z_x_z[i] = -g_0_z_x_z[i] + 2.0 * g_xx_z_x_z[i] * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_0_z_y_x, g_0_z_y_y, g_0_z_y_z, g_x_0_0_0_x_z_y_x, g_x_0_0_0_x_z_y_y, g_x_0_0_0_x_z_y_z, g_xx_z_y_x, g_xx_z_y_y, g_xx_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_z_y_x[i] = -g_0_z_y_x[i] + 2.0 * g_xx_z_y_x[i] * a_exp;

        g_x_0_0_0_x_z_y_y[i] = -g_0_z_y_y[i] + 2.0 * g_xx_z_y_y[i] * a_exp;

        g_x_0_0_0_x_z_y_z[i] = -g_0_z_y_z[i] + 2.0 * g_xx_z_y_z[i] * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_x_0_0_0_x_z_z_x, g_x_0_0_0_x_z_z_y, g_x_0_0_0_x_z_z_z, g_xx_z_z_x, g_xx_z_z_y, g_xx_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_z_z_x[i] = -g_0_z_z_x[i] + 2.0 * g_xx_z_z_x[i] * a_exp;

        g_x_0_0_0_x_z_z_y[i] = -g_0_z_z_y[i] + 2.0 * g_xx_z_z_y[i] * a_exp;

        g_x_0_0_0_x_z_z_z[i] = -g_0_z_z_z[i] + 2.0 * g_xx_z_z_z[i] * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_0_0_0_y_x_x_x, g_x_0_0_0_y_x_x_y, g_x_0_0_0_y_x_x_z, g_xy_x_x_x, g_xy_x_x_y, g_xy_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_x_x_x[i] = 2.0 * g_xy_x_x_x[i] * a_exp;

        g_x_0_0_0_y_x_x_y[i] = 2.0 * g_xy_x_x_y[i] * a_exp;

        g_x_0_0_0_y_x_x_z[i] = 2.0 * g_xy_x_x_z[i] * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_x_0_0_0_y_x_y_x, g_x_0_0_0_y_x_y_y, g_x_0_0_0_y_x_y_z, g_xy_x_y_x, g_xy_x_y_y, g_xy_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_x_y_x[i] = 2.0 * g_xy_x_y_x[i] * a_exp;

        g_x_0_0_0_y_x_y_y[i] = 2.0 * g_xy_x_y_y[i] * a_exp;

        g_x_0_0_0_y_x_y_z[i] = 2.0 * g_xy_x_y_z[i] * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_x_0_0_0_y_x_z_x, g_x_0_0_0_y_x_z_y, g_x_0_0_0_y_x_z_z, g_xy_x_z_x, g_xy_x_z_y, g_xy_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_x_z_x[i] = 2.0 * g_xy_x_z_x[i] * a_exp;

        g_x_0_0_0_y_x_z_y[i] = 2.0 * g_xy_x_z_y[i] * a_exp;

        g_x_0_0_0_y_x_z_z[i] = 2.0 * g_xy_x_z_z[i] * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_0_0_0_y_y_x_x, g_x_0_0_0_y_y_x_y, g_x_0_0_0_y_y_x_z, g_xy_y_x_x, g_xy_y_x_y, g_xy_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_y_x_x[i] = 2.0 * g_xy_y_x_x[i] * a_exp;

        g_x_0_0_0_y_y_x_y[i] = 2.0 * g_xy_y_x_y[i] * a_exp;

        g_x_0_0_0_y_y_x_z[i] = 2.0 * g_xy_y_x_z[i] * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_0_0_0_y_y_y_x, g_x_0_0_0_y_y_y_y, g_x_0_0_0_y_y_y_z, g_xy_y_y_x, g_xy_y_y_y, g_xy_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_y_y_x[i] = 2.0 * g_xy_y_y_x[i] * a_exp;

        g_x_0_0_0_y_y_y_y[i] = 2.0 * g_xy_y_y_y[i] * a_exp;

        g_x_0_0_0_y_y_y_z[i] = 2.0 * g_xy_y_y_z[i] * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_0_0_0_y_y_z_x, g_x_0_0_0_y_y_z_y, g_x_0_0_0_y_y_z_z, g_xy_y_z_x, g_xy_y_z_y, g_xy_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_y_z_x[i] = 2.0 * g_xy_y_z_x[i] * a_exp;

        g_x_0_0_0_y_y_z_y[i] = 2.0 * g_xy_y_z_y[i] * a_exp;

        g_x_0_0_0_y_y_z_z[i] = 2.0 * g_xy_y_z_z[i] * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_0_0_0_y_z_x_x, g_x_0_0_0_y_z_x_y, g_x_0_0_0_y_z_x_z, g_xy_z_x_x, g_xy_z_x_y, g_xy_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_z_x_x[i] = 2.0 * g_xy_z_x_x[i] * a_exp;

        g_x_0_0_0_y_z_x_y[i] = 2.0 * g_xy_z_x_y[i] * a_exp;

        g_x_0_0_0_y_z_x_z[i] = 2.0 * g_xy_z_x_z[i] * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_0_0_0_y_z_y_x, g_x_0_0_0_y_z_y_y, g_x_0_0_0_y_z_y_z, g_xy_z_y_x, g_xy_z_y_y, g_xy_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_z_y_x[i] = 2.0 * g_xy_z_y_x[i] * a_exp;

        g_x_0_0_0_y_z_y_y[i] = 2.0 * g_xy_z_y_y[i] * a_exp;

        g_x_0_0_0_y_z_y_z[i] = 2.0 * g_xy_z_y_z[i] * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_0_0_0_y_z_z_x, g_x_0_0_0_y_z_z_y, g_x_0_0_0_y_z_z_z, g_xy_z_z_x, g_xy_z_z_y, g_xy_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_z_z_x[i] = 2.0 * g_xy_z_z_x[i] * a_exp;

        g_x_0_0_0_y_z_z_y[i] = 2.0 * g_xy_z_z_y[i] * a_exp;

        g_x_0_0_0_y_z_z_z[i] = 2.0 * g_xy_z_z_z[i] * a_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_x_0_0_0_z_x_x_x, g_x_0_0_0_z_x_x_y, g_x_0_0_0_z_x_x_z, g_xz_x_x_x, g_xz_x_x_y, g_xz_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_x_x_x[i] = 2.0 * g_xz_x_x_x[i] * a_exp;

        g_x_0_0_0_z_x_x_y[i] = 2.0 * g_xz_x_x_y[i] * a_exp;

        g_x_0_0_0_z_x_x_z[i] = 2.0 * g_xz_x_x_z[i] * a_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_x_0_0_0_z_x_y_x, g_x_0_0_0_z_x_y_y, g_x_0_0_0_z_x_y_z, g_xz_x_y_x, g_xz_x_y_y, g_xz_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_x_y_x[i] = 2.0 * g_xz_x_y_x[i] * a_exp;

        g_x_0_0_0_z_x_y_y[i] = 2.0 * g_xz_x_y_y[i] * a_exp;

        g_x_0_0_0_z_x_y_z[i] = 2.0 * g_xz_x_y_z[i] * a_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_x_0_0_0_z_x_z_x, g_x_0_0_0_z_x_z_y, g_x_0_0_0_z_x_z_z, g_xz_x_z_x, g_xz_x_z_y, g_xz_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_x_z_x[i] = 2.0 * g_xz_x_z_x[i] * a_exp;

        g_x_0_0_0_z_x_z_y[i] = 2.0 * g_xz_x_z_y[i] * a_exp;

        g_x_0_0_0_z_x_z_z[i] = 2.0 * g_xz_x_z_z[i] * a_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_x_0_0_0_z_y_x_x, g_x_0_0_0_z_y_x_y, g_x_0_0_0_z_y_x_z, g_xz_y_x_x, g_xz_y_x_y, g_xz_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_y_x_x[i] = 2.0 * g_xz_y_x_x[i] * a_exp;

        g_x_0_0_0_z_y_x_y[i] = 2.0 * g_xz_y_x_y[i] * a_exp;

        g_x_0_0_0_z_y_x_z[i] = 2.0 * g_xz_y_x_z[i] * a_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_x_0_0_0_z_y_y_x, g_x_0_0_0_z_y_y_y, g_x_0_0_0_z_y_y_z, g_xz_y_y_x, g_xz_y_y_y, g_xz_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_y_y_x[i] = 2.0 * g_xz_y_y_x[i] * a_exp;

        g_x_0_0_0_z_y_y_y[i] = 2.0 * g_xz_y_y_y[i] * a_exp;

        g_x_0_0_0_z_y_y_z[i] = 2.0 * g_xz_y_y_z[i] * a_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_x_0_0_0_z_y_z_x, g_x_0_0_0_z_y_z_y, g_x_0_0_0_z_y_z_z, g_xz_y_z_x, g_xz_y_z_y, g_xz_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_y_z_x[i] = 2.0 * g_xz_y_z_x[i] * a_exp;

        g_x_0_0_0_z_y_z_y[i] = 2.0 * g_xz_y_z_y[i] * a_exp;

        g_x_0_0_0_z_y_z_z[i] = 2.0 * g_xz_y_z_z[i] * a_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_0_0_0_z_z_x_x, g_x_0_0_0_z_z_x_y, g_x_0_0_0_z_z_x_z, g_xz_z_x_x, g_xz_z_x_y, g_xz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_z_x_x[i] = 2.0 * g_xz_z_x_x[i] * a_exp;

        g_x_0_0_0_z_z_x_y[i] = 2.0 * g_xz_z_x_y[i] * a_exp;

        g_x_0_0_0_z_z_x_z[i] = 2.0 * g_xz_z_x_z[i] * a_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_0_0_0_z_z_y_x, g_x_0_0_0_z_z_y_y, g_x_0_0_0_z_z_y_z, g_xz_z_y_x, g_xz_z_y_y, g_xz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_z_y_x[i] = 2.0 * g_xz_z_y_x[i] * a_exp;

        g_x_0_0_0_z_z_y_y[i] = 2.0 * g_xz_z_y_y[i] * a_exp;

        g_x_0_0_0_z_z_y_z[i] = 2.0 * g_xz_z_y_z[i] * a_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_0_0_0_z_z_z_x, g_x_0_0_0_z_z_z_y, g_x_0_0_0_z_z_z_z, g_xz_z_z_x, g_xz_z_z_y, g_xz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_z_z_x[i] = 2.0 * g_xz_z_z_x[i] * a_exp;

        g_x_0_0_0_z_z_z_y[i] = 2.0 * g_xz_z_z_y[i] * a_exp;

        g_x_0_0_0_z_z_z_z[i] = 2.0 * g_xz_z_z_z[i] * a_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_xy_x_x_x, g_xy_x_x_y, g_xy_x_x_z, g_y_0_0_0_x_x_x_x, g_y_0_0_0_x_x_x_y, g_y_0_0_0_x_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_x_x_x[i] = 2.0 * g_xy_x_x_x[i] * a_exp;

        g_y_0_0_0_x_x_x_y[i] = 2.0 * g_xy_x_x_y[i] * a_exp;

        g_y_0_0_0_x_x_x_z[i] = 2.0 * g_xy_x_x_z[i] * a_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_xy_x_y_x, g_xy_x_y_y, g_xy_x_y_z, g_y_0_0_0_x_x_y_x, g_y_0_0_0_x_x_y_y, g_y_0_0_0_x_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_x_y_x[i] = 2.0 * g_xy_x_y_x[i] * a_exp;

        g_y_0_0_0_x_x_y_y[i] = 2.0 * g_xy_x_y_y[i] * a_exp;

        g_y_0_0_0_x_x_y_z[i] = 2.0 * g_xy_x_y_z[i] * a_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_xy_x_z_x, g_xy_x_z_y, g_xy_x_z_z, g_y_0_0_0_x_x_z_x, g_y_0_0_0_x_x_z_y, g_y_0_0_0_x_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_x_z_x[i] = 2.0 * g_xy_x_z_x[i] * a_exp;

        g_y_0_0_0_x_x_z_y[i] = 2.0 * g_xy_x_z_y[i] * a_exp;

        g_y_0_0_0_x_x_z_z[i] = 2.0 * g_xy_x_z_z[i] * a_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_xy_y_x_x, g_xy_y_x_y, g_xy_y_x_z, g_y_0_0_0_x_y_x_x, g_y_0_0_0_x_y_x_y, g_y_0_0_0_x_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_y_x_x[i] = 2.0 * g_xy_y_x_x[i] * a_exp;

        g_y_0_0_0_x_y_x_y[i] = 2.0 * g_xy_y_x_y[i] * a_exp;

        g_y_0_0_0_x_y_x_z[i] = 2.0 * g_xy_y_x_z[i] * a_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_xy_y_y_x, g_xy_y_y_y, g_xy_y_y_z, g_y_0_0_0_x_y_y_x, g_y_0_0_0_x_y_y_y, g_y_0_0_0_x_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_y_y_x[i] = 2.0 * g_xy_y_y_x[i] * a_exp;

        g_y_0_0_0_x_y_y_y[i] = 2.0 * g_xy_y_y_y[i] * a_exp;

        g_y_0_0_0_x_y_y_z[i] = 2.0 * g_xy_y_y_z[i] * a_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_xy_y_z_x, g_xy_y_z_y, g_xy_y_z_z, g_y_0_0_0_x_y_z_x, g_y_0_0_0_x_y_z_y, g_y_0_0_0_x_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_y_z_x[i] = 2.0 * g_xy_y_z_x[i] * a_exp;

        g_y_0_0_0_x_y_z_y[i] = 2.0 * g_xy_y_z_y[i] * a_exp;

        g_y_0_0_0_x_y_z_z[i] = 2.0 * g_xy_y_z_z[i] * a_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_xy_z_x_x, g_xy_z_x_y, g_xy_z_x_z, g_y_0_0_0_x_z_x_x, g_y_0_0_0_x_z_x_y, g_y_0_0_0_x_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_z_x_x[i] = 2.0 * g_xy_z_x_x[i] * a_exp;

        g_y_0_0_0_x_z_x_y[i] = 2.0 * g_xy_z_x_y[i] * a_exp;

        g_y_0_0_0_x_z_x_z[i] = 2.0 * g_xy_z_x_z[i] * a_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_xy_z_y_x, g_xy_z_y_y, g_xy_z_y_z, g_y_0_0_0_x_z_y_x, g_y_0_0_0_x_z_y_y, g_y_0_0_0_x_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_z_y_x[i] = 2.0 * g_xy_z_y_x[i] * a_exp;

        g_y_0_0_0_x_z_y_y[i] = 2.0 * g_xy_z_y_y[i] * a_exp;

        g_y_0_0_0_x_z_y_z[i] = 2.0 * g_xy_z_y_z[i] * a_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_xy_z_z_x, g_xy_z_z_y, g_xy_z_z_z, g_y_0_0_0_x_z_z_x, g_y_0_0_0_x_z_z_y, g_y_0_0_0_x_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_z_z_x[i] = 2.0 * g_xy_z_z_x[i] * a_exp;

        g_y_0_0_0_x_z_z_y[i] = 2.0 * g_xy_z_z_y[i] * a_exp;

        g_y_0_0_0_x_z_z_z[i] = 2.0 * g_xy_z_z_z[i] * a_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_y, g_0_x_x_z, g_y_0_0_0_y_x_x_x, g_y_0_0_0_y_x_x_y, g_y_0_0_0_y_x_x_z, g_yy_x_x_x, g_yy_x_x_y, g_yy_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_x_x_x[i] = -g_0_x_x_x[i] + 2.0 * g_yy_x_x_x[i] * a_exp;

        g_y_0_0_0_y_x_x_y[i] = -g_0_x_x_y[i] + 2.0 * g_yy_x_x_y[i] * a_exp;

        g_y_0_0_0_y_x_x_z[i] = -g_0_x_x_z[i] + 2.0 * g_yy_x_x_z[i] * a_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_0_x_y_x, g_0_x_y_y, g_0_x_y_z, g_y_0_0_0_y_x_y_x, g_y_0_0_0_y_x_y_y, g_y_0_0_0_y_x_y_z, g_yy_x_y_x, g_yy_x_y_y, g_yy_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_x_y_x[i] = -g_0_x_y_x[i] + 2.0 * g_yy_x_y_x[i] * a_exp;

        g_y_0_0_0_y_x_y_y[i] = -g_0_x_y_y[i] + 2.0 * g_yy_x_y_y[i] * a_exp;

        g_y_0_0_0_y_x_y_z[i] = -g_0_x_y_z[i] + 2.0 * g_yy_x_y_z[i] * a_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_0_x_z_x, g_0_x_z_y, g_0_x_z_z, g_y_0_0_0_y_x_z_x, g_y_0_0_0_y_x_z_y, g_y_0_0_0_y_x_z_z, g_yy_x_z_x, g_yy_x_z_y, g_yy_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_x_z_x[i] = -g_0_x_z_x[i] + 2.0 * g_yy_x_z_x[i] * a_exp;

        g_y_0_0_0_y_x_z_y[i] = -g_0_x_z_y[i] + 2.0 * g_yy_x_z_y[i] * a_exp;

        g_y_0_0_0_y_x_z_z[i] = -g_0_x_z_z[i] + 2.0 * g_yy_x_z_z[i] * a_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_0_y_x_x, g_0_y_x_y, g_0_y_x_z, g_y_0_0_0_y_y_x_x, g_y_0_0_0_y_y_x_y, g_y_0_0_0_y_y_x_z, g_yy_y_x_x, g_yy_y_x_y, g_yy_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_y_x_x[i] = -g_0_y_x_x[i] + 2.0 * g_yy_y_x_x[i] * a_exp;

        g_y_0_0_0_y_y_x_y[i] = -g_0_y_x_y[i] + 2.0 * g_yy_y_x_y[i] * a_exp;

        g_y_0_0_0_y_y_x_z[i] = -g_0_y_x_z[i] + 2.0 * g_yy_y_x_z[i] * a_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_0_y_y_x, g_0_y_y_y, g_0_y_y_z, g_y_0_0_0_y_y_y_x, g_y_0_0_0_y_y_y_y, g_y_0_0_0_y_y_y_z, g_yy_y_y_x, g_yy_y_y_y, g_yy_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_y_y_x[i] = -g_0_y_y_x[i] + 2.0 * g_yy_y_y_x[i] * a_exp;

        g_y_0_0_0_y_y_y_y[i] = -g_0_y_y_y[i] + 2.0 * g_yy_y_y_y[i] * a_exp;

        g_y_0_0_0_y_y_y_z[i] = -g_0_y_y_z[i] + 2.0 * g_yy_y_y_z[i] * a_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_0_y_z_x, g_0_y_z_y, g_0_y_z_z, g_y_0_0_0_y_y_z_x, g_y_0_0_0_y_y_z_y, g_y_0_0_0_y_y_z_z, g_yy_y_z_x, g_yy_y_z_y, g_yy_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_y_z_x[i] = -g_0_y_z_x[i] + 2.0 * g_yy_y_z_x[i] * a_exp;

        g_y_0_0_0_y_y_z_y[i] = -g_0_y_z_y[i] + 2.0 * g_yy_y_z_y[i] * a_exp;

        g_y_0_0_0_y_y_z_z[i] = -g_0_y_z_z[i] + 2.0 * g_yy_y_z_z[i] * a_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_0_z_x_x, g_0_z_x_y, g_0_z_x_z, g_y_0_0_0_y_z_x_x, g_y_0_0_0_y_z_x_y, g_y_0_0_0_y_z_x_z, g_yy_z_x_x, g_yy_z_x_y, g_yy_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_z_x_x[i] = -g_0_z_x_x[i] + 2.0 * g_yy_z_x_x[i] * a_exp;

        g_y_0_0_0_y_z_x_y[i] = -g_0_z_x_y[i] + 2.0 * g_yy_z_x_y[i] * a_exp;

        g_y_0_0_0_y_z_x_z[i] = -g_0_z_x_z[i] + 2.0 * g_yy_z_x_z[i] * a_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_0_z_y_x, g_0_z_y_y, g_0_z_y_z, g_y_0_0_0_y_z_y_x, g_y_0_0_0_y_z_y_y, g_y_0_0_0_y_z_y_z, g_yy_z_y_x, g_yy_z_y_y, g_yy_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_z_y_x[i] = -g_0_z_y_x[i] + 2.0 * g_yy_z_y_x[i] * a_exp;

        g_y_0_0_0_y_z_y_y[i] = -g_0_z_y_y[i] + 2.0 * g_yy_z_y_y[i] * a_exp;

        g_y_0_0_0_y_z_y_z[i] = -g_0_z_y_z[i] + 2.0 * g_yy_z_y_z[i] * a_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_y_0_0_0_y_z_z_x, g_y_0_0_0_y_z_z_y, g_y_0_0_0_y_z_z_z, g_yy_z_z_x, g_yy_z_z_y, g_yy_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_z_z_x[i] = -g_0_z_z_x[i] + 2.0 * g_yy_z_z_x[i] * a_exp;

        g_y_0_0_0_y_z_z_y[i] = -g_0_z_z_y[i] + 2.0 * g_yy_z_z_y[i] * a_exp;

        g_y_0_0_0_y_z_z_z[i] = -g_0_z_z_z[i] + 2.0 * g_yy_z_z_z[i] * a_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_y_0_0_0_z_x_x_x, g_y_0_0_0_z_x_x_y, g_y_0_0_0_z_x_x_z, g_yz_x_x_x, g_yz_x_x_y, g_yz_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_x_x_x[i] = 2.0 * g_yz_x_x_x[i] * a_exp;

        g_y_0_0_0_z_x_x_y[i] = 2.0 * g_yz_x_x_y[i] * a_exp;

        g_y_0_0_0_z_x_x_z[i] = 2.0 * g_yz_x_x_z[i] * a_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_y_0_0_0_z_x_y_x, g_y_0_0_0_z_x_y_y, g_y_0_0_0_z_x_y_z, g_yz_x_y_x, g_yz_x_y_y, g_yz_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_x_y_x[i] = 2.0 * g_yz_x_y_x[i] * a_exp;

        g_y_0_0_0_z_x_y_y[i] = 2.0 * g_yz_x_y_y[i] * a_exp;

        g_y_0_0_0_z_x_y_z[i] = 2.0 * g_yz_x_y_z[i] * a_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_y_0_0_0_z_x_z_x, g_y_0_0_0_z_x_z_y, g_y_0_0_0_z_x_z_z, g_yz_x_z_x, g_yz_x_z_y, g_yz_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_x_z_x[i] = 2.0 * g_yz_x_z_x[i] * a_exp;

        g_y_0_0_0_z_x_z_y[i] = 2.0 * g_yz_x_z_y[i] * a_exp;

        g_y_0_0_0_z_x_z_z[i] = 2.0 * g_yz_x_z_z[i] * a_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_y_0_0_0_z_y_x_x, g_y_0_0_0_z_y_x_y, g_y_0_0_0_z_y_x_z, g_yz_y_x_x, g_yz_y_x_y, g_yz_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_y_x_x[i] = 2.0 * g_yz_y_x_x[i] * a_exp;

        g_y_0_0_0_z_y_x_y[i] = 2.0 * g_yz_y_x_y[i] * a_exp;

        g_y_0_0_0_z_y_x_z[i] = 2.0 * g_yz_y_x_z[i] * a_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_y_0_0_0_z_y_y_x, g_y_0_0_0_z_y_y_y, g_y_0_0_0_z_y_y_z, g_yz_y_y_x, g_yz_y_y_y, g_yz_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_y_y_x[i] = 2.0 * g_yz_y_y_x[i] * a_exp;

        g_y_0_0_0_z_y_y_y[i] = 2.0 * g_yz_y_y_y[i] * a_exp;

        g_y_0_0_0_z_y_y_z[i] = 2.0 * g_yz_y_y_z[i] * a_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_y_0_0_0_z_y_z_x, g_y_0_0_0_z_y_z_y, g_y_0_0_0_z_y_z_z, g_yz_y_z_x, g_yz_y_z_y, g_yz_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_y_z_x[i] = 2.0 * g_yz_y_z_x[i] * a_exp;

        g_y_0_0_0_z_y_z_y[i] = 2.0 * g_yz_y_z_y[i] * a_exp;

        g_y_0_0_0_z_y_z_z[i] = 2.0 * g_yz_y_z_z[i] * a_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_y_0_0_0_z_z_x_x, g_y_0_0_0_z_z_x_y, g_y_0_0_0_z_z_x_z, g_yz_z_x_x, g_yz_z_x_y, g_yz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_z_x_x[i] = 2.0 * g_yz_z_x_x[i] * a_exp;

        g_y_0_0_0_z_z_x_y[i] = 2.0 * g_yz_z_x_y[i] * a_exp;

        g_y_0_0_0_z_z_x_z[i] = 2.0 * g_yz_z_x_z[i] * a_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_y_0_0_0_z_z_y_x, g_y_0_0_0_z_z_y_y, g_y_0_0_0_z_z_y_z, g_yz_z_y_x, g_yz_z_y_y, g_yz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_z_y_x[i] = 2.0 * g_yz_z_y_x[i] * a_exp;

        g_y_0_0_0_z_z_y_y[i] = 2.0 * g_yz_z_y_y[i] * a_exp;

        g_y_0_0_0_z_z_y_z[i] = 2.0 * g_yz_z_y_z[i] * a_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_y_0_0_0_z_z_z_x, g_y_0_0_0_z_z_z_y, g_y_0_0_0_z_z_z_z, g_yz_z_z_x, g_yz_z_z_y, g_yz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_z_z_x[i] = 2.0 * g_yz_z_z_x[i] * a_exp;

        g_y_0_0_0_z_z_z_y[i] = 2.0 * g_yz_z_z_y[i] * a_exp;

        g_y_0_0_0_z_z_z_z[i] = 2.0 * g_yz_z_z_z[i] * a_exp;
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_xz_x_x_x, g_xz_x_x_y, g_xz_x_x_z, g_z_0_0_0_x_x_x_x, g_z_0_0_0_x_x_x_y, g_z_0_0_0_x_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_x_x_x[i] = 2.0 * g_xz_x_x_x[i] * a_exp;

        g_z_0_0_0_x_x_x_y[i] = 2.0 * g_xz_x_x_y[i] * a_exp;

        g_z_0_0_0_x_x_x_z[i] = 2.0 * g_xz_x_x_z[i] * a_exp;
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_xz_x_y_x, g_xz_x_y_y, g_xz_x_y_z, g_z_0_0_0_x_x_y_x, g_z_0_0_0_x_x_y_y, g_z_0_0_0_x_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_x_y_x[i] = 2.0 * g_xz_x_y_x[i] * a_exp;

        g_z_0_0_0_x_x_y_y[i] = 2.0 * g_xz_x_y_y[i] * a_exp;

        g_z_0_0_0_x_x_y_z[i] = 2.0 * g_xz_x_y_z[i] * a_exp;
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_xz_x_z_x, g_xz_x_z_y, g_xz_x_z_z, g_z_0_0_0_x_x_z_x, g_z_0_0_0_x_x_z_y, g_z_0_0_0_x_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_x_z_x[i] = 2.0 * g_xz_x_z_x[i] * a_exp;

        g_z_0_0_0_x_x_z_y[i] = 2.0 * g_xz_x_z_y[i] * a_exp;

        g_z_0_0_0_x_x_z_z[i] = 2.0 * g_xz_x_z_z[i] * a_exp;
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_xz_y_x_x, g_xz_y_x_y, g_xz_y_x_z, g_z_0_0_0_x_y_x_x, g_z_0_0_0_x_y_x_y, g_z_0_0_0_x_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_y_x_x[i] = 2.0 * g_xz_y_x_x[i] * a_exp;

        g_z_0_0_0_x_y_x_y[i] = 2.0 * g_xz_y_x_y[i] * a_exp;

        g_z_0_0_0_x_y_x_z[i] = 2.0 * g_xz_y_x_z[i] * a_exp;
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_xz_y_y_x, g_xz_y_y_y, g_xz_y_y_z, g_z_0_0_0_x_y_y_x, g_z_0_0_0_x_y_y_y, g_z_0_0_0_x_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_y_y_x[i] = 2.0 * g_xz_y_y_x[i] * a_exp;

        g_z_0_0_0_x_y_y_y[i] = 2.0 * g_xz_y_y_y[i] * a_exp;

        g_z_0_0_0_x_y_y_z[i] = 2.0 * g_xz_y_y_z[i] * a_exp;
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_xz_y_z_x, g_xz_y_z_y, g_xz_y_z_z, g_z_0_0_0_x_y_z_x, g_z_0_0_0_x_y_z_y, g_z_0_0_0_x_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_y_z_x[i] = 2.0 * g_xz_y_z_x[i] * a_exp;

        g_z_0_0_0_x_y_z_y[i] = 2.0 * g_xz_y_z_y[i] * a_exp;

        g_z_0_0_0_x_y_z_z[i] = 2.0 * g_xz_y_z_z[i] * a_exp;
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_xz_z_x_x, g_xz_z_x_y, g_xz_z_x_z, g_z_0_0_0_x_z_x_x, g_z_0_0_0_x_z_x_y, g_z_0_0_0_x_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_z_x_x[i] = 2.0 * g_xz_z_x_x[i] * a_exp;

        g_z_0_0_0_x_z_x_y[i] = 2.0 * g_xz_z_x_y[i] * a_exp;

        g_z_0_0_0_x_z_x_z[i] = 2.0 * g_xz_z_x_z[i] * a_exp;
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_xz_z_y_x, g_xz_z_y_y, g_xz_z_y_z, g_z_0_0_0_x_z_y_x, g_z_0_0_0_x_z_y_y, g_z_0_0_0_x_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_z_y_x[i] = 2.0 * g_xz_z_y_x[i] * a_exp;

        g_z_0_0_0_x_z_y_y[i] = 2.0 * g_xz_z_y_y[i] * a_exp;

        g_z_0_0_0_x_z_y_z[i] = 2.0 * g_xz_z_y_z[i] * a_exp;
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_xz_z_z_x, g_xz_z_z_y, g_xz_z_z_z, g_z_0_0_0_x_z_z_x, g_z_0_0_0_x_z_z_y, g_z_0_0_0_x_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_z_z_x[i] = 2.0 * g_xz_z_z_x[i] * a_exp;

        g_z_0_0_0_x_z_z_y[i] = 2.0 * g_xz_z_z_y[i] * a_exp;

        g_z_0_0_0_x_z_z_z[i] = 2.0 * g_xz_z_z_z[i] * a_exp;
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_yz_x_x_x, g_yz_x_x_y, g_yz_x_x_z, g_z_0_0_0_y_x_x_x, g_z_0_0_0_y_x_x_y, g_z_0_0_0_y_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_x_x_x[i] = 2.0 * g_yz_x_x_x[i] * a_exp;

        g_z_0_0_0_y_x_x_y[i] = 2.0 * g_yz_x_x_y[i] * a_exp;

        g_z_0_0_0_y_x_x_z[i] = 2.0 * g_yz_x_x_z[i] * a_exp;
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_yz_x_y_x, g_yz_x_y_y, g_yz_x_y_z, g_z_0_0_0_y_x_y_x, g_z_0_0_0_y_x_y_y, g_z_0_0_0_y_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_x_y_x[i] = 2.0 * g_yz_x_y_x[i] * a_exp;

        g_z_0_0_0_y_x_y_y[i] = 2.0 * g_yz_x_y_y[i] * a_exp;

        g_z_0_0_0_y_x_y_z[i] = 2.0 * g_yz_x_y_z[i] * a_exp;
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_yz_x_z_x, g_yz_x_z_y, g_yz_x_z_z, g_z_0_0_0_y_x_z_x, g_z_0_0_0_y_x_z_y, g_z_0_0_0_y_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_x_z_x[i] = 2.0 * g_yz_x_z_x[i] * a_exp;

        g_z_0_0_0_y_x_z_y[i] = 2.0 * g_yz_x_z_y[i] * a_exp;

        g_z_0_0_0_y_x_z_z[i] = 2.0 * g_yz_x_z_z[i] * a_exp;
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_yz_y_x_x, g_yz_y_x_y, g_yz_y_x_z, g_z_0_0_0_y_y_x_x, g_z_0_0_0_y_y_x_y, g_z_0_0_0_y_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_y_x_x[i] = 2.0 * g_yz_y_x_x[i] * a_exp;

        g_z_0_0_0_y_y_x_y[i] = 2.0 * g_yz_y_x_y[i] * a_exp;

        g_z_0_0_0_y_y_x_z[i] = 2.0 * g_yz_y_x_z[i] * a_exp;
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_yz_y_y_x, g_yz_y_y_y, g_yz_y_y_z, g_z_0_0_0_y_y_y_x, g_z_0_0_0_y_y_y_y, g_z_0_0_0_y_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_y_y_x[i] = 2.0 * g_yz_y_y_x[i] * a_exp;

        g_z_0_0_0_y_y_y_y[i] = 2.0 * g_yz_y_y_y[i] * a_exp;

        g_z_0_0_0_y_y_y_z[i] = 2.0 * g_yz_y_y_z[i] * a_exp;
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_yz_y_z_x, g_yz_y_z_y, g_yz_y_z_z, g_z_0_0_0_y_y_z_x, g_z_0_0_0_y_y_z_y, g_z_0_0_0_y_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_y_z_x[i] = 2.0 * g_yz_y_z_x[i] * a_exp;

        g_z_0_0_0_y_y_z_y[i] = 2.0 * g_yz_y_z_y[i] * a_exp;

        g_z_0_0_0_y_y_z_z[i] = 2.0 * g_yz_y_z_z[i] * a_exp;
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_yz_z_x_x, g_yz_z_x_y, g_yz_z_x_z, g_z_0_0_0_y_z_x_x, g_z_0_0_0_y_z_x_y, g_z_0_0_0_y_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_z_x_x[i] = 2.0 * g_yz_z_x_x[i] * a_exp;

        g_z_0_0_0_y_z_x_y[i] = 2.0 * g_yz_z_x_y[i] * a_exp;

        g_z_0_0_0_y_z_x_z[i] = 2.0 * g_yz_z_x_z[i] * a_exp;
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_yz_z_y_x, g_yz_z_y_y, g_yz_z_y_z, g_z_0_0_0_y_z_y_x, g_z_0_0_0_y_z_y_y, g_z_0_0_0_y_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_z_y_x[i] = 2.0 * g_yz_z_y_x[i] * a_exp;

        g_z_0_0_0_y_z_y_y[i] = 2.0 * g_yz_z_y_y[i] * a_exp;

        g_z_0_0_0_y_z_y_z[i] = 2.0 * g_yz_z_y_z[i] * a_exp;
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_yz_z_z_x, g_yz_z_z_y, g_yz_z_z_z, g_z_0_0_0_y_z_z_x, g_z_0_0_0_y_z_z_y, g_z_0_0_0_y_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_z_z_x[i] = 2.0 * g_yz_z_z_x[i] * a_exp;

        g_z_0_0_0_y_z_z_y[i] = 2.0 * g_yz_z_z_y[i] * a_exp;

        g_z_0_0_0_y_z_z_z[i] = 2.0 * g_yz_z_z_z[i] * a_exp;
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_0_x_x_x, g_0_x_x_y, g_0_x_x_z, g_z_0_0_0_z_x_x_x, g_z_0_0_0_z_x_x_y, g_z_0_0_0_z_x_x_z, g_zz_x_x_x, g_zz_x_x_y, g_zz_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_x_x_x[i] = -g_0_x_x_x[i] + 2.0 * g_zz_x_x_x[i] * a_exp;

        g_z_0_0_0_z_x_x_y[i] = -g_0_x_x_y[i] + 2.0 * g_zz_x_x_y[i] * a_exp;

        g_z_0_0_0_z_x_x_z[i] = -g_0_x_x_z[i] + 2.0 * g_zz_x_x_z[i] * a_exp;
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_0_x_y_x, g_0_x_y_y, g_0_x_y_z, g_z_0_0_0_z_x_y_x, g_z_0_0_0_z_x_y_y, g_z_0_0_0_z_x_y_z, g_zz_x_y_x, g_zz_x_y_y, g_zz_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_x_y_x[i] = -g_0_x_y_x[i] + 2.0 * g_zz_x_y_x[i] * a_exp;

        g_z_0_0_0_z_x_y_y[i] = -g_0_x_y_y[i] + 2.0 * g_zz_x_y_y[i] * a_exp;

        g_z_0_0_0_z_x_y_z[i] = -g_0_x_y_z[i] + 2.0 * g_zz_x_y_z[i] * a_exp;
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_0_x_z_x, g_0_x_z_y, g_0_x_z_z, g_z_0_0_0_z_x_z_x, g_z_0_0_0_z_x_z_y, g_z_0_0_0_z_x_z_z, g_zz_x_z_x, g_zz_x_z_y, g_zz_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_x_z_x[i] = -g_0_x_z_x[i] + 2.0 * g_zz_x_z_x[i] * a_exp;

        g_z_0_0_0_z_x_z_y[i] = -g_0_x_z_y[i] + 2.0 * g_zz_x_z_y[i] * a_exp;

        g_z_0_0_0_z_x_z_z[i] = -g_0_x_z_z[i] + 2.0 * g_zz_x_z_z[i] * a_exp;
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_0_y_x_x, g_0_y_x_y, g_0_y_x_z, g_z_0_0_0_z_y_x_x, g_z_0_0_0_z_y_x_y, g_z_0_0_0_z_y_x_z, g_zz_y_x_x, g_zz_y_x_y, g_zz_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_y_x_x[i] = -g_0_y_x_x[i] + 2.0 * g_zz_y_x_x[i] * a_exp;

        g_z_0_0_0_z_y_x_y[i] = -g_0_y_x_y[i] + 2.0 * g_zz_y_x_y[i] * a_exp;

        g_z_0_0_0_z_y_x_z[i] = -g_0_y_x_z[i] + 2.0 * g_zz_y_x_z[i] * a_exp;
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_0_y_y_x, g_0_y_y_y, g_0_y_y_z, g_z_0_0_0_z_y_y_x, g_z_0_0_0_z_y_y_y, g_z_0_0_0_z_y_y_z, g_zz_y_y_x, g_zz_y_y_y, g_zz_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_y_y_x[i] = -g_0_y_y_x[i] + 2.0 * g_zz_y_y_x[i] * a_exp;

        g_z_0_0_0_z_y_y_y[i] = -g_0_y_y_y[i] + 2.0 * g_zz_y_y_y[i] * a_exp;

        g_z_0_0_0_z_y_y_z[i] = -g_0_y_y_z[i] + 2.0 * g_zz_y_y_z[i] * a_exp;
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_0_y_z_x, g_0_y_z_y, g_0_y_z_z, g_z_0_0_0_z_y_z_x, g_z_0_0_0_z_y_z_y, g_z_0_0_0_z_y_z_z, g_zz_y_z_x, g_zz_y_z_y, g_zz_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_y_z_x[i] = -g_0_y_z_x[i] + 2.0 * g_zz_y_z_x[i] * a_exp;

        g_z_0_0_0_z_y_z_y[i] = -g_0_y_z_y[i] + 2.0 * g_zz_y_z_y[i] * a_exp;

        g_z_0_0_0_z_y_z_z[i] = -g_0_y_z_z[i] + 2.0 * g_zz_y_z_z[i] * a_exp;
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_0_z_x_x, g_0_z_x_y, g_0_z_x_z, g_z_0_0_0_z_z_x_x, g_z_0_0_0_z_z_x_y, g_z_0_0_0_z_z_x_z, g_zz_z_x_x, g_zz_z_x_y, g_zz_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_z_x_x[i] = -g_0_z_x_x[i] + 2.0 * g_zz_z_x_x[i] * a_exp;

        g_z_0_0_0_z_z_x_y[i] = -g_0_z_x_y[i] + 2.0 * g_zz_z_x_y[i] * a_exp;

        g_z_0_0_0_z_z_x_z[i] = -g_0_z_x_z[i] + 2.0 * g_zz_z_x_z[i] * a_exp;
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_0_z_y_x, g_0_z_y_y, g_0_z_y_z, g_z_0_0_0_z_z_y_x, g_z_0_0_0_z_z_y_y, g_z_0_0_0_z_z_y_z, g_zz_z_y_x, g_zz_z_y_y, g_zz_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_z_y_x[i] = -g_0_z_y_x[i] + 2.0 * g_zz_z_y_x[i] * a_exp;

        g_z_0_0_0_z_z_y_y[i] = -g_0_z_y_y[i] + 2.0 * g_zz_z_y_y[i] * a_exp;

        g_z_0_0_0_z_z_y_z[i] = -g_0_z_y_z[i] + 2.0 * g_zz_z_y_z[i] * a_exp;
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_0_z_z_x, g_0_z_z_y, g_0_z_z_z, g_z_0_0_0_z_z_z_x, g_z_0_0_0_z_z_z_y, g_z_0_0_0_z_z_z_z, g_zz_z_z_x, g_zz_z_z_y, g_zz_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_z_z_x[i] = -g_0_z_z_x[i] + 2.0 * g_zz_z_z_x[i] * a_exp;

        g_z_0_0_0_z_z_z_y[i] = -g_0_z_z_y[i] + 2.0 * g_zz_z_z_y[i] * a_exp;

        g_z_0_0_0_z_z_z_z[i] = -g_0_z_z_z[i] + 2.0 * g_zz_z_z_z[i] * a_exp;
    }
}

} // t4c_geom namespace

