#include "GeomDeriv1010OfScalarForSPPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_sppp_0(CSimdArray<double>& buffer_1010_sppp,
                     const CSimdArray<double>& buffer_ppsp,
                     const CSimdArray<double>& buffer_ppdp,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_sppp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_ppsp

    auto g_x_x_0_x = buffer_ppsp[0];

    auto g_x_x_0_y = buffer_ppsp[1];

    auto g_x_x_0_z = buffer_ppsp[2];

    auto g_x_y_0_x = buffer_ppsp[3];

    auto g_x_y_0_y = buffer_ppsp[4];

    auto g_x_y_0_z = buffer_ppsp[5];

    auto g_x_z_0_x = buffer_ppsp[6];

    auto g_x_z_0_y = buffer_ppsp[7];

    auto g_x_z_0_z = buffer_ppsp[8];

    auto g_y_x_0_x = buffer_ppsp[9];

    auto g_y_x_0_y = buffer_ppsp[10];

    auto g_y_x_0_z = buffer_ppsp[11];

    auto g_y_y_0_x = buffer_ppsp[12];

    auto g_y_y_0_y = buffer_ppsp[13];

    auto g_y_y_0_z = buffer_ppsp[14];

    auto g_y_z_0_x = buffer_ppsp[15];

    auto g_y_z_0_y = buffer_ppsp[16];

    auto g_y_z_0_z = buffer_ppsp[17];

    auto g_z_x_0_x = buffer_ppsp[18];

    auto g_z_x_0_y = buffer_ppsp[19];

    auto g_z_x_0_z = buffer_ppsp[20];

    auto g_z_y_0_x = buffer_ppsp[21];

    auto g_z_y_0_y = buffer_ppsp[22];

    auto g_z_y_0_z = buffer_ppsp[23];

    auto g_z_z_0_x = buffer_ppsp[24];

    auto g_z_z_0_y = buffer_ppsp[25];

    auto g_z_z_0_z = buffer_ppsp[26];

    /// Set up components of auxilary buffer : buffer_ppdp

    auto g_x_x_xx_x = buffer_ppdp[0];

    auto g_x_x_xx_y = buffer_ppdp[1];

    auto g_x_x_xx_z = buffer_ppdp[2];

    auto g_x_x_xy_x = buffer_ppdp[3];

    auto g_x_x_xy_y = buffer_ppdp[4];

    auto g_x_x_xy_z = buffer_ppdp[5];

    auto g_x_x_xz_x = buffer_ppdp[6];

    auto g_x_x_xz_y = buffer_ppdp[7];

    auto g_x_x_xz_z = buffer_ppdp[8];

    auto g_x_x_yy_x = buffer_ppdp[9];

    auto g_x_x_yy_y = buffer_ppdp[10];

    auto g_x_x_yy_z = buffer_ppdp[11];

    auto g_x_x_yz_x = buffer_ppdp[12];

    auto g_x_x_yz_y = buffer_ppdp[13];

    auto g_x_x_yz_z = buffer_ppdp[14];

    auto g_x_x_zz_x = buffer_ppdp[15];

    auto g_x_x_zz_y = buffer_ppdp[16];

    auto g_x_x_zz_z = buffer_ppdp[17];

    auto g_x_y_xx_x = buffer_ppdp[18];

    auto g_x_y_xx_y = buffer_ppdp[19];

    auto g_x_y_xx_z = buffer_ppdp[20];

    auto g_x_y_xy_x = buffer_ppdp[21];

    auto g_x_y_xy_y = buffer_ppdp[22];

    auto g_x_y_xy_z = buffer_ppdp[23];

    auto g_x_y_xz_x = buffer_ppdp[24];

    auto g_x_y_xz_y = buffer_ppdp[25];

    auto g_x_y_xz_z = buffer_ppdp[26];

    auto g_x_y_yy_x = buffer_ppdp[27];

    auto g_x_y_yy_y = buffer_ppdp[28];

    auto g_x_y_yy_z = buffer_ppdp[29];

    auto g_x_y_yz_x = buffer_ppdp[30];

    auto g_x_y_yz_y = buffer_ppdp[31];

    auto g_x_y_yz_z = buffer_ppdp[32];

    auto g_x_y_zz_x = buffer_ppdp[33];

    auto g_x_y_zz_y = buffer_ppdp[34];

    auto g_x_y_zz_z = buffer_ppdp[35];

    auto g_x_z_xx_x = buffer_ppdp[36];

    auto g_x_z_xx_y = buffer_ppdp[37];

    auto g_x_z_xx_z = buffer_ppdp[38];

    auto g_x_z_xy_x = buffer_ppdp[39];

    auto g_x_z_xy_y = buffer_ppdp[40];

    auto g_x_z_xy_z = buffer_ppdp[41];

    auto g_x_z_xz_x = buffer_ppdp[42];

    auto g_x_z_xz_y = buffer_ppdp[43];

    auto g_x_z_xz_z = buffer_ppdp[44];

    auto g_x_z_yy_x = buffer_ppdp[45];

    auto g_x_z_yy_y = buffer_ppdp[46];

    auto g_x_z_yy_z = buffer_ppdp[47];

    auto g_x_z_yz_x = buffer_ppdp[48];

    auto g_x_z_yz_y = buffer_ppdp[49];

    auto g_x_z_yz_z = buffer_ppdp[50];

    auto g_x_z_zz_x = buffer_ppdp[51];

    auto g_x_z_zz_y = buffer_ppdp[52];

    auto g_x_z_zz_z = buffer_ppdp[53];

    auto g_y_x_xx_x = buffer_ppdp[54];

    auto g_y_x_xx_y = buffer_ppdp[55];

    auto g_y_x_xx_z = buffer_ppdp[56];

    auto g_y_x_xy_x = buffer_ppdp[57];

    auto g_y_x_xy_y = buffer_ppdp[58];

    auto g_y_x_xy_z = buffer_ppdp[59];

    auto g_y_x_xz_x = buffer_ppdp[60];

    auto g_y_x_xz_y = buffer_ppdp[61];

    auto g_y_x_xz_z = buffer_ppdp[62];

    auto g_y_x_yy_x = buffer_ppdp[63];

    auto g_y_x_yy_y = buffer_ppdp[64];

    auto g_y_x_yy_z = buffer_ppdp[65];

    auto g_y_x_yz_x = buffer_ppdp[66];

    auto g_y_x_yz_y = buffer_ppdp[67];

    auto g_y_x_yz_z = buffer_ppdp[68];

    auto g_y_x_zz_x = buffer_ppdp[69];

    auto g_y_x_zz_y = buffer_ppdp[70];

    auto g_y_x_zz_z = buffer_ppdp[71];

    auto g_y_y_xx_x = buffer_ppdp[72];

    auto g_y_y_xx_y = buffer_ppdp[73];

    auto g_y_y_xx_z = buffer_ppdp[74];

    auto g_y_y_xy_x = buffer_ppdp[75];

    auto g_y_y_xy_y = buffer_ppdp[76];

    auto g_y_y_xy_z = buffer_ppdp[77];

    auto g_y_y_xz_x = buffer_ppdp[78];

    auto g_y_y_xz_y = buffer_ppdp[79];

    auto g_y_y_xz_z = buffer_ppdp[80];

    auto g_y_y_yy_x = buffer_ppdp[81];

    auto g_y_y_yy_y = buffer_ppdp[82];

    auto g_y_y_yy_z = buffer_ppdp[83];

    auto g_y_y_yz_x = buffer_ppdp[84];

    auto g_y_y_yz_y = buffer_ppdp[85];

    auto g_y_y_yz_z = buffer_ppdp[86];

    auto g_y_y_zz_x = buffer_ppdp[87];

    auto g_y_y_zz_y = buffer_ppdp[88];

    auto g_y_y_zz_z = buffer_ppdp[89];

    auto g_y_z_xx_x = buffer_ppdp[90];

    auto g_y_z_xx_y = buffer_ppdp[91];

    auto g_y_z_xx_z = buffer_ppdp[92];

    auto g_y_z_xy_x = buffer_ppdp[93];

    auto g_y_z_xy_y = buffer_ppdp[94];

    auto g_y_z_xy_z = buffer_ppdp[95];

    auto g_y_z_xz_x = buffer_ppdp[96];

    auto g_y_z_xz_y = buffer_ppdp[97];

    auto g_y_z_xz_z = buffer_ppdp[98];

    auto g_y_z_yy_x = buffer_ppdp[99];

    auto g_y_z_yy_y = buffer_ppdp[100];

    auto g_y_z_yy_z = buffer_ppdp[101];

    auto g_y_z_yz_x = buffer_ppdp[102];

    auto g_y_z_yz_y = buffer_ppdp[103];

    auto g_y_z_yz_z = buffer_ppdp[104];

    auto g_y_z_zz_x = buffer_ppdp[105];

    auto g_y_z_zz_y = buffer_ppdp[106];

    auto g_y_z_zz_z = buffer_ppdp[107];

    auto g_z_x_xx_x = buffer_ppdp[108];

    auto g_z_x_xx_y = buffer_ppdp[109];

    auto g_z_x_xx_z = buffer_ppdp[110];

    auto g_z_x_xy_x = buffer_ppdp[111];

    auto g_z_x_xy_y = buffer_ppdp[112];

    auto g_z_x_xy_z = buffer_ppdp[113];

    auto g_z_x_xz_x = buffer_ppdp[114];

    auto g_z_x_xz_y = buffer_ppdp[115];

    auto g_z_x_xz_z = buffer_ppdp[116];

    auto g_z_x_yy_x = buffer_ppdp[117];

    auto g_z_x_yy_y = buffer_ppdp[118];

    auto g_z_x_yy_z = buffer_ppdp[119];

    auto g_z_x_yz_x = buffer_ppdp[120];

    auto g_z_x_yz_y = buffer_ppdp[121];

    auto g_z_x_yz_z = buffer_ppdp[122];

    auto g_z_x_zz_x = buffer_ppdp[123];

    auto g_z_x_zz_y = buffer_ppdp[124];

    auto g_z_x_zz_z = buffer_ppdp[125];

    auto g_z_y_xx_x = buffer_ppdp[126];

    auto g_z_y_xx_y = buffer_ppdp[127];

    auto g_z_y_xx_z = buffer_ppdp[128];

    auto g_z_y_xy_x = buffer_ppdp[129];

    auto g_z_y_xy_y = buffer_ppdp[130];

    auto g_z_y_xy_z = buffer_ppdp[131];

    auto g_z_y_xz_x = buffer_ppdp[132];

    auto g_z_y_xz_y = buffer_ppdp[133];

    auto g_z_y_xz_z = buffer_ppdp[134];

    auto g_z_y_yy_x = buffer_ppdp[135];

    auto g_z_y_yy_y = buffer_ppdp[136];

    auto g_z_y_yy_z = buffer_ppdp[137];

    auto g_z_y_yz_x = buffer_ppdp[138];

    auto g_z_y_yz_y = buffer_ppdp[139];

    auto g_z_y_yz_z = buffer_ppdp[140];

    auto g_z_y_zz_x = buffer_ppdp[141];

    auto g_z_y_zz_y = buffer_ppdp[142];

    auto g_z_y_zz_z = buffer_ppdp[143];

    auto g_z_z_xx_x = buffer_ppdp[144];

    auto g_z_z_xx_y = buffer_ppdp[145];

    auto g_z_z_xx_z = buffer_ppdp[146];

    auto g_z_z_xy_x = buffer_ppdp[147];

    auto g_z_z_xy_y = buffer_ppdp[148];

    auto g_z_z_xy_z = buffer_ppdp[149];

    auto g_z_z_xz_x = buffer_ppdp[150];

    auto g_z_z_xz_y = buffer_ppdp[151];

    auto g_z_z_xz_z = buffer_ppdp[152];

    auto g_z_z_yy_x = buffer_ppdp[153];

    auto g_z_z_yy_y = buffer_ppdp[154];

    auto g_z_z_yy_z = buffer_ppdp[155];

    auto g_z_z_yz_x = buffer_ppdp[156];

    auto g_z_z_yz_y = buffer_ppdp[157];

    auto g_z_z_yz_z = buffer_ppdp[158];

    auto g_z_z_zz_x = buffer_ppdp[159];

    auto g_z_z_zz_y = buffer_ppdp[160];

    auto g_z_z_zz_z = buffer_ppdp[161];

    /// Set up components of integrals buffer : buffer_1010_sppp

    auto g_x_0_x_0_0_x_x_x = buffer_1010_sppp[0];

    auto g_x_0_x_0_0_x_x_y = buffer_1010_sppp[1];

    auto g_x_0_x_0_0_x_x_z = buffer_1010_sppp[2];

    auto g_x_0_x_0_0_x_y_x = buffer_1010_sppp[3];

    auto g_x_0_x_0_0_x_y_y = buffer_1010_sppp[4];

    auto g_x_0_x_0_0_x_y_z = buffer_1010_sppp[5];

    auto g_x_0_x_0_0_x_z_x = buffer_1010_sppp[6];

    auto g_x_0_x_0_0_x_z_y = buffer_1010_sppp[7];

    auto g_x_0_x_0_0_x_z_z = buffer_1010_sppp[8];

    auto g_x_0_x_0_0_y_x_x = buffer_1010_sppp[9];

    auto g_x_0_x_0_0_y_x_y = buffer_1010_sppp[10];

    auto g_x_0_x_0_0_y_x_z = buffer_1010_sppp[11];

    auto g_x_0_x_0_0_y_y_x = buffer_1010_sppp[12];

    auto g_x_0_x_0_0_y_y_y = buffer_1010_sppp[13];

    auto g_x_0_x_0_0_y_y_z = buffer_1010_sppp[14];

    auto g_x_0_x_0_0_y_z_x = buffer_1010_sppp[15];

    auto g_x_0_x_0_0_y_z_y = buffer_1010_sppp[16];

    auto g_x_0_x_0_0_y_z_z = buffer_1010_sppp[17];

    auto g_x_0_x_0_0_z_x_x = buffer_1010_sppp[18];

    auto g_x_0_x_0_0_z_x_y = buffer_1010_sppp[19];

    auto g_x_0_x_0_0_z_x_z = buffer_1010_sppp[20];

    auto g_x_0_x_0_0_z_y_x = buffer_1010_sppp[21];

    auto g_x_0_x_0_0_z_y_y = buffer_1010_sppp[22];

    auto g_x_0_x_0_0_z_y_z = buffer_1010_sppp[23];

    auto g_x_0_x_0_0_z_z_x = buffer_1010_sppp[24];

    auto g_x_0_x_0_0_z_z_y = buffer_1010_sppp[25];

    auto g_x_0_x_0_0_z_z_z = buffer_1010_sppp[26];

    auto g_x_0_y_0_0_x_x_x = buffer_1010_sppp[27];

    auto g_x_0_y_0_0_x_x_y = buffer_1010_sppp[28];

    auto g_x_0_y_0_0_x_x_z = buffer_1010_sppp[29];

    auto g_x_0_y_0_0_x_y_x = buffer_1010_sppp[30];

    auto g_x_0_y_0_0_x_y_y = buffer_1010_sppp[31];

    auto g_x_0_y_0_0_x_y_z = buffer_1010_sppp[32];

    auto g_x_0_y_0_0_x_z_x = buffer_1010_sppp[33];

    auto g_x_0_y_0_0_x_z_y = buffer_1010_sppp[34];

    auto g_x_0_y_0_0_x_z_z = buffer_1010_sppp[35];

    auto g_x_0_y_0_0_y_x_x = buffer_1010_sppp[36];

    auto g_x_0_y_0_0_y_x_y = buffer_1010_sppp[37];

    auto g_x_0_y_0_0_y_x_z = buffer_1010_sppp[38];

    auto g_x_0_y_0_0_y_y_x = buffer_1010_sppp[39];

    auto g_x_0_y_0_0_y_y_y = buffer_1010_sppp[40];

    auto g_x_0_y_0_0_y_y_z = buffer_1010_sppp[41];

    auto g_x_0_y_0_0_y_z_x = buffer_1010_sppp[42];

    auto g_x_0_y_0_0_y_z_y = buffer_1010_sppp[43];

    auto g_x_0_y_0_0_y_z_z = buffer_1010_sppp[44];

    auto g_x_0_y_0_0_z_x_x = buffer_1010_sppp[45];

    auto g_x_0_y_0_0_z_x_y = buffer_1010_sppp[46];

    auto g_x_0_y_0_0_z_x_z = buffer_1010_sppp[47];

    auto g_x_0_y_0_0_z_y_x = buffer_1010_sppp[48];

    auto g_x_0_y_0_0_z_y_y = buffer_1010_sppp[49];

    auto g_x_0_y_0_0_z_y_z = buffer_1010_sppp[50];

    auto g_x_0_y_0_0_z_z_x = buffer_1010_sppp[51];

    auto g_x_0_y_0_0_z_z_y = buffer_1010_sppp[52];

    auto g_x_0_y_0_0_z_z_z = buffer_1010_sppp[53];

    auto g_x_0_z_0_0_x_x_x = buffer_1010_sppp[54];

    auto g_x_0_z_0_0_x_x_y = buffer_1010_sppp[55];

    auto g_x_0_z_0_0_x_x_z = buffer_1010_sppp[56];

    auto g_x_0_z_0_0_x_y_x = buffer_1010_sppp[57];

    auto g_x_0_z_0_0_x_y_y = buffer_1010_sppp[58];

    auto g_x_0_z_0_0_x_y_z = buffer_1010_sppp[59];

    auto g_x_0_z_0_0_x_z_x = buffer_1010_sppp[60];

    auto g_x_0_z_0_0_x_z_y = buffer_1010_sppp[61];

    auto g_x_0_z_0_0_x_z_z = buffer_1010_sppp[62];

    auto g_x_0_z_0_0_y_x_x = buffer_1010_sppp[63];

    auto g_x_0_z_0_0_y_x_y = buffer_1010_sppp[64];

    auto g_x_0_z_0_0_y_x_z = buffer_1010_sppp[65];

    auto g_x_0_z_0_0_y_y_x = buffer_1010_sppp[66];

    auto g_x_0_z_0_0_y_y_y = buffer_1010_sppp[67];

    auto g_x_0_z_0_0_y_y_z = buffer_1010_sppp[68];

    auto g_x_0_z_0_0_y_z_x = buffer_1010_sppp[69];

    auto g_x_0_z_0_0_y_z_y = buffer_1010_sppp[70];

    auto g_x_0_z_0_0_y_z_z = buffer_1010_sppp[71];

    auto g_x_0_z_0_0_z_x_x = buffer_1010_sppp[72];

    auto g_x_0_z_0_0_z_x_y = buffer_1010_sppp[73];

    auto g_x_0_z_0_0_z_x_z = buffer_1010_sppp[74];

    auto g_x_0_z_0_0_z_y_x = buffer_1010_sppp[75];

    auto g_x_0_z_0_0_z_y_y = buffer_1010_sppp[76];

    auto g_x_0_z_0_0_z_y_z = buffer_1010_sppp[77];

    auto g_x_0_z_0_0_z_z_x = buffer_1010_sppp[78];

    auto g_x_0_z_0_0_z_z_y = buffer_1010_sppp[79];

    auto g_x_0_z_0_0_z_z_z = buffer_1010_sppp[80];

    auto g_y_0_x_0_0_x_x_x = buffer_1010_sppp[81];

    auto g_y_0_x_0_0_x_x_y = buffer_1010_sppp[82];

    auto g_y_0_x_0_0_x_x_z = buffer_1010_sppp[83];

    auto g_y_0_x_0_0_x_y_x = buffer_1010_sppp[84];

    auto g_y_0_x_0_0_x_y_y = buffer_1010_sppp[85];

    auto g_y_0_x_0_0_x_y_z = buffer_1010_sppp[86];

    auto g_y_0_x_0_0_x_z_x = buffer_1010_sppp[87];

    auto g_y_0_x_0_0_x_z_y = buffer_1010_sppp[88];

    auto g_y_0_x_0_0_x_z_z = buffer_1010_sppp[89];

    auto g_y_0_x_0_0_y_x_x = buffer_1010_sppp[90];

    auto g_y_0_x_0_0_y_x_y = buffer_1010_sppp[91];

    auto g_y_0_x_0_0_y_x_z = buffer_1010_sppp[92];

    auto g_y_0_x_0_0_y_y_x = buffer_1010_sppp[93];

    auto g_y_0_x_0_0_y_y_y = buffer_1010_sppp[94];

    auto g_y_0_x_0_0_y_y_z = buffer_1010_sppp[95];

    auto g_y_0_x_0_0_y_z_x = buffer_1010_sppp[96];

    auto g_y_0_x_0_0_y_z_y = buffer_1010_sppp[97];

    auto g_y_0_x_0_0_y_z_z = buffer_1010_sppp[98];

    auto g_y_0_x_0_0_z_x_x = buffer_1010_sppp[99];

    auto g_y_0_x_0_0_z_x_y = buffer_1010_sppp[100];

    auto g_y_0_x_0_0_z_x_z = buffer_1010_sppp[101];

    auto g_y_0_x_0_0_z_y_x = buffer_1010_sppp[102];

    auto g_y_0_x_0_0_z_y_y = buffer_1010_sppp[103];

    auto g_y_0_x_0_0_z_y_z = buffer_1010_sppp[104];

    auto g_y_0_x_0_0_z_z_x = buffer_1010_sppp[105];

    auto g_y_0_x_0_0_z_z_y = buffer_1010_sppp[106];

    auto g_y_0_x_0_0_z_z_z = buffer_1010_sppp[107];

    auto g_y_0_y_0_0_x_x_x = buffer_1010_sppp[108];

    auto g_y_0_y_0_0_x_x_y = buffer_1010_sppp[109];

    auto g_y_0_y_0_0_x_x_z = buffer_1010_sppp[110];

    auto g_y_0_y_0_0_x_y_x = buffer_1010_sppp[111];

    auto g_y_0_y_0_0_x_y_y = buffer_1010_sppp[112];

    auto g_y_0_y_0_0_x_y_z = buffer_1010_sppp[113];

    auto g_y_0_y_0_0_x_z_x = buffer_1010_sppp[114];

    auto g_y_0_y_0_0_x_z_y = buffer_1010_sppp[115];

    auto g_y_0_y_0_0_x_z_z = buffer_1010_sppp[116];

    auto g_y_0_y_0_0_y_x_x = buffer_1010_sppp[117];

    auto g_y_0_y_0_0_y_x_y = buffer_1010_sppp[118];

    auto g_y_0_y_0_0_y_x_z = buffer_1010_sppp[119];

    auto g_y_0_y_0_0_y_y_x = buffer_1010_sppp[120];

    auto g_y_0_y_0_0_y_y_y = buffer_1010_sppp[121];

    auto g_y_0_y_0_0_y_y_z = buffer_1010_sppp[122];

    auto g_y_0_y_0_0_y_z_x = buffer_1010_sppp[123];

    auto g_y_0_y_0_0_y_z_y = buffer_1010_sppp[124];

    auto g_y_0_y_0_0_y_z_z = buffer_1010_sppp[125];

    auto g_y_0_y_0_0_z_x_x = buffer_1010_sppp[126];

    auto g_y_0_y_0_0_z_x_y = buffer_1010_sppp[127];

    auto g_y_0_y_0_0_z_x_z = buffer_1010_sppp[128];

    auto g_y_0_y_0_0_z_y_x = buffer_1010_sppp[129];

    auto g_y_0_y_0_0_z_y_y = buffer_1010_sppp[130];

    auto g_y_0_y_0_0_z_y_z = buffer_1010_sppp[131];

    auto g_y_0_y_0_0_z_z_x = buffer_1010_sppp[132];

    auto g_y_0_y_0_0_z_z_y = buffer_1010_sppp[133];

    auto g_y_0_y_0_0_z_z_z = buffer_1010_sppp[134];

    auto g_y_0_z_0_0_x_x_x = buffer_1010_sppp[135];

    auto g_y_0_z_0_0_x_x_y = buffer_1010_sppp[136];

    auto g_y_0_z_0_0_x_x_z = buffer_1010_sppp[137];

    auto g_y_0_z_0_0_x_y_x = buffer_1010_sppp[138];

    auto g_y_0_z_0_0_x_y_y = buffer_1010_sppp[139];

    auto g_y_0_z_0_0_x_y_z = buffer_1010_sppp[140];

    auto g_y_0_z_0_0_x_z_x = buffer_1010_sppp[141];

    auto g_y_0_z_0_0_x_z_y = buffer_1010_sppp[142];

    auto g_y_0_z_0_0_x_z_z = buffer_1010_sppp[143];

    auto g_y_0_z_0_0_y_x_x = buffer_1010_sppp[144];

    auto g_y_0_z_0_0_y_x_y = buffer_1010_sppp[145];

    auto g_y_0_z_0_0_y_x_z = buffer_1010_sppp[146];

    auto g_y_0_z_0_0_y_y_x = buffer_1010_sppp[147];

    auto g_y_0_z_0_0_y_y_y = buffer_1010_sppp[148];

    auto g_y_0_z_0_0_y_y_z = buffer_1010_sppp[149];

    auto g_y_0_z_0_0_y_z_x = buffer_1010_sppp[150];

    auto g_y_0_z_0_0_y_z_y = buffer_1010_sppp[151];

    auto g_y_0_z_0_0_y_z_z = buffer_1010_sppp[152];

    auto g_y_0_z_0_0_z_x_x = buffer_1010_sppp[153];

    auto g_y_0_z_0_0_z_x_y = buffer_1010_sppp[154];

    auto g_y_0_z_0_0_z_x_z = buffer_1010_sppp[155];

    auto g_y_0_z_0_0_z_y_x = buffer_1010_sppp[156];

    auto g_y_0_z_0_0_z_y_y = buffer_1010_sppp[157];

    auto g_y_0_z_0_0_z_y_z = buffer_1010_sppp[158];

    auto g_y_0_z_0_0_z_z_x = buffer_1010_sppp[159];

    auto g_y_0_z_0_0_z_z_y = buffer_1010_sppp[160];

    auto g_y_0_z_0_0_z_z_z = buffer_1010_sppp[161];

    auto g_z_0_x_0_0_x_x_x = buffer_1010_sppp[162];

    auto g_z_0_x_0_0_x_x_y = buffer_1010_sppp[163];

    auto g_z_0_x_0_0_x_x_z = buffer_1010_sppp[164];

    auto g_z_0_x_0_0_x_y_x = buffer_1010_sppp[165];

    auto g_z_0_x_0_0_x_y_y = buffer_1010_sppp[166];

    auto g_z_0_x_0_0_x_y_z = buffer_1010_sppp[167];

    auto g_z_0_x_0_0_x_z_x = buffer_1010_sppp[168];

    auto g_z_0_x_0_0_x_z_y = buffer_1010_sppp[169];

    auto g_z_0_x_0_0_x_z_z = buffer_1010_sppp[170];

    auto g_z_0_x_0_0_y_x_x = buffer_1010_sppp[171];

    auto g_z_0_x_0_0_y_x_y = buffer_1010_sppp[172];

    auto g_z_0_x_0_0_y_x_z = buffer_1010_sppp[173];

    auto g_z_0_x_0_0_y_y_x = buffer_1010_sppp[174];

    auto g_z_0_x_0_0_y_y_y = buffer_1010_sppp[175];

    auto g_z_0_x_0_0_y_y_z = buffer_1010_sppp[176];

    auto g_z_0_x_0_0_y_z_x = buffer_1010_sppp[177];

    auto g_z_0_x_0_0_y_z_y = buffer_1010_sppp[178];

    auto g_z_0_x_0_0_y_z_z = buffer_1010_sppp[179];

    auto g_z_0_x_0_0_z_x_x = buffer_1010_sppp[180];

    auto g_z_0_x_0_0_z_x_y = buffer_1010_sppp[181];

    auto g_z_0_x_0_0_z_x_z = buffer_1010_sppp[182];

    auto g_z_0_x_0_0_z_y_x = buffer_1010_sppp[183];

    auto g_z_0_x_0_0_z_y_y = buffer_1010_sppp[184];

    auto g_z_0_x_0_0_z_y_z = buffer_1010_sppp[185];

    auto g_z_0_x_0_0_z_z_x = buffer_1010_sppp[186];

    auto g_z_0_x_0_0_z_z_y = buffer_1010_sppp[187];

    auto g_z_0_x_0_0_z_z_z = buffer_1010_sppp[188];

    auto g_z_0_y_0_0_x_x_x = buffer_1010_sppp[189];

    auto g_z_0_y_0_0_x_x_y = buffer_1010_sppp[190];

    auto g_z_0_y_0_0_x_x_z = buffer_1010_sppp[191];

    auto g_z_0_y_0_0_x_y_x = buffer_1010_sppp[192];

    auto g_z_0_y_0_0_x_y_y = buffer_1010_sppp[193];

    auto g_z_0_y_0_0_x_y_z = buffer_1010_sppp[194];

    auto g_z_0_y_0_0_x_z_x = buffer_1010_sppp[195];

    auto g_z_0_y_0_0_x_z_y = buffer_1010_sppp[196];

    auto g_z_0_y_0_0_x_z_z = buffer_1010_sppp[197];

    auto g_z_0_y_0_0_y_x_x = buffer_1010_sppp[198];

    auto g_z_0_y_0_0_y_x_y = buffer_1010_sppp[199];

    auto g_z_0_y_0_0_y_x_z = buffer_1010_sppp[200];

    auto g_z_0_y_0_0_y_y_x = buffer_1010_sppp[201];

    auto g_z_0_y_0_0_y_y_y = buffer_1010_sppp[202];

    auto g_z_0_y_0_0_y_y_z = buffer_1010_sppp[203];

    auto g_z_0_y_0_0_y_z_x = buffer_1010_sppp[204];

    auto g_z_0_y_0_0_y_z_y = buffer_1010_sppp[205];

    auto g_z_0_y_0_0_y_z_z = buffer_1010_sppp[206];

    auto g_z_0_y_0_0_z_x_x = buffer_1010_sppp[207];

    auto g_z_0_y_0_0_z_x_y = buffer_1010_sppp[208];

    auto g_z_0_y_0_0_z_x_z = buffer_1010_sppp[209];

    auto g_z_0_y_0_0_z_y_x = buffer_1010_sppp[210];

    auto g_z_0_y_0_0_z_y_y = buffer_1010_sppp[211];

    auto g_z_0_y_0_0_z_y_z = buffer_1010_sppp[212];

    auto g_z_0_y_0_0_z_z_x = buffer_1010_sppp[213];

    auto g_z_0_y_0_0_z_z_y = buffer_1010_sppp[214];

    auto g_z_0_y_0_0_z_z_z = buffer_1010_sppp[215];

    auto g_z_0_z_0_0_x_x_x = buffer_1010_sppp[216];

    auto g_z_0_z_0_0_x_x_y = buffer_1010_sppp[217];

    auto g_z_0_z_0_0_x_x_z = buffer_1010_sppp[218];

    auto g_z_0_z_0_0_x_y_x = buffer_1010_sppp[219];

    auto g_z_0_z_0_0_x_y_y = buffer_1010_sppp[220];

    auto g_z_0_z_0_0_x_y_z = buffer_1010_sppp[221];

    auto g_z_0_z_0_0_x_z_x = buffer_1010_sppp[222];

    auto g_z_0_z_0_0_x_z_y = buffer_1010_sppp[223];

    auto g_z_0_z_0_0_x_z_z = buffer_1010_sppp[224];

    auto g_z_0_z_0_0_y_x_x = buffer_1010_sppp[225];

    auto g_z_0_z_0_0_y_x_y = buffer_1010_sppp[226];

    auto g_z_0_z_0_0_y_x_z = buffer_1010_sppp[227];

    auto g_z_0_z_0_0_y_y_x = buffer_1010_sppp[228];

    auto g_z_0_z_0_0_y_y_y = buffer_1010_sppp[229];

    auto g_z_0_z_0_0_y_y_z = buffer_1010_sppp[230];

    auto g_z_0_z_0_0_y_z_x = buffer_1010_sppp[231];

    auto g_z_0_z_0_0_y_z_y = buffer_1010_sppp[232];

    auto g_z_0_z_0_0_y_z_z = buffer_1010_sppp[233];

    auto g_z_0_z_0_0_z_x_x = buffer_1010_sppp[234];

    auto g_z_0_z_0_0_z_x_y = buffer_1010_sppp[235];

    auto g_z_0_z_0_0_z_x_z = buffer_1010_sppp[236];

    auto g_z_0_z_0_0_z_y_x = buffer_1010_sppp[237];

    auto g_z_0_z_0_0_z_y_y = buffer_1010_sppp[238];

    auto g_z_0_z_0_0_z_y_z = buffer_1010_sppp[239];

    auto g_z_0_z_0_0_z_z_x = buffer_1010_sppp[240];

    auto g_z_0_z_0_0_z_z_y = buffer_1010_sppp[241];

    auto g_z_0_z_0_0_z_z_z = buffer_1010_sppp[242];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_x_0_0_x_x_x, g_x_0_x_0_0_x_x_y, g_x_0_x_0_0_x_x_z, g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_x_x_xx_x, g_x_x_xx_y, g_x_x_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_x_x_x[i] = -2.0 * g_x_x_0_x[i] * a_exp + 4.0 * g_x_x_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_x_y[i] = -2.0 * g_x_x_0_y[i] * a_exp + 4.0 * g_x_x_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_x_z[i] = -2.0 * g_x_x_0_z[i] * a_exp + 4.0 * g_x_x_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_x_0_0_x_y_x, g_x_0_x_0_0_x_y_y, g_x_0_x_0_0_x_y_z, g_x_x_xy_x, g_x_x_xy_y, g_x_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_x_y_x[i] = 4.0 * g_x_x_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_y_y[i] = 4.0 * g_x_x_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_y_z[i] = 4.0 * g_x_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_x_0_0_x_z_x, g_x_0_x_0_0_x_z_y, g_x_0_x_0_0_x_z_z, g_x_x_xz_x, g_x_x_xz_y, g_x_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_x_z_x[i] = 4.0 * g_x_x_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_z_y[i] = 4.0 * g_x_x_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_z_z[i] = 4.0 * g_x_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_0_x_0_0_y_x_x, g_x_0_x_0_0_y_x_y, g_x_0_x_0_0_y_x_z, g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_x_y_xx_x, g_x_y_xx_y, g_x_y_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_y_x_x[i] = -2.0 * g_x_y_0_x[i] * a_exp + 4.0 * g_x_y_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_x_y[i] = -2.0 * g_x_y_0_y[i] * a_exp + 4.0 * g_x_y_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_x_z[i] = -2.0 * g_x_y_0_z[i] * a_exp + 4.0 * g_x_y_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_0_x_0_0_y_y_x, g_x_0_x_0_0_y_y_y, g_x_0_x_0_0_y_y_z, g_x_y_xy_x, g_x_y_xy_y, g_x_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_y_y_x[i] = 4.0 * g_x_y_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_y_y[i] = 4.0 * g_x_y_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_y_z[i] = 4.0 * g_x_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_0_x_0_0_y_z_x, g_x_0_x_0_0_y_z_y, g_x_0_x_0_0_y_z_z, g_x_y_xz_x, g_x_y_xz_y, g_x_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_y_z_x[i] = 4.0 * g_x_y_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_z_y[i] = 4.0 * g_x_y_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_z_z[i] = 4.0 * g_x_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_0_x_0_0_z_x_x, g_x_0_x_0_0_z_x_y, g_x_0_x_0_0_z_x_z, g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_x_z_xx_x, g_x_z_xx_y, g_x_z_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_z_x_x[i] = -2.0 * g_x_z_0_x[i] * a_exp + 4.0 * g_x_z_xx_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_x_y[i] = -2.0 * g_x_z_0_y[i] * a_exp + 4.0 * g_x_z_xx_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_x_z[i] = -2.0 * g_x_z_0_z[i] * a_exp + 4.0 * g_x_z_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_0_x_0_0_z_y_x, g_x_0_x_0_0_z_y_y, g_x_0_x_0_0_z_y_z, g_x_z_xy_x, g_x_z_xy_y, g_x_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_z_y_x[i] = 4.0 * g_x_z_xy_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_y_y[i] = 4.0 * g_x_z_xy_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_y_z[i] = 4.0 * g_x_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_0_x_0_0_z_z_x, g_x_0_x_0_0_z_z_y, g_x_0_x_0_0_z_z_z, g_x_z_xz_x, g_x_z_xz_y, g_x_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_z_z_x[i] = 4.0 * g_x_z_xz_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_z_y[i] = 4.0 * g_x_z_xz_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_z_z[i] = 4.0 * g_x_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_0_y_0_0_x_x_x, g_x_0_y_0_0_x_x_y, g_x_0_y_0_0_x_x_z, g_x_x_xy_x, g_x_x_xy_y, g_x_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_x_x_x[i] = 4.0 * g_x_x_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_x_y[i] = 4.0 * g_x_x_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_x_z[i] = 4.0 * g_x_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_x_0_y_0_0_x_y_x, g_x_0_y_0_0_x_y_y, g_x_0_y_0_0_x_y_z, g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_x_x_yy_x, g_x_x_yy_y, g_x_x_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_x_y_x[i] = -2.0 * g_x_x_0_x[i] * a_exp + 4.0 * g_x_x_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_y_y[i] = -2.0 * g_x_x_0_y[i] * a_exp + 4.0 * g_x_x_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_y_z[i] = -2.0 * g_x_x_0_z[i] * a_exp + 4.0 * g_x_x_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_x_0_y_0_0_x_z_x, g_x_0_y_0_0_x_z_y, g_x_0_y_0_0_x_z_z, g_x_x_yz_x, g_x_x_yz_y, g_x_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_x_z_x[i] = 4.0 * g_x_x_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_z_y[i] = 4.0 * g_x_x_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_z_z[i] = 4.0 * g_x_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_0_y_0_0_y_x_x, g_x_0_y_0_0_y_x_y, g_x_0_y_0_0_y_x_z, g_x_y_xy_x, g_x_y_xy_y, g_x_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_y_x_x[i] = 4.0 * g_x_y_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_x_y[i] = 4.0 * g_x_y_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_x_z[i] = 4.0 * g_x_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_0_y_0_0_y_y_x, g_x_0_y_0_0_y_y_y, g_x_0_y_0_0_y_y_z, g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_x_y_yy_x, g_x_y_yy_y, g_x_y_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_y_y_x[i] = -2.0 * g_x_y_0_x[i] * a_exp + 4.0 * g_x_y_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_y_y[i] = -2.0 * g_x_y_0_y[i] * a_exp + 4.0 * g_x_y_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_y_z[i] = -2.0 * g_x_y_0_z[i] * a_exp + 4.0 * g_x_y_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_0_y_0_0_y_z_x, g_x_0_y_0_0_y_z_y, g_x_0_y_0_0_y_z_z, g_x_y_yz_x, g_x_y_yz_y, g_x_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_y_z_x[i] = 4.0 * g_x_y_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_z_y[i] = 4.0 * g_x_y_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_z_z[i] = 4.0 * g_x_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_0_y_0_0_z_x_x, g_x_0_y_0_0_z_x_y, g_x_0_y_0_0_z_x_z, g_x_z_xy_x, g_x_z_xy_y, g_x_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_z_x_x[i] = 4.0 * g_x_z_xy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_x_y[i] = 4.0 * g_x_z_xy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_x_z[i] = 4.0 * g_x_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_0_y_0_0_z_y_x, g_x_0_y_0_0_z_y_y, g_x_0_y_0_0_z_y_z, g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_x_z_yy_x, g_x_z_yy_y, g_x_z_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_z_y_x[i] = -2.0 * g_x_z_0_x[i] * a_exp + 4.0 * g_x_z_yy_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_y_y[i] = -2.0 * g_x_z_0_y[i] * a_exp + 4.0 * g_x_z_yy_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_y_z[i] = -2.0 * g_x_z_0_z[i] * a_exp + 4.0 * g_x_z_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_0_y_0_0_z_z_x, g_x_0_y_0_0_z_z_y, g_x_0_y_0_0_z_z_z, g_x_z_yz_x, g_x_z_yz_y, g_x_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_z_z_x[i] = 4.0 * g_x_z_yz_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_z_y[i] = 4.0 * g_x_z_yz_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_z_z[i] = 4.0 * g_x_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_x_0_z_0_0_x_x_x, g_x_0_z_0_0_x_x_y, g_x_0_z_0_0_x_x_z, g_x_x_xz_x, g_x_x_xz_y, g_x_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_x_x_x[i] = 4.0 * g_x_x_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_x_y[i] = 4.0 * g_x_x_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_x_z[i] = 4.0 * g_x_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_x_0_z_0_0_x_y_x, g_x_0_z_0_0_x_y_y, g_x_0_z_0_0_x_y_z, g_x_x_yz_x, g_x_x_yz_y, g_x_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_x_y_x[i] = 4.0 * g_x_x_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_y_y[i] = 4.0 * g_x_x_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_y_z[i] = 4.0 * g_x_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_x_0_z_0_0_x_z_x, g_x_0_z_0_0_x_z_y, g_x_0_z_0_0_x_z_z, g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_x_x_zz_x, g_x_x_zz_y, g_x_x_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_x_z_x[i] = -2.0 * g_x_x_0_x[i] * a_exp + 4.0 * g_x_x_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_z_y[i] = -2.0 * g_x_x_0_y[i] * a_exp + 4.0 * g_x_x_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_z_z[i] = -2.0 * g_x_x_0_z[i] * a_exp + 4.0 * g_x_x_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_x_0_z_0_0_y_x_x, g_x_0_z_0_0_y_x_y, g_x_0_z_0_0_y_x_z, g_x_y_xz_x, g_x_y_xz_y, g_x_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_y_x_x[i] = 4.0 * g_x_y_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_x_y[i] = 4.0 * g_x_y_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_x_z[i] = 4.0 * g_x_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_x_0_z_0_0_y_y_x, g_x_0_z_0_0_y_y_y, g_x_0_z_0_0_y_y_z, g_x_y_yz_x, g_x_y_yz_y, g_x_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_y_y_x[i] = 4.0 * g_x_y_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_y_y[i] = 4.0 * g_x_y_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_y_z[i] = 4.0 * g_x_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_x_0_z_0_0_y_z_x, g_x_0_z_0_0_y_z_y, g_x_0_z_0_0_y_z_z, g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_x_y_zz_x, g_x_y_zz_y, g_x_y_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_y_z_x[i] = -2.0 * g_x_y_0_x[i] * a_exp + 4.0 * g_x_y_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_z_y[i] = -2.0 * g_x_y_0_y[i] * a_exp + 4.0 * g_x_y_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_z_z[i] = -2.0 * g_x_y_0_z[i] * a_exp + 4.0 * g_x_y_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_0_z_0_0_z_x_x, g_x_0_z_0_0_z_x_y, g_x_0_z_0_0_z_x_z, g_x_z_xz_x, g_x_z_xz_y, g_x_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_z_x_x[i] = 4.0 * g_x_z_xz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_x_y[i] = 4.0 * g_x_z_xz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_x_z[i] = 4.0 * g_x_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_0_z_0_0_z_y_x, g_x_0_z_0_0_z_y_y, g_x_0_z_0_0_z_y_z, g_x_z_yz_x, g_x_z_yz_y, g_x_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_z_y_x[i] = 4.0 * g_x_z_yz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_y_y[i] = 4.0 * g_x_z_yz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_y_z[i] = 4.0 * g_x_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_0_z_0_0_z_z_x, g_x_0_z_0_0_z_z_y, g_x_0_z_0_0_z_z_z, g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_x_z_zz_x, g_x_z_zz_y, g_x_z_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_z_z_x[i] = -2.0 * g_x_z_0_x[i] * a_exp + 4.0 * g_x_z_zz_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_z_y[i] = -2.0 * g_x_z_0_y[i] * a_exp + 4.0 * g_x_z_zz_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_z_z[i] = -2.0 * g_x_z_0_z[i] * a_exp + 4.0 * g_x_z_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_y_0_x_0_0_x_x_x, g_y_0_x_0_0_x_x_y, g_y_0_x_0_0_x_x_z, g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_y_x_xx_x, g_y_x_xx_y, g_y_x_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_x_x_x[i] = -2.0 * g_y_x_0_x[i] * a_exp + 4.0 * g_y_x_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_x_y[i] = -2.0 * g_y_x_0_y[i] * a_exp + 4.0 * g_y_x_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_x_z[i] = -2.0 * g_y_x_0_z[i] * a_exp + 4.0 * g_y_x_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_y_0_x_0_0_x_y_x, g_y_0_x_0_0_x_y_y, g_y_0_x_0_0_x_y_z, g_y_x_xy_x, g_y_x_xy_y, g_y_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_x_y_x[i] = 4.0 * g_y_x_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_y_y[i] = 4.0 * g_y_x_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_y_z[i] = 4.0 * g_y_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_y_0_x_0_0_x_z_x, g_y_0_x_0_0_x_z_y, g_y_0_x_0_0_x_z_z, g_y_x_xz_x, g_y_x_xz_y, g_y_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_x_z_x[i] = 4.0 * g_y_x_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_z_y[i] = 4.0 * g_y_x_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_z_z[i] = 4.0 * g_y_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_y_0_x_0_0_y_x_x, g_y_0_x_0_0_y_x_y, g_y_0_x_0_0_y_x_z, g_y_y_0_x, g_y_y_0_y, g_y_y_0_z, g_y_y_xx_x, g_y_y_xx_y, g_y_y_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_y_x_x[i] = -2.0 * g_y_y_0_x[i] * a_exp + 4.0 * g_y_y_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_x_y[i] = -2.0 * g_y_y_0_y[i] * a_exp + 4.0 * g_y_y_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_x_z[i] = -2.0 * g_y_y_0_z[i] * a_exp + 4.0 * g_y_y_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_y_0_x_0_0_y_y_x, g_y_0_x_0_0_y_y_y, g_y_0_x_0_0_y_y_z, g_y_y_xy_x, g_y_y_xy_y, g_y_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_y_y_x[i] = 4.0 * g_y_y_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_y_y[i] = 4.0 * g_y_y_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_y_z[i] = 4.0 * g_y_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_y_0_x_0_0_y_z_x, g_y_0_x_0_0_y_z_y, g_y_0_x_0_0_y_z_z, g_y_y_xz_x, g_y_y_xz_y, g_y_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_y_z_x[i] = 4.0 * g_y_y_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_z_y[i] = 4.0 * g_y_y_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_z_z[i] = 4.0 * g_y_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_y_0_x_0_0_z_x_x, g_y_0_x_0_0_z_x_y, g_y_0_x_0_0_z_x_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z, g_y_z_xx_x, g_y_z_xx_y, g_y_z_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_z_x_x[i] = -2.0 * g_y_z_0_x[i] * a_exp + 4.0 * g_y_z_xx_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_x_y[i] = -2.0 * g_y_z_0_y[i] * a_exp + 4.0 * g_y_z_xx_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_x_z[i] = -2.0 * g_y_z_0_z[i] * a_exp + 4.0 * g_y_z_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_y_0_x_0_0_z_y_x, g_y_0_x_0_0_z_y_y, g_y_0_x_0_0_z_y_z, g_y_z_xy_x, g_y_z_xy_y, g_y_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_z_y_x[i] = 4.0 * g_y_z_xy_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_y_y[i] = 4.0 * g_y_z_xy_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_y_z[i] = 4.0 * g_y_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_y_0_x_0_0_z_z_x, g_y_0_x_0_0_z_z_y, g_y_0_x_0_0_z_z_z, g_y_z_xz_x, g_y_z_xz_y, g_y_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_z_z_x[i] = 4.0 * g_y_z_xz_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_z_y[i] = 4.0 * g_y_z_xz_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_z_z[i] = 4.0 * g_y_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_y_0_y_0_0_x_x_x, g_y_0_y_0_0_x_x_y, g_y_0_y_0_0_x_x_z, g_y_x_xy_x, g_y_x_xy_y, g_y_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_x_x_x[i] = 4.0 * g_y_x_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_x_y[i] = 4.0 * g_y_x_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_x_z[i] = 4.0 * g_y_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_y_0_y_0_0_x_y_x, g_y_0_y_0_0_x_y_y, g_y_0_y_0_0_x_y_z, g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_y_x_yy_x, g_y_x_yy_y, g_y_x_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_x_y_x[i] = -2.0 * g_y_x_0_x[i] * a_exp + 4.0 * g_y_x_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_y_y[i] = -2.0 * g_y_x_0_y[i] * a_exp + 4.0 * g_y_x_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_y_z[i] = -2.0 * g_y_x_0_z[i] * a_exp + 4.0 * g_y_x_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_y_0_y_0_0_x_z_x, g_y_0_y_0_0_x_z_y, g_y_0_y_0_0_x_z_z, g_y_x_yz_x, g_y_x_yz_y, g_y_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_x_z_x[i] = 4.0 * g_y_x_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_z_y[i] = 4.0 * g_y_x_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_z_z[i] = 4.0 * g_y_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_y_0_y_0_0_y_x_x, g_y_0_y_0_0_y_x_y, g_y_0_y_0_0_y_x_z, g_y_y_xy_x, g_y_y_xy_y, g_y_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_y_x_x[i] = 4.0 * g_y_y_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_x_y[i] = 4.0 * g_y_y_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_x_z[i] = 4.0 * g_y_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_y_0_y_0_0_y_y_x, g_y_0_y_0_0_y_y_y, g_y_0_y_0_0_y_y_z, g_y_y_0_x, g_y_y_0_y, g_y_y_0_z, g_y_y_yy_x, g_y_y_yy_y, g_y_y_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_y_y_x[i] = -2.0 * g_y_y_0_x[i] * a_exp + 4.0 * g_y_y_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_y_y[i] = -2.0 * g_y_y_0_y[i] * a_exp + 4.0 * g_y_y_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_y_z[i] = -2.0 * g_y_y_0_z[i] * a_exp + 4.0 * g_y_y_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_y_0_y_0_0_y_z_x, g_y_0_y_0_0_y_z_y, g_y_0_y_0_0_y_z_z, g_y_y_yz_x, g_y_y_yz_y, g_y_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_y_z_x[i] = 4.0 * g_y_y_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_z_y[i] = 4.0 * g_y_y_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_z_z[i] = 4.0 * g_y_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_y_0_y_0_0_z_x_x, g_y_0_y_0_0_z_x_y, g_y_0_y_0_0_z_x_z, g_y_z_xy_x, g_y_z_xy_y, g_y_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_z_x_x[i] = 4.0 * g_y_z_xy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_x_y[i] = 4.0 * g_y_z_xy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_x_z[i] = 4.0 * g_y_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_y_0_y_0_0_z_y_x, g_y_0_y_0_0_z_y_y, g_y_0_y_0_0_z_y_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z, g_y_z_yy_x, g_y_z_yy_y, g_y_z_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_z_y_x[i] = -2.0 * g_y_z_0_x[i] * a_exp + 4.0 * g_y_z_yy_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_y_y[i] = -2.0 * g_y_z_0_y[i] * a_exp + 4.0 * g_y_z_yy_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_y_z[i] = -2.0 * g_y_z_0_z[i] * a_exp + 4.0 * g_y_z_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_y_0_y_0_0_z_z_x, g_y_0_y_0_0_z_z_y, g_y_0_y_0_0_z_z_z, g_y_z_yz_x, g_y_z_yz_y, g_y_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_z_z_x[i] = 4.0 * g_y_z_yz_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_z_y[i] = 4.0 * g_y_z_yz_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_z_z[i] = 4.0 * g_y_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_y_0_z_0_0_x_x_x, g_y_0_z_0_0_x_x_y, g_y_0_z_0_0_x_x_z, g_y_x_xz_x, g_y_x_xz_y, g_y_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_x_x_x[i] = 4.0 * g_y_x_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_x_y[i] = 4.0 * g_y_x_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_x_z[i] = 4.0 * g_y_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_y_0_z_0_0_x_y_x, g_y_0_z_0_0_x_y_y, g_y_0_z_0_0_x_y_z, g_y_x_yz_x, g_y_x_yz_y, g_y_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_x_y_x[i] = 4.0 * g_y_x_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_y_y[i] = 4.0 * g_y_x_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_y_z[i] = 4.0 * g_y_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_y_0_z_0_0_x_z_x, g_y_0_z_0_0_x_z_y, g_y_0_z_0_0_x_z_z, g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_y_x_zz_x, g_y_x_zz_y, g_y_x_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_x_z_x[i] = -2.0 * g_y_x_0_x[i] * a_exp + 4.0 * g_y_x_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_z_y[i] = -2.0 * g_y_x_0_y[i] * a_exp + 4.0 * g_y_x_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_z_z[i] = -2.0 * g_y_x_0_z[i] * a_exp + 4.0 * g_y_x_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_y_0_z_0_0_y_x_x, g_y_0_z_0_0_y_x_y, g_y_0_z_0_0_y_x_z, g_y_y_xz_x, g_y_y_xz_y, g_y_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_y_x_x[i] = 4.0 * g_y_y_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_x_y[i] = 4.0 * g_y_y_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_x_z[i] = 4.0 * g_y_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_y_0_z_0_0_y_y_x, g_y_0_z_0_0_y_y_y, g_y_0_z_0_0_y_y_z, g_y_y_yz_x, g_y_y_yz_y, g_y_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_y_y_x[i] = 4.0 * g_y_y_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_y_y[i] = 4.0 * g_y_y_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_y_z[i] = 4.0 * g_y_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_y_0_z_0_0_y_z_x, g_y_0_z_0_0_y_z_y, g_y_0_z_0_0_y_z_z, g_y_y_0_x, g_y_y_0_y, g_y_y_0_z, g_y_y_zz_x, g_y_y_zz_y, g_y_y_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_y_z_x[i] = -2.0 * g_y_y_0_x[i] * a_exp + 4.0 * g_y_y_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_z_y[i] = -2.0 * g_y_y_0_y[i] * a_exp + 4.0 * g_y_y_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_z_z[i] = -2.0 * g_y_y_0_z[i] * a_exp + 4.0 * g_y_y_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_y_0_z_0_0_z_x_x, g_y_0_z_0_0_z_x_y, g_y_0_z_0_0_z_x_z, g_y_z_xz_x, g_y_z_xz_y, g_y_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_z_x_x[i] = 4.0 * g_y_z_xz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_x_y[i] = 4.0 * g_y_z_xz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_x_z[i] = 4.0 * g_y_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_y_0_z_0_0_z_y_x, g_y_0_z_0_0_z_y_y, g_y_0_z_0_0_z_y_z, g_y_z_yz_x, g_y_z_yz_y, g_y_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_z_y_x[i] = 4.0 * g_y_z_yz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_y_y[i] = 4.0 * g_y_z_yz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_y_z[i] = 4.0 * g_y_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_y_0_z_0_0_z_z_x, g_y_0_z_0_0_z_z_y, g_y_0_z_0_0_z_z_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z, g_y_z_zz_x, g_y_z_zz_y, g_y_z_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_z_z_x[i] = -2.0 * g_y_z_0_x[i] * a_exp + 4.0 * g_y_z_zz_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_z_y[i] = -2.0 * g_y_z_0_y[i] * a_exp + 4.0 * g_y_z_zz_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_z_z[i] = -2.0 * g_y_z_0_z[i] * a_exp + 4.0 * g_y_z_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_z_0_x_0_0_x_x_x, g_z_0_x_0_0_x_x_y, g_z_0_x_0_0_x_x_z, g_z_x_0_x, g_z_x_0_y, g_z_x_0_z, g_z_x_xx_x, g_z_x_xx_y, g_z_x_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_x_x_x[i] = -2.0 * g_z_x_0_x[i] * a_exp + 4.0 * g_z_x_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_x_y[i] = -2.0 * g_z_x_0_y[i] * a_exp + 4.0 * g_z_x_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_x_z[i] = -2.0 * g_z_x_0_z[i] * a_exp + 4.0 * g_z_x_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_z_0_x_0_0_x_y_x, g_z_0_x_0_0_x_y_y, g_z_0_x_0_0_x_y_z, g_z_x_xy_x, g_z_x_xy_y, g_z_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_x_y_x[i] = 4.0 * g_z_x_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_y_y[i] = 4.0 * g_z_x_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_y_z[i] = 4.0 * g_z_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_z_0_x_0_0_x_z_x, g_z_0_x_0_0_x_z_y, g_z_0_x_0_0_x_z_z, g_z_x_xz_x, g_z_x_xz_y, g_z_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_x_z_x[i] = 4.0 * g_z_x_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_z_y[i] = 4.0 * g_z_x_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_z_z[i] = 4.0 * g_z_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_z_0_x_0_0_y_x_x, g_z_0_x_0_0_y_x_y, g_z_0_x_0_0_y_x_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z, g_z_y_xx_x, g_z_y_xx_y, g_z_y_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_y_x_x[i] = -2.0 * g_z_y_0_x[i] * a_exp + 4.0 * g_z_y_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_x_y[i] = -2.0 * g_z_y_0_y[i] * a_exp + 4.0 * g_z_y_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_x_z[i] = -2.0 * g_z_y_0_z[i] * a_exp + 4.0 * g_z_y_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_z_0_x_0_0_y_y_x, g_z_0_x_0_0_y_y_y, g_z_0_x_0_0_y_y_z, g_z_y_xy_x, g_z_y_xy_y, g_z_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_y_y_x[i] = 4.0 * g_z_y_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_y_y[i] = 4.0 * g_z_y_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_y_z[i] = 4.0 * g_z_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_z_0_x_0_0_y_z_x, g_z_0_x_0_0_y_z_y, g_z_0_x_0_0_y_z_z, g_z_y_xz_x, g_z_y_xz_y, g_z_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_y_z_x[i] = 4.0 * g_z_y_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_z_y[i] = 4.0 * g_z_y_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_z_z[i] = 4.0 * g_z_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_z_0_x_0_0_z_x_x, g_z_0_x_0_0_z_x_y, g_z_0_x_0_0_z_x_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z, g_z_z_xx_x, g_z_z_xx_y, g_z_z_xx_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_z_x_x[i] = -2.0 * g_z_z_0_x[i] * a_exp + 4.0 * g_z_z_xx_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_x_y[i] = -2.0 * g_z_z_0_y[i] * a_exp + 4.0 * g_z_z_xx_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_x_z[i] = -2.0 * g_z_z_0_z[i] * a_exp + 4.0 * g_z_z_xx_z[i] * a_exp * c_exps[i];
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_z_0_x_0_0_z_y_x, g_z_0_x_0_0_z_y_y, g_z_0_x_0_0_z_y_z, g_z_z_xy_x, g_z_z_xy_y, g_z_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_z_y_x[i] = 4.0 * g_z_z_xy_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_y_y[i] = 4.0 * g_z_z_xy_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_y_z[i] = 4.0 * g_z_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_z_0_x_0_0_z_z_x, g_z_0_x_0_0_z_z_y, g_z_0_x_0_0_z_z_z, g_z_z_xz_x, g_z_z_xz_y, g_z_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_z_z_x[i] = 4.0 * g_z_z_xz_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_z_y[i] = 4.0 * g_z_z_xz_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_z_z[i] = 4.0 * g_z_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_z_0_y_0_0_x_x_x, g_z_0_y_0_0_x_x_y, g_z_0_y_0_0_x_x_z, g_z_x_xy_x, g_z_x_xy_y, g_z_x_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_x_x_x[i] = 4.0 * g_z_x_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_x_y[i] = 4.0 * g_z_x_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_x_z[i] = 4.0 * g_z_x_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_z_0_y_0_0_x_y_x, g_z_0_y_0_0_x_y_y, g_z_0_y_0_0_x_y_z, g_z_x_0_x, g_z_x_0_y, g_z_x_0_z, g_z_x_yy_x, g_z_x_yy_y, g_z_x_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_x_y_x[i] = -2.0 * g_z_x_0_x[i] * a_exp + 4.0 * g_z_x_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_y_y[i] = -2.0 * g_z_x_0_y[i] * a_exp + 4.0 * g_z_x_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_y_z[i] = -2.0 * g_z_x_0_z[i] * a_exp + 4.0 * g_z_x_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_z_0_y_0_0_x_z_x, g_z_0_y_0_0_x_z_y, g_z_0_y_0_0_x_z_z, g_z_x_yz_x, g_z_x_yz_y, g_z_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_x_z_x[i] = 4.0 * g_z_x_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_z_y[i] = 4.0 * g_z_x_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_z_z[i] = 4.0 * g_z_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_z_0_y_0_0_y_x_x, g_z_0_y_0_0_y_x_y, g_z_0_y_0_0_y_x_z, g_z_y_xy_x, g_z_y_xy_y, g_z_y_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_y_x_x[i] = 4.0 * g_z_y_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_x_y[i] = 4.0 * g_z_y_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_x_z[i] = 4.0 * g_z_y_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_z_0_y_0_0_y_y_x, g_z_0_y_0_0_y_y_y, g_z_0_y_0_0_y_y_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z, g_z_y_yy_x, g_z_y_yy_y, g_z_y_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_y_y_x[i] = -2.0 * g_z_y_0_x[i] * a_exp + 4.0 * g_z_y_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_y_y[i] = -2.0 * g_z_y_0_y[i] * a_exp + 4.0 * g_z_y_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_y_z[i] = -2.0 * g_z_y_0_z[i] * a_exp + 4.0 * g_z_y_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_z_0_y_0_0_y_z_x, g_z_0_y_0_0_y_z_y, g_z_0_y_0_0_y_z_z, g_z_y_yz_x, g_z_y_yz_y, g_z_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_y_z_x[i] = 4.0 * g_z_y_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_z_y[i] = 4.0 * g_z_y_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_z_z[i] = 4.0 * g_z_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_z_0_y_0_0_z_x_x, g_z_0_y_0_0_z_x_y, g_z_0_y_0_0_z_x_z, g_z_z_xy_x, g_z_z_xy_y, g_z_z_xy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_z_x_x[i] = 4.0 * g_z_z_xy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_x_y[i] = 4.0 * g_z_z_xy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_x_z[i] = 4.0 * g_z_z_xy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_z_0_y_0_0_z_y_x, g_z_0_y_0_0_z_y_y, g_z_0_y_0_0_z_y_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z, g_z_z_yy_x, g_z_z_yy_y, g_z_z_yy_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_z_y_x[i] = -2.0 * g_z_z_0_x[i] * a_exp + 4.0 * g_z_z_yy_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_y_y[i] = -2.0 * g_z_z_0_y[i] * a_exp + 4.0 * g_z_z_yy_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_y_z[i] = -2.0 * g_z_z_0_z[i] * a_exp + 4.0 * g_z_z_yy_z[i] * a_exp * c_exps[i];
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_z_0_y_0_0_z_z_x, g_z_0_y_0_0_z_z_y, g_z_0_y_0_0_z_z_z, g_z_z_yz_x, g_z_z_yz_y, g_z_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_z_z_x[i] = 4.0 * g_z_z_yz_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_z_y[i] = 4.0 * g_z_z_yz_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_z_z[i] = 4.0 * g_z_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_z_0_z_0_0_x_x_x, g_z_0_z_0_0_x_x_y, g_z_0_z_0_0_x_x_z, g_z_x_xz_x, g_z_x_xz_y, g_z_x_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_x_x_x[i] = 4.0 * g_z_x_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_x_y[i] = 4.0 * g_z_x_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_x_z[i] = 4.0 * g_z_x_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_z_0_z_0_0_x_y_x, g_z_0_z_0_0_x_y_y, g_z_0_z_0_0_x_y_z, g_z_x_yz_x, g_z_x_yz_y, g_z_x_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_x_y_x[i] = 4.0 * g_z_x_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_y_y[i] = 4.0 * g_z_x_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_y_z[i] = 4.0 * g_z_x_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_z_0_z_0_0_x_z_x, g_z_0_z_0_0_x_z_y, g_z_0_z_0_0_x_z_z, g_z_x_0_x, g_z_x_0_y, g_z_x_0_z, g_z_x_zz_x, g_z_x_zz_y, g_z_x_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_x_z_x[i] = -2.0 * g_z_x_0_x[i] * a_exp + 4.0 * g_z_x_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_z_y[i] = -2.0 * g_z_x_0_y[i] * a_exp + 4.0 * g_z_x_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_z_z[i] = -2.0 * g_z_x_0_z[i] * a_exp + 4.0 * g_z_x_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_z_0_z_0_0_y_x_x, g_z_0_z_0_0_y_x_y, g_z_0_z_0_0_y_x_z, g_z_y_xz_x, g_z_y_xz_y, g_z_y_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_y_x_x[i] = 4.0 * g_z_y_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_x_y[i] = 4.0 * g_z_y_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_x_z[i] = 4.0 * g_z_y_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_z_0_z_0_0_y_y_x, g_z_0_z_0_0_y_y_y, g_z_0_z_0_0_y_y_z, g_z_y_yz_x, g_z_y_yz_y, g_z_y_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_y_y_x[i] = 4.0 * g_z_y_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_y_y[i] = 4.0 * g_z_y_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_y_z[i] = 4.0 * g_z_y_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_z_0_z_0_0_y_z_x, g_z_0_z_0_0_y_z_y, g_z_0_z_0_0_y_z_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z, g_z_y_zz_x, g_z_y_zz_y, g_z_y_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_y_z_x[i] = -2.0 * g_z_y_0_x[i] * a_exp + 4.0 * g_z_y_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_z_y[i] = -2.0 * g_z_y_0_y[i] * a_exp + 4.0 * g_z_y_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_z_z[i] = -2.0 * g_z_y_0_z[i] * a_exp + 4.0 * g_z_y_zz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_z_0_z_0_0_z_x_x, g_z_0_z_0_0_z_x_y, g_z_0_z_0_0_z_x_z, g_z_z_xz_x, g_z_z_xz_y, g_z_z_xz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_z_x_x[i] = 4.0 * g_z_z_xz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_x_y[i] = 4.0 * g_z_z_xz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_x_z[i] = 4.0 * g_z_z_xz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_z_0_z_0_0_z_y_x, g_z_0_z_0_0_z_y_y, g_z_0_z_0_0_z_y_z, g_z_z_yz_x, g_z_z_yz_y, g_z_z_yz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_z_y_x[i] = 4.0 * g_z_z_yz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_y_y[i] = 4.0 * g_z_z_yz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_y_z[i] = 4.0 * g_z_z_yz_z[i] * a_exp * c_exps[i];
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_z_0_z_0_0_z_z_x, g_z_0_z_0_0_z_z_y, g_z_0_z_0_0_z_z_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z, g_z_z_zz_x, g_z_z_zz_y, g_z_z_zz_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_z_z_x[i] = -2.0 * g_z_z_0_x[i] * a_exp + 4.0 * g_z_z_zz_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_z_y[i] = -2.0 * g_z_z_0_y[i] * a_exp + 4.0 * g_z_z_zz_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_z_z[i] = -2.0 * g_z_z_0_z[i] * a_exp + 4.0 * g_z_z_zz_z[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace

