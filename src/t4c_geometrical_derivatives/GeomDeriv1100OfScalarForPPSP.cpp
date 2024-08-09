#include "GeomDeriv1100OfScalarForPPSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_ppsp_0(CSimdArray<double>& buffer_1100_ppsp,
                     const CSimdArray<double>& buffer_sssp,
                     const CSimdArray<double>& buffer_sdsp,
                     const CSimdArray<double>& buffer_dssp,
                     const CSimdArray<double>& buffer_ddsp,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_ppsp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sssp

    auto g_0_0_0_x = buffer_sssp[0];

    auto g_0_0_0_y = buffer_sssp[1];

    auto g_0_0_0_z = buffer_sssp[2];

    /// Set up components of auxilary buffer : buffer_sdsp

    auto g_0_xx_0_x = buffer_sdsp[0];

    auto g_0_xx_0_y = buffer_sdsp[1];

    auto g_0_xx_0_z = buffer_sdsp[2];

    auto g_0_xy_0_x = buffer_sdsp[3];

    auto g_0_xy_0_y = buffer_sdsp[4];

    auto g_0_xy_0_z = buffer_sdsp[5];

    auto g_0_xz_0_x = buffer_sdsp[6];

    auto g_0_xz_0_y = buffer_sdsp[7];

    auto g_0_xz_0_z = buffer_sdsp[8];

    auto g_0_yy_0_x = buffer_sdsp[9];

    auto g_0_yy_0_y = buffer_sdsp[10];

    auto g_0_yy_0_z = buffer_sdsp[11];

    auto g_0_yz_0_x = buffer_sdsp[12];

    auto g_0_yz_0_y = buffer_sdsp[13];

    auto g_0_yz_0_z = buffer_sdsp[14];

    auto g_0_zz_0_x = buffer_sdsp[15];

    auto g_0_zz_0_y = buffer_sdsp[16];

    auto g_0_zz_0_z = buffer_sdsp[17];

    /// Set up components of auxilary buffer : buffer_dssp

    auto g_xx_0_0_x = buffer_dssp[0];

    auto g_xx_0_0_y = buffer_dssp[1];

    auto g_xx_0_0_z = buffer_dssp[2];

    auto g_xy_0_0_x = buffer_dssp[3];

    auto g_xy_0_0_y = buffer_dssp[4];

    auto g_xy_0_0_z = buffer_dssp[5];

    auto g_xz_0_0_x = buffer_dssp[6];

    auto g_xz_0_0_y = buffer_dssp[7];

    auto g_xz_0_0_z = buffer_dssp[8];

    auto g_yy_0_0_x = buffer_dssp[9];

    auto g_yy_0_0_y = buffer_dssp[10];

    auto g_yy_0_0_z = buffer_dssp[11];

    auto g_yz_0_0_x = buffer_dssp[12];

    auto g_yz_0_0_y = buffer_dssp[13];

    auto g_yz_0_0_z = buffer_dssp[14];

    auto g_zz_0_0_x = buffer_dssp[15];

    auto g_zz_0_0_y = buffer_dssp[16];

    auto g_zz_0_0_z = buffer_dssp[17];

    /// Set up components of auxilary buffer : buffer_ddsp

    auto g_xx_xx_0_x = buffer_ddsp[0];

    auto g_xx_xx_0_y = buffer_ddsp[1];

    auto g_xx_xx_0_z = buffer_ddsp[2];

    auto g_xx_xy_0_x = buffer_ddsp[3];

    auto g_xx_xy_0_y = buffer_ddsp[4];

    auto g_xx_xy_0_z = buffer_ddsp[5];

    auto g_xx_xz_0_x = buffer_ddsp[6];

    auto g_xx_xz_0_y = buffer_ddsp[7];

    auto g_xx_xz_0_z = buffer_ddsp[8];

    auto g_xx_yy_0_x = buffer_ddsp[9];

    auto g_xx_yy_0_y = buffer_ddsp[10];

    auto g_xx_yy_0_z = buffer_ddsp[11];

    auto g_xx_yz_0_x = buffer_ddsp[12];

    auto g_xx_yz_0_y = buffer_ddsp[13];

    auto g_xx_yz_0_z = buffer_ddsp[14];

    auto g_xx_zz_0_x = buffer_ddsp[15];

    auto g_xx_zz_0_y = buffer_ddsp[16];

    auto g_xx_zz_0_z = buffer_ddsp[17];

    auto g_xy_xx_0_x = buffer_ddsp[18];

    auto g_xy_xx_0_y = buffer_ddsp[19];

    auto g_xy_xx_0_z = buffer_ddsp[20];

    auto g_xy_xy_0_x = buffer_ddsp[21];

    auto g_xy_xy_0_y = buffer_ddsp[22];

    auto g_xy_xy_0_z = buffer_ddsp[23];

    auto g_xy_xz_0_x = buffer_ddsp[24];

    auto g_xy_xz_0_y = buffer_ddsp[25];

    auto g_xy_xz_0_z = buffer_ddsp[26];

    auto g_xy_yy_0_x = buffer_ddsp[27];

    auto g_xy_yy_0_y = buffer_ddsp[28];

    auto g_xy_yy_0_z = buffer_ddsp[29];

    auto g_xy_yz_0_x = buffer_ddsp[30];

    auto g_xy_yz_0_y = buffer_ddsp[31];

    auto g_xy_yz_0_z = buffer_ddsp[32];

    auto g_xy_zz_0_x = buffer_ddsp[33];

    auto g_xy_zz_0_y = buffer_ddsp[34];

    auto g_xy_zz_0_z = buffer_ddsp[35];

    auto g_xz_xx_0_x = buffer_ddsp[36];

    auto g_xz_xx_0_y = buffer_ddsp[37];

    auto g_xz_xx_0_z = buffer_ddsp[38];

    auto g_xz_xy_0_x = buffer_ddsp[39];

    auto g_xz_xy_0_y = buffer_ddsp[40];

    auto g_xz_xy_0_z = buffer_ddsp[41];

    auto g_xz_xz_0_x = buffer_ddsp[42];

    auto g_xz_xz_0_y = buffer_ddsp[43];

    auto g_xz_xz_0_z = buffer_ddsp[44];

    auto g_xz_yy_0_x = buffer_ddsp[45];

    auto g_xz_yy_0_y = buffer_ddsp[46];

    auto g_xz_yy_0_z = buffer_ddsp[47];

    auto g_xz_yz_0_x = buffer_ddsp[48];

    auto g_xz_yz_0_y = buffer_ddsp[49];

    auto g_xz_yz_0_z = buffer_ddsp[50];

    auto g_xz_zz_0_x = buffer_ddsp[51];

    auto g_xz_zz_0_y = buffer_ddsp[52];

    auto g_xz_zz_0_z = buffer_ddsp[53];

    auto g_yy_xx_0_x = buffer_ddsp[54];

    auto g_yy_xx_0_y = buffer_ddsp[55];

    auto g_yy_xx_0_z = buffer_ddsp[56];

    auto g_yy_xy_0_x = buffer_ddsp[57];

    auto g_yy_xy_0_y = buffer_ddsp[58];

    auto g_yy_xy_0_z = buffer_ddsp[59];

    auto g_yy_xz_0_x = buffer_ddsp[60];

    auto g_yy_xz_0_y = buffer_ddsp[61];

    auto g_yy_xz_0_z = buffer_ddsp[62];

    auto g_yy_yy_0_x = buffer_ddsp[63];

    auto g_yy_yy_0_y = buffer_ddsp[64];

    auto g_yy_yy_0_z = buffer_ddsp[65];

    auto g_yy_yz_0_x = buffer_ddsp[66];

    auto g_yy_yz_0_y = buffer_ddsp[67];

    auto g_yy_yz_0_z = buffer_ddsp[68];

    auto g_yy_zz_0_x = buffer_ddsp[69];

    auto g_yy_zz_0_y = buffer_ddsp[70];

    auto g_yy_zz_0_z = buffer_ddsp[71];

    auto g_yz_xx_0_x = buffer_ddsp[72];

    auto g_yz_xx_0_y = buffer_ddsp[73];

    auto g_yz_xx_0_z = buffer_ddsp[74];

    auto g_yz_xy_0_x = buffer_ddsp[75];

    auto g_yz_xy_0_y = buffer_ddsp[76];

    auto g_yz_xy_0_z = buffer_ddsp[77];

    auto g_yz_xz_0_x = buffer_ddsp[78];

    auto g_yz_xz_0_y = buffer_ddsp[79];

    auto g_yz_xz_0_z = buffer_ddsp[80];

    auto g_yz_yy_0_x = buffer_ddsp[81];

    auto g_yz_yy_0_y = buffer_ddsp[82];

    auto g_yz_yy_0_z = buffer_ddsp[83];

    auto g_yz_yz_0_x = buffer_ddsp[84];

    auto g_yz_yz_0_y = buffer_ddsp[85];

    auto g_yz_yz_0_z = buffer_ddsp[86];

    auto g_yz_zz_0_x = buffer_ddsp[87];

    auto g_yz_zz_0_y = buffer_ddsp[88];

    auto g_yz_zz_0_z = buffer_ddsp[89];

    auto g_zz_xx_0_x = buffer_ddsp[90];

    auto g_zz_xx_0_y = buffer_ddsp[91];

    auto g_zz_xx_0_z = buffer_ddsp[92];

    auto g_zz_xy_0_x = buffer_ddsp[93];

    auto g_zz_xy_0_y = buffer_ddsp[94];

    auto g_zz_xy_0_z = buffer_ddsp[95];

    auto g_zz_xz_0_x = buffer_ddsp[96];

    auto g_zz_xz_0_y = buffer_ddsp[97];

    auto g_zz_xz_0_z = buffer_ddsp[98];

    auto g_zz_yy_0_x = buffer_ddsp[99];

    auto g_zz_yy_0_y = buffer_ddsp[100];

    auto g_zz_yy_0_z = buffer_ddsp[101];

    auto g_zz_yz_0_x = buffer_ddsp[102];

    auto g_zz_yz_0_y = buffer_ddsp[103];

    auto g_zz_yz_0_z = buffer_ddsp[104];

    auto g_zz_zz_0_x = buffer_ddsp[105];

    auto g_zz_zz_0_y = buffer_ddsp[106];

    auto g_zz_zz_0_z = buffer_ddsp[107];

    /// Set up components of integrals buffer : buffer_1100_ppsp

    auto g_x_x_0_0_x_x_0_x = buffer_1100_ppsp[0];

    auto g_x_x_0_0_x_x_0_y = buffer_1100_ppsp[1];

    auto g_x_x_0_0_x_x_0_z = buffer_1100_ppsp[2];

    auto g_x_x_0_0_x_y_0_x = buffer_1100_ppsp[3];

    auto g_x_x_0_0_x_y_0_y = buffer_1100_ppsp[4];

    auto g_x_x_0_0_x_y_0_z = buffer_1100_ppsp[5];

    auto g_x_x_0_0_x_z_0_x = buffer_1100_ppsp[6];

    auto g_x_x_0_0_x_z_0_y = buffer_1100_ppsp[7];

    auto g_x_x_0_0_x_z_0_z = buffer_1100_ppsp[8];

    auto g_x_x_0_0_y_x_0_x = buffer_1100_ppsp[9];

    auto g_x_x_0_0_y_x_0_y = buffer_1100_ppsp[10];

    auto g_x_x_0_0_y_x_0_z = buffer_1100_ppsp[11];

    auto g_x_x_0_0_y_y_0_x = buffer_1100_ppsp[12];

    auto g_x_x_0_0_y_y_0_y = buffer_1100_ppsp[13];

    auto g_x_x_0_0_y_y_0_z = buffer_1100_ppsp[14];

    auto g_x_x_0_0_y_z_0_x = buffer_1100_ppsp[15];

    auto g_x_x_0_0_y_z_0_y = buffer_1100_ppsp[16];

    auto g_x_x_0_0_y_z_0_z = buffer_1100_ppsp[17];

    auto g_x_x_0_0_z_x_0_x = buffer_1100_ppsp[18];

    auto g_x_x_0_0_z_x_0_y = buffer_1100_ppsp[19];

    auto g_x_x_0_0_z_x_0_z = buffer_1100_ppsp[20];

    auto g_x_x_0_0_z_y_0_x = buffer_1100_ppsp[21];

    auto g_x_x_0_0_z_y_0_y = buffer_1100_ppsp[22];

    auto g_x_x_0_0_z_y_0_z = buffer_1100_ppsp[23];

    auto g_x_x_0_0_z_z_0_x = buffer_1100_ppsp[24];

    auto g_x_x_0_0_z_z_0_y = buffer_1100_ppsp[25];

    auto g_x_x_0_0_z_z_0_z = buffer_1100_ppsp[26];

    auto g_x_y_0_0_x_x_0_x = buffer_1100_ppsp[27];

    auto g_x_y_0_0_x_x_0_y = buffer_1100_ppsp[28];

    auto g_x_y_0_0_x_x_0_z = buffer_1100_ppsp[29];

    auto g_x_y_0_0_x_y_0_x = buffer_1100_ppsp[30];

    auto g_x_y_0_0_x_y_0_y = buffer_1100_ppsp[31];

    auto g_x_y_0_0_x_y_0_z = buffer_1100_ppsp[32];

    auto g_x_y_0_0_x_z_0_x = buffer_1100_ppsp[33];

    auto g_x_y_0_0_x_z_0_y = buffer_1100_ppsp[34];

    auto g_x_y_0_0_x_z_0_z = buffer_1100_ppsp[35];

    auto g_x_y_0_0_y_x_0_x = buffer_1100_ppsp[36];

    auto g_x_y_0_0_y_x_0_y = buffer_1100_ppsp[37];

    auto g_x_y_0_0_y_x_0_z = buffer_1100_ppsp[38];

    auto g_x_y_0_0_y_y_0_x = buffer_1100_ppsp[39];

    auto g_x_y_0_0_y_y_0_y = buffer_1100_ppsp[40];

    auto g_x_y_0_0_y_y_0_z = buffer_1100_ppsp[41];

    auto g_x_y_0_0_y_z_0_x = buffer_1100_ppsp[42];

    auto g_x_y_0_0_y_z_0_y = buffer_1100_ppsp[43];

    auto g_x_y_0_0_y_z_0_z = buffer_1100_ppsp[44];

    auto g_x_y_0_0_z_x_0_x = buffer_1100_ppsp[45];

    auto g_x_y_0_0_z_x_0_y = buffer_1100_ppsp[46];

    auto g_x_y_0_0_z_x_0_z = buffer_1100_ppsp[47];

    auto g_x_y_0_0_z_y_0_x = buffer_1100_ppsp[48];

    auto g_x_y_0_0_z_y_0_y = buffer_1100_ppsp[49];

    auto g_x_y_0_0_z_y_0_z = buffer_1100_ppsp[50];

    auto g_x_y_0_0_z_z_0_x = buffer_1100_ppsp[51];

    auto g_x_y_0_0_z_z_0_y = buffer_1100_ppsp[52];

    auto g_x_y_0_0_z_z_0_z = buffer_1100_ppsp[53];

    auto g_x_z_0_0_x_x_0_x = buffer_1100_ppsp[54];

    auto g_x_z_0_0_x_x_0_y = buffer_1100_ppsp[55];

    auto g_x_z_0_0_x_x_0_z = buffer_1100_ppsp[56];

    auto g_x_z_0_0_x_y_0_x = buffer_1100_ppsp[57];

    auto g_x_z_0_0_x_y_0_y = buffer_1100_ppsp[58];

    auto g_x_z_0_0_x_y_0_z = buffer_1100_ppsp[59];

    auto g_x_z_0_0_x_z_0_x = buffer_1100_ppsp[60];

    auto g_x_z_0_0_x_z_0_y = buffer_1100_ppsp[61];

    auto g_x_z_0_0_x_z_0_z = buffer_1100_ppsp[62];

    auto g_x_z_0_0_y_x_0_x = buffer_1100_ppsp[63];

    auto g_x_z_0_0_y_x_0_y = buffer_1100_ppsp[64];

    auto g_x_z_0_0_y_x_0_z = buffer_1100_ppsp[65];

    auto g_x_z_0_0_y_y_0_x = buffer_1100_ppsp[66];

    auto g_x_z_0_0_y_y_0_y = buffer_1100_ppsp[67];

    auto g_x_z_0_0_y_y_0_z = buffer_1100_ppsp[68];

    auto g_x_z_0_0_y_z_0_x = buffer_1100_ppsp[69];

    auto g_x_z_0_0_y_z_0_y = buffer_1100_ppsp[70];

    auto g_x_z_0_0_y_z_0_z = buffer_1100_ppsp[71];

    auto g_x_z_0_0_z_x_0_x = buffer_1100_ppsp[72];

    auto g_x_z_0_0_z_x_0_y = buffer_1100_ppsp[73];

    auto g_x_z_0_0_z_x_0_z = buffer_1100_ppsp[74];

    auto g_x_z_0_0_z_y_0_x = buffer_1100_ppsp[75];

    auto g_x_z_0_0_z_y_0_y = buffer_1100_ppsp[76];

    auto g_x_z_0_0_z_y_0_z = buffer_1100_ppsp[77];

    auto g_x_z_0_0_z_z_0_x = buffer_1100_ppsp[78];

    auto g_x_z_0_0_z_z_0_y = buffer_1100_ppsp[79];

    auto g_x_z_0_0_z_z_0_z = buffer_1100_ppsp[80];

    auto g_y_x_0_0_x_x_0_x = buffer_1100_ppsp[81];

    auto g_y_x_0_0_x_x_0_y = buffer_1100_ppsp[82];

    auto g_y_x_0_0_x_x_0_z = buffer_1100_ppsp[83];

    auto g_y_x_0_0_x_y_0_x = buffer_1100_ppsp[84];

    auto g_y_x_0_0_x_y_0_y = buffer_1100_ppsp[85];

    auto g_y_x_0_0_x_y_0_z = buffer_1100_ppsp[86];

    auto g_y_x_0_0_x_z_0_x = buffer_1100_ppsp[87];

    auto g_y_x_0_0_x_z_0_y = buffer_1100_ppsp[88];

    auto g_y_x_0_0_x_z_0_z = buffer_1100_ppsp[89];

    auto g_y_x_0_0_y_x_0_x = buffer_1100_ppsp[90];

    auto g_y_x_0_0_y_x_0_y = buffer_1100_ppsp[91];

    auto g_y_x_0_0_y_x_0_z = buffer_1100_ppsp[92];

    auto g_y_x_0_0_y_y_0_x = buffer_1100_ppsp[93];

    auto g_y_x_0_0_y_y_0_y = buffer_1100_ppsp[94];

    auto g_y_x_0_0_y_y_0_z = buffer_1100_ppsp[95];

    auto g_y_x_0_0_y_z_0_x = buffer_1100_ppsp[96];

    auto g_y_x_0_0_y_z_0_y = buffer_1100_ppsp[97];

    auto g_y_x_0_0_y_z_0_z = buffer_1100_ppsp[98];

    auto g_y_x_0_0_z_x_0_x = buffer_1100_ppsp[99];

    auto g_y_x_0_0_z_x_0_y = buffer_1100_ppsp[100];

    auto g_y_x_0_0_z_x_0_z = buffer_1100_ppsp[101];

    auto g_y_x_0_0_z_y_0_x = buffer_1100_ppsp[102];

    auto g_y_x_0_0_z_y_0_y = buffer_1100_ppsp[103];

    auto g_y_x_0_0_z_y_0_z = buffer_1100_ppsp[104];

    auto g_y_x_0_0_z_z_0_x = buffer_1100_ppsp[105];

    auto g_y_x_0_0_z_z_0_y = buffer_1100_ppsp[106];

    auto g_y_x_0_0_z_z_0_z = buffer_1100_ppsp[107];

    auto g_y_y_0_0_x_x_0_x = buffer_1100_ppsp[108];

    auto g_y_y_0_0_x_x_0_y = buffer_1100_ppsp[109];

    auto g_y_y_0_0_x_x_0_z = buffer_1100_ppsp[110];

    auto g_y_y_0_0_x_y_0_x = buffer_1100_ppsp[111];

    auto g_y_y_0_0_x_y_0_y = buffer_1100_ppsp[112];

    auto g_y_y_0_0_x_y_0_z = buffer_1100_ppsp[113];

    auto g_y_y_0_0_x_z_0_x = buffer_1100_ppsp[114];

    auto g_y_y_0_0_x_z_0_y = buffer_1100_ppsp[115];

    auto g_y_y_0_0_x_z_0_z = buffer_1100_ppsp[116];

    auto g_y_y_0_0_y_x_0_x = buffer_1100_ppsp[117];

    auto g_y_y_0_0_y_x_0_y = buffer_1100_ppsp[118];

    auto g_y_y_0_0_y_x_0_z = buffer_1100_ppsp[119];

    auto g_y_y_0_0_y_y_0_x = buffer_1100_ppsp[120];

    auto g_y_y_0_0_y_y_0_y = buffer_1100_ppsp[121];

    auto g_y_y_0_0_y_y_0_z = buffer_1100_ppsp[122];

    auto g_y_y_0_0_y_z_0_x = buffer_1100_ppsp[123];

    auto g_y_y_0_0_y_z_0_y = buffer_1100_ppsp[124];

    auto g_y_y_0_0_y_z_0_z = buffer_1100_ppsp[125];

    auto g_y_y_0_0_z_x_0_x = buffer_1100_ppsp[126];

    auto g_y_y_0_0_z_x_0_y = buffer_1100_ppsp[127];

    auto g_y_y_0_0_z_x_0_z = buffer_1100_ppsp[128];

    auto g_y_y_0_0_z_y_0_x = buffer_1100_ppsp[129];

    auto g_y_y_0_0_z_y_0_y = buffer_1100_ppsp[130];

    auto g_y_y_0_0_z_y_0_z = buffer_1100_ppsp[131];

    auto g_y_y_0_0_z_z_0_x = buffer_1100_ppsp[132];

    auto g_y_y_0_0_z_z_0_y = buffer_1100_ppsp[133];

    auto g_y_y_0_0_z_z_0_z = buffer_1100_ppsp[134];

    auto g_y_z_0_0_x_x_0_x = buffer_1100_ppsp[135];

    auto g_y_z_0_0_x_x_0_y = buffer_1100_ppsp[136];

    auto g_y_z_0_0_x_x_0_z = buffer_1100_ppsp[137];

    auto g_y_z_0_0_x_y_0_x = buffer_1100_ppsp[138];

    auto g_y_z_0_0_x_y_0_y = buffer_1100_ppsp[139];

    auto g_y_z_0_0_x_y_0_z = buffer_1100_ppsp[140];

    auto g_y_z_0_0_x_z_0_x = buffer_1100_ppsp[141];

    auto g_y_z_0_0_x_z_0_y = buffer_1100_ppsp[142];

    auto g_y_z_0_0_x_z_0_z = buffer_1100_ppsp[143];

    auto g_y_z_0_0_y_x_0_x = buffer_1100_ppsp[144];

    auto g_y_z_0_0_y_x_0_y = buffer_1100_ppsp[145];

    auto g_y_z_0_0_y_x_0_z = buffer_1100_ppsp[146];

    auto g_y_z_0_0_y_y_0_x = buffer_1100_ppsp[147];

    auto g_y_z_0_0_y_y_0_y = buffer_1100_ppsp[148];

    auto g_y_z_0_0_y_y_0_z = buffer_1100_ppsp[149];

    auto g_y_z_0_0_y_z_0_x = buffer_1100_ppsp[150];

    auto g_y_z_0_0_y_z_0_y = buffer_1100_ppsp[151];

    auto g_y_z_0_0_y_z_0_z = buffer_1100_ppsp[152];

    auto g_y_z_0_0_z_x_0_x = buffer_1100_ppsp[153];

    auto g_y_z_0_0_z_x_0_y = buffer_1100_ppsp[154];

    auto g_y_z_0_0_z_x_0_z = buffer_1100_ppsp[155];

    auto g_y_z_0_0_z_y_0_x = buffer_1100_ppsp[156];

    auto g_y_z_0_0_z_y_0_y = buffer_1100_ppsp[157];

    auto g_y_z_0_0_z_y_0_z = buffer_1100_ppsp[158];

    auto g_y_z_0_0_z_z_0_x = buffer_1100_ppsp[159];

    auto g_y_z_0_0_z_z_0_y = buffer_1100_ppsp[160];

    auto g_y_z_0_0_z_z_0_z = buffer_1100_ppsp[161];

    auto g_z_x_0_0_x_x_0_x = buffer_1100_ppsp[162];

    auto g_z_x_0_0_x_x_0_y = buffer_1100_ppsp[163];

    auto g_z_x_0_0_x_x_0_z = buffer_1100_ppsp[164];

    auto g_z_x_0_0_x_y_0_x = buffer_1100_ppsp[165];

    auto g_z_x_0_0_x_y_0_y = buffer_1100_ppsp[166];

    auto g_z_x_0_0_x_y_0_z = buffer_1100_ppsp[167];

    auto g_z_x_0_0_x_z_0_x = buffer_1100_ppsp[168];

    auto g_z_x_0_0_x_z_0_y = buffer_1100_ppsp[169];

    auto g_z_x_0_0_x_z_0_z = buffer_1100_ppsp[170];

    auto g_z_x_0_0_y_x_0_x = buffer_1100_ppsp[171];

    auto g_z_x_0_0_y_x_0_y = buffer_1100_ppsp[172];

    auto g_z_x_0_0_y_x_0_z = buffer_1100_ppsp[173];

    auto g_z_x_0_0_y_y_0_x = buffer_1100_ppsp[174];

    auto g_z_x_0_0_y_y_0_y = buffer_1100_ppsp[175];

    auto g_z_x_0_0_y_y_0_z = buffer_1100_ppsp[176];

    auto g_z_x_0_0_y_z_0_x = buffer_1100_ppsp[177];

    auto g_z_x_0_0_y_z_0_y = buffer_1100_ppsp[178];

    auto g_z_x_0_0_y_z_0_z = buffer_1100_ppsp[179];

    auto g_z_x_0_0_z_x_0_x = buffer_1100_ppsp[180];

    auto g_z_x_0_0_z_x_0_y = buffer_1100_ppsp[181];

    auto g_z_x_0_0_z_x_0_z = buffer_1100_ppsp[182];

    auto g_z_x_0_0_z_y_0_x = buffer_1100_ppsp[183];

    auto g_z_x_0_0_z_y_0_y = buffer_1100_ppsp[184];

    auto g_z_x_0_0_z_y_0_z = buffer_1100_ppsp[185];

    auto g_z_x_0_0_z_z_0_x = buffer_1100_ppsp[186];

    auto g_z_x_0_0_z_z_0_y = buffer_1100_ppsp[187];

    auto g_z_x_0_0_z_z_0_z = buffer_1100_ppsp[188];

    auto g_z_y_0_0_x_x_0_x = buffer_1100_ppsp[189];

    auto g_z_y_0_0_x_x_0_y = buffer_1100_ppsp[190];

    auto g_z_y_0_0_x_x_0_z = buffer_1100_ppsp[191];

    auto g_z_y_0_0_x_y_0_x = buffer_1100_ppsp[192];

    auto g_z_y_0_0_x_y_0_y = buffer_1100_ppsp[193];

    auto g_z_y_0_0_x_y_0_z = buffer_1100_ppsp[194];

    auto g_z_y_0_0_x_z_0_x = buffer_1100_ppsp[195];

    auto g_z_y_0_0_x_z_0_y = buffer_1100_ppsp[196];

    auto g_z_y_0_0_x_z_0_z = buffer_1100_ppsp[197];

    auto g_z_y_0_0_y_x_0_x = buffer_1100_ppsp[198];

    auto g_z_y_0_0_y_x_0_y = buffer_1100_ppsp[199];

    auto g_z_y_0_0_y_x_0_z = buffer_1100_ppsp[200];

    auto g_z_y_0_0_y_y_0_x = buffer_1100_ppsp[201];

    auto g_z_y_0_0_y_y_0_y = buffer_1100_ppsp[202];

    auto g_z_y_0_0_y_y_0_z = buffer_1100_ppsp[203];

    auto g_z_y_0_0_y_z_0_x = buffer_1100_ppsp[204];

    auto g_z_y_0_0_y_z_0_y = buffer_1100_ppsp[205];

    auto g_z_y_0_0_y_z_0_z = buffer_1100_ppsp[206];

    auto g_z_y_0_0_z_x_0_x = buffer_1100_ppsp[207];

    auto g_z_y_0_0_z_x_0_y = buffer_1100_ppsp[208];

    auto g_z_y_0_0_z_x_0_z = buffer_1100_ppsp[209];

    auto g_z_y_0_0_z_y_0_x = buffer_1100_ppsp[210];

    auto g_z_y_0_0_z_y_0_y = buffer_1100_ppsp[211];

    auto g_z_y_0_0_z_y_0_z = buffer_1100_ppsp[212];

    auto g_z_y_0_0_z_z_0_x = buffer_1100_ppsp[213];

    auto g_z_y_0_0_z_z_0_y = buffer_1100_ppsp[214];

    auto g_z_y_0_0_z_z_0_z = buffer_1100_ppsp[215];

    auto g_z_z_0_0_x_x_0_x = buffer_1100_ppsp[216];

    auto g_z_z_0_0_x_x_0_y = buffer_1100_ppsp[217];

    auto g_z_z_0_0_x_x_0_z = buffer_1100_ppsp[218];

    auto g_z_z_0_0_x_y_0_x = buffer_1100_ppsp[219];

    auto g_z_z_0_0_x_y_0_y = buffer_1100_ppsp[220];

    auto g_z_z_0_0_x_y_0_z = buffer_1100_ppsp[221];

    auto g_z_z_0_0_x_z_0_x = buffer_1100_ppsp[222];

    auto g_z_z_0_0_x_z_0_y = buffer_1100_ppsp[223];

    auto g_z_z_0_0_x_z_0_z = buffer_1100_ppsp[224];

    auto g_z_z_0_0_y_x_0_x = buffer_1100_ppsp[225];

    auto g_z_z_0_0_y_x_0_y = buffer_1100_ppsp[226];

    auto g_z_z_0_0_y_x_0_z = buffer_1100_ppsp[227];

    auto g_z_z_0_0_y_y_0_x = buffer_1100_ppsp[228];

    auto g_z_z_0_0_y_y_0_y = buffer_1100_ppsp[229];

    auto g_z_z_0_0_y_y_0_z = buffer_1100_ppsp[230];

    auto g_z_z_0_0_y_z_0_x = buffer_1100_ppsp[231];

    auto g_z_z_0_0_y_z_0_y = buffer_1100_ppsp[232];

    auto g_z_z_0_0_y_z_0_z = buffer_1100_ppsp[233];

    auto g_z_z_0_0_z_x_0_x = buffer_1100_ppsp[234];

    auto g_z_z_0_0_z_x_0_y = buffer_1100_ppsp[235];

    auto g_z_z_0_0_z_x_0_z = buffer_1100_ppsp[236];

    auto g_z_z_0_0_z_y_0_x = buffer_1100_ppsp[237];

    auto g_z_z_0_0_z_y_0_y = buffer_1100_ppsp[238];

    auto g_z_z_0_0_z_y_0_z = buffer_1100_ppsp[239];

    auto g_z_z_0_0_z_z_0_x = buffer_1100_ppsp[240];

    auto g_z_z_0_0_z_z_0_y = buffer_1100_ppsp[241];

    auto g_z_z_0_0_z_z_0_z = buffer_1100_ppsp[242];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_0_0_x, g_0_0_0_y, g_0_0_0_z, g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_x_x_0_0_x_x_0_x, g_x_x_0_0_x_x_0_y, g_x_x_0_0_x_x_0_z, g_xx_0_0_x, g_xx_0_0_y, g_xx_0_0_z, g_xx_xx_0_x, g_xx_xx_0_y, g_xx_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_x_0_x[i] = g_0_0_0_x[i] - 2.0 * g_0_xx_0_x[i] * b_exp - 2.0 * g_xx_0_0_x[i] * a_exp + 4.0 * g_xx_xx_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_0_y[i] = g_0_0_0_y[i] - 2.0 * g_0_xx_0_y[i] * b_exp - 2.0 * g_xx_0_0_y[i] * a_exp + 4.0 * g_xx_xx_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_x_0_z[i] = g_0_0_0_z[i] - 2.0 * g_0_xx_0_z[i] * b_exp - 2.0 * g_xx_0_0_z[i] * a_exp + 4.0 * g_xx_xx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_x_x_0_0_x_y_0_x, g_x_x_0_0_x_y_0_y, g_x_x_0_0_x_y_0_z, g_xx_xy_0_x, g_xx_xy_0_y, g_xx_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_y_0_x[i] = -2.0 * g_0_xy_0_x[i] * b_exp + 4.0 * g_xx_xy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_0_y[i] = -2.0 * g_0_xy_0_y[i] * b_exp + 4.0 * g_xx_xy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_0_z[i] = -2.0 * g_0_xy_0_z[i] * b_exp + 4.0 * g_xx_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_x_x_0_0_x_z_0_x, g_x_x_0_0_x_z_0_y, g_x_x_0_0_x_z_0_z, g_xx_xz_0_x, g_xx_xz_0_y, g_xx_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_z_0_x[i] = -2.0 * g_0_xz_0_x[i] * b_exp + 4.0 * g_xx_xz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_0_y[i] = -2.0 * g_0_xz_0_y[i] * b_exp + 4.0 * g_xx_xz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_0_z[i] = -2.0 * g_0_xz_0_z[i] * b_exp + 4.0 * g_xx_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_x_0_0_y_x_0_x, g_x_x_0_0_y_x_0_y, g_x_x_0_0_y_x_0_z, g_xy_0_0_x, g_xy_0_0_y, g_xy_0_0_z, g_xy_xx_0_x, g_xy_xx_0_y, g_xy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_x_0_x[i] = -2.0 * g_xy_0_0_x[i] * a_exp + 4.0 * g_xy_xx_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_0_y[i] = -2.0 * g_xy_0_0_y[i] * a_exp + 4.0 * g_xy_xx_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_x_0_z[i] = -2.0 * g_xy_0_0_z[i] * a_exp + 4.0 * g_xy_xx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_x_0_0_y_y_0_x, g_x_x_0_0_y_y_0_y, g_x_x_0_0_y_y_0_z, g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_y_0_x[i] = 4.0 * g_xy_xy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_0_y[i] = 4.0 * g_xy_xy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_0_z[i] = 4.0 * g_xy_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_x_0_0_y_z_0_x, g_x_x_0_0_y_z_0_y, g_x_x_0_0_y_z_0_z, g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_z_0_x[i] = 4.0 * g_xy_xz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_0_y[i] = 4.0 * g_xy_xz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_0_z[i] = 4.0 * g_xy_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_x_0_0_z_x_0_x, g_x_x_0_0_z_x_0_y, g_x_x_0_0_z_x_0_z, g_xz_0_0_x, g_xz_0_0_y, g_xz_0_0_z, g_xz_xx_0_x, g_xz_xx_0_y, g_xz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_x_0_x[i] = -2.0 * g_xz_0_0_x[i] * a_exp + 4.0 * g_xz_xx_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_0_y[i] = -2.0 * g_xz_0_0_y[i] * a_exp + 4.0 * g_xz_xx_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_x_0_z[i] = -2.0 * g_xz_0_0_z[i] * a_exp + 4.0 * g_xz_xx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_x_0_0_z_y_0_x, g_x_x_0_0_z_y_0_y, g_x_x_0_0_z_y_0_z, g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_y_0_x[i] = 4.0 * g_xz_xy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_0_y[i] = 4.0 * g_xz_xy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_0_z[i] = 4.0 * g_xz_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_x_0_0_z_z_0_x, g_x_x_0_0_z_z_0_y, g_x_x_0_0_z_z_0_z, g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_z_0_x[i] = 4.0 * g_xz_xz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_0_y[i] = 4.0 * g_xz_xz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_0_z[i] = 4.0 * g_xz_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_x_y_0_0_x_x_0_x, g_x_y_0_0_x_x_0_y, g_x_y_0_0_x_x_0_z, g_xx_xy_0_x, g_xx_xy_0_y, g_xx_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_x_0_x[i] = -2.0 * g_0_xy_0_x[i] * b_exp + 4.0 * g_xx_xy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_0_y[i] = -2.0 * g_0_xy_0_y[i] * b_exp + 4.0 * g_xx_xy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_x_0_z[i] = -2.0 * g_0_xy_0_z[i] * b_exp + 4.0 * g_xx_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_0_0_0_x, g_0_0_0_y, g_0_0_0_z, g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_x_y_0_0_x_y_0_x, g_x_y_0_0_x_y_0_y, g_x_y_0_0_x_y_0_z, g_xx_0_0_x, g_xx_0_0_y, g_xx_0_0_z, g_xx_yy_0_x, g_xx_yy_0_y, g_xx_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_y_0_x[i] = g_0_0_0_x[i] - 2.0 * g_0_yy_0_x[i] * b_exp - 2.0 * g_xx_0_0_x[i] * a_exp + 4.0 * g_xx_yy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_0_y[i] = g_0_0_0_y[i] - 2.0 * g_0_yy_0_y[i] * b_exp - 2.0 * g_xx_0_0_y[i] * a_exp + 4.0 * g_xx_yy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_0_z[i] = g_0_0_0_z[i] - 2.0 * g_0_yy_0_z[i] * b_exp - 2.0 * g_xx_0_0_z[i] * a_exp + 4.0 * g_xx_yy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_x_y_0_0_x_z_0_x, g_x_y_0_0_x_z_0_y, g_x_y_0_0_x_z_0_z, g_xx_yz_0_x, g_xx_yz_0_y, g_xx_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_z_0_x[i] = -2.0 * g_0_yz_0_x[i] * b_exp + 4.0 * g_xx_yz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_0_y[i] = -2.0 * g_0_yz_0_y[i] * b_exp + 4.0 * g_xx_yz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_0_z[i] = -2.0 * g_0_yz_0_z[i] * b_exp + 4.0 * g_xx_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_y_0_0_y_x_0_x, g_x_y_0_0_y_x_0_y, g_x_y_0_0_y_x_0_z, g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_x_0_x[i] = 4.0 * g_xy_xy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_0_y[i] = 4.0 * g_xy_xy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_x_0_z[i] = 4.0 * g_xy_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_y_0_0_y_y_0_x, g_x_y_0_0_y_y_0_y, g_x_y_0_0_y_y_0_z, g_xy_0_0_x, g_xy_0_0_y, g_xy_0_0_z, g_xy_yy_0_x, g_xy_yy_0_y, g_xy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_y_0_x[i] = -2.0 * g_xy_0_0_x[i] * a_exp + 4.0 * g_xy_yy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_0_y[i] = -2.0 * g_xy_0_0_y[i] * a_exp + 4.0 * g_xy_yy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_0_z[i] = -2.0 * g_xy_0_0_z[i] * a_exp + 4.0 * g_xy_yy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_y_0_0_y_z_0_x, g_x_y_0_0_y_z_0_y, g_x_y_0_0_y_z_0_z, g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_z_0_x[i] = 4.0 * g_xy_yz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_0_y[i] = 4.0 * g_xy_yz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_0_z[i] = 4.0 * g_xy_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_y_0_0_z_x_0_x, g_x_y_0_0_z_x_0_y, g_x_y_0_0_z_x_0_z, g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_x_0_x[i] = 4.0 * g_xz_xy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_0_y[i] = 4.0 * g_xz_xy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_x_0_z[i] = 4.0 * g_xz_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_y_0_0_z_y_0_x, g_x_y_0_0_z_y_0_y, g_x_y_0_0_z_y_0_z, g_xz_0_0_x, g_xz_0_0_y, g_xz_0_0_z, g_xz_yy_0_x, g_xz_yy_0_y, g_xz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_y_0_x[i] = -2.0 * g_xz_0_0_x[i] * a_exp + 4.0 * g_xz_yy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_0_y[i] = -2.0 * g_xz_0_0_y[i] * a_exp + 4.0 * g_xz_yy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_0_z[i] = -2.0 * g_xz_0_0_z[i] * a_exp + 4.0 * g_xz_yy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_y_0_0_z_z_0_x, g_x_y_0_0_z_z_0_y, g_x_y_0_0_z_z_0_z, g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_z_0_x[i] = 4.0 * g_xz_yz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_0_y[i] = 4.0 * g_xz_yz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_0_z[i] = 4.0 * g_xz_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_x_z_0_0_x_x_0_x, g_x_z_0_0_x_x_0_y, g_x_z_0_0_x_x_0_z, g_xx_xz_0_x, g_xx_xz_0_y, g_xx_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_x_0_x[i] = -2.0 * g_0_xz_0_x[i] * b_exp + 4.0 * g_xx_xz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_0_y[i] = -2.0 * g_0_xz_0_y[i] * b_exp + 4.0 * g_xx_xz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_x_0_z[i] = -2.0 * g_0_xz_0_z[i] * b_exp + 4.0 * g_xx_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_x_z_0_0_x_y_0_x, g_x_z_0_0_x_y_0_y, g_x_z_0_0_x_y_0_z, g_xx_yz_0_x, g_xx_yz_0_y, g_xx_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_y_0_x[i] = -2.0 * g_0_yz_0_x[i] * b_exp + 4.0 * g_xx_yz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_0_y[i] = -2.0 * g_0_yz_0_y[i] * b_exp + 4.0 * g_xx_yz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_0_z[i] = -2.0 * g_0_yz_0_z[i] * b_exp + 4.0 * g_xx_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_0_0_0_x, g_0_0_0_y, g_0_0_0_z, g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_x_z_0_0_x_z_0_x, g_x_z_0_0_x_z_0_y, g_x_z_0_0_x_z_0_z, g_xx_0_0_x, g_xx_0_0_y, g_xx_0_0_z, g_xx_zz_0_x, g_xx_zz_0_y, g_xx_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_z_0_x[i] = g_0_0_0_x[i] - 2.0 * g_0_zz_0_x[i] * b_exp - 2.0 * g_xx_0_0_x[i] * a_exp + 4.0 * g_xx_zz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_0_y[i] = g_0_0_0_y[i] - 2.0 * g_0_zz_0_y[i] * b_exp - 2.0 * g_xx_0_0_y[i] * a_exp + 4.0 * g_xx_zz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_0_z[i] = g_0_0_0_z[i] - 2.0 * g_0_zz_0_z[i] * b_exp - 2.0 * g_xx_0_0_z[i] * a_exp + 4.0 * g_xx_zz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_x_z_0_0_y_x_0_x, g_x_z_0_0_y_x_0_y, g_x_z_0_0_y_x_0_z, g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_x_0_x[i] = 4.0 * g_xy_xz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_0_y[i] = 4.0 * g_xy_xz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_x_0_z[i] = 4.0 * g_xy_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_x_z_0_0_y_y_0_x, g_x_z_0_0_y_y_0_y, g_x_z_0_0_y_y_0_z, g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_y_0_x[i] = 4.0 * g_xy_yz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_0_y[i] = 4.0 * g_xy_yz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_0_z[i] = 4.0 * g_xy_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_x_z_0_0_y_z_0_x, g_x_z_0_0_y_z_0_y, g_x_z_0_0_y_z_0_z, g_xy_0_0_x, g_xy_0_0_y, g_xy_0_0_z, g_xy_zz_0_x, g_xy_zz_0_y, g_xy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_z_0_x[i] = -2.0 * g_xy_0_0_x[i] * a_exp + 4.0 * g_xy_zz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_0_y[i] = -2.0 * g_xy_0_0_y[i] * a_exp + 4.0 * g_xy_zz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_0_z[i] = -2.0 * g_xy_0_0_z[i] * a_exp + 4.0 * g_xy_zz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_z_0_0_z_x_0_x, g_x_z_0_0_z_x_0_y, g_x_z_0_0_z_x_0_z, g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_x_0_x[i] = 4.0 * g_xz_xz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_0_y[i] = 4.0 * g_xz_xz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_x_0_z[i] = 4.0 * g_xz_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_z_0_0_z_y_0_x, g_x_z_0_0_z_y_0_y, g_x_z_0_0_z_y_0_z, g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_y_0_x[i] = 4.0 * g_xz_yz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_0_y[i] = 4.0 * g_xz_yz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_0_z[i] = 4.0 * g_xz_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_z_0_0_z_z_0_x, g_x_z_0_0_z_z_0_y, g_x_z_0_0_z_z_0_z, g_xz_0_0_x, g_xz_0_0_y, g_xz_0_0_z, g_xz_zz_0_x, g_xz_zz_0_y, g_xz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_z_0_x[i] = -2.0 * g_xz_0_0_x[i] * a_exp + 4.0 * g_xz_zz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_0_y[i] = -2.0 * g_xz_0_0_y[i] * a_exp + 4.0 * g_xz_zz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_0_z[i] = -2.0 * g_xz_0_0_z[i] * a_exp + 4.0 * g_xz_zz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_xy_0_0_x, g_xy_0_0_y, g_xy_0_0_z, g_xy_xx_0_x, g_xy_xx_0_y, g_xy_xx_0_z, g_y_x_0_0_x_x_0_x, g_y_x_0_0_x_x_0_y, g_y_x_0_0_x_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_x_0_x[i] = -2.0 * g_xy_0_0_x[i] * a_exp + 4.0 * g_xy_xx_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_0_y[i] = -2.0 * g_xy_0_0_y[i] * a_exp + 4.0 * g_xy_xx_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_x_0_z[i] = -2.0 * g_xy_0_0_z[i] * a_exp + 4.0 * g_xy_xx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z, g_y_x_0_0_x_y_0_x, g_y_x_0_0_x_y_0_y, g_y_x_0_0_x_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_y_0_x[i] = 4.0 * g_xy_xy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_0_y[i] = 4.0 * g_xy_xy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_0_z[i] = 4.0 * g_xy_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z, g_y_x_0_0_x_z_0_x, g_y_x_0_0_x_z_0_y, g_y_x_0_0_x_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_z_0_x[i] = 4.0 * g_xy_xz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_0_y[i] = 4.0 * g_xy_xz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_0_z[i] = 4.0 * g_xy_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_0_0_0_x, g_0_0_0_y, g_0_0_0_z, g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_y_x_0_0_y_x_0_x, g_y_x_0_0_y_x_0_y, g_y_x_0_0_y_x_0_z, g_yy_0_0_x, g_yy_0_0_y, g_yy_0_0_z, g_yy_xx_0_x, g_yy_xx_0_y, g_yy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_x_0_x[i] = g_0_0_0_x[i] - 2.0 * g_0_xx_0_x[i] * b_exp - 2.0 * g_yy_0_0_x[i] * a_exp + 4.0 * g_yy_xx_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_0_y[i] = g_0_0_0_y[i] - 2.0 * g_0_xx_0_y[i] * b_exp - 2.0 * g_yy_0_0_y[i] * a_exp + 4.0 * g_yy_xx_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_x_0_z[i] = g_0_0_0_z[i] - 2.0 * g_0_xx_0_z[i] * b_exp - 2.0 * g_yy_0_0_z[i] * a_exp + 4.0 * g_yy_xx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_y_x_0_0_y_y_0_x, g_y_x_0_0_y_y_0_y, g_y_x_0_0_y_y_0_z, g_yy_xy_0_x, g_yy_xy_0_y, g_yy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_y_0_x[i] = -2.0 * g_0_xy_0_x[i] * b_exp + 4.0 * g_yy_xy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_0_y[i] = -2.0 * g_0_xy_0_y[i] * b_exp + 4.0 * g_yy_xy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_0_z[i] = -2.0 * g_0_xy_0_z[i] * b_exp + 4.0 * g_yy_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_y_x_0_0_y_z_0_x, g_y_x_0_0_y_z_0_y, g_y_x_0_0_y_z_0_z, g_yy_xz_0_x, g_yy_xz_0_y, g_yy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_z_0_x[i] = -2.0 * g_0_xz_0_x[i] * b_exp + 4.0 * g_yy_xz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_0_y[i] = -2.0 * g_0_xz_0_y[i] * b_exp + 4.0 * g_yy_xz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_0_z[i] = -2.0 * g_0_xz_0_z[i] * b_exp + 4.0 * g_yy_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_y_x_0_0_z_x_0_x, g_y_x_0_0_z_x_0_y, g_y_x_0_0_z_x_0_z, g_yz_0_0_x, g_yz_0_0_y, g_yz_0_0_z, g_yz_xx_0_x, g_yz_xx_0_y, g_yz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_x_0_x[i] = -2.0 * g_yz_0_0_x[i] * a_exp + 4.0 * g_yz_xx_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_0_y[i] = -2.0 * g_yz_0_0_y[i] * a_exp + 4.0 * g_yz_xx_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_x_0_z[i] = -2.0 * g_yz_0_0_z[i] * a_exp + 4.0 * g_yz_xx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_y_x_0_0_z_y_0_x, g_y_x_0_0_z_y_0_y, g_y_x_0_0_z_y_0_z, g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_y_0_x[i] = 4.0 * g_yz_xy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_0_y[i] = 4.0 * g_yz_xy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_0_z[i] = 4.0 * g_yz_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_y_x_0_0_z_z_0_x, g_y_x_0_0_z_z_0_y, g_y_x_0_0_z_z_0_z, g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_z_0_x[i] = 4.0 * g_yz_xz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_0_y[i] = 4.0 * g_yz_xz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_0_z[i] = 4.0 * g_yz_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z, g_y_y_0_0_x_x_0_x, g_y_y_0_0_x_x_0_y, g_y_y_0_0_x_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_x_0_x[i] = 4.0 * g_xy_xy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_0_y[i] = 4.0 * g_xy_xy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_x_0_z[i] = 4.0 * g_xy_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_xy_0_0_x, g_xy_0_0_y, g_xy_0_0_z, g_xy_yy_0_x, g_xy_yy_0_y, g_xy_yy_0_z, g_y_y_0_0_x_y_0_x, g_y_y_0_0_x_y_0_y, g_y_y_0_0_x_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_y_0_x[i] = -2.0 * g_xy_0_0_x[i] * a_exp + 4.0 * g_xy_yy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_0_y[i] = -2.0 * g_xy_0_0_y[i] * a_exp + 4.0 * g_xy_yy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_0_z[i] = -2.0 * g_xy_0_0_z[i] * a_exp + 4.0 * g_xy_yy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z, g_y_y_0_0_x_z_0_x, g_y_y_0_0_x_z_0_y, g_y_y_0_0_x_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_z_0_x[i] = 4.0 * g_xy_yz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_0_y[i] = 4.0 * g_xy_yz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_0_z[i] = 4.0 * g_xy_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_y_y_0_0_y_x_0_x, g_y_y_0_0_y_x_0_y, g_y_y_0_0_y_x_0_z, g_yy_xy_0_x, g_yy_xy_0_y, g_yy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_x_0_x[i] = -2.0 * g_0_xy_0_x[i] * b_exp + 4.0 * g_yy_xy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_0_y[i] = -2.0 * g_0_xy_0_y[i] * b_exp + 4.0 * g_yy_xy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_x_0_z[i] = -2.0 * g_0_xy_0_z[i] * b_exp + 4.0 * g_yy_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_0_0_0_x, g_0_0_0_y, g_0_0_0_z, g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_y_y_0_0_y_y_0_x, g_y_y_0_0_y_y_0_y, g_y_y_0_0_y_y_0_z, g_yy_0_0_x, g_yy_0_0_y, g_yy_0_0_z, g_yy_yy_0_x, g_yy_yy_0_y, g_yy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_y_0_x[i] = g_0_0_0_x[i] - 2.0 * g_0_yy_0_x[i] * b_exp - 2.0 * g_yy_0_0_x[i] * a_exp + 4.0 * g_yy_yy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_0_y[i] = g_0_0_0_y[i] - 2.0 * g_0_yy_0_y[i] * b_exp - 2.0 * g_yy_0_0_y[i] * a_exp + 4.0 * g_yy_yy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_0_z[i] = g_0_0_0_z[i] - 2.0 * g_0_yy_0_z[i] * b_exp - 2.0 * g_yy_0_0_z[i] * a_exp + 4.0 * g_yy_yy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_y_y_0_0_y_z_0_x, g_y_y_0_0_y_z_0_y, g_y_y_0_0_y_z_0_z, g_yy_yz_0_x, g_yy_yz_0_y, g_yy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_z_0_x[i] = -2.0 * g_0_yz_0_x[i] * b_exp + 4.0 * g_yy_yz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_0_y[i] = -2.0 * g_0_yz_0_y[i] * b_exp + 4.0 * g_yy_yz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_0_z[i] = -2.0 * g_0_yz_0_z[i] * b_exp + 4.0 * g_yy_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_y_y_0_0_z_x_0_x, g_y_y_0_0_z_x_0_y, g_y_y_0_0_z_x_0_z, g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_x_0_x[i] = 4.0 * g_yz_xy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_0_y[i] = 4.0 * g_yz_xy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_x_0_z[i] = 4.0 * g_yz_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_y_y_0_0_z_y_0_x, g_y_y_0_0_z_y_0_y, g_y_y_0_0_z_y_0_z, g_yz_0_0_x, g_yz_0_0_y, g_yz_0_0_z, g_yz_yy_0_x, g_yz_yy_0_y, g_yz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_y_0_x[i] = -2.0 * g_yz_0_0_x[i] * a_exp + 4.0 * g_yz_yy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_0_y[i] = -2.0 * g_yz_0_0_y[i] * a_exp + 4.0 * g_yz_yy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_0_z[i] = -2.0 * g_yz_0_0_z[i] * a_exp + 4.0 * g_yz_yy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_y_y_0_0_z_z_0_x, g_y_y_0_0_z_z_0_y, g_y_y_0_0_z_z_0_z, g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_z_0_x[i] = 4.0 * g_yz_yz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_0_y[i] = 4.0 * g_yz_yz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_0_z[i] = 4.0 * g_yz_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z, g_y_z_0_0_x_x_0_x, g_y_z_0_0_x_x_0_y, g_y_z_0_0_x_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_x_0_x[i] = 4.0 * g_xy_xz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_0_y[i] = 4.0 * g_xy_xz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_x_0_z[i] = 4.0 * g_xy_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z, g_y_z_0_0_x_y_0_x, g_y_z_0_0_x_y_0_y, g_y_z_0_0_x_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_y_0_x[i] = 4.0 * g_xy_yz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_0_y[i] = 4.0 * g_xy_yz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_0_z[i] = 4.0 * g_xy_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_xy_0_0_x, g_xy_0_0_y, g_xy_0_0_z, g_xy_zz_0_x, g_xy_zz_0_y, g_xy_zz_0_z, g_y_z_0_0_x_z_0_x, g_y_z_0_0_x_z_0_y, g_y_z_0_0_x_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_z_0_x[i] = -2.0 * g_xy_0_0_x[i] * a_exp + 4.0 * g_xy_zz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_0_y[i] = -2.0 * g_xy_0_0_y[i] * a_exp + 4.0 * g_xy_zz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_0_z[i] = -2.0 * g_xy_0_0_z[i] * a_exp + 4.0 * g_xy_zz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_y_z_0_0_y_x_0_x, g_y_z_0_0_y_x_0_y, g_y_z_0_0_y_x_0_z, g_yy_xz_0_x, g_yy_xz_0_y, g_yy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_x_0_x[i] = -2.0 * g_0_xz_0_x[i] * b_exp + 4.0 * g_yy_xz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_0_y[i] = -2.0 * g_0_xz_0_y[i] * b_exp + 4.0 * g_yy_xz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_x_0_z[i] = -2.0 * g_0_xz_0_z[i] * b_exp + 4.0 * g_yy_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_y_z_0_0_y_y_0_x, g_y_z_0_0_y_y_0_y, g_y_z_0_0_y_y_0_z, g_yy_yz_0_x, g_yy_yz_0_y, g_yy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_y_0_x[i] = -2.0 * g_0_yz_0_x[i] * b_exp + 4.0 * g_yy_yz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_0_y[i] = -2.0 * g_0_yz_0_y[i] * b_exp + 4.0 * g_yy_yz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_0_z[i] = -2.0 * g_0_yz_0_z[i] * b_exp + 4.0 * g_yy_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_0_0_0_x, g_0_0_0_y, g_0_0_0_z, g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_y_z_0_0_y_z_0_x, g_y_z_0_0_y_z_0_y, g_y_z_0_0_y_z_0_z, g_yy_0_0_x, g_yy_0_0_y, g_yy_0_0_z, g_yy_zz_0_x, g_yy_zz_0_y, g_yy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_z_0_x[i] = g_0_0_0_x[i] - 2.0 * g_0_zz_0_x[i] * b_exp - 2.0 * g_yy_0_0_x[i] * a_exp + 4.0 * g_yy_zz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_0_y[i] = g_0_0_0_y[i] - 2.0 * g_0_zz_0_y[i] * b_exp - 2.0 * g_yy_0_0_y[i] * a_exp + 4.0 * g_yy_zz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_0_z[i] = g_0_0_0_z[i] - 2.0 * g_0_zz_0_z[i] * b_exp - 2.0 * g_yy_0_0_z[i] * a_exp + 4.0 * g_yy_zz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_y_z_0_0_z_x_0_x, g_y_z_0_0_z_x_0_y, g_y_z_0_0_z_x_0_z, g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_x_0_x[i] = 4.0 * g_yz_xz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_0_y[i] = 4.0 * g_yz_xz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_x_0_z[i] = 4.0 * g_yz_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_y_z_0_0_z_y_0_x, g_y_z_0_0_z_y_0_y, g_y_z_0_0_z_y_0_z, g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_y_0_x[i] = 4.0 * g_yz_yz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_0_y[i] = 4.0 * g_yz_yz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_0_z[i] = 4.0 * g_yz_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_y_z_0_0_z_z_0_x, g_y_z_0_0_z_z_0_y, g_y_z_0_0_z_z_0_z, g_yz_0_0_x, g_yz_0_0_y, g_yz_0_0_z, g_yz_zz_0_x, g_yz_zz_0_y, g_yz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_z_0_x[i] = -2.0 * g_yz_0_0_x[i] * a_exp + 4.0 * g_yz_zz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_0_y[i] = -2.0 * g_yz_0_0_y[i] * a_exp + 4.0 * g_yz_zz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_0_z[i] = -2.0 * g_yz_0_0_z[i] * a_exp + 4.0 * g_yz_zz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (162-165)

    #pragma omp simd aligned(g_xz_0_0_x, g_xz_0_0_y, g_xz_0_0_z, g_xz_xx_0_x, g_xz_xx_0_y, g_xz_xx_0_z, g_z_x_0_0_x_x_0_x, g_z_x_0_0_x_x_0_y, g_z_x_0_0_x_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_x_0_x[i] = -2.0 * g_xz_0_0_x[i] * a_exp + 4.0 * g_xz_xx_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_0_y[i] = -2.0 * g_xz_0_0_y[i] * a_exp + 4.0 * g_xz_xx_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_x_0_z[i] = -2.0 * g_xz_0_0_z[i] * a_exp + 4.0 * g_xz_xx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (165-168)

    #pragma omp simd aligned(g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z, g_z_x_0_0_x_y_0_x, g_z_x_0_0_x_y_0_y, g_z_x_0_0_x_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_y_0_x[i] = 4.0 * g_xz_xy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_0_y[i] = 4.0 * g_xz_xy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_0_z[i] = 4.0 * g_xz_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (168-171)

    #pragma omp simd aligned(g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z, g_z_x_0_0_x_z_0_x, g_z_x_0_0_x_z_0_y, g_z_x_0_0_x_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_z_0_x[i] = 4.0 * g_xz_xz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_0_y[i] = 4.0 * g_xz_xz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_0_z[i] = 4.0 * g_xz_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (171-174)

    #pragma omp simd aligned(g_yz_0_0_x, g_yz_0_0_y, g_yz_0_0_z, g_yz_xx_0_x, g_yz_xx_0_y, g_yz_xx_0_z, g_z_x_0_0_y_x_0_x, g_z_x_0_0_y_x_0_y, g_z_x_0_0_y_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_x_0_x[i] = -2.0 * g_yz_0_0_x[i] * a_exp + 4.0 * g_yz_xx_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_0_y[i] = -2.0 * g_yz_0_0_y[i] * a_exp + 4.0 * g_yz_xx_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_x_0_z[i] = -2.0 * g_yz_0_0_z[i] * a_exp + 4.0 * g_yz_xx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (174-177)

    #pragma omp simd aligned(g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z, g_z_x_0_0_y_y_0_x, g_z_x_0_0_y_y_0_y, g_z_x_0_0_y_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_y_0_x[i] = 4.0 * g_yz_xy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_0_y[i] = 4.0 * g_yz_xy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_0_z[i] = 4.0 * g_yz_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (177-180)

    #pragma omp simd aligned(g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z, g_z_x_0_0_y_z_0_x, g_z_x_0_0_y_z_0_y, g_z_x_0_0_y_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_z_0_x[i] = 4.0 * g_yz_xz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_0_y[i] = 4.0 * g_yz_xz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_0_z[i] = 4.0 * g_yz_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (180-183)

    #pragma omp simd aligned(g_0_0_0_x, g_0_0_0_y, g_0_0_0_z, g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_z_x_0_0_z_x_0_x, g_z_x_0_0_z_x_0_y, g_z_x_0_0_z_x_0_z, g_zz_0_0_x, g_zz_0_0_y, g_zz_0_0_z, g_zz_xx_0_x, g_zz_xx_0_y, g_zz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_x_0_x[i] = g_0_0_0_x[i] - 2.0 * g_0_xx_0_x[i] * b_exp - 2.0 * g_zz_0_0_x[i] * a_exp + 4.0 * g_zz_xx_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_0_y[i] = g_0_0_0_y[i] - 2.0 * g_0_xx_0_y[i] * b_exp - 2.0 * g_zz_0_0_y[i] * a_exp + 4.0 * g_zz_xx_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_x_0_z[i] = g_0_0_0_z[i] - 2.0 * g_0_xx_0_z[i] * b_exp - 2.0 * g_zz_0_0_z[i] * a_exp + 4.0 * g_zz_xx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (183-186)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_z_x_0_0_z_y_0_x, g_z_x_0_0_z_y_0_y, g_z_x_0_0_z_y_0_z, g_zz_xy_0_x, g_zz_xy_0_y, g_zz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_y_0_x[i] = -2.0 * g_0_xy_0_x[i] * b_exp + 4.0 * g_zz_xy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_0_y[i] = -2.0 * g_0_xy_0_y[i] * b_exp + 4.0 * g_zz_xy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_0_z[i] = -2.0 * g_0_xy_0_z[i] * b_exp + 4.0 * g_zz_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (186-189)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_z_x_0_0_z_z_0_x, g_z_x_0_0_z_z_0_y, g_z_x_0_0_z_z_0_z, g_zz_xz_0_x, g_zz_xz_0_y, g_zz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_z_0_x[i] = -2.0 * g_0_xz_0_x[i] * b_exp + 4.0 * g_zz_xz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_0_y[i] = -2.0 * g_0_xz_0_y[i] * b_exp + 4.0 * g_zz_xz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_0_z[i] = -2.0 * g_0_xz_0_z[i] * b_exp + 4.0 * g_zz_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (189-192)

    #pragma omp simd aligned(g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z, g_z_y_0_0_x_x_0_x, g_z_y_0_0_x_x_0_y, g_z_y_0_0_x_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_x_0_x[i] = 4.0 * g_xz_xy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_0_y[i] = 4.0 * g_xz_xy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_x_0_z[i] = 4.0 * g_xz_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (192-195)

    #pragma omp simd aligned(g_xz_0_0_x, g_xz_0_0_y, g_xz_0_0_z, g_xz_yy_0_x, g_xz_yy_0_y, g_xz_yy_0_z, g_z_y_0_0_x_y_0_x, g_z_y_0_0_x_y_0_y, g_z_y_0_0_x_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_y_0_x[i] = -2.0 * g_xz_0_0_x[i] * a_exp + 4.0 * g_xz_yy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_0_y[i] = -2.0 * g_xz_0_0_y[i] * a_exp + 4.0 * g_xz_yy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_0_z[i] = -2.0 * g_xz_0_0_z[i] * a_exp + 4.0 * g_xz_yy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (195-198)

    #pragma omp simd aligned(g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z, g_z_y_0_0_x_z_0_x, g_z_y_0_0_x_z_0_y, g_z_y_0_0_x_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_z_0_x[i] = 4.0 * g_xz_yz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_0_y[i] = 4.0 * g_xz_yz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_0_z[i] = 4.0 * g_xz_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (198-201)

    #pragma omp simd aligned(g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z, g_z_y_0_0_y_x_0_x, g_z_y_0_0_y_x_0_y, g_z_y_0_0_y_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_x_0_x[i] = 4.0 * g_yz_xy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_0_y[i] = 4.0 * g_yz_xy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_x_0_z[i] = 4.0 * g_yz_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (201-204)

    #pragma omp simd aligned(g_yz_0_0_x, g_yz_0_0_y, g_yz_0_0_z, g_yz_yy_0_x, g_yz_yy_0_y, g_yz_yy_0_z, g_z_y_0_0_y_y_0_x, g_z_y_0_0_y_y_0_y, g_z_y_0_0_y_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_y_0_x[i] = -2.0 * g_yz_0_0_x[i] * a_exp + 4.0 * g_yz_yy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_0_y[i] = -2.0 * g_yz_0_0_y[i] * a_exp + 4.0 * g_yz_yy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_0_z[i] = -2.0 * g_yz_0_0_z[i] * a_exp + 4.0 * g_yz_yy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (204-207)

    #pragma omp simd aligned(g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z, g_z_y_0_0_y_z_0_x, g_z_y_0_0_y_z_0_y, g_z_y_0_0_y_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_z_0_x[i] = 4.0 * g_yz_yz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_0_y[i] = 4.0 * g_yz_yz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_0_z[i] = 4.0 * g_yz_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (207-210)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_z_y_0_0_z_x_0_x, g_z_y_0_0_z_x_0_y, g_z_y_0_0_z_x_0_z, g_zz_xy_0_x, g_zz_xy_0_y, g_zz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_x_0_x[i] = -2.0 * g_0_xy_0_x[i] * b_exp + 4.0 * g_zz_xy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_0_y[i] = -2.0 * g_0_xy_0_y[i] * b_exp + 4.0 * g_zz_xy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_x_0_z[i] = -2.0 * g_0_xy_0_z[i] * b_exp + 4.0 * g_zz_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (210-213)

    #pragma omp simd aligned(g_0_0_0_x, g_0_0_0_y, g_0_0_0_z, g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_z_y_0_0_z_y_0_x, g_z_y_0_0_z_y_0_y, g_z_y_0_0_z_y_0_z, g_zz_0_0_x, g_zz_0_0_y, g_zz_0_0_z, g_zz_yy_0_x, g_zz_yy_0_y, g_zz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_y_0_x[i] = g_0_0_0_x[i] - 2.0 * g_0_yy_0_x[i] * b_exp - 2.0 * g_zz_0_0_x[i] * a_exp + 4.0 * g_zz_yy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_0_y[i] = g_0_0_0_y[i] - 2.0 * g_0_yy_0_y[i] * b_exp - 2.0 * g_zz_0_0_y[i] * a_exp + 4.0 * g_zz_yy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_0_z[i] = g_0_0_0_z[i] - 2.0 * g_0_yy_0_z[i] * b_exp - 2.0 * g_zz_0_0_z[i] * a_exp + 4.0 * g_zz_yy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (213-216)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_z_y_0_0_z_z_0_x, g_z_y_0_0_z_z_0_y, g_z_y_0_0_z_z_0_z, g_zz_yz_0_x, g_zz_yz_0_y, g_zz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_z_0_x[i] = -2.0 * g_0_yz_0_x[i] * b_exp + 4.0 * g_zz_yz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_0_y[i] = -2.0 * g_0_yz_0_y[i] * b_exp + 4.0 * g_zz_yz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_0_z[i] = -2.0 * g_0_yz_0_z[i] * b_exp + 4.0 * g_zz_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (216-219)

    #pragma omp simd aligned(g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z, g_z_z_0_0_x_x_0_x, g_z_z_0_0_x_x_0_y, g_z_z_0_0_x_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_x_0_x[i] = 4.0 * g_xz_xz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_0_y[i] = 4.0 * g_xz_xz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_x_0_z[i] = 4.0 * g_xz_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (219-222)

    #pragma omp simd aligned(g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z, g_z_z_0_0_x_y_0_x, g_z_z_0_0_x_y_0_y, g_z_z_0_0_x_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_y_0_x[i] = 4.0 * g_xz_yz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_0_y[i] = 4.0 * g_xz_yz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_0_z[i] = 4.0 * g_xz_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (222-225)

    #pragma omp simd aligned(g_xz_0_0_x, g_xz_0_0_y, g_xz_0_0_z, g_xz_zz_0_x, g_xz_zz_0_y, g_xz_zz_0_z, g_z_z_0_0_x_z_0_x, g_z_z_0_0_x_z_0_y, g_z_z_0_0_x_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_z_0_x[i] = -2.0 * g_xz_0_0_x[i] * a_exp + 4.0 * g_xz_zz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_0_y[i] = -2.0 * g_xz_0_0_y[i] * a_exp + 4.0 * g_xz_zz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_0_z[i] = -2.0 * g_xz_0_0_z[i] * a_exp + 4.0 * g_xz_zz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (225-228)

    #pragma omp simd aligned(g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z, g_z_z_0_0_y_x_0_x, g_z_z_0_0_y_x_0_y, g_z_z_0_0_y_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_x_0_x[i] = 4.0 * g_yz_xz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_0_y[i] = 4.0 * g_yz_xz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_x_0_z[i] = 4.0 * g_yz_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (228-231)

    #pragma omp simd aligned(g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z, g_z_z_0_0_y_y_0_x, g_z_z_0_0_y_y_0_y, g_z_z_0_0_y_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_y_0_x[i] = 4.0 * g_yz_yz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_0_y[i] = 4.0 * g_yz_yz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_0_z[i] = 4.0 * g_yz_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (231-234)

    #pragma omp simd aligned(g_yz_0_0_x, g_yz_0_0_y, g_yz_0_0_z, g_yz_zz_0_x, g_yz_zz_0_y, g_yz_zz_0_z, g_z_z_0_0_y_z_0_x, g_z_z_0_0_y_z_0_y, g_z_z_0_0_y_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_z_0_x[i] = -2.0 * g_yz_0_0_x[i] * a_exp + 4.0 * g_yz_zz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_0_y[i] = -2.0 * g_yz_0_0_y[i] * a_exp + 4.0 * g_yz_zz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_0_z[i] = -2.0 * g_yz_0_0_z[i] * a_exp + 4.0 * g_yz_zz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (234-237)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_z_z_0_0_z_x_0_x, g_z_z_0_0_z_x_0_y, g_z_z_0_0_z_x_0_z, g_zz_xz_0_x, g_zz_xz_0_y, g_zz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_x_0_x[i] = -2.0 * g_0_xz_0_x[i] * b_exp + 4.0 * g_zz_xz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_0_y[i] = -2.0 * g_0_xz_0_y[i] * b_exp + 4.0 * g_zz_xz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_x_0_z[i] = -2.0 * g_0_xz_0_z[i] * b_exp + 4.0 * g_zz_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (237-240)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_z_z_0_0_z_y_0_x, g_z_z_0_0_z_y_0_y, g_z_z_0_0_z_y_0_z, g_zz_yz_0_x, g_zz_yz_0_y, g_zz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_y_0_x[i] = -2.0 * g_0_yz_0_x[i] * b_exp + 4.0 * g_zz_yz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_0_y[i] = -2.0 * g_0_yz_0_y[i] * b_exp + 4.0 * g_zz_yz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_0_z[i] = -2.0 * g_0_yz_0_z[i] * b_exp + 4.0 * g_zz_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (240-243)

    #pragma omp simd aligned(g_0_0_0_x, g_0_0_0_y, g_0_0_0_z, g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_z_z_0_0_z_z_0_x, g_z_z_0_0_z_z_0_y, g_z_z_0_0_z_z_0_z, g_zz_0_0_x, g_zz_0_0_y, g_zz_0_0_z, g_zz_zz_0_x, g_zz_zz_0_y, g_zz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_z_0_x[i] = g_0_0_0_x[i] - 2.0 * g_0_zz_0_x[i] * b_exp - 2.0 * g_zz_0_0_x[i] * a_exp + 4.0 * g_zz_zz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_0_y[i] = g_0_0_0_y[i] - 2.0 * g_0_zz_0_y[i] * b_exp - 2.0 * g_zz_0_0_y[i] * a_exp + 4.0 * g_zz_zz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_0_z[i] = g_0_0_0_z[i] - 2.0 * g_0_zz_0_z[i] * b_exp - 2.0 * g_zz_0_0_z[i] * a_exp + 4.0 * g_zz_zz_0_z[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

