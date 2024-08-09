#include "GeomDeriv1000OfScalarForPDSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_pdsp_0(CSimdArray<double>& buffer_1000_pdsp,
                     const CSimdArray<double>& buffer_sdsp,
                     const CSimdArray<double>& buffer_ddsp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_pdsp.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_pdsp

    auto g_x_0_0_0_x_xx_0_x = buffer_1000_pdsp[0];

    auto g_x_0_0_0_x_xx_0_y = buffer_1000_pdsp[1];

    auto g_x_0_0_0_x_xx_0_z = buffer_1000_pdsp[2];

    auto g_x_0_0_0_x_xy_0_x = buffer_1000_pdsp[3];

    auto g_x_0_0_0_x_xy_0_y = buffer_1000_pdsp[4];

    auto g_x_0_0_0_x_xy_0_z = buffer_1000_pdsp[5];

    auto g_x_0_0_0_x_xz_0_x = buffer_1000_pdsp[6];

    auto g_x_0_0_0_x_xz_0_y = buffer_1000_pdsp[7];

    auto g_x_0_0_0_x_xz_0_z = buffer_1000_pdsp[8];

    auto g_x_0_0_0_x_yy_0_x = buffer_1000_pdsp[9];

    auto g_x_0_0_0_x_yy_0_y = buffer_1000_pdsp[10];

    auto g_x_0_0_0_x_yy_0_z = buffer_1000_pdsp[11];

    auto g_x_0_0_0_x_yz_0_x = buffer_1000_pdsp[12];

    auto g_x_0_0_0_x_yz_0_y = buffer_1000_pdsp[13];

    auto g_x_0_0_0_x_yz_0_z = buffer_1000_pdsp[14];

    auto g_x_0_0_0_x_zz_0_x = buffer_1000_pdsp[15];

    auto g_x_0_0_0_x_zz_0_y = buffer_1000_pdsp[16];

    auto g_x_0_0_0_x_zz_0_z = buffer_1000_pdsp[17];

    auto g_x_0_0_0_y_xx_0_x = buffer_1000_pdsp[18];

    auto g_x_0_0_0_y_xx_0_y = buffer_1000_pdsp[19];

    auto g_x_0_0_0_y_xx_0_z = buffer_1000_pdsp[20];

    auto g_x_0_0_0_y_xy_0_x = buffer_1000_pdsp[21];

    auto g_x_0_0_0_y_xy_0_y = buffer_1000_pdsp[22];

    auto g_x_0_0_0_y_xy_0_z = buffer_1000_pdsp[23];

    auto g_x_0_0_0_y_xz_0_x = buffer_1000_pdsp[24];

    auto g_x_0_0_0_y_xz_0_y = buffer_1000_pdsp[25];

    auto g_x_0_0_0_y_xz_0_z = buffer_1000_pdsp[26];

    auto g_x_0_0_0_y_yy_0_x = buffer_1000_pdsp[27];

    auto g_x_0_0_0_y_yy_0_y = buffer_1000_pdsp[28];

    auto g_x_0_0_0_y_yy_0_z = buffer_1000_pdsp[29];

    auto g_x_0_0_0_y_yz_0_x = buffer_1000_pdsp[30];

    auto g_x_0_0_0_y_yz_0_y = buffer_1000_pdsp[31];

    auto g_x_0_0_0_y_yz_0_z = buffer_1000_pdsp[32];

    auto g_x_0_0_0_y_zz_0_x = buffer_1000_pdsp[33];

    auto g_x_0_0_0_y_zz_0_y = buffer_1000_pdsp[34];

    auto g_x_0_0_0_y_zz_0_z = buffer_1000_pdsp[35];

    auto g_x_0_0_0_z_xx_0_x = buffer_1000_pdsp[36];

    auto g_x_0_0_0_z_xx_0_y = buffer_1000_pdsp[37];

    auto g_x_0_0_0_z_xx_0_z = buffer_1000_pdsp[38];

    auto g_x_0_0_0_z_xy_0_x = buffer_1000_pdsp[39];

    auto g_x_0_0_0_z_xy_0_y = buffer_1000_pdsp[40];

    auto g_x_0_0_0_z_xy_0_z = buffer_1000_pdsp[41];

    auto g_x_0_0_0_z_xz_0_x = buffer_1000_pdsp[42];

    auto g_x_0_0_0_z_xz_0_y = buffer_1000_pdsp[43];

    auto g_x_0_0_0_z_xz_0_z = buffer_1000_pdsp[44];

    auto g_x_0_0_0_z_yy_0_x = buffer_1000_pdsp[45];

    auto g_x_0_0_0_z_yy_0_y = buffer_1000_pdsp[46];

    auto g_x_0_0_0_z_yy_0_z = buffer_1000_pdsp[47];

    auto g_x_0_0_0_z_yz_0_x = buffer_1000_pdsp[48];

    auto g_x_0_0_0_z_yz_0_y = buffer_1000_pdsp[49];

    auto g_x_0_0_0_z_yz_0_z = buffer_1000_pdsp[50];

    auto g_x_0_0_0_z_zz_0_x = buffer_1000_pdsp[51];

    auto g_x_0_0_0_z_zz_0_y = buffer_1000_pdsp[52];

    auto g_x_0_0_0_z_zz_0_z = buffer_1000_pdsp[53];

    auto g_y_0_0_0_x_xx_0_x = buffer_1000_pdsp[54];

    auto g_y_0_0_0_x_xx_0_y = buffer_1000_pdsp[55];

    auto g_y_0_0_0_x_xx_0_z = buffer_1000_pdsp[56];

    auto g_y_0_0_0_x_xy_0_x = buffer_1000_pdsp[57];

    auto g_y_0_0_0_x_xy_0_y = buffer_1000_pdsp[58];

    auto g_y_0_0_0_x_xy_0_z = buffer_1000_pdsp[59];

    auto g_y_0_0_0_x_xz_0_x = buffer_1000_pdsp[60];

    auto g_y_0_0_0_x_xz_0_y = buffer_1000_pdsp[61];

    auto g_y_0_0_0_x_xz_0_z = buffer_1000_pdsp[62];

    auto g_y_0_0_0_x_yy_0_x = buffer_1000_pdsp[63];

    auto g_y_0_0_0_x_yy_0_y = buffer_1000_pdsp[64];

    auto g_y_0_0_0_x_yy_0_z = buffer_1000_pdsp[65];

    auto g_y_0_0_0_x_yz_0_x = buffer_1000_pdsp[66];

    auto g_y_0_0_0_x_yz_0_y = buffer_1000_pdsp[67];

    auto g_y_0_0_0_x_yz_0_z = buffer_1000_pdsp[68];

    auto g_y_0_0_0_x_zz_0_x = buffer_1000_pdsp[69];

    auto g_y_0_0_0_x_zz_0_y = buffer_1000_pdsp[70];

    auto g_y_0_0_0_x_zz_0_z = buffer_1000_pdsp[71];

    auto g_y_0_0_0_y_xx_0_x = buffer_1000_pdsp[72];

    auto g_y_0_0_0_y_xx_0_y = buffer_1000_pdsp[73];

    auto g_y_0_0_0_y_xx_0_z = buffer_1000_pdsp[74];

    auto g_y_0_0_0_y_xy_0_x = buffer_1000_pdsp[75];

    auto g_y_0_0_0_y_xy_0_y = buffer_1000_pdsp[76];

    auto g_y_0_0_0_y_xy_0_z = buffer_1000_pdsp[77];

    auto g_y_0_0_0_y_xz_0_x = buffer_1000_pdsp[78];

    auto g_y_0_0_0_y_xz_0_y = buffer_1000_pdsp[79];

    auto g_y_0_0_0_y_xz_0_z = buffer_1000_pdsp[80];

    auto g_y_0_0_0_y_yy_0_x = buffer_1000_pdsp[81];

    auto g_y_0_0_0_y_yy_0_y = buffer_1000_pdsp[82];

    auto g_y_0_0_0_y_yy_0_z = buffer_1000_pdsp[83];

    auto g_y_0_0_0_y_yz_0_x = buffer_1000_pdsp[84];

    auto g_y_0_0_0_y_yz_0_y = buffer_1000_pdsp[85];

    auto g_y_0_0_0_y_yz_0_z = buffer_1000_pdsp[86];

    auto g_y_0_0_0_y_zz_0_x = buffer_1000_pdsp[87];

    auto g_y_0_0_0_y_zz_0_y = buffer_1000_pdsp[88];

    auto g_y_0_0_0_y_zz_0_z = buffer_1000_pdsp[89];

    auto g_y_0_0_0_z_xx_0_x = buffer_1000_pdsp[90];

    auto g_y_0_0_0_z_xx_0_y = buffer_1000_pdsp[91];

    auto g_y_0_0_0_z_xx_0_z = buffer_1000_pdsp[92];

    auto g_y_0_0_0_z_xy_0_x = buffer_1000_pdsp[93];

    auto g_y_0_0_0_z_xy_0_y = buffer_1000_pdsp[94];

    auto g_y_0_0_0_z_xy_0_z = buffer_1000_pdsp[95];

    auto g_y_0_0_0_z_xz_0_x = buffer_1000_pdsp[96];

    auto g_y_0_0_0_z_xz_0_y = buffer_1000_pdsp[97];

    auto g_y_0_0_0_z_xz_0_z = buffer_1000_pdsp[98];

    auto g_y_0_0_0_z_yy_0_x = buffer_1000_pdsp[99];

    auto g_y_0_0_0_z_yy_0_y = buffer_1000_pdsp[100];

    auto g_y_0_0_0_z_yy_0_z = buffer_1000_pdsp[101];

    auto g_y_0_0_0_z_yz_0_x = buffer_1000_pdsp[102];

    auto g_y_0_0_0_z_yz_0_y = buffer_1000_pdsp[103];

    auto g_y_0_0_0_z_yz_0_z = buffer_1000_pdsp[104];

    auto g_y_0_0_0_z_zz_0_x = buffer_1000_pdsp[105];

    auto g_y_0_0_0_z_zz_0_y = buffer_1000_pdsp[106];

    auto g_y_0_0_0_z_zz_0_z = buffer_1000_pdsp[107];

    auto g_z_0_0_0_x_xx_0_x = buffer_1000_pdsp[108];

    auto g_z_0_0_0_x_xx_0_y = buffer_1000_pdsp[109];

    auto g_z_0_0_0_x_xx_0_z = buffer_1000_pdsp[110];

    auto g_z_0_0_0_x_xy_0_x = buffer_1000_pdsp[111];

    auto g_z_0_0_0_x_xy_0_y = buffer_1000_pdsp[112];

    auto g_z_0_0_0_x_xy_0_z = buffer_1000_pdsp[113];

    auto g_z_0_0_0_x_xz_0_x = buffer_1000_pdsp[114];

    auto g_z_0_0_0_x_xz_0_y = buffer_1000_pdsp[115];

    auto g_z_0_0_0_x_xz_0_z = buffer_1000_pdsp[116];

    auto g_z_0_0_0_x_yy_0_x = buffer_1000_pdsp[117];

    auto g_z_0_0_0_x_yy_0_y = buffer_1000_pdsp[118];

    auto g_z_0_0_0_x_yy_0_z = buffer_1000_pdsp[119];

    auto g_z_0_0_0_x_yz_0_x = buffer_1000_pdsp[120];

    auto g_z_0_0_0_x_yz_0_y = buffer_1000_pdsp[121];

    auto g_z_0_0_0_x_yz_0_z = buffer_1000_pdsp[122];

    auto g_z_0_0_0_x_zz_0_x = buffer_1000_pdsp[123];

    auto g_z_0_0_0_x_zz_0_y = buffer_1000_pdsp[124];

    auto g_z_0_0_0_x_zz_0_z = buffer_1000_pdsp[125];

    auto g_z_0_0_0_y_xx_0_x = buffer_1000_pdsp[126];

    auto g_z_0_0_0_y_xx_0_y = buffer_1000_pdsp[127];

    auto g_z_0_0_0_y_xx_0_z = buffer_1000_pdsp[128];

    auto g_z_0_0_0_y_xy_0_x = buffer_1000_pdsp[129];

    auto g_z_0_0_0_y_xy_0_y = buffer_1000_pdsp[130];

    auto g_z_0_0_0_y_xy_0_z = buffer_1000_pdsp[131];

    auto g_z_0_0_0_y_xz_0_x = buffer_1000_pdsp[132];

    auto g_z_0_0_0_y_xz_0_y = buffer_1000_pdsp[133];

    auto g_z_0_0_0_y_xz_0_z = buffer_1000_pdsp[134];

    auto g_z_0_0_0_y_yy_0_x = buffer_1000_pdsp[135];

    auto g_z_0_0_0_y_yy_0_y = buffer_1000_pdsp[136];

    auto g_z_0_0_0_y_yy_0_z = buffer_1000_pdsp[137];

    auto g_z_0_0_0_y_yz_0_x = buffer_1000_pdsp[138];

    auto g_z_0_0_0_y_yz_0_y = buffer_1000_pdsp[139];

    auto g_z_0_0_0_y_yz_0_z = buffer_1000_pdsp[140];

    auto g_z_0_0_0_y_zz_0_x = buffer_1000_pdsp[141];

    auto g_z_0_0_0_y_zz_0_y = buffer_1000_pdsp[142];

    auto g_z_0_0_0_y_zz_0_z = buffer_1000_pdsp[143];

    auto g_z_0_0_0_z_xx_0_x = buffer_1000_pdsp[144];

    auto g_z_0_0_0_z_xx_0_y = buffer_1000_pdsp[145];

    auto g_z_0_0_0_z_xx_0_z = buffer_1000_pdsp[146];

    auto g_z_0_0_0_z_xy_0_x = buffer_1000_pdsp[147];

    auto g_z_0_0_0_z_xy_0_y = buffer_1000_pdsp[148];

    auto g_z_0_0_0_z_xy_0_z = buffer_1000_pdsp[149];

    auto g_z_0_0_0_z_xz_0_x = buffer_1000_pdsp[150];

    auto g_z_0_0_0_z_xz_0_y = buffer_1000_pdsp[151];

    auto g_z_0_0_0_z_xz_0_z = buffer_1000_pdsp[152];

    auto g_z_0_0_0_z_yy_0_x = buffer_1000_pdsp[153];

    auto g_z_0_0_0_z_yy_0_y = buffer_1000_pdsp[154];

    auto g_z_0_0_0_z_yy_0_z = buffer_1000_pdsp[155];

    auto g_z_0_0_0_z_yz_0_x = buffer_1000_pdsp[156];

    auto g_z_0_0_0_z_yz_0_y = buffer_1000_pdsp[157];

    auto g_z_0_0_0_z_yz_0_z = buffer_1000_pdsp[158];

    auto g_z_0_0_0_z_zz_0_x = buffer_1000_pdsp[159];

    auto g_z_0_0_0_z_zz_0_y = buffer_1000_pdsp[160];

    auto g_z_0_0_0_z_zz_0_z = buffer_1000_pdsp[161];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_x_0_0_0_x_xx_0_x, g_x_0_0_0_x_xx_0_y, g_x_0_0_0_x_xx_0_z, g_xx_xx_0_x, g_xx_xx_0_y, g_xx_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xx_0_x[i] = -g_0_xx_0_x[i] + 2.0 * g_xx_xx_0_x[i] * a_exp;

        g_x_0_0_0_x_xx_0_y[i] = -g_0_xx_0_y[i] + 2.0 * g_xx_xx_0_y[i] * a_exp;

        g_x_0_0_0_x_xx_0_z[i] = -g_0_xx_0_z[i] + 2.0 * g_xx_xx_0_z[i] * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_x_0_0_0_x_xy_0_x, g_x_0_0_0_x_xy_0_y, g_x_0_0_0_x_xy_0_z, g_xx_xy_0_x, g_xx_xy_0_y, g_xx_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xy_0_x[i] = -g_0_xy_0_x[i] + 2.0 * g_xx_xy_0_x[i] * a_exp;

        g_x_0_0_0_x_xy_0_y[i] = -g_0_xy_0_y[i] + 2.0 * g_xx_xy_0_y[i] * a_exp;

        g_x_0_0_0_x_xy_0_z[i] = -g_0_xy_0_z[i] + 2.0 * g_xx_xy_0_z[i] * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_x_0_0_0_x_xz_0_x, g_x_0_0_0_x_xz_0_y, g_x_0_0_0_x_xz_0_z, g_xx_xz_0_x, g_xx_xz_0_y, g_xx_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_xz_0_x[i] = -g_0_xz_0_x[i] + 2.0 * g_xx_xz_0_x[i] * a_exp;

        g_x_0_0_0_x_xz_0_y[i] = -g_0_xz_0_y[i] + 2.0 * g_xx_xz_0_y[i] * a_exp;

        g_x_0_0_0_x_xz_0_z[i] = -g_0_xz_0_z[i] + 2.0 * g_xx_xz_0_z[i] * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_x_0_0_0_x_yy_0_x, g_x_0_0_0_x_yy_0_y, g_x_0_0_0_x_yy_0_z, g_xx_yy_0_x, g_xx_yy_0_y, g_xx_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yy_0_x[i] = -g_0_yy_0_x[i] + 2.0 * g_xx_yy_0_x[i] * a_exp;

        g_x_0_0_0_x_yy_0_y[i] = -g_0_yy_0_y[i] + 2.0 * g_xx_yy_0_y[i] * a_exp;

        g_x_0_0_0_x_yy_0_z[i] = -g_0_yy_0_z[i] + 2.0 * g_xx_yy_0_z[i] * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_x_0_0_0_x_yz_0_x, g_x_0_0_0_x_yz_0_y, g_x_0_0_0_x_yz_0_z, g_xx_yz_0_x, g_xx_yz_0_y, g_xx_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_yz_0_x[i] = -g_0_yz_0_x[i] + 2.0 * g_xx_yz_0_x[i] * a_exp;

        g_x_0_0_0_x_yz_0_y[i] = -g_0_yz_0_y[i] + 2.0 * g_xx_yz_0_y[i] * a_exp;

        g_x_0_0_0_x_yz_0_z[i] = -g_0_yz_0_z[i] + 2.0 * g_xx_yz_0_z[i] * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_x_0_0_0_x_zz_0_x, g_x_0_0_0_x_zz_0_y, g_x_0_0_0_x_zz_0_z, g_xx_zz_0_x, g_xx_zz_0_y, g_xx_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_zz_0_x[i] = -g_0_zz_0_x[i] + 2.0 * g_xx_zz_0_x[i] * a_exp;

        g_x_0_0_0_x_zz_0_y[i] = -g_0_zz_0_y[i] + 2.0 * g_xx_zz_0_y[i] * a_exp;

        g_x_0_0_0_x_zz_0_z[i] = -g_0_zz_0_z[i] + 2.0 * g_xx_zz_0_z[i] * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_0_0_0_y_xx_0_x, g_x_0_0_0_y_xx_0_y, g_x_0_0_0_y_xx_0_z, g_xy_xx_0_x, g_xy_xx_0_y, g_xy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xx_0_x[i] = 2.0 * g_xy_xx_0_x[i] * a_exp;

        g_x_0_0_0_y_xx_0_y[i] = 2.0 * g_xy_xx_0_y[i] * a_exp;

        g_x_0_0_0_y_xx_0_z[i] = 2.0 * g_xy_xx_0_z[i] * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_0_0_0_y_xy_0_x, g_x_0_0_0_y_xy_0_y, g_x_0_0_0_y_xy_0_z, g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xy_0_x[i] = 2.0 * g_xy_xy_0_x[i] * a_exp;

        g_x_0_0_0_y_xy_0_y[i] = 2.0 * g_xy_xy_0_y[i] * a_exp;

        g_x_0_0_0_y_xy_0_z[i] = 2.0 * g_xy_xy_0_z[i] * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_0_0_0_y_xz_0_x, g_x_0_0_0_y_xz_0_y, g_x_0_0_0_y_xz_0_z, g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_xz_0_x[i] = 2.0 * g_xy_xz_0_x[i] * a_exp;

        g_x_0_0_0_y_xz_0_y[i] = 2.0 * g_xy_xz_0_y[i] * a_exp;

        g_x_0_0_0_y_xz_0_z[i] = 2.0 * g_xy_xz_0_z[i] * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_0_0_0_y_yy_0_x, g_x_0_0_0_y_yy_0_y, g_x_0_0_0_y_yy_0_z, g_xy_yy_0_x, g_xy_yy_0_y, g_xy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yy_0_x[i] = 2.0 * g_xy_yy_0_x[i] * a_exp;

        g_x_0_0_0_y_yy_0_y[i] = 2.0 * g_xy_yy_0_y[i] * a_exp;

        g_x_0_0_0_y_yy_0_z[i] = 2.0 * g_xy_yy_0_z[i] * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_x_0_0_0_y_yz_0_x, g_x_0_0_0_y_yz_0_y, g_x_0_0_0_y_yz_0_z, g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_yz_0_x[i] = 2.0 * g_xy_yz_0_x[i] * a_exp;

        g_x_0_0_0_y_yz_0_y[i] = 2.0 * g_xy_yz_0_y[i] * a_exp;

        g_x_0_0_0_y_yz_0_z[i] = 2.0 * g_xy_yz_0_z[i] * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_x_0_0_0_y_zz_0_x, g_x_0_0_0_y_zz_0_y, g_x_0_0_0_y_zz_0_z, g_xy_zz_0_x, g_xy_zz_0_y, g_xy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_zz_0_x[i] = 2.0 * g_xy_zz_0_x[i] * a_exp;

        g_x_0_0_0_y_zz_0_y[i] = 2.0 * g_xy_zz_0_y[i] * a_exp;

        g_x_0_0_0_y_zz_0_z[i] = 2.0 * g_xy_zz_0_z[i] * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_0_0_0_z_xx_0_x, g_x_0_0_0_z_xx_0_y, g_x_0_0_0_z_xx_0_z, g_xz_xx_0_x, g_xz_xx_0_y, g_xz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xx_0_x[i] = 2.0 * g_xz_xx_0_x[i] * a_exp;

        g_x_0_0_0_z_xx_0_y[i] = 2.0 * g_xz_xx_0_y[i] * a_exp;

        g_x_0_0_0_z_xx_0_z[i] = 2.0 * g_xz_xx_0_z[i] * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_0_0_0_z_xy_0_x, g_x_0_0_0_z_xy_0_y, g_x_0_0_0_z_xy_0_z, g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xy_0_x[i] = 2.0 * g_xz_xy_0_x[i] * a_exp;

        g_x_0_0_0_z_xy_0_y[i] = 2.0 * g_xz_xy_0_y[i] * a_exp;

        g_x_0_0_0_z_xy_0_z[i] = 2.0 * g_xz_xy_0_z[i] * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_0_0_0_z_xz_0_x, g_x_0_0_0_z_xz_0_y, g_x_0_0_0_z_xz_0_z, g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_xz_0_x[i] = 2.0 * g_xz_xz_0_x[i] * a_exp;

        g_x_0_0_0_z_xz_0_y[i] = 2.0 * g_xz_xz_0_y[i] * a_exp;

        g_x_0_0_0_z_xz_0_z[i] = 2.0 * g_xz_xz_0_z[i] * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_0_0_0_z_yy_0_x, g_x_0_0_0_z_yy_0_y, g_x_0_0_0_z_yy_0_z, g_xz_yy_0_x, g_xz_yy_0_y, g_xz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yy_0_x[i] = 2.0 * g_xz_yy_0_x[i] * a_exp;

        g_x_0_0_0_z_yy_0_y[i] = 2.0 * g_xz_yy_0_y[i] * a_exp;

        g_x_0_0_0_z_yy_0_z[i] = 2.0 * g_xz_yy_0_z[i] * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_0_0_0_z_yz_0_x, g_x_0_0_0_z_yz_0_y, g_x_0_0_0_z_yz_0_z, g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_yz_0_x[i] = 2.0 * g_xz_yz_0_x[i] * a_exp;

        g_x_0_0_0_z_yz_0_y[i] = 2.0 * g_xz_yz_0_y[i] * a_exp;

        g_x_0_0_0_z_yz_0_z[i] = 2.0 * g_xz_yz_0_z[i] * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_0_0_0_z_zz_0_x, g_x_0_0_0_z_zz_0_y, g_x_0_0_0_z_zz_0_z, g_xz_zz_0_x, g_xz_zz_0_y, g_xz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_zz_0_x[i] = 2.0 * g_xz_zz_0_x[i] * a_exp;

        g_x_0_0_0_z_zz_0_y[i] = 2.0 * g_xz_zz_0_y[i] * a_exp;

        g_x_0_0_0_z_zz_0_z[i] = 2.0 * g_xz_zz_0_z[i] * a_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_xy_xx_0_x, g_xy_xx_0_y, g_xy_xx_0_z, g_y_0_0_0_x_xx_0_x, g_y_0_0_0_x_xx_0_y, g_y_0_0_0_x_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xx_0_x[i] = 2.0 * g_xy_xx_0_x[i] * a_exp;

        g_y_0_0_0_x_xx_0_y[i] = 2.0 * g_xy_xx_0_y[i] * a_exp;

        g_y_0_0_0_x_xx_0_z[i] = 2.0 * g_xy_xx_0_z[i] * a_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_xy_xy_0_x, g_xy_xy_0_y, g_xy_xy_0_z, g_y_0_0_0_x_xy_0_x, g_y_0_0_0_x_xy_0_y, g_y_0_0_0_x_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xy_0_x[i] = 2.0 * g_xy_xy_0_x[i] * a_exp;

        g_y_0_0_0_x_xy_0_y[i] = 2.0 * g_xy_xy_0_y[i] * a_exp;

        g_y_0_0_0_x_xy_0_z[i] = 2.0 * g_xy_xy_0_z[i] * a_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_xy_xz_0_x, g_xy_xz_0_y, g_xy_xz_0_z, g_y_0_0_0_x_xz_0_x, g_y_0_0_0_x_xz_0_y, g_y_0_0_0_x_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_xz_0_x[i] = 2.0 * g_xy_xz_0_x[i] * a_exp;

        g_y_0_0_0_x_xz_0_y[i] = 2.0 * g_xy_xz_0_y[i] * a_exp;

        g_y_0_0_0_x_xz_0_z[i] = 2.0 * g_xy_xz_0_z[i] * a_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_xy_yy_0_x, g_xy_yy_0_y, g_xy_yy_0_z, g_y_0_0_0_x_yy_0_x, g_y_0_0_0_x_yy_0_y, g_y_0_0_0_x_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yy_0_x[i] = 2.0 * g_xy_yy_0_x[i] * a_exp;

        g_y_0_0_0_x_yy_0_y[i] = 2.0 * g_xy_yy_0_y[i] * a_exp;

        g_y_0_0_0_x_yy_0_z[i] = 2.0 * g_xy_yy_0_z[i] * a_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_xy_yz_0_x, g_xy_yz_0_y, g_xy_yz_0_z, g_y_0_0_0_x_yz_0_x, g_y_0_0_0_x_yz_0_y, g_y_0_0_0_x_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_yz_0_x[i] = 2.0 * g_xy_yz_0_x[i] * a_exp;

        g_y_0_0_0_x_yz_0_y[i] = 2.0 * g_xy_yz_0_y[i] * a_exp;

        g_y_0_0_0_x_yz_0_z[i] = 2.0 * g_xy_yz_0_z[i] * a_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_xy_zz_0_x, g_xy_zz_0_y, g_xy_zz_0_z, g_y_0_0_0_x_zz_0_x, g_y_0_0_0_x_zz_0_y, g_y_0_0_0_x_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_zz_0_x[i] = 2.0 * g_xy_zz_0_x[i] * a_exp;

        g_y_0_0_0_x_zz_0_y[i] = 2.0 * g_xy_zz_0_y[i] * a_exp;

        g_y_0_0_0_x_zz_0_z[i] = 2.0 * g_xy_zz_0_z[i] * a_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_y_0_0_0_y_xx_0_x, g_y_0_0_0_y_xx_0_y, g_y_0_0_0_y_xx_0_z, g_yy_xx_0_x, g_yy_xx_0_y, g_yy_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xx_0_x[i] = -g_0_xx_0_x[i] + 2.0 * g_yy_xx_0_x[i] * a_exp;

        g_y_0_0_0_y_xx_0_y[i] = -g_0_xx_0_y[i] + 2.0 * g_yy_xx_0_y[i] * a_exp;

        g_y_0_0_0_y_xx_0_z[i] = -g_0_xx_0_z[i] + 2.0 * g_yy_xx_0_z[i] * a_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_y_0_0_0_y_xy_0_x, g_y_0_0_0_y_xy_0_y, g_y_0_0_0_y_xy_0_z, g_yy_xy_0_x, g_yy_xy_0_y, g_yy_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xy_0_x[i] = -g_0_xy_0_x[i] + 2.0 * g_yy_xy_0_x[i] * a_exp;

        g_y_0_0_0_y_xy_0_y[i] = -g_0_xy_0_y[i] + 2.0 * g_yy_xy_0_y[i] * a_exp;

        g_y_0_0_0_y_xy_0_z[i] = -g_0_xy_0_z[i] + 2.0 * g_yy_xy_0_z[i] * a_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_y_0_0_0_y_xz_0_x, g_y_0_0_0_y_xz_0_y, g_y_0_0_0_y_xz_0_z, g_yy_xz_0_x, g_yy_xz_0_y, g_yy_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_xz_0_x[i] = -g_0_xz_0_x[i] + 2.0 * g_yy_xz_0_x[i] * a_exp;

        g_y_0_0_0_y_xz_0_y[i] = -g_0_xz_0_y[i] + 2.0 * g_yy_xz_0_y[i] * a_exp;

        g_y_0_0_0_y_xz_0_z[i] = -g_0_xz_0_z[i] + 2.0 * g_yy_xz_0_z[i] * a_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_y_0_0_0_y_yy_0_x, g_y_0_0_0_y_yy_0_y, g_y_0_0_0_y_yy_0_z, g_yy_yy_0_x, g_yy_yy_0_y, g_yy_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yy_0_x[i] = -g_0_yy_0_x[i] + 2.0 * g_yy_yy_0_x[i] * a_exp;

        g_y_0_0_0_y_yy_0_y[i] = -g_0_yy_0_y[i] + 2.0 * g_yy_yy_0_y[i] * a_exp;

        g_y_0_0_0_y_yy_0_z[i] = -g_0_yy_0_z[i] + 2.0 * g_yy_yy_0_z[i] * a_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_y_0_0_0_y_yz_0_x, g_y_0_0_0_y_yz_0_y, g_y_0_0_0_y_yz_0_z, g_yy_yz_0_x, g_yy_yz_0_y, g_yy_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_yz_0_x[i] = -g_0_yz_0_x[i] + 2.0 * g_yy_yz_0_x[i] * a_exp;

        g_y_0_0_0_y_yz_0_y[i] = -g_0_yz_0_y[i] + 2.0 * g_yy_yz_0_y[i] * a_exp;

        g_y_0_0_0_y_yz_0_z[i] = -g_0_yz_0_z[i] + 2.0 * g_yy_yz_0_z[i] * a_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_y_0_0_0_y_zz_0_x, g_y_0_0_0_y_zz_0_y, g_y_0_0_0_y_zz_0_z, g_yy_zz_0_x, g_yy_zz_0_y, g_yy_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_zz_0_x[i] = -g_0_zz_0_x[i] + 2.0 * g_yy_zz_0_x[i] * a_exp;

        g_y_0_0_0_y_zz_0_y[i] = -g_0_zz_0_y[i] + 2.0 * g_yy_zz_0_y[i] * a_exp;

        g_y_0_0_0_y_zz_0_z[i] = -g_0_zz_0_z[i] + 2.0 * g_yy_zz_0_z[i] * a_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_y_0_0_0_z_xx_0_x, g_y_0_0_0_z_xx_0_y, g_y_0_0_0_z_xx_0_z, g_yz_xx_0_x, g_yz_xx_0_y, g_yz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xx_0_x[i] = 2.0 * g_yz_xx_0_x[i] * a_exp;

        g_y_0_0_0_z_xx_0_y[i] = 2.0 * g_yz_xx_0_y[i] * a_exp;

        g_y_0_0_0_z_xx_0_z[i] = 2.0 * g_yz_xx_0_z[i] * a_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_y_0_0_0_z_xy_0_x, g_y_0_0_0_z_xy_0_y, g_y_0_0_0_z_xy_0_z, g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xy_0_x[i] = 2.0 * g_yz_xy_0_x[i] * a_exp;

        g_y_0_0_0_z_xy_0_y[i] = 2.0 * g_yz_xy_0_y[i] * a_exp;

        g_y_0_0_0_z_xy_0_z[i] = 2.0 * g_yz_xy_0_z[i] * a_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_y_0_0_0_z_xz_0_x, g_y_0_0_0_z_xz_0_y, g_y_0_0_0_z_xz_0_z, g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_xz_0_x[i] = 2.0 * g_yz_xz_0_x[i] * a_exp;

        g_y_0_0_0_z_xz_0_y[i] = 2.0 * g_yz_xz_0_y[i] * a_exp;

        g_y_0_0_0_z_xz_0_z[i] = 2.0 * g_yz_xz_0_z[i] * a_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_y_0_0_0_z_yy_0_x, g_y_0_0_0_z_yy_0_y, g_y_0_0_0_z_yy_0_z, g_yz_yy_0_x, g_yz_yy_0_y, g_yz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yy_0_x[i] = 2.0 * g_yz_yy_0_x[i] * a_exp;

        g_y_0_0_0_z_yy_0_y[i] = 2.0 * g_yz_yy_0_y[i] * a_exp;

        g_y_0_0_0_z_yy_0_z[i] = 2.0 * g_yz_yy_0_z[i] * a_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_y_0_0_0_z_yz_0_x, g_y_0_0_0_z_yz_0_y, g_y_0_0_0_z_yz_0_z, g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_yz_0_x[i] = 2.0 * g_yz_yz_0_x[i] * a_exp;

        g_y_0_0_0_z_yz_0_y[i] = 2.0 * g_yz_yz_0_y[i] * a_exp;

        g_y_0_0_0_z_yz_0_z[i] = 2.0 * g_yz_yz_0_z[i] * a_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_y_0_0_0_z_zz_0_x, g_y_0_0_0_z_zz_0_y, g_y_0_0_0_z_zz_0_z, g_yz_zz_0_x, g_yz_zz_0_y, g_yz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_zz_0_x[i] = 2.0 * g_yz_zz_0_x[i] * a_exp;

        g_y_0_0_0_z_zz_0_y[i] = 2.0 * g_yz_zz_0_y[i] * a_exp;

        g_y_0_0_0_z_zz_0_z[i] = 2.0 * g_yz_zz_0_z[i] * a_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_xz_xx_0_x, g_xz_xx_0_y, g_xz_xx_0_z, g_z_0_0_0_x_xx_0_x, g_z_0_0_0_x_xx_0_y, g_z_0_0_0_x_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xx_0_x[i] = 2.0 * g_xz_xx_0_x[i] * a_exp;

        g_z_0_0_0_x_xx_0_y[i] = 2.0 * g_xz_xx_0_y[i] * a_exp;

        g_z_0_0_0_x_xx_0_z[i] = 2.0 * g_xz_xx_0_z[i] * a_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_xz_xy_0_x, g_xz_xy_0_y, g_xz_xy_0_z, g_z_0_0_0_x_xy_0_x, g_z_0_0_0_x_xy_0_y, g_z_0_0_0_x_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xy_0_x[i] = 2.0 * g_xz_xy_0_x[i] * a_exp;

        g_z_0_0_0_x_xy_0_y[i] = 2.0 * g_xz_xy_0_y[i] * a_exp;

        g_z_0_0_0_x_xy_0_z[i] = 2.0 * g_xz_xy_0_z[i] * a_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_xz_xz_0_x, g_xz_xz_0_y, g_xz_xz_0_z, g_z_0_0_0_x_xz_0_x, g_z_0_0_0_x_xz_0_y, g_z_0_0_0_x_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_xz_0_x[i] = 2.0 * g_xz_xz_0_x[i] * a_exp;

        g_z_0_0_0_x_xz_0_y[i] = 2.0 * g_xz_xz_0_y[i] * a_exp;

        g_z_0_0_0_x_xz_0_z[i] = 2.0 * g_xz_xz_0_z[i] * a_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_xz_yy_0_x, g_xz_yy_0_y, g_xz_yy_0_z, g_z_0_0_0_x_yy_0_x, g_z_0_0_0_x_yy_0_y, g_z_0_0_0_x_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yy_0_x[i] = 2.0 * g_xz_yy_0_x[i] * a_exp;

        g_z_0_0_0_x_yy_0_y[i] = 2.0 * g_xz_yy_0_y[i] * a_exp;

        g_z_0_0_0_x_yy_0_z[i] = 2.0 * g_xz_yy_0_z[i] * a_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_xz_yz_0_x, g_xz_yz_0_y, g_xz_yz_0_z, g_z_0_0_0_x_yz_0_x, g_z_0_0_0_x_yz_0_y, g_z_0_0_0_x_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_yz_0_x[i] = 2.0 * g_xz_yz_0_x[i] * a_exp;

        g_z_0_0_0_x_yz_0_y[i] = 2.0 * g_xz_yz_0_y[i] * a_exp;

        g_z_0_0_0_x_yz_0_z[i] = 2.0 * g_xz_yz_0_z[i] * a_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_xz_zz_0_x, g_xz_zz_0_y, g_xz_zz_0_z, g_z_0_0_0_x_zz_0_x, g_z_0_0_0_x_zz_0_y, g_z_0_0_0_x_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_zz_0_x[i] = 2.0 * g_xz_zz_0_x[i] * a_exp;

        g_z_0_0_0_x_zz_0_y[i] = 2.0 * g_xz_zz_0_y[i] * a_exp;

        g_z_0_0_0_x_zz_0_z[i] = 2.0 * g_xz_zz_0_z[i] * a_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_yz_xx_0_x, g_yz_xx_0_y, g_yz_xx_0_z, g_z_0_0_0_y_xx_0_x, g_z_0_0_0_y_xx_0_y, g_z_0_0_0_y_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xx_0_x[i] = 2.0 * g_yz_xx_0_x[i] * a_exp;

        g_z_0_0_0_y_xx_0_y[i] = 2.0 * g_yz_xx_0_y[i] * a_exp;

        g_z_0_0_0_y_xx_0_z[i] = 2.0 * g_yz_xx_0_z[i] * a_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_yz_xy_0_x, g_yz_xy_0_y, g_yz_xy_0_z, g_z_0_0_0_y_xy_0_x, g_z_0_0_0_y_xy_0_y, g_z_0_0_0_y_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xy_0_x[i] = 2.0 * g_yz_xy_0_x[i] * a_exp;

        g_z_0_0_0_y_xy_0_y[i] = 2.0 * g_yz_xy_0_y[i] * a_exp;

        g_z_0_0_0_y_xy_0_z[i] = 2.0 * g_yz_xy_0_z[i] * a_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_yz_xz_0_x, g_yz_xz_0_y, g_yz_xz_0_z, g_z_0_0_0_y_xz_0_x, g_z_0_0_0_y_xz_0_y, g_z_0_0_0_y_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_xz_0_x[i] = 2.0 * g_yz_xz_0_x[i] * a_exp;

        g_z_0_0_0_y_xz_0_y[i] = 2.0 * g_yz_xz_0_y[i] * a_exp;

        g_z_0_0_0_y_xz_0_z[i] = 2.0 * g_yz_xz_0_z[i] * a_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_yz_yy_0_x, g_yz_yy_0_y, g_yz_yy_0_z, g_z_0_0_0_y_yy_0_x, g_z_0_0_0_y_yy_0_y, g_z_0_0_0_y_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yy_0_x[i] = 2.0 * g_yz_yy_0_x[i] * a_exp;

        g_z_0_0_0_y_yy_0_y[i] = 2.0 * g_yz_yy_0_y[i] * a_exp;

        g_z_0_0_0_y_yy_0_z[i] = 2.0 * g_yz_yy_0_z[i] * a_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_yz_yz_0_x, g_yz_yz_0_y, g_yz_yz_0_z, g_z_0_0_0_y_yz_0_x, g_z_0_0_0_y_yz_0_y, g_z_0_0_0_y_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_yz_0_x[i] = 2.0 * g_yz_yz_0_x[i] * a_exp;

        g_z_0_0_0_y_yz_0_y[i] = 2.0 * g_yz_yz_0_y[i] * a_exp;

        g_z_0_0_0_y_yz_0_z[i] = 2.0 * g_yz_yz_0_z[i] * a_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_yz_zz_0_x, g_yz_zz_0_y, g_yz_zz_0_z, g_z_0_0_0_y_zz_0_x, g_z_0_0_0_y_zz_0_y, g_z_0_0_0_y_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_zz_0_x[i] = 2.0 * g_yz_zz_0_x[i] * a_exp;

        g_z_0_0_0_y_zz_0_y[i] = 2.0 * g_yz_zz_0_y[i] * a_exp;

        g_z_0_0_0_y_zz_0_z[i] = 2.0 * g_yz_zz_0_z[i] * a_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_0_xx_0_x, g_0_xx_0_y, g_0_xx_0_z, g_z_0_0_0_z_xx_0_x, g_z_0_0_0_z_xx_0_y, g_z_0_0_0_z_xx_0_z, g_zz_xx_0_x, g_zz_xx_0_y, g_zz_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xx_0_x[i] = -g_0_xx_0_x[i] + 2.0 * g_zz_xx_0_x[i] * a_exp;

        g_z_0_0_0_z_xx_0_y[i] = -g_0_xx_0_y[i] + 2.0 * g_zz_xx_0_y[i] * a_exp;

        g_z_0_0_0_z_xx_0_z[i] = -g_0_xx_0_z[i] + 2.0 * g_zz_xx_0_z[i] * a_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_0_xy_0_x, g_0_xy_0_y, g_0_xy_0_z, g_z_0_0_0_z_xy_0_x, g_z_0_0_0_z_xy_0_y, g_z_0_0_0_z_xy_0_z, g_zz_xy_0_x, g_zz_xy_0_y, g_zz_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xy_0_x[i] = -g_0_xy_0_x[i] + 2.0 * g_zz_xy_0_x[i] * a_exp;

        g_z_0_0_0_z_xy_0_y[i] = -g_0_xy_0_y[i] + 2.0 * g_zz_xy_0_y[i] * a_exp;

        g_z_0_0_0_z_xy_0_z[i] = -g_0_xy_0_z[i] + 2.0 * g_zz_xy_0_z[i] * a_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_0_xz_0_x, g_0_xz_0_y, g_0_xz_0_z, g_z_0_0_0_z_xz_0_x, g_z_0_0_0_z_xz_0_y, g_z_0_0_0_z_xz_0_z, g_zz_xz_0_x, g_zz_xz_0_y, g_zz_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_xz_0_x[i] = -g_0_xz_0_x[i] + 2.0 * g_zz_xz_0_x[i] * a_exp;

        g_z_0_0_0_z_xz_0_y[i] = -g_0_xz_0_y[i] + 2.0 * g_zz_xz_0_y[i] * a_exp;

        g_z_0_0_0_z_xz_0_z[i] = -g_0_xz_0_z[i] + 2.0 * g_zz_xz_0_z[i] * a_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_0_yy_0_x, g_0_yy_0_y, g_0_yy_0_z, g_z_0_0_0_z_yy_0_x, g_z_0_0_0_z_yy_0_y, g_z_0_0_0_z_yy_0_z, g_zz_yy_0_x, g_zz_yy_0_y, g_zz_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yy_0_x[i] = -g_0_yy_0_x[i] + 2.0 * g_zz_yy_0_x[i] * a_exp;

        g_z_0_0_0_z_yy_0_y[i] = -g_0_yy_0_y[i] + 2.0 * g_zz_yy_0_y[i] * a_exp;

        g_z_0_0_0_z_yy_0_z[i] = -g_0_yy_0_z[i] + 2.0 * g_zz_yy_0_z[i] * a_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_0_yz_0_x, g_0_yz_0_y, g_0_yz_0_z, g_z_0_0_0_z_yz_0_x, g_z_0_0_0_z_yz_0_y, g_z_0_0_0_z_yz_0_z, g_zz_yz_0_x, g_zz_yz_0_y, g_zz_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_yz_0_x[i] = -g_0_yz_0_x[i] + 2.0 * g_zz_yz_0_x[i] * a_exp;

        g_z_0_0_0_z_yz_0_y[i] = -g_0_yz_0_y[i] + 2.0 * g_zz_yz_0_y[i] * a_exp;

        g_z_0_0_0_z_yz_0_z[i] = -g_0_yz_0_z[i] + 2.0 * g_zz_yz_0_z[i] * a_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_0_zz_0_x, g_0_zz_0_y, g_0_zz_0_z, g_z_0_0_0_z_zz_0_x, g_z_0_0_0_z_zz_0_y, g_z_0_0_0_z_zz_0_z, g_zz_zz_0_x, g_zz_zz_0_y, g_zz_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_zz_0_x[i] = -g_0_zz_0_x[i] + 2.0 * g_zz_zz_0_x[i] * a_exp;

        g_z_0_0_0_z_zz_0_y[i] = -g_0_zz_0_y[i] + 2.0 * g_zz_zz_0_y[i] * a_exp;

        g_z_0_0_0_z_zz_0_z[i] = -g_0_zz_0_z[i] + 2.0 * g_zz_zz_0_z[i] * a_exp;
    }
}

} // t4c_geom namespace

