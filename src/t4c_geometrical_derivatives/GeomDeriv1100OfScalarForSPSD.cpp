#include "GeomDeriv1100OfScalarForSPSD.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_spsd_0(CSimdArray<double>& buffer_1100_spsd,
                     const CSimdArray<double>& buffer_pssd,
                     const CSimdArray<double>& buffer_pdsd,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_spsd.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pssd

    auto g_x_0_0_xx = buffer_pssd[0];

    auto g_x_0_0_xy = buffer_pssd[1];

    auto g_x_0_0_xz = buffer_pssd[2];

    auto g_x_0_0_yy = buffer_pssd[3];

    auto g_x_0_0_yz = buffer_pssd[4];

    auto g_x_0_0_zz = buffer_pssd[5];

    auto g_y_0_0_xx = buffer_pssd[6];

    auto g_y_0_0_xy = buffer_pssd[7];

    auto g_y_0_0_xz = buffer_pssd[8];

    auto g_y_0_0_yy = buffer_pssd[9];

    auto g_y_0_0_yz = buffer_pssd[10];

    auto g_y_0_0_zz = buffer_pssd[11];

    auto g_z_0_0_xx = buffer_pssd[12];

    auto g_z_0_0_xy = buffer_pssd[13];

    auto g_z_0_0_xz = buffer_pssd[14];

    auto g_z_0_0_yy = buffer_pssd[15];

    auto g_z_0_0_yz = buffer_pssd[16];

    auto g_z_0_0_zz = buffer_pssd[17];

    /// Set up components of auxilary buffer : buffer_pdsd

    auto g_x_xx_0_xx = buffer_pdsd[0];

    auto g_x_xx_0_xy = buffer_pdsd[1];

    auto g_x_xx_0_xz = buffer_pdsd[2];

    auto g_x_xx_0_yy = buffer_pdsd[3];

    auto g_x_xx_0_yz = buffer_pdsd[4];

    auto g_x_xx_0_zz = buffer_pdsd[5];

    auto g_x_xy_0_xx = buffer_pdsd[6];

    auto g_x_xy_0_xy = buffer_pdsd[7];

    auto g_x_xy_0_xz = buffer_pdsd[8];

    auto g_x_xy_0_yy = buffer_pdsd[9];

    auto g_x_xy_0_yz = buffer_pdsd[10];

    auto g_x_xy_0_zz = buffer_pdsd[11];

    auto g_x_xz_0_xx = buffer_pdsd[12];

    auto g_x_xz_0_xy = buffer_pdsd[13];

    auto g_x_xz_0_xz = buffer_pdsd[14];

    auto g_x_xz_0_yy = buffer_pdsd[15];

    auto g_x_xz_0_yz = buffer_pdsd[16];

    auto g_x_xz_0_zz = buffer_pdsd[17];

    auto g_x_yy_0_xx = buffer_pdsd[18];

    auto g_x_yy_0_xy = buffer_pdsd[19];

    auto g_x_yy_0_xz = buffer_pdsd[20];

    auto g_x_yy_0_yy = buffer_pdsd[21];

    auto g_x_yy_0_yz = buffer_pdsd[22];

    auto g_x_yy_0_zz = buffer_pdsd[23];

    auto g_x_yz_0_xx = buffer_pdsd[24];

    auto g_x_yz_0_xy = buffer_pdsd[25];

    auto g_x_yz_0_xz = buffer_pdsd[26];

    auto g_x_yz_0_yy = buffer_pdsd[27];

    auto g_x_yz_0_yz = buffer_pdsd[28];

    auto g_x_yz_0_zz = buffer_pdsd[29];

    auto g_x_zz_0_xx = buffer_pdsd[30];

    auto g_x_zz_0_xy = buffer_pdsd[31];

    auto g_x_zz_0_xz = buffer_pdsd[32];

    auto g_x_zz_0_yy = buffer_pdsd[33];

    auto g_x_zz_0_yz = buffer_pdsd[34];

    auto g_x_zz_0_zz = buffer_pdsd[35];

    auto g_y_xx_0_xx = buffer_pdsd[36];

    auto g_y_xx_0_xy = buffer_pdsd[37];

    auto g_y_xx_0_xz = buffer_pdsd[38];

    auto g_y_xx_0_yy = buffer_pdsd[39];

    auto g_y_xx_0_yz = buffer_pdsd[40];

    auto g_y_xx_0_zz = buffer_pdsd[41];

    auto g_y_xy_0_xx = buffer_pdsd[42];

    auto g_y_xy_0_xy = buffer_pdsd[43];

    auto g_y_xy_0_xz = buffer_pdsd[44];

    auto g_y_xy_0_yy = buffer_pdsd[45];

    auto g_y_xy_0_yz = buffer_pdsd[46];

    auto g_y_xy_0_zz = buffer_pdsd[47];

    auto g_y_xz_0_xx = buffer_pdsd[48];

    auto g_y_xz_0_xy = buffer_pdsd[49];

    auto g_y_xz_0_xz = buffer_pdsd[50];

    auto g_y_xz_0_yy = buffer_pdsd[51];

    auto g_y_xz_0_yz = buffer_pdsd[52];

    auto g_y_xz_0_zz = buffer_pdsd[53];

    auto g_y_yy_0_xx = buffer_pdsd[54];

    auto g_y_yy_0_xy = buffer_pdsd[55];

    auto g_y_yy_0_xz = buffer_pdsd[56];

    auto g_y_yy_0_yy = buffer_pdsd[57];

    auto g_y_yy_0_yz = buffer_pdsd[58];

    auto g_y_yy_0_zz = buffer_pdsd[59];

    auto g_y_yz_0_xx = buffer_pdsd[60];

    auto g_y_yz_0_xy = buffer_pdsd[61];

    auto g_y_yz_0_xz = buffer_pdsd[62];

    auto g_y_yz_0_yy = buffer_pdsd[63];

    auto g_y_yz_0_yz = buffer_pdsd[64];

    auto g_y_yz_0_zz = buffer_pdsd[65];

    auto g_y_zz_0_xx = buffer_pdsd[66];

    auto g_y_zz_0_xy = buffer_pdsd[67];

    auto g_y_zz_0_xz = buffer_pdsd[68];

    auto g_y_zz_0_yy = buffer_pdsd[69];

    auto g_y_zz_0_yz = buffer_pdsd[70];

    auto g_y_zz_0_zz = buffer_pdsd[71];

    auto g_z_xx_0_xx = buffer_pdsd[72];

    auto g_z_xx_0_xy = buffer_pdsd[73];

    auto g_z_xx_0_xz = buffer_pdsd[74];

    auto g_z_xx_0_yy = buffer_pdsd[75];

    auto g_z_xx_0_yz = buffer_pdsd[76];

    auto g_z_xx_0_zz = buffer_pdsd[77];

    auto g_z_xy_0_xx = buffer_pdsd[78];

    auto g_z_xy_0_xy = buffer_pdsd[79];

    auto g_z_xy_0_xz = buffer_pdsd[80];

    auto g_z_xy_0_yy = buffer_pdsd[81];

    auto g_z_xy_0_yz = buffer_pdsd[82];

    auto g_z_xy_0_zz = buffer_pdsd[83];

    auto g_z_xz_0_xx = buffer_pdsd[84];

    auto g_z_xz_0_xy = buffer_pdsd[85];

    auto g_z_xz_0_xz = buffer_pdsd[86];

    auto g_z_xz_0_yy = buffer_pdsd[87];

    auto g_z_xz_0_yz = buffer_pdsd[88];

    auto g_z_xz_0_zz = buffer_pdsd[89];

    auto g_z_yy_0_xx = buffer_pdsd[90];

    auto g_z_yy_0_xy = buffer_pdsd[91];

    auto g_z_yy_0_xz = buffer_pdsd[92];

    auto g_z_yy_0_yy = buffer_pdsd[93];

    auto g_z_yy_0_yz = buffer_pdsd[94];

    auto g_z_yy_0_zz = buffer_pdsd[95];

    auto g_z_yz_0_xx = buffer_pdsd[96];

    auto g_z_yz_0_xy = buffer_pdsd[97];

    auto g_z_yz_0_xz = buffer_pdsd[98];

    auto g_z_yz_0_yy = buffer_pdsd[99];

    auto g_z_yz_0_yz = buffer_pdsd[100];

    auto g_z_yz_0_zz = buffer_pdsd[101];

    auto g_z_zz_0_xx = buffer_pdsd[102];

    auto g_z_zz_0_xy = buffer_pdsd[103];

    auto g_z_zz_0_xz = buffer_pdsd[104];

    auto g_z_zz_0_yy = buffer_pdsd[105];

    auto g_z_zz_0_yz = buffer_pdsd[106];

    auto g_z_zz_0_zz = buffer_pdsd[107];

    /// Set up components of integrals buffer : buffer_1100_spsd

    auto g_x_x_0_0_0_x_0_xx = buffer_1100_spsd[0];

    auto g_x_x_0_0_0_x_0_xy = buffer_1100_spsd[1];

    auto g_x_x_0_0_0_x_0_xz = buffer_1100_spsd[2];

    auto g_x_x_0_0_0_x_0_yy = buffer_1100_spsd[3];

    auto g_x_x_0_0_0_x_0_yz = buffer_1100_spsd[4];

    auto g_x_x_0_0_0_x_0_zz = buffer_1100_spsd[5];

    auto g_x_x_0_0_0_y_0_xx = buffer_1100_spsd[6];

    auto g_x_x_0_0_0_y_0_xy = buffer_1100_spsd[7];

    auto g_x_x_0_0_0_y_0_xz = buffer_1100_spsd[8];

    auto g_x_x_0_0_0_y_0_yy = buffer_1100_spsd[9];

    auto g_x_x_0_0_0_y_0_yz = buffer_1100_spsd[10];

    auto g_x_x_0_0_0_y_0_zz = buffer_1100_spsd[11];

    auto g_x_x_0_0_0_z_0_xx = buffer_1100_spsd[12];

    auto g_x_x_0_0_0_z_0_xy = buffer_1100_spsd[13];

    auto g_x_x_0_0_0_z_0_xz = buffer_1100_spsd[14];

    auto g_x_x_0_0_0_z_0_yy = buffer_1100_spsd[15];

    auto g_x_x_0_0_0_z_0_yz = buffer_1100_spsd[16];

    auto g_x_x_0_0_0_z_0_zz = buffer_1100_spsd[17];

    auto g_x_y_0_0_0_x_0_xx = buffer_1100_spsd[18];

    auto g_x_y_0_0_0_x_0_xy = buffer_1100_spsd[19];

    auto g_x_y_0_0_0_x_0_xz = buffer_1100_spsd[20];

    auto g_x_y_0_0_0_x_0_yy = buffer_1100_spsd[21];

    auto g_x_y_0_0_0_x_0_yz = buffer_1100_spsd[22];

    auto g_x_y_0_0_0_x_0_zz = buffer_1100_spsd[23];

    auto g_x_y_0_0_0_y_0_xx = buffer_1100_spsd[24];

    auto g_x_y_0_0_0_y_0_xy = buffer_1100_spsd[25];

    auto g_x_y_0_0_0_y_0_xz = buffer_1100_spsd[26];

    auto g_x_y_0_0_0_y_0_yy = buffer_1100_spsd[27];

    auto g_x_y_0_0_0_y_0_yz = buffer_1100_spsd[28];

    auto g_x_y_0_0_0_y_0_zz = buffer_1100_spsd[29];

    auto g_x_y_0_0_0_z_0_xx = buffer_1100_spsd[30];

    auto g_x_y_0_0_0_z_0_xy = buffer_1100_spsd[31];

    auto g_x_y_0_0_0_z_0_xz = buffer_1100_spsd[32];

    auto g_x_y_0_0_0_z_0_yy = buffer_1100_spsd[33];

    auto g_x_y_0_0_0_z_0_yz = buffer_1100_spsd[34];

    auto g_x_y_0_0_0_z_0_zz = buffer_1100_spsd[35];

    auto g_x_z_0_0_0_x_0_xx = buffer_1100_spsd[36];

    auto g_x_z_0_0_0_x_0_xy = buffer_1100_spsd[37];

    auto g_x_z_0_0_0_x_0_xz = buffer_1100_spsd[38];

    auto g_x_z_0_0_0_x_0_yy = buffer_1100_spsd[39];

    auto g_x_z_0_0_0_x_0_yz = buffer_1100_spsd[40];

    auto g_x_z_0_0_0_x_0_zz = buffer_1100_spsd[41];

    auto g_x_z_0_0_0_y_0_xx = buffer_1100_spsd[42];

    auto g_x_z_0_0_0_y_0_xy = buffer_1100_spsd[43];

    auto g_x_z_0_0_0_y_0_xz = buffer_1100_spsd[44];

    auto g_x_z_0_0_0_y_0_yy = buffer_1100_spsd[45];

    auto g_x_z_0_0_0_y_0_yz = buffer_1100_spsd[46];

    auto g_x_z_0_0_0_y_0_zz = buffer_1100_spsd[47];

    auto g_x_z_0_0_0_z_0_xx = buffer_1100_spsd[48];

    auto g_x_z_0_0_0_z_0_xy = buffer_1100_spsd[49];

    auto g_x_z_0_0_0_z_0_xz = buffer_1100_spsd[50];

    auto g_x_z_0_0_0_z_0_yy = buffer_1100_spsd[51];

    auto g_x_z_0_0_0_z_0_yz = buffer_1100_spsd[52];

    auto g_x_z_0_0_0_z_0_zz = buffer_1100_spsd[53];

    auto g_y_x_0_0_0_x_0_xx = buffer_1100_spsd[54];

    auto g_y_x_0_0_0_x_0_xy = buffer_1100_spsd[55];

    auto g_y_x_0_0_0_x_0_xz = buffer_1100_spsd[56];

    auto g_y_x_0_0_0_x_0_yy = buffer_1100_spsd[57];

    auto g_y_x_0_0_0_x_0_yz = buffer_1100_spsd[58];

    auto g_y_x_0_0_0_x_0_zz = buffer_1100_spsd[59];

    auto g_y_x_0_0_0_y_0_xx = buffer_1100_spsd[60];

    auto g_y_x_0_0_0_y_0_xy = buffer_1100_spsd[61];

    auto g_y_x_0_0_0_y_0_xz = buffer_1100_spsd[62];

    auto g_y_x_0_0_0_y_0_yy = buffer_1100_spsd[63];

    auto g_y_x_0_0_0_y_0_yz = buffer_1100_spsd[64];

    auto g_y_x_0_0_0_y_0_zz = buffer_1100_spsd[65];

    auto g_y_x_0_0_0_z_0_xx = buffer_1100_spsd[66];

    auto g_y_x_0_0_0_z_0_xy = buffer_1100_spsd[67];

    auto g_y_x_0_0_0_z_0_xz = buffer_1100_spsd[68];

    auto g_y_x_0_0_0_z_0_yy = buffer_1100_spsd[69];

    auto g_y_x_0_0_0_z_0_yz = buffer_1100_spsd[70];

    auto g_y_x_0_0_0_z_0_zz = buffer_1100_spsd[71];

    auto g_y_y_0_0_0_x_0_xx = buffer_1100_spsd[72];

    auto g_y_y_0_0_0_x_0_xy = buffer_1100_spsd[73];

    auto g_y_y_0_0_0_x_0_xz = buffer_1100_spsd[74];

    auto g_y_y_0_0_0_x_0_yy = buffer_1100_spsd[75];

    auto g_y_y_0_0_0_x_0_yz = buffer_1100_spsd[76];

    auto g_y_y_0_0_0_x_0_zz = buffer_1100_spsd[77];

    auto g_y_y_0_0_0_y_0_xx = buffer_1100_spsd[78];

    auto g_y_y_0_0_0_y_0_xy = buffer_1100_spsd[79];

    auto g_y_y_0_0_0_y_0_xz = buffer_1100_spsd[80];

    auto g_y_y_0_0_0_y_0_yy = buffer_1100_spsd[81];

    auto g_y_y_0_0_0_y_0_yz = buffer_1100_spsd[82];

    auto g_y_y_0_0_0_y_0_zz = buffer_1100_spsd[83];

    auto g_y_y_0_0_0_z_0_xx = buffer_1100_spsd[84];

    auto g_y_y_0_0_0_z_0_xy = buffer_1100_spsd[85];

    auto g_y_y_0_0_0_z_0_xz = buffer_1100_spsd[86];

    auto g_y_y_0_0_0_z_0_yy = buffer_1100_spsd[87];

    auto g_y_y_0_0_0_z_0_yz = buffer_1100_spsd[88];

    auto g_y_y_0_0_0_z_0_zz = buffer_1100_spsd[89];

    auto g_y_z_0_0_0_x_0_xx = buffer_1100_spsd[90];

    auto g_y_z_0_0_0_x_0_xy = buffer_1100_spsd[91];

    auto g_y_z_0_0_0_x_0_xz = buffer_1100_spsd[92];

    auto g_y_z_0_0_0_x_0_yy = buffer_1100_spsd[93];

    auto g_y_z_0_0_0_x_0_yz = buffer_1100_spsd[94];

    auto g_y_z_0_0_0_x_0_zz = buffer_1100_spsd[95];

    auto g_y_z_0_0_0_y_0_xx = buffer_1100_spsd[96];

    auto g_y_z_0_0_0_y_0_xy = buffer_1100_spsd[97];

    auto g_y_z_0_0_0_y_0_xz = buffer_1100_spsd[98];

    auto g_y_z_0_0_0_y_0_yy = buffer_1100_spsd[99];

    auto g_y_z_0_0_0_y_0_yz = buffer_1100_spsd[100];

    auto g_y_z_0_0_0_y_0_zz = buffer_1100_spsd[101];

    auto g_y_z_0_0_0_z_0_xx = buffer_1100_spsd[102];

    auto g_y_z_0_0_0_z_0_xy = buffer_1100_spsd[103];

    auto g_y_z_0_0_0_z_0_xz = buffer_1100_spsd[104];

    auto g_y_z_0_0_0_z_0_yy = buffer_1100_spsd[105];

    auto g_y_z_0_0_0_z_0_yz = buffer_1100_spsd[106];

    auto g_y_z_0_0_0_z_0_zz = buffer_1100_spsd[107];

    auto g_z_x_0_0_0_x_0_xx = buffer_1100_spsd[108];

    auto g_z_x_0_0_0_x_0_xy = buffer_1100_spsd[109];

    auto g_z_x_0_0_0_x_0_xz = buffer_1100_spsd[110];

    auto g_z_x_0_0_0_x_0_yy = buffer_1100_spsd[111];

    auto g_z_x_0_0_0_x_0_yz = buffer_1100_spsd[112];

    auto g_z_x_0_0_0_x_0_zz = buffer_1100_spsd[113];

    auto g_z_x_0_0_0_y_0_xx = buffer_1100_spsd[114];

    auto g_z_x_0_0_0_y_0_xy = buffer_1100_spsd[115];

    auto g_z_x_0_0_0_y_0_xz = buffer_1100_spsd[116];

    auto g_z_x_0_0_0_y_0_yy = buffer_1100_spsd[117];

    auto g_z_x_0_0_0_y_0_yz = buffer_1100_spsd[118];

    auto g_z_x_0_0_0_y_0_zz = buffer_1100_spsd[119];

    auto g_z_x_0_0_0_z_0_xx = buffer_1100_spsd[120];

    auto g_z_x_0_0_0_z_0_xy = buffer_1100_spsd[121];

    auto g_z_x_0_0_0_z_0_xz = buffer_1100_spsd[122];

    auto g_z_x_0_0_0_z_0_yy = buffer_1100_spsd[123];

    auto g_z_x_0_0_0_z_0_yz = buffer_1100_spsd[124];

    auto g_z_x_0_0_0_z_0_zz = buffer_1100_spsd[125];

    auto g_z_y_0_0_0_x_0_xx = buffer_1100_spsd[126];

    auto g_z_y_0_0_0_x_0_xy = buffer_1100_spsd[127];

    auto g_z_y_0_0_0_x_0_xz = buffer_1100_spsd[128];

    auto g_z_y_0_0_0_x_0_yy = buffer_1100_spsd[129];

    auto g_z_y_0_0_0_x_0_yz = buffer_1100_spsd[130];

    auto g_z_y_0_0_0_x_0_zz = buffer_1100_spsd[131];

    auto g_z_y_0_0_0_y_0_xx = buffer_1100_spsd[132];

    auto g_z_y_0_0_0_y_0_xy = buffer_1100_spsd[133];

    auto g_z_y_0_0_0_y_0_xz = buffer_1100_spsd[134];

    auto g_z_y_0_0_0_y_0_yy = buffer_1100_spsd[135];

    auto g_z_y_0_0_0_y_0_yz = buffer_1100_spsd[136];

    auto g_z_y_0_0_0_y_0_zz = buffer_1100_spsd[137];

    auto g_z_y_0_0_0_z_0_xx = buffer_1100_spsd[138];

    auto g_z_y_0_0_0_z_0_xy = buffer_1100_spsd[139];

    auto g_z_y_0_0_0_z_0_xz = buffer_1100_spsd[140];

    auto g_z_y_0_0_0_z_0_yy = buffer_1100_spsd[141];

    auto g_z_y_0_0_0_z_0_yz = buffer_1100_spsd[142];

    auto g_z_y_0_0_0_z_0_zz = buffer_1100_spsd[143];

    auto g_z_z_0_0_0_x_0_xx = buffer_1100_spsd[144];

    auto g_z_z_0_0_0_x_0_xy = buffer_1100_spsd[145];

    auto g_z_z_0_0_0_x_0_xz = buffer_1100_spsd[146];

    auto g_z_z_0_0_0_x_0_yy = buffer_1100_spsd[147];

    auto g_z_z_0_0_0_x_0_yz = buffer_1100_spsd[148];

    auto g_z_z_0_0_0_x_0_zz = buffer_1100_spsd[149];

    auto g_z_z_0_0_0_y_0_xx = buffer_1100_spsd[150];

    auto g_z_z_0_0_0_y_0_xy = buffer_1100_spsd[151];

    auto g_z_z_0_0_0_y_0_xz = buffer_1100_spsd[152];

    auto g_z_z_0_0_0_y_0_yy = buffer_1100_spsd[153];

    auto g_z_z_0_0_0_y_0_yz = buffer_1100_spsd[154];

    auto g_z_z_0_0_0_y_0_zz = buffer_1100_spsd[155];

    auto g_z_z_0_0_0_z_0_xx = buffer_1100_spsd[156];

    auto g_z_z_0_0_0_z_0_xy = buffer_1100_spsd[157];

    auto g_z_z_0_0_0_z_0_xz = buffer_1100_spsd[158];

    auto g_z_z_0_0_0_z_0_yy = buffer_1100_spsd[159];

    auto g_z_z_0_0_0_z_0_yz = buffer_1100_spsd[160];

    auto g_z_z_0_0_0_z_0_zz = buffer_1100_spsd[161];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_0_xx, g_x_0_0_xy, g_x_0_0_xz, g_x_0_0_yy, g_x_0_0_yz, g_x_0_0_zz, g_x_x_0_0_0_x_0_xx, g_x_x_0_0_0_x_0_xy, g_x_x_0_0_0_x_0_xz, g_x_x_0_0_0_x_0_yy, g_x_x_0_0_0_x_0_yz, g_x_x_0_0_0_x_0_zz, g_x_xx_0_xx, g_x_xx_0_xy, g_x_xx_0_xz, g_x_xx_0_yy, g_x_xx_0_yz, g_x_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_x_0_xx[i] = -2.0 * g_x_0_0_xx[i] * a_exp + 4.0 * g_x_xx_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_0_xy[i] = -2.0 * g_x_0_0_xy[i] * a_exp + 4.0 * g_x_xx_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_0_xz[i] = -2.0 * g_x_0_0_xz[i] * a_exp + 4.0 * g_x_xx_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_0_yy[i] = -2.0 * g_x_0_0_yy[i] * a_exp + 4.0 * g_x_xx_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_0_yz[i] = -2.0 * g_x_0_0_yz[i] * a_exp + 4.0 * g_x_xx_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_0_zz[i] = -2.0 * g_x_0_0_zz[i] * a_exp + 4.0 * g_x_xx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_x_0_0_0_y_0_xx, g_x_x_0_0_0_y_0_xy, g_x_x_0_0_0_y_0_xz, g_x_x_0_0_0_y_0_yy, g_x_x_0_0_0_y_0_yz, g_x_x_0_0_0_y_0_zz, g_x_xy_0_xx, g_x_xy_0_xy, g_x_xy_0_xz, g_x_xy_0_yy, g_x_xy_0_yz, g_x_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_y_0_xx[i] = 4.0 * g_x_xy_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_0_xy[i] = 4.0 * g_x_xy_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_0_xz[i] = 4.0 * g_x_xy_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_0_yy[i] = 4.0 * g_x_xy_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_0_yz[i] = 4.0 * g_x_xy_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_0_zz[i] = 4.0 * g_x_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_x_0_0_0_z_0_xx, g_x_x_0_0_0_z_0_xy, g_x_x_0_0_0_z_0_xz, g_x_x_0_0_0_z_0_yy, g_x_x_0_0_0_z_0_yz, g_x_x_0_0_0_z_0_zz, g_x_xz_0_xx, g_x_xz_0_xy, g_x_xz_0_xz, g_x_xz_0_yy, g_x_xz_0_yz, g_x_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_z_0_xx[i] = 4.0 * g_x_xz_0_xx[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_0_xy[i] = 4.0 * g_x_xz_0_xy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_0_xz[i] = 4.0 * g_x_xz_0_xz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_0_yy[i] = 4.0 * g_x_xz_0_yy[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_0_yz[i] = 4.0 * g_x_xz_0_yz[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_0_zz[i] = 4.0 * g_x_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_xy_0_xx, g_x_xy_0_xy, g_x_xy_0_xz, g_x_xy_0_yy, g_x_xy_0_yz, g_x_xy_0_zz, g_x_y_0_0_0_x_0_xx, g_x_y_0_0_0_x_0_xy, g_x_y_0_0_0_x_0_xz, g_x_y_0_0_0_x_0_yy, g_x_y_0_0_0_x_0_yz, g_x_y_0_0_0_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_x_0_xx[i] = 4.0 * g_x_xy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_0_xy[i] = 4.0 * g_x_xy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_0_xz[i] = 4.0 * g_x_xy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_0_yy[i] = 4.0 * g_x_xy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_0_yz[i] = 4.0 * g_x_xy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_0_zz[i] = 4.0 * g_x_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_0_xx, g_x_0_0_xy, g_x_0_0_xz, g_x_0_0_yy, g_x_0_0_yz, g_x_0_0_zz, g_x_y_0_0_0_y_0_xx, g_x_y_0_0_0_y_0_xy, g_x_y_0_0_0_y_0_xz, g_x_y_0_0_0_y_0_yy, g_x_y_0_0_0_y_0_yz, g_x_y_0_0_0_y_0_zz, g_x_yy_0_xx, g_x_yy_0_xy, g_x_yy_0_xz, g_x_yy_0_yy, g_x_yy_0_yz, g_x_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_y_0_xx[i] = -2.0 * g_x_0_0_xx[i] * a_exp + 4.0 * g_x_yy_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_0_xy[i] = -2.0 * g_x_0_0_xy[i] * a_exp + 4.0 * g_x_yy_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_0_xz[i] = -2.0 * g_x_0_0_xz[i] * a_exp + 4.0 * g_x_yy_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_0_yy[i] = -2.0 * g_x_0_0_yy[i] * a_exp + 4.0 * g_x_yy_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_0_yz[i] = -2.0 * g_x_0_0_yz[i] * a_exp + 4.0 * g_x_yy_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_0_zz[i] = -2.0 * g_x_0_0_zz[i] * a_exp + 4.0 * g_x_yy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_y_0_0_0_z_0_xx, g_x_y_0_0_0_z_0_xy, g_x_y_0_0_0_z_0_xz, g_x_y_0_0_0_z_0_yy, g_x_y_0_0_0_z_0_yz, g_x_y_0_0_0_z_0_zz, g_x_yz_0_xx, g_x_yz_0_xy, g_x_yz_0_xz, g_x_yz_0_yy, g_x_yz_0_yz, g_x_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_z_0_xx[i] = 4.0 * g_x_yz_0_xx[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_0_xy[i] = 4.0 * g_x_yz_0_xy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_0_xz[i] = 4.0 * g_x_yz_0_xz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_0_yy[i] = 4.0 * g_x_yz_0_yy[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_0_yz[i] = 4.0 * g_x_yz_0_yz[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_0_zz[i] = 4.0 * g_x_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_x_xz_0_xx, g_x_xz_0_xy, g_x_xz_0_xz, g_x_xz_0_yy, g_x_xz_0_yz, g_x_xz_0_zz, g_x_z_0_0_0_x_0_xx, g_x_z_0_0_0_x_0_xy, g_x_z_0_0_0_x_0_xz, g_x_z_0_0_0_x_0_yy, g_x_z_0_0_0_x_0_yz, g_x_z_0_0_0_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_x_0_xx[i] = 4.0 * g_x_xz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_0_xy[i] = 4.0 * g_x_xz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_0_xz[i] = 4.0 * g_x_xz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_0_yy[i] = 4.0 * g_x_xz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_0_yz[i] = 4.0 * g_x_xz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_0_zz[i] = 4.0 * g_x_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_yz_0_xx, g_x_yz_0_xy, g_x_yz_0_xz, g_x_yz_0_yy, g_x_yz_0_yz, g_x_yz_0_zz, g_x_z_0_0_0_y_0_xx, g_x_z_0_0_0_y_0_xy, g_x_z_0_0_0_y_0_xz, g_x_z_0_0_0_y_0_yy, g_x_z_0_0_0_y_0_yz, g_x_z_0_0_0_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_y_0_xx[i] = 4.0 * g_x_yz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_0_xy[i] = 4.0 * g_x_yz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_0_xz[i] = 4.0 * g_x_yz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_0_yy[i] = 4.0 * g_x_yz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_0_yz[i] = 4.0 * g_x_yz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_0_zz[i] = 4.0 * g_x_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_0_0_xx, g_x_0_0_xy, g_x_0_0_xz, g_x_0_0_yy, g_x_0_0_yz, g_x_0_0_zz, g_x_z_0_0_0_z_0_xx, g_x_z_0_0_0_z_0_xy, g_x_z_0_0_0_z_0_xz, g_x_z_0_0_0_z_0_yy, g_x_z_0_0_0_z_0_yz, g_x_z_0_0_0_z_0_zz, g_x_zz_0_xx, g_x_zz_0_xy, g_x_zz_0_xz, g_x_zz_0_yy, g_x_zz_0_yz, g_x_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_z_0_xx[i] = -2.0 * g_x_0_0_xx[i] * a_exp + 4.0 * g_x_zz_0_xx[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_0_xy[i] = -2.0 * g_x_0_0_xy[i] * a_exp + 4.0 * g_x_zz_0_xy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_0_xz[i] = -2.0 * g_x_0_0_xz[i] * a_exp + 4.0 * g_x_zz_0_xz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_0_yy[i] = -2.0 * g_x_0_0_yy[i] * a_exp + 4.0 * g_x_zz_0_yy[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_0_yz[i] = -2.0 * g_x_0_0_yz[i] * a_exp + 4.0 * g_x_zz_0_yz[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_0_zz[i] = -2.0 * g_x_0_0_zz[i] * a_exp + 4.0 * g_x_zz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_y_0_0_xx, g_y_0_0_xy, g_y_0_0_xz, g_y_0_0_yy, g_y_0_0_yz, g_y_0_0_zz, g_y_x_0_0_0_x_0_xx, g_y_x_0_0_0_x_0_xy, g_y_x_0_0_0_x_0_xz, g_y_x_0_0_0_x_0_yy, g_y_x_0_0_0_x_0_yz, g_y_x_0_0_0_x_0_zz, g_y_xx_0_xx, g_y_xx_0_xy, g_y_xx_0_xz, g_y_xx_0_yy, g_y_xx_0_yz, g_y_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_x_0_xx[i] = -2.0 * g_y_0_0_xx[i] * a_exp + 4.0 * g_y_xx_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_0_xy[i] = -2.0 * g_y_0_0_xy[i] * a_exp + 4.0 * g_y_xx_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_0_xz[i] = -2.0 * g_y_0_0_xz[i] * a_exp + 4.0 * g_y_xx_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_0_yy[i] = -2.0 * g_y_0_0_yy[i] * a_exp + 4.0 * g_y_xx_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_0_yz[i] = -2.0 * g_y_0_0_yz[i] * a_exp + 4.0 * g_y_xx_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_0_zz[i] = -2.0 * g_y_0_0_zz[i] * a_exp + 4.0 * g_y_xx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_y_x_0_0_0_y_0_xx, g_y_x_0_0_0_y_0_xy, g_y_x_0_0_0_y_0_xz, g_y_x_0_0_0_y_0_yy, g_y_x_0_0_0_y_0_yz, g_y_x_0_0_0_y_0_zz, g_y_xy_0_xx, g_y_xy_0_xy, g_y_xy_0_xz, g_y_xy_0_yy, g_y_xy_0_yz, g_y_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_y_0_xx[i] = 4.0 * g_y_xy_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_0_xy[i] = 4.0 * g_y_xy_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_0_xz[i] = 4.0 * g_y_xy_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_0_yy[i] = 4.0 * g_y_xy_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_0_yz[i] = 4.0 * g_y_xy_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_0_zz[i] = 4.0 * g_y_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_y_x_0_0_0_z_0_xx, g_y_x_0_0_0_z_0_xy, g_y_x_0_0_0_z_0_xz, g_y_x_0_0_0_z_0_yy, g_y_x_0_0_0_z_0_yz, g_y_x_0_0_0_z_0_zz, g_y_xz_0_xx, g_y_xz_0_xy, g_y_xz_0_xz, g_y_xz_0_yy, g_y_xz_0_yz, g_y_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_z_0_xx[i] = 4.0 * g_y_xz_0_xx[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_0_xy[i] = 4.0 * g_y_xz_0_xy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_0_xz[i] = 4.0 * g_y_xz_0_xz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_0_yy[i] = 4.0 * g_y_xz_0_yy[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_0_yz[i] = 4.0 * g_y_xz_0_yz[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_0_zz[i] = 4.0 * g_y_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_y_xy_0_xx, g_y_xy_0_xy, g_y_xy_0_xz, g_y_xy_0_yy, g_y_xy_0_yz, g_y_xy_0_zz, g_y_y_0_0_0_x_0_xx, g_y_y_0_0_0_x_0_xy, g_y_y_0_0_0_x_0_xz, g_y_y_0_0_0_x_0_yy, g_y_y_0_0_0_x_0_yz, g_y_y_0_0_0_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_x_0_xx[i] = 4.0 * g_y_xy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_0_xy[i] = 4.0 * g_y_xy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_0_xz[i] = 4.0 * g_y_xy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_0_yy[i] = 4.0 * g_y_xy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_0_yz[i] = 4.0 * g_y_xy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_0_zz[i] = 4.0 * g_y_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_y_0_0_xx, g_y_0_0_xy, g_y_0_0_xz, g_y_0_0_yy, g_y_0_0_yz, g_y_0_0_zz, g_y_y_0_0_0_y_0_xx, g_y_y_0_0_0_y_0_xy, g_y_y_0_0_0_y_0_xz, g_y_y_0_0_0_y_0_yy, g_y_y_0_0_0_y_0_yz, g_y_y_0_0_0_y_0_zz, g_y_yy_0_xx, g_y_yy_0_xy, g_y_yy_0_xz, g_y_yy_0_yy, g_y_yy_0_yz, g_y_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_y_0_xx[i] = -2.0 * g_y_0_0_xx[i] * a_exp + 4.0 * g_y_yy_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_0_xy[i] = -2.0 * g_y_0_0_xy[i] * a_exp + 4.0 * g_y_yy_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_0_xz[i] = -2.0 * g_y_0_0_xz[i] * a_exp + 4.0 * g_y_yy_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_0_yy[i] = -2.0 * g_y_0_0_yy[i] * a_exp + 4.0 * g_y_yy_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_0_yz[i] = -2.0 * g_y_0_0_yz[i] * a_exp + 4.0 * g_y_yy_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_0_zz[i] = -2.0 * g_y_0_0_zz[i] * a_exp + 4.0 * g_y_yy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_y_y_0_0_0_z_0_xx, g_y_y_0_0_0_z_0_xy, g_y_y_0_0_0_z_0_xz, g_y_y_0_0_0_z_0_yy, g_y_y_0_0_0_z_0_yz, g_y_y_0_0_0_z_0_zz, g_y_yz_0_xx, g_y_yz_0_xy, g_y_yz_0_xz, g_y_yz_0_yy, g_y_yz_0_yz, g_y_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_z_0_xx[i] = 4.0 * g_y_yz_0_xx[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_0_xy[i] = 4.0 * g_y_yz_0_xy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_0_xz[i] = 4.0 * g_y_yz_0_xz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_0_yy[i] = 4.0 * g_y_yz_0_yy[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_0_yz[i] = 4.0 * g_y_yz_0_yz[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_0_zz[i] = 4.0 * g_y_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_y_xz_0_xx, g_y_xz_0_xy, g_y_xz_0_xz, g_y_xz_0_yy, g_y_xz_0_yz, g_y_xz_0_zz, g_y_z_0_0_0_x_0_xx, g_y_z_0_0_0_x_0_xy, g_y_z_0_0_0_x_0_xz, g_y_z_0_0_0_x_0_yy, g_y_z_0_0_0_x_0_yz, g_y_z_0_0_0_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_x_0_xx[i] = 4.0 * g_y_xz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_0_xy[i] = 4.0 * g_y_xz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_0_xz[i] = 4.0 * g_y_xz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_0_yy[i] = 4.0 * g_y_xz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_0_yz[i] = 4.0 * g_y_xz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_0_zz[i] = 4.0 * g_y_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_y_yz_0_xx, g_y_yz_0_xy, g_y_yz_0_xz, g_y_yz_0_yy, g_y_yz_0_yz, g_y_yz_0_zz, g_y_z_0_0_0_y_0_xx, g_y_z_0_0_0_y_0_xy, g_y_z_0_0_0_y_0_xz, g_y_z_0_0_0_y_0_yy, g_y_z_0_0_0_y_0_yz, g_y_z_0_0_0_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_y_0_xx[i] = 4.0 * g_y_yz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_0_xy[i] = 4.0 * g_y_yz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_0_xz[i] = 4.0 * g_y_yz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_0_yy[i] = 4.0 * g_y_yz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_0_yz[i] = 4.0 * g_y_yz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_0_zz[i] = 4.0 * g_y_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_y_0_0_xx, g_y_0_0_xy, g_y_0_0_xz, g_y_0_0_yy, g_y_0_0_yz, g_y_0_0_zz, g_y_z_0_0_0_z_0_xx, g_y_z_0_0_0_z_0_xy, g_y_z_0_0_0_z_0_xz, g_y_z_0_0_0_z_0_yy, g_y_z_0_0_0_z_0_yz, g_y_z_0_0_0_z_0_zz, g_y_zz_0_xx, g_y_zz_0_xy, g_y_zz_0_xz, g_y_zz_0_yy, g_y_zz_0_yz, g_y_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_z_0_xx[i] = -2.0 * g_y_0_0_xx[i] * a_exp + 4.0 * g_y_zz_0_xx[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_0_xy[i] = -2.0 * g_y_0_0_xy[i] * a_exp + 4.0 * g_y_zz_0_xy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_0_xz[i] = -2.0 * g_y_0_0_xz[i] * a_exp + 4.0 * g_y_zz_0_xz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_0_yy[i] = -2.0 * g_y_0_0_yy[i] * a_exp + 4.0 * g_y_zz_0_yy[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_0_yz[i] = -2.0 * g_y_0_0_yz[i] * a_exp + 4.0 * g_y_zz_0_yz[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_0_zz[i] = -2.0 * g_y_0_0_zz[i] * a_exp + 4.0 * g_y_zz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (108-114)

    #pragma omp simd aligned(g_z_0_0_xx, g_z_0_0_xy, g_z_0_0_xz, g_z_0_0_yy, g_z_0_0_yz, g_z_0_0_zz, g_z_x_0_0_0_x_0_xx, g_z_x_0_0_0_x_0_xy, g_z_x_0_0_0_x_0_xz, g_z_x_0_0_0_x_0_yy, g_z_x_0_0_0_x_0_yz, g_z_x_0_0_0_x_0_zz, g_z_xx_0_xx, g_z_xx_0_xy, g_z_xx_0_xz, g_z_xx_0_yy, g_z_xx_0_yz, g_z_xx_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_x_0_xx[i] = -2.0 * g_z_0_0_xx[i] * a_exp + 4.0 * g_z_xx_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_0_xy[i] = -2.0 * g_z_0_0_xy[i] * a_exp + 4.0 * g_z_xx_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_0_xz[i] = -2.0 * g_z_0_0_xz[i] * a_exp + 4.0 * g_z_xx_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_0_yy[i] = -2.0 * g_z_0_0_yy[i] * a_exp + 4.0 * g_z_xx_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_0_yz[i] = -2.0 * g_z_0_0_yz[i] * a_exp + 4.0 * g_z_xx_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_0_zz[i] = -2.0 * g_z_0_0_zz[i] * a_exp + 4.0 * g_z_xx_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (114-120)

    #pragma omp simd aligned(g_z_x_0_0_0_y_0_xx, g_z_x_0_0_0_y_0_xy, g_z_x_0_0_0_y_0_xz, g_z_x_0_0_0_y_0_yy, g_z_x_0_0_0_y_0_yz, g_z_x_0_0_0_y_0_zz, g_z_xy_0_xx, g_z_xy_0_xy, g_z_xy_0_xz, g_z_xy_0_yy, g_z_xy_0_yz, g_z_xy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_y_0_xx[i] = 4.0 * g_z_xy_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_0_xy[i] = 4.0 * g_z_xy_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_0_xz[i] = 4.0 * g_z_xy_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_0_yy[i] = 4.0 * g_z_xy_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_0_yz[i] = 4.0 * g_z_xy_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_0_zz[i] = 4.0 * g_z_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (120-126)

    #pragma omp simd aligned(g_z_x_0_0_0_z_0_xx, g_z_x_0_0_0_z_0_xy, g_z_x_0_0_0_z_0_xz, g_z_x_0_0_0_z_0_yy, g_z_x_0_0_0_z_0_yz, g_z_x_0_0_0_z_0_zz, g_z_xz_0_xx, g_z_xz_0_xy, g_z_xz_0_xz, g_z_xz_0_yy, g_z_xz_0_yz, g_z_xz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_z_0_xx[i] = 4.0 * g_z_xz_0_xx[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_0_xy[i] = 4.0 * g_z_xz_0_xy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_0_xz[i] = 4.0 * g_z_xz_0_xz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_0_yy[i] = 4.0 * g_z_xz_0_yy[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_0_yz[i] = 4.0 * g_z_xz_0_yz[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_0_zz[i] = 4.0 * g_z_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (126-132)

    #pragma omp simd aligned(g_z_xy_0_xx, g_z_xy_0_xy, g_z_xy_0_xz, g_z_xy_0_yy, g_z_xy_0_yz, g_z_xy_0_zz, g_z_y_0_0_0_x_0_xx, g_z_y_0_0_0_x_0_xy, g_z_y_0_0_0_x_0_xz, g_z_y_0_0_0_x_0_yy, g_z_y_0_0_0_x_0_yz, g_z_y_0_0_0_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_x_0_xx[i] = 4.0 * g_z_xy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_0_xy[i] = 4.0 * g_z_xy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_0_xz[i] = 4.0 * g_z_xy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_0_yy[i] = 4.0 * g_z_xy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_0_yz[i] = 4.0 * g_z_xy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_0_zz[i] = 4.0 * g_z_xy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (132-138)

    #pragma omp simd aligned(g_z_0_0_xx, g_z_0_0_xy, g_z_0_0_xz, g_z_0_0_yy, g_z_0_0_yz, g_z_0_0_zz, g_z_y_0_0_0_y_0_xx, g_z_y_0_0_0_y_0_xy, g_z_y_0_0_0_y_0_xz, g_z_y_0_0_0_y_0_yy, g_z_y_0_0_0_y_0_yz, g_z_y_0_0_0_y_0_zz, g_z_yy_0_xx, g_z_yy_0_xy, g_z_yy_0_xz, g_z_yy_0_yy, g_z_yy_0_yz, g_z_yy_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_y_0_xx[i] = -2.0 * g_z_0_0_xx[i] * a_exp + 4.0 * g_z_yy_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_0_xy[i] = -2.0 * g_z_0_0_xy[i] * a_exp + 4.0 * g_z_yy_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_0_xz[i] = -2.0 * g_z_0_0_xz[i] * a_exp + 4.0 * g_z_yy_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_0_yy[i] = -2.0 * g_z_0_0_yy[i] * a_exp + 4.0 * g_z_yy_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_0_yz[i] = -2.0 * g_z_0_0_yz[i] * a_exp + 4.0 * g_z_yy_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_0_zz[i] = -2.0 * g_z_0_0_zz[i] * a_exp + 4.0 * g_z_yy_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (138-144)

    #pragma omp simd aligned(g_z_y_0_0_0_z_0_xx, g_z_y_0_0_0_z_0_xy, g_z_y_0_0_0_z_0_xz, g_z_y_0_0_0_z_0_yy, g_z_y_0_0_0_z_0_yz, g_z_y_0_0_0_z_0_zz, g_z_yz_0_xx, g_z_yz_0_xy, g_z_yz_0_xz, g_z_yz_0_yy, g_z_yz_0_yz, g_z_yz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_z_0_xx[i] = 4.0 * g_z_yz_0_xx[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_0_xy[i] = 4.0 * g_z_yz_0_xy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_0_xz[i] = 4.0 * g_z_yz_0_xz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_0_yy[i] = 4.0 * g_z_yz_0_yy[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_0_yz[i] = 4.0 * g_z_yz_0_yz[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_0_zz[i] = 4.0 * g_z_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (144-150)

    #pragma omp simd aligned(g_z_xz_0_xx, g_z_xz_0_xy, g_z_xz_0_xz, g_z_xz_0_yy, g_z_xz_0_yz, g_z_xz_0_zz, g_z_z_0_0_0_x_0_xx, g_z_z_0_0_0_x_0_xy, g_z_z_0_0_0_x_0_xz, g_z_z_0_0_0_x_0_yy, g_z_z_0_0_0_x_0_yz, g_z_z_0_0_0_x_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_x_0_xx[i] = 4.0 * g_z_xz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_0_xy[i] = 4.0 * g_z_xz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_0_xz[i] = 4.0 * g_z_xz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_0_yy[i] = 4.0 * g_z_xz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_0_yz[i] = 4.0 * g_z_xz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_0_zz[i] = 4.0 * g_z_xz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (150-156)

    #pragma omp simd aligned(g_z_yz_0_xx, g_z_yz_0_xy, g_z_yz_0_xz, g_z_yz_0_yy, g_z_yz_0_yz, g_z_yz_0_zz, g_z_z_0_0_0_y_0_xx, g_z_z_0_0_0_y_0_xy, g_z_z_0_0_0_y_0_xz, g_z_z_0_0_0_y_0_yy, g_z_z_0_0_0_y_0_yz, g_z_z_0_0_0_y_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_y_0_xx[i] = 4.0 * g_z_yz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_0_xy[i] = 4.0 * g_z_yz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_0_xz[i] = 4.0 * g_z_yz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_0_yy[i] = 4.0 * g_z_yz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_0_yz[i] = 4.0 * g_z_yz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_0_zz[i] = 4.0 * g_z_yz_0_zz[i] * a_exp * b_exp;
    }
    // integrals block (156-162)

    #pragma omp simd aligned(g_z_0_0_xx, g_z_0_0_xy, g_z_0_0_xz, g_z_0_0_yy, g_z_0_0_yz, g_z_0_0_zz, g_z_z_0_0_0_z_0_xx, g_z_z_0_0_0_z_0_xy, g_z_z_0_0_0_z_0_xz, g_z_z_0_0_0_z_0_yy, g_z_z_0_0_0_z_0_yz, g_z_z_0_0_0_z_0_zz, g_z_zz_0_xx, g_z_zz_0_xy, g_z_zz_0_xz, g_z_zz_0_yy, g_z_zz_0_yz, g_z_zz_0_zz  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_z_0_xx[i] = -2.0 * g_z_0_0_xx[i] * a_exp + 4.0 * g_z_zz_0_xx[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_0_xy[i] = -2.0 * g_z_0_0_xy[i] * a_exp + 4.0 * g_z_zz_0_xy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_0_xz[i] = -2.0 * g_z_0_0_xz[i] * a_exp + 4.0 * g_z_zz_0_xz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_0_yy[i] = -2.0 * g_z_0_0_yy[i] * a_exp + 4.0 * g_z_zz_0_yy[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_0_yz[i] = -2.0 * g_z_0_0_yz[i] * a_exp + 4.0 * g_z_zz_0_yz[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_0_zz[i] = -2.0 * g_z_0_0_zz[i] * a_exp + 4.0 * g_z_zz_0_zz[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

