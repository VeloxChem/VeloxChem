#include "GeomDeriv1100OfScalarForSDSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_sdsp_0(CSimdArray<double>& buffer_1100_sdsp,
                     const CSimdArray<double>& buffer_ppsp,
                     const CSimdArray<double>& buffer_pfsp,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_sdsp.number_of_columns();

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

    /// Set up components of auxilary buffer : buffer_pfsp

    auto g_x_xxx_0_x = buffer_pfsp[0];

    auto g_x_xxx_0_y = buffer_pfsp[1];

    auto g_x_xxx_0_z = buffer_pfsp[2];

    auto g_x_xxy_0_x = buffer_pfsp[3];

    auto g_x_xxy_0_y = buffer_pfsp[4];

    auto g_x_xxy_0_z = buffer_pfsp[5];

    auto g_x_xxz_0_x = buffer_pfsp[6];

    auto g_x_xxz_0_y = buffer_pfsp[7];

    auto g_x_xxz_0_z = buffer_pfsp[8];

    auto g_x_xyy_0_x = buffer_pfsp[9];

    auto g_x_xyy_0_y = buffer_pfsp[10];

    auto g_x_xyy_0_z = buffer_pfsp[11];

    auto g_x_xyz_0_x = buffer_pfsp[12];

    auto g_x_xyz_0_y = buffer_pfsp[13];

    auto g_x_xyz_0_z = buffer_pfsp[14];

    auto g_x_xzz_0_x = buffer_pfsp[15];

    auto g_x_xzz_0_y = buffer_pfsp[16];

    auto g_x_xzz_0_z = buffer_pfsp[17];

    auto g_x_yyy_0_x = buffer_pfsp[18];

    auto g_x_yyy_0_y = buffer_pfsp[19];

    auto g_x_yyy_0_z = buffer_pfsp[20];

    auto g_x_yyz_0_x = buffer_pfsp[21];

    auto g_x_yyz_0_y = buffer_pfsp[22];

    auto g_x_yyz_0_z = buffer_pfsp[23];

    auto g_x_yzz_0_x = buffer_pfsp[24];

    auto g_x_yzz_0_y = buffer_pfsp[25];

    auto g_x_yzz_0_z = buffer_pfsp[26];

    auto g_x_zzz_0_x = buffer_pfsp[27];

    auto g_x_zzz_0_y = buffer_pfsp[28];

    auto g_x_zzz_0_z = buffer_pfsp[29];

    auto g_y_xxx_0_x = buffer_pfsp[30];

    auto g_y_xxx_0_y = buffer_pfsp[31];

    auto g_y_xxx_0_z = buffer_pfsp[32];

    auto g_y_xxy_0_x = buffer_pfsp[33];

    auto g_y_xxy_0_y = buffer_pfsp[34];

    auto g_y_xxy_0_z = buffer_pfsp[35];

    auto g_y_xxz_0_x = buffer_pfsp[36];

    auto g_y_xxz_0_y = buffer_pfsp[37];

    auto g_y_xxz_0_z = buffer_pfsp[38];

    auto g_y_xyy_0_x = buffer_pfsp[39];

    auto g_y_xyy_0_y = buffer_pfsp[40];

    auto g_y_xyy_0_z = buffer_pfsp[41];

    auto g_y_xyz_0_x = buffer_pfsp[42];

    auto g_y_xyz_0_y = buffer_pfsp[43];

    auto g_y_xyz_0_z = buffer_pfsp[44];

    auto g_y_xzz_0_x = buffer_pfsp[45];

    auto g_y_xzz_0_y = buffer_pfsp[46];

    auto g_y_xzz_0_z = buffer_pfsp[47];

    auto g_y_yyy_0_x = buffer_pfsp[48];

    auto g_y_yyy_0_y = buffer_pfsp[49];

    auto g_y_yyy_0_z = buffer_pfsp[50];

    auto g_y_yyz_0_x = buffer_pfsp[51];

    auto g_y_yyz_0_y = buffer_pfsp[52];

    auto g_y_yyz_0_z = buffer_pfsp[53];

    auto g_y_yzz_0_x = buffer_pfsp[54];

    auto g_y_yzz_0_y = buffer_pfsp[55];

    auto g_y_yzz_0_z = buffer_pfsp[56];

    auto g_y_zzz_0_x = buffer_pfsp[57];

    auto g_y_zzz_0_y = buffer_pfsp[58];

    auto g_y_zzz_0_z = buffer_pfsp[59];

    auto g_z_xxx_0_x = buffer_pfsp[60];

    auto g_z_xxx_0_y = buffer_pfsp[61];

    auto g_z_xxx_0_z = buffer_pfsp[62];

    auto g_z_xxy_0_x = buffer_pfsp[63];

    auto g_z_xxy_0_y = buffer_pfsp[64];

    auto g_z_xxy_0_z = buffer_pfsp[65];

    auto g_z_xxz_0_x = buffer_pfsp[66];

    auto g_z_xxz_0_y = buffer_pfsp[67];

    auto g_z_xxz_0_z = buffer_pfsp[68];

    auto g_z_xyy_0_x = buffer_pfsp[69];

    auto g_z_xyy_0_y = buffer_pfsp[70];

    auto g_z_xyy_0_z = buffer_pfsp[71];

    auto g_z_xyz_0_x = buffer_pfsp[72];

    auto g_z_xyz_0_y = buffer_pfsp[73];

    auto g_z_xyz_0_z = buffer_pfsp[74];

    auto g_z_xzz_0_x = buffer_pfsp[75];

    auto g_z_xzz_0_y = buffer_pfsp[76];

    auto g_z_xzz_0_z = buffer_pfsp[77];

    auto g_z_yyy_0_x = buffer_pfsp[78];

    auto g_z_yyy_0_y = buffer_pfsp[79];

    auto g_z_yyy_0_z = buffer_pfsp[80];

    auto g_z_yyz_0_x = buffer_pfsp[81];

    auto g_z_yyz_0_y = buffer_pfsp[82];

    auto g_z_yyz_0_z = buffer_pfsp[83];

    auto g_z_yzz_0_x = buffer_pfsp[84];

    auto g_z_yzz_0_y = buffer_pfsp[85];

    auto g_z_yzz_0_z = buffer_pfsp[86];

    auto g_z_zzz_0_x = buffer_pfsp[87];

    auto g_z_zzz_0_y = buffer_pfsp[88];

    auto g_z_zzz_0_z = buffer_pfsp[89];

    /// Set up components of integrals buffer : buffer_1100_sdsp

    auto g_x_x_0_0_0_xx_0_x = buffer_1100_sdsp[0];

    auto g_x_x_0_0_0_xx_0_y = buffer_1100_sdsp[1];

    auto g_x_x_0_0_0_xx_0_z = buffer_1100_sdsp[2];

    auto g_x_x_0_0_0_xy_0_x = buffer_1100_sdsp[3];

    auto g_x_x_0_0_0_xy_0_y = buffer_1100_sdsp[4];

    auto g_x_x_0_0_0_xy_0_z = buffer_1100_sdsp[5];

    auto g_x_x_0_0_0_xz_0_x = buffer_1100_sdsp[6];

    auto g_x_x_0_0_0_xz_0_y = buffer_1100_sdsp[7];

    auto g_x_x_0_0_0_xz_0_z = buffer_1100_sdsp[8];

    auto g_x_x_0_0_0_yy_0_x = buffer_1100_sdsp[9];

    auto g_x_x_0_0_0_yy_0_y = buffer_1100_sdsp[10];

    auto g_x_x_0_0_0_yy_0_z = buffer_1100_sdsp[11];

    auto g_x_x_0_0_0_yz_0_x = buffer_1100_sdsp[12];

    auto g_x_x_0_0_0_yz_0_y = buffer_1100_sdsp[13];

    auto g_x_x_0_0_0_yz_0_z = buffer_1100_sdsp[14];

    auto g_x_x_0_0_0_zz_0_x = buffer_1100_sdsp[15];

    auto g_x_x_0_0_0_zz_0_y = buffer_1100_sdsp[16];

    auto g_x_x_0_0_0_zz_0_z = buffer_1100_sdsp[17];

    auto g_x_y_0_0_0_xx_0_x = buffer_1100_sdsp[18];

    auto g_x_y_0_0_0_xx_0_y = buffer_1100_sdsp[19];

    auto g_x_y_0_0_0_xx_0_z = buffer_1100_sdsp[20];

    auto g_x_y_0_0_0_xy_0_x = buffer_1100_sdsp[21];

    auto g_x_y_0_0_0_xy_0_y = buffer_1100_sdsp[22];

    auto g_x_y_0_0_0_xy_0_z = buffer_1100_sdsp[23];

    auto g_x_y_0_0_0_xz_0_x = buffer_1100_sdsp[24];

    auto g_x_y_0_0_0_xz_0_y = buffer_1100_sdsp[25];

    auto g_x_y_0_0_0_xz_0_z = buffer_1100_sdsp[26];

    auto g_x_y_0_0_0_yy_0_x = buffer_1100_sdsp[27];

    auto g_x_y_0_0_0_yy_0_y = buffer_1100_sdsp[28];

    auto g_x_y_0_0_0_yy_0_z = buffer_1100_sdsp[29];

    auto g_x_y_0_0_0_yz_0_x = buffer_1100_sdsp[30];

    auto g_x_y_0_0_0_yz_0_y = buffer_1100_sdsp[31];

    auto g_x_y_0_0_0_yz_0_z = buffer_1100_sdsp[32];

    auto g_x_y_0_0_0_zz_0_x = buffer_1100_sdsp[33];

    auto g_x_y_0_0_0_zz_0_y = buffer_1100_sdsp[34];

    auto g_x_y_0_0_0_zz_0_z = buffer_1100_sdsp[35];

    auto g_x_z_0_0_0_xx_0_x = buffer_1100_sdsp[36];

    auto g_x_z_0_0_0_xx_0_y = buffer_1100_sdsp[37];

    auto g_x_z_0_0_0_xx_0_z = buffer_1100_sdsp[38];

    auto g_x_z_0_0_0_xy_0_x = buffer_1100_sdsp[39];

    auto g_x_z_0_0_0_xy_0_y = buffer_1100_sdsp[40];

    auto g_x_z_0_0_0_xy_0_z = buffer_1100_sdsp[41];

    auto g_x_z_0_0_0_xz_0_x = buffer_1100_sdsp[42];

    auto g_x_z_0_0_0_xz_0_y = buffer_1100_sdsp[43];

    auto g_x_z_0_0_0_xz_0_z = buffer_1100_sdsp[44];

    auto g_x_z_0_0_0_yy_0_x = buffer_1100_sdsp[45];

    auto g_x_z_0_0_0_yy_0_y = buffer_1100_sdsp[46];

    auto g_x_z_0_0_0_yy_0_z = buffer_1100_sdsp[47];

    auto g_x_z_0_0_0_yz_0_x = buffer_1100_sdsp[48];

    auto g_x_z_0_0_0_yz_0_y = buffer_1100_sdsp[49];

    auto g_x_z_0_0_0_yz_0_z = buffer_1100_sdsp[50];

    auto g_x_z_0_0_0_zz_0_x = buffer_1100_sdsp[51];

    auto g_x_z_0_0_0_zz_0_y = buffer_1100_sdsp[52];

    auto g_x_z_0_0_0_zz_0_z = buffer_1100_sdsp[53];

    auto g_y_x_0_0_0_xx_0_x = buffer_1100_sdsp[54];

    auto g_y_x_0_0_0_xx_0_y = buffer_1100_sdsp[55];

    auto g_y_x_0_0_0_xx_0_z = buffer_1100_sdsp[56];

    auto g_y_x_0_0_0_xy_0_x = buffer_1100_sdsp[57];

    auto g_y_x_0_0_0_xy_0_y = buffer_1100_sdsp[58];

    auto g_y_x_0_0_0_xy_0_z = buffer_1100_sdsp[59];

    auto g_y_x_0_0_0_xz_0_x = buffer_1100_sdsp[60];

    auto g_y_x_0_0_0_xz_0_y = buffer_1100_sdsp[61];

    auto g_y_x_0_0_0_xz_0_z = buffer_1100_sdsp[62];

    auto g_y_x_0_0_0_yy_0_x = buffer_1100_sdsp[63];

    auto g_y_x_0_0_0_yy_0_y = buffer_1100_sdsp[64];

    auto g_y_x_0_0_0_yy_0_z = buffer_1100_sdsp[65];

    auto g_y_x_0_0_0_yz_0_x = buffer_1100_sdsp[66];

    auto g_y_x_0_0_0_yz_0_y = buffer_1100_sdsp[67];

    auto g_y_x_0_0_0_yz_0_z = buffer_1100_sdsp[68];

    auto g_y_x_0_0_0_zz_0_x = buffer_1100_sdsp[69];

    auto g_y_x_0_0_0_zz_0_y = buffer_1100_sdsp[70];

    auto g_y_x_0_0_0_zz_0_z = buffer_1100_sdsp[71];

    auto g_y_y_0_0_0_xx_0_x = buffer_1100_sdsp[72];

    auto g_y_y_0_0_0_xx_0_y = buffer_1100_sdsp[73];

    auto g_y_y_0_0_0_xx_0_z = buffer_1100_sdsp[74];

    auto g_y_y_0_0_0_xy_0_x = buffer_1100_sdsp[75];

    auto g_y_y_0_0_0_xy_0_y = buffer_1100_sdsp[76];

    auto g_y_y_0_0_0_xy_0_z = buffer_1100_sdsp[77];

    auto g_y_y_0_0_0_xz_0_x = buffer_1100_sdsp[78];

    auto g_y_y_0_0_0_xz_0_y = buffer_1100_sdsp[79];

    auto g_y_y_0_0_0_xz_0_z = buffer_1100_sdsp[80];

    auto g_y_y_0_0_0_yy_0_x = buffer_1100_sdsp[81];

    auto g_y_y_0_0_0_yy_0_y = buffer_1100_sdsp[82];

    auto g_y_y_0_0_0_yy_0_z = buffer_1100_sdsp[83];

    auto g_y_y_0_0_0_yz_0_x = buffer_1100_sdsp[84];

    auto g_y_y_0_0_0_yz_0_y = buffer_1100_sdsp[85];

    auto g_y_y_0_0_0_yz_0_z = buffer_1100_sdsp[86];

    auto g_y_y_0_0_0_zz_0_x = buffer_1100_sdsp[87];

    auto g_y_y_0_0_0_zz_0_y = buffer_1100_sdsp[88];

    auto g_y_y_0_0_0_zz_0_z = buffer_1100_sdsp[89];

    auto g_y_z_0_0_0_xx_0_x = buffer_1100_sdsp[90];

    auto g_y_z_0_0_0_xx_0_y = buffer_1100_sdsp[91];

    auto g_y_z_0_0_0_xx_0_z = buffer_1100_sdsp[92];

    auto g_y_z_0_0_0_xy_0_x = buffer_1100_sdsp[93];

    auto g_y_z_0_0_0_xy_0_y = buffer_1100_sdsp[94];

    auto g_y_z_0_0_0_xy_0_z = buffer_1100_sdsp[95];

    auto g_y_z_0_0_0_xz_0_x = buffer_1100_sdsp[96];

    auto g_y_z_0_0_0_xz_0_y = buffer_1100_sdsp[97];

    auto g_y_z_0_0_0_xz_0_z = buffer_1100_sdsp[98];

    auto g_y_z_0_0_0_yy_0_x = buffer_1100_sdsp[99];

    auto g_y_z_0_0_0_yy_0_y = buffer_1100_sdsp[100];

    auto g_y_z_0_0_0_yy_0_z = buffer_1100_sdsp[101];

    auto g_y_z_0_0_0_yz_0_x = buffer_1100_sdsp[102];

    auto g_y_z_0_0_0_yz_0_y = buffer_1100_sdsp[103];

    auto g_y_z_0_0_0_yz_0_z = buffer_1100_sdsp[104];

    auto g_y_z_0_0_0_zz_0_x = buffer_1100_sdsp[105];

    auto g_y_z_0_0_0_zz_0_y = buffer_1100_sdsp[106];

    auto g_y_z_0_0_0_zz_0_z = buffer_1100_sdsp[107];

    auto g_z_x_0_0_0_xx_0_x = buffer_1100_sdsp[108];

    auto g_z_x_0_0_0_xx_0_y = buffer_1100_sdsp[109];

    auto g_z_x_0_0_0_xx_0_z = buffer_1100_sdsp[110];

    auto g_z_x_0_0_0_xy_0_x = buffer_1100_sdsp[111];

    auto g_z_x_0_0_0_xy_0_y = buffer_1100_sdsp[112];

    auto g_z_x_0_0_0_xy_0_z = buffer_1100_sdsp[113];

    auto g_z_x_0_0_0_xz_0_x = buffer_1100_sdsp[114];

    auto g_z_x_0_0_0_xz_0_y = buffer_1100_sdsp[115];

    auto g_z_x_0_0_0_xz_0_z = buffer_1100_sdsp[116];

    auto g_z_x_0_0_0_yy_0_x = buffer_1100_sdsp[117];

    auto g_z_x_0_0_0_yy_0_y = buffer_1100_sdsp[118];

    auto g_z_x_0_0_0_yy_0_z = buffer_1100_sdsp[119];

    auto g_z_x_0_0_0_yz_0_x = buffer_1100_sdsp[120];

    auto g_z_x_0_0_0_yz_0_y = buffer_1100_sdsp[121];

    auto g_z_x_0_0_0_yz_0_z = buffer_1100_sdsp[122];

    auto g_z_x_0_0_0_zz_0_x = buffer_1100_sdsp[123];

    auto g_z_x_0_0_0_zz_0_y = buffer_1100_sdsp[124];

    auto g_z_x_0_0_0_zz_0_z = buffer_1100_sdsp[125];

    auto g_z_y_0_0_0_xx_0_x = buffer_1100_sdsp[126];

    auto g_z_y_0_0_0_xx_0_y = buffer_1100_sdsp[127];

    auto g_z_y_0_0_0_xx_0_z = buffer_1100_sdsp[128];

    auto g_z_y_0_0_0_xy_0_x = buffer_1100_sdsp[129];

    auto g_z_y_0_0_0_xy_0_y = buffer_1100_sdsp[130];

    auto g_z_y_0_0_0_xy_0_z = buffer_1100_sdsp[131];

    auto g_z_y_0_0_0_xz_0_x = buffer_1100_sdsp[132];

    auto g_z_y_0_0_0_xz_0_y = buffer_1100_sdsp[133];

    auto g_z_y_0_0_0_xz_0_z = buffer_1100_sdsp[134];

    auto g_z_y_0_0_0_yy_0_x = buffer_1100_sdsp[135];

    auto g_z_y_0_0_0_yy_0_y = buffer_1100_sdsp[136];

    auto g_z_y_0_0_0_yy_0_z = buffer_1100_sdsp[137];

    auto g_z_y_0_0_0_yz_0_x = buffer_1100_sdsp[138];

    auto g_z_y_0_0_0_yz_0_y = buffer_1100_sdsp[139];

    auto g_z_y_0_0_0_yz_0_z = buffer_1100_sdsp[140];

    auto g_z_y_0_0_0_zz_0_x = buffer_1100_sdsp[141];

    auto g_z_y_0_0_0_zz_0_y = buffer_1100_sdsp[142];

    auto g_z_y_0_0_0_zz_0_z = buffer_1100_sdsp[143];

    auto g_z_z_0_0_0_xx_0_x = buffer_1100_sdsp[144];

    auto g_z_z_0_0_0_xx_0_y = buffer_1100_sdsp[145];

    auto g_z_z_0_0_0_xx_0_z = buffer_1100_sdsp[146];

    auto g_z_z_0_0_0_xy_0_x = buffer_1100_sdsp[147];

    auto g_z_z_0_0_0_xy_0_y = buffer_1100_sdsp[148];

    auto g_z_z_0_0_0_xy_0_z = buffer_1100_sdsp[149];

    auto g_z_z_0_0_0_xz_0_x = buffer_1100_sdsp[150];

    auto g_z_z_0_0_0_xz_0_y = buffer_1100_sdsp[151];

    auto g_z_z_0_0_0_xz_0_z = buffer_1100_sdsp[152];

    auto g_z_z_0_0_0_yy_0_x = buffer_1100_sdsp[153];

    auto g_z_z_0_0_0_yy_0_y = buffer_1100_sdsp[154];

    auto g_z_z_0_0_0_yy_0_z = buffer_1100_sdsp[155];

    auto g_z_z_0_0_0_yz_0_x = buffer_1100_sdsp[156];

    auto g_z_z_0_0_0_yz_0_y = buffer_1100_sdsp[157];

    auto g_z_z_0_0_0_yz_0_z = buffer_1100_sdsp[158];

    auto g_z_z_0_0_0_zz_0_x = buffer_1100_sdsp[159];

    auto g_z_z_0_0_0_zz_0_y = buffer_1100_sdsp[160];

    auto g_z_z_0_0_0_zz_0_z = buffer_1100_sdsp[161];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_x_0_0_0_xx_0_x, g_x_x_0_0_0_xx_0_y, g_x_x_0_0_0_xx_0_z, g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_x_xxx_0_x, g_x_xxx_0_y, g_x_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xx_0_x[i] = -4.0 * g_x_x_0_x[i] * a_exp + 4.0 * g_x_xxx_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_0_y[i] = -4.0 * g_x_x_0_y[i] * a_exp + 4.0 * g_x_xxx_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_xx_0_z[i] = -4.0 * g_x_x_0_z[i] * a_exp + 4.0 * g_x_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_x_0_0_0_xy_0_x, g_x_x_0_0_0_xy_0_y, g_x_x_0_0_0_xy_0_z, g_x_xxy_0_x, g_x_xxy_0_y, g_x_xxy_0_z, g_x_y_0_x, g_x_y_0_y, g_x_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xy_0_x[i] = -2.0 * g_x_y_0_x[i] * a_exp + 4.0 * g_x_xxy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_0_y[i] = -2.0 * g_x_y_0_y[i] * a_exp + 4.0 * g_x_xxy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_0_z[i] = -2.0 * g_x_y_0_z[i] * a_exp + 4.0 * g_x_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_x_0_0_0_xz_0_x, g_x_x_0_0_0_xz_0_y, g_x_x_0_0_0_xz_0_z, g_x_xxz_0_x, g_x_xxz_0_y, g_x_xxz_0_z, g_x_z_0_x, g_x_z_0_y, g_x_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xz_0_x[i] = -2.0 * g_x_z_0_x[i] * a_exp + 4.0 * g_x_xxz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_0_y[i] = -2.0 * g_x_z_0_y[i] * a_exp + 4.0 * g_x_xxz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_0_z[i] = -2.0 * g_x_z_0_z[i] * a_exp + 4.0 * g_x_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_x_0_0_0_yy_0_x, g_x_x_0_0_0_yy_0_y, g_x_x_0_0_0_yy_0_z, g_x_xyy_0_x, g_x_xyy_0_y, g_x_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yy_0_x[i] = 4.0 * g_x_xyy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_0_y[i] = 4.0 * g_x_xyy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_0_z[i] = 4.0 * g_x_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_x_0_0_0_yz_0_x, g_x_x_0_0_0_yz_0_y, g_x_x_0_0_0_yz_0_z, g_x_xyz_0_x, g_x_xyz_0_y, g_x_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_yz_0_x[i] = 4.0 * g_x_xyz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_0_y[i] = 4.0 * g_x_xyz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_0_z[i] = 4.0 * g_x_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_x_0_0_0_zz_0_x, g_x_x_0_0_0_zz_0_y, g_x_x_0_0_0_zz_0_z, g_x_xzz_0_x, g_x_xzz_0_y, g_x_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_zz_0_x[i] = 4.0 * g_x_xzz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_0_y[i] = 4.0 * g_x_xzz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_0_z[i] = 4.0 * g_x_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_xxy_0_x, g_x_xxy_0_y, g_x_xxy_0_z, g_x_y_0_0_0_xx_0_x, g_x_y_0_0_0_xx_0_y, g_x_y_0_0_0_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xx_0_x[i] = 4.0 * g_x_xxy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_0_y[i] = 4.0 * g_x_xxy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_xx_0_z[i] = 4.0 * g_x_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_x_xyy_0_x, g_x_xyy_0_y, g_x_xyy_0_z, g_x_y_0_0_0_xy_0_x, g_x_y_0_0_0_xy_0_y, g_x_y_0_0_0_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xy_0_x[i] = -2.0 * g_x_x_0_x[i] * a_exp + 4.0 * g_x_xyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_0_y[i] = -2.0 * g_x_x_0_y[i] * a_exp + 4.0 * g_x_xyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_0_z[i] = -2.0 * g_x_x_0_z[i] * a_exp + 4.0 * g_x_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_xyz_0_x, g_x_xyz_0_y, g_x_xyz_0_z, g_x_y_0_0_0_xz_0_x, g_x_y_0_0_0_xz_0_y, g_x_y_0_0_0_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xz_0_x[i] = 4.0 * g_x_xyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_0_y[i] = 4.0 * g_x_xyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_0_z[i] = 4.0 * g_x_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_y_0_0_0_yy_0_x, g_x_y_0_0_0_yy_0_y, g_x_y_0_0_0_yy_0_z, g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_x_yyy_0_x, g_x_yyy_0_y, g_x_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yy_0_x[i] = -4.0 * g_x_y_0_x[i] * a_exp + 4.0 * g_x_yyy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_0_y[i] = -4.0 * g_x_y_0_y[i] * a_exp + 4.0 * g_x_yyy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_0_z[i] = -4.0 * g_x_y_0_z[i] * a_exp + 4.0 * g_x_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_x_y_0_0_0_yz_0_x, g_x_y_0_0_0_yz_0_y, g_x_y_0_0_0_yz_0_z, g_x_yyz_0_x, g_x_yyz_0_y, g_x_yyz_0_z, g_x_z_0_x, g_x_z_0_y, g_x_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_yz_0_x[i] = -2.0 * g_x_z_0_x[i] * a_exp + 4.0 * g_x_yyz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_0_y[i] = -2.0 * g_x_z_0_y[i] * a_exp + 4.0 * g_x_yyz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_0_z[i] = -2.0 * g_x_z_0_z[i] * a_exp + 4.0 * g_x_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_x_y_0_0_0_zz_0_x, g_x_y_0_0_0_zz_0_y, g_x_y_0_0_0_zz_0_z, g_x_yzz_0_x, g_x_yzz_0_y, g_x_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_zz_0_x[i] = 4.0 * g_x_yzz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_0_y[i] = 4.0 * g_x_yzz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_0_z[i] = 4.0 * g_x_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_xxz_0_x, g_x_xxz_0_y, g_x_xxz_0_z, g_x_z_0_0_0_xx_0_x, g_x_z_0_0_0_xx_0_y, g_x_z_0_0_0_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xx_0_x[i] = 4.0 * g_x_xxz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_0_y[i] = 4.0 * g_x_xxz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_xx_0_z[i] = 4.0 * g_x_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_xyz_0_x, g_x_xyz_0_y, g_x_xyz_0_z, g_x_z_0_0_0_xy_0_x, g_x_z_0_0_0_xy_0_y, g_x_z_0_0_0_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xy_0_x[i] = 4.0 * g_x_xyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_0_y[i] = 4.0 * g_x_xyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_0_z[i] = 4.0 * g_x_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_x_xzz_0_x, g_x_xzz_0_y, g_x_xzz_0_z, g_x_z_0_0_0_xz_0_x, g_x_z_0_0_0_xz_0_y, g_x_z_0_0_0_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xz_0_x[i] = -2.0 * g_x_x_0_x[i] * a_exp + 4.0 * g_x_xzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_0_y[i] = -2.0 * g_x_x_0_y[i] * a_exp + 4.0 * g_x_xzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_0_z[i] = -2.0 * g_x_x_0_z[i] * a_exp + 4.0 * g_x_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_yyz_0_x, g_x_yyz_0_y, g_x_yyz_0_z, g_x_z_0_0_0_yy_0_x, g_x_z_0_0_0_yy_0_y, g_x_z_0_0_0_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yy_0_x[i] = 4.0 * g_x_yyz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_0_y[i] = 4.0 * g_x_yyz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_0_z[i] = 4.0 * g_x_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_x_yzz_0_x, g_x_yzz_0_y, g_x_yzz_0_z, g_x_z_0_0_0_yz_0_x, g_x_z_0_0_0_yz_0_y, g_x_z_0_0_0_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_yz_0_x[i] = -2.0 * g_x_y_0_x[i] * a_exp + 4.0 * g_x_yzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_0_y[i] = -2.0 * g_x_y_0_y[i] * a_exp + 4.0 * g_x_yzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_0_z[i] = -2.0 * g_x_y_0_z[i] * a_exp + 4.0 * g_x_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_x_z_0_0_0_zz_0_x, g_x_z_0_0_0_zz_0_y, g_x_z_0_0_0_zz_0_z, g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_x_zzz_0_x, g_x_zzz_0_y, g_x_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_zz_0_x[i] = -4.0 * g_x_z_0_x[i] * a_exp + 4.0 * g_x_zzz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_0_y[i] = -4.0 * g_x_z_0_y[i] * a_exp + 4.0 * g_x_zzz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_0_z[i] = -4.0 * g_x_z_0_z[i] * a_exp + 4.0 * g_x_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_y_x_0_0_0_xx_0_x, g_y_x_0_0_0_xx_0_y, g_y_x_0_0_0_xx_0_z, g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_y_xxx_0_x, g_y_xxx_0_y, g_y_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xx_0_x[i] = -4.0 * g_y_x_0_x[i] * a_exp + 4.0 * g_y_xxx_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_0_y[i] = -4.0 * g_y_x_0_y[i] * a_exp + 4.0 * g_y_xxx_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_xx_0_z[i] = -4.0 * g_y_x_0_z[i] * a_exp + 4.0 * g_y_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_y_x_0_0_0_xy_0_x, g_y_x_0_0_0_xy_0_y, g_y_x_0_0_0_xy_0_z, g_y_xxy_0_x, g_y_xxy_0_y, g_y_xxy_0_z, g_y_y_0_x, g_y_y_0_y, g_y_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xy_0_x[i] = -2.0 * g_y_y_0_x[i] * a_exp + 4.0 * g_y_xxy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_0_y[i] = -2.0 * g_y_y_0_y[i] * a_exp + 4.0 * g_y_xxy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_0_z[i] = -2.0 * g_y_y_0_z[i] * a_exp + 4.0 * g_y_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_y_x_0_0_0_xz_0_x, g_y_x_0_0_0_xz_0_y, g_y_x_0_0_0_xz_0_z, g_y_xxz_0_x, g_y_xxz_0_y, g_y_xxz_0_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xz_0_x[i] = -2.0 * g_y_z_0_x[i] * a_exp + 4.0 * g_y_xxz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_0_y[i] = -2.0 * g_y_z_0_y[i] * a_exp + 4.0 * g_y_xxz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_0_z[i] = -2.0 * g_y_z_0_z[i] * a_exp + 4.0 * g_y_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_y_x_0_0_0_yy_0_x, g_y_x_0_0_0_yy_0_y, g_y_x_0_0_0_yy_0_z, g_y_xyy_0_x, g_y_xyy_0_y, g_y_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yy_0_x[i] = 4.0 * g_y_xyy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_0_y[i] = 4.0 * g_y_xyy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_0_z[i] = 4.0 * g_y_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_y_x_0_0_0_yz_0_x, g_y_x_0_0_0_yz_0_y, g_y_x_0_0_0_yz_0_z, g_y_xyz_0_x, g_y_xyz_0_y, g_y_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_yz_0_x[i] = 4.0 * g_y_xyz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_0_y[i] = 4.0 * g_y_xyz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_0_z[i] = 4.0 * g_y_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_y_x_0_0_0_zz_0_x, g_y_x_0_0_0_zz_0_y, g_y_x_0_0_0_zz_0_z, g_y_xzz_0_x, g_y_xzz_0_y, g_y_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_zz_0_x[i] = 4.0 * g_y_xzz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_0_y[i] = 4.0 * g_y_xzz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_0_z[i] = 4.0 * g_y_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_y_xxy_0_x, g_y_xxy_0_y, g_y_xxy_0_z, g_y_y_0_0_0_xx_0_x, g_y_y_0_0_0_xx_0_y, g_y_y_0_0_0_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xx_0_x[i] = 4.0 * g_y_xxy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_0_y[i] = 4.0 * g_y_xxy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_xx_0_z[i] = 4.0 * g_y_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_y_xyy_0_x, g_y_xyy_0_y, g_y_xyy_0_z, g_y_y_0_0_0_xy_0_x, g_y_y_0_0_0_xy_0_y, g_y_y_0_0_0_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xy_0_x[i] = -2.0 * g_y_x_0_x[i] * a_exp + 4.0 * g_y_xyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_0_y[i] = -2.0 * g_y_x_0_y[i] * a_exp + 4.0 * g_y_xyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_0_z[i] = -2.0 * g_y_x_0_z[i] * a_exp + 4.0 * g_y_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_y_xyz_0_x, g_y_xyz_0_y, g_y_xyz_0_z, g_y_y_0_0_0_xz_0_x, g_y_y_0_0_0_xz_0_y, g_y_y_0_0_0_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xz_0_x[i] = 4.0 * g_y_xyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_0_y[i] = 4.0 * g_y_xyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_0_z[i] = 4.0 * g_y_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_y_y_0_0_0_yy_0_x, g_y_y_0_0_0_yy_0_y, g_y_y_0_0_0_yy_0_z, g_y_y_0_x, g_y_y_0_y, g_y_y_0_z, g_y_yyy_0_x, g_y_yyy_0_y, g_y_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yy_0_x[i] = -4.0 * g_y_y_0_x[i] * a_exp + 4.0 * g_y_yyy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_0_y[i] = -4.0 * g_y_y_0_y[i] * a_exp + 4.0 * g_y_yyy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_0_z[i] = -4.0 * g_y_y_0_z[i] * a_exp + 4.0 * g_y_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_y_y_0_0_0_yz_0_x, g_y_y_0_0_0_yz_0_y, g_y_y_0_0_0_yz_0_z, g_y_yyz_0_x, g_y_yyz_0_y, g_y_yyz_0_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_yz_0_x[i] = -2.0 * g_y_z_0_x[i] * a_exp + 4.0 * g_y_yyz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_0_y[i] = -2.0 * g_y_z_0_y[i] * a_exp + 4.0 * g_y_yyz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_0_z[i] = -2.0 * g_y_z_0_z[i] * a_exp + 4.0 * g_y_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_y_y_0_0_0_zz_0_x, g_y_y_0_0_0_zz_0_y, g_y_y_0_0_0_zz_0_z, g_y_yzz_0_x, g_y_yzz_0_y, g_y_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_zz_0_x[i] = 4.0 * g_y_yzz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_0_y[i] = 4.0 * g_y_yzz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_0_z[i] = 4.0 * g_y_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_y_xxz_0_x, g_y_xxz_0_y, g_y_xxz_0_z, g_y_z_0_0_0_xx_0_x, g_y_z_0_0_0_xx_0_y, g_y_z_0_0_0_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xx_0_x[i] = 4.0 * g_y_xxz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_0_y[i] = 4.0 * g_y_xxz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_xx_0_z[i] = 4.0 * g_y_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_y_xyz_0_x, g_y_xyz_0_y, g_y_xyz_0_z, g_y_z_0_0_0_xy_0_x, g_y_z_0_0_0_xy_0_y, g_y_z_0_0_0_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xy_0_x[i] = 4.0 * g_y_xyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_0_y[i] = 4.0 * g_y_xyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_0_z[i] = 4.0 * g_y_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_y_xzz_0_x, g_y_xzz_0_y, g_y_xzz_0_z, g_y_z_0_0_0_xz_0_x, g_y_z_0_0_0_xz_0_y, g_y_z_0_0_0_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xz_0_x[i] = -2.0 * g_y_x_0_x[i] * a_exp + 4.0 * g_y_xzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_0_y[i] = -2.0 * g_y_x_0_y[i] * a_exp + 4.0 * g_y_xzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_0_z[i] = -2.0 * g_y_x_0_z[i] * a_exp + 4.0 * g_y_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_y_yyz_0_x, g_y_yyz_0_y, g_y_yyz_0_z, g_y_z_0_0_0_yy_0_x, g_y_z_0_0_0_yy_0_y, g_y_z_0_0_0_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yy_0_x[i] = 4.0 * g_y_yyz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_0_y[i] = 4.0 * g_y_yyz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_0_z[i] = 4.0 * g_y_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_y_y_0_x, g_y_y_0_y, g_y_y_0_z, g_y_yzz_0_x, g_y_yzz_0_y, g_y_yzz_0_z, g_y_z_0_0_0_yz_0_x, g_y_z_0_0_0_yz_0_y, g_y_z_0_0_0_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_yz_0_x[i] = -2.0 * g_y_y_0_x[i] * a_exp + 4.0 * g_y_yzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_0_y[i] = -2.0 * g_y_y_0_y[i] * a_exp + 4.0 * g_y_yzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_0_z[i] = -2.0 * g_y_y_0_z[i] * a_exp + 4.0 * g_y_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_y_z_0_0_0_zz_0_x, g_y_z_0_0_0_zz_0_y, g_y_z_0_0_0_zz_0_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z, g_y_zzz_0_x, g_y_zzz_0_y, g_y_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_zz_0_x[i] = -4.0 * g_y_z_0_x[i] * a_exp + 4.0 * g_y_zzz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_0_y[i] = -4.0 * g_y_z_0_y[i] * a_exp + 4.0 * g_y_zzz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_0_z[i] = -4.0 * g_y_z_0_z[i] * a_exp + 4.0 * g_y_zzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_z_x_0_0_0_xx_0_x, g_z_x_0_0_0_xx_0_y, g_z_x_0_0_0_xx_0_z, g_z_x_0_x, g_z_x_0_y, g_z_x_0_z, g_z_xxx_0_x, g_z_xxx_0_y, g_z_xxx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xx_0_x[i] = -4.0 * g_z_x_0_x[i] * a_exp + 4.0 * g_z_xxx_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_0_y[i] = -4.0 * g_z_x_0_y[i] * a_exp + 4.0 * g_z_xxx_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_xx_0_z[i] = -4.0 * g_z_x_0_z[i] * a_exp + 4.0 * g_z_xxx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_z_x_0_0_0_xy_0_x, g_z_x_0_0_0_xy_0_y, g_z_x_0_0_0_xy_0_z, g_z_xxy_0_x, g_z_xxy_0_y, g_z_xxy_0_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xy_0_x[i] = -2.0 * g_z_y_0_x[i] * a_exp + 4.0 * g_z_xxy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_0_y[i] = -2.0 * g_z_y_0_y[i] * a_exp + 4.0 * g_z_xxy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_0_z[i] = -2.0 * g_z_y_0_z[i] * a_exp + 4.0 * g_z_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_z_x_0_0_0_xz_0_x, g_z_x_0_0_0_xz_0_y, g_z_x_0_0_0_xz_0_z, g_z_xxz_0_x, g_z_xxz_0_y, g_z_xxz_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xz_0_x[i] = -2.0 * g_z_z_0_x[i] * a_exp + 4.0 * g_z_xxz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_0_y[i] = -2.0 * g_z_z_0_y[i] * a_exp + 4.0 * g_z_xxz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_0_z[i] = -2.0 * g_z_z_0_z[i] * a_exp + 4.0 * g_z_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_z_x_0_0_0_yy_0_x, g_z_x_0_0_0_yy_0_y, g_z_x_0_0_0_yy_0_z, g_z_xyy_0_x, g_z_xyy_0_y, g_z_xyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yy_0_x[i] = 4.0 * g_z_xyy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_0_y[i] = 4.0 * g_z_xyy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_0_z[i] = 4.0 * g_z_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_z_x_0_0_0_yz_0_x, g_z_x_0_0_0_yz_0_y, g_z_x_0_0_0_yz_0_z, g_z_xyz_0_x, g_z_xyz_0_y, g_z_xyz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_yz_0_x[i] = 4.0 * g_z_xyz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_0_y[i] = 4.0 * g_z_xyz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_0_z[i] = 4.0 * g_z_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_z_x_0_0_0_zz_0_x, g_z_x_0_0_0_zz_0_y, g_z_x_0_0_0_zz_0_z, g_z_xzz_0_x, g_z_xzz_0_y, g_z_xzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_zz_0_x[i] = 4.0 * g_z_xzz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_0_y[i] = 4.0 * g_z_xzz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_0_z[i] = 4.0 * g_z_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_z_xxy_0_x, g_z_xxy_0_y, g_z_xxy_0_z, g_z_y_0_0_0_xx_0_x, g_z_y_0_0_0_xx_0_y, g_z_y_0_0_0_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xx_0_x[i] = 4.0 * g_z_xxy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_0_y[i] = 4.0 * g_z_xxy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_xx_0_z[i] = 4.0 * g_z_xxy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_z_x_0_x, g_z_x_0_y, g_z_x_0_z, g_z_xyy_0_x, g_z_xyy_0_y, g_z_xyy_0_z, g_z_y_0_0_0_xy_0_x, g_z_y_0_0_0_xy_0_y, g_z_y_0_0_0_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xy_0_x[i] = -2.0 * g_z_x_0_x[i] * a_exp + 4.0 * g_z_xyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_0_y[i] = -2.0 * g_z_x_0_y[i] * a_exp + 4.0 * g_z_xyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_0_z[i] = -2.0 * g_z_x_0_z[i] * a_exp + 4.0 * g_z_xyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_z_xyz_0_x, g_z_xyz_0_y, g_z_xyz_0_z, g_z_y_0_0_0_xz_0_x, g_z_y_0_0_0_xz_0_y, g_z_y_0_0_0_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xz_0_x[i] = 4.0 * g_z_xyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_0_y[i] = 4.0 * g_z_xyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_0_z[i] = 4.0 * g_z_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_z_y_0_0_0_yy_0_x, g_z_y_0_0_0_yy_0_y, g_z_y_0_0_0_yy_0_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z, g_z_yyy_0_x, g_z_yyy_0_y, g_z_yyy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yy_0_x[i] = -4.0 * g_z_y_0_x[i] * a_exp + 4.0 * g_z_yyy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_0_y[i] = -4.0 * g_z_y_0_y[i] * a_exp + 4.0 * g_z_yyy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_0_z[i] = -4.0 * g_z_y_0_z[i] * a_exp + 4.0 * g_z_yyy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_z_y_0_0_0_yz_0_x, g_z_y_0_0_0_yz_0_y, g_z_y_0_0_0_yz_0_z, g_z_yyz_0_x, g_z_yyz_0_y, g_z_yyz_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_yz_0_x[i] = -2.0 * g_z_z_0_x[i] * a_exp + 4.0 * g_z_yyz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_0_y[i] = -2.0 * g_z_z_0_y[i] * a_exp + 4.0 * g_z_yyz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_0_z[i] = -2.0 * g_z_z_0_z[i] * a_exp + 4.0 * g_z_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_z_y_0_0_0_zz_0_x, g_z_y_0_0_0_zz_0_y, g_z_y_0_0_0_zz_0_z, g_z_yzz_0_x, g_z_yzz_0_y, g_z_yzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_zz_0_x[i] = 4.0 * g_z_yzz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_0_y[i] = 4.0 * g_z_yzz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_0_z[i] = 4.0 * g_z_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_z_xxz_0_x, g_z_xxz_0_y, g_z_xxz_0_z, g_z_z_0_0_0_xx_0_x, g_z_z_0_0_0_xx_0_y, g_z_z_0_0_0_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xx_0_x[i] = 4.0 * g_z_xxz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_0_y[i] = 4.0 * g_z_xxz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_xx_0_z[i] = 4.0 * g_z_xxz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_z_xyz_0_x, g_z_xyz_0_y, g_z_xyz_0_z, g_z_z_0_0_0_xy_0_x, g_z_z_0_0_0_xy_0_y, g_z_z_0_0_0_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xy_0_x[i] = 4.0 * g_z_xyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_0_y[i] = 4.0 * g_z_xyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_0_z[i] = 4.0 * g_z_xyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_z_x_0_x, g_z_x_0_y, g_z_x_0_z, g_z_xzz_0_x, g_z_xzz_0_y, g_z_xzz_0_z, g_z_z_0_0_0_xz_0_x, g_z_z_0_0_0_xz_0_y, g_z_z_0_0_0_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xz_0_x[i] = -2.0 * g_z_x_0_x[i] * a_exp + 4.0 * g_z_xzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_0_y[i] = -2.0 * g_z_x_0_y[i] * a_exp + 4.0 * g_z_xzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_0_z[i] = -2.0 * g_z_x_0_z[i] * a_exp + 4.0 * g_z_xzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_z_yyz_0_x, g_z_yyz_0_y, g_z_yyz_0_z, g_z_z_0_0_0_yy_0_x, g_z_z_0_0_0_yy_0_y, g_z_z_0_0_0_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yy_0_x[i] = 4.0 * g_z_yyz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_0_y[i] = 4.0 * g_z_yyz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_0_z[i] = 4.0 * g_z_yyz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_z_y_0_x, g_z_y_0_y, g_z_y_0_z, g_z_yzz_0_x, g_z_yzz_0_y, g_z_yzz_0_z, g_z_z_0_0_0_yz_0_x, g_z_z_0_0_0_yz_0_y, g_z_z_0_0_0_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_yz_0_x[i] = -2.0 * g_z_y_0_x[i] * a_exp + 4.0 * g_z_yzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_0_y[i] = -2.0 * g_z_y_0_y[i] * a_exp + 4.0 * g_z_yzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_0_z[i] = -2.0 * g_z_y_0_z[i] * a_exp + 4.0 * g_z_yzz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_z_z_0_0_0_zz_0_x, g_z_z_0_0_0_zz_0_y, g_z_z_0_0_0_zz_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z, g_z_zzz_0_x, g_z_zzz_0_y, g_z_zzz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_zz_0_x[i] = -4.0 * g_z_z_0_x[i] * a_exp + 4.0 * g_z_zzz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_0_y[i] = -4.0 * g_z_z_0_y[i] * a_exp + 4.0 * g_z_zzz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_0_z[i] = -4.0 * g_z_z_0_z[i] * a_exp + 4.0 * g_z_zzz_0_z[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

