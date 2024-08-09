#include "GeomDeriv2000OfScalarForPPSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_ppsp_0(CSimdArray<double>& buffer_2000_ppsp,
                     const CSimdArray<double>& buffer_ppsp,
                     const CSimdArray<double>& buffer_fpsp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_ppsp.number_of_columns();

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

    /// Set up components of auxilary buffer : buffer_fpsp

    auto g_xxx_x_0_x = buffer_fpsp[0];

    auto g_xxx_x_0_y = buffer_fpsp[1];

    auto g_xxx_x_0_z = buffer_fpsp[2];

    auto g_xxx_y_0_x = buffer_fpsp[3];

    auto g_xxx_y_0_y = buffer_fpsp[4];

    auto g_xxx_y_0_z = buffer_fpsp[5];

    auto g_xxx_z_0_x = buffer_fpsp[6];

    auto g_xxx_z_0_y = buffer_fpsp[7];

    auto g_xxx_z_0_z = buffer_fpsp[8];

    auto g_xxy_x_0_x = buffer_fpsp[9];

    auto g_xxy_x_0_y = buffer_fpsp[10];

    auto g_xxy_x_0_z = buffer_fpsp[11];

    auto g_xxy_y_0_x = buffer_fpsp[12];

    auto g_xxy_y_0_y = buffer_fpsp[13];

    auto g_xxy_y_0_z = buffer_fpsp[14];

    auto g_xxy_z_0_x = buffer_fpsp[15];

    auto g_xxy_z_0_y = buffer_fpsp[16];

    auto g_xxy_z_0_z = buffer_fpsp[17];

    auto g_xxz_x_0_x = buffer_fpsp[18];

    auto g_xxz_x_0_y = buffer_fpsp[19];

    auto g_xxz_x_0_z = buffer_fpsp[20];

    auto g_xxz_y_0_x = buffer_fpsp[21];

    auto g_xxz_y_0_y = buffer_fpsp[22];

    auto g_xxz_y_0_z = buffer_fpsp[23];

    auto g_xxz_z_0_x = buffer_fpsp[24];

    auto g_xxz_z_0_y = buffer_fpsp[25];

    auto g_xxz_z_0_z = buffer_fpsp[26];

    auto g_xyy_x_0_x = buffer_fpsp[27];

    auto g_xyy_x_0_y = buffer_fpsp[28];

    auto g_xyy_x_0_z = buffer_fpsp[29];

    auto g_xyy_y_0_x = buffer_fpsp[30];

    auto g_xyy_y_0_y = buffer_fpsp[31];

    auto g_xyy_y_0_z = buffer_fpsp[32];

    auto g_xyy_z_0_x = buffer_fpsp[33];

    auto g_xyy_z_0_y = buffer_fpsp[34];

    auto g_xyy_z_0_z = buffer_fpsp[35];

    auto g_xyz_x_0_x = buffer_fpsp[36];

    auto g_xyz_x_0_y = buffer_fpsp[37];

    auto g_xyz_x_0_z = buffer_fpsp[38];

    auto g_xyz_y_0_x = buffer_fpsp[39];

    auto g_xyz_y_0_y = buffer_fpsp[40];

    auto g_xyz_y_0_z = buffer_fpsp[41];

    auto g_xyz_z_0_x = buffer_fpsp[42];

    auto g_xyz_z_0_y = buffer_fpsp[43];

    auto g_xyz_z_0_z = buffer_fpsp[44];

    auto g_xzz_x_0_x = buffer_fpsp[45];

    auto g_xzz_x_0_y = buffer_fpsp[46];

    auto g_xzz_x_0_z = buffer_fpsp[47];

    auto g_xzz_y_0_x = buffer_fpsp[48];

    auto g_xzz_y_0_y = buffer_fpsp[49];

    auto g_xzz_y_0_z = buffer_fpsp[50];

    auto g_xzz_z_0_x = buffer_fpsp[51];

    auto g_xzz_z_0_y = buffer_fpsp[52];

    auto g_xzz_z_0_z = buffer_fpsp[53];

    auto g_yyy_x_0_x = buffer_fpsp[54];

    auto g_yyy_x_0_y = buffer_fpsp[55];

    auto g_yyy_x_0_z = buffer_fpsp[56];

    auto g_yyy_y_0_x = buffer_fpsp[57];

    auto g_yyy_y_0_y = buffer_fpsp[58];

    auto g_yyy_y_0_z = buffer_fpsp[59];

    auto g_yyy_z_0_x = buffer_fpsp[60];

    auto g_yyy_z_0_y = buffer_fpsp[61];

    auto g_yyy_z_0_z = buffer_fpsp[62];

    auto g_yyz_x_0_x = buffer_fpsp[63];

    auto g_yyz_x_0_y = buffer_fpsp[64];

    auto g_yyz_x_0_z = buffer_fpsp[65];

    auto g_yyz_y_0_x = buffer_fpsp[66];

    auto g_yyz_y_0_y = buffer_fpsp[67];

    auto g_yyz_y_0_z = buffer_fpsp[68];

    auto g_yyz_z_0_x = buffer_fpsp[69];

    auto g_yyz_z_0_y = buffer_fpsp[70];

    auto g_yyz_z_0_z = buffer_fpsp[71];

    auto g_yzz_x_0_x = buffer_fpsp[72];

    auto g_yzz_x_0_y = buffer_fpsp[73];

    auto g_yzz_x_0_z = buffer_fpsp[74];

    auto g_yzz_y_0_x = buffer_fpsp[75];

    auto g_yzz_y_0_y = buffer_fpsp[76];

    auto g_yzz_y_0_z = buffer_fpsp[77];

    auto g_yzz_z_0_x = buffer_fpsp[78];

    auto g_yzz_z_0_y = buffer_fpsp[79];

    auto g_yzz_z_0_z = buffer_fpsp[80];

    auto g_zzz_x_0_x = buffer_fpsp[81];

    auto g_zzz_x_0_y = buffer_fpsp[82];

    auto g_zzz_x_0_z = buffer_fpsp[83];

    auto g_zzz_y_0_x = buffer_fpsp[84];

    auto g_zzz_y_0_y = buffer_fpsp[85];

    auto g_zzz_y_0_z = buffer_fpsp[86];

    auto g_zzz_z_0_x = buffer_fpsp[87];

    auto g_zzz_z_0_y = buffer_fpsp[88];

    auto g_zzz_z_0_z = buffer_fpsp[89];

    /// Set up components of integrals buffer : buffer_2000_ppsp

    auto g_xx_0_0_0_x_x_0_x = buffer_2000_ppsp[0];

    auto g_xx_0_0_0_x_x_0_y = buffer_2000_ppsp[1];

    auto g_xx_0_0_0_x_x_0_z = buffer_2000_ppsp[2];

    auto g_xx_0_0_0_x_y_0_x = buffer_2000_ppsp[3];

    auto g_xx_0_0_0_x_y_0_y = buffer_2000_ppsp[4];

    auto g_xx_0_0_0_x_y_0_z = buffer_2000_ppsp[5];

    auto g_xx_0_0_0_x_z_0_x = buffer_2000_ppsp[6];

    auto g_xx_0_0_0_x_z_0_y = buffer_2000_ppsp[7];

    auto g_xx_0_0_0_x_z_0_z = buffer_2000_ppsp[8];

    auto g_xx_0_0_0_y_x_0_x = buffer_2000_ppsp[9];

    auto g_xx_0_0_0_y_x_0_y = buffer_2000_ppsp[10];

    auto g_xx_0_0_0_y_x_0_z = buffer_2000_ppsp[11];

    auto g_xx_0_0_0_y_y_0_x = buffer_2000_ppsp[12];

    auto g_xx_0_0_0_y_y_0_y = buffer_2000_ppsp[13];

    auto g_xx_0_0_0_y_y_0_z = buffer_2000_ppsp[14];

    auto g_xx_0_0_0_y_z_0_x = buffer_2000_ppsp[15];

    auto g_xx_0_0_0_y_z_0_y = buffer_2000_ppsp[16];

    auto g_xx_0_0_0_y_z_0_z = buffer_2000_ppsp[17];

    auto g_xx_0_0_0_z_x_0_x = buffer_2000_ppsp[18];

    auto g_xx_0_0_0_z_x_0_y = buffer_2000_ppsp[19];

    auto g_xx_0_0_0_z_x_0_z = buffer_2000_ppsp[20];

    auto g_xx_0_0_0_z_y_0_x = buffer_2000_ppsp[21];

    auto g_xx_0_0_0_z_y_0_y = buffer_2000_ppsp[22];

    auto g_xx_0_0_0_z_y_0_z = buffer_2000_ppsp[23];

    auto g_xx_0_0_0_z_z_0_x = buffer_2000_ppsp[24];

    auto g_xx_0_0_0_z_z_0_y = buffer_2000_ppsp[25];

    auto g_xx_0_0_0_z_z_0_z = buffer_2000_ppsp[26];

    auto g_xy_0_0_0_x_x_0_x = buffer_2000_ppsp[27];

    auto g_xy_0_0_0_x_x_0_y = buffer_2000_ppsp[28];

    auto g_xy_0_0_0_x_x_0_z = buffer_2000_ppsp[29];

    auto g_xy_0_0_0_x_y_0_x = buffer_2000_ppsp[30];

    auto g_xy_0_0_0_x_y_0_y = buffer_2000_ppsp[31];

    auto g_xy_0_0_0_x_y_0_z = buffer_2000_ppsp[32];

    auto g_xy_0_0_0_x_z_0_x = buffer_2000_ppsp[33];

    auto g_xy_0_0_0_x_z_0_y = buffer_2000_ppsp[34];

    auto g_xy_0_0_0_x_z_0_z = buffer_2000_ppsp[35];

    auto g_xy_0_0_0_y_x_0_x = buffer_2000_ppsp[36];

    auto g_xy_0_0_0_y_x_0_y = buffer_2000_ppsp[37];

    auto g_xy_0_0_0_y_x_0_z = buffer_2000_ppsp[38];

    auto g_xy_0_0_0_y_y_0_x = buffer_2000_ppsp[39];

    auto g_xy_0_0_0_y_y_0_y = buffer_2000_ppsp[40];

    auto g_xy_0_0_0_y_y_0_z = buffer_2000_ppsp[41];

    auto g_xy_0_0_0_y_z_0_x = buffer_2000_ppsp[42];

    auto g_xy_0_0_0_y_z_0_y = buffer_2000_ppsp[43];

    auto g_xy_0_0_0_y_z_0_z = buffer_2000_ppsp[44];

    auto g_xy_0_0_0_z_x_0_x = buffer_2000_ppsp[45];

    auto g_xy_0_0_0_z_x_0_y = buffer_2000_ppsp[46];

    auto g_xy_0_0_0_z_x_0_z = buffer_2000_ppsp[47];

    auto g_xy_0_0_0_z_y_0_x = buffer_2000_ppsp[48];

    auto g_xy_0_0_0_z_y_0_y = buffer_2000_ppsp[49];

    auto g_xy_0_0_0_z_y_0_z = buffer_2000_ppsp[50];

    auto g_xy_0_0_0_z_z_0_x = buffer_2000_ppsp[51];

    auto g_xy_0_0_0_z_z_0_y = buffer_2000_ppsp[52];

    auto g_xy_0_0_0_z_z_0_z = buffer_2000_ppsp[53];

    auto g_xz_0_0_0_x_x_0_x = buffer_2000_ppsp[54];

    auto g_xz_0_0_0_x_x_0_y = buffer_2000_ppsp[55];

    auto g_xz_0_0_0_x_x_0_z = buffer_2000_ppsp[56];

    auto g_xz_0_0_0_x_y_0_x = buffer_2000_ppsp[57];

    auto g_xz_0_0_0_x_y_0_y = buffer_2000_ppsp[58];

    auto g_xz_0_0_0_x_y_0_z = buffer_2000_ppsp[59];

    auto g_xz_0_0_0_x_z_0_x = buffer_2000_ppsp[60];

    auto g_xz_0_0_0_x_z_0_y = buffer_2000_ppsp[61];

    auto g_xz_0_0_0_x_z_0_z = buffer_2000_ppsp[62];

    auto g_xz_0_0_0_y_x_0_x = buffer_2000_ppsp[63];

    auto g_xz_0_0_0_y_x_0_y = buffer_2000_ppsp[64];

    auto g_xz_0_0_0_y_x_0_z = buffer_2000_ppsp[65];

    auto g_xz_0_0_0_y_y_0_x = buffer_2000_ppsp[66];

    auto g_xz_0_0_0_y_y_0_y = buffer_2000_ppsp[67];

    auto g_xz_0_0_0_y_y_0_z = buffer_2000_ppsp[68];

    auto g_xz_0_0_0_y_z_0_x = buffer_2000_ppsp[69];

    auto g_xz_0_0_0_y_z_0_y = buffer_2000_ppsp[70];

    auto g_xz_0_0_0_y_z_0_z = buffer_2000_ppsp[71];

    auto g_xz_0_0_0_z_x_0_x = buffer_2000_ppsp[72];

    auto g_xz_0_0_0_z_x_0_y = buffer_2000_ppsp[73];

    auto g_xz_0_0_0_z_x_0_z = buffer_2000_ppsp[74];

    auto g_xz_0_0_0_z_y_0_x = buffer_2000_ppsp[75];

    auto g_xz_0_0_0_z_y_0_y = buffer_2000_ppsp[76];

    auto g_xz_0_0_0_z_y_0_z = buffer_2000_ppsp[77];

    auto g_xz_0_0_0_z_z_0_x = buffer_2000_ppsp[78];

    auto g_xz_0_0_0_z_z_0_y = buffer_2000_ppsp[79];

    auto g_xz_0_0_0_z_z_0_z = buffer_2000_ppsp[80];

    auto g_yy_0_0_0_x_x_0_x = buffer_2000_ppsp[81];

    auto g_yy_0_0_0_x_x_0_y = buffer_2000_ppsp[82];

    auto g_yy_0_0_0_x_x_0_z = buffer_2000_ppsp[83];

    auto g_yy_0_0_0_x_y_0_x = buffer_2000_ppsp[84];

    auto g_yy_0_0_0_x_y_0_y = buffer_2000_ppsp[85];

    auto g_yy_0_0_0_x_y_0_z = buffer_2000_ppsp[86];

    auto g_yy_0_0_0_x_z_0_x = buffer_2000_ppsp[87];

    auto g_yy_0_0_0_x_z_0_y = buffer_2000_ppsp[88];

    auto g_yy_0_0_0_x_z_0_z = buffer_2000_ppsp[89];

    auto g_yy_0_0_0_y_x_0_x = buffer_2000_ppsp[90];

    auto g_yy_0_0_0_y_x_0_y = buffer_2000_ppsp[91];

    auto g_yy_0_0_0_y_x_0_z = buffer_2000_ppsp[92];

    auto g_yy_0_0_0_y_y_0_x = buffer_2000_ppsp[93];

    auto g_yy_0_0_0_y_y_0_y = buffer_2000_ppsp[94];

    auto g_yy_0_0_0_y_y_0_z = buffer_2000_ppsp[95];

    auto g_yy_0_0_0_y_z_0_x = buffer_2000_ppsp[96];

    auto g_yy_0_0_0_y_z_0_y = buffer_2000_ppsp[97];

    auto g_yy_0_0_0_y_z_0_z = buffer_2000_ppsp[98];

    auto g_yy_0_0_0_z_x_0_x = buffer_2000_ppsp[99];

    auto g_yy_0_0_0_z_x_0_y = buffer_2000_ppsp[100];

    auto g_yy_0_0_0_z_x_0_z = buffer_2000_ppsp[101];

    auto g_yy_0_0_0_z_y_0_x = buffer_2000_ppsp[102];

    auto g_yy_0_0_0_z_y_0_y = buffer_2000_ppsp[103];

    auto g_yy_0_0_0_z_y_0_z = buffer_2000_ppsp[104];

    auto g_yy_0_0_0_z_z_0_x = buffer_2000_ppsp[105];

    auto g_yy_0_0_0_z_z_0_y = buffer_2000_ppsp[106];

    auto g_yy_0_0_0_z_z_0_z = buffer_2000_ppsp[107];

    auto g_yz_0_0_0_x_x_0_x = buffer_2000_ppsp[108];

    auto g_yz_0_0_0_x_x_0_y = buffer_2000_ppsp[109];

    auto g_yz_0_0_0_x_x_0_z = buffer_2000_ppsp[110];

    auto g_yz_0_0_0_x_y_0_x = buffer_2000_ppsp[111];

    auto g_yz_0_0_0_x_y_0_y = buffer_2000_ppsp[112];

    auto g_yz_0_0_0_x_y_0_z = buffer_2000_ppsp[113];

    auto g_yz_0_0_0_x_z_0_x = buffer_2000_ppsp[114];

    auto g_yz_0_0_0_x_z_0_y = buffer_2000_ppsp[115];

    auto g_yz_0_0_0_x_z_0_z = buffer_2000_ppsp[116];

    auto g_yz_0_0_0_y_x_0_x = buffer_2000_ppsp[117];

    auto g_yz_0_0_0_y_x_0_y = buffer_2000_ppsp[118];

    auto g_yz_0_0_0_y_x_0_z = buffer_2000_ppsp[119];

    auto g_yz_0_0_0_y_y_0_x = buffer_2000_ppsp[120];

    auto g_yz_0_0_0_y_y_0_y = buffer_2000_ppsp[121];

    auto g_yz_0_0_0_y_y_0_z = buffer_2000_ppsp[122];

    auto g_yz_0_0_0_y_z_0_x = buffer_2000_ppsp[123];

    auto g_yz_0_0_0_y_z_0_y = buffer_2000_ppsp[124];

    auto g_yz_0_0_0_y_z_0_z = buffer_2000_ppsp[125];

    auto g_yz_0_0_0_z_x_0_x = buffer_2000_ppsp[126];

    auto g_yz_0_0_0_z_x_0_y = buffer_2000_ppsp[127];

    auto g_yz_0_0_0_z_x_0_z = buffer_2000_ppsp[128];

    auto g_yz_0_0_0_z_y_0_x = buffer_2000_ppsp[129];

    auto g_yz_0_0_0_z_y_0_y = buffer_2000_ppsp[130];

    auto g_yz_0_0_0_z_y_0_z = buffer_2000_ppsp[131];

    auto g_yz_0_0_0_z_z_0_x = buffer_2000_ppsp[132];

    auto g_yz_0_0_0_z_z_0_y = buffer_2000_ppsp[133];

    auto g_yz_0_0_0_z_z_0_z = buffer_2000_ppsp[134];

    auto g_zz_0_0_0_x_x_0_x = buffer_2000_ppsp[135];

    auto g_zz_0_0_0_x_x_0_y = buffer_2000_ppsp[136];

    auto g_zz_0_0_0_x_x_0_z = buffer_2000_ppsp[137];

    auto g_zz_0_0_0_x_y_0_x = buffer_2000_ppsp[138];

    auto g_zz_0_0_0_x_y_0_y = buffer_2000_ppsp[139];

    auto g_zz_0_0_0_x_y_0_z = buffer_2000_ppsp[140];

    auto g_zz_0_0_0_x_z_0_x = buffer_2000_ppsp[141];

    auto g_zz_0_0_0_x_z_0_y = buffer_2000_ppsp[142];

    auto g_zz_0_0_0_x_z_0_z = buffer_2000_ppsp[143];

    auto g_zz_0_0_0_y_x_0_x = buffer_2000_ppsp[144];

    auto g_zz_0_0_0_y_x_0_y = buffer_2000_ppsp[145];

    auto g_zz_0_0_0_y_x_0_z = buffer_2000_ppsp[146];

    auto g_zz_0_0_0_y_y_0_x = buffer_2000_ppsp[147];

    auto g_zz_0_0_0_y_y_0_y = buffer_2000_ppsp[148];

    auto g_zz_0_0_0_y_y_0_z = buffer_2000_ppsp[149];

    auto g_zz_0_0_0_y_z_0_x = buffer_2000_ppsp[150];

    auto g_zz_0_0_0_y_z_0_y = buffer_2000_ppsp[151];

    auto g_zz_0_0_0_y_z_0_z = buffer_2000_ppsp[152];

    auto g_zz_0_0_0_z_x_0_x = buffer_2000_ppsp[153];

    auto g_zz_0_0_0_z_x_0_y = buffer_2000_ppsp[154];

    auto g_zz_0_0_0_z_x_0_z = buffer_2000_ppsp[155];

    auto g_zz_0_0_0_z_y_0_x = buffer_2000_ppsp[156];

    auto g_zz_0_0_0_z_y_0_y = buffer_2000_ppsp[157];

    auto g_zz_0_0_0_z_y_0_z = buffer_2000_ppsp[158];

    auto g_zz_0_0_0_z_z_0_x = buffer_2000_ppsp[159];

    auto g_zz_0_0_0_z_z_0_y = buffer_2000_ppsp[160];

    auto g_zz_0_0_0_z_z_0_z = buffer_2000_ppsp[161];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_xx_0_0_0_x_x_0_x, g_xx_0_0_0_x_x_0_y, g_xx_0_0_0_x_x_0_z, g_xxx_x_0_x, g_xxx_x_0_y, g_xxx_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_x_0_x[i] = -6.0 * g_x_x_0_x[i] * a_exp + 4.0 * g_xxx_x_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_0_y[i] = -6.0 * g_x_x_0_y[i] * a_exp + 4.0 * g_xxx_x_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_x_0_z[i] = -6.0 * g_x_x_0_z[i] * a_exp + 4.0 * g_xxx_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_xx_0_0_0_x_y_0_x, g_xx_0_0_0_x_y_0_y, g_xx_0_0_0_x_y_0_z, g_xxx_y_0_x, g_xxx_y_0_y, g_xxx_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_y_0_x[i] = -6.0 * g_x_y_0_x[i] * a_exp + 4.0 * g_xxx_y_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_0_y[i] = -6.0 * g_x_y_0_y[i] * a_exp + 4.0 * g_xxx_y_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_0_z[i] = -6.0 * g_x_y_0_z[i] * a_exp + 4.0 * g_xxx_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_xx_0_0_0_x_z_0_x, g_xx_0_0_0_x_z_0_y, g_xx_0_0_0_x_z_0_z, g_xxx_z_0_x, g_xxx_z_0_y, g_xxx_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_z_0_x[i] = -6.0 * g_x_z_0_x[i] * a_exp + 4.0 * g_xxx_z_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_0_y[i] = -6.0 * g_x_z_0_y[i] * a_exp + 4.0 * g_xxx_z_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_0_z[i] = -6.0 * g_x_z_0_z[i] * a_exp + 4.0 * g_xxx_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_xx_0_0_0_y_x_0_x, g_xx_0_0_0_y_x_0_y, g_xx_0_0_0_y_x_0_z, g_xxy_x_0_x, g_xxy_x_0_y, g_xxy_x_0_z, g_y_x_0_x, g_y_x_0_y, g_y_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_x_0_x[i] = -2.0 * g_y_x_0_x[i] * a_exp + 4.0 * g_xxy_x_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_0_y[i] = -2.0 * g_y_x_0_y[i] * a_exp + 4.0 * g_xxy_x_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_x_0_z[i] = -2.0 * g_y_x_0_z[i] * a_exp + 4.0 * g_xxy_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_xx_0_0_0_y_y_0_x, g_xx_0_0_0_y_y_0_y, g_xx_0_0_0_y_y_0_z, g_xxy_y_0_x, g_xxy_y_0_y, g_xxy_y_0_z, g_y_y_0_x, g_y_y_0_y, g_y_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_y_0_x[i] = -2.0 * g_y_y_0_x[i] * a_exp + 4.0 * g_xxy_y_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_0_y[i] = -2.0 * g_y_y_0_y[i] * a_exp + 4.0 * g_xxy_y_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_0_z[i] = -2.0 * g_y_y_0_z[i] * a_exp + 4.0 * g_xxy_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_xx_0_0_0_y_z_0_x, g_xx_0_0_0_y_z_0_y, g_xx_0_0_0_y_z_0_z, g_xxy_z_0_x, g_xxy_z_0_y, g_xxy_z_0_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_z_0_x[i] = -2.0 * g_y_z_0_x[i] * a_exp + 4.0 * g_xxy_z_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_0_y[i] = -2.0 * g_y_z_0_y[i] * a_exp + 4.0 * g_xxy_z_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_0_z[i] = -2.0 * g_y_z_0_z[i] * a_exp + 4.0 * g_xxy_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_xx_0_0_0_z_x_0_x, g_xx_0_0_0_z_x_0_y, g_xx_0_0_0_z_x_0_z, g_xxz_x_0_x, g_xxz_x_0_y, g_xxz_x_0_z, g_z_x_0_x, g_z_x_0_y, g_z_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_x_0_x[i] = -2.0 * g_z_x_0_x[i] * a_exp + 4.0 * g_xxz_x_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_0_y[i] = -2.0 * g_z_x_0_y[i] * a_exp + 4.0 * g_xxz_x_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_x_0_z[i] = -2.0 * g_z_x_0_z[i] * a_exp + 4.0 * g_xxz_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_xx_0_0_0_z_y_0_x, g_xx_0_0_0_z_y_0_y, g_xx_0_0_0_z_y_0_z, g_xxz_y_0_x, g_xxz_y_0_y, g_xxz_y_0_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_y_0_x[i] = -2.0 * g_z_y_0_x[i] * a_exp + 4.0 * g_xxz_y_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_0_y[i] = -2.0 * g_z_y_0_y[i] * a_exp + 4.0 * g_xxz_y_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_0_z[i] = -2.0 * g_z_y_0_z[i] * a_exp + 4.0 * g_xxz_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_xx_0_0_0_z_z_0_x, g_xx_0_0_0_z_z_0_y, g_xx_0_0_0_z_z_0_z, g_xxz_z_0_x, g_xxz_z_0_y, g_xxz_z_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_z_0_x[i] = -2.0 * g_z_z_0_x[i] * a_exp + 4.0 * g_xxz_z_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_0_y[i] = -2.0 * g_z_z_0_y[i] * a_exp + 4.0 * g_xxz_z_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_0_z[i] = -2.0 * g_z_z_0_z[i] * a_exp + 4.0 * g_xxz_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_xxy_x_0_x, g_xxy_x_0_y, g_xxy_x_0_z, g_xy_0_0_0_x_x_0_x, g_xy_0_0_0_x_x_0_y, g_xy_0_0_0_x_x_0_z, g_y_x_0_x, g_y_x_0_y, g_y_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_x_0_x[i] = -2.0 * g_y_x_0_x[i] * a_exp + 4.0 * g_xxy_x_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_0_y[i] = -2.0 * g_y_x_0_y[i] * a_exp + 4.0 * g_xxy_x_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_x_0_z[i] = -2.0 * g_y_x_0_z[i] * a_exp + 4.0 * g_xxy_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_xxy_y_0_x, g_xxy_y_0_y, g_xxy_y_0_z, g_xy_0_0_0_x_y_0_x, g_xy_0_0_0_x_y_0_y, g_xy_0_0_0_x_y_0_z, g_y_y_0_x, g_y_y_0_y, g_y_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_y_0_x[i] = -2.0 * g_y_y_0_x[i] * a_exp + 4.0 * g_xxy_y_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_0_y[i] = -2.0 * g_y_y_0_y[i] * a_exp + 4.0 * g_xxy_y_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_0_z[i] = -2.0 * g_y_y_0_z[i] * a_exp + 4.0 * g_xxy_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_xxy_z_0_x, g_xxy_z_0_y, g_xxy_z_0_z, g_xy_0_0_0_x_z_0_x, g_xy_0_0_0_x_z_0_y, g_xy_0_0_0_x_z_0_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_z_0_x[i] = -2.0 * g_y_z_0_x[i] * a_exp + 4.0 * g_xxy_z_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_0_y[i] = -2.0 * g_y_z_0_y[i] * a_exp + 4.0 * g_xxy_z_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_0_z[i] = -2.0 * g_y_z_0_z[i] * a_exp + 4.0 * g_xxy_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_xy_0_0_0_y_x_0_x, g_xy_0_0_0_y_x_0_y, g_xy_0_0_0_y_x_0_z, g_xyy_x_0_x, g_xyy_x_0_y, g_xyy_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_x_0_x[i] = -2.0 * g_x_x_0_x[i] * a_exp + 4.0 * g_xyy_x_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_0_y[i] = -2.0 * g_x_x_0_y[i] * a_exp + 4.0 * g_xyy_x_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_x_0_z[i] = -2.0 * g_x_x_0_z[i] * a_exp + 4.0 * g_xyy_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_xy_0_0_0_y_y_0_x, g_xy_0_0_0_y_y_0_y, g_xy_0_0_0_y_y_0_z, g_xyy_y_0_x, g_xyy_y_0_y, g_xyy_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_y_0_x[i] = -2.0 * g_x_y_0_x[i] * a_exp + 4.0 * g_xyy_y_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_0_y[i] = -2.0 * g_x_y_0_y[i] * a_exp + 4.0 * g_xyy_y_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_0_z[i] = -2.0 * g_x_y_0_z[i] * a_exp + 4.0 * g_xyy_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_xy_0_0_0_y_z_0_x, g_xy_0_0_0_y_z_0_y, g_xy_0_0_0_y_z_0_z, g_xyy_z_0_x, g_xyy_z_0_y, g_xyy_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_z_0_x[i] = -2.0 * g_x_z_0_x[i] * a_exp + 4.0 * g_xyy_z_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_0_y[i] = -2.0 * g_x_z_0_y[i] * a_exp + 4.0 * g_xyy_z_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_0_z[i] = -2.0 * g_x_z_0_z[i] * a_exp + 4.0 * g_xyy_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_xy_0_0_0_z_x_0_x, g_xy_0_0_0_z_x_0_y, g_xy_0_0_0_z_x_0_z, g_xyz_x_0_x, g_xyz_x_0_y, g_xyz_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_x_0_x[i] = 4.0 * g_xyz_x_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_0_y[i] = 4.0 * g_xyz_x_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_x_0_z[i] = 4.0 * g_xyz_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_xy_0_0_0_z_y_0_x, g_xy_0_0_0_z_y_0_y, g_xy_0_0_0_z_y_0_z, g_xyz_y_0_x, g_xyz_y_0_y, g_xyz_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_y_0_x[i] = 4.0 * g_xyz_y_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_0_y[i] = 4.0 * g_xyz_y_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_0_z[i] = 4.0 * g_xyz_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_xy_0_0_0_z_z_0_x, g_xy_0_0_0_z_z_0_y, g_xy_0_0_0_z_z_0_z, g_xyz_z_0_x, g_xyz_z_0_y, g_xyz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_z_0_x[i] = 4.0 * g_xyz_z_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_0_y[i] = 4.0 * g_xyz_z_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_0_z[i] = 4.0 * g_xyz_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_xxz_x_0_x, g_xxz_x_0_y, g_xxz_x_0_z, g_xz_0_0_0_x_x_0_x, g_xz_0_0_0_x_x_0_y, g_xz_0_0_0_x_x_0_z, g_z_x_0_x, g_z_x_0_y, g_z_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_x_0_x[i] = -2.0 * g_z_x_0_x[i] * a_exp + 4.0 * g_xxz_x_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_0_y[i] = -2.0 * g_z_x_0_y[i] * a_exp + 4.0 * g_xxz_x_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_x_0_z[i] = -2.0 * g_z_x_0_z[i] * a_exp + 4.0 * g_xxz_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_xxz_y_0_x, g_xxz_y_0_y, g_xxz_y_0_z, g_xz_0_0_0_x_y_0_x, g_xz_0_0_0_x_y_0_y, g_xz_0_0_0_x_y_0_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_y_0_x[i] = -2.0 * g_z_y_0_x[i] * a_exp + 4.0 * g_xxz_y_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_0_y[i] = -2.0 * g_z_y_0_y[i] * a_exp + 4.0 * g_xxz_y_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_0_z[i] = -2.0 * g_z_y_0_z[i] * a_exp + 4.0 * g_xxz_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_xxz_z_0_x, g_xxz_z_0_y, g_xxz_z_0_z, g_xz_0_0_0_x_z_0_x, g_xz_0_0_0_x_z_0_y, g_xz_0_0_0_x_z_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_z_0_x[i] = -2.0 * g_z_z_0_x[i] * a_exp + 4.0 * g_xxz_z_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_0_y[i] = -2.0 * g_z_z_0_y[i] * a_exp + 4.0 * g_xxz_z_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_0_z[i] = -2.0 * g_z_z_0_z[i] * a_exp + 4.0 * g_xxz_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_xyz_x_0_x, g_xyz_x_0_y, g_xyz_x_0_z, g_xz_0_0_0_y_x_0_x, g_xz_0_0_0_y_x_0_y, g_xz_0_0_0_y_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_x_0_x[i] = 4.0 * g_xyz_x_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_0_y[i] = 4.0 * g_xyz_x_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_x_0_z[i] = 4.0 * g_xyz_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_xyz_y_0_x, g_xyz_y_0_y, g_xyz_y_0_z, g_xz_0_0_0_y_y_0_x, g_xz_0_0_0_y_y_0_y, g_xz_0_0_0_y_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_y_0_x[i] = 4.0 * g_xyz_y_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_0_y[i] = 4.0 * g_xyz_y_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_0_z[i] = 4.0 * g_xyz_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_xyz_z_0_x, g_xyz_z_0_y, g_xyz_z_0_z, g_xz_0_0_0_y_z_0_x, g_xz_0_0_0_y_z_0_y, g_xz_0_0_0_y_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_z_0_x[i] = 4.0 * g_xyz_z_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_0_y[i] = 4.0 * g_xyz_z_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_0_z[i] = 4.0 * g_xyz_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_xz_0_0_0_z_x_0_x, g_xz_0_0_0_z_x_0_y, g_xz_0_0_0_z_x_0_z, g_xzz_x_0_x, g_xzz_x_0_y, g_xzz_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_x_0_x[i] = -2.0 * g_x_x_0_x[i] * a_exp + 4.0 * g_xzz_x_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_0_y[i] = -2.0 * g_x_x_0_y[i] * a_exp + 4.0 * g_xzz_x_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_x_0_z[i] = -2.0 * g_x_x_0_z[i] * a_exp + 4.0 * g_xzz_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_xz_0_0_0_z_y_0_x, g_xz_0_0_0_z_y_0_y, g_xz_0_0_0_z_y_0_z, g_xzz_y_0_x, g_xzz_y_0_y, g_xzz_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_y_0_x[i] = -2.0 * g_x_y_0_x[i] * a_exp + 4.0 * g_xzz_y_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_0_y[i] = -2.0 * g_x_y_0_y[i] * a_exp + 4.0 * g_xzz_y_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_0_z[i] = -2.0 * g_x_y_0_z[i] * a_exp + 4.0 * g_xzz_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_xz_0_0_0_z_z_0_x, g_xz_0_0_0_z_z_0_y, g_xz_0_0_0_z_z_0_z, g_xzz_z_0_x, g_xzz_z_0_y, g_xzz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_z_0_x[i] = -2.0 * g_x_z_0_x[i] * a_exp + 4.0 * g_xzz_z_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_0_y[i] = -2.0 * g_x_z_0_y[i] * a_exp + 4.0 * g_xzz_z_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_0_z[i] = -2.0 * g_x_z_0_z[i] * a_exp + 4.0 * g_xzz_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (81-84)

    #pragma omp simd aligned(g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_xyy_x_0_x, g_xyy_x_0_y, g_xyy_x_0_z, g_yy_0_0_0_x_x_0_x, g_yy_0_0_0_x_x_0_y, g_yy_0_0_0_x_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_x_0_x[i] = -2.0 * g_x_x_0_x[i] * a_exp + 4.0 * g_xyy_x_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_0_y[i] = -2.0 * g_x_x_0_y[i] * a_exp + 4.0 * g_xyy_x_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_x_0_z[i] = -2.0 * g_x_x_0_z[i] * a_exp + 4.0 * g_xyy_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (84-87)

    #pragma omp simd aligned(g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_xyy_y_0_x, g_xyy_y_0_y, g_xyy_y_0_z, g_yy_0_0_0_x_y_0_x, g_yy_0_0_0_x_y_0_y, g_yy_0_0_0_x_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_y_0_x[i] = -2.0 * g_x_y_0_x[i] * a_exp + 4.0 * g_xyy_y_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_0_y[i] = -2.0 * g_x_y_0_y[i] * a_exp + 4.0 * g_xyy_y_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_0_z[i] = -2.0 * g_x_y_0_z[i] * a_exp + 4.0 * g_xyy_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (87-90)

    #pragma omp simd aligned(g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_xyy_z_0_x, g_xyy_z_0_y, g_xyy_z_0_z, g_yy_0_0_0_x_z_0_x, g_yy_0_0_0_x_z_0_y, g_yy_0_0_0_x_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_z_0_x[i] = -2.0 * g_x_z_0_x[i] * a_exp + 4.0 * g_xyy_z_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_0_y[i] = -2.0 * g_x_z_0_y[i] * a_exp + 4.0 * g_xyy_z_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_0_z[i] = -2.0 * g_x_z_0_z[i] * a_exp + 4.0 * g_xyy_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (90-93)

    #pragma omp simd aligned(g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_yy_0_0_0_y_x_0_x, g_yy_0_0_0_y_x_0_y, g_yy_0_0_0_y_x_0_z, g_yyy_x_0_x, g_yyy_x_0_y, g_yyy_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_x_0_x[i] = -6.0 * g_y_x_0_x[i] * a_exp + 4.0 * g_yyy_x_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_0_y[i] = -6.0 * g_y_x_0_y[i] * a_exp + 4.0 * g_yyy_x_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_x_0_z[i] = -6.0 * g_y_x_0_z[i] * a_exp + 4.0 * g_yyy_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (93-96)

    #pragma omp simd aligned(g_y_y_0_x, g_y_y_0_y, g_y_y_0_z, g_yy_0_0_0_y_y_0_x, g_yy_0_0_0_y_y_0_y, g_yy_0_0_0_y_y_0_z, g_yyy_y_0_x, g_yyy_y_0_y, g_yyy_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_y_0_x[i] = -6.0 * g_y_y_0_x[i] * a_exp + 4.0 * g_yyy_y_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_0_y[i] = -6.0 * g_y_y_0_y[i] * a_exp + 4.0 * g_yyy_y_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_0_z[i] = -6.0 * g_y_y_0_z[i] * a_exp + 4.0 * g_yyy_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (96-99)

    #pragma omp simd aligned(g_y_z_0_x, g_y_z_0_y, g_y_z_0_z, g_yy_0_0_0_y_z_0_x, g_yy_0_0_0_y_z_0_y, g_yy_0_0_0_y_z_0_z, g_yyy_z_0_x, g_yyy_z_0_y, g_yyy_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_z_0_x[i] = -6.0 * g_y_z_0_x[i] * a_exp + 4.0 * g_yyy_z_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_0_y[i] = -6.0 * g_y_z_0_y[i] * a_exp + 4.0 * g_yyy_z_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_0_z[i] = -6.0 * g_y_z_0_z[i] * a_exp + 4.0 * g_yyy_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (99-102)

    #pragma omp simd aligned(g_yy_0_0_0_z_x_0_x, g_yy_0_0_0_z_x_0_y, g_yy_0_0_0_z_x_0_z, g_yyz_x_0_x, g_yyz_x_0_y, g_yyz_x_0_z, g_z_x_0_x, g_z_x_0_y, g_z_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_x_0_x[i] = -2.0 * g_z_x_0_x[i] * a_exp + 4.0 * g_yyz_x_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_0_y[i] = -2.0 * g_z_x_0_y[i] * a_exp + 4.0 * g_yyz_x_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_x_0_z[i] = -2.0 * g_z_x_0_z[i] * a_exp + 4.0 * g_yyz_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (102-105)

    #pragma omp simd aligned(g_yy_0_0_0_z_y_0_x, g_yy_0_0_0_z_y_0_y, g_yy_0_0_0_z_y_0_z, g_yyz_y_0_x, g_yyz_y_0_y, g_yyz_y_0_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_y_0_x[i] = -2.0 * g_z_y_0_x[i] * a_exp + 4.0 * g_yyz_y_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_0_y[i] = -2.0 * g_z_y_0_y[i] * a_exp + 4.0 * g_yyz_y_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_0_z[i] = -2.0 * g_z_y_0_z[i] * a_exp + 4.0 * g_yyz_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (105-108)

    #pragma omp simd aligned(g_yy_0_0_0_z_z_0_x, g_yy_0_0_0_z_z_0_y, g_yy_0_0_0_z_z_0_z, g_yyz_z_0_x, g_yyz_z_0_y, g_yyz_z_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_z_0_x[i] = -2.0 * g_z_z_0_x[i] * a_exp + 4.0 * g_yyz_z_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_0_y[i] = -2.0 * g_z_z_0_y[i] * a_exp + 4.0 * g_yyz_z_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_0_z[i] = -2.0 * g_z_z_0_z[i] * a_exp + 4.0 * g_yyz_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (108-111)

    #pragma omp simd aligned(g_xyz_x_0_x, g_xyz_x_0_y, g_xyz_x_0_z, g_yz_0_0_0_x_x_0_x, g_yz_0_0_0_x_x_0_y, g_yz_0_0_0_x_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_x_0_x[i] = 4.0 * g_xyz_x_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_0_y[i] = 4.0 * g_xyz_x_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_x_0_z[i] = 4.0 * g_xyz_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (111-114)

    #pragma omp simd aligned(g_xyz_y_0_x, g_xyz_y_0_y, g_xyz_y_0_z, g_yz_0_0_0_x_y_0_x, g_yz_0_0_0_x_y_0_y, g_yz_0_0_0_x_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_y_0_x[i] = 4.0 * g_xyz_y_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_0_y[i] = 4.0 * g_xyz_y_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_0_z[i] = 4.0 * g_xyz_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (114-117)

    #pragma omp simd aligned(g_xyz_z_0_x, g_xyz_z_0_y, g_xyz_z_0_z, g_yz_0_0_0_x_z_0_x, g_yz_0_0_0_x_z_0_y, g_yz_0_0_0_x_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_z_0_x[i] = 4.0 * g_xyz_z_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_0_y[i] = 4.0 * g_xyz_z_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_0_z[i] = 4.0 * g_xyz_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (117-120)

    #pragma omp simd aligned(g_yyz_x_0_x, g_yyz_x_0_y, g_yyz_x_0_z, g_yz_0_0_0_y_x_0_x, g_yz_0_0_0_y_x_0_y, g_yz_0_0_0_y_x_0_z, g_z_x_0_x, g_z_x_0_y, g_z_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_x_0_x[i] = -2.0 * g_z_x_0_x[i] * a_exp + 4.0 * g_yyz_x_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_0_y[i] = -2.0 * g_z_x_0_y[i] * a_exp + 4.0 * g_yyz_x_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_x_0_z[i] = -2.0 * g_z_x_0_z[i] * a_exp + 4.0 * g_yyz_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (120-123)

    #pragma omp simd aligned(g_yyz_y_0_x, g_yyz_y_0_y, g_yyz_y_0_z, g_yz_0_0_0_y_y_0_x, g_yz_0_0_0_y_y_0_y, g_yz_0_0_0_y_y_0_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_y_0_x[i] = -2.0 * g_z_y_0_x[i] * a_exp + 4.0 * g_yyz_y_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_0_y[i] = -2.0 * g_z_y_0_y[i] * a_exp + 4.0 * g_yyz_y_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_0_z[i] = -2.0 * g_z_y_0_z[i] * a_exp + 4.0 * g_yyz_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (123-126)

    #pragma omp simd aligned(g_yyz_z_0_x, g_yyz_z_0_y, g_yyz_z_0_z, g_yz_0_0_0_y_z_0_x, g_yz_0_0_0_y_z_0_y, g_yz_0_0_0_y_z_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_z_0_x[i] = -2.0 * g_z_z_0_x[i] * a_exp + 4.0 * g_yyz_z_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_0_y[i] = -2.0 * g_z_z_0_y[i] * a_exp + 4.0 * g_yyz_z_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_0_z[i] = -2.0 * g_z_z_0_z[i] * a_exp + 4.0 * g_yyz_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (126-129)

    #pragma omp simd aligned(g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_yz_0_0_0_z_x_0_x, g_yz_0_0_0_z_x_0_y, g_yz_0_0_0_z_x_0_z, g_yzz_x_0_x, g_yzz_x_0_y, g_yzz_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_x_0_x[i] = -2.0 * g_y_x_0_x[i] * a_exp + 4.0 * g_yzz_x_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_0_y[i] = -2.0 * g_y_x_0_y[i] * a_exp + 4.0 * g_yzz_x_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_x_0_z[i] = -2.0 * g_y_x_0_z[i] * a_exp + 4.0 * g_yzz_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (129-132)

    #pragma omp simd aligned(g_y_y_0_x, g_y_y_0_y, g_y_y_0_z, g_yz_0_0_0_z_y_0_x, g_yz_0_0_0_z_y_0_y, g_yz_0_0_0_z_y_0_z, g_yzz_y_0_x, g_yzz_y_0_y, g_yzz_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_y_0_x[i] = -2.0 * g_y_y_0_x[i] * a_exp + 4.0 * g_yzz_y_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_0_y[i] = -2.0 * g_y_y_0_y[i] * a_exp + 4.0 * g_yzz_y_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_0_z[i] = -2.0 * g_y_y_0_z[i] * a_exp + 4.0 * g_yzz_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (132-135)

    #pragma omp simd aligned(g_y_z_0_x, g_y_z_0_y, g_y_z_0_z, g_yz_0_0_0_z_z_0_x, g_yz_0_0_0_z_z_0_y, g_yz_0_0_0_z_z_0_z, g_yzz_z_0_x, g_yzz_z_0_y, g_yzz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_z_0_x[i] = -2.0 * g_y_z_0_x[i] * a_exp + 4.0 * g_yzz_z_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_0_y[i] = -2.0 * g_y_z_0_y[i] * a_exp + 4.0 * g_yzz_z_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_0_z[i] = -2.0 * g_y_z_0_z[i] * a_exp + 4.0 * g_yzz_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (135-138)

    #pragma omp simd aligned(g_x_x_0_x, g_x_x_0_y, g_x_x_0_z, g_xzz_x_0_x, g_xzz_x_0_y, g_xzz_x_0_z, g_zz_0_0_0_x_x_0_x, g_zz_0_0_0_x_x_0_y, g_zz_0_0_0_x_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_x_0_x[i] = -2.0 * g_x_x_0_x[i] * a_exp + 4.0 * g_xzz_x_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_0_y[i] = -2.0 * g_x_x_0_y[i] * a_exp + 4.0 * g_xzz_x_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_x_0_z[i] = -2.0 * g_x_x_0_z[i] * a_exp + 4.0 * g_xzz_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (138-141)

    #pragma omp simd aligned(g_x_y_0_x, g_x_y_0_y, g_x_y_0_z, g_xzz_y_0_x, g_xzz_y_0_y, g_xzz_y_0_z, g_zz_0_0_0_x_y_0_x, g_zz_0_0_0_x_y_0_y, g_zz_0_0_0_x_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_y_0_x[i] = -2.0 * g_x_y_0_x[i] * a_exp + 4.0 * g_xzz_y_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_0_y[i] = -2.0 * g_x_y_0_y[i] * a_exp + 4.0 * g_xzz_y_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_0_z[i] = -2.0 * g_x_y_0_z[i] * a_exp + 4.0 * g_xzz_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (141-144)

    #pragma omp simd aligned(g_x_z_0_x, g_x_z_0_y, g_x_z_0_z, g_xzz_z_0_x, g_xzz_z_0_y, g_xzz_z_0_z, g_zz_0_0_0_x_z_0_x, g_zz_0_0_0_x_z_0_y, g_zz_0_0_0_x_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_z_0_x[i] = -2.0 * g_x_z_0_x[i] * a_exp + 4.0 * g_xzz_z_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_0_y[i] = -2.0 * g_x_z_0_y[i] * a_exp + 4.0 * g_xzz_z_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_0_z[i] = -2.0 * g_x_z_0_z[i] * a_exp + 4.0 * g_xzz_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (144-147)

    #pragma omp simd aligned(g_y_x_0_x, g_y_x_0_y, g_y_x_0_z, g_yzz_x_0_x, g_yzz_x_0_y, g_yzz_x_0_z, g_zz_0_0_0_y_x_0_x, g_zz_0_0_0_y_x_0_y, g_zz_0_0_0_y_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_x_0_x[i] = -2.0 * g_y_x_0_x[i] * a_exp + 4.0 * g_yzz_x_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_0_y[i] = -2.0 * g_y_x_0_y[i] * a_exp + 4.0 * g_yzz_x_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_x_0_z[i] = -2.0 * g_y_x_0_z[i] * a_exp + 4.0 * g_yzz_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (147-150)

    #pragma omp simd aligned(g_y_y_0_x, g_y_y_0_y, g_y_y_0_z, g_yzz_y_0_x, g_yzz_y_0_y, g_yzz_y_0_z, g_zz_0_0_0_y_y_0_x, g_zz_0_0_0_y_y_0_y, g_zz_0_0_0_y_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_y_0_x[i] = -2.0 * g_y_y_0_x[i] * a_exp + 4.0 * g_yzz_y_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_0_y[i] = -2.0 * g_y_y_0_y[i] * a_exp + 4.0 * g_yzz_y_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_0_z[i] = -2.0 * g_y_y_0_z[i] * a_exp + 4.0 * g_yzz_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (150-153)

    #pragma omp simd aligned(g_y_z_0_x, g_y_z_0_y, g_y_z_0_z, g_yzz_z_0_x, g_yzz_z_0_y, g_yzz_z_0_z, g_zz_0_0_0_y_z_0_x, g_zz_0_0_0_y_z_0_y, g_zz_0_0_0_y_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_z_0_x[i] = -2.0 * g_y_z_0_x[i] * a_exp + 4.0 * g_yzz_z_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_0_y[i] = -2.0 * g_y_z_0_y[i] * a_exp + 4.0 * g_yzz_z_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_0_z[i] = -2.0 * g_y_z_0_z[i] * a_exp + 4.0 * g_yzz_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (153-156)

    #pragma omp simd aligned(g_z_x_0_x, g_z_x_0_y, g_z_x_0_z, g_zz_0_0_0_z_x_0_x, g_zz_0_0_0_z_x_0_y, g_zz_0_0_0_z_x_0_z, g_zzz_x_0_x, g_zzz_x_0_y, g_zzz_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_x_0_x[i] = -6.0 * g_z_x_0_x[i] * a_exp + 4.0 * g_zzz_x_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_0_y[i] = -6.0 * g_z_x_0_y[i] * a_exp + 4.0 * g_zzz_x_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_x_0_z[i] = -6.0 * g_z_x_0_z[i] * a_exp + 4.0 * g_zzz_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (156-159)

    #pragma omp simd aligned(g_z_y_0_x, g_z_y_0_y, g_z_y_0_z, g_zz_0_0_0_z_y_0_x, g_zz_0_0_0_z_y_0_y, g_zz_0_0_0_z_y_0_z, g_zzz_y_0_x, g_zzz_y_0_y, g_zzz_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_y_0_x[i] = -6.0 * g_z_y_0_x[i] * a_exp + 4.0 * g_zzz_y_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_0_y[i] = -6.0 * g_z_y_0_y[i] * a_exp + 4.0 * g_zzz_y_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_0_z[i] = -6.0 * g_z_y_0_z[i] * a_exp + 4.0 * g_zzz_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (159-162)

    #pragma omp simd aligned(g_z_z_0_x, g_z_z_0_y, g_z_z_0_z, g_zz_0_0_0_z_z_0_x, g_zz_0_0_0_z_z_0_y, g_zz_0_0_0_z_z_0_z, g_zzz_z_0_x, g_zzz_z_0_y, g_zzz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_z_0_x[i] = -6.0 * g_z_z_0_x[i] * a_exp + 4.0 * g_zzz_z_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_0_y[i] = -6.0 * g_z_z_0_y[i] * a_exp + 4.0 * g_zzz_z_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_0_z[i] = -6.0 * g_z_z_0_z[i] * a_exp + 4.0 * g_zzz_z_0_z[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

