#include "GeomDeriv2000OfScalarForPDSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_pdss_0(CSimdArray<double>& buffer_2000_pdss,
                     const CSimdArray<double>& buffer_pdss,
                     const CSimdArray<double>& buffer_fdss,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_pdss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pdss

    auto g_x_xx_0_0 = buffer_pdss[0];

    auto g_x_xy_0_0 = buffer_pdss[1];

    auto g_x_xz_0_0 = buffer_pdss[2];

    auto g_x_yy_0_0 = buffer_pdss[3];

    auto g_x_yz_0_0 = buffer_pdss[4];

    auto g_x_zz_0_0 = buffer_pdss[5];

    auto g_y_xx_0_0 = buffer_pdss[6];

    auto g_y_xy_0_0 = buffer_pdss[7];

    auto g_y_xz_0_0 = buffer_pdss[8];

    auto g_y_yy_0_0 = buffer_pdss[9];

    auto g_y_yz_0_0 = buffer_pdss[10];

    auto g_y_zz_0_0 = buffer_pdss[11];

    auto g_z_xx_0_0 = buffer_pdss[12];

    auto g_z_xy_0_0 = buffer_pdss[13];

    auto g_z_xz_0_0 = buffer_pdss[14];

    auto g_z_yy_0_0 = buffer_pdss[15];

    auto g_z_yz_0_0 = buffer_pdss[16];

    auto g_z_zz_0_0 = buffer_pdss[17];

    /// Set up components of auxilary buffer : buffer_fdss

    auto g_xxx_xx_0_0 = buffer_fdss[0];

    auto g_xxx_xy_0_0 = buffer_fdss[1];

    auto g_xxx_xz_0_0 = buffer_fdss[2];

    auto g_xxx_yy_0_0 = buffer_fdss[3];

    auto g_xxx_yz_0_0 = buffer_fdss[4];

    auto g_xxx_zz_0_0 = buffer_fdss[5];

    auto g_xxy_xx_0_0 = buffer_fdss[6];

    auto g_xxy_xy_0_0 = buffer_fdss[7];

    auto g_xxy_xz_0_0 = buffer_fdss[8];

    auto g_xxy_yy_0_0 = buffer_fdss[9];

    auto g_xxy_yz_0_0 = buffer_fdss[10];

    auto g_xxy_zz_0_0 = buffer_fdss[11];

    auto g_xxz_xx_0_0 = buffer_fdss[12];

    auto g_xxz_xy_0_0 = buffer_fdss[13];

    auto g_xxz_xz_0_0 = buffer_fdss[14];

    auto g_xxz_yy_0_0 = buffer_fdss[15];

    auto g_xxz_yz_0_0 = buffer_fdss[16];

    auto g_xxz_zz_0_0 = buffer_fdss[17];

    auto g_xyy_xx_0_0 = buffer_fdss[18];

    auto g_xyy_xy_0_0 = buffer_fdss[19];

    auto g_xyy_xz_0_0 = buffer_fdss[20];

    auto g_xyy_yy_0_0 = buffer_fdss[21];

    auto g_xyy_yz_0_0 = buffer_fdss[22];

    auto g_xyy_zz_0_0 = buffer_fdss[23];

    auto g_xyz_xx_0_0 = buffer_fdss[24];

    auto g_xyz_xy_0_0 = buffer_fdss[25];

    auto g_xyz_xz_0_0 = buffer_fdss[26];

    auto g_xyz_yy_0_0 = buffer_fdss[27];

    auto g_xyz_yz_0_0 = buffer_fdss[28];

    auto g_xyz_zz_0_0 = buffer_fdss[29];

    auto g_xzz_xx_0_0 = buffer_fdss[30];

    auto g_xzz_xy_0_0 = buffer_fdss[31];

    auto g_xzz_xz_0_0 = buffer_fdss[32];

    auto g_xzz_yy_0_0 = buffer_fdss[33];

    auto g_xzz_yz_0_0 = buffer_fdss[34];

    auto g_xzz_zz_0_0 = buffer_fdss[35];

    auto g_yyy_xx_0_0 = buffer_fdss[36];

    auto g_yyy_xy_0_0 = buffer_fdss[37];

    auto g_yyy_xz_0_0 = buffer_fdss[38];

    auto g_yyy_yy_0_0 = buffer_fdss[39];

    auto g_yyy_yz_0_0 = buffer_fdss[40];

    auto g_yyy_zz_0_0 = buffer_fdss[41];

    auto g_yyz_xx_0_0 = buffer_fdss[42];

    auto g_yyz_xy_0_0 = buffer_fdss[43];

    auto g_yyz_xz_0_0 = buffer_fdss[44];

    auto g_yyz_yy_0_0 = buffer_fdss[45];

    auto g_yyz_yz_0_0 = buffer_fdss[46];

    auto g_yyz_zz_0_0 = buffer_fdss[47];

    auto g_yzz_xx_0_0 = buffer_fdss[48];

    auto g_yzz_xy_0_0 = buffer_fdss[49];

    auto g_yzz_xz_0_0 = buffer_fdss[50];

    auto g_yzz_yy_0_0 = buffer_fdss[51];

    auto g_yzz_yz_0_0 = buffer_fdss[52];

    auto g_yzz_zz_0_0 = buffer_fdss[53];

    auto g_zzz_xx_0_0 = buffer_fdss[54];

    auto g_zzz_xy_0_0 = buffer_fdss[55];

    auto g_zzz_xz_0_0 = buffer_fdss[56];

    auto g_zzz_yy_0_0 = buffer_fdss[57];

    auto g_zzz_yz_0_0 = buffer_fdss[58];

    auto g_zzz_zz_0_0 = buffer_fdss[59];

    /// Set up components of integrals buffer : buffer_2000_pdss

    auto g_xx_0_0_0_x_xx_0_0 = buffer_2000_pdss[0];

    auto g_xx_0_0_0_x_xy_0_0 = buffer_2000_pdss[1];

    auto g_xx_0_0_0_x_xz_0_0 = buffer_2000_pdss[2];

    auto g_xx_0_0_0_x_yy_0_0 = buffer_2000_pdss[3];

    auto g_xx_0_0_0_x_yz_0_0 = buffer_2000_pdss[4];

    auto g_xx_0_0_0_x_zz_0_0 = buffer_2000_pdss[5];

    auto g_xx_0_0_0_y_xx_0_0 = buffer_2000_pdss[6];

    auto g_xx_0_0_0_y_xy_0_0 = buffer_2000_pdss[7];

    auto g_xx_0_0_0_y_xz_0_0 = buffer_2000_pdss[8];

    auto g_xx_0_0_0_y_yy_0_0 = buffer_2000_pdss[9];

    auto g_xx_0_0_0_y_yz_0_0 = buffer_2000_pdss[10];

    auto g_xx_0_0_0_y_zz_0_0 = buffer_2000_pdss[11];

    auto g_xx_0_0_0_z_xx_0_0 = buffer_2000_pdss[12];

    auto g_xx_0_0_0_z_xy_0_0 = buffer_2000_pdss[13];

    auto g_xx_0_0_0_z_xz_0_0 = buffer_2000_pdss[14];

    auto g_xx_0_0_0_z_yy_0_0 = buffer_2000_pdss[15];

    auto g_xx_0_0_0_z_yz_0_0 = buffer_2000_pdss[16];

    auto g_xx_0_0_0_z_zz_0_0 = buffer_2000_pdss[17];

    auto g_xy_0_0_0_x_xx_0_0 = buffer_2000_pdss[18];

    auto g_xy_0_0_0_x_xy_0_0 = buffer_2000_pdss[19];

    auto g_xy_0_0_0_x_xz_0_0 = buffer_2000_pdss[20];

    auto g_xy_0_0_0_x_yy_0_0 = buffer_2000_pdss[21];

    auto g_xy_0_0_0_x_yz_0_0 = buffer_2000_pdss[22];

    auto g_xy_0_0_0_x_zz_0_0 = buffer_2000_pdss[23];

    auto g_xy_0_0_0_y_xx_0_0 = buffer_2000_pdss[24];

    auto g_xy_0_0_0_y_xy_0_0 = buffer_2000_pdss[25];

    auto g_xy_0_0_0_y_xz_0_0 = buffer_2000_pdss[26];

    auto g_xy_0_0_0_y_yy_0_0 = buffer_2000_pdss[27];

    auto g_xy_0_0_0_y_yz_0_0 = buffer_2000_pdss[28];

    auto g_xy_0_0_0_y_zz_0_0 = buffer_2000_pdss[29];

    auto g_xy_0_0_0_z_xx_0_0 = buffer_2000_pdss[30];

    auto g_xy_0_0_0_z_xy_0_0 = buffer_2000_pdss[31];

    auto g_xy_0_0_0_z_xz_0_0 = buffer_2000_pdss[32];

    auto g_xy_0_0_0_z_yy_0_0 = buffer_2000_pdss[33];

    auto g_xy_0_0_0_z_yz_0_0 = buffer_2000_pdss[34];

    auto g_xy_0_0_0_z_zz_0_0 = buffer_2000_pdss[35];

    auto g_xz_0_0_0_x_xx_0_0 = buffer_2000_pdss[36];

    auto g_xz_0_0_0_x_xy_0_0 = buffer_2000_pdss[37];

    auto g_xz_0_0_0_x_xz_0_0 = buffer_2000_pdss[38];

    auto g_xz_0_0_0_x_yy_0_0 = buffer_2000_pdss[39];

    auto g_xz_0_0_0_x_yz_0_0 = buffer_2000_pdss[40];

    auto g_xz_0_0_0_x_zz_0_0 = buffer_2000_pdss[41];

    auto g_xz_0_0_0_y_xx_0_0 = buffer_2000_pdss[42];

    auto g_xz_0_0_0_y_xy_0_0 = buffer_2000_pdss[43];

    auto g_xz_0_0_0_y_xz_0_0 = buffer_2000_pdss[44];

    auto g_xz_0_0_0_y_yy_0_0 = buffer_2000_pdss[45];

    auto g_xz_0_0_0_y_yz_0_0 = buffer_2000_pdss[46];

    auto g_xz_0_0_0_y_zz_0_0 = buffer_2000_pdss[47];

    auto g_xz_0_0_0_z_xx_0_0 = buffer_2000_pdss[48];

    auto g_xz_0_0_0_z_xy_0_0 = buffer_2000_pdss[49];

    auto g_xz_0_0_0_z_xz_0_0 = buffer_2000_pdss[50];

    auto g_xz_0_0_0_z_yy_0_0 = buffer_2000_pdss[51];

    auto g_xz_0_0_0_z_yz_0_0 = buffer_2000_pdss[52];

    auto g_xz_0_0_0_z_zz_0_0 = buffer_2000_pdss[53];

    auto g_yy_0_0_0_x_xx_0_0 = buffer_2000_pdss[54];

    auto g_yy_0_0_0_x_xy_0_0 = buffer_2000_pdss[55];

    auto g_yy_0_0_0_x_xz_0_0 = buffer_2000_pdss[56];

    auto g_yy_0_0_0_x_yy_0_0 = buffer_2000_pdss[57];

    auto g_yy_0_0_0_x_yz_0_0 = buffer_2000_pdss[58];

    auto g_yy_0_0_0_x_zz_0_0 = buffer_2000_pdss[59];

    auto g_yy_0_0_0_y_xx_0_0 = buffer_2000_pdss[60];

    auto g_yy_0_0_0_y_xy_0_0 = buffer_2000_pdss[61];

    auto g_yy_0_0_0_y_xz_0_0 = buffer_2000_pdss[62];

    auto g_yy_0_0_0_y_yy_0_0 = buffer_2000_pdss[63];

    auto g_yy_0_0_0_y_yz_0_0 = buffer_2000_pdss[64];

    auto g_yy_0_0_0_y_zz_0_0 = buffer_2000_pdss[65];

    auto g_yy_0_0_0_z_xx_0_0 = buffer_2000_pdss[66];

    auto g_yy_0_0_0_z_xy_0_0 = buffer_2000_pdss[67];

    auto g_yy_0_0_0_z_xz_0_0 = buffer_2000_pdss[68];

    auto g_yy_0_0_0_z_yy_0_0 = buffer_2000_pdss[69];

    auto g_yy_0_0_0_z_yz_0_0 = buffer_2000_pdss[70];

    auto g_yy_0_0_0_z_zz_0_0 = buffer_2000_pdss[71];

    auto g_yz_0_0_0_x_xx_0_0 = buffer_2000_pdss[72];

    auto g_yz_0_0_0_x_xy_0_0 = buffer_2000_pdss[73];

    auto g_yz_0_0_0_x_xz_0_0 = buffer_2000_pdss[74];

    auto g_yz_0_0_0_x_yy_0_0 = buffer_2000_pdss[75];

    auto g_yz_0_0_0_x_yz_0_0 = buffer_2000_pdss[76];

    auto g_yz_0_0_0_x_zz_0_0 = buffer_2000_pdss[77];

    auto g_yz_0_0_0_y_xx_0_0 = buffer_2000_pdss[78];

    auto g_yz_0_0_0_y_xy_0_0 = buffer_2000_pdss[79];

    auto g_yz_0_0_0_y_xz_0_0 = buffer_2000_pdss[80];

    auto g_yz_0_0_0_y_yy_0_0 = buffer_2000_pdss[81];

    auto g_yz_0_0_0_y_yz_0_0 = buffer_2000_pdss[82];

    auto g_yz_0_0_0_y_zz_0_0 = buffer_2000_pdss[83];

    auto g_yz_0_0_0_z_xx_0_0 = buffer_2000_pdss[84];

    auto g_yz_0_0_0_z_xy_0_0 = buffer_2000_pdss[85];

    auto g_yz_0_0_0_z_xz_0_0 = buffer_2000_pdss[86];

    auto g_yz_0_0_0_z_yy_0_0 = buffer_2000_pdss[87];

    auto g_yz_0_0_0_z_yz_0_0 = buffer_2000_pdss[88];

    auto g_yz_0_0_0_z_zz_0_0 = buffer_2000_pdss[89];

    auto g_zz_0_0_0_x_xx_0_0 = buffer_2000_pdss[90];

    auto g_zz_0_0_0_x_xy_0_0 = buffer_2000_pdss[91];

    auto g_zz_0_0_0_x_xz_0_0 = buffer_2000_pdss[92];

    auto g_zz_0_0_0_x_yy_0_0 = buffer_2000_pdss[93];

    auto g_zz_0_0_0_x_yz_0_0 = buffer_2000_pdss[94];

    auto g_zz_0_0_0_x_zz_0_0 = buffer_2000_pdss[95];

    auto g_zz_0_0_0_y_xx_0_0 = buffer_2000_pdss[96];

    auto g_zz_0_0_0_y_xy_0_0 = buffer_2000_pdss[97];

    auto g_zz_0_0_0_y_xz_0_0 = buffer_2000_pdss[98];

    auto g_zz_0_0_0_y_yy_0_0 = buffer_2000_pdss[99];

    auto g_zz_0_0_0_y_yz_0_0 = buffer_2000_pdss[100];

    auto g_zz_0_0_0_y_zz_0_0 = buffer_2000_pdss[101];

    auto g_zz_0_0_0_z_xx_0_0 = buffer_2000_pdss[102];

    auto g_zz_0_0_0_z_xy_0_0 = buffer_2000_pdss[103];

    auto g_zz_0_0_0_z_xz_0_0 = buffer_2000_pdss[104];

    auto g_zz_0_0_0_z_yy_0_0 = buffer_2000_pdss[105];

    auto g_zz_0_0_0_z_yz_0_0 = buffer_2000_pdss[106];

    auto g_zz_0_0_0_z_zz_0_0 = buffer_2000_pdss[107];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_xx_0_0, g_x_xy_0_0, g_x_xz_0_0, g_x_yy_0_0, g_x_yz_0_0, g_x_zz_0_0, g_xx_0_0_0_x_xx_0_0, g_xx_0_0_0_x_xy_0_0, g_xx_0_0_0_x_xz_0_0, g_xx_0_0_0_x_yy_0_0, g_xx_0_0_0_x_yz_0_0, g_xx_0_0_0_x_zz_0_0, g_xxx_xx_0_0, g_xxx_xy_0_0, g_xxx_xz_0_0, g_xxx_yy_0_0, g_xxx_yz_0_0, g_xxx_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_xx_0_0[i] = -6.0 * g_x_xx_0_0[i] * a_exp + 4.0 * g_xxx_xx_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xy_0_0[i] = -6.0 * g_x_xy_0_0[i] * a_exp + 4.0 * g_xxx_xy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_x_xz_0_0[i] = -6.0 * g_x_xz_0_0[i] * a_exp + 4.0 * g_xxx_xz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yy_0_0[i] = -6.0 * g_x_yy_0_0[i] * a_exp + 4.0 * g_xxx_yy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_x_yz_0_0[i] = -6.0 * g_x_yz_0_0[i] * a_exp + 4.0 * g_xxx_yz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_x_zz_0_0[i] = -6.0 * g_x_zz_0_0[i] * a_exp + 4.0 * g_xxx_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_xx_0_0_0_y_xx_0_0, g_xx_0_0_0_y_xy_0_0, g_xx_0_0_0_y_xz_0_0, g_xx_0_0_0_y_yy_0_0, g_xx_0_0_0_y_yz_0_0, g_xx_0_0_0_y_zz_0_0, g_xxy_xx_0_0, g_xxy_xy_0_0, g_xxy_xz_0_0, g_xxy_yy_0_0, g_xxy_yz_0_0, g_xxy_zz_0_0, g_y_xx_0_0, g_y_xy_0_0, g_y_xz_0_0, g_y_yy_0_0, g_y_yz_0_0, g_y_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_xx_0_0[i] = -2.0 * g_y_xx_0_0[i] * a_exp + 4.0 * g_xxy_xx_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xy_0_0[i] = -2.0 * g_y_xy_0_0[i] * a_exp + 4.0 * g_xxy_xy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_y_xz_0_0[i] = -2.0 * g_y_xz_0_0[i] * a_exp + 4.0 * g_xxy_xz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yy_0_0[i] = -2.0 * g_y_yy_0_0[i] * a_exp + 4.0 * g_xxy_yy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_y_yz_0_0[i] = -2.0 * g_y_yz_0_0[i] * a_exp + 4.0 * g_xxy_yz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_y_zz_0_0[i] = -2.0 * g_y_zz_0_0[i] * a_exp + 4.0 * g_xxy_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_xx_0_0_0_z_xx_0_0, g_xx_0_0_0_z_xy_0_0, g_xx_0_0_0_z_xz_0_0, g_xx_0_0_0_z_yy_0_0, g_xx_0_0_0_z_yz_0_0, g_xx_0_0_0_z_zz_0_0, g_xxz_xx_0_0, g_xxz_xy_0_0, g_xxz_xz_0_0, g_xxz_yy_0_0, g_xxz_yz_0_0, g_xxz_zz_0_0, g_z_xx_0_0, g_z_xy_0_0, g_z_xz_0_0, g_z_yy_0_0, g_z_yz_0_0, g_z_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_xx_0_0[i] = -2.0 * g_z_xx_0_0[i] * a_exp + 4.0 * g_xxz_xx_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xy_0_0[i] = -2.0 * g_z_xy_0_0[i] * a_exp + 4.0 * g_xxz_xy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_z_xz_0_0[i] = -2.0 * g_z_xz_0_0[i] * a_exp + 4.0 * g_xxz_xz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yy_0_0[i] = -2.0 * g_z_yy_0_0[i] * a_exp + 4.0 * g_xxz_yy_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_z_yz_0_0[i] = -2.0 * g_z_yz_0_0[i] * a_exp + 4.0 * g_xxz_yz_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_z_zz_0_0[i] = -2.0 * g_z_zz_0_0[i] * a_exp + 4.0 * g_xxz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_xxy_xx_0_0, g_xxy_xy_0_0, g_xxy_xz_0_0, g_xxy_yy_0_0, g_xxy_yz_0_0, g_xxy_zz_0_0, g_xy_0_0_0_x_xx_0_0, g_xy_0_0_0_x_xy_0_0, g_xy_0_0_0_x_xz_0_0, g_xy_0_0_0_x_yy_0_0, g_xy_0_0_0_x_yz_0_0, g_xy_0_0_0_x_zz_0_0, g_y_xx_0_0, g_y_xy_0_0, g_y_xz_0_0, g_y_yy_0_0, g_y_yz_0_0, g_y_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_xx_0_0[i] = -2.0 * g_y_xx_0_0[i] * a_exp + 4.0 * g_xxy_xx_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xy_0_0[i] = -2.0 * g_y_xy_0_0[i] * a_exp + 4.0 * g_xxy_xy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_x_xz_0_0[i] = -2.0 * g_y_xz_0_0[i] * a_exp + 4.0 * g_xxy_xz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yy_0_0[i] = -2.0 * g_y_yy_0_0[i] * a_exp + 4.0 * g_xxy_yy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_x_yz_0_0[i] = -2.0 * g_y_yz_0_0[i] * a_exp + 4.0 * g_xxy_yz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_x_zz_0_0[i] = -2.0 * g_y_zz_0_0[i] * a_exp + 4.0 * g_xxy_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_xx_0_0, g_x_xy_0_0, g_x_xz_0_0, g_x_yy_0_0, g_x_yz_0_0, g_x_zz_0_0, g_xy_0_0_0_y_xx_0_0, g_xy_0_0_0_y_xy_0_0, g_xy_0_0_0_y_xz_0_0, g_xy_0_0_0_y_yy_0_0, g_xy_0_0_0_y_yz_0_0, g_xy_0_0_0_y_zz_0_0, g_xyy_xx_0_0, g_xyy_xy_0_0, g_xyy_xz_0_0, g_xyy_yy_0_0, g_xyy_yz_0_0, g_xyy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_xx_0_0[i] = -2.0 * g_x_xx_0_0[i] * a_exp + 4.0 * g_xyy_xx_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xy_0_0[i] = -2.0 * g_x_xy_0_0[i] * a_exp + 4.0 * g_xyy_xy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_y_xz_0_0[i] = -2.0 * g_x_xz_0_0[i] * a_exp + 4.0 * g_xyy_xz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yy_0_0[i] = -2.0 * g_x_yy_0_0[i] * a_exp + 4.0 * g_xyy_yy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_y_yz_0_0[i] = -2.0 * g_x_yz_0_0[i] * a_exp + 4.0 * g_xyy_yz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_y_zz_0_0[i] = -2.0 * g_x_zz_0_0[i] * a_exp + 4.0 * g_xyy_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_xy_0_0_0_z_xx_0_0, g_xy_0_0_0_z_xy_0_0, g_xy_0_0_0_z_xz_0_0, g_xy_0_0_0_z_yy_0_0, g_xy_0_0_0_z_yz_0_0, g_xy_0_0_0_z_zz_0_0, g_xyz_xx_0_0, g_xyz_xy_0_0, g_xyz_xz_0_0, g_xyz_yy_0_0, g_xyz_yz_0_0, g_xyz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_xx_0_0[i] = 4.0 * g_xyz_xx_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xy_0_0[i] = 4.0 * g_xyz_xy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_z_xz_0_0[i] = 4.0 * g_xyz_xz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yy_0_0[i] = 4.0 * g_xyz_yy_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_z_yz_0_0[i] = 4.0 * g_xyz_yz_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_z_zz_0_0[i] = 4.0 * g_xyz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_xxz_xx_0_0, g_xxz_xy_0_0, g_xxz_xz_0_0, g_xxz_yy_0_0, g_xxz_yz_0_0, g_xxz_zz_0_0, g_xz_0_0_0_x_xx_0_0, g_xz_0_0_0_x_xy_0_0, g_xz_0_0_0_x_xz_0_0, g_xz_0_0_0_x_yy_0_0, g_xz_0_0_0_x_yz_0_0, g_xz_0_0_0_x_zz_0_0, g_z_xx_0_0, g_z_xy_0_0, g_z_xz_0_0, g_z_yy_0_0, g_z_yz_0_0, g_z_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_xx_0_0[i] = -2.0 * g_z_xx_0_0[i] * a_exp + 4.0 * g_xxz_xx_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xy_0_0[i] = -2.0 * g_z_xy_0_0[i] * a_exp + 4.0 * g_xxz_xy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_x_xz_0_0[i] = -2.0 * g_z_xz_0_0[i] * a_exp + 4.0 * g_xxz_xz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yy_0_0[i] = -2.0 * g_z_yy_0_0[i] * a_exp + 4.0 * g_xxz_yy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_x_yz_0_0[i] = -2.0 * g_z_yz_0_0[i] * a_exp + 4.0 * g_xxz_yz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_x_zz_0_0[i] = -2.0 * g_z_zz_0_0[i] * a_exp + 4.0 * g_xxz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_xyz_xx_0_0, g_xyz_xy_0_0, g_xyz_xz_0_0, g_xyz_yy_0_0, g_xyz_yz_0_0, g_xyz_zz_0_0, g_xz_0_0_0_y_xx_0_0, g_xz_0_0_0_y_xy_0_0, g_xz_0_0_0_y_xz_0_0, g_xz_0_0_0_y_yy_0_0, g_xz_0_0_0_y_yz_0_0, g_xz_0_0_0_y_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_xx_0_0[i] = 4.0 * g_xyz_xx_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xy_0_0[i] = 4.0 * g_xyz_xy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_y_xz_0_0[i] = 4.0 * g_xyz_xz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yy_0_0[i] = 4.0 * g_xyz_yy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_y_yz_0_0[i] = 4.0 * g_xyz_yz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_y_zz_0_0[i] = 4.0 * g_xyz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_x_xx_0_0, g_x_xy_0_0, g_x_xz_0_0, g_x_yy_0_0, g_x_yz_0_0, g_x_zz_0_0, g_xz_0_0_0_z_xx_0_0, g_xz_0_0_0_z_xy_0_0, g_xz_0_0_0_z_xz_0_0, g_xz_0_0_0_z_yy_0_0, g_xz_0_0_0_z_yz_0_0, g_xz_0_0_0_z_zz_0_0, g_xzz_xx_0_0, g_xzz_xy_0_0, g_xzz_xz_0_0, g_xzz_yy_0_0, g_xzz_yz_0_0, g_xzz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_xx_0_0[i] = -2.0 * g_x_xx_0_0[i] * a_exp + 4.0 * g_xzz_xx_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xy_0_0[i] = -2.0 * g_x_xy_0_0[i] * a_exp + 4.0 * g_xzz_xy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_z_xz_0_0[i] = -2.0 * g_x_xz_0_0[i] * a_exp + 4.0 * g_xzz_xz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yy_0_0[i] = -2.0 * g_x_yy_0_0[i] * a_exp + 4.0 * g_xzz_yy_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_z_yz_0_0[i] = -2.0 * g_x_yz_0_0[i] * a_exp + 4.0 * g_xzz_yz_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_z_zz_0_0[i] = -2.0 * g_x_zz_0_0[i] * a_exp + 4.0 * g_xzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_x_xx_0_0, g_x_xy_0_0, g_x_xz_0_0, g_x_yy_0_0, g_x_yz_0_0, g_x_zz_0_0, g_xyy_xx_0_0, g_xyy_xy_0_0, g_xyy_xz_0_0, g_xyy_yy_0_0, g_xyy_yz_0_0, g_xyy_zz_0_0, g_yy_0_0_0_x_xx_0_0, g_yy_0_0_0_x_xy_0_0, g_yy_0_0_0_x_xz_0_0, g_yy_0_0_0_x_yy_0_0, g_yy_0_0_0_x_yz_0_0, g_yy_0_0_0_x_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_xx_0_0[i] = -2.0 * g_x_xx_0_0[i] * a_exp + 4.0 * g_xyy_xx_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xy_0_0[i] = -2.0 * g_x_xy_0_0[i] * a_exp + 4.0 * g_xyy_xy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_x_xz_0_0[i] = -2.0 * g_x_xz_0_0[i] * a_exp + 4.0 * g_xyy_xz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yy_0_0[i] = -2.0 * g_x_yy_0_0[i] * a_exp + 4.0 * g_xyy_yy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_x_yz_0_0[i] = -2.0 * g_x_yz_0_0[i] * a_exp + 4.0 * g_xyy_yz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_x_zz_0_0[i] = -2.0 * g_x_zz_0_0[i] * a_exp + 4.0 * g_xyy_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_y_xx_0_0, g_y_xy_0_0, g_y_xz_0_0, g_y_yy_0_0, g_y_yz_0_0, g_y_zz_0_0, g_yy_0_0_0_y_xx_0_0, g_yy_0_0_0_y_xy_0_0, g_yy_0_0_0_y_xz_0_0, g_yy_0_0_0_y_yy_0_0, g_yy_0_0_0_y_yz_0_0, g_yy_0_0_0_y_zz_0_0, g_yyy_xx_0_0, g_yyy_xy_0_0, g_yyy_xz_0_0, g_yyy_yy_0_0, g_yyy_yz_0_0, g_yyy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_xx_0_0[i] = -6.0 * g_y_xx_0_0[i] * a_exp + 4.0 * g_yyy_xx_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xy_0_0[i] = -6.0 * g_y_xy_0_0[i] * a_exp + 4.0 * g_yyy_xy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_y_xz_0_0[i] = -6.0 * g_y_xz_0_0[i] * a_exp + 4.0 * g_yyy_xz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yy_0_0[i] = -6.0 * g_y_yy_0_0[i] * a_exp + 4.0 * g_yyy_yy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_y_yz_0_0[i] = -6.0 * g_y_yz_0_0[i] * a_exp + 4.0 * g_yyy_yz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_y_zz_0_0[i] = -6.0 * g_y_zz_0_0[i] * a_exp + 4.0 * g_yyy_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_yy_0_0_0_z_xx_0_0, g_yy_0_0_0_z_xy_0_0, g_yy_0_0_0_z_xz_0_0, g_yy_0_0_0_z_yy_0_0, g_yy_0_0_0_z_yz_0_0, g_yy_0_0_0_z_zz_0_0, g_yyz_xx_0_0, g_yyz_xy_0_0, g_yyz_xz_0_0, g_yyz_yy_0_0, g_yyz_yz_0_0, g_yyz_zz_0_0, g_z_xx_0_0, g_z_xy_0_0, g_z_xz_0_0, g_z_yy_0_0, g_z_yz_0_0, g_z_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_xx_0_0[i] = -2.0 * g_z_xx_0_0[i] * a_exp + 4.0 * g_yyz_xx_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xy_0_0[i] = -2.0 * g_z_xy_0_0[i] * a_exp + 4.0 * g_yyz_xy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_z_xz_0_0[i] = -2.0 * g_z_xz_0_0[i] * a_exp + 4.0 * g_yyz_xz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yy_0_0[i] = -2.0 * g_z_yy_0_0[i] * a_exp + 4.0 * g_yyz_yy_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_z_yz_0_0[i] = -2.0 * g_z_yz_0_0[i] * a_exp + 4.0 * g_yyz_yz_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_z_zz_0_0[i] = -2.0 * g_z_zz_0_0[i] * a_exp + 4.0 * g_yyz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_xyz_xx_0_0, g_xyz_xy_0_0, g_xyz_xz_0_0, g_xyz_yy_0_0, g_xyz_yz_0_0, g_xyz_zz_0_0, g_yz_0_0_0_x_xx_0_0, g_yz_0_0_0_x_xy_0_0, g_yz_0_0_0_x_xz_0_0, g_yz_0_0_0_x_yy_0_0, g_yz_0_0_0_x_yz_0_0, g_yz_0_0_0_x_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_xx_0_0[i] = 4.0 * g_xyz_xx_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xy_0_0[i] = 4.0 * g_xyz_xy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_x_xz_0_0[i] = 4.0 * g_xyz_xz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yy_0_0[i] = 4.0 * g_xyz_yy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_x_yz_0_0[i] = 4.0 * g_xyz_yz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_x_zz_0_0[i] = 4.0 * g_xyz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_yyz_xx_0_0, g_yyz_xy_0_0, g_yyz_xz_0_0, g_yyz_yy_0_0, g_yyz_yz_0_0, g_yyz_zz_0_0, g_yz_0_0_0_y_xx_0_0, g_yz_0_0_0_y_xy_0_0, g_yz_0_0_0_y_xz_0_0, g_yz_0_0_0_y_yy_0_0, g_yz_0_0_0_y_yz_0_0, g_yz_0_0_0_y_zz_0_0, g_z_xx_0_0, g_z_xy_0_0, g_z_xz_0_0, g_z_yy_0_0, g_z_yz_0_0, g_z_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_xx_0_0[i] = -2.0 * g_z_xx_0_0[i] * a_exp + 4.0 * g_yyz_xx_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xy_0_0[i] = -2.0 * g_z_xy_0_0[i] * a_exp + 4.0 * g_yyz_xy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_y_xz_0_0[i] = -2.0 * g_z_xz_0_0[i] * a_exp + 4.0 * g_yyz_xz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yy_0_0[i] = -2.0 * g_z_yy_0_0[i] * a_exp + 4.0 * g_yyz_yy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_y_yz_0_0[i] = -2.0 * g_z_yz_0_0[i] * a_exp + 4.0 * g_yyz_yz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_y_zz_0_0[i] = -2.0 * g_z_zz_0_0[i] * a_exp + 4.0 * g_yyz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_y_xx_0_0, g_y_xy_0_0, g_y_xz_0_0, g_y_yy_0_0, g_y_yz_0_0, g_y_zz_0_0, g_yz_0_0_0_z_xx_0_0, g_yz_0_0_0_z_xy_0_0, g_yz_0_0_0_z_xz_0_0, g_yz_0_0_0_z_yy_0_0, g_yz_0_0_0_z_yz_0_0, g_yz_0_0_0_z_zz_0_0, g_yzz_xx_0_0, g_yzz_xy_0_0, g_yzz_xz_0_0, g_yzz_yy_0_0, g_yzz_yz_0_0, g_yzz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_xx_0_0[i] = -2.0 * g_y_xx_0_0[i] * a_exp + 4.0 * g_yzz_xx_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xy_0_0[i] = -2.0 * g_y_xy_0_0[i] * a_exp + 4.0 * g_yzz_xy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_z_xz_0_0[i] = -2.0 * g_y_xz_0_0[i] * a_exp + 4.0 * g_yzz_xz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yy_0_0[i] = -2.0 * g_y_yy_0_0[i] * a_exp + 4.0 * g_yzz_yy_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_z_yz_0_0[i] = -2.0 * g_y_yz_0_0[i] * a_exp + 4.0 * g_yzz_yz_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_z_zz_0_0[i] = -2.0 * g_y_zz_0_0[i] * a_exp + 4.0 * g_yzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_x_xx_0_0, g_x_xy_0_0, g_x_xz_0_0, g_x_yy_0_0, g_x_yz_0_0, g_x_zz_0_0, g_xzz_xx_0_0, g_xzz_xy_0_0, g_xzz_xz_0_0, g_xzz_yy_0_0, g_xzz_yz_0_0, g_xzz_zz_0_0, g_zz_0_0_0_x_xx_0_0, g_zz_0_0_0_x_xy_0_0, g_zz_0_0_0_x_xz_0_0, g_zz_0_0_0_x_yy_0_0, g_zz_0_0_0_x_yz_0_0, g_zz_0_0_0_x_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_xx_0_0[i] = -2.0 * g_x_xx_0_0[i] * a_exp + 4.0 * g_xzz_xx_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xy_0_0[i] = -2.0 * g_x_xy_0_0[i] * a_exp + 4.0 * g_xzz_xy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_x_xz_0_0[i] = -2.0 * g_x_xz_0_0[i] * a_exp + 4.0 * g_xzz_xz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yy_0_0[i] = -2.0 * g_x_yy_0_0[i] * a_exp + 4.0 * g_xzz_yy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_x_yz_0_0[i] = -2.0 * g_x_yz_0_0[i] * a_exp + 4.0 * g_xzz_yz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_x_zz_0_0[i] = -2.0 * g_x_zz_0_0[i] * a_exp + 4.0 * g_xzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_y_xx_0_0, g_y_xy_0_0, g_y_xz_0_0, g_y_yy_0_0, g_y_yz_0_0, g_y_zz_0_0, g_yzz_xx_0_0, g_yzz_xy_0_0, g_yzz_xz_0_0, g_yzz_yy_0_0, g_yzz_yz_0_0, g_yzz_zz_0_0, g_zz_0_0_0_y_xx_0_0, g_zz_0_0_0_y_xy_0_0, g_zz_0_0_0_y_xz_0_0, g_zz_0_0_0_y_yy_0_0, g_zz_0_0_0_y_yz_0_0, g_zz_0_0_0_y_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_xx_0_0[i] = -2.0 * g_y_xx_0_0[i] * a_exp + 4.0 * g_yzz_xx_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xy_0_0[i] = -2.0 * g_y_xy_0_0[i] * a_exp + 4.0 * g_yzz_xy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_y_xz_0_0[i] = -2.0 * g_y_xz_0_0[i] * a_exp + 4.0 * g_yzz_xz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yy_0_0[i] = -2.0 * g_y_yy_0_0[i] * a_exp + 4.0 * g_yzz_yy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_y_yz_0_0[i] = -2.0 * g_y_yz_0_0[i] * a_exp + 4.0 * g_yzz_yz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_y_zz_0_0[i] = -2.0 * g_y_zz_0_0[i] * a_exp + 4.0 * g_yzz_zz_0_0[i] * a_exp * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_z_xx_0_0, g_z_xy_0_0, g_z_xz_0_0, g_z_yy_0_0, g_z_yz_0_0, g_z_zz_0_0, g_zz_0_0_0_z_xx_0_0, g_zz_0_0_0_z_xy_0_0, g_zz_0_0_0_z_xz_0_0, g_zz_0_0_0_z_yy_0_0, g_zz_0_0_0_z_yz_0_0, g_zz_0_0_0_z_zz_0_0, g_zzz_xx_0_0, g_zzz_xy_0_0, g_zzz_xz_0_0, g_zzz_yy_0_0, g_zzz_yz_0_0, g_zzz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_xx_0_0[i] = -6.0 * g_z_xx_0_0[i] * a_exp + 4.0 * g_zzz_xx_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xy_0_0[i] = -6.0 * g_z_xy_0_0[i] * a_exp + 4.0 * g_zzz_xy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_z_xz_0_0[i] = -6.0 * g_z_xz_0_0[i] * a_exp + 4.0 * g_zzz_xz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yy_0_0[i] = -6.0 * g_z_yy_0_0[i] * a_exp + 4.0 * g_zzz_yy_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_z_yz_0_0[i] = -6.0 * g_z_yz_0_0[i] * a_exp + 4.0 * g_zzz_yz_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_z_zz_0_0[i] = -6.0 * g_z_zz_0_0[i] * a_exp + 4.0 * g_zzz_zz_0_0[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

