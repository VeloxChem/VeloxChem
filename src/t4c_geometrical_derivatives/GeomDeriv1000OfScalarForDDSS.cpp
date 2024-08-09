#include "GeomDeriv1000OfScalarForDDSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_ddss_0(CSimdArray<double>& buffer_1000_ddss,
                     const CSimdArray<double>& buffer_pdss,
                     const CSimdArray<double>& buffer_fdss,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_ddss.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_ddss

    auto g_x_0_0_0_xx_xx_0_0 = buffer_1000_ddss[0];

    auto g_x_0_0_0_xx_xy_0_0 = buffer_1000_ddss[1];

    auto g_x_0_0_0_xx_xz_0_0 = buffer_1000_ddss[2];

    auto g_x_0_0_0_xx_yy_0_0 = buffer_1000_ddss[3];

    auto g_x_0_0_0_xx_yz_0_0 = buffer_1000_ddss[4];

    auto g_x_0_0_0_xx_zz_0_0 = buffer_1000_ddss[5];

    auto g_x_0_0_0_xy_xx_0_0 = buffer_1000_ddss[6];

    auto g_x_0_0_0_xy_xy_0_0 = buffer_1000_ddss[7];

    auto g_x_0_0_0_xy_xz_0_0 = buffer_1000_ddss[8];

    auto g_x_0_0_0_xy_yy_0_0 = buffer_1000_ddss[9];

    auto g_x_0_0_0_xy_yz_0_0 = buffer_1000_ddss[10];

    auto g_x_0_0_0_xy_zz_0_0 = buffer_1000_ddss[11];

    auto g_x_0_0_0_xz_xx_0_0 = buffer_1000_ddss[12];

    auto g_x_0_0_0_xz_xy_0_0 = buffer_1000_ddss[13];

    auto g_x_0_0_0_xz_xz_0_0 = buffer_1000_ddss[14];

    auto g_x_0_0_0_xz_yy_0_0 = buffer_1000_ddss[15];

    auto g_x_0_0_0_xz_yz_0_0 = buffer_1000_ddss[16];

    auto g_x_0_0_0_xz_zz_0_0 = buffer_1000_ddss[17];

    auto g_x_0_0_0_yy_xx_0_0 = buffer_1000_ddss[18];

    auto g_x_0_0_0_yy_xy_0_0 = buffer_1000_ddss[19];

    auto g_x_0_0_0_yy_xz_0_0 = buffer_1000_ddss[20];

    auto g_x_0_0_0_yy_yy_0_0 = buffer_1000_ddss[21];

    auto g_x_0_0_0_yy_yz_0_0 = buffer_1000_ddss[22];

    auto g_x_0_0_0_yy_zz_0_0 = buffer_1000_ddss[23];

    auto g_x_0_0_0_yz_xx_0_0 = buffer_1000_ddss[24];

    auto g_x_0_0_0_yz_xy_0_0 = buffer_1000_ddss[25];

    auto g_x_0_0_0_yz_xz_0_0 = buffer_1000_ddss[26];

    auto g_x_0_0_0_yz_yy_0_0 = buffer_1000_ddss[27];

    auto g_x_0_0_0_yz_yz_0_0 = buffer_1000_ddss[28];

    auto g_x_0_0_0_yz_zz_0_0 = buffer_1000_ddss[29];

    auto g_x_0_0_0_zz_xx_0_0 = buffer_1000_ddss[30];

    auto g_x_0_0_0_zz_xy_0_0 = buffer_1000_ddss[31];

    auto g_x_0_0_0_zz_xz_0_0 = buffer_1000_ddss[32];

    auto g_x_0_0_0_zz_yy_0_0 = buffer_1000_ddss[33];

    auto g_x_0_0_0_zz_yz_0_0 = buffer_1000_ddss[34];

    auto g_x_0_0_0_zz_zz_0_0 = buffer_1000_ddss[35];

    auto g_y_0_0_0_xx_xx_0_0 = buffer_1000_ddss[36];

    auto g_y_0_0_0_xx_xy_0_0 = buffer_1000_ddss[37];

    auto g_y_0_0_0_xx_xz_0_0 = buffer_1000_ddss[38];

    auto g_y_0_0_0_xx_yy_0_0 = buffer_1000_ddss[39];

    auto g_y_0_0_0_xx_yz_0_0 = buffer_1000_ddss[40];

    auto g_y_0_0_0_xx_zz_0_0 = buffer_1000_ddss[41];

    auto g_y_0_0_0_xy_xx_0_0 = buffer_1000_ddss[42];

    auto g_y_0_0_0_xy_xy_0_0 = buffer_1000_ddss[43];

    auto g_y_0_0_0_xy_xz_0_0 = buffer_1000_ddss[44];

    auto g_y_0_0_0_xy_yy_0_0 = buffer_1000_ddss[45];

    auto g_y_0_0_0_xy_yz_0_0 = buffer_1000_ddss[46];

    auto g_y_0_0_0_xy_zz_0_0 = buffer_1000_ddss[47];

    auto g_y_0_0_0_xz_xx_0_0 = buffer_1000_ddss[48];

    auto g_y_0_0_0_xz_xy_0_0 = buffer_1000_ddss[49];

    auto g_y_0_0_0_xz_xz_0_0 = buffer_1000_ddss[50];

    auto g_y_0_0_0_xz_yy_0_0 = buffer_1000_ddss[51];

    auto g_y_0_0_0_xz_yz_0_0 = buffer_1000_ddss[52];

    auto g_y_0_0_0_xz_zz_0_0 = buffer_1000_ddss[53];

    auto g_y_0_0_0_yy_xx_0_0 = buffer_1000_ddss[54];

    auto g_y_0_0_0_yy_xy_0_0 = buffer_1000_ddss[55];

    auto g_y_0_0_0_yy_xz_0_0 = buffer_1000_ddss[56];

    auto g_y_0_0_0_yy_yy_0_0 = buffer_1000_ddss[57];

    auto g_y_0_0_0_yy_yz_0_0 = buffer_1000_ddss[58];

    auto g_y_0_0_0_yy_zz_0_0 = buffer_1000_ddss[59];

    auto g_y_0_0_0_yz_xx_0_0 = buffer_1000_ddss[60];

    auto g_y_0_0_0_yz_xy_0_0 = buffer_1000_ddss[61];

    auto g_y_0_0_0_yz_xz_0_0 = buffer_1000_ddss[62];

    auto g_y_0_0_0_yz_yy_0_0 = buffer_1000_ddss[63];

    auto g_y_0_0_0_yz_yz_0_0 = buffer_1000_ddss[64];

    auto g_y_0_0_0_yz_zz_0_0 = buffer_1000_ddss[65];

    auto g_y_0_0_0_zz_xx_0_0 = buffer_1000_ddss[66];

    auto g_y_0_0_0_zz_xy_0_0 = buffer_1000_ddss[67];

    auto g_y_0_0_0_zz_xz_0_0 = buffer_1000_ddss[68];

    auto g_y_0_0_0_zz_yy_0_0 = buffer_1000_ddss[69];

    auto g_y_0_0_0_zz_yz_0_0 = buffer_1000_ddss[70];

    auto g_y_0_0_0_zz_zz_0_0 = buffer_1000_ddss[71];

    auto g_z_0_0_0_xx_xx_0_0 = buffer_1000_ddss[72];

    auto g_z_0_0_0_xx_xy_0_0 = buffer_1000_ddss[73];

    auto g_z_0_0_0_xx_xz_0_0 = buffer_1000_ddss[74];

    auto g_z_0_0_0_xx_yy_0_0 = buffer_1000_ddss[75];

    auto g_z_0_0_0_xx_yz_0_0 = buffer_1000_ddss[76];

    auto g_z_0_0_0_xx_zz_0_0 = buffer_1000_ddss[77];

    auto g_z_0_0_0_xy_xx_0_0 = buffer_1000_ddss[78];

    auto g_z_0_0_0_xy_xy_0_0 = buffer_1000_ddss[79];

    auto g_z_0_0_0_xy_xz_0_0 = buffer_1000_ddss[80];

    auto g_z_0_0_0_xy_yy_0_0 = buffer_1000_ddss[81];

    auto g_z_0_0_0_xy_yz_0_0 = buffer_1000_ddss[82];

    auto g_z_0_0_0_xy_zz_0_0 = buffer_1000_ddss[83];

    auto g_z_0_0_0_xz_xx_0_0 = buffer_1000_ddss[84];

    auto g_z_0_0_0_xz_xy_0_0 = buffer_1000_ddss[85];

    auto g_z_0_0_0_xz_xz_0_0 = buffer_1000_ddss[86];

    auto g_z_0_0_0_xz_yy_0_0 = buffer_1000_ddss[87];

    auto g_z_0_0_0_xz_yz_0_0 = buffer_1000_ddss[88];

    auto g_z_0_0_0_xz_zz_0_0 = buffer_1000_ddss[89];

    auto g_z_0_0_0_yy_xx_0_0 = buffer_1000_ddss[90];

    auto g_z_0_0_0_yy_xy_0_0 = buffer_1000_ddss[91];

    auto g_z_0_0_0_yy_xz_0_0 = buffer_1000_ddss[92];

    auto g_z_0_0_0_yy_yy_0_0 = buffer_1000_ddss[93];

    auto g_z_0_0_0_yy_yz_0_0 = buffer_1000_ddss[94];

    auto g_z_0_0_0_yy_zz_0_0 = buffer_1000_ddss[95];

    auto g_z_0_0_0_yz_xx_0_0 = buffer_1000_ddss[96];

    auto g_z_0_0_0_yz_xy_0_0 = buffer_1000_ddss[97];

    auto g_z_0_0_0_yz_xz_0_0 = buffer_1000_ddss[98];

    auto g_z_0_0_0_yz_yy_0_0 = buffer_1000_ddss[99];

    auto g_z_0_0_0_yz_yz_0_0 = buffer_1000_ddss[100];

    auto g_z_0_0_0_yz_zz_0_0 = buffer_1000_ddss[101];

    auto g_z_0_0_0_zz_xx_0_0 = buffer_1000_ddss[102];

    auto g_z_0_0_0_zz_xy_0_0 = buffer_1000_ddss[103];

    auto g_z_0_0_0_zz_xz_0_0 = buffer_1000_ddss[104];

    auto g_z_0_0_0_zz_yy_0_0 = buffer_1000_ddss[105];

    auto g_z_0_0_0_zz_yz_0_0 = buffer_1000_ddss[106];

    auto g_z_0_0_0_zz_zz_0_0 = buffer_1000_ddss[107];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_0_0_0_xx_xx_0_0, g_x_0_0_0_xx_xy_0_0, g_x_0_0_0_xx_xz_0_0, g_x_0_0_0_xx_yy_0_0, g_x_0_0_0_xx_yz_0_0, g_x_0_0_0_xx_zz_0_0, g_x_xx_0_0, g_x_xy_0_0, g_x_xz_0_0, g_x_yy_0_0, g_x_yz_0_0, g_x_zz_0_0, g_xxx_xx_0_0, g_xxx_xy_0_0, g_xxx_xz_0_0, g_xxx_yy_0_0, g_xxx_yz_0_0, g_xxx_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xx_xx_0_0[i] = -2.0 * g_x_xx_0_0[i] + 2.0 * g_xxx_xx_0_0[i] * a_exp;

        g_x_0_0_0_xx_xy_0_0[i] = -2.0 * g_x_xy_0_0[i] + 2.0 * g_xxx_xy_0_0[i] * a_exp;

        g_x_0_0_0_xx_xz_0_0[i] = -2.0 * g_x_xz_0_0[i] + 2.0 * g_xxx_xz_0_0[i] * a_exp;

        g_x_0_0_0_xx_yy_0_0[i] = -2.0 * g_x_yy_0_0[i] + 2.0 * g_xxx_yy_0_0[i] * a_exp;

        g_x_0_0_0_xx_yz_0_0[i] = -2.0 * g_x_yz_0_0[i] + 2.0 * g_xxx_yz_0_0[i] * a_exp;

        g_x_0_0_0_xx_zz_0_0[i] = -2.0 * g_x_zz_0_0[i] + 2.0 * g_xxx_zz_0_0[i] * a_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_0_0_0_xy_xx_0_0, g_x_0_0_0_xy_xy_0_0, g_x_0_0_0_xy_xz_0_0, g_x_0_0_0_xy_yy_0_0, g_x_0_0_0_xy_yz_0_0, g_x_0_0_0_xy_zz_0_0, g_xxy_xx_0_0, g_xxy_xy_0_0, g_xxy_xz_0_0, g_xxy_yy_0_0, g_xxy_yz_0_0, g_xxy_zz_0_0, g_y_xx_0_0, g_y_xy_0_0, g_y_xz_0_0, g_y_yy_0_0, g_y_yz_0_0, g_y_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xy_xx_0_0[i] = -g_y_xx_0_0[i] + 2.0 * g_xxy_xx_0_0[i] * a_exp;

        g_x_0_0_0_xy_xy_0_0[i] = -g_y_xy_0_0[i] + 2.0 * g_xxy_xy_0_0[i] * a_exp;

        g_x_0_0_0_xy_xz_0_0[i] = -g_y_xz_0_0[i] + 2.0 * g_xxy_xz_0_0[i] * a_exp;

        g_x_0_0_0_xy_yy_0_0[i] = -g_y_yy_0_0[i] + 2.0 * g_xxy_yy_0_0[i] * a_exp;

        g_x_0_0_0_xy_yz_0_0[i] = -g_y_yz_0_0[i] + 2.0 * g_xxy_yz_0_0[i] * a_exp;

        g_x_0_0_0_xy_zz_0_0[i] = -g_y_zz_0_0[i] + 2.0 * g_xxy_zz_0_0[i] * a_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_0_0_0_xz_xx_0_0, g_x_0_0_0_xz_xy_0_0, g_x_0_0_0_xz_xz_0_0, g_x_0_0_0_xz_yy_0_0, g_x_0_0_0_xz_yz_0_0, g_x_0_0_0_xz_zz_0_0, g_xxz_xx_0_0, g_xxz_xy_0_0, g_xxz_xz_0_0, g_xxz_yy_0_0, g_xxz_yz_0_0, g_xxz_zz_0_0, g_z_xx_0_0, g_z_xy_0_0, g_z_xz_0_0, g_z_yy_0_0, g_z_yz_0_0, g_z_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_xz_xx_0_0[i] = -g_z_xx_0_0[i] + 2.0 * g_xxz_xx_0_0[i] * a_exp;

        g_x_0_0_0_xz_xy_0_0[i] = -g_z_xy_0_0[i] + 2.0 * g_xxz_xy_0_0[i] * a_exp;

        g_x_0_0_0_xz_xz_0_0[i] = -g_z_xz_0_0[i] + 2.0 * g_xxz_xz_0_0[i] * a_exp;

        g_x_0_0_0_xz_yy_0_0[i] = -g_z_yy_0_0[i] + 2.0 * g_xxz_yy_0_0[i] * a_exp;

        g_x_0_0_0_xz_yz_0_0[i] = -g_z_yz_0_0[i] + 2.0 * g_xxz_yz_0_0[i] * a_exp;

        g_x_0_0_0_xz_zz_0_0[i] = -g_z_zz_0_0[i] + 2.0 * g_xxz_zz_0_0[i] * a_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_x_0_0_0_yy_xx_0_0, g_x_0_0_0_yy_xy_0_0, g_x_0_0_0_yy_xz_0_0, g_x_0_0_0_yy_yy_0_0, g_x_0_0_0_yy_yz_0_0, g_x_0_0_0_yy_zz_0_0, g_xyy_xx_0_0, g_xyy_xy_0_0, g_xyy_xz_0_0, g_xyy_yy_0_0, g_xyy_yz_0_0, g_xyy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yy_xx_0_0[i] = 2.0 * g_xyy_xx_0_0[i] * a_exp;

        g_x_0_0_0_yy_xy_0_0[i] = 2.0 * g_xyy_xy_0_0[i] * a_exp;

        g_x_0_0_0_yy_xz_0_0[i] = 2.0 * g_xyy_xz_0_0[i] * a_exp;

        g_x_0_0_0_yy_yy_0_0[i] = 2.0 * g_xyy_yy_0_0[i] * a_exp;

        g_x_0_0_0_yy_yz_0_0[i] = 2.0 * g_xyy_yz_0_0[i] * a_exp;

        g_x_0_0_0_yy_zz_0_0[i] = 2.0 * g_xyy_zz_0_0[i] * a_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_x_0_0_0_yz_xx_0_0, g_x_0_0_0_yz_xy_0_0, g_x_0_0_0_yz_xz_0_0, g_x_0_0_0_yz_yy_0_0, g_x_0_0_0_yz_yz_0_0, g_x_0_0_0_yz_zz_0_0, g_xyz_xx_0_0, g_xyz_xy_0_0, g_xyz_xz_0_0, g_xyz_yy_0_0, g_xyz_yz_0_0, g_xyz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_yz_xx_0_0[i] = 2.0 * g_xyz_xx_0_0[i] * a_exp;

        g_x_0_0_0_yz_xy_0_0[i] = 2.0 * g_xyz_xy_0_0[i] * a_exp;

        g_x_0_0_0_yz_xz_0_0[i] = 2.0 * g_xyz_xz_0_0[i] * a_exp;

        g_x_0_0_0_yz_yy_0_0[i] = 2.0 * g_xyz_yy_0_0[i] * a_exp;

        g_x_0_0_0_yz_yz_0_0[i] = 2.0 * g_xyz_yz_0_0[i] * a_exp;

        g_x_0_0_0_yz_zz_0_0[i] = 2.0 * g_xyz_zz_0_0[i] * a_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_x_0_0_0_zz_xx_0_0, g_x_0_0_0_zz_xy_0_0, g_x_0_0_0_zz_xz_0_0, g_x_0_0_0_zz_yy_0_0, g_x_0_0_0_zz_yz_0_0, g_x_0_0_0_zz_zz_0_0, g_xzz_xx_0_0, g_xzz_xy_0_0, g_xzz_xz_0_0, g_xzz_yy_0_0, g_xzz_yz_0_0, g_xzz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_zz_xx_0_0[i] = 2.0 * g_xzz_xx_0_0[i] * a_exp;

        g_x_0_0_0_zz_xy_0_0[i] = 2.0 * g_xzz_xy_0_0[i] * a_exp;

        g_x_0_0_0_zz_xz_0_0[i] = 2.0 * g_xzz_xz_0_0[i] * a_exp;

        g_x_0_0_0_zz_yy_0_0[i] = 2.0 * g_xzz_yy_0_0[i] * a_exp;

        g_x_0_0_0_zz_yz_0_0[i] = 2.0 * g_xzz_yz_0_0[i] * a_exp;

        g_x_0_0_0_zz_zz_0_0[i] = 2.0 * g_xzz_zz_0_0[i] * a_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_xxy_xx_0_0, g_xxy_xy_0_0, g_xxy_xz_0_0, g_xxy_yy_0_0, g_xxy_yz_0_0, g_xxy_zz_0_0, g_y_0_0_0_xx_xx_0_0, g_y_0_0_0_xx_xy_0_0, g_y_0_0_0_xx_xz_0_0, g_y_0_0_0_xx_yy_0_0, g_y_0_0_0_xx_yz_0_0, g_y_0_0_0_xx_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xx_xx_0_0[i] = 2.0 * g_xxy_xx_0_0[i] * a_exp;

        g_y_0_0_0_xx_xy_0_0[i] = 2.0 * g_xxy_xy_0_0[i] * a_exp;

        g_y_0_0_0_xx_xz_0_0[i] = 2.0 * g_xxy_xz_0_0[i] * a_exp;

        g_y_0_0_0_xx_yy_0_0[i] = 2.0 * g_xxy_yy_0_0[i] * a_exp;

        g_y_0_0_0_xx_yz_0_0[i] = 2.0 * g_xxy_yz_0_0[i] * a_exp;

        g_y_0_0_0_xx_zz_0_0[i] = 2.0 * g_xxy_zz_0_0[i] * a_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_x_xx_0_0, g_x_xy_0_0, g_x_xz_0_0, g_x_yy_0_0, g_x_yz_0_0, g_x_zz_0_0, g_xyy_xx_0_0, g_xyy_xy_0_0, g_xyy_xz_0_0, g_xyy_yy_0_0, g_xyy_yz_0_0, g_xyy_zz_0_0, g_y_0_0_0_xy_xx_0_0, g_y_0_0_0_xy_xy_0_0, g_y_0_0_0_xy_xz_0_0, g_y_0_0_0_xy_yy_0_0, g_y_0_0_0_xy_yz_0_0, g_y_0_0_0_xy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xy_xx_0_0[i] = -g_x_xx_0_0[i] + 2.0 * g_xyy_xx_0_0[i] * a_exp;

        g_y_0_0_0_xy_xy_0_0[i] = -g_x_xy_0_0[i] + 2.0 * g_xyy_xy_0_0[i] * a_exp;

        g_y_0_0_0_xy_xz_0_0[i] = -g_x_xz_0_0[i] + 2.0 * g_xyy_xz_0_0[i] * a_exp;

        g_y_0_0_0_xy_yy_0_0[i] = -g_x_yy_0_0[i] + 2.0 * g_xyy_yy_0_0[i] * a_exp;

        g_y_0_0_0_xy_yz_0_0[i] = -g_x_yz_0_0[i] + 2.0 * g_xyy_yz_0_0[i] * a_exp;

        g_y_0_0_0_xy_zz_0_0[i] = -g_x_zz_0_0[i] + 2.0 * g_xyy_zz_0_0[i] * a_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_xyz_xx_0_0, g_xyz_xy_0_0, g_xyz_xz_0_0, g_xyz_yy_0_0, g_xyz_yz_0_0, g_xyz_zz_0_0, g_y_0_0_0_xz_xx_0_0, g_y_0_0_0_xz_xy_0_0, g_y_0_0_0_xz_xz_0_0, g_y_0_0_0_xz_yy_0_0, g_y_0_0_0_xz_yz_0_0, g_y_0_0_0_xz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_xz_xx_0_0[i] = 2.0 * g_xyz_xx_0_0[i] * a_exp;

        g_y_0_0_0_xz_xy_0_0[i] = 2.0 * g_xyz_xy_0_0[i] * a_exp;

        g_y_0_0_0_xz_xz_0_0[i] = 2.0 * g_xyz_xz_0_0[i] * a_exp;

        g_y_0_0_0_xz_yy_0_0[i] = 2.0 * g_xyz_yy_0_0[i] * a_exp;

        g_y_0_0_0_xz_yz_0_0[i] = 2.0 * g_xyz_yz_0_0[i] * a_exp;

        g_y_0_0_0_xz_zz_0_0[i] = 2.0 * g_xyz_zz_0_0[i] * a_exp;
    }
    // integrals block (54-60)

    #pragma omp simd aligned(g_y_0_0_0_yy_xx_0_0, g_y_0_0_0_yy_xy_0_0, g_y_0_0_0_yy_xz_0_0, g_y_0_0_0_yy_yy_0_0, g_y_0_0_0_yy_yz_0_0, g_y_0_0_0_yy_zz_0_0, g_y_xx_0_0, g_y_xy_0_0, g_y_xz_0_0, g_y_yy_0_0, g_y_yz_0_0, g_y_zz_0_0, g_yyy_xx_0_0, g_yyy_xy_0_0, g_yyy_xz_0_0, g_yyy_yy_0_0, g_yyy_yz_0_0, g_yyy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yy_xx_0_0[i] = -2.0 * g_y_xx_0_0[i] + 2.0 * g_yyy_xx_0_0[i] * a_exp;

        g_y_0_0_0_yy_xy_0_0[i] = -2.0 * g_y_xy_0_0[i] + 2.0 * g_yyy_xy_0_0[i] * a_exp;

        g_y_0_0_0_yy_xz_0_0[i] = -2.0 * g_y_xz_0_0[i] + 2.0 * g_yyy_xz_0_0[i] * a_exp;

        g_y_0_0_0_yy_yy_0_0[i] = -2.0 * g_y_yy_0_0[i] + 2.0 * g_yyy_yy_0_0[i] * a_exp;

        g_y_0_0_0_yy_yz_0_0[i] = -2.0 * g_y_yz_0_0[i] + 2.0 * g_yyy_yz_0_0[i] * a_exp;

        g_y_0_0_0_yy_zz_0_0[i] = -2.0 * g_y_zz_0_0[i] + 2.0 * g_yyy_zz_0_0[i] * a_exp;
    }
    // integrals block (60-66)

    #pragma omp simd aligned(g_y_0_0_0_yz_xx_0_0, g_y_0_0_0_yz_xy_0_0, g_y_0_0_0_yz_xz_0_0, g_y_0_0_0_yz_yy_0_0, g_y_0_0_0_yz_yz_0_0, g_y_0_0_0_yz_zz_0_0, g_yyz_xx_0_0, g_yyz_xy_0_0, g_yyz_xz_0_0, g_yyz_yy_0_0, g_yyz_yz_0_0, g_yyz_zz_0_0, g_z_xx_0_0, g_z_xy_0_0, g_z_xz_0_0, g_z_yy_0_0, g_z_yz_0_0, g_z_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_yz_xx_0_0[i] = -g_z_xx_0_0[i] + 2.0 * g_yyz_xx_0_0[i] * a_exp;

        g_y_0_0_0_yz_xy_0_0[i] = -g_z_xy_0_0[i] + 2.0 * g_yyz_xy_0_0[i] * a_exp;

        g_y_0_0_0_yz_xz_0_0[i] = -g_z_xz_0_0[i] + 2.0 * g_yyz_xz_0_0[i] * a_exp;

        g_y_0_0_0_yz_yy_0_0[i] = -g_z_yy_0_0[i] + 2.0 * g_yyz_yy_0_0[i] * a_exp;

        g_y_0_0_0_yz_yz_0_0[i] = -g_z_yz_0_0[i] + 2.0 * g_yyz_yz_0_0[i] * a_exp;

        g_y_0_0_0_yz_zz_0_0[i] = -g_z_zz_0_0[i] + 2.0 * g_yyz_zz_0_0[i] * a_exp;
    }
    // integrals block (66-72)

    #pragma omp simd aligned(g_y_0_0_0_zz_xx_0_0, g_y_0_0_0_zz_xy_0_0, g_y_0_0_0_zz_xz_0_0, g_y_0_0_0_zz_yy_0_0, g_y_0_0_0_zz_yz_0_0, g_y_0_0_0_zz_zz_0_0, g_yzz_xx_0_0, g_yzz_xy_0_0, g_yzz_xz_0_0, g_yzz_yy_0_0, g_yzz_yz_0_0, g_yzz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_zz_xx_0_0[i] = 2.0 * g_yzz_xx_0_0[i] * a_exp;

        g_y_0_0_0_zz_xy_0_0[i] = 2.0 * g_yzz_xy_0_0[i] * a_exp;

        g_y_0_0_0_zz_xz_0_0[i] = 2.0 * g_yzz_xz_0_0[i] * a_exp;

        g_y_0_0_0_zz_yy_0_0[i] = 2.0 * g_yzz_yy_0_0[i] * a_exp;

        g_y_0_0_0_zz_yz_0_0[i] = 2.0 * g_yzz_yz_0_0[i] * a_exp;

        g_y_0_0_0_zz_zz_0_0[i] = 2.0 * g_yzz_zz_0_0[i] * a_exp;
    }
    // integrals block (72-78)

    #pragma omp simd aligned(g_xxz_xx_0_0, g_xxz_xy_0_0, g_xxz_xz_0_0, g_xxz_yy_0_0, g_xxz_yz_0_0, g_xxz_zz_0_0, g_z_0_0_0_xx_xx_0_0, g_z_0_0_0_xx_xy_0_0, g_z_0_0_0_xx_xz_0_0, g_z_0_0_0_xx_yy_0_0, g_z_0_0_0_xx_yz_0_0, g_z_0_0_0_xx_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xx_xx_0_0[i] = 2.0 * g_xxz_xx_0_0[i] * a_exp;

        g_z_0_0_0_xx_xy_0_0[i] = 2.0 * g_xxz_xy_0_0[i] * a_exp;

        g_z_0_0_0_xx_xz_0_0[i] = 2.0 * g_xxz_xz_0_0[i] * a_exp;

        g_z_0_0_0_xx_yy_0_0[i] = 2.0 * g_xxz_yy_0_0[i] * a_exp;

        g_z_0_0_0_xx_yz_0_0[i] = 2.0 * g_xxz_yz_0_0[i] * a_exp;

        g_z_0_0_0_xx_zz_0_0[i] = 2.0 * g_xxz_zz_0_0[i] * a_exp;
    }
    // integrals block (78-84)

    #pragma omp simd aligned(g_xyz_xx_0_0, g_xyz_xy_0_0, g_xyz_xz_0_0, g_xyz_yy_0_0, g_xyz_yz_0_0, g_xyz_zz_0_0, g_z_0_0_0_xy_xx_0_0, g_z_0_0_0_xy_xy_0_0, g_z_0_0_0_xy_xz_0_0, g_z_0_0_0_xy_yy_0_0, g_z_0_0_0_xy_yz_0_0, g_z_0_0_0_xy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xy_xx_0_0[i] = 2.0 * g_xyz_xx_0_0[i] * a_exp;

        g_z_0_0_0_xy_xy_0_0[i] = 2.0 * g_xyz_xy_0_0[i] * a_exp;

        g_z_0_0_0_xy_xz_0_0[i] = 2.0 * g_xyz_xz_0_0[i] * a_exp;

        g_z_0_0_0_xy_yy_0_0[i] = 2.0 * g_xyz_yy_0_0[i] * a_exp;

        g_z_0_0_0_xy_yz_0_0[i] = 2.0 * g_xyz_yz_0_0[i] * a_exp;

        g_z_0_0_0_xy_zz_0_0[i] = 2.0 * g_xyz_zz_0_0[i] * a_exp;
    }
    // integrals block (84-90)

    #pragma omp simd aligned(g_x_xx_0_0, g_x_xy_0_0, g_x_xz_0_0, g_x_yy_0_0, g_x_yz_0_0, g_x_zz_0_0, g_xzz_xx_0_0, g_xzz_xy_0_0, g_xzz_xz_0_0, g_xzz_yy_0_0, g_xzz_yz_0_0, g_xzz_zz_0_0, g_z_0_0_0_xz_xx_0_0, g_z_0_0_0_xz_xy_0_0, g_z_0_0_0_xz_xz_0_0, g_z_0_0_0_xz_yy_0_0, g_z_0_0_0_xz_yz_0_0, g_z_0_0_0_xz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_xz_xx_0_0[i] = -g_x_xx_0_0[i] + 2.0 * g_xzz_xx_0_0[i] * a_exp;

        g_z_0_0_0_xz_xy_0_0[i] = -g_x_xy_0_0[i] + 2.0 * g_xzz_xy_0_0[i] * a_exp;

        g_z_0_0_0_xz_xz_0_0[i] = -g_x_xz_0_0[i] + 2.0 * g_xzz_xz_0_0[i] * a_exp;

        g_z_0_0_0_xz_yy_0_0[i] = -g_x_yy_0_0[i] + 2.0 * g_xzz_yy_0_0[i] * a_exp;

        g_z_0_0_0_xz_yz_0_0[i] = -g_x_yz_0_0[i] + 2.0 * g_xzz_yz_0_0[i] * a_exp;

        g_z_0_0_0_xz_zz_0_0[i] = -g_x_zz_0_0[i] + 2.0 * g_xzz_zz_0_0[i] * a_exp;
    }
    // integrals block (90-96)

    #pragma omp simd aligned(g_yyz_xx_0_0, g_yyz_xy_0_0, g_yyz_xz_0_0, g_yyz_yy_0_0, g_yyz_yz_0_0, g_yyz_zz_0_0, g_z_0_0_0_yy_xx_0_0, g_z_0_0_0_yy_xy_0_0, g_z_0_0_0_yy_xz_0_0, g_z_0_0_0_yy_yy_0_0, g_z_0_0_0_yy_yz_0_0, g_z_0_0_0_yy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yy_xx_0_0[i] = 2.0 * g_yyz_xx_0_0[i] * a_exp;

        g_z_0_0_0_yy_xy_0_0[i] = 2.0 * g_yyz_xy_0_0[i] * a_exp;

        g_z_0_0_0_yy_xz_0_0[i] = 2.0 * g_yyz_xz_0_0[i] * a_exp;

        g_z_0_0_0_yy_yy_0_0[i] = 2.0 * g_yyz_yy_0_0[i] * a_exp;

        g_z_0_0_0_yy_yz_0_0[i] = 2.0 * g_yyz_yz_0_0[i] * a_exp;

        g_z_0_0_0_yy_zz_0_0[i] = 2.0 * g_yyz_zz_0_0[i] * a_exp;
    }
    // integrals block (96-102)

    #pragma omp simd aligned(g_y_xx_0_0, g_y_xy_0_0, g_y_xz_0_0, g_y_yy_0_0, g_y_yz_0_0, g_y_zz_0_0, g_yzz_xx_0_0, g_yzz_xy_0_0, g_yzz_xz_0_0, g_yzz_yy_0_0, g_yzz_yz_0_0, g_yzz_zz_0_0, g_z_0_0_0_yz_xx_0_0, g_z_0_0_0_yz_xy_0_0, g_z_0_0_0_yz_xz_0_0, g_z_0_0_0_yz_yy_0_0, g_z_0_0_0_yz_yz_0_0, g_z_0_0_0_yz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_yz_xx_0_0[i] = -g_y_xx_0_0[i] + 2.0 * g_yzz_xx_0_0[i] * a_exp;

        g_z_0_0_0_yz_xy_0_0[i] = -g_y_xy_0_0[i] + 2.0 * g_yzz_xy_0_0[i] * a_exp;

        g_z_0_0_0_yz_xz_0_0[i] = -g_y_xz_0_0[i] + 2.0 * g_yzz_xz_0_0[i] * a_exp;

        g_z_0_0_0_yz_yy_0_0[i] = -g_y_yy_0_0[i] + 2.0 * g_yzz_yy_0_0[i] * a_exp;

        g_z_0_0_0_yz_yz_0_0[i] = -g_y_yz_0_0[i] + 2.0 * g_yzz_yz_0_0[i] * a_exp;

        g_z_0_0_0_yz_zz_0_0[i] = -g_y_zz_0_0[i] + 2.0 * g_yzz_zz_0_0[i] * a_exp;
    }
    // integrals block (102-108)

    #pragma omp simd aligned(g_z_0_0_0_zz_xx_0_0, g_z_0_0_0_zz_xy_0_0, g_z_0_0_0_zz_xz_0_0, g_z_0_0_0_zz_yy_0_0, g_z_0_0_0_zz_yz_0_0, g_z_0_0_0_zz_zz_0_0, g_z_xx_0_0, g_z_xy_0_0, g_z_xz_0_0, g_z_yy_0_0, g_z_yz_0_0, g_z_zz_0_0, g_zzz_xx_0_0, g_zzz_xy_0_0, g_zzz_xz_0_0, g_zzz_yy_0_0, g_zzz_yz_0_0, g_zzz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_zz_xx_0_0[i] = -2.0 * g_z_xx_0_0[i] + 2.0 * g_zzz_xx_0_0[i] * a_exp;

        g_z_0_0_0_zz_xy_0_0[i] = -2.0 * g_z_xy_0_0[i] + 2.0 * g_zzz_xy_0_0[i] * a_exp;

        g_z_0_0_0_zz_xz_0_0[i] = -2.0 * g_z_xz_0_0[i] + 2.0 * g_zzz_xz_0_0[i] * a_exp;

        g_z_0_0_0_zz_yy_0_0[i] = -2.0 * g_z_yy_0_0[i] + 2.0 * g_zzz_yy_0_0[i] * a_exp;

        g_z_0_0_0_zz_yz_0_0[i] = -2.0 * g_z_yz_0_0[i] + 2.0 * g_zzz_yz_0_0[i] * a_exp;

        g_z_0_0_0_zz_zz_0_0[i] = -2.0 * g_z_zz_0_0[i] + 2.0 * g_zzz_zz_0_0[i] * a_exp;
    }
}

} // t4c_geom namespace

