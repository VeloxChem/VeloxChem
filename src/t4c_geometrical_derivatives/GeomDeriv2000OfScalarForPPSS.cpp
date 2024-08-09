#include "GeomDeriv2000OfScalarForPPSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_ppss_0(CSimdArray<double>& buffer_2000_ppss,
                     const CSimdArray<double>& buffer_ppss,
                     const CSimdArray<double>& buffer_fpss,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_ppss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_ppss

    auto g_x_x_0_0 = buffer_ppss[0];

    auto g_x_y_0_0 = buffer_ppss[1];

    auto g_x_z_0_0 = buffer_ppss[2];

    auto g_y_x_0_0 = buffer_ppss[3];

    auto g_y_y_0_0 = buffer_ppss[4];

    auto g_y_z_0_0 = buffer_ppss[5];

    auto g_z_x_0_0 = buffer_ppss[6];

    auto g_z_y_0_0 = buffer_ppss[7];

    auto g_z_z_0_0 = buffer_ppss[8];

    /// Set up components of auxilary buffer : buffer_fpss

    auto g_xxx_x_0_0 = buffer_fpss[0];

    auto g_xxx_y_0_0 = buffer_fpss[1];

    auto g_xxx_z_0_0 = buffer_fpss[2];

    auto g_xxy_x_0_0 = buffer_fpss[3];

    auto g_xxy_y_0_0 = buffer_fpss[4];

    auto g_xxy_z_0_0 = buffer_fpss[5];

    auto g_xxz_x_0_0 = buffer_fpss[6];

    auto g_xxz_y_0_0 = buffer_fpss[7];

    auto g_xxz_z_0_0 = buffer_fpss[8];

    auto g_xyy_x_0_0 = buffer_fpss[9];

    auto g_xyy_y_0_0 = buffer_fpss[10];

    auto g_xyy_z_0_0 = buffer_fpss[11];

    auto g_xyz_x_0_0 = buffer_fpss[12];

    auto g_xyz_y_0_0 = buffer_fpss[13];

    auto g_xyz_z_0_0 = buffer_fpss[14];

    auto g_xzz_x_0_0 = buffer_fpss[15];

    auto g_xzz_y_0_0 = buffer_fpss[16];

    auto g_xzz_z_0_0 = buffer_fpss[17];

    auto g_yyy_x_0_0 = buffer_fpss[18];

    auto g_yyy_y_0_0 = buffer_fpss[19];

    auto g_yyy_z_0_0 = buffer_fpss[20];

    auto g_yyz_x_0_0 = buffer_fpss[21];

    auto g_yyz_y_0_0 = buffer_fpss[22];

    auto g_yyz_z_0_0 = buffer_fpss[23];

    auto g_yzz_x_0_0 = buffer_fpss[24];

    auto g_yzz_y_0_0 = buffer_fpss[25];

    auto g_yzz_z_0_0 = buffer_fpss[26];

    auto g_zzz_x_0_0 = buffer_fpss[27];

    auto g_zzz_y_0_0 = buffer_fpss[28];

    auto g_zzz_z_0_0 = buffer_fpss[29];

    /// Set up components of integrals buffer : buffer_2000_ppss

    auto g_xx_0_0_0_x_x_0_0 = buffer_2000_ppss[0];

    auto g_xx_0_0_0_x_y_0_0 = buffer_2000_ppss[1];

    auto g_xx_0_0_0_x_z_0_0 = buffer_2000_ppss[2];

    auto g_xx_0_0_0_y_x_0_0 = buffer_2000_ppss[3];

    auto g_xx_0_0_0_y_y_0_0 = buffer_2000_ppss[4];

    auto g_xx_0_0_0_y_z_0_0 = buffer_2000_ppss[5];

    auto g_xx_0_0_0_z_x_0_0 = buffer_2000_ppss[6];

    auto g_xx_0_0_0_z_y_0_0 = buffer_2000_ppss[7];

    auto g_xx_0_0_0_z_z_0_0 = buffer_2000_ppss[8];

    auto g_xy_0_0_0_x_x_0_0 = buffer_2000_ppss[9];

    auto g_xy_0_0_0_x_y_0_0 = buffer_2000_ppss[10];

    auto g_xy_0_0_0_x_z_0_0 = buffer_2000_ppss[11];

    auto g_xy_0_0_0_y_x_0_0 = buffer_2000_ppss[12];

    auto g_xy_0_0_0_y_y_0_0 = buffer_2000_ppss[13];

    auto g_xy_0_0_0_y_z_0_0 = buffer_2000_ppss[14];

    auto g_xy_0_0_0_z_x_0_0 = buffer_2000_ppss[15];

    auto g_xy_0_0_0_z_y_0_0 = buffer_2000_ppss[16];

    auto g_xy_0_0_0_z_z_0_0 = buffer_2000_ppss[17];

    auto g_xz_0_0_0_x_x_0_0 = buffer_2000_ppss[18];

    auto g_xz_0_0_0_x_y_0_0 = buffer_2000_ppss[19];

    auto g_xz_0_0_0_x_z_0_0 = buffer_2000_ppss[20];

    auto g_xz_0_0_0_y_x_0_0 = buffer_2000_ppss[21];

    auto g_xz_0_0_0_y_y_0_0 = buffer_2000_ppss[22];

    auto g_xz_0_0_0_y_z_0_0 = buffer_2000_ppss[23];

    auto g_xz_0_0_0_z_x_0_0 = buffer_2000_ppss[24];

    auto g_xz_0_0_0_z_y_0_0 = buffer_2000_ppss[25];

    auto g_xz_0_0_0_z_z_0_0 = buffer_2000_ppss[26];

    auto g_yy_0_0_0_x_x_0_0 = buffer_2000_ppss[27];

    auto g_yy_0_0_0_x_y_0_0 = buffer_2000_ppss[28];

    auto g_yy_0_0_0_x_z_0_0 = buffer_2000_ppss[29];

    auto g_yy_0_0_0_y_x_0_0 = buffer_2000_ppss[30];

    auto g_yy_0_0_0_y_y_0_0 = buffer_2000_ppss[31];

    auto g_yy_0_0_0_y_z_0_0 = buffer_2000_ppss[32];

    auto g_yy_0_0_0_z_x_0_0 = buffer_2000_ppss[33];

    auto g_yy_0_0_0_z_y_0_0 = buffer_2000_ppss[34];

    auto g_yy_0_0_0_z_z_0_0 = buffer_2000_ppss[35];

    auto g_yz_0_0_0_x_x_0_0 = buffer_2000_ppss[36];

    auto g_yz_0_0_0_x_y_0_0 = buffer_2000_ppss[37];

    auto g_yz_0_0_0_x_z_0_0 = buffer_2000_ppss[38];

    auto g_yz_0_0_0_y_x_0_0 = buffer_2000_ppss[39];

    auto g_yz_0_0_0_y_y_0_0 = buffer_2000_ppss[40];

    auto g_yz_0_0_0_y_z_0_0 = buffer_2000_ppss[41];

    auto g_yz_0_0_0_z_x_0_0 = buffer_2000_ppss[42];

    auto g_yz_0_0_0_z_y_0_0 = buffer_2000_ppss[43];

    auto g_yz_0_0_0_z_z_0_0 = buffer_2000_ppss[44];

    auto g_zz_0_0_0_x_x_0_0 = buffer_2000_ppss[45];

    auto g_zz_0_0_0_x_y_0_0 = buffer_2000_ppss[46];

    auto g_zz_0_0_0_x_z_0_0 = buffer_2000_ppss[47];

    auto g_zz_0_0_0_y_x_0_0 = buffer_2000_ppss[48];

    auto g_zz_0_0_0_y_y_0_0 = buffer_2000_ppss[49];

    auto g_zz_0_0_0_y_z_0_0 = buffer_2000_ppss[50];

    auto g_zz_0_0_0_z_x_0_0 = buffer_2000_ppss[51];

    auto g_zz_0_0_0_z_y_0_0 = buffer_2000_ppss[52];

    auto g_zz_0_0_0_z_z_0_0 = buffer_2000_ppss[53];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_x_0_0, g_x_y_0_0, g_x_z_0_0, g_xx_0_0_0_x_x_0_0, g_xx_0_0_0_x_y_0_0, g_xx_0_0_0_x_z_0_0, g_xxx_x_0_0, g_xxx_y_0_0, g_xxx_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_x_x_0_0[i] = -6.0 * g_x_x_0_0[i] * a_exp + 4.0 * g_xxx_x_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_x_y_0_0[i] = -6.0 * g_x_y_0_0[i] * a_exp + 4.0 * g_xxx_y_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_x_z_0_0[i] = -6.0 * g_x_z_0_0[i] * a_exp + 4.0 * g_xxx_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_xx_0_0_0_y_x_0_0, g_xx_0_0_0_y_y_0_0, g_xx_0_0_0_y_z_0_0, g_xxy_x_0_0, g_xxy_y_0_0, g_xxy_z_0_0, g_y_x_0_0, g_y_y_0_0, g_y_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_y_x_0_0[i] = -2.0 * g_y_x_0_0[i] * a_exp + 4.0 * g_xxy_x_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_y_y_0_0[i] = -2.0 * g_y_y_0_0[i] * a_exp + 4.0 * g_xxy_y_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_y_z_0_0[i] = -2.0 * g_y_z_0_0[i] * a_exp + 4.0 * g_xxy_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_xx_0_0_0_z_x_0_0, g_xx_0_0_0_z_y_0_0, g_xx_0_0_0_z_z_0_0, g_xxz_x_0_0, g_xxz_y_0_0, g_xxz_z_0_0, g_z_x_0_0, g_z_y_0_0, g_z_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_z_x_0_0[i] = -2.0 * g_z_x_0_0[i] * a_exp + 4.0 * g_xxz_x_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_z_y_0_0[i] = -2.0 * g_z_y_0_0[i] * a_exp + 4.0 * g_xxz_y_0_0[i] * a_exp * a_exp;

        g_xx_0_0_0_z_z_0_0[i] = -2.0 * g_z_z_0_0[i] * a_exp + 4.0 * g_xxz_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_xxy_x_0_0, g_xxy_y_0_0, g_xxy_z_0_0, g_xy_0_0_0_x_x_0_0, g_xy_0_0_0_x_y_0_0, g_xy_0_0_0_x_z_0_0, g_y_x_0_0, g_y_y_0_0, g_y_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_x_x_0_0[i] = -2.0 * g_y_x_0_0[i] * a_exp + 4.0 * g_xxy_x_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_x_y_0_0[i] = -2.0 * g_y_y_0_0[i] * a_exp + 4.0 * g_xxy_y_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_x_z_0_0[i] = -2.0 * g_y_z_0_0[i] * a_exp + 4.0 * g_xxy_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_x_0_0, g_x_y_0_0, g_x_z_0_0, g_xy_0_0_0_y_x_0_0, g_xy_0_0_0_y_y_0_0, g_xy_0_0_0_y_z_0_0, g_xyy_x_0_0, g_xyy_y_0_0, g_xyy_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_y_x_0_0[i] = -2.0 * g_x_x_0_0[i] * a_exp + 4.0 * g_xyy_x_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_y_y_0_0[i] = -2.0 * g_x_y_0_0[i] * a_exp + 4.0 * g_xyy_y_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_y_z_0_0[i] = -2.0 * g_x_z_0_0[i] * a_exp + 4.0 * g_xyy_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_xy_0_0_0_z_x_0_0, g_xy_0_0_0_z_y_0_0, g_xy_0_0_0_z_z_0_0, g_xyz_x_0_0, g_xyz_y_0_0, g_xyz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_z_x_0_0[i] = 4.0 * g_xyz_x_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_z_y_0_0[i] = 4.0 * g_xyz_y_0_0[i] * a_exp * a_exp;

        g_xy_0_0_0_z_z_0_0[i] = 4.0 * g_xyz_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_xxz_x_0_0, g_xxz_y_0_0, g_xxz_z_0_0, g_xz_0_0_0_x_x_0_0, g_xz_0_0_0_x_y_0_0, g_xz_0_0_0_x_z_0_0, g_z_x_0_0, g_z_y_0_0, g_z_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_x_x_0_0[i] = -2.0 * g_z_x_0_0[i] * a_exp + 4.0 * g_xxz_x_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_x_y_0_0[i] = -2.0 * g_z_y_0_0[i] * a_exp + 4.0 * g_xxz_y_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_x_z_0_0[i] = -2.0 * g_z_z_0_0[i] * a_exp + 4.0 * g_xxz_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_xyz_x_0_0, g_xyz_y_0_0, g_xyz_z_0_0, g_xz_0_0_0_y_x_0_0, g_xz_0_0_0_y_y_0_0, g_xz_0_0_0_y_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_y_x_0_0[i] = 4.0 * g_xyz_x_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_y_y_0_0[i] = 4.0 * g_xyz_y_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_y_z_0_0[i] = 4.0 * g_xyz_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_x_0_0, g_x_y_0_0, g_x_z_0_0, g_xz_0_0_0_z_x_0_0, g_xz_0_0_0_z_y_0_0, g_xz_0_0_0_z_z_0_0, g_xzz_x_0_0, g_xzz_y_0_0, g_xzz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_z_x_0_0[i] = -2.0 * g_x_x_0_0[i] * a_exp + 4.0 * g_xzz_x_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_z_y_0_0[i] = -2.0 * g_x_y_0_0[i] * a_exp + 4.0 * g_xzz_y_0_0[i] * a_exp * a_exp;

        g_xz_0_0_0_z_z_0_0[i] = -2.0 * g_x_z_0_0[i] * a_exp + 4.0 * g_xzz_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_x_x_0_0, g_x_y_0_0, g_x_z_0_0, g_xyy_x_0_0, g_xyy_y_0_0, g_xyy_z_0_0, g_yy_0_0_0_x_x_0_0, g_yy_0_0_0_x_y_0_0, g_yy_0_0_0_x_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_x_x_0_0[i] = -2.0 * g_x_x_0_0[i] * a_exp + 4.0 * g_xyy_x_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_x_y_0_0[i] = -2.0 * g_x_y_0_0[i] * a_exp + 4.0 * g_xyy_y_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_x_z_0_0[i] = -2.0 * g_x_z_0_0[i] * a_exp + 4.0 * g_xyy_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_y_x_0_0, g_y_y_0_0, g_y_z_0_0, g_yy_0_0_0_y_x_0_0, g_yy_0_0_0_y_y_0_0, g_yy_0_0_0_y_z_0_0, g_yyy_x_0_0, g_yyy_y_0_0, g_yyy_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_y_x_0_0[i] = -6.0 * g_y_x_0_0[i] * a_exp + 4.0 * g_yyy_x_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_y_y_0_0[i] = -6.0 * g_y_y_0_0[i] * a_exp + 4.0 * g_yyy_y_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_y_z_0_0[i] = -6.0 * g_y_z_0_0[i] * a_exp + 4.0 * g_yyy_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_yy_0_0_0_z_x_0_0, g_yy_0_0_0_z_y_0_0, g_yy_0_0_0_z_z_0_0, g_yyz_x_0_0, g_yyz_y_0_0, g_yyz_z_0_0, g_z_x_0_0, g_z_y_0_0, g_z_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_z_x_0_0[i] = -2.0 * g_z_x_0_0[i] * a_exp + 4.0 * g_yyz_x_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_z_y_0_0[i] = -2.0 * g_z_y_0_0[i] * a_exp + 4.0 * g_yyz_y_0_0[i] * a_exp * a_exp;

        g_yy_0_0_0_z_z_0_0[i] = -2.0 * g_z_z_0_0[i] * a_exp + 4.0 * g_yyz_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_xyz_x_0_0, g_xyz_y_0_0, g_xyz_z_0_0, g_yz_0_0_0_x_x_0_0, g_yz_0_0_0_x_y_0_0, g_yz_0_0_0_x_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_x_x_0_0[i] = 4.0 * g_xyz_x_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_x_y_0_0[i] = 4.0 * g_xyz_y_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_x_z_0_0[i] = 4.0 * g_xyz_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_yyz_x_0_0, g_yyz_y_0_0, g_yyz_z_0_0, g_yz_0_0_0_y_x_0_0, g_yz_0_0_0_y_y_0_0, g_yz_0_0_0_y_z_0_0, g_z_x_0_0, g_z_y_0_0, g_z_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_y_x_0_0[i] = -2.0 * g_z_x_0_0[i] * a_exp + 4.0 * g_yyz_x_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_y_y_0_0[i] = -2.0 * g_z_y_0_0[i] * a_exp + 4.0 * g_yyz_y_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_y_z_0_0[i] = -2.0 * g_z_z_0_0[i] * a_exp + 4.0 * g_yyz_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_y_x_0_0, g_y_y_0_0, g_y_z_0_0, g_yz_0_0_0_z_x_0_0, g_yz_0_0_0_z_y_0_0, g_yz_0_0_0_z_z_0_0, g_yzz_x_0_0, g_yzz_y_0_0, g_yzz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_z_x_0_0[i] = -2.0 * g_y_x_0_0[i] * a_exp + 4.0 * g_yzz_x_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_z_y_0_0[i] = -2.0 * g_y_y_0_0[i] * a_exp + 4.0 * g_yzz_y_0_0[i] * a_exp * a_exp;

        g_yz_0_0_0_z_z_0_0[i] = -2.0 * g_y_z_0_0[i] * a_exp + 4.0 * g_yzz_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_x_x_0_0, g_x_y_0_0, g_x_z_0_0, g_xzz_x_0_0, g_xzz_y_0_0, g_xzz_z_0_0, g_zz_0_0_0_x_x_0_0, g_zz_0_0_0_x_y_0_0, g_zz_0_0_0_x_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_x_x_0_0[i] = -2.0 * g_x_x_0_0[i] * a_exp + 4.0 * g_xzz_x_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_x_y_0_0[i] = -2.0 * g_x_y_0_0[i] * a_exp + 4.0 * g_xzz_y_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_x_z_0_0[i] = -2.0 * g_x_z_0_0[i] * a_exp + 4.0 * g_xzz_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_y_x_0_0, g_y_y_0_0, g_y_z_0_0, g_yzz_x_0_0, g_yzz_y_0_0, g_yzz_z_0_0, g_zz_0_0_0_y_x_0_0, g_zz_0_0_0_y_y_0_0, g_zz_0_0_0_y_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_y_x_0_0[i] = -2.0 * g_y_x_0_0[i] * a_exp + 4.0 * g_yzz_x_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_y_y_0_0[i] = -2.0 * g_y_y_0_0[i] * a_exp + 4.0 * g_yzz_y_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_y_z_0_0[i] = -2.0 * g_y_z_0_0[i] * a_exp + 4.0 * g_yzz_z_0_0[i] * a_exp * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_z_x_0_0, g_z_y_0_0, g_z_z_0_0, g_zz_0_0_0_z_x_0_0, g_zz_0_0_0_z_y_0_0, g_zz_0_0_0_z_z_0_0, g_zzz_x_0_0, g_zzz_y_0_0, g_zzz_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_z_x_0_0[i] = -6.0 * g_z_x_0_0[i] * a_exp + 4.0 * g_zzz_x_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_z_y_0_0[i] = -6.0 * g_z_y_0_0[i] * a_exp + 4.0 * g_zzz_y_0_0[i] * a_exp * a_exp;

        g_zz_0_0_0_z_z_0_0[i] = -6.0 * g_z_z_0_0[i] * a_exp + 4.0 * g_zzz_z_0_0[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

