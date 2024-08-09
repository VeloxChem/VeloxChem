#include "GeomDeriv1100OfScalarForSDSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_sdss_0(CSimdArray<double>& buffer_1100_sdss,
                     const CSimdArray<double>& buffer_ppss,
                     const CSimdArray<double>& buffer_pfss,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_sdss.number_of_columns();

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

    /// Set up components of auxilary buffer : buffer_pfss

    auto g_x_xxx_0_0 = buffer_pfss[0];

    auto g_x_xxy_0_0 = buffer_pfss[1];

    auto g_x_xxz_0_0 = buffer_pfss[2];

    auto g_x_xyy_0_0 = buffer_pfss[3];

    auto g_x_xyz_0_0 = buffer_pfss[4];

    auto g_x_xzz_0_0 = buffer_pfss[5];

    auto g_x_yyy_0_0 = buffer_pfss[6];

    auto g_x_yyz_0_0 = buffer_pfss[7];

    auto g_x_yzz_0_0 = buffer_pfss[8];

    auto g_x_zzz_0_0 = buffer_pfss[9];

    auto g_y_xxx_0_0 = buffer_pfss[10];

    auto g_y_xxy_0_0 = buffer_pfss[11];

    auto g_y_xxz_0_0 = buffer_pfss[12];

    auto g_y_xyy_0_0 = buffer_pfss[13];

    auto g_y_xyz_0_0 = buffer_pfss[14];

    auto g_y_xzz_0_0 = buffer_pfss[15];

    auto g_y_yyy_0_0 = buffer_pfss[16];

    auto g_y_yyz_0_0 = buffer_pfss[17];

    auto g_y_yzz_0_0 = buffer_pfss[18];

    auto g_y_zzz_0_0 = buffer_pfss[19];

    auto g_z_xxx_0_0 = buffer_pfss[20];

    auto g_z_xxy_0_0 = buffer_pfss[21];

    auto g_z_xxz_0_0 = buffer_pfss[22];

    auto g_z_xyy_0_0 = buffer_pfss[23];

    auto g_z_xyz_0_0 = buffer_pfss[24];

    auto g_z_xzz_0_0 = buffer_pfss[25];

    auto g_z_yyy_0_0 = buffer_pfss[26];

    auto g_z_yyz_0_0 = buffer_pfss[27];

    auto g_z_yzz_0_0 = buffer_pfss[28];

    auto g_z_zzz_0_0 = buffer_pfss[29];

    /// Set up components of integrals buffer : buffer_1100_sdss

    auto g_x_x_0_0_0_xx_0_0 = buffer_1100_sdss[0];

    auto g_x_x_0_0_0_xy_0_0 = buffer_1100_sdss[1];

    auto g_x_x_0_0_0_xz_0_0 = buffer_1100_sdss[2];

    auto g_x_x_0_0_0_yy_0_0 = buffer_1100_sdss[3];

    auto g_x_x_0_0_0_yz_0_0 = buffer_1100_sdss[4];

    auto g_x_x_0_0_0_zz_0_0 = buffer_1100_sdss[5];

    auto g_x_y_0_0_0_xx_0_0 = buffer_1100_sdss[6];

    auto g_x_y_0_0_0_xy_0_0 = buffer_1100_sdss[7];

    auto g_x_y_0_0_0_xz_0_0 = buffer_1100_sdss[8];

    auto g_x_y_0_0_0_yy_0_0 = buffer_1100_sdss[9];

    auto g_x_y_0_0_0_yz_0_0 = buffer_1100_sdss[10];

    auto g_x_y_0_0_0_zz_0_0 = buffer_1100_sdss[11];

    auto g_x_z_0_0_0_xx_0_0 = buffer_1100_sdss[12];

    auto g_x_z_0_0_0_xy_0_0 = buffer_1100_sdss[13];

    auto g_x_z_0_0_0_xz_0_0 = buffer_1100_sdss[14];

    auto g_x_z_0_0_0_yy_0_0 = buffer_1100_sdss[15];

    auto g_x_z_0_0_0_yz_0_0 = buffer_1100_sdss[16];

    auto g_x_z_0_0_0_zz_0_0 = buffer_1100_sdss[17];

    auto g_y_x_0_0_0_xx_0_0 = buffer_1100_sdss[18];

    auto g_y_x_0_0_0_xy_0_0 = buffer_1100_sdss[19];

    auto g_y_x_0_0_0_xz_0_0 = buffer_1100_sdss[20];

    auto g_y_x_0_0_0_yy_0_0 = buffer_1100_sdss[21];

    auto g_y_x_0_0_0_yz_0_0 = buffer_1100_sdss[22];

    auto g_y_x_0_0_0_zz_0_0 = buffer_1100_sdss[23];

    auto g_y_y_0_0_0_xx_0_0 = buffer_1100_sdss[24];

    auto g_y_y_0_0_0_xy_0_0 = buffer_1100_sdss[25];

    auto g_y_y_0_0_0_xz_0_0 = buffer_1100_sdss[26];

    auto g_y_y_0_0_0_yy_0_0 = buffer_1100_sdss[27];

    auto g_y_y_0_0_0_yz_0_0 = buffer_1100_sdss[28];

    auto g_y_y_0_0_0_zz_0_0 = buffer_1100_sdss[29];

    auto g_y_z_0_0_0_xx_0_0 = buffer_1100_sdss[30];

    auto g_y_z_0_0_0_xy_0_0 = buffer_1100_sdss[31];

    auto g_y_z_0_0_0_xz_0_0 = buffer_1100_sdss[32];

    auto g_y_z_0_0_0_yy_0_0 = buffer_1100_sdss[33];

    auto g_y_z_0_0_0_yz_0_0 = buffer_1100_sdss[34];

    auto g_y_z_0_0_0_zz_0_0 = buffer_1100_sdss[35];

    auto g_z_x_0_0_0_xx_0_0 = buffer_1100_sdss[36];

    auto g_z_x_0_0_0_xy_0_0 = buffer_1100_sdss[37];

    auto g_z_x_0_0_0_xz_0_0 = buffer_1100_sdss[38];

    auto g_z_x_0_0_0_yy_0_0 = buffer_1100_sdss[39];

    auto g_z_x_0_0_0_yz_0_0 = buffer_1100_sdss[40];

    auto g_z_x_0_0_0_zz_0_0 = buffer_1100_sdss[41];

    auto g_z_y_0_0_0_xx_0_0 = buffer_1100_sdss[42];

    auto g_z_y_0_0_0_xy_0_0 = buffer_1100_sdss[43];

    auto g_z_y_0_0_0_xz_0_0 = buffer_1100_sdss[44];

    auto g_z_y_0_0_0_yy_0_0 = buffer_1100_sdss[45];

    auto g_z_y_0_0_0_yz_0_0 = buffer_1100_sdss[46];

    auto g_z_y_0_0_0_zz_0_0 = buffer_1100_sdss[47];

    auto g_z_z_0_0_0_xx_0_0 = buffer_1100_sdss[48];

    auto g_z_z_0_0_0_xy_0_0 = buffer_1100_sdss[49];

    auto g_z_z_0_0_0_xz_0_0 = buffer_1100_sdss[50];

    auto g_z_z_0_0_0_yy_0_0 = buffer_1100_sdss[51];

    auto g_z_z_0_0_0_yz_0_0 = buffer_1100_sdss[52];

    auto g_z_z_0_0_0_zz_0_0 = buffer_1100_sdss[53];

    // integrals block (0-6)

    #pragma omp simd aligned(g_x_x_0_0, g_x_x_0_0_0_xx_0_0, g_x_x_0_0_0_xy_0_0, g_x_x_0_0_0_xz_0_0, g_x_x_0_0_0_yy_0_0, g_x_x_0_0_0_yz_0_0, g_x_x_0_0_0_zz_0_0, g_x_xxx_0_0, g_x_xxy_0_0, g_x_xxz_0_0, g_x_xyy_0_0, g_x_xyz_0_0, g_x_xzz_0_0, g_x_y_0_0, g_x_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_xx_0_0[i] = -4.0 * g_x_x_0_0[i] * a_exp + 4.0 * g_x_xxx_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_0_xy_0_0[i] = -2.0 * g_x_y_0_0[i] * a_exp + 4.0 * g_x_xxy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_0_xz_0_0[i] = -2.0 * g_x_z_0_0[i] * a_exp + 4.0 * g_x_xxz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_0_yy_0_0[i] = 4.0 * g_x_xyy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_0_yz_0_0[i] = 4.0 * g_x_xyz_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_0_zz_0_0[i] = 4.0 * g_x_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (6-12)

    #pragma omp simd aligned(g_x_x_0_0, g_x_xxy_0_0, g_x_xyy_0_0, g_x_xyz_0_0, g_x_y_0_0, g_x_y_0_0_0_xx_0_0, g_x_y_0_0_0_xy_0_0, g_x_y_0_0_0_xz_0_0, g_x_y_0_0_0_yy_0_0, g_x_y_0_0_0_yz_0_0, g_x_y_0_0_0_zz_0_0, g_x_yyy_0_0, g_x_yyz_0_0, g_x_yzz_0_0, g_x_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_xx_0_0[i] = 4.0 * g_x_xxy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_0_xy_0_0[i] = -2.0 * g_x_x_0_0[i] * a_exp + 4.0 * g_x_xyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_0_xz_0_0[i] = 4.0 * g_x_xyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_0_yy_0_0[i] = -4.0 * g_x_y_0_0[i] * a_exp + 4.0 * g_x_yyy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_0_yz_0_0[i] = -2.0 * g_x_z_0_0[i] * a_exp + 4.0 * g_x_yyz_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_0_zz_0_0[i] = 4.0 * g_x_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (12-18)

    #pragma omp simd aligned(g_x_x_0_0, g_x_xxz_0_0, g_x_xyz_0_0, g_x_xzz_0_0, g_x_y_0_0, g_x_yyz_0_0, g_x_yzz_0_0, g_x_z_0_0, g_x_z_0_0_0_xx_0_0, g_x_z_0_0_0_xy_0_0, g_x_z_0_0_0_xz_0_0, g_x_z_0_0_0_yy_0_0, g_x_z_0_0_0_yz_0_0, g_x_z_0_0_0_zz_0_0, g_x_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_xx_0_0[i] = 4.0 * g_x_xxz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_0_xy_0_0[i] = 4.0 * g_x_xyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_0_xz_0_0[i] = -2.0 * g_x_x_0_0[i] * a_exp + 4.0 * g_x_xzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_0_yy_0_0[i] = 4.0 * g_x_yyz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_0_yz_0_0[i] = -2.0 * g_x_y_0_0[i] * a_exp + 4.0 * g_x_yzz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_0_zz_0_0[i] = -4.0 * g_x_z_0_0[i] * a_exp + 4.0 * g_x_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (18-24)

    #pragma omp simd aligned(g_y_x_0_0, g_y_x_0_0_0_xx_0_0, g_y_x_0_0_0_xy_0_0, g_y_x_0_0_0_xz_0_0, g_y_x_0_0_0_yy_0_0, g_y_x_0_0_0_yz_0_0, g_y_x_0_0_0_zz_0_0, g_y_xxx_0_0, g_y_xxy_0_0, g_y_xxz_0_0, g_y_xyy_0_0, g_y_xyz_0_0, g_y_xzz_0_0, g_y_y_0_0, g_y_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_xx_0_0[i] = -4.0 * g_y_x_0_0[i] * a_exp + 4.0 * g_y_xxx_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_0_xy_0_0[i] = -2.0 * g_y_y_0_0[i] * a_exp + 4.0 * g_y_xxy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_0_xz_0_0[i] = -2.0 * g_y_z_0_0[i] * a_exp + 4.0 * g_y_xxz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_0_yy_0_0[i] = 4.0 * g_y_xyy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_0_yz_0_0[i] = 4.0 * g_y_xyz_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_0_zz_0_0[i] = 4.0 * g_y_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (24-30)

    #pragma omp simd aligned(g_y_x_0_0, g_y_xxy_0_0, g_y_xyy_0_0, g_y_xyz_0_0, g_y_y_0_0, g_y_y_0_0_0_xx_0_0, g_y_y_0_0_0_xy_0_0, g_y_y_0_0_0_xz_0_0, g_y_y_0_0_0_yy_0_0, g_y_y_0_0_0_yz_0_0, g_y_y_0_0_0_zz_0_0, g_y_yyy_0_0, g_y_yyz_0_0, g_y_yzz_0_0, g_y_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_xx_0_0[i] = 4.0 * g_y_xxy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_0_xy_0_0[i] = -2.0 * g_y_x_0_0[i] * a_exp + 4.0 * g_y_xyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_0_xz_0_0[i] = 4.0 * g_y_xyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_0_yy_0_0[i] = -4.0 * g_y_y_0_0[i] * a_exp + 4.0 * g_y_yyy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_0_yz_0_0[i] = -2.0 * g_y_z_0_0[i] * a_exp + 4.0 * g_y_yyz_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_0_zz_0_0[i] = 4.0 * g_y_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (30-36)

    #pragma omp simd aligned(g_y_x_0_0, g_y_xxz_0_0, g_y_xyz_0_0, g_y_xzz_0_0, g_y_y_0_0, g_y_yyz_0_0, g_y_yzz_0_0, g_y_z_0_0, g_y_z_0_0_0_xx_0_0, g_y_z_0_0_0_xy_0_0, g_y_z_0_0_0_xz_0_0, g_y_z_0_0_0_yy_0_0, g_y_z_0_0_0_yz_0_0, g_y_z_0_0_0_zz_0_0, g_y_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_xx_0_0[i] = 4.0 * g_y_xxz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_0_xy_0_0[i] = 4.0 * g_y_xyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_0_xz_0_0[i] = -2.0 * g_y_x_0_0[i] * a_exp + 4.0 * g_y_xzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_0_yy_0_0[i] = 4.0 * g_y_yyz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_0_yz_0_0[i] = -2.0 * g_y_y_0_0[i] * a_exp + 4.0 * g_y_yzz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_0_zz_0_0[i] = -4.0 * g_y_z_0_0[i] * a_exp + 4.0 * g_y_zzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (36-42)

    #pragma omp simd aligned(g_z_x_0_0, g_z_x_0_0_0_xx_0_0, g_z_x_0_0_0_xy_0_0, g_z_x_0_0_0_xz_0_0, g_z_x_0_0_0_yy_0_0, g_z_x_0_0_0_yz_0_0, g_z_x_0_0_0_zz_0_0, g_z_xxx_0_0, g_z_xxy_0_0, g_z_xxz_0_0, g_z_xyy_0_0, g_z_xyz_0_0, g_z_xzz_0_0, g_z_y_0_0, g_z_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_xx_0_0[i] = -4.0 * g_z_x_0_0[i] * a_exp + 4.0 * g_z_xxx_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_0_xy_0_0[i] = -2.0 * g_z_y_0_0[i] * a_exp + 4.0 * g_z_xxy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_0_xz_0_0[i] = -2.0 * g_z_z_0_0[i] * a_exp + 4.0 * g_z_xxz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_0_yy_0_0[i] = 4.0 * g_z_xyy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_0_yz_0_0[i] = 4.0 * g_z_xyz_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_0_zz_0_0[i] = 4.0 * g_z_xzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (42-48)

    #pragma omp simd aligned(g_z_x_0_0, g_z_xxy_0_0, g_z_xyy_0_0, g_z_xyz_0_0, g_z_y_0_0, g_z_y_0_0_0_xx_0_0, g_z_y_0_0_0_xy_0_0, g_z_y_0_0_0_xz_0_0, g_z_y_0_0_0_yy_0_0, g_z_y_0_0_0_yz_0_0, g_z_y_0_0_0_zz_0_0, g_z_yyy_0_0, g_z_yyz_0_0, g_z_yzz_0_0, g_z_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_xx_0_0[i] = 4.0 * g_z_xxy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_0_xy_0_0[i] = -2.0 * g_z_x_0_0[i] * a_exp + 4.0 * g_z_xyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_0_xz_0_0[i] = 4.0 * g_z_xyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_0_yy_0_0[i] = -4.0 * g_z_y_0_0[i] * a_exp + 4.0 * g_z_yyy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_0_yz_0_0[i] = -2.0 * g_z_z_0_0[i] * a_exp + 4.0 * g_z_yyz_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_0_zz_0_0[i] = 4.0 * g_z_yzz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (48-54)

    #pragma omp simd aligned(g_z_x_0_0, g_z_xxz_0_0, g_z_xyz_0_0, g_z_xzz_0_0, g_z_y_0_0, g_z_yyz_0_0, g_z_yzz_0_0, g_z_z_0_0, g_z_z_0_0_0_xx_0_0, g_z_z_0_0_0_xy_0_0, g_z_z_0_0_0_xz_0_0, g_z_z_0_0_0_yy_0_0, g_z_z_0_0_0_yz_0_0, g_z_z_0_0_0_zz_0_0, g_z_zzz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_xx_0_0[i] = 4.0 * g_z_xxz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_0_xy_0_0[i] = 4.0 * g_z_xyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_0_xz_0_0[i] = -2.0 * g_z_x_0_0[i] * a_exp + 4.0 * g_z_xzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_0_yy_0_0[i] = 4.0 * g_z_yyz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_0_yz_0_0[i] = -2.0 * g_z_y_0_0[i] * a_exp + 4.0 * g_z_yzz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_0_zz_0_0[i] = -4.0 * g_z_z_0_0[i] * a_exp + 4.0 * g_z_zzz_0_0[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

