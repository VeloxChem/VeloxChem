#include "GeomDeriv1100OfScalarForPPSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_ppss_0(CSimdArray<double>& buffer_1100_ppss,
                     const CSimdArray<double>& buffer_ssss,
                     const CSimdArray<double>& buffer_sdss,
                     const CSimdArray<double>& buffer_dsss,
                     const CSimdArray<double>& buffer_ddss,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_ppss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_ssss

    auto g_0_0_0_0 = buffer_ssss[0];

    /// Set up components of auxilary buffer : buffer_sdss

    auto g_0_xx_0_0 = buffer_sdss[0];

    auto g_0_xy_0_0 = buffer_sdss[1];

    auto g_0_xz_0_0 = buffer_sdss[2];

    auto g_0_yy_0_0 = buffer_sdss[3];

    auto g_0_yz_0_0 = buffer_sdss[4];

    auto g_0_zz_0_0 = buffer_sdss[5];

    /// Set up components of auxilary buffer : buffer_dsss

    auto g_xx_0_0_0 = buffer_dsss[0];

    auto g_xy_0_0_0 = buffer_dsss[1];

    auto g_xz_0_0_0 = buffer_dsss[2];

    auto g_yy_0_0_0 = buffer_dsss[3];

    auto g_yz_0_0_0 = buffer_dsss[4];

    auto g_zz_0_0_0 = buffer_dsss[5];

    /// Set up components of auxilary buffer : buffer_ddss

    auto g_xx_xx_0_0 = buffer_ddss[0];

    auto g_xx_xy_0_0 = buffer_ddss[1];

    auto g_xx_xz_0_0 = buffer_ddss[2];

    auto g_xx_yy_0_0 = buffer_ddss[3];

    auto g_xx_yz_0_0 = buffer_ddss[4];

    auto g_xx_zz_0_0 = buffer_ddss[5];

    auto g_xy_xx_0_0 = buffer_ddss[6];

    auto g_xy_xy_0_0 = buffer_ddss[7];

    auto g_xy_xz_0_0 = buffer_ddss[8];

    auto g_xy_yy_0_0 = buffer_ddss[9];

    auto g_xy_yz_0_0 = buffer_ddss[10];

    auto g_xy_zz_0_0 = buffer_ddss[11];

    auto g_xz_xx_0_0 = buffer_ddss[12];

    auto g_xz_xy_0_0 = buffer_ddss[13];

    auto g_xz_xz_0_0 = buffer_ddss[14];

    auto g_xz_yy_0_0 = buffer_ddss[15];

    auto g_xz_yz_0_0 = buffer_ddss[16];

    auto g_xz_zz_0_0 = buffer_ddss[17];

    auto g_yy_xx_0_0 = buffer_ddss[18];

    auto g_yy_xy_0_0 = buffer_ddss[19];

    auto g_yy_xz_0_0 = buffer_ddss[20];

    auto g_yy_yy_0_0 = buffer_ddss[21];

    auto g_yy_yz_0_0 = buffer_ddss[22];

    auto g_yy_zz_0_0 = buffer_ddss[23];

    auto g_yz_xx_0_0 = buffer_ddss[24];

    auto g_yz_xy_0_0 = buffer_ddss[25];

    auto g_yz_xz_0_0 = buffer_ddss[26];

    auto g_yz_yy_0_0 = buffer_ddss[27];

    auto g_yz_yz_0_0 = buffer_ddss[28];

    auto g_yz_zz_0_0 = buffer_ddss[29];

    auto g_zz_xx_0_0 = buffer_ddss[30];

    auto g_zz_xy_0_0 = buffer_ddss[31];

    auto g_zz_xz_0_0 = buffer_ddss[32];

    auto g_zz_yy_0_0 = buffer_ddss[33];

    auto g_zz_yz_0_0 = buffer_ddss[34];

    auto g_zz_zz_0_0 = buffer_ddss[35];

    /// Set up components of integrals buffer : buffer_1100_ppss

    auto g_x_x_0_0_x_x_0_0 = buffer_1100_ppss[0];

    auto g_x_x_0_0_x_y_0_0 = buffer_1100_ppss[1];

    auto g_x_x_0_0_x_z_0_0 = buffer_1100_ppss[2];

    auto g_x_x_0_0_y_x_0_0 = buffer_1100_ppss[3];

    auto g_x_x_0_0_y_y_0_0 = buffer_1100_ppss[4];

    auto g_x_x_0_0_y_z_0_0 = buffer_1100_ppss[5];

    auto g_x_x_0_0_z_x_0_0 = buffer_1100_ppss[6];

    auto g_x_x_0_0_z_y_0_0 = buffer_1100_ppss[7];

    auto g_x_x_0_0_z_z_0_0 = buffer_1100_ppss[8];

    auto g_x_y_0_0_x_x_0_0 = buffer_1100_ppss[9];

    auto g_x_y_0_0_x_y_0_0 = buffer_1100_ppss[10];

    auto g_x_y_0_0_x_z_0_0 = buffer_1100_ppss[11];

    auto g_x_y_0_0_y_x_0_0 = buffer_1100_ppss[12];

    auto g_x_y_0_0_y_y_0_0 = buffer_1100_ppss[13];

    auto g_x_y_0_0_y_z_0_0 = buffer_1100_ppss[14];

    auto g_x_y_0_0_z_x_0_0 = buffer_1100_ppss[15];

    auto g_x_y_0_0_z_y_0_0 = buffer_1100_ppss[16];

    auto g_x_y_0_0_z_z_0_0 = buffer_1100_ppss[17];

    auto g_x_z_0_0_x_x_0_0 = buffer_1100_ppss[18];

    auto g_x_z_0_0_x_y_0_0 = buffer_1100_ppss[19];

    auto g_x_z_0_0_x_z_0_0 = buffer_1100_ppss[20];

    auto g_x_z_0_0_y_x_0_0 = buffer_1100_ppss[21];

    auto g_x_z_0_0_y_y_0_0 = buffer_1100_ppss[22];

    auto g_x_z_0_0_y_z_0_0 = buffer_1100_ppss[23];

    auto g_x_z_0_0_z_x_0_0 = buffer_1100_ppss[24];

    auto g_x_z_0_0_z_y_0_0 = buffer_1100_ppss[25];

    auto g_x_z_0_0_z_z_0_0 = buffer_1100_ppss[26];

    auto g_y_x_0_0_x_x_0_0 = buffer_1100_ppss[27];

    auto g_y_x_0_0_x_y_0_0 = buffer_1100_ppss[28];

    auto g_y_x_0_0_x_z_0_0 = buffer_1100_ppss[29];

    auto g_y_x_0_0_y_x_0_0 = buffer_1100_ppss[30];

    auto g_y_x_0_0_y_y_0_0 = buffer_1100_ppss[31];

    auto g_y_x_0_0_y_z_0_0 = buffer_1100_ppss[32];

    auto g_y_x_0_0_z_x_0_0 = buffer_1100_ppss[33];

    auto g_y_x_0_0_z_y_0_0 = buffer_1100_ppss[34];

    auto g_y_x_0_0_z_z_0_0 = buffer_1100_ppss[35];

    auto g_y_y_0_0_x_x_0_0 = buffer_1100_ppss[36];

    auto g_y_y_0_0_x_y_0_0 = buffer_1100_ppss[37];

    auto g_y_y_0_0_x_z_0_0 = buffer_1100_ppss[38];

    auto g_y_y_0_0_y_x_0_0 = buffer_1100_ppss[39];

    auto g_y_y_0_0_y_y_0_0 = buffer_1100_ppss[40];

    auto g_y_y_0_0_y_z_0_0 = buffer_1100_ppss[41];

    auto g_y_y_0_0_z_x_0_0 = buffer_1100_ppss[42];

    auto g_y_y_0_0_z_y_0_0 = buffer_1100_ppss[43];

    auto g_y_y_0_0_z_z_0_0 = buffer_1100_ppss[44];

    auto g_y_z_0_0_x_x_0_0 = buffer_1100_ppss[45];

    auto g_y_z_0_0_x_y_0_0 = buffer_1100_ppss[46];

    auto g_y_z_0_0_x_z_0_0 = buffer_1100_ppss[47];

    auto g_y_z_0_0_y_x_0_0 = buffer_1100_ppss[48];

    auto g_y_z_0_0_y_y_0_0 = buffer_1100_ppss[49];

    auto g_y_z_0_0_y_z_0_0 = buffer_1100_ppss[50];

    auto g_y_z_0_0_z_x_0_0 = buffer_1100_ppss[51];

    auto g_y_z_0_0_z_y_0_0 = buffer_1100_ppss[52];

    auto g_y_z_0_0_z_z_0_0 = buffer_1100_ppss[53];

    auto g_z_x_0_0_x_x_0_0 = buffer_1100_ppss[54];

    auto g_z_x_0_0_x_y_0_0 = buffer_1100_ppss[55];

    auto g_z_x_0_0_x_z_0_0 = buffer_1100_ppss[56];

    auto g_z_x_0_0_y_x_0_0 = buffer_1100_ppss[57];

    auto g_z_x_0_0_y_y_0_0 = buffer_1100_ppss[58];

    auto g_z_x_0_0_y_z_0_0 = buffer_1100_ppss[59];

    auto g_z_x_0_0_z_x_0_0 = buffer_1100_ppss[60];

    auto g_z_x_0_0_z_y_0_0 = buffer_1100_ppss[61];

    auto g_z_x_0_0_z_z_0_0 = buffer_1100_ppss[62];

    auto g_z_y_0_0_x_x_0_0 = buffer_1100_ppss[63];

    auto g_z_y_0_0_x_y_0_0 = buffer_1100_ppss[64];

    auto g_z_y_0_0_x_z_0_0 = buffer_1100_ppss[65];

    auto g_z_y_0_0_y_x_0_0 = buffer_1100_ppss[66];

    auto g_z_y_0_0_y_y_0_0 = buffer_1100_ppss[67];

    auto g_z_y_0_0_y_z_0_0 = buffer_1100_ppss[68];

    auto g_z_y_0_0_z_x_0_0 = buffer_1100_ppss[69];

    auto g_z_y_0_0_z_y_0_0 = buffer_1100_ppss[70];

    auto g_z_y_0_0_z_z_0_0 = buffer_1100_ppss[71];

    auto g_z_z_0_0_x_x_0_0 = buffer_1100_ppss[72];

    auto g_z_z_0_0_x_y_0_0 = buffer_1100_ppss[73];

    auto g_z_z_0_0_x_z_0_0 = buffer_1100_ppss[74];

    auto g_z_z_0_0_y_x_0_0 = buffer_1100_ppss[75];

    auto g_z_z_0_0_y_y_0_0 = buffer_1100_ppss[76];

    auto g_z_z_0_0_y_z_0_0 = buffer_1100_ppss[77];

    auto g_z_z_0_0_z_x_0_0 = buffer_1100_ppss[78];

    auto g_z_z_0_0_z_y_0_0 = buffer_1100_ppss[79];

    auto g_z_z_0_0_z_z_0_0 = buffer_1100_ppss[80];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_0_0_0, g_0_xx_0_0, g_0_xy_0_0, g_0_xz_0_0, g_x_x_0_0_x_x_0_0, g_x_x_0_0_x_y_0_0, g_x_x_0_0_x_z_0_0, g_xx_0_0_0, g_xx_xx_0_0, g_xx_xy_0_0, g_xx_xz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_x_x_0_0[i] = g_0_0_0_0[i] - 2.0 * g_0_xx_0_0[i] * b_exp - 2.0 * g_xx_0_0_0[i] * a_exp + 4.0 * g_xx_xx_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_x_y_0_0[i] = -2.0 * g_0_xy_0_0[i] * b_exp + 4.0 * g_xx_xy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_x_z_0_0[i] = -2.0 * g_0_xz_0_0[i] * b_exp + 4.0 * g_xx_xz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_x_0_0_y_x_0_0, g_x_x_0_0_y_y_0_0, g_x_x_0_0_y_z_0_0, g_xy_0_0_0, g_xy_xx_0_0, g_xy_xy_0_0, g_xy_xz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_y_x_0_0[i] = -2.0 * g_xy_0_0_0[i] * a_exp + 4.0 * g_xy_xx_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_y_y_0_0[i] = 4.0 * g_xy_xy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_y_z_0_0[i] = 4.0 * g_xy_xz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_x_0_0_z_x_0_0, g_x_x_0_0_z_y_0_0, g_x_x_0_0_z_z_0_0, g_xz_0_0_0, g_xz_xx_0_0, g_xz_xy_0_0, g_xz_xz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_z_x_0_0[i] = -2.0 * g_xz_0_0_0[i] * a_exp + 4.0 * g_xz_xx_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_z_y_0_0[i] = 4.0 * g_xz_xy_0_0[i] * a_exp * b_exp;

        g_x_x_0_0_z_z_0_0[i] = 4.0 * g_xz_xz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_0_0_0_0, g_0_xy_0_0, g_0_yy_0_0, g_0_yz_0_0, g_x_y_0_0_x_x_0_0, g_x_y_0_0_x_y_0_0, g_x_y_0_0_x_z_0_0, g_xx_0_0_0, g_xx_xy_0_0, g_xx_yy_0_0, g_xx_yz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_x_x_0_0[i] = -2.0 * g_0_xy_0_0[i] * b_exp + 4.0 * g_xx_xy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_x_y_0_0[i] = g_0_0_0_0[i] - 2.0 * g_0_yy_0_0[i] * b_exp - 2.0 * g_xx_0_0_0[i] * a_exp + 4.0 * g_xx_yy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_x_z_0_0[i] = -2.0 * g_0_yz_0_0[i] * b_exp + 4.0 * g_xx_yz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_y_0_0_y_x_0_0, g_x_y_0_0_y_y_0_0, g_x_y_0_0_y_z_0_0, g_xy_0_0_0, g_xy_xy_0_0, g_xy_yy_0_0, g_xy_yz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_y_x_0_0[i] = 4.0 * g_xy_xy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_y_y_0_0[i] = -2.0 * g_xy_0_0_0[i] * a_exp + 4.0 * g_xy_yy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_y_z_0_0[i] = 4.0 * g_xy_yz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_y_0_0_z_x_0_0, g_x_y_0_0_z_y_0_0, g_x_y_0_0_z_z_0_0, g_xz_0_0_0, g_xz_xy_0_0, g_xz_yy_0_0, g_xz_yz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_z_x_0_0[i] = 4.0 * g_xz_xy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_z_y_0_0[i] = -2.0 * g_xz_0_0_0[i] * a_exp + 4.0 * g_xz_yy_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_z_z_0_0[i] = 4.0 * g_xz_yz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_0_0_0_0, g_0_xz_0_0, g_0_yz_0_0, g_0_zz_0_0, g_x_z_0_0_x_x_0_0, g_x_z_0_0_x_y_0_0, g_x_z_0_0_x_z_0_0, g_xx_0_0_0, g_xx_xz_0_0, g_xx_yz_0_0, g_xx_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_x_x_0_0[i] = -2.0 * g_0_xz_0_0[i] * b_exp + 4.0 * g_xx_xz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_x_y_0_0[i] = -2.0 * g_0_yz_0_0[i] * b_exp + 4.0 * g_xx_yz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_x_z_0_0[i] = g_0_0_0_0[i] - 2.0 * g_0_zz_0_0[i] * b_exp - 2.0 * g_xx_0_0_0[i] * a_exp + 4.0 * g_xx_zz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_z_0_0_y_x_0_0, g_x_z_0_0_y_y_0_0, g_x_z_0_0_y_z_0_0, g_xy_0_0_0, g_xy_xz_0_0, g_xy_yz_0_0, g_xy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_y_x_0_0[i] = 4.0 * g_xy_xz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_y_y_0_0[i] = 4.0 * g_xy_yz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_y_z_0_0[i] = -2.0 * g_xy_0_0_0[i] * a_exp + 4.0 * g_xy_zz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_z_0_0_z_x_0_0, g_x_z_0_0_z_y_0_0, g_x_z_0_0_z_z_0_0, g_xz_0_0_0, g_xz_xz_0_0, g_xz_yz_0_0, g_xz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_z_x_0_0[i] = 4.0 * g_xz_xz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_z_y_0_0[i] = 4.0 * g_xz_yz_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_z_z_0_0[i] = -2.0 * g_xz_0_0_0[i] * a_exp + 4.0 * g_xz_zz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_xy_0_0_0, g_xy_xx_0_0, g_xy_xy_0_0, g_xy_xz_0_0, g_y_x_0_0_x_x_0_0, g_y_x_0_0_x_y_0_0, g_y_x_0_0_x_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_x_x_0_0[i] = -2.0 * g_xy_0_0_0[i] * a_exp + 4.0 * g_xy_xx_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_x_y_0_0[i] = 4.0 * g_xy_xy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_x_z_0_0[i] = 4.0 * g_xy_xz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_0_0_0_0, g_0_xx_0_0, g_0_xy_0_0, g_0_xz_0_0, g_y_x_0_0_y_x_0_0, g_y_x_0_0_y_y_0_0, g_y_x_0_0_y_z_0_0, g_yy_0_0_0, g_yy_xx_0_0, g_yy_xy_0_0, g_yy_xz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_y_x_0_0[i] = g_0_0_0_0[i] - 2.0 * g_0_xx_0_0[i] * b_exp - 2.0 * g_yy_0_0_0[i] * a_exp + 4.0 * g_yy_xx_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_y_y_0_0[i] = -2.0 * g_0_xy_0_0[i] * b_exp + 4.0 * g_yy_xy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_y_z_0_0[i] = -2.0 * g_0_xz_0_0[i] * b_exp + 4.0 * g_yy_xz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_y_x_0_0_z_x_0_0, g_y_x_0_0_z_y_0_0, g_y_x_0_0_z_z_0_0, g_yz_0_0_0, g_yz_xx_0_0, g_yz_xy_0_0, g_yz_xz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_z_x_0_0[i] = -2.0 * g_yz_0_0_0[i] * a_exp + 4.0 * g_yz_xx_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_z_y_0_0[i] = 4.0 * g_yz_xy_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_z_z_0_0[i] = 4.0 * g_yz_xz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_xy_0_0_0, g_xy_xy_0_0, g_xy_yy_0_0, g_xy_yz_0_0, g_y_y_0_0_x_x_0_0, g_y_y_0_0_x_y_0_0, g_y_y_0_0_x_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_x_x_0_0[i] = 4.0 * g_xy_xy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_x_y_0_0[i] = -2.0 * g_xy_0_0_0[i] * a_exp + 4.0 * g_xy_yy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_x_z_0_0[i] = 4.0 * g_xy_yz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_0_0_0_0, g_0_xy_0_0, g_0_yy_0_0, g_0_yz_0_0, g_y_y_0_0_y_x_0_0, g_y_y_0_0_y_y_0_0, g_y_y_0_0_y_z_0_0, g_yy_0_0_0, g_yy_xy_0_0, g_yy_yy_0_0, g_yy_yz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_y_x_0_0[i] = -2.0 * g_0_xy_0_0[i] * b_exp + 4.0 * g_yy_xy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_y_y_0_0[i] = g_0_0_0_0[i] - 2.0 * g_0_yy_0_0[i] * b_exp - 2.0 * g_yy_0_0_0[i] * a_exp + 4.0 * g_yy_yy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_y_z_0_0[i] = -2.0 * g_0_yz_0_0[i] * b_exp + 4.0 * g_yy_yz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_y_y_0_0_z_x_0_0, g_y_y_0_0_z_y_0_0, g_y_y_0_0_z_z_0_0, g_yz_0_0_0, g_yz_xy_0_0, g_yz_yy_0_0, g_yz_yz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_z_x_0_0[i] = 4.0 * g_yz_xy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_z_y_0_0[i] = -2.0 * g_yz_0_0_0[i] * a_exp + 4.0 * g_yz_yy_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_z_z_0_0[i] = 4.0 * g_yz_yz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_xy_0_0_0, g_xy_xz_0_0, g_xy_yz_0_0, g_xy_zz_0_0, g_y_z_0_0_x_x_0_0, g_y_z_0_0_x_y_0_0, g_y_z_0_0_x_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_x_x_0_0[i] = 4.0 * g_xy_xz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_x_y_0_0[i] = 4.0 * g_xy_yz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_x_z_0_0[i] = -2.0 * g_xy_0_0_0[i] * a_exp + 4.0 * g_xy_zz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_0_0_0_0, g_0_xz_0_0, g_0_yz_0_0, g_0_zz_0_0, g_y_z_0_0_y_x_0_0, g_y_z_0_0_y_y_0_0, g_y_z_0_0_y_z_0_0, g_yy_0_0_0, g_yy_xz_0_0, g_yy_yz_0_0, g_yy_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_y_x_0_0[i] = -2.0 * g_0_xz_0_0[i] * b_exp + 4.0 * g_yy_xz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_y_y_0_0[i] = -2.0 * g_0_yz_0_0[i] * b_exp + 4.0 * g_yy_yz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_y_z_0_0[i] = g_0_0_0_0[i] - 2.0 * g_0_zz_0_0[i] * b_exp - 2.0 * g_yy_0_0_0[i] * a_exp + 4.0 * g_yy_zz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_y_z_0_0_z_x_0_0, g_y_z_0_0_z_y_0_0, g_y_z_0_0_z_z_0_0, g_yz_0_0_0, g_yz_xz_0_0, g_yz_yz_0_0, g_yz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_z_x_0_0[i] = 4.0 * g_yz_xz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_z_y_0_0[i] = 4.0 * g_yz_yz_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_z_z_0_0[i] = -2.0 * g_yz_0_0_0[i] * a_exp + 4.0 * g_yz_zz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_xz_0_0_0, g_xz_xx_0_0, g_xz_xy_0_0, g_xz_xz_0_0, g_z_x_0_0_x_x_0_0, g_z_x_0_0_x_y_0_0, g_z_x_0_0_x_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_x_x_0_0[i] = -2.0 * g_xz_0_0_0[i] * a_exp + 4.0 * g_xz_xx_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_x_y_0_0[i] = 4.0 * g_xz_xy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_x_z_0_0[i] = 4.0 * g_xz_xz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_yz_0_0_0, g_yz_xx_0_0, g_yz_xy_0_0, g_yz_xz_0_0, g_z_x_0_0_y_x_0_0, g_z_x_0_0_y_y_0_0, g_z_x_0_0_y_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_y_x_0_0[i] = -2.0 * g_yz_0_0_0[i] * a_exp + 4.0 * g_yz_xx_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_y_y_0_0[i] = 4.0 * g_yz_xy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_y_z_0_0[i] = 4.0 * g_yz_xz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_0_0_0_0, g_0_xx_0_0, g_0_xy_0_0, g_0_xz_0_0, g_z_x_0_0_z_x_0_0, g_z_x_0_0_z_y_0_0, g_z_x_0_0_z_z_0_0, g_zz_0_0_0, g_zz_xx_0_0, g_zz_xy_0_0, g_zz_xz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_z_x_0_0[i] = g_0_0_0_0[i] - 2.0 * g_0_xx_0_0[i] * b_exp - 2.0 * g_zz_0_0_0[i] * a_exp + 4.0 * g_zz_xx_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_z_y_0_0[i] = -2.0 * g_0_xy_0_0[i] * b_exp + 4.0 * g_zz_xy_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_z_z_0_0[i] = -2.0 * g_0_xz_0_0[i] * b_exp + 4.0 * g_zz_xz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_xz_0_0_0, g_xz_xy_0_0, g_xz_yy_0_0, g_xz_yz_0_0, g_z_y_0_0_x_x_0_0, g_z_y_0_0_x_y_0_0, g_z_y_0_0_x_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_x_x_0_0[i] = 4.0 * g_xz_xy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_x_y_0_0[i] = -2.0 * g_xz_0_0_0[i] * a_exp + 4.0 * g_xz_yy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_x_z_0_0[i] = 4.0 * g_xz_yz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_yz_0_0_0, g_yz_xy_0_0, g_yz_yy_0_0, g_yz_yz_0_0, g_z_y_0_0_y_x_0_0, g_z_y_0_0_y_y_0_0, g_z_y_0_0_y_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_y_x_0_0[i] = 4.0 * g_yz_xy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_y_y_0_0[i] = -2.0 * g_yz_0_0_0[i] * a_exp + 4.0 * g_yz_yy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_y_z_0_0[i] = 4.0 * g_yz_yz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_0_0_0_0, g_0_xy_0_0, g_0_yy_0_0, g_0_yz_0_0, g_z_y_0_0_z_x_0_0, g_z_y_0_0_z_y_0_0, g_z_y_0_0_z_z_0_0, g_zz_0_0_0, g_zz_xy_0_0, g_zz_yy_0_0, g_zz_yz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_z_x_0_0[i] = -2.0 * g_0_xy_0_0[i] * b_exp + 4.0 * g_zz_xy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_z_y_0_0[i] = g_0_0_0_0[i] - 2.0 * g_0_yy_0_0[i] * b_exp - 2.0 * g_zz_0_0_0[i] * a_exp + 4.0 * g_zz_yy_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_z_z_0_0[i] = -2.0 * g_0_yz_0_0[i] * b_exp + 4.0 * g_zz_yz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_xz_0_0_0, g_xz_xz_0_0, g_xz_yz_0_0, g_xz_zz_0_0, g_z_z_0_0_x_x_0_0, g_z_z_0_0_x_y_0_0, g_z_z_0_0_x_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_x_x_0_0[i] = 4.0 * g_xz_xz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_x_y_0_0[i] = 4.0 * g_xz_yz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_x_z_0_0[i] = -2.0 * g_xz_0_0_0[i] * a_exp + 4.0 * g_xz_zz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_yz_0_0_0, g_yz_xz_0_0, g_yz_yz_0_0, g_yz_zz_0_0, g_z_z_0_0_y_x_0_0, g_z_z_0_0_y_y_0_0, g_z_z_0_0_y_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_y_x_0_0[i] = 4.0 * g_yz_xz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_y_y_0_0[i] = 4.0 * g_yz_yz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_y_z_0_0[i] = -2.0 * g_yz_0_0_0[i] * a_exp + 4.0 * g_yz_zz_0_0[i] * a_exp * b_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_0_0_0_0, g_0_xz_0_0, g_0_yz_0_0, g_0_zz_0_0, g_z_z_0_0_z_x_0_0, g_z_z_0_0_z_y_0_0, g_z_z_0_0_z_z_0_0, g_zz_0_0_0, g_zz_xz_0_0, g_zz_yz_0_0, g_zz_zz_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_z_x_0_0[i] = -2.0 * g_0_xz_0_0[i] * b_exp + 4.0 * g_zz_xz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_z_y_0_0[i] = -2.0 * g_0_yz_0_0[i] * b_exp + 4.0 * g_zz_yz_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_z_z_0_0[i] = g_0_0_0_0[i] - 2.0 * g_0_zz_0_0[i] * b_exp - 2.0 * g_zz_0_0_0[i] * a_exp + 4.0 * g_zz_zz_0_0[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

