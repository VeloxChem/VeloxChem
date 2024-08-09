#include "GeomDeriv1000OfScalarForSDSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_sdsp_0(CSimdArray<double>& buffer_1000_sdsp,
                     const CSimdArray<double>& buffer_pdsp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_sdsp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pdsp

    auto g_x_xx_0_x = buffer_pdsp[0];

    auto g_x_xx_0_y = buffer_pdsp[1];

    auto g_x_xx_0_z = buffer_pdsp[2];

    auto g_x_xy_0_x = buffer_pdsp[3];

    auto g_x_xy_0_y = buffer_pdsp[4];

    auto g_x_xy_0_z = buffer_pdsp[5];

    auto g_x_xz_0_x = buffer_pdsp[6];

    auto g_x_xz_0_y = buffer_pdsp[7];

    auto g_x_xz_0_z = buffer_pdsp[8];

    auto g_x_yy_0_x = buffer_pdsp[9];

    auto g_x_yy_0_y = buffer_pdsp[10];

    auto g_x_yy_0_z = buffer_pdsp[11];

    auto g_x_yz_0_x = buffer_pdsp[12];

    auto g_x_yz_0_y = buffer_pdsp[13];

    auto g_x_yz_0_z = buffer_pdsp[14];

    auto g_x_zz_0_x = buffer_pdsp[15];

    auto g_x_zz_0_y = buffer_pdsp[16];

    auto g_x_zz_0_z = buffer_pdsp[17];

    auto g_y_xx_0_x = buffer_pdsp[18];

    auto g_y_xx_0_y = buffer_pdsp[19];

    auto g_y_xx_0_z = buffer_pdsp[20];

    auto g_y_xy_0_x = buffer_pdsp[21];

    auto g_y_xy_0_y = buffer_pdsp[22];

    auto g_y_xy_0_z = buffer_pdsp[23];

    auto g_y_xz_0_x = buffer_pdsp[24];

    auto g_y_xz_0_y = buffer_pdsp[25];

    auto g_y_xz_0_z = buffer_pdsp[26];

    auto g_y_yy_0_x = buffer_pdsp[27];

    auto g_y_yy_0_y = buffer_pdsp[28];

    auto g_y_yy_0_z = buffer_pdsp[29];

    auto g_y_yz_0_x = buffer_pdsp[30];

    auto g_y_yz_0_y = buffer_pdsp[31];

    auto g_y_yz_0_z = buffer_pdsp[32];

    auto g_y_zz_0_x = buffer_pdsp[33];

    auto g_y_zz_0_y = buffer_pdsp[34];

    auto g_y_zz_0_z = buffer_pdsp[35];

    auto g_z_xx_0_x = buffer_pdsp[36];

    auto g_z_xx_0_y = buffer_pdsp[37];

    auto g_z_xx_0_z = buffer_pdsp[38];

    auto g_z_xy_0_x = buffer_pdsp[39];

    auto g_z_xy_0_y = buffer_pdsp[40];

    auto g_z_xy_0_z = buffer_pdsp[41];

    auto g_z_xz_0_x = buffer_pdsp[42];

    auto g_z_xz_0_y = buffer_pdsp[43];

    auto g_z_xz_0_z = buffer_pdsp[44];

    auto g_z_yy_0_x = buffer_pdsp[45];

    auto g_z_yy_0_y = buffer_pdsp[46];

    auto g_z_yy_0_z = buffer_pdsp[47];

    auto g_z_yz_0_x = buffer_pdsp[48];

    auto g_z_yz_0_y = buffer_pdsp[49];

    auto g_z_yz_0_z = buffer_pdsp[50];

    auto g_z_zz_0_x = buffer_pdsp[51];

    auto g_z_zz_0_y = buffer_pdsp[52];

    auto g_z_zz_0_z = buffer_pdsp[53];

    /// Set up components of integrals buffer : buffer_1000_sdsp

    auto g_x_0_0_0_0_xx_0_x = buffer_1000_sdsp[0];

    auto g_x_0_0_0_0_xx_0_y = buffer_1000_sdsp[1];

    auto g_x_0_0_0_0_xx_0_z = buffer_1000_sdsp[2];

    auto g_x_0_0_0_0_xy_0_x = buffer_1000_sdsp[3];

    auto g_x_0_0_0_0_xy_0_y = buffer_1000_sdsp[4];

    auto g_x_0_0_0_0_xy_0_z = buffer_1000_sdsp[5];

    auto g_x_0_0_0_0_xz_0_x = buffer_1000_sdsp[6];

    auto g_x_0_0_0_0_xz_0_y = buffer_1000_sdsp[7];

    auto g_x_0_0_0_0_xz_0_z = buffer_1000_sdsp[8];

    auto g_x_0_0_0_0_yy_0_x = buffer_1000_sdsp[9];

    auto g_x_0_0_0_0_yy_0_y = buffer_1000_sdsp[10];

    auto g_x_0_0_0_0_yy_0_z = buffer_1000_sdsp[11];

    auto g_x_0_0_0_0_yz_0_x = buffer_1000_sdsp[12];

    auto g_x_0_0_0_0_yz_0_y = buffer_1000_sdsp[13];

    auto g_x_0_0_0_0_yz_0_z = buffer_1000_sdsp[14];

    auto g_x_0_0_0_0_zz_0_x = buffer_1000_sdsp[15];

    auto g_x_0_0_0_0_zz_0_y = buffer_1000_sdsp[16];

    auto g_x_0_0_0_0_zz_0_z = buffer_1000_sdsp[17];

    auto g_y_0_0_0_0_xx_0_x = buffer_1000_sdsp[18];

    auto g_y_0_0_0_0_xx_0_y = buffer_1000_sdsp[19];

    auto g_y_0_0_0_0_xx_0_z = buffer_1000_sdsp[20];

    auto g_y_0_0_0_0_xy_0_x = buffer_1000_sdsp[21];

    auto g_y_0_0_0_0_xy_0_y = buffer_1000_sdsp[22];

    auto g_y_0_0_0_0_xy_0_z = buffer_1000_sdsp[23];

    auto g_y_0_0_0_0_xz_0_x = buffer_1000_sdsp[24];

    auto g_y_0_0_0_0_xz_0_y = buffer_1000_sdsp[25];

    auto g_y_0_0_0_0_xz_0_z = buffer_1000_sdsp[26];

    auto g_y_0_0_0_0_yy_0_x = buffer_1000_sdsp[27];

    auto g_y_0_0_0_0_yy_0_y = buffer_1000_sdsp[28];

    auto g_y_0_0_0_0_yy_0_z = buffer_1000_sdsp[29];

    auto g_y_0_0_0_0_yz_0_x = buffer_1000_sdsp[30];

    auto g_y_0_0_0_0_yz_0_y = buffer_1000_sdsp[31];

    auto g_y_0_0_0_0_yz_0_z = buffer_1000_sdsp[32];

    auto g_y_0_0_0_0_zz_0_x = buffer_1000_sdsp[33];

    auto g_y_0_0_0_0_zz_0_y = buffer_1000_sdsp[34];

    auto g_y_0_0_0_0_zz_0_z = buffer_1000_sdsp[35];

    auto g_z_0_0_0_0_xx_0_x = buffer_1000_sdsp[36];

    auto g_z_0_0_0_0_xx_0_y = buffer_1000_sdsp[37];

    auto g_z_0_0_0_0_xx_0_z = buffer_1000_sdsp[38];

    auto g_z_0_0_0_0_xy_0_x = buffer_1000_sdsp[39];

    auto g_z_0_0_0_0_xy_0_y = buffer_1000_sdsp[40];

    auto g_z_0_0_0_0_xy_0_z = buffer_1000_sdsp[41];

    auto g_z_0_0_0_0_xz_0_x = buffer_1000_sdsp[42];

    auto g_z_0_0_0_0_xz_0_y = buffer_1000_sdsp[43];

    auto g_z_0_0_0_0_xz_0_z = buffer_1000_sdsp[44];

    auto g_z_0_0_0_0_yy_0_x = buffer_1000_sdsp[45];

    auto g_z_0_0_0_0_yy_0_y = buffer_1000_sdsp[46];

    auto g_z_0_0_0_0_yy_0_z = buffer_1000_sdsp[47];

    auto g_z_0_0_0_0_yz_0_x = buffer_1000_sdsp[48];

    auto g_z_0_0_0_0_yz_0_y = buffer_1000_sdsp[49];

    auto g_z_0_0_0_0_yz_0_z = buffer_1000_sdsp[50];

    auto g_z_0_0_0_0_zz_0_x = buffer_1000_sdsp[51];

    auto g_z_0_0_0_0_zz_0_y = buffer_1000_sdsp[52];

    auto g_z_0_0_0_0_zz_0_z = buffer_1000_sdsp[53];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_0_0_0_xx_0_x, g_x_0_0_0_0_xx_0_y, g_x_0_0_0_0_xx_0_z, g_x_xx_0_x, g_x_xx_0_y, g_x_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xx_0_x[i] = 2.0 * g_x_xx_0_x[i] * a_exp;

        g_x_0_0_0_0_xx_0_y[i] = 2.0 * g_x_xx_0_y[i] * a_exp;

        g_x_0_0_0_0_xx_0_z[i] = 2.0 * g_x_xx_0_z[i] * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_0_0_0_xy_0_x, g_x_0_0_0_0_xy_0_y, g_x_0_0_0_0_xy_0_z, g_x_xy_0_x, g_x_xy_0_y, g_x_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xy_0_x[i] = 2.0 * g_x_xy_0_x[i] * a_exp;

        g_x_0_0_0_0_xy_0_y[i] = 2.0 * g_x_xy_0_y[i] * a_exp;

        g_x_0_0_0_0_xy_0_z[i] = 2.0 * g_x_xy_0_z[i] * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_0_0_0_xz_0_x, g_x_0_0_0_0_xz_0_y, g_x_0_0_0_0_xz_0_z, g_x_xz_0_x, g_x_xz_0_y, g_x_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_xz_0_x[i] = 2.0 * g_x_xz_0_x[i] * a_exp;

        g_x_0_0_0_0_xz_0_y[i] = 2.0 * g_x_xz_0_y[i] * a_exp;

        g_x_0_0_0_0_xz_0_z[i] = 2.0 * g_x_xz_0_z[i] * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_0_0_0_0_yy_0_x, g_x_0_0_0_0_yy_0_y, g_x_0_0_0_0_yy_0_z, g_x_yy_0_x, g_x_yy_0_y, g_x_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yy_0_x[i] = 2.0 * g_x_yy_0_x[i] * a_exp;

        g_x_0_0_0_0_yy_0_y[i] = 2.0 * g_x_yy_0_y[i] * a_exp;

        g_x_0_0_0_0_yy_0_z[i] = 2.0 * g_x_yy_0_z[i] * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_0_0_0_0_yz_0_x, g_x_0_0_0_0_yz_0_y, g_x_0_0_0_0_yz_0_z, g_x_yz_0_x, g_x_yz_0_y, g_x_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_yz_0_x[i] = 2.0 * g_x_yz_0_x[i] * a_exp;

        g_x_0_0_0_0_yz_0_y[i] = 2.0 * g_x_yz_0_y[i] * a_exp;

        g_x_0_0_0_0_yz_0_z[i] = 2.0 * g_x_yz_0_z[i] * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_0_0_0_0_zz_0_x, g_x_0_0_0_0_zz_0_y, g_x_0_0_0_0_zz_0_z, g_x_zz_0_x, g_x_zz_0_y, g_x_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_zz_0_x[i] = 2.0 * g_x_zz_0_x[i] * a_exp;

        g_x_0_0_0_0_zz_0_y[i] = 2.0 * g_x_zz_0_y[i] * a_exp;

        g_x_0_0_0_0_zz_0_z[i] = 2.0 * g_x_zz_0_z[i] * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_y_0_0_0_0_xx_0_x, g_y_0_0_0_0_xx_0_y, g_y_0_0_0_0_xx_0_z, g_y_xx_0_x, g_y_xx_0_y, g_y_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xx_0_x[i] = 2.0 * g_y_xx_0_x[i] * a_exp;

        g_y_0_0_0_0_xx_0_y[i] = 2.0 * g_y_xx_0_y[i] * a_exp;

        g_y_0_0_0_0_xx_0_z[i] = 2.0 * g_y_xx_0_z[i] * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_y_0_0_0_0_xy_0_x, g_y_0_0_0_0_xy_0_y, g_y_0_0_0_0_xy_0_z, g_y_xy_0_x, g_y_xy_0_y, g_y_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xy_0_x[i] = 2.0 * g_y_xy_0_x[i] * a_exp;

        g_y_0_0_0_0_xy_0_y[i] = 2.0 * g_y_xy_0_y[i] * a_exp;

        g_y_0_0_0_0_xy_0_z[i] = 2.0 * g_y_xy_0_z[i] * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_y_0_0_0_0_xz_0_x, g_y_0_0_0_0_xz_0_y, g_y_0_0_0_0_xz_0_z, g_y_xz_0_x, g_y_xz_0_y, g_y_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_xz_0_x[i] = 2.0 * g_y_xz_0_x[i] * a_exp;

        g_y_0_0_0_0_xz_0_y[i] = 2.0 * g_y_xz_0_y[i] * a_exp;

        g_y_0_0_0_0_xz_0_z[i] = 2.0 * g_y_xz_0_z[i] * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_y_0_0_0_0_yy_0_x, g_y_0_0_0_0_yy_0_y, g_y_0_0_0_0_yy_0_z, g_y_yy_0_x, g_y_yy_0_y, g_y_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yy_0_x[i] = 2.0 * g_y_yy_0_x[i] * a_exp;

        g_y_0_0_0_0_yy_0_y[i] = 2.0 * g_y_yy_0_y[i] * a_exp;

        g_y_0_0_0_0_yy_0_z[i] = 2.0 * g_y_yy_0_z[i] * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_y_0_0_0_0_yz_0_x, g_y_0_0_0_0_yz_0_y, g_y_0_0_0_0_yz_0_z, g_y_yz_0_x, g_y_yz_0_y, g_y_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_yz_0_x[i] = 2.0 * g_y_yz_0_x[i] * a_exp;

        g_y_0_0_0_0_yz_0_y[i] = 2.0 * g_y_yz_0_y[i] * a_exp;

        g_y_0_0_0_0_yz_0_z[i] = 2.0 * g_y_yz_0_z[i] * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_y_0_0_0_0_zz_0_x, g_y_0_0_0_0_zz_0_y, g_y_0_0_0_0_zz_0_z, g_y_zz_0_x, g_y_zz_0_y, g_y_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_zz_0_x[i] = 2.0 * g_y_zz_0_x[i] * a_exp;

        g_y_0_0_0_0_zz_0_y[i] = 2.0 * g_y_zz_0_y[i] * a_exp;

        g_y_0_0_0_0_zz_0_z[i] = 2.0 * g_y_zz_0_z[i] * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_z_0_0_0_0_xx_0_x, g_z_0_0_0_0_xx_0_y, g_z_0_0_0_0_xx_0_z, g_z_xx_0_x, g_z_xx_0_y, g_z_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xx_0_x[i] = 2.0 * g_z_xx_0_x[i] * a_exp;

        g_z_0_0_0_0_xx_0_y[i] = 2.0 * g_z_xx_0_y[i] * a_exp;

        g_z_0_0_0_0_xx_0_z[i] = 2.0 * g_z_xx_0_z[i] * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_z_0_0_0_0_xy_0_x, g_z_0_0_0_0_xy_0_y, g_z_0_0_0_0_xy_0_z, g_z_xy_0_x, g_z_xy_0_y, g_z_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xy_0_x[i] = 2.0 * g_z_xy_0_x[i] * a_exp;

        g_z_0_0_0_0_xy_0_y[i] = 2.0 * g_z_xy_0_y[i] * a_exp;

        g_z_0_0_0_0_xy_0_z[i] = 2.0 * g_z_xy_0_z[i] * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_z_0_0_0_0_xz_0_x, g_z_0_0_0_0_xz_0_y, g_z_0_0_0_0_xz_0_z, g_z_xz_0_x, g_z_xz_0_y, g_z_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_xz_0_x[i] = 2.0 * g_z_xz_0_x[i] * a_exp;

        g_z_0_0_0_0_xz_0_y[i] = 2.0 * g_z_xz_0_y[i] * a_exp;

        g_z_0_0_0_0_xz_0_z[i] = 2.0 * g_z_xz_0_z[i] * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_z_0_0_0_0_yy_0_x, g_z_0_0_0_0_yy_0_y, g_z_0_0_0_0_yy_0_z, g_z_yy_0_x, g_z_yy_0_y, g_z_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yy_0_x[i] = 2.0 * g_z_yy_0_x[i] * a_exp;

        g_z_0_0_0_0_yy_0_y[i] = 2.0 * g_z_yy_0_y[i] * a_exp;

        g_z_0_0_0_0_yy_0_z[i] = 2.0 * g_z_yy_0_z[i] * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_z_0_0_0_0_yz_0_x, g_z_0_0_0_0_yz_0_y, g_z_0_0_0_0_yz_0_z, g_z_yz_0_x, g_z_yz_0_y, g_z_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_yz_0_x[i] = 2.0 * g_z_yz_0_x[i] * a_exp;

        g_z_0_0_0_0_yz_0_y[i] = 2.0 * g_z_yz_0_y[i] * a_exp;

        g_z_0_0_0_0_yz_0_z[i] = 2.0 * g_z_yz_0_z[i] * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_z_0_0_0_0_zz_0_x, g_z_0_0_0_0_zz_0_y, g_z_0_0_0_0_zz_0_z, g_z_zz_0_x, g_z_zz_0_y, g_z_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_zz_0_x[i] = 2.0 * g_z_zz_0_x[i] * a_exp;

        g_z_0_0_0_0_zz_0_y[i] = 2.0 * g_z_zz_0_y[i] * a_exp;

        g_z_0_0_0_0_zz_0_z[i] = 2.0 * g_z_zz_0_z[i] * a_exp;
    }
}

} // t4c_geom namespace

