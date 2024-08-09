#include "GeomDeriv2000OfScalarForSSPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_sspp_0(CSimdArray<double>& buffer_2000_sspp,
                     const CSimdArray<double>& buffer_sspp,
                     const CSimdArray<double>& buffer_dspp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_sspp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_sspp

    auto g_0_0_x_x = buffer_sspp[0];

    auto g_0_0_x_y = buffer_sspp[1];

    auto g_0_0_x_z = buffer_sspp[2];

    auto g_0_0_y_x = buffer_sspp[3];

    auto g_0_0_y_y = buffer_sspp[4];

    auto g_0_0_y_z = buffer_sspp[5];

    auto g_0_0_z_x = buffer_sspp[6];

    auto g_0_0_z_y = buffer_sspp[7];

    auto g_0_0_z_z = buffer_sspp[8];

    /// Set up components of auxilary buffer : buffer_dspp

    auto g_xx_0_x_x = buffer_dspp[0];

    auto g_xx_0_x_y = buffer_dspp[1];

    auto g_xx_0_x_z = buffer_dspp[2];

    auto g_xx_0_y_x = buffer_dspp[3];

    auto g_xx_0_y_y = buffer_dspp[4];

    auto g_xx_0_y_z = buffer_dspp[5];

    auto g_xx_0_z_x = buffer_dspp[6];

    auto g_xx_0_z_y = buffer_dspp[7];

    auto g_xx_0_z_z = buffer_dspp[8];

    auto g_xy_0_x_x = buffer_dspp[9];

    auto g_xy_0_x_y = buffer_dspp[10];

    auto g_xy_0_x_z = buffer_dspp[11];

    auto g_xy_0_y_x = buffer_dspp[12];

    auto g_xy_0_y_y = buffer_dspp[13];

    auto g_xy_0_y_z = buffer_dspp[14];

    auto g_xy_0_z_x = buffer_dspp[15];

    auto g_xy_0_z_y = buffer_dspp[16];

    auto g_xy_0_z_z = buffer_dspp[17];

    auto g_xz_0_x_x = buffer_dspp[18];

    auto g_xz_0_x_y = buffer_dspp[19];

    auto g_xz_0_x_z = buffer_dspp[20];

    auto g_xz_0_y_x = buffer_dspp[21];

    auto g_xz_0_y_y = buffer_dspp[22];

    auto g_xz_0_y_z = buffer_dspp[23];

    auto g_xz_0_z_x = buffer_dspp[24];

    auto g_xz_0_z_y = buffer_dspp[25];

    auto g_xz_0_z_z = buffer_dspp[26];

    auto g_yy_0_x_x = buffer_dspp[27];

    auto g_yy_0_x_y = buffer_dspp[28];

    auto g_yy_0_x_z = buffer_dspp[29];

    auto g_yy_0_y_x = buffer_dspp[30];

    auto g_yy_0_y_y = buffer_dspp[31];

    auto g_yy_0_y_z = buffer_dspp[32];

    auto g_yy_0_z_x = buffer_dspp[33];

    auto g_yy_0_z_y = buffer_dspp[34];

    auto g_yy_0_z_z = buffer_dspp[35];

    auto g_yz_0_x_x = buffer_dspp[36];

    auto g_yz_0_x_y = buffer_dspp[37];

    auto g_yz_0_x_z = buffer_dspp[38];

    auto g_yz_0_y_x = buffer_dspp[39];

    auto g_yz_0_y_y = buffer_dspp[40];

    auto g_yz_0_y_z = buffer_dspp[41];

    auto g_yz_0_z_x = buffer_dspp[42];

    auto g_yz_0_z_y = buffer_dspp[43];

    auto g_yz_0_z_z = buffer_dspp[44];

    auto g_zz_0_x_x = buffer_dspp[45];

    auto g_zz_0_x_y = buffer_dspp[46];

    auto g_zz_0_x_z = buffer_dspp[47];

    auto g_zz_0_y_x = buffer_dspp[48];

    auto g_zz_0_y_y = buffer_dspp[49];

    auto g_zz_0_y_z = buffer_dspp[50];

    auto g_zz_0_z_x = buffer_dspp[51];

    auto g_zz_0_z_y = buffer_dspp[52];

    auto g_zz_0_z_z = buffer_dspp[53];

    /// Set up components of integrals buffer : buffer_2000_sspp

    auto g_xx_0_0_0_0_0_x_x = buffer_2000_sspp[0];

    auto g_xx_0_0_0_0_0_x_y = buffer_2000_sspp[1];

    auto g_xx_0_0_0_0_0_x_z = buffer_2000_sspp[2];

    auto g_xx_0_0_0_0_0_y_x = buffer_2000_sspp[3];

    auto g_xx_0_0_0_0_0_y_y = buffer_2000_sspp[4];

    auto g_xx_0_0_0_0_0_y_z = buffer_2000_sspp[5];

    auto g_xx_0_0_0_0_0_z_x = buffer_2000_sspp[6];

    auto g_xx_0_0_0_0_0_z_y = buffer_2000_sspp[7];

    auto g_xx_0_0_0_0_0_z_z = buffer_2000_sspp[8];

    auto g_xy_0_0_0_0_0_x_x = buffer_2000_sspp[9];

    auto g_xy_0_0_0_0_0_x_y = buffer_2000_sspp[10];

    auto g_xy_0_0_0_0_0_x_z = buffer_2000_sspp[11];

    auto g_xy_0_0_0_0_0_y_x = buffer_2000_sspp[12];

    auto g_xy_0_0_0_0_0_y_y = buffer_2000_sspp[13];

    auto g_xy_0_0_0_0_0_y_z = buffer_2000_sspp[14];

    auto g_xy_0_0_0_0_0_z_x = buffer_2000_sspp[15];

    auto g_xy_0_0_0_0_0_z_y = buffer_2000_sspp[16];

    auto g_xy_0_0_0_0_0_z_z = buffer_2000_sspp[17];

    auto g_xz_0_0_0_0_0_x_x = buffer_2000_sspp[18];

    auto g_xz_0_0_0_0_0_x_y = buffer_2000_sspp[19];

    auto g_xz_0_0_0_0_0_x_z = buffer_2000_sspp[20];

    auto g_xz_0_0_0_0_0_y_x = buffer_2000_sspp[21];

    auto g_xz_0_0_0_0_0_y_y = buffer_2000_sspp[22];

    auto g_xz_0_0_0_0_0_y_z = buffer_2000_sspp[23];

    auto g_xz_0_0_0_0_0_z_x = buffer_2000_sspp[24];

    auto g_xz_0_0_0_0_0_z_y = buffer_2000_sspp[25];

    auto g_xz_0_0_0_0_0_z_z = buffer_2000_sspp[26];

    auto g_yy_0_0_0_0_0_x_x = buffer_2000_sspp[27];

    auto g_yy_0_0_0_0_0_x_y = buffer_2000_sspp[28];

    auto g_yy_0_0_0_0_0_x_z = buffer_2000_sspp[29];

    auto g_yy_0_0_0_0_0_y_x = buffer_2000_sspp[30];

    auto g_yy_0_0_0_0_0_y_y = buffer_2000_sspp[31];

    auto g_yy_0_0_0_0_0_y_z = buffer_2000_sspp[32];

    auto g_yy_0_0_0_0_0_z_x = buffer_2000_sspp[33];

    auto g_yy_0_0_0_0_0_z_y = buffer_2000_sspp[34];

    auto g_yy_0_0_0_0_0_z_z = buffer_2000_sspp[35];

    auto g_yz_0_0_0_0_0_x_x = buffer_2000_sspp[36];

    auto g_yz_0_0_0_0_0_x_y = buffer_2000_sspp[37];

    auto g_yz_0_0_0_0_0_x_z = buffer_2000_sspp[38];

    auto g_yz_0_0_0_0_0_y_x = buffer_2000_sspp[39];

    auto g_yz_0_0_0_0_0_y_y = buffer_2000_sspp[40];

    auto g_yz_0_0_0_0_0_y_z = buffer_2000_sspp[41];

    auto g_yz_0_0_0_0_0_z_x = buffer_2000_sspp[42];

    auto g_yz_0_0_0_0_0_z_y = buffer_2000_sspp[43];

    auto g_yz_0_0_0_0_0_z_z = buffer_2000_sspp[44];

    auto g_zz_0_0_0_0_0_x_x = buffer_2000_sspp[45];

    auto g_zz_0_0_0_0_0_x_y = buffer_2000_sspp[46];

    auto g_zz_0_0_0_0_0_x_z = buffer_2000_sspp[47];

    auto g_zz_0_0_0_0_0_y_x = buffer_2000_sspp[48];

    auto g_zz_0_0_0_0_0_y_y = buffer_2000_sspp[49];

    auto g_zz_0_0_0_0_0_y_z = buffer_2000_sspp[50];

    auto g_zz_0_0_0_0_0_z_x = buffer_2000_sspp[51];

    auto g_zz_0_0_0_0_0_z_y = buffer_2000_sspp[52];

    auto g_zz_0_0_0_0_0_z_z = buffer_2000_sspp[53];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_0_x_x, g_0_0_x_y, g_0_0_x_z, g_xx_0_0_0_0_0_x_x, g_xx_0_0_0_0_0_x_y, g_xx_0_0_0_0_0_x_z, g_xx_0_x_x, g_xx_0_x_y, g_xx_0_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_0_x_x[i] = -2.0 * g_0_0_x_x[i] * a_exp + 4.0 * g_xx_0_x_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_x_y[i] = -2.0 * g_0_0_x_y[i] * a_exp + 4.0 * g_xx_0_x_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_x_z[i] = -2.0 * g_0_0_x_z[i] * a_exp + 4.0 * g_xx_0_x_z[i] * a_exp * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_0_0_y_x, g_0_0_y_y, g_0_0_y_z, g_xx_0_0_0_0_0_y_x, g_xx_0_0_0_0_0_y_y, g_xx_0_0_0_0_0_y_z, g_xx_0_y_x, g_xx_0_y_y, g_xx_0_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_0_y_x[i] = -2.0 * g_0_0_y_x[i] * a_exp + 4.0 * g_xx_0_y_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_y_y[i] = -2.0 * g_0_0_y_y[i] * a_exp + 4.0 * g_xx_0_y_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_y_z[i] = -2.0 * g_0_0_y_z[i] * a_exp + 4.0 * g_xx_0_y_z[i] * a_exp * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_0_0_z_x, g_0_0_z_y, g_0_0_z_z, g_xx_0_0_0_0_0_z_x, g_xx_0_0_0_0_0_z_y, g_xx_0_0_0_0_0_z_z, g_xx_0_z_x, g_xx_0_z_y, g_xx_0_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_0_z_x[i] = -2.0 * g_0_0_z_x[i] * a_exp + 4.0 * g_xx_0_z_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_z_y[i] = -2.0 * g_0_0_z_y[i] * a_exp + 4.0 * g_xx_0_z_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_0_z_z[i] = -2.0 * g_0_0_z_z[i] * a_exp + 4.0 * g_xx_0_z_z[i] * a_exp * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_xy_0_0_0_0_0_x_x, g_xy_0_0_0_0_0_x_y, g_xy_0_0_0_0_0_x_z, g_xy_0_x_x, g_xy_0_x_y, g_xy_0_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_0_x_x[i] = 4.0 * g_xy_0_x_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_x_y[i] = 4.0 * g_xy_0_x_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_x_z[i] = 4.0 * g_xy_0_x_z[i] * a_exp * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_xy_0_0_0_0_0_y_x, g_xy_0_0_0_0_0_y_y, g_xy_0_0_0_0_0_y_z, g_xy_0_y_x, g_xy_0_y_y, g_xy_0_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_0_y_x[i] = 4.0 * g_xy_0_y_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_y_y[i] = 4.0 * g_xy_0_y_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_y_z[i] = 4.0 * g_xy_0_y_z[i] * a_exp * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_xy_0_0_0_0_0_z_x, g_xy_0_0_0_0_0_z_y, g_xy_0_0_0_0_0_z_z, g_xy_0_z_x, g_xy_0_z_y, g_xy_0_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_0_z_x[i] = 4.0 * g_xy_0_z_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_z_y[i] = 4.0 * g_xy_0_z_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_0_z_z[i] = 4.0 * g_xy_0_z_z[i] * a_exp * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_xz_0_0_0_0_0_x_x, g_xz_0_0_0_0_0_x_y, g_xz_0_0_0_0_0_x_z, g_xz_0_x_x, g_xz_0_x_y, g_xz_0_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_0_x_x[i] = 4.0 * g_xz_0_x_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_x_y[i] = 4.0 * g_xz_0_x_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_x_z[i] = 4.0 * g_xz_0_x_z[i] * a_exp * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_xz_0_0_0_0_0_y_x, g_xz_0_0_0_0_0_y_y, g_xz_0_0_0_0_0_y_z, g_xz_0_y_x, g_xz_0_y_y, g_xz_0_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_0_y_x[i] = 4.0 * g_xz_0_y_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_y_y[i] = 4.0 * g_xz_0_y_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_y_z[i] = 4.0 * g_xz_0_y_z[i] * a_exp * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_xz_0_0_0_0_0_z_x, g_xz_0_0_0_0_0_z_y, g_xz_0_0_0_0_0_z_z, g_xz_0_z_x, g_xz_0_z_y, g_xz_0_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_0_z_x[i] = 4.0 * g_xz_0_z_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_z_y[i] = 4.0 * g_xz_0_z_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_0_z_z[i] = 4.0 * g_xz_0_z_z[i] * a_exp * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_0_0_x_x, g_0_0_x_y, g_0_0_x_z, g_yy_0_0_0_0_0_x_x, g_yy_0_0_0_0_0_x_y, g_yy_0_0_0_0_0_x_z, g_yy_0_x_x, g_yy_0_x_y, g_yy_0_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_0_x_x[i] = -2.0 * g_0_0_x_x[i] * a_exp + 4.0 * g_yy_0_x_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_x_y[i] = -2.0 * g_0_0_x_y[i] * a_exp + 4.0 * g_yy_0_x_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_x_z[i] = -2.0 * g_0_0_x_z[i] * a_exp + 4.0 * g_yy_0_x_z[i] * a_exp * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_0_0_y_x, g_0_0_y_y, g_0_0_y_z, g_yy_0_0_0_0_0_y_x, g_yy_0_0_0_0_0_y_y, g_yy_0_0_0_0_0_y_z, g_yy_0_y_x, g_yy_0_y_y, g_yy_0_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_0_y_x[i] = -2.0 * g_0_0_y_x[i] * a_exp + 4.0 * g_yy_0_y_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_y_y[i] = -2.0 * g_0_0_y_y[i] * a_exp + 4.0 * g_yy_0_y_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_y_z[i] = -2.0 * g_0_0_y_z[i] * a_exp + 4.0 * g_yy_0_y_z[i] * a_exp * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_0_0_z_x, g_0_0_z_y, g_0_0_z_z, g_yy_0_0_0_0_0_z_x, g_yy_0_0_0_0_0_z_y, g_yy_0_0_0_0_0_z_z, g_yy_0_z_x, g_yy_0_z_y, g_yy_0_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_0_z_x[i] = -2.0 * g_0_0_z_x[i] * a_exp + 4.0 * g_yy_0_z_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_z_y[i] = -2.0 * g_0_0_z_y[i] * a_exp + 4.0 * g_yy_0_z_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_0_z_z[i] = -2.0 * g_0_0_z_z[i] * a_exp + 4.0 * g_yy_0_z_z[i] * a_exp * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_yz_0_0_0_0_0_x_x, g_yz_0_0_0_0_0_x_y, g_yz_0_0_0_0_0_x_z, g_yz_0_x_x, g_yz_0_x_y, g_yz_0_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_0_x_x[i] = 4.0 * g_yz_0_x_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_x_y[i] = 4.0 * g_yz_0_x_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_x_z[i] = 4.0 * g_yz_0_x_z[i] * a_exp * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_yz_0_0_0_0_0_y_x, g_yz_0_0_0_0_0_y_y, g_yz_0_0_0_0_0_y_z, g_yz_0_y_x, g_yz_0_y_y, g_yz_0_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_0_y_x[i] = 4.0 * g_yz_0_y_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_y_y[i] = 4.0 * g_yz_0_y_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_y_z[i] = 4.0 * g_yz_0_y_z[i] * a_exp * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_yz_0_0_0_0_0_z_x, g_yz_0_0_0_0_0_z_y, g_yz_0_0_0_0_0_z_z, g_yz_0_z_x, g_yz_0_z_y, g_yz_0_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_0_z_x[i] = 4.0 * g_yz_0_z_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_z_y[i] = 4.0 * g_yz_0_z_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_0_z_z[i] = 4.0 * g_yz_0_z_z[i] * a_exp * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_0_0_x_x, g_0_0_x_y, g_0_0_x_z, g_zz_0_0_0_0_0_x_x, g_zz_0_0_0_0_0_x_y, g_zz_0_0_0_0_0_x_z, g_zz_0_x_x, g_zz_0_x_y, g_zz_0_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_0_x_x[i] = -2.0 * g_0_0_x_x[i] * a_exp + 4.0 * g_zz_0_x_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_x_y[i] = -2.0 * g_0_0_x_y[i] * a_exp + 4.0 * g_zz_0_x_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_x_z[i] = -2.0 * g_0_0_x_z[i] * a_exp + 4.0 * g_zz_0_x_z[i] * a_exp * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_0_0_y_x, g_0_0_y_y, g_0_0_y_z, g_zz_0_0_0_0_0_y_x, g_zz_0_0_0_0_0_y_y, g_zz_0_0_0_0_0_y_z, g_zz_0_y_x, g_zz_0_y_y, g_zz_0_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_0_y_x[i] = -2.0 * g_0_0_y_x[i] * a_exp + 4.0 * g_zz_0_y_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_y_y[i] = -2.0 * g_0_0_y_y[i] * a_exp + 4.0 * g_zz_0_y_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_y_z[i] = -2.0 * g_0_0_y_z[i] * a_exp + 4.0 * g_zz_0_y_z[i] * a_exp * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_0_0_z_x, g_0_0_z_y, g_0_0_z_z, g_zz_0_0_0_0_0_z_x, g_zz_0_0_0_0_0_z_y, g_zz_0_0_0_0_0_z_z, g_zz_0_z_x, g_zz_0_z_y, g_zz_0_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_0_z_x[i] = -2.0 * g_0_0_z_x[i] * a_exp + 4.0 * g_zz_0_z_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_z_y[i] = -2.0 * g_0_0_z_y[i] * a_exp + 4.0 * g_zz_0_z_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_0_z_z[i] = -2.0 * g_0_0_z_z[i] * a_exp + 4.0 * g_zz_0_z_z[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

