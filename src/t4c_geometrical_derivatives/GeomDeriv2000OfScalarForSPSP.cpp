#include "GeomDeriv2000OfScalarForSPSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom2000_spsp_0(CSimdArray<double>& buffer_2000_spsp,
                     const CSimdArray<double>& buffer_spsp,
                     const CSimdArray<double>& buffer_dpsp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_2000_spsp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_spsp

    auto g_0_x_0_x = buffer_spsp[0];

    auto g_0_x_0_y = buffer_spsp[1];

    auto g_0_x_0_z = buffer_spsp[2];

    auto g_0_y_0_x = buffer_spsp[3];

    auto g_0_y_0_y = buffer_spsp[4];

    auto g_0_y_0_z = buffer_spsp[5];

    auto g_0_z_0_x = buffer_spsp[6];

    auto g_0_z_0_y = buffer_spsp[7];

    auto g_0_z_0_z = buffer_spsp[8];

    /// Set up components of auxilary buffer : buffer_dpsp

    auto g_xx_x_0_x = buffer_dpsp[0];

    auto g_xx_x_0_y = buffer_dpsp[1];

    auto g_xx_x_0_z = buffer_dpsp[2];

    auto g_xx_y_0_x = buffer_dpsp[3];

    auto g_xx_y_0_y = buffer_dpsp[4];

    auto g_xx_y_0_z = buffer_dpsp[5];

    auto g_xx_z_0_x = buffer_dpsp[6];

    auto g_xx_z_0_y = buffer_dpsp[7];

    auto g_xx_z_0_z = buffer_dpsp[8];

    auto g_xy_x_0_x = buffer_dpsp[9];

    auto g_xy_x_0_y = buffer_dpsp[10];

    auto g_xy_x_0_z = buffer_dpsp[11];

    auto g_xy_y_0_x = buffer_dpsp[12];

    auto g_xy_y_0_y = buffer_dpsp[13];

    auto g_xy_y_0_z = buffer_dpsp[14];

    auto g_xy_z_0_x = buffer_dpsp[15];

    auto g_xy_z_0_y = buffer_dpsp[16];

    auto g_xy_z_0_z = buffer_dpsp[17];

    auto g_xz_x_0_x = buffer_dpsp[18];

    auto g_xz_x_0_y = buffer_dpsp[19];

    auto g_xz_x_0_z = buffer_dpsp[20];

    auto g_xz_y_0_x = buffer_dpsp[21];

    auto g_xz_y_0_y = buffer_dpsp[22];

    auto g_xz_y_0_z = buffer_dpsp[23];

    auto g_xz_z_0_x = buffer_dpsp[24];

    auto g_xz_z_0_y = buffer_dpsp[25];

    auto g_xz_z_0_z = buffer_dpsp[26];

    auto g_yy_x_0_x = buffer_dpsp[27];

    auto g_yy_x_0_y = buffer_dpsp[28];

    auto g_yy_x_0_z = buffer_dpsp[29];

    auto g_yy_y_0_x = buffer_dpsp[30];

    auto g_yy_y_0_y = buffer_dpsp[31];

    auto g_yy_y_0_z = buffer_dpsp[32];

    auto g_yy_z_0_x = buffer_dpsp[33];

    auto g_yy_z_0_y = buffer_dpsp[34];

    auto g_yy_z_0_z = buffer_dpsp[35];

    auto g_yz_x_0_x = buffer_dpsp[36];

    auto g_yz_x_0_y = buffer_dpsp[37];

    auto g_yz_x_0_z = buffer_dpsp[38];

    auto g_yz_y_0_x = buffer_dpsp[39];

    auto g_yz_y_0_y = buffer_dpsp[40];

    auto g_yz_y_0_z = buffer_dpsp[41];

    auto g_yz_z_0_x = buffer_dpsp[42];

    auto g_yz_z_0_y = buffer_dpsp[43];

    auto g_yz_z_0_z = buffer_dpsp[44];

    auto g_zz_x_0_x = buffer_dpsp[45];

    auto g_zz_x_0_y = buffer_dpsp[46];

    auto g_zz_x_0_z = buffer_dpsp[47];

    auto g_zz_y_0_x = buffer_dpsp[48];

    auto g_zz_y_0_y = buffer_dpsp[49];

    auto g_zz_y_0_z = buffer_dpsp[50];

    auto g_zz_z_0_x = buffer_dpsp[51];

    auto g_zz_z_0_y = buffer_dpsp[52];

    auto g_zz_z_0_z = buffer_dpsp[53];

    /// Set up components of integrals buffer : buffer_2000_spsp

    auto g_xx_0_0_0_0_x_0_x = buffer_2000_spsp[0];

    auto g_xx_0_0_0_0_x_0_y = buffer_2000_spsp[1];

    auto g_xx_0_0_0_0_x_0_z = buffer_2000_spsp[2];

    auto g_xx_0_0_0_0_y_0_x = buffer_2000_spsp[3];

    auto g_xx_0_0_0_0_y_0_y = buffer_2000_spsp[4];

    auto g_xx_0_0_0_0_y_0_z = buffer_2000_spsp[5];

    auto g_xx_0_0_0_0_z_0_x = buffer_2000_spsp[6];

    auto g_xx_0_0_0_0_z_0_y = buffer_2000_spsp[7];

    auto g_xx_0_0_0_0_z_0_z = buffer_2000_spsp[8];

    auto g_xy_0_0_0_0_x_0_x = buffer_2000_spsp[9];

    auto g_xy_0_0_0_0_x_0_y = buffer_2000_spsp[10];

    auto g_xy_0_0_0_0_x_0_z = buffer_2000_spsp[11];

    auto g_xy_0_0_0_0_y_0_x = buffer_2000_spsp[12];

    auto g_xy_0_0_0_0_y_0_y = buffer_2000_spsp[13];

    auto g_xy_0_0_0_0_y_0_z = buffer_2000_spsp[14];

    auto g_xy_0_0_0_0_z_0_x = buffer_2000_spsp[15];

    auto g_xy_0_0_0_0_z_0_y = buffer_2000_spsp[16];

    auto g_xy_0_0_0_0_z_0_z = buffer_2000_spsp[17];

    auto g_xz_0_0_0_0_x_0_x = buffer_2000_spsp[18];

    auto g_xz_0_0_0_0_x_0_y = buffer_2000_spsp[19];

    auto g_xz_0_0_0_0_x_0_z = buffer_2000_spsp[20];

    auto g_xz_0_0_0_0_y_0_x = buffer_2000_spsp[21];

    auto g_xz_0_0_0_0_y_0_y = buffer_2000_spsp[22];

    auto g_xz_0_0_0_0_y_0_z = buffer_2000_spsp[23];

    auto g_xz_0_0_0_0_z_0_x = buffer_2000_spsp[24];

    auto g_xz_0_0_0_0_z_0_y = buffer_2000_spsp[25];

    auto g_xz_0_0_0_0_z_0_z = buffer_2000_spsp[26];

    auto g_yy_0_0_0_0_x_0_x = buffer_2000_spsp[27];

    auto g_yy_0_0_0_0_x_0_y = buffer_2000_spsp[28];

    auto g_yy_0_0_0_0_x_0_z = buffer_2000_spsp[29];

    auto g_yy_0_0_0_0_y_0_x = buffer_2000_spsp[30];

    auto g_yy_0_0_0_0_y_0_y = buffer_2000_spsp[31];

    auto g_yy_0_0_0_0_y_0_z = buffer_2000_spsp[32];

    auto g_yy_0_0_0_0_z_0_x = buffer_2000_spsp[33];

    auto g_yy_0_0_0_0_z_0_y = buffer_2000_spsp[34];

    auto g_yy_0_0_0_0_z_0_z = buffer_2000_spsp[35];

    auto g_yz_0_0_0_0_x_0_x = buffer_2000_spsp[36];

    auto g_yz_0_0_0_0_x_0_y = buffer_2000_spsp[37];

    auto g_yz_0_0_0_0_x_0_z = buffer_2000_spsp[38];

    auto g_yz_0_0_0_0_y_0_x = buffer_2000_spsp[39];

    auto g_yz_0_0_0_0_y_0_y = buffer_2000_spsp[40];

    auto g_yz_0_0_0_0_y_0_z = buffer_2000_spsp[41];

    auto g_yz_0_0_0_0_z_0_x = buffer_2000_spsp[42];

    auto g_yz_0_0_0_0_z_0_y = buffer_2000_spsp[43];

    auto g_yz_0_0_0_0_z_0_z = buffer_2000_spsp[44];

    auto g_zz_0_0_0_0_x_0_x = buffer_2000_spsp[45];

    auto g_zz_0_0_0_0_x_0_y = buffer_2000_spsp[46];

    auto g_zz_0_0_0_0_x_0_z = buffer_2000_spsp[47];

    auto g_zz_0_0_0_0_y_0_x = buffer_2000_spsp[48];

    auto g_zz_0_0_0_0_y_0_y = buffer_2000_spsp[49];

    auto g_zz_0_0_0_0_y_0_z = buffer_2000_spsp[50];

    auto g_zz_0_0_0_0_z_0_x = buffer_2000_spsp[51];

    auto g_zz_0_0_0_0_z_0_y = buffer_2000_spsp[52];

    auto g_zz_0_0_0_0_z_0_z = buffer_2000_spsp[53];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_xx_0_0_0_0_x_0_x, g_xx_0_0_0_0_x_0_y, g_xx_0_0_0_0_x_0_z, g_xx_x_0_x, g_xx_x_0_y, g_xx_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_x_0_x[i] = -2.0 * g_0_x_0_x[i] * a_exp + 4.0 * g_xx_x_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_0_y[i] = -2.0 * g_0_x_0_y[i] * a_exp + 4.0 * g_xx_x_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_x_0_z[i] = -2.0 * g_0_x_0_z[i] * a_exp + 4.0 * g_xx_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_xx_0_0_0_0_y_0_x, g_xx_0_0_0_0_y_0_y, g_xx_0_0_0_0_y_0_z, g_xx_y_0_x, g_xx_y_0_y, g_xx_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_y_0_x[i] = -2.0 * g_0_y_0_x[i] * a_exp + 4.0 * g_xx_y_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_0_y[i] = -2.0 * g_0_y_0_y[i] * a_exp + 4.0 * g_xx_y_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_y_0_z[i] = -2.0 * g_0_y_0_z[i] * a_exp + 4.0 * g_xx_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_xx_0_0_0_0_z_0_x, g_xx_0_0_0_0_z_0_y, g_xx_0_0_0_0_z_0_z, g_xx_z_0_x, g_xx_z_0_y, g_xx_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xx_0_0_0_0_z_0_x[i] = -2.0 * g_0_z_0_x[i] * a_exp + 4.0 * g_xx_z_0_x[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_0_y[i] = -2.0 * g_0_z_0_y[i] * a_exp + 4.0 * g_xx_z_0_y[i] * a_exp * a_exp;

        g_xx_0_0_0_0_z_0_z[i] = -2.0 * g_0_z_0_z[i] * a_exp + 4.0 * g_xx_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_xy_0_0_0_0_x_0_x, g_xy_0_0_0_0_x_0_y, g_xy_0_0_0_0_x_0_z, g_xy_x_0_x, g_xy_x_0_y, g_xy_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_x_0_x[i] = 4.0 * g_xy_x_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_0_y[i] = 4.0 * g_xy_x_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_x_0_z[i] = 4.0 * g_xy_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_xy_0_0_0_0_y_0_x, g_xy_0_0_0_0_y_0_y, g_xy_0_0_0_0_y_0_z, g_xy_y_0_x, g_xy_y_0_y, g_xy_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_y_0_x[i] = 4.0 * g_xy_y_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_0_y[i] = 4.0 * g_xy_y_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_y_0_z[i] = 4.0 * g_xy_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_xy_0_0_0_0_z_0_x, g_xy_0_0_0_0_z_0_y, g_xy_0_0_0_0_z_0_z, g_xy_z_0_x, g_xy_z_0_y, g_xy_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xy_0_0_0_0_z_0_x[i] = 4.0 * g_xy_z_0_x[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_0_y[i] = 4.0 * g_xy_z_0_y[i] * a_exp * a_exp;

        g_xy_0_0_0_0_z_0_z[i] = 4.0 * g_xy_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_xz_0_0_0_0_x_0_x, g_xz_0_0_0_0_x_0_y, g_xz_0_0_0_0_x_0_z, g_xz_x_0_x, g_xz_x_0_y, g_xz_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_x_0_x[i] = 4.0 * g_xz_x_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_0_y[i] = 4.0 * g_xz_x_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_x_0_z[i] = 4.0 * g_xz_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_xz_0_0_0_0_y_0_x, g_xz_0_0_0_0_y_0_y, g_xz_0_0_0_0_y_0_z, g_xz_y_0_x, g_xz_y_0_y, g_xz_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_y_0_x[i] = 4.0 * g_xz_y_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_0_y[i] = 4.0 * g_xz_y_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_y_0_z[i] = 4.0 * g_xz_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_xz_0_0_0_0_z_0_x, g_xz_0_0_0_0_z_0_y, g_xz_0_0_0_0_z_0_z, g_xz_z_0_x, g_xz_z_0_y, g_xz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_xz_0_0_0_0_z_0_x[i] = 4.0 * g_xz_z_0_x[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_0_y[i] = 4.0 * g_xz_z_0_y[i] * a_exp * a_exp;

        g_xz_0_0_0_0_z_0_z[i] = 4.0 * g_xz_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_yy_0_0_0_0_x_0_x, g_yy_0_0_0_0_x_0_y, g_yy_0_0_0_0_x_0_z, g_yy_x_0_x, g_yy_x_0_y, g_yy_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_x_0_x[i] = -2.0 * g_0_x_0_x[i] * a_exp + 4.0 * g_yy_x_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_0_y[i] = -2.0 * g_0_x_0_y[i] * a_exp + 4.0 * g_yy_x_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_x_0_z[i] = -2.0 * g_0_x_0_z[i] * a_exp + 4.0 * g_yy_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_yy_0_0_0_0_y_0_x, g_yy_0_0_0_0_y_0_y, g_yy_0_0_0_0_y_0_z, g_yy_y_0_x, g_yy_y_0_y, g_yy_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_y_0_x[i] = -2.0 * g_0_y_0_x[i] * a_exp + 4.0 * g_yy_y_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_0_y[i] = -2.0 * g_0_y_0_y[i] * a_exp + 4.0 * g_yy_y_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_y_0_z[i] = -2.0 * g_0_y_0_z[i] * a_exp + 4.0 * g_yy_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_yy_0_0_0_0_z_0_x, g_yy_0_0_0_0_z_0_y, g_yy_0_0_0_0_z_0_z, g_yy_z_0_x, g_yy_z_0_y, g_yy_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yy_0_0_0_0_z_0_x[i] = -2.0 * g_0_z_0_x[i] * a_exp + 4.0 * g_yy_z_0_x[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_0_y[i] = -2.0 * g_0_z_0_y[i] * a_exp + 4.0 * g_yy_z_0_y[i] * a_exp * a_exp;

        g_yy_0_0_0_0_z_0_z[i] = -2.0 * g_0_z_0_z[i] * a_exp + 4.0 * g_yy_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_yz_0_0_0_0_x_0_x, g_yz_0_0_0_0_x_0_y, g_yz_0_0_0_0_x_0_z, g_yz_x_0_x, g_yz_x_0_y, g_yz_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_x_0_x[i] = 4.0 * g_yz_x_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_0_y[i] = 4.0 * g_yz_x_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_x_0_z[i] = 4.0 * g_yz_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_yz_0_0_0_0_y_0_x, g_yz_0_0_0_0_y_0_y, g_yz_0_0_0_0_y_0_z, g_yz_y_0_x, g_yz_y_0_y, g_yz_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_y_0_x[i] = 4.0 * g_yz_y_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_0_y[i] = 4.0 * g_yz_y_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_y_0_z[i] = 4.0 * g_yz_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_yz_0_0_0_0_z_0_x, g_yz_0_0_0_0_z_0_y, g_yz_0_0_0_0_z_0_z, g_yz_z_0_x, g_yz_z_0_y, g_yz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_yz_0_0_0_0_z_0_x[i] = 4.0 * g_yz_z_0_x[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_0_y[i] = 4.0 * g_yz_z_0_y[i] * a_exp * a_exp;

        g_yz_0_0_0_0_z_0_z[i] = 4.0 * g_yz_z_0_z[i] * a_exp * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_zz_0_0_0_0_x_0_x, g_zz_0_0_0_0_x_0_y, g_zz_0_0_0_0_x_0_z, g_zz_x_0_x, g_zz_x_0_y, g_zz_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_x_0_x[i] = -2.0 * g_0_x_0_x[i] * a_exp + 4.0 * g_zz_x_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_0_y[i] = -2.0 * g_0_x_0_y[i] * a_exp + 4.0 * g_zz_x_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_x_0_z[i] = -2.0 * g_0_x_0_z[i] * a_exp + 4.0 * g_zz_x_0_z[i] * a_exp * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_zz_0_0_0_0_y_0_x, g_zz_0_0_0_0_y_0_y, g_zz_0_0_0_0_y_0_z, g_zz_y_0_x, g_zz_y_0_y, g_zz_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_y_0_x[i] = -2.0 * g_0_y_0_x[i] * a_exp + 4.0 * g_zz_y_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_0_y[i] = -2.0 * g_0_y_0_y[i] * a_exp + 4.0 * g_zz_y_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_y_0_z[i] = -2.0 * g_0_y_0_z[i] * a_exp + 4.0 * g_zz_y_0_z[i] * a_exp * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_zz_0_0_0_0_z_0_x, g_zz_0_0_0_0_z_0_y, g_zz_0_0_0_0_z_0_z, g_zz_z_0_x, g_zz_z_0_y, g_zz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_zz_0_0_0_0_z_0_x[i] = -2.0 * g_0_z_0_x[i] * a_exp + 4.0 * g_zz_z_0_x[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_0_y[i] = -2.0 * g_0_z_0_y[i] * a_exp + 4.0 * g_zz_z_0_y[i] * a_exp * a_exp;

        g_zz_0_0_0_0_z_0_z[i] = -2.0 * g_0_z_0_z[i] * a_exp + 4.0 * g_zz_z_0_z[i] * a_exp * a_exp;
    }
}

} // t4c_geom namespace

