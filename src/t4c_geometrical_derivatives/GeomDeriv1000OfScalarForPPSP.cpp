#include "GeomDeriv1000OfScalarForPPSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_ppsp_0(CSimdArray<double>& buffer_1000_ppsp,
                     const CSimdArray<double>& buffer_spsp,
                     const CSimdArray<double>& buffer_dpsp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_ppsp.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_ppsp

    auto g_x_0_0_0_x_x_0_x = buffer_1000_ppsp[0];

    auto g_x_0_0_0_x_x_0_y = buffer_1000_ppsp[1];

    auto g_x_0_0_0_x_x_0_z = buffer_1000_ppsp[2];

    auto g_x_0_0_0_x_y_0_x = buffer_1000_ppsp[3];

    auto g_x_0_0_0_x_y_0_y = buffer_1000_ppsp[4];

    auto g_x_0_0_0_x_y_0_z = buffer_1000_ppsp[5];

    auto g_x_0_0_0_x_z_0_x = buffer_1000_ppsp[6];

    auto g_x_0_0_0_x_z_0_y = buffer_1000_ppsp[7];

    auto g_x_0_0_0_x_z_0_z = buffer_1000_ppsp[8];

    auto g_x_0_0_0_y_x_0_x = buffer_1000_ppsp[9];

    auto g_x_0_0_0_y_x_0_y = buffer_1000_ppsp[10];

    auto g_x_0_0_0_y_x_0_z = buffer_1000_ppsp[11];

    auto g_x_0_0_0_y_y_0_x = buffer_1000_ppsp[12];

    auto g_x_0_0_0_y_y_0_y = buffer_1000_ppsp[13];

    auto g_x_0_0_0_y_y_0_z = buffer_1000_ppsp[14];

    auto g_x_0_0_0_y_z_0_x = buffer_1000_ppsp[15];

    auto g_x_0_0_0_y_z_0_y = buffer_1000_ppsp[16];

    auto g_x_0_0_0_y_z_0_z = buffer_1000_ppsp[17];

    auto g_x_0_0_0_z_x_0_x = buffer_1000_ppsp[18];

    auto g_x_0_0_0_z_x_0_y = buffer_1000_ppsp[19];

    auto g_x_0_0_0_z_x_0_z = buffer_1000_ppsp[20];

    auto g_x_0_0_0_z_y_0_x = buffer_1000_ppsp[21];

    auto g_x_0_0_0_z_y_0_y = buffer_1000_ppsp[22];

    auto g_x_0_0_0_z_y_0_z = buffer_1000_ppsp[23];

    auto g_x_0_0_0_z_z_0_x = buffer_1000_ppsp[24];

    auto g_x_0_0_0_z_z_0_y = buffer_1000_ppsp[25];

    auto g_x_0_0_0_z_z_0_z = buffer_1000_ppsp[26];

    auto g_y_0_0_0_x_x_0_x = buffer_1000_ppsp[27];

    auto g_y_0_0_0_x_x_0_y = buffer_1000_ppsp[28];

    auto g_y_0_0_0_x_x_0_z = buffer_1000_ppsp[29];

    auto g_y_0_0_0_x_y_0_x = buffer_1000_ppsp[30];

    auto g_y_0_0_0_x_y_0_y = buffer_1000_ppsp[31];

    auto g_y_0_0_0_x_y_0_z = buffer_1000_ppsp[32];

    auto g_y_0_0_0_x_z_0_x = buffer_1000_ppsp[33];

    auto g_y_0_0_0_x_z_0_y = buffer_1000_ppsp[34];

    auto g_y_0_0_0_x_z_0_z = buffer_1000_ppsp[35];

    auto g_y_0_0_0_y_x_0_x = buffer_1000_ppsp[36];

    auto g_y_0_0_0_y_x_0_y = buffer_1000_ppsp[37];

    auto g_y_0_0_0_y_x_0_z = buffer_1000_ppsp[38];

    auto g_y_0_0_0_y_y_0_x = buffer_1000_ppsp[39];

    auto g_y_0_0_0_y_y_0_y = buffer_1000_ppsp[40];

    auto g_y_0_0_0_y_y_0_z = buffer_1000_ppsp[41];

    auto g_y_0_0_0_y_z_0_x = buffer_1000_ppsp[42];

    auto g_y_0_0_0_y_z_0_y = buffer_1000_ppsp[43];

    auto g_y_0_0_0_y_z_0_z = buffer_1000_ppsp[44];

    auto g_y_0_0_0_z_x_0_x = buffer_1000_ppsp[45];

    auto g_y_0_0_0_z_x_0_y = buffer_1000_ppsp[46];

    auto g_y_0_0_0_z_x_0_z = buffer_1000_ppsp[47];

    auto g_y_0_0_0_z_y_0_x = buffer_1000_ppsp[48];

    auto g_y_0_0_0_z_y_0_y = buffer_1000_ppsp[49];

    auto g_y_0_0_0_z_y_0_z = buffer_1000_ppsp[50];

    auto g_y_0_0_0_z_z_0_x = buffer_1000_ppsp[51];

    auto g_y_0_0_0_z_z_0_y = buffer_1000_ppsp[52];

    auto g_y_0_0_0_z_z_0_z = buffer_1000_ppsp[53];

    auto g_z_0_0_0_x_x_0_x = buffer_1000_ppsp[54];

    auto g_z_0_0_0_x_x_0_y = buffer_1000_ppsp[55];

    auto g_z_0_0_0_x_x_0_z = buffer_1000_ppsp[56];

    auto g_z_0_0_0_x_y_0_x = buffer_1000_ppsp[57];

    auto g_z_0_0_0_x_y_0_y = buffer_1000_ppsp[58];

    auto g_z_0_0_0_x_y_0_z = buffer_1000_ppsp[59];

    auto g_z_0_0_0_x_z_0_x = buffer_1000_ppsp[60];

    auto g_z_0_0_0_x_z_0_y = buffer_1000_ppsp[61];

    auto g_z_0_0_0_x_z_0_z = buffer_1000_ppsp[62];

    auto g_z_0_0_0_y_x_0_x = buffer_1000_ppsp[63];

    auto g_z_0_0_0_y_x_0_y = buffer_1000_ppsp[64];

    auto g_z_0_0_0_y_x_0_z = buffer_1000_ppsp[65];

    auto g_z_0_0_0_y_y_0_x = buffer_1000_ppsp[66];

    auto g_z_0_0_0_y_y_0_y = buffer_1000_ppsp[67];

    auto g_z_0_0_0_y_y_0_z = buffer_1000_ppsp[68];

    auto g_z_0_0_0_y_z_0_x = buffer_1000_ppsp[69];

    auto g_z_0_0_0_y_z_0_y = buffer_1000_ppsp[70];

    auto g_z_0_0_0_y_z_0_z = buffer_1000_ppsp[71];

    auto g_z_0_0_0_z_x_0_x = buffer_1000_ppsp[72];

    auto g_z_0_0_0_z_x_0_y = buffer_1000_ppsp[73];

    auto g_z_0_0_0_z_x_0_z = buffer_1000_ppsp[74];

    auto g_z_0_0_0_z_y_0_x = buffer_1000_ppsp[75];

    auto g_z_0_0_0_z_y_0_y = buffer_1000_ppsp[76];

    auto g_z_0_0_0_z_y_0_z = buffer_1000_ppsp[77];

    auto g_z_0_0_0_z_z_0_x = buffer_1000_ppsp[78];

    auto g_z_0_0_0_z_z_0_y = buffer_1000_ppsp[79];

    auto g_z_0_0_0_z_z_0_z = buffer_1000_ppsp[80];

    // integrals block (0-3)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_x_0_0_0_x_x_0_x, g_x_0_0_0_x_x_0_y, g_x_0_0_0_x_x_0_z, g_xx_x_0_x, g_xx_x_0_y, g_xx_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_x_0_x[i] = -g_0_x_0_x[i] + 2.0 * g_xx_x_0_x[i] * a_exp;

        g_x_0_0_0_x_x_0_y[i] = -g_0_x_0_y[i] + 2.0 * g_xx_x_0_y[i] * a_exp;

        g_x_0_0_0_x_x_0_z[i] = -g_0_x_0_z[i] + 2.0 * g_xx_x_0_z[i] * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_x_0_0_0_x_y_0_x, g_x_0_0_0_x_y_0_y, g_x_0_0_0_x_y_0_z, g_xx_y_0_x, g_xx_y_0_y, g_xx_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_y_0_x[i] = -g_0_y_0_x[i] + 2.0 * g_xx_y_0_x[i] * a_exp;

        g_x_0_0_0_x_y_0_y[i] = -g_0_y_0_y[i] + 2.0 * g_xx_y_0_y[i] * a_exp;

        g_x_0_0_0_x_y_0_z[i] = -g_0_y_0_z[i] + 2.0 * g_xx_y_0_z[i] * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_x_0_0_0_x_z_0_x, g_x_0_0_0_x_z_0_y, g_x_0_0_0_x_z_0_z, g_xx_z_0_x, g_xx_z_0_y, g_xx_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_x_z_0_x[i] = -g_0_z_0_x[i] + 2.0 * g_xx_z_0_x[i] * a_exp;

        g_x_0_0_0_x_z_0_y[i] = -g_0_z_0_y[i] + 2.0 * g_xx_z_0_y[i] * a_exp;

        g_x_0_0_0_x_z_0_z[i] = -g_0_z_0_z[i] + 2.0 * g_xx_z_0_z[i] * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_0_0_0_y_x_0_x, g_x_0_0_0_y_x_0_y, g_x_0_0_0_y_x_0_z, g_xy_x_0_x, g_xy_x_0_y, g_xy_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_x_0_x[i] = 2.0 * g_xy_x_0_x[i] * a_exp;

        g_x_0_0_0_y_x_0_y[i] = 2.0 * g_xy_x_0_y[i] * a_exp;

        g_x_0_0_0_y_x_0_z[i] = 2.0 * g_xy_x_0_z[i] * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_0_0_0_y_y_0_x, g_x_0_0_0_y_y_0_y, g_x_0_0_0_y_y_0_z, g_xy_y_0_x, g_xy_y_0_y, g_xy_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_y_0_x[i] = 2.0 * g_xy_y_0_x[i] * a_exp;

        g_x_0_0_0_y_y_0_y[i] = 2.0 * g_xy_y_0_y[i] * a_exp;

        g_x_0_0_0_y_y_0_z[i] = 2.0 * g_xy_y_0_z[i] * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_0_0_0_y_z_0_x, g_x_0_0_0_y_z_0_y, g_x_0_0_0_y_z_0_z, g_xy_z_0_x, g_xy_z_0_y, g_xy_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_y_z_0_x[i] = 2.0 * g_xy_z_0_x[i] * a_exp;

        g_x_0_0_0_y_z_0_y[i] = 2.0 * g_xy_z_0_y[i] * a_exp;

        g_x_0_0_0_y_z_0_z[i] = 2.0 * g_xy_z_0_z[i] * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_0_0_0_z_x_0_x, g_x_0_0_0_z_x_0_y, g_x_0_0_0_z_x_0_z, g_xz_x_0_x, g_xz_x_0_y, g_xz_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_x_0_x[i] = 2.0 * g_xz_x_0_x[i] * a_exp;

        g_x_0_0_0_z_x_0_y[i] = 2.0 * g_xz_x_0_y[i] * a_exp;

        g_x_0_0_0_z_x_0_z[i] = 2.0 * g_xz_x_0_z[i] * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_0_0_0_z_y_0_x, g_x_0_0_0_z_y_0_y, g_x_0_0_0_z_y_0_z, g_xz_y_0_x, g_xz_y_0_y, g_xz_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_y_0_x[i] = 2.0 * g_xz_y_0_x[i] * a_exp;

        g_x_0_0_0_z_y_0_y[i] = 2.0 * g_xz_y_0_y[i] * a_exp;

        g_x_0_0_0_z_y_0_z[i] = 2.0 * g_xz_y_0_z[i] * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_0_0_0_z_z_0_x, g_x_0_0_0_z_z_0_y, g_x_0_0_0_z_z_0_z, g_xz_z_0_x, g_xz_z_0_y, g_xz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_z_z_0_x[i] = 2.0 * g_xz_z_0_x[i] * a_exp;

        g_x_0_0_0_z_z_0_y[i] = 2.0 * g_xz_z_0_y[i] * a_exp;

        g_x_0_0_0_z_z_0_z[i] = 2.0 * g_xz_z_0_z[i] * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_xy_x_0_x, g_xy_x_0_y, g_xy_x_0_z, g_y_0_0_0_x_x_0_x, g_y_0_0_0_x_x_0_y, g_y_0_0_0_x_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_x_0_x[i] = 2.0 * g_xy_x_0_x[i] * a_exp;

        g_y_0_0_0_x_x_0_y[i] = 2.0 * g_xy_x_0_y[i] * a_exp;

        g_y_0_0_0_x_x_0_z[i] = 2.0 * g_xy_x_0_z[i] * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_xy_y_0_x, g_xy_y_0_y, g_xy_y_0_z, g_y_0_0_0_x_y_0_x, g_y_0_0_0_x_y_0_y, g_y_0_0_0_x_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_y_0_x[i] = 2.0 * g_xy_y_0_x[i] * a_exp;

        g_y_0_0_0_x_y_0_y[i] = 2.0 * g_xy_y_0_y[i] * a_exp;

        g_y_0_0_0_x_y_0_z[i] = 2.0 * g_xy_y_0_z[i] * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_xy_z_0_x, g_xy_z_0_y, g_xy_z_0_z, g_y_0_0_0_x_z_0_x, g_y_0_0_0_x_z_0_y, g_y_0_0_0_x_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_x_z_0_x[i] = 2.0 * g_xy_z_0_x[i] * a_exp;

        g_y_0_0_0_x_z_0_y[i] = 2.0 * g_xy_z_0_y[i] * a_exp;

        g_y_0_0_0_x_z_0_z[i] = 2.0 * g_xy_z_0_z[i] * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_y_0_0_0_y_x_0_x, g_y_0_0_0_y_x_0_y, g_y_0_0_0_y_x_0_z, g_yy_x_0_x, g_yy_x_0_y, g_yy_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_x_0_x[i] = -g_0_x_0_x[i] + 2.0 * g_yy_x_0_x[i] * a_exp;

        g_y_0_0_0_y_x_0_y[i] = -g_0_x_0_y[i] + 2.0 * g_yy_x_0_y[i] * a_exp;

        g_y_0_0_0_y_x_0_z[i] = -g_0_x_0_z[i] + 2.0 * g_yy_x_0_z[i] * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_y_0_0_0_y_y_0_x, g_y_0_0_0_y_y_0_y, g_y_0_0_0_y_y_0_z, g_yy_y_0_x, g_yy_y_0_y, g_yy_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_y_0_x[i] = -g_0_y_0_x[i] + 2.0 * g_yy_y_0_x[i] * a_exp;

        g_y_0_0_0_y_y_0_y[i] = -g_0_y_0_y[i] + 2.0 * g_yy_y_0_y[i] * a_exp;

        g_y_0_0_0_y_y_0_z[i] = -g_0_y_0_z[i] + 2.0 * g_yy_y_0_z[i] * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_y_0_0_0_y_z_0_x, g_y_0_0_0_y_z_0_y, g_y_0_0_0_y_z_0_z, g_yy_z_0_x, g_yy_z_0_y, g_yy_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_y_z_0_x[i] = -g_0_z_0_x[i] + 2.0 * g_yy_z_0_x[i] * a_exp;

        g_y_0_0_0_y_z_0_y[i] = -g_0_z_0_y[i] + 2.0 * g_yy_z_0_y[i] * a_exp;

        g_y_0_0_0_y_z_0_z[i] = -g_0_z_0_z[i] + 2.0 * g_yy_z_0_z[i] * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_y_0_0_0_z_x_0_x, g_y_0_0_0_z_x_0_y, g_y_0_0_0_z_x_0_z, g_yz_x_0_x, g_yz_x_0_y, g_yz_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_x_0_x[i] = 2.0 * g_yz_x_0_x[i] * a_exp;

        g_y_0_0_0_z_x_0_y[i] = 2.0 * g_yz_x_0_y[i] * a_exp;

        g_y_0_0_0_z_x_0_z[i] = 2.0 * g_yz_x_0_z[i] * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_y_0_0_0_z_y_0_x, g_y_0_0_0_z_y_0_y, g_y_0_0_0_z_y_0_z, g_yz_y_0_x, g_yz_y_0_y, g_yz_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_y_0_x[i] = 2.0 * g_yz_y_0_x[i] * a_exp;

        g_y_0_0_0_z_y_0_y[i] = 2.0 * g_yz_y_0_y[i] * a_exp;

        g_y_0_0_0_z_y_0_z[i] = 2.0 * g_yz_y_0_z[i] * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_y_0_0_0_z_z_0_x, g_y_0_0_0_z_z_0_y, g_y_0_0_0_z_z_0_z, g_yz_z_0_x, g_yz_z_0_y, g_yz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_z_z_0_x[i] = 2.0 * g_yz_z_0_x[i] * a_exp;

        g_y_0_0_0_z_z_0_y[i] = 2.0 * g_yz_z_0_y[i] * a_exp;

        g_y_0_0_0_z_z_0_z[i] = 2.0 * g_yz_z_0_z[i] * a_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_xz_x_0_x, g_xz_x_0_y, g_xz_x_0_z, g_z_0_0_0_x_x_0_x, g_z_0_0_0_x_x_0_y, g_z_0_0_0_x_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_x_0_x[i] = 2.0 * g_xz_x_0_x[i] * a_exp;

        g_z_0_0_0_x_x_0_y[i] = 2.0 * g_xz_x_0_y[i] * a_exp;

        g_z_0_0_0_x_x_0_z[i] = 2.0 * g_xz_x_0_z[i] * a_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_xz_y_0_x, g_xz_y_0_y, g_xz_y_0_z, g_z_0_0_0_x_y_0_x, g_z_0_0_0_x_y_0_y, g_z_0_0_0_x_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_y_0_x[i] = 2.0 * g_xz_y_0_x[i] * a_exp;

        g_z_0_0_0_x_y_0_y[i] = 2.0 * g_xz_y_0_y[i] * a_exp;

        g_z_0_0_0_x_y_0_z[i] = 2.0 * g_xz_y_0_z[i] * a_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_xz_z_0_x, g_xz_z_0_y, g_xz_z_0_z, g_z_0_0_0_x_z_0_x, g_z_0_0_0_x_z_0_y, g_z_0_0_0_x_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_x_z_0_x[i] = 2.0 * g_xz_z_0_x[i] * a_exp;

        g_z_0_0_0_x_z_0_y[i] = 2.0 * g_xz_z_0_y[i] * a_exp;

        g_z_0_0_0_x_z_0_z[i] = 2.0 * g_xz_z_0_z[i] * a_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_yz_x_0_x, g_yz_x_0_y, g_yz_x_0_z, g_z_0_0_0_y_x_0_x, g_z_0_0_0_y_x_0_y, g_z_0_0_0_y_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_x_0_x[i] = 2.0 * g_yz_x_0_x[i] * a_exp;

        g_z_0_0_0_y_x_0_y[i] = 2.0 * g_yz_x_0_y[i] * a_exp;

        g_z_0_0_0_y_x_0_z[i] = 2.0 * g_yz_x_0_z[i] * a_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_yz_y_0_x, g_yz_y_0_y, g_yz_y_0_z, g_z_0_0_0_y_y_0_x, g_z_0_0_0_y_y_0_y, g_z_0_0_0_y_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_y_0_x[i] = 2.0 * g_yz_y_0_x[i] * a_exp;

        g_z_0_0_0_y_y_0_y[i] = 2.0 * g_yz_y_0_y[i] * a_exp;

        g_z_0_0_0_y_y_0_z[i] = 2.0 * g_yz_y_0_z[i] * a_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_yz_z_0_x, g_yz_z_0_y, g_yz_z_0_z, g_z_0_0_0_y_z_0_x, g_z_0_0_0_y_z_0_y, g_z_0_0_0_y_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_y_z_0_x[i] = 2.0 * g_yz_z_0_x[i] * a_exp;

        g_z_0_0_0_y_z_0_y[i] = 2.0 * g_yz_z_0_y[i] * a_exp;

        g_z_0_0_0_y_z_0_z[i] = 2.0 * g_yz_z_0_z[i] * a_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_0_x_0_x, g_0_x_0_y, g_0_x_0_z, g_z_0_0_0_z_x_0_x, g_z_0_0_0_z_x_0_y, g_z_0_0_0_z_x_0_z, g_zz_x_0_x, g_zz_x_0_y, g_zz_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_x_0_x[i] = -g_0_x_0_x[i] + 2.0 * g_zz_x_0_x[i] * a_exp;

        g_z_0_0_0_z_x_0_y[i] = -g_0_x_0_y[i] + 2.0 * g_zz_x_0_y[i] * a_exp;

        g_z_0_0_0_z_x_0_z[i] = -g_0_x_0_z[i] + 2.0 * g_zz_x_0_z[i] * a_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_0_y_0_x, g_0_y_0_y, g_0_y_0_z, g_z_0_0_0_z_y_0_x, g_z_0_0_0_z_y_0_y, g_z_0_0_0_z_y_0_z, g_zz_y_0_x, g_zz_y_0_y, g_zz_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_y_0_x[i] = -g_0_y_0_x[i] + 2.0 * g_zz_y_0_x[i] * a_exp;

        g_z_0_0_0_z_y_0_y[i] = -g_0_y_0_y[i] + 2.0 * g_zz_y_0_y[i] * a_exp;

        g_z_0_0_0_z_y_0_z[i] = -g_0_y_0_z[i] + 2.0 * g_zz_y_0_z[i] * a_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_0_z_0_x, g_0_z_0_y, g_0_z_0_z, g_z_0_0_0_z_z_0_x, g_z_0_0_0_z_z_0_y, g_z_0_0_0_z_z_0_z, g_zz_z_0_x, g_zz_z_0_y, g_zz_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_z_z_0_x[i] = -g_0_z_0_x[i] + 2.0 * g_zz_z_0_x[i] * a_exp;

        g_z_0_0_0_z_z_0_y[i] = -g_0_z_0_y[i] + 2.0 * g_zz_z_0_y[i] * a_exp;

        g_z_0_0_0_z_z_0_z[i] = -g_0_z_0_z[i] + 2.0 * g_zz_z_0_z[i] * a_exp;
    }
}

} // t4c_geom namespace

