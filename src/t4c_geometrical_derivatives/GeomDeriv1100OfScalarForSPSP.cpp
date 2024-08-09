#include "GeomDeriv1100OfScalarForSPSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_spsp_0(CSimdArray<double>& buffer_1100_spsp,
                     const CSimdArray<double>& buffer_pssp,
                     const CSimdArray<double>& buffer_pdsp,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_spsp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pssp

    auto g_x_0_0_x = buffer_pssp[0];

    auto g_x_0_0_y = buffer_pssp[1];

    auto g_x_0_0_z = buffer_pssp[2];

    auto g_y_0_0_x = buffer_pssp[3];

    auto g_y_0_0_y = buffer_pssp[4];

    auto g_y_0_0_z = buffer_pssp[5];

    auto g_z_0_0_x = buffer_pssp[6];

    auto g_z_0_0_y = buffer_pssp[7];

    auto g_z_0_0_z = buffer_pssp[8];

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

    /// Set up components of integrals buffer : buffer_1100_spsp

    auto g_x_x_0_0_0_x_0_x = buffer_1100_spsp[0];

    auto g_x_x_0_0_0_x_0_y = buffer_1100_spsp[1];

    auto g_x_x_0_0_0_x_0_z = buffer_1100_spsp[2];

    auto g_x_x_0_0_0_y_0_x = buffer_1100_spsp[3];

    auto g_x_x_0_0_0_y_0_y = buffer_1100_spsp[4];

    auto g_x_x_0_0_0_y_0_z = buffer_1100_spsp[5];

    auto g_x_x_0_0_0_z_0_x = buffer_1100_spsp[6];

    auto g_x_x_0_0_0_z_0_y = buffer_1100_spsp[7];

    auto g_x_x_0_0_0_z_0_z = buffer_1100_spsp[8];

    auto g_x_y_0_0_0_x_0_x = buffer_1100_spsp[9];

    auto g_x_y_0_0_0_x_0_y = buffer_1100_spsp[10];

    auto g_x_y_0_0_0_x_0_z = buffer_1100_spsp[11];

    auto g_x_y_0_0_0_y_0_x = buffer_1100_spsp[12];

    auto g_x_y_0_0_0_y_0_y = buffer_1100_spsp[13];

    auto g_x_y_0_0_0_y_0_z = buffer_1100_spsp[14];

    auto g_x_y_0_0_0_z_0_x = buffer_1100_spsp[15];

    auto g_x_y_0_0_0_z_0_y = buffer_1100_spsp[16];

    auto g_x_y_0_0_0_z_0_z = buffer_1100_spsp[17];

    auto g_x_z_0_0_0_x_0_x = buffer_1100_spsp[18];

    auto g_x_z_0_0_0_x_0_y = buffer_1100_spsp[19];

    auto g_x_z_0_0_0_x_0_z = buffer_1100_spsp[20];

    auto g_x_z_0_0_0_y_0_x = buffer_1100_spsp[21];

    auto g_x_z_0_0_0_y_0_y = buffer_1100_spsp[22];

    auto g_x_z_0_0_0_y_0_z = buffer_1100_spsp[23];

    auto g_x_z_0_0_0_z_0_x = buffer_1100_spsp[24];

    auto g_x_z_0_0_0_z_0_y = buffer_1100_spsp[25];

    auto g_x_z_0_0_0_z_0_z = buffer_1100_spsp[26];

    auto g_y_x_0_0_0_x_0_x = buffer_1100_spsp[27];

    auto g_y_x_0_0_0_x_0_y = buffer_1100_spsp[28];

    auto g_y_x_0_0_0_x_0_z = buffer_1100_spsp[29];

    auto g_y_x_0_0_0_y_0_x = buffer_1100_spsp[30];

    auto g_y_x_0_0_0_y_0_y = buffer_1100_spsp[31];

    auto g_y_x_0_0_0_y_0_z = buffer_1100_spsp[32];

    auto g_y_x_0_0_0_z_0_x = buffer_1100_spsp[33];

    auto g_y_x_0_0_0_z_0_y = buffer_1100_spsp[34];

    auto g_y_x_0_0_0_z_0_z = buffer_1100_spsp[35];

    auto g_y_y_0_0_0_x_0_x = buffer_1100_spsp[36];

    auto g_y_y_0_0_0_x_0_y = buffer_1100_spsp[37];

    auto g_y_y_0_0_0_x_0_z = buffer_1100_spsp[38];

    auto g_y_y_0_0_0_y_0_x = buffer_1100_spsp[39];

    auto g_y_y_0_0_0_y_0_y = buffer_1100_spsp[40];

    auto g_y_y_0_0_0_y_0_z = buffer_1100_spsp[41];

    auto g_y_y_0_0_0_z_0_x = buffer_1100_spsp[42];

    auto g_y_y_0_0_0_z_0_y = buffer_1100_spsp[43];

    auto g_y_y_0_0_0_z_0_z = buffer_1100_spsp[44];

    auto g_y_z_0_0_0_x_0_x = buffer_1100_spsp[45];

    auto g_y_z_0_0_0_x_0_y = buffer_1100_spsp[46];

    auto g_y_z_0_0_0_x_0_z = buffer_1100_spsp[47];

    auto g_y_z_0_0_0_y_0_x = buffer_1100_spsp[48];

    auto g_y_z_0_0_0_y_0_y = buffer_1100_spsp[49];

    auto g_y_z_0_0_0_y_0_z = buffer_1100_spsp[50];

    auto g_y_z_0_0_0_z_0_x = buffer_1100_spsp[51];

    auto g_y_z_0_0_0_z_0_y = buffer_1100_spsp[52];

    auto g_y_z_0_0_0_z_0_z = buffer_1100_spsp[53];

    auto g_z_x_0_0_0_x_0_x = buffer_1100_spsp[54];

    auto g_z_x_0_0_0_x_0_y = buffer_1100_spsp[55];

    auto g_z_x_0_0_0_x_0_z = buffer_1100_spsp[56];

    auto g_z_x_0_0_0_y_0_x = buffer_1100_spsp[57];

    auto g_z_x_0_0_0_y_0_y = buffer_1100_spsp[58];

    auto g_z_x_0_0_0_y_0_z = buffer_1100_spsp[59];

    auto g_z_x_0_0_0_z_0_x = buffer_1100_spsp[60];

    auto g_z_x_0_0_0_z_0_y = buffer_1100_spsp[61];

    auto g_z_x_0_0_0_z_0_z = buffer_1100_spsp[62];

    auto g_z_y_0_0_0_x_0_x = buffer_1100_spsp[63];

    auto g_z_y_0_0_0_x_0_y = buffer_1100_spsp[64];

    auto g_z_y_0_0_0_x_0_z = buffer_1100_spsp[65];

    auto g_z_y_0_0_0_y_0_x = buffer_1100_spsp[66];

    auto g_z_y_0_0_0_y_0_y = buffer_1100_spsp[67];

    auto g_z_y_0_0_0_y_0_z = buffer_1100_spsp[68];

    auto g_z_y_0_0_0_z_0_x = buffer_1100_spsp[69];

    auto g_z_y_0_0_0_z_0_y = buffer_1100_spsp[70];

    auto g_z_y_0_0_0_z_0_z = buffer_1100_spsp[71];

    auto g_z_z_0_0_0_x_0_x = buffer_1100_spsp[72];

    auto g_z_z_0_0_0_x_0_y = buffer_1100_spsp[73];

    auto g_z_z_0_0_0_x_0_z = buffer_1100_spsp[74];

    auto g_z_z_0_0_0_y_0_x = buffer_1100_spsp[75];

    auto g_z_z_0_0_0_y_0_y = buffer_1100_spsp[76];

    auto g_z_z_0_0_0_y_0_z = buffer_1100_spsp[77];

    auto g_z_z_0_0_0_z_0_x = buffer_1100_spsp[78];

    auto g_z_z_0_0_0_z_0_y = buffer_1100_spsp[79];

    auto g_z_z_0_0_0_z_0_z = buffer_1100_spsp[80];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_0_x, g_x_0_0_y, g_x_0_0_z, g_x_x_0_0_0_x_0_x, g_x_x_0_0_0_x_0_y, g_x_x_0_0_0_x_0_z, g_x_xx_0_x, g_x_xx_0_y, g_x_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_x_0_x[i] = -2.0 * g_x_0_0_x[i] * a_exp + 4.0 * g_x_xx_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_0_y[i] = -2.0 * g_x_0_0_y[i] * a_exp + 4.0 * g_x_xx_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_x_0_z[i] = -2.0 * g_x_0_0_z[i] * a_exp + 4.0 * g_x_xx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_x_0_0_0_y_0_x, g_x_x_0_0_0_y_0_y, g_x_x_0_0_0_y_0_z, g_x_xy_0_x, g_x_xy_0_y, g_x_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_y_0_x[i] = 4.0 * g_x_xy_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_0_y[i] = 4.0 * g_x_xy_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_y_0_z[i] = 4.0 * g_x_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_x_0_0_0_z_0_x, g_x_x_0_0_0_z_0_y, g_x_x_0_0_0_z_0_z, g_x_xz_0_x, g_x_xz_0_y, g_x_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_z_0_x[i] = 4.0 * g_x_xz_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_0_y[i] = 4.0 * g_x_xz_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_z_0_z[i] = 4.0 * g_x_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_xy_0_x, g_x_xy_0_y, g_x_xy_0_z, g_x_y_0_0_0_x_0_x, g_x_y_0_0_0_x_0_y, g_x_y_0_0_0_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_x_0_x[i] = 4.0 * g_x_xy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_0_y[i] = 4.0 * g_x_xy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_x_0_z[i] = 4.0 * g_x_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_0_0_x, g_x_0_0_y, g_x_0_0_z, g_x_y_0_0_0_y_0_x, g_x_y_0_0_0_y_0_y, g_x_y_0_0_0_y_0_z, g_x_yy_0_x, g_x_yy_0_y, g_x_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_y_0_x[i] = -2.0 * g_x_0_0_x[i] * a_exp + 4.0 * g_x_yy_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_0_y[i] = -2.0 * g_x_0_0_y[i] * a_exp + 4.0 * g_x_yy_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_y_0_z[i] = -2.0 * g_x_0_0_z[i] * a_exp + 4.0 * g_x_yy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_y_0_0_0_z_0_x, g_x_y_0_0_0_z_0_y, g_x_y_0_0_0_z_0_z, g_x_yz_0_x, g_x_yz_0_y, g_x_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_z_0_x[i] = 4.0 * g_x_yz_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_0_y[i] = 4.0 * g_x_yz_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_z_0_z[i] = 4.0 * g_x_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_xz_0_x, g_x_xz_0_y, g_x_xz_0_z, g_x_z_0_0_0_x_0_x, g_x_z_0_0_0_x_0_y, g_x_z_0_0_0_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_x_0_x[i] = 4.0 * g_x_xz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_0_y[i] = 4.0 * g_x_xz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_x_0_z[i] = 4.0 * g_x_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_yz_0_x, g_x_yz_0_y, g_x_yz_0_z, g_x_z_0_0_0_y_0_x, g_x_z_0_0_0_y_0_y, g_x_z_0_0_0_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_y_0_x[i] = 4.0 * g_x_yz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_0_y[i] = 4.0 * g_x_yz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_y_0_z[i] = 4.0 * g_x_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_0_0_x, g_x_0_0_y, g_x_0_0_z, g_x_z_0_0_0_z_0_x, g_x_z_0_0_0_z_0_y, g_x_z_0_0_0_z_0_z, g_x_zz_0_x, g_x_zz_0_y, g_x_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_z_0_x[i] = -2.0 * g_x_0_0_x[i] * a_exp + 4.0 * g_x_zz_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_0_y[i] = -2.0 * g_x_0_0_y[i] * a_exp + 4.0 * g_x_zz_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_z_0_z[i] = -2.0 * g_x_0_0_z[i] * a_exp + 4.0 * g_x_zz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_y_0_0_x, g_y_0_0_y, g_y_0_0_z, g_y_x_0_0_0_x_0_x, g_y_x_0_0_0_x_0_y, g_y_x_0_0_0_x_0_z, g_y_xx_0_x, g_y_xx_0_y, g_y_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_x_0_x[i] = -2.0 * g_y_0_0_x[i] * a_exp + 4.0 * g_y_xx_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_0_y[i] = -2.0 * g_y_0_0_y[i] * a_exp + 4.0 * g_y_xx_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_x_0_z[i] = -2.0 * g_y_0_0_z[i] * a_exp + 4.0 * g_y_xx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_y_x_0_0_0_y_0_x, g_y_x_0_0_0_y_0_y, g_y_x_0_0_0_y_0_z, g_y_xy_0_x, g_y_xy_0_y, g_y_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_y_0_x[i] = 4.0 * g_y_xy_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_0_y[i] = 4.0 * g_y_xy_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_y_0_z[i] = 4.0 * g_y_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_y_x_0_0_0_z_0_x, g_y_x_0_0_0_z_0_y, g_y_x_0_0_0_z_0_z, g_y_xz_0_x, g_y_xz_0_y, g_y_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_z_0_x[i] = 4.0 * g_y_xz_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_0_y[i] = 4.0 * g_y_xz_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_z_0_z[i] = 4.0 * g_y_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_y_xy_0_x, g_y_xy_0_y, g_y_xy_0_z, g_y_y_0_0_0_x_0_x, g_y_y_0_0_0_x_0_y, g_y_y_0_0_0_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_x_0_x[i] = 4.0 * g_y_xy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_0_y[i] = 4.0 * g_y_xy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_x_0_z[i] = 4.0 * g_y_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_y_0_0_x, g_y_0_0_y, g_y_0_0_z, g_y_y_0_0_0_y_0_x, g_y_y_0_0_0_y_0_y, g_y_y_0_0_0_y_0_z, g_y_yy_0_x, g_y_yy_0_y, g_y_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_y_0_x[i] = -2.0 * g_y_0_0_x[i] * a_exp + 4.0 * g_y_yy_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_0_y[i] = -2.0 * g_y_0_0_y[i] * a_exp + 4.0 * g_y_yy_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_y_0_z[i] = -2.0 * g_y_0_0_z[i] * a_exp + 4.0 * g_y_yy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_y_y_0_0_0_z_0_x, g_y_y_0_0_0_z_0_y, g_y_y_0_0_0_z_0_z, g_y_yz_0_x, g_y_yz_0_y, g_y_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_z_0_x[i] = 4.0 * g_y_yz_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_0_y[i] = 4.0 * g_y_yz_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_z_0_z[i] = 4.0 * g_y_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_y_xz_0_x, g_y_xz_0_y, g_y_xz_0_z, g_y_z_0_0_0_x_0_x, g_y_z_0_0_0_x_0_y, g_y_z_0_0_0_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_x_0_x[i] = 4.0 * g_y_xz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_0_y[i] = 4.0 * g_y_xz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_x_0_z[i] = 4.0 * g_y_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_y_yz_0_x, g_y_yz_0_y, g_y_yz_0_z, g_y_z_0_0_0_y_0_x, g_y_z_0_0_0_y_0_y, g_y_z_0_0_0_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_y_0_x[i] = 4.0 * g_y_yz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_0_y[i] = 4.0 * g_y_yz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_y_0_z[i] = 4.0 * g_y_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_y_0_0_x, g_y_0_0_y, g_y_0_0_z, g_y_z_0_0_0_z_0_x, g_y_z_0_0_0_z_0_y, g_y_z_0_0_0_z_0_z, g_y_zz_0_x, g_y_zz_0_y, g_y_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_z_0_x[i] = -2.0 * g_y_0_0_x[i] * a_exp + 4.0 * g_y_zz_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_0_y[i] = -2.0 * g_y_0_0_y[i] * a_exp + 4.0 * g_y_zz_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_z_0_z[i] = -2.0 * g_y_0_0_z[i] * a_exp + 4.0 * g_y_zz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_z_0_0_x, g_z_0_0_y, g_z_0_0_z, g_z_x_0_0_0_x_0_x, g_z_x_0_0_0_x_0_y, g_z_x_0_0_0_x_0_z, g_z_xx_0_x, g_z_xx_0_y, g_z_xx_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_x_0_x[i] = -2.0 * g_z_0_0_x[i] * a_exp + 4.0 * g_z_xx_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_0_y[i] = -2.0 * g_z_0_0_y[i] * a_exp + 4.0 * g_z_xx_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_x_0_z[i] = -2.0 * g_z_0_0_z[i] * a_exp + 4.0 * g_z_xx_0_z[i] * a_exp * b_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_z_x_0_0_0_y_0_x, g_z_x_0_0_0_y_0_y, g_z_x_0_0_0_y_0_z, g_z_xy_0_x, g_z_xy_0_y, g_z_xy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_y_0_x[i] = 4.0 * g_z_xy_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_0_y[i] = 4.0 * g_z_xy_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_y_0_z[i] = 4.0 * g_z_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_z_x_0_0_0_z_0_x, g_z_x_0_0_0_z_0_y, g_z_x_0_0_0_z_0_z, g_z_xz_0_x, g_z_xz_0_y, g_z_xz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_z_0_x[i] = 4.0 * g_z_xz_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_0_y[i] = 4.0 * g_z_xz_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_z_0_z[i] = 4.0 * g_z_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_z_xy_0_x, g_z_xy_0_y, g_z_xy_0_z, g_z_y_0_0_0_x_0_x, g_z_y_0_0_0_x_0_y, g_z_y_0_0_0_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_x_0_x[i] = 4.0 * g_z_xy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_0_y[i] = 4.0 * g_z_xy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_x_0_z[i] = 4.0 * g_z_xy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_z_0_0_x, g_z_0_0_y, g_z_0_0_z, g_z_y_0_0_0_y_0_x, g_z_y_0_0_0_y_0_y, g_z_y_0_0_0_y_0_z, g_z_yy_0_x, g_z_yy_0_y, g_z_yy_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_y_0_x[i] = -2.0 * g_z_0_0_x[i] * a_exp + 4.0 * g_z_yy_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_0_y[i] = -2.0 * g_z_0_0_y[i] * a_exp + 4.0 * g_z_yy_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_y_0_z[i] = -2.0 * g_z_0_0_z[i] * a_exp + 4.0 * g_z_yy_0_z[i] * a_exp * b_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_z_y_0_0_0_z_0_x, g_z_y_0_0_0_z_0_y, g_z_y_0_0_0_z_0_z, g_z_yz_0_x, g_z_yz_0_y, g_z_yz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_z_0_x[i] = 4.0 * g_z_yz_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_0_y[i] = 4.0 * g_z_yz_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_z_0_z[i] = 4.0 * g_z_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_z_xz_0_x, g_z_xz_0_y, g_z_xz_0_z, g_z_z_0_0_0_x_0_x, g_z_z_0_0_0_x_0_y, g_z_z_0_0_0_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_x_0_x[i] = 4.0 * g_z_xz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_0_y[i] = 4.0 * g_z_xz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_x_0_z[i] = 4.0 * g_z_xz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_z_yz_0_x, g_z_yz_0_y, g_z_yz_0_z, g_z_z_0_0_0_y_0_x, g_z_z_0_0_0_y_0_y, g_z_z_0_0_0_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_y_0_x[i] = 4.0 * g_z_yz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_0_y[i] = 4.0 * g_z_yz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_y_0_z[i] = 4.0 * g_z_yz_0_z[i] * a_exp * b_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_z_0_0_x, g_z_0_0_y, g_z_0_0_z, g_z_z_0_0_0_z_0_x, g_z_z_0_0_0_z_0_y, g_z_z_0_0_0_z_0_z, g_z_zz_0_x, g_z_zz_0_y, g_z_zz_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_z_0_x[i] = -2.0 * g_z_0_0_x[i] * a_exp + 4.0 * g_z_zz_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_0_y[i] = -2.0 * g_z_0_0_y[i] * a_exp + 4.0 * g_z_zz_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_z_0_z[i] = -2.0 * g_z_0_0_z[i] * a_exp + 4.0 * g_z_zz_0_z[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

