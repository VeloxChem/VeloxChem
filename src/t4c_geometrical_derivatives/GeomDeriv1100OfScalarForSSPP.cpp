#include "GeomDeriv1100OfScalarForSSPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_sspp_0(CSimdArray<double>& buffer_1100_sspp,
                     const CSimdArray<double>& buffer_pppp,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_sspp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pppp

    auto g_x_x_x_x = buffer_pppp[0];

    auto g_x_x_x_y = buffer_pppp[1];

    auto g_x_x_x_z = buffer_pppp[2];

    auto g_x_x_y_x = buffer_pppp[3];

    auto g_x_x_y_y = buffer_pppp[4];

    auto g_x_x_y_z = buffer_pppp[5];

    auto g_x_x_z_x = buffer_pppp[6];

    auto g_x_x_z_y = buffer_pppp[7];

    auto g_x_x_z_z = buffer_pppp[8];

    auto g_x_y_x_x = buffer_pppp[9];

    auto g_x_y_x_y = buffer_pppp[10];

    auto g_x_y_x_z = buffer_pppp[11];

    auto g_x_y_y_x = buffer_pppp[12];

    auto g_x_y_y_y = buffer_pppp[13];

    auto g_x_y_y_z = buffer_pppp[14];

    auto g_x_y_z_x = buffer_pppp[15];

    auto g_x_y_z_y = buffer_pppp[16];

    auto g_x_y_z_z = buffer_pppp[17];

    auto g_x_z_x_x = buffer_pppp[18];

    auto g_x_z_x_y = buffer_pppp[19];

    auto g_x_z_x_z = buffer_pppp[20];

    auto g_x_z_y_x = buffer_pppp[21];

    auto g_x_z_y_y = buffer_pppp[22];

    auto g_x_z_y_z = buffer_pppp[23];

    auto g_x_z_z_x = buffer_pppp[24];

    auto g_x_z_z_y = buffer_pppp[25];

    auto g_x_z_z_z = buffer_pppp[26];

    auto g_y_x_x_x = buffer_pppp[27];

    auto g_y_x_x_y = buffer_pppp[28];

    auto g_y_x_x_z = buffer_pppp[29];

    auto g_y_x_y_x = buffer_pppp[30];

    auto g_y_x_y_y = buffer_pppp[31];

    auto g_y_x_y_z = buffer_pppp[32];

    auto g_y_x_z_x = buffer_pppp[33];

    auto g_y_x_z_y = buffer_pppp[34];

    auto g_y_x_z_z = buffer_pppp[35];

    auto g_y_y_x_x = buffer_pppp[36];

    auto g_y_y_x_y = buffer_pppp[37];

    auto g_y_y_x_z = buffer_pppp[38];

    auto g_y_y_y_x = buffer_pppp[39];

    auto g_y_y_y_y = buffer_pppp[40];

    auto g_y_y_y_z = buffer_pppp[41];

    auto g_y_y_z_x = buffer_pppp[42];

    auto g_y_y_z_y = buffer_pppp[43];

    auto g_y_y_z_z = buffer_pppp[44];

    auto g_y_z_x_x = buffer_pppp[45];

    auto g_y_z_x_y = buffer_pppp[46];

    auto g_y_z_x_z = buffer_pppp[47];

    auto g_y_z_y_x = buffer_pppp[48];

    auto g_y_z_y_y = buffer_pppp[49];

    auto g_y_z_y_z = buffer_pppp[50];

    auto g_y_z_z_x = buffer_pppp[51];

    auto g_y_z_z_y = buffer_pppp[52];

    auto g_y_z_z_z = buffer_pppp[53];

    auto g_z_x_x_x = buffer_pppp[54];

    auto g_z_x_x_y = buffer_pppp[55];

    auto g_z_x_x_z = buffer_pppp[56];

    auto g_z_x_y_x = buffer_pppp[57];

    auto g_z_x_y_y = buffer_pppp[58];

    auto g_z_x_y_z = buffer_pppp[59];

    auto g_z_x_z_x = buffer_pppp[60];

    auto g_z_x_z_y = buffer_pppp[61];

    auto g_z_x_z_z = buffer_pppp[62];

    auto g_z_y_x_x = buffer_pppp[63];

    auto g_z_y_x_y = buffer_pppp[64];

    auto g_z_y_x_z = buffer_pppp[65];

    auto g_z_y_y_x = buffer_pppp[66];

    auto g_z_y_y_y = buffer_pppp[67];

    auto g_z_y_y_z = buffer_pppp[68];

    auto g_z_y_z_x = buffer_pppp[69];

    auto g_z_y_z_y = buffer_pppp[70];

    auto g_z_y_z_z = buffer_pppp[71];

    auto g_z_z_x_x = buffer_pppp[72];

    auto g_z_z_x_y = buffer_pppp[73];

    auto g_z_z_x_z = buffer_pppp[74];

    auto g_z_z_y_x = buffer_pppp[75];

    auto g_z_z_y_y = buffer_pppp[76];

    auto g_z_z_y_z = buffer_pppp[77];

    auto g_z_z_z_x = buffer_pppp[78];

    auto g_z_z_z_y = buffer_pppp[79];

    auto g_z_z_z_z = buffer_pppp[80];

    /// Set up components of integrals buffer : buffer_1100_sspp

    auto g_x_x_0_0_0_0_x_x = buffer_1100_sspp[0];

    auto g_x_x_0_0_0_0_x_y = buffer_1100_sspp[1];

    auto g_x_x_0_0_0_0_x_z = buffer_1100_sspp[2];

    auto g_x_x_0_0_0_0_y_x = buffer_1100_sspp[3];

    auto g_x_x_0_0_0_0_y_y = buffer_1100_sspp[4];

    auto g_x_x_0_0_0_0_y_z = buffer_1100_sspp[5];

    auto g_x_x_0_0_0_0_z_x = buffer_1100_sspp[6];

    auto g_x_x_0_0_0_0_z_y = buffer_1100_sspp[7];

    auto g_x_x_0_0_0_0_z_z = buffer_1100_sspp[8];

    auto g_x_y_0_0_0_0_x_x = buffer_1100_sspp[9];

    auto g_x_y_0_0_0_0_x_y = buffer_1100_sspp[10];

    auto g_x_y_0_0_0_0_x_z = buffer_1100_sspp[11];

    auto g_x_y_0_0_0_0_y_x = buffer_1100_sspp[12];

    auto g_x_y_0_0_0_0_y_y = buffer_1100_sspp[13];

    auto g_x_y_0_0_0_0_y_z = buffer_1100_sspp[14];

    auto g_x_y_0_0_0_0_z_x = buffer_1100_sspp[15];

    auto g_x_y_0_0_0_0_z_y = buffer_1100_sspp[16];

    auto g_x_y_0_0_0_0_z_z = buffer_1100_sspp[17];

    auto g_x_z_0_0_0_0_x_x = buffer_1100_sspp[18];

    auto g_x_z_0_0_0_0_x_y = buffer_1100_sspp[19];

    auto g_x_z_0_0_0_0_x_z = buffer_1100_sspp[20];

    auto g_x_z_0_0_0_0_y_x = buffer_1100_sspp[21];

    auto g_x_z_0_0_0_0_y_y = buffer_1100_sspp[22];

    auto g_x_z_0_0_0_0_y_z = buffer_1100_sspp[23];

    auto g_x_z_0_0_0_0_z_x = buffer_1100_sspp[24];

    auto g_x_z_0_0_0_0_z_y = buffer_1100_sspp[25];

    auto g_x_z_0_0_0_0_z_z = buffer_1100_sspp[26];

    auto g_y_x_0_0_0_0_x_x = buffer_1100_sspp[27];

    auto g_y_x_0_0_0_0_x_y = buffer_1100_sspp[28];

    auto g_y_x_0_0_0_0_x_z = buffer_1100_sspp[29];

    auto g_y_x_0_0_0_0_y_x = buffer_1100_sspp[30];

    auto g_y_x_0_0_0_0_y_y = buffer_1100_sspp[31];

    auto g_y_x_0_0_0_0_y_z = buffer_1100_sspp[32];

    auto g_y_x_0_0_0_0_z_x = buffer_1100_sspp[33];

    auto g_y_x_0_0_0_0_z_y = buffer_1100_sspp[34];

    auto g_y_x_0_0_0_0_z_z = buffer_1100_sspp[35];

    auto g_y_y_0_0_0_0_x_x = buffer_1100_sspp[36];

    auto g_y_y_0_0_0_0_x_y = buffer_1100_sspp[37];

    auto g_y_y_0_0_0_0_x_z = buffer_1100_sspp[38];

    auto g_y_y_0_0_0_0_y_x = buffer_1100_sspp[39];

    auto g_y_y_0_0_0_0_y_y = buffer_1100_sspp[40];

    auto g_y_y_0_0_0_0_y_z = buffer_1100_sspp[41];

    auto g_y_y_0_0_0_0_z_x = buffer_1100_sspp[42];

    auto g_y_y_0_0_0_0_z_y = buffer_1100_sspp[43];

    auto g_y_y_0_0_0_0_z_z = buffer_1100_sspp[44];

    auto g_y_z_0_0_0_0_x_x = buffer_1100_sspp[45];

    auto g_y_z_0_0_0_0_x_y = buffer_1100_sspp[46];

    auto g_y_z_0_0_0_0_x_z = buffer_1100_sspp[47];

    auto g_y_z_0_0_0_0_y_x = buffer_1100_sspp[48];

    auto g_y_z_0_0_0_0_y_y = buffer_1100_sspp[49];

    auto g_y_z_0_0_0_0_y_z = buffer_1100_sspp[50];

    auto g_y_z_0_0_0_0_z_x = buffer_1100_sspp[51];

    auto g_y_z_0_0_0_0_z_y = buffer_1100_sspp[52];

    auto g_y_z_0_0_0_0_z_z = buffer_1100_sspp[53];

    auto g_z_x_0_0_0_0_x_x = buffer_1100_sspp[54];

    auto g_z_x_0_0_0_0_x_y = buffer_1100_sspp[55];

    auto g_z_x_0_0_0_0_x_z = buffer_1100_sspp[56];

    auto g_z_x_0_0_0_0_y_x = buffer_1100_sspp[57];

    auto g_z_x_0_0_0_0_y_y = buffer_1100_sspp[58];

    auto g_z_x_0_0_0_0_y_z = buffer_1100_sspp[59];

    auto g_z_x_0_0_0_0_z_x = buffer_1100_sspp[60];

    auto g_z_x_0_0_0_0_z_y = buffer_1100_sspp[61];

    auto g_z_x_0_0_0_0_z_z = buffer_1100_sspp[62];

    auto g_z_y_0_0_0_0_x_x = buffer_1100_sspp[63];

    auto g_z_y_0_0_0_0_x_y = buffer_1100_sspp[64];

    auto g_z_y_0_0_0_0_x_z = buffer_1100_sspp[65];

    auto g_z_y_0_0_0_0_y_x = buffer_1100_sspp[66];

    auto g_z_y_0_0_0_0_y_y = buffer_1100_sspp[67];

    auto g_z_y_0_0_0_0_y_z = buffer_1100_sspp[68];

    auto g_z_y_0_0_0_0_z_x = buffer_1100_sspp[69];

    auto g_z_y_0_0_0_0_z_y = buffer_1100_sspp[70];

    auto g_z_y_0_0_0_0_z_z = buffer_1100_sspp[71];

    auto g_z_z_0_0_0_0_x_x = buffer_1100_sspp[72];

    auto g_z_z_0_0_0_0_x_y = buffer_1100_sspp[73];

    auto g_z_z_0_0_0_0_x_z = buffer_1100_sspp[74];

    auto g_z_z_0_0_0_0_y_x = buffer_1100_sspp[75];

    auto g_z_z_0_0_0_0_y_y = buffer_1100_sspp[76];

    auto g_z_z_0_0_0_0_y_z = buffer_1100_sspp[77];

    auto g_z_z_0_0_0_0_z_x = buffer_1100_sspp[78];

    auto g_z_z_0_0_0_0_z_y = buffer_1100_sspp[79];

    auto g_z_z_0_0_0_0_z_z = buffer_1100_sspp[80];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_x_0_0_0_0_x_x, g_x_x_0_0_0_0_x_y, g_x_x_0_0_0_0_x_z, g_x_x_x_x, g_x_x_x_y, g_x_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_0_x_x[i] = 4.0 * g_x_x_x_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_x_y[i] = 4.0 * g_x_x_x_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_x_z[i] = 4.0 * g_x_x_x_z[i] * a_exp * b_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_x_0_0_0_0_y_x, g_x_x_0_0_0_0_y_y, g_x_x_0_0_0_0_y_z, g_x_x_y_x, g_x_x_y_y, g_x_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_0_y_x[i] = 4.0 * g_x_x_y_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_y_y[i] = 4.0 * g_x_x_y_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_y_z[i] = 4.0 * g_x_x_y_z[i] * a_exp * b_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_x_0_0_0_0_z_x, g_x_x_0_0_0_0_z_y, g_x_x_0_0_0_0_z_z, g_x_x_z_x, g_x_x_z_y, g_x_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_0_z_x[i] = 4.0 * g_x_x_z_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_z_y[i] = 4.0 * g_x_x_z_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_z_z[i] = 4.0 * g_x_x_z_z[i] * a_exp * b_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_y_0_0_0_0_x_x, g_x_y_0_0_0_0_x_y, g_x_y_0_0_0_0_x_z, g_x_y_x_x, g_x_y_x_y, g_x_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_0_x_x[i] = 4.0 * g_x_y_x_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_x_y[i] = 4.0 * g_x_y_x_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_x_z[i] = 4.0 * g_x_y_x_z[i] * a_exp * b_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_y_0_0_0_0_y_x, g_x_y_0_0_0_0_y_y, g_x_y_0_0_0_0_y_z, g_x_y_y_x, g_x_y_y_y, g_x_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_0_y_x[i] = 4.0 * g_x_y_y_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_y_y[i] = 4.0 * g_x_y_y_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_y_z[i] = 4.0 * g_x_y_y_z[i] * a_exp * b_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_y_0_0_0_0_z_x, g_x_y_0_0_0_0_z_y, g_x_y_0_0_0_0_z_z, g_x_y_z_x, g_x_y_z_y, g_x_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_0_z_x[i] = 4.0 * g_x_y_z_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_z_y[i] = 4.0 * g_x_y_z_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_z_z[i] = 4.0 * g_x_y_z_z[i] * a_exp * b_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_z_0_0_0_0_x_x, g_x_z_0_0_0_0_x_y, g_x_z_0_0_0_0_x_z, g_x_z_x_x, g_x_z_x_y, g_x_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_0_x_x[i] = 4.0 * g_x_z_x_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_x_y[i] = 4.0 * g_x_z_x_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_x_z[i] = 4.0 * g_x_z_x_z[i] * a_exp * b_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_z_0_0_0_0_y_x, g_x_z_0_0_0_0_y_y, g_x_z_0_0_0_0_y_z, g_x_z_y_x, g_x_z_y_y, g_x_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_0_y_x[i] = 4.0 * g_x_z_y_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_y_y[i] = 4.0 * g_x_z_y_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_y_z[i] = 4.0 * g_x_z_y_z[i] * a_exp * b_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_z_0_0_0_0_z_x, g_x_z_0_0_0_0_z_y, g_x_z_0_0_0_0_z_z, g_x_z_z_x, g_x_z_z_y, g_x_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_0_z_x[i] = 4.0 * g_x_z_z_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_z_y[i] = 4.0 * g_x_z_z_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_z_z[i] = 4.0 * g_x_z_z_z[i] * a_exp * b_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_y_x_0_0_0_0_x_x, g_y_x_0_0_0_0_x_y, g_y_x_0_0_0_0_x_z, g_y_x_x_x, g_y_x_x_y, g_y_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_0_x_x[i] = 4.0 * g_y_x_x_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_x_y[i] = 4.0 * g_y_x_x_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_x_z[i] = 4.0 * g_y_x_x_z[i] * a_exp * b_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_y_x_0_0_0_0_y_x, g_y_x_0_0_0_0_y_y, g_y_x_0_0_0_0_y_z, g_y_x_y_x, g_y_x_y_y, g_y_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_0_y_x[i] = 4.0 * g_y_x_y_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_y_y[i] = 4.0 * g_y_x_y_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_y_z[i] = 4.0 * g_y_x_y_z[i] * a_exp * b_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_y_x_0_0_0_0_z_x, g_y_x_0_0_0_0_z_y, g_y_x_0_0_0_0_z_z, g_y_x_z_x, g_y_x_z_y, g_y_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_0_z_x[i] = 4.0 * g_y_x_z_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_z_y[i] = 4.0 * g_y_x_z_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_z_z[i] = 4.0 * g_y_x_z_z[i] * a_exp * b_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_y_y_0_0_0_0_x_x, g_y_y_0_0_0_0_x_y, g_y_y_0_0_0_0_x_z, g_y_y_x_x, g_y_y_x_y, g_y_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_0_x_x[i] = 4.0 * g_y_y_x_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_x_y[i] = 4.0 * g_y_y_x_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_x_z[i] = 4.0 * g_y_y_x_z[i] * a_exp * b_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_y_y_0_0_0_0_y_x, g_y_y_0_0_0_0_y_y, g_y_y_0_0_0_0_y_z, g_y_y_y_x, g_y_y_y_y, g_y_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_0_y_x[i] = 4.0 * g_y_y_y_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_y_y[i] = 4.0 * g_y_y_y_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_y_z[i] = 4.0 * g_y_y_y_z[i] * a_exp * b_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_y_y_0_0_0_0_z_x, g_y_y_0_0_0_0_z_y, g_y_y_0_0_0_0_z_z, g_y_y_z_x, g_y_y_z_y, g_y_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_0_z_x[i] = 4.0 * g_y_y_z_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_z_y[i] = 4.0 * g_y_y_z_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_z_z[i] = 4.0 * g_y_y_z_z[i] * a_exp * b_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_y_z_0_0_0_0_x_x, g_y_z_0_0_0_0_x_y, g_y_z_0_0_0_0_x_z, g_y_z_x_x, g_y_z_x_y, g_y_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_0_x_x[i] = 4.0 * g_y_z_x_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_x_y[i] = 4.0 * g_y_z_x_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_x_z[i] = 4.0 * g_y_z_x_z[i] * a_exp * b_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_y_z_0_0_0_0_y_x, g_y_z_0_0_0_0_y_y, g_y_z_0_0_0_0_y_z, g_y_z_y_x, g_y_z_y_y, g_y_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_0_y_x[i] = 4.0 * g_y_z_y_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_y_y[i] = 4.0 * g_y_z_y_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_y_z[i] = 4.0 * g_y_z_y_z[i] * a_exp * b_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_y_z_0_0_0_0_z_x, g_y_z_0_0_0_0_z_y, g_y_z_0_0_0_0_z_z, g_y_z_z_x, g_y_z_z_y, g_y_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_0_z_x[i] = 4.0 * g_y_z_z_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_z_y[i] = 4.0 * g_y_z_z_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_z_z[i] = 4.0 * g_y_z_z_z[i] * a_exp * b_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_z_x_0_0_0_0_x_x, g_z_x_0_0_0_0_x_y, g_z_x_0_0_0_0_x_z, g_z_x_x_x, g_z_x_x_y, g_z_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_0_x_x[i] = 4.0 * g_z_x_x_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_x_y[i] = 4.0 * g_z_x_x_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_x_z[i] = 4.0 * g_z_x_x_z[i] * a_exp * b_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_z_x_0_0_0_0_y_x, g_z_x_0_0_0_0_y_y, g_z_x_0_0_0_0_y_z, g_z_x_y_x, g_z_x_y_y, g_z_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_0_y_x[i] = 4.0 * g_z_x_y_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_y_y[i] = 4.0 * g_z_x_y_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_y_z[i] = 4.0 * g_z_x_y_z[i] * a_exp * b_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_z_x_0_0_0_0_z_x, g_z_x_0_0_0_0_z_y, g_z_x_0_0_0_0_z_z, g_z_x_z_x, g_z_x_z_y, g_z_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_0_z_x[i] = 4.0 * g_z_x_z_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_z_y[i] = 4.0 * g_z_x_z_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_z_z[i] = 4.0 * g_z_x_z_z[i] * a_exp * b_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_z_y_0_0_0_0_x_x, g_z_y_0_0_0_0_x_y, g_z_y_0_0_0_0_x_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_0_x_x[i] = 4.0 * g_z_y_x_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_x_y[i] = 4.0 * g_z_y_x_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_x_z[i] = 4.0 * g_z_y_x_z[i] * a_exp * b_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_z_y_0_0_0_0_y_x, g_z_y_0_0_0_0_y_y, g_z_y_0_0_0_0_y_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_0_y_x[i] = 4.0 * g_z_y_y_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_y_y[i] = 4.0 * g_z_y_y_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_y_z[i] = 4.0 * g_z_y_y_z[i] * a_exp * b_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_z_y_0_0_0_0_z_x, g_z_y_0_0_0_0_z_y, g_z_y_0_0_0_0_z_z, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_0_z_x[i] = 4.0 * g_z_y_z_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_z_y[i] = 4.0 * g_z_y_z_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_z_z[i] = 4.0 * g_z_y_z_z[i] * a_exp * b_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_z_z_0_0_0_0_x_x, g_z_z_0_0_0_0_x_y, g_z_z_0_0_0_0_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_0_x_x[i] = 4.0 * g_z_z_x_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_x_y[i] = 4.0 * g_z_z_x_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_x_z[i] = 4.0 * g_z_z_x_z[i] * a_exp * b_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_z_z_0_0_0_0_y_x, g_z_z_0_0_0_0_y_y, g_z_z_0_0_0_0_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_0_y_x[i] = 4.0 * g_z_z_y_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_y_y[i] = 4.0 * g_z_z_y_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_y_z[i] = 4.0 * g_z_z_y_z[i] * a_exp * b_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_z_z_0_0_0_0_z_x, g_z_z_0_0_0_0_z_y, g_z_z_0_0_0_0_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_0_z_x[i] = 4.0 * g_z_z_z_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_z_y[i] = 4.0 * g_z_z_z_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_z_z[i] = 4.0 * g_z_z_z_z[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

