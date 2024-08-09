#include "GeomDeriv1000OfScalarForSPPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_sppp_0(CSimdArray<double>& buffer_1000_sppp,
                     const CSimdArray<double>& buffer_pppp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_sppp.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_sppp

    auto g_x_0_0_0_0_x_x_x = buffer_1000_sppp[0];

    auto g_x_0_0_0_0_x_x_y = buffer_1000_sppp[1];

    auto g_x_0_0_0_0_x_x_z = buffer_1000_sppp[2];

    auto g_x_0_0_0_0_x_y_x = buffer_1000_sppp[3];

    auto g_x_0_0_0_0_x_y_y = buffer_1000_sppp[4];

    auto g_x_0_0_0_0_x_y_z = buffer_1000_sppp[5];

    auto g_x_0_0_0_0_x_z_x = buffer_1000_sppp[6];

    auto g_x_0_0_0_0_x_z_y = buffer_1000_sppp[7];

    auto g_x_0_0_0_0_x_z_z = buffer_1000_sppp[8];

    auto g_x_0_0_0_0_y_x_x = buffer_1000_sppp[9];

    auto g_x_0_0_0_0_y_x_y = buffer_1000_sppp[10];

    auto g_x_0_0_0_0_y_x_z = buffer_1000_sppp[11];

    auto g_x_0_0_0_0_y_y_x = buffer_1000_sppp[12];

    auto g_x_0_0_0_0_y_y_y = buffer_1000_sppp[13];

    auto g_x_0_0_0_0_y_y_z = buffer_1000_sppp[14];

    auto g_x_0_0_0_0_y_z_x = buffer_1000_sppp[15];

    auto g_x_0_0_0_0_y_z_y = buffer_1000_sppp[16];

    auto g_x_0_0_0_0_y_z_z = buffer_1000_sppp[17];

    auto g_x_0_0_0_0_z_x_x = buffer_1000_sppp[18];

    auto g_x_0_0_0_0_z_x_y = buffer_1000_sppp[19];

    auto g_x_0_0_0_0_z_x_z = buffer_1000_sppp[20];

    auto g_x_0_0_0_0_z_y_x = buffer_1000_sppp[21];

    auto g_x_0_0_0_0_z_y_y = buffer_1000_sppp[22];

    auto g_x_0_0_0_0_z_y_z = buffer_1000_sppp[23];

    auto g_x_0_0_0_0_z_z_x = buffer_1000_sppp[24];

    auto g_x_0_0_0_0_z_z_y = buffer_1000_sppp[25];

    auto g_x_0_0_0_0_z_z_z = buffer_1000_sppp[26];

    auto g_y_0_0_0_0_x_x_x = buffer_1000_sppp[27];

    auto g_y_0_0_0_0_x_x_y = buffer_1000_sppp[28];

    auto g_y_0_0_0_0_x_x_z = buffer_1000_sppp[29];

    auto g_y_0_0_0_0_x_y_x = buffer_1000_sppp[30];

    auto g_y_0_0_0_0_x_y_y = buffer_1000_sppp[31];

    auto g_y_0_0_0_0_x_y_z = buffer_1000_sppp[32];

    auto g_y_0_0_0_0_x_z_x = buffer_1000_sppp[33];

    auto g_y_0_0_0_0_x_z_y = buffer_1000_sppp[34];

    auto g_y_0_0_0_0_x_z_z = buffer_1000_sppp[35];

    auto g_y_0_0_0_0_y_x_x = buffer_1000_sppp[36];

    auto g_y_0_0_0_0_y_x_y = buffer_1000_sppp[37];

    auto g_y_0_0_0_0_y_x_z = buffer_1000_sppp[38];

    auto g_y_0_0_0_0_y_y_x = buffer_1000_sppp[39];

    auto g_y_0_0_0_0_y_y_y = buffer_1000_sppp[40];

    auto g_y_0_0_0_0_y_y_z = buffer_1000_sppp[41];

    auto g_y_0_0_0_0_y_z_x = buffer_1000_sppp[42];

    auto g_y_0_0_0_0_y_z_y = buffer_1000_sppp[43];

    auto g_y_0_0_0_0_y_z_z = buffer_1000_sppp[44];

    auto g_y_0_0_0_0_z_x_x = buffer_1000_sppp[45];

    auto g_y_0_0_0_0_z_x_y = buffer_1000_sppp[46];

    auto g_y_0_0_0_0_z_x_z = buffer_1000_sppp[47];

    auto g_y_0_0_0_0_z_y_x = buffer_1000_sppp[48];

    auto g_y_0_0_0_0_z_y_y = buffer_1000_sppp[49];

    auto g_y_0_0_0_0_z_y_z = buffer_1000_sppp[50];

    auto g_y_0_0_0_0_z_z_x = buffer_1000_sppp[51];

    auto g_y_0_0_0_0_z_z_y = buffer_1000_sppp[52];

    auto g_y_0_0_0_0_z_z_z = buffer_1000_sppp[53];

    auto g_z_0_0_0_0_x_x_x = buffer_1000_sppp[54];

    auto g_z_0_0_0_0_x_x_y = buffer_1000_sppp[55];

    auto g_z_0_0_0_0_x_x_z = buffer_1000_sppp[56];

    auto g_z_0_0_0_0_x_y_x = buffer_1000_sppp[57];

    auto g_z_0_0_0_0_x_y_y = buffer_1000_sppp[58];

    auto g_z_0_0_0_0_x_y_z = buffer_1000_sppp[59];

    auto g_z_0_0_0_0_x_z_x = buffer_1000_sppp[60];

    auto g_z_0_0_0_0_x_z_y = buffer_1000_sppp[61];

    auto g_z_0_0_0_0_x_z_z = buffer_1000_sppp[62];

    auto g_z_0_0_0_0_y_x_x = buffer_1000_sppp[63];

    auto g_z_0_0_0_0_y_x_y = buffer_1000_sppp[64];

    auto g_z_0_0_0_0_y_x_z = buffer_1000_sppp[65];

    auto g_z_0_0_0_0_y_y_x = buffer_1000_sppp[66];

    auto g_z_0_0_0_0_y_y_y = buffer_1000_sppp[67];

    auto g_z_0_0_0_0_y_y_z = buffer_1000_sppp[68];

    auto g_z_0_0_0_0_y_z_x = buffer_1000_sppp[69];

    auto g_z_0_0_0_0_y_z_y = buffer_1000_sppp[70];

    auto g_z_0_0_0_0_y_z_z = buffer_1000_sppp[71];

    auto g_z_0_0_0_0_z_x_x = buffer_1000_sppp[72];

    auto g_z_0_0_0_0_z_x_y = buffer_1000_sppp[73];

    auto g_z_0_0_0_0_z_x_z = buffer_1000_sppp[74];

    auto g_z_0_0_0_0_z_y_x = buffer_1000_sppp[75];

    auto g_z_0_0_0_0_z_y_y = buffer_1000_sppp[76];

    auto g_z_0_0_0_0_z_y_z = buffer_1000_sppp[77];

    auto g_z_0_0_0_0_z_z_x = buffer_1000_sppp[78];

    auto g_z_0_0_0_0_z_z_y = buffer_1000_sppp[79];

    auto g_z_0_0_0_0_z_z_z = buffer_1000_sppp[80];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_0_0_0_x_x_x, g_x_0_0_0_0_x_x_y, g_x_0_0_0_0_x_x_z, g_x_x_x_x, g_x_x_x_y, g_x_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_x_x_x[i] = 2.0 * g_x_x_x_x[i] * a_exp;

        g_x_0_0_0_0_x_x_y[i] = 2.0 * g_x_x_x_y[i] * a_exp;

        g_x_0_0_0_0_x_x_z[i] = 2.0 * g_x_x_x_z[i] * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_0_0_0_x_y_x, g_x_0_0_0_0_x_y_y, g_x_0_0_0_0_x_y_z, g_x_x_y_x, g_x_x_y_y, g_x_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_x_y_x[i] = 2.0 * g_x_x_y_x[i] * a_exp;

        g_x_0_0_0_0_x_y_y[i] = 2.0 * g_x_x_y_y[i] * a_exp;

        g_x_0_0_0_0_x_y_z[i] = 2.0 * g_x_x_y_z[i] * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_0_0_0_x_z_x, g_x_0_0_0_0_x_z_y, g_x_0_0_0_0_x_z_z, g_x_x_z_x, g_x_x_z_y, g_x_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_x_z_x[i] = 2.0 * g_x_x_z_x[i] * a_exp;

        g_x_0_0_0_0_x_z_y[i] = 2.0 * g_x_x_z_y[i] * a_exp;

        g_x_0_0_0_0_x_z_z[i] = 2.0 * g_x_x_z_z[i] * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_0_0_0_0_y_x_x, g_x_0_0_0_0_y_x_y, g_x_0_0_0_0_y_x_z, g_x_y_x_x, g_x_y_x_y, g_x_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_y_x_x[i] = 2.0 * g_x_y_x_x[i] * a_exp;

        g_x_0_0_0_0_y_x_y[i] = 2.0 * g_x_y_x_y[i] * a_exp;

        g_x_0_0_0_0_y_x_z[i] = 2.0 * g_x_y_x_z[i] * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_0_0_0_0_y_y_x, g_x_0_0_0_0_y_y_y, g_x_0_0_0_0_y_y_z, g_x_y_y_x, g_x_y_y_y, g_x_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_y_y_x[i] = 2.0 * g_x_y_y_x[i] * a_exp;

        g_x_0_0_0_0_y_y_y[i] = 2.0 * g_x_y_y_y[i] * a_exp;

        g_x_0_0_0_0_y_y_z[i] = 2.0 * g_x_y_y_z[i] * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_0_0_0_0_y_z_x, g_x_0_0_0_0_y_z_y, g_x_0_0_0_0_y_z_z, g_x_y_z_x, g_x_y_z_y, g_x_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_y_z_x[i] = 2.0 * g_x_y_z_x[i] * a_exp;

        g_x_0_0_0_0_y_z_y[i] = 2.0 * g_x_y_z_y[i] * a_exp;

        g_x_0_0_0_0_y_z_z[i] = 2.0 * g_x_y_z_z[i] * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_0_0_0_0_z_x_x, g_x_0_0_0_0_z_x_y, g_x_0_0_0_0_z_x_z, g_x_z_x_x, g_x_z_x_y, g_x_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_z_x_x[i] = 2.0 * g_x_z_x_x[i] * a_exp;

        g_x_0_0_0_0_z_x_y[i] = 2.0 * g_x_z_x_y[i] * a_exp;

        g_x_0_0_0_0_z_x_z[i] = 2.0 * g_x_z_x_z[i] * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_0_0_0_0_z_y_x, g_x_0_0_0_0_z_y_y, g_x_0_0_0_0_z_y_z, g_x_z_y_x, g_x_z_y_y, g_x_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_z_y_x[i] = 2.0 * g_x_z_y_x[i] * a_exp;

        g_x_0_0_0_0_z_y_y[i] = 2.0 * g_x_z_y_y[i] * a_exp;

        g_x_0_0_0_0_z_y_z[i] = 2.0 * g_x_z_y_z[i] * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_0_0_0_0_z_z_x, g_x_0_0_0_0_z_z_y, g_x_0_0_0_0_z_z_z, g_x_z_z_x, g_x_z_z_y, g_x_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_z_z_x[i] = 2.0 * g_x_z_z_x[i] * a_exp;

        g_x_0_0_0_0_z_z_y[i] = 2.0 * g_x_z_z_y[i] * a_exp;

        g_x_0_0_0_0_z_z_z[i] = 2.0 * g_x_z_z_z[i] * a_exp;
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_y_0_0_0_0_x_x_x, g_y_0_0_0_0_x_x_y, g_y_0_0_0_0_x_x_z, g_y_x_x_x, g_y_x_x_y, g_y_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_x_x_x[i] = 2.0 * g_y_x_x_x[i] * a_exp;

        g_y_0_0_0_0_x_x_y[i] = 2.0 * g_y_x_x_y[i] * a_exp;

        g_y_0_0_0_0_x_x_z[i] = 2.0 * g_y_x_x_z[i] * a_exp;
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_y_0_0_0_0_x_y_x, g_y_0_0_0_0_x_y_y, g_y_0_0_0_0_x_y_z, g_y_x_y_x, g_y_x_y_y, g_y_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_x_y_x[i] = 2.0 * g_y_x_y_x[i] * a_exp;

        g_y_0_0_0_0_x_y_y[i] = 2.0 * g_y_x_y_y[i] * a_exp;

        g_y_0_0_0_0_x_y_z[i] = 2.0 * g_y_x_y_z[i] * a_exp;
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_y_0_0_0_0_x_z_x, g_y_0_0_0_0_x_z_y, g_y_0_0_0_0_x_z_z, g_y_x_z_x, g_y_x_z_y, g_y_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_x_z_x[i] = 2.0 * g_y_x_z_x[i] * a_exp;

        g_y_0_0_0_0_x_z_y[i] = 2.0 * g_y_x_z_y[i] * a_exp;

        g_y_0_0_0_0_x_z_z[i] = 2.0 * g_y_x_z_z[i] * a_exp;
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_y_0_0_0_0_y_x_x, g_y_0_0_0_0_y_x_y, g_y_0_0_0_0_y_x_z, g_y_y_x_x, g_y_y_x_y, g_y_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_y_x_x[i] = 2.0 * g_y_y_x_x[i] * a_exp;

        g_y_0_0_0_0_y_x_y[i] = 2.0 * g_y_y_x_y[i] * a_exp;

        g_y_0_0_0_0_y_x_z[i] = 2.0 * g_y_y_x_z[i] * a_exp;
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_y_0_0_0_0_y_y_x, g_y_0_0_0_0_y_y_y, g_y_0_0_0_0_y_y_z, g_y_y_y_x, g_y_y_y_y, g_y_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_y_y_x[i] = 2.0 * g_y_y_y_x[i] * a_exp;

        g_y_0_0_0_0_y_y_y[i] = 2.0 * g_y_y_y_y[i] * a_exp;

        g_y_0_0_0_0_y_y_z[i] = 2.0 * g_y_y_y_z[i] * a_exp;
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_y_0_0_0_0_y_z_x, g_y_0_0_0_0_y_z_y, g_y_0_0_0_0_y_z_z, g_y_y_z_x, g_y_y_z_y, g_y_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_y_z_x[i] = 2.0 * g_y_y_z_x[i] * a_exp;

        g_y_0_0_0_0_y_z_y[i] = 2.0 * g_y_y_z_y[i] * a_exp;

        g_y_0_0_0_0_y_z_z[i] = 2.0 * g_y_y_z_z[i] * a_exp;
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_y_0_0_0_0_z_x_x, g_y_0_0_0_0_z_x_y, g_y_0_0_0_0_z_x_z, g_y_z_x_x, g_y_z_x_y, g_y_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_z_x_x[i] = 2.0 * g_y_z_x_x[i] * a_exp;

        g_y_0_0_0_0_z_x_y[i] = 2.0 * g_y_z_x_y[i] * a_exp;

        g_y_0_0_0_0_z_x_z[i] = 2.0 * g_y_z_x_z[i] * a_exp;
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_y_0_0_0_0_z_y_x, g_y_0_0_0_0_z_y_y, g_y_0_0_0_0_z_y_z, g_y_z_y_x, g_y_z_y_y, g_y_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_z_y_x[i] = 2.0 * g_y_z_y_x[i] * a_exp;

        g_y_0_0_0_0_z_y_y[i] = 2.0 * g_y_z_y_y[i] * a_exp;

        g_y_0_0_0_0_z_y_z[i] = 2.0 * g_y_z_y_z[i] * a_exp;
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_y_0_0_0_0_z_z_x, g_y_0_0_0_0_z_z_y, g_y_0_0_0_0_z_z_z, g_y_z_z_x, g_y_z_z_y, g_y_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_z_z_x[i] = 2.0 * g_y_z_z_x[i] * a_exp;

        g_y_0_0_0_0_z_z_y[i] = 2.0 * g_y_z_z_y[i] * a_exp;

        g_y_0_0_0_0_z_z_z[i] = 2.0 * g_y_z_z_z[i] * a_exp;
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_z_0_0_0_0_x_x_x, g_z_0_0_0_0_x_x_y, g_z_0_0_0_0_x_x_z, g_z_x_x_x, g_z_x_x_y, g_z_x_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_x_x_x[i] = 2.0 * g_z_x_x_x[i] * a_exp;

        g_z_0_0_0_0_x_x_y[i] = 2.0 * g_z_x_x_y[i] * a_exp;

        g_z_0_0_0_0_x_x_z[i] = 2.0 * g_z_x_x_z[i] * a_exp;
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_z_0_0_0_0_x_y_x, g_z_0_0_0_0_x_y_y, g_z_0_0_0_0_x_y_z, g_z_x_y_x, g_z_x_y_y, g_z_x_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_x_y_x[i] = 2.0 * g_z_x_y_x[i] * a_exp;

        g_z_0_0_0_0_x_y_y[i] = 2.0 * g_z_x_y_y[i] * a_exp;

        g_z_0_0_0_0_x_y_z[i] = 2.0 * g_z_x_y_z[i] * a_exp;
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_z_0_0_0_0_x_z_x, g_z_0_0_0_0_x_z_y, g_z_0_0_0_0_x_z_z, g_z_x_z_x, g_z_x_z_y, g_z_x_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_x_z_x[i] = 2.0 * g_z_x_z_x[i] * a_exp;

        g_z_0_0_0_0_x_z_y[i] = 2.0 * g_z_x_z_y[i] * a_exp;

        g_z_0_0_0_0_x_z_z[i] = 2.0 * g_z_x_z_z[i] * a_exp;
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_z_0_0_0_0_y_x_x, g_z_0_0_0_0_y_x_y, g_z_0_0_0_0_y_x_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_y_x_x[i] = 2.0 * g_z_y_x_x[i] * a_exp;

        g_z_0_0_0_0_y_x_y[i] = 2.0 * g_z_y_x_y[i] * a_exp;

        g_z_0_0_0_0_y_x_z[i] = 2.0 * g_z_y_x_z[i] * a_exp;
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_z_0_0_0_0_y_y_x, g_z_0_0_0_0_y_y_y, g_z_0_0_0_0_y_y_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_y_y_x[i] = 2.0 * g_z_y_y_x[i] * a_exp;

        g_z_0_0_0_0_y_y_y[i] = 2.0 * g_z_y_y_y[i] * a_exp;

        g_z_0_0_0_0_y_y_z[i] = 2.0 * g_z_y_y_z[i] * a_exp;
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_z_0_0_0_0_y_z_x, g_z_0_0_0_0_y_z_y, g_z_0_0_0_0_y_z_z, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_y_z_x[i] = 2.0 * g_z_y_z_x[i] * a_exp;

        g_z_0_0_0_0_y_z_y[i] = 2.0 * g_z_y_z_y[i] * a_exp;

        g_z_0_0_0_0_y_z_z[i] = 2.0 * g_z_y_z_z[i] * a_exp;
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_z_0_0_0_0_z_x_x, g_z_0_0_0_0_z_x_y, g_z_0_0_0_0_z_x_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_z_x_x[i] = 2.0 * g_z_z_x_x[i] * a_exp;

        g_z_0_0_0_0_z_x_y[i] = 2.0 * g_z_z_x_y[i] * a_exp;

        g_z_0_0_0_0_z_x_z[i] = 2.0 * g_z_z_x_z[i] * a_exp;
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_z_0_0_0_0_z_y_x, g_z_0_0_0_0_z_y_y, g_z_0_0_0_0_z_y_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_z_y_x[i] = 2.0 * g_z_z_y_x[i] * a_exp;

        g_z_0_0_0_0_z_y_y[i] = 2.0 * g_z_z_y_y[i] * a_exp;

        g_z_0_0_0_0_z_y_z[i] = 2.0 * g_z_z_y_z[i] * a_exp;
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_z_0_0_0_0_z_z_x, g_z_0_0_0_0_z_z_y, g_z_0_0_0_0_z_z_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_z_z_x[i] = 2.0 * g_z_z_z_x[i] * a_exp;

        g_z_0_0_0_0_z_z_y[i] = 2.0 * g_z_z_z_y[i] * a_exp;

        g_z_0_0_0_0_z_z_z[i] = 2.0 * g_z_z_z_z[i] * a_exp;
    }
}

} // t4c_geom namespace

