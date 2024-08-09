#include "GeomDeriv1010OfScalarForSPSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_spsp_0(CSimdArray<double>& buffer_1010_spsp,
                     const CSimdArray<double>& buffer_pppp,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_spsp.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1010_spsp

    auto g_x_0_x_0_0_x_0_x = buffer_1010_spsp[0];

    auto g_x_0_x_0_0_x_0_y = buffer_1010_spsp[1];

    auto g_x_0_x_0_0_x_0_z = buffer_1010_spsp[2];

    auto g_x_0_x_0_0_y_0_x = buffer_1010_spsp[3];

    auto g_x_0_x_0_0_y_0_y = buffer_1010_spsp[4];

    auto g_x_0_x_0_0_y_0_z = buffer_1010_spsp[5];

    auto g_x_0_x_0_0_z_0_x = buffer_1010_spsp[6];

    auto g_x_0_x_0_0_z_0_y = buffer_1010_spsp[7];

    auto g_x_0_x_0_0_z_0_z = buffer_1010_spsp[8];

    auto g_x_0_y_0_0_x_0_x = buffer_1010_spsp[9];

    auto g_x_0_y_0_0_x_0_y = buffer_1010_spsp[10];

    auto g_x_0_y_0_0_x_0_z = buffer_1010_spsp[11];

    auto g_x_0_y_0_0_y_0_x = buffer_1010_spsp[12];

    auto g_x_0_y_0_0_y_0_y = buffer_1010_spsp[13];

    auto g_x_0_y_0_0_y_0_z = buffer_1010_spsp[14];

    auto g_x_0_y_0_0_z_0_x = buffer_1010_spsp[15];

    auto g_x_0_y_0_0_z_0_y = buffer_1010_spsp[16];

    auto g_x_0_y_0_0_z_0_z = buffer_1010_spsp[17];

    auto g_x_0_z_0_0_x_0_x = buffer_1010_spsp[18];

    auto g_x_0_z_0_0_x_0_y = buffer_1010_spsp[19];

    auto g_x_0_z_0_0_x_0_z = buffer_1010_spsp[20];

    auto g_x_0_z_0_0_y_0_x = buffer_1010_spsp[21];

    auto g_x_0_z_0_0_y_0_y = buffer_1010_spsp[22];

    auto g_x_0_z_0_0_y_0_z = buffer_1010_spsp[23];

    auto g_x_0_z_0_0_z_0_x = buffer_1010_spsp[24];

    auto g_x_0_z_0_0_z_0_y = buffer_1010_spsp[25];

    auto g_x_0_z_0_0_z_0_z = buffer_1010_spsp[26];

    auto g_y_0_x_0_0_x_0_x = buffer_1010_spsp[27];

    auto g_y_0_x_0_0_x_0_y = buffer_1010_spsp[28];

    auto g_y_0_x_0_0_x_0_z = buffer_1010_spsp[29];

    auto g_y_0_x_0_0_y_0_x = buffer_1010_spsp[30];

    auto g_y_0_x_0_0_y_0_y = buffer_1010_spsp[31];

    auto g_y_0_x_0_0_y_0_z = buffer_1010_spsp[32];

    auto g_y_0_x_0_0_z_0_x = buffer_1010_spsp[33];

    auto g_y_0_x_0_0_z_0_y = buffer_1010_spsp[34];

    auto g_y_0_x_0_0_z_0_z = buffer_1010_spsp[35];

    auto g_y_0_y_0_0_x_0_x = buffer_1010_spsp[36];

    auto g_y_0_y_0_0_x_0_y = buffer_1010_spsp[37];

    auto g_y_0_y_0_0_x_0_z = buffer_1010_spsp[38];

    auto g_y_0_y_0_0_y_0_x = buffer_1010_spsp[39];

    auto g_y_0_y_0_0_y_0_y = buffer_1010_spsp[40];

    auto g_y_0_y_0_0_y_0_z = buffer_1010_spsp[41];

    auto g_y_0_y_0_0_z_0_x = buffer_1010_spsp[42];

    auto g_y_0_y_0_0_z_0_y = buffer_1010_spsp[43];

    auto g_y_0_y_0_0_z_0_z = buffer_1010_spsp[44];

    auto g_y_0_z_0_0_x_0_x = buffer_1010_spsp[45];

    auto g_y_0_z_0_0_x_0_y = buffer_1010_spsp[46];

    auto g_y_0_z_0_0_x_0_z = buffer_1010_spsp[47];

    auto g_y_0_z_0_0_y_0_x = buffer_1010_spsp[48];

    auto g_y_0_z_0_0_y_0_y = buffer_1010_spsp[49];

    auto g_y_0_z_0_0_y_0_z = buffer_1010_spsp[50];

    auto g_y_0_z_0_0_z_0_x = buffer_1010_spsp[51];

    auto g_y_0_z_0_0_z_0_y = buffer_1010_spsp[52];

    auto g_y_0_z_0_0_z_0_z = buffer_1010_spsp[53];

    auto g_z_0_x_0_0_x_0_x = buffer_1010_spsp[54];

    auto g_z_0_x_0_0_x_0_y = buffer_1010_spsp[55];

    auto g_z_0_x_0_0_x_0_z = buffer_1010_spsp[56];

    auto g_z_0_x_0_0_y_0_x = buffer_1010_spsp[57];

    auto g_z_0_x_0_0_y_0_y = buffer_1010_spsp[58];

    auto g_z_0_x_0_0_y_0_z = buffer_1010_spsp[59];

    auto g_z_0_x_0_0_z_0_x = buffer_1010_spsp[60];

    auto g_z_0_x_0_0_z_0_y = buffer_1010_spsp[61];

    auto g_z_0_x_0_0_z_0_z = buffer_1010_spsp[62];

    auto g_z_0_y_0_0_x_0_x = buffer_1010_spsp[63];

    auto g_z_0_y_0_0_x_0_y = buffer_1010_spsp[64];

    auto g_z_0_y_0_0_x_0_z = buffer_1010_spsp[65];

    auto g_z_0_y_0_0_y_0_x = buffer_1010_spsp[66];

    auto g_z_0_y_0_0_y_0_y = buffer_1010_spsp[67];

    auto g_z_0_y_0_0_y_0_z = buffer_1010_spsp[68];

    auto g_z_0_y_0_0_z_0_x = buffer_1010_spsp[69];

    auto g_z_0_y_0_0_z_0_y = buffer_1010_spsp[70];

    auto g_z_0_y_0_0_z_0_z = buffer_1010_spsp[71];

    auto g_z_0_z_0_0_x_0_x = buffer_1010_spsp[72];

    auto g_z_0_z_0_0_x_0_y = buffer_1010_spsp[73];

    auto g_z_0_z_0_0_x_0_z = buffer_1010_spsp[74];

    auto g_z_0_z_0_0_y_0_x = buffer_1010_spsp[75];

    auto g_z_0_z_0_0_y_0_y = buffer_1010_spsp[76];

    auto g_z_0_z_0_0_y_0_z = buffer_1010_spsp[77];

    auto g_z_0_z_0_0_z_0_x = buffer_1010_spsp[78];

    auto g_z_0_z_0_0_z_0_y = buffer_1010_spsp[79];

    auto g_z_0_z_0_0_z_0_z = buffer_1010_spsp[80];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_x_0_0_x_0_x, g_x_0_x_0_0_x_0_y, g_x_0_x_0_0_x_0_z, g_x_x_x_x, g_x_x_x_y, g_x_x_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_x_0_x[i] = 4.0 * g_x_x_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_0_y[i] = 4.0 * g_x_x_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_x_0_z[i] = 4.0 * g_x_x_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_x_0_0_y_0_x, g_x_0_x_0_0_y_0_y, g_x_0_x_0_0_y_0_z, g_x_y_x_x, g_x_y_x_y, g_x_y_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_y_0_x[i] = 4.0 * g_x_y_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_0_y[i] = 4.0 * g_x_y_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_y_0_z[i] = 4.0 * g_x_y_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_x_0_0_z_0_x, g_x_0_x_0_0_z_0_y, g_x_0_x_0_0_z_0_z, g_x_z_x_x, g_x_z_x_y, g_x_z_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_z_0_x[i] = 4.0 * g_x_z_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_0_y[i] = 4.0 * g_x_z_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_z_0_z[i] = 4.0 * g_x_z_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_x_0_y_0_0_x_0_x, g_x_0_y_0_0_x_0_y, g_x_0_y_0_0_x_0_z, g_x_x_y_x, g_x_x_y_y, g_x_x_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_x_0_x[i] = 4.0 * g_x_x_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_0_y[i] = 4.0 * g_x_x_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_x_0_z[i] = 4.0 * g_x_x_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_x_0_y_0_0_y_0_x, g_x_0_y_0_0_y_0_y, g_x_0_y_0_0_y_0_z, g_x_y_y_x, g_x_y_y_y, g_x_y_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_y_0_x[i] = 4.0 * g_x_y_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_0_y[i] = 4.0 * g_x_y_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_y_0_z[i] = 4.0 * g_x_y_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_x_0_y_0_0_z_0_x, g_x_0_y_0_0_z_0_y, g_x_0_y_0_0_z_0_z, g_x_z_y_x, g_x_z_y_y, g_x_z_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_z_0_x[i] = 4.0 * g_x_z_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_0_y[i] = 4.0 * g_x_z_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_z_0_z[i] = 4.0 * g_x_z_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_x_0_z_0_0_x_0_x, g_x_0_z_0_0_x_0_y, g_x_0_z_0_0_x_0_z, g_x_x_z_x, g_x_x_z_y, g_x_x_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_x_0_x[i] = 4.0 * g_x_x_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_0_y[i] = 4.0 * g_x_x_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_x_0_z[i] = 4.0 * g_x_x_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_x_0_z_0_0_y_0_x, g_x_0_z_0_0_y_0_y, g_x_0_z_0_0_y_0_z, g_x_y_z_x, g_x_y_z_y, g_x_y_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_y_0_x[i] = 4.0 * g_x_y_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_0_y[i] = 4.0 * g_x_y_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_y_0_z[i] = 4.0 * g_x_y_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_x_0_z_0_0_z_0_x, g_x_0_z_0_0_z_0_y, g_x_0_z_0_0_z_0_z, g_x_z_z_x, g_x_z_z_y, g_x_z_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_z_0_x[i] = 4.0 * g_x_z_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_0_y[i] = 4.0 * g_x_z_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_z_0_z[i] = 4.0 * g_x_z_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (27-30)

    #pragma omp simd aligned(g_y_0_x_0_0_x_0_x, g_y_0_x_0_0_x_0_y, g_y_0_x_0_0_x_0_z, g_y_x_x_x, g_y_x_x_y, g_y_x_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_x_0_x[i] = 4.0 * g_y_x_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_0_y[i] = 4.0 * g_y_x_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_x_0_z[i] = 4.0 * g_y_x_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (30-33)

    #pragma omp simd aligned(g_y_0_x_0_0_y_0_x, g_y_0_x_0_0_y_0_y, g_y_0_x_0_0_y_0_z, g_y_y_x_x, g_y_y_x_y, g_y_y_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_y_0_x[i] = 4.0 * g_y_y_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_0_y[i] = 4.0 * g_y_y_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_y_0_z[i] = 4.0 * g_y_y_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (33-36)

    #pragma omp simd aligned(g_y_0_x_0_0_z_0_x, g_y_0_x_0_0_z_0_y, g_y_0_x_0_0_z_0_z, g_y_z_x_x, g_y_z_x_y, g_y_z_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_z_0_x[i] = 4.0 * g_y_z_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_0_y[i] = 4.0 * g_y_z_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_z_0_z[i] = 4.0 * g_y_z_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (36-39)

    #pragma omp simd aligned(g_y_0_y_0_0_x_0_x, g_y_0_y_0_0_x_0_y, g_y_0_y_0_0_x_0_z, g_y_x_y_x, g_y_x_y_y, g_y_x_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_x_0_x[i] = 4.0 * g_y_x_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_0_y[i] = 4.0 * g_y_x_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_x_0_z[i] = 4.0 * g_y_x_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (39-42)

    #pragma omp simd aligned(g_y_0_y_0_0_y_0_x, g_y_0_y_0_0_y_0_y, g_y_0_y_0_0_y_0_z, g_y_y_y_x, g_y_y_y_y, g_y_y_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_y_0_x[i] = 4.0 * g_y_y_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_0_y[i] = 4.0 * g_y_y_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_y_0_z[i] = 4.0 * g_y_y_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (42-45)

    #pragma omp simd aligned(g_y_0_y_0_0_z_0_x, g_y_0_y_0_0_z_0_y, g_y_0_y_0_0_z_0_z, g_y_z_y_x, g_y_z_y_y, g_y_z_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_z_0_x[i] = 4.0 * g_y_z_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_0_y[i] = 4.0 * g_y_z_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_z_0_z[i] = 4.0 * g_y_z_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (45-48)

    #pragma omp simd aligned(g_y_0_z_0_0_x_0_x, g_y_0_z_0_0_x_0_y, g_y_0_z_0_0_x_0_z, g_y_x_z_x, g_y_x_z_y, g_y_x_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_x_0_x[i] = 4.0 * g_y_x_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_0_y[i] = 4.0 * g_y_x_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_x_0_z[i] = 4.0 * g_y_x_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (48-51)

    #pragma omp simd aligned(g_y_0_z_0_0_y_0_x, g_y_0_z_0_0_y_0_y, g_y_0_z_0_0_y_0_z, g_y_y_z_x, g_y_y_z_y, g_y_y_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_y_0_x[i] = 4.0 * g_y_y_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_0_y[i] = 4.0 * g_y_y_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_y_0_z[i] = 4.0 * g_y_y_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (51-54)

    #pragma omp simd aligned(g_y_0_z_0_0_z_0_x, g_y_0_z_0_0_z_0_y, g_y_0_z_0_0_z_0_z, g_y_z_z_x, g_y_z_z_y, g_y_z_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_z_0_x[i] = 4.0 * g_y_z_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_0_y[i] = 4.0 * g_y_z_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_z_0_z[i] = 4.0 * g_y_z_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (54-57)

    #pragma omp simd aligned(g_z_0_x_0_0_x_0_x, g_z_0_x_0_0_x_0_y, g_z_0_x_0_0_x_0_z, g_z_x_x_x, g_z_x_x_y, g_z_x_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_x_0_x[i] = 4.0 * g_z_x_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_0_y[i] = 4.0 * g_z_x_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_x_0_z[i] = 4.0 * g_z_x_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (57-60)

    #pragma omp simd aligned(g_z_0_x_0_0_y_0_x, g_z_0_x_0_0_y_0_y, g_z_0_x_0_0_y_0_z, g_z_y_x_x, g_z_y_x_y, g_z_y_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_y_0_x[i] = 4.0 * g_z_y_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_0_y[i] = 4.0 * g_z_y_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_y_0_z[i] = 4.0 * g_z_y_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (60-63)

    #pragma omp simd aligned(g_z_0_x_0_0_z_0_x, g_z_0_x_0_0_z_0_y, g_z_0_x_0_0_z_0_z, g_z_z_x_x, g_z_z_x_y, g_z_z_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_z_0_x[i] = 4.0 * g_z_z_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_0_y[i] = 4.0 * g_z_z_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_z_0_z[i] = 4.0 * g_z_z_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (63-66)

    #pragma omp simd aligned(g_z_0_y_0_0_x_0_x, g_z_0_y_0_0_x_0_y, g_z_0_y_0_0_x_0_z, g_z_x_y_x, g_z_x_y_y, g_z_x_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_x_0_x[i] = 4.0 * g_z_x_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_0_y[i] = 4.0 * g_z_x_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_x_0_z[i] = 4.0 * g_z_x_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (66-69)

    #pragma omp simd aligned(g_z_0_y_0_0_y_0_x, g_z_0_y_0_0_y_0_y, g_z_0_y_0_0_y_0_z, g_z_y_y_x, g_z_y_y_y, g_z_y_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_y_0_x[i] = 4.0 * g_z_y_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_0_y[i] = 4.0 * g_z_y_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_y_0_z[i] = 4.0 * g_z_y_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (69-72)

    #pragma omp simd aligned(g_z_0_y_0_0_z_0_x, g_z_0_y_0_0_z_0_y, g_z_0_y_0_0_z_0_z, g_z_z_y_x, g_z_z_y_y, g_z_z_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_z_0_x[i] = 4.0 * g_z_z_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_0_y[i] = 4.0 * g_z_z_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_z_0_z[i] = 4.0 * g_z_z_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (72-75)

    #pragma omp simd aligned(g_z_0_z_0_0_x_0_x, g_z_0_z_0_0_x_0_y, g_z_0_z_0_0_x_0_z, g_z_x_z_x, g_z_x_z_y, g_z_x_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_x_0_x[i] = 4.0 * g_z_x_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_0_y[i] = 4.0 * g_z_x_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_x_0_z[i] = 4.0 * g_z_x_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (75-78)

    #pragma omp simd aligned(g_z_0_z_0_0_y_0_x, g_z_0_z_0_0_y_0_y, g_z_0_z_0_0_y_0_z, g_z_y_z_x, g_z_y_z_y, g_z_y_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_y_0_x[i] = 4.0 * g_z_y_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_0_y[i] = 4.0 * g_z_y_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_y_0_z[i] = 4.0 * g_z_y_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (78-81)

    #pragma omp simd aligned(g_z_0_z_0_0_z_0_x, g_z_0_z_0_0_z_0_y, g_z_0_z_0_0_z_0_z, g_z_z_z_x, g_z_z_z_y, g_z_z_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_z_0_x[i] = 4.0 * g_z_z_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_0_y[i] = 4.0 * g_z_z_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_z_0_z[i] = 4.0 * g_z_z_z_z[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace

