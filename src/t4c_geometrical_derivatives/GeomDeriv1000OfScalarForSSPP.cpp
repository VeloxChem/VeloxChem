#include "GeomDeriv1000OfScalarForSSPP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_sspp_0(CSimdArray<double>& buffer_1000_sspp,
                     const CSimdArray<double>& buffer_pspp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_sspp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_pspp

    auto g_x_0_x_x = buffer_pspp[0];

    auto g_x_0_x_y = buffer_pspp[1];

    auto g_x_0_x_z = buffer_pspp[2];

    auto g_x_0_y_x = buffer_pspp[3];

    auto g_x_0_y_y = buffer_pspp[4];

    auto g_x_0_y_z = buffer_pspp[5];

    auto g_x_0_z_x = buffer_pspp[6];

    auto g_x_0_z_y = buffer_pspp[7];

    auto g_x_0_z_z = buffer_pspp[8];

    auto g_y_0_x_x = buffer_pspp[9];

    auto g_y_0_x_y = buffer_pspp[10];

    auto g_y_0_x_z = buffer_pspp[11];

    auto g_y_0_y_x = buffer_pspp[12];

    auto g_y_0_y_y = buffer_pspp[13];

    auto g_y_0_y_z = buffer_pspp[14];

    auto g_y_0_z_x = buffer_pspp[15];

    auto g_y_0_z_y = buffer_pspp[16];

    auto g_y_0_z_z = buffer_pspp[17];

    auto g_z_0_x_x = buffer_pspp[18];

    auto g_z_0_x_y = buffer_pspp[19];

    auto g_z_0_x_z = buffer_pspp[20];

    auto g_z_0_y_x = buffer_pspp[21];

    auto g_z_0_y_y = buffer_pspp[22];

    auto g_z_0_y_z = buffer_pspp[23];

    auto g_z_0_z_x = buffer_pspp[24];

    auto g_z_0_z_y = buffer_pspp[25];

    auto g_z_0_z_z = buffer_pspp[26];

    /// Set up components of integrals buffer : buffer_1000_sspp

    auto g_x_0_0_0_0_0_x_x = buffer_1000_sspp[0];

    auto g_x_0_0_0_0_0_x_y = buffer_1000_sspp[1];

    auto g_x_0_0_0_0_0_x_z = buffer_1000_sspp[2];

    auto g_x_0_0_0_0_0_y_x = buffer_1000_sspp[3];

    auto g_x_0_0_0_0_0_y_y = buffer_1000_sspp[4];

    auto g_x_0_0_0_0_0_y_z = buffer_1000_sspp[5];

    auto g_x_0_0_0_0_0_z_x = buffer_1000_sspp[6];

    auto g_x_0_0_0_0_0_z_y = buffer_1000_sspp[7];

    auto g_x_0_0_0_0_0_z_z = buffer_1000_sspp[8];

    auto g_y_0_0_0_0_0_x_x = buffer_1000_sspp[9];

    auto g_y_0_0_0_0_0_x_y = buffer_1000_sspp[10];

    auto g_y_0_0_0_0_0_x_z = buffer_1000_sspp[11];

    auto g_y_0_0_0_0_0_y_x = buffer_1000_sspp[12];

    auto g_y_0_0_0_0_0_y_y = buffer_1000_sspp[13];

    auto g_y_0_0_0_0_0_y_z = buffer_1000_sspp[14];

    auto g_y_0_0_0_0_0_z_x = buffer_1000_sspp[15];

    auto g_y_0_0_0_0_0_z_y = buffer_1000_sspp[16];

    auto g_y_0_0_0_0_0_z_z = buffer_1000_sspp[17];

    auto g_z_0_0_0_0_0_x_x = buffer_1000_sspp[18];

    auto g_z_0_0_0_0_0_x_y = buffer_1000_sspp[19];

    auto g_z_0_0_0_0_0_x_z = buffer_1000_sspp[20];

    auto g_z_0_0_0_0_0_y_x = buffer_1000_sspp[21];

    auto g_z_0_0_0_0_0_y_y = buffer_1000_sspp[22];

    auto g_z_0_0_0_0_0_y_z = buffer_1000_sspp[23];

    auto g_z_0_0_0_0_0_z_x = buffer_1000_sspp[24];

    auto g_z_0_0_0_0_0_z_y = buffer_1000_sspp[25];

    auto g_z_0_0_0_0_0_z_z = buffer_1000_sspp[26];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_0_0_0_0_x_x, g_x_0_0_0_0_0_x_y, g_x_0_0_0_0_0_x_z, g_x_0_x_x, g_x_0_x_y, g_x_0_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_0_x_x[i] = 2.0 * g_x_0_x_x[i] * a_exp;

        g_x_0_0_0_0_0_x_y[i] = 2.0 * g_x_0_x_y[i] * a_exp;

        g_x_0_0_0_0_0_x_z[i] = 2.0 * g_x_0_x_z[i] * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_0_0_0_0_y_x, g_x_0_0_0_0_0_y_y, g_x_0_0_0_0_0_y_z, g_x_0_y_x, g_x_0_y_y, g_x_0_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_0_y_x[i] = 2.0 * g_x_0_y_x[i] * a_exp;

        g_x_0_0_0_0_0_y_y[i] = 2.0 * g_x_0_y_y[i] * a_exp;

        g_x_0_0_0_0_0_y_z[i] = 2.0 * g_x_0_y_z[i] * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_0_0_0_0_z_x, g_x_0_0_0_0_0_z_y, g_x_0_0_0_0_0_z_z, g_x_0_z_x, g_x_0_z_y, g_x_0_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_0_z_x[i] = 2.0 * g_x_0_z_x[i] * a_exp;

        g_x_0_0_0_0_0_z_y[i] = 2.0 * g_x_0_z_y[i] * a_exp;

        g_x_0_0_0_0_0_z_z[i] = 2.0 * g_x_0_z_z[i] * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_y_0_0_0_0_0_x_x, g_y_0_0_0_0_0_x_y, g_y_0_0_0_0_0_x_z, g_y_0_x_x, g_y_0_x_y, g_y_0_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_0_x_x[i] = 2.0 * g_y_0_x_x[i] * a_exp;

        g_y_0_0_0_0_0_x_y[i] = 2.0 * g_y_0_x_y[i] * a_exp;

        g_y_0_0_0_0_0_x_z[i] = 2.0 * g_y_0_x_z[i] * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_y_0_0_0_0_0_y_x, g_y_0_0_0_0_0_y_y, g_y_0_0_0_0_0_y_z, g_y_0_y_x, g_y_0_y_y, g_y_0_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_0_y_x[i] = 2.0 * g_y_0_y_x[i] * a_exp;

        g_y_0_0_0_0_0_y_y[i] = 2.0 * g_y_0_y_y[i] * a_exp;

        g_y_0_0_0_0_0_y_z[i] = 2.0 * g_y_0_y_z[i] * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_y_0_0_0_0_0_z_x, g_y_0_0_0_0_0_z_y, g_y_0_0_0_0_0_z_z, g_y_0_z_x, g_y_0_z_y, g_y_0_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_0_z_x[i] = 2.0 * g_y_0_z_x[i] * a_exp;

        g_y_0_0_0_0_0_z_y[i] = 2.0 * g_y_0_z_y[i] * a_exp;

        g_y_0_0_0_0_0_z_z[i] = 2.0 * g_y_0_z_z[i] * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_z_0_0_0_0_0_x_x, g_z_0_0_0_0_0_x_y, g_z_0_0_0_0_0_x_z, g_z_0_x_x, g_z_0_x_y, g_z_0_x_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_0_x_x[i] = 2.0 * g_z_0_x_x[i] * a_exp;

        g_z_0_0_0_0_0_x_y[i] = 2.0 * g_z_0_x_y[i] * a_exp;

        g_z_0_0_0_0_0_x_z[i] = 2.0 * g_z_0_x_z[i] * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_z_0_0_0_0_0_y_x, g_z_0_0_0_0_0_y_y, g_z_0_0_0_0_0_y_z, g_z_0_y_x, g_z_0_y_y, g_z_0_y_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_0_y_x[i] = 2.0 * g_z_0_y_x[i] * a_exp;

        g_z_0_0_0_0_0_y_y[i] = 2.0 * g_z_0_y_y[i] * a_exp;

        g_z_0_0_0_0_0_y_z[i] = 2.0 * g_z_0_y_z[i] * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_z_0_0_0_0_0_z_x, g_z_0_0_0_0_0_z_y, g_z_0_0_0_0_0_z_z, g_z_0_z_x, g_z_0_z_y, g_z_0_z_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_0_z_x[i] = 2.0 * g_z_0_z_x[i] * a_exp;

        g_z_0_0_0_0_0_z_y[i] = 2.0 * g_z_0_z_y[i] * a_exp;

        g_z_0_0_0_0_0_z_z[i] = 2.0 * g_z_0_z_z[i] * a_exp;
    }
}

} // t4c_geom namespace

