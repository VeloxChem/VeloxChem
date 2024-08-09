#include "GeomDeriv1100OfScalarForSSSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_sssp_0(CSimdArray<double>& buffer_1100_sssp,
                     const CSimdArray<double>& buffer_ppsp,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_sssp.number_of_columns();

    /// Set up components of auxilary buffer : buffer_ppsp

    auto g_x_x_0_x = buffer_ppsp[0];

    auto g_x_x_0_y = buffer_ppsp[1];

    auto g_x_x_0_z = buffer_ppsp[2];

    auto g_x_y_0_x = buffer_ppsp[3];

    auto g_x_y_0_y = buffer_ppsp[4];

    auto g_x_y_0_z = buffer_ppsp[5];

    auto g_x_z_0_x = buffer_ppsp[6];

    auto g_x_z_0_y = buffer_ppsp[7];

    auto g_x_z_0_z = buffer_ppsp[8];

    auto g_y_x_0_x = buffer_ppsp[9];

    auto g_y_x_0_y = buffer_ppsp[10];

    auto g_y_x_0_z = buffer_ppsp[11];

    auto g_y_y_0_x = buffer_ppsp[12];

    auto g_y_y_0_y = buffer_ppsp[13];

    auto g_y_y_0_z = buffer_ppsp[14];

    auto g_y_z_0_x = buffer_ppsp[15];

    auto g_y_z_0_y = buffer_ppsp[16];

    auto g_y_z_0_z = buffer_ppsp[17];

    auto g_z_x_0_x = buffer_ppsp[18];

    auto g_z_x_0_y = buffer_ppsp[19];

    auto g_z_x_0_z = buffer_ppsp[20];

    auto g_z_y_0_x = buffer_ppsp[21];

    auto g_z_y_0_y = buffer_ppsp[22];

    auto g_z_y_0_z = buffer_ppsp[23];

    auto g_z_z_0_x = buffer_ppsp[24];

    auto g_z_z_0_y = buffer_ppsp[25];

    auto g_z_z_0_z = buffer_ppsp[26];

    /// Set up components of integrals buffer : buffer_1100_sssp

    auto g_x_x_0_0_0_0_0_x = buffer_1100_sssp[0];

    auto g_x_x_0_0_0_0_0_y = buffer_1100_sssp[1];

    auto g_x_x_0_0_0_0_0_z = buffer_1100_sssp[2];

    auto g_x_y_0_0_0_0_0_x = buffer_1100_sssp[3];

    auto g_x_y_0_0_0_0_0_y = buffer_1100_sssp[4];

    auto g_x_y_0_0_0_0_0_z = buffer_1100_sssp[5];

    auto g_x_z_0_0_0_0_0_x = buffer_1100_sssp[6];

    auto g_x_z_0_0_0_0_0_y = buffer_1100_sssp[7];

    auto g_x_z_0_0_0_0_0_z = buffer_1100_sssp[8];

    auto g_y_x_0_0_0_0_0_x = buffer_1100_sssp[9];

    auto g_y_x_0_0_0_0_0_y = buffer_1100_sssp[10];

    auto g_y_x_0_0_0_0_0_z = buffer_1100_sssp[11];

    auto g_y_y_0_0_0_0_0_x = buffer_1100_sssp[12];

    auto g_y_y_0_0_0_0_0_y = buffer_1100_sssp[13];

    auto g_y_y_0_0_0_0_0_z = buffer_1100_sssp[14];

    auto g_y_z_0_0_0_0_0_x = buffer_1100_sssp[15];

    auto g_y_z_0_0_0_0_0_y = buffer_1100_sssp[16];

    auto g_y_z_0_0_0_0_0_z = buffer_1100_sssp[17];

    auto g_z_x_0_0_0_0_0_x = buffer_1100_sssp[18];

    auto g_z_x_0_0_0_0_0_y = buffer_1100_sssp[19];

    auto g_z_x_0_0_0_0_0_z = buffer_1100_sssp[20];

    auto g_z_y_0_0_0_0_0_x = buffer_1100_sssp[21];

    auto g_z_y_0_0_0_0_0_y = buffer_1100_sssp[22];

    auto g_z_y_0_0_0_0_0_z = buffer_1100_sssp[23];

    auto g_z_z_0_0_0_0_0_x = buffer_1100_sssp[24];

    auto g_z_z_0_0_0_0_0_y = buffer_1100_sssp[25];

    auto g_z_z_0_0_0_0_0_z = buffer_1100_sssp[26];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_x_0_0_0_0_0_x, g_x_x_0_0_0_0_0_y, g_x_x_0_0_0_0_0_z, g_x_x_0_x, g_x_x_0_y, g_x_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_0_0_x[i] = 4.0 * g_x_x_0_x[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_0_y[i] = 4.0 * g_x_x_0_y[i] * a_exp * b_exp;

        g_x_x_0_0_0_0_0_z[i] = 4.0 * g_x_x_0_z[i] * a_exp * b_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_y_0_0_0_0_0_x, g_x_y_0_0_0_0_0_y, g_x_y_0_0_0_0_0_z, g_x_y_0_x, g_x_y_0_y, g_x_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_y_0_0_0_0_0_x[i] = 4.0 * g_x_y_0_x[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_0_y[i] = 4.0 * g_x_y_0_y[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_0_z[i] = 4.0 * g_x_y_0_z[i] * a_exp * b_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_z_0_0_0_0_0_x, g_x_z_0_0_0_0_0_y, g_x_z_0_0_0_0_0_z, g_x_z_0_x, g_x_z_0_y, g_x_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_z_0_0_0_0_0_x[i] = 4.0 * g_x_z_0_x[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_0_y[i] = 4.0 * g_x_z_0_y[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_0_z[i] = 4.0 * g_x_z_0_z[i] * a_exp * b_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_y_x_0_0_0_0_0_x, g_y_x_0_0_0_0_0_y, g_y_x_0_0_0_0_0_z, g_y_x_0_x, g_y_x_0_y, g_y_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_x_0_0_0_0_0_x[i] = 4.0 * g_y_x_0_x[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_0_y[i] = 4.0 * g_y_x_0_y[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_0_z[i] = 4.0 * g_y_x_0_z[i] * a_exp * b_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_y_y_0_0_0_0_0_x, g_y_y_0_0_0_0_0_y, g_y_y_0_0_0_0_0_z, g_y_y_0_x, g_y_y_0_y, g_y_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_y_0_0_0_0_0_x[i] = 4.0 * g_y_y_0_x[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_0_y[i] = 4.0 * g_y_y_0_y[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_0_z[i] = 4.0 * g_y_y_0_z[i] * a_exp * b_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_y_z_0_0_0_0_0_x, g_y_z_0_0_0_0_0_y, g_y_z_0_0_0_0_0_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_z_0_0_0_0_0_x[i] = 4.0 * g_y_z_0_x[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_0_y[i] = 4.0 * g_y_z_0_y[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_0_z[i] = 4.0 * g_y_z_0_z[i] * a_exp * b_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_z_x_0_0_0_0_0_x, g_z_x_0_0_0_0_0_y, g_z_x_0_0_0_0_0_z, g_z_x_0_x, g_z_x_0_y, g_z_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_x_0_0_0_0_0_x[i] = 4.0 * g_z_x_0_x[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_0_y[i] = 4.0 * g_z_x_0_y[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_0_z[i] = 4.0 * g_z_x_0_z[i] * a_exp * b_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_z_y_0_0_0_0_0_x, g_z_y_0_0_0_0_0_y, g_z_y_0_0_0_0_0_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_y_0_0_0_0_0_x[i] = 4.0 * g_z_y_0_x[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_0_y[i] = 4.0 * g_z_y_0_y[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_0_z[i] = 4.0 * g_z_y_0_z[i] * a_exp * b_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_z_z_0_0_0_0_0_x, g_z_z_0_0_0_0_0_y, g_z_z_0_0_0_0_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_z_0_0_0_0_0_x[i] = 4.0 * g_z_z_0_x[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_0_y[i] = 4.0 * g_z_z_0_y[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_0_z[i] = 4.0 * g_z_z_0_z[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

