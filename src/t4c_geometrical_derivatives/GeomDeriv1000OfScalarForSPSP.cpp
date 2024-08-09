#include "GeomDeriv1000OfScalarForSPSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_spsp_0(CSimdArray<double>& buffer_1000_spsp,
                     const CSimdArray<double>& buffer_ppsp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_spsp.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_spsp

    auto g_x_0_0_0_0_x_0_x = buffer_1000_spsp[0];

    auto g_x_0_0_0_0_x_0_y = buffer_1000_spsp[1];

    auto g_x_0_0_0_0_x_0_z = buffer_1000_spsp[2];

    auto g_x_0_0_0_0_y_0_x = buffer_1000_spsp[3];

    auto g_x_0_0_0_0_y_0_y = buffer_1000_spsp[4];

    auto g_x_0_0_0_0_y_0_z = buffer_1000_spsp[5];

    auto g_x_0_0_0_0_z_0_x = buffer_1000_spsp[6];

    auto g_x_0_0_0_0_z_0_y = buffer_1000_spsp[7];

    auto g_x_0_0_0_0_z_0_z = buffer_1000_spsp[8];

    auto g_y_0_0_0_0_x_0_x = buffer_1000_spsp[9];

    auto g_y_0_0_0_0_x_0_y = buffer_1000_spsp[10];

    auto g_y_0_0_0_0_x_0_z = buffer_1000_spsp[11];

    auto g_y_0_0_0_0_y_0_x = buffer_1000_spsp[12];

    auto g_y_0_0_0_0_y_0_y = buffer_1000_spsp[13];

    auto g_y_0_0_0_0_y_0_z = buffer_1000_spsp[14];

    auto g_y_0_0_0_0_z_0_x = buffer_1000_spsp[15];

    auto g_y_0_0_0_0_z_0_y = buffer_1000_spsp[16];

    auto g_y_0_0_0_0_z_0_z = buffer_1000_spsp[17];

    auto g_z_0_0_0_0_x_0_x = buffer_1000_spsp[18];

    auto g_z_0_0_0_0_x_0_y = buffer_1000_spsp[19];

    auto g_z_0_0_0_0_x_0_z = buffer_1000_spsp[20];

    auto g_z_0_0_0_0_y_0_x = buffer_1000_spsp[21];

    auto g_z_0_0_0_0_y_0_y = buffer_1000_spsp[22];

    auto g_z_0_0_0_0_y_0_z = buffer_1000_spsp[23];

    auto g_z_0_0_0_0_z_0_x = buffer_1000_spsp[24];

    auto g_z_0_0_0_0_z_0_y = buffer_1000_spsp[25];

    auto g_z_0_0_0_0_z_0_z = buffer_1000_spsp[26];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_0_0_0_x_0_x, g_x_0_0_0_0_x_0_y, g_x_0_0_0_0_x_0_z, g_x_x_0_x, g_x_x_0_y, g_x_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_x_0_x[i] = 2.0 * g_x_x_0_x[i] * a_exp;

        g_x_0_0_0_0_x_0_y[i] = 2.0 * g_x_x_0_y[i] * a_exp;

        g_x_0_0_0_0_x_0_z[i] = 2.0 * g_x_x_0_z[i] * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_0_0_0_y_0_x, g_x_0_0_0_0_y_0_y, g_x_0_0_0_0_y_0_z, g_x_y_0_x, g_x_y_0_y, g_x_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_y_0_x[i] = 2.0 * g_x_y_0_x[i] * a_exp;

        g_x_0_0_0_0_y_0_y[i] = 2.0 * g_x_y_0_y[i] * a_exp;

        g_x_0_0_0_0_y_0_z[i] = 2.0 * g_x_y_0_z[i] * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_0_0_0_z_0_x, g_x_0_0_0_0_z_0_y, g_x_0_0_0_0_z_0_z, g_x_z_0_x, g_x_z_0_y, g_x_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_z_0_x[i] = 2.0 * g_x_z_0_x[i] * a_exp;

        g_x_0_0_0_0_z_0_y[i] = 2.0 * g_x_z_0_y[i] * a_exp;

        g_x_0_0_0_0_z_0_z[i] = 2.0 * g_x_z_0_z[i] * a_exp;
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_y_0_0_0_0_x_0_x, g_y_0_0_0_0_x_0_y, g_y_0_0_0_0_x_0_z, g_y_x_0_x, g_y_x_0_y, g_y_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_x_0_x[i] = 2.0 * g_y_x_0_x[i] * a_exp;

        g_y_0_0_0_0_x_0_y[i] = 2.0 * g_y_x_0_y[i] * a_exp;

        g_y_0_0_0_0_x_0_z[i] = 2.0 * g_y_x_0_z[i] * a_exp;
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_y_0_0_0_0_y_0_x, g_y_0_0_0_0_y_0_y, g_y_0_0_0_0_y_0_z, g_y_y_0_x, g_y_y_0_y, g_y_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_y_0_x[i] = 2.0 * g_y_y_0_x[i] * a_exp;

        g_y_0_0_0_0_y_0_y[i] = 2.0 * g_y_y_0_y[i] * a_exp;

        g_y_0_0_0_0_y_0_z[i] = 2.0 * g_y_y_0_z[i] * a_exp;
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_y_0_0_0_0_z_0_x, g_y_0_0_0_0_z_0_y, g_y_0_0_0_0_z_0_z, g_y_z_0_x, g_y_z_0_y, g_y_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_z_0_x[i] = 2.0 * g_y_z_0_x[i] * a_exp;

        g_y_0_0_0_0_z_0_y[i] = 2.0 * g_y_z_0_y[i] * a_exp;

        g_y_0_0_0_0_z_0_z[i] = 2.0 * g_y_z_0_z[i] * a_exp;
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_z_0_0_0_0_x_0_x, g_z_0_0_0_0_x_0_y, g_z_0_0_0_0_x_0_z, g_z_x_0_x, g_z_x_0_y, g_z_x_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_x_0_x[i] = 2.0 * g_z_x_0_x[i] * a_exp;

        g_z_0_0_0_0_x_0_y[i] = 2.0 * g_z_x_0_y[i] * a_exp;

        g_z_0_0_0_0_x_0_z[i] = 2.0 * g_z_x_0_z[i] * a_exp;
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_z_0_0_0_0_y_0_x, g_z_0_0_0_0_y_0_y, g_z_0_0_0_0_y_0_z, g_z_y_0_x, g_z_y_0_y, g_z_y_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_y_0_x[i] = 2.0 * g_z_y_0_x[i] * a_exp;

        g_z_0_0_0_0_y_0_y[i] = 2.0 * g_z_y_0_y[i] * a_exp;

        g_z_0_0_0_0_y_0_z[i] = 2.0 * g_z_y_0_z[i] * a_exp;
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_z_0_0_0_0_z_0_x, g_z_0_0_0_0_z_0_y, g_z_0_0_0_0_z_0_z, g_z_z_0_x, g_z_z_0_y, g_z_z_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_z_0_x[i] = 2.0 * g_z_z_0_x[i] * a_exp;

        g_z_0_0_0_0_z_0_y[i] = 2.0 * g_z_z_0_y[i] * a_exp;

        g_z_0_0_0_0_z_0_z[i] = 2.0 * g_z_z_0_z[i] * a_exp;
    }
}

} // t4c_geom namespace

