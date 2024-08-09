#include "GeomDeriv1010OfScalarForSSSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_sssp_0(CSimdArray<double>& buffer_1010_sssp,
                     const CSimdArray<double>& buffer_pspp,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_sssp.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1010_sssp

    auto g_x_0_x_0_0_0_0_x = buffer_1010_sssp[0];

    auto g_x_0_x_0_0_0_0_y = buffer_1010_sssp[1];

    auto g_x_0_x_0_0_0_0_z = buffer_1010_sssp[2];

    auto g_x_0_y_0_0_0_0_x = buffer_1010_sssp[3];

    auto g_x_0_y_0_0_0_0_y = buffer_1010_sssp[4];

    auto g_x_0_y_0_0_0_0_z = buffer_1010_sssp[5];

    auto g_x_0_z_0_0_0_0_x = buffer_1010_sssp[6];

    auto g_x_0_z_0_0_0_0_y = buffer_1010_sssp[7];

    auto g_x_0_z_0_0_0_0_z = buffer_1010_sssp[8];

    auto g_y_0_x_0_0_0_0_x = buffer_1010_sssp[9];

    auto g_y_0_x_0_0_0_0_y = buffer_1010_sssp[10];

    auto g_y_0_x_0_0_0_0_z = buffer_1010_sssp[11];

    auto g_y_0_y_0_0_0_0_x = buffer_1010_sssp[12];

    auto g_y_0_y_0_0_0_0_y = buffer_1010_sssp[13];

    auto g_y_0_y_0_0_0_0_z = buffer_1010_sssp[14];

    auto g_y_0_z_0_0_0_0_x = buffer_1010_sssp[15];

    auto g_y_0_z_0_0_0_0_y = buffer_1010_sssp[16];

    auto g_y_0_z_0_0_0_0_z = buffer_1010_sssp[17];

    auto g_z_0_x_0_0_0_0_x = buffer_1010_sssp[18];

    auto g_z_0_x_0_0_0_0_y = buffer_1010_sssp[19];

    auto g_z_0_x_0_0_0_0_z = buffer_1010_sssp[20];

    auto g_z_0_y_0_0_0_0_x = buffer_1010_sssp[21];

    auto g_z_0_y_0_0_0_0_y = buffer_1010_sssp[22];

    auto g_z_0_y_0_0_0_0_z = buffer_1010_sssp[23];

    auto g_z_0_z_0_0_0_0_x = buffer_1010_sssp[24];

    auto g_z_0_z_0_0_0_0_y = buffer_1010_sssp[25];

    auto g_z_0_z_0_0_0_0_z = buffer_1010_sssp[26];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_x_0_0_0_0_x, g_x_0_x_0_0_0_0_y, g_x_0_x_0_0_0_0_z, g_x_0_x_x, g_x_0_x_y, g_x_0_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_0_0_x[i] = 4.0 * g_x_0_x_x[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_0_y[i] = 4.0 * g_x_0_x_y[i] * a_exp * c_exps[i];

        g_x_0_x_0_0_0_0_z[i] = 4.0 * g_x_0_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_x_0_y_0_0_0_0_x, g_x_0_y_0_0_0_0_y, g_x_0_y_0_0_0_0_z, g_x_0_y_x, g_x_0_y_y, g_x_0_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_y_0_0_0_0_x[i] = 4.0 * g_x_0_y_x[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_0_y[i] = 4.0 * g_x_0_y_y[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_0_z[i] = 4.0 * g_x_0_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_x_0_z_0_0_0_0_x, g_x_0_z_0_0_0_0_y, g_x_0_z_0_0_0_0_z, g_x_0_z_x, g_x_0_z_y, g_x_0_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_z_0_0_0_0_x[i] = 4.0 * g_x_0_z_x[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_0_y[i] = 4.0 * g_x_0_z_y[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_0_z[i] = 4.0 * g_x_0_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (9-12)

    #pragma omp simd aligned(g_y_0_x_0_0_0_0_x, g_y_0_x_0_0_0_0_y, g_y_0_x_0_0_0_0_z, g_y_0_x_x, g_y_0_x_y, g_y_0_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_x_0_0_0_0_x[i] = 4.0 * g_y_0_x_x[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_0_y[i] = 4.0 * g_y_0_x_y[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_0_z[i] = 4.0 * g_y_0_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (12-15)

    #pragma omp simd aligned(g_y_0_y_0_0_0_0_x, g_y_0_y_0_0_0_0_y, g_y_0_y_0_0_0_0_z, g_y_0_y_x, g_y_0_y_y, g_y_0_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_y_0_0_0_0_x[i] = 4.0 * g_y_0_y_x[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_0_y[i] = 4.0 * g_y_0_y_y[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_0_z[i] = 4.0 * g_y_0_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (15-18)

    #pragma omp simd aligned(g_y_0_z_0_0_0_0_x, g_y_0_z_0_0_0_0_y, g_y_0_z_0_0_0_0_z, g_y_0_z_x, g_y_0_z_y, g_y_0_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_z_0_0_0_0_x[i] = 4.0 * g_y_0_z_x[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_0_y[i] = 4.0 * g_y_0_z_y[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_0_z[i] = 4.0 * g_y_0_z_z[i] * a_exp * c_exps[i];
    }
    // integrals block (18-21)

    #pragma omp simd aligned(g_z_0_x_0_0_0_0_x, g_z_0_x_0_0_0_0_y, g_z_0_x_0_0_0_0_z, g_z_0_x_x, g_z_0_x_y, g_z_0_x_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_x_0_0_0_0_x[i] = 4.0 * g_z_0_x_x[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_0_y[i] = 4.0 * g_z_0_x_y[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_0_z[i] = 4.0 * g_z_0_x_z[i] * a_exp * c_exps[i];
    }
    // integrals block (21-24)

    #pragma omp simd aligned(g_z_0_y_0_0_0_0_x, g_z_0_y_0_0_0_0_y, g_z_0_y_0_0_0_0_z, g_z_0_y_x, g_z_0_y_y, g_z_0_y_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_y_0_0_0_0_x[i] = 4.0 * g_z_0_y_x[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_0_y[i] = 4.0 * g_z_0_y_y[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_0_z[i] = 4.0 * g_z_0_y_z[i] * a_exp * c_exps[i];
    }
    // integrals block (24-27)

    #pragma omp simd aligned(g_z_0_z_0_0_0_0_x, g_z_0_z_0_0_0_0_y, g_z_0_z_0_0_0_0_z, g_z_0_z_x, g_z_0_z_y, g_z_0_z_z, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_z_0_0_0_0_x[i] = 4.0 * g_z_0_z_x[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_0_y[i] = 4.0 * g_z_0_z_y[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_0_z[i] = 4.0 * g_z_0_z_z[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace

