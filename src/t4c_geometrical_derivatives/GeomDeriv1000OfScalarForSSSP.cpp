#include "GeomDeriv1000OfScalarForSSSP.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_sssp_0(CSimdArray<double>& buffer_1000_sssp,
                     const CSimdArray<double>& buffer_pssp,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_sssp.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1000_sssp

    auto g_x_0_0_0_0_0_0_x = buffer_1000_sssp[0];

    auto g_x_0_0_0_0_0_0_y = buffer_1000_sssp[1];

    auto g_x_0_0_0_0_0_0_z = buffer_1000_sssp[2];

    auto g_y_0_0_0_0_0_0_x = buffer_1000_sssp[3];

    auto g_y_0_0_0_0_0_0_y = buffer_1000_sssp[4];

    auto g_y_0_0_0_0_0_0_z = buffer_1000_sssp[5];

    auto g_z_0_0_0_0_0_0_x = buffer_1000_sssp[6];

    auto g_z_0_0_0_0_0_0_y = buffer_1000_sssp[7];

    auto g_z_0_0_0_0_0_0_z = buffer_1000_sssp[8];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_0_0_0_0_0_x, g_x_0_0_0_0_0_0_y, g_x_0_0_0_0_0_0_z, g_x_0_0_x, g_x_0_0_y, g_x_0_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_0_0_x[i] = 2.0 * g_x_0_0_x[i] * a_exp;

        g_x_0_0_0_0_0_0_y[i] = 2.0 * g_x_0_0_y[i] * a_exp;

        g_x_0_0_0_0_0_0_z[i] = 2.0 * g_x_0_0_z[i] * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_y_0_0_0_0_0_0_x, g_y_0_0_0_0_0_0_y, g_y_0_0_0_0_0_0_z, g_y_0_0_x, g_y_0_0_y, g_y_0_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_0_0_x[i] = 2.0 * g_y_0_0_x[i] * a_exp;

        g_y_0_0_0_0_0_0_y[i] = 2.0 * g_y_0_0_y[i] * a_exp;

        g_y_0_0_0_0_0_0_z[i] = 2.0 * g_y_0_0_z[i] * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_z_0_0_0_0_0_0_x, g_z_0_0_0_0_0_0_y, g_z_0_0_0_0_0_0_z, g_z_0_0_x, g_z_0_0_y, g_z_0_0_z  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_0_0_x[i] = 2.0 * g_z_0_0_x[i] * a_exp;

        g_z_0_0_0_0_0_0_y[i] = 2.0 * g_z_0_0_y[i] * a_exp;

        g_z_0_0_0_0_0_0_z[i] = 2.0 * g_z_0_0_z[i] * a_exp;
    }
}

} // t4c_geom namespace

