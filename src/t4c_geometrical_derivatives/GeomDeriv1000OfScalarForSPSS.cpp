#include "GeomDeriv1000OfScalarForSPSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_spss_0(CSimdArray<double>& buffer_1000_spss,
                     const CSimdArray<double>& buffer_ppss,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_spss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_ppss

    auto g_x_x_0_0 = buffer_ppss[0];

    auto g_x_y_0_0 = buffer_ppss[1];

    auto g_x_z_0_0 = buffer_ppss[2];

    auto g_y_x_0_0 = buffer_ppss[3];

    auto g_y_y_0_0 = buffer_ppss[4];

    auto g_y_z_0_0 = buffer_ppss[5];

    auto g_z_x_0_0 = buffer_ppss[6];

    auto g_z_y_0_0 = buffer_ppss[7];

    auto g_z_z_0_0 = buffer_ppss[8];

    /// Set up components of integrals buffer : buffer_1000_spss

    auto g_x_0_0_0_0_x_0_0 = buffer_1000_spss[0];

    auto g_x_0_0_0_0_y_0_0 = buffer_1000_spss[1];

    auto g_x_0_0_0_0_z_0_0 = buffer_1000_spss[2];

    auto g_y_0_0_0_0_x_0_0 = buffer_1000_spss[3];

    auto g_y_0_0_0_0_y_0_0 = buffer_1000_spss[4];

    auto g_y_0_0_0_0_z_0_0 = buffer_1000_spss[5];

    auto g_z_0_0_0_0_x_0_0 = buffer_1000_spss[6];

    auto g_z_0_0_0_0_y_0_0 = buffer_1000_spss[7];

    auto g_z_0_0_0_0_z_0_0 = buffer_1000_spss[8];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_0_0_0_x_0_0, g_x_0_0_0_0_y_0_0, g_x_0_0_0_0_z_0_0, g_x_x_0_0, g_x_y_0_0, g_x_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_x_0_0[i] = 2.0 * g_x_x_0_0[i] * a_exp;

        g_x_0_0_0_0_y_0_0[i] = 2.0 * g_x_y_0_0[i] * a_exp;

        g_x_0_0_0_0_z_0_0[i] = 2.0 * g_x_z_0_0[i] * a_exp;
    }
    // integrals block (3-6)

    #pragma omp simd aligned(g_y_0_0_0_0_x_0_0, g_y_0_0_0_0_y_0_0, g_y_0_0_0_0_z_0_0, g_y_x_0_0, g_y_y_0_0, g_y_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_y_0_0_0_0_x_0_0[i] = 2.0 * g_y_x_0_0[i] * a_exp;

        g_y_0_0_0_0_y_0_0[i] = 2.0 * g_y_y_0_0[i] * a_exp;

        g_y_0_0_0_0_z_0_0[i] = 2.0 * g_y_z_0_0[i] * a_exp;
    }
    // integrals block (6-9)

    #pragma omp simd aligned(g_z_0_0_0_0_x_0_0, g_z_0_0_0_0_y_0_0, g_z_0_0_0_0_z_0_0, g_z_x_0_0, g_z_y_0_0, g_z_z_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_z_0_0_0_0_x_0_0[i] = 2.0 * g_z_x_0_0[i] * a_exp;

        g_z_0_0_0_0_y_0_0[i] = 2.0 * g_z_y_0_0[i] * a_exp;

        g_z_0_0_0_0_z_0_0[i] = 2.0 * g_z_z_0_0[i] * a_exp;
    }
}

} // t4c_geom namespace

