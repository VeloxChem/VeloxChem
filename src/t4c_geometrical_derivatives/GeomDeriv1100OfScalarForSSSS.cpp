#include "GeomDeriv1100OfScalarForSSSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1100_ssss_0(CSimdArray<double>& buffer_1100_ssss,
                     const CSimdArray<double>& buffer_ppss,
                     const double a_exp,
                     const double b_exp) -> void
{
    const auto ndims = buffer_1100_ssss.number_of_columns();

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

    /// Set up components of integrals buffer : buffer_1100_ssss

    auto g_x_x_0_0_0_0_0_0 = buffer_1100_ssss[0];

    auto g_x_y_0_0_0_0_0_0 = buffer_1100_ssss[1];

    auto g_x_z_0_0_0_0_0_0 = buffer_1100_ssss[2];

    auto g_y_x_0_0_0_0_0_0 = buffer_1100_ssss[3];

    auto g_y_y_0_0_0_0_0_0 = buffer_1100_ssss[4];

    auto g_y_z_0_0_0_0_0_0 = buffer_1100_ssss[5];

    auto g_z_x_0_0_0_0_0_0 = buffer_1100_ssss[6];

    auto g_z_y_0_0_0_0_0_0 = buffer_1100_ssss[7];

    auto g_z_z_0_0_0_0_0_0 = buffer_1100_ssss[8];

    // integrals block (0-9)

    #pragma omp simd aligned(g_x_x_0_0, g_x_x_0_0_0_0_0_0, g_x_y_0_0, g_x_y_0_0_0_0_0_0, g_x_z_0_0, g_x_z_0_0_0_0_0_0, g_y_x_0_0, g_y_x_0_0_0_0_0_0, g_y_y_0_0, g_y_y_0_0_0_0_0_0, g_y_z_0_0, g_y_z_0_0_0_0_0_0, g_z_x_0_0, g_z_x_0_0_0_0_0_0, g_z_y_0_0, g_z_y_0_0_0_0_0_0, g_z_z_0_0, g_z_z_0_0_0_0_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_x_0_0_0_0_0_0[i] = 4.0 * g_x_x_0_0[i] * a_exp * b_exp;

        g_x_y_0_0_0_0_0_0[i] = 4.0 * g_x_y_0_0[i] * a_exp * b_exp;

        g_x_z_0_0_0_0_0_0[i] = 4.0 * g_x_z_0_0[i] * a_exp * b_exp;

        g_y_x_0_0_0_0_0_0[i] = 4.0 * g_y_x_0_0[i] * a_exp * b_exp;

        g_y_y_0_0_0_0_0_0[i] = 4.0 * g_y_y_0_0[i] * a_exp * b_exp;

        g_y_z_0_0_0_0_0_0[i] = 4.0 * g_y_z_0_0[i] * a_exp * b_exp;

        g_z_x_0_0_0_0_0_0[i] = 4.0 * g_z_x_0_0[i] * a_exp * b_exp;

        g_z_y_0_0_0_0_0_0[i] = 4.0 * g_z_y_0_0[i] * a_exp * b_exp;

        g_z_z_0_0_0_0_0_0[i] = 4.0 * g_z_z_0_0[i] * a_exp * b_exp;
    }
}

} // t4c_geom namespace

