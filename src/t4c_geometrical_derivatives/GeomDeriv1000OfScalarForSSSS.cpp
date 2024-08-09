#include "GeomDeriv1000OfScalarForSSSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1000_ssss_0(CSimdArray<double>& buffer_1000_ssss,
                     const CSimdArray<double>& buffer_psss,
                     const double a_exp) -> void
{
    const auto ndims = buffer_1000_ssss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_psss

    auto g_x_0_0_0 = buffer_psss[0];

    auto g_y_0_0_0 = buffer_psss[1];

    auto g_z_0_0_0 = buffer_psss[2];

    /// Set up components of integrals buffer : buffer_1000_ssss

    auto g_x_0_0_0_0_0_0_0 = buffer_1000_ssss[0];

    auto g_y_0_0_0_0_0_0_0 = buffer_1000_ssss[1];

    auto g_z_0_0_0_0_0_0_0 = buffer_1000_ssss[2];

    // integrals block (0-3)

    #pragma omp simd aligned(g_x_0_0_0, g_x_0_0_0_0_0_0_0, g_y_0_0_0, g_y_0_0_0_0_0_0_0, g_z_0_0_0, g_z_0_0_0_0_0_0_0  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_0_0_0_0_0_0[i] = 2.0 * g_x_0_0_0[i] * a_exp;

        g_y_0_0_0_0_0_0_0[i] = 2.0 * g_y_0_0_0[i] * a_exp;

        g_z_0_0_0_0_0_0_0[i] = 2.0 * g_z_0_0_0[i] * a_exp;
    }
}

} // t4c_geom namespace

