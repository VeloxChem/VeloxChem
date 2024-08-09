#include "GeomDeriv1010OfScalarForSSSS.hpp"

namespace t4c_geom { // t4c_geom namespace

auto
comp_geom1010_ssss_0(CSimdArray<double>& buffer_1010_ssss,
                     const CSimdArray<double>& buffer_psps,
                     const double a_exp,
                     const double* c_exps) -> void
{
    const auto ndims = buffer_1010_ssss.number_of_columns();

    /// Set up components of auxilary buffer : buffer_psps

    auto g_x_0_x_0 = buffer_psps[0];

    auto g_x_0_y_0 = buffer_psps[1];

    auto g_x_0_z_0 = buffer_psps[2];

    auto g_y_0_x_0 = buffer_psps[3];

    auto g_y_0_y_0 = buffer_psps[4];

    auto g_y_0_z_0 = buffer_psps[5];

    auto g_z_0_x_0 = buffer_psps[6];

    auto g_z_0_y_0 = buffer_psps[7];

    auto g_z_0_z_0 = buffer_psps[8];

    /// Set up components of integrals buffer : buffer_1010_ssss

    auto g_x_0_x_0_0_0_0_0 = buffer_1010_ssss[0];

    auto g_x_0_y_0_0_0_0_0 = buffer_1010_ssss[1];

    auto g_x_0_z_0_0_0_0_0 = buffer_1010_ssss[2];

    auto g_y_0_x_0_0_0_0_0 = buffer_1010_ssss[3];

    auto g_y_0_y_0_0_0_0_0 = buffer_1010_ssss[4];

    auto g_y_0_z_0_0_0_0_0 = buffer_1010_ssss[5];

    auto g_z_0_x_0_0_0_0_0 = buffer_1010_ssss[6];

    auto g_z_0_y_0_0_0_0_0 = buffer_1010_ssss[7];

    auto g_z_0_z_0_0_0_0_0 = buffer_1010_ssss[8];

    // integrals block (0-9)

    #pragma omp simd aligned(g_x_0_x_0, g_x_0_x_0_0_0_0_0, g_x_0_y_0, g_x_0_y_0_0_0_0_0, g_x_0_z_0, g_x_0_z_0_0_0_0_0, g_y_0_x_0, g_y_0_x_0_0_0_0_0, g_y_0_y_0, g_y_0_y_0_0_0_0_0, g_y_0_z_0, g_y_0_z_0_0_0_0_0, g_z_0_x_0, g_z_0_x_0_0_0_0_0, g_z_0_y_0, g_z_0_y_0_0_0_0_0, g_z_0_z_0, g_z_0_z_0_0_0_0_0, c_exps  : 64)
    for (size_t i = 0; i < ndims; i++)
    {
        g_x_0_x_0_0_0_0_0[i] = 4.0 * g_x_0_x_0[i] * a_exp * c_exps[i];

        g_x_0_y_0_0_0_0_0[i] = 4.0 * g_x_0_y_0[i] * a_exp * c_exps[i];

        g_x_0_z_0_0_0_0_0[i] = 4.0 * g_x_0_z_0[i] * a_exp * c_exps[i];

        g_y_0_x_0_0_0_0_0[i] = 4.0 * g_y_0_x_0[i] * a_exp * c_exps[i];

        g_y_0_y_0_0_0_0_0[i] = 4.0 * g_y_0_y_0[i] * a_exp * c_exps[i];

        g_y_0_z_0_0_0_0_0[i] = 4.0 * g_y_0_z_0[i] * a_exp * c_exps[i];

        g_z_0_x_0_0_0_0_0[i] = 4.0 * g_z_0_x_0[i] * a_exp * c_exps[i];

        g_z_0_y_0_0_0_0_0[i] = 4.0 * g_z_0_y_0[i] * a_exp * c_exps[i];

        g_z_0_z_0_0_0_0_0[i] = 4.0 * g_z_0_z_0[i] * a_exp * c_exps[i];
    }
}

} // t4c_geom namespace

