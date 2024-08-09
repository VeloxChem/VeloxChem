#include "ElectronRepulsionPrimRecSPSS.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_spss(CSimdArray<double>& prim_buffer_0_spss,
                                  const CSimdArray<double>& prim_buffer_0_ssss,
                                  const CSimdArray<double>& prim_buffer_1_ssss,
                                  const double pb_x,
                                  const double pb_y,
                                  const double pb_z,
                                  const double* wp_x,
                                  const double* wp_y,
                                  const double* wp_z) -> void
{
    const auto ndims = prim_buffer_0_spss.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_ssss

    auto g_0_0_0_0_0 = prim_buffer_0_ssss[0];

    /// Set up components of auxilary buffer : prim_buffer_1_ssss

    auto g_0_0_0_0_1 = prim_buffer_1_ssss[0];

    /// Set up components of targeted buffer : prim_buffer_0_spss

    auto g_0_x_0_0_0 = prim_buffer_0_spss[0];

    auto g_0_y_0_0_0 = prim_buffer_0_spss[1];

    auto g_0_z_0_0_0 = prim_buffer_0_spss[2];

    #pragma omp simd aligned(g_0_0_0_0_0, g_0_0_0_0_1, g_0_x_0_0_0, g_0_y_0_0_0, g_0_z_0_0_0, wp_x, wp_y, wp_z  : 64)
    for (int i = 0; i < ndims; i++)
    {
        g_0_x_0_0_0[i] = g_0_0_0_0_0[i] * pb_x + g_0_0_0_0_1[i] * wp_x[i];

        g_0_y_0_0_0[i] = g_0_0_0_0_0[i] * pb_y + g_0_0_0_0_1[i] * wp_y[i];

        g_0_z_0_0_0[i] = g_0_0_0_0_0[i] * pb_z + g_0_0_0_0_1[i] * wp_z[i];
    }
}

} // erirec namespace

