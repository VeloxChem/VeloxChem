#include "ElectronRepulsionPrimRecSPSP.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_spsp(CSimdArray<double>& prim_buffer_0_spsp,
                                  const CSimdArray<double>& prim_buffer_1_ssss,
                                  const CSimdArray<double>& prim_buffer_0_sssp,
                                  const CSimdArray<double>& prim_buffer_1_sssp,
                                  const double pb_x,
                                  const double pb_y,
                                  const double pb_z,
                                  const double* wp_x,
                                  const double* wp_y,
                                  const double* wp_z,
                                  const double a_exp,
                                  const double b_exp,
                                  const double* c_exps,
                                  const double* d_exps) -> void
{
    const auto ndims = prim_buffer_0_spsp.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_1_ssss

    auto g_0_0_0_0_1 = prim_buffer_1_ssss[0];

    /// Set up components of auxilary buffer : prim_buffer_0_sssp

    auto g_0_0_0_x_0 = prim_buffer_0_sssp[0];

    auto g_0_0_0_y_0 = prim_buffer_0_sssp[1];

    auto g_0_0_0_z_0 = prim_buffer_0_sssp[2];

    /// Set up components of auxilary buffer : prim_buffer_1_sssp

    auto g_0_0_0_x_1 = prim_buffer_1_sssp[0];

    auto g_0_0_0_y_1 = prim_buffer_1_sssp[1];

    auto g_0_0_0_z_1 = prim_buffer_1_sssp[2];

    /// Set up 0-3 components of targeted buffer : prim_buffer_0_spsp

    auto g_0_x_0_x_0 = prim_buffer_0_spsp[0];

    auto g_0_x_0_y_0 = prim_buffer_0_spsp[1];

    auto g_0_x_0_z_0 = prim_buffer_0_spsp[2];

    #pragma omp simd aligned(g_0_0_0_0_1, g_0_0_0_x_0, g_0_0_0_x_1, g_0_0_0_y_0, g_0_0_0_y_1, g_0_0_0_z_0, g_0_0_0_z_1, g_0_x_0_x_0, g_0_x_0_y_0, g_0_x_0_z_0, wp_x  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_x_0_x_0[i] = g_0_0_0_0_1[i] * fi_abcd_0 + g_0_0_0_x_0[i] * pb_x + g_0_0_0_x_1[i] * wp_x[i];

        g_0_x_0_y_0[i] = g_0_0_0_y_0[i] * pb_x + g_0_0_0_y_1[i] * wp_x[i];

        g_0_x_0_z_0[i] = g_0_0_0_z_0[i] * pb_x + g_0_0_0_z_1[i] * wp_x[i];
    }

    /// Set up 3-6 components of targeted buffer : prim_buffer_0_spsp

    auto g_0_y_0_x_0 = prim_buffer_0_spsp[3];

    auto g_0_y_0_y_0 = prim_buffer_0_spsp[4];

    auto g_0_y_0_z_0 = prim_buffer_0_spsp[5];

    #pragma omp simd aligned(g_0_0_0_0_1, g_0_0_0_x_0, g_0_0_0_x_1, g_0_0_0_y_0, g_0_0_0_y_1, g_0_0_0_z_0, g_0_0_0_z_1, g_0_y_0_x_0, g_0_y_0_y_0, g_0_y_0_z_0, wp_y  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_y_0_x_0[i] = g_0_0_0_x_0[i] * pb_y + g_0_0_0_x_1[i] * wp_y[i];

        g_0_y_0_y_0[i] = g_0_0_0_0_1[i] * fi_abcd_0 + g_0_0_0_y_0[i] * pb_y + g_0_0_0_y_1[i] * wp_y[i];

        g_0_y_0_z_0[i] = g_0_0_0_z_0[i] * pb_y + g_0_0_0_z_1[i] * wp_y[i];
    }

    /// Set up 6-9 components of targeted buffer : prim_buffer_0_spsp

    auto g_0_z_0_x_0 = prim_buffer_0_spsp[6];

    auto g_0_z_0_y_0 = prim_buffer_0_spsp[7];

    auto g_0_z_0_z_0 = prim_buffer_0_spsp[8];

    #pragma omp simd aligned(g_0_0_0_0_1, g_0_0_0_x_0, g_0_0_0_x_1, g_0_0_0_y_0, g_0_0_0_y_1, g_0_0_0_z_0, g_0_0_0_z_1, g_0_z_0_x_0, g_0_z_0_y_0, g_0_z_0_z_0, wp_z  : 64)
    for (int i = 0; i < ndims; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_z_0_x_0[i] = g_0_0_0_x_0[i] * pb_z + g_0_0_0_x_1[i] * wp_z[i];

        g_0_z_0_y_0[i] = g_0_0_0_y_0[i] * pb_z + g_0_0_0_y_1[i] * wp_z[i];

        g_0_z_0_z_0[i] = g_0_0_0_0_1[i] * fi_abcd_0 + g_0_0_0_z_0[i] * pb_z + g_0_0_0_z_1[i] * wp_z[i];
    }
}

} // erirec namespace

