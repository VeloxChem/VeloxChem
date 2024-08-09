#include "ElectronRepulsionPrimRecSSSP.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_sssp(CSimdArray<double>& prim_buffer_0_sssp,
                                  const CSimdArray<double>& prim_buffer_0_ssss,
                                  const CSimdArray<double>& prim_buffer_1_ssss,
                                  const double* qd_x,
                                  const double* qd_y,
                                  const double* qd_z,
                                  const double* wq_x,
                                  const double* wq_y,
                                  const double* wq_z) -> void
{
    const auto ndims = prim_buffer_0_sssp.number_of_columns();

    /// Set up components of auxilary buffer : prim_buffer_0_ssss

    auto g_0_0_0_0_0 = prim_buffer_0_ssss[0];

    /// Set up components of auxilary buffer : prim_buffer_1_ssss

    auto g_0_0_0_0_1 = prim_buffer_1_ssss[0];

    /// Set up components of targeted buffer : prim_buffer_0_sssp

    auto g_0_0_0_x_0 = prim_buffer_0_sssp[0];

    auto g_0_0_0_y_0 = prim_buffer_0_sssp[1];

    auto g_0_0_0_z_0 = prim_buffer_0_sssp[2];

    #pragma omp simd aligned(g_0_0_0_0_0, g_0_0_0_0_1, g_0_0_0_x_0, g_0_0_0_y_0, g_0_0_0_z_0, qd_x, qd_y, qd_z, wq_x, wq_y, wq_z  : 64)
    for (int i = 0; i < ndims; i++)
    {
        g_0_0_0_x_0[i] = g_0_0_0_0_0[i] * qd_x[i] + g_0_0_0_0_1[i] * wq_x[i];

        g_0_0_0_y_0[i] = g_0_0_0_0_0[i] * qd_y[i] + g_0_0_0_0_1[i] * wq_y[i];

        g_0_0_0_z_0[i] = g_0_0_0_0_0[i] * qd_z[i] + g_0_0_0_0_1[i] * wq_z[i];
    }
}

} // erirec namespace

