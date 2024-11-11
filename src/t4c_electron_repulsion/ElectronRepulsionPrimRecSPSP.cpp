#include "ElectronRepulsionPrimRecSPSP.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_spsp(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_spsp,
                                  size_t                idx_eri_1_ssss,
                                  size_t                idx_eri_0_sssp,
                                  size_t                idx_eri_1_sssp,
                                  CSimdArray<double>&   factors,
                                  const size_t          idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double          a_exp,
                                  const double          b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SSSS

    auto g_0_0_0_0_1 = pbuffer.data(idx_eri_1_ssss);

    /// Set up components of auxilary buffer : SSSP

    auto g_0_0_0_x_0 = pbuffer.data(idx_eri_0_sssp);

    auto g_0_0_0_y_0 = pbuffer.data(idx_eri_0_sssp + 1);

    auto g_0_0_0_z_0 = pbuffer.data(idx_eri_0_sssp + 2);

    /// Set up components of auxilary buffer : SSSP

    auto g_0_0_0_x_1 = pbuffer.data(idx_eri_1_sssp);

    auto g_0_0_0_y_1 = pbuffer.data(idx_eri_1_sssp + 1);

    auto g_0_0_0_z_1 = pbuffer.data(idx_eri_1_sssp + 2);

    /// Set up 0-3 components of targeted buffer : SPSP

    auto g_0_x_0_x_0 = pbuffer.data(idx_eri_0_spsp);

    auto g_0_x_0_y_0 = pbuffer.data(idx_eri_0_spsp + 1);

    auto g_0_x_0_z_0 = pbuffer.data(idx_eri_0_spsp + 2);

#pragma omp simd aligned( \
        g_0_0_0_0_1, g_0_0_0_x_0, g_0_0_0_x_1, g_0_0_0_y_0, g_0_0_0_y_1, g_0_0_0_z_0, g_0_0_0_z_1, g_0_x_0_x_0, g_0_x_0_y_0, g_0_x_0_z_0, wp_x : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_x_0_x_0[i] = g_0_0_0_0_1[i] * fi_abcd_0 + g_0_0_0_x_0[i] * pb_x + g_0_0_0_x_1[i] * wp_x[i];

        g_0_x_0_y_0[i] = g_0_0_0_y_0[i] * pb_x + g_0_0_0_y_1[i] * wp_x[i];

        g_0_x_0_z_0[i] = g_0_0_0_z_0[i] * pb_x + g_0_0_0_z_1[i] * wp_x[i];
    }

    /// Set up 3-6 components of targeted buffer : SPSP

    auto g_0_y_0_x_0 = pbuffer.data(idx_eri_0_spsp + 3);

    auto g_0_y_0_y_0 = pbuffer.data(idx_eri_0_spsp + 4);

    auto g_0_y_0_z_0 = pbuffer.data(idx_eri_0_spsp + 5);

#pragma omp simd aligned( \
        g_0_0_0_0_1, g_0_0_0_x_0, g_0_0_0_x_1, g_0_0_0_y_0, g_0_0_0_y_1, g_0_0_0_z_0, g_0_0_0_z_1, g_0_y_0_x_0, g_0_y_0_y_0, g_0_y_0_z_0, wp_y : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_y_0_x_0[i] = g_0_0_0_x_0[i] * pb_y + g_0_0_0_x_1[i] * wp_y[i];

        g_0_y_0_y_0[i] = g_0_0_0_0_1[i] * fi_abcd_0 + g_0_0_0_y_0[i] * pb_y + g_0_0_0_y_1[i] * wp_y[i];

        g_0_y_0_z_0[i] = g_0_0_0_z_0[i] * pb_y + g_0_0_0_z_1[i] * wp_y[i];
    }

    /// Set up 6-9 components of targeted buffer : SPSP

    auto g_0_z_0_x_0 = pbuffer.data(idx_eri_0_spsp + 6);

    auto g_0_z_0_y_0 = pbuffer.data(idx_eri_0_spsp + 7);

    auto g_0_z_0_z_0 = pbuffer.data(idx_eri_0_spsp + 8);

#pragma omp simd aligned( \
        g_0_0_0_0_1, g_0_0_0_x_0, g_0_0_0_x_1, g_0_0_0_y_0, g_0_0_0_y_1, g_0_0_0_z_0, g_0_0_0_z_1, g_0_z_0_x_0, g_0_z_0_y_0, g_0_z_0_z_0, wp_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_z_0_x_0[i] = g_0_0_0_x_0[i] * pb_z + g_0_0_0_x_1[i] * wp_z[i];

        g_0_z_0_y_0[i] = g_0_0_0_y_0[i] * pb_z + g_0_0_0_y_1[i] * wp_z[i];

        g_0_z_0_z_0[i] = g_0_0_0_0_1[i] * fi_abcd_0 + g_0_0_0_z_0[i] * pb_z + g_0_0_0_z_1[i] * wp_z[i];
    }
}

}  // namespace erirec
