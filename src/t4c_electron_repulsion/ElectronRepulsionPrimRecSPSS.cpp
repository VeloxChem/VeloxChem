#include "ElectronRepulsionPrimRecSPSS.hpp"

namespace erirec { // erirec namespace

auto
comp_prim_electron_repulsion_spss(CSimdArray<double>& pbuffer,
                                  const size_t idx_eri_0_spss,
                                  size_t idx_eri_0_ssss,
                                  size_t idx_eri_1_ssss,
                                  CSimdArray<double>& factors,
                                  const size_t idx_wp,
                                  const TPoint<double>& r_pb) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

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

    auto g_0_0_0_0_0 = pbuffer.data(idx_eri_0_ssss);

    /// Set up components of auxilary buffer : SSSS

    auto g_0_0_0_0_1 = pbuffer.data(idx_eri_1_ssss);

    /// Set up components of targeted buffer : SPSS

    auto g_0_x_0_0_0 = pbuffer.data(idx_eri_0_spss);

    auto g_0_y_0_0_0 = pbuffer.data(idx_eri_0_spss + 1);

    auto g_0_z_0_0_0 = pbuffer.data(idx_eri_0_spss + 2);

    #pragma omp simd aligned(g_0_0_0_0_0, g_0_0_0_0_1, g_0_x_0_0_0, g_0_y_0_0_0, g_0_z_0_0_0, wp_x, wp_y, wp_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_x_0_0_0[i] = g_0_0_0_0_0[i] * pb_x + g_0_0_0_0_1[i] * wp_x[i];

        g_0_y_0_0_0[i] = g_0_0_0_0_0[i] * pb_y + g_0_0_0_0_1[i] * wp_y[i];

        g_0_z_0_0_0[i] = g_0_0_0_0_0[i] * pb_z + g_0_0_0_0_1[i] * wp_z[i];
    }
}

} // erirec namespace

