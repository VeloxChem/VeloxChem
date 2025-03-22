#include "TwoCenterElectronRepulsionPrimRecSP.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_sp(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_sp,
                                const size_t idx_eri_1_ss,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpb) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up components of auxiliary buffer : SS

    auto g_0_0_1 = pbuffer.data(idx_eri_1_ss);

    // Set up components of targeted buffer : SP

    auto g_0_x_0 = pbuffer.data(idx_eri_0_sp);

    auto g_0_y_0 = pbuffer.data(idx_eri_0_sp + 1);

    auto g_0_z_0 = pbuffer.data(idx_eri_0_sp + 2);

    #pragma omp simd aligned(g_0_0_1, g_0_x_0, g_0_y_0, g_0_z_0, pb_x, pb_y, pb_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_x_0[i] = g_0_0_1[i] * pb_x[i];

        g_0_y_0[i] = g_0_0_1[i] * pb_y[i];

        g_0_z_0[i] = g_0_0_1[i] * pb_z[i];
    }
}

} // t2ceri namespace

