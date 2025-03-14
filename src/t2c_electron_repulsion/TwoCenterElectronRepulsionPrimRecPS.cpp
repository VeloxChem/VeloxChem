#include "TwoCenterElectronRepulsionPrimRecPS.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_ps(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_ps,
                                const size_t idx_eri_1_ss,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : SS

    auto g_0_0_1 = pbuffer.data(idx_eri_1_ss);

    // Set up components of targeted buffer : PS

    auto g_x_0_0 = pbuffer.data(idx_eri_0_ps);

    auto g_y_0_0 = pbuffer.data(idx_eri_0_ps + 1);

    auto g_z_0_0 = pbuffer.data(idx_eri_0_ps + 2);

    #pragma omp simd aligned(g_0_0_1, g_x_0_0, g_y_0_0, g_z_0_0, pa_x, pa_y, pa_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_x_0_0[i] = g_0_0_1[i] * pa_x[i];

        g_y_0_0[i] = g_0_0_1[i] * pa_y[i];

        g_z_0_0[i] = g_0_0_1[i] * pa_z[i];
    }
}

} // t2ceri namespace

