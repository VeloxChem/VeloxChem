#include "ThreeCenterElectronRepulsionPrimRecPSS.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_pss(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_pss,
                                 size_t idx_eri_1_sss,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa) -> void 
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(WA) distances

    auto wa_x = factors.data(idx_wa);

    auto wa_y = factors.data(idx_wa + 1);

    auto wa_z = factors.data(idx_wa + 2);

    /// Set up components of auxilary buffer : SSS

    auto g_0_0_0_1 = pbuffer.data(idx_eri_1_sss);

    /// Set up components of targeted buffer : PSS

    auto g_x_0_0_0 = pbuffer.data(idx_eri_0_pss);

    auto g_y_0_0_0 = pbuffer.data(idx_eri_0_pss + 1);

    auto g_z_0_0_0 = pbuffer.data(idx_eri_0_pss + 2);

    #pragma omp simd aligned(g_0_0_0_1, g_x_0_0_0, g_y_0_0_0, g_z_0_0_0, wa_x, wa_y, wa_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_x_0_0_0[i] = g_0_0_0_1[i] * wa_x[i];

        g_y_0_0_0[i] = g_0_0_0_1[i] * wa_y[i];

        g_z_0_0_0[i] = g_0_0_0_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

