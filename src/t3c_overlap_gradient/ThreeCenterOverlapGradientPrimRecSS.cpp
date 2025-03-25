#include "ThreeCenterOverlapGradientPrimRecSS.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_ss(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_ss,
                              const size_t idx_ss,
                              const CSimdArray<double>& factors,
                              const size_t idx_rgc,
                              const double a_exp,
                              const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up R(GC) distances

    auto gc_x = factors.data(idx_rgc);

    auto gc_y = factors.data(idx_rgc + 1);

    auto gc_z = factors.data(idx_rgc + 2);

    // Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ss);

    // Set up components of targeted buffer : SS

    auto gs_x_0_0 = pbuffer.data(idx_g_ss);

    auto gs_y_0_0 = pbuffer.data(idx_g_ss + 1);

    auto gs_z_0_0 = pbuffer.data(idx_g_ss + 2);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gs_x_0_0, gs_y_0_0, gs_z_0_0, ts_0_0  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        gs_x_0_0[i] = 2.0 * ts_0_0[i] * gc_x[i] * tce_0;

        gs_y_0_0[i] = 2.0 * ts_0_0[i] * gc_y[i] * tce_0;

        gs_z_0_0[i] = 2.0 * ts_0_0[i] * gc_z[i] * tce_0;
    }
}

} // g3ovlrec namespace

