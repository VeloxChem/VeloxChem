#include "ThreeCenterOverlapGradientPrimRecSP.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_sp(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_sp,
                              const size_t idx_ss,
                              const size_t idx_sp,
                              const CSimdArray<double>& factors,
                              const size_t idx_rgc,
                              const double a_exp,
                              const double c_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(GC) distances

    auto gc_x = factors.data(idx_rgc);

    auto gc_y = factors.data(idx_rgc + 1);

    auto gc_z = factors.data(idx_rgc + 2);

    // Set up components of auxiliary buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ss);

    // Set up components of auxiliary buffer : SP

    auto ts_0_x = pbuffer.data(idx_sp);

    auto ts_0_y = pbuffer.data(idx_sp + 1);

    auto ts_0_z = pbuffer.data(idx_sp + 2);

    // Set up components of targeted buffer : SP

    auto gs_x_0_x = pbuffer.data(idx_g_sp);

    auto gs_x_0_y = pbuffer.data(idx_g_sp + 1);

    auto gs_x_0_z = pbuffer.data(idx_g_sp + 2);

    auto gs_y_0_x = pbuffer.data(idx_g_sp + 3);

    auto gs_y_0_y = pbuffer.data(idx_g_sp + 4);

    auto gs_y_0_z = pbuffer.data(idx_g_sp + 5);

    auto gs_z_0_x = pbuffer.data(idx_g_sp + 6);

    auto gs_z_0_y = pbuffer.data(idx_g_sp + 7);

    auto gs_z_0_z = pbuffer.data(idx_g_sp + 8);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gs_x_0_x, gs_x_0_y, gs_x_0_z, gs_y_0_x, gs_y_0_y, gs_y_0_z, gs_z_0_x, gs_z_0_y, gs_z_0_z, ts_0_0, ts_0_x, ts_0_y, ts_0_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_0_x[i] = 2.0 * ts_0_0[i] * gfe_0 * tce_0 + 2.0 * ts_0_x[i] * gc_x[i] * tce_0;

        gs_x_0_y[i] = 2.0 * ts_0_y[i] * gc_x[i] * tce_0;

        gs_x_0_z[i] = 2.0 * ts_0_z[i] * gc_x[i] * tce_0;

        gs_y_0_x[i] = 2.0 * ts_0_x[i] * gc_y[i] * tce_0;

        gs_y_0_y[i] = 2.0 * ts_0_0[i] * gfe_0 * tce_0 + 2.0 * ts_0_y[i] * gc_y[i] * tce_0;

        gs_y_0_z[i] = 2.0 * ts_0_z[i] * gc_y[i] * tce_0;

        gs_z_0_x[i] = 2.0 * ts_0_x[i] * gc_z[i] * tce_0;

        gs_z_0_y[i] = 2.0 * ts_0_y[i] * gc_z[i] * tce_0;

        gs_z_0_z[i] = 2.0 * ts_0_0[i] * gfe_0 * tce_0 + 2.0 * ts_0_z[i] * gc_z[i] * tce_0;
    }
}

} // g3ovlrec namespace

