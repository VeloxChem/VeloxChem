#include "ThreeCenterOverlapGradientPrimRecPS.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_ps(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_ps,
                              const size_t idx_ss,
                              const size_t idx_ps,
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

    // Set up components of auxiliary buffer : PS

    auto ts_x_0 = pbuffer.data(idx_ps);

    auto ts_y_0 = pbuffer.data(idx_ps + 1);

    auto ts_z_0 = pbuffer.data(idx_ps + 2);

    // Set up components of targeted buffer : PS

    auto gs_x_x_0 = pbuffer.data(idx_g_ps);

    auto gs_x_y_0 = pbuffer.data(idx_g_ps + 1);

    auto gs_x_z_0 = pbuffer.data(idx_g_ps + 2);

    auto gs_y_x_0 = pbuffer.data(idx_g_ps + 3);

    auto gs_y_y_0 = pbuffer.data(idx_g_ps + 4);

    auto gs_y_z_0 = pbuffer.data(idx_g_ps + 5);

    auto gs_z_x_0 = pbuffer.data(idx_g_ps + 6);

    auto gs_z_y_0 = pbuffer.data(idx_g_ps + 7);

    auto gs_z_z_0 = pbuffer.data(idx_g_ps + 8);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gs_x_x_0, gs_x_y_0, gs_x_z_0, gs_y_x_0, gs_y_y_0, gs_y_z_0, gs_z_x_0, gs_z_y_0, gs_z_z_0, ts_0_0, ts_x_0, ts_y_0, ts_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_x_0[i] = 2.0 * ts_0_0[i] * gfe_0 * tce_0 + 2.0 * ts_x_0[i] * gc_x[i] * tce_0;

        gs_x_y_0[i] = 2.0 * ts_y_0[i] * gc_x[i] * tce_0;

        gs_x_z_0[i] = 2.0 * ts_z_0[i] * gc_x[i] * tce_0;

        gs_y_x_0[i] = 2.0 * ts_x_0[i] * gc_y[i] * tce_0;

        gs_y_y_0[i] = 2.0 * ts_0_0[i] * gfe_0 * tce_0 + 2.0 * ts_y_0[i] * gc_y[i] * tce_0;

        gs_y_z_0[i] = 2.0 * ts_z_0[i] * gc_y[i] * tce_0;

        gs_z_x_0[i] = 2.0 * ts_x_0[i] * gc_z[i] * tce_0;

        gs_z_y_0[i] = 2.0 * ts_y_0[i] * gc_z[i] * tce_0;

        gs_z_z_0[i] = 2.0 * ts_0_0[i] * gfe_0 * tce_0 + 2.0 * ts_z_0[i] * gc_z[i] * tce_0;
    }
}

} // g3ovlrec namespace

