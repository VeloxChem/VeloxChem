#include "ThreeCenterR2PrimRecPS.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_ps(CSimdArray<double>& pbuffer, 
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

    auto gr_x_0 = pbuffer.data(idx_g_ps);

    auto gr_y_0 = pbuffer.data(idx_g_ps + 1);

    auto gr_z_0 = pbuffer.data(idx_g_ps + 2);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_0, gr_y_0, gr_z_0, ts_0_0, ts_x_0, ts_y_0, ts_z_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_x_0[i] = 2.0 * ts_0_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_x_0[i] * gfe_0 + ts_x_0[i] * rgc2_0;

        gr_y_0[i] = 2.0 * ts_0_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_y_0[i] * gfe_0 + ts_y_0[i] * rgc2_0;

        gr_z_0[i] = 2.0 * ts_0_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_0[i] * gfe_0 + ts_z_0[i] * rgc2_0;
    }
}

} // t3r2rec namespace

