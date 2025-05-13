#include "ThreeCenterRR2PrimRecSS.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_ss(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_ss,
                  const size_t idx_ss,
                  const size_t idx_g_ss,
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

    // Set up components of auxiliary buffer : SS

    auto gr_0_0 = pbuffer.data(idx_g_ss);

    // Set up components of targeted buffer : SS

    auto grr_x_0_0 = pbuffer.data(idx_gr_ss);

    auto grr_y_0_0 = pbuffer.data(idx_gr_ss + 1);

    auto grr_z_0_0 = pbuffer.data(idx_gr_ss + 2);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_0, grr_x_0_0, grr_y_0_0, grr_z_0_0, ts_0_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        grr_x_0_0[i] = 2.0 * ts_0_0[i] * gfe_0 * gc_x[i] + gr_0_0[i] * gc_x[i];

        grr_y_0_0[i] = 2.0 * ts_0_0[i] * gfe_0 * gc_y[i] + gr_0_0[i] * gc_y[i];

        grr_z_0_0[i] = 2.0 * ts_0_0[i] * gfe_0 * gc_z[i] + gr_0_0[i] * gc_z[i];
    }
}

} // t3rr2rec namespace

