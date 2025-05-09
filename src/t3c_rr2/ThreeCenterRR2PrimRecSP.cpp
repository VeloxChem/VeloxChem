#include "ThreeCenterRR2PrimRecSP.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_sp(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_sp,
                  const size_t idx_ss,
                  const size_t idx_g_ss,
                  const size_t idx_sp,
                  const size_t idx_g_sp,
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

    // Set up components of auxiliary buffer : SP

    auto ts_0_x = pbuffer.data(idx_sp);

    auto ts_0_y = pbuffer.data(idx_sp + 1);

    auto ts_0_z = pbuffer.data(idx_sp + 2);

    // Set up components of auxiliary buffer : SP

    auto gr_0_x = pbuffer.data(idx_g_sp);

    auto gr_0_y = pbuffer.data(idx_g_sp + 1);

    auto gr_0_z = pbuffer.data(idx_g_sp + 2);

    // Set up components of targeted buffer : SP

    auto grr_x_0_x = pbuffer.data(idx_gr_sp);

    auto grr_x_0_y = pbuffer.data(idx_gr_sp + 1);

    auto grr_x_0_z = pbuffer.data(idx_gr_sp + 2);

    auto grr_y_0_x = pbuffer.data(idx_gr_sp + 3);

    auto grr_y_0_y = pbuffer.data(idx_gr_sp + 4);

    auto grr_y_0_z = pbuffer.data(idx_gr_sp + 5);

    auto grr_z_0_x = pbuffer.data(idx_gr_sp + 6);

    auto grr_z_0_y = pbuffer.data(idx_gr_sp + 7);

    auto grr_z_0_z = pbuffer.data(idx_gr_sp + 8);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_0, gr_0_x, gr_0_y, gr_0_z, grr_x_0_x, grr_x_0_y, grr_x_0_z, grr_y_0_x, grr_y_0_y, grr_y_0_z, grr_z_0_x, grr_z_0_y, grr_z_0_z, ts_0_0, ts_0_x, ts_0_y, ts_0_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_0_x[i] = ts_0_0[i] * gfe2_0 + gr_0_0[i] * gfe_0 + ts_0_x[i] * gfe_0 * gc_x[i] + gr_0_x[i] * gc_x[i];

        grr_x_0_y[i] = ts_0_y[i] * gfe_0 * gc_x[i] + gr_0_y[i] * gc_x[i];

        grr_x_0_z[i] = ts_0_z[i] * gfe_0 * gc_x[i] + gr_0_z[i] * gc_x[i];

        grr_y_0_x[i] = ts_0_x[i] * gfe_0 * gc_y[i] + gr_0_x[i] * gc_y[i];

        grr_y_0_y[i] = ts_0_0[i] * gfe2_0 + gr_0_0[i] * gfe_0 + ts_0_y[i] * gfe_0 * gc_y[i] + gr_0_y[i] * gc_y[i];

        grr_y_0_z[i] = ts_0_z[i] * gfe_0 * gc_y[i] + gr_0_z[i] * gc_y[i];

        grr_z_0_x[i] = ts_0_x[i] * gfe_0 * gc_z[i] + gr_0_x[i] * gc_z[i];

        grr_z_0_y[i] = ts_0_y[i] * gfe_0 * gc_z[i] + gr_0_y[i] * gc_z[i];

        grr_z_0_z[i] = ts_0_0[i] * gfe2_0 + gr_0_0[i] * gfe_0 + ts_0_z[i] * gfe_0 * gc_z[i] + gr_0_z[i] * gc_z[i];
    }
}

} // t3rr2rec namespace

