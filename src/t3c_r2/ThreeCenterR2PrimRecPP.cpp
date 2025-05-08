#include "ThreeCenterR2PrimRecPP.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_pp(CSimdArray<double>& pbuffer, 
                const size_t idx_g_pp,
                const size_t idx_ss,
                const size_t idx_sp,
                const size_t idx_ps,
                const size_t idx_pp,
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

    // Set up components of auxiliary buffer : PS

    auto ts_x_0 = pbuffer.data(idx_ps);

    auto ts_y_0 = pbuffer.data(idx_ps + 1);

    auto ts_z_0 = pbuffer.data(idx_ps + 2);

    // Set up components of auxiliary buffer : PP

    auto ts_x_x = pbuffer.data(idx_pp);

    auto ts_x_y = pbuffer.data(idx_pp + 1);

    auto ts_x_z = pbuffer.data(idx_pp + 2);

    auto ts_y_x = pbuffer.data(idx_pp + 3);

    auto ts_y_y = pbuffer.data(idx_pp + 4);

    auto ts_y_z = pbuffer.data(idx_pp + 5);

    auto ts_z_x = pbuffer.data(idx_pp + 6);

    auto ts_z_y = pbuffer.data(idx_pp + 7);

    auto ts_z_z = pbuffer.data(idx_pp + 8);

    // Set up 0-3 components of targeted buffer : PP

    auto gr_x_x = pbuffer.data(idx_g_pp);

    auto gr_x_y = pbuffer.data(idx_g_pp + 1);

    auto gr_x_z = pbuffer.data(idx_g_pp + 2);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_x, gr_x_y, gr_x_z, ts_0_0, ts_0_x, ts_0_y, ts_0_z, ts_x_0, ts_x_x, ts_x_y, ts_x_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_x_x[i] = 2.0 * ts_0_0[i] * gfe2_0 + 2.0 * ts_0_x[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_x_x[i] * gfe_0 + ts_x_x[i] * rgc2_0;

        gr_x_y[i] = 2.0 * ts_0_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_x_y[i] * gfe_0 + ts_x_y[i] * rgc2_0;

        gr_x_z[i] = 2.0 * ts_0_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_z[i] * gfe_0 + ts_x_z[i] * rgc2_0;
    }

    // Set up 3-6 components of targeted buffer : PP

    auto gr_y_x = pbuffer.data(idx_g_pp + 3);

    auto gr_y_y = pbuffer.data(idx_g_pp + 4);

    auto gr_y_z = pbuffer.data(idx_g_pp + 5);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_x, gr_y_y, gr_y_z, ts_0_0, ts_0_x, ts_0_y, ts_0_z, ts_y_0, ts_y_x, ts_y_y, ts_y_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_y_x[i] = 2.0 * ts_0_x[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_y_x[i] * gfe_0 + ts_y_x[i] * rgc2_0;

        gr_y_y[i] = 2.0 * ts_0_0[i] * gfe2_0 + 2.0 * ts_0_y[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_y_y[i] * gfe_0 + ts_y_y[i] * rgc2_0;

        gr_y_z[i] = 2.0 * ts_0_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_z[i] * gfe_0 + ts_y_z[i] * rgc2_0;
    }

    // Set up 6-9 components of targeted buffer : PP

    auto gr_z_x = pbuffer.data(idx_g_pp + 6);

    auto gr_z_y = pbuffer.data(idx_g_pp + 7);

    auto gr_z_z = pbuffer.data(idx_g_pp + 8);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_z_x, gr_z_y, gr_z_z, ts_0_0, ts_0_x, ts_0_y, ts_0_z, ts_z_0, ts_z_x, ts_z_y, ts_z_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_z_x[i] = 2.0 * ts_0_x[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_0[i] * gfe_0 * gc_x[i] + 3.0 * ts_z_x[i] * gfe_0 + ts_z_x[i] * rgc2_0;

        gr_z_y[i] = 2.0 * ts_0_y[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_0[i] * gfe_0 * gc_y[i] + 3.0 * ts_z_y[i] * gfe_0 + ts_z_y[i] * rgc2_0;

        gr_z_z[i] = 2.0 * ts_0_0[i] * gfe2_0 + 2.0 * ts_0_z[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_0[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_z[i] * gfe_0 + ts_z_z[i] * rgc2_0;
    }

}

} // t3r2rec namespace

