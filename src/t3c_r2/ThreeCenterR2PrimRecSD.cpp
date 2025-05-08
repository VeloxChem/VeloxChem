#include "ThreeCenterR2PrimRecSD.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_sd(CSimdArray<double>& pbuffer, 
                const size_t idx_g_sd,
                const size_t idx_ss,
                const size_t idx_sp,
                const size_t idx_sd,
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

    // Set up components of auxiliary buffer : SD

    auto ts_0_xx = pbuffer.data(idx_sd);

    auto ts_0_xy = pbuffer.data(idx_sd + 1);

    auto ts_0_xz = pbuffer.data(idx_sd + 2);

    auto ts_0_yy = pbuffer.data(idx_sd + 3);

    auto ts_0_yz = pbuffer.data(idx_sd + 4);

    auto ts_0_zz = pbuffer.data(idx_sd + 5);

    // Set up components of targeted buffer : SD

    auto gr_0_xx = pbuffer.data(idx_g_sd);

    auto gr_0_xy = pbuffer.data(idx_g_sd + 1);

    auto gr_0_xz = pbuffer.data(idx_g_sd + 2);

    auto gr_0_yy = pbuffer.data(idx_g_sd + 3);

    auto gr_0_yz = pbuffer.data(idx_g_sd + 4);

    auto gr_0_zz = pbuffer.data(idx_g_sd + 5);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_xx, gr_0_xy, gr_0_xz, gr_0_yy, gr_0_yz, gr_0_zz, ts_0_0, ts_0_x, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_y, ts_0_yy, ts_0_yz, ts_0_z, ts_0_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gr_0_xx[i] = 2.0 * ts_0_0[i] * gfe_0 + 4.0 * ts_0_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_0_xx[i] * gfe_0 + ts_0_xx[i] * rgc2_0;

        gr_0_xy[i] = 2.0 * ts_0_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_0_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_0_xy[i] * gfe_0 + ts_0_xy[i] * rgc2_0;

        gr_0_xz[i] = 2.0 * ts_0_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_0_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_xz[i] * gfe_0 + ts_0_xz[i] * rgc2_0;

        gr_0_yy[i] = 2.0 * ts_0_0[i] * gfe_0 + 4.0 * ts_0_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_0_yy[i] * gfe_0 + ts_0_yy[i] * rgc2_0;

        gr_0_yz[i] = 2.0 * ts_0_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_0_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_yz[i] * gfe_0 + ts_0_yz[i] * rgc2_0;

        gr_0_zz[i] = 2.0 * ts_0_0[i] * gfe_0 + 4.0 * ts_0_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_0_zz[i] * gfe_0 + ts_0_zz[i] * rgc2_0;
    }
}

} // t3r2rec namespace

