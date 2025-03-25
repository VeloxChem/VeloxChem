#include "ThreeCenterOverlapGradientPrimRecSD.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_sd(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_sd,
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

    auto gs_x_0_xx = pbuffer.data(idx_g_sd);

    auto gs_x_0_xy = pbuffer.data(idx_g_sd + 1);

    auto gs_x_0_xz = pbuffer.data(idx_g_sd + 2);

    auto gs_x_0_yy = pbuffer.data(idx_g_sd + 3);

    auto gs_x_0_yz = pbuffer.data(idx_g_sd + 4);

    auto gs_x_0_zz = pbuffer.data(idx_g_sd + 5);

    auto gs_y_0_xx = pbuffer.data(idx_g_sd + 6);

    auto gs_y_0_xy = pbuffer.data(idx_g_sd + 7);

    auto gs_y_0_xz = pbuffer.data(idx_g_sd + 8);

    auto gs_y_0_yy = pbuffer.data(idx_g_sd + 9);

    auto gs_y_0_yz = pbuffer.data(idx_g_sd + 10);

    auto gs_y_0_zz = pbuffer.data(idx_g_sd + 11);

    auto gs_z_0_xx = pbuffer.data(idx_g_sd + 12);

    auto gs_z_0_xy = pbuffer.data(idx_g_sd + 13);

    auto gs_z_0_xz = pbuffer.data(idx_g_sd + 14);

    auto gs_z_0_yy = pbuffer.data(idx_g_sd + 15);

    auto gs_z_0_yz = pbuffer.data(idx_g_sd + 16);

    auto gs_z_0_zz = pbuffer.data(idx_g_sd + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gs_x_0_xx, gs_x_0_xy, gs_x_0_xz, gs_x_0_yy, gs_x_0_yz, gs_x_0_zz, gs_y_0_xx, gs_y_0_xy, gs_y_0_xz, gs_y_0_yy, gs_y_0_yz, gs_y_0_zz, gs_z_0_xx, gs_z_0_xy, gs_z_0_xz, gs_z_0_yy, gs_z_0_yz, gs_z_0_zz, ts_0_x, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_y, ts_0_yy, ts_0_yz, ts_0_z, ts_0_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_0_xx[i] = 4.0 * ts_0_x[i] * gfe_0 * tce_0 + 2.0 * ts_0_xx[i] * gc_x[i] * tce_0;

        gs_x_0_xy[i] = 2.0 * ts_0_y[i] * gfe_0 * tce_0 + 2.0 * ts_0_xy[i] * gc_x[i] * tce_0;

        gs_x_0_xz[i] = 2.0 * ts_0_z[i] * gfe_0 * tce_0 + 2.0 * ts_0_xz[i] * gc_x[i] * tce_0;

        gs_x_0_yy[i] = 2.0 * ts_0_yy[i] * gc_x[i] * tce_0;

        gs_x_0_yz[i] = 2.0 * ts_0_yz[i] * gc_x[i] * tce_0;

        gs_x_0_zz[i] = 2.0 * ts_0_zz[i] * gc_x[i] * tce_0;

        gs_y_0_xx[i] = 2.0 * ts_0_xx[i] * gc_y[i] * tce_0;

        gs_y_0_xy[i] = 2.0 * ts_0_x[i] * gfe_0 * tce_0 + 2.0 * ts_0_xy[i] * gc_y[i] * tce_0;

        gs_y_0_xz[i] = 2.0 * ts_0_xz[i] * gc_y[i] * tce_0;

        gs_y_0_yy[i] = 4.0 * ts_0_y[i] * gfe_0 * tce_0 + 2.0 * ts_0_yy[i] * gc_y[i] * tce_0;

        gs_y_0_yz[i] = 2.0 * ts_0_z[i] * gfe_0 * tce_0 + 2.0 * ts_0_yz[i] * gc_y[i] * tce_0;

        gs_y_0_zz[i] = 2.0 * ts_0_zz[i] * gc_y[i] * tce_0;

        gs_z_0_xx[i] = 2.0 * ts_0_xx[i] * gc_z[i] * tce_0;

        gs_z_0_xy[i] = 2.0 * ts_0_xy[i] * gc_z[i] * tce_0;

        gs_z_0_xz[i] = 2.0 * ts_0_x[i] * gfe_0 * tce_0 + 2.0 * ts_0_xz[i] * gc_z[i] * tce_0;

        gs_z_0_yy[i] = 2.0 * ts_0_yy[i] * gc_z[i] * tce_0;

        gs_z_0_yz[i] = 2.0 * ts_0_y[i] * gfe_0 * tce_0 + 2.0 * ts_0_yz[i] * gc_z[i] * tce_0;

        gs_z_0_zz[i] = 4.0 * ts_0_z[i] * gfe_0 * tce_0 + 2.0 * ts_0_zz[i] * gc_z[i] * tce_0;
    }
}

} // g3ovlrec namespace

