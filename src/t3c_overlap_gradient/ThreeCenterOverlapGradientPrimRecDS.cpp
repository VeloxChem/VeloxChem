#include "ThreeCenterOverlapGradientPrimRecDS.hpp"

namespace g3ovlrec { // g3ovlrec namespace

auto
comp_prim_overlap_gradient_ds(CSimdArray<double>& pbuffer, 
                              const size_t idx_g_ds,
                              const size_t idx_ps,
                              const size_t idx_ds,
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

    // Set up components of auxiliary buffer : PS

    auto ts_x_0 = pbuffer.data(idx_ps);

    auto ts_y_0 = pbuffer.data(idx_ps + 1);

    auto ts_z_0 = pbuffer.data(idx_ps + 2);

    // Set up components of auxiliary buffer : DS

    auto ts_xx_0 = pbuffer.data(idx_ds);

    auto ts_xy_0 = pbuffer.data(idx_ds + 1);

    auto ts_xz_0 = pbuffer.data(idx_ds + 2);

    auto ts_yy_0 = pbuffer.data(idx_ds + 3);

    auto ts_yz_0 = pbuffer.data(idx_ds + 4);

    auto ts_zz_0 = pbuffer.data(idx_ds + 5);

    // Set up components of targeted buffer : DS

    auto gs_x_xx_0 = pbuffer.data(idx_g_ds);

    auto gs_x_xy_0 = pbuffer.data(idx_g_ds + 1);

    auto gs_x_xz_0 = pbuffer.data(idx_g_ds + 2);

    auto gs_x_yy_0 = pbuffer.data(idx_g_ds + 3);

    auto gs_x_yz_0 = pbuffer.data(idx_g_ds + 4);

    auto gs_x_zz_0 = pbuffer.data(idx_g_ds + 5);

    auto gs_y_xx_0 = pbuffer.data(idx_g_ds + 6);

    auto gs_y_xy_0 = pbuffer.data(idx_g_ds + 7);

    auto gs_y_xz_0 = pbuffer.data(idx_g_ds + 8);

    auto gs_y_yy_0 = pbuffer.data(idx_g_ds + 9);

    auto gs_y_yz_0 = pbuffer.data(idx_g_ds + 10);

    auto gs_y_zz_0 = pbuffer.data(idx_g_ds + 11);

    auto gs_z_xx_0 = pbuffer.data(idx_g_ds + 12);

    auto gs_z_xy_0 = pbuffer.data(idx_g_ds + 13);

    auto gs_z_xz_0 = pbuffer.data(idx_g_ds + 14);

    auto gs_z_yy_0 = pbuffer.data(idx_g_ds + 15);

    auto gs_z_yz_0 = pbuffer.data(idx_g_ds + 16);

    auto gs_z_zz_0 = pbuffer.data(idx_g_ds + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gs_x_xx_0, gs_x_xy_0, gs_x_xz_0, gs_x_yy_0, gs_x_yz_0, gs_x_zz_0, gs_y_xx_0, gs_y_xy_0, gs_y_xz_0, gs_y_yy_0, gs_y_yz_0, gs_y_zz_0, gs_z_xx_0, gs_z_xy_0, gs_z_xz_0, gs_z_yy_0, gs_z_yz_0, gs_z_zz_0, ts_x_0, ts_xx_0, ts_xy_0, ts_xz_0, ts_y_0, ts_yy_0, ts_yz_0, ts_z_0, ts_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double tce_0 = c_exp;

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        gs_x_xx_0[i] = 4.0 * ts_x_0[i] * gfe_0 * tce_0 + 2.0 * ts_xx_0[i] * gc_x[i] * tce_0;

        gs_x_xy_0[i] = 2.0 * ts_y_0[i] * gfe_0 * tce_0 + 2.0 * ts_xy_0[i] * gc_x[i] * tce_0;

        gs_x_xz_0[i] = 2.0 * ts_z_0[i] * gfe_0 * tce_0 + 2.0 * ts_xz_0[i] * gc_x[i] * tce_0;

        gs_x_yy_0[i] = 2.0 * ts_yy_0[i] * gc_x[i] * tce_0;

        gs_x_yz_0[i] = 2.0 * ts_yz_0[i] * gc_x[i] * tce_0;

        gs_x_zz_0[i] = 2.0 * ts_zz_0[i] * gc_x[i] * tce_0;

        gs_y_xx_0[i] = 2.0 * ts_xx_0[i] * gc_y[i] * tce_0;

        gs_y_xy_0[i] = 2.0 * ts_x_0[i] * gfe_0 * tce_0 + 2.0 * ts_xy_0[i] * gc_y[i] * tce_0;

        gs_y_xz_0[i] = 2.0 * ts_xz_0[i] * gc_y[i] * tce_0;

        gs_y_yy_0[i] = 4.0 * ts_y_0[i] * gfe_0 * tce_0 + 2.0 * ts_yy_0[i] * gc_y[i] * tce_0;

        gs_y_yz_0[i] = 2.0 * ts_z_0[i] * gfe_0 * tce_0 + 2.0 * ts_yz_0[i] * gc_y[i] * tce_0;

        gs_y_zz_0[i] = 2.0 * ts_zz_0[i] * gc_y[i] * tce_0;

        gs_z_xx_0[i] = 2.0 * ts_xx_0[i] * gc_z[i] * tce_0;

        gs_z_xy_0[i] = 2.0 * ts_xy_0[i] * gc_z[i] * tce_0;

        gs_z_xz_0[i] = 2.0 * ts_x_0[i] * gfe_0 * tce_0 + 2.0 * ts_xz_0[i] * gc_z[i] * tce_0;

        gs_z_yy_0[i] = 2.0 * ts_yy_0[i] * gc_z[i] * tce_0;

        gs_z_yz_0[i] = 2.0 * ts_y_0[i] * gfe_0 * tce_0 + 2.0 * ts_yz_0[i] * gc_z[i] * tce_0;

        gs_z_zz_0[i] = 4.0 * ts_z_0[i] * gfe_0 * tce_0 + 2.0 * ts_zz_0[i] * gc_z[i] * tce_0;
    }
}

} // g3ovlrec namespace

