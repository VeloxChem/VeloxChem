#include "ThreeCenterRR2PrimRecDS.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_ds(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_ds,
                  const size_t idx_ps,
                  const size_t idx_g_ps,
                  const size_t idx_ds,
                  const size_t idx_g_ds,
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

    // Set up components of auxiliary buffer : PS

    auto gr_x_0 = pbuffer.data(idx_g_ps);

    auto gr_y_0 = pbuffer.data(idx_g_ps + 1);

    auto gr_z_0 = pbuffer.data(idx_g_ps + 2);

    // Set up components of auxiliary buffer : DS

    auto ts_xx_0 = pbuffer.data(idx_ds);

    auto ts_xy_0 = pbuffer.data(idx_ds + 1);

    auto ts_xz_0 = pbuffer.data(idx_ds + 2);

    auto ts_yy_0 = pbuffer.data(idx_ds + 3);

    auto ts_yz_0 = pbuffer.data(idx_ds + 4);

    auto ts_zz_0 = pbuffer.data(idx_ds + 5);

    // Set up components of auxiliary buffer : DS

    auto gr_xx_0 = pbuffer.data(idx_g_ds);

    auto gr_xy_0 = pbuffer.data(idx_g_ds + 1);

    auto gr_xz_0 = pbuffer.data(idx_g_ds + 2);

    auto gr_yy_0 = pbuffer.data(idx_g_ds + 3);

    auto gr_yz_0 = pbuffer.data(idx_g_ds + 4);

    auto gr_zz_0 = pbuffer.data(idx_g_ds + 5);

    // Set up components of targeted buffer : DS

    auto grr_x_xx_0 = pbuffer.data(idx_gr_ds);

    auto grr_x_xy_0 = pbuffer.data(idx_gr_ds + 1);

    auto grr_x_xz_0 = pbuffer.data(idx_gr_ds + 2);

    auto grr_x_yy_0 = pbuffer.data(idx_gr_ds + 3);

    auto grr_x_yz_0 = pbuffer.data(idx_gr_ds + 4);

    auto grr_x_zz_0 = pbuffer.data(idx_gr_ds + 5);

    auto grr_y_xx_0 = pbuffer.data(idx_gr_ds + 6);

    auto grr_y_xy_0 = pbuffer.data(idx_gr_ds + 7);

    auto grr_y_xz_0 = pbuffer.data(idx_gr_ds + 8);

    auto grr_y_yy_0 = pbuffer.data(idx_gr_ds + 9);

    auto grr_y_yz_0 = pbuffer.data(idx_gr_ds + 10);

    auto grr_y_zz_0 = pbuffer.data(idx_gr_ds + 11);

    auto grr_z_xx_0 = pbuffer.data(idx_gr_ds + 12);

    auto grr_z_xy_0 = pbuffer.data(idx_gr_ds + 13);

    auto grr_z_xz_0 = pbuffer.data(idx_gr_ds + 14);

    auto grr_z_yy_0 = pbuffer.data(idx_gr_ds + 15);

    auto grr_z_yz_0 = pbuffer.data(idx_gr_ds + 16);

    auto grr_z_zz_0 = pbuffer.data(idx_gr_ds + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_0, gr_xx_0, gr_xy_0, gr_xz_0, gr_y_0, gr_yy_0, gr_yz_0, gr_z_0, gr_zz_0, grr_x_xx_0, grr_x_xy_0, grr_x_xz_0, grr_x_yy_0, grr_x_yz_0, grr_x_zz_0, grr_y_xx_0, grr_y_xy_0, grr_y_xz_0, grr_y_yy_0, grr_y_yz_0, grr_y_zz_0, grr_z_xx_0, grr_z_xy_0, grr_z_xz_0, grr_z_yy_0, grr_z_yz_0, grr_z_zz_0, ts_x_0, ts_xx_0, ts_xy_0, ts_xz_0, ts_y_0, ts_yy_0, ts_yz_0, ts_z_0, ts_zz_0, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_xx_0[i] = 4.0 * ts_x_0[i] * gfe2_0 + 2.0 * gr_x_0[i] * gfe_0 + 2.0 * ts_xx_0[i] * gfe_0 * gc_x[i] + gr_xx_0[i] * gc_x[i];

        grr_x_xy_0[i] = 2.0 * ts_y_0[i] * gfe2_0 + gr_y_0[i] * gfe_0 + 2.0 * ts_xy_0[i] * gfe_0 * gc_x[i] + gr_xy_0[i] * gc_x[i];

        grr_x_xz_0[i] = 2.0 * ts_z_0[i] * gfe2_0 + gr_z_0[i] * gfe_0 + 2.0 * ts_xz_0[i] * gfe_0 * gc_x[i] + gr_xz_0[i] * gc_x[i];

        grr_x_yy_0[i] = 2.0 * ts_yy_0[i] * gfe_0 * gc_x[i] + gr_yy_0[i] * gc_x[i];

        grr_x_yz_0[i] = 2.0 * ts_yz_0[i] * gfe_0 * gc_x[i] + gr_yz_0[i] * gc_x[i];

        grr_x_zz_0[i] = 2.0 * ts_zz_0[i] * gfe_0 * gc_x[i] + gr_zz_0[i] * gc_x[i];

        grr_y_xx_0[i] = 2.0 * ts_xx_0[i] * gfe_0 * gc_y[i] + gr_xx_0[i] * gc_y[i];

        grr_y_xy_0[i] = 2.0 * ts_x_0[i] * gfe2_0 + gr_x_0[i] * gfe_0 + 2.0 * ts_xy_0[i] * gfe_0 * gc_y[i] + gr_xy_0[i] * gc_y[i];

        grr_y_xz_0[i] = 2.0 * ts_xz_0[i] * gfe_0 * gc_y[i] + gr_xz_0[i] * gc_y[i];

        grr_y_yy_0[i] = 4.0 * ts_y_0[i] * gfe2_0 + 2.0 * gr_y_0[i] * gfe_0 + 2.0 * ts_yy_0[i] * gfe_0 * gc_y[i] + gr_yy_0[i] * gc_y[i];

        grr_y_yz_0[i] = 2.0 * ts_z_0[i] * gfe2_0 + gr_z_0[i] * gfe_0 + 2.0 * ts_yz_0[i] * gfe_0 * gc_y[i] + gr_yz_0[i] * gc_y[i];

        grr_y_zz_0[i] = 2.0 * ts_zz_0[i] * gfe_0 * gc_y[i] + gr_zz_0[i] * gc_y[i];

        grr_z_xx_0[i] = 2.0 * ts_xx_0[i] * gfe_0 * gc_z[i] + gr_xx_0[i] * gc_z[i];

        grr_z_xy_0[i] = 2.0 * ts_xy_0[i] * gfe_0 * gc_z[i] + gr_xy_0[i] * gc_z[i];

        grr_z_xz_0[i] = 2.0 * ts_x_0[i] * gfe2_0 + gr_x_0[i] * gfe_0 + 2.0 * ts_xz_0[i] * gfe_0 * gc_z[i] + gr_xz_0[i] * gc_z[i];

        grr_z_yy_0[i] = 2.0 * ts_yy_0[i] * gfe_0 * gc_z[i] + gr_yy_0[i] * gc_z[i];

        grr_z_yz_0[i] = 2.0 * ts_y_0[i] * gfe2_0 + gr_y_0[i] * gfe_0 + 2.0 * ts_yz_0[i] * gfe_0 * gc_z[i] + gr_yz_0[i] * gc_z[i];

        grr_z_zz_0[i] = 4.0 * ts_z_0[i] * gfe2_0 + 2.0 * gr_z_0[i] * gfe_0 + 2.0 * ts_zz_0[i] * gfe_0 * gc_z[i] + gr_zz_0[i] * gc_z[i];
    }
}

} // t3rr2rec namespace

