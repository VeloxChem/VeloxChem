#include "ThreeCenterRR2PrimRecSD.hpp"

namespace t3rr2rec { // t3rr2rec namespace

auto
comp_prim_r_r2_sd(CSimdArray<double>& pbuffer, 
                  const size_t idx_gr_sd,
                  const size_t idx_sp,
                  const size_t idx_g_sp,
                  const size_t idx_sd,
                  const size_t idx_g_sd,
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

    // Set up components of auxiliary buffer : SP

    auto gr_0_x = pbuffer.data(idx_g_sp);

    auto gr_0_y = pbuffer.data(idx_g_sp + 1);

    auto gr_0_z = pbuffer.data(idx_g_sp + 2);

    // Set up components of auxiliary buffer : SD

    auto ts_0_xx = pbuffer.data(idx_sd);

    auto ts_0_xy = pbuffer.data(idx_sd + 1);

    auto ts_0_xz = pbuffer.data(idx_sd + 2);

    auto ts_0_yy = pbuffer.data(idx_sd + 3);

    auto ts_0_yz = pbuffer.data(idx_sd + 4);

    auto ts_0_zz = pbuffer.data(idx_sd + 5);

    // Set up components of auxiliary buffer : SD

    auto gr_0_xx = pbuffer.data(idx_g_sd);

    auto gr_0_xy = pbuffer.data(idx_g_sd + 1);

    auto gr_0_xz = pbuffer.data(idx_g_sd + 2);

    auto gr_0_yy = pbuffer.data(idx_g_sd + 3);

    auto gr_0_yz = pbuffer.data(idx_g_sd + 4);

    auto gr_0_zz = pbuffer.data(idx_g_sd + 5);

    // Set up components of targeted buffer : SD

    auto grr_x_0_xx = pbuffer.data(idx_gr_sd);

    auto grr_x_0_xy = pbuffer.data(idx_gr_sd + 1);

    auto grr_x_0_xz = pbuffer.data(idx_gr_sd + 2);

    auto grr_x_0_yy = pbuffer.data(idx_gr_sd + 3);

    auto grr_x_0_yz = pbuffer.data(idx_gr_sd + 4);

    auto grr_x_0_zz = pbuffer.data(idx_gr_sd + 5);

    auto grr_y_0_xx = pbuffer.data(idx_gr_sd + 6);

    auto grr_y_0_xy = pbuffer.data(idx_gr_sd + 7);

    auto grr_y_0_xz = pbuffer.data(idx_gr_sd + 8);

    auto grr_y_0_yy = pbuffer.data(idx_gr_sd + 9);

    auto grr_y_0_yz = pbuffer.data(idx_gr_sd + 10);

    auto grr_y_0_zz = pbuffer.data(idx_gr_sd + 11);

    auto grr_z_0_xx = pbuffer.data(idx_gr_sd + 12);

    auto grr_z_0_xy = pbuffer.data(idx_gr_sd + 13);

    auto grr_z_0_xz = pbuffer.data(idx_gr_sd + 14);

    auto grr_z_0_yy = pbuffer.data(idx_gr_sd + 15);

    auto grr_z_0_yz = pbuffer.data(idx_gr_sd + 16);

    auto grr_z_0_zz = pbuffer.data(idx_gr_sd + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_0_x, gr_0_xx, gr_0_xy, gr_0_xz, gr_0_y, gr_0_yy, gr_0_yz, gr_0_z, gr_0_zz, grr_x_0_xx, grr_x_0_xy, grr_x_0_xz, grr_x_0_yy, grr_x_0_yz, grr_x_0_zz, grr_y_0_xx, grr_y_0_xy, grr_y_0_xz, grr_y_0_yy, grr_y_0_yz, grr_y_0_zz, grr_z_0_xx, grr_z_0_xy, grr_z_0_xz, grr_z_0_yy, grr_z_0_yz, grr_z_0_zz, ts_0_x, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_y, ts_0_yy, ts_0_yz, ts_0_z, ts_0_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        grr_x_0_xx[i] = 2.0 * ts_0_x[i] * gfe2_0 + 2.0 * gr_0_x[i] * gfe_0 + ts_0_xx[i] * gfe_0 * gc_x[i] + gr_0_xx[i] * gc_x[i];

        grr_x_0_xy[i] = ts_0_y[i] * gfe2_0 + gr_0_y[i] * gfe_0 + ts_0_xy[i] * gfe_0 * gc_x[i] + gr_0_xy[i] * gc_x[i];

        grr_x_0_xz[i] = ts_0_z[i] * gfe2_0 + gr_0_z[i] * gfe_0 + ts_0_xz[i] * gfe_0 * gc_x[i] + gr_0_xz[i] * gc_x[i];

        grr_x_0_yy[i] = ts_0_yy[i] * gfe_0 * gc_x[i] + gr_0_yy[i] * gc_x[i];

        grr_x_0_yz[i] = ts_0_yz[i] * gfe_0 * gc_x[i] + gr_0_yz[i] * gc_x[i];

        grr_x_0_zz[i] = ts_0_zz[i] * gfe_0 * gc_x[i] + gr_0_zz[i] * gc_x[i];

        grr_y_0_xx[i] = ts_0_xx[i] * gfe_0 * gc_y[i] + gr_0_xx[i] * gc_y[i];

        grr_y_0_xy[i] = ts_0_x[i] * gfe2_0 + gr_0_x[i] * gfe_0 + ts_0_xy[i] * gfe_0 * gc_y[i] + gr_0_xy[i] * gc_y[i];

        grr_y_0_xz[i] = ts_0_xz[i] * gfe_0 * gc_y[i] + gr_0_xz[i] * gc_y[i];

        grr_y_0_yy[i] = 2.0 * ts_0_y[i] * gfe2_0 + 2.0 * gr_0_y[i] * gfe_0 + ts_0_yy[i] * gfe_0 * gc_y[i] + gr_0_yy[i] * gc_y[i];

        grr_y_0_yz[i] = ts_0_z[i] * gfe2_0 + gr_0_z[i] * gfe_0 + ts_0_yz[i] * gfe_0 * gc_y[i] + gr_0_yz[i] * gc_y[i];

        grr_y_0_zz[i] = ts_0_zz[i] * gfe_0 * gc_y[i] + gr_0_zz[i] * gc_y[i];

        grr_z_0_xx[i] = ts_0_xx[i] * gfe_0 * gc_z[i] + gr_0_xx[i] * gc_z[i];

        grr_z_0_xy[i] = ts_0_xy[i] * gfe_0 * gc_z[i] + gr_0_xy[i] * gc_z[i];

        grr_z_0_xz[i] = ts_0_x[i] * gfe2_0 + gr_0_x[i] * gfe_0 + ts_0_xz[i] * gfe_0 * gc_z[i] + gr_0_xz[i] * gc_z[i];

        grr_z_0_yy[i] = ts_0_yy[i] * gfe_0 * gc_z[i] + gr_0_yy[i] * gc_z[i];

        grr_z_0_yz[i] = ts_0_y[i] * gfe2_0 + gr_0_y[i] * gfe_0 + ts_0_yz[i] * gfe_0 * gc_z[i] + gr_0_yz[i] * gc_z[i];

        grr_z_0_zz[i] = 2.0 * ts_0_z[i] * gfe2_0 + 2.0 * gr_0_z[i] * gfe_0 + ts_0_zz[i] * gfe_0 * gc_z[i] + gr_0_zz[i] * gc_z[i];
    }
}

} // t3rr2rec namespace

