#include "ThreeCenterR2PrimRecPD.hpp"

namespace t3r2rec { // t3r2rec namespace

auto
comp_prim_r2_pd(CSimdArray<double>& pbuffer, 
                const size_t idx_g_pd,
                const size_t idx_sp,
                const size_t idx_sd,
                const size_t idx_ps,
                const size_t idx_pp,
                const size_t idx_pd,
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

    // Set up components of auxiliary buffer : PD

    auto ts_x_xx = pbuffer.data(idx_pd);

    auto ts_x_xy = pbuffer.data(idx_pd + 1);

    auto ts_x_xz = pbuffer.data(idx_pd + 2);

    auto ts_x_yy = pbuffer.data(idx_pd + 3);

    auto ts_x_yz = pbuffer.data(idx_pd + 4);

    auto ts_x_zz = pbuffer.data(idx_pd + 5);

    auto ts_y_xx = pbuffer.data(idx_pd + 6);

    auto ts_y_xy = pbuffer.data(idx_pd + 7);

    auto ts_y_xz = pbuffer.data(idx_pd + 8);

    auto ts_y_yy = pbuffer.data(idx_pd + 9);

    auto ts_y_yz = pbuffer.data(idx_pd + 10);

    auto ts_y_zz = pbuffer.data(idx_pd + 11);

    auto ts_z_xx = pbuffer.data(idx_pd + 12);

    auto ts_z_xy = pbuffer.data(idx_pd + 13);

    auto ts_z_xz = pbuffer.data(idx_pd + 14);

    auto ts_z_yy = pbuffer.data(idx_pd + 15);

    auto ts_z_yz = pbuffer.data(idx_pd + 16);

    auto ts_z_zz = pbuffer.data(idx_pd + 17);

    // Set up 0-6 components of targeted buffer : PD

    auto gr_x_xx = pbuffer.data(idx_g_pd);

    auto gr_x_xy = pbuffer.data(idx_g_pd + 1);

    auto gr_x_xz = pbuffer.data(idx_g_pd + 2);

    auto gr_x_yy = pbuffer.data(idx_g_pd + 3);

    auto gr_x_yz = pbuffer.data(idx_g_pd + 4);

    auto gr_x_zz = pbuffer.data(idx_g_pd + 5);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_x_xx, gr_x_xy, gr_x_xz, gr_x_yy, gr_x_yz, gr_x_zz, ts_0_x, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_y, ts_0_yy, ts_0_yz, ts_0_z, ts_0_zz, ts_x_0, ts_x_x, ts_x_xx, ts_x_xy, ts_x_xz, ts_x_y, ts_x_yy, ts_x_yz, ts_x_z, ts_x_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_x_xx[i] = 4.0 * ts_0_x[i] * gfe2_0 + 2.0 * ts_0_xx[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_0[i] * gfe2_0 + 4.0 * ts_x_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_x_xx[i] * gfe_0 + ts_x_xx[i] * rgc2_0;

        gr_x_xy[i] = 2.0 * ts_0_y[i] * gfe2_0 + 2.0 * ts_0_xy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_x_xy[i] * gfe_0 + ts_x_xy[i] * rgc2_0;

        gr_x_xz[i] = 2.0 * ts_0_z[i] * gfe2_0 + 2.0 * ts_0_xz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_xz[i] * gfe_0 + ts_x_xz[i] * rgc2_0;

        gr_x_yy[i] = 2.0 * ts_0_yy[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_0[i] * gfe2_0 + 4.0 * ts_x_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_x_yy[i] * gfe_0 + ts_x_yy[i] * rgc2_0;

        gr_x_yz[i] = 2.0 * ts_0_yz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_x_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_yz[i] * gfe_0 + ts_x_yz[i] * rgc2_0;

        gr_x_zz[i] = 2.0 * ts_0_zz[i] * gfe_0 * gc_x[i] + 2.0 * ts_x_0[i] * gfe2_0 + 4.0 * ts_x_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_x_zz[i] * gfe_0 + ts_x_zz[i] * rgc2_0;
    }

    // Set up 6-12 components of targeted buffer : PD

    auto gr_y_xx = pbuffer.data(idx_g_pd + 6);

    auto gr_y_xy = pbuffer.data(idx_g_pd + 7);

    auto gr_y_xz = pbuffer.data(idx_g_pd + 8);

    auto gr_y_yy = pbuffer.data(idx_g_pd + 9);

    auto gr_y_yz = pbuffer.data(idx_g_pd + 10);

    auto gr_y_zz = pbuffer.data(idx_g_pd + 11);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_y_xx, gr_y_xy, gr_y_xz, gr_y_yy, gr_y_yz, gr_y_zz, ts_0_x, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_y, ts_0_yy, ts_0_yz, ts_0_z, ts_0_zz, ts_y_0, ts_y_x, ts_y_xx, ts_y_xy, ts_y_xz, ts_y_y, ts_y_yy, ts_y_yz, ts_y_z, ts_y_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_y_xx[i] = 2.0 * ts_0_xx[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_0[i] * gfe2_0 + 4.0 * ts_y_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_y_xx[i] * gfe_0 + ts_y_xx[i] * rgc2_0;

        gr_y_xy[i] = 2.0 * ts_0_x[i] * gfe2_0 + 2.0 * ts_0_xy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_y_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_y_xy[i] * gfe_0 + ts_y_xy[i] * rgc2_0;

        gr_y_xz[i] = 2.0 * ts_0_xz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_y_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_xz[i] * gfe_0 + ts_y_xz[i] * rgc2_0;

        gr_y_yy[i] = 4.0 * ts_0_y[i] * gfe2_0 + 2.0 * ts_0_yy[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_0[i] * gfe2_0 + 4.0 * ts_y_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_y_yy[i] * gfe_0 + ts_y_yy[i] * rgc2_0;

        gr_y_yz[i] = 2.0 * ts_0_z[i] * gfe2_0 + 2.0 * ts_0_yz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_yz[i] * gfe_0 + ts_y_yz[i] * rgc2_0;

        gr_y_zz[i] = 2.0 * ts_0_zz[i] * gfe_0 * gc_y[i] + 2.0 * ts_y_0[i] * gfe2_0 + 4.0 * ts_y_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_y_zz[i] * gfe_0 + ts_y_zz[i] * rgc2_0;
    }

    // Set up 12-18 components of targeted buffer : PD

    auto gr_z_xx = pbuffer.data(idx_g_pd + 12);

    auto gr_z_xy = pbuffer.data(idx_g_pd + 13);

    auto gr_z_xz = pbuffer.data(idx_g_pd + 14);

    auto gr_z_yy = pbuffer.data(idx_g_pd + 15);

    auto gr_z_yz = pbuffer.data(idx_g_pd + 16);

    auto gr_z_zz = pbuffer.data(idx_g_pd + 17);

    #pragma omp simd aligned(gc_x, gc_y, gc_z, gr_z_xx, gr_z_xy, gr_z_xz, gr_z_yy, gr_z_yz, gr_z_zz, ts_0_x, ts_0_xx, ts_0_xy, ts_0_xz, ts_0_y, ts_0_yy, ts_0_yz, ts_0_z, ts_0_zz, ts_z_0, ts_z_x, ts_z_xx, ts_z_xy, ts_z_xz, ts_z_y, ts_z_yy, ts_z_yz, ts_z_z, ts_z_zz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];

        const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);

        const double gfe2_0 = gfe_0 * gfe_0;

        gr_z_xx[i] = 2.0 * ts_0_xx[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_0[i] * gfe2_0 + 4.0 * ts_z_x[i] * gfe_0 * gc_x[i] + 3.0 * ts_z_xx[i] * gfe_0 + ts_z_xx[i] * rgc2_0;

        gr_z_xy[i] = 2.0 * ts_0_xy[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_y[i] * gfe_0 * gc_x[i] + 2.0 * ts_z_x[i] * gfe_0 * gc_y[i] + 3.0 * ts_z_xy[i] * gfe_0 + ts_z_xy[i] * rgc2_0;

        gr_z_xz[i] = 2.0 * ts_0_x[i] * gfe2_0 + 2.0 * ts_0_xz[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_z[i] * gfe_0 * gc_x[i] + 2.0 * ts_z_x[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_xz[i] * gfe_0 + ts_z_xz[i] * rgc2_0;

        gr_z_yy[i] = 2.0 * ts_0_yy[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_0[i] * gfe2_0 + 4.0 * ts_z_y[i] * gfe_0 * gc_y[i] + 3.0 * ts_z_yy[i] * gfe_0 + ts_z_yy[i] * rgc2_0;

        gr_z_yz[i] = 2.0 * ts_0_y[i] * gfe2_0 + 2.0 * ts_0_yz[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_z[i] * gfe_0 * gc_y[i] + 2.0 * ts_z_y[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_yz[i] * gfe_0 + ts_z_yz[i] * rgc2_0;

        gr_z_zz[i] = 4.0 * ts_0_z[i] * gfe2_0 + 2.0 * ts_0_zz[i] * gfe_0 * gc_z[i] + 2.0 * ts_z_0[i] * gfe2_0 + 4.0 * ts_z_z[i] * gfe_0 * gc_z[i] + 3.0 * ts_z_zz[i] * gfe_0 + ts_z_zz[i] * rgc2_0;
    }

}

} // t3r2rec namespace

